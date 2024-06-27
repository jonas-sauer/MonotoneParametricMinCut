#pragma once

#include <queue>
#include <vector>

#include "IBFS.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/FlowUtils.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Meta.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class ParametricIBFS {
public:
    using FlowFunction = FLOW_FUNCTION;
    using FlowGraph = ParametricFlowGraph<FlowFunction>;
    using StaticWrapper = ParametricToStaticMaxFlowInstanceWrapper<FlowFunction>;
    using IBFSType = IBFS<StaticWrapper>;

private:
    struct ExcessBuckets {
        ExcessBuckets(const int n) :
            positionOfVertex_(n, -1) {
        }

        inline void addVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex_[vertex] != -1) return;
            if (static_cast<size_t>(dist) >= buckets_.size()) buckets_.resize(dist + 1);
            positionOfVertex_[vertex] = buckets_[dist].size();
            buckets_[dist].emplace_back(vertex);
        }

        inline void removeVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex_[vertex] == -1) return;
            const Vertex other = buckets_[dist].back();
            positionOfVertex_[other] = positionOfVertex_[vertex];
            buckets_[dist][positionOfVertex_[vertex]] = other;
            buckets_[dist].pop_back();
            positionOfVertex_[vertex] = -1;
            while (!buckets_.empty() && buckets_.back().empty()) buckets_.pop_back();
        }

        inline void increaseBucket(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(newDist > oldDist, "Distance has not increased!");
            Assert(static_cast<size_t>(positionOfVertex_[vertex]) < buckets_[oldDist].size(), "Vertex is not in bucket!");
            Assert(buckets_[oldDist][positionOfVertex_[vertex]] == vertex, "Vertex is not in bucket!");
            const Vertex other = buckets_[oldDist].back();
            positionOfVertex_[other] = positionOfVertex_[vertex];
            buckets_[oldDist][positionOfVertex_[vertex]] = other;
            buckets_[oldDist].pop_back();
            if (static_cast<size_t>(newDist) >= buckets_.size()) buckets_.resize(newDist + 1);
            positionOfVertex_[vertex] = buckets_[newDist].size();
            buckets_[newDist].emplace_back(vertex);
        }

        inline bool empty() const noexcept {
            return buckets_.empty();
        }

        inline Vertex pop() noexcept {
            Assert(!empty(), "Buckets are empty!");
            const Vertex vertex = buckets_.back().back();
            buckets_.back().pop_back();
            positionOfVertex_[vertex] = -1;
            while (!buckets_.empty() && buckets_.back().empty()) buckets_.pop_back();
            return vertex;
        }

        std::vector<std::vector<Vertex>> buckets_;
        std::vector<int> positionOfVertex_;
    };

    struct TreeData {
        TreeData(const size_t n) :
            edgeToParent_(n, noEdge),
            children_(n),
            childIndex_(n, -1) {
        }

        void addVertex(const Vertex parent, const Vertex child, const Edge edge) {
            edgeToParent_[child] = edge;
            childIndex_[child] = children_[parent].size();
            children_[parent].emplace_back(child);
        }

        void removeChild(const Vertex parent, const Vertex child) {
            const uint index = childIndex_[child];
            const Vertex back = children_[parent].back();
            childIndex_[back] = index;
            childIndex_[child] = -1;
            children_[parent][index] = back;
            children_[parent].pop_back();
            edgeToParent_[child] = noEdge;
        }

        template<typename FUNCTION>
        void removeChildren(const Vertex parent, const FUNCTION& callback) {
            for (const Vertex child : children_[parent]) {
                callback(child, edgeToParent_[child]);
                edgeToParent_[child] = noEdge;
                childIndex_[child] = -1;
            }
            children_[parent].clear();
        }

        std::vector<Edge> edgeToParent_;
        std::vector<std::vector<Vertex>> children_;
        std::vector<uint> childIndex_;
    };

    struct RootAlphaLabel : public ExternalKHeapElement {
        RootAlphaLabel() : ExternalKHeapElement(), value_(INFTY) {}

        inline bool hasSmallerKey(const RootAlphaLabel* other) const noexcept {
            return value_ < other->value_;
        }
        double value_;
    };

public:
    ParametricIBFS(const ParametricMaxFlowInstance<FlowFunction>& instance) :
        instance_(instance),
        graph_(instance.graph),
        source_(instance.source),
        sink_(instance.sink),
        alphaMin_(instance.alphaMin),
        alphaMax_(instance.alphaMax),
        n(graph_.numVertices()),
        wrapper(instance, instance.alphaMin),
        initialFlow(wrapper),
        residualCapacity_(instance_.getCurrentCapacities()),
        dist_(n, INFTY),
        excessVertices_(n),
        thetaByVertex_(n, INFTY),
        thetaBreakpoints_(1, alphaMin_),
        treeData_(n),
        currentEdge_(n, noEdge),
        rootAlpha_(n),
        alphaQ_(n),
        excess_at_vertex_(n, FlowFunction(0)) {
        for (const Vertex vertex : graph_.vertices()) {
            currentEdge_[vertex] = graph_.beginEdgeFrom(vertex);
        }
    };

    void run() {
        initialize();
        double alpha = alphaMin_;
        while (pmf::doubleLessThanAbs(alpha, alphaMax_)) {
            assert(!alphaQ_.empty());
            assert(alphaQ_.front()->value_ > alpha);
            alpha = alphaQ_.front()->value_;
            if (alpha > alphaMax_) return;
            updateTree(alpha);
            reconnectTree(alpha);
            drainExcess(alpha);
        }
    }

    const std::vector<double>& getBreakpoints() const {
        return thetaBreakpoints_;
    };

    const std::vector<double>& getVertexThetas() const {
        return thetaByVertex_;
    }

private:
    void initialize() {
        initialFlow.run();
        const std::vector<double>& initialResidualCapacity = initialFlow.getResidualCapacities();
        // TODO Obtain source and sink component from initialFlow
        initializeSinkTree(initialResidualCapacity);

        for (const Vertex from : graph_.vertices()) {
            if (from == sink_) continue;
            if (dist_[from] == INFTY) {
                thetaByVertex_[from] = alphaMin_;
            }
            for (const Edge e : graph_.edgesFrom(from)) {
                const Vertex to = graph_.get(ToVertex, e);
                if (dist_[to] == INFTY) continue;
                const Edge revE = graph_.get(ReverseEdge, e);
                if (dist_[from] == INFTY) {
                    if (dist_[to] == INFTY) continue;
                    // Cut edge
                    saturateEdgeInitial(e, revE, to);
                } else if (to == sink_) {
                    if (initialResidualCapacity[e] > 0) {
                        residualCapacity_[e] = FlowFunction(initialResidualCapacity[e]);
                        residualCapacity_[revE] = FlowFunction(initialResidualCapacity[revE]);
                    } else {
                        saturateEdgeInitial(e, revE, to);
                    }
                } else {
                    residualCapacity_[e] = FlowFunction(initialResidualCapacity[e]);
                }
            }
        }

        drainExcess(alphaMin_);
    }

    void saturateEdgeInitial(const Edge e, const Edge revE, const Vertex to) {
        const FlowFunction& capacity = graph_.get(Capacity, e);
        residualCapacity_[e] = FlowFunction(0);
        residualCapacity_[revE] = capacity + graph_.get(Capacity, revE);
        excess_at_vertex_[to] += capacity - FlowFunction(capacity.eval(alphaMin_));
        excessVertices_.addVertex(to, dist_[to]);
    }

    void initializeSinkTree(const std::vector<double>& initialResidualCapacity) {
        std::deque<Vertex> queue(1, sink_);
        dist_[sink_] = 0;

        while (!queue.empty()) {
            const Vertex v = queue.front();
            queue.pop_front();

            for (const Edge e : graph_.edgesFrom(v)) {
                const Edge revE = graph_.get(ReverseEdge, e);
                if (!pmf::doubleIsPositive(initialResidualCapacity[revE])) continue;
                const Vertex w = graph_.get(ToVertex, e);
                if (dist_[w] != INFTY) continue;
                dist_[w] = dist_[v] + 1;
                treeData_.addVertex(v, w, revE);
                alphaQ_.push(rootAlpha_[w]);
                queue.emplace_back(w);
            }
        }
    }

    void updateTree(const double nextAlpha) {
        assert(orphans_.empty());
        while (!alphaQ_.empty() && alphaQ_.front()->value_ == nextAlpha) {
            const Vertex v(alphaQ_.front() - &(rootAlpha_[0]));
            const Edge e = treeData_.edgeToParent_[v];
            assert(e != noEdge);
            assert(!isEdgeResidual(e, nextAlpha));
            const Vertex parent = graph_.get(ToVertex, e);
            removeTreeEdge(e, v, parent, nextAlpha);
            treeData_.removeChild(parent, v);
        }
    }

    bool adoptWithSameDist(const Vertex v, const double nextAlpha) {
        for (Edge e = currentEdge_[v]; e < graph_.endEdgeFrom(v); e++) {
            if (!isEdgeResidual(e, nextAlpha)) continue;
            const Vertex to = graph_.get(ToVertex, e);
            if (!isEdgeAdmissible(v, to)) continue;
            treeData_.addVertex(to, v, e);
            alphaQ_.push(&(rootAlpha_[v]));
            currentEdge_[v] = e;
            return true;
        }
        return false;
    }

    bool adoptWithNewDist(const Vertex v, const double nextAlpha) {
        uint d_min = INFTY;
        Edge e_min = noEdge;
        Vertex v_min = noVertex;

        for (const Edge e : graph_.edgesFrom(v)) {
            if (!isEdgeResidual(e, nextAlpha)) continue;
            const Vertex to = graph_.get(ToVertex, e);
            if (dist_[to] < d_min) {
                e_min = e;
                d_min = dist_[to];
                v_min = to;
            }
        }

        if (d_min > static_cast<uint>(n) - 1) {
            return false;
        }

        assert(d_min >= dist_[v]);
        excessVertices_.increaseBucket(v, dist_[v], d_min + 1); //TODO Only if it has excess
        dist_[v] = d_min + 1;
        treeData_.addVertex(v_min, v, e_min);
        alphaQ_.push(&(rootAlpha_[v]));
        currentEdge_[v] = e_min;
        return true;
    }

    void reconnectTree(const double nextAlpha) {
        while (!orphans_.empty()) {
            const Vertex v = orphans_.back();
            orphans_.pop_back();

            if (adoptWithSameDist(v, nextAlpha)) continue;

            //Children become orphaned, as they may have better options for a parent now
            treeData_.removeChildren(v, [&](const Vertex child, const Edge e) {
                removeTreeEdge(e, child, v, nextAlpha);
            });

            if (adoptWithNewDist(v, nextAlpha)) continue;

            excessVertices_.removeVertex(v, dist_[v]);
            dist_[v] = INFTY;
            thetaByVertex_[v] = nextAlpha;
            if (thetaBreakpoints_.back() != nextAlpha) {
                thetaBreakpoints_.emplace_back(nextAlpha);
            }
        }
    }

    void checkQueue() {
        for (const Vertex v : graph_.vertices()) {
            if (treeData_.edgeToParent_[v] == noEdge) continue;
            assert(alphaQ_.contains(&rootAlpha_[v]));
        }
    }

    void drainExcess(const double nextAlpha) {
        while (!excessVertices_.empty()) {
            const Vertex v = excessVertices_.pop();
            if (v == sink_) continue;
            assert(dist_[v] != INFTY);

            const Edge e = treeData_.edgeToParent_[v];
            const Edge revE = graph_.get(ReverseEdge, e);
            const Vertex w = graph_.get(ToVertex, e);

            residualCapacity_[e] -= excess_at_vertex_[v];
            residualCapacity_[revE] += excess_at_vertex_[v];
            if (w != sink_) {
                excess_at_vertex_[w] += excess_at_vertex_[v];
                excessVertices_.addVertex(w, dist_[w]);
            }
            excess_at_vertex_[v] = FlowFunction(0);
            recalculateRootAlpha(v, e, nextAlpha);
        }
    }

    void removeTreeEdge(const Edge e, const Vertex from, const Vertex to, const double nextAlpha) {
        const FlowFunction& oldResidualCapacity = residualCapacity_[e];
        const FlowFunction newResidualCapacity(residualCapacity_[e].eval(nextAlpha));
        excess_at_vertex_[from] += newResidualCapacity - oldResidualCapacity;
        excessVertices_.addVertex(from, dist_[from]);
        excess_at_vertex_[to] += oldResidualCapacity - newResidualCapacity;
        excessVertices_.addVertex(to, dist_[to]);
        setResidualCapacity(e, newResidualCapacity);
        clearRootAlpha(from);
        orphans_.emplace_back(from);
        assert(dist_[from] != INFTY);
    }

    void setResidualCapacity(const Edge e, const FlowFunction& newResidualCapacity) {
        residualCapacity_[e] = newResidualCapacity;
        const Edge revE = graph_.get(ReverseEdge, e);
        residualCapacity_[revE] = -newResidualCapacity;
    }

    double getNextZeroCrossing(const Edge e, const double alpha) const {
        return residualCapacity_[e].getNextZeroCrossing(alpha);
    }

    void recalculateRootAlpha(const Vertex v, const Edge e, const double alpha) {
        rootAlpha_[v].value_ = getNextZeroCrossing(e, alpha);
        alphaQ_.push(&rootAlpha_[v]);
    }

    void clearRootAlpha(const Vertex v) {
        rootAlpha_[v].value_ = INFTY;
        alphaQ_.remove(&rootAlpha_[v]); //TODO Only if value was not already infty. Ensure that vertices with rootAlpha infty are never in the queue
    }

    bool isEdgeAdmissible(const Vertex from, const Vertex to) const {
        return dist_[from] == dist_[to] + 1;
    }

    bool isEdgeResidual(const Edge edge, const double alpha) const {
        return pmf::doubleIsPositive(residualCapacity_[edge].eval(alpha));
    }

private:
    const ParametricMaxFlowInstance<FlowFunction>& instance_;
    const FlowGraph& graph_;
    const Vertex& source_, sink_;
    const double& alphaMin_, alphaMax_;
    const int n;

    StaticWrapper wrapper;
    IBFSType initialFlow;

    std::vector<FlowFunction> residualCapacity_;
    std::vector<uint> dist_;
    ExcessBuckets excessVertices_;

    std::vector<double> thetaByVertex_;
    std::vector<double> thetaBreakpoints_;
    TreeData treeData_;
    std::vector<Edge> currentEdge_;

    std::vector<RootAlphaLabel> rootAlpha_;
    ExternalKHeap<2, RootAlphaLabel> alphaQ_;

    std::vector<Vertex> orphans_;

    std::vector<FlowFunction> excess_at_vertex_;
};
