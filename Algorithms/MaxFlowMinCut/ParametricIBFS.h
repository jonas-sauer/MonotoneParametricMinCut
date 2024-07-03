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

    struct OrphanBuckets {
        OrphanBuckets(const int n) :
           positionOfVertex_(n, -1), minBucket_(INFTY) {
        }

        inline void assertVertexInBucket(const Vertex vertex, const int dist) const noexcept {
            Assert(positionOfVertex_[vertex] != -1, "Vertex is not in bucket!");
            Assert(static_cast<size_t>(dist) < buckets_.size(), "Vertex is not in bucket!");
            Assert(static_cast<size_t>(positionOfVertex_[vertex]) < buckets_[dist].size(), "Vertex is not in bucket!");
            Assert(buckets_[dist][positionOfVertex_[vertex]] == vertex, "Vertex is not in bucket!");
        }

        inline void addVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex_[vertex] != -1) {
                assertVertexInBucket(vertex, dist);
                return;
            }
            if (static_cast<size_t>(dist) >= buckets_.size()) buckets_.resize(dist + 1);
            positionOfVertex_[vertex] = buckets_[dist].size();
            buckets_[dist].emplace_back(vertex);
            minBucket_ = std::min(minBucket_, dist);
        }

        inline void decreaseBucket(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(newDist < oldDist, "Distance has not decreased!");
            assertVertexInBucket(vertex, oldDist);
            const Vertex other = buckets_[oldDist].back();
            positionOfVertex_[other] = positionOfVertex_[vertex];
            buckets_[oldDist][positionOfVertex_[vertex]] = other;
            buckets_[oldDist].pop_back();
            positionOfVertex_[vertex] = buckets_[newDist].size();
            buckets_[newDist].emplace_back(vertex);
            if (static_cast<size_t>(oldDist) == buckets_.size() - 1) {
                while (!buckets_.empty() && buckets_.back().empty()) buckets_.pop_back();
            }
            minBucket_ = std::min(minBucket_, newDist);
        }

        inline bool empty() const noexcept {
            return buckets_.empty();
        }

        inline Vertex pop() noexcept {
            Assert(!empty(), "Buckets are empty!");
            const Vertex vertex = buckets_[minBucket_].back();
            buckets_[minBucket_].pop_back();
            positionOfVertex_[vertex] = -1;
            while (static_cast<size_t>(minBucket_) < buckets_.size() && buckets_[minBucket_].empty()) minBucket_++;
            if (static_cast<size_t>(minBucket_) == buckets_.size()) {
                buckets_.clear();
                minBucket_ = INFTY;
            }
            return vertex;
        }

        std::vector<std::vector<Vertex>> buckets_;
        std::vector<int> positionOfVertex_;
        int minBucket_;
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
        orphans_(n),
        threePassOrphans_(n),
        excess_at_vertex_(n, FlowFunction(0)) {
        for (const Vertex vertex : graph_.vertices()) {
            currentEdge_[vertex] = graph_.beginEdgeFrom(vertex);
        }
    };

    inline void run() noexcept {
        initialize();
        double alpha = alphaMin_;
        while (pmf::doubleLessThanAbs(alpha, alphaMax_)) {
            if (alphaQ_.empty()) return;
            assert(alphaQ_.front()->value_ > alpha);
            alpha = alphaQ_.front()->value_;
            if (alpha > alphaMax_) return;
            updateTree(alpha);
            //TODO Hybrid adoption
            //reconnectTree(alpha);
            reconnectTreeThreePass(alpha);
            drainExcess(alpha);
        }
    }

    inline const std::vector<double>& getBreakpoints() const noexcept {
        return thetaBreakpoints_;
    };

    inline const std::vector<double>& getVertexThetas() const noexcept {
        return thetaByVertex_;
    }

    inline std::vector<Vertex> getSinkComponent(const double alpha) const noexcept {
        std::vector<Vertex> sinkComponent;
        for (const Vertex vertex : graph_.vertices()) {
            if (thetaByVertex_[vertex] <= alpha) continue;
            sinkComponent.emplace_back(vertex);
        }
        return sinkComponent;
    }

    inline double getFlowValue(const double alpha) const noexcept {
        double flow = 0;
        for (const Vertex from : graph_.vertices()) {
            if (thetaByVertex_[from] > alpha) continue;
            for (const Edge edge : graph_.edgesFrom(from)) {
                const Vertex to = graph_.get(ToVertex, edge);
                if (thetaByVertex_[to] <= alpha) continue;
                flow += instance_.getCapacity(edge, alpha);
            }
        }
        return flow;
    }

private:
    inline void initialize() noexcept {
        initialFlow.run();
        const std::vector<double>& initialResidualCapacity = initialFlow.getCleanResidualCapacities();
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
                        residualCapacity_[e] = graph_.get(Capacity, e) - FlowFunction(graph_.get(Capacity, e).eval(alphaMin_) - initialResidualCapacity[e]);
                        residualCapacity_[revE] = graph_.get(Capacity, revE) - FlowFunction(graph_.get(Capacity, revE).eval(alphaMin_) - initialResidualCapacity[revE]);
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

    inline void initializeSinkTree(const std::vector<double>& initialResidualCapacity) noexcept {
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
                queue.emplace_back(w);
            }
        }
    }

    inline void saturateEdgeInitial(const Edge e, const Edge revE, const Vertex to) noexcept {
        const FlowFunction& capacity = graph_.get(Capacity, e);
        residualCapacity_[e] = FlowFunction(0);
        residualCapacity_[revE] = capacity + graph_.get(Capacity, revE);
        excess_at_vertex_[to] += capacity - FlowFunction(capacity.eval(alphaMin_));
        excessVertices_.addVertex(to, dist_[to]);
    }

    inline void updateTree(const double nextAlpha) noexcept {
        assert(orphans_.empty());
        assert(threePassOrphans_.empty());
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

    inline void reconnectTree(const double nextAlpha) noexcept {
        std::vector<Vertex> moved;
        while (!orphans_.empty()) {
            const Vertex v = orphans_.pop();
            excessVertices_.removeVertex(v, dist_[v]);
            if (adoptWithSameDist(v, nextAlpha)) continue;
            treeData_.removeChildren(v, [&](const Vertex child, const Edge e) {
                removeTreeEdge(e, child, v, nextAlpha);
            });
            excessVertices_.removeVertex(v, dist_[v]);
            if (adoptWithNewDist(v, nextAlpha)) continue;
            moved.emplace_back(v);
            dist_[v] = INFTY;
            thetaByVertex_[v] = nextAlpha;
            if (thetaBreakpoints_.back() != nextAlpha) {
                thetaBreakpoints_.emplace_back(nextAlpha);
            }
        }
    }

    inline bool adoptWithSameDist(const Vertex v, const double nextAlpha) noexcept {
        for (Edge e = currentEdge_[v]; e < graph_.endEdgeFrom(v); e++) {
            if (!isEdgeResidual(e, nextAlpha)) continue;
            const Vertex to = graph_.get(ToVertex, e);
            if (!isEdgeAdmissible(v, to)) continue;
            excessVertices_.addVertex(v, dist_[v]);
            treeData_.addVertex(to, v, e);
            currentEdge_[v] = e;
            return true;
        }
        return false;
    }

    inline bool adoptWithNewDist(const Vertex v, const double nextAlpha) noexcept {
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
        dist_[v] = d_min + 1;
        excessVertices_.addVertex(v, dist_[v]);
        treeData_.addVertex(v_min, v, e_min);
        currentEdge_[v] = e_min;
        return true;
    }

    inline void reconnectTreeThreePass(const double nextAlpha) noexcept {
        // Pass 1: Try to adopt orphans without changing their distance.
        reconnectTreeFirstPass(nextAlpha);
        std::vector<Vertex> moved;
        // Passes 2 and 3: Go through orphans in increasing order of distance.
        while (!threePassOrphans_.empty()) {
            const Vertex v = threePassOrphans_.pop();
            // edgeToParent_[v] is set iff pass 2 was successful.
            if (treeData_.edgeToParent_[v] == noEdge) {
                // Pass 2: Try to find a non-orphan parent for v with minimal distance.
                // If this increases the distance of v, postpone pass 3 until the new bucket of v is scanned.
                //std::cout << "Second pass of " << v << " with distance " << dist_[v] << std::endl;
                if (reconnectTreeSecondPass(v, nextAlpha)) {
                    //std::cout << "Distance increased to " << dist_[v] << std::endl;
                    continue;
                }
            }
            if (treeData_.edgeToParent_[v] != noEdge) {
                // Pass 3 (only if pass 2 succeeded): Finalize the adoption and update the distances of potential children.
                //std::cout << "Third pass of " << v << " with distance " << dist_[v] << std::endl;
                reconnectTreeThirdPass(v, nextAlpha);
            } else {
                // If pass 2 failed, v leaves the sink component.
                moved.emplace_back(v); //TODO: Do this properly
            }
        }
        for (const Vertex v : moved) {
            if (treeData_.edgeToParent_[v] != noEdge) continue;
            excessVertices_.removeVertex(v, dist_[v]);
            dist_[v] = INFTY;
            thetaByVertex_[v] = nextAlpha;
            if (thetaBreakpoints_.back() != nextAlpha) {
                thetaBreakpoints_.emplace_back(nextAlpha);
            }
        }
    }

    // Go through orphans in increasing order of distance.
    // If they can be adopted with the same distance, do so.
    // Otherwise, increment their distance and add them to the orphan buckets for passes 2 and 3.
    inline void reconnectTreeFirstPass(const double nextAlpha) noexcept {
        while (!orphans_.empty()) {
            const Vertex v = orphans_.pop();
            excessVertices_.removeVertex(v, dist_[v]);
            if (adoptWithSameDist(v, nextAlpha)) continue;
            //std::cout << "First pass of " << v << " with distance " << dist_[v] << " failed!" << std::endl;
            treeData_.removeChildren(v, [&](const Vertex child, const Edge e) {
                removeTreeEdge(e, child, v, nextAlpha);
            });
            //TODO Clean up this mess
            excessVertices_.removeVertex(v, dist_[v]);
            dist_[v]++;
            threePassOrphans_.addVertex(v, dist_[v]);
        }
    }

    // Try to find a non-orphan parent with minimal distance for v.
    // This adoption is not necessarily optimal yet - it may be improved once other orphans re-enter the tree.
    // Return true iff the adoption was successful and the distance increased.
    inline bool reconnectTreeSecondPass(const Vertex v, const double nextAlpha) noexcept {
        uint d_min = INFTY;
        Edge e_min = noEdge;
        Vertex v_min = noVertex;

        for (const Edge e : graph_.edgesFrom(v)) {
            if (!isEdgeResidual(e, nextAlpha)) continue;
            const Vertex to = graph_.get(ToVertex, e);
            // Don't allow unadopted orphans as parents.
            // To distinguish between adopted and unadopted orphans, we temporarily set edgeToParent_ during the adoption.
            // This will not necessarily be the final parent edge, so the child relation is not updated yet.
            if (treeData_.edgeToParent_[to] == noEdge) continue;
            if (dist_[to] < d_min) {
                e_min = e;
                d_min = dist_[to];
                v_min = to;
            }
        }

        if (d_min == INFTY) {
            //std::cout << "Adoption failed" << std::endl;
            return false;
        }
        treeData_.edgeToParent_[v] = e_min;
        //std::cout << "Set parent to " << v_min << " with " << dist_[v_min];
        if (d_min + 1 > dist_[v]) {
            dist_[v] = d_min + 1;
            threePassOrphans_.addVertex(v, dist_[v]);
            return true;
        }
        return false;
    }

    inline void reconnectTreeThirdPass(const Vertex v, const double nextAlpha) noexcept {
        // dist_[v] is correct at this point.
        // Find the first admissible parent edge to ensure that currentEdge_ is set correctly.
        for (const Edge e : graph_.edgesFrom(v)) {
            if (!isEdgeResidual(e, nextAlpha)) continue;
            const Vertex to = graph_.get(ToVertex, e);
            if (!isEdgeAdmissible(v, to)) continue;
            if (treeData_.edgeToParent_[to] == noEdge) continue;
            excessVertices_.addVertex(v, dist_[v]);
            treeData_.addVertex(to, v, e);
            currentEdge_[v] = e;
            break;
        }

        // Update the distances of all potential children.
        for (const Edge e : graph_.edgesFrom(v)) {
            const Vertex from = graph_.get(ToVertex, e);
            if (dist_[from] == INFTY || from == sink_) continue;
            const bool isFree = treeData_.edgeToParent_[from] == noEdge;
            const bool isDistanceGreater = dist_[from] > dist_[v] + 1;
            if (!isFree && !isDistanceGreater) continue;
            const Edge rev = graph_.get(ReverseEdge, e);
            if (!isEdgeResidual(rev, nextAlpha)) continue;
            if (isDistanceGreater) {
                //std::cout << "Decrease distance of " << from << " from " << dist_[from] << " to " << dist_[v] + 1 << std::endl;
                threePassOrphans_.decreaseBucket(from, dist_[from], dist_[v] + 1);
            } else {
                assert(isFree);
                //std::cout << "Re-add orphan " << from << " from " << dist_[from] << " to " << dist_[v] + 1 << std::endl;
                if (dist_[from] == dist_[v]) {
                    // TODO This is incorrect but the paper claims that this is done?
                    //threePassOrphans_.increaseBucket(from, dist_[from], dist_[v] + 1);
                    continue;
                } else if (dist_[from] < dist_[v]) {
                    threePassOrphans_.addVertex(from, dist_[v] + 1);
                }
            }
            treeData_.edgeToParent_[from] = rev;
            dist_[from] = dist_[v] + 1;
        }
    }

    inline void drainExcess(const double nextAlpha) noexcept {
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

    inline void removeTreeEdge(const Edge e, const Vertex from, const Vertex to, const double nextAlpha) noexcept {
        const Edge rev = graph_.get(ReverseEdge, e);
        const FlowFunction oldResidualCapacity = residualCapacity_[e];
        const FlowFunction newResidualCapacity(residualCapacity_[e].eval(nextAlpha));
        const FlowFunction add = newResidualCapacity - oldResidualCapacity;
        //Don't add orphan to excessVertices_ yet. Wait until it has been adopted.
        excess_at_vertex_[from] += add;
        excess_at_vertex_[to] -= add;
        residualCapacity_[e] = newResidualCapacity;
        residualCapacity_[rev] -= add;
        excessVertices_.addVertex(to, dist_[to]);
        clearRootAlpha(from);
        orphans_.addVertex(from, dist_[from]);
        assert(dist_[from] != INFTY);
    }

    inline double getNextZeroCrossing(const Edge e, const double alpha) const noexcept {
        return residualCapacity_[e].getNextZeroCrossing(alpha);
    }

    inline void recalculateRootAlpha(const Vertex v, const Edge e, const double alpha) noexcept {
        const double oldValue = rootAlpha_[v].value_;
        rootAlpha_[v].value_ = getNextZeroCrossing(e, alpha);
        //assert(rootAlpha_[v].value_ > alpha);
        if (rootAlpha_[v].value_ == INFTY) {
            if (oldValue < INFTY) {
                alphaQ_.remove(&rootAlpha_[v]);
            }
        }
        else {
            alphaQ_.push(&rootAlpha_[v]);
        }
    }

    inline void clearRootAlpha(const Vertex v) noexcept {
        if (rootAlpha_[v].value_ == INFTY) return;
        rootAlpha_[v].value_ = INFTY;
        alphaQ_.remove(&rootAlpha_[v]);
    }

    inline bool isEdgeAdmissible(const Vertex from, const Vertex to) const noexcept {
        return dist_[from] == dist_[to] + 1;
    }

    inline bool isEdgeResidual(const Edge edge, const double alpha) const noexcept {
        return pmf::doubleIsPositive(residualCapacity_[edge].eval(alpha));
    }

    inline void checkTree(const bool allowOrphans, const double alpha) const noexcept {
        for (const Vertex v : graph_.vertices()) {
            if (dist_[v] == INFTY || v == sink_) continue;
            const Edge edge = treeData_.edgeToParent_[v];
            if (edge == noEdge) {
                assert(allowOrphans);
                continue;
            }
            const Vertex parent = graph_.get(ToVertex, edge);
            assert(dist_[parent] == dist_[v] - 1);
            for (const Edge e : graph_.edgesFrom(v)) {
                const Vertex to = graph_.get(ToVertex, e);
                if (!isEdgeResidual(e, alpha)) continue;
                assert(dist_[to] >= dist_[v] - 1);
            }

        }
    }

    inline void checkQueue() noexcept {
        for (const Vertex v : graph_.vertices()) {
            if (treeData_.edgeToParent_[v] == noEdge) continue;
            if (rootAlpha_[v] == INFTY) continue;
            assert(alphaQ_.contains(&rootAlpha_[v]));
        }
    }

    inline void checkExcessBuckets() noexcept {
        for (size_t i = 0; i < excessVertices_.buckets_.size(); i++) {
            for (const Vertex v : excessVertices_.buckets_[i]) {
                assert(dist_[v] == i);
            }
        }
    }

    inline void checkCapacityConstraints(const double alpha) const noexcept {
        for (const Edge edge : graph_.edges()) {
            Assert(!pmf::isNumberNegative(residualCapacity_[edge].eval(alpha)), "Capacity constraint violated!");
        }
    }

    inline void checkFlowConservation(const bool allowExcesses = false) const noexcept {
        for (const Vertex vertex : graph_.vertices()) {
            checkFlowConservation(vertex, allowExcesses);
        }
    }

    inline void checkFlowConservation(const double alpha, const bool allowExcesses = false) const noexcept {
        for (const Vertex vertex : graph_.vertices()) {
            checkFlowConservation(vertex, alpha, allowExcesses);
        }
    }

    inline FlowFunction getInflow(const Vertex vertex) const noexcept {
        FlowFunction inflow(0);
        for (const Edge edge : graph_.edgesFrom(vertex)) {
            const Edge reverseEdge = graph_.get(ReverseEdge, edge);
            inflow += instance_.getCapacity(reverseEdge) - residualCapacity_[reverseEdge];
        }
        return inflow;
    }

    inline void checkFlowConservation(const Vertex vertex, const bool allowExcesses = false) const noexcept {
        if (vertex == source_ || vertex == sink_) return;
        if (!allowExcesses) Assert(excess_at_vertex_[vertex] == FlowFunction(0), "Vertex " << vertex << " has excess!");
        const FlowFunction inflow = getInflow(vertex);
        const FlowFunction netInflow = inflow - excess_at_vertex_[vertex];
        Assert(netInflow == FlowFunction(0), "Flow conservation not fulfilled!");
    }

    inline void checkFlowConservation(const Vertex vertex, const double alpha, const bool allowExcesses = false) const noexcept {
        if (vertex == source_ || vertex == sink_) return;
        if (!allowExcesses) Assert(excess_at_vertex_[vertex] == FlowFunction(0), "Vertex " << vertex << " has excess!");
        const FlowFunction inflow = getInflow(vertex);
        const FlowFunction netInflow = inflow - excess_at_vertex_[vertex];
        Assert(pmf::doubleEqualAbs(netInflow.eval(alpha), 0), "Flow conservation not fulfilled!");
    }

    inline void checkDisconnected(const Vertex v, const double alpha) const noexcept {
        for (const Edge edge : graph_.edgesFrom(v)) {
            if (!isEdgeResidual(edge, alpha)) continue;
            const Vertex to = graph_.get(ToVertex, edge);
            Assert(dist_[to] == INFTY, "Vertex " << v << " is not disconnected from sink component!");
        }
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

    OrphanBuckets orphans_;
    OrphanBuckets threePassOrphans_;

    std::vector<FlowFunction> excess_at_vertex_;
};
