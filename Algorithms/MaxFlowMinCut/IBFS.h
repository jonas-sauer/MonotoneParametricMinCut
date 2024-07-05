#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<typename MAX_FLOW_INSTANCE>
class IBFS {
public:
    using MaxFlowInstance = MAX_FLOW_INSTANCE;
    using FlowType = MaxFlowInstance::FlowType;
    using GraphType = MaxFlowInstance::GraphType;

private:
    struct TreeData {
        TreeData(const size_t n) :
            parentEdge(n, noEdge),
            parentVertex(n , noVertex),
            firstChild(n, noVertex),
            nextSibling(n, noVertex),
            prevSibling(n, noVertex) {
        }

        inline void addVertex(const Vertex parent, const Vertex child, const Edge edge) noexcept {
            parentEdge[child] = edge;
            parentVertex[child] = parent;
            prevSibling[child] = noVertex;
            nextSibling[child] = firstChild[parent];
            if (nextSibling[child] != noVertex)
                prevSibling[nextSibling[child]] = child;
            firstChild[parent] = child;
        }

        inline void removeVertex(const Vertex vertex) noexcept {
            if (nextSibling[vertex] != noVertex)
                prevSibling[nextSibling[vertex]] = prevSibling[vertex];
            if (prevSibling[vertex] == noVertex)
                firstChild[parentVertex[vertex]] = nextSibling[vertex];
            else
                nextSibling[prevSibling[vertex]] = nextSibling[vertex];
            parentEdge[vertex] = noEdge;
            parentVertex[vertex] = noVertex;
        }

        template<typename FUNCTION>
        inline void removeChildren(const Vertex parent, const FUNCTION& callback) noexcept {
            Vertex child = firstChild[parent];
            while (child != noVertex) {
                const Vertex next = nextSibling[child];
                parentEdge[child] = noEdge;
                parentVertex[child] = noVertex;
                callback(child);
                child = next;
            }
            firstChild[parent] = noVertex;
        }

        // Always directed towards sink: (parent, child) for source tree, (child, parent) for sink tree
        std::vector<Edge> parentEdge;
        std::vector<Vertex> parentVertex;
        std::vector<Vertex> firstChild;
        std::vector<Vertex> nextSibling;
        std::vector<Vertex> prevSibling;
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

    struct Cut {
        Cut(const int n) : inSinkComponent(n, false) {}

        inline void compute(const std::vector<int>& dist) {
            for (size_t i = 0; i < dist.size(); i++) {
                inSinkComponent[i] = (dist[i] < 0);
            }
        }

        inline std::vector<Vertex> getSourceComponent() const noexcept {
            std::vector<Vertex> component;
            for (size_t i = 0; i < inSinkComponent.size(); i++) {
                if (!inSinkComponent[i]) component.emplace_back(Vertex(i));
            }
            return component;
        }

        inline std::vector<Vertex> getSinkComponent() const noexcept {
            std::vector<Vertex> component;
            for (size_t i = 0; i < inSinkComponent.size(); i++) {
                if (inSinkComponent[i]) component.emplace_back(Vertex(i));
            }
            return component;
        }

        std::vector<bool> inSinkComponent;
    };

public:
    explicit IBFS(const MaxFlowInstance& instance) :
        instance(instance),
        graph(instance.graph),
        n(graph.numVertices()),
        terminal{instance.source, instance.sink},
        residualCapacity(instance.getCurrentCapacities()),
        distance(n, 0),
        maxDistance{0, 0},
        currentEdge(n, noEdge),
        treeData(n),
        orphans{n, n},
        cut(n) {
    }

public:
    inline void run() noexcept {
        initialize<FORWARD>();
        initialize<BACKWARD>();
        runAfterInitialize();
        cut.compute(distance);
    }

    inline std::vector<Vertex> getSourceComponent() const noexcept {
        return cut.getSourceComponent();
    }

    inline std::vector<Vertex> getSinkComponent() const noexcept {
        return cut.getSinkComponent();
    }

    inline const std::vector<bool>& getInSinkComponent() const noexcept {
        return cut.inSinkComponent;
    }

    inline std::vector<Edge> getCutEdges() const noexcept {
        std::vector<Edge> edges;
        for (const Vertex vertex : graph.vertices()) {
            if (cut.inSinkComponent[vertex]) continue;
            for (const Edge edge : graph.edgesFrom(vertex)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (!cut.inSinkComponent[to]) continue;
                edges.emplace_back(edge);
            }
        }
        return edges;
    }

    inline FlowType getFlowValue() noexcept {
        FlowType flow = 0;
        for (const Edge edge : graph.edgesFrom(terminal[BACKWARD])) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            flow += instance.getCapacity(reverseEdge) - residualCapacity[reverseEdge];
        }
        return flow;
    }

    inline const std::vector<FlowType>& getResidualCapacities() const noexcept {
        return residualCapacity;
    }

    inline std::vector<FlowType> getCleanResidualCapacities() const noexcept {
        std::vector<FlowType> result = getResidualCapacities();
        for (size_t i = 0; i < result.size(); i++) {
            if (result[i] > INFTY - 1000) result[i] = INFTY;
        }
        return result;
    }

private:
    template<int DIRECTION>
    inline void initialize() noexcept {
        setDistance<DIRECTION>(terminal[DIRECTION], 1);
        maxDistance[DIRECTION] = 1;
        currentEdge[terminal[DIRECTION]] = graph.beginEdgeFrom(terminal[DIRECTION]);
        nextQ[DIRECTION].push(terminal[DIRECTION]);
    }

    inline void runAfterInitialize() noexcept {
        while (true) {
            if (!grow<FORWARD>()) {
                // Continue the backward search to ensure we get the minimal sink component
                while (grow<BACKWARD>());
                break;
            }
            if (!grow<BACKWARD>()) break;
        }
    }

    template<int DIRECTION>
    inline bool grow() noexcept {
        std::swap(Q[DIRECTION], nextQ[DIRECTION]);
        maxDistance[DIRECTION]++;
        while (!Q[DIRECTION].empty()) {
            const Vertex from = Q[DIRECTION].front();
            Q[DIRECTION].pop();
            if (getDistance<DIRECTION>(from) != maxDistance[DIRECTION] - 1) continue;

            for (Edge edge = graph.beginEdgeFrom(from); edge < graph.endEdgeFrom(from); edge++) {
                const Edge edgeTowardsSink = getForwardEdge<DIRECTION>(edge);
                if (!isEdgeResidual(edgeTowardsSink)) continue;
                const Vertex to = graph.get(ToVertex, edge);
                if (distance[to] == 0) {
                    setDistance<DIRECTION>(to, maxDistance[DIRECTION]);
                    currentEdge[to] = graph.beginEdgeFrom(to);
                    treeData.addVertex(from, to, edgeTowardsSink);
                    nextQ[DIRECTION].push(to);
                } else if (!isVertexInTree<DIRECTION>(to)) {
                    const Vertex sourceEndpoint = (DIRECTION == BACKWARD) ? to : from;
                    const Vertex sinkEndpoint = (DIRECTION == BACKWARD) ? from : to;
                    augment(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
                    if (getDistance<DIRECTION>(from) != maxDistance[DIRECTION] - 1) break;
                    if (isEdgeResidual(edgeTowardsSink)) edge--;
                }
            }
        }
        return !nextQ[DIRECTION].empty();
    }

    inline void augment(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) noexcept {
        const FlowType flow = findBottleneckCapacity(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
        augmentPath(sourceEndpoint, sinkEndpoint, edgeTowardsSink, flow);
        adoptOrphans<FORWARD>();
        adoptOrphans<BACKWARD>();
    }

    inline FlowType findBottleneckCapacity(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) const noexcept {
        FlowType bottleneck = residualCapacity[edgeTowardsSink];
        findBottleneckCapacity<FORWARD>(sourceEndpoint, bottleneck);
        findBottleneckCapacity<BACKWARD>(sinkEndpoint, bottleneck);
        return bottleneck;
    }

    template<int DIRECTION>
    inline void findBottleneckCapacity(Vertex vertex, FlowType& bottleneck) const noexcept {
        while (vertex != terminal[DIRECTION]) {
            const Edge edge = treeData.parentEdge[vertex];
            bottleneck = std::min(bottleneck, residualCapacity[edge]);
            vertex = treeData.parentVertex[vertex];
        }
    }

    inline void augmentPath(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink, const FlowType flow) noexcept {
        const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
        residualCapacity[edgeTowardsSink] -= flow;
        residualCapacity[edgeTowardsSource] += flow;
        augmentPath<FORWARD>(sourceEndpoint, flow);
        augmentPath<BACKWARD>(sinkEndpoint, flow);
    }

    template<int DIRECTION>
    inline void augmentPath(Vertex vertex, const FlowType flow) noexcept {
        while (vertex != terminal[DIRECTION]) {
            const Edge edgeTowardsSink = treeData.parentEdge[vertex];
            const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
            const Vertex parentVertex = treeData.parentVertex[vertex];
            residualCapacity[edgeTowardsSink] -= flow;
            residualCapacity[edgeTowardsSource] += flow;
            if (!isEdgeResidual(edgeTowardsSink))
                makeOrphan<DIRECTION>(vertex);
            vertex = parentVertex;
        }
    }

    template<int DIRECTION>
    inline void makeOrphan(const Vertex vertex) noexcept {
        orphans[DIRECTION].addVertex(vertex, getDistance<DIRECTION>(vertex));
        treeData.removeVertex(vertex);
    }

    template<int DIRECTION>
    inline void adoptOrphans() noexcept {
        while (!orphans[DIRECTION].empty()) {
            const Vertex orphan = orphans[DIRECTION].pop();
            if (adoptWithSameDistance<DIRECTION>(orphan)) continue;
            if (getDistance<DIRECTION>(orphan) == maxDistance[DIRECTION]) {
                distance[orphan] = 0;
                continue;
            }
            treeData.removeChildren(orphan, [&](const Vertex child) {
                orphans[DIRECTION].addVertex(child, getDistance<DIRECTION>(child));
            });
            if (!adoptWithNewDistance<DIRECTION>(orphan))
                distance[orphan] = 0;
        }
    }

    template<int DIRECTION>
    inline bool adoptWithSameDistance(const Vertex orphan) noexcept {
        for (Edge edge = currentEdge[orphan]; edge < graph.endEdgeFrom(orphan); edge++) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (!isEdgeAdmissible<DIRECTION>(orphan, from)) continue;
            treeData.addVertex(from, orphan, edgeTowardsSink);
            currentEdge[orphan] = edge;
            return true;
        }
        return false;
    }

    template<int DIRECTION>
    inline bool adoptWithNewDistance(const Vertex orphan) noexcept {
        const int oldDistance = getDistance<DIRECTION>(orphan);
        int newDistance = maxDistance[DIRECTION];
        Edge newEdge = noEdge;
        Edge newEdgeTowardsSink = noEdge;
        Vertex newParent = noVertex;
        for (const Edge edge : graph.edgesFrom(orphan)) {
            const Edge edgeTowardsSink = getBackwardEdge<DIRECTION>(edge);
            if (!isEdgeResidual(edgeTowardsSink)) continue;
            const Vertex from = graph.get(ToVertex, edge);
            if (!isVertexInTree<DIRECTION>(from)) continue;
            const int fromDistance = getDistance<DIRECTION>(from);
            if (fromDistance < newDistance) {
                newDistance = fromDistance;
                newEdge = edge;
                newEdgeTowardsSink = edgeTowardsSink;
                newParent = from;
                if (newDistance == oldDistance) break;
            }
        }
        if (newEdge == noEdge) return false;
        setDistance<DIRECTION>(orphan, newDistance + 1);
        currentEdge[orphan] = newEdge;
        treeData.addVertex(newParent, orphan, newEdgeTowardsSink);
        if (newDistance + 1 == maxDistance[DIRECTION])
            nextQ[DIRECTION].push(orphan);
        return true;
    }

    template<int DIRECTION>
    inline Edge getForwardEdge(const Edge edge) noexcept {
        if (DIRECTION == BACKWARD)
            return graph.get(ReverseEdge, edge);
        else
            return edge;
    }

    template<int DIRECTION>
    inline Edge getBackwardEdge(const Edge edge) noexcept {
        if (DIRECTION == FORWARD)
            return graph.get(ReverseEdge, edge);
        else
            return edge;
    }

    template<int DIRECTION>
    inline int getDistance(const Vertex vertex) const noexcept {
        return (DIRECTION == BACKWARD) ? -distance[vertex] : distance[vertex];
    }

    template<int DIRECTION>
    inline void setDistance(const Vertex vertex, const int dist) noexcept {
        distance[vertex] = (DIRECTION == BACKWARD) ? -dist : dist;
    }

    template<int DIRECTION>
    inline bool isVertexInTree(const Vertex vertex) noexcept {
        return getDistance<DIRECTION>(vertex) > 0;
    }

    template<int DIRECTION>
    inline bool isEdgeAdmissible(const Vertex from, const Vertex to) const noexcept {
        return getDistance<DIRECTION>(from) == getDistance<DIRECTION>(to) + 1;
    }

    inline bool isEdgeResidual(const Edge edge) const noexcept {
        return pmf::isNumberPositive(residualCapacity[edge]);
    }

    inline void checkDistanceInvariants(const bool allowOrphans = false) const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkDistanceInvariants(vertex, allowOrphans);
        }
    }

    inline void checkDistanceInvariants(const Vertex vertex, const bool allowOrphans = false) const noexcept {
        const Vertex parent = treeData.parentVertex[vertex];
        if (abs(distance[vertex]) <= 1) {
            Ensure(parent == noVertex, "Vertex " << vertex << " with distance <= 1 has parent " << parent << "!");
        } else {
            if (!graph.isVertex(parent)) {
                Ensure(allowOrphans, "Vertex " << vertex << " has an invalid parent!");
            } else {
                Ensure(abs(distance[vertex]) == abs(distance[parent]) + 1, "Distance of " << vertex << " is " << distance[vertex] << ", but distance of " << parent << " is " << distance[parent] << "!");
            }
        }
    }

    inline void checkChildrenRelation() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkChildrenRelation(vertex);
        }
    }

    inline void checkChildrenRelation(const Vertex vertex) const noexcept {
        std::vector<bool> isChild(n, false);
        Vertex child = treeData.firstChild[vertex];
        while (child != noVertex) {
            isChild[child] = true;
            child = treeData.nextSibling[child];
        }
        for (const Vertex c : graph.vertices()) {
            const Vertex parent = treeData.parentVertex[c];
            if (isChild[c]) {
                Ensure(parent == vertex, "Child " << c << " of << " << vertex << " has the wrong parent!");
            } else {
                Ensure(parent != vertex, "Vertex " << c << " has parent << " << vertex << " but is not a child!");
            }
        }
    }

private:
    const MaxFlowInstance& instance;
    const GraphType& graph;
    const int n;
    const Vertex terminal[2];
    std::vector<FlowType> residualCapacity;
    std::vector<int> distance; // 1 for source, -1 for sink, 0 for n-vertices
    int maxDistance[2];
    std::vector<Edge> currentEdge;
    TreeData treeData;
    std::queue<Vertex> Q[2];
    std::queue<Vertex> nextQ[2];
    OrphanBuckets orphans[2];
    Cut cut;
};
