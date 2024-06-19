#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

class ExcessesIBFS {

    struct ExcessBuckets {
        ExcessBuckets(const int n) :
            buckets(n), positionOfVertex(n, -1), maxBucket(-1) {
        }

        inline void addVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex[vertex] != -1) return;
            positionOfVertex[vertex] = buckets[dist].size();
            buckets[dist].emplace_back(vertex);
            maxBucket = std::max(maxBucket, dist);
        }

        inline void removeVertex(const Vertex vertex, const int dist) noexcept {
            if (positionOfVertex[vertex] == -1) return;
            const Vertex other = buckets[dist].back();
            positionOfVertex[other] = positionOfVertex[vertex];
            buckets[dist][positionOfVertex[vertex]] = other;
            buckets[dist].pop_back();
            positionOfVertex[vertex] = -1;
            if (maxBucket == dist && buckets[dist].empty()) {
                while (maxBucket >= 0 && buckets[maxBucket].empty()) maxBucket--;
            }
        }

        inline void increaseBucket(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(newDist > oldDist, "Distance has not increased!");
            Assert(static_cast<size_t>(positionOfVertex[vertex]) < buckets[oldDist].size(), "Vertex is not in bucket!");
            Assert(buckets[oldDist][positionOfVertex[vertex]] == vertex, "Vertex is not in bucket!");
            const Vertex other = buckets[oldDist].back();
            positionOfVertex[other] = positionOfVertex[vertex];
            buckets[oldDist][positionOfVertex[vertex]] = other;
            buckets[oldDist].pop_back();
            positionOfVertex[vertex] = buckets[newDist].size();
            buckets[newDist].emplace_back(vertex);
            maxBucket = std::max(maxBucket, newDist);
        }

        inline bool empty() const noexcept {
            return maxBucket < 0;
        }

        inline Vertex front() const noexcept {
            Assert(!empty(), "Buckets are empty!");
            return buckets[maxBucket].back();
        }

        std::vector<std::vector<Vertex>> buckets;
        std::vector<int> positionOfVertex;
        int maxBucket;
    };

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
            if (parentVertex[vertex] == noVertex) return;
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
    explicit ExcessesIBFS(const MaxFlowInstance& instance) :
        instance(instance),
        graph(instance.graph),
        n(graph.numVertices()),
        terminal{instance.source, instance.sink},
        residualCapacity(instance.currentCapacity),
        distance(n, 0),
        excess(n, 0),
        maxDistance{0, 0},
        currentEdge(n, noEdge),
        treeData(n),
        excessVertices{ExcessBuckets(n), ExcessBuckets(n)},
        cut(n) {
    }

public:
    inline void run() noexcept {
        initialize<FORWARD>();
        initialize<BACKWARD>();
        runAfterInitialize();
        cut.compute(distance);
    }

    //TODO: Implement update after capacity change

    inline std::vector<Vertex> getSourceComponent() const noexcept {
        return cut.getSourceComponent();
    }

    inline std::vector<Vertex> getSinkComponent() const noexcept {
        return cut.getSinkComponent();
    }

    //TODO: Maintain the flow value throughout the algorithm.
    inline int getFlowValue() const noexcept {
        int flow = 0;
        for (const auto [edge, from] : graph.edgesWithFromVertex()) {
            const Vertex to = graph.get(ToVertex, edge);
            if (distance[from] >= 0 && distance[to] < 0) flow += instance.getCapacity(edge) - residualCapacity[edge];
        }
        return flow;
    }

private:
    template<int DIRECTION>
    inline void initialize() noexcept {
        setDistance<DIRECTION>(terminal[DIRECTION], 1);
        addExcess<DIRECTION>(terminal[DIRECTION], -INFTY);
        maxDistance[DIRECTION] = 1;
        currentEdge[terminal[DIRECTION]] = graph.beginEdgeFrom(terminal[DIRECTION]);
        nextQ[DIRECTION].push(terminal[DIRECTION]);
    }

    inline void runAfterInitialize() noexcept {
        while (true) {
            //TODO: Clever alternation
            if (!grow<FORWARD>()) break;
            if (!grow<BACKWARD>()) break;
        }
    }

    template<int DIRECTION>
    inline bool grow() {
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
        const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
        const int flow = findBottleneckCapacity(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
        pushFlow<BACKWARD>(sourceEndpoint, sinkEndpoint, edgeTowardsSink, edgeTowardsSource, flow);
        drainExcesses<FORWARD>(sourceEndpoint);
        drainExcesses<BACKWARD>(sinkEndpoint);
    }

    inline int findBottleneckCapacity(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) const noexcept {
        auto[sourceBottleneck, sourceRoot] = findBottleneckCapacity(sourceEndpoint);
        auto[sinkBottleneck, sinkRoot] = findBottleneckCapacity(sinkEndpoint);
        if (sourceRoot == terminal[FORWARD] && sinkRoot == terminal[BACKWARD]) return residualCapacity[edgeTowardsSink];
        else if (sourceRoot == terminal[FORWARD]) return std::min(sourceBottleneck, residualCapacity[edgeTowardsSink]);
        else if (sinkRoot == terminal[BACKWARD]) return std::min(sinkBottleneck, residualCapacity[edgeTowardsSink]);
        return std::min({excess[sourceRoot], sourceBottleneck, residualCapacity[edgeTowardsSink], sinkBottleneck, -excess[sinkRoot]});
    }

    inline std::pair<int, Vertex> findBottleneckCapacity(const Vertex start) const noexcept {
        int bottleneck = INFTY;
        Vertex vertex = start;
        while (treeData.parentVertex[vertex] != noVertex) {
            const Edge edge = treeData.parentEdge[vertex];
            bottleneck = std::min(bottleneck, residualCapacity[edge]);
            vertex = treeData.parentVertex[vertex];
        }
        return {bottleneck, vertex};
    }

    template<int DIRECTION>
    inline void drainExcesses(const Vertex start) noexcept {
        const int exc = getExcess<DIRECTION>(start);
        if (exc < 0) return;
        else if (exc == 0) {
            makeOrphan<DIRECTION>(start);
            adoptOrphans<DIRECTION>();
        } else {
            excessVertices[DIRECTION].addVertex(start, getDistance<DIRECTION>(start));
            while (!excessVertices[DIRECTION].empty()) {
                const Vertex vertex = excessVertices[DIRECTION].front();
                drainExcess<DIRECTION>(vertex);
                adoptOrphans<DIRECTION>();
            }
        }
    }

    template<int DIRECTION>
    inline void drainExcess(Vertex vertex) noexcept {
        while (treeData.parentVertex[vertex] != noVertex) {
            const Vertex parentVertex = treeData.parentVertex[vertex];
            const Edge edgeTowardsSink = treeData.parentEdge[vertex];
            Assert(isEdgeResidual(edgeTowardsSink), "Tree edge is not residual!");
            const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
            const int exc = getExcess<DIRECTION>(vertex);
            const int res = residualCapacity[edgeTowardsSink];
            const int flow = std::min(res, exc);
            pushFlow<DIRECTION>(vertex, parentVertex, edgeTowardsSink, edgeTowardsSource, flow);
            if (flow == res) {
                makeOrphan<DIRECTION>(vertex);
            }
            if (flow == exc) {
                excessVertices[DIRECTION].removeVertex(vertex, getDistance<DIRECTION>(vertex));
            }
            vertex = parentVertex;
            if (getExcess<DIRECTION>(vertex) > 0) {
                excessVertices[DIRECTION].addVertex(vertex, getDistance<DIRECTION>(vertex));
            }
            else Assert(treeData.parentVertex[vertex] == noVertex, "Non-root vertex has zero excess!");
        }

        if (getExcess<DIRECTION>(vertex) >= 0)
            makeOrphan<DIRECTION>(vertex);
    }
    template<int DIRECTION>
    inline void makeOrphan(const Vertex vertex) noexcept {
        orphans[DIRECTION].push(vertex);
        treeData.removeVertex(vertex);
    }

    //TODO: Adopt orphans in increasing order of distance
    //TODO: Three-pass/hybrid adoption
    template<int DIRECTION>
    inline void adoptOrphans() noexcept {
        while (!orphans[DIRECTION].empty()) {
            const Vertex orphan = orphans[DIRECTION].front();
            orphans[DIRECTION].pop();
            if (adoptWithSameDistance<DIRECTION>(orphan)) continue;
            if (getDistance<DIRECTION>(orphan) == maxDistance[DIRECTION]) {
                removeOrphan<DIRECTION>(orphan);
                continue;
            }
            treeData.removeChildren(orphan, [&](const Vertex child) {
                orphans[DIRECTION].push(child);
            });
            if (!adoptWithNewDistance<DIRECTION>(orphan)) {
                removeOrphan<DIRECTION>(orphan);
            }
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
        if (getExcess<DIRECTION>(orphan) > 0)
            excessVertices[DIRECTION].increaseBucket(orphan, oldDistance, newDistance + 1);
        currentEdge[orphan] = newEdge;
        treeData.addVertex(newParent, orphan, newEdgeTowardsSink);
        if (newDistance + 1 == maxDistance[DIRECTION])
            nextQ[DIRECTION].push(orphan);
        return true;
    }

    template<int DIRECTION>
    inline void removeOrphan(const Vertex orphan) noexcept {
        if (getExcess<DIRECTION>(orphan) > 0) {
            excessVertices[DIRECTION].removeVertex(orphan, getDistance<DIRECTION>(orphan));
            setDistance<!DIRECTION>(orphan, maxDistance[!DIRECTION]);
            nextQ[!DIRECTION].push(orphan);
            currentEdge[orphan] = graph.beginEdgeFrom(orphan);
        } else {
            distance[orphan] = 0;
        }
    }

    template<int DIRECTION>
    inline void pushFlow(const Vertex from, const Vertex to, const Edge edge, const Edge reverseEdge, const int flow) noexcept {
        addExcess<DIRECTION>(from, -flow);
        addExcess<DIRECTION>(to, flow);
        residualCapacity[edge] -= flow;
        residualCapacity[reverseEdge] += flow;
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
    inline int getExcess(const Vertex vertex) const noexcept {
        return (DIRECTION == FORWARD) ? -excess[vertex] : excess[vertex];
    }

    template<int DIRECTION>
    inline void addExcess(const Vertex vertex, const int add) noexcept {
        if (DIRECTION == FORWARD)
            excess[vertex] -= add;
        else
            excess[vertex] += add;
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
        return residualCapacity[edge] > 0;
    }

    template<int DIRECTION>
    inline void checkBuckets() const noexcept {
        for (size_t i = 0; static_cast<int>(i) <= excessVertices[DIRECTION].maxBucket; i++) {
            for (const Vertex vertex : excessVertices[DIRECTION].buckets[i]) {
                Assert(getExcess<DIRECTION>(vertex) > 0, "Vertex in bucket has no excess!");
                Assert(getDistance<DIRECTION>(vertex) == int(i), "Vertex is in wrong bucket!");
            }
        }
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

    inline void checkTreeResidual() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            if (treeData.parentEdge[vertex] == noEdge) continue;
            Assert(isEdgeResidual(treeData.parentEdge[vertex]), "Tree edge is not residual!");
        }
    }

    inline void checkFlowConservation() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkFlowConservation(vertex);
        }
    }

    inline int getInflow(const Vertex vertex) const noexcept {
        int inflow = 0;
        for (const Edge edge : graph.edgesFrom(vertex)) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            inflow += instance.getCapacity(reverseEdge) - residualCapacity[reverseEdge];
        }
        return inflow;
    }

    inline void checkFlowConservation(const Vertex vertex) const noexcept {
        if (vertex == terminal[FORWARD] || vertex == terminal[BACKWARD]) return;
        Assert(getInflow(vertex) == excess[vertex], "Flow conservation not fulfilled!");
    }

private:
    const MaxFlowInstance& instance;
    const StaticFlowGraph& graph;
    const int n;
    const Vertex terminal[2];
    std::vector<int> residualCapacity;
    // positive for s-vertices, negative for t-vertices, 0 for n-vertices
    std::vector<int> distance;
    // positive for source roots, negative for sink roots
    // non-positive for other s-vertices, non-negative for other t-vertices
    std::vector<int> excess;
    int maxDistance[2];
    std::vector<Edge> currentEdge; //TODO: Can be represented implicitly, but unclear if that saves time
    TreeData treeData;
    std::queue<Vertex> Q[2];
    std::queue<Vertex> nextQ[2];
    ExcessBuckets excessVertices[2];
    std::queue<Vertex> orphans[2];
    Cut cut;
};
