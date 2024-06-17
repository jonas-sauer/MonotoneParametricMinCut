#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

class IBFS {

    struct TreeData {
        TreeData(const size_t n) :
            parentEdge(n, noEdge),
            parentVertex(n , noVertex),
            firstChild(n, noVertex),
            nextSibling(n, noVertex),
            prevSibling(n, noVertex) {
        }

        inline void clear() noexcept {
            Vector::fill(parentEdge, noEdge);
            Vector::fill(parentVertex, noVertex);
            Vector::fill(firstChild, noVertex);
            Vector::fill(nextSibling, noVertex);
            Vector::fill(prevSibling, noVertex);
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
    explicit IBFS(const StaticFlowGraph& graph) :
        graph(graph),
        n(graph.numVertices()),
        terminal{noVertex, noVertex},
        residualCapacity(graph.get(Capacity)),
        distance(n, 0),
        maxDistance{0, 0},
        currentEdge(n, noEdge),
        treeData(n),
        cut(n) {
    }

public:
    inline void run(const Vertex source, const Vertex sink) noexcept {
        clear();
        terminal[FORWARD] = source;
        terminal[BACKWARD] = sink;
        initialize<FORWARD>();
        initialize<BACKWARD>();
        run();
        cut.compute(distance);
    }

    inline std::vector<Vertex> getSourceComponent() const noexcept {
        return cut.getSourceComponent();
    }

    inline std::vector<Vertex> getSinkComponent() const noexcept {
        return cut.getSinkComponent();
    }

    inline int getFlowValue() noexcept {
        int flow = 0;
        for (const Edge edge : graph.edgesFrom(terminal[BACKWARD])) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            flow += graph.get(Capacity, reverseEdge) - residualCapacity[reverseEdge];
        }
        return flow;
    }

private:
    inline void clear() noexcept {
        //TODO: Too slow?
        residualCapacity = graph.get(Capacity);
        Vector::fill(distance, 0);
        treeData.clear();
        clearDirectional<FORWARD>();
        clearDirectional<BACKWARD>();
    }

    template<int DIRECTION>
    inline void clearDirectional() noexcept {
        maxDistance[DIRECTION] = 0;
        clearQueue(Q[DIRECTION]);
        clearQueue(nextQ[DIRECTION]);
        clearQueue(orphans[DIRECTION]);
    }

    inline static void clearQueue(std::queue<Vertex>& queue) noexcept {
        std::queue<Vertex> emptyQ;
        std::swap(emptyQ, queue);
    }

    template<int DIRECTION>
    inline void initialize() noexcept {
        setDistance<DIRECTION>(terminal[DIRECTION], 1);
        maxDistance[DIRECTION] = 1;
        currentEdge[terminal[DIRECTION]] = graph.beginEdgeFrom(terminal[DIRECTION]);
        nextQ[DIRECTION].push(terminal[DIRECTION]);
    }

    inline void run() noexcept {
        while (true) {
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
        const int flow = findBottleneckCapacity(sourceEndpoint, sinkEndpoint, edgeTowardsSink);
        augmentPath(sourceEndpoint, sinkEndpoint, edgeTowardsSink, flow);
        adoptOrphans<FORWARD>();
        adoptOrphans<BACKWARD>();
    }

    inline int findBottleneckCapacity(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink) const noexcept {
        int bottleneck = residualCapacity[edgeTowardsSink];
        findBottleneckCapacity<FORWARD>(sourceEndpoint, bottleneck);
        findBottleneckCapacity<BACKWARD>(sinkEndpoint, bottleneck);
        return bottleneck;
    }

    template<int DIRECTION>
    inline void findBottleneckCapacity(Vertex vertex, int& bottleneck) const noexcept {
        while (vertex != terminal[DIRECTION]) {
            const Edge edge = treeData.parentEdge[vertex];
            bottleneck = std::min(bottleneck, residualCapacity[edge]);
            vertex = treeData.parentVertex[vertex];
        }
    }

    inline void augmentPath(const Vertex sourceEndpoint, const Vertex sinkEndpoint, const Edge edgeTowardsSink, const int flow) noexcept {
        const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
        residualCapacity[edgeTowardsSink] -= flow;
        residualCapacity[edgeTowardsSource] += flow;
        augmentPath<FORWARD>(sourceEndpoint, flow);
        augmentPath<BACKWARD>(sinkEndpoint, flow);
    }

    template<int DIRECTION>
    inline void augmentPath(Vertex vertex, const int flow) noexcept {
        while (vertex != terminal[DIRECTION]) {
            const Edge edgeTowardsSink = treeData.parentEdge[vertex];
            const Edge edgeTowardsSource = graph.get(ReverseEdge, edgeTowardsSink);
            const Vertex parentVertex = treeData.parentVertex[vertex];
            if (residualCapacity[edgeTowardsSink] == flow)
                makeOrphan<DIRECTION>(vertex);
            residualCapacity[edgeTowardsSink] -= flow;
            residualCapacity[edgeTowardsSource] += flow;
            vertex = parentVertex;
        }
    }

    template<int DIRECTION>
    inline void makeOrphan(const Vertex vertex) noexcept {
        orphans[DIRECTION].push(vertex);
        treeData.removeVertex(vertex);
    }

    //TODO: Adopt orphans in increasing order of distance
    template<int DIRECTION>
    inline void adoptOrphans() noexcept {
        while (!orphans[DIRECTION].empty()) {
            const Vertex orphan = orphans[DIRECTION].front();
            orphans[DIRECTION].pop();
            if (adoptWithSameDistance<DIRECTION>(orphan)) continue;
            //TODO: If getDistance<DIRECTION>(orphan) == maxDistance[DIRECTION], then we can give up
            treeData.removeChildren(orphan, [&](const Vertex child) {
                orphans[DIRECTION].push(child);
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
        return residualCapacity[edge] > 0;
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
    const StaticFlowGraph& graph;
    const int n;
    Vertex terminal[2];
    std::vector<int> residualCapacity;
    std::vector<int> distance; // 1 for source, -1 for target, 0 for n-vertices
    int maxDistance[2];
    std::vector<Edge> currentEdge; //TODO: Can be represented implicitly, but unclear if that saves time
    TreeData treeData;
    std::queue<Vertex> Q[2];
    std::queue<Vertex> nextQ[2];
    std::queue<Vertex> orphans[2];
    Cut cut;
};
