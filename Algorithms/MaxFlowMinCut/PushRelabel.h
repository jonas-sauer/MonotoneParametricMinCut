#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

class PushRelabel {

private:
    inline static constexpr int VertexToEdgeRatio = 12;

public:
    struct VertexBuckets {
        VertexBuckets(const int n) :
            activeVertices(n), inactiveVertices(n), isVertexActive(n, false), positionOfVertex(Vector::id<int>(n)), maxActiveBucket(-1), maxBucket(-1) {
        }

        inline void initialize(const int sink) {
            std::vector<std::vector<Vertex>>(activeVertices.size()).swap(activeVertices);
            std::vector<std::vector<Vertex>>(inactiveVertices.size()).swap(inactiveVertices);
            std::vector<bool>(isVertexActive.size(), false).swap(isVertexActive);
            maxActiveBucket = -1;
            maxBucket = 0;

            for (int i = 0; static_cast<size_t>(i) < isVertexActive.size(); i++) {
                if (i == sink) continue;
                positionOfVertex[i] = inactiveVertices[0].size();
                inactiveVertices[0].emplace_back(i);
            }
        }

        inline void rebuild(const std::vector<int>& distances, const int sink) {
            std::vector<std::vector<Vertex>>(activeVertices.size()).swap(activeVertices);
            std::vector<std::vector<Vertex>>(inactiveVertices.size()).swap(inactiveVertices);
            maxActiveBucket = -1;
            maxBucket = -1;

            for (int i = 0; static_cast<size_t>(i) < isVertexActive.size(); i++) {
                if (i == sink || static_cast<size_t>(distances[i]) == activeVertices.size()) continue;
                maxBucket = std::max(maxBucket, distances[i]);
                if (isVertexActive[i]) {
                    positionOfVertex[i] = activeVertices[distances[i]].size();
                    activeVertices[distances[i]].emplace_back(i);
                    maxActiveBucket = std::max(maxActiveBucket, distances[i]);
                } else {
                    positionOfVertex[i] = inactiveVertices[distances[i]].size();
                    inactiveVertices[distances[i]].emplace_back(i);
                }
            }
        }

        inline void activateVertex(const Vertex vertex, const int dist) noexcept {
            Assert(checkInvariant(vertex, dist), "Vertex position is invalid");
            if (isVertexActive[vertex]) return;
            isVertexActive[vertex] = true;
            deleteInactiveVertex(vertex, dist);
            positionOfVertex[vertex] = activeVertices[dist].size();
            activeVertices[dist].emplace_back(vertex);
            maxActiveBucket = std::max(maxActiveBucket, dist);
        }

        inline void updateDistance(const Vertex vertex, const int oldDist, const int newDist) noexcept {
            Assert(checkInvariant(vertex, oldDist), "Vertex position is invalid");
            Assert(!isVertexActive[vertex], "Cannot update distance of active vertex!");
            deleteInactiveVertex(vertex, oldDist);
            if (static_cast<size_t>(newDist) < inactiveVertices.size()) {
                positionOfVertex[vertex] = inactiveVertices[newDist].size();
                inactiveVertices[newDist].emplace_back(vertex);
                maxBucket = std::max(maxBucket, newDist);
            }
        }

        inline void deleteInactiveVertex(const Vertex vertex, const int dist) noexcept {
            const int pos = positionOfVertex[vertex];
            const Vertex other = inactiveVertices[dist].back();
            inactiveVertices[dist][pos] = other;
            inactiveVertices[dist].pop_back();
            positionOfVertex[other] = pos;
        }

        template<typename FUNCTION>
        inline void pruneGap(const int dist, const FUNCTION& callback) noexcept {
            if (!isBucketEmpty(dist)) return;
            for (int d = dist + 1; d <= maxBucket; d++) {
                for (const Vertex vertex : inactiveVertices[d]) {
                    callback(vertex);
                }
                inactiveVertices[d].clear();
            }
            maxBucket = dist;
        }

        inline Vertex pop() noexcept {
            Assert(!empty(), "All buckets are empty!");
            const Vertex vertex = activeVertices[maxActiveBucket].back();
            activeVertices[maxActiveBucket].pop_back();
            isVertexActive[vertex] = false;
            positionOfVertex[vertex] = inactiveVertices[maxActiveBucket].size();
            inactiveVertices[maxActiveBucket].emplace_back(vertex);
            while (maxActiveBucket >= 0 && activeVertices[maxActiveBucket].empty()) maxActiveBucket--;
            return vertex;
        }

        inline bool empty() const noexcept {
            return maxActiveBucket < 0;
        }

        inline bool isBucketEmpty(const int dist) const noexcept {
            return activeVertices[dist].empty() && inactiveVertices[dist].empty();
        }

        inline bool checkInvariant(const Vertex vertex, const int dist) const noexcept {
            if (static_cast<size_t>(dist) == activeVertices.size()) return true;
            if (isVertexActive[vertex]) {
                if (activeVertices[dist].size() <= static_cast<size_t>(positionOfVertex[vertex])) return false;
                return activeVertices[dist][positionOfVertex[vertex]] == vertex;
            } else {
                if (inactiveVertices[dist].size() <= static_cast<size_t>(positionOfVertex[vertex])) return false;
                return inactiveVertices[dist][positionOfVertex[vertex]] == vertex;
            }
        }

        std::vector<std::vector<Vertex>> activeVertices;
        std::vector<std::vector<Vertex>> inactiveVertices;
        std::vector<bool> isVertexActive;
        std::vector<int> positionOfVertex;
        int maxActiveBucket;
        int maxBucket;
    };

    struct Cut {
        Cut(const int n) : inSinkComponent(n, false) {}

        inline void compute(const std::vector<int>& dist) {
            for (size_t i = 0; i < dist.size(); i++) {
                inSinkComponent[i] = static_cast<size_t>(dist[i]) < inSinkComponent.size();
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
    explicit PushRelabel(const StaticFlowGraph& graph) :
        graph(graph),
        n(graph.numVertices()),
        sourceVertex(noVertex),
        sinkVertex(noVertex),
        residualCapacity(graph.get(Capacity)),
        distance(n, 0),
        excess(n, 0),
        currentEdge(n, noEdge),
        vertexBuckets(n),
        workSinceLastUpdate(0),
        workLimit(VertexToEdgeRatio * graph.numVertices() + graph.numEdges()),
        cut(n) {
    }

public:
    inline void run(const Vertex source, const Vertex sink) noexcept {
        clear();
        sourceVertex = source;
        sinkVertex = sink;
        initialize();
        run();
    }

    /*inline void runWithCapacityUpdate(const std::vector<int>& sourceCapacities, const std::vector<int>& sinkCapacities) noexcept {
        updateCapacities(sourceCapacities, sinkCapacities);
        run();
    }*/

    inline std::vector<Vertex> getSourceComponent() const noexcept {
        return cut.getSourceComponent();
    }

    inline std::vector<Vertex> getSinkComponent() const noexcept {
        return cut.getSinkComponent();
    }

    inline int getFlowValue() noexcept {
        int flow = 0;
        for (const Edge edge : graph.edgesFrom(sinkVertex)) {
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
        Vector::fill(excess, 0);
        for (const Vertex vertex : graph.vertices()) {
            currentEdge[vertex] = graph.beginEdgeFrom(vertex);
        }
    }

    inline void initialize() noexcept {
        vertexBuckets.initialize(sinkVertex);
        distance[sourceVertex] = distance.size();
        for (const Edge edge : graph.edgesFrom(sourceVertex)) {
            const int capacity = graph.get(Capacity, edge);
            if (capacity == 0) continue;
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            residualCapacity[edge] = 0;
            residualCapacity[reverseEdge] += capacity;
            const Vertex to = graph.get(ToVertex, edge);
            excess[to] = capacity;
            makeVertexActive(to);
        }
    }

    //TODO: Cannot change the graph, so handle this differently
    /*inline void updateCapacities(const std::vector<int>& sourceCapacities, const std::vector<int>& sinkCapacities) noexcept {
        Edge edgeFromSource = graph.beginEdgeFrom(sourceVertex);
        for (size_t i = 0; i < sourceCapacities.size(); i++, edgeFromSource++) {
            const int oldCapacity = graph.get(Capacity, edgeFromSource);
            graph.set(Capacity, edgeFromSource, sourceCapacities[i]);
            const Vertex to = graph.get(ToVertex, edgeFromSource);
            if (distance[to] >= n) continue;
            Assert(residualCapacity[edgeFromSource] == 0, "Source-incident cut edge is not saturated!");

            const int newResidualCapacity = sourceCapacities[i] - oldCapacity;
            Assert(newResidualCapacity >= 0, "Capacity of source-incident edge has decreased!");
            if (newResidualCapacity > 0) {
                const Edge edgeToSource = graph.get(ReverseEdge, edgeFromSource);
                pushFlow(sourceVertex, to, edgeFromSource, edgeToSource, newResidualCapacity);
            }
        }

        Edge edgeFromSink = graph.beginEdgeFrom(sinkVertex);
        for (size_t i = 0; i < sinkCapacities.size(); i++, edgeFromSink++) {
            const Edge edgeToSink = graph.get(ReverseEdge, edgeFromSink);
            const int oldCapacity = graph.get(Capacity, edgeToSink);
            graph.set(Capacity, edgeToSink, sinkCapacities[i]);

            const int lostCapacity = oldCapacity - sinkCapacities[i];
            Assert(lostCapacity >= 0, "Capacity of sink-incident edge has increased!");
            const int newResidualCapacity = residualCapacity[edgeToSink] - lostCapacity;
            if (newResidualCapacity >= 0) {
                residualCapacity[edgeToSink] = newResidualCapacity;
            } else {
                const Vertex from = graph.get(ToVertex, edgeFromSink);
                pushFlow(sinkVertex, from, edgeFromSink, edgeToSink, -newResidualCapacity);
            }
        }
    }*/

    inline void run() noexcept {
        workSinceLastUpdate = 0;
        while (!vertexBuckets.empty()) {
            const Vertex vertex = vertexBuckets.pop();
            discharge(vertex);
            if (distance[vertex] < n && excess[vertex] > 0)
                makeVertexActive(vertex);

            if (workSinceLastUpdate > workLimit) {
                globalDistanceUpdate();
                workSinceLastUpdate = 0;
            }
        }
        cut.compute(distance);
    }

    inline void discharge(const Vertex vertex) noexcept {
       while (excess[vertex] > 0) {
            const Edge edge = findPushableEdge(vertex);
            if (edge == noEdge) {
                relabel(vertex);
                break;
            }
            const Vertex to = graph.get(ToVertex, edge);
            pushFlow(vertex, to, edge);
        }
    }

    inline Edge findPushableEdge(const Vertex vertex) noexcept {
        for (Edge edge = currentEdge[vertex]; edge < graph.endEdgeFrom(vertex); edge++) {
            const Vertex to = graph.get(ToVertex, edge);
            if (isEdgeResidual(edge) && isEdgeAdmissible(vertex, to)) {
                currentEdge[vertex] = edge;
                return edge;
            }
        }
        return noEdge;
    }

    inline void pushFlow(const Vertex from, const Vertex to, const Edge edge) noexcept {
        const int flow = std::min(excess[from], residualCapacity[edge]);
        const Edge reverseEdge = graph.get(ReverseEdge, edge);
        pushFlow(from, to, edge, reverseEdge, flow);
    }

    inline void pushFlow(const Vertex from, const Vertex to, const Edge edge, const Edge reverseEdge, const int flow) noexcept {
        residualCapacity[edge] -= flow;
        residualCapacity[reverseEdge] += flow;
        excess[from] -= flow;
        excess[to] += flow;
        makeVertexActive(to);
    }

    inline void relabel(const Vertex vertex) noexcept {
        const int oldDistance = distance[vertex];
        int newDistance = n;
        Edge newAdmissibleEdge = noEdge;
        workSinceLastUpdate += VertexToEdgeRatio + graph.outDegree(vertex);
        for (const Edge edge : graph.edgesFrom(vertex)) {
            if (!isEdgeResidual(edge)) continue;
            const Vertex to = graph.get(ToVertex, edge);
            if (distance[to] + 1 < newDistance) {
                newDistance = distance[to] + 1;
                newAdmissibleEdge = edge;
            }
        }
        Assert(newDistance > distance[vertex], "Relabel did not increase the distance!");
        vertexBuckets.updateDistance(vertex, distance[vertex], newDistance);
        distance[vertex] = newDistance;
        currentEdge[vertex] = newAdmissibleEdge;
        vertexBuckets.pruneGap(oldDistance, [&](const Vertex v){
            distance[v] = n;
        });
    }

    inline bool isEdgeAdmissible(const Vertex from, const Vertex to) const noexcept {
        return distance[from] == distance[to] + 1;
    }

    inline bool isEdgeResidual(const Edge edge) const noexcept {
        return residualCapacity[edge] > 0;
    }

    inline void makeVertexActive(const Vertex vertex) noexcept {
        if (vertex == sinkVertex) return;
        Assert(distance[vertex] < n, "Distance label too high!");
        Assert(excess[vertex] > 0, "Vertex does not have excess!");
        vertexBuckets.activateVertex(vertex, distance[vertex]);
    }

    inline void globalDistanceUpdate() noexcept {
        Vector::fill(distance, n);
        std::queue<Vertex> Q;

        distance[sinkVertex] = 0;
        Q.push(sinkVertex);

        while (!Q.empty()) {
            const Vertex u = Q.front();
            Q.pop();
            for (const Edge e : graph.edgesFrom(u)) {
                const Edge re = graph.get(ReverseEdge, e);
                if (!isEdgeResidual(re)) continue;
                const Vertex v = graph.get(ToVertex, e);
                if (distance[v] < n) continue;
                currentEdge[v] = graph.beginEdgeFrom(v);
                distance[v] = std::min(distance[v], distance[u] + 1);
                Q.push(v);
            }
        }

        vertexBuckets.rebuild(distance, sinkVertex);
    }

    inline bool checkBucketInvariants() noexcept {
        for (const Vertex vertex : graph.vertices()) {
            if (vertex == sourceVertex || vertex == sinkVertex) continue;
            if (!vertexBuckets.checkInvariant(vertex, distance[vertex])) return false;
        }
        return true;
    }

    inline bool checkCurrentEdgeInvariant(const Vertex vertex) const noexcept {
        for (Edge edge = graph.beginEdgeFrom(vertex); edge < currentEdge[vertex]; edge++) {
            const Vertex to = graph.get(ToVertex, edge);
            if (isEdgeResidual(edge) && isEdgeAdmissible(vertex, to)) return false;
        }
        return true;
    }

private:
    const StaticFlowGraph& graph;
    const int n;
    Vertex sourceVertex;
    Vertex sinkVertex;
    std::vector<int> residualCapacity;
    std::vector<int> distance;
    std::vector<int> excess;
    std::vector<Edge> currentEdge;
    VertexBuckets vertexBuckets;
    int workSinceLastUpdate;
    const int workLimit;
    Cut cut;
};
