#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<typename MAX_FLOW_INSTANCE>
class PushRelabel {

public:
    using MaxFlowInstance = MAX_FLOW_INSTANCE;
    using FlowType = MaxFlowInstance::FlowType;
    using GraphType = MaxFlowInstance::GraphType;

private:
    inline static constexpr int VertexToEdgeRatio = 12;

public:
    struct VertexBuckets {
        VertexBuckets(const int n) :
            activeVertices(n), inactiveVertices(n), isVertexActive(n, false), positionOfVertex(Vector::id<int>(n)), maxActiveBucket(-1), maxBucket(-1) {
        }

        inline void initialize(const int sink) {
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
    explicit PushRelabel(const MaxFlowInstance& instance) :
        instance(instance),
        graph(instance.graph),
        n(graph.numVertices()),
        sourceVertex(instance.source),
        sinkVertex(instance.sink),
        residualCapacity(instance.getCurrentCapacities()),
        distance(n, 0),
        excess(n, 0),
        currentEdge(n, noEdge),
        vertexBuckets(n),
        workSinceLastUpdate(0),
        workLimit(VertexToEdgeRatio * graph.numVertices() + graph.numEdges()),
        cut(n) {
        for (const Vertex vertex : graph.vertices()) {
            currentEdge[vertex] = graph.beginEdgeFrom(vertex);
        }
    }

public:
    inline void run() noexcept {
        initialize();
        runAfterInitialize();
    }

    inline void continueAfterUpdate() noexcept {
        updateCapacities();
        runAfterInitialize();
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
        for (const Edge edge : graph.edgesFrom(sinkVertex)) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            flow += instance.getCapacity(reverseEdge) - residualCapacity[reverseEdge];
        }
        return flow;
    }

private:
    inline void initialize() noexcept {
        vertexBuckets.initialize(sinkVertex);
        distance[sourceVertex] = distance.size();
        for (const Edge edge : graph.edgesFrom(sourceVertex)) {
            const FlowType capacity = instance.getCapacity(edge);
            if (capacity == 0) continue;
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            residualCapacity[edge] = 0;
            residualCapacity[reverseEdge] += capacity;
            const Vertex to = graph.get(ToVertex, edge);
            excess[to] = capacity;
            makeVertexActive(to);
        }
    }

    inline void runAfterInitialize() noexcept {
        workSinceLastUpdate = 0;
        while (!vertexBuckets.empty()) {
            const Vertex vertex = vertexBuckets.pop();
            discharge(vertex);
            if (distance[vertex] < n && pmf::isNumberPositive(excess[vertex]))
                makeVertexActive(vertex);

            if (workSinceLastUpdate > workLimit) {
                globalDistanceUpdate();
                workSinceLastUpdate = 0;
            }
        }
        computeCut();
    }

    //TODO: Can this be done faster?
    inline void computeCut() noexcept {
        std::vector<bool>(n, false).swap(cut.inSinkComponent);
        std::queue<Vertex> queue;
        queue.push(sinkVertex);
        cut.inSinkComponent[sinkVertex] = true;
        while (!queue.empty()) {
            const Vertex vertex = queue.front();
            queue.pop();
            for (const Edge edge : graph.edgesFrom(vertex)) {
                const Vertex to = graph.get(ToVertex, edge);
                const Edge edgeToSink = graph.get(ReverseEdge, edge);
                if (!isEdgeResidual(edgeToSink)) continue;
                if (!cut.inSinkComponent[to]) {
                    queue.push(to);
                    cut.inSinkComponent[to] = true;
                }
            }
        }

        Assert(!cut.inSinkComponent[sourceVertex], "No cut found!");
    }

    inline void updateCapacities() noexcept {
        Edge edgeFromSource = graph.beginEdgeFrom(sourceVertex);
        const std::vector<FlowType>& sourceDiff = instance.getSourceDiff();
        for (size_t i = 0; i < sourceDiff.size(); i++, edgeFromSource++) {
            Assert(sourceDiff[i] >= 0, "Capacity of source-incident edge has decreased!");
            residualCapacity[edgeFromSource] += sourceDiff[i];
            const Vertex to = graph.get(ToVertex, edgeFromSource);
            if (distance[to] < n && isEdgeResidual(edgeFromSource)) {
                const FlowType add = residualCapacity[edgeFromSource];
                const Edge edgeToSource = graph.get(ReverseEdge, edgeFromSource);
                residualCapacity[edgeFromSource] = 0;
                residualCapacity[edgeToSource] += add;
                excess[to] += add;
                makeVertexActive(to);
            }
        }

        Edge edgeFromSink = graph.beginEdgeFrom(sinkVertex);
        const std::vector<FlowType>& sinkDiff = instance.getSinkDiff();
        for (size_t i = 0; i < sinkDiff.size(); i++, edgeFromSink++) {
            Assert(sinkDiff[i] <= 0, "Capacity of sink-incident edge has increased!");
            const Edge edgeToSink = graph.get(ReverseEdge, edgeFromSink);
            residualCapacity[edgeToSink] += sinkDiff[i];
            if (pmf::isNumberNegative(residualCapacity[edgeToSink])) {
                const FlowType add = -residualCapacity[edgeToSink];
                const Vertex from = graph.get(ToVertex, edgeFromSink);
                residualCapacity[edgeFromSink] -= add;
                residualCapacity[edgeToSink] = 0;
                excess[from] += add;
                if (distance[from] < n) makeVertexActive(from);
            }
        }
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
        const FlowType flow = std::min(excess[from], residualCapacity[edge]);
        const Edge reverseEdge = graph.get(ReverseEdge, edge);
        pushFlow(from, to, edge, reverseEdge, flow);
    }

    inline void pushFlow(const Vertex from, const Vertex to, const Edge edge, const Edge reverseEdge, const FlowType flow) noexcept {
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
        return pmf::isNumberPositive(residualCapacity[edge]);
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

    inline void checkFlowConservation() const noexcept {
        for (const Vertex vertex : graph.vertices()) {
            checkFlowConservation(vertex);
        }
    }

    inline FlowType getInflow(const Vertex vertex) const noexcept {
        FlowType inflow = 0;
        for (const Edge edge : graph.edgesFrom(vertex)) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            inflow += instance.getCapacity(reverseEdge) - residualCapacity[reverseEdge];
        }
        return inflow;
    }

    inline void checkFlowConservation(const Vertex vertex) const noexcept {
        if (vertex == sourceVertex || vertex == sinkVertex) return;
        Assert(getInflow(vertex) == excess[vertex], "Flow conservation not fulfilled!");
    }

    inline void checkCapacityConstraints() const noexcept {
        for (const Edge edge : graph.edges()) {
            Assert(residualCapacity[edge] >= 0, "Capacity constraint violated!");
        }
    }

private:
    const MaxFlowInstance& instance;
    const GraphType& graph;
    const int n;
    const Vertex sourceVertex;
    const Vertex sinkVertex;
    std::vector<FlowType> residualCapacity;
    std::vector<int> distance;
    std::vector<FlowType> excess;
    std::vector<Edge> currentEdge;
    VertexBuckets vertexBuckets;
    int workSinceLastUpdate;
    const int workLimit;
    Cut cut;
};
