#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>
#include <concepts>

#include "../../../DataStructures/CH/UPGraphs.h"
#include "../CH.h"
#include "../CHUtils.h"

#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/Timer.h"
#include "../../../Helpers/Console/Progress.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/Vector/Vector.h"

#include "../../../DataStructures/Container/ExternalKHeap.h"
#include "../../../DataStructures/Container/Set.h"
#include "BucketBuilder.h"

namespace CH {

//Assumption: Targets include all stops
template<bool USE_TARGET_BUCKETS, bool STALL_ON_DEMAND = true, bool DEBUG = false, size_t MAX_SOURCES = 16>
class ProfileUPQuery {

public:
    constexpr static bool UseTargetBuckets = USE_TARGET_BUCKETS;
    constexpr static bool StallOnDemand = STALL_ON_DEMAND;
    constexpr static bool Debug = DEBUG;
    using TargetGraph = std::conditional_t<UseTargetBuckets, BucketGraph, SweepGraph>;
    constexpr static size_t MaxSources = MAX_SOURCES;
    using Type = ProfileUPQuery<UseTargetBuckets, StallOnDemand, Debug, MaxSources>;
    using BucketBuilderType = BucketBuilder<StallOnDemand, Debug>;

private:
    struct Label {
        Label() {
            reset();
        }

        inline void reset() noexcept {
            std::fill(arrivalTime, arrivalTime + MaxSources, INFTY);
            std::fill(parent, parent + MaxSources, noVertex);
        }

        inline void clear() noexcept {
            std::fill(arrivalTime, arrivalTime + MaxSources, arrivalTime[MaxSources - 1]);
            std::fill(parent, parent + MaxSources, parent[MaxSources - 1]);
        }

        int arrivalTime[MaxSources];
        int parent[MaxSources];
        int timestamp;
    };

    struct Distance : public ExternalKHeapElement {
        Distance() : ExternalKHeapElement(), distance(INFTY) {}
        inline bool hasSmallerKey(const Distance* other) const noexcept {return distance < other->distance;}
        int distance;
    };

public:
    ProfileUPQuery(const CHGraph& forward, const CHGraph& backward, const Order&& order, const Vertex::ValueType numberOfStops, const IndexedSet<false, Vertex>& originalTargets) :
        graph {forward, backward},
        contractionOrder(std::move(order)),
        positionInOrder(Construct::Invert, contractionOrder),
        sweepStart(noVertex),
        stops(graph[FORWARD].numVertices(), Vector::id<Vertex>(numberOfStops)),
        targets(originalTargets),
        Q(graph[FORWARD].numVertices()),
        distance(graph[FORWARD].numVertices()),
        label(graph[FORWARD].numVertices()),
        numDepartureTimes(0),
        timestamp(0),
        bucketSources(graph[FORWARD].numVertices()),
        reachedStops(numberOfStops) {
        std::fill(departureTimes, departureTimes + MaxSources, never);
        reorderVertices();
        buildUpwardSweepGraph();
        if constexpr (UseTargetBuckets) {
            BucketBuilderType bucketBuilder(graph[FORWARD], graph[BACKWARD]);
            targetGraph.initialize(bucketBuilder.build(targets));
        } else {
            targetGraph.build(graph[BACKWARD], targets, false, false);
        }
    }

    ProfileUPQuery(const CH& ch, const Order&& order, const Vertex::ValueType numberOfStops, const IndexedSet<false, Vertex>& targets, const int direction = FORWARD) :
        ProfileUPQuery(ch.getGraph(direction), ch.getGraph(!direction), std::move(order), numberOfStops, targets) {
    }

    inline void initialize() noexcept {
        clear();
    }

    inline void runInitialTransfers(const Vertex externalSource) noexcept {
        const Vertex internalSource = originalToInternal(externalSource);
        distance[internalSource].distance = 0;
        Q.update(&distance[internalSource]);

        upwardSearch();
        if constexpr (UseTargetBuckets) {
            evaluateInitialTargetBuckets();
        } else {
            initialDownwardSweep();
        }
    }

    inline void clearDepartureTimes() noexcept {
        timestamp++;
        if constexpr (UseTargetBuckets) {
            bucketSources.clear();
        }
        std::fill(departureTimes, departureTimes + MaxSources, never);
        numDepartureTimes = 0;
    }

    inline void addDepartureTime(const int departureTime) noexcept {
        numDepartureTimes++;
        departureTimes[numDepartureTimes - 1] = departureTime;
        Assert(numDepartureTimes <= MaxSources, "Exceeded maximum number of sources!");
    }

    inline void addSource(const Vertex vertex, const int initialArrivalTime, const Vertex parentVertex) noexcept {
        addSourceInternal(originalToInternal(vertex), initialArrivalTime, parentVertex);
    }

    inline void runFinalTransfers() noexcept {
        upwardSweep();
        if constexpr (UseTargetBuckets) {
            evaluateFinalTargetBuckets();
        } else {
            finalDownwardSweep();
        }
    }

    inline int getInitialDistance(const Vertex vertex) noexcept {
        return distance[originalToInternal(vertex)].distance;
    }

    inline int getArrivalTime(const size_t departure, const Vertex vertex) noexcept {
        return getLabel(originalToInternal(vertex)).arrivalTime[departure];
    }

    inline Vertex getParent(const size_t departure, const Vertex vertex) noexcept {
        return Vertex(getLabel(originalToInternal(vertex)).parent[departure]);
    }

    inline size_t numDepartures() const noexcept {
        return numDepartureTimes;
    }

    inline const IndexedSet<false, Vertex>& getReachedStops() const noexcept {
        return reachedStops;
    }

private:
    inline void reorderVertices() noexcept {
        reorder(graph[FORWARD]);
        reorder(graph[BACKWARD]);
        stops.applyPermutation(positionInOrder);
        targets.applyPermutation(positionInOrder);
    }

    inline void reorder(CHGraph& graph) noexcept {
        graph.applyVertexPermutation(positionInOrder);
        graph.sortEdges(ToVertex);
    }

    inline void buildUpwardSweepGraph() noexcept {
        upwardSweepGraph.build(graph[FORWARD], stops, true, true);
        sweepStartOf.resize(upwardSweepGraph.graph.numVertices(), noVertex);
        for (const Vertex to : upwardSweepGraph.graph.vertices()) {
            for (const Edge edge : upwardSweepGraph.graph.edgesFrom(to)) {
                const Vertex from = upwardSweepGraph.graph.get(ToVertex, edge);
                if (sweepStartOf[from] != noVertex) continue;
                sweepStartOf[from] = to;
            }
        }
    }

    inline void clear() noexcept {
        Q.clear();
        std::vector<Distance>(graph[FORWARD].numVertices()).swap(distance);
        std::vector<Label>(graph[FORWARD].numVertices()).swap(label);
        timestamp = 0;
        bucketSources.clear();
        reachedStops.clear();
        sweepStart = noVertex;
        std::fill(departureTimes, departureTimes + MaxSources, never);
        numDepartureTimes = 0;
    }

    inline void addSourceInternal(const Vertex vertex, const int initialArrivalTime, const Vertex parentVertex) noexcept {
        addSourceInternal(vertex, initialArrivalTime, parentVertex, numDepartureTimes - 1);
    }

    inline void addSourceInternal(const Vertex vertex, const int initialArrivalTime, const Vertex parentVertex, const size_t departure) noexcept {
        Label& vertexLabel = getLabel(vertex);
        if (initialArrivalTime >= vertexLabel.arrivalTime[departure]) return;
        vertexLabel.arrivalTime[departure] = initialArrivalTime;
        vertexLabel.parent[departure] = parentVertex;
        Assert(upwardSweepGraph.externalToInternal(vertex) < upwardSweepGraph.graph.numVertices(), "Vertex is not in sweep graph! (original: "<< internalToOriginal(vertex) << ", CH: " << vertex << ", sweep: " << upwardSweepGraph.externalToInternal(vertex) << ")");
        sweepStart = std::min(sweepStart, sweepStartOf[upwardSweepGraph.externalToInternal(vertex)]);
        if constexpr (UseTargetBuckets) {
            bucketSources.insert(vertex);
        }
    }

    inline void upwardSearch() noexcept {
        if constexpr (Debug) {
            std::cout << "Running upward search with " << Q.size() << " queue elements" << std::endl;
            timer.restart();
        }

        while (!Q.empty()) {
            settle();
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline void settle() noexcept {
        Distance* distanceU = Q.extractFront();
        const Vertex u = Vertex(distanceU - &(distance[0]));
        Assert(u < graph[FORWARD].numVertices(), u << " is not a valid vertex!");
        if constexpr (StallOnDemand) {
            for (Edge edge : graph[BACKWARD].edgesFrom(u)) {
                const Vertex v = graph[BACKWARD].get(ToVertex, edge);
                if (distance[v].distance < distance[u].distance - graph[BACKWARD].get(Weight, edge)) return;
            }
        }
        for (Edge edge : graph[FORWARD].edgesFrom(u)) {
            const Vertex v = graph[FORWARD].get(ToVertex, edge);
            const int newDistance = distanceU->distance + graph[FORWARD].get(Weight, edge);
            if (distance[v].distance > newDistance) {
                distance[v].distance = newDistance;
                Q.update(&distance[v]);
            }
        }
        if constexpr (UseTargetBuckets) {
            bucketSources.insert(u);
        }
        if (stops.contains(u)) {
            reachedStops.insert(internalToOriginal(u));
        }
    }

    inline void upwardSweep() noexcept {
        if constexpr (Debug) {
            std::cout << "Running upward sweep from " << sweepStart << "/" << Vertex(upwardSweepGraph.graph.numVertices()) << std::endl;
            timer.restart();
        }

        for (Vertex sweepV = sweepStart; sweepV < upwardSweepGraph.graph.numVertices(); sweepV++) {
            const Vertex v = upwardSweepGraph.internalToExternal(sweepV);
            Label& vLabel = getLabel(v);
            bool reached = false;
            for (const Edge edge : upwardSweepGraph.graph.edgesFrom(sweepV)) {
                const Vertex u = upwardSweepGraph.toVertex[edge];
                const Label& uLabel = getLabel(u);
                const int weight = upwardSweepGraph.graph.get(Weight, edge);
                for (size_t i = 0; i < numDepartureTimes; i++) {
                    const int newArrivalTime = uLabel.arrivalTime[i] + weight;
                    const bool update = newArrivalTime < vLabel.arrivalTime[i];
                    vLabel.arrivalTime[i] = branchlessConditional(update, newArrivalTime, vLabel.arrivalTime[i]);
                    vLabel.parent[i] = branchlessConditional(update, uLabel.parent[i], vLabel.parent[i]);
                    if constexpr (UseTargetBuckets) reached |= update;
                }
            }
            if constexpr (UseTargetBuckets) {
                if (reached) bucketSources.insert(v);
            } else {
                suppressUnusedParameterWarning(reached);
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline void initialDownwardSweep() noexcept requires (!UseTargetBuckets) {
        if constexpr (Debug) {
            std::cout << "Running downward sweep for initial transfers" << std::endl;
            timer.restart();
        }

        for (const Vertex sweepV : targetGraph.graph.vertices()) {
            const Vertex v = targetGraph.internalToExternal(sweepV);
            bool reached = false;
            for (const Edge edge : targetGraph.graph.edgesFrom(sweepV)) {
                const Vertex u = targetGraph.toVertex[edge];
                const int newDistance = distance[u].distance + targetGraph.graph.get(Weight, edge);
                const bool update = newDistance < distance[v].distance;
                distance[v].distance = branchlessConditional(update, newDistance, distance[v].distance);
                reached |= update;
            }
            if (stops.contains(v)) {
                reachedStops.insert(internalToOriginal(v));
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline void finalDownwardSweep() noexcept requires (!UseTargetBuckets) {
        if constexpr (Debug) {
            std::cout << "Running downward sweep for final transfers" << std::endl;
            timer.restart();
        }

        for (const Vertex sweepV : targetGraph.graph.vertices()) {
            const Vertex v = targetGraph.internalToExternal(sweepV);
            Label& vLabel = getLabel(v);
            for (const Edge edge : targetGraph.graph.edgesFrom(sweepV)) {
                const Vertex u = targetGraph.toVertex[edge];
                const Label& uLabel = label[u]; //Already known to be up to date
                const int weight = targetGraph.graph.get(Weight, edge);
                for (size_t i = 0; i < numDepartureTimes; i++) {
                    const int newArrivalTime = uLabel.arrivalTime[i] + weight;
                    const bool update = newArrivalTime < vLabel.arrivalTime[i];
                    vLabel.arrivalTime[i] = branchlessConditional(update, newArrivalTime, vLabel.arrivalTime[i]);
                    vLabel.parent[i] = branchlessConditional(update, uLabel.parent[i], vLabel.parent[i]);
                }
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline void evaluateInitialTargetBuckets() noexcept requires UseTargetBuckets {
        if constexpr (Debug) {
            std::cout << "Evaluating target buckets for initial transfers" << std::endl;
            timer.restart();
        }

        for (const Vertex vertex : bucketSources) {
            for (const Edge edge : targetGraph.graph.edgesFrom(vertex)) {
                const int newDistance = distance[vertex].distance + targetGraph.graph.get(Weight, edge);
                const Vertex poi = targetGraph.graph.get(ToVertex, edge);
                distance[poi].distance = std::min(distance[poi].distance, newDistance);
                if (stops.contains(poi)) {
                    reachedStops.insert(internalToOriginal(poi));
                }
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline void evaluateFinalTargetBuckets() noexcept requires UseTargetBuckets {
        if constexpr (Debug) {
            std::cout << "Evaluating target buckets for final transfers" << std::endl;
            timer.restart();
        }

        for (const Vertex vertex : bucketSources) {
            const Label& vertexLabel = label[vertex];
            for (const Edge edge : targetGraph.graph.edgesFrom(vertex)) {
                for (size_t i = 0; i < numDepartureTimes; i++) {
                    if (vertexLabel.arrivalTime[i] == INFTY) continue;
                    Assert(label[vertex].parent[i] != int(noVertex), "Invalid parent!");
                    const int newArrivalTime = vertexLabel.arrivalTime[i] + targetGraph.graph.get(Weight, edge);
                    const Vertex poi = targetGraph.graph.get(ToVertex, edge);
                    Label& poiLabel = getLabel(poi);
                    if (newArrivalTime < poiLabel.arrivalTime[i]) {
                        poiLabel.arrivalTime[i] = newArrivalTime;
                        poiLabel.parent[i] = vertexLabel.parent[i];
                    }
                }
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    inline Vertex originalToInternal(const Vertex vertex) const noexcept {
        return positionInOrder.permutate(vertex);
    }

    inline Vertex internalToOriginal(const Vertex vertex) const noexcept {
        return Vertex(contractionOrder[vertex]);
    }

    inline Label& getLabel(const Vertex vertex) noexcept {
        Label& vertexLabel = label[vertex];
        if (vertexLabel.timestamp != timestamp) {
            vertexLabel.clear();
            vertexLabel.timestamp = timestamp;
        }
        return vertexLabel;
    }

private:
    CHGraph graph[2];
    SweepGraph upwardSweepGraph;
    TargetGraph targetGraph;

    const Order contractionOrder;
    const Permutation positionInOrder;

    Vertex sweepStart;
    std::vector<Vertex> sweepStartOf;
    IndexedSet<false, Vertex> stops;
    IndexedSet<false, Vertex> targets;

    ExternalKHeap<2, Distance> Q;
    std::vector<Distance> distance;
    std::vector<Label> label;
    int departureTimes[MaxSources];
    size_t numDepartureTimes;
    int timestamp;

    IndexedSet<false, Vertex> bucketSources;
    IndexedSet<false, Vertex> reachedStops;

    Timer timer;
};

}
