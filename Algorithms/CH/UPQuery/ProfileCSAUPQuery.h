#pragma once

#include <iostream>
#include <vector>
#include <string>

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
template<bool USE_TARGET_BUCKETS, bool STALL_ON_DEMAND = true, bool DEBUG = false>
class ProfileCSAUPQuery {

public:
    constexpr static bool UseTargetBuckets = USE_TARGET_BUCKETS;
    constexpr static bool StallOnDemand = STALL_ON_DEMAND;
    constexpr static bool Debug = DEBUG;
    using TargetGraph = Meta::IF<UseTargetBuckets, BucketGraph, SweepGraph>;
    using Type = ProfileCSAUPQuery<UseTargetBuckets, StallOnDemand, Debug>;
    using BucketBuilderType = BucketBuilder<StallOnDemand, Debug>;

private:
    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel(int* const distance) :
            ExternalKHeapElement(),
            distance(distance) {
        }

        inline int getDistance() const noexcept {
            return *distance;
        }

        inline void setDistance(const int newDistance) noexcept {
            *distance = newDistance;
        }

        inline bool hasSmallerKey(const DijkstraLabel* other) const noexcept {return getDistance() < other->getDistance();}

        int* const distance;
    };

public:
    ProfileCSAUPQuery(const CHGraph& forward, const CHGraph& backward, const Order&& order, const Vertex::ValueType numberOfStops, const IndexedSet<false, Vertex>& originalTargets) :
        graph {forward, backward},
        contractionOrder(std::move(order)),
        positionInOrder(Construct::Invert, contractionOrder),
        sweepStart(noVertex),
        stops(graph[FORWARD].numVertices(), Vector::id<Vertex>(numberOfStops)),
        targets(originalTargets),
        Q(graph[FORWARD].numVertices()),
        initialDistance(graph[FORWARD].numVertices(), INFTY),
        departureTime(graph[FORWARD].numVertices(), INFTY),
        timestamp(graph[FORWARD].numVertices(), 0),
        currentTimestamp(0),
        currentArrivalTime(never),
        bucketSources(graph[FORWARD].numVertices()),
        reachedTargets(graph[FORWARD].numVertices()) {
        reorderVertices();
        upwardSweepGraph.build(graph[FORWARD], stops, true, false);
        BucketBuilderType bucketBuilder(graph[FORWARD], graph[BACKWARD]);
        if constexpr (UseTargetBuckets) {
            targetGraph.initialize(bucketBuilder.build(targets));
        } else {
            targetGraph.build(graph[BACKWARD], targets, false, true);
        }
        for (const Vertex vertex : graph[FORWARD].vertices()) {
            dijkstraLabel.emplace_back(&initialDistance[vertex]);
        }
    }

    ProfileCSAUPQuery(const CH& ch, const Order&& order, const Vertex::ValueType numberOfStops, const IndexedSet<false, Vertex>& targets, const int direction = FORWARD) :
        ProfileCSAUPQuery(ch.getGraph(direction), ch.getGraph(!direction), std::move(order), numberOfStops, targets) {
    }

    inline void initialize() noexcept {
        clear();
    }

    inline void runInitialTransfers(const Vertex externalSource) noexcept {
        const Vertex internalSource = originalToInternal(externalSource);
        initialDistance[internalSource] = 0;
        Q.update(&dijkstraLabel[internalSource]);

        upwardSearch();
        if constexpr (UseTargetBuckets) {
            evaluateInitialTargetBuckets();
        } else {
            initialDownwardSweep();
        }
    }

    inline void setArrivalTime(const int arrivalTime) noexcept {
        currentTimestamp++;
        currentArrivalTime = arrivalTime;
        if constexpr (UseTargetBuckets) {
            bucketSources.clear();
        }
        reachedTargets.clear();
    }

    inline void addSource(const Vertex vertex, const int time) noexcept {
        addSourceInternal(originalToInternal(vertex), time);
    }

    inline void runFinalTransfers() noexcept {
        upwardSweep();
        if constexpr (UseTargetBuckets) {
            evaluateFinalTargetBuckets();
        } else {
            finalDownwardSweep();
        }
    }

    inline int getInitialDistance(const Vertex vertex) const noexcept {
        const Vertex internalVertex = originalToInternal(vertex);
        return initialDistance[internalVertex];
    }

    inline int getDepartureTime(const Vertex vertex) const noexcept {
        const Vertex internalVertex = originalToInternal(vertex);
        return departureTime[internalVertex];
    }

    inline const IndexedSet<false, Vertex> getReachedTargets() const noexcept {
        return reachedTargets;
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

    inline void clear() noexcept {
        Q.clear();
        Vector::fill(initialDistance, INFTY);
        Vector::fill(departureTime, INFTY);
        Vector::fill(timestamp, 0);
        bucketSources.clear();
        reachedTargets.clear();
        sweepStart = noVertex;
        currentArrivalTime = never;
    }

    inline void addSourceInternal(const Vertex vertex, const int time) noexcept {
        if (time >= departureTime[vertex]) return;
        departureTime[vertex] = time;
        AssertMsg(upwardSweepGraph.externalToInternal(vertex) < upwardSweepGraph.graph.numVertices(), "Vertex is not in sweep graph! (original: "<< internalToOriginal(vertex) << ", CH: " << vertex << ", sweep: " << upwardSweepGraph.externalToInternal(vertex) << ")");
        sweepStart = std::min(sweepStart, upwardSweepGraph.externalToInternal(vertex));
        if constexpr (UseTargetBuckets) {
            bucketSources.insert(vertex);
        }
        timestamp[vertex] = currentTimestamp;
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
        const Vertex u = Vertex(Q.extractFront() - &(dijkstraLabel[0]));
        AssertMsg(u < graph[FORWARD].numVertices(), u << " is not a valid vertex!");
        if constexpr (StallOnDemand) {
            for (Edge edge : graph[BACKWARD].edgesFrom(u)) {
                const Vertex v = graph[BACKWARD].get(ToVertex, edge);
                if (initialDistance[v] < initialDistance[u] - graph[BACKWARD].get(Weight, edge)) return;
            }
        }
        for (Edge edge : graph[FORWARD].edgesFrom(u)) {
            const Vertex v = graph[FORWARD].get(ToVertex, edge);
            const int newDistance = initialDistance[u] + graph[FORWARD].get(Weight, edge);
            if (initialDistance[v] > newDistance) {
                Q.update(&dijkstraLabel[v]);
                initialDistance[v] = newDistance;
            }
        }
        bucketSources.insert(u);
    }

    inline void upwardSweep() noexcept {
        if constexpr (Debug) {
            std::cout << "Running upward sweep from " << sweepStart << "/" << Vertex(upwardSweepGraph.graph.numVertices()) << std::endl;
            timer.restart();
        }

        for (Vertex sweepV = sweepStart; sweepV < upwardSweepGraph.graph.numVertices(); sweepV++) {
            const Vertex v = upwardSweepGraph.internalToExternal(sweepV);
            if (timestamp[v] != currentTimestamp) continue;
            if constexpr (UseTargetBuckets) {
                bucketSources.insert(v);
            }
            if (targets.contains(v)) {
                reachedTargets.insert(internalToOriginal(v));
            }
            for (const Edge edge : upwardSweepGraph.graph.edgesFrom(sweepV)) {
                const Vertex u = upwardSweepGraph.toVertex[edge];
                const int weight = upwardSweepGraph.graph.get(Weight, edge);
                const int newDistance = departureTime[v] + weight;
                const bool update = newDistance < departureTime[u];
                departureTime[u] = branchlessConditional(update, newDistance, departureTime[u]);
                timestamp[u] = branchlessConditional(update, currentTimestamp, timestamp[u]);
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    template<bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && !T>>
    inline void initialDownwardSweep() noexcept {
        downwardSweep<false>(initialDistance);

    }

    template<bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && !T>>
    inline void finalDownwardSweep() noexcept {
        downwardSweep<true>(departureTime);
    }

    template<bool FINAL, bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && !T>>
    inline void downwardSweep(std::vector<int>& distance) noexcept {
        if constexpr (Debug) {
            std::cout << "Running " << (FINAL ? "initial" : "final") << " downward sweep" << std::endl;
            timer.restart();
        }

        for (const Vertex sweepV : targetGraph.graph.vertices()) {
            const Vertex v = targetGraph.internalToExternal(sweepV);
            if constexpr (FINAL) {
                if (timestamp[v] != currentTimestamp) continue;
            }
            if constexpr (FINAL) {
                if (targets.contains(v)) {
                    reachedTargets.insert(internalToOriginal(v));
                }
            }
            for (const Edge edge : targetGraph.graph.edgesFrom(sweepV)) {
                const Vertex u = targetGraph.toVertex[edge];
                const int weight = targetGraph.graph.get(Weight, edge);
                const int newDistance = distance[v] + weight;
                const bool update = newDistance < distance[u];
                distance[u] = branchlessConditional(update, newDistance, distance[u]);
                if constexpr (FINAL) {
                    timestamp[u] = branchlessConditional(update, currentTimestamp, timestamp[u]);
                }
            }
        }

        if constexpr (Debug) {
            std::cout << "Time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        }
    }

    template<bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && T>>
    inline void evaluateInitialTargetBuckets() noexcept {
        evaluateTargetBuckets<false>(initialDistance);
    }

    template<bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && T>>
    inline void evaluateFinalTargetBuckets() noexcept {
        evaluateTargetBuckets<true>(departureTime);
    }

    template<bool FINAL, bool T = UseTargetBuckets, typename = std::enable_if_t<T == UseTargetBuckets && T>>
    inline void evaluateTargetBuckets(std::vector<int>& distance) noexcept {
        if constexpr (Debug) {
            std::cout << "Evaluating " << (FINAL ? "initial" : "final") << " target buckets" << std::endl;
            timer.restart();
        }

        for (const Vertex vertex : bucketSources) {
            for (const Edge edge : targetGraph.graph.edgesFrom(vertex)) {
                const int newDistance = distance[vertex] + targetGraph.graph.get(Weight, edge);
                const Vertex poi = targetGraph.graph.get(ToVertex, edge);
                if (newDistance < distance[poi]) {
                    distance[poi] = newDistance;
                    if constexpr (FINAL) {
                        if (targets.contains(poi)) {
                            reachedTargets.insert(internalToOriginal(poi));
                        }
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

private:
    CHGraph graph[2];
    SweepGraph upwardSweepGraph;
    TargetGraph targetGraph;

    const Order contractionOrder;
    const Permutation positionInOrder;

    Vertex sweepStart;
    IndexedSet<false, Vertex> stops;
    IndexedSet<false, Vertex> targets;

    ExternalKHeap<2, DijkstraLabel> Q;
    std::vector<DijkstraLabel> dijkstraLabel;
    std::vector<int> initialDistance;
    std::vector<int> departureTime;
    std::vector<int> timestamp;
    int currentTimestamp;
    int currentArrivalTime;

    IndexedSet<false, Vertex> bucketSources;
    IndexedSet<false, Vertex> reachedTargets;

    Timer timer;
};

}
