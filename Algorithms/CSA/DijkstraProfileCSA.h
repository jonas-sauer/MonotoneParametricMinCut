#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "../RAPTOR/InitialTransfers.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Journey.h"
#include "../../DataStructures/CSA/Profile/ContiguousProfileVector.h"
#include "../../DataStructures/CSA/Profile/Profile.h"
#include "../../DataStructures/CSA/Profile/ProfileDijkstraLabel.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"
#include "Profiler.h"

namespace CSA {

template<typename INITIAL_TRANSFERS, bool SOURCE_PRUNING, typename PROFILER = NoProfiler, bool CONTIGUOUS_PROFILES = true>
class DijkstraProfileCSA {

public:
    using InitialTransferType = INITIAL_TRANSFERS;
    using InitialTransferGraph = typename InitialTransferType::Graph;
    constexpr static bool SourcePruning = SOURCE_PRUNING;
    using Profiler = PROFILER;
    constexpr static bool ContiguousProfiles = CONTIGUOUS_PROFILES;
    using ProfileVectorType = Meta::IF<ContiguousProfiles, ContiguousProfileVector, ProfileVectorWrapper>;
    using Type = DijkstraProfileCSA<InitialTransferType, SourcePruning, Profiler, ContiguousProfiles>;

public:
    template<typename ATTRIBUTE>
    DijkstraProfileCSA(const Data& data, const TransferGraph& reverseGraph, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        initialTransfers(forwardGraph, backwardGraph, data.numberOfStops(), weight),
        sourceVertex(noVertex),
        sourceStop(noStop),
        targetVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never),
        tripArrivalTime(data.numberOfTrips(), never),
        profiles(data.numberOfStops() + 1),
        dijkstraLabels(data.transferGraph.numVertices()),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    template<typename T = CHGraph, typename = std::enable_if_t<Meta::Equals<T, CHGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    DijkstraProfileCSA(const Data& data, const TransferGraph& reverseGraph, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        DijkstraProfileCSA(data, reverseGraph, chData.forward, chData.backward, Weight, profilerTemplate) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    DijkstraProfileCSA(const Data& data, const TransferGraph& forwardGraph, const TransferGraph& backwardGraph,  const Profiler& profilerTemplate = Profiler()) :
        DijkstraProfileCSA(data, backwardGraph, forwardGraph, backwardGraph, TravelTime, profilerTemplate) {
    }

    inline void run(const Vertex source, const Vertex target, const int minTime = 0, const int maxTime = 24 * 60 * 60) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceVertex = source;
        sourceStop = data.isStop(source) ? StopId(source) : StopId(data.numberOfStops());
        targetVertex = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        initialTransfers.run(sourceVertex, targetVertex);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections();
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.done();
    }

    inline RAPTOR::Profile getProfile(const Vertex vertex) noexcept {
        const StopId stop = (vertex == sourceVertex) ? (sourceStop) : (StopId(vertex));
        RAPTOR::Profile result;
        const int directTransferTime = getDirectWalkingTime(vertex);
        for (size_t i = 1; i < profiles.size(stop); i++) {
            const ProfileEntry& entry = profiles.entry(stop, i);
            if (entry.departureTime > maxDepartureTime) {
                if (entry.arrivalTime - maxDepartureTime >= directTransferTime) continue;
                if (i+1 < profiles.size(stop) && profiles.entry(stop, i+1).departureTime >= maxDepartureTime) continue;
            }
            result.emplace_back(entry.departureTime, entry.arrivalTime);
        }
        std::stable_sort(result.begin(), result.end());
        return result;
    }

    inline RAPTOR::ProfileHandle getProfileHandle(const Vertex vertex) noexcept {
        return RAPTOR::ProfileHandle(getProfile(vertex), minDepartureTime, maxDepartureTime, getDirectWalkingTime(vertex));
    }

    inline int getDirectWalkingTime(const Vertex vertex) noexcept {
        return (vertex == sourceVertex) ? initialTransfers.getDistance() : initialTransfers.getBackwardDistance(vertex);
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        sourceVertex = noVertex;
        sourceStop = noStop;
        targetVertex = noVertex;
        minDepartureTime = never;
        maxDepartureTime = never;
        Vector::fill(tripArrivalTime, never);
        profiles.clear();
        Vector::fill(dijkstraLabels, ProfileDijkstraLabel());
        queue.clear();
    }

    inline int getArrivalTimeViaWalking(const StopId stop, const int time) noexcept {
        const int distance = initialTransfers.getBackwardDistance(stop);
        if (distance == never) return distance;
        return time + distance;
    }

    inline void updateProfiles(const Connection& connection, const int arrivalTime) noexcept {
        const int departureTime = connection.departureTime - data.minTransferTime(connection.departureStopId);
        if (departureTime < minDepartureTime) return;
        if (arrivalTime - departureTime >= initialTransfers.getBackwardDistance(connection.departureStopId)) return;
        if constexpr (SourcePruning) {
            if (arrivalTime - departureTime >= initialTransfers.getDistance()) return;
            if (arrivalTime >= profiles.getEarliestArrivalTime(sourceStop, departureTime)) return;
        }
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        if (!profiles.addEntry(connection.departureStopId, departureTime, arrivalTime)) return;
        arrivalByEdge(connection.departureStopId, departureTime, arrivalTime);
        if (initialTransfers.getForwardDistance(connection.departureStopId) != INFTY) {
            profiler.countMetric(METRIC_EDGES);
            const int sourceDepartureTime = departureTime - initialTransfers.getForwardDistance(connection.departureStopId);
            if (sourceDepartureTime < minDepartureTime) return;
            if (arrivalTime - sourceDepartureTime >= initialTransfers.getDistance()) return;
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            profiles.addEntry(sourceStop, sourceDepartureTime, arrivalTime);
        }
    }

    inline void arrivalByEdge(const Vertex vertex, const int departureTime, const int arrivalTime) noexcept {
        if (!dijkstraLabels[vertex].insert(ProfileEntry(departureTime, arrivalTime))) return;
        queue.update(&dijkstraLabels[vertex]);
    }

    inline void scanConnections() noexcept {
        for (size_t i = data.numberOfConnections() - 1; i != size_t(-1); i--) {
            const Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;
            runDijkstra(connection.departureTime);
            profiler.countMetric(METRIC_CONNECTIONS);

            const int arrivalTimeViaWalking = getArrivalTimeViaWalking(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTimeViaTrip = tripArrivalTime[connection.tripId];
            const int arrivalTimeViaTransfer = profiles.getEarliestArrivalTime(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTime = std::min(arrivalTimeViaWalking, std::min(arrivalTimeViaTrip, arrivalTimeViaTransfer));

            tripArrivalTime[connection.tripId] = arrivalTime;
            updateProfiles(connection, arrivalTime);
        }
    }

    inline void runDijkstra(const int nextDepartureTime) noexcept {
        while ((!queue.empty()) && (queue.min().key() >= nextDepartureTime)) {
            ProfileDijkstraLabel* uLabel = queue.extractFront();
            const ProfileEntry& entry = uLabel->settle();
            if (uLabel->hasUnsettledEntries()) {
                queue.update(uLabel);
            }
            if constexpr (SourcePruning) {
                if (entry.arrivalTime - entry.departureTime >= initialTransfers.getDistance()) continue;
                if (entry.arrivalTime >= profiles.getEarliestArrivalTime(sourceStop, entry.departureTime)) continue;
            }
            const Vertex u = Vertex(uLabel - &(dijkstraLabels[0]));
            for (const Edge edge : reverseGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = reverseGraph.get(ToVertex, edge);
                if (v == sourceVertex) continue;
                const int newDepartureTime = entry.departureTime - reverseGraph.get(TravelTime, edge);
                if (newDepartureTime < minDepartureTime) continue;
                if (entry.arrivalTime - newDepartureTime >= initialTransfers.getBackwardDistance(v)) continue;
                arrivalByEdge(v, newDepartureTime, entry.arrivalTime);
            }
            if (data.isStop(u)) {
                profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
                profiles.addEntry(u, entry.departureTime, entry.arrivalTime);
            }
        }
    }

private:
    const Data& data;
    const TransferGraph& reverseGraph;
    InitialTransferType initialTransfers;

    Vertex sourceVertex;
    StopId sourceStop;
    Vertex targetVertex;
    int minDepartureTime;
    int maxDepartureTime;

    std::vector<int> tripArrivalTime;
    ProfileVectorType profiles;
    std::vector<ProfileDijkstraLabel> dijkstraLabels;
    ExternalKHeap<2, ProfileDijkstraLabel> queue;

    Profiler profiler;
};
}
