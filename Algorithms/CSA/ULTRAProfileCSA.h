#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Journey.h"
#include "../../DataStructures/CSA/Profile/Profile.h"
#include "../../DataStructures/CSA/Profile/ContiguousProfileVector.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"
#include "../RAPTOR/InitialTransfers.h"
#include "Profiler.h"

namespace CSA {

template<bool SOURCE_PRUNING, typename PROFILER = NoProfiler, bool CONTIGUOUS_PROFILES = true>
class ULTRAProfileCSA {

public:
    using InitialTransferGraph = CHGraph;
    constexpr static bool SourcePruning = SOURCE_PRUNING;
    using Profiler = PROFILER;
    constexpr static bool ContiguousProfiles = CONTIGUOUS_PROFILES;
    using ProfileVectorType = Meta::IF<ContiguousProfiles, ContiguousProfileVector, ProfileVectorWrapper>;
    using Type = ULTRAProfileCSA<SourcePruning, Profiler, ContiguousProfiles>;

    ULTRAProfileCSA(const Data& data, const TransferGraph& reverseGraph, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        initialTransfers(chData, FORWARD, data.numberOfStops()),
        sourceVertex(noVertex),
        sourceStop(noStop),
        targetVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never),
        tripArrivalTime(data.numberOfTrips(), never),
        profiles(data.numberOfStops() + 1),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
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
        for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (stop == sourceVertex) continue;
            profiler.countMetric(METRIC_EDGES);
            const int stopDepartureTime = departureTime - reverseGraph.get(TravelTime, edge);
            if (stopDepartureTime < minDepartureTime) continue;
            if (arrivalTime - stopDepartureTime >= initialTransfers.getBackwardDistance(stop)) continue;
            if constexpr (SourcePruning) {
                if (arrivalTime - stopDepartureTime >= initialTransfers.getDistance()) continue;
                if (arrivalTime >= profiles.getEarliestArrivalTime(sourceStop, stopDepartureTime)) continue;
            }
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            profiles.addEntry(stop, stopDepartureTime, arrivalTime);
        }
        if (initialTransfers.getForwardDistance(connection.departureStopId) != INFTY) {
            profiler.countMetric(METRIC_EDGES);
            const int sourceDepartureTime = departureTime - initialTransfers.getForwardDistance(connection.departureStopId);
            if (sourceDepartureTime >= minDepartureTime && arrivalTime - sourceDepartureTime < initialTransfers.getDistance()) {
                profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
                profiles.addEntry(sourceStop, sourceDepartureTime, arrivalTime);
            }
        }
    }

    inline void scanConnections() noexcept {
        for (size_t i = data.numberOfConnections() - 1; i != size_t(-1); i--) {
            const Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;
            profiler.countMetric(METRIC_CONNECTIONS);

            const int arrivalTimeViaWalking = getArrivalTimeViaWalking(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTimeViaTrip = tripArrivalTime[connection.tripId];
            const int arrivalTimeViaTransfer = profiles.getEarliestArrivalTime(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTime = std::min(arrivalTimeViaWalking, std::min(arrivalTimeViaTrip, arrivalTimeViaTransfer));

            tripArrivalTime[connection.tripId] = arrivalTime;
            updateProfiles(connection, arrivalTime);
        }
    }

private:
    const Data& data;
    const TransferGraph& reverseGraph;
    RAPTOR::BucketCHInitialTransfers initialTransfers;

    Vertex sourceVertex;
    StopId sourceStop;
    Vertex targetVertex;
    int minDepartureTime;
    int maxDepartureTime;

    std::vector<int> tripArrivalTime;
    ProfileVectorType profiles;

    Profiler profiler;
};
}
