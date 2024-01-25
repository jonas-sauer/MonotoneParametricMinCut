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
#include "Profiler.h"

namespace CSA {

template<bool SOURCE_PRUNING, typename PROFILER = NoProfiler, bool CONTIGUOUS_PROFILES = true>
class ProfileCSA {

public:
    constexpr static bool SourcePruning = SOURCE_PRUNING;
    using Profiler = PROFILER;
    constexpr static bool ContiguousProfiles = CONTIGUOUS_PROFILES;
    using ProfileVectorType = Meta::IF<ContiguousProfiles, ContiguousProfileVector, ProfileVectorWrapper>;
    using Type = ProfileCSA<SourcePruning, Profiler, ContiguousProfiles>;

    ProfileCSA(const Data& data, const TransferGraph& reverseGraph, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        sourceStop(noStop),
        targetStop(noStop),
        minDepartureTime(never),
        maxDepartureTime(never),
        tripArrivalTime(data.numberOfTrips(), never),
        profiles(data.numberOfStops()),
        transferDistanceToTarget(data.numberOfStops(), never),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline void run(const StopId source, const StopId target, const int minTime = 0, const int maxTime = 24 * 60 * 60) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceStop = source;
        targetStop = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        computeFinalTransfers();
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections();
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.done();
    }

    inline RAPTOR::Profile getProfile(const StopId stop) const noexcept {
        RAPTOR::Profile result;
        const int directTransferTime = transferDistanceToTarget[stop];
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

    inline RAPTOR::ProfileHandle getProfileHandle(const StopId stop) const noexcept {
        return RAPTOR::ProfileHandle(getProfile(stop), minDepartureTime, maxDepartureTime, transferDistanceToTarget[stop]);
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        sourceStop = noStop;
        targetStop = noStop;
        minDepartureTime = never;
        maxDepartureTime = never;
        Vector::fill(tripArrivalTime, never);
        profiles.clear();
        Vector::fill(transferDistanceToTarget, never);
    }

    inline void computeFinalTransfers() noexcept {
        transferDistanceToTarget[targetStop] = 0;
        for (const Edge edge : reverseGraph.edgesFrom(targetStop)) {
            profiler.countMetric(METRIC_EDGES);
            const Vertex vertex = reverseGraph.get(ToVertex, edge);
            transferDistanceToTarget[vertex] = reverseGraph.get(TravelTime, edge);
        }
    }

    inline int getArrivalTimeViaWalking(const StopId stop, const int time) const noexcept {
        const int distance = transferDistanceToTarget[stop];
        if (distance == never) return distance;
        return time + distance;
    }

    inline void updateProfiles(const Connection& connection, const int arrivalTime) noexcept {
        const int departureTime = connection.departureTime - data.minTransferTime(connection.departureStopId);
        if (departureTime < minDepartureTime) return;
        if (arrivalTime - departureTime >= transferDistanceToTarget[connection.departureStopId]) return;
        if constexpr (SourcePruning) {
            if (arrivalTime - departureTime >= transferDistanceToTarget[sourceStop]) return;
            if (arrivalTime > profiles.getEarliestArrivalTime(sourceStop, departureTime)) return;
        }
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        if (!profiles.addEntry(connection.departureStopId, departureTime, arrivalTime)) return;
        for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
            profiler.countMetric(METRIC_EDGES);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            const int stopDepartureTime = departureTime - reverseGraph.get(TravelTime, edge);
            if (stopDepartureTime < minDepartureTime) continue;
            if (arrivalTime - stopDepartureTime >= transferDistanceToTarget[stop]) continue;
            if constexpr (SourcePruning) {
                if (arrivalTime - stopDepartureTime >= transferDistanceToTarget[sourceStop]) continue;
                if (arrivalTime > profiles.getEarliestArrivalTime(sourceStop, stopDepartureTime)) continue;
            }
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            profiles.addEntry(stop, stopDepartureTime, arrivalTime);
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

    StopId sourceStop;
    StopId targetStop;
    int minDepartureTime;
    int maxDepartureTime;

    std::vector<int> tripArrivalTime;
    ProfileVectorType profiles;
    std::vector<int> transferDistanceToTarget;

    Profiler profiler;
};
}
