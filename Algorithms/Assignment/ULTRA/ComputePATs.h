#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "StopLabel.h"

#include "../../CH/Query/BucketQuery.h"

#include "../../../DataStructures/CSA/Data.h"

#include "../../../Helpers/BinarySearch.h"
#include "../../../Helpers/Vector/Vector.h"

namespace Assignment::ULTRA {

template<typename PROFILER>
class ComputePATs {

public:
    using Profiler = PROFILER;
    using Type = ComputePATs<Profiler>;

    struct ConnectionLabel {
        PerceivedTime tripPAT{Unreachable};
        PerceivedTime transferPAT{Unreachable};
        PerceivedTime skipPAT{Unreachable};
    };

public:
    ComputePATs(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, CH::BucketQuery<>& bucketQuery, const Profiler& profiler = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        bucketQuery(bucketQuery),
        connectionLabels(data.numberOfConnections()),
        tripPAT(data.numberOfTrips()),
        stopLabels(data.numberOfStops()),
        profiler(profiler) {
    }

    inline void run(const Vertex target, const int maxDelay, const int transferCost, const double walkingCosts = 0.0, const double waitingCosts = 0.0, const double finalWalkingCosts = 1.0, const int minDepartureTime = -never) noexcept {
        profiler.startInitialization();
        Vector::fill(stopLabels);
        Vector::fill(tripPAT, Unreachable);
        bucketQuery.clear();
        bucketQuery.addTarget(target);
        bucketQuery.run();
        for (const StopId stop : data.stops()) {
            if (bucketQuery.getBackwardDistance(stop) >= INFTY) continue;
            stopLabels[stop].setWalkingPTT(bucketQuery.getBackwardDistance(stop) * finalWalkingCosts);
        }
        profiler.doneInitialization();
        for (ConnectionId i = ConnectionId(data.numberOfConnections() - 1); i < data.numberOfConnections(); i--) {
            profiler.scanConnection(i);
            const CSA::Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;

            Assert(tripPAT[connection.tripId] >= connection.arrivalTime, "TripPAT is to low (Connection: " << i << ", Trip: " << connection.tripId << ", PAT: " << tripPAT[connection.tripId] << ", ArrivalTime: " << connection.arrivalTime << ", Target: " << target << ")!");
            connectionLabels[i].tripPAT = tripPAT[connection.tripId];

            StopLabel& arrivalStop = stopLabels[connection.arrivalStopId];
            connectionLabels[i].transferPAT = arrivalStop.evaluateWithDelay(connection.arrivalTime, maxDelay, waitingCosts) + transferCost;

            StopLabel& departureStop = stopLabels[connection.departureStopId];
            connectionLabels[i].skipPAT = departureStop.evaluateWithoutDelay(connection.departureTime, waitingCosts);
            profiler.evaluateProfile();

            const PerceivedTime pat = std::min(std::min(connectionLabels[i].tripPAT, targetPAT(connection)), connectionLabels[i].transferPAT);
            if (pat >= Unreachable) continue;

            tripPAT[connection.tripId] = pat;
            if (pat >= connectionLabels[i].skipPAT) continue;

            departureStop.addEntry(connection.departureTime, pat, data.minTransferTime(connection.departureStopId), waitingCosts);
            for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
                const Vertex from = reverseGraph.get(ToVertex, edge);
                if (!data.isStop(from)) continue;
                stopLabels[from].addEntry(connection.departureTime, pat, data.minTransferTime(connection.departureStopId), waitingCosts, reverseGraph.get(TravelTime, edge), walkingCosts);
                profiler.relaxEdge(edge);
            }
        }
    }

    inline const ConnectionLabel& connectionLabel(const ConnectionId i) const noexcept {
        return connectionLabels[i];
    }

    inline PerceivedTime targetPAT(const CSA::Connection& connection) const noexcept {
        return (stopLabels[connection.arrivalStopId].getWalkingPTT() < Unreachable) ? (connection.arrivalTime + stopLabels[connection.arrivalStopId].getWalkingPTT()) : (Unreachable);
    }

    inline const StopLabel& getProfile(const StopId stop) const noexcept {
        return stopLabels[stop];
    }

    inline StopLabel& getProfile(const StopId stop) noexcept {
        return stopLabels[stop];
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;

    CH::BucketQuery<>& bucketQuery;

    std::vector<ConnectionLabel> connectionLabels;
    std::vector<PerceivedTime> tripPAT;
    std::vector<StopLabel> stopLabels;

    Profiler profiler;

};

}
