#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../DataStructures/Assignment/Profile.h"
#include "../../DataStructures/Assignment/StopLabel.h"
#include "../../DataStructures/CSA/Data.h"

#include "../../Helpers/Vector/Vector.h"

#include "Profiler.h"

namespace Assignment {

template<typename PROFILER = NoPATProfiler, bool USE_TRANSFER_BUFFER_TIMES = false>
class ComputePATs {

public:
    using Profiler = PROFILER;
    static constexpr bool UseTransferBufferTimes = USE_TRANSFER_BUFFER_TIMES;
    using Type = ComputePATs<Profiler, UseTransferBufferTimes>;

    struct ConnectionLabel {
        PerceivedTime tripPAT{Unreachable};
        PerceivedTime transferPAT{Unreachable};
        PerceivedTime skipPAT{Unreachable};
    };

public:
    ComputePATs(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Profiler& profiler = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        connectionLabels(data.numberOfConnections()),
        tripPAT(data.numberOfTrips(), Unreachable),
        stopLabels(data.numberOfStops()),
        transferDistanceToTarget(data.numberOfStops(), INFTY),
        targetVertex(noVertex),
        profiler(profiler) {
    }

    inline void run(const Vertex target, const int maxDelay, const int transferCost, const double walkingCosts = 0.0, const double waitingCosts = 0.0) noexcept {
        profiler.startInitialization();
        clear();
        initialize(target, walkingCosts);
        profiler.doneInitialization();
        for (ConnectionId i = ConnectionId(data.numberOfConnections() - 1); i < data.numberOfConnections(); i--) {
            profiler.scanConnection(i);
            const CSA::Connection& connection = data.connections[i];
            const ProfileEntry& skipEntry = stopLabels[connection.departureStopId].getSkipEntry();

            Assert(skipEntry.departureTime >= connection.departureTime, "Connections are scanned out of order (" << skipEntry.departureTime << " before " << connection.departureTime << ", index: " << i << ")!");
            connectionLabels[i].tripPAT = tripPAT[connection.tripId];
            connectionLabels[i].transferPAT = stopLabels[connection.arrivalStopId].evaluateWithDelay(connection.arrivalTime, maxDelay, waitingCosts) + transferCost;
            connectionLabels[i].skipPAT = skipEntry.evaluate(connection.departureTime, waitingCosts);
            profiler.evaluateProfile();

            const PerceivedTime pat = std::min(std::min(connectionLabels[i].tripPAT, targetPAT(connection)), connectionLabels[i].transferPAT);
            tripPAT[connection.tripId] = pat;
            if (pat >= connectionLabels[i].skipPAT) continue;

            Assert(pat < Unreachable, "Adding infinity PAT = " << pat << "!");
            stopLabels[connection.departureStopId].addWaitingEntry(ProfileEntry(connection.departureTime, i, pat, waitingCosts));
            profiler.addToProfile();
            const int bufferTime = data.minTransferTime(connection.departureStopId);
            profiler.insertToProfile();
            stopLabels[connection.departureStopId].addTransferEntry(ProfileEntry(connection.departureTime, i, pat, 0, bufferTime, walkingCosts, waitingCosts), profiler);
            for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
                const Vertex from = reverseGraph.get(ToVertex, edge);
                if (!data.isStop(from)) continue;
                profiler.insertToProfile();
                stopLabels[from].addTransferEntry(ProfileEntry(connection.departureTime, i, pat, reverseGraph.get(TravelTime, edge), UseTransferBufferTimes ? bufferTime : 0, walkingCosts, waitingCosts), profiler);
                profiler.relaxEdge(edge);
            }
        }
    }

    inline const ConnectionLabel& connectionLabel(const ConnectionId i) const noexcept {
        return connectionLabels[i];
    }

    inline PerceivedTime targetPAT(const CSA::Connection& connection) const noexcept {
        const int distance = transferDistanceToTarget[connection.arrivalStopId];
        return (distance < INFTY) ? (connection.arrivalTime + distance) : Unreachable;
    }

    inline const Profile& getProfile(const StopId source) const noexcept {
        return stopLabels[source].getWaitingProfile();
    }

    inline Profiler& getProfiler() noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        Vector::fill(tripPAT, Unreachable);
        Vector::fill(stopLabels);
        if (reverseGraph.isVertex(targetVertex)) cleanUp();
    }

    inline void initialize(const Vertex target, const double walkingCosts = 0.0) noexcept {
        targetVertex = target;
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            transferDistanceToTarget[stop] = (walkingCosts + 1) * reverseGraph.get(TravelTime, edge);
        }
        if (data.isStop(targetVertex)) transferDistanceToTarget[targetVertex] = 0;
    }

    inline void cleanUp() noexcept {
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            transferDistanceToTarget[stop] = INFTY;
        }
        if (data.isStop(targetVertex)) transferDistanceToTarget[targetVertex] = INFTY;
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;

    std::vector<ConnectionLabel> connectionLabels;
    std::vector<PerceivedTime> tripPAT;
    std::vector<StopLabel> stopLabels;
    std::vector<int> transferDistanceToTarget;
    Vertex targetVertex;

    Profiler profiler;

};

}
