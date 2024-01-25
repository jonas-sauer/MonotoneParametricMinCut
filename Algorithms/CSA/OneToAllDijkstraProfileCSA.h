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

template<typename PROFILER = NoProfiler, bool CONTIGUOUS_PROFILES = true>
class OneToAllDijkstraProfileCSA {

public:
    using Profiler = PROFILER;
    constexpr static bool ContiguousProfiles = CONTIGUOUS_PROFILES;
    using ProfileVectorType = Meta::IF<ContiguousProfiles, ContiguousProfileVector, ProfileVectorWrapper>;
    using Type = OneToAllDijkstraProfileCSA<Profiler, ContiguousProfiles>;

private:
    struct InitialTransferLabel : public ExternalKHeapElement {
        InitialTransferLabel() : ExternalKHeapElement(), distance(INFTY) {}
        inline bool hasSmallerKey(const InitialTransferLabel* other) const noexcept {return distance < other->distance;}
        int distance;
    };

    public:
    OneToAllDijkstraProfileCSA(const Data& data, const TransferGraph& reverseGraph, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        targetVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never),
        tripArrivalTime(data.numberOfTrips(), never),
        profiles(data.transferGraph.numVertices()),
        initialTransferDistance(data.transferGraph.numVertices()),
        dijkstraLabels(data.transferGraph.numVertices()),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN, PHASE_FINAL_TRANSFERS});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline void run(const Vertex target, const int minTime = 0, const int maxTime = 24 * 60 * 60) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        targetVertex = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        runInitialTransfers();
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections();
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.startPhase();
        runDijkstra(-never);
        profiler.donePhase(PHASE_FINAL_TRANSFERS);

        profiler.done();
    }

    inline RAPTOR::Profile getProfile(const Vertex vertex) noexcept {
        RAPTOR::Profile result;
        const int directTransferTime = getDirectWalkingTime(vertex);
        for (size_t i = 1; i < profiles.size(vertex); i++) {
            const ProfileEntry& entry = profiles.entry(vertex, i);
            if (entry.departureTime > maxDepartureTime) {
                if (entry.arrivalTime - maxDepartureTime >= directTransferTime) continue;
                if (i+1 < profiles.size(vertex) && profiles.entry(vertex, i+1).departureTime >= maxDepartureTime) continue;
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
        return initialTransferDistance[vertex].distance;
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        targetVertex = noVertex;
        minDepartureTime = never;
        maxDepartureTime = never;
        Vector::fill(tripArrivalTime, never);
        profiles.clear();
        Vector::fill(initialTransferDistance, InitialTransferLabel());
        initialQueue.clear();
        Vector::fill(dijkstraLabels, ProfileDijkstraLabel());
        queue.clear();
    }

    inline int getArrivalTimeViaWalking(const StopId stop, const int time) noexcept {
        const int distance = initialTransferDistance[stop].distance;
        if (distance == never) return distance;
        return time + distance;
    }

    inline void updateProfiles(const Connection& connection, const int arrivalTime) noexcept {
        const int departureTime = connection.departureTime - data.minTransferTime(connection.departureStopId);
        if (departureTime < minDepartureTime) return;
        if (arrivalTime - departureTime >= initialTransferDistance[connection.departureStopId].distance) return;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        if (!profiles.addEntry(connection.departureStopId, departureTime, arrivalTime)) return;
        arrivalByEdge(connection.departureStopId, departureTime, arrivalTime);
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

    inline void runInitialTransfers() noexcept {
        initialTransferDistance[targetVertex].distance = 0;
        initialQueue.update(&initialTransferDistance[targetVertex]);
        while (!initialQueue.empty()) {
            InitialTransferLabel* uLabel = initialQueue.extractFront();
            const Vertex u = Vertex(uLabel - &(initialTransferDistance[0]));
            for (const Edge edge : reverseGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = reverseGraph.get(ToVertex, edge);
                const int newDistance = uLabel->distance + reverseGraph.get(TravelTime, edge);
                if (newDistance < initialTransferDistance[v].distance) {
                    initialTransferDistance[v].distance = newDistance;
                    initialQueue.update(&initialTransferDistance[v]);
                }
            }
        }
    }

    inline void runDijkstra(const int nextDepartureTime) noexcept {
        while ((!queue.empty()) && (queue.min().key() >= nextDepartureTime)) {
            ProfileDijkstraLabel* uLabel = queue.extractFront();
            const ProfileEntry& entry = uLabel->settle();
            if (uLabel->hasUnsettledEntries()) {
                queue.update(uLabel);
            }
            const Vertex u = Vertex(uLabel - &(dijkstraLabels[0]));
            for (const Edge edge : reverseGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = reverseGraph.get(ToVertex, edge);
                const int newDepartureTime = entry.departureTime - reverseGraph.get(TravelTime, edge);
                if (newDepartureTime < minDepartureTime) continue;
                if (entry.arrivalTime - newDepartureTime >= initialTransferDistance[v].distance) continue;
                arrivalByEdge(v, newDepartureTime, entry.arrivalTime);
            }
            if (data.isStop(u)) {
                profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            }
            profiles.addEntry(u, entry.departureTime, entry.arrivalTime);
        }
    }

private:
    const Data& data;
    const TransferGraph& reverseGraph;

    Vertex targetVertex;
    int minDepartureTime;
    int maxDepartureTime;

    std::vector<int> tripArrivalTime;
    ProfileVectorType profiles;
    std::vector<InitialTransferLabel> initialTransferDistance;
    ExternalKHeap<2, InitialTransferLabel> initialQueue;
    std::vector<ProfileDijkstraLabel> dijkstraLabels;
    ExternalKHeap<2, ProfileDijkstraLabel> queue;

    Profiler profiler;
};
}
