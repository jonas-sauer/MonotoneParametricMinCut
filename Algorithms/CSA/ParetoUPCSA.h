#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <type_traits>
#include <concepts>

#include "../CH/CH.h"
#include "../RAPTOR/InitialTransfers.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Journey.h"
#include "Profiler.h"

namespace CSA {

template<bool USE_DFS_ORDER, bool PATH_RETRIEVAL = true, typename PROFILER = NoProfiler, size_t MAX_TRIPS = 8>
class ParetoUPCSA {

public:
    constexpr static bool UseDFSOrder = USE_DFS_ORDER;
    constexpr static bool PathRetrieval = PATH_RETRIEVAL;
    using Profiler = PROFILER;
    constexpr static bool Debug = std::is_same_v<Profiler, SimpleProfiler>;
    constexpr static size_t MaxTrips = MAX_TRIPS;
    using Type = ParetoUPCSA<UseDFSOrder, PathRetrieval, Profiler, MaxTrips>;
    using InitialAndFinalTransfers = RAPTOR::GroupedParetoInitialAndFinalTransfers<Debug, MaxTrips>;
    using TripFlag = std::conditional_t<PathRetrieval, ConnectionId, bool>;

private:
    struct QueryData {
        QueryData(const Data& oldData, const CH::CH& oldCHData, const IndexedSet<false, Vertex>& oldTargets, const bool reorder) :
            data(oldData),
            chData(oldCHData),
            targets(oldTargets) {
            const Order chOrder = vertexOrder(chData);
            if (reorder) {
                const Order fullOrder = chOrder.splitAt(data.numberOfStops());
                Order stopOrder = fullOrder;
                stopOrder.resize(data.numberOfStops());

                data.applyStopOrder(stopOrder);
                chData.applyVertexOrder(fullOrder);
                targets.applyPermutation(Permutation(Construct::Invert, fullOrder));

                size_t numStops = 0;
                size_t numVertices = data.numberOfStops();
                for (const size_t i : chOrder) {
                    if (i < data.numberOfStops()) {
                        phastOrder.emplace_back(numStops++);
                    } else {
                        phastOrder.emplace_back(numVertices++);
                    }
                }
            } else {
                phastOrder = chOrder;
            }
        }

        Data data;
        CH::CH chData;
        Order phastOrder;
        IndexedSet<false, Vertex> targets;
    };

    struct ParentLabel {
        ParentLabel(const Vertex parent = noVertex, const bool reachedByTransfer = false, const TripId tripId = noTripId) :
            parent(parent),
            reachedByTransfer(reachedByTransfer),
            tripId(tripId) {
        }

        Vertex parent;
        bool reachedByTransfer;
        union {
            TripId tripId;
            Edge transferId;
        };
    };

public:
    ParetoUPCSA(const Data& oldData, const CH::CH& oldCHData, const IndexedSet<false, Vertex>& oldTargets, const bool reorder, const Profiler& profilerTemplate = Profiler()) :
        queryData(oldData, oldCHData, oldTargets, reorder),
        data(queryData.data),
        initialAndFinalTransfers(data.transferGraph, queryData.chData, std::move(queryData.phastOrder), data.numberOfStops(), queryData.targets),
        sourceVertex(noVertex),
        sourceDepartureTime(never),
        targetVertices(queryData.targets),
        targetStopId(queryData.chData.numVertices(), noStop),
        tripReached(data.numberOfTrips() * MaxTrips, TripFlag()),
        arrivalTime(computeTargetMapping() * MaxTrips, never),
        parentLabel(PathRetrieval ? arrivalTime.size() : 0),
        finalTransferSet(arrivalTime.size(), false),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        Assert(!Graph::hasLoops(data.transferGraph), "Shortcut graph may not have loops!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN, PHASE_UPWARD_SWEEP, PHASE_DOWNWARD_SEARCH});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline static Order vertexOrder(const CH::CH& chData) noexcept {
        if constexpr (UseDFSOrder) {
            return Order(Vector::reverse(CH::getOrder(chData)));
        } else {
            return Order(CH::getLevelOrderTopDown(chData));
        }
    }

    inline void run(const Vertex source, const int departureTime) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        initialAndFinalTransfers.initialize();
        if (data.isStop(source)) {
            arrivalTime[source * MaxTrips] = departureTime;
        }
        runInitialTransfers();
        const ConnectionId firstConnection = firstReachableConnection(departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections(firstConnection, ConnectionId(data.connections.size()));
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.startPhase();
        initialAndFinalTransfers.upwardSweep();
        profiler.donePhase(PHASE_UPWARD_SWEEP);
        profiler.startPhase();
        initialAndFinalTransfers.downwardSearchToTargets();
        profiler.donePhase(PHASE_DOWNWARD_SEARCH);

        profiler.done();
    }

    inline bool reachable(const Vertex vertex) noexcept {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        return getEarliestArrivalTime(vertex) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) noexcept {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        const StopId stop = targetStopId[vertex];
        for (size_t numTrips = MaxTrips - 1; numTrips != size_t(-1); numTrips--) {
            checkFinalTransfer(vertex, numTrips);
            const size_t i = stop * MaxTrips + numTrips;
            if (arrivalTime[i] < never) return arrivalTime[i];
        }
        return never;
    }

    inline std::vector<Journey> getJourneys(const Vertex vertex) noexcept requires PathRetrieval {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        std::vector<Journey> journeys;
        for (size_t i = 0; i < MaxTrips; i++) {
            getJourney(journeys, i, vertex);
        }
        return journeys;
    }

    inline Journey getEarliestJourney(const Vertex vertex) noexcept requires PathRetrieval {
        std::vector<Journey> journeys = getJourneys(vertex);
        return journeys.empty() ? Journey() : journeys.back();
    }

    inline std::vector<int> getArrivalTimes(const Vertex vertex) noexcept {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        const StopId stop = targetStopId[vertex];
        std::vector<int> arrivalTimes;
        for (size_t i = 0; i < MaxTrips; i++) {
            checkFinalTransfer(vertex, i);
            arrivalTimes.emplace_back(std::min(arrivalTime[stop * MaxTrips + i], (arrivalTimes.empty()) ? (never) : (arrivalTimes.back())));
        }
        return arrivalTimes;
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        sourceVertex = noVertex;
        sourceDepartureTime = never;
        Vector::fill(arrivalTime, never);
        Vector::fill(tripReached, TripFlag());
        if constexpr (PathRetrieval) {
            Vector::fill(parentLabel, ParentLabel());
        }
        Vector::fill(finalTransferSet, false);
    }

    inline ConnectionId firstReachableConnection(const int departureTime) const noexcept {
        return ConnectionId(Vector::lowerBound(data.connections, departureTime, [](const Connection& connection, const int time) {
            return connection.departureTime < time;
        }));
    }

    inline void scanConnections(const ConnectionId begin, const ConnectionId end) noexcept {
        for (ConnectionId i = begin; i < end; i++) {
            const Connection& connection = data.connections[i];
            for (size_t numTrips = 1; numTrips < MaxTrips; numTrips++) {
                if (connectionIsReachableFromTrip(connection, numTrips)) {
                    profiler.countMetric(METRIC_CONNECTIONS);
                    arrivalByTrip(connection.arrivalStopId, connection.arrivalTime, numTrips, connection.tripId);
                    break;
                }
                if (connectionIsReachableFromStop(connection, numTrips - 1)) {
                    profiler.countMetric(METRIC_CONNECTIONS);
                    useTrip(connection.tripId, i, numTrips);
                    arrivalByTrip(connection.arrivalStopId, connection.arrivalTime, numTrips, connection.tripId);
                    break;
                }
            }
        }
    }

    inline bool connectionIsReachableFromStop(const Connection& connection, const size_t numTrips) const noexcept {
        return arrivalTime[connection.departureStopId * MaxTrips + numTrips] <= connection.departureTime - data.minTransferTime(connection.departureStopId);
    }

    inline bool connectionIsReachableFromTrip(const Connection& connection, const size_t numTrips) const noexcept {
        return tripReached[connection.tripId * MaxTrips + numTrips] != TripFlag();
    }

    inline void useTrip(const TripId tripId, const ConnectionId connectionId, const size_t numTrips) noexcept {
        const size_t i = tripId * MaxTrips + numTrips;
        if constexpr (PathRetrieval) {
            tripReached[i] = connectionId;
        } else {
            suppressUnusedParameterWarning(connectionId);
            tripReached[i] = true;
        }
    }

    inline void arrivalByTrip(const StopId stop, const int time, const size_t numTrips, const TripId trip) noexcept {
        const size_t i = stop * MaxTrips + numTrips;
        if (arrivalTime[i] <= time) return;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        arrivalTime[i] = time;
        if constexpr (PathRetrieval) {
            parentLabel[i].parent = data.connections[tripReached[trip * MaxTrips + numTrips]].departureStopId;
            parentLabel[i].reachedByTransfer = false;
            parentLabel[i].tripId = trip;
        }

        for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
            profiler.countMetric(METRIC_EDGES);
            const StopId toStop = StopId(data.transferGraph.get(ToVertex, edge));
            const int newArrivalTime = time + data.transferGraph.get(TravelTime, edge);
            arrivalByTransfer(toStop, newArrivalTime, numTrips, stop, edge);
        }

        initialAndFinalTransfers.template addSource<true>(stop, time - sourceDepartureTime, stop, numTrips);
    }

    inline void runInitialTransfers() noexcept {
        initialAndFinalTransfers.startNewRound();
        initialAndFinalTransfers.template addSource<false>(sourceVertex, 0, sourceVertex, 0);
        initialAndFinalTransfers.initialUpwardSearch();
        initialAndFinalTransfers.downwardSearchToStops();
        for (const StopId stop : data.stops()) {
            const int newArrivalTime = sourceDepartureTime + initialAndFinalTransfers.getDistance(0, stop);
            arrivalByTransfer(StopId(stop), newArrivalTime, 0, sourceVertex, noEdge);
        }
    }

    inline void arrivalByTransfer(const StopId stop, const int time, const size_t numTrips, const Vertex parent, const Edge edge) noexcept {
        const size_t i = stop * MaxTrips + numTrips;
        if (arrivalTime[i] <= time) return;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        arrivalTime[i] = time;
        if constexpr (PathRetrieval) {
            parentLabel[i].parent = parent;
            parentLabel[i].reachedByTransfer = true;
            parentLabel[i].transferId = edge;
        } else {
            suppressUnusedParameterWarning(parent);
            suppressUnusedParameterWarning(edge);
        }
    }

    inline size_t computeTargetMapping() noexcept {
        StopId stopId(data.numberOfStops());
        for (const Vertex target : targetVertices) {
            targetStopId[target] = data.isStop(target) ? StopId(target) : stopId;
            stopId++;
        }
        return stopId;
    }

    inline void checkFinalTransfer(const Vertex vertex, const size_t numTrips) noexcept {
        if (!targetVertices.contains(vertex)) return;
        const StopId targetStop = targetStopId[vertex];
        const size_t i = targetStop * MaxTrips + numTrips;
        if (finalTransferSet[i]) return;

        const int distance = initialAndFinalTransfers.getDistance(numTrips, vertex);
        if (distance == INFTY) return;

        const int finalTransferArrivalTime = sourceDepartureTime + distance;
        if (finalTransferArrivalTime >= arrivalTime[i]) return;

        arrivalTime[i] = finalTransferArrivalTime;
        if constexpr (PathRetrieval) {
            parentLabel[i].parent = initialAndFinalTransfers.getParent(numTrips, vertex);
            parentLabel[i].reachedByTransfer = true;
            parentLabel[i].transferId = noEdge;

        }
        finalTransferSet[i] = true;
    }

    inline void getJourney(std::vector<Journey>& journeys, size_t numTrips, Vertex vertex) noexcept requires PathRetrieval {
        StopId stop = targetStopId[vertex];
        checkFinalTransfer(vertex, numTrips);
        if (arrivalTime[stop * MaxTrips + numTrips] >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        Journey journey;
        while (stop != sourceVertex) {
            checkFinalTransfer(stop, numTrips);
            const ParentLabel& label = parentLabel[stop * MaxTrips + numTrips];
            if (label.reachedByTransfer) {
                const int parentDepartureTime = (label.parent == sourceVertex) ? sourceDepartureTime : arrivalTime[label.parent * MaxTrips + numTrips];
                journey.emplace_back(label.parent, stop, parentDepartureTime, arrivalTime[stop * MaxTrips + numTrips], label.transferId);
            } else {
                journey.emplace_back(label.parent, stop, data.connections[tripReached[label.tripId * MaxTrips + numTrips]].departureTime, arrivalTime[stop * MaxTrips + numTrips], label.tripId);
                numTrips--;
            }
            stop = StopId(label.parent);
        }
        journeys.emplace_back(Vector::reverse(journey));
    }

private:
    const QueryData queryData;
    const Data& data;

    InitialAndFinalTransfers initialAndFinalTransfers;

    Vertex sourceVertex;
    int sourceDepartureTime;
    const IndexedSet<false, Vertex>& targetVertices;
    std::vector<StopId> targetStopId;

    std::vector<TripFlag> tripReached;
    std::vector<int> arrivalTime;
    std::vector<ParentLabel> parentLabel;
    std::vector<bool> finalTransferSet;

    Profiler profiler;

};
}
