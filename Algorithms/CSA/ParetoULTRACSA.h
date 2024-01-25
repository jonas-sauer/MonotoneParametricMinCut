#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

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

template<bool PATH_RETRIEVAL = true, typename PROFILER = NoProfiler, size_t MAX_TRIPS = 8>
class ParetoULTRACSA {

public:
    using InitialTransferGraph = CHGraph;
    constexpr static bool PathRetrieval = PATH_RETRIEVAL;
    using Profiler = PROFILER;
    constexpr static size_t MaxTrips = MAX_TRIPS;
    using Type = ParetoULTRACSA<PathRetrieval, Profiler, MaxTrips>;
    using TripFlag = Meta::IF<PathRetrieval, ConnectionId, bool>;

private:
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
    ParetoULTRACSA(const Data& data, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialTransfers(chData, FORWARD, data.numberOfStops()),
        sourceVertex(noVertex),
        sourceDepartureTime(never),
        targetVertex(noVertex),
        tripReached(data.numberOfTrips() * MaxTrips, TripFlag()),
        arrivalTime((data.numberOfStops() + 1) * MaxTrips, never),
        parentLabel(PathRetrieval ? (data.numberOfStops() + 1) * MaxTrips : 0),
        profiler(profilerTemplate) {
        AssertMsg(!Graph::hasLoops(data.transferGraph), "Shortcut graph may not have loops!");
        AssertMsg(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline void run(const Vertex source, const int departureTime, const Vertex target = noVertex) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        targetVertex = target;
        if (target == noVertex) {
            targetStop = noStop;
        } else {
            targetStop = data.isStop(target) ? StopId(target) : StopId(data.numberOfStops());
        }
        if (data.isStop(sourceVertex)) {
            arrivalTime[sourceVertex * MaxTrips] = departureTime;
        }
        runInitialTransfers();
        const ConnectionId firstConnection = firstReachableConnection(departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections(firstConnection, ConnectionId(data.connections.size()));
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.done();
    }

    inline bool reachable(const Vertex vertex) const noexcept {
        const StopId stop = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return getEarliestArrivalTime(stop) < never;
    }

    inline int getArrivalTime(const Vertex vertex, const size_t numTrips) const noexcept {
        const StopId stop = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return arrivalTime[stop * MaxTrips + numTrips];
    }

    inline int getEarliestArrivalTime(const Vertex vertex) const noexcept {
        const StopId stop = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        for (size_t i = (stop + 1) * MaxTrips - 1; i != size_t(stop * MaxTrips - 1); i--) {
            if (arrivalTime[i] < never) return arrivalTime[i];
        }
        return never;
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline std::vector<Journey> getJourneys() const noexcept {
        return getJourneys(targetStop);
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline std::vector<Journey> getJourneys(const Vertex vertex) const noexcept {
        std::vector<Journey> journeys;
        for (size_t i = 0; i < MaxTrips; i++) {
            getJourney(journeys, i, vertex);
        }
        return journeys;
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline Journey getEarliestJourney(const Vertex vertex) const noexcept {
        std::vector<Journey> journeys = getJourneys(vertex);
        return journeys.empty() ? Journey() : journeys.back();
    }

    inline std::vector<int> getArrivalTimes() const noexcept {
        return getArrivalTimes(targetStop);
    }

    inline std::vector<int> getArrivalTimes(const Vertex vertex) const noexcept {
        const StopId stop = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<int> arrivalTimes;
        for (size_t i = 0; i < MaxTrips; i++) {
            arrivalTimes.emplace_back(std::min(arrivalTime[stop * MaxTrips + i], (arrivalTimes.empty()) ? (never) : (arrivalTimes.back())));
        }
        return arrivalTimes;
    }

    inline std::vector<Vertex> getPath(const Vertex vertex) noexcept {
        return journeyToPath(getJourneys(vertex).back());
    }

    inline std::vector<std::string> getRouteDescription(const Vertex vertex) noexcept {
        return data.journeyToText(getJourneys(vertex).back());
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() {
        sourceVertex = noVertex;
        sourceDepartureTime = never;
        targetVertex = noVertex;
        targetStop = noStop;
        Vector::fill(arrivalTime, never);
        Vector::fill(tripReached, TripFlag());
        if constexpr (PathRetrieval) {
            Vector::fill(parentLabel, ParentLabel());
        }
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
                if (targetStop != noStop && connection.departureTime > arrivalTime[targetStop * MaxTrips + numTrips]) break;
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

        if (initialTransfers.getBackwardDistance(stop) != INFTY) {
            profiler.countMetric(METRIC_EDGES);
            const int newArrivalTime = time + initialTransfers.getBackwardDistance(stop);
            arrivalByTransfer(targetStop, newArrivalTime, numTrips, stop, noEdge);
        }
    }

    inline void runInitialTransfers() noexcept {
        initialTransfers.run(sourceVertex, targetVertex);
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            AssertMsg(data.isStop(stop), "Reached POI " << stop << " is not a stop!");
            AssertMsg(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            profiler.countMetric(METRIC_EDGES);
            const int newArrivalTime = sourceDepartureTime + initialTransfers.getForwardDistance(stop);
            arrivalByTransfer(StopId(stop), newArrivalTime, 0, sourceVertex, noEdge);
        }
        if (initialTransfers.getDistance() != INFTY) {
            const int newArrivalTime = sourceDepartureTime + initialTransfers.getDistance();
            profiler.countMetric(METRIC_EDGES);
            arrivalByTransfer(targetStop, newArrivalTime, 0, sourceVertex, noEdge);
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
        }
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline void getJourney(std::vector<Journey>& journeys, size_t numTrips, Vertex vertex) const noexcept {
        StopId stop = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        if (arrivalTime[stop * MaxTrips + numTrips] >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        Journey journey;
        while (stop != sourceVertex) {
            const ParentLabel& label = parentLabel[stop * MaxTrips + numTrips];
            if (label.reachedByTransfer) {
                const int parentDepartureTime = (label.parent == sourceVertex) ? sourceDepartureTime : arrivalTime[label.parent * MaxTrips + numTrips];
                journey.emplace_back(label.parent, stop, parentDepartureTime, arrivalTime[stop * MaxTrips + numTrips], label.transferId);
            } else {
                journey.emplace_back(label.parent, stop, data.connections[tripReached[label.tripId]].departureTime, arrivalTime[stop * MaxTrips + numTrips], label.tripId);
                numTrips--;
            }
            stop = StopId(label.parent);
        }
        journeys.emplace_back(Vector::reverse(journey));
    }

private:
    const Data& data;
    RAPTOR::BucketCHInitialTransfers initialTransfers;

    Vertex sourceVertex;
    int sourceDepartureTime;
    Vertex targetVertex;
    StopId targetStop;

    std::vector<TripFlag> tripReached;
    std::vector<int> arrivalTime;
    std::vector<ParentLabel> parentLabel;

    Profiler profiler;

};
}
