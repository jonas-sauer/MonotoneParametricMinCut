#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../../Helpers/Vector/Vector.h"

#include "../InitialTransfers.h"

#include "../../../DataStructures/RAPTOR/Data.h"
#include "../../../DataStructures/RAPTOR/Entities/Journey.h"
#include "../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../DataStructures/Intermediate/Data.h"
#include "../../../DataStructures/Container/Set.h"
#include "../../../DataStructures/Container/Map.h"
#include "../../../DataStructures/Container/ExternalKHeap.h"

namespace RAPTOR::RangeRAPTOR {

template<typename PROFILER, bool COLLECT_REACHED_VERTICES = false, bool TRIP_PRUNING = true>
class ULTRARAPTORModule {

public:
    using InitialTransferGraph = CHGraph;
    using Profiler = PROFILER;
    static constexpr bool CollectReachedVertices = COLLECT_REACHED_VERTICES;
    static constexpr bool TripPruning = TRIP_PRUNING;
    using Type = ULTRARAPTORModule<Profiler, CollectReachedVertices, TripPruning>;

public:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parentDepartureTime(never), parent(noVertex), usesRoute(false), routeId(noRouteId), timestamp(0) {}
        int arrivalTime;
        int parentDepartureTime;
        Vertex parent;
        bool usesRoute;
        RouteId routeId;
        int timestamp;
    };
    using Round = std::vector<EarliestArrivalLabel>;

public:
    template<typename ATTRIBUTE>
    ULTRARAPTORModule(const Data& data, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialTransfers(forwardGraph, backwardGraph, data.numberOfStops(), weight),
        roundIndex(-1),
        stopsUpdatedByRoute(data.numberOfStops() + 1),
        stopsUpdatedByTransfer(data.numberOfStops() + 1),
        routesServingUpdatedStops(data.numberOfRoutes()),
        sourceVertex(noVertex),
        targetVertex(noVertex),
        targetStop(noStop),
        sourceDepartureTime(intMax),
        timestamp(0),
        profiler(profilerTemplate),
        reachedVertices(data.numberOfStops() + 1) {
        Assert(data.hasImplicitBufferTimes(), "Departure buffer times have to be implicit!");
        if constexpr (TripPruning) resetLastTrip();
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    ULTRARAPTORModule(const Data& data, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        ULTRARAPTORModule(data, chData.forward, chData.backward, Weight, profilerTemplate) {
    }

    template<bool CLEAR = true>
    inline void run(const Vertex source, const int departureTime, const Vertex target = noVertex, const size_t maxRounds = INFTY) noexcept {
        runInitialize<CLEAR>(source, departureTime, target);
        runInitialTransfers();
        evaluateInitialTransfers();
        runRounds(maxRounds);
    }

    template<bool CLEAR = true>
    inline void runInitialize(const Vertex source, const int departureTime, const Vertex target = noVertex) noexcept {
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear<CLEAR>();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime, target);
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.doneRound();
    }

    inline void runInitialTransfers() noexcept {
        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        relaxInitialTransfers();
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.doneRound();
    }

    inline void runAddSource(const StopId source, const int departureTime) noexcept {
        EarliestArrivalLabel& label = currentRoundLabel(source);
        label.arrivalTime = departureTime;
        label.parentDepartureTime = sourceDepartureTime;
        label.usesRoute = false;
        stopsUpdatedByTransfer.insert(source);
    }

    inline void runRounds(const size_t maxRounds = INFTY) noexcept {
        profiler.startRound();
        profiler.startPhase();
        startNewRound();
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.startPhase();
        collectRoutesServingUpdatedStops();
        profiler.donePhase(PHASE_COLLECT);
        profiler.startPhase();
        scanRoutes();
        profiler.donePhase(PHASE_SCAN);
        for (size_t i = 0; (i < maxRounds) && (!stopsUpdatedByRoute.empty()); i++) {
            profiler.startPhase();
            relaxIntermediateTransfers();
            profiler.donePhase(PHASE_TRANSFERS);
            profiler.doneRound();
            profiler.startRound();
            profiler.startPhase();
            startNewRound();
            profiler.donePhase(PHASE_INITIALIZATION);
            profiler.startPhase();
            collectRoutesServingUpdatedStops();
            profiler.donePhase(PHASE_COLLECT);
            profiler.startPhase();
            scanRoutes();
            profiler.donePhase(PHASE_SCAN);
        }
        profiler.doneRound();
        profiler.done();
    }

    inline std::vector<Journey> getJourneys() noexcept {
        return getJourneys(targetStop);
    }

    inline std::vector<Journey> getJourneys(const Vertex vertex) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<Journey> journeys;
        for (size_t i = 0; i <= roundIndex; i++) {
            getJourney(journeys, i, target);
        }
        return journeys;
    }

    inline bool getJourney(const Vertex vertex, const size_t round, Journey& journey) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<Journey> journeys;
        getJourney(journeys, round, target);
        if (journeys.empty()) return false;
        journey = journeys[0];
        return true;
    }

    inline Journey getEarliestJourney(const Vertex vertex) noexcept {
        std::vector<Journey> journeys = getJourneys(vertex);
        return journeys.empty() ? Journey() : journeys.back();
    }

    inline std::vector<ArrivalLabel> getArrivals() noexcept {
        return getArrivals(targetStop);
    }

    inline std::vector<ArrivalLabel> getArrivals(const Vertex vertex) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<ArrivalLabel> labels;
        for (size_t i = 0; i <= roundIndex; i++) {
            getArrival(labels, i, target);
        }
        return labels;
    }

    inline bool reachable(const Vertex vertex) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return getEarliestArrivalTime(target) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return roundLabel(rounds.size() - 1, target).arrivalTime;
    }

    inline int getWalkingArrivalTime() const noexcept {
        return sourceDepartureTime + initialTransfers.getDistance();
    }

    inline int getWalkingArrivalTime(const Vertex vertex) noexcept {
        return sourceDepartureTime + initialTransfers.getForwardDistance(vertex);
    }

    inline int getWalkingTravelTime() const noexcept {
        return initialTransfers.getDistance();
    }

    inline int getWalkingTravelTime(const Vertex vertex) noexcept {
        return initialTransfers.getForwardDistance(vertex);
    }

    inline std::vector<Vertex> getPath(const Vertex vertex) {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return journeyToPath(getJourneys(target).back());
    }

    inline std::vector<std::string> getRouteDescription(const Vertex vertex) {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return data.journeyToText(getJourneys(target).back());
    }

    inline int getArrivalTime(const Vertex vertex, const size_t numberOfTrips) noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        Assert(roundLabel(numberOfTrips, target).arrivalTime < never, "No label found for stop " << target << " in round " << numberOfTrips << "!");
        return roundLabel(numberOfTrips, target).arrivalTime;
    }

    inline const IndexedSet<false, Vertex>& getReachedVertices() const noexcept {
        return reachedVertices;
    }

    inline bool targetReached() const noexcept {
        return reachedVertices.contains(targetStop);
    }

    template<bool RESET = true>
    inline void clear() noexcept {
        roundIndex = -1;
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        sourceVertex = noVertex;
        targetVertex = noVertex;
        targetStop = StopId(data.numberOfStops());
        reachedVertices.clear();
        if constexpr (RESET) {
            if constexpr (TripPruning) resetLastTrip();
            timestamp = 1;
            std::vector<Round>().swap(rounds);
        } else {
            timestamp++;
        }
    }

    inline void reset() noexcept {
        clear<true>();
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    template<bool T = TripPruning, typename = std::enable_if_t<T == TripPruning && T>>
    inline void resetLastTrip() noexcept {
        std::vector<std::vector<std::vector<uint16_t>>>(1).swap(lastTrip);
        std::vector<std::vector<uint16_t>>(data.numberOfRoutes()).swap(lastTrip[0]);
        for (const RouteId route : data.routes()) {
            std::vector<uint16_t>(data.numberOfStopsInRoute(route) - 1, data.numberOfTripsInRoute(route) - 1).swap(lastTrip[0][route]);
        }
    }

    inline void initialize(const Vertex source, const int departureTime, const Vertex target) noexcept {
        sourceVertex = source;
        targetVertex = target;
        if (data.isStop(target)) {
            targetStop = StopId(target);
        }
        sourceDepartureTime = departureTime;
        startNewRound();
        if (data.isStop(source)) {
            stopsUpdatedByRoute.insert(StopId(source));
            EarliestArrivalLabel& label = currentRoundLabel(StopId(source));
            label.arrivalTime = departureTime;
            label.parent = source;
            label.parentDepartureTime = departureTime;
            label.usesRoute = false;
        }
    }

    inline void evaluateInitialTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            if (stop == targetStop) continue;
            Assert(data.isStop(stop), "Reached POI " << stop << " is not a stop!");
            Assert(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            //The initial transfers are evaluated automatically when the label is updated
            roundLabel(0, StopId(stop));
            stopsUpdatedByTransfer.insert(StopId(stop));
            if constexpr (CollectReachedVertices) reachedVertices.insert(stop);
        }
        if (initialTransfers.getDistance() != INFTY) {
            //The initial transfers are evaluated automatically when the label is updated
            roundLabel(0, targetStop);
            if constexpr (CollectReachedVertices) reachedVertices.insert(targetStop);
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            for (const RouteSegment& route : data.routesContainingStop(stop)) {
                Assert(data.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                Assert(data.stopIds[data.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == data.numberOfStopsInRoute(route.routeId)) continue;
                if constexpr (TripPruning) {
                    const uint16_t tripNum = lastTrip[roundIndex][route.routeId][route.stopIndex];
                    if (tripNum == uint16_t(-1)) continue;
                    const StopEvent* trip = data.tripOfRoute(route.routeId, tripNum) + route.stopIndex;
                    if (trip->departureTime < previousRoundLabel(stop).arrivalTime) continue;
                } else {
                    if (data.lastTripOfRoute(route.routeId)[route.stopIndex].departureTime < previousRoundLabel(stop).arrivalTime) continue;
                }
                if (routesServingUpdatedStops.contains(route.routeId)) {
                    routesServingUpdatedStops[route.routeId] = std::min(routesServingUpdatedStops[route.routeId], route.stopIndex);
                } else {
                    routesServingUpdatedStops.insert(route.routeId, route.stopIndex);
                }
            }
        }
    }

    inline void scanRoutes() noexcept {
        stopsUpdatedByRoute.clear();
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            profiler.countMetric(METRIC_ROUTES);
            StopIndex stopIndex = routesServingUpdatedStops[route];
            const size_t tripSize = data.numberOfStopsInRoute(route);
            Assert(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ", RoundIndex: " << roundIndex << ")!");

            const StopId* stops = data.stopArrayOfRoute(route);
            const StopEvent* firstTrip = data.firstTripOfRoute(route);
            const StopEvent* trip;
            if constexpr (TripPruning) {
                trip = data.tripOfRoute(route, lastTrip[roundIndex][route][stopIndex]);
                if (trip < firstTrip) continue;
            } else {
                trip = data.lastTripOfRoute(route);
            }
            StopId stop = stops[stopIndex];
            Assert(trip[stopIndex].departureTime >= previousRoundLabel(stop).arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousRoundLabel(stop).arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ", RoundIndex: " << roundIndex << ")!");

            StopIndex parentIndex = stopIndex;
            while (stopIndex < tripSize - 1) {
                while ((trip > firstTrip) && ((trip - tripSize + stopIndex)->departureTime >= previousRoundLabel(stop).arrivalTime)) {
                    trip -= tripSize;
                    parentIndex = stopIndex;
                }
                if constexpr (TripPruning) lastTrip[roundIndex][route][stopIndex] = uint16_t(((trip - firstTrip) / tripSize) - 1);
                stopIndex++;
                stop = stops[stopIndex];
                profiler.countMetric(METRIC_ROUTE_SEGMENTS);
                if (arrivalByRoute(stop, trip[stopIndex].arrivalTime)) {
                    if constexpr (CollectReachedVertices) reachedVertices.insert(stop);
                    EarliestArrivalLabel& label = currentRoundLabel(stop);
                    label.parent = stops[parentIndex];
                    label.parentDepartureTime = trip[parentIndex].departureTime;
                    label.usesRoute = true;
                    label.routeId = route;
                }
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        initialTransfers.run(sourceVertex, targetVertex);
    }

    inline void relaxIntermediateTransfers() noexcept {
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        for (const StopId stop : stopsUpdatedByRoute) {
            const int earliestArrivalTime = currentRoundLabel(stop).arrivalTime;
            for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
                const StopId toStop = StopId(data.transferGraph.get(ToVertex, edge));
                if (toStop == targetStop) continue;
                profiler.countMetric(METRIC_EDGES);
                const int arrivalTime = earliestArrivalTime + data.transferGraph.get(TravelTime, edge);
                Assert(data.isStop(toStop), "Graph contains edges to non stop vertices!");
                arrivalByTransfer(toStop, arrivalTime, stop, earliestArrivalTime);
            }
            if (initialTransfers.getBackwardDistance(stop) != INFTY) {
                const int arrivalTime = earliestArrivalTime + initialTransfers.getBackwardDistance(stop);
                arrivalByTransfer(targetStop, arrivalTime, stop, earliestArrivalTime);
            }
            stopsUpdatedByTransfer.insert(stop);
        }
    }

    inline EarliestArrivalLabel& roundLabel(const size_t round, const StopId stop) noexcept {
        EarliestArrivalLabel& result = rounds[round][stop];
        if (result.timestamp != timestamp) {
            if (round > 0) {
                result.arrivalTime = std::min(result.arrivalTime, roundLabel(round - 1, stop).arrivalTime);
            } else {
                const int distance = (stop == targetStop) ? initialTransfers.getDistance() : initialTransfers.getForwardDistance(stop);
                if (distance != INFTY) {
                    profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
                    result.arrivalTime = sourceDepartureTime + distance;
                    result.parent = sourceVertex;
                    result.parentDepartureTime = sourceDepartureTime;
                    result.usesRoute = false;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline EarliestArrivalLabel& currentRoundLabel(const StopId stop) noexcept {
        Assert(roundIndex < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        return roundLabel(roundIndex, stop);
    }

    inline EarliestArrivalLabel& previousRoundLabel(const StopId stop) noexcept {
        Assert(roundIndex - 1 < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        Assert(roundIndex > 0, "Cannot return previous round, because no round exists!");
        return roundLabel(roundIndex - 1, stop);
    }

    inline void startNewRound() noexcept {
        Assert(roundIndex + 1 <= rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        roundIndex++;
        if (roundIndex == rounds.size()) {
            if (rounds.empty()) {
                rounds.emplace_back(data.numberOfStops() + 1);
            } else {
                rounds.emplace_back(rounds.back());
                if constexpr (TripPruning) lastTrip.emplace_back(lastTrip.back());
            }
        }
    }

    inline bool arrivalByRoute(const StopId stop, const int arrivalTime) noexcept {
        Assert(roundIndex > 0, "arrivalByRoute cannot be used in the first round!");
        Assert(data.isStop(stop), "Stop " << stop << " is out of range!");
        Assert(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], stop: " << stop << ")!");
        if ((currentRoundLabel(targetStop).arrivalTime <= arrivalTime) || (previousRoundLabel(targetStop).arrivalTime <= arrivalTime)) return false;
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if ((label.arrivalTime <= arrivalTime) || (previousRoundLabel(stop).arrivalTime <= arrivalTime)) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        label.arrivalTime = arrivalTime;
        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    inline void arrivalByTransfer(const StopId stop, const int arrivalTime, const Vertex parent, const int parentDepartureTime) noexcept {
        Assert(data.isStop(stop) || stop == targetStop, "Stop " << stop << " is out of range!");
        Assert(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if (data.isStop(stop)) stopsUpdatedByTransfer.insert(stop);
        if (label.arrivalTime <= arrivalTime) return;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        label.arrivalTime = arrivalTime;
        label.parent = parent;
        label.parentDepartureTime = parentDepartureTime;
        label.usesRoute = false;
        if constexpr (CollectReachedVertices) reachedVertices.insert(stop);
    }

    inline void getJourney(std::vector<Journey>& journeys, size_t round, Vertex vertex) noexcept {
        if (roundLabel(round, StopId(vertex)).arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        Journey journey;
        do {
            Assert(round != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = roundLabel(round, StopId(vertex));
            journey.emplace_back(label.parent, vertex, label.parentDepartureTime, label.arrivalTime, label.usesRoute, label.routeId);
            Assert(data.isStop(label.parent) || label.parent == sourceVertex, "Backtracking parent pointers reached a vertex (" << label.parent << ")!");
            vertex = label.parent;
            if (label.usesRoute) round--;
        } while (journey.back().from != sourceVertex);
        journeys.emplace_back(Vector::reverse(journey));
    }

    inline void getArrival(std::vector<ArrivalLabel>& labels, size_t round, const Vertex vertex) noexcept {
        const EarliestArrivalLabel& label = roundLabel(round, StopId(vertex));
        if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
        labels.emplace_back(label.arrivalTime, round);
    }

    inline void printRoundsForStop(const StopId stop) noexcept {
        Assert(data.isStop(stop), stop << " is not a valid stop!");
        std::cout << "Raptor Label for stop " << stop << ":" << std::endl;
        std::cout << std::setw(10) << "Round" << std::setw(14) << "arrivalTime" << std::setw(14) << "parent" << std::endl;
        for (size_t i = 0; i <= roundIndex; i++) {
            EarliestArrivalLabel& label = roundLabel(i, stop);
            std::cout << std::setw(10) << i << std::setw(14) << label.arrivalTime << std::setw(14) << label.parent << std::endl;
        }
    }

private:
    const Data& data;

    BucketCHInitialTransfers initialTransfers;

    std::vector<Round> rounds;
    size_t roundIndex;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    Vertex targetVertex;
    StopId targetStop;
    int sourceDepartureTime;
    int timestamp;

    std::vector<std::vector<std::vector<uint16_t>>> lastTrip;

    Profiler profiler;

    IndexedSet<false, Vertex> reachedVertices;

};

}
