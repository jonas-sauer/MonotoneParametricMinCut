#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <concepts>

#include "../../../Helpers/Vector/Vector.h"

#include "../InitialTransfers.h"
#include "../Profiler.h"

#include "../../../DataStructures/RAPTOR/Data.h"
#include "../../../DataStructures/RAPTOR/Entities/Journey.h"
#include "../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../DataStructures/Intermediate/Data.h"
#include "../../../DataStructures/Container/Set.h"
#include "../../../DataStructures/Container/Map.h"
#include "../../../DataStructures/Container/ExternalKHeap.h"

namespace RAPTOR::RangeRAPTOR {

template<typename INITIAL_TRANSFERS, typename PROFILER, bool ONE_TO_ONE = false, bool COLLECT_REACHED_VERTICES = false, bool USE_MIN_TRANSFER_TIMES = false, bool TRIP_PRUNING = true>
class DijkstraRAPTORModule {

public:
    using InitialTransferType = INITIAL_TRANSFERS;
    using InitialTransferGraph = typename InitialTransferType::Graph;
    using Profiler = PROFILER;
    static constexpr bool OneToOne = ONE_TO_ONE;
    static constexpr bool CollectReachedVertices = COLLECT_REACHED_VERTICES;
    static constexpr bool UseMinTransferTimes = USE_MIN_TRANSFER_TIMES;
    static constexpr int RoundFactor = UseMinTransferTimes ? 2 : 1;
    static constexpr size_t InitialTransferRound = UseMinTransferTimes ? 1 : 0;
    static constexpr bool TripPruning = TRIP_PRUNING;
    using Type = DijkstraRAPTORModule<InitialTransferType, Profiler, OneToOne, CollectReachedVertices, UseMinTransferTimes, TripPruning>;

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

    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel() : arrivalTime(never), parent(noVertex), timestamp(0) {}
        int arrivalTime;
        Vertex parent;
        int timestamp;
        inline bool hasSmallerKey(const DijkstraLabel* const other) const noexcept {
            return arrivalTime < other->arrivalTime;
        }
    };

public:
    template<typename ATTRIBUTE>
    DijkstraRAPTORModule(const Data& data, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        initialTransfers(forwardGraph, backwardGraph, data.numberOfStops(), weight),
        roundIndex(-1),
        stopsUpdatedByRoute(data.numberOfStops()),
        stopsUpdatedByTransfer(data.numberOfStops()),
        routesServingUpdatedStops(data.numberOfRoutes()),
        sourceVertex(noVertex),
        targetVertex(noVertex),
        targetStop(noStop),
        sourceDepartureTime(intMax),
        timestamp(0),
        labelByNumberOfTrips(1, std::vector<DijkstraLabel>(data.transferGraph.numVertices())),
        profiler(profilerTemplate),
        reachedVertices(forwardGraph.numVertices()) {
        if constexpr (UseMinTransferTimes) {
            minChangeTimeGraph = data.minChangeTimeGraph();
            Assert(!data.hasImplicitBufferTimes(), "Either min transfer times have to be used OR departure buffer times have to be implicit!");
        } else {
            Assert(data.hasImplicitBufferTimes(), "Either min transfer times have to be used OR departure buffer times have to be implicit!");
        }
        if constexpr (TripPruning) resetLastTrip();
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_VERTICES, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    DijkstraRAPTORModule(const Data& data, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) requires std::same_as<InitialTransferGraph, CHGraph> :
        DijkstraRAPTORModule(data, chData.forward, chData.backward, Weight, profilerTemplate) {
    }

    DijkstraRAPTORModule(const Data& data, const TransferGraph& forwardGraph, const TransferGraph& backwardGraph, const Profiler& profilerTemplate = Profiler()) requires std::same_as<InitialTransferGraph, TransferGraph> :
        DijkstraRAPTORModule(data, forwardGraph, backwardGraph, TravelTime, profilerTemplate) {
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
        dijkstraLabel(roundIndex / RoundFactor, source).arrivalTime = departureTime;
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
            if constexpr (UseMinTransferTimes) {
                profiler.startPhase();
                startNewRound();
                profiler.donePhase(PHASE_INITIALIZATION);
            }
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

    inline std::vector<Journey> getJourneys() noexcept requires OneToOne {
        return getJourneys(targetVertex);
    }

    inline std::vector<Journey> getJourneys(const Vertex vertex) noexcept {
        std::vector<Journey> journeys;
        for (size_t i = 0; i <= roundIndex; i += RoundFactor) {
            getJourney(journeys, i, vertex);
        }
        return journeys;
    }

    inline bool getJourney(const Vertex vertex, const size_t round, Journey& journey) noexcept {
        std::vector<Journey> journeys;
        getJourney(journeys, round, vertex);
        if (journeys.empty()) return false;
        journey = journeys[0];
        return true;
    }

    inline Journey getEarliestJourney(const Vertex vertex) noexcept {
        std::vector<Journey> journeys = getJourneys(vertex);
        return journeys.empty() ? Journey() : journeys.back();
    }

    inline std::vector<ArrivalLabel> getArrivals() noexcept requires OneToOne {
        return getArrivals(targetVertex);
    }

    inline std::vector<ArrivalLabel> getArrivals(const Vertex vertex) noexcept {
        std::vector<ArrivalLabel> labels;
        for (size_t i = 0; i <= roundIndex; i += RoundFactor) {
            getArrival(labels, i, vertex);
        }
        return labels;
    }

    inline bool reachable(const Vertex vertex) noexcept {
        Assert(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return getEarliestArrivalTime(vertex) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) noexcept {
        Assert(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                return std::min(roundLabel(rounds.size() - 1, StopId(vertex)).arrivalTime, roundLabel(rounds.size() - 2, StopId(vertex)).arrivalTime);
            } else {
                return roundLabel(rounds.size() - 1, StopId(vertex)).arrivalTime;
            }
        } else {
            return dijkstraLabel(labelByNumberOfTrips.size() - 1, vertex).arrivalTime;
        }
    }

    inline int getWalkingArrivalTime() const noexcept requires OneToOne {
        return sourceDepartureTime + initialTransfers.getDistance();
    }

    inline int getWalkingArrivalTime(const Vertex vertex) noexcept {
        return sourceDepartureTime + initialTransfers.getForwardDistance(vertex);
    }

    inline int getWalkingTravelTime() const noexcept requires OneToOne {
        return initialTransfers.getDistance();
    }

    inline int getWalkingTravelTime(const Vertex vertex) noexcept {
        return initialTransfers.getForwardDistance(vertex);
    }

    inline std::vector<Vertex> getPath(const Vertex vertex) {
        Assert(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return journeyToPath(getJourneys(vertex).back());
    }

    inline std::vector<std::string> getRouteDescription(const Vertex vertex) {
        Assert(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return data.journeyToText(getJourneys(vertex).back());
    }

    inline int getArrivalTime(const Vertex vertex, const size_t numberOfTrips) noexcept {
        if (data.isStop(vertex)) {
            size_t round = numberOfTrips * RoundFactor;
            if constexpr (UseMinTransferTimes) {
                if ((round < roundIndex) && (roundLabel(round + 1, StopId(vertex)).arrivalTime < roundLabel(round, StopId(vertex)).arrivalTime)) round++;
            }
            Assert(roundLabel(round, StopId(vertex)).arrivalTime < never, "No label found for stop " << vertex << " in round " << round << "!");
            return roundLabel(round, StopId(vertex)).arrivalTime;
        } else {
            return dijkstraLabel(numberOfTrips, vertex).arrivalTime;
        }
    }

    inline const IndexedSet<false, Vertex>& getReachedVertices() const noexcept {
        return reachedVertices;
    }

    inline bool targetReached() const noexcept {
        return reachedVertices.contains(targetVertex);
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
        queue.clear();
        reachedVertices.clear();
        if constexpr (RESET) {
            if constexpr (TripPruning) resetLastTrip();
            timestamp = 1;
            std::vector<Round>().swap(rounds);
            std::vector<std::vector<DijkstraLabel>>(1).swap(labelByNumberOfTrips);
            std::vector<DijkstraLabel>(data.transferGraph.numVertices()).swap(labelByNumberOfTrips[0]);
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
    inline void resetLastTrip() noexcept requires TripPruning {
        std::vector<std::vector<std::vector<uint16_t>>>(1).swap(lastTrip);
        std::vector<std::vector<uint16_t>>(data.numberOfRoutes()).swap(lastTrip[0]);
        for (const RouteId route : data.routes()) {
            std::vector<uint16_t>(data.numberOfStopsInRoute(route) - 1, data.numberOfTripsInRoute(route) - 1).swap(lastTrip[0][route]);
        }
    }

    inline void initialize(const Vertex source, const int departureTime, const Vertex target) noexcept {
        sourceVertex = source;
        targetVertex = target;
        if constexpr (OneToOne) {
            if (data.isStop(target)) {
                targetStop = StopId(target);
            }
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
        if constexpr (UseMinTransferTimes) startNewRound();
    }

    inline void evaluateInitialTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            if constexpr (OneToOne) {
                if (stop == targetVertex) continue;
            }
            Assert(data.isStop(stop), "Reached POI " << stop << " is not a stop!");
            Assert(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            //The initial transfers are evaluated automatically when the label is updated
            roundLabel(0, StopId(stop));
            stopsUpdatedByTransfer.insert(StopId(stop));
            if constexpr (CollectReachedVertices) reachedVertices.insert(stop);
        }
        if constexpr (OneToOne) {
            if (initialTransfers.getDistance() != INFTY) {
                //The initial transfers are evaluated automatically when the label is updated
                dijkstraLabel(0, targetVertex);
                if constexpr (CollectReachedVertices) reachedVertices.insert(targetVertex);
            }
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            for (const RouteSegment& route : data.routesContainingStop(stop)) {
                Assert(data.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                Assert(data.stopIds[data.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == data.numberOfStopsInRoute(route.routeId)) continue;
                if constexpr (TripPruning) {
                    const uint16_t tripNum = lastTrip[roundIndex / RoundFactor][route.routeId][route.stopIndex];
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
                trip = data.tripOfRoute(route, lastTrip[roundIndex / RoundFactor][route][stopIndex]);
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
                if constexpr (TripPruning) lastTrip[roundIndex / RoundFactor][route][stopIndex] = uint16_t(((trip - firstTrip) / tripSize) - 1);
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
        if constexpr (OneToOne) {
            initialTransfers.run(sourceVertex, targetVertex);
        } else {
            initialTransfers.template run<FORWARD, BACKWARD>(sourceVertex);
        }
    }

    inline void relaxIntermediateTransfers() noexcept {
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        Assert(queue.empty(), "Queue still has " << queue.size() << " elements!");
        for (const StopId stop : stopsUpdatedByRoute) {
            const int arrivalTime = UseMinTransferTimes ? previousRoundLabel(stop).arrivalTime : currentRoundLabel(stop).arrivalTime;
            if constexpr (OneToOne) {
                if (initialTransfers.getBackwardDistance(stop) != INFTY) {
                    const int targetArrivalTime = arrivalTime + initialTransfers.getBackwardDistance(stop);
                    if (arrivalByEdge<false>(targetVertex, targetArrivalTime, stop)) {
                        arrivalByTransfer(targetStop, targetArrivalTime, stop, arrivalTime);
                    }
                }
            }
            if constexpr (UseMinTransferTimes) {
                for (Edge edge : data.transferGraph.edgesFrom(stop)) {
                    const Vertex to = data.transferGraph.get(ToVertex, edge);
                    if constexpr (OneToOne) {
                        if (to == targetVertex) continue;
                    }
                    profiler.countMetric(METRIC_EDGES);
                    arrivalByEdge<true>(to, arrivalTime + data.transferGraph.get(TravelTime, edge), stop);
                }
                for (Edge edge : minChangeTimeGraph.edgesFrom(stop)) {
                    const Vertex to = minChangeTimeGraph.get(ToVertex, edge);
                    if constexpr (OneToOne) {
                        if (to == targetVertex) continue;
                    }
                    profiler.countMetric(METRIC_EDGES);
                    arrivalByEdge<true>(to, arrivalTime + minChangeTimeGraph.get(TravelTime, edge), stop);
                }
                arrivalByEdge<true>(stop, arrivalTime + data.stopData[stop].minTransferTime, stop);
            } else {
                arrivalByEdge<true>(stop, arrivalTime, stop);
            }
        }
        dijkstra();
    }

    inline void dijkstra() noexcept {
        DijkstraLabel& targetLabel = dijkstraLabel(roundIndex / RoundFactor, OneToOne ? targetVertex : StopId(0));
        while (!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            if constexpr (OneToOne) {if (uLabel->arrivalTime > targetLabel.arrivalTime) break;}
            const Vertex u = Vertex(uLabel - &(labelByNumberOfTrips[roundIndex / RoundFactor][0]));
            for (Edge edge : data.transferGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = data.transferGraph.get(ToVertex, edge);
                if (v == targetVertex || v == uLabel->parent) continue;
                arrivalByEdge<true>(v, uLabel->arrivalTime + data.transferGraph.get(TravelTime, edge), uLabel->parent);
            }
            if constexpr (UseMinTransferTimes) {
                for (Edge edge : minChangeTimeGraph.edgesFrom(u)) {
                    profiler.countMetric(METRIC_EDGES);
                    const Vertex v = minChangeTimeGraph.get(ToVertex, edge);
                    if (v == targetVertex || v == uLabel->parent) continue;
                    arrivalByEdge<true>(v, uLabel->arrivalTime + minChangeTimeGraph.get(TravelTime, edge), uLabel->parent);
                }
            }
            if (data.isStop(u)) {
                arrivalByTransfer(StopId(u), uLabel->arrivalTime, uLabel->parent, getParentDepartureTime(uLabel->parent, UseMinTransferTimes ? roundIndex - 1 : roundIndex));
            }
            profiler.countMetric(METRIC_VERTICES);
        }
        if (!queue.empty()) {
            queue.clear();
        }
    }

    inline EarliestArrivalLabel& roundLabel(const size_t round, const StopId stop) noexcept {
        EarliestArrivalLabel& result = rounds[round][stop];
        if (result.timestamp != timestamp) {
            if (round >= RoundFactor) {
                result.arrivalTime = std::min(result.arrivalTime, roundLabel(round - RoundFactor, stop).arrivalTime);
            } else if (round == InitialTransferRound) {
                int distance = (stop == targetStop) ? initialTransfers.getDistance() : initialTransfers.getForwardDistance(stop);
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

    inline DijkstraLabel& dijkstraLabel(const size_t numTrips, const Vertex vertex) noexcept {
        DijkstraLabel& result = labelByNumberOfTrips[numTrips][vertex];
        if (result.timestamp != timestamp) {
            if (numTrips >= 1) {
                result.arrivalTime = std::min(result.arrivalTime, dijkstraLabel(numTrips - 1, vertex).arrivalTime);
            } else if (numTrips == 0) {
                int distance = (vertex == targetVertex) ? initialTransfers.getDistance() : initialTransfers.getForwardDistance(vertex);
                if (distance != INFTY) {
                    result.arrivalTime = sourceDepartureTime + distance;
                    result.parent = sourceVertex;
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
            if (rounds.size() < RoundFactor) {
                rounds.emplace_back(OneToOne ? data.numberOfStops() + 1 : data.numberOfStops());
            } else {
                rounds.emplace_back(rounds[rounds.size() - RoundFactor]);
            }
        }
        if (roundIndex / RoundFactor == labelByNumberOfTrips.size()) {
            labelByNumberOfTrips.emplace_back(labelByNumberOfTrips.back());
            if constexpr (TripPruning) lastTrip.emplace_back(lastTrip.back());
        }
    }

    inline bool arrivalByRoute(const StopId stop, const int arrivalTime) noexcept {
        Assert(roundIndex > 0, "arrivalByRoute cannot be used in the first round!");
        Assert(data.isStop(stop), "Stop " << stop << " is out of range!");
        Assert(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], stop: " << stop << ")!");
        if constexpr (OneToOne) {if ((currentRoundLabel(targetStop).arrivalTime <= arrivalTime) || (previousRoundLabel(targetStop).arrivalTime <= arrivalTime)) return false;}
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if ((label.arrivalTime <= arrivalTime) || (previousRoundLabel(stop).arrivalTime <= arrivalTime)) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        label.arrivalTime = arrivalTime;
        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    template<bool ADD_TO_QUEUE>
    inline bool arrivalByEdge(const Vertex vertex, const int arrivalTime, const Vertex parent) noexcept {
        Assert(data.isStop(parent), "Parent vertex (" << parent << ") is not a stop");
        Assert(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        DijkstraLabel& label = dijkstraLabel(roundIndex / RoundFactor, vertex);
        if (label.arrivalTime <= arrivalTime) return false;

        label.arrivalTime = arrivalTime;
        label.parent = parent;
        if constexpr (ADD_TO_QUEUE) queue.update(&label);
        if constexpr (CollectReachedVertices) reachedVertices.insert(vertex);
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
    }

    inline int getParentDepartureTime(const Vertex vertex, const size_t round) noexcept {
        Assert(data.isStop(vertex) || vertex == sourceVertex, "Vertex " << vertex << " should be a stop or the source vertex!");
        if (vertex == sourceVertex) {
            return sourceDepartureTime;
        } else {
            return roundLabel(round, StopId(vertex)).arrivalTime;
        }
    }

    inline void getJourney(std::vector<Journey>& journeys, size_t round, Vertex vertex) noexcept {
        if constexpr (UseMinTransferTimes) {
            Assert(round % 2 == 0, "Journeys can only be computed for even rounds! (" << round << ")!");
        }
        Journey journey;
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                if ((round < roundIndex) && (roundLabel(round + 1, StopId(vertex)).arrivalTime < roundLabel(round, StopId(vertex)).arrivalTime)) round++;
            }
            if (roundLabel(round, StopId(vertex)).arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        } else {
            const DijkstraLabel& label = dijkstraLabel(round / RoundFactor, vertex);
            if (label.arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
            Assert(label.parent != noVertex, "Vertex " << vertex << ": Label with " << round / RoundFactor << " trips has no parent! (arrival time: " << label.arrivalTime << ")");
            journey.emplace_back(label.parent, (vertex == targetStop) ? targetVertex : vertex, getParentDepartureTime(label.parent, round), label.arrivalTime, false, noRouteId);
            vertex = label.parent;
        }
        while (vertex != sourceVertex) {
            Assert(round != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = roundLabel(round, StopId(vertex));
            journey.emplace_back(label.parent, vertex, label.parentDepartureTime, label.arrivalTime, label.usesRoute, label.routeId);
            vertex = label.parent;
            if constexpr (UseMinTransferTimes) {
                round--;
            } else {
                if (label.usesRoute) round--;
            }
        }
        journeys.emplace_back(Vector::reverse(journey));
    }

    inline void getArrival(std::vector<ArrivalLabel>& labels, size_t round, const Vertex vertex) noexcept {
        if constexpr (UseMinTransferTimes) {
            Assert(round % 2 == 0, "Arrivals can only be computed for even rounds! (" << round << ")!");
        }
        const size_t numberOfTrips = round / RoundFactor;
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                if ((round < roundIndex) && (roundLabel(round + 1, StopId(vertex)).arrivalTime < roundLabel(round, StopId(vertex)).arrivalTime)) round++;
            }
            const EarliestArrivalLabel& label = roundLabel(round, StopId(vertex));
            if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
            labels.emplace_back(label.arrivalTime, numberOfTrips);
        } else {
            const DijkstraLabel& label = dijkstraLabel(numberOfTrips, vertex);
            if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
            labels.emplace_back(label.arrivalTime, numberOfTrips);
        }
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
    TransferGraph minChangeTimeGraph;

    InitialTransferType initialTransfers;

    std::vector<Round> rounds;
    size_t roundIndex;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    Vertex targetVertex;
    StopId targetStop; //One-to-one only
    int sourceDepartureTime;
    int timestamp;

    std::vector<std::vector<DijkstraLabel>> labelByNumberOfTrips;
    ExternalKHeap<2, DijkstraLabel> queue;

    std::vector<std::vector<std::vector<uint16_t>>> lastTrip;

    Profiler profiler;

    IndexedSet<false, Vertex> reachedVertices;

};

}
