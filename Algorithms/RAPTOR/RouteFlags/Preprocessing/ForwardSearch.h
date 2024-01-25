#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>

#include "../../../../Helpers/Vector/Vector.h"

#include "../../InitialTransfers.h"

#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/PreprocessingJourney.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../DataStructures/Intermediate/Data.h"
#include "../../../../DataStructures/Container/Set.h"
#include "../../../../DataStructures/Container/Map.h"
#include "../../../../DataStructures/Container/ExternalKHeap.h"

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE, typename PROFILER, bool DIRECTION = FORWARD>
class ForwardSearch {

public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    using Profiler = PROFILER;
    constexpr static bool Direction = DIRECTION;
    using Type = ForwardSearch<Flags, FinePartitionType, Profiler, Direction>;
    using RouteFlagsDataType = RouteFlagsData<Flags, FinePartitionType>;

public:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parent(noVertex), departureEvent(NULL), arrivalEvent(NULL), timestamp(0), addedTimestamp(-1) {}
        int arrivalTime;
        Vertex parent;
        const StopEvent* departureEvent;
        const StopEvent* arrivalEvent;
        int timestamp;
        int addedTimestamp;

        inline bool usesRoute() const noexcept {
            return departureEvent != NULL;
        }
    };
    using Round = std::vector<EarliestArrivalLabel>;

    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel() : arrivalTime(never), parent(noVertex), timestamp(0), addedTimestamp(-1) {}
        int arrivalTime;
        Vertex parent;
        int timestamp;
        int addedTimestamp;

        inline bool hasSmallerKey(const DijkstraLabel* const other) const noexcept {
            return arrivalTime < other->arrivalTime;
        }
    };

public:
    ForwardSearch(const RouteFlagsDataType& flagsData, const Profiler& profilerTemplate = Profiler()) :
        raptorData(flagsData.getRaptorData(Direction)),
        flagsData(flagsData),
        shortcutGraph(flagsData.getShortcutGraph(Direction)),
        initialTransfers(flagsData.getRaptorData(Direction).transferGraph, flagsData.getRaptorData(!Direction).transferGraph, raptorData.numberOfStops(), TravelTime),
        roundIndex(-1),
        stopsUpdatedByRoute(raptorData.numberOfStops()),
        stopsUpdatedByTransfer(raptorData.numberOfStops()),
        routesServingUpdatedStops(raptorData.numberOfRoutes()),
        sourceVertex(noVertex),
        sourceDepartureTime(intMax),
        forceInitialTrip(false),
        initialTrip(NULL),
        timestamp(0),
        labelByNumberOfTrips(1, std::vector<DijkstraLabel>(raptorData.transferGraph.numVertices())),
        tempLabelByNumberOfTrips(1, std::vector<DijkstraLabel>(raptorData.transferGraph.numVertices())),
        profiler(profilerTemplate),
        reachedTargetVertices(flagsData.getFinePartition().highestBoundaryVertexId() + 1) {
        AssertMsg(raptorData.hasImplicitBufferTimes(), "Route-Flags requires implicit buffer times!");
        resetLastTrip();
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_VERTICES, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    template<bool CLEAR = true>
    inline void run(const Vertex source, const int departureTime, const bool force = false, const RouteSegment routeSegment = RouteSegment(), const StopEvent* trip = NULL, const size_t maxRounds = INFTY) noexcept {
        runInitialize<CLEAR>(source, departureTime, force, routeSegment, trip);
        runInitialTransfers();
        evaluateInitialTransfers();
        runRounds(maxRounds);
    }

    template<bool CLEAR = true>
    inline void runInitialize(const Vertex source, const int departureTime, const bool force = false, const RouteSegment routeSegment = RouteSegment(), const StopEvent* trip = NULL) noexcept {
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear<CLEAR>();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime, force, routeSegment, trip);
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

    inline void runAddSource(const StopId source, const int) noexcept {
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
        finalize();
        profiler.doneRound();
        profiler.done();
    }

    template<bool NEW_JOURNEYS_ONLY>
    inline std::vector<PreprocessingJourney> getJourneys(const Vertex vertex) noexcept {
        std::vector<PreprocessingJourney> journeys;
        int bestArrivalTime = never;
        for (size_t i = 0; i <= roundIndex; i++) {
            getJourney<NEW_JOURNEYS_ONLY>(journeys, i, vertex, bestArrivalTime);
        }
        return journeys;
    }

    template<bool NEW_JOURNEYS_ONLY>
    inline bool getJourney(const Vertex vertex, const size_t round, PreprocessingJourney& journey) noexcept {
        std::vector<PreprocessingJourney> journeys;
        int bestArrivalTime = never;
        getJourney<NEW_JOURNEYS_ONLY>(journeys, round, vertex, bestArrivalTime);
        if (journeys.empty()) return false;
        journey = journeys[0];
        return true;
    }

    template<bool NEW_JOURNEYS_ONLY>
    inline std::vector<ArrivalLabel> getArrivals(const Vertex vertex) noexcept {
        std::vector<ArrivalLabel> labels;
        for (size_t i = 0; i <= roundIndex; i++) {
            getArrival<NEW_JOURNEYS_ONLY>(labels, i, vertex);
        }
        return labels;
    }

    inline bool reachable(const Vertex vertex) noexcept {
        AssertMsg(raptorData.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return getEarliestArrivalTime(vertex) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) noexcept {
        AssertMsg(raptorData.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return dijkstraLabel(labelByNumberOfTrips.size() - 1, vertex).arrivalTime;
    }

    inline int getWalkingArrivalTime(const Vertex vertex) noexcept {
        return sourceDepartureTime + initialTransfers.getForwardDistance(vertex);
    }

    inline int getWalkingTravelTime(const Vertex vertex) noexcept {
        return initialTransfers.getForwardDistance(vertex);
    }

    inline int getArrivalTime(const Vertex vertex, const size_t numberOfTrips) noexcept {
        return dijkstraLabel(numberOfTrips, vertex).arrivalTime;
    }

    inline const IndexedSet<false, Vertex>& getReachedTargetVertices() const noexcept {
        return reachedTargetVertices;
    }

    inline const std::vector<Vertex>& getStopsReachedByWalking() const noexcept {
        return initialTransfers.getForwardPOIs();
    }

    template<bool RESET = true>
    inline void clear() noexcept {
        roundIndex = -1;
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        sourceVertex = noVertex;
        queue.clear();
        reachedTargetVertices.clear();
        if constexpr (RESET) {
            resetLastTrip();
            timestamp = 1;
            std::vector<Round>().swap(rounds);
            std::vector<Round>().swap(tempRounds);
            std::vector<std::vector<DijkstraLabel>>(1).swap(labelByNumberOfTrips);
            std::vector<DijkstraLabel>(raptorData.transferGraph.numVertices()).swap(labelByNumberOfTrips[0]);
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

    inline int getSourceDepartureTime() const noexcept {
        return sourceDepartureTime;
    }

    inline size_t maxNumberOfTrips() const noexcept {
        return roundIndex;
    }

private:
    inline void resetLastTrip() noexcept {
        std::vector<std::vector<std::vector<uint16_t>>>(1).swap(lastTrip);
        std::vector<std::vector<uint16_t>>(raptorData.numberOfRoutes()).swap(lastTrip[0]);
        for (const RouteId route : raptorData.routes()) {
            std::vector<uint16_t>(raptorData.numberOfStopsInRoute(route) - 1, raptorData.numberOfTripsInRoute(route) - 1).swap(lastTrip[0][route]);
        }
    }

    inline void initialize(const Vertex source, const int departureTime, const bool force, const RouteSegment routeSegment, const StopEvent* trip) noexcept {
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        forceInitialTrip = force;
        initialRouteSegment = routeSegment;
        initialTrip = trip;
        startNewRound();
        if (raptorData.isStop(source)) {
            stopsUpdatedByRoute.insert(StopId(source));
            EarliestArrivalLabel& label = currentRoundLabel(StopId(source));
            label.arrivalTime = departureTime;
            label.parent = source;
            label.addedTimestamp = timestamp;
        }
    }

    inline void evaluateInitialTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            AssertMsg(raptorData.isStop(stop), "Reached POI " << stop << " is not a stop!");
            AssertMsg(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            //The initial transfers are evaluated automatically when the label is updated
            roundLabel(0, StopId(stop));
            stopsUpdatedByTransfer.insert(StopId(stop));
            reachedTargetVertices.insert(stop);
        }
    }

    inline void finalize() noexcept {
        if (raptorData.isStop(sourceVertex)) {
            roundLabel(0, StopId(sourceVertex)).arrivalTime = never;
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        if (forceInitialTrip && roundIndex == 1) {
            routesServingUpdatedStops.insert(initialRouteSegment.routeId, initialRouteSegment.stopIndex);
        } else {
            const size_t numberOfTrips = forceInitialTrip ? roundIndex - 1 : roundIndex;
            for (const StopId stop : stopsUpdatedByTransfer) {
                for (const RouteSegment& route : raptorData.routesContainingStop(stop)) {
                    AssertMsg(raptorData.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                    AssertMsg(raptorData.stopIds[raptorData.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                    if (route.stopIndex + 1 == raptorData.numberOfStopsInRoute(route.routeId)) continue;
                    const uint16_t tripNum = lastTrip[numberOfTrips][route.routeId][route.stopIndex];
                    if (tripNum == uint16_t(-1)) continue;
                    const StopEvent* trip = raptorData.tripOfRoute(route.routeId, tripNum) + route.stopIndex;
                    if (trip->departureTime < previousRoundLabel(stop).arrivalTime) continue;
                    if (routesServingUpdatedStops.contains(route.routeId)) {
                        routesServingUpdatedStops[route.routeId] = std::min(routesServingUpdatedStops[route.routeId], route.stopIndex);
                    } else {
                        routesServingUpdatedStops.insert(route.routeId, route.stopIndex);
                    }
                }
            }
        }
    }

    inline void scanRoutes() noexcept {
        stopsUpdatedByRoute.clear();
        const bool forcedTrip = forceInitialTrip && roundIndex == 1;
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            profiler.countMetric(METRIC_ROUTES);
            StopIndex stopIndex = routesServingUpdatedStops[route];
            const size_t tripSize = raptorData.numberOfStopsInRoute(route);
            AssertMsg(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ", RoundIndex: " << roundIndex << ")!");

            const StopId* stops = raptorData.stopArrayOfRoute(route);
            const StopEvent* firstTrip = raptorData.firstTripOfRoute(route);
            const StopEvent* trip;
            if (forcedTrip) {
                trip = initialTrip;
            } else {
                const size_t numberOfTrips = forceInitialTrip ? roundIndex - 1 : roundIndex;
                trip = raptorData.tripOfRoute(route, lastTrip[numberOfTrips][route][stopIndex]);
            }
            if (trip < firstTrip) continue;
            StopId stop = stops[stopIndex];
            AssertMsg(forcedTrip || trip[stopIndex].departureTime >= previousRoundLabel(stop).arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousRoundLabel(stop).arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ", RoundIndex: " << roundIndex << ")!");

            StopIndex parentIndex = stopIndex;
            while (stopIndex < tripSize - 1) {
                if ((!forcedTrip || stopIndex != initialRouteSegment.stopIndex) && stopsUpdatedByTransfer.contains(stop)) {
                    while ((trip > firstTrip) && ((trip - tripSize + stopIndex)->departureTime >= previousRoundLabel(stop).arrivalTime)) {
                        trip -= tripSize;
                        parentIndex = stopIndex;
                    }
                }
                if (!forceInitialTrip) lastTrip[roundIndex][route][stopIndex] = uint16_t(((trip - firstTrip) / tripSize) - 1);
                stopIndex++;
                stop = stops[stopIndex];
                profiler.countMetric(METRIC_ROUTE_SEGMENTS);
                if (arrivalByRoute(stop, trip[stopIndex].arrivalTime)) {
                    EarliestArrivalLabel& label = currentRoundLabel(stop);
                    label.parent = stops[parentIndex];
                    label.departureEvent = &trip[parentIndex];
                    label.arrivalEvent = &trip[stopIndex];
                }
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        initialTransfers.template run<FORWARD, BACKWARD>(sourceVertex);
    }

    inline void relaxIntermediateTransfers() noexcept {
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        AssertMsg(queue.empty(), "Queue still has " << queue.size() << " elements!");
        for (const StopId stop : stopsUpdatedByRoute) {
            stopsUpdatedByTransfer.insert(stop);
            arrivalByEdge<true>(stop, currentRoundLabel(stop).arrivalTime, stop);
        }
        dijkstra();
        for (const StopId fromStop : stopsUpdatedByRoute) {
            const int departureTime = currentRoundLabel(fromStop).arrivalTime;
            for (const Edge edge : shortcutGraph.edgesFrom(fromStop)) {
                const int arrivalTime = departureTime + shortcutGraph.get(TravelTime, edge);
                const StopId toStop = StopId(shortcutGraph.get(ToVertex, edge));
                profiler.countMetric(METRIC_EDGES);
                arrivalByShortcut(toStop, arrivalTime, fromStop);
            }
        }
    }

    inline void dijkstra() noexcept {
        while (!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(dijkstraLabel(roundIndex, Vertex(0))));
            for (Edge edge : raptorData.transferGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = raptorData.transferGraph.get(ToVertex, edge);
                if (v == uLabel->parent) continue;
                arrivalByEdge<true>(v, uLabel->arrivalTime + raptorData.transferGraph.get(TravelTime, edge), uLabel->parent);
            }
            profiler.countMetric(METRIC_VERTICES);
        }
        if (!queue.empty()) {
            queue.clear();
        }
    }

    inline EarliestArrivalLabel& roundLabel(const size_t round, const StopId stop) noexcept {
        return forceInitialTrip ? tempRoundLabel(round, stop) : realRoundLabel(round, stop);
    }

    inline EarliestArrivalLabel& tempRoundLabel(const size_t round, const StopId stop) noexcept {
        EarliestArrivalLabel& result = tempRounds[round][stop];
        if (result.timestamp != timestamp) {
            if (round == 0) {
                result = realRoundLabel(0, stop);
            } else {
                //TODO: Is this too expensive?
                EarliestArrivalLabel& realLabel = realRoundLabel(round - 1, stop);
                EarliestArrivalLabel& tempLabel = tempRoundLabel(round - 1, stop);
                if (realLabel.arrivalTime < tempLabel.arrivalTime) {
                    result = realLabel;
                } else {
                    result = tempLabel;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline EarliestArrivalLabel& realRoundLabel(const size_t round, const StopId stop) noexcept {
        EarliestArrivalLabel& result = rounds[round][stop];
        if (result.timestamp != timestamp) {
            if (round > 0) {
                result.arrivalTime = std::min(result.arrivalTime, realRoundLabel(round - 1, stop).arrivalTime);
            } else {
                int distance = initialTransfers.getForwardDistance(stop);
                if (distance != INFTY) {
                    profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
                    result.arrivalTime = sourceDepartureTime + distance;
                    result.parent = sourceVertex;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline DijkstraLabel& dijkstraLabel(const size_t numTrips, const Vertex vertex) noexcept {
        return forceInitialTrip ? tempDijkstraLabel(numTrips, vertex) : realDijkstraLabel(numTrips, vertex);
    }

    inline DijkstraLabel& tempDijkstraLabel(const size_t numTrips, const Vertex vertex) noexcept {
        DijkstraLabel& result = tempLabelByNumberOfTrips[numTrips][vertex];
        if (result.timestamp != timestamp) {
            if (numTrips == 0) {
                result = realDijkstraLabel(0, vertex);
            } else {
                //TODO: Is this too expensive?
                DijkstraLabel& realLabel = realDijkstraLabel(numTrips - 1, vertex);
                DijkstraLabel& tempLabel = tempDijkstraLabel(numTrips - 1, vertex);
                if (realLabel.arrivalTime < tempLabel.arrivalTime) {
                    result = realLabel;
                } else {
                    result = tempLabel;
                }
            }
            result.timestamp = timestamp;
        }
        return result;
    }

    inline DijkstraLabel& realDijkstraLabel(const size_t numTrips, const Vertex vertex) noexcept {
        DijkstraLabel& result = labelByNumberOfTrips[numTrips][vertex];
        if (result.timestamp != timestamp) {
            if (numTrips >= 1) {
                result.arrivalTime = std::min(result.arrivalTime, realDijkstraLabel(numTrips - 1, vertex).arrivalTime);
            } else if (numTrips == 0) {
                int distance = initialTransfers.getForwardDistance(vertex);
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
        AssertMsg(roundIndex < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        return roundLabel(roundIndex, stop);
    }

    inline EarliestArrivalLabel& previousRoundLabel(const StopId stop) noexcept {
        AssertMsg(roundIndex - 1 < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        AssertMsg(roundIndex > 0, "Cannot return previous round, because no round exists!");
        return roundLabel(roundIndex - 1, stop);
    }

    inline void startNewRound() noexcept {
        AssertMsg(roundIndex + 1 <= rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        roundIndex++;
        if (roundIndex == rounds.size()) {
            if (rounds.empty()) {
                rounds.emplace_back(raptorData.numberOfStops());
                tempRounds.emplace_back(raptorData.numberOfStops());
            } else {
                rounds.emplace_back(rounds.back());
                tempRounds.emplace_back(tempRounds.back());
            }
        }
        if (roundIndex == labelByNumberOfTrips.size()) {
            labelByNumberOfTrips.emplace_back(labelByNumberOfTrips.back());
            tempLabelByNumberOfTrips.emplace_back(tempLabelByNumberOfTrips.back());
            lastTrip.emplace_back(lastTrip.back());
        }
    }

    inline bool dominates(const int arrivalTimeA, const int arrivalTimeB, const int addedTimestamp) const noexcept {
        if constexpr (Direction == BACKWARD) {
            if (forceInitialTrip) {
                return (arrivalTimeA < arrivalTimeB || (arrivalTimeA == arrivalTimeB && addedTimestamp == timestamp));
            }
        } else {
            suppressUnusedParameterWarning(addedTimestamp);
        }
        return (arrivalTimeA <= arrivalTimeB);
    }

    inline bool arrivalByRoute(const StopId stop, const int arrivalTime) noexcept {
        AssertMsg(roundIndex > 0, "arrivalByRoute cannot be used in the first round!");
        AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], stop: " << stop << ")!");
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if (dominates(label.arrivalTime, arrivalTime, label.addedTimestamp)) return false;
        if (previousRoundLabel(stop).arrivalTime <= arrivalTime) return false;

        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        label.arrivalTime = arrivalTime;
        label.addedTimestamp = timestamp;
        stopsUpdatedByRoute.insert(stop);
        reachedTargetVertices.insert(stop);
        return true;
    }

    template<bool ADD_TO_QUEUE>
    inline bool arrivalByEdge(const Vertex vertex, const int arrivalTime, const Vertex parent) noexcept {
        AssertMsg(raptorData.isStop(parent), "Parent vertex (" << parent << ") is not a stop");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        DijkstraLabel& label = dijkstraLabel(roundIndex, vertex);
        if (dominates(label.arrivalTime, arrivalTime, label.addedTimestamp)) return false;

        label.arrivalTime = arrivalTime;
        label.parent = parent;
        label.addedTimestamp = timestamp;
        if constexpr (ADD_TO_QUEUE) queue.update(&label);
        if (flagsData.isStopOrFineBoundaryVertex(vertex)) {
            AssertMsg(vertex < reachedTargetVertices.capacity(), "Vertex " << vertex << " is not within capacity of reachedTargetVertices: " << reachedTargetVertices.capacity());
            reachedTargetVertices.insert(vertex);
        }
        return true;
    }

    inline void arrivalByShortcut(const StopId stop, const int arrivalTime, const StopId parent) noexcept {
        AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");
        EarliestArrivalLabel& label = currentRoundLabel(stop);
        if (dominates(label.arrivalTime, arrivalTime, label.addedTimestamp)) return;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        label.arrivalTime = arrivalTime;
        label.parent = parent;
        label.departureEvent = NULL;
        label.arrivalEvent = NULL;
        label.addedTimestamp = timestamp;
        stopsUpdatedByTransfer.insert(stop);
    }

    template<bool NEW_JOURNEYS_ONLY>
    inline void getJourney(std::vector<PreprocessingJourney>& journeys, size_t round, Vertex vertex, int& bestArrivalTime) noexcept {
        PreprocessingJourney journey;
        const DijkstraLabel& label = dijkstraLabel(round, vertex);
        if (label.arrivalTime >= bestArrivalTime) return;
        bestArrivalTime = label.arrivalTime;
        if constexpr (NEW_JOURNEYS_ONLY) {
            if (label.addedTimestamp != timestamp) return;
        }
        AssertMsg(label.parent != noVertex, "Vertex " << vertex << ": Label with " << round << " trips has no parent! (arrival time: " << label.arrivalTime << ")");
        vertex = label.parent;

        while (vertex != sourceVertex) {
            AssertMsg(raptorData.isStop(vertex), "Vertex is not a stop!");
            AssertMsg(round != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = roundLabel(round, StopId(vertex));
            if (label.usesRoute()) {
                AssertMsg(label.arrivalEvent, "Invalid arrival event for stop " << vertex << " in round " << round);
                journey.departureStopEvents.emplace_back(label.departureEvent);
                journey.arrivalStopEvents.emplace_back(label.arrivalEvent);
                round--;
            } else if (label.parent != vertex) {
                const int parentArrivalTime = dijkstraLabel(round, label.parent).arrivalTime;
                journey.shortcuts.emplace_back(Shortcut{label.parent, vertex, label.arrivalTime - parentArrivalTime});
            }
            vertex = label.parent;
        }
        if (!journey.empty()) {
            if (forceInitialTrip && !journey.departureStopEvents.empty()) journey.departureStopEvents.pop_back();
            journeys.emplace_back(journey);
        }
    }

    template<bool NEW_JOURNEYS_ONLY>
    inline void getArrival(std::vector<ArrivalLabel>& labels, size_t round, const Vertex vertex) noexcept {
        const DijkstraLabel& label = dijkstraLabel(round, vertex);
        if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
        if constexpr (NEW_JOURNEYS_ONLY) {
            if (label.addedTimestamp != timestamp) return;
        }
        labels.emplace_back(label.arrivalTime, round);
    }

    inline void printRoundsForStop(const StopId stop) noexcept {
        AssertMsg(raptorData.isStop(stop), stop << " is not a valid stop!");
        std::cout << "Raptor Label for stop " << stop << ":" << std::endl;
        std::cout << std::setw(10) << "Round" << std::setw(14) << "arrivalTime" << std::setw(14) << "parent" << std::endl;
        for (size_t i = 0; i <= roundIndex; i++) {
            EarliestArrivalLabel& label = roundLabel(i, stop);
            std::cout << std::setw(10) << i << std::setw(14) << label.arrivalTime << std::setw(14) << label.parent << std::endl;
        }
    }

private:
    const Data& raptorData;
    const RouteFlagsDataType& flagsData;
    const TransferGraph& shortcutGraph;

    DijkstraInitialTransfers initialTransfers;

    std::vector<Round> rounds;
    std::vector<Round> tempRounds;
    size_t roundIndex;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    int sourceDepartureTime;
    bool forceInitialTrip;
    RouteSegment initialRouteSegment;
    const StopEvent* initialTrip;
    int timestamp;

    std::vector<std::vector<DijkstraLabel>> labelByNumberOfTrips;
    std::vector<std::vector<DijkstraLabel>> tempLabelByNumberOfTrips;
    ExternalKHeap<2, DijkstraLabel> queue;

    std::vector<std::vector<std::vector<uint16_t>>> lastTrip;

    Profiler profiler;

    IndexedSet<false, Vertex> reachedTargetVertices;

};

}
