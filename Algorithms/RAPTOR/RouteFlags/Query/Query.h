#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../../CH/CH.h"

#include "../../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/Shortcuts.h"
#include "../../../../DataStructures/Container/Set.h"
#include "../../../../DataStructures/Container/Map.h"
#include "../../../../DataStructures/Graph/Utils/Intersection.h"
#include "../../../../DataStructures/Graph/Utils/Union.h"

#include "../../Profiler.h"
#include "../../InitialTransfers.h"

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE, bool TARGET_PRUNING, typename PROFILER = NoProfiler, bool DIRECTION = FORWARD>
class Query {

public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    static constexpr bool TargetPruning = TARGET_PRUNING;
    using Profiler = PROFILER;
    static constexpr bool Direction = DIRECTION;
    static constexpr bool QueryDeparture = Direction;
    static constexpr bool QueryArrival = !Direction;
    using Type = Query<Flags, FinePartitionType, TargetPruning, Profiler, Direction>;
    using RouteFlagsDataType = RouteFlagsData<Flags, FinePartitionType>;

private:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parentDepartureTime(never), parent(noVertex), usesRoute(false), routeId(noRouteId) {}
        int arrivalTime;
        int parentDepartureTime;
        Vertex parent;
        bool usesRoute;
        union {
            RouteId routeId;
            Edge transferId;
        };
    };
    using Round = std::vector<EarliestArrivalLabel>;

    struct RouteInformation {
        RouteInformation() : stopIndex(noStopIndex), earliestArrival(never) {}

        RouteInformation(const StopIndex stopIndex, const int earliestArrival) :
            stopIndex(stopIndex),
            earliestArrival(earliestArrival) {}

        StopIndex stopIndex;
        int earliestArrival;

        inline void update(const StopIndex newStopIndex, const int newEarliestArrival) noexcept {
            stopIndex = std::min(stopIndex, newStopIndex);
            earliestArrival = std::min(earliestArrival, newEarliestArrival);
        }
    };

public:
    Query(const RouteFlagsDataType& flagsData, const CH::CH& chData, const Profiler& profilerTemplate = Profiler()) :
        raptorData(flagsData.getRaptorData(Direction)),
        flagsData(flagsData),
        forwardShortcuts(flagsData.getShortcuts(Direction)),
        backwardShortcuts(flagsData.getShortcuts(!Direction)),
        forwardCoarseShortcuts(NULL),
        forwardFineShortcuts(NULL),
        backwardFineShortcuts(NULL),
        intersectedShortcuts(forwardShortcuts.fineCellShortcuts[0][0], forwardShortcuts.fineCellShortcuts[0][0]),
        unifiedShortcuts(forwardShortcuts.fineCellShortcuts[0][0], forwardShortcuts.fineCellShortcuts[0][0]),
        initialTransfers(chData, Direction, raptorData.numberOfStops()),
        earliestArrival(raptorData.numberOfStops() + 1, never),
        stopsUpdatedByRoute(raptorData.numberOfStops() + 1),
        stopsUpdatedByTransfer(raptorData.numberOfStops() + 1),
        routesServingUpdatedStops(raptorData.numberOfRoutes()),
        sourceDepartureTime(intMax),
        sourceVertex(noVertex),
        mainSourceCoarseCell(-1),
        targetVertex(noVertex),
        targetStop(noStop),
        mainTargetCoarseCell(-1),
        profiler(profilerTemplate) {
        AssertMsg(raptorData.hasImplicitBufferTimes(), "Route-Flags requires implicit buffer times!");
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_VERTICES, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
        if constexpr (Direction == BACKWARD) {
            forwardShortcuts.revert();
            backwardShortcuts.revert();
        }
    }

    inline void run(const Vertex source, const int departureTime, const Vertex target, const size_t maxRounds = INFTY) noexcept {
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime, target);
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.startPhase();
        relaxInitialTransfers();
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.doneRound();

        for (size_t i = 0; i < maxRounds; i++) {
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
            if (stopsUpdatedByRoute.empty()) {
                profiler.doneRound();
                break;
            }
            profiler.startPhase();
            relaxIntermediateTransfers();
            profiler.donePhase(PHASE_TRANSFERS);
            profiler.doneRound();
        }
        profiler.doneRound();
        profiler.done();
    }

    inline std::vector<Journey> getJourneys() const noexcept {
        return getJourneys(targetStop);
    }

    inline std::vector<Journey> getJourneys(const Vertex vertex) const noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<Journey> journeys;
        for (size_t i = 0; i < rounds.size(); i++) {
            getJourney(journeys, i, target);
        }
        return journeys;
    }

    inline std::vector<ArrivalLabel> getArrivals() const noexcept {
        return getArrivals(targetStop);
    }

    inline std::vector<ArrivalLabel> getArrivals(const Vertex vertex) const noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<ArrivalLabel> labels;
        for (size_t i = 0; i < rounds.size(); i++) {
            getArrival(labels, i, target);
        }
        return labels;
    }

    inline std::vector<int> getArrivalTimes() const noexcept {
        return getArrivalTimes(targetStop);
    }

    inline std::vector<int> getArrivalTimes(const Vertex vertex) const noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        std::vector<int> arrivalTimes;
        for (size_t i = 0; i < rounds.size(); i++) {
            getArrivalTime(arrivalTimes, i, target);
        }
        return arrivalTimes;
    }

    inline bool reachable(const Vertex vertex) const noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return earliestArrival[target] < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) const noexcept {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return earliestArrival[target];
    }

    inline std::vector<Vertex> getPath(const Vertex vertex) const {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return journeyToPath(getJourneys(target).back());
    }

    inline std::vector<std::string> getRouteDescription(const Vertex vertex) const {
        const StopId target = (vertex == targetVertex) ? (targetStop) : (StopId(vertex));
        return raptorData.journeyToText(getJourneys(target).back());
    }

    inline const TransferGraph& getShortcuts() const noexcept {
        return forwardShortcuts.allShortcuts;
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    template<bool RESET_CAPACITIES = false>
    inline void clear() noexcept {
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        targetStop = StopId(raptorData.numberOfStops());
        if constexpr (RESET_CAPACITIES) {
            std::vector<Round>().swap(rounds);
            std::vector<int>(earliestArrival.size(), never).swap(earliestArrival);
        } else {
            rounds.clear();
            Vector::fill(earliestArrival, never);
        }
    }

    inline void reset() noexcept {
        clear<true>();
    }

    inline void initialize(const Vertex source, const int departureTime, const Vertex target) noexcept {
        sourceDepartureTime = departureTime;
        sourceVertex = source;
        sourceCoarseCells = flagsData.getCoarsePartition().cellsOfVertex(source);
        mainSourceCoarseCell = sourceCoarseCells[0];
        sourceFineCells = flagsData.getFinePartition().cellsOfVertex(source);
        targetVertex = target;
        if (raptorData.isStop(target)) {
            targetStop = StopId(target);
        }
        targetCoarseCells = flagsData.getCoarsePartition().cellsOfVertex(target);
        mainTargetCoarseCell = targetCoarseCells[0];
        targetFineCells = flagsData.getFinePartition().cellsOfVertex(target);

        forwardCoarseShortcuts = &(forwardShortcuts.coarseCellShortcuts[mainTargetCoarseCell]);
        forwardFineShortcuts = &(forwardShortcuts.fineCellShortcuts[mainTargetCoarseCell][sourceFineCells[0]]);
        backwardFineShortcuts = &(backwardShortcuts.fineCellShortcuts[mainSourceCoarseCell][targetFineCells[0]]);
        intersectedShortcuts.setGraphs(*forwardFineShortcuts, *backwardFineShortcuts);
        unifiedShortcuts.setGraphs(*forwardFineShortcuts, *backwardFineShortcuts);

        startNewRound();
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        std::vector<StopId> stopsToRemove;
        for (const StopId stop : stopsUpdatedByTransfer) {
            AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
            AssertMsg(flagsData.isNecessaryStop(Direction, QueryDeparture, stop, sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells), "Unnecessary stop " << stop << " was updated!");
            const int arrivalTime = previousRound()[stop].arrivalTime;
            AssertMsg(arrivalTime < never, "Updated stop has arrival time = never!");
            if constexpr (TargetPruning) {
                if (earliestArrival[targetStop] <= arrivalTime) {
                    stopsToRemove.push_back(stop);
                    continue;
                }
            }
            for (const RouteSegment& route : raptorData.routesContainingStop(stop)) {
                AssertMsg(raptorData.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                AssertMsg(raptorData.stopIds[raptorData.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (routesServingUpdatedStops.contains(route.routeId)) {
                    routesServingUpdatedStops[route.routeId].update(route.stopIndex, arrivalTime);
                } else {
                    routesServingUpdatedStops.insert(route.routeId, RouteInformation(route.stopIndex, arrivalTime));
                }
            }
        }

        for (const StopId stop : stopsToRemove) {
            stopsUpdatedByTransfer.remove(stop);
        }

        std::vector<RouteId> routesToRemove;
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            TripIterator tripIterator = raptorData.getTripIterator(route, routesServingUpdatedStops[route].stopIndex);
            routesServingUpdatedStops[route].stopIndex = noStopIndex;
            while (tripIterator.hasFurtherStops()) {
                const StopId stop = tripIterator.stop();
                if (!stopsUpdatedByTransfer.contains(stop)) {
                    tripIterator.nextStop();
                    continue;
                }
                const StopIndex stopIndex = tripIterator.getStopIndex();
                if (!flagsData.isNecessaryRouteSegment(Direction, QueryDeparture, route, stopIndex, sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) {
                    tripIterator.nextStop();
                    continue;
                }

                if (hasNecessaryTripAfter(route, stopIndex, previousRound()[stop].arrivalTime)) {
                    routesServingUpdatedStops[route].stopIndex = stopIndex;
                    break;
                }
                tripIterator.nextStop();
            }

            if (routesServingUpdatedStops[route].stopIndex == noStopIndex) {
                routesToRemove.push_back(route);
            }
        }
        for (const RouteId route : routesToRemove) {
            routesServingUpdatedStops.remove(route);
        }
    }

    inline void scanRoutes() noexcept {
        stopsUpdatedByRoute.clear();
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            scanRoute(route, routesServingUpdatedStops[route]);
        }
    }

    inline void scanRoute(const RouteId route, const RouteInformation& info) noexcept {
        if constexpr (TargetPruning) {if (earliestArrival[targetStop] <= info.earliestArrival) return;}
        profiler.countMetric(METRIC_ROUTES);
        StopIndex stopIndex = info.stopIndex;
        const size_t tripSize = raptorData.numberOfStopsInRoute(route);
        AssertMsg(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ")!");

        const StopId* stops = raptorData.stopArrayOfRoute(route);
        const StopEvent* trip = raptorData.lastTripOfRoute(route);
        StopId stop = stops[stopIndex];
        AssertMsg(trip[stopIndex].departureTime >= previousRound()[stop].arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousRound()[stop].arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ")!");

        StopIndex parentIndex = stopIndex;
        const StopEvent* firstTrip = raptorData.firstTripOfRoute(route);
        bool skip = false;

        while (stopIndex < tripSize - 1) {
            const StopEvent* newTrip = trip;
            const StopEvent* tripCounter = trip;
            if (stopsUpdatedByTransfer.contains(stop) && flagsData.isNecessaryRouteSegment(Direction, QueryDeparture, route, stopIndex, sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) {
                while ((tripCounter > firstTrip) && ((tripCounter - tripSize + stopIndex)->departureTime >= previousRound()[stop].arrivalTime)) {
                    tripCounter -= tripSize;
                    if (flagsData.isNecessaryStopEvent(Direction, QueryDeparture, tripCounter + stopIndex - &(raptorData.stopEvents[0]), sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) newTrip = tripCounter;
                }
                if (newTrip != trip) {
                    trip = newTrip;
                    if constexpr (TargetPruning) {skip = earliestArrival[targetStop] <= trip[stopIndex].departureTime;}
                    parentIndex = stopIndex;
                }
            }
            stopIndex++;
            stop = stops[stopIndex];
            if (skip) continue;
            if (!flagsData.isNecessaryStopEvent(Direction, QueryArrival, trip + stopIndex - &(raptorData.stopEvents[0]), sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) continue;
            profiler.countMetric(METRIC_ROUTE_SEGMENTS);
            if (arrivalByRoute(stop, trip[stopIndex].arrivalTime)) {
                EarliestArrivalLabel& label = currentRound()[stop];
                label.parent = stops[parentIndex];
                label.parentDepartureTime = trip[parentIndex].departureTime;
                label.usesRoute = true;
                label.routeId = route;
            }
        }
    }

    template<typename GRAPH>
    inline void relaxEdges(const StopId stop, const GRAPH* graph) noexcept {
        for (const Edge edge : graph->edgesFrom(stop)) {
            const StopId toStop = StopId(graph->get(ToVertex, edge));
            if (toStop == targetStop) continue;
            const int arrivalTime = currentRound()[stop].arrivalTime + graph->get(TravelTime, edge);
            profiler.countMetric(METRIC_EDGES);
            if (arrivalByTransfer(toStop, arrivalTime)) {
                EarliestArrivalLabel& label = currentRound()[toStop];
                label.parent = stop;
                label.parentDepartureTime = currentRound()[stop].arrivalTime;
                label.usesRoute = false;
                label.transferId = edge;
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        initialTransfers.run(sourceVertex, targetVertex);
        for (const Vertex stop : initialTransfers.getForwardPOIs()) {
            if (stop == targetStop) continue;
            AssertMsg(raptorData.isStop(stop), "Reached POI " << stop << " is not a stop!");
            AssertMsg(initialTransfers.getForwardDistance(stop) != INFTY, "Vertex " << stop << " was not reached!");
            const int arrivalTime = sourceDepartureTime + initialTransfers.getForwardDistance(stop);
            if (arrivalByTransfer(StopId(stop), arrivalTime)) {
                EarliestArrivalLabel& label = currentRound()[stop];
                label.parent = sourceVertex;
                label.parentDepartureTime = sourceDepartureTime;
                label.usesRoute = false;
                label.transferId = noEdge;
            }
        }
        if (initialTransfers.getDistance() != INFTY) {
            const int arrivalTime = sourceDepartureTime + initialTransfers.getDistance();
            if (arrivalByTransfer(targetStop, arrivalTime)) {
                EarliestArrivalLabel& label = currentRound()[targetStop];
                label.parent = sourceVertex;
                label.parentDepartureTime = sourceDepartureTime;
                label.usesRoute = false;
                label.transferId = noEdge;
            }
        }
    }

    inline void relaxIntermediateTransfers() noexcept {
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        for (const StopId stop : stopsUpdatedByRoute) {
            if (flagsData.isNecessaryStop(Direction, QueryDeparture, stop, sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) {
                stopsUpdatedByTransfer.insert(stop);
            }
            if (mainSourceCoarseCell == mainTargetCoarseCell && flagsData.getCoarsePartition().isInCell(stop, mainSourceCoarseCell)) {
                relaxEdges(stop, forwardCoarseShortcuts);
                relaxEdges(stop, &intersectedShortcuts);
            } else {
                const bool inSource = flagsData.getCoarsePartition().isInCell(stop, mainSourceCoarseCell);
                const bool inTarget = flagsData.getCoarsePartition().isInCell(stop, mainTargetCoarseCell);
                if (inSource && inTarget) {
                    relaxEdges(stop, &unifiedShortcuts);
                } else if (inSource) {
                    relaxEdges(stop, forwardFineShortcuts);
                } else if (inTarget) {
                    relaxEdges(stop, backwardFineShortcuts);
                } else {
                    relaxEdges(stop, &intersectedShortcuts);
                }
            }

            if (initialTransfers.getBackwardDistance(stop) != INFTY) {
                const int arrivalTime = currentRound()[stop].arrivalTime + initialTransfers.getBackwardDistance(stop);
                if (arrivalByTransfer(targetStop, arrivalTime)) {
                    EarliestArrivalLabel& label = currentRound()[targetStop];
                    label.parent = stop;
                    label.parentDepartureTime = currentRound()[stop].arrivalTime;
                    label.usesRoute = false;
                    label.transferId = noEdge;
                }
            }
        }
    }

    inline Round& currentRound() noexcept {
        AssertMsg(!rounds.empty(), "Cannot return current round, because no round exists!");
        return rounds.back();
    }

    inline Round& previousRound() noexcept {
        AssertMsg(rounds.size() >= 2, "Cannot return previous round, because less than two rounds exist!");
        return rounds[rounds.size() - 2];
    }

    inline void startNewRound() noexcept {
        rounds.emplace_back(raptorData.numberOfStops() + 1);
    }

    inline bool arrivalByRoute(const StopId stop, const int time) noexcept {
        AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
        if constexpr (TargetPruning) if (earliestArrival[targetStop] <= time) return false;
        if (earliestArrival[stop] <= time) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        currentRound()[stop].arrivalTime = time;
        earliestArrival[stop] = time;
        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    inline bool arrivalByTransfer(const StopId stop, const int time) noexcept {
        AssertMsg(raptorData.isStop(stop) || stop == targetStop, "Stop " << stop << " is out of range!");
        if constexpr (TargetPruning) if (earliestArrival[targetStop] <= time) return false;
        if (earliestArrival[stop] <= time) return false;
        profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
        currentRound()[stop].arrivalTime = time;
        earliestArrival[stop] = time;
        if (raptorData.isStop(stop) && flagsData.isNecessaryStop(Direction, QueryDeparture, stop, sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) {
            stopsUpdatedByTransfer.insert(stop);
        }
        return true;
    }

    inline void getJourney(std::vector<Journey>& journeys, size_t round, StopId stop) const noexcept {
        if (rounds[round][stop].arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        Journey journey;
        do {
            AssertMsg(round != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = rounds[round][stop];
            journey.emplace_back(label.parent, stop, label.parentDepartureTime, label.arrivalTime, label.usesRoute, label.routeId);
            AssertMsg(raptorData.isStop(label.parent) || label.parent == sourceVertex, "Backtracking parent pointers reached a vertex (" << label.parent << ")!");
            stop = StopId(label.parent);
            if (label.usesRoute) round--;
        } while (journey.back().from != sourceVertex);
        journeys.emplace_back(Vector::reverse(journey));
    }

    inline void getArrival(std::vector<ArrivalLabel>& labels, size_t round, const StopId stop) const noexcept {
        if (rounds[round][stop].arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
        labels.emplace_back(rounds[round][stop].arrivalTime, round);
    }

    inline void getArrivalTime(std::vector<int>& labels, size_t round, const StopId stop) const noexcept {
        labels.emplace_back(std::min(rounds[round][stop].arrivalTime, (labels.empty()) ? (never) : (labels.back())));
    }

    inline bool hasNecessaryTripAfter(const RouteId route, const StopIndex stopIndex, const int arrivalTime) const noexcept {
        TripIterator tripIterator = raptorData.getTripIterator(route, stopIndex);
        do {
            if (tripIterator.departureTime() < arrivalTime) return false;
            if (flagsData.isNecessaryStopEvent(Direction, QueryDeparture, tripIterator.stopEvent() - &(raptorData.stopEvents[0]), sourceCoarseCells, targetCoarseCells, sourceFineCells, targetFineCells)) return true;

        } while (tripIterator.decreaseTrip());
        return false;
    }

private:
    const Data& raptorData;
    const RouteFlagsDataType& flagsData;
    Shortcuts forwardShortcuts;
    Shortcuts backwardShortcuts;
    const TransferGraph* forwardCoarseShortcuts;
    const TransferGraph* forwardFineShortcuts;
    const TransferGraph* backwardFineShortcuts;
    Graph::Intersection<TransferGraph, TransferGraph> intersectedShortcuts;
    Graph::Union<TransferGraph, TransferGraph> unifiedShortcuts;

    BucketCHInitialTransfers initialTransfers;

    std::vector<Round> rounds;

    std::vector<int> earliestArrival;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<RouteInformation, false, RouteId> routesServingUpdatedStops;

    int sourceDepartureTime;

    Vertex sourceVertex;
    std::vector<size_t> sourceCoarseCells;
    size_t mainSourceCoarseCell;
    std::vector<size_t> sourceFineCells;
    Vertex targetVertex;
    StopId targetStop;
    std::vector<size_t> targetCoarseCells;
    size_t mainTargetCoarseCell;
    std::vector<size_t> targetFineCells;

    Profiler profiler;

};

}
