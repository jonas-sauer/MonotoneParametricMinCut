#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Vector/Vector.h"
#include "Profiler.h"

#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/Entities/Journey.h"
#include "../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Container/Map.h"
#include "../../DataStructures/Container/ExternalKHeap.h"

namespace RAPTOR {

template<bool TARGET_PRUNING, typename PROFILER, bool COLLECT_PARENT_POINTERS, bool USE_MIN_TRANSFER_TIMES = false>
class DijkstraInitializedRAPTOR {

public:
    static constexpr bool TargetPruning = TARGET_PRUNING;
    using Profiler = PROFILER;
    static constexpr bool CollectParentPointers = COLLECT_PARENT_POINTERS;
    static constexpr bool UseMinTransferTimes = USE_MIN_TRANSFER_TIMES;
    static constexpr int RoundFactor = UseMinTransferTimes ? 2 : 1;
    using Type = DijkstraInitializedRAPTOR<TargetPruning, Profiler, CollectParentPointers, UseMinTransferTimes>;

public:
    struct EarliestArrivalLabel {
        EarliestArrivalLabel() : arrivalTime(never), parentDepartureTime(never), parent(noVertex), usesRoute(false), routeId(noRouteId) {}
        int arrivalTime;
        int parentDepartureTime;
        Vertex parent;
        bool usesRoute;
        RouteId routeId;
    };
    using Round = std::vector<EarliestArrivalLabel>;

    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel() : arrivalTime(never), rootParent(noVertex), directParent(noVertex) {}
        int arrivalTime;
        Vertex rootParent;
        Vertex directParent;
        inline bool hasSmallerKey(const DijkstraLabel* const other) const noexcept {
            return arrivalTime < other->arrivalTime;
        }
    };

    struct ParentPointer {
        ParentPointer(const Vertex from, const Vertex to) :
            from(from),
            to(to) {
        }
        Vertex from;
        Vertex to;
    };

public:
    DijkstraInitializedRAPTOR(const Data& data, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        roundIndex(-1),
        stopsUpdatedByRoute(data.numberOfStops()),
        stopsUpdatedByTransfer(data.numberOfStops()),
        routesServingUpdatedStops(data.numberOfRoutes()),
        sourceVertex(noVertex),
        targetVertex(noVertex),
        sourceDepartureTime(intMax),
        targetArrivalTime(intMax),
        labelByNumberOfTrips(1, std::vector<DijkstraLabel>(data.transferGraph.numVertices())),
        profiler(profilerTemplate) {
        if constexpr (UseMinTransferTimes) {
            minChangeTimeGraph = data.minChangeTimeGraph();
            AssertMsg(!data.hasImplicitBufferTimes(), "Either min transfer times have to be used OR departure buffer times have to be implicit!");
        } else {
            AssertMsg(data.hasImplicitBufferTimes(), "Either min transfer times have to be used OR departure buffer times have to be implicit!");
        }
        profiler.registerExtraRounds({EXTRA_ROUND_CLEAR, EXTRA_ROUND_INITIALIZATION});
        profiler.registerPhases({PHASE_INITIALIZATION, PHASE_COLLECT, PHASE_SCAN, PHASE_TRANSFERS});
        profiler.registerMetrics({METRIC_ROUTES, METRIC_ROUTE_SEGMENTS, METRIC_VERTICES, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline void run(const Vertex source, const int departureTime, const Vertex target = noVertex, const size_t maxRounds = INFTY) noexcept {
        AssertMsg(data.transferGraph.isVertex(source), "source (" << source << ") is not a vertex!");
        profiler.start();
        profiler.startExtraRound(EXTRA_ROUND_CLEAR);
        clear();
        profiler.doneRound();

        profiler.startExtraRound(EXTRA_ROUND_INITIALIZATION);
        profiler.startPhase();
        initialize(source, departureTime, target);
        profiler.donePhase(PHASE_INITIALIZATION);
        profiler.startPhase();
        dijkstra();
        profiler.donePhase(PHASE_TRANSFERS);
        profiler.doneRound();
        profiler.startPhase();
        collectRoutesServingUpdatedStops();
        profiler.donePhase(PHASE_COLLECT);

        for (size_t i = 0; i < maxRounds; i++) {
            profiler.startRound();
            profiler.startPhase();
            startNewRound();
            profiler.donePhase(PHASE_INITIALIZATION);
            profiler.startPhase();
            scanRoutes();
            profiler.donePhase(PHASE_SCAN);
            if (stopsUpdatedByRoute.empty()) {
                profiler.doneRound();
                break;
            }
            if constexpr (!UseMinTransferTimes) {
                profiler.startPhase();
                startNewRound();
                profiler.donePhase(PHASE_INITIALIZATION);
            }
            profiler.startPhase();
            relaxTransfers();
            profiler.donePhase(PHASE_TRANSFERS);
            profiler.startPhase();
            collectRoutesServingUpdatedStops();
            profiler.donePhase(PHASE_COLLECT);
            profiler.doneRound();
            if (routesServingUpdatedStops.empty()) break;
        }
        profiler.done();
    }

    inline std::vector<Journey> getJourneys() const noexcept {
        AssertMsg(data.transferGraph.isVertex(targetVertex), "No target vertex has been specified!");
        return getJourneys(targetVertex);
    }

    inline std::vector<Journey> getJourneys(const Vertex vertex) const noexcept {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        std::vector<Journey> journeys;
        for (size_t i = 0; i <= roundIndex; i += RoundFactor) {
            getJourney(journeys, i, vertex);
        }
        return journeys;
    }

    inline std::vector<ArrivalLabel> getArrivals() const noexcept {
        AssertMsg(data.transferGraph.isVertex(targetVertex), "No target vertex has been specified!");
        return getArrivals(targetVertex);
    }

    inline std::vector<ArrivalLabel> getArrivals(const Vertex vertex) const noexcept {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        std::vector<ArrivalLabel> labels;
        for (size_t i = 0; i <= roundIndex; i += RoundFactor) {
            getArrival(labels, i, vertex);
        }
        return labels;
    }

    inline bool reachable(const Vertex vertex) const noexcept {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return getEarliestArrivalTime(vertex) < never;
    }

    inline int getEarliestArrivalTime(const Vertex vertex) const noexcept {
        AssertMsg(roundIndex > 0, "The algorithm has not been executed yet!");
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                return std::min(rounds[roundIndex][vertex].arrivalTime, rounds[roundIndex - 1][vertex].arrivalTime);
            } else {
                return rounds[roundIndex][vertex].arrivalTime;
            }
        } else {
            return labelByNumberOfTrips.back()[vertex].arrivalTime;
        }
    }

    inline int getWalkingArrivalTime() const noexcept {
        AssertMsg(roundIndex > 0, "The algorithm has not been executed yet!");
        return labelByNumberOfTrips[0][targetVertex].arrivalTime;
    }

    inline int getWalkingTravelTime() const noexcept {
        AssertMsg(roundIndex > 0, "The algorithm has not been executed yet!");
        return labelByNumberOfTrips[0][targetVertex].arrivalTime - sourceDepartureTime;
    }

    inline int getWalkingTravelTime(const Vertex vertex) const noexcept {
        AssertMsg(roundIndex > 0, "The algorithm has not been executed yet!");
        return labelByNumberOfTrips[0][vertex].arrivalTime - sourceDepartureTime;
    }

    inline std::vector<Vertex> getPath(const Vertex vertex) const {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return journeyToPath(getJourneys(vertex).back());
    }

    inline std::vector<std::string> getRouteDescription(const Vertex vertex) const {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        return data.journeyToText(getJourneys(vertex).back());
    }

    inline int getArrivalTime(const Vertex vertex, const size_t numberOfTrips) const noexcept {
        AssertMsg(data.transferGraph.isVertex(vertex), "The Id " << vertex << " does not correspond to any vertex!");
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                return std::min(rounds[numberOfTrips * 2][vertex].arrivalTime, rounds[(numberOfTrips * 2) + 1][vertex].arrivalTime);
            } else {
                return rounds[numberOfTrips][vertex].arrivalTime;
            }
        } else {
            return labelByNumberOfTrips[numberOfTrips][vertex].arrivalTime;
        }
    }

    inline const std::vector<ParentPointer>& getParentPointers() const noexcept {
        return parentPointers;
    }

    template<bool RESET_CAPACITIES = true>
    inline void clear() noexcept {
        roundIndex = -1;
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        sourceVertex = noVertex;
        targetVertex = noVertex;
        sourceDepartureTime = intMax;
        targetArrivalTime = intMax;
        queue.clear();

        if constexpr (RESET_CAPACITIES) {
            std::vector<Round>().swap(rounds);
            for (std::vector<DijkstraLabel>& label : labelByNumberOfTrips) {
                std::vector<DijkstraLabel>(label.size()).swap(label);
            }
            std::vector<ParentPointer>().swap(parentPointers);
        } else {
            rounds.clear();
            for (std::vector<DijkstraLabel>& label : labelByNumberOfTrips) {
                Vector::fill(label);
            }
            parentPointers.clear();
        }
    }

    inline void reset() noexcept {
        clear<true>();
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void initialize(const Vertex source, const int departureTime, const Vertex target) noexcept {
        sourceVertex = source;
        targetVertex = target;
        sourceDepartureTime = departureTime;
        startNewRound();
        if (data.isStop(source)) arrivalByRoute(StopId(source), departureTime, source, source, departureTime, noRouteId);
        if constexpr (UseMinTransferTimes) startNewRound();
        arrivalByEdge(source, departureTime, source, source);
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            for (const RouteSegment& route : data.routesContainingStop(stop)) {
                AssertMsg(data.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                AssertMsg(data.stopIds[data.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == data.numberOfStopsInRoute(route.routeId)) continue;
                if (data.lastTripOfRoute(route.routeId)[route.stopIndex].departureTime < previousLabel(stop).arrivalTime) continue;
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
            AssertMsg(stopIndex < tripSize - 1, "Cannot scan a route starting at/after the last stop (Route: " << route << ", StopIndex: " << stopIndex << ", TripSize: " << tripSize << ", RoundIndex: " << roundIndex << ")!");

            const StopId* stops = data.stopArrayOfRoute(route);
            const StopEvent* trip = data.lastTripOfRoute(route);
            StopId stop = stops[stopIndex];
            AssertMsg(trip[stopIndex].departureTime >= previousLabel(stop).arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << stop << ", StopIndex: " << stopIndex << ", Time: " << previousLabel(stop).arrivalTime << ", LastDeparture: " << trip[stopIndex].departureTime << ", RoundIndex: " << roundIndex << ")!");

            StopIndex parentIndex = stopIndex;
            const StopEvent* firstTrip = data.firstTripOfRoute(route);
            while (stopIndex < tripSize - 1) {
                const int earliestDepartureTime = previousLabel(stop).arrivalTime;
                while ((trip > firstTrip) && ((trip - tripSize + stopIndex)->departureTime >= earliestDepartureTime)) {
                    trip -= tripSize;
                    parentIndex = stopIndex;
                }
                stopIndex++;
                stop = stops[stopIndex];
                profiler.countMetric(METRIC_ROUTE_SEGMENTS);
                arrivalByRoute(stop, trip[stopIndex].arrivalTime, stops[parentIndex], stops[stopIndex - 1], trip[parentIndex].departureTime, route);
            }
        }
    }

    inline void relaxTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        Assert(queue.empty());
        for (const StopId stop : stopsUpdatedByRoute) {
            const EarliestArrivalLabel& label = UseMinTransferTimes ? previousLabel(stop) : currentLabel(stop);
            if constexpr (UseMinTransferTimes) {
                for (Edge edge : data.transferGraph.edgesFrom(stop)) {
                    arrivalByEdge(data.transferGraph.get(ToVertex, edge), label.arrivalTime + data.transferGraph.get(TravelTime, edge), stop, stop);
                }
                for (Edge edge : minChangeTimeGraph.edgesFrom(stop)) {
                    arrivalByEdge(minChangeTimeGraph.get(ToVertex, edge), label.arrivalTime + minChangeTimeGraph.get(TravelTime, edge), stop, stop);
                }
                arrivalByEdge(stop, label.arrivalTime + data.minTransferTime(stop), stop, stop);
            } else {
                arrivalByEdge(stop, label.arrivalTime, stop, stop);
            }
        }
        dijkstra();
    }

    inline void dijkstra() noexcept {
        while (!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(labelByNumberOfTrips[roundIndex / RoundFactor][0]));
            for (Edge edge : data.transferGraph.edgesFrom(u)) {
                profiler.countMetric(METRIC_EDGES);
                const Vertex v = data.transferGraph.get(ToVertex, edge);
                if (v == uLabel->rootParent) continue;
                arrivalByEdge(v, uLabel->arrivalTime + data.transferGraph.get(TravelTime, edge), uLabel->rootParent, u);
            }
            if constexpr (UseMinTransferTimes) {
                for (Edge edge : minChangeTimeGraph.edgesFrom(u)) {
                    profiler.countMetric(METRIC_EDGES);
                    const Vertex v = minChangeTimeGraph.get(ToVertex, edge);
                    if (v == uLabel->rootParent) continue;
                    arrivalByEdge(v, uLabel->arrivalTime + minChangeTimeGraph.get(TravelTime, edge), uLabel->rootParent, u);
                }
            }
            arrivalByTransfer(u, *uLabel);
            profiler.countMetric(METRIC_VERTICES);
        }
    }

private:
    inline Round& currentRound() noexcept {
        AssertMsg(roundIndex < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        return rounds[roundIndex];
    }

    inline EarliestArrivalLabel& currentLabel(const StopId stop) noexcept {
        AssertMsg(roundIndex < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        return rounds[roundIndex][stop];
    }

    inline const Round& previousRound() const noexcept {
        AssertMsg(roundIndex - 1 < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        AssertMsg(roundIndex > 0, "Cannot return previous round, because no round exists!");
        return rounds[roundIndex - 1];
    }

    inline const EarliestArrivalLabel& previousLabel(const StopId stop) const noexcept {
        AssertMsg(roundIndex - 1 < rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        AssertMsg(roundIndex > 0, "Cannot return previous round, because no round exists!");
        return rounds[roundIndex - 1][stop];
    }

    inline void startNewRound() noexcept {
        AssertMsg(roundIndex + 1 <= rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        roundIndex++;
        if (roundIndex == rounds.size()) {
            if (rounds.size() < RoundFactor) {
                rounds.emplace_back(data.numberOfStops());
            } else {
                rounds.emplace_back(rounds[rounds.size() - RoundFactor]);
            }
        }
        if (roundIndex / RoundFactor == labelByNumberOfTrips.size()) {
            labelByNumberOfTrips.emplace_back(labelByNumberOfTrips.back());
        }
    }

    inline void arrivalByRoute(const StopId stop, const int arrivalTime, const Vertex rootParent, const Vertex directParent, const int parentDepartureTime, const RouteId route) noexcept {
        AssertMsg(data.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], Stop: " << stop << ")!");
        if constexpr (TargetPruning) {if (arrivalTime >= targetArrivalTime) return;}
        EarliestArrivalLabel& label = currentLabel(stop);
        if ((arrivalTime >= label.arrivalTime) || ((roundIndex > 0) && (arrivalTime >= previousLabel(stop).arrivalTime))) return;
        label.arrivalTime = arrivalTime;
        label.parentDepartureTime = parentDepartureTime;
        label.parent = rootParent;
        label.usesRoute = true;
        label.routeId = route;
        for (size_t i = roundIndex + RoundFactor; i < rounds.size(); i += RoundFactor) {
            if (rounds[i][stop].arrivalTime <= arrivalTime) break;
            rounds[i][stop].arrivalTime = arrivalTime;
        }
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        stopsUpdatedByRoute.insert(stop);
        if constexpr (CollectParentPointers) {
            collectParentPointer(directParent, stop);
        } else {
            suppressUnusedParameterWarning(directParent);
        }
        if constexpr (TargetPruning) improveTargetArrivalTime(stop, arrivalTime);
    }

    inline void arrivalByTransfer(const Vertex vertex, const DijkstraLabel& dijkstraLabel) noexcept {
        if constexpr (UseMinTransferTimes) {
            AssertMsg(roundIndex > 0, "arrivalByTransfer can not be used in the first round!");
        }
        AssertMsg(dijkstraLabel.arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(dijkstraLabel.arrivalTime) << " [" << dijkstraLabel.arrivalTime << "], Vertex: " << vertex << ")!");
        if (data.isStop(vertex)) {
            if constexpr (TargetPruning) {if (dijkstraLabel.arrivalTime >= targetArrivalTime) return;}
            EarliestArrivalLabel& label = currentLabel(StopId(vertex));
            if (dijkstraLabel.arrivalTime >= label.arrivalTime) return;
            label.arrivalTime = dijkstraLabel.arrivalTime;
            label.parentDepartureTime = getParentDepartureTime(dijkstraLabel.rootParent, roundIndex - 1);
            label.parent = dijkstraLabel.rootParent;
            label.usesRoute = false;
            for (size_t i = roundIndex + RoundFactor; i < rounds.size(); i += RoundFactor) {
                if (rounds[i][vertex].arrivalTime <= dijkstraLabel.arrivalTime) break;
                rounds[i][vertex].arrivalTime = dijkstraLabel.arrivalTime;
            }
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            stopsUpdatedByTransfer.insert(StopId(vertex));
        }
        if constexpr (CollectParentPointers) collectParentPointer(dijkstraLabel.directParent, vertex);
        if constexpr (TargetPruning) improveTargetArrivalTime(vertex, dijkstraLabel.arrivalTime);
    }

    inline void arrivalByEdge(const Vertex vertex, const int arrivalTime, const Vertex rootParent, const Vertex directParent) noexcept {
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], Vertex: " << vertex << ")!");
        if constexpr (TargetPruning) {if (arrivalTime >= targetArrivalTime) return;}
        DijkstraLabel& label = labelByNumberOfTrips[roundIndex / RoundFactor][vertex];
        if (arrivalTime >= label.arrivalTime) return;
        label.rootParent = rootParent;
        label.directParent = directParent;
        for (size_t i = roundIndex / RoundFactor; i < labelByNumberOfTrips.size(); i++) {
            if (labelByNumberOfTrips[i][vertex].arrivalTime <= arrivalTime) break;
            labelByNumberOfTrips[i][vertex].arrivalTime = arrivalTime;
        }
        queue.update(&label);
        if constexpr (TargetPruning) improveTargetArrivalTime(vertex, arrivalTime);
    }

    inline int getParentDepartureTime(const Vertex vertex, const size_t round) const noexcept {
        AssertMsg(data.isStop(vertex) || vertex == sourceVertex, "Vertex " << vertex << " should not be a root parent!");
        if (vertex == sourceVertex) {
            return sourceDepartureTime;
        } else {
            return rounds[round][vertex].arrivalTime;
        }
    }

    inline void improveTargetArrivalTime(const Vertex vertex, const int arrivalTime) noexcept {
        if (vertex != targetVertex) return;
        if (targetArrivalTime <= arrivalTime) return;
        targetArrivalTime = arrivalTime;
    }

    inline void collectParentPointer(const Vertex from, const Vertex to) noexcept {
        if (from == to) return;
        parentPointers.emplace_back(from, to);
    }

    inline void getJourney(std::vector<Journey>& journeys, size_t round, Vertex vertex) const noexcept {
        if constexpr (UseMinTransferTimes) {
            AssertMsg(round % 2 == 0, "Journeys can only be computed for even rounds! (" << round << ")!");
        }
        Journey journey;
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                if ((round < roundIndex) && (rounds[round + 1][vertex].arrivalTime < rounds[round][vertex].arrivalTime)) round++;
            }
            if (rounds[round][vertex].arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        } else {
            const DijkstraLabel& label = labelByNumberOfTrips[round / RoundFactor][vertex];
            if (label.arrivalTime >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
            journey.emplace_back(label.rootParent, vertex, getParentDepartureTime(label.rootParent, round), label.arrivalTime, false, noRouteId);
            vertex = label.rootParent;
        }
        while (vertex != sourceVertex) {
            AssertMsg(round != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = rounds[round][vertex];
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

    inline void getArrival(std::vector<ArrivalLabel>& labels, size_t round, const Vertex vertex) const noexcept {
        if constexpr (UseMinTransferTimes) {
            AssertMsg(round % 2 == 0, "Arrivals can only be computed for even rounds! (" << round << ")!");
        }
        const size_t numberOfTrips = round / RoundFactor;
        if (data.isStop(vertex)) {
            if constexpr (UseMinTransferTimes) {
                if ((round < roundIndex) && (rounds[round + 1][vertex].arrivalTime < rounds[round][vertex].arrivalTime)) round++;
            }
            const EarliestArrivalLabel& label = rounds[round][vertex];
            if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
            labels.emplace_back(label.arrivalTime, numberOfTrips);
        } else {
            const DijkstraLabel& label = labelByNumberOfTrips[numberOfTrips][vertex];
            if (label.arrivalTime >= (labels.empty() ? never : labels.back().arrivalTime)) return;
            labels.emplace_back(label.arrivalTime, numberOfTrips);
        }
    }

    inline void printRoundsForStop(const StopId stop) const noexcept {
        Assert(data.isStop(stop));
        std::cout << "Raptor Label for stop " << stop << ":" << std::endl;
        std::cout << std::setw(10) << "Round" << std::setw(14) << "arrivalTime" << std::setw(14) << "parent" << std::endl;
        for (size_t i = 0; i <= roundIndex; i++) {
            std::cout << std::setw(10) << i << std::setw(14) << rounds[i][stop].arrivalTime << std::setw(14) << rounds[i][stop].parent << std::endl;
        }
    }

private:
    const Data& data;
    TransferGraph minChangeTimeGraph;

    std::vector<Round> rounds;
    size_t roundIndex;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    Vertex targetVertex;
    int sourceDepartureTime;
    int targetArrivalTime;

    std::vector<std::vector<DijkstraLabel>> labelByNumberOfTrips;
    ExternalKHeap<2, DijkstraLabel> queue;

    Profiler profiler;

    std::vector<ParentPointer> parentPointers;

};

}
