#pragma once

#include <vector>

#include "../../../Dijkstra/Dijkstra.h"
#include "../../../../DataStructures/Container/ExternalKHeap.h"
#include "../../../../DataStructures/Container/Set.h"
#include "../../../../DataStructures/Container/Map.h"
#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/Entities/TripIterator.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/PreprocessingJourney.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../Helpers/Vector/Vector.h"
#include "ForwardSearch.h"

namespace RAPTOR::RouteFlags {

template<typename FORWARD_SEARCH, typename FLAGS_MARKER>
class BackwardSearch {

public:
    using ForwardSearchType = FORWARD_SEARCH;
    using FlagsMarkerType = FLAGS_MARKER;
    using Type = BackwardSearch<ForwardSearchType, FlagsMarkerType>;
    using RouteFlagsDataType = typename ForwardSearchType::RouteFlagsDataType;

public:
    struct StackElement {
        Vertex vertex;
        size_t tripsLeft;
        int arrivalTime;
        bool usesRoute;
    };

    struct ChildPointer {
        int departureTime;
        int travelTime;
        Vertex child;
        const StopEvent* departureEvent;
        const StopEvent* arrivalEvent;

        inline bool usesRoute() const noexcept {
            return departureEvent != NULL;
        }

        friend bool operator>(const ChildPointer& pointerA, const ChildPointer& pointerB) noexcept {
            return pointerA.departureTime < pointerB.departureTime;
        }

        friend bool operator>=(const ChildPointer& pointerA, const ChildPointer& pointerB) noexcept {
            return pointerA.departureTime <= pointerB.departureTime;
        }

        friend bool isSameConnection(const ChildPointer& pointerA, const ChildPointer& pointerB) noexcept {
            return pointerA.child == pointerB.child && pointerA.departureEvent == pointerB.departureEvent;
        }

    };

    using Children = std::vector<ChildPointer>;
    using ChildrenRound = std::vector<Children>;

    struct EarliestArrivalLabel {
        EarliestArrivalLabel() :
            arrivalTime(never),
            parent(noVertex),
            route(noRouteId),
            departureEvent(NULL),
            arrivalEvent(NULL),
            timestamp(0) {
        }

        inline void reset(const int time) noexcept {
            arrivalTime = never;
            parent = noVertex;
            route = noRouteId;
            departureEvent = NULL;
            arrivalEvent = NULL;
            timestamp = time;
        }

        inline bool usesRoute() const noexcept {
            return departureEvent != NULL;
        }

        int arrivalTime;
        Vertex parent;
        RouteId route;
        const StopEvent* departureEvent;
        const StopEvent* arrivalEvent;
        int timestamp;
    };
    using Round = std::vector<EarliestArrivalLabel>;

    struct DijkstraLabel : public ExternalKHeapElement {
        DijkstraLabel() : ExternalKHeapElement(), distance(intMax), timestamp(-1) {}

        inline void reset(int time) {
            distance = intMax;
            timestamp = time;
        }

        inline bool hasSmallerKey(const DijkstraLabel* other) const {
            return distance < other->distance;
        }

        int distance;
        int timestamp;
    };

    struct FinishedStop {
        StopId stop;
        size_t completedTrips;
    };

public:
    BackwardSearch(const RouteFlagsDataType& flagsData, ForwardSearchType& forwardSearch, FlagsMarkerType& flagsMarker) :
        raptorData(flagsData.getRaptorData(FORWARD)),
        shortcutGraph(flagsData.getShortcutGraph(FORWARD)),
        flagsData(flagsData),
        forwardSearch(forwardSearch),
        flagsMarker(flagsMarker),
        roundIndex(-1),
        queue(raptorData.transferGraph.numVertices()),
        dijkstraLabels(raptorData.transferGraph.numVertices()),
        stopsUpdatedByRoute(raptorData.numberOfStops()),
        stopsUpdatedByTransfer(raptorData.numberOfStops()),
        routesServingUpdatedStops(raptorData.numberOfRoutes()),
        sourceVertex(noVertex),
        targetVertex(noVertex),
        sourceDepartureTime(intMax),
        maxTrips(-1),
        forceFinalTrip(false),
        finalStopEvent(-1),
        timestamp(0) {
        AssertMsg(raptorData.hasImplicitBufferTimes(), "Route-Flags requires implicit buffer times!");
    }

    inline void setTarget(const Vertex target) noexcept {
        targetVertex = target;
        usedStopEvent.clear();
        tempUsedStopEvent.clear();
        rounds.clear();
        children.clear();
        tempChildren.clear();
        forceFinalTrip = false;
        finalStopEvent = -1;
    }

    inline void setForceFinalTrip(const bool force, const size_t stopEvent) noexcept {
        forceFinalTrip = force;
        if (forceFinalTrip) {
            tempUsedStopEvent.clear();
            tempChildren.clear();
        }
        finalStopEvent = stopEvent;
    }

    inline void run(const Vertex source, const int departureTime, const size_t numTrips) noexcept {
        clear();
        initialize(source, departureTime, numTrips);
        relaxInitialTransfers();
        for (size_t i = 0; (i < maxTrips) && (!stopsUpdatedByTransfer.empty()); i++) {
            startNewRound();
            collectRoutesServingUpdatedStops();
            scanRoutes();
            if (stopsUpdatedByRoute.empty() || i == maxTrips - 1) break;
            relaxIntermediateTransfers();
        }
        for (const StopId stop : stopsUpdatedByRoute) {
            AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
            finishedStops.emplace_back(FinishedStop{stop, roundIndex});
        }

        for (const FinishedStop& finishedStop : finishedStops) {
            flagsMarker.markJourney(getJourney(finishedStop), flagsData.getFinePartition().cellsOfVertex(source));
        }
    }

private:
    inline void clear() noexcept {
        roundIndex = -1;
        sourceVertex = noVertex;
        sourceDepartureTime = intMax;
        maxTrips = -1;
        stopsUpdatedByRoute.clear();
        stopsUpdatedByTransfer.clear();
        routesServingUpdatedStops.clear();
        finishedStops.clear();
        timestamp++;
    }

    inline void initialize(const Vertex source, const int departureTime, const size_t numTrips) noexcept {
        sourceVertex = source;
        sourceDepartureTime = departureTime;
        maxTrips = numTrips;
        while (getChildren().size() <= numTrips) {
            getChildren().emplace_back(raptorData.transferGraph.numVertices());
        }
        while (getUsedStopEvent().size() <= numTrips) {
            getUsedStopEvent().emplace_back(raptorData.numberOfStopEvents());
        }
        startNewRound();
        if (raptorData.isStop(source)) {
            stopsUpdatedByRoute.insert(StopId(source));
            currentRoundLabel(StopId(source)).parent = source;
        }
    }

    inline void collectRoutesServingUpdatedStops() noexcept {
        for (const StopId stop : stopsUpdatedByTransfer) {
            AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
            const int arrivalTime = previousRoundLabel(stop).arrivalTime;
            AssertMsg(arrivalTime < never, "Updated stop " << stop << ", has arrival time = never!");
            bool hasUsedDepartingTrip = false;
            for (const RouteSegment& route : raptorData.routesContainingStop(stop)) {
                AssertMsg(raptorData.isRoute(route.routeId), "Route " << route.routeId << " is out of range!");
                AssertMsg(raptorData.stopIds[raptorData.firstStopIdOfRoute[route.routeId] + route.stopIndex] == stop, "RAPTOR data contains invalid route segments!");
                if (route.stopIndex + 1 == raptorData.numberOfStopsInRoute(route.routeId)) continue;
                const StopEvent* firstTrip = findFirstTripAfter(route.routeId, route.stopIndex, arrivalTime);
                if (!firstTrip || firstTrip->departureTime > getForwardArrivalTime(stop, roundIndex - 1)) continue;
                if (isUsedStopEvent(firstTrip - &(raptorData.stopEvents[0]), maxTrips + 1 - roundIndex)) {
                    hasUsedDepartingTrip = true;
                    continue;
                }
                if (routesServingUpdatedStops.contains(route.routeId)) {
                    routesServingUpdatedStops[route.routeId] = std::min(routesServingUpdatedStops[route.routeId], route.stopIndex);
                } else {
                    routesServingUpdatedStops.insert(route.routeId, route.stopIndex);
                }
            }
            if (hasUsedDepartingTrip) {
                finishedStops.emplace_back(FinishedStop{stop, roundIndex - 1});
            }
        }
    }

    inline const StopEvent* findFirstTripAfter(const RouteId route, const StopIndex stopIndex, const int arrivalTime) const noexcept {
        TripIterator tripIterator = raptorData.getTripIterator(route, stopIndex);
        if (tripIterator.departureTime() < arrivalTime) return NULL;
        while (tripIterator.hasEarlierTrip() && tripIterator.previousDepartureTime() >= arrivalTime) {
            tripIterator.previousTrip();
        }
        return tripIterator.stopEvent();
    }

    inline void scanRoutes() noexcept {
        stopsUpdatedByRoute.clear();
        for (const RouteId route : routesServingUpdatedStops.getKeys()) {
            TripIterator tripIterator = raptorData.getTripIterator(route, routesServingUpdatedStops[route]);
            AssertMsg(tripIterator.departureTime() >= previousRoundLabel(tripIterator.stop()).arrivalTime, "Cannot scan a route after the last trip has departed (Route: " << route << ", Stop: " << tripIterator.stop() << ", StopIndex: " << tripIterator.getStopIndex() << ", Time: " << previousRoundLabel(tripIterator.stop()).arrivalTime << ", LastDeparture: " << tripIterator.departureTime() << ")!");

            StopIndex parentIndex = tripIterator.getStopIndex();
            bool start = true;
            while (tripIterator.hasFurtherStops()) {
                const int arrivalTime = previousRoundLabel(tripIterator.stop()).arrivalTime;
                if (tripIterator.hasEarlierTrip() && tripIterator.previousDepartureTime() >= arrivalTime) {
                    do {
                        tripIterator.previousTrip();
                    } while (tripIterator.hasEarlierTrip() && tripIterator.previousDepartureTime() >= arrivalTime);
                    useStopEvent(tripIterator.stopEvent() - &(raptorData.stopEvents[0]), maxTrips + 1 - roundIndex);
                    parentIndex = tripIterator.getStopIndex();
                } else if (start) {
                    useStopEvent(tripIterator.stopEvent() - &(raptorData.stopEvents[0]), maxTrips + 1 - roundIndex);
                }
                start = false;
                tripIterator.nextStop();
                const bool isFinalStopEvent = forceFinalTrip && size_t(tripIterator.stopEvent() - &(raptorData.stopEvents[0])) == finalStopEvent;
                if (!isFinalStopEvent && getTargetArrivalTime() < tripIterator.arrivalTime()) continue;
                if (arrivalByRoute(tripIterator.stop(), tripIterator.arrivalTime(), isFinalStopEvent)) {
                    EarliestArrivalLabel& label = currentRoundLabel(tripIterator.stop());
                    label.parent = tripIterator.stop(parentIndex);
                    label.route = route;
                    label.departureEvent = tripIterator.stopEvent(parentIndex);
                    label.arrivalEvent = tripIterator.stopEvent();
                }
            }
        }
    }

    inline void relaxInitialTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        AssertMsg(stopsUpdatedByRoute.size() == (raptorData.isStop(sourceVertex) ? 1 : 0), "stopsUpdatedByRoute should contain " << (raptorData.isStop(sourceVertex) ? 1 : 0) << " stops but contains " << stopsUpdatedByRoute.size());
        AssertMsg(queue.empty(), "Queue still has " << queue.size() << " elements!");
        DijkstraLabel& sourceLabel = getDijkstraLabel(sourceVertex);
        sourceLabel.distance = sourceDepartureTime;
        queue.update(&sourceLabel);
        while(!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(dijkstraLabels[0]));
            if (raptorData.isStop(u)) {
                stopsUpdatedByTransfer.insert(StopId(u));
                EarliestArrivalLabel& label = currentRoundLabel(StopId(u));
                label.arrivalTime = uLabel->distance;
                label.parent = sourceVertex;
                label.departureEvent = NULL;
                label.arrivalEvent = NULL;
            }
            if (u == targetVertex) {
                queue.clear();
                break;
            }
            for (const Edge edge : raptorData.transferGraph.edgesFrom(u)) {
                const Vertex v = raptorData.transferGraph.get(ToVertex, edge);
                DijkstraLabel& vLabel = getDijkstraLabel(v);
                const int distance = uLabel->distance + raptorData.transferGraph.get(TravelTime, edge);
                if (vLabel.distance > distance && distance <= getForwardArrivalTime(v, 0)) {
                    AssertMsg(distance == getForwardArrivalTime(v, 0), "Found shorter path! Source: " << sourceVertex << ", target: " << targetVertex << ", stop: " << v << ", Source departure time: " << sourceDepartureTime << ", Target arrival time: " << getTargetArrivalTime() << ", max trips: " << maxTrips << ", Distance: " << distance << ", forward: " << getForwardArrivalTime(v, 0));
                    vLabel.distance = distance;
                    queue.update(&vLabel);
                }
            }
        }
    }

    inline void relaxIntermediateTransfers() noexcept {
        routesServingUpdatedStops.clear();
        stopsUpdatedByTransfer.clear();
        for (const StopId stop : stopsUpdatedByRoute) {
            stopsUpdatedByTransfer.insert(stop);
            const int departureTime = currentRoundLabel(stop).arrivalTime;
            for (const Edge edge : shortcutGraph.edgesFrom(stop)) {
                const int arrivalTime = departureTime + shortcutGraph.get(TravelTime, edge);
                AssertMsg(raptorData.isStop(shortcutGraph.get(ToVertex, edge)), "Graph contains edges to non stop vertices!");
                const StopId toStop = StopId(shortcutGraph.get(ToVertex, edge));
                if (arrivalByTransfer(toStop, arrivalTime)) {
                    EarliestArrivalLabel& label = currentRoundLabel(toStop);
                    label.parent = stop;
                    label.departureEvent = NULL;
                    label.arrivalEvent = NULL;
                }
            }
        }
    }

    inline bool arrivalByRoute(const StopId stop, const int arrivalTime, const bool isFinalStopEvent) noexcept {
        AssertMsg(roundIndex > 0, "arrivalByRoute cannot be used in the first round!");
        AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by route BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "], stop: " << stop << ")!");

        if (previousRoundLabel(stop).arrivalTime <= arrivalTime) return false;
        EarliestArrivalLabel& currentLabel = currentRoundLabel(stop);
        if (!isFinalStopEvent) {
            if (currentLabel.arrivalTime <= arrivalTime) return false;
            if (getTargetArrivalTime() < arrivalTime) return false;
            if (getForwardArrivalTime(stop, roundIndex) < arrivalTime) return false;
        }
        currentLabel.arrivalTime = arrivalTime;

        stopsUpdatedByRoute.insert(stop);
        return true;
    }

    inline bool arrivalByTransfer(const StopId stop, const int arrivalTime) noexcept {
        AssertMsg(raptorData.isStop(stop), "Stop " << stop << " is out of range!");
        AssertMsg(arrivalTime >= sourceDepartureTime, "Arriving by edge BEFORE departing from the source (source departure time: " << String::secToTime(sourceDepartureTime) << " [" << sourceDepartureTime << "], arrival time: " << String::secToTime(arrivalTime) << " [" << arrivalTime << "])!");

        EarliestArrivalLabel& currentLabel = currentRoundLabel(stop);
        if (currentLabel.arrivalTime <= arrivalTime || getForwardArrivalTime(stop, roundIndex) < arrivalTime) return false;
        if (getTargetArrivalTime() < arrivalTime) return false;
        currentLabel.arrivalTime = arrivalTime;

        stopsUpdatedByTransfer.insert(stop);
        return true;
    }

    inline EarliestArrivalLabel& currentRoundLabel(const StopId stop) noexcept {
        return getRoundLabel(roundIndex, stop);
    }

    inline EarliestArrivalLabel& previousRoundLabel(const StopId stop) noexcept {
        return getRoundLabel(roundIndex - 1, stop);
    }

    inline EarliestArrivalLabel& getRoundLabel(const size_t round, const StopId stop) noexcept {
        AssertMsg(round < rounds.size(), "Round index is out of bounds (round = " << round << ", rounds.size() = " << rounds.size() << ")!");
        EarliestArrivalLabel& label = rounds[round][stop];
        if (label.timestamp != timestamp) {
            if (round > 0) {
                label = getRoundLabel(round - 1, stop);
            } else {
                label = EarliestArrivalLabel();
            }
            label.timestamp = timestamp;
        }
        return label;
    }

    inline DijkstraLabel& getDijkstraLabel(const Vertex vertex) noexcept {
        DijkstraLabel& label = dijkstraLabels[vertex];
        if (label.timestamp != timestamp) {
            label.reset(timestamp);
        }
        return label;
    }

    inline void startNewRound() noexcept {
        AssertMsg(roundIndex + 1 <= rounds.size(), "Round index is out of bounds (roundIndex = " << roundIndex << ", rounds.size() = " << rounds.size() << ")!");
        roundIndex++;
        if (roundIndex == rounds.size()) {
            if (rounds.empty()) {
                rounds.emplace_back(raptorData.numberOfStops());
            } else {
                rounds.emplace_back(rounds.back());
            }
        }
    }

    inline int getForwardArrivalTime(const Vertex vertex, const size_t numTrips) const noexcept {
        return -forwardSearch.getArrivalTime(vertex, maxTrips - numTrips);
    }

    inline int getTargetArrivalTime() const noexcept {
        return -forwardSearch.getSourceDepartureTime();
    }

    inline void addChild(const Vertex parent, const size_t tripsLeft, const int previousArrivalTime, const int arrivalTime, const Vertex child, const StopEvent* departureEvent, const StopEvent* arrivalEvent) noexcept {
        const int departureTime = (departureEvent) ? departureEvent->departureTime : previousArrivalTime;
        const int travelTime = arrivalTime - departureTime;
        ChildPointer newPointer{departureTime, travelTime, child, departureEvent, arrivalEvent};
        bool add = true;
        Vector::removeIf(getChildren()[tripsLeft][parent], [&](const ChildPointer& childPointer) {
           if (isSameConnection(childPointer, newPointer)) {
               if (childPointer.departureTime >= newPointer.departureTime) {
                   add = false;
                   return false;
               } else {
                   return true;
               }
           }
           return false;
        });

        if (add) {
            Vector::insertSorted<ChildPointer, true>(getChildren()[tripsLeft][parent], newPointer);
        }
    }

    //TODO: Remove child pointers that have become obsolete?
    inline bool extractChildren(PreprocessingJourney& journey, const Vertex vertex, const size_t tripsLeft, const int arrivalTime) noexcept {
        if (vertex == targetVertex) return true;
        if (tripsLeft == 0) {
            journey.shortcuts.emplace_back(Shortcut{targetVertex, vertex, getTargetArrivalTime() - getForwardArrivalTime(vertex, maxTrips)});
            return true;
        }
        if (getChildren()[tripsLeft][vertex].empty()) return false;
        std::vector<StackElement> stack;
        stack.emplace_back(StackElement{vertex, tripsLeft, arrivalTime, false});

        while (!stack.empty()) {
            StackElement s = stack.back();
            stack.pop_back();
            if (s.vertex == targetVertex) continue;
            if (s.tripsLeft == 0) {
                journey.shortcuts.emplace_back(Shortcut{targetVertex, s.vertex, getTargetArrivalTime() - getForwardArrivalTime(s.vertex, maxTrips)});
                continue;
            }
            for (const ChildPointer& childPointer : getChildren()[s.tripsLeft][s.vertex]) {
                if (!childPointer.usesRoute() && !s.usesRoute) continue;
                if (childPointer.departureTime < s.arrivalTime) break;
                size_t newTripsLeft = s.tripsLeft;
                int newArrivalTime = never;
                if (childPointer.usesRoute()) {
                    journey.departureStopEvents.emplace_back(childPointer.departureEvent);
                    if (!forceFinalTrip || newTripsLeft > 1 || childPointer.child != targetVertex) {
                        journey.arrivalStopEvents.emplace_back(childPointer.arrivalEvent);
                    }
                    newArrivalTime = childPointer.arrivalEvent->arrivalTime;
                    newTripsLeft--;
                } else if (s.vertex != childPointer.child && s.vertex != sourceVertex) {
                    journey.shortcuts.emplace_back(Shortcut{childPointer.child, s.vertex, childPointer.travelTime});
                    newArrivalTime = s.arrivalTime + childPointer.travelTime;
                }
                stack.emplace_back(StackElement{childPointer.child, newTripsLeft, newArrivalTime, childPointer.usesRoute()});
            }
        }
        return true;
    }

    inline std::vector<ChildrenRound>& getChildren() noexcept {
        return forceFinalTrip ? tempChildren : children;
    }

    inline bool isUsedStopEvent(const size_t stopEventId, const size_t numTrips) const noexcept {
        return getUsedStopEvent()[numTrips][stopEventId];
    }

    inline void useStopEvent(const size_t stopEventId, const size_t numTrips) noexcept {
        getUsedStopEvent()[numTrips][stopEventId] = true;
    }

    inline std::vector<std::vector<bool>>& getUsedStopEvent() noexcept {
        return forceFinalTrip ? tempUsedStopEvent : usedStopEvent;
    }

    inline const std::vector<std::vector<bool>>& getUsedStopEvent() const noexcept {
        return forceFinalTrip ? tempUsedStopEvent : usedStopEvent;
    }

    inline PreprocessingJourney getJourney(const FinishedStop& finishedStop) noexcept {
        PreprocessingJourney journeyToTarget;
        Vertex vertex = finishedStop.stop;
        size_t numTrips = finishedStop.completedTrips;

        if (!extractChildren(journeyToTarget, vertex, maxTrips - numTrips, getRoundLabel(numTrips, StopId(vertex)).arrivalTime)) return journeyToTarget;

        PreprocessingJourney journeyFromSource;
        while (vertex != sourceVertex) {
            AssertMsg(raptorData.isStop(vertex), "Vertex is not a stop!");
            AssertMsg(numTrips != size_t(-1), "Backtracking parent pointers did not pass through the source stop!");
            const EarliestArrivalLabel& label = getRoundLabel(numTrips, StopId(vertex));
            AssertMsg(label.parent != noVertex, "Vertex " << vertex << " does not have a parent!");
            const size_t previousNumTrips = label.usesRoute() ? numTrips - 1 : numTrips;
            const int previousArrivalTime = (label.parent == sourceVertex) ? sourceDepartureTime : getRoundLabel(previousNumTrips, StopId(label.parent)).arrivalTime;
            addChild(label.parent, maxTrips - previousNumTrips, previousArrivalTime, label.arrivalTime, vertex, label.departureEvent, label.arrivalEvent);
            if (label.usesRoute()) {
                journeyFromSource.departureStopEvents.emplace_back(label.departureEvent);
                journeyFromSource.arrivalStopEvents.emplace_back(label.arrivalEvent);
            } else if (label.parent != vertex && label.parent != sourceVertex) {
                journeyFromSource.shortcuts.emplace_back(Shortcut{vertex, label.parent, label.arrivalTime - previousArrivalTime});
            }
            vertex = label.parent;
            numTrips = previousNumTrips;
        }

        journeyFromSource.reverse();
        if (forceFinalTrip && journeyToTarget.empty() && !journeyFromSource.arrivalStopEvents.empty()) {
            journeyFromSource.arrivalStopEvents.pop_back();
        }
        journeyFromSource += journeyToTarget;
        return journeyFromSource;
    }

private:
    const Data& raptorData;
    const TransferGraph& shortcutGraph;
    const RouteFlagsDataType& flagsData;
    ForwardSearchType& forwardSearch;
    FlagsMarkerType& flagsMarker;

    std::vector<Round> rounds;
    std::vector<ChildrenRound> children;
    std::vector<ChildrenRound> tempChildren;
    size_t roundIndex;
    ExternalKHeap<2, DijkstraLabel> queue;
    std::vector<DijkstraLabel> dijkstraLabels;

    IndexedSet<false, StopId> stopsUpdatedByRoute;
    IndexedSet<false, StopId> stopsUpdatedByTransfer;
    IndexedMap<StopIndex, false, RouteId> routesServingUpdatedStops;

    Vertex sourceVertex;
    Vertex targetVertex;
    int sourceDepartureTime;
    size_t maxTrips;
    bool forceFinalTrip;
    size_t finalStopEvent;

    int timestamp;

    std::vector<std::vector<bool>> usedStopEvent;
    std::vector<std::vector<bool>> tempUsedStopEvent;
    std::vector<FinishedStop> finishedStops;
};
}
