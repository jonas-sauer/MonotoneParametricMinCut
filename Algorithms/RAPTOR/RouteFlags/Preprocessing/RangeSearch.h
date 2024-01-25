#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <utility>

#include "../../../Dijkstra/Dijkstra.h"

#include "../../../../Helpers/Meta.h"
#include "../../../../Helpers/Helpers.h"
#include "../../../../Helpers/Vector/Vector.h"
#include "../../../../Helpers/Vector/Permutation.h"

#include "../../../../DataStructures/Container/Heap.h"
#include "../../../../DataStructures/Container/Set.h"
#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/ShortcutFlags.h"
#include "../../../../DataStructures/Geometry/Rectangle.h"
#include "../../Profiler.h"
#include "AlternativeBackwardSearch.h"
#include "BackwardSearch.h"
#include "FlagsMarker.h"
#include "ForwardSearch.h"

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE, bool DIRECTION = FORWARD, bool USE_ALTERNATIVE_BACKWARD_SEARCH = false, bool DEBUG = false>
class RangeSearch {

public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    constexpr static bool Direction = DIRECTION;
    constexpr static bool UseAlternativeBackwardSearch = USE_ALTERNATIVE_BACKWARD_SEARCH;
    constexpr static bool Debug = DEBUG;
    using Type = RangeSearch<Flags, FinePartitionType, Direction, UseAlternativeBackwardSearch, Debug>;
    using ForwardSearchType = ForwardSearch<Flags, FinePartitionType, Meta::IF<Debug, SimpleProfiler, NoProfiler>, Direction>;
    using RouteFlagsDataType = RouteFlagsData<Flags, FinePartitionType>;
    using FlagsMarkerType = FlagsMarker<Flags, FinePartitionType, Direction>;
    using BackwardSearchType = Meta::IF<UseAlternativeBackwardSearch, AlternativeBackwardSearch<ForwardSearchType, FlagsMarkerType>, BackwardSearch<ForwardSearchType, FlagsMarkerType>>;

private:
    struct DepartureLabel {
        DepartureLabel(const int departureTime, const StopId departureStop, const bool forceInitialTrip = false, const RouteSegment initialRouteSegment = RouteSegment(), const StopEvent* initialTrip = NULL) :
            departureTime(departureTime),
            departureStop(departureStop),
            forceInitialTrip(forceInitialTrip),
            initialRouteSegment(initialRouteSegment),
            initialTrip(initialTrip) {
        }
        inline bool operator<(const DepartureLabel& other) const noexcept {
            return departureTime > other.departureTime ||
                    (departureTime == other.departureTime && forceInitialTrip < other.forceInitialTrip);
        }

        inline friend std::ostream& operator<<(std::ostream& out, const DepartureLabel& label) noexcept {
            out << "Departure time: " << label.departureTime << ", departureStop: " << label.departureStop << ", force initial trip? " << label.forceInitialTrip << ", route segment: " << label.initialRouteSegment;
            return out;
        }

        inline const StopEvent* initialStopEvent() const noexcept {
            if (!forceInitialTrip) return NULL;
            return initialTrip + initialRouteSegment.stopIndex;
        }

        int departureTime;
        StopId departureStop;
        bool forceInitialTrip;
        RouteSegment initialRouteSegment;
        const StopEvent* initialTrip;
    };

    struct ArrivalAtVertex {
        Vertex vertex;
        int arrivalTime;
        size_t numberOfTrips;

        inline bool operator<(const ArrivalAtVertex& other) const noexcept {
            return arrivalTime < other.arrivalTime;
        }

    };

public:
    RangeSearch(RouteFlagsDataType& flagsData) :
        flagsData(flagsData),
        raptorData(flagsData.getRaptorData(Direction)),
        forwardSearch(flagsData),
        flagsMarker(flagsData),
        backwardSearch(flagsData, forwardSearch, flagsMarker),
        sourceVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never) {
    }

    inline void run(const Vertex source, const int minTime = -never, const int maxTime = never, const size_t maxRounds = INFTY) noexcept {
        clear();
        sourceVertex = source;
        sourceCells = flagsData.getCoarsePartition().cellsOfVertex(sourceVertex);
        flagsMarker.setSource(source, sourceCells);
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        if constexpr (!UseAlternativeBackwardSearch && Direction == BACKWARD) {
            backwardSearch.setTarget(source);
        }
        forwardSearch.template runInitialize<true>(source, 0);
        forwardSearch.runInitialTransfers();
        collectDepartures();
        for (size_t i = 0; i < departures.size(); i++) {
            forwardSearch.template runInitialize<false>(source, departures[i].departureTime, departures[i].forceInitialTrip, departures[i].initialRouteSegment, departures[i].initialTrip);
            forwardSearch.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardSearch.getWalkingTravelTime(departures[i].departureStop));
            if (!departures[i].forceInitialTrip) {
                while (i + 1 < departures.size() && departures[i].departureTime == departures[i + 1].departureTime && !departures[i+1].forceInitialTrip) {
                    i++;
                    forwardSearch.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardSearch.getWalkingTravelTime(departures[i].departureStop));
                }
            }
            if constexpr (Debug) std::cout << "Departure time: " << departures[i].departureTime << std::endl;
            forwardSearch.runRounds(maxRounds);
            if constexpr (Direction == BACKWARD) {
                const StopEvent* initialStopEvent = departures[i].initialStopEvent();
                const size_t stopEventId = initialStopEvent ? flagsData.permutateStopEvent(initialStopEvent - &(raptorData.stopEvents[0])) : -1;
                backwardSearch.setForceFinalTrip(departures[i].forceInitialTrip, stopEventId);
                if constexpr (UseAlternativeBackwardSearch) {
                    collectArrivalsBackwardAlternative();
                } else {
                    collectArrivalsBackward();
                }
            } else {
                collectArrivals();
            }
        }
        flagsMarker.reconstructShortcutsIntoCell();
    }

    inline int getWalkingTime(const StopId stop) noexcept {
        return forwardSearch.getWalkingTravelTime(stop);
    }

    inline void reset() noexcept {
        forwardSearch.reset();
        flagsMarker.clear();
    }

    inline const ShortcutFlags& getShortcutFlags() const noexcept {
        return flagsMarker.getShortcutFlags();
    }

private:
    inline void clear() noexcept {
        flagsMarker.clear();
    }

    inline void collectDepartures() noexcept {
        departures.clear();

        for (const Vertex vertex : forwardSearch.getStopsReachedByWalking()) {
            const StopId stop(vertex);
            for (const RouteSegment routeSegment : raptorData.routesContainingStop(stop)) {
                const RouteId route = routeSegment.routeId;
                const size_t numberOfStops = raptorData.numberOfStopsInRoute(route);
                if (routeSegment.stopIndex + 1 == numberOfStops) continue;
                const StopEvent* stopEvents = raptorData.firstTripOfRoute(route);
                const StopId* stops = raptorData.stopArrayOfRoute(route);
                for (size_t stopEventIndex = routeSegment.stopIndex; stopEventIndex < raptorData.numberOfStopEventsInRoute(route); stopEventIndex += numberOfStops) {
                    const int departureTime = stopEvents[stopEventIndex].departureTime - getWalkingTime(stop);
                    if (departureTime < minDepartureTime || departureTime > maxDepartureTime) continue;
                    departures.emplace_back(departureTime, stop);
                    if (stop == sourceVertex && routeSegment.stopIndex > 0 && flagsData.getCoarsePartition().isNonBoundaryIntraCellEdge(stops[routeSegment.stopIndex - 1], stop)) {
                        departures.emplace_back(std::max(departureTime, stopEvents[stopEventIndex].arrivalTime), stop, true, routeSegment, stopEvents + stopEventIndex - routeSegment.stopIndex);
                    }
                }
            }
        }
        sort(departures);
        if constexpr (Debug) std::cout << "Number of departures: " << departures.size() << std::endl;
    }

    inline void collectArrivalsBackward() {
        std::vector<ArrivalAtVertex> arrivals;
        for (const Vertex target : forwardSearch.getReachedTargetVertices()) {
            AssertMsg(forwardSearch.reachable(target), "Unreachable vertex " << target << " was marked as reachable!");
            for (const ArrivalLabel arrival : forwardSearch.template getArrivals<true>(target)) {
                arrivals.emplace_back(ArrivalAtVertex{target, arrival.arrivalTime, arrival.numberOfTrips});
            }
        }
        std::sort(arrivals.begin(), arrivals.end());
        for (const ArrivalAtVertex& arrival : arrivals) {
            backwardSearch.run(arrival.vertex, -arrival.arrivalTime, arrival.numberOfTrips);
        }
    }

    inline void collectArrivalsBackwardAlternative() {
        if (forwardSearch.getReachedTargetVertices().empty()) return;
        backwardSearch.runInitialize(sourceVertex, forwardSearch.maxNumberOfTrips());
        for (const Vertex target : forwardSearch.getReachedTargetVertices()) {
            AssertMsg(forwardSearch.reachable(target), "Unreachable vertex " << target << " was marked as reachable!");
            for (const ArrivalLabel arrival : forwardSearch.template getArrivals<true>(target)) {
                backwardSearch.addSource(target, arrival.numberOfTrips, -arrival.arrivalTime);
            }
        }
        backwardSearch.run();
    }

    inline void collectArrivals() noexcept {
        for (const Vertex target : forwardSearch.getReachedTargetVertices()) {
            AssertMsg(forwardSearch.reachable(target), "Unreachable vertex " << target << " was marked as reachable!");
            const std::vector<size_t> targetCells = flagsData.getFinePartition().cellsOfVertex(target);
            for (const PreprocessingJourney journey : forwardSearch.template getJourneys<true>(target)) {
                flagsMarker.markJourney(journey, targetCells);
            }
        }
    }

private:
    RouteFlagsDataType& flagsData;
    const Data& raptorData;
    ForwardSearchType forwardSearch;
    FlagsMarkerType flagsMarker;
    BackwardSearchType backwardSearch;

    Vertex sourceVertex;
    std::vector<size_t> sourceCells;
    int minDepartureTime;
    int maxDepartureTime;
    std::vector<DepartureLabel> departures;
};

}
