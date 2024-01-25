#pragma once

#include <vector>

#include "../../../Dijkstra/Dijkstra.h"
#include "../../../../DataStructures/Container/ExternalKHeap.h"
#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "FlagsMarker.h"
#include "ForwardSearch.h"

//TODO: forceFinalTrip: earliestArrivals?
namespace RAPTOR::RouteFlags {

template<typename FORWARD_SEARCH, typename FLAGS_MARKER>
class AlternativeBackwardSearch {

public:
    using ForwardSearchType = FORWARD_SEARCH;
    using FlagsMarkerType = FLAGS_MARKER;
    using Type = AlternativeBackwardSearch<ForwardSearchType, FlagsMarkerType>;
    using RouteFlagsDataType = typename ForwardSearchType::RouteFlagsDataType;

public:
    struct EarliestArrival {
        EarliestArrival() :
            arrivalTime(never),
            timestamp(0) {
        }

        inline void reset(const int time) {
            arrivalTime = never;
            timestamp = time;
        }

        int arrivalTime;
        int timestamp;
    };
    using EarliestArrivals = std::vector<EarliestArrival>;
    using Round = std::vector<EarliestArrivals>;

    struct Departure {
        Departure(const bool isArrival, const RouteId route, const StopIndex stop, const size_t tripNumber, const int departureTime) :
            isArrival(isArrival),
            route(route),
            stop(stop),
            tripNumber(tripNumber),
            departureTime(departureTime) {
        }

        bool isArrival;
        RouteId route;
        StopIndex stop;
        size_t tripNumber;
        int departureTime;

        inline bool operator<(const Departure& other) const noexcept {
            return departureTime < other.departureTime;
        }
    };

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

public:
    AlternativeBackwardSearch(const RouteFlagsDataType& flagsData, ForwardSearchType& forwardSearch, FlagsMarkerType& flagsMarker) :
        raptorData(flagsData.getRaptorData(FORWARD)),
        shortcutGraph(flagsData.getShortcutGraph(FORWARD)),
        flagsData(flagsData),
        forwardSearch(forwardSearch),
        flagsMarker(flagsMarker),
        queue(raptorData.transferGraph.numVertices()),
        dijkstraLabels(raptorData.transferGraph.numVertices()),
        targetVertex(noVertex),
        maxTrips(-1),
        forceFinalTrip(false),
        finalStopEvent(-1),
        timestamp(0),
        dijkstraTimestamp(0),
        reachedStop(raptorData.numberOfStops(), false) {
        AssertMsg(raptorData.hasImplicitBufferTimes(), "Route-Flags requires implicit buffer times!");
        buildDepartures();
    }

    inline void setForceFinalTrip(const bool force, const size_t stopEvent) noexcept {
        forceFinalTrip = force;
        finalStopEvent = stopEvent;
    }

    inline void runInitialize(const Vertex target, const size_t numTrips) noexcept {
        clear();
        initialize(target, numTrips);
    }

    inline void addSource(const Vertex source, const size_t numTrips, const int departureTime) noexcept {
        if (raptorData.isStop(source)) {
            reachedStop[source] = true;
            for (const size_t cell : flagsData.getFinePartition().cellsOfVertex(source)) {
                getEarliestArrival(numTrips, StopId(source), cell).arrivalTime = departureTime;
            }
        }
        if (numTrips == 0) {
            if (raptorData.isStop(source)) {
                markShortcut(targetVertex, source, getTargetArrivalTime() - getLatestDepartureTime(source, 0), flagsData.getFinePartition().cellsOfVertex(source));
            }
        } else {
            dijkstra(source, numTrips, departureTime);
        }
    }

    inline void run() noexcept {
        for (const Departure& departure : departures) {
            if (departure.departureTime > getTargetArrivalTime()) break;
            TripIterator tripIterator = raptorData.getTripIterator(departure.route, departure.stop, departure.tripNumber);
            if (!reachedStop[tripIterator.stop()]) continue;
            if (departure.isArrival) {
                relaxTransfers(tripIterator.stop(), tripIterator.arrivalTime());
            } else {
                scanTrip(tripIterator);
            }
        }
    }

private:
    inline void buildDepartures() noexcept {
        for (const RouteId route : raptorData.routes()) {
            TripIterator tripIterator = raptorData.getTripIterator(route);
            do {
                while (tripIterator.hasFurtherStops()) {
                    departures.emplace_back(false, route, tripIterator.getStopIndex(), tripIterator.getCurrentTripNumber(), tripIterator.departureTime());
                    tripIterator.nextStop();
                    departures.emplace_back(true, route, tripIterator.getStopIndex(), tripIterator.getCurrentTripNumber(), tripIterator.arrivalTime());
                }
                tripIterator.setStopIndex(StopIndex(0));
            } while (tripIterator.decreaseTrip());
        }
        std::sort(departures.begin(), departures.end());
    }

    inline void clear() noexcept {
        targetVertex = noVertex;
        maxTrips = -1;
        timestamp++;
        Vector::fill(reachedStop, false);
    }

    inline void initialize(const Vertex target, const size_t numTrips) noexcept {
        targetVertex = target;
        maxTrips = numTrips;
        while (earliestArrivals.size() < numTrips) {
            earliestArrivals.emplace_back(Round(raptorData.numberOfStops(), EarliestArrivals(flagsData.getFinePartition().numberOfCells())));
        }
    }

    inline void dijkstra(const Vertex source, const size_t numTrips, const int departureTime) noexcept {
        dijkstraTimestamp++;
        AssertMsg(queue.empty(), "Queue still has " << queue.size() << " elements!");
        DijkstraLabel& sourceLabel = getDijkstraLabel(source);
        sourceLabel.distance = departureTime;
        queue.update(&sourceLabel);
        while(!queue.empty()) {
            DijkstraLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(dijkstraLabels[0]));
            if (raptorData.isStop(u)) {
                reachedStop[u] = true;
                for (const size_t cell : flagsData.getFinePartition().cellsOfVertex(source)) {
                    EarliestArrival& arrival = getEarliestArrival(numTrips, StopId(u), cell);
                    if (uLabel->distance < arrival.arrivalTime) {
                        arrival.arrivalTime = uLabel->distance;
                    }
                }
            }
            if (u == targetVertex) {
                queue.clear();
                break;
            }
            for (const Edge edge : raptorData.transferGraph.edgesFrom(u)) {
                const Vertex v = raptorData.transferGraph.get(ToVertex, edge);
                DijkstraLabel& vLabel = getDijkstraLabel(v);
                const int distance = uLabel->distance + raptorData.transferGraph.get(TravelTime, edge);
                if (vLabel.distance > distance && distance <= getLatestDepartureTime(v, numTrips)) {
                    AssertMsg(distance == getLatestDepartureTime(v, numTrips), "Found shorter path! Source: " << source << ", target: " << targetVertex << ", stop: " << v << ", Source departure time: " << departureTime << ", Target arrival time: " << getTargetArrivalTime() << ", trips: " << numTrips << ", Distance: " << distance << ", forward: " << getLatestDepartureTime(v, numTrips));
                    vLabel.distance = distance;
                    queue.update(&vLabel);
                }
            }
        }
    }

    inline void relaxTransfers(const StopId from, const int departureTime) noexcept {
        for (size_t numTrips = 0; numTrips < maxTrips; numTrips++) {
            for (size_t cell = 0; cell < earliestArrivals[numTrips][from].size(); cell++) {
                if (getEarliestArrival(numTrips, from, cell).arrivalTime != departureTime) continue;
                if (numTrips == 0) {
                    markShortcut(targetVertex, from, getTargetArrivalTime() - getLatestDepartureTime(from, 0), cell);
                } else {
                    for (const Edge edge : shortcutGraph.edgesFrom(from)) {
                        const Vertex to = shortcutGraph.get(ToVertex, edge);
                        const int travelTime = shortcutGraph.get(TravelTime, edge);
                        const int arrivalTime = departureTime + travelTime;
                        if (arrivalTime > getLatestDepartureTime(to, numTrips)) continue;
                        reachedStop[to] = true;
                        markShortcut(to, from, travelTime, cell);
                        EarliestArrival& arrival = getEarliestArrival(numTrips, StopId(to), cell);
                        if (arrivalTime < arrival.arrivalTime) {
                            arrival.arrivalTime = arrivalTime;
                        }
                    }
                }
                EarliestArrival& arrival = getEarliestArrival(numTrips, from, cell);
                if (departureTime < arrival.arrivalTime) {
                    arrival.arrivalTime = departureTime;
                }
            }
        }
    }

    inline void scanTrip(TripIterator& tripIterator) noexcept {
        const StopId fromStop = tripIterator.stop();
        const int departureTime = tripIterator.departureTime();
        const StopEvent* departureEvent = tripIterator.stopEvent();
        while (tripIterator.hasFurtherStops()) {
            tripIterator.nextStop();
            const StopId toStop = tripIterator.stop();
            const int arrivalTime = tripIterator.arrivalTime();
            for (size_t numTrips = 1; numTrips < maxTrips; numTrips++) {
                const bool isFinalStopEvent = forceFinalTrip && size_t(tripIterator.stopEvent() - &(raptorData.stopEvents[0])) == finalStopEvent;
                if (!isFinalStopEvent && arrivalTime > getLatestDepartureTime(toStop, numTrips - 1)) continue;
                for (size_t cell = 0; cell < earliestArrivals[numTrips][fromStop].size(); cell++) {
                    if (getEarliestArrival(numTrips, fromStop, cell).arrivalTime > departureTime) continue;
                    reachedStop[toStop] = true;
                    markStopEvent(DEPARTURE, departureEvent, cell);
                    if (!isFinalStopEvent) markStopEvent(ARRIVAL, tripIterator.stopEvent(), cell);
                    EarliestArrival& arrival = getEarliestArrival(numTrips - 1, toStop, cell);
                    if (arrivalTime < arrival.arrivalTime) {
                        arrival.arrivalTime = arrivalTime;
                    }
                }
            }
        }
    }

    inline void markShortcut(const Vertex from, const Vertex to, const int travelTime, const std::vector<size_t> cells) noexcept {
        flagsMarker.markShortcut(Shortcut{from, to, travelTime}, cells);
    }

    inline void markShortcut(const Vertex from, const Vertex to, const int travelTime, const size_t cell) noexcept {
        flagsMarker.markShortcut(Shortcut{from, to, travelTime}, cell);
    }

    inline void markStopEvent(const bool departureOrArrival, const StopEvent* stopEvent, const size_t cell) noexcept {
        flagsMarker.markStopEvent(departureOrArrival, stopEvent, cell);
    }

    inline DijkstraLabel& getDijkstraLabel(const Vertex vertex) noexcept {
        DijkstraLabel& label = dijkstraLabels[vertex];
        if (label.timestamp != dijkstraTimestamp) {
            label.reset(dijkstraTimestamp);
        }
        return label;
    }

    inline EarliestArrival& getEarliestArrival(const size_t round, const StopId stop, const size_t cell) noexcept {
        EarliestArrival& arrival = earliestArrivals[round][stop][cell];
        if (arrival.timestamp != timestamp) {
            arrival.reset(timestamp);
        }
        return arrival;
    }

    inline int getLatestDepartureTime(const Vertex vertex, const size_t numTrips) const noexcept {
        return -forwardSearch.getArrivalTime(vertex, numTrips);
    }

    inline int getTargetArrivalTime() const noexcept {
        return -forwardSearch.getSourceDepartureTime();
    }

private:
    const Data& raptorData;
    const TransferGraph& shortcutGraph;
    const RouteFlagsDataType& flagsData;
    ForwardSearchType& forwardSearch;
    FlagsMarkerType& flagsMarker;

    ExternalKHeap<2, DijkstraLabel> queue;
    std::vector<DijkstraLabel> dijkstraLabels;

    std::vector<Departure> departures;
    std::vector<Round> earliestArrivals;
    Vertex targetVertex;
    size_t maxTrips;
    bool forceFinalTrip;
    size_t finalStopEvent;

    int timestamp;
    int dijkstraTimestamp;

    std::vector<bool> reachedStop;
};
}
