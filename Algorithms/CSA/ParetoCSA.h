#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Journey.h"
#include "Profiler.h"

namespace CSA {

template<bool PATH_RETRIEVAL = true, typename PROFILER = NoProfiler, size_t MAX_TRIPS = 8>
class ParetoCSA {

public:
    constexpr static bool PathRetrieval = PATH_RETRIEVAL;
    using Profiler = PROFILER;
    constexpr static size_t MaxTrips = MAX_TRIPS;
    using Type = ParetoCSA<PathRetrieval, Profiler, MaxTrips>;
    using TripFlag = Meta::IF<PathRetrieval, ConnectionId, bool>;

private:
    struct ParentLabel {
        ParentLabel(const StopId parent = noStop, const bool reachedByTransfer = false, const TripId tripId = noTripId) :
            parent(parent),
            reachedByTransfer(reachedByTransfer),
            tripId(tripId) {
        }

        StopId parent;
        bool reachedByTransfer;
        union {
            TripId tripId;
            Edge transferId;
        };
    };

public:
    ParetoCSA(const Data& data, const Profiler& profilerTemplate = Profiler()) :
        data(data),
        sourceStop(noStop),
        targetStop(noStop),
        tripReached(data.numberOfTrips() * MaxTrips, TripFlag()),
        arrivalTime(data.numberOfStops() * MaxTrips, never),
        parentLabel(PathRetrieval ? data.numberOfStops() * MaxTrips : 0),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline void run(const StopId source, const int departureTime, const StopId target = noStop) noexcept {
        profiler.start();

        profiler.startPhase();
        Assert(data.isStop(source), "Source stop " << source << " is not a valid stop!");
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        sourceStop = source;
        targetStop = target;
        arrivalTime[sourceStop * MaxTrips] = departureTime;
        relaxEdges(sourceStop, departureTime, 0);
        const ConnectionId firstConnection = firstReachableConnection(departureTime);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections(firstConnection, ConnectionId(data.connections.size()));
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.done();
    }

    inline bool reachable(const StopId stop) const noexcept {
        return getEarliestArrivalTime(stop) < never;
    }

    inline int getArrivalTime(const StopId stop, const size_t numTrips) const noexcept {
        return arrivalTime[stop * MaxTrips + numTrips];
    }

    inline int getEarliestArrivalTime(const StopId stop) const noexcept {
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
    inline std::vector<Journey> getJourneys(const StopId stop) const noexcept {
        std::vector<Journey> journeys;
        for (size_t i = 0; i < MaxTrips; i++) {
            getJourney(journeys, i, stop);
        }
        return journeys;
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline Journey getEarliestJourney(const StopId stop) const noexcept {
        std::vector<Journey> journeys = getJourneys(stop);
        return journeys.empty() ? Journey() : journeys.back();
    }

    inline std::vector<int> getArrivalTimes() const noexcept {
        return getArrivalTimes(targetStop);
    }

    inline std::vector<int> getArrivalTimes(const StopId stop) const noexcept {
        std::vector<int> arrivalTimes;
        for (size_t i = 0; i < MaxTrips; i++) {
            arrivalTimes.emplace_back(std::min(arrivalTime[stop * MaxTrips + i], (arrivalTimes.empty()) ? (never) : (arrivalTimes.back())));
        }
        return arrivalTimes;
    }

    inline std::vector<Vertex> getPath(const StopId stop) const noexcept {
        return journeyToPath(getJourneys(stop).back());
    }

    inline std::vector<std::string> getRouteDescription(const StopId stop) const noexcept {
        return data.journeyToText(getJourneys(stop).back());
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() {
        sourceStop = noStop;
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
        relaxEdges(stop, time, numTrips);
    }

    inline void relaxEdges(const StopId stop, const int time, const size_t numTrips) noexcept {
        for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
            profiler.countMetric(METRIC_EDGES);
            const StopId toStop = StopId(data.transferGraph.get(ToVertex, edge));
            const int newArrivalTime = time + data.transferGraph.get(TravelTime, edge);
            arrivalByTransfer(toStop, newArrivalTime, numTrips, stop, edge);
        }
    }

    inline void arrivalByTransfer(const StopId stop, const int time, const size_t numTrips, const StopId parent, const Edge edge) noexcept {
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
    inline void getJourney(std::vector<Journey>& journeys, size_t numTrips, StopId stop) const noexcept {
        if (arrivalTime[stop * MaxTrips + numTrips] >= (journeys.empty() ? never : journeys.back().back().arrivalTime)) return;
        Journey journey;
        while (stop != sourceStop) {
            const ParentLabel& label = parentLabel[stop * MaxTrips + numTrips];
            const int time = arrivalTime[stop * MaxTrips + numTrips];
            if (label.reachedByTransfer) {
                const int travelTime = data.transferGraph.get(TravelTime, label.transferId);
                journey.emplace_back(label.parent, stop, time - travelTime, time, label.transferId);
            } else {
                journey.emplace_back(label.parent, stop, data.connections[tripReached[label.tripId * MaxTrips + numTrips]].departureTime, time, label.tripId);
                numTrips--;
            }
            stop = label.parent;
        }
        journeys.emplace_back(Vector::reverse(journey));
    }

private:
    const Data& data;

    StopId sourceStop;
    StopId targetStop;

    std::vector<TripFlag> tripReached;
    std::vector<int> arrivalTime;
    std::vector<ParentLabel> parentLabel;

    Profiler profiler;
};
}
