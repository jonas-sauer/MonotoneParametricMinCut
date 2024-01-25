#pragma once

#include "IntraCellFlagsCounter.h"
#include "../Data.h"
#include "../../../Helpers/Assert.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/IO/Serialization.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/Vector/Vector.h"

#include <iostream>
#include <vector>

namespace RAPTOR::RouteFlags {

class TripFlags {

public:
    TripFlags(const size_t fromCells, const size_t toCells, const Data& raptorData) :
        raptorData(raptorData),
        tripOfStopEvent(raptorData.mapStopEventsToTripId()),
        tripFlags(fromCells, std::vector<std::vector<bool>>(toCells, std::vector<bool>(raptorData.numberOfTrips(), false))),
        routeFlags(fromCells, std::vector<std::vector<bool>>(toCells, std::vector<bool>(raptorData.numberOfRoutes(), false))),
        stopFlags(fromCells, std::vector<std::vector<bool>>(toCells, std::vector<bool>(raptorData.numberOfStops(), false))) {
    }

    TripFlags(const Data& raptorData, const std::string& filename) :
        raptorData(raptorData),
        tripOfStopEvent(raptorData.mapStopEventsToTripId()) {
        deserialize(filename);
    }

    inline bool isNecessaryStopEvent(const bool, const size_t stopEvent, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryStopEvent(stopEvent, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryStopEvent(const size_t stopEvent, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < tripFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < tripFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stopEvent < tripOfStopEvent.size(), "Invalid stop event index: " << stopEvent);
        return isNecessaryTrip(tripOfStopEvent[stopEvent], fromCell, toCell);
    }

    inline bool isNecessaryRouteSegment(const bool, const RouteId route, const StopIndex stopIndex, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryRouteSegment(route, stopIndex, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryRouteSegment(const bool, const RouteId route, const StopIndex, const size_t fromCell, const size_t toCell) const noexcept {
        return isNecessaryRoute(route, fromCell, toCell);
    }

    inline bool isNecessaryRouteSegment(const RouteId route, const StopIndex, const size_t fromCell, const size_t toCell) const noexcept {
        return isNecessaryRoute(route, fromCell, toCell);
    }

    inline bool isNecessaryTrip(const TripId trip, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryTrip(trip, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryTrip(const TripId trip, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < tripFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < tripFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(trip < tripFlags[fromCell][toCell].size(), "Invalid trip index: " << trip);
        return tripFlags[fromCell][toCell][trip];
    }

    inline bool isNecessaryRoute(const RouteId route, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryRoute(route, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryRoute(const RouteId route, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < routeFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < routeFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(route < routeFlags[fromCell][toCell].size(), "Invalid route index: " << route);
        return routeFlags[fromCell][toCell][route];
    }

    inline bool isNecessaryStop(const bool, const StopId stop, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryStop(stop, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryStop(const StopId stop, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < stopFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < stopFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stop < stopFlags[fromCell][toCell].size(), "Invalid stop index: " << stop);
        return stopFlags[fromCell][toCell][stop];
    }

    inline void markStopEvent(const bool, const size_t stopEvent, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                markStopEvent(stopEvent, fromCell, toCell);
            }
        }
    }

    inline void markStopEvent(const size_t stopEvent, const size_t fromCell, const size_t toCell) noexcept {
        AssertMsg(fromCell < tripFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < tripFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stopEvent < tripOfStopEvent.size(), "Invalid stop event index: " << stopEvent);
        markTrip(tripOfStopEvent[stopEvent], fromCell, toCell);
    }

    inline void markTrip(const TripId trip, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                markTrip(trip, fromCell, toCell);
            }
        }
    }

    inline void markTrip(const TripId trip, const size_t fromCell, const size_t toCell) noexcept {
        AssertMsg(fromCell < tripFlags.size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < tripFlags[fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(trip < tripFlags[fromCell][toCell].size(), "Invalid trip index: " << trip);
        tripFlags[fromCell][toCell][trip] = true;
    }

    inline void markConnection(const bool, const RouteId route, const size_t, const std::vector<size_t>& fromCells) noexcept {
        for (const StopEvent* trip = raptorData.firstTripOfRoute(route); trip <= raptorData.lastTripOfRoute(route); trip += raptorData.numberOfStopsInRoute(route)) {
            for (const size_t fromCell : fromCells) {
                for (size_t toCell = 0; toCell < tripFlags[fromCell].size(); toCell++) {
                    markTrip(tripOfStopEvent[trip - &(raptorData.stopEvents[0])], fromCell, toCell);
                }
            }
        }
    }

    inline void computeExtraFlags() noexcept {
        for (size_t fromCell = 0; fromCell < tripFlags.size(); fromCell++) {
            for (size_t toCell = 0; toCell < tripFlags[fromCell].size(); toCell++) {
                for (const RouteId route : raptorData.routes()) {
                    const size_t tripSize = raptorData.numberOfStopsInRoute(route);
                    for (const StopEvent* trip = raptorData.firstTripOfRoute(route); trip <= raptorData.lastTripOfRoute(route); trip += tripSize) {
                        const size_t stopEventId = trip - &(raptorData.stopEvents[0]);
                        const TripId tripId = tripOfStopEvent[stopEventId];
                        if (tripFlags[fromCell][toCell][tripId]) {
                            routeFlags[fromCell][toCell][route] = true;
                            break;
                        }
                    }

                    if (routeFlags[fromCell][toCell][route]) {
                        for (const StopId stop : raptorData.stopsOfRoute(route)) {
                            stopFlags[fromCell][toCell][stop] = true;
                        }
                    }
                }
            }
        }
    }

    inline void serialize(const std::string& filename) const noexcept {
        IO::serialize(filename, tripFlags, routeFlags, stopFlags);
    }

    inline void deserialize(const std::string& filename) noexcept {
        IO::deserialize(filename, tripFlags, routeFlags, stopFlags);
    }

    inline long long byteSize() const noexcept {
        long long result = Vector::byteSize(tripOfStopEvent);
        result += Vector::byteSize(tripFlags);
        result += Vector::byteSize(routeFlags);
        result += Vector::byteSize(stopFlags);
        return result;
    }

    inline void printCellStatistics(const IntraCellFlagsCounter& intraCellFlagsCounter, const size_t fromCell, const size_t toCell) const noexcept {
        printStatistics([&](const auto& flags) {
            return Vector::count(flags[fromCell][toCell], true);
        }, intraCellFlagsCounter.intraCellTrips(fromCell), intraCellFlagsCounter.intraCellRoutes(fromCell), intraCellFlagsCounter.intraCellStops(DEPARTURE, fromCell));
    }


    inline void printAverageStatistics(const IntraCellFlagsCounter& intraCellFlagsCounter) const noexcept {
        printStatistics([&](const auto& flags) {
            return Vector::mean(flags, [&](const auto& v1) {
                return Vector::mean(v1, [&](const auto& v2) {
                    return Vector::count(v2, true);
                });
             });
        }, intraCellFlagsCounter.averageIntraCellTrips(), intraCellFlagsCounter.averageIntraCellRoutes(), intraCellFlagsCounter.averageIntraCellStops(DEPARTURE));
    }

private:
    template<typename COUNT_FLAGS, typename VALUE_TYPE>
    inline void printStatistics(const COUNT_FLAGS& countFlags, const VALUE_TYPE intraCellTrips, const VALUE_TYPE intraCellRoutes, const VALUE_TYPE intraCellStops) const noexcept {
        std::cout << "    Trips:" << std::endl;
        printCount(countFlags(tripFlags), intraCellTrips, raptorData.numberOfTrips());
        std::cout << "    Routes:" << std::endl;
        printCount(countFlags(routeFlags), intraCellRoutes, raptorData.numberOfRoutes());
        std::cout << "    Stops:" << std::endl;
        printCount(countFlags(stopFlags), intraCellStops, raptorData.numberOfStops());
    }

    template<typename VALUE_TYPE>
    inline static void printCount(const VALUE_TYPE flagged, const VALUE_TYPE intraCell, const size_t total) noexcept {
        std::cout << "        Total: " << flagged << " (" << String::percent((double)flagged / total) << ")" << std::endl;
        std::cout << "        Intra-cell: " << intraCell << " (" << String::percent((double)intraCell / total) << ")" << std::endl;
        std::cout << "        Inter-cell: " << (flagged - intraCell) << " (" << String::percent((double)(flagged - intraCell) / total) << ")" << std::endl;
    }

private:
    const Data& raptorData;
    std::vector<TripId> tripOfStopEvent;

    std::vector<std::vector<std::vector<bool>>> tripFlags;
    std::vector<std::vector<std::vector<bool>>> routeFlags;
    std::vector<std::vector<std::vector<bool>>> stopFlags;
};

}
