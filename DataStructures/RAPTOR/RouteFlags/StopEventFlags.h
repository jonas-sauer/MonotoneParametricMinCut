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

class StopEventFlags {

public:
    StopEventFlags(const size_t fromCells, const size_t toCells, const Data& raptorData) :
        raptorData(raptorData),
        stopEventFlags({buildFlags(fromCells, toCells, raptorData.numberOfStopEvents()), buildFlags(fromCells, toCells, raptorData.numberOfStopEvents())}),
        routeSegmentFlags({buildFlags(fromCells, toCells, raptorData.numberOfRouteSegments()), buildFlags(fromCells, toCells, raptorData.numberOfRouteSegments())}),
        stopFlags({buildFlags(fromCells, toCells, raptorData.numberOfStops()), buildFlags(fromCells, toCells, raptorData.numberOfStops())}) {
    }

    StopEventFlags(const Data& raptorData, const std::string& filename) :
        raptorData(raptorData) {
        deserialize(filename);
    }

    inline bool isNecessaryStopEvent(const bool departureOrArrival, const size_t stopEvent, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryStopEvent(departureOrArrival, stopEvent, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryStopEvent(const bool departureOrArrival, const size_t stopEvent, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < stopEventFlags[departureOrArrival].size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < stopEventFlags[departureOrArrival][fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stopEvent < stopEventFlags[departureOrArrival][fromCell][toCell].size(), "Invalid stop event index: " << stopEvent);
        return stopEventFlags[departureOrArrival][fromCell][toCell][stopEvent];
    }

    inline bool isNecessaryRouteSegment(const bool departureOrArrival, const RouteId route, const StopIndex stopIndex, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryRouteSegment(departureOrArrival, route, stopIndex, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryRouteSegment(const bool departureOrArrival, const RouteId route, const StopIndex stopIndex, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < routeSegmentFlags[departureOrArrival].size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < routeSegmentFlags[departureOrArrival][fromCell].size(), "Invalid cell index: " << toCell);
        const size_t routeSegment = raptorData.getRouteSegmentNum(route, stopIndex);
        AssertMsg(routeSegment < routeSegmentFlags[departureOrArrival][fromCell][toCell].size(), "Invalid route segment index: " << routeSegment);
        return routeSegmentFlags[departureOrArrival][fromCell][toCell][routeSegment];
    }

    inline bool isNecessaryStop(const bool departureOrArrival, const StopId stop, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) const noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                if (!isNecessaryStop(departureOrArrival, stop, fromCell, toCell)) return false;
            }
        }
        return true;
    }

    inline bool isNecessaryStop(const bool departureOrArrival, const StopId stop, const size_t fromCell, const size_t toCell) const noexcept {
        AssertMsg(fromCell < stopFlags[departureOrArrival].size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < stopFlags[departureOrArrival][fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stop < stopFlags[departureOrArrival][fromCell][toCell].size(), "Invalid stop index: " << stop);
        return stopFlags[departureOrArrival][fromCell][toCell][stop];
    }

    inline void markStopEvent(const bool departureOrArrival, const size_t stopEvent, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) noexcept {
        for (const size_t fromCell : fromCells) {
            for (const size_t toCell : toCells) {
                markStopEvent(departureOrArrival, stopEvent, fromCell, toCell);
            }
        }
    }

    inline void markStopEvent(const bool departureOrArrival, const size_t stopEvent, const size_t fromCell, const size_t toCell) noexcept {
        AssertMsg(fromCell < stopEventFlags[departureOrArrival].size(), "Invalid cell index: " << fromCell);
        AssertMsg(toCell < stopEventFlags[departureOrArrival][fromCell].size(), "Invalid cell index: " << toCell);
        AssertMsg(stopEvent < stopEventFlags[departureOrArrival][fromCell][toCell].size(), "Invalid stop event index: " << stopEvent);
        stopEventFlags[departureOrArrival][fromCell][toCell][stopEvent] = true;
    }

    inline void markConnection(const bool direction, const RouteId route, const StopIndex stopId, const std::vector<size_t>& fromCells) noexcept {
        for (const StopEvent* trip = raptorData.firstTripOfRoute(route); trip <= raptorData.lastTripOfRoute(route); trip += raptorData.numberOfStopsInRoute(route)) {
            for (const size_t fromCell : fromCells) {
                for (size_t toCell = 0; toCell < stopEventFlags[DEPARTURE][fromCell].size(); toCell++) {
                    markStopEvent(direction, (size_t)(trip + stopId - &(raptorData.stopEvents[0])), fromCell, toCell);
                    markStopEvent(!direction, (size_t)(trip + stopId + 1 - &(raptorData.stopEvents[0])), fromCell, toCell);
                }
            }
        }
    }

    inline void computeExtraFlags() noexcept {
        computeExtraFlags(DEPARTURE);
        computeExtraFlags(ARRIVAL);
    }

    inline void computeExtraFlags(const bool departureOrArrival) noexcept {
        for (size_t fromCell = 0; fromCell < stopEventFlags[departureOrArrival].size(); fromCell++) {
            for (size_t toCell = 0; toCell < stopEventFlags[departureOrArrival][fromCell].size(); toCell++) {
                for (const RouteSegment& segment : raptorData.routeSegments) {
                    const size_t segmentNum = raptorData.getRouteSegmentNum(segment.routeId, segment.stopIndex);
                    const StopId stop = raptorData.stopsOfRoute(segment.routeId)[segment.stopIndex];
                    for (const StopEvent* trip = raptorData.firstTripOfRoute(segment.routeId); trip <= raptorData.lastTripOfRoute(segment.routeId); trip += raptorData.numberOfStopsInRoute(segment.routeId)) {
                        if (stopEventFlags[departureOrArrival][fromCell][toCell][trip + segment.stopIndex - &(raptorData.stopEvents[0])]) {
                            routeSegmentFlags[departureOrArrival][fromCell][toCell][segmentNum] = true;
                            stopFlags[departureOrArrival][fromCell][toCell][stop] = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    inline void serialize(const std::string& filename) const noexcept {
        IO::serialize(filename, stopEventFlags[DEPARTURE], stopEventFlags[ARRIVAL], routeSegmentFlags[DEPARTURE], routeSegmentFlags[ARRIVAL], stopFlags[DEPARTURE], stopFlags[ARRIVAL]);
    }

    inline void deserialize(const std::string& filename) noexcept {
        IO::deserialize(filename, stopEventFlags[DEPARTURE], stopEventFlags[ARRIVAL], routeSegmentFlags[DEPARTURE], routeSegmentFlags[ARRIVAL], stopFlags[DEPARTURE], stopFlags[ARRIVAL]);
    }

    inline long long byteSize() const noexcept {
        long long result = Vector::byteSize(stopEventFlags[DEPARTURE]);
        result += Vector::byteSize(stopEventFlags[ARRIVAL]);
        result += Vector::byteSize(routeSegmentFlags[DEPARTURE]);
        result += Vector::byteSize(routeSegmentFlags[ARRIVAL]);
        result += Vector::byteSize(stopFlags[DEPARTURE]);
        result += Vector::byteSize(stopFlags[ARRIVAL]);
        return result;
    }

    inline void printCellStatistics(const IntraCellFlagsCounter& intraCellFlagsCounter, const size_t fromCell, const size_t toCell) const noexcept {
        printStatistics([&](const bool departureOrArrival, const auto& flags) {
            return Vector::count(flags[departureOrArrival][fromCell][toCell], true);
        }, intraCellFlagsCounter.intraCellStopEvents(fromCell), intraCellFlagsCounter.intraCellRouteSegments(fromCell), [&](const bool departureOrArrival) {
            return intraCellFlagsCounter.intraCellStops(departureOrArrival, fromCell);
        });
    }


    inline void printAverageStatistics(const IntraCellFlagsCounter& intraCellFlagsCounter) const noexcept {
        printStatistics([&](const bool departureOrArrival, const auto& flags) {
            return Vector::mean(flags[departureOrArrival], [&](const auto& v1) {
                return Vector::mean(v1, [&](const auto& v2) {
                    return Vector::count(v2, true);
                });
             });
        }, intraCellFlagsCounter.averageIntraCellStopEvents(), intraCellFlagsCounter.averageIntraCellRouteSegments(), [&](const bool departureOrArrival) {
            return intraCellFlagsCounter.averageIntraCellStops(departureOrArrival);
        });
    }

private:
    inline std::vector<std::vector<std::vector<bool>>> buildFlags(const size_t fromCells, const size_t toCells, const size_t numFlags) const noexcept {
        return std::vector<std::vector<std::vector<bool>>>(fromCells, std::vector<std::vector<bool>>(toCells, std::vector<bool>(numFlags, false)));
    }

    template<typename COUNT_FLAGS, typename VALUE_TYPE, typename INTRA_CELL_STOPS>
    inline void printStatistics(const COUNT_FLAGS& countFlags, const VALUE_TYPE intraCellStopEvents, const VALUE_TYPE intraCellRouteSegments, const INTRA_CELL_STOPS& intraCellStops) const noexcept {
        std::cout << "    Departure:" << std::endl;
        printStatistics(DEPARTURE, countFlags, intraCellStopEvents, intraCellRouteSegments, intraCellStops(DEPARTURE));
        std::cout << "    Arrival:" << std::endl;
        printStatistics(ARRIVAL, countFlags, intraCellStopEvents, intraCellRouteSegments, intraCellStops(ARRIVAL));
    }

    template<typename COUNT_FLAGS, typename VALUE_TYPE>
    inline void printStatistics(const bool departureOrArrival, const COUNT_FLAGS& countFlags, const VALUE_TYPE intraCellStopEvents, const VALUE_TYPE intraCellRouteSegments, const VALUE_TYPE intraCellStops) const noexcept {
        std::cout << "        Stop events:" << std::endl;
        printCount(countFlags(departureOrArrival, stopEventFlags), intraCellStopEvents, raptorData.numberOfStopEvents());
        std::cout << "        Route segments:" << std::endl;
        printCount(countFlags(departureOrArrival, routeSegmentFlags), intraCellRouteSegments, raptorData.numberOfRouteSegments());
        std::cout << "        Stops:" << std::endl;
        printCount(countFlags(departureOrArrival, stopFlags), intraCellStops, raptorData.numberOfStops());
    }

    template<typename VALUE_TYPE>
    inline static void printCount(const VALUE_TYPE flagged, const VALUE_TYPE intraCell, const size_t total) noexcept {
        std::cout << "            Total: " << flagged << " (" << String::percent((double)flagged / total) << ")" << std::endl;
        std::cout << "            Intra-cell: " << intraCell << " (" << String::percent((double)intraCell / total) << ")" << std::endl;
        std::cout << "            Inter-cell: " << (flagged - intraCell) << " (" << String::percent((double)(flagged - intraCell) / total) << ")" << std::endl;
    }

private:
    const Data& raptorData;

    std::vector<std::vector<std::vector<bool>>> stopEventFlags[2];
    std::vector<std::vector<std::vector<bool>>> routeSegmentFlags[2];
    std::vector<std::vector<std::vector<bool>>> stopFlags[2];
};

}
