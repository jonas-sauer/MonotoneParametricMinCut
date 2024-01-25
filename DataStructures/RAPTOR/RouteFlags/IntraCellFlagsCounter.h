#pragma once

#include "../Data.h"
#include "VertexSeparatorPartition.h"
#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/Types.h"

#include <vector>

namespace RAPTOR::RouteFlags {

class IntraCellFlagsCounter {
public:
    IntraCellFlagsCounter(const VertexSeparatorPartition& partition, const Data& raptorData) :
        raptorData(raptorData),
        intraCellSegmentsPerRoute(partition.numberOfCells(), std::vector<size_t>(raptorData.numberOfRoutes(), 0)),
        isIntraCellStop{buildIsIntraCellStop(partition), buildIsIntraCellStop(partition)}{
        for (const RouteId route : raptorData.routes()) {
            const StopId* stops = raptorData.stopArrayOfRoute(route);
            const size_t tripSize = raptorData.numberOfStopsInRoute(route);
            for (size_t s = 1; s < tripSize; s++) {
                const StopId from = stops[s-1];
                const StopId to = stops[s];
                const size_t cell = partition.getNonBoundaryIntraCellEdgeCellId(from, to);
                if (cell == size_t(-1)) continue;
                AssertMsg(cell < partition.numberOfCells(), "Invalid cell: " << cell);
                intraCellSegmentsPerRoute[cell][route]++;
                isIntraCellStop[DEPARTURE][cell][from] = true;
                isIntraCellStop[ARRIVAL][cell][from] = true;
            }
        }
    }

    inline size_t intraCellStopEvents(const size_t cell) const noexcept {
        return intraCellStopEvents(intraCellSegmentsPerRoute[cell]);
    }

    inline double averageIntraCellStopEvents() const noexcept {
        return Vector::mean(intraCellSegmentsPerRoute, [&](const auto& segments) {
            return intraCellStopEvents(segments);
        });
    }

    inline size_t intraCellRouteSegments(const size_t cell) const noexcept {
        return intraCellRouteSegments(intraCellSegmentsPerRoute[cell]);
    }

    inline double averageIntraCellRouteSegments() const noexcept {
        return Vector::mean(intraCellSegmentsPerRoute, [&](const auto& segments) {
            return intraCellRouteSegments(segments);
        });
    }

    inline size_t intraCellTrips(const size_t cell) const noexcept {
        return intraCellTrips(intraCellSegmentsPerRoute[cell]);
    }

    inline double averageIntraCellTrips() const noexcept {
        return Vector::mean(intraCellSegmentsPerRoute, [&](const auto& segments) {
            return intraCellTrips(segments);
        });
    }

    inline size_t intraCellRoutes(const size_t cell) const noexcept {
        return intraCellRoutes(intraCellSegmentsPerRoute[cell]);
    }

    inline double averageIntraCellRoutes() const noexcept {
        return Vector::mean(intraCellSegmentsPerRoute, [&](const auto& segments) {
            return intraCellRoutes(segments);
        });
    }

    inline size_t intraCellStops(const size_t cell, const bool departureOrArrival) const noexcept {
        return intraCellStops(isIntraCellStop[departureOrArrival][cell]);
    }

    inline double averageIntraCellStops(const bool departureOrArrival) const noexcept {
        return Vector::mean(isIntraCellStop[departureOrArrival], [&](const auto& stops) {
            return intraCellStops(stops);
        });
    }

private:
    inline std::vector<std::vector<bool>> buildIsIntraCellStop(const VertexSeparatorPartition& partition) const noexcept {
        return std::vector<std::vector<bool>>(partition.numberOfCells(), std::vector<bool>(raptorData.numberOfStops(), false));
    }

    inline size_t intraCellStopEvents(const std::vector<size_t>& cellSegments) const noexcept {
        size_t intraCellStopEvents = 0;
        for (const RouteId route : raptorData.routes()) {
            intraCellStopEvents += raptorData.numberOfTripsInRoute(route) * cellSegments[route];
        }
        return intraCellStopEvents;
    }

    inline size_t intraCellRouteSegments(const std::vector<size_t>& cellSegments) const noexcept {
        return Vector::sum(cellSegments);
    }

    inline size_t intraCellTrips(const std::vector<size_t>& cellSegments) const noexcept {
        size_t intraCellTrips = 0;
        for (const RouteId route : raptorData.routes()) {
            if (cellSegments[route] > 0) {
                intraCellTrips += raptorData.numberOfTripsInRoute(route);
            }
        }
        return intraCellTrips;
    }

    inline size_t intraCellRoutes(const std::vector<size_t>& cellSegments) const noexcept {
        return Vector::count(cellSegments, [&](const size_t i) {
           return i > 0;
        });
    }

    inline size_t intraCellStops(const std::vector<bool>& cellStops) const noexcept {
        return Vector::count(cellStops, true);
    }

    const Data& raptorData;
    std::vector<std::vector<size_t>> intraCellSegmentsPerRoute;
    std::vector<std::vector<bool>> isIntraCellStop[2];
};
}
