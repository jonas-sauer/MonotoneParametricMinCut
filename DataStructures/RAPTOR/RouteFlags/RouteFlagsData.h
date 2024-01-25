#pragma once

#include "Shortcuts.h"
#include "ShortcutFlags.h"
#include "StopEventFlags.h"
#include "TripFlags.h"
#include "RouteFlags.h"
#include "VertexSeparatorPartition.h"
#include "EdgeSeparatorPartition.h"
#include "IntraCellFlagsCounter.h"

#include "../../Graph/Graph.h"
#include "../../Partition/NestedDissection.h"
#include "../../Partition/VertexPartition.h"
#include "../Data.h"
#include "../../../Helpers/IO/Serialization.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/Vector/Permutation.h"
#include "../../../Helpers/Vector/Vector.h"

#include <iostream>
#include <string>
#include <vector>

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE>
class RouteFlagsData {
public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    using Type = RouteFlagsData<Flags, FinePartitionType>;

    template<typename T = FinePartitionType, typename = std::enable_if_t<Meta::Equals<T, FinePartitionType>() && Meta::Equals<T, EdgeSeparatorPartition>()>>
    RouteFlagsData(const Data& forwardRaptorData, const Data& backwardRaptorData, const Permutation& stopEventPermutation, const TransferGraph& forwardShortcuts, const TransferGraph& backwardShortcuts, const NestedDissection& coarsePartition, const size_t level, const VertexPartition& finePartition) :
        raptorData{forwardRaptorData, backwardRaptorData},
        stopEventPermutation(stopEventPermutation),
        coarsePartition(coarsePartition, forwardRaptorData, level),
        finePartition(finePartition, forwardRaptorData),
        shortcuts{createShortcuts(forwardShortcuts), createShortcuts(backwardShortcuts)},
        flags{createFlags(FORWARD), createFlags(BACKWARD)} {
        AssertMsg(level < coarsePartition.numberOfLevels(), "Invalid level " << level << "! (#levels: " << coarsePartition.numberOfLevels() << ")");
    }

        template<typename T = FinePartitionType, typename = std::enable_if_t<Meta::Equals<T, FinePartitionType>() && Meta::Equals<T, VertexSeparatorPartition>()>>
    RouteFlagsData(const Data& forwardRaptorData, const Data& backwardRaptorData, const Permutation& stopEventPermutation, const TransferGraph& forwardShortcuts, const TransferGraph& backwardShortcuts, const NestedDissection& coarsePartition, const size_t coarseLevel, const NestedDissection& finePartition, const size_t fineLevel) :
        raptorData{forwardRaptorData, backwardRaptorData},
        stopEventPermutation(stopEventPermutation),
        coarsePartition(coarsePartition, forwardRaptorData, coarseLevel),
        finePartition(finePartition, forwardRaptorData, fineLevel),
        shortcuts{createShortcuts(forwardShortcuts), createShortcuts(backwardShortcuts)},
        flags{createFlags(FORWARD), createFlags(BACKWARD)} {
        AssertMsg(coarseLevel < coarsePartition.numberOfLevels(), "Invalid level " << coarseLevel << "! (#levels: " << coarsePartition.numberOfLevels() << ")");
        AssertMsg(fineLevel < finePartition.numberOfLevels(), "Invalid level " << fineLevel << "! (#levels: " << finePartition.numberOfLevels() << ")");
    }

    RouteFlagsData(const std::string filename) :
        raptorData{Data::FromBinary(filename + ".raptor.forward"), Data::FromBinary(filename + ".raptor.backward")},
        coarsePartition(filename + ".coarse"),
        finePartition(filename + ".fine"),
        shortcuts{Shortcuts(filename + ".shortcuts.forward"), Shortcuts(filename + ".shortcuts.backward")},
        flags{Flags(raptorData[FORWARD], filename + ".flags.forward"), Flags(raptorData[BACKWARD], filename + ".flags.backward")} {
        deserialize(filename);
    }

    inline const Data& getRaptorData(const bool direction) const noexcept {
        return raptorData[direction];
    }

    inline const VertexSeparatorPartition& getCoarsePartition() const noexcept {
        return coarsePartition;
    }

    inline const FinePartitionType& getFinePartition() const noexcept {
        return finePartition;
    }

    inline size_t highestBoundaryVertexId() const noexcept {
        return std::max(coarsePartition.highestBoundaryVertexId(), finePartition.highestBoundaryVertexId());
    }

    inline bool isStopOrFineBoundaryVertex(const Vertex vertex) const noexcept {
        return raptorData[FORWARD].isStop(vertex) || finePartition.isBoundaryVertex(vertex);
    }

    inline const TransferGraph& getShortcutGraph(const bool direction) const noexcept {
        return shortcuts[direction].allShortcuts;
    }

    inline const Shortcuts& getShortcuts(const bool direction) const noexcept {
        return shortcuts[direction];
    }

    inline bool hasShortcut(const bool direction, const Vertex from, const Vertex to) const noexcept {
        return shortcuts[direction].allShortcuts.hasEdge(from, to);
    }

    inline size_t permutateStopEvent(const size_t originalStopEvent) const noexcept {
        return stopEventPermutation[originalStopEvent];
    }

    inline bool isNecessaryStopEvent(const bool direction, const bool departureOrArrival, const size_t stopEvent, const std::vector<size_t>& sourceCoarseCells, const std::vector<size_t>& targetCoarseCells, const std::vector<size_t>& sourceFineCells, const std::vector<size_t>& targetFineCells) const noexcept {
        return flags[direction].isNecessaryStopEvent(departureOrArrival, stopEvent, targetCoarseCells, sourceFineCells) && flags[!direction].isNecessaryStopEvent(departureOrArrival, stopEventPermutation[stopEvent], sourceCoarseCells, targetFineCells);
    }

    inline bool isNecessaryStopEvent(const bool direction, const bool departureOrArrival, const size_t stopEvent, const size_t fromCell, const size_t toCell) const noexcept {
        return flags[direction].isNecessaryStopEvent(departureOrArrival, stopEvent, fromCell, toCell);
    }

    inline bool isNecessaryRouteSegment(const bool direction, const bool departureOrArrival, const RouteId route, const StopIndex stopIndex, const std::vector<size_t>& sourceCoarseCells, const std::vector<size_t>& targetCoarseCells, const std::vector<size_t>& sourceFineCells, const std::vector<size_t>& targetFineCells) const noexcept {
        return flags[direction].isNecessaryRouteSegment(departureOrArrival, route, stopIndex, targetCoarseCells, sourceFineCells) && flags[!direction].isNecessaryRouteSegment(departureOrArrival, route, StopIndex(raptorData[direction].numberOfStopsInRoute(route) - stopIndex - 1), sourceCoarseCells, targetFineCells);
    }

    inline bool isNecessaryRouteSegment(const bool direction, const bool departureOrArrival, const RouteId route, const StopIndex stopIndex, const size_t fromCell, const size_t toCell) const noexcept {
        return flags[direction].isNecessaryRouteSegment(departureOrArrival, route, stopIndex, fromCell, toCell);
    }

    inline bool isNecessaryStop(const bool direction, const bool departureOrArrival, const StopId stop, const std::vector<size_t>& sourceCoarseCells, const std::vector<size_t> targetCoarseCells, const std::vector<size_t>& sourceFineCells, const std::vector<size_t> targetFineCells) const noexcept {
        return flags[direction].isNecessaryStop(departureOrArrival, stop, targetCoarseCells, sourceFineCells) && flags[!direction].isNecessaryStop(departureOrArrival, stop, sourceCoarseCells, targetFineCells);
    }

    inline void markStopEvent(const bool direction, const bool departureOrArrival, const size_t stopEvent, const bool permutate, const std::vector<size_t>& fromCells, const std::vector<size_t>& toCells) noexcept {
        AssertMsg(stopEvent < stopEventPermutation.size(), "Invalid stop event index!");
        flags[direction].markStopEvent(departureOrArrival, permutate ? stopEventPermutation[stopEvent] : stopEvent, fromCells, toCells);
    }

    inline void markIntraCellStopEvents() noexcept {
        markIntraCellStopEvents(FORWARD);
        markIntraCellStopEvents(BACKWARD);
    }

    inline void markIntraCellStopEvents(const bool direction) noexcept {
        for (const RouteId route : raptorData[direction].routes()) {
            const StopId* stops = raptorData[direction].stopArrayOfRoute(route);
            const size_t tripSize = raptorData[direction].numberOfStopsInRoute(route);
            for (StopIndex s = StopIndex(0); s < tripSize - 1; s++) {
                const StopId from = stops[s];
                const StopId to = stops[s+1];
                if (!coarsePartition.areInCommonCell(from, to) || coarsePartition.isBoundaryVertex(to)) continue;
                flags[direction].markConnection(direction, route, s, coarsePartition.incidentCells[from]);
            }
        }
    }

    inline void buildShortcuts(const ShortcutFlags& forwardFlags, const ShortcutFlags& backwardFlags) noexcept {
        const auto commonCells = [&](const Vertex u, const Vertex v) {
            return coarsePartition.getCommonCells(u, v);
        };
        shortcuts[FORWARD].buildFromFlags(commonCells, forwardFlags.shortcutFlags, false);
        shortcuts[BACKWARD].buildFromFlags(commonCells, backwardFlags.shortcutFlags, true);
    }

    inline void computeExtraFlags() noexcept {
        flags[FORWARD].computeExtraFlags();
        flags[BACKWARD].computeExtraFlags();
    }

    inline void sortShortcuts() noexcept {
        shortcuts[FORWARD].sort();
        shortcuts[BACKWARD].sort();
    }

    void deserialize(const std::string filename, const bool verbose = false) noexcept {
        if (verbose) std::cout << "Reading Route-Flags topology data from " << filename << "... " << std::flush;
        IO::deserialize(filename, stopEventPermutation);
        if (verbose) std::cout << "done." << std::endl;
    }

    void serialize(const std::string filename, const bool verbose = false) const noexcept {
        if (verbose) std::cout << "Writing Route-Flags topology data to " << filename << "... " << std::flush;
        IO::serialize(filename, stopEventPermutation);
        raptorData[FORWARD].serialize(filename + ".raptor.forward");
        raptorData[BACKWARD].serialize(filename + ".raptor.backward");
        coarsePartition.serialize(filename + ".coarse");
        finePartition.serialize(filename + ".fine");
        shortcuts[FORWARD].serialize(filename + ".shortcuts.forward");
        shortcuts[BACKWARD].serialize(filename + ".shortcuts.backward");
        flags[FORWARD].serialize(filename + ".flags.forward");
        flags[BACKWARD].serialize(filename + ".flags.backward");
        if (verbose) std::cout << "done." << std::endl;
    }

    inline long long byteSize() const noexcept {
        long long result = raptorData[FORWARD].byteSize();
        result += raptorData[BACKWARD].byteSize();
        result += Vector::byteSize(stopEventPermutation);
        result += coarsePartition.byteSize();
        result += finePartition.byteSize();
        result += shortcuts[FORWARD].byteSize();
        result += shortcuts[BACKWARD].byteSize();
        result += flags[FORWARD].byteSize();
        result += flags[BACKWARD].byteSize();
        return result;
    }

    inline void printCellStatistics(const size_t fromCell, const size_t toCell) const noexcept {
        printNetworkAndPartitionStatistics();
        printShortcutCellStatistics(fromCell, toCell);
        printFlagCellStatistics(fromCell, toCell);
    }

    inline void printAverageStatistics() const noexcept {
        printNetworkAndPartitionStatistics();
        printShortcutAverageStatistics();
        printFlagAverageStatistics();
    }

    inline void printNetworkAndPartitionStatistics() const noexcept {
        std::cout << std::endl;
        std::cout << "RAPTOR network:" << std::endl;
        std::cout << "    #Stops: " << raptorData[FORWARD].numberOfStops() << std::endl;
        std::cout << "    #Routes: " << raptorData[FORWARD].numberOfRoutes() << std::endl;
        std::cout << "    #Stop events: " << raptorData[FORWARD].numberOfStopEvents() << std::endl;
        std::cout << "    #Route segments: " << raptorData[FORWARD].numberOfRouteSegments() << std::endl;
        std::cout << "Footpath graph:" << std::endl;
        std::cout << "    #Vertices: " << raptorData[FORWARD].transferGraph.numVertices() << std::endl;
        std::cout << "    #Edges: " << raptorData[FORWARD].transferGraph.numEdges() << std::endl;
        std::cout << "    #Shortcuts: " << shortcuts[FORWARD].allShortcuts.numEdges() << std::endl;
        std::cout << "Coarse partition:" << std::endl;
        coarsePartition.printStatistics();
        std::cout << "Fine partition:" << std::endl;
        finePartition.printStatistics();
        std::cout << std::endl;
    }

    inline void printShortcutCellStatistics(const size_t fromCell, const size_t toCell) const noexcept {
        std::cout << "Shortcuts from " << fromCell << " to " << toCell << std::endl;
        std::cout << "Forward:" << std::endl;
        shortcuts[FORWARD].printCellStatistics(fromCell, toCell);
        std::cout << "Backward:" << std::endl;
        shortcuts[BACKWARD].printCellStatistics(fromCell, toCell);
        std::cout << std::endl;
    }

    inline void printShortcutAverageStatistics() const noexcept {
        std::cout << "Shortcuts" << std::endl;
        std::cout << "Forward:" << std::endl;
        shortcuts[FORWARD].printAverageStatistics();
        std::cout << "Backward:" << std::endl;
        shortcuts[BACKWARD].printAverageStatistics();
        std::cout << std::endl;
    }

    inline void printFlagCellStatistics(const size_t fromCell, const size_t toCell) const noexcept {
        std::cout << "Shortcuts from " << fromCell << " to " << toCell << std::endl;
        std::cout << "Forward:" << std::endl;
        IntraCellFlagsCounter forwardCounter(coarsePartition, raptorData[FORWARD]);
        flags[FORWARD].printCellStatistics(forwardCounter, fromCell, toCell);
        std::cout << "Backward:" << std::endl;
        IntraCellFlagsCounter backwardCounter(coarsePartition, raptorData[BACKWARD]);
        flags[BACKWARD].printCellStatistics(backwardCounter, fromCell, toCell);
        std::cout << std::endl;
    }

    inline void printFlagAverageStatistics() const noexcept {
        std::cout << "Flags" << std::endl;
        std::cout << "Forward:" << std::endl;
        IntraCellFlagsCounter forwardCounter(coarsePartition, raptorData[FORWARD]);
        flags[FORWARD].printAverageStatistics(forwardCounter);
        std::cout << "Backward:" << std::endl;
        IntraCellFlagsCounter backwardCounter(coarsePartition, raptorData[BACKWARD]);
        flags[BACKWARD].printAverageStatistics(backwardCounter);
        std::cout << std::endl;
    }

    inline void printVertexCells(const Vertex v, std::ostream& out = std::cout) const noexcept {
        out << "Cells of " << v << ": " << std::endl;
        out << "    From: ";
        Vector::printConcise(coarsePartition.cellsOfVertex(v), out);
        out << std::endl;
        out << "    To: ";
        Vector::printConcise(finePartition.cellsOfVertex(v), out);
        out << std::endl;
    }

    inline void printFlags(const Vertex source, const Vertex target, const RouteId route, const StopId stop) const noexcept {
        StopIndex stopIndex = noStopIndex;
        for (const RouteSegment& routeSegment : raptorData[FORWARD].routesContainingStop(stop)) {
            if (routeSegment.routeId == route) {
                stopIndex = routeSegment.stopIndex;
                break;
            }
        }
        if (stopIndex == noStopIndex) {
            std::cout << "Stop " << stop << " does not lie on route " << route << "!";
            return;
        } else {
            std::cout << "Stop index " << stopIndex << std::endl;
        }
        StopIndex reverseStopIndex(raptorData[FORWARD].numberOfStopsInRoute(route) - stopIndex - 1);

        std::vector<size_t> sourceCoarseCells = coarsePartition.cellsOfVertex(source);
        std::vector<size_t> sourceFineCells = finePartition.cellsOfVertex(source);
        std::vector<size_t> targetCoarseCells = coarsePartition.cellsOfVertex(target);
        std::vector<size_t> targetFineCells = finePartition.cellsOfVertex(target);
        printVertexCells(source);
        printVertexCells(target);

        std::cout << std::setw(16) << "Route segment:";
        std::cout << " " << flags[FORWARD].isNecessaryRouteSegment(DEPARTURE, route, stopIndex, targetCoarseCells, sourceFineCells);
        std::cout << " " << flags[FORWARD].isNecessaryRouteSegment(ARRIVAL, route, stopIndex, targetCoarseCells, sourceFineCells);
        std::cout << " " << flags[BACKWARD].isNecessaryRouteSegment(DEPARTURE, route, reverseStopIndex, sourceCoarseCells, targetFineCells);
        std::cout << " " << flags[BACKWARD].isNecessaryRouteSegment(ARRIVAL, route, reverseStopIndex, sourceCoarseCells, targetFineCells);
        std::cout << std::endl;

        TripIterator tripIterator = raptorData[FORWARD].getTripIterator(route, stopIndex);
        do {
            const size_t stopEventNumber = tripIterator.stopEvent() - &(raptorData[FORWARD].stopEvents[0]);
            const size_t reverseStopEventNumber = stopEventPermutation[stopEventNumber];
            std::cout << "(" << std::setw(6) << tripIterator.arrivalTime() << "," << std::setw(6) << tripIterator.departureTime() << "):";
            std::cout << " " << flags[FORWARD].isNecessaryStopEvent(DEPARTURE, stopEventNumber, targetCoarseCells, sourceFineCells);
            std::cout << " " << flags[FORWARD].isNecessaryStopEvent(ARRIVAL, stopEventNumber, targetCoarseCells, sourceFineCells);
            std::cout << " " << flags[BACKWARD].isNecessaryStopEvent(DEPARTURE, reverseStopEventNumber, sourceCoarseCells, targetFineCells);
            std::cout << " " << flags[BACKWARD].isNecessaryStopEvent(ARRIVAL, reverseStopEventNumber, sourceCoarseCells, targetFineCells);
            std::cout << std::endl;
        } while (tripIterator.decreaseTrip());
    }

    inline void printShortcuts(const Vertex source, const Vertex target, const StopId from, const StopId to) const noexcept {
        size_t sourceCoarseCell = coarsePartition.cellsOfVertex(source)[0];
        size_t sourceFineCell = finePartition.cellsOfVertex(source)[0];
        size_t targetCoarseCell = coarsePartition.cellsOfVertex(target)[0];
        size_t targetFineCell = finePartition.cellsOfVertex(target)[0];
        printVertexCells(source);
        printVertexCells(target);
        printVertexCells(from);
        printVertexCells(to);

        std::cout << "Forward coarse: " << shortcuts[FORWARD].hasCoarseShortcut(targetCoarseCell, from, to) << std::endl;
        std::cout << "Forward fine: " << shortcuts[FORWARD].hasFineShortcut(targetCoarseCell, sourceFineCell, from, to) << std::endl;
        std::cout << "Backward coarse: " << shortcuts[BACKWARD].hasCoarseShortcut(sourceCoarseCell, from, to) << std::endl;
        std::cout << "Backward fine: " << shortcuts[BACKWARD].hasFineShortcut(sourceCoarseCell, targetFineCell, from, to) << std::endl;
    }

private:
    inline Shortcuts createShortcuts(const TransferGraph& shortcutGraph) const noexcept {
        return Shortcuts(shortcutGraph, coarsePartition.numberOfCells(), finePartition.numberOfCells());
    }

    inline Flags createFlags(const bool direction) const noexcept {
        return Flags(coarsePartition.numberOfCells(), finePartition.numberOfCells(), raptorData[direction]);
    }

private:
    const Data raptorData[2];
    Permutation stopEventPermutation;
    const VertexSeparatorPartition coarsePartition;
    const FinePartitionType finePartition;
    Shortcuts shortcuts[2];
    Flags flags[2];
};

}
