#pragma once

#include "../../../Dijkstra/Dijkstra.h"

#include "../../../../DataStructures/Container/Map.h"
#include "../../../../DataStructures/RAPTOR/Data.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/PreprocessingJourney.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../../../DataStructures/RAPTOR/RouteFlags/ShortcutFlags.h"

namespace RAPTOR::RouteFlags {

template<typename FLAGS, typename FINE_PARTITION_TYPE, bool DIRECTION = FORWARD>
class FlagsMarker {

public:
    using Flags = FLAGS;
    using FinePartitionType = FINE_PARTITION_TYPE;
    constexpr static bool Direction = DIRECTION;
    using Type = FlagsMarker<Flags, FinePartitionType, Direction>;
    using RouteFlagsDataType = RouteFlagsData<Flags, FinePartitionType>;

private:
    struct ShortcutOriginData {
        int travelTime;
        std::vector<size_t> targetCells;
    };

public:
    FlagsMarker(RouteFlagsDataType& flagsData) :
        flagsData(flagsData),
        raptorData(flagsData.getRaptorData(Direction)),
        shortcutFlags(flagsData.getCoarsePartition().numberOfCells(), flagsData.getFinePartition().numberOfCells(), flagsData.getShortcutGraph(!Direction)),
        sourceVertex(noVertex),
        dijkstra(raptorData.transferGraph, raptorData.transferGraph[TravelTime]),
        shortcutOrigins(raptorData.numberOfStops()) {
    }

    inline void setSource(const Vertex source, const std::vector<size_t>& cells) noexcept {
        sourceVertex = source;
        sourceCells = cells;
    }

    inline void clear() noexcept {
        shortcutOrigins.clear();
    }

    inline const ShortcutFlags& getShortcutFlags() const noexcept {
        return shortcutFlags;
    }

    inline void markJourney(const PreprocessingJourney& journey, const std::vector<size_t>& targetCells) noexcept {
        for (const StopEvent* stopEvent : journey.departureStopEvents) {
            markStopEvent(DEPARTURE, stopEvent, targetCells);
        }
        for (const StopEvent* stopEvent : journey.arrivalStopEvents) {
            markStopEvent(ARRIVAL, stopEvent, targetCells);
        }

        for (const Shortcut& shortcut : journey.shortcuts) {
            markShortcut(shortcut, targetCells);
        }
    }

    inline void markStopEvent(const bool departureOrArrival, const StopEvent* stopEvent, const size_t targetCell) noexcept {
        markStopEvent(departureOrArrival, stopEvent, std::vector<size_t>{targetCell});
    }

    inline void markStopEvent(const bool departureOrArrival, const StopEvent* stopEvent, const std::vector<size_t>& targetCells) noexcept {
        const size_t stopEventId = stopEvent - &(flagsData.getRaptorData(FORWARD).stopEvents[0]);
        flagsData.markStopEvent(!Direction, departureOrArrival, stopEventId, Direction == FORWARD, sourceCells, targetCells);
    }

    inline void markShortcut(const Shortcut& shortcut, const size_t targetCell) noexcept {
        markShortcut(shortcut, std::vector<size_t>{targetCell});
    }

    inline void markShortcut(const Shortcut& shortcut, const std::vector<size_t>& targetCells) noexcept {
        AssertMsg(raptorData.isStop(shortcut.destination), "Destination " << shortcut.destination << " of shortcut must be a stop!");
        AssertMsg(shortcut.origin == sourceVertex || raptorData.isStop(shortcut.origin), "Origin " << shortcut.origin << " of shortcut must be source vertex or a stop!");
        if (raptorData.isStop(shortcut.origin)) {
            shortcutFlags.markShortcut(sourceCells, targetCells, shortcut.destination, shortcut.origin);
        }
        if (shortcut.origin == sourceVertex && !Vector::equals(flagsData.getCoarsePartition().cellsOfVertex(shortcut.destination), sourceCells)) {
            const StopId destination = StopId(shortcut.destination);
            if (shortcutOrigins.contains(destination)) {
                AssertMsg(shortcutOrigins[destination].travelTime == shortcut.travelTime, "Different travel time! Shortcut " << shortcut.origin << " -> " << shortcut.destination << ", old travel time: " << shortcutOrigins[destination].travelTime << ", new travel time: " << shortcut.travelTime);
                shortcutOrigins[destination].targetCells = Vector::sortedUnion(shortcutOrigins[destination].targetCells, targetCells);
            } else {
                shortcutOrigins.insert(destination, ShortcutOriginData{shortcut.travelTime, targetCells});
            }
        }
    }

    inline void reconstructShortcutsIntoCell() noexcept {
        IndexedSet<false, Vertex> shortcutDestinationCandidates(raptorData.transferGraph.numVertices());
        std::vector<Shortcut> shortcutCandidates;
        for (const StopId origin : shortcutOrigins.getKeys()) {
            for (const Edge edge : shortcutFlags.shortcuts.edgesFrom(origin)) {
                const Vertex destination = shortcutFlags.shortcuts.get(ToVertex, edge);
                const int travelTime = shortcutFlags.shortcuts.get(TravelTime, edge);
                AssertMsg(raptorData.isStop(destination), "Shortcut destination " << destination << " is not a stop!");
                if (flagsData.getCoarsePartition().isNonBoundaryIntraCellEdge(destination, sourceVertex)) {
                    shortcutCandidates.emplace_back(Shortcut{origin, destination, travelTime});
                    shortcutDestinationCandidates.insert(destination);
                }
            }
        }
        dijkstra.run(sourceVertex, shortcutDestinationCandidates, NoOperation, NoOperation, [&](const Vertex, const Edge edge) {
            const Vertex to = raptorData.transferGraph.get(ToVertex, edge);
            return !flagsData.getCoarsePartition().areInCommonCell(to, sourceVertex);
        });

        for (const Shortcut& shortcut : shortcutCandidates) {
            if (shortcutOrigins[StopId(shortcut.origin)].travelTime + dijkstra.getDistance(shortcut.destination) == shortcut.travelTime) {
                for (const size_t cell : sourceCells) {
                    if (flagsData.getCoarsePartition().isInCell(shortcut.origin, cell)) continue;
                    if (!flagsData.getCoarsePartition().isInCell(shortcut.destination, cell)) continue;
                    shortcutFlags.markShortcut(cell, shortcutOrigins[StopId(shortcut.origin)].targetCells, shortcut.origin, shortcut.destination);
                }
            }
        }
    }

private:
    RouteFlagsDataType& flagsData;
    const Data& raptorData;
    ShortcutFlags shortcutFlags;

    Vertex sourceVertex;
    std::vector<size_t> sourceCells;

    Dijkstra<TransferGraph> dijkstra;
    IndexedMap<ShortcutOriginData, false, StopId> shortcutOrigins;
};

}
