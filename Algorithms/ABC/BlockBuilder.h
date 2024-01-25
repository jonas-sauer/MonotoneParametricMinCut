#pragma once

#include "../UnionFind.h"
#include "../DepthFirstSearch.h"

#include "../../Helpers/Types.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Container/Map.h"
#include "../../DataStructures/Container/Heap.h"

namespace ABC {

class BlockBuilder {

private:
    static inline constexpr uint32_t undefined = -1;

    struct ColumnLable {
        ColumnLable(const size_t numberOfStops, const size_t numberOfTrips) :
            arrivalTimeAtStop(numberOfStops),
            arrivalInTrip(numberOfTrips) {
        }
        inline int getArrivalTimeAtStop(const StopId stop) noexcept {
            return arrivalTimeAtStop.contains(stop) ? arrivalTimeAtStop[stop] : never;
        }
        inline void setArrivalTimeAtStop(const StopId stop, const int time) noexcept {
            if (arrivalTimeAtStop.contains(stop)) {
                arrivalTimeAtStop[stop] = std::min(arrivalTimeAtStop[stop], time);
            } else {
                arrivalTimeAtStop.insert(stop, time);
            }
        }
        inline bool getArrivalInTrip(const TripId trip) noexcept {
            return arrivalInTrip.contains(trip);
        }
        inline void setArrivalInTrip(const TripId trip) noexcept {
            arrivalInTrip.insert(trip);
        }
        IndexedMap<int, false, StopId> arrivalTimeAtStop;
        IndexedSet<false, TripId> arrivalInTrip;
    };

    struct ArrivalLabel {
        ArrivalLabel(const ConnectionId connection, const int time) :
            connection(connection),
            time(time) {
        }
        inline bool operator<(const ArrivalLabel& other) const noexcept {
            return time < other.time;
        }
        ConnectionId connection;
        int time;
    };

    struct PathCover {
        PathCover(const uint32_t numberOfBlocks = 0) :
            numberOfPaths(0),
            pathOfBlock(numberOfBlocks, numberOfBlocks),
            indexOfBlock(numberOfBlocks, numberOfBlocks) {
        }
        uint32_t numberOfPaths;
        std::vector<uint32_t> pathOfBlock;
        std::vector<uint32_t> indexOfBlock;
    };

public:
    BlockBuilder(CSA::Data& csa, const std::vector<uint32_t>& columnOfStop, const uint32_t targetBlockSize, const int maxBlockTimeDelta, const bool verbose = true, const bool = false) :
        csa(csa),
        columnOfStop(columnOfStop),
        targetBlockSize(targetBlockSize),
        maxBlockTimeDelta(maxBlockTimeDelta),
        blockOfConnection(csa.numberOfConnections(), undefined),
        reachedColumnsOfStop(csa.numberOfStops()),
        reachedColumnsOfTrip(csa.numberOfTrips()) {
        csa.sortConnectionsAscendingByDepartureTime();
        Timer timer;
        if (verbose) std::cout << "   Building blocks from columns... " << std::setw(14) << std::flush;
        findBlocks();
        if (verbose) std::cout << String::msToString(timer.elapsedMilliseconds()) << "\n   Building block graph...   " << std::setw(20) << std::flush;
        if (verbose) timer.restart();
        buildBlockGraph(false);
        if (verbose) std::cout << String::msToString(timer.elapsedMilliseconds()) << "\n   Removing transitive edges... " << std::setw(17) << std::flush;
        if (verbose) timer.restart();
        removeTransitiveEdges(verbose);
        if (verbose) std::cout << String::msToString(timer.elapsedMilliseconds()) << std::endl;
    }

    inline size_t numberOfConnections() const noexcept {
        size_t result = 0;
        for (const std::vector<ConnectionId>& connections : connectionsOfBlock) {
            result += connections.size();
        }
        return result;
    }

    inline size_t numberOfColumns() const noexcept {
        return blocksOfColumn.size();
    }

    inline const std::vector<ConnectionId>& getConnectionsOfBlock(const uint32_t block) const noexcept {
        return connectionsOfBlock[block];
    }

    inline const std::vector<uint32_t> getColumnOfBlocks() const noexcept {
        std::vector<uint32_t> columnOfBlock(numberOfBlocks());
        for (size_t column = 0; column < numberOfColumns(); column++) {
            for (const uint32_t block : blocksOfColumn[column]) {
                columnOfBlock[block] = column;
            }
        }
        return columnOfBlock;
    }

    inline size_t numberOfBlocks() const noexcept {
        return connectionsOfBlock.size();
    }

    inline std::vector<size_t> blockSizes() const noexcept {
        std::vector<size_t> result;
        for (const std::vector<ConnectionId>& connections : connectionsOfBlock) {
            result.emplace_back(connections.size());
        }
        return result;
    }

    inline std::vector<size_t> blockDegrees() const noexcept {
        std::vector<size_t> result;
        for (const Vertex block : dynamicBlockGraph.vertices()) {
            result.emplace_back(dynamicBlockGraph.outDegree(block));
        }
        return result;
    }

    inline const std::vector<std::vector<uint32_t>>& getBlocksOfColumn() const noexcept {
        return blocksOfColumn;
    }

    inline const SimpleDynamicGraph& getDynamicBlockGraph() const noexcept {
        return dynamicBlockGraph;
    }

    inline SimpleDynamicGraph getForwardBlockGraph() const noexcept {
        return dynamicBlockGraph;
    }

    inline SimpleDynamicGraph getBackwardBlockGraph() const noexcept {
        SimpleDynamicGraph result;
        result.addVertices(dynamicBlockGraph.numVertices());
        for (const Vertex from : dynamicBlockGraph.vertices()) {
            for (const Edge edge : dynamicBlockGraph.edgesFrom(from)) {
                const Vertex to = dynamicBlockGraph.get(ToVertex, edge);
                result.addEdge(to, from);
            }
        }
        return result;
    }

    inline uint32_t getPathCoverSize() const noexcept {
        return pathCover.numberOfPaths;
    }

private:
    inline void findBlocks() noexcept {
        for (ConnectionId i = ConnectionId(0); i < csa.numberOfConnections(); i++) {
            const CSA::Connection& connection = csa.connections[i];
            const uint32_t departureColumn = columnOfStop[connection.departureStopId];
            while (departureColumn >= currentBlockOfColumn.size()) {
                currentBlockOfColumn.emplace_back(undefined);
                blocksOfColumn.emplace_back();
                labelOfColumn.emplace_back(csa.numberOfStops(), csa.numberOfTrips());
            }
            if (blockIsFull(currentBlockOfColumn[departureColumn], connection.departureTime)) {
                currentBlockOfColumn[departureColumn] = connectionsOfBlock.size();
                blocksOfColumn[departureColumn].emplace_back(currentBlockOfColumn[departureColumn]);
                connectionsOfBlock.emplace_back();
                ColumnLable& label = labelOfColumn[departureColumn];
                for (const StopId stop : label.arrivalTimeAtStop.getKeys()) {
                    reachedColumnsOfStop[stop].erase(departureColumn);
                }
                label.arrivalTimeAtStop.clear();
                label.arrivalInTrip.clear();
            }
            blockOfConnection[i] = currentBlockOfColumn[departureColumn];

            Set<uint32_t> departureColumns;
            departureColumns.insert(departureColumn);
            for (const uint32_t column : reachedColumnsOfStop[connection.departureStopId]) {
                if (labelOfColumn[column].getArrivalTimeAtStop(connection.departureStopId) <= connection.departureTime) {
                    departureColumns.insert(column);
                }
            }
            for (const uint32_t column : reachedColumnsOfTrip[connection.tripId]) {
                if (labelOfColumn[column].getArrivalInTrip(connection.tripId)) {
                    departureColumns.insert(column);
                }
            }

            for (const uint32_t column : departureColumns) {
                const uint32_t block = currentBlockOfColumn[column];
                if (blockIsFull(block, connection.departureTime)) continue;
                connectionsOfBlock[block].emplace_back(i);
                labelOfColumn[column].setArrivalTimeAtStop(connection.arrivalStopId, connection.arrivalTime);
                reachedColumnsOfStop[connection.arrivalStopId].insert(column);
                labelOfColumn[column].setArrivalInTrip(connection.tripId);
                reachedColumnsOfTrip[connection.tripId].insert(column);
                for (const Edge edge : csa.transferGraph.edgesFrom(connection.arrivalStopId)) {
                    const StopId arrivalStopId = StopId(csa.transferGraph.get(ToVertex, edge));
                    labelOfColumn[column].setArrivalTimeAtStop(arrivalStopId, connection.arrivalTime + csa.transferGraph.get(TravelTime, edge));
                    reachedColumnsOfStop[arrivalStopId].insert(column);
                }
            }
        }
    }

    inline void buildBlockGraph(const bool useLatestTargetBlock = false, const bool addColumnEdges = true) noexcept {
        std::vector<std::vector<uint32_t>> blocksOfConnection(csa.numberOfConnections());
        std::vector<uint32_t> latestBlockOfConnection = blockOfConnection;
        for (uint32_t block = 0; block < connectionsOfBlock.size(); block++) {
            for (const ConnectionId connection : connectionsOfBlock[block]) {
                blocksOfConnection[connection].emplace_back(block);
                if (useLatestTargetBlock && (minTimeOfBlock(latestBlockOfConnection[connection]) < minTimeOfBlock(block))) {
                    latestBlockOfConnection[connection] = block;
                }
            }
        }
        dynamicBlockGraph.clear();
        dynamicBlockGraph.addVertices(connectionsOfBlock.size());
        std::vector<ConnectionId> previousDepartureOfStop(csa.numberOfStops(), noConnection);
        std::vector<ConnectionId> previousConnectionOfTrip(csa.numberOfTrips(), noConnection);
        std::vector<Heap<ArrivalLabel>> previousArrivalsOfStop(csa.numberOfStops());
        for (ConnectionId i = ConnectionId(0); i < csa.numberOfConnections(); i++) {
            const CSA::Connection& connection = csa.connections[i];
            if (previousDepartureOfStop[connection.departureStopId] != noConnection) {
                for (const uint32_t block : blocksOfConnection[previousDepartureOfStop[connection.departureStopId]]) {
                    if (!Vector::contains(blocksOfConnection[i], block)) {
                        dynamicBlockGraph.findOrAddEdge(Vertex(block), Vertex(latestBlockOfConnection[i]));
                    }
                }
            }
            previousDepartureOfStop[connection.departureStopId] = i;
            if (previousConnectionOfTrip[connection.tripId] != noConnection) {
                for (const uint32_t block : blocksOfConnection[previousConnectionOfTrip[connection.tripId]]) {
                    if (!Vector::contains(blocksOfConnection[i], block)) {
                        dynamicBlockGraph.findOrAddEdge(Vertex(block), Vertex(latestBlockOfConnection[i]));
                    }
                }
            }
            previousConnectionOfTrip[connection.tripId] = i;
            while ((!previousArrivalsOfStop[connection.departureStopId].empty()) && (previousArrivalsOfStop[connection.departureStopId].min().time <= connection.departureTime)) {
                for (const uint32_t block : blocksOfConnection[previousArrivalsOfStop[connection.departureStopId].min().connection]) {
                    if (!Vector::contains(blocksOfConnection[i], block)) {
                        dynamicBlockGraph.findOrAddEdge(Vertex(block), Vertex(latestBlockOfConnection[i]));
                    }
                }
                previousArrivalsOfStop[connection.departureStopId].remove_min();
            }
            previousArrivalsOfStop[connection.arrivalStopId].emplace_back(i, connection.arrivalTime);
            for (const Edge edge : csa.transferGraph.edgesFrom(connection.arrivalStopId)) {
                previousArrivalsOfStop[csa.transferGraph.get(ToVertex, edge)].emplace_back(i, connection.arrivalTime + csa.transferGraph.get(TravelTime, edge));
            }
        }
        if (addColumnEdges) {
            for (uint32_t column = 0; column < blocksOfColumn.size(); column++) {
                for (size_t i = 1; i < blocksOfColumn[column].size(); i++) {
                    dynamicBlockGraph.findOrAddEdge(Vertex(blocksOfColumn[column][i - 1]), Vertex(blocksOfColumn[column][i]));
                }
            }
        }
    }

    inline void removeTransitiveEdges(const bool debug = true, const bool keepMinTimeEdges = true) noexcept {
        SimpleDynamicGraph backup;
        if (debug) {
            backup = dynamicBlockGraph;
            std::cout << "\n      Original Number of Edges: " << String::prettyInt(backup.numEdges()) << std::flush;
        }

        std::vector<Vertex> topologicalOrder;
        std::vector<uint32_t> topologicalIndex(numberOfBlocks(), 0);
        dfs(dynamicBlockGraph, NoOperation, NoOperation, [&](const Vertex v){
            topologicalIndex[v] = topologicalOrder.size();
            topologicalOrder.emplace_back(v);
        });
        pathCover = findPathCover(topologicalOrder, topologicalIndex);
        std::vector<std::vector<uint32_t>> exclusiveTransitiveClosureOfBlock(numberOfBlocks(), std::vector<uint32_t>(pathCover.numberOfPaths, numberOfBlocks()));
        std::vector<bool> deleteEdge(dynamicBlockGraph.edgeLimit(), false);
        for (const Vertex block : topologicalOrder) {
            for (const Edge edge : dynamicBlockGraph.edgesFrom(block)) {
                const Vertex otherBlock = dynamicBlockGraph.get(ToVertex, edge);
                for (size_t path = 0; path < pathCover.numberOfPaths; path++) {
                    exclusiveTransitiveClosureOfBlock[block][path] = std::min(exclusiveTransitiveClosureOfBlock[block][path], exclusiveTransitiveClosureOfBlock[otherBlock][path]);
                }
            }
            for (const Edge edge : dynamicBlockGraph.edgesFrom(block)) {
                const Vertex otherBlock = dynamicBlockGraph.get(ToVertex, edge);
                const uint32_t path = pathCover.pathOfBlock[otherBlock];
                if (exclusiveTransitiveClosureOfBlock[block][path] > pathCover.indexOfBlock[otherBlock]) {
                    exclusiveTransitiveClosureOfBlock[block][path] = pathCover.indexOfBlock[otherBlock];
                } else {
                    if ((!keepMinTimeEdges) || (minTimeOfBlock(block) < minTimeOfBlock(otherBlock))) {
                        deleteEdge[edge] = true;
                    }
                }
            }
        }
        dynamicBlockGraph.deleteEdges(deleteEdge, true);
        dynamicBlockGraph.packEdges();

        if (debug) {
            std::cout << "\n      Path cover size: " << String::prettyInt(pathCover.numberOfPaths) << std::flush;
            std::cout << "\n      Testing correctness... " << std::setw(20) << std::flush;
            for (const Vertex block : topologicalOrder) {
                std::vector<bool> reachable(topologicalOrder.size(), false);
                reachable[block] = true;
                for (size_t i = topologicalOrder.size() - 1; i < topologicalOrder.size(); i--) {
                    const Vertex v = Vertex(topologicalOrder[i]);
                    if (!reachable[v]) continue;
                    for (const Edge edge : dynamicBlockGraph.edgesFrom(v)) {
                        reachable[dynamicBlockGraph.get(ToVertex, edge)] = true;
                    }
                }
                for (const Edge edge : backup.edgesFrom(block)) {
                    if (!reachable[backup.get(ToVertex, edge)]) {
                        warning("ERROR: Edge from: ", block, " to: ", backup.get(ToVertex, edge), " is missing!");
                    }
                }
            }
        }
    }

    inline PathCover findPathCover(const std::vector<Vertex>& topologicalOrder, const std::vector<uint32_t> topologicalIndex) const noexcept {
        PathCover result(numberOfBlocks());
        for (size_t i = topologicalOrder.size() - 1; i < topologicalOrder.size(); i--) {
            if (result.pathOfBlock[topologicalOrder[i]] <= result.numberOfPaths) continue;
            uint32_t currentIndex = 0;
            uint32_t block = topologicalOrder[i];
            while (result.pathOfBlock[block] > result.numberOfPaths) {
                result.pathOfBlock[block] = result.numberOfPaths;
                result.indexOfBlock[block] = currentIndex;
                currentIndex++;
                uint32_t nextBlock = block;
                for (const Edge edge : dynamicBlockGraph.edgesFrom(Vertex(block))) {
                    const uint32_t otherBlock = dynamicBlockGraph.get(ToVertex, edge);
                    if (result.pathOfBlock[otherBlock] <= result.numberOfPaths) continue;
                    if ((nextBlock == block) || (topologicalIndex[otherBlock] > topologicalIndex[nextBlock])) {
                        nextBlock = otherBlock;
                    }
                }
                block = nextBlock;
            }
            result.numberOfPaths++;
        }
        return result;
    }

    inline bool blockIsFull(const uint32_t block, const int time) const noexcept {
        if (block >= connectionsOfBlock.size()) return true;
        if (connectionsOfBlock[block].empty()) return false;
        if (csa.connections[connectionsOfBlock[block].back()].departureTime >= time) return false;
        if (connectionsOfBlock[block].size() >= targetBlockSize) return true;
        if (time - minTimeOfBlock(block) > maxBlockTimeDelta) return true;
        return false;
    }

    inline int minTimeOfBlock(const uint32_t block) const noexcept {
        AssertMsg(block < connectionsOfBlock.size(), "The block id " << block << " is out of bounds (0, " << connectionsOfBlock.size() << ")!");
        AssertMsg(!connectionsOfBlock[block].empty(), "The block " << block << " is empty!");
        return csa.connections[connectionsOfBlock[block][0]].departureTime;
    }

private:
    CSA::Data& csa;

    const std::vector<uint32_t>& columnOfStop;
    const uint32_t targetBlockSize;
    const int maxBlockTimeDelta;

    std::vector<std::vector<ConnectionId>> connectionsOfBlock;
    std::vector<uint32_t> blockOfConnection;
    std::vector<uint32_t> currentBlockOfColumn;
    std::vector<std::vector<uint32_t>> blocksOfColumn;

    std::vector<ColumnLable> labelOfColumn;
    std::vector<Set<uint32_t>> reachedColumnsOfStop;
    std::vector<Set<uint32_t>> reachedColumnsOfTrip;

    SimpleDynamicGraph dynamicBlockGraph;
    PathCover pathCover;

};

}
