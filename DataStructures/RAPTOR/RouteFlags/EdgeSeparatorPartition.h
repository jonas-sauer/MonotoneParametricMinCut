#pragma once

#include "../../Container/Set.h"
#include "../../Graph/Graph.h"
#include "../../Partition/VertexPartition.h"
#include "../Data.h"
#include "../../../Helpers/IO/Serialization.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/Vector/Vector.h"

#include <iostream>
#include <vector>

namespace RAPTOR::RouteFlags {

struct EdgeSeparatorPartition {
    EdgeSeparatorPartition(const VertexPartition& partition, const Data& raptorData) :
        partition(partition),
        borderVertices(raptorData.transferGraph.numVertices(), partition.getBorderVertices(raptorData.minTravelTimeGraph())) {
    }

    EdgeSeparatorPartition(const std::string& filename) {
        deserialize(filename);
    }

    inline size_t highestBoundaryVertexId() const noexcept {
        return Vector::max(boundaryVertices());
    }

    inline size_t numberOfBoundaryVertices() const noexcept {
        return borderVertices.size();
    }

    inline const std::vector<Vertex>& boundaryVertices() const noexcept {
        return borderVertices.getValues();
    }

    inline bool isBoundaryVertex(const Vertex vertex) const noexcept {
        return borderVertices.contains(vertex);
    }

    inline size_t numberOfCells() const noexcept {
        return partition.numberOfCells();
    }

    inline size_t cellOfVertex(const Vertex u) const noexcept {
        AssertMsg(u < partition.numberOfVertices(), "Invalid vertex: " << u);
        return partition.getCellIdOfVertex(u);
    }

    inline const std::vector<size_t> cellsOfVertex(const Vertex u) const noexcept {
        AssertMsg(u < partition.numberOfVertices(), "Invalid vertex: " << u);
        return std::vector<size_t>{(size_t)partition.getCellIdOfVertex(u)};
    }

    inline bool isInCell(const Vertex u, const size_t cell) const noexcept {
        return cellOfVertex(u) == cell;
    }

    inline std::vector<size_t> getCommonCells(const Vertex u, const Vertex v) const noexcept {
        const size_t commonCell = getFirstCommonCell(u, v);
        std::vector<size_t> result;
        if (commonCell != size_t(-1)) result.emplace_back(commonCell);
        return result;
    }

    inline size_t getFirstCommonCell(const Vertex u, const Vertex v) const noexcept {
        const size_t uCell = cellOfVertex(u);
        const size_t vCell = cellOfVertex(v);
        return (uCell == vCell) ? uCell : -1;
    }

    inline bool areInCommonCell(const Vertex u, const Vertex v) const noexcept {
        return cellOfVertex(u) == cellOfVertex(v);
    }

    inline void serialize(const std::string& filename) const noexcept {
        IO::serialize(filename, borderVertices);
        partition.serialize(filename + ".partition");
    }

    inline void deserialize(const std::string& filename) noexcept {
        IO::deserialize(filename, borderVertices);
        partition.deserialize(filename + ".partition");
    }

    inline long long byteSize() const noexcept {
        return partition.byteSize() + borderVertices.byteSize();
    }

    inline void printStatistics() const noexcept {
        std::cout << "    Vertex partition with edge separators" << std::endl;
        std::cout << "    #Cells: " << numberOfCells() << std::endl;
        std::cout << "    #Boundary vertices: " << borderVertices.size() << std::endl;
    }

    VertexPartition partition;
    IndexedSet<false, Vertex> borderVertices;
};

}
