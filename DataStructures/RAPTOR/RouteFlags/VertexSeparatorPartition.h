#pragma once

#include "../../Container/Set.h"
#include "../../Graph/Graph.h"
#include "../../Partition/NestedDissection.h"
#include "../Data.h"
#include "../../../Helpers/IO/Serialization.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/Vector/Vector.h"

#include <iostream>
#include <vector>

namespace RAPTOR::RouteFlags {

struct VertexSeparatorPartition {
    VertexSeparatorPartition(const NestedDissection& partition, const Data& raptorData, const size_t level) :
        partition(partition),
        level(level),
        incidentCells(partition.getIncidentCells(level, true)),
        separatorVertices(raptorData.transferGraph.numVertices(), partition.getSeparatorOfLevel(level)) {
    }

    VertexSeparatorPartition(const std::string& filename) {
        deserialize(filename);
    }

    inline size_t highestBoundaryVertexId() const noexcept {
        return Vector::max(boundaryVertices());
    }

    inline size_t numberOfBoundaryVertices() const noexcept {
        return separatorVertices.size();
    }

    inline const std::vector<Vertex>& boundaryVertices() const noexcept {
        return separatorVertices.getValues();
    }

    inline bool isBoundaryVertex(const Vertex vertex) const noexcept {
        return separatorVertices.contains(vertex);
    }

    inline size_t numberOfCells() const noexcept {
        return partition.cellsOfLevel(level).size();
    }

    inline const std::vector<size_t>& cellsOfVertex(const Vertex u) const noexcept {
        AssertMsg(u < partition.numberOfVertices(), "Invalid vertex: " << u);
        return incidentCells[u];
    }

    inline bool isInCell(const Vertex u, const size_t cell) const noexcept {
        return Vector::contains(cellsOfVertex(u), cell);
    }

    inline std::vector<size_t> getCommonCells(const Vertex u, const Vertex v) const noexcept {
        return Vector::sortedIntersection(cellsOfVertex(u), cellsOfVertex(v));
    }

    inline size_t getFirstCommonCell(const Vertex u, const Vertex v) const noexcept {
        const std::vector<size_t>& uCells = cellsOfVertex(u);
        const std::vector<size_t>& vCells = cellsOfVertex(v);
        size_t j = 0;
        for (size_t i = 0; i < uCells.size(); i++) {
            while (j < vCells.size() && vCells[j] < uCells[i]) j++;
            if (j == vCells.size()) return -1;
            if (vCells[j] == uCells[i]) return vCells[j];
        }
        return -1;
    }

    inline bool areInCommonCell(const Vertex u, const Vertex v) const noexcept {
        return getFirstCommonCell(u, v) != size_t(-1);
    }

    inline size_t getNonBoundaryIntraCellEdgeCellId(const Vertex u, const Vertex v) const noexcept {
        const std::vector<size_t>& uCells = cellsOfVertex(u);
        const std::vector<size_t>& vCells = cellsOfVertex(v);
        if (uCells.size() == 1 && Vector::contains(vCells, uCells[0])) return uCells[0];
        if (vCells.size() == 1 && Vector::contains(uCells, vCells[0])) return vCells[0];
        return -1;
    }


    inline bool isNonBoundaryIntraCellEdge(const Vertex u, const Vertex v) const noexcept {
        return getNonBoundaryIntraCellEdgeCellId(u, v) != size_t(-1);
    }

    inline void serialize(const std::string& filename) const noexcept {
        IO::serialize(filename, level, incidentCells, separatorVertices);
        partition.serialize(filename + ".partition");
    }

    inline void deserialize(const std::string& filename) noexcept {
        IO::deserialize(filename, level, incidentCells, separatorVertices);
        partition.deserialize(filename + ".partition");
    }

    inline long long byteSize() const noexcept {
        long long result = partition.byteSize();
        result += sizeof(int);
        result += Vector::byteSize(incidentCells);
        result += separatorVertices.byteSize();
        return result;
    }

    inline void printStatistics() const noexcept {
        std::cout << "    Nested dissection with vertex separators" << std::endl;
        std::cout << "    Level: " << level << std::endl;
        std::cout << "    #Cells: " << numberOfCells() << std::endl;
        std::cout << "    #Boundary vertices: " << separatorVertices.size() << std::endl;
    }

    NestedDissection partition;
    size_t level;
    std::vector<std::vector<size_t>> incidentCells;
    IndexedSet<false, Vertex> separatorVertices;
};

}
