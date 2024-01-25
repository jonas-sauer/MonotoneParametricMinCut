#pragma once

#include "HyperGraph.h"

#include "../Data.h"

#include "../../Container/Set.h"

#include "../../../Helpers/String/String.h"

namespace RAPTOR {

class Partition {

public:
    Partition(const HyperGraph& hyperGraph, const std::string& filename, const bool verbose) {
        if (verbose) std::cout << "Reading hyper-graph partition from: " << filename << std::endl << std::flush;
        std::ifstream partIs(filename);
        AssertMsg(partIs.is_open(), "cannot open file: " << filename);
        cellOfVertice.reserve(hyperGraph.numVertices());
        while (!partIs.eof()) {
            std::string line;
            getline(partIs, line);
            if (line == "") continue;
            const size_t cell = String::lexicalCast<size_t>(line);
            if (cell >= hyperVerticesByCell.size()) hyperVerticesByCell.resize(cell + 1);
            hyperVerticesByCell[cell].emplace_back(cellOfVertice.size());
            cellOfVertice.emplace_back(cell);
        }
        AssertMsg(cellOfVertice.size() == hyperGraph.vertices.size(), "Partition contains " << cellOfVertice.size() << " lines, but should contain " << hyperGraph.vertices.size() << " lines!");
        hyperEdgesByCell.resize(hyperVerticesByCell.size());
        cellsOfEdge.resize(hyperGraph.numEdges());
        for (size_t hyperEdgeId = 0; hyperEdgeId < hyperGraph.edges.size(); hyperEdgeId++) {
            Set<size_t> cells;
            for (const HyperVertexId hyperVertexId : hyperGraph.edges[hyperEdgeId].hyperVertices) {
                cells.insert(cellOfVertice[hyperVertexId]);
            }
            for (const size_t cell : cells) {
                cellsOfEdge[hyperEdgeId].emplace_back(cell);
            }
            AssertMsg(cells.size() > 0, "Hyper-edge " << hyperEdgeId << " is incident to no Cell of the partition!");
            if (cells.size() == 1) {
                hyperEdgesByCell[cellsOfEdge[hyperEdgeId][0]].emplace_back(hyperEdgeId);
            } else {
                cutHyperEdges.emplace_back(hyperEdgeId);
            }
        }
    }

public:
    inline size_t numberOfCells() const noexcept {
        return hyperVerticesByCell.size();
    }

    inline void printInfo() const noexcept {
        std::cout << "Hyper graph partition:" << std::endl;
        std::cout << "   Number of Cells:          " << std::setw(12) << String::prettyInt(numberOfCells()) << std::endl;
        std::cout << "   Number of cut edges:      " << std::setw(12) << String::prettyInt(cutHyperEdges.size()) << std::endl;
    }

public:
    std::vector<std::vector<size_t>> hyperVerticesByCell;
    std::vector<size_t> cellOfVertice;

    std::vector<std::vector<size_t>> hyperEdgesByCell;
    std::vector<size_t> cutHyperEdges;
    std::vector<std::vector<size_t>> cellsOfEdge;

};

}