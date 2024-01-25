#pragma once

#include "../../Graph/Graph.h"
#include "../Data.h"
#include "../../../Helpers/IO/Serialization.h"
#include "../../../Helpers/Types.h"

#include <string>
#include <vector>

namespace RAPTOR::RouteFlags {

struct Shortcuts {
    Shortcuts(const TransferGraph& shortcuts, const size_t fromCells, const size_t toCells) :
        allShortcuts(shortcuts),
        coarseCellShortcuts(fromCells),
        fineCellShortcuts(fromCells, std::vector<TransferGraph>(toCells)) {
        for (size_t from = 0; from < fromCells; from++) {
            coarseCellShortcuts[from].addVertices(shortcuts.numVertices());
            for (size_t to = 0; to < toCells; to++) {
                fineCellShortcuts[from][to].addVertices(shortcuts.numVertices());
            }
        }
    }

    Shortcuts(const std::string filename) {
        deserialize(filename);
    }

    Shortcuts(const Shortcuts& other) :
        allShortcuts(other.allShortcuts),
        coarseCellShortcuts(other.coarseCellShortcuts),
        fineCellShortcuts(other.fineCellShortcuts) {
    }

    template<typename COMMON_CELLS>
    inline void buildFromFlags(const COMMON_CELLS& getCommonCells, const std::vector<std::vector<std::vector<bool>>>& flags, const bool revert) noexcept {
        std::vector<DynamicTransferGraph> tempCoarse(coarseCellShortcuts.size());
        for (size_t i = 0; i < coarseCellShortcuts.size(); i++) {
            tempCoarse[i].addVertices(coarseCellShortcuts[i].numVertices());
        }

        for (const auto [edge, from] : allShortcuts.edgesWithFromVertex()) {
            const Vertex to = allShortcuts.get(ToVertex, edge);
            std::vector<size_t> commonCells = getCommonCells(from, to);
            if (commonCells.size() != 1) continue;
            tempCoarse[commonCells[0]].addEdge(revert ? to : from, revert ? from : to, allShortcuts.edgeRecord(edge));
        }

        for (size_t i = 0; i < coarseCellShortcuts.size(); i++) {
            Graph::move(std::move(tempCoarse[i]), coarseCellShortcuts[i]);
        }

        std::vector<std::vector<DynamicTransferGraph>> tempFine(fineCellShortcuts.size(), std::vector<DynamicTransferGraph>(fineCellShortcuts[0].size()));
        for (size_t i = 0; i < fineCellShortcuts.size(); i++) {
            for (size_t j = 0; j < fineCellShortcuts[i].size(); j++) {
                tempFine[i][j].addVertices(fineCellShortcuts[i][j].numVertices());
            }
        }

        for (const auto [edge, from] : allShortcuts.edgesWithFromVertex()) {
            const Vertex to = allShortcuts.get(ToVertex, edge);
            for (size_t i = 0; i < fineCellShortcuts.size(); i++) {
                if (coarseCellShortcuts[i].hasEdge(revert ? to : from, revert ? from : to)) continue;
                for (size_t j = 0; j < fineCellShortcuts[i].size(); j++) {
                    if (flags[i][j][edge]) {
                        tempFine[i][j].addEdge(revert ? to : from, revert ? from : to, allShortcuts.edgeRecord(edge));
                    }
                }
            }
        }

        for (size_t i = 0; i < fineCellShortcuts.size(); i++) {
            for (size_t j = 0; j < fineCellShortcuts[i].size(); j++) {
                Graph::move(std::move(tempFine[i][j]), fineCellShortcuts[i][j]);
            }
        }
    }

    inline void revert() noexcept {
        for (TransferGraph& graph : coarseCellShortcuts) {
            graph.revert();
        }
        for (std::vector<TransferGraph>& v : fineCellShortcuts) {
            for (TransferGraph& graph : v) {
                graph.revert();
            }
        }
    }

    static inline EdgeFlagsTransferGraph packFineCellShortcuts(const std::vector<TransferGraph>& graphs) noexcept{
        EdgeFlagsTransferGraph result;
        result.addVertices(graphs[0].numVertices());
        const size_t numCells = graphs.size();
        for (size_t cell = 0; cell < numCells; cell++) {
            const TransferGraph& graph = graphs[cell];
            for (const auto [edge, from] : graph.edgesWithFromVertex()) {
                const Vertex to = graph.get(ToVertex, edge);
                Edge newEdge = result.findEdge(from, to);
                if (newEdge == noEdge) newEdge = result.addEdge(from, to, graph.edgeRecord(edge));
                std::vector<bool>& flags = result.get(EdgeFlags, newEdge);
                flags.resize(numCells);
                flags[cell] = true;
            }
        }
        return result;
    }

    static inline void unpackFineCellShortcuts(const EdgeFlagsTransferGraph& packedGraph, std::vector<TransferGraph>& unpackedGraphs) noexcept {
        const size_t numCells = unpackedGraphs.size();
        std::vector<DynamicTransferGraph> tempGraphs(numCells);
        for (DynamicTransferGraph& graph : tempGraphs) {
            graph.addVertices(packedGraph.numVertices());
        }
        for (const auto [edge, from] : packedGraph.edgesWithFromVertex()) {
            const Vertex to = packedGraph.get(ToVertex, edge);
            const std::vector<bool>& flags = packedGraph.get(EdgeFlags, edge);
            AssertMsg(flags.size() == numCells, "Flags vector has wrong size, should have " << numCells << " but has " << flags.size() << "!");
            for (size_t cell = 0; cell < numCells; cell++) {
                if (!flags[cell]) continue;
                tempGraphs[cell].addEdge(from, to, packedGraph.edgeRecord(edge));
            }
        }
        for (size_t cell = 0; cell < numCells; cell++) {
            Graph::move(std::move(tempGraphs[cell]), unpackedGraphs[cell]);
        }
    }

    inline void sort() noexcept {
        for (TransferGraph& graph : coarseCellShortcuts) {
            graph.sortEdges(ToVertex);
        }
        for (std::vector<TransferGraph>& graphs : fineCellShortcuts) {
            for (TransferGraph& graph : graphs) {
                graph.sortEdges(ToVertex);
            }
        }
    }

    inline bool hasCoarseShortcut(const size_t coarseCell, const StopId from, const StopId to) const noexcept {
        return coarseCellShortcuts[coarseCell].hasEdge(from, to);
    }

    inline bool hasFineShortcut(const size_t coarseCell, const size_t fineCell, const StopId from, const StopId to) const noexcept {
        return fineCellShortcuts[coarseCell][fineCell].hasEdge(from, to);
    }

    inline void serialize(const std::string filename) const noexcept {
        allShortcuts.writeBinary(filename + ".all");
        IO::serialize(filename, coarseCellShortcuts.size(), fineCellShortcuts[0].size());
        for (size_t i = 0; i < fineCellShortcuts.size(); i++) {
            coarseCellShortcuts[i].writeBinary(filename + ".coarse." + std::to_string(i));
            packFineCellShortcuts(fineCellShortcuts[i]).writeBinary(filename + ".fine." + std::to_string(i));
        }
    }

    inline void deserialize(const std::string filename) noexcept {
        allShortcuts.readBinary(filename + ".all");
        size_t fromCells, toCells;
        IO::deserialize(filename, fromCells, toCells);
        coarseCellShortcuts.resize(fromCells);
        fineCellShortcuts.resize(fromCells);
        for (size_t i = 0; i < fromCells; i++) {
            coarseCellShortcuts[i].readBinary(filename + ".coarse." + std::to_string(i));
            fineCellShortcuts[i].resize(toCells);
            EdgeFlagsTransferGraph packedGraph;
            packedGraph.readBinary(filename + ".fine." + std::to_string(i));
            unpackFineCellShortcuts(packedGraph, fineCellShortcuts[i]);
        }
    }

    inline long long byteSize() const noexcept {
        long long result = allShortcuts.byteSize();
        for (const TransferGraph& graph : coarseCellShortcuts) {
            result += graph.byteSize();
        }
        for (const std::vector<TransferGraph>& v : fineCellShortcuts) {
            for (const TransferGraph& graph : v) {
                result += graph.byteSize();
            }
        }
        return result;
    }

    inline void printCellStatistics(const size_t fromCell, const size_t toCell) const noexcept {
        std::cout << "    #Coarse: " << coarseCellShortcuts[fromCell].numEdges() << std::endl;
        std::cout << "    #Fine:   " << fineCellShortcuts[fromCell][toCell].numEdges() << std::endl;
    }

    inline void printAverageStatistics() const noexcept {
        double coarseAverage = Vector::mean(coarseCellShortcuts, [&](const TransferGraph& graph) {
            return graph.numEdges();
        });
        double fineAverage = Vector::mean(fineCellShortcuts, [&](const std::vector<TransferGraph>& graphs) {
            return Vector::mean(graphs, [&](const TransferGraph& graph) {
                return graph.numEdges();
            });
        });
        std::cout << "    #Coarse: " << coarseAverage << std::endl;
        std::cout << "    #Fine: " << fineAverage << std::endl;
    }

    TransferGraph allShortcuts;
    std::vector<TransferGraph> coarseCellShortcuts;
    std::vector<std::vector<TransferGraph>> fineCellShortcuts;
};

}
