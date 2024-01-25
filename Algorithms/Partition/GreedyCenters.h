#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../Dijkstra/Dijkstra.h"

#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Console/ProgressBar.h"
#include "../../Helpers/MultiThreading.h"

#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/Partition/VertexPartition.h"

template<bool DEBUG>
class GreedyCenters {

public:
    static constexpr bool Debug = DEBUG;
    using Type = GreedyCenters<Debug>;

public:
    GreedyCenters(const RAPTOR::Data& data, const size_t numberOfCenters = 0) :
        graphCopy(data.minTravelTimeGraph()),
        graph(graphCopy),
        dijkstra(graph, graph[TravelTime]),
        stopValuesCopy(getHeuristicStopValues(data)),
        stopValues(stopValuesCopy) {
        run(numberOfCenters);
    }

    GreedyCenters(const RAPTOR::TransferGraph& transferGraph, const std::vector<int>& stopValues, const size_t numberOfCenters = 0) :
        graph(transferGraph),
        dijkstra(graph, graph[TravelTime]),
        stopValues(stopValues) {
        run(numberOfCenters);
    }

public:
    inline void run(const size_t numberOfCenters, const double cellCoverage = 0.5) noexcept {
        if (Debug) std::cout << "Computing " << numberOfCenters << " Greedy Centers..." << std::endl;
        ProgressBar bar(numberOfCenters, Debug);
        centers.clear();
        std::vector<bool> covered(stopValues.size(), false);
        while(centers.size() < numberOfCenters) {
            centers.emplace_back(bestUncoveredStop(covered));
            bar++;
            covered = currentCover(covered.size() * cellCoverage);
        }
        if (Debug) std::cout << "Done! (" << centers.size() << " centers)" << std::endl;
    }

    inline void run() noexcept {
        run(stopValues.size());
    }

    inline VertexPartition getVertexPartition(const size_t numberOfCells, const double cellCoverage = 1.0, const double borderCoverage = 0.0) {
        std::vector<int> cellOfVertex(graph.numVertices(), -1);
        for (size_t i = 0; i < numberOfCells; i++) {
            cellOfVertex[centers[i]] = i;
        }
        const size_t coverSize = cellOfVertex.size() * cellCoverage;
        if (Debug) std::cout << "cellCoverage = " << cellCoverage << std::endl;
        if (Debug) std::cout << "coverSize = " << coverSize << std::endl;
        size_t coverCount = 0;
        dijkstra.run(centers, noVertex, [&](const Vertex vertex) {
            if (cellOfVertex[vertex] == -1) {
                cellOfVertex[vertex] = cellOfVertex[dijkstra.getParent(vertex)];
            }
            coverCount++;
        }, [&]() {
            return (coverCount > coverSize);
        });
        if (Debug) std::cout << "CoverCount = " << coverCount << std::endl;
        int minCellRadius = intMax;
        std::vector<Vertex> borderVertices;
        for (const Vertex from : graph.vertices()) {
            for (const Edge edge : graph.edgesFrom(from)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (cellOfVertex[from] == cellOfVertex[to]) continue;
                if (dijkstra.visited(from)) {
                    minCellRadius = std::min(minCellRadius, dijkstra.getDistance(from));
                }
                if (dijkstra.visited(to)) {
                    minCellRadius = std::min(minCellRadius, dijkstra.getDistance(to));
                }
                borderVertices.emplace_back(from);
                borderVertices.emplace_back(to);
            }
        }
        if (Debug) std::cout << "MinCellRadius = " << minCellRadius << std::endl;
        if (Debug) std::cout << "BorderVertices = " << borderVertices.size() << std::endl;
        dijkstra.run(borderVertices, noVertex, [&](const Vertex vertex) {
            cellOfVertex[vertex] = -1;
            coverCount--;
        }, [&]() {
            return dijkstra.getDistance(dijkstra.getQFront()) > minCellRadius * borderCoverage;
        });
        if (Debug) std::cout << "CoverCount = " << coverCount << std::endl;
        return VertexPartition(cellOfVertex);
    }

    inline static std::vector<int> getHeuristicStopValues(const RAPTOR::Data& data) noexcept {
        if (Debug) std::cout << "Computing Heuristic Stop Values..." << std::endl;
        const size_t numberOfStops = data.numberOfStops();
        std::vector<int> stopValues(numberOfStops);

        omp_set_num_threads(4);
        #pragma omp parallel
        {
            const int threadId = omp_get_thread_num();
            pinThreadToCoreId(threadId);
            AssertMsg(omp_get_num_threads() == 4, "Number of threads is " << omp_get_num_threads() << ", but should be " << 4 << "!");

            #pragma omp for
            for (size_t stopId = 0; stopId < numberOfStops; stopId++) {
                const StopId stop = StopId(stopId);
                stopValues[stop] = data.maxRouteSpeedTimesDistance(stop) * data.numberOfRoutesContainingStop(stop);
            }
        }

        if (Debug) std::cout << "Done!" << std::endl;
        return stopValues;
    }

private:
    inline StopId bestUncoveredStop(const std::vector<bool>& covered) const noexcept {
        StopId bestStop = StopId(0);
        while (covered[bestStop]) {
            bestStop++;
        }
        for (StopId stop = StopId(bestStop + 1); stop < covered.size(); stop++) {
            if (covered[stop]) continue;
            if (stopValues[stop] <= stopValues[bestStop]) continue;
            bestStop = stop;
        }
        return bestStop;
    }

    inline std::vector<bool> currentCover(const size_t coverSize) noexcept {
        std::vector<bool> covered(stopValues.size(), false);
        size_t coverCount = 0;
        dijkstra.run(centers, noVertex, [&](const Vertex vertex) {
            if (vertex >= covered.size()) return;
            covered[vertex] = true;
            coverCount++;
        }, [&]() {
            return (coverCount > coverSize);
        });
        return covered;
    }

private:
    RAPTOR::TransferGraph graphCopy;
    const RAPTOR::TransferGraph& graph;
    Dijkstra<RAPTOR::TransferGraph, false> dijkstra;
    std::vector<int> stopValuesCopy;

public:
    const std::vector<int>& stopValues;
    std::vector<StopId> centers;

};
