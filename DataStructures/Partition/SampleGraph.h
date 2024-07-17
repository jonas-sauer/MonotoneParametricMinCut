#pragma once

#include "NestedDissection.h"

#include "../RAPTOR/Data.h"
#include "../MaxFlowMinCut/FlowGraphs.h"

#include "../../Algorithms/RAPTOR/DijkstraInitializedRAPTOR.h"
#include "../../Algorithms/RAPTOR/Profiler.h"
#include "../../Helpers/Console/Progress.h"
#include "../../Helpers/MultiThreading.h"

namespace Graph {

template<typename RAPTOR_TYPE>
inline DynamicFlowGraph generateSampleGraph(const ThreadPinning& threadPinning, const RAPTOR::Data& data, const NestedDissection& nd, const size_t level, const size_t numberOfSamples, const bool weighted, const std::vector<int>& sampleTimes = {6 * 60 * 60, 14 * 60 * 60, 22 * 60 * 60}) noexcept {
    DynamicFlowGraph result;
    result.addVertices(data.transferGraph.numVertices());
    std::vector<Vertex> separator = nd.getSeparatorOfLevel(level);
    std::cout << "Separator size: " << String::prettyInt(separator.size()) << std::endl;
    std::vector<double> separatorDistances(separator.size(), intMax);
    std::vector<Vertex> sampleVertices(1, separator[0]);
    while (sampleVertices.size() < numberOfSamples) {
        size_t bestIndex = 0;
        for (size_t i = 0; i < separator.size(); i++) {
            const double dist = Geometry::geoDistanceInCM(data.transferGraph.get(Coordinates, sampleVertices.back()), data.transferGraph.get(Coordinates, separator[i]));
            if (separatorDistances[i] > dist) {
                separatorDistances[i] = dist;
            }
            if (separatorDistances[bestIndex] < separatorDistances[i]) {
                bestIndex = i;
            }
        }
        Assert(!Vector::contains(sampleVertices, separator[bestIndex]), "Sample Vertex is a duplicate!");
        sampleVertices.emplace_back(separator[bestIndex]);
    }
    Progress progress(sampleVertices.size() * sampleTimes.size());
    omp_set_num_threads(threadPinning.numberOfThreads);
    #pragma omp parallel
    {
        threadPinning.pinThread();

        RAPTOR_TYPE raptor(data);
        DynamicFlowGraph localResult = result;

        #pragma omp for
        for (size_t i = 0; i < sampleVertices.size(); i++) {
            const Vertex vertex = sampleVertices[i];
            for (const int time : sampleTimes) {
                raptor.run(vertex, time);
                for (const typename RAPTOR_TYPE::ParentPointer& parentPointer : raptor.getParentPointers()) {
                    const Edge edge = localResult.findOrAddEdge(parentPointer.from, parentPointer.to);
                    if (weighted) {
                        localResult.set(Capacity, edge, localResult.get(Capacity, edge) + 1);
                    } else {
                        localResult.set(Capacity, edge, 1);
                    }
                }
                progress++;
            }
        }

        #pragma omp critical
        {
            for (const Vertex from : localResult.vertices()) {
                for (const Edge localEdge : localResult.edgesFrom(from)) {
                    const Vertex to = localResult.get(ToVertex, localEdge);
                    const Edge edge = result.findOrAddEdge(from, to);
                    if (weighted) {
                        result.set(Capacity, edge, result.get(Capacity, edge) + localResult.get(Capacity, edge));
                    } else {
                        result.set(Capacity, edge, 1);
                    }
                }
            }
        }
    }
    std::cout << std::endl;
    return result;
}

inline DynamicFlowGraph generateSampleGraph(const ThreadPinning& threadPinning, RAPTOR::Data& data, const NestedDissection& nd, const size_t level, const size_t numberOfSamples, const bool weighted, const bool useMinTransferTimes, const std::vector<int>& sampleTimes = {6 * 60 * 60, 14 * 60 * 60, 22 * 60 * 60}) noexcept {
    if (useMinTransferTimes) {
        return generateSampleGraph<RAPTOR::DijkstraInitializedRAPTOR<false, RAPTOR::NoProfiler, true, true>>(threadPinning, data, nd, level, numberOfSamples, weighted, sampleTimes);
    } else {
        data.useImplicitDepartureBufferTimes();
        return generateSampleGraph<RAPTOR::DijkstraInitializedRAPTOR<false, RAPTOR::NoProfiler, true, false>>(threadPinning, data, nd, level, numberOfSamples, weighted, sampleTimes);
    }
}

}
