#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../../Helpers/Types.h"
#include "../../../Helpers/Timer.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/Vector/Vector.h"

#include "../../../DataStructures/Container/ExternalKHeap.h"

namespace CH {

template<typename GRAPH, typename PROFILER, typename WITNESS_SEARCH_A, typename WITNESS_SEARCH_B>
class ComparingWitnessSearch {

public:
    using Graph = GRAPH;
    using Profiler = PROFILER;
    using WitnessSearchA = WITNESS_SEARCH_A;
    using WitnessSearchB = WITNESS_SEARCH_B;
    using Type = ComparingWitnessSearch<Graph, Profiler, WitnessSearchA, WitnessSearchB>;

public:
    ComparingWitnessSearch() : errorCount(0) {}

    inline void initialize(const GRAPH* graph, const std::vector<int>* weight, Profiler* profiler) noexcept {
        witnessSearchA.initialize(graph, weight, profiler);
        witnessSearchB.initialize(graph, weight, profiler);
    }

    inline bool shortcutIsNecessary(const Vertex from, const Vertex to, const Vertex via, const int shortcutDistance) noexcept {
        const bool resultA = witnessSearchA.shortcutIsNecessary(from, to, via, shortcutDistance);
        const bool resultB = witnessSearchB.shortcutIsNecessary(from, to, via, shortcutDistance);
        if (resultA != resultB) {
            std::cout << "Error!" << std::endl;
            std::cout << "  from: " << from << std::endl;
            std::cout << "  to:   " << to << std::endl;
            std::cout << "  via:  " << via << std::endl;
            std::cout << "  dist: " << shortcutDistance << std::endl;
            std::cout << "  resultA: " << resultA << std::endl;
            std::cout << "  resultB: " << resultB << std::endl;
            std::cout << "  errorCount: " << errorCount << std::endl;
            witnessSearchA.shortcutIsNecessary(from, to, via, shortcutDistance, true);
            witnessSearchB.shortcutIsNecessary(from, to, via, shortcutDistance, true);
            errorCount++;
        }
        return resultA && resultB;
    }

private:
    WitnessSearchA witnessSearchA;
    WitnessSearchB witnessSearchB;
    int errorCount;

};

}
