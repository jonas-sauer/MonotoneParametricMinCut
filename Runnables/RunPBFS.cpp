#include <iostream>
#include <fstream>

#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlow/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/ParametricMinCut/ParametricIBFS.h"

inline void usage() noexcept {
    std::cout << "Runs PBFS to compute a parametric min-cut for the given graph. Arguments:" << std::endl;
    std::cout << "\t-i:   Parametric min-cut instance in DIMACS format." << std::endl;
    std::cout << "\t-o:   Output CSV file to which the solution is written." << std::endl;
    std::cout << "\t-inf: All edge weights with this value or higher are set to infinity (default: intMax/2)." << std::endl;
}

int main(int argc, char **argv) {
    CommandLineParser clp(argc, argv);

    const auto inputFileName = clp.value<std::string>("i");
    const auto outputFileName = clp.value<std::string>("o");
    const auto infinity = clp.value<int>("inf", INFTY);

    if (inputFileName.empty() || outputFileName.empty()) {
        usage();
        return EXIT_SUCCESS;
    }

    ParametricMaxFlowInstance<LinearFlowFunction> instance;
    instance.fromDimacs(inputFileName, infinity);
    Graph::printInfo(instance.graph);
    instance.graph.printAnalysis();
    std::ofstream out(outputFileName);

    ParametricIBFS<LinearFlowFunction, true> algo(instance);
    algo.run();
    const std::vector<double> solution = algo.getVertexBreakpoints();

    out << "vertex,breakpoint" << std::endl;
    for (size_t i = 0; i < solution.size(); i++) {
        out << i << "," << solution[i] << std::endl;
    }

    return EXIT_SUCCESS;
}