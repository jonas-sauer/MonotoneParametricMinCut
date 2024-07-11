#include <iostream>
#include <fstream>

#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/MaxFlowMinCut/ExcessesIBFS.h"
#include "../Algorithms/MaxFlowMinCut/IBFS.h"
#include "../Algorithms/MaxFlowMinCut/ParametricIBFS.h"
#include "../Algorithms/MaxFlowMinCut/PushRelabel.h"
#include "../Algorithms/MaxFlowMinCut/RestartableIBFS.h"
#include "../Algorithms/MaxFlowMinCut/ChordScheme.h"

using FlowEdgeList = ParametricFlowGraphEdgeList<pmf::linearFlowFunction>;
using FlowGraph = ParametricFlowGraph<pmf::linearFlowFunction>;
using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<pmf::linearFlowFunction>;

/**
 * Runs one benchmark and returns the results as a csv line
 * @param instance The instance file
 * @param algorithm the name of the algorithm
 * @param mode the mode in which the algorithm is to be executed
 * @return
 */
std::string runExperiment(std::string instance, std::string algorithm, std::string mode) {
    ParametricInstance graph;

    // Check if file exists as a binary or if we need to read in the .max file
    std::ifstream file(instance);

    if (file.good()) {
        graph = ParametricInstance(instance);
    } else {
        graph.fromDimacs(instance);
    }

    file.close();

    double runtime;
    uint numBreakpoints;

    if (algorithm == "parametricIBFS") {
        if (mode == "whole") {
            ParametricIBFS<pmf::linearFlowFunction, false> algo(graph);
            Timer timer;
            algo.run();
            runtime = timer.elapsedMicroseconds();
            numBreakpoints = algo.getBreakpoints().size();
        }
        else if (mode == "specific") {
            //TODO actually write the code for specific analysis
            throw std::runtime_error("specific analysis of parametric IBFS has not yet been implemented");
        }
    }
    else {
        throw std::runtime_error("No valid algorithm was selected");
    }

    return algorithm + "," + instance + "," + std::to_string(graph.graph.numVertices()) + "," + std::to_string(graph.graph.numEdges()) + "," + std::to_string(runtime) + "," + std::to_string(numBreakpoints) + "\n";
}

/**
 * Runs a benchmark for parametric max flow.
 * Takes an input file with a graph instance (either in DIMACS format or as a binary) with -i
 * Takes an output file to which the result is appended as an appropriate CSV line with -o
 * Takes an algorithm to be run with -a. Options are [INSERT ALGORITHMS]
 * Takes a mode to be run with -m. Options are 'whole' for measuring time over the whole run and 'specific' for measuring algorithm specific detail
 * Output of specific results should be appended to a CSV file specific for this, as they have unique formatting
 */
int main(int argc, char** argv) {
    CommandLineParser parser(argc, argv);

    std::string inputFileName = parser.value<std::string>("i");
    std::string outputFileName = parser.value<std::string>("o");
    std::string algorithm = parser.value<std::string>("a");
    std::string mode = parser.value<std::string>("m");

    std::string results = runExperiment(inputFileName, algorithm, mode);

    std::ofstream outputFile(outputFileName, std::ios::app);
    outputFile << results;
}