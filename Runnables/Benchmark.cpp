#include <iostream>
#include <fstream>
#include <sstream>

#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlow/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/IBFS.h"
#include "../Algorithms/ParametricIBFS.h"
#include "../Algorithms/PushRelabel.h"
#include "../Algorithms/RestartableIBFS.h"
#include "../Algorithms/ChordScheme.h"
#include "../Algorithms/ChordSchemeNoContraction.h"

using FlowEdgeList = ParametricFlowGraphEdgeList<LinearFlowFunction>;
using FlowGraph = ParametricFlowGraph<LinearFlowFunction>;
using ParametricInstance = ParametricMaxFlowInstance<LinearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<LinearFlowFunction>;
using ChordSchemeWrapper = ChordSchemeMaxFlowWrapper<LinearFlowFunction>;

inline void runParametricIBFS(const ParametricInstance& instance, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    ParametricIBFS<LinearFlowFunction, true> algo(instance);
    algo.run();
    const double runtime = timer.elapsedMicroseconds();
    out << headerPrefix + "breakpoints,runtime,iterations,bottlenecks,adoptions,avgDistance,drains,initTime,updateTime,reconnectTime,drainTime\n";
    out << rowPrefix << std::to_string(algo.getBreakpoints().size()) + "," + std::to_string(runtime) + "," +
           std::to_string(algo.getNumIterations()) + "," +
           std::to_string(algo.getNumBottlenecks()) + "," +
           std::to_string(algo.getNumAdoptions()) + "," +
           std::to_string(algo.getAvgDistance()) + "," +
           std::to_string(algo.getNumDrains()) + "," +
           std::to_string(algo.getInitTime()) + "," + std::to_string(algo.getUpdateTime()) + "," +
           std::to_string(algo.getReconnectTime()) + "," + std::to_string(algo.getDrainTime()) + "\n";
}

template<typename ALGO>
inline void runChordScheme(const ParametricInstance& instance, const double epsilon, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    ChordScheme<LinearFlowFunction, ALGO, true> algo(instance, epsilon);
    algo.run();
    const double runtime = timer.elapsedMicroseconds();
    out << headerPrefix + "breakpoints,runtime,contractionTime,flowTime,totalVertices\n";
    out << rowPrefix << std::to_string(algo.getBreakpoints().size()) + "," + std::to_string(runtime) + "," +
           std::to_string(algo.getContractionTime()) + "," + std::to_string(algo.getFlowTime()) + "," +
           std::to_string(algo.getTotalVertices()) + "\n";
}

template<typename ALGO>
inline void runChordSchemeNoContraction(const ParametricInstance& instance, const double epsilon, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    ChordSchemeNoContraction<LinearFlowFunction, ALGO> algo(instance, epsilon);
    algo.run();
    const double runtime = timer.elapsedMicroseconds();
    out << headerPrefix + "breakpoints,runtime\n";
    out << rowPrefix << std::to_string(algo.getBreakpoints().size()) + "," + std::to_string(runtime) + "\n";
}

template<typename ALGO>
inline void runRestartableAlgorithm(const ParametricInstance& instance, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    ParametricIBFS<LinearFlowFunction, false> breakpointGetter(instance);
    breakpointGetter.run();

    Timer timer;
    ParametricWrapper wrapper(instance);
    ALGO algo(wrapper);
    algo.run();

    for (uint i = 1; i < breakpointGetter.getBreakpoints().size(); i++) {
        wrapper.setAlpha(breakpointGetter.getBreakpoints()[i]);
        algo.continueAfterUpdate();
    }

    const double runtime = timer.elapsedMicroseconds();

    out << headerPrefix + "breakpoints,runtime,updateTime,flowTime\n";
    out << rowPrefix << std::to_string(breakpointGetter.getBreakpoints().size()) + "," + std::to_string(runtime) + "," +
           std::to_string(algo.getUpdateTime()) + "," + std::to_string(algo.getFlowTime()) + "\n";
}

inline void usage() noexcept {
    std::cout << "Benchmarks a parametric max-flow algorithm. Arguments:" << std::endl;
    std::cout << "\t-i:   Parametric max-flow instance in binary format." << std::endl;
    std::cout << "\t-o:   Output CSV file to which the statistics are written." << std::endl;
    std::cout << "\t-e:   Chord scheme only: Desired precision (default: 0)" << std::endl;
    std::cout << "\t-a:   Algorithm. Options are:" << std::endl;
    std::cout << "\t\t    parametricIBFS" << std::endl;
    std::cout << "\t\t    chordScheme[IBFS]" << std::endl;
    std::cout << "\t\t    chordScheme[PushRelabel]" << std::endl;
    std::cout << "\t\t    chordSchemeNoContraction[IBFS]" << std::endl;
    std::cout << "\t\t    chordSchemeNoContraction[PushRelabel]" << std::endl;
    std::cout << "\t\t    restartableIBFS" << std::endl;
    std::cout << "\t\t    restartablePushRelabel" << std::endl;
}

int main(int argc, char **argv) {
    CommandLineParser clp(argc, argv);

    const auto inputFileName = clp.value<std::string>("i");
    const auto outputFileName = clp.value<std::string>("o");
    const auto algorithm = clp.value<std::string>("a");
    const auto epsilon = clp.value<double>("e");

    if (inputFileName.empty() || outputFileName.empty()) {
        usage();
        return EXIT_SUCCESS;
    }

    const ParametricInstance instance(inputFileName);
    std::ofstream out(outputFileName, std::ios::app);

    std::stringstream epsilonHelper;
    epsilonHelper << epsilon;
    const std::string epsilonPrecise = epsilonHelper.str();

    const std::string headerPrefix = "algorithm,instance,vertices,edges,epsilon,";
    const std::string rowPrefix = algorithm + "," + String::split(inputFileName, '/').back() + "," + std::to_string(instance.graph.numVertices())
            + "," + std::to_string(instance.graph.numEdges()) + "," + epsilonPrecise + ",";

    if (algorithm == "parametricIBFS") {
        runParametricIBFS(instance, out, headerPrefix, rowPrefix);
    } else if (algorithm == "chordScheme[IBFS]") {
        runChordScheme<IBFS<ChordSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "chordScheme[PushRelabel]") {
        runChordScheme<PushRelabel<ChordSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "chordSchemeNoContraction[IBFS]") {
        runChordSchemeNoContraction<IBFS<ChordSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "chordSchemeNoContraction[PushRelabel]") {
        runChordSchemeNoContraction<PushRelabel<ChordSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "restartableIBFS") {
        runRestartableAlgorithm<RestartableIBFS<ParametricWrapper, true>>(instance, out, headerPrefix, rowPrefix);
    } else if (algorithm == "restartablePushRelabel") {
        runRestartableAlgorithm<PushRelabel<ParametricWrapper, true>>(instance, out, headerPrefix, rowPrefix);
    }
    else {
        usage();
    }

    return EXIT_SUCCESS;
}