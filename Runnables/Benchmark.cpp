#include <iostream>
#include <fstream>
#include <sstream>

#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlow/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"

#include "../Algorithms/StaticMaxFlow/IBFS.h"
#include "../Algorithms/ParametricMaxFlow/ParametricIBFS.h"
#include "../Algorithms/StaticMaxFlow/PushRelabel.h"
#include "../Algorithms/ParametricMaxFlow/RestartableIBFS.h"
#include "../Algorithms/ParametricMaxFlow/DichotomicScheme.h"
#include "../Algorithms/ParametricMaxFlow/DichotomicSchemeNoContraction.h"

using FlowEdgeList = ParametricFlowGraphEdgeList<LinearFlowFunction>;
using FlowGraph = ParametricFlowGraph<LinearFlowFunction>;
using ParametricInstance = ParametricMaxFlowInstance<LinearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<LinearFlowFunction>;
using DichotomicSchemeWrapper = DichotomicSchemeMaxFlowWrapper<LinearFlowFunction>;

inline void runParametricIBFS(const ParametricInstance& instance, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    ParametricIBFS<LinearFlowFunction, true> algo(instance);
    algo.run();
    const double runtime = timer.elapsedMicroseconds();
    out << headerPrefix + "breakpoints,runtime,iterations,bottlenecks,adoptions,avgDistance,drains,initTime,updateTime,reconnectTime,drainTime\n";
    out << rowPrefix << std::to_string(algo.getBreakpoints().size()) + "," + std::to_string(runtime) + "," +
           algo.getMeasurementsCSV() + "\n";
}

template<typename ALGO>
inline void runChordScheme(const ParametricInstance& instance, const double epsilon, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    DichotomicScheme<LinearFlowFunction, ALGO, true> algo(instance, epsilon);
    algo.run();
    const double runtime = timer.elapsedMicroseconds();
    out << headerPrefix + "breakpoints,runtime,contractionTime,flowTime,totalVertices\n";
    out << rowPrefix << std::to_string(algo.getBreakpoints().size()) + "," + std::to_string(runtime) + "," +
           algo.getMeasurementsCSV() + "\n";
}

template<typename ALGO>
inline void runChordSchemeNoContraction(const ParametricInstance& instance, const double epsilon, std::ofstream& out, const std::string& headerPrefix, const std::string& rowPrefix) noexcept {
    Timer timer;
    DichotomicSchemeNoContraction<LinearFlowFunction, ALGO> algo(instance, epsilon);
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
           algo.getMeasurementsCSV() + "\n";
}

inline void usage() noexcept {
    std::cout << "Benchmarks a parametric max-flow algorithm. Arguments:" << std::endl;
    std::cout << "\t-i:   Parametric max-flow instance in binary format." << std::endl;
    std::cout << "\t-o:   Output CSV file to which the statistics are written." << std::endl;
    std::cout << "\t-e:   Chord scheme only: Desired precision (default: 0)" << std::endl;
    std::cout << "\t-a:   Algorithm. Options are:" << std::endl;
    std::cout << "\t\t    PBFS" << std::endl;
    std::cout << "\t\t    DS[IBFS]" << std::endl;
    std::cout << "\t\t    DS[PRF]" << std::endl;
    std::cout << "\t\t    DSNoContraction[IBFS]" << std::endl;
    std::cout << "\t\t    DSNoContraction[PRF]" << std::endl;
    std::cout << "\t\t    restartableIBFS" << std::endl;
    std::cout << "\t\t    restartablePRF" << std::endl;
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

    if (algorithm == "PBFS") {
        runParametricIBFS(instance, out, headerPrefix, rowPrefix);
    } else if (algorithm == "DS[IBFS]") {
        runChordScheme<IBFS<DichotomicSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "DS[PRF]") {
        runChordScheme<PushRelabel<DichotomicSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "DSNoContraction[IBFS]") {
        runChordSchemeNoContraction<IBFS<DichotomicSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "DSNoContraction[PRF]") {
        runChordSchemeNoContraction<PushRelabel<DichotomicSchemeWrapper>>(instance, epsilon, out, headerPrefix, rowPrefix);
    } else if (algorithm == "restartableIBFS") {
        runRestartableAlgorithm<RestartableIBFS<ParametricWrapper, true>>(instance, out, headerPrefix, rowPrefix);
    } else if (algorithm == "restartablePRF") {
        runRestartableAlgorithm<PushRelabel<ParametricWrapper, true>>(instance, out, headerPrefix, rowPrefix);
    }
    else {
        usage();
    }

    return EXIT_SUCCESS;
}