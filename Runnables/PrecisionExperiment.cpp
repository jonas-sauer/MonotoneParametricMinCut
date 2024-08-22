#include "../Algorithms/IBFS.h"
#include "../Algorithms/ParametricIBFS.h"
#include "../Algorithms/ChordScheme.h"
#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlow/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"
#include "../Helpers/Console/Progress.h"

using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ChordSchemeWrapper = ChordSchemeMaxFlowWrapper<pmf::linearFlowFunction>;
using PBFSAlgo = ParametricIBFS<pmf::linearFlowFunction, false>;
using ChordAlgo = ChordScheme<pmf::linearFlowFunction, IBFS<ChordSchemeWrapper>, false>;

template<typename TRUTH_ALGO, typename COMP_ALGO>
inline void compare(const TRUTH_ALGO& truthAlgo, const COMP_ALGO& compAlgo, const double tolerance, std::ofstream& out) noexcept {
    const std::vector<double>& groundTruth = truthAlgo.getBreakpoints();
    Progress progress(groundTruth.size());
    double cumulativeError = 0;
    size_t numErrors = 0;
    for (const double breakpoint : groundTruth) {
        const double actualFlow = truthAlgo.getFlowValue(breakpoint);
        const double resultFlow = compAlgo.getFlowValue(breakpoint);
        progress++;
        if (resultFlow <= actualFlow + tolerance) continue;
        cumulativeError += (resultFlow - actualFlow)/actualFlow;
        numErrors++;
    }
    progress.finished();
    const double avgError = numErrors == 0 ? 0 : cumulativeError/static_cast<double>(numErrors);
    const double accuracy = cumulativeError/groundTruth.size();
    out << std::to_string(groundTruth.size()) << "," << std::to_string(numErrors) << "," << std::to_string(cumulativeError)  << "," << std::to_string(avgError) << "," << std::to_string(accuracy) << "\n";
}

inline void runPrecisionExperiment(const std::string& instanceFile, std::ofstream& out, const double tolerance) noexcept {
    ParametricInstance instance(instanceFile);
    PBFSAlgo parametricIBFS(instance);
    ChordAlgo chordScheme(instance, 0);
    parametricIBFS.run();
    chordScheme.run();

    const std::vector<double>& parametricBreakpoints = parametricIBFS.getBreakpoints();
    const std::vector<double>& chordBreakpoints = chordScheme.getBreakpoints();
    std::cout << "Parametric IBFS: " << parametricBreakpoints.size() << " breakpoints" << std::endl;
    std::cout << "Chord scheme: " << chordBreakpoints.size() << " breakpoints" << std::endl;

    std::stringstream toleranceHelper;
    toleranceHelper << tolerance;
    const std::string tolerancePrecise = toleranceHelper.str();

    out << "instance,tolerance,algorithm,groundTruthBreakpoints,errors,cumulativeError,avgError,accuracy\n";
    std::cout << "Evaluate chord scheme:" << std::endl;
    const std::string instanceName = String::split(instanceFile, '/').back();
    out << instanceName << "," << tolerancePrecise << ",chordScheme,";
    compare(parametricIBFS, chordScheme, tolerance, out);
    std::cout << "Evaluate parametric IBFS:" << std::endl;
    out << instanceName << "," << tolerancePrecise << ",parametricIBFS,";
    compare(chordScheme, parametricIBFS, tolerance, out);
}

inline void usage() noexcept {
    std::cout << "Evaluates the precision of the results computed by Parametric IBFS and chord scheme with IBFS. For each "
              << "breakpoint computed by either algorithm, the flow values returned by the two algorithms are compared."
              << "Arguments:" << std::endl;
    std::cout << "\t-i:   Parametric max-flow instance in binary format." << std::endl;
    std::cout << "\t-o:   Output CSV file to which the statistics are written." << std::endl;
    std::cout << "\t-t:   Maximum tolerance for absolute difference between flow values (default: 1e-06)." << std::endl;
}

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    const auto instanceFile = clp.value<std::string>("i");
    const auto outputFile = clp.value<std::string>("o");
    const auto tolerance = clp.value<double>("t", 1e-06);

    if (instanceFile.empty() || outputFile.empty()) {
        usage();
        return EXIT_SUCCESS;
    }
    std::ofstream out(outputFile, std::ios::app);
    runPrecisionExperiment(instanceFile, out, tolerance);
    return EXIT_SUCCESS;
}
