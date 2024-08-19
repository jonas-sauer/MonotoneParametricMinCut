#include "../Algorithms/IBFS.h"
#include "../Algorithms/ParametricIBFS.h"
#include "../Algorithms/ChordScheme.h"
#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlow/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"
#include "../Helpers/Console/Progress.h"

using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ChordSchemeWrapper = ChordSchemeMaxFlowWrapper<pmf::linearFlowFunction>;
using PBFSAlgo = ParametricIBFS<pmf::linearFlowFunction, true>;
using ChordAlgo = ChordScheme<pmf::linearFlowFunction, IBFS<ChordSchemeWrapper>, true>;

template<typename TRUTH_ALGO, typename COMP_ALGO>
inline void compare(const TRUTH_ALGO& truthAlgo, const COMP_ALGO& compAlgo) noexcept {
    const std::vector<double>& groundTruth = truthAlgo.getBreakpoints();
    Progress progress(groundTruth.size());
    double cumulativeError = 0;
    size_t numErrors = 0;
    for (const double breakpoint : groundTruth) {
        const double actualFlow = truthAlgo.getFlowValue(breakpoint);
        const double resultFlow = compAlgo.getFlowValue(breakpoint);
        progress++;
        if (resultFlow <= actualFlow + 1e-06) continue;
        std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        std::cout << actualFlow << " vs. " << resultFlow << " ( " << resultFlow - actualFlow << ")" << std::endl;
        cumulativeError += (resultFlow - actualFlow)/actualFlow;
        numErrors++;
    }
    progress.finished();
    std::cout << "Errors: " << numErrors << "/" << groundTruth.size() << std::endl;
    std::cout << "Cumulative error: " << cumulativeError << std::endl;
    std::cout << "Average error: " << (numErrors == 0 ? 0 : cumulativeError/static_cast<double>(numErrors)) << std::endl;
    std::cout << "Accuracy: " << cumulativeError/groundTruth.size() << std::endl;
}

inline void runPrecisionExperiment(const std::string& instanceFile) noexcept {
    ParametricInstance instance(instanceFile);
    PBFSAlgo parametricIBFS(instance);
    ChordAlgo chordScheme(instance, 0);
    parametricIBFS.run();
    chordScheme.run();

    const std::vector<double>& parametricBreakpoints = parametricIBFS.getBreakpoints();
    const std::vector<double>& chordBreakpoints = chordScheme.getBreakpoints();
    std::cout << "Parametric IBFS: " << parametricBreakpoints.size() << " breakpoints" << std::endl;
    std::cout << "Chord scheme: " << chordBreakpoints.size() << " breakpoints" << std::endl;

    std::cout << "Evaluate chord scheme:" << std::endl;
    compare(parametricIBFS, chordScheme);
    std::cout << "Evaluate parametric IBFS:" << std::endl;
    compare(chordScheme, parametricIBFS);
}

inline void usage() noexcept {
    std::cout << "Runs an experiment to compare the precision of Parametric IBFS and chord scheme with IBFS. Arguments:" << std::endl;
    std::cout << "\t-i:   Parametric max-flow instance in binary format." << std::endl;
}

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    const auto instanceFile = clp.value<std::string>("i");
    if (instanceFile.empty()) {
        usage();
        return EXIT_SUCCESS;
    }
    runPrecisionExperiment(instanceFile);
    return EXIT_SUCCESS;
}
