#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"
#include "../Helpers/Console/CommandLineParser.h"

using StaticInstance = StaticMaxFlowInstance<int>;
using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;

inline void loadStaticFromDimacs(const std::string& inputFile, const std::string& outputFile, const int infinity = INFTY) noexcept {
    StaticInstance instance;
    instance.fromDimacs(inputFile, infinity);
    Graph::printInfo(instance.graph);
    instance.graph.printAnalysis();
    instance.serialize(outputFile);
}

inline void staticToParametric(const std::string& inputFile, const std::string& outputFile, const double sourceProbability = 1.0, const double sinkProbability = 0.0) noexcept {
    const StaticInstance staticInstance(inputFile);
    ParametricInstance instance(staticInstance, sourceProbability, sinkProbability);
    Graph::printInfo(instance.graph);
    instance.graph.printAnalysis();
    instance.serialize(outputFile);
}

inline void loadParametricFromDimacs(const std::string& inputFile, const std::string& outputFile, const int infinity = INFTY) noexcept {
    ParametricInstance instance;
    instance.fromDimacs(inputFile, infinity);
    Graph::printInfo(instance.graph);
    instance.graph.printAnalysis();
    instance.serialize(outputFile);
}

inline void usage() noexcept {
    std::cout << "Parse static or parametric max-flow instances. Arguments:" << std::endl;
    std::cout << "\t-m:   Mode:" << std::endl;
    std::cout << "\t\tstatic:     Converts static max-flow instance in DIMACS format to binary format." << std::endl;
    std::cout << "\t\tparametric: Converts parametric max-flow instance in DIMACS format to binary format." << std::endl;
    std::cout << "\t\tconvert:    Converts static to parametric max-flow instance, both in binary format." << std::endl;
    std::cout << "\t\t            By default, edges retain their static weights. For source- and sink-incident edges," << std::endl;
    std::cout << "\t\t            the parameters -p_src and -p_snk control the probability that the static weight" << std::endl;
    std::cout << "\t\t            is replaced with a random parametric one. In this case, both coefficients are" << std::endl;
    std::cout << "\t\t            chosen uniformly at random from [1, max_weight], where max_weight is the highest" << std::endl;
    std::cout << "\t\t            edge weight in the input. Both parameters default to 0." << std::endl;
    std::cout << "\t-i:   Input file." << std::endl;
    std::cout << "\t-o:   Output file." << std::endl;
    std::cout << "\t-inf: All edge weights with this value or higher are set to infinity (default: intMax/2)." << std::endl;
}

int main(int argc, char** argv) {
    CommandLineParser clp(argc, argv);
    const auto inputFileName = clp.value<std::string>("i");
    const auto outputFileName = clp.value<std::string>("o");
    const auto infinity = clp.value<int>("inf", INFTY);
    const auto mode = clp.value<std::string>("m");

    if (inputFileName.empty() || outputFileName.empty()) {
        usage();
        return EXIT_SUCCESS;
    }

    if (mode == "static") {
        loadStaticFromDimacs(inputFileName, outputFileName, infinity);
    } else if (mode == "convert") {
        const auto sourceProbability = clp.value<double>("p_src", 0.0);
        const auto sinkProbability = clp.value<double>("p_snk", 0.0);
        staticToParametric(inputFileName, outputFileName, sourceProbability, sinkProbability);
    } else if (mode == "parametric") {
        loadParametricFromDimacs(inputFileName, outputFileName, infinity);
    } else {
        usage();
    }

    return EXIT_SUCCESS;
}
