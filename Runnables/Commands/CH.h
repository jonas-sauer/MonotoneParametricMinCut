#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "../../Helpers/Debug.h"
#include "../../Helpers/MultiThreading.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"

#include "../../Algorithms/CH/Preprocessing/CHBuilder.h"
#include "../../Algorithms/CH/Preprocessing/MinLevelKey.h"
#include "../../Algorithms/CH/Preprocessing/BlockingKey.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"
#include "../../Algorithms/CH/Query/BidirectionalRPHAST.h"
#include "../../Algorithms/CH/Query/CHQuery.h"
#include "../../Algorithms/CH/Query/BucketQuery.h"
#include "../../Shell/Shell.h"
using namespace Shell;

inline constexpr int ShortcutWeight = 1024;
inline constexpr int DegreeWeight = 0;
inline constexpr int UnidirectionalPopLimit = 500;
inline constexpr int BidirectionalPopLimit = 200;

template<typename PROFILER>
using UnidirectionalWitnessSearch = CH::WitnessSearch<CHCoreGraph, PROFILER, UnidirectionalPopLimit>;
template<typename PROFILER>
using BidirectionalWitnessSearch = CH::BidirectionalWitnessSearch<CHCoreGraph, PROFILER, BidirectionalPopLimit>;

template<typename WITNESS_SEARCH>
using GreedyKey = CH::GreedyKey<WITNESS_SEARCH>;
template<typename WITNESS_SEARCH>
using OrderKey = CH::OrderKey<WITNESS_SEARCH>;
template<typename WITNESS_SEARCH>
using FactorKey = CH::FactorKey<WITNESS_SEARCH, GreedyKey<WITNESS_SEARCH>>;
template<typename WITNESS_SEARCH>
using MinLevelKey = CH::MinLevelKey<WITNESS_SEARCH, GreedyKey<WITNESS_SEARCH>>;
template<typename WITNESS_SEARCH>
using PartialKey = CH::PartialKey<WITNESS_SEARCH, GreedyKey<WITNESS_SEARCH>>;
template<typename WITNESS_SEARCH>
using BlockingKey = CH::BlockingKey<WITNESS_SEARCH, GreedyKey<WITNESS_SEARCH>>;
using StopCriterion = CH::NoStopCriterion;

template<typename CH_BUILDER>
inline CH::CH finalizeCH(CH_BUILDER&& chBuilder, const std::string& orderOutputFile, const std::string& chOutputFile) noexcept {
    chBuilder.copyCoreToCH();
    Order order;
    for (const Vertex vertex : chBuilder.getOrder()) {
        order.emplace_back(vertex);
    }
    order.serialize(orderOutputFile);
    std::cout << "Obtaining CH" << std::endl;
    CH::CH ch(std::move(chBuilder));
    ch.writeBinary(chOutputFile);
    std::cout << std::endl;
    return ch;
}

template<typename PROFILER, typename WITNESS_SEARCH, typename GRAPH, typename KEY_FUNCTION, typename STOP_CRITERION = StopCriterion>
inline CH::CH buildCH(GRAPH& originalGraph, const std::string& orderOutputFile, const std::string& chOutputFile, const KEY_FUNCTION& keyFunction, const STOP_CRITERION& stopCriterion = StopCriterion()) noexcept {
    TravelTimeGraph graph;
    Graph::copy(originalGraph, graph);
    Graph::printInfo(graph);
    CH::Builder<PROFILER, WITNESS_SEARCH, KEY_FUNCTION, STOP_CRITERION, false, false> chBuilder(std::move(graph), graph[TravelTime], keyFunction, stopCriterion);
    chBuilder.run();
    return finalizeCH(chBuilder, orderOutputFile, chOutputFile);
}

template<typename PROFILER, typename WITNESS_SEARCH, typename KEY_FUNCTION, typename STOP_CRITERION = StopCriterion>
inline CH::CH resumeCH(CH::Data&& data, const std::string& orderOutputFile, const std::string& chOutputFile, const KEY_FUNCTION& keyFunction, const STOP_CRITERION& stopCriterion = StopCriterion()) noexcept {
    CH::Builder<PROFILER, WITNESS_SEARCH, KEY_FUNCTION, STOP_CRITERION, false, false> chBuilder(std::move(data), keyFunction, stopCriterion);
    chBuilder.resume();
    return finalizeCH(chBuilder, orderOutputFile, chOutputFile);
}

class BuildCH : public ParameterizedCommand {

public:
    BuildCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildCH", "Computes a CH with greedy key for the input graph.") {
        addParameter("Graph binary");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        TransferGraph graph(getParameter("Graph binary"));
        GreedyKey<WITNESS_SEARCH> keyFunction(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        buildCH<PROFILER, WITNESS_SEARCH>(graph, getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

class ResumeCH : public ParameterizedCommand {

public:
    ResumeCH(BasicShell& shell) :
        ParameterizedCommand(shell, "resumeCH", "Completes the CH computation for the partially contracted input CH.") {
        addParameter("CH input file");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        const CH::CH ch(getParameter("CH input file"));
        CH::Data data = CH::expandData(ch);
        GreedyKey<WITNESS_SEARCH> keyFunction(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        resumeCH<PROFILER, WITNESS_SEARCH>(std::move(data), getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

class BuildCHFromOrder : public ParameterizedCommand {

public:
    BuildCHFromOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "buildCHFromOrder", "Computes a CH from a given order for the input graph.") {
        addParameter("Graph binary");
        addParameter("Order input file");
        addParameter("CH output file");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        TransferGraph graph(getParameter("Graph binary"));
        const std::string orderFile = getParameter("Order input file");
        const Order order(orderFile);
        OrderKey<WITNESS_SEARCH> keyFunction(order);
        buildCH<PROFILER, WITNESS_SEARCH>(graph, orderFile, getParameter("CH output file"), keyFunction);
    }
};

class BuildCoreCH : public ParameterizedCommand {

public:
    BuildCoreCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildCoreCH", "Computes a core-CH for the input network, where all stops are kept uncontracted.") {
        addParameter("Network input file");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Network output file");
        addParameter("Max core degree", "14");
        addParameter("Network type", "raptor", {"intermediate", "csa", "raptor"});
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            return chooseWitnessSearch<CH::FullProfiler>();
        } else {
            return chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() noexcept {
        const std::string witnessSearchType = getParameter("Witness search type");
        if (witnessSearchType == "normal") {
            chooseNetworkType<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            chooseNetworkType<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void chooseNetworkType() noexcept {
        const std::string networkType = getParameter("Network type");
        if (networkType == "raptor") {
            build<RAPTOR::Data, PROFILER, WITNESS_SEARCH>();
        } else if (networkType == "csa") {
            build<CSA::Data, PROFILER, WITNESS_SEARCH>();
        } else {
            build<Intermediate::Data, PROFILER, WITNESS_SEARCH>();
        }
    }

    template<typename NETWORK_TYPE, typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        NETWORK_TYPE data(getParameter("Network input file"));
        data.printInfo();

        std::vector<bool> contractable(data.numberOfStops(), false);
        contractable.resize(data.transferGraph.numVertices(), true);

        const double maxCoreDegree = getParameter<double>("Max core degree");
        std::cout << "Min. core size: " << String::prettyInt(data.numberOfStops()) << std::endl;
        std::cout << "Max. core degree: " << String::prettyInt(maxCoreDegree) << std::endl;
        GreedyKey<WITNESS_SEARCH> greedyKey(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        PartialKey<WITNESS_SEARCH> keyFunction(contractable, data.transferGraph.numVertices(), greedyKey);
        CH::CoreCriterion stopCriterion(data.numberOfStops(), maxCoreDegree);
        const CH::CH ch = buildCH<PROFILER, WITNESS_SEARCH>(data.transferGraph, getParameter("Order output file"), getParameter("CH output file"), keyFunction, stopCriterion);

        Intermediate::TransferGraph resultGraph;
        resultGraph.addVertices(data.transferGraph.numVertices());
        resultGraph[Coordinates] = data.transferGraph[Coordinates];
        for (const Vertex vertex : resultGraph.vertices()) {
            if (ch.isCoreVertex(vertex)) {
                for (const Edge edge : ch.forward.edgesFrom(vertex)) {
                    resultGraph.addEdge(vertex, ch.forward.get(ToVertex, edge)).set(TravelTime, ch.forward.get(Weight, edge));
                }
            }
        }
        Graph::move(std::move(resultGraph), data.transferGraph);
        data.serialize(getParameter("Network output file"));
    }
};

class BuildFactorKeyCH : public ParameterizedCommand {

public:
    BuildFactorKeyCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildFactorKeyCH", "Computes a CH with factor key for the input graph.") {
        addParameter("RAPTOR binary");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Factor", "3");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        std::vector<float> factorOfVertex(data.numberOfStops(), getParameter<float>("Factor"));
        factorOfVertex.resize(data.transferGraph.numVertices(), 1.0);
        GreedyKey<WITNESS_SEARCH> greedyKey(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        FactorKey<WITNESS_SEARCH> keyFunction(factorOfVertex);
        buildCH<PROFILER, WITNESS_SEARCH>(data.transferGraph, getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

class BuildMinLevelKeyCH : public ParameterizedCommand {

public:
    BuildMinLevelKeyCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildMinLevelKeyCH", "Computes a CH with min level key for the input graph.") {
        addParameter("RAPTOR binary");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Min level", "30");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        std::vector<uint16_t> minLevelOfVertex(data.numberOfStops(), getParameter<uint16_t>("Min level"));
        minLevelOfVertex.resize(data.transferGraph.numVertices(), 0);
        GreedyKey<WITNESS_SEARCH> greedyKey(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        MinLevelKey<WITNESS_SEARCH> keyFunction(minLevelOfVertex, greedyKey);
        buildCH<PROFILER, WITNESS_SEARCH>(data.transferGraph, getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

class BuildPartialKeyCH : public ParameterizedCommand {

public:
    BuildPartialKeyCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildPartialKeyCH", "Computes a CH with partial key for the input graph.") {
        addParameter("RAPTOR binary");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Partial factor", "2.0");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        const size_t numStops = data.numberOfStops();
        const size_t numVertices = data.transferGraph.numVertices();
        std::vector<bool> contractable(numStops, false);
        contractable.resize(numVertices, true);
        const size_t minOrderIndex = numVertices - (numStops * getParameter<double>("Partial factor"));
        GreedyKey<WITNESS_SEARCH> greedyKey(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        PartialKey<WITNESS_SEARCH> keyFunction(contractable, minOrderIndex, greedyKey);
        buildCH<PROFILER, WITNESS_SEARCH>(data.transferGraph, getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

class BuildBlockingKeyCH : public ParameterizedCommand {

public:
    BuildBlockingKeyCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildBlockingKeyCH", "Computes a CH with blocking key for the input graph.") {
        addParameter("Graph binary");
        addParameter("Order output file");
        addParameter("CH output file");
        addParameter("Blocking factor", "0.6");
        addParameter("Final core", "1000");
        addParameter("Use full profiler?", "true");
        addParameter("Witness search type", "bidirectional", { "normal", "bidirectional" });
        addParameter("Level weight", "256");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Use full profiler?")) {
            chooseWitnessSearch<CH::FullProfiler>();
        } else {
            chooseWitnessSearch<CH::TimeProfiler>();
        }
    }

private:
    template<typename PROFILER>
    inline void chooseWitnessSearch() const noexcept {
        if (getParameter("Witness search type") == "normal") {
            build<PROFILER, UnidirectionalWitnessSearch<PROFILER>>();
        } else {
            build<PROFILER, BidirectionalWitnessSearch<PROFILER>>();
        }
    }

    template<typename PROFILER, typename WITNESS_SEARCH>
    inline void build() const noexcept {
        TransferGraph graph(getParameter("Graph binary"));
        GreedyKey<WITNESS_SEARCH> greedyKey(ShortcutWeight, getParameter<int>("Level weight"), DegreeWeight);
        BlockingKey<WITNESS_SEARCH> keyFunction(getParameter<double>("Blocking factor"), getParameter<size_t>("Final core"), greedyKey);
        buildCH<PROFILER, WITNESS_SEARCH>(graph, getParameter("Order output file"), getParameter("CH output file"), keyFunction);
    }
};

struct Query {
    Query(const Vertex source, const Vertex target) :
        source(source),
        target(target) {
    }

    Vertex source;
    Vertex target;
};

inline std::vector<Query> generateRandomQueries(const size_t numVertices, const size_t numQueries) noexcept {
    std::mt19937 randomGenerator(42);
    std::uniform_int_distribution<> vertexDistribution(0, numVertices - 1);
    std::vector<Query> queries;
    for (size_t i = 0; i < numQueries; i++) {
        queries.emplace_back(Vertex(vertexDistribution(randomGenerator)), Vertex(vertexDistribution(randomGenerator)));
    }
    return queries;
}

class RunCHQueries : public ParameterizedCommand {

public:
    RunCHQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "runCHQueries", "Runs random CH queries.") {
        addParameter("RAPTOR binary");
        addParameter("CH file");
        addParameter("Number of test queries", "0");
        addParameter("Collect POIs?", "true");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Collect POIs?")) {
            runQueries<true>();
        } else {
            runQueries<false>();
        }
    }

private:
    template<bool COLLECT_POIS>
    inline void runQueries() const noexcept {
        const size_t numberOfQueries = getParameter<size_t>("Number of test queries");
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        data.printInfo();
        CH::CH ch(getParameter("CH file"));
        std::cout << "Starting " << numberOfQueries << " random CH queries" << std::endl;
        const std::vector<Query> queries = generateRandomQueries(ch.numVertices(), numberOfQueries);
        CH::Query<CHGraph, true, false, COLLECT_POIS> algorithm(ch, FORWARD, data.numberOfStops());
        u_int64_t result = 0;
        size_t poiCount = 0;
        Timer timer;
        for (const Query& query : queries) {
            algorithm.run(query.source, query.target);
            result += algorithm.getDistance();
            if constexpr (COLLECT_POIS) {
                poiCount += algorithm.getForwardPOIs().size();
                poiCount += algorithm.getBackwardPOIs().size();
            } else {
                suppressUnusedParameterWarning(poiCount);
            }
        }
        const double time = timer.elapsedMilliseconds();
        std::cout << "Executed " << numberOfQueries << " random queries in " << String::msToString(time) << " (checksum = " << result;
        if constexpr (COLLECT_POIS) {
            std::cout << ", poiCount: " << poiCount;
        }
        std::cout << ")" << std::endl;

    }
};

class RunBucketCHQueries : public ParameterizedCommand {

public:
    RunBucketCHQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "runBucketCHQueries", "Runs random Bucket-CH queries.") {
        addParameter("RAPTOR binary");
        addParameter("CH file");
        addParameter("Number of test queries", "0");
        addParameter("Build only?", "false");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Build only?")) {
            runQueries<true>();
        } else {
            runQueries<false>();
        }
    }

    template<bool BUILD_ONLY>
    inline void runQueries() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        data.printInfo();
        CH::CH ch(getParameter("CH file"));
        CH::BucketQuery<CHGraph, true, BUILD_ONLY> algorithm(ch, FORWARD, data.numberOfStops());

        if constexpr (!BUILD_ONLY) {
            const size_t numberOfQueries = getParameter<size_t>("Number of test queries");
            std::cout << "Starting " << numberOfQueries << " random CH queries" << std::endl;
            const std::vector<Query> queries = generateRandomQueries(ch.numVertices(), numberOfQueries);
            u_int64_t result = 0;
            size_t poiCount = 0;
            Timer timer;
            for (const Query& query : queries) {
                algorithm.run(query.source, query.target);
                result += algorithm.getDistance();
                poiCount += algorithm.getForwardPOIs().size();
                poiCount += algorithm.getBackwardPOIs().size();
            }
            const double time = timer.elapsedMilliseconds();
            std::cout << "Executed " << numberOfQueries << " random queries in " << String::msToString(time) << " (checksum = " << result << ", poiCount: " << poiCount << ")" << std::endl;
        }
    }
};

class RunBidirectionalRPHASTQueries : public ParameterizedCommand {

public:
    RunBidirectionalRPHASTQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "runBidirectionalRPHASTQueries", "Runs random bidirectional PHAST queries.") {
        addParameter("RAPTOR binary");
        addParameter("CH file");
        addParameter("Number of test queries", "0");
        addParameter("Build only?", "false");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Build only?")) {
            runQueries<true>();
        } else {
            runQueries<false>();
        }
    }

    template<bool BUILD_ONLY>
    inline void runQueries() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        data.printInfo();
        CH::CH ch(getParameter("CH file"));
        CH::BidirectionalRPHAST<CHGraph, true, BUILD_ONLY> algorithm(ch, FORWARD, data.numberOfStops());

        if constexpr (!BUILD_ONLY) {
            const size_t numberOfQueries = getParameter<size_t>("Number of test queries");
            std::cout << "Starting " << numberOfQueries << " random CH queries" << std::endl;
            const std::vector<Query> queries = generateRandomQueries(ch.numVertices(), numberOfQueries);
            u_int64_t result = 0;
            size_t poiCount = 0;
            Timer timer;
            for (const Query& query : queries) {
                algorithm.run(query.source, query.target);
                result += algorithm.getDistance();
                poiCount += algorithm.getForwardPOIs().size();
                poiCount += algorithm.getBackwardPOIs().size();
            }
            const double time = timer.elapsedMilliseconds();
            std::cout << "Executed " << numberOfQueries << " random queries in " << String::msToString(time) << " (checksum = " << result << ", poiCount: " << poiCount << ")" << std::endl;
        }
    }
};
