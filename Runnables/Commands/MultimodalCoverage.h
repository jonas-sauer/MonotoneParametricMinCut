#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <random>

#include "../../Shell/Shell.h"

#include "../../Helpers/MultiThreading.h"

#include "../../DataStructures/Queries/Queries.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/MultimodalData.h"
#include "../../DataStructures/RAPTOR/Entities/Fuzzy.h"
#include "../../DataStructures/TripBased/Data.h"
#include "../../DataStructures/TripBased/DelayData.h"
#include "../../DataStructures/TripBased/MultimodalData.h"

#include "../../Algorithms/CSA/CSA.h"
#include "../../Algorithms/CSA/DijkstraCSA.h"
#include "../../Algorithms/CSA/DijkstraProfileCSA.h"
#include "../../Algorithms/CSA/HLCSA.h"
#include "../../Algorithms/CSA/OneToAllDijkstraCSA.h"
#include "../../Algorithms/CSA/OneToAllDijkstraProfileCSA.h"
#include "../../Algorithms/CSA/ParetoCSA.h"
#include "../../Algorithms/CSA/ParetoULTRACSA.h"
#include "../../Algorithms/CSA/ParetoUPCSA.h"
#include "../../Algorithms/CSA/ProfileCSA.h"
#include "../../Algorithms/CSA/ULTRACSA.h"
#include "../../Algorithms/CSA/ULTRAProfileCSA.h"
#include "../../Algorithms/CSA/UPCSA.h"
#include "../../Algorithms/HL/HLQuery.h"
#include "../../Algorithms/RAPTOR/AlternatingRAPTOR.h"
#include "../../Algorithms/RAPTOR/Bounded/BoundedMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRAMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/HLRAPTOR.h"
#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/InitialTransfers.h"
#include "../../Algorithms/RAPTOR/McRAPTOR.h"
#include "../../Algorithms/RAPTOR/MultimodalMCR.h"
#include "../../Algorithms/RAPTOR/MultimodalULTRAMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/DijkstraRAPTORModule.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/RangeRAPTOR.h"
#include "../../Algorithms/RAPTOR/RAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"
#include "../../Algorithms/RAPTOR/MCR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMHydRA.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMHydRA.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMTBRAPTOR.h"
#include "../../Algorithms/TripBased/BoundedMcQuery/BoundedMcQuery.h"
#include "../../Algorithms/TripBased/Query/HLQuery.h"
#include "../../Algorithms/TripBased/Query/McQuery.h"
#include "../../Algorithms/TripBased/Query/Query.h"
#include "../../Algorithms/TripBased/Query/TransitiveQuery.h"

using namespace Shell;

inline constexpr int fuzzyX[5] = {60, 1, 300, 300, 300};
inline constexpr double fuzzyY[5] = {0.8, 0.1, 0.8, 0.8, 0.8};
inline constexpr int bucketFuzzyX[5] = {60, 1, 1, 1, 1};
inline constexpr double bucketFuzzyY[5] = {0.8, 0.1, 0.1, 0.1, 0.1};
inline constexpr int bucketFactor = 10;

template<typename LABEL>
inline RAPTOR::FuzzyEvaluation<LABEL> getFuzzyEvaluation(const int x[], const double y[]) noexcept {
    static constexpr int NumberOfCriteria = LABEL::NumberOfCriteria;
    int xValues[NumberOfCriteria];
    double yValues[NumberOfCriteria];
    for (size_t i = 0; i < NumberOfCriteria; i++) {
        xValues[i] = x[i];
        yValues[i] = y[i];
    }
    return RAPTOR::FuzzyEvaluation<LABEL>(xValues, yValues);
}

template<typename LABEL>
inline RAPTOR::FuzzyEvaluation<LABEL> getFuzzyEvaluation() noexcept {
    return getFuzzyEvaluation<LABEL>(fuzzyX, fuzzyY);
}

template<typename LABEL>
inline RAPTOR::FuzzyEvaluation<LABEL> getBucketFuzzyEvaluation() noexcept {
    return getFuzzyEvaluation<LABEL>(bucketFuzzyX, bucketFuzzyY);
}

template<typename LABEL>
inline double computeFuzzyCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels, const RAPTOR::FuzzyEvaluation<LABEL>& fuzzyEvaluation) noexcept {
    if (fullLabels.empty()) return 1;
    double containedScore = 0;
    for (const LABEL& fullLabel : fullLabels) {
        containedScore += fuzzyEvaluation.getMaxLEQSimilarity(boundedLabels, fullLabel);
    }
    return containedScore/fullLabels.size();
}

template<typename LABEL>
inline double computeFuzzyCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels) noexcept {
    return computeFuzzyCoverage(fullLabels, boundedLabels, getFuzzyEvaluation<LABEL>());
}

template<typename LABEL>
inline double computeBucketFuzzyCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels) noexcept {
    return computeFuzzyCoverage(fullLabels, boundedLabels, getBucketFuzzyEvaluation<LABEL>());
}

template<typename LABEL>
inline double computeExactCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& algorithmLabels) noexcept {
    if (fullLabels.empty()) return 1;
    size_t contained = 0;
    for (const LABEL& fullLabel : fullLabels) {
        if (Vector::contains(algorithmLabels, fullLabel)) contained++;
    }
    return contained/static_cast<double>(fullLabels.size());
}

struct BucketCoverage {
    BucketCoverage() : coverage(0), queries(0) {}

    inline void add(const double value) noexcept {
        coverage += value;
        queries++;
    }

    inline double operator()() const noexcept {
        return coverage/queries;
    }

    inline friend std::ostream& operator<<(std::ostream& out, const BucketCoverage& c) noexcept {
        if (c.queries == 0) {
            out << "    nan";
        } else {
            out << String::percent(c());
        }
        return out;
    }

    double coverage;
    size_t queries;
};

class BucketCoverageData {
public:
    BucketCoverageData() :
        data(RAPTOR::NUM_TRANSFER_TIME_BUCKETS, std::vector<BucketCoverage>(RAPTOR::NUM_TRANSFER_TIME_BUCKETS)) {
    }

    template<typename LABEL>
    inline void addExact(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& algorithmLabels) noexcept {
        add(fullLabels, algorithmLabels, [&](const std::vector<LABEL>& a, const std::vector<LABEL>& b) {
            return computeExactCoverage(a, b);
        });
    }

    template<typename LABEL>
    inline void addFuzzy(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& algorithmLabels) noexcept {
        add(fullLabels, algorithmLabels, [&](const std::vector<LABEL>& a, const std::vector<LABEL>& b) {
            return computeBucketFuzzyCoverage(a, b);
        });
    }

    inline friend std::ostream& operator<<(std::ostream& out, const BucketCoverageData& data) noexcept {
        for (size_t i = 0; i < data.data.size(); i++) {
            for (size_t j = 0; j < data.data[i].size(); j++) {
                out << "\t" << i << "," << j << ": " << std::setw(7) << data.data[i][j];
            }
            out << std::endl;
        }
        return out;
    }

private:
    template<typename LABEL, typename FUNCTION>
    inline void add(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& algorithmLabels, const FUNCTION& compute) noexcept {
        const std::vector<std::vector<std::vector<LABEL>>> fullBucketLabels = divideByBucket(fullLabels);
        const std::vector<std::vector<std::vector<LABEL>>> algorithmBucketLabels = divideByBucket(algorithmLabels);
        for (size_t i = 0; i < data.size(); i++) {
            for (size_t j = 0; j < data[i].size(); j++) {
                if (fullBucketLabels[i][j].empty()) continue;
                data[i][j].add(compute(fullBucketLabels[i][j], algorithmBucketLabels[i][j]));
            }
        }
    }

    template<typename LABEL>
    inline static std::vector<std::vector<std::vector<LABEL>>> divideByBucket(const std::vector<LABEL>& labels) noexcept {
        std::vector<std::vector<std::vector<LABEL>>> result;
        for (size_t i = 0; i < RAPTOR::NUM_TRANSFER_TIME_BUCKETS; i++) {
            result.emplace_back(std::vector<std::vector<LABEL>>(RAPTOR::NUM_TRANSFER_TIME_BUCKETS));
        }
        for (const LABEL& label : labels) {
            result[label.transferTime[0]/bucketFactor][label.transferTime[1]/bucketFactor].emplace_back(label);
        }
        return result;
    }

    std::vector<std::vector<BucketCoverage>> data;
};

struct CoverageData {
public:
    CoverageData() {}

    template<typename LABEL>
    CoverageData(const std::vector<std::vector<LABEL>>& fullResults, const std::vector<std::vector<LABEL>>& algorithmResults) {
        for (size_t i = 0; i < fullResults.size(); i++) {
            add(fullResults[i], algorithmResults[i]);
        }
    }

    template<typename LABEL>
    inline void add(const std::vector<LABEL>& fullResults, const std::vector<LABEL>& algorithmResults) noexcept {
        coverages.emplace_back(computeExactCoverage(fullResults, algorithmResults));
        fuzzyCoverages.emplace_back(computeFuzzyCoverage(fullResults, algorithmResults));
        const std::vector<LABEL> fullBucketResults = bucketize(fullResults);
        const std::vector<LABEL> algorithmBucketResults = bucketize(algorithmResults);
        bucketCoverages.emplace_back(computeExactCoverage(fullBucketResults, algorithmBucketResults));
        coveragePerBucket.addExact(fullBucketResults, algorithmBucketResults);
        fuzzyBucketCoverages.emplace_back(computeBucketFuzzyCoverage(fullBucketResults, algorithmBucketResults));
        fuzzyCoveragePerBucket.addFuzzy(fullBucketResults, algorithmBucketResults);
    }

    inline friend std::ostream& operator<<(std::ostream& out, const CoverageData& data) noexcept {
        std::cout << "Full coverage:" << std::endl;
        printCoverage(data.coverages);
        std::cout << "Fuzzy coverage:" << std::endl;
        printCoverage(data.fuzzyCoverages);
        std::cout << "Bucket coverage:" << std::endl;
        printCoverage(data.bucketCoverages);
        std::cout << "Coverage per bucket:" << std::endl;
        std::cout << data.coveragePerBucket;
        std::cout << "Fuzzy bucket coverage:" << std::endl;
        printCoverage(data.fuzzyBucketCoverages);
        std::cout << "Fuzzy coverage per bucket:" << std::endl;
        std::cout << data.fuzzyCoveragePerBucket;
        return out;
    }

private:
    template<typename LABEL>
    inline static std::vector<LABEL> bucketize(const std::vector<LABEL>& labels) noexcept {
        RAPTOR::Bag<LABEL> bag;
        for (const LABEL& label : labels) {
            bag.merge(label.bucketize(bucketFactor));
        }
        return bag.labels;
    }

    inline static void printCoverage(std::vector<double> c) noexcept {
        std::sort(c.begin(), c.end());
        std::cout << "\t1st percentile: " << String::percent(Vector::percentile(c, 0.01)) << std::endl;
        std::cout << "\t5th percentile: " << String::percent(Vector::percentile(c, 0.05)) << std::endl;
        std::cout << "\tMedian: " << String::percent(Vector::median(c)) << std::endl;
        std::cout << "\tMean: " << String::percent(Vector::mean(c)) << std::endl;
    }

    std::vector<double> coverages;
    std::vector<double> fuzzyCoverages;
    std::vector<double> bucketCoverages;
    BucketCoverageData coveragePerBucket;
    std::vector<double> fuzzyBucketCoverages;
    BucketCoverageData fuzzyCoveragePerBucket;
};

template<typename QUERY_TYPE, typename ALGORITHM, typename LABEL>
inline void runFullQueries(const ThreadPinning& threadPinning, const std::vector<QUERY_TYPE>& queries, ALGORITHM& algorithm, std::vector<std::vector<LABEL>>& result) noexcept {
    std::vector<double> numJourneys(queries.size());
    result.resize(queries.size());
    Progress progress(queries.size(), true);
    omp_set_num_threads(threadPinning.numberOfThreads);
    std::cout << "Running queries with " << threadPinning.numberOfThreads << " threads." << std::endl;
    #pragma omp parallel
    {
        threadPinning.pinThread();
        ALGORITHM localAlgorithm = algorithm;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < queries.size(); i++) {
            localAlgorithm.run(queries[i].source, queries[i].departureTime, queries[i].target);
            result[i] = localAlgorithm.getResults();
            numJourneys[i] = result[i].size();
            progress++;
        }
    }
    std::cout << "Avg. journeys: " << String::prettyDouble(Vector::mean(numJourneys)) << std::endl;
}

template<typename QUERY_TYPE, typename ALGORITHM_BUILDER, typename LABEL>
inline void runBoundedQueries(const ThreadPinning& threadPinning, const std::vector<QUERY_TYPE>& queries, const ALGORITHM_BUILDER& algorithmBuilder, std::vector<std::vector<LABEL>>& result, const double arrivalSlack, const double tripSlack) noexcept {
    std::vector<double> numJourneys(queries.size());
    result.resize(queries.size());
    Progress progress(queries.size(), true);
    omp_set_num_threads(threadPinning.numberOfThreads);
    std::cout << "Running queries with " << threadPinning.numberOfThreads << " threads." << std::endl;
    #pragma omp parallel
    {
        threadPinning.pinThread();
        typename ALGORITHM_BUILDER::Algorithm algorithm = algorithmBuilder.getAlgorithm();
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < queries.size(); i++) {
            algorithm.run(queries[i].source, queries[i].departureTime, queries[i].target, arrivalSlack, tripSlack);
            result[i] = algorithm.getResults();
            numJourneys[i] = result[i].size();
            progress++;
        }
    }
    std::cout << "Avg. journeys: " << String::prettyDouble(Vector::mean(numJourneys)) << std::endl;
}

inline size_t getMaxTrips(const std::vector<RAPTOR::ArrivalLabel>& anchorLabels, const double tripSlack) noexcept {
    size_t maxTrips = 0;
    for (const RAPTOR::ArrivalLabel& anchorLabel : anchorLabels) {
        maxTrips = std::max(anchorLabel.numberOfTrips, maxTrips);
    }
    return std::ceil(maxTrips * tripSlack);
}

template<typename LABEL>
inline std::vector<RAPTOR::ArrivalLabel> computeAnchorSet(const std::vector<LABEL>& labels) noexcept {
    RAPTOR::Bag<RAPTOR::ArrivalLabel> bag;
    for (const LABEL& label : labels) {
        bag.merge(RAPTOR::ArrivalLabel(label.arrivalTime, label.numberOfTrips));
    }
    Vector::reverse(bag.labels);
    return bag.labels;
}

template<typename LABEL>
inline std::vector<LABEL> computeRestrictedParetoSet(const std::vector<LABEL>& fullLabels, const std::vector<RAPTOR::ArrivalLabel>& anchorLabels, const int departureTime, const double arrivalSlack, const double tripSlack) noexcept {
    std::vector<LABEL> result;
    for (const LABEL& label : fullLabels) {
        if (label.isWithinSlack(anchorLabels, departureTime, arrivalSlack, tripSlack)) {
            result.emplace_back(label);
        }
    }
    return result;
}

template<typename LABEL>
inline std::vector<LABEL> computeRestrictedParetoSet(const std::vector<LABEL>& fullLabels, const int departureTime, const double arrivalSlack, const double tripSlack) noexcept {
    return computeRestrictedParetoSet(fullLabels, computeAnchorSet(fullLabels), departureTime, arrivalSlack, tripSlack);
}

template<typename QUERY_TYPE, typename ANCHOR_ALGORITHM, typename ALGORITHM, typename LABEL>
inline void runBoundedMCRQueries(const ThreadPinning& threadPinning, const std::vector<QUERY_TYPE>& queries, ANCHOR_ALGORITHM& anchorAlgorithm, ALGORITHM& algorithm, std::vector<std::vector<LABEL>>& result, const double arrivalSlack, const double tripSlack) noexcept {
    std::vector<double> numJourneys(queries.size());
    result.resize(queries.size());
    Progress progress(queries.size(), true);
    omp_set_num_threads(threadPinning.numberOfThreads);
    std::cout << "Running queries with " << threadPinning.numberOfThreads << " threads." << std::endl;
    #pragma omp parallel
    {
        threadPinning.pinThread();
        ANCHOR_ALGORITHM localAnchorAlgorithm = anchorAlgorithm;
        ALGORITHM localAlgorithm = algorithm;
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < queries.size(); i++) {
            localAnchorAlgorithm.run(queries[i].source, queries[i].departureTime, queries[i].target);
            std::vector<RAPTOR::ArrivalLabel> anchorResults = localAnchorAlgorithm.getArrivals();
            Vector::reverse(anchorResults);
            const size_t maxTrips = getMaxTrips(anchorResults, tripSlack);
            localAlgorithm.run(queries[i].source, queries[i].departureTime, queries[i].target, maxTrips);
            result[i] = computeRestrictedParetoSet(localAlgorithm.getResults(), anchorResults, queries[i].departureTime, arrivalSlack, tripSlack);
            numJourneys[i] = result[i].size();
            progress++;
        }
    }
    std::cout << "Avg. journeys: " << String::prettyDouble(Vector::mean(numJourneys)) << std::endl;
}

template<size_t NUM_MODES, typename INITIAL_TRANSFERS>
struct MultimodalInitialTransferData {
    inline static constexpr size_t NumTransferModes = NUM_MODES;
    using InitialTransferType = INITIAL_TRANSFERS;
    using MultimodalInitialTransferType = RAPTOR::MultimodalInitialTransfers<NumTransferModes, InitialTransferType>;

    MultimodalInitialTransferData(const std::vector<CH::CH>& chData, const RAPTOR::TransferGraph& forwardTransitiveGraph, const RAPTOR::TransferGraph& backwardTransitiveGraph, const std::vector<size_t>& modes, const size_t numStops) :
        transitiveInitialTransfers(forwardTransitiveGraph, backwardTransitiveGraph),
        multimodalInitialTransfers(transitiveInitialTransfers, initialTransfers, modes, chData[0].numVertices(), numStops) {
        for (const CH::CH& ch : chData) {
            initialTransfers.emplace_back(ch, FORWARD, numStops);
        }
        AssertMsg(modes.size() == NumTransferModes, "Wrong number of modes");
        AssertMsg(initialTransfers.size() == NumTransferModes, "Wrong number of modes");
    }

    RAPTOR::TransitiveInitialTransfers transitiveInitialTransfers;
    std::vector<InitialTransferType> initialTransfers;
    MultimodalInitialTransferType multimodalInitialTransfers;
};

class ComputeMultimodalParetoSets : public ParameterizedCommand {

public:
    ComputeMultimodalParetoSets(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalParetoSets", "Computes full Pareto sets for the given number of queries in a multimodal network.") {
        addParameter("MCR input file");
        addParameter("Core-CH directory");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Number of threads", "max");
        addParameter("Pin multiplier", "1");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData mcrData(getParameter("MCR input file"));
        mcrData.useImplicitDepartureBufferTimes();
        mcrData.printInfo();
        switch (mcrData.modes.size()) {
            case 2:
                run<2>(mcrData);
                break;
            case 3:
                run<3>(mcrData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    inline size_t getNumberOfThreads() const noexcept {
        if (getParameter("Number of threads") == "max") {
            return numberOfCores();
        } else {
            return getParameter<int>("Number of threads");
        }
    }

    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& mcrData) const noexcept {
        const std::string coreCHDirectory(getParameter("Core-CH directory"));
        std::vector<CH::CH> coreCHData;
        for (const size_t mode : mcrData.modes) {
            coreCHData.emplace_back(coreCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }

        const std::vector<VertexQuery> queries = generateRandomVertexQueries(coreCHData[0].numVertices(), getParameter<size_t>("Number of queries"));
        RAPTOR::MultimodalMCR<true, NUM_MODES, RAPTOR::AggregateProfiler> mcr(mcrData, coreCHData);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> mcrResults;
        const size_t numberOfThreads = getNumberOfThreads();
        const size_t pinMultiplier = getParameter<size_t>("Pin multiplier");
        runFullQueries(ThreadPinning(numberOfThreads, pinMultiplier), queries, mcr, mcrResults);
        IO::serialize(getParameter("Output file"), queries, mcrResults);
    }
};

class ComputeMultimodalRestrictedParetoSets : public ParameterizedCommand {

public:
    ComputeMultimodalRestrictedParetoSets(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalRestrictedParetoSets", "Computes restricted Pareto sets for the given full Pareto sets.") {
        addParameter("Number of modes");
        addParameter("Input file");
        addParameter("Arrival slack");
        addParameter("Trip slack");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        switch (getParameter<size_t>("Number of modes")) {
            case 2:
                run<2>();
                break;
            case 3:
                run<3>();
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    template<size_t NUM_MODES>
    inline void run() const noexcept {
        std::vector<VertexQuery> queries;
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> fullResults;
        IO::deserialize(getParameter("Input file"), queries, fullResults);

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> boundedResults;
        for (size_t i = 0; i < queries.size(); i++) {
            boundedResults.emplace_back(computeRestrictedParetoSet(fullResults[i], queries[i].departureTime, arrivalSlack, tripSlack));
        }
        IO::serialize(getParameter("Output file"), queries, boundedResults, arrivalSlack, tripSlack);
    }
};

class ComputeMultimodalULTRAMcRAPTORCoverage : public ParameterizedCommand {

public:
    ComputeMultimodalULTRAMcRAPTORCoverage(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalULTRAMcRAPTORCoverage", "Computes coverage of multimodal ULTRA-McRAPTOR results compared to MCR.") {
        addParameter("ULTRA-RAPTOR input file");
        addParameter("Bucket-CH directory");
        addParameter("Results input file");
        addParameter("Number of threads", "max");
        addParameter("Pin multiplier", "1");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData ultraRaptorData(getParameter("ULTRA-RAPTOR input file"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        switch (ultraRaptorData.modes.size()) {
            case 2:
                run<2>(ultraRaptorData);
                break;
            case 3:
                run<3>(ultraRaptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    inline size_t getNumberOfThreads() const noexcept {
        if (getParameter("Number of threads") == "max") {
            return numberOfCores();
        } else {
            return getParameter<int>("Number of threads");
        }
    }

    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& ultraRaptorData) const noexcept {
        const std::string bucketCHDirectory(getParameter("Bucket-CH directory"));
        std::vector<CH::CH> bucketCHData;
        for (const size_t mode : ultraRaptorData.modes) {
            bucketCHData.emplace_back(bucketCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }

        std::vector<VertexQuery> queries;
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> mcrResults;
        IO::deserialize(getParameter("Results input file"), queries, mcrResults);

        RAPTOR::MultimodalULTRAMcRAPTOR<NUM_MODES, RAPTOR::AggregateProfiler> ultraMcRaptor(ultraRaptorData, bucketCHData);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> ultraResults;
        const size_t numberOfThreads = getNumberOfThreads();
        const size_t pinMultiplier = getParameter<size_t>("Pin multiplier");
        runFullQueries(ThreadPinning(numberOfThreads, pinMultiplier), queries, ultraMcRaptor, ultraResults);

        const CoverageData coverageData(mcrResults, ultraResults);
        std::cout << coverageData << std::endl;
    }
};

template<size_t NUM_MODES, bool ADD_INTERMEDIATE_OVERHEAD, typename INITIAL_TRANSFERS, typename BACKWARD_GRAPH_TYPE>
struct UBMRAPTORBuilder {
    using Algorithm = RAPTOR::MultimodalUBMRAPTOR<NUM_MODES, RAPTOR::NoProfiler, INITIAL_TRANSFERS, ADD_INTERMEDIATE_OVERHEAD>;

    UBMRAPTORBuilder(const RAPTOR::MultimodalData& ultraRaptorData, const RAPTOR::Data& pruningData, const RAPTOR::Data& reversePruningData, const RAPTOR::TransferGraph& backwardTransitiveGraph, const std::vector<BACKWARD_GRAPH_TYPE>& backwardGraphs) :
        ultraRaptorData(ultraRaptorData),
        pruningData(pruningData),
        reversePruningData(reversePruningData),
        backwardTransitiveGraph(backwardTransitiveGraph),
        backwardGraphs(backwardGraphs) {
    }

    inline Algorithm getAlgorithm() const noexcept {
        return Algorithm(ultraRaptorData, pruningData, reversePruningData, backwardTransitiveGraph, backwardGraphs);
    }

    const RAPTOR::MultimodalData& ultraRaptorData;
    const RAPTOR::Data& pruningData;
    const RAPTOR::Data& reversePruningData;
    const RAPTOR::TransferGraph& backwardTransitiveGraph;
    const std::vector<BACKWARD_GRAPH_TYPE>& backwardGraphs;
};

class ComputeMultimodalUBMRAPTORCoverage : public ParameterizedCommand {

public:
    ComputeMultimodalUBMRAPTORCoverage(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalUBMRAPTORCoverage", "Computes coverage of multimodal UBM-RAPTOR results compared to MCR.") {
        addParameter("ULTRA-RAPTOR input file");
        addParameter("Bucket-CH directory");
        addParameter("Results input file");
        addParameter("Add intermediate overhead?", "false");
        addParameter("Number of threads", "max");
        addParameter("Pin multiplier", "1");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData ultraRaptorData(getParameter("ULTRA-RAPTOR input file"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        switch (ultraRaptorData.modes.size()) {
            case 2:
                chooseAddIntermediateOverhead<2>(ultraRaptorData);
                break;
            case 3:
                chooseAddIntermediateOverhead<3>(ultraRaptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    inline size_t getNumberOfThreads() const noexcept {
        if (getParameter("Number of threads") == "max") {
            return numberOfCores();
        } else {
            return getParameter<int>("Number of threads");
        }
    }

    template<size_t NUM_MODES>
    inline void chooseAddIntermediateOverhead(const RAPTOR::MultimodalData& ultraRaptorData) const noexcept {
        if (getParameter<bool>("Add intermediate overhead?")) {
            run<NUM_MODES, true>(ultraRaptorData);
        } else {
            run<NUM_MODES, false>(ultraRaptorData);
        }
    }

    template<size_t NUM_MODES, bool ADD_INTERMEDIATE_OVERHEAD>
    inline void run(const RAPTOR::MultimodalData& ultraRaptorData) const noexcept {
        const RAPTOR::Data pruningData = ultraRaptorData.getPruningData();
        const RAPTOR::Data reversePruningData = pruningData.reverseNetwork();

        const std::string bucketCHDirectory(getParameter("Bucket-CH directory"));
        std::vector<CH::CH> bucketCHData;
        for (const size_t mode : ultraRaptorData.modes) {
            bucketCHData.emplace_back(bucketCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }

        RAPTOR::TransferGraph backwardTransitiveGraph = ultraRaptorData.raptorData.transferGraph;
        backwardTransitiveGraph.revert();

        std::vector<VertexQuery> queries;
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> boundedResults;
        double arrivalSlack, tripSlack;
        IO::deserialize(getParameter("Results input file"), queries, boundedResults, arrivalSlack, tripSlack);

        UBMRAPTORBuilder<NUM_MODES, ADD_INTERMEDIATE_OVERHEAD, RAPTOR::BucketCHInitialTransfers, CH::CH> builder(ultraRaptorData, pruningData, reversePruningData, backwardTransitiveGraph, bucketCHData);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> ultraResults;
        const size_t numberOfThreads = getNumberOfThreads();
        const size_t pinMultiplier = getParameter<size_t>("Pin multiplier");
        runBoundedQueries(ThreadPinning(numberOfThreads, pinMultiplier), queries, builder, ultraResults, arrivalSlack, tripSlack);

        const CoverageData coverageData(boundedResults, ultraResults);
        std::cout << coverageData << std::endl;
    }
};

class ComputeMultimodalBoundedMcRAPTORCoverage : public ParameterizedCommand {

public:
    ComputeMultimodalBoundedMcRAPTORCoverage(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalBoundedMcRAPTORCoverage", "Computes coverage of multimodal transitive Bounded ULTRA-McRAPTOR results compared to MCR.") {
        addParameter("RAPTOR input file");
        addParameter("MCR input file");
        addParameter("Core-CH directory");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
        addParameter("Number of threads", "max");
        addParameter("Pin multiplier", "1");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData ultraRaptorData(getParameter("RAPTOR input file"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        switch (ultraRaptorData.modes.size()) {
            case 2:
                run<2>(ultraRaptorData);
                break;
            case 3:
                run<3>(ultraRaptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    inline size_t getNumberOfThreads() const noexcept {
        if (getParameter("Number of threads") == "max") {
            return numberOfCores();
        } else {
            return getParameter<int>("Number of threads");
        }
    }

    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& ultraRaptorData) const noexcept {
        std::vector<RAPTOR::TransferGraph> backwardGraphs;
        for (const size_t mode : ultraRaptorData.modes) {
            backwardGraphs.emplace_back(ultraRaptorData.getTransferGraph(mode));
            backwardGraphs.back().revert();
        }
        const RAPTOR::Data pruningData = ultraRaptorData.getPruningData();
        const RAPTOR::Data reversePruningData = pruningData.reverseNetwork();
        RAPTOR::MultimodalData mcrData(getParameter("MCR input file"));
        Ensure(mcrData.modes == ultraRaptorData.modes, "Different transfer modes!");
        mcrData.useImplicitDepartureBufferTimes();
        mcrData.printInfo();
        const std::string coreCHDirectory(getParameter("Core-CH directory"));
        std::vector<CH::CH> coreCHData;
        for (const size_t mode : ultraRaptorData.modes) {
            coreCHData.emplace_back(coreCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }
        RAPTOR::TransferGraph backwardTransitiveGraph = ultraRaptorData.raptorData.transferGraph;
        backwardTransitiveGraph.revert();

        const std::vector<StopQuery> queries = generateRandomStopQueries(pruningData.numberOfStops(), getParameter<size_t>("Number of queries"));
        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::MultimodalMCR<true, NUM_MODES, RAPTOR::AggregateProfiler> mcr(mcrData, coreCHData);
        MultimodalInitialTransferData<NUM_MODES, RAPTOR::CoreCHInitialTransfers> multimodalInitialTransferData(coreCHData, ultraRaptorData.raptorData.transferGraph, backwardTransitiveGraph, ultraRaptorData.modes, ultraRaptorData.raptorData.numberOfStops());
        using InitialTransferType = RAPTOR::MultimodalInitialTransfers<NUM_MODES, RAPTOR::CoreCHInitialTransfers>;
        RAPTOR::DijkstraRAPTOR<InitialTransferType, RAPTOR::AggregateProfiler> dijkstraRaptor(pruningData, multimodalInitialTransferData.multimodalInitialTransfers);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> boundedResults;
        const size_t numberOfThreads = getNumberOfThreads();
        const size_t pinMultiplier = getParameter<size_t>("Pin multiplier");
        const ThreadPinning threadPinning(numberOfThreads, pinMultiplier);
        runBoundedMCRQueries(threadPinning, queries, dijkstraRaptor, mcr, boundedResults, arrivalSlack, tripSlack);

        UBMRAPTORBuilder<NUM_MODES, true, RAPTOR::TransitiveInitialTransfers, RAPTOR::TransferGraph> builder(ultraRaptorData, pruningData, reversePruningData, backwardTransitiveGraph, backwardGraphs);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> ultraResults;
        runBoundedQueries(threadPinning, queries, builder, ultraResults, arrivalSlack, tripSlack);

        const CoverageData coverageData(boundedResults, ultraResults);
        std::cout << coverageData << std::endl;
    }
};

template<size_t NUM_MODES>
struct UBMHydRABuilder {
    using Algorithm = RAPTOR::MultimodalUBMHydRA<NUM_MODES, RAPTOR::NoProfiler>;

    UBMHydRABuilder(const TripBased::MultimodalData& ultraTripBasedData, const TripBased::Data& forwardPruningData, const TripBased::Data& backwardPruningData, const RAPTOR::TransferGraph& backwardTransitiveGraph, const std::vector<CH::CH>& backwardGraphs) :
        ultraTripBasedData(ultraTripBasedData),
        forwardPruningData(forwardPruningData),
        backwardPruningData(backwardPruningData),
        backwardTransitiveGraph(backwardTransitiveGraph),
        backwardGraphs(backwardGraphs) {
    }

    inline Algorithm getAlgorithm() const noexcept {
        return Algorithm(ultraTripBasedData, forwardPruningData, backwardPruningData, backwardTransitiveGraph, backwardGraphs);
    }

    const TripBased::MultimodalData& ultraTripBasedData;
    const TripBased::Data& forwardPruningData;
    const TripBased::Data& backwardPruningData;
    const RAPTOR::TransferGraph& backwardTransitiveGraph;
    const std::vector<CH::CH>& backwardGraphs;
};

class ComputeMultimodalUBMHydRACoverage : public ParameterizedCommand {

public:
    ComputeMultimodalUBMHydRACoverage(BasicShell& shell) :
        ParameterizedCommand(shell, "computeMultimodalUBMHydRACoverage", "Computes coverage of multimodal UBM-HydRA results compared to MCR.") {
        addParameter("ULTRA-Trip-Based input file");
        addParameter("Bounded forward Trip-Based input file");
        addParameter("Bounded backward Trip-Based input file");
        addParameter("Bucket-CH directory");
        addParameter("Results input file");
        addParameter("Number of threads", "max");
        addParameter("Pin multiplier", "1");
    }

    virtual void execute() noexcept {
        const TripBased::MultimodalData ultraTripBasedData(getParameter("ULTRA-Trip-Based input file"));
        ultraTripBasedData.printInfo();
        switch (ultraTripBasedData.modes.size()) {
            case 2:
                run<2>(ultraTripBasedData);
                break;
            case 3:
                run<3>(ultraTripBasedData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

private:
    inline size_t getNumberOfThreads() const noexcept {
        if (getParameter("Number of threads") == "max") {
            return numberOfCores();
        } else {
            return getParameter<int>("Number of threads");
        }
    }

    template<size_t NUM_MODES>
    inline void run(const TripBased::MultimodalData& ultraTripBasedData) const noexcept {
        const TripBased::MultimodalData forwardBoundedData(getParameter("Bounded forward Trip-Based input file"));
        Ensure(forwardBoundedData.modes == ultraTripBasedData.modes, "Different transfer modes!");
        forwardBoundedData.printInfo();
        const TripBased::Data forwardPruningData = forwardBoundedData.getPruningData();
        const TripBased::MultimodalData backwardBoundedData(getParameter("Bounded backward Trip-Based input file"));
        Ensure(backwardBoundedData.modes == ultraTripBasedData.modes, "Different transfer modes!");
        backwardBoundedData.printInfo();
        const TripBased::Data backwardPruningData = backwardBoundedData.getPruningData();
        const std::string bucketCHDirectory(getParameter("Bucket-CH directory"));
        std::vector<CH::CH> bucketCHData;
        for (const size_t mode : ultraTripBasedData.modes) {
            bucketCHData.emplace_back(bucketCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }
        RAPTOR::TransferGraph backwardTransitiveGraph = ultraTripBasedData.tripData.raptorData.transferGraph;
        backwardTransitiveGraph.revert();

        std::vector<VertexQuery> queries;
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> boundedResults;
        double arrivalSlack, tripSlack;
        IO::deserialize(getParameter("Results input file"), queries, boundedResults, arrivalSlack, tripSlack);

        UBMHydRABuilder<NUM_MODES> builder(ultraTripBasedData, forwardPruningData, backwardPruningData, backwardTransitiveGraph, bucketCHData);
        std::vector<std::vector<RAPTOR::MultimodalParetoLabel<NUM_MODES>>> ultraResults;
        const size_t numberOfThreads = getNumberOfThreads();
        const size_t pinMultiplier = getParameter<size_t>("Pin multiplier");
        runBoundedQueries(ThreadPinning(numberOfThreads, pinMultiplier), queries, builder, ultraResults, arrivalSlack, tripSlack);

        const CoverageData coverageData(boundedResults, ultraResults);
        std::cout << coverageData << std::endl;
    }
};
