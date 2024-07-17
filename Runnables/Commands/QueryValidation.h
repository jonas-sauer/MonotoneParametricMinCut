#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <random>

#include "../../Shell/Shell.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/MultimodalData.h"
#include "../../DataStructures/RAPTOR/Entities/Fuzzy.h"
#include "../../DataStructures/TripBased/Data.h"
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
#include "../../Algorithms/CSA/UPProfileCSA.h"
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
#include "../../Algorithms/RAPTOR/OneToAllDijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/DijkstraRAPTORModule.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/RangeRAPTOR.h"
#include "../../Algorithms/RAPTOR/RAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"
#include "../../Algorithms/RAPTOR/MCR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMHydRA.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMHydRA.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/NewUBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/UBMTBRAPTOR.h"
#include "../../Algorithms/RAPTOR/UPRAPTOR.h"
#include "../../Algorithms/TripBased/BoundedMcQuery/BoundedMcQuery.h"
#include "../../Algorithms/TripBased/Query/HLQuery.h"
#include "../../Algorithms/TripBased/Query/McQuery.h"
#include "../../Algorithms/TripBased/Query/Query.h"
#include "../../Algorithms/TripBased/Query/TransitiveQuery.h"
#include "../../Algorithms/TripBased/Query/UPQuery.h"

using namespace Shell;

inline constexpr int fuzzyX[5] = {60, 1, 300, 300, 300};
inline constexpr double fuzzyY[5] = {0.8, 0.1, 0.8, 0.8, 0.8};
inline constexpr size_t K = 5;

template<typename LABEL>
inline RAPTOR::FuzzyEvaluation<LABEL> getFuzzyEvaluation() noexcept {
    static constexpr int NumberOfCriteria = LABEL::NumberOfCriteria;
    int xValues[NumberOfCriteria];
    double yValues[NumberOfCriteria];
    for (size_t i = 0; i < NumberOfCriteria; i++) {
        xValues[i] = fuzzyX[i];
        yValues[i] = fuzzyY[i];
    }
    return RAPTOR::FuzzyEvaluation<LABEL>(xValues, yValues);
}

template<typename LABEL>
inline double computeTopKFuzzyCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels) noexcept {
    using ScoredLabel = typename RAPTOR::ScoredLabel<LABEL>;
    if (fullLabels.empty()) return 1;
    const RAPTOR::FuzzyEvaluation<LABEL> fuzzyEvaluation = getFuzzyEvaluation<LABEL>();
    const std::vector<ScoredLabel> topFullLabels = fuzzyEvaluation.getTopKLabels(fullLabels, K);
    const std::vector<ScoredLabel> topBoundedLabels = fuzzyEvaluation.getTopKLabels(boundedLabels, K);
    double containedScore = 0;
    double totalScore = 0;
    for (const ScoredLabel& fullLabel : topFullLabels) {
        totalScore += fullLabel.score;
        if (Vector::contains(topBoundedLabels, fullLabel)) {
            containedScore += fullLabel.score;
        }
    }
    return containedScore/totalScore;
}

template<typename LABEL>
inline double computeFullFuzzyCoverage(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels) noexcept {
    using ScoredLabel = typename RAPTOR::ScoredLabel<LABEL>;
    if (fullLabels.empty()) return 1;
    const RAPTOR::FuzzyEvaluation<LABEL> fuzzyEvaluation = getFuzzyEvaluation<LABEL>();
    const std::vector<ScoredLabel> topFullLabels = fuzzyEvaluation.scoreLabels(fullLabels);
    const std::vector<ScoredLabel> topBoundedLabels = fuzzyEvaluation.scoreLabels(boundedLabels);
    double containedScore = 0;
    double totalScore = 0;
    for (const ScoredLabel& fullLabel : topFullLabels) {
        totalScore += fullLabel.score;
        containedScore += fuzzyEvaluation.getMaxSimilarity(fullLabel, topBoundedLabels) * fullLabel.score;
    }
    return containedScore/totalScore;
}

template<typename LABEL>
inline bool compareBoundedResults(const std::vector<LABEL>& fullLabels, const std::vector<LABEL>& boundedLabels, const std::vector<RAPTOR::ArrivalLabel>& anchorLabels, const int departureTime, const double arrivalSlack, const double tripSlack) noexcept {
    for (const LABEL& label : boundedLabels) {
        if (!Vector::contains(fullLabels, label)) {
            std::cout << "Bounded label is not contained in full labels!" << std::endl;
            std::cout << label << std::endl;
            return false;
        }
    }
    for (const LABEL& label : fullLabels) {
        const bool contained = Vector::contains(boundedLabels, label);
        const bool withinSlack = label.isWithinSlack(anchorLabels, departureTime, arrivalSlack, tripSlack);
        if (contained != withinSlack) {
            std::cout << "Full label is " << (contained ? "" : "not ") << "contained but " << (withinSlack ? "" : "not ") << "within slack!" << std::endl;
            std::cout << label << std::endl;
            return false;
        }
    }
    return true;
}

inline bool compareBoundedResults(const std::vector<RAPTOR::WalkingParetoLabel>& fullLabels, const std::vector<RAPTOR::WalkingParetoLabel>& boundedLabels, const std::vector<RAPTOR::WalkingParetoLabel>& anchorLabels, const double arrivalFactor, const double tripsPerHour) noexcept {
    for (const RAPTOR::WalkingParetoLabel& label : boundedLabels) {
        if (!Vector::contains(fullLabels, label)) {
            std::cout << "Bounded label is not contained in full labels!" << std::endl;
            std::cout << label << std::endl;
            return false;
        }
    }
    for (const RAPTOR::WalkingParetoLabel& label : fullLabels) {
        const bool contained = Vector::contains(boundedLabels, label);
        const bool withinSlack = label.isWithinSlack(anchorLabels, arrivalFactor, tripsPerHour);
        if (contained != withinSlack) {
            std::cout << "Full label is " << (contained ? "" : "not ") << "contained but " << (withinSlack ? "" : "not ") << "within slack!" << std::endl;
            std::cout << label << std::endl;
            return false;
        }
    }
    return true;
}

class ValidateOneToAllDijkstraRAPTOR : public ParameterizedCommand {

public:
    ValidateOneToAllDijkstraRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateOneToAllDijkstraRAPTOR", "Validates journeys computed by one-to-all Dijkstra-RAPTOR by comparing to the Range-RAPTOR implementation on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();

        RAPTOR::TransferGraph backwardGraph = raptorData.transferGraph;
        backwardGraph.revert();
        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::AggregateProfiler> raptor(raptorData, raptorData.transferGraph, backwardGraph);
        RAPTOR::RangeRAPTOR::DijkstraRAPTORModule<RAPTOR::DijkstraInitialTransfers, RAPTOR::AggregateProfiler, false, false, false, false> rangeRaptor(raptorData, raptorData.transferGraph, backwardGraph);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % raptorData.transferGraph.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime);
            rangeRaptor.run(source, departureTime);
            for (const Vertex target : raptorData.transferGraph.vertices()) {
                const std::vector<RAPTOR::ArrivalLabel> raptorArrivals = raptor.getArrivals(target);
                const std::vector<RAPTOR::ArrivalLabel> rangeRaptorArrivals = rangeRaptor.getArrivals(target);
                if (!Vector::equals(raptorArrivals, rangeRaptorArrivals)) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "RAPTOR arrivals:" << std::endl;
                    std::cout << raptorArrivals << std::endl;
                    std::cout << "Range-RAPTOR arrivals:" << std::endl;
                    std::cout << rangeRaptorArrivals << std::endl;
                    std::cout << "RAPTOR journeys:" << std::endl;
                    std::cout << raptor.getJourneys(target) << std::endl;
                    std::cout << "Range-RAPTOR journeys:" << std::endl;
                    std::cout << rangeRaptor.getJourneys(target) << std::endl;
                    return;
                }
            }
            progress++;
        }
        raptor.getProfiler().printStatistics();
        rangeRaptor.getProfiler().printStatistics();
    }
};

class ValidateOneToManyDijkstraRAPTOR : public ParameterizedCommand {

public:
    ValidateOneToManyDijkstraRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateOneToManyDijkstraRAPTOR", "Validates journeys computed by one-to-many Dijkstra-RAPTOR by comparing to the Range-RAPTOR implementation on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Core CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("Core CH data"));

        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler> raptor(raptorData, ch);
        RAPTOR::RangeRAPTOR::DijkstraRAPTORModule<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler, false, false, false, false> rangeRaptor(raptorData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime);
            rangeRaptor.run(source, departureTime);
            for (const Vertex target : raptorData.stops()) {
                const std::vector<RAPTOR::ArrivalLabel> raptorArrivals = raptor.getArrivals(target);
                const std::vector<RAPTOR::ArrivalLabel> rangeRaptorArrivals = rangeRaptor.getArrivals(target);
                if (!Vector::equals(raptorArrivals, rangeRaptorArrivals)) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "RAPTOR arrivals:" << std::endl;
                    std::cout << raptorArrivals << std::endl;
                    std::cout << "Range-RAPTOR arrivals:" << std::endl;
                    std::cout << rangeRaptorArrivals << std::endl;
                    std::cout << "RAPTOR journeys:" << std::endl;
                    std::cout << raptor.getJourneys(target) << std::endl;
                    std::cout << "Range-RAPTOR journeys:" << std::endl;
                    std::cout << rangeRaptor.getJourneys(target) << std::endl;
                    return;
                }
            }
            progress++;
        }
        raptor.getProfiler().printStatistics();
        rangeRaptor.getProfiler().printStatistics();
    }
};

class ValidateULTRARAPTOR : public ParameterizedCommand {

public:
    ValidateULTRARAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRARAPTOR", "Validates journeys computed by ULTRA-RAPTOR by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("Dijkstra RAPTOR input file");
        addParameter("Shortcut RAPTOR input file");
        addParameter("Core CH data");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data dijkstraRaptorData = RAPTOR::Data::FromBinary(getParameter("Dijkstra RAPTOR input file"));
        dijkstraRaptorData.useImplicitDepartureBufferTimes();
        dijkstraRaptorData.printInfo();
        CH::CH coreCH(getParameter("Core CH data"));

        RAPTOR::Data shortcutRaptorData = RAPTOR::Data::FromBinary(getParameter("Shortcut RAPTOR input file"));
        shortcutRaptorData.useImplicitDepartureBufferTimes();
        shortcutRaptorData.printInfo();
        CH::CH bucketCH(getParameter("Bucket CH data"));

        RAPTOR::DijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler> dijkstraRaptor(dijkstraRaptorData, coreCH);
        RAPTOR::ULTRARAPTOR<RAPTOR::AggregateProfiler> ultraRaptor(shortcutRaptorData, bucketCH);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % coreCH.numVertices());
            const Vertex target(rand() % coreCH.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            dijkstraRaptor.run(source, departureTime, target);
            ultraRaptor.run(source, departureTime, target);
            const std::vector<RAPTOR::ArrivalLabel> dijkstraArrivals = dijkstraRaptor.getArrivals();
            const std::vector<RAPTOR::ArrivalLabel> ultraArrivals = ultraRaptor.getArrivals();
            if (!Vector::equals(dijkstraArrivals, ultraArrivals)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "Dijkstra arrivals:" << std::endl;
                std::cout << dijkstraArrivals << std::endl;
                std::cout << "ULTRA arrivals:" << std::endl;
                std::cout << ultraArrivals << std::endl;
                std::cout << "Dijkstra journeys:" << std::endl;
                std::cout << dijkstraRaptor.getJourneys() << std::endl;
                std::cout << "ULTRA journeys:" << std::endl;
                std::cout << ultraRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        dijkstraRaptor.getProfiler().printStatistics();
        ultraRaptor.getProfiler().printStatistics();
    }
};

class ValidateUPRAPTOR : public ParameterizedCommand {

public:
    ValidateUPRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUPRAPTOR", "Validates UP-RAPTOR by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("Regular RAPTOR data");
        addParameter("ULTRA RAPTOR data");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
        addParameter("Targets", { "Vertices", "Stops" });
        addParameter("Grouped rounds");
    }

    virtual void execute() noexcept {
        switch(getParameter<size_t>("Grouped rounds")) {
            case 0:
                run<0>();
                break;
            case 4:
                run<4>();
                break;
            case 6:
                run<6>();
                break;
            case 8:
                run<8>();
                break;
        }
    }

private:
    template<size_t GROUPED_ROUNDS>
    inline void run() const noexcept {
        RAPTOR::Data regularRaptorData = RAPTOR::Data::FromBinary(getParameter("Regular RAPTOR data"));
        regularRaptorData.useImplicitDepartureBufferTimes();
        regularRaptorData.printInfo();
        RAPTOR::Data ultraRaptorData = RAPTOR::Data::FromBinary(getParameter("ULTRA RAPTOR data"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        IndexedSet<false, Vertex> targetSet = getTargetSet(regularRaptorData);

        RAPTOR::Data reverseRegularRaptorData = regularRaptorData.reverseNetwork();
        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::SimpleProfiler> dijkstraRaptor(regularRaptorData, regularRaptorData.transferGraph, reverseRegularRaptorData.transferGraph);
        const bool reorder = getParameter<bool>("Reorder network?");
        using UPRAPTOR = RAPTOR::UPRAPTOR<GROUPED_ROUNDS, RAPTOR::SimpleProfiler>;
        UPRAPTOR upRaptor(ultraRaptorData, regularRaptorData.transferGraph, ch, targetSet, reorder, getParameter("Order") == "DFS");

        srand(42);
        for (size_t i = 0; i < getParameter<size_t>("Number of queries"); i++) {
            const Vertex source(rand() % regularRaptorData.transferGraph.numVertices());
            const int departureTime(rand() % 24 * 60 * 60);
            dijkstraRaptor.run(source, departureTime);
            upRaptor.run(source, departureTime);
            for (const Vertex target : targetSet) {
                std::vector<RAPTOR::ArrivalLabel> dijkstraArrivals = dijkstraRaptor.getArrivals(target);
                std::vector<RAPTOR::ArrivalLabel> ultraArrivals = upRaptor.getArrivals(target);
                if (!Vector::equals(dijkstraArrivals, ultraArrivals)) {
                    std::cout << "Query " << i << ": Arrivals not identical! " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "Dijkstra journeys:" << std::endl;
                    std::cout << dijkstraRaptor.getJourneys(target);
                    std::cout << "ULTRA-PHAST journeys:" << std::endl;
                    std::cout << upRaptor.getJourneys(target);
                    std::cout << std::endl;
                    return;
                }
            }
        }
    }

    inline IndexedSet<false, Vertex> getTargetSet(const RAPTOR::Data& data) const noexcept {
        if (getParameter("Targets") == "Stops") {
            IndexedSet<false, Vertex> result(data.transferGraph.numVertices());
            for (const StopId stop : data.stops()) {
                result.insert(stop);
            }
            return result;
        } else {
            return IndexedSet<false, Vertex>(Construct::Complete, data.transferGraph.numVertices());
        }
    }
};

class ValidateHL : public ParameterizedCommand {

public:
    ValidateHL(BasicShell& shell) :
        ParameterizedCommand(shell, "validateHL", "Validates paths computed by HL by comparing to Dijkstra on random queries.") {
        addParameter("Full graph");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        const TransferGraph fullGraph(getParameter("Full graph"));
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));

        Dijkstra<TransferGraph> dijkstra(fullGraph);
        HL::HLQuery<TransferGraph> hl(outHubs, inHubs);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % fullGraph.numVertices());
            const Vertex target(rand() % fullGraph.numVertices());
            dijkstra.run(source, target);
            hl.run(source, target);
            int dijkstraDistance = dijkstra.getDistance(target);
            if (dijkstraDistance == -1) dijkstraDistance = INFTY;
            const int hlDistance = hl.getDistance();
            if (dijkstraDistance != hlDistance) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << dijkstraDistance << " vs. " << hlDistance << std::endl;
                std::cout << "Dijkstra path:" << std::endl;
                Vector::printConcise(dijkstra.getPath(target));
                std::cout << std::endl;
                std::cout << "Source hubs:" << std::endl;
                const std::vector<Vertex> sourceHubs = hl.getHubs<FORWARD>(source);
                Vector::printConcise(sourceHubs);
                std::cout << std::endl;
                std::cout << "Target hubs:" << std::endl;
                const std::vector<Vertex> targetHubs = hl.getHubs<BACKWARD>(target);
                Vector::printConcise(targetHubs);
                std::cout << std::endl;
                std::cout << "Shared hubs:" << std::endl;
                Vector::printConcise(Vector::sortedIntersection(sourceHubs, targetHubs));
                std::cout << std::endl;
                return;
            }
            progress++;
        }
    }
};

class ValidateHLRAPTOR : public ParameterizedCommand {

public:
    ValidateHLRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateHLRAPTOR", "Validates journeys computed by HL-RAPTOR by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Core CH data");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH coreCH(getParameter("Core CH data"));
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));

        RAPTOR::DijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler> dijkstraRaptor(raptorData, coreCH);
        RAPTOR::HLRAPTOR<RAPTOR::AggregateProfiler> hlRaptor(raptorData, outHubs, inHubs);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % coreCH.numVertices());
            const Vertex target(rand() % coreCH.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            dijkstraRaptor.run(source, departureTime, target);
            hlRaptor.run(source, departureTime, target);
            const std::vector<RAPTOR::ArrivalLabel> dijkstraArrivals = dijkstraRaptor.getArrivals();
            const std::vector<RAPTOR::ArrivalLabel> ultraArrivals = hlRaptor.getArrivals();
            if (!Vector::equals(dijkstraArrivals, ultraArrivals)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "Dijkstra arrivals:" << std::endl;
                std::cout << dijkstraArrivals << std::endl;
                std::cout << "HL arrivals:" << std::endl;
                std::cout << ultraArrivals << std::endl;
                std::cout << "Dijkstra journeys:" << std::endl;
                std::cout << dijkstraRaptor.getJourneys() << std::endl;
                std::cout << "HL journeys:" << std::endl;
                std::cout << hlRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        dijkstraRaptor.getProfiler().printStatistics();
        hlRaptor.getProfiler().printStatistics();
    }
};

class ValidateULTRATripBased : public ParameterizedCommand {

public:
    ValidateULTRATripBased(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRATripBased", "Validates journeys computed by ULTRA-TripBased by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("Dijkstra RAPTOR input file");
        addParameter("Trip-Based input file");
        addParameter("Core CH data");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data dijkstraRaptorData = RAPTOR::Data::FromBinary(getParameter("Dijkstra RAPTOR input file"));
        dijkstraRaptorData.useImplicitDepartureBufferTimes();
        dijkstraRaptorData.printInfo();
        CH::CH coreCH(getParameter("Core CH data"));

        TripBased::Data tripBasedData(getParameter("Trip-Based input file"));
        tripBasedData.printInfo();
        CH::CH bucketCH(getParameter("Bucket CH data"));

        RAPTOR::DijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler> dijkstraRaptor(dijkstraRaptorData, coreCH);
        TripBased::Query<TripBased::AggregateProfiler> ultraTripBased(tripBasedData, bucketCH);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % coreCH.numVertices());
            const Vertex target(rand() % coreCH.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            dijkstraRaptor.run(source, departureTime, target);
            ultraTripBased.run(source, departureTime, target);
            const std::vector<RAPTOR::ArrivalLabel> dijkstraArrivals = dijkstraRaptor.getArrivals();
            const std::vector<RAPTOR::ArrivalLabel> ultraArrivals = ultraTripBased.getArrivals();
            if (!Vector::equals(dijkstraArrivals, ultraArrivals)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "Dijkstra arrivals:" << std::endl;
                std::cout << dijkstraArrivals << std::endl;
                std::cout << "ULTRA arrivals:" << std::endl;
                std::cout << ultraArrivals << std::endl;
                std::cout << "Dijkstra journeys:" << std::endl;
                std::cout << dijkstraRaptor.getJourneys() << std::endl;
                std::cout << "ULTRA journeys:" << std::endl;
                std::cout << ultraTripBased.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        dijkstraRaptor.getProfiler().printStatistics();
        ultraTripBased.getProfiler().printStatistics();
    }
};

class ValidateHLULTRATripBased : public ParameterizedCommand {

public:
    ValidateHLULTRATripBased(BasicShell& shell) :
        ParameterizedCommand(shell, "validateHLULTRATripBased", "Validates journeys computed by HL-ULTRA-TripBased by comparing to ULTRA-TripBased on random queries.") {
        addParameter("Trip-Based input file");
        addParameter("Bucket CH data");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        TripBased::Data tripBasedData(getParameter("Trip-Based input file"));
        tripBasedData.printInfo();
        CH::CH bucketCH(getParameter("Bucket CH data"));
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));

        TripBased::Query<TripBased::AggregateProfiler> ultraTripBased(tripBasedData, bucketCH);
        TripBased::HLQuery<TripBased::AggregateProfiler> hlUltraTripBased(tripBasedData, outHubs, inHubs);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % bucketCH.numVertices());
            const Vertex target(rand() % bucketCH.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            ultraTripBased.run(source, departureTime, target);
            hlUltraTripBased.run(source, departureTime, target);
            const std::vector<RAPTOR::ArrivalLabel> ultraArrivals = ultraTripBased.getArrivals();
            const std::vector<RAPTOR::ArrivalLabel> hlArrivals = hlUltraTripBased.getArrivals();
            if (!Vector::equals(ultraArrivals, hlArrivals)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "ULTRA arrivals:" << std::endl;
                std::cout << ultraArrivals << std::endl;
                std::cout << "HL-ULTRA arrivals:" << std::endl;
                std::cout << hlArrivals << std::endl;
                std::cout << "ULTRA journeys:" << std::endl;
                std::cout << ultraTripBased.getJourneys() << std::endl;
                std::cout << "HL-ULTRA journeys:" << std::endl;
                std::cout << hlUltraTripBased.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        ultraTripBased.getProfiler().printStatistics();
        hlUltraTripBased.getProfiler().printStatistics();
    }
};

class ValidateTransitiveProfileRAPTOR : public ParameterizedCommand {

public:
    ValidateTransitiveProfileRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateTransitiveProfileRAPTOR", "Compares results of transitive Alternating RAPTOR and Range-RAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();

        RAPTOR::RangeRAPTOR::RangeTransitiveRAPTOR<false, true, true, true> rangeRaptor(raptorData, reverseRaptorData);
        RAPTOR::AlternatingTransitiveRAPTOR<false, false> alternatingRaptor(raptorData, reverseRaptorData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % raptorData.numberOfStops());
            const StopId target(rand() % raptorData.numberOfStops());
            rangeRaptor.run(source, target, 0, 24 * 60 * 60);
            alternatingRaptor.run(source, target, 0, 24 * 60 * 60);
            RAPTOR::ProfileHandle rangeProfile = rangeRaptor.getProfileHandle(target);
            RAPTOR::ProfileHandle alternatingProfile = alternatingRaptor.getProfileHandle(target);
            rangeProfile.sort();
            alternatingProfile.sort();
            if (rangeProfile != alternatingProfile) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "Range-RAPTOR profile:" << std::endl;
                std::cout << rangeProfile << std::endl;
                std::cout << "Alternating RAPTOR profile:" << std::endl;
                std::cout << alternatingProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateDijkstraProfileRAPTOR : public ParameterizedCommand {

public:
    ValidateDijkstraProfileRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateDijkstraProfileRAPTOR", "Compares results of Alternating RAPTOR and Range-RAPTOR with unrestricted walking on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        CH::CH ch(getParameter("CH data"));

        RAPTOR::RangeRAPTOR::RangeDijkstraRAPTOR<RAPTOR::BucketCHInitialTransfers, false, true, false, true, true> rangeRaptor(raptorData, reverseRaptorData, ch);
        RAPTOR::AlternatingDijkstraRAPTOR<RAPTOR::BucketCHInitialTransfers, false, false> alternatingRaptor(raptorData, reverseRaptorData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            rangeRaptor.run(source, target, 0, 24 * 60 * 60);
            alternatingRaptor.run(source, target, 0, 24 * 60 * 60);
            RAPTOR::ProfileHandle rangeProfile = rangeRaptor.getProfileHandle(target);
            RAPTOR::ProfileHandle alternatingProfile = alternatingRaptor.getProfileHandle(target);
            rangeProfile.sort();
            alternatingProfile.sort();
            if (rangeProfile != alternatingProfile) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "Range-RAPTOR profile:" << std::endl;
                std::cout << rangeProfile << std::endl;
                std::cout << "Alternating RAPTOR profile:" << std::endl;
                std::cout << alternatingProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateULTRAProfileRAPTOR : public ParameterizedCommand {

public:
    ValidateULTRAProfileRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRAProfileRAPTOR", "Compares results of ULTRA variants Alternating RAPTOR and Range-RAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        CH::CH ch(getParameter("CH data"));

        RAPTOR::RangeRAPTOR::RangeULTRARAPTOR<false, true, true> rangeRaptor(raptorData, reverseRaptorData, ch);
        RAPTOR::AlternatingULTRARAPTOR<false> alternatingRaptor(raptorData, reverseRaptorData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            rangeRaptor.run(source, target, 0, 24 * 60 * 60);
            alternatingRaptor.run(source, target, 0, 24 * 60 * 60);
            RAPTOR::ProfileHandle rangeProfile = rangeRaptor.getProfileHandle(target);
            RAPTOR::ProfileHandle alternatingProfile = alternatingRaptor.getProfileHandle(target);
            rangeProfile.sort();
            alternatingProfile.sort();
            if (rangeProfile != alternatingProfile) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "Range-RAPTOR profile:" << std::endl;
                std::cout << rangeProfile << std::endl;
                std::cout << "Alternating RAPTOR profile:" << std::endl;
                std::cout << alternatingProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateBoundedMcRAPTOR : public ParameterizedCommand {

public:
    ValidateBoundedMcRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateBoundedMcRAPTOR", "Compares results of Bounded McRAPTOR and McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Transitive?");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Transitive?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

private:

    template<bool TRANSITIVE>
    inline void run() const noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        const RAPTOR::Data reverseData = raptorData.reverseNetwork();

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::McRAPTOR<true, TRANSITIVE, RAPTOR::AggregateProfiler> mcRaptor(raptorData);
        RAPTOR::BoundedMcRAPTOR<RAPTOR::AggregateProfiler> boundedMcRaptor(raptorData, reverseData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % raptorData.numberOfStops());
            const StopId target(rand() % raptorData.numberOfStops());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcRaptor.run(source, departureTime, target, arrivalSlack, tripSlack);
            boundedMcRaptor.verify(arrivalSlack, tripSlack, departureTime);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcRaptor.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, boundedMcRaptor.getAnchorLabels(), departureTime, arrivalSlack, tripSlack)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Bounded McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : boundedMcRaptor.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                std::cout << "McRAPTOR journeys:" << std::endl;
                std::cout << mcRaptor.getJourneys() << std::endl;
                std::cout << "Bounded McRAPTOR journeys:" << std::endl;
                std::cout << boundedMcRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        boundedMcRaptor.getProfiler().printStatistics();
    }
};

class ValidateUBMRAPTOR : public ParameterizedCommand {

public:
    ValidateUBMRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUBMRAPTOR", "Compares results of UBM-RAPTOR and ULTRA-McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));
        const RAPTOR::Data reverseData = raptorData.reverseNetwork();

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        RAPTOR::UBMRAPTOR<RAPTOR::AggregateProfiler> boundedMcRaptor(raptorData, reverseData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcRaptor.run(source, departureTime, target, arrivalSlack, tripSlack);
            boundedMcRaptor.verify(arrivalSlack, tripSlack, departureTime);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcRaptor.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, boundedMcRaptor.getAnchorLabels(), departureTime, arrivalSlack, tripSlack)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Bounded McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : boundedMcRaptor.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                std::cout << "McRAPTOR journeys:" << std::endl;
                std::cout << mcRaptor.getJourneys() << std::endl;
                std::cout << "Bounded McRAPTOR journeys:" << std::endl;
                std::cout << boundedMcRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        boundedMcRaptor.getProfiler().printStatistics();
    }
};

class ValidateNewUBMRAPTOR : public ParameterizedCommand {

public:
    ValidateNewUBMRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateNewUBMRAPTOR", "Compares results of new UBM-RAPTOR and ULTRA-McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival factor");
        addParameter("Trips per hour");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));
        const RAPTOR::Data reverseData = raptorData.reverseNetwork();

        const double arrivalFactor = getParameter<double>("Arrival factor");
        const double tripsPerHour = getParameter<double>("Trips per hour");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        RAPTOR::NewUBMRAPTOR<RAPTOR::AggregateProfiler> boundedMcRaptor(raptorData, reverseData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcRaptor.run(source, departureTime, target, arrivalFactor, tripsPerHour);
            boundedMcRaptor.verify(arrivalFactor, tripsPerHour);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcRaptor.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            const std::vector<RAPTOR::WalkingParetoLabel> fullHullLabels = getPreferenceOptimalLabels(fullBag.labels, arrivalFactor, tripsPerHour);
            const std::vector<RAPTOR::WalkingParetoLabel> boundedHullLabels = getPreferenceOptimalLabels(boundedBag.labels, arrivalFactor, tripsPerHour);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, boundedMcRaptor.getAnchorLabels(), arrivalFactor, tripsPerHour) || fullHullLabels != boundedHullLabels) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << "Hull labels:" << std::endl;
                std::cout << fullHullLabels << std::endl;
                std::cout << std::endl;
                std::cout << "Bounded McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << "Hull labels:" << std::endl;
                std::cout << boundedHullLabels << std::endl;
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedMcRaptor.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                return;
            }

            progress++;
        }
        mcRaptor.getProfiler().printStatistics();
        boundedMcRaptor.getProfiler().printStatistics();
    }
};

class ValidateUBMTBRAPTOR : public ParameterizedCommand {

public:
    ValidateUBMTBRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUBMTBRAPTOR", "Compares results of UBM-TB-RAPTOR and ULTRA-McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Bounded forward Trip-Based input file");
        addParameter("Bounded backward Trip-Based input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        const TripBased::Data forwardBoundedData(getParameter("Bounded forward Trip-Based input file"));
        forwardBoundedData.printInfo();
        const TripBased::Data backwardBoundedData(getParameter("Bounded backward Trip-Based input file"));
        backwardBoundedData.printInfo();
        const CH::CH ch(getParameter("CH data"));

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        RAPTOR::UBMTBRAPTOR<RAPTOR::AggregateProfiler> boundedMcRaptor(raptorData, forwardBoundedData, backwardBoundedData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcRaptor.run(source, departureTime, target, arrivalSlack, tripSlack);
            boundedMcRaptor.verify(arrivalSlack, tripSlack, departureTime);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcRaptor.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, boundedMcRaptor.getAnchorLabels(), departureTime, arrivalSlack, tripSlack)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Bounded McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : boundedMcRaptor.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                std::cout << "McRAPTOR journeys:" << std::endl;
                std::cout << mcRaptor.getJourneys() << std::endl;
                std::cout << "Bounded McRAPTOR journeys:" << std::endl;
                std::cout << boundedMcRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        boundedMcRaptor.getProfiler().printStatistics();
    }
};

class ValidateUBMHydRA : public ParameterizedCommand {

public:
    ValidateUBMHydRA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUBMHydRA", "Compares results of UBM-HydRA and ULTRA-McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Trip-Based input file");
        addParameter("Bounded forward Trip-Based input file");
        addParameter("Bounded backward Trip-Based input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        const TripBased::Data tripBasedData(getParameter("Trip-Based input file"));
        tripBasedData.printInfo();
        const TripBased::Data forwardBoundedData(getParameter("Bounded forward Trip-Based input file"));
        forwardBoundedData.printInfo();
        const TripBased::Data backwardBoundedData(getParameter("Bounded backward Trip-Based input file"));
        backwardBoundedData.printInfo();
        const CH::CH ch(getParameter("CH data"));

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        RAPTOR::UBMHydRA<RAPTOR::AggregateProfiler> hydra(tripBasedData, forwardBoundedData, backwardBoundedData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            hydra.run(source, departureTime, target, arrivalSlack, tripSlack);
            hydra.verify(arrivalSlack, tripSlack, departureTime);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = hydra.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, hydra.getAnchorLabels(), departureTime, arrivalSlack, tripSlack)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "HydRA:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : hydra.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                std::cout << "McRAPTOR journeys:" << std::endl;
                std::cout << mcRaptor.getJourneys() << std::endl;
                std::cout << "HydRA journeys:" << std::endl;
                std::cout << hydra.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        hydra.getProfiler().printStatistics();
    }
};

class ValidateULTRAMcRAPTOR : public ParameterizedCommand {

public:
    ValidateULTRAMcRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRAMcRAPTOR", "Compares results of MCR and ULTRA-McRAPTOR on random queries.") {
        addParameter("MCR input file");
        addParameter("ULTRA-RAPTOR input file");
        addParameter("Core CH data");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data mcrData = RAPTOR::Data::FromBinary(getParameter("MCR input file"));
        mcrData.useImplicitDepartureBufferTimes();
        mcrData.printInfo();
        CH::CH coreCH(getParameter("Core CH data"));

        RAPTOR::Data ultraRaptorData = RAPTOR::Data::FromBinary(getParameter("ULTRA-RAPTOR input file"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();
        CH::CH bucketCH(getParameter("Bucket CH data"));

        RAPTOR::MCR<true, RAPTOR::AggregateProfiler> mcr(mcrData, coreCH);
        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> ultraMcRaptor(ultraRaptorData, bucketCH);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % coreCH.numVertices());
            const StopId target(rand() % coreCH.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcr.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> mcrResults = mcr.getResults();
            sort(mcrResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> mcrBag(mcrResults);
            ultraMcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> ultraResults = ultraMcRaptor.getResults();
            sort(ultraResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> ultraBag(ultraResults);
            if (!Vector::equals(mcrBag.labels, ultraBag.labels)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "MCR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : mcrBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "ULTRAMcRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : ultraBag) {
                    std::cout << label << std::endl;
                }
                std::cout << "MCR journeys:" << std::endl;
                std::cout << mcr.getJourneys() << std::endl;
                std::cout << "ULTRAMcRAPTOR journeys:" << std::endl;
                std::cout << ultraMcRaptor.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        mcr.getProfiler().printStatistics();
        ultraMcRaptor.getProfiler().printStatistics();
    }
};

class ValidateBoundedULTRAMcTripBased : public ParameterizedCommand {

public:
    ValidateBoundedULTRAMcTripBased(BasicShell& shell) :
        ParameterizedCommand(shell, "validateBoundedULTRAMcTripBased", "Compares results of Bounded ULTRA-McRAPTOR and ULTRA-McRAPTOR on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("Trip-Based input file");
        addParameter("Bounded forward Trip-Based input file");
        addParameter("Bounded backward Trip-Based input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        const RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();

        TripBased::Data tripBasedData(getParameter("Trip-Based input file"));
        tripBasedData.printInfo();
        TripBased::Data forwardBoundedData(getParameter("Bounded forward Trip-Based input file"));
        forwardBoundedData.printInfo();
        TripBased::Data backwardBoundedData(getParameter("Bounded backward Trip-Based input file"));
        backwardBoundedData.printInfo();
        CH::CH ch(getParameter("CH data"));

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        TripBased::BoundedMcQuery<TripBased::AggregateProfiler> boundedMcTripBased(tripBasedData, forwardBoundedData, backwardBoundedData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcTripBased.run(source, departureTime, target, arrivalSlack, tripSlack);
            boundedMcTripBased.verify(arrivalSlack, tripSlack, departureTime);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcTripBased.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            if (!compareBoundedResults(fullBag.labels, boundedBag.labels, boundedMcTripBased.getAnchorLabels(), departureTime, arrivalSlack, tripSlack)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "McRAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : fullBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Bounded McTripBased:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : boundedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "Anchor labels:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : boundedMcTripBased.getAnchorLabels()) {
                    std::cout << label << std::endl;
                }
                std::cout << "McRAPTOR journeys:" << std::endl;
                std::cout << mcRaptor.getJourneys() << std::endl;
                std::cout << "Bounded McTripBased journeys:" << std::endl;
                std::cout << boundedMcTripBased.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        boundedMcTripBased.getProfiler().printStatistics();
    }
};

class ValidateULTRAMcTripBased : public ParameterizedCommand {

public:
    ValidateULTRAMcTripBased(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRAMcTripBased", "Compares results of ULTRA-McRAPTOR and ULTRA-McTripBased on random queries.") {
        addParameter("RAPTOR input file");
        addParameter("TripBased input file");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();

        TripBased::Data tripBasedData(getParameter("TripBased input file"));
        tripBasedData.printInfo();
        CH::CH ch(getParameter("Bucket CH data"));

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> raptor(raptorData, ch);
        TripBased::McQuery<TripBased::AggregateProfiler> tripBased(tripBasedData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % ch.numVertices());
            const StopId target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> raptorResults = raptor.getResults();
            sort(raptorResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> raptorBag(raptorResults);
            tripBased.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> tripBasedResults = tripBased.getResults();
            sort(tripBasedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> tripBasedBag(tripBasedResults);
            if (!Vector::equals(raptorBag.labels, tripBasedBag.labels)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "RAPTOR:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : raptorBag) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "TripBased:" << std::endl;
                for (const RAPTOR::WalkingParetoLabel& label : tripBasedBag) {
                    std::cout << label << std::endl;
                }
                std::cout << "RAPTOR journeys:" << std::endl;
                std::cout << raptor.getJourneys() << std::endl;
                std::cout << "TripBased journeys:" << std::endl;
                std::cout << tripBased.getJourneys() << std::endl;
                return;
            }
            progress++;
        }
        raptor.getProfiler().printStatistics();
        tripBased.getProfiler().printStatistics();
    }
};


class ValidateTransitiveCSA : public ParameterizedCommand {

public:
    ValidateTransitiveCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateTransitiveCSA", "Validates journeys computed by CSA by comparing to RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::RAPTOR<false, RAPTOR::NoProfiler, true, false, false> raptor(raptorData);
        CSA::CSA<true, CSA::NoProfiler> csa(csaData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % csaData.numberOfStops());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime);
            csa.run(source, departureTime);
            for (const StopId target : csaData.stops()) {
                const int csaArrivalTime = csa.getEarliestArrivalTime(target);
                const int raptorArrivalTime = raptor.getEarliestArrivalTime(target);
                if (csaArrivalTime != raptorArrivalTime) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "CSA journey:" << std::endl;
                    std::cout << csa.getJourney(target) << std::endl;
                    std::cout << "RAPTOR journey:" << std::endl;
                    std::cout << raptor.getEarliestJourney(target) << std::endl;
                }
            }
            progress++;
        }
    }
};

inline bool profilesEqual(const RAPTOR::ProfileHandle& csaProfile, const RAPTOR::ProfileHandle& raptorProfile) noexcept {
    if (csaProfile.size() != raptorProfile.size()) return false;
    for (size_t i = 0; i < csaProfile.size(); i++) {
        if (csaProfile[i].departureTime >= 24 * 60 * 60) {
            if (csaProfile[i].departureTime < raptorProfile[i].departureTime) return false;
        } else {
            if (csaProfile[i].departureTime != raptorProfile[i].departureTime) return false;
        }
        if (csaProfile[i].arrivalTime != raptorProfile[i].arrivalTime) return false;
    }
    return true;
}

class ValidateTransitiveProfileCSA : public ParameterizedCommand {

public:
    ValidateTransitiveProfileCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateTransitiveProfileCSA", "Validates journeys computed by Profile-CSA by comparing to Profile-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("Contiguous profiles?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Contiguous profiles?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

private:
    template<bool CONTIGUOUS_PROFILES>
    inline void run() const noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        CSA::TransferGraph reverseCsaGraph = csaData.transferGraph;
        reverseCsaGraph.revert();

        RAPTOR::RangeRAPTOR::RangeTransitiveRAPTOR<false, true, true, true> raptor(raptorData, reverseRaptorData);
        CSA::ProfileCSA<true, CSA::NoProfiler, CONTIGUOUS_PROFILES> csa(csaData, reverseCsaGraph);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % csaData.numberOfStops());
            const StopId target(rand() % csaData.numberOfStops());
            raptor.run(source, target);
            csa.run(source, target);

            const RAPTOR::ProfileHandle csaProfile = csa.getProfileHandle(source);
            RAPTOR::ProfileHandle raptorProfile = raptor.getProfileHandle(target).reduceToEarliestArrivalTime();
            if (!profilesEqual(csaProfile, raptorProfile)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "CSA profile:" << std::endl;
                std::cout << csaProfile << std::endl;
                std::cout << "RAPTOR profile:" << std::endl;
                std::cout << raptorProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateParetoCSA : public ParameterizedCommand {

public:
    ValidateParetoCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateParetoCSA", "Validates journeys computed by Pareto-CSA by comparing to RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("Number of trips", { "8", "16" });
    }

    virtual void execute() noexcept {
        if (getParameter<size_t>("Number of trips") == 8) {
            run<8>();
        } else {
            run<16>();
        }
    }

private:
    template<size_t MAX_TRIPS>
    inline void run() const noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::RAPTOR<false, RAPTOR::AggregateProfiler, true, false, false> raptor(raptorData);
        CSA::ParetoCSA<true, CSA::AggregateProfiler, MAX_TRIPS> csa(csaData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % csaData.numberOfStops());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, noStop, MAX_TRIPS - 1);
            csa.run(source, departureTime);
            for (const StopId target : csaData.stops()) {
                std::vector<int> csaArrivalTimes = csa.getArrivalTimes(target);
                std::vector<int> raptorArrivalTimes = raptor.getArrivalTimes(target);
                if (!arrivalsEqual(csaArrivalTimes, raptorArrivalTimes)) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "CSA arrivals: ";
                    Vector::printConcise(csaArrivalTimes);
                    std::cout << std::endl;
                    std::cout << "RAPTOR arrivals: ";
                    Vector::printConcise(raptorArrivalTimes);
                    std::cout << std::endl;
                    std::cout << "CSA journeys:" << std::endl;
                    std::cout << csa.getJourneys(target) << std::endl;
                    std::cout << "RAPTOR journeys:" << std::endl;
                    std::cout << raptor.getJourneys(target) << std::endl;
                    return;
                }
            }
            progress++;
        }
        raptor.getProfiler().printStatistics();
        csa.getProfiler().printStatistics();
    }

    inline static bool arrivalsEqual(const std::vector<int>& csaArrivalTimes, const std::vector<int>& raptorArrivalTimes) noexcept {
        if (csaArrivalTimes.size() < raptorArrivalTimes.size()) return false;
        int minimum = never;
        for (size_t i = 0; i < raptorArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] != raptorArrivalTimes[i]) return false;
            minimum = std::min(csaArrivalTimes[i], minimum);
        }
        for (size_t i = raptorArrivalTimes.size(); i < csaArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] < minimum) return false;
        }
        return true;
    }
};


class ValidateDijkstraCSA : public ParameterizedCommand {

public:
    ValidateDijkstraCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateDijkstraCSA", "Validates journeys computed by Dijkstra-CSA by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("CH data");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        CH::CH ch(getParameter("CH data"));
        RAPTOR::DijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::NoProfiler, true, false, false> raptor(raptorData, ch);
        CSA::DijkstraCSA<RAPTOR::CoreCHInitialTransfers, true, CSA::NoProfiler> csa(csaData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, target);
            csa.run(source, departureTime, target);
            const int csaArrivalTime = csa.getEarliestArrivalTime(target);
            const int raptorArrivalTime = raptor.getEarliestArrivalTime(target);
            if (csaArrivalTime != raptorArrivalTime) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << csaArrivalTime << " vs. " << raptorArrivalTime << std::endl;
                std::cout << "CSA journey:" << std::endl;
                std::cout << csa.getJourney(target) << std::endl;
                std::cout << "RAPTOR journey:" << std::endl;
                std::cout << raptor.getEarliestJourney(target) << std::endl;
                return;
            }
            progress++;
        }
    }
};

class ValidateOneToAllDijkstraCSA : public ParameterizedCommand {

public:
    ValidateOneToAllDijkstraCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateOneToAllDijkstraCSA", "Validates journeys computed by one-to-all Dijkstra-CSA by comparing to one-to-all Dijkstra-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::TransferGraph backwardGraph = raptorData.transferGraph;
        backwardGraph.revert();
        RAPTOR::RangeRAPTOR::DijkstraRAPTORModule<RAPTOR::DijkstraInitialTransfers, RAPTOR::NoProfiler, false, false, false, false> raptor(raptorData, raptorData.transferGraph, backwardGraph);
        CSA::OneToAllDijkstraCSA<CSA::TransferGraph, true, CSA::NoProfiler> csa(csaData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % raptorData.transferGraph.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime);
            csa.run(source, departureTime);
            for (const StopId target : csaData.stops()) {
                const int csaArrivalTime = csa.getEarliestArrivalTime(target);
                const int raptorArrivalTime = raptor.getEarliestArrivalTime(target);
                if (csaArrivalTime != raptorArrivalTime) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << csaArrivalTime << " vs. " << raptorArrivalTime << std::endl;
                    std::cout << "CSA journey:" << std::endl;
                    std::cout << csa.getJourney(target) << std::endl;
                    std::cout << "RAPTOR journey:" << std::endl;
                    std::cout << raptor.getEarliestJourney(target) << std::endl;
                    return;
                }
            }
            progress++;
        }
    }
};

class ValidateDijkstraProfileCSA : public ParameterizedCommand {

public:
    ValidateDijkstraProfileCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateDijkstraProfileCSA", "Validates journeys computed by Dijkstra Profile-CSA by comparing to Dijkstra Profile-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Contiguous profiles?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Contiguous profiles?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

private:
    template<bool CONTIGUOUS_PROFILES>
    inline void run() const noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        CSA::TransferGraph reverseCsaGraph = csaData.transferGraph;
        reverseCsaGraph.revert();

        CH::CH ch(getParameter("CH data"));

        RAPTOR::RangeRAPTOR::RangeDijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, false, true, false, true, true> raptor(raptorData, reverseRaptorData, ch);
        CSA::DijkstraProfileCSA<RAPTOR::CoreCHInitialTransfers, true, CSA::NoProfiler, CONTIGUOUS_PROFILES> csa(csaData, reverseCsaGraph, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            raptor.run(source, target);
            csa.run(source, target);

            const RAPTOR::ProfileHandle csaProfile = csa.getProfileHandle(source);
            RAPTOR::ProfileHandle raptorProfile = raptor.getProfileHandle(target).reduceToEarliestArrivalTime();
            if (!profilesEqual(csaProfile, raptorProfile)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "CSA profile:" << std::endl;
                std::cout << csaProfile << std::endl;
                std::cout << "RAPTOR profile:" << std::endl;
                std::cout << raptorProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateOneToAllDijkstraProfileCSA : public ParameterizedCommand {

public:
    ValidateOneToAllDijkstraProfileCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateOneToAllDijkstraProfileCSA", "Validates journeys computed by one-to-all Dijkstra-Profile-CSA by comparing to one-to-one Dijkstra-Profile-CSA on random queries.") {
        addParameter("CSA input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CSA::TransferGraph reverseCsaGraph = csaData.transferGraph;
        reverseCsaGraph.revert();

        CSA::DijkstraProfileCSA<RAPTOR::DijkstraInitialTransfers, true, CSA::NoProfiler, true> oneToOneCSA(csaData, csaData.transferGraph, reverseCsaGraph);
        CSA::OneToAllDijkstraProfileCSA<CSA::NoProfiler, true> oneToAllCSA(csaData, reverseCsaGraph);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % csaData.transferGraph.numVertices());
            const Vertex target(rand() % csaData.transferGraph.numVertices());
            oneToOneCSA.run(source, target);
            oneToAllCSA.run(target);

            const RAPTOR::ProfileHandle oneToOneProfile = oneToOneCSA.getProfileHandle(source);
            const RAPTOR::ProfileHandle oneToAllProfile = oneToAllCSA.getProfileHandle(source);
            if (!profilesEqual(oneToOneProfile, oneToAllProfile)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "One-to-one profile:" << std::endl;
                std::cout << oneToOneProfile << std::endl;
                std::cout << "One-to-all profile:" << std::endl;
                std::cout << oneToAllProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateUPCSA : public ParameterizedCommand {

public:
    ValidateUPCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUPCSA", "Validates UP-CSA by comparing to Dijkstra-CSA on random queries.") {
        addParameter("Regular CSA data");
        addParameter("Regular upward graph");
        addParameter("ULTRA CSA data");
        addParameter("ULTRA CH data");
        addParameter("Number of queries");
        addParameter("Initial transfers", { "Bucket", "PHAST" });
        addParameter("Final transfers", { "Bucket", "PHAST" });
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
        addParameter("Targets", { "Vertices", "Stops" });
    }

    virtual void execute() noexcept {
        if (getParameter("Initial transfers") == "Bucket") {
            chooseFinalTransfers<true>();
        } else {
            chooseFinalTransfers<false>();
        }
    }

private:
    template<bool USE_STOP_BUCKETS>
    inline void chooseFinalTransfers() const noexcept {
        if (getParameter("Final transfers") == "Bucket") {
            run<USE_STOP_BUCKETS, true>();
        } else {
            run<USE_STOP_BUCKETS, false>();
        }
    }

    template<bool USE_STOP_BUCKETS, bool USE_TARGET_BUCKETS>
    inline void run() const noexcept {
        CSA::Data regularCSAData = CSA::Data::FromBinary(getParameter("Regular CSA data"));
        regularCSAData.sortConnectionsAscending();
        regularCSAData.printInfo();
        CHGraph upwardGraph(getParameter("Regular upward graph"));
        CSA::Data ultraCSAData = CSA::Data::FromBinary(getParameter("ULTRA CSA data"));
        ultraCSAData.sortConnectionsAscending();
        ultraCSAData.printInfo();
        CH::CH ch(getParameter("ULTRA CH data"));

        IndexedSet<false, Vertex> targetSet = getTargetSet(regularCSAData);

        CSA::OneToAllDijkstraCSA<CHGraph, true, CSA::SimpleProfiler> dijkstraCSA(regularCSAData, upwardGraph, Weight);
        using UPCSA = CSA::UPCSA<USE_STOP_BUCKETS, USE_TARGET_BUCKETS, true, CSA::SimpleProfiler>;
        const bool reorder = getParameter<bool>("Reorder network?");
        UPCSA upCSA(ultraCSAData, ch, targetSet, reorder, getParameter("Order") == "DFS");

        srand(42);
        for (size_t i = 0; i < getParameter<size_t>("Number of queries"); i++) {
            const Vertex source(rand() % ch.numVertices());
            const int departureTime(rand() % 24 * 60 * 60);
            dijkstraCSA.run(source, departureTime);
            upCSA.run(source, departureTime);
            for (const Vertex target : targetSet) {
                const int dijkstraArrivalTime = dijkstraCSA.getEarliestArrivalTime(target);
                const int ultraArrivalTime = upCSA.getEarliestArrivalTime(target);
                if (dijkstraArrivalTime != ultraArrivalTime) {
                    std::cout << "Query " << i << ": Arrivals not identical! " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "Dijkstra journey:" << std::endl;
                    std::cout << dijkstraCSA.getJourney(target);
                    std::cout << "UP journey:" << std::endl;
                    std::cout << upCSA.getJourney(target);
                    std::cout << std::endl;
                    return;
                }
            }
        }
    }

    inline IndexedSet<false, Vertex> getTargetSet(const CSA::Data& data) const noexcept {
        if (getParameter("Targets") == "Stops") {
            IndexedSet<false, Vertex> result(data.transferGraph.numVertices());
            for (const StopId stop : data.stops()) {
                result.insert(stop);
            }
            return result;
        } else {
            return IndexedSet<false, Vertex>(Construct::Complete, data.transferGraph.numVertices());
        }
    }
};

class ValidateULTRACSA : public ParameterizedCommand {

public:
    ValidateULTRACSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRACSA", "Validates journeys computed by ULTRA-CSA by comparing to ULTRA-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("CH data");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        CH::CH ch(getParameter("CH data"));
        RAPTOR::ULTRARAPTOR<RAPTOR::NoProfiler, false> raptor(raptorData, ch);
        CSA::ULTRACSA<true, CSA::NoProfiler> csa(csaData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, target);
            csa.run(source, departureTime, target);
            const int csaArrivalTime = csa.getEarliestArrivalTime(target);
            const int raptorArrivalTime = raptor.getEarliestArrivalTime(target);
            if (csaArrivalTime != raptorArrivalTime) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "CSA journey:" << std::endl;
                std::cout << csa.getJourney(target) << std::endl;
                std::cout << "RAPTOR journey:" << std::endl;
                std::cout << raptor.getEarliestJourney(target) << std::endl;
                return;
            }
            progress++;
        }
    }
};

class ValidateHLCSA : public ParameterizedCommand {

public:
    ValidateHLCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateHLCSA", "Validates journeys computed by HL-CSA by comparing to Dijkstra-CSA on random queries.") {
        addParameter("CSA input file");
        addParameter("Core CH data");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CH::CH ch(getParameter("Core CH data"));
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));

        CSA::DijkstraCSA<RAPTOR::CoreCHInitialTransfers, true, CSA::AggregateProfiler> dijkstraCSA(csaData, ch);
        CSA::HLCSA<CSA::AggregateProfiler> hlCSA(csaData, outHubs, inHubs);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            dijkstraCSA.run(source, departureTime, target);
            hlCSA.run(source, departureTime, target);
            const int dijkstraArrivalTime = dijkstraCSA.getEarliestArrivalTime(target);
            const int hlArrivalTime = hlCSA.getEarliestArrivalTime(target);
            if (dijkstraArrivalTime != hlArrivalTime) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "Dijkstra-CSA journey:" << std::endl;
                std::cout << dijkstraCSA.getJourney(target) << std::endl;
                std::cout << "HL-CSA journey:" << std::endl;
                std::cout << hlCSA.getJourney(target) << std::endl;
                return;
            }
            progress++;
        }

        dijkstraCSA.getProfiler().printStatistics();
        hlCSA.getProfiler().printStatistics();
    }
};

class ValidateULTRAProfileCSA : public ParameterizedCommand {

public:
    ValidateULTRAProfileCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateULTRAProfileCSA", "Validates journeys computed by ULTRA Profile-CSA by comparing to ULTRA Profile-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Contiguous profiles?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Contiguous profiles?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

private:
    template<bool CONTIGUOUS_PROFILES>
    inline void run() const noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        CSA::TransferGraph reverseCsaGraph = csaData.transferGraph;
        reverseCsaGraph.revert();

        CH::CH ch(getParameter("CH data"));

        RAPTOR::RangeRAPTOR::RangeULTRARAPTOR<false, true, true> raptor(raptorData, reverseRaptorData, ch);
        CSA::ULTRAProfileCSA<true, CSA::NoProfiler, CONTIGUOUS_PROFILES> csa(csaData, reverseCsaGraph, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            raptor.run(source, target);
            csa.run(source, target);

            const RAPTOR::ProfileHandle csaProfile = csa.getProfileHandle(source);
            RAPTOR::ProfileHandle raptorProfile = raptor.getProfileHandle(target).reduceToEarliestArrivalTime();
            if (!profilesEqual(csaProfile, raptorProfile)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << std::endl;
                std::cout << "CSA profile:" << std::endl;
                std::cout << csaProfile << std::endl;
                std::cout << "RAPTOR profile:" << std::endl;
                std::cout << raptorProfile << std::endl;
            }
            progress++;
        }
    }
};

class ValidateUPProfileCSA : public ParameterizedCommand {

public:
    ValidateUPProfileCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUPProfileCSA", "Validates UP-Profile-CSA by comparing to Dijkstra-Profile-CSA on random queries.") {
        addParameter("Regular CSA data");
        addParameter("ULTRA CSA data");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Final transfers", { "Bucket", "PHAST" });
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
        addParameter("Sources", { "Vertices", "Stops" });
    }

    virtual void execute() noexcept {
        if (getParameter("Final transfers") == "Bucket") {
            chooseOrder<true>();
        } else {
            chooseOrder<false>();
        }
    }

private:
    template<bool USE_TARGET_BUCKETS>
    inline void chooseOrder() const noexcept {
        if (getParameter("Order") == "DFS") {
            run<USE_TARGET_BUCKETS, true>();
        } else {
            run<USE_TARGET_BUCKETS, false>();
        }
    }

    template<bool USE_TARGET_BUCKETS, bool USE_DFS_ORDER>
    inline void run() const noexcept {
        CSA::Data regularCSAData = CSA::Data::FromBinary(getParameter("Regular CSA data"));
        regularCSAData.sortConnectionsAscending();
        regularCSAData.printInfo();
        CSA::Data ultraCSAData = CSA::Data::FromBinary(getParameter("ULTRA CSA data"));
        ultraCSAData.sortConnectionsAscending();
        ultraCSAData.printInfo();
        CH::CH ch(getParameter("CH data"));

        IndexedSet<false, Vertex> sourceSet = getSourceSet(regularCSAData);

        const bool reorder = getParameter<bool>("Reorder network?");
        CSA::TransferGraph regularReverseGraph = regularCSAData.transferGraph;
        regularReverseGraph.revert();
        CSA::OneToAllDijkstraProfileCSA<CSA::SimpleProfiler, true> dijkstraCSA(regularCSAData, regularReverseGraph);
        CSA::TransferGraph ultraReverseGraph = ultraCSAData.transferGraph;
        ultraReverseGraph.revert();
        using UPProfileCSA = CSA::UPProfileCSA<USE_TARGET_BUCKETS, USE_DFS_ORDER, CSA::SimpleProfiler, true>;
        UPProfileCSA upCSA(ultraCSAData, ultraReverseGraph, ch, sourceSet, reorder);

        const Order order = reorder ? UPProfileCSA::vertexOrder(ch).splitAt(regularCSAData.numberOfStops()) : Order(Construct::Id, ch.numVertices());
        const Permutation permutation(Construct::Invert, order);

        srand(42);
        const size_t numTargets = (getParameter("Sources") == "Stops") ? regularCSAData.numberOfStops() : ch.numVertices();
        for (size_t i = 0; i < getParameter<size_t>("Number of queries"); i++) {
            const Vertex target(rand() % numTargets);
            dijkstraCSA.run(target);
            upCSA.run(permutation.permutate(target));
            for (const Vertex source : sourceSet) {
                if (regularCSAData.isStop(source)) continue;
                const RAPTOR::ProfileHandle dijkstraProfile = dijkstraCSA.getProfileHandle(source);
                const RAPTOR::ProfileHandle ultraProfile = upCSA.getProfileHandle(permutation.permutate(source));
                if (!Vector::equals(dijkstraProfile, ultraProfile)) {
                    std::cout << "Query " << i << ": Profiles not identical! " << source << " -> " << target << " @ " << std::endl;
                    std::cout << "Dijkstra profile:" << std::endl;
                    std::cout << dijkstraProfile << std::endl;
                    std::cout << "ULTRA-PHAST profile:" << std::endl;
                    std::cout << ultraProfile << std::endl;
                    std::cout << std::endl;
                    return;
                }
            }
        }
    }

    inline IndexedSet<false, Vertex> getSourceSet(const CSA::Data& data) const noexcept {
        if (getParameter("Sources") == "Stops") {
            IndexedSet<false, Vertex> result(data.transferGraph.numVertices());
            for (const StopId stop : data.stops()) {
                result.insert(stop);
            }
            return result;
        } else {
            return IndexedSet<false, Vertex>(Construct::Complete, data.transferGraph.numVertices());
        }
    }
};

class ValidateParetoULTRACSA : public ParameterizedCommand {

public:
    ValidateParetoULTRACSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateParetoULTRACSA", "Validates journeys computed by Pareto-ULTRA-CSA by comparing to ULTRA-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("CH data");
        addParameter("Number of trips", { "8", "16" });
    }

    virtual void execute() noexcept {
        if (getParameter<size_t>("Number of trips") == 8) {
            run<8>();
        } else {
            run<16>();
        }
    }

private:
    template<size_t MAX_TRIPS>
    inline void run() const noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        raptorData.printInfo();

        CH::CH ch(getParameter("CH data"));
        RAPTOR::ULTRARAPTOR<RAPTOR::AggregateProfiler, false> raptor(raptorData, ch);
        CSA::ParetoULTRACSA<true, CSA::AggregateProfiler, MAX_TRIPS> csa(csaData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, target, MAX_TRIPS - 1);
            csa.run(source, departureTime, target);
            std::vector<int> csaArrivalTimes = csa.getArrivalTimes(target);
            std::vector<int> raptorArrivalTimes = raptor.getArrivalTimes(target);
            if (!arrivalsEqual(csaArrivalTimes, raptorArrivalTimes)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "CSA arrivals: ";
                Vector::printConcise(csaArrivalTimes);
                std::cout << std::endl;
                std::cout << "RAPTOR arrivals: ";
                Vector::printConcise(raptorArrivalTimes);
                std::cout << std::endl;
                std::cout << "CSA journeys:" << std::endl;
                std::cout << csa.getJourneys(target) << std::endl;
                std::cout << "RAPTOR journeys:" << std::endl;
                std::cout << raptor.getJourneys(target) << std::endl;
                return;
            }
            progress++;
        }
        raptor.getProfiler().printStatistics();
        csa.getProfiler().printStatistics();
    }

    inline static bool arrivalsEqual(const std::vector<int>& csaArrivalTimes, const std::vector<int>& raptorArrivalTimes) noexcept {
        if (csaArrivalTimes.size() < raptorArrivalTimes.size()) return false;
        int minimum = never;
        for (size_t i = 0; i < raptorArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] != raptorArrivalTimes[i]) return false;
            minimum = std::min(csaArrivalTimes[i], minimum);
        }
        for (size_t i = raptorArrivalTimes.size(); i < csaArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] < minimum) return false;
        }
        return true;
    }
};

class ValidateParetoUPCSA : public ParameterizedCommand {

public:
    ValidateParetoUPCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "validateParetoUPCSA", "Validates journeys computed by Pareto-UP-CSA by comparing to UP-RAPTOR on random queries.") {
        addParameter("CSA input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
        addParameter("CH data");
        addParameter("Number of trips", { "8", "16" });
    }

    virtual void execute() noexcept {
        if (getParameter<size_t>("Number of trips") == 8) {
            run<8>();
        } else {
            run<16>();
        }
    }

private:
    template<size_t MAX_TRIPS>
    inline void run() const noexcept {
        CSA::Data csaData(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        RAPTOR::Data raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        IndexedSet<false, Vertex> targetSet(Construct::Complete, ch.numVertices());

        using UPRAPTOR = RAPTOR::UPRAPTOR<8, RAPTOR::NoProfiler>;
        UPRAPTOR raptor(raptorData, raptorData.transferGraph, ch, targetSet, false, true);
        using UPCSA = CSA::ParetoUPCSA<true, true, CSA::NoProfiler, MAX_TRIPS>;
        UPCSA csa(csaData, ch, targetSet, false);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            raptor.run(source, departureTime, MAX_TRIPS - 1);
            csa.run(source, departureTime);
            for (const StopId target : csaData.stops()) {
                std::vector<int> csaArrivalTimes = csa.getArrivalTimes(target);
                std::vector<int> raptorArrivalTimes = raptor.getArrivalTimes(target);
                if (!arrivalsEqual(csaArrivalTimes, raptorArrivalTimes)) {
                    std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "CSA arrivals: ";
                    Vector::printConcise(csaArrivalTimes);
                    std::cout << std::endl;
                    std::cout << "RAPTOR arrivals: ";
                    Vector::printConcise(raptorArrivalTimes);
                    std::cout << std::endl;
                    std::cout << "CSA journeys:" << std::endl;
                    std::cout << csa.getJourneys(target) << std::endl;
                    std::cout << "RAPTOR journeys:" << std::endl;
                    std::cout << raptor.getJourneys(target) << std::endl;
                    return;
                }
            }
            progress++;
        }
    }

    inline static bool arrivalsEqual(const std::vector<int>& csaArrivalTimes, const std::vector<int>& raptorArrivalTimes) noexcept {
        if (csaArrivalTimes.size() < raptorArrivalTimes.size()) return false;
        int minimum = never;
        for (size_t i = 0; i < raptorArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] != raptorArrivalTimes[i]) return false;
            minimum = std::min(csaArrivalTimes[i], minimum);
        }
        for (size_t i = raptorArrivalTimes.size(); i < csaArrivalTimes.size(); i++) {
            if (csaArrivalTimes[i] < minimum) return false;
        }
        return true;
    }
};

class ValidateTransitiveTripBased : public ParameterizedCommand {

public:
    ValidateTransitiveTripBased(BasicShell& shell) :
        ParameterizedCommand(shell, "validateTransitiveTripBased", "Validates journeys computed by TripBased by comparing to RAPTOR on random queries.") {
        addParameter("TripBased input file");
        addParameter("RAPTOR input file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        TripBased::Data tripBasedData(getParameter("TripBased input file"));
        tripBasedData.printInfo();
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();

        TripBased::TransitiveQuery<TripBased::NoProfiler> tripBased(tripBasedData);
        RAPTOR::RAPTOR<false, RAPTOR::NoProfiler, true, false, false> raptor(raptorData);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        for (size_t i = 0; i < n; i++) {
            const StopId source(rand() % raptorData.numberOfStops());
            const StopId target(rand() % raptorData.numberOfStops());
            const int departureTime(rand() % (24 * 60 * 60));
            tripBased.run(source, departureTime, target);
            raptor.run(source, departureTime, target);
            const std::vector<RAPTOR::ArrivalLabel> raptorArrivals = raptor.getArrivals();
            const std::vector<RAPTOR::ArrivalLabel> tripBasedArrivals = tripBased.getArrivals();
            if (!Vector::equals(raptorArrivals, tripBasedArrivals)) {
                std::cout << "Query " << i << ": " << source << " -> " << target << " @ " << departureTime << std::endl;
                std::cout << "RAPTOR arrivals:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : raptorArrivals) {
                    std::cout << label << std::endl;
                }
                std::cout << std::endl;
                std::cout << "TripBased arrivals:" << std::endl;
                for (const RAPTOR::ArrivalLabel& label : tripBasedArrivals) {
                    std::cout << label << std::endl;
                }
                std::cout << "RAPTOR journeys:" << std::endl;
                for (const RAPTOR::Journey& journey : raptor.getJourneys()) {
                    std::cout << journey << std::endl;
                }
                std::cout << "TripBased journeys:" << std::endl;
                for (const RAPTOR::Journey& journey : tripBased.getJourneys()) {
                    std::cout << journey << std::endl;
                }
            }
            progress++;
        }
    }
};

class ValidateUPTB : public ParameterizedCommand {

public:
    ValidateUPTB(BasicShell& shell) :
        ParameterizedCommand(shell, "validateUPTB", "Validates UP-TB by comparing to Dijkstra-RAPTOR on random queries.") {
        addParameter("RAPTOR data");
        addParameter("ULTRA-TB data");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Targets", { "Vertices", "Stops" });
        addParameter("Reorder network?");
        addParameter("Grouped rounds");
    }

    virtual void execute() noexcept {
        switch(getParameter<size_t>("Grouped rounds")) {
            case 0:
                run<0>();
                break;
            case 4:
                run<4>();
                break;
            case 6:
                run<6>();
                break;
            case 8:
                run<8>();
                break;
        }
    }

private:
    template<size_t GROUPED_ROUNDS>
    inline void run() const noexcept {
        RAPTOR::Data raptorData(getParameter("RAPTOR data"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        TripBased::Data tripBasedData(getParameter("ULTRA-TB data"));
        tripBasedData.printInfo();
        CH::CH ch(getParameter("CH data"));

        IndexedSet<false, Vertex> targetSet = getTargetSet(raptorData);

        RAPTOR::Data reverseRaptorData = raptorData.reverseNetwork();
        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::AggregateProfiler> dijkstraRaptor(raptorData, raptorData.transferGraph, reverseRaptorData.transferGraph);
        using UPTB = TripBased::UPQuery<GROUPED_ROUNDS, TripBased::AggregateProfiler>;
        const bool reorder = getParameter<bool>("Reorder network?");
        UPTB upTB(tripBasedData, raptorData.transferGraph, ch, targetSet, reorder);

        srand(42);
        for (size_t i = 0; i < getParameter<size_t>("Number of queries"); i++) {
            const Vertex source(rand() % raptorData.transferGraph.numVertices());
            const int departureTime(rand() % 24 * 60 * 60);
            dijkstraRaptor.run(source, departureTime);
            upTB.run(source, departureTime);
            for (const Vertex target : targetSet) {
                std::vector<RAPTOR::ArrivalLabel> dijkstraArrivals = dijkstraRaptor.getArrivals(target);
                std::vector<RAPTOR::ArrivalLabel> ultraArrivals = upTB.getArrivals(target);
                if (!Vector::equals(dijkstraArrivals, ultraArrivals)) {
                    std::cout << "Query " << i << ": Arrivals not identical! " << source << " -> " << target << " @ " << departureTime << std::endl;
                    std::cout << "Dijkstra journeys:" << std::endl;
                    std::cout << dijkstraRaptor.getJourneys(target);
                    std::cout << "ULTRA-PHAST journeys:" << std::endl;
                    std::cout << upTB.getJourneys(target);
                    std::cout << std::endl;
                    return;
                }
            }
        }
    }

    inline IndexedSet<false, Vertex> getTargetSet(const RAPTOR::Data& data) const noexcept {
        if (getParameter("Targets") == "Stops") {
            IndexedSet<false, Vertex> result(data.transferGraph.numVertices());
            for (const StopId stop : data.stops()) {
                result.insert(stop);
            }
            return result;
        } else {
            return IndexedSet<false, Vertex>(Construct::Complete, data.transferGraph.numVertices());
        }
    }
};

class ComputeBoundedMcRAPTORCoverage : public ParameterizedCommand {

public:
    ComputeBoundedMcRAPTORCoverage(BasicShell& shell) :
        ParameterizedCommand(shell, "computeBoundedMcRAPTORCoverage", "Computes coverage of Bounded McRAPTOR (ULTRA or transitive) results compared to ULTRA-McRAPTOR.") {
        addParameter("RAPTOR input file");
        addParameter("Limited RAPTOR input file");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Arrival slack");
        addParameter("Trip slack");
        addParameter("Transfer score", "900");
        addParameter("Walking factor", "2");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        RAPTOR::Data limitedRaptorData = RAPTOR::Data::FromBinary(getParameter("Limited RAPTOR input file"));
        limitedRaptorData.useImplicitDepartureBufferTimes();
        limitedRaptorData.printInfo();
        CH::CH ch(getParameter("CH data"));
        const RAPTOR::Data reverseData = limitedRaptorData.reverseNetwork();

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");
        transferScore = getParameter<double>("Transfer score");
        walkingScore = 1/getParameter<double>("Walking factor");

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> mcRaptor(raptorData, ch);
        RAPTOR::UBMRAPTOR<RAPTOR::AggregateProfiler> boundedMcRaptor(limitedRaptorData, reverseData, ch);

        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);
        Progress progress(n);
        double topKCoverage = 0;
        double fullCoverage = 0;
        double dominationCoverage = 0;
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            const Vertex target(rand() % ch.numVertices());
            const int departureTime(rand() % (24 * 60 * 60));
            mcRaptor.run(source, departureTime, target);
            std::vector<RAPTOR::WalkingParetoLabel> fullResults = mcRaptor.getResults();
            sort(fullResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> fullBag(fullResults);
            boundedMcRaptor.run(source, departureTime, target, arrivalSlack, tripSlack);
            std::vector<RAPTOR::WalkingParetoLabel> boundedResults = boundedMcRaptor.getResults();
            sort(boundedResults);
            RAPTOR::Bag<RAPTOR::WalkingParetoLabel> boundedBag(boundedResults);
            topKCoverage += computeTopKFuzzyCoverage(fullResults, boundedResults);
            fullCoverage += computeFullFuzzyCoverage(fullResults, boundedResults);
            dominationCoverage += computeDominationCoverage(fullResults, boundedResults, departureTime);
            progress++;
        }
        topKCoverage /= n;
        fullCoverage /= n;
        dominationCoverage /= n;
        std::cout << "Top " << K << " coverage: " << String::percent(topKCoverage) << std::endl;
        std::cout << "Full coverage: " << String::percent(fullCoverage) << std::endl;
        std::cout << "Domination coverage: " << String::percent(dominationCoverage) << std::endl;
        mcRaptor.getProfiler().printStatistics();
        boundedMcRaptor.getProfiler().printStatistics();
    }

private:
    double transferScore;
    double walkingScore;
    const RAPTOR::FuzzyComparison combinedEvaluation = RAPTOR::FuzzyComparison(180, 0.8);

    inline double getWeightedCost(const RAPTOR::WalkingParetoLabel& label, const int departureTime) const noexcept {
        return (label.arrivalTime - departureTime) + label.numberOfTrips * transferScore + label.walkingDistance * walkingScore;
    }

    inline double computeDominationCoverage(const std::vector<RAPTOR::WalkingParetoLabel>& fullLabels, const std::vector<RAPTOR::WalkingParetoLabel>& boundedLabels, const int departureTime) const noexcept {
        double bestBoundedCost = INFTY;
        double bestUnboundedCost = INFTY;
        for (const RAPTOR::WalkingParetoLabel& label : fullLabels) {
            const double cost = getWeightedCost(label, departureTime);
            if (Vector::contains(boundedLabels, label)) {
                bestBoundedCost = std::min(bestBoundedCost, cost);
            } else {
                bestUnboundedCost = std::min(bestUnboundedCost, cost);
            }
        }
        return combinedEvaluation(bestBoundedCost, bestUnboundedCost).lt();
    }
};
