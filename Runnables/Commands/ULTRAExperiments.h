#pragma once

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../../Algorithms/StronglyConnectedComponents.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"
#include "../../Algorithms/CH/CH.h"
#include "../../Algorithms/CSA/CSA.h"
#include "../../Algorithms/CSA/DijkstraCSA.h"
#include "../../Algorithms/CSA/ULTRACSA.h"
#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/InitialTransfers.h"
#include "../../Algorithms/RAPTOR/MCR.h"
#include "../../Algorithms/RAPTOR/RAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRAMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRA/Builder.h"
#include "../../Algorithms/RAPTOR/ULTRA/McBuilder.h"
#include "../../Algorithms/TripBased/Preprocessing/McULTRABuilder.h"
#include "../../Algorithms/TripBased/Preprocessing/ULTRABuilder.h"
#include "../../Algorithms/TripBased/Query/McQuery.h"
#include "../../Algorithms/TripBased/Query/Query.h"
#include "../../Algorithms/TripBased/Query/TransitiveQuery.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/Queries/Queries.h"

#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/File.h"

#include "../../Shell/Shell.h"

using namespace Shell;

using Profiler = CH::FullProfiler;
using WitnessSearch = CH::BidirectionalWitnessSearch<CHCoreGraph, Profiler, 200>;
inline static constexpr int ShortcutWeight = 1024;
inline static constexpr int LevelWeight = 256;
inline static constexpr int DegreeWeight = 0;
using RegularKeyFunction = CH::GreedyKey<WitnessSearch>;
using RegularStopCriterion = CH::NoStopCriterion;
using RegularCHBuilder = CH::Builder<Profiler, WitnessSearch, RegularKeyFunction, RegularStopCriterion, false, false>;
using CoreKeyFunction = CH::PartialKey<WitnessSearch, RegularKeyFunction>;
using CoreStopCriterion = CH::CoreCriterion;
using CoreCHBuilder = CH::Builder<Profiler, WitnessSearch, CoreKeyFunction, CoreStopCriterion, false, false>;

template<bool COUNT_OPTIMAL_CANDIDATES, bool IGNORE_ISOLATED_CANDIDATES>
using StopShortcutBuilder = RAPTOR::ULTRA::Builder<false, COUNT_OPTIMAL_CANDIDATES, IGNORE_ISOLATED_CANDIDATES>;
template<bool IGNORE_ISOLATED_CANDIDATES>
using EventShortcutBuilder = TripBased::ULTRABuilder<false, IGNORE_ISOLATED_CANDIDATES>;

template<typename PROFILER>
using TransitiveRAPTOR = RAPTOR::RAPTOR<true, PROFILER, true, false, false>;
template<typename PROFILER>
using DijkstraRAPTOR = RAPTOR::DijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, PROFILER, true, false, false>;
template<typename PROFILER>
using ULTRARAPTOR = RAPTOR::ULTRARAPTOR<PROFILER, false>;

template<typename PROFILER>
using TransitiveCSA = CSA::CSA<true, PROFILER>;
template<typename PROFILER>
using DijkstraCSA = CSA::DijkstraCSA<RAPTOR::CoreCHInitialTransfers, true, PROFILER>;
template<typename PROFILER>
using ULTRACSA = CSA::ULTRACSA<true, PROFILER>;

inline static constexpr size_t TimeFactor = 100;

struct CoreCHData {
    CH::CH ch;
    RAPTOR::Data coreData;
    double contractionTime;
};

inline static CoreCHData buildCoreCH(const RAPTOR::Data& raptorData, const size_t coreDegree) noexcept {
    CHCoreGraph graph;
    Graph::copy(raptorData.transferGraph, graph, Weight << TravelTime);
    const size_t numberOfStops = raptorData.numberOfStops();
    std::vector<bool> isNormalVertex(numberOfStops, false);
    isNormalVertex.resize(graph.numVertices(), true);

    CoreCHBuilder chBuilder(std::move(graph), CoreKeyFunction(isNormalVertex, graph.numVertices(), RegularKeyFunction(ShortcutWeight, LevelWeight, DegreeWeight)), CoreStopCriterion(numberOfStops, coreDegree));
    Timer timer;
    chBuilder.run();
    const double time = timer.elapsedMilliseconds();
    chBuilder.copyCoreToCH();
    CH::CH ch(std::move(chBuilder));

    Intermediate::TransferGraph coreGraph;
    coreGraph.addVertices(raptorData.transferGraph.numVertices());
    coreGraph[Coordinates] = raptorData.transferGraph[Coordinates];
    for (const Vertex vertex : coreGraph.vertices()) {
        if (ch.isCoreVertex(vertex)) {
            for (const Edge edge : ch.forward.edgesFrom(vertex)) {
                coreGraph.addEdge(vertex, ch.forward.get(ToVertex, edge)).set(TravelTime, ch.forward.get(Weight, edge));
            }
        }
    }

    CoreCHData result{ch, raptorData, time};
    Graph::move(std::move(coreGraph), result.coreData.transferGraph);
    return result;
}

inline static CH::CH buildCH(TravelTimeGraph&& graph) noexcept {
    RegularCHBuilder chBuilder(std::move(graph), graph[TravelTime]);
    chBuilder.run();
    chBuilder.copyCoreToCH();
    return CH::CH(std::move(chBuilder));
}

struct StopShortcutData {
    RAPTOR::Data raptorData;
    double preprocessingTime;

    inline size_t numShortcuts() const noexcept {
        return raptorData.transferGraph.numEdges();
    }
};

struct EventShortcutData {
    TripBased::Data tripBasedData;
    double preprocessingTime;

    inline size_t numShortcuts() const noexcept {
        return tripBasedData.stopEventGraph.numEdges();
    }
};

template<bool COUNT_OPTIMAL_CANDIDATES, bool IGNORE_ISOLATED_CANDIDATES = false>
inline static StopShortcutData computeStopToStopShortcuts(const RAPTOR::Data& raptorData, const size_t witnessLimit, const size_t numberOfThreads, const size_t pinMultiplier) noexcept {
    RAPTOR::Data shortcutData = raptorData;
    shortcutData.useImplicitDepartureBufferTimes();
    StopShortcutBuilder<COUNT_OPTIMAL_CANDIDATES, IGNORE_ISOLATED_CANDIDATES> shortcutGraphBuilder(shortcutData);
    Timer timer;
    shortcutGraphBuilder.computeShortcuts(ThreadPinning(numberOfThreads, pinMultiplier), witnessLimit);
    const double time = timer.elapsedMilliseconds();
    Graph::move(std::move(shortcutGraphBuilder.getShortcutGraph()), shortcutData.transferGraph);
    shortcutData.dontUseImplicitDepartureBufferTimes();
    return StopShortcutData{shortcutData, time};
}

template<bool IGNORE_ISOLATED_CANDIDATES = false>
inline static EventShortcutData computeEventToEventShortcuts(const RAPTOR::Data& raptorData, const size_t witnessLimit, const size_t numberOfThreads, const size_t pinMultiplier) noexcept {
    TripBased::Data tripBasedData(raptorData);
    EventShortcutBuilder<IGNORE_ISOLATED_CANDIDATES> shortcutGraphBuilder(tripBasedData);
    Timer timer;
    shortcutGraphBuilder.computeShortcuts(ThreadPinning(numberOfThreads, pinMultiplier), witnessLimit);
    const double time = timer.elapsedMilliseconds();
    Graph::move(std::move(shortcutGraphBuilder.getStopEventGraph()), tripBasedData.stopEventGraph);
    return EventShortcutData{tripBasedData, time};
}

inline static StopShortcutData computeStopToStopMcShortcuts(const RAPTOR::Data& raptorData, const size_t numberOfThreads, const size_t pinMultiplier) noexcept {
    RAPTOR::Data shortcutData = raptorData;
    shortcutData.useImplicitDepartureBufferTimes();
    RAPTOR::ULTRA::McBuilder<false, true, false> shortcutGraphBuilder(shortcutData);
    Timer timer;
    shortcutGraphBuilder.computeShortcuts(ThreadPinning(numberOfThreads, pinMultiplier), 0, 3600);
    const double time = timer.elapsedMilliseconds();
    Graph::move(std::move(shortcutGraphBuilder.getShortcutGraph()), shortcutData.transferGraph);
    shortcutData.dontUseImplicitDepartureBufferTimes();
    return StopShortcutData{shortcutData, time};
}

inline static EventShortcutData computeEventToEventMcShortcuts(const RAPTOR::Data& raptorData, const size_t numberOfThreads, const size_t pinMultiplier) noexcept {
    TripBased::Data tripBasedData(raptorData);
    TripBased::McULTRABuilder<false, true, false> shortcutGraphBuilder(tripBasedData);
    Timer timer;
    shortcutGraphBuilder.computeShortcuts(ThreadPinning(numberOfThreads, pinMultiplier), 0, 3600);
    const double time = timer.elapsedMilliseconds();
    Graph::move(std::move(shortcutGraphBuilder.getStopEventGraph()), tripBasedData.stopEventGraph);
    return EventShortcutData{tripBasedData, time};
}

struct RAPTORStatistics {
    inline static void printHeader(IO::OFStream& statistics, const std::string prefix = "") noexcept {
        statistics << "," << prefix << "ClearTime";
        statistics << "," << prefix << "InitTime";
        statistics << "," << prefix << "InitialTransfersTime";
        statistics << "," << prefix << "CollectRoutesTime";
        statistics << "," << prefix << "ScanRoutesTime";
        statistics << "," << prefix << "IntermediateTransfersTime";
        statistics << "," << prefix << "TotalTime";
        statistics << "," << prefix << "ScannedRoutes";
        statistics << "," << prefix << "ScannedRouteSegments";
        statistics << "," << prefix << "RelaxedEdges";
        statistics << "," << prefix << "UpdatedStopsByTrip";
        statistics << "," << prefix << "UpdatedStopsByTransfer";
    }

    inline static void printHeaderLine(IO::OFStream& statistics, const std::string& linePrefix, const std::string columnPrefix = "") noexcept {
        statistics << linePrefix;
        printHeader(statistics, columnPrefix);
        statistics << "\n";
    }

    template<typename PROFILER>
    inline static void print(const PROFILER& profiler, IO::OFStream& statistics) noexcept {
        statistics << "," << profiler.getExtraRoundTime(RAPTOR::EXTRA_ROUND_CLEAR);
        statistics << "," << profiler.getPhaseTime(RAPTOR::PHASE_INITIALIZATION);
        const double initialTransfersTime = profiler.getPhaseTimeInExtraRound(RAPTOR::PHASE_TRANSFERS, RAPTOR::EXTRA_ROUND_INITIALIZATION);
        statistics << "," << initialTransfersTime;
        statistics << "," << profiler.getPhaseTime(RAPTOR::PHASE_COLLECT);
        statistics << "," << profiler.getPhaseTime(RAPTOR::PHASE_SCAN);
        statistics << "," << profiler.getPhaseTime(RAPTOR::PHASE_TRANSFERS) - initialTransfersTime;
        statistics << "," << profiler.getTotalTime();
        statistics << "," << profiler.getMetric(RAPTOR::METRIC_ROUTES);
        statistics << "," << profiler.getMetric(RAPTOR::METRIC_ROUTE_SEGMENTS);
        statistics << "," << profiler.getMetric(RAPTOR::METRIC_EDGES);
        statistics << "," << profiler.getMetric(RAPTOR::METRIC_STOPS_BY_TRIP);
        statistics << "," << profiler.getMetric(RAPTOR::METRIC_STOPS_BY_TRANSFER);
    }

    template<typename PROFILER>
    inline static void printLine(const PROFILER& profiler, IO::OFStream& statistics, const std::string& prefix) noexcept {
        statistics << prefix;
        print(profiler, statistics);
        statistics << "\n";
    }
};

struct CSAStatistics {
    inline static void printHeader(IO::OFStream& statistics, const std::string prefix = "") noexcept {
        statistics << "," << prefix << "ClearTime";
        statistics << "," << prefix << "InitializationTime";
        statistics << "," << prefix << "ConnectionScanTime";
        statistics << "," << prefix << "TotalTime";
        statistics << "," << prefix << "ScannedConnections";
        statistics << "," << prefix << "UpdatedStopsByTransfer";
        statistics << "," << prefix << "RelaxedEdges";
    }

    inline static void printHeaderLine(IO::OFStream& statistics, const std::string& linePrefix, const std::string columnPrefix = "") noexcept {
        statistics << linePrefix;
        printHeader(statistics, columnPrefix);
        statistics << "\n";
    }

    template<typename PROFILER>
    inline static void print(const PROFILER& profiler, IO::OFStream& statistics) noexcept {
        statistics << "," << profiler.getPhaseTime(CSA::PHASE_CLEAR);
        statistics << "," << profiler.getPhaseTime(CSA::PHASE_INITIALIZATION);
        statistics << "," << profiler.getPhaseTime(CSA::PHASE_CONNECTION_SCAN);
        statistics << "," << profiler.getTotalTime();
        statistics << "," << profiler.getMetric(CSA::METRIC_CONNECTIONS);
        statistics << "," << profiler.getMetric(CSA::METRIC_STOPS_BY_TRANSFER);
        statistics << "," << profiler.getMetric(CSA::METRIC_EDGES);
    }

    template<typename PROFILER>
    inline static void printLine(const PROFILER& profiler, IO::OFStream& statistics, const std::string& prefix) noexcept {
        statistics << prefix;
        print(profiler, statistics);
        statistics << "\n";
    }
};

struct TBStatistics {
    inline static void printHeader(IO::OFStream& statistics, const std::string prefix = "") noexcept {
        statistics << "," << prefix << "ScanInitialTime";
        statistics << "," << prefix << "EvaluateInitialTime";
        statistics << "," << prefix << "ScanTime";
        statistics << "," << prefix << "TotalTime";
        statistics << "," << prefix << "ScannedTrips";
        statistics << "," << prefix << "ScannedStops";
        statistics << "," << prefix << "RelaxedTransfers";
        statistics << "," << prefix << "EnqueuedTrips";
    }

    inline static void printHeaderLine(IO::OFStream& statistics, const std::string& linePrefix, const std::string columnPrefix = "") noexcept {
        statistics << linePrefix;
        printHeader(statistics, columnPrefix);
        statistics << "\n";
    }

    template<typename PROFILER>
    inline static void print(const PROFILER& profiler, IO::OFStream& statistics) noexcept {
        statistics << "," << profiler.getPhaseTime(TripBased::PHASE_SCAN_INITIAL);
        statistics << "," << profiler.getPhaseTime(TripBased::PHASE_EVALUATE_INITIAL);
        statistics << "," << profiler.getPhaseTime(TripBased::PHASE_SCAN_TRIPS);
        statistics << "," << profiler.getTotalTime();
        statistics << "," << profiler.getMetric(TripBased::METRIC_SCANNED_TRIPS);
        statistics << "," << profiler.getMetric(TripBased::METRIC_SCANNED_STOPS);
        statistics << "," << profiler.getMetric(TripBased::METRIC_RELAXED_TRANSFERS);
        statistics << "," << profiler.getMetric(TripBased::METRIC_ENQUEUES);
    }

    template<typename PROFILER>
    inline static void printLine(const PROFILER& profiler, IO::OFStream& statistics, const std::string& prefix) noexcept {
        statistics << prefix;
        print(profiler, statistics);
        statistics << "\n";
    }
};

std::vector<size_t> witnessLimits{0, 30 * 60, 60 * 60, 2 * 60 * 60, 4 * 60 * 60, 48 * 60 * 60};
std::vector<size_t> transferSpeeds{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140};

class CoreDegreeWitnessLimitExperiment : public ParameterizedCommand {

public:
    CoreDegreeWitnessLimitExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "coreDegreeWitnessLimitExperiment", "Tests different configurations of contraction degree and witness limit for ULTRA preprocessing.") {
        addParameter("Directory");
        addParameter("Min core degree");
        addParameter("Max core degree");
        addParameter("Core degree step");
    }

    virtual void execute() noexcept {
        const std::string directory = getParameter("Directory");
        const RAPTOR::Data originalData(directory + "Full/raptor.binary");
        const size_t minCoreDegree = getParameter<size_t>("Min core degree");
        const size_t maxCoreDegree = getParameter<size_t>("Max core degree");
        const size_t coreDegreeStep = getParameter<size_t>("Core degree step");

        //Warmup
        const CoreCHData contractionData = buildCoreCH(originalData, minCoreDegree);
        computeStopToStopShortcuts<false>(contractionData.coreData, 0, numberOfCores(), 1);

        IO::OFStream fullStatistics(directory + "Shortcuts/preprocessing_experiment.csv");
        IO::OFStream stopStatistics(directory + "Shortcuts/preprocessing_experiment_stop.csv");
        IO::OFStream eventStatistics(directory + "Shortcuts/preprocessing_experiment_event.csv");
        fullStatistics << "CoreDegree,WitnessLimit,CHTime,StopShortcutTime,StopTotalTime,StopShortcuts,EventShortcutTime,EventTotalTime,EventShortcuts\n";
        stopStatistics << "CoreDegree";
        eventStatistics << "CoreDegree";
        for (const size_t witnessLimit : witnessLimits) {
            stopStatistics << "\tTime" << witnessLimit << "\tShortcuts" << witnessLimit;
            eventStatistics << "\tTime" << witnessLimit << "\tShortcuts" << witnessLimit;
        }
        stopStatistics << "\n";
        eventStatistics << "\n";
        for (size_t coreDegree = minCoreDegree; coreDegree <= maxCoreDegree; coreDegree += coreDegreeStep) {
            stopStatistics << coreDegree;
            eventStatistics << coreDegree;
            const CoreCHData contractionData = buildCoreCH(originalData, coreDegree);
            const std::string chDirectory = directory + "Contracted/Contracted_" + std::to_string(coreDegree) + "/";
            contractionData.ch.writeBinary(chDirectory + "ch");
            contractionData.coreData.serialize(chDirectory + "raptor.binary");
            for (const size_t witnessLimit : witnessLimits) {
                const StopShortcutData stopShortcutData = computeStopToStopShortcuts<false>(contractionData.coreData, witnessLimit, numberOfCores(), 1);
                stopShortcutData.raptorData.serialize(directory + "Shortcuts/Shortcuts_" + std::to_string(coreDegree) + "_" + std::to_string(witnessLimit) + "/raptor.binary");
                const double stopTotalTime = contractionData.contractionTime + stopShortcutData.preprocessingTime;
                const size_t stopShortcuts = stopShortcutData.numShortcuts();
                stopStatistics << "\t" << stopTotalTime / 1000 << "\t" << static_cast<double>(stopShortcuts) / 1000;

                const EventShortcutData eventShortcutData = computeEventToEventShortcuts(contractionData.coreData, witnessLimit, numberOfCores(), 1);
                eventShortcutData.tripBasedData.serialize(directory + "Shortcuts/Shortcuts_" + std::to_string(coreDegree) + "_" + std::to_string(witnessLimit) + "/tripBased.binary");
                const double eventTotalTime = contractionData.contractionTime + eventShortcutData.preprocessingTime;
                const size_t eventShortcuts = eventShortcutData.numShortcuts();
                eventStatistics << "\t" << eventTotalTime / 1000 << "\t" << static_cast<double>(eventShortcuts) / 1000;

                fullStatistics << coreDegree << "," << witnessLimit << "," << contractionData.contractionTime;
                fullStatistics << "," << stopShortcutData.preprocessingTime << "," << stopTotalTime << "," << stopShortcuts;
                fullStatistics << "," << eventShortcutData.preprocessingTime << "," << eventTotalTime << "," << eventShortcuts << "\n";
                fullStatistics.flush();
            }
            stopStatistics << "\n";
            stopStatistics.flush();
            eventStatistics << "\n";
            eventStatistics.flush();
        }
    }
};

class ParallelizationExperiment : public ParameterizedCommand {

public:
    ParallelizationExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "parallelizationExperiment", "Measures ULTRA preprocessing time depending on number of threads used.") {
        addParameter("RAPTOR network");
        addParameter("Core degree");
        addParameter("Witness limit");
        addParameter("Max threads");
        addParameter("Count optimal candidates?");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const RAPTOR::Data originalData(getParameter("RAPTOR network"));
        const size_t coreDegree = getParameter<size_t>("Core degree");
        const size_t witnessLimit = getParameter<size_t>("Witness limit");
        const size_t maxThreads = getParameter<size_t>("Max threads");
        const bool countOptimalCandidates = getParameter<bool>("Count optimal candidates?");
        const std::string outputFile = getParameter("Output file");

        IO::OFStream statistics(outputFile);
        statistics << "Threads,CHTime,ShortcutTime,ShortcutSpeedup,TotalTime,TotalSpeedup\n";
        const CoreCHData contractionData = buildCoreCH(originalData, coreDegree);
        const double chTime = contractionData.contractionTime;
        const double sequentialShortcutTime = chooseCount(contractionData.coreData, witnessLimit, 1, 1, countOptimalCandidates).preprocessingTime;
        statistics << 1 << "," << chTime << "," << sequentialShortcutTime << ",1.0," << (chTime + sequentialShortcutTime) << ",1.0\n";
        for (size_t threads = 2; threads <= maxThreads; threads *= 2) {
            double shortcutTime = chooseCount(contractionData.coreData, witnessLimit, threads, numberOfCores() / threads, countOptimalCandidates).preprocessingTime;
            statistics << threads << "," << chTime << "," << shortcutTime << "," << (sequentialShortcutTime / shortcutTime) << "," << (chTime + shortcutTime) << "," << ((chTime + sequentialShortcutTime) / (chTime + shortcutTime)) << "\n";
            statistics.flush();
        }
    }

private:
    inline static StopShortcutData chooseCount(const RAPTOR::Data& raptorData, const size_t witnessLimit, const size_t numberOfThreads, const size_t threadPinning, const bool countOptimalCandidates) noexcept {
        if (countOptimalCandidates) {
            return computeStopToStopShortcuts<true>(raptorData, witnessLimit, numberOfThreads, threadPinning);
        } else {
            return computeStopToStopShortcuts<false>(raptorData, witnessLimit, numberOfThreads, threadPinning);
        }
    }
};

class TransferSpeedPreprocessingExperiment : public ParameterizedCommand {

public:
    TransferSpeedPreprocessingExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferSpeedPreprocessingExperiment", "Runs ULTRA preprocessing for different transfer speeds.") {
        addParameter("Directory");
    }

    virtual void execute() noexcept {
        const std::string directory = getParameter("Directory");
        const Intermediate::Data originalData(directory + "Original_OSM/intermediate.binary");
        std::vector<ExperimentData> data(transferSpeeds.size());

        const std::string experimentDirectory = directory + "SpeedExperiment/";
        IO::OFStream fullStatistics(experimentDirectory + "experiment.csv");
        fullStatistics << "Speed,ObeySpeedLimits,RequireDirectTransfer,CHTime,StopShortcutTime,StopTotalTime,StopShortcuts,EventShortcutTime,EventTotalTime,EventShortcuts\n";

        for (const bool obeySpeedLimits : {true, false}) {
            for (const bool ignoreIsolatedCandidates : {false, true}) {
                const std::string configurationDirectory = experimentDirectory + (obeySpeedLimits ? "Limited" : "Unlimited") + (ignoreIsolatedCandidates ? "Without" : "With") + "Isolated/";
                for (size_t i = 0; i < transferSpeeds.size(); i++) {
                    const size_t speed = transferSpeeds[i];
                    Intermediate::Data intermediateData = originalData;
                    intermediateData.scaleTimes(TimeFactor);
                    Graph::computeTravelTimes(intermediateData.transferGraph, speed, obeySpeedLimits, TimeFactor);
                    intermediateData.contractDegreeTwoVertices();
                    const RAPTOR::Data raptorData = RAPTOR::Data::FromIntermediate(intermediateData);
                    const std::string currentDirectory = configurationDirectory + "Speed_" + std::to_string(speed) + "/";
                    raptorData.serialize(currentDirectory + "Full/raptor.binary");
                    const CoreCHData contractionData = buildCoreCH(raptorData, CoreDegree);
                    contractionData.ch.writeBinary(currentDirectory + "Contracted/ch");
                    contractionData.coreData.serialize(currentDirectory + "Contracted/raptor.binary");

                    const StopShortcutData stopShortcutData = chooseIgnoreIsolatedStop(contractionData.coreData, ignoreIsolatedCandidates);
                    stopShortcutData.raptorData.serialize(currentDirectory + "Shortcuts/raptor.binary");
                    const double stopPreprocessingTime = contractionData.contractionTime + stopShortcutData.preprocessingTime;
                    data[i].stopTimes.emplace_back(stopPreprocessingTime / 1000);
                    data[i].stopShortcuts.emplace_back(static_cast<double>(stopShortcutData.numShortcuts()) / 1000);

                    const EventShortcutData eventShortcutData = chooseIgnoreIsolatedEvent(contractionData.coreData, ignoreIsolatedCandidates);
                    eventShortcutData.tripBasedData.serialize(currentDirectory + "Shortcuts/tripBased.binary");
                    const double eventPreprocessingTime = contractionData.contractionTime + eventShortcutData.preprocessingTime;
                    data[i].eventTimes.emplace_back(eventPreprocessingTime / 1000);
                    data[i].eventShortcuts.emplace_back(static_cast<double>(eventShortcutData.numShortcuts()) / 1000);

                    fullStatistics << speed << "," << obeySpeedLimits << "," << ignoreIsolatedCandidates << "," << contractionData.contractionTime;
                    fullStatistics << "," << stopShortcutData.preprocessingTime << "," << stopPreprocessingTime << "," << stopShortcutData.numShortcuts();
                    fullStatistics << "," << eventShortcutData.preprocessingTime << "," << eventPreprocessingTime << "," << eventShortcutData.numShortcuts() << "\n";
                    fullStatistics.flush();

                }
            }
        }

        IO::OFStream stopStatistics(experimentDirectory + "experiment_stop.csv");
        stopStatistics << "Speed\tTimeLimitedWithIsolated\tShortcutsLimitedWithIsolated\tTimeLimitedWithoutIsolated\tShortcutsLimitedWithoutIsolated\tTimeUnlimitedWithIsolated\tShortcutsUnlimitedWithIsolated\tTimeUnlimitedWithoutIsolated\tShortcutsUnlimitedWithoutIsolated\n";
        for (size_t i = 0; i < transferSpeeds.size(); i++) {
            stopStatistics << transferSpeeds[i];
            for (size_t j = 0; j < data[i].stopTimes.size(); j++) {
                stopStatistics << "\t" << data[i].stopTimes[j] << "\t" << data[i].stopShortcuts[j];
            }
            stopStatistics << "\n";
        }

        IO::OFStream eventStatistics(experimentDirectory + "experiment_event.csv");
        eventStatistics << "Speed\tTimeLimitedWithIsolated\tShortcutsLimitedWithIsolated\tTimeLimitedWithoutIsolated\tShortcutsLimitedWithoutIsolated\tTimeUnlimitedWithIsolated\tShortcutsUnlimitedWithIsolated\tTimeUnlimitedWithoutIsolated\tShortcutsUnlimitedWithoutIsolated\n";
        for (size_t i = 0; i < transferSpeeds.size(); i++) {
            eventStatistics << transferSpeeds[i];
            for (size_t j = 0; j < data[i].eventTimes.size(); j++) {
                eventStatistics << "\t" << data[i].eventTimes[j] << "\t" << data[i].eventShortcuts[j];
            }
            eventStatistics << "\n";
        }
    }

private:
    struct ExperimentData {
        std::vector<double> stopTimes;
        std::vector<double> stopShortcuts;
        std::vector<double> eventTimes;
        std::vector<double> eventShortcuts;
    };

    inline static constexpr size_t CoreDegree = 14;
    inline static constexpr size_t WitnessLimit = 0;

    inline static StopShortcutData chooseIgnoreIsolatedStop(const RAPTOR::Data& raptorData, const bool ignoreIsolatedCandidates) noexcept {
        if (ignoreIsolatedCandidates) {
            return computeStopToStopShortcuts<false, true>(raptorData, WitnessLimit, numberOfCores(), 1);
        } else {
            return computeStopToStopShortcuts<false, false>(raptorData, WitnessLimit, numberOfCores(), 1);
        }
    }

    inline static EventShortcutData chooseIgnoreIsolatedEvent(const RAPTOR::Data& raptorData, const bool ignoreIsolatedCandidates) noexcept {
        if (ignoreIsolatedCandidates) {
            return computeEventToEventShortcuts<true>(raptorData, WitnessLimit, numberOfCores(), 1);
        } else {
            return computeEventToEventShortcuts<false>(raptorData, WitnessLimit, numberOfCores(), 1);
        }
    }

};

class TransferSpeedMcPreprocessingExperiment : public ParameterizedCommand {

public:
    TransferSpeedMcPreprocessingExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferSpeedMcPreprocessingExperiment", "Runs McULTRA preprocessing for different transfer speeds.") {
        addParameter("Directory");
    }

    virtual void execute() noexcept {
        const std::string directory = getParameter("Directory");
        const Intermediate::Data originalData(directory + "Original_OSM/intermediate.binary");
        std::vector<ExperimentData> data(transferSpeeds.size());

        const std::string experimentDirectory = directory + "McSpeedExperiment/";
        IO::OFStream fullStatistics(experimentDirectory + "experiment.csv");
        fullStatistics << "Speed,ObeySpeedLimits,CHTime,StopShortcutTime,StopTotalTime,StopShortcuts,EventShortcutTime,EventTotalTime,EventShortcuts\n";

        for (const bool obeySpeedLimits : {true, false}) {
            const std::string configurationDirectory = experimentDirectory + (obeySpeedLimits ? "Limited" : "Unlimited") + ("With") + "Isolated/";
            for (size_t i = 0; i < transferSpeeds.size(); i++) {
                const size_t speed = transferSpeeds[i];
                Intermediate::Data intermediateData = originalData;
                intermediateData.scaleTimes(TimeFactor);
                Graph::computeTravelTimes(intermediateData.transferGraph, speed, obeySpeedLimits, TimeFactor);
                intermediateData.contractDegreeTwoVertices();
                const RAPTOR::Data raptorData = RAPTOR::Data::FromIntermediate(intermediateData);
                const std::string currentDirectory = configurationDirectory + "Speed_" + std::to_string(speed) + "/";
                raptorData.serialize(currentDirectory + "Full/raptor.binary");
                const CoreCHData contractionData = buildCoreCH(raptorData, CoreDegree);
                contractionData.ch.writeBinary(currentDirectory + "Contracted/ch");
                contractionData.coreData.serialize(currentDirectory + "Contracted/raptor.binary");

                const StopShortcutData stopShortcutData = computeStopToStopMcShortcuts(contractionData.coreData, numberOfCores(), 1);
                stopShortcutData.raptorData.serialize(currentDirectory + "Shortcuts/raptor.binary");
                const double stopPreprocessingTime = contractionData.contractionTime + stopShortcutData.preprocessingTime;
                data[i].stopTimes.emplace_back(stopPreprocessingTime / 1000);
                data[i].stopShortcuts.emplace_back(static_cast<double>(stopShortcutData.numShortcuts()) / 1000);

                const EventShortcutData eventShortcutData = computeEventToEventMcShortcuts(contractionData.coreData, numberOfCores(), 1);
                eventShortcutData.tripBasedData.serialize(currentDirectory + "Shortcuts/tripBased.binary");
                const double eventPreprocessingTime = contractionData.contractionTime + eventShortcutData.preprocessingTime;
                data[i].eventTimes.emplace_back(eventPreprocessingTime / 1000);
                data[i].eventShortcuts.emplace_back(static_cast<double>(eventShortcutData.numShortcuts()) / 1000);

                fullStatistics << speed << "," << obeySpeedLimits << "," << "," << contractionData.contractionTime;
                fullStatistics << "," << stopShortcutData.preprocessingTime << "," << stopPreprocessingTime << "," << stopShortcutData.numShortcuts();
                fullStatistics << "," << eventShortcutData.preprocessingTime << "," << eventPreprocessingTime << "," << eventShortcutData.numShortcuts() << "\n";
                fullStatistics.flush();
            }
        }

        IO::OFStream stopStatistics(experimentDirectory + "experiment_stop.csv");
        stopStatistics << "Speed\tTimeLimitedWithIsolated\tShortcutsLimitedWithIsolated\tTimeUnlimitedWithIsolated\tShortcutsUnlimitedWithIsolated\n";
        for (size_t i = 0; i < transferSpeeds.size(); i++) {
            stopStatistics << transferSpeeds[i];
            for (size_t j = 0; j < data[i].stopTimes.size(); j++) {
                stopStatistics << "\t" << data[i].stopTimes[j] << "\t" << data[i].stopShortcuts[j];
            }
            stopStatistics << "\n";
        }

        IO::OFStream eventStatistics(experimentDirectory + "experiment_event.csv");
        eventStatistics << "Speed\tTimeLimitedWithIsolated\tShortcutsLimitedWithIsolated\tTimeUnlimitedWithIsolated\tShortcutsUnlimitedWithIsolated\n";
        for (size_t i = 0; i < transferSpeeds.size(); i++) {
            eventStatistics << transferSpeeds[i];
            for (size_t j = 0; j < data[i].eventTimes.size(); j++) {
                eventStatistics << "\t" << data[i].eventTimes[j] << "\t" << data[i].eventShortcuts[j];
            }
            eventStatistics << "\n";
        }
    }

private:
    struct ExperimentData {
        std::vector<double> stopTimes;
        std::vector<double> stopShortcuts;
        std::vector<double> eventTimes;
        std::vector<double> eventShortcuts;
    };

    inline static constexpr size_t CoreDegree = 14;
};

class RAPTORQueryExperiment : public ParameterizedCommand {

public:
    RAPTORQueryExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "raptorQueryExperiment", "Runs ULTRA-RAPTOR, MR and transitive RAPTOR on random queries.") {
        addParameter("Transitive RAPTOR network");
        addParameter("Dijkstra RAPTOR network");
        addParameter("ULTRA-RAPTOR network");
        addParameter("Core CH data");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        RAPTOR::Data transitiveRaptorData(getParameter("Transitive RAPTOR network"));
        transitiveRaptorData.useImplicitDepartureBufferTimes();
        transitiveRaptorData.printInfo();

        RAPTOR::Data dijkstraRaptorData(getParameter("Dijkstra RAPTOR network"));
        dijkstraRaptorData.useImplicitDepartureBufferTimes();
        dijkstraRaptorData.printInfo();

        RAPTOR::Data ultraRaptorData(getParameter("ULTRA-RAPTOR network"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        ultraRaptorData.printInfo();

        const CH::CH coreCH(getParameter("Core CH data"));
        const CH::CH bucketCH(getParameter("Bucket CH data"));

        const size_t numQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");
        const std::vector<StopQuery> stopQueries = generateRandomStopQueries(transitiveRaptorData.numberOfStops(), numQueries, seed);
        const std::vector<VertexQuery> vertexQueries = generateRandomVertexQueries(dijkstraRaptorData.transferGraph.numVertices(), numQueries, seed);

        TransitiveRAPTOR<RAPTOR::AggregateProfiler> transitiveRaptor(transitiveRaptorData);
        DijkstraRAPTOR<RAPTOR::AggregateProfiler> dijkstraRaptor(dijkstraRaptorData, coreCH);
        ULTRARAPTOR<RAPTOR::AggregateProfiler> ultraRaptor(ultraRaptorData, bucketCH);

        for (const StopQuery& q : stopQueries) {
            transitiveRaptor.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            dijkstraRaptor.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            ultraRaptor.run(q.source, q.departureTime, q.target);
        }

        IO::OFStream statistics(getParameter("Output file"));
        RAPTORStatistics::printHeaderLine(statistics, "Algorithm");
        RAPTORStatistics::printLine(transitiveRaptor.getProfiler(), statistics, "Transitive");
        RAPTORStatistics::printLine(dijkstraRaptor.getProfiler(), statistics, "Dijkstra");
        RAPTORStatistics::printLine(ultraRaptor.getProfiler(), statistics, "ULTRA");
        statistics.flush();
    }
};

class CSAQueryExperiment : public ParameterizedCommand {

public:
    CSAQueryExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "csaQueryExperiment", "Runs ULTRA-CSA, MCSA and transitive CSA on random queries.") {
        addParameter("Transitive CSA network");
        addParameter("Dijkstra CSA network");
        addParameter("ULTRA-CSA network");
        addParameter("Core CH data");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        CSA::Data transitiveCSAData(getParameter("Transitive CSA network"));
        transitiveCSAData.sortConnectionsAscending();
        transitiveCSAData.printInfo();

        CSA::Data dijkstraCSAData(getParameter("Dijkstra CSA network"));
        dijkstraCSAData.sortConnectionsAscending();
        dijkstraCSAData.printInfo();

        CSA::Data ultraCSAData(getParameter("ULTRA-CSA network"));
        ultraCSAData.sortConnectionsAscending();
        ultraCSAData.printInfo();

        const CH::CH coreCH(getParameter("Core CH data"));
        const CH::CH bucketCH(getParameter("Bucket CH data"));

        const size_t numQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");
        const std::vector<StopQuery> stopQueries = generateRandomStopQueries(transitiveCSAData.numberOfStops(), numQueries, seed);
        const std::vector<VertexQuery> vertexQueries = generateRandomVertexQueries(dijkstraCSAData.transferGraph.numVertices(), numQueries, seed);

        TransitiveCSA<CSA::AggregateProfiler> transitiveCSA(transitiveCSAData);
        DijkstraCSA<CSA::AggregateProfiler> dijkstraCSA(dijkstraCSAData, coreCH);
        ULTRACSA<CSA::AggregateProfiler> ultraCSA(ultraCSAData, bucketCH);

        for (const StopQuery& q : stopQueries) {
            transitiveCSA.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            dijkstraCSA.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            ultraCSA.run(q.source, q.departureTime, q.target);
        }

        IO::OFStream statistics(getParameter("Output file"));
        CSAStatistics::printHeaderLine(statistics, "Algorithm");
        CSAStatistics::printLine(transitiveCSA.getProfiler(), statistics, "Transitive");
        CSAStatistics::printLine(dijkstraCSA.getProfiler(), statistics, "Dijkstra");
        CSAStatistics::printLine(ultraCSA.getProfiler(), statistics, "ULTRA");
        statistics.flush();
    }
};

class TBQueryExperiment : public ParameterizedCommand {

public:
    TBQueryExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "tbQueryExperiment", "Runs ULTRA-TB and transitive TB on random queries.") {
        addParameter("Transitive TB network");
        addParameter("Sequential ULTRA-TB network");
        addParameter("Integrated ULTRA-TB network");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const TripBased::Data transitiveTBData(getParameter("Transitive TB network"));
        transitiveTBData.printInfo();
        const TripBased::Data sequentialUltraTBData(getParameter("Sequential ULTRA-TB network"));
        sequentialUltraTBData.printInfo();
        const TripBased::Data integratedUltraTBData(getParameter("Integrated ULTRA-TB network"));
        sequentialUltraTBData.printInfo();
        const CH::CH bucketCH(getParameter("Bucket CH data"));

        const size_t numQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");
        const std::vector<StopQuery> stopQueries = generateRandomStopQueries(transitiveTBData.numberOfStops(), numQueries, seed);
        const std::vector<VertexQuery> vertexQueries = generateRandomVertexQueries(bucketCH.numVertices(), numQueries, seed);

        TripBased::TransitiveQuery<TripBased::AggregateProfiler> transitiveTB(transitiveTBData);
        TripBased::Query<TripBased::AggregateProfiler> sequentialUltraTB(sequentialUltraTBData, bucketCH);
        TripBased::Query<TripBased::AggregateProfiler> integratedUltraTB(integratedUltraTBData, bucketCH);

        for (const StopQuery& q : stopQueries) {
            transitiveTB.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            sequentialUltraTB.run(q.source, q.departureTime, q.target);
        }

        for (const VertexQuery& q : vertexQueries) {
            integratedUltraTB.run(q.source, q.departureTime, q.target);
        }

        IO::OFStream statistics(getParameter("Output file"));
        TBStatistics::printHeaderLine(statistics, "Algorithm");
        TBStatistics::printLine(transitiveTB.getProfiler(), statistics, "Transitive");
        TBStatistics::printLine(sequentialUltraTB.getProfiler(), statistics, "SequentialULTRA");
        TBStatistics::printLine(integratedUltraTB.getProfiler(), statistics, "IntegratedULTRA");
        statistics.flush();
    }
};

class TransferSpeedQueryExperiment : public ParameterizedCommand {

public:
    TransferSpeedQueryExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferSpeedQueryExperiment", "Compares Dijkstra RAPTOR, ULTRA-RAPTOR and ULTRA-TB on random queries for different transfer speeds.") {
        addParameter("Base directory");
        addParameter("Number of queries");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string baseDirectory = getParameter("Base directory");
        IO::OFStream statistics(baseDirectory + "speed_experiment.csv");
        statistics << "TransferSpeed";
        RAPTORStatistics::printHeader(statistics, "DijkstraRAPTOR");
        RAPTORStatistics::printHeader(statistics, "ULTRARAPTOR");
        TBStatistics::printHeader(statistics, "ULTRATB");
        statistics << "\n";

        const RAPTOR::TransferGraph sampleGraph(baseDirectory + "Speed_1/Full/raptor.binary.graph");
        const size_t numberOfQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");
        const std::vector<VertexQuery> queries = generateRandomVertexQueries(sampleGraph.numVertices(), numberOfQueries, seed, TimeFactor);

        for (const size_t transferSpeed : transferSpeeds) {
            const std::string speedDirectory = baseDirectory + "Speed_" + std::to_string(transferSpeed) + "/";
            TravelTimeGraph transferGraph(speedDirectory + "Full/raptor.binary.graph");
            CH::CH bucketCH = buildCH(std::move(transferGraph));
            bucketCH.writeBinary(speedDirectory + "CH/ch");

            RAPTOR::Data dijkstraRaptorData(speedDirectory + "Contracted/raptor.binary");
            dijkstraRaptorData.useImplicitDepartureBufferTimes();
            CH::CH coreCH(speedDirectory + "Contracted/ch");
            RAPTOR::Data ultraRaptorData(speedDirectory + "Shortcuts/raptor.binary");
            ultraRaptorData.useImplicitDepartureBufferTimes();
            TripBased::Data ultraTBData(speedDirectory + "Shortcuts/tripBased.binary");

            DijkstraRAPTOR<RAPTOR::AggregateProfiler> dijkstraRaptor(dijkstraRaptorData, coreCH);
            ULTRARAPTOR<RAPTOR::AggregateProfiler> ultraRaptor(ultraRaptorData, bucketCH);
            TripBased::Query<TripBased::AggregateProfiler> ultraTB(ultraTBData, bucketCH);

            for (const VertexQuery& q : queries) {
                dijkstraRaptor.run(q.source, q.departureTime, q.target);
            }

            for (const VertexQuery& q : queries) {
                ultraRaptor.run(q.source, q.departureTime, q.target);
            }

            for (const VertexQuery& q : queries) {
                ultraTB.run(q.source, q.departureTime, q.target);
            }

            statistics << transferSpeed;
            RAPTORStatistics::print(dijkstraRaptor.getProfiler(), statistics);
            RAPTORStatistics::print(ultraRaptor.getProfiler(), statistics);
            TBStatistics::print(ultraTB.getProfiler(), statistics);
            statistics << "\n";
            statistics.flush();
        }
    }
};

class TransferSpeedMcQueryExperiment : public ParameterizedCommand {

public:
    TransferSpeedMcQueryExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferSpeedMcQueryExperiment", "Compares MCR, ULTRA-McRAPTOR and ULTRA-McTB on random queries for different transfer speeds.") {
        addParameter("Base directory");
        addParameter("Number of queries");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string baseDirectory = getParameter("Base directory");
        IO::OFStream statistics(baseDirectory + "mc_speed_experiment.csv");
        statistics << "TransferSpeed";
        RAPTORStatistics::printHeader(statistics, "MCR");
        RAPTORStatistics::printHeader(statistics, "ULTRAMcRAPTOR");
        TBStatistics::printHeader(statistics, "ULTRAMcTB");
        statistics << "\n";

        const RAPTOR::TransferGraph sampleGraph(baseDirectory + "Speed_1/Full/raptor.binary.graph");
        const size_t numberOfQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");
        const std::vector<VertexQuery> queries = generateRandomVertexQueries(sampleGraph.numVertices(), numberOfQueries, seed, TimeFactor);

        for (const size_t transferSpeed : transferSpeeds) {
            const std::string speedDirectory = baseDirectory + "Speed_" + std::to_string(transferSpeed) + "/";
            TravelTimeGraph transferGraph(speedDirectory + "Full/raptor.binary.graph");
            const CH::CH bucketCH = buildCH(std::move(transferGraph));
            bucketCH.writeBinary(speedDirectory + "CH/ch");

            RAPTOR::Data mcrData(speedDirectory + "Contracted/raptor.binary");
            mcrData.useImplicitDepartureBufferTimes();
            const CH::CH coreCH(speedDirectory + "Contracted/ch");
            RAPTOR::Data ultraRaptorData(speedDirectory + "Shortcuts/raptor.binary");
            ultraRaptorData.useImplicitDepartureBufferTimes();
            TripBased::Data ultraTBData(speedDirectory + "Shortcuts/tripBased.binary");

            RAPTOR::MCR<true, RAPTOR::AggregateProfiler> mcr(mcrData, coreCH);
            RAPTOR::ULTRAMcRAPTOR<RAPTOR::AggregateProfiler> ultraMcRaptor(ultraRaptorData, bucketCH);
            TripBased::McQuery<TripBased::AggregateProfiler> ultraMcTB(ultraTBData, bucketCH);

            for (const VertexQuery& q : queries) {
                mcr.run(q.source, q.departureTime, q.target);
            }

            for (const VertexQuery& q : queries) {
                ultraMcRaptor.run(q.source, q.departureTime, q.target);
            }

            for (const VertexQuery& q : queries) {
                ultraMcTB.run(q.source, q.departureTime, q.target);
            }

            statistics << transferSpeed;
            RAPTORStatistics::print(mcr.getProfiler(), statistics);
            RAPTORStatistics::print(ultraMcRaptor.getProfiler(), statistics);
            TBStatistics::print(ultraMcTB.getProfiler(), statistics);
            statistics << "\n";
            statistics.flush();
        }
    }
};

class TransferSpeedTravelTimeExperiment : public ParameterizedCommand {

public:
    TransferSpeedTravelTimeExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferSpeedTravelTimeExperiment", "Measures how the average travel time changes with transfer speed.") {
        addParameter("Base directory");
        addParameter("Number of queries");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const std::string baseDirectory = getParameter("Base directory");
        const std::string outputFile = baseDirectory + "travel_time_experiment.csv";
        const size_t numberOfQueries = getParameter<size_t>("Number of queries");
        const int seed = getParameter<int>("Seed");

        IO::OFStream statistics(outputFile);
        statistics << "TransferSpeed,BestTravelTime,NumberOfTrips,InitialTransferTime,IntermediateTransferTime,TotalTransferTime,DirectTransferTime\n";

        std::vector<VertexQuery> queries;
        const CH::CH tempCH(baseDirectory + "Speed_" + std::to_string(transferSpeeds[0]) + "/CH/ch");
        std::mt19937 randomGenerator(seed);
        std::uniform_int_distribution<> vertexDistribution(0, tempCH.numVertices() - 1);
        std::uniform_int_distribution<> timeDistribution(0, (24 * 60 * 60) - 1);
        RAPTOR::BucketCHInitialTransfers chQuery(tempCH);
        while (queries.size() < numberOfQueries) {
            const Vertex source = Vertex(vertexDistribution(randomGenerator));
            const Vertex target = Vertex(vertexDistribution(randomGenerator));
            chQuery.run(source, target);
            if (chQuery.reachable()) {
                const int departureTime = (timeDistribution(randomGenerator)) * TimeFactor;
                queries.emplace_back(source, target, departureTime);
            }
        }

        for (const size_t transferSpeed : transferSpeeds) {
            const std::string speedDirectory = baseDirectory + "Speed_" + std::to_string(transferSpeed) + "/";
            RAPTOR::Data shortcutRaptorData(speedDirectory + "Shortcuts/raptor.binary");
            shortcutRaptorData.useImplicitDepartureBufferTimes();
            const CH::CH bucketCH(speedDirectory + "CH/ch");

            ULTRARAPTOR<RAPTOR::NoProfiler> ultraRaptor(shortcutRaptorData, bucketCH);
            double averageTravelTime = 0;
            double averageNumberOfTrips = 0;
            double averageInitialTransferTime = 0;
            double averageIntermediateTransferTime = 0;
            double averageTotalTransferTime = 0;
            double averageDirectTransferTime = 0;
            for (const VertexQuery& q : queries) {
                ultraRaptor.run(q.source, q.departureTime, q.target);
                const int arrivalTime = ultraRaptor.getEarliestArrivalTime(q.target);
                averageTravelTime += (arrivalTime - q.departureTime) / TimeFactor;
                const RAPTOR::Journey journey = ultraRaptor.getEarliestJourney(q.target);
                averageNumberOfTrips += countTrips(journey);
                averageInitialTransferTime += initialTransferTime(journey) / TimeFactor;
                averageIntermediateTransferTime += intermediateTransferTime(journey) / TimeFactor;
                averageTotalTransferTime += totalTransferTime(journey) / TimeFactor;
                averageDirectTransferTime += ultraRaptor.getDirectTransferTime() / TimeFactor;
            }
            averageTravelTime /= numberOfQueries;
            averageNumberOfTrips /= numberOfQueries;
            averageInitialTransferTime /= numberOfQueries;
            averageIntermediateTransferTime /= numberOfQueries;
            averageTotalTransferTime /= numberOfQueries;
            averageDirectTransferTime /= numberOfQueries;
            statistics << transferSpeed << "," << averageTravelTime << "," << averageNumberOfTrips << "," << averageInitialTransferTime << "," << averageIntermediateTransferTime << "," << averageTotalTransferTime << "," << averageDirectTransferTime << "\n";
            statistics.flush();
        }
    }
};

class GeoRankExperiment : public ParameterizedCommand {

public:
    GeoRankExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "geoRankExperiment", "Compares Dijkstra RAPTOR, ULTRA-RAPTOR and ULTRA-TB on geo-rank queries.") {
        addParameter("Transfer graph");
        addParameter("Dijkstra RAPTOR network");
        addParameter("Core CH data");
        addParameter("ULTRA-RAPTOR network");
        addParameter("ULTRA-TB network");
        addParameter("Bucket CH data");
        addParameter("Number of queries");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        const RAPTOR::TransferGraph graph(getParameter("Transfer graph"));
        RAPTOR::Data dijkstraRaptorData(getParameter("Dijkstra RAPTOR network"));
        dijkstraRaptorData.useImplicitDepartureBufferTimes();
        const CH::CH coreCH(getParameter("Core CH data"));
        RAPTOR::Data ultraRaptorData(getParameter("ULTRA-RAPTOR network"));
        ultraRaptorData.useImplicitDepartureBufferTimes();
        const TripBased::Data ultraTBData(getParameter("ULTRA-TB network"));
        const CH::CH bucketCH(getParameter("Bucket CH data"));
        const size_t numberOfQueries = getParameter<size_t>("Number of queries");
        const std::string outputFile = getParameter("Output file");
        const int seed = getParameter<int>("Seed");

        std::mt19937 randomGenerator(seed);
        std::uniform_int_distribution<> vertexDistribution(0, graph.numVertices() - 1);
        std::uniform_int_distribution<> timeDistribution(0, (24 * 60 * 60) - 1);

        std::vector<RankQuery> queries;
        std::vector<TargetDistance> targets;
        size_t numRanks = 0;
        for (size_t i = 0; i < numberOfQueries; i++) {
            const Vertex source = Vertex(vertexDistribution(randomGenerator));
            const int departureTime = timeDistribution(randomGenerator);
            targets.clear();
            for (const Vertex vertex : graph.vertices()) {
                targets.emplace_back(vertex, Geometry::geoDistanceInCM(graph.get(Coordinates, source), graph.get(Coordinates, vertex)));
            }
            std::sort(targets.begin(), targets.end());
            size_t queryRank = 0;
            for (size_t rank = 1; rank < targets.size(); rank = rank * 2) {
                queries.emplace_back(RankQuery{queryRank, source, targets[rank].target, departureTime});
                queryRank++;
                numRanks = std::max(numRanks, queryRank);
            }
        }

        std::uniform_int_distribution<> queryDistribution(0, queries.size() - 1);
        for (size_t i = 0; i < queries.size(); i++) {
            const size_t j = queryDistribution(randomGenerator);
            if (i == j) continue;
            std::swap(queries[i], queries[j]);
        }

        DijkstraRAPTOR<RAPTOR::NoProfiler> dijkstraRaptor(dijkstraRaptorData, coreCH);
        ULTRARAPTOR<RAPTOR::NoProfiler> ultraRaptor(ultraRaptorData, bucketCH);
        TripBased::Query<TripBased::NoProfiler> ultraTB(ultraTBData, bucketCH);

        std::vector<std::vector<double>> dijkstraTimes(numRanks);
        Timer timer;
        for (const RankQuery& q : queries) {
            timer.restart();
            dijkstraRaptor.run(q.source, q.departureTime, q.target);
            dijkstraTimes[q.rank].emplace_back(timer.elapsedMicroseconds());
        }

        std::vector<std::vector<double>> ultraRaptorTimes(numRanks);
        for (const RankQuery& q : queries) {
            timer.restart();
            ultraRaptor.run(q.source, q.departureTime, q.target);
            ultraRaptorTimes[q.rank].emplace_back(timer.elapsedMicroseconds());
        }

        std::vector<std::vector<double>> ultraTBTimes(numRanks);
        for (const RankQuery& q : queries) {
            timer.restart();
            ultraTB.run(q.source, q.departureTime, q.target);
            ultraTBTimes[q.rank].emplace_back(timer.elapsedMicroseconds());
        }

        IO::OFStream dijkstraStatistics(outputFile + "_dijkstra.csv");
        IO::OFStream dijkstraFullStatistics(outputFile + "_dijkstra_full.csv");
        IO::OFStream ultraRaptorStatistics(outputFile + "_ultra_raptor.csv");
        IO::OFStream ultraRaptorFullStatistics(outputFile + "_ultra_raptor_full.csv");
        IO::OFStream ultraTBStatistics(outputFile + "_ultra_tb.csv");
        IO::OFStream ultraTBFullStatistics(outputFile + "_ultra_tb_full.csv");
        printStatistics(dijkstraTimes, dijkstraFullStatistics, dijkstraStatistics);
        printStatistics(ultraRaptorTimes, ultraRaptorFullStatistics, ultraRaptorStatistics);
        printStatistics(ultraTBTimes, ultraTBFullStatistics, ultraTBStatistics);
    }

private:
    struct RankQuery {
        size_t rank;
        Vertex source;
        Vertex target;
        int departureTime;
    };

    struct TargetDistance {
        TargetDistance(const Vertex target = noVertex, const double distance = 0) :
            target(target),
            distance(distance) {
        }

        inline bool operator<(const TargetDistance& other) const noexcept {
            return distance < other.distance;
        }

        Vertex target;
        double distance;
    };

    inline static void printStatistics(std::vector<std::vector<double>>& timesPerRank, IO::OFStream& fullStatistics, IO::OFStream& statistics) noexcept {
        size_t rank = 0;
        fullStatistics << "Rank\tQuery\tValue\n";
        statistics << "Rank\tMean\tMedian\tMin\tMax\tLowerQuartile\tUpperQuartile\tLowerWhisker\tUpperWhisker\n";
        for (std::vector<double>& values : timesPerRank) {
            std::sort(values.begin(), values.end());
            for (size_t i = 0; i < values.size(); i++) {
                fullStatistics << rank << "\t" << i << "\t" << values[i] << "\n";
            }
            fullStatistics.flush();
            const double mean = Vector::mean(values);
            const double median = Vector::median(values);
            const double min = Vector::min(values);
            const double max = Vector::max(values);
            const double lowerQuartile = Vector::percentile(values, 0.25);
            const double upperQuartile = Vector::percentile(values, 0.75);
            const double lowerWhiskerLimit = lowerQuartile - 1.5 * (upperQuartile - lowerQuartile);
            const double lowerWhisker = *std::lower_bound(values.begin(), values.end(), lowerWhiskerLimit);
            const double upperWhiskerLimit = upperQuartile + 1.5 * (upperQuartile - lowerQuartile);
            const double upperWhisker = *(std::upper_bound(values.begin(), values.end(), upperWhiskerLimit) - 1);
            statistics << rank;
            statistics << "\t" << mean / 1000;
            statistics << "\t" << median / 1000;
            statistics << "\t" << min / 1000;
            statistics << "\t" << max / 1000;
            statistics << "\t" << lowerQuartile / 1000;
            statistics << "\t" << upperQuartile / 1000;
            statistics << "\t" << lowerWhisker / 1000;
            statistics << "\t" << upperWhisker / 1000;
            statistics << "\n";
            statistics.flush();
            rank++;
        }
    }
};

class EdgeLengthHistogram : public ParameterizedCommand {

public:
    EdgeLengthHistogram(BasicShell& shell) :
        ParameterizedCommand(shell, "edgeLengthHistogram", "Creates histogram of edge lengths.") {
        addParameter("Graph binary");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const RAPTOR::TransferGraph graph(getParameter("Graph binary"));

        std::vector<int> travelTimes;
        for (const Edge edge : graph.edges()) {
            travelTimes.emplace_back(graph.get(TravelTime, edge));
        }
        std::sort(travelTimes.begin(), travelTimes.end());
        std::vector<size_t> histogram(1, 0);
        for (size_t i = 0; i < travelTimes.size(); i++) {
            while (travelTimes[i] >= (1 << histogram.size())) {
                histogram.emplace_back(0);
            }
            histogram.back()++;
        }

        IO::OFStream statistics(getParameter("Output file"));
        statistics << "EdgeLength,#Shortcuts\n";
        for (size_t i = 0; i < histogram.size(); i++) {
            statistics << (1 << i) << "," << histogram[i] << "\n";
        }
        statistics.flush();
    }
};

class EdgeGeoDistanceHistogram : public ParameterizedCommand {

public:
    EdgeGeoDistanceHistogram(BasicShell& shell) :
        ParameterizedCommand(shell, "edgeGeoDistanceHistogram", "Creates histogram of edge geo distances.") {
        addParameter("Graph binary");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const RAPTOR::TransferGraph graph(getParameter("Graph binary"));

        std::vector<double> geoDistances;
        for (const Vertex fromVertex : graph.vertices()) {
            for (const Edge edge : graph.edgesFrom(fromVertex)) {
                const Vertex toVertex = graph.get(ToVertex, edge);
                geoDistances.emplace_back(geoDistanceInCM(graph.get(Coordinates, fromVertex), graph.get(Coordinates, toVertex)) / 100);
            }
        }
        std::sort(geoDistances.begin(), geoDistances.end());
        std::vector<size_t> histogram(1, 0);
        for (size_t i = 0; i < geoDistances.size(); i++) {
            while (geoDistances[i] >= (1 << histogram.size())) {
                histogram.emplace_back(0);
            }
            histogram.back()++;
        }

        IO::OFStream statistics(getParameter("Output file"));
        statistics << "Distance,#Shortcuts\n";
        for (size_t i = 0; i < histogram.size(); i++) {
            statistics << (1 << i) << "," << histogram[i] << "\n";
        }
        statistics.flush();
    }
};

class StopDegreeHistogram : public ParameterizedCommand {

public:
    StopDegreeHistogram(BasicShell& shell) :
        ParameterizedCommand(shell, "stopDegreeHistogram", "Creates histogram of out-/in-degree of stops.") {
        addParameter("Graph binary");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph staticGraph(getParameter("Graph binary"));
        DynamicTransferGraph graph;
        Graph::move(std::move(staticGraph), graph);

        std::vector<size_t> outDegrees;
        std::vector<size_t> inDegrees;
        for (const Vertex vertex : graph.vertices()) {
            outDegrees.emplace_back(graph.outDegree(vertex));
            inDegrees.emplace_back(graph.inDegree(vertex));
        }
        std::sort(outDegrees.begin(), outDegrees.end());
        std::sort(inDegrees.begin(), inDegrees.end());

        std::vector<size_t> outDegreeHistogram(1, 0);
        std::vector<size_t> inDegreeHistogram(1, 0);
        for (size_t i = 0; i < outDegrees.size(); i++) {
            if (outDegrees[i] >= size_t(1 << outDegreeHistogram.size())) {
                outDegreeHistogram.emplace_back(0);
            }
            outDegreeHistogram.back()++;
            if (inDegrees[i] >= size_t(1 << inDegreeHistogram.size())) {
                inDegreeHistogram.emplace_back(0);
            }
            inDegreeHistogram.back()++;
        }

        IO::OFStream statistics(getParameter("Output file"));
        statistics << "Degree,OutVertices,InVertices\n";
        for (size_t i = 0; i < std::max(outDegreeHistogram.size(), inDegreeHistogram.size()); i++) {
            const size_t outVertices = (i < outDegreeHistogram.size()) ? outDegreeHistogram[i] : 0;
            const size_t inVertices = (i < inDegreeHistogram.size()) ? inDegreeHistogram[i] : 0;
            statistics << (1 << i) << "," << outVertices << "," << inVertices << "\n";
        }
        statistics.flush();
    }
};

class ComputeSizeOfTransitiveClosure : public ParameterizedCommand {

public:
    ComputeSizeOfTransitiveClosure(BasicShell& shell) :
        ParameterizedCommand(shell, "computeSizeOfTransitiveClosure", "Computes size of the graph's transitive closure and compares it to current size.") {
        addParameter("Graph binary");
    }

    virtual void execute() noexcept {
        const RAPTOR::TransferGraph graph(getParameter("Graph binary"));

        size_t transitiveEdges = 0;
        CondensationGraph condensation = condense(graph);
        Graph::makeTransitivelyClosed(condensation, TravelTime);
        for (const Vertex u : condensation.vertices()) {
            const size_t componentSize = condensation.get(Size, u);
            transitiveEdges += componentSize * (componentSize - 1);
            for (const Edge e : condensation.edgesFrom(u)) {
                const Vertex v = condensation.get(ToVertex, e);
                transitiveEdges += componentSize * condensation.get(Size, v);
            }
        }
        std::cout << "Transitive closure: " << transitiveEdges << std::endl;
        std::cout << "Original edges: " << graph.numEdges() << " (" << String::percent(static_cast<double>(graph.numEdges()) / transitiveEdges) << ")" << std::endl;
    }

private:
    inline static CondensationGraph condense(const RAPTOR::TransferGraph& graph) noexcept {
        CondensationGraph condensation;
        StronglyConnectedComponents<RAPTOR::TransferGraph> scc(graph);
        scc.run();
        const std::vector<size_t> componentSizes = scc.getComponentSizes();
        condensation.addVertices(scc.numComponents());
        for (const Vertex v : condensation.vertices()) {
            condensation.set(Size, v, componentSizes[v]);
        }
        for (const Vertex from : graph.vertices()) {
            for (const Edge edge : graph.edgesFrom(from)) {
                const Vertex to = graph.get(ToVertex, edge);
                const Vertex fromComponent(scc.getComponent(from));
                const Vertex toComponent(scc.getComponent(to));
                if (fromComponent != toComponent && !condensation.hasEdge(fromComponent, toComponent)) {
                    condensation.addEdge(fromComponent, toComponent);
                }
            }
        }
        return condensation;
    }
};

class ComputeTransitiveClosure : public ParameterizedCommand {

public:
    ComputeTransitiveClosure(BasicShell& shell) :
        ParameterizedCommand(shell, "computeTransitiveClosure", "Computes the transitive closure.") {
        addParameter("Graph binary");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Graph binary"));
        DynamicTransferGraph transitiveClosure;
        Graph::copy(graph, transitiveClosure);
        Graph::makeTransitivelyClosed(transitiveClosure, TravelTime);
        std::cout << "Transitive closure: " << transitiveClosure.numEdges() << std::endl;
        std::cout << "Original edges: " << graph.numEdges() << " (" << String::percent(static_cast<double>(graph.numEdges()) / transitiveClosure.numEdges()) << ")" << std::endl;
        transitiveClosure.writeBinary(getParameter("Output file"));
    }
};

class CompareGraphs : public ParameterizedCommand {

public:
    CompareGraphs(BasicShell& shell) :
        ParameterizedCommand(shell, "compareGraphs", "Compares the two given graphs.") {
        addParameter("Graph 1");
        addParameter("Graph 2");
    }

    virtual void execute() noexcept {
        const RAPTOR::TransferGraph graph1(getParameter("Graph 1"));
        graph1.printAnalysis();
        Graph::printInfo(graph1);
        StronglyConnectedComponents<RAPTOR::TransferGraph> components1(graph1);
        components1.run();
        std::cout << "Graph1:" << std::endl;
        std::cout << "\tNumber of components: " << components1.numComponents() << std::endl;
        std::cout << "\tSize of largest connected component: " << components1.maxComponentSize() << std::endl;

        const RAPTOR::TransferGraph graph2(getParameter("Graph 2"));
        graph2.printAnalysis();
        Graph::printInfo(graph2);
        StronglyConnectedComponents<RAPTOR::TransferGraph> components2(graph2);
        components2.run();
        std::cout << "Graph2:" << std::endl;
        std::cout << "\tNumber of components: " << components2.numComponents() << std::endl;
        std::cout << "\tSize of largest connected component: " << components2.maxComponentSize() << std::endl;

        size_t sharedEdges = 0;
        for (const Vertex u : graph1.vertices()) {
            for (const Edge e : graph1.edgesFrom(u)) {
                const Vertex v = graph1.get(ToVertex, e);
                const Edge f = graph2.findEdge(u, v);
                if (f != noEdge && graph2.get(TravelTime, f) == graph1.get(TravelTime, e)) sharedEdges++;
            }
        }

        std::cout << "Shared edges: " << sharedEdges << std::endl;
        std::cout << String::percent(sharedEdges/static_cast<double>(graph1.numEdges())) << " of edges in graph1" << std::endl;
        std::cout << String::percent(sharedEdges/static_cast<double>(graph2.numEdges())) << " of edges in graph2" << std::endl;
    }
};
