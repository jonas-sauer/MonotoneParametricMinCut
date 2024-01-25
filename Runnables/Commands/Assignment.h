#pragma once

#include "../../Algorithms/Assignment/Capacities/Assignment.h"
#include "../../Algorithms/Assignment/Capacities/ConnectionBasedLockstepAssignment.h"
#include "../../Algorithms/Assignment/Capacities/FixedAssignment.h"
#include "../../Algorithms/Assignment/Capacities/LockstepAssignment.h"
#include "../../Algorithms/Assignment/SequentialAssignment.h"

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Debug.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Calendar.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/File.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Ranges/Range.h"

#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"
#include "../../Helpers/String/TextFileUtils.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Vector/Permutation.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/IdVertexDemand.h"
#include "../../DataStructures/Demand/PassengerData.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Container/Set.h"

#include "../../Algorithms/UnionFind.h"
#include "../../Algorithms/Assignment/GroupAssignment.h"
#include "../../Algorithms/Assignment/MRUM/LogitNormal.h"
#include "../../Algorithms/Assignment/MRUM/Beta.h"
#include "../../Algorithms/Assignment/MRUM/PERT.h"
#include "../../Algorithms/Assignment/MRUM/Kumaraswamy.h"
#include "../../Algorithms/Assignment/MRUM/Assignment.h"
#include "../../Algorithms/Assignment/ULTRA/Assignment.h"

#include "../../Algorithms/Assignment/ComputePATs.h"
#include "../../Algorithms/Assignment/MRUM/SimpleSamplePATComputation.h"
#include "../../Algorithms/CH/Preprocessing/CHBuilder.h"
#include "../../Algorithms/CH/Preprocessing/Profiler.h"
#include "../../Algorithms/CH/Preprocessing/PathAwareKey.h"
#include "../../Algorithms/CH/Preprocessing/WitnessSearch.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"
#include "../../Algorithms/CH/Preprocessing/ComparingWitnessSearch.h"
#include "../../Algorithms/CH/Preprocessing/StopCriterion.h"
#include "../../Algorithms/CH/Query/CHQuery.h"
#include "../../Algorithms/CH/CHUtils.h"
#include "../../Algorithms/CH/CH.h"
#include "../../Algorithms/CSA/Profiler.h"
#include "../../Algorithms/Dijkstra/Dijkstra.h"

#include "../../Algorithms/DecisionModels/Kirchhoff.h"
#include "../../Algorithms/DecisionModels/Linear.h"
#include "../../Algorithms/DecisionModels/Logit.h"
#include "../../Algorithms/DecisionModels/Optimal.h"
#include "../../Algorithms/DecisionModels/RelativeLogit.h"

#include "../../Visualization/MapVisualization.h"
#include "../../Visualization/PDF.h"
#include "../../Visualization/Borders.h"
#include "../../Visualization/Color.h"

using namespace Shell;

class DemandHistogram : public ParameterizedCommand {

public:
    DemandHistogram(BasicShell& shell) :
        ParameterizedCommand(shell, "demandHistogram", "Computes a histogram for the origin-destination relation.") {
        addParameter("CSA binary");
        addParameter("Demand file");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand demand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph);
        demand.sortByDestination();
        Vertex destination = Vertex(0);
        std::vector<size_t> numberOfOriginsHistogram;
        IndexedSet<true, Vertex> origins;
        for (const AccumulatedVertexDemand::Entry& entry : demand.entries) {
            if ((entry.destinationVertex != destination) && (!origins.empty())) {
                if (numberOfOriginsHistogram.size() <= origins.size()) {
                    numberOfOriginsHistogram.resize(origins.size() + 1, 0);
                }
                numberOfOriginsHistogram[origins.size()]++;
                origins.clear();
                destination = entry.destinationVertex;
            }
            origins.insert(entry.originVertex);
        }
        std::cout << "NumberOfOrigins Histogram." << std::endl;
        for (size_t i = 0; i < numberOfOriginsHistogram.size(); i++) {
            std::cout << std::setw(4) << i << " " << std::setw(4) << numberOfOriginsHistogram[i] << std::endl;
        }
        demand.sortByOrigin();
        Vertex origin = Vertex(0);
        std::vector<size_t> numberOfDestinationsHistogram;
        IndexedSet<true, Vertex> destinations;
        for (const AccumulatedVertexDemand::Entry& entry : demand.entries) {
            if ((entry.originVertex != origin) && (!destinations.empty())) {
                if (numberOfDestinationsHistogram.size() <= destinations.size()) {
                    numberOfDestinationsHistogram.resize(destinations.size() + 1, 0);
                }
                numberOfDestinationsHistogram[destinations.size()]++;
                destinations.clear();
                origin = entry.originVertex;
            }
            destinations.insert(entry.destinationVertex);
        }
        std::cout << "NumberOfDestinations Histogram." << std::endl;
        for (size_t i = 0; i < numberOfDestinationsHistogram.size(); i++) {
            std::cout << std::setw(4) << i << " " << std::setw(4) << numberOfDestinationsHistogram[i] << std::endl;
        }
    }

};

class BuildMultiModalAssignmentData : public ParameterizedCommand {

public:
    BuildMultiModalAssignmentData(BasicShell& shell) :
        ParameterizedCommand(shell, "buildMultiModalAssignmentData", "Combines the CSA network is a transfer graph. Zone vertices are preserve as sinks and sources.") {
        addParameter("CSA binary");
        addParameter("Road graph");
        addParameter("Output file");
        addParameter("Walking speed", "4.0");
        addParameter("Max degree", "16.0");
        addParameter("MaxUnificationDistanceInCM", "500");
        addParameter("MaxConnectingDistanceInCM", "10000");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string graphFileName = getParameter("Road graph");
        const std::string outputFileName = getParameter("Output file");
        const double walkingSpeed = getParameter<double>("Walking speed");
        const double maxCoreDegree = getParameter<double>("Max degree");
        const int maxUnificationDistanceInCM = getParameter<int>("MaxUnificationDistanceInCM");
        const int maxConnectingDistanceInCM = getParameter<int>("MaxConnectingDistanceInCM");

        std::cout << blue("Reading CSA input data") << std::endl;
        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;

        std::cout << blue("Converting data to intermediate format") << std::endl;
        Intermediate::Data interData = Intermediate::Data::FromCSA(csaData);
        interData.printInfo();

        std::cout << blue("Creating sources and sinks from zone vertices") << std::endl;
        interData.makeImpassableVertices();
        interData.printInfo();

        std::cout << blue("Writing source and sink zone networks") << std::endl;
        writeNetworks(interData, outputFileName + "SourceAndSinkZones/");

        std::cout << blue("Original number of vertices: ") << interData.transferGraph.numVertices() << std::endl;
        const size_t originalVertexCount = interData.transferGraph.numVertices();

        if (graphFileName != "-") {
            std::cout << blue("Loading walking graph") << std::endl;
            std::vector<Geometry::Point> coordinates = interData.transferGraph[Coordinates];
            for (Vertex v = Vertex(interData.numberOfStops()); v < originalVertexCount; v++) {
                interData.transferGraph.set(Coordinates, v, Geometry::Point(Construct::XY, 10000, 10000));
            }
            Intermediate::TransferGraph graph(graphFileName);
            Graph::printInfo(graph);
            graph.printAnalysis();
            Graph::applyBoundingBox(graph, interData.boundingBox());

            std::cout << blue("Applying max Walking speed of: ") << walkingSpeed << "km/h" << std::endl;
            Graph::computeTravelTimes(graph, walkingSpeed, true);
            Graph::printInfo(graph);
            graph.printAnalysis();

            std::cout << blue("Adding walking graph") << std::endl;
            interData.addTransferGraph(graph, maxUnificationDistanceInCM, maxConnectingDistanceInCM, walkingSpeed, true);
            for (Vertex v = Vertex(interData.numberOfStops()); v < originalVertexCount; v++) {
                interData.transferGraph.set(Coordinates, v, coordinates[v]);
            }
            interData.printInfo();
        }

        std::cout << blue("Writing multi-modal networks") << std::endl;
        writeNetworks(interData, outputFileName + "FullWalkinGraph/");

        buildCoreCH(interData, originalVertexCount, maxCoreDegree, outputFileName);
    }

private:
    inline void buildCoreCH(Intermediate::Data& interData, const size_t minCoreSize, const double maxCoreDegree, const std::string& outputFileName) const noexcept {
        std::cout << blue("Building Core CH") << std::endl;
        std::cout << "   min Core Size: " << String::prettyInt(minCoreSize) << std::endl;
        std::cout << "   max Core Degree: " << String::prettyInt(maxCoreDegree) << std::endl;
        Intermediate::TransferGraph resultGraph;
        resultGraph.addVertices(interData.transferGraph.numVertices());
        resultGraph[Coordinates] = interData.transferGraph[Coordinates];
        std::vector<bool> isContractable(minCoreSize, false);
        isContractable.resize(interData.transferGraph.numVertices(), true);
        using PROFILER = CH::FullProfiler;
        using WITNESS_SEARCH = CH::BidirectionalWitnessSearch<CHCoreGraph, PROFILER, 100>;
        using GREEDY_KEY_FUNCTION = CH::GreedyKey<WITNESS_SEARCH>;
        using KEY_FUNCTION = CH::PartialKey<WITNESS_SEARCH, GREEDY_KEY_FUNCTION>;
        using STOP_CRITERION = CH::CoreCriterion;
        KEY_FUNCTION keyFunction(isContractable, isContractable.size(), GREEDY_KEY_FUNCTION(1024, 256, 0));
        CH::Builder<PROFILER, WITNESS_SEARCH, KEY_FUNCTION, STOP_CRITERION, false, false> chBuilder(std::move(interData.transferGraph), interData.transferGraph[TravelTime], keyFunction, STOP_CRITERION(minCoreSize, maxCoreDegree));
        chBuilder.run();
        chBuilder.copyCoreToCH();
        std::cout << "   obtaining CH" << std::endl;
        CH::CH ch(std::move(chBuilder));
        for (const Vertex vertex : resultGraph.vertices()) {
            if (ch.isCoreVertex(vertex) || vertex < minCoreSize) {
                for (const Edge edge : ch.forward.edgesFrom(vertex)) {
                    resultGraph.addEdge(vertex, ch.forward.get(ToVertex, edge)).set(TravelTime, ch.forward.get(Weight, edge));
                }
            }
        }
        resultGraph.deleteVertices([&](const Vertex vertex){
            return (vertex >= minCoreSize) && (resultGraph.isIsolated(vertex));
        });
        Graph::move(std::move(resultGraph), interData.transferGraph);
        std::cout << blue("Writing core networks") << std::endl;
        writeNetworks(interData, outputFileName + "CoreCH/");
    }

    inline void writeNetworks(const Intermediate::Data& interData, const std::string& outputFileName) const noexcept {
        interData.serialize(outputFileName + "intermediate.binary");
        CSA::Data csaData = CSA::Data::FromIntermediate(interData);
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        csaData.serialize(outputFileName + "csa.binary");
        RAPTOR::Data raptorData = RAPTOR::Data::FromIntermediate(interData);
        raptorData.printInfo();
        raptorData.transferGraph.printAnalysis();
        raptorData.serialize(outputFileName + "raptor.binary");
    }

};

class BuildMultiModalDemandData : public ParameterizedCommand {

public:
    BuildMultiModalDemandData(BasicShell& shell) :
        ParameterizedCommand(shell, "buildMultiModalDemandData", "Adjusts the IDs of the origin and destination zones, in order to take impassable zones into account.") {
        addParameter("Demand file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");

        std::cout << blue("Adjusting demand for source and sink zones") << std::endl;
        AccumulatedVertexDemand::MakeImpassableZones(demandFileName, outputFileName);
        std::cout << std::endl;
    }

};

class RandomDemand : public ParameterizedCommand {

public:
    RandomDemand(BasicShell& shell) :
        ParameterizedCommand(shell, "randomDemand", "Adjusts the IDs of the origin and destination zones, in order to take impassable zones into account.") {
        addParameter("CSA binary");
        addParameter("Output file");
        addParameter("Size");
        addParameter("Stop based");
        addParameter("Vertex based");
        addParameter("Min departure time", "00:00:00");
        addParameter("Max departure time", "24:00:00");
        addParameter("Min departure window", "00:10:00");
        addParameter("Max departure window", "01:00:00");
        addParameter("Min group size", "1");
        addParameter("Max group size", "10");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string outputFileName = getParameter("Output file");
        const size_t size = getParameter<size_t>("Size");
        const bool stopBased = getParameter<bool>("Stop based");
        const bool vertexBased = getParameter<bool>("Vertex based");
        const int minDepartureTime = String::parseSeconds(getParameter("Min departure time"));
        const int maxDepartureTime = String::parseSeconds(getParameter("Max departure time"));
        const int minDepartureWindow = String::parseSeconds(getParameter("Min departure window"));
        const int maxDepartureWindow = String::parseSeconds(getParameter("Max departure window"));
        const int minGroupSize = getParameter<int>("Min group size");
        const int maxGroupSize = getParameter<int>("Max group size");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;

        AccumulatedVertexDemand demand = AccumulatedVertexDemand::Random(csaData, size, stopBased, vertexBased, minDepartureTime, maxDepartureTime, minDepartureWindow, maxDepartureWindow, minGroupSize, maxGroupSize);
        demand.toCSV(outputFileName);
    }

};

class TransferShortcutsForAssignment : public ParameterizedCommand {

public:
    TransferShortcutsForAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "transferShortcutsForAssignment", "Combines RAPTOR transfer shortcuts with an intermediate network including zones into a single CSA network.") {
        addParameter("Intermediate binary");
        addParameter("RAPTOR binary");
        addParameter("CSA binary");
        addParameter("verbose", "true");
    }

    virtual void execute() noexcept {
        const std::string intermediateFileName = getParameter("Intermediate binary");
        const std::string raptorFileName = getParameter("RAPTOR binary");
        const std::string csaFileName = getParameter("CSA binary");
        const bool verbose = getParameter<bool>("verbose");

        if (verbose) std::cout << blue("Reading intermediate input data") << std::endl;
        Intermediate::Data interData = Intermediate::Data::FromBinary(intermediateFileName);
        if (verbose) interData.printInfo();
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Reading RAPTOR input data") << std::endl;
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(raptorFileName);
        if (verbose) raptorData.printInfo();
        if (verbose) Graph::printInfo(raptorData.transferGraph);
        if (verbose) raptorData.transferGraph.printAnalysis();
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Combining edges") << std::endl;
        for (const Vertex from : raptorData.transferGraph.vertices()) {
            for (const Edge shortcut : raptorData.transferGraph.edgesFrom(from)) {
                const Vertex to = raptorData.transferGraph.get(ToVertex, shortcut);
                if (interData.transferGraph.hasEdge(from, to)) {
                    const Edge edge = interData.transferGraph.findEdge(from, to);
                    interData.transferGraph.set(TravelTime, edge, std::min<int>(raptorData.transferGraph.get(TravelTime, shortcut), interData.transferGraph.get(TravelTime, edge)));
                } else {
                    interData.transferGraph.addEdge(from, to).set(TravelTime, raptorData.transferGraph.get(TravelTime, shortcut));
                }
            }
        }
        if (verbose) interData.printInfo();
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Building CSA network") << std::endl;
        CSA::Data csaData = CSA::Data::FromIntermediate(interData);
        if (verbose) csaData.printInfo();
        if (verbose) csaData.transferGraph.printAnalysis();
        csaData.serialize(csaFileName);
    }

};

class SanitizeDemand : public ParameterizedCommand {

public:
    SanitizeDemand(BasicShell& shell) :
        ParameterizedCommand(shell, "sanitizeDemand", "Parses, sanitizes and writes given demand data.") {
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();

        AccumulatedVertexDemand demand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        demand.toZoneIDs(csaData);
        demand.lexicographicalSort();
        demand.sanitize();
        demand.toCSV(outputFileName);
    }
};

class GroupAssignment : public ParameterizedCommand {

public:
    GroupAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "groupAssignment", "Computes a public transit traffic assignment for zone based demand.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "0");
        addParameter("Thread offset", "1");
        addParameter("Aggregate file", "-");
        addParameter("Aggregate prefix", "-");
        addParameter("Use transfer buffer times", "false");
        addParameter("Demand output file", "-");
        addParameter("Demand output size", "-1");
    }

    virtual void execute() noexcept {
        const bool useTransferBufferTimes = getParameter<bool>("Use transfer buffer times");
        if (useTransferBufferTimes) {
            chooseProfiler<true>();
        } else {
            chooseProfiler<false>();
        }
    }

private:
    template<bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseProfiler() {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                chooseDecisionModel<Assignment::NoProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 1: {
                chooseDecisionModel<Assignment::TimeProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 2: {
                chooseDecisionModel<Assignment::DecisionProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        switch (settings.decisionModel) {
            case 0: {
                computeApportionment<Assignment::GroupAssignment<DecisionModels::Linear, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 1: {
                computeApportionment<Assignment::GroupAssignment<DecisionModels::Logit, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 2: {
                computeApportionment<Assignment::GroupAssignment<DecisionModels::Kirchhoff, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::GroupAssignment<DecisionModels::RelativeLogit, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::GroupAssignment<DecisionModels::Optimal, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
        }
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");
        const std::string aggregateFileName = getParameter("Aggregate file");
        const std::string aggregatePrefix = getParameter("Aggregate prefix");
        const std::string demandOutputFileName = getParameter("Demand output file");
        const size_t demandOutputSize = getParameter<size_t>("Demand output size");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }

        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        ma.getProfiler().printStatistics();
        if (demandOutputFileName != "-") {
            AccumulatedVertexDemand outputDemand = originalDemand;
            ma.filterDemand(outputDemand, demandOutputSize);
            outputDemand.toZoneIDs(csaData);
            outputDemand.sanitize();
            outputDemand.toCSV(demandOutputFileName);
        }
        // PassengerData result = ma.getPassengerData(originalDemand);
        // std::cout << result << std::endl;
        // result.serialize(FileSystem::ensureExtension(outputFileName, ".binary"));
        // IdVertexDemand passengers = IdVertexDemand::FromAccumulatedVertexDemand(originalDemand, settings.passengerMultiplier, 100000000);
        // result.writePassengerConnectionPairs(csaData, passengers, FileSystem::ensureExtension(outputFileName, ".csv"));
        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
        if (aggregateFileName != "-") {
            ma.writeConnectionStatistics(aggregateFileName, aggregatePrefix);
        }
    }
};

class SequentialAssignment : public ParameterizedCommand {

public:
    SequentialAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "sequentialAssignment", "Computes a public transit traffic assignment with sequential logit model for zone-based demand.") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "1");
        addParameter("Thread offset", "1");
        addParameter("Aggregate file", "-");
        addParameter("Aggregate prefix", "-");
    }

    virtual void execute() noexcept {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                computeAssignment<Assignment::NoProfiler>(settings);
                break;
            }
            case 1: {
                computeAssignment<Assignment::TimeProfiler>(settings);
                break;
            }
            case 2: {
                computeAssignment<Assignment::DecisionProfiler>(settings);
                break;
            }
        }
    }

private:
    template<typename PROFILER>
    inline void computeAssignment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");
        const std::string aggregateFileName = getParameter("Aggregate file");
        const std::string aggregatePrefix = getParameter("Aggregate prefix");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        Assignment::SequentialAssignment<PROFILER> assignment(csaData, reverseGraph, settings);
        Timer timer;
        const int numCores(numberOfCores());
        std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
        assignment.run(demand, numThreads, pinMultiplier);

        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        assignment.getProfiler().printStatistics();
        assignment.printStatistics(originalDemand, outputFileName);
        assignment.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        assignment.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        assignment.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        assignment.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
        if (aggregateFileName != "-") {
            assignment.writeConnectionStatistics(aggregateFileName, aggregatePrefix);
        }
    }
};

class CapacityAssignment : public ParameterizedCommand {

public:
    CapacityAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "capacityAssignment", "Computes a public transit traffic assignment with vehicle capacities for zone-based demand.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "1");
        addParameter("Thread offset", "1");
        addParameter("Use transfer buffer times", "false");
    }

    virtual void execute() noexcept {
        /*const bool useTransferBufferTimes = getParameter<bool>("Use transfer buffer times");
        if (useTransferBufferTimes) {
            chooseProfiler<true>();
        } else {
            chooseProfiler<false>();
        }*/
        chooseProfiler<false>();
    }

private:
    template<bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseProfiler() {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        /*switch (settings.profilerType) {
            case 0: {*/
                chooseDecisionModel<Assignment::NoProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
            /*    break;
            }
            case 1: {
                chooseDecisionModel<Assignment::TimeProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 2: {
                chooseDecisionModel<Assignment::DecisionProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
        }*/
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        /*switch (settings.decisionModel) {
            case 0: {*/
                computeApportionment<Assignment::Capacities::Assignment<DecisionModels::Linear, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
            /*    break;
            }
            case 1: {
                computeApportionment<Assignment::Capacities::Assignment<DecisionModels::Logit, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 2: {
                computeApportionment<Assignment::Capacities::Assignment<DecisionModels::Kirchhoff, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::Capacities::Assignment<DecisionModels::RelativeLogit, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::Capacities::Assignment<DecisionModels::Optimal, PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
        }*/
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        std::vector<double> connectionCapacity(csaData.numberOfConnections(), 800); //TODO: Dummy capacity

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, connectionCapacity, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        ma.getProfiler().printStatistics();

        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
    }
};

class LockstepCapacityAssignment : public ParameterizedCommand {

public:
    LockstepCapacityAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "lockstepCapacityAssignment", "Computes a public transit traffic assignment with vehicle capacities for zone-based demand.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Capacity file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "1");
        addParameter("Thread offset", "1");
        addParameter("Use transfer buffer times", "false");
    }

    virtual void execute() noexcept {
        /*const bool useTransferBufferTimes = getParameter<bool>("Use transfer buffer times");
        if (useTransferBufferTimes) {
            chooseProfiler<true>();
        } else {
            chooseProfiler<false>();
        }*/
        chooseProfiler<false>();
    }

private:
    template<bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseProfiler() {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                chooseConnectionProfiler<Assignment::NoProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 1: {
                chooseConnectionProfiler<Assignment::TimeProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 2: {
                chooseConnectionProfiler<Assignment::DecisionProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseConnectionProfiler(const Assignment::Settings& settings) {
        if (settings.useConnectionProfiler) {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, true>(settings);
        } else {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, false>(settings);
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES, bool USE_CONNECTION_PROFILER>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        /*switch (settings.decisionModel) {
            case 0: {
                computeApportionment<Assignment::Capacities::LockstepAssignment<DecisionModels::Linear, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 1: {*/
                computeApportionment<Assignment::Capacities::LockstepAssignment<DecisionModels::Logit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
            /*    break;
            }
            case 2: {
                computeApportionment<Assignment::Capacities::LockstepAssignment<DecisionModels::Kirchhoff, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::Capacities::LockstepAssignment<DecisionModels::RelativeLogit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::Capacities::LockstepAssignment<DecisionModels::Optimal, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
        }*/
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string capacityFile = getParameter("Capacity file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        std::vector<double> connectionCapacity = parseCapacities(capacityFile, csaData.numberOfConnections());

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, connectionCapacity, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        //ma.getProfiler().printStatistics();

        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
    }

    inline static std::vector<double> parseCapacities(const std::string filename, const size_t numConnections) noexcept {
        IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(filename);
        in.readHeader("connection_id", "capacity");
        std::vector<double> result(numConnections, 0);
        size_t id;
        double value;
        while (in.readRow(id, value)) {
            Ensure(id < result.size(), "Capacity file has too many entries!");
            result[id] = value;
        }
        return result;
    }
};

class FixedCapacityAssignment : public ParameterizedCommand {

public:
    FixedCapacityAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "fixedCapacityAssignment", "Computes a public transit traffic assignment with vehicle capacities for zone-based demand.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Capacity file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "1");
        addParameter("Thread offset", "1");
        addParameter("Use transfer buffer times", "false");
    }

    virtual void execute() noexcept {
        /*const bool useTransferBufferTimes = getParameter<bool>("Use transfer buffer times");
        if (useTransferBufferTimes) {
            chooseProfiler<true>();
        } else {
            chooseProfiler<false>();
        }*/
        chooseProfiler<false>();
    }

private:
    template<bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseProfiler() {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                chooseConnectionProfiler<Assignment::NoProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 1: {
                chooseConnectionProfiler<Assignment::TimeProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 2: {
                chooseConnectionProfiler<Assignment::DecisionProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseConnectionProfiler(const Assignment::Settings& settings) {
        if (settings.useConnectionProfiler) {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, true>(settings);
        } else {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, false>(settings);
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES, bool USE_CONNECTION_PROFILER>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        /*switch (settings.decisionModel) {
            case 0: {
                computeApportionment<Assignment::Capacities::FixedAssignment<DecisionModels::Linear, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 1: {*/
                computeApportionment<Assignment::Capacities::FixedAssignment<DecisionModels::Logit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
            /*    break;
            }
            case 2: {
                computeApportionment<Assignment::Capacities::FixedAssignment<DecisionModels::Kirchhoff, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::Capacities::FixedAssignment<DecisionModels::RelativeLogit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::Capacities::FixedAssignment<DecisionModels::Optimal, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
        }*/
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string capacityFile = getParameter("Capacity file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        std::vector<double> connectionCapacity = parseCapacities(capacityFile, csaData.numberOfConnections());

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, connectionCapacity, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        //ma.getProfiler().printStatistics();

        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
    }

    inline static std::vector<double> parseCapacities(const std::string filename, const size_t numConnections) noexcept {
        IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(filename);
        in.readHeader("connection_id", "capacity");
        std::vector<double> result(numConnections, 0);
        size_t id;
        double value;
        while (in.readRow(id, value)) {
            Ensure(id < result.size(), "Capacity file has too many entries!");
            result[id] = value;
        }
        return result;
    }
};

class ConnectionBasedLockstepCapacityAssignment : public ParameterizedCommand {

public:
    ConnectionBasedLockstepCapacityAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "connectionBasedLockstepCapacityAssignment", "Computes a public transit traffic assignment with vehicle capacities for zone-based demand.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Capacity file");
        addParameter("Output file");
        addParameter("Demand multiplier", "1");
        addParameter("Num threads", "1");
        addParameter("Thread offset", "1");
        addParameter("Use transfer buffer times", "false");
    }

    virtual void execute() noexcept {
        /*const bool useTransferBufferTimes = getParameter<bool>("Use transfer buffer times");
        if (useTransferBufferTimes) {
            chooseProfiler<true>();
        } else {
            chooseProfiler<false>();
        }*/
        chooseProfiler<false>();
    }

private:
    template<bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseProfiler() {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                chooseConnectionProfiler<Assignment::NoProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 1: {
                chooseConnectionProfiler<Assignment::TimeProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
            case 2: {
                chooseConnectionProfiler<Assignment::DecisionProfiler, USE_TRANSFER_BUFFER_TIMES>(settings);
                break;
            }
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES>
    inline void chooseConnectionProfiler(const Assignment::Settings& settings) {
        if (settings.useConnectionProfiler) {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, true>(settings);
        } else {
            chooseDecisionModel<PROFILER, USE_TRANSFER_BUFFER_TIMES, false>(settings);
        }
    }

    template<typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES, bool USE_CONNECTION_PROFILER>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        /*switch (settings.decisionModel) {
            case 0: {
                computeApportionment<Assignment::Capacities::ConnectionBasedLockstepAssignment<DecisionModels::Linear, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 1: {*/
                computeApportionment<Assignment::Capacities::ConnectionBasedLockstepAssignment<DecisionModels::Logit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
            /*    break;
            }
            case 2: {
                computeApportionment<Assignment::Capacities::ConnectionBasedLockstepAssignment<DecisionModels::Kirchhoff, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::Capacities::ConnectionBasedLockstepAssignment<DecisionModels::RelativeLogit, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::Capacities::ConnectionBasedLockstepAssignment<DecisionModels::Optimal, PROFILER, USE_CONNECTION_PROFILER, USE_TRANSFER_BUFFER_TIMES>>(settings);
                break;
            }
        }*/
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string capacityFile = getParameter("Capacity file");
        const std::string outputFileName = getParameter("Output file");
        const size_t demandMultiplier = getParameter<size_t>("Demand multiplier");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph, demandMultiplier);
        AccumulatedVertexDemand demand = originalDemand;

        std::vector<double> connectionCapacity = parseCapacities(capacityFile, csaData.numberOfConnections());

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, connectionCapacity, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        //ma.getProfiler().printStatistics();

        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"), demand);
    }

    inline static std::vector<double> parseCapacities(const std::string filename, const size_t numConnections) noexcept {
        IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(filename);
        in.readHeader("connection_id", "capacity");
        std::vector<double> result(numConnections, 0);
        size_t id;
        double value;
        while (in.readRow(id, value)) {
            Ensure(id < result.size(), "Capacity file has too many entries!");
            result[id] = value;
        }
        return result;
    }
};

class ParseCapacities : public ParameterizedCommand {
public:
    ParseCapacities(BasicShell& shell) :
        ParameterizedCommand(shell, "parseCapacities", "Parses capacity information from CSV network data.") {
        addParameter("Input file");
        addParameter("CSA binary");
        addParameter("Output file");
        addParameter("Scaling factor", "1.0");
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("Input file");
        const double scalingFactor = getParameter<double>("Scaling factor");
        std::vector<double> tripCapacities;
        IO::readFile(inputFile, "Trips", [&](){
            size_t count = 0;
            IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(inputFile);
            in.readHeader(IO::IGNORE_EXTRA_COLUMN | IO::IGNORE_MISSING_COLUMN, "trip_id", "total_cap");
            TripId tripID;
            double capacity;
            while (in.readRow(tripID, capacity)) {
                if (tripID >= tripCapacities.size()) tripCapacities.resize(tripID + 1, 0);
                tripCapacities[tripID] = capacity * scalingFactor;
                count++;
            }
            return count;
        }, true);

        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA binary"));
        Ensure(csaData.numberOfTrips() == tripCapacities.size(), "CSA data has " << csaData.numberOfTrips() << " trips, but input file has " << tripCapacities.size() << "!");
        std::vector<double> connectionCapacities;
        for (const CSA::Connection& connection : csaData.connections) {
            const TripId trip = connection.tripId;
            connectionCapacities.emplace_back(tripCapacities[trip]);
        }

        IO::OFStream outputFile(getParameter("Output file"));
        outputFile << "connection_id,capacity\n";
        for (size_t i = 0; i < connectionCapacities.size(); i++) {
            outputFile << i << "," << connectionCapacities[i] << "\n";
        }
    }
};

class BuildCapacityAssignmentTestNetwork : public ParameterizedCommand {
public:
    BuildCapacityAssignmentTestNetwork(BasicShell& shell) :
        ParameterizedCommand(shell, "buildCapacityAssignmentTestNetwork", "Builds a test network for capacity assignment.") {
        addParameter("Directory");
    }

    virtual void execute() noexcept {
        const std::string directory = getParameter("Directory");

        EdgeList<WithCoordinates, WithTravelTime> transferGraph;
        transferGraph.addVertices(7);
        transferGraph.addEdge(Vertex(3), Vertex(0)).set(TravelTime, 300);
        transferGraph.addEdge(Vertex(4), Vertex(0)).set(TravelTime, 150);
        transferGraph.addEdge(Vertex(1), Vertex(5)).set(TravelTime, 300);
        transferGraph.addEdge(Vertex(2), Vertex(6)).set(TravelTime, 300);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates,  Vertex(0)),  0);
        stops.emplace_back("A", transferGraph.get(Coordinates,  Vertex(1)),  0);
        stops.emplace_back("T", transferGraph.get(Coordinates,  Vertex(2)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> T", "R1", 1);
        trips.emplace_back("S -> A", "R2", 1);
        trips.emplace_back("A -> T", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 600, 900, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 900, 1200, TripId(0));
        connections.emplace_back(StopId(0), StopId(1), 600, 1020, TripId(1));
        connections.emplace_back(StopId(1), StopId(2), 1080, 1500, TripId(2));

        CSA::Data data = CSA::Data::FromInput(stops, connections, trips, transferGraph);
        data.serialize(directory + "network");

        AccumulatedVertexDemand demand;
        demand.entries.emplace_back(300, 300, Vertex(0), Vertex(2), 100);
        demand.entries.emplace_back(300, 300, Vertex(0), Vertex(3), 100);
        demand.entries.emplace_back(300, 300, Vertex(1), Vertex(2), 100);
        demand.entries.emplace_back(300, 300, Vertex(1), Vertex(3), 100);
        for (size_t i = 0; i < demand.entries.size(); i++) {
            demand.entries[i].demandIndex = i;
        }
        demand.numberOfPassengers = 400;
        demand.toCSV(directory + "demand");
    }
};

class BuildStations : public ParameterizedCommand {

private:
    struct EdgeData {
        EdgeData(const StopId from, const StopId to, const int length) : from(from), to(to), length(length) {}
        inline bool operator<(const EdgeData& other) const noexcept {
            return length < other.length;
        }
        StopId from;
        StopId to;
        int length;
    };

public:
    BuildStations(BasicShell& shell) :
        ParameterizedCommand(shell, "buildStations", "Partitions the stops of a CSA network into stations. If <Max Station Size> is '-', than transitive stations are constructed.") {
        addParameter("CSA binary");
        addParameter("Station file");
        addParameter("Max Station Size", "-");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string stationFileName = getParameter("Station file");

        std::cout << blue("Building stations") << std::endl;
        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        size_t numberOfStations = 0;
        std::vector<int> stationOfStop(csaData.numberOfStops());
        if (getParameter("Max Station Size") == "-") {
            for (const StopId stop : csaData.stops()) {
                stationOfStop[stop] = stop;
                for (const Edge edge : csaData.transferGraph.edgesFrom(stop)) {
                    stationOfStop[stop] = std::min<int>(stationOfStop[stop], csaData.transferGraph.get(ToVertex, edge));
                }
            }
        } else {
            const int maxStationSize = getParameter<int>("Max Station Size");
            std::vector<EdgeData> edges;
            Dijkstra<CSA::TransferGraph> dijkstra(csaData.transferGraph);
            for (const StopId stop : csaData.stops()) {
                dijkstra.run(stop, noVertex, [&](const Vertex u) {
                    if (csaData.isStop(u)) {
                        edges.emplace_back(stop, StopId(u), dijkstra.getDistance(u));
                    }
                }, [&](){
                    return dijkstra.getDistance(dijkstra.getQFront()) > maxStationSize;
                });
            }
            sort(edges);
            UnionFind stations(csaData.numberOfStops());
            std::vector<int> stationSize(csaData.numberOfStops(), 0);
            for (const EdgeData& edge : edges) {
                if (stations(edge.from) == stations(edge.to)) continue;
                const int combinedSize = edge.length + stationSize[stations(edge.from)] + stationSize[stations(edge.to)];
                if (combinedSize > maxStationSize) continue;
                stations.unite(edge.from, edge.to);
                stationSize[stations(edge.from)] = combinedSize;
            }
            for (const StopId stop : csaData.stops()) {
                stationOfStop[stop] = stations(stop);
            }
        }
        for (const StopId stop : csaData.stops()) {
            if (StopId(stationOfStop[stop]) == stop) numberOfStations++;
        }
        IO::serialize(stationFileName, stationOfStop);
        shell << "   Number of stations: " << numberOfStations << newLine;
    }

};

class UltraGroupAssignment : public ParameterizedCommand {

public:
    UltraGroupAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "ultraGroupAssignment", "Computes a public transit traffic assignment for a multi-modal network using ULTRA based PAT computation.", "num threads:", "    positive number  - parallel execution with <num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("CH file");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Num threads", "0");
        addParameter("Thread offset", "1");
        addParameter("Aggregate file", "-");
        addParameter("Aggregate head", "-");
        addParameter("verbose", "false");
    }

    virtual void execute() noexcept {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        chooseProfiler(settings);
    }

private:
    inline void chooseProfiler(const Assignment::Settings& settings) {
        switch (settings.profilerType) {
            case 0: {
                chooseDecisionModel<Assignment::NoProfiler>(settings);
                break;
            }
            case 1: {
                chooseDecisionModel<Assignment::TimeProfiler>(settings);
                break;
            }
            case 2: {
                chooseDecisionModel<Assignment::DecisionProfiler>(settings);
                break;
            }
        }
    }

    template<typename PROFILER>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        switch (settings.decisionModel) {
            case 0: {
                computeApportionment<Assignment::ULTRA::Assignment<DecisionModels::Linear, PROFILER>>(settings);
                break;
            }
            case 1: {
                computeApportionment<Assignment::ULTRA::Assignment<DecisionModels::Logit, PROFILER>>(settings);
                break;
            }
            case 2: {
                computeApportionment<Assignment::ULTRA::Assignment<DecisionModels::Kirchhoff, PROFILER>>(settings);
                break;
            }
            case 3: {
                computeApportionment<Assignment::ULTRA::Assignment<DecisionModels::RelativeLogit, PROFILER>>(settings);
                break;
            }
            case 4: {
                computeApportionment<Assignment::ULTRA::Assignment<DecisionModels::Optimal, PROFILER>>(settings);
                break;
            }
        }
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeApportionment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string chFileName = getParameter("CH file");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");
        const std::string aggregateFileName = getParameter("Aggregate file");
        const std::string aggregateHead = getParameter("Aggregate head");
        const bool verbose = getParameter<bool>("verbose");

        if (verbose) std::cout << blue("Loading CSA network") << std::endl;
        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        if (verbose) csaData.printInfo();
        if (verbose) csaData.transferGraph.printAnalysis();
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Building reverse graph") << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        if (verbose) Graph::printInfo(reverseGraph);
        if (verbose) reverseGraph.printAnalysis();
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Loading CH") << std::endl;
        CH::CH ch(chFileName);
        if (verbose) CH::analyze(ch);
        if (verbose) std::cout << std::endl;

        if (verbose) std::cout << blue("Loading Demand") << std::endl;
        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph);
        AccumulatedVertexDemand demand = originalDemand;
        if (verbose) std::cout << std::endl;

        if (settings.demandIntervalSplitTime >= 0) {
            demand.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
        }
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, ch, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            std::cout << "Running in sequential mode!" << std::endl;
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        ma.getProfiler().printStatistics();
        // PassengerData result = ma.getPassengerData(originalDemand);
        // std::cout << result << std::endl;
        // result.serialize(FileSystem::ensureExtension(outputFileName, ".binary"));
        // result.writePassengerConnectionPairs(csaData, passengers, FileSystem::ensureExtension(outputFileName, ".csv"));
        ma.printStatistics(originalDemand, outputFileName);
        ma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"));
        ma.writeAssignment(FileSystem::ensureExtension(outputFileName, "_assignment.csv"));
        ma.writeGroups(FileSystem::ensureExtension(outputFileName, "_groups.csv"));
        // ma.writeAssignedJourneys(FileSystem::ensureExtension(outputFileName, "_journeys.csv"));
        if (aggregateFileName != "-") {
            ma.writeAggregateStatistics(aggregateFileName, aggregateHead);
        }
        /*if (csaData.numberOfConnections() == 24) {
            std::cout << std::endl;
            std::cout << ma.getPassengerCountForConnection(0) << std::endl;
            std::cout << ma.getPassengerCountForConnection(7) << std::endl;
            std::cout << ma.getPassengerCountForConnection(13) << std::endl;
            std::cout << std::endl;
            std::cout << ma.getPassengerCountForConnection(2) << std::endl;
            std::cout << ma.getPassengerCountForConnection(9) << std::endl;
            std::cout << ma.getPassengerCountForConnection(18) << std::endl;
            std::cout << std::endl;
            std::cout << ma.getPassengerCountForConnection(1) << std::endl;
            std::cout << ma.getPassengerCountForConnection(10) << std::endl;
            std::cout << ma.getPassengerCountForConnection(16) << std::endl;
        }*/
    }

};

class AssignmentTest : public ParameterizedCommand {

public:
    AssignmentTest(BasicShell& shell) :
        ParameterizedCommand(shell, "assignmentTest", "Loads binary CSA data and computes assignments with increasing passenger multiplier.") {
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Num iterations");
        addParameter("Settings file");
        addParameter("Assignment type", {"group", "ultra"});
        addParameter("CH file", "-");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const int numIterations = getParameter<int>("Num iterations");
        const std::string settingsFileName = getParameter("Settings file");
        const std::string assignmentType = getParameter("Assignment type");
        const std::string chFileName = getParameter("CH file");

        pinThreadToCoreId(0);
        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand demand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph);
        std::cout.precision(10);
        const std::vector<int> passengerMultipliers({1, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500});
        if (assignmentType == "group") {
            Assignment::GroupAssignment<DecisionModels::Linear, Assignment::TimeProfiler> ma(csaData, reverseGraph, settings);
            Assignment::TimeProfiler& d = ma.getProfiler();
            for (const int passengerMultiplier : passengerMultipliers) {
                AccumulatedVertexDemand passengers = demand;
                passengers.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
                settings.passengerMultiplier = passengerMultiplier;
                for (int i = 0; i < numIterations; i++) {
                    d.initialize();
                    ma.run(passengers);
                    d.printStatistics();
                    std::cout << "Passenger multiplier: " << passengerMultiplier << std::endl;
                }
            }
        } else {
            CH::CH ch(chFileName);
            Assignment::ULTRA::Assignment<DecisionModels::Kirchhoff, Assignment::TimeProfiler> ma(csaData, reverseGraph, ch, settings);
            Assignment::TimeProfiler& d = ma.getProfiler();
            for (const int passengerMultiplier : passengerMultipliers) {
                AccumulatedVertexDemand passengers = demand;
                if (settings.demandIntervalSplitTime >= 0) {
                    passengers.discretize(settings.demandIntervalSplitTime, settings.keepDemandIntervals, settings.includeIntervalBorder);
                }
                settings.passengerMultiplier = passengerMultiplier;
                for (int i = 0; i < numIterations; i++) {
                    d.initialize();
                    ma.run(passengers);
                    d.printStatistics();
                    std::cout << "Passenger multiplier: " << passengerMultiplier << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }

};

class AssignmentStatistics : public ParameterizedCommand {

public:
    AssignmentStatistics(BasicShell& shell) :
        ParameterizedCommand(shell, "assignmentStatistics", "Loads binary CSA data, a demand, and an assignment and computes statistics.") {
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Assignment file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string assignmentFileName = getParameter("Assignment file");
        const std::string outputFileName = getParameter("Output file");

        pinThreadToCoreId(0);
        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand demand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph);

        std::cout << blue("Reading assignment data") << std::endl;
        std::vector<ConnectionData> connectionData(csaData.numberOfConnections());
        std::vector<PassengerData> passengerData(demand.numberOfPassengers);
        IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(assignmentFileName);
        in.readHeader("connection_id", "passenger_id");
        ConnectionId connection_id;
        size_t passenger_id;
        while (in.readRow(connection_id, passenger_id)) {
            AssertMsg(connection_id < connectionData.size(), "connection_id < connectionData.size()" << connection_id);
            AssertMsg(passenger_id < passengerData.size(), "passenger_id < passengerData.size()" << passenger_id);
            passengerData[passenger_id].connections.emplace_back(connection_id);
            connectionData[connection_id].passengerCount++;
        }

        std::cout << blue("Generating passengers from demand") << std::endl;
        size_t passenger = 0;
        for (const AccumulatedVertexDemand::Entry& entry : demand.entries) {
            for (size_t i = 0; i < entry.numberOfPassengers; i++) {
                passengerData[passenger].origin = entry.originVertex;
                passengerData[passenger].destination = entry.destinationVertex;
                passengerData[passenger].departureTime = entry.earliestDepartureTime;
                passenger++;
            }
        }

        std::cout << blue("Computing statistics") << std::endl;
        for (const size_t i : indices(passengerData)) {
            PassengerData& data = passengerData[i];
            if (data.origin == data.destination) {
                continue;
            } else if (data.connections.empty()) {
                if (csaData.transferGraph.hasEdge(data.origin, data.destination)) {
                    int time = csaData.transferGraph.get(TravelTime, csaData.transferGraph.findEdge(data.origin, data.destination));
                    data.traveltimeDemand = time;
                    data.traveltimeJourney = time;
                    data.walkingTime = time;
                    data.waitingTime = 0;
                } else {
                    data.origin = noVertex;
                    continue;
                }
            } else {
                data.trips++;
                const CSA::Connection& c = csaData.connections[data.connections.front()];
                if (data.origin == c.departureStopId) {
                    int time = c.departureTime - data.departureTime;
                    data.traveltimeDemand = time;
                    data.traveltimeJourney = 0;
                    data.walkingTime = 0;
                    data.waitingTime = time;
                    if (time < 0) warning(4);
                    if (data.waitingTime < 0) warning(1);
                    if (data.waitingTime < 0) warning("               " + std::to_string(time));
                } else {
                    if (csaData.transferGraph.hasEdge(data.origin, c.departureStopId)) {
                        int time = csaData.transferGraph.get(TravelTime, csaData.transferGraph.findEdge(data.origin, c.departureStopId));
                        int time2 = c.departureTime - data.departureTime - time;
                        data.traveltimeDemand = time + time2;
                        data.traveltimeJourney = time;
                        data.walkingTime = time;
                        data.waitingTime = time2;
                        if (time2 < 0) warning(5);
                        if (data.waitingTime < 0) warning(1);
                        if (data.waitingTime < 0) warning("               " + std::to_string(time));
                        if (data.waitingTime < 0) warning("               " + std::to_string(time2));
                    }
                }
                const CSA::Connection& c2 = csaData.connections[data.connections.back()];
                if (data.destination != c2.arrivalStopId) {
                    if (csaData.transferGraph.hasEdge(c2.arrivalStopId, data.destination)) {
                        int time = csaData.transferGraph.get(TravelTime, csaData.transferGraph.findEdge(c2.arrivalStopId, data.destination));
                        data.traveltimeDemand += time;
                        data.traveltimeJourney += time;
                        data.walkingTime += time;
                    }
                }
                data.traveltimeDemand += c.travelTime();
                data.traveltimeJourney += c.travelTime();
                for (size_t j = 1; j < data.connections.size(); j++) {
                    const CSA::Connection& c3 = csaData.connections[data.connections[j - 1]];
                    const CSA::Connection& c4 = csaData.connections[data.connections[j]];
                    if (c3.tripId == c4.tripId) {
                        int time = c4.arrivalTime - c3.arrivalTime;
                        data.traveltimeDemand += time;
                        data.traveltimeJourney += time;
                    } else {
                        data.trips++;
                        data.traveltimeDemand += c4.travelTime();
                        data.traveltimeJourney += c4.travelTime();
                        if (c3.arrivalStopId == c4.departureStopId) {
                            int time = c4.departureTime - c3.arrivalTime;
                            data.traveltimeDemand += time;
                            data.traveltimeJourney += time;
                            data.waitingTime = time;
                            if (data.waitingTime < 0) warning(2);
                            if (data.waitingTime < 0) warning("               " + std::to_string(time));
                        } else {
                            if (csaData.transferGraph.hasEdge(c3.arrivalStopId, c4.departureStopId)) {
                                int time = csaData.transferGraph.get(TravelTime, csaData.transferGraph.findEdge(c3.arrivalStopId, c4.departureStopId));
                                int time2 = c4.departureTime - c3.arrivalTime - time;
                                data.traveltimeDemand += time + time2;
                                data.traveltimeJourney += time + time2;
                                data.walkingTime += time;
                                data.waitingTime += time2;
                                if (data.waitingTime < 0) warning(3);
                                if (data.waitingTime < 0) warning("               " + std::to_string(time));
                                if (data.waitingTime < 0) warning("               " + std::to_string(time2));
                            }
                        }
                    }
                }
            }
        }
        write(outputFileName, connectionData);
        write(outputFileName, passengerData);
        std::cout << std::endl;
    }

private:
    struct PassengerData {
        Vertex origin{noVertex};
        Vertex destination{noVertex};
        int departureTime{-1};
        int traveltimeDemand{0};
        int traveltimeJourney{0};
        int walkingTime{0};
        int waitingTime{0};
        int trips{0};
        std::vector<int> connections;
        inline int vehicleTime() const noexcept {
            return traveltimeDemand - walkingTime - waitingTime;
        }
    };
    struct ConnectionData {
        int passengerCount{0};
    };

    inline void write(const std::string& filename, const std::vector<PassengerData>& data) const noexcept {
        std::string file = filename + "passenger.csv";
        std::ofstream os(file);
        AssertMsg(os, "cannot open file: " << file);
        AssertMsg(os.is_open(), "cannot open file: " << file);
        os << "id,traveltimeDemand,traveltimeJourney,vehicleTime,walkingTime,waitingTime,numConnections,trips\n";
        for (const size_t i : indices(data)) {
            if (data[i].origin < 1) continue;
            os << i << "," << data[i].traveltimeDemand << "," << data[i].traveltimeJourney << ", " << data[i].vehicleTime() << "," << data[i].walkingTime << "," << data[i].waitingTime << "," << data[i].connections.size() << "," << data[i].trips << "\n";
        }
    }

    inline void write(const std::string& filename, const std::vector<ConnectionData>& data) const noexcept {
        std::string file = filename + "_connections.csv";
        std::ofstream os(file);
        AssertMsg(os, "cannot open file: " << file);
        AssertMsg(os.is_open(), "cannot open file: " << file);
        os << "connection,numPassengers\n";
        for (const size_t i : indices(data)) {
            os << i << "," << data[i].passengerCount << "\n";
        }
    }

};

class AssignmentGraph : public ParameterizedCommand {

public:
    AssignmentGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "assignmentGraph", "Loads binary CSA data and computes an assignment. Afterwards, the assignment is visualized as graph.", "num threads:", "    positive number  - parallel execution with <num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Destination vertex");
        addParameter("Output file");
        addParameter("Min departure time", "0");
        addParameter("Max departure time", "86400");
        addParameter("Num threads", "0");
        addParameter("Thread offset", "1");
    }

    virtual void execute() noexcept {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        chooseProfiler(settings);
    }

private:
    inline void chooseProfiler(const Assignment::Settings& settings) {
        switch (settings.profilerType) {
            case 0: {
                chooseDecisionModel<Assignment::NoProfiler>(settings);
                break;
            }
            case 1: {
                chooseDecisionModel<Assignment::TimeProfiler>(settings);
                break;
            }
            case 2: {
                chooseDecisionModel<Assignment::DecisionProfiler>(settings);
                break;
            }
        }
    }

    template<typename PROFILER>
    inline void chooseDecisionModel(const Assignment::Settings& settings) {
        switch (settings.decisionModel) {
            case 0: {
                computeAssignmentGraph<Assignment::GroupAssignment<DecisionModels::Linear, PROFILER>>(settings);
                break;
            }
            case 1: {
                computeAssignmentGraph<Assignment::GroupAssignment<DecisionModels::Logit, PROFILER>>(settings);
                break;
            }
            case 2: {
                computeAssignmentGraph<Assignment::GroupAssignment<DecisionModels::Kirchhoff, PROFILER>>(settings);
                break;
            }
            case 3: {
                computeAssignmentGraph<Assignment::GroupAssignment<DecisionModels::RelativeLogit, PROFILER>>(settings);
                break;
            }
            case 4: {
                computeAssignmentGraph<Assignment::GroupAssignment<DecisionModels::Optimal, PROFILER>>(settings);
                break;
            }
        }
    }

    template<typename APPORTIONMENT_TYPE>
    inline void computeAssignmentGraph(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const Vertex destination = Vertex(getParameter<int>("Destination vertex"));
        const std::string outputFileName = getParameter("Output file");
        const int minDepartureTime = getParameter<int>("Min departure time");
        const int maxDepartureTime = getParameter<int>("Max departure time");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();
        AccumulatedVertexDemand demand = AccumulatedVertexDemand::ForDestination(csaData, destination, minDepartureTime, maxDepartureTime);
        APPORTIONMENT_TYPE ma(csaData, reverseGraph, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            ma.run(demand, numThreads, pinMultiplier);
        } else {
            ma.run(demand);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(ma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(ma.getRemovedCycles()) << std::endl;
        ma.getProfiler().printStatistics();
        SimpleDynamicGraph assignmentGraph;
        assignmentGraph.addVertices(csaData.numberOfStops());
        for (const std::vector<ConnectionId>& connectionIds : ma.getAssignmentData().connectionsPerGroup) {
            if (connectionIds.empty()) continue;
            const CSA::Connection& firstConnection = csaData.connections[connectionIds[0]];
            assignmentGraph.findOrAddEdge(firstConnection.departureStopId, firstConnection.arrivalStopId);
            for (size_t i = 1; i < connectionIds.size(); i++) {
                const CSA::Connection& connectionA = csaData.connections[connectionIds[i - 1]];
                const CSA::Connection& connectionB = csaData.connections[connectionIds[i]];
                if (connectionA.arrivalStopId != connectionB.departureStopId) assignmentGraph.findOrAddEdge(connectionA.arrivalStopId, connectionB.departureStopId);
                assignmentGraph.findOrAddEdge(connectionB.departureStopId, connectionB.arrivalStopId);
            }
        }
        Graph::printInfo(assignmentGraph);
        assignmentGraph.printAnalysis();
        Graph::toGML(outputFileName + ".complete", assignmentGraph);
        IndexedSet<false, Vertex> vertices = IndexedSet<false, Vertex>(Construct::Complete, assignmentGraph.numVertices());
        while (!vertices.empty()) {
            Vertex vertex = vertices.back();
            vertices.remove(vertex);
            if (assignmentGraph.inDegree(vertex) > 0 && assignmentGraph.outDegree(vertex) > 0) continue;
            for (const Vertex neighbor : assignmentGraph.neighbors(vertex)) {
                vertices.insert(neighbor);
            }
            assignmentGraph.isolateVertex(vertex);
        }
        Graph::printInfo(assignmentGraph);
        assignmentGraph.printAnalysis();
        Graph::toGML(outputFileName + ".reduced", assignmentGraph);

        MapVisualization<PDF> vis(outputFileName, csaData.boundingBox());
        for (const Border& border : Borders::GermanStates) {
            vis.drawLine(border, Color::Grey, 5, false);
        }
        vis.drawLine(Borders::Germany, Color::Black, 8, false);
        for (const Vertex from : assignmentGraph.vertices()) {
            for (const Edge edge : assignmentGraph.edgesFrom(from)) {
                const Vertex to = assignmentGraph.get(ToVertex, edge);
                if (assignmentGraph.hasEdge(to, from)) {
                    if (to < from) continue;
                    vis.drawLine(csaData.stopData[from].coordinates, csaData.stopData[to].coordinates, Color::KITseablue, 8.0);
                } else {
                    vis.drawLine(csaData.stopData[from].coordinates, csaData.stopData[to].coordinates, Color::KITseablue >> 0.5, 8.0);
                }
            }
        }
        vis.drawPoint(csaData.stopData[destination].coordinates, Color::KITred, Icon::MalteseCross, 50);

        for (const Vertex from : assignmentGraph.vertices()) {
            for (const Vertex to : assignmentGraph.outgoingNeighbors(from)) {
                if (assignmentGraph.hasEdge(to, from)) {
                    assignmentGraph.deleteEdge(from, to);
                    assignmentGraph.deleteEdge(to, from);
                }
            }
        }
        vertices = IndexedSet<false, Vertex>(Construct::Complete, assignmentGraph.numVertices());
        while (!vertices.empty()) {
            Vertex vertex = vertices.back();
            vertices.remove(vertex);
            if (assignmentGraph.inDegree(vertex) > 0 && assignmentGraph.outDegree(vertex) > 0) continue;
            for (const Vertex neighbor : assignmentGraph.neighbors(vertex)) {
                vertices.insert(neighbor);
            }
            assignmentGraph.isolateVertex(vertex);
        }
        Graph::printInfo(assignmentGraph);
        assignmentGraph.printAnalysis();

        vis.newPage();
        for (const Border& border : Borders::GermanStates) {
            vis.drawLine(border, Color::Grey, 5, false);
        }
        vis.drawLine(Borders::Germany, Color::Black, 8, false);
        for (const Vertex from : assignmentGraph.vertices()) {
            for (const Edge edge : assignmentGraph.edgesFrom(from)) {
                const Vertex to = assignmentGraph.get(ToVertex, edge);
                if (assignmentGraph.hasEdge(to, from)) {
                    if (to < from) continue;
                    vis.drawLine(csaData.stopData[from].coordinates, csaData.stopData[to].coordinates, Color::KITseablue, 8.0);
                } else {
                    vis.drawLine(csaData.stopData[from].coordinates, csaData.stopData[to].coordinates, Color::KITseablue >> 0.5, 8.0);
                }
            }
        }
        vis.drawPoint(csaData.stopData[destination].coordinates, Color::KITred, Icon::MalteseCross, 50);
        vis.drawPoint(csaData.stopData[destination].coordinates, Color::KITred, Icon::Circle, 50);

    }

};

class TestRandomDistribution : public ParameterizedCommand {

public:
    TestRandomDistribution(BasicShell& shell) :
        ParameterizedCommand(shell, "testRandomDistribution", "tests a random distribution.") {
        addParameter("Distribution", {"Beta", "PERT", "Logit-normal", "Kumaraswamy"});
        addParameter("Number of samples");
        addParameter("Min");
        addParameter("Max");
        addParameter("Mean");
        addParameter("Standard deviation");
        addParameter("Seed", "42");
        addParameter("Sampling factor", "100");
    }

    virtual void execute() noexcept {
        const std::string distribution = getParameter("Distribution");
        const size_t numberOfSamples = getParameter<size_t>("Number of samples");
        const double min = getParameter<double>("Min");
        const double max = getParameter<double>("Max");
        const double mean = getParameter<double>("Mean");
        const double standardDeviation = getParameter<double>("Standard deviation");
        const uint32_t seed = getParameter<uint32_t>("Seed");
        const uint32_t samplingFactor = getParameter<uint32_t>("Sampling factor");

        size_t size = (samplingFactor * max) + 1;
        std::vector<size_t> counts(size, 0);

        if (distribution == "Beta") {
            Assignment::MRUM::Beta beta(min, max, seed);
            for (size_t i = 0; i < numberOfSamples; i++) {
                counts[samplingFactor * beta()]++;
            }
        } else if (distribution == "PERT") {
            Assignment::MRUM::PERT pert(min, max, mean, standardDeviation, seed);
            for (size_t i = 0; i < numberOfSamples; i++) {
                counts[samplingFactor * pert()]++;
            }
        } else if (distribution == "Logit-normal") {
            Assignment::MRUM::LogitNormal logitNormal(min, max, mean, standardDeviation, seed);
            for (size_t i = 0; i < numberOfSamples; i++) {
                counts[samplingFactor * logitNormal()]++;
            }
        } else {
            Assignment::MRUM::Kumaraswamy kumaraswamy(min, max, mean, standardDeviation, seed);
            for (size_t i = 0; i < numberOfSamples; i++) {
                counts[samplingFactor * kumaraswamy()]++;
            }
        }

        for (size_t i = 0; i < size; i++) {
            std::cout << i << "\t" << counts[i] << "\n";
        }
    }

};

class SampleAssignment : public ParameterizedCommand {

public:
    SampleAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "sampleAssignment", "Computes a public transit traffic assignment using a Monte-Carlo approach.",  "Num threads:", "    positive number  - parallel execution with <Num threads> threads", "    otherwise        - sequential execution") {
        addParameter("Settings file");
        addParameter("CSA binary");
        addParameter("Demand file");
        addParameter("Output file");
        addParameter("Num threads", "0");
        addParameter("Thread offset", "1");
    }

    virtual void execute() noexcept {
        const std::string settingsFileName = getParameter("Settings file");

        ConfigFile configFile(settingsFileName, true);
        Assignment::Settings settings(configFile);
        configFile.writeIfModified(false);
        switch (settings.profilerType) {
            case 0: {
                computeAssignment<Assignment::NoProfiler>(settings);
                break;
            }
            case 1: {
                computeAssignment<Assignment::TimeProfiler>(settings);
                break;
            }
            case 2: {
                computeAssignment<Assignment::DecisionProfiler>(settings);
                break;
            }
        }
    }

private:
    template<typename PROFILER>
    inline void computeAssignment(const Assignment::Settings& settings) {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string demandFileName = getParameter("Demand file");
        const std::string outputFileName = getParameter("Output file");
        const int numThreads = getParameter<int>("Num threads");
        const int pinMultiplier = getParameter<int>("Thread offset");

        CSA::Data csaData = CSA::Data::FromBinary(csaFileName);
        csaData.sortConnectionsAscendingByDepartureTime();
        csaData.printInfo();
        csaData.transferGraph.printAnalysis();
        std::cout << std::endl;
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();

        AccumulatedVertexDemand originalDemand = AccumulatedVertexDemand::FromZoneCSV(demandFileName, csaData, reverseGraph);
        IdVertexDemand passengers = IdVertexDemand::FromAccumulatedVertexDemand(originalDemand, settings.passengerMultiplier, settings.demandIntervalSplitTime, settings.includeIntervalBorder);
        std::cout << "Number of passengers: " << String::prettyInt(passengers.numberOfPassengers) << std::endl;
        std::cout << "Number of entires: " << String::prettyInt(passengers.entries.size()) << std::endl;
        std::cout << "Number of ids: " << String::prettyInt(passengers.numIds) << std::endl;

        // Assignment::MRUM::SimpleSampleMeat<100, PROFILER> ssm(csaData, reverseGraph, settings);
        // ssm.sampleRandomUtilities();
        // ssm.run(Vertex(42));

        Assignment::MRUM::Assignment<5, PROFILER> mruma(csaData, reverseGraph, settings);
        Timer timer;
        if (numThreads > 0) {
            const int numCores(numberOfCores());
            std::cout << "Using " << numThreads << " threads on " << numCores << " cores!" << std::endl;
            mruma.run(passengers, numThreads, pinMultiplier);
        } else {
            mruma.run(passengers);
        }
        std::cout << "done in " << String::msToString(timer.elapsedMilliseconds()) << "." << std::endl;
        std::cout << "   removed cycle connections: " << String::prettyInt(mruma.getRemovedCycleConnections()) << std::endl;
        std::cout << "   removed cycles: " << String::prettyInt(mruma.getRemovedCycles()) << std::endl;
        mruma.getProfiler().printStatistics();
        PassengerData result = PassengerData::FromApportionment(csaData, passengers, mruma.getPassengersPerConnection(), mruma.getUnassignedPassengers(), mruma.getWalkingPassengers());
        std::cout << result << std::endl;
        result.serialize(FileSystem::ensureExtension(outputFileName, ".binary"));
        // result.writePassengerConnectionPairs(csaData, passengers, FileSystem::ensureExtension(outputFileName, ".csv"));
        // result.writeCumulativeStopDemand(csaData, FileSystem::ensureExtension(outputFileName, "stopDemand.csv"));
        // mruma.writeConnectionsWithLoad(FileSystem::ensureExtension(outputFileName, "_connections.csv"), passengers.passengerMultiplier);
        // passengers.toCSV(FileSystem::ensureExtension(outputFileName, "_demand.csv"));
    }

};

class CSVStatistic : public ParameterizedCommand {

public:
    CSVStatistic(BasicShell& shell) :
        ParameterizedCommand(shell, "csvStatistic", "Analyzes each column of a .csv file.") {
        addParameter("Input file");
        addParameter("Output file", "");
    }

    virtual void execute() noexcept {
        IO::CSVData<double, IO::TrimChars<>, IO::NoQuoteEscape<'\t'>> inputData(getParameter("Input file"));
        IO::CSVData<std::string, IO::TrimChars<>, IO::NoQuoteEscape<'\t'>> result;
        result.columnNames = std::vector<std::string>{"column", "min", "mean", "median", "max", "sum", "#-1", "#0"};
        result.columnData = std::vector<std::vector<std::string>>(result.columnNames.size(), std::vector<std::string>());
        for (const std::string& column : inputData.columnNames) {
            std::vector<double> data = inputData.getColumn(column);
            std::sort(data.begin(), data.end());
            result.columnData[0].emplace_back(column);
            result.columnData[1].emplace_back(std::to_string(data.front()));
            result.columnData[2].emplace_back(std::to_string(Vector::mean(data)));
            result.columnData[3].emplace_back(std::to_string(Vector::median(data)));
            result.columnData[4].emplace_back(std::to_string(data.back()));
            result.columnData[5].emplace_back(std::to_string(Vector::sum(data)));
            result.columnData[6].emplace_back(std::to_string(Vector::count(data, double(-1))));
            result.columnData[7].emplace_back(std::to_string(Vector::count(data, double(0))));
        }
        result.print();
        const std::string outputFile = getParameter("Output file");
        if (outputFile != "") result.write(outputFile);
    }
};
