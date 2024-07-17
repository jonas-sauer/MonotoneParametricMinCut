#pragma once

#include <string>
#include <vector>

#include "../../Shell/Shell.h"

#include "../../Helpers/String/String.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/MultimodalData.h"
#include "../../DataStructures/TripBased/MultimodalData.h"

#include "../../Algorithms/CH/Preprocessing/CHData.h"
#include "../../Algorithms/CH/Query/PHAST.h"
#include "../../Algorithms/CSA/CSA.h"
#include "../../Algorithms/CSA/DijkstraCSA.h"
#include "../../Algorithms/CSA/OneToAllDijkstraCSA.h"
#include "../../Algorithms/CSA/ProfileCSA.h"
#include "../../Algorithms/CSA/ULTRACSA.h"
#include "../../Algorithms/CSA/UPCSA.h"
#include "../../Algorithms/RAPTOR/AlternatingRAPTOR.h"
#include "../../Algorithms/RAPTOR/Profiler.h"
#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/HLRAPTOR.h"
#include "../../Algorithms/RAPTOR/InitialTransfers.h"
#include "../../Algorithms/RAPTOR/OneToAllDijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/DijkstraRAPTORModule.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/RangeRAPTOR.h"
#include "../../Algorithms/RAPTOR/RangeRAPTOR/UPRangeRAPTOR.h"
#include "../../Algorithms/RAPTOR/RAPTOR.h"
#include "../../Algorithms/RAPTOR/McRAPTOR.h"
#include "../../Algorithms/RAPTOR/MCR.h"
#include "../../Algorithms/RAPTOR/MultimodalMCR.h"
#include "../../Algorithms/RAPTOR/MultimodalULTRAMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMHydRA.h"
#include "../../Algorithms/RAPTOR/ULTRABounded/MultimodalUBMRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRAMcRAPTOR.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"
#include "../../Algorithms/RAPTOR/UPRAPTOR.h"
#include "../../Algorithms/TripBased/Query/HLQuery.h"
#include "../../Algorithms/TripBased/Query/Query.h"
#include "../../Algorithms/TripBased/Query/McQuery.h"
#include "../../Algorithms/TripBased/Query/TransitiveQuery.h"

using namespace Shell;

class RunUnrestrictedRAPTORQuery : public ParameterizedCommand {

public:
    RunUnrestrictedRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUnrestrictedRAPTORQuery", "Runs the given RAPTOR query with unrestricted walking.") {
        addParameter("RAPTOR input file");
        addParameter("Use min transfer times?");
        addParameter("Algorithm", { "regular", "range"});
        addParameter("Initial transfers", { "dijkstra", "core", "bucket"});
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
        addParameter("CH data", "");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        const bool useMinTransferTimes = getParameter<bool>("Use min transfer times?");
        if (!useMinTransferTimes) raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        chooseInitialTransfers(raptorData, useMinTransferTimes);
    }

private:
    struct DijkstraData {
        const RAPTOR::Data& raptorData;
        const RAPTOR::TransferGraph& forwardGraph;
        const RAPTOR::TransferGraph& backwardGraph;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(raptorData, forwardGraph, backwardGraph);
        }
    };

    struct CHData {
        const RAPTOR::Data& raptorData;
        const CH::CH& ch;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(raptorData, ch);
        }
    };

    using ProfilerType = RAPTOR::SimpleProfiler;
    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using RangeAlgorithm = RAPTOR::RangeRAPTOR::DijkstraRAPTORModule<INITIAL_TRANSFERS, ProfilerType, true, false, USE_MIN_TRANSFER_TIMES, true>;
    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using RegularAlgorithm = RAPTOR::DijkstraRAPTOR<INITIAL_TRANSFERS, ProfilerType, true, USE_MIN_TRANSFER_TIMES, false>;

    inline void chooseInitialTransfers(const RAPTOR::Data& raptorData, const bool useMinTransferTimes) const noexcept {
        const std::string initialTransferType = getParameter("Initial transfers");
        if (initialTransferType == "dijkstra") {
            RAPTOR::TransferGraph reverseGraph = raptorData.transferGraph;
            reverseGraph.revert();
            DijkstraData data{ raptorData, raptorData.transferGraph, reverseGraph };
            chooseMinTransferTimes<RAPTOR::DijkstraInitialTransfers>(data, useMinTransferTimes);
        } else {
            CH::CH ch(getParameter("CH data"));
            CHData data{ raptorData, ch };
            if (initialTransferType == "core") {
                chooseMinTransferTimes<RAPTOR::CoreCHInitialTransfers>(data, useMinTransferTimes);
            } else {
                chooseMinTransferTimes<RAPTOR::BucketCHInitialTransfers>(data, useMinTransferTimes);
            }
        }
    }

    template<typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseMinTransferTimes(const DATA_TYPE& data, const bool useMinTransferTimes) const noexcept {
        if (useMinTransferTimes) {
            chooseAlgorithm<true, INITIAL_TRANSFERS>(data);
        } else {
            chooseAlgorithm<false, INITIAL_TRANSFERS>(data);
        }
    }

    template<bool USE_MIN_TRANSFER_TIMES, typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseAlgorithm(const DATA_TYPE& data) const noexcept {
        const std::string algorithmType = getParameter("Algorithm");
        if (algorithmType == "range") {
            run<RangeAlgorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>(data);
        } else {
            run<RegularAlgorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>(data);
        }
    }

    template<typename ALGORITHM, typename DATA_TYPE>
    inline void run(const DATA_TYPE& data) const noexcept {
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        ALGORITHM algorithm = data.template createAlgorithm<ALGORITHM>();
        algorithm.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
            std::cout << journey << std::endl;
        }
    }
};

class RunOneToAllDijkstraRAPTORQuery : public ParameterizedCommand {

public:
    RunOneToAllDijkstraRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runOneToAllDijkstraRAPTORQuery", "Runs the given query with one-to-all Dijkstra-RAPTOR.") {
        addParameter("RAPTOR input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();

        RAPTOR::TransferGraph backwardGraph = raptorData.transferGraph;
        backwardGraph.revert();
        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::AggregateProfiler> algorithm(raptorData, raptorData.transferGraph, backwardGraph);

        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");

        algorithm.run(source, departureTime);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunOneToManyDijkstraRAPTORQuery : public ParameterizedCommand {

public:
    RunOneToManyDijkstraRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runOneToManyDijkstraRAPTORQuery", "Runs the given query with one-to-many Dijkstra-RAPTOR.") {
        addParameter("RAPTOR input file");
        addParameter("Core CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("Core CH data"));
        RAPTOR::OneToAllDijkstraRAPTOR<RAPTOR::CoreCHInitialTransfers, RAPTOR::AggregateProfiler> algorithm(raptorData, ch);

        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");

        algorithm.run(source, departureTime);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunTransitiveRAPTORQuery : public ParameterizedCommand {

public:
    RunTransitiveRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runTransitiveRAPTORQuery", "Runs the given transitive RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("Use min transfer times?");
        addParameter("Transitive?");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        const bool useMinTransferTimes = getParameter<bool>("Use min transfer times?");
        if (!useMinTransferTimes) raptorData.useImplicitDepartureBufferTimes();
        const bool transitive = getParameter<bool>("Transitive?");
        raptorData.printInfo();

        if (useMinTransferTimes) {
            if (transitive) {
                run<true, true>(raptorData);
            }
            else {
                run<true, false>(raptorData);
            }
        } else {
            if (transitive) {
                run<false, true>(raptorData);
            }
            else {
                run<false, false>(raptorData);
            }
        }
    }

private:
    template<bool USE_MIN_TRANSFER_TIMES, bool TRANSITIVE>
    inline void run(const RAPTOR::Data& raptorData) const noexcept {
        RAPTOR::RAPTOR<true, RAPTOR::SimpleProfiler, TRANSITIVE, USE_MIN_TRANSFER_TIMES, false> raptor(raptorData);
        const StopId source = getParameter<StopId>("Source");
        const StopId target = getParameter<StopId>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunTransitiveMcRAPTORQuery : public ParameterizedCommand {

public:
    RunTransitiveMcRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runTransitiveMcRAPTORQuery", "Runs the given transitive McRAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();

        RAPTOR::McRAPTOR<true, true, RAPTOR::SimpleProfiler> raptor(raptorData);
        const StopId source = getParameter<StopId>("Source");
        const StopId target = getParameter<StopId>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunMCRQuery : public ParameterizedCommand {

public:
    RunMCRQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runMCRQuery", "Runs the given MCR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        RAPTOR::MCR<true, RAPTOR::SimpleProfiler> raptor(raptorData, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunMultimodalMCRQuery : public ParameterizedCommand {

public:
    RunMultimodalMCRQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runMultimodalMCRQuery", "Runs the given multimodal MCR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH directory");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        switch (raptorData.modes.size()) {
            case 2:
                run<2>(raptorData);
                break;
            case 3:
                run<3>(raptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

    private:
    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& raptorData) const noexcept {
        const std::string chDirectory(getParameter("CH directory"));
        std::vector<CH::CH> chData;
        for (const size_t mode : raptorData.modes) {
            chData.emplace_back(chDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }

        RAPTOR::MultimodalMCR<true, NUM_MODES, RAPTOR::SimpleProfiler> raptor(raptorData, chData);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunULTRARAPTORQuery : public ParameterizedCommand {

public:
    RunULTRARAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runULTRARAPTORQuery", "Runs the given ULTRA-RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        RAPTOR::ULTRARAPTOR<RAPTOR::SimpleProfiler, false> raptor(raptorData, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunUPRAPTORQuery : public ParameterizedCommand {

public:
    RunUPRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUPRAPTORQuery", "Runs the given UP-RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("Dijkstra graph");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Targets");
        addParameter("Departure time");
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
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
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        TransferGraph dijkstraGraph(getParameter("Dijkstra graph"));
        CH::CH ch(getParameter("CH data"));

        const Vertex source = getParameter<Vertex>("Source");
        const std::string targetsString = getParameter("Targets");
        const std::vector<std::string> targetStrings = String::split(targetsString, ',');
        std::vector<Vertex> targets;
        for (const std::string& targetString : targetStrings) {
            targets.emplace_back(String::lexicalCast<Vertex>(targetString));
        }
        IndexedSet<false, Vertex> targetSet(ch.numVertices(), targets);
        const int departureTime = getParameter<int>("Departure time");

        using Algorithm = RAPTOR::UPRAPTOR<GROUPED_ROUNDS, RAPTOR::SimpleProfiler>;
        const bool reorder = getParameter<bool>("Reorder network?");
        Algorithm raptor(raptorData, dijkstraGraph, ch, targetSet, reorder, getParameter("Order") == "DFS");

        raptor.run(source, departureTime);
        for (const Vertex target : targetSet) {
            std::cout << "Journeys to " << target << ":" << std::endl;
            for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
                std::cout << journey << std::endl;
            }
            std::cout << std::endl;
        }
    }
};

class RunULTRAMcRAPTORQuery : public ParameterizedCommand {

public:
    RunULTRAMcRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runULTRAMcRAPTORQuery", "Runs the given ULTRA-McRAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        RAPTOR::ULTRAMcRAPTOR<RAPTOR::SimpleProfiler> raptor(raptorData, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunMultimodalULTRAMcRAPTORQuery : public ParameterizedCommand {

public:
    RunMultimodalULTRAMcRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runMultimodalULTRAMcRAPTORQuery", "Runs the given multimodal ULTRA-McRAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH directory");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        switch (raptorData.modes.size()) {
            case 2:
                run<2>(raptorData);
                break;
            case 3:
                run<3>(raptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

    private:
    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& raptorData) const noexcept {
        const std::string chDirectory(getParameter("CH directory"));
        std::vector<CH::CH> chData;
        for (const size_t mode : raptorData.modes) {
            chData.emplace_back(chDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }

        RAPTOR::MultimodalULTRAMcRAPTOR<NUM_MODES, RAPTOR::SimpleProfiler> raptor(raptorData, chData);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunMultimodalUBMRAPTORQuery : public ParameterizedCommand {

public:
    RunMultimodalUBMRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runMultimodalUBMRAPTORQuery", "Runs the given multimodal UBM-RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("Bucket-CH directory");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        RAPTOR::MultimodalData raptorData(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        switch (raptorData.modes.size()) {
            case 2:
                run<2>(raptorData);
                break;
            case 3:
                run<3>(raptorData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

    private:
    template<size_t NUM_MODES>
    inline void run(const RAPTOR::MultimodalData& raptorData) const noexcept {
        const RAPTOR::Data pruningData = raptorData.getPruningData();
        const RAPTOR::Data reversePruningData = pruningData.reverseNetwork();
        const std::string bucketCHDirectory(getParameter("Bucket-CH directory"));
        std::vector<CH::CH> bucketCHData;
        for (const size_t mode : raptorData.modes) {
            bucketCHData.emplace_back(bucketCHDirectory + RAPTOR::TransferModeNames[mode] + "CH");
        }
        RAPTOR::TransferGraph backwardTransitiveGraph = raptorData.raptorData.transferGraph;
        backwardTransitiveGraph.revert();

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::MultimodalUBMRAPTOR<NUM_MODES, RAPTOR::AggregateProfiler> algorithm(raptorData, pruningData, reversePruningData, backwardTransitiveGraph, bucketCHData);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        algorithm.run(source, departureTime, target, arrivalSlack, tripSlack);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunMultimodalUBMHydRAQuery : public ParameterizedCommand {

public:
    RunMultimodalUBMHydRAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runMultimodalUBMHydRAQuery", "Runs the given multimodal UBM-HydRA query.") {
        addParameter("ULTRA-Trip-Based input file");
        addParameter("Bounded forward Trip-Based input file");
        addParameter("Bounded backward Trip-Based input file");
        addParameter("Bucket-CH directory");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
        addParameter("Arrival slack");
        addParameter("Trip slack");
    }

    virtual void execute() noexcept {
        const TripBased::MultimodalData tripBasedData(getParameter("ULTRA-Trip-Based input file"));
        tripBasedData.printInfo();
        switch (tripBasedData.modes.size()) {
            case 2:
                run<2>(tripBasedData);
                break;
            case 3:
                run<3>(tripBasedData);
                break;
            default:
                Ensure(false, "Unsupported number of modes!");
                break;
        }
    }

    private:
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

        const double arrivalSlack = getParameter<double>("Arrival slack");
        const double tripSlack = getParameter<double>("Trip slack");

        RAPTOR::MultimodalUBMHydRA<NUM_MODES, RAPTOR::AggregateProfiler> hydra(ultraTripBasedData, forwardPruningData, backwardPruningData, backwardTransitiveGraph, bucketCHData);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        hydra.run(source, departureTime, target, arrivalSlack, tripSlack);
        for (const RAPTOR::Journey& journey : hydra.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunHLRAPTORQuery : public ParameterizedCommand {

public:
    RunHLRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runHLRAPTORQuery", "Run the given HL-RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));

        RAPTOR::HLRAPTOR<RAPTOR::AggregateProfiler> raptor(raptorData, outHubs, inHubs);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        raptor.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : raptor.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunTransitiveTripBasedQuery : public ParameterizedCommand {

public:
    RunTransitiveTripBasedQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runTransitiveTripBasedQuery", "Runs the given transitive Trip-Based query.") {
        addParameter("Trip-Based input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        TripBased::Data data(getParameter("Trip-Based input file"));
        data.printInfo();
        TripBased::TransitiveQuery<TripBased::AggregateProfiler> algorithm(data);
        const StopId source = getParameter<StopId>("Source");
        const StopId target = getParameter<StopId>("Target");
        const int departureTime = getParameter<int>("Departure time");
        algorithm.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
            std::cout << journey << std::endl;
        }
    }
};

class RunULTRATripBasedQuery : public ParameterizedCommand {

public:
    RunULTRATripBasedQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runULTRATripBasedQuery", "Runs the given ULTRA-Trip-Based query.") {
        addParameter("Trip-Based input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        TripBased::Data data(getParameter("Trip-Based input file"));
        data.printInfo();
        CH::CH ch(getParameter("CH data"));
        TripBased::Query<TripBased::AggregateProfiler> algorithm(data, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        algorithm.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
            std::cout << journey << std::endl;
        }
    }
};

class RunHLULTRATripBasedQuery : public ParameterizedCommand {

public:
    RunHLULTRATripBasedQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runHLULTRATripBasedQuery", "Runs the given HL-ULTRA-Trip-Based query.") {
        addParameter("Trip-Based input file");
        addParameter("Out-hub file");
        addParameter("In-hub file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        TripBased::Data data(getParameter("Trip-Based input file"));
        data.printInfo();
        const TransferGraph outHubs(getParameter("Out-hub file"));
        const TransferGraph inHubs(getParameter("In-hub file"));
        TripBased::HLQuery<TripBased::AggregateProfiler> algorithm(data, outHubs, inHubs);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        algorithm.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
            std::cout << journey << std::endl;
        }
    }
};

class RunULTRAMcTripBasedQuery : public ParameterizedCommand {

public:
    RunULTRAMcTripBasedQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runULTRAMcTripBasedQuery", "Runs the given ULTRA-McTrip-Based query.") {
        addParameter("Trip-Based input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        TripBased::Data data(getParameter("Trip-Based input file"));
        data.printInfo();
        CH::CH ch(getParameter("CH data"));
        TripBased::McQuery<TripBased::AggregateProfiler> algorithm(data, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        algorithm.run(source, departureTime, target);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
            std::cout << journey << std::endl;
        }
    }
};

class RunUnrestrictedProfileRAPTORQuery : public ParameterizedCommand {

public:
    RunUnrestrictedProfileRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUnrestrictedProfileRAPTORQuery", "Runs the given Alternating-/Range-RAPTOR query with unrestricted walking.") {
        addParameter("RAPTOR input file");
        addParameter("Use min transfer times?");
        addParameter("Algorithm", { "alternating", "range" });
        addParameter("Initial transfers", { "dijkstra", "core", "bucket"});
        addParameter("Source");
        addParameter("Target");
        addParameter("Min departure time", "0");
        addParameter("Max departure time", "86400");
        addParameter("CH data", "");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        const bool useMinTransferTimes = getParameter<bool>("Use min transfer times?");
        if (!useMinTransferTimes) raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        chooseInitialTransfers(raptorData, useMinTransferTimes);
    }

private:
    struct DijkstraData {
        const RAPTOR::Data& forwardRaptorData;
        const RAPTOR::Data& backwardRaptorData;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(forwardRaptorData, backwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);
        }
    };

    struct CHData {
        const RAPTOR::Data& forwardRaptorData;
        const RAPTOR::Data& backwardRaptorData;
        const CH::CH& ch;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(forwardRaptorData, backwardRaptorData, ch);
        }
    };

    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using RangeRAPTOR = RAPTOR::RangeRAPTOR::RangeDijkstraRAPTOR<INITIAL_TRANSFERS, true, true, USE_MIN_TRANSFER_TIMES, true, true>;
    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using AlternatingRAPTOR = RAPTOR::AlternatingDijkstraRAPTOR<INITIAL_TRANSFERS, true, USE_MIN_TRANSFER_TIMES>;

    inline void chooseInitialTransfers(const RAPTOR::Data& raptorData, const bool useMinTransferTimes) const noexcept {
        const std::string initialTransferType = getParameter("Initial transfers");
        RAPTOR::Data reverseData = raptorData.reverseNetwork();
        if (initialTransferType == "dijkstra") {
            DijkstraData data{ raptorData, reverseData };
            chooseMinTransferTimes<RAPTOR::DijkstraInitialTransfers>(data, useMinTransferTimes);
        } else {
            CH::CH ch(getParameter("CH data"));
            CHData data{ raptorData, reverseData, ch };
            if (initialTransferType == "core") {
                chooseMinTransferTimes<RAPTOR::CoreCHInitialTransfers>(data, useMinTransferTimes);
            } else {
                chooseMinTransferTimes<RAPTOR::BucketCHInitialTransfers>(data, useMinTransferTimes);
            }
        }
    }

    template<typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseMinTransferTimes(const DATA_TYPE& data, const bool useMinTransferTimes) const noexcept {
        if (useMinTransferTimes) {
            chooseAlgorithm<true, INITIAL_TRANSFERS>(data);
        } else {
            chooseAlgorithm<false, INITIAL_TRANSFERS>(data);
        }
    }

    template<bool USE_MIN_TRANSFER_TIMES, typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseAlgorithm(const DATA_TYPE& data) const noexcept {
        const std::string algorithmType = getParameter("Algorithm");
        if (algorithmType == "range") {
            run<RangeRAPTOR<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>(data);
        } else {
            run<AlternatingRAPTOR<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>(data);
        }
    }

    template<typename ALGORITHM, typename DATA_TYPE>
    inline void run(const DATA_TYPE& data) const noexcept {
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int minDepartureTime = getParameter<int>("Min departure time");
        const int maxDepartureTime = getParameter<int>("Max departure time");
        ALGORITHM algorithm = data.template createAlgorithm<ALGORITHM>();

        Timer timer;
        timer.restart();
        algorithm.run(source, target, minDepartureTime, maxDepartureTime);
        std::cout << "Total time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        RAPTOR::ProfileHandle profile = algorithm.getProfileHandle();
        profile.sort();
        std::cout << profile;
    }
};

class RunOneToAllRAPTORQuery : public ParameterizedCommand {

public:
    RunOneToAllRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runOneToAllRAPTORQuery", "Runs the given one-to-all RAPTOR query and outputs the journeys to target.") {
        addParameter("RAPTOR input file");
        addParameter("Use min transfer times?");
        addParameter("Initial transfers", { "dijkstra", "core", "bucket"});
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
        addParameter("CH data", "");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        const bool useMinTransferTimes = getParameter<bool>("Use min transfer times?");
        if (!useMinTransferTimes) raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        chooseInitialTransfers(raptorData, useMinTransferTimes);
    }

private:
    struct DijkstraData {
        const RAPTOR::Data& forwardRaptorData;
        const RAPTOR::Data& backwardRaptorData;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(forwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);
        }
    };

    struct CHData {
        const RAPTOR::Data& raptorData;
        const CH::CH& ch;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(raptorData, ch);
        }
    };

    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using Algorithm = RAPTOR::RangeRAPTOR::DijkstraRAPTORModule<INITIAL_TRANSFERS, RAPTOR::SimpleProfiler, false, false, USE_MIN_TRANSFER_TIMES, false>;

    inline void chooseInitialTransfers(const RAPTOR::Data& raptorData, const bool useMinTransferTimes) const noexcept {
        const std::string initialTransferType = getParameter("Initial transfers");
        if (initialTransferType == "dijkstra") {
            RAPTOR::Data reverseData = raptorData.reverseNetwork();
            DijkstraData data{ raptorData, reverseData };
            chooseMinTransferTimes<RAPTOR::DijkstraInitialTransfers>(data, useMinTransferTimes);
        } else {
            CH::CH ch(getParameter("CH data"));
            CHData data{ raptorData, ch };
            if (initialTransferType == "core") {
                chooseMinTransferTimes<RAPTOR::CoreCHInitialTransfers>(data, useMinTransferTimes);
            } else {
                chooseMinTransferTimes<RAPTOR::BucketCHInitialTransfers>(data, useMinTransferTimes);
            }
        }
    }

    template<typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseMinTransferTimes(const DATA_TYPE& data, const bool useMinTransferTimes) const noexcept {
        if (useMinTransferTimes) {
            run<true, INITIAL_TRANSFERS>(data);
        } else {
            run<false, INITIAL_TRANSFERS>(data);
        }
    }

    template<bool USE_MIN_TRANSFER_TIMES, typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void run(const DATA_TYPE& data) const noexcept {
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        Algorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES> algorithm = data.template createAlgorithm<Algorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>();

        algorithm.run(source, departureTime);
        for (const RAPTOR::Journey& journey : algorithm.getJourneys(target)) {
            std::cout << journey << std::endl;
        }
    }
};

class RunOneToAllRangeRAPTORQuery : public ParameterizedCommand {

public:
    RunOneToAllRangeRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runOneToAllRangeRAPTORQuery", "Runs the given one-to-all Range-RAPTOR query and outputs the profile at target.") {
        addParameter("RAPTOR input file");
        addParameter("Use min transfer times?");
        addParameter("Initial transfers", { "dijkstra", "core", "bucket"});
        addParameter("Source");
        addParameter("Target");
        addParameter("Min departure time", "0");
        addParameter("Max departure time", "86400");
        addParameter("CH data", "");
    }

    virtual void execute() noexcept {
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        const bool useMinTransferTimes = getParameter<bool>("Use min transfer times?");
        if (!useMinTransferTimes) raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        chooseInitialTransfers(raptorData, useMinTransferTimes);
    }

private:
    struct DijkstraData {
        const RAPTOR::Data& forwardRaptorData;
        const RAPTOR::Data& backwardRaptorData;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(forwardRaptorData, backwardRaptorData);
        }
    };

    struct CHData {
        const RAPTOR::Data& forwardRaptorData;
        const RAPTOR::Data& backwardRaptorData;
        const CH::CH& ch;

        template<typename ALGORITHM>
        inline ALGORITHM createAlgorithm() const noexcept {
            return ALGORITHM(forwardRaptorData, backwardRaptorData, ch);
        }
    };

    template<typename INITIAL_TRANSFERS, bool USE_MIN_TRANSFER_TIMES>
    using Algorithm = RAPTOR::RangeRAPTOR::RangeDijkstraRAPTOR<INITIAL_TRANSFERS, false, false, USE_MIN_TRANSFER_TIMES, true, false>;

    inline void chooseInitialTransfers(const RAPTOR::Data& raptorData, const bool useMinTransferTimes) const noexcept {
        const std::string initialTransferType = getParameter("Initial transfers");
        RAPTOR::Data reverseData = raptorData.reverseNetwork();
        if (initialTransferType == "dijkstra") {
            DijkstraData data{ raptorData, reverseData };
            chooseMinTransferTimes<RAPTOR::DijkstraInitialTransfers>(data, useMinTransferTimes);
        } else {
            CH::CH ch(getParameter("CH data"));
            CHData data{ raptorData, reverseData, ch };
            if (initialTransferType == "core") {
                chooseMinTransferTimes<RAPTOR::CoreCHInitialTransfers>(data, useMinTransferTimes);
            } else {
                chooseMinTransferTimes<RAPTOR::BucketCHInitialTransfers>(data, useMinTransferTimes);
            }
        }
    }

    template<typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void chooseMinTransferTimes(const DATA_TYPE& data, const bool useMinTransferTimes) const noexcept {
        if (useMinTransferTimes) {
            run<true, INITIAL_TRANSFERS>(data);
        } else {
            run<false, INITIAL_TRANSFERS>(data);
        }
    }

    template<bool USE_MIN_TRANSFER_TIMES, typename INITIAL_TRANSFERS, typename DATA_TYPE>
    inline void run(const DATA_TYPE& data) const noexcept {
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int minDepartureTime = getParameter<int>("Min departure time");
        const int maxDepartureTime = getParameter<int>("Max departure time");
        Algorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES> algorithm = data.template createAlgorithm<Algorithm<INITIAL_TRANSFERS, USE_MIN_TRANSFER_TIMES>>();

        Timer timer;
        timer.restart();
        algorithm.runOneToAll(source, minDepartureTime, maxDepartureTime);
        std::cout << "Total time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        RAPTOR::ProfileHandle profile = algorithm.getProfileHandle(target);
        profile.sort();
        std::cout << profile;
    }
};

class RunUPRangeRAPTORQuery : public ParameterizedCommand {

public:
    RunUPRangeRAPTORQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUPRangeRAPTORQuery", "Runs the given UP-Range-RAPTOR query.") {
        addParameter("RAPTOR input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Targets");
        addParameter("Min departure time");
        addParameter("Max departure time");
        addParameter("Final transfers", { "Bucket", "PHAST" });
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
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
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(getParameter("RAPTOR input file"));
        raptorData.useImplicitDepartureBufferTimes();
        raptorData.printInfo();
        CH::CH ch(getParameter("CH data"));

        const Vertex source = getParameter<Vertex>("Source");
        IndexedSet<false, Vertex> targetSet(Construct::Complete, ch.numVertices());

        const std::string targetsString = getParameter("Targets");
        const std::vector<std::string> targetStrings = String::split(targetsString, ',');
        std::vector<Vertex> targets;
        for (const std::string& targetString : targetStrings) {
            targets.emplace_back(String::lexicalCast<Vertex>(targetString));
        }

        const int minDepartureTime = getParameter<int>("Min departure time");
        const int maxDepartureTime = getParameter<int>("Max departure time");

        using Algorithm = RAPTOR::RangeRAPTOR::UPRangeRAPTOR<USE_TARGET_BUCKETS, USE_DFS_ORDER, true, true>;
        Algorithm raptor(raptorData, ch, targetSet, getParameter<bool>("Reorder network?"));
        raptor.run(source, minDepartureTime, maxDepartureTime);
        for (const Vertex target : targets) {
            std::cout << "Profile of " << target << ":" << std::endl;
            std::cout << raptor.getProfile(target) << std::endl;
        }
    }
};

class RunUnrestrictedCSAQuery : public ParameterizedCommand {

public:
    RunUnrestrictedCSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUnrestrictedCSAQuery", "Runs the given CSA query with unrestricted walking.") {
        addParameter("CSA input file");
        addParameter("Initial transfers", { "dijkstra", "core", "bucket"});
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
        addParameter("CH data", "");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        std::string initialTransferType;
        chooseInitialTransfers(csaData);
    }

private:
    template<typename INITIAL_TRANSFERS>
    using Algorithm = CSA::DijkstraCSA<INITIAL_TRANSFERS, true, CSA::SimpleProfiler>;

    inline void chooseInitialTransfers(const CSA::Data& csaData) const noexcept {
        const std::string initialTransferType = getParameter("Initial transfers");
        if (initialTransferType == "dijkstra") {
            CSA::TransferGraph reverseGraph = csaData.transferGraph;
            reverseGraph.revert();
            Algorithm<RAPTOR::DijkstraInitialTransfers> csa(csaData, csaData.transferGraph, reverseGraph);
            run(csa);
        } else {
            CH::CH ch(getParameter("CH data"));
            if (initialTransferType == "core") {
                Algorithm<RAPTOR::CoreCHInitialTransfers> csa(csaData, ch);
                run(csa);
            } else {
                Algorithm<RAPTOR::BucketCHInitialTransfers> csa(csaData, ch);
                run(csa);
            }
        }
    }

    template<typename ALGORITHM>
    inline void run(ALGORITHM& csa) const noexcept {
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        csa.run(source, departureTime, target);
        std::cout << csa.getJourney() << std::endl;
    }
};

class RunTransitiveCSAQuery : public ParameterizedCommand {

public:
    RunTransitiveCSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runTransitiveCSAQuery", "Runs the given transitive CSA query.") {
        addParameter("CSA input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CSA::CSA<true, CSA::SimpleProfiler> csa(csaData);
        const StopId source = getParameter<StopId>("Source");
        const StopId target = getParameter<StopId>("Target");
        const int departureTime = getParameter<int>("Departure time");
        csa.run(source, departureTime, target);
        std::cout << csa.getJourney();
    }
};

class RunTransitiveProfileCSAQuery : public ParameterizedCommand {

public:
    RunTransitiveProfileCSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runTransitiveProfileCSAQuery", "Runs the given transitive Profile-CSA query.") {
        addParameter("CSA input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Contiguous profiles?");
        addParameter("Min departure time", "0");
        addParameter("Max departure time", "86400");
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
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CSA::TransferGraph reverseGraph = csaData.transferGraph;
        reverseGraph.revert();

        const StopId source = getParameter<StopId>("Source");
        const StopId target = getParameter<StopId>("Target");
        const int minDepartureTime = getParameter<int>("Min departure time");
        const int maxDepartureTime = getParameter<int>("Max departure time");
        CSA::ProfileCSA<true, CSA::SimpleProfiler, CONTIGUOUS_PROFILES> algorithm(csaData, reverseGraph);

        Timer timer;
        timer.restart();
        algorithm.run(source, target, minDepartureTime, maxDepartureTime);
        std::cout << "Total time: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
        RAPTOR::ProfileHandle profile = algorithm.getProfileHandle(source);
        profile.sort();
        std::cout << profile;
    }
};

class RunULTRACSAQuery : public ParameterizedCommand {

public:
    RunULTRACSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runULTRACSAQuery", "Runs the given ULTRA-CSA query.") {
        addParameter("CSA input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CH::CH ch(getParameter("CH data"));

        CSA::ULTRACSA<true, CSA::SimpleProfiler> csa(csaData, ch);
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        csa.run(source, departureTime, target);
        std::cout << csa.getJourney();
    }
};

class RunUPCSAQuery : public ParameterizedCommand {

public:
    RunUPCSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runUPCSAQuery", "Runs the given UP-CSA query.") {
        addParameter("CSA input file");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Targets");
        addParameter("Departure time");
        addParameter("Initial transfers", { "Bucket", "PHAST" });
        addParameter("Final transfers", { "Bucket", "PHAST" });
        addParameter("Reorder network?");
        addParameter("Order", { "DFS", "Level" });
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
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        CH::CH ch(getParameter("CH data"));

        const Vertex source = getParameter<Vertex>("Source");
        const std::string targetsString = getParameter("Targets");
        const std::vector<std::string> targetStrings = String::split(targetsString, ',');
        std::vector<Vertex> targets;
        for (const std::string& targetString : targetStrings) {
            targets.emplace_back(String::lexicalCast<Vertex>(targetString));
        }
        IndexedSet<false, Vertex> targetSet(ch.numVertices(), targets);
        const int departureTime = getParameter<int>("Departure time");

        const bool reorder = getParameter<bool>("Reorder network?");
        using Algorithm = CSA::UPCSA<USE_STOP_BUCKETS, USE_TARGET_BUCKETS, true, CSA::SimpleProfiler>;
        Algorithm csa(csaData, ch, targetSet, reorder, getParameter("Order") == "DFS");

        csa.run(source, departureTime);
        for (const Vertex target : targetSet) {
            std::cout << "Journey to " << target << ":" << std::endl;
            std::cout << csa.getJourney(target) << std::endl;
        }
    }
};

class RunOneToAllDijkstraCSAQuery : public ParameterizedCommand {

public:
    RunOneToAllDijkstraCSAQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runOneToAllDijkstraCSAQuery", "Runs the given one-to-all Dijkstra CSA query and outputs the journey to target.") {
        addParameter("CSA input file");
        addParameter("Source");
        addParameter("Target");
        addParameter("Departure time");
    }

    virtual void execute() noexcept {
        CSA::Data csaData = CSA::Data::FromBinary(getParameter("CSA input file"));
        csaData.sortConnectionsAscending();
        csaData.printInfo();
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure time");
        CSA::OneToAllDijkstraCSA<CSA::TransferGraph, true, CSA::SimpleProfiler> algorithm(csaData);
        algorithm.run(source, departureTime);
        std::cout << algorithm.getJourney(target) << std::endl;
    }
};

class RunDijkstraQuery : public ParameterizedCommand {

public:
    RunDijkstraQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runDijkstraQuery", "Runs the given Dijkstra query.") {
        addParameter("Graph");
        addParameter("Source");
        addParameter("Target");
    }

    virtual void execute() noexcept {
        CSA::TransferGraph graph(getParameter("Graph"));
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");

        Dijkstra<CSA::TransferGraph> dijkstra(graph);
        dijkstra.run(source, target);
        std::cout << "Distance: " << dijkstra.getDistance(target) << std::endl;
        std::cout << "Path: ";
        Vector::printConcise(dijkstra.getPath(target));
        std::cout << std::endl;
    }
};


class RunCHQuery : public ParameterizedCommand {

public:
    RunCHQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runCHQuery", "Runs the given CH query.") {
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
    }

    virtual void execute() noexcept {
        CH::CH ch(getParameter("CH data"));
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");

        CH::Query<CHGraph, true, true, false> chQuery(ch);
        chQuery.run(source, target);
        std::cout << "Distance: " << chQuery.getDistance() << std::endl;
        std::cout << "Intersecting vertex: " << chQuery.getIntersectingVertex() << std::endl;
        std::cout << "Path: ";
        Vector::printConcise(chQuery.getPath());
        std::cout << std::endl;
        std::cout << "Forward leg: ";
        Vector::printConcise(chQuery.getPackedForwardLeg());
        std::cout << std::endl;
        std::cout << "Backward leg: ";
        Vector::printConcise(chQuery.getPackedBackwardLeg());
        std::cout << std::endl;
    }
};

class RunPHASTQuery : public ParameterizedCommand {

public:
    RunPHASTQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runPHASTQuery", "Runs the given PHAST query.") {
        addParameter("CH data");
        addParameter("Source");
        addParameter("Order", { "Contraction", "Level"});
        addParameter("Reorder vertices?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Reorder vertices?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

    template<bool REORDER_VERTICES>
    inline void run() const noexcept {
        CH::CH ch(getParameter("CH data"));
        CH::PHAST<REORDER_VERTICES, true> phastQuery(ch, getOrder(ch));
        phastQuery.run(getParameter<Vertex>("Source"));
    }

private:
    inline Order getOrder(const CH::CH& ch) const noexcept {
        if (getParameter("Order") == "Contraction") {
            return Order(Vector::reverse(CH::getOrder(ch)));
        } else {
            return Order(CH::getLevelOrderTopDown(ch));
        }
    }
};
