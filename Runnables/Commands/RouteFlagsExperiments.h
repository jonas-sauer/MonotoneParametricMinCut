#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../../DataStructures/Partition/NestedDissection.h"
#include "../../DataStructures/Partition/VertexPartition.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"
#include "../../Algorithms/CH/Preprocessing/CHBuilder.h"
#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/RouteFlags/Preprocessing/Preprocessing.h"
#include "../../Algorithms/RAPTOR/RouteFlags/Query/Query.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Permutation.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Shell/Shell.h"
#include "../../Visualization/MapVisualization.h"
#include "../../Visualization/PDF.h"
#include "../../Visualization/PNG.h"

using namespace Shell;

const Color highlightColor = Color::Red;
const Color regularColor = Color::Grey;

inline static Color getColor(const bool highlight) noexcept {
    return highlight ? highlightColor : regularColor;
}

template<typename FLAGS, typename FINE_PARTITION>
using RouteFlagsDataType = RAPTOR::RouteFlags::RouteFlagsData<FLAGS, FINE_PARTITION>;

inline constexpr int ShortcutWeight = 1024;
inline constexpr int LevelWeight = 256;
inline constexpr int DegreeWeight = 0;
inline constexpr int BidirectionalPopLimit = 200;

using CHProfiler = CH::FullProfiler;
using WitnessSearch = CH::BidirectionalWitnessSearch<CHCoreGraph, CHProfiler, BidirectionalPopLimit>;
using GreedyKey = CH::GreedyKey<WitnessSearch>;
using PartialKey = CH::PartialKey<WitnessSearch, GreedyKey>;
using StopCriterion = CH::CoreCriterion;
using CHBuilder = CH::Builder<CHProfiler, WitnessSearch, PartialKey, StopCriterion, false, false>;

class BuildPartitionCoreCH : public ParameterizedCommand {

public:
    BuildPartitionCoreCH(BasicShell& shell) :
        ParameterizedCommand(shell, "buildPartitionCoreCH", "Computes a core-CH for the input network and partition, where all stops and separator vertices are kept uncontracted.") {
        addParameter("RAPTOR binary");
        addParameter("Coarse partition");
        addParameter("Coarse level");
        addParameter("Fine partition");
        addParameter("Fine level");
        addParameter("Max core degree");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        RAPTOR::Data data(getParameter("RAPTOR binary"));
        data.printInfo();
        const size_t numberOfVertices = data.transferGraph.numVertices();
        Intermediate::TransferGraph transferGraph;
        transferGraph.addVertices(numberOfVertices);
        for (const Vertex vertex : transferGraph.vertices()) {
            transferGraph.set(Coordinates, vertex, data.transferGraph.get(Coordinates, vertex));
        }
        std::vector<Vertex> coarseSeparator = getCoarseSeparator();
        std::cout << "Coarse separator size: " << String::prettyInt(coarseSeparator.size()) << std::endl;
        std::vector<Vertex> fineSeparator = getFineSeparator(data);
        std::cout << "Fine separator size: " << String::prettyInt(fineSeparator.size()) << std::endl;
        std::vector<bool> isContractable(numberOfVertices, true);
        std::vector<int> group(numberOfVertices, 2);
        for (const Vertex vertex : coarseSeparator) {
            group[vertex] = 1;
            isContractable[vertex] = false;
        }
        for (const Vertex vertex : fineSeparator) {
            group[vertex] = 1;
            isContractable[vertex] = false;
        }
        for (const StopId stop : data.stops()) {
            group[stop] = 0;
            isContractable[stop] = false;
        }
        const size_t minCoreSize = Vector::count(isContractable, false);
        const double maxCoreDegree = getParameter<double>("Max core degree");
        std::cout << "Min. core size: " << String::prettyInt(minCoreSize) << std::endl;
        std::cout << "Max. core degree: " << String::prettyInt(maxCoreDegree) << std::endl;

        GreedyKey greedyKey(ShortcutWeight, LevelWeight, DegreeWeight);
        PartialKey keyFunction(isContractable, numberOfVertices, greedyKey);
        StopCriterion stopCriterion(minCoreSize, maxCoreDegree);
        CHBuilder chBuilder(std::move(data.transferGraph), data.transferGraph[TravelTime], keyFunction, stopCriterion);
        chBuilder.run();
        chBuilder.copyCoreToCH();
        std::cout << "Obtaining CH" << std::endl;
        CH::CH ch(std::move(chBuilder));
        std::cout << "Obtaining Order" << std::endl;
        size_t coreSize = numberOfVertices;
        for (const Vertex vertex : chBuilder.getOrder()) {
            AssertMsg(group[vertex] == 2, "Vertex " << vertex << " is not contractable but was contracted!");
            group[vertex] = 3;
            coreSize--;
        }
        std::cout << "Final core size: " << String::prettyInt(coreSize) << std::endl;
        std::cout << "Replacing RAPTOR transfer graph" << std::endl;
        for (const Vertex vertex : transferGraph.vertices()) {
            for (const Edge edge : ch.forward.edgesFrom(vertex)) {
                const Edge newEdge = transferGraph.findOrAddEdge(vertex, ch.forward.get(ToVertex, edge));
                transferGraph.set(TravelTime, newEdge, ch.forward.get(Weight, edge));
            }
        }
        Graph::move(std::move(transferGraph), data.transferGraph);
        std::cout << "Reordering vertices" << std::endl;
        Order order(Construct::Id, numberOfVertices);
        sort(order, [&](const size_t a, const size_t b){
            return (group[a] < group[b]) || ((group[a] == group[b]) && (a < b));
        });
        ch.applyVertexOrder(order);
        data.transferGraph.applyVertexOrder(order);
        data.transferGraph.deleteVertices([&](const Vertex vertex){return vertex >= coreSize;});
        CH::analyze(ch);
        data.printInfo();
        Graph::printInfo(data.transferGraph);
        data.transferGraph.printAnalysis();
        std::cout << "Writing results" << std::endl;
        const std::string outputFile = getParameter("Output file");
        order.serialize(outputFile + ".vertexOrder");
        ch.writeBinary(outputFile + ".ch");
        data.serialize(outputFile + ".raptor.binary");
        writeCoarseSeparator(order);
        writeFineSeparator(order);
    }

private:
    inline std::vector<Vertex> getCoarseSeparator() noexcept {
        coarseND.deserialize(getParameter("Coarse partition"));
        return coarseND.getSeparatorOfLevel(getParameter<int>("Coarse level"));
    }

    inline std::vector<Vertex> getFineSeparator(const RAPTOR::Data& data) noexcept {
        const std::string finePartition = getParameter("Fine partition");

        if (getParameter("Fine level") == "-") {
            fineVP.deserialize(finePartition);
            return fineVP.getBorderVertices(data.transferGraph);
        } else {
            fineND.deserialize(finePartition);
            return fineND.getSeparatorOfLevel(getParameter<int>("Fine level"));
        }
    }
    inline void writeCoarseSeparator(const Order& order) noexcept {
        coarseND.applyVertexOrder(order);
        coarseND.serialize(getParameter("Output file") + ".coarse.nd");
    }

    inline void writeFineSeparator(const Order& order) noexcept {
        const std::string outputFile = getParameter("Output file");

        if (getParameter("Fine level") == "-") {
            fineVP.applyVertexOrder(order);
            fineVP.serialize(outputFile + ".fine.vp");
        } else {
            fineND.applyVertexOrder(order);
            fineND.serialize(outputFile + ".fine.nd");
        }
    }

private:
    NestedDissection coarseND;
    NestedDissection fineND;
    VertexPartition fineVP;

};

class BuildRouteFlagsData : public ParameterizedCommand {

public:
    BuildRouteFlagsData(BasicShell& shell) :
        ParameterizedCommand(shell, "buildRouteFlagsData", "Builds Route-Flags data structure from the given data.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("RAPTOR network");
        addParameter("Shortcut graph");
        addParameter("Coarse partition data");
        addParameter("Coarse level");
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Fine partition data");
        addParameter("Output file");
        addParameter("Fine level", "0");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            readData<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            readData<RAPTOR::RouteFlags::TripFlags>();
        } else {
            readData<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void readData() const noexcept {
        const std::string raptorPath = getParameter("RAPTOR network");
        RAPTOR::Data forwardRaptorData = RAPTOR::Data::FromBinary(raptorPath);
        forwardRaptorData.useImplicitDepartureBufferTimes();
        Permutation stopEventPermutation;
        RAPTOR::Data backwardRaptorData = forwardRaptorData.reverseNetwork(stopEventPermutation);

        const std::string shortcutGraphPath = getParameter("Shortcut graph");
        RAPTOR::TransferGraph forwardShortcutGraph;
        forwardShortcutGraph.readBinary(shortcutGraphPath);
        RAPTOR::TransferGraph backwardShortcutGraph = forwardShortcutGraph;
        backwardShortcutGraph.revert();

        const std::string coarsePartitionPath = getParameter("Coarse partition data");
        NestedDissection coarsePartition(coarsePartitionPath);
        const size_t coarseLevel = getParameter<size_t>("Coarse level");

        const std::string finePartitionType = getParameter("Fine partition type");
        std::string finePartitionPath = getParameter("Fine partition data");
        const size_t fineLevel = getParameter<size_t>("Fine level");
        if (finePartitionType == "vp") {
            VertexPartition finePartition(finePartitionPath);
            RouteFlagsDataType<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition> routeFlagsData(forwardRaptorData, backwardRaptorData, stopEventPermutation, forwardShortcutGraph, backwardShortcutGraph, coarsePartition, coarseLevel, finePartition);
            write(routeFlagsData);
        } else {
            NestedDissection finePartition(finePartitionPath);
            RouteFlagsDataType<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition> routeFlagsData(forwardRaptorData, backwardRaptorData, stopEventPermutation, forwardShortcutGraph, backwardShortcutGraph, coarsePartition, coarseLevel, finePartition, fineLevel);
            write(routeFlagsData);
        }
    }

    template<typename ROUTE_FLAGS_DATA>
    inline void write(const ROUTE_FLAGS_DATA& data) const noexcept {
        const std::string outputFilename = getParameter("Output file");
        data.serialize(outputFilename);
    }
};


class RunRouteFlagsPreprocessing : public ParameterizedCommand {

public:
    RunRouteFlagsPreprocessing(BasicShell& shell) :
        ParameterizedCommand(shell, "runRouteFlagsPreprocessing", "Runs preprocessing on the given Route-Flags data.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Threads");
        addParameter("Alternative backward search?", "false");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            chooseAlternativeBackwardSearch<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            chooseAlternativeBackwardSearch<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void chooseAlternativeBackwardSearch() const noexcept {
        const bool useAlternativeBackwardSearch = getParameter<bool>("Alternative backward search?");
        if (useAlternativeBackwardSearch) {
            run<FLAGS, PARTITION, true>();
        } else {
            run<FLAGS, PARTITION, false>();
        }
    }

    template<typename FLAGS, typename PARTITION, bool ALTERNATIVE>
    inline void run() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        const size_t numThreads = getParameter<size_t>("Threads");
        Timer timer;
        ThreadPinning threadPinner(numThreads, 1);
        RAPTOR::RouteFlags::Preprocessing<FLAGS, PARTITION, ALTERNATIVE> preprocessing(routeFlagsData);
        timer.restart();
        preprocessing.run(threadPinner, 0);
        std::cout << "Total time: " << timer.elapsedMilliseconds() << " ms." << std::endl;
        routeFlagsData.serialize(routeFlagsFilename);
    }
};

class VisualizeIntraCellConnections : public ParameterizedCommand {

public:
    VisualizeIntraCellConnections(BasicShell& shell) :
        ParameterizedCommand(shell, "visualizeIntraCellConnections", "Visualizes intra-cell connections in the given Route-Flags data.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            readData<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            readData<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void readData() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        const std::string outputFilename = getParameter("Output file");
        const RAPTOR::Data& raptorData = routeFlagsData.getRaptorData(FORWARD);
        visualizePartition(routeFlagsData.getCoarsePartition(), raptorData, outputFilename + "_coarse");
        visualizePartition(routeFlagsData.getFinePartition(), raptorData, outputFilename + "_fine");
    }

    template<typename PARTITION>
    inline static void visualizePartition(const PARTITION& partition, const RAPTOR::Data& raptorData, std::string filename) noexcept {
        MapVisualization<PDF> visualization(filename, BoundingBoxes::bestMatch(raptorData.transferGraph[Coordinates]));
        for (size_t cell = 0; cell < partition.numberOfCells(); cell++) {
            for (const RouteId route : raptorData.routes()) {
                RAPTOR::TripIterator tripIterator = raptorData.getTripIterator(route);
                while (tripIterator.hasFurtherStops()) {
                    const StopId from = tripIterator.stop();
                    tripIterator.nextStop();
                    const StopId to = tripIterator.stop();
                    const bool isIntraCell = partition.isInCell(from, cell) && partition.isInCell(to, cell);
                    visualization.drawLine(raptorData.stopData[from].coordinates, raptorData.stopData[to].coordinates, getColor(isIntraCell), isIntraCell ? 2 : 1);
                }
            }
            visualization.newPage();
        }
    }
};

class VisualizeTransferShortcuts : public ParameterizedCommand {

public:
    VisualizeTransferShortcuts(BasicShell& shell) :
        ParameterizedCommand(shell, "visualizeTransferShortcuts", "Visualizes transfer shortcuts.") {
        addParameter("RAPTOR network");
        addParameter("Shortcut graph");
        addParameter("Output file");
        addParameter("Walking limit", "-1");
    }

    virtual void execute() noexcept {
        const std::string raptorPath = getParameter("RAPTOR network");
        RAPTOR::Data raptorData = RAPTOR::Data::FromBinary(raptorPath);
        raptorData.useImplicitDepartureBufferTimes();
        const std::string shortcutGraphPath = getParameter("Shortcut graph");
        RAPTOR::TransferGraph shortcutGraph;
        shortcutGraph.readBinary(shortcutGraphPath);
        const std::string outputFilename = getParameter("Output file");
        const int walkingLimit = getParameter<int>("Walking limit");

        MapVisualization<PNG> visualization(outputFilename, BoundingBoxes::bestMatch(raptorData.transferGraph[Coordinates]));
        for (const auto [edge, from] : shortcutGraph.edgesWithFromVertex()) {
            const Vertex to = shortcutGraph.get(ToVertex, edge);
            if (shortcutGraph.get(TravelTime, edge) <= walkingLimit) continue;
            visualization.drawLine(raptorData.stopData[from].coordinates, raptorData.stopData[to].coordinates, highlightColor, 2);
        }
    }
};


class VisualizeRouteFlagsPreprocessing : public ParameterizedCommand {

public:
    VisualizeRouteFlagsPreprocessing(BasicShell& shell) :
        ParameterizedCommand(shell, "visualizeRouteFlagsPreprocessing", "Visualizes flags for the given pair of cells in the given Route-Flags data.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Coarse cell");
        addParameter("Fine cell");
        addParameter("Departure/arrival?", { "departure", "arrival" });
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            readData<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            readData<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void readData() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        const size_t coarseCell = getParameter<size_t>("Coarse cell");
        if (coarseCell >= routeFlagsData.getCoarsePartition().numberOfCells()) return;
        const size_t fineCell = getParameter<size_t>("Fine cell");
        if (fineCell >= routeFlagsData.getFinePartition().numberOfCells()) return;
        const std::string departureOrArrivalString = getParameter("Departure/arrival?");
        bool departureOrArrival = (departureOrArrivalString == "arrival");
        const std::string outputFilename = getParameter("Output file");

        const std::string baseFilename = outputFilename + "_" + String::lexicalCast<std::string, int>(coarseCell) + "_" + String::lexicalCast<std::string, int>(fineCell) + "_" + departureOrArrivalString;
        draw(routeFlagsData, coarseCell, fineCell, FORWARD, departureOrArrival, baseFilename + "_forward");
        draw(routeFlagsData, coarseCell, fineCell, BACKWARD, departureOrArrival, baseFilename + "_backward");
    }

    template<typename ROUTE_FLAGS_DATA>
    inline static void draw(const ROUTE_FLAGS_DATA& flagsData, const size_t coarseCell, const size_t fineCell, const bool direction, const bool departureOrArrival, const std::string filename) noexcept {
        const RAPTOR::Data& raptorData = flagsData.getRaptorData(direction);
        MapVisualization<PDF> visualization(filename, BoundingBoxes::bestMatch(raptorData.transferGraph[Coordinates]));
        drawRoutes(visualization, raptorData, flagsData, coarseCell, fineCell, direction, departureOrArrival, true);
        visualization.newPage();
        drawRoutes(visualization, raptorData, flagsData, coarseCell, fineCell, direction, departureOrArrival, false);
        visualization.newPage();
        drawGraph(visualization, raptorData, flagsData.getShortcuts(direction).coarseCellShortcuts[coarseCell]);
        drawGraph(visualization, raptorData, flagsData.getShortcuts(direction).fineCellShortcuts[coarseCell][fineCell]);
    }

    template<typename VISUALIZATION, typename ROUTE_FLAGS_DATA>
    inline static void drawRoutes(VISUALIZATION& visualization, const RAPTOR::Data& raptorData, const ROUTE_FLAGS_DATA& flagsData, const size_t coarseCell, const size_t fineCell, const bool direction, const bool departureOrArrival, const bool drawUnflagged) noexcept {
        for (const RouteId route : raptorData.routes()) {
            RAPTOR::TripIterator tripIterator = raptorData.getTripIterator(route);
            while (tripIterator.hasFurtherStops()) {
                const StopId from = tripIterator.stop();
                tripIterator.nextStop();
                const bool hasNecessary = flagsData.isNecessaryRouteSegment(direction, departureOrArrival, route, tripIterator.getStopIndex(), coarseCell, fineCell);
                const StopId to = tripIterator.stop();
                if (drawUnflagged || hasNecessary) {
                    visualization.drawLine(raptorData.stopData[from].coordinates, raptorData.stopData[to].coordinates, getColor(hasNecessary), hasNecessary ? 2 : 1);
                }
            }
        }
    }

    template<typename VISUALIZATION>
    inline static void drawGraph(VISUALIZATION& visualization, const RAPTOR::Data& raptorData, const RAPTOR::TransferGraph& graph) noexcept {
        for (const auto [edge, from] : graph.edgesWithFromVertex()) {
            const Vertex to = graph.get(ToVertex, edge);
            visualization.drawLine(raptorData.stopData[from].coordinates, raptorData.stopData[to].coordinates, highlightColor, 2);
        }
    }
};

class RunRouteFlagsQuery : public ParameterizedCommand {

public:
    RunRouteFlagsQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "runRouteFlagsQuery", "Runs a single Route-Flags query.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("CH data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Backward?");
        addParameter("Departure/arrival time");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            chooseDirection<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            chooseDirection<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void chooseDirection() const noexcept {
        const bool direction = getParameter<bool>("Backward?");
        if (direction == BACKWARD) {
            run<BACKWARD, FLAGS, PARTITION>();
        } else {
            run<FORWARD, FLAGS, PARTITION>();
        }
    }

    template<bool DIRECTION, typename FLAGS_TYPE, typename FINE_PARTITION_TYPE>
    inline void run() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS_TYPE, FINE_PARTITION_TYPE> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const std::string chFilename = getParameter("CH data");
        CH::CH chData(chFilename);

        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const int departureTime = getParameter<int>("Departure/arrival time");
        std::cout << "Running " << (DIRECTION ? "backward" : "forward") << " query from " << source << " to " << target << std::endl;
        std::cout << (DIRECTION ? "Arrival" : "Departure") << " time: " << departureTime << std::endl;

        const RAPTOR::Data& forwardRaptorData = routeFlagsData.getRaptorData(DIRECTION);
        const RAPTOR::Data& backwardRaptorData = routeFlagsData.getRaptorData(!DIRECTION);
        RAPTOR::RouteFlags::Query<FLAGS_TYPE, FINE_PARTITION_TYPE, true, RAPTOR::SimpleProfiler, DIRECTION> routeFlagsQuery(routeFlagsData, chData);
        RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::SimpleProfiler, true, false, false> raptorQuery(forwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);
        routeFlagsQuery.run(source, departureTime, target);
        raptorQuery.run(source, departureTime, target);
        std::vector<RAPTOR::ArrivalLabel> routeFlagsArrivals = routeFlagsQuery.getArrivals();
        std::vector<RAPTOR::ArrivalLabel> raptorArrivals = raptorQuery.getArrivals();
        if (!Vector::equals(routeFlagsArrivals, raptorArrivals)) {
            std::cout << "Arrivals not identical!" << std::endl;
        }
        std::cout << "Route-Flags journeys:" << std::endl;
        std::cout << routeFlagsQuery.getJourneys();
        std::cout << "RAPTOR journeys:" << std::endl;
        std::cout << raptorQuery.getJourneys();
    }

    template<typename ROUTE_FLAGS_DATA>
    inline static void printRoute(const RAPTOR::Data& raptorData, const ROUTE_FLAGS_DATA& routeFlagsData, const RouteId route) noexcept {
        for (const StopId stop : raptorData.stopsOfRoute(route)) {
            printCells(routeFlagsData, stop);
        }
        std::cout << std::endl;
    }

    template<typename ROUTE_FLAGS_DATA>
    inline static void printFootpath(Dijkstra<RAPTOR::TransferGraph>& dijkstra, const ROUTE_FLAGS_DATA& routeFlagsData, const Vertex from, const Vertex to) noexcept {
        dijkstra.run(from, to);
        for (const Vertex v : dijkstra.getPath(to)) {
            printCells(routeFlagsData, v);
        }
        std::cout << std::endl;
    }

    template<typename ROUTE_FLAGS_DATA>
    inline static void printCells(const ROUTE_FLAGS_DATA& routeFlagsData, const Vertex v) noexcept {
        std::cout << v << ": " << std::endl;
        Vector::printConcise(routeFlagsData.getCoarsePartition().cellsOfVertex(v));
        std::cout << std::endl;
        Vector::printConcise(routeFlagsData.getFinePartition().cellsOfVertex(v));
        std::cout << std::endl;
    }
};

/*class VisualizeRouteFlagsQuery : public ParameterizedCommand {

public:
    VisualizeRouteFlagsQuery(BasicShell& shell) :
        ParameterizedCommand(shell, "visualizeRouteFlagsQuery", "Visualizes the search space for the given Route-Flags query.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("CH data");
        addParameter("Source (-1 for random)");
        addParameter("Target (-1 for random)");
        addParameter("Backward?");
        addParameter("Departure/arrival time");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            readData<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            readData<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void readData() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const std::string chFilename = getParameter("CH data");
        CH::CH chData(chFilename);

        srand(42);

        int sourceId = getParameter<int>("Source (-1 for random)");
        if (sourceId == -1) sourceId = rand() % routeFlagsData.getRaptorData(FORWARD).transferGraph.numVertices();
        const Vertex source(sourceId);
        int targetId = getParameter<int>("Target (-1 for random)");
        if (targetId == -1) targetId = rand() % routeFlagsData.getRaptorData(FORWARD).transferGraph.numVertices();
        const Vertex target(targetId);
        const bool direction = getParameter<bool>("Backward?");
        const int departureTime = getParameter<int>("Departure/arrival time");

        const std::string outputFilename = getParameter("Output file");
        MapVisualization<PDF> visualization(outputFilename, BoundingBoxes::bestMatch(routeFlagsData.getRaptorData(FORWARD).transferGraph[Coordinates]));

        std::cout << "Running " << (direction ? "backward" : "forward") << " query from " << source << " to " << target << std::endl;
        std::cout << (direction ? "Arrival" : "Departure") << " time: " << departureTime << std::endl;
        std::cout << "Source coarse cells: ";
        Vector::printConcise(routeFlagsData.getCoarsePartition().cellsOfVertex(source));
        std::cout << std::endl;
        std::cout << "Source fine cells: ";
        Vector::printConcise(routeFlagsData.getFinePartition().cellsOfVertex(source));
        std::cout << std::endl;
        std::cout << "Target coarse cells: ";
        Vector::printConcise(routeFlagsData.getCoarsePartition().cellsOfVertex(target));
        std::cout << std::endl;
        std::cout << "Target fine cells: ";
        Vector::printConcise(routeFlagsData.getFinePartition().cellsOfVertex(target));
        std::cout << std::endl;

        if (direction == BACKWARD) {
            visualize<BACKWARD>(visualization, routeFlagsData, chData, source, target, departureTime);
        } else {
            visualize<FORWARD>(visualization, routeFlagsData, chData, source, target, departureTime);
        }

    }

    template<bool DIRECTION, typename VISUALIZATION, typename FLAGS_TYPE, typename FINE_PARTITION_TYPE>
    inline static void visualize(VISUALIZATION& visualization, const RouteFlagsDataType<FLAGS_TYPE, FINE_PARTITION_TYPE>& routeFlagsData, const CH::CH& chData, const Vertex source, const Vertex target, const int departureTime) noexcept {
        const RAPTOR::Data& forwardRaptorData = routeFlagsData.getRaptorData(DIRECTION);
        const RAPTOR::Data& backwardRaptorData = routeFlagsData.getRaptorData(!DIRECTION);
        RAPTOR::RouteFlags::Query<FLAGS_TYPE, FINE_PARTITION_TYPE, true, RAPTOR::SearchSpaceProfiler, DIRECTION> routeFlagsQuery(routeFlagsData, chData);
        RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::SearchSpaceProfiler, true, false, false> raptorQuery(forwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);

        routeFlagsQuery.run(source, departureTime, target);
        raptorQuery.run(source, departureTime, target);

        const std::vector<int> routeFlagsShortcutScans = routeFlagsQuery.getProfiler().getScansPerShortcut();
        const std::vector<int> raptorEdgeScans = raptorQuery.getProfiler().getScansPerEdge();
        const std::vector<bool> routeFlagsScannedRoutes = routeFlagsQuery.getProfiler().getScannedRoutes();
        const std::vector<bool> raptorScannedRoutes = raptorQuery.getProfiler().getScannedRoutes();
        const std::vector<int> routeFlagsRouteSegmentScans = routeFlagsQuery.getProfiler().getScansPerRouteSegment();
        const std::vector<int> raptorRouteSegmentScans = raptorQuery.getProfiler().getScansPerRouteSegment();

        const int maxEdgeScans = std::max(Vector::max(routeFlagsShortcutScans), Vector::max(raptorEdgeScans));
        const int maxRouteSegmentScans = std::max(Vector::max(routeFlagsRouteSegmentScans), Vector::max(raptorRouteSegmentScans));

        drawRoutes(visualization, forwardRaptorData, routeFlagsScannedRoutes, routeFlagsRouteSegmentScans, maxRouteSegmentScans);
        drawSourceAndTarget(visualization, forwardRaptorData, source, target);
        visualization.newPage();
        drawRoutes(visualization, forwardRaptorData, raptorScannedRoutes, raptorRouteSegmentScans, maxRouteSegmentScans);
        drawSourceAndTarget(visualization, forwardRaptorData, source, target);
        visualization.newPage();
        drawFootpaths(visualization, routeFlagsQuery.getShortcuts(), routeFlagsShortcutScans, maxEdgeScans);
        drawSourceAndTarget(visualization, forwardRaptorData, source, target);
        visualization.newPage();
        drawFootpaths(visualization, forwardRaptorData.transferGraph, raptorEdgeScans, maxEdgeScans);
        drawSourceAndTarget(visualization, forwardRaptorData, source, target);
    }

    template<typename VISUALIZATION>
    inline static void drawSourceAndTarget(VISUALIZATION& visualization, const RAPTOR::Data& raptorData, const Vertex source, const Vertex target) noexcept {
        visualization.drawPoint(raptorData.stopData[source].coordinates, Color::Black, 100);
        visualization.write("S", raptorData.stopData[source].coordinates, Color::Black, XAlignment::Center, YAlignment::Center);
        visualization.drawPoint(raptorData.stopData[target].coordinates, Color::Black, 100);
        visualization.write("T", raptorData.stopData[target].coordinates, Color::Black, XAlignment::Center, YAlignment::Center);

    }

    template<typename VISUALIZATION>
    inline static void drawRoutes(VISUALIZATION& visualization, const RAPTOR::Data& raptorData, const std::vector<bool>& scannedRoutes, const std::vector<int>& routeSegmentScans, const int maxRouteSegmentScans) noexcept {
        for (const RouteId route : raptorData.routes()) {
            if (!scannedRoutes[route]) continue;
            RAPTOR::TripIterator tripIterator = raptorData.getTripIterator(route);
            bool routeEntered = false;
            while (tripIterator.hasFurtherStops()) {
                const StopId from = tripIterator.stop();
                tripIterator.nextStop();
                const StopId to = tripIterator.stop();
                const int scans = routeSegmentScans[raptorData.getRouteSegmentNum(route, tripIterator.getStopIndex())];
                routeEntered |= scans > 0;
                if (routeEntered) {
                    visualization.drawLine(raptorData.stopData[from].coordinates, raptorData.stopData[to].coordinates, getGradientColor(maxRouteSegmentScans, scans), 1);
                }
            }
        }
    }

    template<typename VISUALIZATION>
    inline static void drawFootpaths(VISUALIZATION& visualization, const RAPTOR::TransferGraph& graph, const std::vector<int>& edgeScans, const int maxEdgeScans) noexcept {
        for (const auto [edge, from] : graph.edgesWithFromVertex()) {
            if (edgeScans[edge] == 0) continue;
            const Vertex to = graph.get(ToVertex, edge);
            visualization.drawLine(graph.get(Coordinates, from), graph.get(Coordinates, to), getGradientColor(maxEdgeScans, edgeScans[edge]), 1);
        }
    }

    inline static Color getGradientColor(const int maximum, const int value) noexcept {
        AssertMsg(value <= maximum, "Value is higher than maximum!");
        if (value == 0) return Color::LightGrey;
        return Color::getGradientColor(Color::Green, Color::Red, value / (double) maximum);
    }
};*/

class RunRouteFlagsQueries : public ParameterizedCommand {

public:
    RunRouteFlagsQueries(BasicShell& shell) :
        ParameterizedCommand(shell, "runRouteFlagsQueries", "Runs random Route-Flags queries and writes results to a .csv file.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("CH data");
        addParameter("Backward?");
        addParameter("Output file");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            chooseDirection<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            chooseDirection<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void chooseDirection() const noexcept {
        const bool direction = getParameter<bool>("Backward?");
        if (direction == BACKWARD) {
            run<BACKWARD, FLAGS, PARTITION>();
        } else {
            run<FORWARD, FLAGS, PARTITION>();
        }
    }

    template<bool DIRECTION, typename FLAGS_TYPE, typename FINE_PARTITION_TYPE>
    inline void run() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS_TYPE, FINE_PARTITION_TYPE> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const std::string chFilename = getParameter("CH data");
        CH::CH chData(chFilename);

        const std::string outputFilename = getParameter("Output file");
        IO::OFStream result(outputFilename + ".csv");
        result << "queryId,sourceStop,targetStop,departureTime,routeFlagsTime,raptorTime,speedup\n";
        const size_t n = getParameter<size_t>("Number of queries");
        srand(42);

        const RAPTOR::Data& forwardRaptorData = routeFlagsData.getRaptorData(DIRECTION);
        const RAPTOR::Data& backwardRaptorData = routeFlagsData.getRaptorData(!DIRECTION);
        RAPTOR::RouteFlags::Query<FLAGS_TYPE, FINE_PARTITION_TYPE, true, RAPTOR::NoProfiler, DIRECTION> routeFlagsQuery(routeFlagsData, chData);
        RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::NoProfiler, true, false, false> raptorQuery(forwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);
        Progress progress(n);
        progress.SetCheckTimeStep(10);
        Timer timer;
        double totalRouteFlagsTime = 0;
        double totalRaptorTime = 0;
        for (size_t i = 0; i < n; i++) {
            const StopId source = StopId(rand() % forwardRaptorData.numberOfStops());
            const StopId target = StopId(rand() % forwardRaptorData.numberOfStops());
            const int departureTime = (rand() % (16 * 60 * 60)) + (5 * 60 * 60);
            timer.restart();
            routeFlagsQuery.run(source, departureTime, target);
            const double routeFlagsTime = timer.elapsedMilliseconds();
            totalRouteFlagsTime += routeFlagsTime;
            timer.restart();
            raptorQuery.run(source, departureTime, target);
            const double raptorTime = timer.elapsedMilliseconds();
            totalRaptorTime += raptorTime;
            result << i << "," << int(source) << "," << int(target) << "," << departureTime << "," << routeFlagsTime << "," << raptorTime << "," << raptorTime / routeFlagsTime << "\n";
            std::vector<RAPTOR::ArrivalLabel> routeFlagsArrivals = routeFlagsQuery.getArrivals();
            std::vector<RAPTOR::ArrivalLabel> raptorArrivals = raptorQuery.getArrivals();
            if (!Vector::equals(routeFlagsArrivals, raptorArrivals)) {
                std::cout << "Arrivals not identical for query " << i << ": " << source << " -> " << target << " @ " << departureTime << "!" << std::endl;
                std::cout << "Route-Flags journeys:" << std::endl;
                std::cout << routeFlagsQuery.getJourneys();
                std::cout << "RAPTOR journeys:" << std::endl;
                std::cout << raptorQuery.getJourneys();
                break;
            }
            progress++;
        }
        result << "-1,-,-,-," << totalRouteFlagsTime << "," << totalRaptorTime << "," << totalRaptorTime / totalRouteFlagsTime << "\n";
        result.flush();
    }
};

class CheckFlags : public ParameterizedCommand {

public:
    CheckFlags(BasicShell& shell) :
        ParameterizedCommand(shell, "checkFlags", "Outputs which flags are set for the given stop event.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Source");
        addParameter("Target");
        addParameter("Route");
        addParameter("Stop");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            checkFlags<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            checkFlags<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void checkFlags() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const RouteId route = getParameter<RouteId>("Route");
        const StopId stop = getParameter<StopId>("Stop");
        routeFlagsData.printFlags(source, target, route, stop);
    }
};

class CheckShortcut : public ParameterizedCommand {

public:
    CheckShortcut(BasicShell& shell) :
        ParameterizedCommand(shell, "checkShortcut", "Checks if a shortcut between the given stops exists for the given source and target cells.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Source");
        addParameter("Target");
        addParameter("From");
        addParameter("To");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            checkShortcut<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            checkShortcut<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void checkShortcut() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        const StopId from = getParameter<StopId>("From");
        const StopId to = getParameter<StopId>("To");
        routeFlagsData.printShortcuts(source, target, from, to);
    }
};

class CellsOfRoute : public ParameterizedCommand {

public:
    CellsOfRoute(BasicShell& shell) :
        ParameterizedCommand(shell, "cellsOfRoute", "Outputs the cells for each stop on the given route.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Route");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            printCells<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            printCells<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void printCells() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const RouteId route = getParameter<RouteId>("Route");
        for (const StopId stop : routeFlagsData.getRaptorData(FORWARD).stopsOfRoute(route)) {
            routeFlagsData.printVertexCells(stop);
        }
    }
};

class CellsOfFootpath : public ParameterizedCommand {

public:
    CellsOfFootpath(BasicShell& shell) :
        ParameterizedCommand(shell, "cellsOfFootpath", "Outputs the cells for each vertex on the shortest path between source and target.") {
        addParameter("Flags type", { "event", "trip", "route" });
        addParameter("Fine partition type", { "vp", "nd" });
        addParameter("Route-Flags data");
        addParameter("Source");
        addParameter("Target");
    }

    virtual void execute() noexcept {
        const std::string flagsType = getParameter("Flags type");
        if (flagsType == "event") {
            chooseFinePartition<RAPTOR::RouteFlags::StopEventFlags>();
        } else if (flagsType == "trip") {
            chooseFinePartition<RAPTOR::RouteFlags::TripFlags>();
        } else {
            chooseFinePartition<RAPTOR::RouteFlags::RouteFlags>();
        }
    }

private:
    template<typename FLAGS>
    inline void chooseFinePartition() const noexcept {
        const std::string finePartitionType = getParameter("Fine partition type");
        if (finePartitionType == "vp") {
            printCells<FLAGS, RAPTOR::RouteFlags::EdgeSeparatorPartition>();
        } else {
            printCells<FLAGS, RAPTOR::RouteFlags::VertexSeparatorPartition>();
        }
    }

    template<typename FLAGS, typename PARTITION>
    inline void printCells() const noexcept {
        const std::string routeFlagsFilename = getParameter("Route-Flags data");
        RouteFlagsDataType<FLAGS, PARTITION> routeFlagsData(routeFlagsFilename);
        routeFlagsData.sortShortcuts();
        const Vertex source = getParameter<Vertex>("Source");
        const Vertex target = getParameter<Vertex>("Target");
        Dijkstra<RAPTOR::TransferGraph> dijkstra(routeFlagsData.getRaptorData(FORWARD).transferGraph);
        dijkstra.run(source, target);
        for (const Vertex vertex : dijkstra.getPath(target)) {
            routeFlagsData.printVertexCells(vertex);
        }
    }
};
