#pragma once

#include "../UnitTests.h"

#include "../../Algorithms/CH/Preprocessing/CHBuilder.h"
#include "../../Algorithms/CH/Preprocessing/MinLevelKey.h"
#include "../../Algorithms/CH/Preprocessing/BlockingKey.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"

#include "../../Algorithms/RAPTOR/ULTRA/Builder.h"
#include "../../Algorithms/RAPTOR/ULTRARAPTOR.h"

#include "../../Algorithms/TripBased/Preprocessing/ULTRABuilder.h"
#include "../../Algorithms/TripBased/Query/Query.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

namespace UnitTests {

class ULTRATests {

public:
    inline void checkWeakDominationStopToStop() {
        const RAPTOR::Data data = buildWeakDominationNetworkRAPTOR();
        RAPTOR::ULTRA::Builder<false> shortcutBuilder(data);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        RAPTOR::Data shortcutData = data;
        Graph::move(std::move(shortcutBuilder.getShortcutGraph()), shortcutData.transferGraph);

        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(1), StopId(2)), "Shortcut A -> B is missing!");
        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(3), StopId(5)), "Shortcut C -> E is missing!");
        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(4), StopId(5)), "Shortcut D -> E is missing!");
        UnitTests::check(shortcutData.transferGraph.numEdges() == 3, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        RAPTOR::ULTRARAPTOR<> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(6));
        UnitTests::check(query.getEarliestArrivalTime() == 10800, "Query answered incorrectly!");
    }

    inline void checkWeakDominationEventToEvent() {
        const RAPTOR::Data data = buildWeakDominationNetworkRAPTOR();
        TripBased::Data shortcutData(data);
        TripBased::ULTRABuilder<false> shortcutBuilder(shortcutData);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        Graph::move(std::move(shortcutBuilder.getStopEventGraph()), shortcutData.stopEventGraph);

        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(1), Vertex(2)), "Shortcut R1 -> R2 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(1), Vertex(4)), "Shortcut R1 -> R3 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(3), Vertex(6)), "Shortcut R2 -> R5 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(9), Vertex(6)), "Shortcut R4 -> R5 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.numEdges() == 4, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        TripBased::Query<TripBased::NoProfiler> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(6));
        UnitTests::check(query.getEarliestArrivalTime() == 10800, "Query answered incorrectly!");
    }

    inline void checkCyclicalDominationStopToStop() {
        const RAPTOR::Data data = buildCyclicalDominationNetworkRAPTOR();
        RAPTOR::ULTRA::Builder<false> shortcutBuilder(data);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        RAPTOR::Data shortcutData = data;
        Graph::move(std::move(shortcutBuilder.getShortcutGraph()), shortcutData.transferGraph);

        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(1), StopId(2)), "Shortcut A -> B is missing!");
        UnitTests::check(shortcutData.transferGraph.numEdges() == 1, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        RAPTOR::ULTRARAPTOR<> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(3));
        UnitTests::check(query.getEarliestArrivalTime() == 6000, "Query answered incorrectly!");
    }

    inline void checkCyclicalDominationEventToEvent() {
        const RAPTOR::Data data = buildCyclicalDominationNetworkRAPTOR();
        TripBased::Data shortcutData(data);
        TripBased::ULTRABuilder<false> shortcutBuilder(shortcutData);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        Graph::move(std::move(shortcutBuilder.getStopEventGraph()), shortcutData.stopEventGraph);

        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(1), Vertex(2)), "Shortcut R1 -> R2 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.numEdges() == 1, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        TripBased::Query<TripBased::NoProfiler> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(3));
        UnitTests::check(query.getEarliestArrivalTime() == 6000, "Query answered incorrectly!");
    }

    inline void checkRouteOrderStopToStop() {
        const RAPTOR::Data data = buildRouteOrderNetworkRAPTOR();
        RAPTOR::ULTRA::Builder<false> shortcutBuilder(data);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        RAPTOR::Data shortcutData = data;
        Graph::move(std::move(shortcutBuilder.getShortcutGraph()), shortcutData.transferGraph);

        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(1), StopId(2)), "Shortcut A -> B is missing!");
        UnitTests::check(shortcutData.transferGraph.hasEdge(StopId(4), StopId(5)), "Shortcut D -> E is missing!");
        UnitTests::check(shortcutData.transferGraph.numEdges() == 2, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        RAPTOR::ULTRARAPTOR<> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(6));
        UnitTests::check(query.getEarliestArrivalTime() == 10800, "Query answered incorrectly!");
    }

    inline void checkRouteOrderEventToEvent() {
        const RAPTOR::Data data = buildRouteOrderNetworkRAPTOR();
        TripBased::Data shortcutData(data);
        TripBased::ULTRABuilder<false> shortcutBuilder(shortcutData);
        shortcutBuilder.computeShortcuts(ThreadPinning(1, 1), INFTY);
        Graph::move(std::move(shortcutBuilder.getStopEventGraph()), shortcutData.stopEventGraph);


        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(1), Vertex(2)), "Shortcut R1 -> R2 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.hasEdge(Vertex(3), Vertex(12)), "Shortcut R2 -> R4 is missing!");
        UnitTests::check(shortcutData.stopEventGraph.numEdges() == 2, "Number of shortcuts is incorrect!");

        CH::CH chData = buildCH(data);
        TripBased::Query<TripBased::NoProfiler> query(shortcutData, chData);
        query.run(StopId(0), 0, StopId(6));
        UnitTests::check(query.getEarliestArrivalTime() == 10800, "Query answered incorrectly!");
    }

protected:
    inline RAPTOR::Data buildWeakDominationNetworkRAPTOR() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(9);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.0));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.0));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 2.0, 0.0));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 3.0, 0.0));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 4.0, 0.0));
        transferGraph.set(Coordinates, Vertex(5), Geometry::Point(Construct::XY, 5.0, 0.0));
        transferGraph.set(Coordinates, Vertex(6), Geometry::Point(Construct::XY, 6.0, 0.0));
        transferGraph.set(Coordinates, Vertex(7), Geometry::Point(Construct::XY, 2.0, 2.0));
        transferGraph.set(Coordinates, Vertex(8), Geometry::Point(Construct::XY, 3.0, 1.0));

        transferGraph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 600);
        transferGraph.addEdge(Vertex(2), Vertex(8)).set(TravelTime, 600);
        transferGraph.addEdge(Vertex(3), Vertex(4)).set(TravelTime, 600);
        transferGraph.addEdge(Vertex(4), Vertex(5)).set(TravelTime, 600);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  0);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(4)),  0);
        stops.emplace_back("E", transferGraph.get(Coordinates, Vertex(5)),  0);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(6)),  0);
        stops.emplace_back("X", transferGraph.get(Coordinates, Vertex(7)),  0);
        stops.emplace_back("BB", transferGraph.get(Coordinates, Vertex(8)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("B -> C", "R2", 1);
        trips.emplace_back("B -> X", "R3", 1);
        trips.emplace_back("BB -> D", "R4", 1);
        trips.emplace_back("E -> T", "R5", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 2400, TripId(0));
        connections.emplace_back(StopId(2), StopId(3), 3600, 6000, TripId(1));
        connections.emplace_back(StopId(2), StopId(7), 4200, 7200, TripId(2));
        connections.emplace_back(StopId(8), StopId(4), 4800, 7200, TripId(3));
        connections.emplace_back(StopId(5), StopId(6), 8400, 10800, TripId(4));

        CSA::Data csa = CSA::Data::FromInput(stops, connections, trips, transferGraph);
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        RAPTOR::Data raptor = RAPTOR::Data::FromIntermediate(inter);
        raptor.useImplicitDepartureBufferTimes();
        return raptor;
    }

    inline RAPTOR::Data buildCyclicalDominationNetworkRAPTOR() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(5);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.0));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.0));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 2.0, 0.0));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 3.0, 0.0));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 0.0, 1.0));

        transferGraph.addEdge(Vertex(0), Vertex(4)).set(TravelTime, 0);
        transferGraph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 600);
        transferGraph.addEdge(Vertex(4), Vertex(0)).set(TravelTime, 0);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  0);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("S2", transferGraph.get(Coordinates, Vertex(4)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("B -> T", "R2", 1);
        trips.emplace_back("S2 -> A", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 3000, TripId(0));
        connections.emplace_back(StopId(2), StopId(3), 3600, 6000, TripId(1));
        connections.emplace_back(StopId(4), StopId(1), 0, 3000, TripId(2));

        CSA::Data csa = CSA::Data::FromInput(stops, connections, trips, transferGraph);
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        RAPTOR::Data raptor = RAPTOR::Data::FromIntermediate(inter);
        raptor.useImplicitDepartureBufferTimes();
        return raptor;
    }

    inline RAPTOR::Data buildRouteOrderNetworkRAPTOR() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(7);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.0));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.0));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 2.0, 0.0));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 2.0, -1.0));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 3.0, 0.0));
        transferGraph.set(Coordinates, Vertex(5), Geometry::Point(Construct::XY, 4.0, 0.0));
        transferGraph.set(Coordinates, Vertex(6), Geometry::Point(Construct::XY, 5.0, 0.0));

        transferGraph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 600);
        transferGraph.addEdge(Vertex(1), Vertex(3)).set(TravelTime, 300);
        transferGraph.addEdge(Vertex(4), Vertex(5)).set(TravelTime, 0);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  0);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(4)),  0);
        stops.emplace_back("E", transferGraph.get(Coordinates, Vertex(5)),  0);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(6)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("B -> D", "R2", 1);
        trips.emplace_back("C -> D", "R3", 1);
        trips.emplace_back("C -> D", "R3", 1);
        trips.emplace_back("E -> T", "R4", 1);
        trips.emplace_back("B -> T", "R5", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 3000, TripId(0));
        connections.emplace_back(StopId(2), StopId(4), 3600, 7200, TripId(1));
        connections.emplace_back(StopId(3), StopId(2), 2400, 3600, TripId(2));
        connections.emplace_back(StopId(2), StopId(4), 3600, 7200, TripId(2));
        connections.emplace_back(StopId(3), StopId(2), 6000, 7200, TripId(3));
        connections.emplace_back(StopId(2), StopId(4), 7200, 10800, TripId(3));
        connections.emplace_back(StopId(5), StopId(6), 7200, 10800, TripId(4));
        connections.emplace_back(StopId(3), StopId(6), 2400, 7200, TripId(5));

        CSA::Data csa = CSA::Data::FromInput(stops, connections, trips, transferGraph);
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        RAPTOR::Data raptor = RAPTOR::Data::FromIntermediate(inter);
        raptor.useImplicitDepartureBufferTimes();
        return raptor;
    }

    inline CH::CH buildCH(const RAPTOR::Data& data) const noexcept {
        TravelTimeGraph graph;
        Graph::copy(data.transferGraph, graph);

        using Profiler = CH::NoProfiler;
        using WitnessSearch = CH::BidirectionalWitnessSearch<CHCoreGraph, Profiler, 200>;
        using KeyFunction = CH::GreedyKey<WitnessSearch>;
        using StopCriterion = CH::NoStopCriterion;
        using CHBuilder = CH::Builder<Profiler, WitnessSearch, KeyFunction, StopCriterion, false, false>;

        KeyFunction keyFunction(1024, 1024, 0);
        CHBuilder chBuilder(std::move(graph), graph[TravelTime], keyFunction, StopCriterion());
        chBuilder.run();
        chBuilder.copyCoreToCH();
        return CH::CH(std::move(chBuilder));
    }
};

}
