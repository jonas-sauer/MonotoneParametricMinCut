#pragma once

#include "../../Algorithms/Dijkstra/Dijkstra.h"
#include "../../Algorithms/RAPTOR/RouteFlags/Preprocessing/RangeSearch.h"
#include "../UnitTests.h"

#include "../../DataStructures/Partition/NestedDissection.h"
#include "../../DataStructures/Partition/VertexPartition.h"
#include "../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

namespace UnitTests {

class RouteFlagsPreprocessing {

public:
    using FlagsType = RAPTOR::RouteFlags::StopEventFlags;
    using PartitionType = RAPTOR::RouteFlags::EdgeSeparatorPartition;
    using RouteFlagsDataType = RAPTOR::RouteFlags::RouteFlagsData<FlagsType, PartitionType>;
    using ForwardPreprocessingRaptorType = RAPTOR::RouteFlags::RangeSearch<FlagsType, PartitionType, FORWARD, false, false>;
    using BackwardPreprocessingRaptorType = RAPTOR::RouteFlags::RangeSearch<FlagsType, PartitionType, BACKWARD, false, false>;

    inline void check() {
        checkPreprocessing1();
        checkPreprocessing2();
        checkPreprocessing3();
        checkPreprocessing4();
        checkPreprocessing5();
        checkPreprocessing7();
        checkPreprocessing8();
    }

protected:
    inline void checkPreprocessing1() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA1(), buildNestedDissection1(), 1, buildVertexPartition1());

        ForwardPreprocessingRaptorType forwardPreprocessing(flagsData);
        forwardPreprocessing.run(Vertex(1), 0, 10);
        checkFlags(1, flagsData, BACKWARD, DEPARTURE, 0, 1, { false, false, false, true,  false, false, false });
        checkFlags(1, flagsData, BACKWARD, ARRIVAL,   0, 1, { false, false, true,  false, true,  false, false });

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(3), -10, 0);
        checkFlags(1, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, false, true,  false, false, false });
        checkFlags(1, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  true,  false, true,  false, false });
    }

    inline void checkPreprocessing2() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA2(), buildNestedDissection2(), 1, buildVertexPartition2());

        ForwardPreprocessingRaptorType forwardPreprocessing(flagsData);
        forwardPreprocessing.run(Vertex(2), 0, 12);
        checkFlags(2, flagsData, BACKWARD, DEPARTURE, 0, 1, { false, false, false, false, false, false, true,  false });
        checkFlags(2, flagsData, BACKWARD, ARRIVAL,   0, 1, { false, false, true,  false, false, true,  false, true  });

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(3), -12, 0);
        checkFlags(2, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, false, true,  false, false, true,  false });
        checkFlags(2, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  true,  false, false, true,  false, true  });
    }

    inline void checkPreprocessing3() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA3(), buildNestedDissection3(), 1, buildVertexPartition3());

        ForwardPreprocessingRaptorType forwardPreprocessing(flagsData);
        forwardPreprocessing.run(Vertex(1), 99, 99);
        checkFlags(3, flagsData, BACKWARD, DEPARTURE, 0, 1, { false, false, false, false, false, false, true,  false, false, false });
        checkFlags(3, flagsData, BACKWARD, ARRIVAL,   0, 1, { false, false, true,  false, false, true,  false, true,  false, false });

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(3), -111, -108);
        checkFlags(3, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, false, true,  false, false, true,  false, true,  false });
        checkFlags(3, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  false, false, false, true,  false, true,  false, true  });
    }

    inline void checkPreprocessing4() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA4(), buildNestedDissection4(), 1, buildVertexPartition4());

        ForwardPreprocessingRaptorType forwardPreprocessing(flagsData);
        forwardPreprocessing.run(Vertex(1), 99, 99);
        checkFlags(4, flagsData, BACKWARD, DEPARTURE, 0, 1, { false, false, false, true,  false });
        checkFlags(4, flagsData, BACKWARD, ARRIVAL,   0, 1, { false, false, true,  false, true  });

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(2), -120, -115);
        checkFlags(4, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, false, true,  false });
        checkFlags(4, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  true,  false, true  });
    }

    inline void checkPreprocessing5() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA5(), buildNestedDissection5(), 2, buildVertexPartition5());

        ForwardPreprocessingRaptorType forwardPreprocessing(flagsData);
        forwardPreprocessing.run(Vertex(0), 0, 30);
        checkFlags(5, flagsData, BACKWARD, DEPARTURE, 1, 4, { true,  false, false, false, true,  false, false });
        checkFlags(5, flagsData, BACKWARD, ARRIVAL,   1, 4, { false, true,  false, false, false, false, true  });

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(2), -30, 0);
        checkFlags(5, flagsData, FORWARD,  DEPARTURE, 3, 0, { true,  false, true,  false, true,  false, false });
        checkFlags(5, flagsData, FORWARD,  ARRIVAL,   3, 0, { false, true,  false, true,  false, false, false });

        flagsData.markIntraCellStopEvents();
        checkFlags(6, flagsData, FORWARD,  DEPARTURE, 3, 0, { true,  false, true,  false, true,  true,  false });
        checkFlags(6, flagsData, FORWARD,  ARRIVAL,   3, 0, { false, true,  false, true,  false, false, true  });
        checkFlags(6, flagsData, BACKWARD, DEPARTURE, 1, 4, { true,  false, false, false, true,  false, false });
        checkFlags(6, flagsData, BACKWARD, ARRIVAL,   1, 4, { false, true,  false, false, false, false, true  });
    }

    inline void checkPreprocessing7() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA7(), buildNestedDissection7(), 1, buildVertexPartition7());

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(4), -30, -30);
        checkFlags(7, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, true,  false, false, false, false, false });
        checkFlags(7, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  false, false, true,  false, false, false });
    }

    inline void checkPreprocessing8() {
        RouteFlagsDataType flagsData = buildRouteFlagsData(buildNetworkCSA8(), buildNestedDissection8(), 1, buildVertexPartition8());

        BackwardPreprocessingRaptorType backwardPreprocessing(flagsData);
        backwardPreprocessing.run(Vertex(4), -31, -30);
        checkFlags(8, flagsData, FORWARD,  DEPARTURE, 1, 0, { true,  false, true,  false, true,  false, true,  false, false, true,  false });
        checkFlags(8, flagsData, FORWARD,  ARRIVAL,   1, 0, { false, true,  false, true,  false, true,  false, false, false, false, true  });
    }

    inline void checkFlags(const int testNum, const RouteFlagsDataType& flagsData, const bool direction, const bool departureOrArrival, const size_t fromCell, const size_t toCell, const std::vector<bool>& expectedFlags) const noexcept {
        const RAPTOR::Data& raptorData = flagsData.getRaptorData(direction);
        const std::string directionString = (direction) ? "Backward" : "Forward";
        const std::string departureOrArrivalString = (departureOrArrival) ? "arrival" : "departure";
        for (size_t stopEvent = 0; stopEvent < raptorData.numberOfStopEvents(); stopEvent++) {
            const size_t actualStopEvent = (direction == BACKWARD) ? flagsData.permutateStopEvent(stopEvent) : stopEvent;
            const int departureTime = (direction == BACKWARD) ? -raptorData.stopEvents[actualStopEvent].arrivalTime : raptorData.stopEvents[actualStopEvent].departureTime;
            const int arrivalTime = (direction == BACKWARD) ? -raptorData.stopEvents[actualStopEvent].departureTime : raptorData.stopEvents[actualStopEvent].arrivalTime;
            UnitTests::check(flagsData.isNecessaryStopEvent(direction, departureOrArrival, actualStopEvent, fromCell, toCell) == expectedFlags[stopEvent], "checkPreprocessing", testNum, ": ", directionString, " ", departureOrArrivalString, " stop event ", stopEvent, " (departure time: ", departureTime, ", arrival time: ", arrivalTime, ") is not flagged as expected (expected: ", expectedFlags[stopEvent], ")");
        }
    }

    inline void addVertices(TransferEdgeList& transferGraph, const size_t numVertices) const noexcept {
        transferGraph.addVertices(numVertices);
        for (const Vertex v : transferGraph.vertices()) {
            transferGraph.set(Coordinates, v, Geometry::Point(Construct::XY, v, v));
        }
    }

    inline CSA::Data buildNetworkCSA1() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 5);
        transferGraph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(3), Vertex(4)).set(TravelTime, 2);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  1);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  1);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  1);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(3)),  1);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> T", "R1", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("B -> T", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 1, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 1, 6, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 2, 4, TripId(1));
        connections.emplace_back(StopId(2), StopId(3), 3, 5, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA2() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 4);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  4);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  4);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  4);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(3)),  4);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> D", "R1", 1);
        trips.emplace_back("B -> D", "R2", 1);
        trips.emplace_back("C -> D", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(2), 0, 5, TripId(0));
        connections.emplace_back(StopId(2), StopId(3), 5, 12, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 0, 7, TripId(1));
        connections.emplace_back(StopId(2), StopId(3), 7, 12, TripId(1));
        connections.emplace_back(StopId(2), StopId(3), 10, 11, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA3() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 4);
        transferGraph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 5);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  1);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  1);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  1);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(3)),  1);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> D", "R1", 1);
        trips.emplace_back("A -> D", "R2", 1);
        trips.emplace_back("B -> D", "R3", 1);
        trips.emplace_back("C -> D", "R4", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 90, 95, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 100, 111, TripId(0));
        connections.emplace_back(StopId(0), StopId(1), 95, 100, TripId(1));
        connections.emplace_back(StopId(1), StopId(3), 100, 110, TripId(1));
        connections.emplace_back(StopId(1), StopId(3), 100, 108, TripId(2));
        connections.emplace_back(StopId(2), StopId(3), 106, 111, TripId(3));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA4() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 3);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  1);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  1);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  1);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> C", "R1", 1);
        trips.emplace_back("B -> C", "R2", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 90, 95, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 100, 120, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 100, 115, TripId(1));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA5() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 8);
        transferGraph.addEdge(Vertex(3), Vertex(2)).set(TravelTime, 0);
        transferGraph.addEdge(Vertex(0), Vertex(5)).set(TravelTime, 0);
        transferGraph.addEdge(Vertex(5), Vertex(6)).set(TravelTime, 0);
        transferGraph.addEdge(Vertex(6), Vertex(7)).set(TravelTime, 0);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  5);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(4)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("S -> C", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 10, TripId(0));
        connections.emplace_back(StopId(0), StopId(3), 0, 10, TripId(1));
        connections.emplace_back(StopId(1), StopId(2), 10, 20, TripId(2));
        connections.emplace_back(StopId(2), StopId(4), 20, 30, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA7() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 5);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  0);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("E", transferGraph.get(Coordinates, Vertex(4)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> B", "R1", 1);
        trips.emplace_back("B -> E", "R2", 1);
        trips.emplace_back("B -> E", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 10, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 10, 20, TripId(1));
        connections.emplace_back(StopId(2), StopId(4), 20, 30, TripId(1));
        connections.emplace_back(StopId(1), StopId(3), 10, 20, TripId(2));
        connections.emplace_back(StopId(3), StopId(4), 20, 30, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline CSA::Data buildNetworkCSA8() const noexcept {
        TransferEdgeList transferGraph;
        addVertices(transferGraph, 6);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  0);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  0);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  0);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(3)),  0);
        stops.emplace_back("E", transferGraph.get(Coordinates, Vertex(4)),  5);
        stops.emplace_back("F", transferGraph.get(Coordinates, Vertex(5)),  0);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> B", "R1", 1);
        trips.emplace_back("B -> C", "R2", 1);
        trips.emplace_back("B -> D", "R3", 1);
        trips.emplace_back("C -> F", "R4", 1);
        trips.emplace_back("D -> E", "R5", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 10, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 10, 20, TripId(1));
        connections.emplace_back(StopId(1), StopId(3), 10, 20, TripId(2));
        connections.emplace_back(StopId(2), StopId(4), 20, 31, TripId(3));
        connections.emplace_back(StopId(4), StopId(5), 31, 40, TripId(3));
        connections.emplace_back(StopId(3), StopId(4), 20, 30, TripId(4));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR(CSA::Data csa) const noexcept {
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        RAPTOR::Data result = RAPTOR::Data::FromIntermediate(inter);
        result.useImplicitDepartureBufferTimes();
        return result;
    }

    inline RAPTOR::TransferGraph buildShortcutGraph(const RAPTOR::Data& raptorData) const noexcept {
        DynamicTransferGraph shortcutGraph;
        shortcutGraph.addVertices(raptorData.numberOfStops());
        Dijkstra<RAPTOR::TransferGraph> dijkstra(raptorData.transferGraph);
        for (const Vertex u : shortcutGraph.vertices()) {
            dijkstra.run(u);
            for (const Vertex v : shortcutGraph.vertices()) {
                if (u == v) continue;
                if (dijkstra.reachable(v)) {
                    shortcutGraph.addEdge(u, v).set(TravelTime, dijkstra.getDistance(v));
                }
            }
        }
        RAPTOR::TransferGraph result;
        Graph::move(std::move(shortcutGraph), result);
        return result;
    }

    inline ::NestedDissection buildNestedDissection1() const noexcept {
        ::NestedDissection partition(5);
        partition.divideCell(0, { Vertex(1) }, { Vertex(0) }, { Vertex(2) , Vertex(3), Vertex(4) });
        return partition;
    }

    inline VertexPartition buildVertexPartition1() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 0, 1, 1, 1});
    }

    inline ::NestedDissection buildNestedDissection2() const noexcept {
        ::NestedDissection partition(4);
        partition.divideCell(0, { Vertex(2) }, { Vertex(0), Vertex(1) }, { Vertex(3) });
        return partition;
    }

    inline VertexPartition buildVertexPartition2() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 0, 0, 1 });
    }

    inline ::NestedDissection buildNestedDissection3() const noexcept {
        ::NestedDissection partition(4);
        partition.divideCell(0, { Vertex(1) }, { Vertex(0) }, { Vertex(2), Vertex(3) });
        return partition;
    }

    inline VertexPartition buildVertexPartition3() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 0, 1, 1 });
    }

    inline ::NestedDissection buildNestedDissection4() const noexcept {
        ::NestedDissection partition(3);
        partition.divideCell(0, { Vertex(1) }, { Vertex(0) }, { Vertex(2) });
        return partition;
    }

    inline VertexPartition buildVertexPartition4() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 0, 1 });
    }

    inline ::NestedDissection buildNestedDissection5() const noexcept {
        ::NestedDissection partition(8);
        partition.divideCell(0, { Vertex(0) }, { Vertex(5), Vertex(6), Vertex(7) }, { Vertex(1), Vertex(2), Vertex(3), Vertex(4) });
        partition.divideCell(1, { Vertex(6) }, { Vertex(7) }, { Vertex(5) });
        partition.divideCell(2, { Vertex(2) }, { Vertex(1), Vertex(3) }, { Vertex(4) });
        return partition;
    }

    inline VertexPartition buildVertexPartition5() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 1, 2, 3, 4, 5, 6, 7 });
    }

    inline ::NestedDissection buildNestedDissection7() const noexcept {
        ::NestedDissection partition(5);
        partition.divideCell(0, { Vertex(1) }, { Vertex(0) }, { Vertex(2), Vertex(3), Vertex(4) });
        return partition;
    }

    inline VertexPartition buildVertexPartition7() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 1, 2, 3, 4 });
    }

    inline ::NestedDissection buildNestedDissection8() const noexcept {
        ::NestedDissection partition(6);
        partition.divideCell(0, { Vertex(4) }, { Vertex(0), Vertex(1), Vertex(2), Vertex(3) }, { Vertex(5) });
        return partition;
    }

    inline VertexPartition buildVertexPartition8() const noexcept {
        return VertexPartition(std::vector<int>{ 0, 1, 2, 3, 4, 5 });
    }

    inline RouteFlagsDataType buildRouteFlagsData(const CSA::Data network, ::NestedDissection coarsePartition, const size_t level, const VertexPartition finePartition) const noexcept {
        const RAPTOR::Data forwardRaptorData = buildNetworkRAPTOR(network);
        ::Permutation stopEventPermutation;
        const RAPTOR::Data backwardRaptorData = forwardRaptorData.reverseNetwork(stopEventPermutation);
        coarsePartition.computeIncidentCells(forwardRaptorData.minTravelTimeGraph());
        RAPTOR::TransferGraph forwardShortcutGraph = buildShortcutGraph(forwardRaptorData);
        RAPTOR::TransferGraph backwardShortcutGraph = buildShortcutGraph(backwardRaptorData);
        return RouteFlagsDataType(forwardRaptorData, backwardRaptorData, stopEventPermutation, forwardShortcutGraph, backwardShortcutGraph, coarsePartition, level, finePartition);
    }
};
}
