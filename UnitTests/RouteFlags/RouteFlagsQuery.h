#pragma once

#include "../../Algorithms/RAPTOR/DijkstraRAPTOR.h"
#include "../../Algorithms/RAPTOR/RouteFlags/Preprocessing/Preprocessing.h"
#include "../../Algorithms/CH/CH.h"
#include "../../Algorithms/CH/Preprocessing/BidirectionalWitnessSearch.h"
#include "../../Algorithms/CH/Preprocessing/CHBuilder.h"
#include "../../Algorithms/RAPTOR/RouteFlags/Query/Query.h"
#include "../UnitTests.h"

#include "../../DataStructures/Partition/NestedDissection.h"
#include "../../DataStructures/Partition/VertexPartition.h"
#include "../../DataStructures/RAPTOR/RouteFlags/RouteFlagsData.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

#include "../../Helpers/MultiThreading.h"

namespace UnitTests {

template<typename FLAGS>
class RouteFlagsQuery {

public:
    using Flags = FLAGS;
    using Type = RouteFlagsQuery<Flags>;
    //using FinePartitionType = RAPTOR::RouteFlags::VertexSeparatorPartition;
    using FinePartitionType = RAPTOR::RouteFlags::EdgeSeparatorPartition;
    using RouteFlagsDataType = RAPTOR::RouteFlags::RouteFlagsData<Flags, FinePartitionType>;
    using PreprocessingType = RAPTOR::RouteFlags::Preprocessing<Flags, FinePartitionType>;

    template<bool DIRECTION>
    using RouteFlagsQueryType = RAPTOR::RouteFlags::Query<Flags, FinePartitionType, true, RAPTOR::NoProfiler, DIRECTION>;
    using RaptorQueryType = RAPTOR::DijkstraRAPTOR<RAPTOR::DijkstraInitialTransfers, RAPTOR::NoProfiler, true, false, false>;

    inline void check() {
        RouteFlagsDataType flagsData = buildRouteFlagsData();
        PreprocessingType preprocessing(flagsData);
        ThreadPinning pinner(1, 1);
        preprocessing.run(pinner, -never, never, false);

        runQuery<FORWARD>(flagsData, 0, 120);
        runQuery<BACKWARD>(flagsData, -120, 0);
    }

protected:
    template<bool DIRECTION>
    inline void runQuery(const RouteFlagsDataType& flagsData, const int minTime, const int maxTime) const noexcept {
        const std::string directionString = (DIRECTION == BACKWARD) ? "backward" : "forward";
        const RAPTOR::Data& forwardRaptorData = flagsData.getRaptorData(DIRECTION);
        const RAPTOR::Data& backwardRaptorData = flagsData.getRaptorData(!DIRECTION);
        const CH::CH chData = computeCH(forwardRaptorData);
        RouteFlagsQueryType<DIRECTION> forwardRouteFlagsQuery(flagsData, chData);
        RaptorQueryType forwardRaptorQuery(forwardRaptorData, forwardRaptorData.transferGraph, backwardRaptorData.transferGraph);

        for (const StopId u : forwardRaptorData.stops()) {
            for (const StopId v : forwardRaptorData.stops()) {
                for (int departureTime = minTime; departureTime <= maxTime; departureTime += 10) {
                    forwardRouteFlagsQuery.run(u, departureTime, v);
                    forwardRaptorQuery.run(u, departureTime, v);
                    std::vector<RAPTOR::ArrivalLabel> routeFlagsArrivals = forwardRouteFlagsQuery.getArrivals(v);
                    std::vector<RAPTOR::ArrivalLabel> raptorArrivals = forwardRaptorQuery.getArrivals(v);
                    UnitTests::check(routeFlagsArrivals.size() == raptorArrivals.size(), "RouteFlagsQuery (", u, " -> ", v, " @ ", departureTime, ", ", directionString, "): Number of arrivals is different!\n", printArrivals(routeFlagsArrivals, raptorArrivals));
                    for (size_t i = 0; i < routeFlagsArrivals.size(); i++) {
                        if (i >= raptorArrivals.size()) break;
                        UnitTests::check(routeFlagsArrivals[i] == raptorArrivals[i], "RouteFlagsQuery (", u, " -> ", v, " @ ", departureTime, ", ", directionString, "): Arrival ", i, " is different!\n", printArrivals(routeFlagsArrivals, raptorArrivals));
                    }
                }
            }
        }
    }

    inline static std::string printArrivals(const std::vector<RAPTOR::ArrivalLabel>& routeFlagsArrivals, const std::vector<RAPTOR::ArrivalLabel>& raptorArrivals) noexcept {
        std::stringstream result;
        result << "Route-Flags: (" << routeFlagsArrivals.size() << ")" << std::endl;
        for (const RAPTOR::ArrivalLabel& arrival : routeFlagsArrivals) {
            result << "    " << arrival << std::endl;
        }
        result << "RAPTOR: (" << raptorArrivals.size() << ")" << std::endl;
        for (const RAPTOR::ArrivalLabel& arrival : raptorArrivals) {
            result << "    " << arrival << std::endl;
        }
        return result.str();
    }

    inline CSA::Data buildNetworkCSA() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(9);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.0));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 0.5, 0.0));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 1.0, 0.0));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 0.5, 0.5));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 1.0, 0.5));
        transferGraph.set(Coordinates, Vertex(5), Geometry::Point(Construct::XY, 0.0, 1.0));
        transferGraph.set(Coordinates, Vertex(6), Geometry::Point(Construct::XY, 0.5, 1.0));
        transferGraph.set(Coordinates, Vertex(7), Geometry::Point(Construct::XY, 1.0, 1.0));
        transferGraph.set(Coordinates, Vertex(8), Geometry::Point(Construct::XY, 0.0, 0.5));

        transferGraph.addEdge(Vertex(0), Vertex(8)).set(TravelTime, 5);
        transferGraph.addEdge(Vertex(8), Vertex(0)).set(TravelTime, 5);
        transferGraph.addEdge(Vertex(8), Vertex(5)).set(TravelTime, 5);
        transferGraph.addEdge(Vertex(5), Vertex(8)).set(TravelTime, 5);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(0)),  1);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(1)),  1);
        stops.emplace_back("C", transferGraph.get(Coordinates, Vertex(2)),  1);
        stops.emplace_back("D", transferGraph.get(Coordinates, Vertex(3)),  1);
        stops.emplace_back("E", transferGraph.get(Coordinates, Vertex(4)),  1);
        stops.emplace_back("F", transferGraph.get(Coordinates, Vertex(5)),  1);
        stops.emplace_back("G", transferGraph.get(Coordinates, Vertex(6)),  1);
        stops.emplace_back("H", transferGraph.get(Coordinates, Vertex(7)),  1);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("A -> H", "R1", 1);
        trips.emplace_back("A -> H", "R1", 1);
        trips.emplace_back("A -> E", "R2", 1);
        trips.emplace_back("A -> E", "R2", 1);
        trips.emplace_back("B -> F", "R3", 1);
        trips.emplace_back("B -> F", "R3", 1);
        trips.emplace_back("H -> A", "R4", 1);
        trips.emplace_back("H -> A", "R4", 1);
        trips.emplace_back("E -> A", "R5", 1);
        trips.emplace_back("E -> A", "R5", 1);
        trips.emplace_back("F -> B", "R6", 1);
        trips.emplace_back("F -> B", "R6", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 0, 5, TripId(0));
        connections.emplace_back(StopId(1), StopId(2), 5, 10, TripId(0));
        connections.emplace_back(StopId(2), StopId(4), 10, 15, TripId(0));
        connections.emplace_back(StopId(4), StopId(7), 15, 20, TripId(0));
        connections.emplace_back(StopId(0), StopId(1), 60, 65, TripId(1));
        connections.emplace_back(StopId(1), StopId(2), 65, 70, TripId(1));
        connections.emplace_back(StopId(2), StopId(4), 70, 75, TripId(1));
        connections.emplace_back(StopId(4), StopId(7), 75, 80, TripId(1));
        connections.emplace_back(StopId(0), StopId(3), 10, 15, TripId(2));
        connections.emplace_back(StopId(3), StopId(4), 15, 20, TripId(2));
        connections.emplace_back(StopId(0), StopId(3), 70, 75, TripId(3));
        connections.emplace_back(StopId(3), StopId(4), 75, 80, TripId(3));
        connections.emplace_back(StopId(1), StopId(3), 10, 15, TripId(4));
        connections.emplace_back(StopId(3), StopId(6), 15, 20, TripId(4));
        connections.emplace_back(StopId(6), StopId(5), 20, 25, TripId(4));
        connections.emplace_back(StopId(1), StopId(3), 70, 75, TripId(5));
        connections.emplace_back(StopId(3), StopId(6), 75, 80, TripId(5));
        connections.emplace_back(StopId(6), StopId(5), 80, 85, TripId(5));
        connections.emplace_back(StopId(7), StopId(4), 30, 35, TripId(6));
        connections.emplace_back(StopId(4), StopId(2), 35, 40, TripId(6));
        connections.emplace_back(StopId(2), StopId(1), 40, 45, TripId(6));
        connections.emplace_back(StopId(1), StopId(0), 45, 50, TripId(6));
        connections.emplace_back(StopId(7), StopId(4), 90, 95, TripId(7));
        connections.emplace_back(StopId(4), StopId(2), 95, 100, TripId(7));
        connections.emplace_back(StopId(2), StopId(1), 100, 105, TripId(7));
        connections.emplace_back(StopId(1), StopId(0), 105, 110, TripId(7));
        connections.emplace_back(StopId(4), StopId(3), 40, 45, TripId(8));
        connections.emplace_back(StopId(3), StopId(0), 45, 50, TripId(8));
        connections.emplace_back(StopId(4), StopId(3), 100, 105, TripId(9));
        connections.emplace_back(StopId(3), StopId(0), 105, 110, TripId(9));
        connections.emplace_back(StopId(5), StopId(6), 40, 45, TripId(10));
        connections.emplace_back(StopId(6), StopId(3), 45, 50, TripId(10));
        connections.emplace_back(StopId(3), StopId(1), 50, 55, TripId(10));
        connections.emplace_back(StopId(5), StopId(6), 100, 105, TripId(11));
        connections.emplace_back(StopId(6), StopId(3), 105, 110, TripId(11));
        connections.emplace_back(StopId(3), StopId(1), 110, 115, TripId(11));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR(CSA::Data csa) const noexcept {
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        RAPTOR::Data result = RAPTOR::Data::FromIntermediate(inter);
        result.useImplicitDepartureBufferTimes();
        return result;
    }

    inline RAPTOR::TransferGraph buildShortcutGraph() const noexcept {
        Intermediate::TransferGraph shortcutGraph;
        shortcutGraph.addVertices(9);
        shortcutGraph.addEdge(Vertex(0), Vertex(5)).set(TravelTime, 10);
        shortcutGraph.addEdge(Vertex(5), Vertex(0)).set(TravelTime, 10);
        RAPTOR::TransferGraph result;
        Graph::move(std::move(shortcutGraph), result);
        return result;
    }

    inline ::NestedDissection buildNestedDissection() const noexcept {
        ::NestedDissection partition(9);
        partition.divideCell(0, { Vertex(1), Vertex(3) }, { Vertex(0), Vertex(5), Vertex(6), Vertex(8) }, { Vertex(2) , Vertex(4), Vertex(7) });
        partition.divideCell(1, { Vertex(8) }, { Vertex(0) }, { Vertex(5), Vertex(6) });
        partition.divideCell(2, { Vertex(4) }, { Vertex(2) }, { Vertex(7) });
        return partition;
    }

    inline VertexPartition buildVertexPartition() const noexcept {
        VertexPartition partition(Vector::id<int>(9));
        return partition;
    }

    inline RouteFlagsDataType buildRouteFlagsData() const noexcept {
        const RAPTOR::Data forwardRaptorData = buildNetworkRAPTOR(buildNetworkCSA());
        ::Permutation stopEventPermutation;
        const RAPTOR::Data backwardRaptorData = forwardRaptorData.reverseNetwork(stopEventPermutation);
        const RAPTOR::TransferGraph shortcutGraph = buildShortcutGraph();
        ::NestedDissection coarsePartition = buildNestedDissection();
        coarsePartition.computeIncidentCells(forwardRaptorData.minTravelTimeGraph());
        const size_t coarseLevel = 2;
        const VertexPartition finePartition = buildVertexPartition();
        RouteFlagsDataType flagsData(forwardRaptorData, backwardRaptorData, stopEventPermutation, shortcutGraph, shortcutGraph, coarsePartition, coarseLevel, finePartition);
        return flagsData;
    }

    /*inline RouteFlagsDataType buildRouteFlagsData() const noexcept {
        const RAPTOR::Data forwardRaptorData = buildNetworkRAPTOR(buildNetworkCSA());
        ::Permutation stopEventPermutation;
        const RAPTOR::Data backwardRaptorData = forwardRaptorData.reverseNetwork(stopEventPermutation);
        const RAPTOR::TransferGraph shortcutGraph = buildShortcutGraph();
        ::NestedDissection partition = buildNestedDissection();
        partition.computeIncidentCells(forwardRaptorData.minTravelTimeGraph());
        const size_t coarseLevel = 2;
        const size_t fineLevel = 2;
        RouteFlagsDataType flagsData(forwardRaptorData, backwardRaptorData, stopEventPermutation, shortcutGraph, shortcutGraph, partition, coarseLevel, partition, fineLevel);
        return flagsData;
    }*/

    inline CH::CH computeCH(RAPTOR::Data raptorData) const noexcept {
        std::vector<bool> contractable(raptorData.transferGraph.numVertices(), false);
        using PROFILER = CH::NoProfiler;
        using WITNESS_SEARCH = CH::BidirectionalWitnessSearch<CHCoreGraph, PROFILER, -1>;
        using GREEDY_KEY_FUNCTION = CH::GreedyKey<WITNESS_SEARCH>;
        using KEY_FUNCTION = CH::PartialKey<WITNESS_SEARCH, GREEDY_KEY_FUNCTION>;
        using STOP_CRITERION = CH::MinCoreSize;
        KEY_FUNCTION keyFunction(contractable, contractable.size(), GREEDY_KEY_FUNCTION(1024, 1024, 0));
        CH::Builder<PROFILER, WITNESS_SEARCH, KEY_FUNCTION, STOP_CRITERION, false, false> chBuilder(raptorData.transferGraph, raptorData.transferGraph[TravelTime], keyFunction, STOP_CRITERION(9));
        chBuilder.run();
        chBuilder.copyCoreToCH();
        return CH::CH(std::move(chBuilder));
    }
};
}
