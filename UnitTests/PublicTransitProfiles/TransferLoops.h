#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

#include "../../DataStructures/CSA/Data.h"

namespace UnitTests {

class ShortTransferLoops {

public:
    template<typename ALGORITHM>
    inline void check(ALGORITHM& algorithm, const std::string& algorithmName) {
        algorithm.run(StopId(0), StopId(3));
        RAPTOR::Profile result = algorithm.getProfile();
        UnitTests::check(result.size() == 1, "ShortTransferLoops (", algorithmName, "): Profile size should be 1 but is ", result.size());
        if (result.size() >= 1) {
            UnitTests::check(result[0].arrivalTime == 120, "ShortTransferLoops (", algorithmName, "): Arrival time should be 120 but is ", result[0].arrivalTime);
            UnitTests::check(result[0].departureTime == 100, "ShortTransferLoops (", algorithmName, "): Departure time should be 100 but is ", result[0].departureTime);
            UnitTests::check(result[0].numberOfTrips == 2, "ShortTransferLoops (", algorithmName, "): Number of trips should be 2 but is ", result[0].numberOfTrips);
        }
    }

    template<typename ALGORITHM>
    inline void checkCSA(const std::string& algorithmName) {
        CSA::Data data = buildNetworkCSA();
        data.sortConnectionsAscendingByDepartureTime();
        ALGORITHM algorithm(data, data.transferGraph);
        check(algorithm, algorithmName);
    }

    template<typename ALGORITHM>
    inline void checkRAPTOR(const std::string& algorithmName) {
        RAPTOR::Data data = buildNetworkRAPTOR();
        RAPTOR::Data reverseData = data.reverseNetwork();
        ALGORITHM algorithm(data, reverseData);
        check(algorithm, algorithmName);
    }

protected:
    inline CSA::Data buildNetworkCSA() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(5);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.5));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.5));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 1.6, 0.2));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 2.0, 0.5));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 1.3, 0.2));

        transferGraph.addEdge(Vertex(1), Vertex(4)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(4), Vertex(1)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(4), Vertex(2)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(2), Vertex(4)).set(TravelTime, 1);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  5);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  5);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  5);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(3)),  5);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("B -> T", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 100, 105, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 108, 115, TripId(1));
        connections.emplace_back(StopId(2), StopId(3), 111, 120, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR() const noexcept {
        CSA::Data csa = buildNetworkCSA();
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        return RAPTOR::Data::FromIntermediate(inter);
    }

};

class LongTransferLoops {

public:
    template<typename ALGORITHM>
    inline void check(ALGORITHM& algorithm, const std::string& algorithmName) {
        algorithm.run(StopId(0), StopId(3));
        RAPTOR::Profile result = algorithm.getProfile();
        UnitTests::check(result.size() == 1, "LongTransferLoops (", algorithmName, "): Profile size should be 1 but is ", result.size());
        if (result.size() >= 1) {
            UnitTests::check(result[0].arrivalTime == 120, "LongTransferLoops (", algorithmName, "): Arrival time should be 120 but is ", result[0].arrivalTime);
            UnitTests::check(result[0].departureTime == 100, "LongTransferLoops (", algorithmName, "): Departure time should be 100 but is ", result[0].departureTime);
            UnitTests::check(result[0].numberOfTrips == 2, "LongTransferLoops (", algorithmName, "): Number of trips should be 2 but is ", result[0].numberOfTrips);
        }
    }

    template<typename ALGORITHM>
    inline void checkCSA(const std::string& algorithmName) {
        CSA::Data data = buildNetworkCSA();
        data.sortConnectionsAscendingByDepartureTime();
        ALGORITHM algorithm(data, data.transferGraph);
        check(algorithm, algorithmName);
    }

    template<typename ALGORITHM>
    inline void checkRAPTOR(const std::string& algorithmName) {
        RAPTOR::Data data = buildNetworkRAPTOR();
        RAPTOR::Data reverseData = data.reverseNetwork();
        ALGORITHM algorithm(data, reverseData);
        check(algorithm, algorithmName);
    }

protected:
    inline CSA::Data buildNetworkCSA() const noexcept {
        TransferEdgeList transferGraph;
        transferGraph.addVertices(6);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.5));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.5));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 1.6, 0.2));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 2.0, 0.5));
        transferGraph.set(Coordinates, Vertex(4), Geometry::Point(Construct::XY, 1.3, 0.2));
        transferGraph.set(Coordinates, Vertex(5), Geometry::Point(Construct::XY, 1.3, 0.8));

        transferGraph.addEdge(Vertex(1), Vertex(4)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(4), Vertex(1)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(4), Vertex(2)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(2), Vertex(4)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(1), Vertex(5)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(5), Vertex(1)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(4), Vertex(5)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(5), Vertex(4)).set(TravelTime, 1);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  5);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  5);
        stops.emplace_back("B", transferGraph.get(Coordinates, Vertex(2)),  5);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(3)),  5);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("B -> T", "R3", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 100, 105, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 108, 115, TripId(1));
        connections.emplace_back(StopId(2), StopId(3), 111, 120, TripId(2));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR() const noexcept {
        CSA::Data csa = buildNetworkCSA();
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        return RAPTOR::Data::FromIntermediate(inter);
    }

};

class TransferLoops {

public:
    template<typename ALGORITHM>
    inline void checkCSA(const std::string& algorithmName) {
        ShortTransferLoops().checkCSA<ALGORITHM>(algorithmName);
        LongTransferLoops().checkCSA<ALGORITHM>(algorithmName);
    }

    template<typename ALGORITHM>
    inline void checkRAPTOR(const std::string& algorithmName) {
        ShortTransferLoops().checkRAPTOR<ALGORITHM>(algorithmName);
        LongTransferLoops().checkRAPTOR<ALGORITHM>(algorithmName);
    }

};

}
