#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

#include "../../DataStructures/CSA/Data.h"

namespace UnitTests {

class UnrestrictedWalking {

public:
    template<typename ALGORITHM>
    inline void check(ALGORITHM& algorithm, const std::string& algorithmName) {
        algorithm.run(StopId(0), StopId(2), 0, 140);
        RAPTOR::Profile result = algorithm.getProfile();
        UnitTests::check(result.size() == 2, "UnrestrictedWalking (", algorithmName, "): Profile size should be 2 but is ", result.size());
        if (result.size() >= 1) {
            UnitTests::check(result[0].arrivalTime == 115, "UnrestrictedWalking (", algorithmName, "): Arrival time should be 115 but is ", result[0].arrivalTime);
            UnitTests::check(result[0].departureTime == 100, "UnrestrictedWalking (", algorithmName, "): Departure time should be 100 but is ", result[0].departureTime);
            UnitTests::check(result[0].numberOfTrips == 2, "UnrestrictedWalking (", algorithmName, "): Number of trips should be 2 but is ", result[0].numberOfTrips);
        }
        if (result.size() >= 2) {
            UnitTests::check(result[1].arrivalTime == 135, "UnrestrictedWalking (", algorithmName, "): Arrival time should be 135 but is ", result[1].arrivalTime);
            UnitTests::check(result[1].departureTime == 120, "UnrestrictedWalking (", algorithmName, "): Departure time should be 120 but is ", result[1].departureTime);
            UnitTests::check(result[1].numberOfTrips == 2, "UnrestrictedWalking (", algorithmName, "): Number of trips should be 2 but is ", result[1].numberOfTrips);
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
        transferGraph.addVertices(4);

        transferGraph.set(Coordinates, Vertex(0), Geometry::Point(Construct::XY, 0.0, 0.5));
        transferGraph.set(Coordinates, Vertex(1), Geometry::Point(Construct::XY, 1.0, 0.8));
        transferGraph.set(Coordinates, Vertex(2), Geometry::Point(Construct::XY, 2.0, 0.5));
        transferGraph.set(Coordinates, Vertex(3), Geometry::Point(Construct::XY, 1.0, 0.2));

        transferGraph.addEdge(Vertex(0), Vertex(3)).set(TravelTime, 10);
        transferGraph.addEdge(Vertex(3), Vertex(0)).set(TravelTime, 10);
        transferGraph.addEdge(Vertex(3), Vertex(2)).set(TravelTime, 10);
        transferGraph.addEdge(Vertex(2), Vertex(3)).set(TravelTime, 10);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates, Vertex(0)),  5);
        stops.emplace_back("A", transferGraph.get(Coordinates, Vertex(1)),  5);
        stops.emplace_back("T", transferGraph.get(Coordinates, Vertex(2)),  5);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("S -> A", "R1", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("A -> T", "R2", 1);
        trips.emplace_back("A -> T", "R2", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 100, 105, TripId(0));
        connections.emplace_back(StopId(0), StopId(1), 110, 115, TripId(1));
        connections.emplace_back(StopId(0), StopId(1), 120, 125, TripId(2));
        connections.emplace_back(StopId(0), StopId(1), 130, 135, TripId(3));
        connections.emplace_back(StopId(0), StopId(1), 150, 155, TripId(4));
        connections.emplace_back(StopId(1), StopId(2), 110, 115, TripId(5));
        connections.emplace_back(StopId(1), StopId(2), 130, 135, TripId(6));
        connections.emplace_back(StopId(1), StopId(2), 150, 155, TripId(7));
        connections.emplace_back(StopId(1), StopId(2), 170, 175, TripId(8));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR() const noexcept {
        CSA::Data csa = buildNetworkCSA();
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        return RAPTOR::Data::FromIntermediate(inter);
    }

};

}
