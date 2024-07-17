#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"

namespace UnitTests {

class ReverseNetwork {

public:
    template<typename ALGORITHM>
    inline void check() {
        RAPTOR::Data data = buildNetworkRAPTOR();
        RAPTOR::Data reverseData = data.reverseNetwork();
        ALGORITHM algorithm(data, reverseData);
        ALGORITHM reverseAlgorithm(reverseData, data);

        for (const StopId u : data.stops()) {
            for (const StopId v : data.stops()) {
                if (u == v) continue;
                algorithm.run(u, v, 100, 225);
                reverseAlgorithm.run(v, u, -225, -100);
                RAPTOR::Profile result = algorithm.getProfile();
                RAPTOR::Profile reverseResult = reverseAlgorithm.getProfile();
                std::sort(reverseResult.begin(), reverseResult.end(), [&](const auto& a, const auto& b) {
                   return a.arrivalTime > b.arrivalTime || (a.arrivalTime == b.arrivalTime && a.departureTime > b.departureTime);
                });
                UnitTests::check(result.size() == reverseResult.size(), "ReverseNetwork (", u, " -> ", v, "): Profile size should be ", result.size(), " but is ", reverseResult.size());
                for (size_t i = 0; i < result.size(); i++) {
                    if (i >= reverseResult.size()) break;
                    UnitTests::check(result[i].arrivalTime == -reverseResult[i].departureTime, "ReverseNetwork (", u, " -> ", v, "): Departure time should be ", result[i].arrivalTime, " but is ", -reverseResult[i].departureTime);
                    UnitTests::check(result[i].departureTime == -reverseResult[i].arrivalTime, "ReverseNetwork (", u, " -> ", v, "): Arrival time should be ", result[i].departureTime, " but is ", -reverseResult[i].arrivalTime);
                    UnitTests::check(result[i].numberOfTrips == reverseResult[i].numberOfTrips, "ReverseNetwork (", u, " -> ", v, "): Number of trips should be ", result[i].numberOfTrips, " but is ", reverseResult[i].numberOfTrips);
                }
            }
        }
    }

protected:
    inline void addVertices(TransferEdgeList& transferGraph, const size_t numVertices) const noexcept {
        transferGraph.addVertices(numVertices);
        for (const Vertex v : transferGraph.vertices()) {
            transferGraph.set(Coordinates, v, Geometry::Point(Construct::XY, v, v));
        }
    }

    inline CSA::Data buildNetworkCSA() const noexcept {
        EdgeList<WithCoordinates, WithTravelTime> transferGraph;
        addVertices(transferGraph, 6);
        transferGraph.addVertices(6);
        transferGraph.addEdge(Vertex(1), Vertex(5)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(5), Vertex(1)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(5), Vertex(2)).set(TravelTime, 1);
        transferGraph.addEdge(Vertex(2), Vertex(5)).set(TravelTime, 1);

        std::vector<CSA::Stop> stops;
        stops.emplace_back("S", transferGraph.get(Coordinates,  Vertex(0)),  5);
        stops.emplace_back("A", transferGraph.get(Coordinates,  Vertex(1)),  5);
        stops.emplace_back("B", transferGraph.get(Coordinates,  Vertex(2)),  5);
        stops.emplace_back("C", transferGraph.get(Coordinates,  Vertex(3)),  5);
        stops.emplace_back("T", transferGraph.get(Coordinates,  Vertex(4)),  5);

        std::vector<CSA::Trip> trips;
        trips.emplace_back("S -> T", "R1", 1);
        trips.emplace_back("S -> T", "R1", 1);
        trips.emplace_back("S -> T", "R2", 1);
        trips.emplace_back("S -> T", "R2", 1);

        std::vector<CSA::Connection> connections;
        connections.emplace_back(StopId(0), StopId(1), 101, 105, TripId(0));
        connections.emplace_back(StopId(1), StopId(3), 108, 110, TripId(0));
        connections.emplace_back(StopId(3), StopId(4), 118, 120, TripId(0));
        connections.emplace_back(StopId(0), StopId(1), 201, 205, TripId(1));
        connections.emplace_back(StopId(1), StopId(3), 208, 210, TripId(1));
        connections.emplace_back(StopId(3), StopId(4), 218, 220, TripId(1));
        connections.emplace_back(StopId(0), StopId(2), 100, 105, TripId(2));
        connections.emplace_back(StopId(2), StopId(4), 108, 125, TripId(2));
        connections.emplace_back(StopId(0), StopId(2), 200, 205, TripId(3));
        connections.emplace_back(StopId(2), StopId(4), 208, 225, TripId(3));

        return CSA::Data::FromInput(stops, connections, trips, transferGraph);
    }

    inline RAPTOR::Data buildNetworkRAPTOR() const noexcept {
        CSA::Data csa = buildNetworkCSA();
        Intermediate::Data inter = Intermediate::Data::FromCSA(csa);
        return RAPTOR::Data::FromIntermediate(inter);
    }

};
}
