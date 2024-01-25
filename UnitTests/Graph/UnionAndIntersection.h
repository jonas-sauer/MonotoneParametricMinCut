#pragma once

#include <vector>

#include "../UnitTests.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Graph/Utils/Union.h"
#include "../../DataStructures/Graph/Utils/Intersection.h"

namespace UnitTests {

class UnionAndIntersection {

public:
    inline void check() {
        DynamicTransferGraph graphA = buildGraphA();
        TransferGraph graphB;
        Graph::move(buildGraphB(), graphB);
        TransferGraph graphC;
        Graph::move(buildGraphC(), graphC);

        auto unionAB = Graph::makeUnion(graphA, graphB);
        auto intersectionABC = Graph::makeIntersection(unionAB, graphC);

        UnitTests::check(unionAB.numVertices() == 6, "Number of vertices in the union should be 6, but is ", unionAB.numVertices(), "!");
        UnitTests::check(intersectionABC.numVertices() == 6, "Number of vertices in the intersection should be 6, but is ", intersectionABC.numVertices(), "!");

        std::vector<std::vector<int>> travelTimes(6, std::vector<int>(6, 0));

        for (const Vertex from : unionAB.vertices()) {
            for (const Edge edge : unionAB.edgesFrom(from)) {
                travelTimes[from][unionAB.get(ToVertex, edge)] = unionAB.get(TravelTime, edge);
            }
        }

        checkTravelTime(travelTimes, 0, 0, 0);
        checkTravelTime(travelTimes, 0, 1, 0);
        checkTravelTime(travelTimes, 0, 2, 0);
        checkTravelTime(travelTimes, 0, 3, 0);
        checkTravelTime(travelTimes, 0, 4, 0);
        checkTravelTime(travelTimes, 0, 5, 0);
        checkTravelTime(travelTimes, 1, 0, 0);
        checkTravelTime(travelTimes, 1, 1, 0);
        checkTravelTime(travelTimes, 1, 2, 102);
        checkTravelTime(travelTimes, 1, 3, 103);
        checkTravelTime(travelTimes, 1, 4, 104);
        checkTravelTime(travelTimes, 1, 5, 105);
        checkTravelTime(travelTimes, 2, 0, 200);
        checkTravelTime(travelTimes, 2, 1, 201);
        checkTravelTime(travelTimes, 2, 2, 202);
        checkTravelTime(travelTimes, 2, 3, 203);
        checkTravelTime(travelTimes, 2, 4, 204);
        checkTravelTime(travelTimes, 2, 5, 0);
        checkTravelTime(travelTimes, 3, 0, 300);
        checkTravelTime(travelTimes, 3, 1, 301);
        checkTravelTime(travelTimes, 3, 2, 302);
        checkTravelTime(travelTimes, 3, 3, 0);
        checkTravelTime(travelTimes, 3, 4, 0);
        checkTravelTime(travelTimes, 3, 5, 305);
        checkTravelTime(travelTimes, 4, 0, 0);
        checkTravelTime(travelTimes, 4, 1, 401);
        checkTravelTime(travelTimes, 4, 2, 0);
        checkTravelTime(travelTimes, 4, 3, 0);
        checkTravelTime(travelTimes, 4, 4, 404);
        checkTravelTime(travelTimes, 4, 5, 0);
        checkTravelTime(travelTimes, 5, 0, 0);
        checkTravelTime(travelTimes, 5, 1, 0);
        checkTravelTime(travelTimes, 5, 2, 0);
        checkTravelTime(travelTimes, 5, 3, 0);
        checkTravelTime(travelTimes, 5, 4, 504);
        checkTravelTime(travelTimes, 5, 5, 0);

        for (std::vector<int>& vector : travelTimes) {
            for (int& value : vector) {
                value = 0;
            }
        }

        for (const Vertex from : intersectionABC.vertices()) {
            for (const Edge edge : intersectionABC.edgesFrom(from)) {
                travelTimes[from][intersectionABC.get(ToVertex, edge)] = intersectionABC.get(TravelTime, edge);
            }
        }

        checkTravelTime(travelTimes, 0, 0, 0);
        checkTravelTime(travelTimes, 0, 1, 0);
        checkTravelTime(travelTimes, 0, 2, 0);
        checkTravelTime(travelTimes, 0, 3, 0);
        checkTravelTime(travelTimes, 0, 4, 0);
        checkTravelTime(travelTimes, 0, 5, 0);
        checkTravelTime(travelTimes, 1, 0, 0);
        checkTravelTime(travelTimes, 1, 1, 0);
        checkTravelTime(travelTimes, 1, 2, 0);
        checkTravelTime(travelTimes, 1, 3, 0);
        checkTravelTime(travelTimes, 1, 4, 0);
        checkTravelTime(travelTimes, 1, 5, 105);
        checkTravelTime(travelTimes, 2, 0, 0);
        checkTravelTime(travelTimes, 2, 1, 201);
        checkTravelTime(travelTimes, 2, 2, 0);
        checkTravelTime(travelTimes, 2, 3, 203);
        checkTravelTime(travelTimes, 2, 4, 0);
        checkTravelTime(travelTimes, 2, 5, 0);
        checkTravelTime(travelTimes, 3, 0, 0);
        checkTravelTime(travelTimes, 3, 1, 0);
        checkTravelTime(travelTimes, 3, 2, 302);
        checkTravelTime(travelTimes, 3, 3, 0);
        checkTravelTime(travelTimes, 3, 4, 0);
        checkTravelTime(travelTimes, 3, 5, 0);
        checkTravelTime(travelTimes, 4, 0, 0);
        checkTravelTime(travelTimes, 4, 1, 0);
        checkTravelTime(travelTimes, 4, 2, 0);
        checkTravelTime(travelTimes, 4, 3, 0);
        checkTravelTime(travelTimes, 4, 4, 0);
        checkTravelTime(travelTimes, 4, 5, 0);
        checkTravelTime(travelTimes, 5, 0, 0);
        checkTravelTime(travelTimes, 5, 1, 0);
        checkTravelTime(travelTimes, 5, 2, 0);
        checkTravelTime(travelTimes, 5, 3, 0);
        checkTravelTime(travelTimes, 5, 4, 0);
        checkTravelTime(travelTimes, 5, 5, 0);
    }

protected:
    inline void checkTravelTime(const std::vector<std::vector<int>>& travelTimes, const size_t from, const size_t to, const int travelTime) const noexcept {
        UnitTests::check(travelTimes[from][to] == travelTime, "Travel time from vertex ", from, " to vertex ", to, " should be ", travelTime, ", but is ", travelTimes[from][to], "!");
    }

    inline DynamicTransferGraph buildGraphA() const noexcept {
        DynamicTransferGraph graph;
        graph.addVertices(6);
        graph.addEdge(Vertex(1), Vertex(5)).set(TravelTime, 105);
        graph.addEdge(Vertex(1), Vertex(3)).set(TravelTime, 103);
        graph.addEdge(Vertex(2), Vertex(2)).set(TravelTime, 202);
        graph.addEdge(Vertex(2), Vertex(3)).set(TravelTime, 203);
        graph.addEdge(Vertex(2), Vertex(4)).set(TravelTime, 204);
        graph.addEdge(Vertex(3), Vertex(0)).set(TravelTime, 300);
        graph.addEdge(Vertex(3), Vertex(2)).set(TravelTime, 302);
        graph.addEdge(Vertex(3), Vertex(5)).set(TravelTime, 305);
        graph.addEdge(Vertex(4), Vertex(1)).set(TravelTime, 401);
        graph.addEdge(Vertex(5), Vertex(4)).set(TravelTime, 504);
        graph.sortEdges(ToVertex);
        return graph;
    }

    inline DynamicTransferGraph buildGraphB() const noexcept {
        DynamicTransferGraph graph;
        graph.addVertices(5);
        graph.addEdge(Vertex(1), Vertex(4)).set(TravelTime, 104);
        graph.addEdge(Vertex(1), Vertex(2)).set(TravelTime, 102);
        graph.addEdge(Vertex(2), Vertex(1)).set(TravelTime, 201);
        graph.addEdge(Vertex(2), Vertex(2)).set(TravelTime, 202);
        graph.addEdge(Vertex(2), Vertex(0)).set(TravelTime, 200);
        graph.addEdge(Vertex(3), Vertex(2)).set(TravelTime, 302);
        graph.addEdge(Vertex(3), Vertex(1)).set(TravelTime, 301);
        graph.addEdge(Vertex(4), Vertex(4)).set(TravelTime, 404);
        graph.addEdge(Vertex(4), Vertex(1)).set(TravelTime, 401);
        graph.sortEdges(ToVertex);
        return graph;
    }

    inline DynamicTransferGraph buildGraphC() const noexcept {
        DynamicTransferGraph graph;
        graph.addVertices(8);
        graph.addEdge(Vertex(1), Vertex(0)).set(TravelTime, 100);
        graph.addEdge(Vertex(1), Vertex(1)).set(TravelTime, 101);
        graph.addEdge(Vertex(1), Vertex(6)).set(TravelTime, 106);
        graph.addEdge(Vertex(1), Vertex(5)).set(TravelTime, 105);
        graph.addEdge(Vertex(2), Vertex(1)).set(TravelTime, 201);
        graph.addEdge(Vertex(2), Vertex(7)).set(TravelTime, 207);
        graph.addEdge(Vertex(2), Vertex(3)).set(TravelTime, 203);
        graph.addEdge(Vertex(2), Vertex(5)).set(TravelTime, 205);
        graph.addEdge(Vertex(3), Vertex(2)).set(TravelTime, 302);
        graph.addEdge(Vertex(3), Vertex(7)).set(TravelTime, 307);
        graph.addEdge(Vertex(5), Vertex(1)).set(TravelTime, 501);
        graph.addEdge(Vertex(5), Vertex(6)).set(TravelTime, 506);
        graph.addEdge(Vertex(5), Vertex(2)).set(TravelTime, 502);
        graph.addEdge(Vertex(4), Vertex(2)).set(TravelTime, 402);
        graph.addEdge(Vertex(4), Vertex(3)).set(TravelTime, 403);
        graph.sortEdges(ToVertex);
        return graph;
    }

};

}
