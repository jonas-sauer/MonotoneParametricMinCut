#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Graph/Graph.h"

namespace UnitTests {

class StaticGraph {

private:
    using SimpleGraphType = ::StaticGraph<NoVertexAttributes, NoEdgeAttributes>;
    using SimpleDynamicGraphType = ::DynamicGraph<NoVertexAttributes, NoEdgeAttributes>;
    using TravelTimeAndDistanceGraphType = ::StaticGraph<NoVertexAttributes, WithTravelTimeAndDistance>;
    using TravelTimeAndDistanceDynamicGraphType = ::DynamicGraph<NoVertexAttributes, WithTravelTimeAndDistance>;
    using ReverseEdgeGraphType = ::StaticGraph<NoVertexAttributes, WithReverseEdges>;
    using ReverseEdgeDynamicGraphType = ::DynamicGraph<NoVertexAttributes, WithReverseEdges>;
    using ViaVertexGraphType = ::StaticGraph<WithCoordinates, WithReverseEdgesAndViaVertex>;
    using ViaVertexDynamicGraphType = ::DynamicGraph<WithCoordinates, WithReverseEdgesAndViaVertex>;
    using EdgeListType = ::EdgeList<NoVertexAttributes, WithReverseEdgesAndViaVertex>;

public:
    inline void check() {
        checkAddAndDeleteVertices();
        checkAddEdges();
        checkPermutation();
        checkRecords();
        checkReadWrite();
        checkDimacs();
        checkMoveBetweenStatic();
        checkMoveToEdgeList();
        checkReverse();
    }

protected:
    inline void checkAddAndDeleteVertices() {
        SimpleGraphType graph;
        graph.addVertex();
        UnitTests::check(graph.numVertices() == 1, "Expected 1 vertex, got ", graph.numVertices());
        graph.addVertices(5);
        UnitTests::check(graph.numVertices() == 6, "Expected 6 vertices, got ", graph.numVertices());
        UnitTests::check(graph.vertices().size() == graph.numVertices(), "graph.vertices() has ", graph.vertices().size(), " vertices, expected ", graph.numVertices());
        for (const Vertex v : graph.vertices()) {
            UnitTests::check(graph.isVertex(v), "Vertex ", v, " is not a vertex!");
            UnitTests::check(graph.outDegree(v) == 0, "Vertex ", v, " has outgoing edges!");
        }
        std::vector<Vertex> vertices {Vertex(0), Vertex(2), Vertex(5)};
        graph.deleteVertices([&](const Vertex v) {
            return Vector::contains(vertices, v);
        });
        UnitTests::check(graph.numVertices() == 3, "Expected 3 vertices, got ", graph.numVertices());
        UnitTests::check(!graph.isVertex(Vertex(3)), "Vertex 3 should not exist!");
    }

    inline void checkAddEdges() {
        ReverseEdgeDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(5);
        dynamicGraph.addEdge(Vertex(0), Vertex(1));
        dynamicGraph.addEdge(Vertex(1), Vertex(0));
        dynamicGraph.addEdge(Vertex(2), Vertex(4));

        ReverseEdgeGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);
        UnitTests::check(graph.numEdges() == 3, "Expected 3 edges, got ", graph.numEdges());
        UnitTests::check(graph.findEdge(Vertex(0), Vertex(1)) == 0, "Edge between 0 and 1 should be 0, is ", graph.findEdge(Vertex(0), Vertex(1)));
        UnitTests::check(graph.findEdge(Vertex(1), Vertex(0)) == 1, "Edge between 1 and 0 should be 0, is ", graph.findEdge(Vertex(1), Vertex(0)));
        UnitTests::check(graph.isEdge(Edge(2)), "Edge 2 is not an edge!");
        UnitTests::check(graph.outDegree(Vertex(0)) == 1, "Vertex 0 has out degree ", graph.outDegree(Vertex(0)), ", expected 1!");

        graph.addEdge(Vertex(4), Vertex(2));
        UnitTests::check(graph.findEdge(Vertex(4), Vertex(2)) == 3, "Edge between 4 and 2 should be 3, is ", graph.findEdge(Vertex(4), Vertex(2)));
        UnitTests::check(graph.get(ReverseEdge, Edge(3)) == 2, "Reverse edge of 3 should be 2 but is ", graph.get(ReverseEdge, Edge(3)));
        UnitTests::check(graph.get(ReverseEdge, Edge(2)) == 3, "Reverse edge of 2 should be 3 but is ", graph.get(ReverseEdge, Edge(2)));

        graph.redirectEdge(Edge(3), Vertex(4), Vertex(3));
        UnitTests::check(graph.findEdge(Vertex(4), Vertex(2)) == noEdge, "Edge between 4 and 2 should not exist!");
        UnitTests::check(graph.findEdge(Vertex(4), Vertex(3)) == 3, "Edge between 4 and 3 should be 3, is ", graph.findEdge(Vertex(4), Vertex(3)));
        UnitTests::check(graph.get(ReverseEdge, Edge(3)) == noEdge, "Edge 3 should not have a reverse edge");
        UnitTests::check(graph.get(ReverseEdge, Edge(2)) == noEdge, "Edge 2 should not have a reverse edge");

        graph.deleteVertices([&](const Vertex v) {
            return v % 2 == 1;
        });
        UnitTests::check(graph.numEdges() == 1, "Expected 1 edge, got ", graph.numEdges());
        UnitTests::check(graph.findEdge(Vertex(1), Vertex(2)) == 0, "Edge between 1 and 2 should be 0, is ", graph.findEdge(Vertex(1), Vertex(2)));
    }

    inline void checkPermutation() {
        const int numVertices = 10;

        ViaVertexDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(numVertices);
        for (const Vertex v : dynamicGraph.vertices()) {
            dynamicGraph.set(Coordinates, v, Geometry::Point(Construct::XY, double(v), -double(v)));
            dynamicGraph.addEdge(v, Vertex((v + 2) % numVertices)).set(ViaVertex, Vertex((v + 1) % numVertices));
        }

        ViaVertexGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);
        ViaVertexGraphType copy = graph;

        Permutation perm(Construct::Random, numVertices);
        graph.applyVertexPermutation(perm);

        for (const Vertex v : copy.vertices()) {
            const Vertex newV(perm[v]);
            UnitTests::check(graph.get(Coordinates, newV).x == copy.get(Coordinates, v).x, "X coordinate of vertex ",  newV, " should be ", copy.get(Coordinates, v).x, " but is ", graph.get(Coordinates, newV).x);
            UnitTests::check(graph.get(Coordinates, newV).y == copy.get(Coordinates, v).y, "Y coordinate of vertex ", newV, " should be ", copy.get(Coordinates, v).y, " but is ", graph.get(Coordinates, newV).y);

            for (const Edge e : copy.edgesFrom(v)) {
                const Vertex to = copy.get(ToVertex, e);
                const Vertex newTo(perm[to]);
                const Edge newE = graph.findEdge(newV, newTo);
                UnitTests::check(newE != noEdge, "Could not find edge from ", newV, " to ", newTo);
                UnitTests::check(graph.get(ToVertex, newE) == newTo, "To vertex of edge ", newE, " should be ", newTo, " but is ", graph.get(ToVertex, newE));
                UnitTests::check(graph.get(ViaVertex, newE) == perm[copy.get(ViaVertex, e)], "To vertex of edge ", newE, " should be ", perm[copy.get(ViaVertex, e)], " but is ", graph.get(ViaVertex, newE));
            }
        }
    }

    inline void checkRecords() {
        TravelTimeAndDistanceDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(5);
        for (const Vertex v : dynamicGraph.vertices()) {
            dynamicGraph.addEdge(v, Vertex((v + 3) % 5)).set(TravelTime, v).set(Distance, v * 10);
        }

        TravelTimeAndDistanceGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);

        TravelTimeAndDistanceGraphType::EdgeRecord record = graph.edgeRecord(Edge(2));
        record.forEach([&](auto& value) {
            value++;
        });
        UnitTests::check(record.get(TravelTime) == 3, "Travel time should be 3 but is ", record.get(TravelTime));
        UnitTests::check(record.get(Distance) == 21, "Travel time should be 21 but is ", record.get(Distance));

        TravelTimeAndDistanceGraphType::EdgeHandle e = graph.addEdge(Vertex(4), Vertex(1), record);
        UnitTests::check(e.get(ToVertex) == Vertex(1), "To vertex should be ", 1, " but is ", e.get(ToVertex));
        UnitTests::check(e.get(TravelTime) == 3, "Travel time should be 3 but is ", e.get(TravelTime));
        UnitTests::check(e.get(Distance) == 21, "Travel time should be 21 but is ", e.get(Distance));

        const Edge edge = Edge(0);
        graph.setEdgeAttributes(edge, record);
        UnitTests::check(graph.get(ToVertex, edge) == Vertex(3), "To vertex should be ", 3, " but is ", graph.get(ToVertex, edge));
        UnitTests::check(graph.get(TravelTime, edge) == 3, "Travel time should be 3 but is ", graph.get(TravelTime, edge));
        UnitTests::check(graph.get(Distance, edge) == 21, "Travel time should be 21 but is ", graph.get(Distance, edge));
        UnitTests::check(Graph::edgeToString(graph, edge) == "id: 0, to: 3 | Distance: 21, TravelTime: 3", "Edge string should be \"id: 0, to: 3 | Distance: 21, TravelTime: 3\" but is ", Graph::edgeToString(graph, edge));
    }

    inline void checkReadWrite() {
        TravelTimeAndDistanceDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(20);
        for (const Vertex v : dynamicGraph.vertices()) {
            dynamicGraph.addEdge(v, Vertex((v + 3) % 20)).set(TravelTime, v).set(Distance, v * 10);
            dynamicGraph.addEdge(v, Vertex((v + 7) % 20)).set(TravelTime, v * 2).set(Distance, v * 15);
        }
        dynamicGraph.deleteEdges([&](const Edge e) {
            return e % 2;
        });

        TravelTimeAndDistanceGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);
        graph.writeBinary("UnitTestStaticGraphTemp");
        TravelTimeAndDistanceGraphType copiedGraph;
        copiedGraph.readBinary("UnitTestStaticGraphTemp", ".", false);
        std::remove("UnitTestStaticGraphTemp.distance");
        std::remove("UnitTestStaticGraphTemp.beginOut");
        std::remove("UnitTestStaticGraphTemp.toVertex");
        std::remove("UnitTestStaticGraphTemp.travelTime");
        std::remove("UnitTestStaticGraphTemp.statistics.txt");
        std::remove("UnitTestStaticGraphTemp.attributesSize");
        UnitTests::check(graph.numVertices() == copiedGraph.numVertices(), "Number of vertices is different! (", graph.numVertices(), " vs. ", copiedGraph.numVertices(), ")");
        UnitTests::check(graph.numEdges() == copiedGraph.numEdges(), "Number of vertices is different! (", graph.numEdges(), " vs. ", copiedGraph.numEdges(), ")");
        for (const Vertex u : graph.vertices()) {
            for (const Edge e : graph.edgesFrom(u)) {
                const Vertex v = graph.get(ToVertex, e);
                UnitTests::check(copiedGraph.findEdge(u, v) == e, "Could not find edge from ", u, " to ", v);
                UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(TravelTime, e), "Travel time is different! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(TravelTime, e), ")");
                UnitTests::check(graph.get(Distance, e) == copiedGraph.get(Distance, e), "Distance is different! (", graph.get(Distance, e), " vs. ", copiedGraph.get(Distance, e), ")");
            }
        }
    }

    inline void checkDimacs() {
        TravelTimeAndDistanceDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(20);
        for (const Vertex v : dynamicGraph.vertices()) {
            dynamicGraph.addEdge(v, Vertex((v + 3) % 20)).set(TravelTime, v).set(Distance, v * 10);
            dynamicGraph.addEdge(v, Vertex((v + 7) % 20)).set(TravelTime, v * 2).set(Distance, v * 15);
        }
        dynamicGraph.deleteEdges([&](const Edge e) {
            return e % 2;
        });

        TravelTimeAndDistanceGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);
        Graph::toDimacs("temp", graph, graph.get(TravelTime));
        TravelTimeAndDistanceGraphType copiedGraph;
        Graph::fromDimacs("temp", copiedGraph);
        std::remove("temp.gr");
        std::remove("temp.info");
        UnitTests::check(graph.numVertices() == copiedGraph.numVertices(), "Number of vertices is different! (", graph.numVertices(), " vs. ", copiedGraph.numVertices(), ")");
        UnitTests::check(graph.numEdges() == copiedGraph.numEdges(), "Number of vertices is different! (", graph.numEdges(), " vs. ", copiedGraph.numEdges(), ")");
        for (const Vertex u : graph.vertices()) {
            for (const Edge e : graph.edgesFrom(u)) {
                const Vertex v = graph.get(ToVertex, e);
                UnitTests::check(copiedGraph.findEdge(u, v) == e, "Could not find edge from ", u, " to ", v);
                UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(TravelTime, e), "Travel time is different! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(TravelTime, e), ")");
                UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(Distance, e), "Distance is different from travel time! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(Distance, e), ")");
            }
        }
    }

    inline void checkMoveBetweenStatic() {
        SimpleDynamicGraphType simpleDynamicGraph;
        simpleDynamicGraph.addVertices(10);
        for (const Vertex v : simpleDynamicGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            simpleDynamicGraph.addEdge(v, next);
            simpleDynamicGraph.addEdge(next, v);
        }
        simpleDynamicGraph.deleteEdge(Edge(1));
        SimpleGraphType simpleGraph;
        Graph::move(std::move(simpleDynamicGraph), simpleGraph);
        ViaVertexGraphType augmentedGraph;
        SimpleGraphType copy = simpleGraph;
        Graph::move(std::move(simpleGraph), augmentedGraph);
        UnitTests::check(augmentedGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", augmentedGraph.numVertices());
        UnitTests::check(augmentedGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", augmentedGraph.numEdges());

        for (const Vertex v : augmentedGraph.vertices()) {
            UnitTests::check(augmentedGraph.outDegree(v) == copy.outDegree(v), "Out degree should be ", copy.outDegree(v), " but is ", augmentedGraph.outDegree(v));
            for (const Edge e : augmentedGraph.edgesFrom(v)) {
                const Vertex u = augmentedGraph.get(ToVertex, e);
                UnitTests::check(augmentedGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
                UnitTests::check(u == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", u);
                UnitTests::check(augmentedGraph.get(ReverseEdge, e) == copy.findEdge(u, v), "Reverse edge should be ", copy.findEdge(u, v), " but is ", augmentedGraph.get(ReverseEdge, e));
            }
        }
    }

    inline void checkMoveToEdgeList() {
        SimpleDynamicGraphType simpleDynamicGraph;
        simpleDynamicGraph.addVertices(10);
        for (const Vertex v : simpleDynamicGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            simpleDynamicGraph.addEdge(v, next);
            simpleDynamicGraph.addEdge(next, v);
        }
        simpleDynamicGraph.deleteEdge(Edge(1));
        SimpleGraphType simpleGraph;
        Graph::move(std::move(simpleDynamicGraph), simpleGraph);

        EdgeListType edgeList;
        SimpleGraphType copy = simpleGraph;
        Graph::move(std::move(simpleGraph), edgeList);
        UnitTests::check(edgeList.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", edgeList.numVertices());
        UnitTests::check(edgeList.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", edgeList.numEdges());

        for (const Vertex v : copy.vertices()) {
            for (const Edge e : copy.edgesFrom(v)) {
                const Vertex u = copy.get(ToVertex, e);
                const Edge edgeListEdge = edgeList.findEdge(v, u);
                UnitTests::check(edgeListEdge != noEdge, "Could not find edge from ", v, " to ", u);
                UnitTests::check(edgeList.get(ViaVertex, edgeListEdge) == noVertex, "Edge ", edgeListEdge, " should not have via vertex!");
                UnitTests::check(edgeList.get(ReverseEdge, edgeListEdge) == copy.findEdge(u, v), "Reverse edge should be ", copy.findEdge(u, v), " but is ", edgeList.get(ReverseEdge, edgeListEdge));
            }
        }

        ViaVertexGraphType viaGraph;
        Graph::move(std::move(edgeList), viaGraph);
        UnitTests::check(viaGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", viaGraph.numVertices());
        UnitTests::check(viaGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", viaGraph.numEdges());

        for (const Vertex v : viaGraph.vertices()) {
            UnitTests::check(viaGraph.outDegree(v) == copy.outDegree(v), "Out degree should be ", copy.outDegree(v), " but is ", viaGraph.outDegree(v));
            for (const Edge e : viaGraph.edgesFrom(v)) {
                const Vertex u = viaGraph.get(ToVertex, e);
                const Edge originalEdge = copy.findEdge(v, u);
                UnitTests::check(originalEdge != noEdge, "Could not find original edge from ", v, " to ", u);
                UnitTests::check(viaGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
                const Edge reverseEdge = viaGraph.get(ReverseEdge, e);
                if (reverseEdge == noEdge) {
                    UnitTests::check(copy.findEdge(u, v) == noEdge, "Edge from ", u, " to ", v, " was not copied!");
                } else {
                    UnitTests::check(copy.findEdge(u, v) != noEdge, "Edge from ", u, " to ", v, " does not exist in original graph!");
                }
            }
        }
    }

    inline void checkReverse() {
        ReverseEdgeDynamicGraphType dynamicGraph;
        dynamicGraph.addVertices(10);
        for (const Vertex v : dynamicGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            dynamicGraph.addEdge(v, next);
            dynamicGraph.addEdge(next, v);
        }
        ReverseEdgeGraphType graph;
        Graph::move(std::move(dynamicGraph), graph);
        ReverseEdgeGraphType reverseGraph = graph;
        reverseGraph.revert();
        UnitTests::check(reverseGraph.numVertices() == graph.numVertices(), "Expected ", graph.numVertices(), " vertices but got ", reverseGraph.numVertices());
        UnitTests::check(reverseGraph.numEdges() == graph.numEdges(), "Expected ", graph.numEdges(), " edges but got ", reverseGraph.numEdges());

        for (const Vertex from : reverseGraph.vertices()) {
            for (const Edge e : reverseGraph.edgesFrom(from)) {
                const Vertex to = reverseGraph.get(ToVertex, e);
                const Edge originalEdge = graph.findEdge(to, from);
                UnitTests::check(originalEdge != noEdge, "Could not find original edge!");
                const Edge reverseEdge = reverseGraph.get(ReverseEdge, e);
                UnitTests::check(reverseGraph.get(ToVertex, reverseEdge) == from, "To vertex of reverse edge should be ", from, " but is ", reverseGraph.get(ToVertex, reverseEdge));
            }
        }
    }
};

}
