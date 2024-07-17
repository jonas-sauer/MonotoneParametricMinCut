#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Graph/Graph.h"

namespace UnitTests {

class EdgeList {

private:
    using SimpleGraphType = ::EdgeList<NoVertexAttributes, NoEdgeAttributes>;
    using TravelTimeAndDistanceGraphType = ::EdgeList<NoVertexAttributes, WithTravelTimeAndDistance>;
    using ReverseEdgeGraphType = ::EdgeList<NoVertexAttributes, WithReverseEdges>;
    using ViaVertexGraphType = ::EdgeList<WithCoordinates, WithReverseEdgesAndViaVertex>;

public:
    inline void check() {
        checkAddAndDeleteVertices();
        checkAddAndDeleteEdges();
        checkVertexPermutation();
        checkEdgePermutation();
        checkRecords();
        checkReadWrite();
        checkDimacs();
        checkReverseEdges();
        checkMoveBetweenEdgeList();
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
        }
        std::vector<Vertex> vertices {Vertex(0), Vertex(2), Vertex(5)};
        graph.deleteVertices([&](const Vertex v) {
            return Vector::contains(vertices, v);
        });
        UnitTests::check(graph.numVertices() == 3, "Expected 3 vertices, got ", graph.numVertices());
        UnitTests::check(!graph.isVertex(Vertex(3)), "Vertex 3 should not exist!");
    }

    inline void checkAddAndDeleteEdges() {
        SimpleGraphType graph;
        graph.addVertices(5);
        UnitTests::check(graph.findEdge(Vertex(0), Vertex(1)) == noEdge, "Edge between 0 and 1 should not exist!");
        graph.addEdge(Vertex(0), Vertex(1));
        UnitTests::check(graph.findEdge(Vertex(0), Vertex(1)) == 0, "Edge between 0 and 1 should be 0, is ", graph.findEdge(Vertex(0), Vertex(1)));
        graph.addEdge(Vertex(1), Vertex(0));
        UnitTests::check(graph.findEdge(Vertex(1), Vertex(0)) == 1, "Edge between 1 and 0 should be 0, is ", graph.findEdge(Vertex(1), Vertex(0)));
        UnitTests::check(graph.findReverseEdge(Edge(0)) == 1, "Reverse edge of 0 should be 1, is ", graph.findReverseEdge(Edge(0)));
        UnitTests::check(graph.findReverseEdge(Edge(1)) == 0, "Reverse edge of 1 should be 0, is ", graph.findReverseEdge(Edge(1)));
        graph.addEdge(Vertex(2), Vertex(4));
        UnitTests::check(graph.numEdges() == 3, "Expected 3 edges, got ", graph.numEdges());
        UnitTests::check(graph.isEdge(Edge(2)), "Edge 2 is not an edge!");
        graph.deleteEdge(Edge(1));
        UnitTests::check(graph.numEdges() == 2, "Expected 2 edges, got ", graph.numEdges());
        UnitTests::check(graph.findEdge(Vertex(2), Vertex(4)) == 1, "Edge between 2 and 4 should be 1, is ", graph.findEdge(Vertex(2), Vertex(4)));
        graph.clear();
        UnitTests::check(graph.numVertices() == 0, "Graph still has vertices!");
        UnitTests::check(graph.numEdges() == 0, "Graph still has edges!");

        graph.addVertices(10);
        for (const Vertex u : graph.vertices()) {
            for (const Vertex v : graph.vertices()) {
                graph.addEdge(u, v);
            }
        }
        graph.deleteEdges([&](const Edge e) {
           return graph.get(ToVertex, e) >= 5;
        });
        UnitTests::check(graph.numEdges() == 50, "Expected 50 edges, got ", graph.numEdges());
        graph.addEdge(Vertex(0), Vertex(5));
        UnitTests::check(graph.hasEdge(Vertex(0), Vertex(5)), "Reinsertion failed!");

        graph.deleteVertices([&](const Vertex v) {
            return v == 0;
        });
        UnitTests::check(graph.numVertices() == 9, "Expected 9 vertices, got ", graph.numVertices());
        UnitTests::check(graph.numEdges() == 36, "Expected 36 edges, got ", graph.numEdges());
        UnitTests::check(!graph.hasEdge(Vertex(0), Vertex(4)), "Edge between 0 and 4 (formerly 1 and 5) should not exist!");

        graph.addVertices(10);
        graph.deleteIsolatedVertices();
        UnitTests::check(graph.numVertices() == 9, "Expected 9 vertices, got ", graph.numVertices());
    }

    inline void checkVertexPermutation() {
        const int numVertices = 10;

        ViaVertexGraphType graph;
        graph.addVertices(numVertices);
        for (const Vertex v : graph.vertices()) {
            graph.set(Coordinates, v, Geometry::Point(Construct::XY, double(v), -double(v)));
            graph.addEdge(v, Vertex((v + 2) % numVertices)).set(ViaVertex, Vertex((v + 1) % numVertices));
        }

        ViaVertexGraphType copy = graph;

        Permutation perm(Construct::Random, numVertices);
        graph.applyVertexPermutation(perm);

        for (const Vertex v : copy.vertices()) {
            const Vertex newV(perm[v]);
            UnitTests::check(graph.get(Coordinates, newV).x == copy.get(Coordinates, v).x, "X coordinate of vertex ",  newV, " should be ", copy.get(Coordinates, v).x, " but is ", graph.get(Coordinates, newV).x);
            UnitTests::check(graph.get(Coordinates, newV).y == copy.get(Coordinates, v).y, "Y coordinate of vertex ", newV, " should be ", copy.get(Coordinates, v).y, " but is ", graph.get(Coordinates, newV).y);
        }

        for (const Edge e : graph.edges()) {
            UnitTests::check(graph.get(FromVertex, e) == perm[copy.get(FromVertex, e)], "From vertex of edge ", e, " should be ", perm[copy.get(FromVertex, e)], " but is ", graph.get(FromVertex, e));
            UnitTests::check(graph.get(ToVertex, e) == perm[copy.get(ToVertex, e)], "To vertex of edge ", e, " should be ", perm[copy.get(ToVertex, e)], " but is ", graph.get(ToVertex, e));
            UnitTests::check(graph.get(ViaVertex, e) == perm[copy.get(ViaVertex, e)], "Via vertex of edge ", e, " should be ", perm[copy.get(ViaVertex, e)], " but is ", graph.get(ViaVertex, e));
        }
    }

    inline void checkEdgePermutation() {
        const int numVertices = 10;

        ViaVertexGraphType graph;
        graph.addVertices(numVertices);
        for (const Vertex v : graph.vertices()) {
            graph.set(Coordinates, v, Geometry::Point(Construct::XY, double(v), -double(v)));
            graph.addEdge(v, Vertex((v + 2) % numVertices)).set(ViaVertex, Vertex((v + 1) % numVertices));
        }

        ViaVertexGraphType copy = graph;

        Permutation perm(Construct::Random, numVertices);
        graph.applyEdgePermutation(perm);

        for (const Vertex v : copy.vertices()) {
            UnitTests::check(graph.get(Coordinates, v).x == copy.get(Coordinates, v).x, "X coordinate of vertex ",  v, " should be ", copy.get(Coordinates, v).x, " but is ", graph.get(Coordinates, v).x);
            UnitTests::check(graph.get(Coordinates, v).y == copy.get(Coordinates, v).y, "Y coordinate of vertex ", v, " should be ", copy.get(Coordinates, v).y, " but is ", graph.get(Coordinates, v).y);
        }

        for (const Edge e : copy.edges()) {
            Edge newE(perm[e]);
            UnitTests::check(graph.get(FromVertex, newE) == copy.get(FromVertex, e), "From vertex of edge ", newE, " should be ", copy.get(FromVertex, e), " but is ", graph.get(FromVertex, newE));
            UnitTests::check(graph.get(ToVertex, newE) == copy.get(ToVertex, e), "To vertex of edge ", newE, " should be ", copy.get(ToVertex, e), " but is ", graph.get(ToVertex, newE));
            UnitTests::check(graph.get(ViaVertex, newE) == copy.get(ViaVertex, e), "Via vertex of edge ", newE, " should be ", copy.get(ViaVertex, e), " but is ", graph.get(ViaVertex, newE));
        }
    }


    inline void checkRecords() {
        TravelTimeAndDistanceGraphType graph;
        graph.addVertices(5);
        for (const Vertex v : graph.vertices()) {
            graph.addEdge(v, Vertex((v + 3) % 5)).set(TravelTime, v).set(Distance, v * 10);
        }
        TravelTimeAndDistanceGraphType::EdgeRecord record = graph.edgeRecord(Edge(2));
        record.forEach([&](auto& value) {
            value++;
        });
        UnitTests::check(record.get(TravelTime) == 3, "Travel time should be 3 but is ", record.get(TravelTime));
        UnitTests::check(record.get(Distance) == 21, "Travel time should be 21 but is ", record.get(Distance));

        graph.addVertices(5);
        for (Vertex v = Vertex(5); v < Vertex(graph.numVertices()); v++) {
            graph.addEdge(v, Vertex((v + 1) % graph.numVertices()), record);
        }

        for (Edge e = Edge(5); e < Edge(graph.numEdges()); e++) {
            UnitTests::check(graph.get(FromVertex, e) == Vertex(e), "From vertex should be ", e, " but is ", graph.get(FromVertex, e));
            UnitTests::check(graph.get(ToVertex, e) == Vertex((e + 1) % graph.numVertices()), "To vertex should be ", (e + 1) % graph.numVertices(), " but is ", graph.get(ToVertex, e));
            UnitTests::check(graph.get(TravelTime, e) == 3, "Travel time should be 3 but is ", graph.get(TravelTime, e));
            UnitTests::check(graph.get(Distance, e) == 21, "Travel time should be 21 but is ", graph.get(Distance, e));
        }

        const Edge e = Edge(0);
        graph.setEdgeAttributes(e, record);
        UnitTests::check(graph.get(FromVertex, e) == Vertex(0), "From vertex should be ", 0, " but is ", graph.get(FromVertex, e));
        UnitTests::check(graph.get(ToVertex, e) == Vertex(3), "To vertex should be ", 3, " but is ", graph.get(ToVertex, e));
        UnitTests::check(graph.get(TravelTime, e) == 3, "Travel time should be 3 but is ", graph.get(TravelTime, e));
        UnitTests::check(graph.get(Distance, e) == 21, "Travel time should be 21 but is ", graph.get(Distance, e));
        UnitTests::check(Graph::edgeToString(graph, e) == "id: 0, to: 3 | Distance: 21, TravelTime: 3", "Edge string should be \"id: 0, to: 3 | Distance: 21, TravelTime: 3\" but is ", Graph::edgeToString(graph, e));
    }

    inline void checkReadWrite() {
        TravelTimeAndDistanceGraphType graph;
        graph.addVertices(20);
        for (const Vertex v : graph.vertices()) {
            graph.addEdge(v, Vertex((v + 3) % 20)).set(TravelTime, v).set(Distance, v * 10);
            graph.addEdge(v, Vertex((v + 7) % 20)).set(TravelTime, v * 2).set(Distance, v * 15);
        }
        graph.deleteEdges([&](const Edge e) {
            return e % 2;
        });
        graph.writeBinary("UnitTestEdgeListTemp");
        TravelTimeAndDistanceGraphType copiedGraph;
        copiedGraph.readBinary("UnitTestEdgeListTemp", ".", false);
        std::remove("UnitTestEdgeListTemp.distance");
        std::remove("UnitTestEdgeListTemp.fromVertex");
        std::remove("UnitTestEdgeListTemp.info");
        std::remove("UnitTestEdgeListTemp.toVertex");
        std::remove("UnitTestEdgeListTemp.travelTime");
        std::remove("UnitTestEdgeListTemp.statistics.txt");
        std::remove("UnitTestEdgeListTemp.attributesSize");
        UnitTests::check(graph.numVertices() == copiedGraph.numVertices(), "Number of vertices is different! (", graph.numVertices(), " vs. ", copiedGraph.numVertices(), ")");
        UnitTests::check(graph.numEdges() == copiedGraph.numEdges(), "Number of vertices is different! (", graph.numEdges(), " vs. ", copiedGraph.numEdges(), ")");
        for (const Edge e : graph.edges()) {
            UnitTests::check(graph.get(FromVertex, e) == copiedGraph.get(FromVertex, e), "From vertex is different! (", graph.get(FromVertex, e), " vs. ", copiedGraph.get(FromVertex, e), ")");
            UnitTests::check(graph.get(ToVertex, e) == copiedGraph.get(ToVertex, e), "To vertex is different! (", graph.get(ToVertex, e), " vs. ", copiedGraph.get(ToVertex, e), ")");
            UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(TravelTime, e), "Travel time is different! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(TravelTime, e), ")");
            UnitTests::check(graph.get(Distance, e) == copiedGraph.get(Distance, e), "Distance is different! (", graph.get(Distance, e), " vs. ", copiedGraph.get(Distance, e), ")");
        }
    }

    inline void checkDimacs() {
        TravelTimeAndDistanceGraphType graph;
        graph.addVertices(20);
        for (const Vertex v : graph.vertices()) {
            graph.addEdge(v, Vertex((v + 3) % 20)).set(TravelTime, v).set(Distance, v * 10);
            graph.addEdge(v, Vertex((v + 7) % 20)).set(TravelTime, v * 2).set(Distance, v * 15);
        }
        graph.deleteEdges([&](const Edge e) {
            return e % 2;
        });
        Graph::toDimacs("temp", graph, graph.get(TravelTime));
        TravelTimeAndDistanceGraphType copiedGraph;
        Graph::fromDimacs("temp", copiedGraph);
        std::remove("temp.gr");
        std::remove("temp.info");
        UnitTests::check(graph.numVertices() == copiedGraph.numVertices(), "Number of vertices is different! (", graph.numVertices(), " vs. ", copiedGraph.numVertices(), ")");
        UnitTests::check(graph.numEdges() == copiedGraph.numEdges(), "Number of vertices is different! (", graph.numEdges(), " vs. ", copiedGraph.numEdges(), ")");
        for (const Edge e : graph.edges()) {
            UnitTests::check(graph.get(FromVertex, e) == copiedGraph.get(FromVertex, e), "From vertex is different! (", graph.get(FromVertex, e), " vs. ", copiedGraph.get(FromVertex, e), ")");
            UnitTests::check(graph.get(ToVertex, e) == copiedGraph.get(ToVertex, e), "To vertex is different! (", graph.get(ToVertex, e), " vs. ", copiedGraph.get(ToVertex, e), ")");
            UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(TravelTime, e), "Travel time is different! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(TravelTime, e), ")");
            UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(Distance, e), "Distance is different from travel time! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(Distance, e), ")");
        }
    }

    inline void checkReverseEdges() {
        ReverseEdgeGraphType graph;
        graph.addVertices(10);
        for (const Vertex v : graph.vertices()) {
            graph.addEdge(v, Vertex((v + 1) % 10));
        }
        for (const Edge e : graph.edges()) {
            UnitTests::check(graph.get(ReverseEdge, e) == noEdge, "Edge ", e, " should not have reverse edge!");
        }

        ReverseEdgeGraphType::EdgeHandle newEdge1 = graph.addEdge(Vertex(1), Vertex(0));
        UnitTests::check(newEdge1.get(ReverseEdge) == Edge(0), "Reverse edge of ", newEdge1, " should be 0 but is ", newEdge1.get(ReverseEdge));
        UnitTests::check(graph.get(ReverseEdge, Edge(0)) == newEdge1, "Reverse edge of 0 should be ", newEdge1, " but is ", graph.get(ReverseEdge, Edge(0)));

        ReverseEdgeGraphType::EdgeHandle newEdge2 = graph.addReverseEdge(Edge(2));
        UnitTests::check(newEdge2.get(ReverseEdge) == Edge(2), "Reverse edge of ", newEdge2, " should be 2 but is ", newEdge2.get(ReverseEdge));
        UnitTests::check(graph.get(ReverseEdge, Edge(2)) == newEdge2, "Reverse edge of 2 should be ", newEdge2, " but is ", graph.get(ReverseEdge, Edge(2)));

        graph.redirectEdge(Edge(2), Vertex(4));
        UnitTests::check(Edge(graph.findEdge(Vertex(2), Vertex(4))) == Edge(2), "Redirection failed!");
        UnitTests::check(newEdge2.get(ReverseEdge) == noEdge, newEdge2, " should not have reverse edge!");
    }

    inline void checkMoveBetweenEdgeList() {
        SimpleGraphType simpleGraph;
        ViaVertexGraphType augmentedGraph;
        simpleGraph.addVertices(10);
        for (const Vertex v : simpleGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            simpleGraph.addEdge(v, next);
            simpleGraph.addEdge(next, v);
        }
        simpleGraph.deleteEdge(Edge(1));
        SimpleGraphType copy = simpleGraph;
        Graph::move(std::move(simpleGraph), augmentedGraph);
        UnitTests::check(augmentedGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", augmentedGraph.numVertices());
        UnitTests::check(augmentedGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", augmentedGraph.numEdges());

        for (const Edge e : augmentedGraph.edges()) {
            UnitTests::check(augmentedGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            UnitTests::check(augmentedGraph.get(FromVertex, e) == copy.get(FromVertex, e), "From vertex should be ", copy.get(FromVertex, e), " but is ", augmentedGraph.get(FromVertex, e));
            UnitTests::check(augmentedGraph.get(ToVertex, e) == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", augmentedGraph.get(ToVertex, e));
            UnitTests::check(augmentedGraph.get(ReverseEdge, e) == copy.findReverseEdge(e), "Reverse edge should be ", copy.findReverseEdge(e), " but is ", augmentedGraph.get(ReverseEdge, e));
        }
    }

    inline void checkReverse() {
        ReverseEdgeGraphType graph;
        graph.addVertices(10);
        for (const Vertex v : graph.vertices()) {
            const Vertex next((v + 1) % 10);
            graph.addEdge(v, next);
            graph.addEdge(next, v);
        }
        ReverseEdgeGraphType reverseGraph = graph;
        reverseGraph.revert();
        UnitTests::check(reverseGraph.numVertices() == graph.numVertices(), "Expected ", graph.numVertices(), " vertices but got ", reverseGraph.numVertices());
        UnitTests::check(reverseGraph.numEdges() == graph.numEdges(), "Expected ", graph.numEdges(), " edges but got ", reverseGraph.numEdges());

        for (const Edge e : reverseGraph.edges()) {
            const Vertex from = reverseGraph.get(FromVertex, e);
            const Vertex to = reverseGraph.get(ToVertex, e);
            const Edge originalEdge = graph.findEdge(to, from);
            UnitTests::check(originalEdge != noEdge, "Could not find original edge!");
            const Edge reverseEdge = reverseGraph.get(ReverseEdge, e);
            UnitTests::check(reverseGraph.get(FromVertex, reverseEdge) == to, "From vertex of reverse edge should be ", to, " but is ", reverseGraph.get(FromVertex, reverseEdge) == to);
            UnitTests::check(reverseGraph.get(ToVertex, reverseEdge) == from, "To vertex of reverse edge should be ", from, " but is ", reverseGraph.get(ToVertex, reverseEdge) == from);
        }
    }
};

}
