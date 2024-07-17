#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Graph/Graph.h"

namespace UnitTests {

class DynamicGraph {

private:
    using SimpleGraphType = ::DynamicGraph<NoVertexAttributes, NoEdgeAttributes>;
    using TravelTimeAndDistanceGraphType = ::DynamicGraph<NoVertexAttributes, WithTravelTimeAndDistance>;
    using ReverseEdgeGraphType = ::DynamicGraph<NoVertexAttributes, WithReverseEdges>;
    using ViaVertexGraphType = ::DynamicGraph<WithCoordinates, WithReverseEdgesAndViaVertex>;
    using StaticGraphType = ::StaticGraph<NoVertexAttributes, WithReverseEdgesAndViaVertex>;
    using EdgeListType = ::EdgeList<NoVertexAttributes, WithReverseEdgesAndViaVertex>;

public:
    inline void check() {
        checkAddAndDeleteVertices();
        checkAddAndDeleteEdges();
        checkPermutation();
        checkRecords();
        checkReadWrite();
        checkDimacs();
        checkReverseEdges();
        checkMultiEdges();
        checkRemoveDegreeOneVertices();
        checkContractDegreeTwoVertices();
        checkConnectedComponents();
        checkMoveBetweenDynamic();
        checkMoveToStatic();
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
            UnitTests::check(graph.inDegree(v) == 0, "Vertex ", v, " has incoming edges!");
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
        UnitTests::check(graph.isIsolated(Vertex(0)), "Vertex 0 should be isolated!");
        UnitTests::check(graph.isIsolated(Vertex(1)), "Vertex 1 should be isolated!");
        UnitTests::check(graph.satisfiesInvariants(), "Graph data structure is corrupted!");
        graph.addEdge(Vertex(0), Vertex(1));
        UnitTests::check(!graph.isIsolated(Vertex(0)), "Vertex 0 should not be isolated!");
        UnitTests::check(!graph.isIsolated(Vertex(1)), "Vertex 1 should not be isolated!");
        UnitTests::check(graph.findEdge(Vertex(0), Vertex(1)) == 0, "Edge between 0 and 1 should be 0, is ", graph.findEdge(Vertex(0), Vertex(1)));
        UnitTests::check(graph.satisfiesInvariants(), "Graph data structure is corrupted!");
        graph.addEdge(Vertex(1), Vertex(0));
        UnitTests::check(graph.findEdge(Vertex(1), Vertex(0)) == 1, "Edge between 1 and 0 should be 0, is ", graph.findEdge(Vertex(1), Vertex(0)));
        UnitTests::check(graph.findReverseEdge(Edge(0)) == 1, "Reverse edge of 0 should be 1, is ", graph.findReverseEdge(Edge(0)));
        UnitTests::check(graph.findReverseEdge(Edge(1)) == 0, "Reverse edge of 1 should be 0, is ", graph.findReverseEdge(Edge(1)));
        UnitTests::check(graph.degree(Vertex(0)) == 1, "Vertex 0 has degree ", graph.degree(Vertex(0)), ", expected 1!");
        UnitTests::check(graph.satisfiesInvariants(), "Graph data structure is corrupted!");
        graph.addEdge(Vertex(2), Vertex(4));
        UnitTests::check(graph.numEdges() == 3, "Expected 3 edges, got ", graph.numEdges());
        UnitTests::check(graph.isEdge(Edge(2)), "Edge 2 is not an edge!");
        UnitTests::check(graph.satisfiesInvariants(), "Graph data structure is corrupted!");
        graph.deleteEdge(Edge(1));
        UnitTests::check(!graph.isEdge(Edge(1)), "Edge 1 still exists after deletion!");
        UnitTests::check(graph.numEdges() == 2, "Expected 2 edges, got ", graph.numEdges());
        UnitTests::check(graph.findEdge(Vertex(2), Vertex(4)) == 2, "Edge between 2 and 4 should be 2, is ", graph.findEdge(Vertex(2), Vertex(4)));
        UnitTests::check(graph.satisfiesInvariants(), "Graph data structure is corrupted!");
        graph.packEdges();
        UnitTests::check(graph.numEdges() == 2, "Expected 2 edges, got ", graph.numEdges());
        UnitTests::check(graph.findEdge(Vertex(2), Vertex(4)) == 1, "Edge between 2 and 4 should be 1, is ", graph.findEdge(Vertex(2), Vertex(4)));
        graph.addEdge(Vertex(2), Vertex(1));
        UnitTests::check(graph.outDegree(Vertex(2)) == 2, "Vertex 2 has out degree ", graph.outDegree(Vertex(2)), ", expected 2!");
        UnitTests::check(graph.inDegree(Vertex(1)) == 2, "Vertex 1 has in degree ", graph.inDegree(Vertex(1)), ", expected 1!");
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
        for (const Vertex v : graph.vertices()) {
            UnitTests::check(graph.outDegree(v) == 5, "Vertex ", v, " has out degree ", graph.outDegree(v), ", expected 5!");
            UnitTests::check(graph.inDegree(v) == ((v < 5) ? 10 : 0), "Vertex ", v, " has in degree ", graph.inDegree(v), ", expected ", (v < 5) ? 10 : 0, "!");
            for (const Edge e : graph.edgesFrom(v)) {
                UnitTests::check(graph.isEdge(e), "Vertex ", v, " has invalid outgoing edge ", e, " ");
            }
            for (const Edge e : graph.edgesTo(v)) {
                UnitTests::check(graph.isEdge(e), "Vertex ", v, " has invalid incoming edge ", e, " ");
            }
        }

        graph.addEdge(Vertex(0), Vertex(5));
        UnitTests::check(graph.hasEdge(Vertex(0), Vertex(5)), "Reinsertion failed!");
    }

    inline void checkPermutation() {
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

            for (const Edge e : copy.edgesFrom(v)) {
                const Vertex to = copy.get(ToVertex, e);
                const Vertex newTo(perm[to]);
                const Edge newE = graph.findEdge(newV, newTo);
                UnitTests::check(newE != noEdge, "Could not find edge from ", newV, " to ", newTo);
                UnitTests::check(graph.get(FromVertex, newE) == newV, "From vertex of edge ", newE, " should be ", newV, " but is ", graph.get(FromVertex, newE));
                UnitTests::check(graph.get(ToVertex, newE) == newTo, "To vertex of edge ", newE, " should be ", newTo, " but is ", graph.get(ToVertex, newE));
                UnitTests::check(graph.get(ViaVertex, newE) == perm[copy.get(ViaVertex, e)], "To vertex of edge ", newE, " should be ", perm[copy.get(ViaVertex, e)], " but is ", graph.get(ViaVertex, newE));
            }
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
        UnitTests::check(Graph::edgeToString(graph, e) == "id: 0, to: 3 | Distance: 21, TravelTime: 3", "Edge string should be \"id: 0, to: 3 | Distance: 21, TravelTime: 3\" but is \"", Graph::edgeToString(graph, e), "\"");
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
        graph.writeBinary("UnitTestDynamicGraphTemp");
        TravelTimeAndDistanceGraphType copiedGraph;
        copiedGraph.readBinary("UnitTestDynamicGraphTemp", ".", false);
        std::remove("UnitTestDynamicGraphTemp.distance");
        std::remove("UnitTestDynamicGraphTemp.beginOut");
        std::remove("UnitTestDynamicGraphTemp.fromVertex");
        std::remove("UnitTestDynamicGraphTemp.incomingEdgePointer");
        std::remove("UnitTestDynamicGraphTemp.incomingEdges");
        std::remove("UnitTestDynamicGraphTemp.outDegree");
        std::remove("UnitTestDynamicGraphTemp.toVertex");
        std::remove("UnitTestDynamicGraphTemp.travelTime");
        std::remove("UnitTestDynamicGraphTemp.valid");
        std::remove("UnitTestDynamicGraphTemp.statistics.txt");
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
        for (const Vertex u : graph.vertices()) {
            for (const Edge e : graph.edgesFrom(u)) {
                const Vertex v = graph.get(ToVertex, e);
                const Edge newEdge = copiedGraph.findEdge(u, v);
                UnitTests::check(newEdge != noEdge, "Could not find edge from ", u, " to ", v);
                UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(TravelTime, newEdge), "Travel time is different! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(TravelTime, newEdge), ")");
                UnitTests::check(graph.get(TravelTime, e) == copiedGraph.get(Distance, newEdge), "Distance is different from travel time! (", graph.get(TravelTime, e), " vs. ", copiedGraph.get(Distance, newEdge), ")");
            }
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

    inline void checkMultiEdges() {
        TravelTimeAndDistanceGraphType graph;
        graph.addVertices(3);
        for (const Vertex v : graph.vertices()) {
            graph.addEdge(v, Vertex((v + 1) % 3)).set(TravelTime, v);
            graph.addEdge(v, Vertex((v + 1) % 3)).set(TravelTime, v + 3);
            graph.addEdge(v, Vertex((v + 1) % 3)).set(TravelTime, v + 6);
        }
        UnitTests::check(graph.numEdges() == 9, "Expected 9 edges but got ", graph.numEdges());
        graph.reduceMultiEdgesBy(TravelTime);
        UnitTests::check(graph.numEdges() == 3, "Expected 3 edges but got ", graph.numEdges());
        graph.packEdges();
        for (const Edge e : graph.edges()) {
            UnitTests::check(graph.get(TravelTime, e) == int(e), "Travel time should be ", e, " but is ", graph.get(TravelTime, e));
        }
    }

    inline void checkRemoveDegreeOneVertices() {
        SimpleGraphType graph;
        graph.addVertices(10);
        for (Vertex v(0); v < graph.numVertices() - 1; v++) {
            graph.addEdge(v, Vertex(v+1));
        }
        graph.removeDegreeOneVertices();
        UnitTests::check(graph.empty(), "Graph is not empty!");
    }

    inline void checkContractDegreeTwoVertices() {
        TravelTimeAndDistanceGraphType graph;
        graph.addVertices(10);
        for (Vertex v(0); v < graph.numVertices() - 1; v++) {
            graph.addEdge(v, Vertex(v+1)).set(TravelTime, v).set(Distance, 1);
        }
        graph.contractDegreeTwoVertices();
        graph.packEdges();
        UnitTests::check(graph.numEdges() == 1, "Graph should have 1 edge, has ", graph.numEdges());
        UnitTests::check(graph.get(TravelTime, Edge(0)) == 36, "Travel time should be 36, is ", graph.get(TravelTime, Edge(0)));
        UnitTests::check(graph.get(Distance, Edge(0)) == 9, "Distance should be 9, is ", graph.get(Distance, Edge(0)));
    }

    inline void checkConnectedComponents() {
        SimpleGraphType graph;
        graph.addVertices(10);
        for (Vertex v(0); v < 7; v++) {
            graph.addEdge(v, Vertex(v+1));
            graph.addEdge(Vertex(v+1), v);
        }
        Graph::reduceToBiggestStronglyConnectedComponent<SimpleGraphType, false>(graph);
        UnitTests::check(graph.numVertices() == 8, "Graph should have 8 vertices, has ", graph.numVertices());
        UnitTests::check(graph.numEdges() == 14, "Graph should have 14 edges, has ", graph.numEdges());
    }

    inline void checkMoveBetweenDynamic() {
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

        for (const Vertex v : augmentedGraph.vertices()) {
            UnitTests::check(augmentedGraph.get(BeginOut, v) == copy.get(BeginOut, v), "BeginOut should be ", copy.get(BeginOut, v), " but is ", augmentedGraph.get(BeginOut, v));
            UnitTests::check(augmentedGraph.get(OutDegree, v) == copy.get(OutDegree, v), "OutDegree should be ", copy.get(OutDegree, v), " but is ", augmentedGraph.get(OutDegree, v));
        }

        for (const Edge e : augmentedGraph.edges()) {
            UnitTests::check(augmentedGraph.get(Valid, e) == copy.get(Valid, e), "Validity should be ", copy.get(Valid, e), " but is ", augmentedGraph.get(Valid, e));
            UnitTests::check(augmentedGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            UnitTests::check(augmentedGraph.get(FromVertex, e) == copy.get(FromVertex, e), "From vertex should be ", copy.get(FromVertex, e), " but is ", augmentedGraph.get(FromVertex, e));
            UnitTests::check(augmentedGraph.get(ToVertex, e) == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", augmentedGraph.get(ToVertex, e));
            UnitTests::check(augmentedGraph.get(ReverseEdge, e) == copy.findReverseEdge(e), "Reverse edge should be ", copy.findReverseEdge(e), " but is ", augmentedGraph.get(ReverseEdge, e));
        }
    }

    inline void checkMoveToStatic() {
        SimpleGraphType simpleGraph;
        StaticGraphType staticGraph;
        simpleGraph.addVertices(10);
        for (const Vertex v : simpleGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            simpleGraph.addEdge(v, next);
            simpleGraph.addEdge(next, v);
        }
        simpleGraph.deleteEdge(Edge(1));
        SimpleGraphType copy = simpleGraph;
        copy.packEdges();

        Graph::move(std::move(simpleGraph), staticGraph);
        UnitTests::check(staticGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", staticGraph.numVertices());
        UnitTests::check(staticGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", staticGraph.numEdges());

        for (const Vertex v : staticGraph.vertices()) {
            UnitTests::check(staticGraph.outDegree(v) == copy.outDegree(v), "Out degree should be ", copy.outDegree(v), " but is ", staticGraph.outDegree(v));
        }

        for (const Edge e : staticGraph.edges()) {
            UnitTests::check(staticGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            UnitTests::check(staticGraph.get(ToVertex, e) == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", staticGraph.get(ToVertex, e));
            UnitTests::check(staticGraph.get(ReverseEdge, e) == copy.findReverseEdge(e), "Reverse edge should be ", copy.findReverseEdge(e), " but is ", staticGraph.get(ReverseEdge, e));
        }

        ViaVertexGraphType viaGraph;
        Graph::move(std::move(staticGraph), viaGraph);
        UnitTests::check(viaGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", viaGraph.numVertices());
        UnitTests::check(viaGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", viaGraph.numEdges());

        for (const Vertex v : viaGraph.vertices()) {
            UnitTests::check(viaGraph.get(BeginOut, v) == copy.get(BeginOut, v), "BeginOut should be ", copy.get(BeginOut, v), " but is ", viaGraph.get(BeginOut, v));
        }

        for (const Edge e : viaGraph.edges()) {
            UnitTests::check(viaGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            UnitTests::check(viaGraph.get(ToVertex, e) == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", viaGraph.get(ToVertex, e));
            UnitTests::check(viaGraph.get(ReverseEdge, e) == copy.findReverseEdge(e), "Reverse edge should be ", copy.findReverseEdge(e), " but is ", viaGraph.get(ReverseEdge, e));
        }
    }

    inline void checkMoveToEdgeList() {
        SimpleGraphType simpleGraph;
        EdgeListType edgeList;
        simpleGraph.addVertices(10);
        for (const Vertex v : simpleGraph.vertices()) {
            const Vertex next((v + 1) % 10);
            simpleGraph.addEdge(v, next);
            simpleGraph.addEdge(next, v);
        }
        simpleGraph.deleteEdge(Edge(1));
        SimpleGraphType copy = simpleGraph;
        copy.packEdges();

        Graph::move(std::move(simpleGraph), edgeList);
        UnitTests::check(edgeList.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", edgeList.numVertices());
        UnitTests::check(edgeList.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", edgeList.numEdges());

        for (const Edge e : edgeList.edges()) {
            UnitTests::check(edgeList.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            UnitTests::check(edgeList.get(FromVertex, e) == copy.get(FromVertex, e), "From vertex should be ", copy.get(FromVertex, e), " but is ", edgeList.get(FromVertex, e));
            UnitTests::check(edgeList.get(ToVertex, e) == copy.get(ToVertex, e), "To vertex should be ", copy.get(ToVertex, e), " but is ", edgeList.get(ToVertex, e));
            UnitTests::check(edgeList.get(ReverseEdge, e) == copy.findReverseEdge(e), "Reverse edge should be ", copy.findReverseEdge(e), " but is ", edgeList.get(ReverseEdge, e));
        }

        ViaVertexGraphType viaGraph;
        Graph::move(std::move(edgeList), viaGraph);
        UnitTests::check(viaGraph.numVertices() == copy.numVertices(), "Expected ", copy.numVertices(), " vertices but got ", viaGraph.numVertices());
        UnitTests::check(viaGraph.numEdges() == copy.numEdges(), "Expected ", copy.numEdges(), " edges but got ", viaGraph.numEdges());

        for (const Vertex v : viaGraph.vertices()) {
            UnitTests::check(viaGraph.get(BeginOut, v) == copy.get(BeginOut, v), "BeginOut should be ", copy.get(BeginOut, v), " but is ", viaGraph.get(BeginOut, v));
        }

        for (const Edge e : viaGraph.edges()) {
            const Vertex from = viaGraph.get(FromVertex, e);
            const Vertex to = viaGraph.get(ToVertex, e);
            const Edge originalEdge = copy.findEdge(from, to);
            UnitTests::check(originalEdge != noEdge, "Could not find original edge from ", from, " to ", to);
            UnitTests::check(viaGraph.get(ViaVertex, e) == noVertex, "Edge ", e, " should not have via vertex!");
            const Edge reverseEdge = viaGraph.get(ReverseEdge, e);
            if (reverseEdge == noEdge) {
                UnitTests::check(copy.findEdge(to, from) == noEdge, "Edge from ", to, " to ", from, " was not copied!");
            } else {
                UnitTests::check(copy.findEdge(to, from) != noEdge, "Edge from ", to, " to ", from, " does not exist in original graph!");
            }
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
            UnitTests::check(reverseGraph.get(Valid, e) == graph.get(Valid, originalEdge), "Validity should be ", graph.get(Valid, originalEdge), " but is ", reverseGraph.get(Valid, e));
            const Edge reverseEdge = reverseGraph.get(ReverseEdge, e);
            UnitTests::check(reverseGraph.get(FromVertex, reverseEdge) == to, "From vertex of reverse edge should be ", to, " but is ", reverseGraph.get(FromVertex, reverseEdge) == to);
            UnitTests::check(reverseGraph.get(ToVertex, reverseEdge) == from, "To vertex of reverse edge should be ", from, " but is ", reverseGraph.get(ToVertex, reverseEdge) == from);
        }
    }
};

}
