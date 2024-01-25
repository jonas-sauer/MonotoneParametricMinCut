#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/MaxFlowMinCut/FlowGraphs.h"
#include "../../Algorithms/MaxFlowMinCut/Dinic.h"

namespace UnitTests {

class Dinic {

public:
    inline void check() {
        checkA();
        checkB();
    }

protected:
    inline void checkA() {
        DynamicFlowGraph graph = buildGraphA();
        ::Dinic dinic(graph);

        dinic.run(Vertex(0), Vertex(24));
        UnitTests::check(dinic.getFlow() == 3, "Flow from {0} to {24} should be 3 but is ", dinic.getFlow());

        dinic.run(Vertex(0), Vertex(24));
        UnitTests::check(dinic.getFlow() == 3, "Flow from {0} to {24} should be 3 but is ", dinic.getFlow());

        dinic.run({Vertex(0), Vertex(5), Vertex(10), Vertex(15), Vertex(20)}, {Vertex(4), Vertex(9), Vertex(14), Vertex(19), Vertex(24)});
        UnitTests::check(dinic.getFlow() == 8, "Flow from {0, 5, 10, 15, 20} to {4, 9, 14, 19, 24} should be 8 but is ", dinic.getFlow());

        dinic.run(Vertex(0), Vertex(24));
        UnitTests::check(dinic.getFlow() == 3, "Flow from {0} to {24} should be 3 but is ", dinic.getFlow());
    }

    inline void checkB() {
        DynamicFlowGraph graph = buildGraphB();
        ::Dinic dinic(graph);

        dinic.run(Vertex(6), Vertex(7));
        UnitTests::check(dinic.getFlow() == 2, "Flow from {6} to {7} should be 2 but is ", dinic.getFlow());
    }

    inline DynamicFlowGraph buildGraphA() const noexcept {
        DynamicFlowGraph graph;

        graph.addVertices(25);
        addEdge(graph, 0, 1, 2);
        addEdge(graph, 1, 2, 1);
        addEdge(graph, 2, 3, 2);
        addEdge(graph, 3, 4, 3);

        addEdge(graph, 0, 5, 1);
        addEdge(graph, 1, 6, 1);
        addEdge(graph, 2, 7, 2);
        addEdge(graph, 3, 8, 1);
        addEdge(graph, 4, 9, 2);

        addEdge(graph, 5, 6, 1);
        addEdge(graph, 6, 7, 2);
        addEdge(graph, 7, 8, 12);
        addEdge(graph, 8, 9, 1);

        addEdge(graph, 5, 10, 1);
        addEdge(graph, 6, 11, 1);
        addEdge(graph, 7, 12, 2);
        addEdge(graph, 8, 13, 2);
        addEdge(graph, 9, 14, 1);

        addEdge(graph, 10, 11, 2);
        addEdge(graph, 11, 12, 7);
        addEdge(graph, 12, 13, 3);
        addEdge(graph, 13, 14, 1);

        addEdge(graph, 10, 15, 1);
        addEdge(graph, 11, 16, 1);
        addEdge(graph, 12, 17, 2);
        addEdge(graph, 13, 18, 1);
        addEdge(graph, 14, 19, 1);

        addEdge(graph, 15, 16, 1);
        addEdge(graph, 16, 17, 2);
        addEdge(graph, 17, 18, 1);
        addEdge(graph, 18, 19, 1);

        addEdge(graph, 15, 20, 2);
        addEdge(graph, 16, 21, 1);
        addEdge(graph, 17, 22, 2);
        addEdge(graph, 18, 23, 1);
        addEdge(graph, 19, 24, 2);

        addEdge(graph, 20, 21, 2);
        addEdge(graph, 21, 22, 1);
        addEdge(graph, 22, 23, 3);
        addEdge(graph, 23, 24, 2);

        return graph;
    }

    inline DynamicFlowGraph buildGraphB() const noexcept {
        SimpleDynamicGraph graph;

        graph.addVertices(6);
        graph.addEdge(Vertex(1), Vertex(0));
        graph.addEdge(Vertex(2), Vertex(0));
        graph.addEdge(Vertex(3), Vertex(0));
        graph.addEdge(Vertex(1), Vertex(4));
        graph.addEdge(Vertex(1), Vertex(5));

        graph.addVertices(1);
        graph.addEdge(Vertex(6), Vertex(1));
        graph.addEdge(Vertex(6), Vertex(2));
        graph.addEdge(Vertex(6), Vertex(3));

        graph.addVertices(1);
        graph.addEdge(Vertex(0), Vertex(7));
        graph.addEdge(Vertex(4), Vertex(7));
        graph.addEdge(Vertex(5), Vertex(7));

        return Graph::generateFlowGraph(graph);
    }

    inline void addEdge(DynamicFlowGraph& graph, const int from, const int to, const int capacity) const noexcept {
        graph.findOrAddEdge(Vertex(from), Vertex(to)).set(Capacity, capacity);
        graph.findOrAddEdge(Vertex(to), Vertex(from)).set(Capacity, capacity);
    }

};

}
