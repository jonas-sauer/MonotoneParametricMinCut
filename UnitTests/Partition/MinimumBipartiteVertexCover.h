#pragma once

#include "../UnitTests.h"

#include "../../Algorithms/BipartiteGraphAlgorithms.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../Helpers/Vector/Vector.h"

namespace UnitTests {

class MinimumBipartiteVertexCover {

public:
    inline void check() {
        checkGraphA();
        checkGraphB();
        checkGraphC();
    }

protected:
    inline void checkGraphA() const noexcept {
        SimpleDynamicGraph graph;
        graph.addVertices(6);
        graph.addEdge(Vertex(0), Vertex(1));
        graph.addEdge(Vertex(0), Vertex(2));
        graph.addEdge(Vertex(0), Vertex(3));
        graph.addEdge(Vertex(1), Vertex(4));
        graph.addEdge(Vertex(1), Vertex(5));
        std::vector<Vertex> vertexCover = minimumBipartiteVertexCover(graph);
        UnitTests::check(vertexCover.size() == 2, "Vertex Cover should have size 2 but has size ", vertexCover.size());
        UnitTests::check(Vector::contains(vertexCover, Vertex(0)), "Vertex Cover should contain vertex 0 but does not!");
        UnitTests::check(Vector::contains(vertexCover, Vertex(1)), "Vertex Cover should contain vertex 1 but does not!");
    }

    inline void checkGraphB() const noexcept {
        SimpleDynamicGraph graph;
        graph.addVertices(8);
        for (Vertex from = Vertex(0); from < 5; from++) {
            for (Vertex to = Vertex(5); to < 8; to++) {
                graph.addEdge(from, to);
            }
        }
        graph.addEdge(Vertex(7), Vertex(1));
        graph.addEdge(Vertex(5), Vertex(3));
        std::vector<Vertex> vertexCover = minimumBipartiteVertexCover(graph);
        UnitTests::check(vertexCover.size() == 3, "Vertex Cover should have size 3 but has size ", vertexCover.size());
        UnitTests::check(Vector::contains(vertexCover, Vertex(5)), "Vertex Cover should contain vertex 5 but does not!");
        UnitTests::check(Vector::contains(vertexCover, Vertex(6)), "Vertex Cover should contain vertex 6 but does not!");
        UnitTests::check(Vector::contains(vertexCover, Vertex(7)), "Vertex Cover should contain vertex 7 but does not!");
    }

    inline void checkGraphC() const noexcept {
        SimpleDynamicGraph graph;
        graph.addVertices(14);
        graph.addEdge(Vertex(0), Vertex(1));
        graph.addEdge(Vertex(1), Vertex(2));
        graph.addEdge(Vertex(2), Vertex(3));
        graph.addEdge(Vertex(3), Vertex(4));
        graph.addEdge(Vertex(4), Vertex(5));
        graph.addEdge(Vertex(5), Vertex(0));
        graph.addEdge(Vertex(6), Vertex(7));
        graph.addEdge(Vertex(8), Vertex(7));
        graph.addEdge(Vertex(8), Vertex(9));
        graph.addEdge(Vertex(6), Vertex(9));
        graph.addEdge(Vertex(10), Vertex(11));
        graph.addEdge(Vertex(10), Vertex(12));
        graph.addEdge(Vertex(10), Vertex(6));
        graph.addEdge(Vertex(10), Vertex(8));
        std::vector<Vertex> vertexCover = minimumBipartiteVertexCover(graph);
        UnitTests::check(vertexCover.size() == 6, "Vertex Cover should have size 3 but has size ", vertexCover.size());
        UnitTests::check(Vector::contains(vertexCover, Vertex(10)), "Vertex Cover should contain vertex 10 but does not!");
    }

};

}
