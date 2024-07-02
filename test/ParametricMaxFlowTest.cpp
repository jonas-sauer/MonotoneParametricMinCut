#include <gtest/gtest.h>

#include <random>

#include "../Algorithms/MaxFlowMinCut/IBFS.h"
#include "../Algorithms/MaxFlowMinCut/ParametricIBFS.h"
#include "../DataStructures/Graph/Utils/Conversion.h"

const uint bigNum = 1000000000;

using EdgeListTemp = EdgeList<NoVertexAttributes, List<Attribute<FromVertex, Vertex>, Attribute<ToVertex, Vertex>, Attribute<Capacity, pmf::linearFlowFunction>, Attribute<pmf::Flow, pmf::linearFlowFunction>, Attribute<ReverseEdge, Edge>>>;

/*
void
addEdge(EdgeListTemp& edgeList, uint v, uint w, const pmf::linearFlowFunction& cap) {
    Edge e, er;

    e = edgeList.addEdge(Vertex(v), Vertex(w));
    er = edgeList.addEdge(Vertex(w), Vertex(v));
    edgeList.set(Capacity, e, cap);
    edgeList.set(Capacity, er, pmf::linearFlowFunction(0));
    edgeList.set(ReverseEdge, e, er);
    edgeList.set(ReverseEdge, er, e);
    edgeList.set(pmf::Flow, e, pmf::linearFlowFunction(0));
    edgeList.set(pmf::Flow, er, pmf::linearFlowFunction(0));
}

pmf::linearFlowGraph
edgeVectorToStaticGraph(uint n, std::vector<std::pair<std::pair<uint, uint>, pmf::linearFlowFunction>> edges) {
    EdgeListTemp edgeList;

    edgeList.addVertices(n);

    for (const auto &e: edges)
        addEdge(edgeList, Vertex(e.first.first), Vertex(e.first.second), e.second);

    pmf::linearFlowGraph staticGraph;

    Graph::copy(edgeList, staticGraph);

    return std::move(staticGraph);
}
 */

TEST(parametricMaxFlow, smallTest) {
    ParametricMaxFlowInstance<pmf::linearFlowFunction> graph;

    graph.fromDimacs("../../test/smallTest");

    ParametricIBFS<pmf::linearFlowFunction> algo(graph);

    algo.run();

    std::vector<double> vertexThetas = algo.getVertexThetas();

    for (uint i = 0; i < 5; i++) {
        std::cout << "Vertex " << i + 1 << " has value " << vertexThetas[i] << std::endl;
    }

    EXPECT_DOUBLE_EQ(vertexThetas[0], 0.0);
    EXPECT_DOUBLE_EQ(vertexThetas[1], 0.5);
    EXPECT_DOUBLE_EQ(vertexThetas[2], 0.3333333333333333);
    EXPECT_DOUBLE_EQ(vertexThetas[3], 0.0);
    EXPECT_DOUBLE_EQ(vertexThetas[4], INFTY);
}

/*
TEST(parametricMaxFlow, largeTest) {
    uint n = 10000;

    std::mt19937 rng(42);
    std::uniform_int_distribution<int> disN(2, bigNum);
    std::uniform_real_distribution<double> disR(0.0, 20.0);

    std::vector<std::pair<std::pair<uint, uint>, pmf::linearFlowFunction>> edges;

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            edges.emplace_back(
                    std::pair<std::pair<uint, uint>, pmf::linearFlowFunction>{{0, i}, pmf::linearFlowFunction(a, b)});
        }
    }

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            edges.emplace_back(std::pair<std::pair<uint, uint>, pmf::linearFlowFunction>{{i, 1},
                                                                                         pmf::linearFlowFunction(-a, a +
                                                                                                                     b)});
        }
    }

    for (uint i = 2; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            uint v = i, w = j;
            if (disN(rng) % 2 == 0)
                std::swap(v, w);

            if (disN(rng) % (n * n) < 20) {
                double b = disR(rng);

                edges.emplace_back(
                        std::pair<std::pair<uint, uint>, pmf::linearFlowFunction>{{v, 2}, pmf::linearFlowFunction(b)});
            }
        }
    }

    pmf::linearFlowGraph graph = edgeVectorToStaticGraph(n, edges);

    ParametricIBFS algo(graph, 0.0, 1.0, Vertex(0), Vertex(1));

    algo.run();

    for (double breakPoint: algo.getBreakpoints())
        std::cout << "Breakpoint at " << breakPoint << std::endl;

    for (uint i = 0; i <= 10; i++) {
        IBFS staticAlgo(graph, Vertex(0), Vertex(1), static_cast<double >(i) * 0.1);

        staticAlgo.run();

        for (Vertex v: staticAlgo.getSourceComponent()) {
            std::cout << "For alpha = " << static_cast<double >(i) * 0.1 << " vertex " << v
                      << " is in the source component" << std::endl;
            EXPECT_LE(algo.getVertexThetas()[v], static_cast<double >(i) * 0.1);
        }

        for (Vertex v: staticAlgo.getSinkComponent()) {
            std::cout << "For alpha = " << static_cast<double >(i) * 0.1 << " vertex " << v
                      << " is in the sink component" << std::endl;
            EXPECT_GT(algo.getVertexThetas()[v], static_cast<double >(i) * 0.1);
        }
    }
}
*/
