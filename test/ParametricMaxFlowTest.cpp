#include <gtest/gtest.h>

#include <random>
#include <fstream>

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

    graph.fromDimacs("../../test/instances/smallTest");

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

void writeParametricTestFile(std::string filename, uint n) {
    std::ifstream file("../../test/instances/" + filename);

    if (file.good()) {
        std::cout << "Test has already been generated" << std::endl;
        return;
    }

    file.close();

    std::ofstream outFile("../../test/instances/" + filename);

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

    outFile << "p pmax " << n << " " << edges.size() << std::endl;
    outFile << "n 1 s" << std::endl;
    outFile << "n 2 t" << std::endl;
    for (auto e : edges) {
        outFile << "a " + std::to_string(e.first.first + 1) + " " + std::to_string(e.first.second + 1) + " " + e.second.toString() << std::endl;
    }
}

TEST(parametricMaxFlow, largeTest) {
    uint n = 1000;

    writeParametricTestFile("randomInstance_" + std::to_string(n) + ".max", n);

    ParametricMaxFlowInstance<pmf::linearFlowFunction> graph;
    graph.fromDimacs("../../test/instances/randomInstance_" + std::to_string(n) + ".max");

    ParametricIBFS<pmf::linearFlowFunction> algo(graph);

    algo.run();

    for (double breakPoint: algo.getBreakpoints())
        std::cout << "Breakpoint at " << breakPoint << std::endl;

    for (double breakPoint : algo.getBreakpoints()) {
        if (breakPoint > 1)
            break;

        std::cout << "Manually checking breakpoint " << breakPoint << std::endl;

        ParametricFlowGraphEdgeList<double> temp;
        temp.addVertices(graph.graph.numVertices());

        for (Vertex v : graph.graph.vertices()) {
            for (Edge e : graph.graph.edgesFrom(v)) {
                Vertex w = graph.graph.get(ToVertex, e);
                if (v < w) {
                    Edge eNew = temp.addEdge(v, w);
                    Edge eRev = temp.addEdge(w, v);
                    temp.set(Capacity, eNew, graph.getCapacity(e).eval(breakPoint));
                    temp.set(Capacity, eRev, graph.getCapacity(graph.graph.get(ReverseEdge, e)).eval(breakPoint));
                    temp.set(ReverseEdge, eNew, eRev);
                    temp.set(ReverseEdge, eRev, eNew);
                }
            }
        }

        StaticMaxFlowInstance<double> staticGraph;
        staticGraph.fromEdgeList(temp, Vertex(0), Vertex(1));

        IBFS<StaticMaxFlowInstance<double>> staticAlgo(staticGraph);

        staticAlgo.run();

        for (Vertex v: staticAlgo.getSourceComponent()) {
            EXPECT_LE(algo.getVertexThetas()[v], breakPoint);
        }

        for (Vertex v: staticAlgo.getSinkComponent()) {
            EXPECT_GT(algo.getVertexThetas()[v], breakPoint);
        }
    }
}


TEST(parametricMaxFlow, profilingTest) {
#ifdef NDEBUG
    uint n = 100000;

    writeParametricTestFile("randomInstance_" + std::to_string(n) + ".max", n);

    ParametricMaxFlowInstance<pmf::linearFlowFunction> graph;
    graph.fromDimacs("../../test/instances/randomInstance_" + std::to_string(n) + ".max");

    ParametricIBFS<pmf::linearFlowFunction> algo(graph);

    algo.run();

#endif
}