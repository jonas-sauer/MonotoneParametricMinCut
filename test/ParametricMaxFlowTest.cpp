#include <gtest/gtest.h>

#include <random>
#include <fstream>

#include "../Algorithms/MaxFlowMinCut/ExcessesIBFS.h"
#include "../Algorithms/MaxFlowMinCut/IBFS.h"
#include "../Algorithms/MaxFlowMinCut/ParametricIBFS.h"
#include "../Algorithms/MaxFlowMinCut/PushRelabel.h"
#include "../Algorithms/MaxFlowMinCut/RestartableIBFS.h"

using FlowEdgeList = ParametricFlowGraphEdgeList<pmf::linearFlowFunction>;
using FlowGraph = ParametricFlowGraph<pmf::linearFlowFunction>;
using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<pmf::linearFlowFunction>;

TEST(parametricMaxFlow, smallTest) {
    ParametricInstance instance;
    instance.fromDimacs("../../test/instances/smallTest");

    ParametricIBFS<pmf::linearFlowFunction> algo(instance);
    algo.run();

    const std::vector<double>& vertexThetas = algo.getVertexThetas();

    for (uint i = 0; i < 5; i++) {
        std::cout << "Vertex " << i + 1 << " has value " << vertexThetas[i] << std::endl;
    }

    EXPECT_DOUBLE_EQ(vertexThetas[0], 0.0);
    EXPECT_DOUBLE_EQ(vertexThetas[1], 0.5);
    EXPECT_DOUBLE_EQ(vertexThetas[2], 0.3333333333333333);
    EXPECT_DOUBLE_EQ(vertexThetas[3], 0.0);
    EXPECT_DOUBLE_EQ(vertexThetas[4], INFTY);
}

const uint bigNum = 1000000000;

void addEdge(FlowEdgeList& edgeList, const Vertex from, const Vertex to, const pmf::linearFlowFunction& cap) {
    edgeList.addEdge(from, to).set(Capacity, cap);
    edgeList.addEdge(to, from).set(Capacity, pmf::linearFlowFunction(0));
}

ParametricInstance createRandomParametricInstance(const uint n) {
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> disN(2, bigNum);
    std::uniform_real_distribution<double> disR(0.0, 20.0);

    FlowEdgeList edgeList;
    edgeList.addVertices(n);

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            addEdge(edgeList, Vertex(0), Vertex(i), pmf::linearFlowFunction(a, b));
        }
    }

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            addEdge(edgeList, Vertex(i), Vertex(1), pmf::linearFlowFunction(-a, a+b));
        }
    }

    for (uint i = 2; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            uint v = i, w = j;
            if (disN(rng) % 2 == 0)
                std::swap(v, w);

            if (disN(rng) % (n * n) < 20) {
                double b = disR(rng);
                addEdge(edgeList, Vertex(v), Vertex(2), pmf::linearFlowFunction(b));
            }
        }
    }

    FlowGraph graph;
    Graph::move(std::move(edgeList), graph);
    return {graph, Vertex(0), Vertex(1), 0.0, 1.0};
}

TEST(parametricMaxFlow, largeTestStatic) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    ParametricWrapper wrapper(instance);
    const size_t steps = 10;
    for (size_t i = 0; i <= steps; i++) {
        const double alpha = instance.alphaMin + static_cast<double>(i) * (instance.alphaMax - instance.alphaMin)/steps;
        wrapper.setAlpha(alpha);
        IBFS<ParametricWrapper> ibfs(wrapper);
        ibfs.run();
        ExcessesIBFS<ParametricWrapper> eibfs(wrapper);
        eibfs.run();
        PushRelabel<ParametricWrapper> prf(wrapper);
        prf.run();

        EXPECT_NEAR(ibfs.getFlowValue(), eibfs.getFlowValue(), pmf::epsilon);
        EXPECT_NEAR(ibfs.getFlowValue(), prf.getFlowValue(), pmf::epsilon);
        EXPECT_EQ(ibfs.getSinkComponent(), eibfs.getSinkComponent());
        EXPECT_EQ(ibfs.getSinkComponent(), prf.getSinkComponent());
    }
}

TEST(parametricMaxFlow, largeTestRestartable) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    ParametricWrapper wrapper(instance);
    PushRelabel<ParametricWrapper> prf(wrapper);
    RestartableIBFS<ParametricWrapper> restartableIbfs(wrapper);
    const size_t steps = 10;
    for (size_t i = 0; i <= steps; i++) {
        const double alpha = instance.alphaMin + static_cast<double>(i) * (instance.alphaMax - instance.alphaMin)/steps;
        wrapper.setAlpha(alpha);
        IBFS<ParametricWrapper> ibfs(wrapper);
        ibfs.run();
        if (i == 0) {
            restartableIbfs.run();
            prf.run();
        } else {
            restartableIbfs.continueAfterUpdate();
            prf.continueAfterUpdate();
        }

        EXPECT_NEAR(ibfs.getFlowValue(), restartableIbfs.getFlowValue(), pmf::epsilon);
        EXPECT_NEAR(ibfs.getFlowValue(), prf.getFlowValue(), pmf::epsilon);
        EXPECT_EQ(ibfs.getSinkComponent(), restartableIbfs.getSinkComponent());
        EXPECT_EQ(ibfs.getSinkComponent(), prf.getSinkComponent());
    }
}

template<typename STATIC_ALGO, typename RESTARTABLE_ALGO>
void largeTest(const uint n) {
    const ParametricInstance instance = createRandomParametricInstance(n);

    ParametricIBFS<pmf::linearFlowFunction> algo(instance);
    algo.run();

    for (const double breakPoint : algo.getBreakpoints())
        std::cout << "Breakpoint at " << breakPoint << std::endl;

    ParametricWrapper wrapper(instance);
    RESTARTABLE_ALGO restartableAlgo(wrapper);

    for (const double breakPoint : algo.getBreakpoints()) {
        std::cout << "Checking breakpoint " << breakPoint << std::endl;

        wrapper.setAlpha(breakPoint);
        STATIC_ALGO staticAlgo(wrapper);
        staticAlgo.run();

        if (breakPoint == instance.alphaMin) {
            restartableAlgo.run();
        } else {
            restartableAlgo.continueAfterUpdate();
        }

        EXPECT_NEAR(algo.getFlowValue(breakPoint), staticAlgo.getFlowValue(), pmf::epsilon);
        EXPECT_NEAR(algo.getFlowValue(breakPoint), restartableAlgo.getFlowValue(), pmf::epsilon);
        EXPECT_EQ(algo.getSinkComponent(breakPoint), staticAlgo.getSinkComponent());
        EXPECT_EQ(algo.getSinkComponent(breakPoint), restartableAlgo.getSinkComponent());

        for (const Vertex v : staticAlgo.getSourceComponent()) {
            EXPECT_LE(algo.getVertexThetas()[v], breakPoint);
        }

        for (const Vertex v : staticAlgo.getSinkComponent()) {
            EXPECT_GT(algo.getVertexThetas()[v], breakPoint);
        }
    }
}

TEST(parametricMaxFlow, largeTestParametricPushRelabel) {
    largeTest<PushRelabel<ParametricWrapper>, PushRelabel<ParametricWrapper>>(1000);
}

TEST(parametricMaxFlow, largeTestParametricRestartableIBFS) {
    largeTest<IBFS<ParametricWrapper>, RestartableIBFS<ParametricWrapper>>(1000);
}