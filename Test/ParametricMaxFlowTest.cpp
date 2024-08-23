#include <gtest/gtest.h>

#include <random>
#include <fstream>

#include "../Helpers/Console/Progress.h"

#include "../Algorithms/ExcessesIBFS.h"
#include "../Algorithms/IBFS.h"
#include "../Algorithms/ParametricIBFS.h"
#include "../Algorithms/PushRelabel.h"
#include "../Algorithms/RestartableIBFS.h"
#include "../Algorithms/DichotomicScheme.h"

using FlowEdgeList = ParametricFlowGraphEdgeList<LinearFlowFunction>;
using FlowGraph = ParametricFlowGraph<LinearFlowFunction>;
using ParametricInstance = ParametricMaxFlowInstance<LinearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<LinearFlowFunction>;


TEST(parametricMaxFlow, smallTest) {
    ParametricInstance instance;
    instance.fromDimacs("../../Test/Instances/smallTest");

    ParametricIBFS<LinearFlowFunction> algo(instance);
    algo.run();

    const std::vector<double>& vertexBreakpoints = algo.getVertexBreakpoints();

    for (uint i = 0; i < 5; i++) {
        std::cout << "Vertex " << i + 1 << " has value " << vertexBreakpoints[i] << std::endl;
    }

    EXPECT_DOUBLE_EQ(vertexBreakpoints[0], 0.0);
    EXPECT_DOUBLE_EQ(vertexBreakpoints[1], 0.5);
    EXPECT_DOUBLE_EQ(vertexBreakpoints[2], 0.3333333333333333);
    EXPECT_DOUBLE_EQ(vertexBreakpoints[3], 0.0);
    EXPECT_DOUBLE_EQ(vertexBreakpoints[4], INFTY);
}

inline void addEdge(FlowEdgeList& edgeList, const Vertex from, const Vertex to, const LinearFlowFunction& cap) noexcept {
    edgeList.addEdge(from, to).set(Capacity, cap);
    edgeList.addEdge(to, from).set(Capacity, LinearFlowFunction(0));
}

inline ParametricInstance createRandomParametricInstance(const uint n) noexcept {
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> disN(2, 1000000000);
    std::uniform_real_distribution<double> disR(0.0, 20.0);

    FlowEdgeList edgeList;
    edgeList.addVertices(n);

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            addEdge(edgeList, Vertex(0), Vertex(i), LinearFlowFunction(a, b));
        }
    }

    for (uint i = 2; i < n; i++) {
        if (disN(rng) % 5 == 0) {
            double a, b;
            a = disR(rng);
            b = disR(rng);
            addEdge(edgeList, Vertex(i), Vertex(1), LinearFlowFunction(-a, a + b));
        }
    }

    for (uint i = 2; i < n; i++) {
        for (uint j = i + 1; j < n; j++) {
            uint v = i, w = j;
            if (disN(rng) % 2 == 0)
                std::swap(v, w);

            if (disN(rng) % (n * n) < 20) {
                double b = disR(rng);
                addEdge(edgeList, Vertex(v), Vertex(2), LinearFlowFunction(b));
            }
        }
    }

    FlowGraph graph;
    Graph::move(std::move(edgeList), graph);
    return {graph, Vertex(0), Vertex(1), 0.0, 1.0};
}

template<typename ALGORITHM1, typename ALGORITHM2>
inline void compareAlgorithmResults(const ALGORITHM1& algo1, const ALGORITHM2& algo2, const double tolerance = 1e-9) noexcept {
    EXPECT_NEAR(algo1.getFlowValue(), algo2.getFlowValue(), tolerance);
    EXPECT_EQ(algo1.getSinkComponent(), algo2.getSinkComponent());
}

inline void validateStaticAlgorithms(const ParametricInstance& instance, const size_t steps) noexcept {
    ParametricWrapper wrapper(instance);
    for (size_t i = 0; i <= steps; i++) {
        const double alpha = instance.alphaMin + static_cast<double>(i) * (instance.alphaMax - instance.alphaMin)/steps;
        wrapper.setAlpha(alpha);
        IBFS<ParametricWrapper> ibfs(wrapper);
        ibfs.run();
        ExcessesIBFS<ParametricWrapper> eibfs(wrapper);
        eibfs.run();
        PushRelabel<ParametricWrapper> prf(wrapper);
        prf.run();
        compareAlgorithmResults(ibfs, eibfs);
        compareAlgorithmResults(ibfs,prf);
    }
}

inline void validateRestartableAlgorithms(const ParametricInstance& instance, const size_t steps) noexcept {
    ParametricWrapper wrapper(instance);
    PushRelabel<ParametricWrapper> restartablePRF(wrapper);
    RestartableIBFS<ParametricWrapper> restartableIBFS(wrapper);
    for (size_t i = 0; i <= steps; i++) {
        const double alpha = instance.alphaMin + static_cast<double>(i) * (instance.alphaMax - instance.alphaMin)/steps;
        wrapper.setAlpha(alpha);
        IBFS<ParametricWrapper> ibfs(wrapper);
        ibfs.run();
        if (i == 0) {
            restartableIBFS.run();
            restartablePRF.run();
        } else {
            restartableIBFS.continueAfterUpdate();
            restartablePRF.continueAfterUpdate();
        }
        compareAlgorithmResults(ibfs, restartableIBFS);
        compareAlgorithmResults(ibfs,restartablePRF);
    }
}

template<typename PARAMETRIC_ALGO, typename STATIC_ALGO>
inline void compareParametricAlgorithmResults(const PARAMETRIC_ALGO& parametricAlgo, const STATIC_ALGO& staticAlgo, const double alpha, const double tolerance) noexcept {
    EXPECT_NEAR(parametricAlgo.getFlowValue(alpha), staticAlgo.getFlowValue(), tolerance);
    EXPECT_EQ(parametricAlgo.getSinkComponent(alpha), staticAlgo.getSinkComponent());
}

template<typename STATIC_ALGO, typename RESTARTABLE_ALGO>
inline void validateParametricIBFS(const ParametricInstance& instance, const double tolerance) {
    ParametricIBFS<LinearFlowFunction> algo(instance);
    algo.run();
    ParametricWrapper wrapper(instance);
    RESTARTABLE_ALGO restartableAlgo(wrapper);

    const std::vector<double>& breakpoints = algo.getBreakpoints();
    Progress progress(breakpoints.size());
    for (const double breakpoint : breakpoints) {
        wrapper.setAlpha(breakpoint);
        STATIC_ALGO staticAlgo(wrapper);
        staticAlgo.run();
        if (breakpoint == instance.alphaMin) {
            restartableAlgo.run();
        } else {
            restartableAlgo.continueAfterUpdate();
        }
        compareAlgorithmResults(staticAlgo, restartableAlgo, tolerance);
        compareParametricAlgorithmResults(algo, staticAlgo, breakpoint, tolerance);
        compareParametricAlgorithmResults(algo, restartableAlgo, breakpoint, tolerance);
        progress++;
    }
    progress.finished();
}

template<typename RESTARTABLE_ALGO>
inline void validateParametricIBFSFast(const ParametricInstance& instance, const double tolerance) {
    ParametricIBFS<LinearFlowFunction> algo(instance);
    algo.run();
    ParametricWrapper wrapper(instance);
    RESTARTABLE_ALGO restartableAlgo(wrapper);

    const std::vector<double>& breakpoints = algo.getBreakpoints();
    Progress progress(breakpoints.size());
    for (const double breakpoint : breakpoints) {
        wrapper.setAlpha(breakpoint);
        if (breakpoint == instance.alphaMin) {
            restartableAlgo.run();
        } else {
            restartableAlgo.continueAfterUpdate();
        }
        compareParametricAlgorithmResults(algo, restartableAlgo, breakpoint, tolerance);
        progress++;
    }
    progress.finished();
}

template<typename STATIC_ALGO>
inline void validateDichotomicScheme(const ParametricInstance& instance, const double precision, const double tolerance) {
    using SearchAlgorithm = IBFS<DichotomicSchemeMaxFlowWrapper<LinearFlowFunction>>;
    using Dichotomic = DichotomicScheme<LinearFlowFunction, SearchAlgorithm>;
    Dichotomic algo(instance, precision);
    algo.run();
    ParametricWrapper wrapper(instance);

    const std::vector<double>& breakpoints = algo.getBreakpoints();
    Progress progress(breakpoints.size());
    for (const double breakpoint : breakpoints) {
        wrapper.setAlpha(breakpoint);
        STATIC_ALGO staticAlgo(wrapper);
        staticAlgo.run();
        EXPECT_NEAR(algo.getFlowValue(breakpoint), staticAlgo.getFlowValue(), tolerance);
        EXPECT_EQ(algo.getSinkComponent(breakpoint), staticAlgo.getSinkComponent());
        progress++;
    }
    progress.finished();
}

inline void compareParametricAlgorithms(const ParametricInstance& instance, const double precision, const double tolerance) {
    using SearchAlgorithm = IBFS<DichotomicSchemeMaxFlowWrapper<LinearFlowFunction>>;
    using Dichotomic = DichotomicScheme<LinearFlowFunction, SearchAlgorithm>;
    Dichotomic dichotomicScheme(instance, precision);
    dichotomicScheme.run();
    ParametricIBFS<LinearFlowFunction> parametricIBFS(instance);
    parametricIBFS.run();

    const std::vector<double>& parametricBreakpoints = parametricIBFS.getBreakpoints();
    const std::vector<double>& dichotomicBreakpoints = dichotomicScheme.getBreakpoints();
    std::cout << "Parametric IBFS: " << parametricBreakpoints.size() << " breakpoints" << std::endl;
    std::cout << "Dichotomic scheme: " << dichotomicBreakpoints.size() << " breakpoints" << std::endl;

    std::cout << "Check validity of parametric IBFS" << std::endl;
    Progress progress(dichotomicBreakpoints.size());
    for (const double breakpoint : dichotomicBreakpoints) {
        EXPECT_LE(parametricIBFS.getFlowValue(breakpoint), dichotomicScheme.getFlowValue(breakpoint) + tolerance);
        progress++;
    }

    std::cout << "Check validity of dichotomic scheme" << std::endl;
    progress.init(parametricBreakpoints.size());
    size_t i = 0;
    for (const double breakpoint : parametricBreakpoints) {
        while (i < dichotomicBreakpoints.size() && dichotomicBreakpoints[i] < breakpoint) i++;
        double dichotomicValue = INFTY;
        if (i < dichotomicBreakpoints.size()) {
            dichotomicValue = std::min(dichotomicValue, dichotomicScheme.getFlowValue(dichotomicBreakpoints[i]));
        }
        if (i >= 1) {
            dichotomicValue = std::min(dichotomicValue, dichotomicScheme.getFlowValue(dichotomicBreakpoints[i - 1]));
        }
        const double ibfsValue = parametricIBFS.getFlowValue(dichotomicBreakpoints[i]);
        EXPECT_LE(dichotomicValue, (1 + precision) * ibfsValue);
        progress++;
    }
}

TEST(parametricMaxFlow, randomStatic) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    validateStaticAlgorithms(instance, 100);
}

TEST(parametricMaxFlow, randomRestartable) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    validateRestartableAlgorithms(instance, 100);
}

TEST(parametricMaxFlow, randomParametricIBFSPushRelabel) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    validateParametricIBFS<PushRelabel<ParametricWrapper>, PushRelabel<ParametricWrapper>>(instance, 1e-9);
}

TEST(parametricMaxFlow, randomParametricIBFSRestartableIBFS) {
    const ParametricInstance instance = createRandomParametricInstance(1000);
    validateParametricIBFS<IBFS<ParametricWrapper>, RestartableIBFS<ParametricWrapper>>(instance, 1e-9);
}

const std::string ahremPath = "../../Data/Instances/Aggregation/small/ahrem";
const std::string bonnPath = "../../Data/Instances/Aggregation/large/bonn";

TEST(parametricMaxFlow, ahremParametricIBFS) {
    const ParametricInstance instance(ahremPath);
    validateParametricIBFS<PushRelabel<ParametricWrapper>, PushRelabel<ParametricWrapper>>(instance, 1e-4);
}

TEST(parametricMaxFlow, ahremDichotomic) {
    const ParametricInstance instance(ahremPath);
    validateDichotomicScheme<PushRelabel<ParametricWrapper>>(instance, 1e-16, 1e-5);
}

TEST(parametricMaxFlow, ahremCompare) {
    const ParametricInstance instance(ahremPath);
    compareParametricAlgorithms(instance, 1e-16, 1e-5);
}

TEST(parametricMaxFlow, bonnParametricIBFS) {
    const ParametricInstance instance(bonnPath);
    validateParametricIBFSFast<PushRelabel<ParametricWrapper>>(instance, 1e-3);
}

TEST(parametricMaxFlow, bonnDichotomic) {
    const ParametricInstance instance(bonnPath);
    validateDichotomicScheme<PushRelabel<ParametricWrapper>>(instance, 1e-16, 1e-4);
}

TEST(parametricMaxFlow, bonnCompare) {
    const ParametricInstance instance(bonnPath);
    compareParametricAlgorithms(instance, 1e-16, 1e-5);
}