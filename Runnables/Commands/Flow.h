#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Console/ProgressBar.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Timer.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"
#include "../../Algorithms/MaxFlowMinCut/ChordScheme.h"
#include "../../Algorithms/MaxFlowMinCut/PushRelabel.h"
#include "../../Algorithms/MaxFlowMinCut/IBFS.h"
#include "../../Algorithms/MaxFlowMinCut/ExcessesIBFS.h"
#include "../../Algorithms/MaxFlowMinCut/RestartableIBFS.h"
#include "../../Algorithms/MaxFlowMinCut/ParametricIBFS.h"
#include "../../Helpers/Console/Progress.h"

using namespace Shell;

using StaticInstance = StaticMaxFlowInstance<int>;
using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ParametricWrapper = RestartableMaxFlowWrapper<pmf::linearFlowFunction>;

class LoadMaxFlowInstanceFromDimacs : public ParameterizedCommand {

public:
    LoadMaxFlowInstanceFromDimacs(BasicShell& shell) :
        ParameterizedCommand(shell, "loadMaxFlowInstanceFromDimacs", "Load the given max-flow instance in DIMACS format.") {
        addParameter("Input file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        StaticInstance instance;
        instance.fromDimacs(getParameter("Input file"));
        Graph::printInfo(instance.graph);
        instance.graph.printAnalysis();
        instance.serialize(getParameter("Output file"));
    }
};

class MakeStaticMaxFlowInstanceParametric : public ParameterizedCommand {

public:
    MakeStaticMaxFlowInstanceParametric(BasicShell& shell) :
        ParameterizedCommand(shell, "makeStaticMaxFlowInstanceParametric", "Converts the given static max-flow instance to a parametric one.") {
        addParameter("Input file");
        addParameter("Output file");
        addParameter("Steps", "0");
    }

    virtual void execute() noexcept {
        StaticInstance staticInstance(getParameter("Input file"));
        ParametricInstance instance(staticInstance, getParameter<int>("Steps"));
        Graph::printInfo(instance.graph);
        instance.graph.printAnalysis();
        instance.serialize(getParameter("Output file"));
    }
};

class LoadParametricMaxFlowInstanceFromDimacs : public ParameterizedCommand {

public:
    LoadParametricMaxFlowInstanceFromDimacs(BasicShell& shell) :
    ParameterizedCommand(shell, "loadParametricMaxFlowInstanceFromDimacs", "Load the given parametric max-flow instance in DIMACS format.") {
        addParameter("Input file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        ParametricInstance instance;
        instance.fromDimacs(getParameter("Input file"));
        Graph::printInfo(instance.graph);
        instance.graph.printAnalysis();
        instance.serialize(getParameter("Output file"));
    }
};

class RunPushRelabel : public ParameterizedCommand {

public:
    RunPushRelabel(BasicShell& shell) :
        ParameterizedCommand(shell, "runPushRelabel", "Computes a minimum s-t-cut on the given graph with push-relabel.") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        StaticInstance instance(getParameter("Instance file"));
        PushRelabel<StaticInstance> algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "Flow value: " << algorithm.getFlowValue() << std::endl;
    }
};

class TestParametricPushRelabel : public ParameterizedCommand {

public:
    TestParametricPushRelabel(BasicShell& shell) :
    ParameterizedCommand(shell, "testParametricPushRelabel", "Compares restartable push-relabel to regular push-relabel on the given graph.") {
        addParameter("Instance file");
        addParameter("Steps");
    }

    virtual void execute() noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        ParametricWrapper wrapper(instance);
        PushRelabel<ParametricWrapper> algorithm(wrapper);
        run(algorithm, false);

        const int steps = getParameter<int>("Steps");
        for (int i = 1; i <= steps; i++) {
            const double alpha = instance.alphaMin + i * (instance.alphaMax - instance.alphaMin) / steps;
            std::cout << "Alpha: " << alpha << std::endl;
            wrapper.setAlpha(alpha);
            run(algorithm, true);
            PushRelabel<ParametricWrapper> newAlgorithm(wrapper);
            run(newAlgorithm, false);
        }
    }

private:
    inline void run(PushRelabel<ParametricWrapper>& algorithm, const bool update) const noexcept {
        Timer timer;
        if (update) {
            algorithm.continueAfterUpdate();
        } else {
            algorithm.run();
        }
        std::cout << "\tTime: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "\t#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "\t#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "\tFlow value: " << algorithm.getFlowValue() << std::endl;
    }
};

class RunIBFS : public ParameterizedCommand {

public:
    RunIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "runIBFS", "Computes a minimum s-t-cut on the given graph with IBFS.") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        StaticInstance instance(getParameter("Instance file"));
        IBFS<StaticInstance> algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "Flow value: " << algorithm.getFlowValue() << std::endl;
    }
};

class RunExcessesIBFS : public ParameterizedCommand {

public:
    RunExcessesIBFS(BasicShell& shell) :
    ParameterizedCommand(shell, "runExcessesIBFS", "Computes a minimum s-t-cut on the given graph with Excesses IBFS.") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        StaticInstance instance(getParameter("Instance file"));
        ExcessesIBFS<StaticInstance> algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "Flow value: " << algorithm.getFlowValue() << std::endl;
    }
};

class TestRestartableIBFS : public ParameterizedCommand {

public:
    TestRestartableIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "testRestartableIBFS", "Compares restartable IBFS to regular IBFS on the given graph.") {
        addParameter("Instance file");
        addParameter("Steps");
    }

    virtual void execute() noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        ParametricWrapper wrapper(instance);
        RestartableIBFS<ParametricWrapper> algorithm(wrapper);
        run(algorithm, false);

        const int steps = getParameter<int>("Steps");
        for (int i = 1; i <= steps; i++) {
            const double alpha = instance.alphaMin + i * (instance.alphaMax - instance.alphaMin) / steps;
            std::cout << "Alpha: " << alpha << std::endl;
            wrapper.setAlpha(alpha);
            run(algorithm, true);
            RestartableIBFS<ParametricWrapper> newAlgorithm(wrapper);
            run(newAlgorithm, false);
        }
    }

private:
    inline void run(RestartableIBFS<ParametricWrapper>& algorithm, const bool update) const noexcept {
        Timer timer;
        if (update) {
            algorithm.continueAfterUpdate();
        } else {
            algorithm.run();
        }
        std::cout << "\tTime: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "\t#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "\t#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "\tFlow value: " << algorithm.getFlowValue() << std::endl;
    }
};

class RunParametricIBFS : public ParameterizedCommand {

public:
    RunParametricIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "runParametricIBFS", "Computes a parametric minimum s-t-cut on the given graph with Parametric IBFS.") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        ParametricIBFS<ParametricInstance::FlowFunction> algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Breakpoints: " << algorithm.getBreakpoints().size() << std::endl;
    }
};


inline void compareSinkComponents(const double breakpoint, const std::vector<Vertex>& a, const std::vector<Vertex>& b, const size_t n) noexcept {
    bool header = false;
    std::vector<bool> inA(n, false);
    std::vector<bool> inB(n, false);
    for (const Vertex v : a) {
        inA[v] = true;
    }
    for (const Vertex v : b) {
        inB[v] = true;
    }
    for (const Vertex v : a) {
        if (inB[v]) continue;
        if (!header) {
            std::cout << "Breakpoint " << breakpoint << ":" << std::endl;
            header = true;
        }
        std::cout << "\tVertex " << v << " is in parametric but not restartable" << std::endl;
    }
    for (const Vertex v : b) {
        if (inA[v]) continue;
        if (!header) {
            std::cout << "Breakpoint " << breakpoint << ":" << std::endl;
            header = true;
        }
        std::cout << "Vertex " << v << " is in restartable but not parametric" << std::endl;
    }
}

inline bool areFlowValuesEqual(const double a, const double b) noexcept {
    return abs(a - b) <= 0.01;
}

class TestParametricIBFS : public ParameterizedCommand {

public:
    TestParametricIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "testParametricIBFS", "Compares Parametric IBFS to a restartable algorithm (Push-Relabel or IBFS) on the given graph.") {
        addParameter("Instance file");
        addParameter("Restartable algorithm", {"Push-Relabel", "IBFS"});
    }

    virtual void execute() noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        ParametricIBFS<ParametricInstance::FlowFunction> algorithm(instance);
        Timer timer;
        algorithm.run();
        const std::vector<double>& breakpoints = algorithm.getBreakpoints();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Breakpoints: " << breakpoints.size() << std::endl;
        if (getParameter("Restartable algorithm") == "Push-Relabel") {
            compare<PushRelabel<ParametricWrapper>>(instance, algorithm, breakpoints);
        } else {
            compare<RestartableIBFS<ParametricWrapper>>(instance, algorithm, breakpoints);
        }

    }

private:
    template<typename RESTARTABLE_ALGORITHM>
    inline void compare(const ParametricInstance& instance, const ParametricIBFS<ParametricInstance::FlowFunction>& parametricAlgorithm, const std::vector<double>& breakpoints) const noexcept {
        ParametricWrapper wrapper(instance);
        RESTARTABLE_ALGORITHM restartableAlgorithm(wrapper);
        Progress progress(breakpoints.size());
        Timer timer;
        double restartableTime = 0;
        for (size_t i = 0; i < breakpoints.size(); i++) {
            timer.restart();
            if (i == 0) {
                restartableAlgorithm.run();
            } else {
                wrapper.setAlpha(breakpoints[i]);
                restartableAlgorithm.continueAfterUpdate();
            }
            restartableTime += timer.elapsedMicroseconds();
            if (!areFlowValuesEqual(parametricAlgorithm.getFlowValue(breakpoints[i]), restartableAlgorithm.getFlowValue())) {
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
                std::cout << "Flow values for breakpoint " << breakpoints[i] << " are not equal! ";
                std::cout << "Parametric: " << parametricAlgorithm.getFlowValue(breakpoints[i]) << ", restartable: " << restartableAlgorithm.getFlowValue() << std::endl;
            }
            compareSinkComponents(breakpoints[i], parametricAlgorithm.getSinkComponent(breakpoints[i]), restartableAlgorithm.getSinkComponent(), instance.graph.numVertices());
            progress++;
        }
        std::cout << "Restartable time: " << String::musToString(restartableTime) << std::endl;
    }
};

class TestChordScheme : public ParameterizedCommand {

public:
    TestChordScheme(BasicShell& shell) :
        ParameterizedCommand(shell, "testChordScheme", "Solves the given parametric max-flow instance with the chord scheme.") {
        addParameter("Instance file");
        addParameter("Precision");
        addParameter("Chord scheme algorithm", {"Push-Relabel", "IBFS"});
        addParameter("Restartable algorithm", {"Push-Relabel", "IBFS"});
    }

    using ChordSchemeWrapper = ChordSchemeMaxFlowWrapper<pmf::linearFlowFunction>;

    virtual void execute() noexcept {
        if (getParameter("Chord scheme algorithm") == "Push-Relabel") {
            run<PushRelabel<ChordSchemeWrapper>>();
        } else {
            run<IBFS<ChordSchemeWrapper>>();
        }
    }

private:
    template<typename SEARCH_ALGORITHM>
    inline void run() const noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        const double precision = std::pow(10, -getParameter<int>("Precision"));
        ChordScheme<pmf::linearFlowFunction, SEARCH_ALGORITHM> chordScheme(instance, precision);
        Timer timer;
        chordScheme.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "Solutions: " << chordScheme.getSolutions().size() << std::endl;
        if (getParameter("Restartable algorithm") == "Push-Relabel") {
            compare<PushRelabel<ParametricWrapper>>(instance, chordScheme);
        } else {
            compare<RestartableIBFS<ParametricWrapper>>(instance, chordScheme);
        }
    }

    template<typename RESTARTABLE_ALGORITHM, typename CHORD_SCHEME>
    inline void compare(const ParametricInstance& instance, const CHORD_SCHEME& chordScheme) const noexcept {
        using Solution = CHORD_SCHEME::Solution;
        const std::vector<Solution>& solutions = chordScheme.getSolutions();
        ParametricWrapper wrapper(instance);
        RESTARTABLE_ALGORITHM restartableAlgorithm(wrapper);
        Progress progress(solutions.size());
        Timer timer;
        double restartableTime = 0;
        for (size_t i = 0; i < solutions.size(); i++) {
            const double alpha = solutions[i].breakpoint;
            timer.restart();
            if (i == 0) {
                restartableAlgorithm.run();
            } else {
                wrapper.setAlpha(alpha);
                restartableAlgorithm.continueAfterUpdate();
            }
            restartableTime += timer.elapsedMicroseconds();
            const double flowValue = solutions[i].getFlowValue();
            if (!areFlowValuesEqual(flowValue, restartableAlgorithm.getFlowValue())) {
                std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
                std::cout << "Flow values for breakpoint " << alpha << " are not equal! ";
                std::cout << "Parametric: " << flowValue << ", restartable: " << restartableAlgorithm.getFlowValue() << std::endl;
            }
            compareSinkComponents(alpha, chordScheme.getSinkComponent(alpha), restartableAlgorithm.getSinkComponent(), instance.graph.numVertices());
            progress++;
        }
        std::cout << "Restartable time: " << String::musToString(restartableTime) << std::endl;
    }
};