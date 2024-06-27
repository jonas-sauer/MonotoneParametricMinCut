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
#include "../../Algorithms/MaxFlowMinCut/PushRelabel.h"
#include "../../Algorithms/MaxFlowMinCut/IBFS.h"
#include "../../Algorithms/MaxFlowMinCut/ExcessesIBFS.h"
#include "../../Algorithms/MaxFlowMinCut/RestartableIBFS.h"
#include "../../Algorithms/MaxFlowMinCut/ParametricIBFS.h"

using namespace Shell;

using StaticInstance = StaticMaxFlowInstance<int>;
using ParametricInstance = ParametricMaxFlowInstance<pmf::linearFlowFunction>;
using ParametricWrapper = ParametricToStaticMaxFlowInstanceWrapper<pmf::linearFlowFunction>;

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

//TODO: Do this properly once we have true parametric instances
class TestParametricPushRelabel : public ParameterizedCommand {

public:
    TestParametricPushRelabel(BasicShell& shell) :
    ParameterizedCommand(shell, "testParametricPushRelabel", "temp") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        StaticInstance staticInstance(getParameter("Instance file"));
        ParametricInstance instance(staticInstance, 10);
        ParametricWrapper wrapper(instance, 0);
        PushRelabel<ParametricWrapper> algorithm(wrapper);
        run(algorithm, false);

        for (int i = 1; i <= instance.alphaMax; i++) {
            std::cout << "Instance " << i << std::endl;
            wrapper.setAlpha(i);
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

//TODO: Do this properly once we have true parametric instances
class TestRestartableIBFS : public ParameterizedCommand {

public:
    TestRestartableIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "testRestartableIBFS", "temp") {
        addParameter("Instance file");
    }

    virtual void execute() noexcept {
        StaticInstance staticInstance(getParameter("Instance file"));
        ParametricInstance instance(staticInstance, 10);
        ParametricWrapper wrapper(instance, 0);
        RestartableIBFS<ParametricWrapper> algorithm(wrapper);
        run(algorithm, false);

        for (int i = 1; i <= instance.alphaMax; i++) {
            std::cout << "Instance " << i << std::endl;
            wrapper.setAlpha(i);
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
        ParametricIBFS algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Breakpoints: " << algorithm.getBreakpoints() << std::endl;
    }
};