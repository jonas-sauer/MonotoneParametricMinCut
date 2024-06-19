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

using namespace Shell;

class LoadMaxFlowInstanceFromDimacs : public ParameterizedCommand {

public:
    LoadMaxFlowInstanceFromDimacs(BasicShell& shell) :
        ParameterizedCommand(shell, "loadMaxFlowInstanceFromDimacs", "Load the given max-flow instance in DIMACS format.") {
        addParameter("Input file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        MaxFlowInstance instance;
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
        MaxFlowInstance instance(getParameter("Instance file"));
        PushRelabel algorithm(instance);
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
        MaxFlowInstance instance(getParameter("Instance file"));
        instance.maxParameter = 2;
        Edge edgeFromSource = instance.graph.beginEdgeFrom(instance.source);
        for (size_t i = 0; i < instance.graph.outDegree(instance.source); edgeFromSource++, i++) {
            instance.sourceEdgeSlopes[i] = instance.getCapacity(edgeFromSource) / instance.maxParameter;
            instance.graph.set(Capacity, edgeFromSource, 0);
            instance.currentCapacity[edgeFromSource] = 0;
        }
        Edge edgeFromSink = instance.graph.beginEdgeFrom(instance.sink);
        for (size_t i = 0; i < instance.graph.outDegree(instance.sink); edgeFromSink++, i++) {
            const Edge edgeToSink = instance.graph.get(ReverseEdge, edgeFromSink);
            instance.sinkEdgeSlopes[i] = -instance.getCapacity(edgeToSink) / instance.maxParameter;
        }

        PushRelabel algorithm(instance);
        run(algorithm, false);

        for (int i = 1; i <= instance.maxParameter; i++) {
            instance.setCurrentParameter(i);
            run(algorithm, true);
            PushRelabel newAlgorithm(instance);
            run(newAlgorithm, false);
        }
    }

private:
    inline void run(PushRelabel& algorithm, const bool update) const noexcept {
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
        MaxFlowInstance instance(getParameter("Instance file"));
        IBFS algorithm(instance);
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
        MaxFlowInstance instance(getParameter("Instance file"));
        ExcessesIBFS algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Source component: " << algorithm.getSourceComponent().size() << std::endl;
        std::cout << "#Sink component: " << algorithm.getSinkComponent().size() << std::endl;
        std::cout << "Flow value: " << algorithm.getFlowValue() << std::endl;
    }
};