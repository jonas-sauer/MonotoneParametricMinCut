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

class RunParametricIBFS : public ParameterizedCommand {

public:
    RunParametricIBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "runParametricIBFS", "Computes a parametric minimum s-t-cut on the given graph with Parametric IBFS.") {
        addParameter("Instance file");
        addParameter("With measurements?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("With measurements?")) {
            run<true>();
        } else {
            run<false>();
        }
    }

private:
    template<bool MEASUREMENTS>
    inline void run() noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        ParametricIBFS<ParametricInstance::FlowFunction, MEASUREMENTS> algorithm(instance);
        Timer timer;
        algorithm.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "#Breakpoints: " << algorithm.getBreakpoints().size() << std::endl;
    }
};

class RunChordScheme : public ParameterizedCommand {

public:
    RunChordScheme(BasicShell& shell) :
        ParameterizedCommand(shell, "runChordScheme", "Computes a parametric minimum s-t-cut on the given graph with the chord scheme.") {
        addParameter("Instance file");
        addParameter("Precision");
        addParameter("Flow algorithm", {"Push-Relabel", "IBFS"});
    }

    using ChordSchemeWrapper = ChordSchemeMaxFlowWrapper<pmf::linearFlowFunction>;

    virtual void execute() noexcept {
        if (getParameter("Flow algorithm") == "Push-Relabel") {
            run<PushRelabel<ChordSchemeWrapper>>();
        } else {
            run<IBFS<ChordSchemeWrapper>>();
        }
    }

private:
    template<typename SEARCH_ALGORITHM>
    inline void run() const noexcept {
        ParametricInstance instance(getParameter("Instance file"));
        const int exponent = getParameter<int>("Precision");
        const double precision = (exponent < 0) ? 0 : std::pow(10, -exponent);
        ChordScheme<pmf::linearFlowFunction, SEARCH_ALGORITHM, true> chordScheme(instance, precision);
        Timer timer;
        chordScheme.run();
        std::cout << "Time: " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
        std::cout << "Solutions: " << chordScheme.getBreakpoints().size() << std::endl;
    }
};

class PrecisionExperiment : public ParameterizedCommand {

public:
    PrecisionExperiment(BasicShell& shell) :
        ParameterizedCommand(shell, "precisionExperiment", "Compares the solution of Parametric IBFS and chord scheme.") {
        addParameter("Instance file");
        addParameter("Chord scheme algorithm", {"Push-Relabel", "IBFS"});
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
        ParametricIBFS<ParametricInstance::FlowFunction, true> parametricIBFS(instance);
        ChordScheme<pmf::linearFlowFunction, SEARCH_ALGORITHM, true> chordScheme(instance, 0);
        parametricIBFS.run();
        chordScheme.run();

        const std::vector<double>& parametricBreakpoints = parametricIBFS.getBreakpoints();
        const std::vector<double>& chordBreakpoints = chordScheme.getBreakpoints();
        std::cout << "Parametric IBFS: " << parametricBreakpoints.size() << " breakpoints" << std::endl;
        std::cout << "Chord scheme: " << chordBreakpoints.size() << " breakpoints" << std::endl;

        std::cout << "Evaluate chord scheme:" << std::endl;
        compare(parametricIBFS, chordScheme);
        std::cout << "Evaluate parametric IBFS:" << std::endl;
        compare(chordScheme, parametricIBFS);
    }

    template<typename TRUTH_ALGO, typename COMP_ALGO>
    inline void compare(const TRUTH_ALGO& truthAlgo, const COMP_ALGO& compAlgo) const noexcept {
        const std::vector<double>& groundTruth = truthAlgo.getBreakpoints();
        Progress progress(groundTruth.size());
        double cumulativeError = 0;
        size_t numErrors = 0;
        for (const double breakpoint : groundTruth) {
            const double actualFlow = truthAlgo.getFlowValue(breakpoint);
            const double resultFlow = compAlgo.getFlowValue(breakpoint);
            progress++;
            if (resultFlow <= actualFlow + 1e-06) continue;
            std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << actualFlow << " vs. " << resultFlow << " ( " << resultFlow - actualFlow << ")" << std::endl;
            cumulativeError += (resultFlow - actualFlow)/actualFlow;
            numErrors++;
        }
        progress.finished();
        std::cout << "Errors: " << numErrors << "/" << groundTruth.size() << std::endl;
        std::cout << "Cumulative error: " << cumulativeError << std::endl;
        std::cout << "Average error: " << (numErrors == 0 ? 0 : cumulativeError/numErrors) << std::endl;
        std::cout << "Accuracy: " << cumulativeError/groundTruth.size() << std::endl;
    }
};