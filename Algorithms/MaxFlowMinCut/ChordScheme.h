#pragma once

#include <vector>

#include "PushRelabel.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlowMinCut/FlowUtils.h"
#include "../../DataStructures/MaxFlowMinCut/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Meta.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION, typename SEARCH_ALGORITHM, bool MEASUREMENTS = false>
class ChordScheme {
public:
    using FlowFunction = FLOW_FUNCTION;
    using FlowType = FlowFunction::FlowType;
    using ParametricInstance = ParametricMaxFlowInstance<FlowFunction>;
    using ParametricWrapper = ChordSchemeMaxFlowWrapper<FlowFunction>;
    using SearchAlgorithm = SEARCH_ALGORITHM;

    struct Solution {
        Solution(const ParametricWrapper& wrapper, const SearchAlgorithm& search, const double alpha) :
            breakpoint(alpha),
            inSinkComponent(search.getInSinkComponent()),
            flowFunction(wrapper.getCapacity(search.getCutEdges())) {
        }

        double breakpoint;
        std::vector<bool> inSinkComponent;
        FlowFunction flowFunction;

        inline FlowType getFlowValue() const noexcept {
            return flowFunction.eval(breakpoint);
        }

        inline bool operator<(const Solution& other) const noexcept {
            return breakpoint < other.breakpoint;
        }
    };

    ChordScheme(const ParametricInstance& instance, const double epsilon) : instance(instance), epsilon(epsilon), breakpointOfVertex(instance.graph.numVertices(), INFTY) {}

    inline void run() noexcept {
        ParametricWrapper wrapper(instance);
        const Solution solMin = runSearch(wrapper, instance.alphaMin);
        addSolution(instance.alphaMin, solMin, wrapper);
        const Solution solMax = runSearch(wrapper, instance.alphaMax);
        if (instance.alphaMax < INFTY) addSolution(instance.alphaMax, solMax, wrapper);
        for (const Vertex vertex : instance.graph.vertices()) {
            if (!solMin.inSinkComponent[vertex]) {
                breakpointOfVertex[vertex] = instance.alphaMin;
            }
        }
        timer.restart();
        ParametricWrapper contractedWrapper = wrapper.contractSourceAndSinkComponents(solMin.inSinkComponent, solMax.inSinkComponent);
        if constexpr (MEASUREMENTS) contractionTime += timer.elapsedMicroseconds();
        recurse(instance.alphaMin, instance.alphaMax, solMin, solMax, contractedWrapper, wrapper);
        if constexpr (MEASUREMENTS) {
            std::cout << "Contraction time: " << String::musToString(contractionTime) << std::endl;
            std::cout << "Flow time: " << String::musToString(flowTime) << std::endl;
            std::cout << "#Bad splits: " << badSplits << std::endl;
        }
    }

    inline const std::vector<double>& getBreakpoints() const noexcept {
        return breakpoints;
    }

    inline const std::vector<double>& getVertexBreakpoints() const noexcept {
        return breakpointOfVertex;
    }

    inline std::vector<Vertex> getSinkComponent(const double alpha) const noexcept {
        std::vector<Vertex> sinkComponent;
        for (const Vertex vertex : instance.graph.vertices()) {
            if (breakpointOfVertex[vertex] <= alpha) continue;
            sinkComponent.emplace_back(vertex);
        }
        return sinkComponent;
    }

    inline double getFlowValue(const double alpha) const noexcept {
        double flow = 0;
        for (const Vertex from : instance.graph.vertices()) {
            if (breakpointOfVertex[from] > alpha) continue;
            for (const Edge edge : instance.graph.edgesFrom(from)) {
                const Vertex to = instance.graph.get(ToVertex, edge);
                if (breakpointOfVertex[to] <= alpha) continue;
                flow += instance.getCapacity(edge, alpha);
            }
        }
        return flow;
    }

    inline double getContractionTime() const noexcept {
        // if (!MEASUREMENTS) throw std::runtime_error("Detailed measurements are only done if template parameter MEASUREMENTS is true");
        return contractionTime;
    }

    inline double getFlowTime() const noexcept {
        // if (!MEASUREMENTS) throw std::runtime_error("Detailed measurements are only done if template parameter MEASUREMENTS is true");
        return flowTime;
    }

    inline size_t getNumBadSplits() const noexcept {
        // if (!MEASUREMENTS) throw std::runtime_error("Detailed measurements are only done if template parameter MEASUREMENTS is true");
        return badSplits;
    }

private:
    inline Solution runSearch(ParametricWrapper& wrapper, const double alpha) noexcept {
        if constexpr (MEASUREMENTS) timer.restart();
        wrapper.setAlpha(alpha);
        SearchAlgorithm search(wrapper);
        search.run();
        Solution result(wrapper, search, alpha);
        if constexpr (MEASUREMENTS) flowTime += timer.elapsedMicroseconds();
        return result;
    }

    inline void recurse(const double left, const double right, const Solution& solLeft, const Solution& solRight, ParametricWrapper& wrapper, const ParametricWrapper& evalWrapper) noexcept {
        if (right <= left) {
            addSolution(left, solRight, evalWrapper);
            return;
        }
        const double mid = findIntersectionPoint(solLeft.flowFunction, solRight.flowFunction);
        if (mid <= left || mid >= right) {
            addSolution(left, solRight, evalWrapper);
            return;
        }

        const Solution solMid = runSearch(wrapper, mid);
        const double oldVal = solLeft.flowFunction.eval(mid);
        const double newVal = solMid.flowFunction.eval(mid);
        if (oldVal <= (1 + epsilon) * newVal) {
            addSolution(mid, solRight, evalWrapper);
            return;
        }

        if constexpr (MEASUREMENTS) timer.restart();
        ParametricWrapper wrapperLeft = wrapper.contractSinkComponent(solMid.inSinkComponent);
        ParametricWrapper wrapperRight = wrapper.contractSourceComponent(solMid.inSinkComponent);
        if constexpr (MEASUREMENTS) contractionTime += timer.elapsedMicroseconds();
        if constexpr (MEASUREMENTS) {
            const size_t small = std::min(wrapperLeft.graph.numVertices(), wrapperRight.graph.numVertices());
            const size_t large = std::max(wrapperLeft.graph.numVertices(), wrapperRight.graph.numVertices());
            if (large >= 9 * small) badSplits++;
        }
        recurse(left, mid, solLeft, solMid, wrapperLeft, wrapper);
        recurse(mid, right, solMid, solRight, wrapperRight, evalWrapper);
    }

    inline void addSolution(const double breakpoint, const Solution& solution, const ParametricWrapper& wrapper) noexcept {
        bool addSolution = false;
        for (const Vertex vertex : wrapper.graph.vertices()) {
            if (solution.inSinkComponent[vertex] || vertex == wrapper.source) continue;
            const Vertex realVertex = wrapper.newToOldVertex[vertex];
            if (breakpoint < breakpointOfVertex[realVertex]) {
                breakpointOfVertex[realVertex] = breakpoint;
                addSolution = true;
            }
        }
        if (addSolution) breakpoints.emplace_back(breakpoint);
    }

private:
    const ParametricInstance& instance;
    const double epsilon;
    std::vector<double> breakpoints;
    std::vector<double> breakpointOfVertex;

    double contractionTime = 0;
    double flowTime = 0;
    size_t badSplits = 0;
    Timer timer;
};