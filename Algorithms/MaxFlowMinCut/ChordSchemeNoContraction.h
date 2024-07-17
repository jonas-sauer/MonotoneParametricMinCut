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

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION, typename SEARCH_ALGORITHM>
class ChordSchemeNoContraction {
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

    ChordSchemeNoContraction(const ParametricInstance& instance, const double epsilon) : instance(instance), wrapper(instance), epsilon(epsilon), breakpointOfVertex(instance.graph.numVertices(), INFTY) {}

    inline void run() noexcept {
        const Solution solMin = runSearch(instance.alphaMin);
        const Solution solMax = runSearch(instance.alphaMax);
        for (const Vertex vertex : instance.graph.vertices()) {
            if (!solMin.inSinkComponent[vertex]) {
                breakpointOfVertex[vertex] = instance.alphaMin;
            }
        }
        breakpoints.emplace_back(instance.alphaMin);
        recurse(instance.alphaMin, instance.alphaMax, solMin, solMax);
        if (instance.alphaMax < INFTY) addSolution(instance.alphaMax, solMax);
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

private:
    inline Solution runSearch(const double alpha) noexcept {
        wrapper.setAlpha(alpha);
        SearchAlgorithm search(wrapper);
        search.run();
        Solution result(wrapper, search, alpha);
        return result;
    }

    inline void recurse(const double left, const double right, const Solution& solLeft, const Solution& solRight) noexcept {
        if (right <= left) {
            addSolution(left, solRight);
            return;
        }
        const double mid = findIntersectionPoint(solLeft.flowFunction, solRight.flowFunction);
        if (mid <= left || mid >= right) {
            addSolution(left, solRight);
            return;
        }

        const Solution solMid = runSearch(mid);
        const double oldVal = solLeft.flowFunction.eval(mid);
        const double newVal = solMid.flowFunction.eval(mid);
        if (oldVal <= (1 + epsilon) * newVal) {
            addSolution(mid, solRight);
            return;
        }

        recurse(left, mid, solLeft, solMid);
        recurse(mid, right, solMid, solRight);
    }

    inline void addSolution(const double breakpoint, const Solution& solution) noexcept {
        bool addSolution = false;
        for (const Vertex vertex : wrapper.graph.vertices()) {
            if (solution.inSinkComponent[vertex] || vertex == wrapper.source) continue;
            if (breakpoint < breakpointOfVertex[vertex]) {
                breakpointOfVertex[vertex] = breakpoint;
                addSolution = true;
            }
        }
        if (addSolution) breakpoints.emplace_back(breakpoint);
    }

private:
    const ParametricInstance& instance;
    ParametricWrapper wrapper;
    const double epsilon;
    std::vector<double> breakpoints;
    std::vector<double> breakpointOfVertex;
};