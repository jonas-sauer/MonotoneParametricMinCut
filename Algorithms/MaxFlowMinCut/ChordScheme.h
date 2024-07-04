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

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class ChordScheme {
public:
    using FlowFunction = FLOW_FUNCTION;
    using ParametricInstance = ParametricMaxFlowInstance<FlowFunction>;
    using ParametricWrapper = ChordSchemeMaxFlowWrapper<FlowFunction>;
    using SearchAlgorithm = PushRelabel<ParametricWrapper>;

    struct Solution {
        double breakpoint;
        std::vector<Vertex> sinkComponent;
        FlowFunction flowFunction;
        double flowValue;

        inline bool operator<(const Solution& other) const noexcept {
            return breakpoint < other.breakpoint;
        }
    };

    ChordScheme(const ParametricInstance& instance, const double epsilon) : instance(instance), wrapper(instance), epsilon(epsilon) {}

    inline void run() noexcept {
        solutions.clear();
        const Solution solMin = runSearch(instance.alphaMin);
        solutions.push_back(solMin);
        const Solution solMax = runSearch(instance.alphaMax);
        solutions.push_back(solMax);
        recurse(instance.alphaMin, instance.alphaMax, solMin.flowFunction, solMax.flowFunction);
        std::sort(solutions.begin(), solutions.end());
    }

    [[nodiscard]] const std::vector<Solution>& getSolutions() const noexcept {
        return solutions;
    }

private:
    inline Solution runSearch(const double alpha) noexcept {
        SearchAlgorithm search(wrapper);
        wrapper.setAlpha(alpha);
        search.run();
        const FlowFunction flowFunction = wrapper.getCapacity(search.getCutEdges());
        return Solution{alpha, search.getSinkComponent(), flowFunction, flowFunction.eval(alpha)};
    }

    inline void recurse(const double left, const double right, const FlowFunction& flowLeft, const FlowFunction& flowRight) noexcept {
        if (right <= left) return;
        const double mid = findIntersectionPoint(flowLeft, flowRight);
        assert(mid > left && mid < right);

        const Solution& solMid = runSearch(mid);
        const double oldVal = flowLeft.eval(mid);
        const double newVal = solMid.flowFunction.eval(mid);
        if (oldVal <= (1 + epsilon) * newVal) return;

        solutions.emplace_back(solMid);
        recurse(left, mid, flowLeft, solMid.flowFunction);
        recurse(mid, right, solMid.flowFunction, flowRight);
    }

private:
    const ParametricInstance& instance;
    ParametricWrapper wrapper;
    const double epsilon;
    std::vector<Solution> solutions;
};