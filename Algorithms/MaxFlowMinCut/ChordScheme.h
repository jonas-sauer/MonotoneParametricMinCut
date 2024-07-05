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
    using FlowGraph = ParametricInstance::GraphType;

    struct Solution {
        Solution(const ParametricWrapper& wrapper, const SearchAlgorithm& search, const double alpha) :
            breakpoint(alpha),
            sinkComponent(wrapper.translate(search.getSinkComponent())),
            inSinkComponent(search.getInSinkComponent()),
            flowFunction(wrapper.getCapacity(search.getCutEdges())),
            flowValue(flowFunction.eval(alpha)) {
        }

        double breakpoint;
        std::vector<Vertex> sinkComponent;
        std::vector<bool> inSinkComponent;
        FlowFunction flowFunction;
        double flowValue;

        inline bool operator<(const Solution& other) const noexcept {
            return breakpoint < other.breakpoint;
        }
    };

    ChordScheme(const ParametricInstance& instance, const double epsilon) : instance(instance), epsilon(epsilon) {}

    inline void run() noexcept {
        ParametricWrapper wrapper(instance);
        solutions.clear();
        const Solution solMin = runSearch(wrapper, instance.alphaMin);
        solutions.push_back(solMin);
        const Solution solMax = runSearch(wrapper, instance.alphaMax);
        solutions.push_back(solMax);
        ParametricWrapper contractedWrapper = wrapper.contractSourceAndSinkComponents(solMin.inSinkComponent, solMax.inSinkComponent);
        recurse(instance.alphaMin, instance.alphaMax, solMin.flowFunction, solMax.flowFunction, contractedWrapper);
        std::sort(solutions.begin(), solutions.end());
    }

    [[nodiscard]] const std::vector<Solution>& getSolutions() const noexcept {
        return solutions;
    }

private:
    inline Solution runSearch(ParametricWrapper& wrapper, const double alpha) noexcept {
        wrapper.setAlpha(alpha);
        SearchAlgorithm search(wrapper);
        search.run();
        return Solution(wrapper, search, alpha);
    }

    inline void recurse(const double left, const double right, const FlowFunction& flowLeft, const FlowFunction& flowRight, ParametricWrapper& wrapper) noexcept {
        if (right <= left) return;
        const double mid = findIntersectionPoint(flowLeft, flowRight);
        assert(mid >= left && mid <= right);

        const Solution& solMid = runSearch(wrapper, mid);
        const double oldVal = flowLeft.eval(mid);
        const double newVal = solMid.flowFunction.eval(mid);
        if (oldVal <= (1 + epsilon) * newVal) return;

        solutions.emplace_back(solMid);
        ParametricWrapper wrapperLeft = wrapper.contractSinkComponent(solMid.inSinkComponent);
        ParametricWrapper wrapperRight = wrapper.contractSourceComponent(solMid.inSinkComponent);
        recurse(left, mid, flowLeft, solMid.flowFunction, wrapperLeft);
        recurse(mid, right, solMid.flowFunction, flowRight, wrapperRight);
    }

private:
    const ParametricInstance& instance;
    const double epsilon;
    std::vector<Solution> solutions;
};