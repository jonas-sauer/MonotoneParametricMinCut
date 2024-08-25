#pragma once

#include <vector>
#include <type_traits>

#include "../StaticMinCut/PushRelabel.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/MaxFlow/FlowFunction.h"
#include "../../DataStructures/MaxFlow/MaxFlowInstance.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Meta.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

template<Meta::Derived<FlowFunction> FLOW_FUNCTION, typename SEARCH_ALGORITHM, bool MEASUREMENTS = false>
class DichotomicScheme {
public:
    using FlowFunction = FLOW_FUNCTION;
    using FlowType = FlowFunction::FlowType;
    using ParametricInstance = ParametricMaxFlowInstance<FlowFunction>;
    using ParametricWrapper = DichotomicSchemeMaxFlowWrapper<FlowFunction>;
    using SearchAlgorithm = SEARCH_ALGORITHM;

    struct NoMeasurements {
        inline void startTimer() noexcept {}
        inline void measureContractionTime() noexcept {}
        inline void measureFlowTime() noexcept {}
        inline void addToTotalVertices(const size_t) noexcept {}
        inline void print() const noexcept {}
        [[nodiscard]] inline std::string getCSV() const noexcept { return ""; }
    };

    struct Measurements {
        inline void startTimer() noexcept {
            timer.restart();
        }

        inline void measureContractionTime() noexcept {
            contractionTime += timer.elapsedMicroseconds();
        }

        inline void measureFlowTime() noexcept {
            flowTime += timer.elapsedMicroseconds();
        }

        inline void addToTotalVertices(const size_t value) noexcept {
            totalVertices += value;
        }

        inline void print() const noexcept {
            std::cout << "Contraction time: " << String::musToString(contractionTime) << std::endl;
            std::cout << "Flow time: " << String::musToString(flowTime) << std::endl;
            std::cout << "#Vertices (total): " << totalVertices << std::endl;
        }

        [[nodiscard]] inline std::string getCSV() const noexcept {
            return std::to_string(contractionTime) + "," +
                   std::to_string(flowTime) + "," +
                   std::to_string(totalVertices);
        }

        double contractionTime = 0;
        double flowTime = 0;
        long long totalVertices = 0;
        Timer timer;
    };

    using MeasurementsType = std::conditional_t<MEASUREMENTS, Measurements, NoMeasurements>;

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

    DichotomicScheme(const ParametricInstance& instance, const double epsilon) : instance(instance), epsilon(epsilon), breakpointOfVertex(instance.graph.numVertices(), INFTY) {}

    inline void run() noexcept {
        ParametricWrapper wrapper(instance);
        const Solution solMin = runSearch(wrapper, instance.alphaMin);
        const Solution solMax = runSearch(wrapper, instance.alphaMax);
        for (const Vertex vertex : instance.graph.vertices()) {
            if (!solMin.inSinkComponent[vertex]) {
                breakpointOfVertex[vertex] = instance.alphaMin;
            }
        }
        breakpoints.emplace_back(instance.alphaMin);
        measurements.startTimer();
        ParametricWrapper contractedWrapper = wrapper.contractSourceAndSinkComponents(solMin.inSinkComponent, solMax.inSinkComponent);
        measurements.measureContractionTime();
        recurse(instance.alphaMin, instance.alphaMax, solMin, solMax, contractedWrapper, wrapper);
        if (instance.alphaMax < INFTY) addSolution(instance.alphaMax, solMax, wrapper);
        measurements.print();
    }

    [[nodiscard]] inline const std::vector<double>& getBreakpoints() const noexcept {
        return breakpoints;
    }

    [[nodiscard]] inline const std::vector<double>& getVertexBreakpoints() const noexcept {
        return breakpointOfVertex;
    }

    [[nodiscard]] inline std::vector<Vertex> getSinkComponent(const double alpha) const noexcept {
        std::vector<Vertex> sinkComponent;
        for (const Vertex vertex : instance.graph.vertices()) {
            if (breakpointOfVertex[vertex] <= alpha) continue;
            sinkComponent.emplace_back(vertex);
        }
        return sinkComponent;
    }

    [[nodiscard]] inline double getFlowValue(const double alpha) const noexcept {
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

    inline std::string getMeasurementsCSV() const noexcept {
        return measurements.getCSV();
    }

private:
    inline Solution runSearch(ParametricWrapper& wrapper, const double alpha) noexcept {
        measurements.startTimer();
        measurements.addToTotalVertices(wrapper.graph.numVertices());
        wrapper.setAlpha(alpha);
        SearchAlgorithm search(wrapper);
        search.run();
        Solution result(wrapper, search, alpha);
        measurements.measureFlowTime();
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

        measurements.startTimer();
        ParametricWrapper wrapperLeft = wrapper.contractSinkComponent(solMid.inSinkComponent);
        ParametricWrapper wrapperRight = wrapper.contractSourceComponent(solMid.inSinkComponent);
        measurements.measureContractionTime();
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
    MeasurementsType measurements;
};