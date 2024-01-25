#pragma once

#include <cmath>
#include <vector>

#include "../../DataStructures/Assignment/Settings.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"

namespace DecisionModels {

struct NodeData {
    PerceivedTime expectedPAT{Unreachable};
    int arrivalTime{INFTY};
};

class DistanceLogit {
public:
    DistanceLogit(const Assignment::Settings& settings) :
        beta(settings.distanceBeta),
        referenceTravelTime(settings.referenceTravelTime),
        probabilityFactor(std::log((1 - settings.probabilityCutoff)/settings.probabilityCutoff)) {
    }

    inline void merge(NodeData& result, const std::vector<NodeData>& options, const int departureTime) const noexcept {
        result = NodeData();
        if (options.empty()) return;
        if (options.size() == 1) {
            result = options[0];
            return;
        }
        PerceivedTime minExpectedPAT = INFTY;
        for (const NodeData& option : options) {
            if (option.arrivalTime == departureTime) {
                result = option;
                return;
            }
            result.arrivalTime = std::min(option.arrivalTime, result.arrivalTime);
            minExpectedPAT = std::min(option.expectedPAT, minExpectedPAT);
        }
        const double scalingFactor = getScalingFactor(result.arrivalTime - departureTime);
        PerceivedTime skew = 0;
        for (const NodeData& option : options) {
            AssertMsg(option.expectedPAT < Unreachable, "Invalid option!");
            if (option.expectedPAT - minExpectedPAT > scalingFactor * probabilityFactor) continue;
            skew += std::exp((minExpectedPAT - option.expectedPAT)/scalingFactor);
        }
        skew = scalingFactor * std::log(skew);
        result.expectedPAT = minExpectedPAT - skew;
    }

    inline void merge(NodeData& result, const NodeData& a, const NodeData& b, const int departureTime) const noexcept {
        if (a.arrivalTime >= INFTY || b.arrivalTime == departureTime) {
            result = b;
            return;
        }
        if (b.arrivalTime >= INFTY || a.arrivalTime == departureTime) {
            result = a;
            return;
        }
        result.arrivalTime = std::min(a.arrivalTime, b.arrivalTime);
        const double scalingFactor = getScalingFactor(result.arrivalTime - departureTime);
        const PerceivedTime minPAT = std::min(a.expectedPAT, b.expectedPAT);
        const PerceivedTime maxPAT = std::max(a.expectedPAT, b.expectedPAT);
        if (maxPAT - minPAT > scalingFactor * probabilityFactor) {
            result.expectedPAT = minPAT;
            return;
        }
        const PerceivedTime maxLogitValue = std::exp((minPAT - maxPAT)/scalingFactor);
        const PerceivedTime skew = scalingFactor * std::log(maxLogitValue + 1);
        result.expectedPAT = minPAT - skew;
    }

    inline std::vector<int> distribution(const std::vector<PerceivedTime> pats, const int travelTime) const noexcept {
        const PerceivedTime minPAT = Vector::min(pats);
        AssertMsg(minPAT < Unreachable, "All options are invalid!");
        std::vector<int> result(pats.size() + 1, 0);
        if (travelTime == 0) {
            for (size_t i = 0; i < pats.size(); i++) {
                if (pats[i] == minPAT) {
                    result[i] = 1;
                    result.back()++;
                }
            }
        } else {
            const double scalingFactor = getScalingFactor(travelTime);
            for (size_t i = 0; i < pats.size(); i++) {
                if (pats[i] >= Unreachable) continue;
                if (pats[i] - minPAT > scalingFactor * probabilityFactor) continue;
                result[i] = logitValue(pats[i], minPAT, scalingFactor);
                AssertMsg(result[i] >= 0, "Logit value is negative ( " << String::prettyInt(result[i]) << ")!");
                result.back() += result[i];
            }
        }
        AssertMsg(result.back() > 0, "Probability of all options cannot be zero!");
        return result;
    }

    inline std::array<int, 3> distribution(const PerceivedTime a, const PerceivedTime b, const int travelTime) const noexcept {
        if (travelTime == 0) {
            if (a == b) return std::array<int, 3>{1, 1, 2};
            else if (a < b) return std::array<int, 3>{1, 0, 1};
            else return std::array<int, 3>{0, 1, 1};
        }
        AssertMsg(a < Unreachable || b < Unreachable, "Both options are invalid!");
        const PerceivedTime minPAT = std::min(a, b);
        const double scalingFactor = getScalingFactor(travelTime);
        if (b > a + scalingFactor * probabilityFactor) return std::array<int, 3>{ 1, 0, 1 };
        if (a > b + scalingFactor * probabilityFactor) return std::array<int, 3>{ 0, 1, 1 };
        const int valueA = logitValue(a, minPAT, scalingFactor);
        const int valueB = logitValue(b, minPAT, scalingFactor);
        AssertMsg(valueA + valueB > 0, "Probability of all options cannot be zero (" << a << ", " << b << ")!");
        return std::array<int, 3>{valueA, valueB, valueA + valueB};
    }

private:
    inline int logitValue(const PerceivedTime pat, const PerceivedTime minPAT, const double scalingFactor) const noexcept {
        if (pat >= Unreachable) return 0;
        return std::exp(10 + ((minPAT - pat)/scalingFactor));
    }

    inline double getScalingFactor(const int travelTime) const noexcept {
        return std::sqrt(6 * beta * referenceTravelTime * travelTime)/M_PI;
    }

    const double beta;
    const int referenceTravelTime;
    const double probabilityFactor;
};

}
