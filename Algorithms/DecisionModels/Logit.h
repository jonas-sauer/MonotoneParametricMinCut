#pragma once

#include <cmath>
#include <array>
#include <vector>

#include "../../DataStructures/Assignment/Settings.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Vector.h"

namespace DecisionModels {

class Logit {

public:
    Logit(const int delayTolerance, const double beta = 1) :
        delayTolerance(delayTolerance),
        beta(beta) {
    }
    Logit(const Assignment::Settings& settings) :
        Logit(settings.delayTolerance, settings.beta) {
    }

    inline void cumulativeDistribution(const std::vector<int>& values, std::vector<int>& result) const noexcept {
        if (values.empty()) return;
        const std::array<int, 2> minValues = Vector::twoSmallestValues(values);
        result.resize(values.size());
        if (minValues[1] - minValues[0] > delayTolerance) {
            result[0] = (values[0] == minValues[0]) ? 1 : 0;
            for (size_t i = 1; i < values.size(); i++) {
                result[i] = (values[i] == minValues[0]) ? 1 : result[i - 1];
            }
        } else {
            result[0] = (values[0] - minValues[0] <= delayTolerance) ? logitValue(values[0], minValues[0]) : 0;
            for (size_t i = 1; i < values.size(); i++) {
                result[i] = result[i - 1] + ((values[i] - minValues[0] <= delayTolerance) ? logitValue(values[i], minValues[0]) : 0);
                AssertMsg(result[i] >= result[i - 1], "Accumulated values are decreasing from " << String::prettyInt(result[i - 1]) << " to " << String::prettyInt(result[i]) << "!");
            }
        }
        AssertMsg(result.back() > 0, "Probability of all options cannot be zero!");
    }

    inline std::vector<int> cumulativeDistribution(const std::vector<int>& values) const noexcept {
        std::vector<int> result;
        cumulativeDistribution(values, result);
        return result;
    }

    inline std::array<int, 2> cumulativeDistribution(const int a, const int b) const noexcept {
        if (b - a > delayTolerance) {
            return std::array<int, 2>{1, 1};
        } else if (a - b > delayTolerance) {
            return std::array<int, 2>{0, 1};
        } else {
            const int minValue = std::min(a, b);
            const int valueA = logitValue(a, minValue);
            const int valueB = logitValue(b, minValue);
            AssertMsg(valueA + valueB > 0, "Probability of all options cannot be zero (" << a << ", " << b << ")!");
            return std::array<int, 2>{valueA, valueA + valueB};
        }
    }

    inline void distribution(const std::vector<int>& values, std::vector<int>& result) const noexcept {
        if (values.empty()) return;
        const std::array<int, 2> minValues = Vector::twoSmallestValues(values);
        result.resize(values.size() + 1);
        result.back() = 0;
        if (minValues[1] - minValues[0] > delayTolerance) {
            for (size_t i = 0; i < values.size(); i++) {
                result[i] = (values[i] == minValues[0]) ? 1 : 0;
                result.back() += result[i];
            }
        } else {
            for (size_t i = 0; i < values.size(); i++) {
                result[i] = (values[i] - minValues[0] <= delayTolerance) ? logitValue(values[i], minValues[0]) : 0;
                AssertMsg(result[i] >= 0, "Logit value is negative ( " << String::prettyInt(result[i]) << ")!");
                result.back() += result[i];
            }
        }
        AssertMsg(result.back() > 0, "Probability of all options cannot be zero!");
    }

    inline std::vector<int> distribution(const std::vector<int>& values) const noexcept {
        std::vector<int> result;
        distribution(values, result);
        return result;
    }

    inline std::array<int, 3> distribution(const int a, const int b) const noexcept {
        if (b - a > delayTolerance) {
            return std::array<int, 3>{1, 0, 1};
        } else if (a - b > delayTolerance) {
            return std::array<int, 3>{0, 1, 1};
        } else {
            const int minValue = std::min(a, b);
            const int valueA = logitValue(a, minValue);
            const int valueB = logitValue(b, minValue);
            AssertMsg(valueA + valueB > 0, "Probability of all options cannot be zero (" << a << ", " << b << ")!");
            return std::array<int, 3>{valueA, valueB, valueA + valueB};
        }
    }

private:
    inline int logitValue(const int value, const int minValue) const noexcept {
        return std::exp(10 + (beta * (minValue - value)));
    }

private:
    const int delayTolerance;
    const double beta;

};

}
