#pragma once

#include <array>
#include <vector>

#include "../../DataStructures/Assignment/Settings.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Vector.h"

namespace DecisionModels {

class Linear {

public:
    Linear(const int delayTolerance) :
        delayTolerance(delayTolerance),
        delayValue(delayTolerance) {
    }
    Linear(const int delayTolerance, const int delayValue) :
        delayTolerance(delayTolerance),
        delayValue(delayValue) {
        AssertMsg(delayValue >= delayTolerance, "The delayValue has to be greater or equal to delayTolerance!");
    }
    Linear(const Assignment::Settings& settings) :
        Linear(settings.delayTolerance, settings.delayValue) {
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
            result[0] = gain(values[0], minValues);
            for (size_t i = 1; i < values.size(); i++) {
                result[i] = result[i - 1] + gain(values[i], minValues);
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
            return std::array<int, 2>{b - a + delayValue, delayValue * 2};
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
                result[i] = gain(values[i], minValues);
                AssertMsg(result[i] >= 0, "Gain is negative ( " << String::prettyInt(result[i]) << ")!");
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
            return std::array<int, 3>{b - a + delayValue, a - b + delayValue, 2 * delayValue};
        }
    }

private:
    inline int gain(const int value, const std::array<int, 2>& minValues) const noexcept {
        if (value - minValues[0] > delayTolerance) return 0;
        const int minValue = (value == minValues[0]) ? minValues[1] : minValues[0];
        return minValue - value + delayValue;
    }

private:
    const int delayTolerance;
    const int delayValue;

};

}
