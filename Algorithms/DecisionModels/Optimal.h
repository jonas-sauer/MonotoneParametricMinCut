#pragma once

#include <array>
#include <vector>

#include "../../DataStructures/Assignment/Settings.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Vector.h"

namespace DecisionModels {

class Optimal {

public:
    Optimal(const int delayTolerance) :
        deterministic(delayTolerance > 0) {
    }
    Optimal(const Assignment::Settings& settings) :
        Optimal(settings.delayTolerance) {
    }

    inline void cumulativeDistribution(const std::vector<int>& values, std::vector<int>& result) const noexcept {
        if (values.empty()) return;
        const int minValue = Vector::min(values);
        result.resize(values.size());
        if (deterministic) {
            int currentValue = 0;
            for (size_t i = 0; i < values.size(); i++) {
                if (values[i] == minValue) {
                    currentValue = 1;
                }
                result[i] = currentValue;
            }
        } else {
            int currentValue = 0;
            for (size_t i = 0; i < values.size(); i++) {
                if (values[i] == minValue) {
                    currentValue++;
                }
                result[i] = currentValue;
            }
        }
    }

    inline std::vector<int> cumulativeDistribution(const std::vector<int>& values) const noexcept {
        std::vector<int> result;
        cumulativeDistribution(values, result);
        return result;
    }

    inline std::array<int, 2> cumulativeDistribution(const double a, const double b) const noexcept {
        if ((deterministic) || (a != b)) {
            if (a <= b) {
                return std::array<int, 2>{1, 1};
            } else {
                return std::array<int, 2>{0, 1};
            }
        } else {
            return std::array<int, 2>{1, 2};
        }
    }

    inline void distribution(const std::vector<int>& values, std::vector<int>& result) const noexcept {
        if (values.empty()) return;
        const int minValue = Vector::min(values);
        result.resize(values.size() + 1);
        result.back() = 0;
        if (deterministic) {
            bool minimumFound = false;
            for (size_t i = 0; i < values.size(); i++) {
                if (!minimumFound && values[i] == minValue) {
                    minimumFound = true;
                    result[i] = 1;
                } else {
                    result[i] = 0;
                }
            }
            result.back() = 1;
        } else {
            for (size_t i = 0; i < values.size(); i++) {
                result[i] = (values[i] == minValue) ? 1 : 0;
                result.back() += result[i];
            }
        }
    }

    inline std::vector<int> distribution(const std::vector<int>& values) const noexcept {
        std::vector<int> result;
        distribution(values, result);
        return result;
    }

    inline std::array<int, 3> distribution(const double a, const double b) const noexcept {
        if ((deterministic) || (a != b)) {
            if (a <= b) {
                return std::array<int, 3>{1, 0, 1};
            } else {
                return std::array<int, 3>{0, 1, 1};
            }
        } else {
            return std::array<int, 3>{1, 1, 2};
        }
    }

private:
    const bool deterministic;

};

}
