#pragma once

#include <concepts>
#include <numeric>
#include <cstdlib>

#include "../../Helpers/Types.h"
#include "../../Helpers/IO/Serialization.h"

namespace pmf {
    inline constexpr double epsilon = 1E-9;

    /**
     * Function to check if two double values are almost equal
     * @param a
     * @param b
     * @return true if |a - b| < epsilon
     */
    inline bool doubleEqualAbs(const double a, const double b) noexcept {
        return std::abs(a - b) < epsilon;
    }

    /**
     * Function to check if two double values are almost equal taking into account growing imprecision with higher numbers
     * @param a
     * @param b
     * @return true if |a - b| < |a| * epsilon
     */
    inline bool doubleEqualRelative(const double a, const double b) noexcept {
        return std::abs(a - b) < std::abs(a * epsilon);
    }

    /**
     * Function to check a < b while accounting for double inequality
     * @param a
     * @param b
     * @return true iff a + epsilon < b
     */
    inline bool doubleLessThanAbs(const double a, const double b) noexcept {
        return a + epsilon < b;
    }

    /**
     * Function to check a > 0 while accounting for double imprecision
     * @param a
     * @param b
     * @return true iff a > epsilon
     */
    inline bool doubleIsPositive(const double a) noexcept {
        return a > epsilon;
    }

    template<typename T>
    inline bool isNumberPositive(const T a) noexcept {
        if constexpr (std::is_floating_point<T>::value) {
            return a > epsilon;
        } else {
            return a > 0;
        }
    }

    template<typename T>
    inline bool isNumberNegative(const T a) noexcept {
        if constexpr (std::is_floating_point<T>::value) {
            return a < -epsilon;
        } else {
            return a < 0;
        }
    }

    template<typename T>
    inline bool areNumbersEqual(const T a, const T b) noexcept {
        if constexpr (std::is_floating_point<T>::value) {
            return std::abs(a - b) < epsilon;
        } else {
            return a == b;
        }
    }
} // namespace pmf