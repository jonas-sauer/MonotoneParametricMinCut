#pragma once

#include <concepts>
#include <numeric>
#include <cstdlib>

#include "Types.h"
#include "IO/Serialization.h"

inline constexpr double FLOAT_EPS = 1E-9;

/**
 * Check if a > 0 while accounting for floating-point imprecision
 */
template<typename T>
inline bool isNumberPositive(const T a) noexcept {
    if constexpr (std::is_floating_point<T>::value) {
        return a > FLOAT_EPS;
    } else {
        return a > 0;
    }
}

/**
 * Check if a < 0 while accounting for floating-point imprecision
 */
template<typename T>
inline bool isNumberNegative(const T a) noexcept {
    if constexpr (std::is_floating_point<T>::value) {
        return a < -FLOAT_EPS;
    } else {
        return a < 0;
    }
}

/**
 * Check if a == b while accounting for floating-point imprecision
 */
template<typename T>
inline bool areNumbersEqualAbsolute(const T a, const T b) noexcept {
    if constexpr (std::is_floating_point<T>::value) {
        return std::abs(a - b) < FLOAT_EPS;
    } else {
        return a == b;
    }
}

/**
 * Check if a == b while accounting for growing floating-point imprecision with larger numbers
 */
template<typename T>
inline bool areNumbersEqualRelative(const T a, const T b) noexcept {
    if constexpr (std::is_floating_point<T>::value) {
        return std::abs(a - b) < std::abs(a * FLOAT_EPS);
    } else {
        return a == b;
    }
}

/**
 * Check if a < b while accounting for floating-point imprecision
 */
template<typename T>
inline bool isNumberLessThanAbsolute(const T a, const T b) noexcept {
    if constexpr (std::is_floating_point<T>::value) {
        return a + FLOAT_EPS < b;
    } else {
        return a < b;
    }
}