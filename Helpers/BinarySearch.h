#pragma once

#include <iostream>
#include <algorithm>

#include "Assert.h"
#include "Vector/Vector.h"

template<typename FUNCTION>
inline int firstFalseIndex(int minI, int maxI, const FUNCTION& f, const int maxIntervalSize = 16) {
    while (maxI - minI >= maxIntervalSize) {
        Assert(minI >= 0, "minI is negative!");
        Assert(minI < maxI, "minI exceeds maxI!");
        const int j = minI + ((maxI - minI) >> 1);
        if (f(j)) {
            minI = j + 1;
        } else {
            maxI = j - 1;
        }
    }
    for (int i = minI; i <= maxI; i++) {
        if (!f(i)) return i;
    }
    return maxI + 1;
}

template<typename VECTOR, typename FUNCTION>
inline int firstFalseIndex(const VECTOR& vector, const FUNCTION& f, const int maxIntervalSize = 16) {
    return firstFalseIndex(0, vector.size() - 1, f, maxIntervalSize);
}

template<typename FUNCTION>
inline int firstFalseIndex(const int size, const FUNCTION& f, const int maxIntervalSize = 16) {
    return firstFalseIndex(0, size - 1, f, maxIntervalSize);
}

template<typename FUNCTION>
inline int lastTrueIndex(const int minI, const int maxI, const FUNCTION& f, const int maxIntervalSize = 16) {
    return firstFalseIndex(minI, maxI, f, maxIntervalSize) - 1;
}

template<typename VECTOR, typename FUNCTION>
inline int lastTrueIndex(const VECTOR& vector, const FUNCTION& f, const int maxIntervalSize = 16) {
    return firstFalseIndex(0, vector.size() - 1, f, maxIntervalSize) - 1;
}

template<typename FUNCTION>
inline int lastTrueIndex(const int size, const FUNCTION& f, const int maxIntervalSize = 16) {
    return firstFalseIndex(0, size - 1, f, maxIntervalSize) - 1;
}

template<typename VALUE, typename VECTOR, typename LESS>
inline int firstIndexGreaterOrEqual(const int size, const VECTOR& vector, const VALUE& value, const LESS& less, const int maxIntervalSize = 30) {
    int minI = 0;
    int maxI = size - 1;
    while (maxI - minI >= maxIntervalSize) {
        Assert(minI >= 0, "minI is negative!");
        Assert(maxI < size, "maxI exceeds size!");
        Assert(minI < maxI, "minI exceeds maxI!");
        const VALUE& minValue = vector(minI);
        if (!less(minValue, value)) return minI;
        const VALUE& maxValue = vector(maxI);
        if (less(maxValue, value)) return maxI + 1;
        Assert(!less(maxValue, minValue), "minValue exceeds maxValue!");
        const int i = minI + ((maxI - minI) * ((value - minValue) / static_cast<double>(maxValue - minValue)));
        Assert(i >= minI, "i must be in [minI, maxI]!");
        Assert(i <= maxI, "i must be in [minI, maxI]!");
        if (less(vector(i), value)) {
            minI = i + 1;
        } else {
            maxI = i - 1;
        }
        Assert(minI >= 0, "minI is negative!");
        Assert(maxI < size, "maxI exceeds size!");
        Assert(minI < maxI, "minI exceeds maxI!");
        const int j = minI + ((maxI - minI) >> 1);
        if (less(vector(j), value)) {
            minI = j + 1;
        } else {
            maxI = j - 1;
        }
    }
    for (int i = minI; i <= maxI; i++) {
        if (!less(vector(i), value)) return i;
    }
    return maxI + 1;
}

template<typename VALUE, typename VECTOR, typename LESS>
inline int firstIndexGreaterOrEqual(const VECTOR& vector, const VALUE& value, const LESS& less, const int maxIntervalSize = 30) {
    return firstIndexGreaterOrEqual(vector.size(), [&](const int i){return vector[i];}, value, less, maxIntervalSize);
}

template<typename VALUE, typename VECTOR>
inline int firstIndexGreaterOrEqual(const VECTOR& vector, const VALUE& value, const int maxIntervalSize = 100) {
    return firstIndexGreaterOrEqual(vector.size(), [&](const int i){return vector[i];}, value, [](const VALUE& a, const VALUE& b){return a < b;}, maxIntervalSize);
}

template<typename VALUE, typename VECTOR, typename LESS>
inline int lastIndexLess(const int size, const VECTOR& vector, const VALUE& value, const LESS& less, const int maxIntervalSize = 30) {
    return firstIndexGreaterOrEqual(size, vector, value, less, maxIntervalSize) - 1;
}

template<typename VALUE, typename VECTOR, typename LESS>
inline int lastIndexLess(const VECTOR& vector, const VALUE& value, const LESS& less, const int maxIntervalSize = 30) {
    return firstIndexGreaterOrEqual(vector.size(), [&](const int i){return vector[i];}, value, less, maxIntervalSize) - 1;
}

template<typename VALUE, typename VECTOR>
inline int lastIndexLess(const VECTOR& vector, const VALUE& value, const int maxIntervalSize = 100) {
    return firstIndexGreaterOrEqual(vector.size(), [&](const int i){return vector[i];}, value, [](const VALUE& a, const VALUE& b){return a < b;}, maxIntervalSize) - 1;
}
