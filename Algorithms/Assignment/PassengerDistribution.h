#pragma once

#include <array>
#include <random>
#include <vector>

#include "../../Helpers/Types.h"

namespace Assignment {

using Random = std::minstd_rand;

struct SampleElement {
    SampleElement() : index(-1), value(1) {}

    SampleElement(const size_t index, const size_t weight, Random& random) :
        index(index),
        value(weight == 0 ? 1 : -std::log(random()) / static_cast<double>(weight)) {
    }

    size_t index;
    double value;

    inline bool operator<(const SampleElement& other) const noexcept {
        return value < other.value;
    }
};

inline std::vector<size_t> getGroupSizesCumulative(const std::vector<int>& values, const size_t numberOfPassengers, Random& random, std::vector<size_t>& groupSizes, std::vector<SampleElement>& sampleElements) noexcept {
    const int64_t valueSum = values.back();
    size_t rejectedPassengers = 0;
    for (size_t i = values.size() - 1; i > 0; i--) {
        const int64_t value = values[i] - values[i - 1];
        groupSizes[i] = (value * numberOfPassengers) / valueSum;
        rejectedPassengers += groupSizes[i];
        sampleElements[i] = SampleElement(i, (values[i] * numberOfPassengers) % valueSum, random);
    }
    if (!groupSizes.empty()) {
        const int64_t value = values[0];
        groupSizes[0] = (value * numberOfPassengers) / valueSum;
        rejectedPassengers += groupSizes[0];
        sampleElements[0] = SampleElement(0, (value * numberOfPassengers) % valueSum, random);
    }
    const size_t remainingPassengers = numberOfPassengers - rejectedPassengers;
    if (remainingPassengers > 0) {
        std::nth_element(sampleElements.begin(), sampleElements.begin() + remainingPassengers - 1, sampleElements.end());
        for (size_t i = 0; i < remainingPassengers; i++) {
            groupSizes[sampleElements[i].index]++;
        }
    }
    AssertMsg(Vector::sum(groupSizes) == numberOfPassengers, "New groups should comprise " << numberOfPassengers << " passengers, but is " << Vector::sum(groupSizes) << "!");
    return groupSizes;
}

inline std::vector<size_t> getGroupSizesCumulative(const std::vector<int>& values, const size_t numberOfPassengers, Random& random) noexcept {
    std::vector<size_t> groupSizes(values.size(), 0);
    std::vector<SampleElement> sampleElements(values.size() - 1);
    getGroupSizesCumulative(values, numberOfPassengers, random, groupSizes, sampleElements);
    return groupSizes;
}

inline std::vector<size_t> getGroupSizes(const std::vector<int>& values, const size_t numberOfPassengers, Random& random, std::vector<size_t>& groupSizes, std::vector<SampleElement>& sampleElements) noexcept {
    const int64_t valueSum = values.back();
    size_t rejectedPassengers = 0;
    for (size_t i = values.size() - 2; i != size_t(-1); i--) {
        groupSizes[i] = (values[i] * numberOfPassengers) / valueSum;
        rejectedPassengers += groupSizes[i];
        sampleElements[i] = SampleElement(i, (values[i] * numberOfPassengers) % valueSum, random);
    }
    const size_t remainingPassengers = numberOfPassengers - rejectedPassengers;
    if (remainingPassengers > 0) {
        std::nth_element(sampleElements.begin(), sampleElements.begin() + remainingPassengers - 1, sampleElements.end());
        for (size_t i = 0; i < remainingPassengers; i++) {
            groupSizes[sampleElements[i].index]++;
        }
    }
    AssertMsg(Vector::sum(groupSizes) == numberOfPassengers, "New groups should comprise " << numberOfPassengers << " passengers, but is " << Vector::sum(groupSizes) << "!");
    return groupSizes;
}

inline std::vector<size_t> getGroupSizes(const std::vector<int>& values, const size_t numberOfPassengers, Random& random) noexcept {
    std::vector<size_t> groupSizes(values.size() - 1, 0);
    std::vector<SampleElement> sampleElements(values.size() - 1);
    getGroupSizes(values, numberOfPassengers, random, groupSizes, sampleElements);
    return groupSizes;
}

inline std::array<int, 2> getGroupSizesCumulative(const std::array<int, 2> values, const int64_t numberOfPassengers, Random& random) noexcept {
    const int64_t valueSum = values[1];
    std::array<int, 2> groupSizes = {static_cast<int>((values[0] * numberOfPassengers) / valueSum), static_cast<int>(((values[1] - values[0]) * numberOfPassengers) / valueSum)};
    std::array<int, 2> remainderValues = {static_cast<int>((values[0] * numberOfPassengers) % valueSum), static_cast<int>(((values[1] - values[0]) * numberOfPassengers) % valueSum)};
    AssertMsg(remainderValues[0] + remainderValues[1] == 0 || remainderValues[0] + remainderValues[1] == valueSum, "Remainder sum is " << (remainderValues[0] + remainderValues[1]) << ", but should be " << valueSum << "!");
    if (remainderValues[0] != 0) {
        if (static_cast<int>(random() % (remainderValues[0] + remainderValues[1])) < remainderValues[0]) {
            groupSizes[0]++;
        } else {
            groupSizes[1]++;
        }
    }
    AssertMsg(groupSizes[0] + groupSizes[1] == numberOfPassengers, "New groups should comprise " << numberOfPassengers << " passengers, but is " << (groupSizes[0] + groupSizes[1]) << "!");
    return groupSizes;
}

inline std::array<int, 2> getGroupSizes(const std::array<int, 3> values, const int64_t numberOfPassengers, Random& random) noexcept {
    const int64_t valueSum = values[2];
    std::array<int, 2> groupSizes = {static_cast<int>((values[0] * numberOfPassengers) / valueSum), static_cast<int>((values[1] * numberOfPassengers) / valueSum)};
    std::array<int, 2> remainderValues = {static_cast<int>((values[0] * numberOfPassengers) % valueSum), static_cast<int>((values[1] * numberOfPassengers) % valueSum)};
    AssertMsg(remainderValues[0] + remainderValues[1] == 0 || remainderValues[0] + remainderValues[1] == valueSum, "Remainder sum is " << (remainderValues[0] + remainderValues[1]) << ", but should be " << valueSum << "!");
    if (remainderValues[0] != 0) {
        if (static_cast<int>(random() % (remainderValues[0] + remainderValues[1])) < remainderValues[0]) {
            groupSizes[0]++;
        } else {
            groupSizes[1]++;
        }
    }
    AssertMsg(groupSizes[0] + groupSizes[1] == numberOfPassengers, "New groups should comprise " << numberOfPassengers << " passengers, but is " << (groupSizes[0] + groupSizes[1]) << "!");
    return groupSizes;
}

}
