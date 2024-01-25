#pragma once

#include <vector>
#include <algorithm>

#include "../Helpers/Assert.h"
#include "../Helpers/Helpers.h"
#include "../Helpers/IO/Serialization.h"

class ArithmeticProgression {

public:
    struct Iterator {
        Iterator(const int value, const size_t increment) : value(value), increment(increment) {}
        inline bool operator!=(const Iterator& other) const {return value != other.value;}
        inline int operator*() const {return value;}
        inline Iterator& operator++() {value += increment; return *this;}
        inline Iterator& operator+(const int n) {value += (n * increment); return *this;}
        int value;
        const size_t increment;
    };

public:
    ArithmeticProgression(const int minimum = 0, const size_t numberOfElements = 0, const size_t increment = 1) :
        minimum(minimum),
        numberOfElements(numberOfElements),
        increment(increment) {
        Assert(increment > 0, "Step has to be greater than 0!");
    }
    ArithmeticProgression(IO::Deserialization& deserialize) {
        this->deserialize(deserialize);
    }

    inline operator std::vector<int>() const noexcept {
        std::vector<int> result;
        for (int i = minimum; i != firstNotIncluded(); i += increment) {
            result.emplace_back(i);
        }
        return result;
    }

    inline int firstIncluded() const noexcept {
        return minimum;
    }

    inline int firstNotIncluded() const noexcept {
        return minimum + (numberOfElements * increment);
    }

    inline int lastIncluded() const noexcept {
        return firstNotIncluded() - increment;
    }

    inline int lastNotIncluded() const noexcept {
        return firstIncluded() - increment;
    }

    inline int front() const noexcept {
        return firstIncluded();
    }

    inline int back() const noexcept {
        return lastIncluded();
    }

    inline bool empty() const noexcept {
        return numberOfElements == 0;
    }

    inline size_t size() const noexcept {
        return numberOfElements;
    }

    inline bool contains(const int number) const noexcept {
        if ((number < minimum) || (((number - minimum) % increment) != 0)) return false;
        return number < firstNotIncluded();
    }

    inline Iterator begin() const noexcept {
        return Iterator(minimum, increment);
    }

    inline Iterator end() const noexcept {
        return Iterator(firstNotIncluded(), increment);
    }

    inline int operator[](const size_t i) const noexcept {
        return minimum + (i * increment);
    }

    inline void shift(const int amount) noexcept {
        minimum += amount;
    }

    inline ArithmeticProgression shifted(const int amount) const noexcept {
        return ArithmeticProgression(minimum + amount, numberOfElements, increment);
    }

    inline void prepend(const size_t n = 1) noexcept {
        Assert(numberOfElements > 0, "N has to be greater than 0!");
        minimum -= (n * increment);
        numberOfElements += n;
    }

    inline void append(const size_t n = 1) noexcept {
        Assert(numberOfElements > 0, "N has to be greater than 0!");
        numberOfElements += n;
    }

    friend std::ostream& operator<<(std::ostream& out, const ArithmeticProgression& ap) {
        return out << "ArithmeticProgression(" << ap.minimum  << ", " << ap.numberOfElements  << ", " << ap.increment << ")";
    }

    inline void serialize(IO::Serialization& serialize) const noexcept {
        serialize(minimum, numberOfElements, increment);
    }

    inline void deserialize(IO::Deserialization& deserialize) noexcept {
        deserialize(minimum, numberOfElements, increment);
    }

public:
    int minimum;
    size_t numberOfElements;
    size_t increment;
};

class CoverByArithmeticProgressions {

public:
    CoverByArithmeticProgressions(const std::vector<int>& numbers) : input(numbers) {
        sort(input);
        for (size_t i = 1; i < input.size(); i++) {
            Assert(input[i - 1] < input[i], "The vector of numbers contains duplicates (" << input[i - 1] << ", " << input[i] << ")!");
        }
    }

    inline std::vector<ArithmeticProgression> greedyCover() const noexcept {
        std::vector<ArithmeticProgression> result;
        std::vector<bool> isCovered(input.size(), false);
        for (size_t i = 0; i < input.size(); i++) {
            if (isCovered[i]) continue;
            ArithmeticProgression bestCover(input[i], 1, 1);
            size_t bestCoverCount = 1;
            for (size_t j = i + 1; j < input.size(); j++) {
                ArithmeticProgression cover(input[i], 2, input[j] - input[i]);
                size_t coverCount = isCovered[j] ? 1 : 2;
                size_t coverSize = coverCount;
                size_t k = j + 1;
                while (k < input.size() && cover.firstNotIncluded() >= input[k]) {
                    if (cover.firstNotIncluded() == input[k]) {
                        cover.append();
                        if (!isCovered[k]) {
                            coverCount++;
                            coverSize = cover.numberOfElements;
                        }
                    }
                    k++;
                }
                cover.numberOfElements = coverSize;
                if ((coverCount > bestCoverCount) || ((coverCount == bestCoverCount) && (cover.size() < bestCover.size()))) {
                    bestCover = cover;
                    bestCoverCount = coverCount;
                }
            }
            for (size_t j = i; j < input.size(); j++) {
                if (input[j] > bestCover.lastIncluded()) break;
                if (bestCover.contains(input[j])) isCovered[j] = true;
            }
            result.emplace_back(bestCover);
        }
        return result;
    }

private:
    std::vector<int> input;

};