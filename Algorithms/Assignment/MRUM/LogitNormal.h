#pragma once

#include <cmath>
#include <random>

namespace Assignment::MRUM {

class LogitNormal {

public:
    LogitNormal(const double min, const double max, const double mean, const double standardDeviation, const uint32_t seed) :
        min(min),
        span(max - min),
        randomGenerator(seed),
        normalDistribution(mean - ((max + min) / 2.0), standardDeviation) {
    }

    inline double operator()() const noexcept {
        return min + (span / (1 + std::exp(-normalDistribution(randomGenerator))));
    }

    const double min;
    const double span;

    mutable std::mt19937 randomGenerator;
    mutable std::normal_distribution<double> normalDistribution;

};

}
