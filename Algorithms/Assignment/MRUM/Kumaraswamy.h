#pragma once

#include <cmath>
#include <random>

namespace Assignment::MRUM {

class Kumaraswamy {

public:
    Kumaraswamy(const double min, const double max, const double a, const double b, const uint32_t seed) :
        min(min),
        span(max - min),
        exponentA(1 / a),
        exponentB(1 / b),
        randomGenerator(seed),
        uniformDistribution(0.0, 1.0) {
    }

    inline double operator()() const noexcept {
        return min + (span  * std::pow(1 - std::pow(1 - uniformDistribution(randomGenerator), exponentB), exponentA));
    }

    const double min;
    const double span;
    const double exponentA;
    const double exponentB;

    mutable std::mt19937 randomGenerator;
    mutable std::uniform_real_distribution<double> uniformDistribution;

};

}
