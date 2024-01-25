#pragma once

#include <cmath>
#include <random>

#include "Beta.h"

namespace Assignment::MRUM {

class PERT {

public:
    PERT(const double min, const double max, const double mean, const double lambda, const uint32_t seed) :
        min(min),
        span(max - min),
        beta(1 + (lambda * ((mean - min) / span)), 1 + (lambda * ((max - mean) / span)), seed) {
    }

    inline double operator()() const noexcept {
        return min + (span * beta());
    }

    const double min;
    const double span;

    mutable Beta beta;

};

}
