#pragma once

#include <random>

namespace Assignment::MRUM {

class Beta {

public:
    Beta(const double alpha, const double beta, const uint32_t seed) :
        randomGenerator(seed),
        gammaDistributionX(alpha, 1.0),
        gammaDistributionY(beta, 1.0) {
    }

    inline double operator()() const noexcept {
        const double x = gammaDistributionX(randomGenerator);
        return x / (x + gammaDistributionY(randomGenerator));
    }

    mutable std::mt19937 randomGenerator;
    mutable std::gamma_distribution<double> gammaDistributionX;
    mutable std::gamma_distribution<double> gammaDistributionY;

};

}
