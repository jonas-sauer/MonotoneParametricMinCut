#pragma once

#include <array>
#include <cmath>

#include "../../../DataStructures/Assignment/AssignmentData.h"
#include "../../../DataStructures/Assignment/ConnectionLoadData.h"

namespace Assignment::Capacities {

inline double applyIteration(const double oldValue, const double newValue, const std::array<double, 3>& stepSizes) noexcept {
    return (stepSizes[0] * oldValue + stepSizes[1] * newValue)/stepSizes[2];
}

class StepIterator {
public:
    StepIterator() : iteration(1.0) {}

    virtual ~StepIterator() {}

    inline virtual std::array<double, 3> getStepSizes() noexcept {
        return std::array<double, 3>{iteration - 1, 1.0, iteration};
    }

    inline void operator++() noexcept {
        iteration++;
    }

    inline virtual void reset() noexcept {
        iteration = 1;
    }

    inline size_t getIteration() const noexcept {
        return iteration;
    }

protected:
    double iteration;
};

struct WeightedStepIterator : public StepIterator {
    WeightedStepIterator(const double stepSize) : StepIterator(), stepSize(stepSize), sum(0) {}

    inline std::array<double, 3> getStepSizes() noexcept override {
        const double newTerm = std::pow(iteration, stepSize);
        sum += newTerm;
        return std::array<double, 3>{sum - newTerm, newTerm, sum};
    }

    inline void reset() noexcept {
        StepIterator::reset();
        sum = 0;
    }

    const double stepSize;
    double sum;
};

struct SelfRegulatingStepIterator : public StepIterator {
public:
    SelfRegulatingStepIterator(const AssignmentData& currentAssignmentData, const std::vector<ConnectionLoadData>& loadData, const size_t numberOfConnections, const double passengerMultiplier, const double speedupStepSize, const double slowdownStepSize) :
        StepIterator(),
        currentAssignmentData(currentAssignmentData),
        loadData(loadData),
        numberOfConnections(numberOfConnections),
        passengerMultiplier(passengerMultiplier),
        speedupStepSize(speedupStepSize),
        slowdownStepSize(slowdownStepSize),
        iteration(1),
        sum(0),
        diff(0) {
    }

    inline std::array<double, 3> getStepSizes() noexcept override {
        const double newDiff = calculateDiff();
        sum += (newDiff >= diff) ? speedupStepSize : slowdownStepSize;
        diff = newDiff;
        return std::array<double, 3>{sum - 1, 1.0, sum};
    }

    inline void reset() noexcept {
        StepIterator::reset();
        sum = 0;
        diff = 0;
    }

private:
    inline double calculateDiff() const noexcept {
        double maxDiff = 0.0;
        for (ConnectionId connection(0); connection <  numberOfConnections; connection++) {
            const double currentLoad = currentAssignmentData.getConnectionLoad(connection) / passengerMultiplier;
            const double diff = fabs(currentLoad - loadData[connection].load)/loadData[connection].capacity;
            maxDiff = std::max(maxDiff, diff);
        }
        return maxDiff;
    }

    const AssignmentData& currentAssignmentData;
    const std::vector<ConnectionLoadData>& loadData;
    const size_t numberOfConnections;
    const double passengerMultiplier;

    const double speedupStepSize;
    const double slowdownStepSize;

    double iteration;
    double sum;
    double diff;
};

}
