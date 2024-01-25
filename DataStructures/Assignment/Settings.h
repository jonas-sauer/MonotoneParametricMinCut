#pragma once

#include <vector>
#include <string>

#include "../../Helpers/ConfigFile.h"

namespace Assignment {

enum CycleMode {
    KeepCycles,
    RemoveStopCycles,
    RemoveStationCycles
};

enum DepartureTimeChoice {
    DecisionModelWithoutAdaption,
    DecisionModelWithAdaption,
    Uniform,
    Rooftop,
    DecisionModelWithBoxCox,
};

class Settings {

public:
    Settings() {}
    Settings(ConfigFile& config) {
        cycleMode = config.get("cycleMode", cycleMode);
        profilerType = config.get("profilerType", profilerType);
        useConnectionProfiler = config.get("useConnectionProfiler", useConnectionProfiler);
        connectionProfilerFrequency = config.get("connectionProfilerFrequency", connectionProfilerFrequency);
        randomSeed = config.get("randomSeed", randomSeed);
        passengerMultiplier = config.get("passengerMultiplier", passengerMultiplier);
        allowDepartureStops = config.get("allowDepartureStops", allowDepartureStops);
        transferCosts = config.get("transferCosts", transferCosts);
        walkingCosts = config.get("walkingCosts", walkingCosts);
        waitingCosts = config.get("waitingCosts", waitingCosts);
        failureCosts = config.get("failureCosts", failureCosts);
        congestionEnterCosts = config.get("congestionEnterCosts", congestionEnterCosts);
        congestionTravelCosts = config.get("congestionTravelCosts", congestionTravelCosts);
        congestionExitCosts = config.get("congestionExitCosts", congestionExitCosts);
        strandingWaitingTime = config.get("strandingWaitingTime", strandingWaitingTime);
        decisionModel = config.get("decisionModel", decisionModel);
        beta = config.get("beta", beta);
        distanceBeta = config.get("distanceBeta", distanceBeta);
        referenceTravelTime = config.get("referenceTravelTime", referenceTravelTime);
        delayTolerance = config.get("delayTolerance", delayTolerance);
        delayValue = config.get("delayValue", delayValue);
        probabilityCutoff = config.get("probabilityCutoff", probabilityCutoff);
        maxDelay = config.get("maxDelay", maxDelay);
        loadFactorCutoff = config.get("loadFactorCutoff", loadFactorCutoff);
        loadFactorSwitchPoint = config.get("loadFactorSwitchPoint", loadFactorSwitchPoint);
        loadFactorCoefficient1 = config.get("loadFactorCoefficient1", loadFactorCoefficient1);
        loadFactorCoefficient2 = config.get("loadFactorCoefficient2", loadFactorCoefficient2);
        iterationType = config.get("iterationType", iterationType);
        stepExponent = config.get("stepExponent", stepExponent);
        speedupStepSize = config.get("speedupStepSize", speedupStepSize);
        slowdownStepSize = config.get("slowdownStepSize", slowdownStepSize);
        convergenceLimit = config.get("convergenceLimit", convergenceLimit);
        demandIntervalSplitTime = config.get("demandIntervalSplitTime", demandIntervalSplitTime);
        keepDemandIntervals = config.get("keepDemandIntervals", keepDemandIntervals);
        includeIntervalBorder = config.get("includeIntervalBorder", includeIntervalBorder);
        departureTimeChoice = config.get("departureTimeChoice", departureTimeChoice);
        maxAdaptationTime = config.get("maxAdaptationTime", maxAdaptationTime);
        adaptationCost = config.get("adaptationCost", adaptationCost);
        adaptationOffset = config.get("adaptationOffset", adaptationOffset);
        adaptationBeta = config.get("adaptationBeta", adaptationBeta);
        adaptationLambda = config.get("adaptationLambda", adaptationLambda);
        randomUtilityDistribution = config.get("randomUtilityDistribution", randomUtilityDistribution);
        networkWalkingCostsMin = config.get("networkWalkingCostsMin", networkWalkingCostsMin);
        networkWalkingCostsMax = config.get("networkWalkingCostsMax", networkWalkingCostsMax);
        networkWalkingCostsA = config.get("networkWalkingCostsA", networkWalkingCostsA);
        networkWalkingCostsB = config.get("networkWalkingCostsB", networkWalkingCostsB);
        networkWaitingCostsMin = config.get("networkWaitingCostsMin", networkWaitingCostsMin);
        networkWaitingCostsMax = config.get("networkWaitingCostsMax", networkWaitingCostsMax);
        networkWaitingCostsA = config.get("networkWaitingCostsA", networkWaitingCostsA);
        networkWaitingCostsB = config.get("networkWaitingCostsB", networkWaitingCostsB);
        networkDrivingCostsMin = config.get("networkDrivingCostsMin", networkDrivingCostsMin);
        networkDrivingCostsMax = config.get("networkDrivingCostsMax", networkDrivingCostsMax);
        networkDrivingCostsA = config.get("networkDrivingCostsA", networkDrivingCostsA);
        networkDrivingCostsB = config.get("networkDrivingCostsB", networkDrivingCostsB);
        networkTransferCostsMin = config.get("networkTransferCostsMin", networkTransferCostsMin);
        networkTransferCostsMax = config.get("networkTransferCostsMax", networkTransferCostsMax);
        networkTransferCostsA = config.get("networkTransferCostsA", networkTransferCostsA);
        networkTransferCostsB = config.get("networkTransferCostsB", networkTransferCostsB);
        passengerWalkingCostsMin = config.get("passengerWalkingCostsMin", passengerWalkingCostsMin);
        passengerWalkingCostsMax = config.get("passengerWalkingCostsMax", passengerWalkingCostsMax);
        passengerWalkingCostsA = config.get("passengerWalkingCostsA", passengerWalkingCostsA);
        passengerWalkingCostsB = config.get("passengerWalkingCostsB", passengerWalkingCostsB);
        passengerWaitingCostsMin = config.get("passengerWaitingCostsMin", passengerWaitingCostsMin);
        passengerWaitingCostsMax = config.get("passengerWaitingCostsMax", passengerWaitingCostsMax);
        passengerWaitingCostsA = config.get("passengerWaitingCostsA", passengerWaitingCostsA);
        passengerWaitingCostsB = config.get("passengerWaitingCostsB", passengerWaitingCostsB);
        passengerDrivingCostsMin = config.get("passengerDrivingCostsMin", passengerDrivingCostsMin);
        passengerDrivingCostsMax = config.get("passengerDrivingCostsMax", passengerDrivingCostsMax);
        passengerDrivingCostsA = config.get("passengerDrivingCostsA", passengerDrivingCostsA);
        passengerDrivingCostsB = config.get("passengerDrivingCostsB", passengerDrivingCostsB);
        passengerTransferCostsMin = config.get("passengerTransferCostsMin", passengerTransferCostsMin);
        passengerTransferCostsMax = config.get("passengerTransferCostsMax", passengerTransferCostsMax);
        passengerTransferCostsA = config.get("passengerTransferCostsA", passengerTransferCostsA);
        passengerTransferCostsB = config.get("passengerTransferCostsB", passengerTransferCostsB);
    }

    ConfigFile toConfigFile(const std::string& fileName) const noexcept {
        ConfigFile config(fileName);
        config.set("cycleMode", cycleMode);
        config.set("profilerType", profilerType);
        config.set("useConnectionProfiler", useConnectionProfiler);
        config.set("connectionProfilerFrequency", connectionProfilerFrequency);
        config.set("randomSeed", randomSeed);
        config.set("passengerMultiplier", passengerMultiplier);
        config.set("allowDepartureStops", allowDepartureStops);
        config.set("transferCosts", transferCosts);
        config.set("walkingCosts", walkingCosts);
        config.set("waitingCosts", waitingCosts);
        config.set("failureCosts", failureCosts);
        config.set("congestionEnterCosts", congestionEnterCosts);
        config.set("congestionTravelCosts", congestionTravelCosts);
        config.set("congestionExitCosts", congestionExitCosts);
        config.set("strandingWaitingTime", strandingWaitingTime);
        config.set("decisionModel", decisionModel);
        config.set("beta", beta);
        config.set("distanceBeta", distanceBeta);
        config.set("referenceTravelTime", referenceTravelTime);
        config.set("delayTolerance", delayTolerance);
        config.set("delayValue", delayValue);
        config.set("probabilityCutoff", probabilityCutoff);
        config.set("maxDelay", maxDelay);
        config.set("loadFactorCutoff", loadFactorCutoff);
        config.set("loadFactorSwitchPoint", loadFactorSwitchPoint);
        config.set("loadFactorCoefficient1", loadFactorCoefficient1);
        config.set("loadFactorCoefficient2", loadFactorCoefficient2);
        config.set("iterationType", iterationType);
        config.set("stepExponent", stepExponent);
        config.set("speedupStepSize", speedupStepSize);
        config.set("slowdownStepSize", slowdownStepSize);
        config.set("convergenceLimit", convergenceLimit);
        config.set("demandIntervalSplitTime", demandIntervalSplitTime);
        config.set("keepDemandIntervals", keepDemandIntervals);
        config.set("includeIntervalBorder", includeIntervalBorder);
        config.set("departureTimeChoice", departureTimeChoice);
        config.set("maxAdaptationTime", maxAdaptationTime);
        config.set("adaptationCost", adaptationCost);
        config.set("adaptationOffset", adaptationOffset);
        config.set("adaptationBeta", adaptationBeta);
        config.set("adaptationLambda", adaptationLambda);
        config.set("randomUtilityDistribution", randomUtilityDistribution);
        config.set("networkWalkingCostsMin", networkWalkingCostsMin);
        config.set("networkWalkingCostsMax", networkWalkingCostsMax);
        config.set("networkWalkingCostsA", networkWalkingCostsA);
        config.set("networkWalkingCostsB", networkWalkingCostsB);
        config.set("networkWaitingCostsMin", networkWaitingCostsMin);
        config.set("networkWaitingCostsMax", networkWaitingCostsMax);
        config.set("networkWaitingCostsA", networkWaitingCostsA);
        config.set("networkWaitingCostsB", networkWaitingCostsB);
        config.set("networkDrivingCostsMin", networkDrivingCostsMin);
        config.set("networkDrivingCostsMax", networkDrivingCostsMax);
        config.set("networkDrivingCostsA", networkDrivingCostsA);
        config.set("networkDrivingCostsB", networkDrivingCostsB);
        config.set("networkTransferCostsMin", networkTransferCostsMin);
        config.set("networkTransferCostsMax", networkTransferCostsMax);
        config.set("networkTransferCostsA", networkTransferCostsA);
        config.set("networkTransferCostsB", networkTransferCostsB);
        config.set("passengerWalkingCostsMin", passengerWalkingCostsMin);
        config.set("passengerWalkingCostsMax", passengerWalkingCostsMax);
        config.set("passengerWalkingCostsA", passengerWalkingCostsA);
        config.set("passengerWalkingCostsB", passengerWalkingCostsB);
        config.set("passengerWaitingCostsMin", passengerWaitingCostsMin);
        config.set("passengerWaitingCostsMax", passengerWaitingCostsMax);
        config.set("passengerWaitingCostsA", passengerWaitingCostsA);
        config.set("passengerWaitingCostsB", passengerWaitingCostsB);
        config.set("passengerDrivingCostsMin", passengerDrivingCostsMin);
        config.set("passengerDrivingCostsMax", passengerDrivingCostsMax);
        config.set("passengerDrivingCostsA", passengerDrivingCostsA);
        config.set("passengerDrivingCostsB", passengerDrivingCostsB);
        config.set("passengerTransferCostsMin", passengerTransferCostsMin);
        config.set("passengerTransferCostsMax", passengerTransferCostsMax);
        config.set("passengerTransferCostsA", passengerTransferCostsA);
        config.set("passengerTransferCostsB", passengerTransferCostsB);
        return config;
    }

    int cycleMode{RemoveStationCycles}; // Cycle removal (CycleMode)
    int profilerType{0}; // 0 = NoProfiler, 1 = TimeProfiler, 2 = DecisionProfiler
    bool useConnectionProfiler{false}; // Runtime analysis per connection?
    int connectionProfilerFrequency{10000}; // Frequency of connection profiler output

    int randomSeed{42}; // random seed of the Monte Carlo simulation
    int passengerMultiplier{100}; // multiplier for the demand
    bool allowDepartureStops{true}; // Can demand use stops as origins?

    int transferCosts{5 * 60}; // PAT overhead for changing vehicles
    double walkingCosts{2.0}; // cost factor for the walking time in the PAT (must be >= 0, walking is counted 1 + walkingCosts times)
    double waitingCosts{0.0}; // cost factor for the waiting time in the PAT (must be >= 0, waiting is counted 1 + waitingCosts times)
    double failureCosts{5 * 60}; // PAT overhead for failing to board a vehicle
    int congestionEnterCosts{6 * 60}; // cost factor for entering a congested connection
    double congestionTravelCosts{1.2}; // cost factor for traveling with a congested connection
    int congestionExitCosts{2 * 60}; // cost factor for existing a congested connection
    int strandingWaitingTime{8 * 60 * 60}; // assumed waiting time if there is no alternative in failure-to-board case

    int decisionModel{0}; // 0 = Linear, 1 = Logit, 2 = Kirchhoff, 3 = RelativeLogit, 4 = Optimal
    double beta{1.0}; // Adjustment parameter for Logit & Kirchhoff
    double distanceBeta{1.0}; // Adjust parameter for DistanceLogit
    int referenceTravelTime{30 * 60}; // Typical travel time, used to calibrate DistanceLogit
    int delayTolerance{5 * 60}; // maximum difference a journey PAT can have to the optimal PAT in order to be considered for assignment of passengers
    //TODO: Can this be removed?
    int delayValue{5 * 60}; // Linear: PAT overhead for non-optimal journeys
    double probabilityCutoff{0.01}; //Sequential: Minimum non-zero probability that an option may have of being chosen over the best option

    int maxDelay{0}; // max delay of vehicles in the MEAT model

    // Capacity assignment
    double loadFactorCutoff{0.3}; // Maximum load value at which congestion is not punished
    double loadFactorSwitchPoint{0.5}; // Load value at which congestion punishment switches from quadratic to exponential
    double loadFactorCoefficient1{5}; // Coefficient for quadratic congestion punishment
    double loadFactorCoefficient2{2}; // Coefficient for exponential congestion punishment
    int iterationType{0}; // 0 = Successive averages, 1 = Successive weighted averages, 2 = Self-regulating averages
    double stepExponent{1}; // Weighted averages: Exponent of step size increase
    double speedupStepSize{1.8}; //Self-regulating averages: Step size increase when iterate diverges
    double slowdownStepSize{0.3}; //Self-regulating averages: Step size increase when iterate converges
    double convergenceLimit{0.1}; // Maximum allowed load difference between current iteration and average of previous iterations

    int demandIntervalSplitTime{86400}; // Time interval size for discretization of the demand input time intervals (negative value indicates no discretization)
    bool keepDemandIntervals{true}; // false = collapse demand departure time intervals to their minimal value, true = keep full intervals
    bool includeIntervalBorder{false}; // true = intervals before discretization are interpreted as (min <= x <= max), false intervals before discretization are interpreted as (min <= x < max)

    int departureTimeChoice{DecisionModelWithoutAdaption}; // Handling of departure times in demand (DepartureTimeChoice)
    int maxAdaptationTime{0}; // Maximum amount that passengers are willing to adjust their departure time
    double adaptationCost{2.0}; // DecisionModelWithAdaptation/Rooftop: Cost factor for adjusting departure time
    int adaptationOffset{0}; // DecisionModelWithAdaptation: Maximum amount of adaptation that is allowed without incurring costs in
    double adaptationBeta{0.1}; // DecisionModelWithBoxCox: Beta value for Box-Cox transformation
    double adaptationLambda{2.0}; // DecisionModelWithBoxCox: Lambda value for Box-Cox transformation

    // MRUM
    int randomUtilityDistribution{0}; // 0 = Logit-normal, 1 = PERT, 2 = Kumaraswamy, 3 = Beta
    double networkWalkingCostsMin{1.0};
    double networkWalkingCostsMax{2.0};
    double networkWalkingCostsA{0.5};
    double networkWalkingCostsB{1.0};
    double networkWaitingCostsMin{1.0};
    double networkWaitingCostsMax{2.0};
    double networkWaitingCostsA{0.5};
    double networkWaitingCostsB{1.0};
    double networkDrivingCostsMin{1.0};
    double networkDrivingCostsMax{2.0};
    double networkDrivingCostsA{0.5};
    double networkDrivingCostsB{1.0};
    double networkTransferCostsMin{1.0};
    double networkTransferCostsMax{2.0};
    double networkTransferCostsA{0.5};
    double networkTransferCostsB{1.0};
    double passengerWalkingCostsMin{1.0};
    double passengerWalkingCostsMax{2.0};
    double passengerWalkingCostsA{0.5};
    double passengerWalkingCostsB{1.0};
    double passengerWaitingCostsMin{1.0};
    double passengerWaitingCostsMax{2.0};
    double passengerWaitingCostsA{0.5};
    double passengerWaitingCostsB{1.0};
    double passengerDrivingCostsMin{1.0};
    double passengerDrivingCostsMax{2.0};
    double passengerDrivingCostsA{0.5};
    double passengerDrivingCostsB{1.0};
    double passengerTransferCostsMin{1.0};
    double passengerTransferCostsMax{2.0};
    double passengerTransferCostsA{0.5};
    double passengerTransferCostsB{1.0};

};

}
