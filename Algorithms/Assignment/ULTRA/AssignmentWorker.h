#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "../../../DataStructures/Assignment/AssignmentData.h"
#include "../../../DataStructures/Assignment/ChoiceSet.h"
#include "../../../DataStructures/Assignment/GroupData.h"
#include "../../../DataStructures/Assignment/GroupTrackingData.h"
#include "../../../DataStructures/Assignment/Settings.h"
#include "../../../DataStructures/CSA/Data.h"
#include "../../../DataStructures/Container/Heap.h"
#include "../../../DataStructures/Demand/AccumulatedVertexDemand.h"
#include "../../../DataStructures/Graph/Graph.h"

#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/Vector/Vector.h"

#include "../../CH/CH.h"
#include "../../CH/Query/BucketQuery.h"

#include "../CycleRemoval.h"
#include "../PassengerDistribution.h"
#include "../Profiler.h"

#include "ComputePATs.h"
#include "StopLabel.h"

namespace Assignment::ULTRA {

template<typename DECISION_MODEL, typename PROFILER>
class AssignmentWorker {

public:
    using DecisionModel = DECISION_MODEL;
    using Profiler = PROFILER;
    using Type = AssignmentWorker<DecisionModel, Profiler>;
    using ConnectionLabel = typename ComputePATs<NoPATProfiler>::ConnectionLabel;

    struct DistanceLabel {
        DistanceLabel(const StopId stop, const int distance) :
            stop(stop),
            distance(distance) {
        }
        StopId stop;
        int distance;
    };

public:
    AssignmentWorker(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings, const DecisionModel& decisionModel, const std::vector<std::vector<DistanceLabel>>& zoneToStop, const CH::BucketQuery<>& bucketQuery) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        decisionModel(decisionModel),
        zoneToStop(zoneToStop),
        bucketQuery(bucketQuery),
        patComputation(data, reverseGraph, this->bucketQuery),
        groupTrackingData(data.numberOfStops(), data.numberOfTrips()),
        assignmentData(data.numberOfConnections()),
        cycleRemoval(data, settings.cycleMode, assignmentData) {
        profiler.initialize(data);
    }

    inline void run(const Vertex destinationVertex, std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        Assert(!demand.empty(), "Demand for destination vertex " << destinationVertex << " is empty!");
        Assert(data.isStop(destinationVertex) || reverseGraph.outDegree(destinationVertex) > 0, "Destination vertex " << destinationVertex << " is isolated!");
        profiler.startAssignmentForDestination(destinationVertex);

        sort(demand, [](const AccumulatedVertexDemand::Entry& a, const AccumulatedVertexDemand::Entry& b){return a.earliestDepartureTime < b.earliestDepartureTime;});

        profiler.startPATComputation();
        patComputation.run(destinationVertex, settings.maxDelay, settings.transferCosts, settings.walkingCosts, settings.waitingCosts, settings.walkingCosts + 1.0, demand.front().earliestDepartureTime);
        profiler.donePATComputation();

        profiler.startInitialWalking();
        initializeAssignment(demand);
        profiler.doneInitialWalking();
        profiler.startAssignment();
        for (const ConnectionId i : data.connectionIds()) {
            profiler.assignConnection(i);
            processConnection(i);
        }
        profiler.doneAssignment();

        profiler.doneAssignmentForDestination(destinationVertex);
    }

    inline void runCycleRemoval() noexcept {
        profiler.startCycleElimination();
        cycleRemoval.run();
        profiler.doneCycleElimination();
    }

    inline const AssignmentData& getAssignmentData() const noexcept {
        return assignmentData;
    }

    inline u_int64_t getRemovedCycleConnections() const noexcept {
        return cycleRemoval.getRemovedCycleConnections();
    }

    inline u_int64_t getRemovedCycles() const noexcept {
        return cycleRemoval.getRemovedCycles();
    }

    inline Profiler& getProfiler() noexcept {
        return profiler;
    }

private:
    inline void initializeAssignment(std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        groupTrackingData.validate();
        for (const StopId stop : data.stops()) {
            patComputation.getProfile(stop).resetIndex();
        }
        walkToInitialStops(demand);
        for (const StopId stop : data.stops()) {
            patComputation.getProfile(stop).resetIndex();
            sort(groupTrackingData.groupsOriginatingAtStop[stop], [](const GroupArrivalLabel& a, const GroupArrivalLabel& b){return a.arrivalTime > b.arrivalTime;});
        }
    }

    inline void walkToInitialStops(std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        sort(demand, [](const AccumulatedVertexDemand::Entry& a, const AccumulatedVertexDemand::Entry& b){return a.originVertex < b.originVertex;});
        switch(settings.departureTimeChoice) {
            case DecisionModelWithoutAdaption:
                walkToInitialStops<DecisionModelWithoutAdaption>(demand);
                break;
            case DecisionModelWithAdaption:
                walkToInitialStops<DecisionModelWithAdaption>(demand);
                break;
            case Uniform:
                walkToInitialStops<Uniform>(demand);
                break;
            case Rooftop:
                walkToInitialStops<Rooftop>(demand);
                break;
            case DecisionModelWithBoxCox:
                walkToInitialStops<DecisionModelWithBoxCox>(demand);
                break;
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline void walkToInitialStops(const std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        const std::vector<int> targetDistance = bucketQuery.getBackwardDistance();
        for (const AccumulatedVertexDemand::Entry& demandEntry : demand) {
            Assert(demandEntry.originVertex != demandEntry.destinationVertex, "Origin and destination vertex of demand are identical (" << demandEntry.originVertex << ")!");
            Assert(settings.allowDepartureStops || !data.isStop(demandEntry.originVertex), "Demand is originating from a stop (" << demandEntry.originVertex << ")!");
            Assert(data.isStop(demandEntry.originVertex) || data.transferGraph.outDegree(demandEntry.originVertex) > 0, "Origin vertex " << demandEntry.originVertex << " of demand is isolated!");
            ChoiceSet choiceSet;
            int minPAT;
            if (targetDistance[demandEntry.originVertex] < INFTY) {
                const int pat = targetDistance[demandEntry.originVertex] * (settings.walkingCosts + 1);
                choiceSet.addChoice(StopId(demandEntry.destinationVertex), demandEntry.earliestDepartureTime, pat);
                Assert(targetDistance[demandEntry.originVertex] <= pat, "Overflow!");
            }
            for (const DistanceLabel& label : zoneToStop[demandEntry.originVertex]) {
                const int walkingPTT = label.distance * (settings.walkingCosts + 1);
                Assert(label.distance <= walkingPTT, "Overflow!");
                if (walkingPTT - settings.delayTolerance > minPAT) break;
                if (walkingPTT + patComputation.getProfile(label.stop).getMinimumPTT() - settings.delayTolerance > minPAT) continue;
                evaluateInitialStop<DEPARTURE_TIME_CHOICE>(demandEntry, label.stop, label.distance, choiceSet, minPAT);
            }
            const GroupId originalGroup = assignmentData.createNewGroup(demandEntry, settings.passengerMultiplier);
            if (choiceSet.empty()) {
                assignmentData.unassignedGroups.emplace_back(originalGroup);
            } else if (choiceSet.size() == 1) {
                if (data.isStop(Vertex(choiceSet.options[0]))) {
                    groupTrackingData.groupsOriginatingAtStop[choiceSet.options[0]].emplace_back(originalGroup, choiceSet.departureTimes[0]);
                } else {
                    Assert(Vertex(choiceSet.options[0]) == demandEntry.destinationVertex, "Passengers a walking to a vertex that is not their destination!");
                    assignmentData.directWalkingGroups.emplace_back(originalGroup);
                }
            } else {
                profiler.distributePassengersPATs(choiceSet.pats, choiceSet.departureTimes);
                std::vector<int> distribution;
                if (DEPARTURE_TIME_CHOICE == Rooftop) {
                    distribution = choiceSet.rooftopDistribution(demandEntry, settings.adaptationCost);
                } else {
                    distribution = decisionModel.distribution(choiceSet.pats);
                }
                profiler.distributePassengersProbabilities(distribution);
                const std::vector<size_t> groupSizes = getGroupSizes(distribution, demandEntry.numberOfPassengers * settings.passengerMultiplier, random);
                profiler.distributePassengersSizes(groupSizes);
                size_t originalGroupIndex = choiceSet.size();
                for (size_t i = 0; i < choiceSet.size(); i++) {
                    if (groupSizes[i] == 0) continue;
                    GroupId group;
                    if (originalGroupIndex < i) {
                        group = assignmentData.splitGroup(originalGroup, groupSizes[i]);
                    } else {
                        group = originalGroup;
                        originalGroupIndex = i;
                    }
                    if (data.isStop(Vertex(choiceSet.options[i]))) {
                        groupTrackingData.groupsOriginatingAtStop[choiceSet.options[i]].emplace_back(group, choiceSet.departureTimes[i]);
                    } else {
                        Assert(Vertex(choiceSet.options[i]) == demandEntry.destinationVertex, "Passengers a walking to a vertex that is not their destination!");
                        assignmentData.directWalkingGroups.emplace_back(group);
                    }
                }
                Assert(originalGroupIndex < choiceSet.size(), "No groups have been assigned!");
                Assert(assignmentData.groups[originalGroup].groupSize == groupSizes[originalGroupIndex], "Original group has wrong size (size should be: " << groupSizes[originalGroupIndex] << ", size is: " << assignmentData.groups[originalGroup].groupSize << ")!");
            }
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline void evaluateInitialStop(const AccumulatedVertexDemand::Entry& demandEntry, const StopId stop, const int transferTime, ChoiceSet& choiceSet, int& minPAT) noexcept {
        int departureTime = demandEntry.earliestDepartureTime - getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        Assert(demandEntry.earliestDepartureTime <= departureTime, "Overflow!");
        const int latestDepartureTime = demandEntry.latestDepartureTime + getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        patComputation.getProfile(stop).findIndexFast(departureTime);
        while (departureTime <= latestDepartureTime) {
            const int value = patComputation.getProfile(stop).evaluateWithoutWaiting(departureTime, settings.waitingCosts);
            if constexpr ((DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) | (DEPARTURE_TIME_CHOICE == DecisionModelWithBoxCox)) {
                if (departureTime > latestDepartureTime) return;
            }
            if (value >= INFTY) return;
            const int pat = value - departureTime + (transferTime * (1 + settings.walkingCosts)) + getAdaptationCost<DEPARTURE_TIME_CHOICE>(demandEntry, departureTime - transferTime);
            choiceSet.addChoice(stop, departureTime, pat);
            minPAT = std::min(minPAT, pat);
            departureTime++;

        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline int getMaxAdaptationTime() const noexcept {
        if constexpr ((DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) | (DEPARTURE_TIME_CHOICE == DecisionModelWithBoxCox)) {
            return settings.maxAdaptationTime;
        } else {
            return 0;
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline int getAdaptationCost(const AccumulatedVertexDemand::Entry& demandEntry, const int departureTime) const noexcept {
        if constexpr ((DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) | (DEPARTURE_TIME_CHOICE == DecisionModelWithBoxCox)) {
            const int adaptationTime = std::max(0, std::max(demandEntry.earliestDepartureTime - departureTime, departureTime - demandEntry.latestDepartureTime));
            int adaptationCost;
            if constexpr (DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) {
                adaptationCost = std::max(0, adaptationTime - settings.adaptationOffset) * settings.adaptationCost;
            } else {
                adaptationCost = 60 * settings.adaptationBeta * (std::pow(adaptationTime / 60, settings.adaptationLambda) - 1) / settings.adaptationLambda;
            }
            return adaptationCost;
        } else {
            return 0;
        }
    }

    inline void processConnection(const ConnectionId i) noexcept {
        const CSA::Connection& connection = data.connections[i];
        groupTrackingData.processOriginatingGroups(connection);
        groupTrackingData.processWalkingGroups(connection);
        if (groupTrackingData.groupsWaitingAtStop[connection.departureStopId].empty() && groupTrackingData.groupsInTrip[connection.tripId].empty()) return;
        const ConnectionLabel& label = patComputation.connectionLabel(i);
        const double targetPAT = patComputation.targetPAT(connection);
        const double hopOffPAT = std::min(targetPAT, label.transferPAT);
        const double hopOnPAT = std::min(hopOffPAT, label.tripPAT);
        Assert(hopOnPAT >= connection.arrivalTime, "TripPAT is to low (Connection: " << i << ", Trip: " << connection.tripId << ", PAT: " << label.tripPAT << ", ArrivalTime: " << connection.arrivalTime << ")!");
        moveGroups(groupTrackingData.groupsWaitingAtStop[connection.departureStopId], groupTrackingData.groupsInTrip[connection.tripId], label.skipPAT, hopOnPAT, "skip", "board");
        for (const GroupId group : groupTrackingData.groupsInTrip[connection.tripId]) {
            Assert(group < assignmentData.connectionsPerGroup.size(), "Group " << group << " is out of bounds (0, " << assignmentData.connectionsPerGroup.size() << ")");
            assignmentData.connectionsPerGroup[group].emplace_back(i);
        }
        GroupList groupsHoppingOff;
        moveGroups(groupTrackingData.groupsInTrip[connection.tripId], groupsHoppingOff, label.tripPAT, hopOffPAT, "continue", "alight");
        if (groupsHoppingOff.empty()) return;
        moveGroups(groupsHoppingOff, groupTrackingData.groupsAtTarget, label.transferPAT, targetPAT, "walk", "target");
        if (groupsHoppingOff.empty()) return;
        Assert(label.transferPAT - targetPAT <= settings.delayTolerance, "Groups are not walking straight to the target (transferPAT = " << label.transferPAT << ", targetPAT= " << targetPAT << ")!");
        distributePassengers(connection, groupsHoppingOff);
    }

    template<typename TO_PASSENGER_LIST>
    inline void moveGroups(GroupList& from, TO_PASSENGER_LIST& to, const double fromPAT, const double toPAT, const std::string& fromName, const std::string& toName) noexcept {
        if (from.empty()) return;
        profiler.moveGroups(fromName, toName);
        profiler.moveGroupsPATs(fromPAT, toPAT);
        const std::array<int, 3> values = decisionModel.distribution(fromPAT, toPAT);
        profiler.moveGroupsProbabilities(values);
        for (size_t i = 0; i < from.size(); i++) {
            const std::array<int, 2> groupSizes = getGroupSizes(values, assignmentData.groups[from[i]].groupSize, random);
            profiler.moveGroupsSizes(groupSizes);
            if (groupSizes[0] == 0) {
                to.emplace_back(from[i]);
                from[i] = from.back();
                from.pop_back();
                i--;
            } else if (groupSizes[1] != 0) {
                to.emplace_back(assignmentData.splitGroup(from[i], groupSizes[1]));
            }
        }
    }

    inline void distributePassengers(const CSA::Connection& connection, GroupList& groupList) noexcept {
        if (data.transferGraph.outDegree(connection.arrivalStopId) == 0) {
            groupTrackingData.groupsWalkingToStop[connection.arrivalStopId].push_back(GroupArrivalLabel(groupList, connection.arrivalTime + data.minTransferTime(connection.arrivalStopId)));
        } else {
            walkToNextStop(connection.arrivalStopId, groupList, connection.arrivalTime);
        }
    }

    inline void walkToNextStop(const StopId from, GroupList& groupList, const int time) noexcept {
        Assert(data.transferGraph.outDegree(from) != 0, "Cannot walk to next stop from a vertex with out degree zero!");
        ChoiceSet choiceSet;
        for (const Edge edge : data.transferGraph.edgesFrom(from)) {
            const Vertex intermediateStop = data.transferGraph.get(ToVertex, edge);
            if (!data.isStop(intermediateStop)) continue;
            evaluateIntermediateStop(StopId(intermediateStop), time, data.transferGraph.get(TravelTime, edge), choiceSet);
        }
        if (data.isStop(from)) {
            evaluateIntermediateStop(StopId(from), time, 0, choiceSet);
        }
        Assert(!choiceSet.empty(), "" << groupList.size() << " groups arrived at stop " << from << " but have nowhere to go!");
        if (choiceSet.size() == 1) {
            groupTrackingData.groupsWalkingToStop[choiceSet.options[0]].emplace_back(groupList, choiceSet.departureTimes[0]);
        } else {
            profiler.distributePassengersPATs(choiceSet.pats, choiceSet.departureTimes);
            const std::vector<int> distribution = decisionModel.distribution(choiceSet.pats);
            profiler.distributePassengersProbabilities(distribution);
            std::vector<GroupList> groupListsByIndex(choiceSet.size());
            for (size_t i = 0; i < groupList.size(); i++) {
                const std::vector<size_t> groupSizes = getGroupSizes(distribution, assignmentData.groups[groupList[i]].groupSize, random);
                profiler.distributePassengersSizes(groupSizes);
                bool movedOriginalGroup = false;
                for (size_t j = 0; j < groupSizes.size(); j++) {
                    if (groupSizes[j] == 0) continue;
                    GroupId group;
                    if (movedOriginalGroup) {
                        group = assignmentData.splitGroup(groupList[i], groupSizes[j]);
                    } else {
                        group = groupList[i];
                        movedOriginalGroup = true;
                    }
                    groupListsByIndex[j].emplace_back(group);
                }
                Assert(movedOriginalGroup, "Group has not moved to the next stop (Group: " << assignmentData.groups[groupList[i]] << ")");
            }
            for (size_t i = 0; i < groupListsByIndex.size(); i++) {
                if (groupListsByIndex[i].empty()) continue;
                groupTrackingData.groupsWalkingToStop[choiceSet.options[i]].emplace_back(groupListsByIndex[i], choiceSet.departureTimes[i]);
            }
        }
    }

    inline void evaluateIntermediateStop(const StopId stop, const int time, const int transferTime, ChoiceSet& choiceSet) noexcept {
        const int bufferTime = data.minTransferTime(stop);
        const int departureTime = time + transferTime + bufferTime;
        const int value = patComputation.getProfile(stop).evaluateWithWaitingCosts(departureTime - bufferTime, settings.waitingCosts);
        if (value >= INFTY) return;
        const int pat = value + (transferTime * settings.walkingCosts);
        choiceSet.addChoice(stop, departureTime, pat);
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;
    const DecisionModel& decisionModel;
    const std::vector<std::vector<DistanceLabel>>& zoneToStop;

    CH::BucketQuery<> bucketQuery;
    ComputePATs<NoPATProfiler> patComputation;

    GroupTrackingData groupTrackingData;
    AssignmentData assignmentData;

    CycleRemoval cycleRemoval;
    Random random;
    Profiler profiler;

};

}
