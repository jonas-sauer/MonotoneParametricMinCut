#pragma once

#include <vector>

#include "../../DataStructures/Assignment/AssignmentData.h"
#include "../../DataStructures/Assignment/ChoiceSet.h"
#include "../../DataStructures/Assignment/GroupData.h"
#include "../../DataStructures/Assignment/GroupTrackingData.h"
#include "../../DataStructures/Assignment/Settings.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/AccumulatedVertexDemand.h"
#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/Helpers.h"
#include "../../Helpers/Vector/Vector.h"

#include "ComputePATs.h"
#include "CycleRemoval.h"
#include "PassengerDistribution.h"
#include "Profiler.h"

namespace Assignment {

template<typename DECISION_MODEL, typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES = false>
class AssignmentWorker {

public:
    using DecisionModel = DECISION_MODEL;
    using Profiler = PROFILER;
    constexpr static inline bool UseTransferBufferTimes = USE_TRANSFER_BUFFER_TIMES;
    using Type = AssignmentWorker<DecisionModel, Profiler, UseTransferBufferTimes>;

    using PATComputationType = ComputePATs<NoPATProfiler, UseTransferBufferTimes>;
    using ConnectionLabel = typename PATComputationType::ConnectionLabel;

    struct ProfileReader {
    public:
        ProfileReader() : i(-1), data(nullptr) {}

        ProfileReader(const Profile& profile) : i(profile.size() - 1), data(&profile) {}

        inline void initialize(const Profile& profile) noexcept {
            i = profile.size() - 1;
            data = &profile;
        }

        inline void reset() noexcept {
            i = data->size() - 1;
        }

        inline const ProfileEntry& findEntry(const int time) noexcept {
            while (i + 1 < data->size() && (*data)[i + 1].departureTime >= time) i++;
            AssertMsg((i + 1 >= data->size()) || ((*data)[i + 1].departureTime < time), "Profile is not scanned monotonously (current time: " << time << " previous time: " << (*data)[i + 1].departureTime << ")!");
            while ((*data)[i].departureTime < time) {
                i--;
            }
            return (*data)[i];
        }

    private:
        size_t i;
        const Profile* data;
    };

public:
    AssignmentWorker(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings, const DecisionModel& decisionModel) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        decisionModel(decisionModel),
        patComputation(data, reverseGraph),
        profiles(data.numberOfStops()),
        groupTrackingData(data.numberOfStops(), data.numberOfTrips()),
        assignmentData(data.numberOfConnections()),
        cycleRemoval(data, settings.cycleMode, assignmentData) {
        profiler.initialize(data);
    }

    inline void run(const Vertex destinationVertex, std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        AssertMsg(!demand.empty(), "Demand for destination vertex " << destinationVertex << " is empty!");
        AssertMsg(data.isStop(destinationVertex) || reverseGraph.outDegree(destinationVertex) > 0, "Destination vertex " << destinationVertex << " is isolated!");
        profiler.startAssignmentForDestination(destinationVertex);

        sort(demand, [](const AccumulatedVertexDemand::Entry& a, const AccumulatedVertexDemand::Entry& b){return a.earliestDepartureTime < b.earliestDepartureTime;});

        profiler.startPATComputation();
        patComputation.run(destinationVertex, settings.maxDelay, settings.transferCosts, settings.walkingCosts, settings.waitingCosts);
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
    inline void initializeAssignment(const std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        groupTrackingData.validate();
        for (const StopId stop : data.stops()) {
            profiles[stop].initialize(patComputation.getProfile(stop));
        }
        walkToInitialStops(demand);
        for (const StopId stop : data.stops()) {
            profiles[stop].reset();
            sort(groupTrackingData.groupsOriginatingAtStop[stop], [](const GroupArrivalLabel& a, const GroupArrivalLabel& b){return a.arrivalTime > b.arrivalTime;});
        }
    }

    inline void walkToInitialStops(const std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
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
        for (const AccumulatedVertexDemand::Entry& demandEntry : demand) {
            AssertMsg(demandEntry.originVertex != demandEntry.destinationVertex, "Origin and destination vertex of demand are identical (" << demandEntry.originVertex << ")!");
            AssertMsg(settings.allowDepartureStops || !data.isStop(demandEntry.originVertex), "Demand is originating from a stop (" << demandEntry.originVertex << ")!");
            AssertMsg(data.isStop(demandEntry.originVertex) || data.transferGraph.outDegree(demandEntry.originVertex) > 0, "Origin vertex " << demandEntry.originVertex << " of demand is isolated!");
            ChoiceSet choiceSet = collectInitialWalkingChoices<DEPARTURE_TIME_CHOICE>(demandEntry);
            const GroupId originalGroup = assignmentData.createNewGroup(demandEntry, settings.passengerMultiplier);
            if (choiceSet.empty()) {
                assignmentData.unassignedGroups.emplace_back(originalGroup);
            } else if (choiceSet.size() == 1) {
                groupTrackingData.groupsOriginatingAtStop[choiceSet.options[0]].emplace_back(originalGroup, choiceSet.departureTimes[0]);
            } else {
                profiler.distributePassengersPATs(choiceSet.pats, choiceSet.departureTimes);
                std::vector<int> distribution;
                if constexpr (DEPARTURE_TIME_CHOICE == Rooftop) {
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
                    groupTrackingData.groupsOriginatingAtStop[choiceSet.options[i]].emplace_back(group, choiceSet.departureTimes[i]);
                    if (choiceSet.options[i] == demandEntry.destinationVertex) {
                        assignmentData.directWalkingGroups.emplace_back(group);
                    }
                }
                AssertMsg(originalGroupIndex < choiceSet.size(), "No groups have been assigned!");
                AssertMsg(assignmentData.groups[originalGroup].groupSize == groupSizes[originalGroupIndex], "Original group has wrong size (size should be: " << groupSizes[originalGroupIndex] << ", size is: " << assignmentData.groups[originalGroup].groupSize << ")!");
            }
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline ChoiceSet collectInitialWalkingChoices(const AccumulatedVertexDemand::Entry& demandEntry) noexcept {
        ChoiceSet choiceSet;
        bool foundInitialStop = false;
        for (const Edge edge : data.transferGraph.edgesFrom(demandEntry.originVertex)) {
            const Vertex initialStop = data.transferGraph.get(ToVertex, edge);
            if (!data.isStop(initialStop)) continue;
            evaluateInitialStop<DEPARTURE_TIME_CHOICE>(demandEntry, initialStop, data.transferGraph.get(TravelTime, edge), choiceSet);
            foundInitialStop = true;
        }
        if (data.isStop(demandEntry.originVertex)) {
            evaluateInitialStop<DEPARTURE_TIME_CHOICE>(demandEntry, demandEntry.originVertex, 0, choiceSet);
            foundInitialStop = true;
        }
        AssertMsg(foundInitialStop, "Demand is originating from a vertex that is not connected to a stop (" << demandEntry.originVertex << ")!");
        return choiceSet;
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline void evaluateInitialStop(const AccumulatedVertexDemand::Entry& demandEntry, const Vertex stop, const int transferTime, ChoiceSet& choiceSet) noexcept {
        int departureTime = demandEntry.earliestDepartureTime - getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        const int latestDepartureTime = demandEntry.latestDepartureTime + getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        while (departureTime <= latestDepartureTime) {
            const ProfileEntry& entry = profiles[stop].findEntry(departureTime);
            departureTime = entry.departureTime;
            if constexpr ((DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) | (DEPARTURE_TIME_CHOICE == DecisionModelWithBoxCox)) {
                if (departureTime > latestDepartureTime) return;
            }
            const int value = entry.evaluate(departureTime, settings.waitingCosts);
            if (value >= INFTY) return;
            const int pat = value - departureTime + (transferTime * (1 + settings.walkingCosts)) + getAdaptationCost<DEPARTURE_TIME_CHOICE>(demandEntry, departureTime - transferTime);
            choiceSet.addChoice(StopId(stop), departureTime, pat);
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
        const ConnectionLabel& label = patComputation.connectionLabel(i);
        const double targetPAT = patComputation.targetPAT(connection);
        const double hopOffPAT = std::min(targetPAT, label.transferPAT);
        const double hopOnPAT = std::min(hopOffPAT, label.tripPAT);
        moveGroups(groupTrackingData.groupsWaitingAtStop[connection.departureStopId], groupTrackingData.groupsInTrip[connection.tripId], label.skipPAT, hopOnPAT, "skip", "board");
        for (const GroupId group : groupTrackingData.groupsInTrip[connection.tripId]) {
            AssertMsg(group < assignmentData.connectionsPerGroup.size(), "Group " << group << " is out of bounds (0, " << assignmentData.connectionsPerGroup.size() << ")");
            assignmentData.connectionsPerGroup[group].emplace_back(i);
        }
        GroupList groupsHoppingOff;
        moveGroups(groupTrackingData.groupsInTrip[connection.tripId], groupsHoppingOff, label.tripPAT, hopOffPAT, "continue", "alight");
        if (groupsHoppingOff.empty()) return;
        moveGroups(groupsHoppingOff, groupTrackingData.groupsAtTarget, label.transferPAT, targetPAT, "walk", "target");
        if (groupsHoppingOff.empty()) return;
        AssertMsg(label.transferPAT - targetPAT <= settings.delayTolerance, "Groups are not walking straight to the target (transferPAT = " << label.transferPAT << ", targetPAT = " << targetPAT << ")!");
        walkToNextStop(connection.arrivalStopId, groupsHoppingOff, connection.arrivalTime);
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

    inline void walkToNextStop(const StopId from, GroupList& groupList, const int time) noexcept {
        if (data.transferGraph.outDegree(from) == 0) {
            groupTrackingData.groupsWalkingToStop[from].emplace_back(groupList, time + data.minTransferTime(from));
            return;
        }
        ChoiceSet choiceSet = collectIntermediateWalkingChoices(from, time);
        AssertMsg(!choiceSet.empty(), "" << groupList.size() << " groups arrived at stop " << from << " but have nowhere to go!");
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
                AssertMsg(movedOriginalGroup, "Group has not moved to the next stop (Group: " << assignmentData.groups[groupList[i]] << ")");
            }
            for (size_t i = 0; i < groupListsByIndex.size(); i++) {
                if (groupListsByIndex[i].empty()) continue;
                groupTrackingData.groupsWalkingToStop[choiceSet.options[i]].emplace_back(groupListsByIndex[i], choiceSet.departureTimes[i]);
            }
        }
    }

    inline ChoiceSet collectIntermediateWalkingChoices(const StopId from, const int time) noexcept {
        ChoiceSet choiceSet;
        for (const Edge edge : data.transferGraph.edgesFrom(from)) {
            const Vertex intermediateStop = data.transferGraph.get(ToVertex, edge);
            if (!data.isStop(intermediateStop)) continue;
            const int bufferTime = UseTransferBufferTimes ? data.minTransferTime(StopId(intermediateStop)) : 0;
            evaluateIntermediateStop(intermediateStop, time, data.transferGraph.get(TravelTime, edge), bufferTime, choiceSet);
        }
        if (data.isStop(from)) {
            evaluateIntermediateStop(from, time, 0, data.minTransferTime(from), choiceSet);
        }
        return choiceSet;
    }

    inline void evaluateIntermediateStop(const Vertex stop, const int time, const int transferTime, const int bufferTime, ChoiceSet& choiceSet) noexcept {
        const int departureTime = time + transferTime + bufferTime;
        const ProfileEntry& entry = profiles[stop].findEntry(departureTime);
        const int value = entry.evaluate(departureTime - bufferTime, settings.waitingCosts);
        if (value >= INFTY) return;
        const int pat = value + (transferTime * settings.walkingCosts);
        choiceSet.addChoice(StopId(stop), departureTime, pat);
    }

private:
    //Input
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;
    const DecisionModel& decisionModel;

    //PAT computation
    PATComputationType patComputation;
    std::vector<ProfileReader> profiles;

    GroupTrackingData groupTrackingData;
    AssignmentData assignmentData;

    CycleRemoval cycleRemoval;
    Random random;
    Profiler profiler;
};

}
