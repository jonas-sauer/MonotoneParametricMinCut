#pragma once

#include <vector>

#include "../../DataStructures/Assignment/AssignmentData.h"
#include "../../DataStructures/Assignment/ChoiceSet.h"
#include "../../DataStructures/Assignment/GroupData.h"
#include "../../DataStructures/Assignment/GroupTrackingData.h"
#include "../../DataStructures/Assignment/Settings.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/TimeExpandedNetwork.h"
#include "../../DataStructures/Demand/AccumulatedVertexDemand.h"
#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/Helpers.h"

#include "ComputeSequentialPATs.h"
#include "PassengerDistribution.h"
#include "Profiler.h"

namespace Assignment {

//TODO: Only explore efficient arcs
//TODO: Round all PATs to integers?
template<typename PROFILER>
class SequentialAssignmentWorker {

public:
    using Profiler = PROFILER;
    using Type = SequentialAssignmentWorker<Profiler>;

    using PATComputationType = ComputeSequentialPATs;
    using NodeData = typename PATComputationType::NodeData;
    using ConnectionLabel = typename PATComputationType::ConnectionLabel;
    using DecisionModel = typename PATComputationType::DecisionModel;

    struct ConnectionReader {
        ConnectionReader(const CSA::Data& data) : data(data), connectionsPerStop(data.numberOfStops()) {
            for (const ConnectionId i : data.connectionIds()) {
                const CSA::Connection& connection = data.connections[i];
                connectionsPerStop[connection.departureStopId].emplace_back(i);
            }
        }

        //TODO: Find range faster?
        template<typename FUNCTION>
        inline void scanConnectionsInRange(const Vertex stop, const int begin, const int end, const FUNCTION& function) const noexcept {
            for (size_t i = 0; i < connectionsPerStop[stop].size(); i++) {
                const ConnectionId connectionId = connectionsPerStop[stop][i];
                const CSA::Connection& connection = data.connections[connectionId];
                if (connection.departureTime < begin) continue;
                if (connection.departureTime > end) return;
                function(connectionId);
            }
        }

    private:
        const CSA::Data& data;
        std::vector<std::vector<ConnectionId>> connectionsPerStop;
    };

public:
    SequentialAssignmentWorker(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const CSA::TimeExpandedNetwork& timeExpandedNetwork, const Settings& settings, const DecisionModel& decisionModel) :
        data(data),
        reverseGraph(reverseGraph),
        timeExpandedNetwork(timeExpandedNetwork),
        settings(settings),
        decisionModel(decisionModel),
        patComputation(data, reverseGraph, timeExpandedNetwork, settings, decisionModel),
        connectionReader(data),
        groupTrackingData(data.numberOfConnections(), data.numberOfTrips()),
        assignmentData(data.numberOfConnections()) {
        profiler.initialize(data);
    }

    inline void run(const Vertex destinationVertex, std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        AssertMsg(!demand.empty(), "Demand for destination vertex " << destinationVertex << " is empty!");
        AssertMsg(data.isStop(destinationVertex) || reverseGraph.outDegree(destinationVertex) > 0, "Destination vertex " << destinationVertex << " is isolated!");
        profiler.startAssignmentForDestination(destinationVertex);

        sort(demand, [](const AccumulatedVertexDemand::Entry& a, const AccumulatedVertexDemand::Entry& b){return a.earliestDepartureTime < b.earliestDepartureTime;});

        profiler.startPATComputation();
        patComputation.run(destinationVertex);
        profiler.donePATComputation();

        profiler.startInitialWalking();
        groupTrackingData.validate();
        walkToInitialStops(demand);
        profiler.doneInitialWalking();
        for (const ConnectionId i : data.connectionIds()) {
            profiler.assignConnection(i);
            processConnection(i);
        }
        profiler.doneAssignment();

        profiler.doneAssignmentForDestination(destinationVertex);
    }

    inline void addGroupsToConnections() noexcept {
        assignmentData.addGroupsToConnections();
    }

    inline const AssignmentData& getAssignmentData() const noexcept {
        return assignmentData;
    }

    inline Profiler& getProfiler() noexcept {
        return profiler;
    }

private:
    inline void walkToInitialStops(const std::vector<AccumulatedVertexDemand::Entry>& demand) noexcept {
        for (const AccumulatedVertexDemand::Entry& demandEntry : demand) {
            AssertMsg(demandEntry.originVertex != demandEntry.destinationVertex, "Origin and destination vertex of demand are identical (" << demandEntry.originVertex << ")!");
            AssertMsg(settings.allowDepartureStops || !data.isStop(demandEntry.originVertex), "Demand is originating from a stop (" << demandEntry.originVertex << ")!");
            AssertMsg(data.isStop(demandEntry.originVertex) || data.transferGraph.outDegree(demandEntry.originVertex) > 0, "Origin vertex " << demandEntry.originVertex << " of demand is isolated!");
            int minTravelTime = INFTY;
            SequentialChoiceSet choiceSet = collectInitialWalkingChoices(demandEntry, minTravelTime);
            const GroupId originalGroup = assignmentData.createNewGroup(demandEntry, settings.passengerMultiplier);
            if (choiceSet.empty()) {
                assignmentData.unassignedGroups.emplace_back(originalGroup);
            } else if (choiceSet.size() == 1) {
                spawnInitialWalkingGroup(originalGroup, choiceSet.options[0]);
            } else {
                profiler.distributePassengersPATs(choiceSet.pats, choiceSet.options);
                std::vector<int> distribution;
                distribution = decisionModel.distribution(choiceSet.pats, minTravelTime);
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
                    spawnInitialWalkingGroup(group, choiceSet.options[i]);
                }
                AssertMsg(originalGroupIndex < choiceSet.size(), "No groups have been assigned!");
                AssertMsg(assignmentData.groups[originalGroup].groupSize == groupSizes[originalGroupIndex], "Original group has wrong size (size should be: " << groupSizes[originalGroupIndex] << ", size is: " << assignmentData.groups[originalGroup].groupSize << ")!");
            }
        }
    }

    inline SequentialChoiceSet collectInitialWalkingChoices(const AccumulatedVertexDemand::Entry& demandEntry, int& minTravelTime) noexcept {
        SequentialChoiceSet choiceSet;
        for (const Edge edge : data.transferGraph.edgesFrom(demandEntry.originVertex)) {
            const Vertex initialStop = data.transferGraph.get(ToVertex, edge);
            if (initialStop == demandEntry.destinationVertex) {
                const int directWalkingTime = data.transferGraph.get(TravelTime, edge);
                minTravelTime = std::min(minTravelTime, directWalkingTime);
                choiceSet.addChoice(noConnection, (settings.walkingCosts + 1) * directWalkingTime);
                continue;
            } else {
                if (!data.isStop(initialStop)) continue;
            }
            evaluateInitialStop(demandEntry, initialStop, data.transferGraph.get(TravelTime, edge), choiceSet, minTravelTime);
        }
        if (data.isStop(demandEntry.originVertex)) {
            evaluateInitialStop(demandEntry, demandEntry.originVertex, 0, choiceSet, minTravelTime);
        }
        return choiceSet;
    }

    inline void evaluateInitialStop(const AccumulatedVertexDemand::Entry& demandEntry, const Vertex stop, const int transferTime, SequentialChoiceSet& choiceSet, int& minTravelTime) noexcept {
        int earliestDepartureTime = demandEntry.earliestDepartureTime - settings.maxAdaptationTime + transferTime;
        int latestDepartureTime = demandEntry.latestDepartureTime + settings.maxAdaptationTime + transferTime;
        connectionReader.scanConnectionsInRange(stop, earliestDepartureTime, latestDepartureTime, [&](const ConnectionId i) {
            const NodeData& node = patComputation.connectionLabel(i).entryNode;
            const PerceivedTime value = node.expectedPAT;
            if (value >= Unreachable) return;
            const CSA::Connection& connection = data.connections[i];
            const int departureTime = connection.departureTime - transferTime;
            minTravelTime = std::min(minTravelTime, node.arrivalTime - departureTime);
            const PerceivedTime ptt = value - departureTime + (transferTime * (1 + settings.walkingCosts)) + getAdaptationCost(demandEntry, departureTime);
            choiceSet.addChoice(i, ptt);
        });
    }

    inline PerceivedTime getAdaptationCost(const AccumulatedVertexDemand::Entry& demandEntry, const int departureTime) const noexcept {
        const int adaptationTime = std::max(0, std::max(demandEntry.earliestDepartureTime - departureTime, departureTime - demandEntry.latestDepartureTime));
        return 60 * settings.adaptationBeta * (std::pow(adaptationTime / 60, settings.adaptationLambda) - 1) / settings.adaptationLambda;
    }

    inline void spawnInitialWalkingGroup(const GroupId group, const ConnectionId connection) noexcept {
        if (connection == noConnection) {
            assignmentData.directWalkingGroups.emplace_back(group);
        } else {
            groupTrackingData.groupsWaitingAtConnection[connection].emplace_back(group);
        }
    }

    inline void processConnection(const ConnectionId i) noexcept {
        makeBoardingDecision(i);

        const CSA::Connection& connection = data.connections[i];
        GroupList& groupsInTrip = groupTrackingData.groupsInTrip[connection.tripId];
        if (groupsInTrip.empty()) return;

        for (const GroupId group : groupsInTrip) {
            AssertMsg(group < assignmentData.connectionsPerGroup.size(), "Group " << group << " is out of bounds (0, " << assignmentData.connectionsPerGroup.size() << ")");
            assignmentData.connectionsPerGroup[group].emplace_back(i);
        }

        GroupList groupsHoppingOff;
        const ConnectionLabel& label = patComputation.connectionLabel(i);
        const int exitTravelTime = label.exitNode.arrivalTime - connection.arrivalTime;
        moveGroups(groupsInTrip, groupsHoppingOff, label.tripPAT, label.transferNode.expectedPAT, exitTravelTime, "continue", "alight");
        if (groupsHoppingOff.empty()) return;

        walkToNextStop(i, groupsHoppingOff, exitTravelTime);
    }

    inline void makeBoardingDecision(const ConnectionId i) noexcept {
        GroupList& waitingGroups = groupTrackingData.groupsWaitingAtConnection[i];
        if (waitingGroups.empty()) return;
        const CSA::Connection& connection = data.connections[i];
        const ConnectionLabel& label = patComputation.connectionLabel(i);
        const int boardingTravelTime = label.entryNode.arrivalTime - connection.departureTime;
        GroupList& groupsInTrip = groupTrackingData.groupsInTrip[connection.tripId];
        const PerceivedTime boardingPAT = label.exitNode.expectedPAT + settings.transferCosts;
        moveGroups(waitingGroups, groupsInTrip, label.skipPAT, boardingPAT, boardingTravelTime, "skip", "board");
        if (waitingGroups.empty()) return;
        const ConnectionId nextAtStop = timeExpandedNetwork.getNextDepartingConnection(i);
        if (nextAtStop == noConnection) return;
        groupTrackingData.groupsWaitingAtConnection[nextAtStop] += std::move(waitingGroups);
    }

    template<typename TO_PASSENGER_LIST>
    inline void moveGroups(GroupList& from, TO_PASSENGER_LIST& to, const PerceivedTime fromPAT, const PerceivedTime toPAT, const int travelTime, const std::string& fromName, const std::string& toName) noexcept {
        profiler.moveGroups(fromName, toName);
        profiler.moveGroupsPATs(fromPAT, toPAT);
        const std::array<int, 3> values = decisionModel.distribution(fromPAT, toPAT, travelTime);
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

    inline void walkToNextStop(const ConnectionId from, GroupList& groupList, const int travelTime) noexcept {
        SequentialChoiceSet choiceSet = collectIntermediateWalkingChoices(from);
        AssertMsg(!choiceSet.empty(), "" << groupList.size() << " groups arrived at stop " << from << " but have nowhere to go!");
        if (choiceSet.size() == 1) {
            if (choiceSet.options[0] == noConnection) return;
            groupTrackingData.groupsWaitingAtConnection[choiceSet.options[0]] += std::move(groupList);
        } else {
            profiler.distributePassengersPATs(choiceSet.pats, choiceSet.options);
            const std::vector<int> distribution = decisionModel.distribution(choiceSet.pats, travelTime);
            profiler.distributePassengersProbabilities(distribution);
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
                    if (choiceSet.options[j] == noConnection) continue;
                    groupTrackingData.groupsWaitingAtConnection[choiceSet.options[j]].emplace_back(group);
                }
                AssertMsg(movedOriginalGroup, "Group has not moved to the next stop (Group: " << assignmentData.groups[groupList[i]] << ")");
            }
        }
    }

    inline SequentialChoiceSet collectIntermediateWalkingChoices(const ConnectionId from) noexcept {
        SequentialChoiceSet choiceSet;
        const CSA::Connection& connection = data.connections[from];
        const PerceivedTime targetPAT = patComputation.targetPAT(connection);
        if (targetPAT < Unreachable) {
            choiceSet.addChoice(noConnection, targetPAT);
        }
        for (const LinkId id : timeExpandedNetwork.getOutgoingLinkIDs(from)) {
            const CSA::Link& link = timeExpandedNetwork.getLink(id);
            const PerceivedTime expectedPAT = patComputation.connectionLabel(link.headConnection).entryNode.expectedPAT;
            if (expectedPAT >= INFTY) continue;
            const PerceivedTime pat = expectedPAT + link.travelTime * settings.walkingCosts + link.waitingTime * settings.waitingCosts;
            choiceSet.addChoice(link.headConnection, pat);
        }
        return choiceSet;
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const CSA::TimeExpandedNetwork& timeExpandedNetwork;
    const Settings& settings;
    const DecisionModel& decisionModel;

    PATComputationType patComputation;
    ConnectionReader connectionReader;

    SequentialGroupTrackingData groupTrackingData;
    AssignmentData assignmentData;

    Random random;
    Profiler profiler;
};

}
