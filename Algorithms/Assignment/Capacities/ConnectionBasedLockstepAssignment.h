#pragma once

#include <cmath>
#include <random>
#include <vector>

#include "../../../DataStructures/Assignment/AssignmentData.h"
#include "../../../DataStructures/Assignment/ChoiceSet.h"
#include "../../../DataStructures/Assignment/ConnectionLoadData.h"
#include "../../../DataStructures/Assignment/GroupData.h"
#include "../../../DataStructures/Assignment/GroupTrackingData.h"
#include "../../../DataStructures/Assignment/Settings.h"
#include "../../../DataStructures/Container/Heap.h"
#include "../../../DataStructures/CSA/Data.h"
#include "../../../DataStructures/CSA/TimeExpandedNetwork.h"
#include "../../../DataStructures/Demand/AccumulatedVertexDemand.h"

#include "../../../Helpers/MultiThreading.h"
#include "../../../Helpers/Vector/Vector.h"
#include "ComputeLockstepPATs.h"
#include "StepIterator.h"

#include "../PassengerDistribution.h"
#include "../CycleRemoval.h"

namespace Assignment::Capacities {

template<typename DECISION_MODEL, typename PROFILER, bool USE_CONNECTION_PROFILER, bool USE_TRANSFER_BUFFER_TIMES = false>
class ConnectionBasedLockstepAssignment {

public:
    using DecisionModel = DECISION_MODEL;
    using Profiler = PROFILER;
    constexpr static inline bool UseConnectionProfiler = USE_CONNECTION_PROFILER;
    using ConnectionProfilerType = Meta::IF<UseConnectionProfiler, ConnectionProfiler, NoConnectionProfiler>;
    constexpr static inline bool UseTransferBufferTimes = USE_TRANSFER_BUFFER_TIMES;
    using Type = ConnectionBasedLockstepAssignment<DecisionModel, Profiler, UseConnectionProfiler, UseTransferBufferTimes>;

    using PATComputationType = ComputeLockstepPATs<NoPATProfiler, UseTransferBufferTimes>;
    using PATData = typename PATComputationType::PATData;
    using ConnectionLabel = typename PATComputationType::ConnectionLabel;

    using ConnectionChoiceSet = ParameterizedChoiceSet<ConnectionId>;
    using LinkChoiceSet = ParameterizedChoiceSet<LinkId>;

    struct DecisionPATs {
        double targetPAT;
        double failureTargetPAT;
        double hopOffPAT;
    };

public:
    ConnectionBasedLockstepAssignment(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const std::vector<double>& connectionCapacity, const Settings& settings) :
        data(data),
        timeExpandedNetwork(data, UseTransferBufferTimes),
        reverseGraph(reverseGraph),
        settings(settings),
        decisionModel(settings),
        destinationsWithDemand(data.transferGraph.numVertices()),
        destinationIndex(data.transferGraph.numVertices(), -1),
        groupTrackingData(data.numberOfConnections(), timeExpandedNetwork.numberOfLinks()),
        assignmentData(data.numberOfConnections()),
        currentBoardingData(data.numberOfConnections()),
        cycleRemoval(data, settings.cycleMode, assignmentData) {
        profiler.initialize(data);
        stepIterator = initializeStepIterator();
        for (ConnectionId connection : data.connectionIds()) {
            loadData.emplace_back(connectionCapacity[connection]);
        }
    }

    ~ConnectionBasedLockstepAssignment() {
        delete stepIterator;
    }

    inline void run(const AccumulatedVertexDemand& inputDemand, const int numberOfThreads = 1, const int pinMultiplier = 1) noexcept {
        clear();
        initializeDemand(inputDemand);

        const int numCores = numberOfCores();
        omp_set_num_threads(numberOfThreads);
        #pragma omp parallel
        {
            srand(settings.randomSeed);
        }

        do {
            std::cout << "Iteration " << stepIterator->getIteration() << std::endl;
            assignmentData.clear();
            Vector::fill(currentBoardingData);

            #pragma omp parallel
            {
                const int threadId = omp_get_thread_num();
                pinThreadToCoreId((threadId * pinMultiplier) % numCores);
                AssertMsg(omp_get_num_threads() == numberOfThreads, "Number of threads is " << omp_get_num_threads() << ", but should be " << numberOfThreads << "!");

                PATComputationType patComputation(data, reverseGraph, settings, loadData, patData);
                Profiler threadProfiler;

                #pragma omp for schedule(guided,1)
                for (size_t i = 0; i < destinationsWithDemand.size(); i++) {
                    const Vertex destinationVertex = destinationsWithDemand[i];
                    destinationIndex[destinationVertex] = i;
                    threadProfiler.startPATComputation();
                    patComputation.run(destinationVertex, i);
                    threadProfiler.donePATComputation();
                }

                #pragma omp critical
                profiler += threadProfiler;
            }

            profiler.startInitialWalking();
            initializeAssignment();
            profiler.doneInitialWalking();
            profiler.startAssignment();
            connectionProfiler.initialize(settings.connectionProfilerFrequency);
            for (const ConnectionId i : data.connectionIds()) {
                profiler.assignConnection(i);
                processConnection(i);
                connectionProfiler.doneConnection();
            }
            connectionProfiler.print();
            profiler.doneAssignment();

            profiler.startCycleElimination();
            cycleRemoval.run();
            profiler.doneCycleElimination();
            profiler.doneAssignmentForDestination(0);
        } while(updateConnectionLoad());
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

    inline long long byteSize() const noexcept {
        long long result = Vector::byteSize(demand);
        result += destinationsWithDemand.byteSize();
        result += Vector::byteSize(patData);
        result += Vector::byteSize(decisionPATs);
        result += Vector::byteSize(destinationIndex);
        result += sizeof(size_t);
        result += assignmentData.byteSize();
        result += Vector::byteSize(currentBoardingData);
        result += Vector::byteSize(loadData);
        return result;
    }

    inline double getPassengerCountForConnection(const ConnectionId connectionId) const noexcept {
        return assignmentData.getConnectionLoad(connectionId) / static_cast<double>(settings.passengerMultiplier);
    }

    inline std::vector<double> getPassengerCountsPerConnection() const noexcept {
        std::vector<double> passengerCounts(data.numberOfConnections());
        for (const ConnectionId i : data.connectionIds()) {
            passengerCounts[i] = getPassengerCountForConnection(i);
        }
        return passengerCounts;
    }

    inline void writeConnectionsWithLoad(const std::string& fileName) const noexcept {
        IO::OFStream file(fileName);
        file << CSA::Connection::CSV_HEADER << ",connectionId,load,boardingProbability\n";
        for (const ConnectionId i : data.connectionIds()) {
            data.connections[i].toCSV(file) << "," << i.value() << "," << getPassengerCountForConnection(i) << "," << loadData[i].boardingProbability() << "\n";
        }
    }

    inline void writeAssignment(const std::string& fileName) const noexcept {
        assignmentData.writeAssignment(fileName);
    }

    inline void writeGroups(const std::string& fileName) const noexcept {
        assignmentData.writeGroups(fileName);
    }

    inline void writeAssignedJourneys(const std::string& fileName, const AccumulatedVertexDemand& demand) const noexcept {
        JourneyWriter journeyWriter(data, settings, demand, assignmentData);
        journeyWriter.write(fileName);
    }

    inline void writeConnectionStatistics(const std::string& fileName, const std::string& prefix) const noexcept {
        ConnectionStatistics statistics(data, settings, assignmentData, getPassengerCountsPerConnection());
        statistics.write(fileName, prefix);
    }

    inline void printStatistics(const AccumulatedVertexDemand& demand, const std::string& fileName) const noexcept {
        const std::string textFileName = fileName + ".statistics.txt";
        const std::string binaryFileName = fileName + ".statistics.binary";
        GroupAssignmentStatistic stats(data, demand, assignmentData, settings.passengerMultiplier);
        std::cout << stats << std::endl;
        std::ofstream statistics(textFileName);
        AssertMsg(statistics, "Cannot create output stream for: " << textFileName);
        AssertMsg(statistics.is_open(), "Cannot open output stream for: " << textFileName);
        statistics << stats << std::endl;
        statistics.close();
        stats.serialize(binaryFileName);
    }

private:
    inline void clear() noexcept {
        demand.clear();
        for (ConnectionLoadData& load : loadData) {
            load.clear();
        }
        Vector::fill(destinationIndex, size_t(-1));
        assignmentData.clear();
        Vector::fill(currentBoardingData);
    }

    inline void initializeDemand(const AccumulatedVertexDemand& inputDemand) noexcept {
        originalDemand = &inputDemand.entries;
        for (const AccumulatedVertexDemand::Entry& entry : inputDemand.entries) {
            if (entry.originVertex == entry.destinationVertex) continue;
            if (!settings.allowDepartureStops && data.isStop(entry.originVertex)) continue;
            if (!data.isStop(entry.originVertex) && data.transferGraph.outDegree(entry.originVertex) == 0) continue;
            if (!data.isStop(entry.destinationVertex) && reverseGraph.outDegree(entry.destinationVertex) == 0) continue;
            demand.push_back(entry);
            destinationsWithDemand.insert(entry.destinationVertex);
        }
        sort(demand, [](const AccumulatedVertexDemand::Entry& a, const AccumulatedVertexDemand::Entry& b){return a.earliestDepartureTime < b.earliestDepartureTime;});
        patData.resize(destinationsWithDemand.size(), PATData(data.numberOfStops(), data.numberOfConnections()));
        decisionPATs.resize(patData.size());
    }

    inline void initializeAssignment() noexcept {
        groupTrackingData.validate();
        walkToInitialStops();
    }

    inline void walkToInitialStops() noexcept {
        switch(settings.departureTimeChoice) {
            case DecisionModelWithoutAdaption:
                walkToInitialStops<DecisionModelWithoutAdaption>();
                break;
            case DecisionModelWithAdaption:
                walkToInitialStops<DecisionModelWithAdaption>();
                break;
            case Uniform:
                walkToInitialStops<Uniform>();
                break;
            case Rooftop:
                walkToInitialStops<Rooftop>();
                break;
            case DecisionModelWithBoxCox:
                walkToInitialStops<DecisionModelWithBoxCox>();
                break;
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline void walkToInitialStops() noexcept {
        for (const AccumulatedVertexDemand::Entry& demandEntry : demand) {
            AssertMsg(demandEntry.originVertex != demandEntry.destinationVertex, "Origin and destination vertex of demand are identical (" << demandEntry.originVertex << ")!");
            AssertMsg(settings.allowDepartureStops || !data.isStop(demandEntry.originVertex), "Demand is originating from a stop (" << demandEntry.originVertex << ")!");
            AssertMsg(data.isStop(demandEntry.originVertex) || data.transferGraph.outDegree(demandEntry.originVertex) > 0, "Origin vertex " << demandEntry.originVertex << " of demand is isolated!");
            AssertMsg(data.isStop(demandEntry.destinationVertex) || reverseGraph.outDegree(demandEntry.destinationVertex) > 0, "Destination vertex " << demandEntry.destinationVertex << " of demand is isolated!");
            ConnectionChoiceSet choiceSet = collectInitialWalkingChoices<DEPARTURE_TIME_CHOICE>(demandEntry);
            const GroupId originalGroup = assignmentData.createNewGroup(demandEntry, settings.passengerMultiplier);
            if (choiceSet.empty()) {
                assignmentData.unassignedGroups.emplace_back(originalGroup);
            } else if (choiceSet.size() == 1) {
                groupTrackingData.groupsWaitingAtConnection[choiceSet.options[0]].emplace_back(originalGroup);
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
                    groupTrackingData.groupsWaitingAtConnection[choiceSet.options[i]].emplace_back(group);
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
    inline ConnectionChoiceSet collectInitialWalkingChoices(const AccumulatedVertexDemand::Entry& demandEntry) noexcept {
        ConnectionChoiceSet choiceSet;
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
    inline void evaluateInitialStop(const AccumulatedVertexDemand::Entry& demandEntry, const Vertex stop, const int transferTime, ConnectionChoiceSet& choiceSet) noexcept {
        int departureTime = demandEntry.earliestDepartureTime - getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        const int latestDepartureTime = demandEntry.latestDepartureTime + getMaxAdaptationTime<DEPARTURE_TIME_CHOICE>() + transferTime;
        while (departureTime <= latestDepartureTime) {
            const ProfileEntry& entry = patData[destinationIndex[demandEntry.destinationVertex]].profiles.findEntry(stop, departureTime);
            departureTime = entry.departureTime;
            if constexpr ((DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) | (DEPARTURE_TIME_CHOICE == DecisionModelWithBoxCox)) {
                if (departureTime > latestDepartureTime) return;
            }
            const int value = entry.evaluate(departureTime, settings.waitingCosts);
            if (value >= INFTY) return;
            const int pat = value - departureTime + (transferTime * (1 + settings.walkingCosts)) + getAdaptationCost<DEPARTURE_TIME_CHOICE>(demandEntry, departureTime - transferTime);
            choiceSet.addChoice(entry.connectionId, departureTime, pat);
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
        connectionProfiler.startCollectIncomingGroups();
        collectIncomingGroups(i);
        connectionProfiler.stopCollectIncomingGroups();
        connectionProfiler.startCalculateDecisionPATs();
        calculateDecisionPATs(i);
        connectionProfiler.stopCalculateDecisionPATs();

        connectionProfiler.startBoarding();
        GroupList& waitingGroups = groupTrackingData.groupsWaitingAtConnection[i];
        GroupList& groupsInTrip = groupTrackingData.groupsInTrip[i];
        const size_t inTripOffset = groupsInTrip.size();
        const size_t currentLoad = Vector::sum(groupsInTrip, [&](const GroupId group) {
            return assignmentData.groups[group].groupSize;
        });
        currentBoardingData[i].boardingCapacity = loadData[i].capacity * settings.passengerMultiplier - currentLoad;
        GroupList walkingGroups;
        currentBoardingData[i].boardingDemand = moveWaitingGroupsIntoTrip(waitingGroups, groupsInTrip, walkingGroups, currentBoardingData[i].boardingCapacity, [&](const size_t destination) {
            return patData[destination].connectionLabels[i].skipPAT;
        }, [&](const size_t destination) {
            return patData[destination].connectionLabels[i].hopOnPAT;
        });
        mergeNewGroups<false>(groupsInTrip, inTripOffset);
        connectionProfiler.stopBoarding();
        if (!walkingGroups.empty()) {
            processWalkingDecisions<true>(i, walkingGroups);
        }

        connectionProfiler.startAssignment();
        for (const GroupId group : groupsInTrip) {
            AssertMsg(group < assignmentData.connectionsPerGroup.size(), "Group " << group << " is out of bounds (0, " << assignmentData.connectionsPerGroup.size() << ")");
            assignmentData.connectionsPerGroup[group].emplace_back(i);
        }
        connectionProfiler.stopAssignment();

        connectionProfiler.startAlighting();
        GroupList groupsHoppingOff;
        moveGroups(groupsInTrip, groupsHoppingOff, [&](const size_t destination) {
            return patData[destination].connectionLabels[i].tripPAT;
        }, [&](const size_t destination) {
            return decisionPATs[destination].hopOffPAT;
        }, "continue", "alight");
        connectionProfiler.stopAlighting();
        if (groupsHoppingOff.empty()) return;

        processWalkingDecisions<false>(i, groupsHoppingOff);
    }

    template<bool FAILED_TO_BOARD>
    inline void processWalkingDecisions(const ConnectionId i, GroupList& groups) noexcept {
        connectionProfiler.startTarget();
        moveGroups(groups, groupTrackingData.groupsAtTarget, [&](const size_t destination) {
            if constexpr (FAILED_TO_BOARD) {
                return patData[destination].connectionLabels[i].failureTransferPAT;
            } else {
                return patData[destination].connectionLabels[i].transferPAT;
            }
        }, [&](const size_t destination) {
            if constexpr (FAILED_TO_BOARD) {
                return decisionPATs[destination].failureTargetPAT;
            } else {
                return decisionPATs[destination].targetPAT;
            }
        }, "walk", "target");
        connectionProfiler.stopTarget();
        if (groups.empty()) return;

        connectionProfiler.startWalking();
        walkToNextStop<FAILED_TO_BOARD>(i, groups);
        connectionProfiler.stopWalking();
    }


    inline void collectIncomingGroups(const ConnectionId i) noexcept {
        GroupList& waitingGroups = groupTrackingData.groupsWaitingAtConnection[i];
        const size_t waitingOffset = waitingGroups.size();
        for (const LinkId link : timeExpandedNetwork.getIncomingLinkIDs(i)) {
            waitingGroups += std::move(groupTrackingData.groupsWalkingOnLink[link]);
        }
        const ConnectionId previousAtStop = timeExpandedNetwork.getPreviousDepartingConnection(i);
        if (previousAtStop != noConnection) {
            waitingGroups += std::move(groupTrackingData.groupsWaitingAtConnection[previousAtStop]);
        }
        mergeNewGroups<true>(waitingGroups, waitingOffset);

        GroupList& groupsInTrip = groupTrackingData.groupsInTrip[i];
        const ConnectionId previousInTrip = timeExpandedNetwork.getPreviousConnectionInTrip(i);
        if (previousInTrip != noConnection) {
            groupTrackingData.groupsInTrip[previousInTrip].swap(groupsInTrip);
        }
    }

    template<bool SORT_NEW_GROUPS>
    inline void mergeNewGroups(GroupList& groupList, const size_t newOffset) noexcept {
        auto comp = [&](const GroupId a, const GroupId b) {
            return getDestinationOfGroup(a) < getDestinationOfGroup(b);
        };
        if constexpr (SORT_NEW_GROUPS) sort(groupList.begin() + newOffset, groupList.end(), comp);
        std::inplace_merge(groupList.begin(), groupList.begin() + newOffset, groupList.end(), comp);
    }

    inline Vertex getDestinationOfGroup(const GroupId group) const noexcept {
        return (*originalDemand)[assignmentData.groups[group].demandIndex].destinationVertex;
    }

    inline void calculateDecisionPATs(const ConnectionId i) noexcept {
        const CSA::Connection& connection = data.connections[i];
        IndexedSet<false, size_t> destinations(patData.size()); //TODO: Reuse between connections?
        for (const GroupId group : groupTrackingData.groupsWaitingAtConnection[i]) {
            destinations.insert(destinationIndex[getDestinationOfGroup(group)]);
        }
        for (const GroupId group : groupTrackingData.groupsInTrip[i]) {
            destinations.insert(destinationIndex[getDestinationOfGroup(group)]);
        }
        for (const size_t destination : destinations) {
            const ConnectionLabel& label = patData[destination].connectionLabels[i];
            const double hopOffLoadCost = label.loadFactor * settings.congestionExitCosts;
            decisionPATs[destination].targetPAT = patData[destination].targetPAT(connection);
            decisionPATs[destination].failureTargetPAT = patData[destination].failureTargetPAT(connection);
            decisionPATs[destination].hopOffPAT = std::min(decisionPATs[destination].targetPAT, label.transferPAT) + hopOffLoadCost;
        }
    }

    template<typename STAY_PAT, typename MOVE_PAT>
    inline size_t moveWaitingGroupsIntoTrip(GroupList& waitingGroups, GroupList& boardingGroups, GroupList& walkingGroups, const size_t remainingSpace, const STAY_PAT& stayPAT, const MOVE_PAT& movePAT) noexcept {
        if (waitingGroups.empty()) return 0;
        profiler.moveGroups("skip", "board");
        size_t destination = -1;
        std::array<int, 3> values;
        std::vector<int> waitingGroupSizes(waitingGroups.size());
        std::vector<int> boardingGroupSizes(waitingGroups.size() + 1);
        std::vector<int> walkingGroupSizes(waitingGroups.size(), 0);
        for (size_t i = 0; i < waitingGroups.size(); i++) {
            const size_t newDestination = destinationIndex[getDestinationOfGroup(waitingGroups[i])];
            if (newDestination != destination) {
                destination = newDestination;
                profiler.moveGroupsDestination(destination);
                profiler.moveGroupsPATs(stayPAT(destination), movePAT(destination));
                values = decisionModel.distribution(stayPAT(destination), movePAT(destination));
                profiler.moveGroupsProbabilities(values);
            }
            const std::array<int, 2> groupSizes = getGroupSizes(values, assignmentData.groups[waitingGroups[i]].groupSize, random);
            profiler.moveGroupsSizes(groupSizes);
            waitingGroupSizes[i] = groupSizes[0];
            boardingGroupSizes[i] = groupSizes[1];
            boardingGroupSizes.back() += groupSizes[1];
        }
        const size_t desiredLoad = boardingGroupSizes.back();
        if (desiredLoad > remainingSpace) {
            const std::vector<size_t> rejectedPassengers = getGroupSizes(boardingGroupSizes, desiredLoad - remainingSpace, random);
            profiler.distributePassengersSizes(rejectedPassengers);
            for (size_t i = 0; i < waitingGroupSizes.size(); i++) {
                walkingGroupSizes[i] += rejectedPassengers[i];
                boardingGroupSizes[i] -= rejectedPassengers[i];
            }
        }
        size_t offset = 0;
        for (size_t i = 0; i < waitingGroups.size(); i++) {
            const size_t newI = i - offset;
            if (offset > 0) {
                waitingGroups[newI] = waitingGroups[i];
            }
            if (waitingGroupSizes[i] == 0) {
                if (boardingGroupSizes[i] == 0) {
                    walkingGroups.emplace_back(waitingGroups[newI]);
                } else {
                    boardingGroups.emplace_back(waitingGroups[newI]);
                    if (walkingGroupSizes[i] != 0) {
                        walkingGroups.emplace_back(assignmentData.splitGroup(boardingGroups.back(), walkingGroupSizes[i]));
                    }
                }
                offset++;
            } else {
                if (boardingGroupSizes[i] != 0) {
                    boardingGroups.emplace_back(assignmentData.splitGroup(waitingGroups[newI], boardingGroupSizes[i]));
                }
                if (walkingGroupSizes[i] != 0) {
                    walkingGroups.emplace_back(assignmentData.splitGroup(waitingGroups[newI], walkingGroupSizes[i]));
                }
            }
        }
        if (offset > 0) {
            waitingGroups.resize(waitingGroups.size() - offset);
        }
        return desiredLoad;
    }

    template<typename TO_PASSENGER_LIST, typename STAY_PAT, typename MOVE_PAT>
    inline void moveGroups(GroupList& from, TO_PASSENGER_LIST& to, const STAY_PAT& stayPAT, const MOVE_PAT& movePAT, const std::string& fromName, const std::string& toName) noexcept {
        if (from.empty()) return;
        profiler.moveGroups(fromName, toName);
        size_t destination = -1;
        std::array<int, 3> values;
        size_t offset = 0;
        for (size_t i = 0; i < from.size(); i++) {
            const size_t newDestination = destinationIndex[getDestinationOfGroup(from[i])];
            if (newDestination != destination) {
                destination = newDestination;
                profiler.moveGroupsDestination(destination);
                profiler.moveGroupsPATs(stayPAT(destination), movePAT(destination));
                values = decisionModel.distribution(stayPAT(destination), movePAT(destination));
                profiler.moveGroupsProbabilities(values);
            }
            const std::array<int, 2> groupSizes = getGroupSizes(values, assignmentData.groups[from[i]].groupSize, random);
            profiler.moveGroupsSizes(groupSizes);
            const size_t newI = i - offset;
            if (offset > 0) {
                from[newI] = from[i];
            }
            if (groupSizes[0] == 0) {
                to.emplace_back(from[newI]);
                offset++;
            } else if (groupSizes[1] != 0) {
                to.emplace_back(assignmentData.splitGroup(from[newI], groupSizes[1]));
            }
        }
        if (offset > 0) {
            from.resize(from.size() - offset);
        }
    }

    template<bool FAILED_TO_BOARD>
    inline void walkToNextStop(const ConnectionId from, GroupList& groupList) noexcept {
        size_t destination = -1;
        LinkChoiceSet choiceSet;
        std::vector<int> distribution;
        std::vector<size_t> groupSizes;
        std::vector<SampleElement> sampleElements;
        for (size_t i = 0; i < groupList.size(); i++) {
            const size_t newDestination = destinationIndex[getDestinationOfGroup(groupList[i])];
            if (newDestination != destination) {
                destination = newDestination;
                choiceSet = collectIntermediateWalkingChoices<FAILED_TO_BOARD>(from, destination);
                AssertMsg(FAILED_TO_BOARD || !choiceSet.empty(), "" << groupList.size() << " groups arrived at connection " << from << " but have nowhere to go!");
                if (choiceSet.size() > 1) {
                    profiler.distributePassengersPATs(choiceSet.pats, choiceSet.departureTimes);
                    decisionModel.distribution(choiceSet.pats, distribution);
                    profiler.distributePassengersProbabilities(distribution);
                }
                groupSizes.resize(choiceSet.size());
                sampleElements.resize(choiceSet.size());
            }
            if constexpr (FAILED_TO_BOARD) {
                if (choiceSet.empty()) {
                    assignmentData.unassignedGroups.emplace_back(groupList[i]);
                    continue;
                }
            }
            if (choiceSet.size() == 1) {
                groupTrackingData.groupsWalkingOnLink[choiceSet.options[0]].emplace_back(groupList[i]);
            } else {
                getGroupSizes(distribution, assignmentData.groups[groupList[i]].groupSize, random, groupSizes, sampleElements);
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
                    groupTrackingData.groupsWalkingOnLink[choiceSet.options[j]].emplace_back(group);
                }
                AssertMsg(movedOriginalGroup, "Group has not moved to the next stop (Group: " << assignmentData.groups[groupList[i]] << ")");
            }
        }
    }

    template<bool FAILED_TO_BOARD>
    inline LinkChoiceSet collectIntermediateWalkingChoices(const ConnectionId from, const size_t destination) noexcept {
        LinkChoiceSet choiceSet;
        const CSA::Connection& connection = data.connections[from];
        const StopId fromStop = FAILED_TO_BOARD ? connection.departureStopId : connection.arrivalStopId;
        const int time = FAILED_TO_BOARD ? connection.departureTime : connection.arrivalTime;
        for (const LinkId id : FAILED_TO_BOARD ? timeExpandedNetwork.getOutgoingFailureLinkIDs(from) : timeExpandedNetwork.getOutgoingLinkIDs(from)) {
            const CSA::Link& link = timeExpandedNetwork.getLink(id);
            const CSA::Connection& nextConnection = data.connections[link.headConnection];
            const StopId nextStop = nextConnection.departureStopId;
            const int bufferTime = (nextStop == fromStop || UseTransferBufferTimes) ? data.minTransferTime(nextStop) : 0;
            const int bestPAT = patData[destination].connectionLabels[link.headConnection].bestPAT();
            if (bestPAT >= INFTY) continue;
            const int departureTime = time + link.travelTime + bufferTime;
            const int waitingTime = bufferTime + (nextConnection.departureTime - departureTime);
            const int pat = bestPAT + (link.travelTime * settings.walkingCosts) + (waitingTime * settings.waitingCosts);
            choiceSet.addChoice(id, departureTime, pat);
        }
        return choiceSet;
    }

    inline bool updateConnectionLoad() noexcept {
        size_t unfinishedConnections = 0;
        double maxDiff = 0.0;
        const std::array<double, 3> stepSizes = stepIterator->getStepSizes();
        for (const ConnectionId connection : data.connectionIds()) {
            const double currentLoad = getPassengerCountForConnection(connection);
            AssertMsg(currentLoad <= loadData[connection].capacity, "Connection " << connection << " has a load of " << currentLoad << ", but only holds " << loadData[connection].capacity << " passengers!");
            const double newLoad = applyIteration(loadData[connection].load, currentLoad, stepSizes);
            BoardingData& boardingData = loadData[connection].boardingData;
            BoardingData newBoardingData;
            newBoardingData.boardingDemand = applyIteration(boardingData.boardingDemand, currentBoardingData[connection].boardingDemand, stepSizes);
            newBoardingData.boardingCapacity = applyIteration(boardingData.boardingCapacity, currentBoardingData[connection].boardingCapacity, stepSizes);
            profiler.printConnectionStatistics(connection, loadData[connection].load, currentLoad, newLoad, boardingData, currentBoardingData[connection], newBoardingData);
            const double diff = fabs(newBoardingData.boardingDemand - boardingData.boardingDemand)/(loadData[connection].capacity * settings.passengerMultiplier);
            maxDiff = std::max(maxDiff, diff);
            if (diff >= settings.convergenceLimit) unfinishedConnections++;
            loadData[connection].load = newLoad;
            boardingData = newBoardingData;
        }
        if (unfinishedConnections == 0) return false;
        ++(*stepIterator);
        std::cout << "\tUnfinished connections: " << unfinishedConnections << "/" << data.numberOfConnections() << std::endl;
        std::cout << "\tMaximum relative load difference: " << maxDiff << std::endl;
        return true;
    }

    inline StepIterator* initializeStepIterator() const noexcept {
        switch (settings.iterationType) {
            case 0:
                return new StepIterator();
            case 1:
                return new WeightedStepIterator(settings.stepExponent);
            case 2:
                return new SelfRegulatingStepIterator(assignmentData, loadData, data.numberOfConnections(), settings.passengerMultiplier, settings.speedupStepSize, settings.slowdownStepSize);
            default:
                return new StepIterator();
        }
    }

private:
    //Input
    const CSA::Data& data;
    CSA::CapacityTimeExpandedNetwork timeExpandedNetwork;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;
    const DecisionModel decisionModel;

    //Demand
    const std::vector<AccumulatedVertexDemand::Entry>* originalDemand;
    std::vector<AccumulatedVertexDemand::Entry> demand;
    IndexedSet<false, Vertex> destinationsWithDemand;

    //Data
    std::vector<PATData> patData;
    std::vector<DecisionPATs> decisionPATs;
    std::vector<size_t> destinationIndex;
    ConnectionBasedGroupTrackingData groupTrackingData;
    StepIterator* stepIterator;

    //Output
    AssignmentData assignmentData;
    std::vector<BoardingData> currentBoardingData;
    std::vector<ConnectionLoadData> loadData;

    CycleRemoval cycleRemoval;

    Random random;
    Profiler profiler;
    ConnectionProfilerType connectionProfiler;
};

}
