#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../../DataStructures/Assignment/AssignmentData.h"
#include "../../DataStructures/Assignment/ConnectionStatistics.h"
#include "../../DataStructures/Assignment/GroupAssignmentStatistic.h"
#include "../../DataStructures/Assignment/GroupData.h"
#include "../../DataStructures/Assignment/JourneyWriter.h"
#include "../../DataStructures/Assignment/Settings.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/AccumulatedVertexDemand.h"
#include "../../DataStructures/Demand/IdVertexDemand.h"
#include "../../DataStructures/Demand/SplitDemand.h"
#include "../../DataStructures/Demand/PassengerData.h"
#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/Helpers.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/Vector/Vector.h"

#include "AssignmentWorker.h"
#include "CycleRemoval.h"
#include "Profiler.h"

namespace Assignment {

template<typename DECISION_MODEL, typename PROFILER, bool USE_TRANSFER_BUFFER_TIMES = false>
class GroupAssignment {

public:
    using DecisionModel = DECISION_MODEL;
    using Profiler = PROFILER;
    constexpr static inline bool UseTransferBufferTimes = USE_TRANSFER_BUFFER_TIMES;
    using Type = GroupAssignment<DecisionModel, Profiler, UseTransferBufferTimes>;
    using WorkerType = AssignmentWorker<DecisionModel, Profiler, UseTransferBufferTimes>;

public:
    GroupAssignment(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        decisionModel(settings),
        assignmentData(data.numberOfConnections()),
        removedCycleConnections(0),
        removedCycles(0) {
        profiler.initialize(data);
    }

    inline void run(const AccumulatedVertexDemand& demand, const int numberOfThreads = 1, const int pinMultiplier = 1) noexcept {
        profiler.start();
        clear();
        SplitDemand<AccumulatedVertexDemand::Entry> demandByDestination(Construct::SplitByDestination, data, reverseGraph, demand.entries, settings.allowDepartureStops);

        const int numCores = numberOfCores();
        omp_set_num_threads(numberOfThreads);
        #pragma omp parallel
        {
            srand(settings.randomSeed);
            int threadId = omp_get_thread_num();
            pinThreadToCoreId((threadId * pinMultiplier) % numCores);
            Assert(omp_get_num_threads() == numberOfThreads, "Number of threads is " << omp_get_num_threads() << ", but should be " << numberOfThreads << "!");

            WorkerType worker(data, reverseGraph, settings, decisionModel);

            #pragma omp for schedule(guided,1)
            for (size_t i = 0; i < demandByDestination.size(); i++) {
                const Vertex destinationVertex = demandByDestination.vertexAtIndex(i);
                worker.run(destinationVertex, demandByDestination[destinationVertex]);
            }

            worker.runCycleRemoval();

            #pragma omp critical
            {
                finalize(worker);
            }
        }
        profiler.done();
    }

    inline const AssignmentData& getAssignmentData() const noexcept {
        return assignmentData;
    }

    inline u_int64_t getRemovedCycleConnections() const noexcept {
        return removedCycleConnections;
    }

    inline u_int64_t getRemovedCycles() const noexcept {
        return removedCycles;
    }

    inline Profiler& getProfiler() noexcept {
        return profiler;
    }

    inline long long byteSize() const noexcept {
        long long result = assignmentData.byteSize();
        result += 2*sizeof(u_int64_t);
        return result;
    }

    inline double getPassengerCountForConnection(const ConnectionId connectionId) const noexcept {
        return (assignmentData.getConnectionLoad(connectionId) / static_cast<double>(settings.passengerMultiplier));
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
        file << CSA::Connection::CSV_HEADER << ",connectionId,load\n";
        for (const ConnectionId i : data.connectionIds()) {
            data.connections[i].toCSV(file) << "," << i.value() << "," << getPassengerCountForConnection(i) << "\n";
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
        Assert(statistics, "Cannot create output stream for: " << textFileName);
        Assert(statistics.is_open(), "Cannot open output stream for: " << textFileName);
        statistics << stats << std::endl;
        statistics.close();
        stats.serialize(binaryFileName);
    }

    inline PassengerData getPassengerData(const AccumulatedVertexDemand& demand) const noexcept {
        IdVertexDemand idVertexDemand = IdVertexDemand::FromAccumulatedVertexDemand(demand, settings.passengerMultiplier, 100000000);
        std::vector<GlobalPassengerList> globalPassengerListByDemandIndex;
        size_t idVertexDemandIndex = 0;
        for (const AccumulatedVertexDemand::Entry& demandEntry : demand.entries) {
            Assert(demandEntry.demandIndex + 1 >= globalPassengerListByDemandIndex.size(), "AccumulatedVertexDemand is not sorted by index, " << demandEntry.demandIndex << " comes after " << globalPassengerListByDemandIndex.size() << "!");
            globalPassengerListByDemandIndex.resize(demandEntry.demandIndex + 1);
            int passengerCount = demandEntry.numberOfPassengers * settings.passengerMultiplier;
            while (passengerCount > 0) {
                Assert(idVertexDemandIndex < idVertexDemand.entries.size(), "IdVertexDemandIndex is out of bounds (IdVertexDemandIndex: " << idVertexDemandIndex << ", Size: " << idVertexDemand.entries.size() << ")!");
                Assert(idVertexDemand.entries[idVertexDemandIndex].destinationVertex == demandEntry.destinationVertex, "DestinationVertex of AccumulatedVertexDemand does not match IdVertexDemand (" << idVertexDemand.entries[idVertexDemandIndex].destinationVertex << " != " << demandEntry.destinationVertex << ")!");
                Assert(idVertexDemand.entries[idVertexDemandIndex].originVertex == demandEntry.originVertex, "OriginVertex of AccumulatedVertexDemand does not match IdVertexDemand (" << idVertexDemand.entries[idVertexDemandIndex].originVertex << " != " << demandEntry.originVertex << ")!");
                Assert(idVertexDemand.entries[idVertexDemandIndex].departureTime == demandEntry.earliestDepartureTime, "DepartureTime of AccumulatedVertexDemand does not match IdVertexDemand (" << idVertexDemand.entries[idVertexDemandIndex].departureTime << " != " << demandEntry.earliestDepartureTime << ")!");
                for (const DestinationSpecificPassengerId destinationSpecificPassengerId : idVertexDemand.entries[idVertexDemandIndex].ids) {
                    globalPassengerListByDemandIndex[demandEntry.demandIndex].emplace_back(getGlobalPassengerId(demandEntry.destinationVertex, destinationSpecificPassengerId));
                    passengerCount--;
                }
                idVertexDemandIndex++;
            }
            Assert(passengerCount == 0, "Did not find IdVertexDemand for every passenger (demand index: " << demandEntry.demandIndex << ", idVertexIndex: " << idVertexDemandIndex << ", passengerCount: " << passengerCount << ")!");
        }
        std::vector<GlobalPassengerList> globalPassengerListByGroupId(assignmentData.groups.size());
        for (size_t i = assignmentData.groups.size() - 1; i < assignmentData.groups.size(); i--) {
            const GroupData& group = assignmentData.groups[i];
            Assert(group.groupSize <= globalPassengerListByDemandIndex[group.demandIndex].size(), "Not enough passengers for group (GroupSize: " << group.groupSize << ", Available passengers: " << globalPassengerListByDemandIndex[group.demandIndex].size() << ", demand index: " << group.demandIndex << ")!");
            for (size_t j = 0; j < group.groupSize; j++) {
                globalPassengerListByGroupId[i].emplace_back(globalPassengerListByDemandIndex[group.demandIndex].back());
                globalPassengerListByDemandIndex[group.demandIndex].pop_back();
            }
        }
        for (const GlobalPassengerList& globalPassengerList : globalPassengerListByDemandIndex) {
            Assert(globalPassengerList.empty(), "Passengers have not been assigned to group!");
        }
        std::vector<GlobalPassengerList> passengersInConnection(assignmentData.groupsPerConnection.size());
        for (size_t connection = 0; connection < assignmentData.groupsPerConnection.size(); connection++) {
            for (const GroupId& groupId : assignmentData.groupsPerConnection[connection]) {
                for (const GlobalPassengerId globalPassengerId : globalPassengerListByGroupId[groupId]) {
                    passengersInConnection[connection].emplace_back(globalPassengerId);
                }
            }
        }
        GlobalPassengerList unassignedPassengers;
        for (const GroupId& groupId : assignmentData.unassignedGroups) {
            for (const GlobalPassengerId globalPassengerId : globalPassengerListByGroupId[groupId]) {
                unassignedPassengers.emplace_back(globalPassengerId);
            }
        }
        GlobalPassengerList walkingPassengers;
        for (const GroupId& groupId : assignmentData.directWalkingGroups) {
            for (const GlobalPassengerId globalPassengerId : globalPassengerListByGroupId[groupId]) {
                walkingPassengers.emplace_back(globalPassengerId);
            }
        }
        return PassengerData::FromApportionment(data, idVertexDemand, passengersInConnection, unassignedPassengers, walkingPassengers, ((settings.departureTimeChoice == DecisionModelWithAdaption) || (settings.departureTimeChoice == Rooftop)));
    }

    inline void filterDemand(AccumulatedVertexDemand& demand, const size_t maxSize = -1) const noexcept {
        assignmentData.filterDemand(demand, maxSize);
    }

private:
    inline void clear() noexcept {
        assignmentData.clear();
        removedCycleConnections = 0;
        removedCycles = 0;
    }

    inline void finalize(WorkerType& worker) noexcept {
        assignmentData += worker.getAssignmentData();
        removedCycleConnections += worker.getRemovedCycleConnections();
        removedCycles += worker.getRemovedCycles();
        profiler += worker.getProfiler();
    }

private:
    //Input
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;
    const DecisionModel decisionModel;

    //Output
    AssignmentData assignmentData;
    u_int64_t removedCycleConnections;
    u_int64_t removedCycles;

    Profiler profiler;

};

}
