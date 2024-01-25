#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../../../DataStructures/Assignment/AssignmentData.h"
#include "../../../DataStructures/Assignment/ConnectionStatistics.h"
#include "../../../DataStructures/Assignment/GroupAssignmentStatistic.h"
#include "../../../DataStructures/Assignment/GroupData.h"
#include "../../../DataStructures/Assignment/Settings.h"
#include "../../../DataStructures/CSA/Data.h"
#include "../../../DataStructures/Demand/AccumulatedVertexDemand.h"
#include "../../../DataStructures/Demand/IdVertexDemand.h"
#include "../../../DataStructures/Demand/SplitDemand.h"
#include "../../../DataStructures/Demand/PassengerData.h"
#include "../../../DataStructures/Graph/Graph.h"

#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/MultiThreading.h"
#include "../../../Helpers/Vector/Vector.h"

#include "../../CH/CH.h"
#include "../../CH/Query/BucketQuery.h"

#include "../CycleRemoval.h"
#include "../Profiler.h"

#include "AssignmentWorker.h"

namespace Assignment::ULTRA {

template<typename DECISION_MODEL, typename PROFILER>
class Assignment {

public:
    using DecisionModel = DECISION_MODEL;
    using Profiler = PROFILER;
    using Type = Assignment<DecisionModel, Profiler>;
    using WorkerType = AssignmentWorker<DecisionModel, Profiler>;
    using DistanceLabel = typename WorkerType::DistanceLabel;

public:
    Assignment(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const CH::CH& ch, const Settings& settings) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        decisionModel(settings),
        bucketQuery(ch, FORWARD, data.numberOfStops()),
        zoneToStop(reverseGraph.numVertices()),
        assignmentData(data.numberOfConnections()),
        removedCycleConnections(0),
        removedCycles(0) {
        initializeZoneToStation();
        profiler.initialize(data);
    }

    inline void run(const AccumulatedVertexDemand& demand) noexcept {
        warning("Sequential");
        profiler.start();
        clear();
        srand(settings.randomSeed);
        SplitDemand<AccumulatedVertexDemand::Entry> demandByDestination(Construct::SplitByDestination, data, reverseGraph, demand.entries, settings.allowDepartureStops);

        WorkerType worker(data, reverseGraph, settings, decisionModel, zoneToStop, bucketQuery);

        for (size_t i = 0; i < demandByDestination.size(); i++) {
            const Vertex destinationVertex = demandByDestination.vertexAtIndex(i);
            worker.run(destinationVertex, demandByDestination[destinationVertex]);
        }

        worker.runCycleRemoval();

        finalize(worker);
        profiler.done();
        profiler.setPathsPerPassenger(assignmentData.groups.size() / (double)(demand.numberOfPassengers));
    }

    inline void run(const AccumulatedVertexDemand& demand, const int numberOfThreads, const int pinMultiplier = 1) noexcept {
        warning("Parallel");
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

            WorkerType worker(data, reverseGraph, settings, decisionModel, zoneToStop, bucketQuery);

            #pragma omp for schedule(guided,1)
            for (size_t i = 0; i < demandByDestination.size(); i++) {
                const Vertex destinationVertex = demandByDestination.vertexAtIndex(i);
                worker.run(destinationVertex, demandByDestination[destinationVertex]);
            }

            worker.runCycleRemoval();

            size_t groupOffset = 0;
            const AssignmentData& workerData = worker.getAssignmentData();
            #pragma omp critical
            {
                groupOffset = assignmentData.groups.size();
                for (const GroupData& group : workerData.groups) {
                    Assert(group.groupId + groupOffset == assignmentData.groups.size(), "Current group id is " << (group.groupId + groupOffset) << ", but should be " << assignmentData.groups.size() << "!");
                    assignmentData.groups.emplace_back(assignmentData.groups.size(), group.demandIndex, group.groupSize);
                }
                for (const GroupId group : workerData.unassignedGroups) {
                    assignmentData.unassignedGroups.emplace_back(group + groupOffset);
                }
                for (const GroupId group : workerData.directWalkingGroups) {
                    assignmentData.directWalkingGroups.emplace_back(group + groupOffset);
                }
                removedCycleConnections += worker.getRemovedCycleConnections();
                removedCycles += worker.getRemovedCycles();
                profiler += worker.getProfiler();
            }

            const size_t stepSize = (data.numberOfConnections() / numberOfThreads) + 1;
            for (size_t step = 0; step < (size_t)numberOfThreads; step++) {
                const size_t begin = ((step + threadId) % numberOfThreads) * stepSize;
                const size_t end = std::min(data.numberOfConnections(), begin + stepSize);
                for (size_t i = begin; i < end; i++) {
                    for (const GroupId group : workerData.groupsPerConnection[i]) {
                        assignmentData.groupsPerConnection[i].emplace_back(group + groupOffset);
                    }
                }
                #pragma omp barrier
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
        result += Vector::byteSize(zoneToStop);
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

    inline void writeAssignedJourneys(const std::string& fileName) const noexcept {
        IO::OFStream file(fileName);
        file << "journey,numberOfPassengers,demandIndex,journeyDescription\n";
        for (const size_t group : indices(assignmentData.groups)) {
            file << data.journeyToShortText(assignmentData.connectionsPerGroup[group]) << ",";
            file << (assignmentData.groups[group].groupSize / settings.passengerMultiplier) << ",";
            file << assignmentData.groups[group].demandIndex << ",";
            file << data.journeyToText(assignmentData.connectionsPerGroup[group]) << "\n";
        }
    }

    inline void writeAggregateStatistics(const std::string& fileName, const std::string& prefix) const noexcept {
        ConnectionStatistics statistics(data, settings, assignmentData, getPassengerCountsPerConnection());
        statistics.writeAggregateText(fileName, prefix);
    }

    inline void printStatistics(const AccumulatedVertexDemand& demand, const std::string& fileName) const noexcept {
        const std::string textFileName = fileName + ".statistics.txt";
        const std::string binaryFileName = fileName + ".statistics.binary";
        GroupAssignmentStatistic stats(data, bucketQuery, demand, assignmentData, settings.passengerMultiplier);
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
                globalPassengerListByDemandIndex[demandEntry.demandIndex].shrink_to_fit();
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
            globalPassengerListByGroupId[i].shrink_to_fit();
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
        std::vector<GlobalPassengerList>().swap(globalPassengerListByGroupId);
        std::vector<GlobalPassengerList>().swap(globalPassengerListByDemandIndex);
        for (GlobalPassengerList& list : passengersInConnection) {
            list.shrink_to_fit();
        }
        passengersInConnection.shrink_to_fit();
        unassignedPassengers.shrink_to_fit();
        walkingPassengers.shrink_to_fit();
        return PassengerData::FromApportionment(data, idVertexDemand, passengersInConnection, unassignedPassengers, walkingPassengers, ((settings.departureTimeChoice == DecisionModelWithAdaption) || (settings.departureTimeChoice == Rooftop)), false);
    }

private:
    inline void initializeZoneToStation() noexcept {
        for (const Vertex zone : reverseGraph.vertices()) {
            if (data.isStop(zone)) continue;
            if (data.transferGraph.outDegree(zone) == 0) continue;
            bucketQuery.clear();
            bucketQuery.addSource(zone);
            bucketQuery.run();
            for (const Vertex vertex : bucketQuery.getForwardPOIs()) {
                if (!data.isStop(vertex)) continue;
                zoneToStop[zone].emplace_back(StopId(vertex), bucketQuery.getForwardDistance(vertex));
            }
        }
    }

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

    //Data
    CH::BucketQuery<> bucketQuery;
    std::vector<std::vector<DistanceLabel>> zoneToStop;

    //Output
    AssignmentData assignmentData;
    u_int64_t removedCycleConnections;
    u_int64_t removedCycles;

    Profiler profiler;

};

}
