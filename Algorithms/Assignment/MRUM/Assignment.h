#pragma once

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <string>

#include "../Profiler.h"

#include "../../CH/CH.h"
#include "../../CH/Query/BucketQuery.h"

#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/Vector/Permutation.h"
#include "../../../Helpers/MultiThreading.h"
#include "../../../Helpers/ConfigFile.h"
#include "../../../Helpers/Ranges/Range.h"
#include "../../../Helpers/BinarySearch.h"

#include "../../../DataStructures/Assignment/GroupAssignmentStatistic.h"
#include "../../../DataStructures/Assignment/GroupData.h"
#include "../../../DataStructures/Assignment/Settings.h"
#include "../../../DataStructures/CSA/Data.h"
#include "../../../DataStructures/Graph/Graph.h"
#include "../../../DataStructures/Geometry/Point.h"
#include "../../../DataStructures/Container/Heap.h"
#include "../../../DataStructures/Container/Set.h"
#include "../../../DataStructures/Demand/IdVertexDemand.h"
#include "../../../DataStructures/Demand/SplitDemand.h"
#include "../../../DataStructures/Demand/PassengerData.h"
#include "SimpleSamplePATComputation.h"

namespace Assignment::MRUM {

template<uint32_t NUMBER_OF_SAMPLES, typename PROFILER>
class Assignment {

public:
    inline static constexpr uint32_t NumberOfSamples = NUMBER_OF_SAMPLES;
    using Profiler = PROFILER;
    using Type = Assignment<NumberOfSamples, Profiler>;

    using ComputePATs = SimpleSamplePATComputation<NumberOfSamples, NoPATProfiler>;
    using ConnectionLabel = typename ComputePATs::ConnectionLabel;
    using ParentSample = typename ComputePATs::ParentSample;
    using Profile = typename ComputePATs::VertexLabel;

    struct PathLabel {
        PathLabel(const CSA::Connection& c, const std::vector<int>& stationByStop) :
            time(c.departureTime),
            trip(c.tripId),
            stop(c.departureStopId),
            station(stationByStop[stop]){
        }
        inline void update(const CSA::Connection& c, int arrivalStation) noexcept {
            time = c.arrivalTime;
            trip = c.tripId;
            stop = c.arrivalStopId;
            station = arrivalStation;
        }
        int time;
        TripId trip;
        StopId stop;
        int station;
    };

    struct AlgorithmData {
        AlgorithmData(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings, const size_t numberOfIds) :
            patComputation(data, reverseGraph, settings),
            connectionsByPassenger(numberOfIds),
            stops(data.numberOfStops(), size_t(-1)),
            passengersInConnection(data.numberOfConnections()),
            localPassengersInConnection(data.numberOfConnections()),
            removedCycleConnections(0),
            removedCycles(0) {
            patComputation.sampleRandomUtilities();
        }
        ComputePATs patComputation;
        std::vector<std::vector<int>> connectionsByPassenger;
        std::vector<size_t> stops;
        std::vector<GlobalPassengerList> passengersInConnection;
        std::vector<DestinationSpecificPassengerList> localPassengersInConnection;
        GlobalPassengerList strandedPassengers;
        GlobalPassengerList walkingPassengers;
        u_int64_t removedCycleConnections;
        u_int64_t removedCycles;
    };

public:
    Assignment(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        passengersPerConnection(data.numberOfConnections()),
        removedCycleConnections(0),
        removedCycles(0) {
        stationByStop.resize(data.numberOfStops());
        for (const StopId stop : data.stops()) {
            stationByStop[stop] = stop;
            for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
                stationByStop[stop] = std::min<int>(stationByStop[stop], data.transferGraph.get(ToVertex, edge));
            }
        }
        profiler.initialize(data);
    }

    inline void run(const AccumulatedVertexDemand& demand) noexcept {
        run(IdVertexDemand::FromAccumulatedVertexDemand(demand, NumberOfSamples, settings.demandIntervalSplitTime));
    }

    inline void run(const IdVertexDemand& demand) noexcept {
        warning("Sequential");
        profiler.start();
        clear();
        SplitDemand<IdVertexDemand::Entry> demandByDestination(Construct::SplitByDestination, data, reverseGraph, demand.entries, settings.allowDepartureStops);

        AlgorithmData algorithmData(data, reverseGraph, settings, demand.numIds);

        for (size_t i = 0; i < demandByDestination.size(); i++) {
            const Vertex destinationVertex = demandByDestination.vertexAtIndex(i);
            Assert(!demandByDestination[destinationVertex].empty(), "Demand for destination vertex " << destinationVertex << " is empty!");
            Assert(data.isStop(destinationVertex) || reverseGraph.outDegree(destinationVertex) > 0, "Destination vertex " << destinationVertex << " is isolated!");
            profiler.startAssignmentForDestination(destinationVertex);

            sort(demandByDestination[destinationVertex], [](const IdVertexDemand::Entry& a, const IdVertexDemand::Entry& b){return a.departureTime < b.departureTime;});

            profiler.startPATComputation();
            algorithmData.patComputation.run(destinationVertex, demandByDestination[destinationVertex].front().departureTime);
            profiler.donePATComputation();

            profiler.startInitialWalking();
            for (const Vertex vertex : reverseGraph.vertices()) {
                algorithmData.patComputation.getProfile(vertex).setIndexFront();
            }
            walkToInitialStops(algorithmData, demandByDestination[destinationVertex], destinationVertex);
            profiler.doneInitialWalking();
            profiler.startAssignment();
            for (const ConnectionId i : data.connectionIds()) {
                const ConnectionLabel& label = algorithmData.patComputation.getConnectionLabel(i);
                for (const DestinationSpecificPassengerId id : algorithmData.localPassengersInConnection[i]) {
                    algorithmData.connectionsByPassenger[id].emplace_back(i);
                    const uint32_t sample = id % NumberOfSamples;
                    if (data.isConnection(label.parents[sample])) {
                        algorithmData.localPassengersInConnection[label.parents[sample]].emplace_back(id);
                    }
                }
                algorithmData.localPassengersInConnection[i].clear();
            }
            profiler.doneAssignment();

            profiler.startCycleElimination();
            if (settings.cycleMode == KeepCycles) keepCycles(algorithmData, destinationVertex);
            if (settings.cycleMode == RemoveStopCycles) removeStopCycles(algorithmData, destinationVertex);
            if (settings.cycleMode == RemoveStationCycles) removeStationCycles(algorithmData, destinationVertex);
            profiler.doneCycleElimination();

            profiler.doneAssignmentForDestination(destinationVertex);
        }

        removedCycleConnections += algorithmData.removedCycleConnections;
        unassignedPassengers += algorithmData.strandedPassengers;
        walkingPassengers += algorithmData.walkingPassengers;
        removedCycles += algorithmData.removedCycles;
        for (const ConnectionId i : data.connectionIds()) {
            passengersPerConnection[i] += algorithmData.passengersInConnection[i];
        }
        profiler.done();

    }

    inline void run(const IdVertexDemand& demand, const int numberOfThreads, const int pinMultiplier = 1) noexcept {
        warning("Parallel");
        profiler.start();
        clear();
        SplitDemand<IdVertexDemand::Entry> demandByDestination(Construct::SplitByDestination, data, reverseGraph, demand.entries, settings.allowDepartureStops);

        const int numCores = numberOfCores();
        omp_set_num_threads(numberOfThreads);
        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            pinThreadToCoreId((threadId * pinMultiplier) % numCores);
            Assert(omp_get_num_threads() == numberOfThreads, "Number of threads is " << omp_get_num_threads() << ", but should be " << numberOfThreads << "!");

            AlgorithmData algorithmData(data, reverseGraph, settings, demand.numIds);

            #pragma omp for schedule(guided,1)
            for (size_t i = 0; i < demandByDestination.size(); i++) {
                const Vertex destinationVertex = demandByDestination.vertexAtIndex(i);
                if (demandByDestination[destinationVertex].empty()) continue;
                if (!data.isStop(destinationVertex) && reverseGraph.outDegree(destinationVertex) == 0) continue;

                sort(demandByDestination[destinationVertex], [](const IdVertexDemand::Entry& a, const IdVertexDemand::Entry& b){return a.departureTime < b.departureTime;});

                algorithmData.patComputation.run(destinationVertex, demandByDestination[destinationVertex].front().departureTime);

                for (const Vertex vertex : reverseGraph.vertices()) {
                    algorithmData.patComputation.getProfile(vertex).setIndexFront();
                }
                walkToInitialStops(algorithmData, demandByDestination[destinationVertex], destinationVertex);
                for (const ConnectionId i : data.connectionIds()) {
                    const ConnectionLabel& label = algorithmData.patComputation.getConnectionLabel(i);
                    for (const DestinationSpecificPassengerId id : algorithmData.localPassengersInConnection[i]) {
                        algorithmData.connectionsByPassenger[id].emplace_back(i);
                        const uint32_t sample = id % NumberOfSamples;
                        if (data.isConnection(label.parents[sample])) {
                            algorithmData.localPassengersInConnection[label.parents[sample]].emplace_back(id);
                        }
                    }
                    algorithmData.localPassengersInConnection[i].clear();
                }

                if (settings.cycleMode == KeepCycles) keepCycles(algorithmData, destinationVertex);
                if (settings.cycleMode == RemoveStopCycles) removeStopCycles(algorithmData, destinationVertex);
                if (settings.cycleMode == RemoveStationCycles) removeStationCycles(algorithmData, destinationVertex);
            }

            #pragma omp critical
            {
                removedCycleConnections += algorithmData.removedCycleConnections;
                unassignedPassengers += algorithmData.strandedPassengers;
                walkingPassengers += algorithmData.walkingPassengers;
                removedCycles += algorithmData.removedCycles;
            }

            const size_t stepSize = (data.numberOfConnections() / numberOfThreads) + 1;
            for (size_t step = 0; step < (size_t)numberOfThreads; step++) {
                const size_t begin = ((step + threadId) % numberOfThreads) * stepSize;
                const size_t end = std::min(data.numberOfConnections(), begin + stepSize);
                for (size_t i = begin; i < end; i++) {
                    passengersPerConnection[i] += algorithmData.passengersInConnection[i];
                }
                #pragma omp barrier
            }
        }
        profiler.done();
    }

    inline const std::vector<GlobalPassengerList>& getPassengersPerConnection() const noexcept {
        return passengersPerConnection;
    }

    inline const GlobalPassengerList& getUnassignedPassengers() const noexcept {
        return unassignedPassengers;
    }

    inline const GlobalPassengerList& getWalkingPassengers() const noexcept {
        return walkingPassengers;
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
        long long result = Vector::byteSize(stationByStop);
        result += Vector::byteSize(passengersPerConnection);
        result += Vector::byteSize(unassignedPassengers);
        result += Vector::byteSize(walkingPassengers);
        result += sizeof(Type);
        return result;
    }

    inline void writeConnectionsWithLoad(const std::string& fileName, const int = 10) const noexcept {
        IO::OFStream file(fileName);
        file << CSA::Connection::CSV_HEADER << ",load\n";
        for (size_t connectionId = 0; connectionId < data.numberOfConnections(); connectionId++) {
            data.connections[connectionId].toCSV(file) << "," << (passengersPerConnection[connectionId].size() / static_cast<double>(NumberOfSamples)) << "\n";
        }
    }

private:
    inline void clear() noexcept {
        Vector::fill(passengersPerConnection);
        unassignedPassengers.clear();
        walkingPassengers.clear();
        removedCycleConnections = 0;
        removedCycles = 0;
    }

    inline void walkToInitialStops(AlgorithmData& algorithmData, const std::vector<IdVertexDemand::Entry>& demand, const int destinationVertex) const noexcept {
        switch(settings.departureTimeChoice) {
            case DecisionModelWithoutAdaption : {
                walkToInitialStop<DecisionModelWithoutAdaption>(algorithmData, demand, destinationVertex);
                break;
            }
            case DecisionModelWithAdaption : {
                walkToInitialStop<DecisionModelWithAdaption>(algorithmData, demand, destinationVertex);
                break;
            }
        }
    }

    template<int DEPARTURE_TIME_CHOICE>
    inline void walkToInitialStop(AlgorithmData& algorithmData, const std::vector<IdVertexDemand::Entry>& demand, const int destinationVertex) const noexcept {
        for (const IdVertexDemand::Entry& passengers : demand) {
            Assert(passengers.originVertex != passengers.destinationVertex, "Origin and destination vertex of p are identical (" << passengers.originVertex << ")!");
            Assert(settings.allowDepartureStops || !data.isStop(passengers.originVertex), "Demand is originating from a stop (" << passengers.originVertex << ")!");
            Assert(data.isStop(passengers.originVertex) || data.transferGraph.outDegree(passengers.originVertex) > 0, "Origin vertex " << passengers.originVertex << " of demand is isolated!");
            ParentSample parents;
            if (DEPARTURE_TIME_CHOICE == DecisionModelWithAdaption) {
                parents = algorithmData.patComputation.getProfile(passengers.originVertex).evaluate(passengers.departureTime, passengers.departureTime, settings.maxAdaptationTime, settings.adaptationCost);
            } else {
                parents = algorithmData.patComputation.getProfile(passengers.originVertex).evaluate(passengers.departureTime, passengers.departureTime);
            }
            for (const DestinationSpecificPassengerId id : passengers.ids) {
                const uint32_t sample = id % NumberOfSamples;
                if (data.isConnection(parents[sample])) {
                    algorithmData.localPassengersInConnection[parents[sample]].emplace_back(id);
                } else {
                    algorithmData.strandedPassengers.emplace_back(getGlobalPassengerId(destinationVertex, id));
                }
            }
        }
    }

    inline void keepCycles(AlgorithmData& algorithmData, const int destinationVertex) const noexcept {
        for (const size_t passenger : indices(algorithmData.connectionsByPassenger)) {
            std::vector<int>& connections = algorithmData.connectionsByPassenger[passenger];
            if (connections.empty()) continue;
            for (const int connection : connections) {
                algorithmData.passengersInConnection[connection].emplace_back(getGlobalPassengerId(destinationVertex, passenger));
            }
            connections.clear();
        }
    }

    inline void removeStopCycles(AlgorithmData& algorithmData, const int destinationVertex) const noexcept {
        for (const size_t passenger : indices(algorithmData.connectionsByPassenger)) {
            std::vector<int>& connections = algorithmData.connectionsByPassenger[passenger];
            if (connections.empty()) continue;
            for (int i = connections.size() - 1; i >= 0; i--) {
                const CSA::Connection& connection = data.connections[connections[i]];
                algorithmData.stops[connection.departureStopId] = i;
                algorithmData.stops[connection.arrivalStopId] = i + 1;
            }
            size_t usedConnectionsCount = 0;
            for (size_t i = connections.size() - 1; i < connections.size(); i--) {
                Assert(algorithmData.stops[data.connections[connections[i]].arrivalStopId] - 1 <= i, "Increasing path index at arrival stop from " << i << " to " << (algorithmData.stops[data.connections[connections[i]].arrivalStopId] - 1) << "!");
                i = algorithmData.stops[data.connections[connections[i]].arrivalStopId] - 1;
                //if (i < 0) break;
                algorithmData.passengersInConnection[connections[i]].emplace_back(getGlobalPassengerId(destinationVertex, passenger));
                usedConnectionsCount++;
                Assert(algorithmData.stops[data.connections[connections[i]].departureStopId] <= i, "Increasing path index at departure stop from " << i << " to " << (algorithmData.stops[data.connections[connections[i]].departureStopId]) << "!");
                i = algorithmData.stops[data.connections[connections[i]].departureStopId];
            }
            if (usedConnectionsCount == 0) {
                algorithmData.walkingPassengers.emplace_back(getGlobalPassengerId(destinationVertex, passenger));
            }
            if (connections.size() != usedConnectionsCount) {
                algorithmData.removedCycleConnections += (connections.size() - usedConnectionsCount);
                algorithmData.removedCycles++;
            }
            connections.clear();
        }
    }

    inline void removeStationCycles(AlgorithmData& algorithmData, const int destinationVertex) const noexcept {
        std::vector<int> path;
        for (const size_t passenger : indices(algorithmData.connectionsByPassenger)) {
            Assert(path.empty(), "Path contains stations from a previous iteration!");
            std::vector<int>& connections = algorithmData.connectionsByPassenger[passenger];
            if (connections.empty()) continue;
            PathLabel label(data.connections[connections.front()], stationByStop);
            path.emplace_back(label.station);
            for (const size_t i : indices(connections)) {
                algorithmData.stops[path.back()] = i; // Integer array of |stations| size, Contains last appearance index of every station in the journey-path
                path.emplace_back(stationByStop[data.connections[connections[i]].arrivalStopId]); // journey represent as sequence of stations
            }
            size_t i = 0; // real (cycle free) starting index of the journey-path
            if (algorithmData.stops[label.station] > i) { // first station of the journey-path is part of a cycle
                size_t j = algorithmData.stops[label.station]; // potential real start index of the (cycle free) journey
                while (j > i) {
                    if (path[j] == path[i]) { // check if skipping the cycle yelds a valid journey
                        const CSA::Connection& nextConnection = data.connections[connections[j]];
                        if ((nextConnection.tripId != label.trip) && (data.isCombinable<false>(label.stop, label.time, nextConnection))) break;
                    }
                    j--;
                }
                i = j; // increase journey-path start index to the end of leading cycles
            }
            size_t usedConnectionsCount = 0;
            while (i < connections.size()) {
                const CSA::Connection& connection = data.connections[connections[i]];
                if ((label.station == path.back()) && (label.trip != connection.tripId)) break; // if you can reach the destination by walking => do not board a new trip
                algorithmData.passengersInConnection[connections[i]].emplace_back(getGlobalPassengerId(destinationVertex, passenger));
                usedConnectionsCount++;
                i++;
                if (i >= connections.size()) break;
                label.update(connection, path[i]);
                if (algorithmData.stops[label.station] > i) {
                    size_t j = algorithmData.stops[label.station];
                    while (j > i) {
                        if (path[j] == path[i]) {
                            const CSA::Connection& nextConnection = data.connections[connections[j]];
                            if ((nextConnection.tripId != label.trip) && (data.isCombinable<true>(label.stop, label.time, nextConnection))) break;
                        }
                        j--;
                    }
                    i = j;
                }
            }
            if (usedConnectionsCount == 0) {
                algorithmData.walkingPassengers.emplace_back(getGlobalPassengerId(destinationVertex, passenger));
            }
            if (connections.size() != usedConnectionsCount) {
                algorithmData.removedCycleConnections += (connections.size() - usedConnectionsCount);
                algorithmData.removedCycles++;
            }
            connections.clear();
            path.clear();
        }
    }





private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;

    std::vector<int> stationByStop;

    std::vector<GlobalPassengerList> passengersPerConnection;
    GlobalPassengerList unassignedPassengers;
    GlobalPassengerList walkingPassengers;
    u_int64_t removedCycleConnections;
    u_int64_t removedCycles;

    Profiler profiler;

};

}
