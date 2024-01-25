#pragma once

#include <vector>

#include "../../DataStructures/Assignment/AssignmentData.h"
#include "../../DataStructures/Assignment/GroupData.h"
#include "../../DataStructures/Assignment/Settings.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Connection.h"


namespace Assignment {

class CycleRemoval {
private:
    struct PathLabel {
        PathLabel(const CSA::Connection& c, const std::vector<StopId>& stationByStop) :
            time(c.departureTime),
            trip(c.tripId),
            stop(c.departureStopId),
            station(stationByStop[stop]) {
        }

        inline void update(const CSA::Connection& c, StopId arrivalStation) noexcept {
            time = c.arrivalTime;
            trip = c.tripId;
            stop = c.arrivalStopId;
            station = arrivalStation;
        }

        int time;
        TripId trip;
        StopId stop;
        StopId station;
    };

public:
    CycleRemoval(const CSA::Data& data, const int mode, AssignmentData& assignmentData) :
        data(data),
        mode(mode),
        assignmentData(assignmentData),
        stationByStop(data.numberOfStops(), noStop),
        stopCycleIndex(data.numberOfStops(), size_t(-1)),
        removedCycleConnections(0),
        removedCycles(0) {
        if (mode == RemoveStationCycles) {
            for (const StopId stop : data.stops()) {
                stationByStop[stop] = stop;
                for (const Edge edge : data.transferGraph.edgesFrom(stop)) {
                    stationByStop[stop] = std::min<StopId>(stationByStop[stop], StopId(data.transferGraph.get(ToVertex, edge)));
                }
            }
        }
    }

    inline void run() noexcept {
        removedCycleConnections = 0;
        removedCycles = 0;

        switch (mode) {
            case KeepCycles:
                keepCycles();
                break;
            case RemoveStopCycles:
                removeStopCycles();
                break;
            case RemoveStationCycles:
                removeStationCycles();
                break;
        }
    }

    inline long long byteSize() const noexcept {
        return Vector::byteSize(stationByStop) + Vector::byteSize(stopCycleIndex) + 2*sizeof(u_int64_t);
    }

    inline u_int64_t getRemovedCycleConnections() const noexcept {
        return removedCycleConnections;
    }

    inline u_int64_t getRemovedCycles() const noexcept {
        return removedCycles;
    }

private:
    inline void keepCycles() noexcept {
        assignmentData.addGroupsToConnections();
    }

    inline void removeStopCycles() noexcept {
        for (GroupId group = 0; group < assignmentData.connectionsPerGroup.size(); group++) {
            std::vector<ConnectionId>& connections = assignmentData.connectionsPerGroup[group];
            if (connections.empty()) continue;
            const size_t size = connections.size();
            for (size_t i = connections.size() - 1; i < size; i--) {
                const CSA::Connection& connection = data.connections[connections[i]];
                stopCycleIndex[connection.departureStopId] = i;
                stopCycleIndex[connection.arrivalStopId] = i + 1;
            } // stopCycleIndex[x] holds the index of the first connection used after visiting x <=> (stopCycleIndex[x] - 1) is the index of the first connection reaching x
            std::vector<ConnectionId> usedConnections;
            for (size_t i = connections.size() - 1; i < size; i--) {
                AssertMsg(stopCycleIndex[data.connections[connections[i]].arrivalStopId] - 1 <= i, "Increasing path index at arrival stop from " << i << " to " << (stopCycleIndex[data.connections[connections[i]].arrivalStopId] - 1) << "!");
                i = stopCycleIndex[data.connections[connections[i]].arrivalStopId] - 1;
                if (i >= size) break;
                assignmentData.groupsPerConnection[connections[i]].emplace_back(group);
                usedConnections.emplace_back(connections[i]);
                AssertMsg(stopCycleIndex[data.connections[connections[i]].departureStopId] <= i, "Increasing path index at departure stop from " << i << " to " << (stopCycleIndex[data.connections[connections[i]].departureStopId]) << "!");
                i = stopCycleIndex[data.connections[connections[i]].departureStopId];
            }
            Vector::reverse(usedConnections);
            if (usedConnections.empty()) {
                assignmentData.directWalkingGroups.emplace_back(group);
            }
            if (connections.size() != usedConnections.size()) {
                removedCycleConnections += (connections.size() - usedConnections.size());
                removedCycles++;
            }
            usedConnections.swap(connections);
        }
    }

    inline void removeStationCycles() noexcept {
        std::vector<StopId> path;
        for (GroupId group = 0; group < assignmentData.connectionsPerGroup.size(); group++) {
            AssertMsg(path.empty(), "Path contains stations from a previous iteration!");
            std::vector<ConnectionId>& connections = assignmentData.connectionsPerGroup[group];
            if (connections.empty()) continue;
            PathLabel label(data.connections[connections.front()], stationByStop);
            path.emplace_back(label.station);
            for (const size_t i : indices(connections)) {
                stopCycleIndex[path.back()] = i; // Integer array of |stations| size, Contains last appearance index of every station in the journey-path
                path.emplace_back(stationByStop[data.connections[connections[i]].arrivalStopId]); // journey represented as sequence of stations
            }
            size_t i = 0; // real (cycle free) starting index of the journey-path
            if (stopCycleIndex[label.station] > i) { // first station of the journey-path is part of a cycle
                size_t j = stopCycleIndex[label.station]; // potential real start index of the (cycle free) journey
                while (j > i) {
                    if (path[j] == path[i]) { // check if skipping the cycle yields a valid journey
                        const CSA::Connection& nextConnection = data.connections[connections[j]];
                        if ((nextConnection.tripId != label.trip) && (data.isCombinable<false>(label.stop, label.time, nextConnection))) break;
                    }
                    j--;
                }
                i = j; // increase journey-path start index to the end of leading cycles
            }
            std::vector<ConnectionId> usedConnections;
            while (i < connections.size()) {
                const CSA::Connection& connection = data.connections[connections[i]];
                if ((label.station == path.back()) && (label.trip != connection.tripId)) break; // if you can reach the destination by walking => do not board a new trip
                assignmentData.groupsPerConnection[connections[i]].emplace_back(group);
                usedConnections.emplace_back(connections[i]);
                i++;
                if (i >= connections.size()) break;
                label.update(connection, path[i]);
                if (stopCycleIndex[label.station] > i) {
                    size_t j = stopCycleIndex[label.station];
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
            if (usedConnections.empty()) {
                assignmentData.directWalkingGroups.emplace_back(group);
            }
            if (connections.size() != usedConnections.size()) {
                removedCycleConnections += (connections.size() - usedConnections.size());
                removedCycles++;
            }
            usedConnections.swap(connections);
            path.clear();
        }
    }

private:
    const CSA::Data& data;
    const int mode;
    AssignmentData& assignmentData;

    std::vector<StopId> stationByStop;
    std::vector<size_t> stopCycleIndex;
    u_int64_t removedCycleConnections;
    u_int64_t removedCycles;
};

}
