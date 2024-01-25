#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

#include "../Demand/Passenger.h"
#include "../Demand/AccumulatedVertexDemand.h"

#include "../CSA/Data.h"
#include "../Geometry/Point.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Console/Progress.h"

#include "AssignmentData.h"
#include "GroupData.h"

#include "../Container/Set.h"
#include "../Container/Map.h"

#include "../../Algorithms/Dijkstra/Dijkstra.h"
#include "../../Algorithms/CH/Query/BucketQuery.h"

namespace Assignment {

class GroupAssignmentStatistic {

private:
    template<typename INT_TYPE = int32_t>
    struct Entry {
        Entry(const INT_TYPE value = 0) :
            numberOfTrips(value),
            numberOfConnections(value),
            timeInVehicle(value),
            walkingTime(value),
            travelTimeWithoutInitialWaiting(value),
            travelTimeWithInitialWaiting(value),
            groupSize(value) {
        }
        inline std::ostream& toCSV(std::ostream& out) const {
            out << numberOfTrips << ", " << numberOfConnections << ", ";
            out << timeInVehicle << ", " << walkingTime << ", ";
            out << travelTimeWithoutInitialWaiting << ", " << travelTimeWithInitialWaiting << ", " << groupSize << "\n";
            return out;
        }
        inline void print() const noexcept {
            std::cout << "PassengerData:\n";
            std::cout << "numberOfTrips:                   " << std::setw(12) << String::prettyInt(numberOfTrips)                     << "\n";
            std::cout << "numberOfConnections:             " << std::setw(12) << String::prettyInt(numberOfConnections)               << "\n";
            std::cout << "timeInVehicle:                   " << std::setw(12) << String::secToString(timeInVehicle)                   << "\n";
            std::cout << "walkingTime:                     " << std::setw(12) << String::secToString(walkingTime)                     << "\n";
            std::cout << "travelTimeWithoutInitialWaiting: " << std::setw(12) << String::secToString(travelTimeWithoutInitialWaiting) << "\n";
            std::cout << "travelTimeWithInitialWaiting:    " << std::setw(12) << String::secToString(travelTimeWithInitialWaiting)    << "\n";
            std::cout << "groupSize:                       " << std::setw(12) << String::prettyDouble(groupSize)                      << "\n";
        }
        template<typename FUNCTION, typename OTHER>
        inline void apply(const OTHER& other, const FUNCTION& function) noexcept {
            numberOfTrips = function(numberOfTrips, other.numberOfTrips);
            numberOfConnections = function(numberOfConnections, other.numberOfConnections);
            timeInVehicle = function(timeInVehicle, other.timeInVehicle);
            walkingTime = function(walkingTime, other.walkingTime);
            travelTimeWithoutInitialWaiting = function(travelTimeWithoutInitialWaiting, other.travelTimeWithoutInitialWaiting);
            travelTimeWithInitialWaiting = function(travelTimeWithInitialWaiting, other.travelTimeWithInitialWaiting);
            groupSize = function(groupSize, other.groupSize);
        }
        template<typename OTHER>
        inline void maximize(const OTHER& other) noexcept {
            apply(other, [](const auto& a, const auto& b){return std::max(a, b);});
        }
        template<typename OTHER>
        inline void minimize(const OTHER& other) noexcept {
            apply(other, [](const auto& a, const auto& b){return std::min(a, b);});
        }
        template<typename OTHER>
        inline void add(const OTHER& other) noexcept {
            apply(other, [&](const auto& a, const auto& b){return a + (b * other.groupSize);});
        }
        INT_TYPE numberOfTrips{0};
        INT_TYPE numberOfConnections{0};
        INT_TYPE timeInVehicle{0};
        INT_TYPE walkingTime{0};
        INT_TYPE travelTimeWithoutInitialWaiting{0};
        INT_TYPE travelTimeWithInitialWaiting{0};
        double groupSize{0};
    };

    struct Path {
        Path(const std::vector<ConnectionId>& p) {
            data.reserve(p.size());
            for (const ConnectionId s : p) {
                data.emplace_back(s);
            }
        }
        inline bool operator==(const Path& other) const noexcept {
            if (data.size() != other.data.size()) return false;
            for (size_t i = 0; i < data.size(); i++) {
                if (data[i] != other.data[i]) return false;
            }
            return true;
        }
        inline bool operator<(const Path& other) const noexcept {
            if (data.size() == other.data.size()) {
                for (size_t i = 0; i < data.size(); i++) {
                    if (data[i] < other.data[i]) return true;
                    if (data[i] > other.data[i]) return false;
                }
                return false;
            } else {
                return data.size() < other.data.size();
            }
        }
        std::vector<int32_t> data;
    };

    struct DistanceLabel {
        DistanceLabel(const int initialDistance = -1, const int finalDistance = -1) :
            initialDistance(initialDistance),
            finalDistance(finalDistance) {
        }
        int initialDistance;
        int finalDistance;
    };

    struct DijkstraDistance {
        DijkstraDistance(const TransferGraph& graph) :
            dijkstra(graph),
            r(false),
            d(0) {
        }
        inline void run(const Vertex from, const Vertex to) noexcept {
            dijkstra.run(from, to);
            r = dijkstra.reachable(to);
            if (r) {
                d = dijkstra.getDistance(to);
            }
        }
        inline bool reachable() const noexcept {
            return r;
        }
        inline int getDistance() const noexcept {
            return d;
        }
        Dijkstra<TransferGraph> dijkstra;
        bool r;
        int d;
    };

    struct BucketDistance {
        BucketDistance(const TransferGraph& graph, const size_t numberOfStops, const CH::BucketQuery<>& bucketQuery) :
            numberOfStops(numberOfStops),
            bucketQuery(bucketQuery),
            dijkstra(graph),
            r(false),
            d(0) {
        }
        inline void run(const Vertex from, const Vertex to) noexcept {
            if (from >= numberOfStops || to >= numberOfStops) {
                bucketQuery.run(from, to);
                r = bucketQuery.reachable();
                if (r) {
                    d = bucketQuery.getDistance();
                }
            } else {
                dijkstra.run(from, to);
                r = dijkstra.reachable(to);
                if (r) {
                    d = dijkstra.getDistance(to);
                }
            }
        }
        inline bool reachable() const noexcept {
            return r;
        }
        inline int getDistance() const noexcept {
            return d;
        }
        size_t numberOfStops;
        CH::BucketQuery<> bucketQuery;
        Dijkstra<TransferGraph> dijkstra;
        bool r;
        int d;
    };

public:
    GroupAssignmentStatistic() {}

    GroupAssignmentStatistic(const std::string& filename) {
        deserialize(filename);
    }

    template<typename DIST>
    GroupAssignmentStatistic(DIST dist, const CSA::Data& data, const AccumulatedVertexDemand& demands, const AssignmentData& assignmentData, const int passengerMultiplier) :
        min(INFTY),
        max(0),
        sum(0),
        totalNumberOfTrips(data.numberOfTrips()),
        totalNumberOfGroups(assignmentData.groups.size()),
        totalNumberOfPassengers(0),
        totalNumberOfConnections(data.numberOfConnections()),
        numberOfWalkingPassengers(0),
        numberOfUnassignedPassengers(0),
        numberOfEmptyConnections(0),
        numberOfEmptyTrips(0),
        numberOfUsedPaths(0),
        passengerMultiplier(passengerMultiplier) {
        std::vector<bool> tripIsUsed(totalNumberOfTrips, false);
        for (size_t i = 0; i < assignmentData.groupsPerConnection.size(); i++) {
            if (assignmentData.groupsPerConnection[i].empty()) {
                numberOfEmptyConnections++;
            } else {
                tripIsUsed[data.connections[i].tripId] = true;
            }
        }
        numberOfEmptyTrips = Vector::count(tripIsUsed, false);
        std::vector<bool> groupIsWalking(assignmentData.groups.size(), false);
        for (const GroupId group : assignmentData.directWalkingGroups) {
            groupIsWalking[group] = true;
        }
        std::vector<bool> groupIsUnassigned(assignmentData.groups.size(), false);
        for (const GroupId group : assignmentData.unassignedGroups) {
            groupIsUnassigned[group] = true;
        }
        entries.reserve(assignmentData.groups.size());
        std::vector<GroupList> groupsPerDemand(demands.entries.size());
        Progress progress(assignmentData.groups.size());
        progress.SetCheckTimeStep(1000);
        for (size_t i = 0; i < assignmentData.groups.size(); i++) {
            const GroupData& group = assignmentData.groups[i];
            totalNumberOfPassengers += group.groupSize;
            if (groupIsUnassigned[i]) {
                numberOfUnassignedPassengers += group.groupSize;
            } else {
                Entry<int32_t> newEntry(0);
                if (demands.entries.size() <= group.demandIndex) {
                    warning("Group ", i, " (demand ", group.demandIndex, ") corresponds to no demand!");
                    continue;
                }
                const AccumulatedVertexDemand::Entry& demand = demands.entries[group.demandIndex];
                if (groupsPerDemand.size() <= group.demandIndex) {
                    groupsPerDemand.resize(group.demandIndex + 1);
                }
                groupsPerDemand[group.demandIndex].emplace_back(i);
                if (groupIsWalking[i]) {
                    numberOfWalkingPassengers += group.groupSize;
                    dist.run(demand.originVertex, demand.destinationVertex);
                    if (!dist.reachable()) {
                        warning("Group ", i, " (demand ", group.demandIndex, ") walks, but the destination is not reachable!");
                        continue;
                    }
                    newEntry.numberOfTrips = 0;
                    newEntry.numberOfConnections = 0;
                    newEntry.timeInVehicle = 0;
                    newEntry.walkingTime = dist.getDistance();
                    newEntry.travelTimeWithoutInitialWaiting = newEntry.walkingTime;
                    newEntry.travelTimeWithInitialWaiting = newEntry.walkingTime;
                    newEntry.groupSize = group.groupSize;
                } else {
                    if (assignmentData.connectionsPerGroup[i].empty()) {
                        warning("Group ", i, " (demand ", group.demandIndex, ") drives, but no connections were used!");
                        continue;
                    }
                    Set<size_t> usedTrips;
                    Vertex currentPosition = demand.destinationVertex;
                    int currentTime = -never;
                    for (size_t j = assignmentData.connectionsPerGroup[i].size() - 1; j < assignmentData.connectionsPerGroup[i].size(); j--) {
                        const CSA::Connection& connection = data.connections[assignmentData.connectionsPerGroup[i][j]];
                        if ((j > 0) && (data.connections[assignmentData.connectionsPerGroup[i][j - 1]].arrivalTime > connection.departureTime)) {
                            warning("Group ", i, " (demand ", group.demandIndex, ", connection ", assignmentData.connectionsPerGroup[i][j], ") uses connections out of order!");
                        }
                        usedTrips.insert(connection.tripId);
                        newEntry.numberOfConnections++;
                        newEntry.timeInVehicle += connection.travelTime();
                        newEntry.travelTimeWithoutInitialWaiting += connection.travelTime();
                        newEntry.travelTimeWithInitialWaiting += connection.travelTime();
                        if (connection.arrivalStopId != currentPosition) {
                            dist.run(connection.arrivalStopId, currentPosition);
                            if (!dist.reachable()) {
                                warning("Group ", i, " (demand ", group.demandIndex, ", connection ", assignmentData.connectionsPerGroup[i][j], ", from ", connection.arrivalStopId, ", to ", currentPosition, ") walks intermediate, but the next stop is not reachable!");
                                continue;
                            }
                            const int walkingTime = dist.getDistance();
                            newEntry.walkingTime += walkingTime;
                            newEntry.travelTimeWithoutInitialWaiting += walkingTime;
                            newEntry.travelTimeWithInitialWaiting += walkingTime;
                            if (connection.arrivalTime + walkingTime < currentTime) {
                                const int waitingTime = currentTime - walkingTime - connection.arrivalTime;
                                newEntry.travelTimeWithoutInitialWaiting += waitingTime;
                                newEntry.travelTimeWithInitialWaiting += waitingTime;
                            }
                        }
                        currentPosition = connection.departureStopId;
                        currentTime = connection.departureTime;
                    }
                    if (demand.originVertex != currentPosition) {
                        dist.run(demand.originVertex, currentPosition);
                        if (!dist.reachable()) {
                            warning("Group ", i, " (demand ", group.demandIndex, ", connection ", assignmentData.connectionsPerGroup[i][0], ", from ", demand.originVertex, ", to ", currentPosition, ") walks initially, but the next stop is not reachable!");
                            continue;
                        }
                        const int walkingTime = dist.getDistance();
                        newEntry.walkingTime += walkingTime;
                        newEntry.travelTimeWithoutInitialWaiting += walkingTime;
                        newEntry.travelTimeWithInitialWaiting += walkingTime;
                        if (demand.latestDepartureTime + walkingTime < currentTime) {
                            const int waitingTime = currentTime - walkingTime - demand.latestDepartureTime;
                            newEntry.travelTimeWithInitialWaiting += waitingTime;
                        }
                    }
                    newEntry.numberOfTrips = usedTrips.size();
                    newEntry.groupSize = group.groupSize;
                }
                min.minimize(newEntry);
                max.maximize(newEntry);
                sum.add(newEntry);
                entries.emplace_back(newEntry);
            }
            progress++;
        }
        for (const GroupList& list : groupsPerDemand) {
            Set<Path> paths;
            for (const GroupId group : list) {
                paths.insert(assignmentData.connectionsPerGroup[group]);
            }
            numberOfUsedPaths += paths.size();
        }
        std::cout << std::endl;
    }

    GroupAssignmentStatistic(const CSA::Data& data, const AccumulatedVertexDemand& demands, const AssignmentData& assignmentData, const int passengerMultiplier) :
        GroupAssignmentStatistic(DijkstraDistance(data.transferGraph), data, demands, assignmentData, passengerMultiplier) {
    }

    GroupAssignmentStatistic(const CSA::Data& data, const CH::BucketQuery<>& bucketQuery, const AccumulatedVertexDemand& demands, const AssignmentData& assignmentData, const int passengerMultiplier) :
        GroupAssignmentStatistic(BucketDistance(data.transferGraph, data.numberOfStops(), bucketQuery), data, demands, assignmentData, passengerMultiplier) {
    }

    inline void printInfo() const noexcept {
        std::cout << *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const GroupAssignmentStatistic& d) {
        const double n = d.totalNumberOfPassengers - d.numberOfUnassignedPassengers;
        out << "GroupAssignmentStatistic (" << String::prettyDouble(d.totalNumberOfPassengers) << " passengers in " << String::prettyInt(d.totalNumberOfGroups) << " groups):\n";
        out << "Value                            " << std::setw(12) << "Min"                                                      << std::setw(14) << "Mean"                                                            << std::setw(14) << "Max" << "\n";
        out << "numberOfTrips:                   " << std::setw(12) << String::prettyInt(d.min.numberOfTrips)                     << std::setw(14) << String::prettyDouble(d.sum.numberOfTrips / n)                     << std::setw(14) << String::prettyInt(d.max.numberOfTrips) << "\n";
        out << "numberOfConnections:             " << std::setw(12) << String::prettyInt(d.min.numberOfConnections)               << std::setw(14) << String::prettyDouble(d.sum.numberOfConnections / n)               << std::setw(14) << String::prettyInt(d.max.numberOfConnections) << "\n";
        out << "timeInVehicle:                   " << std::setw(12) << String::secToString(d.min.timeInVehicle)                   << std::setw(14) << String::secToString(d.sum.timeInVehicle / n)                      << std::setw(14) << String::secToString(d.max.timeInVehicle) << "\n";
        out << "walkingTime:                     " << std::setw(12) << String::secToString(d.min.walkingTime)                     << std::setw(14) << String::secToString(d.sum.walkingTime / n)                        << std::setw(14) << String::secToString(d.max.walkingTime) << "\n";
        out << "travelTimeWithoutInitialWaiting: " << std::setw(12) << String::secToString(d.min.travelTimeWithoutInitialWaiting) << std::setw(14) << String::secToString(d.sum.travelTimeWithoutInitialWaiting / n)    << std::setw(14) << String::secToString(d.max.travelTimeWithoutInitialWaiting) << "\n";
        out << "travelTimeWithInitialWaiting:    " << std::setw(12) << String::secToString(d.min.travelTimeWithInitialWaiting)    << std::setw(14) << String::secToString(d.sum.travelTimeWithInitialWaiting / n)       << std::setw(14) << String::secToString(d.max.travelTimeWithInitialWaiting) << "\n";
        out << "groupSize:                       " << std::setw(12) << String::prettyInt(d.min.groupSize)                         << std::setw(14) << String::prettyDouble(d.sum.groupSize / n)                         << std::setw(14) << String::prettyInt(d.max.groupSize) << "\n";
        out << "number of walking passengers:    " << String::prettyDouble(d.numberOfWalkingPassengers) << " (" << String::percent(d.numberOfWalkingPassengers / d.totalNumberOfPassengers) << ")" << std::endl;
        out << "number of unassigned passengers: " << String::prettyDouble(d.numberOfUnassignedPassengers) << " (" << String::percent(d.numberOfUnassignedPassengers / d.totalNumberOfPassengers) << ")" << std::endl;
        out << "number of empty connections:     " << String::prettyInt(d.numberOfEmptyConnections) << " (" << String::percent(d.numberOfEmptyConnections / static_cast<double>(d.totalNumberOfConnections)) << ")" << std::endl;
        out << "number of empty trips:           " << String::prettyInt(d.numberOfEmptyTrips) << " (" << String::percent(d.numberOfEmptyTrips / static_cast<double>(d.totalNumberOfTrips)) << ")" << std::endl;
        out << "number of paths:                 " << String::prettyInt(d.numberOfUsedPaths) << " (" << String::prettyDouble(d.numberOfUsedPaths / (d.totalNumberOfPassengers / static_cast<double>(d.passengerMultiplier))) << " p.P.)" << std::endl;
        return out;
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, entries, min, max, sum, totalNumberOfTrips, totalNumberOfGroups, totalNumberOfPassengers, totalNumberOfConnections, numberOfWalkingPassengers, numberOfUnassignedPassengers, numberOfEmptyConnections, numberOfEmptyTrips, numberOfUsedPaths, passengerMultiplier);
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, entries, min, max, sum, totalNumberOfTrips, totalNumberOfGroups, totalNumberOfPassengers, totalNumberOfConnections, numberOfWalkingPassengers, numberOfUnassignedPassengers, numberOfEmptyConnections, numberOfEmptyTrips, numberOfUsedPaths, passengerMultiplier);
    }

private:
    std::vector<Entry<int32_t>> entries;
    Entry<int32_t> min;
    Entry<int32_t> max;
    Entry<int64_t> sum;

    size_t totalNumberOfTrips;
    size_t totalNumberOfGroups;
    double totalNumberOfPassengers;
    size_t totalNumberOfConnections;
    double numberOfWalkingPassengers;
    double numberOfUnassignedPassengers;
    size_t numberOfEmptyConnections;
    size_t numberOfEmptyTrips;
    size_t numberOfUsedPaths;
    int passengerMultiplier;

};

}
