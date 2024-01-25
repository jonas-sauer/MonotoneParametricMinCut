#pragma once

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>

#include "Passenger.h"
#include "IdVertexDemand.h"

#include "../CSA/Data.h"
#include "../Geometry/Point.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Console/Progress.h"

#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Container/Map.h"

#include "../../Algorithms/Dijkstra/Dijkstra.h"

class PassengerData {

public:
    static const std::string CSV_HEADER;

    template<typename TIME = int32_t>
    struct Entry {
        Entry(const int value = 0) :
            index(value),
            id(value),
            originVertex(value),
            destinationVertex(value),
            firstStop(value),
            lastStop(value),
            numberOfTrips(value),
            numberOfConnections(value),
            departureTime(value),
            arrivalTime(value),
            timeInVehicle(value),
            travelTimeWithoutInitialWaiting(value),
            travelTimeWithInitialWaiting(value),
            beelineDistanceST(value),
            beelineDistanceOD(value),
            pathDistance(value) {
        }
        inline std::ostream& toCSV(std::ostream& out) const {
            out << index << ", " << id << ", " << originVertex << ", " << destinationVertex << ", ";
            out << firstStop << ", " << lastStop << ", ";
            out << numberOfTrips << ", " << numberOfConnections << ", ";
            out << departureTime << ", " << arrivalTime << ", " << timeInVehicle << ", ";
            out << travelTimeWithoutInitialWaiting << ", " << travelTimeWithInitialWaiting << ", ";
            out << beelineDistanceST << ", " << beelineDistanceOD << ", " << pathDistance << "\n";
            return out;
        }
        inline void print() const noexcept {
            std::cout << "PassengerData:\n";
            std::cout << "index:                           " << std::setw(12) << String::prettyInt(index)                             << "\n";
            std::cout << "originVertex:                    " << std::setw(12) << originVertex                                         << "\n";
            std::cout << "destinationVertex:               " << std::setw(12) << destinationVertex                                    << "\n";
            std::cout << "firstStop:                       " << std::setw(12) << firstStop                                            << "\n";
            std::cout << "lastStop:                        " << std::setw(12) << lastStop                                             << "\n";
            std::cout << "numberOfTrips:                   " << std::setw(12) << String::prettyInt(numberOfTrips)                     << "\n";
            std::cout << "numberOfConnections:             " << std::setw(12) << String::prettyInt(numberOfConnections)               << "\n";
            std::cout << "departureTime:                   " << std::setw(12) << String::secToTime(departureTime)                     << "\n";
            std::cout << "arrivalTime:                     " << std::setw(12) << String::secToTime(arrivalTime)                       << "\n";
            std::cout << "timeInVehicle:                   " << std::setw(12) << String::secToString(timeInVehicle)                   << "\n";
            std::cout << "travelTimeWithoutInitialWaiting: " << std::setw(12) << String::secToString(travelTimeWithoutInitialWaiting) << "\n";
            std::cout << "travelTimeWithInitialWaiting:    " << std::setw(12) << String::secToString(travelTimeWithInitialWaiting)    << "\n";
            std::cout << "beelineDistanceST:               " << std::setw(12) << String::prettyInt(beelineDistanceST)                 << "\n";
            std::cout << "beelineDistanceOD:               " << std::setw(12) << String::prettyInt(beelineDistanceOD)                 << "\n";
            std::cout << "pathDistance:                    " << std::setw(12) << String::prettyInt(pathDistance)                      << "\n";
        }
        template<typename FUNCTION, typename OTHER>
        inline void apply(const OTHER& other, const FUNCTION& function) noexcept {
            index = function(index, other.index);
            id = function(id, other.id);
            originVertex = function(originVertex, other.originVertex);
            destinationVertex = function(destinationVertex, other.destinationVertex);
            firstStop = function(firstStop, other.firstStop);
            lastStop = function(lastStop, other.lastStop);
            numberOfTrips = function(numberOfTrips, other.numberOfTrips);
            numberOfConnections = function(numberOfConnections, other.numberOfConnections);
            departureTime = function(departureTime, other.departureTime);
            arrivalTime = function(arrivalTime, other.arrivalTime);
            timeInVehicle = function(timeInVehicle, other.timeInVehicle);
            travelTimeWithoutInitialWaiting = function(travelTimeWithoutInitialWaiting, other.travelTimeWithoutInitialWaiting);
            travelTimeWithInitialWaiting = function(travelTimeWithInitialWaiting, other.travelTimeWithInitialWaiting);
            beelineDistanceST = function(beelineDistanceST, other.beelineDistanceST);
            beelineDistanceOD = function(beelineDistanceOD, other.beelineDistanceOD);
            pathDistance = function(pathDistance, other.pathDistance);
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
            apply(other, [](const auto& a, const auto& b){return a + b;});
        }
        u_int32_t index{0};
        GlobalPassengerId id{0};
        Vertex originVertex{0};
        Vertex destinationVertex{0};
        StopId firstStop{0};
        StopId lastStop{0};
        int numberOfTrips{0};
        int numberOfConnections{0};
        TIME departureTime{0};
        TIME arrivalTime{0};
        TIME timeInVehicle{0};
        TIME travelTimeWithoutInitialWaiting{0};
        TIME travelTimeWithInitialWaiting{0};
        double beelineDistanceST{0};
        double beelineDistanceOD{0};
        double pathDistance{0};
    };

    struct Path {
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

private:
    PassengerData() {}

public:
    inline static PassengerData FromBinary(const std::string& filename) noexcept {
        PassengerData result;
        result.deserialize(filename);
        return result;
    }

    inline static PassengerData FromApportionment(const CSA::Data& data, const AccumulatedVertexDemand& demand, const std::vector<GlobalPassengerList>& passengersInConnection, const GlobalPassengerList& unassignedPassengers, const GlobalPassengerList& walkingPassengers, const bool allowEarlyDeparture = false, const bool checkPaths = true) noexcept {
        return FromApportionment(data, IdVertexDemand::FromAccumulatedVertexDemand(demand, 10, 300), passengersInConnection, unassignedPassengers, walkingPassengers, allowEarlyDeparture, checkPaths);
    }

    inline static PassengerData FromApportionment(const CSA::Data& data, const IdVertexDemand& demand, const std::vector<GlobalPassengerList>& passengersInConnection, const GlobalPassengerList& unassignedPassengers, const GlobalPassengerList& walkingPassengers, const bool allowEarlyDeparture = false, const bool checkPaths = true) noexcept {
        PassengerData result;
        size_t entryCount = 0;
        for (const IdVertexDemand::Entry& entry : demand.entries) {
            entryCount += entry.ids.size();
        }
        result.entries.resize(entryCount);
        result.entries.clear();
        result.passengersInConnection = passengersInConnection;
        result.unassignedPassengers = unassignedPassengers;
        result.walkingPassengers = walkingPassengers;
        result.min = Entry<>(intMax);
        result.max = Entry<>(0);
        result.sum = Entry<int64_t>(0);
        result.numberOfPaths = 0;
        result.numberOfConnections = data.numberOfConnections();
        std::vector<std::vector<std::vector<ConnectionId>>> connectionsByPassengerByDestination = result.getConnectionsByPassengerByDestination(data);
        Set<GlobalPassengerId> passengersWithoutConnection;
        for (const GlobalPassengerId id : unassignedPassengers) {
            passengersWithoutConnection.insert(id);
        }
        for (const GlobalPassengerId id : walkingPassengers) {
            passengersWithoutConnection.insert(id);
        }
        Dijkstra<TransferGraph> dijkstra(data.transferGraph);
        std::set<Path> paths;
        Progress progress(entryCount);
        progress.SetCheckTimeStep(1000);
        for (const IdVertexDemand::Entry& entry : demand.entries) {
            const Vertex destinationVertex = entry.destinationVertex;
            Assert((destinationVertex < connectionsByPassengerByDestination.size()), "There is no data for destination vertex " << destinationVertex << " (last destination: " << (connectionsByPassengerByDestination.size() - 1) << ")!");
            std::vector<std::vector<ConnectionId>>& connectionsByPassenger = connectionsByPassengerByDestination[destinationVertex];
            for (const DestinationSpecificPassengerId passenger : entry.ids) {
                Entry<> newEntry;
                newEntry.index = result.entries.size();
                newEntry.id = getGlobalPassengerId(destinationVertex, passenger);
                newEntry.originVertex = entry.originVertex;
                newEntry.destinationVertex = destinationVertex;
                newEntry.beelineDistanceOD = Geometry::geoDistanceInCM(data.transferGraph.get(Coordinates, entry.originVertex), data.transferGraph.get(Coordinates, destinationVertex)) / 100.0;
                if ((passenger >= connectionsByPassenger.size()) || (connectionsByPassenger[passenger].empty())) {
                    Assert(passengersWithoutConnection.contains(getGlobalPassengerId(destinationVertex, passenger)), "Passenger (destination vertex: " << destinationVertex << ", passenger id: " << passenger << ", global id: " << getGlobalPassengerId(destinationVertex, passenger) << ") is neither unassigned nor does he use any connections!");
                    newEntry.numberOfConnections = 0;
                    newEntry.departureTime = entry.departureTime;
                    newEntry.arrivalTime = entry.departureTime + getTravelTime<false>(data, dijkstra, entry.originVertex, destinationVertex);
                    newEntry.beelineDistanceST = 0;
                    newEntry.timeInVehicle = 0;
                    newEntry.pathDistance = 0;
                    newEntry.numberOfTrips = 0;
                    Path path;
                    path.data.emplace_back(entry.originVertex);
                    path.data.emplace_back(destinationVertex);
                    paths.insert(path);
                } else {
                    if (checkPaths) Assert(result.isValidPath(data, connectionsByPassenger[passenger], entry, allowEarlyDeparture), "Passenger (destination vertex: " << destinationVertex << ", passenger id: " << passenger << ", global id: " << getGlobalPassengerId(destinationVertex, passenger) << ") has a path that does not comply with his demand!");
                    Assert(!passengersWithoutConnection.contains(getGlobalPassengerId(destinationVertex, passenger)), "Passenger (destination vertex: " << destinationVertex << ", passenger id: " << passenger << ", global id: " << getGlobalPassengerId(destinationVertex, passenger) << ") is unassigned and uses some connections!");
                    std::vector<ConnectionId>& connections = connectionsByPassenger[passenger];
                    newEntry.firstStop = data.connections[connections.front()].departureStopId;
                    newEntry.lastStop = data.connections[connections.back()].arrivalStopId;
                    newEntry.numberOfConnections = connections.size();
                    newEntry.departureTime = data.connections[connections.front()].departureTime - getTravelTime<true>(data, dijkstra, entry.originVertex, data.connections[connections.front()].departureStopId);
                    newEntry.arrivalTime = data.connections[connections.back()].arrivalTime + getTravelTime<true>(data, dijkstra, data.connections[connections.back()].arrivalStopId, destinationVertex);
                    newEntry.travelTimeWithoutInitialWaiting = data.connections[connections.back()].arrivalTime - data.connections[connections.front()].departureTime;
                    newEntry.travelTimeWithInitialWaiting = data.connections[connections.back()].arrivalTime - entry.departureTime;
                    newEntry.beelineDistanceST = Geometry::geoDistanceInCM(data.stopData[newEntry.firstStop].coordinates, data.stopData[newEntry.lastStop].coordinates) / 100.0;
                    newEntry.timeInVehicle = 0;
                    newEntry.pathDistance = 0;
                    std::set<TripId> trips;
                    Path path;
                    path.data.emplace_back(entry.originVertex);
                    for (const ConnectionId i : connections) {
                        path.data.emplace_back(i);
                        trips.insert(data.connections[i].tripId);
                        newEntry.timeInVehicle += data.connections[i].travelTime();
                        newEntry.pathDistance += Geometry::geoDistanceInCM(data.stopData[data.connections[i].departureStopId].coordinates, data.stopData[data.connections[i].arrivalStopId].coordinates) / 100.0;
                    }
                    path.data.emplace_back(destinationVertex);
                    path.data.shrink_to_fit();
                    paths.insert(path);
                    newEntry.numberOfTrips = trips.size();
                }
                Assert(newEntry.arrivalTime >= newEntry.departureTime, "newEntry.arrivalTime >= newEntry.departureTime (" << newEntry.arrivalTime << " >= " << newEntry.departureTime << ")!");
                newEntry.travelTimeWithoutInitialWaiting = newEntry.arrivalTime - newEntry.departureTime;
                Assert(newEntry.arrivalTime >= entry.departureTime, "newEntry.arrivalTime >= entry.departureTime (" << newEntry.arrivalTime << " >= " << entry.departureTime << ")!");
                newEntry.travelTimeWithInitialWaiting = newEntry.arrivalTime - entry.departureTime;
                result.min.minimize(newEntry);
                result.max.maximize(newEntry);
                result.sum.add(newEntry);
                result.entries.emplace_back(newEntry);
                progress++;
            }
            result.numberOfPaths += paths.size();
            paths.clear();
        }
        std::cout << std::endl;
        for (const ConnectionId connection : range(ConnectionId(passengersInConnection.size()))) {
            if (passengersInConnection[connection].empty()) {
                result.emptyConnections.push_back(connection);
            }
        }
        return result;
    }

private:
    template<bool ENSURE_EDGE_EXISTS>
    inline static int getTravelTime(const CSA::Data& data, const Vertex from, const Vertex to) noexcept {
        Assert(data.transferGraph.isVertex(from), "Invalid vertex id: " << from << "!");
        Assert(data.transferGraph.isVertex(to), "Invalid vertex id: " << to << "!");
        if (from == to) {
            return 0;
        } else {
            const Edge edge = data.transferGraph.findEdge(from, to);
            if (ENSURE_EDGE_EXISTS) {
                Assert(data.transferGraph.isEdge(edge), "The edge from " << from << " to " << to << " is missing!");
                return data.transferGraph.get(TravelTime, edge);
            } else if (data.transferGraph.isEdge(edge)) {
                return data.transferGraph.get(TravelTime, edge);
            } else {
                return 0;
            }
        }
    }

    template<bool ENSURE_EDGE_EXISTS>
    inline static int getTravelTime(const CSA::Data& data, Dijkstra<TransferGraph>& dijkstra, const Vertex from, const Vertex to) noexcept {
        suppressUnusedParameterWarning(dijkstra);
        Assert(data.transferGraph.isVertex(from), "Invalid vertex id: " << from << "!");
        Assert(data.transferGraph.isVertex(to), "Invalid vertex id: " << to << "!");
        if (from == to) {
            return 0;
        } else {
            const Edge edge = data.transferGraph.findEdge(from, to);
            if (data.transferGraph.isEdge(edge)) {
                return data.transferGraph.get(TravelTime, edge);
            } else {
                dijkstra.run(from, to);
                if constexpr (ENSURE_EDGE_EXISTS) {
                    Assert(dijkstra.reachable(to), "The path from " << from << " to " << to << " is missing!");
                    return dijkstra.getDistance(to);
                } else {
                    if (dijkstra.reachable(to)) {
                        return dijkstra.getDistance(to);
                    } else {
                        return 0;
                    }
                }
            }
        }
    }

public:
    inline std::vector<std::vector<std::vector<ConnectionId>>> getConnectionsByPassengerByDestination(const CSA::Data& data) const noexcept {
        std::vector<std::vector<std::vector<ConnectionId>>> connectionsByPassengerByDestination(data.transferGraph.numVertices());
        for (const ConnectionId i : range(ConnectionId(data.numberOfConnections()))) {
            for (const GlobalPassengerId id : passengersInConnection[i]) {
                const u_int32_t destination = getDestination(id);
                if (destination >= connectionsByPassengerByDestination.size()) connectionsByPassengerByDestination.resize(destination + 1);
                std::vector<std::vector<ConnectionId>>& connectionsByPassenger = connectionsByPassengerByDestination[destination];
                const DestinationSpecificPassengerId passenger = getDestinationSpecificPassengerId(id);
                if (passenger >= connectionsByPassenger.size()) connectionsByPassenger.resize(passenger + 1);
                connectionsByPassenger[passenger].emplace_back(i);
            }
        }
        connectionsByPassengerByDestination.shrink_to_fit();
        for (std::vector<std::vector<ConnectionId>>& connectionsByPassenger : connectionsByPassengerByDestination) {
            connectionsByPassenger.shrink_to_fit();
            for (std::vector<ConnectionId>& connections : connectionsByPassenger) {
                connections.shrink_to_fit();
            }
        }
        return connectionsByPassengerByDestination;
    }

    inline std::vector<StopId> getPath(const CSA::Data& data, const u_int32_t index) const noexcept {
        std::vector<StopId> path;
        const GlobalPassengerId id = entries[index].id;
        for (const ConnectionId i : data.connectionIds()) {
            if (!Vector::contains(passengersInConnection[i], id)) continue;
            const CSA::Connection& connection = data.connections[i];
            path.emplace_back(connection.departureStopId);
            path.emplace_back(connection.arrivalStopId);
        }
        return path;
    }

    inline static bool isValidPath(const CSA::Data& data, const std::vector<ConnectionId>& path, const IdVertexDemand::Entry& demand, const bool  allowEarlyDeparture = false) noexcept {
        if (path.empty()) {
            return true;
        } else {
            Assert(data.isConnection(path.front()), "First connection id " << path.front() << " does not represent a connection!");
            Assert(data.isConnection(path.back()), "Last connection id " << path.back() << " does not represent a connection!");
            if (!data.isCombinable(demand.originVertex, ((allowEarlyDeparture) ? (-intMax) : (demand.departureTime)), data.connections[path.front()])) return false;
            if (!data.isCombinable(data.connections[path.back()], demand.destinationVertex)) return false;
            for (const size_t i : range(path.size() - 1)) {
                Assert(data.isConnection(path[i]), "" << i << "th connection id " << path[i] << " does not represent a connection!");
                Assert(data.isConnection(path[i + 1]), "" << i << "th connection id " << path[i] << " does not represent a connection!");
                if (!data.isCombinable(data.connections[path[i]], data.connections[path[i + 1]])) return false;
            }
            return true;
        }
    }

    inline static void printPath(const CSA::Data& data, const std::vector<ConnectionId>& path, const IdVertexDemand::Entry& demand) noexcept {
        std::cout << "starting at vertex " << demand.originVertex << " at " << String::secToTime(demand.departureTime) << std::endl;
        if (path.empty()) {
            std::cout << "   walk direct to the destination" << std::endl;
        } else {
            std::cout << "   walk to stop " << data.connections[path.front()].departureStopId << std::endl;
            int departureTime = demand.departureTime;
            if (demand.originVertex == data.connections[path.front()].departureStopId) {
                std::cout << "      which is the origin!" << std::endl;
            } else {
                Edge edge = data.transferGraph.findEdge(demand.originVertex, data.connections[path.front()].departureStopId);
                if (!data.transferGraph.isEdge(edge)) {
                    std::cout << "      which is impossible since the edge does not exist!" << std::endl;
                } else {
                    departureTime += data.transferGraph.get(TravelTime, edge);
                    std::cout << "      taking      " << std::setw(13) << String::secToString(data.transferGraph.get(TravelTime, edge)) << std::endl;
                    std::cout << "      arriving at " << std::setw(13) << String::secToTime(departureTime) << std::endl;
                }
            }
            std::cout << "   take connection " << std::setw(13) << path.front() << std::endl;
            std::cout << "      to stop      " << std::setw(13) << data.connections[path.front()].arrivalStopId << std::endl;
            std::cout << "      departing at " << std::setw(13) << String::secToTime(data.connections[path.front()].departureTime) << std::endl;
            std::cout << "      arriving at  " << std::setw(13) << String::secToTime(data.connections[path.front()].arrivalTime) << std::endl;
            if (data.connections[path.front()].departureTime < departureTime) {
                std::cout << "      which is impossible since the connection departs to early!" << std::endl;
            }
            if (!data.isCombinable(demand.originVertex, demand.departureTime, data.connections[path.front()])) {
                std::cout << "      which is impossible since the connection is not combinable!" << std::endl;
            }
            for (const size_t i : range(path.size() - 1)) {
                if (data.connections[path[i]].tripId == data.connections[path[i + 1]].tripId) {
                } else if (data.connections[path[i]].arrivalStopId == data.connections[path[i + 1]].departureStopId) {
                    std::cout << "   wait at stop " << data.connections[path[i + 1]].departureStopId << std::endl;
                    std::cout << "      taking      " << std::setw(13) << String::secToString(data.minTransferTime(data.connections[path[i + 1]].departureStopId)) << std::endl;
                    std::cout << "      arriving at " << std::setw(13) << String::secToTime(data.connections[path[i]].arrivalTime + data.minTransferTime(data.connections[path[i + 1]].departureStopId)) << std::endl;
                } else {
                    std::cout << "   walk to stop " << data.connections[path[i + 1]].departureStopId << std::endl;
                    Edge edge = data.transferGraph.findEdge(data.connections[path[i]].arrivalStopId, data.connections[path[i + 1]].departureStopId);
                    if (!data.transferGraph.isEdge(edge)) {
                        std::cout << "      which is impossible since the edge does not exist!" << std::endl;
                    } else {
                        std::cout << "      taking      " << std::setw(13) << String::secToString(data.transferGraph.get(TravelTime, edge)) << std::endl;
                        std::cout << "      arriving at " << std::setw(13) << String::secToTime(data.connections[path[i]].arrivalTime + data.transferGraph.get(TravelTime, edge)) << std::endl;
                    }
                }
                std::cout << "   take connection " << std::setw(13) << path[i + 1] << std::endl;
                std::cout << "      to stop      " << std::setw(13) << data.connections[path[i + 1]].arrivalStopId << std::endl;
                std::cout << "      departing at " << std::setw(13) << String::secToTime(data.connections[path[i + 1]].departureTime) << std::endl;
                std::cout << "      arriving at  " << std::setw(13) << String::secToTime(data.connections[path[i + 1]].arrivalTime) << std::endl;
                if (!data.isCombinable(data.connections[path[i]], data.connections[path[i + 1]])) {
                    std::cout << "      which is impossible since the connection is not combinable!" << std::endl;
                }
            }
            std::cout << "   walk to destination" << std::endl;
            if (data.connections[path.back()].arrivalStopId == demand.destinationVertex) {
                std::cout << "      which is the last stop!" << std::endl;
            } else {
                Edge edge = data.transferGraph.findEdge(data.connections[path.back()].arrivalStopId, demand.destinationVertex);
                if (!data.transferGraph.isEdge(edge)) {
                    std::cout << "      which is impossible since the edge does not exist!" << std::endl;
                } else {
                    std::cout << "      taking       " << std::setw(13) << String::secToString(data.transferGraph.get(TravelTime, edge)) << std::endl;
                    std::cout << "      arriving at  " << std::setw(13) << String::secToTime(data.connections[path.back()].arrivalTime + data.transferGraph.get(TravelTime, edge)) << std::endl;
                }
            }
            if (!data.isCombinable(data.connections[path.back()], demand.destinationVertex)) {
                std::cout << "      which is impossible since the connection is not combinable!" << std::endl;
            }
        }
        std::cout << "arrive at destination vertex " << demand.destinationVertex << std::endl;
    }

    inline bool usesPublicTransit(const GlobalPassengerId id) const noexcept {
        return (!Vector::contains(walkingPassengers, id)) && (!Vector::contains(unassignedPassengers, id));
    }

    inline void printInfo() const noexcept {
        std::cout << *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const PassengerData& d) {
        const double n = d.entries.size();
        out << "PassengerData (" << String::prettyInt(n) << " entries):\n";
        out << "Value                            " << std::setw(12) << "Min"                                                      << std::setw(15) << "Mean"                                                            << std::setw(12) << "Max" << "\n";
        out << "index:                           " << std::setw(12) << String::prettyInt(d.min.index)                             << std::setw(15) << String::prettyDouble(d.sum.index / n)                             << std::setw(12) << String::prettyInt(d.max.index) << "\n";
        out << "originVertex:                    " << std::setw(12) << d.min.originVertex                                         << std::setw(15) << String::prettyDouble(d.sum.originVertex / n)                      << std::setw(12) << d.max.originVertex << "\n";
        out << "destinationVertex:               " << std::setw(12) << d.min.destinationVertex                                    << std::setw(15) << String::prettyDouble(d.sum.destinationVertex / n)                 << std::setw(12) << d.max.destinationVertex << "\n";
        out << "firstStop:                       " << std::setw(12) << d.min.firstStop                                            << std::setw(15) << String::prettyDouble(d.sum.firstStop / n)                         << std::setw(12) << d.max.firstStop << "\n";
        out << "lastStop:                        " << std::setw(12) << d.min.lastStop                                             << std::setw(15) << String::prettyDouble(d.sum.lastStop / n)                          << std::setw(12) << d.max.lastStop << "\n";
        out << "numberOfTrips:                   " << std::setw(12) << String::prettyInt(d.min.numberOfTrips)                     << std::setw(15) << String::prettyDouble(d.sum.numberOfTrips / n)                     << std::setw(12) << String::prettyInt(d.max.numberOfTrips) << "\n";
        out << "numberOfConnections:             " << std::setw(12) << String::prettyInt(d.min.numberOfConnections)               << std::setw(15) << String::prettyDouble(d.sum.numberOfConnections / n)               << std::setw(12) << String::prettyInt(d.max.numberOfConnections) << "\n";
        out << "departureTime:                   " << std::setw(12) << String::secToTime(d.min.departureTime)                     << std::setw(15) << String::secToTime(d.sum.departureTime / n)                        << std::setw(12) << String::secToTime(d.max.departureTime) << "\n";
        out << "arrivalTime:                     " << std::setw(12) << String::secToTime(d.min.arrivalTime)                       << std::setw(15) << String::secToTime(d.sum.arrivalTime / n)                          << std::setw(12) << String::secToTime(d.max.arrivalTime) << "\n";
        out << "timeInVehicle:                   " << std::setw(12) << String::secToString(d.min.timeInVehicle)                   << std::setw(15) << String::secToString(d.sum.timeInVehicle / n)                      << std::setw(12) << String::secToString(d.max.timeInVehicle) << "\n";
        out << "travelTimeWithoutInitialWaiting: " << std::setw(12) << String::secToString(d.min.travelTimeWithoutInitialWaiting) << std::setw(15) << String::secToString(d.sum.travelTimeWithoutInitialWaiting / n)    << std::setw(12) << String::secToString(d.max.travelTimeWithoutInitialWaiting) << "\n";
        out << "travelTimeWithInitialWaiting:    " << std::setw(12) << String::secToString(d.min.travelTimeWithInitialWaiting)    << std::setw(15) << String::secToString(d.sum.travelTimeWithInitialWaiting / n)       << std::setw(12) << String::secToString(d.max.travelTimeWithInitialWaiting) << "\n";
        out << "beelineDistanceST:               " << std::setw(12) << String::prettyInt(d.min.beelineDistanceST)                 << std::setw(15) << String::prettyDouble(d.sum.beelineDistanceST / n)                 << std::setw(12) << String::prettyInt(d.max.beelineDistanceST) << "\n";
        out << "beelineDistanceOD:               " << std::setw(12) << String::prettyInt(d.min.beelineDistanceOD)                 << std::setw(15) << String::prettyDouble(d.sum.beelineDistanceOD / n)                 << std::setw(12) << String::prettyInt(d.max.beelineDistanceOD) << "\n";
        out << "pathDistance:                    " << std::setw(12) << String::prettyInt(d.min.pathDistance)                      << std::setw(15) << String::prettyDouble(d.sum.pathDistance / n)                      << std::setw(12) << String::prettyInt(d.max.pathDistance) << "\n";
        out << "number of walking passengers:    " << String::prettyInt(d.walkingPassengers.size()) << " (" << String::percent(d.walkingPassengers.size() / (double)(d.entries.size())) << ")" << std::endl;
        out << "number of unassigned passengers: " << String::prettyInt(d.unassignedPassengers.size()) << " (" << String::percent(d.unassignedPassengers.size() / (double)(d.entries.size())) << ")" << std::endl;
        out << "number of empty connections:     " << String::prettyInt(d.emptyConnections.size()) << " (" << String::percent(d.emptyConnections.size() / (double)(d.numberOfConnections)) << ")" << std::endl;
        out << "number of paths:                 " << String::prettyInt(d.numberOfPaths) << std::endl;
        return out;
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, passengersInConnection, unassignedPassengers, walkingPassengers, entries, min, max, sum, emptyConnections, numberOfConnections, numberOfPaths);
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, passengersInConnection, unassignedPassengers, walkingPassengers, entries, min, max, sum, emptyConnections, numberOfConnections, numberOfPaths);
    }

    inline void writePassengerConnectionPairs(const CSA::Data& data, const IdVertexDemand& demand, const std::string& fileName) const noexcept {
        std::ofstream os(fileName);
        Assert(os, "cannot open file: " << fileName);
        Assert(os.is_open(), "cannot open file: " << fileName);
        os << "connection_id,passenger_id\n";
        int passengerId = 0;
        std::vector<std::vector<std::vector<ConnectionId>>> connectionsByPassengerByDestination = getConnectionsByPassengerByDestination(data);
        for (const IdVertexDemand::Entry& entry : demand.entries) {
            std::vector<std::vector<ConnectionId>>& connectionsByPassenger = connectionsByPassengerByDestination[entry.destinationVertex];
            for (const DestinationSpecificPassengerId passenger : entry.ids) {
                if ((passenger < connectionsByPassenger.size()) && (!connectionsByPassenger[passenger].empty())) {
                    for (const ConnectionId connection : connectionsByPassenger[passenger]) {
                        os << connection.value() << "," << passengerId << "\n";
                    }
                }
                passengerId++;
            }
        }
    }

    inline void writeCumulativeStopDemand(const CSA::Data& data, const std::string& fileName) const noexcept {
        std::vector<Map<StopId, size_t>> demandBySourceStop(data.numberOfStops());
        for (const Entry<>& entry : entries) {
            if (!data.isStop(entry.firstStop) || !data.isStop(entry.lastStop)) continue;
            if (!demandBySourceStop[entry.firstStop].contains(entry.lastStop)) {
                demandBySourceStop[entry.firstStop][entry.lastStop] = 1;
            } else {
                demandBySourceStop[entry.firstStop][entry.lastStop]++;
            }
        }
        std::ofstream os(fileName);
        Assert(os, "cannot open file: " << fileName);
        Assert(os.is_open(), "cannot open file: " << fileName);
        os << "sourceId,targetId,sourceLat,sourceLon,targeLat,targeLon,passengerCount\n";
        for (const StopId source : data.stops()) {
            for (std::pair<StopId, size_t> entry : demandBySourceStop[source]) {
                os << source << "," << entry.first;
                os << "," << data.stopData[source].coordinates.latitude << "," << data.stopData[source].coordinates.longitude;
                os << "," << data.stopData[entry.first].coordinates.latitude << "," << data.stopData[entry.first].coordinates.longitude;
                os << "," << entry.second;
                os << "\n";
            }

        }
    }

    inline std::ostream& toCSV(std::ostream& out) const {
        out << CSV_HEADER << "\n";
        for (const Entry<>& entry : entries) {
            entry.toCSV(out);
        }
        return out;
    }

    inline void toCSV(const std::string& fileName) const {
        std::ofstream os(fileName);
        Assert(os, "cannot open file: " << fileName);
        Assert(os.is_open(), "cannot open file: " << fileName);
        toCSV(os);
    }

public:
    std::vector<GlobalPassengerList> passengersInConnection;
    GlobalPassengerList unassignedPassengers;
    GlobalPassengerList walkingPassengers;
    std::vector<Entry<>> entries;
    Entry<> min;
    Entry<> max;
    Entry<int64_t> sum;
    std::vector<ConnectionId> emptyConnections;
    size_t numberOfConnections;
    size_t numberOfPaths;

};

const std::string PassengerData::CSV_HEADER = "index,id,origin,destination,firstStop,lastStop,numberOfTrips,numberOfConnections,timeInVehicle,travelTimeWithoutInitialWaiting,travelTimeWithInitialWaiting,beelineDistance,pathDistance";
