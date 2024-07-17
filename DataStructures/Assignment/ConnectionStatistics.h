#pragma once

#include <string>
#include <vector>

#include "../CSA/Data.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Ranges/Range.h"
#include "AssignmentData.h"
#include "GroupData.h"
#include "Settings.h"

namespace Assignment {

struct ConnectionAggregateData {
    ConnectionAggregateData(const double passengers, const int tripLength, const int tripTime, const int time) :
        passengers(passengers),
        tripLength(tripLength),
        tripTime(tripTime),
        time(time) {
    }

    inline bool operator<(const ConnectionAggregateData& other) const noexcept {
        if (tripLength > other.tripLength) return true;
        if (tripLength < other.tripLength) return false;
        if (tripTime < other.tripTime) return true;
        if (tripTime > other.tripTime) return false;
        return (time < other.time);
    }

    double passengers;
    int tripLength;
    int tripTime;
    int time;
};

struct JourneyAggregateData {
    JourneyAggregateData(const double passengers = 0.0) :
        passengers(passengers),
        travelTime(0.0),
        waitingTime(0.0),
        perceivedTravelTime(0.0),
        numberOfTransfers(0) {
    }

    inline JourneyAggregateData& operator+=(const JourneyAggregateData& other) noexcept {
        passengers += other.passengers;
        travelTime += other.travelTime;
        waitingTime += other.waitingTime;
        perceivedTravelTime += other.perceivedTravelTime;
        return *this;
    }

    inline void print(std::stringstream& ss, const int passengerMultiplier) const noexcept {
        ss << (passengers / passengerMultiplier) << ",";
        const double scalingFactor = ((passengers == 0.0) ? 1 : passengers) * 60;
        ss << (travelTime / scalingFactor) << ",";
        ss << (perceivedTravelTime / scalingFactor) << ",";
        ss << (waitingTime / scalingFactor);
    }

    double passengers;
    double travelTime;
    double waitingTime;
    double perceivedTravelTime;
    size_t numberOfTransfers;
};

class ConnectionStatistics {
public:
    ConnectionStatistics(const CSA::Data& data, const Settings& settings, const AssignmentData& assignmentData, const std::vector<double> passengerCountsPerConnection) :
        data(data),
        settings(settings),
        assignmentData(assignmentData) {
        std::vector<int> tripLength(data.numberOfTrips(), 0);
        std::vector<int> tripTime(data.numberOfTrips(), INFTY);
        for (const ConnectionId i : data.connectionIds()) {
            const CSA::Connection& connection = data.connections[i];
            tripLength[connection.tripId]++;
            tripTime[connection.tripId] = std::min(tripTime[connection.tripId], connection.departureTime);
        }
        for (const ConnectionId i : data.connectionIds()) {
            const CSA::Connection& connection = data.connections[i];
            connectionAggregateData.emplace_back(passengerCountsPerConnection[i], tripLength[connection.tripId], tripTime[connection.tripId], connection.departureTime);
        }
        std::sort(connectionAggregateData.begin(), connectionAggregateData.end());
    }

    inline void write(const std::string& fileName, const std::string& prefix) const noexcept {
        const bool headerNeeded = !FileSystem::isFile(fileName);
        std::ofstream resultFile(fileName, std::ios_base::app);
        if (headerNeeded) resultFile << getPrefixHeader() << "," << getAggregateHeader() << "," << getConnectionHeader() << "\n";
        resultFile << prefix << "," << getAggregateText() << "," << getConnectionText() << "\n";
        resultFile.close();
    }

    inline void writeAggregateText(const std::string& fileName, const std::string& prefix) const noexcept {
        std::ofstream resultFile(fileName, std::ios_base::app);
        resultFile << prefix << "," << getAggregateText() << "\n";
        resultFile.close();
    }

private:
    inline std::string getPrefixHeader() const noexcept {
        return "S,FPZ,TI,dAbf";
    }

    inline std::string getAggregateHeader() const noexcept {
        return "BusPers,BusRZ,BusERZ,BusUWZ,BusZugPers,BusZugRZ,BusZugERZ,BusZugUWZ";
    }

    inline JourneyAggregateData getJourneyAggregateData(const size_t group) const noexcept {
        JourneyAggregateData journeyAggregateData(assignmentData.groups[group].groupSize);
        const std::vector<ConnectionId>& connections = assignmentData.connectionsPerGroup[group];
        if (!connections.empty()) {
            journeyAggregateData.travelTime = data.connections[connections.back()].arrivalTime - data.connections[connections.front()].departureTime;
        }
        for (size_t i = 0; i < connections.size() - 1; i++) {
            const CSA::Connection& current = data.connections[connections[i]];
            const CSA::Connection& next = data.connections[connections[i + 1]];
            if (current.tripId == next.tripId) continue;
            journeyAggregateData.numberOfTransfers++;
            journeyAggregateData.waitingTime += next.departureTime - current.arrivalTime;
        }
        journeyAggregateData.numberOfTransfers = std::min(size_t(2), journeyAggregateData.numberOfTransfers);
        journeyAggregateData.perceivedTravelTime = journeyAggregateData.travelTime + (settings.transferCosts * journeyAggregateData.numberOfTransfers) + (settings.waitingCosts * journeyAggregateData.waitingTime);
        journeyAggregateData.travelTime *= journeyAggregateData.passengers;
        journeyAggregateData.waitingTime *= journeyAggregateData.passengers;
        journeyAggregateData.perceivedTravelTime *= journeyAggregateData.passengers;
        return journeyAggregateData;
    }

    inline std::string getAggregateText() const noexcept {
        std::stringstream result;
        std::vector<JourneyAggregateData> journeyAggregateData(3);
        for (const size_t group : indices(assignmentData.groups)) {
            const JourneyAggregateData groupData = getJourneyAggregateData(group);
            journeyAggregateData[groupData.numberOfTransfers] += groupData;
        }
        journeyAggregateData[0].print(result, settings.passengerMultiplier);
        result << ",";
        journeyAggregateData[1].print(result, settings.passengerMultiplier);
        return result.str();
    }

    inline std::string getConnectionHeader() const noexcept {
        std::stringstream result;
        Vector::printConciseMapped(connectionAggregateData, [&](const ConnectionAggregateData& data) {
            return "C" + String::secToTime(data.time);
        }, result, ",");
        return result.str();
    }

    inline std::string getConnectionText() const noexcept {
        std::stringstream result;
        Vector::printConciseMapped(connectionAggregateData, [&](const ConnectionAggregateData& data) {
            return data.passengers;
        }, result, ",");
        return result.str();
    }

private:
    const CSA::Data& data;
    const Settings& settings;
    const AssignmentData& assignmentData;
    std::vector<ConnectionAggregateData> connectionAggregateData;
};

}
