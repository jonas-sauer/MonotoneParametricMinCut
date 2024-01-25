#pragma once

#include <vector>

#include "AssignmentData.h"
#include "GroupData.h"
#include "Settings.h"
#include "../CSA/Data.h"
#include "../Demand/AccumulatedVertexDemand.h"
#include "../../Helpers/Types.h"

namespace Assignment {

struct JourneyLeg {
    JourneyLeg() :
        trip(noTripId),
        firstConnection(noConnection),
        lastConnection(noConnection),
        departureStop(noStop),
        arrivalStop(noStop),
        departureTime(never),
        arrivalTime(never) {
    }

    TripId trip;
    ConnectionId firstConnection;
    ConnectionId lastConnection;
    StopId departureStop;
    StopId arrivalStop;
    int departureTime;
    int arrivalTime;

    inline bool operator<(const JourneyLeg& other) const noexcept {
        if (trip != other.trip) return trip < other.trip;
        if (firstConnection != other.firstConnection) return firstConnection < other.firstConnection;
        return lastConnection < other.lastConnection;
    }

    inline bool operator==(const JourneyLeg& other) const noexcept {
        return trip == other.trip && firstConnection == other.firstConnection && lastConnection == other.lastConnection;
    }
};

struct JourneyInfo {
    JourneyInfo(const double groupSize = 0.0) :
        groupSize(groupSize),
        departureTime(0),
        arrivalTime(0),
        travelTime(0),
        transferWaitingTime(0),
        perceivedTravelTime(0) {
    }

    double groupSize;
    std::vector<JourneyLeg> legs;
    std::vector<std::string> tripNames;
    std::vector<std::string> departureStopNames;
    std::vector<std::string> arrivalStopNames;
    int departureTime;
    int arrivalTime;
    int travelTime;
    int transferWaitingTime;
    int perceivedTravelTime;

    inline size_t numberOfTransfers() const noexcept {
        return legs.empty() ? 0 : legs.size() - 1;
    }

    inline friend bool operator<(const JourneyInfo& a, const JourneyInfo& b) noexcept {
        return std::lexicographical_compare(a.legs.begin(), a.legs.end(), b.legs.begin(), b.legs.end());
    }
};

class JourneyWriter {
public:
    JourneyWriter(const CSA::Data& data, const Settings& settings, const AccumulatedVertexDemand& demand, const AssignmentData& assignmentData) :
        data(data),
        settings(settings),
        demand(demand),
        assignmentData(assignmentData) {
    }

    inline void write(const std::string& fileName) const noexcept {
        std::vector<JourneyInfo> journeys;
        for (const size_t group : indices(assignmentData.groups)) {
            journeys.emplace_back(getJourneyInfo(group));

        }
        std::sort(journeys.begin(), journeys.end());
        std::vector<JourneyInfo> mergedJourneys;
        for (JourneyInfo& journey : journeys) {
            if (!mergedJourneys.empty() && Vector::equals(journey.legs, mergedJourneys.back().legs)) {
                mergedJourneys.back().groupSize += journey.groupSize;
            } else {
                mergedJourneys.emplace_back(journey);
            }
        }

        IO::OFStream file(fileName);
        file << "numberOfPassengers,trips,departureStops,arrivalStops,numberOfTrips,departureTime,arrivalTime,travelTime,transferWaitingTime,perceivedTravelTime\n";
        for (const JourneyInfo& journey : mergedJourneys) {
            file << (journey.groupSize / static_cast<double>(settings.passengerMultiplier)) << ",";
            Vector::printConcise(journey.tripNames, file, "|");
            file << ",";
            Vector::printConcise(journey.departureStopNames, file, "|");
            file << ",";
            Vector::printConcise(journey.arrivalStopNames, file, "|");
            file << "," << journey.legs.size() << ",";
            file << String::secToTime(journey.departureTime) << ",";
            file << String::secToTime(journey.arrivalTime) << ",";
            file << journey.travelTime << ",";
            file << journey.transferWaitingTime << ",";
            file << journey.perceivedTravelTime << "\n";
        }
    }

private:
    inline int getWalkingTime(const Vertex from, const Vertex to) const noexcept {
        if (from == to) return 0;
        const Edge edge = data.transferGraph.findEdge(from, to);
        if (edge == noEdge) return INFTY;
        return data.transferGraph.get(TravelTime, edge);
    }

    inline JourneyInfo getJourneyInfo(const size_t group) const noexcept {
        JourneyInfo journey(assignmentData.groups[group].groupSize);
        JourneyLeg leg;
        for (const ConnectionId connection : assignmentData.connectionsPerGroup[group]) {
            if (data.connections[connection].tripId != leg.trip) {
                if (leg.trip != noTripId) {
                    journey.legs.emplace_back(leg);
                }
                leg.trip = data.connections[connection].tripId;
                leg.firstConnection = ConnectionId(connection);
            }
            leg.lastConnection = ConnectionId(connection);
        }
        if (leg.trip != noTripId) {
            journey.legs.emplace_back(leg);
        }

        for (JourneyLeg& leg : journey.legs) {
            leg.departureStop = data.connections[leg.firstConnection].departureStopId;
            leg.arrivalStop = data.connections[leg.lastConnection].arrivalStopId;
            leg.departureTime = data.connections[leg.firstConnection].departureTime;
            leg.arrivalTime = data.connections[leg.lastConnection].arrivalTime;
            journey.tripNames.emplace_back(data.tripData[leg.trip].tripName);
            journey.departureStopNames.emplace_back(data.stopData[leg.departureStop].name);
            journey.arrivalStopNames.emplace_back(data.stopData[leg.arrivalStop].name);
        }

        const AccumulatedVertexDemand::Entry& demandEntry = demand.entries[assignmentData.groups[group].demandIndex];
        if (journey.legs.empty()) {
            journey.departureTime = demandEntry.earliestDepartureTime;
            journey.travelTime = getWalkingTime(demandEntry.originVertex, demandEntry.destinationVertex);
            journey.arrivalTime = journey.departureTime + journey.travelTime;
        } else {
            const int firstDepartureTime = journey.legs[0].departureTime;
            const int initialWalkingTime = getWalkingTime(demandEntry.originVertex, journey.legs[0].departureStop);
            journey.departureTime = firstDepartureTime - initialWalkingTime;

            const int lastArrivalTime = journey.legs.back().arrivalTime;
            const int finalWalkingTime = getWalkingTime(journey.legs.back().arrivalStop, demandEntry.destinationVertex);
            journey.arrivalTime = lastArrivalTime + finalWalkingTime;
            journey.travelTime = journey.arrivalTime - journey.departureTime;

            for (size_t i = 1; i < journey.legs.size(); i++) {
                const int arrivalTime = journey.legs[i-1].arrivalTime;
                const int departureTime = journey.legs[i].departureTime;
                const int walkingTime = getWalkingTime(journey.legs[i-1].arrivalStop, journey.legs[i].departureStop);
                journey.transferWaitingTime += departureTime - arrivalTime - walkingTime;
            }
            journey.perceivedTravelTime =  journey.travelTime + (settings.transferCosts * journey.numberOfTransfers()) + (settings.waitingCosts * journey.transferWaitingTime);
        }

        return journey;
    }

    inline std::vector<JourneyInfo> mergeIdenticalJourneys(std::vector<JourneyInfo>& journeys) const noexcept {
        std::sort(journeys.begin(), journeys.end());
        std::vector<JourneyInfo> mergedJourneys;
        for (JourneyInfo& journey : journeys) {
            if (!mergedJourneys.empty() && Vector::equals(journey.legs, mergedJourneys.back().legs)) {
                mergedJourneys.back().groupSize += journey.groupSize;
            } else {
                mergedJourneys.emplace_back(journey);
            }
        }
        return mergedJourneys;
    }

private:
    const CSA::Data& data;
    const Settings& settings;
    const AccumulatedVertexDemand& demand;
    const AssignmentData& assignmentData;
};

}
