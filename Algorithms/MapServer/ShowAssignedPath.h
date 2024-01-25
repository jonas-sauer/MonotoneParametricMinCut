#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/PassengerData.h"

class ShowAssignedPath : public Algorithm {

public:
    ShowAssignedPath(const CSA::Data& csaData, const PassengerData& assignmentData) : Algorithm(), csaData(csaData), assignmentData(assignmentData) {
        addParameter(Parameter("int", "minIndex", "Min. Passenger Index", "0", true));
        addParameter(Parameter("int", "numResults", "Number of Results", "1", true));
    }

    virtual std::vector<Result> run(const Vertex, const Vertex) noexcept {
        const int minIndex = getParameter<int>("minIndex");
        const int numResults = getParameter<int>("numResults");
        const int maxIndex = std::min<int>(minIndex + numResults, assignmentData.entries.size());
        std::vector<Result> results;
        if (minIndex < 0) return results;
        for (int index = minIndex; index < maxIndex; index++) {
            results.emplace_back();
            results.back().color = KIT::green;
            results.back().nodes.push_back(Result::Node(Vertex(assignmentData.entries[index].originVertex), "O", KIT::red));
            results.back().nodes.push_back(Result::Node(Vertex(assignmentData.entries[index].destinationVertex), "D", KIT::red));
            results.back().info = "Index: " + String::prettyInt(index) + "</br>"
                                + "Id: " + String::prettyInt(getDestinationSpecificPassengerId(assignmentData.entries[index].id)) + "</br></br>"
                                + "Origin: " + assignmentData.entries[index].originVertex + "</br>"
                                + "Destination: " + assignmentData.entries[index].destinationVertex + "</br>";
            if (assignmentData.entries[index].numberOfConnections == 0) {
                results.back().info += "</br> No Connections are used";
            } else {
                results.back().nodes.push_back(Result::Node(assignmentData.entries[index].firstStop, "F", KIT::blue));
                results.back().nodes.push_back(Result::Node(assignmentData.entries[index].lastStop, "L", KIT::blue));
                results.back().pathes.push_back(Result::Path(assignmentData.getPath(csaData, index), "A Path", KIT::green));
                results.back().info += "First stop: " + assignmentData.entries[index].firstStop + "</br>"
                                     + "Last stop: " + assignmentData.entries[index].lastStop + "</br></br>"
                                     + "Travel time: " + String::secToString(assignmentData.entries[index].travelTimeWithoutInitialWaiting) + "</br>"
                                     + "#Connections: " + String::prettyInt(assignmentData.entries[index].numberOfConnections) + "</br>"
                                     + "#Trips: " + String::prettyInt(assignmentData.entries[index].numberOfTrips);
            }
        }
        return results;
    }

private:
    const CSA::Data& csaData;
    const PassengerData& assignmentData;

};
