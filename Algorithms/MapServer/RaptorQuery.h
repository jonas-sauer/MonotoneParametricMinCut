#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/RAPTOR/Data.h"

#include "../../Visualization/TimeTableVisualization.h"
#include "../../Visualization/PDF.h"

#include "../../Helpers/Meta.h"

template<typename ALGORITHM, bool STOPS_ONLY = false>
class RaptorQuery : public Algorithm {

public:
    using Algo = ALGORITHM;
    static constexpr bool StopsOnly = STOPS_ONLY;
    using Type = RaptorQuery<Algo, StopsOnly>;
    using VertexType = Meta::IF<STOPS_ONLY, StopId, Vertex>;

    RaptorQuery(Algo& algorithm, const RAPTOR::Data& data) :
        Algorithm(),
        algorithm(algorithm),
        data(data) {
        addParameter(Parameter("time", "departureTime", "Departure time", "28800", true, 0));
    }

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetVertexId) noexcept {
        if constexpr (StopsOnly) {
            AssertMsg(data.isStop(sourceVertexId), "The sourceVertexId " << sourceVertexId << " does not represent a stop!");
            AssertMsg(data.isStop(targetVertexId), "The targetVertexId " << targetVertexId << " does not represent a stop!");
        }
        const VertexType sourceVertex(sourceVertexId);
        const VertexType targetVertex(targetVertexId);
        const int departureTime = getParameter<int>("departureTime");
        std::cout << "Run RAPTOR algorithm!" << std::endl;
        algorithm.run(sourceVertex, departureTime, targetVertex);
        std::cout << "Done!" << std::endl << std::endl;
        std::vector<Result> results(1);
        results.back().color = KIT::green;
        results.back().nodes.push_back(Result::Node(sourceVertexId, "Source", KIT::blue));
        results.back().nodes.push_back(Result::Node(targetVertexId, "Target", KIT::blue));
        results.back().info = "Source: " + sourceVertexId + "</br>Target: " + targetVertexId;
        TimeTableVisualization<PDF> pdf = TimeTableVisualization<PDF>::FromRAPTOR("RaptorQueryResult.pdf", data);
        if (algorithm.reachable(targetVertex)) {
            for (const RAPTOR::Journey& journey : algorithm.getJourneys()) {
                results.emplace_back();
                pdf.drawJourney(journey);
                pdf.newPage();
                std::vector<std::string> directions = data.journeyToText(journey);
                std::stringstream info;
                info << "Arrival time: " + String::secToTime(journey.back().arrivalTime, true);
                for (std::string& s : directions) {
                    std::cout << s << std::endl;
                    info << "</br></br>" << String::trim(String::replaceAll(s, '_', " "));
                }
                std::cout << std::endl;
                results.back().color = KIT::green;
                results.back().info = info.str();
                results.back().pathes.push_back(Result::Path(RAPTOR::journeyToPath(journey), "A Path", KIT::green));
            }
        } else {
            results.back().info += "</br>No Journeys found.";
        }
        return results;
    }

private:
    Algo& algorithm;
    const RAPTOR::Data& data;

};
