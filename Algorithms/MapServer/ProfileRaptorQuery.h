#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"

#include "../../Visualization/TimeTableVisualization.h"
#include "../../Visualization/FunctionVisualization.h"
#include "../../Visualization/PDF.h"

#include "../../Helpers/Meta.h"
#include "../../Helpers/Timer.h"

template<typename ALGORITHM, bool STOPS_ONLY = false>
class ProfileRaptorQuery : public Algorithm {

private:
    using Point = Geometry::Point;
    using Rectangle = Geometry::Rectangle;
    using Algo = ALGORITHM;
    static constexpr bool StopsOnly = STOPS_ONLY;
    using Type = ProfileRaptorQuery<Algo, StopsOnly>;
    using VertexType = Meta::IF<STOPS_ONLY, StopId, Vertex>;

public:
    ProfileRaptorQuery(Algo& algorithm, const RAPTOR::Data& data) :
        Algorithm(),
        algorithm(algorithm),
        data(data) {
        addParameter(Parameter("time", "minDepartureTime", "Min Departure time", "28800", true, 0));
        addParameter(Parameter("time", "maxDepartureTime", "Max Departure time", "32400", true, 0));
        addParameter(Parameter("bool", "travelTimeProfile", "Travel Time Profile", "false", true, 0));
    }

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetVertexId) noexcept {
        if constexpr (StopsOnly) {
            Assert(data.isStop(sourceVertexId), "The sourceVertexId " << sourceVertexId << " does not represent a stop!");
            Assert(data.isStop(targetVertexId), "The targetVertexId " << targetVertexId << " does not represent a stop!");
        }
        const int minDepartureTime = getParameter<int>("minDepartureTime");
        const int maxDepartureTime = getParameter<int>("maxDepartureTime");
        const bool travelTimeProfile = getParameter<bool>("travelTimeProfile");
        std::cout << "Run Profile-RAPTOR algorithm!" << std::endl;
        Timer timer;
        algorithm.run(VertexType(sourceVertexId), VertexType(targetVertexId), minDepartureTime, maxDepartureTime);
        const double time = timer.elapsedMilliseconds();
        std::cout << "Done (" << String::msToString(time) << ")!" << std::endl << std::endl;
        std::stringstream info;
        info << "Source: " << sourceVertexId;
        info << "</br>Target: " << targetVertexId;
        info << "</br>Computation Time: " << String::msToString(time);
        RAPTOR::ProfileHandle profile = algorithm.getProfileHandle();
        profile.sort();
        if (profile.empty()) {
            info << "</br></br>No Journeys found!";
        } else {
            const std::string pdfFileName = FileSystem::extendPath(FileSystem::getWorkingDirectory(), "profile_" + String::timeString() + ".pdf");
            info << "</br>Profile entries: " << profile.size();
            info << "</br>Profile : <a href='" << pdfFileName << "'>PDF</a>";
            const std::vector<std::vector<Point>> functions = (travelTimeProfile) ? (algorithm.getTravelTimeFunctions()) : (algorithm.getArrivalTimeFunctions());
            Rectangle bb1 = Rectangle::BoundingBox(Point(Construct::XY, minDepartureTime, std::min(0, minDepartureTime)), Point(Construct::XY, maxDepartureTime, 1000000));
            Rectangle bb2 = Rectangle::BoundingBox(functions);
            const double dt = findTimeStep(bb2.dy() / 16);
            bb2.discretize(dt, dt);
            FunctionVisualization<PDF> pdf = FunctionVisualization<PDF>::FromBoundingBox(pdfFileName, Rectangle::Intersection(bb1, bb2));
            pdf.drawGrid(dt, dt);
            for (size_t i = 0; i < functions.size(); i++) {
                pdf.drawFunction(functions[i], colors[i]);
            }
            pdf.clearOverflow();
            pdf.drawBox();
            pdf.drawXAxis(dt, [](const double x){return String::secToTime(x);});
            pdf.drawYAxis(dt, [](const double x){return String::secToTime(x);});
            std::cout << profile << std::endl;
        }
        std::vector<Result> results(1);
        results.back().color = KIT::green;
        results.back().nodes.push_back(Result::Node(sourceVertexId, "Source", KIT::blue));
        results.back().nodes.push_back(Result::Node(targetVertexId, "Target", KIT::blue));
        results.back().info = info.str();
        return results;
    }

private:
    inline double findTimeStep(const double dt) const noexcept {
        double result = 240 * 60;
        for (const double d : {180 * 60, 120 * 60, 60 * 60, 30 * 60, 20 * 60, 15 * 60, 10 * 60, 5 * 60}) {
            if (d > dt && d < result) result = d;
        }
        return result;
    }

private:
    Algo& algorithm;
    const RAPTOR::Data& data;

};

