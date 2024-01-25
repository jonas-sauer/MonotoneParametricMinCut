#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/RAPTOR/Data.h"

#include "../../Visualization/TimeTableVisualization.h"
#include "../../Visualization/FunctionVisualization.h"
#include "../../Visualization/PDF.h"

#include "../../Helpers/Timer.h"

template<typename ALGORITHM>
class MeatQuery : public Algorithm {

private:
    using Point = Geometry::Point;
    using Rectangle = Geometry::Rectangle;

public:
    MeatQuery(ALGORITHM& algorithm) :
        Algorithm(),
        algorithm(algorithm) {
        addParameter(Parameter("time", "minDepartureTime", "Min Departure time", "28800", true, 0));
        addParameter(Parameter("time", "maxDepartureTime", "Max Departure time", "32400", true, 1));
        addParameter(Parameter("time", "maxDelay", "Max Delay", "1800", true, 2));
        addParameter(Parameter("time", "transferCost", "Transfer Costs", "0", true, 3));
        addParameter(Parameter("bool", "travelTimeProfile", "Travel Time Profile", "false", true, 4));
    }

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetVertexId) noexcept {
        const int minDepartureTime = getParameter<int>("minDepartureTime");
        const int maxDepartureTime = getParameter<int>("maxDepartureTime");
        const int maxDelay = getParameter<int>("maxDelay");
        const int transferCost = getParameter<int>("transferCost");
        const bool travelTimeProfile = getParameter<bool>("travelTimeProfile");
        std::cout << "Run MEAT algorithm!" << std::endl;
        Timer timer;
        algorithm.run(targetVertexId, maxDelay, transferCost);
        const double time = timer.elapsedMilliseconds();
        std::cout << "Done (" << String::msToString(time) << ")!" << std::endl << std::endl;
        std::stringstream info;
        info << "Source: " << sourceVertexId;
        info << "</br>Target: " << targetVertexId;
        info << "</br>Computation Time: " << String::msToString(time);
        const std::vector<std::vector<Point>> functions = (travelTimeProfile) ? (algorithm.getTravelTimeFunctions(sourceVertexId)) : (algorithm.getArrivalTimeFunctions(sourceVertexId));
        if (functions[0].empty()) {
            info << "</br></br>No Journeys found!";
        } else {
            const std::string pdfFileName = FileSystem::extendPath(FileSystem::getWorkingDirectory(), "profile_" + String::timeString() + ".pdf");
            info << "</br>Profile entries: " << functions[0].size();
            info << "</br>Profile : <a href='" << pdfFileName << "'>PDF</a>";
            Rectangle bb1 = Rectangle::BoundingBox(Point(Construct::XY, minDepartureTime, std::min(0, minDepartureTime)), Point(Construct::XY, maxDepartureTime, 1000000));
            Rectangle bb2 = Rectangle::Empty();
            for (const std::vector<Point>& function : functions) {
                for (const Point& p : function) {
                    if (p.x >= minDepartureTime && p.x <= maxDepartureTime) bb2.extend(p);
                }
            }
            bb2.min.x = minDepartureTime;
            bb2.max.x = maxDepartureTime;
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
    ALGORITHM& algorithm;

};

