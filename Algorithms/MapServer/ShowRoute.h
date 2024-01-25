#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Geometry/Rectangle.h"
#include "../../DataStructures/RAPTOR/Data.h"

#include "../../Helpers/String/String.h"

class ShowRoute : public Algorithm {

public:
    ShowRoute(const RAPTOR::Data& data) : Algorithm(), data(data) {
        addParameter(Parameter("int", "routeId", "Route Id", "0", true));
    }

    virtual std::vector<Result> run(const Vertex, const Vertex) noexcept {
        const RouteId routeId = RouteId(getParameter<int>("routeId"));
        std::vector<Result> results;
        if (data.isRoute(routeId)) {
            data.printRoute(routeId);
            std::vector<Vertex> vertices;
            for (const StopId stop : data.stopsOfRoute(routeId)) {
                vertices.push_back(stop);
            }
            results.emplace_back(data.routeData[routeId].name);
            results.back().nodes.emplace_back(Result::Node(vertices.front(), "S", KIT::green));
            results.back().nodes.emplace_back(Result::Node(vertices.back(), "T", KIT::green));
            results.back().pathes.emplace_back(Result::Path(vertices, "Route", KIT::red));
        } else {
            results.emplace_back("Route not found!");
        }
        return results;
    }

private:
    const RAPTOR::Data& data;

};

class ShowRoutes : public Algorithm {

public:
    ShowRoutes(const RAPTOR::Data& data) : Algorithm(), data(data) {
    }

    virtual std::vector<Result> run(const Vertex source, const Vertex) noexcept {
        Assert(data.isStop(source), "The vertex " << source << " does not represent a stop!");
        std::vector<Result> results;
        for (const RAPTOR::RouteSegment route : data.routesContainingStop(StopId(source))) {
            data.printRoute(route.routeId);
            std::cout << "   maxRouteSpeed: " << data.maxRouteSpeed(route) << std::endl;
            std::cout << "   maxRouteDistance: " << data.maxRouteDistance(route) << std::endl;
            std::cout << "   maxRouteSpeedTimesDistance: " << data.maxRouteSpeedTimesDistance(route) << std::endl;
            std::vector<Vertex> vertices;
            for (const StopId stop : data.stopsOfRoute(route.routeId)) {
                vertices.push_back(stop);
            }
            results.emplace_back(data.routeData[route.routeId].name + " [" + std::to_string(route.routeId) + "]");
            results.back().nodes.emplace_back(Result::Node(vertices.front(), "S", KIT::green));
            results.back().nodes.emplace_back(Result::Node(vertices.back(), "T", KIT::green));
            results.back().pathes.emplace_back(Result::Path(vertices, "Route", KIT::red));
        }
        return results;
    }

private:
    const RAPTOR::Data& data;

};
