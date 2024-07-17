#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Geometry/Rectangle.h"
#include "../../DataStructures/Graph/Graph.h"

#include "../../Helpers/String/String.h"

template<typename GRAPH>
class ShowEdges : public Algorithm {

public:
    ShowEdges(const GRAPH& graph) : Algorithm(), graph(graph), coordinates(graph[Coordinates]), bendEdges(true) {
        addParameter(Parameter("int", "travelTime", "Travel Time", "0", true));
        addParameter(Parameter("bool", "bendEdges", "Bend Edges", "true", true));
    }

    virtual std::vector<Result> run(const Vertex, const Vertex) noexcept {
        const int travelTime = getParameter<int>("travelTime");
        bendEdges = getParameter<bool>("bendEdges");
        const std::string label = "Edges";
        Result result(label);
        result.color = KIT::green;
        size_t count = 0;
        for (Vertex from : graph.vertices()) {
            for (Edge edge : graph.edgesFrom(from)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (graph.get(TravelTime, edge) != travelTime) continue;
                result.polyLines.push_back(getPolyLine(from, to));
                count++;
            }
        }
        result.info = "Found " + std::to_string(count) +  " edges.";
        return std::vector<Result>({result});
    }

    inline Result::PolyLine getPolyLine(const Vertex from, const Vertex to) noexcept {
        if (bendEdges) {
            Geometry::Point d = coordinates[to] - coordinates[from];
            Geometry::Point o = Geometry::Point(Construct::XY, d.y, -d.x);
            Geometry::Point m = coordinates[from] + (0.5 * d) + (0.02 * o);
            return Result::PolyLine(std::vector<Geometry::Point>{coordinates[from], m, coordinates[to]}, "", KIT::red);
        } else {
            return Result::PolyLine(std::vector<Geometry::Point>{coordinates[from], coordinates[to]}, "", KIT::red);
        }
    }

private:
    const GRAPH& graph;
    const std::vector<Geometry::Point>& coordinates;
    bool bendEdges;

};
