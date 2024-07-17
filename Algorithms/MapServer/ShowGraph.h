#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Geometry/Rectangle.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/String/String.h"

template<typename GRAPH>
class ShowGraph : public Algorithm {

public:
    ShowGraph(const GRAPH& graph) : Algorithm(), graph(graph), coordinates(graph[Coordinates]), bendEdges(true) {
        addParameter(Parameter("bool", "bendEdges", "Bend Edges", "true", true));
        addParameter(Parameter("bool", "showVertices", "Show Vertices", "false", true));
    }

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetVertexId) noexcept {
        bendEdges = getParameter<bool>("bendEdges");
        showVertices = getParameter<bool>("showVertices");
        const std::string label = "Corner vertices:</br>" + sourceVertexId + "</br>" + targetVertexId;
        Result edges(label);
        Result vertices("Vertices");
        edges.color = KIT::green;
        vertices.color = KIT::orange;
        edges.nodes.push_back(Result::Node(sourceVertexId, "S", KIT::green));
        edges.nodes.push_back(Result::Node(targetVertexId, "T", KIT::green));
        Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(coordinates[sourceVertexId], coordinates[targetVertexId]);
        for (Vertex from : graph.vertices()) {
            if (!boundingBox.contains(coordinates[from])) continue;
            for (Edge edge : graph.edgesFrom(from)) {
                Vertex to = graph.get(ToVertex, edge);
                if (!boundingBox.contains(coordinates[to])) continue;
                edges.polyLines.push_back(getPolyLine(from, to));
            }
            if (showVertices) {
                vertices.nodes.push_back(Result::Node(from, "V", KIT::orange));
            }
        }
        if (showVertices) {
            return std::vector<Result>({edges, vertices});
        } else {
            return std::vector<Result>({edges});
        }
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
    bool showVertices;

};
