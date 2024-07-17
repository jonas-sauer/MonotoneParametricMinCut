#pragma once

#include <vector>
#include <string>

#include "Geometry/Point.h"

#include "../Helpers/Types.h"
#include "../Helpers/JSONString.h"
#include "../Helpers/Vector/Vector.h"

namespace KIT {
    const std::string green  = "009682";
    const std::string blue   = "4664AA";
    const std::string orange = "DF9B1B";
    const std::string red    = "A22223";
    const std::string violet = "A3107C";
    const std::string cyan   = "23A1E0";
}

class Result {

public:
    struct PolyLine {
        PolyLine(const std::vector<Geometry::Point>& points = std::vector<Geometry::Point>(), const std::string& info = "", const std::string& color = "FF0000") :
            info(info),
            color(color),
            points(points) {
        }
        std::string info;
        std::string color;
        std::vector<Geometry::Point> points;
    };

    struct Path {
        Path(const std::vector<Vertex>& nodeIDs = std::vector<Vertex>(), const std::string& info = "", const std::string& color = "FF0000") :
            info(info),
            color(color),
            nodeIDs(nodeIDs) {
        }
        Path(const std::vector<StopId>& stopIDs, const std::string& info = "", const std::string& color = "FF0000") :
            info(info),
            color(color) {
            Vector::assign(nodeIDs, stopIDs);
        }
        std::string info;
        std::string color;
        std::vector<Vertex> nodeIDs;
    };

    struct Node {
        Node(const Vertex id, const std::string& info = "", const std::string& color = KIT::green) :
            info(info),
            color(color),
            id(id) {
        }
        std::string info;
        std::string color;
        Vertex id;
    };

    Result(const std::string& info = "", const std::string& color = "FF0000") :
        info(info),
        color(color) {
    }

    static Result fromPath(const std::vector<Vertex>& path) noexcept {
        Result result;
        result.pathes.push_back(Path(path));
        return result;
    }

    static Result fromNodeSet(const std::vector<Vertex>& nodeSet) noexcept {
        Result result;
        for (const Vertex v : nodeSet) {
            result.nodes.push_back(Node(v));
        }
        return result;
    }

    std::string info;
    std::string color;
    std::vector<Path> pathes;
    std::vector<Node> nodes;
    std::vector<PolyLine> polyLines;
    std::vector<Geometry::Point> function;
    std::vector<std::pair<std::string, std::string>> parameter;

};

inline std::string resultsToJSON(const std::vector<Result>& results, const std::vector<Geometry::Point>& coordinates) {
    JSONString json;
    json.startNewObject();
        json.addQualifier("ResultSet");
        json.startNewArray();
            for (size_t i = 0; i < results.size(); i++) {
                const Result& result = results[i];
                json.startNewObject();
                    json.addQualifierWithValue("Info", result.info);
                    json.addQualifierWithValue("Color", result.color);
                    json.addQualifier("Pathes");
                    json.startNewArray();
                        for (size_t j = 0; j < result.pathes.size(); j++) {
                            const Result::Path& path = result.pathes[j];
                            json.startNewObject();
                                json.addQualifierWithValue("Info", path.info);
                                json.addQualifierWithValue("Color", path.color);
                                json.addQualifier("Nodes");
                                json.startNewArray();
                                    for (size_t k = 0; k < path.nodeIDs.size(); k++) {
                                        const Geometry::Point& pos = coordinates[path.nodeIDs[k]];
                                        json.startNewObject();
                                            json.addQualifierWithValue("Lat", pos.latitude);
                                            json.addQualifierWithValue("Lon", pos.longitude);
                                        json.endObject();
                                    }
                                json.endArray();
                            json.endObject();
                        }
                        for (size_t j = 0; j < result.polyLines.size(); j++) {
                            const Result::PolyLine& polyLine = result.polyLines[j];
                            json.startNewObject();
                                json.addQualifierWithValue("Info", polyLine.info);
                                json.addQualifierWithValue("Color", polyLine.color);
                                json.addQualifier("Nodes");
                                json.startNewArray();
                                    for (size_t k = 0; k < polyLine.points.size(); k++) {
                                        const Geometry::Point& pos = polyLine.points[k];
                                        json.startNewObject();
                                            json.addQualifierWithValue("Lat", pos.latitude);
                                            json.addQualifierWithValue("Lon", pos.longitude);
                                        json.endObject();
                                    }
                                json.endArray();
                            json.endObject();
                        }
                    json.endArray();
                    json.addQualifier("Nodes");
                    json.startNewArray();
                        for (size_t j = 0; j < result.nodes.size(); j++) {
                            const Result::Node& node = result.nodes[j];
                            json.startNewObject();
                                json.addQualifierWithValue("Info", node.info);
                                json.addQualifierWithValue("Color", node.color);
                                const Geometry::Point& pos = coordinates[node.id];
                                json.addQualifier("Pos");
                                json.startNewObject();
                                    json.addQualifierWithValue("Lat", pos.latitude);
                                    json.addQualifierWithValue("Lon", pos.longitude);
                                json.endObject();
                            json.endObject();
                        }
                    json.endArray();
                    json.addQualifier("Function");
                    json.startNewArray();
                        for (const Geometry::Point& point : result.function) {
                            json.startNewObject();
                                json.addQualifierWithValue("X", point.x);
                                json.addQualifierWithValue("Y", point.y);
                            json.endObject();
                        }
                    json.endArray();
                    json.addQualifier("Parameter");
                    json.startNewObject();
                        for (const std::pair<std::string, std::string>& parameter : result.parameter) {
                            json.addQualifierWithValue(parameter.first, parameter.second);
                        }
                    json.endObject();
                json.endObject();
            }
        json.endArray();
    json.endObject();
    return json.str();
}
