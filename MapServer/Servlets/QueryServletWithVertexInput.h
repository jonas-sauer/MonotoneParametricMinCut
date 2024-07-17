#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Servlet.h"

#include "../QueryStringParser.h"
#include "../Responses/ResponseJSON.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/Helpers.h"
#include "../../Helpers/JSONString.h"

#include "../../Algorithms/MapServer/Algorithm.h"

#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"

template<typename COORDINATE_TREE_TYPE>
class QueryServletWithVertexInput : public MapServer::Servlet {

public:

    QueryServletWithVertexInput(std::string handler, Algorithm& algorithm, const std::vector<Geometry::Point>& coordinates, const COORDINATE_TREE_TYPE& cs) :
        Servlet(handler),
        coordinates(coordinates),
        cs(cs),
        algorithm(algorithm),
        parameters(algorithm.getParameters()) {
        parameters["overwriteVertices"] = Parameter("bool", "overwriteVertices", "Use manual vertex IDs", "false", false, 1000);
        parameters["sourceId"] = Parameter("size_t", "sourceId", "Source Id", "-1", true, 1001);
        parameters["targetId"] = Parameter("size_t", "targetId", "Target Id", "-1", true, 1002);
    }

    virtual ~QueryServletWithVertexInput() {}

    virtual void Process(mg_connection* conn, const mg_request_info* requestInfo) {
        MapServer::QueryStringParser qsp(requestInfo);
        readParameters(qsp);
        const Vertex sourceNode = getSourceNodeId(qsp);
        const Vertex targetNode = getTargetNodeId(qsp);
        std::cout << "   Source = " << sourceNode << std::endl;
        std::cout << "   Target = " << targetNode  << std::endl;
        if (targetNode == noVertex) {
            sendResponse(conn, algorithm.runSourceOnly(sourceNode), sourceNode, targetNode);
        } else if (sourceNode == noVertex) {
            sendResponse(conn, algorithm.runTargetOnly(targetNode), sourceNode, targetNode);
        } else {
            sendResponse(conn, algorithm.run(sourceNode, targetNode), sourceNode, targetNode);
        }
        std::cout << "The calculation has been completed." << std::endl << std::endl << std::flush;
    }

    inline const std::map<std::string, Parameter>& getParameters() const {
        return parameters;
    }

protected:
    Vertex getSourceNodeId(MapServer::QueryStringParser& qsp) {
        if (parameters.find("overwriteVertices")->second.template getValue<bool>()) {
            return Vertex(parameters.find("sourceId")->second.template getValue<size_t>());
        } else {
            double lat = qsp.getValue<double>("SourceLat", 1000.0);
            double lon = qsp.getValue<double>("SourceLon", 1000.0);
            return getNodeId(Geometry::Point(Construct::LatLong, lat, lon));
        }
    }

    Vertex getTargetNodeId(MapServer::QueryStringParser& qsp) {
        if (parameters.find("overwriteVertices")->second.template getValue<bool>()) {
            return Vertex(parameters.find("targetId")->second.template getValue<size_t>());
        } else {
            double lat = qsp.getValue<double>("TargetLat", 1000.0);
            double lon = qsp.getValue<double>("TargetLon", 1000.0);
            return getNodeId(Geometry::Point(Construct::LatLong, lat, lon));
        }
    }

    virtual Vertex getNodeId(Geometry::Point pos) {
        if (pos.latitude == 1000.0 || pos.longitude == 1000.0) return noVertex;
        return cs.getNearestNeighbor(pos);
    }

    virtual void readParameters(MapServer::QueryStringParser& qsp) {
        for (std::map<std::string, Parameter>::iterator iterator = algorithm.getParameters().begin(); iterator != algorithm.getParameters().end(); iterator++) {
            iterator->second.value = qsp.getValue<std::string>(iterator->second.name, iterator->second.defaultValue);
        }
        for (std::map<std::string, Parameter>::iterator iterator = parameters.begin(); iterator != parameters.end(); iterator++) {
            iterator->second.value = qsp.getValue<std::string>(iterator->second.name, iterator->second.defaultValue);
            std::cout << "   Parameter(" << iterator->second.name << ") = " << iterator->second.value << std::endl;
        }
    }

    virtual void sendResponse(mg_connection* conn, std::vector<Result>&& response, Vertex sourceId = noVertex, Vertex targetId = noVertex) {
        if (!response.empty()) {
            for (std::map<std::string, Parameter>::iterator iterator = algorithm.getParameters().begin(); iterator != algorithm.getParameters().end(); iterator++) {
                response.front().parameter.push_back(std::make_pair(iterator->second.name, iterator->second.value));
            }
            response.front().parameter.push_back(std::make_pair("sourceId", std::to_string(sourceId)));
            response.front().parameter.push_back(std::make_pair("targetId", std::to_string(targetId)));
        }
        MapServer::ResponseJSON jsonresponse(resultsToJSON(response, coordinates));
        jsonresponse.Write(conn);
    }

private:
    const std::vector<Geometry::Point>& coordinates;
    const COORDINATE_TREE_TYPE& cs;
    Algorithm& algorithm;

    std::map<std::string, Parameter> parameters;

};
