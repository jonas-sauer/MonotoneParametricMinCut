#pragma once

#include <vector>

#include "Servlet.h"

#include "../Responses/ResponseJSON.h"

#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Geometry/Rectangle.h"

class BoundingBoxServlet : public MapServer::Servlet {

public:
    BoundingBoxServlet(const Geometry::Rectangle& boundingBox) :
        Servlet("bounding_box"),
        boundingBox(boundingBox) {
    }
    BoundingBoxServlet(const std::vector<Geometry::Point>& coordinates) :
        Servlet("bounding_box"),
        boundingBox(Geometry::Rectangle::BoundingBox(coordinates)) {
    }

    virtual ~BoundingBoxServlet() {}

    virtual void Process(mg_connection* conn, const mg_request_info*) {
        std::stringstream json;
        json << "{" << std::endl
             << "\"Min\": { \"Latitude\":" << boundingBox.min.latitude << ", \"Longitude\":" << boundingBox.min.longitude << "}, " << std::endl
             << "\"Max\": { \"Latitude\":" << boundingBox.max.latitude << ", \"Longitude\":" << boundingBox.max.longitude << "}" << std::endl
             << "}" << std::endl;

        MapServer::ResponseJSON response(json.str());
        response.Write(conn);
    }

private:
    Geometry::Rectangle boundingBox;

};
