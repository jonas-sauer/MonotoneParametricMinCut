#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "MapVisualization.h"

#include "../DataStructures/CSA/Data.h"
#include "../DataStructures/Graph/Graph.h"
#include "../DataStructures/GTFS/Entities/Route.h"

template<typename DOCUMENT_TYPE>
class CsaVisualization : public MapVisualization<DOCUMENT_TYPE> {

public:
    using DocumentType = DOCUMENT_TYPE;
    using Type = CsaVisualization<DocumentType>;
    using Point = Geometry::Point;
    using Rectangle = Geometry::Rectangle;

private:
    using Super = MapVisualization<DocumentType>;

private:
    CsaVisualization(const std::string& fileName, const CSA::Data& data, const Rectangle& boundingBox) :
        Super(fileName, boundingBox),
        data(data) {
    }
    CsaVisualization(const std::string& fileName, const CSA::Data& data, const std::vector<Point>& coordinates, const double sideBarWidth) :
        Super(fileName, coordinates, sideBarWidth),
        data(data) {
    }

public:
    inline static Type FromBoundingBox(const std::string& fileName, const CSA::Data& data, Rectangle boundingBox) noexcept {
        return Type(fileName, data, boundingBox);
    }
    inline static Type FromCSA(const std::string& fileName, const CSA::Data& data, const double sideBarWidth = 0.0) noexcept {
        return Type(fileName, data, data.getCoordinates(), sideBarWidth);
    }

    inline void drawConnection(const size_t connectionId, const Color& color, const double stroke, const bool removeIfClipped = true) noexcept {
        const CSA::Connection& connection = data.connections[connectionId];
        Super::drawLine(data.transferGraph.get(Coordinates, connection.departureStopId), data.transferGraph.get(Coordinates, connection.arrivalStopId), color, stroke, removeIfClipped);
    }
    inline void drawConnection(const size_t connectionId, const Color& color, const bool removeIfClipped = true) noexcept {drawLine(connectionId, color, Super::defaultStroke, removeIfClipped);}
    inline void drawConnection(const size_t connectionId, const double stroke, const bool removeIfClipped = true) noexcept {drawLine(connectionId, Super::defaultColor, stroke, removeIfClipped);}
    inline void drawConnection(const size_t connectionId, const bool removeIfClipped = true) noexcept {drawLine(connectionId, Super::defaultColor, Super::defaultStroke, removeIfClipped);}

    inline void drawConnections(const Color& color, const double stroke, const std::string& caption = "Connections") noexcept {
        for (const CSA::Connection& connection : data.connections) {
            Super::drawLine(data.transferGraph.get(Coordinates, connection.departureStopId), data.transferGraph.get(Coordinates, connection.arrivalStopId), color, stroke);
        }
        if (caption != "") Super::write(caption + ": " + String::prettyInt(data.numberOfConnections()) + "\n", color);
    }
    inline void drawConnections(const Color& color, const std::string& caption = "Connections") noexcept {drawConnections(color, Super::defaultStroke, caption);}
    inline void drawConnections(const double stroke, const std::string& caption = "Connections") noexcept {drawConnections(Super::defaultColor, stroke, caption);}
    inline void drawConnections(const std::string& caption = "Connections") noexcept {drawConnections(Super::defaultColor, Super::defaultStroke, caption);}

    inline void drawConnectionsOnce(const Color& color, const double stroke, const std::string& caption, const int minTravelTime = 0) noexcept {
        SimpleDynamicGraph graph;
        graph.addVertices(data.numberOfStops());
        for (const CSA::Connection& connection : data.connections) {
            if (connection.travelTime() < minTravelTime) continue;
            if (connection.departureStopId < connection.arrivalStopId) {
                graph.findOrAddEdge(connection.departureStopId, connection.arrivalStopId);
            } else {
                graph.findOrAddEdge(connection.arrivalStopId, connection.departureStopId);
            }
        }
        size_t lineCount = 0;
        for (const Vertex from : graph.vertices()) {
            for (const Edge edge : graph.edgesFrom(from)) {
                Super::drawLine(data.transferGraph.get(Coordinates, from), data.transferGraph.get(Coordinates, graph.get(ToVertex, edge)), color, stroke);
                lineCount++;
            }
        }
        std::cout << "Used " << String::prettyInt(lineCount) << " lines to draw " << String::prettyInt(data.numberOfConnections()) << " connections!" << std::endl;
        if (caption != "") Super::write(caption + ": " + String::prettyInt(data.numberOfConnections()) + "\n", color);
    }
    inline void drawConnectionsOnce(const Color& color, const double stroke, const int minTravelTime = 0) noexcept {drawConnectionsOnce(color, stroke, "Connections", minTravelTime);}
    inline void drawConnectionsOnce(const Color& color, const std::string& caption, const int minTravelTime = 0) noexcept {drawConnectionsOnce(color, Super::defaultStroke, caption, minTravelTime);}
    inline void drawConnectionsOnce(const double stroke, const std::string& caption, const int minTravelTime = 0) noexcept {drawConnectionsOnce(Super::defaultColor, stroke, caption, minTravelTime);}
    inline void drawConnectionsOnce(const Color& color, const int minTravelTime = 0) noexcept {drawConnectionsOnce(color, Super::defaultStroke, "Connections", minTravelTime);}
    inline void drawConnectionsOnce(const double stroke, const int minTravelTime = 0) noexcept {drawConnectionsOnce(Super::defaultColor, stroke, "Connections", minTravelTime);}
    inline void drawConnectionsOnce(const std::string& caption, const int minTravelTime = 0) noexcept {drawConnectionsOnce(Super::defaultColor, Super::defaultStroke, caption, minTravelTime);}
    inline void drawConnectionsOnce(const int minTravelTime = 0) noexcept {drawConnectionsOnce(Super::defaultColor, Super::defaultStroke, "Connections", minTravelTime);}

    inline void drawEdges(const Color& color, const double stroke, const std::string& caption = "Transfers") noexcept {
        for (const Vertex fromVertex : data.transferGraph.vertices()) {
            for (const Edge edge : data.transferGraph.edgesFrom(fromVertex)) {
                const Vertex toVertex = data.transferGraph.get(ToVertex, edge);
                if (toVertex < fromVertex) continue;
                Super::drawLine(data.transferGraph.get(Coordinates, toVertex), data.transferGraph.get(Coordinates, fromVertex), color, stroke);
            }
        }
        if (caption != "") Super::write(caption + ": " + String::prettyInt(data.transferGraph.numEdges() / 2) + "\n", color);
    }
    inline void drawEdges(const Color& color, const std::string& caption = "Transfers") noexcept {drawEdges(color, Super::defaultStroke, caption);}
    inline void drawEdges(const double stroke, const std::string& caption = "Transfers") noexcept {drawEdges(Super::defaultColor, stroke, caption);}
    inline void drawEdges(const std::string& caption = "Transfers") noexcept {drawEdges(Super::defaultColor, Super::defaultStroke, caption);}

    inline void drawStops(const Color& color, const Icon icon, const double size, const std::string& caption = "Stops") noexcept {
        for (const StopId stop : data.stops()) {
            Super::drawPoint(data.stopData[stop].coordinates, color, icon, size);
        }
        if (caption != "") {
            Super::write(color, icon, size);
            Super::write(caption + ": " + String::prettyInt(data.numberOfStops()) + "\n", color);
        }
    }
    inline void drawStops(const Color& color, const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(color, icon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const Color& color, const double size, const std::string& caption = "Stops") noexcept {drawStops(color, Super::defaultIcon, size, caption);}
    inline void drawStops(const Color& color, const std::string& caption = "Stops") noexcept {drawStops(color, Super::defaultIcon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const Icon icon, const double size, const std::string& caption = "Stops") noexcept {drawStops(Super::defaultColor, icon, size, caption);}
    inline void drawStops(const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(Super::defaultColor, icon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const double size, const std::string& caption = "Stops") noexcept {drawStops(Super::defaultColor, Super::defaultIcon, size, caption);}
    inline void drawStops(const std::string& caption = "Stops") noexcept {drawStops(Super::defaultColor, Super::defaultIcon, Super::defaultStroke * 5, caption);}

    inline void drawStops(const std::vector<StopId>& stops, const Color& color, const Icon icon, const double size, const std::string& caption = "Stops") noexcept {
        for (const StopId stop : stops) {
            Super::drawPoint(data.stopData[stop].coordinates, color, icon, size);
        }
        if (caption != "") {
            Super::write(color, icon);
            Super::write(caption + " : " + String::prettyInt(stops.size()) + "\n", color);
        }
    }
    inline void drawStops(const std::vector<StopId>& stops, const Color& color, const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(stops, color, icon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const Color& color, const double size, const std::string& caption = "Stops") noexcept {drawStops(stops, color, Super::defaultIcon, size, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const Color& color, const std::string& caption = "Stops") noexcept {drawStops(stops, color, Super::defaultIcon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const Icon icon, const double size, const std::string& caption = "Stops") noexcept {drawStops(stops, Super::defaultColor, icon, size, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(stops, Super::defaultColor, icon, Super::defaultStroke * 5, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const double size, const std::string& caption = "Stops") noexcept {drawStops(stops, Super::defaultColor, Super::defaultIcon, size, caption);}
    inline void drawStops(const std::vector<StopId>& stops, const std::string& caption = "Stops") noexcept {drawStops(stops, Super::defaultColor, Super::defaultIcon, Super::defaultStroke * 5, caption);}

    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Color& color, const Icon icon, const double size, const std::string& caption = "Stops") noexcept {
        size_t stopCount = 0;
        for (const StopId stop : data.stops()) {
            if (draw(stop)) {
                Super::drawPoint(data.stopData[stop].coordinates, color, icon, size);
                stopCount++;
            }
        }
        if (caption != "") {
            Super::write(color, icon, size);
            Super::write(caption + ": " + String::prettyInt(stopCount) + "\n", color);
        }
    }
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Color& color, const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(draw, color, icon, Super::defaultStroke * 5, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Color& color, const double size, const std::string& caption = "Stops") noexcept {drawStops(draw, color, Super::defaultIcon, size, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Color& color, const std::string& caption = "Stops") noexcept {drawStops(draw, color, Super::defaultIcon, Super::defaultStroke * 5, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Icon icon, const double size, const std::string& caption = "Stops") noexcept {drawStops(draw, Super::defaultColor, icon, size, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const Icon icon, const std::string& caption = "Stops") noexcept {drawStops(draw, Super::defaultColor, icon, Super::defaultStroke * 5, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const double size, const std::string& caption = "Stops") noexcept {drawStops(draw, Super::defaultColor, Super::defaultIcon, size, caption);}
    template<typename DRAW>
    inline void drawStops(const DRAW& draw, const std::string& caption = "Stops") noexcept {drawStops(draw, Super::defaultColor, Super::defaultIcon, Super::defaultStroke * 5, caption);}

private:
    const CSA::Data& data;

};
