#pragma once

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <math.h>

#include "../../Shell/Shell.h"

#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"
#include "../../Helpers/String/TextFileUtils.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Vector/Permutation.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Debug.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Calendar.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/File.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/Vector/Vector.h"

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/GTFS/Entities/Calendar.h"
#include "../../DataStructures/GTFS/Data.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Demand/IdVertexDemand.h"
#include "../../DataStructures/Demand/PassengerData.h"
#include "../../DataStructures/BikeSharing/Stations.h"

#include "../../Visualization/Color.h"
#include "../../Visualization/PDF.h"
#include "../../Visualization/PNG.h"
#include "../../Visualization/SVG.h"
#include "../../Visualization/TimeTableVisualization.h"
#include "../../Visualization/CsaVisualization.h"
#include "../../Visualization/Borders.h"

using namespace Shell;

class DrawRAPTOR : public ParameterizedCommand {

public:
    DrawRAPTOR(BasicShell& shell) :
        ParameterizedCommand(shell, "drawRAPTOR", "Draws binary RAPTOR data in a graphics file.") {
        addParameter("RAPTOR binary");
        addParameter("Output file");
        addParameter("Format", { "pdf", "png", "svg"});
    }

    virtual void execute() noexcept {
        const std::string format = getParameter("Format");
        if (format == "pdf") {
            draw<PDF>();
        } else if (format == "png") {
            draw<PNG>();
        } else if (format == "svg") {
            draw<SVG>();
        }
    }

private:
    template<typename FORMAT>
    inline void draw() const noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        data.printInfo();
        draw(TimeTableVisualization<FORMAT>::FromRAPTOR(getParameter("Output file"), data, 0.3));
    }

    template<typename DOCUMENT>
    inline static void draw(DOCUMENT doc) noexcept {
        doc.drawEdges(Color::Grey);
        doc.drawRoutesByType();
        for (const int type : GTFS::Types) {
            doc.newPage();
            doc.drawRoutesByType(type);
        }
        doc.newPage();
        doc.drawEdges(Color::Grey);
        doc.newPage();
        doc.drawStops(Color::Black);
        doc.newPage();
        doc.drawRoutesByType();
        doc.drawStops(Color::Black, true);
        doc.newPage();
        doc.drawEdges(Color::LightGrey);
        doc.drawRoutes(Color::KITseablue, 20.0);
    }
};

class DrawCSA : public ParameterizedCommand {

public:
    DrawCSA(BasicShell& shell) :
        ParameterizedCommand(shell, "drawCSA", "Draws binary CSA data to a graphics file.") {
        addParameter("CSA binary");
        addParameter("Output file");
        addParameter("Format", {"pdf", "png", "svg"});
        addParameter("Draw edges", "true", {"true", "false"});
        addParameter("Draw connections", "true", {"true", "false"});
        addParameter("Draw stops", "true", {"true", "false"});
        addParameter("Draw borders", "false", {"true", "false"});
        addParameter("Draw connections once", "false", {"true", "false"});
        addParameter("Draw legend", "true", {"true", "false"});
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string outputFile = getParameter("Output file");
        const std::string format = getParameter("Format");
        const bool drawLegend = getParameter<bool>("Draw legend");

        CSA::Data data = CSA::Data::FromBinary(inputFile);
        data.printInfo();

        const double legendWidth = drawLegend ? 0.3 : 0.0;
        if (format == "pdf") draw(CsaVisualization<PDF>::FromCSA(outputFile, data, legendWidth));
        if (format == "png") draw(CsaVisualization<PNG>::FromCSA(outputFile, data, legendWidth));
        if (format == "svg") draw(CsaVisualization<SVG>::FromCSA(outputFile, data, legendWidth));
    }

private:
    template<typename DOCUMENT>
    inline void draw(DOCUMENT doc) const noexcept {
        const bool drawEdges = getParameter<bool>("Draw edges");
        const bool drawConnections = getParameter<bool>("Draw connections");
        const bool drawStops = getParameter<bool>("Draw stops");
        const bool drawBorders = getParameter<bool>("Draw borders");
        const bool drawConnectionsOnce = getParameter<bool>("Draw connections once");
        const bool drawLegend = getParameter<bool>("Draw legend");

        if (drawBorders) {
            for (const Border& border : Borders::GermanStates) {
                doc.drawLine(border, Color::Grey, 5, false);
            }
            doc.drawLine(Borders::Germany, Color::Black, 8, false);
        }
        if (drawEdges) doc.drawEdges(Color::Grey, drawLegend ? "Edges" : "");
        if (drawConnections) doc.drawConnections(Color::KITred, drawLegend ? "Connections" : "");
        if (drawConnectionsOnce) doc.drawConnectionsOnce(Color::KITblue, drawLegend ? "Connections" : "");
        if (drawStops) doc.drawStops(Color::Black, drawLegend ? "Stops" : "");
    }

};

class DrawAssignment : public ParameterizedCommand {

private:
    struct ConnectionLoad {
        ConnectionLoad(const ConnectionId connectionId = noConnection, const double load = 0) :
            connectionId(connectionId),
            load(load) {
        }
        ConnectionId connectionId;
        double load;
    };

public:
    DrawAssignment(BasicShell& shell) :
        ParameterizedCommand(shell, "drawAssignment", "Draws a CSA network according to loads computed by an assignment algorithm.") {
        addParameter("CSA binary");
        addParameter("Load file");
        addParameter("Output file");
        addParameter("Output format", "png", {"pdf", "png", "svg"});
        addParameter("Min time", "00:00:00");
        addParameter("Max time", "24:00:00");
        addParameter("Min latitude", "0");
        addParameter("Min longitude", "0");
        addParameter("Max latitude", "90");
        addParameter("Max longitude", "90");
        addParameter("Max travel time", "86400");
        addParameter("Max connection length in cm", "1000000");
    }

    virtual void execute() noexcept {
        const std::string csaFileName = getParameter("CSA binary");
        const std::string outputFileName = getParameter("Output file");
        const std::string outputFormat = getParameter("Output format");
        const double minLatitude = getParameter<double>("Min latitude");
        const double minLongitude = getParameter<double>("Min longitude");
        const double maxLatitude = getParameter<double>("Max latitude");
        const double maxLongitude = getParameter<double>("Max longitude");
        const Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(Geometry::Point(Construct::LatLong, minLatitude, minLongitude), Geometry::Point(Construct::LatLong, maxLatitude, maxLongitude));

        CSA::Data data = CSA::Data::FromBinary(csaFileName);
        data.printInfo();
        if (outputFormat == "pdf") draw(CsaVisualization<PDF>::FromBoundingBox(outputFileName, data, boundingBox), data);
        if (outputFormat == "png") draw(CsaVisualization<PNG>::FromBoundingBox(outputFileName, data, boundingBox), data);
        if (outputFormat == "svg") draw(CsaVisualization<SVG>::FromBoundingBox(outputFileName, data, boundingBox), data);
    }

private:
    template<typename DOCUMENT>
    inline void draw(DOCUMENT doc, const CSA::Data& data) const noexcept {
        const std::string loadFileName = getParameter("Load file");
        const int minTime = String::parseSeconds(getParameter("Min time"));
        const int maxTime = String::parseSeconds(getParameter("Max time"));
        const int maxTravelTime = String::parseSeconds(getParameter("Max travel time"));
        const int maxConnectionLength = String::parseSeconds(getParameter("Max connection length in cm"));

        std::vector<ConnectionLoad> connections;
        ConnectionLoad cl;
        double minLoad = 10000000.0;
        double maxLoad = 0;
        IO::CSVReader<2, IO::TrimChars<>, IO::DoubleQuoteEscape<',','"'>> in(loadFileName);
        in.readHeader("connectionId", "load");
        while (in.readRow(cl.connectionId, cl.load)) {
            if (cl.connectionId > data.connections.size()) continue;
            const CSA::Connection& connection = data.connections[cl.connectionId];
            if (connection.arrivalTime < minTime) continue;
            if (connection.departureTime > maxTime) continue;
            if (connection.travelTime() > maxTravelTime) continue;
            if (geoDistanceInCM(data.transferGraph.get(Coordinates, connection.departureStopId), data.transferGraph.get(Coordinates, connection.arrivalStopId)) > maxConnectionLength) continue;
            if ((!doc.contains(data.transferGraph.get(Coordinates, connection.departureStopId))) && (!doc.contains(data.transferGraph.get(Coordinates, connection.arrivalStopId)))) continue;
            connections.emplace_back(cl);
            if (minLoad > cl.load) minLoad = cl.load;
            if (maxLoad < cl.load) maxLoad = cl.load;
        }
        std::cout << "minLoad: " << minLoad << ", maxLoad: " << maxLoad << std::endl;

        drawEdges(doc, data);
        drawConnections(doc, connections, minLoad, maxLoad);
        doc.newPage();
        drawConnections(doc, connections, minLoad, maxLoad);
        doc.newPage();
        sort(connections, [](const ConnectionLoad& a, const ConnectionLoad& b){return a.load < b.load;});
        drawConnections(doc, connections, minLoad, maxLoad);
        doc.newPage();
        sort(connections, [](const ConnectionLoad& a, const ConnectionLoad& b){return a.load > b.load;});
        drawConnections(doc, connections, minLoad, maxLoad);
    }

    template<typename DOCUMENT>
    inline void drawEdges(DOCUMENT& doc, const CSA::Data& data) const noexcept {
        for (const Vertex from : data.transferGraph.vertices()) {
            std::vector<Edge> edges;
            for (const Edge edge : data.transferGraph.edgesFrom(from)) {
                edges.emplace_back(edge);
            }
            sort(edges, [&](const Edge a, const Edge b){
                return data.transferGraph.get(TravelTime, a) < data.transferGraph.get(TravelTime, b);
            });
            if (edges.size() > 4) {
                edges.resize(4);
            }
            for (const Edge edge : edges) {
                const Vertex to = data.transferGraph.get(ToVertex, edge);
                doc.drawLine(data.transferGraph.get(Coordinates, from), data.transferGraph.get(Coordinates, to), Color::Grey >> 0.5, 3, false);
            }
        }
    }

    template<typename DOCUMENT>
    inline void drawConnections(DOCUMENT& doc, const std::vector<ConnectionLoad>& connections, const double minLoad, const double maxLoad) const noexcept {
        const std::string loadFileName = getParameter("Load file");

        for (const ConnectionLoad& cl : connections) {
            const double relativeLoad = log(cl.load - minLoad + 1) / log(maxLoad - minLoad + 1);
            doc.drawConnection(cl.connectionId, Color::getGradientColor(Color::KITgreen, Color::KITorange, Color::KITred, relativeLoad), 5 + (40 * relativeLoad), false);
        }
    }

};
