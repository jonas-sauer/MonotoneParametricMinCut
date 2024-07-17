#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../../Algorithms/UnionFind.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/TripBased/Data.h"

#include "../../Helpers/HighlightText.h"

#include "../../Shell/Shell.h"

using namespace Shell;

class AnalyzeNetwork : public ParameterizedCommand {

public:
    AnalyzeNetwork(BasicShell& shell) :
        ParameterizedCommand(shell, "analyzeNetwork", "Analyzes and prints statistics for a network.") {
        addParameter("Network binary");
        addParameter("Network type", { "intermediate", "raptor", "csa" });
    }

    virtual void execute() noexcept {
        const std::string networkFile = getParameter("Network binary");
        const std::string networkType = getParameter("Network type");

        if (networkType == "intermediate") {
            analyze<Intermediate::Data>(networkFile);
        } else if (networkType == "raptor") {
            analyze<RAPTOR::Data>(networkFile);
        } else {
            analyze<CSA::Data>(networkFile);
        }
    }

private:
    template<typename DATA_TYPE>
    inline static void analyze(const std::string& networkFile) noexcept {
        DATA_TYPE network = DATA_TYPE::FromBinary(networkFile);
        network.printInfo();
    }
};

class AnalyzeConnectivity : public ParameterizedCommand {

private:
    struct EdgeEntry {
        Vertex from;
        Vertex to;
        int length;
        inline bool operator<(const EdgeEntry& other) const noexcept {return length < other.length;}
    };

public:
    AnalyzeConnectivity(BasicShell& shell) :
        ParameterizedCommand(shell, "analyzeConnectivity", "Analyzes the connectivity of the transfer graph.") {
        addParameter("Intermediate binary");
    }

    virtual void execute() noexcept {
        const std::string intermediateFile = getParameter("Intermediate binary");

        Intermediate::Data inter = Intermediate::Data::FromBinary(intermediateFile);
        inter.printInfo();
        std::vector<int64_t> sizeOfComponent(inter.numberOfStops(), 1);
        sizeOfComponent.resize(inter.transferGraph.numVertices(), 0);
        int64_t numberOfComponents = inter.numberOfStops();
        int64_t numberOfIsolatedStops = inter.numberOfStops();
        int64_t maxComponentSize = 0;
        int64_t numberOfEdges = 0;
        UnionFind components(inter.transferGraph.numVertices());
        std::vector<EdgeEntry> edges;
        for (const Vertex v : inter.transferGraph.vertices()) {
            for (const Edge e : inter.transferGraph.edgesFrom(v)) {
                edges.emplace_back(EdgeEntry{v, inter.transferGraph.get(ToVertex, e), inter.transferGraph.get(TravelTime, e)});
            }
        }
        std::sort(edges.begin(), edges.end());
        std::cout << std::setw(10) << "Index" << std::setw(14) << "Edge length" << std::setw(14) << "#Components" << std::setw(14) << "#Isolated" << std::setw(14) << "max size" << std::setw(14) << "#edges" << std::endl;
        int64_t lastLength = -1;
        for (const size_t i : indices(edges)) {
            if (components.find(edges[i].from) == components.find(edges[i].to)) continue;
            const int64_t sizeA = sizeOfComponent[components.find(edges[i].from)];
            const int64_t sizeB = sizeOfComponent[components.find(edges[i].to)];
            const int64_t newSize = sizeA + sizeB;
            const int64_t newEdges = sizeA * sizeB;
            components.unite(edges[i].from, edges[i].to);
            sizeOfComponent[components.find(edges[i].from)] = newSize;
            if (newEdges == 0) continue;
            numberOfComponents--;
            if (sizeA == 1) numberOfIsolatedStops--;
            if (sizeB == 1) numberOfIsolatedStops--;
            if (maxComponentSize < newSize) maxComponentSize = newSize;
            numberOfEdges += newEdges;
            if (numberOfEdges < 0) break;
            if (edges[i].length == lastLength) continue;
            lastLength = edges[i].length;
            std::cout << std::setw(10) << String::prettyInt(i) << std::setw(14) << String::prettyInt(edges[i].length) << std::setw(14) << String::prettyInt(numberOfComponents);
            std::cout << std::setw(14) << String::prettyInt(numberOfIsolatedStops) << std::setw(14) << String::prettyInt(maxComponentSize) << std::setw(14) << String::prettyInt(numberOfEdges) << std::endl;
        }
    }

};

class AnalyzeConnectedComponents : public ParameterizedCommand {

public:
    AnalyzeConnectedComponents(BasicShell& shell) :
        ParameterizedCommand(shell, "analyzeConnectedComponents", "Analyzes the connected components of a RAPTOR network.") {
        addParameter("RAPTOR binary");
    }

    virtual void execute() noexcept {
        const std::string raptorFile = getParameter("RAPTOR binary");

        RAPTOR::Data data = RAPTOR::Data::FromBinary(raptorFile);
        data.printInfo();
        RAPTOR::TransferGraph graph = data.minTravelTimeGraph();
        Graph::printInfo(graph);
        graph.printAnalysis();
        std::vector<int> componentSize;
        std::vector<int> component(graph.numVertices(), -1);
        for (const Vertex vertex : graph.vertices()) {
            if (component[vertex] != -1) continue;
            component[vertex] = componentSize.size();
            std::vector<Vertex> queue(1, vertex);
            componentSize.emplace_back(1);
            while (!queue.empty()) {
                const Vertex from = queue.back();
                queue.pop_back();
                for (const Edge edge : graph.edgesFrom(from)) {
                    const Vertex to = graph.get(ToVertex, edge);
                    if (component[to] != -1) continue;
                    component[to] = componentSize.size() - 1;
                    queue.emplace_back(to);
                    componentSize.back()++;
                }
            }
        }
        std::cout << "Connected components:" << std::endl;
        std::cout << std::setw(10) << "Index" << std::setw(14) << "Size" << std::endl;
        for (size_t i = 0; i < componentSize.size(); i++) {
            std::cout << std::setw(10) << String::prettyInt(i) << std::setw(14) << String::prettyInt(componentSize[i]) << std::endl;
        }

    }

};

class AnalyzeMinTransferTimes : public ParameterizedCommand {

public:
    AnalyzeMinTransferTimes(BasicShell& shell) :
        ParameterizedCommand(shell, "analyzeMinTransferTimes", "Analyzes minimum transfer times for a RAPTOR network.") {
        addParameter("RAPTOR binary");
    }

    virtual void execute() noexcept {
        const RAPTOR::Data data(getParameter("RAPTOR binary"));
        std::vector<size_t> histogram;
        for (const StopId stop : data.stops()) {
            const int minTransferTime = data.minTransferTime(stop);
            if (histogram.size() <= static_cast<size_t>(minTransferTime)) {
                histogram.resize(minTransferTime + 1, 0);
            }
            histogram[minTransferTime]++;
        }
        for (size_t i = 0; i < histogram.size(); i++) {
            if (histogram[i] == 0) continue;
            std::cout << "#Stops with min transfer time " << i << ": " << histogram[i] << std::endl;
        }
    }
};

class PrintRoute : public ParameterizedCommand {

public:
    PrintRoute(BasicShell& shell) :
        ParameterizedCommand(shell, "printRoute", "Prints all stop events of the given route.") {
        addParameter("RAPTOR binary");
        addParameter("Route ID");
    }

    virtual void execute() noexcept {
        RAPTOR::Data data = RAPTOR::Data::FromBinary(getParameter("RAPTOR binary"));
        data.printInfo();
        data.printRoute(getParameter<RouteId>("Route ID"));
    }
};

class PrintTrip : public ParameterizedCommand {

public:
    PrintTrip(BasicShell& shell) :
        ParameterizedCommand(shell, "printTrip", "Prints all stop events of the given trip.") {
        addParameter("Intermediate binary");
        addParameter("Trip ID");
    }

    virtual void execute() noexcept {
        Intermediate::Data data = Intermediate::Data::FromBinary(getParameter("Intermediate binary"));
        data.printInfo();
        data.printTrip(getParameter<TripId>("Trip ID"));
    }
};

class PrintStopEventID : public ParameterizedCommand {

public:
    PrintStopEventID(BasicShell& shell) :
        ParameterizedCommand(shell, "printStopEventID", "Finds the specified stop event and prints its ID.") {
        addParameter("Trip-Based file");
        addParameter("Route");
        addParameter("Stop");
        addParameter("Time");
    }

    virtual void execute() noexcept {
        const std::string tripBasedFile = getParameter("Trip-Based file");
        const RouteId route = getParameter<RouteId>("Route");
        const StopId stop = getParameter<StopId>("Stop");
        const int time = getParameter<int>("Time");

        TripBased::Data tripBasedData(tripBasedFile);
        tripBasedData.printInfo();

        tripBasedData.raptorData.printRoute(route);

        for (size_t i = tripBasedData.raptorData.firstStopEventOfRoute[route]; i < tripBasedData.raptorData.firstStopEventOfRoute[route + 1]; i++) {
            if (tripBasedData.raptorData.stopsOfRoute(route)[tripBasedData.indexOfStopEvent[i]] == stop) {
                if (tripBasedData.raptorData.stopEvents[i].arrivalTime == time) {
                    warning("Arrival even ", i);
                }
                if (tripBasedData.raptorData.stopEvents[i].departureTime == time) {
                    warning("Departure even ", i);
                }
            }
        }

    }

};

class PrintOutgoingTransfers : public ParameterizedCommand {

public:
    PrintOutgoingTransfers(BasicShell& shell) :
        ParameterizedCommand(shell, "printOutgoingTransfers", "Prints all transfers outgoing from the specified stop event.") {
        addParameter("Trip-Based file");
        addParameter("Stop event ID");
    }

    virtual void execute() noexcept {
        const std::string tripBasedFile = getParameter("Trip-Based file");
        const int stopEvent = getParameter<int>("Stop event ID");

        TripBased::Data tripBasedData(tripBasedFile);
        tripBasedData.printInfo();

        for (const Edge edge : tripBasedData.stopEventGraph.edgesFrom(Vertex(stopEvent))) {
            std::cout << tripBasedData.stopEventGraph.get(ToVertex, edge) << std::endl;
        }

    }

};
