#pragma once

#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Intermediate/Data.h"

#include "../../Algorithms/CH/CH.h"

#include "../../Visualization/MapVisualization.h"

#include "../../Shell/Shell.h"
using namespace Shell;

class LoadDimacsGraph : public ParameterizedCommand {

public:
    LoadDimacsGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "loadDimacsGraph", "Converts DIMACS graph data to our transfer graph format.") {
        addParameter("Input file");
        addParameter("Output file");
        addParameter("Graph type", "dynamic", { "static", "dynamic" });
        addParameter("Coordinate factor", "0.000001");
    }

    virtual void execute() noexcept {
        std::string graphType = getParameter("Graph type");
        if (graphType == "static") {
            load<TransferGraph>();
        } else {
            load<DynamicTransferGraph>();
        }
    }

private:
    template<typename GRAPH_TYPE>
    inline void load() const noexcept {
        DimacsGraphWithCoordinates dimacs;
        dimacs.fromDimacs<true>(getParameter("Input file"), getParameter<double>("Coordinate factor"));
        //Graph::applyBoundingBox(dimacs, BoundingBoxes::Europe);
        Graph::printInfo(dimacs);
        dimacs.printAnalysis();
        GRAPH_TYPE graph;
        Graph::move(std::move(dimacs), graph);
        Graph::printInfo(graph);
        graph.printAnalysis();
        graph.writeBinary(getParameter("Output file"));
    }
};

class DynamicToStaticTransferGraph : public ParameterizedCommand {

public:
    DynamicToStaticTransferGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "dynamicToStaticTransferGraph", "Converts a dynamic to a static transfer graph.") {
        addParameter("Input file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        DynamicTransferGraph dynamicGraph(getParameter("Input file"));
        TransferGraph staticGraph;
        Graph::move(std::move(dynamicGraph), staticGraph);
        staticGraph.writeBinary(getParameter("Output file"));
    }
};

class ReverseTransferGraph : public ParameterizedCommand {

public:
    ReverseTransferGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "reverseTransferGraph", "Reverses a static transfer graph.") {
        addParameter("Input file");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        TransferGraph graph(getParameter("Input file"));
        graph.revert();
        graph.writeBinary(getParameter("Output file"));
    }
};

class SaveDimacsGraph : public ParameterizedCommand {

public:
    SaveDimacsGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "saveDimacsGraph", "Saves a static transfer graph in DIMACS format.") {
        addParameter("Static graph");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        TransferGraph graph(getParameter("Static graph"));
        Graph::printInfo(graph);
        graph.printAnalysis();
        Graph::toDimacs(getParameter("Output file"), graph, graph[TravelTime]);
    }
};
