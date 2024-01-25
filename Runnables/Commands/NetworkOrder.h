#pragma once

#include <string>
#include <vector>

#include "../../Shell/Shell.h"

#include "../../Algorithms/CH/CH.h"
#include "../../Algorithms/BreadthFirstSearch.h"
#include "../../Algorithms/CH/Query/PHAST.h"
#include "../../Algorithms/Dijkstra/SimpleDijkstra.h"
#include "../../Algorithms/GraphHilbertCurveGenerator.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"

#include "../../Visualization/MapVisualization.h"
#include "../../Visualization/PDF.h"
#include "../../Visualization/PNG.h"
#include "../../Visualization/SVG.h"

using namespace Shell;

class MeasureDijkstra : public ParameterizedCommand {

public:
    MeasureDijkstra(BasicShell& shell) :
        ParameterizedCommand(shell, "measureDijkstra", "Measures the performance of Dijkstra's algorithm on the given graph.") {
        addParameter("Graph");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Graph"));

        SimpleDijkstra<RAPTOR::TransferGraph> dijkstra(graph);
        Timer timer;
        double totalTime = 0;
        const size_t n = getParameter<size_t>("Number of queries");
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % graph.numVertices());
            timer.restart();
            dijkstra.run(source, TravelTime);
            totalTime += timer.elapsedMilliseconds();
        }
        totalTime /= n;
        std::cout << "Dijkstra took " << String::msToString(totalTime) << std::endl;
    }
};

class MeasureBFS : public ParameterizedCommand {

public:
    MeasureBFS(BasicShell& shell) :
        ParameterizedCommand(shell, "measureBFS", "Measures the performance of breadth-first search on the given graph.") {
        addParameter("Graph");
        addParameter("Number of queries");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Graph"));
        srand(42);

        BreadthFirstSearch<RAPTOR::TransferGraph> bfs(graph);
        Timer timer;
        double totalTime = 0;
        const size_t n = getParameter<size_t>("Number of queries");
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % graph.numVertices());
            timer.restart();
            bfs.run(source);
            totalTime += timer.elapsedMilliseconds();
        }
        totalTime /= n;
        std::cout << "BFS took " << String::msToString(totalTime) << std::endl;
    }
};

class MeasurePHAST : public ParameterizedCommand {

public:
    MeasurePHAST(BasicShell& shell) :
        ParameterizedCommand(shell, "measurePHAST", "Measures the performance of PHAST on the given graph.") {
        addParameter("Graph");
        addParameter("CH data");
        addParameter("Number of queries");
        addParameter("Order");
        addParameter("Reorder vertices?");
        addParameter("Path retrieval?");
    }

    virtual void execute() noexcept {
        if (getParameter<bool>("Reorder vertices?")) {
            choosePathRetrieval<true>();
        } else {
            choosePathRetrieval<false>();
        }
    }

    template<bool REORDER_VERTICES>
    inline void choosePathRetrieval() const noexcept {
        if (getParameter<bool>("Path retrieval?")) {
            run<REORDER_VERTICES, true>();
        } else {
            run<REORDER_VERTICES, false>();
        }
    }

    template<bool REORDER_VERTICES, bool PATH_RETRIEVAL>
    inline void run() const noexcept {
        RAPTOR::TransferGraph graph(getParameter("Graph"));
        CH::CH ch(getParameter("CH data"));

        const Order order = getOrder(ch);
        if constexpr (REORDER_VERTICES) graph.applyVertexOrder(order);
        SimpleDijkstra<RAPTOR::TransferGraph> dijkstra(graph);
        CH::PHAST<REORDER_VERTICES, true, PATH_RETRIEVAL> phast(ch, order);

        Timer timer;
        double phastTime = 0;
        double dijkstraTime = 0;
        const size_t n = getParameter<size_t>("Number of queries");
        for (size_t i = 0; i < n; i++) {
            const Vertex source(rand() % ch.numVertices());
            timer.restart();
            phast.run(source);
            phastTime += timer.elapsedMilliseconds();
            timer.restart();
            dijkstra.run(source, TravelTime);
            dijkstraTime += timer.elapsedMilliseconds();
            for (const Vertex vertex : graph.vertices()) {
                const int phastDistance = phast.getDistance(vertex);
                const int dijkstraDistance = dijkstra.getDistance(vertex);
                if (phastDistance != dijkstraDistance) {
                    std::cout << "Vertex " << vertex << std::endl;
                    std::cout << "PHAST distance: " << phastDistance << std::endl;
                    std::cout << "Dijkstra distance: " << dijkstraDistance << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        phastTime /= n;
        dijkstraTime /= n;
        std::cout << "PHAST took " << String::msToString(phastTime) << std::endl;
        std::cout << "Dijkstra took " << String::msToString(dijkstraTime) << std::endl;
    }

private:
    inline Order getOrder(const CH::CH& ch) const noexcept {
        const std::string orderType = getParameter("Order");
        if (orderType == "Contraction") {
            return Order(Vector::reverse(CH::getOrder(ch)));
        } else if (orderType == "LevelTopDown") {
            return Order(CH::getLevelOrderTopDown(ch));
        } else if (orderType == "LevelBottomUp") {
            return Order(Vector::reverse(CH::getLevelOrderBottomUp(ch)));
        } else {
            Order order;
            IO::deserialize(orderType, order);
            Ensure(order.size() == ch.numVertices(), "Loaded order has size " << order.size() << ", but graph has " << ch.numVertices() << " vertices!");
            Vector::reverse(order);
            return order;
        }
    }
};

class VisualizeVertexOrder : public ParameterizedCommand {

public:
    VisualizeVertexOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "visualizeVertexOrder", "Draws the given graph with the vertices color-coded by order.") {
        addParameter("Input graph");
        addParameter("Visualization file");
        addParameter("Format", { "pdf", "png", "svg"});
        addParameter("Connect vertices?");
    }

    virtual void execute() noexcept {
        const std::string format = getParameter("Format");
        if (format == "pdf") {
            visualize<PDF>();
        } else if (format == "png") {
            visualize<PNG>();
        } else if (format == "svg") {
            visualize<SVG>();
        }
    }

private:
    template<typename FORMAT>
    inline void visualize() const noexcept {
        const std::string graphPath = getParameter("Input graph");
        RAPTOR::TransferGraph graph(graphPath);
        MapVisualization<FORMAT> visualization(getParameter("Visualization file"), Graph::boundingBox(graph));

        if (getParameter<bool>("Connect vertices?")) {
            for (size_t i = 1; i < graph.numVertices(); i++) {
                visualization.drawLine(graph.get(Coordinates, Vertex(i-1)), graph.get(Coordinates, Vertex(i)), Color::getGradientColor(Color::KITgreen, Color::KITorange, Color::KITred, (double)i/graph.numVertices()));
            }
        } else {
            for (const Vertex vertex : graph.vertices()) {
                visualization.drawPoint(graph.get(Coordinates, vertex), Color::getGradientColor(Color::KITgreen, Color::KITorange, Color::KITred, (double)vertex/graph.numVertices()), Icon::Dot, 16);
            }
        }
    }
};

class ReorderGraph : public ParameterizedCommand {

public:
    ReorderGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "reorderGraph", "Reorders the given graph according to the given vertex order.") {
        addParameter("Input graph");
        addParameter("Order");
        addParameter("Output graph");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Input graph"));
        Order order(getParameter("Order"));
        if (order.size() != graph.numVertices()) {
            std::cout << "Loaded order has size " << order.size() << ", but graph has " << graph.numVertices() << " vertices!" << std::endl;
            return;
        }
        graph.applyVertexOrder(order);
        graph.writeBinary(getParameter("Output graph"));
    }
};

class ReorderNetwork : public ParameterizedCommand {

public:
    ReorderNetwork(BasicShell& shell) :
        ParameterizedCommand(shell, "reorderNetwork", "Reorders the given public transit network according to the given vertex order.") {
        addParameter("Input network");
        addParameter("Network type", { "RAPTOR", "CSA", "Intermediate" });
        addParameter("Order");
        addParameter("Reorder stops?");
        addParameter("Output network");
    }

    virtual void execute() noexcept {
        const std::string networkType = getParameter("Network type");
        if (networkType == "RAPTOR") {
            reorder<RAPTOR::Data>();
        } else if (networkType == "CSA") {
            reorder<CSA::Data>();
        } else if (networkType == "Intermediate") {
            reorder<Intermediate::Data>();
        }
    }

private:
    template<typename DATA_TYPE>
    inline void reorder() const noexcept {
        DATA_TYPE data = DATA_TYPE::FromBinary(getParameter("Input network"));
        Order order(getParameter("Order"));
        if (order.size() != data.transferGraph.numVertices()) {
            std::cout << "Loaded order has size " << order.size() << ", but network has " << data.transferGraph.numVertices() << " vertices!" << std::endl;
            return;
        }
        data.applyVertexOrder(order, getParameter<bool>("Reorder stops?"));
        data.serialize(getParameter("Output network"));
    }
};

class ReorderNetworkStops : public ParameterizedCommand {

public:
    ReorderNetworkStops(BasicShell& shell) :
        ParameterizedCommand(shell, "reorderNetworkStops", "Reorders the stops of the given public transit network according to the given stop order.") {
        addParameter("Input network");
        addParameter("Network type", { "RAPTOR", "CSA", "Intermediate" });
        addParameter("Order");
        addParameter("Output network");
    }

    virtual void execute() noexcept {
        const std::string networkType = getParameter("Network type");
        if (networkType == "RAPTOR") {
            reorder<RAPTOR::Data>();
        } else if (networkType == "CSA") {
            reorder<CSA::Data>();
        } else if (networkType == "Intermediate") {
            reorder<Intermediate::Data>();
        }
    }

private:
    template<typename DATA_TYPE>
    inline void reorder() const noexcept {
        DATA_TYPE data = DATA_TYPE::FromBinary(getParameter("Input network"));
        Order order(getParameter("Order"));
        if (order.size() != data.numberOfStops()) {
            std::cout << "Loaded order has size " << order.size() << ", but network has " << data.numberOfStops() << " stops!" << std::endl;
            return;
        }
        data.applyStopOrder(order);
        data.serialize(getParameter("Output network"));
    }
};

class CreateRandomVertexOrder : public ParameterizedCommand {

public:
    CreateRandomVertexOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "createRandomVertexOrder", "Creates a random vertex order of the given size.") {
        addParameter("Size");
        addParameter("Output file");
        addParameter("Seed", "42");
    }

    virtual void execute() noexcept {
        Order order(Construct::Random, getParameter<size_t>("Size"), getParameter<size_t>("Seed"));
        order.serialize(getParameter("Output file"));
    }
};

class CreateDFSVertexOrder : public ParameterizedCommand {

public:
    CreateDFSVertexOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "createDFSVertexOrder", "Creates a vertex order from depth-first search on the given graph.") {
        addParameter("Input graph");
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Input graph"));
        getDFSVertexOrder(graph).serialize(getParameter("Output file"));
    }
};

class CreateDFSStopOrder : public ParameterizedCommand {

public:
    CreateDFSStopOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "createDFSStopOrder", "Creates a stop order from depth-first search on the given public transit network.") {
        addParameter("Input network");
        addParameter("Network type", { "RAPTOR", "CSA", "Intermediate" });
        addParameter("Output file");
    }

    virtual void execute() noexcept {
        const std::string networkType = getParameter("Network type");
        if (networkType == "RAPTOR") {
            createOrder<RAPTOR::Data>();
        } else if (networkType == "CSA") {
            createOrder<CSA::Data>();
        } else if (networkType == "Intermediate") {
            createOrder<Intermediate::Data>();
        }
    }

private:
    template<typename DATA_TYPE>
    inline void createOrder() const noexcept {
        DATA_TYPE data = DATA_TYPE::FromBinary(getParameter("Input network"));
        data.transferGraph.clear();
        getDFSVertexOrder(data.minTravelTimeGraph()).serialize(getParameter("Output file"));
    }
};

class CreateHilbertCurveVertexOrder : public ParameterizedCommand {

public:
    CreateHilbertCurveVertexOrder(BasicShell& shell) :
        ParameterizedCommand(shell, "createHilbertCurveVertexOrder", "Creates a vertex order from a Hilbert curve on the given graph.") {
        addParameter("Input graph");
        addParameter("Output file");
        addParameter("Max recursion depth", "32");
        addParameter("Visualization file", "");
    }

    virtual void execute() noexcept {
        RAPTOR::TransferGraph graph(getParameter("Input graph"));
        GraphHilbertCurveGenerator<RAPTOR::TransferGraph> curveGenerator(graph);
        curveGenerator.run(getParameter<size_t>("Max recursion depth"), true);
        curveGenerator.getVertexOrder().serialize(getParameter("Output file"));
        const std::string visualizationFile = getParameter("Visualization file");
        if (visualizationFile != "") {
            curveGenerator.template visualize<PDF>(visualizationFile);
        }
    }
};
