#pragma once

#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <optional>

#include "../../Shell/Shell.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Debug.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Calendar.h"
#include "../../Helpers/MultiThreading.h"
#include "../../Helpers/IO/File.h"
#include "../../Helpers/IO/CSVData.h"
#include "../../Helpers/IO/ParserCSV.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"
#include "../../Helpers/String/TextFileUtils.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Vector/Permutation.h"
#include "../../Helpers/Ranges/Range.h"

#include "../../DataStructures/Geometry/Rectangle.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Intermediate/Data.h"
#include "../../DataStructures/RAPTOR/Data.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/ABC/Data.h"

#include "../../Algorithms/CSA/CSA.h"
#include "../../Algorithms/CSA/ShortestPathCover.h"
#include "../../Algorithms/ABC/BlockBuilder.h"
#include "../../Algorithms/ABC/KMeans.h"

#include "../../Visualization/PDF.h"
#include "../../Visualization/Color.h"
#include "../../Visualization/MapVisualization.h"
#include "../../Visualization/CsaVisualization.h"

using namespace Shell;

class KMeansOnStops : public ParameterizedCommand {

public:
    KMeansOnStops(BasicShell& shell) :
        ParameterizedCommand(shell, "kMeansOnStops", "Uses kMeans to partition the stops of a CSA network.") {
        addParameter("CSA binary");
        addParameter("k");
        addParameter("Output file");
        addParameter("iterations");
        addParameter("dist iterations");
        addParameter("seed", "42");
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string outputFile = getParameter("Output file");
        const uint32_t k = getParameter<uint32_t>("k");
        const uint32_t iterations = getParameter<uint32_t>("iterations");
        const uint32_t distiterations = getParameter<uint32_t>("dist iterations");
        const uint32_t seed = getParameter<uint32_t>("seed");

        CSA::Data csa = CSA::Data::FromBinary(inputFile);
        csa.printInfo(false);
        ABC::KMeans partition(csa, k, seed, iterations, distiterations);
        CsaVisualization<PDF> doc(CsaVisualization<PDF>::FromCSA(outputFile, csa, 0.3));
        std::vector<std::vector<StopId>> stopsOfColumn(k);
        for (const StopId stop : csa.stops()) {
            if (partition.columnOfVertex[stop] < k) {
                stopsOfColumn[partition.columnOfVertex[stop]].emplace_back(stop);
            }
        }
        for (size_t i = 0; i < k; i++) {
            std::stringstream ss;
            ss << "Column " << i << ", weight: " << partition.labelOfColumn[i].weight << ", stops";
            doc.drawStops(stopsOfColumn[i], cyclicColor(i), cyclicIcon(i / 8), 20.0, ss.str());
        }
        // for (const StopId stop : csa.stops()) {
        //     for (const Edge edge : partition.graph.edgesFrom(stop)) {
        //         if (partition.graph.get(TravelTime, edge) > 600) continue;
        //         const Vertex to = partition.graph.get(ToVertex, edge);
        //         if (partition.columnOfVertex[stop] == partition.columnOfVertex[to]) continue;
        //         doc.drawLine(partition.graph.get(Coordinates, stop), partition.graph.get(Coordinates, to), Color::Black, 10.0);
        //     }
        // }
    }

};

class BuildMetisStopGraph : public ParameterizedCommand {

public:
    BuildMetisStopGraph(BasicShell& shell) :
        ParameterizedCommand(shell, "buildMetisStopGraph", "Exports the stop graph of a CSA network to metis format.") {
        addParameter("CSA binary");
        addParameter("Max edge time");
        addParameter("Output file");
        addParameter("Make connected", "false", {"true", "false"});
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string outputFile = getParameter("Output file");
        const int maxEdgeTime = String::parseSeconds(getParameter("Max edge time"));
        const bool makeConnected = getParameter<bool>("Make connected");

        CSA::Data csa = CSA::Data::FromBinary(inputFile);
        csa.printInfo(false);
        Ensure(csa.numberOfStops() >= csa.transferGraph.numVertices(), "Transfer graph contains to many vertices!");

        std::vector<int> connectionsPerStop(csa.numberOfStops(), 0);
        DynamicGraph<NoVertexAttributes, WithWeight> graph;
        graph.addVertices(csa.numberOfStops());
        for (const CSA::Connection& connection : csa.connections) {
            connectionsPerStop[connection.departureStopId]++;
            if (connection.travelTime() <= maxEdgeTime) {
                if (connection.departureStopId != connection.arrivalStopId) {
                    Edge forwardEdge = graph.findEdge(connection.departureStopId, connection.arrivalStopId);
                    if (graph.isEdge(forwardEdge)) {
                        graph.set(Weight, forwardEdge, graph.get(Weight, forwardEdge) + 1);
                    } else {
                        graph.addEdge(connection.departureStopId, connection.arrivalStopId).set(Weight, 1);
                    }
                    Edge backwardEdge = graph.findEdge(connection.arrivalStopId, connection.departureStopId);
                    if (graph.isEdge(backwardEdge)) {
                        graph.set(Weight, backwardEdge, graph.get(Weight, backwardEdge) + 1);
                    } else {
                        graph.addEdge(connection.arrivalStopId, connection.departureStopId).set(Weight, 1);
                    }
                }
            }
            for (const Edge transfer : csa.transferGraph.edgesFrom(connection.departureStopId)) {
                if (connection.travelTime() + csa.transferGraph.get(TravelTime, transfer) <= maxEdgeTime) {
                    const Vertex to = csa.transferGraph.get(ToVertex, transfer);
                    if (connection.departureStopId != to) {
                        Edge forwardEdge = graph.findEdge(connection.departureStopId, to);
                        if (graph.isEdge(forwardEdge)) {
                            graph.set(Weight, forwardEdge, graph.get(Weight, forwardEdge) + 1);
                        } else {
                            graph.addEdge(connection.departureStopId, to).set(Weight, 1);
                        }
                        Edge backwardEdge = graph.findEdge(to, connection.departureStopId);
                        if (graph.isEdge(backwardEdge)) {
                            graph.set(Weight, backwardEdge, graph.get(Weight, backwardEdge) + 1);
                        } else {
                            graph.addEdge(to, connection.departureStopId).set(Weight, 1);
                        }
                    }
                }
            }
        }
        if (makeConnected) {
            std::cout << "Analyzing components of the graph:" << std::endl;
            uint32_t maxComponent = 0;
            std::vector<uint64_t> sizeOfComponent(1, 0);
            std::vector<uint32_t> componentOfVertex(graph.numVertices(), 0);
            Dijkstra<DynamicGraph<NoVertexAttributes, WithWeight>, false> dijkstra(graph, graph[Weight]);
            for (const Vertex vertex : graph.vertices()) {
                if (componentOfVertex[vertex] != 0) continue;
                sizeOfComponent.emplace_back(0);
                dijkstra.run(vertex, noVertex, [&](const Vertex v){
                    componentOfVertex[v] = sizeOfComponent.size() - 1;
                    sizeOfComponent.back()++;
                });
                if (sizeOfComponent[maxComponent] < sizeOfComponent.back()) {
                    maxComponent = sizeOfComponent.size() - 1;
                }
            }
            std::cout << "   Number of components: " << String::prettyInt(sizeOfComponent.size() - 1) << std::endl;
            std::cout << "   Max component size: " << String::prettyInt(sizeOfComponent[maxComponent]) << std::endl;
            std::vector<Geometry::Point> coordinates = csa.transferGraph[Coordinates];
            for (const Vertex vertex : graph.vertices()) {
                if (componentOfVertex[vertex] == maxComponent) continue;
                coordinates[vertex] = Geometry::Point(Construct::XY, 100.0, 100.0);
            }
            Geometry::Rectangle boundingBox = Geometry::Rectangle::BoundingBox(coordinates);
            Geometry::GeoMetricAproximation metric = Geometry::GeoMetricAproximation::ComputeCorrection(boundingBox.center());
            CoordinateTree<Geometry::GeoMetricAproximation> ct(metric, coordinates);
            for (const Vertex vertex : graph.vertices()) {
                if (componentOfVertex[vertex] == maxComponent) continue;
                if (sizeOfComponent[componentOfVertex[vertex]] == 0) continue;
                const Vertex other = ct.getNearestNeighbor(csa.transferGraph.get(Coordinates, vertex));
                graph.addEdge(other, vertex).set(Weight, 1);
                graph.addEdge(vertex, other).set(Weight, 1);
                sizeOfComponent[componentOfVertex[vertex]] = 0;
            }
        }
        Graph::printInfo(graph);
        graph.printAnalysis();

        std::ofstream metis(outputFile);
        Assert(metis);
        Assert(metis.is_open());
        metis << graph.numVertices() << " " << (graph.numEdges() / 2) << " 011";
        for (const Vertex vertex : graph.vertices()) {
            metis << "\n" << std::max<int>(1, sqrt(connectionsPerStop[vertex]) + 0.5);
            for (const Edge edge : graph.edgesFrom(vertex)) {
                metis << " " << (graph.get(ToVertex, edge) + 1) << " " << graph.get(Weight, edge);
            }
        }
        metis.close();
    }

};

class ComputeShortestPathCover : public ParameterizedCommand {

public:
    ComputeShortestPathCover(BasicShell& shell) :
        ParameterizedCommand(shell, "computeShortestPathCover", "Draws the main stops of a CSA network.") {
        addParameter("CSA binary");
        addParameter("Output file");
        addParameter("Number of stops");
        addParameter("Number of samples");
        addParameter("Number of sample groups");
        addParameter("Max station transfer time");
        addParameter("Max tree time");
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string outputFile = getParameter("Output file");
        const uint32_t numberOfStops = getParameter<uint32_t>("Number of stops");
        const uint32_t numberOfSamples = getParameter<uint32_t>("Number of samples");
        const uint32_t numberOfSampleGroups = getParameter<uint32_t>("Number of sample groups");
        const int maxStationTransferTime = getParameter<int>("Max station transfer time");
        const int maxTreeTime = getParameter<int>("Max tree time");

        CSA::Data csa = CSA::Data::FromBinary(inputFile);
        csa.printInfo(false);

        CSA::ShortestPathCover spc(csa, numberOfStops, numberOfSamples, maxStationTransferTime, maxTreeTime, numberOfSampleGroups);

        IO::serialize(outputFile + ".binary", spc.getCover());

        CsaVisualization<PDF> doc(CsaVisualization<PDF>::FromCSA(outputFile + ".pdf", csa, 0.3));
        doc.drawConnectionsOnce(Color::Grey, 200);
        doc.drawStops(spc.getCover(), Color::KITred, Icon::MalteseCross, 60, "Shortest Path Cover");
    }

};

class LoadMetisPartition : public ParameterizedCommand {

public:
    LoadMetisPartition(BasicShell& shell) :
        ParameterizedCommand(shell, "loadMetisPartition", "Loads a metis partition of a CSA network.") {
        addParameter("CSA binary");
        addParameter("Metis file");
        addParameter("Max block size");
        addParameter("Max block time");
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string metisFile = getParameter("Metis file");
        const uint32_t maxBlockSize = getParameter<uint32_t>("Max block size");
        const int maxBlockTime = String::parseSeconds(getParameter("Max block time"));

        CSA::Data csa = CSA::Data::FromBinary(inputFile);
        csa.printInfo(false);
        Ensure(csa.numberOfStops() >= csa.transferGraph.numVertices(), "Transfer graph contains to many vertices!");

        std::vector<int> connectionsPerStop(csa.numberOfStops(), 0);
        for (const CSA::Connection& connection : csa.connections) {
            connectionsPerStop[connection.departureStopId]++;
        }

        std::vector<uint32_t> columnOfStop;
        std::vector<uint64_t> weightOfColumn;
        std::vector<std::vector<StopId>> stopsOfColumn;
        std::ifstream mateisData(metisFile);
        AssertMsg(mateisData.is_open(), "cannot open file: " << metisFile);
        while (!mateisData.eof()) {
            std::string line;
            getline(mateisData, line);
            if (line == "") continue;
            const uint32_t column = String::lexicalCast<uint32_t>(line);
            while (column >= stopsOfColumn.size()) {
                stopsOfColumn.emplace_back();
                weightOfColumn.emplace_back(0);
            }
            stopsOfColumn[column].emplace_back(StopId(columnOfStop.size()));
            weightOfColumn[column] += connectionsPerStop[columnOfStop.size()];
            columnOfStop.emplace_back(column);
        }

        std::cout << "Building block graph:" << std::endl;
        ABC::BlockBuilder bb(csa, columnOfStop, maxBlockSize, maxBlockTime);
        std::cout << std::endl;
        std::vector<size_t> blockSizes = bb.blockSizes();
        sort(blockSizes);
        std::cout << "   Number of components: " << String::prettyInt(bb.numberOfConnections()) << std::endl;
        std::cout << "   Number of blocks: " << String::prettyInt(bb.numberOfBlocks()) << std::endl;
        std::cout << "      Min connections per block:  " << String::prettyInt(Vector::min(blockSizes)) << std::endl;
        std::cout << "      Mean connections per block: " << String::prettyInt(Vector::mean(blockSizes)) << std::endl;
        std::cout << "      Max connections per block:  " << String::prettyInt(Vector::max(blockSizes)) << std::endl;
        std::cout << "      25% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.25)) << std::endl;
        std::cout << "      50% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.5)) << std::endl;
        std::cout << "      75% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.75)) << std::endl;
        std::vector<size_t> blockDegrees = bb.blockDegrees();
        sort(blockDegrees);
        std::cout << "   Number of edges: " << String::prettyInt(bb.getDynamicBlockGraph().numEdges()) << std::endl;
        std::cout << "      Min edges per block:        " << String::prettyInt(Vector::min(blockDegrees)) << std::endl;
        std::cout << "      Mean edges per block:       " << String::prettyInt(Vector::mean(blockDegrees)) << std::endl;
        std::cout << "      Max edges per block:        " << String::prettyInt(Vector::max(blockDegrees)) << std::endl;
        std::cout << "      25% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.25)) << std::endl;
        std::cout << "      50% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.5)) << std::endl;
        std::cout << "      75% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.75)) << std::endl;
        if (Graph::isAcyclic(bb.getDynamicBlockGraph())) {
            std::cout << "   Block graph is a DAG." << std::endl;
        } else {
            std::cout << "   Block graph is not a DAG!" << std::endl;
        }
        Graph::printInfo(bb.getDynamicBlockGraph());
        bb.getDynamicBlockGraph().printAnalysis();

        CsaVisualization<PDF> doc(CsaVisualization<PDF>::FromCSA(metisFile + ".pdf", csa, 0.3));
        for (size_t i = 0; i < stopsOfColumn.size(); i++) {
            std::stringstream ss;
            ss << "Column " << i << ", weight: " << weightOfColumn[i] << ", stops";
            doc.drawStops(stopsOfColumn[i], cyclicColor(i), cyclicIcon(i / 8), 20.0, ss.str());
        }
    }

};

class LoadMetisPartitions : public ParameterizedCommand {

public:
    LoadMetisPartitions(BasicShell& shell) :
        ParameterizedCommand(shell, "loadMetisPartitions", "Loads multiple metis partitions of a CSA network.") {
        addParameter("CSA binary");
        addParameter("Metis file");
        addParameter("Results file");
    }

    virtual void execute() noexcept {
        const std::string inputFile = getParameter("CSA binary");
        const std::string metisFile = getParameter("Metis file");
        const std::string resultsFile = getParameter("Results file");

        const std::vector<int> edgeTimes{10, 20, 30, 40};
        const std::vector<int> parts{25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
        const std::vector<int> maxBlockTimes{10, 20, 30, 40, 50, 60, 90, 120};
        const std::vector<uint32_t> maxBlockSizes{512, 1024, 1536, 2048, 3072, 4096};
        const std::vector<bool> blockEdgeModes{true, false};

        CSA::Data csa = CSA::Data::FromBinary(inputFile);
        csa.printInfo(false);
        Ensure(csa.numberOfStops() >= csa.transferGraph.numVertices(), "Transfer graph contains to many vertices!");

        std::vector<int> connectionsPerStop(csa.numberOfStops(), 0);
        for (const CSA::Connection& connection : csa.connections) {
            connectionsPerStop[connection.departureStopId]++;
        }

        std::ofstream result(resultsFile);
        Assert(result);
        Assert(result.is_open());
        result << "EdgeTime,Parts,BlockTime,BlockSize,BlockMode,PathCover,Blocks,Connections,ConnectionsMean,ConnectionsMax,ConnectionsMedian,Edges,EdgesMean,EdgesMax,EdgesMedian,Time\n";

        for (const int edgeTime : edgeTimes) {
            for (const int part : parts) {
                std::stringstream fullMetisFile;
                fullMetisFile << metisFile << edgeTime << "min.metis.part." << part;
                const std::string fullMetisFileName = fullMetisFile.str();
                std::vector<uint32_t> columnOfStop;
                std::vector<uint64_t> weightOfColumn;
                std::vector<std::vector<StopId>> stopsOfColumn;
                std::ifstream mateisData(fullMetisFileName);
                AssertMsg(mateisData.is_open(), "cannot open file: " << fullMetisFileName);
                while (!mateisData.eof()) {
                    std::string line;
                    getline(mateisData, line);
                    if (line == "") continue;
                    const uint32_t column = String::lexicalCast<uint32_t>(line);
                    while (column >= stopsOfColumn.size()) {
                        stopsOfColumn.emplace_back();
                        weightOfColumn.emplace_back(0);
                    }
                    stopsOfColumn[column].emplace_back(StopId(columnOfStop.size()));
                    weightOfColumn[column] += connectionsPerStop[columnOfStop.size()];
                    columnOfStop.emplace_back(column);
                }
                for (const int maxBlockTime : maxBlockTimes) {
                    for (const int maxBlockSize : maxBlockSizes) {
                        for (const bool blockEdgeMode : blockEdgeModes) {
                            std::cout << std::endl << "Building block graph:" << std::endl;
                            Timer timer;
                            ABC::BlockBuilder bb(csa, columnOfStop, maxBlockSize, maxBlockTime * 60, true, blockEdgeMode);
                            const double time = timer.elapsedMilliseconds();
                            std::vector<size_t> blockSizes = bb.blockSizes();
                            sort(blockSizes);
                            std::cout << "   Number of connections: " << String::prettyInt(bb.numberOfConnections()) << std::endl;
                            std::cout << "   Number of blocks: " << String::prettyInt(bb.numberOfBlocks()) << std::endl;
                            std::cout << "      Min connections per block:  " << String::prettyInt(Vector::min(blockSizes)) << std::endl;
                            std::cout << "      Mean connections per block: " << String::prettyInt(Vector::mean(blockSizes)) << std::endl;
                            std::cout << "      Max connections per block:  " << String::prettyInt(Vector::max(blockSizes)) << std::endl;
                            std::cout << "      25% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.25)) << std::endl;
                            std::cout << "      50% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.5)) << std::endl;
                            std::cout << "      75% connections per block:  " << String::prettyInt(Vector::percentile(blockSizes, 0.75)) << std::endl;
                            std::vector<size_t> blockDegrees = bb.blockDegrees();
                            sort(blockDegrees);
                            std::cout << "   Number of edges: " << String::prettyInt(bb.getDynamicBlockGraph().numEdges()) << std::endl;
                            std::cout << "      Min edges per block:        " << String::prettyInt(Vector::min(blockDegrees)) << std::endl;
                            std::cout << "      Mean edges per block:       " << String::prettyInt(Vector::mean(blockDegrees)) << std::endl;
                            std::cout << "      Max edges per block:        " << String::prettyInt(Vector::max(blockDegrees)) << std::endl;
                            std::cout << "      25% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.25)) << std::endl;
                            std::cout << "      50% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.5)) << std::endl;
                            std::cout << "      75% edges per block:        " << String::prettyInt(Vector::percentile(blockDegrees, 0.75)) << std::endl;
                            if (Graph::isAcyclic(bb.getDynamicBlockGraph())) {
                                std::cout << "   Block graph is a DAG." << std::endl;
                            } else {
                                std::cout << "   Block graph is not a DAG!" << std::endl;
                            }
                            result << edgeTime << "," << part << "," << maxBlockTime << "," << maxBlockSize << "," << blockEdgeMode << ",";
                            result << bb.getPathCoverSize() << "," << bb.numberOfBlocks() << ",";
                            result << bb.numberOfConnections() << "," << Vector::mean(blockSizes) << "," << Vector::max(blockSizes) << "," << Vector::median(blockSizes) << ",";
                            result << bb.getDynamicBlockGraph().numEdges() << "," << Vector::mean(blockDegrees) << "," << Vector::max(blockDegrees) << "," << Vector::median(blockDegrees) << ",";
                            result << time << "\n" << std::flush;
                        }
                    }
                }
                CsaVisualization<PDF> doc(CsaVisualization<PDF>::FromCSA(fullMetisFileName + ".pdf", csa, 0.3));
                for (size_t i = 0; i < stopsOfColumn.size(); i++) {
                    std::stringstream ss;
                    ss << "Column " << i << ", weight: " << weightOfColumn[i] << ", stops";
                    doc.drawStops(stopsOfColumn[i], cyclicColor(i), cyclicIcon(i / 8), 20.0, ss.str());
                }
            }
        }
    }

};
