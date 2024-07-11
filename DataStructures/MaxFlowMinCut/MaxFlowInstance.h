#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Console/ProgressBar.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/Meta.h"

#include "../Graph/Graph.h"

#include "FlowUtils.h"

template<typename FLOW_TYPE = int>
class StaticMaxFlowInstance {
public:
    using FlowType = FLOW_TYPE;
    using GraphType = ParametricFlowGraph<FlowType>;

    StaticMaxFlowInstance() : source(noVertex), sink(noVertex) {}

    explicit StaticMaxFlowInstance(const std::string& fileName) {
        deserialize(fileName);
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, source, sink);
        graph.writeBinary(fileName + ".graph");
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, source, sink);
        graph.readBinary(fileName + ".graph");
    }

    inline FlowType getCapacity(const Edge edge) const noexcept {
        return graph.get(Capacity, edge);
    }

    inline const std::vector<FlowType>& getCurrentCapacities() const noexcept {
        return graph.get(Capacity);
    }

    inline const std::vector<FlowType>& getSourceDiff() const noexcept {
        throw std::runtime_error("getSourceDiff() not available in static max-flow instance!");
    }

    inline const std::vector<FlowType>& getSinkDiff() const noexcept {
        throw std::runtime_error("getSinkDiff() not available in static max-flow instance!");
    }

    template<bool VERBOSE = true>
    inline void fromDimacs(const std::string& fileName) noexcept {
        const std::string fileNameWithExtension = FileSystem::ensureExtension(fileName, ".max");
        if constexpr (VERBOSE) std::cout << "Reading DIMACS max-flow graph from: " << fileNameWithExtension << std::endl << std::flush;
        std::ifstream is(fileNameWithExtension);
        Assert(is.is_open(), "Cannot open file: " << fileNameWithExtension);
        ProgressBar bar(1);
        std::cout << "\r                     \r" << std::flush;

        ParametricFlowGraphEdgeList<FlowType> temp;
        std::vector<bool> hasEdgeFromSource;
        std::vector<bool> hasEdgeFromSink;
        size_t vertexCount = -1;
        size_t edgeCount = -1;
        source = noVertex;
        sink = noVertex;

        while (!is.eof()) {
            std::string line;
            getline(is, line);
            line = String::trim(line);
            if (line.empty() || line[0] == 'c') continue;
            const std::vector<std::string> tokens = String::split(line, ' ');
            if (vertexCount == size_t(-1)) {
                if (tokens.size() != 4 || tokens[0] != "p" || tokens[1] != "max") {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else {
                    vertexCount = String::lexicalCast<size_t>(tokens[2]);
                    hasEdgeFromSource.resize(vertexCount, false);
                    hasEdgeFromSink.resize(vertexCount, false);
                    edgeCount = String::lexicalCast<size_t>(tokens[3]);
                    temp.reserve(vertexCount, edgeCount);
                    temp.addVertices(vertexCount);
                    if constexpr (VERBOSE) bar.init(edgeCount);
                }
            } else if (source == noVertex || sink == noVertex) {
                if (tokens.size() != 3 || tokens[0] != "n" || (tokens[2] != "s" && tokens[2] != "t")) {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else if (tokens[2] == "s") {
                    source = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    hasEdgeFromSource[source] = true;
                    if (!temp.isVertex(source)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                } else {
                    sink = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    hasEdgeFromSink[sink] = true;
                    if (!temp.isVertex(sink)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                }
            } else {
                if (tokens.size() != 4 || tokens[0] != "a") {
                    std::cout << "WARNING, ignoring line in .max-file: " << line << std::endl;
                    continue;
                } else {
                    const Vertex from(String::lexicalCast<size_t>(tokens[1]) - 1);
                    const Vertex to(String::lexicalCast<size_t>(tokens[2]) - 1);
                    const FlowType capacity = String::lexicalCast<FlowType>(tokens[3]);
                    if (from == source) hasEdgeFromSource[to] = true;
                    else if (from == sink) hasEdgeFromSink[to] = true;
                    if (!temp.isVertex(from)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    } else if (!temp.isVertex(to)) {
                        std::cout << "ERROR, " << tokens[2] << " does not name a vertex!" << std::endl;
                        break;
                    }
                    temp.addEdge(from, to).set(Capacity, capacity);
                    temp.addEdge(to, from).set(Capacity, 0);
                    if constexpr (VERBOSE) bar++;
                }
            }
        }
        is.close();
        if constexpr (VERBOSE) std::cout << std::endl;
        if (temp.numEdges() / 2 != edgeCount) {
            std::cout << "WARNING, found " << temp.numEdges() / 2 << " edges, but " << edgeCount << " edges were declared." << std::endl;
        }

        for (const Vertex vertex: temp.vertices()) {
            if (!hasEdgeFromSource[vertex]) {
                temp.addEdge(source, vertex).set(Capacity, 0);
                temp.addEdge(vertex, source).set(Capacity, 0);
            }
            if (!hasEdgeFromSink[vertex]) {
                temp.addEdge(sink, vertex).set(Capacity, 0);
                temp.addEdge(vertex, sink).set(Capacity, 0);
            }
        }

        temp.sortEdges();
        ParametricFlowGraphEdgeList<FlowType> temp2;
        temp2.reserve(vertexCount, edgeCount);
        temp2.addVertices(vertexCount);
        Vertex lastFrom = noVertex;
        Vertex lastTo = noVertex;
        Edge lastEdge = noEdge;
        for (const Edge edge : temp.edges()) {
            const Vertex from = temp.get(FromVertex, edge);
            const Vertex to = temp.get(ToVertex, edge);
            if (from == lastFrom && to == lastTo) {
                temp2.get(Capacity, lastEdge) += temp.get(Capacity, edge);
            } else {
                lastFrom = from;
                lastTo = to;
                lastEdge = temp2.addEdge(from, to).set(Capacity, temp.get(Capacity, edge));
            }
        }
        Graph::move(std::move(temp2), graph);
    }

public:
    GraphType graph;
    Vertex source;
    Vertex sink;
};

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class ParametricMaxFlowInstance {
public:
    using FlowFunction = FLOW_FUNCTION;
    using FlowType = FlowFunction::FlowType;
    using GraphType = ParametricFlowGraph<FlowFunction>;

    ParametricMaxFlowInstance() : source(noVertex), sink(noVertex), alphaMin(0), alphaMax(INFTY) {}

    ParametricMaxFlowInstance(const GraphType& graph, const Vertex source, const Vertex sink, const double alphaMin, const double alphaMax) :
        graph(graph), source(source), sink(sink), alphaMin(alphaMin), alphaMax(alphaMax) {
    }

    explicit ParametricMaxFlowInstance(const std::string& fileName) {
        deserialize(fileName);
    }

    template<typename T>
    explicit ParametricMaxFlowInstance(const StaticMaxFlowInstance<T>& staticInstance) : source(staticInstance.source), sink(staticInstance.sink), alphaMin(0), alphaMax(2) {
        Graph::copy(staticInstance.graph, graph);
        for (const Edge e : graph.edges()) {
            graph.set(Capacity, e, FlowFunction(staticInstance.graph.get(Capacity, e)));
        }
        FlowType minSourceCapacity = INFTY;
        FlowType maxSourceCapacity = 0;
        for (const Edge edge : graph.edgesFrom(source)) {
            const FlowType capacity = staticInstance.graph.get(Capacity, edge);
            minSourceCapacity = std::min(capacity, minSourceCapacity);
            maxSourceCapacity = std::max(capacity, maxSourceCapacity);
            graph.set(Capacity, edge, FlowFunction(capacity, 0));
        }
        std::cout << "Min source capacity: " << minSourceCapacity << std::endl;
        std::cout << "Max source capacity: " << maxSourceCapacity << std::endl;
        FlowType minSinkCapacity = INFTY;
        FlowType maxSinkCapacity = 0;
        for (const Edge edge : graph.edgesFrom(sink)) {
            const Edge reverseEdge = graph.get(ReverseEdge, edge);
            const FlowType capacity = staticInstance.graph.get(Capacity, reverseEdge);
            minSinkCapacity = std::min(capacity, minSinkCapacity);
            maxSinkCapacity = std::max(capacity, maxSinkCapacity);
            graph.set(Capacity, reverseEdge, FlowFunction(-capacity, 2 * capacity));
        }
        std::cout << "Min sink capacity: " << minSinkCapacity << std::endl;
        std::cout << "Max sink capacity: " << maxSinkCapacity << std::endl;
    }

    ParametricMaxFlowInstance(const StaticMaxFlowInstance<int>& staticInstance, const double edgeProbability, const bool sinkEdges) :
        source(staticInstance.source),
        sink(staticInstance.sink),
        alphaMin(0),
        alphaMax(INFTY) {
        Graph::copy(staticInstance.graph, graph);
        int maxCapacity = 0;
        for (const Edge e : graph.edges()) {
            const int capacity = staticInstance.graph.get(Capacity, e);
            if (capacity < INFTY) maxCapacity = std::max(maxCapacity, capacity);
            graph.set(Capacity, e, FlowFunction(capacity));
        }
        std::cout << "Max capacity: " << maxCapacity << std::endl;
        std::mt19937 randomGenerator;
        std::uniform_int_distribution<> distribution(1, maxCapacity);
        std::uniform_real_distribution<> probDist(0, 1);
        for (const Edge edge : graph.edgesFrom(source)) {
            if (probDist(randomGenerator) > edgeProbability) {
                graph.set(Capacity, edge, FlowFunction(0));
            } else {
                const int slope = distribution(randomGenerator);
                const int minValue = distribution(randomGenerator);
                graph.set(Capacity, edge, FlowFunction(slope, minValue));
            }
        }
        if (sinkEdges) {
            alphaMax = 1;
            for (const Edge edge : graph.edgesFrom(sink)) {
                const Edge reverseEdge = graph.get(ReverseEdge, edge);
                if (probDist(randomGenerator) > edgeProbability) {
                    graph.set(Capacity, reverseEdge, FlowFunction(0));
                } else {
                    const int slope = distribution(randomGenerator);
                    const int maxValue = distribution(randomGenerator);
                    const double limit = 1 + static_cast<double>(maxValue)/slope;
                    alphaMax = std::min(limit, alphaMax);
                    graph.set(Capacity, reverseEdge, FlowFunction(-slope, maxValue + slope));
                }
            }
        }
    }

    inline void serialize(const std::string& fileName) const noexcept {
        IO::serialize(fileName, source, sink, alphaMin, alphaMax);
        graph.writeBinary(fileName + ".graph");
    }

    inline void deserialize(const std::string& fileName) noexcept {
        IO::deserialize(fileName, source, sink, alphaMin, alphaMax);
        graph.readBinary(fileName + ".graph");
    }

    inline const FlowFunction& getCapacity(const Edge edge) const noexcept {
        return graph.get(Capacity, edge);
    }

    inline FlowType getCapacity(const Edge edge, const double alpha) const noexcept {
        return graph.get(Capacity, edge).eval(alpha);
    }

    inline const std::vector<FlowFunction>& getCurrentCapacities() const noexcept {
        return graph.get(Capacity);
    }

    template<bool VERBOSE = true>
    inline void fromDimacs(const std::string& fileName) noexcept {
        const std::string fileNameWithExtension = FileSystem::ensureExtension(fileName, ".max");
        if constexpr (VERBOSE) std::cout << "Reading DIMACS max-flow graph from: " << fileNameWithExtension << std::endl << std::flush;
        std::ifstream is(fileNameWithExtension);
        Assert(is.is_open(), "Cannot open file: " << fileNameWithExtension);
        ProgressBar bar(1);
        std::cout << "\r                     \r" << std::flush;

        ParametricFlowGraphEdgeList<FlowFunction> temp;
        std::vector<bool> hasEdgeFromSource;
        std::vector<bool> hasEdgeFromSink;
        size_t vertexCount = -1;
        size_t edgeCount = -1;
        source = noVertex;
        sink = noVertex;

        while (!is.eof()) {
            std::string line;
            getline(is, line);
            line = String::trim(line);
            if (line.empty() || line[0] == 'c') continue;
            const std::vector<std::string> tokens = String::split(line, ' ');
            if (vertexCount == size_t(-1)) {
                if (tokens.size() != 4 || tokens[0] != "p" || tokens[1] != "pmax") {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else {
                    vertexCount = String::lexicalCast<size_t>(tokens[2]);
                    hasEdgeFromSource.resize(vertexCount, false);
                    hasEdgeFromSink.resize(vertexCount, false);
                    edgeCount = String::lexicalCast<size_t>(tokens[3]);
                    temp.reserve(vertexCount, edgeCount);
                    temp.addVertices(vertexCount);
                    if constexpr (VERBOSE) bar.init(edgeCount);
                }
            } else if (source == noVertex || sink == noVertex) {
                if (tokens.size() != 3 || tokens[0] != "n" || (tokens[2] != "s" && tokens[2] != "t")) {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else if (tokens[2] == "s") {
                    source = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    hasEdgeFromSource[source] = true;
                    if (!temp.isVertex(source)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                } else {
                    sink = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    hasEdgeFromSource[sink] = true;
                    if (!temp.isVertex(sink)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                }
            } else {
                if (tokens.size() != 5 || tokens[0] != "a") {
                    std::cout << "WARNING, ignoring line in .max-file: " << line << std::endl;
                    continue;
                } else {
                    const Vertex from(String::lexicalCast<size_t>(tokens[1]) - 1);
                    const Vertex to(String::lexicalCast<size_t>(tokens[2]) - 1);
                    const FlowType capacityA = String::lexicalCast<FlowType>(tokens[3]);
                    const FlowType capacityB = String::lexicalCast<FlowType>(tokens[4]);
                    if (from == source) hasEdgeFromSource[to] = true;
                    else if (from == sink) hasEdgeFromSink[to] = true;
                    if (!temp.isVertex(from)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    } else if (!temp.isVertex(to)) {
                        std::cout << "ERROR, " << tokens[2] << " does not name a vertex!" << std::endl;
                        break;
                    }
                    temp.addEdge(from, to).set(Capacity, FlowFunction(capacityA, capacityB));
                    temp.addEdge(to, from).set(Capacity, FlowFunction(0));
                    if constexpr (VERBOSE) bar++;
                }
            }
        }
        is.close();
        if constexpr (VERBOSE) std::cout << std::endl;
        if (temp.numEdges() / 2 != edgeCount) {
            std::cout << "WARNING, found " << temp.numEdges() / 2 << " edges, but " << edgeCount << " edges were declared." << std::endl;
        }

        for (const Vertex vertex: temp.vertices()) {
            if (!hasEdgeFromSource[vertex]) {
                temp.addEdge(source, vertex).set(Capacity, FlowFunction(0));
                temp.addEdge(vertex, source).set(Capacity, FlowFunction(0));
            }
            if (!hasEdgeFromSink[vertex]) {
                temp.addEdge(sink, vertex).set(Capacity, FlowFunction(0));
                temp.addEdge(vertex, sink).set(Capacity, FlowFunction(0));
            }
        }

        temp.sortEdges();
        ParametricFlowGraphEdgeList<FlowFunction> temp2;
        temp2.reserve(vertexCount, edgeCount);
        temp2.addVertices(vertexCount);
        Vertex lastFrom = noVertex;
        Vertex lastTo = noVertex;
        Edge lastEdge = noEdge;
        for (const Edge edge : temp.edges()) {
            const Vertex from = temp.get(FromVertex, edge);
            const Vertex to = temp.get(ToVertex, edge);
            if (from == lastFrom && to == lastTo) {
                temp2.get(Capacity, lastEdge) += temp.get(Capacity, edge);
            } else {
                lastFrom = from;
                lastTo = to;
                lastEdge = temp2.addEdge(from, to).set(Capacity, temp.get(Capacity, edge));
            }
        }
        Graph::move(std::move(temp2), graph);
    }

public:
    GraphType graph;
    Vertex source;
    Vertex sink;
    double alphaMin;
    double alphaMax;
};

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class RestartableMaxFlowWrapper {
public:
    using FlowFunction = FLOW_FUNCTION;
    using InstanceType = ParametricMaxFlowInstance<FlowFunction>;
    using GraphType = InstanceType::GraphType;
    using FlowType = InstanceType::FlowType;

public:
    RestartableMaxFlowWrapper(const InstanceType& instance, const double alpha) :
        instance(instance),
        graph(instance.graph),
        source(instance.source),
        sink(instance.sink),
        alpha(alpha),
        currentCapacity(graph.numEdges()),
        sourceDiff(graph.outDegree(source), 0),
        sinkDiff(graph.outDegree(sink), 0) {
        for (const Edge edge : instance.graph.edges()) {
            currentCapacity[edge] = instance.getCapacity(edge, alpha);
        }
    }

    RestartableMaxFlowWrapper(const InstanceType& instance) :
        RestartableMaxFlowWrapper(instance, instance.alphaMin) {}


    inline const FlowType& getCapacity(const Edge edge) const noexcept {
        return currentCapacity[edge];
    }

    inline const std::vector<FlowType>& getCurrentCapacities() const noexcept {
        return currentCapacity;
    }

    inline const std::vector<FlowType>& getSourceDiff() const noexcept {
        return sourceDiff;
    }

    inline const std::vector<FlowType>& getSinkDiff() const noexcept {
        return sinkDiff;
    }

    inline void setAlpha(const double newAlpha) noexcept {
        Assert(newAlpha >= alpha, "Parameter did not increase!");
        Assert(newAlpha <= instance.alphaMax, "Parameter out of range!");
        alpha = newAlpha;
        Edge edgeFromSource = graph.beginEdgeFrom(source);
        for (size_t i = 0; i < sourceDiff.size(); i++, edgeFromSource++) {
            const FlowType newCapacity = instance.getCapacity(edgeFromSource).eval(alpha);
            sourceDiff[i] = newCapacity - currentCapacity[edgeFromSource];
            Assert(sourceDiff[i] >= 0, "Capacity of source-incident edge has decreased!");
            currentCapacity[edgeFromSource] = newCapacity;
        }
        Edge edgeFromSink = graph.beginEdgeFrom(sink);
        for (size_t i = 0; i < sinkDiff.size(); i++, edgeFromSink++) {
            const Edge edgeToSink = graph.get(ReverseEdge, edgeFromSink);
            const FlowType newCapacity = instance.getCapacity(edgeToSink).eval(alpha);
            Assert(newCapacity >= 0, "Negative capacity!");
            sinkDiff[i] = newCapacity - currentCapacity[edgeToSink];
            Assert(sinkDiff[i] <= 0, "Capacity of sink-incident edge has increased!");
            currentCapacity[edgeToSink] = newCapacity;
        }
    }

public:
    const InstanceType& instance;
    const GraphType& graph;
    const Vertex& source;
    const Vertex& sink;
    double alpha;
    std::vector<FlowType> currentCapacity;
    std::vector<FlowType> sourceDiff;
    std::vector<FlowType> sinkDiff;
};

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class ChordSchemeMaxFlowWrapper {
public:
    using FlowFunction = FLOW_FUNCTION;
    using Type = ChordSchemeMaxFlowWrapper<FlowFunction>;
    using InstanceType = ParametricMaxFlowInstance<FlowFunction>;
    using GraphType = InstanceType::GraphType;
    using FlowType = InstanceType::FlowType;

public:
    ChordSchemeMaxFlowWrapper(const InstanceType& instance, const double alpha) :
        graph(instance.graph),
        source(instance.source),
        sink(instance.sink),
        alpha(alpha),
        currentCapacity(instance.graph.numEdges()),
        newToOldVertex(Vector::id<Vertex>(graph.numVertices())) {
        for (const Edge edge : instance.graph.edges()) {
            currentCapacity[edge] = instance.getCapacity(edge, alpha);
        }
    }

    ChordSchemeMaxFlowWrapper(const InstanceType& instance) :
        ChordSchemeMaxFlowWrapper(instance, instance.alphaMin) {}

    ChordSchemeMaxFlowWrapper(const Type& parent, const GraphType& graph, const Vertex source, const Vertex sink, const std::vector<Vertex>& newToOldVertex) :
        graph(graph),
        source(source),
        sink(sink),
        alpha(parent.alpha),
        currentCapacity(graph.numEdges()),
        newToOldVertex(newToOldVertex) {
    }

    inline const FlowType& getCapacity(const Edge edge) const noexcept {
        return currentCapacity[edge];
    }

    inline FlowFunction getCapacity(const std::vector<Edge>& edges) const noexcept {
        FlowFunction result(0);
        for (const Edge edge : edges) {
            result += graph.get(Capacity, edge);
        }
        return result;
    }

    inline const std::vector<FlowType>& getCurrentCapacities() const noexcept {
        return currentCapacity;
    }

    inline void setAlpha(const double newAlpha) noexcept {
        alpha = newAlpha;
        for (const Edge edge: graph.edges()) {
            currentCapacity[edge] = graph.get(Capacity, edge).eval(alpha);
        }
    }

    inline Type contractSourceComponent(const std::vector<bool>& inSinkComponent) const noexcept {
        GraphType contractedGraph;
        std::vector<Vertex> newToOldMapping;
        Vertex newSource = noVertex;
        Vertex newSink = noVertex;
        Graph::copy(graph, contractedGraph);
        std::vector<Edge> fromSourceEdge(contractedGraph.numVertices(), noEdge);
        std::vector<Edge> toSourceEdge(contractedGraph.numVertices(), noEdge);
        for (const Edge edge : contractedGraph.edgesFrom(source)) {
            const Vertex to = contractedGraph.get(ToVertex, edge);
            if (!inSinkComponent[to]) continue;
            fromSourceEdge[to] = edge;
            toSourceEdge[to] = contractedGraph.get(ReverseEdge, edge);
        }
        for (const Vertex vertex : contractedGraph.vertices()) {
            if (inSinkComponent[vertex] || vertex == source) {
                if (vertex == source) newSource = Vertex(newToOldMapping.size());
                else if (vertex == sink) newSink = Vertex(newToOldMapping.size());
                newToOldMapping.emplace_back(newToOldVertex[vertex]);
                continue;
            }
            for (const Edge edge : contractedGraph.edgesFrom(vertex)) {
                const Vertex to = contractedGraph.get(ToVertex, edge);
                if (!inSinkComponent[to]) continue;
                contractedGraph.get(Capacity, fromSourceEdge[to]) += contractedGraph.get(Capacity, edge);
                const Edge reverseEdge = contractedGraph.get(ReverseEdge, edge);
                contractedGraph.get(Capacity, toSourceEdge[to]) += contractedGraph.get(Capacity, reverseEdge);
            }
        }
        contractedGraph.deleteVertices([&](const Vertex vertex) {
            return !inSinkComponent[vertex] && vertex != source;
        });
        return {*this, contractedGraph, newSource, newSink, newToOldMapping};
    }

    inline Type contractSinkComponent(const std::vector<bool>& inSinkComponent) const noexcept {
        GraphType contractedGraph;
        std::vector<Vertex> newToOldMapping;
        Vertex newSource = noVertex;
        Vertex newSink = noVertex;
        Graph::copy(graph, contractedGraph);
        std::vector<Edge> fromSinkEdge(contractedGraph.numVertices(), noEdge);
        std::vector<Edge> toSinkEdge(contractedGraph.numVertices(), noEdge);
        for (const Edge edge : contractedGraph.edgesFrom(sink)) {
            const Vertex to = contractedGraph.get(ToVertex, edge);
            if (inSinkComponent[to]) continue;
            fromSinkEdge[to] = edge;
            toSinkEdge[to] = contractedGraph.get(ReverseEdge, edge);
        }
        for (const Vertex vertex : contractedGraph.vertices()) {
            if (!inSinkComponent[vertex] || vertex == sink) {
                if (vertex == source) newSource = Vertex(newToOldMapping.size());
                else if (vertex == sink) newSink = Vertex(newToOldMapping.size());
                newToOldMapping.emplace_back(newToOldVertex[vertex]);
                continue;
            }
            for (const Edge edge : contractedGraph.edgesFrom(vertex)) {
                const Vertex to = contractedGraph.get(ToVertex, edge);
                if (inSinkComponent[to]) continue;
                contractedGraph.get(Capacity, fromSinkEdge[to]) += contractedGraph.get(Capacity, edge);
                const Edge reverseEdge = contractedGraph.get(ReverseEdge, edge);
                contractedGraph.get(Capacity, toSinkEdge[to]) += contractedGraph.get(Capacity, reverseEdge);
            }
        }
        contractedGraph.deleteVertices([&](const Vertex vertex) {
            return inSinkComponent[vertex] && vertex != sink;
        });
        return {*this, contractedGraph, newSource, newSink, newToOldMapping};
    }

    inline Type contractSourceAndSinkComponents(const std::vector<bool>& inSinkComponent1, const std::vector<bool>& inSinkComponent2) const noexcept {
        GraphType contractedGraph;
        std::vector<Vertex> newToOldMapping;
        Vertex newSource = noVertex;
        Vertex newSink = noVertex;
        Graph::copy(graph, contractedGraph);
        std::vector<Edge> fromSourceEdge(contractedGraph.numVertices(), noEdge);
        std::vector<Edge> toSourceEdge(contractedGraph.numVertices(), noEdge);
        for (const Edge edge : contractedGraph.edgesFrom(source)) {
            const Vertex to = contractedGraph.get(ToVertex, edge);
            if (!inSinkComponent1[to]) continue;
            fromSourceEdge[to] = edge;
            toSourceEdge[to] = contractedGraph.get(ReverseEdge, edge);
        }
        std::vector<Edge> fromSinkEdge(contractedGraph.numVertices(), noEdge);
        std::vector<Edge> toSinkEdge(contractedGraph.numVertices(), noEdge);
        for (const Edge edge : contractedGraph.edgesFrom(sink)) {
            const Vertex to = contractedGraph.get(ToVertex, edge);
            if (inSinkComponent2[to]) continue;
            fromSinkEdge[to] = edge;
            toSinkEdge[to] = contractedGraph.get(ReverseEdge, edge);
        }
        for (const Vertex vertex : contractedGraph.vertices()) {
            if (!inSinkComponent1[vertex] && vertex != source) {
                for (const Edge edge : contractedGraph.edgesFrom(vertex)) {
                    const Vertex to = contractedGraph.get(ToVertex, edge);
                    if (!inSinkComponent1[to]) continue;
                    contractedGraph.get(Capacity, fromSourceEdge[to]) += contractedGraph.get(Capacity, edge);
                    const Edge reverseEdge = contractedGraph.get(ReverseEdge, edge);
                    contractedGraph.get(Capacity, toSourceEdge[to]) += contractedGraph.get(Capacity, reverseEdge);
                }
            } else if (inSinkComponent2[vertex] && vertex != sink) {
                for (const Edge edge : contractedGraph.edgesFrom(vertex)) {
                    const Vertex to = contractedGraph.get(ToVertex, edge);
                    if (inSinkComponent2[to]) continue;
                    contractedGraph.get(Capacity, fromSinkEdge[to]) += contractedGraph.get(Capacity, edge);
                    const Edge reverseEdge = contractedGraph.get(ReverseEdge, edge);
                    contractedGraph.get(Capacity, toSinkEdge[to]) += contractedGraph.get(Capacity, reverseEdge);
                }
            } else {
                if (vertex == source) newSource = Vertex(newToOldMapping.size());
                else if (vertex == sink) newSink = Vertex(newToOldMapping.size());
                newToOldMapping.emplace_back(newToOldVertex[vertex]);
            }
        }
        contractedGraph.deleteVertices([&](const Vertex vertex) {
            return (!inSinkComponent1[vertex] && vertex != source) || (inSinkComponent2[vertex] && vertex != sink);
        });
        return {*this, contractedGraph, newSource, newSink, newToOldMapping};
    }

public:
    GraphType graph;
    const Vertex source;
    const Vertex sink;
    double alpha;
    std::vector<FlowType> currentCapacity;
    std::vector<Vertex> newToOldVertex;
};