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
                    edgeCount = String::lexicalCast<size_t>(tokens[3]);
                    temp.addVertices(vertexCount);
                    if constexpr (VERBOSE) bar.init(edgeCount);
                }
            } else if (source == noVertex || sink == noVertex) {
                if (tokens.size() != 3 || tokens[0] != "n" || (tokens[2] != "s" && tokens[2] != "t")) {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else if (tokens[2] == "s") {
                    source = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    if (!temp.isVertex(source)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                } else {
                    sink = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
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
                    if (!temp.isVertex(from)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    } else if (!temp.isVertex(to)) {
                        std::cout << "ERROR, " << tokens[2] << " does not name a vertex!" << std::endl;
                        break;
                    }
                    temp.addEdge(from, to).set(Capacity, capacity);
                    if constexpr (VERBOSE) bar++;
                }
            }
        }
        is.close();
        if constexpr (VERBOSE) std::cout << std::endl;
        if (temp.numEdges() != edgeCount) {
            std::cout << "WARNING, found " << temp.numEdges() << " edges, but " << edgeCount << " edges were declared." << std::endl;
        }

        for (const auto [edge, from] : temp.edgesWithFromVertex()) {
            if (temp.hasReverseEdge(edge)) continue;
            temp.addReverseEdge(edge).set(Capacity, 0);
        }
        Graph::move(std::move(temp), graph);
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

    ParametricMaxFlowInstance() : source(noVertex), sink(noVertex), alphaMin(0), alphaMax(0) {}

    explicit ParametricMaxFlowInstance(const std::string& fileName) {
        deserialize(fileName);
    }

    template<typename T>
    explicit ParametricMaxFlowInstance(const StaticMaxFlowInstance<T>& staticInstance, const int steps = 0) : source(staticInstance.source), sink(staticInstance.sink), alphaMin(0), alphaMax(steps) {
        Graph::copy(staticInstance.graph, graph);
        for (const Edge e : graph.edges()) {
            graph.set(Capacity, e, FlowFunction(staticInstance.graph.get(Capacity, e)));
        }
        if (steps > 0) {
            for (const Edge edge : graph.edgesFrom(source)) {
                const FlowType capacity = staticInstance.graph.get(Capacity, edge);
                graph.set(Capacity, edge, FlowFunction(0, capacity / steps));
            }
            for (const Edge edge : graph.edgesFrom(sink)) {
                const Edge reverseEdge = graph.get(ReverseEdge, edge);
                const FlowType capacity = staticInstance.graph.get(Capacity, reverseEdge);
                graph.set(Capacity, reverseEdge, FlowFunction(capacity, -capacity / steps));
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
                    edgeCount = String::lexicalCast<size_t>(tokens[3]);
                    temp.addVertices(vertexCount);
                    if constexpr (VERBOSE) bar.init(edgeCount);
                }
            } else if (source == noVertex || sink == noVertex) {
                if (tokens.size() != 3 || tokens[0] != "n" || (tokens[2] != "s" && tokens[2] != "t")) {
                    std::cout << "ERROR, invalid DIMACS .max-file header: " << line << std::endl;
                    break;
                } else if (tokens[2] == "s") {
                    source = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
                    if (!temp.isVertex(source)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    }
                } else {
                    sink = Vertex(String::lexicalCast<size_t>(tokens[1]) - 1);
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
                    if (!temp.isVertex(from)) {
                        std::cout << "ERROR, " << tokens[1] << " does not name a vertex!" << std::endl;
                        break;
                    } else if (!temp.isVertex(to)) {
                        std::cout << "ERROR, " << tokens[2] << " does not name a vertex!" << std::endl;
                        break;
                    }
                    temp.addEdge(from, to).set(Capacity, FlowFunction(capacityA, capacityB));
                    if constexpr (VERBOSE) bar++;
                }
            }
        }
        is.close();
        if constexpr (VERBOSE) std::cout << std::endl;
        if (temp.numEdges() != edgeCount) {
            std::cout << "WARNING, found " << temp.numEdges() << " edges, but " << edgeCount << " edges were declared." << std::endl;
        }

        for (const auto [edge, from] : temp.edgesWithFromVertex()) {
            if (temp.hasReverseEdge(edge)) continue;
            temp.addReverseEdge(edge).set(Capacity, FlowFunction(0));
        }
        Graph::move(std::move(temp), graph);
    }

public:
    GraphType graph;
    Vertex source;
    Vertex sink;
    double alphaMin;
    double alphaMax;
};

template<Meta::Derived<pmf::flowFunction> FLOW_FUNCTION>
class ParametricToStaticMaxFlowInstanceWrapper {
public:
    using FlowFunction = FLOW_FUNCTION;
    using InstanceType = ParametricMaxFlowInstance<FlowFunction>;
    using GraphType = InstanceType::GraphType;
    using FlowType = InstanceType::FlowType;

public:
    ParametricToStaticMaxFlowInstanceWrapper(const InstanceType& instance, const double alpha) :
        instance(instance),
        graph(instance.graph),
        source(instance.source),
        sink(instance.sink),
        alpha(alpha),
        sourceDiff(graph.outDegree(source), 0),
        sinkDiff(graph.outDegree(sink), 0) {
        for (const Edge edge : instance.graph.edges()) {
            currentCapacity[edge] = instance.getCapacity(edge, alpha);
        }
    }

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
        Assert(newAlpha > alpha, "Parameter did not increase!");
        Assert(newAlpha <= instance.alphaMax, "Parameter out of range!");
        alpha = newAlpha;
        Edge edgeFromSource = graph.beginEdgeFrom(source);
        for (size_t i = 0; i < sourceDiff.size(); i++, edgeFromSource++) {
            const int newCapacity = instance.getCapacity(edgeFromSource).eval(alpha);
            sourceDiff[i] = newCapacity - currentCapacity[edgeFromSource];
            currentCapacity[edgeFromSource] = newCapacity;
        }
        Edge edgeFromSink = graph.beginEdgeFrom(sink);
        for (size_t i = 0; i < sinkDiff.size(); i++, edgeFromSink++) {
            const Edge edgeToSink = graph.get(ReverseEdge, edgeFromSink);
            const int newCapacity = instance.getCapacity(edgeToSink).eval(alpha);
            Assert(newCapacity >= 0, "Negative capacity!");
            sinkDiff[i] = newCapacity - currentCapacity[edgeToSink];
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