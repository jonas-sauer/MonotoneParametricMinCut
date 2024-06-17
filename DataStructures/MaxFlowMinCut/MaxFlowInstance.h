#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Shell/Shell.h"

#include "../../Helpers/Console/ProgressBar.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/String/String.h"

#include "../../DataStructures/Graph/Graph.h"

//TODO: Capacities should be double
//TODO: Delete source-/sink-incident edges and maintain them implictly?
class MaxFlowInstance {
public:
    MaxFlowInstance() : source(noVertex), sink(noVertex) {}

    explicit MaxFlowInstance(const std::string& fileName) {
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

    template<bool VERBOSE = true>
    inline void fromDimacs(const std::string& fileName) noexcept {
        const std::string fileNameWithExtension = FileSystem::ensureExtension(fileName, ".max");
        if constexpr (VERBOSE) std::cout << "Reading DIMACS max-flow graph from: " << fileNameWithExtension << std::endl << std::flush;
        std::ifstream is(fileNameWithExtension);
        Assert(is.is_open(), "Cannot open file: " << fileNameWithExtension);
        ProgressBar bar(1);
        std::cout << "\r                     \r" << std::flush;

        FlowGraphEdgeList temp;
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
                    const int capacity = String::lexicalCast<int>(tokens[3]);
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
    StaticFlowGraph graph;
    Vertex source;
    Vertex sink;
};