#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Assert.h"
#include "../../Helpers/FileSystem/FileSystem.h"

namespace HyperGraph {

template<typename VERTEX, typename EDGE>
class EdgeList {

public:
    using HyperVertex = VERTEX;
    using HyperEdge = EDGE;
    using Type = EdgeList<HyperVertex, HyperEdge>;

public:
    inline size_t numVertices() const noexcept {
        return vertices.size();
    }

    inline size_t numEdges() const noexcept {
        return edges.size();
    }

    inline HyperVertex& addVertex(const HyperVertex& vertex) noexcept {
        vertices.emplace_back(vertex);
        vertices.back().id = HyperVertexId(vertices.size() - 1);
        return vertices.back();
    }

    inline HyperEdge& addEdge(const HyperEdge& edge) noexcept {
        edges.emplace_back(edge);
        return edges.back();
    }

    inline void toHMetis(const std::string& fileBaseName) const {
        std::ofstream hgrOs(FileSystem::ensureExtension(fileBaseName, ".hgr"));
        Assert(hgrOs);
        Assert(hgrOs.is_open());
        hgrOs << numEdges() << " " << numVertices() << " 11" << std::endl;
        for (const HyperEdge& edge : edges) {
            hgrOs << edge.getWeight();
            for (const HyperVertexId vertex : edge.hyperVertices) {
                hgrOs << " " << (vertex.value() + 1);
            }
            hgrOs << std::endl;
        }
        for (const HyperVertex& vertex : vertices) {
            hgrOs << vertex.getWeight() << std::endl;
        }
        hgrOs.close();
    }

    inline void printInfo() const noexcept {
        std::cout << "Hyper graph:" << std::endl;
        std::cout << "   Number of vertices:   " << std::setw(12) << String::prettyInt(vertices.size()) << std::endl;
        std::cout << "   Number of edges:      " << std::setw(12) << String::prettyInt(edges.size()) << std::endl;
    }

public:
    std::vector<HyperVertex> vertices;
    std::vector<HyperEdge> edges;

};

}
