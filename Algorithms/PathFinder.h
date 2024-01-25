#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "DepthFirstSearch.h"

#include "../Helpers/Types.h"

template<typename GRAPH>
class PathFinder {

private:
    class TraverseTreeEdge {

    public:
        TraverseTreeEdge(PathFinder<GRAPH>& super) :
            super(super),
            lastVertexIndex(-1) {
        }

        inline void operator() (const Edge edge, const Vertex) {
            Vertex vertex = super.graph.get(ToVertex, edge);
            if (super.index[vertex] >= 0) {
                super.numVerticesOnPaths++;
                if (lastVertexIndex < 0) super.numPaths++;
                super.index[vertex] = lastVertexIndex + 1;
            }
            lastVertexIndex = super.index[vertex];
        }

    private:
        PathFinder<GRAPH>& super;
        int lastVertexIndex;

    };

public:
    PathFinder(const GRAPH& graph) :
        graph(graph),
        numVerticesOnPaths(0),
        numPaths(0) {
    }

    inline void run() {
        initializeIndices();
        numVerticesOnPaths = 0;
        numPaths = 0;

        DepthFirstSearch<GRAPH, TraverseTreeEdge> dfs(graph, TraverseTreeEdge(*this));
        dfs.initialize();
        for (Vertex vertex : graph.vertices()) {
            if (!dfs.hasSettled(vertex) && index[vertex] < 0) {
                dfs.run(vertex);
            }
        }
    }

    inline bool isOnPath(Vertex vertex) {
        return index[vertex] >= 0;
    }

    inline int getIndex(Vertex vertex) {
        return index[vertex];
    }

    inline size_t getNumVerticesOnPaths() {
        return numVerticesOnPaths;
    }

    inline size_t getNumPaths() {
        return numPaths;
    }

private:
    inline void initializeIndices() {
        std::vector<int>(graph.numVertices(), 0).swap(index);
        std::vector<Vertex> firstNeighbor(graph.numVertices(), noVertex);
        std::vector<Vertex> secondNeighbor(graph.numVertices(), noVertex);

        for (Vertex v : graph.vertices()) {
            for (Edge edge : graph.edgesFrom(v)) {
                Vertex w = graph.get(ToVertex, edge);
                addNeighbor(w, firstNeighbor[v], secondNeighbor[v], index[v]);
                addNeighbor(v, firstNeighbor[w], secondNeighbor[w], index[w]);
            }
        }
        for (Vertex vertex : graph.vertices()) {
            if (secondNeighbor[vertex] == noVertex) index[vertex] = -1;
        }
    }

    inline void addNeighbor(const Vertex newNeighbor, Vertex& firstNeighbor, Vertex& secondNeighbor, int& index) {
        if (index < 0) return;
        if (firstNeighbor == noVertex || firstNeighbor == newNeighbor) {
            firstNeighbor = newNeighbor;
        } else if (secondNeighbor == noVertex || secondNeighbor == newNeighbor) {
            secondNeighbor = newNeighbor;
        } else {
            index = -1;
        }
    }

private:
    const GRAPH& graph;

    std::vector<int> index;
    size_t numVerticesOnPaths;
    size_t numPaths;

};
