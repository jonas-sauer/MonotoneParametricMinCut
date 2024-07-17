#pragma once

#include <queue>
#include <vector>

#include "../DataStructures/Graph/Classes/GraphInterface.h"
#include "../Helpers/Vector/Vector.h"

template<typename GRAPH>
class BreadthFirstSearch {

public:
    BreadthFirstSearch(const GRAPH& graph) :
        graph(graph),
        found(graph.numVertices(), false),
        parent(graph.numVertices(), noVertex),
        sourceVertex(noVertex) {
    }

    inline void run(const Vertex source) noexcept {
        clear();

        sourceVertex = source;
        found[source] = true;
        queue.push(source);

        while (!queue.empty()) {
            const Vertex u = queue.front();
            queue.pop();
            for (const Edge e : graph.edgesFrom(u)) {
                const Vertex v = graph.get(ToVertex, e);
                if (!found[v]) {
                    found[v] = true;
                    parent[v] = u;
                    queue.push(v);
                }
            }
        }
    }

    inline bool reached(const Vertex vertex) const noexcept {
        return found[vertex];
    }

    inline std::vector<Vertex> getPath(const Vertex target) const noexcept {
        Vertex v = target;
        std::vector<Vertex> path;
        path.emplace_back(v);
        while (v != sourceVertex) {
            v = parent[v];
            path.emplace_back(v);
        }
        return Vector::reverse(path);
    }

private:
    inline void clear() noexcept {
        std::vector<bool>(found.size(), false).swap(found);
        std::vector<Vertex>(parent.size(), noVertex).swap(parent);
        std::queue<Vertex>().swap(queue);
        sourceVertex = noVertex;
    }

    const GRAPH& graph;
    std::vector<bool> found;
    std::vector<Vertex> parent;
    std::queue<Vertex> queue;
    Vertex sourceVertex;
};
