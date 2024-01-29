#pragma once

#include <queue>
#include <vector>

#include "../../DataStructures/Container/ExternalKHeap.h"
#include "../../DataStructures/Graph/Classes/GraphInterface.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/Types.h"

template<typename GRAPH>
class SimpleDijkstra {

private:
    struct VertexLabel : public ExternalKHeapElement {
        VertexLabel() : distance(INFTY) {}
        int distance;

        inline bool hasSmallerKey(const VertexLabel* other) const noexcept {
            return distance < other->distance;
        }
    };

public:
    SimpleDijkstra(const GRAPH& graph) :
        graph(graph),
        label(graph.numVertices()),
        parent(graph.numVertices(), noVertex),
        sourceVertex(noVertex) {
    }

    template<typename ATTRIBUTE>
    inline void run(const Vertex source, const ATTRIBUTE weight) noexcept {
        clear();
        sourceVertex = source;
        label[source].distance = 0;
        queue.update(&label[source]);

        while (!queue.empty()) {
            const VertexLabel* uLabel = queue.extractFront();
            const Vertex u = Vertex(uLabel - &(label[0]));
            const int uDistance = uLabel->distance;
            for (const Edge e : graph.edgesFrom(u)) {
                const Vertex v = graph.get(ToVertex, e);
                const int vDistance = uDistance + graph.get(weight, e);
                if (vDistance < label[v].distance) {
                    label[v].distance = vDistance;
                    parent[v] = u;
                    queue.update(&label[v]);
                }
            }
        }
    }

    inline int getDistance(const Vertex vertex) const noexcept {
        return label[vertex].distance;
    }

    inline bool reached(const Vertex vertex) const noexcept {
        return label[vertex].distance != INFTY;
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
        std::vector<VertexLabel>(label.size()).swap(label);
        std::vector<Vertex>(parent.size(), noVertex).swap(parent);
        queue.clear();
        sourceVertex = noVertex;
    }

    const GRAPH& graph;
    std::vector<VertexLabel> label;
    std::vector<Vertex> parent;
    ExternalKHeap<2, VertexLabel> queue;
    Vertex sourceVertex;
};
