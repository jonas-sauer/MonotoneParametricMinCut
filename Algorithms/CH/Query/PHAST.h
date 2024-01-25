#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../CH.h"
#include "../CHUtils.h"

#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/Vector/Permutation.h"

namespace CH {

template<bool REORDER_VERTICES, bool DEBUG = false, bool PATH_RETRIEVAL = false>
class PHAST {

public:
    using Graph = CHGraph;
    constexpr static bool ReorderVertices = REORDER_VERTICES;
    constexpr static bool Debug = DEBUG;
    constexpr static bool PathRetrieval = PATH_RETRIEVAL;
    using Type = PHAST<ReorderVertices, Debug, PathRetrieval>;

private:
    struct Label : public ExternalKHeapElement {
        Label(int* const distance) :
            ExternalKHeapElement(),
            distance(distance) {
        }

        inline int getDistance() const noexcept {
            return *distance;
        }

        inline void setDistance(const int newDistance) noexcept {
            *distance = newDistance;
        }

        inline bool hasSmallerKey(const Label* other) const noexcept {return getDistance() < other->getDistance();}

        int* const distance;
    };

    struct Data {
        Data(const Graph& forwardGraph, const Graph& backwardGraph, const Order& chOrder) :
            forward(forwardGraph),
            backward(backwardGraph),
            order(chOrder),
            fromVertex(backwardGraph.numEdges()) {
            if constexpr (REORDER_VERTICES) {
                forward.applyVertexOrder(order);
                forward.sortEdges(ToVertex);
                backward.applyVertexOrder(order);
                backward.sortEdges(ToVertex);

                for (const Vertex v : backward.vertices()) {
                    for (const Edge e : backward.edgesFrom(v)) {
                        fromVertex[e] = v;
                    }
                }
            }
        }

        inline Vertex internalVertex(const Vertex externalVertex) const noexcept {
            if constexpr (REORDER_VERTICES) {
                return externalVertex;
            } else {
                return Vertex(order[externalVertex]);
            }
        }

        Graph forward;
        Graph backward;
        Order order;
        std::vector<Vertex> fromVertex;
    };

public:
    PHAST(const CH& ch, const Order& order) :
        data(ch.forward, ch.backward, order),
        Q(ch.numVertices()),
        distance(ch.numVertices()),
        parent(PathRetrieval ? ch.numVertices() : 0, noVertex) {
        for (const Vertex vertex : ch.vertices()) {
            label.emplace_back(&distance[vertex]);
        }
    }

    inline void run(const Vertex source) noexcept {
        if constexpr (Debug) {
            std::cout << "Starting PHAST query" << std::endl;
            timer.restart();
        }
        clear();
        if constexpr (Debug) {
            std::cout << "Clear: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
            timer.restart();
        }
        upwardSearch(source);
        if constexpr (Debug) {
            std::cout << "Upward search: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
            timer.restart();
        }
        downwardSweep();
        if constexpr (Debug) {
            std::cout << "Downward sweep: " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
            timer.restart();
        }
    }

    inline bool reachable(const Vertex vertex) const noexcept {
        return getDistance(vertex) != INFTY;
    }

    inline int getDistance(const Vertex vertex) const noexcept {
        return distance[vertex];
    }

    template<bool T = PathRetrieval, typename = std::enable_if_t<T == PathRetrieval && T>>
    inline std::vector<Vertex> getPath(Vertex vertex) const noexcept {
        std::vector<Vertex> path{vertex};
        while (parent[vertex] != noVertex) {
            vertex = Vertex(parent[vertex]);
            path.emplace_back(vertex);
        }
        return path;
    }

private:
    inline void clear() noexcept {
        Vector::fill(distance, INFTY);
        if constexpr (PathRetrieval) Vector::fill(parent, int(noVertex));
    }

    inline void upwardSearch(const Vertex source) noexcept {
        distance[source] = 0;
        Q.update(&label[source]);

        while (!Q.empty()) {
            const Vertex u = Vertex(Q.extractFront() - &(label[0]));
            for (const Edge edge : data.forward.edgesFrom(u)) {
                const Vertex v = data.forward.get(ToVertex, edge);
                const int newDistance = distance[u] + data.forward.get(Weight, edge);
                if (newDistance < distance[v]) {
                    distance[v] = newDistance;
                    if constexpr (PathRetrieval) parent[v] = u;
                    Q.update(&label[v]);
                }
            }
        }
    }

    inline void downwardSweep() noexcept {
        if constexpr (ReorderVertices) {
            for (const Edge edge : data.backward.edges()) {
                const Vertex fromVertex = data.fromVertex[edge];
                const Vertex toVertex = data.backward.get(ToVertex, edge);
                const int newDistance = distance[toVertex] + data.backward.get(Weight, edge);
                const bool update = newDistance < distance[fromVertex];
                distance[fromVertex] = branchlessConditional(update, newDistance, distance[fromVertex]);
                if constexpr (PathRetrieval) parent[fromVertex] = branchlessConditional(update, toVertex, parent[fromVertex]);
            }
        } else {
            for (Vertex i(0); i < data.backward.numVertices(); i++) {
                const Vertex vertex = data.internalVertex(i);
                for (const Edge edge : data.backward.edgesFrom(vertex)) {
                    const Vertex toVertex = data.backward.get(ToVertex, edge);
                    const int newDistance = distance[toVertex] + data.backward.get(Weight, edge);
                    const bool update = newDistance < distance[vertex];
                    distance[vertex] = branchlessConditional(update, newDistance, distance[vertex]);
                    if constexpr (PathRetrieval) parent[vertex] = branchlessConditional(update, toVertex, parent[vertex]);
                }
            }
        }
    }

private:
    Data data;
    ExternalKHeap<2, Label> Q;
    std::vector<Label> label;
    std::vector<int> distance;
    std::vector<int> parent;
    Timer timer;
};

}
