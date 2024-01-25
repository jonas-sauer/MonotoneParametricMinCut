#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "CHQuery.h"

#include "../CHUtils.h"

#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/Vector/Permutation.h"

namespace CH {

template<typename GRAPH = CHGraph, bool STALL_ON_DEMAND = true, bool DEBUG = false>
class BidirectionalRPHAST {

public:
    using Graph = GRAPH;
    constexpr static bool StallOnDemand = STALL_ON_DEMAND;
    constexpr static bool Debug = DEBUG;
    using Type = BidirectionalRPHAST<Graph, StallOnDemand, Debug>;

    using BaseQuery = Query<Graph, StallOnDemand, false, true>;

private:
    struct Data {
        Data(const Graph& forwardGraph, const Graph& backwardGraph, const std::vector<Vertex>& chOrder, const Vertex::ValueType endOfPOIs) :
            vertexCount(0),
            forward(forwardGraph),
            backward(backwardGraph) {
            if (Debug) std::cout << "Computing restricted vertex set" << std::endl;
            std::vector<bool> isRequired(forward.numVertices(), false);
            std::vector<Vertex> stack;
            for (Vertex vertex = Vertex(0); vertex < endOfPOIs; vertex++) {
                isRequired[vertex] = true;
                stack.emplace_back(vertex);
                vertexCount++;
            }
            while (!stack.empty()) {
                const Vertex vertex = stack.back();
                stack.pop_back();
                for (const Edge edge : forward.edgesFrom(vertex)) {
                    const Vertex other = forward.get(ToVertex, edge);
                    if (isRequired[other]) continue;
                    isRequired[other] = true;
                    stack.emplace_back(other);
                    vertexCount++;
                }
                for (const Edge edge : backward.edgesFrom(vertex)) {
                    const Vertex other = backward.get(ToVertex, edge);
                    if (isRequired[other]) continue;
                    isRequired[other] = true;
                    stack.emplace_back(other);
                    vertexCount++;
                }
            }
            if (Debug) std::cout << "Reordering vertices" << std::endl;
            vertexOrder = Order(chOrder);
            vertexPermutation = Permutation(Construct::Invert, vertexOrder);
            sort(vertexOrder, [&](const size_t a, const size_t b) {
                return (isRequired[a] > isRequired[b]) || ((isRequired[a] == isRequired[b]) && (vertexPermutation[a] > vertexPermutation[b]));
            });
            forward.applyVertexOrder(vertexOrder);
            forward.sortEdges(ToVertex);
            backward.applyVertexOrder(vertexOrder);
            backward.sortEdges(ToVertex);
            vertexPermutation = Permutation(Construct::Invert, vertexOrder);
            if (Debug) std::cout << "Computing restricted graph" << std::endl;
            CHConstructionGraph temp;
            temp.addVertices(vertexCount * 3);
            for (Vertex vertex = Vertex(0); vertex < vertexCount; vertex++) {
                for (const Edge edge : forward.edgesFrom(vertex)) {
                    const Vertex other = forward.get(ToVertex, edge);
                    if (other >= vertexCount) continue;
                    const Edge reverse = backward.findEdge(vertex, other);
                    if ((!backward.isEdge(reverse)) || (forward.get(Weight, edge) != backward.get(Weight, reverse))) {
                        temp.addEdge(Vertex(other * 3), vertex).set(Weight, forward.get(Weight, edge));
                    } else {
                        temp.addEdge(Vertex((other * 3) + 1), vertex).set(Weight, forward.get(Weight, edge));
                    }
                }
                for (const Edge edge : backward.edgesFrom(vertex)) {
                    const Vertex other = backward.get(ToVertex, edge);
                    if (other >= vertexCount) continue;
                    const Edge reverse = backward.findEdge(vertex, other);
                    if ((!forward.isEdge(reverse)) || (backward.get(Weight, edge) != forward.get(Weight, reverse))) {
                        temp.addEdge(Vertex((other * 3) + 2), vertex).set(Weight, backward.get(Weight, edge));
                    }
                }
            }
            ::Graph::move(std::move(temp), restricted);
            restricted.sortEdges(ToVertex);
            if (Debug) ::Graph::printInfo(restricted);
            if (Debug) restricted.printAnalysis();
        }
        inline Vertex internalVertex(const Vertex external) const noexcept {
            return Vertex(vertexPermutation[external]);
        }
        inline Vertex externalVertex(const Vertex internal) const noexcept {
            return Vertex(vertexOrder[internal]);
        }
        size_t vertexCount;
        Order vertexOrder;
        Permutation vertexPermutation;
        Graph forward;
        Graph backward;
        CHGraph restricted;
    };

public:
    BidirectionalRPHAST(const Graph& forward, const Graph& backward, const std::vector<Vertex>& chOrder, const Vertex::ValueType endOfPOIs) :
        data(forward, backward, chOrder, endOfPOIs),
        baseQuery(data.forward, data.backward, data.forward.numVertices(), Weight),
        endOfPOIs(endOfPOIs),
        distance {std::vector<int>(data.vertexCount, INFTY), std::vector<int>(data.vertexCount, INFTY)},
        reachedPOIs {std::vector<Vertex>(), std::vector<Vertex>()} {
    }

    BidirectionalRPHAST(const CH& ch, const int direction = FORWARD, const Vertex::ValueType endOfPOIs = 0) :
        BidirectionalRPHAST(ch.getGraph(direction), ch.getGraph(!direction), getOrder(ch), endOfPOIs) {
    }

    inline void run(const Vertex from, const Vertex to) noexcept {
        if (Debug) std::cout << "Starting bidirectional RPHAST query" << std::endl;
        if (Debug) timer.restart();

        clear<FORWARD>();
        clear<BACKWARD>();

        baseQuery.run(data.internalVertex(from), data.internalVertex(to));

        collectPOIs<FORWARD>();
        collectPOIs<BACKWARD>();

        const int maxDistance = baseQuery.getDistance();
        for (Vertex vertex = Vertex(0); vertex < data.vertexCount; vertex++) {
            const bool relaxForward = distance[FORWARD][vertex] < maxDistance;
            const bool relaxBackward = distance[BACKWARD][vertex] < maxDistance;

            if (relaxForward) {
                if (relaxBackward) {
                    relaxEdges<true, true>(vertex);
                } else {
                    relaxEdges<true, false>(vertex);
                }
            } else {
                if (relaxBackward) {
                    relaxEdges<false, true>(vertex);
                }
            }
        }

        if (Debug) std::cout << "   Time = " << String::msToString(timer.elapsedMilliseconds()) << std::endl;
    }

    inline bool reachable() const noexcept {
        return baseQuery.reachable();
    }

    inline bool visited(const Vertex vertex) const noexcept {
        return baseQuery.visited(data.internalVertex(vertex));
    }

    inline int getDistance(const Vertex = noVertex) const noexcept {
        return baseQuery.getDistance();
    }

    inline int getForwardDistance(const Vertex vertex) noexcept {
        return distance[FORWARD][data.internalVertex(vertex)];
    }

    inline int getBackwardDistance(const Vertex vertex) noexcept {
        return distance[BACKWARD][data.internalVertex(vertex)];
    }

    inline const std::vector<Vertex>& getForwardPOIs() const noexcept {
        return reachedPOIs[FORWARD];
    }

    inline const std::vector<Vertex>& getBackwardPOIs() const noexcept {
        return reachedPOIs[BACKWARD];
    }

    inline int getSettleCount() const noexcept {
        return baseQuery.getSettleCount();
    }

    inline int getStallCount() const noexcept {
        return baseQuery.getStallCount();
    }

private:
    template<int DIRECTION>
    inline void clear() noexcept {
        Vector::fill(distance[DIRECTION], INFTY);
        reachedPOIs[DIRECTION].clear();
    }

    template<int DIRECTION>
    inline void collectPOIs() noexcept {
        const int maxDistance = baseQuery.getDistance();
        for (const Vertex vertex : baseQuery.template getPOIs<DIRECTION>()) {
            if (baseQuery.template getDistanceToPOI<DIRECTION>(vertex) > maxDistance) break;
            if (vertex >= distance[DIRECTION].size()) continue;
            distance[DIRECTION][vertex] = baseQuery.template getDistanceToPOI<DIRECTION>(vertex);
        }
    }

    template<bool RELAX_FORWARD, bool RELAX_BACKWARD>
    inline void relaxEdges(const Vertex vertex) noexcept {
        if constexpr (!RELAX_FORWARD && !RELAX_BACKWARD) {
            suppressUnusedParameterWarning(vertex);
            return;
        }

        const Vertex external = data.externalVertex(vertex);

        if constexpr (RELAX_BACKWARD) {
            if (external < endOfPOIs) {
                reachedPOIs[BACKWARD].emplace_back(external);
            }
            for (const Edge edge : data.restricted.edgesFrom(Vertex(vertex * 3))) {
                updateDistance<BACKWARD>(vertex, edge);
            }
        }
        for (const Edge edge : data.restricted.edgesFrom(Vertex((vertex * 3) + 1))) {
            if constexpr (RELAX_BACKWARD) updateDistance<BACKWARD>(vertex, edge);
            if constexpr (RELAX_FORWARD) updateDistance<FORWARD>(vertex, edge);
        }
        if constexpr (RELAX_FORWARD) {
            if (external < endOfPOIs) {
                reachedPOIs[FORWARD].emplace_back(external);
            }
            for (const Edge edge : data.restricted.edgesFrom(Vertex((vertex * 3) + 2))) {
                updateDistance<FORWARD>(vertex, edge);
            }
        }
    }

    template<int DIRECTION>
    inline void updateDistance(const Vertex fromVertex, const Edge edge) noexcept {
        const int newDistance = distance[DIRECTION][fromVertex] + data.restricted.get(Weight, edge);
        const Vertex toVertex = data.restricted.get(ToVertex, edge);
        distance[DIRECTION][toVertex] = std::min(distance[DIRECTION][toVertex], newDistance);
    }

private:
    Data data;
    BaseQuery baseQuery;
    Vertex endOfPOIs;

    std::vector<int> distance[2];
    std::vector<Vertex> reachedPOIs[2];

    Timer timer;

};

}
