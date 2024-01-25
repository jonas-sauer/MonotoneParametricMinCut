#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "HLQuery.h"

#include "../../DataStructures/Container/Set.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../Helpers/Console/Progress.h"
#include "../../Helpers/Timer.h"

namespace HL {

template<typename GRAPH = TransferGraph, bool DEBUG = false>
class BucketHL {

public:
    using Graph = GRAPH;
    constexpr static bool Debug = DEBUG;
    using Type = BucketHL<Graph, Debug>;
    using BaseQuery = HLQuery<Graph, Debug>;

public:
    BucketHL(const Graph& outHubs, const Graph& inHubs, const Vertex::ValueType endOfPOIs) :
        baseQuery(outHubs, inHubs),
        hubs {outHubs, inHubs},
        bucketGraph {TransferGraph(), TransferGraph()},
        distance {std::vector<int>(outHubs.numVertices(), INFTY), std::vector<int>(inHubs.numVertices(), INFTY)},
        root {noVertex, noVertex},
        endOfPOIs(endOfPOIs),
        reachedPOIs {IndexedSet<false, Vertex>(endOfPOIs), IndexedSet<false, Vertex>(endOfPOIs)} {
        hubs[FORWARD].sortEdges(TravelTime);
        hubs[BACKWARD].sortEdges(TravelTime);
        buildBucketGraph<FORWARD, BACKWARD>();
        buildBucketGraph<BACKWARD, FORWARD>();
    }

    inline void run(const Vertex from, const Vertex to) noexcept {
        if (root[FORWARD] == from && root[BACKWARD] == to) return;
        if constexpr (Debug) {
            std::cout << "Starting bucket query" << std::endl;
            timer.restart();
        }

        root[FORWARD] = from;
        root[BACKWARD] = to;
        baseQuery.run(from, to);
        run<FORWARD>();
        run<BACKWARD>();

        if constexpr (Debug) std::cout << "   Time = " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
    }

    inline int getDistance(const Vertex = noVertex) const noexcept {
        return baseQuery.getDistance();
    }

    inline const std::vector<int>& getForwardDistance() const noexcept {
        return distance[FORWARD];
    }

    inline int getForwardDistance(const Vertex vertex) const noexcept {
        return distance[FORWARD][vertex];
    }

    inline const std::vector<int>& getBackwardDistance() const noexcept {
        return distance[BACKWARD];
    }

    inline int getBackwardDistance(const Vertex vertex) const noexcept {
        return distance[BACKWARD][vertex];
    }

    inline const IndexedSet<false, Vertex>& getForwardPOIs() const noexcept {
        return reachedPOIs[FORWARD];
    }

    inline const IndexedSet<false, Vertex>& getBackwardPOIs() const noexcept {
        return reachedPOIs[BACKWARD];
    }

private:
    template<int I, int J>
    inline void buildBucketGraph() noexcept {
        if constexpr (Debug) std::cout << "Building " << ((I == FORWARD) ? ("forward") : ("backward")) << " bucket graph" << std::endl;
        TransferEdgeList temp;
        temp.addVertices(distance[I].size());
        Progress progress(endOfPOIs, Debug);
        for (Vertex vertex = Vertex(0); vertex < endOfPOIs; vertex++) {
            for (const Edge edge : hubs[J].edgesFrom(vertex)) {
                const Vertex hub = hubs[J].get(ToVertex, edge);
                temp.addEdge(hub, vertex).set(TravelTime, hubs[J].get(TravelTime, edge));
            }
            progress++;
        }
        ::Graph::move(std::move(temp), bucketGraph[I]);
        bucketGraph[I].sortEdges(TravelTime);
        if constexpr (Debug) {
            std::cout << std::endl;
            ::Graph::printInfo(bucketGraph[I]);
            bucketGraph[I].printAnalysis();
        }
    }

    template<int I>
    inline void run() noexcept {
        clear<I>();
        scanHubs<I>();
        scanBuckets<I>();
    }

    template<int I>
    inline void clear() noexcept {
        for (const Vertex vertex : updatedHubs[I]) {
            distance[I][vertex] = INFTY;
        }
        for (const Vertex vertex : reachedPOIs[I]) {
            distance[I][vertex] = INFTY;
        }
        updatedHubs[I].clear();
        reachedPOIs[I].clear();
    }

    template<int I>
    inline void scanHubs() noexcept {
        distance[I][root[I]] = 0;
        updatedHubs[I].emplace_back(root[I]);
        for (const Edge edge : hubs[I].edgesFrom(root[I])) {
            const Vertex hub = hubs[I].get(ToVertex, edge);
            const int newDistance = hubs[I].get(TravelTime, edge);
            if (newDistance > baseQuery.getDistance()) break;
            distance[I][hub] = newDistance;
            updatedHubs[I].emplace_back(hub);
        }
    }

    template<int I>
    inline void scanBuckets() noexcept {
        for (const Vertex hub : updatedHubs[I]) {
            if (hub < endOfPOIs) reachedPOIs[I].insert(hub);
            for (const Edge edge : bucketGraph[I].edgesFrom(hub)) {
                const Vertex poi = bucketGraph[I].get(ToVertex, edge);
                const int newDistance = distance[I][hub] + bucketGraph[I].get(TravelTime, edge);
                if (newDistance > baseQuery.getDistance()) break;
                if (newDistance < distance[I][poi]) {
                    distance[I][poi] = newDistance;
                    reachedPOIs[I].insert(poi);
                }
            }
        }
    }

private:
    BaseQuery baseQuery;

    TransferGraph hubs[2];
    TransferGraph bucketGraph[2];
    std::vector<int> distance[2];

    Vertex root[2];

    Vertex endOfPOIs;
    std::vector<Vertex> updatedHubs[2];
    IndexedSet<false, Vertex> reachedPOIs[2];

    Timer timer;

};

}
