#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../DataStructures/Graph/Graph.h"
#include "../../Helpers/Timer.h"

namespace HL {

template<typename GRAPH = TransferGraph, bool DEBUG = false>
class HLQuery {

public:
    using Graph = GRAPH;
    constexpr static bool Debug = DEBUG;
    using Type = HLQuery<Graph, Debug>;

public:
    HLQuery(const Graph& outHubs, const Graph& inHubs) :
        hubs {outHubs, inHubs},
        distance(INFTY) {
        hubs[FORWARD].sortEdges(ToVertex);
        hubs[BACKWARD].sortEdges(ToVertex);
    }

    inline void run(const Vertex from, const Vertex to) noexcept {
        if constexpr (Debug) {
            std::cout << "Starting HL query" << std::endl;
            std::cout << "   Source vertex: " << from << std::endl;
            std::cout << "   Target vertex: " << to << std::endl;
            timer.restart();
        }

        distance = INFTY;
        runInternal(from, to);

        if constexpr (Debug) std::cout << "   Time = " << String::musToString(timer.elapsedMicroseconds()) << std::endl;
    }

    inline int getDistance() const noexcept {
        return distance;
    }

    template<int DIRECTION>
    inline std::vector<Vertex> getHubs(const Vertex vertex) const noexcept {
        std::vector<Vertex> result;
        for (const Edge edge : hubs[DIRECTION].edgesFrom(vertex)) {
            result.emplace_back(hubs[DIRECTION].get(ToVertex, edge));
        }
        return result;
    }

private:
    inline void runInternal(const Vertex from, const Vertex to) noexcept {
        size_t j = 0;
        const Range<Edge>& outEdges = hubs[FORWARD].edgesFrom(from);
        const Range<Edge>& inEdges = hubs[BACKWARD].edgesFrom(to);

        for (size_t i = 0; i < outEdges.size(); i++) {
            while (hubs[FORWARD].get(ToVertex, outEdges[i]) > hubs[BACKWARD].get(ToVertex, inEdges[j])) {
                j++;
                if (j >= inEdges.size()) return;
            }
            if (hubs[FORWARD].get(ToVertex, outEdges[i]) == hubs[BACKWARD].get(ToVertex, inEdges[j])) {
                const int newDistance = hubs[FORWARD].get(TravelTime, outEdges[i]) + hubs[BACKWARD].get(TravelTime, inEdges[j]);
                distance = std::min(distance, newDistance);
            }
        }
    }

    TransferGraph hubs[2];
    int distance;
    Timer timer;

};

}
