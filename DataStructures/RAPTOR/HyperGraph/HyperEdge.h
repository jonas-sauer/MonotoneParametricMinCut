#pragma once

#include "../Data.h"

#include "../../../Algorithms/Dijkstra/Dijkstra.h"

namespace RAPTOR {

class HyperEdge {

private:
    HyperEdge(const Vertex vertex, const int weight) :
        dualVertex(vertex),
        weight(weight) {
    }

public:
    inline static HyperEdge FromVertexUnitWeight(const Vertex vertex, const Data&, const int weight = 1) {
        return HyperEdge(vertex, weight);
    }
    inline static HyperEdge FromEdgeTripWeight(const Vertex vertex, const Data& data, const int walkingDistance = 300) {
        int weight = 0;
        Dijkstra<TransferGraph, false> dijkstra(data.transferGraph, data.transferGraph[TravelTime]);
        dijkstra.run(vertex, noVertex, [&](const Vertex u) {
            if (data.isStop(u)) {
                weight += data.numberOfTripsContainingStop(StopId(u));
            }
        }, [&]() {
            return dijkstra.getDistance(dijkstra.getQFront()) > walkingDistance;
        });
        return HyperEdge(vertex, weight);
    }

public:
    inline int getWeight() const {
        return weight;
    }

public:
    Vertex dualVertex;
    int weight;
    std::vector<HyperVertexId> hyperVertices;

};

}
