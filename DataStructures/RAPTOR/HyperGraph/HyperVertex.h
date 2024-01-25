#pragma once

#include "../Data.h"

namespace RAPTOR {

class HyperVertex {

private:
    HyperVertex(const RouteId route, const int weight) :
        id(noHyperVertex),
        isRouteVertex(true),
        weight(weight),
        route(route) {
    }

    HyperVertex(const Edge edge, const int weight) :
        id(noHyperVertex),
        isRouteVertex(false),
        weight(weight),
        dualEdges(1, edge) {
    }

public:
    inline static HyperVertex FromRouteUnitWeight(const RouteId route, const Data& data) {
        return HyperVertex(VertexIdFromRoute(route, data), 1);
    }
    inline static HyperVertex FromRouteTripWeight(const RouteId route, const Data& data) {
        return HyperVertex(VertexIdFromRoute(route, data), data.numberOfTripsInRoute(route));
    }
    inline static HyperVertex FromRouteStopWeight(const RouteId route, const Data& data) {
        return HyperVertex(VertexIdFromRoute(route, data), data.numberOfStopsInRoute(route));
    }
    inline static HyperVertex FromRouteProductWeight(const RouteId route, const Data& data) {
        return HyperVertex(VertexIdFromRoute(route, data), data.numberOfTripsInRoute(route) * data.numberOfStopsInRoute(route));
    }
    inline static HyperVertex FromRouteSumWeight(const RouteId route, const Data& data) {
        return HyperVertex(VertexIdFromRoute(route, data), data.numberOfTripsInRoute(route) + data.numberOfStopsInRoute(route));
    }

    inline static HyperVertex FromEdgeUnitWeight(const Edge edge, const int weight = 1) {
        return HyperVertex(edge, weight);
    }

    inline static RouteId VertexIdFromRoute(const RouteId route, const Data&) {
        return route;
    }

public:
    // inline void setRoute(const RouteId route) noexcept {
    //     this->route = route;
    //     isRouteVertex = true;
    // }

    // inline void setEdge(const Edge) noexcept {
    //     //this->edges.emplace_back(edge);
    //     isRouteVertex = false;
    // }

    inline int getWeight() const {
        return weight;
    }

public:
    HyperVertexId id;
    bool isRouteVertex;
    int weight;
    RouteId route;
    std::vector<Edge> dualEdges;

};

}
