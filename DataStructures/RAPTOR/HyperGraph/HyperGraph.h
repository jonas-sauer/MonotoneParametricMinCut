#pragma once

#include "HyperVertex.h"
#include "HyperEdge.h"

#include "../Data.h"

#include "../../HyperGraph/EdgeList.h"

namespace RAPTOR {

class HyperGraph : public ::HyperGraph::EdgeList<HyperVertex, HyperEdge> {

public:
    using Super = HyperGraph::EdgeList<HyperVertex, HyperEdge>;

public:
    HyperGraph(const Data& data, const int routeWeight = 0, const int edgeWeight = 1, const int vertexWeight = 1) :
        vertexIdByEdge(data.transferGraph.numEdges(), data.numberOfRoutes() + data.transferGraph.numEdges()) {
        Super& hyperGraph = *this;
        for (const Vertex vertex : data.transferGraph.vertices()) {
            if (vertexWeight < 0) {
                hyperGraph.addEdge(HyperEdge::FromEdgeTripWeight(vertex, data));
            } else {
                hyperGraph.addEdge(HyperEdge::FromVertexUnitWeight(vertex, data, vertexWeight));
            }
            if (data.isStop(vertex)) {
                for (const RouteSegment& route : data.routesContainingStop(StopId(vertex))) {
                    hyperGraph.edges.back().hyperVertices.emplace_back(HyperVertexId(route.routeId));
                }
            }
        }
        for (const RouteId route : data.routes()) {
            switch(routeWeight) {
            case 0:
                hyperGraph.addVertex(HyperVertex::FromRouteUnitWeight(route, data));
                break;
            case 1:
                hyperGraph.addVertex(HyperVertex::FromRouteTripWeight(route, data));
                break;
            case 2:
                hyperGraph.addVertex(HyperVertex::FromRouteStopWeight(route, data));
                break;
            case 3:
                hyperGraph.addVertex(HyperVertex::FromRouteProductWeight(route, data));
                break;
            case 4:
                hyperGraph.addVertex(HyperVertex::FromRouteSumWeight(route, data));
                break;
            default:
                AssertMsg(false, routeWeight << " is not a valid route weight!");
            }
            AssertMsg(hyperGraph.vertices.back().id == route, "Route " << route << " has invalid hyper vertex id: " << hyperGraph.vertices.back().id);
        }
        for (const Vertex vertex : data.transferGraph.vertices()) {
            for (const Edge edge : data.transferGraph.edgesFrom(vertex)) {
                if (hyperGraph.numVertices() >= vertexIdByEdge[edge]) continue;
                HyperVertex hyperVertex = hyperGraph.addVertex(HyperVertex::FromEdgeUnitWeight(edge, edgeWeight));
                vertexIdByEdge[edge] = hyperVertex.id;
                hyperGraph.edges[vertex].hyperVertices.emplace_back(hyperVertex.id);
                hyperGraph.edges[data.transferGraph.get(ToVertex, edge)].hyperVertices.emplace_back(hyperVertex.id);
                for (const Edge reverseEdge : data.transferGraph.edgesFrom(data.transferGraph.get(ToVertex, edge))) {
                    if (data.transferGraph.get(ToVertex, reverseEdge) != vertex) continue;
                    hyperVertex.dualEdges.emplace_back(reverseEdge);
                    vertexIdByEdge[reverseEdge] = hyperVertex.id;
                    break;
                }
            }
        }
    }

public:
    std::vector<size_t> vertexIdByEdge;

};

}
