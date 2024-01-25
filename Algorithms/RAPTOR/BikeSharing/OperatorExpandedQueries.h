#pragma once

#include "OperatorDependentQueries.h"

#include "../Profiler.h"
#include "../ULTRARAPTOR.h"
#include "../DijkstraRAPTOR.h"
#include "../InitialTransfers.h"

#include "../../../DataStructures/BikeSharing/ExtendedNetwork.h"

namespace RAPTOR {

class BikeSharingProfiler : public NoProfiler {

public:
    BikeSharingProfiler() :
        stopCount(0),
        vertexCount(0),
        routeCount(0),
        roundCount(0) {
    }

public:
    inline void start() noexcept {
        roundCount++;
    }

    inline void newRound() noexcept {
        roundCount++;
    }

    inline void scanRoute(const RouteId) noexcept {
        routeCount++;
    }

    inline void updateStopByRoute(const StopId, const int) noexcept {
        stopCount++;
    }

    inline void updateStopByTransfer(const StopId, const int) noexcept {
        stopCount++;
    }

    inline void settleVertex(const Vertex) noexcept {
        vertexCount++;
    }

    inline void printData(const double f = 1.0) noexcept {
        std::cout << "Number of scanned routes: " << String::prettyDouble(routeCount / f, 0) << std::endl;
        std::cout << "Number of settled vertices: " << String::prettyDouble(vertexCount / f, 0) << std::endl;
        std::cout << "Number of rounds: " << String::prettyDouble(roundCount / f, 2) << std::endl;
        stopCount = 0;
        vertexCount = 0;
        routeCount = 0;
        roundCount = 0;
    }

public:
    size_t stopCount;
    size_t vertexCount;
    size_t routeCount;
    size_t roundCount;

};

template<bool DEBUG = false>
class OperatorExpandedRAPTOR {

public:
    static constexpr bool Debug = DEBUG;
    using Type = OperatorExpandedRAPTOR<Debug>;

public:
    OperatorExpandedRAPTOR(const BikeSharing::ExtendedNetwork& extendedNetwork) :
        extendedNetwork(extendedNetwork),
        raptor(extendedNetwork.data, extendedNetwork.coreCH) {
    }

public:
    inline void run(const Vertex source, const int departureTime, const Vertex target) noexcept {
        raptor.run(extendedNetwork.oldToNewVertex[source], departureTime, extendedNetwork.oldToNewVertex[target]);
    }

    inline int getEarliestArrivalTime() const noexcept {
        return raptor.getEarliestArrivalTime();
    }

    inline int getEarliestArrivalNumberOfTrips() const noexcept {
        return raptor.getEarliestArrivalNumberOfTrips();
    }

    inline void debug(const double f = 1.0) noexcept {
        raptor.getProfiler().printData(f);
    }

private:
    const BikeSharing::ExtendedNetwork& extendedNetwork;

    DijkstraRAPTOR<CoreCHInitialTransfers, BikeSharingProfiler> raptor;

};

template<bool DEBUG = false>
class OperatorExpandedULTRA {

public:
    static constexpr bool Debug = DEBUG;
    using Type = OperatorExpandedULTRA<Debug>;

public:
    OperatorExpandedULTRA(const BikeSharing::ExtendedNetwork& extendedNetwork) :
        extendedNetwork(extendedNetwork),
        ultra(extendedNetwork.data, extendedNetwork.ch) {
    }

public:
    inline void run(const Vertex source, const int departureTime, const Vertex target) noexcept {
        ultra.run(extendedNetwork.oldToNewVertex[source], departureTime, extendedNetwork.oldToNewVertex[target]);
    }

    inline int getEarliestArrivalTime() const noexcept {
        return ultra.getEarliestArrivalTime();
    }

    inline int getEarliestArrivalNumberOfTrips() const noexcept {
        return ultra.getEarliestArrivalNumberOfTrips();
    }

    inline void debug(const double f = 1.0) noexcept {
        ultra.getProfiler().printData(f);
    }

    inline void printPath() const noexcept {
        std::vector<Vertex> newToOldVertex(extendedNetwork.graph.numVertices(), noVertex);
        for (size_t i = 0; i < extendedNetwork.oldToNewVertex.size(); i++) {
            newToOldVertex[extendedNetwork.oldToNewVertex[i]] = Vertex(i);
        }
        RAPTOR::Journey journey = ultra.getJourneys().back();
        for (size_t i = journey.size() - 1; i < journey.size(); i --) {
            std::cout << "from: " << newToOldVertex[journey[i].from] << ", to: " << newToOldVertex[journey[i].to] << ", departureTime: " << journey[i].departureTime << ", arrivalTime: " << journey[i].arrivalTime << ", " << (journey[i].usesRoute ? journey[i].routeId : noRouteId) << std::endl;
        }
        std::vector<Vertex> path = ultra.getPath();
        for (size_t i = path.size() - 1; i < path.size(); i --) {
            std::cout << newToOldVertex[path[i]] << std::endl;
        }
    }

private:
    const BikeSharing::ExtendedNetwork& extendedNetwork;

    ULTRARAPTOR<BikeSharingProfiler> ultra;

};

template<bool DEBUG = false>
class OperatorExpandedOperatorDependentRAPTOR {

public:
    static constexpr bool Debug = DEBUG;
    using Type = OperatorExpandedULTRA<Debug>;

public:
    OperatorExpandedOperatorDependentRAPTOR(const BikeSharing::ExtendedNetwork& extendedNetwork) :
        extendedNetwork(extendedNetwork),
        extendedData(extendedNetwork.data, extendedNetwork.graph, extendedNetwork.coreCH),
        raptor(extendedData) {
    }

public:
    inline void run(const Vertex source, const int departureTime, const Vertex target) noexcept {
        raptor.run(extendedNetwork.oldToNewVertex[source], departureTime, extendedNetwork.oldToNewVertex[target]);
    }

    inline int getEarliestArrivalTime() const noexcept {
        return raptor.getEarliestArrivalTime();
    }

    inline int getEarliestArrivalNumberOfTrips() const noexcept {
        return raptor.getEarliestArrivalNumberOfTrips();
    }

    inline void debug(const double f = 1.0) noexcept {
        raptor.debug(f);
    }

    inline void printPath() const noexcept {
        raptor.printPath();
    }

private:
    const BikeSharing::ExtendedNetwork& extendedNetwork;
    const BikeSharing::Data extendedData;

    OperatorDependentRAPTOR<false, Debug, false> raptor;

};

}
