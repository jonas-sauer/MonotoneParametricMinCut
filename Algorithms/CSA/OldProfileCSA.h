#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>

#include "../Dijkstra/Dijkstra.h"

#include "../../Helpers/Assert.h"
#include "../../Helpers/Console/Progress.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"
#include "../../DataStructures/Container/ExternalKHeap.h"
#include "../../DataStructures/Profile/Profile.h"
#include "../../DataStructures/Profile/Entities/ArrivalArrayEntry.h"
#include "../../DataStructures/Profile/Entities/ArrivalAndParentArrayEntry.h"
#include "../../DataStructures/Profile/Entities/EarliestArrivalEntry.h"

namespace CSA {

template<size_t MAX_NUMBER_OF_TRIPS, bool ALLOW_MULTI_HOP_TRANSFERS = true, bool TARGET_PRUNING = false, bool DEBUG = false>
class OldProfileCSA {

public:
    constexpr static size_t MaxNumberOfTrips = MAX_NUMBER_OF_TRIPS;
    constexpr static bool AllowMultiHopTransfers = ALLOW_MULTI_HOP_TRANSFERS;
    constexpr static bool TargetPruning = TARGET_PRUNING;
    constexpr static bool Debug = DEBUG;
    using Type = OldProfileCSA<MaxNumberOfTrips, AllowMultiHopTransfers, TargetPruning, Debug>;
    using ArrivalEntry = std::conditional_t<AllowMultiHopTransfers, ArrivalAndParentArrayEntry<MaxNumberOfTrips, int>, ArrivalArrayEntry<MaxNumberOfTrips, int>>;
    using Label = Profile<ArrivalEntry>;

public:
    OldProfileCSA(const Data& data, const TransferGraph& reverseGraph) :
        data(data),
        reverseGraph(reverseGraph),
        sourceStop(noStop),
        targetStop(noStop),
        minDepartureTime(never),
        maxDepartureTime(never),
        arrivalByTrip(data.numberOfTrips()),
        profileWithoutInitialTransfers(data.numberOfStops()),
        labels(reverseGraph.numVertices()),
        queue(reverseGraph.numVertices()),
        dijkstra(reverseGraph, reverseGraph[TravelTime]) {
        if (AllowMultiHopTransfers) buildMinChangeTimeGraph();
    }

    OldProfileCSA(const Data& data) : OldProfileCSA(data, data.transferGraph) {}

    inline void run(const StopId source, const StopId target, const int minTime = 0, const int maxTime = 60 * 60 * 24) noexcept {
        clear();
        sourceStop = source;
        targetStop = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        initialize();
        Progress progress(data.numberOfConnections(), Debug);
        for (int i = data.numberOfConnections() - 1; i >= 0; i--) {
            const Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;

            if (AllowMultiHopTransfers) runDijkstra(connection.departureTime);

            Label& departureProfile = profileWithoutInitialTransfers[connection.departureStopId];
            Label& departureLabel = labels[connection.departureStopId];
            const Label& arrivalLabel = labels[connection.arrivalStopId];
            Assert(departureProfile.earliestDeparture() >= connection.departureTime, "Connections are scanned out of order (" << departureProfile.earliestDeparture() << " before " << connection.departureTime << ")!\nTry: data.sortConnectionsAscendingByDepartureTime();");

            ArrivalEntry arrival = arrivalByTrip[connection.tripId];
            if (arrivalLabel.isConnectedToTarget()) arrival.minimize(connection.arrivalTime + arrivalLabel.getWalkingTime());
            arrival.minimizeShifted(arrivalLabel.evaluateLinear(connection.arrivalTime));
            arrivalByTrip[connection.tripId] = arrival;

            if (!(arrival < departureProfile.earliestArrival())) continue;
            departureProfile.append(connection.departureTime, arrival);
            const int newDepartureTime = connection.departureTime - data.stopData[connection.departureStopId].minTransferTime;
            if (TargetPruning) {
                if (labels[sourceStop].dominates(newDepartureTime, arrival)) continue;
            }
            arrival.setParent(connection.departureStopId);
            if (departureLabel.insertIfImprovement(newDepartureTime, arrival)) {
                if (AllowMultiHopTransfers) queue.update(&departureLabel);
            }
            settle(connection.departureStopId, connection.departureTime, arrival);

            if (Debug) progress++;
        }
        if (AllowMultiHopTransfers) runDijkstra(minDepartureTime);
        mergeProfiles();
        if (Debug) std::cout << labels[sourceStop] << std::endl;
    }

    inline RAPTOR::Profile getProfile(const size_t maxTrips = Label::EntrySize) const noexcept {
        RAPTOR::Profile fullProfile;
        ArrivalEntry arrival = ArrivalEntry::Max();
        const Label& label = labels[sourceStop];
        for (size_t i = 1; i < label.size(); i++) {
            int time = never;
            for (size_t j = 0; j < arrival.Size && j < maxTrips; j++) {
                if (label.arrivalTime(i)(j) < time && label.arrivalTime(i)(j) < arrival(j)) {
                    time = label.arrivalTime(i)(j);
                    fullProfile.emplace_back(label.departureTime(i), time, j + 1);
                }
            }
            arrival.minimize(label.arrivalTime(i));
        }
        std::stable_sort(fullProfile.begin(), fullProfile.end());
        RAPTOR::Profile result;
        std::vector<int> arrivalTimes(maxTrips + 1, std::numeric_limits<int>::max());
        const int walkingTravelTime = labels[sourceStop].getWalkingTime();
        for (const RAPTOR::ProfileEntry& entry : fullProfile) {
            if (entry.departureTime < minDepartureTime) continue;
            if (entry.travelTime() >= walkingTravelTime) continue;
            if (entry.arrivalTime - maxDepartureTime >= walkingTravelTime) continue;
            if (entry.departureTime < maxDepartureTime) {
                result.emplace_back(entry);
            } else if (entry.arrivalTime < arrivalTimes[entry.numberOfTrips]) {
                result.emplace_back(entry);
                arrivalTimes[entry.numberOfTrips] = entry.arrivalTime;
                for (size_t i = entry.numberOfTrips + 1; i < arrivalTimes.size(); i++) {
                    arrivalTimes[i] = std::min(arrivalTimes[i], arrivalTimes[i - 1]);
                }
            }
        }
        return result;
    }

    inline RAPTOR::ProfileHandle getProfileHandle(const int maxTrips = Label::EntrySize) const noexcept {
        return RAPTOR::ProfileHandle(getProfile(maxTrips), minDepartureTime, maxDepartureTime, labels[sourceStop].getWalkingTime());
    }

    inline std::vector<Geometry::Point> getArrivalTimeFunction() const noexcept {
        return getProfileHandle().getArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getTravelTimeFunction() const noexcept {
        return getProfileHandle().getTravelTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction() const noexcept {
        return getProfileHandle().getWalkingArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction() const noexcept {
        return getProfileHandle().getWalkingTravelTimeFunction();
    }

    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions() const noexcept {
        return getProfileHandle().getArrivalTimeFunctions();
    }

    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions() const noexcept {
        return getProfileHandle().getTravelTimeFunctions();
    }

    inline void reset() noexcept {
        warning("Resetting...");
        queue.clear();
        for (Label& profile : profileWithoutInitialTransfers) {
            profile.reset_xxxxx();
        }
        for (Label& label : labels) {
            label.reset_();
        }
    }

private:
    inline void clear() noexcept {
        queue.clear();
        for (Label& profile : profileWithoutInitialTransfers) {
            profile.clear();
        }
        for (Label& label : labels) {
            label.clear();
        }
        for (ArrivalEntry& arrival : arrivalByTrip) {
            arrival.maximize();
        }
    }

    inline void initialize() noexcept {
        dijkstra.run(targetStop, noVertex, [&](const Vertex u){
            labels[u].setWalkingTime(dijkstra.getDistance(u));
        });
    }

    inline void runDijkstra(const int keyTime) noexcept {
        while ((!queue.empty()) && (queue.min().key() >= keyTime)) {
            Label* const uLabel = queue.extractFront();
            settle(Vertex(uLabel - &(labels[0])), uLabel->key(), uLabel->lastUnsettledEntry());
            uLabel->settle();
            if (uLabel->hasUnsettledEntries()) {
                queue.update(uLabel);
            }
        }
    }

    inline void settle(const Vertex u, const int departureTime, const ArrivalEntry arrival) noexcept {
        for (const Edge edge : reverseGraph.edgesFrom(u)) {
            Label& vLabel = labels[reverseGraph.get(ToVertex, edge)];
            const int newDepartureTime = departureTime - reverseGraph.get(TravelTime, edge);
            if (TargetPruning) {
                if (labels[sourceStop].dominates(newDepartureTime, arrival)) continue;
            }
            if (arrival.hasParent(reverseGraph.get(ToVertex, edge))) {
                ArrivalEntry partialArrival = arrival;
                partialArrival.removeParent(reverseGraph.get(ToVertex, edge));
                if (vLabel.insertIfImprovement(newDepartureTime, partialArrival)) {
                    if (AllowMultiHopTransfers) queue.update(&vLabel);
                }
            } else {
                if (vLabel.insertIfImprovement(newDepartureTime, arrival)) {
                    if (AllowMultiHopTransfers) queue.update(&vLabel);
                }
            }
        }
        if (!AllowMultiHopTransfers) return;
        for (const Edge edge : minChangeTimeGraph.edgesFrom(u)) {
            Label& vLabel = labels[minChangeTimeGraph.get(ToVertex, edge)];
            const int newDepartureTime = departureTime - minChangeTimeGraph.get(TravelTime, edge);
            if (TargetPruning) {
                if (labels[sourceStop].dominates(newDepartureTime, arrival)) continue;
            }
            if (arrival.hasParent(minChangeTimeGraph.get(ToVertex, edge))) {
                ArrivalEntry partialArrival = arrival;
                partialArrival.removeParent(minChangeTimeGraph.get(ToVertex, edge));
                if (vLabel.insertIfImprovement(newDepartureTime, partialArrival)) {
                    if (AllowMultiHopTransfers) queue.update(&vLabel);
                }
            } else {
                if (vLabel.insertIfImprovement(newDepartureTime, arrival)) {
                    if (AllowMultiHopTransfers) queue.update(&vLabel);
                }
            }
        }
    }

    inline void mergeProfiles() noexcept {
        for (size_t i = 0; i < data.numberOfStops(); i++) {
            labels[i] = Label::Min(labels[i], profileWithoutInitialTransfers[i]);
        }
    }

    inline void buildMinChangeTimeGraph() noexcept {
        Intermediate::TransferGraph graph;
        graph.addVertices(reverseGraph.numVertices());
        Dijkstra<TransferGraph, false> dijkstra(reverseGraph, reverseGraph[TravelTime]);
        for (const StopId stop : data.stops()) {
            const int range = data.stopData[stop].minTransferTime / 2;
            dijkstra.run(stop, noVertex, [&](const Vertex u) {
                if (u == stop || dijkstra.getParent(u) == stop) return;
                const int travelTime = dijkstra.getDistance(u);
                if (travelTime > range || data.isStop(u)) {
                    graph.addEdge(stop, u).set(TravelTime, travelTime);
                    graph.addEdge(u, stop).set(TravelTime, travelTime);
                }
            }, [&]() {
                return false;
            }, [&](const Vertex from, const Edge) {
                return dijkstra.getDistance(from) > range;
            });
        }
        graph.reduceMultiEdgesBy(TravelTime);
        Graph::move(std::move(graph), minChangeTimeGraph);
    }

private:
    const Data& data;
    const TransferGraph& reverseGraph;
    TransferGraph minChangeTimeGraph;

    StopId sourceStop;
    StopId targetStop;
    int minDepartureTime;
    int maxDepartureTime;

    std::vector<ArrivalEntry> arrivalByTrip;
    std::vector<Label> profileWithoutInitialTransfers;
    std::vector<Label> labels;
    ExternalKHeap<2, Label> queue;

    Dijkstra<TransferGraph, false> dijkstra;

};

}
