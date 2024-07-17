#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <type_traits>

#include "../../Helpers/Assert.h"
#include "../../Helpers/Timer.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Vector.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/CSA/Entities/Journey.h"
#include "../../DataStructures/CSA/Profile/Profile.h"
#include "../../DataStructures/CSA/Profile/ContiguousProfileVector.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"
#include "../RAPTOR/InitialTransfers.h"
#include "Profiler.h"

namespace CSA {

template<bool USE_TARGET_BUCKETS, bool USE_DFS_ORDER, typename PROFILER = NoProfiler, bool CONTIGUOUS_PROFILES = true>
class UPProfileCSA {

public:
    static constexpr bool UseTargetBuckets = USE_TARGET_BUCKETS;
    static constexpr bool UseDFSOrder = USE_DFS_ORDER;
    using Profiler = PROFILER;
    constexpr static bool ContiguousProfiles = CONTIGUOUS_PROFILES;
    using ProfileVectorType = std::conditional_t<ContiguousProfiles, ContiguousProfileVector, ProfileVectorWrapper>;
    using InitialAndFinalTransfers = RAPTOR::ProfileCSAInitialAndFinalTransfers<false, UseTargetBuckets>;
    using Type = UPProfileCSA<UseTargetBuckets, UseDFSOrder, Profiler, ContiguousProfiles>;

private:
    struct QueryData {
        QueryData(const Data& oldData, const TransferGraph& oldReverseGraph, const CH::CH& oldCHData, const IndexedSet<false, Vertex>& oldSources, const bool reorder) :
            data(oldData),
            reverseGraph(oldReverseGraph),
            chData(oldCHData),
            sources(oldSources) {
            const Order chOrder = vertexOrder(chData);
            if (reorder) {
                const Order fullOrder = chOrder.splitAt(data.numberOfStops());
                Order stopOrder = fullOrder;
                stopOrder.resize(data.numberOfStops());

                data.applyStopOrder(stopOrder);
                reverseGraph.applyVertexOrder(stopOrder);
                chData.applyVertexOrder(fullOrder);
                sources.applyPermutation(Permutation(Construct::Invert, fullOrder));

                size_t numStops = 0;
                size_t numVertices = data.numberOfStops();
                for (const size_t i : chOrder) {
                    if (i < data.numberOfStops()) {
                        phastOrder.emplace_back(numStops++);
                    } else {
                        phastOrder.emplace_back(numVertices++);
                    }
                }
            } else {
                phastOrder = chOrder;
            }
        }

        Data data;
        TransferGraph reverseGraph;
        CH::CH chData;
        Order phastOrder;
        IndexedSet<false, Vertex> sources;
    };

    struct FinalTransferLabel {
        FinalTransferLabel(const int departureTime = never, const StopId sourceStop = noStop) :
            departureTime(departureTime),
            sourceStop(sourceStop) {
        }

        int departureTime;
        StopId sourceStop;
    };

public:
    UPProfileCSA(const Data& oldData, const TransferGraph& oldReverseGraph, const CH::CH& oldCHData, const IndexedSet<false, Vertex>& oldSources, const bool reorder, const Profiler& profilerTemplate = Profiler()) :
        queryData(oldData, oldReverseGraph, oldCHData, oldSources, reorder),
        data(queryData.data),
        reverseGraph(queryData.reverseGraph),
        initialAndFinalTransfers(queryData.chData, std::move(queryData.phastOrder), data.numberOfStops(), queryData.sources, BACKWARD),
        sourceVertices(queryData.sources),
        sourceStopId(queryData.chData.numVertices(), noStop),
        targetVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never),
        arrivalTimeOffset(computeArrivalTimeOffset()),
        arrivalTimeMap(48 * 60 * 60 + arrivalTimeOffset),
        tripArrivalTime(data.numberOfTrips(), never),
        profiles(computeSourceMapping()),
        profiler(profilerTemplate) {
        Assert(Vector::isSorted(data.connections), "Connections must be sorted in ascending order!");
        profiler.registerPhases({PHASE_CLEAR, PHASE_INITIALIZATION, PHASE_CONNECTION_SCAN, PHASE_FINAL_TRANSFERS});
        profiler.registerMetrics({METRIC_CONNECTIONS, METRIC_EDGES, METRIC_STOPS_BY_TRIP, METRIC_STOPS_BY_TRANSFER});
        profiler.initialize();
    }

    inline static Order vertexOrder(const CH::CH& chData) noexcept {
        if constexpr (UseDFSOrder) {
            return Order(Vector::reverse(CH::getOrder(chData)));
        } else {
            return Order(CH::getLevelOrderTopDown(chData));
        }
    }

    inline void run(const Vertex target, const int minTime = 0, const int maxTime = 24 * 60 * 60) noexcept {
        profiler.start();

        profiler.startPhase();
        clear();
        profiler.donePhase(PHASE_CLEAR);

        profiler.startPhase();
        initialAndFinalTransfers.initialize();
        targetVertex = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        initialAndFinalTransfers.runInitialTransfers(targetVertex);
        profiler.donePhase(PHASE_INITIALIZATION);

        profiler.startPhase();
        scanConnections();
        profiler.donePhase(PHASE_CONNECTION_SCAN);

        profiler.startPhase();
        relaxFinalTransfers();
        profiler.donePhase(PHASE_FINAL_TRANSFERS);

        profiler.done();
    }

    inline RAPTOR::Profile getProfile(const Vertex vertex) noexcept {
        const StopId stop = sourceStopId[vertex];
        RAPTOR::Profile result;
        const int directTransferTime = getDirectWalkingTime(vertex);
        for (size_t i = 1; i < profiles.size(stop); i++) {
            const ProfileEntry& entry = profiles.entry(stop, i);
            if (entry.departureTime > maxDepartureTime) {
                if (entry.arrivalTime - maxDepartureTime >= directTransferTime) continue;
                if (i+1 < profiles.size(stop) && profiles.entry(stop, i+1).departureTime >= maxDepartureTime) continue;
            }
            result.emplace_back(entry.departureTime, entry.arrivalTime);
        }
        std::stable_sort(result.begin(), result.end());
        return result;
    }

    inline RAPTOR::ProfileHandle getProfileHandle(const Vertex vertex) noexcept {
        return RAPTOR::ProfileHandle(getProfile(vertex), minDepartureTime, maxDepartureTime, getDirectWalkingTime(vertex));
    }

    inline int getDirectWalkingTime(const Vertex vertex) noexcept {
        return initialAndFinalTransfers.getInitialDistance(vertex);
    }

    inline const Profiler& getProfiler() const noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        targetVertex = noVertex;
        minDepartureTime = never;
        maxDepartureTime = never;
        arrivalTimeMap.clear();
        Vector::fill(tripArrivalTime, never);
        profiles.clear();
    }

    inline int getArrivalTimeViaWalking(const StopId stop, const int time) noexcept {
        const int distance = initialAndFinalTransfers.getInitialDistance(stop);
        if (distance == never) return distance;
        return time + distance;
    }

    inline void updateProfiles(const Connection& connection, const int arrivalTime) noexcept {
        const int departureTime = connection.departureTime - data.minTransferTime(connection.departureStopId);
        if (departureTime < minDepartureTime) return;
        if (arrivalTime - departureTime >= initialAndFinalTransfers.getInitialDistance(connection.departureStopId)) return;
        profiler.countMetric(METRIC_STOPS_BY_TRIP);
        if (!profiles.addEntry(connection.departureStopId, departureTime, arrivalTime)) return;
        for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            profiler.countMetric(METRIC_EDGES);
            const int stopDepartureTime = departureTime - reverseGraph.get(TravelTime, edge);
            if (stopDepartureTime < minDepartureTime) continue;
            if (arrivalTime - stopDepartureTime >= initialAndFinalTransfers.getInitialDistance(stop)) continue;
            profiler.countMetric(METRIC_STOPS_BY_TRANSFER);
            profiles.addEntry(stop, stopDepartureTime, arrivalTime);
        }
    }

    inline void scanConnections() noexcept {
        for (size_t i = data.numberOfConnections() - 1; i != size_t(-1); i--) {
            const Connection& connection = data.connections[i];
            if (connection.departureTime < minDepartureTime) break;
            profiler.countMetric(METRIC_CONNECTIONS);

            const int arrivalTimeViaWalking = getArrivalTimeViaWalking(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTimeViaTrip = tripArrivalTime[connection.tripId];
            const int arrivalTimeViaTransfer = profiles.getEarliestArrivalTime(connection.arrivalStopId, connection.arrivalTime);
            const int arrivalTime = std::min(arrivalTimeViaWalking, std::min(arrivalTimeViaTrip, arrivalTimeViaTransfer));

            tripArrivalTime[connection.tripId] = arrivalTime;
            updateProfiles(connection, arrivalTime);
        }
    }

    inline size_t computeSourceMapping() noexcept {
        StopId stopId(data.numberOfStops());
        for (const Vertex source : sourceVertices) {
            sourceStopId[source] = data.isStop(source) ? StopId(source) : stopId;
            stopId++;
        }
        return stopId;
    }

    inline void relaxFinalTransfers() noexcept {
        std::vector<std::vector<FinalTransferLabel>> finalTransferLabels;
        std::vector<int> arrivalTimes;
        for (const StopId stop : data.stops()) {
            for (size_t i = 1; i < profiles.size(stop); i++) {
                const ProfileEntry& entry = profiles.entry(stop, i);
                const int mapIndex = entry.arrivalTime + arrivalTimeOffset;
                if (!arrivalTimeMap.contains(mapIndex)) {
                    arrivalTimeMap.insert(mapIndex, finalTransferLabels.size());
                    arrivalTimes.emplace_back(entry.arrivalTime);
                    finalTransferLabels.emplace_back();
                }
                finalTransferLabels[arrivalTimeMap[mapIndex]].emplace_back(entry.departureTime, stop);
            }
        }
        std::sort(arrivalTimes.begin(), arrivalTimes.end());
        for (const int arrivalTime : arrivalTimes) {
            initialAndFinalTransfers.setArrivalTime(-arrivalTime);
            for (const FinalTransferLabel& label : finalTransferLabels[arrivalTimeMap[arrivalTime + arrivalTimeOffset]]) {
                initialAndFinalTransfers.addSource(label.sourceStop, -label.departureTime);
            }
            initialAndFinalTransfers.runFinalTransfers();
            for (const Vertex sourceVertex : initialAndFinalTransfers.getReachedTargets()) {
                Assert(sourceVertices.contains(sourceVertex), "Reached non-source vertex " << sourceVertex);
                const int departureTime = -initialAndFinalTransfers.getDepartureTime(sourceVertex);
                if (departureTime < minDepartureTime) continue;
                if (arrivalTime - departureTime >= initialAndFinalTransfers.getInitialDistance(sourceVertex)) continue;
                profiles.addEntry(sourceStopId[sourceVertex], departureTime, arrivalTime);
            }
        }
    }

    inline int computeArrivalTimeOffset() const noexcept {
        int result = 0;
        for (const Connection& connection : data.connections) {
            result = std::max(result, -connection.arrivalTime);
        }
        return result;
    }

private:
    const QueryData queryData;
    const Data& data;
    const TransferGraph& reverseGraph;
    InitialAndFinalTransfers initialAndFinalTransfers;

    const IndexedSet<false, Vertex>& sourceVertices;
    std::vector<StopId> sourceStopId;
    Vertex targetVertex;
    int minDepartureTime;
    int maxDepartureTime;

    int arrivalTimeOffset;
    IndexedMap<size_t, true, int> arrivalTimeMap;

    std::vector<int> tripArrivalTime;
    ProfileVectorType profiles;

    Profiler profiler;
};
}
