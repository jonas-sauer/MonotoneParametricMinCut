#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../../DataStructures/Assignment/ConnectionLoadData.h"
#include "../../../DataStructures/Assignment/Profile.h"
#include "../../../DataStructures/Assignment/Settings.h"
#include "../../../DataStructures/Assignment/StopLabel.h"
#include "../../../DataStructures/CSA/Data.h"

#include "../../../Helpers/Vector/Vector.h"

#include "../Profiler.h"

namespace Assignment::Capacities {

template<typename PROFILER = NoPATProfiler, bool USE_TRANSFER_BUFFER_TIMES = false>
class ComputeReassignmentPATs {

public:
    using Profiler = PROFILER;
    static constexpr bool UseTransferBufferTimes = USE_TRANSFER_BUFFER_TIMES;
    using Type = ComputeReassignmentPATs<Profiler, UseTransferBufferTimes>;

    struct PATProfileContainer {
        PATProfileContainer() :
            firstEntry({0}) {
        }

        inline size_t begin(const Vertex stop) const noexcept {
            return firstEntry[stop];
        }

        inline size_t end(const Vertex stop) const noexcept {
            return firstEntry[stop + 1];
        }

        inline size_t size(const Vertex stop) const noexcept {
            return end(stop) - begin(stop);
        }

        inline void resetScanIndex(const Vertex stop) noexcept {
            scanIndex[stop] = end(stop) - 1;
        }

        inline void resetScanIndices() noexcept {
            for (const size_t i : indices(scanIndex)) {
                resetScanIndex(Vertex(i));
            }
        }

        inline void addProfile(const Profile& profile) noexcept {
            entries += profile;
            scanIndex.push_back(entries.size() - 1);
            firstEntry.push_back(entries.size());
        }

        inline const ProfileEntry& findEntry(const Vertex stop, const int time) noexcept {
            while (scanIndex[stop] + 1 < end(stop) && entries[scanIndex[stop] + 1].departureTime >= time) scanIndex[stop]++;
            AssertMsg((scanIndex[stop] + 1 >= end(stop) || entries[scanIndex[stop] + 1].departureTime < time), "Profile is not scanned monotonously (current time: " << time << " previous time: " << entries[scanIndex[stop] + 1] << ")!");
            while (entries[scanIndex[stop]].departureTime < time) {
                scanIndex[stop]--;
                AssertMsg(scanIndex[stop] >= begin(stop), "There seems to be no profile entry for time = " << time << "!");
            }
            return entries[scanIndex[stop]];
        }

        inline void clear() noexcept {
            std::vector<size_t>{0}.swap(firstEntry);
            entries.clear();
            scanIndex.clear();
        }

        std::vector<size_t> firstEntry;
        Profile entries;
        std::vector<size_t> scanIndex;
    };

    struct ConnectionLabel {
        PerceivedTime tripPAT{Unreachable};
        PerceivedTime transferPAT{Unreachable};
        PerceivedTime hopOnPAT{Unreachable};
        PerceivedTime skipPAT{Unreachable};

        inline PerceivedTime bestPAT() const noexcept {
            return std::min(hopOnPAT, skipPAT);
        }
    };

    struct PATData {
        PATData(const size_t numberOfStops, const size_t numberOfConnections) :
            transferDistanceToTarget(numberOfStops, INFTY),
            connectionLabels(numberOfConnections) {
        }

        PATProfileContainer profiles;
        std::vector<int> transferDistanceToTarget;
        std::vector<ConnectionLabel> connectionLabels;

        inline PerceivedTime targetPAT(const CSA::Connection& connection) const noexcept {
            return (transferDistanceToTarget[connection.arrivalStopId] < INFTY) ? (connection.arrivalTime + transferDistanceToTarget[connection.arrivalStopId]) : Unreachable;
        }
    };

public:
    ComputeReassignmentPATs(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const Settings& settings, const std::vector<ConnectionLoadData>& loadData, std::vector<PATData>& patData, const Profiler& profiler = Profiler()) :
        data(data),
        reverseGraph(reverseGraph),
        settings(settings),
        loadData(loadData),
        patData(patData),
        patDataIndex(-1),
        tripPAT(data.numberOfTrips(), Unreachable),
        stopLabels(data.numberOfStops()),
        targetVertex(noVertex),
        profiler(profiler) {
    }

    inline void run(const Vertex target, const size_t index) noexcept {
        profiler.startInitialization();
        patDataIndex = index;
        clear();
        initialize(target);
        profiler.doneInitialization();
        for (ConnectionId i = ConnectionId(data.numberOfConnections() - 1); i < data.numberOfConnections(); i--) {
            profiler.scanConnection(i);
            const CSA::Connection& connection = data.connections[i];
            StopLabel& departureStop = stopLabels[connection.departureStopId];
            const ProfileEntry& skipEntry = departureStop.getSkipEntry();
            ConnectionLabel& label = patData[patDataIndex].connectionLabels[i];

            AssertMsg(skipEntry.departureTime >= connection.departureTime, "Connections are scanned out of order (" << skipEntry.departureTime << " before " << connection.departureTime << ", index: " << i << ")!");
            label.tripPAT = tripPAT[connection.tripId];
            label.transferPAT = stopLabels[connection.arrivalStopId].evaluateWithDelay(connection.arrivalTime, settings.maxDelay, settings.waitingCosts) + settings.transferCosts;
            label.skipPAT = skipEntry.evaluate(connection.departureTime, settings.waitingCosts);
            profiler.evaluateProfile();

            const PerceivedTime travelPAT = label.tripPAT;
            const PerceivedTime walkingPAT = targetPAT(connection);
            const PerceivedTime transferPAT = label.transferPAT;
            //TODO: Allow kick-out?
            const PerceivedTime pat = loadData[i].isFull() ? Unreachable : std::min(std::min(travelPAT, walkingPAT), transferPAT);
            tripPAT[connection.tripId] = pat;
            label.hopOnPAT = pat;
            if (pat >= Unreachable) continue;

            AssertMsg(pat < Unreachable, "Adding infinity PAT = " << pat << "!");
            if (pat >= label.skipPAT) continue;
            stopLabels[connection.departureStopId].addWaitingEntry(ProfileEntry(connection.departureTime, i, pat, settings.waitingCosts));
            profiler.addToProfile();
            const int bufferTime = data.minTransferTime(connection.departureStopId);
            stopLabels[connection.departureStopId].addTransferEntry(ProfileEntry(connection.departureTime, i, pat, 0, bufferTime, settings.walkingCosts, settings.waitingCosts), profiler);
            for (const Edge edge : reverseGraph.edgesFrom(connection.departureStopId)) {
                const Vertex from = reverseGraph.get(ToVertex, edge);
                if (!data.isStop(from)) continue;
                stopLabels[from].addTransferEntry(ProfileEntry(connection.departureTime, i, pat, reverseGraph.get(TravelTime, edge), UseTransferBufferTimes ? bufferTime : 0, settings.walkingCosts, settings.waitingCosts), profiler);
                profiler.relaxEdge(edge);
            }
        }
        patData[patDataIndex].profiles.clear();
        size_t numEntries = 0;
        for (const StopId stop : data.stops()) {
            numEntries += stopLabels[stop].getWaitingProfile().size();
        }
        patData[patDataIndex].profiles.entries.reserve(numEntries);
        for (const StopId stop : data.stops()) {
            patData[patDataIndex].profiles.addProfile(stopLabels[stop].getWaitingProfile());
        }
    }

    inline Profiler& getProfiler() noexcept {
        return profiler;
    }

private:
    inline void clear() noexcept {
        Vector::fill(tripPAT, Unreachable);
        Vector::fill(stopLabels);
        if (reverseGraph.isVertex(targetVertex)) cleanUp();
    }

    inline void initialize(const Vertex target) noexcept {
        targetVertex = target;
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            patData[patDataIndex].transferDistanceToTarget[stop] = (settings.walkingCosts + 1) * reverseGraph.get(TravelTime, edge);
        }
        if (data.isStop(targetVertex)) patData[patDataIndex].transferDistanceToTarget[targetVertex] = 0;
    }

    inline void cleanUp() noexcept {
        for (const Edge edge : reverseGraph.edgesFrom(targetVertex)) {
            profiler.relaxEdge(edge);
            const Vertex stop = reverseGraph.get(ToVertex, edge);
            if (!data.isStop(stop)) continue;
            patData[patDataIndex].transferDistanceToTarget[stop] = INFTY;
        }
        if (data.isStop(targetVertex)) patData[patDataIndex].transferDistanceToTarget[targetVertex] = INFTY;
    }

    inline PerceivedTime targetPAT(const CSA::Connection& connection) const noexcept {
        const int distance = patData[patDataIndex].transferDistanceToTarget[connection.arrivalStopId];
        return (distance < INFTY) ? (connection.arrivalTime + distance) : Unreachable;
    }

private:
    const CSA::Data& data;
    const CSA::TransferGraph& reverseGraph;
    const Settings& settings;
    const std::vector<ConnectionLoadData>& loadData;
    std::vector<PATData>& patData;
    size_t patDataIndex;

    std::vector<PerceivedTime> tripPAT;
    std::vector<StopLabel> stopLabels;
    Vertex targetVertex;

    Profiler profiler;

};

}
