#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../InitialTransfers.h"
#include "../../Dijkstra/Dijkstra.h"

#include "../../../Helpers/Meta.h"
#include "../../../Helpers/Helpers.h"
#include "../../../Helpers/Vector/Vector.h"
#include "../../../Helpers/Vector/Permutation.h"

#include "../../../DataStructures/Container/Heap.h"
#include "../../../DataStructures/Container/Set.h"
#include "../../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../../DataStructures/RAPTOR/Entities/Profile.h"
#include "../../../DataStructures/Geometry/Rectangle.h"
#include "../Profiler.h"
#include "UPRAPTORModule.h"

namespace RAPTOR::RangeRAPTOR {

template<bool USE_TARGET_BUCKETS, bool USE_DFS_ORDER, bool DEBUG = false, bool TRIP_PRUNING = true>
class UPRangeRAPTOR {

public:
    static constexpr bool UseTargetBuckets = USE_TARGET_BUCKETS;
    static constexpr bool UseDFSOrder = USE_DFS_ORDER;
    static constexpr bool Debug = DEBUG;
    using Profiler = Meta::IF<Debug, SimpleProfiler, NoProfiler>;
    static constexpr bool TripPruning = TRIP_PRUNING;
    using RaptorModule = UPRAPTORModule<UseTargetBuckets, UseDFSOrder, Profiler, TripPruning>;
    using InitialAndFinalTransfers = typename RaptorModule::InitialAndFinalTransfers;
    using Type = UPRangeRAPTOR<UseTargetBuckets, UseDFSOrder, Debug, TripPruning>;

private:
    struct DepartureLabel {
        DepartureLabel(const int departureTime, const StopId departureStop) :
            departureTime(departureTime),
            departureStop(departureStop) {
        }
        inline bool operator<(const DepartureLabel& other) const noexcept {
            return departureTime > other.departureTime;
        }
        int departureTime;
        StopId departureStop;
    };

    struct FinalTransferLabel {
        FinalTransferLabel(const int departureTime, const int arrivalTime, const StopId stop) :
            departureTime(departureTime),
            arrivalTime(arrivalTime),
            stop(stop) {
        }

        int departureTime;
        int arrivalTime;
        StopId stop;
    };

    struct SegmentedProfileEntry {
        SegmentedProfileEntry(const int departureTime = never, const int arrivalTime = never) :
            departureTime(departureTime),
            arrivalTime(arrivalTime) {
        }

        inline int travelTime() const noexcept {
            return arrivalTime - departureTime;
        }

        inline bool operator<(const SegmentedProfileEntry& other) const noexcept {
            return (departureTime < other.departureTime) || ((departureTime == other.departureTime) && (arrivalTime < other.arrivalTime));
        }

        inline bool operator==(const SegmentedProfileEntry& other) const noexcept {
            return (departureTime == other.departureTime) && (arrivalTime == other.arrivalTime);
        }

        inline bool dominates(const SegmentedProfileEntry& other) const noexcept {
            return (departureTime >= other.departureTime) && (arrivalTime <= other.arrivalTime);
        }

        inline friend std::ostream& operator<<(std::ostream& out, const SegmentedProfileEntry& entry) noexcept {
            return out << "Departure: " << String::secToTime(entry.departureTime) << ", Arrival: " << String::secToTime(entry.arrivalTime);
        }

        int departureTime;
        int arrivalTime;
    };

    class ProfileSegment : public std::vector<SegmentedProfileEntry> {
    public:
        inline bool add(const int departureTime, const int arrivalTime) noexcept {
            if (!this->empty() && this->back().arrivalTime <= arrivalTime) return false;
            this->emplace_back(departureTime, arrivalTime);
            return true;
        }

        inline void reverse() noexcept {
            Vector::reverse(*this);
        }

        inline void sort() noexcept {
            std::sort(this->begin(), this->end(), [&](const SegmentedProfileEntry& a, const SegmentedProfileEntry& b) {
               return b < a;
            });
        }
    };

    class SegmentedProfile : public std::vector<ProfileSegment> {
    public:
        inline bool add(const ArrivalLabel arrival, const int departureTime) noexcept {
            Assert(arrival.numberOfTrips != size_t(-1), "Arrival label has invalid number of trips!");
            if (arrival.numberOfTrips >= this->size()) this->resize(arrival.numberOfTrips + 1);
            return (*this)[arrival.numberOfTrips].add(departureTime, arrival.arrivalTime);
        }

        inline void reverse() noexcept {
            for (ProfileSegment& segment : *this) {
                segment.reverse();
            }
        }

        inline void sort() noexcept {
            for (ProfileSegment& segment : *this) {
                segment.sort();
            }
        }

        inline void prune(const int maxDepartureTime) noexcept {
            int bestArrivalTime = never;
            for (ProfileSegment& segment : *this) {
                size_t size = segment.size();
                while (size > 0 && segment[size - 1].departureTime >= maxDepartureTime) {
                    size--;
                }
                if (size == segment.size()) continue;
                if (segment[size].arrivalTime < bestArrivalTime) {
                    bestArrivalTime = segment[size].arrivalTime;
                    size++;
                }
                segment.resize(size);
            }
        }

        inline SegmentedProfile getCutoffProfile(const size_t maxTrips) const noexcept {
            SegmentedProfile result = *this;
            result.resize(maxTrips);
            return result;
        }

        inline Profile getDesegmentedProfile() const noexcept {
            Profile result;
            std::vector<size_t> nextEntry(this->size(), 0);
            while (true) {
                const Profile::iterator resultEnd = result.end();
                for (size_t i = 0; i < nextEntry.size(); i++) {
                    if (nextEntry[i] >= (*this)[i].size()) continue;
                    const SegmentedProfileEntry& entry = (*this)[i][nextEntry[i]];
                    result.push_back(ProfileEntry(entry.departureTime, entry.arrivalTime, i));
                    nextEntry[i]++;
                }
                if (result.end() == resultEnd) break;
                std::sort(resultEnd, result.end(), [&](const ProfileEntry& a, const ProfileEntry& b) {
                   return b < a;
                });
            }
            return result;
        }

        inline Profile getDesegmentedProfile(const size_t maxTrips) const noexcept {
            return getCutoffProfile(maxTrips).getDesegmentedProfile();
        }

        inline size_t numEntries() const noexcept {
            size_t result = 0;
            for (const ProfileSegment& segment : (*this)) {
                result += segment.size();
            }
            return result;
        }
    };

    inline static RaptorModule createRaptorModule(Data& data, CH::CH& chData, IndexedSet<false, Vertex>& targets, const bool reorder) noexcept {
        if (reorder) {
            return RaptorModule::Reordered(data, chData, targets);
        } else {
            return RaptorModule::NotReordered(data, chData, targets);
        }
    }

public:
    UPRangeRAPTOR(Data& data, CH::CH& chData, IndexedSet<false, Vertex>& targets, const bool reorder) :
        raptorModule(createRaptorModule(data, chData, targets, reorder)),
        initialAndFinalTransfers(raptorModule.getInitialAndFinalTransfers()),
        data(data),
        sourceVertex(noVertex),
        targets(targets),
        minDepartureTime(never),
        maxDepartureTime(never),
        profileByVertex(chData.numVertices()) {
    }

    inline void run(const Vertex source, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        clear();
        sourceVertex = source;
        targetVertices = targets;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;

        raptorModule.template runInitialize<true>(source, maxDepartureTime);
        raptorModule.template runInitialTransfers();

        collectDepartures();
        for (size_t i = 0; i < departures.size(); i++) {
            raptorModule.template runInitialize<false>(source, departures[i].departureTime);
            raptorModule.runAddSource(departures[i].departureStop, departures[i].departureTime + raptorModule.getWalkingTravelTime(departures[i].departureStop));
            while (i + 1 < departures.size() && departures[i].departureTime == departures[i + 1].departureTime) {
                i++;
                raptorModule.runAddSource(departures[i].departureStop, departures[i].departureTime + raptorModule.getWalkingTravelTime(departures[i].departureStop));
            }
            if constexpr (Debug) std::cout << "Departure Time: " << departures[i].departureTime << std::endl;
            raptorModule.template runRounds(maxRounds);
            collectArrivals(departures[i].departureTime, raptorModule.getReachedStops());
        }

        relaxFinalTransfers();

        for (const Vertex target : targets) {
            profileByVertex[target].reverse();
            profileByVertex[target].prune(maxDepartureTime);
        }
    }

    inline Profile getProfile(const Vertex vertex) const noexcept {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        return profileByVertex[vertex].getDesegmentedProfile();
    }

    inline Profile getProfile(const Vertex vertex, const size_t maxTrips) const noexcept {
        return profileByVertex[vertex].getDesegmentedProfile(maxTrips);
    }

    inline ProfileHandle getProfileHandle(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return ProfileHandle(getProfile(vertex, maxTrips), minDepartureTime, maxDepartureTime, raptorModule.getWalkingTravelTime(vertex));
    }

    inline std::vector<Geometry::Point> getArrivalTimeFunction(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getTravelTimeFunction(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getTravelTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction(const Vertex vertex) noexcept {
        return getProfileHandle(vertex).getWalkingArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction(const Vertex vertex) noexcept {
        return getProfileHandle(vertex).getWalkingTravelTimeFunction();
    }

    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getArrivalTimeFunctions();
    }

    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getTravelTimeFunctions();
    }

    inline void reset() noexcept {
        raptorModule.reset();
        std::vector<SegmentedProfile>(profileByVertex.size()).swap(profileByVertex);
    }

private:
    inline void clear() noexcept {
        for (SegmentedProfile& profile : profileByVertex) {
            profile.clear();
        }
        finalTransferLabels.clear();
    }

    inline void collectDepartures() noexcept {
        departures.clear();
        for (const RouteId route : data.routes()) {
            const size_t numberOfStops = data.numberOfStopsInRoute(route);
            const StopId* stops = data.stopArrayOfRoute(route);
            const StopEvent* stopEvents = data.firstTripOfRoute(route);
            for (uint32_t stopEventIndex = 0; stopEventIndex < data.numberOfStopEventsInRoute(route); stopEventIndex++) {
                if ((stopEventIndex + 1) % numberOfStops == 0) continue;
                const StopId stop = stops[stopEventIndex % numberOfStops];
                const int walkingTime = raptorModule.getWalkingTravelTime(stop);
                if (walkingTime == intMax) continue;
                const int departureTime = stopEvents[stopEventIndex].departureTime - walkingTime;
                if (departureTime < minDepartureTime) continue;
                departures.emplace_back(departureTime, stop);
            }
        }
        sort(departures);
        if constexpr (Debug) std::cout << "Number of departures: " << departures.size() << std::endl;
    }

    inline void collectArrivals(const int departureTime, const IndexedSet<false, StopId>& reachedStops) noexcept {
        for (const StopId stop : reachedStops) {
            for (const ArrivalLabel arrival : raptorModule.getArrivals(stop)) {
                if (arrival.numberOfTrips == 0) continue;
                if (profileByVertex[stop].add(arrival, departureTime)) {
                    if (finalTransferLabels.size() < arrival.numberOfTrips) finalTransferLabels.resize(arrival.numberOfTrips);
                    finalTransferLabels[arrival.numberOfTrips - 1].emplace_back(departureTime, arrival.arrivalTime, stop);
                }
            }
        }
    }

    //TODO: How can later rounds be pruned with earlier rounds?
    inline void relaxFinalTransfers() noexcept {
        for (size_t numTrips = 1; numTrips - 1 < finalTransferLabels.size(); numTrips++) {
            const std::vector<FinalTransferLabel>& labels = finalTransferLabels[numTrips - 1];
            size_t i = 0;
            std::vector<int> departureTimes(RaptorModule::InitialAndFinalTransfers::MaxSources, never);
            while (i < labels.size()) {
                initialAndFinalTransfers.clearDepartureTimes();
                for (size_t j = 0; j < RaptorModule::InitialAndFinalTransfers::MaxSources; j++) {
                    if (i >= labels.size()) break;
                    departureTimes[j] = labels[i].departureTime;
                    initialAndFinalTransfers.addDepartureTime(departureTimes[j]);
                    while (labels[i].departureTime == departureTimes[j]) {
                        initialAndFinalTransfers.addSource(labels[i].stop, labels[i].arrivalTime, labels[i].stop);
                        i++;
                    }
                }
                initialAndFinalTransfers.runFinalTransfers();
                //TODO: Only retrieve entries that were actually updated
                for (const Vertex target : targets) {
                    for (size_t j = 0; j < initialAndFinalTransfers.numDepartures(); j++) {
                        const ArrivalLabel label(initialAndFinalTransfers.getArrivalTime(j, target), numTrips);
                        profileByVertex[target].add(label, departureTimes[j]);
                    }
                }
            }
        }
    }

private:
    RaptorModule raptorModule;
    InitialAndFinalTransfers& initialAndFinalTransfers;
    const Data& data;

    Vertex sourceVertex;
    const IndexedSet<false, Vertex> targets;
    int minDepartureTime;
    int maxDepartureTime;
    std::vector<DepartureLabel> departures;

    std::vector<SegmentedProfile> profileByVertex;
    IndexedSet<false, Vertex> targetVertices;

    std::vector<std::vector<FinalTransferLabel>> finalTransferLabels;
};

}
