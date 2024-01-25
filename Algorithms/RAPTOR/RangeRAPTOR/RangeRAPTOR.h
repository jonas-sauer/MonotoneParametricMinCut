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
#include "DijkstraRAPTORModule.h"
#include "TransitiveRAPTORModule.h"
#include "ULTRARAPTORModule.h"

namespace RAPTOR::RangeRAPTOR {

template<typename RAPTOR_MODULE, bool ONE_TO_ONE = false, bool DEBUG = false, bool CORRECT_DEPARTURE_TIMES = true>
class RangeRAPTOR {

public:
    using RaptorModule = RAPTOR_MODULE;
    static constexpr bool OneToOne = ONE_TO_ONE;
    static constexpr bool Debug = DEBUG;
    static constexpr bool CorrectDepartureTimes = CORRECT_DEPARTURE_TIMES;
    using InitialTransferGraph = typename RaptorModule::InitialTransferGraph;
    using Type = RangeRAPTOR<RaptorModule, OneToOne, Debug, CorrectDepartureTimes>;

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

    struct ProfileWrapper {
        inline void add(const ArrivalLabel arrival, const int departureTime) noexcept {
            Assert(arrival.numberOfTrips != size_t(-1), "Arrival label has invalid number of trips!");
            if (arrival.numberOfTrips >= arrivalTimeByNumberOfTrips.size()) {
                arrivalTimeByNumberOfTrips.resize(arrival.numberOfTrips + 1, never);
            }
            if (arrivalTimeByNumberOfTrips[arrival.numberOfTrips] <= arrival.arrivalTime) return;
            arrivalTimeByNumberOfTrips[arrival.numberOfTrips] = arrival.arrivalTime;
            profile.emplace_back(arrival, departureTime);
        }

        inline void sort() noexcept {
            std::sort(profile.begin(), profile.end(), [&](const ProfileEntry& a, const ProfileEntry& b) {
               return b < a;
            });
        }

        inline void prune(const int maxDepartureTime) noexcept {
            size_t removedEntriesCount = 0;
            std::vector<int> arrivalTimes(arrivalTimeByNumberOfTrips.size(), never);
            for (size_t i = 0; i < profile.size(); i++) {
                const ProfileEntry& entry = profile[i];
                if (entry.departureTime < maxDepartureTime) continue;
                if (entry.arrivalTime >= arrivalTimes[entry.numberOfTrips]) {
                    removedEntriesCount++;
                    continue;
                } else {
                    arrivalTimes[entry.numberOfTrips] = entry.arrivalTime;
                    for (size_t i = entry.numberOfTrips + 1; i < arrivalTimes.size(); i++) {
                        arrivalTimes[i] = std::min(arrivalTimes[i], arrivalTimes[i - 1]);
                    }
                }
                if (removedEntriesCount > 0) {
                    profile[i - removedEntriesCount] = entry;
                }
            }
            if (removedEntriesCount > 0) {
                profile.resize(profile.size() - removedEntriesCount);
            }
        }

        inline void clear() noexcept {
            arrivalTimeByNumberOfTrips.clear();
            profile.clear();
        }

        inline const Profile& getProfile() const noexcept {
            return profile;
        }

        inline Profile getProfile(const size_t maxTrips) const noexcept {
            Profile result;
            for (const ProfileEntry& entry : profile) {
                if (entry.numberOfTrips > maxTrips) continue;
                result.emplace_back(entry);
            }
            return result;
        }

        std::vector<int> arrivalTimeByNumberOfTrips;
        Profile profile;
    };

public:
    template<typename ATTRIBUTE>
    RangeRAPTOR(const Data& data, const Data& reverseData, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight) :
        forwardRaptor(data, forwardGraph, backwardGraph, weight),
        backwardRaptor(reverseData, backwardGraph, forwardGraph, weight),
        data(data),
        reverseData(reverseData),
        sourceVertex(noVertex),
        minDepartureTime(never),
        maxDepartureTime(never),
        profileByVertex(OneToOne ? 0 : forwardGraph.numVertices()) {
    }

    template<typename T = CHGraph, typename = std::enable_if_t<Meta::Equals<T, CHGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    RangeRAPTOR(const Data& data, const Data& reverseData, const CH::CH& chData) :
        RangeRAPTOR(data, reverseData, chData.forward, chData.backward, Weight) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    RangeRAPTOR(const Data& data, const Data& reverseData, const TransferGraph& forwardGraph, const TransferGraph& backwardGraph) :
        RangeRAPTOR(data, reverseData, forwardGraph, backwardGraph, TravelTime) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    RangeRAPTOR(const Data& data, const Data& reverseData) :
        RangeRAPTOR(data, reverseData, data.transferGraph, reverseData.transferGraph) {
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline void runOneToAll(const Vertex source, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        run(source, IndexedSet<false, Vertex>(Construct::Complete, data.transferGraph.numVertices()), minTime, maxTime, maxRounds);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline void runOneToAllStops(const Vertex source, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        run(source, IndexedSet<false, Vertex>(Construct::Complete, data.numberOfStops()), minTime, maxTime, maxRounds);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline void run(const Vertex source, const std::vector<Vertex>& targets, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        run(source, IndexedSet<false, Vertex>(data.transferGraph.numVertices(), targets), minTime, maxTime, maxRounds);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline void run(const Vertex source, const IndexedSet<false, Vertex>& targets, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        clear();
        sourceVertex = source;
        targetVertices = targets;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;

        if constexpr (CorrectDepartureTimes) {
            forwardRaptor.template run<true>(source, maxDepartureTime, noVertex, maxRounds);
            collectArrivals<true>(maxDepartureTime, targets, forwardRaptor.getReachedVertices());
        } else {
            forwardRaptor.template runInitialize<true>(source, maxDepartureTime);
            forwardRaptor.template runInitialTransfers();
        }

        collectDepartures();
        for (size_t i = 0; i < departures.size(); i++) {
            forwardRaptor.template runInitialize<false>(source, departures[i].departureTime);
            forwardRaptor.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardRaptor.getWalkingTravelTime(departures[i].departureStop));
            while (i + 1 < departures.size() && departures[i].departureTime == departures[i + 1].departureTime) {
                i++;
                forwardRaptor.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardRaptor.getWalkingTravelTime(departures[i].departureStop));
            }
            if constexpr (Debug) std::cout << "Departure Time: " << departures[i].departureTime << std::endl;
            forwardRaptor.template runRounds(maxRounds);
            collectArrivals(departures[i].departureTime, targets, forwardRaptor.getReachedVertices());
        }

        for (const Vertex target : targets) {
            Vector::reverse(profileByVertex[target].profile);
            if constexpr (!CorrectDepartureTimes) profileByVertex[target].prune(maxDepartureTime);
        }
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline void run(const Vertex source, const Vertex target, const int minTime = 0, const int maxTime = 60 * 60 * 24, const size_t maxRounds = INFTY) noexcept {
        clear();
        sourceVertex = source;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;

        if constexpr (CorrectDepartureTimes) {
            forwardRaptor.template run<true>(source, maxDepartureTime, target, maxRounds);
            collectArrivals<true>(maxDepartureTime, target, targetProfile);
        } else {
            forwardRaptor.template runInitialize<true>(source, 0, target);
            forwardRaptor.template runInitialTransfers();
        }

        collectDepartures();
        for (size_t i = 0; i < departures.size(); i++) {
            forwardRaptor.template runInitialize<false>(source, departures[i].departureTime, target);
            forwardRaptor.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardRaptor.getWalkingTravelTime(departures[i].departureStop));
            while (i + 1 < departures.size() && departures[i].departureTime == departures[i + 1].departureTime) {
                i++;
                forwardRaptor.runAddSource(departures[i].departureStop, departures[i].departureTime + forwardRaptor.getWalkingTravelTime(departures[i].departureStop));
            }
            if constexpr (Debug) std::cout << "Departure Time: " << departures[i].departureTime << std::endl;
            forwardRaptor.template runRounds(maxRounds);
            if (forwardRaptor.targetReached()) {
                collectArrivals(departures[i].departureTime, target, targetProfile);
            }
        }

        Vector::reverse(targetProfile.profile);
        if constexpr (!CorrectDepartureTimes) targetProfile.prune(maxDepartureTime);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline const Profile& getProfile() const noexcept {
        return targetProfile.getProfile();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline const Profile& getProfile(const Vertex vertex) const noexcept {
        Assert(targetVertices.contains(vertex), "Vertex " << vertex << " is not a target!");
        return profileByVertex[vertex].getProfile();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline Profile getProfile(const size_t maxTrips) const noexcept {
        return targetProfile.getProfile(maxTrips);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline Profile getProfile(const Vertex vertex, const size_t maxTrips) const noexcept {
        return profileByVertex[vertex].getProfile(maxTrips);
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline ProfileHandle getProfileHandle(const size_t maxTrips = 100) noexcept {
        return ProfileHandle(getProfile(maxTrips), minDepartureTime, maxDepartureTime, forwardRaptor.getWalkingTravelTime());
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline ProfileHandle getProfileHandle(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return ProfileHandle(getProfile(vertex, maxTrips), minDepartureTime, maxDepartureTime, forwardRaptor.getWalkingTravelTime(vertex));
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<Geometry::Point> getArrivalTimeFunction(const size_t maxTrips = 100) noexcept {
        return getProfileHandle(maxTrips).getArrivalTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<Geometry::Point> getArrivalTimeFunction(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getArrivalTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<Geometry::Point> getTravelTimeFunction(const size_t maxTrips = 100) noexcept {
        return getProfileHandle(maxTrips).getTravelTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<Geometry::Point> getTravelTimeFunction(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getTravelTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction() noexcept {
        return getProfileHandle().getWalkingArrivalTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction(const Vertex vertex) noexcept {
        return getProfileHandle(vertex).getWalkingArrivalTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction() noexcept {
        return getProfileHandle().getWalkingTravelTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction(const Vertex vertex) noexcept {
        return getProfileHandle(vertex).getWalkingTravelTimeFunction();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions(const size_t maxTrips = 100) noexcept {
        return getProfileHandle(maxTrips).getArrivalTimeFunctions();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getArrivalTimeFunctions();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && T>>
    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions(const size_t maxTrips = 100) noexcept {
        return getProfileHandle(maxTrips).getTravelTimeFunctions();
    }

    template<bool T = OneToOne, typename = std::enable_if_t<T == OneToOne && !T>>
    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions(const Vertex vertex, const size_t maxTrips = 100) noexcept {
        return getProfileHandle(vertex, maxTrips).getTravelTimeFunctions();
    }

    inline void reset() noexcept {
        forwardRaptor.reset();
        backwardRaptor.reset();
        if constexpr (OneToOne) {
            ProfileWrapper().swap(targetProfile);
        } else {
            std::vector<ProfileWrapper>(profileByVertex.size()).swap(profileByVertex);
        }
    }

private:
    inline void clear() noexcept {
        if constexpr (OneToOne) {
            targetProfile.clear();
        } else {
            for (ProfileWrapper& profile : profileByVertex) {
                profile.clear();
            }
        }
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
                const int walkingTime = forwardRaptor.getWalkingTravelTime(stop);
                if (walkingTime == intMax) continue;
                const int departureTime = stopEvents[stopEventIndex].departureTime - walkingTime;
                if (departureTime < minDepartureTime) continue;
                if constexpr (CorrectDepartureTimes) {
                    if (departureTime >= maxDepartureTime) continue;
                }
                departures.emplace_back(departureTime, stop);
            }
        }
        sort(departures);
        if constexpr (Debug) std::cout << "Number of departures: " << departures.size() << std::endl;
    }

    template<bool RUN_BACKWARD_QUERIES = false>
    inline void collectArrivals(const int departureTime, const IndexedSet<false, Vertex>& targets, const IndexedSet<false, Vertex>& reachedVertices) noexcept {
        const IndexedSet<false, Vertex>& smallTargetSet = (targets.size() < reachedVertices.size()) ? targets : reachedVertices;
        const IndexedSet<false, Vertex>& largeTargetSet = (targets.size() < reachedVertices.size()) ? reachedVertices : targets;
        for (const Vertex target : smallTargetSet) {
            if (!largeTargetSet.contains(target)) continue;
            collectArrivals<RUN_BACKWARD_QUERIES>(departureTime, target, profileByVertex[target]);
        }
    }

    template<bool RUN_BACKWARD_QUERIES = false>
    inline void collectArrivals(const int departureTime, const Vertex target, ProfileWrapper& profile) noexcept {
        if constexpr (RUN_BACKWARD_QUERIES) {
            backwardRaptor.reset();
            suppressUnusedParameterWarning(departureTime);
            const std::vector<ArrivalLabel>& arrivals = forwardRaptor.getArrivals(target);
            for (size_t i = arrivals.size() - 1; i != size_t(-1); i--) {
                if (arrivals[i].numberOfTrips == 0) continue;
                if constexpr (OneToOne && Debug) std::cout << arrivals[i] << std::endl;
                backwardRaptor.template run<true>(target, -arrivals[i].arrivalTime, sourceVertex, arrivals[i].numberOfTrips);
                const int realDepartureTime = -backwardRaptor.getArrivalTime(sourceVertex, arrivals[i].numberOfTrips);
                Assert(realDepartureTime >= departureTime, "Backward query found an earlier departure time! Original: " << departureTime << ", New: " << realDepartureTime);
                profile.add(arrivals[i], realDepartureTime);
            }
            profile.sort();
        } else {
            for (const ArrivalLabel arrival : forwardRaptor.getArrivals(target)) {
                if (arrival.numberOfTrips == 0) continue;
                if constexpr (OneToOne && Debug) std::cout << arrival << std::endl;
                profile.add(arrival, departureTime);
            }
        }
    }

private:
    RaptorModule forwardRaptor;
    RaptorModule backwardRaptor;
    const Data& data;
    const Data& reverseData;

    Vertex sourceVertex;
    int minDepartureTime;
    int maxDepartureTime;
    std::vector<DepartureLabel> departures;

    //One-to-one only
    ProfileWrapper targetProfile;

    //One-to-many only
    std::vector<ProfileWrapper> profileByVertex;
    IndexedSet<false, Vertex> targetVertices;
};


template<bool DEBUG = false, bool ONE_TO_ONE = false, bool TRIP_PRUNING = true, bool CORRECT_DEPARTURE_TIMES = true>
using RangeTransitiveRAPTOR = RangeRAPTOR<TransitiveRAPTORModule<Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, ONE_TO_ONE, true, TRIP_PRUNING>, ONE_TO_ONE, DEBUG, CORRECT_DEPARTURE_TIMES>;

template<typename INITIAL_TRANSFERS, bool DEBUG = false, bool ONE_TO_ONE = false, bool USE_MIN_TRANSFER_TIMES = false, bool TRIP_PRUNING = true, bool CORRECT_DEPARTURE_TIMES = true>
using RangeDijkstraRAPTOR = RangeRAPTOR<DijkstraRAPTORModule<INITIAL_TRANSFERS, Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, ONE_TO_ONE, true, USE_MIN_TRANSFER_TIMES, TRIP_PRUNING>, ONE_TO_ONE, DEBUG, CORRECT_DEPARTURE_TIMES>;

template<bool DEBUG = false, bool TRIP_PRUNING = true, bool CORRECT_DEPARTURE_TIMES = true>
using RangeULTRARAPTOR = RangeRAPTOR<ULTRARAPTORModule<Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, true, TRIP_PRUNING>, true, DEBUG, CORRECT_DEPARTURE_TIMES>;

}
