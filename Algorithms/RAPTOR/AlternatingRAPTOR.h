#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "../../Helpers/Helpers.h"
#include "../../Helpers/Meta.h"

#include "../../DataStructures/Container/Heap.h"
#include "../../DataStructures/RAPTOR/Entities/ArrivalLabel.h"
#include "../../DataStructures/RAPTOR/Entities/Profile.h"
#include "../../DataStructures/Geometry/Rectangle.h"
#include "DijkstraRAPTOR.h"
#include "Profiler.h"
#include "RAPTOR.h"
#include "ULTRARAPTOR.h"

namespace RAPTOR {

template<typename RAPTOR_MODULE, bool DEBUG = false>
class AlternatingRAPTOR {

public:
    using RaptorModule = RAPTOR_MODULE;
    static constexpr bool Debug = DEBUG;
    using InitialTransferGraph = typename RaptorModule::InitialTransferGraph;
    using Type = AlternatingRAPTOR<RaptorModule, Debug>;
    using SourceType = typename RaptorModule::SourceType;

public:
    template<typename ATTRIBUTE>
    AlternatingRAPTOR(const Data& data, const Data& reverseData, const InitialTransferGraph& forwardGraph, const InitialTransferGraph& backwardGraph, const ATTRIBUTE weight) :
        data(data),
        reverseData(reverseData),
        forwardRaptor(data, forwardGraph, backwardGraph, weight),
        backwardRaptor(reverseData, backwardGraph, forwardGraph, weight),
        sourceVertex(SourceType()),
        targetVertex(SourceType()),
        minDepartureTime(never),
        maxDepartureTime(never),
        currentArrivalTime(-never) {
    }

    template<typename T = CHGraph, typename = std::enable_if_t<Meta::Equals<T, CHGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    AlternatingRAPTOR(const Data& data, const Data& reverseData, const CH::CH& chData) :
        AlternatingRAPTOR(data, reverseData, chData.forward, chData.backward, Weight) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    AlternatingRAPTOR(const Data& data, const Data& reverseData, const TransferGraph& forwardGraph, const TransferGraph& backwardGraph) :
        AlternatingRAPTOR(data, reverseData, forwardGraph, backwardGraph, TravelTime) {
    }

    template<typename T = TransferGraph, typename = std::enable_if_t<Meta::Equals<T, TransferGraph>() && Meta::Equals<T, InitialTransferGraph>()>>
    AlternatingRAPTOR(const Data& data, const Data& reverseData) :
        AlternatingRAPTOR(data, reverseData, data.transferGraph, reverseData.transferGraph) {
    }

    inline void run(const SourceType source, const SourceType target, const int minTime = 0, const int maxTime = 60 * 60 * 24) noexcept {
        clear();
        sourceVertex = source;
        targetVertex = target;
        minDepartureTime = minTime;
        maxDepartureTime = maxTime;
        currentArrivalTime = minTime;
        runForward(minTime);
        while (!queue.empty()) {
            ArrivalLabel label = queue.pop_min();
            Assert(currentArrivalTime < label.arrivalTime, "Queue key is decreasing: old key = " << currentArrivalTime << ", new key = " << label.arrivalTime << "!");
            currentArrivalTime = label.arrivalTime;
            if constexpr (Debug) {
                std::cout << "Current queue label: arrivalTime = " << String::secToTime(label.arrivalTime) << " (" << label.arrivalTime << "), numberOfTrips = " << label.numberOfTrips << std::endl;
                std::cout << "   New arrival time..." << std::endl;
            }
            runBackward(label);
            if constexpr (Debug) std::cout << "   departureTime = " << String::secToTime(-backwardRaptor.getArrivalTime(sourceVertex, label.numberOfTrips)) << " (" << -backwardRaptor.getArrivalTime(sourceVertex, label.numberOfTrips) << ")" << std::endl;
            addProfileEntry(label, -backwardRaptor.getArrivalTime(sourceVertex, label.numberOfTrips));
            while ((!queue.empty()) && (queue.min().arrivalTime == label.arrivalTime)) {
                if constexpr (Debug) std::cout << "   Repeated arrival time..." << std::endl;
                label = queue.pop_min();
                addProfileEntry(label, -backwardRaptor.getArrivalTime(sourceVertex, label.numberOfTrips));
            }
        }
        pruneProfileWithWalking();
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

    inline ProfileHandle getProfileHandle(const size_t maxTrips = 100) const noexcept {
        return ProfileHandle(getProfile(maxTrips), minDepartureTime, maxDepartureTime, forwardRaptor.getWalkingTravelTime());
    }

    inline std::vector<Geometry::Point> getArrivalTimeFunction(const size_t maxTrips = 100) const noexcept {
        return getProfileHandle(maxTrips).getArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getTravelTimeFunction(const size_t maxTrips = 100) const noexcept {
        return getProfileHandle(maxTrips).getTravelTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingArrivalTimeFunction() const noexcept {
        return getProfileHandle().getWalkingArrivalTimeFunction();
    }

    inline std::vector<Geometry::Point> getWalkingTravelTimeFunction() const noexcept {
        return getProfileHandle().getWalkingTravelTimeFunction();
    }

    inline std::vector<std::vector<Geometry::Point>> getArrivalTimeFunctions(const size_t maxTrips = 100) const noexcept {
        return getProfileHandle(maxTrips).getArrivalTimeFunctions();
    }

    inline std::vector<std::vector<Geometry::Point>> getTravelTimeFunctions(const size_t maxTrips = 100) const noexcept {
        return getProfileHandle(maxTrips).getTravelTimeFunctions();
    }

    inline void reset() noexcept {
        clear<true>();
    }

private:
    template<bool RESET_CAPACITIES = false>
    inline void clear() noexcept {
        forwardRaptor.template clear<RESET_CAPACITIES>();
        backwardRaptor.template clear<RESET_CAPACITIES>();
        queue.clear();
        if constexpr (RESET_CAPACITIES) {
            Profile().swap(profile);
            std::vector<int>().swap(forwardDepartures);
        } else {
            profile.clear();
            forwardDepartures.clear();
        }
    }

    inline void addProfileEntry(const ArrivalLabel& label, const int departureTime) noexcept {
        while ((!queue.empty()) && (queue.min() == label)) {
            queue.pop_min();
        }
        if constexpr (Debug) std::cout << "   NEW profile entry: departure time = " << String::secToTime(departureTime) << " (" << departureTime << "), arrival time = " << String::secToTime(label.arrivalTime) << " (" << label.arrivalTime << "), number of trips = " << label.numberOfTrips << std::endl;
        profile.emplace_back(label, departureTime);
        runForward(departureTime + 1);
    }

    inline void runForward(const int departureTime) noexcept {
        for (int i = forwardDepartures.size() - 1; i >= 0; i--) {
            if (forwardDepartures[i] == departureTime) return;
        }
        if (departureTime > maxDepartureTime) return;
        if constexpr (Debug) std::cout << "Running forward query, departureTime = " << String::secToTime(departureTime) << " (" << departureTime << ")" << std::endl;
        forwardRaptor.run(sourceVertex, departureTime, targetVertex);
        std::vector<ArrivalLabel> labels = forwardRaptor.getArrivals();
        for (int i = labels.size() - 1; i >= 0; i--) {
            if (labels[i].arrivalTime < currentArrivalTime) {
                if constexpr (Debug) std::cout << "   NOT adding label to queue: arrivalTime = " << String::secToTime(labels[i].arrivalTime) << " (" << labels[i].arrivalTime << "), numberOfTrips = " << labels[i].numberOfTrips << std::endl;
                continue;
            } else {
                if constexpr (Debug) std::cout << "   Adding label to queue: arrivalTime = " << String::secToTime(labels[i].arrivalTime) << " (" << labels[i].arrivalTime << "), numberOfTrips = " << labels[i].numberOfTrips << std::endl;
                queue.push_back(labels[i]);
            }
        }
        forwardDepartures.emplace_back(departureTime);
    }

    inline void runBackward(const ArrivalLabel& label) noexcept {
        if constexpr (Debug) std::cout << "Running backward query, arrivalTime = " << String::secToTime(label.arrivalTime) << " (" << label.arrivalTime << "), max Round = " << ((label.numberOfTrips + 1) * 2) << std::endl;
        backwardRaptor.run(targetVertex, -label.arrivalTime, sourceVertex);
    }

    inline void pruneProfileWithWalking() noexcept {
        const int walkingTravelTime = forwardRaptor.getWalkingTravelTime();
        size_t numberOfRemovedEntries = 0;
        for (size_t i = 0; i < profile.size(); i++) {
            if ((profile[i].travelTime() < walkingTravelTime) && (profile[i].arrivalTime - maxDepartureTime < walkingTravelTime)) {
                if (numberOfRemovedEntries > 0) {
                    profile[i - numberOfRemovedEntries] = profile[i];
                }
            } else {
                numberOfRemovedEntries++;
            }
        }
        if (numberOfRemovedEntries > 0) {
            profile.resize(profile.size() - numberOfRemovedEntries);
        }
    }

private:
    const Data& data;
    const Data& reverseData;
    RaptorModule forwardRaptor;
    RaptorModule backwardRaptor;

    SourceType sourceVertex;
    SourceType targetVertex;
    int minDepartureTime;
    int maxDepartureTime;

    Profile profile;
    Heap<ArrivalLabel> queue;

    int currentArrivalTime;
    std::vector<int> forwardDepartures;

};

template<bool DEBUG = false, bool USE_MIN_TRANSFER_TIMES = false>
using AlternatingTransitiveRAPTOR = AlternatingRAPTOR<RAPTOR<true, Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, true, USE_MIN_TRANSFER_TIMES, true>, DEBUG>;

template<typename INITIAL_TRANSFERS, bool DEBUG = false, bool USE_MIN_TRANSFER_TIMES = false>
using AlternatingDijkstraRAPTOR = AlternatingRAPTOR<DijkstraRAPTOR<INITIAL_TRANSFERS, Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, true, USE_MIN_TRANSFER_TIMES, true>, DEBUG>;

template<bool DEBUG = false>
using AlternatingULTRARAPTOR = AlternatingRAPTOR<ULTRARAPTOR<Meta::IF<DEBUG, SimpleProfiler, NoProfiler>, true>, DEBUG>;

}
