#pragma once

#include <vector>

#include "GroupData.h"
#include "../Container/Heap.h"
#include "../CSA/Data.h"
#include "../../Helpers/Types.h"

namespace Assignment {

struct ConnectionBasedGroupTrackingData {
    ConnectionBasedGroupTrackingData(const size_t numberOfConnections, size_t numberOfLinks) :
        groupsWalkingOnLink(numberOfLinks),
        groupsWaitingAtConnection(numberOfConnections),
        groupsInTrip(numberOfConnections) {
    }

    inline void validate() const noexcept {
        for (const GroupList& groupList : groupsWalkingOnLink) {
            AssertMsg(groupList.empty(), "There are groups walking to a stop from the last iteration!");
        }
        for (const GroupList& groupList : groupsWaitingAtConnection) {
            AssertMsg(groupList.empty(), "There are groups waiting at a stop from the last iteration!");
        }
        for (const GroupList& groupList : groupsInTrip) {
            AssertMsg(groupList.empty(), "There are groups in a trip from the last iteration!");
        }
    }

    std::vector<GroupList> groupsWalkingOnLink;
    std::vector<GroupList> groupsWaitingAtConnection;
    std::vector<GroupList> groupsInTrip;
    DummyGroupList groupsAtTarget;
};

struct SequentialGroupTrackingData {
    SequentialGroupTrackingData(const size_t numberOfConnections, const size_t numberOfTrips) :
        groupsWaitingAtConnection(numberOfConnections),
        groupsInTrip(numberOfTrips) {
    }

    inline void validate() const noexcept {
        for (const GroupList& groupList : groupsWaitingAtConnection) {
            AssertMsg(groupList.empty(), "There are groups waiting at a stop from the last iteration!");
        }
        for (const GroupList& groupList : groupsInTrip) {
            AssertMsg(groupList.empty(), "There are groups in a trip from the last iteration!");
        }
    }

    std::vector<GroupList> groupsWaitingAtConnection;
    std::vector<GroupList> groupsInTrip;
    DummyGroupList groupsAtTarget;
};

struct GroupTrackingData {
    GroupTrackingData(const size_t numberOfStops, const size_t numberOfTrips) :
        groupsOriginatingAtStop(numberOfStops),
        groupsWalkingToStop(numberOfStops),
        groupsWaitingAtStop(numberOfStops),
        groupsInTrip(numberOfTrips) {
    }

    inline void validate() const noexcept {
        for (const GroupList& groupList : groupsInTrip) {
            AssertMsg(groupList.empty(), "There are groups in a trip from the last iteration!");
        }
        for (size_t stop = 0; stop < groupsOriginatingAtStop.size(); stop++) {
            AssertMsg(groupsOriginatingAtStop[stop].empty(), "There are groups originating at a stop from the last iteration!");
            AssertMsg(groupsWalkingToStop[stop].empty(), "There are groups walking to a stop from the last iteration!");
            AssertMsg(groupsWaitingAtStop[stop].empty(), "There are groups waiting at a stop from the last iteration!");
        }
    }

    inline void processOriginatingGroups(const CSA::Connection& connection) noexcept {
        std::vector<GroupArrivalLabel>& originatingGroups = groupsOriginatingAtStop[connection.departureStopId];
        GroupList& waitingGroups = groupsWaitingAtStop[connection.departureStopId];
        size_t i = originatingGroups.size() - 1;
        size_t newGroups = 0;
        while (i != size_t(-1) && originatingGroups[i].arrivalTime <= connection.departureTime) {
            AssertMsg(!originatingGroups[i].ids.empty(), "There is an empty set of passengers originating!");
            newGroups += originatingGroups[i].ids.size();
            i--;
        }
        waitingGroups.reserve(waitingGroups.size() + newGroups);
        for (size_t j = originatingGroups.size() - 1; j != i; j--) {
            waitingGroups.insert(waitingGroups.end(), originatingGroups[j].ids.begin(), originatingGroups[j].ids.end());
        }
        originatingGroups.resize(i + 1);
    }

    inline void processWalkingGroups(const CSA::Connection& connection) noexcept {
        Heap<GroupArrivalLabel>& walkingGroups = groupsWalkingToStop[connection.departureStopId];
        GroupList& waitingGroups = groupsWaitingAtStop[connection.departureStopId];
        std::vector<GroupArrivalLabel> removedLabels;
        size_t newGroups = 0;
        while ((!walkingGroups.empty()) && (walkingGroups.min().arrivalTime <= connection.departureTime)) {
            AssertMsg(!walkingGroups.min().ids.empty(), "There is an empty set of passengers walking!");
            removedLabels.emplace_back(std::move(walkingGroups.pop_min()));
            newGroups += removedLabels.back().ids.size();
        }
        waitingGroups.reserve(waitingGroups.size() + newGroups);
        for (const GroupArrivalLabel& label : removedLabels) {
            waitingGroups.insert(waitingGroups.end(), label.ids.begin(), label.ids.end());
        }
    }

    std::vector<std::vector<GroupArrivalLabel>> groupsOriginatingAtStop;
    std::vector<Heap<GroupArrivalLabel>> groupsWalkingToStop;
    std::vector<GroupList> groupsWaitingAtStop;
    std::vector<GroupList> groupsInTrip;
    DummyGroupList groupsAtTarget;
};

}
