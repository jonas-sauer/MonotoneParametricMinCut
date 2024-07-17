#pragma once

#include <iostream>
#include <vector>

#include "../../Helpers/Types.h"

namespace Assignment {

using GroupId = size_t;
using GroupList = std::vector<GroupId>;

struct DummyGroupList {
    inline void emplace_back(const int) const noexcept {}
};

struct GroupData {
    GroupData(const GroupId groupId = -1, const size_t demandIndex = -1, const double groupSize = 0.0) :
        groupId(groupId),
        demandIndex(demandIndex),
        groupSize(groupSize) {
    }

    friend std::ostream& operator<<(std::ostream& out, const GroupData& g) {
        return out << "Assignment::GroupData{" << g.groupId << ", " << g.demandIndex  << ", " << g.groupSize << "}";
    }

    GroupId groupId;
    size_t demandIndex;
    double groupSize;
};

struct GroupArrivalLabel {
    GroupArrivalLabel() : ids(), arrivalTime(INFTY) {}
    GroupArrivalLabel(const int id, const int arrivalTime = INFTY) : ids(1, id), arrivalTime(arrivalTime) {}
    GroupArrivalLabel(GroupList& groupList, const int arrivalTime) : ids(), arrivalTime(arrivalTime) {ids.swap(groupList);}

    inline bool operator<(const GroupArrivalLabel& other) const noexcept {
        return arrivalTime < other.arrivalTime;
    }

    inline friend void swap(GroupArrivalLabel& a, GroupArrivalLabel& b) noexcept {
        std::swap(a.ids, b.ids);
        std::swap(a.arrivalTime, b.arrivalTime);
    }

    GroupList ids;
    int arrivalTime;
};

}
