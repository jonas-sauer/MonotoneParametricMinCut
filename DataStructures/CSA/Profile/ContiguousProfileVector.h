#pragma once

#include <iterator>
#include <utility>
#include <vector>

#include "ProfileEntry.h"

namespace CSA {

struct ContiguousProfileVector {
    ContiguousProfileVector(const size_t size) :
        profileBegin(size + 1),
        predecessor(size),
        successor(size),
        lastEntry(size - 1),
        profileSize(size, 1),
        entries(size) {
        for (size_t i = 0; i < size; i++) {
            profileBegin[i] = i;
            predecessor[i] = i-1;
            successor[i] = i+1;
        }
        profileBegin[size] = size;
    }

    inline int getEarliestArrivalTime(const size_t vertex, const int departureTime) const noexcept {
        for (size_t i = profileEnd(vertex) - 1; i != size_t(profileBegin[vertex] - 1); i--) {
            if (departureTime <= entries[i].departureTime) {
                return entries[i].arrivalTime;
            }
        }
        Assert(false, "Could not find a breakpoint after " << departureTime << "!");
        return entry(vertex, 0).arrivalTime;
    }

    inline bool pushEntry(const size_t vertex, const int departureTime, const int arrivalTime) noexcept {
        ProfileEntry& lastEntry = entries[profileEnd(vertex) - 1];
        Assert(departureTime <= lastEntry.departureTime, "Entry does not belong at the end of the profile!");
        if (arrivalTime >= lastEntry.arrivalTime) return false;
        if (departureTime == lastEntry.departureTime) {
            lastEntry.arrivalTime = arrivalTime;
        } else {
            if (profileCapacity(vertex) == profileSize[vertex]) {
                moveToEnd(vertex, 2*profileSize[vertex] + 2);
            }
            entries[profileEnd(vertex)].departureTime = departureTime;
            entries[profileEnd(vertex)].arrivalTime = arrivalTime;
            profileSize[vertex]++;
        }
        return true;
    }

    inline bool addEntry(const size_t vertex, const int departureTime, const int arrivalTime) noexcept {
        size_t i = profileEnd(vertex);
        while (i > profileBegin[vertex] && entries[i-1].departureTime < departureTime) {
            i--;
        }
        if (arrivalTime >= entries[i-1].arrivalTime) return false;
        std::vector<ProfileEntry> temp;
        temp.insert(temp.end(), std::make_move_iterator(entries.begin() + i), std::make_move_iterator(end(vertex)));
        const size_t oldSize = profileSize[vertex];
        profileSize[vertex] = i - profileBegin[vertex];
        if (departureTime == entries[i-1].departureTime) {
            entries[i-1].arrivalTime = arrivalTime;
        } else {
            if (profileCapacity(vertex) == oldSize) {
                moveToEnd(vertex, 2*oldSize + 2);
            }
            entries[profileEnd(vertex)].departureTime = departureTime;
            entries[profileEnd(vertex)].arrivalTime = arrivalTime;
            profileSize[vertex]++;
        }
        for (const ProfileEntry& entry : temp) {
            if (entry.arrivalTime < entries[profileEnd(vertex) - 1].arrivalTime) {
                entries[profileEnd(vertex)] = entry;
                profileSize[vertex]++;
            }
        }
        return true;
    }

    inline void moveToEnd(const size_t vertex, const size_t newCapacity) noexcept {
        const size_t newBegin = entries.size();
        entries.resize(entries.size() + newCapacity);
        std::move(begin(vertex), end(vertex), entries.begin() + newBegin);
        profileBegin[vertex] = newBegin;
        if (successor[vertex] != entries.size()) {
            predecessor[successor[vertex]] = predecessor[vertex];
        }
        if (predecessor[vertex] != size_t(-1)) {
            successor[predecessor[vertex]] = successor[vertex];
        }
        successor[lastEntry] = vertex;
        predecessor[vertex] = lastEntry;
        lastEntry = vertex;
    }

    inline std::vector<ProfileEntry>::iterator begin(const size_t vertex) noexcept {
        return entries.begin() + profileBegin[vertex];
    }

    inline std::vector<ProfileEntry>::const_iterator begin(const size_t vertex) const noexcept {
        return entries.begin() + profileBegin[vertex];
    }

    inline std::vector<ProfileEntry>::iterator end(const size_t vertex) noexcept {
        return begin(vertex) + profileSize[vertex];
    }

    inline std::vector<ProfileEntry>::const_iterator end(const size_t vertex) const noexcept {
        return begin(vertex) + profileSize[vertex];
    }

    inline size_t profileEnd(const size_t vertex) const noexcept {
        return profileBegin[vertex] + profileSize[vertex];
    }

    inline size_t profileCapacity(const size_t vertex) const noexcept {
        return profileBegin[successor[vertex]] - profileBegin[vertex];
    }

    inline size_t size(const size_t vertex) const noexcept {
        return profileSize[vertex];
    }

    inline ProfileEntry& entry(const size_t vertex, const size_t i) noexcept {
        return entries[profileBegin[vertex] + i];
    }

    inline const ProfileEntry& entry(const size_t vertex, const size_t i) const noexcept {
        return entries[profileBegin[vertex] + i];
    }

    inline void clear() noexcept {
        for (size_t i = 0; i < predecessor.size(); i++) {
            profileBegin[i] = i;
            predecessor[i] = i-1;
            successor[i] = i+1;
        }
        profileBegin.back() = profileBegin.size();
        lastEntry = predecessor.size() - 1;
        Vector::fill(profileSize, size_t(1));
        std::vector<ProfileEntry>(predecessor.size()).swap(entries);
    }

    std::vector<size_t> profileBegin;
    std::vector<size_t> predecessor;
    std::vector<size_t> successor;
    size_t lastEntry;
    std::vector<size_t> profileSize;
    std::vector<ProfileEntry> entries;
};

}
