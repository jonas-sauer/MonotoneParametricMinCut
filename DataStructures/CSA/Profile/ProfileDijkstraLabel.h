#pragma once

#include <deque>
#include <iterator>

#include "ProfileEntry.h"

#include "../../Container/ExternalKHeap.h"

namespace CSA {

//TODO: Implement contiguous variant
struct ProfileDijkstraLabel : public ExternalKHeapElement {
    ProfileDijkstraLabel() : entries(1) {}

    inline bool hasUnsettledEntries() const noexcept {
        return entries.size() > 1;
    }

    inline int key() const noexcept {
        AssertMsg(hasUnsettledEntries(), "Label does not have entries and should not be in queue!");
        return entries[1].departureTime;
    }

    inline ProfileEntry settle() noexcept {
        AssertMsg(hasUnsettledEntries(), "Label does not have entries to settle!");
        entries.pop_front();
        return entries[0];
    }

    inline bool insert(const ProfileEntry& entry) noexcept {
        AssertMsg(entry.departureTime <= entries[0].departureTime, "Entry is inserted out of order! Last departure time: " << entries[0].departureTime << ", New departure time: " << entry.departureTime);
        size_t i = entries.size();
        while (i > 0 && entries[i-1].departureTime < entry.departureTime) {
            i--;
        }
        if (entry.arrivalTime >= entries[i-1].arrivalTime) return false;
        std::vector<ProfileEntry> temp;
        temp.insert(temp.end(), std::make_move_iterator(entries.begin() + i), std::make_move_iterator(entries.end()));
        entries.resize(i);
        if (entry.departureTime == entries.back().departureTime && hasUnsettledEntries()) {
            entries.back().arrivalTime = entry.arrivalTime;
        } else {
            entries.emplace_back(entry.departureTime, entry.arrivalTime);
        }
        for (const ProfileEntry& entry : temp) {
            if (entry.arrivalTime < entries.back().arrivalTime) {
                entries.push_back(entry);
            }
        }
        AssertMsg(hasUnsettledEntries(), "Label does not have entries to settle!");
        return true;
    }

    inline bool hasSmallerKey(const ProfileDijkstraLabel* const other) const noexcept {
        return key() > other->key();
    }

    //TODO: Use dynamic circular buffer?
    std::deque<ProfileEntry> entries;
};

}
