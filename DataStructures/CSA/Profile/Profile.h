#pragma once

#include <iterator>
#include <vector>

#include "ProfileEntry.h"

#include "../../../Helpers/Vector/Vector.h"

namespace CSA {

struct Profile {
    Profile() : entries(1) { }

    inline int getEarliestArrivalTime(const int departureTime) const noexcept {
        for (size_t i = entries.size() - 1; i != size_t(-1); i--) {
            if (departureTime <= entries[i].departureTime) {
                return entries[i].arrivalTime;
            }
        }
        AssertMsg(false, "Could not find a breakpoint after " << departureTime << "!");
        return entries[0].arrivalTime;
    }

    inline bool pushEntry(const int departureTime, const int arrivalTime) noexcept {
        AssertMsg(departureTime <= entries.back().departureTime, "Entry does not belong at the end of the profile!");
        if (arrivalTime >= entries.back().arrivalTime) return false;
        if (departureTime == entries.back().departureTime) {
            entries.back().arrivalTime = arrivalTime;
        } else {
            entries.emplace_back(departureTime, arrivalTime);
        }
        return true;
    }

    inline bool addEntry(const int departureTime, const int arrivalTime) noexcept {
        size_t i = entries.size();
        while (i > 0 && entries[i-1].departureTime < departureTime) {
            i--;
        }
        if (arrivalTime >= entries[i-1].arrivalTime) return false;
        std::vector<ProfileEntry> temp;
        temp.insert(temp.end(), std::make_move_iterator(entries.begin() + i), std::make_move_iterator(entries.end()));
        entries.resize(i);
        if (departureTime == entries.back().departureTime) {
            entries.back().arrivalTime = arrivalTime;
        } else {
            entries.emplace_back(departureTime, arrivalTime);
        }
        for (const ProfileEntry& entry : temp) {
            if (entry.arrivalTime < entries.back().arrivalTime) {
                entries.push_back(entry);
            }
        }
        return true;
    }

    inline bool isSorted() const noexcept {
        return Vector::isSorted(entries, [&](const ProfileEntry& a, const ProfileEntry& b) {
            return a.departureTime >= b.departureTime || a.arrivalTime >= b.arrivalTime;
        });
    }

    std::vector<ProfileEntry> entries;
};

struct ProfileVectorWrapper {
    ProfileVectorWrapper(const size_t size) :
        profile(size) {
    }

    inline int getEarliestArrivalTime(const size_t vertex, const int departureTime) const noexcept {
        return profile[vertex].getEarliestArrivalTime(departureTime);
    }

    inline bool pushEntry(const size_t vertex, const int departureTime, const int arrivalTime) noexcept {
        return profile[vertex].pushEntry(departureTime, arrivalTime);
    }

    inline bool addEntry(const size_t vertex, const int departureTime, const int arrivalTime) noexcept {
        return profile[vertex].addEntry(departureTime, arrivalTime);
    }

    inline size_t size(const size_t vertex) const noexcept {
        return profile[vertex].entries.size();
    }

    inline ProfileEntry& entry(const size_t vertex, const size_t i) noexcept {
        return profile[vertex].entries[i];
    }

    inline const ProfileEntry& entry(const size_t vertex, const size_t i) const noexcept {
        return profile[vertex].entries[i];
    }

    inline void clear() noexcept {
        Vector::fill(profile, Profile());
    }

    std::vector<Profile> profile;
};

}
