#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <string>

#include "../Geometry/Point.h"

#include "../../Helpers/Types.h"
#include "../../Helpers/BinarySearch.h"
#include "../../Helpers/String/String.h"
#include "../../Helpers/String/Enumeration.h"

#include "../../DataStructures/Container/ExternalKHeap.h"

template<typename ENTRY_TYPE>
class Profile : public ExternalKHeapElement {

public:
    using EntryType = ENTRY_TYPE;
    using Type = Profile<EntryType>;
    static constexpr size_t EntrySize = EntryType::Size;

    using DepartureType = std::vector<int>;
    using ArrivalType = std::vector<EntryType>;

public:
    Profile() :
        walkingTime(intMax),
        departureTimes(1, intMax),
        arrivalTimes(1, EntryType::Max()),
        lastUnsettledIndex(1) {
    }

    static Type Min(const Type& p1, const Type& p2) noexcept {
        Type result;
        size_t i1 = 1;
        size_t i2 = 1;
        while ((i1 < p1.size()) && (i2 < p2.size())) {
            if (p1.departureTimes[i1] > p2.departureTimes[i2]) {
                result.append(p1.departureTimes[i1], p1.arrivalTimes[i1]);
                i1++;
            } else {
                result.append(p2.departureTimes[i2], p2.arrivalTimes[i2]);
                i2++;
            }
        }
        while (i1 < p1.size()) {
            result.append(p1.departureTimes[i1], p1.arrivalTimes[i1]);
            i1++;
        }
        while (i2 < p2.size()) {
            result.append(p2.departureTimes[i2], p2.arrivalTimes[i2]);
            i2++;
        }
        result.walkingTime = std::min(p1.walkingTime, p2.walkingTime);
        return result;
    }

    inline bool empty() const noexcept {
        return departureTimes.empty();
    }

    inline size_t size() const noexcept {
        return departureTimes.size();
    }

    inline int departureTime(const int i) const noexcept {
        return departureTimes[i];
    }

    inline const EntryType& arrivalTime(const int i) const noexcept {
        return arrivalTimes[i];
    }

    inline void clear() noexcept {
        walkingTime = intMax;
        DepartureType(1, intMax).swap(departureTimes);
        ArrivalType(1, EntryType::Max()).swap(arrivalTimes);
        lastUnsettledIndex = 1;
    }

    inline void reset() noexcept {
        walkingTime = intMax;
        DepartureType(1, intMax).swap(departureTimes);
        ArrivalType(1, EntryType::Max()).swap(arrivalTimes);
        lastUnsettledIndex = 1;
    }

    inline const EntryType& back() const noexcept {
        return arrivalTimes.front();
    }

    inline bool isConnectedToTarget() const noexcept {
        return walkingTime < intMax;
    }

    inline void setWalkingTime(const int time) noexcept {
        walkingTime = time;
    }

    inline int getWalkingTime() const noexcept {
        return walkingTime;
    }

    inline bool hasUnsettledEntries() const noexcept {
        return lastUnsettledIndex < departureTimes.size();
    }

    inline int key() const noexcept {
        Assert(hasUnsettledEntries(), "Profile key is compared, but Profile has no unsettled entry");
        return departureTimes[lastUnsettledIndex];
    }

    inline const EntryType& lastUnsettledEntry() const noexcept {
        Assert(hasUnsettledEntries(), "Unsettled entry is requested, but Profile has no unsettled entry");
        return arrivalTimes[lastUnsettledIndex];
    }

    inline bool hasSmallerKey(const Profile* const other) const noexcept {
        return key() > other->key();
    }

    inline void settle() noexcept {
        lastUnsettledIndex++;
    }

    inline int earliestDeparture() const noexcept {
        return departureTimes.back();
    }

    inline const EntryType& earliestArrival() const noexcept {
        return arrivalTimes.back();
    }

    inline bool isImprovedBy(const int departureTime, const EntryType& entry) const noexcept {
        return entry < evaluateLinear(departureTime);
    }

    inline bool dominates(const int departureTime, const EntryType& entry) const noexcept {
        return !isImprovedBy(departureTime, entry);
    }

    inline void insert(const int departureTime, const EntryType& arrivalTime) noexcept {
        Assert(isImprovedBy(departureTime, arrivalTime), "Cannot insert an entry that does not improve the profile!");
        insertIfImprovement(departureTime, arrivalTime);
    }

    inline bool insertIfImprovement(const int departureTime, const EntryType& arrivalTime) noexcept {
        for (size_t i = departureTimes.size() - 1; i < departureTimes.size(); i--) {
            if (departureTimes[i] < departureTime) continue;
            if (!(arrivalTime < arrivalTimes[i])) return false;
            if (departureTimes[i] == departureTime) {
                for (size_t j = departureTimes.size() - 1; j >= i; j--) {
                    arrivalTimes[j].minimize(arrivalTime);
                }
            } else {
                i++;
                if (i == departureTimes.size()) {
                    departureTimes.emplace_back(departureTime);
                    arrivalTimes.emplace_back(arrivalTime);
                    arrivalTimes[i].minimize(arrivalTimes[i - 1]);
                } else {
                    departureTimes.emplace_back(departureTimes.back());
                    arrivalTimes.emplace_back(arrivalTimes.back());
                    for (size_t j = departureTimes.size() - 3; j >= i; j--) {
                        arrivalTimes[j + 2].minimize(arrivalTime);
                        departureTimes[j + 1] = departureTimes[j];
                        arrivalTimes[j + 1] = arrivalTimes[j];
                    }
                    arrivalTimes[i + 1].minimize(arrivalTime);
                    departureTimes[i] = departureTime;
                    arrivalTimes[i] = arrivalTime;
                    arrivalTimes[i].minimize(arrivalTimes[i - 1]);
                }
            }
            if (lastUnsettledIndex > i) {
                lastUnsettledIndex = i;
            }
            return true;
        }
        Assert(false, "Departure time " << departureTime << " did not evaluate the sentinel entry!");
        return false;
    }

    inline void append(const int departureTime, const EntryType& arrivalTime) noexcept {
        Assert(departureTimes.back() >= departureTime, "Adding profile entries in increasing order (last Entry: " << String::secToTime(departureTimes.back()) << ", new Entry: " << String::secToTime(departureTime) << ")!");
        if (!(arrivalTime < arrivalTimes.back())) return;
        if (departureTimes.back() == departureTime) {
            arrivalTimes.back().minimize(arrivalTime);
        } else {
            departureTimes.emplace_back(departureTime);
            arrivalTimes.emplace_back(arrivalTime);
            arrivalTimes.back().minimize(arrivalTimes[arrivalTimes.size() - 2]);
        }
        if (lastUnsettledIndex > departureTimes.size() - 1) {
            lastUnsettledIndex = departureTimes.size() - 1;
        }
    }

    inline const EntryType& evaluateBinary(int departureTime) const {
        const int i = lastTrueIndex(departureTimes.size(), [&](const int i){return departureTimes[i] >= departureTime;}, 4);
        Assert(i >= 0, "Departure time " << departureTime << " did not evaluate the sentinel entry!");
        Assert(i < arrivalTimes.size(), "Index i = " << i << " is out of bounds (" << 0 << ", " << arrivalTimes.size() << ")!");
        Assert(departureTimes[i] >= departureTime, "Departure time " << departureTime << " evaluated at " << departureTimes[i] << "!");
        return arrivalTimes[i];
    }

    inline const EntryType& evaluateLinear(int departureTime) const {
        for (int i = departureTimes.size() - 1; i >= 0; i--) {
            if (departureTimes[i] >= departureTime) return arrivalTimes[i];
        }
        Assert(false, "Departure time " << departureTime << " did not evaluate the sentinel entry!");
        return arrivalTimes.front();
    }

    inline long long byteSize() const noexcept {
        long long result = 2 * sizeof(int);
        result += Vector::byteSize(departureTimes);
        result += sizeof(arrivalTimes);
        for (const EntryType& entry : arrivalTimes) {
            result += entry.byteSize();
        }
        result += sizeof(int) * (departureTimes.capacity() - departureTimes.size());
        result += sizeof(EntryType) * (arrivalTimes.capacity() - arrivalTimes.size());
        return result;
    }

    friend inline std::ostream& operator<<(std::ostream& out, const Type& p) noexcept {
        out << "Walking time: " << String::secToString(p.walkingTime) << " (" << p.walkingTime << ")" << std::endl;
        for (int i = p.departureTimes.size() - 1; i > 0; i--) {
            out << "   " << std::setw(15) << String::secToString(p.departureTimes[i]) << " " << std::setw(7) << p.departureTimes[i] << "  |  ";
            out << p.arrivalTimes[i] << std::endl;
        }
        return out;
    }

private:
    int walkingTime;
    DepartureType departureTimes;
    ArrivalType arrivalTimes;

    size_t lastUnsettledIndex;

};
