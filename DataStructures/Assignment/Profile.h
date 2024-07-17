#pragma once

#include <iostream>
#include <vector>

#include "../../Helpers/Assert.h"
#include "../../Helpers/Types.h"

namespace Assignment {

struct ProfileEntry {
public:
    ProfileEntry() :
        departureTime(INFTY),
        connectionId(noConnection),
        normalizedPAT(Unreachable) {
    }

    ProfileEntry(const int departureTime, const ConnectionId connectionId, const PerceivedTime pat, const double waitingCosts) :
        departureTime(departureTime),
        connectionId(connectionId),
        normalizedPAT(pat + (departureTime * waitingCosts)) {
    }

    ProfileEntry(const int departureTime, const ConnectionId connectionId, const PerceivedTime originalPAT, const int transferTime, const int waitingTime, const double walkingCosts, const double waitingCosts) :
        departureTime(departureTime - transferTime - waitingTime),
        connectionId(connectionId),
        normalizedPAT(originalPAT + ((departureTime - transferTime) * waitingCosts) + (transferTime * walkingCosts)) {
    }

    inline bool operator<(const ProfileEntry& other) const noexcept {
        return (departureTime <= other.departureTime) && (normalizedPAT < other.normalizedPAT);
    }

    inline bool patDominates(const ProfileEntry& other) const noexcept {
        return normalizedPAT <= other.normalizedPAT;
    }

    inline PerceivedTime evaluate(const int time, const double waitingCosts) const noexcept {
        Assert(time <= departureTime, "Evaluation time lies after departureTime!");
        if (normalizedPAT == Unreachable) return Unreachable;
        return normalizedPAT - time * waitingCosts;
    }

    inline void print(const double waitingCosts) const noexcept {
        std::cout << "(" << departureTime << ", " << evaluate(departureTime, waitingCosts) << ", " << connectionId << ")" << std::endl;
    }

    inline friend std::ostream& operator<<(std::ostream& out, const ProfileEntry& p) {
        return out << "(" << p.departureTime << ", " << p.normalizedPAT << ", " << p.connectionId << ")";
    }

    int departureTime;
    ConnectionId connectionId;

private:
    PerceivedTime normalizedPAT;
};

using Profile = std::vector<ProfileEntry>;

}
