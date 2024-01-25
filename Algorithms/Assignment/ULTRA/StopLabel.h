#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "../../../Helpers/BinarySearch.h"
#include "../../../Helpers/Vector/Vector.h"

namespace Assignment::ULTRA {

class StopLabel {

public:
    using Type = StopLabel;

    struct ProfileEntry {
        ProfileEntry() :
            departureTime(Unreachable),
            normalizedPAT(Unreachable) {
        }

        ProfileEntry(const int departureTime, const PerceivedTime normalizedPAT) :
            departureTime(departureTime),
            normalizedPAT(normalizedPAT) {
        }

        inline friend std::ostream& operator<<(std::ostream& out, const ProfileEntry& p) {
            return out << "(" << p.departureTime << ", " << p.normalizedPAT << ")";
        }

        int departureTime;
        PerceivedTime normalizedPAT;
    };

    using Profile = std::vector<ProfileEntry>;

public:
    StopLabel() :
        walkingPTT(Unreachable),
        minimumPTT(Unreachable),
        profile(1),
        index(0),
        transferProfile(1) {
    }

    inline void addEntry(const int departureTime, const PerceivedTime pat, const int waitingTime, const double waitingCosts) noexcept {
        Assert(!profile.empty(), "The profile is not valid!");
        Assert(departureTime <= profile.back().departureTime, "New entry (" << departureTime << ") departs later than " << profile.back() << "!");
        if (profile.back().departureTime == departureTime) {
            profile.back().normalizedPAT = pat + (departureTime * waitingCosts);
        } else {
            profile.emplace_back(departureTime, pat + (departureTime * waitingCosts));
        }
        minimumPTT = std::min(minimumPTT, pat - departureTime);
        addEntry(departureTime, pat, waitingTime, waitingCosts, 0, 0);
    }

    inline void addEntry(const int departureTime, const PerceivedTime pat, const int waitingTime, const double waitingCosts, const int transferTime, const double transferCosts) noexcept {
        addProfileEntry(ProfileEntry(departureTime - waitingTime - transferTime, pat + (departureTime * waitingCosts) + (transferTime * transferCosts) - (transferTime * waitingCosts)));
    }

    inline void addProfileEntry(const ProfileEntry& entry) noexcept {
        if (transferProfile.size() <= 1) {
            transferProfile.emplace_back(entry);
        } else {
            size_t insertionIndex = transferProfile.size() - 1;
            int shift = -1;
            while (transferProfile[insertionIndex].departureTime < entry.departureTime) {
                Assert(insertionIndex > 0, "Insertion index reached sentinel (sentinel time: " << transferProfile[0].departureTime << ", entry time: " << entry.departureTime << ")!");
                if (transferProfile[insertionIndex].normalizedPAT >= entry.normalizedPAT) shift++;
                insertionIndex--;
            }
            if (transferProfile[insertionIndex].normalizedPAT <= entry.normalizedPAT) return;
            if (transferProfile[insertionIndex].departureTime == entry.departureTime) {
                Assert(insertionIndex > 0, "Insertion index reached sentinel (sentinel time: " << transferProfile[0].departureTime << ", entry time: " << entry.departureTime << ")!");
                shift++;
                insertionIndex--;
            }
            if (shift == 0) {
                transferProfile[insertionIndex + 1] = entry;
            } else if (shift == -1) {
                transferProfile.emplace_back(transferProfile.back());
                for (size_t i = transferProfile.size() - 3; i > insertionIndex; i--) {
                    transferProfile[i + 1] = transferProfile[i];
                }
                transferProfile[insertionIndex + 1] = entry;
            } else {
                transferProfile[insertionIndex + 1] = entry;
                for (size_t i = insertionIndex + 2; i + shift < transferProfile.size(); i++) {
                    transferProfile[i] = transferProfile[i + shift];
                }
                transferProfile.resize(transferProfile.size() - shift);
            }
        }
    }

    inline void setWalkingPTT(const int newWalkingPTT) noexcept {
        walkingPTT = newWalkingPTT;
        minimumPTT = std::min<int>(minimumPTT, newWalkingPTT);
    }

    inline PerceivedTime evaluateWithDelay(const int time, const int maxDelay, const double waitingCosts) const noexcept {
        PerceivedTime pat = 0.0;
        double probability = 0.0;
        for (size_t i = transferProfile.size() - 1; i > 0; i--) {
            if (transferProfile[i].departureTime < time) continue;
            const double newProbability = delayProbability(transferProfile[i].departureTime - time, maxDelay);
            Assert((newProbability >= probability && newProbability <= 1.0), "delayProbability (" << newProbability << ") is not a probability! (x: " << (transferProfile[i].departureTime - time) << ", maxDelay: " << maxDelay << ")");
            pat += (newProbability - probability) * (transferProfile[i].normalizedPAT - (time * waitingCosts));
            Assert(pat < Unreachable, "PAT has reached infinity (time: " << time << ", transferProfile[i].departureTime: " << transferProfile[i].departureTime << ", probability: " << delayProbability(transferProfile[i].departureTime - time, maxDelay) << ")!");
            probability = newProbability;
            if (probability >= 1) break;
        }
        if (probability < 1.0) {
            pat = (probability > 0.0000001) ? (pat / probability) : (Unreachable);
        }
        Assert(pat == pat, "PAT calculation failed (result = " << pat << ")!");
        return pat;
    }

    inline PerceivedTime evaluateWithoutDelay(const int time, const double waitingCosts) const noexcept {
        Assert(profile.back().departureTime >= time, "Profile contains entries that are too early!");
        if (profile.size() <= 1) {
            return Unreachable;
        } else {
            return profile.back().normalizedPAT - (profile.back().departureTime * waitingCosts);
        }
    }

    inline static double delayProbability(const double time, const double maxDelay) noexcept {
        if (time < 0) return 0.0;
        if (time >= maxDelay) return 1.0;
        return (31.0/30.0) - ((11.0/30.0) * (maxDelay / ((10.0 * time) + maxDelay)));
    }

    inline void findIndex(const int time) noexcept {
        Assert(index < profile.size(), "Index = " << index  << " is out of bounds (0, " << profile.size() << ")!");
        if (profile.back().departureTime >= time) {
            index = profile.size() - 1;
        } else if (profile[index].departureTime >= time) {
            do {
                index++;
                Assert(index < profile.size(), "Index = " << index  << " is out of bounds (0, " << profile.size() << ")!");
            } while (profile[index].departureTime >= time);
            index--;
        } else {
            do {
                index--;
                Assert(index < profile.size(), "Index = " << index  << " is out of bounds (0, " << profile.size() << ")!");
            } while (profile[index].departureTime < time);
        }
        Assert(index < profile.size(), "Index = " << index  << " is out of bounds (0, " << profile.size() << ")!");
        Assert(profile[index].departureTime >= time, "Departure time is " << time << " but profile was evaluated for " << profile[index].departureTime << "!");
        Assert((index == profile.size() - 1) || (profile[index + 1].departureTime < time), "Departure time is " << time << " but profile was not evaluated for " << profile[index + 1].departureTime << "!");
    }

    inline void findIndexFast(const int time) noexcept {
        index = firstFalseIndex(0, profile.size() - 1, [&](const int i){return profile[i].departureTime >= time;}, 16) - 1;
        Assert(index < profile.size(), "Index = " << index  << " is out of bounds (0, " << profile.size() << ")!");
        Assert(profile[index].departureTime >= time, "Departure time is " << time << " but profile was evaluated for " << profile[index].departureTime << "!");
        Assert((index == profile.size() - 1) || (profile[index + 1].departureTime < time), "Departure time is " << time << " but profile was not evaluated for " << profile[index + 1].departureTime << "!");
    }

    inline double evaluateIndex(const size_t j, const int time, const double waitingCosts) noexcept {
        Assert(j < profile.size(), "Index = " << j  << " is out of bounds (0, " << profile.size() << ")!");
        Assert(profile[j].departureTime >= time, "Departure time is " << time << " but profile was evaluated for " << profile[j].departureTime << "!");
        Assert((j == profile.size() - 1) || (profile[j + 1].departureTime < time), "Departure time is " << time << " but profile was not evaluated for " << profile[j + 1].departureTime << "!");
        if (index == 0) {
            return Unreachable;
        } else {
            return profile[j].normalizedPAT - (time * waitingCosts);
        }
    }

    inline double evaluateWithWaitingCosts(const int time, const double waitingCosts) noexcept {
        findIndex(time);
        return evaluateIndex(index, time, waitingCosts);
    }

    inline double evaluateWithoutWaiting(int& time, const double waitingCosts) noexcept {
        findIndex(time);
        time = profile[index].departureTime;
        return evaluateIndex(index, time, waitingCosts);
    }

    inline void resetIndex() noexcept {
        index = profile.size() - 1;
    }

    inline PerceivedTime getWalkingPTT() const noexcept {
        return walkingPTT;
    }

    inline PerceivedTime getMinimumPTT() const noexcept {
        return minimumPTT;
    }

    inline friend std::ostream& operator<<(std::ostream& out, const StopLabel& s) {
        return out << s.profile;
    }

private:
    PerceivedTime walkingPTT;
    PerceivedTime minimumPTT;

    Profile profile;
    size_t index;

    Profile transferProfile;
};

}
