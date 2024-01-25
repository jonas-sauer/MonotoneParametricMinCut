#pragma once

#include <limits>
#include <vector>

#include "Profile.h"

#include "../../Helpers/Assert.h"

namespace Assignment {

class StopLabel {
public:
    StopLabel() : waitingProfile(1), transferProfile(1) {}

    inline void addWaitingEntry(const ProfileEntry& entry) noexcept {
        Assert(!waitingProfile.empty(), "Missing sentinel entry!");
        Assert(entry < waitingProfile.back(), "New entry " << entry << " is dominated by " << waitingProfile.back() << "!");
        if (entry.departureTime == waitingProfile.back().departureTime) {
            waitingProfile.back() = entry;
        } else {
            waitingProfile.emplace_back(entry);
        }
    }

    template<typename PROFILER>
    inline void addTransferEntry(const ProfileEntry& entry, const PROFILER& profiler) noexcept {
        if (transferProfile.size() <= 1) {
            transferProfile.emplace_back(entry);
        } else {
            size_t insertionIndex = transferProfile.size() - 1;
            int shift = -1;
            while (transferProfile[insertionIndex].departureTime < entry.departureTime) {
                Assert(insertionIndex > 0, "Insertion index reached sentinel (sentinel time: " << transferProfile[0].departureTime << ", entry time: " << entry.departureTime << ")!");
                if (entry.patDominates(transferProfile[insertionIndex])) shift++;
                insertionIndex--;
            }
            if (transferProfile[insertionIndex].patDominates(entry)) return;
            if (transferProfile[insertionIndex].departureTime == entry.departureTime) {
                Assert(insertionIndex > 0, "Insertion index reached sentinel (sentinel time: " << transferProfile[0].departureTime << ", entry time: " << entry.departureTime << ")!");
                shift++;
                insertionIndex--;
            }
            if (shift == 0) {
                transferProfile[insertionIndex + 1] = entry;
            } else if (shift == -1) {
                profiler.copyEntry();
                transferProfile.emplace_back(transferProfile.back());
                for (size_t i = transferProfile.size() - 3; i > insertionIndex; i--) {
                    profiler.copyEntry();
                    transferProfile[i + 1] = transferProfile[i];
                }
                transferProfile[insertionIndex + 1] = entry;
            } else {
                transferProfile[insertionIndex + 1] = entry;
                for (size_t i = insertionIndex + 2; i + shift < transferProfile.size(); i++) {
                    profiler.copyEntry();
                    transferProfile[i] = transferProfile[i + shift];
                }
                transferProfile.resize(transferProfile.size() - shift);
            }
        }
        Assert(checkTransferProfile(), "Profile is not monotone!");
    }

    inline PerceivedTime evaluateWithDelay(const int time, const int maxDelay, const double waitingCosts) const noexcept {
        PerceivedTime pat = 0.0;
        double probability = 0.0;
        for (size_t i = transferProfile.size() - 1; i > 0; i--) {
            if (transferProfile[i].departureTime < time) continue;
            const double newProbability = delayProbability(transferProfile[i].departureTime - time, maxDelay);
            Assert((newProbability >= probability && newProbability <= 1.0), "delayProbability (" << newProbability << ") is not a probability! (x: " << (transferProfile[i].departureTime - time) << ", maxDelay: " << maxDelay << ")");
            pat += (newProbability - probability) * transferProfile[i].evaluate(time, waitingCosts);
            Assert(pat < INFTY, "PAT has reached infinity (time: " << time << ", entry.departureTime: " << transferProfile[i].departureTime << ", probability: " << newProbability << ")!");
            probability = newProbability;
            if (probability >= 1) break;
        }
        if (probability < 1.0) {
            pat = (probability > 0.0000001) ? (pat / probability) : (std::numeric_limits<PerceivedTime>::infinity());
        }
        Assert(pat == pat, "PAT calculation failed (result = " << pat << ")!");
        return pat;
    }

    inline const ProfileEntry& getSkipEntry() const noexcept {
        Assert(!waitingProfile.empty(), "Missing sentinel entry!");
        return waitingProfile.back();
    }

    inline const ProfileEntry& getFailureEntry(const int time) const noexcept {
        Assert(!transferProfile.empty(), "Missing sentinel entry!");
        size_t i = transferProfile.size() - 1;
        while (transferProfile[i].departureTime < time) i--;
        return transferProfile[i];
    }

    inline const Profile& getWaitingProfile() const noexcept {
        return waitingProfile;
    }

private:
    inline static double delayProbability(const double time, const double maxDelay) noexcept {
        if (time < 0) return 0.0;
        if (time >= maxDelay) return 1.0;
        return (31.0/30.0) - ((11.0/30.0) * (maxDelay / ((10 * time) + maxDelay)));
    }

    inline bool checkTransferProfile() const noexcept {
        for (size_t i = 0; i < transferProfile.size() - 1; i++) {
            Assert(transferProfile[i + 1] < transferProfile[i], "Profile is not monotone! (" << transferProfile[i + 1] << " >= " << transferProfile[i] << ")");
        }
        return true;
    }

    Profile waitingProfile;
    Profile transferProfile;
};
}
