#pragma once

#include <vector>

#include "../Demand/AccumulatedVertexDemand.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Vector/Permutation.h"

namespace Assignment {

template<typename OPTION_TYPE>
struct ParameterizedChoiceSet {
    using OptionType = OPTION_TYPE;

    ParameterizedChoiceSet() {}

    ParameterizedChoiceSet(const OptionType option, const int departureTime, const int pat) :
        options(1, option),
        departureTimes(1, departureTime),
        pats(1, pat) {
    }

    inline void addChoice(const OptionType option, const int departureTime, const int pat) noexcept {
        options.push_back(option);
        departureTimes.push_back(departureTime);
        pats.push_back(pat);
    }

    inline void sort() noexcept {
        Order order(Construct::Sort, departureTimes);
        order.order(options);
        order.order(departureTimes);
        order.order(pats);
    }

    inline bool empty() const noexcept {
        Assert(options.size() == pats.size(), "options and pats have different size!");
        Assert(departureTimes.size() == pats.size(), "departureTimes and pats have different size!");
        return pats.empty();
    }

    inline size_t size() const noexcept {
        Assert(options.size() == pats.size(), "options and pats have different size!");
        Assert(departureTimes.size() == pats.size(), "departureTimes and pats have different size!");
        return pats.size();
    }

    inline std::vector<int> rooftopDistribution(const AccumulatedVertexDemand::Entry& demandEntry, const int adaptationCost) noexcept {
        std::vector<int> result = rooftopRawDistribution(demandEntry, adaptationCost);
        result.emplace_back(Vector::sum(result));
        return result;
    }

    inline std::vector<int> rooftopCumulativeDistribution(const AccumulatedVertexDemand::Entry& demandEntry, const int adaptationCost) noexcept {
        std::vector<int> result = rooftopRawDistribution(demandEntry, adaptationCost);
        for (size_t i = 1; size(); i++) {
            result[i] += result[i - 1];
        }
        return result;
    }

private:
    inline std::vector<size_t> getOptimalRooftopChoices(const int adaptationCost) noexcept {
        sort();
        size_t bestPrevious = 0;
        std::vector<bool> isDominated(size(), false);
        for (size_t i = 1; i < size(); i++) {
            const int prevPAT = pats[bestPrevious] + adaptationCost * (departureTimes[i] - departureTimes[bestPrevious]);
            if (prevPAT <= pats[i]) {
                isDominated[i] = true;
            } else {
                bestPrevious = i;
            }
        }

        size_t bestNext = size() - 1;
        while (isDominated[bestNext]) bestNext--;
        for (size_t i = bestNext - 1; i != size_t(-1); i--) {
            if (isDominated[i]) continue;
            const int postPAT = pats[bestNext] + adaptationCost * (departureTimes[bestNext] - departureTimes[i]);
            if (postPAT <= pats[i]) {
                isDominated[i] = true;
            } else {
                bestNext = i;
            }
        }

        std::vector<size_t> result;
        for (const size_t i : indices(*this)) {
            if (!isDominated[i]) result.push_back(i);
        }
        return result;
    }

    inline std::vector<int> rooftopRawDistribution(const AccumulatedVertexDemand::Entry& demandEntry, const int adaptationCost) noexcept {
        std::vector<int> distribution(size(), 0);
        std::vector<size_t> relevantChoices = getOptimalRooftopChoices(adaptationCost);
        distribution[relevantChoices[0]] += departureTimes[relevantChoices[0]] - demandEntry.earliestDepartureTime;
        distribution[relevantChoices.back()] += demandEntry.latestDepartureTime - departureTimes[relevantChoices.back()];

        for (size_t i = 1; i < relevantChoices.size(); i++) {
            const size_t currentChoice = relevantChoices[i];
            const size_t previousChoice = relevantChoices[i - 1];
            const int adaptationPenalty = adaptationCost * (departureTimes[currentChoice] - departureTimes[previousChoice]);
            distribution[currentChoice] += (pats[previousChoice] + adaptationPenalty - pats[currentChoice]) / (2 * adaptationCost);
            distribution[previousChoice] += (pats[currentChoice] + adaptationPenalty - pats[previousChoice]) / (2 * adaptationCost);
        }
        return distribution;
    }

public:
    std::vector<OptionType> options;
    std::vector<int> departureTimes;
    std::vector<int> pats;
};

using ChoiceSet = ParameterizedChoiceSet<StopId>;

struct SequentialChoiceSet {

    SequentialChoiceSet() {}

    inline void addChoice(const ConnectionId option, const PerceivedTime pat) noexcept {
        options.push_back(option);
        pats.push_back(pat);
    }

    inline bool empty() const noexcept {
        Assert(options.size() == pats.size(), "options and pats have different size!");
        return pats.empty();
    }

    inline size_t size() const noexcept {
        Assert(options.size() == pats.size(), "options and pats have different size!");
        return pats.size();
    }

    std::vector<ConnectionId> options;
    std::vector<PerceivedTime> pats;

};

}
