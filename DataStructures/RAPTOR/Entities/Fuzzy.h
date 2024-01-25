#pragma once

#include <array>
#include <vector>

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/Helpers.h"

#include "ArrivalLabel.h"

namespace RAPTOR {

template<typename LABEL>
struct ScoredLabel {
    using Label = LABEL;
    using Type = ScoredLabel<LABEL>;

    ScoredLabel(const Label label = Label(), const double score = 0.0) :
        label(label),
        score(score) {
    }

    inline bool operator==(const Type& other) const noexcept {
        return label == other.label;
    }

    inline bool operator<(const Type& other) const noexcept {
        return score > other.score;
    }

    inline friend std::ostream& operator<<(std::ostream& out, const Type& scoredLabel) noexcept {
        return out << scoredLabel.label << ", score: " << scoredLabel.score;
    }

    Label label;
    double score;
};

using ScoredWalkingLabel = ScoredLabel<WalkingParetoLabel>;

struct FuzzyValues {
    FuzzyValues(const double lt = 0, const double eq = 0, const double gt = 0) :
        values{lt, eq, gt} {
        }

    inline FuzzyValues operator+(const FuzzyValues& other) const noexcept {
        return FuzzyValues(values[0] + other.values[0], values[1] + other.values[1], values[2] + other.values[2]);
    }

    inline FuzzyValues& operator+=(const FuzzyValues& other) noexcept {
        values[0] += other.values[0];
        values[1] += other.values[1];
        values[2] += other.values[2];
        return *this;
    }

    inline double lt() const noexcept {
        return values[0];
    }

    inline double eq() const noexcept {
        return values[1];
    }

    inline double gt() const noexcept {
        return values[2];
    }

    std::array<double, 3> values;
};

struct FuzzyComparison {
public:
    FuzzyComparison() : scalingFactor(0) {}

    //x = eps, y = chi
    FuzzyComparison(const int x, const double y) :
        scalingFactor(std::log(y)/(x * x)) {
    }

    inline FuzzyValues operator()(const int valueA, const int valueB) const noexcept {
        const int value = valueA - valueB;
        if (value == 0) return FuzzyValues(0, 1, 0);
        const double eq = std::exp(scalingFactor * value * value);
        if (value < 0) return FuzzyValues(1 - eq, eq, 0);
        return FuzzyValues(0, eq, 1 - eq);
    }

private:
    double scalingFactor;
};

template<typename LABEL>
class FuzzyEvaluation {
public:
    using Label = LABEL;
    using LabelWithScore = ScoredLabel<Label>;
    inline static constexpr int M = Label::NumberOfCriteria;

    FuzzyEvaluation(const int x[M], const double y[M]) {
        for (size_t i = 0; i < M; i++) {
            comparison[i] = FuzzyComparison(x[i], y[i]);
        }
    }

    inline double getDegreeOfDomination(const Label& a, const Label& b) const noexcept {
        const FuzzyValues counts = getCriteriaCounts(a, b);
        const double lt = counts.lt();
        const double eq = counts.eq();
        const double gt = counts.gt();
        AssertMsg(std::abs(lt + eq + gt - M) < 0.001, "Fuzzy sum (" << lt + eq + gt << ") does not match number of criteria (" << M << ")!");
        if (eq == M) return 0; //To simplify fuzzy scoring, exactly equal labels don't dominate each other
        if (lt < gt) return 0;
        return 1 - gt/lt;
    }

    inline double getScore(const Label& label, const std::vector<Label>& labels) const noexcept {
        double maxDomination = 0;
        for (const Label& other : labels) {
            maxDomination = std::max(maxDomination, getDegreeOfDomination(other, label));
        }
        return 1 - maxDomination;
    }

    inline LabelWithScore scoreLabel(const Label& label, const std::vector<Label>& labels) const noexcept {
        return LabelWithScore(label, getScore(label, labels));
    }

    inline std::vector<LabelWithScore> scoreLabels(const std::vector<Label>& labels) const noexcept {
        std::vector<LabelWithScore> scoredLabels;
        for (const Label& label : labels) {
            scoredLabels.emplace_back(scoreLabel(label, labels));
        }
        sort(scoredLabels);
        return scoredLabels;
    }

    inline std::vector<LabelWithScore> getTopKLabels(const std::vector<Label>& labels, const size_t k) const noexcept {
        std::vector<LabelWithScore> scoredLabels = scoreLabels(labels);
        if (scoredLabels.size() > k) scoredLabels.resize(k);
        return scoredLabels;
    }

    inline double getMaxSimilarity(const LabelWithScore& referenceLabel, const std::vector<LabelWithScore>& labels) const noexcept {
        double similarity = 0;
        for (const LabelWithScore& label : labels) {
            similarity = std::max(similarity, computeSimilarity(referenceLabel.label, label.label));
        }
        return similarity;
    }

    inline double getMaxSimilarity(const Label& referenceLabel, const std::vector<Label>& labels) const noexcept {
        double similarity = 0;
        for (const Label& label : labels) {
            similarity = std::max(similarity, computeSimilarity(referenceLabel, label));
        }
        return similarity;
    }

    inline double getMaxLEQSimilarity(const std::vector<Label>& labels, const Label& referenceLabel) const noexcept {
        double similarity = 0;
        for (const Label& label : labels) {
            similarity = std::max(similarity, computeLEQSimilarity(label, referenceLabel));
        }
        return similarity;
    }

private:
    inline double computeSimilarity(const Label& a, const Label& b) const noexcept {
        double result = INFTY;
        for (size_t i = 0; i < M; i++) {
            result = std::min(result, getComparison(a, b, i).eq());
        }
        return result;
    }

    inline double computeLEQSimilarity(const Label& a, const Label& b) const noexcept {
        double result = INFTY;
        for (size_t i = 0; i < M; i++) {
            result = std::min(result, 1 - getComparison(a, b, i).gt());
        }
        return result;
    }

    inline FuzzyValues getCriteriaCounts(const Label& a, const Label& b) const noexcept {
        FuzzyValues result;
        for (size_t i = 0; i < M; i++) {
            result += getComparison(a, b, i);
        }
        return result;
    }

    inline FuzzyValues getComparison(const Label& a, const Label& b, size_t i) const noexcept {
        return comparison[i](a.getCriterion(i), b.getCriterion(i));
    }

    FuzzyComparison comparison[M];
};

}
