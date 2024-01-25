#pragma once

namespace Assignment {

struct BoardingData {
    double boardingDemand{0.0};
    double boardingCapacity{0.0};

    inline void reset() noexcept {
        boardingDemand = 0.0;
        boardingCapacity = 0.0;
    }

    inline double boardingProbability() const noexcept {
        if (boardingDemand <= boardingCapacity) return 1.0;
        return boardingCapacity/boardingDemand;
    }
};

struct ConnectionLoadData {
public:
    ConnectionLoadData(const double capacity) : capacity(capacity), load(0) {}

    const double capacity;
    double load;
    BoardingData boardingData;

    inline void clear() noexcept {
        load = 0.0;
        boardingData.reset();
    }

    inline bool isFull() const noexcept {
        return load == capacity;
    }

    inline double relativeLoad() const noexcept {
        return load/capacity;
    }

    inline double boardingProbability() const noexcept {
        return boardingData.boardingProbability();
    }
};

}
