#pragma once

#include <chrono>

class Timer {

public:
    using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;
    using Microseconds = std::chrono::duration<double, std::micro>;
    using Milliseconds = std::chrono::duration<double, std::milli>;

    Timer() : start(timestamp()) {}

    inline void restart() noexcept {
        start = timestamp();
    }

    [[nodiscard]] inline double elapsedMicroseconds() const noexcept {
        return std::chrono::duration_cast<Microseconds>(timestamp() - start).count();
    }

    [[nodiscard]] inline double elapsedMilliseconds() const noexcept {
        return std::chrono::duration_cast<Milliseconds>(timestamp() - start).count();
    }

private:
    inline static TimePoint timestamp() noexcept {
        return std::chrono::steady_clock::now();
    }

private:
    TimePoint start;

};

