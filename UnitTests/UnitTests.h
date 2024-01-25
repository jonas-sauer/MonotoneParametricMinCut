#pragma once

#include "../Helpers/Assert.h"
#include "../Helpers/Debug.h"
#include "../Helpers/HighlightText.h"

namespace UnitTests {

int testCount = 0;
int failureCount = 0;

template<typename... MSG>
inline void check(const bool condition, const MSG&... msg) noexcept {
    testCount++;
    if (!condition) {
        error(msg...);
        printStackTrace();
        std::cout << std::endl;
        failureCount++;
    }
}

inline void evaluate() noexcept {
    if (failureCount == 0) {
        green(testCount, " out of ", testCount, " unit tests were successful!\n");
    } else {
        red(failureCount, " out of ", testCount, " unit tests failed!\n");
    }
}

}
