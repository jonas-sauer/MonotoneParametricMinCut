#pragma once

#include <limits>

#include "TaggedInteger.h"

using Vertex = TaggedInteger<0, u_int32_t, -u_int32_t(1)>;
constexpr Vertex noVertex(Vertex::InvalidValue);

using Edge = TaggedInteger<1, u_int32_t, -u_int32_t(1)>;
constexpr Edge noEdge(Edge::InvalidValue);

inline constexpr int intMax = std::numeric_limits<int>::max();
inline constexpr int INFTY = std::numeric_limits<int>::max() / 2;
inline constexpr double doubleMax = std::numeric_limits<double>::max();

inline constexpr int FORWARD = 0;
inline constexpr int BACKWARD = 1;

struct NO_OPERATION {
    template<typename... ARGS>
    constexpr inline bool operator() (ARGS...) const noexcept {return false;}
};

NO_OPERATION NoOperation;
