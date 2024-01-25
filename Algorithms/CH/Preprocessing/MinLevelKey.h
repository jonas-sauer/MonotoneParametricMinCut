#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "CHData.h"
#include "KeyFunction.h"
#include "../../PathFinder.h"

namespace CH {

template<typename WITNESS_SEARCH, typename KEY_FUNCTION>
class MinLevelKey {

public:
    using WitnessSearch = WITNESS_SEARCH;
    using KeyFunction = KEY_FUNCTION;
    using Type = MinLevelKey<WitnessSearch, KeyFunction>;

    using KeyType = typename KeyFunction::KeyType;
    constexpr static int offset = (1 << 30) - 1;

public:
    MinLevelKey(const std::vector<uint16_t>& minLevel, const KeyFunction& keyFunction = KeyFunction()) :
        data(nullptr),
        keyFunction(keyFunction),
        minLevel(minLevel) {
    }

    inline KeyType operator() (Vertex vertex) {
        int key = keyFunction(vertex);
        if (data->level[vertex] < minLevel[vertex]) {
            key = std::min(key, offset);
            key += offset;
        }
        return key;
    }

    template<typename T> inline void update(T&) noexcept {}

    inline void initialize(const Data* data, WitnessSearch* witnessSearch) noexcept {
        this->data = data;
        keyFunction.initialize(data, witnessSearch);
    }

private:
    const Data* data;
    KeyFunction keyFunction;
    std::vector<uint16_t> minLevel;

};

}
