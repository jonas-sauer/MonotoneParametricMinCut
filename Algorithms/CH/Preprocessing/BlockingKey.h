#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "CHData.h"
#include "WitnessSearch.h"

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/Ranges/Range.h"
#include "../../../Helpers/Vector/Permutation.h"

#include "../../../DataStructures/Container/Set.h"

namespace CH {

template<typename WITNESS_SEARCH, typename KEY_FUNCTION = GreedyKey<WITNESS_SEARCH>>
class BlockingKey {

public:
    using WitnessSearch = WITNESS_SEARCH;
    using KeyFunction = KEY_FUNCTION;
    using KeyType = typename KeyFunction::KeyType;
    using Type = PartialKey<WitnessSearch, KeyFunction>;

public:
    BlockingKey(const double blockingFactor, const size_t finalCore, const KeyFunction& keyFunction = KeyFunction()) :
        data(nullptr),
        keyFunction(keyFunction),
        blocking(false),
        initialized(false),
        blockingFactor(blockingFactor),
        finalCore(finalCore) {
    }

    inline KeyType operator() (const Vertex vertex) noexcept {
        if (blocking) {
            blockedVertices.insert(vertex);
            return intMax;
        } else {
            return keyFunction(vertex);
        }
    }

    template<typename T>
    inline void update(T& t) noexcept {
        if (initialized) {
            if (blockedVertices.size() > (data->coreSize() * blockingFactor)) {
                blocking = false;
                for (Vertex vertex : blockedVertices) {
                    t.reKey(vertex);
                }
                blockedVertices.clear();
                if (data->coreSize() > finalCore) {
                    blocking = true;
                }
            }
        } else {
            blocking = true;
            initialized = true;
        }
    }

    inline void initialize(const Data* data, WitnessSearch* witnessSearch) noexcept {
        this->data = data;
        keyFunction.initialize(data, witnessSearch);
        blockedVertices = IndexedSet<false, Vertex>(data->numVertices);
        blocking = false;
        initialized = false;
    }

private:
    const Data* data;
    KeyFunction keyFunction;

    IndexedSet<false, Vertex> blockedVertices;
    bool blocking;
    bool initialized;

    double blockingFactor;
    size_t finalCore;

};

}
