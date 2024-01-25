#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <strings.h>

#include "CHData.h"
#include "KeyFunction.h"
#include "../../PathFinder.h"

#include "../../../Helpers/Helpers.h"

namespace CH {

template<typename WITNESS_SEARCH, typename KEY_FUNCTION>
class PathAwareKey {

public:
    using WitnessSearch = WITNESS_SEARCH;
    using KeyFunction = KEY_FUNCTION;
    using Type = PathAwareKey<WitnessSearch, KeyFunction>;

    using KeyType = int;

public:
    PathAwareKey(const KeyFunction& keyFunction = KeyFunction()) : data(0), paths(0), keyFunction(keyFunction) {}
    ~PathAwareKey() {delete paths;}

    inline KeyType operator() (Vertex vertex) {
        if (paths->isOnPath(vertex)) {
            return leastSignificantBit(paths->getIndex(vertex) + 1) - 100000;
        } else {
            return keyFunction(vertex);
        }
    }

    template<typename T> inline void update(T&) {}

    inline void initialize(const Data* data, WitnessSearch* witnessSearch) {
        keyFunction.initialize(data, witnessSearch);
        this->data = data;
        delete paths;
        paths = new PathFinder<CHCoreGraph>(data->core);
        paths->run();
    }

private:
    KeyFunction keyFunction;
    const Data* data;

    PathFinder<CHCoreGraph>* paths;

};

}
