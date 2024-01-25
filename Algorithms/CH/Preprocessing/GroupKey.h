#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "CHData.h"
#include "KeyFunction.h"
#include "../../PathFinder.h"

namespace CH {

template<typename WITNESS_SEARCH, typename KEY_FUNCTION>
class GroupKey {

public:
    using WitnessSearch = WITNESS_SEARCH;
    using KeyFunction = KEY_FUNCTION;
    using Type = GroupKey<WitnessSearch, KeyFunction>;

private:
    struct GroupAndKeyData {
        GroupAndKeyData(const int32_t group, const int32_t key) :
            key(key),
            group(group) {
        }
        int32_t key;
        int32_t group;
    };

public:
    union KeyType {
        KeyType(const int32_t group, const int32_t key) : data(group, key) {value += key;}
        inline bool operator<(const KeyType& other) const noexcept {return value < other.value;}
        inline bool operator==(const KeyType& other) const noexcept {return value == other.value;}
        GroupAndKeyData data;
        int64_t value;
    };

public:
    GroupKey(const std::vector<int>& group, const KeyFunction& keyFunction = KeyFunction()) :
        group(group),
        keyFunction(keyFunction) {
        assert(KeyType(-1, 1) < KeyType(1, -1));
    }
    GroupKey(const int numVertices, const KeyFunction& keyFunction = KeyFunction()) :
        group(numVertices, 0),
        keyFunction(keyFunction) {
        assert(KeyType(-1, 1) < KeyType(1, -1));
    }

    inline KeyType operator() (const Vertex vertex) noexcept {
        return KeyType(group[vertex], keyFunction(vertex));
    }

    template<typename T> inline void update(T&) noexcept {}

    inline void initialize(const Data* data, WitnessSearch* witnessSearch) {
        keyFunction.initialize(data, witnessSearch);
    }

    inline int& getGroup(const Vertex vertex) noexcept {return group[vertex];}
    inline int getGroup(const Vertex vertex) const noexcept {return group[vertex];}
    inline void setGroup(const Vertex vertex, const int g) const noexcept {group[vertex] = g;}

private:
    std::vector<int> group;
    KeyFunction keyFunction;

};

}
