#pragma once

#include <vector>
#include <algorithm>

#include "../Helpers/Assert.h"

class UnionFind {

public:
    UnionFind(const int n) : parent(n, n), n(n) {}

    inline int find(const int i) noexcept {
        if (parent[i] >= n) {
            return i;
        } else {
            parent[i] = find(parent[i]);
            return parent[i];
        }
    }
    inline int operator()(const int i) noexcept {
        return find(i);
    }

    inline void unite(const int i, const int j) noexcept {
        if (find(i) != find(j)) {
            link(find(i), find(j));
        }
    }
    inline void operator()(const int i, const int j) noexcept {
        unite(i, j);
    }

protected:
    inline void link(const int i, const int j) noexcept {
        Assert(parent[i] >= n);
        Assert(parent[j] >= n);
        Assert(i != j);
        if (parent[i] < parent[j]) {
            parent[i] = j;
        } else if (parent[j] < parent[i]) {
            parent[j] = i;
        } else {
            parent[i] = j;
            parent[j]++;
        }
    }

    std::vector<int> parent;
    int n;

};