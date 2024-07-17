#pragma once

#include <array>
#include <vector>

#include "CSA.h"

#include "../UnionFind.h"

#include "../../DataStructures/CSA/Data.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/Console/Progress.h"

namespace CSA {

class ShortestPathCover {

private:
    struct Node {

    public:
        Node(const uint32_t stopId = noStop, const uint32_t parentIndex = 0, const uint32_t subTreeSize = 0) :
            stopId(stopId),
            parentIndex(parentIndex),
            subTreeSize(subTreeSize) {
        }

    public:
        uint32_t stopId;
        uint32_t parentIndex;
        uint32_t subTreeSize;

    };

    struct Tree {

    public:
        Tree(const StopId root, const std::vector<StopId>& parent, std::vector<uint64_t>& stopWeights) :
            stopWeights(stopWeights) {
            std::vector<std::vector<StopId>> children(parent.size());
            size_t treeSize = 1;
            for (StopId stop = StopId(0); stop < parent.size(); stop++) {
                if (parent[stop] == noStop) continue;
                treeSize++;
                children[parent[stop]].emplace_back(stop);
            }
            nodes.reserve(treeSize + 2);
            nodes.emplace_back();
            std::vector<Node> stack;
            stack.emplace_back(root, 0, 1);
            while (!stack.empty()) {
                nodes.emplace_back(stack.back());
                stack.pop_back();
                for (const StopId child : children[nodes.back().stopId]) {
                    stack.emplace_back(child, nodes.size() - 1, 1);
                }
            }
            nodes.emplace_back();
            nodes.shrink_to_fit();
            for (size_t i = nodes.size() - 2; i > 0; i--) {
                nodes[nodes[i].parentIndex].subTreeSize += nodes[i].subTreeSize;
                stopWeights[nodes[i].stopId] += nodes[i].subTreeSize;
            }
        }

        inline void remove(const StopId stop) noexcept {
            const uint32_t beginOfSubTree = findBeginOfSubTree(stop);
            if (beginOfSubTree + 1 == nodes.size()) return;
            const uint32_t subTreeSize = nodes[beginOfSubTree].subTreeSize;
            for (uint32_t i = beginOfSubTree; i > 0; i = nodes[i].parentIndex) {
                stopWeights[nodes[i].stopId] -= subTreeSize;
                nodes[i].subTreeSize -= subTreeSize;
            }
            uint32_t endOfSubTree = beginOfSubTree + 1;
            while (nodes[endOfSubTree].parentIndex >= beginOfSubTree) {
                if (stopWeights[nodes[endOfSubTree].stopId] < nodes[endOfSubTree].subTreeSize) warning("B");
                stopWeights[nodes[endOfSubTree].stopId] -= nodes[endOfSubTree].subTreeSize;
                endOfSubTree++;
            }
            const uint32_t shift = endOfSubTree - beginOfSubTree;
            for (uint32_t i = endOfSubTree; i < nodes.size(); i++) {
                nodes[i - shift].stopId = nodes[i].stopId;
                nodes[i - shift].parentIndex = (nodes[i].parentIndex < beginOfSubTree) ? (nodes[i].parentIndex) : (nodes[i].parentIndex - shift);
                nodes[i - shift].subTreeSize = nodes[i].subTreeSize;
            }
            nodes.resize(nodes.size() - shift);
        }

    private:
        inline uint32_t findBeginOfSubTree(const StopId stop) noexcept {
            nodes.back().stopId = stop;
            uint32_t i = 1;
            while (nodes[i].stopId != stop) {
                i++;
            }
            return i;
        }

    public:
        std::vector<uint64_t>& stopWeights;
        std::vector<Node> nodes;

    };

public:
    ShortestPathCover(const Data& csaData, const uint32_t coverSize, const uint32_t numberOfTrees, const int maxStationTransferTime, const int maxTreeTime, const uint32_t numberOfTreeGroups) :
        stopWeights(numberOfTreeGroups, std::vector<uint64_t>(csaData.numberOfStops(), 0)) {
        UnionFind unionFind(csaData.numberOfStops());
        for (const StopId stop : csaData.stops()) {
            for (const Edge edge : csaData.transferGraph.edgesFrom(stop)) {
                if (csaData.transferGraph.get(TravelTime, edge) > maxStationTransferTime) continue;
                unionFind(stop, csaData.transferGraph.get(ToVertex, edge));
            }
        }
        std::vector<StopId> stationOfStop;
        for (const StopId stop : csaData.stops()) {
            //if (stop < 100) warning(stop, " -> ", unionFind(stop));
            //stationOfStop.emplace_back(StopId(unionFind(stop)));
            stationOfStop.emplace_back(stop);
        }
        std::vector<StopId> parents(csaData.numberOfStops(), noStop);
        CSA<true, false, false> csaQuery(csaData);
        std::cout << "Building shortest path trees... " << std::endl;
        Progress treeProgress(numberOfTrees);
        for (uint32_t i = 0; i < numberOfTrees; i++) {
            const StopId source = stationOfStop[rand() % csaData.numberOfStops()];
            const int time = (rand() % (16 * 60 * 60)) + (4 * 60 * 60);
            csaQuery.run(source, time, noStop, time + maxTreeTime);
            for (const StopId stop : csaData.stops()) {
                parents[stationOfStop[stop]] = (csaQuery.reachable(stop) & stop != source) ? (stationOfStop[csaQuery.getParent(stop)]) : (noStop);
            }
            trees.emplace_back(source, parents, stopWeights[i % numberOfTreeGroups]);
            treeProgress++;
        }
        std::cout << "Selecting shortest path cover stops... " << std::endl;
        Progress coverProgress(coverSize);
        while (cover.size() < coverSize) {
            StopId bestStop = StopId(0);
            uint32_t bestWeight = getWeight<3>(bestStop);
            for (const StopId stop : csaData.stops()) {
                const uint32_t newWeight = getWeight<3>(stop);
                if (bestWeight >= newWeight) continue;
                bestWeight = newWeight;
                bestStop = stop;
            }
            std::cout << "   Selecting stop " << bestStop << ", with weight " << String::prettyInt(bestWeight) << std::endl;
            cover.emplace_back(bestStop);
            for (Tree& tree : trees) {
                tree.remove(bestStop);
            }
            coverProgress++;
        }
    }

    inline const std::vector<StopId>& getCover() const noexcept {
        return cover;
    }

private:
    template<size_t N>
    inline uint32_t getWeight(const StopId stop) const noexcept {
        std::array<uint32_t, N> maxWeights;
        for (size_t i = 0; i < N; i++) {
            maxWeights[i] = 0;
        }
        for (const std::vector<uint64_t>& weights : stopWeights) {
            if (weights[stop] <= maxWeights[0]) continue;
            maxWeights[0] = weights[stop];
            for (size_t i = 1; i < N; i++) {
                if (maxWeights[i - 1] <= maxWeights[i]) break;
                std::swap(maxWeights[i - 1], maxWeights[i]);
            }
        }
        return maxWeights[0];
    }

private:
    std::vector<std::vector<uint64_t>> stopWeights;
    std::vector<Tree> trees;
    std::vector<StopId> cover;

};

}
