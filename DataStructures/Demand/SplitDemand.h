#pragma once

#include <vector>

#include "../Container/Set.h"
#include "../CSA/Data.h"
#include "../../Helpers/Types.h"

template<typename DEMAND_TYPE>
class SplitDemand {
public:
    SplitDemand(const Construct::SplitByOriginTag, const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const std::vector<DEMAND_TYPE>& demand, const bool allowDepartureStops) :
        SplitDemand(data, reverseGraph, demand, allowDepartureStops, [&](const DEMAND_TYPE& entry) {
        return entry.originVertex;
    }) {
    }

    SplitDemand(const Construct::SplitByDestinationTag, const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const std::vector<DEMAND_TYPE>& demand, const bool allowDepartureStops) :
        SplitDemand(data, reverseGraph, demand, allowDepartureStops, [&](const DEMAND_TYPE& entry) {
        return entry.destinationVertex;
    }) {
    }

    inline size_t size() const noexcept {
        return verticesWithDemand.size();
    }

    inline Vertex vertexAtIndex(const size_t i) const noexcept {
        return verticesWithDemand[i];
    }

    inline std::vector<DEMAND_TYPE>& operator[](const Vertex vertex) noexcept {
        return entries[vertex];
    }

    inline const std::vector<DEMAND_TYPE>& operator[](const Vertex vertex) const noexcept {
        return entries[vertex];
    }

private:
    template<typename SPLIT_VERTEX>
    SplitDemand(const CSA::Data& data, const CSA::TransferGraph& reverseGraph, const std::vector<DEMAND_TYPE>& demand, const bool allowDepartureStops, const SPLIT_VERTEX& splitVertex) {
        IndexedSet<false, Vertex> set(data.transferGraph.numVertices());
        for (const DEMAND_TYPE& entry : demand) {
            if (entry.originVertex == entry.destinationVertex) continue;
            if (!allowDepartureStops && data.isStop(entry.originVertex)) continue;
            if (!data.isStop(entry.originVertex) && data.transferGraph.outDegree(entry.originVertex) == 0) continue;
            if (!data.isStop(entry.destinationVertex) && reverseGraph.outDegree(entry.destinationVertex) == 0) continue;

            const Vertex vertex = splitVertex(entry);
            set.insert(vertex);
            if (vertex >= entries.size()) entries.resize(vertex + 1);
            entries[vertex].emplace_back(entry);
        }

        verticesWithDemand = set.getValues();
    }

    std::vector<std::vector<DEMAND_TYPE>> entries;
    std::vector<Vertex> verticesWithDemand;
};
