#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"
#include "../../DataStructures/RAPTOR/Data.h"

class IsolatedStops : public Algorithm {

public:
    IsolatedStops(const RAPTOR::Data& data) : Algorithm(), data(data) {}

    virtual std::vector<Result> run(const Vertex, const Vertex) noexcept {
        Result result;
        result.color = KIT::green;
        size_t isolatedStops = 0;
        for (const StopId stop : data.stops()) {
            if (data.transferGraph.outDegree(stop) == 0) {
                isolatedStops++;
                result.nodes.push_back(Result::Node(stop, "S", KIT::blue));
            }
            result.info = "Isolated Stops: " + std::to_string(isolatedStops);
        }
        return std::vector<Result>({result});
    }

    virtual std::vector<Result> runSourceOnly(const Vertex) noexcept {
        return run(noVertex, noVertex);
    }

    virtual std::vector<Result> runTargetOnly(const Vertex) noexcept {
        return run(noVertex, noVertex);
    }

private:
    const RAPTOR::Data& data;

};

