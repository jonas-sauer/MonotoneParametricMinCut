#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"

template<typename ALGORITHM, bool TIME = false>
class ShortestPath : public Algorithm {

public:
    ShortestPath(ALGORITHM& algorithm, const double secondsPerTimeUnit = 1) :
        Algorithm(),
        algorithm(algorithm),
        secondsPerTimeUnit(secondsPerTimeUnit) {
    }

    virtual std::vector<Result> run(const Vertex sourceVertexId, const Vertex targetVertexId) noexcept {
        algorithm.run(sourceVertexId, targetVertexId);
        const std::vector<Vertex> path(algorithm.getPath(targetVertexId));
        Result result;
        if (!path.empty()) {
            if constexpr (TIME) {
                result.info = "Distance: " + String::secToString(algorithm.getDistance(targetVertexId) * secondsPerTimeUnit);
            } else {
                result.info = "Distance: " + std::to_string(algorithm.getDistance(targetVertexId));
            }
        } else {
            result.info = "Target not reachable.";
        }
        result.color = KIT::red;
        result.nodes.push_back(Result::Node(sourceVertexId, "Source", KIT::blue));
        result.nodes.push_back(Result::Node(targetVertexId, "Target", KIT::blue));
        result.pathes.push_back(Result::Path(path, "A Path", KIT::red));
        return std::vector<Result>({result});
    }

private:
    ALGORITHM& algorithm;
    const double secondsPerTimeUnit;

};
