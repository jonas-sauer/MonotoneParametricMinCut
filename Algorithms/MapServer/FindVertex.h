#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "Algorithm.h"

#include "../../Helpers/Types.h"
#include "../../DataStructures/Parameter.h"
#include "../../DataStructures/Result.h"

#include "../../Helpers/String/String.h"

class FindVertex : public Algorithm {

public:
    FindVertex() : Algorithm() {
        addParameter(Parameter("size_t", "vertex", "Vertex ID", "0", true));
    }

    virtual std::vector<Result> run(const Vertex, const Vertex) noexcept {
        Vertex vertex = Vertex(getParameter<size_t>("vertex"));
        std::string label = "Vertex " + vertex;
        Result result(label);
        result.color = KIT::green;
        result.nodes.push_back(Result::Node(vertex, "V", KIT::green));
        return std::vector<Result>({result});
    }

    virtual std::vector<Result> runSourceOnly(const Vertex) noexcept {
        return run(noVertex, noVertex);
    }

    virtual std::vector<Result> runTargetOnly(const Vertex) noexcept {
        return run(noVertex, noVertex);
    }

};
