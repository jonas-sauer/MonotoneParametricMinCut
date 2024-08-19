#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <algorithm>

#include "../Classes/GraphInterface.h"

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/Vector/Permutation.h"

namespace Graph {
    template<typename GRAPH>
    inline void printInfo(const GRAPH& graph, std::ostream& out = std::cout) noexcept {
        const std::string typeString = graphType(graph);
        const std::string vertexData = cleanGraphType(attributeListToString<typename GRAPH::ListOfRecordVertexAttributes>());
        const std::string edgeData = cleanGraphType(attributeListToString<typename GRAPH::ListOfRecordEdgeAttributes>());
        out << typeString.substr(0, typeString.find_first_of('<')) << " with " << String::prettyInt(graph.numVertices()) << " vertices and " << String::prettyInt(graph.numEdges()) << " edges"
            << " (" << String::bytesToString(graph.byteSize()) << " on disk)." << std::endl;
        if (!vertexData.empty()) out << "    Vertices contain: " << vertexData << "." << std::endl;
        if (!edgeData.empty()) out << "    Edges contain: " << edgeData << "." << std::endl;
        if (vertexData.empty() && edgeData.empty()) out << "      no additional data exists." << std::endl;
    }

    template<typename GRAPH>
    inline void writeStatisticsFile(const GRAPH& graph, const std::string& fileNameBase, const std::string& separator = ".") noexcept {
        const std::string fileName = fileNameBase + separator + "statistics.txt";
        std::ofstream statistics(fileName);
        Assert(statistics, "Cannot create output stream for: " << fileName);
        Assert(statistics.is_open(), "Cannot open output stream for: " << fileName);
        printInfo(graph, statistics);
        graph.printAnalysis(statistics);
        statistics.close();
    }
}
