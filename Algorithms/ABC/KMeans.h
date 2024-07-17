#pragma once

#include <random>

#include "../Dijkstra/Dijkstra.h"

#include "../../Helpers/Types.h"
#include "../../DataStructures/CSA/Data.h"
#include "../../DataStructures/Graph/Graph.h"
#include "../../DataStructures/Container/ExternalKHeap.h"

namespace ABC {

class KMeans {

private:
    static inline constexpr uint32_t undefined = -1;

    struct VertexLabel : public ExternalKHeapElement {
        VertexLabel() :
            ExternalKHeapElement(),
            distance(intMax),
            timeStamp(0) {
        }
        inline bool hasSmallerKey(const VertexLabel* other) const {
            return distance < other->distance;
        }
        int distance;
        int timeStamp;
    };

    struct ColumnLabel : public ExternalKHeapElement {
        ColumnLabel(const size_t numberOfVertices) :
            ExternalKHeapElement(),
            labelOfVertex(numberOfVertices),
            vertexQueue(numberOfVertices),
            weight(0) {
        }
        inline void setSource(const Vertex vertex, const int timeStamp) noexcept {
            weight = 0;
            vertexQueue.clear();
            checkTimeStamp(vertex, timeStamp);
            labelOfVertex[vertex].distance = 0;
            update(vertex);
        }
        inline void update(const Vertex vertex) noexcept {
            vertexQueue.update(&(labelOfVertex[vertex]));
        }
        inline void checkTimeStamp(const Vertex vertex, const int timeStamp) noexcept {
            if (labelOfVertex[vertex].timeStamp < timeStamp) {
                labelOfVertex[vertex].timeStamp = timeStamp;
                labelOfVertex[vertex].distance = intMax;
            }
        }
        inline bool hasSmallerKey(const ColumnLabel* other) const {
            return weight < other->weight;
        }
        std::vector<VertexLabel> labelOfVertex;
        ExternalKHeap<2, VertexLabel> vertexQueue;
        double weight;
    };

public:
    KMeans(const CSA::Data& csa, const uint32_t k, const uint32_t seed, const int iterations = 1, const int distColumns = 1) :
        graph(csa.minTravelTimeGraph()),
        vertexWeight(graph.numVertices(), 0),
        columnOfVertex(graph.numVertices(), undefined),
        labelOfColumn(k, ColumnLabel(graph.numVertices())),
        columnQueue(k),
        timeStamp(0),
        randomGenerator(seed),
        randomDistribution(0, graph.numVertices() - 1) {
        for (const CSA::Connection connection : csa.connections) {
            vertexWeight[connection.departureStopId] += 1;
        }
        assignColumns(randomCenters(k));
        // assignColumns(distCenters(k));
        // assignColumns(weightCenters(k));
        for (int i = 0; i < iterations; i++) {
            assignColumns(geoCenters());
        }
        for (int i = 0; i < distColumns; i++) {
            assignColumnsByDist(geoCenters());
        }
    }

private:
    inline std::vector<Vertex> randomCenters(const uint32_t k) noexcept {
        std::vector<bool> picked(graph.numVertices(), false);
        std::vector<Vertex> result;
        for (uint32_t i = 0; i < k; i++) {
            Vertex vertex = randomVertex();
            while (picked[vertex]) {
                vertex = randomVertex();
            }
            picked[vertex] = true;
            result.emplace_back(vertex);
        }
        return result;
    }

    inline std::vector<Vertex> weightCenters(const uint32_t k) noexcept {
        std::vector<Vertex> vertices;
        std::vector<double> weights;
        for (const Vertex vertex : graph.vertices()) {
            vertices.emplace_back(vertex);
            weights.emplace_back(vertexWeight[vertex]);
        }
        Dijkstra<CSA::TransferGraph, false> dijkstra(graph);
        std::vector<Vertex> result;
        for (uint32_t i = 0; i < k; i++) {
            std::sort(vertices.begin(), vertices.end(), [&](const Vertex a, const Vertex b){return weights[a] > weights[b];});
            result.emplace_back(vertices.front());
            int count = 0;
            dijkstra.run(vertices.front(), noVertex, [&](const Vertex v){
                weights[v] = 0;
                count++;
            }, [&](){return count > 30;});
        }
        return result;
    }

    inline std::vector<Vertex> distCenters(const uint32_t k) noexcept {
        std::vector<Vertex> result;
        result.emplace_back(randomVertex());
        Dijkstra<CSA::TransferGraph, false> dijkstra(graph);
        for (uint32_t i = 1; i < k; i++) {
            Vertex nextVertex = randomVertex();
            dijkstra.run(result, noVertex, [&](const Vertex v){
                nextVertex = v;
            });
            result.emplace_back(nextVertex);
        }
        return result;
    }

    inline void assign(const std::vector<Vertex> centers) noexcept {
        for (const Vertex vertex : graph.vertices()) {
            columnOfVertex[vertex] = 0;
        }
        for (size_t i = 0; i < centers.size(); i++) {
            warning(centers[i], "   ", i);
            columnOfVertex[centers[i]] = i;
        }
    }

    inline void assignColumns(const std::vector<Vertex> centers) noexcept {
        // warning("ZZZ");
        timeStamp++;
        columnQueue.clear();
        for (uint32_t column = 0; column < centers.size(); column++) {
            labelOfColumn[column].setSource(centers[column], timeStamp);
            columnQueue.update(&(labelOfColumn[column]));
        }
        for (const Vertex vertex : graph.vertices()) {
            columnOfVertex[vertex] = undefined;
        }
        size_t assignedVertices = 0;
        while (assignedVertices < columnOfVertex.size()) {
            if (columnQueue.empty()) break;
            uint32_t column = columnQueue.front() - &(labelOfColumn[0]);
            if (labelOfColumn[column].vertexQueue.empty()) {
                columnQueue.pop();
            } else {
                VertexLabel* currentLabel = labelOfColumn[column].vertexQueue.extractFront();
                const Vertex currentVertex = Vertex(currentLabel - &(labelOfColumn[column].labelOfVertex[0]));
                if (columnOfVertex[currentVertex] != undefined) continue;
                labelOfColumn[column].weight += vertexWeight[currentVertex];
                columnQueue.update(&(labelOfColumn[column]));
                columnOfVertex[currentVertex] = column;
                assignedVertices++;
                for (const Edge edge : graph.edgesFrom(currentVertex)) {
                    const Vertex nextVertex = graph.get(ToVertex, edge);
                    if (columnOfVertex[nextVertex] != undefined) continue;
                    labelOfColumn[column].checkTimeStamp(nextVertex, timeStamp);
                    VertexLabel* nextLabel = &(labelOfColumn[column].labelOfVertex[nextVertex]);
                    const int distance = currentLabel->distance + graph.get(TravelTime, edge);
                    if (nextLabel->distance > distance) {
                        nextLabel->distance = distance;
                        labelOfColumn[column].vertexQueue.update(nextLabel);
                    }
                }
            }
        }
        for (const Vertex vertex : graph.vertices()) {
            if (columnOfVertex[vertex] == undefined) columnOfVertex[vertex] = 0;
        }
        // warning("ZZZ");
    }

    inline void assignColumnsByDist(const std::vector<Vertex> centers) noexcept {
        // warning("YYY");
        timeStamp++;
        columnQueue.clear();
        for (uint32_t column = 0; column < centers.size(); column++) {
            labelOfColumn[column].setSource(centers[column], timeStamp);
            columnQueue.update(&(labelOfColumn[column]));
        }
        for (const Vertex vertex : graph.vertices()) {
            columnOfVertex[vertex] = undefined;
        }
        size_t assignedVertices = 0;
        while (assignedVertices < columnOfVertex.size()) {
            if (columnQueue.empty()) break;
            uint32_t column = columnQueue.front() - &(labelOfColumn[0]);
            if (labelOfColumn[column].vertexQueue.empty()) {
                columnQueue.pop();
            } else {
                columnQueue.pop();
                VertexLabel* currentLabel = labelOfColumn[column].vertexQueue.extractFront();
                const Vertex currentVertex = Vertex(currentLabel - &(labelOfColumn[column].labelOfVertex[0]));
                if (!labelOfColumn[column].vertexQueue.empty()) {
                    labelOfColumn[column].weight = labelOfColumn[column].vertexQueue.front()->distance;
                    columnQueue.update(&(labelOfColumn[column]));
                }
                if (columnOfVertex[currentVertex] != undefined) continue;
                columnOfVertex[currentVertex] = column;
                assignedVertices++;
                for (const Edge edge : graph.edgesFrom(currentVertex)) {
                    const Vertex nextVertex = graph.get(ToVertex, edge);
                    if (columnOfVertex[nextVertex] != undefined) continue;
                    labelOfColumn[column].checkTimeStamp(nextVertex, timeStamp);
                    VertexLabel* nextLabel = &(labelOfColumn[column].labelOfVertex[nextVertex]);
                    const int distance = currentLabel->distance + graph.get(TravelTime, edge);
                    if (nextLabel->distance > distance) {
                        nextLabel->distance = distance;
                        labelOfColumn[column].vertexQueue.update(nextLabel);
                    }
                }
                if (!labelOfColumn[column].vertexQueue.empty()) {
                    labelOfColumn[column].weight = labelOfColumn[column].vertexQueue.front()->distance;
                    columnQueue.update(&(labelOfColumn[column]));
                }
            }
        }
        for (const Vertex vertex : graph.vertices()) {
            if (columnOfVertex[vertex] == undefined) columnOfVertex[vertex] = 0;
        }
        // warning("YYY");
    }

    inline std::vector<Vertex> geoCenters() noexcept {
        // warning("XXX");
        timeStamp++;
        for (const Vertex vertex : graph.vertices()) {
            labelOfColumn[columnOfVertex[vertex]].checkTimeStamp(vertex, timeStamp);
        }
        // warning("XXX1");
        for (uint32_t column = 0; column < labelOfColumn.size(); column++) {
            labelOfColumn[column].vertexQueue.clear();
        }
        // warning("XXX2");
        for (const Vertex from : graph.vertices()) {
            for (const Edge edge : graph.edgesFrom(from)) {
                const Vertex to = graph.get(ToVertex, edge);
                if (columnOfVertex[from] != columnOfVertex[to]) {
                    labelOfColumn[columnOfVertex[from]].labelOfVertex[from].distance = std::min(labelOfColumn[columnOfVertex[from]].labelOfVertex[from].distance, graph.get(TravelTime, edge));
                    labelOfColumn[columnOfVertex[to]].labelOfVertex[to].distance = std::min(labelOfColumn[columnOfVertex[to]].labelOfVertex[to].distance, graph.get(TravelTime, edge));
                    labelOfColumn[columnOfVertex[from]].vertexQueue.update(&(labelOfColumn[columnOfVertex[from]].labelOfVertex[from]));
                    labelOfColumn[columnOfVertex[to]].vertexQueue.update(&(labelOfColumn[columnOfVertex[to]].labelOfVertex[to]));
                }
            }
        }
        // warning("XXX3");
        std::vector<Vertex> result(labelOfColumn.size());
        for (uint32_t column = 0; column < result.size(); column++) {
            while (!labelOfColumn[column].vertexQueue.empty()) {
                VertexLabel* currentLabel = labelOfColumn[column].vertexQueue.extractFront();
                const Vertex currentVertex = Vertex(currentLabel - &(labelOfColumn[column].labelOfVertex[0]));
                result[column] = currentVertex;
                for (const Edge edge : graph.edgesFrom(currentVertex)) {
                    const Vertex nextVertex = graph.get(ToVertex, edge);
                    if (columnOfVertex[nextVertex] != column) continue;
                    VertexLabel* nextLabel = &(labelOfColumn[column].labelOfVertex[nextVertex]);
                    const int distance = currentLabel->distance + graph.get(TravelTime, edge);
                    if (nextLabel->distance > distance) {
                        nextLabel->distance = distance;
                        labelOfColumn[column].vertexQueue.update(nextLabel);
                    }
                }
            }
        }
        // warning("XXX");
        return result;
    }

    inline Vertex randomVertex() noexcept {
        return Vertex(randomDistribution(randomGenerator));
    }

public:
    const CSA::TransferGraph graph;
    std::vector<double> vertexWeight;

    std::vector<uint32_t> columnOfVertex;
    std::vector<ColumnLabel> labelOfColumn;
    ExternalKHeap<2, ColumnLabel> columnQueue;
    int timeStamp;

    std::mt19937 randomGenerator;
    std::uniform_int_distribution<uint32_t> randomDistribution;

};

}
