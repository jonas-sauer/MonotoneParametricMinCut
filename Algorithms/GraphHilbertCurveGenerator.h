#pragma once

#include <iostream>
#include <vector>

#include "../DataStructures/Geometry/Point.h"
#include "../DataStructures/Geometry/Rectangle.h"
#include "../DataStructures/Graph/Graph.h"
#include "../Visualization/MapVisualization.h"

template<typename GRAPH>
class GraphHilbertCurveGenerator {
private:
    struct VertexPoint {
        VertexPoint(const Vertex vertex, const Geometry::Point& coordinates) :
            vertices({vertex}),
            coordinates(coordinates) {
        }

        inline bool operator<(const VertexPoint& other) const noexcept {
            return coordinates.x < other.coordinates.x || (coordinates.x == other.coordinates.x && coordinates.y < other.coordinates.y);
        }

        std::vector<Vertex> vertices;
        Geometry::Point coordinates;
    };

    struct Cell {
        Cell(const Geometry::Rectangle& box, const size_t x, const size_t y) :
            boundingBox(box),
            x(x),
            y(y),
            curvePosition(0) {
        }

        Cell(const GRAPH& graph) :
            boundingBox(Graph::boundingBox(graph)),
            x(0),
            y(0),
            curvePosition(0) {
            std::vector<VertexPoint> vertices;
            for (const Vertex v : graph.vertices()) {
                vertices.emplace_back(v, graph.get(Coordinates, v));
            }
            std::sort(vertices.begin(), vertices.end());
            for (const VertexPoint& v : vertices) {
                addPoint(v.vertices[0], v.coordinates);
            }
        }

        inline void addPoint(const Vertex v, const Geometry::Point& coordinates) noexcept {
            if (points.empty() || points.back().coordinates != coordinates) {
                points.emplace_back(v, coordinates);
            } else {
                points.back().vertices.emplace_back(v);
            }
        }

        inline std::vector<Cell> divide() const noexcept {
            if (points.size() == 1) {
                Cell result(boundingBox, 2*x, 2*y);
                result.points = points;
                return std::vector<Cell>{ result };
            }

            const Geometry::Point topLeft = boundingBox.topLeft;
            const Geometry::Point bottomRight = boundingBox.bottomRight;

            const double xmin = topLeft.x;
            const double xmax = bottomRight.x;
            const double xmid = (xmax + xmin)/2;
            const double ymin = bottomRight.y;
            const double ymax = topLeft.y;
            const double ymid = (ymax + ymin)/2;

            const Geometry::Point topMid(Construct::XY, xmid, ymax);
            const Geometry::Point midMid(Construct::XY, xmid, ymid);
            const Geometry::Point midRight(Construct::XY, xmax, ymid);
            const Geometry::Point midLeft(Construct::XY, xmin, ymid);
            const Geometry::Point bottomMid(Construct::XY, xmid, ymin);

            std::vector<Cell> cells;
            cells.emplace_back(Geometry::Rectangle::BoundingBox(midLeft, bottomMid), 2*x, 2*y);
            cells.emplace_back(Geometry::Rectangle::BoundingBox(topLeft, midMid), 2*x, 2*y+1);
            cells.emplace_back(Geometry::Rectangle::BoundingBox(midMid, bottomRight), 2*x+1, 2*y);
            cells.emplace_back(Geometry::Rectangle::BoundingBox(topMid, midRight), 2*x+1, 2*y+1);

            for (const VertexPoint& v : points) {
                for (Cell& cell : cells) {
                    if (cell.boundingBox.contains(v.coordinates)) {
                        cell.points.emplace_back(v);
                        break;
                    }
                }
            }

            Vector::removeIf(cells, [&](const Cell& cell) {
               return cell.points.empty();
            });

            return cells;
        }

        inline void computeCurvePosition(const size_t recursionDepth) noexcept {
            curvePosition = 0;
            size_t currentX = x;
            size_t currentY = y;
            for (size_t s = recursionDepth/2; s > 0; s /= 2) {
                const size_t rx = (currentX & s) > 0;
                const size_t ry = (currentY & s) > 0;
                curvePosition += s * s * ((3 * rx) ^ ry);
                if (ry == 0) {
                    if (rx == 1) {
                        currentX = recursionDepth-currentX-1;
                        currentY = recursionDepth-currentY-1;
                    }
                    std::swap(currentX, currentY);
                }
            }
        }

        inline bool operator<(const Cell& other) const noexcept {
            return curvePosition < other.curvePosition;
        }

        Geometry::Rectangle boundingBox;
        size_t x;
        size_t y;
        std::vector<VertexPoint> points;
        size_t curvePosition;
    };

public:
    GraphHilbertCurveGenerator(const GRAPH& graph) :
        graph(graph),
        cells{ Cell(graph) } {
    }

    inline void run(const size_t maxRecursionDepth, const bool verbose = false) noexcept {
        if (verbose) std::cout << "Computing Hilbert curve with max. recursion depth " << maxRecursionDepth << "." << std::endl;

        size_t recursionDepth = 1;
        for (; recursionDepth <= maxRecursionDepth; recursionDepth++) {
            std::vector<Cell> newCells;
            for (const Cell& cell : cells) {
                newCells += cell.divide();
            }
            if (cells.size() == newCells.size()) break;
            cells = newCells;
            if (verbose) std::cout << "Cells after step " << recursionDepth << ": " << cells.size() << "/" << graph.numVertices() << std::endl;
        }

        for (Cell& cell : cells) {
            cell.computeCurvePosition(1 << recursionDepth);
        }
        std::sort(cells.begin(), cells.end());
    }

    inline Order getVertexOrder() const noexcept {
        Order order;
        for (const Cell& cell : cells) {
            for (const VertexPoint& point : cell.points) {
                for (const Vertex v : point.vertices) {
                    order.emplace_back(v);
                }
            }
        }
        return order;
    }

    template<typename FORMAT>
    inline void visualize(const std::string& visualizationFile) const noexcept {
        MapVisualization<FORMAT> visualization(visualizationFile, Graph::boundingBox(graph));
        for (size_t i = 1; i < cells.size(); i++) {
            visualization.drawLine(cells[i-1].boundingBox.center(), cells[i].boundingBox.center(), Color::getGradientColor(Color::KITgreen, Color::KITorange, Color::KITred, (double)i/cells.size()));
        }
    }

private:
    const GRAPH& graph;
    std::vector<Cell> cells;
};
