#pragma once

#include "../UnitTests.h"

#include "../../DataStructures/Geometry/CoordinateTree.h"

namespace UnitTests {

class CoordinateTree {

public:
    inline void check() {
        std::vector<Geometry::Point> coordinates;
        coordinates.emplace_back(Geometry::Point(Construct::XY, 0, 0));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 0, 1));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 0));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 1));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 2));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 5));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 8));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 1, 10));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 3, 0));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 3, 2));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 3, 10));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 4, 2));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 4, 5));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 6, 3));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 6, 7));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 8, 5));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 10, 0));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 10, 3));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 10, 7));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 13, 0));
        coordinates.emplace_back(Geometry::Point(Construct::XY, 13, 3));
        ::CoordinateTree<Geometry::EuclideanMetric> ct(Geometry::EuclideanMetric(), coordinates, 3);
        for (size_t i = 0; i < coordinates.size(); i++) {
            UnitTests::check(ct.getNearestNeighbor(coordinates[i]) == static_cast<Vertex>(i), "Point i = ", coordinates[i], " has nearest neighbor ", ct.getNearestNeighbor(coordinates[i]), " = ", coordinates[ct.getNearestNeighbor(coordinates[i])], "!");
        }
        UnitTests::check(ct.getNearestNeighbor(Geometry::Point(Construct::XY, 20, 0)) == 19, "Point (20, 0) has nearest neighbor ", ct.getNearestNeighbor(Geometry::Point(Construct::XY, 20, 0)), " != 19!");
        UnitTests::check(ct.getNearestNeighbor(Geometry::Point(Construct::XY, -1, -1)) == 0, "Point (-1, -1) has nearest neighbor ", ct.getNearestNeighbor(Geometry::Point(Construct::XY, -1, -1)), " != 0!");
        UnitTests::check(ct.getNearestNeighbor(Geometry::Point(Construct::XY, 0.5, 0.5)) < 4, "Point (0.5, 0.5) has nearest neighbor ", ct.getNearestNeighbor(Geometry::Point(Construct::XY, 0.5, 0.5)), " >= 4!");
        UnitTests::check(ct.getNearestNeighbor(Geometry::Point(Construct::XY, 2, 7)) == 6, "Point (2, 7) has nearest neighbor ", ct.getNearestNeighbor(Geometry::Point(Construct::XY, 2, 7)), " != 6!");
        UnitTests::check(ct.getNearestNeighbor(Geometry::Point(Construct::XY, 6.00001, 5)) == 15, "Point (6.00001, 5) has nearest neighbor ", ct.getNearestNeighbor(Geometry::Point(Construct::XY, 6.00001, 5)), " != 15!");
        std::vector<Vertex> n1 = ct.getNeighbors(Geometry::Point(Construct::XY, 6, 5), 2);
        UnitTests::check(n1.size() == 4, "Point (6, 5) has ", n1.size(), " neighbors != 4");
        UnitTests::check(Vector::contains(n1, Vertex(12)), "Neighbors of point (6, 5) is missing point 12 = ", coordinates[12], "!");
        UnitTests::check(Vector::contains(n1, Vertex(13)), "Neighbors of point (6, 5) is missing point 13 = ", coordinates[13], "!");
        UnitTests::check(Vector::contains(n1, Vertex(14)), "Neighbors of point (6, 5) is missing point 14 = ", coordinates[14], "!");
        UnitTests::check(Vector::contains(n1, Vertex(15)), "Neighbors of point (6, 5) is missing point 15 = ", coordinates[15], "!");
        std::vector<Vertex> n2 = ct.getNeighbors(Geometry::Point(Construct::XY, 1, 1.5), 1.5);
        UnitTests::check(n2.size() == 4, "Point (1, 1.5) has ", n2.size(), " neighbors != 4");
        UnitTests::check(Vector::contains(n2, Vertex(1)), "Neighbors of point (1, 1.5) is missing point 1 = ", coordinates[1], "!");
        UnitTests::check(Vector::contains(n2, Vertex(2)), "Neighbors of point (1, 1.5) is missing point 2 = ", coordinates[2], "!");
        UnitTests::check(Vector::contains(n2, Vertex(3)), "Neighbors of point (1, 1.5) is missing point 3 = ", coordinates[3], "!");
        UnitTests::check(Vector::contains(n2, Vertex(4)), "Neighbors of point (1, 1.5) is missing point 4 = ", coordinates[4], "!");
    }

};

}
