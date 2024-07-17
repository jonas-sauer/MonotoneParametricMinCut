#pragma once

#include <deque>
#include <vector>

#include "Point.h"

namespace Geometry {

inline constexpr double EPS = 1e-6;

struct Halfplane {

public:
    Halfplane(const Point& p, const Point& v) :
        p(p),
        v(v) {
    }

    Halfplane(const double slope, const double offset, const bool direction) :
        p(Construct::XY, 0, offset),
        v(Construct::XY, 1, slope) {
        if (!direction) v *= -1;
    }

    inline double getAngle() const noexcept {
        return atan2(v.y, v.x);
    }

    inline bool operator<(const Halfplane& other) const noexcept {
        return getAngle() < other.getAngle();
    }

    inline bool isPointRightOf(const Point& o) const noexcept {
        return crossProduct(v, o - p) <= 0;
    }

    friend Point getIntersection(const Halfplane& a, const Halfplane& b) {
        const double x = crossProduct(b.p - a.p, b.v) / crossProduct(a.v, b.v);
        return a.p + x*a.v;
    }

public:
    Point p;
    Point v;
};

inline std::vector<Halfplane> createBoundingBox(const double arrivalFactor = INFTY, const double tripsPerHour = INFTY) noexcept {
    std::vector<Halfplane> result;
    const Point p1(Construct::XY, 0, 0);
    const Point p2(Construct::XY, 1, 0);
    const Point p3(Construct::XY, 0, 1);
    result.emplace_back(p1, p2 - p1);
    result.emplace_back(p2, p3 - p2);
    result.emplace_back(p3, p1 - p3);
    if (arrivalFactor < INFTY) result.emplace_back(arrivalFactor, 0, false);
    if (tripsPerHour < INFTY) {
        const double tripFactor = tripsPerHour/(tripsPerHour + 1);
        result.emplace_back(-tripFactor, tripFactor, false);
    }
    return result;
}

//Based on https://cp-algorithms.com/geometry/halfplane-intersection.html
inline std::vector<Point> intersectHalfplanes(std::vector<Halfplane> halfplanes, const double arrivalFactor = INFTY, const double tripsPerHour = INFTY) noexcept {
    halfplanes += createBoundingBox(arrivalFactor, tripsPerHour);
    std::sort(halfplanes.begin(), halfplanes.end());
    std::deque<Halfplane> deque;

    size_t qSize = 0;
    for (size_t i = 0; i < halfplanes.size(); i++) {
        //Remove redundant halfplanes from back of queue
        while (qSize > 1 && halfplanes[i].isPointRightOf(getIntersection(deque[qSize-1], deque[qSize-2]))) {
            deque.pop_back();
            qSize--;
        }

        //Remove redundant halfplanes from front of queue
        while (qSize > 1 && halfplanes[i].isPointRightOf(getIntersection(deque[0], deque[1]))) {
            deque.pop_front();
            qSize--;
        }

        //Check for parallel halfplanes
        if (qSize > 0 && fabs(crossProduct(halfplanes[i].v, deque[qSize-1].v)) < EPS) {
            // Opposite direction -> Intersection is empty
            if (dotProduct(halfplanes[i].v, deque[qSize-1].v) < 0.0)
                return std::vector<Point>();

            // Same direction -> redundant
            if (halfplanes[i].isPointRightOf(deque[qSize-1].p)) {
                deque.pop_back();
                qSize--;
            }
            else continue;
        }

        deque.push_back(halfplanes[i]);
        qSize++;
    }

    //Finalize: Remove redundant halfplanes
    while (qSize > 2 && deque[0].isPointRightOf(getIntersection(deque[qSize-1], deque[qSize-2]))) {
        deque.pop_back();
        qSize--;
    }

    while (qSize > 2 && deque[qSize-1].isPointRightOf(getIntersection(deque[0], deque[1]))) {
        deque.pop_front();
        qSize--;
    }

    if (qSize < 3) return std::vector<Point>();

    //Construct polygon
    std::vector<Point> ret(qSize);
    for (size_t i = 0; i+1 < qSize; i++) {
        ret[i] = getIntersection(deque[i], deque[i+1]);
    }
    ret.back() = getIntersection(deque[qSize-1], deque[0]);

    return ret;
}

}
