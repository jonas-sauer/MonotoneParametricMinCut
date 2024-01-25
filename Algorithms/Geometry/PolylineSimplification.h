#pragma once

#include <vector>

#include "../../Helpers/Vector/Vector.h"

#include "../../DataStructures/Geometry/Point.h"
#include "../../DataStructures/Geometry/Rectangle.h"

namespace Geometry {

class ParableLine {

public:
    ParableLine(const double f = 5) : f(f) {}

    inline bool operator()(const Point& lineA, const Point& lineB, const Point& point) const noexcept {
        const Point m = (lineA + lineB) / 2;
        return (point.distanceToPoint(m) + (f * point.distanceToLine(lineA, lineB))) <= (lineA.distanceToPoint(m));
    }

private:
    double f;

};

class AbsoluteParableLine {

public:
    AbsoluteParableLine(const double abs, const double f = 5) : abs(abs), f(f) {}

    inline bool operator()(const Point& lineA, const Point& lineB, const Point& point) const noexcept {
        const Point m = (lineA + lineB) / 2;
        return (point.distanceToPoint(m) + (f * point.distanceToLine(lineA, lineB))) <= abs;
    }

private:
    double abs;
    double f;

};

class BoxLine {

public:
    BoxLine(const double f = 12) : f(f) {}

    inline bool operator()(const Point& lineA, const Point& lineB, const Point& point) const noexcept {
        return (f * point.distanceToLine(lineA, lineB)) <= lineA.distanceToPoint(lineB);
    }

private:
    double f;

};

class AbsoluteBoxLine {

public:
    AbsoluteBoxLine(const double abs, const double f = 12) : abs(abs), f(f) {}

    inline bool operator()(const Point& lineA, const Point& lineB, const Point& point) const noexcept {
        return (f * point.distanceToLine(lineA, lineB)) <= abs;
    }

private:
    double abs;
    double f;

};

class PolyLine {

public:
    PolyLine(const std::vector<Point>& points) :
        points(points),
        important(points.size(), false),
        boundingBox(Rectangle::BoundingBox(points)),
        cache(0) {
        important.front() = true;
        important.back() = true;
    }

    inline void markCorners(PolyLine& other) noexcept {
        if (!boundingBox.intersects(other.boundingBox)) return;
        for (size_t i = 0; i < points.size(); i++) {
            const size_t j = other.index(points[i]);
            if (j == other.points.size()) continue;
            if (important[i]) {
                other.important[j] = true;
            } else if (other.important[j]) {
                important[i] = true;
            } else if ((points[i - 1] != other.points[j - 1]) && (points[i - 1] != other.points[j + 1])) {
                other.important[j] = true;
                important[i] = true;
            } else if ((points[i + 1] != other.points[j - 1]) && (points[i + 1] != other.points[j + 1])) {
                other.important[j] = true;
                important[i] = true;
            }
        }
    }

    inline size_t index(const Point& p) const noexcept {
        if (!boundingBox.contains(p)) {
            return points.size();
        }
        if (cache > 0 && points[cache - 1] == p) {
            cache--;
            return cache;
        }
        if (cache < points.size() - 1 && points[cache + 1] == p) {
            cache++;
            return cache;
        }
        for (size_t i = 0; i < points.size(); i++) {
            if (points[i] == p) {
                cache = i;
                return cache;
            }
        }
        return points.size();
    }

    inline size_t countImportant() const noexcept {
        size_t result = 0;
        for (size_t i = 0; i < important.size(); i++) {
            if (important[i]) result++;
        }
        return result;
    }

public:
    std::vector<Point> points;
    std::vector<bool> important;
    Rectangle boundingBox;
    mutable size_t cache;

};

template<typename IS_LINE>
inline std::vector<Point> simplify(const std::vector<Point>& points, const IS_LINE& isLine) noexcept {
    if (points.size() <= 2) return points;
    std::vector<size_t> indices{0, 1};
    for (size_t i = 2; i < points.size(); i++) {
        bool replace = true;
        const size_t k = indices[indices.size() - 2];
        for (size_t j = k + 1; j < i; j++) {
            if (!isLine(points[k], points[i], points[j])) {
                replace = false;
                break;
            }
        }
        if (replace) {
            indices.back() = i;
        } else {
            indices.emplace_back(i);
        }
    }
    std::vector<Point> result;
    result.reserve(indices.size());
    for (const size_t i : indices) {
        result.emplace_back(points[i]);
    }
    return result;
}

template<typename IS_LINE>
inline std::vector<Point> simplify(const PolyLine& polyline, const IS_LINE& isLine) noexcept {
    if (polyline.points.size() <= 2) return polyline.points;
    std::vector<Point> result;
    std::vector<Point> temp;
    for (size_t i = 0; i < polyline.points.size(); i++) {
        if (polyline.important[i] && !temp.empty()) {
            temp.emplace_back(polyline.points[i]);
            if (temp.front().x < temp.back().x || ((temp.front().x == temp.back().x) && (temp.front().y < temp.back().y))) {
                Vector::reverse(temp);
                temp = simplify(temp, isLine);
                Vector::reverse(temp);
            } else {
                temp = simplify(temp, isLine);
            }
            temp.pop_back();
            result += temp;
            temp.clear();
        }
        temp.emplace_back(polyline.points[i]);
    }
    result.emplace_back(polyline.points.back());
    return result;
}

}
