#pragma once

#include <cmath>
#include <vector>
#include <limits>

#include "Point.h"

namespace Geometry::Function {

template<typename FUNCTION>
class Inverse {

public:
    using Function = FUNCTION;
    using Type = Inverse<Function>;

public:
    Inverse(const Function* const function, double* const variable) :
        function(function),
        variable(variable) {
    }

    inline void operator=(const double y) const noexcept {
        (*variable) = function->evaluateInverse(y);
    }

private:
    const Function* function;
    double* variable;

};

class Linear {

public:
    Linear(const double a, const double b) : a(a), b(b) {}
    Linear(const double a, const Point& p) : a(a),
        b(p.y - (a * p.x)) {
        // std::cout << (*this) << std::endl;
    }
    Linear(const Point& p, const Point& q) :
        Linear((p.y - q.y) / (p.x - q.x), p) {
    }

    inline double operator()(const double x) const noexcept {
        return evaluate(x);
    }

    inline Inverse<Linear> operator()(double* const x) const noexcept {
        return Inverse<Linear>(this, x);
    }

    inline double evaluate(const double x) const noexcept {
        return (a * x) + b;
    }

    inline double evaluateInverse(const double y) const noexcept {
        Assert((a != 0), "No inverse value for " << y << " does exist!");
        return (y - b) / a;
    }

    friend std::ostream& operator<<(std::ostream& out, const Linear& f) {
        return out << "f(x) = (" << f.a << " * x) + " << f.b;
    }

private:
    double a;
    double b;

};

class Exponential {

public:
    Exponential(const double a, const double b, const double c, const double d, const double e) : a(a), b(b), c(c), d(d), e(e) {
        Assert(c != 0.0, "Invalid function");
    }
    Exponential(const Point& p, const Point& q, const Point& r) :
        a(Linear((p.y - r.y) / (p.x - r.x), q)((p.x + r.x) / 2.0) - p.y),
        b((r.y - p.y - a) / a),
        c((r.x - p.x) / 2.0),
        d(p.x),
        e(p.y) {
        Assert(c != 0.0, "Invalid function");
        Assert((p.x < q.x) && (q.x < r.x), "Points are not sorted by their x coordinate!");
        Assert(((p.y < q.y) ? (q.y < r.y) : (q.y > r.y)), "Points are not sorted by their y coordinate!");
        // std::cout << p << ", " << q << ", " << r << "\na: " << a << ", b: " << b << ", c: " << c << "\n" << (*this) << std::endl;
    }

    inline double operator()(const double x) const noexcept {
        return evaluate(x);
    }

    inline Inverse<Exponential> operator()(double* const x) const noexcept {
        return Inverse<Exponential>(this, x);
    }

    inline double evaluate(const double x) const noexcept {
        if (b == 1.0) {
            return ((a / c) * (x - d)) + e;
        } else {
            return ((a * (std::pow(b, (x - d) / c) - 1.0)) / (b - 1.0)) + e;
        }
    }

    inline double evaluateInverse(const double y) const noexcept {
        if (b == 1.0) {
            return ((y - e) / (a / c)) + d;
        } else {
            Assert(((b < 1.0) ? (y < limit()) : (y > limit())), "No inverse value for " << y << " does exist!");
            return ((c * std::log1p(((b - 1.0) * (y - e)) / a)) / std::log(b)) + d;
        }
    }

    inline bool hasInverse(const double y) const noexcept {
        if (b == 1.0) {
            return true;
        } else {
            return (b < 1.0) ? (y < limit()) : (y > limit());
        }
    }

    inline double limit() const noexcept {
        return (a / (1 - b)) + e;
    }

    inline void straighten(const double factor) noexcept {
        const double oldOffset = evaluate(d + (2.0 * a));
        b = 1.0 - ((1.0 - b) * (1.0 - factor));
        const double newOffset = evaluate(d + (2.0 * a));
        e = e - newOffset + oldOffset;
    }

    friend std::ostream& operator<<(std::ostream& out, const Exponential& f) {
        return out << "f(x) = ((" << f.a << " * ((" << f.b << "^((x - " << f.d << ") / " << f.c << ")) - 1)) / (" << f.b << " - 1)) + " << f.e;
    }

private:
    double a;
    double b;
    double c;
    double d;
    double e;

};

class Root {

public:
    Root(const double a, const double b, const double c) : a(a), b(b), c(c) {}
    Root(const Point& p, const Point& q, const Point& r) : a(0.0), b(0.0), c(0.0) {
        const double ypq = p.y - q.y;
        const double yqr = q.y - r.y;
        const double ypr = p.y - r.y;
        const double y = ypq * yqr * ypr;
        Assert(y != 0, "Invalid function");
        a = ((r.x * ypq) + (p.x * yqr) - (q.x * ypr)) / y;
        b = ((r.x * q.y * q.y) - (r.x * p.y * p.y) + (q.x * p.y * p.y) - (q.x * r.y * r.y) + (p.x * r.y * r.y) - (p.x * q.y * q.y)) / y;
        c = ((r.x * p.y * q.y * ypq) + (p.x * q.y * r.y * yqr) - (q.x * p.y * r.y * ypr)) / y;
        // std::cout << "A=" << p << ", B=" << q << ", C=" << r << "\na=" << a << ", b=" << b << ", c=" << c << "\n" << (*this) << std::endl;
        // std::cout << "A=(" << p.y << ", " << p.x << "), B=(" << q.y << ", " << q.x << "), C=(" << r.y << ", " << r.x << ")\na=" << a << ", b=" << b << ", c=" << c << "\n(" << a << " * y * y) + (" << b << " * y) + " << c << std::endl;
    }

    inline double operator()(const double x) const noexcept {
        return evaluate(x);
    }

    inline Inverse<Root> operator()(double* const x) const noexcept {
        return Inverse<Root>(this, x);
    }

    inline double evaluate(const double x) const noexcept {
        return (sqrt((4.0 * a * (x - c)) + (b * b)) - b) / (2.0 * a);
    }

    inline double evaluateInverse(const double y) const noexcept {
        return (a * y * y) + (b * y) + c;
    }

    friend std::ostream& operator<<(std::ostream& out, const Root& f) {
        return out << "(sqrt((4.0 * " << f.a << " * (x - " << f.c << ")) + (" << f.b << " * " << f.b << ")) - " << f.b << ") / (2.0 * " << f.a << ")";
    }

private:
    double a;
    double b;
    double c;

};

class Piecewise {

public:
    Piecewise() {}
    Piecewise(const std::vector<Point>& data) : data(data) {}
    Piecewise(std::vector<Point>&& data) : data(std::move(data)) {}

    inline operator std::vector<Point>() const noexcept {
        return data;
    }

    inline double operator()(const double x) const noexcept {
        return evaluate(x);
    }

    inline Inverse<Piecewise> operator()(double* const x) const noexcept {
        return Inverse<Piecewise>(this, x);
    }

    inline double evaluate(const double x) const noexcept {
        if (data.empty() || x < data.front().x || x > data.back().x) return std::numeric_limits<double>::quiet_NaN();
        size_t i = 0;
        while (data[i].x < x) {
            i++;
        }
        if (data[i].x == x) {
            double result = data[i].y;
            while (i < data.size() && data[i].x == x) {
                if (result > data[i].y) {
                    result = data[i].y;
                }
                i++;
            }
            return result;
        }
        return Linear(data[i - 1], data[i])(x);
    }

    inline double evaluateInverse(const double y) const noexcept {
        if (y == data.front().y) {
            return data.front().x;
        }
        size_t i = 1;
        while (i < data.size()) {
            if (y == data[i].y) {
                return data[i].x;
            }
            if (((data[i].y > y) && (y > data[i - 1].y)) || ((data[i].y < y) && (y < data[i - 1].y))) {
                return Linear(data[i - 1], data[i]).evaluateInverse(y);
            }
            i++;
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    inline Linear extrapolateLinear(const size_t minIndex, const size_t maxIndex) const noexcept {
        Assert(data.size() >= 2, "To few data points for extrapolation!");
        const size_t i = std::min(data.size() - 2, minIndex);
        const size_t j = std::max(i + 1, std::min(data.size() - 1, maxIndex));
        return Linear(data[i], data[j]);
    }

    inline Linear extrapolateLinear(const size_t minIndex = 0) const noexcept {
        return extrapolateLinear(minIndex, data.size() - 1);
    }

    inline Exponential extrapolateExponential(const size_t minIndex, const size_t maxIndex) const noexcept {
        Assert(data.size() >= 3, "To few data points for extrapolation!");
        const size_t i = std::min(data.size() - 3, minIndex);
        const size_t j = std::max(i + 1, std::min(data.size() - 2, (minIndex + maxIndex) / 2));
        const size_t k = std::max(j + 1, std::min(data.size() - 1, maxIndex));
        return Exponential(data[i], data[j], data[k]);
    }

    inline Exponential extrapolateExponential(const size_t minIndex = 0) const noexcept {
        return extrapolateExponential(minIndex, data.size() - 1);
    }

    inline Root extrapolateRoot(const size_t minIndex, const size_t maxIndex) const noexcept {
        Assert(data.size() >= 3, "To few data points for extrapolation!");
        const size_t i = std::min(data.size() - 3, minIndex);
        const size_t j = std::max(i + 1, std::min(data.size() - 2, (minIndex + maxIndex) / 2));
        const size_t k = std::max(j + 1, std::min(data.size() - 1, maxIndex));
        return Root(data[i], data[j], data[k]);
    }

    inline Root extrapolateRoot(const size_t minIndex = 0) const noexcept {
        return extrapolateRoot(minIndex, data.size() - 1);
    }

    inline void append(const Point& p) noexcept {
        Assert(((data.empty()) || (data.back().x <= p.x)), "The function is already defined for " << p.x << "!");
        data.emplace_back(p);
    }

    inline void append(const double x, const double y) noexcept {
        Assert(((data.empty()) || (data.back().x <= x)), "The function is already defined for " << x << "!");
        data.emplace_back(Point(Construct::XY, x, y));
    }

    inline size_t size() noexcept {
        return data.size();
    }

    inline void clear() noexcept {
        data.clear();
    }

private:
    std::vector<Point> data;

};

}
