#pragma once

#include <cstdlib>

#include "../../Helpers/FloatingPointMath.h"
#include "../../Helpers/Types.h"
#include "../../Helpers/IO/Serialization.h"

/**
 * Abstract class for flow functions in a single parameter
 */
class FlowFunction {
    /**
     * Returns the zero crossing of f with minimum x coordinate for which x >= minVal, or INFTY if none exists
     */
    [[nodiscard]] inline virtual double getNextZeroCrossing(const double /*minVal*/) const {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    [[nodiscard]] inline virtual double eval(const double) const {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    friend inline double findIntersectionPoint(const FlowFunction &, const FlowFunction &) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    friend inline FlowFunction &operator+(const FlowFunction &, const FlowFunction &) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    friend inline FlowFunction &operator-(const FlowFunction &, const FlowFunction &) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    inline virtual FlowFunction operator+=(const FlowFunction &) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    inline virtual FlowFunction operator-=(const FlowFunction &) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    /**
     * Equality that takes floating-point inaccuracy into account
     */
    inline virtual bool operator==(const FlowFunction &) const {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    inline virtual void serialize(IO::Serialization&) const {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    };

    inline virtual void deserialize(IO::Deserialization&) {
        throw std::runtime_error("Cannot use FlowFunction directly, use a derived class instead!");
    }
};

/**
* Linear flow function of the form f(x) = ax + b
*/
class LinearFlowFunction : public FlowFunction {
public:
    using FlowType = double;
    inline static constexpr double CONST_INF = INFTY;

    LinearFlowFunction(double a, double b) : a(std::min(a, CONST_INF)), b(std::min(b, CONST_INF)) {};

    explicit LinearFlowFunction(double b) : a(0), b(std::min(b, CONST_INF)) {};

    explicit LinearFlowFunction() : a(0), b(0) {};

    [[nodiscard]] inline double getNextZeroCrossing(const double minVal) const noexcept override {
        if (a == 0) {
            if (b <= FLOAT_EPS) {
                return minVal;
            } else {
                return CONST_INF;
            }
        }
        if (b >= CONST_INF) return CONST_INF;
        const double crossing = -(b / a);
        if (minVal > crossing) return CONST_INF;
        else return crossing;
    }

    [[nodiscard]] inline double eval(const double x) const noexcept override {
        if (x == CONST_INF) return a;
        if (b == CONST_INF) return CONST_INF;
        return x * a + b;
    }

    friend inline double findIntersectionPoint(const LinearFlowFunction& lhs, const LinearFlowFunction& rhs) noexcept {
        if (lhs.a == rhs.a) return CONST_INF;
        return (rhs.b - lhs.b)/(lhs.a - rhs.a);
    }

    friend inline LinearFlowFunction operator+(const LinearFlowFunction& lhs, const LinearFlowFunction &rhs) noexcept {
        LinearFlowFunction result = lhs;
        result.a += rhs.a;
        if (lhs.b == CONST_INF || rhs.b == CONST_INF) {
            result.b = CONST_INF;
        } else {
            result.b += rhs.b;
        }
        return result;
    }

    friend inline LinearFlowFunction operator-(const LinearFlowFunction& lhs, const LinearFlowFunction &rhs) noexcept {
        LinearFlowFunction result = lhs;
        result.a -= rhs.a;
        if (lhs.b != CONST_INF) {
            result.b -= rhs.b;
        }
        return result;
    }

    inline LinearFlowFunction operator-() const noexcept {
        return {-this->a, -this->b};
    }

    inline LinearFlowFunction& operator=(const int rhs) noexcept {
        a = 0;
        b = rhs;
        return *this;
    }

    inline LinearFlowFunction operator+=(const LinearFlowFunction &rhs) noexcept {
        a += rhs.a;
        if (b == CONST_INF || rhs.b == CONST_INF) {
            b = CONST_INF;
        } else {
            b += rhs.b;
        }
        return *this;
    }

    inline LinearFlowFunction operator-=(const LinearFlowFunction &rhs) noexcept {
        a -= rhs.a;
        if (b != CONST_INF) {
            b -= rhs.b;
        }
        return *this;
    }

    inline bool operator==(const LinearFlowFunction &rhs) const noexcept {
        return areNumbersEqualAbsolute(a, rhs.a) && areNumbersEqualAbsolute(b, rhs.b);
    }

    inline void serialize(IO::Serialization& serialize) const noexcept override {
        serialize(a, b);
    }

    inline void deserialize(IO::Deserialization& deserialize) noexcept override {
        deserialize(a, b);
    }

    inline friend std::ostream& operator<<(std::ostream& out, const LinearFlowFunction& f) noexcept {
        out << f.a << "*x+" << f.b;
        return out;
    }

private:
    double a, b;
};