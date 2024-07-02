#pragma once

#include <concepts>
#include <numeric>
#include <cstdlib>

namespace pmf {

    constexpr double epsilon = 1E-9;

    /**
     * Function to check if two double values are almost equal
     * @param a
     * @param b
     * @return true if |a - b| < epsilon
     */
    bool doubleEqualAbs(double a, double b) {
        return std::abs(a - b) < epsilon;
    }

    /**
     * Function to check if two double values are almost equal taking into account growing imprecision with higher numbers
     * @param a
     * @param b
     * @return true if |a - b| < |a| * epsilon
     */
    bool doubleEqualRelative(double a, double b) {
        return std::abs(a - b) < std::abs(a * epsilon);
    }

    /**
     * Function to check a < b while accounting for double inequality
     * @param a
     * @param b
     * @return true iff a + epsilon < b
     */
    bool doubleLessThanAbs(double a, double b) {
        return a + epsilon < b;
    }

    /**
     * Function to check a > 0 while accounting for double imprecision
     * @param a
     * @param b
     * @return true iff a > epsilon
     */
    bool doubleIsPositive(double a) {
        return a > epsilon;
    }

    template<typename T>
    bool isNumberPositive(const T a) {
        if constexpr (std::is_floating_point<T>::value) {
            return a > epsilon;
        } else {
            return a > 0;
        }
    }

    template<typename T>
    bool isNumberNegative(const T a) {
        if constexpr (std::is_floating_point<T>::value) {
            return a < -epsilon;
        } else {
            return a < 0;
        }
    }

    template<typename T>
    bool areNumbersEqual(const T a, const T b) {
        if constexpr (std::is_floating_point<T>::value) {
            return std::abs(a - b) < epsilon;
        } else {
            return a == b;
        }
    }

    /**
     * Abstract class for functions to be used for parametric maxFlow
     */
    class flowFunction {
        /**
         * Returns the zero crossing of f with minimum x coordinate for which x >= minVal
         * @param minVal minimum value for the zero crossing
         * @return the zero crossing, std::numeric_limits<double>::infinity() if none exists
         */
        virtual double getNextZeroCrossing(double) const {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };

        [[nodiscard]] virtual double eval(double) const {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };

        friend flowFunction &operator+(flowFunction, const flowFunction&) {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };

        virtual flowFunction operator+=(const flowFunction&) {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };

        virtual flowFunction operator-=(const flowFunction&) {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };

        /**
         * When overriding this operator float inaccuracy should be considered
         */
        virtual bool operator==(const flowFunction &&) const {
            throw std::runtime_error("Flow function used directly, use class implementing flow function instead");
        };
    };


/**
 * Class to wrap "linear" functions of form f(x) = ax + b
 */
    class linearFlowFunction : public flowFunction {
    public:
        using FlowType = double;
        inline static constexpr double CONST_INF = INFTY;

        linearFlowFunction(double a, double b) : a(std::min(a, CONST_INF)), b(std::min(b, CONST_INF)) {};

        explicit linearFlowFunction(double b) : a(0), b(std::min(b, CONST_INF)) {};

        explicit linearFlowFunction() : a(0), b(0) {};

        double getNextZeroCrossing(double minVal) const override {
            if (a == 0) {
                if (b == 0) {
                    return minVal;
                } else {
                    return CONST_INF;
                }
            }
            if (b >= CONST_INF) return CONST_INF;
            double crossing = -(b / a);
            if (minVal > crossing) return CONST_INF;
            else return crossing;
        }

        std::string toString() {
            return std::to_string(b) + " " + std::to_string(a);
        }

        [[nodiscard]] double eval(double x) const override {
            if (b == CONST_INF) return CONST_INF;
            return x * a + b;
        }

        friend linearFlowFunction operator+(linearFlowFunction lhs, const linearFlowFunction &rhs) {
            lhs.a += rhs.a;
            if (lhs.b == CONST_INF || rhs.b == CONST_INF) {
                lhs.b = CONST_INF;
            } else {
                lhs.b += rhs.b;
            }
            return lhs;
        }

        friend linearFlowFunction operator-(linearFlowFunction lhs, const linearFlowFunction &rhs) {
            lhs.a -= rhs.a;
            if (lhs.b != CONST_INF) {
                lhs.b -= rhs.b;
            }
            return lhs;
        }

        linearFlowFunction operator-() const {
            return {-this->a, -this->b};
        }

        linearFlowFunction& operator=(const int rhs) {
            a = 0;
            b = rhs;
            return *this;
        }

        linearFlowFunction operator+=(const linearFlowFunction &rhs) {
            a += rhs.a;
            if (b == CONST_INF || rhs.b == CONST_INF) {
                b = CONST_INF;
            } else {
                b += rhs.b;
            }
            return *this;
        }

        linearFlowFunction operator-=(const linearFlowFunction &rhs) {
            a -= rhs.a;
            if (b != CONST_INF) {
                b -= rhs.b;
            }
            return *this;
        }

        bool operator==(const linearFlowFunction &rhs) const {
            return (std::abs(this->a - rhs.a) < epsilon) && (std::abs(this->b - rhs.b) < epsilon);
        }

        inline void serialize(IO::Serialization& serialize) const noexcept {
            serialize(a, b);
        }

        inline void deserialize(IO::Deserialization& deserialize) noexcept {
            deserialize(a, b);
        }

        inline friend std::ostream& operator<<(std::ostream& out, const linearFlowFunction& f) noexcept {
            out << f.a << "*x+" << f.b;
            return out;
        }

    private:
        double a, b;
    };


    constexpr ImplementationDetail::DistanceType Flow;
} // namespace pmf