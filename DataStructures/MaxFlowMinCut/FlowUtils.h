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

        linearFlowFunction(double a, double b) : a(a), b(b) {};

        explicit linearFlowFunction(double b) : a(0), b(b) {};

        explicit linearFlowFunction() : a(0), b(0) {};

        double getNextZeroCrossing(double minVal) const override {
            if (a == 0) {
                if (b == 0) {
                    return minVal;
                } else {
                    return std::numeric_limits<double>::infinity();
                }
            }
            double crossing = -(b / a);
            if (minVal > crossing) return std::numeric_limits<double>::infinity();
            else return crossing;
        }

        std::string toString() {
            return "f(x) = " + std::to_string(a) + "*x + " + std::to_string(b);
        }

        [[nodiscard]] double eval(double x) const override {
            return x * a + b;
        }

        friend linearFlowFunction operator+(linearFlowFunction lhs, const linearFlowFunction &rhs) {
            lhs.a += rhs.a;
            lhs.b += rhs.b;
            return lhs;
        }

        friend linearFlowFunction operator-(linearFlowFunction lhs, const linearFlowFunction &rhs) {
            lhs.a -= rhs.a;
            lhs.b -= rhs.b;
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
            b += rhs.b;
            return *this;
        }

        linearFlowFunction operator-=(const linearFlowFunction &rhs) {
            a -= rhs.a;
            b -= rhs.b;
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

    private:
        double a, b;
    };


    constexpr ImplementationDetail::DistanceType Flow;
} // namespace pmf