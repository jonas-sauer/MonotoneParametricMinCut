#pragma once

#include <initializer_list>
#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "../../../Helpers/Assert.h"

template<typename VALUE = std::uint32_t>
class EarliestArrivalEntry {

public:
    static constexpr size_t Size = 1;
    using Value = VALUE;
    using Type = EarliestArrivalEntry<Value>;

private:
    EarliestArrivalEntry(const Value value) :
        value(value) {
    }

public:
    EarliestArrivalEntry() :
        value(std::numeric_limits<Value>::max()) {
    }
    EarliestArrivalEntry(std::initializer_list<Value> defaultValues) :
        value(*defaultValues.begin()) {
    }

    static inline Type Max() noexcept {
        return Type();
    }

    inline Value operator()(const size_t i = 0) const noexcept {
        AssertMsg(i == 0, "Index i = " << i << " out of bounds (0, 0)!");
        return value;
    }

    inline Type& operator+=(const Type& other) noexcept {
        value += other.value;
        return *this;
    }

    inline Type& operator+=(const Value other) noexcept {
        value += other;
        return *this;
    }

    inline Type operator+(const Value other) const noexcept {
        return Type(value + other);
    }

    inline Type& maximize() noexcept {
        value = std::numeric_limits<Value>::max();
        return *this;
    }

    inline Type& minimize() noexcept {
        value = std::numeric_limits<Value>::min();
        return *this;
    }

    inline Type& minimize(const Type& other) noexcept {
        value = (value < other.value) ? (value) : (other.value);
        return *this;
    }

    inline Type& minimize(const Value other) noexcept {
        value = (value < other) ? (value) : (other);
        return *this;
    }

    inline Type& minimizeShifted(const Type& other) noexcept {
        value = (value < other.value) ? (value) : (other.value);
        return *this;
    }

    inline Value min(const int i = Size - 1) const noexcept {
        AssertMsg(i == 0, "Index i = " << i << " out of bounds (0, 0)!");
        return value;
    }

    inline Type& prepend(Value v = std::numeric_limits<Value>::max()) noexcept {
        value = value < v ? value : v;
        return *this;
    }

    inline Type& shift() noexcept {
        return prepend();
    }

    inline bool operator==(const Type& other) const noexcept {
        return value == other.value;
    }

    inline bool operator<(const Type& other) const noexcept {
        return value < other.value;
    }

    inline bool operator>(const Type& other) const noexcept {
        return other < *this;
    }

    inline long long byteSize() const noexcept {
        return sizeof(Type);
    }

    friend inline std::ostream& operator<<(std::ostream& out, const Type& l) noexcept {
        return out << "(" << l.value << ")";
    }

private:
    Value value;

};
