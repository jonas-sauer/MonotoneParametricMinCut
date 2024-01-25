#pragma once

#include <initializer_list>
#include <iostream>
#include <vector>
#include <string>
#include <limits>

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/String/Enumeration.h"

template<size_t SIZE = 8, typename VALUE = std::uint32_t>
class ArrivalArrayEntry {

public:
    static constexpr size_t Size = SIZE;
    using Value = VALUE;
    using Parent = VALUE;
    using Type = ArrivalArrayEntry<SIZE, Value>;

private:
    ArrivalArrayEntry(const Type& base, const Value addend) {
        for (size_t i = 0; i < Size; i++) {
            values[i] = base.values[i] + addend;
        }
    }

public:
    ArrivalArrayEntry() {
        std::fill(values, values + Size, std::numeric_limits<Value>::max());
    }
    ArrivalArrayEntry(std::initializer_list<Value> defaultValues) {
        std::copy(defaultValues.begin(), defaultValues.end(), values);
    }

    inline void setParent(const Parent) noexcept {}

    inline Parent parent(const int) const noexcept {return std::numeric_limits<Parent>::max();}

    inline bool hasParent(const Parent) const noexcept {return false;}

    inline void removeParent(const Parent) noexcept {}

    static inline Type Max() noexcept {
        return Type();
    }

    inline Value operator()(const size_t i = 0) const noexcept {
        AssertMsg(i < Size, "Index i = " << i << " out of bounds (0, " << Size << ")!");
        return values[i];
    }

    inline Type& operator+=(const Type& other) noexcept {
        for (size_t i = 0; i < Size; i++) {
            values[i] += other.values[i];
        }
        return *this;
    }

    inline Type& operator+=(const Value other) noexcept {
        for (size_t i = 0; i < Size; i++) {
            values[i] += other;
        }
        return *this;
    }

    inline Type operator+(const Value other) const noexcept {
        return Type(*this, other);
    }

    inline Type& maximize() noexcept {
        std::fill(values, values + Size, std::numeric_limits<Value>::max());
        return *this;
    }

    inline Type& minimize() noexcept {
        std::fill(values, values + Size, std::numeric_limits<Value>::min());
        return *this;
    }

    inline Type& minimize(const Type& other) noexcept {
        for (size_t i = 0; i < Size; i++) {
            values[i] = (values[i] < other.values[i]) ? (values[i]) : (other.values[i]);
        }
        return *this;
    }

    inline Type& minimize(const Value other) noexcept {
        for (size_t i = 0; i < Size; i++) {
            values[i] = (values[i] < other) ? (values[i]) : (other);
        }
        return *this;
    }

    inline Type& minimizeShifted(const Type& other) noexcept {
        for (size_t i = 1; i < Size; i++) {
            values[i] = (values[i] < other.values[i - 1]) ? (values[i]) : (other.values[i - 1]);
        }
        return *this;
    }

    inline Value min(const size_t i = Size - 1) const noexcept {
        AssertMsg(i < Size, "Index i = " << i << " out of bounds (0, " << Size << ")!");
        Value minValue = values[0];
        for (size_t j = 1; j <= i; j++) {
            if (minValue > values[i]) {
                minValue = values[i];
            }
        }
        return minValue;
    }

    inline Type& prepend(Value v = std::numeric_limits<Value>::max()) noexcept {
        memmove(values + 1, values, (Size - 1) * sizeof(Value));
        values[0] = v;
        return *this;
    }

    inline Type& shift() noexcept {
        return prepend();
    }

    inline bool operator==(const Type& other) const noexcept {
        for (size_t i = 0; i < Size; i++) {
            if (values[i] != other.values[i]) {
                return false;
            }
        }
        return true;
    }

    inline bool operator<(const Type& other) const noexcept {
        for (size_t i = 0; i < Size; i++) {
            if (values[i] < other.values[i]) {
                return true;
            }
        }
        return false;
    }

    inline bool operator>(const Type& other) const noexcept {
        return other < *this;
    }

    inline long long byteSize() const noexcept {
        return sizeof(Type);
    }

    friend inline std::ostream& operator<<(std::ostream& out, const Type& l) noexcept {
        Enumeration e;
        for (const Value i : l.values) {
            e << String::secToTime(i) << sep;
        }
        return out << "(" << e << ")";
    }

private:
    Value values[Size];

};
