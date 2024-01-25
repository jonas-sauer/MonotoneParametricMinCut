#pragma once

#include <initializer_list>
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <cstdint>

#include <stdio.h>
#include <string.h>
#include <emmintrin.h>
#include <smmintrin.h>

#include "../../../Helpers/Assert.h"
#include "../../../Helpers/Types.h"
#include "../../../Helpers/String/String.h"
#include "../../../Helpers/String/Enumeration.h"

template<size_t SIZE = 8>
class ArrivalSIMDEntry {

public:
    static constexpr size_t N = SIZE >> 2;
    static constexpr size_t Size = N << 2;
    using Value = std::uint32_t;
    using Type = ArrivalSIMDEntry<SIZE>;

public:
    ArrivalSIMDEntry() {
        std::fill(values, values + Size, std::numeric_limits<Value>::max());
    }
    ArrivalSIMDEntry(std::initializer_list<Value> defaultValues) {
        std::copy(defaultValues.begin(), defaultValues.end(), values);
    }

    static inline Type Max() noexcept {
        return Type();
    }

    inline Value operator()(const size_t i = 0) const noexcept {
        AssertMsg(i < Size, "Index i = " << i << " out of bounds (0, " << Size << ")!");
        return values[i];
    }

    inline Type& operator+=(const Type& other) noexcept {
        for (size_t i = 0; i < N; i++) {
            value(i) = _mm_add_epi32(value(i), other.value(i));
        }
        return *this;
    }

    inline Type& operator+=(const Value other) noexcept {
        __m128i o = _mm_set1_epi32(other);
        for (size_t i = 0; i < N; i++) {
            value(i) = _mm_add_epi32(value(i), o);
        }
        return *this;
    }

    inline Type operator+(const Value other) const noexcept {
        return (Type(*this) += other);
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
        for (size_t i = 0; i < N; i++) {
            value(i) = _mm_min_epu32(value(i), other.value(i));
        }
        return *this;
    }

    inline Type& minimize(const Value other) noexcept {
        __m128i o = _mm_set1_epi32(other);
        for (size_t i = 0; i < N; i++) {
            value(i) = _mm_min_epu32(value(i), o);
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

    inline Type& prepend(const Value v = std::numeric_limits<Value>::max()) noexcept {
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
    Value values[Size] __attribute__ ((aligned (32)));

    inline __m128i& value(const size_t i) {return reinterpret_cast<__m128i*>(values)[i];}
    inline const __m128i& value(const size_t i) const {return reinterpret_cast<const __m128i*>(values)[i];}

};
