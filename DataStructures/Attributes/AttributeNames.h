#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

#include "../../Helpers/Meta.h"

// Type used in templates, to identify an attribute
using AttributeNameType = uint32_t;

// Helper template class, so that the AttributeNameType of a value can be inferred from a parameter
template<AttributeNameType NAME, typename VALUE_TYPE>
class AttributeValueWrapper {

public:
    constexpr static AttributeNameType Name = NAME;
    using ValueType = VALUE_TYPE;
    using Type = AttributeValueWrapper<Name, ValueType>;

public:
    AttributeValueWrapper(const ValueType& value) : internalValue(value) {}

    inline operator const ValueType&() const noexcept {
        return internalValue;
    }

    inline const ValueType& value() const noexcept {
        return internalValue;
    }

private:
    const ValueType& internalValue;

};

// Helper template class, so that the AttributeNameType can be inferred from a parameter
template<AttributeNameType NAME>
class AttributeNameWrapper {

public:
    constexpr static AttributeNameType Name = NAME;
    using Type = AttributeNameWrapper<Name>;

public:
    constexpr AttributeNameWrapper() {}

    // Automatic conversion back to plain AttributeNameType, which will be used as template parameter
    constexpr inline operator AttributeNameType() const noexcept {
        return Name;
    }

    template<typename TYPE>
    inline AttributeValueWrapper<Name, TYPE> operator()(const TYPE& value) const noexcept {
        return AttributeValueWrapper<Name, TYPE>(value);
    }

};

// List of possible Attribute names, so that we can refer to Attributes by a meaningful name.
// For new Attributes you have to follow three simple steps:
//    1.  Add an instantiation of the template AttributeNameWrapper<> to the ImplementationDetail namespace
//    2.  Add a global constexpr variable, which will be used to refer to the attribute
//    3.  Add the name of the attribute as string to the array ImplementationDetail::AttributeNameStrings.

namespace ImplementationDetail {

    using WeightType = AttributeNameWrapper<0>;
    using LengthType = AttributeNameWrapper<1>;
    using DistanceType = AttributeNameWrapper<2>;
    using TravelTimeType = AttributeNameWrapper<3>;
    using FromVertexType = AttributeNameWrapper<4>;
    using ToVertexType = AttributeNameWrapper<5>;
    using ViaVertexType = AttributeNameWrapper<6>;
    using ReverseEdgeType = AttributeNameWrapper<7>;
    using FunctionType = AttributeNameWrapper<8>;
    using CapacityType = AttributeNameWrapper<9>;
    using ValidType = AttributeNameWrapper<10>;
    using IncomingEdgePointerType = AttributeNameWrapper<11>;
    using CoordinatesType = AttributeNameWrapper<12>;
    using BeginOutType = AttributeNameWrapper<13>;
    using OutDegreeType = AttributeNameWrapper<14>;
    using IncomingEdgesType = AttributeNameWrapper<15>;
    using SizeType = AttributeNameWrapper<16>;
    using UnknownType = AttributeNameWrapper<17>;
    // Ensure that Unknown is the last entry!

}

constexpr ImplementationDetail::WeightType Weight;
constexpr ImplementationDetail::LengthType Length;
constexpr ImplementationDetail::DistanceType Distance;
constexpr ImplementationDetail::TravelTimeType TravelTime;
constexpr ImplementationDetail::FromVertexType FromVertex;
constexpr ImplementationDetail::ToVertexType ToVertex;
constexpr ImplementationDetail::ViaVertexType ViaVertex;
constexpr ImplementationDetail::ReverseEdgeType ReverseEdge;
constexpr ImplementationDetail::FunctionType Function;
constexpr ImplementationDetail::CapacityType Capacity;
constexpr ImplementationDetail::ValidType Valid;
constexpr ImplementationDetail::IncomingEdgePointerType IncomingEdgePointer;
constexpr ImplementationDetail::CoordinatesType Coordinates;
constexpr ImplementationDetail::BeginOutType BeginOut;
constexpr ImplementationDetail::OutDegreeType OutDegree;
constexpr ImplementationDetail::IncomingEdgesType IncomingEdges;
constexpr ImplementationDetail::SizeType Size;
constexpr ImplementationDetail::UnknownType Unknown;
// Ensure that Unknown is the last entry!

namespace ImplementationDetail {

    constexpr const char* AttributeNameStrings[] = {
        /*  0 */ "Weight",
        /*  1 */ "Length",
        /*  2 */ "Distance",
        /*  3 */ "TravelTime",
        /*  4 */ "FromVertex",
        /*  5 */ "ToVertex",
        /*  6 */ "ViaVertex",
        /*  7 */ "ReverseEdge",
        /*  8 */ "Function",
        /*  9 */ "Capacity",
        /* 10 */ "Valid",
        /* 11 */ "IncomingEdgePointer",
        /* 12 */ "Coordinates",
        /* 13 */ "BeginOut",
        /* 14 */ "OutDegree",
        /* 15 */ "IncomingEdges",
        /* 16 */ "Size",
        /* 17 */ "Unknown"
        // Ensure that Unknown is the last entry!
    };

}

// Number of attribute names currently available. Used for automatic AttributeNameType to string conversion
constexpr std::size_t NumberOfAttributeNames = sizeof(ImplementationDetail::AttributeNameStrings) / sizeof(char*);

// Converts an attribute (as runtime parameter) to the corresponding string
inline constexpr const char* attributeToString(const AttributeNameType name) noexcept {
    return (name >= NumberOfAttributeNames) ? "Unknown" : ImplementationDetail::AttributeNameStrings[name];
}

// Converts an attribute (as template parameter) to the corresponding string
template<AttributeNameType ATTRIBUTE_NAME>
inline constexpr const char* attributeToString() noexcept {
    return (ATTRIBUTE_NAME >= NumberOfAttributeNames) ? "Unknown" : ImplementationDetail::AttributeNameStrings[ATTRIBUTE_NAME];
}

template<AttributeNameType ATTRIBUTE_NAME, typename TYPE>
inline TYPE getAttributeValue(const AttributeNameWrapper<ATTRIBUTE_NAME>, const TYPE& defaultValue) noexcept {
    return defaultValue;
}

template<AttributeNameType ATTRIBUTE_NAME, typename TYPE, AttributeNameType HEAD_NAME, typename HEAD_TYPE, typename... TAIL>
inline TYPE getAttributeValue(const AttributeNameWrapper<ATTRIBUTE_NAME>, const TYPE& defaultValue, const AttributeValueWrapper<HEAD_NAME, HEAD_TYPE>& head, const TAIL&... tail) noexcept {
    if constexpr (HEAD_NAME == ATTRIBUTE_NAME) {
        return head.value();
    } else {
        return getAttributeValue(AttributeNameWrapper<ATTRIBUTE_NAME>(), defaultValue, tail...);
    }
}

// Helper template class for changing the name of an attribute
template<AttributeNameType OLD_NAME, AttributeNameType NEW_NAME>
class NameChange {

public:
    constexpr static AttributeNameType OldName = OLD_NAME;
    constexpr static AttributeNameType NewName = NEW_NAME;
    using Type = NameChange<OldName, NewName>;

public:
    constexpr NameChange() {}

};

template<AttributeNameType OLD_NAME, AttributeNameType NEW_NAME>
constexpr inline NameChange<OLD_NAME, NEW_NAME> operator<<(const AttributeNameWrapper<NEW_NAME>, const AttributeNameWrapper<OLD_NAME>) noexcept {
    return NameChange<OLD_NAME, NEW_NAME>();
}

namespace ImplementationDetail {

    template<AttributeNameType NAME, typename LIST, typename... NAME_CHANGES>
    struct GetOldAttributeName;

    template<AttributeNameType NAME>
    struct GetOldAttributeName<NAME, Meta::List<>> : AttributeNameWrapper<NAME> {};

    template<AttributeNameType NAME, AttributeNameType OLD_NAME, AttributeNameType NEW_NAME, typename... NAME_CHANGES>
    struct GetOldAttributeName<NAME, Meta::List<NameChange<OLD_NAME, NEW_NAME>, NAME_CHANGES...>> : std::conditional_t<(NAME == OLD_NAME),
        UnknownType,
        GetOldAttributeName<NAME, Meta::List<NAME_CHANGES...>>> {
    };

    template<AttributeNameType NAME, AttributeNameType OLD_NAME, AttributeNameType NEW_NAME, typename... NAME_CHANGES, typename... LIST>
    struct GetOldAttributeName<NAME, Meta::List<LIST...>, NameChange<OLD_NAME, NEW_NAME>, NAME_CHANGES...> : std::conditional_t<(NAME == NEW_NAME),
        AttributeNameWrapper<OLD_NAME>,
        GetOldAttributeName<NAME, Meta::List<NameChange<OLD_NAME, NEW_NAME>, LIST...>, NAME_CHANGES...>> {
    };

}

template<AttributeNameType NAME, typename... NAME_CHANGES>
inline constexpr ImplementationDetail::GetOldAttributeName<NAME, Meta::List<>, NAME_CHANGES...> getOldAttributeName(const AttributeNameWrapper<NAME>, const NAME_CHANGES...) noexcept {
    return ImplementationDetail::GetOldAttributeName<NAME, Meta::List<>, NAME_CHANGES...>();
}
