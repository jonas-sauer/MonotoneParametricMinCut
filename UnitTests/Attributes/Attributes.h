#pragma once

#include <vector>

#include "../UnitTests.h"

#include "../../DataStructures/Attributes/Attributes.h"

namespace UnitTests {

class Attributes {

private:
    using IntWeight = Attribute<Weight, int>;
    using DoubleWeight = Attribute<Weight, double>;
    using VectorCapacity = Attribute<Capacity, std::vector<int>>;
    using IntCapacity = Attribute<Capacity, int>;
    using DoubleCapacity = Attribute<Capacity, double>;
    using BoolValid = Attribute<Valid, bool>;
    using DoubleValid = Attribute<Valid, double>;
    using StringValid = Attribute<Valid, std::string>;
    using DoubleFunction = Attribute<Function, double>;
    using IntFunction = Attribute<Function, int>;

public:
    inline void check() {
        checkAccess();
        checkForEach();
        checkHandle();
        checkAssignment();
        checkRecord();
        checkNameChanges();
    }

private:
    inline void checkAccess() {
        using AttributesType = ::Attributes<Meta::SortAttributes<List<IntWeight, VectorCapacity, DoubleWeight, BoolValid>>>;
        UnitTests::check(Meta::Equals<AttributesType::AttributeType<Weight>, int>(), "Attribute Weight should have type int, but has type ", Meta::type<AttributesType::AttributeType<Weight>>());
        UnitTests::check(Meta::Equals<AttributesType::AttributeType<Valid>, bool>(), "Attribute Valid should have type bool, but has type ", Meta::type<AttributesType::AttributeType<Valid>>());
        UnitTests::check(Meta::Equals<AttributesType::AttributeType<Capacity>, std::vector<int>>(), "Attribute Capacity should have type std::vector<int>, but has type ", Meta::type<AttributesType::AttributeType<Capacity>>());
        AttributesType attributes;
        UnitTests::check(attributes.NumberOfAttributes == 3, "Attributes should have 3 attributes, but has ", attributes.NumberOfAttributes);
        attributes.resize(2);
        UnitTests::check(attributes.size() == 2, "Attributes should have size 2, but has size ", attributes.size());
        UnitTests::check(attributes[Weight].size() == 2, "Weight vector should have size 2, but has size ", attributes[Weight].size());
        UnitTests::check(attributes[Weight][0] == 0, "Weight[0] should be 0 , but is", attributes[Weight][0]);
        UnitTests::check(attributes[Capacity][0].empty(), "Capacity[0] should be empty, but is not!");
        UnitTests::check(attributes.AttributeNames[0] == Weight, "The first Attribute should be named 'Weight' but is named '", attributeToString(attributes.AttributeNames[0]), "'!");
        attributes.setDefaultValue(Valid, true);
        attributes.resize(3);
        UnitTests::check(attributes[Weight][2] == 0, "Weight[2] should be 0 , but is", attributes[Weight][2]);
        UnitTests::check(attributes[Valid][1] == false, "Valid[1] should be false , but is", attributes[Valid][1]);
        UnitTests::check(attributes[Valid][2] == true, "Valid[2] should be true , but is", attributes[Valid][2]);
        attributes.set(Weight, 1, 42);
        UnitTests::check(attributes[Weight][1] == 42, "Weight[1] should be 42 , but is", attributes[Weight][1]);
        UnitTests::check(attributes[Weight][2] == 0, "Weight[2] should be 0 , but is", attributes[Weight][2]);
        attributes[Weight][1] = 43;
        UnitTests::check(attributes[Weight][1] == 43, "Weight[1] should be 43 , but is", attributes[Weight][1]);
        UnitTests::check(attributes[Weight][2] == 0, "Weight[2] should be 0 , but is", attributes[Weight][2]);
        UnitTests::check(attributes.get(Weight)[1] == 43, "Weight[1] should be 43 , but is", attributes.get(Weight)[1]);
        UnitTests::check(attributes.get(Weight, 1) == 43, "Weight[1] should be 43 , but is", attributes.get(Weight, 1));
        attributes.emplaceBack(Weight(77));
        UnitTests::check(attributes.size() == 4, "Attributes should have size 4, but has size ", attributes.size());
        UnitTests::check(attributes.get(Weight, 3) == 77, "Weight[3] should be 77, but is ", attributes.get(Weight, 3));
        UnitTests::check(attributes.back().get(Capacity).empty(), "Capacity.back() should be empty, but is not!");
        UnitTests::check(attributes[Valid][3], "Valid[3] should be true, but is not!");
        attributes.setDefaultValue(Weight, 42);
        attributes.emplaceBack(Valid(false), Capacity(std::vector<int>(3)));
        UnitTests::check(attributes.size() == 5, "Attributes should have size 5, but has size ", attributes.size());
        UnitTests::check(attributes.back().get(Weight) == 42, "Weight.back() should be 42, but is ", attributes.back().get(Weight));
        UnitTests::check(attributes.back().get(Capacity).size() == 3, "Capacity.back() should have size 3, but has size ", attributes.back().get(Capacity).size());
        UnitTests::check(!attributes[Valid][4], "Valid[4] should be false, but is not!");
        attributes.setDefaultValue(Capacity, std::vector<int>(2));
        attributes.emplaceBack(AnyAttribute(0));
        UnitTests::check(attributes.size() == 6, "Attributes should have size 6, but has size ", attributes.size());
        UnitTests::check(attributes.back().get(Weight) == 0, "Weight.back() should be 0, but is ", attributes.back().get(Weight));
        UnitTests::check(attributes.back().get(Capacity).size() == 2, "Capacity.back() should have size 2, but has size ", attributes.back().get(Capacity).size());
        UnitTests::check(!attributes[Valid][5], "Valid[5] should be false, but is not!");
    }

    inline void checkForEach() {
        using AttributesType = ::Attributes<Meta::SortAttributes<List<IntWeight, IntCapacity, DoubleWeight, DoubleValid>>>;
        AttributesType attributes;
        UnitTests::check(attributes.NumberOfAttributes == 3, "Attributes should have 3 attributes, but has ", attributes.NumberOfAttributes);
        UnitTests::check(attributes.empty(), "Attributes should be empty, but is not!");
        attributes.emplaceBack();
        UnitTests::check(!attributes.empty(), "Attributes should not be empty, but is!");
        UnitTests::check(attributes.size() == 1, "Attributes should have size 1, but has size ", attributes.size());
        attributes.forEach([](const auto& data){
            UnitTests::check(data.size() == 1, "Data vector should have size 1, but has size ", data.size());
            UnitTests::check(data[0] == 0, "Data[0] should be 0, but is ", data[0]);
        });
        attributes.forEach([](std::vector<double>& data){
            data[0] = 42.0;
        });
        UnitTests::check(attributes[Weight][0] == 0, "Weight[0] should be 0, but is ", attributes[Weight][0]);
        UnitTests::check(attributes[Capacity][0] == 0, "Capacity[0] should be 0, but is ", attributes[Capacity][0]);
        UnitTests::check(attributes[Valid][0] == 42, "Valid[0] should be 42, but is ", attributes[Valid][0]);
        attributes.forEach([](auto& data, const AttributeNameType name){
            data[0] = name;
        });
        UnitTests::check(attributes[Weight][0] == static_cast<int>(Weight), "Weight[0] should be 'Weight', but is ", attributes[Weight][0]);
        UnitTests::check(attributes[Capacity][0] == static_cast<int>(Capacity), "Capacity[0] should be 'Capacity', but is ", attributes[Capacity][0]);
        UnitTests::check(attributes[Valid][0] == static_cast<double>(Valid), "Valid[0] should be 'Valid', but is ", attributes[Valid][0]);
        const AttributesType& constAttributes = attributes;
        int sum = 0;
        constAttributes.forEach([&](const auto& data){
            sum += data[0];
        });
        UnitTests::check(sum == static_cast<int>(Weight + Capacity + Valid), "sum should be ", (Weight + Capacity + Valid), ", but is ", sum);
        attributes.emplaceBack(AnyAttribute(99.9));
        UnitTests::check(attributes.size() == 2, "Attributes should have size 2, but has size ", attributes.size());
        UnitTests::check(attributes[Weight][1] == 99, "Weight[1] should be 99, but is ", attributes[Weight][1]);
        UnitTests::check(attributes[Capacity][1] == 99, "Capacity[1] should be 99, but is ", attributes[Capacity][1]);
        UnitTests::check(attributes[Valid][1] == 99.9, "Valid[1] should be 99.9, but is ", attributes[Valid][1]);
    }

    inline void checkHandle() {
        using AttributesType = ::Attributes<Meta::SortAttributes<List<IntWeight, IntCapacity, DoubleWeight, DoubleValid>>>;
        using HandleType = AttributeHandle<AttributesType, uint32_t>;
        UnitTests::check(Meta::Equals<HandleType::AttributeType<Weight>, int>(), "Attribute Weight should have type int, but has type ", Meta::type<HandleType::AttributeType<Weight>>());
        AttributesType attributes(4);
        UnitTests::check(attributes.size() == 4, "Attributes should have size 4, but has size ", attributes.size());
        HandleType handle(attributes, 1);
        handle[Weight] = 42;
        UnitTests::check(attributes[Weight][1] == 42, "Weight[1] should be 42, but is ", attributes[Weight][1]);
        handle.set(Capacity, 10).set(Valid, 0.5);
        UnitTests::check(attributes[Capacity][1] == 10, "Capacity[1] should be 10, but is ", attributes[Capacity][1]);
        UnitTests::check(attributes[Valid][1] == 0.5, "Weight[1] should be 0.5, but is ", attributes[Valid][1]);
        attributes.setDefaultValue(Capacity, 1234);
        attributes.emplaceBack()[Valid] = 0.9876;
        UnitTests::check(attributes.back()[Weight] == 0, "Weight[4] should be 0, but is ", attributes.back()[Weight]);
        UnitTests::check(attributes.back()[Valid] == 0.9876, "Valid[4] should be 0.9876, but is ", attributes.back()[Valid]);
        UnitTests::check(attributes.back()[Capacity] == 1234, "Weight[4] should be 1234, but is ", attributes.back()[Capacity]);
        UnitTests::check(attributes[Weight][4] == 0, "Weight[4] should be 0, but is ", attributes[Weight][4]);
        UnitTests::check(attributes[Valid][4] == 0.9876, "Valid[4] should be 0.9876, but is ", attributes[Valid][4]);
        UnitTests::check(attributes[Capacity][4] == 1234, "Weight[4] should be 1234, but is ", attributes[Capacity][4]);
    }

    inline void checkAssignment() {
        using AttributesTypeA = ::Attributes<List<DoubleWeight, IntCapacity, DoubleValid>>;
        using AttributesTypeB = ::Attributes<List<IntCapacity, DoubleFunction, IntWeight>>;
        AttributesTypeA attributesA(4);
        attributesA[Weight][1] = 42.5;
        attributesA[Capacity][2] = 13;
        attributesA[Valid][3] = 17.2;
        AttributesTypeB attributesB(attributesA);
        UnitTests::check(attributesA.size() == 4, "AttributesA should have size 4, but has size ", attributesA.size());
        UnitTests::check(attributesB.size() == 4, "AttributesB should have size 4, but has size ", attributesB.size());
        UnitTests::check(attributesA[Weight][1] == 42.5, "A Weight[1] should be 42, but is ", attributesA[Weight][1]);
        UnitTests::check(attributesB[Weight][1] == 42, "B Weight[1] should be 42, but is ", attributesB[Weight][1]);
        UnitTests::check(attributesA[Capacity][2] == 13, "A Capacity[2] should be 13, but is ", attributesA[Capacity][2]);
        UnitTests::check(attributesB[Capacity][2] == 13, "B Capacity[2] should be 13, but is ", attributesB[Capacity][2]);
        UnitTests::check(attributesA[Valid][3] == 17.2, "A Valid[3] should be 17.2, but is ", attributesA[Valid][3]);
        UnitTests::check(attributesB[Function][3] == 0, "B Function[3] should be 0, but is ", attributesB[Function][3]);
    }

    inline void checkRecord() {
        using ListA = List<DoubleWeight, IntCapacity, DoubleValid>;
        using ListB = List<IntCapacity, DoubleFunction, IntWeight>;
        using AttributesTypeA = ::Attributes<ListA>;
        using AttributesTypeB = ::Attributes<ListB>;
        using AttributeRecordA = AttributeRecord<ListA>;
        UnitTests::check(AttributeRecordA::NumberOfAttributes == 3, "AttributeRecordA should have 3 attributes, but has ", AttributeRecordA::NumberOfAttributes);
        UnitTests::check(Meta::Equals<AttributeRecordA::AttributeType<Weight>, double>(), "AttributeRecordA Weight should have type double, but has type ", Meta::type<AttributeRecordA::AttributeType<Weight>>());
        UnitTests::check(Meta::Equals<AttributeRecordA::AttributeType<Valid>, double>(), "AttributeRecordA Valid should have type double, but has type ", Meta::type<AttributeRecordA::AttributeType<Valid>>());
        UnitTests::check(Meta::Equals<AttributeRecordA::AttributeType<Capacity>, int>(), "AttributeRecordA Capacity should have type int, but has type ", Meta::type<AttributeRecordA::AttributeType<Capacity>>());
        AttributesTypeA attributesA(4);
        attributesA[Weight][1] = 42.5;
        attributesA[Capacity][2] = 13;
        attributesA[Valid][3] = 17.2;
        AttributeRecordA recordA;
        UnitTests::check(recordA.NumberOfAttributes == 3, "recordA should have 3 attributes, but has ", recordA.NumberOfAttributes);
        UnitTests::check(recordA[Weight] == 0, "recordA Weight should be 0, but is ", recordA[Weight]);
        UnitTests::check(recordA[Capacity] == 0, "recordA Capacity should be 0, but is ", recordA[Capacity]);
        UnitTests::check(recordA.get(Valid) == 0, "recordA Valid should be 0, but is ", recordA.get(Valid));
        recordA[Weight] = 55.5;
        recordA[Capacity] = 7;
        recordA.set(Valid, 0.99);
        UnitTests::check(recordA[Weight] == 55.5, "recordA Weight should be 55.5, but is ", recordA[Weight]);
        UnitTests::check(recordA[Capacity] == 7, "recordA Capacity should be 7, but is ", recordA[Capacity]);
        UnitTests::check(recordA.get(Valid) == 0.99, "recordA Valid should be 0.99, but is ", recordA.get(Valid));
        attributesA.resize(8, recordA);
        UnitTests::check(recordA[Weight] == 55.5, "recordA Weight should be 55.5, but is ", recordA[Weight]);
        UnitTests::check(recordA[Capacity] == 7, "recordA Capacity should be 7, but is ", recordA[Capacity]);
        UnitTests::check(recordA.get(Valid) == 0.99, "recordA Valid should be 0.99, but is ", recordA.get(Valid));
        UnitTests::check(attributesA[Weight][3] == 0, "attributesA Weight[3] should be 0, but is ", attributesA[Weight][3]);
        UnitTests::check(attributesA[Capacity][3] == 0, "attributesA Capacity[3] should be 0, but is ", attributesA[Capacity][3]);
        UnitTests::check(attributesA.get(Valid)[3] == 17.2, "attributesA Valid[3] should be 17.2, but is ", attributesA.get(Valid)[3]);
        UnitTests::check(attributesA[Weight][4] == 55.5, "attributesA Weight[4] should be 55.5, but is ", attributesA[Weight][4]);
        UnitTests::check(attributesA[Capacity][4] == 7, "attributesA Capacity[4] should be 7, but is ", attributesA[Capacity][4]);
        UnitTests::check(attributesA.get(Valid)[4] == 0.99, "attributesA Valid[4] should be 0.99, but is ", attributesA.get(Valid)[4]);
        UnitTests::check(attributesA[Weight][7] == 55.5, "attributesA Weight[7] should be 55.5, but is ", attributesA[Weight][7]);
        UnitTests::check(attributesA[Capacity][7] == 7, "attributesA Capacity[7] should be 7, but is ", attributesA[Capacity][7]);
        UnitTests::check(attributesA.get(Valid)[7] == 0.99, "attributesA Valid[7] should be 0.99, but is ", attributesA.get(Valid)[7]);
        UnitTests::check(recordA.toString() == "Weight: 55.5, Capacity: 7, Valid: 0.99", "recordA to string should be 'Weight: 55.5, Capacity: 7, Valid: 0.99', but is '", recordA.toString(), "'!");
        AttributeRecordA recordB(attributesA, 1);
        UnitTests::check(recordB[Weight] == 42.5, "recordB Weight should be 42.5, but is ", recordB[Weight]);
        UnitTests::check(recordB[Capacity] == 0, "recordB Capacity should be 0, but is ", recordB[Capacity]);
        UnitTests::check(recordB.get(Valid) == 0, "recordB Valid should be 0, but is ", recordB.get(Valid));
        AttributesTypeB attributesB;
        attributesB.emplaceBack(recordB);
        attributesB.emplaceBack(recordA);
        attributesB.resize(4, recordA);
        attributesB.set(2, recordB);
        UnitTests::check(attributesB[Weight][0] == 42, "attributesB Weight[0] should be 42, but is ", attributesB[Weight][0]);
        UnitTests::check(attributesB[Weight][1] == 55, "attributesB Weight[1] should be 55, but is ", attributesB[Weight][1]);
        UnitTests::check(attributesB[Weight][2] == 42, "attributesB Weight[2] should be 42, but is ", attributesB[Weight][2]);
        UnitTests::check(attributesB[Weight][3] == 55, "attributesB Weight[3] should be 55, but is ", attributesB[Weight][3]);
        UnitTests::check(attributesB[Capacity][0] == 0, "attributesB Capacity[0] should be 0, but is ", attributesB[Capacity][0]);
        UnitTests::check(attributesB[Capacity][1] == 7, "attributesB Capacity[1] should be 7, but is ", attributesB[Capacity][1]);
        UnitTests::check(attributesB[Capacity][2] == 0, "attributesB Capacity[2] should be 0, but is ", attributesB[Capacity][2]);
        UnitTests::check(attributesB[Capacity][3] == 7, "attributesB Capacity[3] should be 7, but is ", attributesB[Capacity][3]);
        UnitTests::check(attributesB[Function][0] == 0, "attributesB Function[0] should be 0, but is ", attributesB[Function][0]);
        UnitTests::check(attributesB[Function][1] == 0, "attributesB Function[1] should be 0, but is ", attributesB[Function][1]);
        UnitTests::check(attributesB[Function][2] == 0, "attributesB Function[2] should be 0, but is ", attributesB[Function][2]);
        UnitTests::check(attributesB[Function][3] == 0, "attributesB Function[3] should be 0, but is ", attributesB[Function][3]);
        AttributeRecord<List<IntWeight>> AttributeRecordB(12345);
        UnitTests::check(AttributeRecordB[Weight] == 12345, "AttributeRecordB Weight should be 12345, but is ", AttributeRecordB[Weight]);
    }

    inline void checkNameChanges() {
        using ListA = List<DoubleWeight, StringValid, IntCapacity, IntFunction>;
        using ListB = List<StringValid, DoubleFunction, IntWeight, DoubleCapacity>;
        using AttributesTypeA = ::Attributes<ListA>;
        using AttributesTypeB = ::Attributes<ListB>;
        AttributesTypeA attributesA;
        attributesA.emplaceBack(Weight(17.4), Valid("Foo"), Capacity(42), Function(99));
        attributesA.emplaceBack(Function(98), Weight(23.9), Valid("Bar"), Capacity(43));
        attributesA.emplaceBack(Valid("XXX"), Capacity(24), Function(97), Weight(19.5));
        attributesA.emplaceBack(Weight(12.3), Function(96), Capacity(23), Valid("___"));
        const int* capacityPointer = &attributesA[Capacity][0];
        UnitTests::check(attributesA.size() == 4, "AttributesA should have size 4, but has size ", attributesA.size());
        UnitTests::check(attributesA[Valid][2] == "XXX", "attributesA Valid[2] should be 'XXX', but is '", attributesA[Valid][2], "'");
        UnitTests::check(attributesA[Capacity][2] == 24, "attributesA Capacity[2] should be 24, but is ", attributesA[Capacity][2]);
        UnitTests::check(attributesA[Function][2] == 97, "attributesA Function[2] should be 97, but is ", attributesA[Function][2]);
        UnitTests::check(attributesA[Weight][2] == 19.5, "attributesA Weight[2] should be 19.5, but is ", attributesA[Weight][2]);
        AttributesTypeB attributesB(std::move(attributesA), Capacity << Function, Weight << Capacity);
        UnitTests::check(attributesB.size() == 4, "AttributesA should have size 4, but has size ", attributesA.size());
        UnitTests::check(attributesB[Valid][2] == "XXX", "attributesB Valid[2] should be 'XXX', but is '", attributesB[Valid][2], "'");
        UnitTests::check(attributesB[Weight][2] == 24, "attributesB Weight[2] should be 24, but is ", attributesB[Weight][2]);
        UnitTests::check(attributesB[Capacity][2] == 97, "attributesB Capacity[2] should be 97, but is ", attributesB[Capacity][2]);
        UnitTests::check(attributesB[Function][0] == 0, "attributesB Function[0] should be 0, but is ", attributesB[Function][0]);
        UnitTests::check(attributesB[Function][1] == 0, "attributesB Function[1] should be 0, but is ", attributesB[Function][1]);
        UnitTests::check(attributesB[Function][2] == 0, "attributesB Function[2] should be 0, but is ", attributesB[Function][2]);
        UnitTests::check(attributesB[Function][3] == 0, "attributesB Function[3] should be 0, but is ", attributesB[Function][3]);
        const int* weightPointer = &attributesB[Weight][0];
        UnitTests::check(capacityPointer == weightPointer, "Moving Attributes of the same Type should not change memory addresses!");
    }

};

}
