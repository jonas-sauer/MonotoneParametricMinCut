#pragma once

#include <iostream>
#include <vector>
#include <array>

#include "../UnitTests.h"

#include "../../Helpers/Vector/Vector.h"
#include "../../Helpers/IO/Serialization.h"
#include "../../Helpers/FileSystem/FileSystem.h"

namespace UnitTests {

class IO {

private:
    class Serializable {

    public:
        inline void serialize(::IO::Serialization& serialize) const noexcept {
            serialize(i, d, s, v, a);
        }

        inline void deserialize(::IO::Deserialization& deserialize) noexcept {
            deserialize(i, d, s, v, a);
        }

        int i;
        double d;
        std::string s;
        std::vector<std::vector<bool>> v;
        std::array<bool, 7> a;

    };

    class NotSerializable {

    public:
        int i;
        float f;
        uint8_t u;
        bool b;

    };

public:
    inline void check() {
        std::vector<Serializable> v1;
        v1.resize(3);
        v1[0].i = 17;
        v1[1].i = 42;
        v1[2].i = 999;
        v1[0].d = 0.99;
        v1[1].d = 42.5;
        v1[1].s = "Hallo";
        v1[2].s = " World!\n\tFoo\nBar";
        v1[0].v = std::vector<std::vector<bool>>{{false, true, false}, {false, false, true, false}, {true, true, false, true, false}};
        v1[1].v = std::vector<std::vector<bool>>{{false, true, true, false}, {true, true, true, true, false, true}, {false, true, false}, {true, false, true}, {true, true}};
        v1[2].v = std::vector<std::vector<bool>>{{false, false}, {false, true}, {true, false}, {true, true}, {true, false, false}, {true, false, true}, {true, true, false}, {true, true, true}};
        v1[1].a = std::array<bool, 7>{false, true, false, false, true, true, false};
        ::IO::serialize("UnitTest.testFile", v1);
        std::vector<Serializable> v2;
        ::IO::deserialize("UnitTest.testFile", v2);
        UnitTests::check(v2[0].i == 17, "v2[0].i should be 17, but is ", v2[0].i);
        UnitTests::check(v2[1].i == 42, "v2[1].i should be 42, but is ", v2[1].i);
        UnitTests::check(v2[2].i == 999, "v2[2].i should be 999, but is ", v2[2].i);
        UnitTests::check(v2[0].d == 0.99, "v2[0].i should be 0.99, but is ", v2[0].d);
        UnitTests::check(v2[1].d == 42.5, "v2[1].i should be 42.5, but is ", v2[1].d);
        UnitTests::check(v2[1].s == "Hallo", "v2[1].i should be '', but is ", v2[1].s);
        UnitTests::check(v2[2].s == " World!\n\tFoo\nBar", "v2[2].i should be ' World!\n\tFoo\nBar', but is ", v2[2].s);
        UnitTests::check(Vector::equals(v2[0].v, v1[0].v), "v2[0].v should be equal to v1[0].v, but is not!");
        UnitTests::check(Vector::equals(v2[1].v, v1[1].v), "v2[1].v should be equal to v1[1].v, but is not!");
        UnitTests::check(Vector::equals(v2[2].v, v1[2].v), "v2[2].v should be equal to v1[2].v, but is not!");
        UnitTests::check(v2[1].a[0] == false, "v2[1].a[0] should be false, but is ", v2[1].a[0]);
        UnitTests::check(v2[1].a[1] == true, "v2[1].a[1] should be false, but is ", v2[1].a[1]);
        UnitTests::check(v2[1].a[2] == false, "v2[1].a[2] should be false, but is ", v2[1].a[2]);
        UnitTests::check(v2[1].a[3] == false, "v2[1].a[3] should be false, but is ", v2[1].a[3]);
        UnitTests::check(v2[1].a[4] == true, "v2[1].a[4] should be false, but is ", v2[1].a[4]);
        UnitTests::check(v2[1].a[5] == true, "v2[1].a[5] should be false, but is ", v2[1].a[5]);
        UnitTests::check(v2[1].a[6] == false, "v2[1].a[6] should be false, but is ", v2[1].a[6]);

        std::vector<NotSerializable> v3;
        v3.resize(4);
        v3[0].i = 19;
        v3[1].i = 199;
        v3[2].i = 19999;
        v3[1].f = -1.1111;
        v3[2].f = 0.0909090909;
        v3[3].f = 10000000000;
        v3[0].u = 42;
        v3[2].u = 84;
        v3[1].b = true;
        v3[3].b = true;
        ::IO::serialize("UnitTest.testFile", v3);
        std::vector<NotSerializable> v4;
        ::IO::deserialize("UnitTest.testFile", v4);
        UnitTests::check(v4[0].i == 19, "v4[0].i should be 19, but is ", v4[0].i);
        UnitTests::check(v4[1].i == 199, "v4[1].i should be 199, but is ", v4[1].i);
        UnitTests::check(v4[2].i == 19999, "v4[2].i should be 19999, but is ", v4[2].i);
        UnitTests::check(v4[1].f == v3[1].f, "v4[1].i should be -1.1111, but is ", v4[1].f);
        UnitTests::check(v4[2].f == v3[2].f, "v4[2].i should be 0.0909090909, but is ", v4[2].f);
        UnitTests::check(v4[3].f == v3[3].f, "v4[3].i should be 10000000000, but is ", v4[3].f);
        UnitTests::check(v4[0].u == 42, "v4[0].i should be 42, but is ", v4[0].u);
        UnitTests::check(v4[2].u == 84, "v4[2].i should be 84, but is ", v4[2].u);
        UnitTests::check(v4[1].b == true, "v4[1].i should be true, but is ", v4[1].b);
        UnitTests::check(v4[3].b == true, "v4[3].i should be true, but is ", v4[3].b);

        FileSystem::deleteFile("UnitTest.testFile");
    }

};

}
