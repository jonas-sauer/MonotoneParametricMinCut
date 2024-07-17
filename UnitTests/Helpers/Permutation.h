#pragma once

#include <vector>
#include <string>
#include <sstream>

#include "../UnitTests.h"

#include "../../Helpers/Vector/Permutation.h"
#include "../../Helpers/String/Enumeration.h"

namespace UnitTests {

class Permutation {

public:
    inline void check() {
        std::vector<char> chars1({'e', 'd', 'f', 'b', 'g', 'a', 'c'});

        Order order1(Construct::Sort, chars1);
        UnitTests::check(getString(order1) == "5 3 6 1 0 2 4", "The Order 1 should be '5 3 6 1 0 2 4', but is '", getString(order1), "'!");

        ::Permutation permutation1(Construct::Invert, order1);
        UnitTests::check(getString(permutation1) == "4 3 5 1 6 0 2", "The Permutation 1 should be '4 3 5 1 6 0 2', but is '", getString(permutation1), "'!");

        Order order2(Construct::Invert, permutation1);
        UnitTests::check(getString(order2) == "5 3 6 1 0 2 4", "The Order 2 should be '5 3 6 1 0 2 4', but is '", getString(order2), "'!");

        ::Permutation permutation2(Construct::Invert, std::move(order2));
        UnitTests::check(getString(permutation2) == "4 3 5 1 6 0 2", "The Permutation 2 should be '4 3 5 1 6 0 2', but is '", getString(permutation2), "'!");
        UnitTests::check(order2.empty(), "The Order 2 should be empty but is '", getString(order2), "'!");

        std::vector<char> chars2 = permutation1.getPermuted(chars1);
        UnitTests::check(getString(permutation1) == "4 3 5 1 6 0 2", "The Permutation 1 should be '4 3 5 1 6 0 2', but is '", getString(permutation1), "'!");
        UnitTests::check(getString(chars1) == "e d f b g a c", "The Chars 1 should be 'e d f b g a c', but is '", getString(chars1), "'!");
        UnitTests::check(getString(chars2) == "a b c d e f g", "The Chars 2 should be 'a b c d e f g', but is '", getString(chars2), "'!");

        std::vector<char> chars3 = chars1;
        permutation1.permutate(chars3);
        UnitTests::check(getString(permutation1) == "4 3 5 1 6 0 2", "The Permutation 1 should be '4 3 5 1 6 0 2', but is '", getString(permutation1), "'!");
        UnitTests::check(getString(chars3) == "a b c d e f g", "The Chars 3 should be 'a b c d e f g', but is '", getString(chars3), "'!");

        std::vector<char> chars4 = order1.getOrdered(chars1);
        UnitTests::check(getString(order1) == "5 3 6 1 0 2 4", "The Order 1 should be '5 3 6 1 0 2 4', but is '", getString(order1), "'!");
        UnitTests::check(getString(chars1) == "e d f b g a c", "The Chars 1 should be 'e d f b g a c', but is '", getString(chars1), "'!");
        UnitTests::check(getString(chars4) == "a b c d e f g", "The Chars 4 should be 'a b c d e f g', but is '", getString(chars4), "'!");

        std::vector<char> chars5 = chars1;
        order1.order(chars5);
        UnitTests::check(getString(order1) == "5 3 6 1 0 2 4", "The Order 1 should be '5 3 6 1 0 2 4', but is '", getString(order1), "'!");
        UnitTests::check(getString(chars5) == "a b c d e f g", "The Chars 5 should be 'a b c d e f g', but is '", getString(chars5), "'!");

        std::vector<int> ints = Vector::id<int>(chars1.size());
        for (int& c : ints) {
            c = permutation1.permutate(c);
        }
        UnitTests::check(getString(permutation1) == "4 3 5 1 6 0 2", "The Permutation 1 should be '4 3 5 1 6 0 2', but is '", getString(permutation1), "'!");
        UnitTests::check(getString(ints) == "4 3 5 1 6 0 2", "The Chars 6 should be '4 3 5 1 6 0 2', but is '", getString(ints), "'!");
    }

private:

    template<typename T>
    inline std::string getString(const std::vector<T>& vec) {
        Enumeration result;
        for (const T& t : vec) {
            result << t << Sep(" ");
        }
        return result;
    }

};

}
