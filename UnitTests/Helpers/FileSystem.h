#pragma once

#include <iostream>
#include <vector>
#include <array>

#include "../UnitTests.h"

#include "../../Helpers/FileSystem/FileSystem.h"

namespace UnitTests {

class FileSystem {

public:
    inline void check() {
        UnitTests::check(::FileSystem::isAbsolutePath("/temp"), "Path '/temp' should be absolute, but is not!");
        UnitTests::check(::FileSystem::isAbsolutePath("/temp/foo"), "Path '/temp/foo' should be absolute, but is not!");
        UnitTests::check(::FileSystem::isAbsolutePath("/temp/"), "Path '/temp/' should be absolute, but is not!");
        UnitTests::check(::FileSystem::isAbsolutePath("/"), "Path '/' should be absolute, but is not!");
        UnitTests::check(::FileSystem::isAbsolutePath("//"), "Path '//' should be absolute, but is not!");
        UnitTests::check(!::FileSystem::isAbsolutePath(""), "Path '' should not be absolute, but is!");
        UnitTests::check(!::FileSystem::isAbsolutePath("temp"), "Path 'temp' should not be absolute, but is!");
        UnitTests::check(!::FileSystem::isAbsolutePath("a/b"), "Path 'a/b' should not be absolute, but is!");
        UnitTests::check(!::FileSystem::isAbsolutePath("./"), "Path './' should not be absolute, but is!");

        checkGetParentDirectory("/temp/foo", "/temp");
        checkGetParentDirectory("/temp/foo/", "/temp");
        checkGetParentDirectory("/", "");
        checkGetParentDirectory("/foo", "");
        checkGetParentDirectory("temp/foo", "temp");
        checkGetParentDirectory("temp/foo/", "temp");
        checkGetParentDirectory("temp", "");
        checkGetParentDirectory("", "");

        checkExtendPath("/temp/", "foo/bar", "/temp/foo/bar");
        checkExtendPath("/temp/", "/foo/bar", "/foo/bar");
        checkExtendPath("/temp/", "./", "/temp/");
        checkExtendPath("/temp/", ".", "/temp/");
        checkExtendPath("/temp", "foo/bar", "/temp/foo/bar");
        checkExtendPath("/temp", "/foo/bar", "/foo/bar");
        checkExtendPath("/temp", "./", "/temp");
        checkExtendPath("/temp", ".", "/temp");
        checkExtendPath("temp", "foo/bar", "temp/foo/bar");
        checkExtendPath("temp", "/foo/bar", "/foo/bar");
        checkExtendPath("temp", "./", "temp");
        checkExtendPath("temp", ".", "temp");
        checkExtendPath("/temp/", "./././foo", "/temp/foo");
        checkExtendPath("/temp/foo/", "..", "/temp");
        checkExtendPath("/temp/foo", "..", "/temp");
        checkExtendPath("/temp/", "..", "");
        checkExtendPath("temp/", "..", "");
        checkExtendPath("/temp/foo/", "../", "/temp");
        checkExtendPath("/temp/foo", "../", "/temp");
        checkExtendPath("/temp/", "../", "");
        checkExtendPath("temp/", "../", "");
        checkExtendPath("/temp/foo/x/", "./../.", "/temp/foo");
        checkExtendPath("/temp/foo/x/", "./.././", "/temp/foo");
        checkExtendPath("/temp/foo/x/", "./.././bar", "/temp/foo/bar");
        checkExtendPath("/temp/foo/x", "./../.", "/temp/foo");
        checkExtendPath("/temp/foo/x", "./.././", "/temp/foo");
        checkExtendPath("/temp/foo/x", "./.././bar", "/temp/foo/bar");
        checkExtendPath("/temp/foo/", "./../.", "/temp");
        checkExtendPath("/temp/foo/", "./.././", "/temp");
        checkExtendPath("/temp/foo/", "./.././bar", "/temp/bar");
        checkExtendPath("/temp/foo", "./../.", "/temp");
        checkExtendPath("/temp/foo", "./.././", "/temp");
        checkExtendPath("/temp/foo", "./.././bar", "/temp/bar");
        checkExtendPath("", "", "");
    }

private:
    inline void checkGetParentDirectory(const std::string& path, const std::string& parent) const noexcept {
        UnitTests::check(::FileSystem::getParentDirectory(path) == parent, "Parent directory of path '", path, "' should be '", parent, "', but is '", ::FileSystem::getParentDirectory(path), "'!");
    }

    inline void checkExtendPath(const std::string& base, const std::string& file, const std::string& result) const noexcept {
        UnitTests::check(::FileSystem::extendPath(base, file) == result, "Extending path '", base, "' with '", file, "' should be '", result, "', but is '", ::FileSystem::extendPath(base, file), "'!");
    }

};

}
