cmake_minimum_required(VERSION 3.22)
project(ParametricMaxFlow)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

#Compiler settings
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "-pipe -march=native -Wfatal-errors")
set(CMAKE_CXX_FLAGS_DEBUG "-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math -DNDEBUG -O3")

add_executable(Benchmark Runnables/Benchmark.cpp)
add_executable(InstanceLoader Runnables/InstanceLoader.cpp)
add_executable(PrecisionExperiment Runnables/PrecisionExperiment.cpp)

# tests
option(TESTS "Build the tests" ON)
if (TESTS)
    enable_testing()
    add_subdirectory(libraries/googletest)
    include(GoogleTest)
    add_subdirectory(Test)
endif ()