cmake_minimum_required(VERSION 3.22)
project(RoutingFramework)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

#Compiler settings
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "-pipe -march=native -Wfatal-errors")
set(CMAKE_CXX_FLAGS_DEBUG "-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math -DNDEBUG -O3")

#Libraries
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_executable(Benchmark Runnables/Benchmark.cpp)

#Shells
function(add_shell name)
    add_executable(${ARGV})
endfunction()

set(SHELLS FlowShell InstanceLoader)
foreach(f ${SHELLS})
    add_shell(${f} Runnables/${f}.cpp)
endforeach()

add_custom_target(Shells)
add_dependencies(Shells ${SHELLS})

# tests
option(TESTS "Build the tests" ON)
if (TESTS)
    enable_testing()
    add_subdirectory(libraries/googletest)
    include(GoogleTest)
    add_subdirectory(test)
endif ()

# benchmark
option(BENCHMARK "Additional Code" OFF)
if (BENCHMARK)
    add_subdirectory(libraries/googlebenchmark)
    add_executable(test_benchmark)
    target_include_directories(test_benchmark
            PRIVATE
            ${PROJECT_SOURCE_DIR}/include
            )
    target_link_libraries(test_benchmark
            benchmark::benchmark)
endif ()
