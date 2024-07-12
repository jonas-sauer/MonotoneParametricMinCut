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
find_package(Cairo REQUIRED)
include_directories(${CAIRO_INCLUDE_DIR})
find_package(OpenMP REQUIRED)
add_subdirectory(MapServer/mongoose)

add_executable(ParametricMaxFlowBenchmark Runnables/ParametricMaxFlowBenchmark.cpp)

#HalloWorldMapServer
add_executable(HalloWorldMapServer MapServer/HalloWorldMapServer.cpp)
target_link_libraries(HalloWorldMapServer PUBLIC OpenMP::OpenMP_CXX PRIVATE Threads::Threads dl mongoose)

#Shells
function(add_shell name)
    add_executable(${ARGV})
    target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX ${CAIRO_LIBRARIES})
endfunction()

set(SHELLS Assignment BenchmarkShell BikeSharing CH DelayExperimentsShell FlowShell IsochroneShell MultimodalCoverageShell NetworkBuilder NetworkOrderShell PartitionShell QueriesShell ULTRA UltraAssignment ULTRAExperimentsShell UnrestrictedWalkingShell ValidationShell)
foreach(f ${SHELLS})
    add_shell(${f} Runnables/${f}.cpp)
endforeach()

add_custom_target(Shells)
add_dependencies(Shells ${SHELLS})

#MapServers
set(MAP_SERVERS DemandMapServer RaptorMapServer)
function(add_server name)
    add_executable(${ARGV})
    target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX PRIVATE Threads::Threads dl mongoose ${CAIRO_LIBRARIES})
endfunction()

foreach(f ${MAP_SERVERS})
    add_server(${f} Runnables/${f}.cpp)
endforeach()

add_custom_target(MapServers)
add_dependencies(MapServers ${MAP_SERVERS})

#UnitTests
add_executable(UnitTests Runnables/UnitTests.cpp)
target_link_libraries(UnitTests PUBLIC OpenMP::OpenMP_CXX)

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
