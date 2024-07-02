cmake_minimum_required(VERSION 3.25)
project(Mongoose)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS " -march=native -Wno-deprecated-register")
set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math -DNDEBUG -O3")
set(CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_RELEASE})

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

set_source_files_properties(mongoose.c PROPERTIES LANGUAGE CXX)
add_library(mongoose SHARED mongoose.c)
target_include_directories(mongoose PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(mongoose PUBLIC OpenMP::OpenMP_CXX PRIVATE Threads::Threads)