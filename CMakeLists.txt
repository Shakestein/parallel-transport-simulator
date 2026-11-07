cmake_minimum_required(VERSION 3.16)
project(ParallelTransportSimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Main executable — header-only design, only main.cpp needs compiling
add_executable(parallel_transport
    src/main.cpp
)

target_include_directories(parallel_transport PRIVATE include)

# Optimization for Release builds
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    target_compile_options(parallel_transport PRIVATE /O2 /W4)
else()
    target_compile_options(parallel_transport PRIVATE -O2 -Wall -Wextra)
endif()

# Define _USE_MATH_DEFINES for MSVC to get M_PI
if(MSVC)
    target_compile_definitions(parallel_transport PRIVATE _USE_MATH_DEFINES)
endif()
