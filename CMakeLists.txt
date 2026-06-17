cmake_minimum_required(VERSION 3.20)

option(FWI_WITH_CUDA "Build with CUDA acceleration" ON)
option(FWI_WITH_PYTHON "Build Python bindings" ON)
option(FWI_BUILD_TESTS "Build test programs" ON)
option(FWI_BUILD_EXAMPLES "Build example programs" ON)

if(FWI_WITH_CUDA)
    project(seismic_fwi_engine LANGUAGES CXX CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
else()
    project(seismic_fwi_engine LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CUDA enabled: ${FWI_WITH_CUDA}")
message(STATUS "Python bindings: ${FWI_WITH_PYTHON}")
message(STATUS "Build tests: ${FWI_BUILD_TESTS}")
message(STATUS "Build examples: ${FWI_BUILD_EXAMPLES}")

include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)

if(FWI_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(FWI_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
