cmake_minimum_required(VERSION 3.11)

project(PolynomialCalc)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/external/cereal/include)
include_directories(${CMAKE_SOURCE_DIR}/external/cereal/include/cereal)

SET(CMAKE_CXX_FLAGS "-O3 -Wall -frtti -Wextra")
include_directories("/opt/homebrew/include")
link_directories("/opt/homebrew/lib/")

add_subdirectory(src)

# SET(RUN_HAVE_STD_REGEX 0)
# SET(RUN_HAVE_POSIX_REGEX 0)
# SET(HAVE_STEADY_CLOCK 0)
# SET(HAVE_POSIX_REGEX 0)
# SET(CMAKE_CROSSCOMPILING 1)

include(FetchContent)
find_package(Git REQUIRED)
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        benchmark
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG aa59d40f8822358d012e854488c2e52857b8a1f0
        PATCH_COMMAND ""
)
FetchContent_MakeAvailable(benchmark)

add_subdirectory(benchmark)


file (GLOB EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
foreach (app ${EXAMPLES_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DERICTORY}/examples)
    target_link_libraries(${exe} PUBLIC polymath)
endforeach()