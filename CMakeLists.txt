cmake_minimum_required(VERSION 3.11)

project(PolynomialCalc)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)

# add_executable(polycalc main.cpp)

# targer_link_libtlraries(polycalc PUBLIC polymath)

file (GLOB EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
foreach (app ${EXAMPLES_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    target_link_libraries(${exe} PUBLIC polymath)
endforeach()