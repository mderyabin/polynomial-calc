cmake_minimum_required(VERSION 3.11)

project(PolynomialCalc)
include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)

option(USE_SERIAL "Use serialization" OFF) 

set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/bin)

if(USE_SERIAL)
ExternalProject_Add(cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    SOURCE_DIR "${CMAKE_BINARY_DIR}/cereal"
    BUILD_COMMAND cmake -E echo "Skipping build step for cereal."
    INSTALL_COMMAND cmake -E echo "Skipping install step for cereal."
)
include_directories(${CMAKE_BINARY_DIR}/cereal/include)
include_directories(${CMAKE_BINARY_DIR}/cereal/include/cereal)
endif()

SET(CMAKE_CXX_FLAGS "-O3 -Wall -frtti -Wextra")

add_subdirectory(src)

file (GLOB EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
foreach (app ${EXAMPLES_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    target_link_libraries(${exe} PUBLIC polymath)
    if(!USE_SERIAL)
        target_compile_definitions(${exe} PUBLIC "N_USE_SERIAL")
    endif()
endforeach()