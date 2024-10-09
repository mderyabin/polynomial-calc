file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp)

add_library(polymath STATIC ${SRC_FILES})

target_include_directories(polymath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DERICTORY ${CMAKE_BINARY_DIR}/lib)

target_link_libraries(polymath PUBLIC gmp)
target_link_libraries(polymath PUBLIC m)
target_link_libraries(polymath PUBLIC ntl)
target_link_libraries(polymath PUBLIC fftw3)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)