file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp)

add_library(polymath STATIC ${SRC_FILES})

target_include_directories(polymath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_compile_options(polymath PRIVATE -Wall -Werror)
