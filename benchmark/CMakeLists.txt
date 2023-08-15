
file (GLOB EXAMPLES_SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
foreach (app ${EXAMPLES_SRC_FILES})
    get_filename_component ( exe ${app} NAME_WE )
    add_executable ( ${exe} ${app} )
    target_include_directories( ${exe} PUBLIC ${benchmark_SOURCE_DIR}/include )
    set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DERICTORY}/benchmark)
    target_link_libraries(${exe} PRIVATE polymath)
    target_link_libraries(${exe} PRIVATE benchmark::benchmark)

endforeach()