cmake_minimum_required(VERSION 3.14)

add_library(pydis MODULE)

add_subdirectory(calforce)
add_subdirectory(collision)
add_subdirectory(remesh)
add_subdirectory(util)
add_subdirectory(tests)

target_include_directories(pydis PRIVATE include)

install(TARGETS pydis DESTINATION ${CMAKE_SOURCE_DIR}/lib)
