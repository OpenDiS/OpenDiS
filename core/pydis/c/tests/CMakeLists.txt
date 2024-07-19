cmake_minimum_required(VERSION 3.14)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
add_executable(test_calforce)

target_sources(test_calforce PRIVATE test_calforce.c)

add_compile_definitions("SOURCE_ROOT=${CMAKE_SOURCE_DIR}")

target_link_libraries(test_calforce ${PYTHON_LIBRARIES})