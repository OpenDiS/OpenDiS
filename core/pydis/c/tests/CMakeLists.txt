cmake_minimum_required(VERSION 3.14)

find_package(PythonLibs)
if(PYTHONLIBS_FOUND)
    add_executable(test_calforce)
    include_directories(${PYTHON_INCLUDE_DIRS})

    target_sources(test_calforce PRIVATE test_calforce.c)

    add_compile_definitions("SOURCE_ROOT=${CMAKE_SOURCE_DIR}")

    target_link_libraries(test_calforce ${PYTHON_LIBRARIES})
else()
    message("skips building core/pydis/c/tests/")
endif(PYTHONLIBS_FOUND)


