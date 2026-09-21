# The tests must run on the interpreter that pyexadis was built for. The module
# carries that ABI in its name (pyexadis.cpython-39-x86_64-linux-gnu.so), so any
# other python reports it as missing, which reads as a broken build rather than a
# mismatch. pybind11 leaves the interpreter it built against in PYTHON_EXECUTABLE,
# so follow that one. A plain find_package(Python3) does not: it defaults to the
# newest python on the machine rather than the first on PATH, and on a cluster
# that is some unrelated system anaconda.
#
# Override with -DPython3_EXECUTABLE=/path/to/python3.
foreach(var Python3_EXECUTABLE PYTHON_EXECUTABLE)
    if(NOT OPENDIS_TEST_PYTHON AND ${var})
        set(OPENDIS_TEST_PYTHON "${${var}}")
    endif()
endforeach()

if(NOT OPENDIS_TEST_PYTHON)
    set(Python3_FIND_STRATEGY LOCATION)   # first on PATH, not newest on the machine
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    set(OPENDIS_TEST_PYTHON "${Python3_EXECUTABLE}")
endif()

message(STATUS "tests: interpreter ${OPENDIS_TEST_PYTHON}")

add_subdirectory(unit_tests)
add_subdirectory(full_runs)

# unit_tests are fast and self-contained; full_runs drive the scripts under examples/ and are
# slow. They carry the labels "unit" and "full_run" respectively:
#     ctest -L unit
#     ctest -L full_run
