# Register every tests/unit_tests/<case>/test_*.py with CTest.
#
# The scripts locate the repository and their own input_data/ and ref_data/ relative to
# __file__, so no WORKING_DIRECTORY is needed and CTest can run them from the build tree.
# Each script exits 0 on success and 1 on failure, which is what CTest reads.
#
# Note the interpreter matters: it must be the one with numpy that the build targets. If
# CMake picks the wrong one, configure with
#     -DPython3_EXECUTABLE=/path/to/python3

find_package(Python3 REQUIRED COMPONENTS Interpreter)

file(GLOB UNIT_TEST_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test_*.py")

foreach(script IN LISTS UNIT_TEST_SCRIPTS)
    file(RELATIVE_PATH rel "${CMAKE_CURRENT_SOURCE_DIR}" "${script}")
    string(REGEX REPLACE "\\.py$" "" name "${rel}")
    string(REPLACE "/" "." name "${name}")
    add_test(NAME unit.${name} COMMAND ${Python3_EXECUTABLE} "${script}")
    set_tests_properties(unit.${name} PROPERTIES LABELS unit)
endforeach()

list(LENGTH UNIT_TEST_SCRIPTS n_unit_tests)
message(STATUS "unit_tests: registered ${n_unit_tests} CTest case(s)")
