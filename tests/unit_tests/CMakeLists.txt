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

# Extra invocations, for scripts that cover more than one case through an argument. The glob
# above registers each script once, with no arguments, so anything reached only by a flag
# needs its own entry here or ctest would cover less than the folder makefiles do.
#
# test3_remesh_rule compares one coarsening algorithm per run, chosen by --coarsen-mode, and
# its default is 1. Mode 0 is a different algorithm, not a variation, so it is registered
# separately; test3_remesh_rule/makefile runs both for the same reason.
add_test(NAME unit.test3_remesh_rule.test_remesh_lengthbased_pydis_exadis.coarsen0
         COMMAND ${Python3_EXECUTABLE}
                 "${CMAKE_CURRENT_SOURCE_DIR}/test3_remesh_rule/test_remesh_lengthbased_pydis_exadis.py"
                 --coarsen-mode 0)
set_tests_properties(unit.test3_remesh_rule.test_remesh_lengthbased_pydis_exadis.coarsen0
                     PROPERTIES LABELS unit)
math(EXPR n_unit_tests "${n_unit_tests} + 1")

message(STATUS "unit_tests: registered ${n_unit_tests} CTest case(s)")
