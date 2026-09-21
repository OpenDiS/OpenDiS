# Register every tests/unit_tests/<case>/test_*.py with CTest.
#
# The scripts locate the repository and their own input_data/ and ref_data/ relative to
# __file__, so those follow the script rather than the working directory. Each script exits
# 0 on success and 1 on failure, which is what CTest reads.
#
# The interpreter is OPENDIS_TEST_PYTHON, resolved in tests/CMakeLists.txt.

# One working directory per case, inside the build tree, as tests/full_runs does and for
# the same reason: a test writes its output into the directory it runs in, which is what
# you want running it by hand from its own folder, and means they must not share one under
# ctest. Keyed on the test name, so two cases of one script separate.
function(add_unit_case name script)
    add_test(NAME ${name} COMMAND ${OPENDIS_TEST_PYTHON} "${script}" ${ARGN})
    set(rundir "${CMAKE_CURRENT_BINARY_DIR}/run/${name}")
    file(MAKE_DIRECTORY "${rundir}")   # ctest does not create WORKING_DIRECTORY
    set_tests_properties(${name} PROPERTIES LABELS unit WORKING_DIRECTORY "${rundir}")
endfunction()

file(GLOB UNIT_TEST_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test_*.py")

# Scripts that compare one coarsening algorithm per invocation. They are registered once
# per mode below, with --coarsen-mode always given, and left out of the plain registration
# here. Passing it in both cases rather than letting one of them fall through to the
# script's default keeps ctest from depending on what that default is, and keeps mode 0
# from looking like the special one. test3_remesh_rule/makefile runs both modes the same
# way, so make and ctest cover the same set.
set(COARSEN_SCRIPTS
    "${CMAKE_CURRENT_SOURCE_DIR}/test3_remesh_rule/test_remesh_lengthbased_pydis_exadis.py")
set(COARSEN_MODES 0 1)

foreach(script IN LISTS COARSEN_SCRIPTS)
    if(NOT EXISTS "${script}")
        message(FATAL_ERROR "unit_tests: ${script} is listed in COARSEN_SCRIPTS but does "
                            "not exist. Renaming it silently drops its coverage.")
    endif()
endforeach()

set(n_unit_tests 0)

foreach(script IN LISTS UNIT_TEST_SCRIPTS)
    file(RELATIVE_PATH rel "${CMAKE_CURRENT_SOURCE_DIR}" "${script}")
    string(REGEX REPLACE "\\.py$" "" name "${rel}")
    string(REPLACE "/" "." name "${name}")
    if(script IN_LIST COARSEN_SCRIPTS)
        foreach(mode IN LISTS COARSEN_MODES)
            add_unit_case(unit.${name}.coarsen_mode${mode} "${script}"
                          --coarsen-mode ${mode})
            math(EXPR n_unit_tests "${n_unit_tests} + 1")
        endforeach()
    else()
        add_unit_case(unit.${name} "${script}")
        math(EXPR n_unit_tests "${n_unit_tests} + 1")
    endif()
endforeach()

message(STATUS "unit_tests: registered ${n_unit_tests} CTest case(s)")
