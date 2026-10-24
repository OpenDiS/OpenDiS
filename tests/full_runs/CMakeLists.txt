# Register every tests/full_runs/<case>/test_*.py with CTest.
#
# These are wrappers: they run the scripts under examples/ as subprocesses and compare the
# results, so the examples remain the single source of truth for how a simulation is set up.
# They are far slower than the unit tests, hence the separate "full_run" label:
#     ctest -L unit        fast
#     ctest -L full_run    slow
#
# The test list is a glob, matching the $(wildcard test_*.py) used by the makefiles, so make
# and ctest cover the same set.

# One working directory per case, inside the build tree. A test writes its output into the
# directory it runs in, which is what you want running it by hand from its own folder, and
# the reason they must not share one under ctest: 02's clear_recordings sweeps
# output/ordered_pydis, the names 03 writes to, so run together in one directory 02 deletes
# 03's recordings mid-run. Keyed on the test name, so two cases of one script separate.
function(add_full_run_case name script)
    add_test(NAME ${name} COMMAND ${OPENDIS_TEST_PYTHON} "${script}" ${ARGN})
    set(rundir "${CMAKE_CURRENT_BINARY_DIR}/run/${name}")
    file(MAKE_DIRECTORY "${rundir}")   # ctest does not create WORKING_DIRECTORY
    set_tests_properties(${name} PROPERTIES
        LABELS full_run
        TIMEOUT 1800
        WORKING_DIRECTORY "${rundir}")
endfunction()

file(GLOB FULL_RUN_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test_*.py")

# Scripts that compare one coarsening branch per invocation. They are registered once per
# mode below, with --coarsen-mode always given, and left out of the plain registration
# here. Passing it in both cases rather than letting one of them fall through to the
# script's default keeps ctest from depending on what that default is, and keeps mode 0
# from looking like the special one.
set(COARSEN_SCRIPTS
    "${CMAKE_CURRENT_SOURCE_DIR}/02_frank_read_src/test_frank_read_src_pydis_exadis_elast.py")
set(COARSEN_MODES 0 1)

foreach(script IN LISTS COARSEN_SCRIPTS)
    if(NOT EXISTS "${script}")
        message(FATAL_ERROR "full_runs: ${script} is listed in COARSEN_SCRIPTS but does "
                            "not exist. Renaming it silently drops its coverage.")
    endif()
endforeach()

set(n_full_runs 0)

foreach(script IN LISTS FULL_RUN_SCRIPTS)
    file(RELATIVE_PATH rel "${CMAKE_CURRENT_SOURCE_DIR}" "${script}")
    string(REGEX REPLACE "\\.py$" "" name "${rel}")
    string(REPLACE "/" "." name "${name}")
    if(script IN_LIST COARSEN_SCRIPTS)
        foreach(mode IN LISTS COARSEN_MODES)
            add_full_run_case(full_run.${name}.coarsen_mode${mode} "${script}"
                              --coarsen-mode ${mode})
            math(EXPR n_full_runs "${n_full_runs} + 1")
        endforeach()
    else()
        add_full_run_case(full_run.${name} "${script}")
        math(EXPR n_full_runs "${n_full_runs} + 1")
    endif()
endforeach()

message(STATUS "full_runs: registered ${n_full_runs} CTest case(s)")
