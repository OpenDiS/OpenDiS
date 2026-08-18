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

find_package(Python3 REQUIRED COMPONENTS Interpreter)

file(GLOB FULL_RUN_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test_*.py")

foreach(script IN LISTS FULL_RUN_SCRIPTS)
    file(RELATIVE_PATH rel "${CMAKE_CURRENT_SOURCE_DIR}" "${script}")
    string(REGEX REPLACE "\\.py$" "" name "${rel}")
    string(REPLACE "/" "." name "${name}")
    add_test(NAME full_run.${name} COMMAND ${Python3_EXECUTABLE} "${script}")
    # each case writes its output into the directory it runs in, kept inside the build tree
    set_tests_properties(full_run.${name} PROPERTIES
        LABELS full_run
        TIMEOUT 1800)
endforeach()

list(LENGTH FULL_RUN_SCRIPTS n_full_runs)
message(STATUS "full_runs: registered ${n_full_runs} CTest case(s)")
