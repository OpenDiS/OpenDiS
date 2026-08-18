add_subdirectory(unit_tests)
add_subdirectory(full_runs)

# unit_tests are fast and self-contained; full_runs drive the scripts under examples/ and are
# slow. They carry the labels "unit" and "full_run" respectively:
#     ctest -L unit
#     ctest -L full_run
