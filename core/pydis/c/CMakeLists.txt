cmake_minimum_required(VERSION 3.14)

add_library(pydis MODULE)

add_subdirectory(calforce)
add_subdirectory(collision)
add_subdirectory(remesh)
add_subdirectory(util)
add_subdirectory(tests)

target_include_directories(pydis PRIVATE include)

# PYDIS_BITREPRO_MATH redirects calforce/SegSegForce.c's log()/atan() calls
# to the portable pydis_log()/pydis_atan() (calforce/portable_math.c) via
# -include rather than by editing the file; see calforce/portable_math_shim.h
# for why. Every SYS whose name ends in _repro (mc3_cpu_repro, mac_repro,
# ...) sets this ON, so that its build agrees bitwise with every other
# _repro build on this kernel regardless of which machine's libm each would
# otherwise have used; every other SYS leaves it unset. This buys bitwise
# reproducibility across machines, not extra accuracy, so it stays opt-in
# rather than the default: pydis_log and pydis_atan agree with the system
# libm to about 1ulp, not exactly.
#
# Must be set here, not in calforce/CMakeLists.txt: source file properties
# are only visible to targets in the same directory as where they were set
# (or its subdirectories), and pydis is created in this, the parent,
# directory.
if(PYDIS_BITREPRO_MATH)
  set_source_files_properties(calforce/SegSegForce.c PROPERTIES
    COMPILE_FLAGS "-include ${CMAKE_CURRENT_SOURCE_DIR}/calforce/portable_math_shim.h")
endif()

install(TARGETS pydis DESTINATION ${CMAKE_SOURCE_DIR}/lib)
