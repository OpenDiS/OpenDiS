cmake_minimum_required(VERSION 3.16)

option(Kokkos_ENABLE_SERIAL "kokkos enable serial" On)
option(Kokkos_ENABLE_OPENMP "kokkos enable openmp" On)
option(Kokkos_ENABLE_CUDA "kokkos enable cuda" Off)
option(Kokkos_ENABLE_CUDA_LAMBDA "kokkos enable cuda lambda" Off)
option(Kokkos_ENABLE_HIP "kokkos enable hip" Off)

message("***************************************")
if(DEFINED SYS)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/sys.cmake.${SYS}")
    message(" options set by cmake/sys.cmake.${SYS}")
    include(cmake/sys.cmake.${SYS})
  else()
    message(" cmake/sys.cmake.${SYS} does not exist")
  endif()
else()
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/sys.cmake.ext")
    message("options set by cmake/sys.cmake.ext")
    # the options can be modified by the following file if it exists
    include(cmake/sys.cmake.ext OPTIONAL)
  else()
    message(" cmake/sys.cmake.ext does not exist")
  endif()
endif()

if(Kokkos_ENABLE_CUDA AND NOT DEFINED CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER ${CMAKE_CURRENT_SOURCE_DIR}/core/exadis/kokkos/bin/nvcc_wrapper)
endif()

if(DEFINED CMAKE_CXX_COMPILER)
    message(" CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
endif()
message("***************************************")


project(OpenDiS LANGUAGES CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX install)

# A _repro build must agree bitwise with every other _repro build, so it must
# round every floating-point operation the same way. Compilers may instead fuse
# a*b + c into one multiply-add, rounded once rather than twice, wherever the
# target has the instruction: clang does by default (-ffp-contract=on), so on
# arm64, which always has it, every compiled kernel moved in its last bits
# against x86_64, which has it only with -mfma or a -march that includes it.
# So contraction is off for every C and C++ target of such a build, the pydis
# library, ExaDiS and Kokkos alike, all added below. On x86_64 without FMA this
# changes nothing; on x86_64 with -mfma it gives back those same results.
if(PYDIS_BITREPRO_MATH)
    add_compile_options($<$<COMPILE_LANGUAGE:C,CXX>:-ffp-contract=off>)
    message(" floating-point contraction off (-ffp-contract=off)")
endif()

add_subdirectory(core/pydis)
add_subdirectory(core/exadis)

# Bitwise-reproducible log/atan for ExaDiS' segment-segment kernel, on the same
# flag the pydis half uses, so a _repro build turns both on together. See
# cmake/exadis_bitrepro/force_common.h for how the redirection is reached
# without touching the submodule, and why it is not done with -include.
#
# portable_math.c is compiled into the ExaDiS targets (through a one-line .cpp
# wrapper, so no C language is needed) rather than duplicated:
# one implementation, so the two codes cannot drift to different roundings of
# the same function. That does couple the two trees at the file level, which is
# acceptable for a local build and would not be for an upstream contribution.
if(PYDIS_BITREPRO_MATH)
    set(EXADIS_BITREPRO_INC ${CMAKE_CURRENT_SOURCE_DIR}/cmake/exadis_bitrepro)
    set(EXADIS_BITREPRO_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/exadis_bitrepro/portable_math_cxx.cpp)
    foreach(_t pyexadis exadis)
        if(TARGET ${_t})
            target_include_directories(${_t} BEFORE PRIVATE ${EXADIS_BITREPRO_INC})
            target_sources(${_t} PRIVATE ${EXADIS_BITREPRO_SRC})
            message(" ${_t}: bitwise-reproducible log/atan ON")
        endif()
    endforeach()
endif()

if(DEFINED cpython_lib_fix_cmd)
    add_custom_target(cpython_lib_fix DEPENDS pyexadis)
    add_custom_command(TARGET cpython_lib_fix POST_BUILD COMMAND ${cpython_lib_fix_cmd})
endif()

enable_testing()
add_subdirectory(tests)
