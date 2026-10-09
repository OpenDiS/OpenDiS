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
