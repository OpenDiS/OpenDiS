cmake_minimum_required(VERSION 3.16.8)

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

if(Kokkos_ENABLE_CUDA STREQUAL "On")
    set(CMAKE_CXX_COMPILER ${CMAKE_CURRENT_SOURCE_DIR}/core/exadis/kokkos/bin/nvcc_wrapper)
endif()

message(" CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}")
message("***************************************")


project(OpenDiS LANGUAGES CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX install)

add_subdirectory(core/exadis)
#add_subdirectory(core/exadis/kokkos)


#enable_testing()
#add_subdirectory(ctests)
