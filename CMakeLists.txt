cmake_minimum_required(VERSION 3.16.8)

project(OpenDiS LANGUAGES CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX install)

option(Kokkos_ENABLE_SERIAL "kokkos enable serial" On)
option(Kokkos_ENABLE_OPENMP "kokkos enable openmp" On)
option(Kokkos_ENABLE_CUDA "kokkos enable cuda" Off)
option(Kokkos_ENABLE_CUDA_LAMBDA "kokkos enable cuda lambda" Off)
option(Kokkos_ENABLE_HIP "kokkos enable hip" Off)

include(cmake/sys.cmake.ext OPTIONAL)

#add_subdirectory(core/exadis)
add_subdirectory(core/exadis/kokkos)

#enable_testing()
#add_subdirectory(ctests)
