cmake_minimum_required(VERSION 3.16.8)

project(OpenDiS LANGUAGES CXX)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX install)

option(Kokkos_ENABLE_SERIAL "kokkos enable serial" ON)
option(Kokkos_ENABLE_OPENMP "kokkos enable openmp" ON)
option(Kokkos_ENABLE_CUDA "kokkos enable cuda" OFF)
option(Kokkos_ENABLE_CUDA_LAMBDA "kokkos enable cuda lambda" OFF)
option(Kokkos_ENABLE_HIP "kokkos enable hip" OFF)


#add_subdirectory(core)
add_subdirectory(core/exadis/kokkos)

#enable_testing()
#add_subdirectory(ctests)
