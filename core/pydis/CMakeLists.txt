cmake_minimum_required(VERSION 3.14)

project(PyDiS LANGUAGES C)
set(LIB_PYDIS_SO libpydis.so)

add_subdirectory(c)

add_subdirectory(python)