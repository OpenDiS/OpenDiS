# Pin the OpenMP runtime to the one shipped by the interpreter this build
# targets, when that interpreter ships one.
#
# WHY THIS EXISTS
#
# ExaDiS links an OpenMP runtime, and so does the numpy inside whichever
# python imports pyexadis. If those are two different copies of the LLVM
# OpenMP library, the process has two runtimes and Kokkos::initialize()
# segfaults. The crash is order dependent, so it does not look like a
# library conflict: importing numpy before pyexadis.initialize() dies,
# the reverse order survives, and every pydis-only test passes either way.
#
# The trap is that which interpreter gets targeted is not something the
# SYS file chooses. CMake's FindPython honours CONDA_PREFIX by default
# (Python3_FIND_VIRTUALENV=FIRST), so configuring with a conda
# environment activated builds the module for that environment whatever
# SYS says. -DSYS=mac run inside an activated env therefore produced a
# module built FOR conda python but linked AGAINST the system OpenMP,
# which is exactly the combination that crashes.
#
# Rather than document a rule about which shell to configure from, the
# OpenMP choice follows the interpreter. Then -DSYS=mac is correct in
# either shell.
#
# Only the library is overridden, never OpenMP_CXX_INCLUDE_DIR. Dropping
# the system include path takes fftw3.h with it and the build fails on a
# missing header that has nothing to do with OpenMP. The omp.h that comes
# with the system libomp is ABI compatible with a conda one, both being
# LLVM OpenMP, so leaving the include path alone is correct as well as
# convenient.

# Which interpreter will this build target? An explicit setting wins,
# under any of the three spellings the various finders consult. Otherwise
# ask CMake, whose answer accounts for both PATH and CONDA_PREFIX and is
# therefore the same one pybind11 will arrive at.
set(_omp_target_python "")
foreach(_var Python_EXECUTABLE PYTHON_EXECUTABLE Python3_EXECUTABLE)
  if(NOT _omp_target_python AND DEFINED ${_var})
    set(_omp_target_python "${${_var}}")
  endif()
endforeach()

if(NOT _omp_target_python)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    set(_omp_target_python "${Python3_EXECUTABLE}")
  endif()
endif()

if(_omp_target_python AND EXISTS "${_omp_target_python}")
  get_filename_component(_omp_py_bin "${_omp_target_python}" DIRECTORY)
  get_filename_component(_omp_py_prefix "${_omp_py_bin}" DIRECTORY)
  set(_omp_candidate "${_omp_py_prefix}/lib/libomp.dylib")
  if(EXISTS "${_omp_candidate}")
    set(OpenMP_libomp_LIBRARY "${_omp_candidate}" CACHE FILEPATH "" FORCE)
    message(" OpenMP follows ${_omp_target_python}")
    message("   OpenMP_libomp_LIBRARY = ${OpenMP_libomp_LIBRARY}")
  else()
    # Nothing to match, so the system runtime is the only one in play and
    # is the right choice. Said out loud because the alternative is a
    # segfault whose cause is invisible.
    message(" OpenMP: ${_omp_target_python} ships no libomp.dylib, "
            "using the system runtime")
  endif()
else()
  message(" OpenMP: could not identify the target interpreter, "
          "using the system runtime")
endif()
