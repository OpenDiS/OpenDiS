# Compilation Overview

OpenDiS is a python-based framework that relies on core libraries to implement calculation modules, see [Project Overview](../../code_structure/project_overview.md). While some OpenDiS scripts can be run without requiring compilation of the code, most of the simulations would require compilation to enjoy high-performance modules.

```{hint}
All scripts/simulations that use modules from the `pyexadis` library require the code to be compiled before they can be run. The `pyexadis` library is built during compilation of the ExaDiS core library, included in the OpenDiS compilation tree.
```


## CMake build

OpenDiS is compiled using CMake. Core libraries and modules that need to be compiled are included in the OpenDiS compilation tree. At the moment, a build includes compilation of the ExaDiS core library (including the `pyexadis` library), and some modules of the PyDiS library.

A typical installation follows the steps below:

### Step 1: Prepare the build directory

Make sure any previous build directory is removed:
```bash
cd ${OPENDIS_DIR}
rm -rf build/
```

### Step 2: Configure the build

Configure the build for your system by passing build options to the `configure.sh` script, including options to build for CPU or GPU architectures. (See list of options in the [Build Options](#build-options) section below.)

Below are various examples of how OpenDiS can be configured:

* Example: default build with ExaDiS `SERIAL` and `OPENMP` backends
```bash
./configure.sh
```

* Example: GPU build with ExaDiS `CUDA` backend and device architecture `AMPERE80` (e.g. A100)
```bash
./configure.sh \
    -DKokkos_ENABLE_CUDA=On \
    -DKokkos_ARCH_AMPERE80=On
```

* Example: use a pre-defined build options defined in `cmake/sys.cmake.<system>`, and then pass build argument `-DSYS=<system>`. E.g., to build for `SYS=mac` (i.e. using options set in file `cmake/sys.cmake.mac`):
```bash
./configure.sh -DSYS=mac
```
For more detailed instructions on pre-defined build options available in OpenDiS, see the next pages of this section.

* Example: create your own build options by setting options specific to your system in file e.g. `cmake/sys.cmake.mysystem`, and then pass build argument `-DSYS=mysystem`:
```bash
./configure.sh -DSYS=mysystem
```

````{tip}
If no `SYS` option is passed to `configure.sh`, e.g. as for a default build
```bash
./configure.sh
```
then by default the configuration will attempt to read build options from file `cmake/sys.cmake.ext` if it exists. Thus, file `cmake/sys.cmake.ext` can be used as a sandbox to test different build configurations.
````
````{Important}
In general, it is recommended to test the installation of ExaDiS by making sure that important unit tests execute correctly and return the right results. To enable the testing, add build option `-DEXADIS_BUILD_TESTS=On`, e.g.
```bash
./configure.sh -DSYS=mysystem -DEXADIS_BUILD_TESTS=On
```
````


### Step 3: Build and install the codes

Once configuration is completed, build and install the codes using the following command:
```bash
cmake --build build -j 8 ; cmake --build build --target install
```
A successful compilation will produce various files in `${OPENDIS_DIR}/lib`, and shared library file `${OPENDIS_DIR}/core/exadis/python/pyexadis.cpython-*.so`.

```{note}
Building for GPU (e.g. with `nvcc` or `hipcc`) may be pretty slow, please be patient!
```

### Step 4: Test your installation by running some examples

Once compilation is completed, run examples from the `examples/` folder, e.g.:

```bash
cd examples/02_frank_read_src
python3 test_frank_read_src_pydis.py
python3 test_frank_read_src_exadis.py
```
More information about these two examples can be found in the [Frank-Read Source](../../tutorials/frank_read_src/index) section, and more examples are detailed in the [Tutorials](../../tutorials/index) section.

````{Note}
When the code was built with option `-DEXADIS_BUILD_TESTS=On` (recommended), ExaDiS unit tests are performed as follows:
```bash
cd core/exadis/tests/unit_tests
python3 run_tests.py
```
All tests should pass.
````
If errors occur when trying to run the examples or unit tests, please refer to the [Troubleshooting](#troubleshooting) section.


## Build options

Besides the standard `CMAKE` options and the `SYS` option, most OpenDiS build options are related to the ExaDiS/Kokkos core library. See the list of [ExaDiS build options](../../core_libraries/exadis_documentation/user_guide/obtaining.md#detailed-build-instructions).


## Troubleshooting

This section reports common errors encountered by OpenDiS users. If you encounter other issues or bugs, please also browse the [Issues](https://github.com/OpenDiS/OpenDiS/issues) section and/or open a new issue.


### pyexadis import errors

Issue: Python is not being able to import `pyexadis` (python binding to the ExaDiS core library), with errors such as

```bash
import pyexadis
ModuleNotFoundError: No module named 'pyexadis'
```
or
```bash
raise ImportError('Cannot import pyexadis')
ImportError: Cannot import pyexadis
```
or
```bash
[Test_Pyexadis_Import] FAIL
[Test_Pyexadis_Init] FAIL
```

This typically happens when the python installation used for compilation is not the same as the python installation used to run the code.

When running the `./configure.sh` script, the python executable to be used for compilation should be indicated, e.g.:
```bash
-- Found PythonInterp: /path/to/python3
```
Then after compilation one should be able to run the scripts using the same python executable, e.g.
```bash
/path/to/python3 test_frank_read_src_exadis.py
```
Alternatively, one can explicitly specify the python executable to be used by passing build option `-DPYTHON_EXECUTABLE=...` to the `configure.sh` script, e.g.
```bash
./configure.sh \
    -DKokkos_ENABLE_CUDA=On \
    -DKokkos_ARCH_AMPERE80=On \
    -DPYTHON_EXECUTABLE=$(which python3)
```
in which case one should be able to run the scripts with, e.g.
```bash
python3 test_frank_read_src_exadis.py
```

### Hanging simulations

Issue: Simulations appear to be hanging after initialization, e.g. after diplaying
```
Run for X steps
```

The cause of this issue is generally an inter-operability issue between the compilers, or between python and the compiler(s).

First, identify whether this issue is specific to `pyexadis` or also occurs for regular ExaDiS simulations. To do so, rebuild the code with option `-DEXADIS_BUILD_TESTS=On`, e.g.
```bash
cd ${OPENDIS_DIR}
rm -rf build/
./configure.sh \
    -DKokkos_ENABLE_CUDA=On \
    -DKokkos_ARCH_AMPERE80=On \
    -DEXADIS_BUILD_TESTS=On
cmake --build build -j 8 ; cmake --build build --target install
```
This should create a few executables in directory `${OPENDIS_DIR}/build/core/exadis/tests`. Run these executables to help identify the source of the issue:
```bash
cd ${OPENDIS_DIR}/build/core/exadis/tests
./test_cuda # if compiled with -DKokkos_ENABLE_CUDA=On
./test_kokkos
./test_system
./test_exadis 0
./test_exadis 1
...
./test_exadis 8
```
* If none of the tests fail but the hanging still happens when running the python-based examples, then the issue is with the compilation of `pyexadis`. First make sure that the python installation used to run is the same as the python installation used for compilation, e.g. see [previous issue](#pyexadis-import-errors). If the issue remains, try using a different combination of compilers and python installations.
* If `test_cuda` fails then the issue is likely with the CUDA compiler
* If `test_kokkos` fails then the issue is with the compilation of Kokkos
* If `test_system` fails then the issue is likely with the use of unified memory
* If `test_exadis` fails then the issue is with the compilation of ExaDiS

Try using a different combination of compilers to resolve compatibility issues. Some known issues are:
* Compatibility issues between gcc11 and some versions of cuda11
* Compatibility issues have been reported between Kokkos and gcc12
* Newest versions of cuda (e.g. 12.5+) have not been fully tested with Kokkos and pybind11 and may trigger compatibility issues
