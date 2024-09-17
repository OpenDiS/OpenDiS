# Compilation Overview

OpenDiS is a python-based framework that relies on core libraries to implement calculation modules, see [Project Overview](../../code_structure/project_overview.md). While some OpenDiS scripts can be run without requiring compilation of the code, most of the simulations would require compilation to enjoy high-performance modules.


## CMake build

OpenDiS is compiled using CMake. Core libraries and modules that need to be compiled are included in the OpenDiS compilation tree. At the moment, a build includes compilation of the ExaDiS core library, and some modules of the PyDiS library.

A typical installation follows the steps below:

### Step 1: Prepare the build directory

Make sure any previous build directory is removed:
```bash
cd ${OPENDIS_DIR}
rm -rf build/
```

### Step 2: Configure the build

Configure the build for your system by passing build options to the `configure.sh` script. (See list of options in the Build Options section below.)

Below are various examples of how OpenDiS can be configured:

* Example: use a pre-defined build options defined in `cmake/sys.cmake.<system>`, and then pass build argument `-DSYS=<system>`. E.g., to build for `SYS=mac` (i.e. using options set in file `cmake/sys.cmake.mac`):
```bash
./configure.sh -DSYS=mac
```
For more detailed instructions on pre-defined build options available in OpenDiS, see the next sections of this documentation.

* Example: create your own build options by setting options specific to your system in file e.g. `cmake/sys.cmake.mysystem`, and then pass build argument `-DSYS=mysystem`:
```bash
./configure.sh -DSYS=mysystem
```

* Example: if no `SYS` option is passed to `configure.sh`, e.g.
```bash
./configure.sh
```
then by default the configuration will attempt to read build options from file `cmake/sys.cmake.ext` if it exists. Thus, file `cmake/sys.cmake.ext` can be used as a sandbox to test different build configurations.

* Example: for more advanced users, you can also pass all build options directly to the `configure.sh` script, e.g.
```bash
./configure.sh \
    -DKokkos_ENABLE_CUDA=On \
    -DKokkos_ENABLE_CUDA_LAMBDA=On \
    -DKokkos_ARCH_VOLTA70=On \
    -DEXADIS_BUILD_TESTS=ON \
    -DPYTHON_EXECUTABLE=$(which python)
```

### Step 3: Build and install the codes

Once configuration is completed, build and install the codes using the following command:
```bash
cmake --build build -j 8 ; cmake --build build --target install
```
```{note}
Building for GPU (e.g. with `nvcc` or `hipcc`) may be pretty slow, please be patient!
```


## Build options

Besides the standard `CMAKE` options and the `SYS` option, most OpenDiS build options are related to the ExaDiS/Kokkos core library. See the list of [ExaDiS build options](../../core_libraries/exadis_documentation/user_guide/obtaining.md#detailed-build-instructions).
