## Obtaining, building and running the code

ExaDiS is implemented using the [Kokkos](https://kokkos.org) framework and built using the CMake build system. The code is available as a standalone, or as part of the OpenDiS framework where the `ExaDiS` repository is included as a submodule.

### Obtaining the code from OpenDiS
```{note}
This is the preferred way of obtaining the code.
```
You can obtain the code from OpenDiS, in which ExaDiS is included as a submodule. To do so, please follow the steps described in the [Installation section](../../../installation/index).


### Obtaining the code as a standalone

Alternatively, you can obtain the code as a standalone from the `ExaDiS` repository. A typical installation of the code follows the steps below:

* Step 1: Clone the repository and submodules
```
git clone --recursive https://github.com/LLNL/exadis.git
cd exadis
```
Alternatively, you can use the following commands to achieve the same
```
git clone https://github.com/LLNL/exadis.git
cd exadis
git submodule init
git submodule update
```

* Step 2: Configure the build for your system by passing build options to the `configure.sh` script. (See list of options in the Build Options section below.)
    * Example: default build with `SERIAL` and `OPENMP` backends
    ```
    ./configure.sh
    ```
    * Example: build with `CUDA` backend and device architecture `VOLTA70`
    ```
    ./configure.sh \
        -DKokkos_ENABLE_CUDA=On \
        -DKokkos_ENABLE_CUDA_LAMBDA=On \
        -DKokkos_ARCH_VOLTA70=On  
    ```
    * You can also use pre-defined build options and/or create your own build options by setting the options in files `cmake/sys.cmake.<mysystem>`, and then passing build argument `-DSYS=<mysystem>`. E.g., to build for `SYS=lassen` (i.e. using options set in file `cmake/sys.cmake.lassen`):
    ```
    ./configure.sh -DSYS=lassen
    ```

* Step 3: Build the code
```
cmake --build build -j8
```

```{note}
Building for GPU (e.g. with `nvcc` or `hipcc`) may be pretty slow, please be patient!
```
For additional building options and troubleshooting see section Detailed build instructions below.

* Step 4: Test your installation by running an example (assuming `-DEXADIS_PYTHON_BINDING=On`)
```
cd examples/02_frank_read_src
python test_frank_read_src.py
```


### Detailed build instructions

#### Dependencies

* Kokkos:
    * ExaDiS is implemented using the Kokkos framework. Kokkos is included as a submodule to the repository and will be automatically cloned to the `kokkos/` folder when using the git submodule commands or cloning with the `--recursive` option (see Step 1 of Quick Start section). By default, Kokkos will be built in-tree while building ExaDiS. ExaDiS will be compiled for the backend(s) selected to build Kokkos. For instance, if Kokkos is built to run on GPUs (e.g. with build option `-DKokkos_ENABLE_CUDA=ON`), then ExaDiS will be compiled to run on GPUs. If a prior Kokkos installation exists on the machine, its installation path can be provided with ExaDiS build option `-DKokkos_ROOT`, in which case Kokkos will not be built in-tree. Instructions on how to configure/install Kokkos are found at https://github.com/kokkos/kokkos.
    
* FFT libraries
    * ExaDiS uses FFT libraries to compute long-range elastic interactions. To compile ExaDiS without this module (e.g. if no FFT library is available) use build option `-DEXADIS_FFT=Off`. Otherwise (default), different FFT libraries are invoked depending on the target backend:
        * Serial/OpenMP backend: uses FFTW. Include and library directories can be specified with build options `FFTW_INC_DIR` and `FFTW_LIB_DIR`, respectively.
        * Cuda backend: uses cuFFT
        * HIP backend: uses hipFFT
        
* pybind11
    * ExaDiS uses [pybind11](https://github.com/pybind/pybind11) for the python binding module. pybind11 is included as a submodule to the repository and will be automatically cloned to the `python/pybind11` folder when using the git submodule commands or cloning with the `--recursive` option (see Step 1 of Quick Start section).
    To use a specific python version/executable, use build option `PYTHON_EXECUTABLE`. If needed, the include path to the `python-dev` package (containing file `Python.h`) can be provided with build option `PYTHON_DEV_INC_DIR`.
    To compile ExaDiS without this module, use build option `-DEXADIS_PYTHON_BINDING=Off`.


#### Build options

Below is a list of the various CMake build option specific to ExaDiS. The build options are passed as arguments to the cmake command as `-D<BUILD_OPTION_NAME>=<value>`.

* `EXADIS_PYTHON_BINDING` (optional, default=`On`): enable/disable compilation of the python module
* `PYTHON_EXECUTABLE` (optional, default=''): specifies the path of a specific python version to be used
* `PYTHON_DEV_INC_DIR` (optional, default=''): specifies the path to the python-dev include directory
* `EXADIS_FFT` (optional, default=`On`): enable/disable compilation of the FFT-based long-range force calculation module
* `FFTW_INC_DIR` (optional, default=''): specifies the path of the FFTW include directory
* `FFTW_LIB_DIR` (optional, default=''): specifies the path of the FFTW library directory
* `EXADIS_BUILD_EXAMPLES` (optional, default=`Off`): builds examples that are in the `examples/` folder
* `EXADIS_BUILD_TESTS` (optional, default=`Off`): builds test cases that are in the `tests/` folder

Kokkos related main build options: (see full list [here](https://kokkos.org/kokkos-core-wiki/keywords.html))
* `Kokkos_ENABLE_SERIAL` (optional, default=`On`): enable/disable compilation with the serial (CPU) backend
* `Kokkos_ENABLE_OPENMP` (optional, default=`On`): enable/disable compilation with the OpenMP backend
* `Kokkos_ENABLE_CUDA` (optional, default=`Off`): enable/disable compilation with the CUDA backend. If `On`, option `-DKokkos_ENABLE_CUDA_LAMBDA=On` is also required, and a device architecture must be provided, e.g. `-DKokkos_ARCH_VOLTA70=On`.
* `Kokkos_ENABLE_HIP` (optional, default=`Off`): enable/disable compilation with the HIP backend.
* `Kokkos_ROOT` (optional, default=none) : specifies the path to a pre-existing Kokkos installation. Do not specify any of the above Kokkos options if this option is used; ExaDiS will be built with the backends that the pre-existing Kokkos installation was built for.
