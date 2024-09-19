## System Requirements

### Running Python examples in OpenDiS

You can run some OpenDiS examples that do not require any compilation, as long as you have Python3 installed. Refer to [Python tutorial](../tutorials/frank_read_src/frank_read_src_by_python.md) for more details.

### Requirements for high performance
 
To run high performance simulations (CPU and GPU), OpenDiS relies on the core [ExaDiS](https://github.com/LLNL/exadis) library included as a submodule to the project. ExaDiS is a C++ code built on [Kokkos](https://kokkos.org) that needs to be compiled before it can be used. Before you begin the installation of ExaDiS/Kokkos, make sure the following software packages are installed on your system:

- **CMake**
  - Version: 3.16 or higher
  - [Download CMake](https://cmake.org/download/)

- **GCC**
  - Version: 8.2 or higher
  - GCC can be installed via your package manager on Linux systems (e.g., `sudo apt install gcc` for Ubuntu).

- **CUDA**
  - Version: 11.5 or higher (optional, for NVIDIA GPU support)
  - [Download CUDA](https://developer.nvidia.com/cuda-downloads)
  
- **ROCM**
  - Version: 5.2 or higher (optional, for AMD GPU support)
  - [ROCM documentation](https://rocm.docs.amd.com)

- **Python Development Package**
  - Version: 3.6 or higher
  - Ensure Python and the development headers are installed (e.g., `sudo apt install python3-dev` for Ubuntu).
  - Some of the Python packages required include [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/).
  - Optional Python package [networkx](https://networkx.org/), needed only if you call the DisNet to_networkx(), from_networkx() functions

- **FFTW**
  - Install [FFTW](https://www.fftw.org/) for handling discrete Fourier transforms (only required for CPU build)
  - It can typically be installed via your package manager (e.g., `sudo apt install libfftw3-dev`).
