## System Requirements

### Running Python examples in OpenDiS

You can run some OpenDiS examples that do not require any compilation, as long as you have Python3 installed. Refer to [Python tutorial](../tutorials/frank_read_src/frank_read_src_by_python.md) for more details.

### Requirements for high performance
 
Before you begin the installation of ExaDiS/KOKKOS, make sure the following software packages are installed on your system:

- **CMake**
  - Version: 3.16 or higher
  - [Download CMake](https://cmake.org/download/)

- **GCC**
  - Version: 8.3 or higher
  - GCC can be installed via your package manager on Linux systems (e.g., `sudo apt install gcc` for Ubuntu).

- **CUDA**
  - Version: 11.5 or higher (optional, for GPU support)
  - [Download CUDA](https://developer.nvidia.com/cuda-downloads)

- **Python Development Package**
  - Version: 3.6 or higher
  - Ensure Python and the development headers are installed (e.g., `sudo apt install python3-dev` for Ubuntu).
  - Some of the Python packages required include [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/).
  - Optional Python package [networkx](https://networkx.org/), needed only if you call the DisNet to_networkx()/from_networkx() functions

- **FFTW**
  - Install [FFTW](https://www.fftw.org/) for handling discrete Fourier transforms.
  - It can typically be installed via your package manager (e.g., `sudo apt install libfftw3-dev`).
