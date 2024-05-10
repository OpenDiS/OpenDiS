## System Requirements

### Running Python examples in OpenDiS

You can run some OpenDiS examples that do not require any compilation, as long as you have Python3 installed. Refer to :doc:`../tutorials/frank_read_src/frank_read_src_by_python` for more details.

### Requirements for high performance
 
Before you begin the installation of ExaDiS/KOKKOS, make sure the following software packages are installed on your system:

- **CMake**
  - Version: 3.18 or higher
  - [Download CMake](https://cmake.org/download/)

- **GCC**
  - Version: 9.0 or higher
  - GCC can be installed via your package manager on Linux systems (e.g., `sudo apt install gcc` for Ubuntu).

- **CUDA**
  - Version: 11.0 or higher (optional, for GPU support)
  - [Download CUDA](https://developer.nvidia.com/cuda-downloads)

- **Python Development Package**
  - Ensure Python and the development headers are installed (e.g., `sudo apt install python3-dev` for Ubuntu).

- **FFTW**
  - Install FFTW for handling discrete Fourier transforms.
  - It can typically be installed via your package manager (e.g., `sudo apt install libfftw3-dev`).
