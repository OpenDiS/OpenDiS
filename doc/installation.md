# Installation 


# How to get the code

## Include ExaDiS as submodule

To include ExaDiS as a submodule of OpenDiS, use the following steps:

```bash
cd ~/Codes
git clone --recurse-submodule https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
```

Alternatively, you can use the following commands to achieve the same:

```bash
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
git submodule update --init --recursive
```

## Getting them separately

If you prefer to get the submodules separately, follow these steps:

```bash
git clone https://gitlab.com/micronano/OpenDiS.git OpenDiS.git
cd OpenDiS.git
cd core
git clone https://github.com/LLNL/exadis exadis
cd exadis
git clone https://github.com/kokkos/kokkos.git --branch 4.2.00
cd kokkos
```

# System Requirements

## Running Python examples in OpenDiS

You can run some OpenDiS examples that do not require any compilation, as long as you have Python3 installed. Refer to *Section 3.1* for more details.

## Requirements for high performance

For high-performance applications, you will need to compile ExaDiS and KOKKOS. Ensure the following packages are installed on your system before beginning to build ExaDiS/KOKKOS:

- **cmake**: version ?? or higher
- **gcc**: version ?? or higher
- **cuda**: version ?? or higher (if you want to run on GPU)
- **Python development package**
- **FFTW**

# Compile and Install

## Compile on Mac

### Build ExaDiS/KOKKOS

```bash
cd ~/Codes/OpenDiS.git/
rm -rf build/; ./configure.sh -DSYS=mac
cmake --build build -j 8 ; cmake --build build --target install
```

Alternatively, you can also copy the ``cmake/sys.cmake.ubuntu`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

```bash
cp cmake/sys.cmake.ubuntu cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
```

When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.

# Run test case (OMP version)


``` code-block:: bash

   export OMP_NUM_THREADS=8
   cd examples/02_frank_read_src
   python3 -i test_frank_read_src_exadis.py
```

# Compile on Sherlock.stanford.edu (Linux cluster - CentOS)


Put the following lines in your ``~/.bash_profile`` file, exit and login again:

``` bash

   module load cmake/3.24.2
   module load cuda/11.0.3
   module load gcc/9.1.0
   module load fftw/3.3.10 
   module load py-numpy/1.24.2_py39
```
