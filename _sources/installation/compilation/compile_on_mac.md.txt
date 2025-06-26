### Compile on Mac

#### Install required packages
If CMake is having problem finding the [FFTW](https://www.fftw.org/) package on your system, you can install them manually and specify its location in the ``cmake/sys.cmake.ext`` file.  For example, assuming that you have installed FFTW in your home directory, you may add the following lines in your ``cmake/sys.cmake.ext`` file (and then configure without -DSYS , see below).
```cmake
set(FFTW_LIB_DIR $ENV{HOME}/usr/lib)
set(FFTW_INC_DIR $ENV{HOME}/usr/include)

message("FFTW_LIB_DIR = ${FFTW_LIB_DIR}")
message("FFTW_INC_DIR = ${FFTW_INC_DIR}")
```


#### Build ExaDiS/KOKKOS
```bash
cd ${OPENDIS_DIR}
rm -rf build/; ./configure.sh -DSYS=mac
cmake --build build -j 8 ; cmake --build build --target install
```

Alternatively, you can also copy ``cmake/sys.cmake.ubuntu`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS .  The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

```bash
cd ${OPENDIS_DIR}
cp cmake/sys.cmake.ubuntu cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
```

When compilation is successful, you should see a file like ``pyexadis.cpython*.so``  in the ``core/exadis/python`` folder.  If you encounter errors, the [following section](#compiling-problems) may help you.

#### Run test case (OMP version)

```bash
export OMP_NUM_THREADS=8
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
python3 -i test_frank_read_src_exadis.py
```


#### Compiling Problems

If your Mac does not have OpenMP enabled, you can try to install it using brew.
```bash
brew list
brew update
xcode-select --install
```

Here is an alternative way to install OpenMP on Mac.
```bash
brew install llvm
brew install libomp
```

If you cannot install OpenMP on your Mac, you can turn it off during compilation (configure) as follows.
```bash
cd ${OPENDIS_DIR}
rm -rf build/; ./configure.sh -DSYS=mac -DKokkos_ENABLE_OPENMP=off
cmake --build build -j 8 ; cmake --build build --target install
```

