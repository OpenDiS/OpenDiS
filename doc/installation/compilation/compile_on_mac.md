### Compile on Mac

```{Error}
If CMake is having problem finding the FFTW package, you can add the following line to your ``cmake/sys.cmake.ext`` file. In this example, we are assuming that you installed FFTW in your home directory.
```cmake
set(FFTW_LIB_DIR $ENV{HOME}/usr/lib)
set(FFTW_INC_DIR $ENV{HOME}/usr/include)

message("FFTW_LIB_DIR = ${FFTW_LIB_DIR}")
message("FFTW_INC_DIR = ${FFTW_INC_DIR}")
```

#### Build ExaDiS/KOKKOS

```bash
cd ~/Codes/OpenDiS.git/
rm -rf build/; ./configure.sh -DSYS=mac
cmake --build build -j 8 ; cmake --build build --target install
```

```{Error}
the error is on this line ```./configure.sh -DSYS=mac```

```Could NOT find OpenMP_CXX (missing: OpenMP_CXX_FLAGS OpenMP_CXX_LIB_NAMES)```

Need to write ``` rm -rf build/; ./configure.sh -DSYS=mac -DKokkos_ENABLE_OPENMP=off``` instead

additionally make sure to ```conda deactivate```
```


```{Hint}
Alternatively, you can also copy the ``cmake/sys.cmake.mac`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.
```bash
cp cmake/sys.cmake.mac cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
```
```{Error}
CMake Error at /usr/local/Cellar/cmake/3.29.3/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:230 (message):
  Could NOT find OpenMP_CXX (missing: OpenMP_CXX_FLAGS OpenMP_CXX_LIB_NAMES)
```
```{Hint}
If your Mac does not have OpenMP enabled, you can try this.
```bash
brew list
brew update
xcode-select --install
```

```{Hint}
An alternative way to install OpenMP on Mac.
```bash
brew install llvm
brew install libomp
```

When compilation is successful, you should see a file like ```pyexadis.cpython*.so``` in the ```core/exadis/python``` folder.

#### Run test case (OMP version)

```bash
export OMP_NUM_THREADS=8
cd examples/02_frank_read_src
python3 -i test_frank_read_src_exadis.py
```
