### Compile on Mac

#### Build ExaDiS/KOKKOS

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

#### Run test case (OMP version)

```bash
export OMP_NUM_THREADS=8
cd examples/02_frank_read_src
python3 -i test_frank_read_src_exadis.py
```
