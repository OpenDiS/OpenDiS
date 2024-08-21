### Compile on Android Phone (Termux)

[Termux](https://termux.dev/en/) is a terminal emulator and Linux environment app that runs on Android phones and tablets.

#### Install required packages
```bash
pkg install git
pkg install libx11 clang
pkg install xorgproto
pkg install binutils
pkg install cmake
pkg install vim
pkg install python3
pkg install python-numpy
python3 -m pip install pillow
python3 -m pip install networkx
pkg install matplotlib
pkg install pyqt5
pkg install tigervnc
# run the following command and then launch RVNC app if you want to see graphics
vncserver -localhost
```

#### Build ExaDiS/KOKKOS
```bash
cd ${OPENDIS_DIR}
rm -rf build/; ./configure.sh -DSYS=termux
cmake --build build -j 8 ; cmake --build build --target cpython_lib_fix
```

Alternatively, you can also copy ``cmake/sys.cmake.termux`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS .  The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

```bash
cp cmake/sys.cmake.termux cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target cpython_lib_fix
```

When compilation is successful, you should see a file like ``pyexadis.cpython*.so``  in the ``core/exadis/python`` folder. 

#### Run test case (OMP version)

```bash
export OMP_NUM_THREADS=8
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
python3 -i test_frank_read_src_exadis.py
```

```{figure} frank_read_src_on_termux.png
:alt: Screenshot of OpenDiS running on Android phone
:width: 552px
```
