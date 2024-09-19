### Compile on Sherlock

[Sherlock](https://www.sherlock.stanford.edu/).stanford.edu is a Linux (CentOS) cluster.

#### Load modules 

Put the following lines in your ``~/.bash_profile file``, exit and login again
````bash
module load cmake/3.24.2
module load cuda/11.0.3
module load gcc/9.1.0
module load fftw/3.3.10
module load python/3.9.0
module load py-numpy/1.24.2_py39
````

Install Matplotlib python module
````bash
python3 -m pip install matplotlib --user
````
  
#### Build ExaDiS/KOKKOS (OMP version)

````bash
cd ${OPENDIS_DIR}
rm -rf build/; ./configure.sh -DSYS=sherlock
cmake --build build -j 8 ; cmake --build build --target install
````

Alternatively, you can also copy the ``cmake/sys.cmake.sherlock`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

````bash
cp cmake/sys.cmake.sherlock cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
````

When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.

#### Run test case (OMP version)

````bash
export OMP_NUM_THREADS=8
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
python3 -i test_frank_read_src_exadis.py
````
