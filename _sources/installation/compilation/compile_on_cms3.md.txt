### Compile on CMS3

Cms3-fast.tamu.edu is a high performance computing (HPC) cluster at Texas A&M University. It has NVIDIA T4 GPUs. OpenDiS has been used in the teaching at the summer school of [CMS3](https://cms3.tamu.edu/) in July 2024.

#### Load modules 

Put the following lines in your ``~/.bash_profile file``, exit and login again
````bash
module load gcc/8.3.0
module load cuda/12.2.2

````
  
#### Build ExaDiS/KOKKOS (GPU version)

````bash
cd ${OPENDIS_DIR}
rm -rf build/; ./configure.sh -DSYS=cms3-fast
cmake --build build -j 8 ; cmake --build build --target install
````

Alternatively, you can also copy the ``cmake/sys.cmake.cms3-fast`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

````bash
cp cmake/sys.cmake.cms3-fast cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
````

When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.


#### Run test case (GPU version)

````bash
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````
