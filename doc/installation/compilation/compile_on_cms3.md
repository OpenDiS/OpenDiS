### Compile on MC3

Cms3-fast.tamu.edu is a high performance computing (HPC) cluster at Texas A&M University. It has NVIDIA T4 GPUs. OpenDiS has been used in the teaching at the summer school of [CMS3](https://cms3.tamu.edu/) in July 2024.

#### Load modules 

Put the following lines in your ``~/.bash_profile file``, exit and login again
````bash
module load gcc/8.3.0
module load cuda/12.2.2

````

#### Install FFTW
Follow the instructions at this [wiki page](http://micro.stanford.edu/wiki/Install_FFTW3) to install FFW in your home directory.  The library and include files should be put in paths consistent with what’s specified in the ``cmake/sys.cmake.mc3_cpu`` file.
  
#### Build ExaDiS/KOKKOS (GPU version)

````bash
cd ~/Codes/OpenDiS.git
rm -rf build/; ./configure.sh -DSYS=cms3-fast
cmake --build build -j 8 ; cmake --build build --target install
````

Alternatively, you can also copy the ``cmake/sys.cmake.cms3-fast`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

````bash
cp cmake/sys.cmake.mc3_cpu cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
````

When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.


#### Run test case (GPU version)

````bash
export OMP_NUM_THREADS=8
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````
