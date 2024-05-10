### Compile on MC3

Mc3.stanford.edu is a Linux (CentOS) cluster.

#### Create conda environment (opendis) 
This is for installing Python ``matplotlib`` needed to run test case
```bash
module load intel/python
conda init
conda create --name opendis python=3.7.7
conda activate opendis
conda install matplotlib
````

#### Load modules 

Put the following lines in your ``~/.bash_profile file``, exit and login again
````bash
module load cmake/3.24.1
module load gnu8/8.3.0
module load cuda/11.8
````

#### Install FFTW
Follow the instructions at this [wiki page](http://micro.stanford.edu/wiki/Install_FFTW3) to install FFW in your home directory.  The library and include files should be put in paths consistent with what’s specified in the ``cmake/sys.cmake.mc3_cpu`` file.
  
#### Build ExaDiS/KOKKOS (OMP version)

````bash
cd ~/Codes/OpenDiS.git
rm -rf build/; ./configure.sh -DSYS=mc3_cpu
cmake --build build -j 8 ; cmake --build build --target install
````

Alternatively, you can also copy the ``cmake/sys.cmake.mc3`` file to ``cmake/sys.cmake.ext`` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.

````bash
cp cmake/sys.cmake.mc3_cpu cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
````

When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.

#### Run test case (OMP version)

````bash
export OMP_NUM_THREADS=8
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````

#### Build ExaDiS/KOKKOS (GPU version)

````bash
# first we need to get on a GPU node in order to use nvcc
srun -p gpu-tesla -n 1 --x11 --pty bash
cd ~/Codes/OpenDiS.git
rm -rf build/; ./configure.sh -DSYS=mc3_gpu
cmake --build build -j 8 ; cmake --build build --target install
````

#### Run test case (GPU version)

````bash
export OMP_NUM_THREADS=8
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````
