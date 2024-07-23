### Compile on MC3

Mc3.stanford.edu is a Linux (CentOS) cluster.

#### Create conda environment (opendis) 
This is for installing Python3 executable as well as modules (such as ``matplotlib``) needed to run test case
````bash
eval "$(/opt/ohpc/pub/compiler/anaconda3/2024.02-1/bin/conda shell.bash hook)"
conda init
conda create --name opendis python=3.9
conda activate opendis
conda install matplotlib
conda install networkx
````

To ensure ``matplotlib`` runs properly on MC3, we need to specify its backend to be tkagg in the ``matplotlibrc`` file
````bash
vi ~/.config/matplotlib/matplotlibrc
backend: tkagg
````


#### Load modules 

Put the following lines in your ``~/.bash_profile file``, exit and login again
````bash
module load cmake/3.24.2
module load gnu9/9.4.0
module load cuda/12.5
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

```{hint}
Make sure you are in the ``opendis`` conda environment using the following command.
```bash
conda activate opendis
```

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
