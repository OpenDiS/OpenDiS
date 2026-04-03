### Compile on MC3

[Mc3](https://hpcc-intranet.stanford.edu/resources/mc3-cluster/).stanford.edu is a Linux ([Rocky Linux](https://rockylinux.org/) 8) cluster.

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
module load cuda/nvhpc/23.7
````

#### Install FFTW
Follow the instructions at this [wiki page](http://micro.stanford.edu/wiki/Install_FFTW3) to install FFTW in your home directory.  The library and include files should be put in paths consistent with what’s specified in the ``cmake/sys.cmake.mc3_cpu`` file.
```{note}
FFTW is only required when compiling for CPU/OMP. When compiling for cuda GPU, FFT transforms will be performed using the `cuFFT` library.  
```
  
#### Build ExaDiS/KOKKOS (OMP version)

````bash
cd ${OPENDIS_DIR}
conda activate opendis
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
If you encounter the error ```pyexadis not found```, it is likely that the Python used for compiling and the Python used for running are different.  The Python executable used in compiling is specified in the ```cmake/sys.cmake.mc3_cpu``` file.  Check which Python version you are running by e.g. ```python3 --version```. Make sure you are in the ``opendis`` conda environment using the following command before both compiling and running.
```bash
conda activate opendis
```
When compilation is successful, you should see a file like ``pyexadis.cpython*.so`` in the ``core/exadis/python`` folder.

#### Run test case (OMP version)

````bash
export OMP_NUM_THREADS=8
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````

#### Build ExaDiS/KOKKOS (GPU version)

````bash
# first we need to get on a GPU node (gpu-ampere) in order to use nvcc
srun -p gpu-ampere -n 1 --x11 --pty bash
cd ${OPENDIS_DIR}
conda activate opendis
rm -rf build/; ./configure.sh -DSYS=mc3_ampere
cmake --build build -j 8 ; cmake --build build --target install
````

```{hint}
The Python executable used in compiling is ```${HOME}/.conda/envs/opendis/bin/python```, as specified in the ```cmake/sys.cmake.mc3_tesla``` file.  Make sure you use the same Python version when running as that used for compiling by using the following command.
```bash
conda activate opendis
```

#### Run test case (GPU version)

````bash
export OMP_NUM_THREADS=8
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
conda activate opendis
python3 -i test_frank_read_src_exadis.py
````
