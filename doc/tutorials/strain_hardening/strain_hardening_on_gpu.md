### Strain Hardening Simulation on GPU
To run strain hardening simulation on GPU, execute the same commands as those written in [Sec. 3.3.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/strain_hardening/strain_hardening_on_cpu.html) in the same example directory. The outputs are also the same as those calculated on CPU.

To run the simulation, execute the following commands in the corresponding example directory:
```bash
cd ~/Codes/OpenDiS.git/examples/10_strain_hardening/
export OMP_NUM_THREADS=8
python3 test_strain_hardening_exadis.py
```
</br>

**Explanation**

The information of all dislocation nodes is stored in the data file “180chains_16.10e.data”, which can be read using
```bash
G = ExaDisNet()
G.read_paradis('180chains_16.10e.data')
```

The simulation settings are assigned via params, calforce, mobility, timeint, collision, topology and remesh. Then, they are used in SimulateNetworkPerf(). The simulation is conducted using 
```bash
sim.run(net)
```
This is a simulation for max_step=100. After simulation, we can see a folder “output_fcc_Cu_15um_1e3”, where there is a data file for every 100 steps and a file called “stress_strain_dens.dat” to store step, strain, stress (Pa), dislocation density (1/m^2) and walltime (sec) in five columns for each step, respectively.

**Initial dislocation configuration**

The initial dislocation configuration is visualized below:
```{figure} initial_con figuration.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

**Stress-strain curve**

The stress-strain curve for a simulation of 10000 steps is shown below:
The initial dislocation configuration is visualized below:
```{figure} Stress_strain.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

**Density-strain curve**

The density-strain curve for a simulation of 10000 steps is shown below:
The initial dislocation configuration is visualized below:
```{figure} Density_strain.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

