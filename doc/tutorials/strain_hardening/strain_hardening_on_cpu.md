### Strain Hardening Simulation on CPU

**Introduction**

Strain hardening, also known as work hardening, in single crystals is a process that increases the strength and hardness of a material by plastic deformation. This phenomenon occurs due to the movement and interaction of dislocations within the crystal structure. This is an example demonstrating how to simulate strain hardening with ExaDiS on CPU.

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
sim.run(net, state)
```
This is a simulation for max_step=100. After simulation, we can see a folder “output_fcc_Cu_15um_1e3”, where there is a data file for every 100 steps and a file called “stress_strain_dens.dat” to store step, strain, stress (Pa), dislocation density (1/m^{2}) and walltime (sec) in five columns for each step, respectively.
</br>

**Initial dislocation configuration**

The dimensions for initial dislocation configuration (180chains_16.10e.data) are ∼ 15 𝜇m × 15 𝜇m × 15 𝜇m. Periodic boundary condition is applied along all three dimensions. The dislocation density is 𝜌0 ≈ 1.2 × 10^{12} m^{−2}. It is visualized below using Ovito and the Ovito settings are also shown below.
```{figure} initial_configuration_Ovito.png
:alt: Screenshot of the final configuration
:width: 552px
```

```{figure} Ovito_settings.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

**Final dislocation configuration**

The final dislocation configuration (config.1600.data) in this simulation is visualized below using Ovito. The Ovito settings are the same as the above screenshot.
```{figure} CPU_final_configuration_Ovito.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

**Stress-strain curve**

A simulation of 10000 steps is conducted as an example, which costs 2 days.
The stress-strain curve is shown below:
```{figure} Stress_strain_GPU.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

**Density-strain curve**

A simulation of 10000 steps is conducted as an example, which costs 2 days.
The density-strain curve is shown below:
```{figure} Density_strain_GPU.png
:alt: Screenshot of the final configuration
:width: 552px
```
</br>

The below python script is used to plot the above curves:
```bash
import numpy as np
import matplotlib.pyplot as plt

step, strain, stress, density, walltime = np.loadtxt('stress_strain_dens.dat', usecols=(0,1,2,3,4), unpack=True)

plt.figure(figsize=(8, 8))
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24
plt.rcParams['mathtext.fontset'] = 'stix'
plt.plot(strain*100, stress/1000000, linewidth=1.5, color='b')
plt.xlabel('Strain (%)')
plt.ylabel('Stress (MPa)')

plt.figure(figsize=(8, 8))
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24
plt.rcParams['mathtext.fontset'] = 'stix'
plt.plot(strain*100, density, linewidth=1.5, color='b', label='OpenDiS')
plt.xlabel('Strain (%)')
plt.ylabel('Density (1/m$^{2}$)')

print("Figure has been done!")

plt.show()
```
