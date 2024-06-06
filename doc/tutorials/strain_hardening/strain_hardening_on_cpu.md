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
sim.run(net)
```
This is a simulation for max_step=100. After simulation, we can see a folder “output_fcc_Cu_15um_1e3”, where there is a data file for every 100 steps and a file called “stress_strain_dens.dat” to store step, strain, stress, dislocation density and walltime (sec) in five columns for each step, respectively.
</br>
