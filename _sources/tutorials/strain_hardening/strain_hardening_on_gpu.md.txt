### Strain Hardening Simulation on GPU
We can run a tensile test simulation of a single-crystal Cu using the following commands (running ExaDiS on GPU).  GPU allows us to run the simulation more efficiently and hence reach a greater strain to see the strain-hardening behavior more clearly.
Execute the following commands to enter the directory containing the example.
```bash
cd ${OPENDIS_DIR}
cd examples/10_strain_hardening/
```
Before running the simulation, we shall edit the ```test_strain_hardening_exadis.py``` file to change the constructor of ```sim``` from ```maxstep=100``` to ```max_step=10000```.
Execute the following commands to run the simulation.
```bash
export OMP_NUM_THREADS=8
python3 test_strain_hardening_exadis.py
```

#### Simulation Behavior
The simulation creates a folder called ```output_fcc_Cu_15um_1e3``` to store the results files.  On MC3.stanford.edu (gpu-ampere), it takes about 13.6 hours to run 10000 steps of the simulation.  The simulation will write a data file to the output folder for every 100 steps.  The ```stress_strain_dens.dat``` file stores certain essential information of the tensile test --- it contains 5 columns corresponding to step, strain, stress (Pa), dislocation density (m<sup>-2</sup>) and wall-clock time (sec), respectively.

The final dislocation configuration (config.10000.data) after 10000 steps is shown below.
```{figure} GPU_final_configuration_Ovito.png
:alt: Screenshot of the final configuration
:width: 552px
```

The predicted stress-strain curve is shown below.
```{figure} Stress_strain_ampere.png
:alt: stress-strain curve
:width: 352px
```

Here is how the total dislocation density changes with strain.  The increase of dislocation density (i.e. dislocation multiplication) with strain is a key mechanism for strain-hardening.
```{figure} Density_strain_ampere.png
:alt: dislocation density-strain curve
:width: 352px
```
