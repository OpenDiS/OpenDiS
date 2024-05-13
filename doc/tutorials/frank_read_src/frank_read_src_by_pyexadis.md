### Frank-Read Source by Python calling ExaDiS


Run the Frank-Read Source example with ExaDiS, open your OpenDiS, and open the corresponding example directory to execute the following commands:

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src/
export OMP_NUM_THREADS=8
python3 -i test_frank_read_src_exadis.py
```
Similarly, one obtains the simulation visualization snapshots as shown in [Sec. 3.1.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/frank_read_src/frank_read_src_by_python.html). Compared with the serial execution in [3.1.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/frank_read_src/frank_read_src_by_python.html), by setting ```OMP_NUM_THREADS=8```, the simulation speed is significantly increased. In our test on a Stanford supercomputer, the simulation time for running with pure Python is 34.267s, and using the provided OMP threads the simulation time is 11.142s. You can also test it on your own workstation with ```time python3 test_frank_read_src_exadis.py```. By increasing the MPI threads, one expects the simulation speed to increase.

**Explore dislocation data:**

In the examples provided, both ```pydis``` and ```exadis``` utilized the length-based remesh rule. For the pure Python case, the rule is implemented via 
```{python}
remesh    = Remesh(remesh_rule='LengthBased', Lmin=0.1, Lmax=0.3)
```
and in ```exadis``` it's implemented from 
```{python}
remesh    = Remesh(remesh_rule='LengthBased', params=params)
```
Note that in the ```exadis``` case, the parameters can be specified in detail as in line #79:
```{python}
params = {"burgmag": 3e-10, "mu": 160e9, "nu": 0.31, "a": 0.01, "maxseg": 0.3, "minseg": 0.1, "rann": 0.02}
```
Note that this simulation is different from [3.1.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/frank_read_src/frank_read_src_by_python.html) as it turns the PBC option off
```{python}
G = init_frank_read_src_loop(pbc=False)
```
By setting it back in, the two simulations are identical w.r.t. initial and boundary conditions.

```{figure} frank_read_schematic.png
Final simulation snapshot of the Frank-Read Source.
:width: 552px
```

 

```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.

2. can we make the remesh rules consistent between pydis and exadis, no mesh_refine on segments with both ends fixed? Or at least have an option to get to this behavior?

3. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.  How is this simulation different (e.g. time) from 3.1.1?

4. how to examine the data structure of the dislocation network (DisNet) in python interactive


```
