### Frank-Read Source by Python calling ExaDiS


Run the Frank-Read Source example with ExaDiS, open your OpenDiS, and open the corresponding example directory to execute the following commands:

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src/
export OMP_NUM_THREADS=8
python3 -i test_frank_read_src_exadis.py
```
Similarly, one obtains the simulation visualization snapshots as shown in [Sec. 3.1.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/frank_read_src/frank_read_src_by_python.html). Compared with the serial execution in [3.1.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/frank_read_src/frank_read_src_by_python.html), by setting ```OMP_NUM_THREADS=8```, the simulation speed is significantly increased. In our test on a Stanford supercomputer, the simulation time for running with pure Python is 34.267s, and using the provided OMP threads the simulation time is 11.142s. You can also test it on your own workstation with ```time python3 test_frank_read_src_exadis.py```. By increasing the MPI threads, one expects the simulation speed to increase.

**Explore dislocation data:**

Similarly, you can also check the ```export_data()``` of the dislocation graph structure with DisNetManager. 

```python
net.G.export_data()
```
You will see something like:
```python

```

```python
net.G.export_data().get("nodes")
```

```python
net.G.export_data().get("segs")
```

 
```{figure} frank_read_schematic.png
frank-read source (fix caption)
```

 

```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.

2. can we make the remesh rules consistent between pydis and exadis, no mesh_refine on segments with both ends fixed? Or at least have an option to get to this behavior?

3. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.  How is this simulation different (e.g. time) from 3.1.1?

4. how to examine the data structure of the dislocation network (DisNet) in python interactive


```
