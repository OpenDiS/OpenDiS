### Frank-Read Source by Pure Python
Run OpenDiS in pure python mode (PyDiS) -- no compilation needed

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src
python3 -i test_frank_read_src.py
```

```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.
2. how to examine the data structure of the dislocation network (DisNet) in python interactive
```

**Initial condition:** 

The initial condition is described in the function ```init_frank_read_src_loop()```. The initial nodes are described in the ```cell_list```, where the ```CellList``` function is being imported from the DisNet manager to initialize the simulation cell. One defines the simulation cell (```cell = ```). The initial nodal positions are introduced in the variable ```rn```.


**Boundary condition:**

In this example, the periodic boundary condition (PBC) is being implemented. One can turn the ```pbc``` option on or off by setting ```init_frank_read_src_loop(pbc=True)```.

**Explore dislocation data:**

```python
net.G.nodes
```
```python
net.G.edges
```

```python
net.export_data()
```

```{figure} frank_read_schematic.png
:alt: Screenshot of the final configuration
:width: 552px
```

```{figure} frank_read_vid.gif
:alt: Video of the whole simulation
:width: 552px
```



```{hint}
In case the window for the visualization is not showing up, please do ```ssh -Y <your_account>@your.address```
```
