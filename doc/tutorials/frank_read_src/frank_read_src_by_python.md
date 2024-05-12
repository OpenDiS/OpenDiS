### Frank-Read Source by Pure Python

Frank-Read Source is a mechanism of the generation of dislocations during the plastic deformation of solids. It explains the increase of dislocation density during the material's deformation process. This simple example demonstrates how to simulate Frank-Read Source in pure python mode (PyDiS) using OpenDiS.

To run the simulation, simply execute:

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src
python3 -i test_frank_read_src.py
```

In the simulation, there are two dislocation jogs pinning the four dislocation nodes. The bottom two nodes are fixed so that the dislocation line glides on the top plane. The propagation of the dislocation subsumes force calculation, collision, mobility law, remeshing rules, time integration and topolgical operation, as defined in ```pydis``` ([Bulatov & Cai](https://core.ac.uk/reader/44178170)).

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
