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

The pydis stores the simulated dislocation networks as graphs. One can explore the properties of this graph by checking its nodes

```python
net.G.nodes
```
and obtain
```python
NodeView(((0, 0), (0, 2), (0, 3), (0, 4), (0, 39), (0, 40), (0, 106), (0, 48), (0, 47), (0, 125), (0, 8), (0, 55), (0, 56), (0, 13), (0, 65), (0, 64), (0, 111), (0, 89), (0, 116), (0, 103), (0, 114), (0, 60), (0, 7)))
```
which represents all the nodes in the dislocation graph.

One can also obtain the graph's edges:
```python
net.G.edges
```
and obtain
```python
OutEdgeView([((0, 0), (0, 4)), ((0, 0), (0, 48)), ((0, 2), (0, 3)), ((0, 2), (0, 47)), ((0, 3), (0, 2)), ((0, 3), (0, 4)), ((0, 4), (0, 3)), ((0, 4), (0, 0)), ((0, 39), (0, 56)), ((0, 39), (0, 13)), ((0, 40), (0, 65)), ((0, 40), (0, 116)), ((0, 106), (0, 114)), ((0, 106), (0, 60)), ((0, 48), (0, 0)), ((0, 48), (0, 13)), ((0, 47), (0, 2)), ((0, 47), (0, 65)), ((0, 125), (0, 89)), ((0, 125), (0, 103)), ((0, 8), (0, 64)), ((0, 8), (0, 116)), ((0, 55), (0, 89)), ((0, 55), (0, 7)), ((0, 56), (0, 39)), ((0, 56), (0, 103)), ((0, 13), (0, 48)), ((0, 13), (0, 39)), ((0, 65), (0, 47)), ((0, 65), (0, 40)), ((0, 64), (0, 8)), ((0, 64), (0, 114)), ((0, 111), (0, 60)), ((0, 111), (0, 7)), ((0, 89), (0, 125)), ((0, 89), (0, 55)), ((0, 116), (0, 40)), ((0, 116), (0, 8)), ((0, 103), (0, 125)), ((0, 103), (0, 56)), ((0, 114), (0, 106)), ((0, 114), (0, 64)), ((0, 60), (0, 111)), ((0, 60), (0, 106)), ((0, 7), (0, 111)), ((0, 7), (0, 55))])
```
representing all the edges in the dislocation graph.

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
