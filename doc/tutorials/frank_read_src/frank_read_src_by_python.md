### Frank-Read Source by Pure Python

Frank-Read Source is a mechanism of the generation of dislocations during the plastic deformation of solids. It explains the increase of dislocation density during the material's deformation process. This simple example demonstrates how to simulate Frank-Read Source in pure python mode (PyDiS) using OpenDiS.

To run the simulation, simply execute:

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src
python3 -i test_frank_read_src.py
```

In the simulation, there are two dislocation jogs pinning the four dislocation nodes. The bottom two nodes are fixed so that the dislocation line glides on the top plane. The propagation of the dislocation subsumes force calculation, collision, mobility law, re-meshing rules, time integration, and topological operation, as defined in ```pydis``` (Chap. 10 of [Bulatov & Cai, 2006](https://core.ac.uk/reader/44178170)). The dislocation source propagates and interacts with itself, resulting in a dislocation loop and regeneration of a dislocation that can repeat the sequence ([MAE 324, Princeton University](https://www.princeton.edu/~maelabs/mae324/glos324/frankreed.htm))

```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.
2. how to examine the data structure of the dislocation network (DisNet) in python interactive
```

**Initial condition:** 

The initial condition is described in the function ```init_frank_read_src_loop()```. The initial nodes are described in the ```cell_list```, where the ```CellList``` function is being imported from the DisNet manager to initialize the simulation cell. One defines the simulation cell (```cell = ```). The initial nodal positions are introduced in the variable ```rn```.


**Boundary condition:**

In this example, the periodic boundary condition (PBC) is being implemented. One can turn the ```pbc``` option on or off by setting ```init_frank_read_src_loop(pbc=True)```.

**Explore dislocation data:**

The ```pydis``` module stores the simulated dislocation networks as graphs. One can explore the properties of this graph by checking its nodes

```python
net.G.nodes
```
and obtain
```python
NodeView(((0, 0), (0, 2), (0, 3), (0, 4), (0, 39), (0, 40), (0, 106), (0, 48), (0, 47), (0, 125), (0, 8), (0, 55), (0, 56), (0, 13), (0, 65), (0, 64), (0, 111), (0, 89), (0, 116), (0, 103), (0, 114), (0, 60), (0, 7)))
```
which represents all the nodes in the dislocation graph. Note that each node is labeled as ```(0, i)``` where ```i``` is an integer but not continuously labeled.

One can also obtain the graph's edges:
```python
net.G.edges
```
and obtain
```python
OutEdgeView([((0, 0), (0, 4)), ((0, 0), (0, 48)), ((0, 2), (0, 3)), ((0, 2), (0, 47)), ((0, 3), (0, 2)), ((0, 3), (0, 4)), ((0, 4), (0, 3)), ((0, 4), (0, 0)), ((0, 39), (0, 56)), ((0, 39), (0, 13)), ((0, 40), (0, 65)), ((0, 40), (0, 116)), ((0, 106), (0, 114)), ((0, 106), (0, 60)), ((0, 48), (0, 0)), ((0, 48), (0, 13)), ((0, 47), (0, 2)), ((0, 47), (0, 65)), ((0, 125), (0, 89)), ((0, 125), (0, 103)), ((0, 8), (0, 64)), ((0, 8), (0, 116)), ((0, 55), (0, 89)), ((0, 55), (0, 7)), ((0, 56), (0, 39)), ((0, 56), (0, 103)), ((0, 13), (0, 48)), ((0, 13), (0, 39)), ((0, 65), (0, 47)), ((0, 65), (0, 40)), ((0, 64), (0, 8)), ((0, 64), (0, 114)), ((0, 111), (0, 60)), ((0, 111), (0, 7)), ((0, 89), (0, 125)), ((0, 89), (0, 55)), ((0, 116), (0, 40)), ((0, 116), (0, 8)), ((0, 103), (0, 125)), ((0, 103), (0, 56)), ((0, 114), (0, 106)), ((0, 114), (0, 64)), ((0, 60), (0, 111)), ((0, 60), (0, 106)), ((0, 7), (0, 111)), ((0, 7), (0, 55))])
```
representing all the edges in the dislocation graph that denote the connection between nodes.

One can export the data that stores the information of this graph:
```python
net.G.export_data()
```
where the data is stored as
```python
{'cell': {'h': array([[8., 0., 0.],
       [0., 8., 0.],
       [0., 0., 8.]]), 'is_periodic': [True, True, True]}, 'nodes': [[0.0, -0.5, 0.0, 7], [0.0, 0.5, 0.0, 7], [0.0, 0.5, -1.0, 7], [0.0, -0.5, -1.0, 7], [0.6171914542157114, -0.9203767203483139, 0.0, 0], [0.5896056897739462, 0.9205821223444571, 0.0, 0], [1.527654949085545, 0.37859632672384724, 0.0, 0], [0.10738792556990989, -0.7427245725072399, 0.0, 0], [0.10444403949640058, 0.7398113693714066, 0.0, 0], [1.2394214646843567, -0.7124467613152135, 0.0, 0], [1.038784432520676, 0.8285147317619693, 0.0, 0], [1.5716894184507821, -0.28790456076296594, 0.0, 0], [0.9117097060289435, -0.8758517344507986, 0.0, 0], [0.34637940319295846, -0.8795321491305716, 0.0, 0], [0.3378369332510002, 0.8766402879318167, 0.0, 0], [1.2872062053182376, 0.6724966429085716, 0.0, 0], [1.6233434358862746, 0.05363704629680839, 0.0, 0], [1.429249934659248, -0.5312531663550655, 0.0, 0], [0.816660983695957, 0.9004808461025454, 0.0, 0], [1.0825050529428142, -0.806130063119192, 0.0, 0], [1.4213149895828676, 0.5379245440205022, 0.0, 0], [1.5920293249425574, 0.2218589104669993, 0.0, 0], [1.6172723087695884, -0.11733006641199845, 0.0, 0]], 'segs': [[0, 3, -1.0, -0.0, -0.0, 0.0, -1.0, 0.0], [0, 7, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [1, 2, 1.0, 0.0, 0.0, -0.0, 1.0, 0.0], [1, 8, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [2, 3, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0], [4, 12, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [4, 13, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [5, 14, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [5, 18, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [6, 20, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [6, 21, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [7, 13, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [8, 14, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [9, 17, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [9, 19, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [10, 15, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [10, 18, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [11, 17, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [11, 22, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [12, 19, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [15, 20, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0], [16, 21, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [16, 22, -1.0, -0.0, -0.0, 0.0, 0.0, 1.0]]}
```
One may pass this to a variable that surrogates the dislocation structures after the simulation.

If being ran successfully, one may expect the final simulation configuration to look something like this:
```{figure} frank_read_schematic.png
:alt: Screenshot of the final configuration
:width: 552px
```

The simulation process follows this animation:
```{figure} frank_read_vid.gif
:alt: Video of the whole simulation
:width: 552px
```


```{hint}
In case the window for the visualization is not showing up, please do ```ssh -Y <your_account>@your.address```
```
