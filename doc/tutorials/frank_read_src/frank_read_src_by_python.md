### Frank-Read Source by Pure Python

We can run the following test case without having to compile OpenDiS.

To run the simulation, simply execute:

```bash
cd ~/Codes/OpenDiS.git
cd examples/02_frank_read_src
python3 -i test_frank_read_src_pydis.py
```

#### Initial Condition
The initial configuration of this simulation is a rectangular (prismatic) dislocation loop, where all four segments are of pure edge type.  Their Burgers vectors are normal to the plane of the rectangular loop.  All four corner nodes have their constraint == FIXED.  This is sufficient to pin the bottom and two side arms of the prismatic loop.  However, the top arm contains a node (in the middle) that is free to move.  This allows the top arm to bow out (under the action of applied stress) and act as a Frank-Read source.

#### Boundary Condition

Periodic boundary condition (PBC) is applied in all three directions, as specified in the following line of the test script.
```python
   net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)
```

#### Simulation Behavior
If the simulation runs successfully, a (Python Matplotlib) window should open and the final dislocation configuration at the end of the simulation should look something like this.
```{figure} frank_read_final_config.png
:alt: Screenshot of the final configuration
:width: 552px
```

Here is what you should see during the simulation.
```{figure} frank_read_opendis.gif
:alt: Video of the whole simulation
:width: 552px
```

```{hint}
If you do not see a new window displaying the dislocation configuration, it is possible because you are running on a remote server and didn't have the X11 forwarding set up correctly. 
 Add the ```-Y``` option in your ```ssh``` command, such as  ```ssh -Y <your_account>@your.cluster.address```.  You can also try the ```xeyes``` command (which will open a new window) to verify that your X11 channel is set up correctly.
```




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

One can export the data that stores the information in this graph:
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

