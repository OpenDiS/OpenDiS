### Data Structure


#### DisNet Class
DisNet is the fundamental data structure used in PyDiS to represent the dislocation network.  The following is a sketch of the DisNet data structure, using the initial condition of the [Frank-Read source](../tutorials/frank_read_src/frank_read_src_by_python.md) as an example.  DisNet represents the dislocation network as an undirected graph, meaning that if node A is connected to node B, then node B is also connected to node A.  However, each edge carries a Burgers vector (as an attribute) that depends on from which node we refer to the edge (see below).
```{figure} DisNet_data_struct.png
:alt: Data structure of DisNet
:width: 640px
```
In this example, the DisNet represents a (prismatic) dislocation loop with 5 nodes and 5 edges (called segments).  

#### ExaDisNet and DisNetManager
ExaDisNet is what ExaDis uses to represent the dislocation network.

DisNetManager provides a convenient way to convert the data between DisNet and ExaDisNet.

