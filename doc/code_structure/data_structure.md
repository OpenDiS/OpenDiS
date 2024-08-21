### Data Structure

#### DisNetManager Class

A `DisNetManager` object is a container that allows different core implementations of dislocation network data structures to co-exist and interact with each other within the OpenDiS framework. This allows for inter-operability between modules coming from different core libraries (e.g. PyDiS and ExaDiS). Before any core network object, e.g. `DisNet` or `ExaDisNet`, can be used in modules, they need to be wrapped into a `DisNetManager` object, e.g.:
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```
Then, by specification, each OpenDiS module receives a `DisNetManager` object as input, and can convert it to its desired internal data structure via the `get_disnet()` method as needed, e.g. for an `ExaDisNet`-based module:
```python
class MyExaDisModule():
    def Foo(self, N: DisNetManager, state: dict):
        G = N.get_disnet(ExaDisNet) # convert network to ExaDisNet object
        # perform some operations on the ExaDisNet object
```
Internally, `DisNetManager` retains the type of the active network instance to minimize the number of conversion and memory transfer operations.

The `DisNetManager` provides the `export_data()` method that returns the raw network data stored in a dictionary:
```python
data = N.export_data()
```
where `data` is a dictionary containing the following entries:
* `data["cell"]`: dictionary defining the simulation cell, containing the following entries:
    - `data["cell"]["h"]`: simulation cell matrix (rows are the cell vectors)
    - `data["cell"]["origin"]`: simulation cell origin
    - `data["cell"]["is_periodic"]`: periodic flag along the three dimensions
* `data["nodes"]`: dictionary containing the nodes attributes within the following entries:
    - `data["nodes"]["tags"]`: array of nodes tags (domain,index), size=(Nnodes,2)
    - `data["nodes"]["positions"]`: array of nodes positions, size=(Nnodes,3)
    - `data["nodes"]["constrains"]`: array of nodes constraints, size=(Nnodes,1)
* `data["segs"]`: dictionary containing the segments attributes within the following entries:
    - `data["segs"]["nodeids"]`: array of indices of segments end-nodes (node1,node2), size=(Nsegs,2)
    - `data["segs"]["burgers"]`: array of segments Burgers vectors, size=(Nsegs,3)
    - `data["segs"]["plane"]`: array of segments plane normal, size=(Nsegs,3)




#### DisNet Class
DisNet is the fundamental data structure used in PyDiS to represent the dislocation network.  The following is a sketch of the DisNet data structure, using the initial condition of the [Frank-Read source](../tutorials/frank_read_src/frank_read_src_by_python.md) simulation as an example.  DisNet represents the dislocation network as an undirected graph, meaning that if node A is connected to node B, then node B is also connected to node A.  However, each edge carries a Burgers vector (as an attribute) that depends on from which node we refer to the edge (see below).
```{figure} DisNet_data_struct.png
:alt: Data structure of DisNet
:width: 640px
```
In this example, the DisNet represents a (prismatic) dislocation loop with 5 nodes and 5 edges (called segments).  In DisNet, each node is uniquely specified by a ```Tag```, which is a tuple of two integers, ```(domainID, index)```, following the convention of ParaDiS.  (In this example, the ```domainID``` of all nodes is zero.)  The Tags of the 5 nodes are: (0,0), (0,1), (0,2), (0,3), (0,4).  Each node also carries attributes: ```R``` (3D vector for its position) and ```constraint``` (integer, 0 for unconstraint, 7 for fixed).  The tags and attributes for every node are represented by the yellow boxes in the figure above.

Each segment in DisNet is specified by a tuple of tags specifying the nodes it connects with.  For example, segment ((0,0), (0,1)) represents the segment connecting nodes (0,0) and (0,1).  Since DisNet is an undirected graph, ((0,1), (0,0)) would refer to the same segment.  Each segment has its attributes: ```burg_vec``` (3D vector for Burgers vector) and ```plane_normal``` (3D vector for glide plane normal).  Importantly, the sign of the Burgers vector of a segment depends on the order of the nodes we use to refer to the segment.  For example, the Burger vector of the above segment going from node (0,0) is (1,0,0).  The Burgers vector of the same segment going from node (0,1) is (-1,0,0).  (The Python function to retrieve its Burgers vector of a segment is burg_vec_from(self, from_tag), where ```from_tag``` specifies the node we are going from.)  The tag pairs and attributes for every segment are represented by the green boxes in the figure above.

Section [Explore Dislocation Network](../tutorials/frank_read_src/frank_read_src_by_python.md#explore-dislocation-network) provides an example of how to interact with a DisNet object in Python.

#### ExaDisNet Class
In `pyexadis`, dislocation networks are defined as instances of the `ExaDisNet` class from `python/pyexadis_base.py`. `ExaDisNet` is a wrapper class around the internal data structure representation of the dislocation network in ExaDiS, which internally handles memory movements between the different execution spaces (e.g. CPU to GPU).

An `ExaDisNet` network can be instantiated in several ways. The native way is to provide a cell, an array of nodes, and an array of segments as arguments:
```python
G = ExaDisNet(cell, nodes, segs)
```
with
* `cell`: `pyexadis.Cell` object defining the simulation cell. Constructor arguments:
    - `h` (required): simulation cell matrix (rows are the cell vectors)
    - `origin` (optional): cell origin. Default: (0,0,0)
    - `is_periodic` (optional): periodic flag along the three dimensions. Default: [true,true,true].
* `nodes`: array of nodes where each row contains a node attributes.
    * Attributes for all nodes must be of the following formats:
        - x, y, z
        - x, y, z, constraint
        - domain, local_id, x, y, z
        - domain, local_id, x, y, z, constraint
    * where x, y, z are the nodes coordinates, constraint is the node constraint (`pyexadis_base.NodeConstraints.UNCONSTRAINED` or `pyexadis_base.NodeConstraints.PINNED_NODE`), domain is the simulation domain, local_id is the local id of the node in the domain.
* `segs`: array of segments defining the directed dislocation graph.
    * Segments must be defined only once, e.g. if a segment from node i to node j is defined, then the segment from node j to node i must not be defined.
    * Attributes of the segments must be of the following formats:
        - n1, n2, bx, by, bz
        - n1, n2, bx, by, bz, nx, ny, nz
    * where n1, n2 are the end nodes indices (index in the `nodes` array), bx, by, bz are components of the Burgers vector when going from node n1 to node n2, nx, ny, nz are components of the segment slip plane normal.
        
Note that all lengths (e.g. cell size, nodes coordinates, Burgers vectors) are defined in units of global parameter `burgmag`. Burgers vectors and plane normals are defined in the global frame of the simulation.

For instance, we can define a dislocation line of length 100b along [1,0,0] lying on plane [0,1,0] and discretized with 3 nodes:
```python
Ldis = 100.0
Lbox = 2*Ldis
cell = pyexadis.Cell(h=Lbox*np.eye(3), origin=-0.5*Lbox*np.ones(3))
nodes = np.array([[-0.5*Ldis, 0.0, 0.0, NodeConstraints.PINNED_NODE],
                  [ 0.0*Ldis, 0.0, 0.0, NodeConstraints.UNCONSTRAINED],
                  [ 0.5*Ldis, 0.0, 0.0, NodeConstraints.PINNED_NODE]])
segs = np.array([[0, 1, b[0], b[1], b[2], 0.0, 1.0, 0.0],
                 [1, 2, b[0], b[1], b[2], 0.0, 1.0, 0.0]])
G = ExaDisNet(cell, nodes, segs)
```
Some utility functions are also provided in file `python/pyexadis_utils.py` to generate basic dislocation graphs (Frank-Read source, infinite lines, etc.).

Another convenient method is to initialize a `ExaDisNet` object by reading a dislocation network in legacy ParaDiS format from file using built-in method `read_paradis()`:
```python
G = ExaDisNet()
G.read_paradis('config.data')
```

Then, the dislocation network defined in a `ExaDisNet` object must be wrapped into a `DisNetManager` object before it can be used within modules, e.g.
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```

DisNetManager provides a convenient way to convert the data between DisNet and ExaDisNet.
Section [Graph Data Conversion](../tutorials/frank_read_src/graph_data_conversion.md) describes how to convert between DisNet and ExaDiSNet objects.
