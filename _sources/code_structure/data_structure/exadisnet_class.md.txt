#### ExaDisNet Class
In `pyexadis`, dislocation networks are defined as instances of the `ExaDisNet` class from `python/pyexadis_base.py`. `ExaDisNet` is a wrapper class around the internal data structure representation of the dislocation network in ExaDiS, which internally handles memory movements between the different execution spaces (e.g. CPU to GPU).

An `ExaDisNet` network can be instantiated in several ways. The native way is to provide a cell, an array of nodes, and an array of segments as arguments:
```python
G = ExaDisNet(cell, nodes, segs)
```
with
* `cell`: `pyexadis.Cell` object defining the simulation cell. Constructor arguments:
    - `h` (required): simulation cell matrix (columns are the cell vectors)
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
Section [Graph Data Conversion](../../tutorials/frank_read_src/graph_data_conversion.md) describes how to convert between DisNet and ExaDiSNet objects.
