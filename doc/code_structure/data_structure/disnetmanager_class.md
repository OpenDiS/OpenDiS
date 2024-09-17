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
    - `data["cell"]["h"]`: simulation cell matrix (columns are the cell vectors)
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

DisNetManager provides a convenient way to convert the data between DisNet and ExaDisNet.
Section [Graph Data Conversion](../../tutorials/frank_read_src/graph_data_conversion.md) describes how to convert between DisNet and ExaDiSNet objects.
