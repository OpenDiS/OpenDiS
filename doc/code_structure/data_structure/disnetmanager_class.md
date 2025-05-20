## DisNetManager Class

### Description

`DisNetManager` is a container class that allows different core implementations of dislocation network data structures to co-exist and interact with each other within the OpenDiS framework. This allows for inter-operability between modules coming from different core libraries (e.g. PyDiS and ExaDiS). Before any core network object, e.g. `DisNet` or `ExaDisNet`, can be used in modules, they need to be wrapped into a `DisNetManager` object, e.g.:
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```
Then, by specification, each OpenDiS module receives a `DisNetManager` object as input, and can convert it to its desired internal data structure via the `get_disnet()` method as needed. For instance, an ExaDiS-based module can request an instance of an `ExaDisNet` object from the input `DisNetManager` object:
```python
class MyExaDisModule():
    def Foo(self, N: DisNetManager, state: dict):
        G = N.get_disnet(ExaDisNet) # convert network to ExaDisNet object
        # perform some operations on the ExaDisNet object
```
Internally, `DisNetManager` performs the conversion between different network object types by invoking the `export_data()` / `import_data()` methods implemented in the core network classes. It also retains the type of the active network instance to minimize the number of conversion and memory transfer operations. See section [Graph Data Conversion](../../tutorials/frank_read_src/graph_data_conversion.md) for an example of how to convert between DisNet and ExaDiSNet objects.

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


Core network classes [`DisNet`](disnet_class.md) from PyDiS and [`ExaDisNet`](exadisnet_class.md) from ExaDiS are the native data structures used to handle dislocation networks in OpenDiS. However, if desired, a user can use or implement its own data structure to interact with dislocation networks in OpenDiS, as long as the new data structure (or wrapper to it) abides by the [DisNet_Base](https://github.com/OpenDiS/OpenDiS/blob/main/python/framework/disnet_base.py) abstract class.


### Properties
- `DisNetManager.cell`: network cell object


### Methods
- `DisNetManager.get_disnet(disnet_type)`: Converts the dislocation network into the requested `disnet_type` object.
- `DisNetManager.get_active_type()`: Returns the type of network object that is currently active.
- `DisNetManager.export_data()`: Exports the raw network data into a `data` dictionary.
- `DisNetManager.import_data(data)`: Set the content of the network by importing it from a `data` dictionary. Argument `data` must be the output of an `export_data()` method.
- `DisNetManager.write_json(filename)`: Writes the DisNetManager data to a JSON file.
- `DisNetManager.read_json(filename)`: Reads the DisNetManager data from a JSON file.
- `DisNetManager.num_nodes()`: Returns the number of nodes in the network.
- `DisNetManager.num_segments()`: Returns the number of segments in the network.
- `DisNetManager.is_sane()`: Checks if the network is sane.
