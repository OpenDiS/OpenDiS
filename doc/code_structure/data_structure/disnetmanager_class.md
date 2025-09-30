## DisNetManager Class

### Description

`DisNetManager` is a container class that allows different core implementations of dislocation network data structures to co-exist and interact with each other within the OpenDiS framework. This allows for inter-operability between modules coming from different core libraries (e.g. PyDiS and ExaDiS). Before any core network object, e.g. `DisNet` or `ExaDisNet`, can be used in modules, they need to be wrapped into a `DisNetManager` object, e.g.:
```python
from framework.disnet_manager import DisNetManager

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

Core network classes [`DisNet`](disnet_class.md) from PyDiS and [`ExaDisNet`](exadisnet_class.md) from ExaDiS are the native data structures used to handle dislocation networks in OpenDiS. However, if desired, a user can use or implement its own data structure to interact with dislocation networks in OpenDiS, as long as the new data structure (or wrapper to it) adheres to the [DisNet_Base](https://github.com/OpenDiS/OpenDiS/blob/main/python/framework/disnet_base.py) specification class.


### Accessing raw data

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
    - `data["segs"]["planes"]`: array of segments plane normals, size=(Nsegs,3)


### Attributes and methods

```{eval-rst}
.. py:class:: DisNetManager()
    :noindex:

    Class for managing multiple implementations of dislocation network.
    Implements synchronization between different implementations of DisNet.

    **Attributes**
    
    .. py:property:: G
        :noindex:

        Get the active network instance.
    
    .. py:property:: cell
        :noindex:

        Get the simulation cell object.

    **Constructor**

    .. py:method:: __init__(disnet)
        :noindex:

        Constructs a DisNetManager from an existing dislocation network object.
        
        :param disnet: existing dislocation network object

    **Methods**
    
    .. py:method:: get_disnet(disnet_type)
        :noindex:

        Converts the dislocation network into the requested `disnet_type` object.

        :rtype: disnet_type network object
        
    .. py:method:: get_active_type()
        :noindex:

        Returns the type of network object that is currently active.

        :rtype: disnet_type

    .. py:method:: import_data(data)
        :noindex:

        Imports raw network data from a `data` dictionary. Argument `data` must be the output of an `export_data()` method.

        :param data: Dictionary with keys "cell", "nodes", "segs"

    .. py:method:: export_data()
        :noindex:

        Exports the raw network data into a `data` dictionary.

        :returns: Dictionary with keys "cell", "nodes", "segs"

    .. py:method:: write_json(filename)
        :noindex:

        Writes the DisNetManager data to a JSON file.

        :param filename: Path to JSON file (str)

    .. py:method:: read_json(filename)
        :noindex:

        Reads the DisNetManager data from a JSON file.

        :param filename: Path to JSON file (str)

    .. py:method:: num_nodes()
        :noindex:

        Returns the number of nodes in the network.

        :rtype: int

    .. py:method:: num_segments()
        :noindex:

        Returns the number of segments in the network.

        :rtype: int

    .. py:method:: is_sane()
        :noindex:

        Checks if the network is sane.

        :rtype: bool
       
```
