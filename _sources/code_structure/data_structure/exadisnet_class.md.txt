## ExaDisNet Class

### Description

`ExaDisNet` from `pyexadis_base` is the base class used in the OpenDiS python interface of ExaDiS to represent the dislocation network. `ExaDisNet` is a wrapper class around the internal data structure representation of the dislocation network in ExaDiS, which internally handles memory movements between the different execution spaces (e.g. CPU to GPU).

An `ExaDisNet` network can be instantiated in several ways. The native way is to provide a cell, an array of nodes, and an array of segments as arguments:
```python
G = ExaDisNet(cell, nodes, segs)
```
with
* `cell`: [`pyexadis.Cell`](../../core_libraries/exadis_documentation/user_guide/python_interface/python_binding.rst#pyexadis.Cell) object defining the simulation cell. Constructor arguments:
    - `h` (required): simulation cell matrix (columns are the cell vectors)
    - `origin` (optional): cell origin. Default: (0,0,0)
    - `is_periodic` (optional): periodic flag along the three dimensions. Default: [true,true,true].
* `nodes`: array of nodes where each row contains a node attributes.
    * Attributes for all nodes must be of the following formats:
        - x, y, z
        - x, y, z, constraint
        - domain_id, local_id, x, y, z
        - domain_id, local_id, x, y, z, constraint
    * where x, y, z are the nodes coordinates, constraint is the node constraint (`pyexadis_base.NodeConstraints.UNCONSTRAINED` or `pyexadis_base.NodeConstraints.PINNED_NODE`), domain_id is the simulation domain index, local_id is the local index of the node in the domain.
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
Some utility functions are also provided in file `python/pyexadis_utils.py` to generate basic dislocation graphs (Frank-Read source, infinite lines, etc.). See [Creating initial dislocation configurations](../../../tutorials/initial_configuration) for more information.

Another convenient method is to initialize a `ExaDisNet` object by reading a dislocation network in legacy ParaDiS format from file using built-in method `read_paradis()`:
```python
G = ExaDisNet().read_paradis('config.data')
```

A dislocation network defined in a `ExaDisNet` object must be wrapped into a [`DisNetManager`](disnetmanager_class.md) object before it can be used within modules, e.g.
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```


### Attributes and Methods

```{eval-rst}
.. py:class:: pyexadis_base.ExaDisNet(DisNet_Base)
    :noindex:

    Wrapper class for pyexadis dislocation network.  
    Provides basic functions to manipulate the network.
    The underlying dislocation network object is stored in attribute :attr:`net`.

    **Attributes**
    
    .. py:attribute:: net
       :type: pyexadis.ExaDisNet
       :noindex:
       
       Get the pointer to the ExaDiS network binding object
    
    .. py:property:: cell
        :noindex:

        Get the simulation cell.

        :rtype: pyexadis.Cell

    **Constructor**

    .. py:method:: __init__(*args)
        :noindex:

        Construct ExaDisNet object by instantiating the :attr:`net` attribute.

        - If 3 arguments: `cell`, `nodes`, `segs`
        :param cell: pyexadis.Cell object
        :param nodes: node array
        :param segs: segment array
        
        - If 1 argument: `net`
        :param net: pyexadis.ExaDisNet object
        
        - If 0 arguments:
        Creates empty pyexadis.ExaDisNet
        
        :raises ValueError: if number of arguments is invalid

    **Import/Export methods**

    .. py:method:: import_data(data)
        :noindex:

        Import network data from a dictionary. Argument data must be the output of an export_data() method.

        :param data: Dictionary with keys "cell", "nodes", "segs"
        :returns: self

    .. py:method:: export_data()
        :noindex:

        Export network data as a dictionary.

        :returns: Dictionary with keys "cell", "nodes", "segs"

    .. py:method:: read_paradis(datafile)
        :noindex:

        Read network from a Paradis-format data file.

        :param datafile: Path to data file (str)
        :returns: self

    .. py:method:: read_data(datafile)
        :noindex:

        Alias for :meth:`read_paradis`.

    .. py:method:: write_data(datafile)
        :noindex:

        Write network to a Paradis-format data file.

        :param datafile: Path to data file (str)

    **Generation methods**

    .. py:method:: generate_prismatic_config(crystal, Lbox, num_loops, radius, maxseg=-1, Rorient=None, seed=1234, uniform=False)
        :noindex:

        Generate a prismatic loop configuration.

        :param crystal: Crystal object
        :param Lbox: Box dimensions or cell object
        :param num_loops: Number of loops (int)
        :param radius: Loop radius (float)
        :param maxseg: Maximum segment size (int, default -1)
        :param Rorient: Orientation matrix (optional)
        :param seed: RNG seed (int, default 1234)
        :param uniform: Uniform distribution (bool, default False)
        :returns: self

    .. py:method:: generate_line_config(crystal, Lbox, num_lines, theta=None, maxseg=-1, Rorient=None, seed=-1, verbose=True)
        :noindex:

        Generate a line configuration.

        :param crystal: Crystal object
        :param Lbox: Box dimensions or cell object
        :param num_lines: Number of lines (int)
        :param theta: Line character angle (optional)
        :param maxseg: Maximum segment size (int, default -1)
        :param Rorient: Orientation matrix (optional)
        :param seed: RNG seed (int, default -1)
        :param verbose: Print info (bool, default True)
        :returns: self

    **Network information and access methods**

    .. py:method:: num_nodes()
        :noindex:

        Number of nodes in the network.

        :rtype: int

    .. py:method:: num_segments()
        :noindex:

        Number of segments in the network.

        :rtype: int

    .. py:method:: is_sane()
        :noindex:

        Check network sanity.

        :rtype: bool

    .. py:method:: get_segs_data()
        :noindex:

        Get segment data as a dictionary.

        :returns: Dictionary with keys "nodeids", "burgers", "planes"

    .. py:method:: get_nodes_data()
        :noindex:

        Get node data as a dictionary.

        :returns: Dictionary with keys "tags", "positions", "constraints"

    .. py:method:: get_tags()
        :noindex:

        Get node tags (domain,index).

        :rtype: ndarray, size=(num_nodes,2)

    .. py:method:: get_positions()
        :noindex:

        Get node positions.

        :rtype: ndarray, size=(num_nodes,3)

    .. py:method:: get_forces()
        :noindex:

        Get node forces.

        :rtype: ndarray, size=(num_nodes,3)

    .. py:method:: get_velocities()
        :noindex:

        Get node velocities.

        :rtype: ndarray, size=(num_nodes,3)

    .. py:method:: set_positions(pos)
        :noindex:

        Set node positions.

        :param pos: ndarray of positions, size=(num_nodes,3)
       
```
