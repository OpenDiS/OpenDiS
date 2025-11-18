## Code architecture and data structures

### Overview

ExaDiS is implemented based on a modular design in which each functionality can be viewed as an independent module. This design is essential to make the code extensible and enable its coupling with various other modules within the [OpenDiS](https://github.com/OpenDiS/) framework.

ExaDiS core backend is written in modern C++ and built based on the [Kokkos](https://kokkos.org) framework. This allows to maintain a unique version of the source code that is largely agnostic of the target hardware: the same source code can be compiled for a wide variety of systems and hardwares, e.g. for serial CPU machines, using OpenMP, or for GPU machines using cuda or hip architectures, for instance.

ExaDiS classes and functions are implemented with high-performance and high-parallelism kernel execution in mind. Yet, the code is designed to abstract away as much complexity as possible, so as to lower the entry barrier for developing new functionalities. For prototyping or first implementation pass, ExaDiS also provides convenient helper functions and data structures that are easier to manipulate than Kokkos-based implementations. An illustration of this is the `Topology` class implementation, where class `TopologySerial` implements the split-multi-node procedure in a serial fashion on the host (CPU), while `TopologyParallel` is an implementation of the same split-multi-node procedure but executed with highly-parallel kernels on the device (GPU).


### Project structure

Here is a brief description of the directory structure of the ExaDiS code:

- `cmake/`: **cmake** related files, including pre-defined build system options
- `examples/`: examples of scripts and simulations
    - each example is numbered and placed in a dedicated subfolder
- `kokkos/`: directory containing the **Kokkos** submodule
- `python/`: files related to the python binding implementation
    - `pybind11/` directory containing the **pybind11** submodule
    - `exadis_pybind.cpp`: implementation of the C++ classes / functions binding
    - `pyexadis_base.py`: interface enabling the use of **pyexadis** within OpenDiS
    - `pyexadis_utils.py`: **pyexadis** utility functions
- `src/` : C++ source and header files (`*.cpp`, `*.h`)
    - base class files and common functions are at the root of the `src` directory
    - individual module implementations are placed in subfolders
    - `collision_types/`: implementation of collision modules
    - `force_types/`: implementation of force modules
    - `integrator_types/`: implementation of integrator modules
    - `mobility_types/`: implementation of mobility modules
    - `neighbor_types/`: implementation of neighbor modules
    - `topology_types/`: implementation of topology modules
- `tests/`: files for testing and debugging
    - `benchmark/`: benchmark tests to evaluate the performance of the code
    - `debug.h`: debug utility functions
    - `test_exadis.cpp`: simple ExaDiS test simulations
    


### Dislocation network and data structures

Dislocation networks in ExaDiS are stored using two major classes, `SerialDisNet` and `DeviceDisNet`. Class `DisNetManager` is used as a container to synchronize dislocation networks between the two classes for use in the different execution spaces (host CPU / device GPU).

#### SerialDisNet class
`SerialDisNet` is a STL-based class that allows to easily create, manipulate, and modify dislocation networks on the host: a `SerialDisNet` instance is marked for a `Kokkos::Serial` execution space associated with a `Kokkos::HostSpace` memory space. The `SerialDisNet` class implements all low-level topological operations, e.g. `split_seg()`, `split_node()`, and `merge_nodes()` methods. It is designed to create initial dislocation networks and perform topological operations on existing networks.

##### Constructors
- `SerialDisNet()`: instantiate an empty dislocation network
- `SerialDisNet(double Lbox)`: instantiate an empty dislocation network within a cubic cell of side `Lbox`
- `SerialDisNet(const Cell& cell)`: instantiate an empty dislocation network within a cell `Cell`

##### Properties
- `Cell cell`: dislocation network cell
- `std::vector<DisNode> nodes`: array of nodes
- `std::vector<DisSeg> segs`: array of segments
- `std::vector<Conn> conn`: array of node connectivity
- `int Nnodes_local`: number of local nodes
- `int Nsegs_local`: number of local segments

##### Methods
- `int number_of_nodes()`: returns the number of nodes in the network
- `int number_of_segs()`: returns the number of segments in the network
- `NodeTag get_new_tag()`: returns a new available node tag
- `void free_tag(NodeTag& tag)`: release an existing node tag
- `void update_ptr()`: updates pointers to the nodes, segments and connectivity data as STL arrays may be reallocated under-the-hood, e.g. due to resize during topology changes
- `DisNode* get_nodes()`: returns a pointer to the nodes array data
- `DisSeg* get_segs()` returns a pointer to the segments array data
- `Conn* get_conn()` returns a pointer to the connectivity array data
- `void add_node(const Vec3& pos)`: add a node by specifying its position `pos`
- `void add_node(const Vec3& pos, int constraint)`: add a node by specifying its position `pos` and its constraint `constraint`
- `void add_node(const NodeTag& tag, const Vec3& pos)`: add a node by specifying its tag `tag` and position `pos`
- `void add_node(const NodeTag& tag, const Vec3& pos, int constraint)`: add a node by specifying its tag `tag`, position `pos`, and its constraint `constraint`
- `void add_seg(int n1, int n2, const Vec3& b)`: add a segment with Burgers vector `b` connecting node at index `n1` with node at index `n2`
- `void add_seg(int n1, int n2, const Vec3& b, const Vec3& p)`: add a segment with Burgers vector `b` and plane nornal `p` connecting node at index `n1` with node at index `n2`
- `int find_connection(int n1, int n2)`: finds in the node connectivity array of `n1` if there exists a connection (segment) to node at index `n2`
- `void generate_connectivity()`: generate the connectivity array for all nodes
- `double seg_length(int i)`: returns the length of segment `i`
- `bool constrained_node(int i)`: returns if node `i` is a constrained node
- `bool discretization_node(int i)`: returns if node `i` is a discretization node
- `void update_node_plastic_strain(int i, const Vec3& pold, const Vec3& pnew, Mat33& dEp)`: update the plastic strain tensor `dEp` with the increment of moving node `i` from position `pold` to position `pnew`
- `int split_seg(int i, const Vec3& pos, bool update_conn=true)`: splits segment `i` by inserting a new node at position `pos` and returns the index of the new node
- `int split_node(int i, std::vector<int> arms)`: splits node `i` into a new node that contains the subset of arms indices `arms`, and returns the index of the new node
- `bool merge_nodes(int n1, int n2, Mat33& dEp)`: merges nodes `n1` and `n2` into node `n1`, updates the corresponding plastic strain tensor `dEp`, and returns whether the merge succeeded
- `bool merge_nodes_position(int n1, int n2, const Vec3 &pos, Mat33& dEp)`: merges nodes `n1` and `n2` into node `n1` at position `pos`, updates the corresponding plastic strain tensor `dEp`, and returns whether the merge succeeded
- `void remove_segs(std::vector<int> seglist)`: removes all segments with indices specified in `seglist` from the network
- `void remove_nodes(std::vector<int> nodelist)`: removes all nodes with indices specified in `nodelist` from the network
- `void purge_network()`: purge the network from unconnected nodes and zero Burgers vectors segments
- `std::vector<std::vector<int> > physical_links()`: decompose the network into a set of dislocation links connecting physical network nodes
- `double dislocation_density(double burgmag)`: computes the dislocation density in units of 1/m^2 given the Burgers vector scaling factor `burgmag`
- `void write_data(std::string filename)`: write the dislocation network into a file using the legacy ParaDiS `.data` format
- `SaveNode save_node(int i)`: save a dislocation node and all its connections
- `void restore_node(SaveNode& saved_node)`: restore a dislocation node and all its connections from a save node `saved_node`


#### DeviceDisNet class
`DeviceDisNet` is a class that uses Kokkos views to store and access dislocation nodes, segments and connections for device execution (GPU): a `DeviceDisNet` instance is marked for a `Kokkos::DefaultExecutionSpace` execution space associated with corresponding default memory space. These will default to the highest available execution/memory spaces available at compile time (e.g. device spaces when compiling for GPU). As of now, the `DeviceDisNet` class does not implement topological operations.

##### Constructor
- `DeviceDisNet(const Cell& cell)`: instantiates an empty dislocation network within a cell `Cell`

##### Properties
- `Cell cell`: dislocation network cell
- `T_nodes nodes`: Kokkos view of nodes
- `T_segs segs`: Kokkos view of segments
- `T_conn conn`: Kokkos view of node connectivity
- `int Nnodes_local`: number of local nodes
- `int Nsegs_local`: number of local segments

##### Methods
- `void update_ptr()`: updates network pointers (dummy for `DeviceDisNet`)
- `T_nodes::pointer_type get_nodes()`: returns a pointer to the nodes view data
- `T_segs::pointer_type get_segs()` returns a pointer to the segments view data
- `T_conn::pointer_type get_conn()` returns a pointer to the connectivity view data


#### DisNetManager class
Class `DisNetManager` is used as a container to synchronize dislocation networks between the two classes `SerialDisNet` and `DeviceDisNet` for use in the different execution spaces. A `DisNetManager` object is instantiated by providing a `SerialDisNet` (resp. `DeviceDisNet`) object, and a mirror `DeviceDisNet` (resp. `SerialDisNet`) object is automatically created. A given instance of a network type is requested by using functions `get_serial_network()` and `get_device_network()`. When calling these functions, a memory copy is triggered only if the requested network type instance is not marked as active to minimize memory transfers between host and device memory spaces. `DisNetManager` is used as the fundamental network class associated with the ExaDiS `System` object that is being propagated through the different ExaDiS modules.

##### Constructors
- `DisNetManager(SerialDisNet* n)`: instantiates a `DisNetManager` object from a `SerialDisNet` dislocation network `n`
- `DisNetManager(DeviceDisNet* d)`: instantiates a `DisNetManager` object from a `DeviceDisNet` dislocation network `d`

##### Methods
- `SerialDisNet* get_serial_network()`: returns a pointer to a `SerialDisNet` instance of the dislocation network
- `DeviceDisNet* get_device_network()`: returns a pointer to a `DeviceDisNet` instance of the dislocation network
- `int get_active()`: returns the currently active network type (`SERIAL_ACTIVE` or `DEVICE_ACTIVE`)
- `void set_active(int a)`: sets as active the network type `a = SERIAL_ACTIVE` or `a = DEVICE_ACTIVE`
- `int Nnodes_local()`: returns the number of local nodes in the network
- `int Nsegs_local()`: returns the number of local segments in the network


### System class
`System` is the base class in ExaDiS that contains all information about the simulated dislocation system, including the parameters, the crystal instance, and the dislocation network object. A `System` object is the fundamental data structure that is being propagated from modules to modules.

```{important}
A `System` object must be allocated using the `exadis_new()` or the `make_system()` helper functions to ensure it is placed on a memory space accessible to all execution spaces.
```

##### Constructors
- `System()`: instantiates an empty `System` object
- `System* make_system(SerialDisNet* net, Crystal crystal, Params params)`: creates a `System` object from an initial dislocation network, crystal instance, and parameters object.

##### Properties
- `DisNetManager* net_mngr`: dislocation network manager of the system
- `double neighbor_cutoff`: neighbor cutoff for the simulated system
- `Params params`: parameters of the simulated system
- `Crystal crystal`: crystal object of the simulated system
- `T_x xold`: Kokkos view to store old nodal positions
- `Mat33 extstress`: external/applied stress tensor
- `double realdt`: current global time step size of the simulation
- `Mat33 dEp`: current increment of plastic strain for the time step
- `Mat33 dWp`: current increment of plastic spin for the time step
- `double density`: current dislocation density in the system
- `SystemTimer timer[]`: array of timers of the system
- `int numdevtimer`: number of development timers of the system
- `SystemTimer devtimer[]`: array of development timers of the system

##### Methods
- `void initialize(Params _params, Crystal _crystal, SerialDisNet* network)`: initialize a `System` object with a parameters object, a crystal instance, and an initial dislocation network
- `void register_neighbor_cutoff(double cutoff)`: register a minimum cutoff distance to be used in the simulated system
- `SerialDisNet* get_serial_network()`: returns a pointer to a `SerialDisNet` instance of the dislocation network of the system
- `DeviceDisNet* get_device_network()`: returns a pointer to a `DeviceDisNet` instance of the dislocation network of the system
- `int Nnodes_local()`: returns the local number of dislocation nodes in the system
- `int Nsegs_local()`: returns the local number of dislocation segments in the system
- `int Nnodes_total()`: returns the total number of dislocation nodes in the system
- `int Nsegs_total()`: returns the total number of dislocation segments in the system
- `void plastic_strain()`: computes the plastic strain resulting from the motion of all dislocation lines between the nodal positions stored in `xold` and the current positions
- `void reset_glide_planes()`: reset and select appropriate glide planes for all segments in the network
- `void write_config(std::string filename)`: writes the network configuration into a file (in legacy ParaDiS format for now)
- `void print_timers(bool dev=false)`: prints system timers
