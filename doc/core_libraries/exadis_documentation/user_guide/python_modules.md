## Python modules: pyexadis
This section documents the various ExaDiS modules available through the python interface, `pyexadis`. These modules are bindings to the backend C++ modules implemented in ExaDiS. For documentation about the backend C++ modules please see the [Developer guide](../developer_guide/index) section of the documentation.

### Dislocation network

#### ExaDisNet
Dislocation networks are defined as instances of the `ExaDisNet` class from `python/pyexadis_base.py`. `ExaDisNet` is a wrapper class around the internal data structure representation of the dislocation network in ExaDiS, which internally handles memory movements between the different execution spaces (e.g. CPU to GPU).

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
    * Attributes for all nodes must be one of the following formats:
        - x, y, z
        - x, y, z, constraint
        - domain, local_id, x, y, z
        - domain, local_id, x, y, z, constraint
    * where x, y, z are the nodes coordinates, constraint is the node constraint (`pyexadis_base.NodeConstraints.UNCONSTRAINED` or `pyexadis_base.NodeConstraints.PINNED_NODE`), domain is the simulation domain, local_id is the local index of the node in the domain.
* `segs`: array of segments defining the directed dislocation graph.
    * Segments must be defined only once, e.g. if a segment from node i to node j is defined, then the segment from node j to node i must not be defined.
    * Attributes of the segments must be one of the following formats:
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

Utility method `generate_prismatic_config()` is also provided to create an initial configuration made of prismatic loops:
```python
G = ExaDisNet()
G.generate_prismatic_config(crystal, Lbox, numsources, radius, maxseg, seed)
```

Another convenient method is to initialize a `ExaDisNet` object by reading a dislocation network in legacy ParaDiS format from file using built-in method `read_paradis()`:
```python
G = ExaDisNet()
G.read_paradis('config.data')
```

Important: When used in ExaDiS or OpenDiS modules, the dislocation network defined in a `ExaDisNet` object must first be wrapped into a `DisNetManager` object before it can be used within modules (see next section), e.g.
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```

##### Properties
- `ExaDisNet.net`: pointer to the ExaDiS network binding object
- `ExaDisNet.cell`: network cell object

##### Methods
- `ExaDisNet.import_data(data)`: Set the content of the `ExaDisNet` object by importing it from a `data` dictionary. Argument `data` must be the output of an `export_data()` method.
- `data = ExaDisNet.export_data()`: Export the `ExaDisNet` object into a `data` dictionary.
- `ExaDisNet.read_paradis(datafile)`: Set the content of the `ExaDisNet` object by reading a legacy ParaDiS data file.
- `ExaDisNet.write_data(datafile)`: Write the network into a legacy ParaDiS data file.
- `ExaDisNet.generate_prismatic_config(crystal, Lbox, numsources, radius, maxseg, seed)`: Set the content of the `ExaDisNet` object by generating a configuration made of prismatic dislocation loops.
- `nodes_data = ExaDisNet.get_nodes_data()`: Returns a dictionary of nodes data, containing entries `tags`, `positions`, and `constraints`.
- `ExaDisNet.get_tags()`: Returns an array of the nodes tags (domain,index), size=(Nnodes,2).
- `ExaDisNet.get_positions()`: Returns an array of the nodes positions, size=(Nnodes,3).
- `ExaDisNet.get_forces()`: Returns an array of the nodes forces, size=(Nnodes,3).
- `ExaDisNet.get_velocities()`: Returns an array of the nodes velocities, size=(Nnodes,3).
- `ExaDisNet.get_segs_data()`: Returns a dictionary of segments data, containing entries `nodeids`, `burgers`, and `planes`.


#### DisNetManager

A `DisNetManager` object is a container that allows different core implementations of dislocation network data structures to co-exist and interact with each other within the OpenDiS framework. This allows for inter-operability between modules coming from different core libraries (e.g. PyDiS and ExaDiS). Before any `ExaDisNet` network object can be used in modules, it needs to be wrapped into a `DisNetManager` object:
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```
Then, by specification, each module receives a `DisNetManager` object as input, and can convert it to its desired/internal data structure via the `get_disnet()` method as needed, e.g.
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

##### Methods
- `G = DisNetManager.get_disnet(disnet_type)`: Converts the dislocation network into the desired `disnet_type` object.
- `data = DisNetManager.export_data()`: Export the network into a `data` dictionary.


### Forces

#### CalForce modules
A force module is declared using the `CalForce` class from `python/pyexadis_base.py`. `CalForce` is a wrapper class for computing forces at dislocation nodes.

##### Usage
A `CalForce` module is instantiated by providing the global `state` dictionary, the force mode, and specific parameters to the force mode:
```python
calforce = CalForce(state=state, force_mode='ForceName', ...)
```
Available force modes are:

* `force_mode='DDD_FFT_MODEL'`: DDD-FFT model where forces are computed by summing (i) external PK forces, (ii) core forces, (iii) short-range elastic interactions computed with the non-singular model, and (iv) long-range forces computed using FFT. Specific mode parameters are:
    - `Ngrid` (required): size of the FFT grid along each dimension.
    - `cell` (required): simulation cell object.
    - `Ec` (optional): core energy factor in Pa. If not provided or negative, it defaults to the standard value Ec=mu/4/pi*ln(a/0.1).
    - `Ecore_junc_fact` (optional): ratio of the core energy of junction to native Burgers vectors. Default: 1.0.


* `force_mode='SUBCYCLING_MODEL'`: Force mode required to be used when using  subcycling time-integration. It uses the DDD-FFT model described above. Specific mode parameters are:
    - `Ngrid` (required): size of the FFT grid along each dimension.
    - `cell` (required): simulation cell object.
    - `drift` (optional): toggle to use the drift subcycling scheme. If drift=0 the original subcycling scheme is used, which is not compatible with non-linear mobilities. Default: 1.
    - `Ec` (optional): core energy factor in Pa. If not provided or negative, it defaults to the standard value Ec=mu/4/pi*ln(a/0.1).
    - `Ecore_junc_fact` (optional): ratio of the core energy of junction to native Burgers vectors. Default: 1.0.


* `force_mode='CUTOFF_MODEL'`: Cutoff model where forces are computed by summing (i) external PK forces, (ii) core forces, and (iii) short-range elastic interactions up to a given cutoff distance between segment pairs, computed with the non-singular model. Specific mode parameters are:
    - `cutoff` (required): cutoff for segment pair interaction distance.
    - `Ec` (optional): core energy factor in Pa. If not provided or negative, it defaults to the standard value Ec=mu/4/pi*ln(a/0.1).
    - `Ecore_junc_fact` (optional): ratio of the core energy of junction to native Burgers vectors. Default: 1.0.


* `force_mode='LINE_TENSION_MODEL'` or `force_mode='LineTension'`: line-tension model where forces are computed from (i) the external PK force and (ii) the line energy (core force), and no pairwise elastic interaction is accounted for. Specific mode parameters are:
    - `Ec` (optional): core energy factor in Pa. If not provided or negative, it defaults to the standard value Ec=mu/4/pi*ln(a/0.1).
    - `Ecore_junc_fact` (optional): ratio of the core energy of junction to native Burgers vectors. Default: 1.0.


Nodal forces are computed using the `NodeForce()` method
```python
state = calforce.NodeForce(N, state)
```
The `NodeForce()` method takes a dislocation network `DisNetManager` object and the `state` dictionary as arguments. The `state` dictionary must contain an `applied_stress` entry that specifies the value of the external stress tensor used to compute the external Peach-Koehler force contribution to the nodal forces. The `applied_stress` value must be an array of size 6 that contains the external stress tensor components in standard Voigt notation (xx,yy,zz,yz,xz,xy). The method returns an updated `state` dictionary with entry `nodeforces` that contains the array of nodal forces and corresponding node tags stored in entry `nodeforcetags`:
```python
forces = state["nodeforces"]
tags = state["nodeforcetags"]
print(forces[0]) # print force on node tags[0]
```
The node tags are used to keep track of the nodes across modules, as different modules (e.g. in OpenDiS) may use different data structures and may shuffle the nodes ordering when being converted from one data structure to another via the `DisNetManager` object.

##### Properties
- `CalForce.force_mode`: name of the force mode
- `CalForce.params`: ExaDiS global parameters used to setup the force model
- `CalForce.force`: pointer to the ExaDiS force binding object

##### Methods
- `state = CalForce.NodeForce(DisNetManager, state)`: Computes nodal forces for all nodes and place the result in the `state` dictionary.
- `state = CalForce.PreCompute(DisNetManager, state)`: Call to the force pre-compute operation, which may be required before calls to the `OneNodeForce()` method. This is because some force modules (e.g. FFT) requires certain values to be pre-computed and up-to-date before forces can be evaluated efficiently. There is no need to call `PreCompute()` before a call to the global `NodeForce()` method, as the pre-compute step is always called internally in `NodeForce()`.
- `f = CalForce.OneNodeForce(DisNetManager, state, tag, update_state)`: Compute the nodal force `f` on a single node specified by its `tag`. The `PreCompute()` operation may need to be called before calling `OneNodeForce()`.


#### Force calculation wrappers
In addition to `CalForce` modules, `pyexadis` provides wrappers to some internal ExaDiS force calculation functions, which can be called from python with `pyexadis.<function>(...)`. The following functions are currently available:
- `f = pyexadis.compute_force_n2(ExaDisNet, mu, nu, a)`: compute elastic interaction forces using a brute-force N^2 segment pair summation.
- `f = pyexadis.compute_force_cutoff(ExaDisNet, mu, nu, a, cutoff, maxseg)`: compute elastic forces including interactions up to a given cutoff distance between segment pairs. Argument `maxseg` is needed to ensure that no segment pair is missed when building the neighbor list.
- `f = pyexadis.compute_force_segseglist(ExaDisNet, mu, nu, a, segseglist)`: compute elastic forces provided a list user-defined list of segment pairs. Argument `segseglist` must be an array of size=(Npair,2) indicating the indices of each pair of segments to consider in the calculation. Each pair must be unique, e.g. if pair (2,3) is included, then pair (3,2) must not be included or its contribution will be double-counted. Self-interactions, e.g. pair (2,2), must not be included.



### Mobility laws
A mobility module is declared using the `MobilityLaw` class from `python/pyexadis_base.py`. `MobilityLaw` is a wrapper class for computing nodal velocities.

#### Usage
A `MobilityLaw` module is instantiated by providing the global `state` parameters dictionary, the mobility type, and specific parameters to the mobility type:
```python
mobility = MobilityLaw(state=state, mobility_law='MobilityName', ...)
```
Available mobility types are:

* `mobility_law='SimpleGlide'`: simple linear mobility law where dislocation segments glide on the planes defined by their normal (nx,ny,nz). It does not require a crystal type and thus does not check whether the planes are crystallographically relevant. All Burgers vectors are treated equally. Is provided for testing purposes only. Specific mobility parameters:
    - `mob` (optional): isotropic mobility coefficient used for all dislocation character angles in 1/(Pa.s). Default: 1.0
    - `Medge` (optional): mobility coefficient used for the edge dislocation component in 1/(Pa.s). Must be used in pair with `Mscrew`. Mixed dislocation mobility coefficient will be linearly interpolated between edge and screw values. Default value: none.
    - `Mscrew` (optional): mobility coefficient used for the screw dislocation component in 1/(Pa.s). Must be used in pair with `Medge`. Default value: none.


* `mobility_law='FCC_0'`: generic linear mobility law for FCC crystals. Requires `crystal` to be set to `'fcc'` in the global parameters dictionary. Motion is only allowed on {111} planes. Junction nodes are restricted to move along their line direction (glide plane intersections). Mobility coefficients are linearly interpolated between edge and screw values. Specific mobility parameters:
    - `Medge` (required): mobility coefficient used for edge dislocation component in 1/(Pa.s)
    - `Mscrew` (required): mobility coefficient used for screw dislocation component in 1/(Pa.s)
    - `vmax` (optional): maximum dislocation velocity in m/s. A relativistic capping of the velocity is applied if specified. Default: none.


* `mobility_law='BCC_0b'`: generic linear mobility law for BCC crystals. Requires `crystal` to be set to `'bcc'` in the global parameters dictionary. The mobility matrix on glissile segments is constructed by summing a contribution from a pencil glide behavior of the screw segment component with a contribution from a planar behavior of the edge segment component. Specific mobility parameters:
    - `Medge` (required): mobility coefficient used for edge dislocation component in 1/(Pa.s)
    - `Mscrew` (required): mobility coefficient used for screw dislocation component in 1/(Pa.s)
    - `Mclimb` (required): mobility coefficient used for the climb component in 1/(Pa.s)
    - `vmax` (optional): maximum dislocation velocity in m/s. A relativistic capping of the velocity is applied if specified. Default: none.


Nodal velocities are computed using the `Mobility()` method
```python
state = mobility.Mobility(N, state)
```
The `Mobility()` method takes a dislocation network `DisNetManager` object and the `state` dictionary as arguments. The `state` dictionary must contain nodal force values stored in entries `nodeforces` and `nodeforcetags`, e.g. as would be stored after a call to a `CalForce.NodeForce()` method. The method returns an updated `state` dictionary with entry `nodevels` that contains the array of nodal velocities and corresponding node tags stored in entry `nodeveltags`:
```python
vels = state["nodevels"]
tags = state["nodeveltags"]
print(vels[0]) # print velocity of node tags[0]
```
The node tags are used to keep track of the nodes across modules, as different modules (e.g. in OpenDiS) may use different data structures and may shuffle the nodes ordering when being converted from one data structure to another via the `DisNetManager` object.

#### Properties
- `MobilityLaw.mobility_law`: name of the mobility law
- `MobilityLaw.mobility`: pointer to the ExaDiS mobility binding object

#### Methods
- `state = MobilityLaw.Mobility(DisNetManager, state)`: Computes nodal velocities for all nodes and place the result in the `state` dictionary.
- `v = MobilityLaw.OneNodeMobility(DisNetManager, state, tag, f)`: Compute the nodal velocity `v` on a single node specified by its `tag`. The nodal force must be provided with argument `f`.



### Time-Integrators
A time-integration module is declared using the `TimeIntegration` class from `python/pyexadis_base.py`. `TimeIntegration` is a wrapper class for advancing node positions during a time step.

#### Usage
A `TimeIntegration` module is instantiated by providing the global `state` dictionary, the time-integrator type, and specific parameters to the time-integrator type:
```python
timeint = TimeIntegration(state=state, integrator='IntegratorName', ...)
```
Most integrators requires `force` and `mobility` arguments pointing to force/mobility modules to be provided, as nodal forces and mobilities generally need to be evaluated internally during an integration step.

Available time-integrator types are:

* `integrator='EulerForward'`: Euler forward time-integrator that increments the node positions by multiplying all nodal velocities by a constant time step size. Specific integrator parameters:
    - `dt`: constant time step size in units of s. Default: 1e-8.
    

* `integrator='Trapezoid'`: Trapezoid time-integrator that uses an error-controlled scheme to select the time step size. Force and mobility modules must be provided as arguments. Multi time-stepping is enabled using parameter `multi`. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `state["rtol"]`: integrator tolerance in `burgmag` units.
    - `multi` (optional): multi time-stepping parameter that defines the number of sub time steps to execute the integrator for during a single global time step. Default value: 1 (no multi-step).
    - `state["nexdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


* `integrator='RKF'`: Runge–Kutta–Fehlberg (RKF45) time-integrator that uses a fifth-order accurate error-controlled scheme to select the time step size. Force and mobility modules must be provided as arguments. Multi time-stepping is enabled using parameter `multi`. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `state["rtol"]`: absolute tolerance in `burgmag` units.
    - `rtolrel` (optional): relative tolerance in `burgmag` units. Default: 0.1. 
    - `rtolth` (optional): threshold tolerance in `burgmag` units. Default: 1.0.
    - `multi` (optional): multi time-stepping parameter that defines the number of sub time steps to execute the integrator for during a single global time step. Default value: 1 (no multi-step).
    - `state["nexdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


* `integrator='Subcycling'`: Subcycling time integrator where force contributions are separated in various groups (5 groups) based on segment pair distances and integrated in turn in an asynchronous fashion. Force and mobility modules must be provided as arguments. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `rgroups`: list of group radii in increasing order. Must be a list of size 4.
    - `state["rtol"]`: absolute tolerance in `burgmag` units.
    - `rtolrel` (optional): relative tolerance in `burgmag` units. Default: 0.1. 
    - `rtolth` (optional): threshold tolerance in `burgmag` units. Default: 1.0.
    - `state["nexdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


#### Properties
- `TimeIntegration.integrator_type`: name of the integrator type
- `TimeIntegration.dt`: current time step size
- `TimeIntegration.integrator`: pointer to the ExaDiS integrator binding object

#### Methods
- `state = TimeIntegration.Update(DisNetManager, state)`: Perform a time-integration step and update the nodal positions accordingly.



### Collision
A collision module is declared using the `Collision` class from `python/pyexadis_base.py`. `Collision` is a wrapper class for handling intersections of dislocation segments.

#### Usage
A `Collision` module is instantiated by providing the global `state` dictionary, the collision mode, and specific parameters to the collision mode:
```python
collision = Collision(state=state, collision_mode='CollisionName', ...)
```
Available collision modes are:

* `collision_mode='Retroactive'`: Retroactive collision procedure in which all potential collisions that may have happened during a time-interval are tested for, even retroactively, as well as collisions within a proximity criterion. In order for the retroactive algorithm to be effective, entries `oldnodes_dict` containing a dictionary of nodes positions at the start of the time-interval must be provided in the `state` dictionary. Specific collision parameters:
    - `state["rann"]` (optional): annihilation/capture radius for proximity collisions. Default: 2*`state["rtol"]`.


#### Properties
- `Collision.collision_mode`: name of the collision mode
- `Collision.collision`: pointer to the ExaDiS collision binding object

#### Methods
- `state = Collision.HandleCol(DisNetManager, state)`: Handle collisions and modify the topology of the dislocation network accordingly.


### Topology
A topology module is declared using the `Topology` class from `python/pyexadis_base.py`. `Topology` is a wrapper class for handling topological events such as the split-multi-nodes procedure.

#### Usage
A `Topology` module is instantiated by providing the global `state` dictionary, the topology mode, `force` and `mobility` modules, and specific parameters to the topology mode:
```python
topology = Topology(state=state, topology_mode='TopologyName', force=calforce, mobility=mobilitylaw, ...)
```
Available topology modes are:

* `topology_mode='TopologySerial'`: Topology class that performs the split-multi-nodes procedure in a serial fashion on the host. One should prefer the much more efficient `topology_mode='TopologyParallel'` for production runs on GPU. Specific topology parameters:
    - `state["rann"]` (optional): annihilation/capture radius used to set the trial node splitting distance. Default: 2*`state["rtol"]`.
    - `state["split3node"]` (optional): flag to enable the splitting of 3-nodes (nodal cross-slip). Only operational for bcc crystals for now. Default: 1.
    - `splitMultiNodeAlpha` (optional): noise coefficient for the split-multi-nodes procedure. Default: 1e-3.


* `topology_mode='TopologyParallel'`: Topology class that performs the split-multi-nodes procedure in parallel fashion on the device (GPU). Specific topology parameters:
    - `state["rann"]` (optional): annihilation/capture radius used to set the trial node splitting distance. Default: 2*`state["rtol"]`.
    - `params["split3node"]` (optional): flag to enable the splitting of 3-nodes (nodal cross-slip). Only operational for bcc crystals for now. Default: 1.
    - `splitMultiNodeAlpha` (optional): noise coefficient for the split-multi-nodes procedure. Default: 1e-3.


#### Properties
- `Topology.topology_mode`: name of the topology mode
- `Topology.topology`: pointer to the ExaDiS topology binding object

#### Methods
- `state = Topology.Handle(DisNetManager, state)`: Handle topological changes (e.g. split-multi-nodes) and modify the topology of the dislocation network accordingly.


### Remesh
A remesh module is declared using the `Remesh` class from `python/pyexadis_base.py`. `Collision` is a wrapper class for performing remeshing of the dislocation segments.

#### Usage
A `Remesh` module is instantiated by providing the global `state` dictionary, the remesh rule, and specific parameters to the remesh rule:
```python
remesh = Remesh(state=state, remesh_rule='RemeshName', ...)
```
Available remesh rules are:

* `remesh_rule='LengthBased'`: Remeshing procedure based on minimum and maximum segment length parameters. Segments longer than `state["maxseg"]` will be bisected. End-nodes of segments shorter than `state["maxseg"]` will be merged when allowed. Specific remesh parameters:
    - `state["maxseg"]`: maximum segment discretization size length in units of `burgmag`.
    - `state["minseg"]`: minimum segment discretization size length in units of `burgmag`.


#### Properties
- `Remesh.remesh_rule`: name of the remesh rule
- `Remesh.remesh`: pointer to the ExaDiS remesh binding object

#### Methods
- `state = Remesh.Remesh(DisNetManager, state)`: Remesh the dislocation network and modify its topology accordingly.


### Simulation driver
A simulation driver is declared using the `SimulateNetwork` or `SimulateNetworkPerf` class from `python/pyexadis_base.py`. `SimulateNetwork` and `SimulateNetworkPerf` are base classes to drive a traditional DDD simulation that invoke the different stages of the simulation cycle. It must be defined by passing all base modules to be used for the simulation (`CalForce`, `MobilityLaw`, etc.) plus some additional parameters defining the simulation control, e.g. strain rate, number of steps, etc. For instance:
```python
# Default driver optimized for inter-operability and prototyping
sim = SimulateNetwork(
    calforce=calforce, mobility=mobility, timeint=timeint, 
    collision=collision, topology=topology, remesh=remesh, vis=vis,
    loading_mode=loading_mode, erate=erate, edir=edir,
    max_step=max_step, burgmag=state["burgmag"], state=state,
    print_freq=print_freq, plot_freq=plot_freq, plot_pause_seconds=0.0001,
    write_freq=write_freq, write_dir=output_dir
)
# ExaDiS driver optimized for performance (production runs)
sim = SimulateNetworkPerf(
    calforce=calforce, mobility=mobility, timeint=timeint, 
    collision=collision, topology=topology, remesh=remesh, vis=vis,
    loading_mode=loading_mode, erate=erate, edir=edir,
    max_strain=max_strain, burgmag=state["burgmag"], state=state,
    print_freq=print_freq, write_freq=write_freq, write_dir=output_dir,
    restart=restart
)
```
Once the driver object is instantiated, the simulation is run with method `run()` passing the initial dislocation network `DisNetManager` object and `state` dictionary as its arguments
```python
sim.run(N, state)
```

When driving a pure ExaDiS simulation (i.e. all modules used are ExaDiS modules) and high performance is desired, it is advised to use the `SimulateNetworkPerf` driver. The `SimulateNetworkPerf` performance driver is a direct wrapper around the C++ ExaDiS driver and thus removes data movement through the python interface, in contrast with `SimulateNetwork`, which is a more flexible but slightly less efficient class intended to work with arbitrary OpenDiS modules. `SimulateNetworkPerf` cannot be used with external modules (must use ExaDiS modules only), and it does not allow for interactive data visualization or manipulation. This driver is recommended for production runs of large systems, e.g. for strain-hardening simulations.

`SimulateNetworkPerf` also provides a restart mechanism to resume a simulation from a previous run, using the `restart` argument (e.g. see example `examples/22_fcc_Cu_15um_1e3/example_fcc_Cu_15um_1e3.py` for how to use the restart feature).

Similar to the backend C++ driver, `SimulateNetworkPerf` provides several ways for setting the stopping criterion for a simulation via the following input arguments:
- `num_steps`: number of simulation steps to reach for the current loading instance
- `max_steps`: total number of simulation steps to reach across all loading instances (e.g. including all restarts)
- `max_strain`: maximum (total) strain to reach for the simulation
- `max_time`: maximum physical simulation time (cumulative time step sizes) in s to reach for the simulation
- `max_walltime`: maximum wall-clock (execution) time in s to reach for the simulation
