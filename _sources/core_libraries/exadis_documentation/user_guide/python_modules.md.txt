## Python modules: pyexadis
This section documents the various ExaDiS modules available through the python interface, `pyexadis`. These modules are bindings to the backend C++ modules implemented in ExaDiS. For documentation about the backend C++ modules please see the [Developer guide](../developer_guide/index) section of the documentation.

```{hint}
Before going through this section, make sure you have read the [tutorial section](tutorials.md#setting-up-a-python-driven-simulation) on how to set up a script/simulation using the `pyexadis` interface to ExaDiS.
```

### Dislocation network

#### `ExaDisNet`
Dislocation networks are defined as instances of the `ExaDisNet` class from `python/pyexadis_base.py`. `ExaDisNet` is a wrapper class around the internal data structure representation of the dislocation network in ExaDiS, which internally handles memory movements between the different execution spaces (e.g. CPU to GPU).

See section [ExaDisNet Class](../../../code_structure/data_structure/exadisnet_class.md) for more information.

An `ExaDisNet` network can be instantiated in several ways, see section [Creating initial dislocation configurations](../../../tutorials/initial_configuration).

````{Important}
When used in ExaDiS or OpenDiS modules, the dislocation network defined in a `ExaDisNet` object must first be wrapped into a `DisNetManager` object before it can be used within modules (see next section), e.g.
```python
G = ExaDisNet(...)
N = DisNetManager(G)
```
````

#### `DisNetManager`

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

See section [DisNetManager Class](../../../code_structure/data_structure/disnetmanager_class.md) to get more information.


### Forces

#### `CalForce` modules
A force module is declared using the `CalForce` class from `python/pyexadis_base.py`. `CalForce` is a wrapper class for computing forces at dislocation nodes.

##### Usage
A `CalForce` module is instantiated by providing the global `state` dictionary, the force mode, and specific parameters to the force mode:
```python
calforce = CalForce(state=state, force_mode='ForceName', ...)
```
Available force modes are:

* `force_mode='DDD_FFT_MODEL'`: DDD-FFT model where forces are computed by summing (i) external PK forces from the applied stress, (ii) core forces, (iii) short-range elastic interactions computed with the non-singular model, and (iv) long-range forces computed using FFT. Specific mode parameters are:
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


* `force_mode='CUTOFF_MODEL'`: Cutoff model where forces are computed by summing (i) external PK forces from the applied stress, (ii) core forces, and (iii) short-range elastic interactions up to a given cutoff distance between segment pairs, computed with the non-singular model. Specific mode parameters are:
    - `cutoff` (required): cutoff for segment pair interaction distance.
    - `Ec` (optional): core energy factor in Pa. If not provided or negative, it defaults to the standard value Ec=mu/4/pi*ln(a/0.1).
    - `Ecore_junc_fact` (optional): ratio of the core energy of junction to native Burgers vectors. Default: 1.0.


* `force_mode='LINE_TENSION_MODEL'` or `force_mode='LineTension'`: line-tension model where forces are computed from (i) the external PK force from the applied stress and (ii) the line energy (core force), and no pairwise elastic interaction is accounted for. Specific mode parameters are:
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
print(forces[0]) # print force at node "tags[0]"
```
The node tags are used to keep track of the nodes across modules, as different modules (e.g. in OpenDiS) may use different data structures and may shuffle the nodes ordering when being converted from one data structure to another via the `DisNetManager` object.

##### Properties
- `CalForce.force_mode`: name of the force mode
- `CalForce.params`: ExaDiS global parameters used to setup the force model
- `CalForce.force`: pointer to the ExaDiS force binding object

##### Methods
- `state = CalForce.NodeForce(N: DisNetManager, state: dict, pre_compute=True)`: Computes nodal forces for all nodes and places the result in the `state` dictionary. By default, option `pre_compute=True`, meaning that a call to the `PreCompute()` operation will be made before evaluating the forces.
- `state = CalForce.PreCompute(N: DisNetManager, state: dict)`: Call to the force pre-compute operation, e.g. as some force modules (e.g. FFT) requires certain values to be pre-computed and up-to-date before forces can be evaluated efficiently. In general, there is no need to call `PreCompute()` before a call to the global `NodeForce()` method, as by default the pre-compute step is called internally in `NodeForce()`. However, a call to `PreCompute()` may be required before calling `OneNodeForce()` (e.g. if the network has been updated since the last call to `NodeForce()`). 
- `f = CalForce.OneNodeForce(N: DisNetManager, state: state, tag, update_state=True)`: Compute the nodal force `f` on a single node specified by its `tag`. The `PreCompute()` operation may need to be called before calling `OneNodeForce()`.


#### Force calculation wrappers
In addition to `CalForce` modules, `pyexadis` provides wrappers to some internal ExaDiS force calculation functions, which can be called from python with `pyexadis.<function>(...)`. The following functions are currently available:
- `f = pyexadis.compute_force_n2(G: ExaDisNet, mu, nu, a)`: compute elastic interaction forces using a brute-force N^2 segment pair summation.
- `f = pyexadis.compute_force_cutoff(G: ExaDisNet, mu, nu, a, cutoff, maxseg)`: compute elastic forces including interactions up to a given cutoff distance between segment pairs. Argument `maxseg` is needed to ensure that no segment pair is missed when building the neighbor list.
- `f = pyexadis.compute_force_segseglist(G: ExaDisNet, mu, nu, a, segseglist)`: compute elastic forces provided a user-defined list of segment pairs. Argument `segseglist` must be an array of size=(Npair,2) indicating the indices of each pair of segments to consider in the calculation. Each pair must be unique, e.g. if pair (2,3) is included, then pair (3,2) must not be included or its contribution will be double-counted. Self-interactions, e.g. pair (2,2), must not be included.



### Mobility laws
A mobility module is declared using the `MobilityLaw` class from `python/pyexadis_base.py`. `MobilityLaw` is a wrapper class for computing nodal velocities.

#### Usage
A `MobilityLaw` module is instantiated by providing the global `state` parameters dictionary, the mobility type, and specific parameters to the mobility type:
```python
mobility = MobilityLaw(state=state, mobility_law='MobilityName', ...)
```
Available mobility types are:

* `mobility_law='SimpleGlide'`: simple linear mobility law where dislocation segments glide on the planes defined by their normal (nx,ny,nz). It does not require a crystal type and thus does not check whether the planes are crystallographically relevant. All Burgers vectors are treated equally. It is provided for testing purposes only. Specific mobility parameters:
    - `mob` (optional): isotropic mobility coefficient used for all dislocation character angles in 1/(Pa.s). Default: 1.0
    - `Medge` (optional): mobility coefficient used for the edge dislocation component in 1/(Pa.s). Must be used in pair with `Mscrew`. Mixed dislocation mobility coefficient will be linearly interpolated between edge and screw values. Default value: none.
    - `Mscrew` (optional): mobility coefficient used for the screw dislocation component in 1/(Pa.s). Must be used in pair with `Medge`. Default value: none.


* `mobility_law='FCC_0'`: generic planar, linear mobility law for FCC crystals. Requires `crystal` to be set to `'fcc'` in the global parameters dictionary. Motion is only strictly allowed on {111} planes. Junction nodes are restricted to move along their line direction (glide plane intersections). By default, glide constraints are enforced by systematically projecting node velocities onto their glide planes. Mobility coefficients are linearly interpolated between edge and screw values. Specific mobility parameters:
    - `Medge` (required): mobility coefficient used for edge dislocation component in 1/(Pa.s)
    - `Mscrew` (required): mobility coefficient used for screw dislocation component in 1/(Pa.s)
    - `vmax` (optional): maximum dislocation velocity in m/s. A relativistic capping of the velocity is applied if specified. Default: none.


* `mobility_law='FCC_0_FRIC'`: same mobility law as `FCC_0` but with the possibility to define friction stresses and spatially varying mobility/friction fields. Specific mobility parameters:
    - `Medge` (required): mobility coefficient used for edge dislocation component in 1/(Pa.s)
    - `Mscrew` (required): mobility coefficient used for screw dislocation component in 1/(Pa.s)
    - `vmax` (optional): maximum dislocation velocity in m/s. A relativistic capping of the velocity is applied if specified. Default: none.
    - `Fedge` (optional): friction stress used for edge dislocation component in Pa. Default: 0.0.
    - `Fscrew` (optional): friction stress used for screw dislocation component in Pa. Default: 0.0.
    - `mobility_field` (optional): path to a file defining field values by which the mobility coefficient will be scaled throughout the simulation volume. File format must be the following: the 3 first rows specify the grid size in the 3 directions (Nx, Ny, Nz), and the remaining rows are the grid values (Nx x Ny x Nz values) in C-like index order (last axis index changing the fastest). If one of the grid dimensions is 0 or 1, then the field will be treated as 2D and constant along this dimension.
    - `friction_field` (optional): path to a file defining field values by which the friction stress will be scaled throughout the simulation volume. File format is the same as for the `mobility_field` parameter.


* `mobility_law='FCC_0B'`: generic non-planar, linear mobility law for FCC crystals . Requires `crystal` to be set to `'fcc'` in the global parameters dictionary. In contrast to `FCC_0`, it is modeled after the `BCC_0B` mobility and allows for a climb component of the velocity. To allow for climb, option `state["enforce_glide_planes"] = 0` must be specified, otherwise node velocities will be projected back to glide planes. Specific mobility parameters:
    - `Medge` (required): mobility coefficient used for edge dislocation component in 1/(Pa.s)
    - `Mscrew` (required): mobility coefficient used for screw dislocation component in 1/(Pa.s)
    - `Mclimb` (required): mobility coefficient used for the climb component in 1/(Pa.s)
    - `vmax` (optional): maximum dislocation velocity in m/s. A relativistic capping of the velocity is applied if specified. Default: none.


* `mobility_law='BCC_0B'`: generic linear mobility law for BCC crystals. Requires `crystal` to be set to `'bcc'` in the global parameters dictionary. The mobility matrix on glissile segments is constructed by summing a contribution from a pencil glide behavior of the screw segment component with a contribution from a planar behavior of the edge segment component. By default, no explicit glide planes are defined (non-planar mobility). To use and enforce glide planes (planar mobility), option `state["use_glide_planes"] = 1` and `state["enforce_glide_planes"] = 1` must be specified. Specific mobility parameters:
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
print(vels[0]) # print velocity of node "tags[0]"
```
The node tags are used to keep track of the nodes across modules, as different modules (e.g. in OpenDiS) may use different data structures and may shuffle the nodes ordering when being converted from one data structure to another via the `DisNetManager` object.

#### Properties
- `MobilityLaw.mobility_law`: name of the mobility law
- `MobilityLaw.mobility`: pointer to the ExaDiS mobility binding object

#### Methods
- `state = MobilityLaw.Mobility(N: DisNetManager, state: dict)`: Computes nodal velocities for all nodes and places the result in the `state` dictionary.
- `v = MobilityLaw.OneNodeMobility(N: DisNetManager, state: dict, tag, f)`: Compute the nodal velocity `v` on a single node specified by its `tag`. The nodal force must be provided with argument `f`.



### Time-Integrators
A time-integration module is declared using the `TimeIntegration` class from `python/pyexadis_base.py`. `TimeIntegration` is a wrapper class for advancing node positions during a time step.

#### Usage
A `TimeIntegration` module is instantiated by providing the global `state` dictionary, the time-integrator type, and specific parameters to the time-integrator type:
```python
timeint = TimeIntegration(state=state, integrator='IntegratorName', ...)
```
Most integrators require `force` and `mobility` arguments pointing to force/mobility modules to be provided, as nodal forces and mobilities generally need to be evaluated internally during an integration step.

Available time-integrator types are:

* `integrator='EulerForward'`: Euler forward time-integrator that increments the node positions by multiplying all nodal velocities by a constant time step size. Specific integrator parameters:
    - `dt`: constant time step size in units of s. Default: 1e-8.
    

* `integrator='Trapezoid'`: Trapezoid time-integrator that uses an error-controlled scheme to select the time step size. Force and mobility modules must be provided as arguments. Multi time-stepping is enabled using parameter `multi`. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `state["rtol"]`: integrator tolerance in `burgmag` units.
    - `multi` (optional): multi time-stepping parameter that defines the number of sub time steps to execute the integrator for during a single global time step. Default value: 1 (no multi-step).
    - `state["nextdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


* `integrator='RKF'`: Runge–Kutta–Fehlberg (RKF45) time-integrator that uses a fifth-order accurate error-controlled scheme to select the time step size. Force and mobility modules must be provided as arguments. Multi time-stepping is enabled using parameter `multi`. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `state["rtol"]`: absolute tolerance in `burgmag` units.
    - `rtolrel` (optional): relative tolerance in `burgmag` units. Default: 0.1. 
    - `rtolth` (optional): threshold tolerance in `burgmag` units. Default: 1.0.
    - `multi` (optional): multi time-stepping parameter that defines the number of sub time steps to execute the integrator for during a single global time step. Default value: 1 (no multi-step).
    - `state["nextdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


* `integrator='Subcycling'`: Subcycling time integrator where force contributions are separated in various groups (5 groups) based on segment pair distances and integrated in turn in an asynchronous fashion. The force model must be `SUBCYCLING_MODEL`. Specific integrator parameters:
    - `force`: `CalForce` module to be used for nodal force calculations. Must be a `CalForce` `SUBCYCLING_MODEL` model.
    - `mobility`: `MobilityLaw` module to be used for nodal velocities calculations.
    - `rgroups`: list of group radii in increasing order. Must be a list of size 4.
    - `state["rtol"]`: absolute tolerance in `burgmag` units.
    - `rtolrel` (optional): relative tolerance in `burgmag` units. Default: 0.1. 
    - `rtolth` (optional): threshold tolerance in `burgmag` units. Default: 1.0.
    - `state["nextdt"]` (optional): initial trial value for the time step size in s. Default: 1e-12.
    - `state["maxdt"]` (optional): maximum value for the time step size in s. Default: 1e-7.


#### Properties
- `TimeIntegration.integrator_type`: name of the integrator type
- `TimeIntegration.dt`: current time step size
- `TimeIntegration.integrator`: pointer to the ExaDiS integrator binding object

#### Methods
- `state = TimeIntegration.Update(N: DisNetManager, state: dict)`: Perform a time-integration step and update the nodal positions accordingly.



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


* `collision_mode='Proximity'`: Will default to the retroactive collision procedure above. Specific collision parameters:
    - `state["rann"]` (optional): annihilation/capture radius for proximity collisions. Default: 2*`state["rtol"]`.


* `collision_mode='None'`: No collision procedure.


#### Properties
- `Collision.collision_mode`: name of the collision mode
- `Collision.collision`: pointer to the ExaDiS collision binding object

#### Methods
- `state = Collision.HandleCol(N: DisNetManager, state: dict)`: Handle collisions and modify the topology of the dislocation network accordingly.


### Topology
A topology module is declared using the `Topology` class from `python/pyexadis_base.py`. `Topology` is a wrapper class for handling topological events such as the split-multi-nodes procedure.

#### Usage
A `Topology` module is instantiated by providing the global `state` dictionary, the topology mode, `force` and `mobility` modules, and specific parameters to the topology mode:
```python
topology = Topology(state=state, topology_mode='TopologyName', force=calforce, mobility=mobilitylaw, ...)
```
Available topology modes are:

* `topology_mode='TopologySerial'`: Topology class that performs the split-multi-nodes procedure in a serial fashion on the host. Note: one should prefer the much more efficient `'TopologyParallel'` mode for production runs on GPU. Specific topology parameters:
    - `state["rann"]` (optional): annihilation/capture radius used to set the trial node splitting distance. Default: 2*`state["rtol"]`.
    - `state["split3node"]` (optional): flag to enable the splitting of 3-nodes (nodal cross-slip). Only operational for bcc crystals for now. Default: 1.
    - `splitMultiNodeAlpha` (optional): noise coefficient for the split-multi-nodes procedure. Default: 1e-3.


* `topology_mode='TopologyParallel'`: Topology class that performs the split-multi-nodes procedure in parallel fashion on the device (GPU). Specific topology parameters:
    - `state["rann"]` (optional): annihilation/capture radius used to set the trial node splitting distance. Default: 2*`state["rtol"]`.
    - `params["split3node"]` (optional): flag to enable the splitting of 3-nodes (nodal cross-slip). Only operational for bcc crystals for now. Default: 1.
    - `splitMultiNodeAlpha` (optional): noise coefficient for the split-multi-nodes procedure. Default: 1e-3.


* `topology_mode='None'`: No topology procedure.


#### Properties
- `Topology.topology_mode`: name of the topology mode
- `Topology.topology`: pointer to the ExaDiS topology binding object

#### Methods
- `state = Topology.Handle(N: DisNetManager, state: dict)`: Handle topological changes (e.g. split-multi-nodes) and modify the topology of the dislocation network accordingly.


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


* `remesh_rule='None'`: No remesh procedure.


#### Properties
- `Remesh.remesh_rule`: name of the remesh rule
- `Remesh.remesh`: pointer to the ExaDiS remesh binding object

#### Methods
- `state = Remesh.Remesh(N: DisNetManager, state: dict)`: Remesh the dislocation network and modify its topology accordingly.



### Cross-slip
A cross-slip module is declared using the `CrossSlip` class from `python/pyexadis_base.py`. `CrossSlip` is a wrapper class for handling cross-slip of dislocation segments.

#### Usage
A `CrossSlip` module is instantiated by providing the global `state` dictionary, the cross-slip mode, and specific parameters to the cross-slip mode:
```python
cross_slip = CrossSlip(state=state, cross_slip_mode='CrossSlipName', ...)
```
Available cross-slip modes are:

* `cross_slip_mode='ForceBasedSerial'`: Cross-slip procedure where near screw segments are cross-slipped to a new plane when the resolved force on the new plane exceeds that of the original plane by some threshold. This module is executed in serial fashion on the host. Specific cross-slip parameters:
    - `force`: `CalForce` force module to evaluate the force-based cross-slip criterion.


* `cross_slip_mode='ForceBasedParallel'`: Same cross-slip procedure as `ForceBasedSerial` but where the module is executed in parallel fashion on the device (GPU). Specific cross-slip parameters:
    - `force`: `CalForce` force module to evaluate the force-based cross-slip criterion.


* `cross_slip_mode='None'`: No cross-slip procedure.


#### Properties
- `CrossSlip.cross_slip_mode`: name of the cross-slip mode
- `CrossSlip.cross_slip`: pointer to the ExaDiS cross-slip binding object

#### Methods
- `state = CrossSlip.Handle(N: DisNetManager, state: dict)`: Handle cross-slip operations and update the dislocation network accordingly.



### Simulation driver
A simulation driver is declared using the `SimulateNetwork` or `SimulateNetworkPerf` class from `python/pyexadis_base.py`. `SimulateNetwork` and `SimulateNetworkPerf` are base classes to drive a traditional DDD simulation that invoke the different stages of the simulation cycle. It must be defined by passing all base modules to be used for the simulation (`CalForce`, `MobilityLaw`, etc.) plus some additional parameters defining the simulation control, e.g. strain rate, number of steps, etc. For instance:
```python
# Default driver optimized for inter-operability and prototyping
sim = SimulateNetwork(
    calforce=calforce, mobility=mobility, timeint=timeint, collision=collision, 
    topology=topology, remesh=remesh, cross_slip=cross_slip, vis=vis,
    loading_mode=loading_mode, erate=erate, edir=edir,
    max_step=max_step, burgmag=state["burgmag"], state=state,
    print_freq=print_freq, plot_freq=plot_freq, plot_pause_seconds=0.0001,
    write_freq=write_freq, write_dir=output_dir
)

# ExaDiS driver optimized for performance (e.g. production runs on GPU)
sim = SimulateNetworkPerf(
    calforce=calforce, mobility=mobility, timeint=timeint, collision=collision, 
    topology=topology, remesh=remesh, cross_slip=cross_slip, vis=vis,
    loading_mode=loading_mode, erate=erate, edir=edir,
    max_strain=max_strain, burgmag=state["burgmag"], state=state,
    print_freq=print_freq, write_freq=write_freq, write_dir=output_dir,
    out_props=['step', 'strain', 'stress', 'density'],
    restart=restart
)
```
Once the driver object is instantiated, the simulation is run with method `run()` passing the initial dislocation network `DisNetManager` object and `state` dictionary as its arguments
```python
sim.run(N, state)
```

#### `SimulateNetwork`

`SimulateNetwork` is a flexible, python-based driver that allows for the use and inter-operability of OpenDiS modules (e.g. `pyexadis` modules can be mixed with other OpenDiS modules), and for interactive data visualization and manipulation.

Note that when using `SimulateNetwork` however, ExaDiS simulation modules may suffer some performance hit (compared to a native C++ driven simulation) due to data flow through the python pipeline. This may be particularly apparent when running on GPU. This overhead is generally reasonable (e.g. 20%) and allows for a great flexibility of the code. If maximum performance is desired, it is advised to use the `SimulateNetworkPerf` driver instead.

Currently, the two following stopping criteria are available:
- `max_step`: maximum simulation step number to reach (across all loading instances)
- `num_steps`: number of simulation steps to perform for the current loading instance


#### `SimulateNetworkPerf`

`SimulateNetworkPerf` is a driver intended for maximum performance: when driving a pure ExaDiS simulation (i.e. all modules used are `pyexadis` modules) and high performance is desired, it is advised to use the `SimulateNetworkPerf` driver. `SimulateNetworkPerf` is a direct wrapper around the C++ ExaDiS driver and thus removes data movement through the python interface, in contrast to `SimulateNetwork`, which is a more flexible but slightly less efficient class intended to work with arbitrary OpenDiS modules. Conversely, `SimulateNetworkPerf` cannot be used with external modules (must use `pyexadis` modules only), and it does not allow for interactive data visualization or manipulation. This driver is recommended for production runs of large systems, e.g. for strain-hardening simulations.

`SimulateNetworkPerf` also provides a restart mechanism to resume a simulation from a previous run, using the `restart` argument (e.g. see example `examples/22_fcc_Cu_15um_1e3/example_fcc_Cu_15um_1e3.py` for how to use the restart feature).

Similar to the backend C++ driver, `SimulateNetworkPerf` provides several ways for setting the stopping criterion for a simulation via the following input arguments:
- `num_steps`: number of simulation steps to perform for the current loading instance
- `max_step`: maximum simulation step number to reach across all loading instances (e.g. including all restarts)
- `max_strain`: maximum (total) strain to reach for the simulation
- `max_time`: maximum physical simulation time (cumulative time step sizes) in s to reach for the simulation
- `max_walltime`: maximum wall-clock (execution) time in s to reach for the simulation

`SimulateNetworkPerf` provides an option to specify the simulation properties to be written in the output property file `{write_dir}/stress_strain_dens.dat` at frequency `propfreq`. This can be specified by passing a list of properties to `SimulateNetworkPerf` using argument `out_props=[...]`. Available properties are:
- `step`: simulation step number
- `strain`: uniaxial strain
- `stress`: uniaxial stress
- `density`: total dislocation density
- `Nnodes`: number of dislocation nodes
- `Nsegs`: number of dislocation segments
- `dt`: current timestep size
- `time`: simulation time
- `walltime`: wall-clock time
- `edir`: current loading direction
- `Rorient`: current crystal orientation matrix
- `allstress`: all components of the applied stress tensor

If not specified, `out_props=['step', 'strain', 'stress', 'density']` by default.


### Utility functions

A set of utility functions are provided in file `python/pyexadis_utils.py`. Available functions are:

#### Network generation

* `insert_frank_read_src(cell, nodes, segs, burg, plane, length, center, theta=0.0, linedir=None, numnodes=10)`: Insert a Frank-Read source into the list of nodes and segments.
    - `cell`: network cell object
    - `nodes`: list of nodes
    - `segs`: list of segments
    - `burg`: Burgers vector of the source
    - `plane`: habit plane normal of the source
    - `length`: length of the source
    - `center`: center position of the source
    - `theta`: character angle of the source in degrees
    - `linedir`: line direction of the source
    - `numnodes`: number of discretization nodes for the source
    
* `insert_infinite_line(cell, nodes, segs, burg, plane, origin, theta=0.0, linedir=None, maxseg=-1, trial=False)`: Insert an infinite line into the list of nodes and segments
    - `cell`: network cell object
    - `nodes`: list of nodes
    - `segs`: list of segments
    - `burg`: Burgers vector of the line
    - `plane`: habit plane normal of the source
    - `origin`: origin position of the line
    - `theta`: character angle of the line in degrees
    - `linedir`: line direction
    - `maxseg`: maximum discretization length of the line
    - `trial`: do a trial insertion only (to test if insertion is possible)
    
* `generate_line_config(crystal, Lbox, num_lines, theta=None, maxseg=-1, Rorient=None, seed=-1, verbose=True)`: Generate a configuration made of straight, infinite dislocation lines into an `ExaDisNet` object
    - `crystal`: crystal structure
    - `Lbox`: box size or network cell object
    - `num_lines`: number of dislocation lines
    - `theta`: list of possible character angles in degrees
    - `maxseg`: maximum discretization length of the lines
    - `Rorient`: crystal orientation matrix
    - `seed`: seed for random number generation
    - `verbose`: print information

* `generate_prismatic_config(crystal, Lbox, num_loops, radius, maxseg=-1, Rorient=None, seed=-1, uniform=False)`: Generate a configuration made of prismatic dislocation loops into an `ExaDisNet` object
    - `crystal`: crystal structure
    - `Lbox`: box size or network cell object
    - `num_loops`: number of dislocation loops
    - `radius`: loop radius or [min_radius, max_radius] range
    - `maxseg`: maximum discretization length of the lines
    - `Rorient`: crystal orientation matrix
    - `seed`: seed for random number generation
    - `uniform`: make the spatial loop distribution close to uniform

#### Network manipulation

* `combine_networks(Nlist)`: Combine several DisNetManager into a single network
    - `Nlist`: list of `DisNetManager` objects
    
* `extract_segments(N: DisNetManager, seglist)`: Return a new network that contains a subset of segments from the input network
    - `N`: `DisNetManager` object
    - `seglist`: list of segment indices to extract

* `delete_segments(N: DisNetManager, seglist)`: Return a new network in which segments have been deleted from the input network
    - `N`: `DisNetManager` object
    - `seglist`: list of segment indices to delete

#### Network properties

* `get_segments_length(N: DisNetManager)`: Return the list of dislocation segment lenghts of the network
    - `N`: `DisNetManager` object

* `dislocation_density(N: DisNetManager, burgmag: float)`: Return the dislocation density of the network
    - `N`: `DisNetManager` object
    - `burgmag`: Burgers vector scale to return the density in units of 1/m^2.
    
* `dislocation_charge(N: DisNetManager)`: Return the dislocation charge (net Nye's tensor) of the network
    - `N`: `DisNetManager` object

#### Input / output

* `read_paradis(datafile: str)`: Read dislocation network in ParaDiS format into a `DisNetManager` object
    - `datafile`: path of the data file to read

* `write_data(N: DisNetManager, datafile: str)`: Write dislocation network in ParaDiS format
    - `N`: `DisNetManager` object
    - `datafile`: path of the data file to write

* `write_vtk(N: DisNetManager, vtkfile: str, segprops={}, pbc_wrap=True)`: Write dislocation network in vtk format
    - `N`: `DisNetManager` object
    - `vtkfile`: path of the vtk file to write
    - `segprops`: dictionary of additional segments property fields ("name": values) to write
    - `pbc_wrap`: wrap dislocation segments into primary volume
