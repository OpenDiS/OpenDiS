### Setting up a simulation

This section details the various steps to set up and run a DDD simulation with OpenDiS using a python script.

#### 1. Importing modules and core libraries

First, import the OpenDiS framework modules, e.g.
```python
import os, sys

opendis_path = '/path/to/OpenDiS/python/'
if not opendis_path in sys.path: sys.path.append(opendis_path)
from framework.disnet_manager import DisNetManager
```

Then, for importing `pydis` modules, e.g.:
```python
# Importing pydis
pydis_paths = ['/path/to/OpenDiS/core/pydis/python', '/path/to/OpenDiS/lib']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
from pydis import DisNet, Cell, DisNode
from pydis import CalForce, MobilityLaw, TimeIntegration, Topology, Collision, Remesh
```

For importing `pyexadis` modules, e.g.:
```python
# Importing pyexadis
pyexadis_path = '/path/to/OpenDiS/core/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
try:
    import pyexadis
    from pyexadis_base import ExaDisNet, SimulateNetwork, VisualizeNetwork
    from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
except ImportError:
    raise ImportError('Cannot import pyexadis')
```

````{Warning}
Standard modules (e.g. `CalForce`, `MobilityLaw`, etc.) have the same names in core libraries, so watch out for name conflicts when importing them. For instance, doing
```python
from pydis import CalForce
from pyexadis_base import CalForce
```
will result in only the `CalForce` module from `pyexadis` to be imported. If needing to import both modules simultaneously, use aliases, e.g.
```python
from pydis import CalForce as CalForce_pydis
from pyexadis_base import CalForce as CalForce_pyexadis
```
````

#### 2. Initializing pyexadis

When using `pyexadis` modules, `pyexadis` must be explicitly initialized before it can be used. For this, all the simulation script must be inserted in between calls to the `pyexadis.initialize()` and `pyexadis.finalize()` functions:
```python
# Importing stuff

pyexadis.initialize()

# Code to setup and run the simulation here...

pyexadis.finalize()
```
This is to make sure that Kokkos is initialized and terminated properly on the chosen execution spaces. The `pyexadis.initialize()` can take the following arguments:
* `num_threads`: specifies the number of OpenMP threads to be used. For instance, use `pyexadis.initialize(num_threads=8)` to run a simulation using 8 threads. If not specified (default), ExaDiS will use the maximum number of threads available on the machine or defined from the `OMP_NUM_THREADS` environment variable if set.
* `device_id`: specifies the id of the device (GPU) to be used. If not specified (default), ExaDiS will select the first entry in the list of available devices.

If `pyexadis.initialize()` is not called before execution of a `pyexadis` function, the code will throw an error of the following type:
```
Constructing View and initializing data with uninitialized execution space
```

#### 3. Defining the global state dictionary

The global `state` dictionary is used to hold a collection of variables and arrays defining the state of a simulation, see [State dictionary section](../code_structure/data_structure/state_dictionary.md). It is propagated through the modules during a simulation: modules operate and communicate with each other by reading and writing entries to the `state` dictionary.

The `state` dictionary must be initialized with the global simulation parameters, e.g.

```python
state = {
    "crystal": 'fcc',
    "burgmag": 2.55e-10,
    "mu": 54.6e9,
    "nu": 0.324,
    "a": 6.0,
    "maxseg": 2000.0,
    "minseg": 300.0,
    "rtol": 10.0,
    "rann": 10.0,
    "nextdt": 1e-10,
    "maxdt": 1e-9,
}
```

Global simulation parameters are those that must be consistent across the different modules used in a simulation. Some parameters are required and other optional, although it is good practice to specify a maximum number of these. Here is a list of typical parameters:
* `crystal` (optional): Crystal structure. Optional but most mobility laws will require it to be provided. Supported options are `'fcc'` and `'bcc'` for now. Default: none.
* `burgmag` (required): Burgers vector magnitude in m. All lengths used in a simulation will be scaled by this value.
* `mu` (required): Shear modulus in Pa.
* `nu` (required): Poisson's ratio.
* `a` (required): Dislocation non-singular core radius in `burgmag` units.
* `maxseg` (required): Maximum discretization segment length in `burgmag` units.
* `minseg` (required): Minimum discretization segment length in `burgmag` units.
* `Rorient` (optional): Crystal orientation matrix. Default: Identity matrix.
* `rtol` (optional): Integration tolerance in `burgmag` units. Default: 0.25*`a`.
* `rann` (optional): Annihilation radius in `burgmag` units. Default: 2*`rtol`.
* `maxdt` (optional): Maximum timestep size. Default: 1e-7 s.
* `nextdt` (optional): Initial timestep size. Default: 1e-12 s.

For ExaDiS modules, the following global parameters can be specified:
* `split3node` (optional): Toggle to enable splitting of 3-nodes (only for BCC crystals). Default: 1.
* `use_glide_planes` (optional): Toggle to enable the use of glide planes. Default: 1 for `'fcc'`, 0 for `'bcc'`.
* `enforce_glide_planes` (optional): Toggle to enforce glide constraints by systematically projecting node velocities onto glide planes. Default: value of `use_glide_planes`.

By specification, the `state` dictionary needs to be later passed as an argument when initializing simulation modules.


#### 4. Creating the initial dislocation configuration

Ways to create initial dislocation configurations are detailed in the [Creating initial dislocation configurations](initial_configuration.md) section. For instance, this can be done by reading a dislocation configuration from a ParaDiS format file:
```python
G = ExaDisNet()
G.read_paradis('my_config.data')
```

The network object then needs to be wrapped into a `DisNetManager` object:
```python
N = DisNetManager(G)
```
This step is required to ensure inter-operability between the different core implementations and data structures.


#### 5. Defining simulation modules

Next we define the simulation modules to be used in the simulation. A full DDD simulation typically uses the following modules corresponding to the different elementary stages of the DDD simulation cycle:
* `CalForce`: module that computes forces at dislocation nodes
* `MobilityLaw`: module that computes velocities of dislocation nodes
* `TimeIntegration`: module that performs time-integration of the system
* `Collision`: module that implements dislocation collisions
* `Topology`: module that implements dislocation topological operations
* `CrossSlip`: module that performs cross-slip events
* `Remesh`: module that performs remeshing of the dislocation lines

Each module must be instantiated by passing the global state dictionary `state` and module specific arguments, e.g.
```python
mobility = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
```
* Available modules and parameters from `pyexadis` are documented in the [ExaDiS Python Modules](../core_libraries/exadis_documentation/user_guide/python_modules) section.
* Available modules and parameters from `pydis` will be documented in the [PyDiS documentation](../core_libraries/pydis_documentation/index.rst) section.

Finally, special modules such as `SimulateNetwork` and `VisualizeNetwork` can be defined.

Module `SimulateNetwork` (or `SimulateNetworkPerf`, see [Performance considerations](#performance-considerations)) defines the simulation driver used to run a simulation. It loops through the base modules at each timestep of a simulation to evolve the dislocation network. It must be defined by passing all base modules defined above (`CalForce`, `MobilityLaw`, etc.) plus some additional parameters defining the simulation run, e.g. strain rate, number of steps, etc. For documentation of the `SimulateNetwork` from `pyexadis`, see the [ExaDiS Simulation Driver](../core_libraries/exadis_documentation/user_guide/python_modules.md#simulation-driver) section.

Once a `SimulateNetwork` object is instantiated, the simulation is run with method `run()` passing the initial dislocation network `DisNetManager` object and `state` dictionary as its arguments
```python
sim = SimulateNetwork(...)
sim.run(N, state)
```

Module `VisualizeNetwork` can be used to visualize the evolution of the dislocation network during a simulation, e.g. through a matplotlib window.


#### 6. Running the simulation

After the above python setup is saved in script, e.g. `my_script.py`, it is run by executing the python script
```shell
python my_script.py
```
For debugging or interacting with the simulation, we can also use the interactive mode of python
```shell
python -i my_script.py
```
In this mode, after the script has been executed, one can manipulate the data, e.g. to investigate the dislocation network or perform additional analyses, see [Dislocation network examination](frank_read_src/frank_read_src_by_python.md#dislocation-network-examination) for an example.

````{Hint}
When using `pyexadis`, one must not call the `pyexadis.finalize()` instruction when using the interactive mode, otherwise the ExaDiS memory will be freed and thus inaccessible. A good option can be to call the `pyexadis.finalize()` function only when a non-interactive python mode is detected, e.g.
```python
if not sys.flags.interactive:
    pyexadis.finalize()
```
````

#### 7. Performance considerations

* `pydis` vs. `pyexadis` modules performance:
    * As a pure python code, `pydis` modules are generally slow. They are usually used for learning, prototyping, or running small-scale simulations.
    * ExaDiS is implemented as a HPC code designed to leverage high efficiency of modern computing architectures (including GPUs). Thus `pyexadis` modules are very efficient and intended to be used for production runs and large-scale simulations.

* `pyexadis` provides the two following simulation drivers to run a simulation:
    * `SimulateNetwork`: flexible simulation driver that allows modules from different core libraries to be mixed together, and supports user-defined enhancements. To allow for this, data must flow through the python pipeline, which may result in a small overhead (e.g. 20%), especially if running modules on GPUs.
    * `SimulateNetworkPerf`: high-performance simulation drivers that removes data movement through the python pipeline during a run. This is a direct wrapper to the C++ backend driver implemented in ExaDiS. The downsides is that it cannot be used with external modules (must use `pyexadis` modules only), and it does not allow for interactive data visualization or user-defined enhancements. This driver is recommended for production runs of large systems, e.g. for strain-hardening simulations.
    * See the [ExaDiS Simulation Driver](../core_libraries/exadis_documentation/user_guide/python_modules.md#simulation-driver) section for more information. 
