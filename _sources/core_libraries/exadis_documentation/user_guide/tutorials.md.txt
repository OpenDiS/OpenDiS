## Tutorials

### Basics

ExaDiS is a modular DDD code. Contrary to other DDD codes such as the legacy ParaDiS code, standard DDD simulations in ExaDiS are not driven via a set of input files provided as argument to the command line code. Instead, simulations are driven either through the python or C++ interfaces via scripts. ExaDiS modules can also be used as independent building blocks (e.g. within the OpenDiS framework) to create various DDD applications.

```{important}
All paths, e.g. `examples/02_frank_read_src`, used in this section of the documentation refer to paths within the ExaDiS directory. If ExaDiS is installed as a submodule within OpenDiS, the corresponding path is

`${OPENDIS_DIR}/core/exadis/examples/02_frank_read_src`
```


### Running a Frank-Read source simulation

As a starting example, folder `examples/02_frank_read_src` contains an example of a simulation of a simple Frank-Read source, either driven through the python or C++ interfaces. To run the simulation using the python interface, the code needs to be built with flag `-DEXADIS_PYTHON_BINDING`. Then the simulation can be run with:
```shell
cd examples/02_frank_read_src
python test_frank_read_src.py
```
A successful run should display a matplotlib window showing the evolution of the Frank-Read source as the simulations proceeds. 

To run the simulation with the C++ interface, the code must be built with flag `-DEXADIS_BUILD_TESTS`. Then the simulation can be run with:
```shell
cd build/examples/02_frank_read_src
./test_frank_read_src
```
In contrast to the python run, no visualization window will be displayed. Simulation results will be written to the `output` directory, including dislocation configurations in the `.data` format.


### Simulation examples
Folder `examples/` contains a number of ExaDiS simulation examples:
* `02_frank_read_src`: Frank-Read source simulation under constant stress
* `03_collision`: Collision test case between two dislocations
* `04_bcc_junction`: Formation of a binary junction in a BCC crystal
* `05_fcc_junctions`: Formation of several types of binary reactions in a FCC crystal
* `06_nodal_xslip`: Example of nodal x-slip mechanisms (zipping and unzipping) in an elementary BCC network
* `20_compare_forces`: Example of a script to compare force calculation between pydis and ExaDiS modules. This example requires OpenDiS.
* `21_bcc_Ta_100nm_2e8`: Example of a large-scale strength simulation of a BCC tantalum crystal loaded at ultra-high strain rate
* `22_fcc_Cu_15um_1e3`: Example of a large-scale strain-hardening simulation of a FCC copper crystal loaded at high strain rate


### Setting up a python driven simulation
This section shows how to setup and drive an ExaDiS simulation using a python script.

Note: a major advantage of driving a simulation through the python interface is that is does not need compilation of the simulation file, as is required when driving a simulation through the C++ interface.

#### Importing the pyexadis modules
First, we need to make sure that the path to the `pyexadis` library is included in the python path. A successful compilation of ExaDiS produces library file `python/pyexadis.cpython-*.so`. To add the path of this library to python, we can add the following lines at the beginning of our script, for instance:
```python
pyexadis_path = '/path/to/your/exadis/python/'
if not pyexadis_path in sys.path: sys.path.append(pyexadis_path)
```
When the path is properly setup, we should now be able to import `pyexadis` from our script.
```python
import pyexadis
```
Library `pyexadis` contains binding to the backends C++ functions implemented in ExaDiS. To facilitate access to the ExaDiS functions through python, wrapper classes and functions are implemented in the `pyexadis_base.py` and `pyexadis_utils.py` files located in the `python/` folder. ExaDiS data structures (e.g. dislocation network class) and modules (e.g. force, mobility, etc...) can be directly imported from these files:
```python
from pyexadis_base import ExaDisNet, DisNetManager, SimulateNetwork, VisualizeNetwork
from pyexadis_base import CalForce, MobilityLaw, TimeIntegration, Collision, Topology, Remesh
```
Here we have imported the following base classes and modules:
* `ExaDisNet`: ExaDiS dislocation network (graph) object. Handles synchronization between the different execution spaces (e.g. memory movements between CPU and GPU).
* `DisNetManager`: dislocation network manager to convert different instances of a dislocation network between various data structures used in OpenDiS. Any network passed to the `pyexadis_base` wrapper functions must be wrapped into a `DisNetManager` object.
* `SimulateNetwork`: python simulation driver
* `VisualizeNetwork`: python visualization toggle
* `CalForce`: wrapper class to instantiate a model to calculate dislocation forces
* `MobilityLaw`: wrapper class to instantiate a dislocation mobility law
* `TimeIntegration`: wrapper class to instantiate a time-integration scheme
* `Collision`: wrapper class to instantiate a dislocation collision method
* `Topology`: wrapper class to instantiate a dislocation topology method
* `Remesh`: wrapper class to instantiate a dislocation line remeshing method

#### Initializing Kokkos
Now that we have imported all the base classes we need, we can start setting up the simulation. All the simulation script must be inserted in between calls to the `pyexadis.initialize()` and `pyexadis.finalize()` functions:
```python
pyexadis.initialize()

# Code to setup and run the simulation here...

pyexadis.finalize()
```
This is to make sure that Kokkos is initialized and terminated properly on the chosen execution spaces. The `pyexadis.initialize()` can take the following arguments:
* `num_threads`: specifies the number of OpenMP threads to be used. For instance, use `pyexadis.initialize(num_threads=8)` to run a simulation using 8 threads. If not specified (default), ExaDiS will use the maximum number of threads available on the machine or defined from the `OMP_NUM_THREADS` environment variable if set.
* `device_id`: specifies the id of the device (GPU) to be used. If not specified (default), ExaDiS will select the first entry in the list of available devices.

#### Setting up the initial dislocation configuration
After this, the first step is generally to setup the initial dislocation configuration using an `ExaDisNet` object. This can be done by reading a dislocation configuration file, e.g. in ParaDiS format
```python
G = ExaDisNet()
G.read_paradis('my_config.data')
```
or using a user-defined python function that creates an initial configuration, e.g. as in the `examples/02_frank_read_src\test_frank_read_src.py` example
```python
G = init_frank_read_src_loop(pbc=False)
```

Then, we need to wrap the `ExaDisNet` object into a `DisNetManager` object
```python
N = DisNetManager(G)
```
so that it can be used in ExaDiS modules. This step is necessary to make the ExaDiS modules and data structures compliant with the OpenDiS framework.

#### Global state dictionary
The rest of the simulation script consists in defining the global `state` dictionary and the different modules to be used. The global `state` dictionary holds a collection of variables and arrays used during a simulation. Modules operate and communicate with each other by reading and writing entries to the `state` dictionary. The `state` dictionary must be initialized with the global simulation parameters, e.g.
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
Global simulation parameters are used internally in ExaDiS and must be consistent across the different modules used in a simulation. Some parameters are required and other optional, although it is good practice to specify a maximum number of these. Here is a full list of possible parameters:
* `burgmag` (required): Burgers vector magnitude in m. All lengths used in a simulation will be scaled by this value.
* `mu` (required): Shear modulus in Pa.
* `nu` (required): Poisson's ratio.
* `a` (required): Dislocation non-singular core radius in `burgmag` units.
* `maxseg` (required): Maximum discretization segment length in `burgmag` units.
* `minseg` (required): Minimum discretization segment length in `burgmag` units.
* `crystal` (optional): Crystal structure. Optional but most mobility laws will require it to be provided. Supported options are `'fcc'` and `'bcc'` for now. Default: none.
* `Rorient` (optional): Crystal orientation matrix. Default: Identity matrix.
* `rtol` (optional): Integration tolerance in `burgmag` units. Default: 0.25*`a`.
* `rann` (optional): Annihilation radius in `burgmag` units. Default: 2*`rtol`.
* `maxdt` (optional): Maximum timestep size. Default: 1e-7 s.
* `nextdt` (optional): Initial timestep size. Default: 1e-12 s.
* `split3node` (optional): Toggle to enable splitting of 3-nodes (only for BCC crystals). Default: 1.

#### Simulation modules
Next we define the simulation modules to be used in the simulation. A full simulation typically uses the following modules corresponding to the different elementary stages of the DDD simulation cycle: `CalForce`, `MobilityLaw`, `TimeIntegration`, `Collision`, `Topology`, `Remesh`. Each ExaDiS module must be instantiated by passing the global state dictionary `state` and module specific arguments, e.g.
```python
mobility = MobilityLaw(mobility_law='FCC_0', state=state, Medge=64103.0, Mscrew=64103.0, vmax=4000.0)
```
Modules parameters and options are documented in the Modules section of this guide.

Finally, special modules `SimulateNetwork` (or `SimulateNetworkPerf`) and `VisualizeNetwork` can be defined. Module `SimulateNetwork` (or `SimulateNetworkPerf`) is the simulation driver used to run a simulation. It must be defined by passing all base modules defined above (`CalForce`, `MobilityLaw`, etc.) plus some additional parameters defining the simulation run, e.g. strain rate, number of steps, etc. Once the `SimulateNetwork` object is instantiated, the simulation is run with method `run()` passing the initial dislocation network `DisNetManager` object and `state` dictionary as its arguments
```python
sim = SimulateNetwork(...)
sim.run(N, state)
```

Module `VisualizeNetwork` is used to visualize the evolution of the dislocation network during a simulation, e.g. through a matplotlib window.


#### Running the simulation
After the above python driven simulation setup is saved in script, e.g. `my_script.py`, it is run by executing the python script
```shell
python my_script.py
```
For debugging or interacting with the simulation, we can also use the interactive mode of python
```shell
python -i my_script.py
```
In this mode, after the script has been executed, we can manipulate the data, e.g. to investigate the dislocation network or perform additional analyses. Note that we must not call the `pyexadis.finalize()` instruction when using the interactive mode, otherwise the ExaDiS data will be freed and thus inaccessible. A good option is to call the `pyexadis.finalize()` function only when a non-interactive python mode is detected, e.g.
```python
if not sys.flags.interactive:
    pyexadis.finalize()
```

#### Performance considerations
ExaDiS is implemented as a HPC code designed to leverage high efficiency of modern computing architectures. However, when using the python interface, ExaDiS may suffer some performance hit (compared to a native C++ driven simulation) due to data flow through the python pipeline. This overhead is generally reasonable (e.g. 20%) and allows for a great flexibility of the code, e.g. ExaDiS modules can be used within OpenDiS and interfaced with other DDD implementations.

However, when using the python interface to drive a standalone ExaDiS simulation (i.e. outside of OpenDiS or when all modules used are ExaDiS modules) and high performance is desired, then we can simply use the `SimulateNetworkPerf` driver instead of the default `SimulateNetwork` driver. The `SimulateNetworkPerf` performance driver removes data movement through the python interface during the run and should results in simulation performance identical to when using the C++ driving interface. The downsides are that it cannot be used with external modules (must use ExaDiS modules only), and it does not allow for interactive data visualization or manipulation. This driver is recommended for production runs of large systems, e.g. for strain-hardening simulations.


### Using ExaDiS modules to build applications
The python interface and modular design of ExaDiS allows for modules to be used as independent building blocks to build various DDD applications. A simple example is given in `examples/20_compare_forces`, where ExaDiS is only used to evaluate the nodal forces on a given network configuration (as opposed to running a full simulation).


### Using ExaDiS modules in OpenDiS
ExaDiS modules can be imported and used within the OpenDiS python framework. This allows modules to be mixed with external modules coming from other DDD simulations (e.g. pydis).

To use ExaDiS modules within OpenDiS, the only required change is to use the `DisNetManager` class from the OpenDiS framework instead of the one provided in `pyexadis_base.py`, i.e. replace
```python
from pyexadis_base import DisNetManager
```
with
```python
from framework.disnet_manager import DisNetManager
```
where `framework` designates the `python/framework` directory of the OpenDiS code that must be included in the python path. Then ExaDiS modules are essentially used the same way than has been explained in the previous section. Please refer to the OpenDiS documentation for more information and examples on how to use ExaDiS from the OpenDiS framework.


### Setting up a C++ driven simulation
As ExaDiS is written in C++, it also gives the option to setup and drive simulations directly through the C++ interface. This can be useful for certain purposes and allows for a more direct control of the simulations via access to the native implementation. However, for standard applications such as running a typical DDD simulation, it is generally recommended to use the python interface, which is generally easier to use. Beyond syntaxic differences, a major difference is that a simulation driven though C++ must be compiled before it can be executed.

Overall, the way simulations are setup in the C++ interface is similar to the python approach. A good starting point is to look at the Frank-Read simulation example in `examples/02_frank_read_src` and compare the python driven setup in `test_frank_read_src.py` with the C++ driven setup in `test_frank_read_src.cpp`. The typical workflow consists in instantiating a simulation object, declaring the simulation parameters, initializing the dislocation configuration, setting up the various modules, and running the simulation.

#### Includes and namespace
First, we need to include the required header files for ExaDiS and the simulation driver functionalities. We can also use the ExaDiS namespace to simplify code syntax in the file.
```cpp
#include "exadis.h"
#include "driver.h"

using namespace ExaDiS;
```

#### Initializing the application
Then we need to initialize the application, by first initializing Kokkos (here using the `ScopeGuard` approach). An instance of a DDD simulation application is then created using the `ExaDiSApp` driver class from `driver.h`.
```cpp
Kokkos::ScopeGuard guard(argc, argv);
ExaDiS::ExaDiSApp exadis(argc, argv);
```

#### Setting up the simulation parameters
Global simulation parameters are defined using the `Params` struct. The global parameters are essentially the same as those detailed in the python section, expect for the crystal information which are defined separately using the `Crystal` class. The constructor of the `Params` object takes in the six required simulation parameters (`burgmag`, `MU`, `NU`, `a`, `maxseg`, `minseg`), and other parameters can be specified manually afterwards, e.g.
```cpp
Params params(burgmag, MU, NU, a, maxseg, minseg);
params.nextdt = dt;
params.rann = rann;
params.rtol = rtol;
```
Another set of parameters are the simulation control parameters that specify the control settings such as the number of simulation steps, loading type, applied stress, output frequencies, etc. Each loading instance is defined using the `ExaDiSApp::Control` struct, for instance
```cpp
ExaDiSApp::Control ctrl;
ctrl.nsteps = 200;
ctrl.loading = ExaDiSApp::STRESS_CONTROL;
ctrl.appstress = Mat33().symmetric(0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0);
ctrl.erate = 0.0;
ctrl.edir = Vec3(0.0, 0.0, 1.0);
ctrl.rotation = 0;
ctrl.printfreq = 1;
ctrl.propfreq = 10;
ctrl.outfreq = 100;
using Prop = ExaDiSApp::Prop;
ctrl.props = {Prop::STEP, Prop::STRAIN, Prop::STRESS, Prop::DENSITY};
```
where the control properties are:
* `nsteps` (required): number of steps for the loading instance. Valid options are:
    - `N`: the number of steps N to run from this instant
    - `ExaDiSApp::NUM_STEPS(N)`: same as above
    - `ExaDiSApp::MAX_STEPS(N)`: the maximum number of steps N to reach
    - `ExaDiSApp::MAX_STRAIN(strain)`: the maximum strain to reach
    - `ExaDiSApp::MAX_TIME(time)`: the maximum physical simulation time in s to reach
    - `ExaDiSApp::MAX_WALLTIME(time)`: the maximum wall-clock time in s to reach
* `loading` (required): loading type. Available options are:
    - `ExaDiSApp::STRESS_CONTROL`: constant stress loading. The stress value is defined using the `appstress` control property.
    - `ExaDiSApp::STRAIN_RATE_CONTROL`: constant, uniaxial strain rate loading. The strain rate is specified using property `erate`. The loading direction is specfified using property `edir`. The simulation can be started at a specific initial value of stress defined using the `appstress` control property.
* `erate` (optional): loading strain rate in 1/s for the `STRAIN_RATE_CONTROL` loading mode. Unused in `STRESS_CONTROL` mode. Default: 1e3/s.
* `edir` (optional): uniaxial loading direction for the `STRAIN_RATE_CONTROL` mode. Unused in `STRESS_CONTROL` mode. Default: (0,0,1)
* `appstress` (optional): initial value of the applied stress tensor. Default: zero matrix.
* `rotation` (optional): flag to enable crystal rotation. Default: 0
* `printfreq` (optional): step frequency to print console information
* `propfreq` (optional): step frequency to print simulation properties
* `outfreq` (optional): step frequency to output dislocation configurations and restart files
* `props` (optional): list of simulation properties to be printed in the output property file at frequency `propfreq`. Default: {STEP, STRAIN, STRESS, DENSITY}. Available items are:
    - `ExaDiSApp::Prop::STEP`: simulation step number
    - `ExaDiSApp::Prop::STRAIN`: uniaxial strain
    - `ExaDiSApp::Prop::STRESS`: uniaxial stress
    - `ExaDiSApp::Prop::DENSITY`: total dislocation density
    - `ExaDiSApp::Prop::NNODES`: number of dislocation nodes
    - `ExaDiSApp::Prop::NSEGS`: number of dislocation segments
    - `ExaDiSApp::Prop::DT`: timestep size
    - `ExaDiSApp::Prop::TIME`: simulation time
    - `ExaDiSApp::Prop::WALLTIME`: wall-clock time
    - `ExaDiSApp::Prop::EDIR`: current loading direction
    - `ExaDiSApp::Prop::ALLSTRESS`: full stress components

#### Setting up the initial dislocation configuration
The initial dislocation network is set up using an instance of a `SerialDisNet` object. The `SerialDisNet` class is a simple STL-based implementation of an ExaDiS dislocation network that is used to create, hold and manipulate dislocation networks on the host execution space (CPU). User can create their dislocation network easily by implementing user-defined functions (e.g. see function `init_frank_read_src_loop()` in `test_frank_read_src.cpp`), or using network generation functions implemented in ExaDiS. For instance, we can create an initial configuration made of prismatic loops with
```cpp
Crystal crystal(BCC_CRYSTAL);
SerialDisNet* config = generate_prismatic_config(crystal, Lbox, numloops, radius, maxseg, seed);
```
or load a configuration form a ParaDiS `.data` file
```cpp
SerialDisNet* config = read_paradis(datafile);
```
Once we have an initial configuration, we need to initialize the ExaDiS `System` object of the driver
```cpp
System* system = exadis->system;
system->initialize(params, crystal, config);
```
`System` is the base class of ExaDiS that holds all the information about the dislocation system being simulated (global parameters, crystal, dislocation network, etc.). All ExaDiS modules and functions usually require the system's pointer as one of their argument. Once the dislocation configuration is passed to the system `initialize()` function, it becomes available to use on kernels of all execution spaces (e.g. GPU).

#### Modules initialization
Next we need to define and initialize the various modules needed in the simulation such as the force model, mobility, integrator, collision, topology, and remeshing components. Each of these six modules corresponds to a stage of the traditional DDD cycle, and is required to be defined, e.g.
```cpp
exadis->force = new ForceType::LINE_TENSION_MODEL(system);
exadis->mobility = new MobilityType::GLIDE(system, MobilityType::GLIDE::Params(Mob));
exadis->integrator = new IntegratorEuler(system);
exadis->collision = new CollisionRetroactive(system);
exadis->topology = new TopologySerial(system, exadis->force, exadis->mobility);
exadis->remesh = new RemeshSerial(system);
```

#### Simulation setup and run
Finally, we need to call the `set_simulation()` function of the driver to make sure everything is set up properly and perform some sanity checks
```cpp
exadis->outputdir = outputdir;
exadis->set_simulation();
```
We need to make sure that the `outputdir` parameter of the driver is properly specified before calling `set_simulation()`, because this function also handles the creation of directories necessary for the simulation.

Once everything is setup, we can call the `run()` function to perform the simulation
```cpp
exadis->run(ctrl);
```
The `run()` function takes as an argument an instance of a simulation control `ExaDiSApp::Control` object, or a vector of such objects (allowing to perform a sequence of various loadings).

#### Compilation
After the above C++ driven simulation setup is saved in a cpp file, e.g. `my_simulation.cpp`, it must be compiled to produce a simulation executable. This can be done by adding a new executable to the CMake build tree of ExaDiS and linking against the ExaDiS library, e.g.
```cmake
add_executable(my_simulation my_simulation.cpp)
target_include_directories(my_simulation PRIVATE ${EXADIS_INCLUDE_DIRS})
target_link_libraries(my_simulation PRIVATE exadis ${EXADIS_EXTERN_LIBS})
```

For more details about the backend C++ classes and modules options, please see the [Developer guide](../developer_guide/index) section of the documentation.
