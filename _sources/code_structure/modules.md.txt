## Modules Structure

### Specifications

OpenDiS relies on a modular architecture: each functionality (e.g. force calculation, mobility law, post-processing analysis, etc.) is implemented as a module, exposed as a python class to the framework.

Below is a typical prototype for an OpenDiS module:

```python
class MyModule:
    """ Prototype of an OpenDiS module """
    def __init__(self, state: dict, **kwargs) -> None:
        # do some initialization
        pass
        
    def Compute(self, N: DisNetManager, state: dict) -> dict:
        
        # Get dislocation network in module-specific format
        G = N.get_disnet(MyLibraryDisNet)
        
        # Perform some computations on the network
        nodevalues = compute_values(G, state)
        
        # Make some modifications to the network
        G, mymoduleflags = modify_network(G, state)
        
        # Update state dictionary
        state["nodevalues"] = nodevalues
        state["mymoduleflags"] = mymoduleflags
        
        return state
```

- Each method implemented in a module (e.g. the `Compute()` method) must adhere to the following specifications:
    * Receives a DisNetManager and the state dictionary as input arguments
    * Can perform computation(s) on the network
    * Can modify the network (e.g. topological operations)
    * Returns an updated state dictionary

- Special core modules such as Force calculation and Mobility law modules have additional specifications to guarantee consistency of the implementations:
    * Force calculations modules: see the [CalForce_Base](https://github.com/OpenDiS/OpenDiS/blob/main/python/framework/calforce_base.py) specification class.
    * Mobility law modules: see the [MobilityLaw_Base](https://github.com/OpenDiS/OpenDiS/blob/main/python/framework/mobility_base.py) specification class.

- Compatibility and inter-operability between modules is maintained via the two following objects:
    * `DisNetManager`: container class that allows for co-existence and conversion between various dislocation network data structures. This allows for different modules to use their own data structures to manipulate the networks. See [DisNetManager class](data_structure/disnetmanager_class.md).
    * `state` dictionary: dictionary containing all relevant parameters (e.g. materials parameters, simulation settings, nodal forces, etc.) pertaining to the state of a simulation / analysis. Modules can communicate by writing/reading items into the state dictionary.


### Building a simulation

In the general sense, a simulation (or analysis pipeline) in OpenDiS can be constructed as a sequence of modules, e.g. as illustrated in the figure below:

```{figure} modules_workflow.png
:alt: OpenDiS workflow
:width: 400px
```

where individual modules may come from different core libraries or implementations.

For running traditional DDD simulations, simulation drivers are provided in core libraries (e.g. see [Simulation driver section](../core_libraries/exadis_documentation/user_guide/python_modules.md#simulation-driver)). They loop over a sequence of modules that constitute the traditional DDD cycle, e.g.
1. Force module
2. Mobility module
3. Time-integration module
4. Collision module
5. Topology module
6. Cross-slip module
7. Remesh module
8. Update stress module
9. Output module
