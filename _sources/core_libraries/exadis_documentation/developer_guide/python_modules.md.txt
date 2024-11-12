
## Calling user-defined python modules

ExaDiS python interface allows the code to interact with user-defined, python-based modules implemented within the specifications of OpenDiS. For instance, a generic module that performs calculations and modifications on the dislocation network can be created as:
```Python
class MyModule:
    """ Prototype of an python-based OpenDiS module """
    def __init__(self, state: dict, **kwargs) -> None:
        # do some initialization
        pass
        
    def some_internal_function(self, G, param):
        # compute some stuff
        # values = ...
        return values
        
    def Compute(self, N: DisNetManager, state: dict) -> dict:
        
        # Get dislocation network in module-specific format
        G = N.get_disnet(MyLibraryDisNet)
        
        # Make some calculations on the network
        values = self.some_internal_function(G, state["param"])
        
        # Make some changes to the network
        G, flag = change_network(G)
        
        # Update state dictionary
        state["nodevalues"] = np.array(values)
        state["param"] = param
        state["flag"] = flag
        # other state values...
        
        return state
```
and the results can then be passed to ExaDiS modules, e.g.:
```Python
# Create initial network
G = ExaDisNet(...)
N = DisNetManager(G)
# Compute forces
forces = CalForce(...)
forces.Compute(N, state)
# Perform some user-defined computations/modifications on the network
mymodule = MyModule(...)
mymodule.Compute(N, state)
# Recompute forces on the modified network
forces.Compute(N, state)
```

It is also possible to create user-defined, python-based `CalForce`, `MobilityLaw`, etc, modules. For instance, we can create a python-based `CalForce` module following the OpenDiS specification template below

```Python
class MyCalForceModule:
    """ Template for a CalForce python-based module: """
    def __init__(self, state: dict, **kwargs) -> None:
        # Do some initialization
        pass
        
    def PreCompute(self, N: DisNetManager, state: dict) -> dict:
        # Do some force pre-computation
        # ...
        return state
    
    def NodeForce(self, N: DisNetManager, state: dict) -> dict:
        # Get dislocation network in module-specific format
        G = N.get_disnet(MyLibraryDisNet)
        # Compute forces on network G
        # f = ...
        # Store forces in state dictionary
        state["nodeforces"] = np.array(f)
        state["nodeforcetags"] = data.get("nodes")["tags"]
        return state
    
    def OneNodeForce(self, N: DisNetManager, state: dict, tag, update_state=True) -> np.array:
        # Get dislocation network in module-specific format
        G = N.get_disnet(MyLibraryDisNet)
        tags = G.get_tags()
        # Find node index for which to compute the force
        ind = np.where((tags[:,0]==tag[0])&(tags[:,1]==tag[1]))[0]
        if ind.size != 1:
            raise ValueError("Cannot find node tag (%d,%d) in OneNodeForce" % tuple(tag))
        # Compute force on node tag and return it
        # fnode = ...
        # Update node force in state dictionary if needed:
        if update_state:
            if "nodeforces" in state and "nodeforcetags" in state:
                nodeforcetags = state["nodeforcetags"]
                ind = np.where((nodeforcetags[:,0]==tag[0])&(nodeforcetags[:,1]==tag[1]))[0]
                if ind.size == 1:
                    state["nodeforces"][ind[0]] = fnode
                else:
                    state["nodeforces"] = np.vstack((state["nodeforces"], fnode))
                    state["nodeforcetags"] = np.vstack((state["nodeforcetags"], tag))
            else:
                state["nodeforces"] = np.array([fnode])
                state["nodeforcetags"] = np.array([tag])
        return fnode
```
and call this module from within ExaDiS, even when compiled/running on GPU, e.g.

```Python
# Declare modules
mycalforce = MyCalForceModule(...) # user-defined class
exadis_mobility = MobilityLaw(...) # imported from pyexadis_base
exadis_timeint = TimeIntegration(..., force=mycalforce, mobility=exadis_mobility) # imported from pyexadis_base
# Integrate the system
mycalforce.NodeForce(N, state)
exadis_mobility.Mobility(N, state)
exadis_timeint.Update(N, state)
```

This usage is allowed by the `CalForcePython`/`ForcePython` wrappers implemented in ExaDiS enabling python objects to be called from the C++ code. The same mechanism allows for mixing ExaDiS and PyDiS modules in an OpenDiS simulation, e.g. as exemplified in the `examples/02_frank_read_src/test_frank_read_src_pydis_exadis.py` example file in the [OpenDiS](https://github.com/OpenDiS/OpenDiS) repository.
