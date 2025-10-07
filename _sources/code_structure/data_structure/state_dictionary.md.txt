## State Dictionary

The global `state` dictionary is used to hold a collection of variables and arrays defining the state of a simulation. It is propagated through the modules during a simulation: modules operate and communicate with each other by reading and writing entries to the `state` dictionary.

The `state` dictionary thus serves two purposes:
* Define the global simulation parameters (e.g. materials parameters) that are common to several simulation modules
* Serve as a memory buffer for modules to store / retrieve information (e.g. force calculation modules write nodal forces into the `state` dictionary so that other modules can fetch them)

### Global simulation parameters

Before initializing modules and running a simulation, the `state` dictionary must be initialized with the global simulation parameters, e.g.

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

The global simulation parameters are generally those parameters that are shared by several modules and defining them in the `state` dictionary ensures that their value will be consistent across the different modules used in a simulation. Some parameters are required and other are optional, although it is good practice to specify a maximum number of these. A description of the typical parameters is given at [Setting up a simulation section](../../tutorials/setting_simulation.md#defining-the-global-state-dictionary).

By specification, the `state` dictionary needs to be passed as an argument to modules constructors.


### Module state variables

Modules can write and read entries in the `state` dictionary. By specification, some modules are required to populate or look up for certain items in the `state` dictionary. For instance, the `NodeForce()` method of a `CalForce` module is required to write node forces information into the `state` dictionary while the `Mobility()` method of a `MobilityLaw` module is required to look up for the nodal forces values stored in the `state` dictionary, e.g.

```python
calforce = CalForce(state=state, ...)
mobility = MobilityLaw(state=state, ...)

# Compute nodal forces and store them into the state dictionary
state = calforce.NodeForce(N, state)
print('nodal forces', state["nodeforces"])

# Compute nodal velocities from the nodal forces stored into the state dictionary
state = mobility.Mobility(N, state)
print('nodal velocities', state["nodevels"])
```

In general, each module is free to dump variables into the `state` dictionary, e.g. to save information or communicate with other modules. A good rule for determining variables that need to be saved by a module is to consider the `state` dictionary as information required to restart a simulation from the current state.
