### Inter-operability between PyDiS and ExaDiS Modules

To run the test case, simply execute:

```bash
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
python3 -i test_frank_read_src_pydis_exadis.py
```

#### Simulation Setup (default)
In this example, two sets of classes (such as ```CalForce```, ```MobilityLaw```, ```Topology```) are imported from both ```PyDiS``` and ```ExaDiS``` modules.
They can be combined in various ways to construct a dislocation dynamics simulation.
This Python script demonstrates two ways to do so, as specified by the optional command line argument, which can be either ```1``` (default) or ```2```.

```{hint}
The last line in the block above has the same effect as
```bash
python3 -i test_frank_read_src_pydis_exadis.py 1
```
In this Python script, the ```sim``` object is still constructed from the ```SimulateNetwork``` class of ```PyDiS```.
Here with ```option``` equal to ```1```, ```sim``` uses ```MobilityLaw```, ```TimeIntegration```, ```Topology```, ```VisualizeNetwork``` classes imported from ```PyDiS```, and uses ```CalForce```, ```Collision```, ```Remesh``` classes imported from ```ExaDiS```.

#### Simulation Setup (option 2)

Alternatively, we can run the python script using option 2 as follows.
```bash
python3 -i test_frank_read_src_pydis_exadis.py 2
```

Now ```sim``` uses ```MobilityLaw```,  ```VisualizeNetwork``` classes imported from ```PyDiS```, and uses ```CalForce```, ```TimeIntegration```, ```Topology```, ```Collision```, ```Remesh``` classes imported from ```ExaDiS```.

It is interesting to note that some of the objects make use of other objects that may be from either ```PyDiS``` or ```ExaDiS``` modules.  For example, ```exadis_timeint``` and ```exadis_topology``` both use ```CalForce``` from ```PyDiS``` and ```Mobility``` from ```ExaDiS```. 

