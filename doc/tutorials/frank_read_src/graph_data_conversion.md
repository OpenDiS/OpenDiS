### Graph Data Conversion

To run the test case, simply execute:

```bash
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
python3 -i test_data_convert.py
```

#### Expected Behavior
This test case first run the same simulation as in ```test_frank_read_src_exadis.py```, and then converts the dislocation network (graph) to various formats and check for self-consistency.  If all tests are successful, you should see


```bash
1. test G1 sanity check PASSED
2. test G1_eq_G2 PASSED
3. test G1_eq_G4 PASSED
```

#### Conversion between ExaDiS and PyDiS using DisNetManager

At the end of the simulation, the final dislocation network is managed by the ```DisNetManager``` called ```net```.  The following line in the python code
```python
G0 = net.get_disnet(ExaDisNet)
```
extracts the dislocation graph ```G0``` in the ```ExaDisNet``` format (which is the native format for the ```ExaDiS``` library).

The following line in the python code
```python
G1 = net.get_disnet(DisNet)
```
converts the dislocation graph to ```G1``` in the ```DisNet``` format, which is more convenient to interact with by a human.  For example, the ```DisNet``` has a function for sanity checks (e.g. the Burgers vectors for all dislocation lines going out of any node sum up to zero), which can be called as follows.
```python
G1_sanity_check = G1.is_sane()
```
If ```G1_sanity_check``` is ```True```, then the first test passes.

#### Conversion using export_data and import_data Functions

Both ```DisNet``` and ```ExaDisNet``` implement the ```export_data()``` and ```import_data()``` functions that can be used to convert data to different formats.  As a self-consistency check, we first create an empty ```DisNet``` object ```G2```, and then import the data exported by ```G1``` as follows.
```python
G2 = DisNet()
G2.import_data(G1.export_data())
```
The result ```G2``` should be a copy that is identical to the original graph ```G1```.  We can verify this by calling the ```is_equivalent()``` function implemented in ```DisNet```.
```python
G1_eq_G2 = G1.is_equivalent(G2) and G2.is_equivalent(G1)
```
If ```G1_eq_G2``` is ```True```, then the second test passes.

#### Conversion to and from Networkx
[Networkx](https://networkx.org/) is a widely used network library.  ```DisNet``` has utility functions to export its graph to ```networkx``` format (```DiGraph```) and to import from ```networkx```. As a self-consistency check, we first export the dislocation graph to ```networkx``` object ```G3```, and then import data from ```G3``` into a new ```DisNet``` object ```G4``` as follows.
```python
G3 = G1.to_networkx()
G4 = DisNet()
G4.from_networkx(G3)
```
The result ```G4``` should again be a copy that is identical to the original graph ```G1```.  We can verify this by calling the ```is_equivalent()``` function implemented in ```DisNet```.
```python
G1_eq_G4 = G1.is_equivalent(G4) and G4.is_equivalent(G1)
```
If ```G1_eq_G4``` is ```True```, then the third test passes.
