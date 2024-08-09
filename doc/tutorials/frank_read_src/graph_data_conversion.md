### Graph Data Conversion

To run the test case, simply execute:

```bash
cd ~/Codes/OpenDiS.git
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
converts the dislocation graph to ```G1``` in the ```DisNet``` format, which is more convenient to interact with by a human.
