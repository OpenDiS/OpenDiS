### Inter-operability between PyDiS and ExaDiS Modules

To run the test case, simply execute:

```bash
cd ~/Codes/OpenDiS.git
cd examples/02_frank_read_src
python3 -i test_frank_read_src_pydis_exadis.py
```

#### Simulation Setup
In this test case, two sets of classes (such as ```CalForce```, ```MobilityLaw```, ```Topology```) are imported from both ```PyDiS``` and ```ExaDiS``` modules.
They can be combined in various ways to construct a dislocation dynamics simulation.
This Python script demonstrates two ways to do so, as specified by the optional command line argument, which can be either ```1``` (default) or ```2```.

```{hint}
The last line in the block above has the same effect as
```bash
python3 -i test_frank_read_src_pydis_exadis.py 1
```
