### Frank-Read Source by Python calling ExaDiS

We can run the same Frank-Read Source example using ExaDiS (after compilation) using the following commands.

```bash
cd ~/Codes/OpenDiS.git
cd examples/02_frank_read_src
export OMP_NUM_THREADS=8
python3 -i test_frank_read_src_exadis.py
```
This is assuming ExaDis has been compiled using OpenMP, and ```OMP_NUM_THREADS=8``` specifies the number of threads used for the OpenMP run.  The behavior should be similar to the previous test case (using pydis), except that the simulation should run much faster.

#### Explore Dislocation Network

After the simulation has finished, we can examine the data stored in object ```G```.
```python
G
```
The output is
```python
<pyexadis_base.ExaDisNet object at 0x10342fd30>
```
We can see that ```G``` is an ```ExaDiSNet``` object (so that we cannot interact with it the same way as a ```DisNet``` object).

We can use the following command to see the content of ```G```.
```python
G.export_data()
```

We can use the following command to convert an ```ExaDiSNet``` object to a ```DisNet``` object, so that we can interact with it using the same approach as in the previous test case.
```python
from pydis import DisNet
G1=DisNet()
G1.import_data(G.export_data())
```

```{hint}
This can also be accomplished using the function of DisNetManger ```net```.
```python
from pydis import DisNet
G1=net.get_disnet(DisNet)
```

Now the object ```G1``` is a DisNet object.  We can use the same command as in the previous test case to interact with ```G1```, such as,
```python
G1.all_nodes_tags()
G1.nodes((0,0)).view()
list(G1.all_segments_tags())
G1.segments(((0, 0), (0, 4))).view()
G1.segments(((0, 0), (0, 4))).burg_vec_from((0,4))
G1.segments(((0, 0), (0, 4))).burg_vec_from((0,0))
```
