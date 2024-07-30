### Binary Junction by Python calling ExaDiS

We can run the same Binary Junction example using ExaDiS (after compilation) using the following commands.

```bash
cd ~/Codes/OpenDiS.git
cd examples/03_binary_junction
export OMP_NUM_THREADS=8
python3 -i test_binary_junction_exadis.py
```
This is assuming ExaDis has been compiled using OpenMP, and ```OMP_NUM_THREADS=8``` specifies the number of threads used for the OpenMP run.  The behavior should be similar to the previous test case (using pydis), except that the simulation should run much faster.

#### Explore Dislocation Network

After the simulation has finished, we can examine the data stored in object ```G``` by typing the following command in the Python terminal.
```python
G
```
The output is
```python
<pyexadis_base.ExaDisNet object at 0xffff82a5bc70>
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
This can also be accomplished using the ```get_disnet``` function of DisNetManger ```net```.
```python
from pydis import DisNet
G1=net.get_disnet(DisNet)
```

Now the object ```G1``` is a DisNet object.  We can use the same command as in the previous test case to interact with ```G1```, such as,
```python
G1.all_nodes_tags()
>> dict_keys([(0, 0), (0, 5), (0, 2), (0, 3), (0, 32), (0, 1), (0, 6), (0, 7), (0, 8), (0, 9), (0, 30), (0, 11), (0, 28), (0, 13), (0, 14), (0, 27), (0, 16), (0, 17), (0, 4), (0, 19), (0, 20), (0, 31), (0, 22), (0, 23), (0, 24), (0, 25), (0, 26), (0, 21)])

G1.nodes((0,0)).view()
>> {'R': array([4., 3., 3.]), 'constraint': 7}

list(G1.all_segments_tags())
>> [((0, 0), (0, 17)), ((0, 32), (0, 5)), ((0, 3), (0, 19)), ((0, 27), (0, 7)), ((0, 1), (0, 13)), ((0, 6), (0, 22)), ((0, 7), (0, 23)), ((0, 8), (0, 24)), ((0, 9), (0, 25)), ((0, 20), (0, 26)), ((0, 11), (0, 27)), ((0, 20), (0, 28)), ((0, 13), (0, 21)), ((0, 14), (0, 30)), ((0, 28), (0, 8)), ((0, 16), (0, 32)), ((0, 17), (0, 9)), ((0, 4), (0, 31)), ((0, 19), (0, 11)), ((0, 30), (0, 2)), ((0, 31), (0, 20)), ((0, 22), (0, 14)), ((0, 23), (0, 21)), ((0, 24), (0, 16)), ((0, 25), (0, 1)), ((0, 26), (0, 6)), ((0, 21), (0, 4))]

G1.segments(((0, 3), (0, 19))).view()
>> {'source_tag': (0, 3), 'target_tag': (0, 19), 'burg_vec': array([ 1., -1.,  1.]), 'plane_normal': array([-0.70710678,  0.        ,  0.70710678])}

G1.segments(((0, 3), (0, 19))).burg_vec_from((0, 3))
>> array([ 1., -1.,  1.])

G1.segments(((0, 3), (0, 19))).burg_vec_from((0, 19))
>> array([-1.,  1., -1.])
```
