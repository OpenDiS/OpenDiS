### Frank-Read Source by Python calling ExaDiS

We can run the same Frank-Read Source example using ExaDiS (after compilation) using the following commands.

```bash
cd ${OPENDIS_DIR}
cd examples/02_frank_read_src
export OMP_NUM_THREADS=8
python3 -i test_frank_read_src_exadis.py
```
This is assuming ExaDis has been compiled using OpenMP, and ```OMP_NUM_THREADS=8``` specifies the number of threads used for the OpenMP run.  The behavior should be similar to the previous test case (using pydis), except that the simulation should run much faster.

One difference between this test case and the previous one is that here we have set
```python
topology  = None
```
The purpose is just to demonstrate that, without a topology module (which handles split_multi_node), the Frank-Read source test case can still behave properly.

#### Explore Dislocation Network

After the simulation has finished, we can examine the data stored in object ```G``` by typing the following command in the Python terminal.
```python
G
```
The output is
```python
<pyexadis_base.ExaDisNet object at 0x10342fd30>
```
We can see that ```G``` is an ```ExaDiSNet``` object (so that we cannot interact with it the same way as a ```DisNet``` object).


We can use the following command to see the content of G.

```python
G.export_data()
```

<details>
  <summary>
   Here is the output for this command.
  </summary>

 ```python
{'cell': {'h': array([[1000.,    0.,    0.],
       [   0., 1000.,    0.],
       [   0.,    0., 1000.]]), 'origin': array([0., 0., 0.]), 'is_periodic': [1, 1, 1]}, 'nodes': {'tags': array([[  0,   0],
       [  0, 107],
       [  0,   2],
       [  0,   3],
       [  0,   4],
       ...,
       [  0,  63],
       [  0,  88],
       [  0,   8],
       [  0,  82],
       [  0, 114]]), 'positions': array([[500.        , 437.5       , 500.        ],
       [504.8985802 , 405.30790288, 500.        ],
       [500.        , 562.5       , 500.        ],
       [500.        , 562.5       , 375.        ],
       [500.        , 437.5       , 375.        ],
       ...,
       [ 98.42129319, 129.36609803, 500.        ],
       [239.65231892,  21.11140067, 500.        ],
       [571.97968178, 379.25113058, 500.        ],
       [529.45283533, 613.42150792, 500.        ],
       [ 83.49577868, 877.26034318, 500.        ]]), 'constraints': array([[7],
       [0],
       [7],
       [7],
       [7],
       ...,
       [0],
       [0],
       [0],
       [0],
       [0]])}, 'segs': {'nodeids': array([[ 0,  1],
       [16, 52],
       [ 2,  3],
       [ 3,  4],
       [ 4,  0],
       ...,
       [29, 49],
       [25, 17],
       [32,  8],
       [ 8, 13],
       [19, 32]]), 'burgers': array([[1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.],
       ...,
       [1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.],
       [1., 0., 0.]]), 'planes': array([[ 0.,  0.,  1.],
       [ 0.,  0.,  1.],
       [-0.,  1.,  0.],
       [ 0.,  0., -1.],
       [ 0., -1.,  0.],
       ...,
       [ 0.,  0.,  1.],
       [ 0.,  0.,  1.],
       [ 0.,  0.,  1.],
       [ 0.,  0.,  1.],
       [ 0.,  0.,  1.]])}}
```
</details>


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
G1.nodes((0,0)).view()
list(G1.all_segments_tags())
G1.segments(((0, 0), (0, 4))).view()
G1.segments(((0, 0), (0, 4))).burg_vec_from((0,4))
G1.segments(((0, 0), (0, 4))).burg_vec_from((0,0))
```
