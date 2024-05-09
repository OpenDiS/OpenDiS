# Frank-Read Docs

## Frank-Read Source by Pure Python
Run OpenDiS in pure python mode (PyDiS) -- no compilation needed

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src
python3 -i test_frank_read_src.py
```

```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.
2. how to examine the data structure of the dislocation network (DisNet) in python interactive
```

## Frank-Read Source by Python calling ExaDiS

```bash
cd ~/Codes/OpenDiS.git/examples/02_frank_read_src/
export OMP_NUM_THREADS=8
python3 -i test_frank_read_src_exadis.py
```


Markdown (MyST)
```{figure} frank_read_ex.png
frank-read dislocation (fix caption)
```

This is from Markdown.


```{important}
1. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.

2. can we make the remesh rules consistent between pydis and exadis, no mesh_refine on segments with both ends fixed? Or at least have an option to get to this behavior?

3. describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.  How is this simulation different (e.g. time) from 3.1.1?

4. how to examine the data structure of the dislocation network (DisNet) in python interactive


```
