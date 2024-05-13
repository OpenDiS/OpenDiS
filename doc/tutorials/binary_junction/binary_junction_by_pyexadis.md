### Binary Junction by Python calling ExaDiS

**Run OMP version**

By calling exadis, you can run OMP version where you are able to control the number of OMP threads.

```bash
cd ~/Codes/OpenDiS.git/examples/03_binary_junction
export OMP_NUM_THREADS=8
python3 -i test_binary_junction_exadis.py (to be created)
```

The initial conditions, boundary conditions, and expected dislocation behaviors are identical to those written in [Sec. 3.2.1](https://caiwei-stanford.github.io/opendis-doc/tutorials/binary_junction/binary_junction_by_python.html). In fact, you will be able to observe significant speed increase since nunber of computing threads has increased to 8 via ```export OMP_NUM_THREADS=8```.
