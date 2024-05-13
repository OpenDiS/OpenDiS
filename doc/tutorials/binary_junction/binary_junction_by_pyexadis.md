### Binary Junction by Python calling ExaDiS

**Run OMP version**

By calling exadis, you can run OMP version where you are able to control the number of OMP threads.

```bash
cd ~/Codes/OpenDiS.git/examples/03_binary_junction
export OMP_NUM_THREADS=8
python3 -i test_binary_junction_exadis.py (to be created)
```

The initial conditions, boundary conditions, and expected dislocation behaviors are identical to those written in section 
<mark>To do: create this test case – copy from the similar folder in core/exadis/examples</mark>
