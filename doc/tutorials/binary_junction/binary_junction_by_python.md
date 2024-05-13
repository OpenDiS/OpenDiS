### Binary Junction by Pure Python

**Introduction**

This example shows how to simulate a binary junction formation process. Binary junctions are formed when two different dislocations overlap each other. It results from an energy relaxation process where the two dislocations with distinct burgers vectors are merged. More information about binary junctions can be found in [Bulatov & Cai, 2006](https://core.ac.uk/reader/44178170). The code can be run by PyDis in pure Python mode, where no compilation is needed.

```bash
cd ~/Codes/OpenDiS.git/examples/03_binary_junction
python3 -i test_binary_junction.py
```

</br>

**Initial conditions**

<img src=./figures/binary_junction_python_init.png alt="" width="500" />

In the beginning, two dislocations with distinct burgers vectors are created. Both of those dislocations are initially represented by three nodes. The initial positions of the nodes of the dislocations are assigned inside of ```init_two_disl_lines()``` function to variable ```rn``` as

```python
rn    = np.array([[0.0, -z0, -z0,  DisNode.Constraints.PINNED_NODE],
                  [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                  [0.0,  z0,  z0,  DisNode.Constraints.PINNED_NODE],
                  [-z0,  0.0,-z0,  DisNode.Constraints.PINNED_NODE],
                  [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                  [ z0,  0.0, z0,  DisNode.Constraints.PINNED_NODE]])
```

Specific values related to the positions of the nodes, simulation box size, and boundary condition are set inside the ```main()``` function by

```python
net = init_two_disl_lines(z0=4000, box_length=3.5e4, pbc=False)
```

Also, other initial conditions such as Young's modulus, Poisson ratio, or external stress conditions can be set within the ```main()``` function as

```python
calforce = CalForce(mu=160e9, nu=0.31, a=0.1, Ec=1.0e6,
                    applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                    force_mode='LineTension')
```

</br>

**Boundary conditions**

By default, the simulation is set to ```pbc=False```, which does not use the periodic boundary condition. You can easily set ```pbc=True``` by

```python
net = init_two_disl_lines(z0=4000, box_length=3.5e4, pbc=False)
```

inside the ```main()``` function. Moreover, it is possible to choose boundary conditions for each direction (x, y, z) by modifying 

```python
cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])
```

inside the ```init_two_disl_lines()``` function.

</br>

**Expected behavior of the dislocation**

<img src=./figures/binary_junction_python.png alt="" width="500" />

Two different dislocations form a binary junction near the region where they intersect each other. As a result, three different dislocation lines are formed, which have distinct burgers vectors. The figure above shows the simulation result after ```max_step=200``` with time step size ```dt0 = 1.0e-8```.

 

The code will produce an animation of the binary junction formation process. The following two figures show the simulation's result.



The first figure shows the initial configuration of the simulation. 


```{attention}
Climate change is real.
```
