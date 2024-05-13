### Binary Junction by Pure Python

**Introduction**

This example shows a binary junction formation simulation. Binary junctions are formed when two different dislocations overlap each other. It is a result of an energy relaxation process where the two dislocations with disticnt burgers vector are merged. More information about binary junctions can be found in ([Bulatov & Cai, 2006](https://core.ac.uk/reader/44178170)). 

it produces additional dislocation line with distinct burgers vector

The code can be simply run by PyDis in pure Python mode where no compilation is needed.

Two dislocation lines becomes into three dislocation lines.
This is part of the energy minimization process.
Two dislocation is merged when they are in contact --> junction.

```bash
cd ~/Codes/OpenDiS.git/examples/03_binary_junction
python3 -i test_binary_junction.py
```



**Initial condition**

<img src=./figures/binary_junction_python_init.png alt="" width="500" />

Initial node positions are assigned to variable ```rn``` as
```
rn    = np.array([[0.0, -z0, -z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  z0,  z0,  DisNode.Constraints.PINNED_NODE],
                      [-z0,  0.0,-z0,  DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0, 0.0, DisNode.Constraints.UNCONSTRAINED],
                      [ z0,  0.0, z0,  DisNode.Constraints.PINNED_NODE]])
```


**Boundary condition**

By default, the boundary condition is set to ```pbc=False``` which disables the periodic boundary condition. You can easily set ```pbc=True``` by

```
net = init_two_disl_lines(z0=4000, box_length=3.5e4, pbc=False)
```

while initialization process. 




<mark>To do: describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.</mark>
<mark>Done: describe this simulation.  Initial condition.  Boundary condition.  Behavior of the dislocation.</mark>



**Expected result**

<img src=./figures/binary_junction_python.png alt="" width="500" />

Two different dislocations forms a binary junction where they intersect with each other, which results in three distinct burgers vecctors. Simulation result after ```max_step=200``` with time step size ```dt0 = 1.0e-8``` is shown in the figure above.

it has an independent burgers vector. 

The code will produce animation of binary junction formation process. The following two figures shows the result of the simulation.



The first figure shows the initial configuration of the simulation. 


```{attention}
Climate change is real.
```
