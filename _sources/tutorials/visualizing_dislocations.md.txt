### Visualizing dislocation evolution

Dislocation evolution can be visualized in several ways:

#### Interactive matplotlib visualization module

Module `VisualizeNetwork` can be used to visualize the evolution of the dislocation network during a simulation through a matplotlib window, e.g.

```python
from pyexadis_base import VisualizeNetwork

vis = VisualizeNetwork()

sim = SimulateNetwork(
    ...,
    vis=vis, plot_freq=10
)
```

```{Important}
Interactive plotting of the dislocation network through matplotlib module can induce a significant overhead to the simulation run time.
```

```{Note}
The interactive visualization window will not display when using the performance driver [`SimulateNetworkPerf`](../core_libraries/exadis_documentation/user_guide/python_interface/python_modules.md#simulatenetworkperf) from `pyexadis`, as this driver is intended for production runs in HPC environments.
```

#### Visualizing .data output files

Simulation drivers [`SimulateNetwork`](../core_libraries/exadis_documentation/user_guide/python_interface/python_modules.md#simulation-driver) and [`SimulateNetworkPerf`](../core_libraries/exadis_documentation/user_guide/python_interface/python_modules.md#simulation-driver) from `pyexadis` regularly output the dislocation network in `.data` files (with frequency defined by `write_freq`) during a simulation run. There are several ways the `.data` dislocation configurations can be visualized:

* [VisIt](https://visit-dav.github.io/visit-website/index.html): the VisIt visualization software supports `.data` files.
```{Note}
Some users have reported issues with animations in the latest versions of VisIt.
```

* [Ovito Pro](https://www.ovito.org/): you will need to use a custom file reader, following [this thread](https://github.com/OpenDiS/OpenDiS/issues/3).
```{Note}
This requires Ovito Pro, and will not work with Ovito Basic.
```

* [ParaView](https://www.paraview.org/): `.data` files can be converted to `.vtk` files for visualization in ParaView using functions from `pyexadis_utils.py`. Here is an example of an OpenDiS script to convert `.data` files in batch:

```python
import os, glob
import pyexadis
from pyexadis_utils import read_paradis, write_vtk
pyexadis.initialize()

output_path = '/path/to/opendis/simulation/output'
for f in glob.glob(output_path+'/*.data'):
    N = read_paradis(f)
    write_vtk(N, os.path.splitext(f)[0]+'.vtk')
```
