+++
label = "Uniaxial tension"
engine = "exadis"
title = "Uniaxial tension · ExaDiS"
order = 250
script_name = "uniaxial_tension_exadis.py"
final_json = "uniaxial_tension_exadis_final.json"
tutorial_script = "examples/05_uniaxial_tension/test_uniaxial_tension_exadis.py"
run_hint = "The 500 steps take some seconds with ExaDiS."

[config]
kind = "loop_sources"
box_length = 10000.0
num_loops = 4
radius_ratio = 0.21
seed = 123
uniform = true

[material]
burgmag = 2.85e-10
mu = 55.0e9
nu = 0.339
a = 6.0
maxseg_ratio = 0.1
minseg_ratio = 0.02
rann = 1.0
crystal = "bcc"
rtol = 1.0
nextdt = 1e-12

[forces]
force_mode = "LineTension"
mobility_law = "BCC_0B"
Medge = 2600.0
Mscrew = 2600.0
Mclimb = 1e-4
vmax = 3400.0

[timeint]
integrator = "Trapezoid"
multi = 2

[topology]
topology_mode = "TopologyParallel"
collision_mode = "Retroactive"
remesh_rule = "LengthBased"


[loading]
loading_mode = "strain_rate"
erate = 3e3
edir = [0.0, 0.0, 1.0]

[run]
max_step = 500
print_freq = 1
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.01
+++

Prismatic loops in bcc tantalum expand as the box is pulled along z at a constant strain rate. The loops are generated as in the tutorial's [ensemble of prismatic loops](https://opendis.github.io/OpenDiS/tutorials/initial_configuration.html#ensemble-of-prismatic-loops). The history view plots the stress against the strain.
