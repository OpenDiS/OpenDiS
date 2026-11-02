+++
label = "Precipitates"
engine = "exadis"
title = "Precipitates · ExaDiS"
order = 290
script_name = "precipitates_exadis.py"
final_json = "precipitates_exadis_final.json"
tutorial_script = "tutorials/06_simple_precipitates/example_simple_precipitates.py"
tutorial_url = "https://github.com/OpenDiS/OpenDiS/tree/main/tutorials/06_simple_precipitates"
tutorial_title = "06_simple_precipitates"
run_hint = "The 200 steps take about 5 s with ExaDiS. They cover about 2 ns."

[config]
kind = "precipitates"
box_length = 50000.0
theta = 90.0
num_precipitates = 10
radius_min = 500.0
radius_max = 2000.0
seed = 1234

[material]
burgmag = 2.5e-10
mu = 50e9
nu = 0.3
a = 1.0
maxseg_ratio = 0.04
minseg_ratio = 0.004
rann = 2.0
crystal = "fcc"
rtol = 1.0
nextdt = 1e-11
maxdt = 1e-10

[forces]
force_mode = "CUTOFF_MODEL"
cutoff = 5000.0
mobility_law = "FCC_0"
Medge = 60000.0
Mscrew = 60000.0
vmax = 3400.0

[timeint]
integrator = "Trapezoid"

[topology]
topology_mode = "TopologySerial"
collision_mode = "Retroactive"
remesh_rule = "LengthBased"
cross_slip = "ForceBasedSerial"

[loading]
loading_mode = "stress"
applied_stress = [0.0, 0.0, 0.0, 0.0, 1e8, 0.0]

[run]
max_step = 2000
print_freq = 10
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.0001
+++

This example illustrates unshearable precipitates in a DDD simulation. A dislocation line glides through spherical precipitates, which push back its nodes that enter them.
