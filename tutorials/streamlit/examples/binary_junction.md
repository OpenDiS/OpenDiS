+++
label = "Binary junction"
title = "Binary junction · PyDiS"
order = 200
script_name = "binary_junction_pydis.py"
final_json = "two_disl_lines_final.json"
tutorial_script = "examples/03_binary_junction/test_binary_junction_pydis.py"
tutorial_url = "https://opendis.github.io/OpenDiS/tutorials/binary_junction/binary_junction_by_python.html"
run_hint = "The tutorial's 200 steps take about 8 s."

[config]
kind = "two_disl_lines"
box_length = 8.0
pbc = false
z0 = 1.0
b1 = [-1.0, 1.0, 1.0]
b2 = [1.0, -1.0, 1.0]

[material]
burgmag = 3e-10
mu = 160e9
nu = 0.31
a = 0.1
maxseg_ratio = 0.04
minseg_ratio = 0.01
rann = 0.0316

[forces]
force_mode = "LineTension"
Ec = 1.0e6
mobility_law = "SimpleGlide"

[timeint]
integrator = "EulerForward"
dt = 1.0e-9

[topology]
split_mode = "MaxDiss"
collision_mode = "Proximity"
remesh_rule = "LengthBased"
n_div = [4, 4, 4]

[loading]
loading_mode = "stress"
applied_stress = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

[run]
max_step = 200
print_freq = 10
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.1
+++

The initial configuration is made of two dislocation lines that intersect at a single point. Each line has three nodes: its two end nodes are pinned, and its middle node is free. This lets the two lines zip into a binary junction under zero stress.

The boundary conditions are not periodic.
