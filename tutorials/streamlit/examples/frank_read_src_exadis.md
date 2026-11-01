+++
label = "Frank-Read src"
engine = "exadis"
title = "Frank-Read source · ExaDiS"
order = 100
script_name = "frank_read_src_exadis.py"
final_json = "frank_read_src_exadis_final.json"
tutorial_script = "examples/02_frank_read_src/test_frank_read_src_exadis.py"
tutorial_url = "https://opendis.github.io/OpenDiS/tutorials/frank_read_src/frank_read_src_by_pyexadis.html"
run_hint = "The tutorial's 300 steps take about 5 s."

[config]
kind = "frank_read_src_loop"
box_length = 1000.0
pbc = true
arm_ratio = 0.125
burg_vec = [1.0, 0.0, 0.0]

[material]
burgmag = 3e-10
mu = 50e9
nu = 0.3
a = 1.0
maxseg_ratio = 0.04
minseg_ratio = 0.01
rann = 3.0

[forces]
force_mode = "LineTension"
mobility_law = "SimpleGlide"

[timeint]
integrator = "EulerForward"
dt = 1.0e-8

[topology]
topology_mode = "None"
collision_mode = "Retroactive"
remesh_rule = "LengthBased"

[loading]
loading_mode = "stress"
applied_stress = [0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]

[run]
max_step = 300
print_freq = 10
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.01
+++

A prismatic loop is pinned at its four corners. Its top arm bows out under the applied stress, as a Frank-Read source.
