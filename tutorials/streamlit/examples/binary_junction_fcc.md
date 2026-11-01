+++
label = "Binary junction"
variant = "FCC"
variant_label = "Crystal structure"
title = "Binary junction · FCC · PyDiS"
order = 200
script_name = "fcc_junctions_pydis.py"
final_json = "fcc_junctions_pydis_final.json"
tutorial_script = "examples/03_binary_junction/test_fcc_junctions_pydis.py"
run_hint = "The tutorial's 200 steps take about 30 s."

[config]
kind = "fcc_junction"
box_length = 1000.0
junction_type = "lomer"
phi1 = 30.0
phi2 = 30.0

[material]
burgmag = 2.55e-10
mu = 54.6e9
nu = 0.324
a = 6.0
maxseg_ratio = 0.1
minseg_ratio = 0.02
rann = 10.0
crystal = "fcc"
rtol = 5.0
nextdt = 1e-12
maxdt = 1e-10
use_glide_planes = true

[forces]
force_mode = "Elasticity_SBA"
cutoff = 300.0
mobility_law = "FCC_0"
Medge = 64103.0
Mscrew = 64103.0
vmax = 4000.0

[timeint]
integrator = "Trapezoid"

[topology]
split_mode = "Serial"
collision_mode = "Retroactive"
remesh_rule = "LengthBased"
n_div = [8, 8, 8]
split_alpha = 1e-3
coarsen_mode = 1

[loading]
loading_mode = "stress"
applied_stress = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

[run]
max_step = 200
print_freq = 1
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.0001
+++

Two lines lie on two fcc slip systems, each pinned at its two ends. Under zero stress, they zip into a binary junction of the chosen type. It follows ExaDiS's example [05_fcc_junctions](https://github.com/llnl/exadis/blob/main/examples/05_fcc_junctions/test_fcc_junctions.py), with a cutoff of 300 b instead of its FFT model.
