+++
label = "User stress"
title = "User stress · PyDiS"
order = 280
script_name = "user_stress_pydis.py"
final_json = "user_stress_pydis_final.json"
tutorial_script = "tutorials/05_user_stress/example_internal_stress_pydis.py"
run_hint = "The 300 steps take about a minute with PyDiS. They cover about 17 ns."

[config]
kind = "user_stress_line"
box_length = 58824.0
theta = 0.0
seed = 0

[material]
burgmag = 2.55e-10
mu = 54.6e9
nu = 0.324
a = 6.0
maxseg_ratio = 0.017
minseg_ratio = 0.005
rann = 10.0
crystal = "fcc"
rtol = 10.0
nextdt = 1e-11
maxdt = 1e-9
use_glide_planes = true

[forces]
force_mode = "Elasticity_SBA"
cutoff = 5000.0
mobility_law = "FCC_0"
Medge = 64103.0
Mscrew = 64103.0
vmax = 4000.0
stress_field = true
stress_kind = "frozen"
stress_url = "https://www.dropbox.com/scl/fi/q8025nk18zrovpfupo1vt/rs0084.stress.npz?rlkey=z1ertzr41fn68c6sjzptrysyq&st=5kvnt6qi&dl=1"

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
tau = 5e6

[run]
max_step = 300
print_freq = 10
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.0001
+++

This example illustrates how to add a user-defined stress field to a DDD simulation, as in [user_stress.py](https://github.com/OpenDiS/OpenDiS/blob/main/tutorials/05_user_stress/user_stress.py). A dislocation moves through that field, driven by the applied shear stress τ.
