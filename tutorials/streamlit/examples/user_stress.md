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
force_mode = "LineTension"
mobility_law = "FCC_0"
Medge = 64103.0
Mscrew = 64103.0
vmax = 4000.0
stress_field = true
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
tau = 3e6

[run]
max_step = 300
print_freq = 10
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.0001
+++

A single dislocation moves through the internal stress of a dislocation network, frozen in place. The network is the ParaDiS restart rs0084, of 29,576 nodes, in a periodic box of 58824 b, about 15 µm. ExaDiS computed its stress on a 128³ grid of the box. The moving line does not change this field.

The line lies on the fcc slip system B4: b along [-1 0 1], on (1 1 1). At first it is straight and screw, and it runs through a random point of the box (the seed). A shear stress τ on its slip system drives it, which step 6 sets. Alone, τ would move it at the drag velocity M b τ, 16.35 m/s per MPa. The frozen field slows it, the more so at a low stress.

The force module CalForceUserStress (step 3) adds the field's Peach-Koehler force to the line tension. It integrates the stress along each segment. It follows the tutorial [05_user_stress](https://github.com/OpenDiS/OpenDiS/tree/main/tutorials/05_user_stress). Step 3 can leave the field out: the line then moves at M b τ.

The script downloads the field once, a file of 55 MB, into tutorials/05_user_stress/. Step 3 shows its Dropbox link, which can be changed to another field's. That field must be on the same box, in the same layout, which the link's help gives. The history view plots the line's mean displacement against time: the area that the line has swept, from the plastic strain, over its initial length. It gives the line's velocity over the second half of the run, next to M b τ.

The state takes three of ExaDiS's defaults, so that PyDiS runs the same simulation: the glide planes, splitMultiNodeAlpha = 1e-3 and coarsen_mode = 1.
