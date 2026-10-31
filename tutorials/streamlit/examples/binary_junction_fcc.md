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

Two Frank-Read sources lie on two fcc slip systems, in a periodic box of 1000 b. Each is as long as the box and pinned at its two ends. The first has b1 along [0 1 1], on (1 1 -1); the junction type sets the second. The lines cross near the box's centre, at the angles phi1 and phi2 from the intersection of their planes. Under zero stress they zip into a binary junction: the lomer one, b1 + b2 along [1 1 0], lies off the {111} planes and can only move along its line.

The setup is that of ExaDiS's example 05_fcc_junctions. Its force model, DDD_FFT_MODEL, is missing from PyDiS, so the segments here interact within a cutoff of 300 b, as in ExaDiS's CUTOFF_MODEL. The mobility law is FCC_0, and the time step adapts (Trapezoid).

The state takes three of ExaDiS's defaults, so that PyDiS runs the same simulation: the glide planes, splitMultiNodeAlpha = 1e-3 and coarsen_mode = 1. The lomer and glissile junctions then agree with ExaDiS for all 200 steps.
