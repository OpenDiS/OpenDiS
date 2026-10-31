+++
label = "Uniaxial tension"
title = "Uniaxial tension · PyDiS"
order = 250
script_name = "uniaxial_tension_pydis.py"
final_json = "uniaxial_tension_pydis_final.json"
tutorial_script = "examples/05_uniaxial_tension/test_uniaxial_tension_pydis.py"
run_hint = "The 50 steps take about 15 s with PyDiS, and 500 about 20 minutes."

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
split_mode = "Serial"
collision_mode = "Retroactive"
remesh_rule = "LengthBased"
n_div = [8, 8, 8]

[loading]
loading_mode = "strain_rate"
erate = 3e3
edir = [0.0, 0.0, 1.0]

[run]
max_step = 50
print_freq = 1
write_freq = 10
plot_freq = 10
plot_pause_seconds = 0.01
+++

A bcc crystal of tantalum's constants holds four prismatic loops, of radius 0.21 Lbox, in a periodic box of 10000 b, about 2.9 µm. It is pulled along z at a constant strain rate of 3000 /s. The loops expand on their glide cylinders, and the applied stress follows from the plastic strain that they make.

The mobility law is BCC_0B: an edge and a screw mobility of 2600 /(Pa s), a climb mobility of 1e-4 /(Pa s), and a speed capped near 3400 m/s. The time step adapts: the Trapezoid integrator starts from nextdt and keeps the error below rtol, two steps at a time (multi = 2). The split mode is Serial, which moves the two new nodes apart, as ParaDiS and ExaDiS do. With MaxDiss, a junction node splits and merges again at every step, and the time step stalls.

The setup is that of uniaxial_strain_rate.py. It follows ExaDiS's example 21_bcc_Ta_100nm_2e8, with a larger box, fewer loops and a much lower strain rate. The history view plots the stress against the strain.
