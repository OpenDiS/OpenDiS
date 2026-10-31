+++
label = "Initial configurations"
engine = "exadis"
title = "Initial configurations · ExaDiS"
order = 300
script_name = "initial_configurations_exadis.py"
final_json = "initial_configurations_exadis_final.json"
tutorial_script = "tutorials/01_initial_configurations/initial_configurations.py"
tutorial_title = "Creating initial dislocation configurations"
tutorial_url = "https://opendis.github.io/OpenDiS/tutorials/initial_configuration.html"
run_hint = "Run builds this network with ExaDiS, and writes it. It runs no simulation."

[config]
kind = "simple_source"
Ldis = 100.0
burg = [1.0, 1.0, 0.0]
plane = [-1.0, 1.0, 1.0]
+++

The tutorial's script builds a dislocation network in one of five ways, and plots it. Its first two are one here: the Core menu runs them with PyDiS or with ExaDiS. Each way runs with either core.

- **simple dislocation source**: a straight screw dislocation along its Burgers vector, pinned at its two ends, with a free node in its middle (the tutorial's options 1 and 2).
- **ParaDiS file**: a network read from a ParaDiS data file of the repository (option 3).
- **infinite lines**: straight lines through the periodic box, on the crystal's slip systems (option 4). The script also writes them as `output/infinite_lines.data` and `output/infinite_lines.vtk`.
- **prismatic loops**: loops of the crystal's Burgers vectors, spread through the box (option 5).

The script runs no simulation: Run builds the network with the chosen core, and writes it.
