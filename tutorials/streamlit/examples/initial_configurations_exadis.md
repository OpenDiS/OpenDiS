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

This example illustrates several ways to build an initial dislocation configuration. It runs no simulation.
