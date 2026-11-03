# Mobility-law tutorials

This directory contains three tutorials that progress from evaluating built-in
mobility laws, to fitting an existing mobility law to reference data, to implementing a custom PyDiS mobility law.

Run the examples from this directory so that their relative paths and input
files resolve correctly:

```bash
cd tutorials/02_mobility_law
```

## Tutorial 1: Calculate stress-velocity curves

**File:** `test_stress_velocity.py`

This tutorial calculates the velocity of an infinite straight dislocation over a
range of applied shear stresses. It compares several dislocation character
angles and plots the resulting stress-velocity curves. Linear fits provide an
estimate of the effective mobility.

The script demonstrates how to:

- construct a straight periodic dislocation;
- select and configure a built-in mobility law;
- calculate nodal forces and velocities; and
- compare screw, mixed, and edge dislocation behavior.

Run it with:

```bash
python3 test_stress_velocity.py
```

The mobility law, its parameters, character angles, and stress range can be
changed near the bottom of the script.

## Tutorial 2: Fit a mobility law to reference data

**Files:**

- `mobility_data.json` contains stress-velocity MD data for screw and edge
  dislocations at several temperatures.
- `plot_mobility_data.py` loads the dataset and plots its temperature-dependent
  edge and screw responses.
- `fit_mobility.py` selects a built-in ExaDiS mobility law, chooses which
  parameters may vary, fits the model to the reference data, and compares the
  fitted curves with the data. Fitted mobility parameters can the be tested
  by simulating the evolution of a dislocation glissile loop under shear.

First inspect the input data:

```bash
python3 plot_mobility_data.py
```

Then run the fit:

```bash
python3 fit_mobility.py
```

In `fit_mobility.py`, change `mobility_law` and `dislocation_type` to explore
other model and dataset combinations. The `MOBILITY_PARAMS` dictionary defines
the available fitting profiles, initial values, parameter bounds, and flags.
The supporting `lm_fitter.py` module supplies the nonlinear least-squares
fitter used by this tutorial.

## Tutorial 3: Implement a custom PyDiS mobility law

**Files:**

- `user_mob_bccglide.py` implements `My_MobilityLaw`, a character-dependent BCC
  glide mobility law, along with a compatible custom topology handler.
- `test_mob_frank_read_src_pydis.py` applies the custom law to a Frank-Read
  source simulation.

Read `user_mob_bccglide.py` first to see how a user mobility law derives from
`MobilityLaw_Base` and provides whole-network and single-node mobility methods.
Then run the example with:

```bash
python3 test_mob_frank_read_src_pydis.py
```

The simulation visualizes the evolving network and writes its final state to
`output/frank_read_src_pydis_final.json`.
