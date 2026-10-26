# 04_fcc_slip: plastic strain on an FCC slip system

A Frank-Read source on slip system A3 = (-1,1,1)[101]/2 under uniaxial
loading along [001], with the plastic strain that dislocation motion produces
and the constant-strain-rate stress update that consumes it.

pydis `SimulateNetwork` computes both, with the arguments of ExaDiS'
driver: `loading_mode='strain_rate'` with `erate` and `edir`, or `'stress'`
and `'stress_rate'` (`srate`), and `rotation`, `max_strain` and `burgmag`.
The plastic strain is computed in every mode and written, with the stress
and the density, to `stress_strain_dens.dat` when `print_freq` is set.
`sim_uniaxial.py` only adds the per-step history the tests read.

```
make PYTHON=/path/to/python3     # every test but the slow one, about a minute
make leak                        # test_leak_measurement.py, about 8 minutes
```

The three tests named `test_*_pydis_exadis.py` need `pyexadis`, built with
the rest of OpenDiS. The other tests need pydis only. The merge and split
tests among them need ExaDiS `32082bb` or later. An older ExaDiS miscounts
the plastic strain they check.

## Tests

| File | What it checks |
|---|---|
| `test_fcc_slip_sys.py` | `pydis.FCC_slip_sys`: the 12 systems, lookups both ways, Schmid factors, A3 on [001] = 1/sqrt(6) |
| `test_frank_read_src_A3_pydis.py` | the source at a constant 400 MPa: it operates, and every moving segment stays on A3 |
| `test_plastic_strain_analytic.py` | the swept-area formula against a closed form, across box faces, under sign reversal, and against MLmat's own implementation; the loading-axis rotation |
| `test_uniaxial_A3_pydis.py` | constant strain rate, 1e3 /s: sigma = E (eps - eps_p) every step, each step's dEp a pure A3 glide increment, the ExaDiS output format |
| `test_loading_modes_pydis.py` | `SimulateNetwork`'s `'stress'` and `'stress_rate'` modes, `max_strain`, and `stress_strain_dens.dat`, on the source with every node pinned |
| `test_plastic_strain_pydis_exadis.py` | against ExaDiS: one prescribed step through PyDiS, ExaDiS' python and C++ plastic strain; the A3 strain-rate run through both drivers, step by step, with the loading axis fixed and turned |
| `test_merge_plastic_strain_pydis_exadis.py` | a merge of two connected nodes at a point off their segment's line, by ExaDiS' own collision. ExaDiS and PyDiS both count that segment's swept area once. ExaDiS before `32082bb` counted it twice |
| `test_split_plastic_strain_pydis_exadis.py` | a split that creates a segment, its two products moved apart off one line, replayed through ExaDiS' own `execute_split`. ExaDiS and PyDiS both count the area that segment sweeps. ExaDiS before `32082bb` left it out |
| `test_topological_plastic_strain.py` | plastic strain swept by node insertion, removal, merge, and split; a loop annihilated by remesh sweeps its whole area |
| `test_leak_measurement.py` | the A3 run with the topological accounting on and off |

## Conventions

- Voigt order `[xx, yy, zz, yz, xz, xy]`; spins as `[yz, xz, xy]`.
- A segment 1->2 with Burgers vector b moving from `(r3, r4)` to `(r1, r2)`
  contributes `P = outer(0.5 (r2 - r3) x (r1 - r4), b) / V`, with
  `dEp = 0.5 (P + P^T)` and `dWp = 0.5 (P - P^T)`. This is ParaDiS
  DeltaPlasticStrain and ExaDiS System::plastic_strain. The spin sign is
  ExaDiS'; ParaDiS and MLmat store `-dWp` and flip their rotation angles to
  match.
- The stress update is ParaDiS LoadCurve, loadType 1:
  `dsigma = E (erate dt - e.dEp.e)` along the loading axis e.
- `rotation=True` turns the axis and the applied stress by the plastic spin,
  as ExaDiS does. ParaDiS turns only the axis. Off by default.

## Topological operations

Collision, topology, and remesh move and remove nodes outside time
integration. Without accounting for that, the area those moves sweep is lost.
`SimulateNetwork` keeps a `pydis.util.plastic_strain.PlasticStrainAccum`
in `state["plast_acc"]`, and the `DisNet` operations add to it:
`move_node`, `insert_node_between`, `remove_two_arm_node`, `merge_node`, and
`split_node` all take `plast_acc`. Left `None`, the default, they behave
exactly as before. Collision, topology, or remesh done by ExaDiS components
inside the pydis driver are not accounted.

An operation that moves several nodes is accounted as one move, each segment
counted once with both ends moving. ExaDiS accounts the move of each node in
turn. Since `32082bb` it gives the same plastic strain. Before that, it
counted the collapsing segment of a merge twice. It also left out the new
segment of a split. ExaDiS accounts no plastic spin for topological
operations. PyDiS does.

The check that nothing leaks is the identity `sum P = F(start) - F(now)`,
where `F = sum over segments of outer(0.5 r1 x r2, b) / V` is the network's
area functional (`plast_strain_util.network_area_dyad`). It holds to
round-off for a network that conserves Burgers vector at every node, until a
loop reaches its own periodic image: the A3 plane holds the lattice vector
(1000, 1000, 0), and the A3 source reaches its image near step 594. After
that F is not defined, so `test_leak_measurement.py` checks only before it.

## Notes

- PyDiS' `LineTension` force uses the core-energy line tension
  `mu / (4 pi) ln(a / 0.1)`, so this 500 b source operates from about
  100 to 130 MPa axial.
- One source in a 1000 b box cannot keep up with 1e3 /s: after 300 steps the
  plastic strain is about 1.3% of the total, and the curve sits about 5 MPa
  below the elastic line. A yield plateau needs a far lower rate and thousands
  of steps.
- The source bows from the first loaded step, so there is no stage where a
  free source gives exactly zero plastic strain; the purely elastic check
  uses a network with every node pinned.
- `fcc_slip_sys.py` and `disl_analysis_on_dict_util.py` are reference copies
  from the MLmat analysis code. pydis does not import them;
  `test_plastic_strain_analytic.py` imports the second to compare against it.
