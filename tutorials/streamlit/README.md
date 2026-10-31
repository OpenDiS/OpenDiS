# OpenDiS — Streamlit app

A web app that sets up, runs and shows the PyDiS and ExaDiS simulations of the
[OpenDiS tutorials](https://opendis.github.io/OpenDiS/tutorials/index.html), step by step, and
the networks that their initial configurations build. **Run** runs the tutorial's script,
with your values, as a separate job. Its screen output shows as it runs, and its dislocation
network in 3-D.

There are two ways to use the app:

1. **In the browser, at https://opendis.streamlit.app.** No installation is needed. How the
   app runs there is in [On Streamlit Community Cloud](#on-streamlit-community-cloud), at the
   end.
2. **On your own computer.** The sections from [Installation](#installation) to
   [Tests](#tests) describe this way.

## Installation

Create the conda environment `opendis-vis` once, and install the app's packages into it, from
the root of the OpenDiS repository:

```
conda create -n opendis-vis python=3.12
conda activate opendis-vis
pip install -r requirements.txt
```

Another environment with Python 3.12 works too, for example a venv:

```
python3.12 -m venv ~/venvs/opendis
source ~/venvs/opendis/bin/activate
pip install -r requirements.txt
```

The root [`requirements.txt`](../../requirements.txt) holds numpy, matplotlib, Streamlit,
Plotly and watchdog, and CMake on Linux. Streamlit Community Cloud installs the same
packages. The app is tested with Streamlit 1.61 and 1.66, and with Plotly 6.9 and 7.1. PyDiS
needs its git submodule `core/pydis/python/pydis/graph`. If a clone lacks it, the app fetches
it when it starts.

PyDiS runs in pure Python, but its C library, `lib/libpydis.so`, makes it faster, and the
scripts use it if it is there, as the tutorials do. If a clone lacks it, the app builds it
when it starts, once, in the background (`build_pydis_library()` in `common/paths.py`). It
builds it as CMake's `_repro` builds do, without CMake: gcc compiles the sources that
`core/pydis/c`'s `CMakeLists.txt` files list, with `-O3 -fPIC -ffp-contract=off` and the
portable log/atan shim. Then ctypesgen writes `lib/pydis_lib.py`; ctypesgen is the git
submodule `core/pydis/python/ctypesgen`, which the app fetches first if need be. Its runs then
agree byte for byte with those of a `_repro` build made with CMake. Until the library is
there, and if gcc is missing, PyDiS uses its Python versions, and its output starts with
"pydis_lib not found".

The Frank-Read source also runs with ExaDiS, through its Python module `pyexadis`. On Linux,
if a clone lacks the module for the app's Python, the app builds it after PyDiS's library, in
the background (`build_exadis()` in `common/paths.py`). It uses
OpenDiS's CMake project with its `ubuntu` preset, `-DSYS=ubuntu`. Two options differ from a
usual build. `-DEXADIS_FFT=Off` leaves out ExaDiS's FFT force, since the Cloud has no FFTW.
`-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF` turns off link-time optimization, whose last step
would compile all the files at once. The build compiles one file at a time, to fit in the
Cloud's memory, and takes some minutes. ExaDiS is the git submodule `core/exadis`, with Kokkos
and pybind11 as submodules of its own. The app fetches them first if need be, without their
history. The module goes to `core/exadis/python/`, where the scripts look for it. Until it is
there, an example run by ExaDiS cannot run: its Run is off. On macOS, build ExaDiS with CMake,
as for the tutorials.

## Launch

```
cd tutorials/streamlit
./launch_app.sh
```

The app is served at http://localhost:8505, or at the next free port if another program uses
that one. No browser window opens: open the address that the script prints. The ports from
8501, Streamlit's default, are left to other Streamlit apps. `launch_app.sh` finds the Python
of `opendis-vis` by itself, so that environment need not be activated. For another
environment, set `OPENDIS_PYTHON` to its Python, as in
`OPENDIS_PYTHON=~/venvs/opendis/bin/python ./launch_app.sh`.

## The page

- **Banner:** "OpenDiS" with the logo of the OpenDiS website, pinned at the top.
- **Sidebar:** the example, in a drop-down menu, and below it a second one, **Core**: PyDiS or
  ExaDiS. Every example has both. Then the
  steps, where a • marks a step whose values differ from the example's, and ‹ prev and next ›
  buttons. A simulation has seven steps; the initial configurations, which run no simulation,
  have two.
- **2nd left column:** the chosen step's selections. « hides it, » shows it again, and
  choosing a step shows it too. The binary junction's step 1 starts with a third menu,
  **Crystal structure**: BCC or FCC. Each is a sub-example, with its own values in every step,
  and choosing one goes to its step 1.
- **Center:** the 3-D window and the code, PyDiS's or ExaDiS's.
  - The 3-D window shows the initial configuration, as PyDiS builds it from your values
    whichever code runs the example; the
    simulation, which during a run is its newest snapshot, and after it a replay of all its
    snapshots, with ▶ play and a slider (it starts on the final network); and the history, the
    number of nodes and the dislocation density of each snapshot against the time step; for a
    strain-rate run, also the stress against the strain, from the run's
    `output/stress_strain_dens.dat`. For the initial configurations, the second view, **run**,
    shows the network that the run built, and there is no history.
  - The cube first shows its x axis along its front-left edge and its y axis along its
    front-right edge, seen a little from the side rather than along its diagonal. Turn it
    with the mouse: the view you turn to stays when another network of the same example is
    drawn, as when a value changes, a replay plays or you switch between the views; another
    example starts in the first view, which is the same for every example, and so does an
    example that you come back to. Plotly's "Reset camera to default" goes back to the first
    view. A replay holds while the
    mouse button is down, so that you can turn the view as it plays. So does a run: its newest
    snapshot keeps the view you turned to, and waits until the mouse button is up.
  - The networks are drawn by a small component of the app's own, `components/network3d/`,
    since Streamlit's own Plotly chart puts the camera back whenever it draws a new figure.
    During a run, the app writes the newest snapshot into the component's `live/` folder, and
    the component fetches it there every 0.5 s, so that one figure lasts the whole run.
  - The code panel shows this step's code, or the full script, with a download button.
  - Its third mode, **Python console**, is the browser tab's own Python session, to look into
    the run's objects. `net`, `sim`, `state` and `G` are as after `python -i` on the full
    script, and `snapshot(step)` loads one of the run's snapshots (`snapshot_steps` lists
    them). The session runs the run's script in a process of its own, with
    `SimulateNetwork.run` replaced by a load of the run's final network, so nothing is
    simulated again. It works on a copy of the run's results, so it changes neither the runs
    nor the other panels. Guard rails refuse what looking into the objects does not need:
    other programs, threads, the network, native libraries, and changes to files outside the
    session's folder. A command stops after 30 s, and the session ends after 15 min without
    one (`common/console.py`, `common/console_child.py`).
- **Right:** the simulation output: what the script prints, as it runs, and after a run its
  `log.txt`, which ends with a line on how the job ended. Below it, **🗑 clear all results**
  deletes the results of all past runs, after a confirmation: those of every example, for
  everyone who uses the app. A run in progress is kept. Below it is ExaDiS's status: ready,
  or a bar while the app builds it, which refreshes itself every 5 s.
- **Run (the last step)** runs the full script as a separate job, in its run folder. With
  ExaDiS, it is off until ExaDiS is built, and the step says why. The app stops a run only if
  you press **Stop** or change any selection during the run, if it prints nothing for 5 min,
  or if it runs for an hour. After a run, the step offers its results for download.

## Examples

An example, run by one code, is one Markdown file, `examples/<name>.md`. The Example menu offers
one entry per `label`, and the Core menu one per file of that label. The files of a label may
also be sub-examples, each with its `variant`; a menu at the top of step 1 picks one, and the
Core menu then offers its files. The file has two parts:

- **Front matter**, between two lines of `+++`, in TOML: the `label` that the menu shows, the
  `variant` if any (`BCC`) and the name of step 1's menu (`variant_label`: `Crystal
  structure`), the code that runs it (`engine`: `pydis`, the default, or `exadis`), the
  `title`, the script's file name (`script_name`) and that of its final network
  (`final_json`), the tutorial's script and page (`tutorial_script`, `tutorial_url`, and the
  link's text `tutorial_title` if not the default), a `run_hint`, and one table of values for
  each step: `[config]`, `[material]`, `[forces]`, `[timeint]`, `[topology]`, `[loading]` and
  `[run]` (see `common/spec.py`). A file that builds a network and runs no simulation has
  `[config]` only.
- **Body**, in Markdown: the example in words, which step 1 shows.

| Example | Tutorial's script | Initial configuration |
| --- | --- | --- |
| `examples/frank_read_src.md` | `examples/02_frank_read_src/test_frank_read_src_pydis.py` | `frank_read_src_loop` |
| `examples/frank_read_src_exadis.md` | `examples/02_frank_read_src/test_frank_read_src_exadis.py` | `frank_read_src_loop` |
| `examples/binary_junction.md` | `examples/03_binary_junction/test_binary_junction_pydis.py` | `two_disl_lines` |
| `examples/binary_junction_exadis.md` | `examples/03_binary_junction/test_binary_junction_exadis.py` | `two_disl_lines` |
| `examples/binary_junction_fcc.md` | `examples/03_binary_junction/test_fcc_junctions_pydis.py` | `fcc_junction` |
| `examples/binary_junction_fcc_exadis.md` | `examples/03_binary_junction/test_fcc_junctions_exadis.py` | `fcc_junction` |
| `examples/uniaxial_tension.md` | `examples/05_uniaxial_tension/test_uniaxial_tension_pydis.py` | `loop_sources` |
| `examples/uniaxial_tension_exadis.md` | `examples/05_uniaxial_tension/test_uniaxial_tension_exadis.py` | `loop_sources` |
| `examples/initial_configurations.md` | `tutorials/01_initial_configurations/initial_configurations.py` | `simple_source` |
| `examples/initial_configurations_exadis.md` | the same | `simple_source` |

The binary junction has two crystal structures. BCC is the tutorial's two lines, now run by
ExaDiS too. FCC follows ExaDiS's example `05_fcc_junctions`: two Frank-Read sources, each as
long as a periodic box of 1000 b, on two slip systems that cross near its centre. Step 1 offers
its four junction types (lomer, glissile, hirth and colinear) and the lines' two angles. Its two
tutorials in `examples/03_binary_junction/` are new.

ExaDiS's example uses `DDD_FFT_MODEL`. PyDiS lacks it, and the Cloud's ExaDiS, built without
FFT, does too. So the segments interact within a cutoff of 300 b, with PyDiS's `Elasticity_SBA`
or ExaDiS's `CUTOFF_MODEL`, which select the same segment pairs. PyDiS had no `FCC_0`. It is
now in `core/pydis/python/pydis/mobility/mobility_disnet.py`, translated from ParaDiS's
`MobilityLaw_FCC_0.c`, with ExaDiS's speed cap. PyDiS had no `insert_frank_read_src` either:
`core/pydis/python/pydis/util/network_util.py` is a port of ExaDiS's.

PyDiS's state takes three of ExaDiS's defaults, which the app offers in steps 2 and 5:
`use_glide_planes`, `splitMultiNodeAlpha` and `coarsen_mode` = 1. `splitMultiNodeAlpha` = 1e-3
is a factor that ExaDiS puts on the noise floor of its split test; PyDiS's `Serial` split now
reads it too, and its default, 1, is ParaDiS's. With these, the lomer and glissile junctions
agree with ExaDiS for all 200 steps, node for node. The hirth junction parts at step 149, and
the colinear one at step 7: ExaDiS keeps a collision and a remesh to the glide planes, and
PyDiS does not. PyDiS's 200 steps take about 30 s, and ExaDiS's about 3 s. ExaDiS warns that
the Burgers vector is not conserved at 4 nodes: they are the lines' pinned ends.

The uniaxial tension is the setup of `uniaxial_strain_rate.py`: four prismatic loops in a
periodic bcc crystal, pulled along z at a constant strain rate. It follows ExaDiS's example
`21_bcc_Ta_100nm_2e8`, with a larger box, fewer loops and a much lower strain rate. Its two
tutorials in `examples/05_uniaxial_tension/` run it with each core. It takes options that the
other examples do not: the crystal, `rtol` and `nextdt` in step 2; the mobility law `BCC_0B`,
with `Medge`, `Mscrew`, `Mclimb` and `vmax`; the integrator `Trapezoid`, whose time step
adapts, with `multi` steps in a row; and the loading `strain_rate`, with `erate` and `edir`. PyDiS had neither `BCC_0B` nor `Trapezoid`. They are
now in `core/pydis/python/pydis/mobility/mobility_disnet.py` and `timeint/timeint_disnet.py`,
translated from ParaDiS's `MobilityLaw_BCC_0b.c` and `TrapezoidIntegrator.c`. Compared with
ExaDiS, they follow it where ExaDiS adds to ParaDiS: the junction nodes and the speed cap of
`BCC_0B`, the prediction from the current velocities, and the steps in a row. PyDiS splits the
multi-nodes in its `Serial` mode, which the app now offers next to `MaxDiss`. With `MaxDiss`,
the two new nodes stay on top of each other, a junction node splits and merges again at every
step, and the time step stalls near 1e-11 s. PyDiS runs 50 steps by default, in about 15 s;
its 500 steps take about 20 minutes, ExaDiS's about 11 s.

The initial configurations follow the tutorial's script, which builds a network in one of five
ways and plots it. Its options 1 (PyDiS) and 2 (ExaDiS) are one network here, the simple
dislocation source; the Core menu picks the code. Step 1's **Network** menu offers the four:
`simple_source`, `paradis_file` (option 3), `line_config` (option 4) and `prismatic_config`
(option 5). Another network starts from the tutorial's values. Each runs with either core:
options 3 to 5 used ExaDiS only, so PyDiS's `DisNet` has ExaDisNet's `read_paradis`,
`write_data`, `generate_line_config` and `generate_prismatic_config` too, and `write_vtk`
(`core/pydis/python/pydis/disnet.py`). They are ports of ExaDiS's code, with the same random
numbers: for the same values they give ExaDiS's network, and the same `.data` and `.vtk`
files. The ParaDiS file is one of the repository's, from a fixed list, not any path on the
server.

Another example of a known initial configuration is another file, and the app's code does not
change; a new or edited file shows after you reload the page. A new kind of initial
configuration needs its function in `common/codegen.py` (`INIT_FUNCTIONS`, one for each code,
and `KIND_IMPORTS` for what it imports) and its values in `common/spec.py` (`CONFIGS`,
`ENGINE_KINDS`). Each code offers its own
choices in the steps (`CHOICES`). ExaDiS has a topology mode where PyDiS has a split mode, and
no cell list. A file with problems is not offered: the app lists what is wrong
with it, under **Input files with problems**.

## The generated script

The full script is the tutorial's script with your values. For an example's own values, it
runs the same simulation as the tutorial's script, to the last bit (`tests/test_codegen.py`).
It differs from that script in four ways. It finds PyDiS and ExaDiS from `$OPENDIS_DIR`, else
from the repository that the app runs in, so that it runs from any folder. Each step's code is
one block, so a few lines come in another order. An optional value (Ec, mob, a Burgers vector)
appears only when it is not the default. And ExaDiS's script also writes its snapshots as
PyDiS's does, `output/disnet_<step>.json`, for the app's views (`SimulateNetworkJSON`).

The script of the initial configurations is the tutorial's option for that network and core.
It builds the network, plots it and writes it as its final network; its other files go into
`output/`. The simple source's Burgers vector is a direction, which the script makes of unit
length: `[1, 1, 0]` gives the tutorial's `1/sqrt(2)*[1,1,0]` to the last bit.

You can download it and run it outside the app, with `python <script name>`: it plots with
matplotlib, as the tutorials do, and writes `output/` in the current folder. In the app,
matplotlib draws off screen, and the app shows its own 3-D view.

## Data and outputs

| Path | Role | Git |
| --- | --- | --- |
| `app.py`, `common/` | the app | committed |
| `examples/<name>.md` | one example: its input file | committed |
| `assets/opendis_logo.png` | the banner's logo, scaled from `doc/nodebox_logo.png` (the `document` branch) | committed |
| `assets/opendis_favicon.png` | the browser tab's icon, the website's favicon, scaled from `doc/_static/favicon_node.ico` | committed |
| `components/network3d/` | the 3-D window's component: its page and `network3d.js` | committed |
| `components/network3d/plotly.min.js` | plotly.js, which the app copies there from the Plotly package when it starts | ignored |
| `components/network3d/live/` | a running job's newest snapshot, `<run_id>.json`, for the 3-D window; removed when the run ends | ignored |
| `launch_app.sh` | the launcher | committed |
| `tests/` | the tests, and `record_demo.py`, which records `recordings/app_demo.mp4`, a demo video of the app | committed |
| `runs/<example>/<run_id>/` | one folder per run: the script, `spec.json`, `log.txt`, and `output/`, with a snapshot of the network every write_freq steps (`disnet_<step>.json`, and ExaDiS's `config.<step>.data`), the final network, and any other file that the script writes there | ignored |
| `$TMPDIR/opendis_console_*/` | a Python console's session: its copy of the run's script and `output/`, removed when the session ends | outside the repository |

`run_id` is a hash of the script, so the same selections reuse their results without running
again. To run them again, delete their run folder, or press **🗑 clear all results**.
`OPENDIS_APP_RUNS` points the app at another
runs folder, and `OPENDIS_DIR` at another OpenDiS repository.

## Launch details

- `launch_app.sh` runs the app with `OPENDIS_PYTHON` if set, else
  `/opt/anaconda3/envs/opendis-vis/bin/python` if that file exists, else the environment that
  `conda` knows as `opendis-vis`, else the `python3` on the PATH. It serves the app on port 8505,
  or the next free port up to 8599; `STREAMLIT_SERVER_PORT` sets the port instead (and the
  launch fails if that port is taken), and `STREAMLIT_SERVER_ADDRESS` the address.
- Instead of the script, you can run `python -m streamlit run app.py --server.headless true`
  in the activated environment.
- **After you edit the code, restart the app** (Ctrl-C, then launch again). Reloading the tab
  is not enough.
- **The 3-D window needs WebGL.** If it says "WebGL is not supported", turn on hardware
  acceleration in the browser.

## Tests

Run each test from `tutorials/streamlit/`, with the environment's Python:

```
python -m tests.test_codegen          # the generated scripts against the tutorials' own (60 steps)
python -m tests.test_codegen --full   # ... at the tutorials' full length (about a minute)
python -m tests.test_job              # the job runner: a run, a stop, a time limit, leftovers
python -m tests.test_views            # networks from PyDiS's JSON, and the figures
python -m tests.test_app_smoke        # the whole page, and 20-step runs with PyDiS and ExaDiS
python -m tests.test_console          # the Python console, on a run of each example
python -m tests.test_build            # PyDiS's library and ExaDiS's module (some minutes)
node tests/test_network3d.js          # the 3-D window's script (needs Node.js)
```

Each test prints `PASS` or `FAIL` for every check, and exits with a non-zero code on a failure.

- **`test_codegen`**: for each example, the script that the app writes for the example's own
  values must run the same simulation as the tutorial's script: every snapshot of the network,
  and the final one, equal to the last bit, and the same printed lines but for the wall time.
  For ExaDiS, the snapshots are its `config.<step>.data` files, and the app's script must also
  write a JSON snapshot for each. The tutorial's script runs from a copy whose only change is
  that its paths are absolute (they are relative to its own folder). Also checked: the numbers
  in the scripts, each step's code inside the full script, the spec's JSON round trip, and
  selections that must be refused (a zero time step, minseg over half of maxseg, a NaN stress,
  ...). For the initial configurations, each of the four networks with each core: the app's
  job and the tutorial's option for it give the same network, the 3-D window's own build is
  that network, and the infinite lines' `.data` and `.vtk` files are the tutorial's. The
  tutorial's plot is replaced by a write of its network. Its options 3 to 5 run ExaDiS, so
  PyDiS's network is checked against ExaDiS's.
- **`test_job`**: a run that ends well, one that fails, one that writes no final network, a
  stop from the caller (as Streamlit stops a script run), a time limit, a job left by an
  earlier session, the run folders, and their deletion, which keeps a folder whose job runs;
  each with a small stand-in script in a temporary folder, and none left running at the end.
- **`test_views`**: the initial configurations, as PyDiS builds them; a network written by
  `DisNetManager.write_json` and read back without PyDiS, across a periodic boundary; the
  periodic images, against PyDiS's own `Cell.closest_image`; the figures (their traces, the
  first view, the camera kept, the stress against the strain of a strain-rate run); and a run's
  results (300 snapshots played as 120, the history, the zip).
- **`test_app_smoke`**: the banner, the two menus in the sidebar, no chat, ExaDiS's status,
  every step of each example, sub-example and core, the binary junction's Crystal structure
  menu, a change and a reset, and a 20-step
  run of the binary junction: its output, status and downloads, and the networks that the 3-D
  window's component gets (AppTest does not run the component's script); then the Python
  console, and **🗑 clear all results**. Then a 20-step run of the Frank-Read source with
  ExaDiS, its Run off while ExaDiS is not built, and the build's status bar. Then the initial
  configurations: their two steps, each network drawn in step 1, and a run that builds the
  prismatic loops with ExaDiS. Its runs go to a temporary folder (`OPENDIS_APP_RUNS`).
- **`test_console`**: for each example and code, a 20-step run (the initial configurations
  build their network), then a session whose `G`, `sim` and `net` are the run's (no `sim`
  for a build); then, in one session on the binary junction's run: `state`, code of several
  lines, `print`, an error's traceback, names that stay, `snapshot()`, a file written to the
  session's own copy, the app's environment kept out, `input()`, the guard rails (shell commands, subprocess,
  sockets, threads, signals, native libraries, writes and removals outside, the repository's
  git files, and a replaced `os.path.abspath`), PyDiS's own computations, `import pyexadis`
  where ExaDiS is built for this Python, the limit on sessions, and a command that runs too
  long, which ends the session and removes its folder.
- **`test_build`**: PyDiS's C library, built as the app builds it, into a copy of the
  repository's Python code; a 20-step run of the binary junction with it; and, if the
  repository has its own library, the same run with that one, which must agree byte for byte
  (as it does when that library is a `_repro` build). Then ExaDiS's module, built as the app
  builds it, into a temporary folder; on macOS without OpenMP, which Apple's clang lacks. With
  it, the binary junction's ExaDiS tutorial runs as it is, but for its plots.
- **`test_network3d.js`**: the component's script in Node.js, with stand-ins for Plotly, the
  page and the browser: a run starts on its final network in the first view; the camera that
  Plotly reports after a drag is the one every later network of the example is drawn with, in
  its next view too, while another example, or a return to the first, starts in the first view;
  play, the slider, a hold while the mouse button is down, the end, and new networks; and
  during a run, the newest snapshot fetched and drawn with the user's camera, held during a
  drag, and the fetching stopped when the replay comes.

## On Streamlit Community Cloud

The app runs at https://opendis.streamlit.app, from the `main` branch of
[OpenDiS/OpenDiS](https://github.com/OpenDiS/OpenDiS), with `tutorials/streamlit/app.py` as its
entry point.

- The Cloud installs the packages of the root `requirements.txt`, with uv. It looks for a
  dependency file next to `app.py` first, and an `environment.yml` there would make it use
  conda. So the app's folder has none.
- The Cloud clones the repository without its git submodules. So the app fetches PyDiS's
  `graph` submodule when it starts, with `git submodule update --init` (`fetch_graph()` in
  `common/paths.py`).
- All visitors share the app's `runs/` folder. So a run's results show for everyone who picks
  the same selections, until **🗑 clear all results** or a reboot of the app.
- The Cloud's clone has no `lib/libpydis.so`, so the app builds it when it starts. The Cloud
  has gcc, g++ and make, and CMake only from pip. PyDiS's library takes some seconds. Runs
  that start before it is done use PyDiS's Python versions. ExaDiS's module comes next, in
  some minutes, one file at a time. Until then, ExaDiS's Run is off, and the status bar below
  **🗑 clear all results** shows the build. The app's log, under Manage app, says when each
  was built, or why not. The Cloud runs Python 3.14.
- A push updates the app in place, and keeps what the app built. After a change to PyDiS's C
  code or to ExaDiS, reboot the app, so that it builds them anew.
- The Python console works there too, for every visitor. It is meant for looking into a
  run's objects. Each session has its own process, its own copy of the results and an
  environment without the app's variables. It has the guard rails above, limits on time,
  CPU, memory, file size and open files, and at most 4 sessions run at once. It is not a
  sandbox against a determined expert: Python's audit hooks cannot be one, and code typed
  there runs on the app's machine.
