# OpenDiS — Streamlit app

A web app that sets up, runs and shows the PyDiS simulations of the
[OpenDiS tutorials](https://opendis.github.io/OpenDiS/tutorials/index.html), step by step.
**Run** runs the tutorial's script, with your values, as a separate PyDiS job. Its screen
output shows as it runs, and its dislocation network in 3-D.

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
Plotly and watchdog. Streamlit Community Cloud installs the same packages. The app is tested
with Streamlit 1.61 and 1.66, and with Plotly 6.9 and 7.1. PyDiS is pure Python and needs no
build; if `lib/libpydis.so` has been built, the scripts use it, as the tutorials do. PyDiS
needs its git submodule `core/pydis/python/pydis/graph`. If a clone lacks it, the app fetches
it when it starts.

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

- **Banner:** "OpenDiS" with the logo of the OpenDiS website, pinned at the top, and on its
  right the example, in one drop-down menu.
- **Sidebar:** the steps, where a • marks a step whose values differ from the example's, and
  ‹ prev and next › buttons.
- **2nd left column:** the chosen step's selections. « hides it, » shows it again, and
  choosing a step shows it too.
- **Center:** the 3-D window and the PyDiS code.
  - The 3-D window shows the initial configuration, as PyDiS builds it from your values; the
    simulation, which during a run is its newest snapshot, and after it a replay of all its
    snapshots, with ▶ play and a slider (it starts on the final network); and the history, the
    number of nodes and the dislocation density of each snapshot against the time step.
  - The cube first shows its x axis along its front-left edge and its y axis along its
    front-right edge, seen a little from the side rather than along its diagonal. Turn it
    with the mouse: the view you turn to stays when another network of the same example is
    drawn, as when a value changes, a replay plays or you switch between the views; another
    example starts in the first view, which is the same for every example, and so does an
    example that you come back to. Plotly's "Reset camera to default" goes back to the first
    view. A replay holds while the
    mouse button is down, so that you can turn the view as it plays. (During a run, its newest
    snapshot starts again in the first view.)
  - The networks are drawn by a small component of the app's own, `components/network3d/`,
    since Streamlit's own Plotly chart puts the camera back whenever it draws a new figure.
  - The code panel shows this step's code, or the full script, with a download button.
- **Right:** the simulation output: what the script prints, as it runs, and after a run its
  `log.txt`, which ends with a line on how the job ended. Below it, **🗑 clear all results**
  deletes the results of all past runs, after a confirmation: those of every example, for
  everyone who uses the app. A run in progress is kept.
- **Run (step 7)** runs the full script as a separate job, in its run folder. The app stops it
  only if you press **Stop** or change any selection during the run, if it prints nothing for
  5 min, or if it runs for an hour. After a run, step 7 offers its results for download.

## Examples

An example is one Markdown file, `examples/<name>.md`; the drop-down menu offers one entry per
file. The file has two parts:

- **Front matter**, between two lines of `+++`, in TOML: the `label` that the menu shows, the
  `title`, the script's file name (`script_name`) and that of its final network
  (`final_json`), the tutorial's script and page (`tutorial_script`, `tutorial_url`), a
  `run_hint`, and one table of values for each step: `[config]`, `[material]`, `[forces]`,
  `[timeint]`, `[topology]`, `[loading]` and `[run]` (see `common/spec.py`).
- **Body**, in Markdown: the example in words, which step 1 shows.

| Example | Tutorial's script | Initial configuration |
| --- | --- | --- |
| `examples/frank_read_src.md` | `examples/02_frank_read_src/test_frank_read_src_pydis.py` | `frank_read_src_loop` |
| `examples/binary_junction.md` | `examples/03_binary_junction/test_binary_junction_pydis.py` | `two_disl_lines` |

Another example of a known initial configuration is another file, and the app's code does not
change; a new or edited file shows after you reload the page. A new kind of initial
configuration needs its function in `common/codegen.py` (`INIT_FUNCTIONS`) and its values in
`common/spec.py` (`CONFIGS`). A file with problems is not offered: the app lists what is wrong
with it, under **Input files with problems**.

## The generated script

The full script is the tutorial's script with your values. For an example's own values, it
runs the same simulation as the tutorial's script, to the last bit (`tests/test_codegen.py`).
It differs from that script in three ways: it finds PyDiS from `$OPENDIS_DIR`, else from the
repository that the app runs in, so that it runs from any folder; each step's code is one
block, so a few lines come in another order; and an optional value (Ec, mob, a Burgers vector)
appears only when it is not the default.

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
| `launch_app.sh` | the launcher | committed |
| `tests/` | the tests, and `record_demo.py`, which records `recordings/app_demo.mp4`, a demo video of the app | committed |
| `runs/<example>/<run_id>/` | one folder per run: the script, `spec.json`, `log.txt`, and `output/`, with a snapshot of the network every write_freq steps (`disnet_<step>.json`), the final network, and any other file that the script writes there | ignored |

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
python -m tests.test_app_smoke        # the whole page, and a 20-step run of the binary junction
node tests/test_network3d.js          # the 3-D window's script (needs Node.js)
```

Each test prints `PASS` or `FAIL` for every check, and exits with a non-zero code on a failure.

- **`test_codegen`**: for each example, the script that the app writes for the example's own
  values must run the same simulation as the tutorial's script: every snapshot of the network,
  and the final one, equal to the last bit, and the same printed lines but for the wall time.
  The tutorial's script runs from a copy whose only change is that its PyDiS paths are
  absolute (they are relative to its own folder). Also checked: the numbers in the scripts,
  each step's code inside the full script, the spec's JSON round trip, and selections that
  must be refused (a zero time step, minseg over half of maxseg, a NaN stress, ...).
- **`test_job`**: a run that ends well, one that fails, one that writes no final network, a
  stop from the caller (as Streamlit stops a script run), a time limit, a job left by an
  earlier session, the run folders, and their deletion, which keeps a folder whose job runs;
  each with a small stand-in script in a temporary folder, and none left running at the end.
- **`test_views`**: the initial configurations, as PyDiS builds them; a network written by
  `DisNetManager.write_json` and read back without PyDiS, across a periodic boundary; the
  periodic images, against PyDiS's own `Cell.closest_image`; the figures (their traces, the
  first view, the camera kept); and a run's results (300 snapshots played as 120, the history,
  the zip).
- **`test_app_smoke`**: the banner, the one drop-down menu, no chat, every step of both
  examples, a change and a reset, and a 20-step run of the binary junction: its output, status
  and downloads, and the networks that the 3-D window's component gets (AppTest does not run
  the component's script); then **🗑 clear all results**. Its runs go to a temporary folder
  (`OPENDIS_APP_RUNS`).
- **`test_network3d.js`**: the component's script in Node.js, with stand-ins for Plotly, the
  page and the browser: a run starts on its final network in the first view; the camera that
  Plotly reports after a drag is the one every later network of the example is drawn with, in
  its next view too, while another example, or a return to the first, starts in the first view;
  play, the slider, a hold while the mouse button is down, the end, and new networks.

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
- The Cloud has no compiled `lib/libpydis.so`, so PyDiS uses its Python version of the
  functions that the library holds. Its output starts with "pydis_lib not found".
