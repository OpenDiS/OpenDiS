# engines/

The 50 engines of the OpenDiS tutorial, `opendis-00.streamlit.app` to `opendis-49.streamlit.app`. Each engine runs the app of [`app.py`](../app.py) on a server of its own, for one user. So the users of a tutorial do not slow each other down, and do not share results. The front page, [opendis.streamlit.app](https://opendis.streamlit.app), links to the engines. It shows how many users each one has at that moment ([`common/front_page.py`](../common/front_page.py)).

## Files

| File | Purpose |
|------|---------|
| `opendis_00.py` to `opendis_49.py` | The entry files of the 50 engines. Streamlit Community Cloud deploys one app per repository, branch and entry file. Each file calls `main()` of `engine.py`. |
| [`engine.py`](engine.py) | Runs the app as `app.py` does. It puts the app's folder first on `sys.path` and on `PYTHONPATH`. Streamlit reloads a changed module only in the entry file's folder or on `PYTHONPATH`. So a push updates an engine in place, as it updates `app.py`. It also resets the engine when `RESET` changes. It puts the engine's number in the banner, with its users, and signs off an idle session. |
| [`RESET`](RESET) | One line. A push that changes it resets every engine. At its next page run, each engine deletes the results of its past runs and clears its caches. Its builds of PyDiS and ExaDiS stay. |
| [`manage.py`](manage.py) | The engines from a terminal: `status`, `open`, `links` and `reset`. |
| [`announcements.md`](announcements.md) | The announcements of the front page, in Markdown. |

## Deploy the engines, once

No script drives the Cloud's dashboard: section 7 of its [terms of use](https://streamlit.io/terms-of-use) bars robots. So the engines are deployed by hand, in the [workspace](https://share.streamlit.io).

1. `python manage.py links` prints each engine's name and the GitHub link of its entry file.
2. For each engine: **Create app**, deploy a public app from GitHub, and paste the link of its entry file. The App URL is the engine's name, such as `opendis-07`. In **Advanced settings**, choose the Python of opendis.streamlit.app (3.14 on 2026-11-01).
3. Deploy `opendis-00` first, and check it with `python manage.py status`. Then deploy the other 49.

## The front page

opendis.streamlit.app shows the front page when its secret `OPENDIS_FRONT_PAGE` is true. Set it once, in the app's **Settings**, under **Secrets**: `OPENDIS_FRONT_PAGE = true`. Without it, the deployment shows the app.

The page reads each engine's gauge `active_sessions` at `https://opendis-NN.streamlit.app/~/+/_stcore/metrics`. Streamlit publishes it without a login. A session is one open browser tab. The table refreshes itself every 10 s, for 30 minutes after the page opened. The 50 engines are read at most once in 10 s, for all visitors together. An engine whose gauge does not come is shown as asleep: it may also be starting, or failing. An engine that the Cloud does not know is shown as not deployed.

Below the announcements, a short guide says how to use an engine. Everyone who opens the same engine shares it. Closing the engine's tab releases it. An idle engine signs its user off (below). After a sign-off, the front page says which engine signed the visitor off.

## Announcements

The front page shows the text of [`announcements.md`](announcements.md) in a framed box, below the banner and above the table. To post an announcement, write it there in Markdown, commit it, and push it to GitHub's `main`. An open front page shows the new text within 10 s. An empty file shows no box.

## In each engine

- **The banner** reads "OpenDiS 00" in `opendis-00`, and so on, as does the browser tab. At its right end, the engine shows its users at that moment, refreshed every 10 s: its open browser tabs, the count of the front page. With more than one, the count turns amber.
- **The idle sign-off.** A session without a click for 60 minutes is signed off: its tab goes back to the front page, which frees the engine. A warning comes 5 minutes before, with a **Stay signed in** button. A run in progress holds the check until it ends. A push counts as a click, since it reruns the page. Use that stays in the browser counts too: a press, a key or the wheel anywhere in the tab, the 3-D window and the charts included. Such use in the warning's 5 minutes closes the warning by itself. Pointer moves alone do not count: the browser may send them by itself.
- **Its settings** are `OPENDIS_IDLE_MINUTES` (60; 0 turns the sign-off off) and `OPENDIS_IDLE_WARNING_MINUTES` (5). Each comes from the environment, else from the engine's secrets, else the default. `OPENDIS_FRONT_URL` is where a signed-off tab goes, https://opendis.streamlit.app by default.
- **On the Cloud,** the app runs in a sandboxed frame, which may not move the tab. So the sign-off follows a link in the tab's own page, of the same origin. If that fails, the frame leaves the app, which ends the session all the same.

## Manage the engines

| Command | What it does |
|---------|--------------|
| `python manage.py status` | Prints the front page's table. |
| `python manage.py open` | Opens in the browser each engine that does not answer. Press the button of each asleep one. Opening an engine also starts its builds, about 6 minutes for ExaDiS. So do it about 10 minutes before a tutorial. `--all` opens every engine. |
| `python manage.py reset` | Writes a new line in `RESET`. Commit it, and push it to GitHub's `main`. Each engine then resets at its next page run. Do it when nobody uses the engines. |

An engine sleeps after 12 hours without visitors. A push to `main` updates every engine in place. A full reboot of an engine, from a fresh clone with fresh dependencies, is done by hand: the engine's menu in the workspace, **Reboot**. It also rebuilds ExaDiS.

## Test on this computer

`python -m tests.test_engines`, from `tutorials/streamlit/`, tests the entry files, `engine.py`, the front page and `manage.py`. To see the front page read live engines, start two engines and the front page from `tutorials/streamlit/`, each in a terminal:

```
OPENDIS_APP_RUNS=runs_00 python -m streamlit run engines/opendis_00.py --server.port 8700
OPENDIS_APP_RUNS=runs_01 python -m streamlit run engines/opendis_01.py --server.port 8701
OPENDIS_FRONT_PAGE=1 OPENDIS_ENGINE_PAGE="http://localhost:87{n:02d}" \
  OPENDIS_ENGINE_METRICS="http://localhost:87{n:02d}/_stcore/metrics" \
  python -m streamlit run app.py --server.port 8750
```

Engines that do not run are shown as asleep. For a quick look at the idle sign-off, start an engine with `OPENDIS_IDLE_MINUTES=1 OPENDIS_IDLE_WARNING_MINUTES=0.5 OPENDIS_FRONT_URL=http://localhost:8750`.
