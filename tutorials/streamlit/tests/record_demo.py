"""tests/record_demo.py: record the app's demo video, recordings/app_demo.mp4, with Playwright and Chrome."""
import argparse
import base64
import json
import os
import shutil
import signal
import subprocess
import tempfile
import time
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image
from playwright.sync_api import Frame, Page, sync_playwright  # pylint: disable=import-error

APP = Path(__file__).resolve().parents[1]
PORT = 8597
VIEWPORT = {"width": 1600, "height": 1000}
COMPONENT = 'iframe[title*="network3d"]'
STEPS = ("2. Material and lengths", "3. Forces and mobility", "4. Time integration",
         "5. Topology, collision and remesh", "6. Loading", "7. Run")
CODE_STEP = "5. Topology, collision and remesh"
NOTE = "Press Run in step 7 to simulate these selections."  # in the 3-D window before a run
DEG_PER_PX = 0.32  # about how far a drag turns the view
INK = 0.01  # a frame shows the 3-D view when this share of the view's pixels is not white
MARGIN_S = 0.3  # the screencast's frame times and the page's events agree to within this
# A visible cursor (headless Chrome records none), in the top page; the 3-D view's frame moves it too.
# It listens on the window: a frame keeps its window when its first, blank document is replaced,
# so listeners on that document would hear nothing.
CURSOR_JS = """
(() => {
  const dot = () => { try { return window.top.document.getElementById('demo-cursor'); }
                      catch (e) { return null; } };
  if (window === window.top) {
    window.addEventListener('DOMContentLoaded', () => {
      const d = document.createElement('div');
      d.id = 'demo-cursor';
      Object.assign(d.style, {position: 'fixed', zIndex: 2147483647, width: '18px',
        height: '18px', borderRadius: '50%', background: 'rgba(200, 30, 30, 0.55)',
        border: '2px solid white', boxShadow: '0 0 5px rgba(0, 0, 0, 0.5)',
        pointerEvents: 'none', transform: 'translate(-50%, -50%)', left: '-50px', top: '-50px'});
      document.body.appendChild(d);
    });
  }
  const offset = () => {  // this frame's corner in the top page
    let x = 0, y = 0, w = window;
    try { while (w !== window.top) { const r = w.frameElement.getBoundingClientRect();
                                     x += r.left; y += r.top; w = w.parent; } } catch (e) { }
    return [x, y];
  };
  window.addEventListener('pointermove', e => {
    const d = dot(), [x, y] = offset();
    if (d) { d.style.left = e.clientX + x + 'px'; d.style.top = e.clientY + y + 'px'; }
  }, true);
  const press = down => () => {
    const d = dot();
    if (d) { d.style.width = d.style.height = down ? '28px' : '18px';
             d.style.background = down ? 'rgba(200, 30, 30, 0.9)' : 'rgba(200, 30, 30, 0.55)'; }
  };
  window.addEventListener('pointerdown', press(true), true);
  window.addEventListener('pointerup', press(false), true);
})();
"""
# Scroll the app's main pane smoothly, to its 'bottom' or its 'top'
SCROLL_JS = """([to, duration]) => new Promise(resolve => {
  const main = document.querySelector('section[data-testid="stMain"]');
  const start = main.scrollTop, end = to === 'bottom' ? main.scrollHeight - main.clientHeight : 0;
  const t0 = performance.now();
  const frame = now => {
    const f = Math.min((now - t0) / duration, 1);
    main.scrollTop = start + (end - start) * (f < 0.5 ? 2 * f * f : 1 - (2 - 2 * f) ** 2 / 2);
    if (f < 1) requestAnimationFrame(frame); else resolve();
  };
  requestAnimationFrame(frame);
})"""
# The 3-D view's camera azimuth, in degrees, as Plotly's scene has it now
AZIMUTH_JS = """() => {
  const eye = document.getElementById('plot')._fullLayout.scene._scene.getCamera().eye;
  return Math.atan2(eye.y, eye.x) * 180 / Math.PI;
}"""


def glide_click(page: Page, target) -> None:
    """Move the visible cursor to the target's centre, pause, and click there."""
    target.scroll_into_view_if_needed()
    box = target.bounding_box()
    page.mouse.move(box["x"] + box["width"] / 2, box["y"] + box["height"] / 2, steps=28)
    page.wait_for_timeout(450)
    page.mouse.down()
    page.wait_for_timeout(90)
    page.mouse.up()


def view(page: Page) -> Frame:
    """The 3-D window's frame, once its figure is drawn."""
    frame = page.wait_for_selector(COMPONENT, timeout=120000).content_frame()
    frame.wait_for_selector("#plot .main-svg", timeout=120000)
    return frame


def turn(page: Page, frame: Frame, azimuth: float, step_px: int = 6) -> None:
    """Turn the 3-D view about its vertical axis to ``azimuth`` (degrees) by dragging across it; a
    drag to the right turns the cube's near side to the right, which lowers the azimuth."""
    def remaining() -> float:
        return (azimuth - frame.evaluate(AZIMUTH_JS) + 180) % 360 - 180

    sign = 1 if remaining() < 0 else -1  # +1: drag to the right
    length = abs(remaining()) / DEG_PER_PX
    box = page.locator(COMPONENT).first.bounding_box()
    x, y = box["x"] + box["width"] / 2 - sign * length / 2, box["y"] + box["height"] * 0.55
    page.mouse.move(x, y, steps=20)
    page.wait_for_timeout(300)
    page.mouse.down()
    for _ in range(int(3 * length / step_px)):
        if sign * remaining() >= 0:
            break
        x += sign * step_px
        page.mouse.move(x, y)
        page.wait_for_timeout(40)
    page.mouse.up()


def walk(page: Page, show_code: bool) -> tuple[list[float], dict]:
    """One example: turn its initial configuration and back, every step by next ›, (the code of
    step 5), Run, play the replay, and turn the view at its end. Returns the times in which the 3-D
    view shows the run (from its first snapshot until the replay is drawn) and the view's place."""
    frame = view(page)
    region = page.locator(COMPONENT).first.bounding_box()
    page.wait_for_timeout(1500)
    start = frame.evaluate(AZIMUTH_JS)
    turn(page, frame, start - 55)
    page.wait_for_timeout(1200)
    turn(page, frame, start)
    page.wait_for_timeout(1500)
    for step in STEPS:
        glide_click(page, page.get_by_role("button", name="next ›"))
        page.locator(".opendis-panel-title").filter(has_text=step).first.wait_for(timeout=60000)
        page.wait_for_timeout(2600)
        if show_code and step == CODE_STEP:
            page.evaluate(SCROLL_JS, ["bottom", 1600])
            page.wait_for_timeout(4500)
            page.evaluate(SCROLL_JS, ["top", 1600])
            page.wait_for_timeout(1000)
    glide_click(page, page.get_by_role("button", name="▶ Run"))
    page.get_by_text(NOTE, exact=True).wait_for(state="detached", timeout=120000)
    run = [time.time() - MARGIN_S]  # the screencast's frames have times on this clock
    page.get_by_text("Download the results").wait_for(timeout=600000)
    frame = view(page)
    run.append(time.time() + MARGIN_S)
    page.wait_for_timeout(2000)
    glide_click(page, frame.locator("#play"))
    frame.wait_for_function("() => document.getElementById('play').textContent.includes('play')",
                            timeout=300000)  # the end of the replay
    page.wait_for_timeout(1500)
    turn(page, frame, frame.evaluate(AZIMUTH_JS) + 40)
    page.wait_for_timeout(2500)
    return run, region


class Screencast:  # pylint: disable=too-few-public-methods
    """The tab's frames from Chrome's screencast, each with its time."""

    def __init__(self, page: Page, folder: Path) -> None:
        self.folder, self.frames = folder, []
        folder.mkdir(parents=True, exist_ok=True)
        self.cdp = page.context.new_cdp_session(page)
        self.cdp.on("Page.screencastFrame", self._frame)
        self.cdp.send("Page.startScreencast", {"format": "jpeg", "quality": 85,
                                               "maxWidth": VIEWPORT["width"],
                                               "maxHeight": VIEWPORT["height"]})

    def _frame(self, event: dict) -> None:
        """Keep a frame, and acknowledge it so that Chrome sends the next."""
        path = self.folder / f"{len(self.frames):06d}.jpg"
        path.write_bytes(base64.b64decode(event["data"]))
        self.frames.append((float(event["metadata"]["timestamp"]), path))
        self.cdp.send("Page.screencastFrameAck", {"sessionId": event["sessionId"]})

    def stop(self) -> list[tuple[float, Path]]:
        """Stop the screencast; returns the frames, oldest first."""
        self.cdp.send("Page.stopScreencast")
        return self.frames


def ink(path: Path, region: dict) -> float:
    """The share of the 3-D view's pixels, below its title and legend, that are not white."""
    x, y, w, h = (region[k] for k in ("x", "y", "width", "height"))
    with Image.open(path) as image:
        pixels = np.asarray(image.convert("RGB").crop((x + 20, y + 80, x + w - 20, y + h - 20)))
    return float((pixels.min(axis=2) < 225).mean())


def steady(frames: list[tuple[float, Path]], runs: list[list[float]], region: dict
           ) -> list[tuple[float, Path]]:
    """The frames without those of a run in which the 3-D view is blank (Plotly redraws it at each
    snapshot, so it blinks): a frame without the view is dropped, and the one before stays on
    screen."""
    dropped = {path for time_s, path in frames
               if any(start <= time_s <= end for start, end in runs) and ink(path, region) < INK}
    print(f"dropped {len(dropped)} of {len(frames)} frames, where a run's 3-D view was blank")
    return [(time_s, path) for time_s, path in frames if path not in dropped]


def to_mp4(frames: list[tuple[float, Path]], mp4: Path, ffmpeg: str, hold_s: float = 1.0) -> None:
    """Encode the frames, each shown until the next one's time, as an H.264 MP4 at 30 fps."""
    listing = frames[0][1].parent / "frames.txt"
    lines = []
    for (time_s, path), (next_s, _) in zip(frames, frames[1:] + [(frames[-1][0] + hold_s, None)]):
        lines += [f"file '{path.name}'", f"duration {max(next_s - time_s, 0.001):.4f}"]
    listing.write_text("\n".join(lines + [f"file '{frames[-1][1].name}'"]) + "\n",
                       encoding="utf-8")
    mp4.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-f", "concat", "-safe", "0",
                    "-i", str(listing), "-vf", "fps=30,scale=out_range=tv,format=yuv420p",
                    "-color_range", "tv", "-c:v", "libx264", "-crf", "22", "-preset", "medium",
                    "-movflags", "+faststart", str(mp4)],
                   check=True)


def record(work: Path) -> None:
    """Start the app on a spare port with an empty runs folder, and record the demo into ``work``:
    the frames, and ``frames.json`` with their times, the runs' times and the 3-D view's place."""
    env = {**os.environ, "STREAMLIT_SERVER_PORT": str(PORT), "OPENDIS_APP_RUNS": str(work / "runs")}
    with open(work / "app.log", "w", encoding="utf-8") as log:
        app = subprocess.Popen([str(APP / "launch_app.sh")], cwd=APP, env=env,  # pylint: disable=consider-using-with
                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    try:
        for _ in range(240):
            try:
                with urllib.request.urlopen(f"http://localhost:{PORT}/_stcore/health",
                                            timeout=2) as reply:
                    if reply.read().strip() == b"ok":
                        break
            except OSError:
                time.sleep(0.5)
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(channel="chrome", headless=True)
            context = browser.new_context(viewport=VIEWPORT)
            context.add_init_script(CURSOR_JS)
            page = context.new_page()
            page.goto(f"http://localhost:{PORT}")
            view(page)
            page.wait_for_timeout(1500)  # the first frame shows the app ready
            cast, runs, region = Screencast(page, work / "frames"), [], None
            try:
                page.mouse.move(800, 520, steps=10)
                page.wait_for_timeout(2500)
                run, region = walk(page, show_code=True)  # the default example, Frank-Read
                runs.append(run)
                glide_click(page, page.locator(
                    '.st-key-opendis_picker [data-testid="stSelectbox"]').first)
                glide_click(page, page.get_by_role("option", name="Binary junction"))
                page.locator(".opendis-panel-title").filter(
                    has_text="Binary junction").first.wait_for(timeout=60000)
                page.wait_for_timeout(2500)
                runs.append(walk(page, show_code=False)[0])
                page.wait_for_timeout(1500)
            finally:
                frames = cast.stop()
                browser.close()
    finally:
        os.killpg(app.pid, signal.SIGTERM)
        app.wait(timeout=10)
    (work / "frames.json").write_text(json.dumps({
        "frames": [(time_s, path.name) for time_s, path in frames], "runs": runs,
        "region": region}), encoding="utf-8")


def main() -> None:
    """Record the demo and write the MP4; ``--work`` keeps the frames, ``--encode`` encodes the
    frames of an earlier ``--work`` again, without recording."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, help="keep the frames in this folder")
    parser.add_argument("--encode", action="store_true", help="only encode the --work frames")
    args = parser.parse_args()
    try:
        import imageio_ffmpeg  # pylint: disable=import-outside-toplevel
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        ffmpeg = shutil.which("ffmpeg")
    out = APP / "recordings" / "app_demo.mp4"
    with tempfile.TemporaryDirectory(prefix="opendis_demo_") as tmp:
        work = (args.work or Path(tmp)).resolve()
        work.mkdir(parents=True, exist_ok=True)
        if not args.encode:
            record(work)
        data = json.loads((work / "frames.json").read_text(encoding="utf-8"))
        frames = [(time_s, work / "frames" / name) for time_s, name in data["frames"]]
        frames = steady(frames, data["runs"], data["region"])
        to_mp4(frames, out, ffmpeg)
    print(f"wrote {out}: {len(frames)} frames, {frames[-1][0] - frames[0][0]:.0f} s")


if __name__ == "__main__":
    main()
