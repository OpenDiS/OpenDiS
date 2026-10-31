// tests/test_network3d.js: the 3-D window's script, with stand-ins for Plotly and the browser (README.md, Tests).
"use strict";
const fs = require("fs");
const path = require("path");
const vm = require("vm");
const source = fs.readFileSync(path.join(__dirname, "..", "components", "network3d", "network3d.js"), "utf8");
let failures = 0;
function check(name, ok, detail) {
  console.log(`${name}: ${ok ? "PASS" : "FAIL — " + JSON.stringify(detail)}`);
  if (!ok) failures += 1;
}

const storage = new Map();  // sessionStorage, shared by the instances of one "tab"

function instance() {
  const handlers = {};       // window's
  const elements = {};
  const posted = [];
  const draws = [];          // Plotly.react calls: the camera and the title
  const relayouts = [];
  const timers = {next: 1, live: new Map()};
  const files = new Map();   // what the page can fetch, by URL
  const fetched = [];
  function element(id) {
    const listeners = {};
    return elements[id] = {
      id, textContent: "", value: "", max: "", className: "", scrollHeight: 600, handlers: {},
      addEventListener(type, fn) { (listeners[type] = listeners[type] || []).push(fn); },
      fire(type, event) { (listeners[type] || []).forEach((fn) => fn(event || {})); },
      on(type, fn) { (this.handlers[type] = this.handlers[type] || []).push(fn); },
    };
  }
  ["plot", "controls", "play", "slider", "label"].forEach(element);
  const plot = elements.plot;
  const Plotly = {
    react(gd, data, layout) {
      draws.push({camera: layout.scene.camera, title: layout.title.text, n: data.length});
      return Promise.resolve(gd);
    },
    relayout(gd, update) {
      relayouts.push(update);
      (gd.handlers.plotly_relayout || []).forEach((fn) => fn(update));
      return Promise.resolve(gd);
    },
  };
  const window = {
    parent: {postMessage(message) { posted.push(message); }},
    addEventListener(type, fn) { (handlers[type] = handlers[type] || []).push(fn); },
    sessionStorage: {
      getItem(k) { return storage.has(k) ? storage.get(k) : null; },
      setItem(k, v) { storage.set(k, String(v)); },
      removeItem(k) { storage.delete(k); },
    },
  };
  const context = {
    window, Plotly, JSON, Number, String,
    document: {getElementById: (id) => elements[id], body: {scrollHeight: 600}},
    setInterval(fn) { const id = timers.next++; timers.live.set(id, fn); return id; },
    clearInterval(id) { timers.live.delete(id); },
    setTimeout(fn) { setImmediate(fn); },
    Date,
    fetch(url) {
      fetched.push(url);
      const body = files.get(url.split("?")[0]);
      return Promise.resolve(body === undefined ? {ok: false}
                             : {ok: true, json: () => Promise.resolve(JSON.parse(body))});
    },
  };
  vm.runInNewContext(source, context);
  return {
    posted, draws, relayouts, elements, plot, timers, files, fetched,
    render(args) { (handlers.message || []).forEach((fn) => fn({data: {type: "streamlit:render", args}})); },
    tick() { [...timers.live.values()].forEach((fn) => fn()); },
    pointer(type) { (type === "down" ? plot : {fire: (t) => (handlers[t] || []).forEach((fn) => fn({}))}).fire(type === "down" ? "pointerdown" : "pointerup"); },
  };
}

const same = (x, y) => JSON.stringify(x) === JSON.stringify(y);
const settle = () => new Promise((resolve) => setImmediate(resolve));
const CAMERA = {eye: {x: 1.138, y: -1.484, z: 1.091}};
const USER = {eye: {x: 0.2, y: -2.0, z: 0.5}};
function frames(labels) { return labels.map((l) => ({label: l, title: l, data: [{}, {}, {}, {}]})); }
function args(labels, id, extra) {
  return Object.assign({frames: frames(labels), layout: {scene: {camera: CAMERA}, title: {}},
                        config: {}, start: "last", frame_ms: 300, data_id: id,
                        scope: "frank_read_src"}, extra || {});
}

(async () => {
  const a = instance();
  check("it tells Streamlit that it is ready", a.posted[0] && a.posted[0].type === "streamlit:componentReady"
        && a.posted[0].apiVersion === 1, a.posted[0]);
  a.render(args(["step 0", "step 10", "final"], "run1"));
  await settle(); await settle();
  check("a run starts on its final network, in the app's first view",
        a.draws.length === 1 && a.draws[0].title === "final" && same(a.draws[0].camera, CAMERA), a.draws);
  check("the controls show, for more than one snapshot", a.elements.controls.className === "on");
  check("its height is sent", a.posted.some((m) => m.type === "streamlit:setFrameHeight" && m.height === 600));
  a.plot.handlers.plotly_relayout.forEach((fn) => fn({"scene.camera": USER}));  // the user turns it
  a.elements.play.fire("click");
  check("play at the end goes back to step 0, with the user's camera",
        a.draws.at(-1).title === "step 0" && a.draws.at(-1).camera === USER, a.draws.at(-1));
  a.tick();
  check("a tick draws the next snapshot, with the user's camera",
        a.draws.at(-1).title === "step 10" && a.draws.at(-1).camera === USER, a.draws.at(-1));
  const before = a.draws.length;
  a.pointer("down");
  a.tick();
  check("no snapshot is drawn while the mouse button is down", a.draws.length === before);
  a.pointer("up");
  a.tick();
  check("after the drag it goes on", a.draws.at(-1).title === "final", a.draws.at(-1));
  a.tick();
  check("at the end it stops", a.timers.live.size === 0 && a.elements.play.textContent === "▶ play");
  a.render(args(["step 0", "step 10", "final"], "run1"));
  await settle();
  check("the same networks again (a rerun) are not drawn again", a.draws.length === before + 1, a.draws.length - before);
  a.render(args(["step 0", "final"], "run2", {start: "first"}));
  await settle();
  check("new networks are drawn with the user's camera, from the first if asked",
        a.draws.at(-1).title === "step 0" && a.draws.at(-1).camera === USER, a.draws.at(-1));
  a.elements.slider.value = "1";
  a.elements.slider.fire("input");
  check("the slider draws its snapshot", a.draws.at(-1).title === "final");

  const b = instance();  // the next view, in the same tab
  b.render(args(["initial configuration"], "init1"));
  await settle(); await settle();
  check("the next view starts in the app's first view, then takes the tab's camera",
        same(b.draws[0].camera, CAMERA) && b.relayouts.length === 1 &&
        JSON.stringify(b.relayouts[0]["scene.camera"]) === JSON.stringify(USER), b.relayouts);
  check("one network: no controls", b.elements.controls.className === "");
  b.render(args(["initial configuration"], "init2"));
  await settle();
  check("a changed value keeps the camera",
        JSON.stringify(b.draws.at(-1).camera) === JSON.stringify(USER), b.draws.at(-1));
  const c = instance();  // another example, in the same tab
  c.render(args(["initial configuration"], "bj1", {scope: "binary_junction"}));
  await settle(); await settle();
  check("another example starts in the app's first view, and drops the kept camera",
        same(c.draws[0].camera, CAMERA) && c.relayouts.length === 0 && !storage.has("opendis-app-camera"),
        [c.draws[0], c.relayouts]);
  const d = instance();  // back to the first example
  d.render(args(["initial configuration"], "fr1"));
  await settle(); await settle();
  check("back to the first example: the first view again", same(d.draws[0].camera, CAMERA)
        && d.relayouts.length === 0, d.relayouts);

  // during a run: the job's newest snapshot, which the app writes into live/
  const settled = async () => { for (let i = 0; i < 6; i++) await settle(); };
  const snapshot = (label, id) => JSON.stringify({frames: frames([label]), data_id: id,
                                                  layout: {scene: {camera: CAMERA}, title: {}}});
  const DRAG = {eye: {x: -1.5, y: 0.3, z: 1.2}};
  const e = instance();
  e.render({live: "live/r1.json", live_ms: 500, config: {}, scope: "binary_junction"});
  await settled();
  check("during a run, nothing is drawn before the first snapshot, and the page fetches again",
        e.draws.length === 0 && e.timers.live.size === 1 && e.fetched[0].startsWith("live/r1.json?t="),
        [e.draws, e.fetched]);
  e.files.set("live/r1.json", snapshot("step 0", "s0"));
  e.tick();
  await settled();
  check("a run's first snapshot is drawn in the app's first view, without controls",
        e.draws.length === 1 && e.draws[0].title === "step 0" && same(e.draws[0].camera, CAMERA)
        && e.elements.controls.className === "", e.draws);
  e.plot.handlers.plotly_relayout.forEach((fn) => fn({"scene.camera": USER}));
  e.files.set("live/r1.json", snapshot("step 10", "s10"));
  e.tick();
  await settled();
  check("the run's next snapshot keeps the user's camera",
        e.draws.at(-1).title === "step 10" && e.draws.at(-1).camera === USER, e.draws.at(-1));
  e.tick();
  await settled();
  check("the same snapshot is not drawn again", e.draws.length === 2, e.draws.length);
  e.pointer("down");
  e.plot.handlers.plotly_relayouting.forEach((fn) => fn({"scene.camera": DRAG}));
  e.files.set("live/r1.json", snapshot("step 20", "s20"));
  e.tick();
  await settled();
  check("a new snapshot waits while the mouse button is down", e.draws.length === 2, e.draws.length);
  e.pointer("up");
  await settled();
  check("then it is drawn, with the camera that the drag got to",
        e.draws.at(-1).title === "step 20" && e.draws.at(-1).camera === DRAG, e.draws.at(-1));
  e.render(args(["step 0", "step 10", "final"], "run3", {scope: "binary_junction"}));
  await settled();
  const count = e.fetched.length;
  e.tick();
  check("the replay stops the fetching, and keeps the camera",
        e.timers.live.size === 0 && e.fetched.length === count && e.draws.at(-1).title === "final"
        && e.draws.at(-1).camera === DRAG, [e.timers.live.size, e.draws.at(-1)]);
  const f = instance();  // the replay's own frame, after the run
  f.render(args(["step 0", "final"], "run3", {scope: "binary_junction"}));
  await settled();
  check("the replay's frame takes the camera that the run's view had",
        f.relayouts.length === 1 && same(f.relayouts[0]["scene.camera"], DRAG), f.relayouts);
  console.log(failures ? `\n${failures} FAILED` : "\nALL PASSED");
  process.exit(failures ? 1 : 0);
})();
