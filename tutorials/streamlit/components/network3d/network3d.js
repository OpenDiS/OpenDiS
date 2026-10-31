// The 3-D window of the OpenDiS app: a Streamlit component that draws the networks it is given,
// one network or the snapshots of a run, as a Plotly figure. Unlike st.plotly_chart, it keeps the
// camera that the user turned to: Plotly reports the camera after each drag, and every network is
// drawn with it, in this view and the next of the same example (the browser tab keeps it, with
// the example's name, its "scope"); another example starts in the app's first view. A run's
// replay plays here, in the browser; a snapshot is not drawn while the mouse button is down, so a
// drag goes on undisturbed. During a run, the job's newest snapshot, which the app writes next to
// this page, is fetched and drawn the same way. Plotly's "Reset camera to default" goes back to
// the app's first view.
(function () {
  "use strict";
  var STORE = "opendis-app-camera";  // {scope, camera}, kept for the browser tab (sessionStorage)
  var plot = document.getElementById("plot");
  var controls = document.getElementById("controls");
  var playButton = document.getElementById("play");
  var slider = document.getElementById("slider");
  var label = document.getElementById("label");
  var state = {
    dataId: null, scope: null, frames: [], layout: null, config: {}, index: 0, frameMs: 300,
    camera: null, timer: null, dragging: false, started: false,
    live: null, liveArgs: null, liveTimer: null, held: null
  };

  function send(type, data) {
    var message = {isStreamlitMessage: true, type: type};
    for (var name in data) { message[name] = data[name]; }
    window.parent.postMessage(message, "*");
  }

  function setHeight() {
    send("streamlit:setFrameHeight", {height: document.body.scrollHeight});
  }

  // The camera that the tab kept for this example; another example's is dropped
  function loadCamera() {
    try {
      var kept = JSON.parse(window.sessionStorage.getItem(STORE));
      if (kept && kept.scope === state.scope) { return kept.camera; }
      window.sessionStorage.removeItem(STORE);
    } catch (e) { /* no storage: the first view */ }
    return null;
  }

  function saveCamera() {
    try {
      window.sessionStorage.setItem(STORE, JSON.stringify({scope: state.scope,
                                                           camera: state.camera}));
    } catch (e) { /* no storage */ }
  }

  // Draw the current snapshot, with the user's camera once there is one
  function draw() {
    var frame = state.frames[state.index];
    var layout = JSON.parse(JSON.stringify(state.layout));
    if (state.camera) { layout.scene.camera = state.camera; }
    layout.title = layout.title || {};
    layout.title.text = frame.title;
    slider.value = String(state.index);
    label.textContent = frame.label;
    return Plotly.react(plot, frame.data, layout, state.config);
  }

  // The first draw is in the app's first view, so that Plotly's reset goes back to it
  function firstDraw() {
    state.started = true;
    return draw().then(function () {
      var turned = function (update) {  // during a drag too, which may end outside the plot
        if (update && update["scene.camera"]) {
          state.camera = update["scene.camera"];
          saveCamera();
        }
      };
      plot.on("plotly_relayout", turned);
      plot.on("plotly_relayouting", turned);
    });
  }

  function stop() {
    if (state.timer !== null) { clearInterval(state.timer); state.timer = null; }
    playButton.textContent = "▶ play";
  }

  function tick() {
    if (state.dragging) { return; }  // hold while the user turns the view
    if (state.index + 1 >= state.frames.length) { stop(); return; }
    state.index += 1;
    draw();
  }

  function play() {
    if (state.index + 1 >= state.frames.length) {  // at the end: from the start
      state.index = 0;
      draw();
    }
    playButton.textContent = "❚❚ pause";
    state.timer = setInterval(tick, state.frameMs);
  }

  playButton.addEventListener("click", function () {
    if (state.timer !== null) { stop(); } else { play(); }
  });
  slider.addEventListener("input", function () {
    state.index = Number(slider.value);
    draw();
  });
  function released() {
    state.dragging = false;
    setTimeout(showHeld, 0);  // after Plotly has reported the drag's last camera
  }
  plot.addEventListener("pointerdown", function () { state.dragging = true; });
  window.addEventListener("pointerup", released);
  window.addEventListener("pointercancel", released);

  // During a run: the job's newest snapshot, from the file that the app writes; a new one is
  // drawn at once, or when the mouse button is up again
  function poll() {
    var url = state.live;
    fetch(url + "?t=" + Date.now(), {cache: "no-store"}).then(function (response) {
      return response.ok ? response.json() : null;
    }).then(function (snapshot) {
      if (!snapshot || url !== state.live || snapshot.data_id === state.dataId) { return; }
      state.held = Object.assign({}, state.liveArgs, snapshot);
      showHeld();
    }).catch(function () { /* not written yet, or the run has ended */ });
  }

  function showHeld() {
    if (state.held && !state.dragging) {
      var args = state.held;
      state.held = null;
      show(args);
    }
  }

  function startLive(args) {
    stopLive();
    state.live = args.live;
    state.liveArgs = {config: args.config, scope: args.scope, start: "last"};
    poll();
    state.liveTimer = setInterval(poll, args.live_ms || 500);
  }

  function stopLive() {
    if (state.liveTimer !== null) { clearInterval(state.liveTimer); state.liveTimer = null; }
    state.live = null;
    state.held = null;
  }

  // Streamlit sends the arguments on every rerun; only new networks are drawn again, and another
  // example's in the app's first view. During a run, they name the file to fetch (live)
  function onRender(args) {
    if (args.live) {
      if (args.live !== state.live) { startLive(args); }
      return;
    }
    stopLive();
    if (!args.frames || (args.data_id === state.dataId && args.scope === state.scope)) { return; }
    show(args);
  }

  function show(args) {
    stop();
    var kept = null;  // a camera to apply after the draw, which Plotly reports like a drag
    if (args.scope !== state.scope) {
      state.scope = args.scope;
      state.camera = null;   // the app's first view, ...
      kept = loadCamera();   // ... or the camera that the tab kept for this example
    }
    state.dataId = args.data_id;
    state.frames = args.frames;
    state.layout = args.layout;
    state.config = args.config || {};
    state.frameMs = args.frame_ms || 300;
    state.index = args.start === "first" ? 0 : state.frames.length - 1;
    slider.max = String(state.frames.length - 1);
    controls.className = state.frames.length > 1 ? "on" : "";
    (state.started ? draw() : firstDraw()).then(function () {
      return kept ? Plotly.relayout(plot, {"scene.camera": kept}) : null;
    }).then(setHeight);
  }

  window.addEventListener("message", function (event) {
    if (event.data && event.data.type === "streamlit:render") {
      onRender(event.data.args || {});
    }
  });
  send("streamlit:componentReady", {apiVersion: 1});
})();
