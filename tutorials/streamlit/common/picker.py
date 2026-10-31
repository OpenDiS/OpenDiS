"""common.picker: the navigation, in the sidebar: the example and its code, the steps, prev / next.

The picker is two drop-down menus, drawn inside the app's keyed container ``opendis_picker``:
the examples' labels, then the core libraries that run the chosen one, PyDiS or ExaDiS (the
keyed container ``opendis_core``). Each core of an example is an input file of its own
(``common/catalog.py``). The steps and the prev / next buttons come below. An example with
sub-examples, such as the binary junction's crystal structures, has a third menu, in step 1.
"""
from collections.abc import Collection, Sequence
from dataclasses import dataclass

import streamlit as st

from common.catalog import Example
from common.spec import ENGINE_NAMES, ENGINES, STEPS

PICKER_KEY = "picker_example"
CORE_KEY = "picker_core"  # one per label: "picker_core/<label>"
# The chosen sub-example, one per label: "picker_variant/<label>". Not a widget's key, so that
# it stays while step 1, which has the menu ("picker_variant_menu/<label>"), is not shown.
VARIANT_KEY = "picker_variant"


@dataclass(frozen=True)
class Step:
    """One step of the workflow, as the sidebar lists it."""

    key: str
    label: str


STEP_LIST = tuple(Step(key, f"{number}. {title}")
                  for number, (key, title) in enumerate(STEPS, start=1))  # a simulation's
BUILD_STEPS = (STEP_LIST[0], Step("run", "2. Run"))  # a script that only builds a network


def step_list(example: Example) -> tuple[Step, ...]:
    """The example's steps: seven for a simulation, two for a script that builds a network."""
    return STEP_LIST if example.default.config.simulates else BUILD_STEPS


def run_step(example: Example) -> Step:
    """The example's last step, Run."""
    return step_list(example)[-1]


def _ensure_valid(key: str, options: list[str]) -> None:
    """Reset a keyed widget's stored value when it is not among its options.

    This is called before the widget is created, because Streamlit raises when a
    keyed widget's value is not one of its options, e.g. after an example's file is renamed.
    """
    if options and st.session_state.get(key) not in options:
        st.session_state[key] = options[0]


def _variants(examples: Sequence[Example], label: str) -> list[str]:
    """The sub-examples of a label, in menu order; [""] when it has none."""
    return list(dict.fromkeys(e.variant for e in examples if e.label == label))


def example_picker(examples: Sequence[Example]) -> Example:
    """A label, a select box of the examples' labels, then one of the chosen example's cores;
    returns the chosen example, of the label's chosen sub-example."""
    st.markdown("**Example**")  # bold like "Steps"
    labels = list(dict.fromkeys(example.label for example in examples))
    _ensure_valid(PICKER_KEY, labels)
    label = st.selectbox("Example", labels, key=PICKER_KEY, label_visibility="collapsed")
    variants = _variants(examples, label)
    _ensure_valid(f"{VARIANT_KEY}/{label}", variants)
    variant = st.session_state.get(f"{VARIANT_KEY}/{label}", variants[0])
    cores = sorted((e for e in examples if e.label == label and e.variant == variant),
                   key=lambda e: ENGINES.index(e.engine))
    names = [ENGINE_NAMES[example.engine] for example in cores]
    key = f"{CORE_KEY}/{label}"
    _ensure_valid(key, names)
    with st.container(key="opendis_core"):  # its label is smaller (common/style.py)
        name = st.selectbox("Core", names, key=key, disabled=len(names) == 1,
                            help="The core library that runs the example's script: PyDiS, in "
                                 "Python, or ExaDiS, in C++ through its Python module pyexadis.")
    return cores[names.index(name)]


def _choose_variant(label: str, menu: str, steps_keys: list[str]) -> None:
    """The menu's choice becomes the label's sub-example, at its step 1 (an on_change callback)."""
    st.session_state[f"{VARIANT_KEY}/{label}"] = st.session_state[menu]
    for key in steps_keys:
        st.session_state[key] = STEP_LIST[0].label


def variant_picker(example: Example, examples: Sequence[Example]) -> None:
    """Step 1's menu of the example's sub-examples, such as its crystal structures; nothing if
    it has none."""
    variants = _variants(examples, example.label)
    if len(variants) < 2:
        return
    menu = f"{VARIANT_KEY}_menu/{example.label}"
    steps_keys = [_steps_key(e) for e in examples if e.label == example.label]
    st.selectbox(example.variant_label or "Variant", variants,
                 index=variants.index(example.variant), key=menu, on_change=_choose_variant,
                 args=(example.label, menu, steps_keys),
                 help="Each is a sub-example, with its own values in every step.")


def _steps_key(example: Example) -> str:
    """The session-state key of an example's steps radio (one per example)."""
    return f"steps_{example.name}"


def current_step(example: Example) -> Step:
    """The chosen step, read before the steps radio is drawn (it is drawn later in the run)."""
    steps = step_list(example)
    labels = [step.label for step in steps]
    key = _steps_key(example)
    _ensure_valid(key, labels)
    return steps[labels.index(st.session_state.get(key, labels[0]))]


def steps_radio(example: Example, marked: Collection[str] = ()) -> Step:
    """The steps as a radio, with " •" after the changed ones; returns the chosen step.

    ``marked`` holds the keys of the steps whose selections differ from the example's.
    """
    st.markdown("**Steps**")
    steps = step_list(example)
    labels = [step.label for step in steps]
    marks = {step.label for step in steps if step.key in marked}
    key = _steps_key(example)
    _ensure_valid(key, labels)
    label = st.radio("Steps", labels, key=key, label_visibility="collapsed",
                     format_func=lambda text: f"{text} •" if text in marks else text)
    return steps[labels.index(label)]


def _step_by(key: str, labels: list[str], delta: int) -> None:
    """Move the steps radio by ``delta``. This is an on_click callback, run before the rerun."""
    current = st.session_state.get(key, labels[0])
    i = labels.index(current) if current in labels else 0
    st.session_state[key] = labels[min(max(i + delta, 0), len(labels) - 1)]


def go_to(example: Example, step: Step) -> None:
    """Choose a step (an on_click callback)."""
    st.session_state[_steps_key(example)] = step.label


def prev_next(example: Example) -> None:
    """The ‹ prev and next › buttons that walk through the steps."""
    labels = [step.label for step in step_list(example)]
    key = _steps_key(example)
    current = st.session_state.get(key, labels[0])
    previous, forward = st.columns(2)
    previous.button("‹ prev", key=f"{key}_prev", width="stretch", on_click=_step_by,
                    args=(key, labels, -1), disabled=current == labels[0])
    forward.button("next ›", key=f"{key}_next", width="stretch", on_click=_step_by,
                   args=(key, labels, 1), disabled=current == labels[-1])
