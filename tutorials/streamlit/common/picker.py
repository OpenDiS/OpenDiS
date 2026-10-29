"""common.picker: the navigation: the example picker (in the banner), the steps and prev / next.

The picker, one drop-down menu of the examples' labels, is drawn inside the app's keyed container
``opendis_picker``, which the CSS in ``common/style.py`` pins into the banner. The steps and the
prev / next buttons are in the sidebar.
"""
from collections.abc import Collection, Sequence
from dataclasses import dataclass

import streamlit as st

from common.catalog import Example
from common.spec import STEPS

PICKER_KEY = "picker_example"


@dataclass(frozen=True)
class Step:
    """One step of the workflow, as the sidebar lists it."""

    key: str
    label: str


STEP_LIST = tuple(Step(key, f"{number}. {title}")
                  for number, (key, title) in enumerate(STEPS, start=1))
RUN_STEP = next(step for step in STEP_LIST if step.key == "run")


def _ensure_valid(key: str, options: list[str]) -> None:
    """Reset a keyed widget's stored value when it is not among its options.

    This is called before the widget is created, because Streamlit raises when a
    keyed widget's value is not one of its options, e.g. after an example's file is renamed.
    """
    if options and st.session_state.get(key) not in options:
        st.session_state[key] = options[0]


def example_picker(examples: Sequence[Example]) -> Example:
    """A label and one select box of the examples; returns the chosen example."""
    st.markdown("**Example**")  # bold like the sidebar's "Steps"
    labels = [example.label for example in examples]
    _ensure_valid(PICKER_KEY, labels)
    label = st.selectbox("Example", labels, key=PICKER_KEY, label_visibility="collapsed")
    return examples[labels.index(label)]


def _steps_key(example: Example) -> str:
    """The session-state key of an example's steps radio (one per example)."""
    return f"steps_{example.name}"


def current_step(example: Example) -> Step:
    """The chosen step, read before the steps radio is drawn (it is drawn later in the run)."""
    labels = [step.label for step in STEP_LIST]
    key = _steps_key(example)
    _ensure_valid(key, labels)
    return STEP_LIST[labels.index(st.session_state.get(key, labels[0]))]


def steps_radio(example: Example, marked: Collection[str] = ()) -> Step:
    """The steps as a radio, with " •" after the changed ones; returns the chosen step.

    ``marked`` holds the keys of the steps whose selections differ from the example's.
    """
    st.markdown("**Steps**")
    labels = [step.label for step in STEP_LIST]
    marks = {step.label for step in STEP_LIST if step.key in marked}
    key = _steps_key(example)
    _ensure_valid(key, labels)
    label = st.radio("Steps", labels, key=key, label_visibility="collapsed",
                     format_func=lambda text: f"{text} •" if text in marks else text)
    return STEP_LIST[labels.index(label)]


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
    labels = [step.label for step in STEP_LIST]
    key = _steps_key(example)
    current = st.session_state.get(key, labels[0])
    previous, forward = st.columns(2)
    previous.button("‹ prev", key=f"{key}_prev", width="stretch", on_click=_step_by,
                    args=(key, labels, -1), disabled=current == labels[0])
    forward.button("next ›", key=f"{key}_next", width="stretch", on_click=_step_by,
                   args=(key, labels, 1), disabled=current == labels[-1])
