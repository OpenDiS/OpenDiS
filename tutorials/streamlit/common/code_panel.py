"""common.code_panel: the code panel's content, this step's code or the full script.

``app.py`` draws the panel's title and its mode radio, and this module fills the panel.
"""
import streamlit as st

from common import codegen
from common.catalog import Example
from common.picker import Step
from common.spec import Spec

MODES = ("this step", "full script")
HEIGHT = 460  # px; the code scrolls inside


def show(mode: str, example: Example, spec: Spec, step: Step) -> None:
    """Fill the panel for its mode."""
    errors = spec.validate()
    if errors:
        st.error("Fix these selections first:\n\n" + "\n".join(f"- {e}" for e in errors))
        return
    if mode == "this step":
        st.code(codegen.step_code(spec, step.key, example.final_json), language="python",
                line_numbers=True, height=HEIGHT)
        return
    script = codegen.full_script(spec, example.final_json)
    st.download_button(f"⬇ download {example.script_name}", script,
                       file_name=example.script_name, mime="text/x-python", on_click="ignore",
                       key=f"download/{example.name}")
    st.code(script, language="python", line_numbers=True, height=HEIGHT)
