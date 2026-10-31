"""common.style: the page CSS, the pinned banner and panel titles.

The rules come from the example app that this app follows, a Streamlit app for FEniCSx, which
tuned them for Streamlit 1.61: re-check the ``data-testid`` selectors after an upgrade. Its chat
column became the output column. This module does not import Streamlit; ``app.py`` passes these
strings to ``st.markdown``.
"""
import base64
import html
from pathlib import Path

BANNER_TITLE = "OpenDiS"

# A plain string, not an f-string, so that the CSS braces stay single. Two CSS
# lines (the font URL and one selector) are longer than 99 characters.
# pylint: disable=line-too-long
CSS = """<style>
@import url('https://fonts.googleapis.com/css2?family=STIX+Two+Text:ital,wght@0,400;0,600;0,700;1,400;1,600&display=swap');
:root {
  --opendis-banner-h: 2.9rem;
  --opendis-serif: "STIX Two Text", STIXGeneral, "Times New Roman", Times, serif;
}
/* Pinned top banner — light chrome, centered title */
.opendis-banner {
  position: fixed !important; top: 0 !important; left: 0 !important; right: 0 !important;
  width: 100vw !important; height: var(--opendis-banner-h) !important;
  box-sizing: border-box !important; z-index: 2147483000 !important;
  background: #e6eef6; color: #1a2332;
  padding: 0.35rem 1rem !important;
  border-bottom: 1px solid #b8c5d4;
  box-shadow: 0 1px 4px rgba(26, 35, 50, 0.08);
  text-align: center !important;
  display: flex !important; align-items: center !important; justify-content: center !important;
}
.opendis-banner h1 {
  font-family: var(--opendis-serif) !important; font-size: 1.35rem !important;
  font-weight: 650; margin: 0 !important; line-height: 1.15 !important;
  letter-spacing: 0.01em; color: #1a2332 !important; text-align: center !important;
}
/* Push sidebar under the banner (do not cover title) */
section[data-testid="stSidebar"] {
  z-index: 100 !important;
  top: var(--opendis-banner-h) !important;
  height: calc(100vh - var(--opendis-banner-h)) !important;
  max-height: calc(100vh - var(--opendis-banner-h)) !important;
  padding-top: 0.85rem !important;
}
/* Streamlit 1.61: expand control is stExpandSidebarButton in the header
   toolbar — the fixed banner was covering it. Pin it below the banner. */
header[data-testid="stHeader"] {
  background: transparent !important;
  top: var(--opendis-banner-h) !important;
  z-index: 2147483001 !important;
}
[data-testid="stToolbar"] { z-index: 2147483001 !important; }
/* The header's empty band below the banner lets clicks through; only its buttons take
   clicks */
header[data-testid="stHeader"],
header[data-testid="stHeader"] [data-testid="stToolbar"] { pointer-events: none !important; }
header[data-testid="stHeader"] button,
header[data-testid="stHeader"] a,
header[data-testid="stHeader"] [role="button"],
[data-testid="stExpandSidebarButton"] { pointer-events: auto !important; }
[data-testid="stExpandSidebarButton"] {
  position: fixed !important;
  top: calc(var(--opendis-banner-h) + 0.4rem) !important; left: 0.45rem !important;
  z-index: 2147483002 !important;
  display: flex !important; align-items: center !important; justify-content: center !important;
  gap: 0.25rem !important; width: auto !important; min-width: 2.6rem !important;
  height: 2.6rem !important; padding: 0.35rem 0.65rem !important;
  background: #f7fafc !important; border: 1px solid #9aa8b8 !important;
  border-radius: 0.5rem !important; box-shadow: 0 1px 3px rgba(26, 35, 50, 0.08) !important;
  color: #1a2332 !important;
}
[data-testid="stExpandSidebarButton"]::after {
  content: "Example · steps"; color: #1a2332 !important; font-size: 1.05rem !important;
  font-weight: 650 !important; margin-left: 0.15rem;
}
[data-testid="stExpandSidebarButton"] svg {
  width: 1.4rem !important; height: 1.4rem !important;
  stroke: #1a2332 !important; fill: none !important; color: #1a2332 !important;
}
[data-testid="stExpandSidebarButton"]:hover {
  background: #eef3f8 !important; border-color: #5b6b7c !important;
}
/* Legacy selectors (older Streamlit) */
[data-testid="collapsedControl"],
[data-testid="stSidebarCollapsedControl"] {
  top: calc(var(--opendis-banner-h) + 0.45rem) !important; left: 0.35rem !important;
  z-index: 2147483002 !important; background: #f7fafc !important;
  border: 1px solid #9aa8b8 !important; color: #1a2332 !important;
}
/* Main pane clears the banner */
div[data-testid="stAppViewContainer"] { padding-top: 0 !important; }
div[data-testid="stAppViewContainer"] > .main {
  padding-top: calc(var(--opendis-banner-h) + 0.65rem) !important;
}
.opendis-banner .opendis-logo { height: 2.3rem; width: auto; margin-right: 0.6rem; }
/* In-sidebar collapse — native Material « double-arrow (soft gray, not black) */
[data-testid="stSidebarCollapseButton"] {
  background: transparent !important;
  border: 1px solid #c5cad1 !important;
  border-radius: 0.35rem !important;
  padding: 0.2rem 0.35rem !important;
  z-index: 102 !important;
  display: inline-flex !important;
  align-items: center !important;
  justify-content: center !important;
  min-width: 2.0rem !important;
  min-height: 1.85rem !important;
  color: #8b949e !important;
  box-shadow: none !important;
}
[data-testid="stSidebarCollapseButton"]::before {
  content: none !important;
}
[data-testid="stSidebarCollapseButton"] svg {
  display: block !important;
  width: 1.35rem !important;
  height: 1.35rem !important;
  fill: #8b949e !important;
  color: #8b949e !important;
  stroke: none !important;
}
[data-testid="stSidebarCollapseButton"] span,
[data-testid="stSidebarCollapseButton"] i {
  display: inline-flex !important;
  color: #8b949e !important;
  font-size: 1.35rem !important;
  line-height: 1 !important;
}
[data-testid="stSidebarCollapseButton"]:hover {
  border-color: #9aa8b8 !important;
  background: #eef3f8 !important;
  color: #5b6b7c !important;
}
[data-testid="stSidebarCollapseButton"]:hover svg,
[data-testid="stSidebarCollapseButton"]:hover span,
[data-testid="stSidebarCollapseButton"]:hover i {
  fill: #5b6b7c !important;
  color: #5b6b7c !important;
}
/* Larger app chrome fonts */
html, body, [data-testid="stAppViewContainer"],
[data-testid="stMarkdownContainer"],
[data-testid="stWidgetLabel"] {
  font-size: 1.12rem !important;
}
div[data-testid="stRadio"] label,
div[data-testid="stRadio"] p,
div[data-testid="stMarkdownContainer"] p,
div[data-testid="stMarkdownContainer"] span,
div[data-testid="stExpander"] summary {
  font-size: 1.2rem !important;
}
div[data-testid="stMarkdownContainer"] h2,
div[data-testid="stMarkdownContainer"] h3,
[data-testid="stHeadingWithActionElements"] h2,
[data-testid="stHeadingWithActionElements"] h3,
h2, h3 {
  font-size: 1.55rem !important;
  margin-bottom: 0.45rem !important;
}
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] p,
section[data-testid="stSidebar"] div[data-testid="stMarkdownContainer"] strong,
section[data-testid="stSidebar"] label,
section[data-testid="stSidebar"] p {
  font-size: 1.22rem !important;
}
/* Tighten main-column padding so plots use more horizontal space */
div[data-testid="stMainBlockContainer"] {
  padding-left: 0.45rem !important;
  padding-right: 0.45rem !important;
  max-width: 100% !important;
}
div[data-testid="stHorizontalBlock"] {
  gap: 0.35rem !important;
}
div[data-testid="stPlotlyChart"] {
  width: 100% !important;
  display: flex !important;
  justify-content: center !important;
}
div[data-testid="stPlotlyChart"] > div {
  margin-left: auto !important;
  margin-right: auto !important;
}
/* Centred panel headings */
div[data-testid="stHeadingWithActionElements"] {
  justify-content: center !important;
  text-align: center !important;
  width: 100% !important;
}
div[data-testid="stHeadingWithActionElements"] h2,
div[data-testid="stHeadingWithActionElements"] h3,
[data-testid="stVerticalBlockBorderWrapper"] h2,
[data-testid="stVerticalBlockBorderWrapper"] h3 {
  text-align: center !important;
  width: 100% !important;
}
.opendis-panel-title {
  text-align: center !important;
  font-family: var(--opendis-serif) !important;
  font-size: 1.45rem !important;
  font-weight: 650 !important;
  margin: 0.05rem 0 0.25rem 0 !important;
}
/* Centred mode radios (the 3-D window's views, the code panel's modes) */
div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"] div[data-testid="stRadio"],
div[data-testid="stRadio"] {
  text-align: center !important;
}
div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"] div[data-testid="stRadio"] > div,
div[data-testid="stRadio"] > div[role="radiogroup"],
div[data-testid="stRadio"] div[role="radiogroup"] {
  display: flex !important;
  justify-content: center !important;
  align-items: center !important;
  flex-wrap: wrap !important;
  gap: 0.75rem !important;
  width: 100% !important;
}
div[data-testid="stRadio"] label {
  justify-content: center !important;
}
/* The 2nd left column (the keyed container opendis_selections): smaller labels, titles,
   buttons and text. Captions and inputs are already small and keep their size. */
.st-key-opendis_selections [data-testid="stMarkdownContainer"],
.st-key-opendis_selections [data-testid="stMarkdownContainer"] p,
.st-key-opendis_selections [data-testid="stMarkdownContainer"] span,
.st-key-opendis_selections [data-testid="stWidgetLabel"] p,
.st-key-opendis_selections [data-testid="stExpander"] summary p,
.st-key-opendis_selections [data-testid="stCheckbox"] label p,
.st-key-opendis_selections [role="radiogroup"] label p,
.st-key-opendis_selections button p {
  font-size: 0.95rem !important;
}
.st-key-opendis_selections [data-testid="stCaptionContainer"] p {
  font-size: 0.875rem !important;
}
/* The output column (the keyed container opendis_output): the simulation's screen output in a
   smaller monospace font, so that a line of PyDiS's fits */
.st-key-opendis_output [data-testid="stCode"] code,
.st-key-opendis_output [data-testid="stCode"] pre {
  font-size: 0.74rem !important; line-height: 1.35 !important;
}
.st-key-opendis_output [data-testid="stMarkdownContainer"] p,
.st-key-opendis_output button p {
  font-size: 0.95rem !important;
}
.st-key-opendis_output [data-testid="stCaptionContainer"] p {
  font-size: 0.875rem !important;
}
/* The steps in the sidebar: flush left, one line each, smaller text (after the rules above) */
section[data-testid="stSidebar"] div[data-testid="stRadio"],
section[data-testid="stSidebar"] div[data-testid="stRadio"] label {
  text-align: left !important; justify-content: flex-start !important;
}
section[data-testid="stSidebar"] div[data-testid="stRadio"] div[role="radiogroup"] {
  flex-direction: column !important; flex-wrap: nowrap !important;
  align-items: flex-start !important; justify-content: flex-start !important;
  gap: 0.2rem !important;
}
/* The Core menu's label, at the size of the steps (the keyed container opendis_core) */
section[data-testid="stSidebar"] .st-key-opendis_core [data-testid="stWidgetLabel"] p {
  font-size: 0.92rem !important;
}
/* The status of ExaDiS's build, at the size of the simulation output's text */
.st-key-opendis_build_status div[data-testid="stMarkdownContainer"] p,
.st-key-opendis_build_status [data-testid="stCaptionContainer"] p {
  font-size: 0.74rem !important; line-height: 1.35 !important;
}
section[data-testid="stSidebar"] div[data-testid="stRadio"] label p,
section[data-testid="stSidebar"] div[data-testid="stRadio"] label span {
  font-size: 0.92rem !important; white-space: nowrap !important;
}
/* The prev / next buttons below the steps, in the steps' size */
section[data-testid="stSidebar"] [data-testid="stButton"] button p {
  font-size: 0.92rem !important;
}
</style>
"""
# pylint: enable=line-too-long


def css() -> str:
    """The page CSS: pinned banner, sidebar below it, fonts, centred titles and radios."""
    return CSS


def png_data_uri(path: Path) -> str:
    """A PNG file as a data URI, for an ``<img>`` tag (``st.markdown`` cannot load local files)."""
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def banner_html(title: str, logo: Path) -> str:
    """The pinned top banner: the logo, then the title, centred together."""
    return (
        '<div class="opendis-banner">'
        f'<img class="opendis-logo" src="{png_data_uri(logo)}" alt="OpenDiS logo">'
        f"<h1>{html.escape(title)}</h1>"
        "</div>"
    )


def panel_title_html(text: str) -> str:
    """A centred serif panel title, one for each panel."""
    return f'<div class="opendis-panel-title">{html.escape(text)}</div>'
