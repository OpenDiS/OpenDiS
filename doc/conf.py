# Configuration file for the Sphinx documentation builder .
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

#sys.path.extend([os.path.abspath('../python'),os.path.abspath('../lib'),os.path.abspath('../core/pydis/python')])



# -- Project information -----------------------------------------------------

project = 'OpenDiS'
author = 'OpenDiS'


from datetime import datetime
year = datetime.now().year
if year > 2024:
    copyright = f'2024-{year}, {author}'
else:
    copyright = f'2024, {author}'

# does not work
#copyright = '2024, OpenDiS | [contributors](contributors.md)'

# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser','sphinx_inline_tabs','sphinx_copybutton','sphinx_favicon',]


# templates_path = ['_templates']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme_options = {
#     "announcement": " ",
# }

# # change theme color 
# html_theme_options = {
#     "light_css_variables": {
#         "color-brand-primary": "#7C4DFF",
#         "color-brand-content": "#7C4DFF",
#     },
#     "show_breadcrumbs": True
# }

html_theme_options = {
    "footer_icons": [
        {
            "name": "Contributors",
            "url": "https://caiwei-stanford.github.io/opendis-doc/contributors.html",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
}

html_static_path = ['_static']

def setup(app):
    app.add_css_file('custom.css')
    app.add_js_file('custom.js')
    

html_logo = "nodebox_logo.png"
html_theme =  "furo"

#favicons = ["logo_box_upright.png",]
#favicons = [{"href": "favicon.ico"},]
favicons = ["favicon_node.ico",]
