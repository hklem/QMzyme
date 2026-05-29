# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime
import subprocess
# Ensure that modules can be imported without installing aqme
sys.path.insert(0, os.path.abspath('..')) 



# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

now = datetime.datetime.now()
project = 'QMzyme Documentation'
copyright = f'{now.year}, Heidi Klem'
author = 'Heidi Klem'
packageversion = subprocess.check_output(
        ["git", "describe", "--tags", "--abbrev=0"],
        text=True
    ).strip()
release = packageversion
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', 'nbsphinx', 'sphinx_copybutton']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Disable  smartquotes which might transform '--' into a different character
smartquotes = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': True,
}
html_static_path = ['_static']
html_logo = "_static/QMzyme_logo.png"
autoclass_content = 'both'
toc_object_entries = False
autodoc_member_order = 'bysource'

# Custom sidebar templates, maps document names to template names.
# alabaster sidebars
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',
        'searchbox.html',
    ]
}

rst_prolog = """
.. role:: specialcode
   :class: specialcode
"""

def setup(app):
    app.add_css_file('custom.css')

# Set the syntax highlighting style to the Pygments default (matches standard Jupyter)
pygments_style = 'default'
nbsphinx_codecell_lexer = 'python3'

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True