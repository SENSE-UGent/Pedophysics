# Configuration file for the Sphinx documentation builder.

from __future__ import annotations

import os
import sys

# -- Path setup --------------------------------------------------------------
# docs/source/conf.py -> add project src/ so autodoc can import pedophysics
sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------------
project = "Pedophysics"
copyright = "2024, Gaston Mendoza Veirana"
author = "Gaston Matias Mendoza Veirana"

# The full version, including alpha/beta/rc tags
release = "0.1"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
html_theme = "alabaster"

# If you don't have these folders yet, Sphinx may warn; we'll create them below.
html_static_path = ["_static"]
