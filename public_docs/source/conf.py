# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

# Add project root to path for autodoc
sys.path.insert(0, os.path.abspath("../../quasix"))

# Mock imports for modules that aren't available on ReadTheDocs
# This allows autodoc to work without building the Rust extension
autodoc_mock_imports = [
    "quasix.quasix",  # The Rust extension module
    "pyscf",
    "h5py",
    "matplotlib",
]

# -- Project information -----------------------------------------------------
project = "QuasiX"
copyright = f"2024-{datetime.now().year}, QuasiX Development Team"
author = "Viwat Vchirawongkwin"
release = "0.6.0"
version = "0.6.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.ifconfig",
    "sphinx.ext.githubpages",
    "nbsphinx",
    "myst_parser",  # For including markdown files
]

# Add myst-parser for markdown support
myst_enable_extensions = [
    "dollarmath",  # Support $...$ and $$...$$ math syntax
    "colon_fence",
    "deflist",
    "html_image",
]

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

# The master document
master_doc = "index"

# Language
language = "en"

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
}

# HTML context for GitHub integration
html_context = {
    "display_github": True,
    "github_user": "ExaPsi",
    "github_repo": "QuasiX",
    "github_version": "main",
    "conf_py_path": "/public_docs/source/",
}

# Custom sidebar templates
html_sidebars = {
    "**": [
        "relations.html",
        "searchbox.html",
        "globaltoc.html",
    ]
}

# -- Options for autodoc -----------------------------------------------------
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "special-members": "__init__",
    "undoc-members": True,
    "exclude-members": "__weakref__",
}

# -- Options for napoleon (Google/NumPy docstrings) --------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# -- Options for intersphinx -------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pyscf": ("https://pyscf.org/", None),
}

# -- Options for MathJax -----------------------------------------------------
mathjax3_config = {
    "tex": {
        "macros": {
            "GW": r"\mathrm{GW}",
            "BSE": r"\mathrm{BSE}",
            "QP": r"\mathrm{QP}",
            "Ha": r"\mathrm{Ha}",
            "eV": r"\mathrm{eV}",
        },
        "inlineMath": [["$", "$"], ["\\(", "\\)"]],
        "displayMath": [["$$", "$$"], ["\\[", "\\]"]],
    }
}

# -- Options for nbsphinx ----------------------------------------------------
nbsphinx_execute = "never"  # Don't execute notebooks during build
nbsphinx_allow_errors = True

# -- Options for todo extension ----------------------------------------------
todo_include_todos = False  # Set to True during development

# -- Options for linkcheck ---------------------------------------------------
linkcheck_ignore = [
    r"http://localhost:\d+/",
]

# -- Custom CSS --------------------------------------------------------------
def setup(app):
    # Add custom CSS if _static/custom.css exists
    if os.path.exists(os.path.join(app.srcdir, "_static", "custom.css")):
        app.add_css_file("custom.css")
