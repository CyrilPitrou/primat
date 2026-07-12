# -*- coding: utf-8 -*-
"""
Sphinx configuration for the primat documentation site.

This is the *scaffolding* half of FABLEADVICE O-3 (Opus): it fixes the stack
and information architecture so the bulk content migration (Sonnet's half) has
a clean, ``sphinx-build -W``-passing skeleton to fill in.

Stack (all pinned in the ``docs`` extra of ``pyproject.toml``):

* ``furo``              -- clean, responsive, light/dark HTML theme.
* ``myst-parser``       -- lets us author pages in Markdown (the repo's existing
                           README/EXTENDING/notebook-README content is all .md),
                           so migration is copy-paste rather than reStructuredText
                           rewriting.
* ``myst-nb``           -- renders the ``notebooks/`` gallery as tutorial pages.
                           Execution is *off* by default here (see
                           ``nb_execution_mode``); the notebooks ship with stored
                           outputs and are only re-executed in CI's nightly lane.
* ``sphinx.ext.autodoc``+``napoleon`` -- API reference straight from primat's
                           (already excellent) NumPy/Google-style docstrings.
* ``sphinx.ext.intersphinx`` -- cross-links to numpy/scipy/python object pages.
* ``sphinx-copybutton`` -- one-click copy on every code block.
* ``sphinxarg.ext``     -- auto-generates the CLI reference from
                           ``primat.cli._build_parser`` (no hand-maintained flag
                           list that can drift from ``--help``).

Build locally with::

    pip install -e ".[docs]"
    sphinx-build -W -b html docs docs/_build/html

The ``-W`` (warnings-as-errors) invocation is what CI runs, so keep the tree
warning-clean: every ``toctree`` entry must resolve to an existing page, and
autodoc must be able to import (or mock) every documented module.
"""

from __future__ import annotations

from importlib.metadata import version as _pkg_version

# -- Project information -----------------------------------------------------

project = "primat"
author = "Cyril Pitrou, Alain Coc, Jean-Philippe Uzan, Elisabeth Vangioni"
copyright = "2018-2026, the primat authors"

# Single source of truth for the version is ``pyproject.toml`` (read back via
# the installed package metadata), so the docs never drift from the package.
release = _pkg_version("primat")
# The short X.Y version shown in the sidebar.
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    # NOTE: load only ``myst_nb`` here, NOT ``myst_parser`` as well -- myst-nb
    # already sets up myst-parser internally, and listing both double-registers
    # its roles/transforms (which errors out on newer Sphinx). myst-nb gives us
    # the full MyST Markdown parser plus notebook rendering in one extension.
    "myst_nb",
    "sphinx_copybutton",
    "sphinxarg.ext",
]

# Source files: Markdown (MyST) for prose, notebooks for tutorials. Note we do
# NOT list ``.rst`` here because the whole site is authored in Markdown; add it
# back if a hand-written .rst page is ever needed.
source_suffix = {
    ".md": "myst-nb",   # let myst-nb own .md too so MyST directives work uniformly
    ".ipynb": "myst-nb",
}

master_doc = "index"
language = "en"

# Files/dirs Sphinx should ignore when scanning ``docs/``.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "**.ipynb_checkpoints",
    # Development-process artifacts (plans/specs written by the "superpowers"
    # Claude Code skill while implementing a feature) that got committed
    # under docs/ by mistake -- they document the *process*, not the
    # product, and were never meant to be toctree pages. Excluding the
    # directory (rather than deleting the files, which are still useful
    # implementation history) keeps `sphinx-build -W` warning-clean without
    # touching their content.
    "superpowers",
]

# -- MyST (Markdown) options -------------------------------------------------

# Enable the commonly-needed MyST extensions. ``dollarmath``/``amsmath`` let us
# write the BBN equations inline ($...$) and in display blocks; ``colon_fence``
# lets ``:::{note}`` admonitions work; ``deflist`` renders the parameter tables
# nicely; ``linkify`` turns bare URLs into links.
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
    "fieldlist",
    "linkify",
    "substitution",
    "tasklist",
]
# Auto-generate anchors for headings up to level 3 so intra-doc links like
# ``[...](networks.md#choosing-amax)`` resolve.
myst_heading_anchors = 3

# -- myst-nb (notebook) options ----------------------------------------------

# Do NOT execute notebooks during the ordinary ``-W`` docs build: they are slow
# and some need optional extras (numba/vegas/papermill). The committed notebooks
# ship with stored outputs, which myst-nb renders as-is. CI's nightly lane is
# the place to re-execute them (set ``NB_EXECUTION_MODE=cache`` there).
import os as _os

nb_execution_mode = _os.environ.get("NB_EXECUTION_MODE", "off")
nb_execution_timeout = 900  # seconds, only relevant when execution is enabled
nb_execution_allow_errors = False

# -- autodoc / autosummary ---------------------------------------------------

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_rtype = False
# Render a class's ``Attributes`` section as inline ``:ivar:`` fields inside the
# class docstring rather than as separate ``py:attribute`` directives. Without
# this, a dataclass documents each attribute twice -- once from the annotated
# field (``:members:``) and once from the napoleon ``Attributes`` block -- which
# is a duplicate-object-description warning (fatal under ``-W``).
napoleon_use_ivar = True

# Optional/heavy third-party deps that primat imports lazily. Mock them so
# ``autodoc`` can import primat's modules on a docs-only environment (RTD, the
# CI docs job) without pulling in numba/vegas/plotly/streamlit toolchains.
autodoc_mock_imports = [
    "numba",
    "vegas",
    "joblib",
    "plotly",
    "streamlit",
    "pandas",
    "matplotlib",
    "papermill",
]

# -- intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- HTML output -------------------------------------------------------------

html_theme = "furo"
html_title = f"primat {version}"
html_static_path = ["_static"]
# Point the "edit this page" / repo links at GitHub.
html_theme_options = {
    "source_repository": "https://github.com/CyrilPitrou/primat",
    "source_branch": "master",
    "source_directory": "docs/",
}

# -- copybutton --------------------------------------------------------------

# Strip common shell/REPL prompts when copying so pasted commands are clean.
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
