# -*- coding: utf-8 -*-
"""
primat — core package for the PRIMAT BBN solver.

Public API::

    from primat import PRIMAT
    result = PRIMAT({"Omegabh2": 0.02242}).solve()
"""

from importlib.metadata import version as _version, PackageNotFoundError

from .main import PRIMAT, mc_uncertainty
from .background import Background, StandardBackground
from .nuclear_network import NuclearNetwork
from .network_data import nuclide_latex
from .credits import CITATION_BIBTEX as __citation__

# Discoverability aliases: `primat.backend` remains the
# canonical import used throughout docs/notebooks/examples, but IDE users
# exploring the top-level `primat` package get these for free too. Safe to
# import eagerly here since `backend.py` only imports `.main` lazily inside
# its functions, so there is no import cycle.
from .backend import run_bbn, run_mc, HAS_C_BACKEND
from .sensitivity import sensitivity_table, SensitivityTable, SensTarget

# Single source of truth for the version is pyproject.toml; we read it back
# from the installed distribution metadata so the number is never duplicated.
try:
    __version__ = _version("primat")
except PackageNotFoundError:
    # Running from a source checkout that was never installed (e.g. no
    # `pip install -e .`): metadata is absent, so fall back to a sentinel.
    __version__ = "0.0.0+unknown"

__all__ = ["PRIMAT", "mc_uncertainty", "Background", "StandardBackground",
           "NuclearNetwork", "nuclide_latex", "__version__", "__citation__",
           "run_bbn", "run_mc", "HAS_C_BACKEND",
           "sensitivity_table", "SensitivityTable", "SensTarget"]
