"""
Notebook smoke test.

Executes demonstration notebooks end-to-end with ``papermill`` and checks
they run without raising. Two tiers:

* ``FAST_NOTEBOOKS`` -- no Monte Carlo, run as-is:
    - ``AbundanceEvolution.ipynb`` -- small/large(amax=8)/large solves, ~5 s.
    - ``CompareSmallNetworks.ipynb`` -- two small-network solves, ~4 s.
    - ``ReactionRates.ipynb`` -- one small-network build; plots the weak and
      nuclear reaction rates vs the Hubble rate, no solve(), ~4 s.
* ``MC_NOTEBOOKS`` -- normally run a Monte Carlo scan at publication-quality
  sample counts (``num_mc``/``N_MC`` ~100-500, sometimes over a parameter
  grid too, e.g. ``StandardPlots.ipynb``'s 20 eta points x 100 MC samples).
  The cell that sets the sample count is tagged ``parameters`` in each of
  these notebooks, so papermill overrides it down to ``MC_NOTEBOOK_NUM_MC``
  (3) here -- enough to exercise the MC code path (`run_bbn`/`run_mc`
  wiring, plotting of central value + error band) without paying for
  publication-quality statistics:
    - ``AbundancesNrelat.ipynb``, ``AbundancesXi.ipynb`` -- ~21/11-point
      grids x 3 MC samples.
    - ``PosteriorBaryons.ipynb`` -- ~17-point grid x 3 MC samples.
    - ``StandardPlots.ipynb`` -- 20-point grid x 3 MC samples.
    - ``MonteCarloRates.ipynb`` -- 3 full-BBN MC samples (no grid).

This is a regression guard against import-path bugs (the notebooks still
imported the pre-reorganisation ``pypr`` package name)
and against API drift in ``primat.main.PRIMAT``/``primat.backend.run_mc``: a
renamed/removed attribute that the notebooks rely on (``r.A``,
``r.abundance_names``, ``r[name](t)``, ``run_mc(...)``'s return shape, ...)
makes one of these cells raise, and papermill re-raises that as a
``CellExecutionError`` here.

Two notebooks live outside ``notebooks/`` and derive their root from the
working directory, so each is executed from a throwaway copy of its own
directory: ``generate_rates/thermal_average.ipynb`` (the rate-generation
how-to) and ``manual/primat_doc_figures.ipynb`` (the manual's figure set).
Every notebook in the repository is therefore covered; ``AbundancesXi``/
``AbundancesNrelat``'s MC scans run at reduced sample counts, as above.

Requires ``papermill`` (an optional ``notebooks`` extra, see
``pyproject.toml``); skipped if not installed.
"""
import shutil
from pathlib import Path

import pytest

pytestmark = [pytest.mark.slow, pytest.mark.notebook]

NOTEBOOKS_DIR = Path(__file__).resolve().parents[1] / "notebooks"

FAST_NOTEBOOKS = [
    "AbundanceEvolution.ipynb",
    "CompareSmallNetworks.ipynb",
    "ReactionRates.ipynb",
    # No Monte Carlo either, so no parameter override is needed: Sensitivity
    # runs the derivative scan (~6 s) and AnimatedAbundances builds its frames
    # from 40 solves (21 for the DeltaNeff sweep, 19 for Omega_b h^2).
    "Sensitivity.ipynb",
    "AnimatedAbundances.ipynb",
]

# Notebook name -> papermill parameter dict overriding its MC sample count
# (the cell tagged "parameters" in each of these notebooks).
MC_NOTEBOOK_NUM_MC = 3
MC_NOTEBOOKS = {
    "AbundancesNrelat.ipynb": {"num_mc": MC_NOTEBOOK_NUM_MC},
    "AbundancesXi.ipynb": {"num_mc": MC_NOTEBOOK_NUM_MC},
    "PosteriorBaryons.ipynb": {"num_mc": MC_NOTEBOOK_NUM_MC},
    "StandardPlots.ipynb": {"num_mc": MC_NOTEBOOK_NUM_MC},
    "MonteCarloRates.ipynb": {"N_MC": MC_NOTEBOOK_NUM_MC},
}


def _run_notebook(name, tmp_path, monkeypatch, parameters=None):
    """Execute one notebook with papermill from a throwaway copy of notebooks/."""
    papermill = pytest.importorskip("papermill")

    # Headless plotting backend: notebooks call plt.savefig()/plt.show(),
    # which would otherwise try (and fail) to open a GUI window in CI.
    monkeypatch.setenv("MPLBACKEND", "Agg")

    # Run from a throwaway copy of notebooks/, so plt.savefig('plots/...')
    # (a relative path resolved against cwd) writes into tmp_path instead
    # of overwriting the tracked PDFs in notebooks/plots/.
    work_dir = tmp_path / "notebooks"
    shutil.copytree(NOTEBOOKS_DIR, work_dir)

    # The notebooks save figures with plt.savefig('plots/<name>') -- a path
    # relative to cwd (== work_dir). notebooks/plots/ is .gitignore'd (its PDFs
    # are regenerable output, not tracked), so a fresh CI checkout has no
    # plots/ subdir for copytree to copy, and savefig() raises
    # FileNotFoundError. Create it here so the test does not depend on the
    # author's local, untracked plots/ directory being present.
    (work_dir / "plots").mkdir(exist_ok=True)

    papermill.execute_notebook(
        str(work_dir / name), str(work_dir / f"out_{name}"),
        cwd=str(work_dir),
        progress_bar=False,
        parameters=parameters or {},
    )


@pytest.mark.parametrize("name", FAST_NOTEBOOKS)
def test_fast_notebook_executes(name, tmp_path, monkeypatch):
    """Run a fast demo notebook with papermill; fail if any cell raises."""
    _run_notebook(name, tmp_path, monkeypatch)


@pytest.mark.parametrize("name", MC_NOTEBOOKS)
def test_mc_notebook_executes_with_few_samples(name, tmp_path, monkeypatch):
    """Run an MC demo notebook with its sample count cut to 3, via papermill
    parameter injection into the notebook's tagged "parameters" cell."""
    _run_notebook(name, tmp_path, monkeypatch, parameters=MC_NOTEBOOKS[name])


def test_thermal_average_notebook_executes(tmp_path, monkeypatch):
    """Run ``generate_rates/thermal_average.ipynb``, the rate-generation how-to.

    It is the only notebook whose self-checks are assertions rather than plots: three analytic thermal
    averages (constant sigma, 1/v, and the two d+d S-factor channels) that fail
    the cell if the quadrature drifts. It reads primat's shipped rate tables
    and reruns a small BBN solve with the table it just wrote, so an API rename
    in ``network_data``/``backend`` breaks it exactly as it breaks the others.

    The notebook derives the repo root from ``Path.cwd().parent`` and writes
    into ``<root>/generate_rates/rate_tables_out/``, so it is run from a
    throwaway ``generate_rates/`` whose parent links back to the real package
    -- output lands in tmp_path, inputs come from the checkout.
    """
    papermill = pytest.importorskip("papermill")
    monkeypatch.setenv("MPLBACKEND", "Agg")

    repo = Path(__file__).resolve().parents[1]
    work_root = tmp_path / "repo"
    (work_root / "generate_rates").mkdir(parents=True)
    try:
        (work_root / "primat").symlink_to(repo / "primat", target_is_directory=True)
    except (OSError, NotImplementedError):
        # Windows needs a privilege for symlinks; a copy costs ~24 MB and
        # works everywhere.
        shutil.copytree(repo / "primat", work_root / "primat")
    work_dir = work_root / "generate_rates"
    shutil.copy(repo / "generate_rates" / "thermal_average.ipynb", work_dir)

    papermill.execute_notebook(
        str(work_dir / "thermal_average.ipynb"),
        str(work_dir / "out_thermal_average.ipynb"),
        cwd=str(work_dir),
        progress_bar=False,
    )


def test_manual_figure_notebook_executes(tmp_path, monkeypatch):
    """Run ``manual/primat_doc_figures.ipynb``, which regenerates every figure
    in the LaTeX manual.

    It is the only notebook that reaches into a solved run's *intermediate*
    state -- ``background.Tg_vec``, ``plasma.rho_e``, the weak-rate
    interpolants -- so it fails on an attribute rename that the result-dict
    based notebooks would not notice. Run from a throwaway ``manual/`` whose
    parent links back to the real package, so its figures land in tmp_path.
    """
    papermill = pytest.importorskip("papermill")
    monkeypatch.setenv("MPLBACKEND", "Agg")

    repo = Path(__file__).resolve().parents[1]
    work_root = tmp_path / "repo"
    (work_root / "manual").mkdir(parents=True)
    try:
        (work_root / "primat").symlink_to(repo / "primat", target_is_directory=True)
    except (OSError, NotImplementedError):
        # Windows needs a privilege for symlinks; a copy works everywhere.
        shutil.copytree(repo / "primat", work_root / "primat")
    work_dir = work_root / "manual"
    shutil.copy(repo / "manual" / "primat_doc_figures.ipynb", work_dir)

    papermill.execute_notebook(
        str(work_dir / "primat_doc_figures.ipynb"),
        str(work_dir / "out_primat_doc_figures.ipynb"),
        cwd=str(work_dir),
        progress_bar=False,
    )
