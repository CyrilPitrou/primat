"""Guard the tracked documentation against staling relative to the code.

`README.md`, `tests/README.md`, `notebooks/README.md`, `docs/tutorials/index.md`
and the two generated parameter templates quote specific PRIMATConfig
defaults, specific `runfiles/primat_reference_run.py` parameter names/values,
and specific reference numbers. None of that is machine-checked by anything
else, so a config refactor can silently leave them wrong (this happened: the
docs used to cite an `n_temperature_table`/`sampling_nTOp` that no longer
exist). These tests assert the quoted facts still hold, so a future change
that breaks them fails a test instead of just leaving stale prose.

Scope note: `CLAUDE.md` is deliberately *not* read here. It is a local,
gitignored file (`.gitignore`), so it is not present in a public clone or in
CI and cannot be asserted against. Anything that must be enforced has to live
in a tracked file -- which is why the "Validation reference" numbers were
moved into `tests/README.md` and `tests/reference_values.py`, where the two
`*_matches_reference_constants` tests below pin them.

Every test in this module is a static file read -- no solve -- so the whole
file stays in the fast (`-m "not slow"`) lane.
"""
import ast
import os
import re

import pytest

from primat.config import PRIMATConfig

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def _read_text(path):
    """Read a text file as UTF-8, explicitly.

    README.md, tests/README.md, config.py and the run templates all contain
    non-ASCII physics characters (ν, ↔, →, σ, …). ``open()`` with no
    ``encoding`` uses the platform's locale default, which on Windows is
    cp1252 ("charmap") and raises ``UnicodeDecodeError`` on those bytes --
    so every doc-consistency read must pin UTF-8 to stay portable.
    """
    with open(path, encoding="utf-8") as f:
        return f.read()


def test_cprimat_version_matches_pyproject():
    """primat-c/include/config.h's CPRIMAT_VERSION must track pyproject.toml's version.

    The sync is performed by hand in the same commit as a version bump, with
    no other automated check. Parse both files directly
    (no import of a built/installed package) so this test only depends on
    the two source files staying in the same commit.
    """
    pyproject_path = os.path.join(REPO_ROOT, "pyproject.toml")
    pyproject_text = _read_text(pyproject_path)
    pyproject_match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', pyproject_text)
    assert pyproject_match, "version field not found in pyproject.toml"
    pyproject_version = pyproject_match.group(1)

    config_h_path = os.path.join(REPO_ROOT, "primat-c", "include", "config.h")
    config_h_text = _read_text(config_h_path)
    config_h_match = re.search(r'#define\s+CPRIMAT_VERSION\s+"([^"]+)"', config_h_text)
    assert config_h_match, "CPRIMAT_VERSION macro not found in primat-c/include/config.h"
    config_h_version = config_h_match.group(1)

    assert config_h_version == pyproject_version, (
        f"CPRIMAT_VERSION ({config_h_version!r}) in primat-c/include/config.h "
        f"is out of sync with pyproject.toml's version ({pyproject_version!r}); "
        "update both in the same commit."
    )


def test_save_nTOp_defaults_match_readme():
    """README's n<->p weak-rate section states save_nTOp/save_nTOp_thermal default True."""
    cfg = PRIMATConfig()
    assert cfg.save_nTOp is True
    assert cfg.save_nTOp_thermal is True


def _reference_run_options():
    """Parse MyOptions out of primat_reference_run.py without running it.

    The script performs an expensive multi-minute solve as a side effect of
    import, so we extract the literal dict via the AST instead of importing
    the module.
    """
    path = os.path.join(REPO_ROOT, "runfiles", "primat_reference_run.py")
    tree = ast.parse(_read_text(path), filename=path)
    # MyOptions references module-level names (e.g. "Omegabh2": omegabh2), so
    # literal_eval alone can't resolve it; evaluate against a namespace built
    # from this module's own simple top-level literal assignments instead of
    # importing the module (which would trigger its expensive solve()).
    namespace = {}
    my_options_node = None
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)):
            continue
        if node.targets[0].id == 'MyOptions':
            my_options_node = node.value
            continue
        try:
            namespace[node.targets[0].id] = ast.literal_eval(node.value)
        except ValueError:
            pass

    if my_options_node is None:
        raise AssertionError("MyOptions dict not found in primat_reference_run.py")
    code = compile(ast.Expression(body=my_options_node), filename=path, mode='eval')
    return eval(code, {}, namespace)


@pytest.mark.parametrize("key,expected", [
    ("sampling_temperature_per_decade", 2000),
    ("numerical_precision", 1e-10),
    ("sampling_nTOp_per_decade", 125),
    ("T_start_cosmo_MeV", 100.0),
])
def test_reference_run_params_match_published_settings(key, expected):
    """The reference-run settings tests/README.md publishes must exist verbatim.

    tests/README.md's "Validation reference" section states the exact
    high-precision settings its numbers were produced with; if
    primat_reference_run.py stops setting one of them, those numbers are no
    longer reproducible by the documented command."""
    options = _reference_run_options()
    assert key in options, f"{key!r} no longer in primat_reference_run.py's MyOptions"
    assert options[key] == expected


def test_reference_run_params_are_known_to_config():
    """Every MyOptions key must be a real PRIMATConfig field (catches silent typos)."""
    options = _reference_run_options()
    with _no_warning_context():
        PRIMATConfig(options)


@pytest.mark.parametrize("key", [
    "output_mc_file_prefix",
    "output_mc_covariance",
    "output_mc_correlation",
])
def test_readme_mc_key_names_are_real_params(key):
    """README's MC section quotes these DEFAULT_PARAMS
    key names verbatim; a rename that isn't mirrored in README would leave
    users following a documented option that raises an 'unknown parameter'
    warning instead of doing anything."""
    from primat.config import DEFAULT_PARAMS
    assert key in DEFAULT_PARAMS, f"{key!r} no longer a DEFAULT_PARAMS key"
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = _read_text(readme_path)
    assert key in readme_text, f"{key!r} no longer documented in README.md"


def test_readme_does_not_reference_old_mc_file_key():
    """output_mc_file was hard-renamed to output_mc_file_prefix (no deprecated
    alias, author decision -- primat is not on PyPI yet); README must not
    keep referencing the old name."""
    from primat.config import DEFAULT_PARAMS
    assert "output_mc_file" not in DEFAULT_PARAMS
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = _read_text(readme_path)
    # Match the old key as a whole word so "output_mc_file_prefix" (the
    # correct, current name) does not trip this assertion.
    assert not re.search(r'\boutput_mc_file\b(?!_prefix)', readme_text)


def test_streamlit_wheel_matches_pyproject_version():
    """requirements.txt's last line (the Streamlit Cloud deployment chain,
    see wheels/README.md) must point at a
    wheel file that (a) actually exists under wheels/ and (b) has the same
    version as pyproject.toml, or the public demo silently keeps serving an
    old build after a version bump."""
    pyproject_path = os.path.join(REPO_ROOT, "pyproject.toml")
    pyproject_text = _read_text(pyproject_path)
    pyproject_version = re.search(r'(?m)^version\s*=\s*"([^"]+)"', pyproject_text).group(1)

    requirements_path = os.path.join(REPO_ROOT, "requirements.txt")
    lines = [l.strip() for l in _read_text(requirements_path).splitlines() if l.strip()]
    wheel_line = lines[-1]
    assert wheel_line.startswith("./wheels/") and wheel_line.endswith(".whl"), (
        f"requirements.txt's last line is expected to be the Streamlit-Cloud "
        f"wheel path, got {wheel_line!r}"
    )
    wheel_filename = wheel_line[len("./wheels/"):]
    wheel_path = os.path.join(REPO_ROOT, "wheels", wheel_filename)
    assert os.path.isfile(wheel_path), f"{wheel_line!r} does not exist on disk"

    # Wheel filenames are "primat-<version>-<tag...>.whl" (PEP 427).
    assert wheel_filename.startswith(f"primat-{pyproject_version}-"), (
        f"wheels/{wheel_filename} does not match pyproject.toml's version "
        f"{pyproject_version!r} -- rebuild via build_linux.yml and update "
        f"requirements.txt (see wheels/README.md)."
    )


def test_readme_set_syntax_is_key_equals_value():
    """README's CLI section documents `--set KEY=VALUE`; primat/cli.py actually
    requires the '=' form (argparse splits on it), not the 'KEY VALUE' form
    README used to show -- which primat --set rejects."""
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = _read_text(readme_path)
    assert "--set KEY=VALUE" in readme_text
    assert "--set tau_n=880.1" in readme_text
    # The old (wrong) space-separated form must not reappear.
    assert "--set tau_n 880.1" not in readme_text


def test_readme_python_only_features_list_matches_backend():
    """README's 'Python-only features' list must match
    primat/backend.py's actual auto-fallback gate. The only
    inherently-Python run_bbn feature is background= (a custom Background
    object); extra_rho, decay_era, custom_network and output_time_evolution
    are all supported on the C backend (only MC prev is additionally
    Python-only, on the run_mc side)."""
    backend_path = os.path.join(REPO_ROOT, "primat", "backend.py")
    backend_text = _read_text(backend_path)
    # The feature actually gated in run_bbn()'s fallback logic -- if this
    # string disappears from backend.py, the module was refactored and
    # README's list needs re-verifying against the new code.
    assert "python_only_feature = background is not None" in backend_text

    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = _read_text(readme_path)
    assert "background=" in readme_text
    assert "MC `prev`" in readme_text
    # extra_rho/decay_era/custom_network/output_time_evolution must NOT be
    # listed as Python-only any more -- all supported on the C backend.
    python_only_section = readme_text[readme_text.index("Python-only features"):]
    python_only_section = python_only_section[:python_only_section.index("### Using primat-c directly")]
    assert "custom_network` (GUI" not in python_only_section
    assert "output_time_evolution=True (write full time series)" not in python_only_section
    # ...they appear only in the sentence explicitly excluding them from the
    # Python-only list.
    assert "extra_rho" in python_only_section and "decay_era" in python_only_section
    assert "all four are supported on the C" in python_only_section


def test_notebooks_readme_lists_every_notebook():
    """notebooks/README.md's table must mention every notebooks/*.ipynb file --
    a new notebook silently missing from the README
    is undiscoverable from the folder's own index."""
    notebooks_dir = os.path.join(REPO_ROOT, "notebooks")
    readme_text = _read_text(os.path.join(notebooks_dir, "README.md"))
    ipynb_files = sorted(f for f in os.listdir(notebooks_dir) if f.endswith(".ipynb"))
    assert ipynb_files, "no notebooks found in notebooks/"
    missing = [f for f in ipynb_files if f not in readme_text]
    assert not missing, f"notebooks/README.md is missing rows for: {missing}"


class _no_warning_context:
    """Fail the test if PRIMATConfig(options) emits an 'unknown parameter' warning."""

    def __enter__(self):
        import warnings
        self._cw = warnings.catch_warnings(record=True)
        self._records = self._cw.__enter__()
        warnings.simplefilter("always")
        return self

    def __exit__(self, *exc):
        self._cw.__exit__(*exc)
        unknown = [r for r in self._records if "unknown parameter" in str(r.message)]
        assert not unknown, f"PRIMATConfig reported unknown keys: {unknown}"


# ---------------------------------------------------------------------------
# DEFAULT_PARAMS three-file sync (config.py's DEFAULT_PARAMS <-> the
# generated primat_run_explanatory.py and run_basic.ini). This used to be a
# purely manual rule; these tests make drift a test failure.
# ---------------------------------------------------------------------------

_TEMPLATE_PY = os.path.join(REPO_ROOT, "runfiles", "primat_run_explanatory.py")
_TEMPLATE_INI = os.path.join(REPO_ROOT, "primat-c", "examples", "run_basic.ini")


def _template_py_keys():
    """Keys listed in the Python template's commented-out cfg dict."""
    text = _read_text(_TEMPLATE_PY)
    # Entries look like "    # numerical_precision=1e-7,  # comment" -- key
    # immediately followed by '=' (prose comment lines never match this).
    return set(re.findall(r"^\s*#\s*([A-Za-z_][A-Za-z0-9_]*)=", text, re.M))


def _template_ini_keys():
    """Keys listed in the INI template (commented or active settings)."""
    text = _read_text(_TEMPLATE_INI)
    # Entries look like "# verbose = false" or "numerical_precision = 1e-7".
    # The INI template ALSO carries prose lines that explain a boolean's
    # meaning in the same "# word = text" shape, e.g.
    # "# false = instantaneous decoupling" / "# true = Born approximation";
    # those capture the literal value word, never a parameter, so drop the
    # boolean/None literals (no DEFAULT_PARAMS key is ever named that).
    keys = set(re.findall(r"^#?\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", text, re.M))
    return keys - {"true", "false", "none"}


def test_explanatory_template_lists_every_default_param():
    """runfiles/primat_run_explanatory.py lists every DEFAULT_PARAMS key, and no
    stale ones -- it is the primary user-facing reference for what can be set."""
    from primat.config import DEFAULT_PARAMS
    found, expected = _template_py_keys(), set(DEFAULT_PARAMS)
    assert expected - found == set(), \
        f"primat_run_explanatory.py is missing keys: {sorted(expected - found)}"
    assert found - expected == set(), \
        f"primat_run_explanatory.py lists stale/unknown keys: {sorted(found - expected)}"


def test_ini_template_lists_every_default_param():
    """primat-c/examples/run_basic.ini lists every DEFAULT_PARAMS key, and no
    stale ones -- the C backend's equivalent of the template above."""
    from primat.config import DEFAULT_PARAMS
    found, expected = _template_ini_keys(), set(DEFAULT_PARAMS)
    assert expected - found == set(), \
        f"run_basic.ini is missing keys: {sorted(expected - found)}"
    assert found - expected == set(), \
        f"run_basic.ini lists stale/unknown keys: {sorted(found - expected)}"


def test_param_templates_match_generator():
    """Both templates are build products of primat/tools/gen_param_templates.py
    (see that module's docstring): this asserts the committed files are
    byte-for-byte what the generator would produce right now, so a
    DEFAULT_PARAMS/PARAM_GROUPS change that isn't followed by
    `python -m primat.tools.gen_param_templates` fails here instead of
    silently drifting (the "Keeping DEFAULT_PARAMS ... in sync" chore, now
    enforced rather than merely documented)."""
    from primat.tools.gen_param_templates import (
        generate_run_explanatory, generate_run_basic_ini)

    assert _read_text(_TEMPLATE_PY) == generate_run_explanatory(), \
        "runfiles/primat_run_explanatory.py is stale -- regenerate with " \
        "`python -m primat.tools.gen_param_templates`"
    assert _read_text(_TEMPLATE_INI) == generate_run_basic_ini(), \
        "primat-c/examples/run_basic.ini is stale -- regenerate with " \
        "`python -m primat.tools.gen_param_templates`"


def test_param_count_comments_match_default_params():
    """The two templates' count comments must quote len(DEFAULT_PARAMS)
    exactly. (Any '(currently NN keys)' count in the untracked CLAUDE.md is
    NOT asserted -- see this module's docstring; keep it updated by hand.)"""
    from primat.config import DEFAULT_PARAMS
    n = len(DEFAULT_PARAMS)
    assert f"All {n} DEFAULT_PARAMS keys are listed" in _read_text(_TEMPLATE_PY)
    assert f"all {n} keys round-trip" in _read_text(_TEMPLATE_INI)


def _validation_reference_section():
    """The text of tests/README.md's 'Validation reference' section."""
    readme = _read_text(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "README.md"))
    return readme[readme.index("## Validation reference"):]


def test_validation_reference_table_matches_reference_constants():
    """tests/README.md's 'Validation reference' observable tables and
    tests/reference_values.py must quote the same numbers, so neither can be
    updated without the other (the regression tier asserts the constants
    against actual solves)."""
    from reference_values import (REF_SMALL_YPBBN, REF_SMALL_DOH,
                                  REF_LARGE8_YPBBN, REF_LARGE8_DOH)
    section = _validation_reference_section()
    yp = re.findall(r"\|\s*YP \(BBN\)\s*\|\s*([0-9.eE+-]+)\s*\|", section)
    dh = re.findall(r"\|\s*D/H\s*\|\s*([0-9.eE+-]+)\s*\|", section)
    # First table is the small network, second is large+amax=8.
    assert [float(v) for v in yp[:2]] == [REF_SMALL_YPBBN, REF_LARGE8_YPBBN]
    assert [float(v) for v in dh[:2]] == [REF_SMALL_DOH, REF_LARGE8_DOH]


def test_per_nuclide_reference_table_matches_reference_constants():
    """tests/README.md's per-nuclide table and reference_values.NUCLIDE_REFERENCE
    must quote the same numbers, in the same column order.

    GOAL: close the half of the "Validation reference" section that nothing
    parsed. Until now only the two observable tables above were checked, so
    the 21-cell per-nuclide table could -- and did -- drift silently: by
    2026-08-05 every one of its rows was stale in the 5th significant figure
    (up to 6.2e-05 relative on the free neutron), while being advertised as
    test-pinned.

    This is the static half (README text vs. the constants);
    tests/test_regression.py::test_per_nuclide_abundances_match_the_reference_table
    is the live half (the constants vs. an actual solve).
    """
    from reference_values import NUCLIDE_REFERENCE, NUCLIDE_COLUMNS
    section = _validation_reference_section()
    for nuclide, expected in NUCLIDE_REFERENCE.items():
        m = re.search(rf"^\|\s*{nuclide}\s*\|([^\n]*)\|\s*$", section, re.M)
        assert m, f"tests/README.md's per-nuclide table has no {nuclide!r} row"
        cells = [float(c) for c in m.group(1).split("|")]
        assert len(cells) == len(NUCLIDE_COLUMNS), (
            f"{nuclide!r} row has {len(cells)} columns, "
            f"expected {NUCLIDE_COLUMNS}")
        assert cells == list(expected), (
            f"tests/README.md's {nuclide!r} row {cells} disagrees with "
            f"reference_values.NUCLIDE_REFERENCE {list(expected)}")


def test_readme_documents_the_unified_evolution_schema():
    """README's Output section must describe the real unified TSV header
    (t_s/a/T_gamma_MeV/...), not the pre-unification legacy column list."""
    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    for col in ("`t_s`", "T_gamma_MeV", "T_nutau_MeV"):
        assert col in readme, f"README lost the unified-schema column {col}"
    assert "n_to_p_weak_rate" not in readme  # the legacy column list


def test_readme_result_dict_table_is_complete():
    """Every key of a real run_bbn() result dict appears in README."""
    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    # Keep in sync with cpr_assemble_results/_python_solve; cheap static list
    # (a live solve here would drag this file into the slow tier).
    for key in ("YPBBN", "YPCMB", "He4oH", "DoH", "He3oH", "He3oHe4", "Li7oH",
                "Neff", "Omeganurel", "OneOverOmeganunr", "Y_final"):
        assert f"`{key}`" in readme, f"README result-dict table is missing {key}"


def test_readme_gui_is_not_called_source_only():
    """primat.gui ships in the wheel ([tool.setuptools] packages); only
    runfiles/ is source-only."""
    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    assert "Streamlit app (optional, source-only)" not in readme


def test_readme_rate_columns_match_implementation():
    """README's output_rates_time_evolution claims must describe the
    implemented per-reaction forward-rate column block, not the
    historical no-op."""
    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    assert "_frwrd" in readme
    assert "no-op" not in readme


def test_tests_readme_lists_every_test_file():
    """tests/README.md's structure table must mention every tests/test_*.py
    -- the suite's README documents the goal of every test group.
    Mirrors test_notebooks_readme_lists_every_notebook."""
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    readme_text = _read_text(os.path.join(tests_dir, "README.md"))
    files = sorted(f for f in os.listdir(tests_dir)
                   if f.startswith("test_") and f.endswith(".py"))
    assert files, "no test files found"
    missing = [f for f in files if f not in readme_text]
    assert not missing, f"tests/README.md is missing rows for: {missing}"


def test_docs_tutorials_gallery_complete():
    """Every notebook in notebooks/ is rendered in the docs tutorial
    gallery: symlinked into docs/tutorials/ and listed in its toctree.

    Guards against the drift found 2026-07-10 where ReactionRates.ipynb and
    AnimatedAbundances.ipynb were added to notebooks/ but never surfaced on
    the docs site."""
    notebooks_dir = os.path.join(REPO_ROOT, "notebooks")
    tutorials_dir = os.path.join(REPO_ROOT, "docs", "tutorials")
    notebooks = {os.path.splitext(f)[0] for f in os.listdir(notebooks_dir)
                 if f.endswith(".ipynb")}
    index = _read_text(os.path.join(tutorials_dir, "index.md"))
    opt_out = set()  # notebooks intentionally kept off the docs site
    for stem in sorted(notebooks - opt_out):
        assert os.path.exists(os.path.join(tutorials_dir, f"{stem}.ipynb")), stem
        assert stem in index, f"{stem} missing from docs/tutorials/index.md"


@pytest.mark.parametrize("page", ["docs/index.md", "README.md"])
def test_quick_start_numbers_match_the_reference_constants(page):
    """Both landing pages' quick starts must quote the tracked reference values.

    GOAL: each prints YP and D/H for both backends, to the decimal counts the
    project reports at. docs/index.md went stale twice, by 8.6e-09 in D/H once
    and by 2.1e-07 in YP after that was hand-corrected; the README carried the
    same 8.6e-09 error for longer, because only one of the two was checked.
    The two pages quote the same numbers on purpose -- the README is meant to
    stand alone -- so both are pinned here. Compared against
    tests/reference_values.py rather than a fresh solve, for the same reason
    the tests/README.md table is: the eighth decimal is not portable. Measured
    on one commit across the CI matrix, the C backend's YPBBN is 0.24699907 on
    macOS/arm64 and 0.24699900 on Linux and Windows/x86_64 -- a 6.7e-08 spread,
    only 3x below the drift this test exists to catch, so a live-run bound
    could not separate the two. The constants themselves are pinned to live
    solves by test_regression.py.
    """
    from reference_values import (DOH_REFERENCE, PY_DOH_REFERENCE,
                                  PY_YPBBN_REFERENCE, YPBBN_REFERENCE)

    text = _read_text(os.path.join(REPO_ROOT, *page.split("/")))
    yp_c = re.search(r"YPBBN'\]:\.8f\}\"\)\s*#\s*([0-9.]+)", text).group(1)
    doh_c = re.search(r"DoH'\]:\.7e\}\"\)\s*#\s*([0-9.e+-]+)", text).group(1)
    yp_py, doh_py = re.search(
        r"pure-Python backend gives `([0-9.]+)` / `([0-9.e+-]+)`", text).groups()

    # Exact strings, not approx: both sides are tracked files, so a mismatch
    # means one was edited without the other.
    assert yp_c == f"{YPBBN_REFERENCE:.8f}"
    assert doh_c == f"{DOH_REFERENCE:.7e}"
    assert yp_py == f"{PY_YPBBN_REFERENCE:.8f}"
    assert doh_py == f"{PY_DOH_REFERENCE:.7e}"


def test_development_notes_only_point_at_files_that_exist():
    """docs/development.md's cross-references must all resolve.

    GOAL: apply that file's own "no unverifiable claims about other files"
    rule to itself. It routes the reader to the test that enforces each
    project rule, so a renamed or deleted test would send them nowhere -- the
    exact decay the rule exists to prevent.
    """
    text = _read_text(os.path.join(REPO_ROOT, "docs", "development.md"))

    referenced = set(re.findall(r"`(tests/[\w./]+\.py)`", text))
    referenced |= set(re.findall(r"`(primat(?:-c)?/[\w./-]+\.(?:py|c|h))`", text))
    referenced |= {"README.md", "tests/README.md"}
    assert len(referenced) >= 8, f"suspiciously few references parsed: {referenced}"

    missing = [r for r in sorted(referenced)
               if not os.path.exists(os.path.join(REPO_ROOT, r))]
    assert not missing, (
        f"docs/development.md points at files that do not exist: {missing}")

    # The two section titles it sends the reader to must still be there.
    assert "Backend parity contract" in _read_text(
        os.path.join(REPO_ROOT, "README.md"))
    assert "Known cross-backend divergences" in _read_text(
        os.path.join(REPO_ROOT, "tests", "README.md"))
