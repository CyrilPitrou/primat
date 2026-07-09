"""Guard against README/CLAUDE.md documentation staling relative to the code.

Both docs quote specific PRIMATConfig defaults and specific
runfiles/primat_reference_run.py parameter names/values (CLAUDE.md's
"Validation before committing" section says references were produced with
particular settings). Neither file is machine-checked by anything else, so a
config refactor can silently leave them wrong (this happened: CLAUDE.md used
to cite a `n_temperature_table`/`sampling_nTOp` that no longer exist). These
tests assert the quoted facts still hold, so a future config change that
breaks them fails a test instead of just leaving stale prose.
"""
import ast
import os
import re

import pytest

from primat.config import PRIMATConfig

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def test_cprimat_version_matches_pyproject():
    """primat-c/include/config.h's CPRIMAT_VERSION must track pyproject.toml's version.

    CLAUDE.md documents this sync as manual ("update CPRIMAT_VERSION by hand
    in the same commit") with no automated check. Parse both files directly
    (no import of a built/installed package) so this test only depends on
    the two source files staying in the same commit.
    """
    pyproject_path = os.path.join(REPO_ROOT, "pyproject.toml")
    pyproject_text = open(pyproject_path).read()
    pyproject_match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', pyproject_text)
    assert pyproject_match, "version field not found in pyproject.toml"
    pyproject_version = pyproject_match.group(1)

    config_h_path = os.path.join(REPO_ROOT, "primat-c", "include", "config.h")
    config_h_text = open(config_h_path).read()
    config_h_match = re.search(r'#define\s+CPRIMAT_VERSION\s+"([^"]+)"', config_h_text)
    assert config_h_match, "CPRIMAT_VERSION macro not found in primat-c/include/config.h"
    config_h_version = config_h_match.group(1)

    assert config_h_version == pyproject_version, (
        f"CPRIMAT_VERSION ({config_h_version!r}) in primat-c/include/config.h "
        f"is out of sync with pyproject.toml's version ({pyproject_version!r}); "
        "update both in the same commit (see CLAUDE.md)."
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
    tree = ast.parse(open(path).read(), filename=path)
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
def test_reference_run_params_match_claude_md(key, expected):
    """The param names/values CLAUDE.md quotes for the reference run must exist verbatim."""
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
    """README's MC section (FABLEADVICE F-1/F-2) quotes these DEFAULT_PARAMS
    key names verbatim; a rename that isn't mirrored in README would leave
    users following a documented option that raises an 'unknown parameter'
    warning instead of doing anything."""
    from primat.config import DEFAULT_PARAMS
    assert key in DEFAULT_PARAMS, f"{key!r} no longer a DEFAULT_PARAMS key"
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = open(readme_path).read()
    assert key in readme_text, f"{key!r} no longer documented in README.md"


def test_readme_does_not_reference_old_mc_file_key():
    """output_mc_file was hard-renamed to output_mc_file_prefix (no deprecated
    alias, author decision -- primat is not on PyPI yet); README must not
    keep referencing the old name."""
    from primat.config import DEFAULT_PARAMS
    assert "output_mc_file" not in DEFAULT_PARAMS
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = open(readme_path).read()
    # Match the old key as a whole word so "output_mc_file_prefix" (the
    # correct, current name) does not trip this assertion.
    assert not re.search(r'\boutput_mc_file\b(?!_prefix)', readme_text)


def test_streamlit_wheel_matches_pyproject_version():
    """requirements.txt's last line (the Streamlit Cloud deployment chain,
    see CLAUDE.md/wheels/README.md -- FABLEADVICE.md S-5) must point at a
    wheel file that (a) actually exists under wheels/ and (b) has the same
    version as pyproject.toml, or the public demo silently keeps serving an
    old build after a version bump."""
    pyproject_path = os.path.join(REPO_ROOT, "pyproject.toml")
    pyproject_text = open(pyproject_path).read()
    pyproject_version = re.search(r'(?m)^version\s*=\s*"([^"]+)"', pyproject_text).group(1)

    requirements_path = os.path.join(REPO_ROOT, "requirements.txt")
    lines = [l.strip() for l in open(requirements_path).read().splitlines() if l.strip()]
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
    README used to show (FABLEADVICE.md S-3) -- which primat --set rejects."""
    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = open(readme_path).read()
    assert "--set KEY=VALUE" in readme_text
    assert "--set tau_n=880.1" in readme_text
    # The old (wrong) space-separated form must not reappear.
    assert "--set tau_n 880.1" not in readme_text


def test_readme_python_only_features_list_matches_backend():
    """README's 'Python-only features' list (FABLEADVICE.md S-3) must match
    primat/backend.py's actual auto-fallback gate. After O-8 the only
    inherently-Python run_bbn feature is background= (a custom Background
    object); extra_rho, decay_era, custom_network and output_time_evolution
    are all supported on the C backend (only MC prev is additionally
    Python-only, on the run_mc side)."""
    backend_path = os.path.join(REPO_ROOT, "primat", "backend.py")
    backend_text = open(backend_path).read()
    # The feature actually gated in run_bbn()'s fallback logic -- if this
    # string disappears from backend.py, the module was refactored and
    # README's list needs re-verifying against the new code.
    assert "python_only_feature = background is not None" in backend_text

    readme_path = os.path.join(REPO_ROOT, "README.md")
    readme_text = open(readme_path).read()
    assert "background=" in readme_text
    assert "MC `prev`" in readme_text
    # extra_rho/decay_era/custom_network/output_time_evolution must NOT be
    # listed as Python-only any more -- all supported on the C backend (O-8).
    python_only_section = readme_text[readme_text.index("Python-only features"):]
    python_only_section = python_only_section[:python_only_section.index("### Using primat-c directly")]
    assert "custom_network` (GUI" not in python_only_section
    assert "output_time_evolution=True (write full time series)" not in python_only_section
    # ...they appear only in the sentence explicitly excluding them from the
    # Python-only list.
    assert "extra_rho" in python_only_section and "decay_era" in python_only_section
    assert "all four are supported on the C" in python_only_section


def test_notebooks_readme_lists_every_notebook():
    """notebooks/README.md's table must mention every notebooks/*.ipynb file
    (FABLEADVICE.md S-8) -- a new notebook silently missing from the README
    is undiscoverable from the folder's own index."""
    notebooks_dir = os.path.join(REPO_ROOT, "notebooks")
    readme_text = open(os.path.join(notebooks_dir, "README.md")).read()
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
# DEFAULT_PARAMS three-file sync (CLAUDE.md's "Keeping DEFAULT_PARAMS,
# primat_run_explanatory.py, and run_basic.ini in sync"). Until now this was
# a purely manual rule; these tests make drift a test failure.
# ---------------------------------------------------------------------------

_TEMPLATE_PY = os.path.join(REPO_ROOT, "runfiles", "primat_run_explanatory.py")
_TEMPLATE_INI = os.path.join(REPO_ROOT, "primat-c", "examples", "run_basic.ini")


def _template_py_keys():
    """Keys listed in the Python template's commented-out cfg dict."""
    text = open(_TEMPLATE_PY).read()
    # Entries look like "    # numerical_precision=1e-7,  # comment" -- key
    # immediately followed by '=' (prose comment lines never match this).
    return set(re.findall(r"^\s*#\s*([A-Za-z_][A-Za-z0-9_]*)=", text, re.M))


def _template_ini_keys():
    """Keys listed in the INI template (commented or active settings)."""
    text = open(_TEMPLATE_INI).read()
    # Entries look like "# verbose = false" or "numerical_precision = 1e-7".
    # The INI template ALSO carries prose lines that explain a boolean's
    # meaning in the same "# word = text" shape, e.g.
    # "# false = instantaneous decoupling" / "# true = Born approximation";
    # those capture the literal value word, never a parameter, so drop the
    # boolean/None literals (no DEFAULT_PARAMS key is ever named that).
    keys = set(re.findall(r"^#?\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", text, re.M))
    return keys - {"true", "false", "none"}


def test_explanatory_template_lists_every_default_param():
    from primat.config import DEFAULT_PARAMS
    found, expected = _template_py_keys(), set(DEFAULT_PARAMS)
    assert expected - found == set(), \
        f"primat_run_explanatory.py is missing keys: {sorted(expected - found)}"
    assert found - expected == set(), \
        f"primat_run_explanatory.py lists stale/unknown keys: {sorted(found - expected)}"


def test_ini_template_lists_every_default_param():
    from primat.config import DEFAULT_PARAMS
    found, expected = _template_ini_keys(), set(DEFAULT_PARAMS)
    assert expected - found == set(), \
        f"run_basic.ini is missing keys: {sorted(expected - found)}"
    assert found - expected == set(), \
        f"run_basic.ini lists stale/unknown keys: {sorted(found - expected)}"


def test_param_count_comments_match_default_params():
    """The two templates' count comments must quote len(DEFAULT_PARAMS)
    exactly. (CLAUDE.md's '(currently NN keys)' count is NOT asserted --
    CLAUDE.md is an untracked local file since Task 0; keep it updated by
    hand.)"""
    from primat.config import DEFAULT_PARAMS
    n = len(DEFAULT_PARAMS)
    assert f"All {n} DEFAULT_PARAMS keys are listed" in open(_TEMPLATE_PY).read()
    assert f"all {n} keys round-trip" in open(_TEMPLATE_INI).read()


def test_validation_reference_table_matches_reference_constants():
    """tests/README.md's 'Validation reference' tables (moved there from
    CLAUDE.md by Task 0) and tests/reference_values.py must quote the same
    numbers, so neither can be updated without the other (the regression
    tier asserts the constants against actual solves)."""
    from reference_values import (REF_SMALL_YPBBN, REF_SMALL_DOH,
                                  REF_LARGE8_YPBBN, REF_LARGE8_DOH)
    readme = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "README.md")).read()
    section = readme[readme.index("## Validation reference"):]
    yp = re.findall(r"\|\s*YP \(BBN\)\s*\|\s*([0-9.eE+-]+)\s*\|", section)
    dh = re.findall(r"\|\s*D/H\s*\|\s*([0-9.eE+-]+)\s*\|", section)
    # First table is the small network, second is large+amax=8.
    assert [float(v) for v in yp[:2]] == [REF_SMALL_YPBBN, REF_LARGE8_YPBBN]
    assert [float(v) for v in dh[:2]] == [REF_SMALL_DOH, REF_LARGE8_DOH]
