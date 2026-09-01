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
CI and cannot be asserted against. Anything that must be enforced lives in a
tracked file instead: the "Validation reference" numbers are in
`tests/README.md` and `tests/reference_values.py`, where the two
`*_matches_reference_constants` tests below pin them.

No test in this module runs a solve -- most are static file reads, the rest
load a network or collect the suite -- so the whole file stays in the fast
(`-m "not slow"`) lane.
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


def _pyproject_version():
    """pyproject.toml's version -- the single source of truth for the release."""
    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"',
                      _read_text(os.path.join(REPO_ROOT, "pyproject.toml")))
    assert match, "version field not found in pyproject.toml"
    return match.group(1)


def test_citation_cff_version_matches_pyproject():
    """CITATION.cff's version is what GitHub's "Cite this repository" and Zenodo
    report, so a bump that misses it hands out the wrong release number.

    The release date is checked against CHANGELOG.md's heading for the same
    version rather than asserted literally: the two are written in different
    commits and drifted once, leaving a 0.3.1 date on a 0.3.2 record.
    """
    cff_text = _read_text(os.path.join(REPO_ROOT, "CITATION.cff"))
    version = _pyproject_version()

    cff_version = re.search(r'(?m)^version:\s*(\S+)\s*$', cff_text)
    assert cff_version, "version field not found in CITATION.cff"
    assert cff_version.group(1).strip('"\'') == version, (
        f"CITATION.cff's version ({cff_version.group(1)}) is out of sync with "
        f"pyproject.toml's ({version!r}); update it in the same commit "
        "(see PyPiGuide.md, Step 1)."
    )

    cff_date = re.search(r'(?m)^date-released:\s*"?([0-9]{4}-[0-9]{2}-[0-9]{2})"?',
                         cff_text)
    assert cff_date, "date-released field not found in CITATION.cff"

    changelog = _read_text(os.path.join(REPO_ROOT, "CHANGELOG.md"))
    heading = re.search(r'(?m)^## \[%s\] - ([0-9]{4}-[0-9]{2}-[0-9]{2})\s*$'
                        % re.escape(version), changelog)
    if heading is None:
        # The version is still unreleased (no dated CHANGELOG heading yet);
        # nothing to compare the date against.
        return
    assert cff_date.group(1) == heading.group(1), (
        f"CITATION.cff says {version} was released on {cff_date.group(1)}, "
        f"CHANGELOG.md says {heading.group(1)}."
    )


def test_manual_declares_the_current_version():
    """manual/ declares the version it documents in four tracked places plus its
    own filenames, none of which a version bump touches automatically.

    The .tex title page and the intro paragraph are the ones a reader sees; the
    two filenames and manual/README.md's four references to them are what a bump
    has to rename. A stale manual claims to describe a release it predates.
    """
    version = _pyproject_version()
    manual_dir = os.path.join(REPO_ROOT, "manual")

    tex_name = "primat_documentation_v%s.tex" % version
    pdf_name = "primat_documentation_v%s.pdf" % version
    for name in (tex_name, pdf_name):
        assert os.path.isfile(os.path.join(manual_dir, name)), (
            f"manual/{name} does not exist -- rename the manual's .tex/.pdf to "
            f"the current version {version!r} and update manual/README.md "
            "(see PyPiGuide.md, Step 1)."
        )

    tex_text = _read_text(os.path.join(manual_dir, tex_name))
    assert f"Documentation for primat version {version}" in tex_text, (
        f"manual/{tex_name}'s title page does not declare version {version!r}."
    )

    readme_text = _read_text(os.path.join(manual_dir, "README.md"))
    assert f"**primat version {version}**" in readme_text, (
        f"manual/README.md does not declare version {version!r}."
    )
    assert tex_name in readme_text and pdf_name in readme_text, (
        f"manual/README.md does not reference {tex_name} and {pdf_name}."
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
    """Every tracked prose copy of the parameter count must quote
    len(DEFAULT_PARAMS) exactly -- see this module's docstring for why only
    tracked files are asserted."""
    from primat.config import DEFAULT_PARAMS
    n = len(DEFAULT_PARAMS)
    assert f"All {n} DEFAULT_PARAMS keys are listed" in _read_text(_TEMPLATE_PY)
    assert f"all {n} keys round-trip" in _read_text(_TEMPLATE_INI)
    # docs/index.md's is the one a reader meets first, and unlike the two
    # templates it is not regenerated, so nothing else would catch it.
    index_md = os.path.join(REPO_ROOT, "docs", "index.md")
    assert f"lists all {n} keys" in _read_text(index_md)


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


# --- The README's two console blocks -------------------------------------
#
# The README prints two blocks of real console output that nothing used to
# compare against: the CLI summary under "Command-line interface" and the
# --mc block under "Computing the uncertainty". Both are pinned below.


def _block_with(text, needle, offset=0):
    """The fenced code block containing ``needle``, or ``offset`` blocks later.

    ``offset=1`` picks up a command block's output block, which the README
    prints as the next fence down.
    """
    blocks = re.findall(r"^```[^\n]*\n(.*?)^```", text, re.S | re.M)
    for i, body in enumerate(blocks):
        if needle in body:
            return blocks[i + offset]
    raise AssertionError(f"no fenced block contains {needle!r}")


def _summary_labels_and_values(block):
    """``{label: value}`` for every ``LABEL = VALUE`` line of a CLI summary.

    Reads the block as printed, tolerating the leading ``# `` the README puts
    on the --mc block's lines (it sits inside a ``bash`` fence). The label must
    start the line, which is what leaves out the centred banner.
    """
    out = {}
    for line in block.splitlines():
        m = re.match(r"^(?:# )?([A-Za-z0-9][A-Za-z0-9 ()/]*?)\s+=\s+(\S+)",
                     line)
        if m:
            out[m.group(1).strip()] = m.group(2)
    return out


def _cli_summary_labels():
    """The observable labels ``primat`` prints, in the order it prints them."""
    src = _read_text(os.path.join(REPO_ROOT, "primat", "cli.py"))
    return re.findall(r'print\(f"([A-Za-z0-9 ()/]+?)\s*= \{results\[', src)


def test_readme_cli_summary_block_matches_the_cli_and_the_constants():
    """README's `--network large --amax 8` output block is a real run.

    GOAL: the block is copied console output, so it goes stale silently. Its
    labels are checked against the ones primat/cli.py actually prints, and its
    two published observables against tests/reference_values.py -- which
    test_regression.py's reference tier pins to a live solve, so a physics
    change moves the constants and fails this too. Static on purpose: this
    file stays in the fast lane (see the module docstring).
    """
    from reference_values import (REF_LARGE8_DOH, REF_LARGE8_YPBBN,
                                  ROUTINE_RUN_DOH_ABS_TOL,
                                  ROUTINE_RUN_YPBBN_ABS_TOL)

    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    block = _block_with(readme, "primat --Omegabh2 0.02242 --network large "
                                "--amax 8", offset=1)
    documented = _summary_labels_and_values(block)
    printed = _cli_summary_labels()
    assert printed, ("could not read any observable label out of "
                     "primat/cli.py -- its print idiom changed, so this test "
                     "is no longer reading the contract it thinks it is")

    # Every label the block shows must be one the CLI prints, in the CLI's
    # order. The block omits "CNO (mass)", which this network does not report.
    assert set(documented) <= set(printed), (
        f"README's CLI block shows labels primat never prints: "
        f"{sorted(set(documented) - set(printed))}")
    assert [lb for lb in printed if lb in documented] == list(documented), (
        "README's CLI block lists the observables in a different order than "
        "primat/cli.py prints them")
    for required in ("Neff", "YP (BBN)", "D/H", "Li7/H"):
        assert required in documented, f"README's CLI block lost {required!r}"

    # The block runs `large, amax=8` at the default precision, so it is held
    # to the routine bound, not the reference run's tighter one.
    assert abs(float(documented["YP (BBN)"]) - REF_LARGE8_YPBBN) < \
        ROUTINE_RUN_YPBBN_ABS_TOL
    assert abs(float(documented["D/H"]) - REF_LARGE8_DOH) < \
        ROUTINE_RUN_DOH_ABS_TOL


def test_readme_mc_block_matches_the_cli_and_the_constants():
    """README's `--mc 300` block quotes real central values and a real matrix.

    GOAL: same reason as the CLI block above, plus one claim of its own -- the
    README states that the value before each `+/-` is the unperturbed central
    solve and so does not depend on N or the seed. That makes those two
    numbers the default-run constants, and pins them exactly; the sigmas and
    the correlations are finite-sample and are checked for shape, not digits.
    """
    from reference_values import DOH_REFERENCE, YPBBN_REFERENCE

    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    block = _block_with(readme, "primat --Omegabh2 0.02242 --mc 300")
    documented = _summary_labels_and_values(block)

    # Exact strings: both sides are tracked files, so a mismatch means one was
    # edited without the other (the idiom of the quick-start test above).
    assert documented["YP (BBN)"] == f"{YPBBN_REFERENCE:.8f}"
    assert documented["D/H"] == f"{DOH_REFERENCE:.7e}"
    assert "+/-" in block, "README's MC block lost the +/- sigma column"

    # The correlation matrix: the four main products the CLI pairs, in the
    # CLI's order, with a unit diagonal and symmetric off-diagonals.
    from primat.cli import _MC_MAIN_PRODUCTS
    rows = re.findall(r"^#\s*(" + "|".join(_MC_MAIN_PRODUCTS) + r")\s+"
                      r"([-\d. ]+)$", block, re.M)
    assert [name for name, _ in rows] == list(_MC_MAIN_PRODUCTS), (
        "README's correlation matrix does not list primat --mc's four main "
        f"products in order; found {[n for n, _ in rows]}")
    matrix = [[float(x) for x in cells.split()] for _, cells in rows]
    n = len(_MC_MAIN_PRODUCTS)
    for i in range(n):
        assert matrix[i][i] == 1.0, "correlation diagonal is not 1"
        for j in range(n):
            assert matrix[i][j] == matrix[j][i], "correlation is not symmetric"
            assert -1.0 <= matrix[i][j] <= 1.0, "correlation outside [-1, 1]"


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


@pytest.mark.parametrize("page", ["docs/index.md", "README.md",
                                  "docs/tutorials/first-run.md"])
def test_quick_start_numbers_match_the_reference_constants(page):
    """Every page quoting the quick start must use the tracked reference values.

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
    m_py = re.search(
        r"pure-Python backend gives `([0-9.]+)` / `([0-9.e+-]+)`", text)

    # Exact strings, not approx: both sides are tracked files, so a mismatch
    # means one was edited without the other.
    assert yp_c == f"{YPBBN_REFERENCE:.8f}"
    assert doh_c == f"{DOH_REFERENCE:.7e}"
    # The tutorial quotes only the default backend's pair; the two landing
    # pages quote both, so that they can stand alone.
    if m_py is not None:
        yp_py, doh_py = m_py.groups()
        assert yp_py == f"{PY_YPBBN_REFERENCE:.8f}"
        assert doh_py == f"{PY_DOH_REFERENCE:.7e}"
    else:
        assert page.endswith("first-run.md"), f"{page} lost its Python-backend pair"


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


# --- The methods paper's citation, and two counts nothing used to pin ------

# Every place that prints the Physics Reports reference for a human to copy.
# The manual's own front matter and README use the "volume (year) first-page"
# short form instead, so they are not in this list.
_CITATION_PAGES = ("README.md", "docs/index.md", "docs/citing.md",
                   "docs/physics.md", "docs/extending.md")


def _citation_cff_reference():
    """(volume, start, end) of CITATION.cff's `preferred-citation` entry."""
    text = _read_text(os.path.join(REPO_ROOT, "CITATION.cff"))
    block = text.split("preferred-citation:", 1)[1]
    grab = lambda k: re.search(rf'^\s+{k}:\s*"?(\d+)"?\s*$', block, re.M).group(1)
    return grab("volume"), grab("start"), grab("end")


@pytest.mark.parametrize("page", _CITATION_PAGES)
def test_paper_citation_matches_citation_cff(page):
    """Every human-readable citation of the methods paper gives the volume and
    pages CITATION.cff gives.

    GOAL: the one reference this project asks a reader to copy must resolve.
    These pages once printed "Physics Reports 04 (2018) 005" -- the DOI suffix
    read as a volume and a page -- next to a BibTeX block saying 754 (2018)
    1-66.
    """
    volume, start, end = _citation_cff_reference()
    # Collapsed to one line first: docs/extending.md wraps its citation
    # mid-reference, so a line-by-line scan would see only half of it.
    text = " ".join(_read_text(os.path.join(REPO_ROOT, *page.split("/"))).split())
    cites = re.findall(r"Physics Reports[^.]{0,80}?\(2018\)[^.]{0,20}", text)
    assert cites, f"{page} no longer carries a Physics Reports citation"
    for cite in cites:
        assert volume in cite, f"{page}: citation lacks volume {volume}: {cite}"
        assert re.search(rf"{start}[-\u2013]{end}", cite), (
            f"{page}: citation lacks pages {start}-{end}: {cite}")


def test_no_document_prints_the_doi_suffix_as_a_volume():
    """"04 (2018) 005" is the DOI's month and article number, not a citation.

    GOAL: guard the whole prose tree, not only the pages listed above, against
    the mis-rendering `test_paper_citation_matches_citation_cff` fixed.
    """
    bad = re.compile(r"\b04\b[^\n]{0,4}\(2018\)[^\n]{0,4}\b005\b")
    # CHANGELOG.md quotes the mis-rendering in the entry that records fixing
    # it, which is what a changelog is for.
    paths = [os.path.join(REPO_ROOT, f) for f in os.listdir(REPO_ROOT)
             if f.endswith(".md") and f != "CHANGELOG.md"]
    for sub in ("docs", "manual", "primat"):
        for dirpath, dirnames, filenames in os.walk(os.path.join(REPO_ROOT, sub)):
            dirnames[:] = [d for d in dirnames
                           if d not in ("_build", "__pycache__", "data")]
            paths += [os.path.join(dirpath, n) for n in filenames
                      if n.endswith((".md", ".py", ".tex"))]
    offenders = [os.path.relpath(p, REPO_ROOT) for p in paths
                 if bad.search(_read_text(p))]
    assert not offenders, (
        f"these files print the DOI suffix as a volume/page: {sorted(set(offenders))}")


def test_documented_rate_column_counts_match_the_networks():
    """README and the C header quote the real `<reaction>_frwrd` column count.

    GOAL: the block has one column per LT reaction *except* n<->p, which has
    no rate table, so the count is `n_reac - 1` -- 12 / 67 / 428, not the LT
    network sizes 13 / 68 / 429 that both documents used to quote.
    """
    from primat.network_data import load_network

    counts = {}
    for label, params in (("small", {"network": "small"}),
                          ("amax8", {"network": "large", "amax": 8}),
                          ("large", {"network": "large"})):
        counts[label] = load_network(PRIMATConfig(dict(params))).n_reac - 1
    assert counts == {"small": 12, "amax8": 67, "large": 428}, counts

    readme = _read_text(os.path.join(REPO_ROOT, "README.md"))
    header = _read_text(os.path.join(
        REPO_ROOT, "primat-c", "include", "nuclear_network.h"))
    for text, where in ((readme, "README.md"),
                        (header, "primat-c/include/nuclear_network.h")):
        # Both documents state the exception, then the three counts.
        assert "no rate table" in text, f"{where}: no column-count sentence found"
        window = text[text.index("no rate table"):][:220]
        quoted = [int(n) for n in re.findall(r"\b\d{2,3}\b", window)[:3]]
        assert quoted == [counts["small"], counts["amax8"], counts["large"]], (
            f"{where} quotes {quoted} as the column counts, "
            f"but they are {[counts['small'], counts['amax8'], counts['large']]}")


def test_tests_readme_test_counts_match_the_collected_suite():
    """tests/README.md's two totals are what pytest actually collects.

    GOAL: those totals are the only numbers in that file nothing checked, and
    they drifted three tests behind the suite. One collection reports both --
    "<fast>/<total> tests collected (<slow> deselected)".
    """
    import subprocess
    import sys

    out = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/", "--collect-only", "-q",
         "-m", "not slow"],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=600).stdout
    m = re.search(r"(\d+)/(\d+) tests collected \((\d+) deselected\)", out)
    assert m, f"could not parse the collection summary from:\n{out[-2000:]}"
    total, deselected = m.group(2), m.group(3)

    readme = _read_text(os.path.join(REPO_ROOT, "tests", "README.md"))
    assert f"everything: {total} tests" in readme, (
        f"tests/README.md's total is stale; the suite collects {total} tests")
    assert f"marker excludes {deselected} tests" in readme, (
        f"tests/README.md's slow count is stale; `slow` excludes {deselected}")
