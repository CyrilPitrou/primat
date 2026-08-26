# -*- coding: utf-8 -*-
"""
primat.tools.gen_param_templates
=================================
Generator for the three user-facing parameter listings that used to require
hand-editing on every ``DEFAULT_PARAMS`` change:

- ``runfiles/primat_run_explanatory.py`` -- Python template
- ``primat-c/examples/run_basic.ini``    -- INI template for the C backend
- ``docs/parameters.md``                 -- the documentation site's reference

All three are generated from the same three sources of truth: ``DEFAULT_PARAMS``
(keys + default values), ``PARAM_GROUPS`` (section grouping/order, both in
``primat/config.py``), and ``_TEMPLATE_DESCRIPTIONS`` below (one curated
one-line description per key -- kept here, not derived from config.py's own
inline comments, because those are written for a physicist reading the
source and are often multi-line/discursive, whereas a template description
must be a single line short enough to sit next to a commented-out
``key=value``).

Usage::

    python -m primat.tools.gen_param_templates          # write all three files
    python -m primat.tools.gen_param_templates --check   # exit 1 on drift

``tests/test_docs_consistency.py`` already checks that both templates list
every ``DEFAULT_PARAMS`` key and quote the right key count; running this
module (or importing :func:`check` from it) additionally guarantees the
*content* -- not just the key set -- exactly matches what ``DEFAULT_PARAMS``/
``PARAM_GROUPS``/``_TEMPLATE_DESCRIPTIONS`` would produce, so a stale
docstring paraphrase can no longer hide from review.
"""
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from primat.config import DEFAULT_PARAMS, PARAM_GROUPS  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TEMPLATE_PY = os.path.join(REPO_ROOT, "runfiles", "primat_run_explanatory.py")
_TEMPLATE_INI = os.path.join(REPO_ROOT, "primat-c", "examples", "run_basic.ini")
_DOCS_PARAMS = os.path.join(REPO_ROOT, "docs", "parameters.md")


# ---------------------------------------------------------------------------
# One curated one-line description per DEFAULT_PARAMS key, shared verbatim by
# both templates. Keep in sync with config.py's own (more detailed) inline
# comments when the physics they describe changes; test_config.py only
# checks completeness (every key present), not wording, so this needs a
# human's attention on review, same as any other doc.
# ---------------------------------------------------------------------------
_TEMPLATE_DESCRIPTIONS = {
    "verbose": "print primat's own progress messages",
    "debug": "print extra debug messages",
    "show_progress": "print compact stderr progress indicators (HT./MT./LT./done., MC counter) when verbose=False",
    "numerical_precision": "rtol for all solve_ivp/ODE calls",
    "use_numba": "a request, not a fact: set False to run without the JIT kernels even where numba is installed; True is cleared at runtime if it is not (unused on the C backend, kept for round-trip parity)",
    "strict_params": "True: raise on an unknown param key (typo); False: warn with a \"did you mean ...?\" hint",
    "incomplete_decoupling": "True: NEVO non-instantaneous decoupling table; False: instantaneous decoupling",
    "QED_corrections": "QED corrections to the EM plasma equation of state",
    "n_electron_table": "grid points for the electron-thermo tables",
    "recompute_electron_thermo": "force recomputation of the electron-thermo cache",
    "recompute_qed_corrections": "force recomputation (and overwrite) of cache_plasma_weak/plasma/QED_*.txt",
    "spectral_distortions": "n<->p rate corrections from non-Fermi-Dirac neutrino spectra",
    "analytic_distortions": "True: analytic y_SZ/y_gray distortion instead of the NEVO spectral table",
    "y_SZ": "amplitude of the y-type (Compton/SZ-like) distortion",
    "y_gray": "amplitude of the gray-type (temperature-rescaling) distortion",
    "nevo_file": "override the 6/7-column thermo table",
    "nevo_spectral_file": "override the 86-column spectral-distortion table",
    "nevo_grid_file": "override the y-grid for nevo_spectral_file",
    "nevo_file_prefix": "base filename for the default NEVO thermo/spectral tables",
    "data_dir": "replace the entire primat/data/ tree (NEVO/, cache_plasma_weak/, nuclear/, csv/)",
    "user_nuclear_dir": "additive overlay for nuclear networks & rate tables only (primat/data/nuclear/ equivalent)",
    "external_scale_factor": "read a(T_gamma) directly from the NEVO table's x column",
    "custom_background": "path to a user-supplied background file (T, t, a columns)",
    "GN": "Newton's constant, SI [m^3 kg^-1 s^-2] (CODATA literal)",
    "alphaem": "fine-structure constant (CODATA 2018)",
    "GF": "Fermi constant [MeV^-2] (PDG 2020)",
    "mZ": "Z boson mass [MeV] (PDG 2020)",
    "me": "electron mass [MeV] (CODATA 2018)",
    "mn": "neutron mass [MeV] (PDG/CODATA 2018)",
    "mp": "proton mass [MeV] (CODATA 2018)",
    "T0CMB": "CMB photon temperature today [K] (Fixsen 2009)",
    "gA": "nucleon axial coupling (PDG 2018)",
    "Vud": "CKM matrix element V_ud (PDG 2018); drops out when tau_n_normalization=True",
    "kappa_p": "proton anomalous magnetic moment (CODATA 2018)",
    "kappa_n": "neutron anomalous magnetic moment (CODATA 2018)",
    "radproton": "proton charge radius [cm] (CODATA 2018)",
    "ma": "unified atomic mass unit [MeV] (CODATA 2010)",
    "He4Overma": "M(He4) / u (AME2020)",
    "HOverma": "M(H) / u (AME2016)",
    "Neff_SM": "Standard-Model prediction for Neff; only used to normalise the EDE era",
    "T_start_cosmo_MeV": "starting temperature [MeV]",
    "T_end_MeV": "end temperature for nuclear integration [MeV]",
    "sampling_temperature_per_decade": "points per decade of T for the background a(T)/t(T) grid",
    "radiative_corrections": "Coulomb + T=0 resummed radiative corrections (CCR); False: Born approximation",
    "finite_mass_corrections": "Fokker-Planck finite-nucleon-mass correction",
    "thermal_corrections": "finite-temperature radiative corrections (CCRTh)",
    "cache_dir": "writable dir for ALL regenerable caches (weak-rate + plasma); unset = <data_dir>/cache_plasma_weak/. Set on read-only installs: caches are written to <cache_dir>/{weak,plasma}/ and read from there first, falling back to the shipped caches (overlay). Not part of any fingerprint.",
    "weak_rate_cache": "if False, never load the non-thermal weak-rate cache (the thermal one is loaded whenever its file exists)",
    "save_nTOp": "save computed n<->p rates to cache_plasma_weak/weak/ (or the cache_dir redirect)",
    "sampling_nTOp_per_decade": "points per decade of T in the n<->p rate grid",
    "save_nTOp_thermal": "save computed thermal n<->p rates to cache_plasma_weak/weak/ (or the cache_dir redirect)",
    "sampling_nTOp_thermal_per_decade": "points per decade of T for the thermal-correction table",
    "tau_n_normalization": "normalise weak rates using the neutron lifetime tau_n",
    "tau_n": "neutron lifetime [s]",
    "std_tau_n": "1-sigma uncertainty on tau_n [s] (used for MC sampling)",
    "vegas_n_eval": "vegas: evaluations per iteration (thermal-correction integral; unused on the C backend, deterministic quadrature only)",
    "vegas_n_itn": "vegas: number of iterations (unused on the C backend)",
    "epsrel_thermal": "dblquad fallback relative tolerance",
    "output_time_evolution": "write the unified time-evolution TSV (see evolution.py's module docstring for the schema)",
    "output_rates_time_evolution": "append per-reaction forward-rate columns (<reaction>_frwrd) to the time-evolution TSV; one per reaction in the active LT network (~12 small, ~429 full large)",
    "output_n_points": "number of points in the time-evolution TSV",
    "output_file": "path for output_time_evolution",
    "output_final_result": "write a two-column (nuclide, Y) final-abundances file",
    "output_final_file": "path for output_final_result",
    "output_background_evolution": "write the cosmological background's own time-evolution TSV",
    "output_background_file": "path for output_background_evolution",
    "output_mc_samples": "write every MC sample to <output_mc_file_prefix>_samples.tsv (one column per quantity)",
    "output_mc_covariance": "write the (n_q,n_q) sample covariance matrix (ddof=1) to <output_mc_file_prefix>_covariance.tsv",
    "output_mc_correlation": "write the matching correlation matrix to <output_mc_file_prefix>_correlation.tsv",
    "output_mc_file_prefix": "filename stem for the three MC output files above",
    "rate_grid_npts": "points in the master T9 grid used to resample every rate table",
    "rate_grid_T9_min": "minimum T9 [GK] of the master rate grid",
    "rate_grid_T9_max": "maximum T9 [GK] of the master rate grid",
    "network": '"small" / "small_parthenope" / "large" / custom network filename',
    "amax": "filter any network to reactions with A <= amax",
    "atol_LT": "solve_ivp absolute tolerance for the LT era of every network (the name is historical)",
    "mc_rate_rescale_cap": "clamp the MC rate variation factor to [1/cap, cap]; None disables the cap",
    "nuclear_qed_corrections": "QED correction to select radiative-capture rates (Pitrou & Pospelov 2020)",
    "Omegabh2": "baryon density Omega_b h^2 (Planck 2018 default)",
    "Omegach2": "cold dark matter density Omega_c h^2 (Planck 2018)",
    "h": "reduced Hubble constant h = H_0 / (100 km/s/Mpc) (Planck 2018)",
    "DeltaNeff": "extra relativistic species beyond SM neutrinos",
    "munuOverTnu": "reduced neutrino chemical potential xi = mu/T (common default for all 3 flavours)",
    "munuOverTnu_e": "per-flavour xi_e of nu_e; None = inherit munuOverTnu (only xi_e shifts the n<->p weak rates)",
    "munuOverTnu_mu": "per-flavour xi_mu of nu_mu; None = inherit munuOverTnu (gravitates only, via Neff)",
    "munuOverTnu_tau": "per-flavour xi_tau of nu_tau; None = inherit munuOverTnu (gravitates only, via Neff)",
    "decay_reverse_rates": "compute detailed-balance reverse rates for radioactive decays",
    "decay_era": 'run a 4th "Decay Time" era after LT (network="large" only)',
    "t_decay_end": "DT era duration [s] (default: 1 Gyr)",
    "decay_n_points": "log-spaced output points in the DT era",
    "output_decay_evolution": "write a TSV of the DT-era abundance time evolution",
    "output_decay_file": "path for output_decay_evolution",
    "fEDE": "EDE fraction at peak; 0 = disabled",
    "zcEDE": "redshift of EDE peak",
    "wnEDE": "EDE equation-of-state parameter",
}


def _check_descriptions_complete():
    missing = set(DEFAULT_PARAMS) - set(_TEMPLATE_DESCRIPTIONS)
    stale = set(_TEMPLATE_DESCRIPTIONS) - set(DEFAULT_PARAMS)
    if missing or stale:
        raise AssertionError(
            f"_TEMPLATE_DESCRIPTIONS out of sync with DEFAULT_PARAMS: "
            f"missing={sorted(missing)} stale={sorted(stale)}")


def _py_literal(value) -> str:
    """Render a DEFAULT_PARAMS value as a Python literal for the ``.py``
    template (double-quoted strings, to match the surrounding file style)."""
    if value is None or isinstance(value, bool):
        return repr(value)
    if isinstance(value, str):
        return f'"{value}"'
    return repr(value)


def _ini_literal(value):
    """Render a DEFAULT_PARAMS value for the ``.ini`` template, or ``None``
    if the value has no direct INI representation (the ``None`` sentinel
    itself -- ``cpr_ini_load`` has no null literal), in which case the
    caller emits an ``<unset>`` placeholder line instead."""
    if value is None:
        return None
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    return repr(value)


def generate_run_explanatory() -> str:
    """Return the full text of ``runfiles/primat_run_explanatory.py``."""
    _check_descriptions_complete()
    n = len(DEFAULT_PARAMS)
    lines = [
        '# -*- coding: utf-8 -*-',
        '"""',
        'primat_run_explanatory.py',
        '==========================',
        'Minimal, heavily-commented template for a standalone BBN run. Copy this file',
        'and uncomment/edit whichever options you need; every option shown below is at',
        'its default value, so running this file unmodified reproduces the standard',
        'run (see tests/README.md\'s "Validation reference (authoritative copy)"',
        'for the expected YPBBN/D-H values and the tolerance that applies).',
        '',
        'Run from the repo root so that the shipped ``data/`` data resolve correctly:',
        '',
        '    python runfiles/primat_run_explanatory.py',
        '',
        'Most runs change only four of the keys below: ``network`` (which reaction',
        'set to integrate), ``Omegabh2`` (the baryon density), ``DeltaNeff`` (extra',
        'relativistic species) and ``tau_n`` (the neutron lifetime). The rest are',
        'listed so that nothing is hidden, not because a first run needs them.',
        'docs/parameters.md is the same list as a reference table, and',
        'docs/tutorials/first-run.md walks through what each printed number means.',
        '',
        'Generated by primat/tools/gen_param_templates.py -- do not hand-edit the cfg',
        'block below; regenerate with `python -m primat.tools.gen_param_templates`.',
        '"""',
        'import os',
        'import sys',
        '',
        'sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))',
        '',
        'from primat.backend import run_bbn  # , run_mc, dump_mc_samples (see bottom of file)',
        '',
        'cfg = dict(',
        '    # Every key below is shown at its DEFAULT_PARAMS default (see',
        '    # primat/config.py for the authoritative, more detailed comments this',
        '    # file summarises); uncomment and edit whichever you need to override.',
        f'    # All {n} DEFAULT_PARAMS keys are listed, grouped exactly as in config.py.',
    ]
    for group, keys in PARAM_GROUPS.items():
        lines.append('')
        lines.append(f'    # ---- {group} ----')
        for key in keys:
            value = _py_literal(DEFAULT_PARAMS[key])
            desc = _TEMPLATE_DESCRIPTIONS[key]
            lines.append(f'    # {key}={value},  # {desc}')
    lines.append(')')
    lines.extend([
        '# force_backend: None/"auto" (default: C extension if built, else pure Python),',
        '# "c", or "python" -- see primat/backend.py\'s module docstring for exactly',
        '# which features (extra_rho/custom_network/background=, output_time_evolution)',
        '# always fall back to "python" regardless of this setting.',
        'result = run_bbn(cfg, force_backend="auto")',
        'print("Neff  =", result.get("Neff"))',
        'print("YPBBN =", result["YPBBN"])',
        'print("D/H   =", result["DoH"])',
        '',
        '# Monte-Carlo nuclear-rate/tau_n uncertainty propagation (uncomment to run):',
        '# the same dispatch story as run_bbn -- C backend when available, else',
        '# pure-Python joblib (see primat/backend.py\'s module docstring for the RNG',
        '# caveat: C and Python samples are statistically, not bit-for-bit, equal).',
        '#',
        '# from primat.backend import (run_mc, dump_mc_samples,',
        '#                              dump_mc_covariance, dump_mc_correlation)',
        "# mc = run_mc(50, ['YPBBN', 'DoH'], params=cfg, force_backend=\"auto\")",
        '# print("YPBBN =", mc[\'YPBBN\'].mean, "+/-", mc[\'YPBBN\'].std)',
        '# # Joint uncertainty (needed for a multi-abundance likelihood): the full',
        '# # covariance/correlation matrices, plus scalar access by name.',
        "# print(\"corr(YPBBN, DoH) =\", mc.corr('YPBBN', 'DoH'))",
        '# with open("results/output_mc_samples.tsv", "w") as f:',
        '#     f.write(dump_mc_samples(mc))',
        '# with open("results/output_mc_covariance.tsv", "w") as f:',
        '#     f.write(dump_mc_covariance(mc))',
        '# with open("results/output_mc_correlation.tsv", "w") as f:',
        '#     f.write(dump_mc_correlation(mc))',
        '',
    ])
    return "\n".join(lines)


def generate_run_basic_ini() -> str:
    """Return the full text of ``primat-c/examples/run_basic.ini``."""
    _check_descriptions_complete()
    n = len(DEFAULT_PARAMS)
    lines = [
        '# run_basic.ini -- minimal, heavily-commented template for a standalone',
        '# primat-c run (KEY=VALUE syntax, one setting per line -- trailing',
        '# same-line comments are NOT supported by cpr_ini_load, see src/ini.c).',
        '# Copy this file and uncomment/edit whichever options you need; every',
        '# option shown below is at its default value, so running this file',
        '# unmodified reproduces the standard run (see ../../tests/README.md\'s',
        '# "Validation reference (authoritative copy)" for the expected YPBBN/D-H',
        '# values and the tolerance that applies). The shipped data/ tree is found',
        '# automatically next to the binary; pass --data_dir PATH to point elsewhere:',
        '#',
        '#   cd primat-c && make && ./build/primat-c --ini examples/run_basic.ini',
        '#',
        '# Two shorter starting points live beside this file: run_small.ini (the',
        '# default 12-reaction network) and run_large_amax8.ini (68 reactions).',
        '#',
        '# Every key below mirrors a DEFAULT_PARAMS entry in ../../primat/config.py',
        f'# (see that file for more detailed physics comments); all {n} keys round-trip',
        "# through cpr_config_set_by_name (primat-c/src/config.c's FIELD_TABLE), so",
        '# any of them may be uncommented here. A value of "<unset>" marks a key whose',
        '# Python default is None (INI has no null literal); replace it with a real',
        '# value to use it, or leave the line commented to keep the field unset.',
        '#',
        '# Generated by primat/tools/gen_param_templates.py -- do not hand-edit the',
        '# key list below; regenerate with `python -m primat.tools.gen_param_templates`.',
    ]
    for group, keys in PARAM_GROUPS.items():
        lines.append('')
        lines.append(f'# ---- {group} ----')
        for key in keys:
            desc = _TEMPLATE_DESCRIPTIONS[key]
            ini_value = _ini_literal(DEFAULT_PARAMS[key])
            lines.append(f'# {key}: {desc}')
            if ini_value is None:
                lines.append(f'# {key} = <unset>')
            else:
                lines.append(f'# {key} = {ini_value}')
    lines.append('')
    return "\n".join(lines)


# A bare ``evolution.py`` or ``_samples.tsv`` in a description reads as a
# hostname to MyST's linkify extension, which turns it into a dead
# ``http://evolution.py`` the docs link checker then fails on. Backticks stop
# that, and are what such a token should be in prose anyway.
_FILENAME_TOKEN = re.compile(r"(?<![`\w/])([A-Za-z][\w-]*\.(?:py|ini|txt|csv|tsv|md))\b")


def _md_cell(text: str) -> str:
    """Escape a description or default for a Markdown table cell.

    ``<`` and ``>`` would be swallowed as an HTML tag by MyST (several
    descriptions contain ``<data_dir>``, ``<reaction>_frwrd``, ``n<->p``), and
    a bare ``|`` would end the cell.
    """
    text = _FILENAME_TOKEN.sub(r"`\1`", text)
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;").replace("|", "\\|"))


def generate_parameters_page() -> str:
    """Return the full text of ``docs/parameters.md``."""
    _check_descriptions_complete()
    n = len(DEFAULT_PARAMS)
    lines = [
        "<!-- Generated by primat/tools/gen_param_templates.py -- do not hand-edit;",
        "     regenerate with `python -m primat.tools.gen_param_templates`. -->",
        "",
        "# Parameters",
        "",
        f"Every one of primat's {n} run parameters, grouped as `primat/config.py`",
        "groups them. All are optional: a run overrides only what it names and takes",
        "the rest from the defaults below.",
        "",
        "The same key works everywhere a parameter can be set:",
        "",
        "```python",
        'run_bbn({"Omegabh2": 0.02242, "network": "large", "amax": 8})',
        "```",
        "",
        "```bash",
        "primat --Omegabh2 0.02242 --network large --amax 8   # named flag, if it has one",
        "primat --set T_end_MeV=1e-4                          # any key at all",
        "```",
        "",
        "```ini",
        "; an --ini file for the C backend or the standalone C CLI",
        "Omegabh2 = 0.02242",
        "network = large",
        "```",
        "",
        "This table is generated from `DEFAULT_PARAMS`, so it cannot drift from the",
        "code. `primat --list-params` prints the same list in the terminal, and",
        "`runfiles/primat_run_explanatory.py` and `primat-c/examples/run_basic.ini`",
        "are the same list again as copy-and-edit templates.",
        "",
        "Per-reaction rate knobs (`p_<reaction>`, `delta_<reaction>`) are not listed",
        "here: they depend on the chosen network, so there is no fixed set of them.",
        "`primat --list-reactions` prints the names a given `--network`/`--amax`",
        "accepts — see {doc}`howto/rate-variation-mc`.",
        "",
        "A `None` default means the feature is off or the path is unset.",
    ]
    for group, keys in PARAM_GROUPS.items():
        lines.append("")
        lines.append(f"## {group}")
        lines.append("")
        lines.append("| Parameter | Default | What it does |")
        lines.append("|---|---|---|")
        for key in keys:
            desc = _md_cell(_TEMPLATE_DESCRIPTIONS[key])
            default = _md_cell(_py_literal(DEFAULT_PARAMS[key]))
            lines.append(f"| `{key}` | `{default}` | {desc} |")
    lines.append("")
    return "\n".join(lines)


# The generated files, in the order the module writes them.
_OUTPUTS = ((_TEMPLATE_PY, generate_run_explanatory),
            (_TEMPLATE_INI, generate_run_basic_ini),
            (_DOCS_PARAMS, generate_parameters_page))


def check() -> bool:
    """Compare the committed generated files against freshly generated text.

    Returns True (and prints nothing) if all three match; otherwise prints
    which file(s) are stale and returns False. Used by both the ``--check``
    CLI flag and ``tests/test_docs_consistency.py``.
    """
    ok = True
    for path, generate in _OUTPUTS:
        with open(path, encoding="utf-8") as f:
            committed = f.read()
        fresh = generate()
        if committed != fresh:
            print(f"STALE: {path} does not match generated output "
                  f"-- run `python -m primat.tools.gen_param_templates`.")
            ok = False
    return ok


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                         help="check for drift without writing files; exit 1 if stale")
    args = parser.parse_args(argv)

    if args.check:
        sys.exit(0 if check() else 1)

    for path, generate in _OUTPUTS:
        with open(path, "w", encoding="utf-8") as f:
            f.write(generate())
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
