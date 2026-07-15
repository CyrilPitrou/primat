# -*- coding: utf-8 -*-
"""
primat.cli
============

Command-line entry point for primat.

This wraps the same "build a params dict and call PRIMAT" pattern used by
``runfiles/primat_run.py``, exposing the handful of options most users
need (baryon density, extra relativistic species, network choice) so a
``pip install``-ed user can get BBN abundances without writing any Python::

    primat --Omegabh2 0.02242 --network large --amax 8

The output filenames are exposed as named flags as well
(``--output_file``, ``--output_final_file``, ``--output_background_file``,
``--output_mc_file_prefix``) so they show up in ``primat --help`` alongside
the other basic options. Monte-Carlo output is controlled by three standard
config flags (``--output_mc_samples`` / ``--output_mc_covariance`` /
``--output_mc_correlation``); each writes ``<prefix>_samples.tsv`` /
``<prefix>_covariance.tsv`` / ``<prefix>_correlation.tsv`` respectively, with
the stem taken from ``--output_mc_file_prefix``.

Anything not exposed as a named flag here can still be set without writing a
script, via the (intentionally undocumented in ``--help``, to keep the
printed help short) ``--set KEY=VALUE`` escape hatch, repeatable for any
``PRIMATConfig`` key (including ``p_<reaction>``/``delta_<reaction>``
rate-variation keys), e.g.::

    primat --set T_end_MeV=1e-4 --set decay_era=True --set network=large

Values are parsed with ``ast.literal_eval`` (so ``True``/``False``/``None``,
numbers, and quoted strings all work); anything that fails to parse as a
Python literal is kept as a plain string (e.g. ``--set network=large``).
"""
import argparse
import ast
import json
import os
import sys
import time

from . import PRIMAT, __version__
from .credits import cli_credits_text
from .backend import (HAS_C_BACKEND, dump_mc_correlation, dump_mc_covariance,
                      dump_mc_samples, run_bbn, run_mc)
from .cache_utils import clear_weak_cache, list_weak_cache_files, weak_cache_dir
from .config import (DEFAULT_PARAMS, PARAM_GROUPS, PRIMATConfig,
                     _default_params_comments, _rates_overlay_notice)


def _parse_set_value(raw: str):
    """Parse the value half of a ``--set KEY=VALUE`` CLI argument.

    Tries ``ast.literal_eval`` first, so numeric, boolean, ``None``, and
    quoted-string values are converted to the right Python type (e.g.
    ``"True"`` -> ``True``, ``"1e-4"`` -> ``1e-4``). Falls back to the raw
    string unchanged when it is not a valid Python literal (e.g. an
    unquoted network name like ``large``), since ``PRIMATConfig`` string
    parameters (``network``, ``custom_background``, ...) are passed this way.

    Example
    -------
        >>> _parse_set_value("1e-4")
        0.0001
        >>> _parse_set_value("large")
        'large'
    """
    try:
        return ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return raw


# The four main BBN products whose joint MC uncertainty the CLI prints as
# aligned correlation/covariance matrices.
_MC_MAIN_PRODUCTS = ("YPBBN", "DoH", "He3oHe4", "Li7oH")


def _print_mc_matrices(mc):
    """Print the 4x4 correlation and covariance matrices of the four main BBN
    products (``YPBBN``, ``DoH``, ``He3oHe4``, ``Li7oH``) from an
    :class:`primat.main.MCResult`, aligned and labelled.

    The correlation matrix is dimensionless (printed to 3 decimals); the
    covariance matrix uses ``%.3e``.  Together they give the *joint* nuclear-
    rate/tau_n uncertainty a user needs to build a likelihood over several
    abundances at once (the off-diagonal YP-D/H term is not recoverable from
    the per-observable sigmas alone).  Only the products this network actually
    produced are shown; nothing is printed if fewer than two are present.

    The layout (fixed-width fields) is byte-for-byte identical to the C
    backend's ``print_mc_matrices`` (``primat-c/src/cli.c``), per CLAUDE.md's
    verbose/output-parity mandate.
    """
    names = mc.quantity_names()
    labels = [q for q in _MC_MAIN_PRODUCTS if q in names]
    if len(labels) < 2:
        return
    # A sample covariance/correlation needs at least 2 samples; skip the
    # matrices entirely for a degenerate single-sample run (they'd be all-NaN).
    if mc.samples_array().shape[0] < 2:
        return
    title = ", ".join(labels)
    # Correlation: 8-wide row-label column, 9-wide value columns (%9.3f); the
    # header names are right-justified to the same 9 so they sit over values.
    # A leading blank line separates the matrix block from the preceding
    # per-observable summary for readability (mirrored in the C backend's
    # print_mc_matrices, per the byte-for-byte output-parity mandate).
    print()
    print(f"Correlation matrix ({title}):")
    print(" " * 8 + "".join(f"{q:>9}" for q in labels))
    for a in labels:
        print(f"{a:>8}" + "".join(f"{mc.corr(a, b):9.3f}" for b in labels))
    # Covariance: same 8-wide labels, 13-wide value columns (%13.3e).
    # Leading blank line, as for the correlation block above.
    print()
    print(f"Covariance matrix ({title}):")
    print(" " * 8 + "".join(f"{q:>13}" for q in labels))
    for a in labels:
        print(f"{a:>8}" + "".join(f"{mc.cov(a, b):13.3e}" for b in labels))


def _print_list_params():
    """Print every ``PRIMATConfig`` parameter (``DEFAULT_PARAMS`` key), its
    default value, and a one-line description, for ``primat --list-params``
    (CLI discoverability). The ``--set KEY=VALUE`` escape hatch is
    deliberately left undocumented in ``--help`` to keep the printed help
    short (see the module docstring), so this is the intended way power
    users discover every parameter it accepts -- descriptions are parsed
    straight out of ``config.py``'s own inline comments
    (:func:`primat.config._default_params_comments`) rather than duplicated
    here, so they cannot drift out of sync with the source of truth. Keys
    are printed grouped and ordered by :data:`primat.config.PARAM_GROUPS`
    (the same taxonomy the GUI and the param-template generator use), rather
    than in raw ``DEFAULT_PARAMS`` dict order, so related flags sit together.
    """
    comments = _default_params_comments()
    key_width = max(len(key) for key in DEFAULT_PARAMS)
    for group, keys in PARAM_GROUPS.items():
        print(f"# ---- {group} ----")
        for key in keys:
            line = f"{key:<{key_width}} = {DEFAULT_PARAMS[key]!r}"
            comment = comments.get(key, "")
            if comment:
                line += f"  # {comment}"
            print(line)
        print()


def _build_parser():
    """Build the ``argparse.ArgumentParser`` for the ``primat`` CLI.

    Only the most commonly varied ``PRIMATConfig`` keys are exposed as flags;
    each flag's ``dest`` matches the corresponding config key so that
    ``main()`` can forward it to ``PRIMAT(params=...)`` unchanged.  Flags
    default to ``None`` (rather than duplicating ``PRIMATConfig``'s defaults)
    so that only options the user actually passed override the config.
    """
    parser = argparse.ArgumentParser(
        prog="primat",
        description="Run a Big Bang Nucleosynthesis computation with "
                     "primat and print the resulting Neff/abundances.",
        epilog="Any other PRIMATConfig parameter (including p_<reaction>/"
               "delta_<reaction> rate variations) can be set with "
               "repeated --set KEY=VALUE, e.g. --set T_end_MeV=1e-4. Use "
               "--list-params to see every parameter's default and a "
               "one-line description.",
    )
    # `version` action prints the string and exits before any computation;
    # the version itself comes from the installed distribution metadata via
    # primat.__version__ (single source of truth in pyproject.toml). Also
    # reports C-backend availability, since "primat 0.3.2" alone doesn't tell
    # a user whether --backend c will work on their install.
    parser.add_argument(
        "--credits", action="store_true",
        help="Print the project credits and exit.",
    )
    backend_status = "available" if HAS_C_BACKEND else "unavailable"
    parser.add_argument(
        "--version", action="version",
        version=f"%(prog)s {__version__} (C backend: {backend_status})",
        help="Print the primat version and C-backend availability, then exit.",
    )
    parser.add_argument(
        "--list-params", action="store_true",
        help="Print every PRIMATConfig parameter, its default value, and a "
             "one-line description, then exit.",
    )
    parser.add_argument(
        "--Omegabh2", type=float, default=None, metavar="VALUE",
        help=f"Baryon density Omega_b h^2 (PRIMATConfig default: "
             f"{DEFAULT_PARAMS['Omegabh2']}).",
    )
    parser.add_argument(
        "--DeltaNeff", type=float, default=None, metavar="VALUE",
        help=f"Extra relativistic degrees of freedom on top of the SM "
             f"neutrino sector (PRIMATConfig default: {DEFAULT_PARAMS['DeltaNeff']}).",
    )
    parser.add_argument(
        "--network", default=None, metavar="NAME",
        help=f"Nuclear reaction network used in the LT era "
             f"(PRIMATConfig default: {DEFAULT_PARAMS['network']!r}). Built-in "
             "choices are 'small', 'small_parthenope' and 'large', but any "
             "name for which data/nuclear/networks/<NAME>.txt exists is "
             "accepted; PRIMATConfig raises a ValueError if no such file is found.",
    )
    parser.add_argument(
        "--amax", type=int, default=None, metavar="A",
        help="Drop reactions involving any nuclide with mass number > A "
             "(must be a positive integer); applies to any --network, not "
             "just 'large'. E.g. --network large --amax 8 reproduces the "
             "old 'medium' network's 68 reactions.",
    )
    parser.add_argument(
        "--numerical_precision", type=float, default=None, metavar="RTOL",
        help=f"Relative tolerance passed to solve_ivp (PRIMATConfig default: "
             f"{DEFAULT_PARAMS['numerical_precision']}).",
    )
    parser.add_argument(
        "--data_dir", default=None, metavar="PATH",
        help="Replace the entire data tree (NEVO/, nuclear/, csv/, "
             "cache_plasma_weak/) with PATH. Must contain those "
             "subdirectories. Default (None): use the package's own "
             "primat/data/ tree.",
    )
    parser.add_argument(
        "--user_nuclear_dir", default=None, metavar="PATH",
        help="Additive overlay for nuclear networks and rate tables only "
             "(the primat/data/nuclear/ equivalent, containing networks/ and "
             "tables/). Checked before the default tree; shipped networks "
             "remain accessible even when this is set. Default: None.",
    )
    parser.add_argument(
        "--munuOverTnu", type=float, default=None, metavar="XI",
        help=f"Reduced neutrino chemical potential mu/T, the common default "
             f"for all flavours (PRIMATConfig default: {DEFAULT_PARAMS['munuOverTnu']}).",
    )
    parser.add_argument(
        "--munuOverTnu_e", type=float, default=None, metavar="XI_E",
        help="Per-flavour reduced chemical potential of nu_e; overrides "
             "--munuOverTnu for the electron neutrino (which alone shifts the "
             "n<->p weak rates). Default: inherit --munuOverTnu.",
    )
    parser.add_argument(
        "--munuOverTnu_mu", type=float, default=None, metavar="XI_MU",
        help="Per-flavour reduced chemical potential of nu_mu (gravitates only). "
             "Default: inherit --munuOverTnu.",
    )
    parser.add_argument(
        "--munuOverTnu_tau", type=float, default=None, metavar="XI_TAU",
        help="Per-flavour reduced chemical potential of nu_tau (gravitates only). "
             "Default: inherit --munuOverTnu.",
    )
    parser.add_argument(
        "--output_file", default=None, metavar="FILE",
        help="Write the full time-evolution TSV to FILE when "
             "--output_time_evolution is enabled.",
    )
    parser.add_argument(
        "--output_final_file", default=None, metavar="FILE",
        help="Write the final-abundance table to FILE when "
             "--output_final_result is enabled.",
    )
    parser.add_argument(
        "--output_background_file", default=None, metavar="FILE",
        help="Write the background time-evolution TSV to FILE when "
             "--output_background_evolution is enabled.",
    )
    parser.add_argument(
        "--output_mc_file_prefix", default=None, metavar="PREFIX",
        help="Filename stem for the Monte-Carlo output files written when --mc "
             "is used: PREFIX_samples.tsv / PREFIX_covariance.tsv / "
             "PREFIX_correlation.tsv, each gated by --output_mc_samples / "
             "--output_mc_covariance / --output_mc_correlation respectively.",
    )
    # Boolean PRIMATConfig flags exposed as --flag/--no-flag pairs (argparse's
    # BooleanOptionalAction), so they don't need the --set escape hatch. Each
    # flag's default is looked up from DEFAULT_PARAMS rather than duplicated
    # as a literal here, so the printed "(PRIMATConfig default: ...)" cannot
    # drift from config.py if that default is ever changed.
    for flag_name, help_text in (
        ("QED_corrections",
         "QED interaction corrections to the EM plasma equation of state."),
        ("nuclear_qed_corrections",
         "QED corrections to radiative-capture nuclear reaction rates "
         "(Pitrou & Pospelov 2020)."),
        ("radiative_corrections",
         "Coulomb + T=0 resummed radiative corrections to n<->p (CCR); "
         "if False, use the Born approximation instead."),
        ("finite_mass_corrections",
         "Finite-nucleon-mass (Fokker-Planck) correction to n<->p."),
        ("thermal_corrections",
         "Finite-temperature radiative corrections to n<->p (CCRTh; "
         "Brown & Sawyer 2001)."),
        ("spectral_distortions",
         "Correct n<->p rates for non-Fermi-Dirac neutrino distributions."),
        ("output_time_evolution",
         "Write the full time-evolution series (in-memory always; to disk "
         "if output_file is set)."),
        ("output_final_result",
         "Write the final results dict to output_file."),
        ("output_background_evolution",
         "Write the cosmological background time series to disk."),
        ("output_mc_samples",
         "Write --mc samples to <output_mc_file_prefix>_samples.tsv."),
        ("output_mc_covariance",
         "Write the --mc sample covariance matrix to "
         "<output_mc_file_prefix>_covariance.tsv."),
        ("output_mc_correlation",
         "Write the --mc sample correlation matrix to "
         "<output_mc_file_prefix>_correlation.tsv."),
        ("show_progress",
         "Print compact stderr progress indicators ('[primat]  HT.  MT.  "
         "LT.  done.' phase markers, '[MC] ...' sample counter)."),
    ):
        parser.add_argument(
            f"--{flag_name}", action=argparse.BooleanOptionalAction, default=None,
            help=f"{help_text} (PRIMATConfig default: {DEFAULT_PARAMS[flag_name]}).",
        )
    parser.add_argument(
        "--backend", choices=("auto", "c", "python"), default="auto",
        help="Which solver implementation to use: 'auto' (default) picks the "
             "compiled C extension when available, 'c' forces it (error if "
             "unavailable), 'python' forces the pure-Python implementation.",
    )
    parser.add_argument(
        "--mc", type=int, default=None, metavar="N",
        help="Also run an N-sample Monte-Carlo nuclear-rate/tau_n uncertainty "
             "propagation (primat.backend.run_mc) and print each observable "
             "as 'value +/- sigma' instead of a bare value. Uses the C "
             "backend when available (--backend auto/c), else joblib in "
             "pure Python; see primat.backend's docstring for the RNG caveat.",
    )
    parser.add_argument(
        "--mc-seed", type=int, default=0, metavar="SEED",
        help="Base RNG seed for --mc (default: 0); sample i uses seed+i.",
    )
    parser.add_argument(
        "--mc-jobs", type=int, default=-1, metavar="N",
        help="Parallel worker count for --mc (default: -1, all CPUs; Python "
             "backend only -- the C backend always uses one pthread per sample).",
    )
    parser.add_argument(
        "--json", action="store_true",
        help="Print the full results dict as JSON instead of a short summary.",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable primat's internal progress messages (timings, cache hits, ...).",
    )
    parser.add_argument(
        "--cache-info", action="store_true",
        help="Print the number of cached n<->p weak-rate files "
             "(cache_plasma_weak/weak/nTOp_*.txt, or the cache_dir redirect) "
             "and exit, without running a solve.",
    )
    parser.add_argument(
        "--cache-clear", action="store_true",
        help="Delete every cached n<->p weak-rate file and exit, without "
             "running a solve. The cache is always safely regenerable: a "
             "later run just pays the one-time recompute cost again.",
    )
    # Generic escape hatch: lets any PRIMATConfig key (including p_<reaction>/
    # delta_<reaction>) be set from the CLI without a dedicated flag.
    # help=SUPPRESS keeps it out of --help, per the handful of named flags
    # above being the only ones intended to show there; see the module
    # docstring for usage.
    parser.add_argument(
        "--set", action="append", dest="set_params", metavar="KEY=VALUE",
        default=[], help=argparse.SUPPRESS,
    )
    return parser


def main(argv=None):
    """Entry point for the ``primat`` console script.

    Parses command-line arguments into a ``PRIMATConfig`` ``params`` dict,
    runs ``PRIMAT(params).primat_results()``, and prints either a short
    human-readable summary (default) or the full results dict as JSON
    (``--json``).

    Parameters
    ----------
    argv : list of str, optional
        Argument vector to parse; defaults to ``sys.argv[1:]``. Exposed as a
        parameter so the CLI can be invoked programmatically (e.g. in tests)
        without spawning a subprocess.

    Returns
    -------
    int
        Process exit code, always ``0`` on success (argparse itself exits
        with code 2 on a bad argument).

    Example
    -------
        $ primat --Omegabh2 0.02242 --network large --amax 8
        Neff       = 3.04397730
        YP (BBN)   = 0.24691900
        ...
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.credits:
        print(cli_credits_text())
        return 0

    if args.list_params:
        _print_list_params()
        return 0

    if args.cache_info or args.cache_clear:
        cfg = PRIMATConfig({})
        if args.cache_clear:
            n = clear_weak_cache(cfg)
            print(f"Removed {n} cached weak-rate file(s) from {weak_cache_dir(cfg)}/.")
        else:
            n = len(list_weak_cache_files(cfg))
            print(f"{n} cached weak-rate file(s) in {weak_cache_dir(cfg)}/.")
        return 0

    # Only forward options the user actually set, so unset flags fall back
    # to PRIMATConfig's own defaults rather than a value duplicated here.
    params = {}
    for key in (
        "Omegabh2", "DeltaNeff", "network", "amax", "numerical_precision",
        "munuOverTnu", "munuOverTnu_e", "munuOverTnu_mu", "munuOverTnu_tau",
        "QED_corrections", "nuclear_qed_corrections",
        "radiative_corrections", "finite_mass_corrections", "thermal_corrections",
        "spectral_distortions", "output_time_evolution", "output_final_result",
        "output_background_evolution", "output_mc_samples",
        "output_mc_covariance", "output_mc_correlation", "show_progress",
        "output_file", "output_final_file", "output_background_file",
        "output_mc_file_prefix", "data_dir", "user_nuclear_dir",
    ):
        value = getattr(args, key)
        if value is not None:
            params[key] = value
    if args.verbose:
        params["verbose"] = True
    for entry in args.set_params:
        if "=" not in entry:
            parser.error(f"--set {entry!r}: expected KEY=VALUE")
        key, _, raw_value = entry.partition("=")
        params[key] = _parse_set_value(raw_value)

    for key in ("data_dir", "user_nuclear_dir"):
        if params.get(key) is not None:
            print(_rates_overlay_notice(key, params[key]), file=sys.stderr)

    start_time = time.time()
    if args.mc is not None:
        # run_mc already computes the central (nominal) solve internally;
        # derive results from mc[q].central to avoid a redundant second solve.
        mc = run_mc(args.mc, params=params, force_backend=args.backend,
                    seed=args.mc_seed, n_jobs=args.mc_jobs)
        # Flat dict: each quantity's central value plus a sigma_<name> entry
        # holding its MC uncertainty (see MCResult.to_flat_dict).
        results = mc.to_flat_dict()
    else:
        results = run_bbn(params=params, force_backend=args.backend)
        mc = None
        # output_mc_samples/output_mc_covariance/output_mc_correlation only
        # have an effect inside the `if mc is not None:` file-writing block
        # below (mc_uncertainty/run_mc is what produces the MCResult they are
        # dumped from); without --mc there is no MCResult, so any of these
        # flags being set is silently a no-op unless we flag it here.
        cfg_check = PRIMATConfig(params)
        requested = [name for name, enabled in (
            ("output_mc_samples", cfg_check.output_mc_samples),
            ("output_mc_covariance", cfg_check.output_mc_covariance),
            ("output_mc_correlation", cfg_check.output_mc_correlation),
        ) if enabled]
        if requested:
            # Plain stderr print (matching _rates_overlay_notice above and
            # the C CLI's fprintf(stderr, "warning: ...") in cli.c), not
            # warnings.warn: this fires inside main() itself, the CLI entry
            # point, so warnings.warn's source-line echo would just quote
            # whichever line of this (multi-line) call happens to sit at the
            # reported lineno -- noise, not information -- for a warning that
            # has no meaningful caller frame to attribute to anyway.
            print(f"warning: {', '.join(requested)} set but --mc was not "
                  "passed; no MC output file(s) will be written.",
                  file=sys.stderr)
    elapsed = time.time() - start_time

    if args.json:
        out = dict(results)
        if mc is not None:
            out["mc"] = {q: {"central": mc[q].central, "mean": mc[q].mean,
                              "std": mc[q].std, "values": list(mc[q].values)}
                         for q in mc.quantity_names()}
        print(json.dumps(out, indent=2))
    else:
        T_end_MeV = params.get("T_end_MeV", 1e-3)
        sep = "─" * 52
        header = f"PRIMAT results at T = {T_end_MeV:g} MeV"
        print(sep)
        print(f"{header:^52}")
        print(sep)
        print(f"Neff       = {results['Neff']:.8f}" +
              (f" +/- {mc['Neff'].std:.8f}" if mc is not None and "Neff" in mc.quantity_names() else ""))
        print(f"YP (BBN)   = {results['YPBBN']:.8f}" +
              (f" +/- {mc['YPBBN'].std:.8f}" if mc is not None and "YPBBN" in mc.quantity_names() else ""))
        print(f"YP (CMB)   = {results['YPCMB']:.8f}" +
              (f" +/- {mc['YPCMB'].std:.8f}" if mc is not None and "YPCMB" in mc.quantity_names() else ""))
        print(f"D/H        = {results['DoH']:.7e}" +
              (f" +/- {mc['DoH'].std:.7e}" if mc is not None and "DoH" in mc.quantity_names() else ""))
        print(f"He3/H      = {results['He3oH']:.7e}" +
              (f" +/- {mc['He3oH'].std:.7e}" if mc is not None and "He3oH" in mc.quantity_names() else ""))
        print(f"He3/He4    = {results['He3oHe4']:.7e}" +
              (f" +/- {mc['He3oHe4'].std:.7e}" if mc is not None and "He3oHe4" in mc.quantity_names() else ""))
        print(f"Li7/H      = {results['Li7oH']:.6e}" +
              (f" +/- {mc['Li7oH'].std:.6e}" if mc is not None and "Li7oH" in mc.quantity_names() else ""))
        if "Li6oLi7" in results:
            print(f"Li6/Li7    = {results['Li6oLi7']:.6e}" +
                  (f" +/- {mc['Li6oLi7'].std:.6e}" if mc is not None and "Li6oLi7" in mc.quantity_names() else ""))
        if "YCNO" in results:
            print(f"CNO (mass) = {results['YCNO']:.6e}" +
                  (f" +/- {mc['YCNO'].std:.6e}" if mc is not None and "YCNO" in mc.quantity_names() else ""))
        # After the per-observable value +/- sigma block, print the joint
        # uncertainty of the four main products (YPBBN/DoH/He3oHe4/Li7oH): the
        # 4x4 correlation matrix (dimensionless, 3 decimals) and covariance
        # matrix (%.3e), computed from the same MC samples. The correlation is
        # what a user needs to build a joint likelihood over several abundances.
        if mc is not None:
            _print_mc_matrices(mc)
        print(f"--- running time: {elapsed:.2f} seconds ---")

    if mc is not None:
        # All three MC files share one filename stem (output_mc_file_prefix),
        # each gated by its own boolean, and are written verbatim from the same
        # dump_mc_* helpers the GUI uses (backend-agnostic MCResult).
        cfg_check = PRIMATConfig(params)
        prefix = cfg_check.output_mc_file_prefix
        for enabled, suffix, writer, label in (
            (cfg_check.output_mc_samples,     "_samples.tsv",     dump_mc_samples,     "samples"),
            (cfg_check.output_mc_covariance,  "_covariance.tsv",  dump_mc_covariance,  "covariance"),
            (cfg_check.output_mc_correlation, "_correlation.tsv", dump_mc_correlation, "correlation"),
        ):
            if not enabled:
                continue
            out_path = os.path.abspath(prefix + suffix)
            out_dir = os.path.dirname(out_path)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            with open(out_path, "w") as f:
                f.write(writer(mc))
            if label == "samples":
                sample_word = "sample" if args.mc == 1 else "samples"
                print(f"[output] MC samples ({args.mc} {sample_word}) written to {out_path}")
            else:
                print(f"[output] MC {label} matrix written to {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
