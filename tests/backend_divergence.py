# -*- coding: utf-8 -*-
"""Cross-backend divergence harness: where C and Python disagree, term by term.

``tests/README.md``'s "Known cross-backend divergences" names the causes;
``tests/test_backend_parity.py`` pins the observables' total gap. This module
is the instrument in between: it measures each *link of the chain* separately
-- background, nuclear rate tables, CCRTh interpolation, per-nuclide
abundances -- so a widened total can be attributed instead of merely noticed.

Run it directly for a report::

    python -m tests.backend_divergence                  # small + large,amax=8
    python -m tests.backend_divergence --precision 1e-9 # converged tolerance

Every number is a relative difference ``(C - Python) / Python``, reported over
the BBN window T = 0.008..1.2 MeV where the light-element abundances are set.
"""
from __future__ import annotations

import numpy as np

from primat.backend import run_bbn

# The temperature window that sets the light-element abundances: D burns in
# around 0.1 MeV and the network freezes out below ~0.01 MeV.
BBN_WINDOW_MEV = (0.008, 1.2)

# Observables reported by both backends, in the order the report prints them.
OBSERVABLES = ("YPBBN", "DoH", "He3oH", "Li7oH", "Neff")


def _rel(c, p):
    """Relative difference (C - Python) / Python, 0 where Python's value is 0."""
    c, p = np.asarray(c, float), np.asarray(p, float)
    out = np.zeros_like(p)
    nz = p != 0.0
    out[nz] = (c[nz] - p[nz]) / p[nz]
    return out


def _stats(rel, weights=None):
    """(max|rel|, median|rel|, signed mean) of a relative-difference array."""
    a = np.abs(rel)
    if a.size == 0:
        return dict(max=0.0, median=0.0, mean=0.0)
    return dict(max=float(a.max()), median=float(np.median(a)),
                mean=float(np.mean(rel)))


def observable_gap(params, precision=None):
    """Cross-backend gap on the five reported observables, plus per-nuclide Y.

    Returns ``{"observables": {name: rel}, "observables_value": {name: python
    value}, "nuclides": {name: rel}}``. Per-nuclide entries below ``1e-25``
    abundance per baryon are skipped as noise.
    """
    p = dict(params)
    if precision is not None:
        p["numerical_precision"] = precision
    r_c = run_bbn(dict(p), force_backend="c")
    r_py = run_bbn(dict(p), force_backend="python")
    obs, val = {}, {}
    for k in OBSERVABLES:
        if k in r_c and k in r_py:
            obs[k] = float(_rel(r_c[k], r_py[k]))
            val[k] = float(r_py[k])
    nuc = {}
    for s, v in r_py.get("Y_final", {}).items():
        if v > 1e-25 and s in r_c["Y_final"]:
            nuc[s] = float(_rel(r_c["Y_final"][s], v))
    return {"observables": obs, "observables_value": val, "nuclides": nuc}


def background_gap(params, tmp_dir, precision=None):
    """Gap in a(T), t(T), H(T) and the neutrino temperature, at matched T.

    Both backends write ``output_background.tsv`` on a grid log-spaced in
    *cosmic time*, so the two T columns do not coincide: the Python table is
    re-interpolated (log-log cubic) onto C's temperatures before comparing.
    ``tmp_dir`` receives the two TSVs.
    """
    from scipy.interpolate import CubicSpline

    tabs = {}
    for backend in ("c", "python"):
        path = f"{tmp_dir}/background_{backend}.tsv"
        p = dict(params, output_background_evolution=True,
                 output_background_file=path)
        if precision is not None:
            p["numerical_precision"] = precision
        run_bbn(p, force_backend=backend)
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        tabs[backend] = (header, np.loadtxt(path, skiprows=1))

    header, c = tabs["c"]
    _, p_tab = tabs["python"]
    Tc, Tp = c[:, 0], p_tab[:, 0]
    order = np.argsort(Tp)
    # Stay strictly inside the Python table so the comparison never rides on
    # spline extrapolation at either end.
    inside = (Tc >= Tp.min() * 1.001) & (Tc <= Tp.max() * 0.999)
    window = inside & (Tc > BBN_WINDOW_MEV[0]) & (Tc < BBN_WINDOW_MEV[1])

    out = {}
    for j, name in enumerate(header):
        key = name.split()[0]
        if key not in ("t", "a", "H", "Tnue"):
            continue
        yc, yp = c[:, j], p_tab[:, j]
        if not (np.all(yc > 0) and np.all(yp > 0)):
            continue
        interp = CubicSpline(np.log(Tp[order]), np.log(yp[order]))
        rel = _rel(yc, np.exp(interp(np.log(Tc))))
        out[key] = _stats(rel[window])
    return out


def rate_column_gap(params, precision=None):
    """Gap in the per-reaction forward nuclear rates, at matched T.

    Uses the optional ``<reaction>_frwrd`` evolution columns both backends
    emit, interpolating Python's onto C's temperatures. Returns
    ``{reaction: stats}`` over the BBN window.
    """
    from scipy.interpolate import interp1d

    p = dict(params, output_time_evolution=True,
             output_rates_time_evolution=True, output_file=None)
    if precision is not None:
        p["numerical_precision"] = precision
    evo_c = run_bbn(dict(p), force_backend="c")["evolution"]
    evo_py = run_bbn(dict(p), force_backend="python")["evolution"]

    Tc, Tp = evo_c.T_gamma, evo_py.T_gamma
    order = np.argsort(Tp)
    window = ((Tc >= max(Tc.min(), Tp.min())) & (Tc <= min(Tc.max(), Tp.max()))
              & (Tc > BBN_WINDOW_MEV[0]) & (Tc < BBN_WINDOW_MEV[1]))
    out = {}
    for name in evo_c.rates:
        ip = interp1d(Tp[order], evo_py.rates[name][order], kind="cubic",
                      fill_value="extrapolate")(Tc)
        # Several shipped tables are exactly 0 over most of the grid and switch
        # on abruptly; interpolating across that hard zero flips between 0 and
        # ~1e-20, which is a threshold artefact, not a rate disagreement. The
        # floor scales with the column (same convention as
        # test_backend_parity.test_rates_columns_backend_parity).
        floor = 1e-12 * float(np.max(np.abs(evo_c.rates[name])))
        keep = window & (np.abs(ip) > floor)
        if keep.any():
            out[name] = _stats(_rel(evo_c.rates[name][keep], ip[keep]))
    return out


def ccrth_interpolant_gap(params=None):
    """Spread between the two backends' CCRTh interpolants.

    The finite-temperature (CCRTh) correction is read from a cache the two
    backends share, so only the curve *between* its nodes can differ. Both fit
    a not-a-knot cubic in linear (T, L) space -- Python
    ``interp1d(kind='cubic')``, C ``cpr_cubic_spline_fit_notaknot``, two
    implementations of the same mathematical object. This evaluates Python's
    against an independent scipy not-a-knot fit at the log-midpoints between
    nodes (the worst case), as a fraction of the n->p rate: the scale both weak
    terms enter the network's right-hand side against.

    It guards the *scheme*. A mismatch here is invisible to any end-to-end
    comparison -- both curves pass through the shared nodes -- and moves the
    observables. Needs no C backend.
    """
    from scipy.interpolate import CubicSpline

    from primat import PRIMAT
    from primat.weak_rates.cache import thermal_fingerprint, fingerprint_hash
    from primat.cache_utils import resolve_cache_file

    run = PRIMAT(dict(params or {"network": "small"}))
    cfg = run.cfg
    if not cfg.thermal_corrections:
        return None
    fname = "nTOp_thermal_%s.txt" % fingerprint_hash(thermal_fingerprint(cfg))
    tab = np.loadtxt(resolve_cache_file(cfg, "weak", fname))
    T, Ln, Lp = tab[:, 0], tab[:, 1], tab[:, 2]
    mid = np.sqrt(T[:-1] * T[1:])
    T_MeV = mid / cfg.MeV_to_Kelvin
    mid = mid[(T_MeV > BBN_WINDOW_MEV[0]) & (T_MeV < BBN_WINDOW_MEV[1])]
    scale = run.background.weak_nTOp(mid) * cfg.tau_n

    # The live closures, so a change of scheme in corrections.py shows up here
    # rather than being re-asserted by this file. They return L/Fn, the units
    # the non-thermal table (and `scale`) are in.
    from primat.weak_rates.corrections import (ComputeFn,
                                                _thermal_correction_interpolants)
    inv_Fn = 1.0 / ComputeFn(cfg)
    fn, fp = _thermal_correction_interpolants(
        [run.background.Tg_vec, run.background.Tnue_vec], cfg)

    out = {}
    for name, L, f in (("n_to_p", Ln, fn), ("p_to_n", Lp, fp)):
        reference = CubicSpline(T, L * inv_Fn, bc_type="not-a-knot")(mid)
        out[name] = _stats((f(mid) - reference) / scale)
    return out


def report(params, tmp_dir, precision=None):
    """Print the full divergence report for one configuration."""
    label = f"{params.get('network', 'small')}"
    if params.get("amax"):
        label += f", amax={params['amax']}"
    prec = precision if precision is not None else "default (1e-7)"
    print(f"\n=== {label}   numerical_precision = {prec}")

    g = observable_gap(params, precision)
    print("  observables (C - py)/py:")
    for k, v in g["observables"].items():
        extra = ""
        if k == "YPBBN":
            # YP is quoted absolutely everywhere else (tests/README.md's
            # validation reference, test_backend_parity's budget table).
            extra = f"   [abs {v * g['observables_value'][k]:+.3e}]"
        print(f"    {k:8s} {v:+.3e}{extra}")
    worst = sorted(g["nuclides"].items(), key=lambda kv: -abs(kv[1]))[:5]
    print("  worst per-nuclide Y_final: "
          + ", ".join(f"{s}={v:+.2e}" for s, v in worst))

    print("  background at matched T (BBN window):")
    for k, s in background_gap(params, tmp_dir, precision).items():
        print(f"    {k:8s} max={s['max']:.2e}  median={s['median']:.2e}  "
              f"mean={s['mean']:+.2e}")

    rates = rate_column_gap(params, precision)
    worst_r = sorted(rates.items(), key=lambda kv: -kv[1]["max"])[:3]
    print(f"  nuclear rate columns ({len(rates)} reactions), worst 3 by max:")
    for name, s in worst_r:
        print(f"    {name:24s} max={s['max']:.2e}  median={s['median']:.2e}")

    th = ccrth_interpolant_gap(params)
    if th:
        print("  CCRTh interpolant vs an independent not-a-knot cubic fit, "
              "as a fraction of the n->p rate:")
        for k, s in th.items():
            print(f"    {k:8s} max={s['max']:.2e}  median={s['median']:.2e}")


def main(argv=None):
    import argparse
    import tempfile

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--precision", type=float, default=None,
                    help="numerical_precision for both backends (default: the "
                         "package default, 1e-7)")
    ap.add_argument("--network", default=None,
                    help="run one network only (default: small and large/amax=8)")
    args = ap.parse_args(argv)

    configs = ([{"network": args.network}] if args.network
               else [{"network": "small"}, {"network": "large", "amax": 8}])
    with tempfile.TemporaryDirectory() as tmp:
        for cfg in configs:
            report(cfg, tmp, args.precision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
