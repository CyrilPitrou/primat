# -*- coding: utf-8 -*-
"""
nuclide_table.py
================
Offline helpers (used only by the rate/network *generation* command, never at
primat run time) that turn the reaction list extracted from AC2024 +
PRIMAT-main.m into:

  1. the **set of nuclides** the network touches, each resolved to its
     (N, Z, A, charge Q, mass excess, spin), and
  2. the **detailed-balance coefficients** (alpha, beta, gamma) of every
     reaction that has a reverse rate.

Why offline.  These quantities never change for a fixed PRIMAT version, and the
NUBASE2020 table is ~760 kB, so we resolve everything once here and bake the
result into small CSV files that primat simply reads at start-up
(``nuclides.csv``, ``detailed_balance.csv``).

Token convention.  Reaction sides come from the AC2024/PRIMAT sources as token
lists that mix spellings: ``a``/``He4``, ``d``/``H2``, ``t``/``H3``, the bare
nucleons ``n``/``p``, ordinary nuclides ``Be9``/``C12``/..., the photon ``g``
(or ``2g``), and the beta-decay leptons ``Bm`` (electron, e^-) and ``Bp``
(positron, e^+).  :func:`resolve_token` maps any of these to a canonical record;
nuclides are keyed by a canonical name (``n``, ``p``, ``H2``, ``H3``, ``He3``,
``He4``, ``Be9``, ...) chosen to match primat's existing ``Nuclides`` keys.
"""
import re
from collections import Counter

# Element symbol -> atomic number Z.  The BBN+ network reaches Na (Z=11); we
# list through Ca (Z=20) so the table comfortably covers any token that appears.
_ELEMENT_Z = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
}
_Z_ELEMENT = {z: sym for sym, z in _ELEMENT_Z.items()}

# Short single-letter aliases used in the sources, mapped to (Z, A).
_SHORT = {"n": (0, 1), "p": (1, 1), "d": (1, 2), "t": (1, 3), "a": (2, 4)}


class Token:
    """Resolved reaction token.

    kind : 'nuclide' | 'photon' | 'lepton'
    Z, A : atomic number and mass number (0 for photon/lepton).
    Q    : electric charge in units of e (Z for nuclei; Bm=-1, Bp=+1; 0 otherwise).
    name : canonical nuclide name (None for photon/lepton).
    """
    __slots__ = ("kind", "Z", "A", "Q", "name")

    def __init__(self, kind, Z, A, Q, name):
        self.kind, self.Z, self.A, self.Q, self.name = kind, Z, A, Q, name


def canonical_name(Z, A):
    """Canonical nuclide name from (Z, A), matching primat's ``Nuclides`` keys.

    The bare nucleons are special-cased to ``n``/``p``; everything else is the
    element symbol followed by the mass number (``H2``, ``H3``, ``He4``,
    ``Be9``, ``C12``, ...).  This deliberately yields ``H2`` (not ``d``) so the
    generated ``nuclides.csv`` lines up with the names primat already uses.
    """
    if (Z, A) == (0, 1):
        return "n"
    if (Z, A) == (1, 1):
        return "p"
    return f"{_Z_ELEMENT[Z]}{A}"


def resolve_token(tok):
    """Resolve one reaction token (e.g. ``'a'``, ``'He4'``, ``'Be9'``, ``'Bm'``,
    ``'g'``) to a :class:`Token`.

    Photons (``g``) carry no baryon number or charge; the beta leptons ``Bm``
    (e^-) and ``Bp`` (e^+) carry charge -1 / +1 and A=0 -- both are needed for
    the formal charge/baryon conservation check but are *not* tracked species.
    """
    if tok in ("g", "gamma"):
        return Token("photon", 0, 0, 0, None)
    if tok == "Bm":                       # beta-minus: emits an electron
        return Token("lepton", 0, 0, -1, None)
    if tok == "Bp":                       # beta-plus: emits a positron
        return Token("lepton", 0, 0, +1, None)
    if tok in _SHORT:
        Z, A = _SHORT[tok]
        return Token("nuclide", Z, A, Z, canonical_name(Z, A))
    m = re.fullmatch(r"([A-Z][a-z]?)(\d+)", tok)
    if not m:
        raise ValueError(f"cannot resolve reaction token {tok!r}")
    sym, A = m.group(1), int(m.group(2))
    if sym not in _ELEMENT_Z:
        raise ValueError(f"unknown element symbol {sym!r} in token {tok!r}")
    Z = _ELEMENT_Z[sym]
    return Token("nuclide", Z, A, Z, canonical_name(Z, A))


# ---------------------------------------------------------------------------
# NUBASE2020 reader (general: keyed by (Z, A), not a fixed nuclide list)
# ---------------------------------------------------------------------------
def _parse_spin(jpi_field):
    """Leading J (integer or fraction) of a NUBASE ``Jpi`` field (``3/2-*`` -> 1.5)."""
    m = re.match(r"\s*\(?(\d+)(?:/(\d+))?", jpi_field)
    if m is None:
        return None
    return float(m.group(1)) / float(m.group(2)) if m.group(2) else float(m.group(1))


def load_nubase_all(nubase_path):
    """Read mass excesses and spins of *every* ground state in a NUBASE2020 file.

    Fixed-width columns (0-indexed): A = ``[0:3]``, Z = ``[4:7]``, isomer index
    ``[7]`` (``0`` = ground state), mass excess [keV] = ``[18:31]``, Jpi =
    ``[88:102]``.  Estimated values carry a trailing ``#`` we strip.

    Returns ``{(Z, A): (mass_excess_keV, spin_J)}``.
    """
    table = {}
    with open(nubase_path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#") or len(line) < 102:
                continue
            if line[7] != "0":                      # keep ground states only
                continue
            try:
                A = int(line[0:3])
                Z = int(line[4:7])
            except ValueError:
                continue
            excess = float(line[18:31].replace("#", ""))
            table[(Z, A)] = (excess, _parse_spin(line[88:102]))
    return table


# Half-life units NUBASE spells out in its 2-character unit field, in seconds.
# The year is the Julian year (365.2422 d), matching the conversion factor
# convert_ac2024_rates.py's coded half-lives use.
_UNIT_TO_S = {
    "ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15,
    "ps": 1e-12, "ns": 1e-9,  "us": 1e-6,  "ms": 1e-3,
    "s":  1.0,   "m":  60.0,  "h":  3600.0,
    "d":  86400.0,
    "y":  86400.0 * 365.2422,
    "ky": 86400.0 * 365.2422 * 1e3,
    "My": 86400.0 * 365.2422 * 1e6,
    "Gy": 86400.0 * 365.2422 * 1e9,
    "Ty": 86400.0 * 365.2422 * 1e12,
}


def load_nubase_halflives(nubase_path):
    """Read every ground state's half-life [s] from a NUBASE2020 file.

    Companion to :func:`load_nubase_all` (which returns mass excesses and
    spins).  Used by ``convert_ac2024_rates._validate_decay_halflives`` to
    cross-check the half-lives hard-coded in ``_ANALYTIC_REACTIONS`` against
    the evaluation, so a copy-paste error in a decay rate is caught at
    generation time rather than shipped.

    Fixed-width layout, quoted verbatim from the format block at the top of
    ``nubase_4.mas20.txt`` (its columns are **1-based**, so each slice below is
    ``[start-1 : stop]``)::

        70: 78   T #         f9.4   Half-life ("stbl", "p-unst", or a value)
        79: 80   unit T        a2   Half-life unit

    i.e. value = ``line[69:78]`` (nine characters), unit = ``line[78:80]``.
    Reading the value one column late silently drops the leading digit of any
    half-life wide enough to fill the field -- 8 ground states in the shipped
    table are, e.g. Ne18's ``"1664.20  ms"``, which then reads as 664.20 ms.

    Args:
        nubase_path: str, path to the NUBASE2020 fixed-width text file
            (``generate_rates/nubase_4.mas20.txt``).

    Returns:
        dict ``{(Z, A): half_life_seconds_or_None}``.  ``None`` means "no
        usable measured half-life": stable, particle-unstable, a bare limit
        (``">912.4 ys"`` -- a bound, not a value, so comparing a coded number
        against it would be meaningless), or an unrecognised unit.

    Example:
        >>> t12 = load_nubase_halflives("generate_rates/nubase_4.mas20.txt")
        >>> round(t12[(10, 18)], 4)          # Ne18 -> F18, 1664.20 ms
        1.6642
    """
    halflives = {}
    with open(nubase_path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#") or len(line) < 82:
                continue
            if line[7] != "0":                      # ground states only
                continue
            try:
                A = int(line[0:3])
                Z = int(line[4:7])
            except ValueError:
                continue
            t_str = line[69:78].strip().replace("#", "")   # '#' = systematics
            # For plain units ("s", "d", "y") the field's first slot is a
            # space; for "ms"/"ky"/"My"/"Gy" the SI prefix sits there.
            unit = line[78:80].strip()
            if t_str in ("stbl", "p-unst", "") or t_str[0] in "<>":
                halflives[(Z, A)] = None
                continue
            try:
                t_val = float(t_str)
            except ValueError:
                halflives[(Z, A)] = None
                continue
            s_per_unit = _UNIT_TO_S.get(unit)
            halflives[(Z, A)] = None if s_per_unit is None else t_val * s_per_unit
    return halflives


def build_nuclide_table(reactions, nubase_path):
    """Deduce the nuclide set from the reaction list and attach NUBASE properties.

    ``reactions`` is a list of dicts with ``reactants``/``products`` token
    lists.  Photons and leptons are dropped; every remaining distinct nuclide is
    resolved to (N, Z, A, Q) and matched to its NUBASE mass excess [keV] and
    ground-state spin.

    Returns an ``Ordered‑ish`` dict ``name -> record`` where record is
    ``dict(name, N, Z, A, Q, excess_keV, spin)``, ordered by increasing (Z, A)
    so the file is deterministic and human-readable.
    """
    nubase = load_nubase_all(nubase_path)
    seen = {}
    for rxn in reactions:
        for tok in rxn["reactants"] + rxn["products"]:
            t = resolve_token(tok)
            if t.kind != "nuclide":
                continue
            if t.name in seen:
                continue
            if (t.Z, t.A) not in nubase:
                raise ValueError(
                    f"nuclide {t.name} (Z={t.Z}, A={t.A}) not found in NUBASE "
                    f"file {nubase_path}")
            excess, spin = nubase[(t.Z, t.A)]
            seen[t.name] = dict(name=t.name, N=t.A - t.Z, Z=t.Z, A=t.A, Q=t.Z,
                                excess_keV=excess, spin=spin)
    # Deterministic order: by (Z, A), so n, p, H2, H3, He3, He4, ... come first.
    return {rec["name"]: rec
            for rec in sorted(seen.values(), key=lambda r: (r["Z"], r["A"]))}


# ---------------------------------------------------------------------------
# Formal conservation check (baryon number A and electric charge Q)
# ---------------------------------------------------------------------------
def conservation_residual(reactants, products):
    """Return ``(dA, dQ)`` = products-minus-reactants of (baryon number, charge).

    Both must be 0 for a physical reaction.  Photons contribute (0, 0); the beta
    leptons contribute (0, -1) for ``Bm`` and (0, +1) for ``Bp``.  This is a
    *formal* (exact integer) check -- no floating point involved.
    """
    def totals(side):
        A = Q = 0
        for tok in side:
            t = resolve_token(tok)
            A += t.A
            Q += t.Q
        return A, Q
    Ar, Qr = totals(reactants)
    Ap, Qp = totals(products)
    return Ap - Ar, Qp - Qr


# ---------------------------------------------------------------------------
# Detailed balance (same physics as primat.network_data.compute_detailed_
# balance_coefficients, reimplemented standalone here so this offline
# generator has no runtime dependency on a live PRIMATConfig instance)
# ---------------------------------------------------------------------------
class _DBConfig:
    """Minimal stand-in for ``PRIMATConfig`` exposing exactly what
    :func:`nuclear_data.detailed_balance` reads: the nuclide property dicts
    (built here for the *whole* large network) and the fundamental constants
    (copied verbatim from ``primat.constants.CONST``, the frozen single
    source of truth -- avoids instantiating a full, throwaway ``PRIMATConfig``
    just to read seven numbers).  This lets the
    offline generator compute detailed balance over an arbitrary nuclide set
    without needing a full PRIMATConfig."""

    def __init__(self, nuclide_table):
        from primat.constants import CONST
        for k in ("keV", "kB", "MeV", "ma", "me", "clight", "hbar"):
            setattr(self, k, getattr(CONST, k))
        self.Nuclides = {n: [r["N"], r["Z"]] for n, r in nuclide_table.items()}
        self.NuclExcessMass = {n: r["excess_keV"] for n, r in nuclide_table.items()}
        self.NuclSpin = {n: r["spin"] for n, r in nuclide_table.items()}


def make_detailed_balance(nuclide_table):
    """Return ``db(reactants, products) -> (Q_keV, alpha, beta, gamma)``.

    ``reactants``/``products`` are token lists (any spelling); photons and
    leptons are dropped, the rest canonicalised, and the result handed to
    :func:`nuclear_data.detailed_balance`.  ``Q_keV`` is the
    energy released (positive = exothermic).  Reactions that emit a lepton
    (decays) have no reverse rate and must not be passed here.
    """
    from nuclear_data import detailed_balance
    cfg = _DBConfig(nuclide_table)

    def to_canonical(side):
        out = []
        for tok in side:
            t = resolve_token(tok)
            if t.kind == "nuclide":           # drop photons (no mass/spin)
                out.append(t.name)
        return out

    def db(reactants, products):
        rc, pc = to_canonical(reactants), to_canonical(products)
        alpha, beta, gamma = detailed_balance(rc, pc, cfg)
        # gamma = -Q/(kB*1e9 K); recover Q in keV for the stored table.
        Q_keV = -gamma * cfg.kB * 1e9 / cfg.keV
        return Q_keV, alpha, beta, gamma

    return db


def is_decay(reactants, products):
    """True if the reaction emits a beta lepton (``Bm``/``Bp``) -- i.e. it is a
    weak decay with no reverse rate, so it has no detailed-balance coefficients."""
    return any(tok in ("Bm", "Bp") for tok in reactants + products)
