# -*- coding: utf-8 -*-
"""
primat.sensitivity — logarithmic sensitivity of BBN observables to parameters.

Referees of BBN papers routinely ask for a *sensitivity table*

.. math::

    S(O, p) \\equiv \\frac{\\partial \\ln O}{\\partial \\ln p},

the dimensionless *elasticity*: :math:`S = 1` means a 1 % variation of
:math:`p` shows up as a 1 % variation of :math:`O`, :math:`S = -0.2` means a
1 % rise in :math:`p` lowers :math:`O` by 0.2 %. **Every** row of the table is
that same ``d ln O / d ln p``, whatever the parameter's flavour — that is what
makes the rows comparable with one another. :func:`sensitivity_table` computes
the whole matrix by symmetric finite-differencing full BBN solves and returns a
:class:`SensitivityTable` dataclass with ``to_markdown()`` / ``to_dataframe()``
views ready to paste into a paper or notebook.

Why symmetric finite differences? For a fractional step :math:`\\delta` the
central estimate

.. math::

    S(O, p) \\approx
      \\frac{\\ln O\\!\\left(p(1+\\delta)\\right) - \\ln O\\!\\left(p(1-\\delta)\\right)}
           {2\\,\\ln(1+\\delta)}

is accurate to :math:`O(\\delta^2)` (the linear error term cancels), so the
default 1 % step (``rel_step=0.01``) already gives ~4 correct digits without
paying for tiny-step round-off. The denominator ``2 ln(1+δ)`` — rather than the
naive ``2δ`` — makes the result an exact *logarithmic* derivative: the two runs
sit at :math:`\\ln p \\pm \\ln(1+\\delta)` in log-parameter space.

Three flavours of parameter need three variation recipes, all expressed through
the :class:`SensTarget` helper (a plain string is auto-classified):

* **Nuclear reaction rates** (e.g. ``"n_p__d_g"``): varied through primat's
  ``rescale_nuclear_rates`` + ``delta_<rxn>`` mechanism, which multiplies the
  median rate by ``(1 + delta)`` — see the *Rate variation* how-to. The
  resulting number is a genuine :math:`\\partial\\ln O/\\partial\\ln(\\text{rate})`.
* **Multiplicative physical parameters** with a non-zero fiducial value
  (``tau_n``, ``GN``, ``Omegabh2``, ...): scaled by ``(1 ± δ)`` about their
  effective value, giving :math:`\\partial\\ln O/\\partial\\ln p` directly.
* **Additive parameters** whose own fiducial value is zero, so a multiplicative
  step is degenerate (``DeltaNeff`` at its SM value 0): varied by an absolute
  ``±step``. Such a knob is normally an *offset* on a physical parameter that is
  perfectly non-zero — ``DeltaNeff`` displaces
  :math:`N_{\\rm eff} = N_{\\rm eff}^{\\rm SM} + \\Delta N_{\\rm eff}`, fiducial
  :math:`\\simeq 3.044` — so name that parameter's fiducial with ``ref`` and the
  row comes out as :math:`\\partial\\ln O/\\partial\\ln N_{\\rm eff}`, the same
  elasticity as every other row (``ref="Neff"`` reads the run's own central
  value). Without ``ref`` the fallback denominator is the linear separation
  ``2*step``, giving :math:`\\partial\\ln O/\\partial p` *per unit of* ``p`` —
  useful, but not an elasticity and therefore not comparable with the other
  rows, since there is no :math:`\\ln p` at :math:`p = 0`.

Example
-------
>>> from primat.sensitivity import sensitivity_table, SensTarget
>>> tab = sensitivity_table(
...     params={"network": "small"},
...     observables=["YPBBN", "DoH"],
...     targets=[
...         "n_p__d_g",                       # a nuclear rate  (auto -> rate)
...         "tau_n",                          # neutron lifetime (auto -> mult.)
...         "Omegabh2",                       # baryon density   (auto -> mult.)
...         # d ln O / d ln Neff: an absolute +-0.1 step on the DeltaNeff
...         # offset, normalised by the run's own central Neff (~3.044).
...         SensTarget("DeltaNeff", kind="additive", step=0.1, ref="Neff",
...                    label=r"$N_{\\rm eff}$"),
...     ],
... )
>>> print(tab.to_markdown())          # doctest: +SKIP
| Parameter | $Y_P$ | D/H |
| --- | --- | --- |
| n_p__d_g | +0.0045 | -0.2024 |
| tau_n | +0.7355 | +0.4230 |
| Omegabh2 | +0.0391 | -1.6501 |
| $N_{\\rm eff}$ | +0.1647 | +0.4104 |

Read every row the same way — the last one as "1 % more radiation density,
measured as :math:`N_{\\rm eff}`, gives 0.16 % more helium and 0.41 % more
deuterium". The familiar :math:`{\\rm D/H} \\propto (\\Omega_b h^2)^{-1.6}`
scaling is the ``Omegabh2`` row. The ``0.1`` step keeps that row a derivative:
``step=1.0`` would be a :math:`\\pm 33\\,\\%` excursion in
:math:`N_{\\rm eff}`, far enough for the secant to depart from the tangent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

# Human-readable column headers for the observables primat's result dict
# exposes. Anything not listed falls back to its raw result-dict key, so this
# is purely cosmetic (it never changes the computed numbers).
DEFAULT_OBS_LABELS: dict[str, str] = {
    "YPBBN": "$Y_P$",
    "YPCMB": "$Y_P^{\\rm CMB}$",
    "He4oH": "$^4$He/H",
    "DoH": "D/H",
    "He3oH": "$^3$He/H",
    "Li7oH": "$^7$Li/H",
    "Neff": "$N_{\\rm eff}$",
}


@dataclass(frozen=True)
class SensTarget:
    """One row of a sensitivity table: *what parameter to vary and how*.

    A :class:`SensTarget` knows how to turn a fractional step ``rel_step`` into
    the two bracketing parameter overrides (``plus``/``minus``) fed to
    :func:`primat.backend.run_bbn`, plus the finite-difference denominator. Most
    callers never build one explicitly: passing a bare string to
    :func:`sensitivity_table` constructs ``SensTarget(param)`` with
    ``kind="auto"``, which classifies the string as a reaction rate (if it names
    a reaction in the effective config's ``p_rxn`` table) or a multiplicative
    physical parameter otherwise.

    Attributes:
        param: config-parameter name (e.g. ``"tau_n"``, ``"Omegabh2"``,
            ``"DeltaNeff"``) *or* a nuclear reaction name (e.g. ``"n_p__d_g"``).
        label: display name for this row in tables. Defaults to ``param``.
        kind: ``"auto"`` (default; classify at resolve time), ``"rate"``
            (force the ``delta_<rxn>`` rescaling mechanism), ``"param"`` (force
            multiplicative ``p(1±δ)``) or ``"additive"`` (absolute ``±step``).
        step: for ``kind="additive"`` only, the absolute half-step (e.g.
            ``1.0`` for ``DeltaNeff``). Ignored otherwise; defaults to
            ``rel_step`` when ``None``.
        ref: for ``kind="additive"`` only — the fiducial value of the *physical*
            parameter the additive step displaces, which turns the row into a
            true elasticity ``dln O/dln(that parameter)``, directly comparable
            with every other row. ``DeltaNeff`` displaces
            ``Neff = Neff_SM + DeltaNeff``, so ``ref="Neff"`` (a string naming a
            key of the central solve's result dict, read from *this* run rather
            than a hard-coded 3.044) or ``ref=3.044`` (an explicit number) both
            give ``dln O/dln Neff``. ``None`` (default) falls back to the
            per-unit form below. Setting it on a non-additive target raises,
            since those already differentiate about their own fiducial.
        denom: override the finite-difference denominator outright. ``None``
            (default) uses ``2 ln(1+rel_step)`` for ``rate``/``param`` rows and,
            for ``additive`` rows, ``2*step/ref`` when ``ref`` is given, else
            ``2*step`` — which makes the row a semi-logarithmic ``dln O/dp``
            *per unit of p*, NOT an elasticity (there is no ``ln p`` at p = 0)
            and therefore not comparable with the other rows. Prefer ``ref``
            over a hand-computed ``denom``.

    Example:
        >>> SensTarget("Omegabh2", label=r"$\\Omega_b h^2$")
        SensTarget(param='Omegabh2', ...)
        >>> # dln O / dln Neff -- the form to prefer for DeltaNeff
        >>> SensTarget("DeltaNeff", kind="additive", step=0.1, ref="Neff")
        SensTarget(param='DeltaNeff', ...)
        >>> # dln O / dDeltaNeff, per unit: not an elasticity, not comparable
        >>> SensTarget("DeltaNeff", kind="additive", step=1.0)
        SensTarget(param='DeltaNeff', ...)
    """

    param: str
    label: str | None = None
    kind: str = "auto"
    step: float | None = None
    denom: float | None = None
    ref: float | str | None = None

    def display_label(self) -> str:
        """Row label to print — the explicit ``label`` or the parameter name."""
        return self.label if self.label is not None else self.param

    def resolve(self, cfg: Any, rel_step: float,
                central: dict | None = None) -> tuple[dict, dict, float]:
        """Return ``(plus_params, minus_params, denom)`` for this target.

        ``cfg`` is the *effective* :class:`~primat.config.PRIMATConfig` built
        from the caller's ``params`` (so fiducial values honour user overrides,
        e.g. a non-default ``tau_n``). ``rel_step`` is the fractional step
        :math:`\\delta`. ``central`` is the shared central solve's result dict,
        needed only to resolve a string :attr:`ref` (e.g. ``ref="Neff"`` reads
        this run's own central ``Neff``). The two returned dicts are merged onto
        the base params for the ``+`` and ``-`` runs; ``denom`` divides
        ``ln O+ - ln O-``.

        Raises:
            ValueError: for a multiplicative target whose fiducial value is 0
                (a proportional step can never move it), an unknown ``kind``,
                a non-positive ``ref`` (there is no ``ln p`` to take), or a
                ``ref`` on a non-additive target (where it would be a silent
                no-op -- the other kinds already have their own fiducial).
            KeyError: for a string ``ref`` naming no key of ``central``.
        """
        d = rel_step
        # Default log-derivative denominator: the two runs sit at
        # ln(p) ± ln(1+δ), so dividing by 2 ln(1+δ) yields d ln O / d ln p.
        default_denom = 2.0 * math.log1p(d)

        kind = self.kind
        if kind == "auto":
            # A reaction name lives in the config's per-reaction rate table;
            # anything else is treated as a scalar physical parameter.
            kind = "rate" if self.param in cfg.p_rxn else "param"

        if self.ref is not None and kind != "additive":
            raise ValueError(
                f"ref={self.ref!r} is only meaningful for kind='additive' "
                f"(target {self.param!r} resolved to kind={kind!r}, which is "
                "already differentiated about its own non-zero fiducial); drop "
                "ref, or pass kind='additive' explicitly."
            )

        if kind == "rate":
            # primat's deterministic rate-rescaling knob: with
            # rescale_nuclear_rates=True and p_<rxn>=0 (default), the rate
            # becomes median*(1 + delta_<rxn>). So delta=±δ brackets the rate
            # multiplicatively and the log-derivative denominator applies.
            key = f"delta_{self.param}"
            plus = {"rescale_nuclear_rates": True, key: +d}
            minus = {"rescale_nuclear_rates": True, key: -d}
            return plus, minus, (self.denom if self.denom is not None else default_denom)

        if kind == "param":
            fid = float(getattr(cfg, self.param))
            if fid == 0.0:
                raise ValueError(
                    f"parameter {self.param!r} has fiducial value 0, so a "
                    f"multiplicative step p*(1±δ) cannot move it; pass "
                    f"SensTarget({self.param!r}, kind='additive', step=...) "
                    f"to vary it by an absolute amount instead."
                )
            plus = {self.param: fid * (1.0 + d)}
            minus = {self.param: fid * (1.0 - d)}
            return plus, minus, (self.denom if self.denom is not None else default_denom)

        if kind == "additive":
            # Absolute step about the base value, for knobs whose own fiducial
            # is 0 so that p*(1±delta) cannot move them (DeltaNeff).
            s = self.step if self.step is not None else d
            base = float(getattr(cfg, self.param))
            plus = {self.param: base + s}
            minus = {self.param: base - s}
            if self.denom is not None:
                return plus, minus, self.denom
            # Denominator = the two runs' separation in the space we are
            # differentiating against.
            #
            # With `ref`: the knob is an OFFSET on a physical parameter P whose
            # fiducial is ref (e.g. Neff = Neff_SM + DeltaNeff, ref ~ 3.044), so
            # the runs sit at ln P = ln(ref) ± ... and the separation is
            # d ln P = 2*step/ref. The cell is then dln(O)/dln(P) -- the SAME
            # dimensionless elasticity as every other row ("+1 means 1 % in P
            # gives 1 % in O"), which is the documented meaning of the whole
            # table. This is the recommended form for an additive knob.
            #
            # Without `ref`: fall back to the linear separation 2*step, giving
            # dln(O)/dp per unit of p. Not an elasticity (there is no ln p at
            # p = 0), so such a row is NOT comparable with the others -- hence
            # the strong preference for ref.
            #
            # Either way the denominator must NOT be the multiplicative rows'
            # 2 ln(1+rel_step): rel_step plays no part in an absolute ±step
            # perturbation, and using it made the reported number scale with a
            # parameter that changes nothing (at DeltaNeff/step=1.0/rel_step=0.01
            # it inflated the cell by 1/ln(1.01) ~ 100.5).
            if self.ref is None:
                return plus, minus, 2.0 * s
            ref = self._resolve_ref(central)
            return plus, minus, 2.0 * s / ref

        raise ValueError(
            f"unknown SensTarget kind {self.kind!r}; expected one of "
            f"'auto', 'rate', 'param', 'additive'"
        )

    def _resolve_ref(self, central: dict | None) -> float:
        """Numeric value of :attr:`ref`: itself when a float, else the named key
        of the central solve's result dict (``ref="Neff"`` -> this run's own
        central ``Neff``, so the elasticity is taken about the actual fiducial
        rather than a hard-coded 3.044).

        Raises:
            KeyError: a string ref naming no result key (or no central solve
                available, which cannot happen through
                :func:`sensitivity_table`).
            ValueError: a non-positive ref, whose logarithm does not exist.
        """
        if isinstance(self.ref, str):
            if central is None or self.ref not in central:
                available = sorted(central) if central else []
                raise KeyError(
                    f"ref={self.ref!r} names no key of the central result dict "
                    f"(available: {available}); pass a number instead."
                )
            ref = float(central[self.ref])
        else:
            ref = float(self.ref)
        if ref <= 0.0:
            raise ValueError(
                f"ref must be strictly positive to take d ln(ref) (got {ref!r} "
                f"for target {self.param!r})."
            )
        return ref


@dataclass
class SensitivityTable:
    """Result of :func:`sensitivity_table`: a (targets × observables) matrix.

    Holds the logarithmic-sensitivity matrix plus enough context to render it.
    ``values[i, j]`` is the dimensionless elasticity
    :math:`\\partial \\ln O_j / \\partial \\ln p_i` — "1 % in :math:`p_i` gives
    that many % in :math:`O_j`" — for ``rate`` and ``param`` rows, and for
    ``additive`` rows given a ``ref`` (where the logarithm is taken of the
    physical parameter the offset displaces, e.g. :math:`N_{\\rm eff}` for
    ``DeltaNeff``). The one exception is an ``additive`` row *without* ``ref``:
    that cell is :math:`\\partial \\ln O_j / \\partial p_i` per unit of
    :math:`p_i`, which is not an elasticity and must not be ranked against the
    others — see :class:`SensTarget`.

    Attributes:
        row_labels: display label of each varied parameter (table rows).
        observables: result-dict keys of each observable (table columns).
        obs_labels: pretty column headers, aligned with ``observables``.
        values: ``(n_targets, n_observables)`` float array of sensitivities.
        fiducial: observable-key → fiducial (unperturbed) value, from the single
            shared central solve.
        rel_step: the fractional step used for the finite differences.

    Example:
        >>> tab.to_dataframe()                       # doctest: +SKIP
        >>> print(tab.to_markdown())                 # doctest: +SKIP
        >>> tab.values[0, 1]                         # dln(D/H)/dln(first param)
    """

    row_labels: list[str]
    observables: list[str]
    obs_labels: list[str]
    values: np.ndarray
    fiducial: dict[str, float]
    rel_step: float = 0.01

    def to_dataframe(self):
        """Return the sensitivity matrix as a pandas ``DataFrame``.

        Rows are indexed by parameter label, columns by observable label.
        Requires pandas (an optional dependency); raised as an ``ImportError``
        with an actionable message if it is missing.

        Example:
            >>> df = tab.to_dataframe()              # doctest: +SKIP
            >>> df.loc[r"$\\tau_n$", "D/H"]          # doctest: +SKIP
        """
        try:
            import pandas as pd
        except ImportError as exc:  # pragma: no cover - trivial guard
            raise ImportError(
                "SensitivityTable.to_dataframe() needs pandas; install it with "
                "`pip install pandas` or use to_markdown() for a text table."
            ) from exc
        return pd.DataFrame(self.values, index=self.row_labels, columns=self.obs_labels)

    def to_markdown(self, fmt: str = "+.4f") -> str:
        """Render the table as a GitHub-flavoured Markdown string.

        Hand-rolled (no pandas/tabulate dependency) so it always works. ``fmt``
        is the per-cell format spec (default ``"+.4f"`` — signed, 4 decimals,
        which resolves the ~1e-2–1e-3 sensitivities BBN cares about).

        Example:
            >>> print(tab.to_markdown())             # doctest: +SKIP
            | Parameter | $Y_P$ | D/H |
            | --- | --- | --- |
            | $\\tau_n$ | +0.7290 | +0.4130 |
        """
        header = "| Parameter | " + " | ".join(self.obs_labels) + " |"
        sep = "| --- | " + " | ".join("---" for _ in self.obs_labels) + " |"
        lines = [header, sep]
        for i, label in enumerate(self.row_labels):
            cells = " | ".join(format(self.values[i, j], fmt)
                               for j in range(len(self.observables)))
            lines.append(f"| {label} | {cells} |")
        return "\n".join(lines)

    def __repr__(self) -> str:  # concise, avoids dumping the whole ndarray
        return (f"SensitivityTable({len(self.row_labels)} params × "
                f"{len(self.observables)} observables, rel_step={self.rel_step})")


def sensitivity_table(
    params: dict[str, Any] | None,
    observables: list[str],
    targets: list[Any],
    rel_step: float = 0.01,
    *,
    obs_labels: list[str] | None = None,
    force_backend: str | None = None,
    progress: bool = False,
) -> SensitivityTable:
    """Compute the logarithmic-sensitivity matrix of BBN observables.

    For each target parameter this runs two full BBN solves bracketing the
    fiducial point (``p(1±δ)`` for multiplicative targets, ``±delta`` for rate
    targets, ``base±step`` for additive targets) and forms the symmetric
    finite-difference logarithmic derivative
    :math:`\\partial\\ln O/\\partial\\ln p` — including for additive targets
    carrying a ``ref``, whose logarithm is taken of the physical parameter the
    offset displaces; only a ``ref``-less additive row is instead
    :math:`\\partial\\ln O/\\partial p` per unit of ``p`` (see
    :class:`SensTarget`). A single *shared* central solve at
    the fiducial parameters is run once and stored in the result for reference
    (the symmetric estimate itself uses only the two bracketing runs).

    Cost is ``2 * len(targets) + 1`` solves; with the default C backend a
    ``small``-network solve is ~sub-second, so a full 12-rate + 4-parameter
    table is a few tens of seconds.

    Args:
        params: base ("fiducial") parameter overrides, exactly as accepted by
            :func:`primat.backend.run_bbn` / ``PRIMAT(params=...)``. ``None`` ==
            all defaults. Fiducial values of multiplicative targets are read
            from the config built from *these* params, so overriding e.g.
            ``tau_n`` here shifts the point about which it is differentiated.
        observables: result-dict keys to differentiate, e.g.
            ``["YPBBN", "DoH", "He3oH", "Li7oH"]``.
        targets: parameters to vary. Each item is either a plain string
            (auto-classified into a rate or multiplicative parameter — see
            :class:`SensTarget`) or an explicit :class:`SensTarget` for full
            control (additive knobs and their ``ref``, custom labels, custom
            denominators).
        rel_step: fractional finite-difference step :math:`\\delta`
            (default 0.01 = 1 %). Small enough for ~4-digit accuracy, large
            enough to stay clear of solver round-off.
        obs_labels: optional pretty column headers aligned with ``observables``;
            defaults to :data:`DEFAULT_OBS_LABELS` (falling back to the raw key).
        force_backend: forwarded to :func:`primat.backend.run_bbn`
            (``None``/``"auto"`` prefers the fast C backend; ``"c"``/``"python"``
            force one).
        progress: forwarded to ``run_bbn`` — leave ``False`` to silence the
            per-solve phase markers during the (many) sensitivity runs.

    Returns:
        SensitivityTable: with ``.values`` the ``(len(targets), len(observables))``
        matrix, ``.to_markdown()`` / ``.to_dataframe()`` views, and ``.fiducial``
        the shared central observable values.

    Example:
        >>> tab = sensitivity_table(
        ...     {"network": "small"},
        ...     observables=["YPBBN", "DoH"],
        ...     targets=["tau_n", "Omegabh2", "n_p__d_g"],
        ... )
        >>> tab.values.shape
        (3, 2)
    """
    from .backend import run_bbn
    from .config import PRIMATConfig

    base = dict(params or {})
    # Effective config: fiducial values of multiplicative targets and the
    # p_rxn membership test (used to auto-classify string targets) both come
    # from here, so any user override in `params` is respected.
    cfg = PRIMATConfig(base)

    # Column headers: caller-supplied, else the pretty default, else raw key.
    if obs_labels is None:
        obs_labels = [DEFAULT_OBS_LABELS.get(o, o) for o in observables]
    elif len(obs_labels) != len(observables):
        raise ValueError(
            f"obs_labels has {len(obs_labels)} entries but there are "
            f"{len(observables)} observables"
        )

    def _solve(overrides: dict) -> dict:
        # One BBN solve with the given overrides merged onto the fiducial
        # params; verbose/progress forced off so a big table stays quiet.
        merged = {**base, **overrides, "verbose": False}
        return run_bbn(merged, force_backend=force_backend, progress=progress)

    # Single shared central solve at the fiducial point, reused by every
    # perturbed sensitivity solve. Used for the reported fiducial values.
    r0 = _solve({})
    missing = [o for o in observables if o not in r0]
    if missing:
        raise KeyError(
            f"observable(s) {missing} are not in the result dict; available "
            f"keys include e.g. YPBBN, DoH, He3oH, Li7oH, Neff"
        )
    fiducial = {o: float(r0[o]) for o in observables}

    # Normalise every target to a SensTarget (bare strings -> kind='auto').
    norm_targets = [t if isinstance(t, SensTarget) else SensTarget(t) for t in targets]

    values = np.zeros((len(norm_targets), len(observables)))
    row_labels: list[str] = []
    for i, t in enumerate(norm_targets):
        row_labels.append(t.display_label())
        # r0 is passed so a string `ref` (e.g. SensTarget(..., ref="Neff"))
        # resolves against THIS run's central observables rather than a
        # hard-coded fiducial.
        plus_over, minus_over, denom = t.resolve(cfg, rel_step, central=r0)
        rp = _solve(plus_over)
        rm = _solve(minus_over)
        for j, o in enumerate(observables):
            # Symmetric logarithmic finite difference; the O(δ²) error term of
            # the naive forward difference cancels between the ± runs.
            values[i, j] = (math.log(rp[o]) - math.log(rm[o])) / denom

    return SensitivityTable(
        row_labels=row_labels,
        observables=list(observables),
        obs_labels=list(obs_labels),
        values=values,
        fiducial=fiducial,
        rel_step=rel_step,
    )
