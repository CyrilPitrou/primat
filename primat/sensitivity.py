# -*- coding: utf-8 -*-
"""
primat.sensitivity — logarithmic sensitivity of BBN observables to parameters.

This module promotes the ad-hoc finite-difference loop that used to live in
``notebooks/Sensitivity.ipynb`` into a reusable public API. Referees of BBN
papers routinely ask for a *sensitivity table*

.. math::

    S(O, p) \\equiv \\frac{\\partial \\ln O}{\\partial \\ln p},

the dimensionless "if parameter :math:`p` rises by 1 %, observable :math:`O`
rises by :math:`S` %". :func:`sensitivity_table` computes the whole matrix by
symmetric finite-differencing full BBN solves and returns a
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
* **Additive parameters** whose fiducial value is zero, so a multiplicative
  step is degenerate (``DeltaNeff`` at its SM value 0): varied by an absolute
  ``±step`` about the base value. The denominator is user-controlled (default
  ``2 ln(1+rel_step)``) because the "natural" normalisation of an additive knob
  is a modelling choice, not a mathematical given.

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
...         SensTarget("DeltaNeff", kind="additive", step=1.0),
...     ],
... )
>>> print(tab.to_markdown())          # doctest: +SKIP
| Parameter | $Y_P$ | D/H |
| --- | --- | --- |
| n_p__d_g  | +0.024 | +0.207 |
| ...
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
        denom: override the finite-difference denominator. ``None`` (default)
            uses ``2 ln(1+rel_step)``, which makes ``rate``/``param`` rows exact
            logarithmic derivatives. Handy for additive rows where you want a
            different normalisation.

    Example:
        >>> SensTarget("Omegabh2", label=r"$\\Omega_b h^2$")
        SensTarget(param='Omegabh2', ...)
        >>> SensTarget("DeltaNeff", kind="additive", step=1.0)   # per unit ΔNeff
        SensTarget(param='DeltaNeff', ...)
    """

    param: str
    label: str | None = None
    kind: str = "auto"
    step: float | None = None
    denom: float | None = None

    def display_label(self) -> str:
        """Row label to print — the explicit ``label`` or the parameter name."""
        return self.label if self.label is not None else self.param

    def resolve(self, cfg: Any, rel_step: float) -> tuple[dict, dict, float]:
        """Return ``(plus_params, minus_params, denom)`` for this target.

        ``cfg`` is the *effective* :class:`~primat.config.PRIMATConfig` built
        from the caller's ``params`` (so fiducial values honour user overrides,
        e.g. a non-default ``tau_n``). ``rel_step`` is the fractional step
        :math:`\\delta`. The two returned dicts are merged onto the base params
        for the ``+`` and ``-`` runs; ``denom`` divides ``ln O+ - ln O-``.

        Raises:
            ValueError: for a multiplicative target whose fiducial value is 0
                (a proportional step can never move it), or an unknown ``kind``.
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
            # Absolute step about the base value (used for knobs whose fiducial
            # is 0, e.g. DeltaNeff). The denominator defaults to the same
            # log-derivative normalisation so an additive row is directly
            # comparable to the multiplicative ones at the chosen rel_step.
            s = self.step if self.step is not None else d
            base = float(getattr(cfg, self.param))
            plus = {self.param: base + s}
            minus = {self.param: base - s}
            return plus, minus, (self.denom if self.denom is not None else default_denom)

        raise ValueError(
            f"unknown SensTarget kind {self.kind!r}; expected one of "
            f"'auto', 'rate', 'param', 'additive'"
        )


@dataclass
class SensitivityTable:
    """Result of :func:`sensitivity_table`: a (targets × observables) matrix.

    Holds the logarithmic-sensitivity matrix plus enough context to render it.
    ``values[i, j]`` is :math:`\\partial \\ln O_j / \\partial \\ln p_i` (or the
    additive analogue for additive targets).

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
    :math:`\\partial\\ln O/\\partial\\ln p`. A single *shared* central solve at
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
            control (additive knobs, custom labels, custom denominators).
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
        plus_over, minus_over, denom = t.resolve(cfg, rel_step)
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
