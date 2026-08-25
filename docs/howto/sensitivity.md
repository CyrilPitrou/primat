# Sensitivity tables (∂ln O / ∂ln p)

Referees of BBN papers routinely ask *"how sensitive is each abundance to each
input?"*. primat answers this in one call with
{func}`primat.sensitivity.sensitivity_table`, which finite-differences full BBN
solves to build the **logarithmic-sensitivity matrix**

$$
S(O, p) \equiv \frac{\partial \ln O}{\partial \ln p}
\approx \frac{\ln O\!\left(p(1+\delta)\right) - \ln O\!\left(p(1-\delta)\right)}
             {2\,\ln(1+\delta)} .
$$

$S(O,p)$ is dimensionless: *"if parameter $p$ rises by 1 %, observable $O$ rises
by $S$ %."* The symmetric difference is accurate to $O(\delta^2)$, so the
default 1 % step (`rel_step=0.01`) already gives ~4 correct digits.

## Quick start

```python
from primat.sensitivity import sensitivity_table

tab = sensitivity_table(
    params={"network": "small", "Omegabh2": 0.02242},
    observables=["YPBBN", "DoH", "He3oH", "Li7oH"],
    targets=[
        "n_p__d_g",   # a nuclear reaction rate  (auto-detected)
        "tau_n",      # neutron lifetime          (multiplicative)
        "Omegabh2",   # baryon density            (multiplicative)
    ],
)

print(tab.to_markdown())     # paste-ready GitHub table
df = tab.to_dataframe()      # pandas view, indexed by parameter / observable
tab.values                   # raw (n_targets, n_observables) numpy array
tab.fiducial                 # unperturbed observable values (shared central solve)
```

The call runs one shared central solve plus two bracketing solves per target
(`2 * len(targets) + 1` solves total) on the fast C backend, so a full
12-rate + 4-parameter table takes a few tens of seconds.

## Three flavours of parameter

Each entry of `targets` is either a bare string (auto-classified) or an explicit
{class}`~primat.sensitivity.SensTarget` for full control.

| Flavour | How it is varied | How to pass it |
|---------|------------------|----------------|
| **Nuclear rate** (`n_p__d_g`, `d_p__He3_g`, …) | `delta_<rxn>=±δ`, i.e. `rate = median·(1±δ)` | bare string, or `SensTarget("n_p__d_g")` |
| **Multiplicative parameter** (`tau_n`, `GN`, `Omegabh2`, …) | scaled by `(1±δ)` about its fiducial value | bare string, or `SensTarget("tau_n", label=r"$\tau_n$")` |
| **Additive parameter** (`DeltaNeff`, fiducial 0) | varied by an absolute `±step` about the base value, normalised by the `ref` parameter it offsets | `SensTarget("DeltaNeff", kind="additive", step=0.1, ref="Neff")` |

A bare string is auto-classified: it becomes a **rate** target if it names a
reaction in the config's rate table, otherwise a **multiplicative** target.
Because a proportional step `p·(1±δ)` can never move a parameter whose fiducial
value is 0, `DeltaNeff` (its Standard-Model value is 0) must be given as an
explicit `additive` target — otherwise `sensitivity_table` raises a `ValueError`
telling you exactly that.

### Keeping an additive row an elasticity: `ref`

Every cell of the table means the same thing — $\partial\ln O/\partial\ln p$,
"a 1 % variation of $p$ appears as this many % of $O$" — and that is what makes
the rows comparable. An additive knob has no $\ln p$ to differentiate against at
$p = 0$, but it is normally just an *offset* on a parameter that does:
$N_{\rm eff} = N_{\rm eff}^{\rm SM} + \Delta N_{\rm eff}$. Name that parameter's
fiducial with `ref` and the row is $\partial\ln O/\partial\ln N_{\rm eff}$ like
any other:

```python
SensTarget("DeltaNeff", kind="additive", step=0.1, ref="Neff",
           label=r"$N_{\rm eff}$")
```

`ref` accepts either a string naming a key of the result dict — read from *this*
run's central solve, so you are not hard-coding 3.044 — or an explicit number.
Keep `step` small (0.1, not 1.0): `step=1.0` is a ±33 % excursion in
$N_{\rm eff}$, whose secant differs from the derivative by about 1 %.

Omitting `ref` is allowed and falls back to dividing by the linear separation
`2·step`, giving $\partial\ln O/\partial p$ **per unit** of $p$. That is a
different quantity — not an elasticity — so such a row must not be ranked
against the others; say so in the caption if you publish one.

### Why the `delta_<rxn>` mechanism for rates?

Reaction rates are varied through primat's deterministic rescaling knob (see
[Rate variation and Monte-Carlo uncertainty](rate-variation-mc.md)), *not* the
log-normal `p_<rxn>` knob: `delta_<rxn>=δ` multiplies the median rate by exactly
`(1+δ)`, independent of the rate's tabulated uncertainty, which is what a clean
$\partial\ln O/\partial\ln(\text{rate})$ requires.

## Custom labels, steps, and observables

```python
from primat.sensitivity import sensitivity_table, SensTarget

tab = sensitivity_table(
    params={"network": "small"},
    observables=["YPBBN", "DoH"],
    obs_labels=[r"$Y_P$", r"D/H"],           # pretty column headers
    rel_step=0.02,                            # 2% step instead of the 1% default
    targets=[
        SensTarget("n_p__d_g", label=r"p+n→d+γ"),
        SensTarget("tau_n",    label=r"$\tau_n$"),
        SensTarget("DeltaNeff", label=r"$N_{\rm eff}$",
                   kind="additive", step=0.1, ref="Neff"),  # dlnO/dlnNeff
    ],
)
```

Fiducial values of multiplicative targets are read from the config built from
your `params`, so overriding e.g. `tau_n` in `params` differentiates about that
shifted point. `force_backend={"auto","c","python"}` selects the solver
(default prefers the C backend).

## Worked example notebook

{doc}`../tutorials/Sensitivity` (`notebooks/Sensitivity.ipynb` in the
repository) is a full demo: it builds the 12-reaction + 4-parameter target
list, calls `sensitivity_table` once, prints the tables, and renders them as a
sensitivity heat-map.

:::{seealso}
- [Rate variation and Monte-Carlo uncertainty](rate-variation-mc.md) — the
  `delta_<rxn>` / `p_<rxn>` knobs and full covariance propagation.
- {func}`primat.sensitivity.sensitivity_table` — full API reference.
:::
