# Custom backgrounds and `extra_rho`

primat offers three mechanisms for driving the nuclear network with a
non-standard expansion history, from lightest-touch to most flexible:
`extra_rho` (add to the standard background), `custom_background` (replace
it with a table), and `background=` (replace it with your own class).

## `extra_rho` — add an energy-density component

The generic plug-in point. Pass a list of callables `rho(Tg) -> MeV^4` to
`PRIMAT(params, extra_rho=[...])`; each is summed into `rho_tot` by
`StandardBackground.Hubble` (`primat.background`) every time the Friedmann
equation is evaluated. This is the right tool for "add an extra
energy-density component to a standard run" — e.g. a constant dark-radiation
density:

```python
from primat import PRIMAT

PRIMAT({"network": "small"}, extra_rho=[lambda Tg: dRho])
```

Early Dark Energy (`cfg.fEDE > 0`) is itself implemented this way
(`StandardBackground._setup_EDE` appends to `self.extra_rho`) — read that
method as a worked example of a temperature-dependent, parametrised
contribution rather than a flat constant.

`extra_rho` is Python-only: any non-`None` value forces the Python backend
regardless of `force_backend="auto"` (and raises if `force_backend="c"`).

## `custom_background` — replace the expansion history with a table

For replacing the expansion history itself rather than adding to it. Set
`cfg.custom_background` to a path to a delimited file with columns `T`
[MeV], `t` [s], `a` (scale factor normalised to 1 today);
`CustomBackground` (`primat.background`) reads `T(t)`/`t`/`a(t)` directly
from the table instead of solving the entropy-conservation ODE, and falls
back to the instantaneous-decoupling approximation for neutrino
temperatures (`incomplete_decoupling`/`spectral_distortions` are forced to
`False` — NEVO tables are not loaded in this mode). `Neff` is estimated
from the Friedmann equation given the supplied `a(t)`.

This is the intended path for investigating modified expansion histories
(early dark energy, non-standard radiation content, time-varying
dark-energy equations of state), alternative neutrino physics (via a table
that encodes a different `T_ν(T_γ)` relationship implicitly through the
supplied `a(t)`), or non-standard temperature-time relationships (modified
gravity, other exotic scenarios).

**How it works:**

1. Create a tab- or comma-delimited text file with at least three columns
   named `T` [MeV], `t` [s], and `a` (scale factor, normalised to 1 today).
   The file must span the full BBN temperature range (typically from
   `T ~ 10 MeV` down to `T ~ 0.001 MeV`).
2. Set `cfg.custom_background` to the path of this file.
3. primat automatically forces `incomplete_decoupling=False` and
   `spectral_distortions=False`, then uses your table directly:
   - `T_of_t`/`t_of_T`/`a_of_T`/`T_of_a`/`t_of_a`/`a_of_t` are all read from
     the supplied table via linear interpolation.
   - Neutrino temperatures use the **instantaneous-decoupling**
     approximation (`T_ν = (4/11)^(1/3) * T_γ` for all flavours), computed
     from your supplied photon temperature `T`.
   - The n↔p weak rates are computed using these instantaneous-decoupling
     neutrino temperatures (no NEVO spectral distortions).
   - `Neff` is estimated at the end of BBN from the Friedmann equation
     `H² = 8πG/3 · ρ_tot` where `ρ_tot = ρ_plasma + ρ_ν`, and `ρ_ν` is
     inferred as the difference between `ρ_tot` (from `H` via your `a(t)`)
     and the known plasma density `ρ_plasma(T_γ)`.

```python
from primat.backend import run_bbn

result = run_bbn({
    "custom_background": "my_cosmology.tsv",
    "network": "small",
    "Omegabh2": 0.022425,
})

print(f"YP (BBN) = {result['YPBBN']:.8f}")
print(f"Neff     = {result['Neff']:.8f}")
```

**Important notes:**

- The scale factor `a` in your table **must** be normalised so that
  `a · T_γ → T_0CMB` as `T_γ → 0` (`a ≈ 1/T_γ` in radiation domination —
  the standard convention — so `a = 1` today when `T_γ = T_0CMB`).
- Rows may be in any order; primat sorts them internally by cosmic time `t`.
- Extra columns in your file are silently ignored.
- Mutually exclusive with {doc}`nevo-tables`'s `external_scale_factor=True`.
- Supported identically on both backends — see `tests/test_custom_background.py`
  for a round-trip example (write a reference background from a standard
  run, re-run through `custom_background`, check observables agree to
  <1e-5 relative).

## `background=` — inject a custom `Background` subclass

A `PRIMAT(background=<Background instance>)` hook for a fully custom
expansion history that needs more than `extra_rho` can express — subclass
`primat.background.Background` (whose docstring lists the compulsory
`T_of_t`/`t_of_T`/`rhoB_BBN`/`weak_nTOp_frwrd`/`weak_nTOp_bkwrd` and
optional methods) and pass an instance. `PRIMAT` takes `self.cfg`/
`self.plasma` from the supplied instance rather than building its own, so
build the instance with your own `PRIMATConfig`/`Plasma` first:

```python
from primat import Background
from primat.config import PRIMATConfig
from primat.plasma import Plasma

cfg = PRIMATConfig({"network": "small"})
plasma = Plasma(cfg)

class MyBackground(Background):
    ...  # implement T_of_t, t_of_T, rhoB_BBN, weak_nTOp_*

PRIMAT(background=MyBackground(cfg, plasma))
```

`None` (default) preserves the standard `cfg.custom_background`-based
dispatch. `background=` is mutually exclusive with `params`/`extra_rho`/
`cfg.custom_background` (only meaningful for the default dispatch);
supplying `background` together with either emits a warning, and the
supplied `background` instance wins. Like `extra_rho`, this is Python-only.

## Which mechanism to reach for

| Mechanism | Use for |
|-----------|---------|
| `extra_rho` | Adding a component on top of the standard background |
| `custom_background` | Replacing the expansion history with an externally-computed table |
| `background=` | Replacing the expansion history with your own `Background` subclass, built in code |
