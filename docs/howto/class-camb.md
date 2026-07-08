# Using primat with CLASS/CAMB tables, and in MCMC chains

## Driving primat from a CLASS/CAMB background

primat does not call CLASS or CAMB directly; instead, feed their
thermodynamics output through the {doc}`backgrounds`'s
`custom_background=` mechanism. Export `(T_γ, t, a)` from your CLASS/CAMB
run as a tab- or comma-delimited file (columns named `T` [MeV], `t` [s], `a`
— scale factor normalised to 1 today) and point `custom_background` at it:

```python
from primat.backend import run_bbn

result = run_bbn({
    "custom_background": "class_output.tsv",
    "network": "small",
    "Omegabh2": 0.02242,
})
```

See {doc}`backgrounds` for the exact table format, normalisation convention
(`a · T_γ → T_0CMB` as `T_γ → 0`), and the caveats around neutrino
temperatures (instantaneous-decoupling approximation in this mode).

## Using primat in MCMC chains (Cobaya)

For embedding BBN directly in MCMC analyses of CMB or other cosmological
data, use the Cobaya wrapper in the separate
[primat_tools](https://github.com/CyrilPitrou/primat_tools) repository (for
use with [Cobaya](https://cobaya.readthedocs.io)). It exposes `Omegabh2`,
`DeltaNeff`, and the nuclear-rate uncertainty parameters (see
{doc}`rate-variation-mc`) as Cobaya theory/likelihood inputs, and returns the
standard BBN observables (`YPBBN`, `DoH`, etc. — see {doc}`output`) for use
in a likelihood.

primat itself has no Cobaya dependency and none is implemented in this
package — keep the `run_bbn`/`PRIMAT.solve()` result-dict keys stable across
releases so the wrapper keeps working (see the primat_tools repo for its own
docs and installation instructions).
