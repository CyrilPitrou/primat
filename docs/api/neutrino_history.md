# `primat.neutrino_history` — neutrino decoupling

Where `T_ν(T_γ)` comes from: the tabulated non-instantaneous decoupling
history (`NEVOTable`), the analytic instantaneous-decoupling limit
(`InstantaneousDecoupling`), and the optional `AnalyticDistortion` wrapper.
`make_neutrino_history(cfg, plasma)` is the dispatch, and the place a new
variant plugs in — see {doc}`../extending`.

To swap the *data* under the existing classes rather than write a new one,
use the `nevo_file`/`nevo_spectral_file`/`nevo_grid_file`/`nevo_file_prefix`
parameters ({doc}`../howto/nevo-tables`).

```{eval-rst}
.. automodule:: primat.neutrino_history
   :members:
   :undoc-members:
   :show-inheritance:
```
