# `primat.backend` — backend dispatch

`run_bbn(params)` is primat's main entry point: it picks the fast C engine
when the compiled extension imported (`HAS_C_BACKEND`) and the pure-Python
implementation otherwise, and returns the same dict either way.

**How closely the two agree.** Every physics formula, correction, clamp,
tolerance and default exists in both, and `tests/test_backend_parity.py` holds
them to `rel=5e-5` on D/H. The gap actually measured is 6.6e-06 relative on
`small` and 6.8e-06 on `large, amax=8`; `Neff` is identical to every digit.
That residual is not noise — it is the two backends integrating the
high-temperature era with different methods — so it does not shrink with
`numerical_precision`. It is also far larger than the same-backend regression
bound (±3e-9 absolute on D/H), which is the practical point: pick one backend
for a comparison and stay on it.

**The one feature gap** is `background=`, a user-supplied `Background`
instance, which cannot cross the C ABI: `force_backend="auto"` falls back to
Python for it and `force_backend="c"` raises. `extra_rho`, `decay_era`,
`custom_network` and the data overlays all work on both.

```{eval-rst}
.. automodule:: primat.backend
   :members:
   :undoc-members:
   :show-inheritance:
```
