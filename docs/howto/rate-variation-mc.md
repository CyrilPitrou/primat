# Rate variation and Monte-Carlo uncertainty

:::{note}
*(stub — FABLEADVICE O-3)* Migrate from the README's rate-variation section and
`runfiles/primat_mc.py`: the per-reaction `p_<reaction>` variation parameters
(`median * exp(p * expsigma)`), `mc_uncertainty` / `run_mc`, the
`MCResult.corr()` / `cov()` covariance and correlation outputs, and parallel
execution via the `mc` / `recommended` extras (joblib).
:::

Each reaction rate has a `p_<reaction>` knob (e.g. `p_n_p__d_g`,
`p_Li7_p__a_a`); setting it to a non-zero float samples the rate at
`median * exp(p * expsigma)`, the hook Monte-Carlo uncertainty propagation
drives.
