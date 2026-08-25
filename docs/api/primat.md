# `primat` — top-level facade

`PRIMAT` is the object the pure-Python backend builds: a configuration, a
background and a nuclear network, with `solve()` returning the result dict
{doc}`../howto/output` describes. Most users reach it through
{doc}`backend`'s `run_bbn()` instead, and come here for the two things
`run_bbn` cannot express — the `background=` hook and direct access to a
solved instance's intermediate quantities (`get_quantity`).

```{eval-rst}
.. automodule:: primat.main
   :members:
   :undoc-members:
   :show-inheritance:
```
