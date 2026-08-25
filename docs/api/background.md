# `primat.background` — expansion history

The cosmological background the nuclear network is integrated against:
`a ↔ t ↔ T`, the Hubble rate, the baryon density, and the n↔p weak rates as
the network sees them. `StandardBackground` is what a default run builds;
`CustomBackground` reads the history from a table instead.

`Background` is the base class to subclass for the `background=` hook — its
docstring lists the methods an implementation must provide and the ones it
may. {doc}`../howto/backgrounds` covers the three ways to change the
expansion history and when to reach for each.

```{eval-rst}
.. automodule:: primat.background
   :members:
   :undoc-members:
   :show-inheritance:
```
