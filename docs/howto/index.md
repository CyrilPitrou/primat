# How-to guides

Task-oriented recipes, one per job. {doc}`../tutorials/first-run` is the
guided introduction if you have not run primat yet.

## Setting up a run

- **{doc}`gui`** — the `primat-gui` browser app: the parameter form, the four
  result tabs, and the two downloads that make a session reproducible
  elsewhere.
- **{doc}`networks`** — pick a reaction network with `network`/`amax`, and
  what each choice costs in time and accuracy.
- **{doc}`custom-networks`** — build your own reaction set: interactively in
  the GUI, or permanently by adding a rate table to the source tree.

## Getting numbers out

- **{doc}`output`** — every key of the result dict, and the time-evolution TSV
  both backends write.
- **{doc}`rate-variation-mc`** — perturb one rate (`p_<reaction>`,
  `delta_<reaction>`), or run the Monte Carlo that turns rate uncertainties
  into an error bar on any observable.
- **{doc}`sensitivity`** — the ∂ln(observable)/∂ln(parameter) table, in one
  call, for every input at once.

## Non-standard physics

- **{doc}`backgrounds`** — three ways to drive the network with a
  non-standard expansion history: `extra_rho`, `custom_background`, and your
  own `Background` subclass.
- **{doc}`nevo-tables`** — swap the neutrino-decoupling tables for your own.
- **{doc}`class-camb`** — feed a CLASS or CAMB thermodynamics table into
  primat, and embed a run in an MCMC chain.

## Installation-level knobs

- **{doc}`data-overlays`** — redirect where data is read and written:
  `data_dir`, `user_nuclear_dir`, `cache_dir`. The one to reach for on a
  read-only install.
- **{doc}`weak-rate-cache`** — how the expensive n↔p rates are fingerprinted
  and cached, and how to control the recompute.

```{toctree}
:maxdepth: 1
:hidden:

gui
networks
custom-networks
output
rate-variation-mc
sensitivity
backgrounds
nevo-tables
class-camb
data-overlays
weak-rate-cache
```
