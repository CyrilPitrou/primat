# The graphical interface

`primat-gui` is a browser front end to the same solver the API and CLI use:
a parameter form on the left, results on the right, and every file a run
produces downloadable. It is the fastest way to see what a parameter does,
and the only interface that can build a network interactively.

```bash
pip install "primat[gui]"
primat-gui
```

A browser tab opens on a local Streamlit server. Nothing leaves your machine;
a public instance of the same app runs at
[primat.streamlit.app](https://primat.streamlit.app).

## The sidebar: setting up a run

The form covers the parameters a physics run usually varies, not all 96 of
them, in five expanders:

| Expander | Holds |
|---|---|
| **Cosmology** | `Omegabh2`, `DeltaNeff`, the neutrino degeneracy `munuOverTnu` |
| **Nuclear reactions** | `network`, `amax`, and the network builder |
| **Physics** | the correction flags, under "Weak rates" / "Plasma physics" / "Nuclear QED" |
| **Constants** | `GN`, `tau_n` and the measured constants |
| **Uncertainty** | the quick Monte Carlo, below |

Each widget matches its parameter's type — a toggle for a flag, a number box
for a number, a menu for a choice — and carries the same one-line explanation
{doc}`../parameters` gives. Anything not in the form is reachable from the API
or the CLI instead.

Two controls are the GUI's own rather than parameters:

- **Manage networks** (under "Nuclear reactions") — the gateway to every
  network action: build one, import one from a `.zip`, rename, remove, or
  select which to run. See {doc}`custom-networks`.
- **Uncertainty → Quick MC uncertainty** — after the main run, draw a small
  Monte Carlo (2–100 samples, default 30) over the nuclear rates and τ_n and
  show ±1σ beside each result. Raising the sample count reuses what has
  already been computed and solves only the difference. It is a noisy
  order-of-magnitude estimate, not a publication error bar — use
  {doc}`rate-variation-mc` for that.

The sidebar footer names the backend in use. To force one, set
`PRIMAT_GUI_BACKEND=c` or `=python` in the environment before launching.

**Run BBN** starts the solve. The button turns primary again whenever an edit
makes the displayed results out of date, so what you are looking at is always
either current or visibly marked stale.

## The four tabs

| Tab | What it shows |
|-----|---------------|
| **Reactions summary** | every reaction in the network as currently configured, with the rate table behind each — built from the sidebar as you edit it, before any run |
| **Final abundances** | the standard ratios (`YP`, `D/H`, `He3/H`, `Li7/H`, …), then every tracked nuclide's `A`, `Z` and final abundance per baryon |
| **Abundance evolution** | `A_i Y_i(t)` for every nuclide, interactive |
| **Output tables** | the run's output files, one download button each |

## Reproducing a GUI run elsewhere

Two downloads make a session's work portable:

- **Reproduce these results** (Final abundances tab) — a `.zip` holding a
  Python script *and* a `primat-c` INI file that reproduce exactly the numbers
  above them, with the backend pinned, plus the custom network if there is
  one.
- **Download network (zip)** (Reactions summary tab) — the network alone:
  every rate table verbatim on its original grid, re-importable through
  "Manage networks" and usable from the API as `custom_network=`. See
  {doc}`custom-networks` for what the archive contains and why it round-trips
  exactly.
