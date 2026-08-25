# Extending primat with new physics

A short guide to the three supported extension points. For the full physics
formalism (rate evaluation, background equations, weak-rate corrections),
the primary reference is the Physics Reports paper linked from
{doc}`physics` (Pitrou, Coc, Uzan, Vangioni, *Physics Reports* 04 (2018)
005); see also the other PDFs in the repository's `biblio/` directory.

## Add a nuclear reaction

See {doc}`howto/custom-networks` — either the GUI's "Manage networks" dialog
or the `custom_network=` API for a run-time toggle, or a permanent addition
to the source tree (drop a rate table under
`primat/data/nuclear/tables/<name>/`, add `<name>` to a network file under
`primat/data/nuclear/networks/`).

For a one-off sensitivity study rather than a permanent addition, use the
existing `p_<reaction>`/`delta_<reaction>` config knobs to perturb a rate
already in the network — see {doc}`howto/rate-variation-mc` — no file
changes needed.

## Add a dark-sector / non-standard background component

See {doc}`howto/backgrounds` for the three mechanisms (`extra_rho`,
`custom_background`, `background=`) and when to reach for each.

## Add a neutrino-history variant

Neutrino-temperature/spectral-distortion history is dispatched by
`make_neutrino_history(cfg, plasma)` (`primat.neutrino_history`):

1. Base regime: `NEVOTable` (`cfg.incomplete_decoupling=True`, the default —
   reads the non-instantaneous-decoupling table) or
   `InstantaneousDecoupling` (analytic `T_ν(T_γ)` from EM entropy
   conservation).
2. Optionally wrapped in `AnalyticDistortion` when
   `cfg.spectral_distortions and cfg.analytic_distortions`, layering an
   analytic μ/y-type spectral distortion on top of the base regime.

A new variant — e.g. a different non-instantaneous-decoupling table, or a
new analytic distortion shape — is a new class exposing the same interface
(`Tnue_of_Tg`/`Tnumu_of_Tg`/`Tnutau_of_Tg`/`N_NEVO_of_Tg`/`dFDneu_func`/
`rho_nu_SD`) plus a new branch in `make_neutrino_history`'s dispatch.

For swapping the *data* underlying the existing `NEVOTable`/
`InstantaneousDecoupling` classes without writing a new class, see
{doc}`howto/nevo-tables` (`nevo_file`/`nevo_spectral_file`/
`nevo_grid_file`/`nevo_file_prefix`).
