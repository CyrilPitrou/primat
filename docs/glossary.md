# Glossary

The shorthand this project uses, in one line each. Terms are expanded at
first use in every file; where a file uses one repeatedly, it links here
instead of repeating the expansion.

```{glossary}
BBN
  Big Bang Nucleosynthesis — the formation of the light nuclei (D, ³He, ⁴He,
  ⁷Li) in the first ~20 minutes of the universe. What this code computes.

Y_i
  The abundance of nuclide *i* per baryon, `Y_i = n_i/n_b` (dimensionless).
  A *number* fraction, not a mass fraction: the mass fraction is `A_i Y_i`,
  and `sum(A_i Y_i) = 1`. The `Y_final` dict and the `Y_<nuclide>` columns of
  the time-evolution output both hold `Y_i`.

YP
  The primordial helium-4 *mass* fraction, `YP = 4 Y_He4` (dimensionless).
  Reported in two conventions: `YPBBN` (BBN convention, what this code
  integrates to) and `YPCMB` (CMB convention, converted via the He4/H ratio).

D/H
  Deuterium-to-hydrogen ratio by number, `Y_H2 / Y_p`. Likewise `He3/H`,
  `Li7/H`, `He3/He4`, `Li6/Li7`; result-dict keys `DoH`, `He3oH`, `Li7oH`,
  `He3oHe4`, `Li6oLi7`.

Neff
  Effective number of neutrino species (dimensionless) — the relativistic
  energy density beyond photons, expressed in units of one instantaneously
  decoupled neutrino flavour. 3.044 in the standard model.

Omegabh2
  The baryon density parameter `Ω_b h²` (dimensionless). primat's default,
  0.02242, is the Planck 2018 + BAO value.

eta_b
  The baryon-to-photon ratio `n_b/n_γ`; `eta10 = 10¹⁰ η_b`. Fixed by
  `Omegabh2` through a conversion constant, and the x-axis of a Schramm
  diagram.

tau_n
  The free-neutron lifetime [s], default 878.4. Normalises the n↔p weak rates
  when `tau_n_normalization=True`.

T9
  Temperature in units of 10⁹ K (equivalently GK). The unit every nuclear
  rate table is tabulated in.

HT / MT / LT
  The three temperature eras the nuclear network is integrated across:
  **HT** (high, T > ~1 MeV) evolves n and p only; **MT** (middle, down to
  ~0.1 MeV) uses the chosen network intersected with a fixed 18-reaction
  list, the full network being too stiff there; **LT** (low, down to
  ~0.001 MeV) uses the chosen network in full.

n↔p
  The weak interconversion of neutrons and protons (six processes: β decay,
  electron and neutrino capture, and their inverses). Its rates set the
  neutron fraction at freeze-out and therefore `YP`.

Born
  The lowest-order n↔p rate, with no radiative, finite-mass or
  spectral-distortion corrections. What `radiative_corrections=False`
  selects.

CCR
  Coulomb and T = 0 resummed radiative corrections to the n↔p rates
  (`radiative_corrections=True`, the default).

FM
  The finite-nucleon-mass correction to the n↔p rates — a Fokker-Planck
  expansion to first order in `T/m_N`, including weak magnetism
  (`finite_mass_corrections=True`).

CCRTh
  The finite-*temperature* radiative correction to the n↔p rates (Brown &
  Sawyer 2001), including bremsstrahlung (`thermal_corrections=True`). The
  expensive one: computing it needs a `vegas` Monte-Carlo integration, which
  is why it has a cache of its own.

SD
  Spectral distortions — the correction for neutrino distributions that are
  not exactly Fermi-Dirac after decoupling (`spectral_distortions=True`).

NEVO
  The pre-computed neutrino-evolution tables under `data/NEVO/`, holding the
  non-instantaneous decoupling history: the per-flavour temperature ratios
  `T_ν(T_γ)`, the heating function, and (in the spectral tables) the
  distorted neutrino spectra. Read when `incomplete_decoupling=True`.

QED corrections
  Finite-temperature QED corrections to the electromagnetic plasma's
  pressure and energy density (`QED_corrections=True`). Tabulated under
  `data/cache_plasma_weak/plasma/`.

amax
  Maximum mass number `A`. Setting it drops every reaction involving a
  nuclide heavier than `A` from *any* network — `network="large", amax=8`
  gives 68 reactions.

p_&lt;reaction&gt;
  Per-reaction log-normal variation parameter, e.g. `p_n_p__d_g`. The rate
  becomes `median × exp(p × σ)`, so `p = 1` is a +1σ rate. Monte Carlo
  samples these from N(0,1).

delta_&lt;reaction&gt;
  Per-reaction *additive* rescaling, active when
  `rescale_nuclear_rates=True`: the rate becomes `median × (1 + delta)`. For
  deterministic sensitivity studies, where `p_<reaction>` is for uncertainty
  propagation.

expsigma
  The log-normal uncertainty width σ of a rate, read from the third column
  of its rate table. The σ in `median × exp(p × σ)`.

MC
  Monte Carlo — repeated solves with randomly sampled rates and `tau_n`,
  used to propagate nuclear-rate uncertainty onto the observables
  (`run_mc()`, or `--mc N`).

master T9 grid
  The single temperature grid every rate table is resampled onto at load
  time (1000 log-spaced points over T9 ∈ [10⁻³, 10], by default), so tables
  tabulated on different grids can be mixed freely.

fingerprint
  The hash of every configuration field that can change a cached table's
  contents, stored in the cache file's header and in its filename. A cache
  whose fingerprint does not match the current configuration is recomputed,
  never used.

overlay
  A directory searched *before* the shipped data tree, falling back to it on
  a miss: `cache_dir` for the regenerable caches, `user_nuclear_dir` for
  networks and rate tables. `data_dir` is not an overlay — it replaces the
  shipped tree outright.

backend
  One of the two interchangeable solver implementations: the compiled **C**
  engine (default, faster) and the pure-**Python** one (fallback, and the
  only one supporting a custom `Background` object). Selected with
  `force_backend` or `--backend`.

network
  The set of nuclear reactions integrated: `small` (12 reactions, 8
  nuclides), `small_parthenope` (the same 12 with Parthenope 3.0 rate
  tables), or `large` (~429 reactions, ~59 nuclides), optionally restricted
  by `amax`.

AC2024
  The 2024 nuclear-reaction rate compilation the `large` network's tables
  are derived from.

Parthenope
  Another BBN code; `small_parthenope` uses its 3.0 rate tables for the same
  12 reactions as `small`, for comparison runs.

detailed balance
  The thermodynamic relation fixing each reverse reaction rate from its
  forward rate and the nuclear Q value, so only forward rates are tabulated.

NSE
  Nuclear statistical equilibrium — the abundance distribution reached when
  every reaction is fast compared with the expansion, used as the initial
  condition and as an independent check on the nuclide data.

Schramm diagram
  The standard BBN figure: primordial abundances plotted against `eta_b` (or
  `Omegabh2`), with observational constraints overlaid.
```
