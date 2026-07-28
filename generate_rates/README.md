# generate_rates/

Offline scripts that produce the rate tables shipped under `primat/data/`.
None of these run at BBN-solve time — they are run once (or whenever an
upstream source changes) to regenerate committed data files.

## Pipeline map

- **`convert_ac2024_rates.py`** — the main entry point. Builds the
  500-point log-uniform-T9 rate-table set from `BBNRatesAC2024.dat` (tabulated
  reactions) plus a hard-coded analytic-rate table (reactions PRIMAT evaluates
  in closed form). Produces:
  - `primat/data/nuclear/tables/*.txt` — per-reaction rate tables.
  - `primat/data/csv/*.csv` — `nuclides.csv`, `reactions_large.csv`,
    `detailed_balance.csv`.
  - `primat/data/nuclear/networks/large.txt` — the large-network reaction
    list.

  It imports `nuclide_table.py` (nuclide property table, NUBASE2020 parsing,
  detailed-balance helper) and, through it, `nuclear_data.py` (the
  `detailed_balance` reverse-rate-coefficient formula).

- **`generate_qed_tables.py`** — independent of the above. Recomputes the
  analytic QED plasma-pressure correction tables (`primat.qed_pressure`)
  and writes them to `primat/data/plasma/`.

- **`parthenope3.0_extract/`** — a separate, self-contained sub-pipeline that
  extracts the 12 primat *small-network* rates directly from the Parthenope
  3.0 Fortran source (verbatim code fragments, not retyped formulas). See its
  own `parthenope3.0_extract/README.md` for the method and how to run it.

- **`migrate_tables_to_folders.py`** — a one-off, already-applied migration
  that flattened `tables/<name>.txt` into `tables/<name>/<name>.txt` (the
  per-reaction-folder layout used today). Kept for reference/reproducibility
  in case a fresh `--keep-source-grid` export ever needs the same treatment.

- **`thermal_average.ipynb`** — the user-facing entry point for adding a rate
  from your *own* cross-section data, as opposed to the bulk regeneration
  above. Give it an astrophysical S-factor `S(E)` (or a cross-section
  `σ(E)`) plus a parameter covariance, and it computes the
  Maxwell–Boltzmann-averaged `N_A<σv>` on the 60-point Wagoner temperature
  grid, propagates the uncertainty by Monte Carlo into the table's
  multiplicative error column, and writes a primat-format table with the
  correct detailed-balance header. All nuclide data (masses, charges, spins,
  Q-values) is read from a live `PRIMATConfig`, so it cannot drift from the
  solver's own nuclear data.

  Output goes to `generate_rates/rate_tables_out/` (gitignored), laid out as a
  ready-made `user_nuclear_dir` overlay — the shipped `primat/data/` tree is
  never modified. The notebook's last section builds a matching network list
  file and runs BBN with the new rate to show the effect. This notebook
  supersedes the Mathematica `Thermal-Average.nb`.

- **`PRIMAT-Main.m`** — the original Mathematica source
  `convert_ac2024_rates.py`'s hard-coded analytic-rate table was extracted
  from. Kept only as the input to `--dump-analytic`, which regenerates that
  hard-coded table if `PRIMAT-Main.m` changes; not needed for a normal
  rate-table rebuild.

- **`nubase_4.mas20.txt`**, **`BBNRatesAC2024*.dat`** — upstream source data
  (NUBASE2020 mass table; the AC2024 tabulated-rate compilation).
