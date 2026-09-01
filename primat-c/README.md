# primat-c

A standalone C99 implementation of the PRIMAT Big Bang Nucleosynthesis (BBN) solver.

Two kinds of reader come here. Most users never do: `pip install primat` builds
this tree into a compiled extension and `primat.backend.run_bbn()` uses it
automatically. You are in the right place if you want to **run BBN without a
Python interpreter** — on a cluster node, inside a C or Fortran pipeline, from a
Makefile — or if you are **changing the solver** and need the C half of a
physics change that must land on both backends.

## Overview

`primat-c/` is a complete C99 port of the PRIMAT BBN solver, providing:

- **Identical physics**: the same observables as the Python backend, held to a
  `rel=5e-5` cross-backend D/H budget by `tests/test_backend_parity.py` (the
  measured gap is ~7e-6; `Neff` agrees to every digit)
- **High performance**: ~10× faster than the pure-Python backend on a warm
  `small` run and ~5× on `large, amax=8` (`examples/baseline_timings.txt`)
- **Standalone usage**: Can be compiled and run independently of the Python package
- **Python integration**: Compiled as an extension to provide the default fast backend for `primat.backend.run_bbn()`

## Compilation

### Prerequisites

- C99-compatible compiler (GCC, Clang, MSVC)
- GNU Make (for the provided Makefile)
- Python development headers (for the Python extension bridge, optional)

### Building the standalone C executable

From the `primat-c/` directory:

```bash
cd primat-c
make
```

This produces the `primat-c` executable in the `build/` directory.

**Available targets:**
- `make` or `make all` - Build the standalone executable (optimized, `-O2` by default)
- `make clean` - Remove build artifacts
- `make debug` - Build with debug symbols and sanitizers instead
- `make test` - Build and run the unit test suite (27 programs). This does
  *not* build the standalone executable; the Python suite's CLI tests need
  it, so run plain `make` too
- `make debug-test` - The same, instrumented with ASan and UBSan
- `make bench` - Build and run the timing benchmark
- `make leak-test` - Build and run the memory-leak check under a sanitizer
- `make fuzz RUNS=N` - Fuzz every parser that reads a user-supplied file, under
  AddressSanitizer and UndefinedBehaviorSanitizer (see `tests/fuzz/README.md`)
- `make fuzz-coverage` - Report what each fuzz target actually reached

### Platform-specific notes

#### Linux (GCC/Clang)

```bash
# Install build essentials on Debian/Ubuntu
sudo apt-get install build-essential

# Then compile as above
make
```

#### macOS (Clang)

```bash
# Install Xcode command line tools if not already present
xcode-select --install

# Then compile
make
```

#### Windows (MinGW/MSYS2)

The standalone executable's Makefile is POSIX-only (no MSVC/`nmake` project is
provided); build it under MSYS2's MinGW toolchain instead:

```bash
pacman -S mingw-w64-x86_64-gcc
make
```

### Building the Python extension

The Python extension is automatically built during `pip install primat` and included in the wheel distribution. To build manually:

```bash
# from the repository root, not from primat-c/
python setup.py build_ext --inplace
```

## Usage

### Standalone C executable

After compilation, run from the `primat-c/` directory:

```bash
# Basic usage with default parameters
./build/primat-c

# With custom parameters
./build/primat-c --Omegabh2 0.02242 --network large --amax 8

# List all available options
./build/primat-c --help
```

**Common options:**
- `--Omegabh2 VALUE` - Baryon density Ω_b h² (default: 0.02242)
- `--DeltaNeff VALUE` - Extra relativistic degrees of freedom (default: 0)
- `--network NAME` - Nuclear reaction network: small, small_parthenope, large (default: small)
- `--amax N` - Maximum mass number A for reactions (filters any network)
- `--numerical_precision VALUE` - ODE solver relative tolerance (default: 1e-7)
- `--output_file PATH` - Where to write the time-evolution TSV
- `--ini FILE` - Read parameters from an INI file instead of flags
- `--data_dir DIR` - Use a different data tree (see "Data directory" below)
- `--json` - Output results as JSON

### Configuration file

Instead of command-line arguments, you can use an INI-style configuration file:

```bash
# Run with a configuration file
./build/primat-c --ini examples/run_basic.ini
```

### Via Python API

The C backend is automatically used as the default by `primat.backend.run_bbn()`:

```python
from primat.backend import run_bbn

# Automatically uses C backend if available
result = run_bbn({"Omegabh2": 0.02242, "network": "large"})

# Force C backend explicitly
result = run_bbn({"Omegabh2": 0.02242}, force_backend="c")
```

### Custom networks

Custom nuclear reaction networks can be used with both the standalone executable and Python API:

```bash
# Using a custom network file
./build/primat-c --network my_custom_network
```

The network file should be placed in `data/nuclear/networks/` or made accessible via `--user_nuclear_dir` (additive overlay) or `--data_dir` (full data-tree replacement).

## Output

The C backend produces identical output to the Python backend:

### Return values (via Python API)

`run_bbn()` returns a dictionary with:

- `YPBBN` - Helium-4 mass fraction (BBN convention)
- `YPCMB` - Helium-4 mass fraction (CMB convention)  
- `He4oH` - He4/H ratio (by number)
- `DoH` - D/H ratio
- `He3oH` - (He3+T)/H ratio
- `Li7oH` - (Li7+Be7)/H ratio
- `Neff` - Effective number of neutrino species
- `Omeganurel` - Ω_ν h² × 10⁶ (relativistic)
- `OneOverOmeganunr` - 1 / (Ω_ν h² × 10⁻⁶) (non-relativistic)
- `evolution` - Time evolution data (when `output_time_evolution=True`)

### Command-line output

Default console output:
```
────────────────────────────────────────────────────
          PRIMAT results at T = 0.001 MeV
────────────────────────────────────────────────────
Neff       = 3.04397730
YP (BBN)   = 0.24699907
YP (CMB)   = 0.24567276
He4/H      = 8.2011454e-02
D/H        = 2.4358767e-05
He3/H      = 1.0399348e-05
He3/He4    = 1.2680361e-04
Li7/H      = 5.557664e-10
--- running time: 0.06 seconds ---
```

`Li6/Li7` and `YCNO` appear only when the network tracks them, i.e. with
`--network large`. The running time is whatever this machine took; the
abundances are the `small` network's, pinned by `tests/README.md`'s validation
reference.

With `--json` flag, full results are output as JSON.

## Data directory structure

`primat-c/` ships no data of its own: both backends read the one tree inside
the Python package, `primat/data/`.

```
primat/data/
  nuclear/
    tables/              # Per-reaction rate tables (one folder per reaction)
    networks/            # Network list files (small.txt, large.txt, etc.)
  csv/                   # Reaction catalog (nuclides.csv, detailed_balance.csv, reactions_large.csv)
  NEVO/                  # Neutrino-decoupling history tables
  cache_plasma_weak/
    plasma/              # Pre-computed QED pressure and electron-thermo tables
    weak/                # Cached n<->p forward/backward rates
```

The standalone executable finds that tree without being told, in this order:

1. `--data_dir DIR`, if given;
2. the `CPRIMAT_DATA_DIR` environment variable;
3. the sibling package next to the executable itself — `build/primat-c` resolves
   `../../primat/data`. This is anchored to where the binary lives, **not** to
   the working directory, so `./build/primat-c` gives the same answer from any
   directory.

Two further knobs redirect parts of the tree: `--user_nuclear_dir DIR` overlays
networks and rate tables only, and `--cache_dir DIR` sends the two regenerable
cache trees somewhere writable. Used as a Python extension, the root comes from
the installed package.

## Code structure

Every header carries the contract for what it declares — ownership of
pointers, units, error returns — so `include/` is the API document and `src/`
need not be read to use the library.

```
primat-c/
  include/            # 27 headers: one per module, plus the two compat shims
    api.h             # cprimat_run: one full BBN solve, and its CPRResults
    config.h          # CPRConfig: every parameter, and the by-name setter
    cli.h             # the executable's entry point
    mc.h              # threaded Monte-Carlo rate/tau_n uncertainty propagation
    background.h      # a<->t<->T, the Friedmann rate, the neutrino sector
    nuclear_network.h # the HT/MT/LT era integration
    network_data.h    # network files, rate tables, the solver-ready network
    weak_rates.h      # n<->p rate tables and their corrections
    plasma.h          # photon/e+-/neutrino thermodynamics
    ...               # and the numerics: ode_bdf, ode_rk, spline, linalg,
                      # quad, vegas, rng, cache, table_io, xalloc, ...

  src/                # 26 files: one per module header, plus main.c
    api.c             # cprimat_run, cpr_assemble_results, the verbose banner
    cli.c             # argument parsing, the printed report, --json, --mc
    mc.c              # the MC worker threads
    config.c          # the field table, validation, path resolution
    network_data.c    # cpr_load_network and the rate-table pipeline
    nuclear_network.c # cpr_nuclear_network_solve and the output writers
    background.c      # cpr_bg_init_standard / _custom and the queries
    weak_rates.c      # the Born/CCR/FM/SD/CCRTh rate integrals

  tests/
    unit/             # 27 test programs, run by `make test`
    fuzz/             # 11 fuzz targets + a coverage-guided engine (`make fuzz`)

  examples/           # run_basic.ini and two reference-run configurations
```

Ownership convention, uniform across the library: every `cpr_*_free` releases
what the struct owns but **never the struct itself**, and none of them accepts
`NULL`. A `char **errmsg` out-parameter is always the caller's to `free`.

## Custom backgrounds and advanced usage

The C backend supports the same extension points as the Python one, all as
ordinary parameters:

- `custom_background FILE` — a user-supplied `(T, t, a)` table replacing the
  standard cosmology. Setting it forces `incomplete_decoupling=False` and
  `spectral_distortions=False` (custom-background mode uses
  instantaneous-decoupling weak rates) and says so on stderr.
- `nevo_file` / `nevo_spectral_file` / `nevo_grid_file` / `nevo_file_prefix` —
  alternative neutrino-decoupling tables.
- `p_<reaction>` / `delta_<reaction>` — per-reaction rate variation, the basis
  of `--mc` and of sensitivity analysis. `--list-reactions` prints the keys the
  selected network accepts.
- `external_scale_factor`, `decay_era`, `fEDE`/`zcEDE`/`wnEDE` — see
  `--list-params` for the full set with its current values, in a spelling an
  INI file accepts back.

`include/` carries the contract for each of these; `docs/howto/` documents the
physics behind them.

## Error handling and debugging

### Verbose output

Use the `--verbose` flag to see detailed progress messages:

```bash
./build/primat-c --verbose
```

This shows timing information, cache hits, and integration progress.

### Debug builds

To compile with debug symbols and no optimizations:

```bash
make debug
```

### Checking for memory leaks

```bash
make leak-test
```

Builds `tests/unit/test_memory_stress.c` (-O0 -g) and runs it under the
platform's leak checker: macOS's `leaks --atExit`, or `valgrind
--leak-check=full` elsewhere. The test cycles `cprimat_run`/
`cpr_mc_uncertainty`/`cpr_config_init_defaults` (success path, error path,
and MC) several times over in one process, each time freeing everything
via the matching `*_free` call, so a leak of even a few hundred bytes per
call accumulates into something the checker reliably flags. Exits nonzero
if any leak is found.

### Checking installation

To verify the C backend is working correctly:

```python
from primat.backend import HAS_C_BACKEND
print(f"C backend available: {HAS_C_BACKEND}")
```

## Version compatibility

The `primat-c` version is kept in sync with the main `primat` package version. The version is defined in `include/config.h` as `CPRIMAT_VERSION` and must match the version in `pyproject.toml`.

## Support and contribution

- Report issues on the main PRIMAT repository
- Contributions to the C backend should maintain parity with the Python implementation
- Any changes to physics or numerics must be mirrored in both backends

## License

Same as the main PRIMAT package - see the repository LICENSE file.