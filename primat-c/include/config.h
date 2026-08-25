/* config.h -- CPRIMAT run-time configuration (port of primat/config.py).
 *
 * Unlike Python's dynamically-typed PRIMATConfig, CPRConfig is a single plain
 * struct with one typed field per DEFAULT_PARAMS entry: C has no convenient
 * dynamic-attribute story, and a struct is both simpler to read and faster
 * to access than threading every physics formula through a generic
 * key/value lookup. The *external* interface (ini file, CLI flags, --set)
 * still goes through a generic tagged-union CPRParam (see cpr_parse_literal
 * / cpr_config_set_by_name below) -- that
 * union is the parsing/dispatch boundary, not the storage representation.
 *
 * Optional ("None"-able) Python values are represented as:
 *   - string-typed param, no value     -> NULL char* (nevo_file, ...)
 *   - amax (int-or-None)               -> -1 sentinel (Python requires a
 *                                          positive int when set, so -1 is
 *                                          unambiguous)
 */
#ifndef CPRIMAT_CONFIG_H
#define CPRIMAT_CONFIG_H

/* MUST be kept in sync with pyproject.toml's `version`; bump this by hand
 * alongside pyproject.toml whenever the package version changes. Checked by
 * tests/test_docs_consistency.py::test_cprimat_version_matches_pyproject. */
#define CPRIMAT_VERSION "0.3.2"

/* Buffer sizes for filesystem paths built by joining cfg->data_dir (see
 * CPRConfig::data_dir below) with one or two relative path components
 * (e.g. "<data_dir>/plasma/QED_pressure_correction_e2.txt"). Sized with
 * headroom beyond data_dir's own capacity plus a generous per-component
 * margin so -Wformat-truncation can prove the snprintf join never
 * truncates, even though real data-dir paths are always far shorter. */
#define CPR_DATA_DIR_LEN 4096
#define CPR_PATH_BUF_LEN  (CPR_DATA_DIR_LEN + 256)  /* data_dir + "/" + one component */
#define CPR_PATH_BUF_LEN2 (CPR_PATH_BUF_LEN + 256)  /* + a second joined component */

#include <stddef.h>

#include "constants.h"

/* The four fingerprinted caches. Each is keyed on its own subset of the
 * physical constants -- the ones it actually reads -- so a constant it does
 * not read cannot invalidate it. Mirrors cache_utils.CACHE_CONSTANTS; the
 * subsets themselves live in cpr_constants_hash (cache.c). */
typedef enum {
    CPR_CONSTS_WEAK = 0,      /* nTOp_<hash>.txt */
    CPR_CONSTS_THERMAL,       /* nTOp_thermal_<hash>.txt */
    CPR_CONSTS_ELEC_THERMO,   /* electron_thermo_<hash>.txt */
    CPR_CONSTS_QED,           /* QED_pressure_correction_e{2,3}.txt */
    CPR_CONSTS_N_CACHES
} CPRConstsCache;

/* ---- Generic tagged-union value, used only at the parsing/CLI/ini
 * boundary (cpr_parse_literal, cpr_config_set_by_name). ---- */
typedef enum { CPR_NONE, CPR_BOOL, CPR_INT, CPR_DOUBLE, CPR_STRING } CPRType;

typedef struct {
    CPRType type;
    union {
        int b;          /* CPR_BOOL: 0/1 */
        long i;         /* CPR_INT */
        double d;       /* CPR_DOUBLE */
        const char *s;  /* CPR_STRING; not owned -- caller-managed lifetime */
    } v;
} CPRParam;

/* A single named (key, value) pair -- the unit ini/cli parsing produces. */
typedef struct {
    const char *key;   /* not owned */
    CPRParam value;
} CPRParamSet;

/* Parses one literal token the same way primat.cli's --set escape hatch
 * does (ast.literal_eval-equivalent): try int, then float, then
 * true/false/none (case-insensitive), else fall back to the literal string
 * (quotes, if any, are stripped).
 *
 * LIFETIME (important): a CPR_STRING result points into the caller's `buf`
 * (of `bufsize` bytes), NOT into `s` and NOT into fresh storage, so it is
 * valid exactly as long as that buffer is. Give each parse whose value you
 * keep its own buffer. Any caller retaining a value past the immediately
 * following cpr_config_set_by_name (which strdup's its own copy) must copy
 * it; cpr_paramlist_add below does exactly that. */
CPRParam cpr_parse_literal(const char *s, char *buf, size_t bufsize);

/* ---- A retained, self-owning list of (key, value) overrides ----
 *
 * The MC driver re-applies the user's overrides to a *fresh* CPRConfig in
 * every worker thread (mc.c's worker_setup), so the CLI/ini front end has to
 * hand it the complete override set as CPRParamSet[] -- long after argv and
 * the callers' parse buffers have gone out of scope. This list copies both
 * halves of every pair into its own storage, so a retained entry can never
 * alias argv, an ini line buffer, or the previous parse.
 *
 * Keys longer than CPR_PARAM_KEY_LEN-1 / string values longer than
 * CPR_PARAM_VAL_LEN-1 are truncated (both bounds are far above any real
 * parameter name or path). */
#define CPR_PARAM_KEY_LEN 256
#define CPR_PARAM_VAL_LEN 1024

typedef struct {
    CPRParamSet *items;                        /* items[i].key -> key_store[i] */
    char (*key_store)[CPR_PARAM_KEY_LEN];
    char (*val_store)[CPR_PARAM_VAL_LEN];      /* CPR_STRING values only */
    size_t n, cap;
} CPRParamList;

/* Appends a copy of (key, value). Grows on demand and re-points every
 * previously stored entry at its (possibly moved) storage, so
 * `pl->items` stays a valid CPRParamSet array of length `pl->n` after any
 * number of adds. */
void cpr_paramlist_add(CPRParamList *pl, const char *key, CPRParam value);

/* Releases the three backing arrays and zeroes the list. */
void cpr_paramlist_free(CPRParamList *pl);

/* Small open dictionary for p_<rxn> / delta_<rxn>, mirroring
 * PRIMATConfig.p_rxn / delta_rxn. Linear-scan array: the reaction count is
 * at most ~430 (the "large" network), so a hash table buys nothing here. */
typedef struct {
    char name[40];
    double value;
} CPRRxnEntry;

typedef struct {
    CPRRxnEntry *entries;
    size_t n, cap;
} CPRRxnMap;

double cpr_rxnmap_get(const CPRRxnMap *map, const char *name); /* 0.0 default */
void cpr_rxnmap_set(CPRRxnMap *map, const char *name, double value);
void cpr_rxnmap_free(CPRRxnMap *map);

/* One nuclide row from data/csv/nuclides.csv. */
typedef struct {
    char name[16];
    int N, Z;
    double mass_excess_keV;
    double spin;
} CPRNuclide;

typedef struct {
    CPRNuclide *items;
    size_t n;
} CPRNuclideTable;

/* ------------------------------------------------------------------------
 * CPRConfig: every DEFAULT_PARAMS entry as a typed field, grouped exactly
 * as in config.py's DEFAULT_PARAMS dict (comments there explain each flag
 * in physics terms; not repeated here -- see config.py).
 * ------------------------------------------------------------------------ */
typedef struct {
    /* ---- general behaviour and numerical settings ---- */
    int verbose;
    int show_progress; /* print [primat] HT./MT./LT./done. phase markers to stderr (default 1; suppressed during MC batch samples) */
    int debug;
    double numerical_precision;
    int use_numba; /* unused in C (no JIT path); kept for CLI/ini parity */
    int strict_params; /* unknown-key policy for the INI/CLI loaders (default 0). Mirrors
                        * PRIMATConfig.strict_params: 0 = warn-and-ignore, 1 = fatal. Stored
                        * for round-trip/wrapper parity; the standalone C loaders already treat
                        * unknown keys as fatal, so this only downgrades that when 0. */

    /* ---- neutrino decoupling ---- */
    int incomplete_decoupling;

    /* ---- electromagnetic plasma ---- */
    int QED_corrections;
    int n_electron_table;
    int recompute_electron_thermo;
    int recompute_qed_corrections;

    /* ---- spectral distortions ----
     * analytic_distortions / y_SZ / y_gray select the closed-form y-type
     * (SZ/Compton) + gray-type distortion (neutrino_history.
     * AnalyticDistortion), an alternative to the default NEVO-spectrum-
     * table distortion; PRIMATConfig pairs analytic_distortions=True with
     * incomplete_decoupling=False (cpr_config_validate enforces this).
     * (There is deliberately no mu-type / delta_xi_nu distortion: a
     * genuine neutrino chemical potential is munuOverTnu, which IS
     * ported -- it shifts the weak rates and, via
     * cpr_rho_nu_chempot_excess, the neutrino energy density / Neff.) */
    int spectral_distortions;
    int analytic_distortions;
    double y_SZ;
    double y_gray;

    /* ---- custom NEVO tables (NULL = unset / use shipped default) ---- */
    char *nevo_file;
    char *nevo_spectral_file;
    char *nevo_grid_file;
    char *nevo_file_prefix; /* never NULL; defaults to "NEVOPRIMAT" */

    /* ---- background mode ---- */
    int external_scale_factor;
    char *custom_background; /* NULL = not set */

    /* ---- fundamental constants (overridable) ---- */
    /* Internal, not a DEFAULT_PARAMS key and not part of any fingerprint:
     * set by the Python bridge (_wrapper.c) so config warnings are printed
     * once, by primat/config.py's _warn_off_default_risks, instead of twice
     * (Python's warning plus this backend's stderr copy). The standalone CLI
     * leaves it 0 and prints them itself. */
    int suppress_config_warnings;

    double GN; /* natural units [MeV^-2]; set/read via cpr_config_set_GN()/
                * cpr_config_get_GN(), which convert to/from the SI units
                * [m^3 kg^-1 s^-2] that DEFAULT_PARAMS["GN"] and
                * cpr_config_set_by_name("GN", ...) use -- never assign this
                * field directly with an SI-unit value. */

    /* Physical constants for THIS run. Seeded from g_const by
     * cpr_config_init_defaults; the 16 measured fields are then settable by
     * name (alphaem, me, gA, ...) like any other parameter, while the ten
     * exact ones have no field-table entry and stay at their defaults. Every
     * solver read of a measured constant goes through here -- never through
     * g_const -- so coexisting configs and MC worker threads cannot see each
     * other's values. */
    CPRConstants consts;

    /* One 16-hex-digit hash per fingerprinted cache, each over only the
     * constants that cache reads (cache_utils.CACHE_CONSTANTS on the Python
     * side), refreshed by cpr_config_refresh_constants whenever a constant is
     * set. Every fingerprint builder embeds its own, and stores this very
     * pointer in a CPRFPField -- which is why they live here, with the
     * config's lifetime, rather than in a builder-local buffer. */
    char consts_hash[CPR_CONSTS_N_CACHES][17];

    /* ---- background thermodynamics ---- */
    double T_start_cosmo_MeV;
    double T_end_MeV;
    int sampling_temperature_per_decade;

    /* ---- n <-> p weak rates ---- */
    int radiative_corrections;
    int finite_mass_corrections;
    int thermal_corrections;
    int weak_rate_cache;
    int save_nTOp;
    int sampling_nTOp_per_decade;
    int save_nTOp_thermal;
    int sampling_nTOp_thermal_per_decade;
    int tau_n_normalization;
    double tau_n;
    double std_tau_n;
    int vegas_n_eval;     /* evaluations per VEGAS iteration, see vegas.h */
    int vegas_n_itn;      /* VEGAS warmup/measure iterations, see vegas.h */
    double epsrel_thermal;

    /* ---- output options ---- */
    int output_time_evolution;
    int output_rates_time_evolution;
    int output_n_points;
    char *output_file;
    int output_final_result;
    char *output_final_file;
    int output_background_evolution;
    char *output_background_file;
    int output_mc_samples;
    int output_mc_covariance;
    int output_mc_correlation;
    char *output_mc_file_prefix; /* stem for <prefix>_samples/_covariance/_correlation.tsv */

    /* ---- nuclear network ---- */
    int rate_grid_npts;
    double rate_grid_T9_min;
    double rate_grid_T9_max;
    char *network;
    int amax; /* -1 = None (no filter); else positive int */
    double atol_LT;
    /* Cap on the MC rate rescaling factor: variation is clamped to [1/cap, cap]
     * before multiplying the median rate.  0.0 = no cap (mirrors Python None). */
    double mc_rate_rescale_cap;
    int nuclear_qed_corrections;

    /* ---- nuclear overlay (mirrors PRIMATConfig.user_nuclear_dir; see
     * docs/howto/data-overlays.md). NULL = unset (shipped data/nuclear/
     * tree only). Wired through cpr_config_resolve_rates_path() at the same
     * two call sites as the Python side: the network-file path
     * (nuclear/networks/<name>.txt) and each reaction's rate-table
     * file (nuclear/tables/<rxn>/<file>) -- NOT the reaction catalog
     * (nuclides.csv/reactions_large.csv/detailed_balance.csv) or decays.txt,
     * which stay on data_dir. Overlay roots are treated as the equivalent of
     * `primat/data/nuclear`, so they should contain `networks/` and `tables/`
     * directly.  The full data-tree takeover (PRIMATConfig.data_dir) is handled
     * at the C level by cpr_config_init_defaults(data_dir): the Python
     * backend.py passes cfg.resolved_data_dir there, so data_dir already
     * reflects any user override before any field is set. */
    char *user_nuclear_dir;  /* additive nuclear overlay, checked before the shipped default */

    /* ---- writable cache redirect (mirrors PRIMATConfig.cache_dir).
     * NULL = unset: both regenerable cache trees (the n<->p weak rates and the
     * plasma electron-thermo/QED tables) live under
     * <data_dir>/cache_plasma_weak/{weak,plasma}/. When set, cache files are
     * READ from <cache_dir>/{weak,plasma}/ first and, on a miss, from the
     * shipped <data_dir>/cache_plasma_weak/ (overlay -- shipped caches never
     * shadowed), and WRITTEN only to <cache_dir> (created on demand). Set it
     * when the install location is read-only. Wired through
     * cpr_config_cache_write_dir()/cpr_config_resolve_cache_file() below. Cache
     * LOCATION only: never part of any fingerprint. */
    char *cache_dir;

    /* ---- cosmological inputs ---- */
    double Omegabh2_; /* backing field; use cpr_config_set_Omegabh2() to set
                          (mirrors the Python @property that recomputes
                          eta0b on assignment) */
    double Omegach2;
    double h;
    double DeltaNeff;
    double munuOverTnu;
    /* Per-flavour neutrino chemical potentials. Each defaults to NAN,
     * the "inherit munuOverTnu" sentinel (mirrors Python's None). Read the
     * EFFECTIVE per-flavour value through cpr_config_xi_nu_e/mu/tau(), which
     * resolve NAN -> munuOverTnu. Only xi_e enters the n<->p weak rates (nu_e
     * appears in n <-> p + e + nu_e); all three enter the neutrino energy
     * density / Neff (cpr_rho_nu_chempot_excess). Set via cpr_config_set_by_name
     * as F_DOUBLE_OR_NAN (None -> NAN). */
    double munuOverTnu_e;
    double munuOverTnu_mu;
    double munuOverTnu_tau;

    /* ---- decay-era options. The decay_era Decay-Time propagation is now
     * implemented (see nuclear_network.h's cpr_nuclear_network_decay_era):
     * decay_era + output_decay_evolution + the `large` network writes the
     * output_decay_file TSV. decay_reverse_rates still only affects rate
     * loading, not a distinct solver era. ---- */
    int decay_reverse_rates;
    int decay_era;
    double t_decay_end;
    int decay_n_points;
    int output_decay_evolution;
    char *output_decay_file;

    /* ---- Early Dark Energy ---- */
    double fEDE;
    double zcEDE;
    double wnEDE;

    /* ------------------------------------------------------------------
     * Derived / non-DEFAULT_PARAMS state
     * ------------------------------------------------------------------ */
    double Omegabh2_to_eta0b;
    double eta0b;

    CPRRxnMap p_rxn;
    CPRRxnMap delta_rxn;

    CPRNuclideTable nuclides;

    /* ---- Tabulated extra energy density (mirrors PRIMAT.__init__'s
     * `extra_rho` list of rho(Tg) callables). Python cannot ship a live
     * callable across the C ABI, so
     * backend.py evaluates the *sum* of the user's extra_rho callables on a
     * dense log-spaced Tg grid once and hands the (Tg[], rho[]) arrays here;
     * cpr_bg_init_standard fits a cubic spline over them (in log10(Tg)) and
     * cpr_bg_Hubble adds rho(Tg) to rho_tot -- the exact place Python's
     * StandardBackground.Hubble sums self.extra_rho. Tg is [MeV], rho is
     * [MeV^4]. NULL/0 = no extra energy density (the common case). Owned
     * (malloc'd by the _wrapper.c bridge); freed by cpr_config_free. Not a
     * DEFAULT_PARAMS key and not routed through cpr_config_set_by_name -- it
     * is array-valued, set directly on the struct by the caller. */
    double *extra_rho_T;    /* Tg grid [MeV], strictly increasing, length extra_rho_n */
    double *extra_rho_val;  /* summed extra rho [MeV^4] at each extra_rho_T node */
    size_t  extra_rho_n;    /* number of grid points (>= 4 required by the cubic spline) */

    char data_dir[CPR_DATA_DIR_LEN]; /* the data folder itself (NEVO/, nuclear/, csv/, cache_plasma_weak/) */
} CPRConfig;

/* True iff cfg->network == "small" (mirrors network_is_small). Compares the
 * name, not the size: "small_parthenope" has the same twelve reactions and
 * answers 0. */
int cpr_config_network_is_small(const CPRConfig *cfg);

/* Recomputes everything that follows from cfg->consts: the cache-key hash
 * cfg->consts_hash, and the eta0b chain built on n0CMB/ma/maOvermB. Called
 * after every successful cpr_config_set_by_name of a measured constant, so no
 * path can leave a config whose derived values describe other constants.
 * Mirrors PRIMATConfig._update_constants. */
void cpr_config_refresh_constants(CPRConfig *cfg);

/* Derived constants depending on overridable params (mirrors the Python
 * @property of the same name). */
double cpr_config_Mpl(const CPRConfig *cfg);
double cpr_config_rhocOverh2(const CPRConfig *cfg);
double cpr_config_T_start_cosmo(const CPRConfig *cfg); /* [K] */
double cpr_config_T_end(const CPRConfig *cfg);         /* [K] */

/* Effective per-flavour neutrino chemical potentials: resolve the NAN
 * "inherit" sentinel of munuOverTnu_e/mu/tau back to the common munuOverTnu.
 * Mirror PRIMATConfig.xi_nu_e / xi_nu_mu / xi_nu_tau. Only xi_nu_e feeds the
 * n<->p weak rates; all three feed the neutrino energy density / Neff. */
double cpr_config_xi_nu_e(const CPRConfig *cfg);
double cpr_config_xi_nu_mu(const CPRConfig *cfg);
double cpr_config_xi_nu_tau(const CPRConfig *cfg);

/* Fills `cfg` with every DEFAULT_PARAMS value (string fields strdup'd so
 * the whole struct can later be freed uniformly by cpr_config_free).
 * `data_dir` is the data folder itself (e.g. .../primat/data, containing
 * NEVO/, nuclear/, csv/, cache_plasma_weak/) -- passed in rather than derived
 * from argv[0], since CPRIMAT supports --data-dir / the CPRIMAT_DATA_DIR
 * env var ahead of the executable-relative default -- see cli.c). Loads
 * nuclides.csv from `data_dir`/csv/. Returns 0 on success, nonzero (with
 * *errmsg set, caller frees) if nuclides.csv is missing or malformed. */
int cpr_config_init_defaults(CPRConfig *cfg, const char *data_dir, char **errmsg);

/* Resolves `relpath` (e.g. "nuclear/networks/large.txt" or
 * "nuclear/tables/<rxn>/<file>.txt") through the overlay chain:
 *   cfg->user_nuclear_dir (additive nuclear overlay, NULL = skip) ->
 *   cfg->data_dir + "/" + relpath (resolved default, tried last so
 *   shipped files are never unreachable when user_nuclear_dir is set).
 * Overlay roots for user_nuclear_dir are treated as the equivalent of
 * `primat/data/nuclear`: the resolver first tries
 * `base/<relpath without a leading "nuclear/">` and then the legacy nested
 * layout `base/<relpath>` for compatibility. The first candidate that
 * exists on disk wins; if none exist, the resolved-default path is written
 * anyway (so callers get a "missing file" error pointing at the expected
 * location). Writes into `out` (size `outsize`,
 * truncated/snprintf-safe like every other path builder in this codebase). */
void cpr_config_resolve_rates_path(const CPRConfig *cfg, const char *relpath,
                                    char *out, size_t outsize);

/* Cache-tree overlay (mirror of primat/cache_utils.py's
 * {cache_write_dir, resolve_cache_file}). `sub` is "weak" or "plasma".
 *
 * cpr_config_cache_write_dir: writes the WRITE directory
 *   <cache_dir>/<sub> if cache_dir is set, else <data_dir>/cache_plasma_weak/<sub>.
 * cpr_config_resolve_cache_file: writes the READ path for `file` through the
 *   overlay -- <cache_dir>/<sub>/<file> if it exists, else the shipped
 *   <data_dir>/cache_plasma_weak/<sub>/<file> if it exists, else the write
 *   path (where the file WILL be written). snprintf-safe like every other
 *   path builder here. Cache LOCATION only: never part of any fingerprint. */
void cpr_config_cache_write_dir(const CPRConfig *cfg, const char *sub,
                                char *out, size_t outsize);
void cpr_config_resolve_cache_file(const CPRConfig *cfg, const char *sub,
                                   const char *file, char *out, size_t outsize);

/* Sets cfg->Omegabh2_ and recomputes Omegabh2_to_eta0b/eta0b (the C
 * equivalent of the Python Omegabh2 property setter). */
void cpr_config_set_Omegabh2(CPRConfig *cfg, double value);
double cpr_config_get_Omegabh2(const CPRConfig *cfg);

/* GN is stored in cfg->GN in natural units [MeV^-2] (as consumed by
 * cpr_config_Mpl() and the Friedmann-equation Hubble helper), but exposed
 * to callers in SI units [m^3 kg^-1 s^-2] (matching
 * primat/config.py's DEFAULT_PARAMS["GN"]) -- always go through these two
 * functions rather than reading/writing cfg->GN directly. */
void cpr_config_set_GN(CPRConfig *cfg, double GN_SI);
double cpr_config_get_GN(const CPRConfig *cfg);

/* cpr_config_set_by_name return codes.
 *
 * The two failure modes need different handling and must not be conflated:
 *
 *   CPR_SET_UNKNOWN_KEY -- the name matches no DEFAULT_PARAMS key. A typo, or
 *       a key from a newer/older version. Python's PRIMATConfig warns and
 *       ignores it by default (strict_params=False), so the loaders do too.
 *   CPR_SET_BAD_VALUE   -- the key is known but the value has the wrong type
 *       (e.g. `network = 3`). Python raises TypeError unconditionally, so this
 *       is fatal on the C side as well: it is a malformed request, not a
 *       forward-compatibility question.
 *
 * A name with prefix "p_" or "delta_" goes into the corresponding CPRRxnMap
 * (value coerced to double). Booleans accept CPR_BOOL or CPR_INT (0/1,
 * mirroring Python's duck-typed bool/int interchangeability in
 * DEFAULT_PARAMS); numeric fields accept CPR_INT for CPR_DOUBLE (widened).
 * On any failure *errmsg is set (caller frees) and the target field is left
 * exactly as it was. */
#define CPR_SET_OK          0
#define CPR_SET_UNKNOWN_KEY 1
#define CPR_SET_BAD_VALUE   2

int cpr_config_set_by_name(CPRConfig *cfg, const char *name, CPRParam value,
                            char **errmsg);

/* ---- Parameter enumeration, for `cprimat --list-params` ----
 *
 * Walks the same internal field table cpr_config_set_by_name dispatches on
 * (plus the three keys routed around it: Omegabh2, GN, data_dir), in
 * DEFAULT_PARAMS order, so the listing cannot drift from what is actually
 * settable. Indices run 0 .. cpr_config_field_count()-1.
 *
 * cpr_config_format_value renders the value `cfg` currently holds for `name`
 * in the same spelling the INI/--set parser accepts back (True/False, an
 * integer, %g for reals, a bare string, None for an unset optional), so a
 * listed line can be pasted straight into an ini file. Returns 0 on success,
 * nonzero if `name` is not a known parameter. */
size_t cpr_config_field_count(void);
const char *cpr_config_field_name(size_t index);   /* NULL past the end */
int cpr_config_format_value(const CPRConfig *cfg, const char *name,
                            char *out, size_t outsize);

/* Validates flag-combination invariants (mirrors the `raise ValueError`
 * blocks in PyPRConfig.__init__, except the ones that require modules not
 * yet ported -- see config.c's top-of-function comment for the current
 * list). Returns 0 if valid, nonzero with *errmsg set (caller frees)
 * otherwise. Call once after all overrides (ini/cli/--set) are applied. */
int cpr_config_validate(CPRConfig *cfg, char **errmsg);

/* Frees every strdup'd string field, the nuclide table, and the two
 * CPRRxnMap dictionaries. Does not free `cfg` itself. */
void cpr_config_free(CPRConfig *cfg);

#endif /* CPRIMAT_CONFIG_H */
