#include "config.h"
#include "xalloc.h"
#include "cache.h"
#include "constants.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compat_posix.h"  /* unistd.h + sys/stat.h + strtok_r, portable */
#ifndef _WIN32
#include <pwd.h>           /* getpwuid/getpwnam for "~"/"~user" expansion  */
#endif

/* ===========================================================================
 * Literal parsing (--set KEY=VALUE / ini values), mirroring the
 * ast.literal_eval-equivalent used by primat/cli.py.
 * ===========================================================================
 */
CPRParam cpr_parse_literal(const char *s)
{
    CPRParam p;
    char *end;

    while (isspace((unsigned char)*s)) s++;
    size_t len = strlen(s);
    while (len > 0 && isspace((unsigned char)s[len - 1])) len--;

    /* Quoted string literal: strip matching quotes, return as-is. */
    if (len >= 2 && ((s[0] == '"' && s[len - 1] == '"') ||
                      (s[0] == '\'' && s[len - 1] == '\''))) {
        static char buf[1024];
        size_t n = len - 2 < sizeof(buf) - 1 ? len - 2 : sizeof(buf) - 1;
        memcpy(buf, s + 1, n);
        buf[n] = '\0';
        p.type = CPR_STRING;
        p.v.s = buf;
        return p;
    }

    if (strncasecmp(s, "none", len) == 0 && len == 4) {
        p.type = CPR_NONE;
        return p;
    }
    if (strncasecmp(s, "true", len) == 0 && len == 4) {
        p.type = CPR_BOOL;
        p.v.b = 1;
        return p;
    }
    if (strncasecmp(s, "false", len) == 0 && len == 5) {
        p.type = CPR_BOOL;
        p.v.b = 0;
        return p;
    }

    /* Try integer (must consume the whole trimmed token). */
    {
        char tmp[256];
        size_t n = len < sizeof(tmp) - 1 ? len : sizeof(tmp) - 1;
        memcpy(tmp, s, n);
        tmp[n] = '\0';
        long iv = strtol(tmp, &end, 10);
        if (end == tmp + n && n > 0) {
            p.type = CPR_INT;
            p.v.i = iv;
            return p;
        }
        double dv = strtod(tmp, &end);
        if (end == tmp + n && n > 0) {
            p.type = CPR_DOUBLE;
            p.v.d = dv;
            return p;
        }
    }

    /* Fall back to literal (unquoted) string. An *empty* token is not a valid
     * literal: Python's ast.literal_eval("") raises and cli.py turns that into
     * a parser.error, so `--set network=` must not quietly become the empty
     * string (which then fails much later, on a nonsensical ".txt" open).
     * CPR_NONE is not usable as the "invalid" signal (it is a legitimate
     * value), so this is reported as an empty *string* the callers reject --
     * see cli.c's --set handler and ini.c. */
    {
        static char buf[1024];
        size_t n = len < sizeof(buf) - 1 ? len : sizeof(buf) - 1;
        memcpy(buf, s, n);
        buf[n] = '\0';
        p.type = CPR_STRING;
        p.v.s = buf;
        return p;
    }
}

/* ===========================================================================
 * CPRParamList: retained, self-owning (key, value) overrides (see config.h).
 * ===========================================================================
 */
void cpr_paramlist_add(CPRParamList *pl, const char *key, CPRParam value)
{
    if (pl->n == pl->cap) {
        pl->cap = pl->cap ? pl->cap * 2 : 32;
        pl->items     = CPR_XREALLOC(pl->items, pl->cap * sizeof(*pl->items));
        pl->key_store = CPR_XREALLOC(pl->key_store, pl->cap * sizeof(*pl->key_store));
        pl->val_store = CPR_XREALLOC(pl->val_store, pl->cap * sizeof(*pl->val_store));
        /* realloc may have moved the two string arenas; every already-stored
         * entry points into them, so re-point them all before returning. */
        for (size_t i = 0; i < pl->n; i++) {
            pl->items[i].key = pl->key_store[i];
            if (pl->items[i].value.type == CPR_STRING)
                pl->items[i].value.v.s = pl->val_store[i];
        }
    }
    snprintf(pl->key_store[pl->n], CPR_PARAM_KEY_LEN, "%s", key);
    pl->items[pl->n].key = pl->key_store[pl->n];
    pl->items[pl->n].value = value;
    if (value.type == CPR_STRING) {
        /* Copy the string too: the caller's may be argv, an ini line buffer,
         * or cpr_parse_literal's static scratch -- none of which survive. */
        snprintf(pl->val_store[pl->n], CPR_PARAM_VAL_LEN, "%s", value.v.s ? value.v.s : "");
        pl->items[pl->n].value.v.s = pl->val_store[pl->n];
    }
    pl->n++;
}

void cpr_paramlist_free(CPRParamList *pl)
{
    free(pl->items);
    free(pl->key_store);
    free(pl->val_store);
    pl->items = NULL;
    pl->key_store = NULL;
    pl->val_store = NULL;
    pl->n = pl->cap = 0;
}

/* ===========================================================================
 * CPRRxnMap: p_<rxn> / delta_<rxn> dictionary.
 * ===========================================================================
 */
double cpr_rxnmap_get(const CPRRxnMap *map, const char *name)
{
    for (size_t i = 0; i < map->n; i++)
        if (strcmp(map->entries[i].name, name) == 0)
            return map->entries[i].value;
    return 0.0;
}

void cpr_rxnmap_set(CPRRxnMap *map, const char *name, double value)
{
    for (size_t i = 0; i < map->n; i++) {
        if (strcmp(map->entries[i].name, name) == 0) {
            map->entries[i].value = value;
            return;
        }
    }
    if (map->n == map->cap) {
        map->cap = map->cap ? map->cap * 2 : 64;
        map->entries = CPR_XREALLOC(map->entries, map->cap * sizeof(CPRRxnEntry));
    }
    strncpy(map->entries[map->n].name, name, sizeof(map->entries[map->n].name) - 1);
    map->entries[map->n].name[sizeof(map->entries[map->n].name) - 1] = '\0';
    map->entries[map->n].value = value;
    map->n++;
}

void cpr_rxnmap_free(CPRRxnMap *map)
{
    free(map->entries);
    map->entries = NULL;
    map->n = map->cap = 0;
}

/* ===========================================================================
 * nuclides.csv loader (mirrors PyPRConfig._load_nuclide_data).
 * ===========================================================================
 */
static int load_nuclides(CPRConfig *cfg, char **errmsg)
{
    char path[4200];
    snprintf(path, sizeof(path), "%s/csv/nuclides.csv", cfg->data_dir);

    FILE *f = fopen(path, "r");
    if (!f) {
        *errmsg = strdup("nuclides.csv not found (data_dir misconfigured?)");
        return 1;
    }

    char line[512];
    /* header: name,N,Z,A,Q,mass_excess_keV,spin -- locate columns by name so
     * a reordering of nuclides.csv doesn't silently break this loader. */
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        *errmsg = strdup("nuclides.csv is empty");
        return 1;
    }
    int col_name = -1, col_N = -1, col_Z = -1, col_mex = -1, col_spin = -1;
    {
        char hdr[512];
        strncpy(hdr, line, sizeof(hdr) - 1);
        hdr[sizeof(hdr) - 1] = '\0';
        int idx = 0;
        char *strtok_state = NULL;
        /* strtok_r, not strtok: mc.c's worker threads each call
         * cpr_config_init_defaults concurrently, and strtok keeps its
         * cursor in a single static buffer shared by every caller in the
         * process -- under threading that corrupts other threads'
         * in-progress parses (observed as spurious "header missing"/
         * dropped-row failures). strtok_r's state lives on this thread's
         * own stack instead. */
        for (char *tok = strtok_r(hdr, ",\r\n", &strtok_state); tok;
             tok = strtok_r(NULL, ",\r\n", &strtok_state), idx++) {
            if (strcmp(tok, "name") == 0) col_name = idx;
            else if (strcmp(tok, "N") == 0) col_N = idx;
            else if (strcmp(tok, "Z") == 0) col_Z = idx;
            else if (strcmp(tok, "mass_excess_keV") == 0) col_mex = idx;
            else if (strcmp(tok, "spin") == 0) col_spin = idx;
        }
    }
    if (col_name < 0 || col_N < 0 || col_Z < 0 || col_mex < 0 || col_spin < 0) {
        fclose(f);
        *errmsg = strdup("nuclides.csv header missing one of name,N,Z,mass_excess_keV,spin");
        return 1;
    }

    size_t cap = 64, n = 0;
    CPRNuclide *items = CPR_XMALLOC(cap * sizeof(CPRNuclide));
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        char row[512];
        strncpy(row, line, sizeof(row) - 1);
        row[sizeof(row) - 1] = '\0';
        char *fields[16] = {0};
        int nf = 0;
        char *strtok_state = NULL;
        for (char *tok = strtok_r(row, ",\r\n", &strtok_state); tok && nf < 16;
             tok = strtok_r(NULL, ",\r\n", &strtok_state))
            fields[nf++] = tok;
        if (nf <= col_name || nf <= col_N || nf <= col_Z || nf <= col_mex || nf <= col_spin)
            continue;

        if (n == cap) {
            cap *= 2;
            items = CPR_XREALLOC(items, cap * sizeof(CPRNuclide));
        }
        CPRNuclide *nuc = &items[n];
        strncpy(nuc->name, fields[col_name], sizeof(nuc->name) - 1);
        nuc->name[sizeof(nuc->name) - 1] = '\0';
        nuc->N = atoi(fields[col_N]);
        nuc->Z = atoi(fields[col_Z]);
        nuc->mass_excess_keV = atof(fields[col_mex]);
        nuc->spin = atof(fields[col_spin]);
        n++;
    }
    fclose(f);

    cfg->nuclides.items = items;
    cfg->nuclides.n = n;
    return 0;
}

/* ===========================================================================
 * Defaults + field table for generic name-based dispatch.
 * ===========================================================================
 */
/* F_DOUBLE_OR_NAN: a double whose "None" (Python) maps to NAN rather than 0.0,
 * used for the per-flavour munuOverTnu_e/mu/tau overrides where NAN is the
 * "inherit munuOverTnu" sentinel (0.0 is a legitimate value there, so it cannot
 * double as the sentinel the way F_DOUBLE_OR_NONE's 0.0 does for the MC cap). */
typedef enum { F_BOOL, F_INT, F_INT_OR_NONE, F_DOUBLE, F_DOUBLE_OR_NONE, F_DOUBLE_OR_NAN, F_STRING } FieldKind;

typedef struct {
    const char *name;
    FieldKind kind;
    size_t offset;
} FieldDesc;

#define FLD(field, kind) { #field, kind, offsetof(CPRConfig, field) }
/* The 16 measured physical constants live in the nested cfg->consts, but are
 * named without the prefix on every user-facing surface (params dict, --set,
 * ini), exactly as in Python's DEFAULT_PARAMS. */
#define FLD_CONST(field) { #field, F_DOUBLE, offsetof(CPRConfig, consts.field) }

static const FieldDesc FIELD_TABLE[] = {
    FLD(verbose, F_BOOL),
    FLD(show_progress, F_BOOL),
    FLD(debug, F_BOOL),
    FLD(numerical_precision, F_DOUBLE),
    FLD(numba_installed, F_BOOL),
    FLD(strict_params, F_BOOL),
    FLD(incomplete_decoupling, F_BOOL),
    FLD(QED_corrections, F_BOOL),
    FLD(n_electron_table, F_INT),
    FLD(recompute_electron_thermo, F_BOOL),
    FLD(recompute_qed_corrections, F_BOOL),
    FLD(spectral_distortions, F_BOOL),
    FLD(analytic_distortions, F_BOOL),
    FLD(y_SZ, F_DOUBLE),
    FLD(y_gray, F_DOUBLE),
    FLD(nevo_file, F_STRING),
    FLD(nevo_spectral_file, F_STRING),
    FLD(nevo_grid_file, F_STRING),
    FLD(nevo_file_prefix, F_STRING),
    FLD(external_scale_factor, F_BOOL),
    FLD(custom_background, F_STRING),
    /* The 16 measured physical constants (constants.OVERRIDABLE_CONSTANTS).
     * The exact ten are deliberately absent: they cannot be set by name on
     * either backend (Python's PRIMATConfig rejects them too). */
    FLD_CONST(alphaem),
    FLD_CONST(GF),
    FLD_CONST(mZ),
    FLD_CONST(me),
    FLD_CONST(mn),
    FLD_CONST(mp),
    FLD_CONST(T0CMB),
    FLD_CONST(gA),
    FLD_CONST(Vud),
    FLD_CONST(kappa_p),
    FLD_CONST(kappa_n),
    FLD_CONST(radproton),
    FLD_CONST(ma),
    FLD_CONST(He4Overma),
    FLD_CONST(HOverma),
    FLD_CONST(Neff_SM),
    FLD(T_start_cosmo_MeV, F_DOUBLE),
    FLD(T_end_MeV, F_DOUBLE),
    FLD(sampling_temperature_per_decade, F_INT),
    FLD(radiative_corrections, F_BOOL),
    FLD(finite_mass_corrections, F_BOOL),
    FLD(thermal_corrections, F_BOOL),
    FLD(weak_rate_cache, F_BOOL),
    FLD(save_nTOp, F_BOOL),
    FLD(sampling_nTOp_per_decade, F_INT),
    FLD(save_nTOp_thermal, F_BOOL),
    FLD(sampling_nTOp_thermal_per_decade, F_INT),
    FLD(tau_n_normalization, F_BOOL),
    FLD(tau_n, F_DOUBLE),
    FLD(std_tau_n, F_DOUBLE),
    FLD(vegas_n_eval, F_INT),
    FLD(vegas_n_itn, F_INT),
    FLD(epsrel_thermal, F_DOUBLE),
    FLD(output_time_evolution, F_BOOL),
    FLD(output_rates_time_evolution, F_BOOL),
    FLD(output_n_points, F_INT),
    FLD(output_file, F_STRING),
    FLD(output_final_result, F_BOOL),
    FLD(output_final_file, F_STRING),
    FLD(output_background_evolution, F_BOOL),
    FLD(output_background_file, F_STRING),
    FLD(output_mc_samples, F_BOOL),
    FLD(output_mc_covariance, F_BOOL),
    FLD(output_mc_correlation, F_BOOL),
    FLD(output_mc_file_prefix, F_STRING),
    FLD(rate_grid_npts, F_INT),
    FLD(rate_grid_T9_min, F_DOUBLE),
    FLD(rate_grid_T9_max, F_DOUBLE),
    FLD(network, F_STRING),
    FLD(amax, F_INT_OR_NONE),
    FLD(atol_large_LT, F_DOUBLE),
    FLD(rescale_nuclear_rates, F_BOOL),
    FLD(mc_rate_rescale_cap, F_DOUBLE_OR_NONE),
    FLD(nuclear_qed_corrections, F_BOOL),
    FLD(user_nuclear_dir, F_STRING),
    FLD(cache_dir, F_STRING),
    FLD(Omegach2, F_DOUBLE),
    FLD(h, F_DOUBLE),
    FLD(DeltaNeff, F_DOUBLE),
    FLD(munuOverTnu, F_DOUBLE),
    FLD(munuOverTnu_e, F_DOUBLE_OR_NAN),
    FLD(munuOverTnu_mu, F_DOUBLE_OR_NAN),
    FLD(munuOverTnu_tau, F_DOUBLE_OR_NAN),
    FLD(decay_reverse_rates, F_BOOL),
    FLD(decay_era, F_BOOL),
    FLD(t_decay_end, F_DOUBLE),
    FLD(decay_n_points, F_INT),
    FLD(output_decay_evolution, F_BOOL),
    FLD(output_decay_file, F_STRING),
    FLD(fEDE, F_DOUBLE),
    FLD(zcEDE, F_DOUBLE),
    FLD(wnEDE, F_DOUBLE),
    /* Omegabh2 deliberately absent: routed to cpr_config_set_Omegabh2()
     * by cpr_config_set_by_name() below, mirroring the Python @property. */
    /* GN deliberately absent: routed to cpr_config_set_GN() by
     * cpr_config_set_by_name() below -- see that function for why a raw
     * FLD() entry is unsafe. */
};
#define FIELD_TABLE_N (sizeof(FIELD_TABLE) / sizeof(FIELD_TABLE[0]))

/* The three DEFAULT_PARAMS keys cpr_config_set_by_name handles ahead of
 * FIELD_TABLE (they need a setter, not a raw field write). Listed here so
 * cpr_config_field_name() can enumerate the *complete* settable surface --
 * `cprimat --list-params` omitting them would be a lie about what --set
 * accepts. Order matches config.py's DEFAULT_PARAMS grouping. */
/* True iff a FIELD_TABLE offset points inside the nested cfg->consts, i.e. the
 * entry was declared with FLD_CONST and writing it invalidates consts_hash. */
static int offset_is_in_consts(size_t offset)
{
    return offset >= offsetof(CPRConfig, consts)
        && offset < offsetof(CPRConfig, consts) + sizeof(CPRConstants);
}

void cpr_config_refresh_constants(CPRConfig *cfg)
{
    cpr_constants_hash(&cfg->consts, cfg->consts_hash);
    /* eta0b is built from n0CMB, ma and maOvermB, so it must follow an
     * override of T0CMB/ma/He4Overma/HOverma -- it was computed once from the
     * defaults when Omegabh2 was first set. Mirrors _update_constants ->
     * _update_derived in primat/config.py. */
    cpr_config_set_Omegabh2(cfg, cfg->Omegabh2_);
}

static const char * const EXTRA_FIELD_NAMES[] = { "Omegabh2", "GN", "data_dir" };
#define EXTRA_FIELD_N (sizeof(EXTRA_FIELD_NAMES) / sizeof(EXTRA_FIELD_NAMES[0]))

size_t cpr_config_field_count(void) { return FIELD_TABLE_N + EXTRA_FIELD_N; }

const char *cpr_config_field_name(size_t index)
{
    if (index < FIELD_TABLE_N)
        return FIELD_TABLE[index].name;
    if (index < FIELD_TABLE_N + EXTRA_FIELD_N)
        return EXTRA_FIELD_NAMES[index - FIELD_TABLE_N];
    return NULL;
}

int cpr_config_format_value(const CPRConfig *cfg, const char *name,
                            char *out, size_t outsize)
{
    /* The three setter-routed keys first (they have no FIELD_TABLE entry). */
    if (strcmp(name, "Omegabh2") == 0) {
        snprintf(out, outsize, "%g", cpr_config_get_Omegabh2(cfg));
        return 0;
    }
    if (strcmp(name, "GN") == 0) {
        snprintf(out, outsize, "%g", cpr_config_get_GN(cfg));
        return 0;
    }
    if (strcmp(name, "data_dir") == 0) {
        snprintf(out, outsize, "%s", cfg->data_dir);
        return 0;
    }

    for (size_t i = 0; i < FIELD_TABLE_N; i++) {
        if (strcmp(FIELD_TABLE[i].name, name) != 0)
            continue;
        const void *field = (const char *)cfg + FIELD_TABLE[i].offset;
        switch (FIELD_TABLE[i].kind) {
        case F_BOOL:
            snprintf(out, outsize, "%s", *(const int *)field ? "True" : "False");
            return 0;
        case F_INT:
            snprintf(out, outsize, "%d", *(const int *)field);
            return 0;
        case F_INT_OR_NONE:
            /* -1 is the "None" sentinel (see config.h's header comment). */
            if (*(const int *)field < 0) snprintf(out, outsize, "None");
            else snprintf(out, outsize, "%d", *(const int *)field);
            return 0;
        case F_DOUBLE:
            snprintf(out, outsize, "%g", *(const double *)field);
            return 0;
        case F_DOUBLE_OR_NONE:
            /* 0.0 is the "no cap" sentinel (Python None). */
            if (*(const double *)field == 0.0) snprintf(out, outsize, "None");
            else snprintf(out, outsize, "%g", *(const double *)field);
            return 0;
        case F_DOUBLE_OR_NAN:
            /* NAN is the "inherit munuOverTnu" sentinel (Python None). */
            if (isnan(*(const double *)field)) snprintf(out, outsize, "None");
            else snprintf(out, outsize, "%g", *(const double *)field);
            return 0;
        case F_STRING:
            snprintf(out, outsize, "%s",
                     *(char * const *)field ? *(char * const *)field : "None");
            return 0;
        }
    }
    return 1;
}

static char *cpr_strdup(const char *s) { return s ? strdup(s) : NULL; }

static int cpr_is_path_field(const char *name)
{
    return strcmp(name, "nevo_file") == 0
        || strcmp(name, "nevo_spectral_file") == 0
        || strcmp(name, "nevo_grid_file") == 0
        || strcmp(name, "custom_background") == 0
        || strcmp(name, "user_nuclear_dir") == 0
        || strcmp(name, "cache_dir") == 0
        || strcmp(name, "output_file") == 0
        || strcmp(name, "output_final_file") == 0
        || strcmp(name, "output_background_file") == 0
        || strcmp(name, "output_mc_file_prefix") == 0
        || strcmp(name, "output_decay_file") == 0;
}

static char *cpr_expanduser_path(const char *path)
{
    if (!path)
        return NULL;
    if (path[0] != '~')
        return strdup(path);

    const char *rest = NULL;
    const char *home = NULL;
    if (path[1] == '\0' || path[1] == '/') {
        /* "~" and "~/" both expand against the current user's home dir. */
        home = getenv("HOME");
#ifdef _WIN32
        /* Windows has no HOME/passwd; fall back to the standard profile
         * environment variable (%USERPROFILE%, e.g. C:\Users\alice). */
        if (!home || !home[0])
            home = getenv("USERPROFILE");
#else
        if (!home || !home[0]) {
            struct passwd *pw = getpwuid(getuid());
            if (pw && pw->pw_dir && pw->pw_dir[0])
                home = pw->pw_dir;
        }
#endif
        rest = path + 1;
    } else {
#ifdef _WIN32
        /* "~user" (another account's home) is not resolvable without the
         * Windows user-profile lookup APIs; leave such paths unexpanded. */
        return strdup(path);
#else
        /* "~user/..." expands against that account's passwd entry. */
        const char *slash = strchr(path + 1, '/');
        size_t user_len = slash ? (size_t)(slash - (path + 1)) : strlen(path + 1);
        char user[256];
        if (user_len == 0 || user_len >= sizeof(user))
            return strdup(path);
        memcpy(user, path + 1, user_len);
        user[user_len] = '\0';
        struct passwd *pw = getpwnam(user);
        if (pw && pw->pw_dir && pw->pw_dir[0]) {
            home = pw->pw_dir;
            rest = slash ? slash : "";
        } else {
            return strdup(path);
        }
#endif
    }

    if (!home || !home[0])
        return strdup(path);

    size_t home_len = strlen(home);
    size_t rest_len = strlen(rest);
    char *out = malloc(home_len + rest_len + 1);
    if (!out)
        return NULL;
    memcpy(out, home, home_len);
    memcpy(out + home_len, rest, rest_len + 1);
    return out;
}

/* Assign cfg->data_dir, expanding a leading "~" and warning on truncation.
 *
 * cfg->data_dir is a fixed CPR_DATA_DIR_LEN buffer (not a malloc'd char*,
 * unlike user_nuclear_dir/cache_dir), so an absurdly long
 * --data_dir/CPRIMAT_DATA_DIR value must be truncated, not overflowed.
 * strncpy alone does not guarantee NUL-termination when the source is >= the
 * copy length, so we terminate explicitly rather than relying on any prior
 * memset; the warning is unconditional (not gated by cfg->verbose, unlike
 * cpr_log) since every subsequent rates lookup under a truncated data_dir
 * would otherwise fail with a confusing "file not found" pointing at a
 * mangled path.
 *
 * Shared by cpr_config_init_defaults (the --data_dir / CPRIMAT_DATA_DIR /
 * Python-extension argument) and cpr_config_set_by_name's "data_dir" case (an
 * INI key or a params-dict entry), so every path normalises identically --
 * including the "~" expansion the Python side applies via _PATH_PARAMS. */
static void cpr_config_assign_data_dir(CPRConfig *cfg, const char *data_dir)
{
    char *expanded = cpr_expanduser_path(data_dir);
    const char *src = expanded ? expanded : data_dir;   /* OOM -> use as-is */
    strncpy(cfg->data_dir, src, sizeof(cfg->data_dir) - 1);
    cfg->data_dir[sizeof(cfg->data_dir) - 1] = '\0';
    if (strlen(src) >= sizeof(cfg->data_dir)) {
        fprintf(stderr,
                "warning: data_dir path (%zu bytes) exceeds the %zu-byte "
                "internal limit and was truncated to %.60s...; rates lookups "
                "will very likely fail. Use a shorter data_dir path.\n",
                strlen(src), sizeof(cfg->data_dir) - 1, cfg->data_dir);
    }
    free(expanded);
}

int cpr_config_init_defaults(CPRConfig *cfg, const char *data_dir, char **errmsg)
{
    memset(cfg, 0, sizeof(*cfg));
    cpr_constants_init();
    /* This run's own copy of the physical constants: the 16 measured ones are
     * then settable by name (see FLD_CONST in FIELD_TABLE), the ten exact ones
     * stay at these values. */
    cfg->consts = g_const;
    cpr_constants_hash(&cfg->consts, cfg->consts_hash);

    cpr_config_assign_data_dir(cfg, data_dir);

    cfg->verbose = 0;
    cfg->show_progress = 1;
    cfg->debug = 0;
    cfg->numerical_precision = 1.e-7;
    cfg->numba_installed = 1;
    cfg->strict_params = 0;

    cfg->incomplete_decoupling = 1;

    cfg->QED_corrections = 1;
    cfg->n_electron_table = 2000;
    cfg->recompute_electron_thermo = 0;
    cfg->recompute_qed_corrections = 0;

    cfg->spectral_distortions = 1;
    cfg->analytic_distortions = 0;
    cfg->y_SZ = 0.;
    cfg->y_gray = 0.;

    cfg->nevo_file = NULL;
    cfg->nevo_spectral_file = NULL;
    cfg->nevo_grid_file = NULL;
    cfg->nevo_file_prefix = cpr_strdup("NEVOPRIMAT");

    cfg->external_scale_factor = 0;
    cfg->custom_background = NULL;

    /* SI-unit default (m^3 kg^-1 s^-2), matching primat/config.py's
     * DEFAULT_PARAMS["GN"] exactly (the CODATA-tabulated 5-significant-figure
     * literal -- do not replace with a natural-units round-trip result, see
     * that file's comment); converted to the natural units (MeV^-2)
     * cfg->GN is stored in by cpr_config_set_GN(). */
    cpr_config_set_GN(cfg, 6.6743e-11);

    cfg->T_start_cosmo_MeV = 40.0;
    cfg->T_end_MeV = 1.e-3;
    cfg->sampling_temperature_per_decade = 600;

    cfg->radiative_corrections = 1;
    cfg->finite_mass_corrections = 1;
    cfg->thermal_corrections = 1;
    cfg->weak_rate_cache = 1;
    cfg->save_nTOp = 1;
    cfg->sampling_nTOp_per_decade = 80;
    cfg->save_nTOp_thermal = 1;
    cfg->sampling_nTOp_thermal_per_decade = 20;
    cfg->tau_n_normalization = 1;
    cfg->tau_n = 878.4;
    cfg->std_tau_n = 0.5;
    cfg->vegas_n_eval = 20000;
    cfg->vegas_n_itn = 20;
    cfg->epsrel_thermal = 1.e-2;

    cfg->output_time_evolution = 0;
    cfg->output_rates_time_evolution = 0;
    cfg->output_n_points = 500;
    cfg->output_file = cpr_strdup("results/output_tables.tsv");
    cfg->output_final_result = 0;
    cfg->output_final_file = cpr_strdup("results/output_final.dat");
    cfg->output_background_evolution = 0;
    cfg->output_background_file = cpr_strdup("results/output_background.tsv");
    cfg->output_mc_samples = 0;
    cfg->output_mc_covariance = 0;
    cfg->output_mc_correlation = 0;
    cfg->output_mc_file_prefix = cpr_strdup("results/output_mc");

    cfg->rate_grid_npts = 1000;
    cfg->rate_grid_T9_min = 1.0e-3;
    cfg->rate_grid_T9_max = 10.0;
    cfg->network = cpr_strdup("small");
    cfg->amax = -1; /* None */
    cfg->atol_large_LT = 1.e-26;
    cfg->rescale_nuclear_rates = 0;
    cfg->mc_rate_rescale_cap = 30.0; /* 0.0 = no cap (mirrors Python None) */
    cfg->nuclear_qed_corrections = 1;
    cfg->user_nuclear_dir = NULL;
    cfg->cache_dir = NULL;

    cfg->Omegabh2_ = 0.02242;
    cfg->Omegach2 = 0.11933;
    cfg->h = 0.6766;
    cfg->DeltaNeff = 0.;
    cfg->munuOverTnu = 0.;
    /* NAN = "inherit munuOverTnu" (mirrors Python None); resolved by
     * cpr_config_xi_nu_e/mu/tau(). */
    cfg->munuOverTnu_e = NAN;
    cfg->munuOverTnu_mu = NAN;
    cfg->munuOverTnu_tau = NAN;

    cfg->decay_reverse_rates = 0;
    cfg->decay_era = 0;
    cfg->t_decay_end = 3.156e16;
    cfg->decay_n_points = 200;
    cfg->output_decay_evolution = 0;
    cfg->output_decay_file = cpr_strdup("results/output_decay_evolution.tsv");

    cfg->fEDE = 0.;
    cfg->zcEDE = 1.e8;
    cfg->wnEDE = 1.;

    /* Tabulated extra_rho: empty by default (set directly by the caller when
     * the Python `extra_rho` constructor argument was used -- see config.h). */
    cfg->extra_rho_T = NULL;
    cfg->extra_rho_val = NULL;
    cfg->extra_rho_n = 0;

    if (load_nuclides(cfg, errmsg))
        return 1;

    /* Omegabh2_to_eta0b / eta0b depend on Omegabh2_, set just above. */
    cpr_config_set_Omegabh2(cfg, cfg->Omegabh2_);
    return 0;
}

int cpr_config_is_small(const CPRConfig *cfg) { return strcmp(cfg->network, "small") == 0; }
int cpr_config_is_large(const CPRConfig *cfg) { return strcmp(cfg->network, "large") == 0; }

double cpr_config_Mpl(const CPRConfig *cfg) { return 1. / sqrt(cfg->GN); }

double cpr_config_rhocOverh2(const CPRConfig *cfg)
{
    double H = cpr_HubbleOverh();
    return 3. / (8. * M_PI * cfg->GN) * H * H;
}

double cpr_config_T_start_cosmo(const CPRConfig *cfg)
{
    return cfg->T_start_cosmo_MeV * cpr_MeV_to_Kelvin();
}

double cpr_config_T_end(const CPRConfig *cfg)
{
    return cfg->T_end_MeV * cpr_MeV_to_Kelvin();
}

/* Effective per-flavour ξ: NAN in the override field means "inherit the
 * common munuOverTnu"; any finite value is that flavour's own ξ. Mirrors
 * PRIMATConfig.xi_nu_e / xi_nu_mu / xi_nu_tau (None -> munuOverTnu). */
double cpr_config_xi_nu_e(const CPRConfig *cfg)
{
    return isnan(cfg->munuOverTnu_e) ? cfg->munuOverTnu : cfg->munuOverTnu_e;
}
double cpr_config_xi_nu_mu(const CPRConfig *cfg)
{
    return isnan(cfg->munuOverTnu_mu) ? cfg->munuOverTnu : cfg->munuOverTnu_mu;
}
double cpr_config_xi_nu_tau(const CPRConfig *cfg)
{
    return isnan(cfg->munuOverTnu_tau) ? cfg->munuOverTnu : cfg->munuOverTnu_tau;
}

static int path_exists(const char *path)
{
    struct stat st;
    return stat(path, &st) == 0;
}

/* snprintf wrapper shared by every path-building helper below, so
 * truncation is detected in exactly one place instead of after each of the
 * ~8 individual snprintf calls that join a (possibly attacker/user-supplied,
 * unbounded-length) data_dir/user_nuclear_dir/cache_dir with a relative
 * path. snprintf itself is always memory-safe (it never writes past
 * outsize and always NUL-terminates when outsize > 0), so a truncated
 * result cannot overflow -- but it silently becomes a syntactically valid,
 * *wrong* path, which then fails a later path_exists()/fopen() with a
 * confusing "file not found" instead of the real "path too long" cause.
 * `what` names the field/call site for that warning; printed unconditionally
 * (unlike cpr_log, not gated by cfg->verbose) since a truncated data path
 * is a misconfiguration a user needs to see regardless. */
static int cpr_snprintf_path(char *out, size_t outsize, const char *what,
                              const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    int n = vsnprintf(out, outsize, fmt, args);
    va_end(args);
    if (n >= 0 && (size_t)n >= outsize) {
        fprintf(stderr,
                "warning: %s path was truncated (%d bytes needed, only a "
                "%zu-byte buffer available); the resulting path is very "
                "likely wrong. Use a shorter data_dir/user_nuclear_dir/"
                "cache_dir.\n", what, n, outsize);
    }
    return n;
}

void cpr_config_resolve_rates_path(const CPRConfig *cfg, const char *relpath,
                                    char *out, size_t outsize)
{
    char candidate[CPR_PATH_BUF_LEN2];

    if (cfg->user_nuclear_dir) {
        /* user_nuclear_dir is the equivalent of primat/data/nuclear, so strip
         * any leading "nuclear/" component before joining, then fall through
         * to the legacy nested layout for compatibility. */
        if (strncmp(relpath, "nuclear/", 8) == 0) {
            cpr_snprintf_path(candidate, sizeof(candidate), "user_nuclear_dir",
                               "%s/%s", cfg->user_nuclear_dir, relpath + 8);
            if (path_exists(candidate)) {
                cpr_snprintf_path(out, outsize, "resolved rates", "%s", candidate);
                return;
            }
        }
        cpr_snprintf_path(candidate, sizeof(candidate), "user_nuclear_dir",
                           "%s/%s", cfg->user_nuclear_dir, relpath);
        if (path_exists(candidate)) {
            cpr_snprintf_path(out, outsize, "resolved rates", "%s", candidate);
            return;
        }
    }
    /* Resolved default (cfg->data_dir, set by cpr_config_init_defaults from
     * the --data-dir flag / CPRIMAT_DATA_DIR env var / the Python backend's
     * cfg._resolved_data_dir), always tried last (and returned even when
     * missing, so the caller's "file not found" error points at the expected
     * location). cfg->data_dir is the data folder itself
     * (e.g. .../primat/data), not its parent. */
    cpr_snprintf_path(out, outsize, "data_dir", "%s/%s", cfg->data_dir, relpath);
}

/* Cache-tree overlay -- mirror of primat/cache_utils.py's
 * {cache_write_dir, resolve_cache_file}. Cache LOCATION only, never part of
 * any fingerprint. `sub` is "weak" or "plasma". */
void cpr_config_cache_write_dir(const CPRConfig *cfg, const char *sub,
                                char *out, size_t outsize)
{
    if (cfg->cache_dir && cfg->cache_dir[0])
        cpr_snprintf_path(out, outsize, "cache_dir", "%s/%s", cfg->cache_dir, sub);
    else
        cpr_snprintf_path(out, outsize, "data_dir", "%s/cache_plasma_weak/%s", cfg->data_dir, sub);
}

void cpr_config_resolve_cache_file(const CPRConfig *cfg, const char *sub,
                                   const char *file, char *out, size_t outsize)
{
    char cand[CPR_PATH_BUF_LEN2];
    if (cfg->cache_dir && cfg->cache_dir[0]) {          /* overlay: redirect first */
        cpr_snprintf_path(cand, sizeof(cand), "cache_dir", "%s/%s/%s", cfg->cache_dir, sub, file);
        if (path_exists(cand)) { cpr_snprintf_path(out, outsize, "resolved cache", "%s", cand); return; }
    }
    /* shipped package copy (always tried, always last) */
    cpr_snprintf_path(cand, sizeof(cand), "data_dir", "%s/cache_plasma_weak/%s/%s",
             cfg->data_dir, sub, file);
    if (path_exists(cand)) { cpr_snprintf_path(out, outsize, "resolved cache", "%s", cand); return; }
    /* miss -> the write path (where it WILL be written) */
    cpr_config_cache_write_dir(cfg, sub, out, outsize);
    size_t n = strlen(out);  /* always < outsize: snprintf NUL-terminates within outsize */
    cpr_snprintf_path(out + n, outsize - n, "cache write path", "/%s", file);
}

void cpr_config_set_Omegabh2(CPRConfig *cfg, double value)
{
    cfg->Omegabh2_ = value;
    /* Omegabh2_to_eta0b = (rhocOverh2 / n0CMB) / (ma / maOvermB); eta0b =
     * Omegabh2_to_eta0b * Omegabh2 (Phys. Rep. baryon-to-photon ratio). */
    cfg->Omegabh2_to_eta0b = (cpr_config_rhocOverh2(cfg) / cpr_n0CMB(&cfg->consts))
                             / (cfg->consts.ma / cpr_maOvermB(&cfg->consts));
    cfg->eta0b = cfg->Omegabh2_to_eta0b * cfg->Omegabh2_;
}

double cpr_config_get_Omegabh2(const CPRConfig *cfg) { return cfg->Omegabh2_; }

/* GN_SI_to_MeV2: conversion factor from Newton's constant in SI units
 * [m^3 kg^-1 s^-2] to the natural units [MeV^-2] that cfg->GN is stored in
 * internally (used by cpr_config_Mpl / the Friedmann-equation Hubble
 * helper). Mirrors primat/constants.py's CONST.GN_SI_to_MeV2 exactly:
 * G_natural[MeV^-2] = G_SI[m^3 kg^-1 s^-2] / (hbar[erg s] * clight[cm/s]^5
 * / MeV[erg]^2 * 1e-3), the 1e-3 converting cm^3 g^-1 s^-2 to
 * m^3 kg^-1 s^-2. */
static double GN_SI_to_MeV2(void)
{
    double GN_MeV2_to_SI = g_const.hbar * pow(g_const.clight, 5)
                            / (g_const.MeV * g_const.MeV) * 1e-3;
    return 1. / GN_MeV2_to_SI;
}

void cpr_config_set_GN(CPRConfig *cfg, double GN_SI)
{
    /* cfg->GN (natural units, MeV^-2) is what cpr_config_Mpl() and the
     * Friedmann-equation Hubble helper actually consume; GN_SI (the
     * user-facing value, matching DEFAULT_PARAMS["GN"] and the GUI's
     * "Constants" panel) must be converted before storing, exactly like
     * the Python-only PRIMATConfig._GN_MeV2 property does. Without this
     * conversion, an SI-unit GN (~6.674e-11) landing in a field meant for
     * natural units (~6.709e-45) is off by ~34 orders of magnitude and
     * produces a meaningless Hubble rate. */
    cfg->GN = GN_SI * GN_SI_to_MeV2();
}

double cpr_config_get_GN(const CPRConfig *cfg)
{
    return cfg->GN / GN_SI_to_MeV2();
}

int cpr_config_set_by_name(CPRConfig *cfg, const char *name, CPRParam value,
                            char **errmsg)
{
    if (strncmp(name, "p_", 2) == 0) {
        double d = value.type == CPR_DOUBLE ? value.v.d
                 : value.type == CPR_INT ? (double)value.v.i
                 : value.type == CPR_BOOL ? (double)value.v.b : 0.0;
        cpr_rxnmap_set(&cfg->p_rxn, name + 2, d);
        return 0;
    }
    if (strncmp(name, "delta_", 6) == 0) {
        double d = value.type == CPR_DOUBLE ? value.v.d
                 : value.type == CPR_INT ? (double)value.v.i
                 : value.type == CPR_BOOL ? (double)value.v.b : 0.0;
        cpr_rxnmap_set(&cfg->delta_rxn, name + 6, d);
        return 0;
    }
    if (strcmp(name, "Omegabh2") == 0) {
        double d = value.type == CPR_DOUBLE ? value.v.d
                 : value.type == CPR_INT ? (double)value.v.i : NAN;
        if (isnan(d)) {
            *errmsg = strdup("Omegabh2 requires a numeric value");
            return CPR_SET_BAD_VALUE;
        }
        cpr_config_set_Omegabh2(cfg, d);
        return CPR_SET_OK;
    }
    if (strcmp(name, "GN") == 0) {
        double d = value.type == CPR_DOUBLE ? value.v.d
                 : value.type == CPR_INT ? (double)value.v.i : NAN;
        if (isnan(d)) {
            *errmsg = strdup("GN requires a numeric value");
            return CPR_SET_BAD_VALUE;
        }
        cpr_config_set_GN(cfg, d);
        return CPR_SET_OK;
    }
    /* data_dir is routed here rather than through FIELD_TABLE because it is a
     * fixed char[] buffer, and F_STRING's handler free()s and strdup()s a
     * char* -- which would corrupt it. It must still round-trip like every
     * other key: primat/backend.py forwards the whole params dict to the
     * extension (data_dir included, on top of the positional argument), and
     * primat-c/examples/run_basic.ini advertises `data_dir` as a settable INI
     * key. Rejecting it here used to abort the Python C-backend call outright
     * and made the INI key a warn-and-ignore no-op.
     *
     * None means "keep the current root" (the Python default, where None
     * selects the shipped primat/data/ tree). Because the nuclide table was
     * already loaded from the previous root by cpr_config_init_defaults, a
     * *changed* root must reload it -- otherwise an INI-supplied data_dir
     * would silently keep the shipped nuclides.csv. */
    if (strcmp(name, "data_dir") == 0) {
        if (value.type == CPR_NONE)
            return CPR_SET_OK;
        if (value.type != CPR_STRING) {
            *errmsg = strdup("data_dir expects a string or None");
            return CPR_SET_BAD_VALUE;
        }
        char previous[CPR_DATA_DIR_LEN];
        snprintf(previous, sizeof(previous), "%s", cfg->data_dir);
        cpr_config_assign_data_dir(cfg, value.v.s);
        if (strcmp(previous, cfg->data_dir) != 0) {
            free(cfg->nuclides.items);
            cfg->nuclides.items = NULL;
            cfg->nuclides.n = 0;
            if (load_nuclides(cfg, errmsg))
                return CPR_SET_BAD_VALUE;
        }
        return CPR_SET_OK;
    }

    for (size_t i = 0; i < FIELD_TABLE_N; i++) {
        if (strcmp(FIELD_TABLE[i].name, name) != 0)
            continue;
        void *field = (char *)cfg + FIELD_TABLE[i].offset;
        switch (FIELD_TABLE[i].kind) {
        case F_BOOL:
            if (value.type != CPR_BOOL && value.type != CPR_INT) {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects a bool", name);
                return CPR_SET_BAD_VALUE;
            }
            *(int *)field = value.type == CPR_BOOL ? value.v.b : (int)value.v.i;
            return CPR_SET_OK;
        case F_INT:
            if (value.type != CPR_INT && value.type != CPR_BOOL) {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects an int", name);
                return CPR_SET_BAD_VALUE;
            }
            *(int *)field = value.type == CPR_INT ? (int)value.v.i : value.v.b;
            return CPR_SET_OK;
        case F_INT_OR_NONE:
            if (value.type == CPR_NONE) {
                *(int *)field = -1;
                return CPR_SET_OK;
            }
            if (value.type != CPR_INT) {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects an int or None", name);
                return CPR_SET_BAD_VALUE;
            }
            *(int *)field = (int)value.v.i;
            return CPR_SET_OK;
        case F_DOUBLE:
            if (value.type == CPR_DOUBLE) *(double *)field = value.v.d;
            else if (value.type == CPR_INT) *(double *)field = (double)value.v.i;
            else {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects a number", name);
                return CPR_SET_BAD_VALUE;
            }
            /* One of the 16 measured constants (an FLD_CONST entry): keep the
             * cached hash in step, so the fingerprints can never describe
             * another set of constants than the one just set. */
            if (offset_is_in_consts(FIELD_TABLE[i].offset))
                cpr_config_refresh_constants(cfg);
            return CPR_SET_OK;
        case F_DOUBLE_OR_NONE:
            /* None → 0.0 (sentinel for "no cap"); any positive number is the cap value. */
            if (value.type == CPR_NONE) {
                *(double *)field = 0.0;
                return CPR_SET_OK;
            }
            if (value.type == CPR_DOUBLE) *(double *)field = value.v.d;
            else if (value.type == CPR_INT) *(double *)field = (double)value.v.i;
            else {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects a number or None", name);
                return CPR_SET_BAD_VALUE;
            }
            return CPR_SET_OK;
        case F_DOUBLE_OR_NAN:
            /* None → NAN (sentinel for "inherit munuOverTnu"); any number is a
             * concrete per-flavour ξ (may be negative). */
            if (value.type == CPR_NONE) {
                *(double *)field = NAN;
                return CPR_SET_OK;
            }
            if (value.type == CPR_DOUBLE) *(double *)field = value.v.d;
            else if (value.type == CPR_INT) *(double *)field = (double)value.v.i;
            else {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects a number or None", name);
                return CPR_SET_BAD_VALUE;
            }
            return CPR_SET_OK;
        case F_STRING: {
            /* Build the replacement FIRST, and only then release the old
             * value. The obvious ordering (free, then switch on the type)
             * leaves the field holding a freed pointer whenever the value is
             * not a string -- and since ini.c/cli.c used to continue past that
             * error, the run went on to *read* it (a garbage network path) and
             * cpr_config_free double-freed it at exit. Nothing here touches
             * `field` until the new value is in hand. */
            char *newval = NULL;
            if (value.type == CPR_STRING) {
                newval = cpr_is_path_field(name)
                    ? cpr_expanduser_path(value.v.s)
                    : strdup(value.v.s);
                if (!newval) {
                    *errmsg = strdup("out of memory while copying string parameter");
                    return CPR_SET_BAD_VALUE;
                }
            } else if (value.type != CPR_NONE) {
                *errmsg = malloc(128);
                snprintf(*errmsg, 128, "%s expects a string or None", name);
                return CPR_SET_BAD_VALUE;   /* field untouched */
            }
            free(*(char **)field);
            *(char **)field = newval;       /* NULL for CPR_NONE */
            return CPR_SET_OK;
        }
        }
    }

    *errmsg = malloc(256);
    snprintf(*errmsg, 256, "unknown parameter key: %s", name);
    return CPR_SET_UNKNOWN_KEY;
}

int cpr_config_validate(CPRConfig *cfg, char **errmsg)
{
    /* custom_background: force instantaneous decoupling / no spectral
     * distortions, mirroring PRIMATConfig.__init__. The forcing is reported on
     * stderr, as the Python side reports it with warnings.warn: silently
     * overriding two flags the user explicitly set leaves them believing the
     * run included non-instantaneous decoupling / spectral distortions when it
     * did not. Printed unconditionally (not via cpr_log, which is gated on
     * cfg->verbose) for the same reason the Python warning is not gated. */
    if (cfg->custom_background != NULL) {
        if (cfg->external_scale_factor) {
            *errmsg = strdup("custom_background and external_scale_factor are mutually exclusive");
            return 1;
        }
        if (cfg->incomplete_decoupling || cfg->spectral_distortions) {
            fprintf(stderr,
                    "warning: custom_background: forcing %s%s%s "
                    "(custom-background mode uses instantaneous-decoupling "
                    "weak rates; spectral distortions are not supported).\n",
                    cfg->incomplete_decoupling ? "incomplete_decoupling=False" : "",
                    (cfg->incomplete_decoupling && cfg->spectral_distortions) ? ", " : "",
                    cfg->spectral_distortions ? "spectral_distortions=False" : "");
        }
        cfg->incomplete_decoupling = 0;
        cfg->spectral_distortions = 0;
    }

    /* NOTE: the network-file-existence check (PRIMATConfig.__init__'s
     * `network must be 'small' or name an existing file...`) and the
     * p_<rxn>/delta_<rxn> typo check against the configured network's
     * reaction list both need network_data.c to enumerate valid reaction
     * names, so they are not performed here. Callers reaching the C solver
     * through primat/backend.py get both checks anyway (it builds a
     * PRIMATConfig from the same params dict first); the standalone
     * primat-c CLI instead reports a missing network later, when the list is
     * opened ("error: cannot open network list '<path>'"). */

    if (cfg->amax != -1 && cfg->amax < 1) {
        *errmsg = strdup("amax must be None (-1) or a positive integer");
        return 1;
    }

    /* user_nuclear_dir must name an existing directory, mirroring
     * PRIMATConfig._validate_dir_field ("user_nuclear_dir=... is not an
     * existing directory", a ValueError). Checked HERE rather than only on
     * cli.c's --user_nuclear_dir flag, so the same typo is caught however the
     * key arrives -- an ini line, --set, or the Python extension's params
     * dict. Without it the overlay silently resolves nothing and every rate
     * table falls back to the shipped tree: plausible numbers computed from
     * the wrong rate set, with no error anywhere.
     * cache_dir is deliberately NOT checked (on either side): it is a write
     * target created on demand, so "does not exist yet" is its normal state. */
    if (cfg->user_nuclear_dir && cfg->user_nuclear_dir[0]) {
        struct stat st;
        if (stat(cfg->user_nuclear_dir, &st) != 0 || !S_ISDIR(st.st_mode)) {
            *errmsg = malloc(CPR_PARAM_VAL_LEN);
            snprintf(*errmsg, CPR_PARAM_VAL_LEN,
                     "user_nuclear_dir='%.700s' is not an existing directory",
                     cfg->user_nuclear_dir);
            return 1;
        }
    }

    /* Physical/numerical range checks, mirroring primat/config.py's
     * _PARAM_RANGE table, for clearer config-validation error messages. Each
     * guards a value that
     * must stay strictly positive (a physical scale, tolerance, count, or
     * time) or non-negative. Emitted here rather than in cpr_config_set_by_name
     * so a value set via any path (INI, --set, wrapper) is caught uniformly. */
#define CPR_REQUIRE(cond, fmt, val) \
    do { if (!(cond)) { *errmsg = malloc(160); \
        snprintf(*errmsg, 160, fmt, val); return 1; } } while (0)

    /* strictly positive doubles */
    CPR_REQUIRE(cpr_config_get_Omegabh2(cfg) > 0,
                "Omegabh2=%.6g is out of range: must be > 0", cpr_config_get_Omegabh2(cfg));
    CPR_REQUIRE(cpr_config_get_GN(cfg) > 0,
                "GN=%.6g is out of range: must be > 0", cpr_config_get_GN(cfg));
    CPR_REQUIRE(cfg->numerical_precision > 0,
                "numerical_precision=%.6g is out of range: must be > 0", cfg->numerical_precision);
    CPR_REQUIRE(cfg->T_start_cosmo_MeV > 0,
                "T_start_cosmo_MeV=%.6g is out of range: must be > 0", cfg->T_start_cosmo_MeV);
    CPR_REQUIRE(cfg->T_end_MeV > 0,
                "T_end_MeV=%.6g is out of range: must be > 0", cfg->T_end_MeV);
    CPR_REQUIRE(cfg->tau_n > 0,
                "tau_n=%.6g is out of range: must be > 0", cfg->tau_n);
    CPR_REQUIRE(cfg->Omegach2 > 0,
                "Omegach2=%.6g is out of range: must be > 0", cfg->Omegach2);
    CPR_REQUIRE(cfg->h > 0,
                "h=%.6g is out of range: must be > 0", cfg->h);
    CPR_REQUIRE(cfg->atol_large_LT > 0,
                "atol_large_LT=%.6g is out of range: must be > 0", cfg->atol_large_LT);
    CPR_REQUIRE(cfg->epsrel_thermal > 0,
                "epsrel_thermal=%.6g is out of range: must be > 0", cfg->epsrel_thermal);
    CPR_REQUIRE(cfg->t_decay_end > 0,
                "t_decay_end=%.6g is out of range: must be > 0", cfg->t_decay_end);
    CPR_REQUIRE(cfg->zcEDE > 0,
                "zcEDE=%.6g is out of range: must be > 0", cfg->zcEDE);
    CPR_REQUIRE(cfg->rate_grid_T9_min > 0,
                "rate_grid_T9_min=%.6g is out of range: must be > 0", cfg->rate_grid_T9_min);
    CPR_REQUIRE(cfg->rate_grid_T9_max > 0,
                "rate_grid_T9_max=%.6g is out of range: must be > 0", cfg->rate_grid_T9_max);
    /* mc_rate_rescale_cap: 0.0 is the "no cap" sentinel (Python None); any
     * value >= 1 is the cap. Below 1 the clamp bounds [1/cap, cap] *cross*
     * (cap = 0.5 -> [2, 0.5]), which silently pins every sampled rate factor
     * to cap instead of capping the variation -- see network_data.c's
     * rate-variation clamp and _validate_ranges in primat/config.py. */
    CPR_REQUIRE(cfg->mc_rate_rescale_cap == 0.0 || cfg->mc_rate_rescale_cap >= 1.0,
                "mc_rate_rescale_cap=%.6g is out of range: must be >= 1 or None "
                "(it clamps the MC rate factor to [1/cap, cap], whose bounds "
                "cross below 1)", cfg->mc_rate_rescale_cap);
    /* Measured physical constants: a mass, a coupling or a temperature that is
     * zero or negative is not a sensitivity study but a typo (the integrands
     * divide by me and take sqrt(E^2 - me^2)). The two anomalous magnetic
     * moments are absent -- kappa_n is negative by nature -- and Neff_SM only
     * needs >= 0. Mirrors primat/config.py's _PARAM_RANGE. */
    CPR_REQUIRE(cfg->consts.alphaem > 0,
                "alphaem=%.6g is out of range: must be > 0", cfg->consts.alphaem);
    CPR_REQUIRE(cfg->consts.GF > 0,
                "GF=%.6g is out of range: must be > 0", cfg->consts.GF);
    CPR_REQUIRE(cfg->consts.mZ > 0,
                "mZ=%.6g is out of range: must be > 0", cfg->consts.mZ);
    CPR_REQUIRE(cfg->consts.me > 0,
                "me=%.6g is out of range: must be > 0", cfg->consts.me);
    CPR_REQUIRE(cfg->consts.mn > 0,
                "mn=%.6g is out of range: must be > 0", cfg->consts.mn);
    CPR_REQUIRE(cfg->consts.mp > 0,
                "mp=%.6g is out of range: must be > 0", cfg->consts.mp);
    CPR_REQUIRE(cfg->consts.T0CMB > 0,
                "T0CMB=%.6g is out of range: must be > 0", cfg->consts.T0CMB);
    CPR_REQUIRE(cfg->consts.gA > 0,
                "gA=%.6g is out of range: must be > 0", cfg->consts.gA);
    CPR_REQUIRE(cfg->consts.Vud > 0,
                "Vud=%.6g is out of range: must be > 0", cfg->consts.Vud);
    CPR_REQUIRE(cfg->consts.radproton > 0,
                "radproton=%.6g is out of range: must be > 0", cfg->consts.radproton);
    CPR_REQUIRE(cfg->consts.ma > 0,
                "ma=%.6g is out of range: must be > 0", cfg->consts.ma);
    CPR_REQUIRE(cfg->consts.He4Overma > 0,
                "He4Overma=%.6g is out of range: must be > 0", cfg->consts.He4Overma);
    CPR_REQUIRE(cfg->consts.HOverma > 0,
                "HOverma=%.6g is out of range: must be > 0", cfg->consts.HOverma);

    /* non-negative doubles (a 1-sigma width may legitimately be 0) */
    CPR_REQUIRE(cfg->std_tau_n >= 0,
                "std_tau_n=%.6g is out of range: must be >= 0", cfg->std_tau_n);
    CPR_REQUIRE(cfg->consts.Neff_SM >= 0,
                "Neff_SM=%.6g is out of range: must be >= 0", cfg->consts.Neff_SM);

    /* strictly positive integer counts */
    CPR_REQUIRE(cfg->n_electron_table >= 1,
                "n_electron_table=%d is out of range: must be a positive integer (>= 1)", cfg->n_electron_table);
    CPR_REQUIRE(cfg->sampling_temperature_per_decade >= 1,
                "sampling_temperature_per_decade=%d is out of range: must be a positive integer (>= 1)", cfg->sampling_temperature_per_decade);
    CPR_REQUIRE(cfg->sampling_nTOp_per_decade >= 1,
                "sampling_nTOp_per_decade=%d is out of range: must be a positive integer (>= 1)", cfg->sampling_nTOp_per_decade);
    CPR_REQUIRE(cfg->sampling_nTOp_thermal_per_decade >= 1,
                "sampling_nTOp_thermal_per_decade=%d is out of range: must be a positive integer (>= 1)", cfg->sampling_nTOp_thermal_per_decade);
    CPR_REQUIRE(cfg->vegas_n_eval >= 1,
                "vegas_n_eval=%d is out of range: must be a positive integer (>= 1)", cfg->vegas_n_eval);
    CPR_REQUIRE(cfg->vegas_n_itn >= 1,
                "vegas_n_itn=%d is out of range: must be a positive integer (>= 1)", cfg->vegas_n_itn);
    CPR_REQUIRE(cfg->output_n_points >= 1,
                "output_n_points=%d is out of range: must be a positive integer (>= 1)", cfg->output_n_points);
    /* >= 2, not >= 1: a single-point grid leaves cpr_network_fill_buffer
     * reading g[ii+1] past the end of a one-element array (and
     * cpr_find_segment returning n-2 = SIZE_MAX). Python's equivalent
     * degenerate grid divides by zero. Two points is the true minimum for
     * the linear interpolant; the not-a-knot >= 4 rule applies to the
     * SOURCE table, not to this master grid. */
    CPR_REQUIRE(cfg->rate_grid_npts >= 2,
                "rate_grid_npts=%d is out of range: must be >= 2 (a one-point grid has no interval to interpolate on)", cfg->rate_grid_npts);
    CPR_REQUIRE(cfg->decay_n_points >= 1,
                "decay_n_points=%d is out of range: must be a positive integer (>= 1)", cfg->decay_n_points);
#undef CPR_REQUIRE

    /* Cross-field consistency, mirroring _validate_ranges in
     * primat/config.py: constraints that involve two parameters at once, each
     * of which otherwise surfaces only as an opaque "cpr_ode_bdf: step size
     * underflowed below machine precision" from deep inside the MT era. */
    if (cfg->rate_grid_T9_min >= cfg->rate_grid_T9_max) {
        *errmsg = malloc(224);
        snprintf(*errmsg, 224,
                 "rate_grid_T9_min=%.6g must be < rate_grid_T9_max=%.6g: they "
                 "bound the log-spaced master T9 grid every nuclear rate table "
                 "is resampled onto, which must be increasing",
                 cfg->rate_grid_T9_min, cfg->rate_grid_T9_max);
        return 1;
    }
    if (cfg->T_end_MeV >= cfg->T_start_cosmo_MeV) {
        *errmsg = malloc(224);
        snprintf(*errmsg, 224,
                 "T_end_MeV=%.6g must be < T_start_cosmo_MeV=%.6g: the "
                 "background and the nuclear network are integrated from "
                 "T_start_cosmo_MeV down to T_end_MeV",
                 cfg->T_end_MeV, cfg->T_start_cosmo_MeV);
        return 1;
    }

    /* fEDE is the EDE fraction of the total energy density at its peak;
     * background.c has (1 - fEDE) in the denominator, so fEDE >= 1 diverges. */
    if (!(cfg->fEDE >= 0.0 && cfg->fEDE < 1.0)) {
        *errmsg = malloc(128);
        snprintf(*errmsg, 128,
                 "fEDE=%.6g is out of range: must satisfy 0 <= fEDE < 1", cfg->fEDE);
        return 1;
    }

    /* wnEDE > 1/3 whenever EDE is on. background.c's _setup_EDE locates the
     * EDE-fraction peak at u^(3(1+wnEDE)) = 4/(3*wnEDE - 1), which has no
     * solution for wnEDE <= 1/3: such a component dilutes no faster than
     * radiation, so its fraction never peaks during radiation domination and
     * fEDE (defined at that peak) is meaningless. Without this check C's pow()
     * quietly produces NaN -- 4.0/0.0 -> +inf at wnEDE = 1/3, or pow(negative,
     * fractional) -> NaN below -- and the whole background silently becomes
     * NaN, whereas Python raises. Mirrors _validate_fEDE in config.py. */
    if (cfg->fEDE != 0.0 && !(cfg->wnEDE > 1.0 / 3.0)) {
        *errmsg = malloc(256);
        snprintf(*errmsg, 256,
                 "wnEDE=%.6g is out of range: must satisfy wnEDE > 1/3 when "
                 "fEDE > 0 (the EDE peak u^(3(1+wnEDE)) = 4/(3*wnEDE - 1) has "
                 "no solution otherwise). For V ~ (1 - cos phi)^n use "
                 "wnEDE = (n-1)/(n+1) with n >= 3.", cfg->wnEDE);
        return 1;
    }

    if (cfg->external_scale_factor && !cfg->incomplete_decoupling) {
        *errmsg = strdup("external_scale_factor=True requires incomplete_decoupling=True");
        return 1;
    }

    /* Validate spectral-distortion flag combination (mirrors
     * PRIMATConfig.__init__'s equivalent block in config.py). */
    if (cfg->spectral_distortions) {
        if (cfg->analytic_distortions) {
            if (cfg->incomplete_decoupling) {
                *errmsg = strdup(
                    "spectral_distortions=True with analytic_distortions=True "
                    "requires instantaneous decoupling (incomplete_decoupling=False).");
                return 1;
            }
        } else {
            if (!cfg->incomplete_decoupling) {
                *errmsg = strdup(
                    "spectral_distortions=True with analytic_distortions=False "
                    "requires incomplete_decoupling=True (the full NEVO spectrum "
                    "file is only available in the non-instantaneous decoupling mode).");
                return 1;
            }
        }
    }

    /* NEVO override existence/shape checks: deferred to neutrino_history.c
     * (Phase 3a), which owns resolve_nevo_path() and the CSV column counts. */

    return 0;
}

void cpr_config_free(CPRConfig *cfg)
{
    for (size_t i = 0; i < FIELD_TABLE_N; i++) {
        if (FIELD_TABLE[i].kind == F_STRING) {
            void *field = (char *)cfg + FIELD_TABLE[i].offset;
            free(*(char **)field);
            *(char **)field = NULL;
        }
    }
    cpr_rxnmap_free(&cfg->p_rxn);
    cpr_rxnmap_free(&cfg->delta_rxn);
    free(cfg->nuclides.items);
    cfg->nuclides.items = NULL;
    cfg->nuclides.n = 0;
    /* Tabulated extra_rho arrays (config.h): owned here, malloc'd by the
     * caller (the _wrapper.c bridge) after cpr_config_init_defaults. */
    free(cfg->extra_rho_T);
    free(cfg->extra_rho_val);
    cfg->extra_rho_T = NULL;
    cfg->extra_rho_val = NULL;
    cfg->extra_rho_n = 0;
}
