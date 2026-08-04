/* test_config.c -- checks cpr_config_init_defaults/cpr_config_set_by_name's
 * generic override dispatch, the C equivalent of
 * ../../tests/test_config.py's PyPRConfig checks: sane defaults, that a
 * named override actually lands on the right typed field (bool/int/double/
 * string), that an unknown key is reported as an error (mirrors Python's
 * "unknown key" warning -- the C port makes it a hard error via *errmsg
 * rather than a warning, since cli.c/ini.c are the ones that decide to
 * downgrade it), and that p_<rxn>/delta_<rxn> keys route into the
 * corresponding CPRRxnMap instead of the fixed-field table.
 */
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

int main(void)
{
    char *err = NULL;
    CPRConfig cfg;
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL cpr_config_init_defaults: %s\n", err);
        return 1;
    }

    /* ---- Defaults sanity (mirrors test_config.py's test_default_construction) ---- */
    CHECK(strcmp(cfg.network, "small") == 0, "default network is 'small'");
    CHECK(cfg.numerical_precision == 1e-7, "default numerical_precision is 1e-7");
    CHECK(cfg.incomplete_decoupling == 1, "default incomplete_decoupling is True");
    CHECK(cfg.QED_corrections == 1, "default QED_corrections is True");

    /* ---- Overrides land on the right typed field, one per kind ---- */
    char *set_err = NULL;
    CPRParam p_bool = { .type = CPR_BOOL, .v.b = 0 };
    CHECK(cpr_config_set_by_name(&cfg, "QED_corrections", p_bool, &set_err) == 0,
          "set_by_name accepts a bool override");
    CHECK(cfg.QED_corrections == 0, "QED_corrections overridden to False");

    CPRParam p_int = { .type = CPR_INT, .v.i = 1234 };
    CHECK(cpr_config_set_by_name(&cfg, "n_electron_table", p_int, &set_err) == 0,
          "set_by_name accepts an int override");
    CHECK(cfg.n_electron_table == 1234, "n_electron_table overridden to 1234");

    CPRParam p_double = { .type = CPR_DOUBLE, .v.d = 1e-9 };
    CHECK(cpr_config_set_by_name(&cfg, "numerical_precision", p_double, &set_err) == 0,
          "set_by_name accepts a double override");
    CHECK(cfg.numerical_precision == 1e-9, "numerical_precision overridden to 1e-9");

    CPRParam p_str = { .type = CPR_STRING, .v.s = "large" };
    CHECK(cpr_config_set_by_name(&cfg, "network", p_str, &set_err) == 0,
          "set_by_name accepts a string override");
    CHECK(strcmp(cfg.network, "large") == 0, "network overridden to 'large'");

    /* An int is a valid CPR_DOUBLE-field value too (numeric widening,
     * mirroring Python's duck-typed int/float interchangeability). */
    CPRParam p_int_for_double = { .type = CPR_INT, .v.i = 2 };
    CHECK(cpr_config_set_by_name(&cfg, "numerical_precision", p_int_for_double, &set_err) == 0,
          "set_by_name widens an int to a double field");
    CHECK(cfg.numerical_precision == 2.0, "numerical_precision widened to 2.0");

    /* A string value on a bool field is a type mismatch -> error, not a
     * silent truthy coercion. */
    CPRParam p_str_for_bool = { .type = CPR_STRING, .v.s = "true" };
    free(set_err); set_err = NULL;
    int rc = cpr_config_set_by_name(&cfg, "QED_corrections", p_str_for_bool, &set_err);
    CHECK(rc != 0, "set_by_name rejects a string value for a bool field");
    CHECK(set_err != NULL, "type-mismatch error message is set");
    free(set_err); set_err = NULL;

    /* ---- Unknown key is an error (mirrors test_unknown_key_warns) ---- */
    CPRParam p_anything = { .type = CPR_INT, .v.i = 1 };
    rc = cpr_config_set_by_name(&cfg, "not_a_real_parameter", p_anything, &set_err);
    CHECK(rc != 0, "unknown key is reported as an error");
    CHECK(set_err != NULL && strstr(set_err, "not_a_real_parameter") != NULL,
          "unknown-key error message names the bad key");
    free(set_err); set_err = NULL;

    /* ---- p_<rxn>/delta_<rxn> route into the reaction-rate maps, not the
     * fixed-field table (mirrors test_p_rxn_typo_warns's sibling,
     * test_p_rxn_valid_reaction_does_not_warn) ---- */
    CPRParam p_rate = { .type = CPR_DOUBLE, .v.d = 0.5 };
    CHECK(cpr_config_set_by_name(&cfg, "p_n_p__d_g", p_rate, &set_err) == 0,
          "p_<rxn> key is accepted");
    CHECK(cpr_rxnmap_get(&cfg.p_rxn, "n_p__d_g") == 0.5,
          "p_n_p__d_g landed in cfg.p_rxn under the bare reaction name");
    CHECK(cpr_rxnmap_get(&cfg.p_rxn, "some_other_reaction") == 0.0,
          "cpr_rxnmap_get defaults to 0.0 for an unset reaction");

    CPRParam delta_rate = { .type = CPR_DOUBLE, .v.d = -0.25 };
    CHECK(cpr_config_set_by_name(&cfg, "delta_n_p__d_g", delta_rate, &set_err) == 0,
          "delta_<rxn> key is accepted");
    CHECK(cpr_rxnmap_get(&cfg.delta_rxn, "n_p__d_g") == -0.25,
          "delta_n_p__d_g landed in cfg.delta_rxn, independent of cfg.p_rxn");

    /* ---- Omegabh2 setter recomputes eta0b (dedicated branch in
     * cpr_config_set_by_name, not the generic field table) ---- */
    double eta0b_before = cfg.eta0b;
    CPRParam p_ombh2 = { .type = CPR_DOUBLE, .v.d = 0.02 };
    CHECK(cpr_config_set_by_name(&cfg, "Omegabh2", p_ombh2, &set_err) == 0,
          "Omegabh2 override is accepted");
    CHECK(cpr_config_get_Omegabh2(&cfg) == 0.02, "Omegabh2 getter reflects the override");
    CHECK(cfg.eta0b != eta0b_before, "eta0b was recomputed after the Omegabh2 override");

    /* ---- Physical/numerical range checks in cpr_config_validate,
     * mirroring test_config.py's range-check tests. Each bad value must make
     * cpr_config_validate return non-zero with a message naming the field.
     * strict_params default (0) and the strict_params field round-trip are
     * checked too. ---- */
    CHECK(cfg.strict_params == 0, "default strict_params is False (0)");
    CPRParam p_strict = { .type = CPR_BOOL, .v.b = 1 };
    CHECK(cpr_config_set_by_name(&cfg, "strict_params", p_strict, &set_err) == 0,
          "strict_params round-trips through the field table");
    CHECK(cfg.strict_params == 1, "strict_params overridden to True (1)");

    /* Reset to a clean valid config, then perturb one field at a time. */
    cpr_config_free(&cfg);
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL cpr_config_init_defaults (range block): %s\n", err);
        return 1;
    }
    char *verr = NULL;
    CHECK(cpr_config_validate(&cfg, &verr) == 0, "pristine defaults pass validation");
    free(verr); verr = NULL;

    /* Negative Omegabh2 -> out of range (routed through the setter). */
    cpr_config_set_Omegabh2(&cfg, -0.02);
    CHECK(cpr_config_validate(&cfg, &verr) != 0 && verr && strstr(verr, "Omegabh2"),
          "negative Omegabh2 is rejected with a naming message");
    free(verr); verr = NULL;
    cpr_config_set_Omegabh2(&cfg, 0.02242); /* restore */

    /* Non-positive tau_n. */
    cfg.tau_n = 0.0;
    CHECK(cpr_config_validate(&cfg, &verr) != 0 && verr && strstr(verr, "tau_n"),
          "tau_n=0 is rejected");
    free(verr); verr = NULL;
    cfg.tau_n = 878.4;

    /* Non-positive integer count. */
    cfg.rate_grid_npts = 0;
    CHECK(cpr_config_validate(&cfg, &verr) != 0 && verr && strstr(verr, "rate_grid_npts"),
          "rate_grid_npts=0 is rejected");
    free(verr); verr = NULL;
    cfg.rate_grid_npts = 1000;

    /* std_tau_n is allowed to be exactly 0 (non-negative, not strictly positive). */
    cfg.std_tau_n = 0.0;
    CHECK(cpr_config_validate(&cfg, &verr) == 0, "std_tau_n=0 is accepted (non-negative)");
    free(verr); verr = NULL;

    /* ---- Per-flavour neutrino chemical potentials. Defaults are the NAN
     * "inherit munuOverTnu" sentinel, resolved by cpr_config_xi_nu_e/mu/tau;
     * None round-trips back to NAN, a number pins that one flavour. Mirrors
     * PRIMATConfig.xi_nu_e / munuOverTnu_e=None. ---- */
    CHECK(isnan(cfg.munuOverTnu_e) && isnan(cfg.munuOverTnu_mu)
              && isnan(cfg.munuOverTnu_tau),
          "per-flavour ξ default to the NAN inherit-sentinel");
    cfg.munuOverTnu = 0.07;
    CHECK(cpr_config_xi_nu_e(&cfg) == 0.07 && cpr_config_xi_nu_mu(&cfg) == 0.07
              && cpr_config_xi_nu_tau(&cfg) == 0.07,
          "unset per-flavour ξ inherit the common munuOverTnu");
    CPRParam p_xe = { .type = CPR_DOUBLE, .v.d = 0.2 };
    CHECK(cpr_config_set_by_name(&cfg, "munuOverTnu_e", p_xe, &set_err) == 0,
          "munuOverTnu_e override is accepted");
    CHECK(cpr_config_xi_nu_e(&cfg) == 0.2, "munuOverTnu_e override wins for ξ_e");
    CHECK(cpr_config_xi_nu_mu(&cfg) == 0.07 && cpr_config_xi_nu_tau(&cfg) == 0.07,
          "the other two flavours still inherit munuOverTnu");
    CPRParam p_none = { .type = CPR_NONE };
    CHECK(cpr_config_set_by_name(&cfg, "munuOverTnu_e", p_none, &set_err) == 0
              && isnan(cfg.munuOverTnu_e),
          "None resets munuOverTnu_e back to the NAN inherit-sentinel");
    CHECK(cpr_config_xi_nu_e(&cfg) == 0.07,
          "after reset ξ_e inherits munuOverTnu again");

    /* ---- cache_dir redirect + cache_plasma_weak/ overlay. Mirrors
     * tests/test_cache_utils.py's cache_dir tests: unset -> the write dir is
     * <data_dir>/cache_plasma_weak/<sub>; set -> <cache_dir>/<sub>; and the
     * READ resolver still finds a shipped file (present only in the package
     * tree) even when cache_dir points elsewhere (overlay, never shadowed). ---- */
    {
        char buf[CPR_PATH_BUF_LEN2];

        /* Unset (default): write dir under <data_dir>/cache_plasma_weak/. */
        cpr_config_cache_write_dir(&cfg, "weak", buf, sizeof(buf));
        CHECK(strcmp(buf, "../primat/data/cache_plasma_weak/weak") == 0,
              "cache_write_dir(weak) defaults to <data_dir>/cache_plasma_weak/weak");
        cpr_config_cache_write_dir(&cfg, "plasma", buf, sizeof(buf));
        CHECK(strcmp(buf, "../primat/data/cache_plasma_weak/plasma") == 0,
              "cache_write_dir(plasma) defaults to <data_dir>/cache_plasma_weak/plasma");

        /* A file present only in the shipped tree resolves there even before
         * cache_dir is set (overlay base = the package copy, always last). */
        cpr_config_resolve_cache_file(&cfg, "weak",
            "nTOp_b8cdcc18d4677cc5.txt", buf, sizeof(buf));
        CHECK(strcmp(buf, "../primat/data/cache_plasma_weak/weak/"
                          "nTOp_b8cdcc18d4677cc5.txt") == 0,
              "resolve_cache_file finds the shipped weak cache (default)");

        /* Set cache_dir: writes redirect, but the shipped file (absent under
         * cache_dir) still resolves to the package copy (overlay fallback). */
        CPRParam p_cdir = { .type = CPR_STRING, .v.s = "/tmp/primat-cache-xyz" };
        CHECK(cpr_config_set_by_name(&cfg, "cache_dir", p_cdir, &set_err) == 0,
              "cache_dir override is accepted (generic string field)");
        free(set_err); set_err = NULL;
        cpr_config_cache_write_dir(&cfg, "weak", buf, sizeof(buf));
        CHECK(strcmp(buf, "/tmp/primat-cache-xyz/weak") == 0,
              "cache_write_dir(weak) redirects under cache_dir");
        cpr_config_resolve_cache_file(&cfg, "weak",
            "nTOp_b8cdcc18d4677cc5.txt", buf, sizeof(buf));
        CHECK(strcmp(buf, "../primat/data/cache_plasma_weak/weak/"
                          "nTOp_b8cdcc18d4677cc5.txt") == 0,
              "resolve_cache_file falls back to the shipped weak cache when "
              "cache_dir lacks it (overlay, never shadowed)");
    }

    /* ---- Absurdly long data_dir/user_nuclear_dir/cache_dir: every
     * path-building helper must truncate safely (NUL-terminated, never
     * write past the caller's buffer) rather than overflow, and warn on
     * stderr instead of silently handing back a mangled path (CODE_REVIEW.md
     * item 6: "verify each caller checks the return or that truncation is
     * detected once centrally"). data_dir is a fixed CPR_DATA_DIR_LEN
     * buffer (the tightest case); user_nuclear_dir/cache_dir are malloc'd
     * char* with no length cap at cpr_config_set_by_name time, so they are
     * the more realistic way a user could supply an oversized path. ---- */
    {
        char huge[CPR_DATA_DIR_LEN + 1000];
        memset(huge, 'a', sizeof(huge) - 1);
        huge[sizeof(huge) - 1] = '\0';

        CPRConfig cfg2;
        char *err2 = NULL;
        /* huge is not a real directory, so init_defaults legitimately fails
         * (nuclides.csv not found under it) -- the point of this block is
         * only that it fails *cleanly* (no crash, no overflow) rather than
         * succeeding; data_dir is set/truncated before that later nuclides.csv
         * check runs, so it is still safe to inspect and to feed to the
         * path-building helpers below even though init "failed". */
        CHECK(cpr_config_init_defaults(&cfg2, huge, &err2) != 0,
              "init_defaults reports failure for a non-existent (truncated) data_dir");
        CHECK(strlen(cfg2.data_dir) == CPR_DATA_DIR_LEN - 1,
              "an oversized data_dir is truncated to exactly the buffer capacity");
        CHECK(cfg2.data_dir[CPR_DATA_DIR_LEN - 1] == '\0',
              "the truncated data_dir is NUL-terminated");

        char buf2[CPR_PATH_BUF_LEN2];
        cpr_config_resolve_rates_path(&cfg2, "nuclear/networks/small.txt", buf2, sizeof(buf2));
        CHECK(strlen(buf2) < sizeof(buf2),
              "resolve_rates_path on a truncated data_dir stays within its output buffer");
        CHECK(buf2[strlen(buf2)] == '\0',
              "resolve_rates_path output is NUL-terminated");

        /* user_nuclear_dir: malloc'd char*, no length cap -- exercise a
         * candidate buffer smaller than the field itself. */
        char *set_err2 = NULL;
        CPRParam p_unc = { .type = CPR_STRING, .v.s = huge };
        CHECK(cpr_config_set_by_name(&cfg2, "user_nuclear_dir", p_unc, &set_err2) == 0,
              "user_nuclear_dir accepts an absurdly long override");
        free(set_err2); set_err2 = NULL;
        cpr_config_resolve_rates_path(&cfg2, "nuclear/networks/small.txt", buf2, sizeof(buf2));
        CHECK(strlen(buf2) < sizeof(buf2),
              "resolve_rates_path with an absurdly long user_nuclear_dir stays within bounds");

        /* cache_dir: same malloc'd-string story, exercised through both
         * cache_write_dir and resolve_cache_file (the latter also builds an
         * "out + n" tail-append that must not underflow outsize-n). */
        CPRParam p_cd = { .type = CPR_STRING, .v.s = huge };
        CHECK(cpr_config_set_by_name(&cfg2, "cache_dir", p_cd, &set_err2) == 0,
              "cache_dir accepts an absurdly long override");
        free(set_err2); set_err2 = NULL;
        cpr_config_cache_write_dir(&cfg2, "weak", buf2, sizeof(buf2));
        CHECK(strlen(buf2) < sizeof(buf2),
              "cache_write_dir with an absurdly long cache_dir stays within bounds");
        cpr_config_resolve_cache_file(&cfg2, "weak", "nTOp_deadbeef.txt", buf2, sizeof(buf2));
        CHECK(strlen(buf2) < sizeof(buf2),
              "resolve_cache_file with an absurdly long cache_dir stays within bounds");

        free(err2);
        cpr_config_free(&cfg2);
    }

    cpr_config_free(&cfg);

    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
