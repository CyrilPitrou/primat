/* test_ini.c -- verifies cpr_ini_load applies KEY=VALUE overrides on top of
 * defaults, that examples/run_large_amax8.ini round-trips amax/network, and
 * that the three policies the loader owes the Python side hold:
 *
 *   1. every applied pair is COLLECTED into the caller's CPRParamList. MC
 *      workers rebuild their config from defaults + that list, so an ini key
 *      missing from it silently vanishes from every MC sample -- the sigmas
 *      would then describe a different model than the central value printed
 *      next to them.
 *   2. a wrong-typed value is FATAL (PRIMATConfig raises TypeError), and the
 *      target field is left untouched. Warning-and-continuing here used to
 *      leave a freed pointer in the config: the run went on to read it, and
 *      cpr_config_free double-freed it at exit.
 *   3. an unknown key is only a WARNING while strict_params is off
 *      (PRIMATConfig's documented default), and fatal once it is on.
 */
#include "config.h"
#include "ini.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

/* Writes `body` to a scratch ini file and returns its path (a static buffer). */
static const char *write_ini(const char *body)
{
    static char path[] = "build/test_ini_tmp.ini";
    FILE *f = fopen(path, "w");
    if (!f) { printf("FAIL: cannot create %s\n", path); exit(1); }
    fputs(body, f);
    fclose(f);
    return path;
}

int main(void)
{
    char *err = NULL;
    CPRConfig cfg;
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL cpr_config_init_defaults: %s\n", err);
        return 1;
    }

    /* ---- 1. shipped example ini, and collection into a param list ---- */
    CPRParamList pl;
    memset(&pl, 0, sizeof(pl));
    if (cpr_ini_load(&cfg, "examples/run_large_amax8.ini", &pl, &err)) {
        printf("FAIL cpr_ini_load: %s\n", err);
        return 1;
    }
    CHECK(strcmp(cfg.network, "large") == 0, "network == large after ini load");
    CHECK(cfg.amax == 8, "amax == 8 after ini load");
    CHECK(cfg.Omegabh2_ == 0.02242, "Omegabh2 == 0.02242 after ini load");

    /* Every key the ini set must be replayable to an MC worker. */
    int saw_network = 0, saw_amax = 0;
    for (size_t i = 0; i < pl.n; i++) {
        if (strcmp(pl.items[i].key, "network") == 0) {
            saw_network = 1;
            CHECK(pl.items[i].value.type == CPR_STRING
                  && strcmp(pl.items[i].value.v.s, "large") == 0,
                  "collected network value survives as its own copy");
        }
        if (strcmp(pl.items[i].key, "amax") == 0) saw_amax = 1;
    }
    CHECK(pl.n > 0, "ini keys are collected for MC replay");
    CHECK(saw_network, "collected list contains 'network'");
    CHECK(saw_amax, "collected list contains 'amax'");
    cpr_paramlist_free(&pl);

    if (cpr_config_validate(&cfg, &err)) {
        printf("FAIL cpr_config_validate: %s\n", err);
        failures++;
    } else {
        printf("ok: validate succeeds after ini load\n");
    }
    cpr_config_free(&cfg);

    /* ---- 2. a wrong-typed value is fatal, and the field is untouched ---- */
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL cpr_config_init_defaults (2): %s\n", err);
        return 1;
    }
    const char *network_before = cfg.network;
    err = NULL;
    int rc = cpr_ini_load(&cfg, write_ini("network = 3\n"), NULL, &err);
    CHECK(rc != 0, "a wrong-typed ini value is fatal");
    CHECK(err != NULL, "the type error carries a message");
    /* Same pointer, same contents: nothing was freed or replaced. Reading
     * cfg.network at all is only safe because of that. */
    CHECK(cfg.network == network_before && strcmp(cfg.network, "small") == 0,
          "the target field is untouched after a type error");
    free(err);
    cpr_config_free(&cfg);   /* must not double-free */

    /* ---- 3. unknown keys: warning by default, fatal under strict_params ---- */
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL cpr_config_init_defaults (3): %s\n", err);
        return 1;
    }
    err = NULL;
    printf("(one 'unknown parameter key' warning is expected next)\n");
    rc = cpr_ini_load(&cfg, write_ini("no_such_key = 1\nnetwork = large\n"),
                      NULL, &err);
    CHECK(rc == 0, "an unknown ini key is only a warning (strict_params=False)");
    CHECK(strcmp(cfg.network, "large") == 0,
          "the load continues past an unknown key");
    free(err);

    err = NULL;
    rc = cpr_ini_load(&cfg, write_ini("strict_params = True\nno_such_key = 1\n"),
                      NULL, &err);
    CHECK(rc != 0, "an unknown ini key is fatal under strict_params=True");
    free(err);

    /* ---- 4. an empty value is not a literal ---- */
    err = NULL;
    rc = cpr_ini_load(&cfg, write_ini("network =\n"), NULL, &err);
    CHECK(rc != 0, "an empty ini value is rejected, not read as \"\"");
    free(err);

    cpr_config_free(&cfg);
    remove("build/test_ini_tmp.ini");

    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
