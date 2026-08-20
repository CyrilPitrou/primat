/* test_fuzz_regressions.c -- replays every fuzzing input the tree keeps.
 *
 * Two kinds of input, both under tests/fuzz (see tests/fuzz/README.md):
 *
 *   regressions/<kind>/  a minimised reproducer for a defect the fuzzer found.
 *                        Each one is asserted to produce the specific
 *                        diagnosis it was fixed to produce, not merely to
 *                        avoid crashing.
 *   corpus/<kind>/       the evolved corpus the campaign ended with. Replayed
 *                        wholesale: under `make debug-test` (ASan+UBSan) this
 *                        is the cheapest guard against a loader regressing on
 *                        input shapes the fuzzer already reached.
 *
 * Run from primat-c/ (the Makefile's `test` target does).
 */
#include "config.h"
#include "ini.h"
#include "network_data.h"
#include "table_io.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

#define DATA_DIR "../primat/data"

/* ------------------------------------------------------------------ */
/* The reproducers, each pinned to the message its fix installed.      */
/* ------------------------------------------------------------------ */

/* A non-nullable string parameter set to None used to leave the field NULL,
 * and cpr_config_validate's strcmp(cfg->network, "small") then died on
 * SIGSEGV printing nothing -- while PRIMATConfig raised a TypeError naming
 * the key. Both doors into that field are checked: the ini loader (where the
 * fuzzer found it) and cpr_config_set_by_name (what --set uses). */
static void test_network_none_is_a_typed_error(void)
{
    const char *expected =
        "network=None has the wrong type: expected str, got NoneType";

    CPRConfig cfg;
    char *err = NULL;
    if (cpr_config_init_defaults(&cfg, DATA_DIR, &err)) {
        printf("FAIL: cannot init defaults: %s\n", err ? err : "?");
        failures++;
        free(err);
        cpr_config_free(&cfg);
        return;
    }

    CPRParam none = { CPR_NONE, {0} };
    char *set_err = NULL;
    int rc = cpr_config_set_by_name(&cfg, "network", none, &set_err);
    CHECK(rc != CPR_SET_OK, "network=None is rejected by cpr_config_set_by_name");
    CHECK(set_err && strcmp(set_err, expected) == 0,
          "the rejection message matches PRIMATConfig's TypeError word for word");
    CHECK(cfg.network != NULL, "the field is left holding its previous value");
    free(set_err);

    /* Same value through the ini door. */
    CPRParamList cp;
    memset(&cp, 0, sizeof(cp));
    err = NULL;
    rc = cpr_ini_load(&cfg, "tests/fuzz/regressions/ini/network_is_none.ini",
                      &cp, &err);
    CHECK(rc != 0, "an ini setting network = none is fatal");
    CHECK(err && strstr(err, expected) != NULL,
          "the ini error carries the same sentence");
    free(err);
    cpr_paramlist_free(&cp);

    /* The nullable path parameters are unaffected: None is their documented
     * "use the default / skip this output" value. */
    err = NULL;
    rc = cpr_config_set_by_name(&cfg, "output_file", none, &err);
    CHECK(rc == CPR_SET_OK && cfg.output_file == NULL,
          "output_file=None is still accepted and clears the field");
    free(err);

    cpr_config_free(&cfg);
}

/* A line longer than a reader's buffer used to come back from fgets in
 * chunks, each treated as a line of its own. For cpr_table_read that turned
 * the tail of an over-long "#" header into a data row -- three rows where
 * numpy.loadtxt reads two, silently, and only for the one input that
 * straddles the buffer. Comment lines of any length are now skipped whole
 * (matching numpy), and an over-long *data* line is a named error rather
 * than a partial parse. */
static void test_an_overlong_line_is_never_parsed_in_part(void)
{
    const size_t line_cap = 8191;          /* cpr_table_read's char line[8192] */
    const char *path = "build/test_fuzz_overlong.txt";

    /* Padded so the split lands exactly on the numbers at the end. */
    FILE *f = fopen(path, "w");
    if (!f) { printf("FAIL: cannot create %s\n", path); failures++; return; }
    fputs("# ref=", f);
    for (size_t i = 0; i < line_cap - 6; i++) fputc('A', f);
    fputs("7.0 8.0 9.0\n1.0 1.0 0.1\n2.0 2.0 0.2\n", f);
    fclose(f);

    CPRTable t;
    char *err = NULL;
    int rc = cpr_table_read(path, 0, &t, &err);
    CHECK(rc == 0, "a table with an over-long comment header still loads");
    if (rc == 0) {
        CHECK(t.n_rows == 2, "its comment contributes no data row");
        CHECK(t.n_rows == 2 && t.cols[0][0] == 1.0,
              "the first data row is the first real one");
        cpr_table_free(&t);
    } else {
        free(err);
    }

    /* The same length in a data line cannot be parsed at all. */
    f = fopen(path, "w");
    if (!f) { printf("FAIL: cannot reopen %s\n", path); failures++; return; }
    fputs("1.0 1.0 0.1\n", f);
    for (size_t i = 0; i < line_cap; i++) fputc(i % 4 == 3 ? ' ' : '1', f);
    fputs(" 2.0\n", f);
    fclose(f);

    err = NULL;
    rc = cpr_table_read(path, 0, &t, &err);
    CHECK(rc != 0, "an over-long data line is rejected");
    CHECK(err && strstr(err, "longer than") != NULL,
          "the rejection names the length limit");
    free(err);
    remove(path);
}

/* ------------------------------------------------------------------ */
/* Corpus replay                                                       */
/* ------------------------------------------------------------------ */

typedef void (*ReplayFn)(const char *path);

static void replay_table(const char *path)
{
    CPRTable t;
    char *err = NULL;
    if (cpr_table_read(path, 0, &t, &err) == 0) {
        if (t.n_cols >= 2) {
            char *verr = NULL;
            if (cpr_validate_rate_table(t.cols[0], t.cols[1],
                                        t.n_cols >= 3 ? t.cols[2] : NULL,
                                        t.n_rows, "replay", &verr) != 0)
                free(verr);
        }
        cpr_table_free(&t);
    } else {
        free(err);
    }
}

static void replay_network_list(const char *path)
{
    CPRNetworkList l;
    char *err = NULL;
    if (cpr_load_network_list(path, &l, &err) == 0) cpr_network_list_free(&l);
    else free(err);
}

static void replay_decays(const char *path)
{
    CPRDecayTable t;
    char *err = NULL;
    if (cpr_load_decays(path, &t, &err) == 0) cpr_decay_table_free(&t);
    else free(err);
}

static void replay_detailed_balance(const char *path)
{
    CPRDetailedBalanceTable t;
    char *err = NULL;
    if (cpr_load_detailed_balance(path, &t, &err) == 0) cpr_detailed_balance_free(&t);
    else free(err);
}

static void replay_reactions_large(const char *path)
{
    CPRReactionTable t;
    char *err = NULL;
    if (cpr_load_reactions_large(path, &t, &err) == 0) cpr_reaction_table_free(&t);
    else free(err);
}

static void replay_ini(const char *path)
{
    CPRConfig cfg;
    char *err = NULL;
    if (cpr_config_init_defaults(&cfg, DATA_DIR, &err) != 0) {
        free(err);
        cpr_config_free(&cfg);
        return;
    }
    CPRParamList cp;
    memset(&cp, 0, sizeof(cp));
    if (cpr_ini_load(&cfg, path, &cp, &err) != 0) {
        free(err);
    } else {
        char *verr = NULL;
        if (cpr_config_validate(&cfg, &verr) != 0) free(verr);
    }
    cpr_paramlist_free(&cp);
    cpr_config_free(&cfg);
}

static const struct { const char *dir; ReplayFn fn; } REPLAY[] = {
    { "rate_table",       replay_table },
    { "nevo",             replay_table },
    { "nevo_spectral",    replay_table },
    { "cache",            replay_table },
    { "network_list",     replay_network_list },
    { "decays",           replay_decays },
    { "detailed_balance", replay_detailed_balance },
    { "reactions_large",  replay_reactions_large },
    { "ini",              replay_ini },
};

/* Replays every file under tests/fuzz/<root>/<dir>/, returning the count. */
static int replay_dir(const char *root, const char *dir, ReplayFn fn)
{
    char path[512];
    snprintf(path, sizeof(path), "tests/fuzz/%s/%s", root, dir);
    DIR *dp = opendir(path);
    if (!dp) return 0;
    struct dirent *de;
    int n = 0;
    while ((de = readdir(dp))) {
        if (de->d_name[0] == '.') continue;
        char file[1024];
        snprintf(file, sizeof(file), "%s/%s", path, de->d_name);
        fn(file);
        n++;
    }
    closedir(dp);
    return n;
}

int main(void)
{
    printf("== test_fuzz_regressions ==\n");
    test_network_none_is_a_typed_error();
    test_an_overlong_line_is_never_parsed_in_part();

    /* The loaders diagnose each malformed input on stderr, which for a few
     * hundred deliberately-malformed files buries the test's own result.
     * Silenced for the replay only, then restored. */
    fflush(stderr);
    int saved_stderr = dup(2);
    FILE *devnull = freopen("/dev/null", "w", stderr);
    int total = 0;
    for (size_t i = 0; i < sizeof(REPLAY) / sizeof(REPLAY[0]); i++) {
        total += replay_dir("corpus", REPLAY[i].dir, REPLAY[i].fn);
        total += replay_dir("regressions", REPLAY[i].dir, REPLAY[i].fn);
    }
    fflush(stderr);
    if (devnull) dup2(saved_stderr, 2);
    close(saved_stderr);
    /* The corpora are committed, so an empty walk means the tree lost them
     * rather than that everything passed. */
    CHECK(total >= 100, "the committed fuzz corpora were found and replayed");
    printf("replayed %d fuzzing inputs\n", total);

    printf(failures ? "FAILURES: %d\n" : "all ok (%d failures)\n", failures);
    return failures ? 1 : 0;
}
