/* fuzz_engine.c -- a small coverage-guided fuzzing engine.
 *
 * Apple's Command Line Tools clang instruments with -fsanitize-coverage but
 * ships no libFuzzer runtime (libclang_rt.fuzzer_osx.a is absent), so this
 * supplies the same loop: an edge map fed by __sanitizer_cov_trace_pc_guard,
 * a corpus that grows whenever an input reaches a new edge, and a mutator.
 * Targets keep libFuzzer's LLVMFuzzerTestOneInput signature, so they build
 * against the real runtime unchanged where it exists.
 *
 * Crashes are not caught in-process: the sanitizers end the process, so the
 * engine writes each input to <artifacts>/.current before running it and
 * run_fuzz.py turns that file into a saved artifact after an abnormal exit.
 */
#include "fuzz.h"

#include <dirent.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/* ------------------------------------------------------------------ */
/* Edge coverage                                                       */
/* ------------------------------------------------------------------ */

#define MAP_SIZE (1u << 18)

static uint8_t cur_map[MAP_SIZE];    /* edges hit by the input being run */
static uint8_t tot_map[MAP_SIZE];    /* edges hit by anything so far */
static uint32_t n_guards = 0;

void __sanitizer_cov_trace_pc_guard_init(uint32_t *start, uint32_t *stop)
{
    static uint32_t next = 1;
    if (start == stop || *start) return;
    for (uint32_t *p = start; p < stop; p++) *p = next++;
    n_guards = next - 1;
}

void __sanitizer_cov_trace_pc_guard(uint32_t *guard)
{
    if (!*guard) return;
    cur_map[*guard & (MAP_SIZE - 1)] = 1;
}

static size_t map_popcount(const uint8_t *m)
{
    size_t n = 0;
    for (size_t i = 0; i < MAP_SIZE; i++) n += m[i] != 0;
    return n;
}

/* ------------------------------------------------------------------ */
/* Scratch directory                                                   */
/* ------------------------------------------------------------------ */

static char g_tmpdir[512];

const char *fuzz_tmpdir(void) { return g_tmpdir; }

const char *fuzz_data_dir(void)
{
    const char *d = getenv("CPRIMAT_FUZZ_DATA_DIR");
    return d ? d : "../primat/data";
}

static void mkdir_p(const char *path)
{
    char buf[1024];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char *p = buf + 1; *p; p++) {
        if (*p == '/') { *p = '\0'; mkdir(buf, 0755); *p = '/'; }
    }
    mkdir(buf, 0755);
}

const char *fuzz_write_file(const char *relpath, const uint8_t *data, size_t size)
{
    static char full[1024];
    snprintf(full, sizeof(full), "%s/%s", g_tmpdir, relpath);
    char dir[1024];
    snprintf(dir, sizeof(dir), "%s", full);
    char *slash = strrchr(dir, '/');
    if (slash) { *slash = '\0'; mkdir_p(dir); }
    FILE *f = fopen(full, "wb");
    if (!f) { perror("fuzz_write_file"); exit(2); }
    if (size) fwrite(data, 1, size, f);
    fclose(f);
    return full;
}

/* ------------------------------------------------------------------ */
/* Corpus                                                              */
/* ------------------------------------------------------------------ */

#define FUZZ_MAX_INPUT (1u << 16)
#define MAX_CORPUS 4096

typedef struct { uint8_t *buf; size_t len; } Input;

static Input corpus[MAX_CORPUS];
static size_t n_corpus = 0;

static void corpus_add(const uint8_t *d, size_t n, const char *dir)
{
    if (n_corpus >= MAX_CORPUS) return;
    Input *in = &corpus[n_corpus++];
    in->buf = malloc(n ? n : 1);
    memcpy(in->buf, d, n);
    in->len = n;
    if (!dir) return;
    /* Persist so a restart after a crash keeps what was learned. */
    char path[1024];
    snprintf(path, sizeof(path), "%s/id_%06zu", dir, n_corpus - 1);
    FILE *f = fopen(path, "wb");
    if (f) { if (n) fwrite(d, 1, n, f); fclose(f); }
}

/* ------------------------------------------------------------------ */
/* Mutation                                                            */
/* ------------------------------------------------------------------ */

static uint64_t rng_state = 0x2545F4914F6CDD1DULL;

static uint64_t rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static size_t rnd_below(size_t n) { return n ? (size_t)(rnd() % n) : 0; }

/* Tokens that matter to line-oriented numeric parsers: the values that break
 * strtod/strtol, the separators these formats key on, and the shapes past
 * passes found by hand (empty, comments-only, duplicate, truncated). */
static const char *const dict[] = {
    "nan", "NaN", "-nan", "inf", "-inf", "INF", "1e999", "-1e999", "1e-999",
    "0x10", "0b1", "1_000", "9223372036854775808", "-9223372036854775809",
    "1.7976931348623157e308", "4.9406564584124654e-324", "0", "-0", "1", "-1",
    "#", "# ", ",", ",,,,", "=", "==", ";", "\n", "\r\n", "\t", "  ",
    "true", "false", "none", "None", "True", "False",
    "name,N,Z,A,Q,mass_excess_keV,spin", "reaction,Q_keV,alpha,beta,gamma",
    "n_p__d_g", "n_p__d_g,n_p__d_g_primat.txt", "n", "p", "H2", "He4", "Li7",
    "network", "small", "large", "amax", "verbose", "data_dir", "cache_dir",
    "# fingerprint_hash: 0123456789abcdef", "# fingerprint: {}",
    "ref=", "%", "\\", "\"", "'", "[", "]", "{", "}",
};
#define N_DICT (sizeof(dict) / sizeof(dict[0]))

static size_t mutate(uint8_t *out, const uint8_t *in, size_t len)
{
    if (len > FUZZ_MAX_INPUT) len = FUZZ_MAX_INPUT;
    memcpy(out, in, len);
    int rounds = 1 + (int)rnd_below(4);
    for (int r = 0; r < rounds; r++) {
        switch (rnd_below(11)) {
        case 0:                                     /* bit flip */
            if (len) out[rnd_below(len)] ^= (uint8_t)(1u << rnd_below(8));
            break;
        case 1:                                     /* random byte */
            if (len) out[rnd_below(len)] = (uint8_t)rnd();
            break;
        case 2: {                                   /* printable byte */
            if (len) {
                static const char pr[] = "0123456789.eE+-, \t\n#,;=abcdefABCDEF_/";
                out[rnd_below(len)] = (uint8_t)pr[rnd_below(sizeof(pr) - 1)];
            }
            break;
        }
        case 3: {                                   /* delete a chunk */
            if (len > 1) {
                size_t pos = rnd_below(len), n = 1 + rnd_below(len - pos);
                memmove(out + pos, out + pos + n, len - pos - n);
                len -= n;
            }
            break;
        }
        case 4: {                                   /* duplicate a chunk */
            if (len && len < FUZZ_MAX_INPUT / 2) {
                size_t pos = rnd_below(len), n = 1 + rnd_below(len - pos);
                if (len + n > FUZZ_MAX_INPUT) n = FUZZ_MAX_INPUT - len;
                memmove(out + pos + n, out + pos, len - pos);
                memcpy(out + pos, out + pos + n, n);
                len += n;
            }
            break;
        }
        case 5: case 6: {                           /* insert a dict token */
            const char *tok = dict[rnd_below(N_DICT)];
            size_t n = strlen(tok);
            size_t pos = rnd_below(len + 1);
            if (len + n > FUZZ_MAX_INPUT) break;
            memmove(out + pos + n, out + pos, len - pos);
            memcpy(out + pos, tok, n);
            len += n;
            break;
        }
        case 7: {                                   /* overwrite with token */
            const char *tok = dict[rnd_below(N_DICT)];
            size_t n = strlen(tok);
            if (len >= n) memcpy(out + rnd_below(len - n + 1), tok, n);
            break;
        }
        case 8: {                                   /* splice from corpus */
            if (n_corpus) {
                Input *o = &corpus[rnd_below(n_corpus)];
                size_t keep = rnd_below(len + 1);
                size_t take = o->len ? rnd_below(o->len) : 0;
                if (keep + take > FUZZ_MAX_INPUT) take = FUZZ_MAX_INPUT - keep;
                memcpy(out + keep, o->buf, take);
                len = keep + take;
            }
            break;
        }
        case 9: {                                   /* repeat a line */
            const uint8_t *nl = memchr(out, '\n', len);
            if (nl && len < FUZZ_MAX_INPUT / 2) {
                size_t lineln = (size_t)(nl - out) + 1;
                size_t reps = 1 + rnd_below(8);
                for (size_t k = 0; k < reps && len + lineln <= FUZZ_MAX_INPUT; k++) {
                    memmove(out + lineln, out, len);
                    memcpy(out, out + lineln, lineln);
                    len += lineln;
                }
            }
            break;
        }
        case 10:                                    /* truncate */
            if (len) len = rnd_below(len);
            break;
        }
    }
    return len;
}

/* ------------------------------------------------------------------ */
/* Timeout                                                             */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t g_timeout_s = 10;

/* The engine's own report channel, duplicated before any target runs: the CLI
 * target reopens stdout and stderr on /dev/null (its console output is not
 * what is under test), which would otherwise swallow the run summary. */
static int g_report_fd = 2;

static void on_alarm(int sig)
{
    (void)sig;
    static const char msg[] = "\n=== TIMEOUT: input exceeded the per-exec budget ===\n";
    ssize_t w = write(g_report_fd, msg, sizeof(msg) - 1);
    (void)w;
    _exit(88);
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

static int run_one(const uint8_t *d, size_t n, const char *cur_path)
{
    if (cur_path) {
        FILE *f = fopen(cur_path, "wb");
        if (f) { if (n) fwrite(d, 1, n, f); fclose(f); }
    }
    memset(cur_map, 0, sizeof(cur_map));
    alarm((unsigned)g_timeout_s);
    LLVMFuzzerTestOneInput(d, n);
    alarm(0);
    int novel = 0;
    for (size_t i = 0; i < MAP_SIZE; i++)
        if (cur_map[i] && !tot_map[i]) { tot_map[i] = 1; novel = 1; }
    return novel;
}

int main(int argc, char **argv)
{
    const char *corpus_dir = NULL, *art_dir = NULL;
    long runs = 10000;
    unsigned long seed = 0;
    const char *single = NULL;

    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "-runs=", 6) == 0) runs = atol(argv[i] + 6);
        else if (strncmp(argv[i], "-seed=", 6) == 0) seed = strtoul(argv[i] + 6, NULL, 10);
        else if (strncmp(argv[i], "-timeout=", 9) == 0) g_timeout_s = atoi(argv[i] + 9);
        else if (strncmp(argv[i], "-artifacts=", 11) == 0) art_dir = argv[i] + 11;
        else if (strncmp(argv[i], "-run_one=", 9) == 0) single = argv[i] + 9;
        else corpus_dir = argv[i];
    }
    if (seed) rng_state = seed * 6364136223846793005ULL + 1442695040888963407ULL;

    snprintf(g_tmpdir, sizeof(g_tmpdir), "/tmp/cprfuzz.%d", (int)getpid());
    mkdir_p(g_tmpdir);
    signal(SIGALRM, on_alarm);
    g_report_fd = dup(1);
    LLVMFuzzerInitialize(&argc, &argv);

    /* -run_one: replay a single saved artifact (used by the C regression
     * test and when minimising). */
    if (single) {
        FILE *f = fopen(single, "rb");
        if (!f) { fprintf(stderr, "cannot open %s\n", single); return 2; }
        static uint8_t buf[FUZZ_MAX_INPUT];
        size_t n = fread(buf, 1, sizeof(buf), f);
        fclose(f);
        run_one(buf, n, NULL);
        dprintf(g_report_fd, "ran %s (%zu bytes) ok\n", single, n);
        return 0;
    }

    char cur_path[1024] = {0};
    if (art_dir) { mkdir_p(art_dir); snprintf(cur_path, sizeof(cur_path), "%s/.current", art_dir); }

    /* Replay the on-disk corpus first: it rebuilds the edge map after a
     * restart, so a crash does not cost the run its accumulated coverage. */
    size_t n_seeds = 0;
    if (corpus_dir) {
        mkdir_p(corpus_dir);
        DIR *dp = opendir(corpus_dir);
        struct dirent *de;
        char names[MAX_CORPUS][256];
        size_t nn = 0;
        while (dp && (de = readdir(dp)) && nn < MAX_CORPUS) {
            if (de->d_name[0] == '.') continue;
            snprintf(names[nn++], 256, "%s", de->d_name);
        }
        if (dp) closedir(dp);
        for (size_t i = 0; i < nn; i++) {
            char p[1024];
            snprintf(p, sizeof(p), "%s/%s", corpus_dir, names[i]);
            FILE *f = fopen(p, "rb");
            if (!f) continue;
            static uint8_t buf[FUZZ_MAX_INPUT];
            size_t n = fread(buf, 1, sizeof(buf), f);
            fclose(f);
            run_one(buf, n, cur_path[0] ? cur_path : NULL);
            corpus_add(buf, n, NULL);
            n_seeds++;
        }
    }
    if (n_corpus == 0) corpus_add((const uint8_t *)"1 1\n", 4, NULL);

    dprintf(g_report_fd, "[fuzz] %s: %zu seeds, %u guards, %zu edges before mutation\n",
            argv[0], n_seeds, n_guards, map_popcount(tot_map));

    static uint8_t buf[FUZZ_MAX_INPUT];
    long found = 0;
    time_t t0 = time(NULL);
    for (long i = 0; i < runs; i++) {
        Input *src = &corpus[rnd_below(n_corpus)];
        size_t n = mutate(buf, src->buf, src->len);
        if (run_one(buf, n, cur_path[0] ? cur_path : NULL)) {
            corpus_add(buf, n, corpus_dir);
            found++;
        }
    }
    if (cur_path[0]) remove(cur_path);

    dprintf(g_report_fd, "[fuzz] %s: runs=%ld edges=%zu/%u corpus=%zu new=%ld elapsed=%.0fs\n",
            argv[0], runs, map_popcount(tot_map), n_guards, n_corpus, found,
            difftime(time(NULL), t0));
    return 0;
}
