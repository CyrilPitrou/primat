/* cli.c -- see cli.h. */
#include "cli.h"
#include "xalloc.h"
#include "api.h"
#include "cache.h"
#include "config.h"
#include "ini.h"
#include "mc.h"
#include "network_data.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>          /* timespec_get() for the "--- running time ---" line */
#include "compat_posix.h"  /* sys/stat.h + unistd.h/getcwd/mkdir, portable */
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

/* Matches "<prefix>*.txt" for one of the hash-named cache families, mirroring
 * primat.cache_utils._CACHE_PREFIXES:
 *   weak/   "nTOp_"            (thermal caches are "nTOp_thermal_*.txt" and
 *                               match the same prefix, as in Python)
 *   plasma/ "electron_thermo_" (hash-named since configurations stopped
 *                               evicting one another)
 * The fixed-name QED pressure tables are deliberately not matched by any
 * prefix: they cannot proliferate, so there is nothing to clean. */
static int is_cache_name(const char *name, const char *prefix)
{
    size_t plen = strlen(prefix), nlen = strlen(name);
    return strncmp(name, prefix, plen) == 0
        && nlen > 4 && strcmp(name + nlen - 4, ".txt") == 0;
}

/* Byte-for-byte twin of primat.credits.cli_credits_text()
 * (_CREDITS_CORE + _CREDITS_CLI_SUFFIX + "\n\n" + CITATION_BIBTEX), which is
 * what `primat --credits` prints. The paragraph breaks and the closing BibTeX
 * entry are part of that text, not decoration: users copy-paste the entry
 * rather than hand-formatting a reference from the arXiv link. Any edit to
 * primat/credits.py must be mirrored here. */
static void print_credits(void)
{
    fputs("primat is developed by Cyril Pitrou (https://www2.iap.fr/users/pitrou/) "
          "with features related to neutrino physics written by Julien Froustey.\n\n",
          stdout);
    fputs("The story started in the 1980s with BBN codes written by Elisabeth "
          "Vangioni and Alain Coc which eventually lead to 'ezbbn', a large "
          "nuclear network FORTRAN code whose nuclear rates tables were maintained "
          "by Alain Coc.\n",
          stdout);
    fputs("PRIMAT, initially a Mathematica code, was based on "
          "'ezbbn' with improved neutrino physics. It is now translated into a "
          "python code, but it also relies on a C backend to improve its "
          "performance.\n\n",
          stdout);
    fputs("For notebooks, examples and documentation, download the source code "
          "(https://github.com/CyrilPitrou/primat).\n",
          stdout);
    fputs("Please cite the publication (https://arxiv.org/abs/1801.08023) if "
          "you use it.\n\n",
          stdout);
    fputs("@article{Pitrou:2018cgg,\n"
          "    author = \"Pitrou, Cyril and Coc, Alain and Uzan, Jean-Philippe and Vangioni, Elisabeth\",\n"
          "    title = \"{Precision big bang nucleosynthesis with improved Helium-4 predictions}\",\n"
          "    eprint = \"1801.08023\",\n"
          "    archivePrefix = \"arXiv\",\n"
          "    primaryClass = \"astro-ph.CO\",\n"
          "    doi = \"10.1016/j.physrep.2018.04.005\",\n"
          "    journal = \"Phys. Rept.\",\n"
          "    volume = \"754\",\n"
          "    pages = \"1--66\",\n"
          "    year = \"2018\"\n"
          "}\n",
          stdout);
}

/* `--list-params`: every settable parameter with the value it currently holds
 * (i.e. its default, since this runs before any override is applied), so a
 * user can discover the full `--set KEY=VALUE` surface without reading
 * config.py.
 *
 * Deliberately WITHOUT the one-line descriptions Python's --list-params
 * prints: those live in config.py's inline comments, and Python parses them
 * out of the source rather than duplicating them (_default_params_comments).
 * Copying ~80 of them into this file would create a third place to keep in
 * sync -- exactly what generating the templates from DEFAULT_PARAMS exists to
 * prevent. The generated examples/run_basic.ini already carries the same
 * descriptions for the C side, so this points there instead. */
static void print_list_params(const CPRConfig *cfg)
{
    printf("# Every parameter settable via --set KEY=VALUE or an --ini file,\n"
           "# with its default value. One-line descriptions for each are in\n"
           "# primat-c/examples/run_basic.ini (or `primat --list-params`).\n");
    size_t n = cpr_config_field_count();
    for (size_t i = 0; i < n; i++) {
        const char *name = cpr_config_field_name(i);
        char value[CPR_PARAM_VAL_LEN];
        if (cpr_config_format_value(cfg, name, value, sizeof(value)) == 0)
            printf("%-32s = %s\n", name, value);
    }
}

/* Counts (and optionally deletes) the hash-named cache files of ONE family:
 * `subdir` is "weak" or "plasma", `prefix` the matching basename prefix.
 * Mirrors primat.cache_utils.list_cache_files/clear_cache restricted to a
 * single subdir, which is how cli.py reports the two counts separately. */
static int list_or_clear_cache(const CPRConfig *cfg, const char *subdir,
                                const char *prefix, int clear)
{
    /* Overlay-aware: the writable cache dir is cache_dir/<subdir> if set,
     * else <data_dir>/cache_plasma_weak/<subdir>. */
    char dir_path[CPR_PATH_BUF_LEN2];
    cpr_config_cache_write_dir(cfg, subdir, dir_path, sizeof(dir_path));

    DIR *d = opendir(dir_path);
    if (!d) {
        fprintf(stderr, "cannot open cache directory '%s'\n", dir_path);
        return 0;
    }
    int n = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (!is_cache_name(ent->d_name, prefix))
            continue;
        n++;
        if (clear) {
            char file_path[CPR_PATH_BUF_LEN2];
            snprintf(file_path, sizeof(file_path), "%s/%s", dir_path, ent->d_name);
            remove(file_path);
        }
    }
    closedir(d);
    return n;
}

/* Boolean PRIMATConfig flags that accept --flag / --no-flag pairs on the CLI,
 * mirroring Python cli.py's BooleanOptionalAction loop. */
static const char * const bool_flags[] = {
    "QED_corrections",
    "nuclear_qed_corrections",
    "radiative_corrections",
    "finite_mass_corrections",
    "thermal_corrections",
    "spectral_distortions",
    "output_time_evolution",
    "output_final_result",
    "output_background_evolution",
    "output_mc_samples",
    "output_mc_covariance",
    "output_mc_correlation",
    "show_progress",
    NULL
};

/* The 16 measured physical constants (primat.constants.OVERRIDABLE_CONSTANTS),
 * each exposed as its own --<name> VALUE flag so a sensitivity study reads the
 * same on both CLIs. The ten exact-by-definition constants are deliberately
 * absent: no config can carry them across the ABI. */
static const char * const CONSTANT_FLAGS[] = {
    "alphaem", "GF", "mZ", "me", "mn", "mp", "T0CMB", "gA", "Vud",
    "kappa_p", "kappa_n", "radproton", "ma", "He4Overma", "HOverma", "Neff_SM"
};
#define CONSTANT_FLAGS_N (sizeof(CONSTANT_FLAGS) / sizeof(CONSTANT_FLAGS[0]))

/* Returns the constant name if `arg` is "--<one of CONSTANT_FLAGS>", else NULL. */
static const char *constant_flag_name(const char *arg)
{
    if (arg[0] != '-' || arg[1] != '-') return NULL;
    for (size_t k = 0; k < CONSTANT_FLAGS_N; k++)
        if (strcmp(arg + 2, CONSTANT_FLAGS[k]) == 0) return CONSTANT_FLAGS[k];
    return NULL;
}

/* `--list-reactions`: every reaction of the configured network -- exactly the
 * bare names the p_<name>/delta_<name> rate-variation keys accept. That family
 * is per-network and unbounded, so it cannot appear in --list-params. Honours
 * network/amax, and skips the leading n<->p weak reaction, which carries no
 * rate table to vary. Mirrors cli.py's _print_list_reactions. */
static int name_cmp(const void *a, const void *b)
{
    return strcmp(*(const char * const *)a, *(const char * const *)b);
}

/* p_<rxn>/delta_<rxn> typo check -- the standalone-CLI half of
 * PRIMATConfig._warn_unknown_rate_variations. A key naming a reaction outside
 * the loaded network is a silent no-op otherwise, and reaction names are long
 * and underscore-heavy enough to be mistyped. Callers arriving through
 * primat/backend.py were already checked by PRIMATConfig, which is why this
 * lives here and not in cprimat_run: it must not report the same typo twice.
 * Runs only when the user actually passed such a key, so the usual path pays
 * no network load. Returns nonzero (having printed) iff strict_params makes
 * it fatal. */
static int check_rate_variation_keys(const CPRConfig *cfg)
{
    if (cfg->p_rxn.n == 0 && cfg->delta_rxn.n == 0) return 0;

    CPRNetworkDef net;
    char *err = NULL;
    if (cpr_load_network(cfg, "LT", NULL, 0, NULL, &net, &err)) {
        /* Not this check's job to report a broken network: the run below
         * fails on the same load with its own message. */
        free(err);
        return 0;
    }

    char bad[512];
    size_t used = 0;
    int n_bad = 0;
    for (int pass = 0; pass < 2; pass++) {
        const CPRRxnMap *map = pass ? &cfg->delta_rxn : &cfg->p_rxn;
        const char *prefix = pass ? "delta_" : "p_";
        for (size_t i = 0; i < map->n; i++) {
            int known = 0;
            for (size_t j = 0; j < net.n_reac && !known; j++)
                known = (strcmp(net.names[j], map->entries[i].name) == 0);
            if (known) continue;
            n_bad++;
            if (used < sizeof(bad) - 64)
                used += (size_t)snprintf(bad + used, sizeof(bad) - used,
                                         "%s'%s%s'", used ? ", " : "",
                                         prefix, map->entries[i].name);
        }
    }
    cpr_network_def_free(&net);
    if (n_bad == 0) return 0;

    /* Same wording as the Python warning/ValueError, per the output-parity
     * mandate; strict_params promotes it to a fatal error there too. */
    fprintf(stderr,
            "%s: PRIMATConfig: rate-variation key%s %s do%s not match any "
            "reaction in network '%s'; %s no effect on the run.%s\n",
            cfg->strict_params ? "error" : "warning",
            n_bad > 1 ? "s" : "", bad, n_bad > 1 ? "" : "es", cfg->network,
            n_bad > 1 ? "they have" : "it has",
            cfg->strict_params ? " [strict_params=True]" : "");
    return cfg->strict_params ? 1 : 0;
}

static int print_list_reactions(const CPRConfig *cfg)
{
    CPRNetworkDef net;
    char *err = NULL;
    if (cpr_load_network(cfg, "LT", NULL, 0, NULL, &net, &err)) {
        fprintf(stderr, "error: %s\n", err ? err : "cannot load network");
        free(err);
        return 2;
    }
    const char **names = CPR_XMALLOC(net.n_reac * sizeof(*names));
    size_t n = 0;
    for (size_t i = 0; i < net.n_reac; i++)
        if (strcmp(net.names[i], "n__p") != 0)
            names[n++] = net.names[i];
    qsort(names, n, sizeof(*names), name_cmp);

    printf("# %zu reactions in network '%s'", n, cfg->network);
    if (cfg->amax != -1) printf(" with amax=%d", cfg->amax);
    printf("\n");
    printf("# Vary any of them with --set p_<name>=<sigmas> (log-normal, in "
           "units of the\n# tabulated 1-sigma factor) or --set "
           "delta_<name>=<fraction> (additive).\n");
    for (size_t i = 0; i < n; i++)
        printf("%s\n", names[i]);

    free(names);
    cpr_network_def_free(&net);
    return 0;
}

static void usage(const char *prog)
{
    printf("usage: %s [-h] [--credits] [--version] [--list-params] [--list-reactions]\n"
           "          [--Omegabh2 VALUE] [--DeltaNeff VALUE] [--network NAME]\n"
           "          [--amax A] [--numerical_precision RTOL] [--munuOverTnu XI]\n"
           "          [--munuOverTnu_e XI_E] [--munuOverTnu_mu XI_MU] [--munuOverTnu_tau XI_TAU]\n"
           "          [--alphaem V] [--GF V] [--mZ V] [--me V] [--mn V] [--mp V]\n"
           "          [--T0CMB V] [--gA V] [--Vud V] [--kappa_p V] [--kappa_n V]\n"
           "          [--radproton V] [--ma V] [--He4Overma V] [--HOverma V] [--Neff_SM V]\n"
           "          [--output_file FILE] [--output_final_file FILE]\n"
           "          [--output_background_file FILE] [--output_mc_file_prefix PREFIX]\n"
           "          [--QED_corrections | --no-QED_corrections]\n"
           "          [--nuclear_qed_corrections | --no-nuclear_qed_corrections]\n"
           "          [--radiative_corrections | --no-radiative_corrections]\n"
           "          [--finite_mass_corrections | --no-finite_mass_corrections]\n"
           "          [--thermal_corrections | --no-thermal_corrections]\n"
           "          [--spectral_distortions | --no-spectral_distortions]\n"
           "          [--output_time_evolution | --no-output_time_evolution]\n"
           "          [--output_final_result | --no-output_final_result]\n"
           "          [--output_background_evolution | --no-output_background_evolution]\n"
           "          [--output_mc_samples | --no-output_mc_samples]\n"
           "          [--output_mc_covariance | --no-output_mc_covariance]\n"
           "          [--output_mc_correlation | --no-output_mc_correlation]\n"
           "          [--show_progress | --no-show_progress]\n"
           "          [--mc N] [--mc-seed SEED] [--mc-jobs N]\n"
           "          [--json] [--verbose] [--cache-info] [--cache-clear]\n"
           "          [--ini PATH] [--data_dir PATH] [--user_nuclear_dir PATH]\n"
           "          [--set KEY=VALUE ...]\n\n"
           "Run a Big Bang Nucleosynthesis computation with primat-c and print the\n"
           "resulting Neff/abundances.\n\n"
           "options:\n"
           "  -h, --help            Show this help message and exit.\n"
           "  --credits             Print the project credits and exit.\n"
           "  --version             Print the primat-c version and exit.\n"
           "  --list-params         Print every parameter settable via --set/--ini\n"
           "  --list-reactions      Print every reaction name of the selected\n"
           "                        --network/--amax (the p_<reaction>/\n"
           "                        delta_<reaction> rate-variation keys), then exit.\n"
           "                        with its default value, then exit. One-line\n"
           "                        descriptions are in examples/run_basic.ini.\n"
           "  --Omegabh2 VALUE      Baryon density Omega_b h^2 (default: 0.02242).\n"
           "  --DeltaNeff VALUE     Extra relativistic degrees of freedom on top of\n"
           "                        the SM neutrino sector (default: 0).\n"
           "  --network NAME        Nuclear reaction network used in the LT era\n"
           "                        (default: small). Built-in choices are 'small',\n"
           "                        'small_parthenope' and 'large', but any name for\n"
           "                        which data/nuclear/networks/<NAME>.txt exists is\n"
           "                        accepted.\n"
           "  --amax A              Drop reactions involving any nuclide with mass\n"
           "                        number > A (positive integer); applies to any\n"
           "                        --network. E.g. --network large --amax 8\n"
           "                        reproduces the old 'medium' network's 68 reactions.\n"
           "  --numerical_precision RTOL\n"
           "                        Relative tolerance passed to the ODE solver\n"
           "                        (default: 1e-7).\n"
           "  --munuOverTnu XI      Reduced neutrino chemical potential mu/T, the\n"
           "                        common default for all flavours (default: 0).\n"
           "  --munuOverTnu_e XI_E  Per-flavour ξ of nu_e; overrides --munuOverTnu\n"
           "                        for the electron neutrino, the only flavour that\n"
           "                        shifts the n<->p weak rates. (default: inherit).\n"
           "  --munuOverTnu_mu XI_MU\n"
           "                        Per-flavour ξ of nu_mu (gravitates only; default:\n"
           "                        inherit --munuOverTnu).\n"
           "  --munuOverTnu_tau XI_TAU\n"
           "                        Per-flavour ξ of nu_tau (gravitates only; default:\n"
           "                        inherit --munuOverTnu).\n"
           "  --alphaem V, --GF V, --mZ V, --me V, --mn V, --mp V, --T0CMB V,\n"
           "  --gA V, --Vud V, --kappa_p V, --kappa_n V, --radproton V, --ma V,\n"
           "  --He4Overma V, --HOverma V, --Neff_SM V\n"
           "                        The 16 measured physical constants, settable for\n"
           "                        sensitivity studies. --list-params prints each\n"
           "                        one's default. The remaining ten constants are\n"
           "                        exact by definition and cannot be set.\n"
           "  --output_file FILE    Write the full time-evolution TSV to FILE when\n"
           "                        --output_time_evolution is enabled.\n"
           "  --output_final_file FILE\n"
           "                        Write the final-abundance table to FILE when\n"
           "                        --output_final_result is enabled.\n"
           "  --output_background_file FILE\n"
           "                        Write the background time-evolution TSV to FILE\n"
           "                        when --output_background_evolution is enabled.\n"
           "  --output_mc_file_prefix PREFIX\n"
           "                        Filename stem for the Monte-Carlo output files\n"
           "                        written when --mc is used: PREFIX_samples.tsv /\n"
           "                        PREFIX_covariance.tsv / PREFIX_correlation.tsv,\n"
           "                        each gated by --output_mc_samples /\n"
           "                        --output_mc_covariance / --output_mc_correlation.\n"
           "  --QED_corrections, --no-QED_corrections\n"
           "                        QED interaction corrections to the EM plasma\n"
           "                        equation of state. (default: True).\n"
           "  --nuclear_qed_corrections, --no-nuclear_qed_corrections\n"
           "                        QED corrections to radiative-capture nuclear\n"
           "                        reaction rates (Pitrou & Pospelov 2020).\n"
           "                        (default: True).\n"
           "  --radiative_corrections, --no-radiative_corrections\n"
           "                        Coulomb + T=0 resummed radiative corrections to\n"
           "                        n<->p (CCR); if False, use Born approximation.\n"
           "                        (default: True).\n"
           "  --finite_mass_corrections, --no-finite_mass_corrections\n"
           "                        Finite-nucleon-mass (Fokker-Planck) correction\n"
           "                        to n<->p. (default: True).\n"
           "  --thermal_corrections, --no-thermal_corrections\n"
           "                        Finite-temperature radiative corrections to\n"
           "                        n<->p (CCRTh; Brown & Sawyer 2001). (default: True).\n"
           "  --spectral_distortions, --no-spectral_distortions\n"
           "                        Correct n<->p rates for non-Fermi-Dirac neutrino\n"
           "                        distributions. (default: True).\n"
           "  --output_time_evolution, --no-output_time_evolution\n"
           "                        Write the full time-evolution series (in-memory\n"
           "                        always; to disk if output_file is set).\n"
           "                        (default: False).\n"
           "  --output_final_result, --no-output_final_result\n"
           "                        Write the final results dict to output_final_file.\n"
           "                        (default: False).\n"
           "  --output_background_evolution, --no-output_background_evolution\n"
           "                        Write the cosmological background time series to\n"
           "                        disk. (default: False).\n"
           "  --output_mc_samples, --no-output_mc_samples\n"
           "                        Write --mc samples to\n"
           "                        <output_mc_file_prefix>_samples.tsv.\n"
           "                        (default: False).\n"
           "  --output_mc_covariance, --no-output_mc_covariance\n"
           "                        Write the --mc sample covariance matrix to\n"
           "                        <output_mc_file_prefix>_covariance.tsv.\n"
           "                        (default: False).\n"
           "  --output_mc_correlation, --no-output_mc_correlation\n"
           "                        Write the --mc sample correlation matrix to\n"
           "                        <output_mc_file_prefix>_correlation.tsv.\n"
           "                        (default: False).\n"
           "  --show_progress, --no-show_progress\n"
           "                        Print compact stderr progress indicators\n"
           "                        ('[primat]  HT.  MT.  LT.  done.' phase markers,\n"
           "                        '[MC] ...' sample counter) when --verbose is not\n"
           "                        used. (default: True).\n"
           "  --mc N                Run an N-sample Monte-Carlo nuclear-rate/tau_n\n"
           "                        uncertainty propagation and print each observable\n"
           "                        as 'value +/- sigma'. Uses all available CPU\n"
           "                        cores unless --mc-jobs says otherwise.\n"
           "  --mc-seed SEED        Base RNG seed for --mc (default: 0); sample i\n"
           "                        uses seed+i.\n"
           "  --mc-jobs N           Worker threads for --mc (default: -1, one per\n"
           "                        available core). Use a small N to leave cores\n"
           "                        free on a shared machine.\n"
           "  --json                Print the full results dict as JSON instead of a\n"
           "                        short summary.\n"
           "  --verbose             Enable internal progress messages (timings,\n"
           "                        cache hits, ...).\n"
           "  --cache-info          Print the number of cached hash-named files --\n"
           "                        n<->p weak-rate and e+- thermodynamic -- and exit,\n"
           "                        without running a solve.\n"
           "  --cache-clear         Delete every cached n<->p weak-rate and e+-\n"
           "                        thermodynamic file and exit, without running a\n"
           "                        solve. The cache is always safely regenerable\n"
           "                        (seconds per configuration for the plain rates and\n"
           "                        the electron-thermo tables, minutes for the thermal\n"
           "                        ones).\n"
           "  --ini PATH            Load parameters from an INI file (applied after\n"
           "                        defaults, before named flags and --set).\n"
           "  --data_dir PATH       Replace the entire data tree (NEVO/, nuclear/,\n"
           "                        csv/, cache_plasma_weak/) with PATH.\n"
           "                        Default: auto-detected from the executable location\n"
           "                        or CPRIMAT_DATA_DIR environment variable.\n"
           "  --user_nuclear_dir PATH\n"
           "                        Additive overlay for nuclear networks and rate\n"
           "                        tables only (primat/data/nuclear/ equivalent).\n"
           "                        Checked before the default tree; shipped networks\n"
           "                        remain accessible even when this is set.\n"
           "  --set KEY=VALUE       Set any CPRConfig parameter (including\n"
           "                        p_<reaction>/delta_<reaction> rate variations),\n"
           "                        e.g. --set T_end_MeV=1e-4. Repeatable; later\n"
           "                        values win.\n",
           prog);
}

/* S_ISDIR, not `st_mode & S_IFDIR`: S_IFDIR (0040000) is a *value* within the
 * S_IFMT field, not a standalone bit, so the bitwise test also accepts
 * S_IFBLK (0060000) and S_IFSOCK (0140000) -- a block device or socket passed
 * as --user_nuclear_dir would have been taken for a directory. */
static int path_is_dir(const char *path)
{
    struct stat st;
    return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
}

/* Best-effort absolute path to the running executable's own directory, so
 * the default data dir can be anchored to where `cprimat` itself lives
 * rather than to the caller's CWD (the old ".." default silently broke
 * whenever invoked from anywhere other than primat-c/, e.g. from
 * primat-c/build/ or the repo root). Returns 0 and fills
 * `out` on success, nonzero if the platform call fails (caller falls back
 * to a CWD-relative guess). */
static int executable_dir(char *out, size_t outsize)
{
    char exe_path[4096];
#if defined(__APPLE__)
    uint32_t sz = sizeof(exe_path);
    if (_NSGetExecutablePath(exe_path, &sz) != 0) return 1;
#elif defined(__linux__)
    ssize_t n = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
    if (n <= 0) return 1;
    exe_path[n] = '\0';
#else
    return 1;
#endif
    char *slash = strrchr(exe_path, '/');
    if (!slash) return 1;
    *slash = '\0';
    snprintf(out, outsize, "%s", exe_path);
    return 0;
}

/* Resolves the default data dir: CPRIMAT_DATA_DIR env var wins outright;
 * otherwise try "<exe_dir>/../primat/data" (works for `cprimat` run from
 * primat-c/, primat-c/build/, or any installed location with that sibling
 * layout); otherwise fall back to the legacy CWD-relative "../primat/data"
 * guess (works only when invoked with CWD == primat-c/). Does not require
 * the resolved directory to exist -- cpr_config_init_defaults reports a
 * clear "nuclides.csv not found" error downstream if it's wrong, and the
 * user can always override with --data_dir. */
static const char *default_data_dir(char *buf, size_t bufsize)
{
    const char *env = getenv("CPRIMAT_DATA_DIR");
    if (env) return env;

    char exe_dir[4096];
    if (executable_dir(exe_dir, sizeof(exe_dir)) == 0) {
        /* The binary normally lives in primat-c/build/cprimat, so the
         * sibling primat/ package is two levels up; also try one level up
         * in case cprimat was copied/symlinked directly into primat-c/. */
        snprintf(buf, bufsize, "%s/../../primat/data", exe_dir);
        if (path_is_dir(buf)) return buf;
        snprintf(buf, bufsize, "%s/../primat/data", exe_dir);
        if (path_is_dir(buf)) return buf;
    }
    snprintf(buf, bufsize, "../primat/data");
    return buf;
}

/* ---- Collected CLI + ini overrides ----
 *
 * Every override the user supplied is recorded in one CPRParamList and
 * forwarded verbatim to cpr_mc_uncertainty as base_params, because MC workers
 * rebuild their CPRConfig from *defaults plus this list* (mc.c's
 * worker_setup) rather than inheriting the main run's cfg. Anything applied to
 * cfg but missing from the list therefore vanishes from every sample: the
 * printed "value +/- sigma" would pair a central value from one model with a
 * sigma from another, silently. That is why the ini file (cpr_ini_load's
 * `collect` argument) feeds the same list, and why CPRParamList copies both
 * key and value instead of retaining pointers into argv or into
 * cpr_parse_literal's static scratch buffer.
 *
 * Apply the param to cfg and, on success, record it. A failure is fatal:
 * returns nonzero, and the caller exits (see cpr_config_set_by_name's
 * CPR_SET_* contract -- an unknown key is only a warning, and only while
 * strict_params is off). */
static int apply_param(CPRConfig *cfg, CPRParamList *cp,
                       const char *key, CPRParam val, const char *flag_label)
{
    char *set_err = NULL;
    int rc = cpr_config_set_by_name(cfg, key, val, &set_err);
    if (rc == CPR_SET_UNKNOWN_KEY && !cfg->strict_params) {
        /* Warn and ignore, mirroring PRIMATConfig's strict_params=False
         * default (and primat/backend.py's _c_params filter, which drops
         * keys unknown to both sides before they reach the extension). */
        fprintf(stderr, "%s: warning: %s\n", flag_label,
                set_err ? set_err : "unknown parameter key");
        free(set_err);
        return 0;
    }
    if (rc != CPR_SET_OK) {
        fprintf(stderr, "error: %s: %s%s\n", flag_label,
                set_err ? set_err : "could not set key",
                rc == CPR_SET_UNKNOWN_KEY ? " [strict_params=True]" : "");
        free(set_err);
        return 1;
    }
    free(set_err);
    cpr_paramlist_add(cp, key, val);
    return 0;
}

/* apply_param + "release everything and leave with status 2", the form every
 * call site below wants (2 = usage error, as for an unrecognised argument). */
#define APPLY_OR_FAIL(cfg, cp, key, val, label)                 \
    do {                                                        \
        if (apply_param((cfg), (cp), (key), (val), (label))) {  \
            cpr_paramlist_free(cp);                             \
            cpr_config_free(cfg);                               \
            return 2;                                           \
        }                                                       \
    } while (0)

/* ---- JSON output ---- */

/* Prints a JSON-safe string (escaping backslash and double-quote). */
static void print_json_str(const char *s)
{
    putchar('"');
    for (; *s; s++) {
        if (*s == '"' || *s == '\\') putchar('\\');
        putchar(*s);
    }
    putchar('"');
}

/* Emits the plain solve() result dict (no --mc): the observables, each guarded
 * exactly as main.py guards its dict key, plus the nested per-nuclide
 * "Y_final". `sep` is the running separator ("" before the first key, ",\n"
 * after); the updated value is returned so the caller can continue the object. */
static const char *print_json_results_body(const CPRResults *results,
                                           const char *sep)
{
    if (results->has_Neff) {
        printf("%s  \"Neff\": %.10g", sep, results->Neff); sep = ",\n";
    }
    printf("%s  \"YPBBN\": %.10g", sep, results->YPBBN);   sep = ",\n";
    printf("%s  \"YPCMB\": %.10g", sep, results->YPCMB);   sep = ",\n";
    printf("%s  \"He4oH\": %.10g", sep, results->He4oH);    sep = ",\n";
    printf("%s  \"DoH\": %.10g",   sep, results->DoH);      sep = ",\n";
    printf("%s  \"He3oH\": %.10g", sep, results->He3oH);    sep = ",\n";
    printf("%s  \"He3oHe4\": %.10g", sep, results->He3oHe4); sep = ",\n";
    printf("%s  \"Li7oH\": %.10g", sep, results->Li7oH);    sep = ",\n";
    if (results->has_Li6oLi7) {
        printf("%s  \"Li6oLi7\": %.10g", sep, results->Li6oLi7); sep = ",\n";
    }
    if (results->has_YCNO) {
        printf("%s  \"YCNO\": %.10g", sep, results->YCNO); sep = ",\n";
    }
    if (results->has_Omeganurel) {
        printf("%s  \"Omeganurel\": %.10g", sep, results->Omeganurel); sep = ",\n";
    }
    if (results->has_OneOverOmeganunr) {
        printf("%s  \"OneOverOmeganunr\": %.10g", sep, results->OneOverOmeganunr); sep = ",\n";
    }

    /* Per-nuclide final abundances. */
    if (results->n_nuclides > 0) {
        printf("%s  \"Y_final\": {", sep); sep = ",\n";
        for (size_t i = 0; i < results->n_nuclides; i++) {
            printf("%s\n    ", i > 0 ? "," : "");
            print_json_str(results->nuclide_names[i]);
            printf(": %.10g", results->Y_final[i]);
        }
        printf("\n  }");
    }
    return sep;
}

/* JSON payload, matching `primat --json` key-for-key.
 *
 * The two --mc/no---mc shapes are NOT the same dict plus extras, and the C
 * side must follow Python's lead in both:
 *
 *  - without --mc, the payload is the solve()'s result dict: the observables,
 *    plus a nested "Y_final" of every tracked nuclide;
 *  - with --mc, cli.py replaces it wholesale with MCResult.to_flat_dict(),
 *    which emits every MC quantity -- nuclides included -- as a FLAT
 *    top-level key alongside its "sigma_<name>", and therefore carries no
 *    "Y_final" (the nuclides are already top-level) and no Omeganurel /
 *    OneOverOmeganunr (not MC quantities). Both then get the "mc" sub-dict.
 *
 * Emitting the no-MC shape with sigmas bolted on, as this used to, meant a
 * script parsing `--json --mc` needed a different key list per binary. */
static void print_json(const CPRResults *results, const CPRMCResult *mc)
{
    printf("{\n");
    const char *sep = "";
    int have_mc = (mc && mc->n > 0);

    if (have_mc) {
        /* to_flat_dict(): every quantity's central, then its sigma_<name>,
         * in MC quantity order (observables first, then nuclides -- the order
         * cpr_cli_main built the `quantities` array in, which mirrors
         * backend.py's run_mc). */
        for (size_t i = 0; i < mc->n; i++) {
            const CPRMCQuantity *q = &mc->items[i];
            char sigma_name[48];
            printf("%s  ", sep);
            print_json_str(q->name);
            printf(": %.10g", q->central);
            sep = ",\n";
            snprintf(sigma_name, sizeof(sigma_name), "sigma_%s", q->name);
            printf("%s  ", sep);
            print_json_str(sigma_name);
            printf(": %.10g", q->std);
        }
    } else {
        sep = print_json_results_body(results, sep);
    }

    /* MC summary (central/mean/std per quantity; not the full sample array). */
    if (mc && mc->n > 0) {
        printf("%s  \"mc\": {", sep);
        for (size_t i = 0; i < mc->n; i++) {
            const CPRMCQuantity *q = &mc->items[i];
            printf("%s\n    ", i > 0 ? "," : "");
            print_json_str(q->name);
            printf(": {\"central\": %.10g, \"mean\": %.10g, \"std\": %.10g}",
                   q->central, q->mean, q->std);
        }
        printf("\n  }");
    }

    printf("\n}\n");
}

/* ---- Monte-Carlo output files (samples / covariance / correlation) ----
 *
 * The standalone C CLI writes the same three MC files as primat/cli.py, with
 * byte-identical header wording and value formatting:
 * <output_mc_file_prefix>_samples.tsv / _covariance.tsv /
 * _correlation.tsv, each gated by its own boolean flag. The samples file is a
 * straight port of primat.backend.dump_mc_samples; the two matrix files port
 * dump_mc_covariance/dump_mc_correlation, computing the ddof=1 sample
 * covariance/correlation here from the CPRMCResult's own per-quantity value
 * arrays. */

/* mkdir -p equivalent (mirrors nuclear_network.c's static mkdir_p and Python's
 * os.makedirs(exist_ok=True)); creates each '/'-delimited component in turn. */
static void cli_mkdir_p(const char *path)
{
    char buf[CPR_PATH_BUF_LEN2];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char *p = buf + 1; *p; p++) {
        if (*p == '/') { *p = '\0'; mkdir(buf, 0755); *p = '/'; }
    }
    mkdir(buf, 0755);
}

/* Resolve `path` to an absolute path (prepending the current working directory
 * when it is relative), mirroring os.path.abspath so the "[output] ..." lines
 * print the same absolute path as cli.py. No "." / ".." normalisation is done
 * (unnecessary for the simple prefix stems used here). */
static void cli_abspath(const char *path, char *out, size_t outsize)
{
    if (path[0] == '/') { snprintf(out, outsize, "%s", path); return; }
    char cwd[CPR_PATH_BUF_LEN];
    if (getcwd(cwd, sizeof(cwd)))
        snprintf(out, outsize, "%s/%s", cwd, path);
    else
        snprintf(out, outsize, "%s", path);
}

/* Wall-clock seconds, for the closing "--- running time: ... ---" line.
 *
 * timespec_get (C11, TIME_UTC) rather than clock(): clock() measures CPU time
 * summed over all threads, so a --mc run on 10 cores would report roughly ten
 * times the elapsed time. cli.py's counterpart is time.time(), i.e. wall
 * clock, and the two numbers have to be comparable. */
static double cli_wall_seconds(void)
{
    struct timespec ts;
    if (timespec_get(&ts, TIME_UTC) != TIME_UTC)
        return 0.0;
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Startup note for a data-tree override, the twin of config.py's
 * _rates_overlay_notice() (printed to stderr by cli.py for both keys). Same
 * wording, same two variants, same absolute-path quoting -- a run whose rate
 * tables come from somewhere other than the shipped tree must say so
 * identically on either backend. */
static void print_overlay_notice(const char *field, const char *path)
{
    if (!path || !path[0])
        return;
    char abs[CPR_PATH_BUF_LEN2];
    cli_abspath(path, abs, sizeof(abs));
    if (strcmp(field, "data_dir") == 0)
        fprintf(stderr, "[init-c] data_dir full-takeover data directory override: "
                        "entire data tree (NEVO/, nuclear/, csv/, cache_plasma_weak/) "
                        "replaced under '%s'.\n", abs);
    else
        fprintf(stderr, "[init-c] user_nuclear_dir additive nuclear overlay override: "
                        "nuclear networks and rate tables under '%s'.\n", abs);
}

/* Create the parent directory of `path` (if any), mirroring cli.py's
 * os.makedirs(os.path.dirname(abspath)). */
static void cli_ensure_parent_dir(const char *path)
{
    char dir[CPR_PATH_BUF_LEN2];
    snprintf(dir, sizeof(dir), "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) { *slash = '\0'; if (dir[0]) cli_mkdir_p(dir); }
}

/* The ddof=1 sample covariance kernel is cpr_mc_sample_cov (mc.h), shared with
 * the covariance/correlation matrix builders below and unit-tested in
 * tests/unit/test_mc.c. */

/* Port of primat.backend.dump_mc_samples: header = tab-joined quantity names,
 * then one %.10e row per sample. */
static int mc_write_samples(const char *path, const CPRMCResult *mc, int num_mc)
{
    cli_ensure_parent_dir(path);
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "warning: cannot write MC samples to %s\n", path); return 1; }
    for (size_t j = 0; j < mc->n; j++)
        fprintf(f, "%s%s", j ? "\t" : "", mc->items[j].name);
    fputc('\n', f);
    for (int s = 0; s < num_mc; s++) {
        for (size_t j = 0; j < mc->n; j++)
            fprintf(f, "%s%.10e", j ? "\t" : "", mc->items[j].values[s]);
        fputc('\n', f);
    }
    fclose(f);
    return 0;
}

/* Port of dump_mc_covariance (is_corr=0) / dump_mc_correlation (is_corr=1):
 * line 1 = a '#' comment (N, seed, estimator convention, byte-identical to the
 * Python writers); line 2 = "quantity\t<name>..."; then one name-first %.10e
 * row per quantity. For correlation, R[i,j]=C[i,j]/(std_i*std_j) with a unit
 * diagonal and NaN off-diagonal for any zero-variance quantity (mirrors
 * MCResult.corr). */
static int mc_write_matrix(const char *path, const CPRMCResult *mc,
                           int num_mc, int seed, int is_corr)
{
    cli_ensure_parent_dir(path);
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "warning: cannot write MC %s to %s\n",
                is_corr ? "correlation" : "covariance", path);
        return 1;
    }
    if (is_corr)
        fprintf(f, "# Correlation matrix of the N=%d primat MC samples (seed=%d): "
                   "R[i,j] = Pearson correlation (ddof=1) of quantities i and j; "
                   "unit diagonal.\n", num_mc, seed);
    else
        fprintf(f, "# Covariance matrix of the N=%d primat MC samples (seed=%d): "
                   "C[i,j] = sample covariance (ddof=1) of quantities i and j.\n",
                num_mc, seed);
    fputs("quantity", f);
    for (size_t j = 0; j < mc->n; j++) fprintf(f, "\t%s", mc->items[j].name);
    fputc('\n', f);
    for (size_t i = 0; i < mc->n; i++) {
        fputs(mc->items[i].name, f);
        for (size_t j = 0; j < mc->n; j++) {
            double v;
            if (is_corr) {
                if (i == j) {
                    v = 1.0;   /* unit diagonal by convention */
                } else {
                    double si = mc->items[i].std, sj = mc->items[j].std;
                    v = (si == 0.0 || sj == 0.0) ? NAN
                        : cpr_mc_sample_cov(mc->items[i].values, mc->items[j].values, num_mc)
                          / (si * sj);
                }
            } else {
                v = cpr_mc_sample_cov(mc->items[i].values, mc->items[j].values, num_mc);
            }
            fprintf(f, "\t%.10e", v);
        }
        fputc('\n', f);
    }
    fclose(f);
    return 0;
}

/* Print the 4x4 correlation and covariance matrices of the four main BBN
 * products (YPBBN/DoH/He3oHe4/Li7oH), aligned; byte-for-byte identical layout
 * to primat/cli.py's _print_mc_matrices (verbose/output-parity). Only the
 * products this network produced are shown; nothing prints if fewer than two
 * are present. */
static void print_mc_matrices(const CPRMCResult *mc, int num_mc)
{
    static const char *main4[] = {"YPBBN", "DoH", "He3oHe4", "Li7oH"};
    size_t idx[4]; size_t m = 0;
    for (size_t k = 0; k < 4; k++) {
        size_t ix = cpr_mc_result_index(mc, main4[k]);
        if (ix < mc->n) idx[m++] = ix;
    }
    if (m < 2) return;
    /* A sample covariance/correlation needs >= 2 samples; skip for a
     * single-sample run (mirrors cli.py's _print_mc_matrices guard). */
    if (num_mc < 2) return;

    /* Title: comma-joined present labels (", " separator, matching Python). */
    char title[80]; title[0] = '\0';
    for (size_t k = 0; k < m; k++) {
        strncat(title, mc->items[idx[k]].name, sizeof(title) - strlen(title) - 1);
        if (k + 1 < m) strncat(title, ", ", sizeof(title) - strlen(title) - 1);
    }

    /* Correlation: 8-wide row labels, 9-wide %9.3f value columns.
     * Leading blank line separates the matrix block from the preceding
     * per-observable summary for readability (mirrors cli.py's
     * _print_mc_matrices, per the byte-for-byte output-parity mandate). */
    putchar('\n');
    printf("Correlation matrix (%s):\n", title);
    printf("%8s", "");
    for (size_t k = 0; k < m; k++) printf("%9s", mc->items[idx[k]].name);
    putchar('\n');
    for (size_t i = 0; i < m; i++) {
        printf("%8s", mc->items[idx[i]].name);
        for (size_t j = 0; j < m; j++) {
            double v;
            if (i == j) {
                v = 1.0;
            } else {
                double si = mc->items[idx[i]].std, sj = mc->items[idx[j]].std;
                v = (si == 0.0 || sj == 0.0) ? NAN
                    : cpr_mc_sample_cov(mc->items[idx[i]].values,
                                  mc->items[idx[j]].values, num_mc) / (si * sj);
            }
            printf("%9.3f", v);
        }
        putchar('\n');
    }
    /* Covariance: same 8-wide labels, 13-wide %13.3e value columns.
     * Leading blank line, as for the correlation block above. */
    putchar('\n');
    printf("Covariance matrix (%s):\n", title);
    printf("%8s", "");
    for (size_t k = 0; k < m; k++) printf("%13s", mc->items[idx[k]].name);
    putchar('\n');
    for (size_t i = 0; i < m; i++) {
        printf("%8s", mc->items[idx[i]].name);
        for (size_t j = 0; j < m; j++) {
            double v = cpr_mc_sample_cov(mc->items[idx[i]].values,
                                   mc->items[idx[j]].values, num_mc);
            printf("%13.3e", v);
        }
        putchar('\n');
    }
}

/* ---- Plain-text report (mirrors cli.py's default output) ---- */

static void print_plain(const CPRConfig *cfg, const CPRResults *results,
                        const CPRMCResult *mc, int mc_n, double elapsed_s)
{
    const char *sep = "────────────────────────────────────────────────────";
    char header[80];
    snprintf(header, sizeof(header), "PRIMAT results at T = %g MeV", cfg->T_end_MeV);
    printf("%s\n", sep);
    /* Python centres with f"{header:^52}", which pads on BOTH sides -- the
     * trailing run of spaces is part of the line. Reproduced here (rather than
     * left-padding only) so the two CLIs' output is byte-identical, as the
     * output-parity mandate requires. */
    int total_pad = 52 - (int)strlen(header);
    if (total_pad < 0) total_pad = 0;
    int left_pad = total_pad / 2;               /* str.center: extra space right */
    printf("%*s%s%*s\n", left_pad, "", header, total_pad - left_pad, "");
    printf("%s\n", sep);

/* Helper: if mc has this quantity, append " +/- std"; else nothing. */
#define MC_STD(name, fmt) do { \
    if (mc) { \
        size_t idx = cpr_mc_result_index(mc, name); \
        if (idx < mc->n) printf(" +/- " fmt, mc->items[idx].std); \
    } \
} while (0)

    if (results->has_Neff) {
        printf("Neff       = %.8f", results->Neff);
        MC_STD("Neff", "%.8f");
        putchar('\n');
    }
    printf("YP (BBN)   = %.8f", results->YPBBN);
    MC_STD("YPBBN", "%.8f"); putchar('\n');
    printf("YP (CMB)   = %.8f", results->YPCMB);
    MC_STD("YPCMB", "%.8f"); putchar('\n');
    printf("He4/H      = %.7e", results->He4oH);
    MC_STD("He4oH", "%.7e"); putchar('\n');
    printf("D/H        = %.7e", results->DoH);
    MC_STD("DoH", "%.7e"); putchar('\n');
    printf("He3/H      = %.7e", results->He3oH);
    MC_STD("He3oH", "%.7e"); putchar('\n');
    printf("He3/He4    = %.7e", results->He3oHe4);
    MC_STD("He3oHe4", "%.7e"); putchar('\n');
    printf("Li7/H      = %.6e", results->Li7oH);
    MC_STD("Li7oH", "%.6e"); putchar('\n');
    if (results->has_Li6oLi7) {
        printf("Li6/Li7    = %.6e", results->Li6oLi7);
        MC_STD("Li6oLi7", "%.6e"); putchar('\n');
    }
    if (results->has_YCNO) {
        printf("CNO (mass) = %.6e", results->YCNO);
        MC_STD("YCNO", "%.6e"); putchar('\n');
    }
#undef MC_STD

    if (mc) {
        /* Joint uncertainty of the four main products: 4x4 correlation +
         * covariance matrices (same layout as cli.py's _print_mc_matrices). */
        print_mc_matrices(mc, mc_n);
    }
    /* Closing line, matching cli.py's final print in both the plain and the
     * --mc case. (The C CLI used to print "--- Monte-Carlo: N samples ---"
     * here instead, and nothing at all without --mc; neither line exists on
     * the Python side, so the two outputs could not be diffed.) */
    printf("--- running time: %.2f seconds ---\n", elapsed_s);
}

int cpr_cli_main(int argc, char **argv)
{
    char data_dir_buf[CPR_PATH_BUF_LEN];
    const char *data_dir = default_data_dir(data_dir_buf, sizeof(data_dir_buf));
    const char *custom_nuclear_dir = NULL;
    const char *ini_path = NULL;
    int cache_info = 0, cache_clear = 0, credits = 0, version = 0;
    int list_params = 0;
    int list_reactions = 0;
    int do_json = 0;
    int mc_n = 0, mc_seed = 0, mc_jobs = -1;   /* -1 = one worker per core */
    /* mc_n == 0 also means "no --mc at all", so a separate flag is needed to
     * tell an omitted --mc from an explicit "--mc 0" (which is an error, as in
     * cli.py: a sigma needs at least 2 samples). */
    int mc_given = 0;

    /* --data_dir, --user_nuclear_dir and --ini must be known before
     * cpr_config_init_defaults runs (the first picks the data directory;
     * the others are applied after defaults), so scan for them first;
     * everything else is applied in a second pass, in the same precedence
     * order as cli.py: defaults, then .ini, then named flags, then --set
     * (later wins). */
    int data_dir_given = 0;   /* only announce an override the user asked for */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--data_dir") == 0 && i + 1 < argc) {
            data_dir = argv[++i];
            data_dir_given = 1;
        } else if (strcmp(argv[i], "--user_nuclear_dir") == 0 && i + 1 < argc) {
            custom_nuclear_dir = argv[++i];
        } else if (strcmp(argv[i], "--ini") == 0 && i + 1 < argc) {
            ini_path = argv[++i];
        } else if (strcmp(argv[i], "--cache-info") == 0) {
            cache_info = 1;
        } else if (strcmp(argv[i], "--cache-clear") == 0) {
            cache_clear = 1;
        } else if (strcmp(argv[i], "--credits") == 0) {
            credits = 1;
        } else if (strcmp(argv[i], "--version") == 0) {
            version = 1;
        } else if (strcmp(argv[i], "--list-params") == 0) {
            list_params = 1;
        } else if (strcmp(argv[i], "--list-reactions") == 0) {
            list_reactions = 1;
        } else if (strcmp(argv[i], "--json") == 0) {
            do_json = 1;
        } else if (strcmp(argv[i], "--mc") == 0 && i + 1 < argc) {
            mc_n = atoi(argv[++i]);
            mc_given = 1;
        } else if (strcmp(argv[i], "--mc-seed") == 0 && i + 1 < argc) {
            mc_seed = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--mc-jobs") == 0 && i + 1 < argc) {
            mc_jobs = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        }
    }

    /* Reject a non-positive explicit --mc before any setup: the sampler sizes
     * its buffers from num_mc, so a negative count used to abort the process in
     * CPR_XMALLOC with an "out of memory (1.8e19 bytes)" message naming neither
     * the flag nor the mistake. Mirrors cli.py's parser.error. */
    if (mc_given && mc_n < 1) {
        fprintf(stderr, "error: --mc must be >= 1 (got %d); a sigma needs at "
                        "least 2 samples.\n", mc_n);
        return 2;
    }

    if (version) {
        printf("primat-c %s\n", CPRIMAT_VERSION);
        return 0;
    }
    if (credits) {
        print_credits();
        return 0;
    }

    CPRConfig cfg;
    char *err = NULL;
    if (cpr_config_init_defaults(&cfg, data_dir, &err)) {
        fprintf(stderr, "error: %s\n", err);
        free(err);
        return 1;
    }

    if (list_params) {
        print_list_params(&cfg);
        cpr_config_free(&cfg);
        return 0;
    }

    if (custom_nuclear_dir) {
        if (!path_is_dir(custom_nuclear_dir)) {
            fprintf(stderr, "--user_nuclear_dir: '%s' is not a directory\n", custom_nuclear_dir);
            cpr_config_free(&cfg);
            return 2;
        }
        free(cfg.user_nuclear_dir);
        cfg.user_nuclear_dir = strdup(custom_nuclear_dir);
    }

    if (cache_info || cache_clear) {
        /* Both hash-named families, broken down per tree so the user can see
         * which one is accumulating -- same wording as cli.py. */
        int n_weak = list_or_clear_cache(&cfg, "weak", "nTOp_", cache_clear);
        int n_plasma = list_or_clear_cache(&cfg, "plasma", "electron_thermo_",
                                            cache_clear);
        char wdir[CPR_PATH_BUF_LEN2], pdir[CPR_PATH_BUF_LEN2];
        cpr_config_cache_write_dir(&cfg, "weak", wdir, sizeof(wdir));
        cpr_config_cache_write_dir(&cfg, "plasma", pdir, sizeof(pdir));
        if (cache_clear)
            printf("Removed %d cached file(s): %d weak-rate from %s/, "
                   "%d electron-thermo from %s/.\n",
                   n_weak + n_plasma, n_weak, wdir, n_plasma, pdir);
        else {
            printf("%d cached weak-rate file(s) in %s/.\n", n_weak, wdir);
            printf("%d cached electron-thermo file(s) in %s/.\n", n_plasma, pdir);
        }
        cpr_config_free(&cfg);
        return 0;
    }

    /* Collect all user-supplied overrides here; forwarded to MC workers.
     * Declared BEFORE the ini load so the ini's own keys land in it too --
     * the workers rebuild from defaults + this list, so anything missing from
     * it is silently absent from every MC sample (see apply_param's comment). */
    CPRParamList cp;
    memset(&cp, 0, sizeof(cp));

    /* Record user_nuclear_dir in base_params so MC workers inherit it. */
    if (custom_nuclear_dir)
        cpr_paramlist_add(&cp, "user_nuclear_dir",
                          (CPRParam){CPR_STRING, .v.s = custom_nuclear_dir});

    if (ini_path) {
        if (cpr_ini_load(&cfg, ini_path, &cp, &err)) {
            fprintf(stderr, "error: %s\n", err);
            free(err);
            cpr_paramlist_free(&cp);
            cpr_config_free(&cfg);
            return 1;
        }
    }

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        int has_val = (i + 1 < argc);

        /* Already handled in first pass or not a solver param. */
        if (strcmp(a, "--data_dir") == 0 || strcmp(a, "--user_nuclear_dir") == 0
            || strcmp(a, "--ini") == 0) { i++; continue; }
        if (strcmp(a, "--cache-info") == 0 || strcmp(a, "--cache-clear") == 0
            || strcmp(a, "--credits") == 0 || strcmp(a, "--version") == 0
            || strcmp(a, "--list-params") == 0 || strcmp(a, "--json") == 0
            || strcmp(a, "--list-reactions") == 0
            || strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0) continue;
        if (strcmp(a, "--mc") == 0 || strcmp(a, "--mc-seed") == 0
            || strcmp(a, "--mc-jobs") == 0) { i++; continue; }

        /* ---- The 16 measured physical constants ---- */
        const char *const_name = constant_flag_name(a);
        if (const_name && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, const_name, p, a);
            continue;
        }

        /* ---- Simple scalar flags (string or numeric) ---- */
        if (strcmp(a, "--Omegabh2") == 0 && has_val) {
            CPRParam p = {CPR_DOUBLE, .v.d = atof(argv[++i])};
            APPLY_OR_FAIL(&cfg, &cp, "Omegabh2", p, "--Omegabh2");
        } else if (strcmp(a, "--DeltaNeff") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "DeltaNeff", p, "--DeltaNeff");
        } else if (strcmp(a, "--network") == 0 && has_val) {
            CPRParam p = {CPR_STRING, .v.s = argv[++i]};
            APPLY_OR_FAIL(&cfg, &cp, "network", p, "--network");
        } else if (strcmp(a, "--amax") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "amax", p, "--amax");
        } else if (strcmp(a, "--numerical_precision") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "numerical_precision", p, "--numerical_precision");
        } else if (strcmp(a, "--munuOverTnu") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "munuOverTnu", p, "--munuOverTnu");
        } else if (strcmp(a, "--munuOverTnu_e") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "munuOverTnu_e", p, "--munuOverTnu_e");
        } else if (strcmp(a, "--munuOverTnu_mu") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "munuOverTnu_mu", p, "--munuOverTnu_mu");
        } else if (strcmp(a, "--munuOverTnu_tau") == 0 && has_val) {
            CPRParam p = cpr_parse_literal(argv[++i]);
            APPLY_OR_FAIL(&cfg, &cp, "munuOverTnu_tau", p, "--munuOverTnu_tau");
        } else if (strcmp(a, "--verbose") == 0) {
            CPRParam p = {CPR_BOOL, .v.b = 1};
            APPLY_OR_FAIL(&cfg, &cp, "verbose", p, "--verbose");

        /* ---- Output file paths ---- */
        } else if (strcmp(a, "--output_file") == 0 && has_val) {
            CPRParam p = {CPR_STRING, .v.s = argv[++i]};
            APPLY_OR_FAIL(&cfg, &cp, "output_file", p, "--output_file");
        } else if (strcmp(a, "--output_final_file") == 0 && has_val) {
            CPRParam p = {CPR_STRING, .v.s = argv[++i]};
            APPLY_OR_FAIL(&cfg, &cp, "output_final_file", p, "--output_final_file");
        } else if (strcmp(a, "--output_background_file") == 0 && has_val) {
            CPRParam p = {CPR_STRING, .v.s = argv[++i]};
            APPLY_OR_FAIL(&cfg, &cp, "output_background_file", p, "--output_background_file");
        } else if (strcmp(a, "--output_mc_file_prefix") == 0 && has_val) {
            CPRParam p = {CPR_STRING, .v.s = argv[++i]};
            APPLY_OR_FAIL(&cfg, &cp, "output_mc_file_prefix", p, "--output_mc_file_prefix");

        /* ---- Boolean --flag / --no-flag pairs ---- */
        } else if (strncmp(a, "--", 2) == 0) {
            /* Check --no-<flag> first (longer prefix), then --<flag>. */
            int matched = 0;
            for (int fi = 0; bool_flags[fi]; fi++) {
                char pos_flag[64], neg_flag[70];
                snprintf(pos_flag, sizeof(pos_flag), "--%s", bool_flags[fi]);
                snprintf(neg_flag, sizeof(neg_flag), "--no-%s", bool_flags[fi]);
                if (strcmp(a, neg_flag) == 0) {
                    CPRParam p = {CPR_BOOL, .v.b = 0};
                    APPLY_OR_FAIL(&cfg, &cp, bool_flags[fi], p, neg_flag);
                    matched = 1; break;
                } else if (strcmp(a, pos_flag) == 0) {
                    CPRParam p = {CPR_BOOL, .v.b = 1};
                    APPLY_OR_FAIL(&cfg, &cp, bool_flags[fi], p, pos_flag);
                    matched = 1; break;
                }
            }
            if (!matched) {
                /* ---- --set KEY=VALUE ---- */
                if (strcmp(a, "--set") == 0 && has_val) {
                    const char *entry = argv[++i];
                    const char *eq = strchr(entry, '=');
                    if (!eq) {
                        fprintf(stderr, "--set %s: expected KEY=VALUE\n", entry);
                        cpr_paramlist_free(&cp);
                        cpr_config_free(&cfg);
                        return 2;
                    }
                    /* An empty value is not a Python literal either: cli.py's
                     * ast.literal_eval("") raises and becomes a parser.error.
                     * Without this, `--set network=` set the network to "" and
                     * failed much later on a nonsensical ".txt" open. */
                    if (eq[1] == '\0') {
                        fprintf(stderr, "--set %s: expected KEY=VALUE "
                                        "(the value is empty)\n", entry);
                        cpr_paramlist_free(&cp);
                        cpr_config_free(&cfg);
                        return 2;
                    }
                    /* The key half of argv's "KEY=VALUE" is not NUL-terminated
                     * at the '=', so it needs its own copy; cpr_paramlist_add
                     * then copies it again into the list's own arena. */
                    char key[CPR_PARAM_KEY_LEN];
                    size_t klen = (size_t)(eq - entry);
                    if (klen >= sizeof(key)) klen = sizeof(key) - 1;
                    memcpy(key, entry, klen);
                    key[klen] = '\0';
                    CPRParam p = cpr_parse_literal(eq + 1);
                    APPLY_OR_FAIL(&cfg, &cp, key, p, entry);
                } else {
                    fprintf(stderr, "unrecognized argument: %s\n", a);
                    usage(argv[0]);
                    cpr_paramlist_free(&cp);
                    cpr_config_free(&cfg);
                    return 2;
                }
            }
        } else {
            fprintf(stderr, "unrecognized argument: %s\n", a);
            usage(argv[0]);
            cpr_paramlist_free(&cp);
            cpr_config_free(&cfg);
            return 2;
        }
    }

    if (cpr_config_validate(&cfg, &err)) {
        fprintf(stderr, "error: %s\n", err);
        free(err);
        cpr_paramlist_free(&cp);
        cpr_config_free(&cfg);
        return 1;
    }

    if (list_reactions) {
        int rc = print_list_reactions(&cfg);
        cpr_paramlist_free(&cp);
        cpr_config_free(&cfg);
        return rc;
    }

    if (check_rate_variation_keys(&cfg)) {
        cpr_paramlist_free(&cp);
        cpr_config_free(&cfg);
        /* 2, the exit status both CLIs already use for every fatal
         * configuration error (a bad range, an unknown network). */
        return 2;
    }

    /* Startup note for an overlay/takeover data directory, byte-identical to
     * cli.py's `print(_rates_overlay_notice(key, ...), file=sys.stderr)` --
     * a run reading rate tables from somewhere other than the shipped tree
     * should say so, on both backends. Printed after validation so a bad path
     * is reported as an error rather than announced as if it worked. */
    if (custom_nuclear_dir || cfg.user_nuclear_dir)
        print_overlay_notice("user_nuclear_dir", cfg.user_nuclear_dir);
    if (data_dir_given)
        print_overlay_notice("data_dir", cfg.data_dir);

    /* Wall clock for the closing "--- running time: X seconds ---" line,
     * matching cli.py's start_time (set just before the solve, so table
     * loading counts but argument parsing does not). */
    double t_start = cli_wall_seconds();

    CPRResults results;
    if (cprimat_run(&cfg, NULL, &results, &err)) {
        fprintf(stderr, "error: %s\n", err);
        free(err);
        cpr_paramlist_free(&cp);
        cpr_config_free(&cfg);
        return 1;
    }

    /* ---- Optional Monte-Carlo uncertainty propagation ---- */
    CPRMCResult mc_result;
    memset(&mc_result, 0, sizeof(mc_result));
    CPRMCResult *mc = NULL;

    if (mc_n > 0) {
        /* MC quantity set: the standard observables (mirroring
         * _DEFAULT_MC_OBSERVABLES in primat/backend.py, filtered to those this
         * network actually produces -- e.g. "Li6oLi7"/"YCNO" only for large)
         * FIRST, then every tracked nuclide's final Y, in nuclide_names order.
         * This matches the order primat/backend.py's run_mc builds
         * (`quantities_plus_observables + nuclides`), so the samples/covariance/
         * correlation files this CLI writes have the same columns, in the same
         * order, as the Python CLI's -- covering all MC quantities, not just
         * the observables. cpr_mc_uncertainty resolves each
         * nuclide name via cpr_results_get_quantity's Y_final fallback. */
        const char *all_quantities[] = {
            "Neff", "YPBBN", "YPCMB", "He4oH", "DoH", "He3oH", "He3oHe4",
            "Li7oH", "Li6oLi7", "YCNO"
        };
        size_t n_obs = sizeof(all_quantities)/sizeof(all_quantities[0]);
        size_t n_q = 0;
        /* observables (present) + every nuclide; observables never collide with
         * nuclide names (ratios vs element names), so no dedup is needed. */
        const char **quantities = CPR_XMALLOC((n_obs + results.n_nuclides)
                                         * sizeof(*quantities));
        for (size_t qi = 0; qi < n_obs; qi++) {
            int found = 0;
            cpr_results_get_quantity(&results, all_quantities[qi], &found);
            if (found) quantities[n_q++] = all_quantities[qi];
        }
        for (size_t k = 0; k < results.n_nuclides; k++)
            quantities[n_q++] = results.nuclide_names[k];

        if (cpr_mc_uncertainty(mc_n, quantities, n_q,
                               data_dir,
                               cp.items, cp.n,
                               mc_seed, mc_jobs, NULL,
                               NULL, NULL, 0, cfg.show_progress,
                               &mc_result, &err)) {
            /* Fatal, as on the Python side: cli.py lets run_mc's exception
             * propagate, so `primat --mc` exits non-zero. Printing the central
             * results and exiting 0 (as this used to) tells a script the run
             * succeeded while silently omitting every sigma it asked for. */
            fprintf(stderr, "error: MC: %s\n", err);
            free(err);
            free(quantities);
            cprimat_results_free(&results);
            cpr_paramlist_free(&cp);
            cpr_config_free(&cfg);
            return 1;
        }
        mc = &mc_result;
        free(quantities);
    } else {
        /* output_mc_samples/output_mc_covariance/output_mc_correlation only
         * have an effect inside the `if (mc)` file-writing block below (a
         * CPRMCResult from cpr_mc_uncertainty is what they are dumped from);
         * without --mc there is no CPRMCResult, so any of these flags being
         * set is silently a no-op unless we flag it here. Mirrors cli.py's
         * equivalent warning. */
        const char *requested[3];
        int n_requested = 0;
        if (cfg.output_mc_samples)     requested[n_requested++] = "output_mc_samples";
        if (cfg.output_mc_covariance)  requested[n_requested++] = "output_mc_covariance";
        if (cfg.output_mc_correlation) requested[n_requested++] = "output_mc_correlation";
        if (n_requested > 0) {
            fprintf(stderr, "warning: ");
            for (int k = 0; k < n_requested; k++)
                fprintf(stderr, "%s%s", k ? ", " : "", requested[k]);
            fprintf(stderr, " set but --mc was not passed; "
                            "no MC output file(s) will be written.\n");
        }
    }

    /* ---- Output ---- */
    double elapsed_s = cli_wall_seconds() - t_start;
    if (do_json) {
        print_json(&results, mc);
    } else {
        print_plain(&cfg, &results, mc, mc_n, elapsed_s);
    }

    /* ---- Optional MC output files (samples / covariance / correlation) ----
     * All three share the output_mc_file_prefix stem, each gated by its own
     * boolean, mirroring cli.py's writer block (byte-identical file format). */
    if (mc) {
        const char *prefix = cfg.output_mc_file_prefix;
        char rel[CPR_PATH_BUF_LEN], path[CPR_PATH_BUF_LEN2];
        if (cfg.output_mc_samples) {
            snprintf(rel, sizeof(rel), "%s_samples.tsv", prefix);
            cli_abspath(rel, path, sizeof(path));
            if (mc_write_samples(path, mc, mc_n) == 0)
                printf("[output] MC samples (%d sample%s) written to %s\n",
                       mc_n, mc_n == 1 ? "" : "s", path);
        }
        if (cfg.output_mc_covariance) {
            snprintf(rel, sizeof(rel), "%s_covariance.tsv", prefix);
            cli_abspath(rel, path, sizeof(path));
            if (mc_write_matrix(path, mc, mc_n, mc_seed, 0) == 0)
                printf("[output] MC covariance matrix written to %s\n", path);
        }
        if (cfg.output_mc_correlation) {
            snprintf(rel, sizeof(rel), "%s_correlation.tsv", prefix);
            cli_abspath(rel, path, sizeof(path));
            if (mc_write_matrix(path, mc, mc_n, mc_seed, 1) == 0)
                printf("[output] MC correlation matrix written to %s\n", path);
        }
    }

    if (mc) cpr_mc_result_free(&mc_result);
    cprimat_results_free(&results);
    cpr_paramlist_free(&cp);
    cpr_config_free(&cfg);
    return 0;
}
