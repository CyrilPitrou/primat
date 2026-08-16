/* plasma.c -- see plasma.h. */
#include "plasma.h"
#include "xalloc.h"
#include "constants.h"
#include "qed_pressure.h"
#include "quad.h"
#include "table_io.h"
#include "cache.h"
#include "log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compat_posix.h"  /* sys/stat.h + mkdir, portable across POSIX & MSVC */

/* Below Tg = me / _ELEC_THERMO_LOWT_RATIO the e+- number density is
 * Boltzmann-suppressed by exp(-me/Tg) < exp(-30) ~ 1e-13 relative to
 * photons, so all four e+- quantities are exactly 0 (Phys. Rep. App. A.2
 * non-relativistic limit) -- avoids integrating a negligible, potentially
 * slow tail. Matches plasma.py's _ELEC_THERMO_LOWT_RATIO. */
#define ELEC_THERMO_LOWT_RATIO 30.0

/* High-T limit of spl/Tg^3 (photon + e+-): Phys. Rep. Eq. 25d/26d give
 * sbar_g = 4 pi^2/45, sbar_e+- = 7/8 sbar_g each, so
 * sbar_pl = sbar_g + 2*(7/8)*sbar_g = (11/4) * 4pi^2/45 = 11 pi^2/45. */
#define SIGMA_INF (11.0 * M_PI * M_PI / 45.0)

/* v2: the e+- quadratures moved from an absolute to a relative tolerance
 * (elec_integrate's two-pass scheme below), which changes the tabulated
 * values at the grid's low-T edge. v3: constants_hash narrowed to me, the one
 * constant this table reads. Must stay equal to plasma.py's
 * ELECTRON_THERMO_FORMAT_VERSION -- the two backends share one cache file. */
#define ELECTRON_THERMO_FORMAT_VERSION 3

/* Fixed (T_min [MeV], T_max [MeV], n_pts) of the QED plasma-pressure
 * correction tables. Unlike the electron-thermo grid, whose upper edge tracks
 * cfg->T_start_cosmo_MeV, this one is deliberately constant (see the
 * extrapolation-domain warning above for why it is not rescaled). Named here
 * because the values must be used twice and the two uses MUST agree: they are
 * both an input to cpr_qed_compute_tables and fields of the fingerprint
 * (cpr_qed_fingerprint) deciding whether a cached table describes the grid we
 * are about to build. Must stay equal to plasma.py's _QED_TABLE_GRID -- the
 * two backends share one cache file. */
#define QED_T_MIN 1e-3
#define QED_T_MAX 1e2
#define QED_N_PTS ((size_t)500)

double cpr_rho_g(double Tg)    { return 2.0 * (M_PI * M_PI / 30.0) * pow(Tg, 4.0); }
double cpr_drho_g_dT(double Tg) { return 4.0 * cpr_rho_g(Tg) / Tg; }
double cpr_rho_nu(double Tnu)  { return 2.0 * (7.0 / 8.0) * (M_PI * M_PI / 30.0) * pow(Tnu, 4.0); }
double cpr_drho_nu_dT(double Tnu) { return 4.0 * cpr_rho_nu(Tnu) / Tnu; }

/* Per-flavour (nu+nubar) energy-density excess from a genuine reduced chemical
 * potential c = mu/Tnu: rho(c) - rho(0) = Tnu^4 (c^2/4 + c^4/(8 pi^2)). Even in
 * c (the antineutrino carries -c). See cpr_rho_nu_chempot_excess in plasma.h. */
double cpr_rho_nu_chempot_excess(double Tnu, double c)
{
    return pow(Tnu, 4.0) * (c * c / 4.0 + pow(c, 4.0) / (8.0 * M_PI * M_PI));
}

/* ---------------------------------------------------------------------
 * CPRInterp1D: not-a-knot cubic spline, however the values were obtained.
 * The !is_spline (raw x/y linear) branch is retained only as the transient
 * state in which qed_load_tables stages a freshly read file before fitting.
 * ------------------------------------------------------------------- */

double cpr_interp1d_eval(const CPRInterp1D *itp, double xq)
{
    if (itp->is_spline) return cpr_cubic_spline_eval(&itp->spl, xq);
    return cpr_interp_linear(itp->x, itp->y, itp->n, xq, CPR_EXTRAP_LINEAR);
}

void cpr_interp1d_free(CPRInterp1D *itp)
{
    if (itp->is_spline) cpr_cubic_spline_free(&itp->spl);
    else { free(itp->x); free(itp->y); }
    memset(itp, 0, sizeof(*itp));
}

static int file_exists(const char *path)
{
    struct stat st;
    return stat(path, &st) == 0;
}

/* ---------------------------------------------------------------------
 * QED interaction-pressure correction tables (plasma.Plasma._load_tables).
 * ------------------------------------------------------------------- */

static int load_qed_tables(CPRPlasma *pl, const CPRConfig *cfg, char **errmsg)
{
    if (!cfg->QED_corrections) {
        pl->qed_active = 0;
        memset(&pl->P_QED, 0, sizeof(pl->P_QED));
        memset(&pl->dP_QED, 0, sizeof(pl->dP_QED));
        memset(&pl->d2P_QED, 0, sizeof(pl->d2P_QED));
        return 0;
    }
    pl->qed_active = 1;

    /* Both the shipped files and the analytic fallback (below) cover a
     * fixed T in [1e-3, 1e2] MeV. Unlike the electron-thermo cache's Tmax
     * (build_electron_tables, scaled with cfg->T_start_cosmo_MeV), this
     * upper bound is NOT rescaled: delta_P_a/delta_P_e3 grow ~T^4 (times a
     * slowly-varying log/alpha prefactor, Phys. Rep. Eq. 47-49), so the
     * linear extrapolation used below/past 100 MeV is only trustworthy
     * close to the boundary, not as a stand-in for the true T^4-ish growth.
     * Mirrors the Python-side warning in plasma.py's _setup_qed_pressure. */
    if (cfg->T_start_cosmo_MeV > 100.0) {
        fprintf(stderr,
                "warning: T_start_cosmo_MeV=%.6g MeV exceeds the QED plasma-pressure "
                "correction table's fixed upper bound (100 MeV); dP_QED will be "
                "linearly extrapolated above 100 MeV, underestimating its true "
                "(~T^4) growth. Neff/YP results above this temperature may be biased.\n",
                cfg->T_start_cosmo_MeV);
    }

    char e2_file[CPR_PATH_BUF_LEN2], e3_file[CPR_PATH_BUF_LEN2], old_file[CPR_PATH_BUF_LEN2];
    /* Legacy 3-file names for backward compat with old cached copies. */
    char p_file_leg[CPR_PATH_BUF_LEN2], dp_file_leg[CPR_PATH_BUF_LEN2], d2p_file_leg[CPR_PATH_BUF_LEN2];
    /* Overlay reads: each QED table is resolved individually through the
     * cache_dir->shipped overlay; any recompute is WRITTEN to the writable
     * base's plasma/ subdir (plasma_wdir). Mirrors plasma.py. */
    char plasma_wdir[CPR_PATH_BUF_LEN];
    cpr_config_cache_write_dir(cfg, "plasma", plasma_wdir, sizeof(plasma_wdir));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_pressure_correction_e2.txt", e2_file, sizeof(e2_file));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_pressure_correction_e3.txt", e3_file, sizeof(e3_file));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_tables.txt", old_file, sizeof(old_file));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_P_int.txt", p_file_leg, sizeof(p_file_leg));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_dP_intdT.txt", dp_file_leg, sizeof(dp_file_leg));
    cpr_config_resolve_cache_file(cfg, "plasma", "QED_d2P_intdT2.txt", d2p_file_leg, sizeof(d2p_file_leg));

    /* --- Fingerprint gate on the current two-file format ------------------
     * The tables are a function of alpha and me (through the integrands) and
     * of the grid; cpr_qed_fingerprint records both. QED_T_MIN/QED_T_MAX/QED_N_PTS
     * are shared with the cpr_qed_compute_tables call below so the hash always
     * describes the grid actually built. Mirrors plasma.py's _load_tables. */
    CPRFPField qed_fields[5];
    size_t n_qed_fp = cpr_qed_fingerprint(QED_T_MIN, QED_T_MAX, QED_N_PTS,
                                          cfg->consts_hash[CPR_CONSTS_QED], qed_fields);
    char *qed_fp_hash = cpr_fingerprint_hash(qed_fields, n_qed_fp);

    int split_on_disk = file_exists(e2_file) && file_exists(e3_file);
    int split_valid = 0;
    if (split_on_disk) {
        char *h2 = cpr_cache_read_fingerprint_hash(e2_file);
        char *h3 = cpr_cache_read_fingerprint_hash(e3_file);
        split_valid = (h2 && h3 && strcmp(h2, qed_fp_hash) == 0
                              && strcmp(h3, qed_fp_hash) == 0);
        free(h2);
        free(h3);
    }
    /* A current-format pair whose header does not match (different constants,
     * different grid, or no header at all) is STALE: rebuild rather than use
     * it -- and rather than falling back to one of the superseded layouts
     * below, which would be staler still. It is NOT overwritten: these files
     * keep fixed names, so one config's tables would replace another's, and
     * alphaem/me are ordinary parameters. Mirrors plasma.py. */
    int split_stale = split_on_disk && !split_valid;

    /* The superseded layouts (single 7-column QED_tables.txt, and the older
     * 3-file trio) predate fingerprinting and can never carry a matching
     * header. They stay readable so an old cached checkout still starts, but
     * only when no current-format pair exists at all. */
    int old_present    = file_exists(old_file);
    int legacy_present = file_exists(p_file_leg) && file_exists(dp_file_leg)
                         && file_exists(d2p_file_leg);
    int split_present  = split_valid;
    int files_present  = split_valid
                         || (!split_stale && (old_present || legacy_present));
    int recompute = cfg->recompute_qed_corrections;
    free(qed_fp_hash);

    if (recompute || !files_present) {
        /* Analytic path: compute on a fresh 500-point grid (~0.3 s) and
         * build not-a-knot cubic-spline interpolants directly from the
         * computed arrays -- smoother than the linear interpolation used
         * when loading from a file (mirrors Python's choice exactly). */
        /* Written only when the user asked for a recompute: neither a miss
         * nor a stale pair saves, because any write to these fixed names
         * replaces whatever configuration's tables are already there.
         * Redirect cache_dir to keep a second configuration's pair.
         * Mirrors plasma.py. */
        int save = recompute;
        cpr_log(cfg, "init", "Computing QED plasma-pressure tables (%s)...",
                 recompute ? "recompute requested"
                           : split_stale ? "cached tables stale (fingerprint mismatch)"
                                         : "files not found");
        CPRQEDTables t;
        if (cpr_qed_compute_tables(QED_T_MIN, QED_T_MAX, QED_N_PTS,
                                    cfg->consts.alphaem, cfg->consts.me, &t, errmsg))
            return 1;
        if (save) {
            /* Non-fatal on a read-only install: the freshly computed
             * tables below are valid, only the disk cache is skipped -- warn
             * and point the user at the cache_dir remedy, do NOT abort. */
            if (cpr_qed_save_tables(&t, plasma_wdir,
                                     QED_T_MIN, QED_T_MAX, QED_N_PTS,
                                     cfg->consts_hash[CPR_CONSTS_QED], errmsg)) {
                cpr_log(cfg, "plasma",
                        "could not write cache to %s: results are unaffected, "
                        "but the next run will recompute. Set the cache_dir "
                        "parameter to redirect the cache to a writable directory.",
                        plasma_wdir);
                if (errmsg && *errmsg) { free(*errmsg); *errmsg = NULL; }
            }
        }
        double *sumP = CPR_XMALLOC(t.n * sizeof(double));
        double *sumdP = CPR_XMALLOC(t.n * sizeof(double));
        double *sumd2P = CPR_XMALLOC(t.n * sizeof(double));
        for (size_t i = 0; i < t.n; i++) {
            sumP[i]   = t.dP_e2[i] + t.dP_e3[i];
            sumdP[i]  = t.d_dP_e2_dT[i] + t.d_dP_e3_dT[i];
            sumd2P[i] = t.d2_dP_e2_dT2[i] + t.d2_dP_e3_dT2[i];
        }
        pl->P_QED.is_spline = pl->dP_QED.is_spline = pl->d2P_QED.is_spline = 1;
        int rc = cpr_cubic_spline_fit_notaknot(t.T, sumP, t.n, &pl->P_QED.spl, errmsg)
              || cpr_cubic_spline_fit_notaknot(t.T, sumdP, t.n, &pl->dP_QED.spl, errmsg)
              || cpr_cubic_spline_fit_notaknot(t.T, sumd2P, t.n, &pl->d2P_QED.spl, errmsg);
        free(sumP); free(sumdP); free(sumd2P);
        cpr_qed_tables_free(&t);
        return rc;
    }

    /* File mode: load and sum the e^2/e^3 columns into (x, y) pairs, which the
     * loop at the end of this function then converts to the SAME not-a-knot
     * cubic spline the analytic path above builds (matches plasma.py's
     * _qed_spline). The columns used to be interpolated linearly here, which
     * left the "loaded" and "computed" paths ~8e-4 apart in delta_P -- enough
     * to move Neff in its 6th decimal depending only on whether the tables
     * happened to be on disk. */
    CPRTable tab;
    if (split_present) {
        /* Current format: two 4-column files, one per order in e
         * (T, dP, d(dP)/dT, d2(dP)/dT2), summed column-by-column. Both
         * files share the same T grid (generated together), so the e2
         * file's T column is used for both interpolants' x-axis. */
        CPRTable tab_e2, tab_e3;
        if (cpr_table_read(e2_file, 4, &tab_e2, errmsg)) return 1;
        if (cpr_table_read(e3_file, 4, &tab_e3, errmsg)) { cpr_table_free(&tab_e2); return 1; }
        CPRInterp1D *targets[3] = { &pl->P_QED, &pl->dP_QED, &pl->d2P_QED };
        for (int k = 0; k < 3; k++) {
            targets[k]->is_spline = 0;
            targets[k]->n = tab_e2.n_rows;
            targets[k]->x = CPR_XMALLOC(tab_e2.n_rows * sizeof(double));
            targets[k]->y = CPR_XMALLOC(tab_e2.n_rows * sizeof(double));
            for (size_t i = 0; i < tab_e2.n_rows; i++) {
                targets[k]->x[i] = tab_e2.cols[0][i];
                targets[k]->y[i] = tab_e2.cols[k + 1][i] + tab_e3.cols[k + 1][i];
            }
        }
        cpr_table_free(&tab_e2);
        cpr_table_free(&tab_e3);
    } else if (old_present) {
        /* Older 7-column format: T, dP_a, dP_e3, d(dP_a)/dT, d(dP_e3)/dT,
         * d2(dP_a)/dT2, d2(dP_e3)/dT2. */
        if (cpr_table_read(old_file, 7, &tab, errmsg)) return 1;
        /* col indices: 0=T, 1=dP_a, 2=dP_e3, 3=ddP_a/dT, 4=ddP_e3/dT,
         *              5=d2dP_a/dT2, 6=d2dP_e3/dT2 */
        int col_pairs[3][2] = { {1,2}, {3,4}, {5,6} };
        CPRInterp1D *targets[3] = { &pl->P_QED, &pl->dP_QED, &pl->d2P_QED };
        for (int k = 0; k < 3; k++) {
            targets[k]->is_spline = 0;
            targets[k]->n = tab.n_rows;
            targets[k]->x = CPR_XMALLOC(tab.n_rows * sizeof(double));
            targets[k]->y = CPR_XMALLOC(tab.n_rows * sizeof(double));
            for (size_t i = 0; i < tab.n_rows; i++) {
                targets[k]->x[i] = tab.cols[0][i];
                targets[k]->y[i] = tab.cols[col_pairs[k][0]][i]
                                  + tab.cols[col_pairs[k][1]][i];
            }
        }
        cpr_table_free(&tab);
    } else {
        /* Legacy 3-file format: backward compat with old cached copies. */
        const char *files[3]    = { p_file_leg, dp_file_leg, d2p_file_leg };
        CPRInterp1D *targets[3] = { &pl->P_QED, &pl->dP_QED, &pl->d2P_QED };
        for (int k = 0; k < 3; k++) {
            if (cpr_table_read(files[k], 3, &tab, errmsg)) return 1;
            targets[k]->is_spline = 0;
            targets[k]->n = tab.n_rows;
            targets[k]->x = CPR_XMALLOC(tab.n_rows * sizeof(double));
            targets[k]->y = CPR_XMALLOC(tab.n_rows * sizeof(double));
            for (size_t i = 0; i < tab.n_rows; i++) {
                targets[k]->x[i] = tab.cols[0][i];
                targets[k]->y[i] = tab.cols[1][i] + tab.cols[2][i];
            }
            cpr_table_free(&tab);
        }
    }

    /* Convert the three freshly loaded (x, y) pairs into not-a-knot cubic
     * splines, so every file format above lands on the same interpolant as
     * the analytic branch (see the file-mode comment). Done once here rather
     * than in each branch because all three formats fill x/y identically. */
    CPRInterp1D *loaded[3] = { &pl->P_QED, &pl->dP_QED, &pl->d2P_QED };
    for (int k = 0; k < 3; k++) {
        CPRCubicSpline spl;
        if (cpr_cubic_spline_fit_notaknot(loaded[k]->x, loaded[k]->y,
                                          loaded[k]->n, &spl, errmsg)) {
            /* Release every target before bailing: those with index < k are
             * already converted splines, those with index >= k still hold the
             * raw x/y pair allocated above. cpr_interp1d_free dispatches on
             * each one's own is_spline flag, so a single loop covers both
             * cases (it also zeroes them, so the caller's cleanup is a
             * no-op rather than a double free). */
            for (int j = 0; j < 3; j++) cpr_interp1d_free(loaded[j]);
            return 1;
        }
        free(loaded[k]->x);
        free(loaded[k]->y);
        loaded[k]->x = NULL;
        loaded[k]->y = NULL;
        loaded[k]->spl = spl;
        loaded[k]->is_spline = 1;
    }
    return 0;
}

/* ---------------------------------------------------------------------
 * e+- exact integrands and quadrature (plasma.Plasma._*_exact).
 *
 * Each integrand is evaluated over the dimensionless energy variable
 * E = eps/Tg, lower bound x = me/Tg, fixed upper bound 100 (well past
 * where exp(-E) makes any further contribution negligible at double
 * precision). As in qed_pressure.c, a single cpr_quad_adaptive call over
 * the full [x, 100] risks missing the E~O(1-5) peak when x is small and
 * the domain is wide (the coarse first-level Simpson sample would then
 * land entirely in the exponentially-suppressed tail) -- so integration
 * is split into breakpoints anchored at x and widening geometrically,
 * exactly the same fix as cpr_qed_I01/I2m1.
 * ------------------------------------------------------------------- */

typedef double (*ElecIntegrand)(double E, double x);

static double rho_e_intgd(double E, double x)     { 
    if (E <= x) return 0.0; 
    return E * E * sqrt(E * E - x * x) / (exp(E) + 1.0); 
}
static double drho_e_dT_intgd(double E, double x) { 
    if (E <= x) return 0.0; 
    return E * E * E * sqrt(E * E - x * x) / pow(cosh(E / 2.0), 2.0); 
}
static double p_e_intgd(double E, double x)       { 
    if (E <= x) return 0.0; 
    return pow(E * E - x * x, 1.5) / (exp(E) + 1.0); 
}
static double dp_e_dT_intgd(double E, double x)   { 
    if (E <= x) return 0.0; 
    return E * pow(E * E - x * x, 1.5) / pow(cosh(E / 2.0), 2.0); 
}

typedef struct { ElecIntegrand fn; double x; } ElecCtx;

static double elec_quad_wrapper(double E, void *ctx)
{
    ElecCtx *c = (ElecCtx *)ctx;
    return c->fn(E, c->x);
}

static const double ELEC_OFFSETS[] = { 0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0 };
#define N_ELEC_OFFSETS (sizeof(ELEC_OFFSETS) / sizeof(ELEC_OFFSETS[0]))

/* One segmented sweep at a given ABSOLUTE per-segment tolerance. */
static double elec_integrate_at_tol(ElecCtx *ctx, double x, double upper, double seg_tol)
{
    double total = 0.0;
    double prev = x;
    for (size_t i = 1; i < N_ELEC_OFFSETS; i++) {
        double next = x + ELEC_OFFSETS[i];
        if (next >= upper) { next = upper; }
        if (next > prev)
            total += cpr_quad_adaptive(elec_quad_wrapper, ctx, prev, next, seg_tol, 30, NULL);
        prev = next;
        if (prev >= upper) break;
    }
    if (prev < upper)
        total += cpr_quad_adaptive(elec_quad_wrapper, ctx, prev, upper, seg_tol, 30, NULL);
    return total;
}

/* Relative accuracy target, matching plasma.py's quad(epsabs=0, epsrel=1e-12). */
#define ELEC_QUAD_EPSREL 1e-12

/* Two-pass so the tolerance is RELATIVE to the integral's own magnitude.
 *
 * cpr_quad_adaptive's `tol` is absolute (it is compared against
 * |refined - whole|), so this used to pass a fixed 1e-12/N_SEGMENTS. At the
 * electron-thermo grid's low-T edge (x = me/Tgamma = 30) the integrals are of
 * order e^-30 ~ 9e-14 -- smaller than that tolerance -- so every segment was
 * accepted at its first Simpson estimate and the result had no significant
 * digits. Python had the identical defect (epsabs=1e-12), which is why the two
 * backends' tables disagreed by ~1e-4 on rho_e/p_e down there despite sharing
 * one fingerprint.
 *
 * Pass 1 gets the magnitude at a cheap loose tolerance; pass 2 redoes the
 * sweep at epsrel * |I_est|, floored at DBL_MIN-ish so a genuinely zero
 * integral cannot ask for an unreachable absolute accuracy. See
 * ELECTRON_THERMO_FORMAT_VERSION (bumped to 2 in both backends). */
static double elec_integrate(ElecIntegrand fn, double x, double upper)
{
    ElecCtx ctx = { fn, x };
    double est = elec_integrate_at_tol(&ctx, x, upper, 1e-9 / (double)N_ELEC_OFFSETS);
    double scale = fabs(est);
    if (!(scale > 0.0)) return est;   /* exactly 0 (or NaN): nothing to refine */
    double seg_tol = ELEC_QUAD_EPSREL * scale / (double)N_ELEC_OFFSETS;
    if (seg_tol < 1e-300) seg_tol = 1e-300;
    return elec_integrate_at_tol(&ctx, x, upper, seg_tol);
}

static double rho_e_exact(double Tg, double me)
{
    if (Tg < me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    double r = elec_integrate(rho_e_intgd, me / Tg, 100.0);
    return 4.0 / (2.0 * M_PI * M_PI) * pow(Tg, 4.0) * r;
}

static double drho_e_dT_exact(double Tg, double me)
{
    if (Tg < me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    double r = elec_integrate(drho_e_dT_intgd, me / Tg, 100.0);
    return 1.0 / (2.0 * M_PI * M_PI) * pow(Tg, 3.0) * r;
}

static double p_e_exact(double Tg, double me)
{
    if (Tg < me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    double r = elec_integrate(p_e_intgd, me / Tg, 100.0);
    return 4.0 / (6.0 * M_PI * M_PI) * pow(Tg, 4.0) * r;
}

static double dp_e_dT_exact(double Tg, double me)
{
    if (Tg < me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    double r = elec_integrate(dp_e_dT_intgd, me / Tg, 100.0);
    return 1.0 / (6.0 * M_PI * M_PI) * pow(Tg, 3.0) * r;
}

/* ---------------------------------------------------------------------
 * e+- pre-tabulation with a fingerprinted on-disk cache
 * (plasma.Plasma._build_electron_tables).
 * ------------------------------------------------------------------- */

static int build_electron_tables(CPRPlasma *pl, const CPRConfig *cfg, char **errmsg)
{
    const double me = cfg->consts.me;
    double Tmin = me / ELEC_THERMO_LOWT_RATIO;
    double Tmax = fmax(cfg->T_start_cosmo_MeV, 100.0) * 1.5;
    size_t npts = (size_t)cfg->n_electron_table;

    CPRFPField fields[4];
    fields[0] = (CPRFPField){ "format_version", { CPR_INT, { .i = ELECTRON_THERMO_FORMAT_VERSION } } };
    fields[1] = (CPRFPField){ "n_electron_table", { CPR_INT, { .i = cfg->n_electron_table } } };
    fields[2] = (CPRFPField){ "T_start_cosmo_MeV", { CPR_DOUBLE, { .d = cfg->T_start_cosmo_MeV } } };
    /* Physical constants: the e+- integrands and the grid's lower edge
     * (Tmin = me/30 above) are functions of the electron mass, so editing `me`
     * changes every row. Mirrors _build_electron_tables in plasma.py. */
    fields[3] = (CPRFPField){ "constants_hash",
                              { CPR_STRING, { .s = cfg->consts_hash[CPR_CONSTS_ELEC_THERMO] } } };
    char *fp_hash = cpr_fingerprint_hash(fields, 4);

    /* The hash goes in the FILENAME, not just the header, mirroring the weak
     * tree's nTOp_<hash>.txt (and plasma.py's _build_electron_tables). With a
     * single fixed name two configurations could not coexist: each run whose
     * fingerprint differed from the file on disk overwrote it, leaving the
     * git-tracked shipped copy modified and making every alternation pay the
     * rebuild. There is deliberately no fallback to the old fixed
     * "electron_thermo_cache.txt" name -- a miss costs one table rebuild,
     * cheaper than carrying compatibility code. */
    char cache_name[64];
    snprintf(cache_name, sizeof(cache_name), "electron_thermo_%s.txt", fp_hash);

    /* Overlay read (cache_dir first, else shipped copy); write to the writable
     * base's plasma/ subdir (cache_dir if set, else the package tree). */
    char cache_read[CPR_PATH_BUF_LEN2];
    cpr_config_resolve_cache_file(cfg, "plasma", cache_name,
                                  cache_read, sizeof(cache_read));

    if (!cfg->recompute_electron_thermo) {
        char *cached_hash = cpr_cache_read_fingerprint_hash(cache_read);
        if (cached_hash && strcmp(cached_hash, fp_hash) == 0) {
            free(cached_hash);
            CPRTable tab;
            if (cpr_table_read(cache_read, 5, &tab, errmsg) == 0) {
                int rc = cpr_cubic_spline_fit_notaknot(tab.cols[0], tab.cols[1], tab.n_rows, &pl->rho_e_tab, errmsg)
                      || cpr_cubic_spline_fit_notaknot(tab.cols[0], tab.cols[2], tab.n_rows, &pl->p_e_tab, errmsg)
                      || cpr_cubic_spline_fit_notaknot(tab.cols[0], tab.cols[3], tab.n_rows, &pl->drho_e_dT_tab, errmsg)
                      || cpr_cubic_spline_fit_notaknot(tab.cols[0], tab.cols[4], tab.n_rows, &pl->dp_e_dT_tab, errmsg);
                cpr_table_free(&tab);
                free(fp_hash);
                if (rc == 0) {
                    cpr_log(cfg, "init", "Electron-thermo tables loaded from cache (%d points).",
                             cfg->n_electron_table);
                } else {
                    /* The || chain short-circuits, so an unknown prefix of the
                     * four fits succeeded and holds allocations. Free all four:
                     * pl was memset to 0 by cpr_plasma_init, so the not-yet-
                     * fitted ones are all-NULL and cpr_cubic_spline_free is a
                     * no-op on them. */
                    cpr_cubic_spline_free(&pl->rho_e_tab);
                    cpr_cubic_spline_free(&pl->p_e_tab);
                    cpr_cubic_spline_free(&pl->drho_e_dT_tab);
                    cpr_cubic_spline_free(&pl->dp_e_dT_tab);
                }
                return rc;
            }
            /* Fall through to recompute if the cache file turned out to
             * be unreadable despite a matching fingerprint header
             * (matches Python's try/except warn-and-recompute path). */
        }
        free(cached_hash);
    }

    double *grid = CPR_XMALLOC(npts * sizeof(double));
    double *rho_e_arr = CPR_XMALLOC(npts * sizeof(double));
    double *p_e_arr = CPR_XMALLOC(npts * sizeof(double));
    double *drho_e_dT_arr = CPR_XMALLOC(npts * sizeof(double));
    double *dp_e_dT_arr = CPR_XMALLOC(npts * sizeof(double));

    double log_min = log10(Tmin), log_max = log10(Tmax);
    for (size_t i = 0; i < npts; i++) {
        double frac = (npts == 1) ? 0.0 : (double)i / (double)(npts - 1);
        grid[i] = pow(10.0, log_min + frac * (log_max - log_min));
        rho_e_arr[i] = rho_e_exact(grid[i], me);
        p_e_arr[i] = p_e_exact(grid[i], me);
        drho_e_dT_arr[i] = drho_e_dT_exact(grid[i], me);
        dp_e_dT_arr[i] = dp_e_dT_exact(grid[i], me);
    }

    int rc = cpr_cubic_spline_fit_notaknot(grid, rho_e_arr, npts, &pl->rho_e_tab, errmsg)
          || cpr_cubic_spline_fit_notaknot(grid, p_e_arr, npts, &pl->p_e_tab, errmsg)
          || cpr_cubic_spline_fit_notaknot(grid, drho_e_dT_arr, npts, &pl->drho_e_dT_tab, errmsg)
          || cpr_cubic_spline_fit_notaknot(grid, dp_e_dT_arr, npts, &pl->dp_e_dT_tab, errmsg);

    if (rc == 0) {
        double *columns[5] = { grid, rho_e_arr, p_e_arr, drho_e_dT_arr, dp_e_dT_arr };
        /* Write to the writable base (cache_dir if set, else the package's
         * cache_plasma_weak/plasma/). A cache-write failure is non-fatal
         * (matches Python's warn-and-continue): the tables we just built in
         * memory are still valid for this run, only the on-disk cache for
         * future runs is skipped -- warn and name the cache_dir remedy. */
        char cache_wdir[CPR_PATH_BUF_LEN];
        cpr_config_cache_write_dir(cfg, "plasma", cache_wdir, sizeof(cache_wdir));
        char cache_write[CPR_PATH_BUF_LEN2];
        snprintf(cache_write, sizeof(cache_write), "%s/%s", cache_wdir, cache_name);
        /* Create the dir tree on demand (a fresh cache_dir has no plasma/). */
        char mkdir_cmd[CPR_PATH_BUF_LEN2];
        snprintf(mkdir_cmd, sizeof(mkdir_cmd), "%s/", cache_wdir);
        for (char *p = mkdir_cmd + 1; *p; p++) {
            if (*p == '/') { *p = '\0'; mkdir(mkdir_cmd, 0755); *p = '/'; }
        }
        if (cpr_cache_write(cache_write, fields, 4, "grid rho_e p_e drho_e_dT dp_e_dT",
                             columns, 5, npts, NULL) != 0) {
            cpr_log(cfg, "plasma",
                    "could not write cache to %s: results are unaffected, but "
                    "the next run will recompute. Set the cache_dir parameter "
                    "to redirect the cache to a writable directory.", cache_write);
        }
    }

    free(grid); free(rho_e_arr); free(p_e_arr); free(drho_e_dT_arr); free(dp_e_dT_arr);
    free(fp_hash);
    if (rc == 0)
        cpr_log(cfg, "init", "Electron-thermo tables built (%d points).", cfg->n_electron_table);
    return rc;
}

/* ---------------------------------------------------------------------
 * Public API.
 * ------------------------------------------------------------------- */

int cpr_plasma_init(CPRPlasma *pl, const CPRConfig *cfg, char **errmsg)
{
    memset(pl, 0, sizeof(*pl));
    pl->cfg = cfg;
    if (load_qed_tables(pl, cfg, errmsg)) return 1;
    if (build_electron_tables(pl, cfg, errmsg)) {
        cpr_interp1d_free(&pl->P_QED);
        cpr_interp1d_free(&pl->dP_QED);
        cpr_interp1d_free(&pl->d2P_QED);
        return 1;
    }
    cpr_log(cfg, "init", "QED pressure corrections tables loaded.");
    return 0;
}

void cpr_plasma_free(CPRPlasma *pl)
{
    if (pl->qed_active) {
        cpr_interp1d_free(&pl->P_QED);
        cpr_interp1d_free(&pl->dP_QED);
        cpr_interp1d_free(&pl->d2P_QED);
    }
    cpr_cubic_spline_free(&pl->rho_e_tab);
    cpr_cubic_spline_free(&pl->p_e_tab);
    cpr_cubic_spline_free(&pl->drho_e_dT_tab);
    cpr_cubic_spline_free(&pl->dp_e_dT_tab);
    memset(pl, 0, sizeof(*pl));
}

static double qed_P(const CPRPlasma *pl, double Tg)   { return pl->qed_active ? cpr_interp1d_eval(&pl->P_QED, Tg) : 0.0; }
static double qed_dP(const CPRPlasma *pl, double Tg)  { return pl->qed_active ? cpr_interp1d_eval(&pl->dP_QED, Tg) : 0.0; }
static double qed_d2P(const CPRPlasma *pl, double Tg) { return pl->qed_active ? cpr_interp1d_eval(&pl->d2P_QED, Tg) : 0.0; }

double cpr_plasma_rho_e(const CPRPlasma *pl, double Tg)
{
    if (Tg < pl->cfg->consts.me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    return cpr_cubic_spline_eval(&pl->rho_e_tab, Tg);
}

double cpr_plasma_drho_e_dT(const CPRPlasma *pl, double Tg)
{
    if (Tg < pl->cfg->consts.me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    return cpr_cubic_spline_eval(&pl->drho_e_dT_tab, Tg);
}

double cpr_plasma_p_e(const CPRPlasma *pl, double Tg)
{
    if (Tg < pl->cfg->consts.me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    return cpr_cubic_spline_eval(&pl->p_e_tab, Tg);
}

double cpr_plasma_dp_e_dT(const CPRPlasma *pl, double Tg)
{
    if (Tg < pl->cfg->consts.me / ELEC_THERMO_LOWT_RATIO) return 0.0;
    return cpr_cubic_spline_eval(&pl->dp_e_dT_tab, Tg);
}

double cpr_plasma_rho_nu_extra(const CPRPlasma *pl, double Tg)
{
    if (pl->cfg->DeltaNeff == 0.0) return 0.0;
    double Tnu_dec = cpr_plasma_T_nu_decoupling(pl, Tg);
    return pl->cfg->DeltaNeff * 2.0 * (7.0 / 8.0) * (M_PI * M_PI / 30.0) * pow(Tnu_dec, 4.0);
}

double cpr_plasma_rho_SM(const CPRPlasma *pl, double Tg, double Tnue, double Tnumu)
{
    double rho_qed = Tg * qed_dP(pl, Tg) - qed_P(pl, Tg);
    return cpr_rho_g(Tg) + cpr_plasma_rho_e(pl, Tg) + rho_qed
         + cpr_rho_nu(Tnue) + 2.0 * cpr_rho_nu(Tnumu)
         + cpr_plasma_rho_nu_extra(pl, Tg);
}

double cpr_plasma_p_SM(const CPRPlasma *pl, double Tg, double Tnue, double Tnumu)
{
    return cpr_rho_g(Tg) / 3.0 + cpr_plasma_p_e(pl, Tg) + qed_P(pl, Tg)
         + (cpr_rho_nu(Tnue) + 2.0 * cpr_rho_nu(Tnumu)) / 3.0
         + cpr_plasma_rho_nu_extra(pl, Tg) / 3.0;
}

double cpr_plasma_spl(const CPRPlasma *pl, double Tg)
{
    double rho_pl = cpr_rho_g(Tg) + cpr_plasma_rho_e(pl, Tg);
    double p_pl   = cpr_rho_g(Tg) / 3.0 + cpr_plasma_p_e(pl, Tg);
    double rho_qed = Tg * qed_dP(pl, Tg) - qed_P(pl, Tg);
    double p_qed   = qed_P(pl, Tg);
    return (rho_pl + p_pl + rho_qed + p_qed) / Tg;
}

void cpr_plasma_spl_and_dspl_dT(const CPRPlasma *pl, double Tg, double *s, double *ds_dT)
{
    double rho_g_val = cpr_rho_g(Tg);
    double rho_e_val = cpr_plasma_rho_e(pl, Tg);
    double p_e_val   = cpr_plasma_p_e(pl, Tg);
    double P_val   = qed_P(pl, Tg);
    double dP_val  = qed_dP(pl, Tg);
    double d2P_val = qed_d2P(pl, Tg);

    double rho_pl = rho_g_val + rho_e_val;
    double p_pl   = rho_g_val / 3.0 + p_e_val;
    double rho_qed = Tg * dP_val - P_val;
    double p_qed   = P_val;
    *s = (rho_pl + p_pl + rho_qed + p_qed) / Tg;

    double drho_g_val = cpr_drho_g_dT(Tg);
    double drho_pl_dT = drho_g_val + cpr_plasma_drho_e_dT(pl, Tg);
    double dp_pl_dT   = drho_g_val / 3.0 + cpr_plasma_dp_e_dT(pl, Tg);
    double drho_qed_dT = Tg * d2P_val;  /* d/dT[T dP/dT - P] = T d^2P/dT^2 */
    double dp_qed_dT   = dP_val;        /* d/dT[P] = dP/dT */
    *ds_dT = (drho_pl_dT + dp_pl_dT + drho_qed_dT + dp_qed_dT) / Tg - *s / Tg;
}

double cpr_plasma_dspl_dT(const CPRPlasma *pl, double Tg)
{
    double s, ds_dT;
    cpr_plasma_spl_and_dspl_dT(pl, Tg, &s, &ds_dT);
    return ds_dT;
}

double cpr_plasma_T_nu_decoupling(const CPRPlasma *pl, double Tg)
{
    return Tg * pow(cpr_plasma_spl(pl, Tg) / (SIGMA_INF * pow(Tg, 3.0)), 1.0 / 3.0);
}
