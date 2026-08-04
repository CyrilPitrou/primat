/* nuclear_network.c -- see nuclear_network.h.
 *
 * Reference: Pitrou, Coc, Uzan & Vangioni, Phys. Rep. 2018 (arXiv:1806.11095),
 * cited below as "Phys. Rep.".
 */
#include "nuclear_network.h"
#include "xalloc.h"
#include "constants.h"
#include "network_builder.h"
#include "ode_rk.h"
#include "ode_bdf.h"
#include "linalg.h"
#include "log.h"
#include "spline.h"   /* cpr_find_segment, for the per-reaction rate columns */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compat_posix.h"  /* sys/stat.h + mkdir, portable across POSIX & MSVC */
#include <time.h>

/* Riemann zeta(3) (Apery's constant) -- see constants.c's identical
 * literal; duplicated here (not exposed via constants.h) since it is
 * needed only by the Saha (YA) equilibrium formula below. */
#define ZETA3 1.2020569031595942854

static const CPRNuclide *find_nuclide(const CPRConfig *cfg, const char *name)
{
    for (size_t i = 0; i < cfg->nuclides.n; i++)
        if (strcmp(cfg->nuclides.items[i].name, name) == 0) return &cfg->nuclides.items[i];
    return NULL;
}

/* Saha (Nuclear Statistical Equilibrium) mass-fraction abundance of
 * nuclide `name`, in equilibrium with free neutrons/protons at
 * temperature T_K [Kelvin] and baryon-to-photon ratio eta_b (port of
 * solve()'s local YA closure). Phys. Rep. SS V.A:
 *
 *   Y_A = g_A zeta(3)^(A-1) pi^((1-A)/2) 2^((3A-5)/2)
 *         x (M_A / mn^N mp^Z)^(3/2)
 *         x (kB T)^(3(A-1)/2) eta_b^(A-1)
 *         x Yn^N Yp^Z exp(B_A / kB T)
 *
 * where A=N+Z, g_A=2J+1 (spin degeneracy), B_A the binding energy, all
 * masses/energies carried in erg (the "natural units" convention shared
 * with every other module here -- g_const.MeV/keV are the MeV/keV ->
 * erg conversion factors). Used only at T = T_weak (the MT-era seed,
 * cpr_nuclear_network_solve), so eta_b is the single value eta_b_weak
 * evaluated once there -- no eta_b(T) interpolant is needed. */
static double saha_YA(const CPRConfig *cfg, double eta_b, const char *name,
                       double Yn, double Yp, double T_K)
{
    const CPRNuclide *nuc = find_nuclide(cfg, name);
    const CPRNuclide *nuc_n = find_nuclide(cfg, "n");
    const CPRNuclide *nuc_p = find_nuclide(cfg, "p");
    double A = (double)(nuc->N + nuc->Z);
    double Z = (double)nuc->Z;
    double N = A - Z;

    double Mass = A * g_const.ma * g_const.MeV
                  + g_const.keV * nuc->mass_excess_keV
                  - Z * g_const.me * g_const.MeV;
    double BindE = N * nuc_n->mass_excess_keV + Z * nuc_p->mass_excess_keV
                   - nuc->mass_excess_keV;
    /* (M_A / mn^N mp^Z)^(3/2): ratio of nuclear to free-nucleon masses. */
    double NormYA = pow(Mass / (pow(g_const.mn * g_const.MeV, A - Z)
                                  * pow(g_const.mp * g_const.MeV, Z)),
                          1.5);

    return (2.0 * nuc->spin + 1.0)
           * pow(ZETA3, A - 1.0) * pow(M_PI, (1.0 - A) / 2.0)
           * pow(2.0, (3.0 * A - 5.0) / 2.0)
           * NormYA
           * pow(g_const.kB * T_K, 1.5 * (A - 1.0))
           * pow(eta_b, A - 1.0)
           * pow(Yp, Z) * pow(Yn, N)
           * exp(BindE * g_const.keV / (g_const.kB * T_K));
}

/* ---- Growable per-era (t, Y) recorder, fed by the ODE integrators'
 * step_cb hook; seeded with the initial point (the integrators only
 * report *accepted steps after* t0, mirroring solve_ivp.t/.y which
 * include the initial condition as their first row). ---- */
typedef struct {
    double *t;
    double *Y;   /* row-major (cap x n_sp) */
    size_t n_sp, n, cap;
} CPRRecorder;

static void recorder_init(CPRRecorder *r, size_t n_sp)
{
    r->n_sp = n_sp; r->n = 0; r->cap = 64;
    r->t = CPR_XMALLOC(r->cap * sizeof(double));
    r->Y = CPR_XMALLOC(r->cap * n_sp * sizeof(double));
}

static void recorder_push(CPRRecorder *r, double t, const double *y)
{
    if (r->n == r->cap) {
        r->cap *= 2;
        r->t = CPR_XREALLOC(r->t, r->cap * sizeof(double));
        r->Y = CPR_XREALLOC(r->Y, r->cap * r->n_sp * sizeof(double));
    }
    r->t[r->n] = t;
    memcpy(&r->Y[r->n * r->n_sp], y, r->n_sp * sizeof(double));
    r->n++;
}

static void recorder_cb(double t, const double *y, size_t n, void *ctx)
{
    (void)n;
    CPRRecorder *r = ctx;
    if (getenv("CPR_NN_DEBUG") && (r->n % 2000 == 0))
        fprintf(stderr, "[nn debug] step=%zu t=%.6e\n", r->n, t);
    recorder_push(r, t, y);
}

static void recorder_free(CPRRecorder *r)
{
    free(r->t); free(r->Y);
}

/* Scatters one era's local abundance row (in_row, named by in_names) into
 * the wider canonical-column row out_row (named by nn->abundance_names);
 * columns absent from in_names are left at out_row's existing value
 * (caller pre-zeros each row, mirroring _embed's np.zeros base). O(n_in *
 * n_out) name matching is fine here: both are at most ~60 (the `large`
 * network). */
static void embed_row(double *out_row, char (*out_names)[16], size_t n_out,
                       const double *in_row, char (*in_names)[16], size_t n_in)
{
    for (size_t j = 0; j < n_in; j++)
        for (size_t k = 0; k < n_out; k++)
            if (strcmp(in_names[j], out_names[k]) == 0) { out_row[k] = in_row[j]; break; }
}

/* ---- ODE right-hand-side / Jacobian glue: each era's CPRODEFunc/
 * CPRODEJacFunc closes over the background + (for MT/LT) compiled rate
 * kernels it needs, matching solve()'s local Y_prime_HT/MT/LT closures. */

typedef struct { CPRBackground *bg; } HTCtx;

static int ht_rhs(double t, const double *Y, double *dY, void *ctx)
{
    HTCtx *c = ctx;
    double T_K = cpr_bg_T_of_t(c->bg, t) * cpr_MeV_to_Kelvin();
    double f = cpr_bg_weak_nTOp_frwrd(c->bg, T_K);
    double b = cpr_bg_weak_nTOp_bkwrd(c->bg, T_K);
    dY[0] = b * Y[1] - f * Y[0];
    dY[1] = f * Y[0] - b * Y[1];
    return 0;
}

typedef struct { CPRBackground *bg; CPRNuclearRates *nucl; } MTLTCtx;

static int mt_rhs(double t, const double *Y, double *dY, void *ctx)
{
    MTLTCtx *c = ctx;
    double rho = cpr_bg_rhoB_BBN(c->bg, t);
    double T_K = cpr_bg_T_of_t(c->bg, t) * cpr_MeV_to_Kelvin();
    double f = cpr_bg_weak_nTOp_frwrd(c->bg, T_K), b = cpr_bg_weak_nTOp_bkwrd(c->bg, T_K);
    cpr_nuclear_rates_rhs_mt(c->nucl, Y, T_K, rho, f, b, dY);
    return 0;
}

static int mt_jac(double t, const double *Y, double *J, void *ctx)
{
    MTLTCtx *c = ctx;
    double rho = cpr_bg_rhoB_BBN(c->bg, t);
    double T_K = cpr_bg_T_of_t(c->bg, t) * cpr_MeV_to_Kelvin();
    double f = cpr_bg_weak_nTOp_frwrd(c->bg, T_K), b = cpr_bg_weak_nTOp_bkwrd(c->bg, T_K);
    cpr_nuclear_rates_jac_mt(c->nucl, Y, T_K, rho, f, b, J);
    return 0;
}

static int lt_rhs(double t, const double *Y, double *dY, void *ctx)
{
    MTLTCtx *c = ctx;
    double rho = cpr_bg_rhoB_BBN(c->bg, t);
    double T_K = cpr_bg_T_of_t(c->bg, t) * cpr_MeV_to_Kelvin();
    double f = cpr_bg_weak_nTOp_frwrd(c->bg, T_K), b = cpr_bg_weak_nTOp_bkwrd(c->bg, T_K);
    cpr_nuclear_rates_rhs_lt(c->nucl, Y, T_K, rho, f, b, dY);
    return 0;
}

static int lt_jac(double t, const double *Y, double *J, void *ctx)
{
    MTLTCtx *c = ctx;
    double rho = cpr_bg_rhoB_BBN(c->bg, t);
    double T_K = cpr_bg_T_of_t(c->bg, t) * cpr_MeV_to_Kelvin();
    double f = cpr_bg_weak_nTOp_frwrd(c->bg, T_K), b = cpr_bg_weak_nTOp_bkwrd(c->bg, T_K);
    cpr_nuclear_rates_jac_lt(c->nucl, Y, T_K, rho, f, b, J);
    return 0;
}

static double *find_in(double *raw_vals, char (*raw_names)[16], size_t n_raw, const char *name)
{
    for (size_t i = 0; i < n_raw; i++)
        if (strcmp(raw_names[i], name) == 0) return &raw_vals[i];
    return NULL;
}

int cpr_nuclear_network_solve(CPRNuclearNetwork *nn, const CPRConfig *cfg,
                                CPRNuclearRates *nucl, CPRBackground *background,
                                char **errmsg)
{
    memset(nn, 0, sizeof(*nn));
    nn->cfg = cfg; nn->background = background; nn->nucl = nucl;

    /* Refresh nuclear rates with the current rate-variation parameters
     * (mirrors solve()'s nucl.apply_variations(cfg) call at the top). */
    cpr_nuclear_rates_apply_variations(nucl, cfg);

    /* ---- Temperature era boundaries [s]. cpr_T_start/T_weak/T_nucl are
     * *fixed* era boundaries in Kelvin (10/1/0.11 MeV respectively,
     * independent of cfg -- see constants.h), unlike T_end which is the
     * user-configurable cfg->T_end_MeV. ---- */
    double T_start_K = cpr_T_start(), T_weak_K = cpr_T_weak(), T_nucl_K = cpr_T_nucl();
    double T_end_K = cpr_config_T_end(cfg);
    /* MeV values of era boundaries, used only in verbose log messages. */
    double T_start_MeV = T_start_K / cpr_MeV_to_Kelvin();
    double T_weak_MeV  = T_weak_K  / cpr_MeV_to_Kelvin();
    double T_nucl_MeV  = T_nucl_K  / cpr_MeV_to_Kelvin();
    double t_start = cpr_bg_t_of_T(background, T_start_K / cpr_MeV_to_Kelvin());
    double t_weak  = cpr_bg_t_of_T(background, T_weak_K  / cpr_MeV_to_Kelvin());
    double t_nucl  = cpr_bg_t_of_T(background, T_nucl_K  / cpr_MeV_to_Kelvin());
    double t_end   = cpr_bg_t_of_T(background, T_end_K   / cpr_MeV_to_Kelvin());
    nn->t_end = t_end;

    /* ---- Baryon-to-photon ratio at T_weak, for the MT-era Saha seed. ---- */
    double nB_weak = cpr_bg_rhoB_BBN(background, t_weak) / (g_const.ma * cpr_MeV4_to_gcmm3());
    double ngamma_weak = (2.0 * ZETA3 / (M_PI * M_PI)) * pow(T_weak_K / cpr_MeV_to_Kelvin(), 3.0);
    double eta_b_weak = nB_weak / ngamma_weak;

    /* ------------------------------------------------------------------
     * HT era: n <-> p only, non-stiff RK45.
     * ------------------------------------------------------------------ */
    double f0 = cpr_bg_weak_nTOp_frwrd(background, T_start_K);
    double b0 = cpr_bg_weak_nTOp_bkwrd(background, T_start_K);
    double Y_ht[2] = { b0 / (b0 + f0), 0.0 };
    Y_ht[1] = 1.0 - Y_ht[0];

    CPRRecorder rec_ht; recorder_init(&rec_ht, 2);
    recorder_push(&rec_ht, t_start, Y_ht);
    HTCtx ht_ctx = { background };
    /* HT integrator: Dormand-Prince RK45 here, LSODA on the Python side
     * (nuclear_network.py's _solve_HT). This is a KNOWN, accepted divergence,
     * not an oversight -- recorded here so it is not "fixed" by accident.
     *
     * Both run at the same rtol (cfg->numerical_precision) and atol (1e-10),
     * and the era is n <-> p only. Measured on the default config: swapping
     * Python's HT to RK45 moves YPBBN by ~5e-07 relative, i.e. the method
     * choice accounts for essentially all of the cross-backend YP gap. But
     * neither method is more accurate -- sweeping numerical_precision 1e-6 to
     * 1e-10 with the HT era set to LSODA, RK45 and BDF in turn shows all three
     * converging to the same YPBBN (spread 1.6e-07 relative at rtol=1e-9),
     * which is well inside the accuracy the default rtol=1e-7 delivers.
     *
     * Aligning them was tried (both on BDF with the exact analytic 2x2
     * Jacobian, the HT system being linear in Y) and did NOT improve
     * cross-backend agreement: YP parity got worse at the default tolerance
     * (5.2e-07 -> 1.1e-06) while D/H improved only ~1.5x, because the residual
     * is dominated by the MT/LT BDF solves walking different step sequences,
     * not by the HT method. Alignment therefore buys reviewability, not
     * numbers; it is available if wanted, but is not currently applied. */
    CPRRKOpts rk_opts = cpr_ode_rk_default_opts();
    rk_opts.rtol = cfg->numerical_precision; rk_opts.atol = 1.0e-10;
    if (cfg->show_progress && !cfg->verbose) {
        fprintf(stderr, "[primat]  HT."); fflush(stderr);
    }
    cpr_log(cfg, "nucl", "Solving neutron decoupling at high temperature era"
                         " (T = %.4g -> %.4g MeV)", T_start_MeV, T_weak_MeV);
    clock_t _t_ht0 = clock();
    if (cpr_ode_rk45(ht_rhs, &ht_ctx, t_start, t_weak, Y_ht, 2, rk_opts,
                      recorder_cb, &rec_ht, errmsg)) {
        recorder_free(&rec_ht);
        return 1;
    }
    cpr_log(cfg, "nucl", "[HT] Finished solve_ivp in %.2f s",
             (double)(clock() - _t_ht0) / CLOCKS_PER_SEC);
    double Yn_HT_f = Y_ht[0], Yp_HT_f = Y_ht[1];

    /* ------------------------------------------------------------------
     * MT era: fixed 18-reaction subset, stiff BDF with analytic Jacobian.
     * ------------------------------------------------------------------ */
    char (*mt_names)[16] = nucl->mt_net.species;
    size_t n_mt = nucl->mt_net.n_species;
    double *Yi_MT = CPR_XMALLOC(n_mt * sizeof(double));
    for (size_t i = 0; i < n_mt; i++) {
        if (strcmp(mt_names[i], "n") == 0) Yi_MT[i] = Yn_HT_f;
        else if (strcmp(mt_names[i], "p") == 0) Yi_MT[i] = Yp_HT_f;
        else Yi_MT[i] = saha_YA(cfg, eta_b_weak, mt_names[i], Yn_HT_f, Yp_HT_f, T_weak_K);
    }

    CPRRecorder rec_mt; recorder_init(&rec_mt, n_mt);
    recorder_push(&rec_mt, t_weak, Yi_MT);
    MTLTCtx mt_ctx = { background, nucl };
    CPRBDFOpts bdf_opts = cpr_ode_bdf_default_opts();
    /* atol 1e-15, matching primat/nuclear_network.py's MT solve_ivp call. It
     * read 1e-16 here, an undocumented tolerance divergence: swapping the two
     * values on the Python side moves YPBBN 0.24699729 -> 0.24699702 and D/H
     * 2.4358985e-05 -> 2.4358951e-05, i.e. 1.4e-06 relative, ~470x the +-3e-9
     * same-backend D/H regression pin (the C's own BDF is far less sensitive
     * here, ~4e-08, which is why it went unnoticed). */
    bdf_opts.rtol = cfg->numerical_precision; bdf_opts.atol = 1.0e-15;
    if (cfg->show_progress && !cfg->verbose) {
        fprintf(stderr, "  MT."); fflush(stderr);
    }
    cpr_log(cfg, "nucl", "Solving nuclear network at mid temperature era"
                         " (T = %.4g -> %.4g MeV)", T_weak_MeV, T_nucl_MeV);
    clock_t _t_mt0 = clock();
    if (cpr_ode_bdf(mt_rhs, mt_jac, &mt_ctx, t_weak, t_nucl, Yi_MT, n_mt, bdf_opts,
                     recorder_cb, &rec_mt, errmsg)) {
        free(Yi_MT); recorder_free(&rec_ht); recorder_free(&rec_mt);
        return 1;
    }
    cpr_log(cfg, "nucl", "[MT] Finished solve_ivp (%s network, %zu species) in %.2f s",
             cfg->network, n_mt, (double)(clock() - _t_mt0) / CLOCKS_PER_SEC);

    /* ------------------------------------------------------------------
     * LT era: the chosen network (small/large, optionally amax-restricted),
     * stiff BDF with analytic Jacobian.
     * ------------------------------------------------------------------ */
    char (*lt_names)[16] = nucl->lt_net.species;
    size_t n_lt = nucl->lt_net.n_species;
    double *Yi_LT = CPR_XMALLOC(n_lt * sizeof(double));
    for (size_t i = 0; i < n_lt; i++) {
        double *v = find_in(Yi_MT, mt_names, n_mt, lt_names[i]);
        Yi_LT[i] = v ? *v : 0.0;
    }

    CPRRecorder rec_lt; recorder_init(&rec_lt, n_lt);
    recorder_push(&rec_lt, t_nucl, Yi_LT);
    MTLTCtx lt_ctx = { background, nucl };
    CPRBDFOpts bdf_opts_lt = cpr_ode_bdf_default_opts();
    bdf_opts_lt.rtol = 10.0 * cfg->numerical_precision;
    /* Universal LT absolute tolerance (cfg->atol_large_LT) for every network,
     * not just "large" -- was `cpr_config_is_large(cfg) ? atol_large_LT :
     * 1e-20`, keyed on the literal network name, which broke bit-for-bit
     * reproduction of a custom network run under a renamed user_nuclear_dir
     * overlay (is_large=False -> looser atol). One atol everywhere removes that
     * name dependence; it only tightens non-large networks. Keep in lockstep
     * with primat/nuclear_network.py's `atol = cfg.atol_large_LT`. */
    bdf_opts_lt.atol = cfg->atol_large_LT;
    if (cfg->show_progress && !cfg->verbose) {
        fprintf(stderr, "  LT."); fflush(stderr);
    }
    cpr_log(cfg, "nucl", "Solving nuclear network at low temperature era"
                         " (T = %.4g -> %.4g MeV)", T_nucl_MeV, cfg->T_end_MeV);
    clock_t _t_lt0 = clock();
    if (cpr_ode_bdf(lt_rhs, lt_jac, &lt_ctx, t_nucl, t_end, Yi_LT, n_lt, bdf_opts_lt,
                     recorder_cb, &rec_lt, errmsg)) {
        free(Yi_MT); free(Yi_LT);
        recorder_free(&rec_ht); recorder_free(&rec_mt); recorder_free(&rec_lt);
        return 1;
    }
    cpr_log(cfg, "nucl", "[LT] Finished solve_ivp (%s network, %zu species) in %.2f s",
             cfg->network, n_lt, (double)(clock() - _t_lt0) / CLOCKS_PER_SEC);
    if (cfg->show_progress && !cfg->verbose) {
        fprintf(stderr, "  done.\n"); fflush(stderr);
    }

    /* ---- Final abundances: the LT species list is the canonical name
     * list for any network (mirrors solve()'s self.abundance_names = species_L). ---- */
    nn->n_species = n_lt;
    nn->abundance_names = CPR_XMALLOC(n_lt * sizeof(*nn->abundance_names));
    memcpy(nn->abundance_names, lt_names, n_lt * sizeof(*nn->abundance_names));
    nn->Y_final = CPR_XMALLOC(n_lt * sizeof(double));
    memcpy(nn->Y_final, Yi_LT, n_lt * sizeof(double));

    /* Mirrors nuclear_network.py's solve() verbose dump of the final
     * abundances (same header/format, no cpr_log tag prefix -- that block
     * isn't tagged on the Python side either). */
    if (cfg->verbose) {
        printf("--------------------------------------------------\n");
        printf("Primordial abundances (%zu nuclides) at T = %.4g MeV\n",
               n_lt, cfg->T_end_MeV);
        printf("--------------------------------------------------\n");
        for (size_t i = 0; i < n_lt; i++)
            printf("  Y%-5s= %.6e\n", nn->abundance_names[i], nn->Y_final[i]);
    }

    /* ---- Concatenated HT+MT+LT history, embedding each era's narrower
     * vector into the common n_lt columns by species name (mirrors
     * solve()'s _embed/Y_of_t construction; DT-era extension not ported,
     * see this module's header comment). MT/LT recorders' first row
     * duplicates the previous era's last time point exactly (both eras
     * are seeded at the boundary time), so it is dropped here exactly as
     * Python's sol_MT.t[1:]/sol_LT.t[1:] does. ---- */
    char ht_names[2][16] = { "n", "p" };
    nn->n_t = rec_ht.n + (rec_mt.n - 1) + (rec_lt.n - 1);
    nn->t_hist = CPR_XMALLOC(nn->n_t * sizeof(double));
    nn->Y_hist = CPR_XCALLOC(nn->n_t * n_lt, sizeof(double));
    size_t row = 0;
    for (size_t i = 0; i < rec_ht.n; i++, row++) {
        nn->t_hist[row] = rec_ht.t[i];
        embed_row(&nn->Y_hist[row * n_lt], nn->abundance_names, n_lt,
                  &rec_ht.Y[i * 2], ht_names, 2);
    }
    for (size_t i = 1; i < rec_mt.n; i++, row++) {
        nn->t_hist[row] = rec_mt.t[i];
        embed_row(&nn->Y_hist[row * n_lt], nn->abundance_names, n_lt,
                  &rec_mt.Y[i * n_mt], mt_names, n_mt);
    }
    for (size_t i = 1; i < rec_lt.n; i++, row++) {
        nn->t_hist[row] = rec_lt.t[i];
        embed_row(&nn->Y_hist[row * n_lt], nn->abundance_names, n_lt,
                  &rec_lt.Y[i * n_lt], lt_names, n_lt);
    }
    nn->t_start = nn->t_hist[0];

    free(Yi_MT); free(Yi_LT);
    recorder_free(&rec_ht); recorder_free(&rec_mt); recorder_free(&rec_lt);

    /* Propagate write failures. These return codes used to be discarded, so a
     * bad output_file (unwritable directory, full disk) left the run reporting
     * success with exit status 0 and no message -- while Python raised OSError
     * for the same input -- and leaked the strdup'd *errmsg nobody would ever
     * read or free (confirmed by `leaks`: 96 bytes). Failing the solve matches
     * Python and makes the caller free the message on the normal error path. */
    /* Unlike the HT/MT/LT failure paths above, these run with `nn` already
     * populated -- and this function's contract is that a nonzero return
     * leaves nothing for the caller to free (mc.c and api.c both just
     * propagate the error without calling cpr_nuclear_network_free). Release
     * it here so failing the run does not leak the abundance/history arrays
     * (~384 KB for the default network). */
    if (cfg->output_time_evolution
        && cpr_nuclear_network_write_time_evolution(nn, cfg->output_n_points, errmsg)) {
        cpr_nuclear_network_free(nn);
        return 1;
    }
    if (cfg->output_final_result
        && cpr_nuclear_network_write_final_result(nn, errmsg)) {
        cpr_nuclear_network_free(nn);
        return 1;
    }

    return 0;
}

void cpr_nuclear_network_free(CPRNuclearNetwork *nn)
{
    free(nn->abundance_names); free(nn->Y_final);
    free(nn->t_hist); free(nn->Y_hist);
    memset(nn, 0, sizeof(*nn));
}

double cpr_nuclear_network_get(const CPRNuclearNetwork *nn, const char *name)
{
    for (size_t i = 0; i < nn->n_species; i++)
        if (strcmp(nn->abundance_names[i], name) == 0) return nn->Y_final[i];
    return 0.0;
}

double cpr_nuclear_network_Y_of_t(const CPRNuclearNetwork *nn, const char *name, double t)
{
    size_t col = nn->n_species;
    for (size_t i = 0; i < nn->n_species; i++)
        if (strcmp(nn->abundance_names[i], name) == 0) { col = i; break; }
    if (col == nn->n_species) return 0.0; /* not tracked */

    if (t <= nn->t_start) return 0.0;                       /* before HT start */
    if (t >= nn->t_hist[nn->n_t - 1]) return nn->Y_final[col]; /* past LT end */

    /* Binary search the bracketing segment, then linear-interpolate --
     * mirrors interp1d's default (linear, here with the explicit
     * fill_value=(0, Y[-1]) clamps applied above instead of bounds_error). */
    size_t lo = 0, hi = nn->n_t - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (nn->t_hist[mid] <= t) lo = mid; else hi = mid;
    }
    double t0 = nn->t_hist[lo], t1 = nn->t_hist[hi];
    double y0 = nn->Y_hist[lo * nn->n_species + col], y1 = nn->Y_hist[hi * nn->n_species + col];
    double frac = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;
    return y0 + frac * (y1 - y0);
}

/* Creates every directory component of `path` in turn, including `path`
 * itself (mkdir -p equivalent without a shell call) -- mirrors
 * os.makedirs(exist_ok=True). The intermediate-component loop alone
 * (walking only embedded '/' characters) never creates the final,
 * slash-free component -- e.g. mkdir_p("results") used to silently do
 * nothing, so a fresh checkout's first write into the (not-yet-existing)
 * "results/" directory failed with ENOENT; the explicit mkdir() below the
 * loop covers that last, most common case. */
static void mkdir_p(const char *path)
{
    char buf[4300];
    snprintf(buf, sizeof(buf), "%s", path);
    for (char *p = buf + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            mkdir(buf, 0755);
            *p = '/';
        }
    }
    mkdir(buf, 0755);
}

int cpr_nuclear_network_write_final_result(const CPRNuclearNetwork *nn, char **errmsg)
{
    /* Resolve relative to the current working directory (matching
     * os.path.abspath's behaviour -- Python resolves against cwd, not
     * data_dir). */
    char path[4300];
    snprintf(path, sizeof(path), "%s", nn->cfg->output_final_file);

    char dir[4300];
    snprintf(dir, sizeof(dir), "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) { *slash = '\0'; mkdir_p(dir); }

    FILE *f = fopen(path, "w");
    if (!f) {
        char buf[4400];
        snprintf(buf, sizeof(buf), "cpr_nuclear_network_write_final_result: cannot open %s", path);
        *errmsg = strdup(buf);
        return 1;
    }
    fprintf(f, "# %-12sY\n", "nuclide");
    for (size_t i = 0; i < nn->n_species; i++)
        fprintf(f, "%-14s%.6e\n", nn->abundance_names[i], nn->Y_final[i]);
    fclose(f);
    printf("[output] Final abundances (%zu nuclides) written to %s\n", nn->n_species, path);
    fflush(stdout);
    return 0;
}

void cpr_nuclear_network_sample_time_evolution(const CPRNuclearNetwork *nn, int n_points,
                                                  double *t_out, double *T_out, double *a_out,
                                                  double *Tnue_out, double *Tnumu_out,
                                                  double *Tnutau_out, double *Y_out)
{
    const CPRConfig *cfg = nn->cfg;
    CPRBackground *bg = nn->background;

    double t_cosmo = cpr_bg_t_of_T(bg, cfg->T_start_cosmo_MeV);
    double t_end = nn->t_end;
    size_t n = (size_t)n_points;
    double logTlo = log10(t_cosmo), logThi = log10(t_end);

    for (size_t i = 0; i < n; i++) {
        double frac = (n == 1) ? 0.0 : (double)i / (double)(n - 1);
        double t = pow(10.0, logTlo + frac * (logThi - logTlo));
        t_out[i] = t;
        T_out[i] = cpr_bg_T_of_t(bg, t);
        a_out[i] = bg->has_scale_factor ? cpr_bg_a_of_t(bg, t) : NAN;

        double Tnue, Tnumu, Tnutau;
        if (cpr_bg_Tnu_of_t(bg, t, &Tnue, &Tnumu, &Tnutau)) {
            Tnue_out[i] = Tnue; Tnumu_out[i] = Tnumu; Tnutau_out[i] = Tnutau;
        } else {
            Tnue_out[i] = Tnumu_out[i] = Tnutau_out[i] = NAN;
        }

        for (size_t s = 0; s < nn->n_species; s++)
            Y_out[i * nn->n_species + s] = cpr_nuclear_network_Y_of_t(nn, nn->abundance_names[s], t);
    }
}

/* Whether this run emits per-reaction forward-rate columns: simply the flag
 * (one column per reaction in the active LT network, whatever the network/amax
 * selects). Mirrors Python's `if cfg.output_rates_time_evolution:` gate in
 * NuclearNetwork._write_time_evolution. */
static int rate_columns_enabled(const CPRConfig *cfg)
{
    return cfg->output_rates_time_evolution;
}

/* Build the sorted per-reaction rate-column list for the active LT network.
 * names[k] = "<reaction>_frwrd", rows[k] = the reaction's thermonuclear row
 * index into lt->fwd (= LT names index - 1, since names[0] is the prepended
 * weak n__p with no fwd row). Returns the column count (n LT reactions - 1).
 * Sorted by column name so the order matches Python's sorted() exactly.
 * Callers pass buffers of at least (lt->n_reac - 1) entries. */
static size_t build_rate_columns(const CPRNuclearNetwork *nn,
                                   char (*names)[64], size_t *rows)
{
    const CPRNetworkDef *lt = &nn->nucl->lt_net;
    size_t count = 0;
    for (size_t i = 1; i < lt->n_reac; i++) {   /* skip n__p at index 0 */
        snprintf(names[count], 64, "%s_frwrd", lt->names[i]);
        rows[count] = i - 1;                     /* fwd row for names[i] */
        count++;
    }
    /* Insertion sort by column name (small n ~ 12), keeping rows aligned --
     * strcmp order == Python's str sort over these ASCII names. */
    for (size_t a = 1; a < count; a++) {
        char kn[64];
        snprintf(kn, sizeof(kn), "%s", names[a]);
        size_t kr = rows[a];
        size_t b = a;
        while (b > 0 && strcmp(names[b - 1], kn) > 0) {
            snprintf(names[b], 64, "%s", names[b - 1]);
            rows[b] = rows[b - 1];
            b--;
        }
        snprintf(names[b], 64, "%s", kn);
        rows[b] = kr;
    }
    return count;
}

/* Active forward reaction rate of LT thermonuclear row `row` at photon
 * temperature `T_MeV`, linearly interpolated on the master T9 grid -- the exact
 * interpolation cpr_network_fill_buffer and Python's <rxn>_frwrd use
 * (T9 = T[K]/1e9, searchsorted-1 clamped, linear weight). */
static double frwrd_at(const CPRNetworkDef *lt, size_t row, double T_MeV)
{
    double T9 = T_MeV * cpr_MeV_to_Kelvin() * 1.0e-9;
    const double *g = lt->grid;
    size_t ii = cpr_find_segment(g, lt->n_grid, T9);
    double w = (T9 - g[ii]) / (g[ii + 1] - g[ii]);
    return lt->fwd[row * lt->n_grid + ii] * (1.0 - w)
         + lt->fwd[row * lt->n_grid + ii + 1] * w;
}

size_t cpr_nuclear_network_rate_columns(const CPRNuclearNetwork *nn,
                                          char (*out_names)[64])
{
    if (!rate_columns_enabled(nn->cfg))
        return 0;
    size_t n = nn->nucl->lt_net.n_reac - 1;   /* one column per thermonuclear reaction */
    if (out_names) {
        size_t *rows = CPR_XMALLOC(n * sizeof(size_t));
        build_rate_columns(nn, out_names, rows);
        free(rows);
    }
    return n;
}

void cpr_nuclear_network_sample_rates(const CPRNuclearNetwork *nn,
                                        const double *T_MeV, int n_points,
                                        double *rates_out)
{
    if (!rate_columns_enabled(nn->cfg))
        return;
    const CPRNetworkDef *lt = &nn->nucl->lt_net;
    size_t n_cols = lt->n_reac - 1;
    char (*names)[64] = CPR_XMALLOC(n_cols * sizeof(*names));
    size_t *rows = CPR_XMALLOC(n_cols * sizeof(size_t));
    build_rate_columns(nn, names, rows);
    for (size_t i = 0; i < (size_t)n_points; i++)
        for (size_t k = 0; k < n_cols; k++)
            rates_out[i * n_cols + k] = frwrd_at(lt, rows[k], T_MeV[i]);
    free(names);
    free(rows);
}

int cpr_nuclear_network_write_time_evolution(const CPRNuclearNetwork *nn, int n_points,
                                                char **errmsg)
{
    /* cfg->output_file == NULL/"" is the in-memory-only escape hatch
     * (mirrors Python's NuclearNetwork._write_time_evolution skipping disk
     * I/O when cfg.output_file is None, e.g. primat-gui/run_bbn's
     * in-memory-only callers via CPRResults's evol_* arrays, populated by
     * cpr_assemble_results regardless of this flag). */
    if (!nn->cfg->output_file || !nn->cfg->output_file[0])
        return 0;

    size_t n = (size_t)n_points;
    double *t_out = CPR_XMALLOC(n * sizeof(double));
    double *T_out = CPR_XMALLOC(n * sizeof(double));
    double *a_out = CPR_XMALLOC(n * sizeof(double));
    double *Tnue_out = CPR_XMALLOC(n * sizeof(double));
    double *Tnumu_out = CPR_XMALLOC(n * sizeof(double));
    double *Tnutau_out = CPR_XMALLOC(n * sizeof(double));
    double *Y_out = CPR_XMALLOC(n * nn->n_species * sizeof(double));
    cpr_nuclear_network_sample_time_evolution(nn, n_points, t_out, T_out, a_out,
                                                Tnue_out, Tnumu_out, Tnutau_out, Y_out);

    /* Optional per-reaction forward-rate columns (flag on; one per active LT
     * reaction), sampled
     * at the same T_gamma grid -- appended after the Y_ block, matching the
     * Python backend's dump_evolution order. */
    size_t n_rate = cpr_nuclear_network_rate_columns(nn, NULL);
    char (*rate_names)[64] = NULL;
    double *rate_out = NULL;
    if (n_rate) {
        rate_names = CPR_XMALLOC(n_rate * sizeof(*rate_names));
        rate_out = CPR_XMALLOC(n * n_rate * sizeof(double));
        cpr_nuclear_network_rate_columns(nn, rate_names);
        cpr_nuclear_network_sample_rates(nn, T_out, n_points, rate_out);
    }

    const char *rel = nn->cfg->output_file;
    char path[4300];
    snprintf(path, sizeof(path), "%s", rel);
    char dir[4300];
    snprintf(dir, sizeof(dir), "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) { *slash = '\0'; mkdir_p(dir); }

    FILE *f = fopen(path, "w");
    if (!f) {
        free(t_out); free(T_out); free(a_out);
        free(Tnue_out); free(Tnumu_out); free(Tnutau_out); free(Y_out);
        free(rate_names); free(rate_out);
        char buf[4400];
        snprintf(buf, sizeof(buf), "cpr_nuclear_network_write_time_evolution: cannot open %s", path);
        *errmsg = strdup(buf);
        return 1;
    }

    /* Unified schema, header-compatible with
     * primat.evolution.dump_evolution/load_evolution: no leading "#",
     * tab-separated, t_s/a/T_*_MeV core block then one Y_<nuclide> column
     * per tracked species, then the optional per-reaction <reaction>_frwrd
     * columns (cfg->output_rates_time_evolution; one per active LT reaction)
     * in the SAME sorted order the Python backend writes. */
    fprintf(f, "t_s\ta\tT_gamma_MeV\tT_nue_MeV\tT_numu_MeV\tT_nutau_MeV");
    for (size_t s = 0; s < nn->n_species; s++) fprintf(f, "\tY_%s", nn->abundance_names[s]);
    for (size_t k = 0; k < n_rate; k++) fprintf(f, "\t%s", rate_names[k]);
    fprintf(f, "\n");

    for (size_t i = 0; i < n; i++) {
        fprintf(f, "%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e",
                t_out[i], a_out[i], T_out[i], Tnue_out[i], Tnumu_out[i], Tnutau_out[i]);
        for (size_t s = 0; s < nn->n_species; s++)
            fprintf(f, "\t%.8e", Y_out[i * nn->n_species + s]);
        for (size_t k = 0; k < n_rate; k++)
            fprintf(f, "\t%.8e", rate_out[i * n_rate + k]);
        fprintf(f, "\n");
    }
    fclose(f);
    free(t_out); free(T_out); free(a_out);
    free(Tnue_out); free(Tnumu_out); free(Tnutau_out); free(Y_out);
    free(rate_names); free(rate_out);
    printf("[output] Time-evolution data (%zu rows) written to %s\n", n, path);
    fflush(stdout);
    return 0;
}

/* =====================================================================
 * Decay Time (DT) era -- port of nuclear_network.py's _build_decay_matrix
 * / _integrate_decay_era / _write_decay_evolution (see CLAUDE.md's
 * backend feature gaps). After BBN ends at t_end, long-lived radioactive
 * isotopes (C14, Be10, Na22, the residual free neutron, ...) keep decaying
 * for years to Myr under the *constant* decay-rate matrix D alone (no Hubble
 * expansion, no thermal production, since T is far too low for any thermal
 * activation): dY/dt = D.Y, solved exactly by Y(t) = exp(D*(t-t_end)) Y0.
 *
 * Gated (mirroring Python) on cfg->decay_era plus the LT network actually
 * carrying a decay reaction -- weak_flags holds n__p (index 0) plus every
 * decay, so "any weak reaction past index 0" is exactly that test. It used to
 * be keyed on the literal network name `large`, which silently skipped the era
 * for a large-equivalent network reproduced under a renamed user_nuclear_dir
 * overlay.
 *
 * DIVERGENCE FROM PYTHON, deliberate: this port additionally requires
 * cfg->output_decay_evolution, i.e. it runs only when the TSV was requested.
 * Y_final and the result dict are the end-of-LT state on both backends, so no
 * BBN *observable* differs. But Python does more than write the file: its
 * solve() also extends the public Y_of_t interpolator across the DT era, so
 * run[species](t) past t_end returns the decayed state there. That accessor
 * has no representation across this ABI -- there is nothing for C to expose it
 * through -- so computing 200 dense matrix exponentials whose result nothing
 * could read would be pure waste. If a future revision surfaces a post-t_end
 * abundance history through the C API, drop the output_decay_evolution term
 * from the gate below so the two triggers match again.
 * ===================================================================== */

/* Dense n x n matrix multiply C = A*B, row-major. n is small (<= ~60, the
 * large network's LT nuclide count), so the naive triple loop is fine. */
static void dt_mat_mul(const double *A, const double *B, double *C, size_t n)
{
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) s += A[i * n + k] * B[k * n + j];
            C[i * n + j] = s;
        }
}

/* Dense matrix exponential E = expm(A) (A n x n row-major, left unchanged),
 * via the scaling-and-squaring + degree-13 Padé algorithm of Higham (2005) --
 * the same method scipy.linalg.expm uses, and the reason Python's
 * _integrate_decay_era can handle D's ~16-decade eigenvalue spread (the
 * fastest residual decay vs. a ~Gyr Delta-t) in milliseconds: the squaring
 * count grows only logarithmically in ||A||. Degree 13 alone (scipy also
 * uses lower degrees for small ||A|| purely as an optimisation) reproduces
 * exp to ~machine precision at every norm once scaled below theta13, far
 * finer than the cross-backend tolerance. Returns 0 on success, nonzero
 * with *errmsg set (caller frees) on OOM or a singular Pade denominator.
 *
 * Reference: N. J. Higham, "The Scaling and Squaring Method for the Matrix
 * Exponential Revisited", SIAM J. Matrix Anal. Appl. 26 (2005) 1179-1193. */
static int dt_expm(const double *A_in, size_t n, double *E, char **errmsg)
{
    /* Degree-13 Pade numerator/denominator coefficients (Higham 2005 Table). */
    static const double b[14] = {
        64764752532480000.0, 32382376266240000.0, 7771770303897600.0,
        1187353796428800.0, 129060195264000.0, 10559470521600.0,
        670442572800.0, 33522128640.0, 1323241920.0, 40840800.0,
        960960.0, 16380.0, 182.0, 1.0
    };
    const double theta13 = 5.371920351148152; /* scaling threshold for degree 13 */

    size_t nn = n * n;
    /* One workspace block, carved into the matrices the algorithm needs. */
    double *A  = malloc(nn * sizeof(double));
    double *A2 = malloc(nn * sizeof(double));
    double *A4 = malloc(nn * sizeof(double));
    double *A6 = malloc(nn * sizeof(double));
    double *U  = malloc(nn * sizeof(double));
    double *V  = malloc(nn * sizeof(double));
    double *W  = malloc(nn * sizeof(double));  /* scratch */
    double *P  = malloc(nn * sizeof(double));  /* V + U (numerator) */
    double *Q  = malloc(nn * sizeof(double));  /* V - U (denominator) */
    size_t *piv = malloc(n * sizeof(size_t));
    double *col = malloc(n * sizeof(double));
    if (!A || !A2 || !A4 || !A6 || !U || !V || !W || !P || !Q || !piv || !col) {
        free(A); free(A2); free(A4); free(A6); free(U); free(V);
        free(W); free(P); free(Q); free(piv); free(col);
        if (errmsg) *errmsg = strdup("dt_expm: out of memory");
        return 1;
    }
    memcpy(A, A_in, nn * sizeof(double));

    /* 1-norm (max absolute column sum) sets the scaling s so ||A/2^s|| <= theta13. */
    double norm1 = 0.0;
    for (size_t j = 0; j < n; j++) {
        double colsum = 0.0;
        for (size_t i = 0; i < n; i++) colsum += fabs(A[i * n + j]);
        if (colsum > norm1) norm1 = colsum;
    }
    int s = 0;
    if (norm1 > theta13) {
        s = (int)ceil(log2(norm1 / theta13));
        if (s < 0) s = 0;
        double scale = ldexp(1.0, -s); /* 2^-s */
        for (size_t i = 0; i < nn; i++) A[i] *= scale;
    }

    /* Even powers of the (scaled) A. */
    dt_mat_mul(A, A, A2, n);
    dt_mat_mul(A2, A2, A4, n);
    dt_mat_mul(A4, A2, A6, n);

    /* U = A * (A6*(b13*A6 + b11*A4 + b9*A2) + b7*A6 + b5*A4 + b3*A2 + b1*I)
     * V =      A6*(b12*A6 + b10*A4 + b8*A2) + b6*A6 + b4*A4 + b2*A2 + b0*I */
    for (size_t i = 0; i < nn; i++)
        W[i] = b[13] * A6[i] + b[11] * A4[i] + b[9] * A2[i];
    dt_mat_mul(A6, W, V, n);                 /* reuse V as scratch for A6*W */
    for (size_t i = 0; i < nn; i++)
        V[i] += b[7] * A6[i] + b[5] * A4[i] + b[3] * A2[i];
    for (size_t i = 0; i < n; i++)
        V[i * n + i] += b[1];                /* + b1*I */
    dt_mat_mul(A, V, U, n);                  /* U = A * (...) */

    for (size_t i = 0; i < nn; i++)
        W[i] = b[12] * A6[i] + b[10] * A4[i] + b[8] * A2[i];
    dt_mat_mul(A6, W, V, n);                 /* V = A6*(...) */
    for (size_t i = 0; i < nn; i++)
        V[i] += b[6] * A6[i] + b[4] * A4[i] + b[2] * A2[i];
    for (size_t i = 0; i < n; i++)
        V[i * n + i] += b[0];                /* + b0*I */

    /* Solve (V - U) R = (V + U). R = E (result before squaring). */
    for (size_t i = 0; i < nn; i++) { P[i] = V[i] + U[i]; Q[i] = V[i] - U[i]; }
    if (cpr_lu_factor(Q, n, piv)) {
        free(A); free(A2); free(A4); free(A6); free(U); free(V);
        free(W); free(P); free(Q); free(piv); free(col);
        if (errmsg) *errmsg = strdup("dt_expm: singular Pade denominator");
        return 1;
    }
    /* Solve column by column: Q * E[:,j] = P[:,j]. */
    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < n; i++) col[i] = P[i * n + j];
        cpr_lu_solve(Q, n, piv, col);
        for (size_t i = 0; i < n; i++) E[i * n + j] = col[i];
    }

    /* Undo the scaling: square E s times (E <- E*E). */
    for (int k = 0; k < s; k++) {
        dt_mat_mul(E, E, W, n);
        memcpy(E, W, nn * sizeof(double));
    }

    free(A); free(A2); free(A4); free(A6); free(U); free(V);
    free(W); free(P); free(Q); free(piv); free(col);
    return 0;
}

/* Build the constant N x N decay-rate matrix D [s^-1] for the DT era (port
 * of _build_decay_matrix). dY/dt = D.Y, D's columns indexing the *parent*:
 * each decay X -> products with rate lambda contributes D[X,X] -= lambda*mult_X
 * and D[P,X] += lambda*mult_P.
 *
 * CONVENTION: Y is the *number* abundance per baryon, Y_s = n_s/n_B with
 * sum_s A_s Y_s = 1 -- not a mass fraction. That is what the LT/MT right-hand
 * side uses (cpr_network_builder applies the bare integer stoichiometry
 * c_prod - c_react, and the conservation check verifies sum_s A_s dY_s = 0).
 * The product gain is therefore the bare multiplicity, with NO A_P/A_X factor.
 * An earlier version carried such a factor *in addition to* mult_P, which
 * broke baryon conservation for every decay whose products differ in mass
 * number from the parent (Li8 -> a+a yielded one alpha instead of two; C9 ->
 * a+a+p lost 4/9 of the baryon number). The 33 ordinary beta decays have
 * A_P == A_X and were unaffected, which is why it went unnoticed. Kept in
 * lockstep with nuclear_network.py's _build_decay_matrix.
 *
 * Decay reactions are the LT network's weak reactions other than n__p
 * (index 0, handled by the HT/MT/LT thermal weak rate); their rate is the
 * T9-independent decays.txt constant, stored uniformly across the master T9
 * grid, so grid index 0 is representative. The free-neutron beta decay
 * n -> p (lambda = 1/tau_n) is added explicitly -- it is the T->0 limit of
 * the thermal n<->p rate, absent from decays.txt, so without it the residual
 * free neutrons at t_end would never decay. D is written into the caller's
 * N*N buffer (row-major, zeroed here first). */
static void dt_build_decay_matrix(const CPRNetworkDef *net, const CPRConfig *cfg, double *D)
{
    size_t N = net->n_species;
    memset(D, 0, N * N * sizeof(double));

    for (size_t rxn = 1; rxn < net->n_reac; rxn++) {   /* skip index 0 (n__p) */
        if (!net->weak_flags[rxn]) continue;           /* decays are the weak reactions */
        double rate = net->fwd_median[(rxn - 1) * net->n_grid + 0]; /* [s^-1], T9-independent */
        if (rate == 0.0) continue;

        const CPRStoichSide *react = &net->network[rxn].reactants;
        const CPRStoichSide *prod  = &net->network[rxn].products;
        for (size_t a = 0; a < react->n; a++) {
            long X = react->species_idx[a];
            double X_mult = (double)react->mult[a];
            D[(size_t)X * N + (size_t)X] -= rate * X_mult;   /* parent loss */
            for (size_t p = 0; p < prod->n; p++) {
                long P = prod->species_idx[p];
                double P_mult = (double)prod->mult[p];
                /* number-abundance gain: the bare multiplicity (see the
                 * CONVENTION note above -- no A_P/A_X weighting). */
                D[(size_t)P * N + (size_t)X] += rate * P_mult;
            }
        }
    }

    /* Free-neutron beta decay n -> p (rate 1/tau_n); A_n = A_p = 1. */
    long n_idx = -1, p_idx = -1;
    for (size_t i = 0; i < N; i++) {
        if (strcmp(net->species[i], "n") == 0) n_idx = (long)i;
        else if (strcmp(net->species[i], "p") == 0) p_idx = (long)i;
    }
    if (n_idx >= 0 && p_idx >= 0) {
        double lam_n = 1.0 / cfg->tau_n;
        D[(size_t)n_idx * N + (size_t)n_idx] -= lam_n;
        D[(size_t)p_idx * N + (size_t)n_idx] += lam_n;
    }
}

int cpr_nuclear_network_decay_era(const CPRNuclearNetwork *nn, char **errmsg)
{
    const CPRConfig *cfg = nn->cfg;
    const CPRNetworkDef *net = &nn->nucl->lt_net;

    /* See this section's top comment for the gate and its one deliberate
     * divergence from Python. "Has decays" = any weak reaction past index 0
     * (index 0 is n__p, handled by the thermal weak rate, not decays.txt). */
    int has_decays = 0;
    for (size_t i = 1; i < net->n_reac; i++)
        if (net->weak_flags[i]) { has_decays = 1; break; }
    if (!cfg->decay_era || !has_decays || !cfg->output_decay_evolution)
        return 0;

    size_t N = net->n_species;
    int M = cfg->decay_n_points;
    if (M < 1) return 0;

    double *D      = malloc(N * N * sizeof(double));
    double *Y0     = malloc(N * sizeof(double));
    double *t_grid = malloc((size_t)M * sizeof(double));
    double *Y_DT   = malloc((size_t)M * N * sizeof(double));
    double *E      = malloc(N * N * sizeof(double));
    double *Ddt    = malloc(N * N * sizeof(double));
    if (!D || !Y0 || !t_grid || !Y_DT || !E || !Ddt) {
        free(D); free(Y0); free(t_grid); free(Y_DT); free(E); free(Ddt);
        if (errmsg) *errmsg = strdup("cpr_nuclear_network_decay_era: out of memory");
        return 1;
    }

    dt_build_decay_matrix(net, cfg, D);

    /* Y0 = end-of-LT abundances, in abundance_names (== lt_net.species) order. */
    for (size_t i = 0; i < N; i++) Y0[i] = nn->Y_final[i];

    /* Time grid log-spaced in the *elapsed* time Delta-t = t - t_end from
     * Delta-t = 1 s to t_decay_end, then offset to absolute cosmic time --
     * mirrors solve()'s `t_end + np.logspace(log10(1), log10(t_decay_end),
     * decay_n_points)`. Spacing in Delta-t (not absolute t) is essential to
     * resolve the fast residual free-neutron decay (tau_n ~ 880 s), ~10
     * decades below t_end (~1.3e6 s). */
    double log_lo = log10(1.0), log_hi = log10(cfg->t_decay_end);
    for (int k = 0; k < M; k++) {
        double frac = (M == 1) ? 0.0 : (double)k / (double)(M - 1);
        t_grid[k] = nn->t_end + pow(10.0, log_lo + frac * (log_hi - log_lo));
    }

    /* Y(t_k) = expm(D * (t_k - t_end)) @ Y0, per output time. */
    for (int k = 0; k < M; k++) {
        double dt = t_grid[k] - nn->t_end;
        for (size_t i = 0; i < N * N; i++) Ddt[i] = D[i] * dt;
        if (dt_expm(Ddt, N, E, errmsg)) {
            free(D); free(Y0); free(t_grid); free(Y_DT); free(E); free(Ddt);
            return 1;
        }
        for (size_t i = 0; i < N; i++) {
            double acc = 0.0;
            for (size_t j = 0; j < N; j++) acc += E[i * N + j] * Y0[j];
            /* Clip tiny negatives from matrix-exp cancellation (mirrors
             * np.clip(..., 0, None) in _integrate_decay_era). */
            Y_DT[(size_t)k * N + i] = acc < 0.0 ? 0.0 : acc;
        }
    }

    /* Write the decay-evolution TSV (port of _write_decay_evolution):
     * header "t\tY<species>...", then one row per time point. np.savetxt's
     * default %.18e format is matched so a loader sees identical precision. */
    char path[4300];
    snprintf(path, sizeof(path), "%s", cfg->output_decay_file);
    char dir[4300];
    snprintf(dir, sizeof(dir), "%s", path);
    char *slash = strrchr(dir, '/');
    if (slash) { *slash = '\0'; mkdir_p(dir); }

    FILE *f = fopen(path, "w");
    if (!f) {
        free(D); free(Y0); free(t_grid); free(Y_DT); free(E); free(Ddt);
        char buf[4400];
        snprintf(buf, sizeof(buf), "cpr_nuclear_network_decay_era: cannot open %s", path);
        if (errmsg) *errmsg = strdup(buf);
        return 1;
    }
    fprintf(f, "t");
    for (size_t s = 0; s < N; s++) fprintf(f, "\tY%s", nn->abundance_names[s]);
    fprintf(f, "\n");
    for (int k = 0; k < M; k++) {
        fprintf(f, "%.18e", t_grid[k]);
        for (size_t s = 0; s < N; s++) fprintf(f, "\t%.18e", Y_DT[(size_t)k * N + s]);
        fprintf(f, "\n");
    }
    fclose(f);
    printf("[output] Decay-era evolution (%d rows) written to %s\n", M, path);
    fflush(stdout);

    free(D); free(Y0); free(t_grid); free(Y_DT); free(E); free(Ddt);
    return 0;
}
