/* vegas.c -- see vegas.h.
 *
 * Implements the core of Lepage's VEGAS algorithm (G.P. Lepage, J. Comput.
 * Phys. 27 (1978) 192): per-dimension adaptive importance sampling, with
 * the grid in each dimension refined after every iteration from that
 * iteration's f^2-weighted histogram, damped by VEGAS_DAMPING and equal-mass
 * rebinned (the textbook formulation, as in Numerical Recipes' `vegas` or
 * GSL's `gsl_monte_vegas`).
 *
 * Simplification relative to full Lepage VEGAS: samples are drawn by
 * picking a bin independently and uniformly in each dimension (not jointly
 * stratified across the ng^ndim grid of cells). This keeps the
 * implementation simple and is still an unbiased, correctly-adapted
 * importance-sampling estimator -- it just forgoes the extra variance
 * reduction joint stratification would give in higher dimensions. At
 * ndim=2 with the evaluation budgets used here (cfg.vegas_n_eval ~ 2e4 per
 * iteration), this is more than sufficient for the ~1e-3 relative accuracy
 * the CCRTh thermal correction already targets (see weak_rates.c).
 */
#include "vegas.h"
#include "rng.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define VEGAS_NDIM 2
#define VEGAS_NBINS 50           /* bins per dimension */
/* Grid-refinement damping exponent. 1.5 is the GSL / Numerical Recipes
 * default; Lepage's own `vegas` Python package -- the one
 * weak_rates/corrections.py imports -- defaults to alpha=0.5 instead, so the
 * two backends adapt their grids at different rates. That affects how fast
 * each converges, not what it converges to, and the CCRTh caches are
 * fingerprinted per backend; VEGAS_DAMPING is deliberately NOT in that fingerprint,
 * so changing it here would silently stale the shipped tables. */
#define VEGAS_DAMPING 1.5
#define VEGAS_TINY 1.0e-300

typedef struct {
    double edges[VEGAS_NDIM][VEGAS_NBINS + 1];   /* bin edges in real (domain) coordinates */
    double hist[VEGAS_NDIM][VEGAS_NBINS];        /* per-iteration f^2-weight accumulator */
} CPRVegasGrid;

static void grid_init(CPRVegasGrid *g, const double lo[VEGAS_NDIM], const double hi[VEGAS_NDIM])
{
    for (int d = 0; d < VEGAS_NDIM; d++) {
        for (int i = 0; i <= VEGAS_NBINS; i++)
            g->edges[d][i] = lo[d] + (hi[d] - lo[d]) * ((double)i / VEGAS_NBINS);
        memset(g->hist[d], 0, sizeof(g->hist[d]));
    }
}

/* Equal-mass rebin of one dimension's grid from its accumulated histogram,
 * following the standard VEGAS smoothing + damping + rebinning recipe. */
static void refine_dim(double edges[VEGAS_NBINS + 1], const double hist_in[VEGAS_NBINS])
{
    double d[VEGAS_NBINS];

    /* Smooth the raw f^2 histogram (three-point running average, with
     * two-point averages at the two edge bins) -- damps single-bin noise
     * spikes that would otherwise overreact the next grid refinement. */
    if (VEGAS_NBINS == 1) {
        d[0] = hist_in[0];
    } else {
        d[0] = (hist_in[0] + hist_in[1]) / 2.0;
        for (int i = 1; i < VEGAS_NBINS - 1; i++)
            d[i] = (hist_in[i - 1] + hist_in[i] + hist_in[i + 1]) / 3.0;
        d[VEGAS_NBINS - 1] = (hist_in[VEGAS_NBINS - 2] + hist_in[VEGAS_NBINS - 1]) / 2.0;
    }

    double sum = 0.0;
    for (int i = 0; i < VEGAS_NBINS; i++) sum += d[i];
    if (sum <= 0.0) return;  /* f was exactly 0 everywhere this iteration: keep the old grid */
    for (int i = 0; i < VEGAS_NBINS; i++) d[i] /= sum;  /* normalise to a probability per bin */

    /* Damped importance weight per old bin: rc[i] grows with how far above
     * (or below) uniform (1/VEGAS_NBINS) that bin's probability mass is, raised to
     * VEGAS_DAMPING -- bins seeing more of the integrand's f^2 mass get more of the
     * VEGAS_NBINS new bins allocated to them next iteration. */
    double rc[VEGAS_NBINS], sum_rc = 0.0;
    for (int i = 0; i < VEGAS_NBINS; i++) {
        if (d[i] < VEGAS_TINY) {
            rc[i] = VEGAS_TINY;
        } else if (d[i] == 1.0) {
            rc[i] = 1.0;  /* (d-1)/log(d) -> 1 in this limit */
        } else {
            rc[i] = pow((d[i] - 1.0) / log(d[i]), VEGAS_DAMPING);
        }
        sum_rc += rc[i];
    }

    /* Equal-mass rebin: walk the old bins accumulating rc, placing each new
     * edge where the cumulative rc mass crosses a multiple of sum_rc/VEGAS_NBINS,
     * linearly interpolating the edge position within whichever old bin
     * that crossing falls in. */
    double old_edges[VEGAS_NBINS + 1];
    memcpy(old_edges, edges, sizeof(old_edges));
    double target = sum_rc / VEGAS_NBINS;
    int old_i = 0;
    double acc = 0.0;  /* rc mass of bin old_i already consumed */
    double new_edges[VEGAS_NBINS + 1];
    new_edges[0] = old_edges[0];
    new_edges[VEGAS_NBINS] = old_edges[VEGAS_NBINS];
    for (int new_i = 1; new_i < VEGAS_NBINS; new_i++) {
        double needed = target;
        while (old_i < VEGAS_NBINS - 1 && rc[old_i] - acc < needed) {
            needed -= (rc[old_i] - acc);
            old_i++;
            acc = 0.0;
        }
        acc += needed;
        double frac = (rc[old_i] > 0.0) ? (acc / rc[old_i]) : 1.0;
        if (frac > 1.0) frac = 1.0;
        new_edges[new_i] = old_edges[old_i] + frac * (old_edges[old_i + 1] - old_edges[old_i]);
    }
    memcpy(edges, new_edges, sizeof(new_edges));
}

/* Draws one sample point, returning its domain coordinates in x[VEGAS_NDIM] and
 * the Jacobian |dx/dy| of the grid mapping at that point (the product,
 * over dimensions, of VEGAS_NBINS * (local bin width) -- see vegas.h). */
static double sample_point(const CPRVegasGrid *g, CPRRng *rng, double x[VEGAS_NDIM], int iy[VEGAS_NDIM])
{
    double jac = 1.0;
    for (int d = 0; d < VEGAS_NDIM; d++) {
        double u = cpr_rng_uniform(rng) * VEGAS_NBINS;
        int i = (int)u;
        if (i >= VEGAS_NBINS) i = VEGAS_NBINS - 1;
        double frac = u - i;
        double lo = g->edges[d][i], hi = g->edges[d][i + 1];
        x[d] = lo + frac * (hi - lo);
        jac *= VEGAS_NBINS * (hi - lo);
        iy[d] = i;
    }
    return jac;
}

/* Runs one VEGAS iteration of n_eval samples: accumulates the histogram for
 * the next refine_dim() call (if adapt is set) and returns this iteration's
 * (mean, variance-of-the-mean) estimate of the integral. */
static void run_iteration(CPRVegasGrid *g, CPRVegasFunc f, void *ctx, CPRRng *rng,
                           int n_eval, int adapt, double *mean_out, double *var_out)
{
    if (adapt) {
        for (int d = 0; d < VEGAS_NDIM; d++) memset(g->hist[d], 0, sizeof(g->hist[d]));
    }

    double sum = 0.0, sum_sq = 0.0;
    for (int n = 0; n < n_eval; n++) {
        double x[VEGAS_NDIM];
        int iy[VEGAS_NDIM];
        double jac = sample_point(g, rng, x, iy);
        double g_val = f(x, ctx) * jac;
        sum += g_val;
        sum_sq += g_val * g_val;
        if (adapt) {
            double w = g_val * g_val;
            for (int d = 0; d < VEGAS_NDIM; d++) g->hist[d][iy[d]] += w;
        }
    }

    double mean = sum / n_eval;
    double var_of_mean = (sum_sq / n_eval - mean * mean) / n_eval;
    if (var_of_mean < 0.0) var_of_mean = 0.0;  /* guards a tiny negative from FP cancellation */

    if (adapt) {
        for (int d = 0; d < VEGAS_NDIM; d++) refine_dim(g->edges[d], g->hist[d]);
    }

    *mean_out = mean;
    *var_out = var_of_mean;
}

CPRVegasResult cpr_vegas_integrate(CPRVegasFunc f, void *ctx,
                                    const double lo[2], const double hi[2],
                                    int n_eval, int n_itn_warmup,
                                    int n_itn_measure, uint64_t seed)
{
    CPRVegasGrid grid;
    grid_init(&grid, lo, hi);

    CPRRng rng;
    cpr_rng_seed(&rng, seed);

    for (int it = 0; it < n_itn_warmup; it++) {
        double mean, var;
        run_iteration(&grid, f, ctx, &rng, n_eval, 1, &mean, &var);
    }

    /* Inverse-variance-weighted combination of the measure-phase
     * iterations, mirroring vegas.Integrator's running weighted average
     * (and corrections.py's reliance on `result['myres'].mean`). */
    double sum_w = 0.0, sum_w_mean = 0.0;
    for (int it = 0; it < n_itn_measure; it++) {
        double mean, var;
        run_iteration(&grid, f, ctx, &rng, n_eval, 1, &mean, &var);
        double w = 1.0 / (var > VEGAS_TINY ? var : VEGAS_TINY);
        sum_w += w;
        sum_w_mean += w * mean;
    }

    CPRVegasResult result;
    if (sum_w > 0.0) {
        result.mean = sum_w_mean / sum_w;
        result.sigma = sqrt(1.0 / sum_w);
    } else {
        result.mean = 0.0;
        result.sigma = 0.0;
    }
    return result;
}
