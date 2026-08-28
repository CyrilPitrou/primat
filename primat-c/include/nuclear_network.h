/* nuclear_network.h -- integrating the reaction network down the three
 * temperature eras, which is where the abundances come from.
 *
 * CPRNuclearNetwork is driven purely through the query surface of a
 * CPRBackground (cpr_bg_T_of_t/t_of_T/rhoB_BBN/weak_nTOp/pTOn, background.h)
 * and a CPRNuclearRates (the compiled MT/LT RHS/Jacobian kernels,
 * network_data.h) -- it knows nothing about *how* the background or the rate
 * tables were built, which is what lets a caller substitute either.
 *
 * cpr_nuclear_network_solve integrates:
 *   HT  (T > T_weak ~ 1 MeV):        n <-> p only, non-stiff RK45 (ode_rk.h).
 *   MT  (T_weak -> T_nucl ~ 0.11 MeV): the network's intersection with
 *                                       CPR_ORDER_MT -- 18 reactions for
 *                                       "large", 13 for "small" --
 *                                       (nr->mt_net/mt_compiled), stiff BDF
 *                                       (ode_bdf.h) with the analytic
 *                                       Jacobian.
 *   LT  (T_nucl -> T_end ~ 0.001 MeV): the chosen network
 *                                       (nr->lt_net/lt_compiled), stiff BDF.
 *
 * A fourth, optional stage runs past T_end: the Decay-Time era, which
 * propagates the long-lived isotopes forward under a constant decay matrix.
 * It changes no BBN observable and exists only to produce its own output file
 * -- see cpr_nuclear_network_decay_era below.
 *
 * Reference: Pitrou, Coc, Uzan & Vangioni, Phys. Rep. 2018 (arXiv:1806.11095).
 */
#ifndef CPRIMAT_NUCLEAR_NETWORK_H
#define CPRIMAT_NUCLEAR_NETWORK_H

#include "config.h"
#include "background.h"
#include "network_data.h"
#include <stddef.h>

typedef struct {
    const CPRConfig *cfg;       /* borrowed */
    CPRBackground *background;  /* borrowed; not const since cpr_bg_* query
                                  * functions are logically const but the
                                  * header declares them on non-const
                                  * CPRBackground* in a couple of spots */
    CPRNuclearRates *nucl;       /* borrowed */

    /* ---- Final state (populated by cpr_nuclear_network_solve). ---- */
    char (*abundance_names)[16]; /* owned copy of nr->lt_net.species, length n_species */
    size_t n_species;
    double *Y_final;             /* length n_species, same order as abundance_names */
    double t_end;                 /* cosmic time [s] at T_end (end of LT era) */

    /* ---- Concatenated HT+MT+LT abundance history, for
     * cpr_nuclear_network_Y_of_t / cpr_nuclear_network_write_time_evolution.
     * Row-major (n_t x n_species), in abundance_names column order; each
     * era's solution is embedded into the common columns by species name
     * (mirrors solve()'s _embed/Y_of_t, minus the DT-era extension). ---- */
    double *t_hist;   /* length n_t, strictly ascending */
    double *Y_hist;   /* length n_t * n_species */
    size_t n_t;
    double t_start;   /* t_hist[0]; cosmic time at T_start (10 MeV), kept
                        * separately since write_time_evolution needs the
                        * *nuclear-network* start time, not T_start_cosmo's
                        * (generally earlier) background start. */
} CPRNuclearNetwork;

/* Integrates the HT->MT->LT eras and populates `nn` (zeroed first).
 * `nucl` must already be cpr_nuclear_rates_init'd and have had
 * cpr_nuclear_rates_apply_variations applied for the current cfg (mirrors
 * solve()'s own nucl.apply_variations(cfg) call at the top -- done here,
 * not by the caller, for exact parity). `background` must already be
 * cpr_bg_init_standard/custom'd. Both are borrowed and must outlive `nn`.
 *
 * Returns 0 on success (caller must cpr_nuclear_network_free), nonzero
 * with *errmsg set (caller frees) on an ODE integration failure (mirrors
 * solve_ivp returning status != 0 in Python, which that code does not
 * itself check -- here a failed era integration is treated as a hard
 * error rather than silently returning a wrong answer). */
int cpr_nuclear_network_solve(CPRNuclearNetwork *nn, const CPRConfig *cfg,
                                CPRNuclearRates *nucl, CPRBackground *background,
                                char **errmsg);
void cpr_nuclear_network_free(CPRNuclearNetwork *nn);

/* Final abundance per baryon Y of nuclide `name` (mirrors Y_final.get(name, 0.0)
 * via PRIMAT's __getitem__); 0.0 if `name` is not tracked by the active network. */
double cpr_nuclear_network_get(const CPRNuclearNetwork *nn, const char *name);

/* Abundance vector Y(t) of nuclide `name` at cosmic time t [s], linearly
 * interpolated over nn->t_hist/Y_hist (mirrors Y_of_t's interp1d); 0
 * before nn->t_start, held at the final value of `name` (cpr_nuclear_network_get)
 * beyond the last grid point (constant extrapolation both sides, matching
 * interp1d's fill_value=(0, Y[-1])). Returns 0.0 (not an error) if `name`
 * is not tracked. */
double cpr_nuclear_network_Y_of_t(const CPRNuclearNetwork *nn, const char *name, double t);

/* Writes a two-column "nuclide  Y" table of final abundances to
 * cfg->output_final_file (port of _write_final_result). Returns 0 on
 * success, nonzero with *errmsg set (caller frees) on a file-write
 * failure. */
int cpr_nuclear_network_write_final_result(const CPRNuclearNetwork *nn, char **errmsg);

/* Samples the unified time-evolution schema (columns
 * t_s/a/T_gamma_MeV/T_nue_MeV/T_numu_MeV/T_nutau_MeV/Y_<nuclide>) at
 * `n_points` log-spaced rows in cosmic time between T_start_cosmo and
 * nn->t_end (mirrors Python's _write_time_evolution's t_out grid). Writes
 * into caller-allocated buffers, all of length n_points except `Y_out`
 * (length n_points * nn->n_species, row-major, one row per time step in
 * nn->abundance_names order). `a_out`/`Tnue_out`/`Tnumu_out`/`Tnutau_out`
 * are filled with NaN wherever the active background has no scale-factor/
 * neutrino-sector tracking, exactly like Python's EvolutionResult does for
 * a minimal/custom background. Shared by
 * cpr_nuclear_network_write_time_evolution (TSV) and the Python-extension
 * bridge (_wrapper.c, in-memory output_time_evolution=True) so the two
 * never drift. */
void cpr_nuclear_network_sample_time_evolution(const CPRNuclearNetwork *nn, int n_points,
                                                  double *t_out, double *T_out, double *a_out,
                                                  double *Tnue_out, double *Tnumu_out,
                                                  double *Tnutau_out, double *Y_out);

/* Writes the unified time-evolution TSV (header-compatible
 * with primat.evolution.dump_evolution's output) to
 * cfg->output_file, via cpr_nuclear_network_sample_time_evolution.
 * `n_points` is the number of log-spaced output rows (mirrors
 * cfg->output_n_points). When cfg->output_rates_time_evolution and the
 * network is small/small_parthenope, the optional per-reaction forward-rate
 * columns (see cpr_nuclear_network_rate_columns) are appended after the
 * Y_<nuclide> block, matching the Python backend. A NULL/empty
 * cfg->output_file is the in-memory-only escape hatch (mirrors Python's
 * output_file=None): no-op, returns 0. Returns 0 on success, nonzero with
 * *errmsg set (caller frees) on a file-write failure. */
int cpr_nuclear_network_write_time_evolution(const CPRNuclearNetwork *nn, int n_points,
                                                char **errmsg);

/* Per-reaction forward-rate columns for output_rates_time_evolution (mirrors
 * Python's EvolutionResult.rates / NuclearNetwork._write_time_evolution).
 * Returns the number of columns: 0 unless cfg->output_rates_time_evolution,
 * else one per LT reaction *except* n<->p, which has no rate table -- 12 for
 * small/small_parthenope, 67 for large+amax=8, 428 for full large. When the
 * return is >0 and `out_names`!=NULL, fills out_names[0..count) with
 * "<reaction>_frwrd" lexicographically sorted -- the IDENTICAL names and order
 * the Python backend produces via sorted().
 * `out_names` must hold at least (n LT reactions - 1) entries; pass NULL to
 * query only the count. */
size_t cpr_nuclear_network_rate_columns(const CPRNuclearNetwork *nn,
                                          char (*out_names)[64]);

/* Samples the per-reaction forward rates (columns from
 * cpr_nuclear_network_rate_columns, same sorted order) at the `n_points`
 * photon temperatures `T_MeV` [MeV] -- typically the T_gamma_MeV column of a
 * prior cpr_nuclear_network_sample_time_evolution call. `rates_out` is
 * caller-allocated, length n_points * n_cols, row-major (one row per time
 * step, one column per reaction). Each value is the active forward reaction
 * rate linearly interpolated on the master T9 grid, exactly as
 * cpr_network_fill_buffer / Python's <rxn>_frwrd do. No-op if there are no
 * rate columns. */
void cpr_nuclear_network_sample_rates(const CPRNuclearNetwork *nn,
                                        const double *T_MeV, int n_points,
                                        double *rates_out);

/* Decay-Time (DT) era: propagate the end-of-LT abundances past T_end under
 * the constant radioactive-decay matrix D (Y(t) = exp(D*(t-t_end)) Y0), via
 * scaling-and-squaring Pade-13 matrix exponentiation -- port of
 * nuclear_network.py's _build_decay_matrix/_integrate_decay_era/
 * _write_decay_evolution. No-op (returns 0) unless
 * cfg->decay_era && the `large` network && cfg->output_decay_evolution: the
 * era changes no BBN observable (Y_final is the end-of-LT state on both
 * backends), so its only effect is the output_decay_evolution TSV (columns
 * t, then Y<nuclide> per species; np.savetxt-compatible %.18e). Returns 0 on
 * success, nonzero with *errmsg set (caller frees) on OOM / file-write /
 * singular-Pade failure. */
int cpr_nuclear_network_decay_era(const CPRNuclearNetwork *nn, char **errmsg);

#endif /* CPRIMAT_NUCLEAR_NETWORK_H */
