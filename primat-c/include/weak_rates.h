/* weak_rates.h -- n<->p weak-rate tables (port of primat/weak_rates/).
 *
 * Always-needed pieces: the non-thermal n<->p
 * rate (Born/CCR + finite-nucleon-mass + spectral-distortion corrections),
 * computed from scratch via fixed Gauss-Legendre quadrature when no cache
 * file matches the current configuration's fingerprint, or loaded directly
 * from data/cache_plasma_weak/weak/nTOp_<hash>.txt otherwise (cache.c already ports the
 * fingerprint/hash/cache-file machinery, see cache.h). The finite-
 * temperature radiative correction (CCRTh, Brown & Sawyer 2001) is S7b:
 * cpr_weak_rates_init *loads* its cache file
 * (data/cache_plasma_weak/weak/nTOp_thermal_<hash>.txt) when cfg->thermal_corrections is set
 * and a matching file exists, and otherwise recomputes it from scratch via
 * the same algorithm Python's `corrections.py` uses -- VEGAS adaptive
 * Monte Carlo (vegas.h) for the three 2D sub-integrals, deterministic 1D
 * quadrature for the one 1D sub-integral -- see weak_rates.c's CCRTh
 * section.
 *
 * The combined spectral-distortion / finite-mass terms apply only in
 * analytic-distortion mode, i.e. when cfg->analytic_distortions,
 * cfg->spectral_distortions and cfg->finite_mass_corrections are all set.
 *
 * Reference: Pitrou, Coc, Uzan & Vangioni, Phys. Rep. 2018
 * (arXiv:1806.11095), cited below as "Phys. Rep.".
 */
#ifndef CPRIMAT_WEAK_RATES_H
#define CPRIMAT_WEAK_RATES_H

#include "config.h"
#include "neutrino_history.h"
#include "spline.h"
#include <stddef.h>

/* Fermi-Coulomb factor F(b) (Phys. Rep. S III.D); see corrections.FermiCoulomb. */
double cpr_fermi_coulomb(double b, const CPRConfig *cfg);

/* Resummed T=0 radiative correction R(b,y,en) (Phys. Rep. Eq. 101-105);
 * see corrections.RadCorrResum. b = v/c, y = E_nu/me, en = E_e/me. */
double cpr_rad_corr_resum(double b, double y, double en, const CPRConfig *cfg);

/* Neutron-decay phase-space integral Fn (Phys. Rep. Eq. 89-91); normalises
 * K = 1/(tau_n * Fn). See corrections.ComputeFn. */
double cpr_compute_fn(const CPRConfig *cfg);

/* The non-thermal rate table plus, when available, the separately-cached
 * thermal (CCRTh) correction -- mirrors RecomputeWeakRates's two-piece
 * return value, but stored as raw arrays plus fitted interpolants rather than
 * closures. Both tables are in units of 1/tau_n (the caller multiplies by
 * 1/cfg->tau_n, or by the corresponding K, to get the physical rate in s^-1).
 *
 * The nonthermal forward/backward rates are evaluated through log10-log10
 * not-a-knot cubic splines (frwrd_sp/bkwrd_sp), matching Python's
 * _weak_rate_loglog_interp and the nuclear rate tables' cpr_resample_rate_table
 * -- see cpr_weak_rate_nTOp. This replaced the earlier linear-space local
 * 3-point quadratic, which differed from Python's interpolant by ~1e-4 in the
 * n/p freeze-out window and was the dominant C-vs-Python D/H parity gap
 * (tests/test_backend_parity.py). The backward rate is exp(-Q/T)-suppressed to
 * exact 0 in a low-T prefix where log10 is undefined, so its interpolant
 * covers only the contiguous positive suffix and the rate is pinned to 0 below
 * that suffix's lowest T (frwrd is positive throughout, so its suffix is the
 * full grid). */

/* One nonthermal rate's log10-log10 interpolant over its positive suffix,
 * mirroring Python's _weak_rate_loglog_interp: a not-a-knot cubic spline of
 * log10(Gamma) vs log10(T) when the suffix has >= 4 knots, else a log-log
 * linear interpolant (the same <4-knot fallback Python uses -- unreachable for
 * the default grids, kept only for a pathologically short custom grid). */
typedef struct {
    double *logT, *logR; /* log10(T[K]), log10(Gamma[1/tau_n]) of the suffix */
    size_t n;            /* suffix length */
    double Tmin;         /* below this raw T[K], rate == 0; -INFINITY when the
                          * whole grid is positive (forward rate: extrapolate) */
    CPRCubicSpline sp;   /* valid iff cubic != 0 */
    int cubic;           /* 1 iff the not-a-knot spline was built (n >= 4) */
} CPRWeakInterp;

typedef struct {
    double *T, *frwrd, *bkwrd; /* nonthermal table: T[K], Gamma_nTOp, Gamma_pTOn */
    size_t n;
    CPRWeakInterp frwrd_i, bkwrd_i; /* log10-log10 interpolants, see above */
    double *T_th, *Lnth, *Lpth; /* thermal correction table, only if has_thermal */
    size_t n_th;
    int has_thermal;
    /* CCRTh interpolants: not-a-knot cubic splines of L vs T in LINEAR space
     * (the n->p correction changes sign, so the log-log scheme the nonthermal
     * rates use is unavailable here), mirroring Python's
     * corrections._L_CCRTh_interpolants. Valid iff th_cubic; below 4 knots
     * both backends fall back to linear interpolation of the same nodes. */
    CPRCubicSpline Lnth_sp, Lpth_sp;
    int th_cubic;
} CPRWeakRates;

/* Builds the n<->p weak-rate tables for the given background, mirroring
 * weak_rates.RecomputeWeakRates([Tg_vec, Tnu_vec], cfg, dFDneu_func=...).
 *
 * Tg_MeV/Tnu_MeV (length n_bg): photon and (electron-flavour) neutrino
 * temperatures in MeV, e.g. StandardBackground._setup_background_and_cosmo's Tg_vec/
 * Tnue_vec -- despite ComputeWeakRates's Python docstring saying "Kelvin",
 * background.py actually passes MeV arrays (_build_rate_context converts
 * via cfg.MeV_to_Kelvin); confirmed by reading the caller in background.py.
 * Used only to build the T_nu(T_gamma)/T_gamma ratio interpolant feeding
 * the rate integrands -- not stored.
 *
 * nh: neutrino history (cpr_neutrino_history_init), supplies the NEVO
 * spectral-distortion correction dFDneu when cfg->spectral_distortions.
 *
 * On a fingerprint cache hit (cfg->weak_rate_cache and a matching
 * data/cache_plasma_weak/weak/nTOp_<hash>.txt exists), the nonthermal table is loaded
 * directly (no integration). Otherwise it is computed via the
 * Gauss-Legendre rate integrals (Born/CCR/FM/SD) and, if cfg->save_nTOp,
 * written to that cache file. The thermal correction is loaded from
 * data/cache_plasma_weak/weak/nTOp_thermal_<hash>.txt when cfg->thermal_corrections is set
 * and that file exists (`has_thermal` is then 1), and otherwise computed
 * from scratch by L_CCRTh_compute.
 *
 * Returns 0 on success (caller must cpr_weak_rates_free), nonzero with
 * *errmsg set (caller frees) otherwise. */
int cpr_weak_rates_init(CPRWeakRates *wr, const double *Tg_MeV, const double *Tnu_MeV,
                         size_t n_bg, const CPRConfig *cfg, const CPRNeutrinoHistory *nh,
                         char **errmsg);

void cpr_weak_rates_free(CPRWeakRates *wr);

/* Reject a tabulated n<->p rate table containing NaN or infinity. `source` is
 * "computed" or a cache-file path, quoted back in the message, which mirrors
 * primat/weak_rates/api.py's validate_weak_rates_finite word for word.
 * Returns 0 if every entry is finite, 1 with *errmsg set (caller frees). */
int cpr_validate_weak_rates_finite(const CPRConfig *cfg, const double *T,
                                    const double *frwrd, const double *bkwrd,
                                    size_t n, const char *source, char **errmsg);

/* Gamma_{n->p}(T)/Gamma_{p->n}(T) in units of 1/tau_n, the sum of the
 * nonthermal rate (log10-log10 cubic-interpolated, see CPRWeakInterp) and,
 * when present, the thermal correction table -- mirrors RecomputeWeakRates's
 * returned closures. */
double cpr_weak_rate_nTOp(const CPRWeakRates *wr, double T_K);
double cpr_weak_rate_pTOn(const CPRWeakRates *wr, double T_K);

#endif /* CPRIMAT_WEAK_RATES_H */
