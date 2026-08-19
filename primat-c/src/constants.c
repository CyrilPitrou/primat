#include "constants.h"
#include <math.h>

/* Riemann zeta(3) (Apery's constant), needed by cpr_n0CMB() below. Python
 * gets this from scipy.special.zeta(3); libm has no zeta function, so the
 * literal (17 significant digits, well beyond double precision) is the
 * simplest faithful port. */
#define ZETA3 1.2020569031595942854
#ifndef M_PI
#  define M_PI  3.141592653589793238462643383279502884
#endif

/* Every field is a literal, so the whole table is initialised by the loader
 * before main() runs and is never written afterwards. That is what makes it
 * safe to read from several threads at once: cpr_config_init_defaults copies
 * it into each run's own cfg->consts, and MC worker threads build their
 * configs concurrently. */
const CPRConstants g_const = {
    /* CGS base units, natural-units convention */
    .Kelvin = 1., .second = 1., .cm = 1., .gram = 1.,

    .kB     = 1.380649e-16,
    .clight = 2.99792458e+10,
    .hbar   = 6.62607015 / (2. * M_PI) * 1e-27,
    .Mpc    = 3.08567758149e+24,
    .MeV    = 1.602176634e-6,
    .keV    = 1.602176634e-9,

    .alphaem = 1. / 137.035999084,
    .GF      = 1.1663787e-5 * 1.e-6,
    .mZ      = 91.1876e3,

    .me = 0.51099895,
    .mn = 939.56542052,
    .mp = 938.27208816,

    .T0CMB   = 2.7255,
    .Neff_SM = 3.044,

    .gA        = 1.2756,
    .kappa_p   = 2.79284734463 - 1.,
    .kappa_n   = -1.91304273,
    .Vud       = 0.9738,
    .radproton = 0.8409e-13,

    .ma        = 931.494061,
    .He4Overma = 4.0026032541,
    .HOverma   = 1.00782503223,
};

double cpr_MeV_to_Kelvin(void) { return g_const.MeV / g_const.kB; }
double cpr_MeV_to_secm1(void)  { return g_const.MeV / g_const.hbar; }
double cpr_MeV_to_g(void)      { return g_const.MeV / (g_const.clight * g_const.clight); }
double cpr_MeV_to_cmm1(void)   { return g_const.MeV / (g_const.hbar * g_const.clight); }

double cpr_MeV4_to_gcmm3(void)
{
    double cmm1 = cpr_MeV_to_cmm1();
    return cpr_MeV_to_g() * cmm1 * cmm1 * cmm1;
}

double cpr_T_start(void) { return 10.0 * cpr_MeV_to_Kelvin(); }
double cpr_T_weak(void)  { return 1.0 * cpr_MeV_to_Kelvin(); }
double cpr_T_nucl(void)  { return 0.11 * cpr_MeV_to_Kelvin(); }

double cpr_sW2(const CPRConstants *c)
{
    /* On-shell relation: sin^2(theta_W) from GF, mZ, alphaem. */
    return 0.5 * (1. - sqrt(1. - 2. * sqrt(2.) * M_PI * c->alphaem
                             / (c->GF * c->mZ * c->mZ)));
}

double cpr_geL(const CPRConstants *c) { return 0.5 + cpr_sW2(c); }
double cpr_geR(const CPRConstants *c) { return cpr_sW2(c); }
double cpr_gmuL(const CPRConstants *c) { return -0.5 + cpr_sW2(c); }
double cpr_gmuR(const CPRConstants *c) { return cpr_sW2(c); }
double cpr_deltakappa(const CPRConstants *c) { return c->kappa_p - c->kappa_n; }

double cpr_s0bar(void)
{
    /* Relativistic boson gas, g=2 (photon): s_gamma = (4 pi^2/45) T^3
     * (Phys. Rep. Eq. 24). */
    return 4. * M_PI * M_PI / 45.;
}

double cpr_n0CMB(const CPRConstants *c)
{
    /* n_gamma = (2 zeta(3)/pi^2) T^3 for a bosonic gas with g=2 (photon). */
    double t = c->T0CMB / cpr_MeV_to_Kelvin();
    return (2. * ZETA3) / (M_PI * M_PI) * t * t * t;
}

double cpr_mB(const CPRConstants *c)
{
    /* Mean baryon mass [MeV] for a 24.7% He4 mass-fraction mixture with H. */
    const double percentHe = 24.7 / 100.;
    return ((1. - percentHe) * c->HOverma
            + percentHe * c->He4Overma / 4.) * c->ma;
}

double cpr_maOvermB(const CPRConstants *c) { return c->ma / cpr_mB(c); }

double cpr_HubbleOverh(void)
{
    /* 100 km/s/Mpc converted to natural (MeV) units via the cm/s/Mpc chain. */
    return (100. * (1.e+5 * g_const.cm * cpr_MeV_to_cmm1()))
           / (g_const.second * cpr_MeV_to_secm1())
           / (g_const.Mpc * cpr_MeV_to_cmm1());
}
