/* constants.h -- physical constants and unit-conversion factors.
 *
 * Direct port of primat/constants.py's `Constants` dataclass: the base fields
 * are literal PDG values (verbatim, no computation differs from Python).
 *
 * The 26 fields split by one question -- does this number have an error bar?
 * The 16 MEASURED ones (alphaem, GF, mZ, me, mn, mp, T0CMB, gA, Vud, kappa_p,
 * kappa_n, radproton, ma, He4Overma, HOverma, Neff_SM) are ordinary
 * parameters, stored per run in `cfg->consts` and settable through
 * cpr_config_set_by_name like any other key. The other 10 are exact by
 * definition and never move, so `g_const` (the defaults) is the right source
 * for them anywhere -- including inside the derived helpers below that take
 * no argument.
 *
 * Derived quantities are functions rather than dataclass `@property`s (C has
 * no lazy per-instance property mechanism, and they are cheap to recompute).
 * Those depending on a measured field take a `const CPRConstants *`; those
 * built from the exact ten take none.
 *
 * Units convention (unchanged from Python): natural units throughout
 * (Kelvin = second = cm = gram = 1); the `cpr_MeV_to_*` functions convert
 * *to* CGS only where needed.
 */
#ifndef CPRIMAT_CONSTANTS_H
#define CPRIMAT_CONSTANTS_H

typedef struct {
    /* ---- CGS base units (dimensionless by convention: natural units) ---- */
    double Kelvin, second, cm, gram;

    /* ---- Fundamental constants (PDG) ---- */
    double kB;      /* Boltzmann constant [erg/K] */
    double clight;  /* speed of light [cm/s] */
    double hbar;    /* reduced Planck constant h/2pi [erg s] */
    double Mpc;     /* megaparsec [cm] */
    double MeV;     /* 1 MeV [erg] */
    double keV;     /* 1 keV [erg] */

    /* ---- Electroweak sector (PDG) ---- */
    double alphaem; /* fine-structure constant */
    double GF;      /* Fermi constant [MeV^-2] */
    double mZ;      /* Z boson mass [MeV] */

    /* ---- Fermion masses [MeV] (PDG) ---- */
    double me, mn, mp;

    /* ---- CMB ---- */
    double T0CMB;   /* photon temperature today [K] */

    /* ---- Standard-model effective neutrino number ---- */
    double Neff_SM; /* 3 instantaneous-decoupling flavours + NEVO/QED heating corrections */

    /* ---- Weak-rate nuclear-structure constants (PDG) ---- */
    double gA;        /* nucleon axial coupling */
    double kappa_p;    /* proton anomalous magnetic moment */
    double kappa_n;    /* neutron anomalous magnetic moment */
    double Vud;        /* CKM matrix element |V_ud| */
    double radproton;  /* proton charge radius [cm] */

    /* ---- Atomic masses ---- */
    double ma;         /* 1 unified atomic mass unit [MeV] */
    double He4Overma;  /* M(He4) / u */
    double HOverma;    /* M(H) / u */
} CPRConstants;

/* The DEFAULT values, initialised statically and never written. Read it
 * directly only for the ten exact constants; a run's 16 measured ones live in
 * `cfg->consts`, which cpr_config_init_defaults seeds from here. Being
 * immutable is what lets several threads build configurations at once. */
extern const CPRConstants g_const;

/* ---- Derived from the exact ten only (never move) ---- */
double cpr_MeV_to_Kelvin(void);
double cpr_MeV_to_secm1(void);
double cpr_MeV_to_g(void);
double cpr_MeV_to_cmm1(void);
double cpr_MeV4_to_gcmm3(void);

/* ---- Fixed temperature eras: quoted in MeV, returned in Kelvin ----
 * cpr_T_start_nucl is where the *nuclear* network starts; the *background*
 * integration starts higher, at the overridable cfg->T_start_cosmo_MeV. */
double cpr_T_start_nucl(void);  /* 10 MeV */
double cpr_T_weak(void);   /* 1 MeV */
double cpr_T_nucl(void);   /* 0.11 MeV */

/* ---- Electroweak mixing angle and effective couplings (alphaem, GF, mZ) ---- */
double cpr_sW2(const CPRConstants *c);
double cpr_geL(const CPRConstants *c);
double cpr_geR(const CPRConstants *c);
double cpr_gmuL(const CPRConstants *c);
double cpr_gmuR(const CPRConstants *c);
double cpr_deltakappa(const CPRConstants *c);  /* kappa_p, kappa_n */

/* ---- High-T plasma entropy/number-density normalisations ---- */
double cpr_s0bar(void);
double cpr_n0CMB(const CPRConstants *c);  /* T0CMB */

/* ---- Mean baryon mass (H + He4 mixture): ma, He4Overma, HOverma ---- */
double cpr_mB(const CPRConstants *c);
double cpr_maOvermB(const CPRConstants *c);

/* ---- Hubble constant in natural units, per unit h ---- */
double cpr_HubbleOverh(void);

#endif /* CPRIMAT_CONSTANTS_H */
