# -*- coding: utf-8 -*-
"""
constants.py
============
Physical constants and unit-conversion factors for primat.

This module is the single source of truth for the constants used throughout
the code: PDG masses and couplings, the CGS-vs-natural-units conversion
factors, and the derived quantities (the Weinberg angle ``sW2``, the mean
baryon mass ``mB``, …) that follow from them.

The 26 fields split by one question — *does this number have an error bar?*
:data:`OVERRIDABLE_CONSTANTS` (16) are measured, so they are ``DEFAULT_PARAMS``
keys settable like any other parameter; :data:`FROZEN_CONSTANTS` (10) are
exact by the 2019 SI redefinition, an IAU definition, or the natural-units
convention, so they stay fixed and ``PRIMATConfig`` rejects an override.

``CONST`` is the frozen default instance.  A ``PRIMATConfig`` exposes its own
values as plain attributes (``cfg.me``, ``cfg.sW2``, …) and its own snapshot as
``cfg.constants``; ``CONST`` remains the right import for code with no config
in hand:

    >>> from primat.constants import CONST
    >>> CONST.me
    0.51099895
    >>> CONST.MeV_to_Kelvin   # 1 MeV in Kelvin
    11604518121.5...

Units convention
-----------------
All "CGS" quantities (``Kelvin``, ``second``, ``cm``, ``gram``) are set to
1: lengths/times/temperatures/masses are expressed in natural (MeV-based)
units throughout the code, and the ``MeV_to_*`` factors below convert *to*
CGS only where needed (e.g. for printing or comparison with CGS-valued
inputs such as ``T0CMB`` [K]).
"""

from dataclasses import dataclass
import numpy as np
from scipy.special import zeta

__all__ = ['Constants', 'CONST', 'OVERRIDABLE_CONSTANTS', 'FROZEN_CONSTANTS',
           'DERIVED_OVERRIDABLE']

# The 16 measured constants, promoted to DEFAULT_PARAMS keys by
# primat.config: each carries an experimental uncertainty, so varying it is a
# sensitivity study rather than a redefinition.
OVERRIDABLE_CONSTANTS = (
    'alphaem', 'GF', 'mZ', 'me', 'mn', 'mp', 'T0CMB', 'gA', 'Vud',
    'kappa_p', 'kappa_n', 'radproton', 'ma', 'He4Overma', 'HOverma', 'Neff_SM',
)

# The 10 that stay fixed: Kelvin/second/cm/gram are 1 by the natural-units
# convention, kB/clight/hbar/MeV/keV are exact by the 2019 SI redefinition,
# and Mpc is an IAU definition. PRIMATConfig rejects an override of any of
# them (see PRIMATConfig.validate_frozen_constants).
FROZEN_CONSTANTS = (
    'Kelvin', 'second', 'cm', 'gram',
    'kB', 'clight', 'hbar', 'MeV', 'keV', 'Mpc',
)

# The derived properties that depend on at least one OVERRIDABLE_CONSTANTS
# field, and so must be recomputed whenever a config overrides one. The rest
# (MeV_to_*, T_start/T_weak/T_nucl, s0bar, HubbleOverh, GN_*_to_*) are
# functions of frozen fields only and never move.
DERIVED_OVERRIDABLE = (
    'sW2', 'geL', 'geR', 'gmuL', 'gmuR',   # alphaem, GF, mZ
    'deltakappa',                          # kappa_p, kappa_n
    'n0CMB',                               # T0CMB
    'mB', 'maOvermB',                      # ma, He4Overma, HOverma
)


@dataclass(frozen=True)
class Constants:
    """Physical constants and unit-conversion factors (frozen dataclass).

    A ``PRIMATConfig`` carries its own instance (``cfg.constants``), built by
    replacing the :data:`OVERRIDABLE_CONSTANTS` fields with that config's
    values; :data:`CONST` is the all-defaults one.

    Grouped as:

    - CGS base units (all set to 1; natural units throughout)
    - Fundamental constants: kB, clight, hbar, Mpc, MeV, keV
    - Electroweak sector: alphaem, GF, mZ
    - Fermion masses [MeV]: me, mn, mp
    - CMB: T0CMB
    - Weak-rate nuclear-structure constants: gA, kappa_p, kappa_n,
      Vud, radproton
    - Atomic masses: ma, He4Overma, HOverma

    Each field below carries the edition it came from (CODATA/PDG/AME year, or
    "SI 2019, exact"); they are not all from one vintage.

    plus derived quantities (unit-conversion factors, fixed temperature
    eras, electroweak mixing-angle couplings, and the high-T plasma
    entropy/number-density normalisations) exposed as read-only
    properties computed from the fields above.
    """

    # ---- CGS base units (dimensionless by convention: natural units) ----
    Kelvin: float = 1.
    second: float = 1.
    cm:     float = 1.
    gram:   float = 1.

    # ------------------------------------------------------------------
    # Provenance. These values were assembled over time and are NOT all from
    # one edition; the per-field tags below record which, so a future refresh
    # can tell what must move together. All are well within primat's target
    # precision as they stand -- this is bookkeeping, not a correction list.
    # The weak-rate sector (gA, Vud, and the two anomalous moments) is the one
    # to watch: it feeds the n<->p rates directly (primat/weak_rates/).
    # ------------------------------------------------------------------

    # ---- Fundamental constants ----
    # The first four are *exact* by the 2019 SI redefinition (no uncertainty).
    kB:     float = 1.380649e-16          # Boltzmann constant [erg/K]      (SI 2019, exact)
    clight: float = 2.99792458e+10        # speed of light [cm/s]           (SI 2019, exact)
    hbar:   float = 6.62607015 / (2 * np.pi) * 1e-27  # Planck constant [erg s]  (SI 2019, exact h)
    Mpc:    float = 3.08567758149e+24     # megaparsec [cm]                 (IAU 2015)
    MeV:    float = 1.602176634e-6        # 1 MeV [erg]                     (SI 2019, exact e)
    keV:    float = 1.602176634e-9        # 1 keV [erg]                     (SI 2019, exact e)

    # ---- Electroweak sector ----
    alphaem: float = 1. / 137.035999084   # fine-structure constant         (CODATA 2018)
    GF:      float = 1.1663787e-5 * 1.e-6 # Fermi constant [MeV^-2]         (PDG 2020)
    mZ:      float = 91.1876e3            # Z boson mass [MeV]              (PDG 2020)

    # ---- Fermion masses [MeV] ----
    me: float = 0.51099895                # electron                        (CODATA 2018)
    mn: float = 939.56542052              # neutron                         (PDG/CODATA 2018)
    mp: float = 938.27208816              # proton                          (CODATA 2018)

    # ---- CMB ----
    T0CMB: float = 2.7255                 # photon temperature today [K]    (Fixsen 2009, ApJ 707, 916)

    # ---- Weak-rate nuclear-structure constants ----
    gA:        float = 1.2756              # nucleon axial coupling         (PDG 2018)
    kappa_p:   float = 2.79284734463 - 1.  # proton anomalous magnetic moment  (CODATA 2018)
    kappa_n:   float = -1.91304273         # neutron anomalous magnetic moment (CODATA 2018)
    Vud:       float = 0.9738              # CKM matrix element |V_ud|      (PDG 2018)
    radproton: float = 0.8409e-13          # proton charge radius [cm]      (CODATA 2018)

    # ---- Atomic masses ----
    ma:        float = 931.494061          # 1 unified atomic mass unit [MeV]  (CODATA 2010;
                                           #   CODATA 2018 is 931.49410242, a 4e-8 relative shift)
    He4Overma: float = 4.0026032541        # M(He4) / u                     (AME2020)
    HOverma:   float = 1.00782503223       # M(H) / u                       (AME2016)

    # ---- Standard-model effective neutrino number ----
    # The Standard Model prediction for N_eff including non-instantaneous
    # decoupling, finite-temperature QED, and flavour-oscillation effects
    # (Bennett et al. 2021, arXiv:2012.02726; de Salas & Pastor numerical
    # value to 3.044 quoted from the more recent re-evaluations of the
    # Mangano/Miele 3.046 result). Used wherever the *standard*-physics value
    # of N_eff is needed as an input, rather than the value primat itself
    # solves for via the NEVO table: numerically only in the EDE-era radiation
    # normalisation, elsewhere as the reference point of the reported
    # ``Neff = Neff_SM + DeltaNeff`` relation.
    Neff_SM:   float = 3.044

    # ------------------------------------------------------------------
    # Derived quantities (pure functions of the constants above)
    # ------------------------------------------------------------------

    @property
    def MeV_to_Kelvin(self) -> float:
        """Conversion factor: 1 MeV / kB, in Kelvin."""
        return self.MeV / self.kB

    @property
    def MeV_to_secm1(self) -> float:
        """Conversion factor: 1 MeV / hbar, in s^-1."""
        return self.MeV / self.hbar

    @property
    def MeV_to_g(self) -> float:
        """Conversion factor: 1 MeV / c^2, in g."""
        return self.MeV / self.clight**2

    @property
    def MeV_to_cmm1(self) -> float:
        """Conversion factor: 1 MeV / (hbar c), in cm^-1."""
        return self.MeV / (self.hbar * self.clight)

    @property
    def MeV4_to_gcmm3(self) -> float:
        """Conversion factor for an energy density [MeV^4] to a mass density [g/cm^3]."""
        return self.MeV_to_g * self.MeV_to_cmm1**3

    @property
    def GN_MeV2_to_SI(self) -> float:
        """Conversion factor: Newton's constant in natural units [MeV^-2] to
        SI units [m^3 kg^-1 s^-2].

        Particle physics conventionally expresses ``G`` in the natural-units
        convention (hbar = c = 1), where it has the dimension of
        [energy]^-2 (``G = 1/m_Pl^2``, with ``m_Pl`` the Planck mass in MeV).
        Restoring hbar and c, ``G`` in CGS units [cm^3 g^-1 s^-2] is
        ``G_natural[MeV^-2] * hbar[erg s] * clight[cm/s]^5 / MeV[erg]^2``
        (the erg/cm/s factors cancel exactly into CGS's
        gravitational-constant units); one more factor of ``1e-3`` converts
        cm^3 g^-1 s^-2 to m^3 kg^-1 s^-2.

        ``cfg.GN`` itself is stored and exposed in SI units throughout
        primat (``DEFAULT_PARAMS["GN"]`` in ``config.py``, the GUI's
        "Constants" panel); this factor (and its inverse,
        :attr:`GN_SI_to_MeV2`) is used internally to convert to/from the
        natural units that the Friedmann equation
        (``PRIMATConfig.Mpl``/``PRIMATConfig.rhocOverh2``) is written in.

        Example
        -------
        >>> CONST.GN_MeV2_to_SI * 6.70883e-45   # doctest: +SKIP
        6.674...e-11
        """
        return self.hbar * self.clight**5 / self.MeV**2 * 1e-3

    @property
    def GN_SI_to_MeV2(self) -> float:
        """Conversion factor: Newton's constant in SI units
        [m^3 kg^-1 s^-2] to natural units [MeV^-2]. Inverse of
        :attr:`GN_MeV2_to_SI` -- see that property for the derivation.
        """
        return 1. / self.GN_MeV2_to_SI

    # ---- Fixed temperature eras [MeV, converted to Kelvin] ----
    @property
    def T_start(self) -> float:
        return 10.0 * self.MeV_to_Kelvin

    @property
    def T_weak(self) -> float:
        return 1.0 * self.MeV_to_Kelvin

    @property
    def T_nucl(self) -> float:
        return 0.11 * self.MeV_to_Kelvin
        
    # ---- Electroweak mixing angle and effective couplings ----
    @property
    def sW2(self) -> float:
        """sin^2(theta_W), from GF, mZ, alphaem (on-shell relation)."""
        return 0.5 * (1. - np.sqrt(1. - 2.*np.sqrt(2.)*np.pi*self.alphaem
                                    / (self.GF * self.mZ**2)))

    @property
    def geL(self) -> float:
        return 0.5 + self.sW2

    @property
    def geR(self) -> float:
        return self.sW2

    @property
    def gmuL(self) -> float:
        return -0.5 + self.sW2

    @property
    def gmuR(self) -> float:
        return self.sW2

    @property
    def deltakappa(self) -> float:
        return self.kappa_p - self.kappa_n

    # ---- High-T plasma entropy/number-density normalisations ----
    @property
    def s0bar(self) -> float:
        """Dimensionless prefactor in the photon entropy density: s_gamma = s0bar T^3.

        For a relativistic boson gas with g=2 (photon) the entropy density is
            s_gamma = (2 pi^2/45) x 2 x T^3 = (4 pi^2/45) T^3  [Phys. Rep. Eq. 24].
        """
        return 4. * np.pi**2 / 45.

    @property
    def n0CMB(self) -> float:
        """Present-day CMB photon number density [MeV^3].

        n_gamma = (2 zeta(3)/pi^2) T^3 for a bosonic gas with g=2 (photon).
        """
        return (2. * zeta(3)) / np.pi**2 * (self.T0CMB / self.MeV_to_Kelvin)**3

    # ---- Mean baryon mass (H + He4 mixture) ----
    @property
    def mB(self) -> float:
        """Mean baryon mass [MeV], for a 24.7% He4 mass-fraction mixture with H."""
        percentHe = 24.7 / 100.
        return ((1. - percentHe) * self.HOverma
                + percentHe * self.He4Overma / 4.) * self.ma

    @property
    def maOvermB(self) -> float:
        return self.ma / self.mB

    # ---- Hubble constant in natural units, per unit h ----
    @property
    def HubbleOverh(self) -> float:
        """H0 / h, converted to natural (MeV) units, in MeV.

        100 km/s/Mpc converted via the cm/s/Mpc -> MeV chain.
        """
        return (100. * (1.e+5 * self.cm * self.MeV_to_cmm1)
                / (self.second * self.MeV_to_secm1)
                / (self.Mpc * self.MeV_to_cmm1))


# Single shared instance: all fields/properties above are pure constants,
# so one frozen object suffices for the whole process.
CONST = Constants()
