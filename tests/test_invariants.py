"""
Physical invariants: what any correct BBN answer must satisfy, whatever the
code computes.

The rest of the suite pins today's numbers as today's numbers. These tests
instead assert what the physics *requires* — conservation laws, monotonicity,
analytic limits — so that a change moving both backends the same wrong way is
still caught. Each bound below is a physical statement with slack, not a
measurement: none of them is a number to be "updated" when a result moves.

Where the invariant can be checked during the integration it is, at every
accepted step of every era, because a violation that cancels by the end is
still a bug.
"""
import numpy as np
import pytest
from scipy.special import zeta

from primat.config import DEFAULT_PARAMS

# Baryon number, positivity and charge are exact statements; the bounds here
# are solver noise floors with two orders of magnitude of slack.
BARYON_TOL = 1e-10
# An abundance is never negative. The floor is the LT integration's own
# absolute tolerance: below it the solver makes no claim about the value at
# all, so an excursion smaller than that is noise, while a real sign error
# lands orders of magnitude above it. Deliberately not a measured number --
# the previous -1e-40 was one, and it failed on scipy 1.11 (where BDF
# undershoots to -2.7e-31 on the large network) while passing on scipy 1.18.
NEGATIVE_FLOOR = -DEFAULT_PARAMS["atol_large_LT"]


def _AZ(cfg, names):
    """Mass number and charge of each nuclide, as abundance-vector arrays."""
    A = np.array([sum(cfg.Nuclides[s]) for s in names], float)
    Z = np.array([cfg.Nuclides[s][1] for s in names], float)
    return A, Z


def _solve_capturing_eras(params):
    """Run a Python-backend solve, returning the three eras' solver results.

    ``solve_ivp`` is called without ``t_eval``, so each returned ``sol.t`` is
    the era's accepted-step sequence — which is what lets the conservation
    tests look inside the integration instead of only at its endpoint.
    """
    from primat.main import PRIMAT
    from primat import nuclear_network as NN

    store = {}
    originals = {n: getattr(NN.NuclearNetwork, n)
                 for n in ("_solve_HT", "_solve_MT", "_solve_LT")}

    def wrap(name):
        original = originals[name]

        def wrapped(self, *a, **k):
            out = original(self, *a, **k)
            store[name[-2:]] = out[0]
            store["net"] = self
            return out
        return wrapped

    for name in originals:
        setattr(NN.NuclearNetwork, name, wrap(name))
    try:
        run = PRIMAT(params=dict(params))
        run.solve(progress=False)
    finally:
        for name, original in originals.items():
            setattr(NN.NuclearNetwork, name, original)
    store["run"] = run
    return store


# ---------------------------------------------------------------------------
# Conservation and positivity
# ---------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.parametrize("network", ["small", "large"])
def test_baryon_number_conserved_at_every_accepted_step(network):
    """Nucleons are conserved: sum_s A_s Y_s = 1 throughout every era.

    Checked at each accepted step of HT, MT and LT rather than only at the end,
    so a violation that cancels before the last step still fails.
    """
    store = _solve_capturing_eras({"network": network})
    net = store["net"]
    cfg = net.cfg
    eras = (("HT", ["n", "p"]),
            ("MT", net.nucl._mt_net.species),
            ("LT", net.nucl.species_large))
    for era, names in eras:
        sol = store[era]
        A, _ = _AZ(cfg, names)
        baryons = A @ sol.y
        worst = float(np.abs(baryons - 1.0).max())
        assert worst < BARYON_TOL, (
            f"{network}/{era}: baryon number departs from 1 by {worst:.3e} "
            f"over {len(sol.t)} accepted steps")


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.parametrize("network", ["small", "large"])
def test_abundances_never_go_negative(network):
    """An abundance is a number of nuclei per baryon: it cannot be negative.

    Checked at every accepted step of every era, before the final clamp in
    ``_solve_LT`` hides any excursion.
    """
    store = _solve_capturing_eras({"network": network})
    net = store["net"]
    for era in ("HT", "MT", "LT"):
        worst = float(store[era].y.min())
        assert worst > NEGATIVE_FLOOR, (
            f"{network}/{era}: abundance reached {worst:.3e}")


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.backend
@pytest.mark.parametrize("backend", ["python", "c"])
def test_evolution_output_conserves_baryons_and_stays_positive(backend):
    """The published time evolution obeys the same two laws as the solver.

    Rows before the network starts integrating are exactly zero by design (see
    ``NuclearNetwork._write_time_evolution``); they are excluded, and every
    populated row must carry a full baryon budget and no negative abundance.
    """
    from primat.backend import run_bbn, HAS_C_BACKEND
    from primat.config import PRIMATConfig
    if backend == "c" and not HAS_C_BACKEND:
        pytest.skip("C extension not built")
    ev = run_bbn({"network": "large", "output_time_evolution": True,
                  "output_file": None}, force_backend=backend)["evolution"]
    cfg = PRIMATConfig()
    names = list(ev.Y.keys())
    A = np.array([sum(cfg.Nuclides[s]) for s in names], float)
    Y = np.array([ev.Y[s] for s in names])
    baryons = A @ Y
    populated = baryons > 0.5
    assert populated.sum() > 100
    assert float(np.abs(baryons[populated] - 1.0).max()) < BARYON_TOL
    assert float(Y[:, populated].min()) > NEGATIVE_FLOOR


@pytest.mark.slow
@pytest.mark.solve
def test_decay_era_conserves_baryons_and_stays_positive():
    """Radioactive decay moves nucleons between nuclides, never destroys them.

    The decay matrix must annihilate the mass-number vector exactly (A . D = 0),
    which makes the whole matrix exponential baryon-conserving by construction.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "large", "decay_era": True})
    run.solve(progress=False)
    nn = run.nuclear
    cfg = run.cfg
    names = nn.abundance_names
    A = np.array([sum(cfg.Nuclides[s]) for s in names], float)
    D = nn._build_decay_matrix(nn.nucl._lt_net)
    assert float(np.abs(A @ D).max()) < 1e-12 * float(np.abs(D).max())

    Y0 = np.array([nn.Y_final.get(s, 0.0) for s in names])
    t_end = nn._lt_t_end_s()
    t_DT = t_end + np.logspace(0.0, np.log10(cfg.t_decay_end), cfg.decay_n_points)
    Y = nn._integrate_decay_era(D, Y0, t_end, t_DT)
    assert float(np.abs(Y @ A - 1.0).max()) < BARYON_TOL
    assert float(Y.min()) > NEGATIVE_FLOOR


# ---------------------------------------------------------------------------
# Monotonicity and direction
# ---------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.solve
def test_background_runs_one_way_in_time_temperature_and_scale_factor():
    """The universe expands and cools: t and a increase, T_gamma decreases."""
    from primat.main import PRIMAT
    bg = PRIMAT(params={"network": "small"}).background
    assert np.all(np.diff(bg.t_vec) > 0)
    assert np.all(np.diff(bg.Tg_vec) < 0)
    assert np.all(np.diff(bg.a_of_T(bg.Tg_vec)) > 0)


@pytest.mark.slow
@pytest.mark.solve
def test_neutrino_temperatures_never_exceed_the_photon_temperature():
    """Photons are heated by e+e- annihilation, neutrinos (mostly) are not.

    So T_nu <= T_gamma for all three flavours at all times, and the electron
    flavour — which keeps exchanging energy with the plasma longest — stays the
    warmest of the three.
    """
    from primat.main import PRIMAT
    bg = PRIMAT(params={"network": "small"}).background
    Tg = bg.Tg_vec
    for flavour, Tnu in (("nue", bg.Tnue_vec), ("numu", bg.Tnumu_vec),
                         ("nutau", bg.Tnutau_vec)):
        assert float((Tnu / Tg).max()) < 1.0 + 1e-9, flavour
    assert float(((bg.Tnumu_vec - bg.Tnue_vec) / Tg).max()) < 1e-9
    # Well after decoupling the three are cooler than the photons by the
    # familiar (4/11)^(1/3) ~ 0.714, corrected upward by incomplete decoupling.
    assert 0.71 < float(bg.Tnue_vec[-1] / Tg[-1]) < 0.72


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.backend
@pytest.mark.parametrize("backend", ["python", "c"])
def test_published_evolution_never_reports_a_neutrino_hotter_than_the_photons(backend):
    """The same law as above, on the file a user actually reads.

    The output grid is not the solver's, so the temperature columns are
    interpolated onto it — and if the neutrino and photon columns were put on
    different interpolation schemes their ratio would exceed 1 between nodes,
    where the physics says the two are equal. Only the shipped table's own
    rounding (a few parts in 1e12) is allowed through.
    """
    from primat.backend import run_bbn, HAS_C_BACKEND
    if backend == "c" and not HAS_C_BACKEND:
        pytest.skip("C extension not built")
    ev = run_bbn({"network": "small", "output_time_evolution": True,
                  "output_file": None}, force_backend=backend)["evolution"]
    Tg = np.asarray(ev.T_gamma)
    for flavour in ("e", "mu", "tau"):
        worst = float((np.asarray(ev.T_nu[flavour]) / Tg).max())
        assert worst < 1.0 + 1e-9, f"{flavour}: T_nu/T_gamma reached {worst!r}"


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.parametrize("rate_grid_T9_max", [10.0, 2.0])
def test_reaction_rates_are_never_negative_outside_the_master_grid(rate_grid_T9_max):
    """A reaction rate is non-negative wherever the solver asks for one.

    The rate buffer interpolates linearly on the master T9 grid and continues
    the end slope outside it, so a temperature far enough beyond the grid can
    carry the extrapolation straight through zero. The MT era starts above the
    grid's top by default, and the grid's span is configurable, so this is
    checked at the era boundary for a grid deliberately too short to cover it.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "small",
                         "rate_grid_T9_max": rate_grid_T9_max})
    cfg = run.cfg
    for net in (run.nucl._mt_net, run.nucl._lt_net):
        buf = np.array(net.fill_buffer(cfg.T_weak, run.background.weak_nTOp_frwrd,
                                       run.background.weak_nTOp_bkwrd, clamp=False))
        assert buf[2::2].min() >= 0.0, "forward rate"
        assert buf[3::2].min() >= 0.0, "reverse rate"


@pytest.mark.slow
@pytest.mark.solve
def test_comoving_entropy_is_conserved_when_nothing_heats_the_neutrinos():
    """With the heating switched off, the a(T) ODE *is* entropy conservation.

    ``incomplete_decoupling=False`` sets the NEVO heating N to zero, and the
    entropy-conservation right-hand side then integrates to s a^3 = const
    exactly. The spread of that product is therefore a direct readout of the
    background ODE's own error, independent of any tolerance it was asked for.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "small", "incomplete_decoupling": False,
                         "QED_corrections": False, "spectral_distortions": False})
    T = run.background.Tg_vec
    S = np.array([run.plasma.spl(t) for t in T]) * run.background.a_of_T(T) ** 3
    assert float(S.max() / S.min() - 1.0) < 1e-4


# ---------------------------------------------------------------------------
# Analytic limits
# ---------------------------------------------------------------------------
@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.backend
@pytest.mark.parametrize("backend", ["python", "c"])
def test_Neff_is_exactly_three_in_the_instantaneous_no_QED_limit(backend):
    """Neff counts neutrino species: with instantaneous decoupling and no QED
    correction to the plasma, three flavours must give exactly 3."""
    from primat.backend import run_bbn, HAS_C_BACKEND
    if backend == "c" and not HAS_C_BACKEND:
        pytest.skip("C extension not built")
    res = run_bbn({"network": "small", "incomplete_decoupling": False,
                   "QED_corrections": False, "spectral_distortions": False},
                  force_backend=backend)
    assert abs(res["Neff"] - 3.0) < 1e-12


@pytest.mark.slow
@pytest.mark.solve
def test_weak_rate_tends_to_free_neutron_decay_at_low_temperature():
    """With no plasma left to interact with, n -> p is just free decay.

    Gamma(n->p) must therefore approach 1/tau_n from below as T -> 0.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "small"})
    cfg = run.cfg
    for T_MeV, tol in ((1e-2, 1e-3), (1e-3, 1e-4), (2e-4, 1e-5)):
        rate = float(run.background.weak_nTOp_frwrd(T_MeV * cfg.MeV_to_Kelvin))
        assert abs(rate * cfg.tau_n - 1.0) < tol, f"T = {T_MeV} MeV"


@pytest.mark.slow
@pytest.mark.solve
def test_weak_rates_obey_detailed_balance_above_neutrino_decoupling():
    """Above decoupling the neutrinos share the photon temperature, so the two
    n<->p rates must sit in thermal equilibrium with each other.

    In the Born approximation that is Gamma(p->n)/Gamma(n->p) = exp(-Q/T). With
    the finite-nucleon-mass correction on (the shipped default) the nucleon
    recoil adds the non-relativistic phase-space factor (mn/mp)^(3/2), which
    the check below divides out.
    """
    from primat.main import PRIMAT
    born = PRIMAT(params={"network": "small", "radiative_corrections": False,
                          "finite_mass_corrections": False,
                          "thermal_corrections": False,
                          "spectral_distortions": False})
    cfg = born.cfg
    Q = cfg.mn - cfg.mp
    for T_MeV, tol in ((10.0, 1e-6), (5.0, 1e-5)):
        T_K = T_MeV * cfg.MeV_to_Kelvin
        ratio = (float(born.background.weak_nTOp_bkwrd(T_K))
                 / float(born.background.weak_nTOp_frwrd(T_K)))
        assert abs(ratio / np.exp(-Q / T_MeV) - 1.0) < tol, f"Born, T = {T_MeV}"

    default = PRIMAT(params={"network": "small"})
    recoil = (cfg.mn / cfg.mp) ** 1.5
    for T_MeV in (5.0, 3.0, 2.0, 1.5):
        T_K = T_MeV * cfg.MeV_to_Kelvin
        ratio = (float(default.background.weak_nTOp_bkwrd(T_K))
                 / float(default.background.weak_nTOp_frwrd(T_K)))
        assert abs(ratio / (recoil * np.exp(-Q / T_MeV)) - 1.0) < 1e-3, \
            f"finite mass, T = {T_MeV}"


@pytest.mark.slow
@pytest.mark.solve
@pytest.mark.parametrize("network", ["small", "large"])
def test_saha_equilibrium_nulls_every_thermonuclear_reaction(network):
    """Two independent encodings of nuclear equilibrium must agree.

    Nuclear Statistical Equilibrium abundances come from mass excesses and
    spins (``nuclides.csv``, via ``_saha_YA``); the reverse reaction rates come
    from the detailed-balance triples (``detailed_balance.csv``). Both describe
    the same physics, so at NSE abundances every thermonuclear reaction's
    forward and backward fluxes must cancel. A wrong Q value or spin on any
    nuclide breaks this for that reaction alone.
    """
    from primat.main import PRIMAT
    from primat.network_builder import compile_network
    run = PRIMAT(params={"network": network})
    cfg = run.cfg
    nn = run.nuclear
    for net in (run.nucl._mt_net, run.nucl._lt_net):
        comp = compile_network(net.network, net.species)
        for T_MeV in (1.0, 0.5):
            T_K = T_MeV * cfg.MeV_to_Kelvin
            rho = run.background.rhoB_BBN(run.background.t_of_T(T_MeV))
            n_B = rho / (cfg.ma * cfg.MeV4_to_gcmm3)
            eta_b = n_B / ((2.0 * zeta(3) / np.pi ** 2) * T_MeV ** 3)
            Y = np.array([0.25 if s == "n" else 0.75 if s == "p"
                          else nn._saha_YA(s, 0.25, 0.75, T_K, eta_b)
                          for s in net.species])
            buf = np.array(net.fill_buffer(T_K, run.background.weak_nTOp_frwrd,
                                           run.background.weak_nTOp_bkwrd,
                                           clamp=False))
            fwd, bwd = _reaction_fluxes(comp, Y, rho, buf)
            checkable = np.isfinite(fwd) & np.isfinite(bwd) & (fwd > 0) & (bwd > 0)
            for i in net.weak_indices:      # weak and beta channels are not in NSE
                checkable[i] = False
            assert checkable.sum() > 5
            worst = float(np.abs(1.0 - bwd[checkable] / fwd[checkable]).max())
            assert worst < 1e-4, f"{network}, T = {T_MeV} MeV: {worst:.3e}"


def _reaction_fluxes(comp, Y, rho, rates):
    """Forward and backward mass-action flux of each reaction, one by one.

    Mirrors ``network_builder._rhs_kernel``'s inner loop but keeps the two
    fluxes separate instead of accumulating their difference, which is what
    lets a caller ask whether they cancel.
    """
    n_rx = comp.ri_len.shape[0]
    fwd = np.zeros(n_rx)
    bwd = np.zeros(n_rx)
    for i in range(n_rx):
        f = rates[2 * i] * rho ** comp.Rm1[i] * comp.invsr[i]
        for k in range(comp.ri_len[i]):
            f *= Y[comp.ri_idx[i, k]] ** comp.ri_pow[i, k]
        b = rates[2 * i + 1] * rho ** comp.Pm1[i] * comp.invsp[i]
        for k in range(comp.pi_len[i]):
            b *= Y[comp.pi_idx[i, k]] ** comp.pi_pow[i, k]
        fwd[i], bwd[i] = f, b
    return fwd, bwd


@pytest.mark.slow
@pytest.mark.solve
def test_deuterium_follows_its_known_power_law_in_the_baryon_density():
    """More baryons burn deuterium away faster: D/H falls as a steep power law.

    The standard result is D/H ~ (Omega_b h^2)^-1.6; the bound below is wide
    enough that only a change to the scaling itself, not to its coefficient,
    can trip it.
    """
    from primat.backend import run_bbn
    grid = np.logspace(np.log10(0.006), np.log10(0.05), 13)
    DoH = np.array([run_bbn({"network": "small", "Omegabh2": float(ob)})["DoH"]
                    for ob in grid])
    assert np.all(np.diff(DoH) < 0)
    slope = float(np.polyfit(np.log(grid), np.log(DoH), 1)[0])
    assert -1.8 < slope < -1.5


@pytest.mark.slow
@pytest.mark.solve
def test_lithium7_has_a_valley_in_the_baryon_density():
    """Li7 is made two ways — directly at low eta, as Be7 at high eta — so its
    abundance has a minimum where the two channels cross.

    The valley sits near eta10 ~ 2.5, i.e. well below the observed baryon
    density: this is the lithium problem, and its *position* is a prediction
    the network must reproduce, not only its depth.
    """
    from primat.backend import run_bbn
    grid = np.logspace(np.log10(0.006), np.log10(0.016), 11)
    Li7 = np.array([run_bbn({"network": "small", "Omegabh2": float(ob)})["Li7oH"]
                    for ob in grid])
    i = int(np.argmin(Li7))
    assert 0 < i < len(grid) - 1, "the valley is not bracketed by the scan"
    eta10 = 273.9 * grid[i]
    assert 2.0 < eta10 < 3.2, f"Li7 valley at eta10 = {eta10:.2f}"


@pytest.mark.slow
@pytest.mark.solve
def test_free_neutrons_decay_with_the_neutron_lifetime_after_BBN():
    """Once nucleosynthesis is over, the leftover free neutrons only decay.

    Their abundance must follow exp(-dt/tau_n) exactly, with tau_n the config's
    own neutron lifetime — the decay era propagates nothing else.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "large", "decay_era": True})
    run.solve(progress=False)
    nn = run.nuclear
    cfg = run.cfg
    names = nn.abundance_names
    D = nn._build_decay_matrix(nn.nucl._lt_net)
    Y0 = np.array([nn.Y_final.get(s, 0.0) for s in names])
    t_end = nn._lt_t_end_s()
    t_DT = t_end + np.logspace(0.0, np.log10(cfg.t_decay_end), cfg.decay_n_points)
    Y = nn._integrate_decay_era(D, Y0, t_end, t_DT)
    k = names.index("n")
    dt = t_DT - t_end
    predicted = Y0[k] * np.exp(-dt / cfg.tau_n)
    resolvable = predicted > 1e-300
    assert resolvable.sum() > 20
    worst = float(np.abs(Y[resolvable, k] / predicted[resolvable] - 1.0).max())
    assert worst < 1e-6


@pytest.mark.slow
@pytest.mark.solve
def test_photon_temperature_stays_positive_past_the_end_of_the_grid():
    """T(a) is a temperature everywhere, including outside its own node range.

    The scale-factor grid ends at T_end_MeV. A linear extrapolation of a
    decaying T(a) crosses zero a finite distance past the last node, so the
    tail follows the power law through the last two nodes instead. Checked out
    to 1e6 times the final scale factor, and against T ~ 1/a, which is what
    radiation domination requires there.
    """
    from primat.main import PRIMAT
    run = PRIMAT(params={"network": "small"})
    bg = run.background
    a_end = float(np.asarray(bg.a_of_T(bg.Tg_vec)).max())
    T_end = float(bg.T_of_a(a_end))
    for factor in (1.5, 2.0, 10.0, 1e3, 1e6):
        T = float(bg.T_of_a(a_end * factor))
        assert T > 0.0, f"T_of_a({factor} x a_end) = {T}"
        assert T == pytest.approx(T_end / factor, rel=1e-3)
