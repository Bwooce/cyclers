"""#905: the Sun-forced periodic-orbit driver (``search/sun_forced_905.py``) and its published
positive controls.

The controls run the driver from a THREE-BODY parent (Sun off) through the Rhouma-Chicone
screens, the Melnikov zeros in the Sun phase and the continuation in eps, and compare the
orbits it finds at eps = 1 with the PRINTED states of the papers. They are not closure checks
of the printed states (those are in tests/core); the driver has to find the orbits.

* Oshima (2022), ASR 70:1325, Tables 1-5 (p1326, p1332): the four 1:1 synodic resonant spatial
  retrograde families in the bicircular model, printed to ten digits at the y = 0 crossings with
  the Sun angle. The paper prints no three-body parent; the parent is obtained by dropping the
  Sun and correcting Table 4 row 1 as a symmetric three-body orbit at period Tg (it moves 0.022
  from the printed state, so a match to 1e-9 is not inherited from the seed).
* Leiva & Briozzo (2008), CMDA 101:225, Table 1 (p234, three-body parents on the section
  x = L1) and Table 2 (p238, the quasi-bicircular periodic orbits at clock t_i, frame rotated by
  pi). Their homotopy H_RTBP + eps (H_QBCP - H_RTBP) is the driver's ``qbcp`` homotopy; mass
  ratio 0.0121505482 (the value their section abscissa implies; see
  tests/core/test_leiva_briozzo_2008_tables.py).
  The printed Table 2 states themselves close only to 5e-6 .. 6e-5 in this module, so that is
  the yardstick for agreement, not their nine printed digits.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.search import sun_forced_905 as sf

FloatArray = np.ndarray

sf.set_threads(4)

# --- the fields against the core models --------------------------------------------------------

_RNG = np.random.default_rng(905)
_STATES = [
    np.array([0.5, 0.1, 0.05, 0.02, 0.3, -0.01]) + 0.2 * _RNG.normal(size=6) for _ in range(4)
]


@pytest.mark.parametrize("state", _STATES)
def test_bicircular_field_matches_core_model_and_three_body_limit(state: FloatArray) -> None:
    model = sf.bcr4bp_model()
    system = bcr4bp.andreu_default()
    for t in (0.0, 1.7, 5.3):
        f1 = sf.vector_field(model, 1.0, t, state)
        assert np.max(np.abs(f1 - bcr4bp.bcr4bp_eom(t, state, system))) < 1e-14
        f0 = sf.vector_field(model, 0.0, t, state)
        assert np.max(np.abs(f0 - cr3bp.cr3bp_eom(t, state, model.mu))) < 1e-14


@pytest.mark.parametrize("state", _STATES)
def test_coherent_field_matches_core_model_and_three_body_limit(state: FloatArray) -> None:
    model = sf.qbcp_model()
    system = qbcp.qbcp_default()
    for t in (0.0, 1.7, 5.3):
        f1 = sf.vector_field(model, 1.0, t, state)
        assert np.max(np.abs(f1 - qbcp.qbcp_eom(t, state, system))) < 1e-13
        assert (
            np.max(
                np.abs(sf.pv_to_model(model, state, t, 1.0) - qbcp.state_pv_to_pm(state, t, system))
            )
            < 1e-15
        )
        # eps = 0 is the three-body problem in canonical variables px = vx - y, py = vy + x.
        can = sf.pv_to_model(model, state, t, 0.0)
        f0 = sf.vector_field(model, 0.0, t, can)
        ref = cr3bp.cr3bp_eom(t, state, model.mu)
        dv = np.array([f0[3] + f0[1], f0[4] - f0[0], f0[5]])
        assert np.max(np.abs(f0[:3] - ref[:3])) < 1e-14
        assert np.max(np.abs(dv - ref[3:])) < 1e-14


@pytest.mark.parametrize("kind", ["bcr4bp", "qbcp"])
def test_stm_and_eps_sensitivity_against_finite_differences(kind: str) -> None:
    model = sf.bcr4bp_model() if kind == "bcr4bp" else sf.qbcp_model()
    s = np.array([0.8, 0.05, 0.02, 0.1, 0.3, 0.01])
    t0, t1, eps = 1.3, 3.1, 0.6
    _, stm, sens = sf.flow(model, eps, s, t0, t1, variational=True)
    assert stm is not None and sens is not None
    h = 1e-6
    fd = np.empty((6, 6))
    for j in range(6):
        d = np.zeros(6)
        d[j] = h
        fd[:, j] = (
            sf.flow(model, eps, s + d, t0, t1)[0][0] - sf.flow(model, eps, s - d, t0, t1)[0][0]
        ) / (2 * h)
    assert np.max(np.abs(fd - stm[0])) < 1e-6 * max(1.0, float(np.max(np.abs(stm[0]))))
    se = (
        sf.flow(model, eps + 1e-6, s, t0, t1)[0][0] - sf.flow(model, eps - 1e-6, s, t0, t1)[0][0]
    ) / 2e-6
    assert np.max(np.abs(se - sens[0])) < 1e-7


@pytest.mark.parametrize("kind", ["bcr4bp", "qbcp"])
def test_numba_integrator_matches_scipy_dop853(kind: str) -> None:
    s = np.array([0.8, 0.05, 0.02, 0.1, 0.3, 0.01])
    if kind == "bcr4bp":
        model = sf.bcr4bp_model()
        ref = solve_ivp(
            bcr4bp.bcr4bp_eom,
            (0.0, 2 * model.tg),
            s,
            args=(bcr4bp.andreu_default(),),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        ).y[:, -1]
        got = sf.flow(model, 1.0, s, 0.0, 2 * model.tg)[0][0]
    else:
        model = sf.qbcp_model()
        system = qbcp.qbcp_default()
        st = sf.pv_to_model(model, s, 0.4, 1.0)
        ref = solve_ivp(
            qbcp.qbcp_eom,
            (0.4, 0.4 + model.tg),
            st,
            args=(system,),
            method="DOP853",
            rtol=1e-13,
            atol=1e-13,
        ).y[:, -1]
        got = sf.flow(model, 1.0, st, 0.4, 0.4 + model.tg)[0][0]
    assert np.max(np.abs(got - ref)) < 1e-9


# --- Oshima (2022): the four spatial 1:1 families found from their three-body parent ----------

_OSHIMA_OMEGA = 0.925195985  # Table 1, printed as -0.925195985 (clockwise Sun)
_OSHIMA = bcr4bp.BCR4BPSystem(
    mu=0.0121506683, mu_sun=328900.541, a_sun_nondim=388.811143, omega_sun_nondim=_OSHIMA_OMEGA
)
# Tables 2-5 (p1332): (x, z, vx, vy, vz, theta_S) at the y = 0 crossing with x > 0 and z < 0.
_OSHIMA_ROWS = {
    "vz0-S0": (1.090174251, -0.204803847, 0.0, -2.061909684, 0.0, 0.0),  # Table 4 #1
    "vz0-Spi": (1.090649738, -0.204909100, 0.0, -2.061914819, 0.0, 3.141592654),  # Table 5 #1
    "z0-S0": (1.111203054, -0.180245301, -0.000134973, -2.062045643, -0.000001655, 4.712253115),
    "z0-Spi": (1.111203054, -0.180245301, 0.000134973, -2.062045643, 0.000001655, 1.570932192),
}


@pytest.fixture(scope="module")
def oshima_found() -> dict[str, object]:
    model = sf.bcr4bp_model(_OSHIMA)
    x, z, _, vy, _, _ = _OSHIMA_ROWS["vz0-S0"]
    member = sf.correct_symmetric(model, np.array([x, z, vy]), model.tg)
    parent = sf.Parent("oshima-1:1", member.state, model.tg, 1, 1)
    scr = sf.screens(model, parent)
    mel = sf.melnikov(model, parent)
    zeros = sf.melnikov_zeros(mel)
    orbits = []
    for zz in zeros["zeros"]:
        prob, xs0 = sf.forced_problem(model, parent, zz["tau"])
        br = sf.continue_in_eps(prob, xs0)
        orbits.append((zz, prob, br))
    return {
        "model": model,
        "member": member,
        "parent": parent,
        "screens": scr,
        "mel": mel,
        "zeros": zeros,
        "orbits": orbits,
    }


def _oshima_miss(model: sf.SunModel, prob: sf.Shooting, nodes: FloatArray, row: str) -> float:
    x, z, vx, vy, vz, theta = _OSHIMA_ROWS[row]
    t = ((model.theta_sun0 - theta) % (2.0 * math.pi)) / model.omega_sun
    pv = sf.value_at_clock(prob, nodes, 1.0, t)
    return float(np.max(np.abs(pv - np.array([x, 0.0, z, vx, vy, vz]))))


def test_oshima_parent_is_a_stable_three_body_orbit_away_from_the_printed_states(
    oshima_found: dict[str, object],
) -> None:
    """The parent (no Sun) is 0.022 from the printed eps = 1 state (measured 2.2e-2), passes every
    screen, and is linearly stable (p1331: "the original orbit in the CR3BP (epsilon = 0) is
    linearly stable")."""
    member = oshima_found["member"]
    scr = oshima_found["screens"]
    assert isinstance(member, sf.SymmetricMember) and isinstance(scr, dict)
    x, z, vx, vy, vz, _ = _OSHIMA_ROWS["vz0-S0"]
    assert np.max(np.abs(member.state - np.array([x, 0.0, z, vx, vy, vz]))) > 1e-2
    assert scr["pass"]
    assert not scr["planar"]
    assert scr["parent_floquet"]["max_abs"] < 1.0 + 1e-4


def test_oshima_melnikov_has_four_simple_zeros_at_the_quarter_sun_phases(
    oshima_found: dict[str, object],
) -> None:
    """A doubly symmetric 1:1 orbit: zeros where the Sun is at 0, pi/2, pi, 3 pi/2 at the
    reference point (the printed Sun angles at the corresponding crossings are 0, 1.5709,
    3.1416, 4.7123). The quadrature agrees with the variational value (both ~0 at a zero)."""
    zeros = oshima_found["zeros"]
    assert isinstance(zeros, dict)
    angles = sorted(z["sun_angle"] % (2 * math.pi) for z in zeros["zeros"])
    angles = sorted(a if a < 2 * math.pi - 1e-6 else 0.0 for a in angles)
    assert len(angles) == 4
    for a, want in zip(angles, (0.0, 0.5 * math.pi, math.pi, 1.5 * math.pi), strict=True):
        assert abs(a - want) < 1e-6
    assert all(z["simple"] for z in zeros["zeros"])
    model = oshima_found["model"]
    parent = oshima_found["parent"]
    assert isinstance(model, sf.SunModel) and isinstance(parent, sf.Parent)
    for z in zeros["zeros"]:
        assert (
            abs(sf.melnikov_variational(model, parent, z["tau"])) < 1e-9 * zeros["amplitude"] * 1e3
        )


def test_oshima_each_zero_continues_to_its_own_printed_family(
    oshima_found: dict[str, object],
) -> None:
    """Every zero reaches eps = 1 with no fold, and the orbit found matches ONE printed family at
    the printed Sun angle to 2e-9 (measured 1.6e-10 .. 1.6e-9) and misses the other three by more
    than 0.3 (the discriminator: the vz0-S0 / vz0-Spi pair differ by 0.41, the z0 pair by 0.36).
    """
    model = oshima_found["model"]
    orbits = oshima_found["orbits"]
    assert isinstance(model, sf.SunModel) and isinstance(orbits, list)
    matched = set()
    for _, prob, br in orbits:
        assert br.stop_reason == "reached_target"
        assert not br.folds
        misses = {row: _oshima_miss(model, prob, br.final_nodes, row) for row in _OSHIMA_ROWS}
        best = min(misses, key=misses.__getitem__)
        assert misses[best] < 2e-9
        assert all(m > 0.3 for row, m in misses.items() if row != best)
        matched.add(best)
    assert matched == set(_OSHIMA_ROWS)


def test_oshima_found_orbits_diagnostics(oshima_found: dict[str, object]) -> None:
    """Closure over one Sun period with the Sun back at its starting phase; the Radau closure
    agrees; the z0 orbits are weakly unstable (Fig. 8: about 1.06) and the vz0 ones stable."""
    orbits = oshima_found["orbits"]
    model = oshima_found["model"]
    assert isinstance(orbits, list) and isinstance(model, sf.SunModel)
    for _, prob, br in orbits:
        d = sf.orbit_diagnostics(prob, br.final_nodes, 1.0)
        assert d["sun_phase_returns"]
        assert d["period_in_sun_periods"] == pytest.approx(1.0, abs=1e-12)
        assert d["closure_dop853"] < 1e-10
        assert d["closure_radau"] < 1e-9
        row = min(_OSHIMA_ROWS, key=lambda r: _oshima_miss(model, prob, br.final_nodes, r))
        if row.startswith("z0"):
            assert 1.05 < d["floquet"]["max_abs"] < 1.07
        else:
            assert d["floquet"]["max_abs"] < 1.0 + 1e-8


# --- Leiva & Briozzo (2008): QBCP orbits found from their Table 1 three-body parents -------------

_LB_MU = 0.0121505482
_LB_SECTION_X = 0.836915310  # printed as -0.836915310 (their frame)
# Table 1 (p234): h, y, ydot, p, q (paper frame).
_LB_TABLE_1 = {
    "013": (-1.58740571, -0.0399746624, -0.0441622383, 4, 1),
    "180A_1": (-1.59005198, 0.00520342002, 0.0479318298, 5, 2),
    "180A_2": (-1.57583831, -0.0283283340, 0.117872065, 5, 2),
    "357": (-1.52791268, 0.129037155, 0.115396082, 5, 2),
}
# Table 2 (p238): t_i, x, xdot, y, ydot (paper frame), state at QBCP clock t_i.
_LB_TABLE_2 = {
    "013_t3": (1.92708674, -0.841058432, -0.0710601802, -0.0415648661, -0.0231934953),
    "013_t4": (5.32268367, -0.841255581, -0.0709901229, -0.0417347404, -0.0224675395),
}


def _lb_model() -> sf.SunModel:
    return sf.qbcp_model(dataclasses.replace(qbcp.qbcp_default(), mu=_LB_MU))


def _lb_parent(model: sf.SunModel, name: str) -> sf.Parent:
    h, y, ydot, p, q = _LB_TABLE_1[name]
    x, y, ydot = _LB_SECTION_X, -y, -ydot
    r1 = math.hypot(x + _LB_MU, y)
    r2 = math.hypot(x - 1.0 + _LB_MU, y)
    xdot = math.sqrt(2 * h + x * x + y * y + 2 * (1 - _LB_MU) / r1 + 2 * _LB_MU / r2 - ydot**2)
    pv, res = sf.correct_parent(model, np.array([x, y, 0.0, xdot, ydot, 0.0]), p / q * model.tg)
    assert res < 1e-11
    assert abs(sf.jacobi_pv(pv, _LB_MU) + 2 * h) < 1e-7
    return sf.Parent(name, pv, p / q * model.tg, q, p)


@pytest.fixture(scope="module")
def lb013() -> dict[str, object]:
    model = _lb_model()
    parent = _lb_parent(model, "013")
    scr = sf.screens(model, parent)
    zeros = sf.melnikov_zeros(sf.melnikov(model, parent))
    orbits = []
    for zz in zeros["zeros"]:
        prob, xs0 = sf.forced_problem(model, parent, zz["tau"])
        orbits.append((zz, prob, sf.continue_in_eps(prob, xs0)))
    return {"model": model, "parent": parent, "screens": scr, "zeros": zeros, "orbits": orbits}


def test_leiva_briozzo_013_found_from_its_three_body_parent(lb013: dict[str, object]) -> None:
    """Parent 013 (Table 1, 4 Sun periods, q = 1) passes the screens; the Melnikov function has
    four simple zeros (their Eq. 31 also gives four phases). All four continue to eps = 1. Two of
    the orbits found are the printed 013_t3 and 013_t4: at the printed clock times they agree
    with the printed states to 2e-6 in position and 1e-5 in velocity (measured 1.1e-6 and
    4.9e-6; the printed states themselves close only to 5e-6 .. 6e-5 in this module), and every
    other orbit misses each printed state by more than 0.05. Their largest multiplier 155.88
    (Table 5 prints |s1| = 155.8 for both, s = lambda + 1/lambda) and their closest approach
    to the Earth 137,125 km (printed 137125 and 137126) are reproduced. The other two phases
    (|lambda| about 125, periselene about 2,350 km) are not in the paper and are reported, not
    claimed."""
    model = lb013["model"]
    scr = lb013["screens"]
    zeros = lb013["zeros"]
    orbits = lb013["orbits"]
    assert isinstance(model, sf.SunModel) and isinstance(scr, dict)
    assert isinstance(zeros, dict) and isinstance(orbits, list)
    assert scr["pass"] and scr["planar"]
    assert len(zeros["zeros"]) == 4 and all(z["simple"] for z in zeros["zeros"])
    hits = {}
    for _, prob, br in orbits:
        assert br.stop_reason == "reached_target"
        for row, (t_i, x, xdot, y, ydot) in _LB_TABLE_2.items():
            pv = sf.value_at_clock(prob, br.final_nodes, 1.0, t_i)
            want = np.array([-x, -y, 0.0, -xdot, -ydot, 0.0])
            dpos = float(np.max(np.abs(pv[:3] - want[:3])))
            dvel = float(np.max(np.abs(pv[3:] - want[3:])))
            if dpos < 2e-6 and dvel < 1e-5:
                hits[row] = (prob, br)
            else:
                assert max(dpos, dvel) > 0.05
    assert set(hits) == set(_LB_TABLE_2)
    for prob, br in hits.values():
        res, _, _, stms = prob.evaluate(br.final_nodes, 1.0)
        assert float(np.max(np.abs(res))) < 1e-10
        lam = sf.floquet_report(sf.monodromy(stms))["max_abs"]
        assert abs(lam + 1.0 / lam - 155.8) < 0.1


def test_leiva_briozzo_phases_are_near_the_melnikov_zeros(lb013: dict[str, object]) -> None:
    """Their t_i (the clock time at the Table 1 section point, from the first-order BCP
    condition) lie within 0.005 TU of the driver's QBCP Melnikov zeros (measured 4.6e-3 for
    both)."""
    zeros = lb013["zeros"]
    assert isinstance(zeros, dict)
    taus = [z["tau"] for z in zeros["zeros"]]
    for t_i, *_ in _LB_TABLE_2.values():
        assert min(abs(t_i - t) for t in taus) < 5e-3


def test_leiva_briozzo_5_2_c32_member_completes_as_a_periodic_orbit() -> None:
    """Leiva & Briozzo (p239) obtained 180A_1 (the project's C32 member at 5/2, C = 3.18010396)
    only as a periodic ARC: single shooting over 5 Sun periods did not converge. With multiple
    shooting the driver closes it as a true periodic orbit of period 5 Sun periods (Sun back at
    its starting phase) at both Melnikov phases (measured: closure 7e-12; largest multipliers
    2.5e3 and 3.9e4; periselene 7,348 km at the phase of their arc 180A_1_t1, whose printed
    d_M is 7,371 km; one branch point recorded at eps 0.52 on the other phase). This test runs
    the first phase only."""
    model = _lb_model()
    parent = _lb_parent(model, "180A_1")
    assert parent.laps == 2 and parent.n_sun == 5
    zeros = sf.melnikov_zeros(sf.melnikov(model, parent))
    assert len(zeros["zeros"]) == 2
    zz = min(zeros["zeros"], key=lambda z: abs(z["tau"] - 1.51986327))  # their t1
    prob, xs0 = sf.forced_problem(model, parent, zz["tau"])
    br = sf.continue_in_eps(prob, xs0)
    assert br.stop_reason == "reached_target" and br.final_nodes is not None
    res, _, _, stms = prob.evaluate(br.final_nodes, 1.0)
    assert float(np.max(np.abs(res))) < 1e-10
    assert prob.period == pytest.approx(5 * model.tg, rel=1e-14)
    lam = sf.floquet_report(sf.monodromy(stms))["max_abs"]
    assert 1e3 < lam < 1e4


def test_oshima_four_families_are_two_symmetry_classes(oshima_found: dict[str, object]) -> None:
    """Oshima's Tables 4 and 5 (and 2 and 3) are images of each other under the reflection
    z -> -z (his eq. 8: Table 5 #3 is Table 4 #1 with z negated, at the same Sun angle), so the
    four orbits found are two classes, one per start point (z0 and vz0), not four."""
    orbits = oshima_found["orbits"]
    model = oshima_found["model"]
    assert isinstance(orbits, list) and isinstance(model, sf.SunModel)
    samples, rows = [], []
    for _, prob, br in orbits:
        d = sf.orbit_diagnostics(prob, br.final_nodes, 1.0, radau=False)
        samples.append(d["_clock_samples"])
        rows.append(min(_OSHIMA_ROWS, key=lambda r: _oshima_miss(model, prob, br.final_nodes, r)))
    classes = sf.symmetry_classes(samples, 128)
    assert len(set(classes)) == 2
    by_row = dict(zip(rows, classes, strict=True))
    assert by_row["vz0-S0"] == by_row["vz0-Spi"]
    assert by_row["z0-S0"] == by_row["z0-Spi"]
    assert by_row["vz0-S0"] != by_row["z0-S0"]


def test_leiva_briozzo_013_refined_distances_match_table_5(lb013: dict[str, object]) -> None:
    """Refined closest approaches of the found 013 orbits against Table 5 (p241): d_E 137125 and
    137126 km (centre distances, 384,400 km units; measured 137,125.05 for both); d_M printed
    2729 and 2733 km, which the #896 addendum found to be printed up to 7.5 km above the true
    minimum (a sampled minimum); measured 2,727.16 km. The parent's own refined passes are in
    the screens."""
    orbits = lb013["orbits"]
    scr = lb013["screens"]
    assert isinstance(orbits, list) and isinstance(scr, dict)
    assert not scr["parent_below_moon_surface"] and scr["parent_periselene_km"] > 1737.4
    peris = []
    for _, prob, br in orbits:
        t_i, x, _, y, _ = _LB_TABLE_2["013_t3"]
        pv = sf.value_at_clock(prob, br.final_nodes, 1.0, t_i)
        if float(np.max(np.abs(pv[:2] - np.array([-x, -y])))) > 1e-5:
            continue
        d = sf.orbit_diagnostics(prob, br.final_nodes, 1.0, radau=False)
        assert abs(d["perigee_km"] - 137125.0) < 1.0
        assert 2729.0 - 7.5 <= d["periselene_km"] <= 2729.0
        peris.append(d["periselene_km"])
    assert len(peris) == 1


def test_sub_surface_parent_is_flagged_by_the_refined_periselene() -> None:
    """Review section 8: the Casoliva 2:1(b) member at C = 0.567 (period Tg; state from the #884
    three-body walk) is itself an impact orbit, periselene 1,415.4 km from the Moon's centre by
    the reviewer's separate code, where a minimum sampled at 400 points per TU read 1,810 km. The
    refined periselene here agrees to 1 km and the parent is excluded."""
    model = sf.bcr4bp_model()
    pv, res = sf.correct_parent(
        model, np.array([-1.7957517461647838, 0.0, 0.0, 0.0, 1.9426441634823564, 0.0]), model.tg
    )
    assert res < 1e-11
    shape = sf.minimal_period_and_planarity(
        model, sf.Parent("casoliva-2-1b-low", pv, model.tg, 1, 1)
    )
    assert abs(shape["parent_periselene_km"] - 1415.4) < 1.0
    assert shape["parent_below_moon_surface"]
