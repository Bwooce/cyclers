"""Tests for the #942/#943 two-working-body generator and date-residual corrector."""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.kepler import propagate
from cyclerfinder.search.cycler_assembly import half_rev_intersection
from cyclerfinder.search.generic_return import RussellModel
from cyclerfinder.search.hollister_menning_1970 import (
    block_types,
    build_cycle,
    load_table3,
    pick_branches,
)
from cyclerfinder.search.two_working_body import (
    Cycle,
    HalfRevLeg,
    LambertLeg,
    MeanElementSystem,
    ResonantLeg,
    correct_dates,
    cycle_flybys,
    encounter_self_consistency,
    half_rev_arrival,
    half_rev_vectors,
    heliocentric_circular,
    resonant_circle,
    resonant_vec,
)

DAY = SECONDS_PER_DAY


def test_circular_periods_and_kepler3() -> None:
    sysm = heliocentric_circular({"E": 1.0, "V": 0.61520})
    r0, _ = sysm.state("V", 0.0)
    r1, _ = sysm.state("V", sysm.period_s("V"))
    assert np.allclose(r0, r1, atol=1e-3)
    # Kepler III: the circular speed of each body is consistent with its period
    for c in "EV":
        r, v = sysm.state(c, 0.0)
        per = 2 * math.pi * float(np.linalg.norm(r)) / float(np.linalg.norm(v))
        assert per == pytest.approx(sysm.period_s(c), rel=1e-12)


@pytest.mark.parametrize("phi", [0.0, 1.0, 2.5, 4.0])
def test_resonant_circle_is_a_true_full_rev_return(phi: float) -> None:
    sysm = heliocentric_circular({"E": 1.0, "V": 0.61520})
    leg = ResonantLeg("V", 2, 1)  # 2 Venus periods, 1 spacecraft revolution
    vinf = 6.0
    circ = resonant_circle(sysm, leg, 0.0, vinf)
    assert circ is not None
    u = resonant_vec(circ, phi)
    assert float(np.linalg.norm(u)) == pytest.approx(vinf, rel=1e-12)
    r, w = sysm.state("V", 0.0)
    dt = 2 * sysm.period_s("V")
    r1, _ = propagate(r, w + u, dt, sysm.mu)
    rb, _ = sysm.state("V", dt)
    assert float(np.linalg.norm(r1 - rb)) < 1.0  # km


def test_half_rev_lands_on_body_and_matches_russell_z() -> None:
    sysm = heliocentric_circular({"E": 1.0, "V": 0.61520})
    leg = HalfRevLeg("E", half_periods=3, k_sc=1, via_peri=True, sign=0)
    vinf_kms = 10.0
    us = half_rev_vectors(sysm, leg, 0.0, vinf_kms)
    assert len(us) == 2
    r, w = sysm.state("E", 0.0)
    for u in us:
        assert float(np.linalg.norm(u)) == pytest.approx(vinf_kms, rel=1e-10)
        dt = 1.5 * sysm.period_s("E")
        r1, _ = propagate(r, w + u, dt, sysm.mu)
        rb, _ = sysm.state("E", dt)
        assert float(np.linalg.norm(r1 - rb)) < 10.0  # km, of 1.5e8
        ua = half_rev_arrival(sysm, leg, 0.0, u)
        assert float(np.linalg.norm(ua)) == pytest.approx(vinf_kms, rel=1e-8)
    # Independent cross-check: Russell's half-rev intersection (Eqs 2.23-2.25)
    # gives the same component along the body velocity for the same conic.
    v_e = float(np.linalg.norm(w))
    v_sc = w + us[0]
    rn = float(np.linalg.norm(r))
    a_km = 1.0 / (2.0 / rn - float(v_sc @ v_sc) / sysm.mu)
    rm = RussellModel()
    au = 149597870.7
    vu = math.sqrt(sysm.mu / au)  # km/s per canonical speed unit (1 AU circle)
    tip, _ = half_rev_intersection(rm, "E", a_km / au, vinf_kms / vu)
    assert float(us[0] @ (w / v_e)) / vu == pytest.approx(float(tip[2]), abs=2e-4)


def test_cycle_rejects_broken_chain() -> None:
    with pytest.raises(ValueError):
        Cycle((LambertLeg("E", "V"), LambertLeg("E", "V")), 1.0)


def test_table3_block_structure_and_transcription_fix() -> None:
    t = load_table3()
    assert len(t) == 15
    assert t[1][11].planet == "E"  # printed "E 3163" (transcribed as V)
    for orbit, rows in t.items():
        assert "".join(r.planet for r in rows) == "EEVVV" * 5 + "E"
        assert len(block_types(orbit, rows)) == 5
    # H&M text: orbits 3-8 have symmetric returns at Earth, orbit 1 has none
    assert all(et == "FR" for et, _ in block_types(1, t[1]))
    assert sum(et == "SY" for et, _ in block_types(3, t[3])) == 5


def test_hm1970_orbit1_corrector_reproduces_printed_dates_and_turns() -> None:
    """Positive control (expected values: H&M 1970 Table 3, orbit 1, the paper).

    The corrector, seeded at the printed dates, converges on the
    inclined-elliptic mean-element model to dates within 3 d and turn angles
    within 3 deg of every printed encounter (minimax full-rev directions)."""
    sysm = MeanElementSystem()
    rows = load_table3()[1]
    branches, _ = pick_branches(sysm, 1, rows)
    cyc, x0 = build_cycle(sysm, 1, rows, branches)
    sol = correct_dates(sysm, cyc, x0)
    assert sol.converged
    assert encounter_self_consistency(sysm, cyc, sol.x) < 1.0
    fl = cycle_flybys(sysm, cyc, sol.x)
    assert fl is not None and len(fl) == 25
    for r in rows[:-1]:
        f = min(
            (f for f in fl if f.body == r.planet),
            key=lambda f: min(abs(f.t_s / DAY - r.date), abs(f.t_s / DAY - 5844 - r.date)),
        )
        d = min(abs(f.t_s / DAY - r.date), abs(f.t_s / DAY - 5844 - r.date))
        assert d < 3.0
        assert f.turn_deg == pytest.approx(r.theta_deg, abs=3.0)


def test_kepler_step_matches_universal_variables() -> None:
    from cyclerfinder.search.two_working_body import kepler_step

    mu = 1.327e11
    r0 = np.array([1.4e8, 2e7, 1e6])
    v0 = np.array([-3.0, 31.0, 0.5])
    for dt in (1e5, 3e6, 2e7, 7e7, -4e6):
        a = kepler_step(r0, v0, dt, mu)
        b = propagate(r0, v0, dt, mu)
        assert float(np.linalg.norm(a[0] - b[0])) < 1e-3
        assert float(np.linalg.norm(a[1] - b[1])) < 1e-9


def test_x1_one_body_control_russell_strange_gancal5() -> None:
    """X1 control (expected values: Russell & Strange 2009 Tables 2, 3, 5).

    GanCal#5: Ganymede generic returns g1.50425 and g3.74691 (both one-rev),
    Callisto massless. The date corrector must find the split, V_inf 3.24 km/s
    and the 328 km minimum flyby altitude, and carry a massless Callisto
    encounter placed on the cycler without changing them."""
    from cyclerfinder.search.two_working_body import CircularSystem, FlybyBody, gate_cycle
    from cyclerfinder.search.two_working_body_enum import Zero, assess, place_massless_target
    from cyclerfinder.verify.turn_gate_closures import (
        RS_BODIES,
        RS_GENERIC_CYCLERS,
        RS_PRIMARY_GM,
    )

    mu = RS_PRIMARY_GM["Jupiter"]
    g, c = RS_BODIES["Ganymede"], RS_BODIES["Callisto"]

    def a_of(p: float) -> float:
        return float((mu * (p / (2 * math.pi)) ** 2) ** (1 / 3))

    sysm = CircularSystem(
        mu,
        {
            "Ganymede": (a_of(g.ideal_period_s), g.ideal_period_s, 0.0),
            "Callisto": (a_of(c.ideal_period_s), c.ideal_period_s, 0.0),
        },
        flyby_overrides={"Ganymede": FlybyBody("Ganymede", g.gm_km3_s2, g.radius_km, 0.0)},
    )
    row = RS_GENERIC_CYCLERS[0]
    assert row.rs_id == "GanCal#5"
    p = g.ideal_period_s
    period = sum(leg.n_body_revs for leg in row.legs) * p
    # Blind over the four branch pairs and five seeds: the gate-passing zeros
    # that reproduce the published row must all be the ("low", "high") pair.
    matches = []
    for b1, b2 in itertools.product(("low", "high"), repeat=2):
        cyc_try = Cycle(
            (LambertLeg("Ganymede", "Ganymede", 1, b1), LambertLeg("Ganymede", "Ganymede", 1, b2)),
            period,
        )
        for split in (1.3, 1.55, 1.8, 2.1, 2.4):
            s_try = correct_dates(sysm, cyc_try, np.array([0.0, split * p]))
            fl_try = cycle_flybys(sysm, cyc_try, s_try.x) if s_try.converged else None
            if not fl_try:
                continue
            rep_try = gate_cycle(sysm, fl_try)
            if (
                rep_try.status == "pass"
                and abs(fl_try[0].vinf_kms - row.vinf_flyby_kms) < 0.015
                and abs(rep_try.gate.min_required_alt_km - row.min_flyby_alt_km) < 2.0
            ):
                matches.append((b1, b2, cyc_try, s_try))
    assert matches
    assert {(m[0], m[1]) for m in matches} == {("low", "high")}
    _, _, cyc, sol = matches[0]
    assert (sol.x[1] - sol.x[0]) / p == pytest.approx(row.legs[0].n_body_revs, abs=1e-4)
    placed = place_massless_target(sysm, cyc, sol.x, "Callisto")
    assert placed
    s2, c2, seed = placed[0]
    z = correct_dates(s2, c2, seed)
    assert z.converged
    a = assess(s2, Zero(c2, z.x, z.max_abs_residual_kms))
    assert a.report is not None
    assert a.report.gate.min_required_alt_km == pytest.approx(row.min_flyby_alt_km, abs=2.0)
    assert a.r_min_km == pytest.approx(row.min_dist_primary_km, rel=1e-3)
    assert a.r_max_km == pytest.approx(row.max_dist_primary_km, rel=1e-3)
    assert a.max_encounter_miss_km < 1.0


def test_hollister_3h_circular_recall() -> None:
    """Recall control (Hollister & Menning 1970 p.1194): in the circular
    coplanar model a 3.2-yr (two E-V synodic periods) orbit with a symmetric
    return at Earth and two full-revolution returns at Venus exists with
    ballistic flybys. Consistency check only: the 1.37-yr symmetric return
    H&M inserted (p.1195) belongs to their inclined-elliptic sequential
    modification, not to the circular orbit."""
    from cyclerfinder.search.two_working_body_enum import assess, solve_structure

    s = heliocentric_circular({"E": 1.0, "V": 0.61520})
    syn = s.synodic_s("E", "V")
    cyc = Cycle(
        (
            LambertLeg("E", "E", 1, "high"),
            LambertLeg("E", "V"),
            ResonantLeg("V"),
            ResonantLeg("V"),
            LambertLeg("V", "E"),
        ),
        2 * syn,
    )
    zs = solve_structure(s, cyc, phase_period_s=syn, n_phase=48, n_split=16, n_refine=300)
    passing = [z for z in zs if assess(s, z).status == "pass"]
    assert passing
    # the symmetric return runs from the first Lambert start to the second
    sy_years = [(z.x[1] - z.x[0]) / (365.25 * DAY) for z in passing]
    assert min(abs(y - 1.37) for y in sy_years) < 0.03


def test_x1_ganeur43_blind_through_the_enumerator() -> None:
    """X1 control (expected: Russell & Strange 2009 Tables 3, 5, row GanEur#43:
    V_inf 1.87 km/s, minimum flyby altitude 8861 km, 564,558-1,072,330 km).

    Every k = 2 Ganymede-Europa structure with Europa massless and up to two
    revolutions per transfer is enumerated at the production settings; every
    gate-passing zero must be this cycler (the split labels are one conic)."""
    import importlib.util
    from pathlib import Path

    from cyclerfinder.search.two_working_body_enum import (
        CatalogueSpec,
        assess,
        solve_structure,
        structures,
    )

    path = Path(__file__).resolve().parents[2] / "scripts" / "run_942_enumerate.py"
    spec = importlib.util.spec_from_file_location("run_942_enumerate", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sysm, a, b = mod.cell_system("ge1")
    syn = sysm.synodic_s(a, b)
    passing = []
    for cyc in structures(
        sysm,
        a,
        b,
        2,
        max_returns={a: 0, b: 0},
        spec={c: CatalogueSpec() for c in (a, b)},
        transfer_revs=(0, 1, 2),
    ):
        for z in solve_structure(
            sysm, cyc, phase_period_s=syn, n_phase=36, n_split=12, n_refine=40
        ):
            ass = assess(sysm, z)
            if ass.status == "pass":
                passing.append(ass)
    assert passing
    for ass in passing:
        assert ass.vinf_kms["Ganymede"] == pytest.approx(1.87, abs=0.015)
        assert ass.report is not None
        assert ass.report.gate.min_required_alt_km == pytest.approx(8861.0, abs=2.0)
        assert ass.r_min_km == pytest.approx(564_558.0, rel=1e-3)
        assert ass.r_max_km == pytest.approx(1_072_330.0, rel=1e-3)


def test_r1c_one_body_control_venmar45_through_the_enumerator() -> None:
    """R1(c) control (expected: Russell & Strange 2007, AAS 07-118, Tables 2, 3, 5,
    row VenMar#45: V_inf 8.22 / 12.96 km/s, minimum Venus flyby altitude
    19,784 km, maximum distance 341,571,371 km, period 2 x 333.9 d).

    Venus-Mars, Mars massless, k = 2, no returns, transfers of 0-1 revolutions:
    every gate-passing zero must be VenMar#45."""
    import importlib.util
    from pathlib import Path

    from cyclerfinder.search.two_working_body_enum import (
        CatalogueSpec,
        assess,
        solve_structure,
        structures,
    )

    path = Path(__file__).resolve().parents[2] / "scripts" / "run_942_enumerate.py"
    spec = importlib.util.spec_from_file_location("run_942_enumerate", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sysm, a, b = mod.cell_system("vm")
    syn = sysm.synodic_s(a, b)
    assert 2 * syn / DAY == pytest.approx(667.8, abs=0.1)
    passing = []
    for cyc in structures(
        sysm,
        a,
        b,
        2,
        max_returns={a: 0, b: 0},
        spec={c: CatalogueSpec() for c in (a, b)},
        transfer_revs=(0, 1),
    ):
        for z in solve_structure(
            sysm, cyc, phase_period_s=syn, n_phase=36, n_split=12, n_refine=40
        ):
            ass = assess(sysm, z)
            if ass.status == "pass":
                passing.append(ass)
    assert passing
    for ass in passing:
        assert ass.vinf_kms["V"] == pytest.approx(8.22, abs=0.01)
        assert ass.vinf_kms["M"] == pytest.approx(12.96, abs=0.01)
        assert ass.report is not None
        assert ass.report.gate.min_required_alt_km == pytest.approx(19_784.0, abs=2.0)
        assert ass.r_max_km == pytest.approx(341_571_371.0, rel=1e-5)
