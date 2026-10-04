"""#888 -- the demanded-turn gate (``verify/turn_gate.py``) and its extractors.

GOLDEN DISCIPLINE. Expected values come from mathematical identities or from
published tables, never from this project's earlier output:

* identities: zero turn on a same-conic pass, a known rotation is recovered,
  the bend formula against a hand calculation, frame invariance, the
  equal-magnitude impulse identity;
* POSITIVE CONTROL 1, McConaghy, Longuski and Byrnes, "Analysis of a Broad
  Class of Earth-Mars Cycler Trajectories", AIAA 2002-4420 (filed in the
  private paper corpus as
  mcconaghy-longuski-byrnes-2002-analysis-broad-class-earth-mars-cycler-trajectories-AIAA-2002-4420.pdf),
  Table 4 p.6: required and maximum turn angles (200 km Earth flyby, p.5) of
  19 nPr cyclers, three of them ballistic (footnote e); and McConaghy,
  Longuski and Byrnes, JSR 41(4) 2004, p.627 (filed as
  mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.11939.pdf):
  the Aldrin (1L1) cycler "required flyby altitude is -1731 km";
* POSITIVE CONTROL 2, Russell and Strange, "Cycler Trajectories in Planetary
  Moon Systems", JGCD 32(1) 2009 (filed as
  russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf),
  Tables 3-4 (pp.149-150) "Min flyby alt." and Tables 5-6 (p.151) leg
  nomenclature, constants from Table 2 (p.148);
* REGRESSION: the six withdrawn Uranian rows (``data/withdrawn/``) fail at
  every encounter. Only the verdict and a conservative bound are asserted.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.verify.turn_gate import (
    Encounter,
    body_constants,
    demanded_turn_gate,
    encounters_from_vinf_nodes,
    impulse_beyond_bend_kms,
    required_periapsis_alt_km,
    rotate_towards,
    to_body_local,
)
from cyclerfinder.verify.turn_gate_closures import (
    RS_GENERIC_CYCLERS,
    RS_TITAN_MIN_ALT_KM,
    mcconaghy_npr_cycler,
    russell_strange_generic_cycler,
    symmetric_closure,
)

REPO = Path(__file__).resolve().parents[2]


def _rot_z(v: np.ndarray, ang: float) -> np.ndarray:
    c, s = math.cos(ang), math.sin(ang)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]])


# ---------------------------------------------------------------------------
# Identities
# ---------------------------------------------------------------------------


def test_same_conic_pass_demands_zero_turn() -> None:
    v = np.array([1.3, -0.4, 0.2])
    rep = demanded_turn_gate([Encounter.for_body("Oberon", v, v)])
    (e,) = rep.encounters
    assert e.demanded_turn_deg == pytest.approx(0.0, abs=1e-6)
    assert e.ratio == 0.0
    assert e.required_alt_km == math.inf
    assert e.impulse_beyond_bend_kms == pytest.approx(0.0, abs=1e-12)
    assert rep.turn_feasible and rep.ballistic


@pytest.mark.parametrize("angle_deg", [1.0, 30.0, 97.0, 179.0])
def test_known_rotation_is_recovered(angle_deg: float) -> None:
    v = np.array([0.9, 0.3, 0.0])
    w = _rot_z(v, math.radians(angle_deg))
    (e,) = demanded_turn_gate([Encounter.for_body("Oberon", v, w)]).encounters
    assert e.demanded_turn_deg == pytest.approx(angle_deg, abs=1e-9)


def test_bend_formula_against_hand_calculation() -> None:
    # Hand calculation, Earth GM 398600.435507 km^3/s^2, R 6378.137 km, 200 km floor,
    # V_inf 6.54 km/s: rp v^2 / GM = 6578.137 * 42.7716 / 398600.435507 = 0.705863,
    # sin(d/2) = 1 / 1.705863 = 0.586213, d = 2 asin(.) = 2 * 35.8888 = 71.7775 deg.
    rp = 6378.137 + 200.0
    k = rp * 6.54**2 / 398600.435507
    expected = 2.0 * math.degrees(math.asin(1.0 / (1.0 + k)))
    assert expected == pytest.approx(71.7775, abs=1e-3)
    v = np.array([6.54, 0.0, 0.0])
    (e,) = demanded_turn_gate([Encounter.for_body("E", v, _rot_z(v, 0.5))]).encounters
    assert e.available_bend_deg == pytest.approx(expected, abs=1e-9)
    assert e.alt_floor_km == PLANETS["E"].safe_alt_km == 200.0


def test_required_altitude_inverts_the_bend_formula() -> None:
    bc = body_constants("Titania")
    v = 1.7
    for alt in (50.0, 300.0, 2000.0):
        rp = bc.radius_km + alt
        turn = 2.0 * math.asin(1.0 / (1.0 + rp * v * v / bc.mu_km3_s2))
        assert required_periapsis_alt_km(bc.mu_km3_s2, bc.radius_km, v, turn) == pytest.approx(
            alt, abs=1e-6
        )
    # A 180 degree turn needs periapsis at the centre.
    assert required_periapsis_alt_km(bc.mu_km3_s2, bc.radius_km, v, math.pi) == pytest.approx(
        -bc.radius_km
    )


def test_frame_invariance() -> None:
    rng = np.random.default_rng(888)
    v_in = np.array([1.1, -0.2, 0.35])
    v_out = np.array([0.4, 1.0, -0.1]) * (np.linalg.norm(v_in) / np.linalg.norm([0.4, 1.0, -0.1]))
    base = demanded_turn_gate([Encounter.for_body("Ganymede", v_in, v_out)]).encounters[0]
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rot = demanded_turn_gate([Encounter.for_body("Ganymede", q @ v_in, q @ v_out)]).encounters[0]
    assert rot.demanded_turn_deg == pytest.approx(base.demanded_turn_deg, abs=1e-10)
    assert rot.available_bend_deg == pytest.approx(base.available_bend_deg, abs=1e-12)
    assert rot.impulse_beyond_bend_kms == pytest.approx(base.impulse_beyond_bend_kms, abs=1e-12)


def test_body_local_frame_round_trip_and_rotation_invariance() -> None:
    r = np.array([7.0e5, 1.0e5, 0.0])
    w = np.array([-1.0, 7.0, 0.0])
    v = np.array([0.3, -1.2, 0.1])
    loc = to_body_local(v, r, w)
    assert np.linalg.norm(loc) == pytest.approx(np.linalg.norm(v))
    ang = 1.234
    loc_rot = to_body_local(_rot_z(v, ang), _rot_z(r, ang), _rot_z(w, ang))
    assert np.allclose(loc, loc_rot, atol=1e-12)


def test_equal_magnitude_impulse_identity() -> None:
    v = np.array([2.0, 0.0, 0.0])
    demanded, bend = math.radians(80.0), math.radians(25.0)
    w = _rot_z(v, demanded)
    expected = 2.0 * 2.0 * math.sin(0.5 * (demanded - bend))
    assert impulse_beyond_bend_kms(v, w, bend) == pytest.approx(expected, abs=1e-12)
    # Within the bend: no impulse.
    assert impulse_beyond_bend_kms(v, w, math.radians(81.0)) == pytest.approx(0.0, abs=1e-12)
    # A pure magnitude mismatch costs exactly the mismatch.
    assert impulse_beyond_bend_kms(v, 1.3 * v, bend) == pytest.approx(0.6, abs=1e-12)


def test_rotate_towards_handles_antiparallel() -> None:
    v = np.array([1.0, 0.0, 0.0])
    out = rotate_towards(v, -v, math.radians(90.0))
    assert np.linalg.norm(out) == pytest.approx(1.0)
    assert float(out @ v) == pytest.approx(0.0, abs=1e-12)


def test_magnitude_mismatch_reported_and_smaller_speed_used() -> None:
    v_in = np.array([1.0, 0.0, 0.0])
    v_out = _rot_z(np.array([1.2, 0.0, 0.0]), math.radians(5.0))
    rep = demanded_turn_gate([Encounter.for_body("Oberon", v_in, v_out)])
    (e,) = rep.encounters
    assert e.magnitude_mismatch_kms == pytest.approx(0.2)
    bc = body_constants("Oberon")
    bmax = 2.0 * math.asin(1.0 / (1.0 + (bc.radius_km + bc.alt_floor_km) * 1.0 / bc.mu_km3_s2))
    assert e.available_bend_deg == pytest.approx(math.degrees(bmax))
    assert e.turn_feasible and not rep.ballistic
    assert math.isnan(e.impulse_periapsis_kms)


def test_vinf_nodes_adapter_with_wrap() -> None:
    v = np.array([5.0, 0.0, 0.0])
    nodes = {
        "b0_out": v,
        "b1_in": v,
        "b1_out": _rot_z(v, math.radians(20.0)),
        "b2_in": _rot_z(v, math.radians(40.0)),
    }
    encs = encounters_from_vinf_nodes(nodes, ("E", "M", "E"), wrap_rotation_rad=math.radians(10.0))
    rep = demanded_turn_gate(encs)
    assert [e.body for e in rep.encounters] == ["M", "E"]
    assert rep.encounters[0].demanded_turn_deg == pytest.approx(20.0)
    assert rep.encounters[1].demanded_turn_deg == pytest.approx(30.0)  # 40 vs b0_out rotated by 10


# ---------------------------------------------------------------------------
# Positive control 1: McConaghy, Longuski and Byrnes (2002) Table 4, p.6
# ---------------------------------------------------------------------------

# (n, P, r, aphelion AU, V_inf at Earth km/s, required turn deg, max turn deg)
MCCONAGHY_TABLE4 = [
    (1, "L", 1, 2.23, 6.54, 84, 72),
    (2, "L", 2, 2.33, 10.06, 134, 44),
    (2, "L", 3, 1.51, 5.65, 135, 82),
    (3, "L", 4, 1.89, 11.78, 167, 35),
    (3, "L", 5, 1.45, 7.61, 167, 62),
    (3, "S", 5, 1.52, 12.27, 167, 33),
    (4, "S", 5, 1.82, 11.23, 167, 38),
    (4, "S", 6, 1.53, 8.51, 167, 54),
    (5, "S", 4, 2.49, 10.62, 134, 41),
    (5, "S", 5, 2.09, 9.08, 134, 50),
    (5, "S", 6, 1.79, 7.51, 135, 62),
    (5, "S", 7, 1.54, 5.86, 135, 79),
    (5, "S", 8, 1.34, 4.11, 136, 103),
    (6, "S", 4, 2.81, 7.93, 83, 59),
    (6, "S", 5, 2.37, 6.94, 84, 68),
    (6, "S", 6, 2.04, 5.96, 84, 78),
    (6, "S", 7, 1.78, 4.99, 85, 90),
    (6, "S", 8, 1.57, 4.02, 85, 104),
    (6, "S", 9, 1.40, 3.04, 86, 120),
]
MCCONAGHY_BALLISTIC = {"6S7", "6S8", "6S9"}  # Table 4 footnote e

# Integer-degree columns: +-1 deg (rounding plus the unstated constants of
# their model); two-decimal V_inf and aphelion columns: +-0.01 / +-0.01.
TURN_TOL_DEG = 1.0


@pytest.mark.parametrize(("n", "p", "r", "ra", "vinf", "req", "mx"), MCCONAGHY_TABLE4)
def test_mcconaghy_table4_turns_and_verdict(
    n: int, p: str, r: int, ra: float, vinf: float, req: float, mx: float
) -> None:
    c = mcconaghy_npr_cycler(n, p, r)  # type: ignore[arg-type]
    assert c is not None
    # Identify the arc by the columns that do not depend on the turn.
    assert c.vinf_earth_kms == pytest.approx(vinf, abs=0.01)
    assert c.aphelion_au == pytest.approx(ra, abs=0.011)
    rep = demanded_turn_gate([c.encounter])
    (e,) = rep.encounters
    assert e.demanded_turn_deg == pytest.approx(req, abs=TURN_TOL_DEG)
    assert e.available_bend_deg == pytest.approx(mx, abs=TURN_TOL_DEG)
    assert rep.turn_feasible == (f"{n}{p}{r}" in MCCONAGHY_BALLISTIC)


def test_mcconaghy_aldrin_required_altitude() -> None:
    """JSR 2004 p.627: the 1L1 required flyby altitude is -1731 km.

    Tolerance: dh/d(turn) is about 135 km/deg near 84 deg at 6.54 km/s, so the
    integer-degree turn column alone allows about +-70 km; 100 km is used.
    """
    c = mcconaghy_npr_cycler(1, "L", 1)
    assert c is not None
    (e,) = demanded_turn_gate([c.encounter]).encounters
    assert e.required_alt_km == pytest.approx(-1731.0, abs=100.0)
    assert not e.turn_feasible
    assert e.impulse_beyond_bend_kms > 0.0


# ---------------------------------------------------------------------------
# Positive control 2: Russell and Strange (2009) Tables 3-6
# ---------------------------------------------------------------------------

# GanIo#403 carries the branch label "Ll" (two candidate branches with V_inf
# 4.28 and 4.23 km/s against a published 4.29); its rebuilt minimum altitude is
# 557 km against 540 km published, so it is held to a looser bound and the
# ambiguity is reported, not hidden.
RS_ALT_TOL_KM = {"GanIo#403": 25.0}


@pytest.mark.parametrize("row", RS_GENERIC_CYCLERS, ids=[r.rs_id for r in RS_GENERIC_CYCLERS])
def test_russell_strange_min_flyby_altitude(row) -> None:  # type: ignore[no-untyped-def]
    floor = RS_TITAN_MIN_ALT_KM if row.flyby_body == "Titan" else 0.0
    rb = russell_strange_generic_cycler(row, alt_floor_km=floor)
    # Identification checks that do not involve the turn.
    assert rb.period_days == pytest.approx(row.period_days, abs=0.06)
    assert min(rb.leg_vinf_kms) == pytest.approx(row.vinf_flyby_kms, abs=0.015)
    assert min(rb.leg_rp_km) == pytest.approx(row.min_dist_primary_km, rel=1e-3)
    assert max(rb.leg_ra_km) == pytest.approx(row.max_dist_primary_km, rel=1e-3)
    rep = demanded_turn_gate(rb.encounters)
    # The golden: the needed altitude reproduces "Min flyby alt. at body A".
    tol = RS_ALT_TOL_KM.get(row.rs_id, 2.0)
    assert rep.min_required_alt_km == pytest.approx(row.min_flyby_alt_km, abs=tol)
    # Published as ballistic in the ideal model: the gate passes it at the
    # source floor (1000 km at Titan, p.144; at the Galilean moons every
    # published altitude is positive, so the surface).
    assert rep.turn_feasible


# ---------------------------------------------------------------------------
# Regression: the six withdrawn Uranian (1,1) rows
# ---------------------------------------------------------------------------

WITHDRAWN = sorted((REPO / "data" / "withdrawn").glob("*-1-1-uranian-quasi-cycler-2026.yaml"))


def _rebuild_withdrawn(path: Path):  # type: ignore[no-untyped-def]
    row = yaml.safe_load(path.read_text())[0]
    seq = row["sequence_canonical"].split("-")
    tof = float(row["legs"][0]["tof_days"])
    want = [float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"]]
    best = None
    # rel_offset is not a field of the rows and the #312 row's legs block lists
    # n_revs 0 against its (1,1) name: choose by reproducing the row's V_inf.
    for rel in (0.0, 180.0):
        for nr in ((0, 0), (1, 1)):
            c = symmetric_closure(
                "Uranus", seq[0], seq[1], tof_days=tof, n_rev=nr, rel_offset_deg=rel
            )
            if c is None:
                continue
            err = max(abs(a - b) for a, b in zip(c.stored_convention_vinf, want, strict=True))
            if best is None or err < best[0]:
                best = (err, c)
    assert best is not None
    return best


def test_six_withdrawn_rows_present() -> None:
    assert len(WITHDRAWN) == 6


@pytest.mark.parametrize("path", WITHDRAWN, ids=[p.stem for p in WITHDRAWN])
def test_withdrawn_rows_fail_at_every_encounter(path: Path) -> None:
    err, c = _rebuild_withdrawn(path)
    assert err < 1e-3  # the row's 4-decimal V_inf are reproduced
    for floor in (None, 50.0):
        rep = demanded_turn_gate(c.encounters, alt_floor_km=floor)
        assert not rep.turn_feasible
        assert len(rep.encounters) == 2
        for e in rep.encounters:
            assert not e.turn_feasible
            assert e.ratio > 1.5
            assert e.required_alt_km < 0.0  # inside the moon


def test_symmetric_closure_wrap_matches_local_frame_when_commensurate() -> None:
    """For tof = n T_syn / 2 the chain repeats in the rotating frame, so the
    directly solved third leg and the rotated first leg give the same turn."""
    path = REPO / "data" / "withdrawn" / "ariel-oberon-1-1-uranian-quasi-cycler-2026.yaml"
    _err, c = _rebuild_withdrawn(path)
    direct = demanded_turn_gate(c.encounters).encounters[1].demanded_turn_deg
    assert c.wrap_local_check_deg == pytest.approx(direct, abs=1e-6)
