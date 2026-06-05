"""Physics-invariant tests for :mod:`cyclerfinder.core.sims_flanagan`.

There is no usable *literature* anchor for the Sims-Flanagan leg model yet (Yam
2010's worked examples are Jupiter rendezvous, out of cycler scope; the Vasile &
Campagnola JBIS tables are transcription-blocked as golden EXPECTED — see
``docs/superpowers/plans/2026-06-05-sims-flanagan-lowthrust.md`` Phase 5). These
tests therefore rest entirely on **physics invariants**, never fabricated
numbers. The load-bearing regression anchor is *zero-thrust reduces to Kepler*:
a zero-Delta-V leg must reproduce :func:`cyclerfinder.core.kepler.propagate`
exactly (it is the same propagator) and close with zero match-point defect.
"""

from __future__ import annotations

from math import exp, pi, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import (
    AU_KM,
    MU_SUN_KM3_S2,
    SECONDS_PER_DAY,
    STANDARD_GRAVITY_KM_S2,
)
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.sims_flanagan import (
    SimsFlanaganError,
    SimsFlanaganLeg,
    final_mass,
    match_point_defect,
    propagate_backward,
    propagate_forward,
    segment_dv_bounds,
)


def _earth_to_mars_leg(
    n_segments: int = 10,
    tof_days: float = 250.0,
    tmax_kn: float = 0.0,
    match_index: int = -1,
) -> SimsFlanaganLeg:
    """Build a leg from Earth's circular state to Mars's, ``tof_days`` later.

    Endpoints come from the circular ephemeris; the start velocity is the true
    heliocentric circular velocity (so a zero-thrust leg is a real Kepler arc).
    ``tmax_kn`` defaults to 0 (coast-only); raise it for thrust-bound tests.
    """
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = tof_days * SECONDS_PER_DAY
    rf, vf = eph.state("M", tof_s)
    return SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n_segments,
        m0_kg=1000.0,
        isp_s=3000.0,
        tmax_kn=tmax_kn,
        match_index=match_index,
    )


def _zero_schedule(n_segments: int) -> np.ndarray:
    return np.zeros((n_segments, 3), dtype=np.float64)


# ---------------------------------------------------------------------------
# Invariant 1: zero-thrust forward propagation == Kepler (the key anchor)
# ---------------------------------------------------------------------------


def test_zero_thrust_forward_matches_kepler() -> None:
    """A zero-Delta-V forward pass to the match point equals direct Kepler."""
    leg = _earth_to_mars_leg(n_segments=8, match_index=8)
    fwd = propagate_forward(leg, _zero_schedule(8))
    # match_index == n_segments => forward covers the whole leg.
    r_kep, v_kep = propagate(leg.r0, leg.v0, leg.tof_s, leg.mu)
    assert float(np.linalg.norm(fwd.r - r_kep)) < 1.0e-3
    assert float(np.linalg.norm(fwd.v - v_kep)) < 1.0e-9


def test_zero_thrust_partial_forward_matches_kepler() -> None:
    """Forward to a mid-leg match point equals Kepler over that sub-interval."""
    n = 10
    leg = _earth_to_mars_leg(n_segments=n, match_index=4)
    fwd = propagate_forward(leg, _zero_schedule(n))
    dt_to_match = leg.match_index * leg.dt_seg_s
    r_kep, v_kep = propagate(leg.r0, leg.v0, dt_to_match, leg.mu)
    assert float(np.linalg.norm(fwd.r - r_kep)) < 1.0e-3
    assert float(np.linalg.norm(fwd.v - v_kep)) < 1.0e-9


# ---------------------------------------------------------------------------
# Invariant 2: zero-thrust leg closes (Lambert / Kepler) with ~0 defect
# ---------------------------------------------------------------------------


def test_zero_thrust_defect_zero_kepler_endpoints() -> None:
    """End state = Kepler-propagated start => zero-Delta-V defect ~ 0."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = 200.0 * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    n = 12
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n,
        m0_kg=500.0,
        isp_s=2500.0,
        tmax_kn=0.0,
    )
    dvs = _zero_schedule(n)
    defect = match_point_defect(leg, dvs, final_mass(leg, dvs))
    assert float(np.linalg.norm(defect[0:3])) < 1.0e-3  # position, km
    assert float(np.linalg.norm(defect[3:6])) < 1.0e-9  # velocity, km/s
    assert abs(defect[6]) < 1.0e-12  # mass, kg (no burns => exactly closed)


def test_zero_thrust_defect_zero_lambert_endpoints() -> None:
    """Lambert-consistent endpoints close ballistically with zero Delta V."""
    eph = Ephemeris(model="circular")
    r0, _ = eph.state("E", 0.0)
    tof_s = 180.0 * SECONDS_PER_DAY
    rf, _ = eph.state("M", tof_s)
    sol = lambert(r0, rf, tof_s)[0]  # single-rev solution
    n = 10
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=sol.v1,
        rf=rf,
        vf=sol.v2,
        tof_s=tof_s,
        n_segments=n,
        m0_kg=750.0,
        isp_s=3200.0,
        tmax_kn=0.0,
    )
    dvs = _zero_schedule(n)
    defect = match_point_defect(leg, dvs, final_mass(leg, dvs))
    assert float(np.linalg.norm(defect[0:3])) < 1.0e-2  # km (Lambert tol scale)
    assert float(np.linalg.norm(defect[3:6])) < 1.0e-8  # km/s


# ---------------------------------------------------------------------------
# Invariant 3: mass monotone non-increasing forward; backward inverts it
# ---------------------------------------------------------------------------


def test_mass_monotone_non_increasing() -> None:
    """final_mass <= m0, strictly less when any segment burns."""
    n = 6
    leg = _earth_to_mars_leg(n_segments=n)
    assert final_mass(leg, _zero_schedule(n)) == pytest.approx(leg.m0_kg)

    dvs = np.zeros((n, 3), dtype=np.float64)
    dvs[2] = np.array([0.05, 0.0, 0.0])  # one real burn
    mf = final_mass(leg, dvs)
    assert mf < leg.m0_kg
    # Tsiolkovsky check on the single burn.
    expected = leg.m0_kg * exp(-0.05 / (STANDARD_GRAVITY_KM_S2 * leg.isp_s))
    assert mf == pytest.approx(expected, rel=1e-12)


def test_backward_mass_inverts_forward_mass() -> None:
    """Backward mass at the match point matches forward mass for feasible mf."""
    n = 8
    leg = _earth_to_mars_leg(n_segments=n, tmax_kn=0.05, match_index=4)
    rng = np.random.default_rng(42)
    dvs = rng.normal(scale=0.01, size=(n, 3))
    mf = final_mass(leg, dvs)
    fwd = propagate_forward(leg, dvs)
    bwd = propagate_backward(leg, dvs, mf)
    assert bwd.mass_kg == pytest.approx(fwd.mass_kg, rel=1e-12)


# ---------------------------------------------------------------------------
# Invariant 4: per-segment Delta V respects the thrust-capability bound
# ---------------------------------------------------------------------------


def test_segment_dv_bound_formula_and_growth() -> None:
    """Bound equals (T_max/m_i)*dt_seg and grows as mass falls."""
    n = 5
    leg = _earth_to_mars_leg(n_segments=n, tmax_kn=0.1)
    dvs = np.full((n, 3), 0.02, dtype=np.float64)  # equal burns each segment
    bounds = segment_dv_bounds(leg, dvs)
    # First-segment bound uses m0 exactly.
    assert bounds[0] == pytest.approx((leg.tmax_kn / leg.m0_kg) * leg.dt_seg_s, rel=1e-12)
    # Mass strictly decreases => bound strictly increases.
    assert np.all(np.diff(bounds) > 0.0)


def test_capability_bound_respected_by_construction() -> None:
    """A schedule built at the per-segment bound satisfies |dv_i| <= bound_i."""
    n = 6
    leg = _earth_to_mars_leg(n_segments=n, tmax_kn=0.08)
    # A realistic low-thrust duty cycle burns a small slice of the (large)
    # ballistic per-segment capability. Size each segment burn off the
    # first-segment capability (avoiding a runaway where shrinking mass inflates
    # later caps); a small absolute Delta V keeps every segment under its bound.
    cap0 = (leg.tmax_kn / leg.m0_kg) * leg.dt_seg_s
    dv_per_seg = 0.02 * cap0  # 2% of the first-segment capability
    dvs = np.zeros((n, 3), dtype=np.float64)
    dvs[:, 0] = dv_per_seg
    bounds = segment_dv_bounds(leg, dvs)
    mags = np.linalg.norm(dvs, axis=1)
    # The bound grows along the leg (mass falls), so every segment is under it.
    assert np.all(mags < bounds)


# ---------------------------------------------------------------------------
# Invariant 6: energy conservation on coast; energy change tracks the impulse
# ---------------------------------------------------------------------------


def _specific_energy(r: np.ndarray, v: np.ndarray, mu: float) -> float:
    return float(np.dot(v, v)) / 2.0 - mu / float(np.linalg.norm(r))


def test_coast_conserves_specific_energy() -> None:
    """A zero-thrust leg conserves specific orbital energy end to end."""
    n = 10
    leg = _earth_to_mars_leg(n_segments=n, match_index=n)
    e0 = _specific_energy(leg.r0, leg.v0, leg.mu)
    fwd = propagate_forward(leg, _zero_schedule(n))
    e1 = _specific_energy(fwd.r, fwd.v, leg.mu)
    assert abs(e1 - e0) / abs(e0) < 1.0e-8


# ---------------------------------------------------------------------------
# Invariant 7: forward-then-backward over a full leg round-trips the start
# ---------------------------------------------------------------------------


def test_full_leg_round_trip() -> None:
    """Forward over [0,N) then backward over [0,N) returns the start state."""
    n = 8
    # match_index = N => forward covers everything; then a backward pass with
    # match_index = 0 covers everything in reverse and must land at the start.
    rng = np.random.default_rng(7)
    dvs = rng.normal(scale=0.005, size=(n, 3))

    leg_fwd = _earth_to_mars_leg(n_segments=n, tmax_kn=0.05, match_index=n)
    fwd = propagate_forward(leg_fwd, dvs)
    mf = final_mass(leg_fwd, dvs)

    # A leg whose end state is the forward result, propagated fully backward.
    leg_bwd = SimsFlanaganLeg(
        r0=leg_fwd.r0,
        v0=leg_fwd.v0,
        rf=fwd.r,
        vf=fwd.v,
        tof_s=leg_fwd.tof_s,
        n_segments=n,
        m0_kg=leg_fwd.m0_kg,
        isp_s=leg_fwd.isp_s,
        tmax_kn=leg_fwd.tmax_kn,
        match_index=0,
    )
    bwd = propagate_backward(leg_bwd, dvs, mf)
    assert float(np.linalg.norm(bwd.r - leg_fwd.r0)) < 1.0e-2
    assert float(np.linalg.norm(bwd.v - leg_fwd.v0)) < 1.0e-8
    assert bwd.mass_kg == pytest.approx(leg_fwd.m0_kg, rel=1e-12)


# ---------------------------------------------------------------------------
# Hyperbolic-regime smoke: zero-thrust still tracks Kepler off a fast state
# ---------------------------------------------------------------------------


def test_zero_thrust_hyperbolic_tracks_kepler() -> None:
    """Zero-Delta-V leg from a hyperbolic state matches Kepler propagation."""
    r0 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    v_esc = sqrt(2.0 * MU_SUN_KM3_S2 / AU_KM)
    v0 = np.array([0.0, 1.3 * v_esc, 0.0], dtype=np.float64)
    tof_s = 40.0 * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    n = 6
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n,
        m0_kg=300.0,
        isp_s=2000.0,
        tmax_kn=0.0,
        match_index=n,
    )
    fwd = propagate_forward(leg, _zero_schedule(n))
    assert float(np.linalg.norm(fwd.r - rf)) < 1.0e-2
    assert float(np.linalg.norm(fwd.v - vf)) < 1.0e-9


# ---------------------------------------------------------------------------
# Defect / mass non-trivial: a burn schedule that does NOT close has nonzero defect
# ---------------------------------------------------------------------------


def test_nonzero_schedule_breaks_ballistic_closure() -> None:
    """Adding burns to a ballistically-closed leg opens the match-point defect."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    tof_s = 200.0 * SECONDS_PER_DAY
    rf, vf = propagate(r0, v0, tof_s, MU_SUN_KM3_S2)
    n = 10
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=tof_s,
        n_segments=n,
        m0_kg=500.0,
        isp_s=2500.0,
        tmax_kn=0.1,
        match_index=5,
    )
    dvs = np.zeros((n, 3), dtype=np.float64)
    dvs[1] = np.array([0.02, 0.0, 0.0])  # a burn only on the forward half
    # mf consistent with the full schedule, so the mass defect stays ~0; the
    # position/velocity defect must open because the forward half now diverges.
    defect = match_point_defect(leg, dvs, final_mass(leg, dvs))
    assert float(np.linalg.norm(defect[0:3])) > 1.0  # km, clearly nonzero
    assert abs(defect[6]) < 1.0e-9  # mass still consistent


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------


def test_default_match_index_is_temporal_midpoint() -> None:
    leg = _earth_to_mars_leg(n_segments=10)
    assert leg.match_index == 5


def test_dt_seg_property() -> None:
    leg = _earth_to_mars_leg(n_segments=8, tof_days=240.0)
    assert leg.dt_seg_s == pytest.approx(240.0 * SECONDS_PER_DAY / 8.0)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"tof_s": 0.0}, "tof_s"),
        ({"n_segments": 0}, "n_segments"),
        ({"m0_kg": 0.0}, "m0_kg"),
        ({"isp_s": -1.0}, "isp_s"),
        ({"tmax_kn": -1.0}, "tmax_kn"),
        ({"match_index": 99}, "match_index"),
    ],
)
def test_config_validation_rejects_bad_params(kwargs: dict[str, float], match: str) -> None:
    base = dict(
        r0=np.array([AU_KM, 0.0, 0.0]),
        v0=np.array([0.0, 30.0, 0.0]),
        rf=np.array([0.0, AU_KM, 0.0]),
        vf=np.array([-30.0, 0.0, 0.0]),
        tof_s=100.0 * SECONDS_PER_DAY,
        n_segments=10,
        m0_kg=1000.0,
        isp_s=3000.0,
        tmax_kn=0.1,
    )
    base.update(kwargs)
    with pytest.raises(SimsFlanaganError, match=match):
        SimsFlanaganLeg(**base)  # type: ignore[arg-type]


def test_schedule_shape_validation() -> None:
    leg = _earth_to_mars_leg(n_segments=5)
    with pytest.raises(SimsFlanaganError, match="shape"):
        propagate_forward(leg, np.zeros((4, 3)))


def test_backward_rejects_nonpositive_mass() -> None:
    leg = _earth_to_mars_leg(n_segments=5)
    with pytest.raises(SimsFlanaganError, match="mf_kg"):
        propagate_backward(leg, _zero_schedule(5), 0.0)


def test_defect_is_length_seven() -> None:
    n = 6
    leg = _earth_to_mars_leg(n_segments=n)
    dvs = _zero_schedule(n)
    defect = match_point_defect(leg, dvs, final_mass(leg, dvs))
    assert defect.shape == (7,)


def test_full_circle_period_round_trip() -> None:
    """A coast leg over one Earth period returns to the start (closes)."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    a_km = float(np.linalg.norm(r0))
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
    n = 12
    rf, vf = propagate(r0, v0, period_s, MU_SUN_KM3_S2)
    leg = SimsFlanaganLeg(
        r0=r0,
        v0=v0,
        rf=rf,
        vf=vf,
        tof_s=period_s,
        n_segments=n,
        m0_kg=400.0,
        isp_s=2200.0,
        tmax_kn=0.0,
        match_index=n,
    )
    fwd = propagate_forward(leg, _zero_schedule(n))
    assert float(np.linalg.norm(fwd.r - r0)) < 1.0  # 1 km, matches kepler.py
    assert float(np.linalg.norm(fwd.v - v0)) < 1.0e-6
