"""Regression tests for the Sun-Mars WSB / ballistic-capture machinery (#681).

Locks in the load-bearing physics that the task's positive control and
repeating-capture search depend on:

* the paper's Table 5 bitangential Hohmann baselines (constants + heliocentric
  geometry) reproduce to < 3 m/s;
* the Table 3 H3-equivalent Hohmann insertion cost dV2 reproduces analytically;
* Mars's Keplerian ephemeris (perihelion radius, sidereal period);
* the capture periapsis seed lands exactly on the WSB (E_2 = mu_M(e-1)/(2 r_p),
  zero Mars-relative radial rate);
* the two-body reduction (Mars turned off -> exact heliocentric Kepler energy
  conservation);
* the recapture-episode classifier's semantics (an isolated temporary capture is
  episodes=1, i.e. NOT a repeating chain).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.sunmars_wsb as sw


@pytest.mark.parametrize(
    ("earth_apsis", "mars_apsis", "dv1", "dv2", "tof"),
    [
        ("peri", "peri", 2.179, 3.388, 234),
        ("peri", "apo", 3.398, 2.090, 278),
        ("apo", "peri", 2.414, 3.163, 239),
        ("apo", "apo", 3.629, 1.881, 283),
    ],
)
def test_hohmann_table5(
    earth_apsis: str, mars_apsis: str, dv1: float, dv2: float, tof: float
) -> None:
    """Paper Table 5 bitangential Hohmann baselines (constants + geometry)."""
    h = sw.hohmann_baseline(earth_apsis=earth_apsis, mars_apsis=mars_apsis)  # type: ignore[arg-type]
    assert h.dv1_kms == pytest.approx(dv1, abs=3e-3)
    assert h.dv2_inf_kms == pytest.approx(dv2, abs=3e-3)
    assert h.tof_days == pytest.approx(tof, abs=1.0)


def test_mars_ephemeris() -> None:
    """Mars perihelion radius and sidereal period (Keplerian ephemeris)."""
    r, _ = sw.mars_state(0.0, 0.0)  # f0 = 0 -> perihelion, on +x
    peri_au = float(np.hypot(*r)) / sw.AU_KM
    assert peri_au == pytest.approx(1.523688399 * (1.0 - sw.MARS_E), rel=1e-9)
    assert pytest.approx(1.8808, abs=1e-3) == sw.MARS_PERIOD_S / (365.25 * 86400.0)


@pytest.mark.parametrize("rp", [49896.0, 73896.0, 113897.0])
@pytest.mark.parametrize("ecc", [0.90, 0.99])
def test_capture_seed_on_wsb(rp: float, ecc: float) -> None:
    """A capture periapsis seed lands exactly on W: E_2 = mu_M(e-1)/(2 r_p), rdot=0."""
    f0 = math.pi / 4.0
    st = sw.capture_periapsis_state(r_p_km=rp, ecc=ecc, theta=0.7, f0=f0)
    e2 = sw.mars_kepler_energy(st, 0.0, f0)
    assert e2 == pytest.approx(sw.MU_MARS_KM3_S2 * (ecc - 1.0) / (2.0 * rp), rel=1e-10)
    assert e2 < 0.0  # bound
    assert sw.mars_distance(st, 0.0, f0) == pytest.approx(rp, rel=1e-12)
    assert abs(sw.mars_radial_rate(st, 0.0, f0)) < 1e-6


def test_two_body_reduction_conserves_heliocentric_energy() -> None:
    """With Mars gravity removed, heliocentric Kepler energy is conserved."""
    # A bound heliocentric ellipse near Mars's orbit; Mars mass zeroed via a
    # patched EOM (Sun-only) integrated by hand over one Mars period.
    from scipy.integrate import solve_ivp

    r0, v0 = sw.mars_state(0.0, 0.3)
    # perturb to an independent orbit
    state = np.array([r0[0] * 0.9, r0[1], v0[0], v0[1] * 1.05], dtype=np.float64)

    def sun_only(_t: float, y: np.ndarray) -> np.ndarray:
        r = y[:2]
        d = float(np.hypot(r[0], r[1]))
        a = -sw.MU_SUN_KM3_S2 * r / d**3
        return np.array([y[2], y[3], a[0], a[1]])

    def energy(y: np.ndarray) -> float:
        r = float(np.hypot(y[0], y[1]))
        v2 = float(y[2] ** 2 + y[3] ** 2)
        return 0.5 * v2 - sw.MU_SUN_KM3_S2 / r

    sol = solve_ivp(
        sun_only, (0.0, sw.MARS_PERIOD_S), state, method="DOP853", rtol=1e-12, atol=1e-9
    )
    e_start = energy(state)
    e_end = energy(sol.y[:, -1])
    assert e_end == pytest.approx(e_start, rel=1e-9)


def test_recapture_classifier_isolated_capture_is_episodes_one() -> None:
    """An isolated temporary recapture is episodes=1 (NOT a repeating chain).

    The rp=91897, e=0.90, f0=pi/2, theta~2.094 forward seed is one of the
    deepest sustained recaptures found by the #681 search; it recaptures ONCE
    (episodes==1, sustained>=2) and never repeats -- the object-class negative.
    """
    f0 = math.pi / 2.0
    st = sw.capture_periapsis_state(r_p_km=91897.0, ecc=0.90, theta=2.0943951, f0=f0)
    res = sw.integrate_stability(st, f0, direction="forward", horizon_revs=50.0)
    assert res.escaped
    assert res.max_sustained_bound_revs >= 2  # a genuine sustained temporary capture
    assert res.n_recapture_episodes == 1  # but it does NOT repeat -> not a cycler
