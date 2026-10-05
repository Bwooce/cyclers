"""Tests for :mod:`cyclerfinder.core.kepler` (universal-variable propagator)."""

from __future__ import annotations

from math import pi, sqrt

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate


def test_kepler_zero_dt_identity() -> None:
    """``dt == 0`` returns the original state exactly (as a copy)."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    r, v = propagate(r0, v0, 0.0)
    assert np.array_equal(r, r0)
    assert np.array_equal(v, v0)
    # And it's a copy, not a view.
    r[0] = 12345.0
    assert r0[0] != 12345.0


def test_kepler_circular_period() -> None:
    """Earth's circular state propagates one Keplerian period back to itself."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    a_km = float(np.linalg.norm(r0))
    period_s = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
    r1, v1 = propagate(r0, v0, period_s)
    # 1 km positional tolerance per plan §4.6.
    assert float(np.linalg.norm(r1 - r0)) < 1.0
    assert float(np.linalg.norm(v1 - v0)) < 1.0e-6


def test_kepler_reversibility() -> None:
    """``propagate(propagate(s, +dt), -dt) == s`` within numerical noise."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("M", 0.0)
    dt = 137.0 * SECONDS_PER_DAY
    r1, v1 = propagate(r0, v0, dt)
    r2, v2 = propagate(r1, v1, -dt)
    assert float(np.linalg.norm(r2 - r0)) < 1.0
    assert float(np.linalg.norm(v2 - v0)) < 1.0e-6


def test_kepler_energy_conservation() -> None:
    """Specific orbital energy is conserved to high relative precision."""
    eph = Ephemeris(model="circular")
    r0, v0 = eph.state("E", 0.0)
    energy0 = float(np.dot(v0, v0)) / 2.0 - MU_SUN_KM3_S2 / float(np.linalg.norm(r0))
    for t_days in (10.0, 100.0, 365.0, 700.0):
        r, v = propagate(r0, v0, t_days * SECONDS_PER_DAY)
        energy = float(np.dot(v, v)) / 2.0 - MU_SUN_KM3_S2 / float(np.linalg.norm(r))
        # Relative tolerance 1e-8 per plan §4.6.
        assert abs(energy - energy0) / abs(energy0) < 1.0e-8


def test_kepler_hyperbolic() -> None:
    """A fabricated hyperbolic state round-trips identity under +dt then -dt."""
    # Start at 1 AU with velocity 1.4x escape (well into the hyperbolic regime).
    r0 = np.array([AU_KM, 0.0, 0.0], dtype=np.float64)
    v_esc = sqrt(2.0 * MU_SUN_KM3_S2 / AU_KM)
    v0 = np.array([0.0, 1.4 * v_esc, 0.0], dtype=np.float64)
    dt = 30.0 * SECONDS_PER_DAY
    r1, v1 = propagate(r0, v0, dt)
    r2, v2 = propagate(r1, v1, -dt)
    assert float(np.linalg.norm(r2 - r0)) < 1.0
    assert float(np.linalg.norm(v2 - v0)) < 1.0e-6


def test_kepler_negative_dt_directly() -> None:
    """Backward propagation from a future state lands at the past state."""
    eph = Ephemeris(model="circular")
    r_at_zero, v_at_zero = eph.state("E", 0.0)
    dt = 200.0 * SECONDS_PER_DAY
    r_at_dt, v_at_dt = eph.state("E", dt)
    # Propagate the future state backward by dt; should land at t=0.
    r_back, v_back = propagate(r_at_dt, v_at_dt, -dt)
    assert float(np.linalg.norm(r_back - r_at_zero)) < 1.0
    assert float(np.linalg.norm(v_back - v_at_zero)) < 1.0e-6


def _dop853_two_body(
    r0: np.ndarray, v0: np.ndarray, dt: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Independent reference: integrate the two-body ODE (no universal variables)."""
    from scipy.integrate import solve_ivp

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        r = y[:3]
        return np.concatenate([y[3:], -mu * r / float(np.linalg.norm(r)) ** 3])

    sol = solve_ivp(
        rhs, (0.0, dt), np.concatenate([r0, v0]), method="DOP853", rtol=1e-13, atol=1e-13
    )
    return sol.y[:3, -1], sol.y[3:, -1]


def test_kepler_near_parabolic_sweep_934() -> None:
    """#934: propagate must not fail near e = 1 (either side), and must be right.

    Sweeps e - 1 over log-spaced magnitudes 1e-12..1e-2 on both sides of 1 and
    exactly 1, against a DOP853 reference, with energy / angular-momentum
    conservation checked on the result.
    """
    mu = 1.0
    deltas = [0.0]
    for mag in np.logspace(-12, -2, 11):
        deltas.extend([float(mag), -float(mag)])
    failures: list[str] = []
    for rp in (1.0, 7.0):
        for delta in deltas:
            e = 1.0 + delta
            r0 = np.array([rp, 0.0, 0.0])
            v0 = np.array([0.0, sqrt(mu * (1.0 + e) / rp), 0.0])
            eps0 = 0.5 * float(v0 @ v0) - mu / rp
            h0 = np.cross(r0, v0)
            for dt in (0.1, 0.7, 2.5, 9.0, -2.5):
                try:
                    r, v = propagate(r0, v0, dt, mu)
                except Exception as exc:
                    failures.append(f"rp={rp} e-1={delta:.3e} dt={dt}: {type(exc).__name__}")
                    continue
                r_ref, v_ref = _dop853_two_body(r0, v0, dt, mu)
                scale = max(1.0, float(np.linalg.norm(r_ref)))
                if float(np.linalg.norm(r - r_ref)) > 1e-8 * scale:
                    failures.append(f"rp={rp} e-1={delta:.3e} dt={dt}: position mismatch")
                if float(np.linalg.norm(v - v_ref)) > 1e-8:
                    failures.append(f"rp={rp} e-1={delta:.3e} dt={dt}: velocity mismatch")
                eps = 0.5 * float(v @ v) - mu / float(np.linalg.norm(r))
                if abs(eps - eps0) > 1e-10 * max(abs(eps0), mu / rp):
                    failures.append(f"rp={rp} e-1={delta:.3e} dt={dt}: energy drift")
                h = np.cross(r, v)
                if float(np.linalg.norm(h - h0)) > 1e-10 * float(np.linalg.norm(h0)):
                    failures.append(f"rp={rp} e-1={delta:.3e} dt={dt}: angular momentum drift")
    assert not failures, f"{len(failures)} failing cases:\n" + "\n".join(failures)


def test_shepperd_stm_near_parabolic_934() -> None:
    """#934: the STM propagator shares the guess code and must also survive e ~ 1."""
    from cyclerfinder.core.kepler_stm import shepperd_stm

    mu = 1.0
    for rp in (1.0, 7.0):
        for delta in (1e-8, -1e-8, 1e-7, 0.0):
            r0 = np.array([rp, 0.0, 0.0])
            v0 = np.array([0.0, sqrt(mu * (2.0 + delta) / rp), 0.0])
            for dt in (0.7, 2.5, 9.0, -2.5):
                r_ref, v_ref = _dop853_two_body(r0, v0, dt, mu)
                r, v, _phi = shepperd_stm(r0, v0, dt, mu)
                assert float(np.linalg.norm(r - r_ref)) < 1e-8 * max(1.0, float(np.linalg.norm(r)))
                assert float(np.linalg.norm(v - v_ref)) < 1e-8


# #963 reproducer (papercut 2026-10-05-twobody-gen-opus-kepler-propagate-nonconvergence): a
# 173-day leg of an elliptic, near-radial heliocentric orbit (a = 0.754 au, perihelion about
# 9e6 km, period 239 d) from a #942 Venus-to-Earth sample. Unguarded Newton on f(chi) diverged
# because f'(chi) = r is small near the first iterate.
_R0_963 = np.array([-88039995.65452953, -62916481.70866319, 0.0])
_V0_963 = np.array([-19.984711843308432, -29.60574045546609, 0.0])
_MU_963 = 132717453059.67786
_DT_963 = 14953799.40036204


def _kepler_equation_state(
    r0: np.ndarray, v0: np.ndarray, dt: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Independent elliptic reference: Kepler's equation M = E - e sin E solved by bisection
    (monotone in E, no starting guess), then the perifocal state. Planar orbits only."""
    rn = float(np.linalg.norm(r0))
    h = float(r0[0] * v0[1] - r0[1] * v0[0])
    a = 1.0 / (2.0 / rn - float(v0 @ v0) / mu)
    ev = ((float(v0 @ v0) - mu / rn) * r0 - float(r0 @ v0) * v0) / mu
    e = float(np.linalg.norm(ev))
    w = float(np.arctan2(ev[1], ev[0]))
    cos_e0 = (1.0 - rn / a) / e
    sin_e0 = float(r0 @ v0) / (e * sqrt(mu * a))
    e0 = float(np.arctan2(sin_e0, cos_e0))
    n = sqrt(mu / a**3)
    m = e0 - e * np.sin(e0) + n * dt
    lo, hi = m - e - 1.0, m + e + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - e * np.sin(mid) < m:
            lo = mid
        else:
            hi = mid
    ea = 0.5 * (lo + hi)
    b = a * sqrt(1.0 - e * e)
    xp, yp = a * (np.cos(ea) - e), b * np.sin(ea)
    edot = n / (1.0 - e * np.cos(ea))
    vxp, vyp = -a * np.sin(ea) * edot, b * np.cos(ea) * edot
    sgn = 1.0 if h >= 0.0 else -1.0
    cw, sw = np.cos(w), np.sin(w)
    r = np.array([cw * xp - sw * sgn * yp, sw * xp + cw * sgn * yp, 0.0])
    v = np.array([cw * vxp - sw * sgn * vyp, sw * vxp + cw * sgn * vyp, 0.0])
    return r, v


def _assert_matches_kepler_equation(
    r: np.ndarray, v: np.ndarray, r0: np.ndarray, v0: np.ndarray, dt: float, mu: float
) -> None:
    r_ref, v_ref = _kepler_equation_state(r0, v0, dt, mu)
    scale_r = float(np.linalg.norm(r_ref))
    a = 1.0 / (2.0 / float(np.linalg.norm(r0)) - float(v0 @ v0) / mu)
    # Both sides lose digits in the near-radial pass: the reference through the angle
    # (eccentric anomaly) and E - e sin E near E = 0, the propagator through f and g.
    assert float(np.linalg.norm(r - r_ref)) < 1e-8 * max(scale_r, a)
    assert float(np.linalg.norm(v - v_ref)) < 1e-8 * float(np.linalg.norm(v_ref)) + 1e-8 * sqrt(
        mu / a
    )


def test_kepler_near_radial_elliptic_reproducer_963() -> None:
    r, v = propagate(_R0_963, _V0_963, _DT_963, _MU_963)
    _assert_matches_kepler_equation(r, v, _R0_963, _V0_963, _DT_963, _MU_963)
    r_b, _v_b = propagate(r, v, -_DT_963, _MU_963)
    # round trip measured 4e-15 relative
    assert float(np.linalg.norm(r_b - _R0_963)) < 1e-12 * float(np.linalg.norm(_R0_963))


def test_kepler_near_radial_elliptic_sweep_963() -> None:
    """Elliptic orbits with e from 0.9 to 0.9999, started at several true anomalies (both
    branches), propagated by fractions of the period up to two revolutions, both signs of dt,
    against Kepler's equation. Every case must converge and match."""
    mu = 1.0
    failures: list[str] = []
    for e in (0.9, 0.99, 0.999, 0.9999):
        a = 1.0
        p = a * (1.0 - e * e)
        period = 2.0 * pi * sqrt(a**3 / mu)
        for nu in (-3.0, -2.0, -0.3, 0.0, 0.3, 2.0, 3.0):
            rn = p / (1.0 + e * np.cos(nu))
            r0 = np.array([rn * np.cos(nu), rn * np.sin(nu), 0.0])
            vf = sqrt(mu / p)
            v0 = np.array([-vf * np.sin(nu), vf * (e + np.cos(nu)), 0.0])
            for frac in (0.1, 0.37, 0.5, 0.72, 0.95, 1.3, 2.0):
                for sign in (1.0, -1.0):
                    dt = sign * frac * period
                    try:
                        r, v = propagate(r0, v0, dt, mu)
                        if e < 0.9995:
                            _assert_matches_kepler_equation(r, v, r0, v0, dt, mu)
                        else:
                            # e = 0.9999 from r0 = 1e-4: the f and g sums cancel by about
                            # a / r0 = 1e4, so the universal-variable state is good to about
                            # 1e-8 a (measured 9.4e-9 against a 40-digit Kepler solve; the
                            # float Kepler-equation reference has 2e-10). Convergence is what
                            # #963 is about; accuracy here is held at 1e-7 a.
                            r_ref, _ = _kepler_equation_state(r0, v0, dt, mu)
                            assert float(np.linalg.norm(r - r_ref)) < 1e-7 * a
                    except Exception as exc:
                        failures.append(f"e={e} nu={nu} dt={dt:.3f}: {type(exc).__name__}")
    assert not failures, f"{len(failures)} failing cases:\n" + "\n".join(failures[:40])


def test_shepperd_stm_near_radial_elliptic_963() -> None:
    from cyclerfinder.core.kepler_stm import shepperd_stm

    r, v, _phi = shepperd_stm(_R0_963, _V0_963, _DT_963, _MU_963)
    _assert_matches_kepler_equation(r, v, _R0_963, _V0_963, _DT_963, _MU_963)
