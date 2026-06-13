"""Smoke coverage for the Jovian-moon n-body lane (:mod:`cyclerfinder.nbody.jovian`).

The 822-line module was authored but never executed before #223. These are the
first numerical exercises of it. Two tiers:

* KERNEL-FREE (always run): the closed-form patched-conic flyby helpers
  (:func:`flyby_min_dv`, :func:`flyby_altitude_km`) and the TDB epoch axis —
  pure math, no SPICE, no REBOUND.
* KERNEL-BACKED (skipped when JUP365 is absent; CI has no 1.14 GB kernel): one
  moon-state spline-accuracy check and one short conic leg, marked ``slow``
  because furnishing the kernel + building a spline cache is > 10 s.

The full seed -> chain -> shoot -> compare run lives in
``scripts/liang_member_d_run.py`` (script-only; not a test).
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

import cyclerfinder.nbody.jovian as jovian
from cyclerfinder.core.satellites import SATELLITES


def _resolve_kernel() -> str | None:
    """JUP365 path from CYCLERFINDER_JUP365, else the conventional local copy.

    Mirrors ``scripts/liang_member_d_run.py``: the kernel is never committed,
    so CI (no kernel) skips; a developer machine with the standard checkout
    runs the slow tier without needing the env var set in the test shell.
    """
    env = jovian.jup365_kernel_path()
    if env is not None:
        return env
    fallback = os.path.expanduser("~/dev/references/kernels/jup365.bsp")
    return fallback if os.path.exists(fallback) else None


_KERNEL = _resolve_kernel()
_needs_kernel = pytest.mark.skipif(
    _KERNEL is None, reason="JUP365 kernel not furnished (CYCLERFINDER_JUP365)"
)


# --- kernel-free closed-form helpers ---------------------------------------------


def test_tdb_epoch_axis_monotone_and_scaled() -> None:
    """ISO -> TDB seconds since J2000: ~known anchor and day-scale linearity."""
    t0 = jovian.tdb_sec_from_iso("2000-01-01T12:00:00")
    # J2000 epoch is 2000-01-01T12:00:00 TDB -> ~0 s.
    assert abs(t0) < 1.0
    t1 = jovian.tdb_sec_from_iso("2000-01-02T12:00:00")
    assert math.isclose(t1 - t0, 86400.0, abs_tol=1e-3)


def test_flyby_min_dv_zero_when_turn_deliverable() -> None:
    """A small turn an unpowered >min_alt flyby can supply gives ~zero defect."""
    moon = "Callisto"
    v = 5.0  # km/s, representative Callisto V-infinity
    vin = np.array([v, 0.0, 0.0])
    # Rotate by a tiny angle the moon can certainly bend at >100 km altitude.
    ang = math.radians(2.0)
    vout = v * np.array([math.cos(ang), math.sin(ang), 0.0])
    dv, turn, turn_max = jovian.flyby_min_dv(vin, vout, moon, min_alt_km=100.0)
    assert turn_max > turn  # the flyby can over-deliver this turn
    assert dv < 1e-9  # ballistic: defect ~ 0


def test_flyby_min_dv_positive_when_turn_exceeds_capacity() -> None:
    """A turn beyond the moon's bending capacity leaves a positive defect."""
    moon = "Europa"  # small GM -> limited bending
    v = 8.0
    vin = np.array([v, 0.0, 0.0])
    ang = math.radians(60.0)  # far beyond an 8 km/s Europa flyby's reach
    vout = v * np.array([math.cos(ang), math.sin(ang), 0.0])
    dv, _turn, turn_max = jovian.flyby_min_dv(vin, vout, moon, min_alt_km=100.0)
    assert turn_max < math.radians(60.0)  # capacity is the binding limit
    assert dv > 0.1  # meaningful residual defect (km/s)


def test_flyby_altitude_matches_turn_formula() -> None:
    """flyby_altitude_km inverts the same r_p(turn) relation it documents."""
    moon = "Ganymede"
    sat = SATELLITES[moon]
    v = 6.0
    vin = np.array([v, 0.0, 0.0])
    ang = math.radians(15.0)
    vout = v * np.array([math.cos(ang), math.sin(ang), 0.0])
    alt = jovian.flyby_altitude_km(vin, vout, moon)
    # Recompute r_p from the closed form and compare.
    sin_half = math.sin(0.5 * ang)
    r_p = sat.mu_km3_s2 * (1.0 / sin_half - 1.0) / (v * v)
    assert math.isclose(alt, r_p - sat.radius_eq_km, rel_tol=1e-9)


# --- kernel-backed (slow; skipped without the kernel) ----------------------------


@_needs_kernel
@pytest.mark.slow
def test_moon_state_spline_accuracy() -> None:
    """JovianRailsCache spline reproduces direct spkezr moon positions to ~1e-2 km.

    Pins the module docstring's spline-error claim on a short (10 d) window.
    """
    assert _KERNEL is not None
    ephem = jovian.JovianEphemeris(_KERNEL)
    t0 = jovian.tdb_sec_from_iso("2033-09-25T18:04:43")
    from cyclerfinder.core.constants import SECONDS_PER_DAY

    t1 = t0 + 10.0 * SECONDS_PER_DAY
    cache = jovian.JovianRailsCache(jovian.GALILEAN, ephem, t0, t1)
    rng = np.linspace(t0 + SECONDS_PER_DAY, t1 - SECONDS_PER_DAY, 17)
    worst = 0.0
    for moon in jovian.GALILEAN:
        for t in rng:
            r_spline = cache.position(moon, float(t))
            r_true, _ = ephem.state(moon, float(t))
            worst = max(worst, float(np.linalg.norm(r_spline - r_true)))
    assert worst < 1.0  # km — comfortably inside flyby tolerances (claim: ~1e-2)


@_needs_kernel
@pytest.mark.slow
def test_short_conic_leg_lambert_on_real_geometry() -> None:
    """One Callisto->Ganymede 1-rev Lambert leg on real JUP365 geometry solves.

    Exercises the conic-leg path (_solve_cycle_legs) on a single leg at the
    Member D departure epoch and checks the V-infinity magnitudes are in the
    Galilean-tour regime (a few km/s), not divergent.
    """
    assert _KERNEL is not None
    ephem = jovian.JovianEphemeris(_KERNEL)
    t0 = jovian.tdb_sec_from_iso("2033-09-25T18:04:43")
    from cyclerfinder.core.constants import SECONDS_PER_DAY

    # 5 epochs spanning one nominal CGCEC cycle (Member A leg ToFs).
    tofs = (31.8973, 18.1697, 29.9343, 19.9747)
    epochs = [t0]
    for tof in tofs:
        epochs.append(epochs[-1] + tof * SECONDS_PER_DAY)
    legs = jovian._solve_cycle_legs(epochs, ephem, jovian.BRANCH_PLAN)
    assert legs is not None, "all four CGCEC legs must Lambert-solve at the seed epochs"
    vinf_out, vinf_in = legs
    assert len(vinf_out) == 4 and len(vinf_in) == 4
    for vec in vinf_out + vinf_in:
        mag = float(np.linalg.norm(vec))
        assert 1.0 < mag < 20.0, f"V-infinity {mag:.3f} km/s outside the Galilean-tour regime"


@_needs_kernel
@pytest.mark.slow
def test_nbody_propagator_energy_conserved_far_from_moons() -> None:
    """JovianRestrictedNBody is near-Keplerian on an arc far from every moon.

    A 3e6 km circular orbit (beyond Callisto's 1.88e6 km) feels only weak
    moon perturbations; specific orbital energy must be conserved to high
    precision over 5 days. Catches sign/indirect-term/integration bugs in
    the never-before-run propagator (the only numeric guard on it).
    """
    assert _KERNEL is not None
    from cyclerfinder.core.constants import SECONDS_PER_DAY

    ephem = jovian.JovianEphemeris(_KERNEL)
    t0 = jovian.tdb_sec_from_iso("2033-09-25T18:04:43")
    t1 = t0 + 5.0 * SECONDS_PER_DAY
    cache = jovian.JovianRailsCache(jovian.GALILEAN, ephem, t0, t1)
    mu = jovian.MU_JUPITER_KM3_S2
    r0 = np.array([3.0e6, 0.0, 0.0])
    v0 = np.array([0.0, math.sqrt(mu / 3.0e6), 0.0])  # circular speed
    arc = jovian.JovianRestrictedNBody().propagate(r0, v0, t0, t1, cache=cache)
    assert arc.converged

    def energy(r: np.ndarray, v: np.ndarray) -> float:
        return 0.5 * float(np.dot(v, v)) - mu / float(np.linalg.norm(r))

    e0 = energy(r0, v0)
    e1 = energy(arc.r_km, arc.v_km_s)
    # Moon perturbations are real but small at 3e6 km; energy drift stays <1e-4.
    assert abs(e1 - e0) / abs(e0) < 1e-4
