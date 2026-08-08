"""N-body Phase A: propagator interface + rails-force shape (plan Phase A)."""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.nbody.forces import rails_acceleration


def test_sun_only_acceleration_is_central_inverse_square() -> None:
    # No perturbers -> pure -mu_sun r_hat / r^2.
    r = np.array([1.5e8, 0.0, 0.0])  # ~1 AU on +x
    a = rails_acceleration(r, t_sec=0.0, bodies=(), ephem=None)
    expected = -MU_SUN_KM3_S2 / (1.5e8**2)
    assert a[0] == pytest.approx(expected, rel=1e-12)
    assert a[1] == 0.0 and a[2] == 0.0


def test_propagator_backend_selectable() -> None:
    pytest.importorskip("rebound")
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    assert prop.backend == "rebound"
    with pytest.raises(ValueError, match="backend"):
        RestrictedNBody("tudat")  # deferred slot: not wired yet (design Q1)


@pytest.mark.slow
def test_rails_acceleration_matches_integrator_in_far_field() -> None:
    """Far-field agreement (C3): the public unsoftened ``rails_acceleration`` and
    the integrator's softened callback are identical in the cruise regime.

    With the spacecraft ~1 AU from the Sun and ~10^7 km from Earth — far outside
    Earth's safe-flyby periapsis, so the callback's ``|d|`` clamp never engages —
    the effective acceleration REBOUND applies over a short integration must equal
    the unsoftened public function. (They diverge ONLY inside the safe periapsis,
    by design; that near-flyby regime is the callback's clamped path.) The applied
    acceleration is recovered by finite-differencing the velocity over a short step
    (REBOUND 5.x exposes no single-eval acceleration accessor).
    """
    pytest.importorskip("rebound")
    import rebound

    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.nbody.forces import RailsEphemerisCache, rails_acceleration
    from cyclerfinder.nbody.propagator import _G_KM3_KG_S2, _install_rails_forces

    ephem = Ephemeris("astropy")
    r0 = np.array([1.496e8, 2.0e7, 0.0])  # ~1 AU from Sun, ~1e7 km off the x-axis
    v0 = np.array([0.0, 29.78, 0.0])
    t0 = 0.0
    bodies = ("E",)

    a_ref = rails_acceleration(r0, t0, bodies, ephem)

    sim = rebound.Simulation()
    sim.G = _G_KM3_KG_S2
    sim.integrator = "ias15"
    sim.add(m=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
    sim.add(
        m=0.0,
        x=float(r0[0]),
        y=float(r0[1]),
        z=float(r0[2]),
        vx=float(v0[0]),
        vy=float(v0[1]),
        vz=float(v0[2]),
    )
    sim.t = t0
    cache = RailsEphemerisCache(bodies, ephem, t0, t0 + 200.0)
    _install_rails_forces(sim, bodies, cache)

    dt = 10.0
    sc = sim.particles[1]
    v_before = np.array([sc.vx, sc.vy, sc.vz])
    sim.integrate(t0 + dt)
    v_after = np.array([sc.vx, sc.vy, sc.vz])
    a_fd = (v_after - v_before) / dt

    # ~1e-6 rel residual is the trajectory curvature over the step, not a model
    # divergence: the two FORCE models agree, the clamp is inactive in far field.
    assert np.linalg.norm(a_fd - a_ref) / np.linalg.norm(a_ref) < 1e-4


def test_rails_cache_batch_samples_match_per_point() -> None:
    """Pin the I5 batching optimisation: the rails spline is built from the
    vectorised ephem.states() batch, whose samples must be BIT-IDENTICAL to the
    per-point ingest_planet_state() calls they replaced (asserted directly below
    via array_equal — the actual byte-for-byte parity claim).

    The spline-evaluation check on top of that uses a few-ULP tolerance, NOT
    exact reproduction, because CubicSpline knot fidelity is not exact at the
    FINAL knot (#803): PPoly's interval lookup is left-closed, so every interior
    knot returns the stored constant coefficient exactly (offset 0), but the
    last knot is evaluated as the last interval's cubic at offset h = 86400 s
    (scipy _ppoly ascending-power accumulation), which reproduces the endpoint
    sample only to floating-point rounding. Whether that rounds bit-exact is a
    knife edge in the spline coefficients (scipy's banded LAPACK solve — BLAS/
    Accelerate-build dependent): Mars' ~1.69e8 km x-coordinate lands 1 ULP
    (2**-25 km ≈ 3e-8) off on macOS/Accelerate. The old fixed 1e-9 km bound sat
    BELOW one ULP of any coordinate in [1.34e8, 2.68e8) km, i.e. it demanded
    exact endpoint rounding and passed only by luck. Honest bound: 4 ULPs of the
    largest coordinate (~1.2e-7 km at Mars) — still ~6 orders of magnitude
    tighter than the cache's real mid-interval interpolation accuracy vs DE440
    (~0.08 km), so the parity pin loses nothing physically meaningful.
    """
    pytest.importorskip("astropy")
    from cyclerfinder.core.constants import SECONDS_PER_DAY
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.nbody.forces import RailsEphemerisCache, ingest_planet_state

    ephem = Ephemeris("astropy")
    bodies = ("V", "E", "M")
    t0, t1 = 0.0, 60.0 * SECONDS_PER_DAY
    cache = RailsEphemerisCache(bodies, ephem, t0, t1)

    # Reconstruct the exact grid the cache builds (defaults: step 1 d, pad 5 d).
    pad, step = 5.0 * SECONDS_PER_DAY, 1.0 * SECONDS_PER_DAY
    lo, hi = min(t0, t1) - pad, max(t0, t1) + pad
    n = max(4, int((hi - lo) / step) + 1)
    grid = np.linspace(lo, hi, n)
    grid_list = [float(t) for t in grid]

    for body in bodies:
        batch = ephem.states([body] * len(grid_list), grid_list)
        for i, t in enumerate(grid_list):
            ref_r, ref_v = ingest_planet_state(body, t, ephem)
            # The I5 pin proper: batched states() ≡ per-point state(), bitwise.
            assert np.array_equal(batch[i][0], ref_r)
            assert np.array_equal(batch[i][1], ref_v)
            # Spline knot fidelity: exact at interior knots, ≤ a few ULPs of the
            # coordinate magnitude at the final knot (see docstring derivation).
            got = cache.position(body, t)
            assert np.max(np.abs(got - ref_r)) <= 4.0 * np.spacing(np.max(np.abs(ref_r)))


@pytest.mark.slow
def test_anchor_err_km_is_reserved_zero_contract() -> None:
    """RestrictedNBody seeds at the exact (r0, v0) and never anchors against an
    ephemeris, so anchor_err_km is a RESERVED slot that is always 0.0. Callers
    that need a closure error compute it themselves (rung.py does). Pin the
    contract so a future backend cannot silently start populating it unannounced.
    """
    pytest.importorskip("rebound")
    from cyclerfinder.nbody.propagator import RestrictedNBody

    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=10.0 * 86400.0, bodies=(), accuracy=1e-10
    )
    assert arc.anchor_err_km == 0.0
