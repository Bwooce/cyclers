"""Tests for the impulsive primer-vector optimality diagnostic.

Golden gates are SOURCED (golden-test discipline): expected verdicts trace to
published primer-vector theory (Lawden; Lion & Handelsman 1968; Prussing &
Conway, *Orbital Mechanics* / *Spacecraft Trajectory Optimization* Ch. 2), not
to magnitudes our own code produced.

Sourced facts used as gates
---------------------------
* The coplanar-circular Hohmann transfer satisfies the primer NECESSARY
  conditions (``|p| <= 1`` on the transfer coast) — the endpoints are the
  local maxima at exactly unity (Prussing & Conway). This holds for the
  symmetric 180-degree Hohmann transfer at every radius ratio (verified below
  at ratio 2 and ratio 20). NB: the published bi-elliptic-vs-Hohmann ΔV
  threshold (~11.94) is a *global* / coast-extension comparison and does NOT
  appear as an interior ``|p| > 1`` bulge on the symmetric Hohmann transfer
  coast — see the results note. We therefore do NOT (and must not) assert an
  interior violation on a Hohmann coast that has none.
* A fixed-time transfer that takes the "long way" (transfer angle > 180°) is a
  classic NON-optimal trajectory: Lion & Handelsman's diagnostic shows ``|p|``
  bulges above unity in the interior, signalling that an added/relocated
  impulse lowers total ΔV (IMPROVABLE). We assert only the qualitative
  ``max|p| > 1`` necessary-condition violation, never an invented magnitude.

Construction invariants (no sourced magnitude needed)
-----------------------------------------------------
* ``|p(t_i)| = 1`` at each impulse (boundary conditions).
* Continuity of ``p`` across a coast endpoint.
* ``G(r)`` symmetry and tracelessness.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.verify.primer import (
    CoastPrimerResult,
    PrimerDiagnostic,
    PrimerVerdict,
    diagnose_impulse_schedule,
    gravity_gradient,
    hohmann_primer_diagnostic,
    primer_on_coast,
)

MU_CANONICAL = 1.0


# ---------------------------------------------------------------------------
# G(r) algebraic invariants
# ---------------------------------------------------------------------------


def test_gravity_gradient_symmetric_and_traceless() -> None:
    """G = (μ/r³)(3 r̂r̂ᵀ - I) is symmetric and traceless by construction."""
    for r in (
        np.array([1.0, 0.0, 0.0]),
        np.array([3.0, -1.0, 2.0]),
        np.array([0.0, 0.0, 5.0]),
    ):
        g = gravity_gradient(r, mu=MU_CANONICAL)
        assert np.allclose(g, g.T, atol=1e-14)
        assert abs(np.trace(g)) < 1e-12


def test_gravity_gradient_radial_eigenvalue() -> None:
    """Along r̂ the eigenvalue is +2μ/r³; transverse directions give -μ/r³."""
    r = np.array([2.0, 0.0, 0.0])
    g = gravity_gradient(r, mu=MU_CANONICAL)
    r_n = 2.0
    assert g[0, 0] == pytest.approx(2.0 * MU_CANONICAL / r_n**3)
    assert g[1, 1] == pytest.approx(-MU_CANONICAL / r_n**3)
    assert g[2, 2] == pytest.approx(-MU_CANONICAL / r_n**3)


# ---------------------------------------------------------------------------
# Boundary-condition / construction invariants
# ---------------------------------------------------------------------------


def _circular_state(r: float, mu: float = MU_CANONICAL) -> tuple[np.ndarray, np.ndarray]:
    return np.array([r, 0.0, 0.0]), np.array([0.0, np.sqrt(mu / r), 0.0])


def test_endpoint_magnitudes_unit_by_construction() -> None:
    """|p(t_i)| = 1 at both impulses bounding a coast (Lawden BC)."""
    r0, v0 = _circular_state(1.0)
    # Mild eccentric coast over a quarter period.
    v0 = v0 * 1.1
    a = 1.0 / (2.0 / 1.0 - (1.1 * np.sqrt(MU_CANONICAL)) ** 2 / MU_CANONICAL)
    tof = 0.25 * 2.0 * np.pi * np.sqrt(a**3 / MU_CANONICAL)
    p0 = np.array([0.3, 0.95, 0.0])
    p1 = np.array([-0.6, 0.8, 0.0])
    res = primer_on_coast(r0, v0, p0, p1, tof, mu=MU_CANONICAL, n_samples=200)
    assert res.endpoint_magnitudes[0] == pytest.approx(1.0, abs=1e-9)
    assert res.endpoint_magnitudes[1] == pytest.approx(1.0, abs=1e-9)


def test_primer_continuity_grid_refinement_invariant() -> None:
    """The recovered max|p| is grid-independent (continuity / convergence).

    The primer is a smooth function of time on a coast; the diagnostic samples
    it on a grid. A coarse and a fine grid must agree on max|p| to high
    precision, confirming the STM reconstruction is continuous and the sampling
    converged (no spurious jumps between adjacent grid points).
    """
    r0, v0 = _circular_state(1.0)
    v0 = v0 * 1.05
    a = 1.0 / (2.0 / 1.0 - (1.05 * np.sqrt(MU_CANONICAL)) ** 2 / MU_CANONICAL)
    tof = 0.4 * 2.0 * np.pi * np.sqrt(a**3 / MU_CANONICAL)
    p0 = np.array([0.0, 1.0, 0.0])
    p1 = np.array([0.2, 0.97, 0.0])
    coarse = primer_on_coast(r0, v0, p0, p1, tof, mu=MU_CANONICAL, n_samples=201)
    fine = primer_on_coast(r0, v0, p0, p1, tof, mu=MU_CANONICAL, n_samples=801)
    assert fine.max_primer_magnitude == pytest.approx(coarse.max_primer_magnitude, rel=2e-3)
    # Both endpoints remain pinned to unity regardless of grid.
    assert fine.endpoint_magnitudes[0] == pytest.approx(1.0, abs=1e-9)
    assert fine.endpoint_magnitudes[1] == pytest.approx(1.0, abs=1e-9)


# ---------------------------------------------------------------------------
# Sourced golden gate 1 — Hohmann satisfies the necessary conditions
# ---------------------------------------------------------------------------


def test_hohmann_ratio_2_optimal() -> None:
    """Hohmann r2/r1 = 2: |p| ≤ 1 on the coast → OPTIMAL (Prussing & Conway).

    Below the ~11.94 bi-elliptic threshold the Hohmann transfer is ΔV-optimal
    and its primer satisfies the necessary conditions.
    """
    diag = hohmann_primer_diagnostic(1.0, 2.0, mu=MU_CANONICAL)
    assert diag.overall_verdict is PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET
    assert diag.max_primer_magnitude <= 1.0 + 1e-6
    (coast,) = diag.coasts
    assert coast.endpoint_magnitudes[0] == pytest.approx(1.0, abs=1e-9)
    assert coast.endpoint_magnitudes[1] == pytest.approx(1.0, abs=1e-9)


def test_hohmann_ratio_20_coast_satisfies_necessary_conditions() -> None:
    """Hohmann r2/r1 = 20: the transfer COAST still has |p| ≤ 1.

    The symmetric (180°) Hohmann primer touches unity only at its endpoints at
    every radius ratio (Prussing & Conway). The bi-elliptic superiority above
    ~11.94 is a global / coast-extension comparison, NOT an interior bulge on
    the Hohmann transfer coast — so the per-coast necessary-conditions check
    correctly returns OPTIMAL here. (See the results note for why this does not
    contradict the published bi-elliptic threshold.)
    """
    diag = hohmann_primer_diagnostic(1.0, 20.0, mu=MU_CANONICAL)
    assert diag.max_primer_magnitude <= 1.0 + 1e-6
    assert diag.overall_verdict is PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET
    # The interior genuinely dips below unity (endpoints are the maxima).
    assert diag.coasts[0].time_of_max_s in (
        pytest.approx(0.0, abs=1e-6),
        pytest.approx(diag.coasts[0].duration_s, rel=1e-6),
    )


# ---------------------------------------------------------------------------
# Sourced golden gate 2 — a non-optimal long-way transfer is IMPROVABLE
# ---------------------------------------------------------------------------


def _long_way_transfer(
    r1: float, r2: float, theta_deg: float, tof_factor: float, mu: float = MU_CANONICAL
) -> PrimerDiagnostic:
    """Build a fixed-time transfer with transfer angle > 180° and diagnose it.

    Departs a circular orbit of radius ``r1`` at angle 0, arrives a circular
    orbit of radius ``r2`` at ``theta_deg`` (> 180°, the long way), over a time
    ``tof_factor x t_Hohmann``. The primer is pinned to the actual departure /
    arrival ΔV directions (Lambert solution minus the local circular velocity).
    Such transfers are classic non-optimal cases (Lion & Handelsman): an
    interior impulse lowers the cost.
    """
    th = np.radians(theta_deg)
    big_r1 = np.array([r1, 0.0, 0.0])
    big_r2 = np.array([r2 * np.cos(th), r2 * np.sin(th), 0.0])
    a_t = 0.5 * (r1 + r2)
    t_hohmann = np.pi * np.sqrt(a_t**3 / mu)
    tof = tof_factor * t_hohmann
    sol = lambert(big_r1, big_r2, tof, mu=mu)[0]
    v_dep = np.asarray(sol.v1, dtype=np.float64)
    v_arr = np.asarray(sol.v2, dtype=np.float64)
    vc1 = np.sqrt(mu / r1) * np.array([0.0, 1.0, 0.0])
    vc2 = np.sqrt(mu / r2) * np.array([-np.sin(th), np.cos(th), 0.0])
    dv1 = v_dep - vc1
    dv2 = vc2 - v_arr
    return diagnose_impulse_schedule(
        [(big_r1, v_dep, tof)],
        [dv1 / np.linalg.norm(dv1), dv2 / np.linalg.norm(dv2)],
        mu=mu,
        n_samples=600,
    )


def test_long_way_transfer_improvable() -> None:
    """A 200°, fixed-time transfer is non-optimal: max|p| > 1 (IMPROVABLE).

    Lion & Handelsman necessary-condition violation: the interior |p| exceeds
    unity, so an added/relocated impulse reduces total ΔV. Only the qualitative
    verdict is asserted (no invented magnitude).
    """
    diag = _long_way_transfer(1.0, 4.0, theta_deg=200.0, tof_factor=1.0)
    assert diag.overall_verdict is PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    assert diag.max_primer_magnitude > 1.0
    # The peak is in the interior (not at an endpoint), which is the
    # signature of a beneficial midcourse impulse.
    coast = diag.coasts[0]
    assert 0.0 < coast.time_of_max_s < coast.duration_s


def test_long_way_transfer_peak_strictly_interior() -> None:
    """A 250°, longer-time transfer bulges further above unity, interior peak."""
    diag = _long_way_transfer(1.0, 4.0, theta_deg=250.0, tof_factor=1.5)
    assert diag.overall_verdict is PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    assert diag.max_primer_magnitude > 1.5
    coast = diag.coasts[0]
    frac = coast.time_of_max_s / coast.duration_s
    assert 0.05 < frac < 0.95


# ---------------------------------------------------------------------------
# Multi-coast schedule plumbing
# ---------------------------------------------------------------------------


def test_multi_coast_schedule_aggregates_verdict() -> None:
    """A schedule with one optimal and one improvable coast reports IMPROVABLE."""
    # Coast A: Hohmann-like optimal arc (reuse the ratio-2 transfer geometry).
    r1, r2 = 1.0, 2.0
    a_t = 0.5 * (r1 + r2)
    r0a = np.array([r1, 0.0, 0.0])
    v0a = np.array([0.0, np.sqrt(MU_CANONICAL * (2.0 / r1 - 1.0 / a_t)), 0.0])
    tofa = np.pi * np.sqrt(a_t**3 / MU_CANONICAL)
    # Coast B: long-way arc that bulges.
    th = np.radians(200.0)
    r0b = np.array([1.0, 0.0, 0.0])
    big_r2 = np.array([4.0 * np.cos(th), 4.0 * np.sin(th), 0.0])
    ab = 0.5 * (1.0 + 4.0)
    tofb = np.pi * np.sqrt(ab**3 / MU_CANONICAL)
    sol = lambert(r0b, big_r2, tofb, mu=MU_CANONICAL)[0]
    v0b = np.asarray(sol.v1, dtype=np.float64)

    # Directions: A endpoints tangential; B endpoints from the Lambert ΔV.
    p_a0 = np.array([0.0, 1.0, 0.0])
    p_a1 = np.array([0.0, -1.0, 0.0])
    vc2 = np.sqrt(MU_CANONICAL / 4.0) * np.array([-np.sin(th), np.cos(th), 0.0])
    dv_b = vc2 - np.asarray(sol.v2, dtype=np.float64)
    p_b1 = dv_b / np.linalg.norm(dv_b)

    diag = diagnose_impulse_schedule(
        [(r0a, v0a, tofa), (r0b, v0b, tofb)],
        [p_a0, p_a1, p_b1],
        mu=MU_CANONICAL,
        n_samples=400,
    )
    assert len(diag.coasts) == 2
    assert diag.overall_verdict is PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    assert isinstance(diag.coasts[0], CoastPrimerResult)
    assert diag.caveat  # non-empty caveat is always present


def test_schedule_rejects_wrong_impulse_count() -> None:
    r0, v0 = _circular_state(1.0)
    with pytest.raises(ValueError, match="impulse directions"):
        diagnose_impulse_schedule(
            [(r0, v0 * 1.05, 1.0)],
            [np.array([0.0, 1.0, 0.0])],  # only 1 dir for 1 coast (need 2)
            mu=MU_CANONICAL,
        )


def test_coast_stm_matches_kepler_reference() -> None:
    """The STM integrator's reference arc matches the project Kepler propagator."""
    from cyclerfinder.verify.primer import _coast_stm

    r0, v0 = _circular_state(1.0)
    v0 = v0 * 1.1
    a = 1.0 / (2.0 / 1.0 - (1.1 * np.sqrt(MU_CANONICAL)) ** 2 / MU_CANONICAL)
    tof = 0.3 * 2.0 * np.pi * np.sqrt(a**3 / MU_CANONICAL)
    _times, ref, _stms = _coast_stm(r0, v0, tof, MU_CANONICAL, n_samples=50)
    r_end_kepler, v_end_kepler = propagate(r0, v0, tof, MU_CANONICAL)
    assert np.allclose(ref[-1, :3], r_end_kepler, rtol=1e-7, atol=1e-7)
    assert np.allclose(ref[-1, 3:], v_end_kepler, rtol=1e-7, atol=1e-7)


# ---------------------------------------------------------------------------
# Aldrin maintenance-schedule application (DIAGNOSTIC, PROVISIONAL)
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_aldrin_maintenance_primer_diagnostic() -> None:
    """Primer diagnostic on the real in-family Aldrin E-M-E maintenance schedule.

    DIAGNOSTIC / PROVISIONAL: the maintenance ΔV is our own value (not sourced),
    and long multi-rev cycler legs are exactly where linearised primer theory is
    flagged as fragile (Guzman 2002, not yet acquired). This test pins only that
    the diagnostic RUNS on the real schedule and returns a verdict; it does not
    assert an optimality claim as a golden.
    """
    from datetime import UTC, datetime

    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv

    ephem = Ephemeris("astropy")
    result = optimise_aldrin_maintenance_dv(
        ephem,
        real_window_priority_date=datetime(1985, 10, 28, tzinfo=UTC),
    )
    assert result.converged, "Aldrin maintenance solve did not converge"

    cycler = result.cycler
    # Build the impulse schedule from the cycler's heliocentric legs: each leg
    # is a coast; at each encounter the spacecraft V∞ is bent (the ΔV).
    coast_states: list[tuple[np.ndarray, np.ndarray, float]] = []
    for i, (leg, enc) in enumerate(zip(cycler.legs, cycler.encounters[:-1], strict=True)):
        dur = float(cycler.encounters[i + 1].t - enc.t)
        coast_states.append(
            (np.asarray(enc.r, dtype=np.float64), np.asarray(leg.v_depart, dtype=np.float64), dur)
        )
    # Impulse directions: the heliocentric ΔV at each encounter (v_depart of the
    # next leg minus v_arrive of the previous), with the closure ΔV at the ends.
    dirs: list[np.ndarray] = []
    n_enc = len(cycler.encounters)
    for i, enc in enumerate(cycler.encounters):
        if i == 0:
            dv = cycler.legs[0].v_depart - (enc.v_planet + enc.vinf_in)
        elif i == n_enc - 1:
            dv = (enc.v_planet + enc.vinf_out) - cycler.legs[-1].v_arrive
        else:
            dv = cycler.legs[i].v_depart - cycler.legs[i - 1].v_arrive
        norm = float(np.linalg.norm(dv))
        # Degenerate (near-ballistic) impulse: fall back to the V∞-bend direction.
        if norm < 1e-9:
            dv = np.asarray(enc.vinf_out - enc.vinf_in, dtype=np.float64)
            norm = float(np.linalg.norm(dv))
        if norm < 1e-12:
            dv = np.array([1.0, 0.0, 0.0])
            norm = 1.0
        dirs.append(np.asarray(dv, dtype=np.float64) / norm)

    diag = diagnose_impulse_schedule(coast_states, dirs, n_samples=300)
    assert diag.overall_verdict in (
        PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET,
        PrimerVerdict.IMPROVABLE_ADD_IMPULSE,
    )
    assert diag.caveat
    # Surface the per-coast diagnostic for the results note (visible with -s).
    print("\nAldrin maintenance primer diagnostic (DIAGNOSTIC, PROVISIONAL):")
    print(f"  maintenance_dv_kms (our value) = {result.maintenance_dv_kms:.4f}")
    for c in diag.coasts:
        print(
            f"  coast {c.coast_index}: max|p|={c.max_primer_magnitude:.4f} "
            f"at t={c.time_of_max_s / 86400.0:.1f} d / {c.duration_s / 86400.0:.1f} d "
            f"-> {c.verdict.value}"
        )
    print(f"  overall: {diag.overall_verdict.value}")
