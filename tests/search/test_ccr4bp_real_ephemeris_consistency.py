"""Tests for the CCR4BP-to-real-ephemeris consistency check (`#704`).

Structure: pure-math frame-conversion unit tests first (no SPICE needed),
then the MANDATORY positive control (`#704`'s own TDD requirement -- does
this module's real-N-body pipeline reduce EXACTLY to the idealized CCR4BP's
own rotating-frame propagation when fed CIRCULAR-COPLANAR-approximated moon
positions, the direct analogue of every other ``mu_gan=0``-style reduction
check in this whole `#689`-`#703` arc), then SPICE-gated integration tests
(skipped if the URA111 kernel is not installed, matching
``tests/data/test_v4_uranus_strict.py``'s own skip convention) that exercise
the module against `#701`'s own real, saved connection.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.ccr4bp_umbriel_titania as ut
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs
import cyclerfinder.search.ccr4bp_manifold_globalize as mg
import cyclerfinder.search.ccr4bp_real_ephemeris_consistency as rec
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.genome.composed_moon_map import resonance_semimajor

_KERNELS_PRESENT = (
    rec.DEFAULT_LSK_PATH.exists()
    and rec.DEFAULT_PCK_PATH.exists()
    and rec.DEFAULT_URA_PATH.exists()
)
_SKIP_REASON = (
    f"URA111 SPICE kernel not installed (looked at {rec.DEFAULT_URA_PATH}); "
    "see docs/notes or scripts/install_uranian_spice.sh"
)


# --------------------------------------------------------------------------- #
# Pure-math frame conversion
# --------------------------------------------------------------------------- #


def test_osculating_frame_circular_orbit_basis_and_rate() -> None:
    """Circular orbit: R's columns are orthonormal, x_hat points at the
    body, omega_inst matches the exact circular mean motion."""
    r_mag = 265_986.0
    v_mag = 4.0
    theta = 1.234
    r = np.array([r_mag * math.cos(theta), r_mag * math.sin(theta), 0.0])
    v = np.array([-v_mag * math.sin(theta), v_mag * math.cos(theta), 0.0])
    rot, omega_inst = rec.osculating_frame(r, v)
    assert np.allclose(rot.T @ rot, np.eye(3), atol=1e-12)
    assert np.allclose(rot[:, 0], r / r_mag, atol=1e-12)
    assert np.allclose(rot[:, 2], np.array([0.0, 0.0, 1.0]), atol=1e-12)
    assert math.isclose(omega_inst, v_mag / r_mag, rel_tol=1e-12)


def test_osculating_frame_matches_finite_difference_angle_rate_eccentric() -> None:
    """Eccentric two-body orbit: omega_inst = h/r^2 matches a direct
    finite-difference estimate of dtheta/dt from a short numerical
    propagation -- confirms the areal-velocity formula holds pointwise at
    nonzero eccentricity (not just the circular special case above)."""
    mu = 5.7945564e6  # GM_Uranus, km^3/s^2 (order-of-magnitude realistic)
    a = 265_986.0
    ecc = 0.30
    r_peri = a * (1.0 - ecc)
    v_peri = math.sqrt(mu * (2.0 / r_peri - 1.0 / a))
    state0 = np.array([r_peri, 0.0, 0.0, 0.0, v_peri, 0.0])

    def kepler_rhs(_t: float, y: np.ndarray) -> np.ndarray:
        r = y[:3]
        r_norm = np.linalg.norm(r)
        return np.concatenate([y[3:], -mu * r / r_norm**3])

    dt = 50.0
    sol = solve_ivp(kepler_rhs, (0.0, dt), state0, method="DOP853", rtol=1e-13, atol=1e-13)
    r0, v0 = state0[:3], state0[3:]
    r1 = sol.y[:3, -1]
    _rot0, omega0 = rec.osculating_frame(r0, v0)
    theta0 = math.atan2(r0[1], r0[0])
    theta1 = math.atan2(r1[1], r1[0])
    dtheta_dt_fd = (theta1 - theta0) / dt
    assert math.isclose(omega0, dtheta_dt_fd, rel_tol=2e-3)


def test_nondim_inertial_roundtrip() -> None:
    """nondim -> inertial -> nondim recovers the original state exactly
    (to floating-point precision) for an arbitrary moon state."""
    rng = np.random.default_rng(704)
    l_km = 265_986.0
    v_unit_km_s = 4.667
    r_moon = np.array([2.5e5, -1.2e5, 3.0e3])
    v_moon = np.array([1.1, 3.3, -0.02])
    state6 = rng.normal(size=6) * np.array([1.6, 0.2, 0.0, 0.03, 0.9, 0.0])
    r_km, v_km_s = rec.nondim_state_to_inertial(state6, r_moon, v_moon, l_km, v_unit_km_s)
    state6_back = rec.inertial_state_to_nondim(r_km, v_km_s, r_moon, v_moon, l_km, v_unit_km_s)
    assert np.allclose(state6, state6_back, atol=1e-9)


def test_tu_to_seconds_matches_umbriel_like_mean_motion() -> None:
    l_km = 265_986.0
    gm_uranus = 5.7945564e6
    gm_umbriel = 85.1
    n1 = math.sqrt((gm_uranus + gm_umbriel) / l_km**3)
    v_unit_km_s = l_km * n1
    tu_s = rec.tu_to_seconds(l_km, v_unit_km_s)
    assert math.isclose(tu_s, 1.0 / n1, rel_tol=1e-12)
    # Sanity: Umbriel's real orbital period is ~4.144 days; 2*pi TU should
    # land in that ballpark (this system's own base-moon mean motion).
    period_days = 2.0 * math.pi * tu_s / 86400.0
    assert 4.0 < period_days < 4.3


# --------------------------------------------------------------------------- #
# MANDATORY positive control: reduction to the idealized CCR4BP
# --------------------------------------------------------------------------- #


def test_reduces_to_idealized_ccr4bp_when_fed_circular_coplanar_positions() -> None:
    """`#704`'s own TDD requirement: feeding :func:`idealized_moon_state_fn`
    (circular-coplanar substitute, NOT real SPICE) into
    :func:`rec.propagate_real` must reproduce -- to tight integrator
    tolerance -- the SAME endpoint as directly propagating in the idealized
    CCR4BP's own rotating frame (:func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp`)
    and transforming the result to inertial via the identical frame
    construction. This is the ``mu_gan=0``-style reduction check every
    other module in this arc uses before trusting a NEW propagator on a
    novel comparison.
    """
    system = ut.uranus_umbriel_titania_default()
    l_km = ut.L_KM
    v_unit_km_s = ut.v_unit_km_s()
    n1 = v_unit_km_s / l_km

    # An arbitrary but fixed nondim departure state near the base resonant
    # orbit's own scale (not on the actual torus -- this positive control
    # only needs to exercise the frame/force-model machinery, not #701's
    # own torus).
    state0_nondim = np.array([1.6, 0.05, 0.0, 0.02, -0.9, 0.0])
    t_flow_tu = 3.7
    t_flow_s = t_flow_tu / n1

    # Path A: idealized rotating-frame propagation, then transform to
    # inertial using the SAME osculating-frame construction, evaluated at
    # Umbriel's IDEALIZED (circular) state at t=0 and at t=t_flow.
    idealized_fn = rec.idealized_moon_state_fn(
        system, l_km, v_unit_km_s, system.mu, system.theta_gan0, system.a_gan, system.omega_gan
    )
    r_umbriel0, v_umbriel0 = idealized_fn("Umbriel", 0.0)
    r0_km, v0_km_s = rec.nondim_state_to_inertial(
        state0_nondim, r_umbriel0, v_umbriel0, l_km, v_unit_km_s, mu=system.mu
    )

    arc = ccr4bp.propagate_ccr4bp(system, state0_nondim, t_flow_tu, t0=0.0, rtol=1e-13, atol=1e-13)
    r_umbriel_f, v_umbriel_f = idealized_fn("Umbriel", t_flow_s)
    r_f_expected_km, v_f_expected_km_s = rec.nondim_state_to_inertial(
        arc.state_f, r_umbriel_f, v_umbriel_f, l_km, v_unit_km_s, mu=system.mu
    )

    # Path B: this module's own real-N-body-shaped propagator, fed the
    # circular-coplanar substitute instead of real SPICE.
    prop = rec.propagate_real(
        r0_km,
        v0_km_s,
        t_flow_s,
        et0=0.0,
        mu_uranus=rec.PRIMARIES["Uranus"],
        perturber_moons=("Umbriel", "Titania"),
        moon_state_fn=idealized_fn,
        rtol=1e-13,
        atol=1e-9,
    )
    assert prop.success

    pos_gap_km = float(np.linalg.norm(prop.r_f_km - r_f_expected_km))
    vel_gap_km_s = float(np.linalg.norm(prop.v_f_km_s - v_f_expected_km_s))
    # Both paths integrate the SAME physics (central Uranus + Umbriel +
    # Titania, only Umbriel's own GM differs from the CCR4BP's exact
    # mass-ratio convention by float round-trip through SATELLITES -- see
    # module docstring); expect near-machine-precision agreement, not just
    # "small".
    assert pos_gap_km < 1e-3, f"reduction check pos_gap={pos_gap_km} km"
    assert vel_gap_km_s < 1e-8, f"reduction check vel_gap={vel_gap_km_s} km/s"


# --------------------------------------------------------------------------- #
# SPICE-gated integration tests
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_spice_moon_state_fn_umbriel_sma_sane() -> None:
    import spiceypy as spice

    from cyclerfinder.data.validation.v4_uranus_strict import (
        _ephemeris_time_seconds,
        _spice_furnsh_all,
    )

    spice.kclear()
    try:
        _spice_furnsh_all(
            (str(rec.DEFAULT_LSK_PATH), str(rec.DEFAULT_PCK_PATH), str(rec.DEFAULT_URA_PATH))
        )
        et = _ephemeris_time_seconds("2030-01-01T00:00:00")
        r, v = rec.spice_moon_state_fn("Umbriel", et)
    finally:
        spice.kclear()
    r_mag = float(np.linalg.norm(r))
    # Umbriel SMA 265,986 km, e=0.0041 -> r in [264,900, 267,072] roughly.
    assert 260_000.0 < r_mag < 272_000.0
    assert np.linalg.norm(v) > 0.0


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_check_connection_survives_real_ephemeris_e2e_runs_and_is_finite() -> None:
    """End-to-end smoke test against `#701`'s OWN actual saved connection
    (deterministically reconstructed the same way `#701`'s own driver
    script builds it), at one fixed epoch. Does not assert genuineness --
    only that the pipeline runs and produces finite, physically-scaled
    numbers (regression floor for the module actually working end to end).
    """
    torus = _rebuild_701_torus()
    departure_u, target_s, t_u_tu = _reconstruct_701_best_robust(torus)

    system = ut.uranus_umbriel_titania_default()
    result = rec.check_connection_survives_real_ephemeris(
        "2030-01-01T00:00:00",
        departure_u,
        target_s,
        t_u_tu,
        ut.L_KM,
        ut.v_unit_km_s(),
        system.mu,
    )
    assert result.propagation_success
    assert math.isfinite(result.pos_gap_km)
    assert math.isfinite(result.vel_gap_km_s)
    # Sanity scale: gap should be well within the system's own physical
    # extent (Titania's SMA, ~4.4e5 km) -- a total blowup would indicate a
    # wiring bug, not a genuine "collapses" physics finding (which would
    # still be bounded by escape/collision, not literally unbounded on a
    # ~13-day arc).
    assert result.pos_gap_km < 5.0e6


def _rebuild_701_torus() -> vt.CCR4BPTorusVariationalResult:
    system = ut.uranus_umbriel_titania_default()
    s0, period, res = _resonant_symmetric_orbit(system.mu, 1, 2)
    assert res < 1e-10
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )


def _reconstruct_701_best_robust(
    torus: vt.CCR4BPTorusVariationalResult,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Reproduce `#701`'s own saved ``best_robust_genuine_connection_corrected``
    departure/target states bit-for-bit from its own committed seed (see
    ``data/found/701_ccr4bp_umbriel_titania_search/result.json``)."""
    seed = hs.ManifoldCandidate(
        theta2_u=3.665191429188092,
        t_u=19.314531534610783,
        theta2_s=3.5604716740684323,
        t_s=18.187850528425155,
        gap_planar=0.0011381326404778097,
    )
    refined = hs.refine_candidate(
        torus, torus, seed, lobe_sign_u=-1.0, lobe_sign_s=-1.0, n_segments_dir=32
    )
    assert refined is not None
    assert math.isclose(refined.residual_norm, 1.1159187446079244e-14, rel_tol=1e-3)
    departure_u = mg.manifold_state_at(
        torus,
        "unstable",
        0.0,
        refined.theta2_u,
        0.0,
        lobe_sign=-1.0,
        ref_vec=refined.ref_vec_u,
        n_segments_dir=32,
    )
    assert departure_u is not None
    return departure_u, refined.state_s, refined.t_u


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to `#701`'s own driver script."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res
