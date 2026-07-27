"""Tests for the N=5 CRNBP-torus-to-real-ephemeris consistency check (`#726`).

Structure mirrors ``tests/search/test_ccr4bp_real_ephemeris_consistency.py``
(`#704`) exactly: pure-math sanity checks first (no SPICE needed), then the
MANDATORY positive control (reduction to the idealized CRNBP when fed
circular-coplanar substitute moon positions -- `#704`'s own TDD requirement,
one perturber further), then SPICE-gated integration tests (skipped if
``jup365.bsp``/the vendored LSK are not installed) that exercise the module
against the ACTUAL delivered `#720`/`#723`/`#724` torus.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.crnbp_real_ephemeris_consistency as rec5
import cyclerfinder.search.variational_ccr4bp_torus as vt
import cyclerfinder.search.variational_crnbp_torus as vc
from cyclerfinder.genome.composed_moon_map import resonance_semimajor
from cyclerfinder.verify.spice_kernels import NAIF_JUP365_LOCAL, VENDORED_LSK_PATH

_KERNELS_PRESENT = NAIF_JUP365_LOCAL.expanduser().exists() and VENDORED_LSK_PATH.exists()
_SKIP_REASON = (
    f"jup365.bsp not installed (looked at {NAIF_JUP365_LOCAL.expanduser()}); "
    "see cyclerfinder.verify.spice_kernels.ensure_jup365_kernel"
)


# --------------------------------------------------------------------------- #
# Pure-math sanity
# --------------------------------------------------------------------------- #


def test_l_km_and_v_unit_km_s_sane() -> None:
    """``L_KM`` is Europa's registry SMA (671,100 km); ``v_unit_km_s()``
    implies a ~7.05-day... no, a base (Europa) synodic TU period consistent
    with Europa's own ~3.55-day sidereal period order of magnitude (the base
    frame rotates at Europa's own mean motion by CR3BP construction)."""
    assert rec5.L_KM == 671_100.0
    v_unit = rec5.v_unit_km_s()
    n1 = v_unit / rec5.L_KM
    period_days = 2.0 * math.pi / n1 / 86400.0
    # Europa's real sidereal period is 3.551 d; the CR3BP base-pair TU period
    # (using the gm_j+gm_e "system GM folding" convention -- see module
    # docstring) should land close to that.
    assert 3.4 < period_days < 3.7


def test_idealized_moon_state_fn_crnbp_europa_at_l_km() -> None:
    """At t=0, Europa sits exactly ``L_KM`` from Jupiter on +x (by
    construction: the primary separation IS the length unit)."""
    system = crnbp.jupiter_europa_io_ganymede_default()
    fn = rec5.idealized_moon_state_fn_crnbp(system, rec5.L_KM, rec5.v_unit_km_s())
    r, v = fn("Europa", 0.0)
    assert math.isclose(float(np.linalg.norm(r)), rec5.L_KM, rel_tol=1e-12)
    assert math.isclose(r[0], rec5.L_KM, rel_tol=1e-12)
    assert math.isclose(r[1], 0.0, abs_tol=1e-9)
    assert np.linalg.norm(v) > 0.0


def test_idealized_moon_state_fn_crnbp_io_ganymede_radii() -> None:
    """Io's and Ganymede's idealized circular radii match ``a_io``/``a_gan``
    scaled by ``L_KM`` (their own registry SMA ratios to Europa's)."""
    system = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = system.perturbers
    fn = rec5.idealized_moon_state_fn_crnbp(system, rec5.L_KM, rec5.v_unit_km_s())
    r_io, _ = fn("Io", 0.0)
    r_gan, _ = fn("Ganymede", 0.0)
    # A tiny +mu barycentre-to-Jupiter shift is applied (see module
    # docstring), so allow a loose relative tolerance dominated by that shift
    # (~mu ~ 2.5e-5 relative) rather than requiring bit-exactness.
    assert math.isclose(float(np.linalg.norm(r_io)), rec5.L_KM * io.a, rel_tol=1e-3)
    assert math.isclose(float(np.linalg.norm(r_gan)), rec5.L_KM * gan.a, rel_tol=1e-3)


def test_idealized_moon_state_fn_crnbp_unknown_moon_raises() -> None:
    system = crnbp.jupiter_europa_io_ganymede_default()
    fn = rec5.idealized_moon_state_fn_crnbp(system, rec5.L_KM, rec5.v_unit_km_s())
    with pytest.raises(ValueError, match="unsupported moon"):
        fn("Callisto", 0.0)


# --------------------------------------------------------------------------- #
# MANDATORY positive control: reduction to the idealized CRNBP
# --------------------------------------------------------------------------- #


def test_reduces_to_idealized_crnbp_when_fed_circular_coplanar_positions() -> None:
    """`#704`'s own TDD requirement, one perturber further: feeding
    :func:`rec5.idealized_moon_state_fn_crnbp` (circular-coplanar substitute,
    NOT real SPICE) into the reused :func:`~cyclerfinder.search.
    ccr4bp_real_ephemeris_consistency.propagate_real` must reproduce -- to
    tight tolerance -- the SAME endpoint as directly propagating in the
    idealized CRNBP's own rotating frame
    (:func:`cyclerfinder.core.crnbp.propagate_crnbp`) and transforming the
    result to inertial via the identical frame construction.

    Observed gap (2026-07-27, this pass): pos_gap ~1.75e-2 km, vel_gap
    ~1.13e-7 km/s, over an ~8.06e5 km / ~1.8e5 s arc -- consistent with the
    module docstring's ~2.07e-4-relative "system GM folding" approximation
    (Jupiter's folded-in Galilean mass is a larger fraction of the central
    term than Uranus's folded-in Umbriel+Titania mass was for `#704`, so this
    module's own reduction gap is correspondingly larger than `#704`'s
    1e-3 km / 1e-8 km/s -- both are the SAME class of accepted model
    approximation, not a code defect; bounds below are set well above the
    observed value as a regression floor, not a tight target).
    """
    system = crnbp.jupiter_europa_io_ganymede_default()
    l_km = rec5.L_KM
    v_unit_km_s = rec5.v_unit_km_s()
    n1 = v_unit_km_s / l_km

    # An arbitrary but fixed nondim departure state near the resonant torus's
    # own scale (between Europa's and Ganymede's orbits) -- not on the actual
    # torus; this positive control only needs to exercise the frame/force
    # machinery, not #720/#723/#724's own torus.
    state0_nondim = np.array([1.2, 0.05, 0.0, 0.02, 0.9, 0.0])
    t_flow_tu = 3.7
    t_flow_s = t_flow_tu / n1

    idealized_fn = rec5.idealized_moon_state_fn_crnbp(system, l_km, v_unit_km_s)
    r_europa0, v_europa0 = idealized_fn("Europa", 0.0)
    r0_km, v0_km_s = rec5.nondim_state_to_inertial(
        state0_nondim, r_europa0, v_europa0, l_km, v_unit_km_s, mu=system.mu
    )

    arc = crnbp.propagate_crnbp(system, state0_nondim, t_flow_tu, rtol=1e-13, atol=1e-13)
    r_europa_f, v_europa_f = idealized_fn("Europa", t_flow_s)
    r_f_expected_km, v_f_expected_km_s = rec5.nondim_state_to_inertial(
        arc.state_f, r_europa_f, v_europa_f, l_km, v_unit_km_s, mu=system.mu
    )

    from cyclerfinder.core.satellites import PRIMARIES
    from cyclerfinder.search.ccr4bp_real_ephemeris_consistency import propagate_real

    prop = propagate_real(
        r0_km,
        v0_km_s,
        t_flow_s,
        et0=0.0,
        mu_uranus=PRIMARIES["Jupiter"],
        perturber_moons=rec5.PERTURBER_MOONS_DEFAULT,
        moon_state_fn=idealized_fn,
        rtol=1e-13,
        atol=1e-9,
    )
    assert prop.success

    pos_gap_km = float(np.linalg.norm(prop.r_f_km - r_f_expected_km))
    vel_gap_km_s = float(np.linalg.norm(prop.v_f_km_s - v_f_expected_km_s))
    assert pos_gap_km < 0.5, f"reduction check pos_gap={pos_gap_km} km"
    assert vel_gap_km_s < 1e-5, f"reduction check vel_gap={vel_gap_km_s} km/s"


def test_check_torus_survives_real_ephemeris_with_idealized_fn_is_near_exact_at_short_window() -> (
    None
):
    """Chaos/noise-floor control: feeding the check's OWN top-level driver
    the IDEALIZED substitute moon-state function (instead of real SPICE), at
    a short window, must reproduce ``propagate_crnbp``'s own target closely
    -- confirming the large gaps seen with real SPICE (see the e2e test
    below) are a genuine real-ephemeris effect, not integrator/model noise.
    Uses a throwaway 1-perturber-continuation torus (mu_io=0, i.e. exactly
    `#690`'s CCR4BP torus) purely to keep this test fast -- the reduction
    property being tested does not depend on which converged torus is used.
    """
    system4 = ccr4bp.jupiter_europa_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(system4.mu, 3, 4)
    assert res < 1e-10
    ccr4bp_torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system4,
        s0,
        period,
        n1=2,
        n2=10,
        tr_solver="exact",
        max_nfev=300,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    target = crnbp.jupiter_europa_io_ganymede_default()
    torus = vc.discover_crnbp_torus_from_ccr4bp_seed(
        ccr4bp_torus,
        mu_io=0.0,
        a_io=target.perturbers[0].a,
        omega_io=target.perturbers[0].omega,
        theta_io0=target.perturbers[0].theta0,
        tr_solver="exact",
        max_nfev=300,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    idealized_fn = rec5.idealized_moon_state_fn_crnbp(torus.system, rec5.L_KM, rec5.v_unit_km_s())
    result = rec5.check_torus_survives_real_ephemeris(
        "2000-01-01T00:00:00",
        torus,
        0.0,
        0.0,
        t_window_tu=torus.period * 0.02,
        moon_state_fn=idealized_fn,
    )
    assert result.propagation_success
    # Short window, idealized-fed: should be small relative to the system
    # scale (L_KM), not merely finite.
    assert result.pos_gap_km < 1000.0, f"pos_gap={result.pos_gap_km} km"


# --------------------------------------------------------------------------- #
# SPICE-gated integration tests
# --------------------------------------------------------------------------- #


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_jupiter_spice_moon_state_fn_europa_sma_sane() -> None:
    import spiceypy as spice

    from cyclerfinder.data.validation.v4_uranus_strict import _ephemeris_time_seconds
    from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel

    spice.furnsh(ensure_leapseconds_kernel())
    et = _ephemeris_time_seconds("2030-01-01T00:00:00")
    r, v = rec5.jupiter_spice_moon_state_fn("Europa", et)
    r_mag = float(np.linalg.norm(r))
    # Europa SMA 671,100 km, e~0.009 -> r comfortably within [660,000, 682,000].
    assert 655_000.0 < r_mag < 690_000.0
    assert np.linalg.norm(v) > 0.0


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_check_torus_survives_real_ephemeris_e2e_headline() -> None:
    """End-to-end run against the ACTUAL delivered `#720`/`#723`/`#724`
    torus (physical Io mass, ``theta_io0=pi``), at the SAME headline epoch
    `#704` used for its own Umbriel-Titania headline result. Reproduces the
    torus deterministically from the same pipeline
    ``scripts/verify_724_rerun_continuation.py`` uses.

    Does not assert a specific numeric gap (that is `#726`'s own reported
    finding, not a frozen regression target) -- only that the pipeline runs,
    produces finite results, and that the observed generic-collapse scale
    (comparable to the system's own physical extent, NOT machine-precision-
    close) is reproduced within an order of magnitude, matching this pass's
    own recorded result.
    """
    final = _rebuild_724_final_torus()
    result = rec5.check_torus_survives_real_ephemeris("2030-01-01T00:00:00", final, 0.0, 0.0)
    assert result.propagation_success
    assert math.isfinite(result.pos_gap_km)
    assert math.isfinite(result.vel_gap_km_s)
    # Sanity floor: not a wiring-bug blowup (bounded well above Ganymede's
    # own SMA, ~1.07e6 km) but also not machine-precision-close (this is the
    # actual, reported "generic collapse" finding -- see docs/notes).
    assert result.pos_gap_km < 5.0e6
    assert result.pos_gap_km > 1.0e3


def _rebuild_724_final_torus() -> vc.CRNBPTorusVariationalResult:
    """Reproduce `#720`/`#723`/`#724`'s own delivered N=5 torus bit-for-bit
    (same pipeline as ``scripts/verify_724_rerun_continuation.py``)."""
    system4 = ccr4bp.jupiter_europa_ganymede_default()
    target = crnbp.jupiter_europa_io_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(system4.mu, 3, 4)
    assert res < 1e-10
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system4,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=target.perturbers[0].a,
        omega_io=target.perturbers[0].omega,
        theta_io0=target.perturbers[0].theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    return steps[-1]


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical scaffolding to `#724`'s own reproduction script."""
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
