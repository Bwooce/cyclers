"""Tests for the #648 deflated seedless-corrector family enumerator.

Combines `deflated_newton.py` (#524, previously only ever aimed at basin-
restricted shooting residuals) with the #606 seedless spectral corrector's
own Fourier-coefficient residual, using a gauge-invariant (phase-minimized
via FFT cross-correlation) distance metric between candidate solutions plus
a mandatory Radau cross-check on every deflated find. See
`src/cyclerfinder/search/deflated_variational_periodic_orbit.py`'s module
docstring for the full design rationale.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.deflated_variational_periodic_orbit import (
    _magnitude_fingerprint,
    enumerate_families_fixed_jacobi,
    gauge_distance,
    gauge_distance_z,
    phase_shift_coeffs,
    same_family,
)
from cyclerfinder.search.variational_periodic_orbit import (
    _n_free_actual,
    _unpack,
    discover_periodic_orbit_fixed_jacobi,
)

# --------------------------------------------------------------------------
# Gauge-invariant distance metric -- pure algebra, no solver noise involved.
# --------------------------------------------------------------------------

_Loop = tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]


def _random_loop(rng: np.random.Generator, n_harmonics: int = 8) -> _Loop:
    dc = np.array([0.8, 0.0, 0.02])
    cosc = rng.normal(scale=0.1, size=(3, n_harmonics))
    sinc = rng.normal(scale=0.1, size=(3, n_harmonics))
    period = 2.7
    return dc, cosc, sinc, period


def _pack_free_anchor_z(
    dc: NDArray[np.float64], cosc: NDArray[np.float64], sinc: NDArray[np.float64], period: float
) -> NDArray[np.float64]:
    """Pack ``(dc, cosc, sinc, period)`` into the raw free-variable vector
    matching :func:`~cyclerfinder.search.variational_periodic_orbit._unpack`'s
    layout with ALL THREE amplitude anchors free (``anchor_x1=anchor_y1=
    anchor_z1=None``, the pattern used throughout the #648 deflated
    enumerator) -- the inverse of ``_unpack`` for that one fixed anchor
    pattern, used here to build synthetic test fixtures directly."""
    parts = [dc[0], dc[1], dc[2], cosc[0, 0]]
    parts.extend(cosc[0, 1:])
    parts.extend(sinc[0, 1:])
    parts.append(sinc[1, 0])
    parts.append(cosc[1, 0])
    parts.extend(cosc[1, 1:])
    parts.extend(sinc[1, 1:])
    parts.append(cosc[2, 0])
    parts.append(sinc[2, 0])
    parts.extend(cosc[2, 1:])
    parts.extend(sinc[2, 1:])
    parts.append(np.log(period))
    return np.array(parts, dtype=np.float64)


def test_gauge_distance_self_consistency_phase_shift() -> None:
    """The load-bearing correctness test (#648 step 4): a loop compared
    against its OWN Fourier representation, phase-shifted by a range of
    delta (including near-2pi wraparound), must register as the SAME orbit
    at near-zero distance -- not a numerically-coincidental single point.
    """
    rng = np.random.default_rng(0)
    dc, cosc, sinc, period = _random_loop(rng)
    for delta in [0.0, 0.3, 1.7, -2.2, 3.14159, 6.0]:
        cosc2, sinc2 = phase_shift_coeffs(cosc, sinc, delta)
        dist, info = gauge_distance(dc, cosc, sinc, period, dc, cosc2, sinc2, period, n_samples=256)
        assert dist < 1e-6, f"delta={delta}: expected near-zero, got {dist}"
        assert info["period_rel_diff"] == 0.0


def test_gauge_distance_different_shapes_far_apart() -> None:
    """Two genuinely different loops (independent random coefficients, same
    period) must NOT collapse to a small distance -- the metric has real
    discriminating power, not just a trivial always-zero degenerate case.
    """
    rng = np.random.default_rng(1)
    dc, cosc_a, sinc_a, period = _random_loop(rng)
    _dc2, cosc_b, sinc_b, _p2 = _random_loop(rng)
    dist, info = gauge_distance(dc, cosc_a, sinc_a, period, dc, cosc_b, sinc_b, period)
    assert dist > 0.1
    assert info["period_rel_diff"] == 0.0


def test_gauge_distance_catches_period_difference() -> None:
    """An identical SHAPE at a different period is a physically different
    orbit -- period difference is gauge-invariant by construction (a time
    shift never changes a loop's own period), so this must never collapse."""
    rng = np.random.default_rng(2)
    dc, cosc, sinc, period_a = _random_loop(rng)
    period_b = period_a * 1.3
    dist, info = gauge_distance(dc, cosc, sinc, period_a, dc, cosc, sinc, period_b)
    assert info["shape_dist"] < 1e-6  # identical shape
    assert dist > 0.5  # but combined distance is dominated by the period gap
    assert info["period_rel_diff"] == pytest.approx(0.3 / 1.3, rel=1e-6)


def test_same_family_threshold_behaviour() -> None:
    """Uses ``delta=pi`` for the z-vector-level shift: the z-packing
    (``_unpack``'s convention, see ``_pack_free_anchor_z``) structurally
    forces ``sinc[0, 0] = 0`` (the phase gauge is not a free variable), so
    an ARBITRARY phase shift is not faithfully representable as a z-vector
    at all -- only shifts that happen to preserve that zero are (a
    structural property of this specific parameterization, not a
    limitation of :func:`same_family`/:func:`phase_shift_coeffs`
    themselves, both already validated at the raw-coefficient level by
    ``test_gauge_distance_self_consistency_phase_shift`` above, which
    sweeps arbitrary deltas without going through z-packing).
    ``delta=pi`` is the natural nontrivial such shift: ``cos(k*pi) =
    (-1)^k``, ``sin(k*pi) = 0`` for every integer ``k``, so
    ``new_sinc[0, 0] = sinc[0, 0] * (-1) = 0`` regardless of ``cosc[0, 0]``
    -- a genuine half-period relabeling (flips every ODD harmonic's sign,
    preserves every EVEN one), which the z format CAN represent.
    """
    rng = np.random.default_rng(3)
    n = 8
    dc, cosc, sinc, period = _random_loop(rng, n)
    n_free = _n_free_actual(n, None, None, None)

    z_a = _pack_free_anchor_z(dc, cosc, sinc, period)
    assert z_a.shape == (n_free,)
    cosc_shift, sinc_shift = phase_shift_coeffs(cosc, sinc, np.pi)
    z_b = _pack_free_anchor_z(dc, cosc_shift, sinc_shift, period)

    is_same, info = same_family(z_a, z_b, n)
    assert is_same
    assert info["shape_dist"] < 1e-6

    _dc2, cosc_c, sinc_c, _p2 = _random_loop(np.random.default_rng(4), n)
    z_c = _pack_free_anchor_z(dc, cosc_c, sinc_c, period)
    is_same_c, info_c = same_family(z_a, z_c, n)
    assert not is_same_c
    assert info_c["shape_dist"] > 0.05


def test_magnitude_fingerprint_is_phase_invariant() -> None:
    """The cheap in-loop deflation invariant (harmonic MAGNITUDES) must be
    exactly unchanged under a phase shift -- the analytic property
    :func:`~cyclerfinder.search.deflated_variational_periodic_orbit.
    _magnitude_fingerprint`'s docstring claims, checked directly.

    Uses ``delta=pi`` for the same z-packing-gauge reason documented in
    ``test_same_family_threshold_behaviour`` above (the only nontrivial
    shift a raw z-vector can faithfully represent); a broader arbitrary-
    delta sweep at the raw-coefficient level (bypassing z-packing entirely)
    already lives in ``test_gauge_distance_self_consistency_phase_shift``.
    """
    rng = np.random.default_rng(5)
    n = 6
    dc, cosc, sinc, period = _random_loop(rng, n)
    n_free = _n_free_actual(n, None, None, None)

    z = _pack_free_anchor_z(dc, cosc, sinc, period)
    assert z.shape == (n_free,)
    fp1 = _magnitude_fingerprint(z, n)
    cosc2, sinc2 = phase_shift_coeffs(cosc, sinc, np.pi)
    z2 = _pack_free_anchor_z(dc, cosc2, sinc2, period)
    fp2 = _magnitude_fingerprint(z2, n)
    assert np.linalg.norm(fp1 - fp2) < 1e-8


def test_unpack_roundtrip_matches_gauge_distance_z() -> None:
    """Sanity: ``gauge_distance_z`` on the SAME raw coefficient vector twice
    is exactly zero (basic self-distance property)."""
    rng = np.random.default_rng(6)
    n = 6
    n_free = _n_free_actual(n, None, None, None)
    z = np.zeros(n_free)
    z[0] = 0.8
    z[3:-1] = rng.normal(scale=0.05, size=n_free - 4)
    z[-1] = np.log(2.7)
    dist, info = gauge_distance_z(z, z, n)
    assert dist < 1e-10
    assert info["shape_dist"] < 1e-10


# --------------------------------------------------------------------------
# discover_periodic_orbit_fixed_jacobi -- single-attempt corrector
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_discover_periodic_orbit_fixed_jacobi_reaches_target_jacobi() -> None:
    """A converged fixed-Jacobi solve must actually land AT the requested
    Jacobi constant (not just satisfy the EOM at some nearby energy).
    Marked ``slow`` -- see the multi-family test's docstring for why."""
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    target = 3.000
    n_ok = 0
    for seed in range(6):
        rng = np.random.default_rng(seed)
        res = discover_periodic_orbit_fixed_jacobi(
            sysm,
            target,
            n_harmonics=8,
            center_guess=(0.8369, 0.0, 0.0),
            period_guess=2.2,
            rng=rng,
            tol=1e-5,
            jacobi_tol=1e-5,
            max_nfev=2500,
            n_continuation_steps=5,
        )
        if res.converged:
            n_ok += 1
            assert abs(res.jacobi - target) < 1e-5
            assert res.residual_rms < 1e-5
            # Independent check: propagate through the TRUE EOM.
            assert res.closure_residual < 1e-2
    assert n_ok >= 1, "expected at least one of 6 cold starts to converge"


def test_discover_periodic_orbit_fixed_jacobi_smoke() -> None:
    """Fast, always-on mechanical smoke test: one bounded-budget attempt
    must run end-to-end without crashing and return internally-consistent
    fields (whether or not it happens to converge -- convergence odds are
    covered by the slow test above)."""
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    rng = np.random.default_rng(0)
    res = discover_periodic_orbit_fixed_jacobi(
        sysm,
        3.000,
        n_harmonics=6,
        center_guess=(0.8369, 0.0, 0.0),
        period_guess=2.2,
        rng=rng,
        tol=1e-5,
        jacobi_tol=1e-5,
        max_nfev=800,
        n_continuation_steps=3,
    )
    assert res.n_harmonics == 6
    assert res.period > 0.0
    assert np.isfinite(res.jacobi)
    assert np.isfinite(res.residual_rms)


def test_discover_periodic_orbit_fixed_jacobi_rejects_bad_n_harmonics() -> None:
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError, match="n_harmonics"):
        discover_periodic_orbit_fixed_jacobi(sysm, 3.0, n_harmonics=0)


# --------------------------------------------------------------------------
# enumerate_families_fixed_jacobi -- integration + positive control (bounded)
# --------------------------------------------------------------------------


def test_enumerate_families_fixed_jacobi_smoke() -> None:
    """Fast mechanical wiring smoke test (always-on, small restart budget):
    the enumerator must run end-to-end without crashing, and whatever it
    finds must be internally consistent (Radau-confirmed, no reported
    duplicate). The full multi-family recovery claim is a much more
    expensive stochastic search -- see the `slow`-marked test below and
    `docs/notes/2026-07-19-648-positive-control.md` for the real run.
    """
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    result = enumerate_families_fixed_jacobi(
        sysm,
        3.000,
        n_harmonics=8,
        n_restarts=3,
        center_guess=(0.8369, 0.0, 0.0),
        period_guess=2.2,
        max_nfev=1500,
        n_continuation_steps=4,
        rng=np.random.default_rng(7),
    )
    assert result.n_attempts == 3
    assert result.n_converged_raw >= 0
    for fam in result.families:
        assert fam.radau_ok
        assert fam.radau_closure_residual < 1e-6


@pytest.mark.slow
def test_enumerate_families_finds_at_least_two_distinct_families_bounded() -> None:
    """Bounded regression of the full #648 positive control: at Earth-Moon
    C=3.0 (a window with >=5 documented distinct JPL SSD families -- halo N,
    halo S, planar Lyapunov, vertical, axial -- see
    `docs/notes/2026-07-19-648-positive-control.md`), a deflated enumeration
    run must recover AT LEAST two genuinely distinct families, none of which
    are gauge-duplicate false-positives of each other. Marked ``slow``
    (deselected by default per this project's convention, see
    `pyproject.toml`) -- each cold-start attempt runs a several-step Jacobi-
    constant continuation and can legitimately take tens of seconds; this is
    a capability-demonstration test, not a catalogue V-gauntlet evidence
    test (which this project's own convention forbids marking slow).
    """
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    rng = np.random.default_rng(7)
    result = enumerate_families_fixed_jacobi(
        sysm,
        3.000,
        n_harmonics=8,
        n_restarts=25,
        center_guess=(0.8369, 0.0, 0.0),
        period_guess=2.2,
        rng=rng,
    )
    assert len(result.families) >= 2
    # Cross-check: every pair of accepted families is genuinely distinct
    # under the SAME gauge-invariant test used internally.
    for i in range(len(result.families)):
        for j in range(i + 1, len(result.families)):
            is_same, _info = same_family(
                result.families[i].result.raw_coeffs,
                result.families[j].result.raw_coeffs,
                8,
            )
            assert not is_same
    # Every accepted family passed the mandatory Radau cross-check.
    for fam in result.families:
        assert fam.radau_ok
        assert fam.radau_closure_residual < 1e-6


def test_enumerate_families_rejects_bad_n_harmonics() -> None:
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    with pytest.raises(ValueError, match="n_harmonics"):
        enumerate_families_fixed_jacobi(sysm, 3.0, n_harmonics=0, n_restarts=1)


def test_unpack_still_matches_n_free_actual() -> None:
    """Basic wiring sanity: the free-anchor pattern used throughout this
    module (anchor_x1=anchor_y1=anchor_z1=None) round-trips through
    ``_unpack``/``_n_free_actual`` without shape mismatches."""
    n = 6
    n_free = _n_free_actual(n, None, None, None)
    z = np.zeros(n_free)
    z[0:3] = [0.8, 0.0, 0.0]
    z[-1] = np.log(2.7)
    dc, cosc, sinc, period = _unpack(z, n, None, None, None)
    assert dc.shape == (3,)
    assert cosc.shape == (3, n)
    assert sinc.shape == (3, n)
    assert period == pytest.approx(2.7)
    assert sinc[0, 0] == 0.0  # phase gauge always pinned
