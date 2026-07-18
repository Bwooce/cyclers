"""Deflated enumeration of distinct CR3BP periodic-orbit families at fixed Jacobi (#648).

Combines two pieces of machinery that had never been combined before this
task: `deflated_newton.py` (#524, previously only ever aimed at
basin-restricted SHOOTING residuals) and the #606 seedless spectral
(harmonic-balance) periodic-orbit corrector
(`variational_periodic_orbit.py`), whose own module docstring documents a
family-SELECTION BIAS -- a free amplitude tends to slide onto the easier,
already-known vertical-Lyapunov family instead of a genuine halo family even
when both exist at the same energy level. This module turns that corrector
into a systematic DISTINCT-FAMILY ENUMERATOR at a fixed Jacobi constant:
repeatedly deflated-solve :func:`~cyclerfinder.search.
variational_periodic_orbit.discover_periodic_orbit_fixed_jacobi`'s residual
from independent random cold starts, pushing each new attempt away from every
family already found, until the restart budget is exhausted.

Gauge invariance (the load-bearing correctness issue)
-------------------------------------------------------
A truncated real Fourier series ``f(theta) = dc + sum_k [c_k cos(k*theta) +
s_k sin(k*theta)]`` parameterizes a closed loop by an ARBITRARY choice of
time origin (``theta=0``). The SAME physical periodic orbit, refit from two
different cold starts, can converge to two Fourier-coefficient vectors that
differ by nothing but a phase shift ``theta -> theta + delta`` -- naive
Euclidean deflation (`deflated_newton.deflation_factor`'s literal distance)
would then either (a) fail to repel a re-discovery of an already-found
family expressed at a different phase (missing the intended repulsion,
re-reporting the SAME family many times as if each phase were a distinct
root), or (b) if by luck it does land near-coincidentally on the SAME phase,
work correctly only by accident. Two invariants are used here, doing
DIFFERENT jobs:

1. **In-loop repulsion (cheap, smooth, exact)**: shifting the time origin by
   ``delta`` rotates each harmonic's ``(c_k, s_k)`` pair by angle
   ``k*delta`` (a standard Fourier shift-theorem fact, derived in
   :func:`phase_shift_coeffs`'s docstring) -- so the per-harmonic COMPLEX
   MAGNITUDE ``sqrt(c_k^2 + s_k^2)`` is EXACTLY invariant under any phase
   shift, with no alignment search needed. :func:`_magnitude_fingerprint`
   builds this invariant descriptor (offsets, period, and per-harmonic
   magnitudes); `deflated_newton.deflation_factor` -- REUSED verbatim, not
   reimplemented -- is then called on these fingerprints instead of the raw
   coefficient vectors, giving a genuinely gauge-invariant Farrell deflation
   multiplier with no FFT/argmax inside the hot loop. This is a NECESSARY
   (fingerprint-equal for a true phase-shift) but not manifestly SUFFICIENT
   invariant (two different-shaped loops could in principle share magnitudes
   while differing in the relative phase BETWEEN coordinates); it is used
   only to steer the optimizer's repulsion, never for the final accept/
   reject decision.
2. **Final classification (expensive, explicit alignment, task-specified)**:
   :func:`gauge_distance` sample-reconstructs both loops on a dense
   ``theta`` grid and finds the SINGLE best-aligning circular shift via FFT
   cross-correlation (``ifft(fft(a) * conj(fft(b)))``, standard phase-
   correlation), applied jointly across all three coordinates (so a
   relative-phase difference between x/y/z that the magnitude fingerprint
   alone cannot see is still caught), then reports the RMS distance between
   the aligned sampled trajectories plus the (gauge-invariant-by-
   construction) period. :func:`same_family` thresholds this into the
   accept/reject decision used both for de-duplicating a genuinely new find
   against the accepted-family list and for the self-consistency test
   (`tests/search/test_deflated_variational_periodic_orbit.py`) that
   directly verifies: phase-shifting a converged orbit's OWN coefficients
   and comparing against the original yields near-zero distance.

Ghost-minima discipline (#620)
-------------------------------
The seedless corrector has NO forward integration inside its own search (the
entire point of #606's method) -- a harmonic-balance residual near zero does
NOT by itself certify a genuine closed trajectory; #620 hit exactly this
failure mode in a related collocation-arc corrector (machine-zero nodal
residual, endpoints matched, yet an independent Radau re-propagation showed
an O(1) loop defect). Every candidate that clears the Fourier residual AND
Jacobi-matching gates here is therefore re-propagated with
``topology_audit.check_periodic_orbit_closure`` -- Radau, a DIFFERENT
integrator from anything used in the search itself -- before being accepted
as a genuine family; candidates failing this are counted and reported
separately (``n_rejected_ghost``), never silently dropped.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares, minimize_scalar

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.deflated_newton import deflation_factor
from cyclerfinder.search.topology_audit import (
    INDEPENDENT_CLOSURE_FLOOR_NONDIM,
    check_periodic_orbit_closure,
)
from cyclerfinder.search.variational_periodic_orbit import (
    VariationalOrbitResult,
    _eval_series,
    _harmonic_balance_residual,
    _n_free_actual,
    _reconstruct_state0,
    _residual_fixed_jacobi,
    _unpack,
)

__all__ = [
    "EnumeratedFamily",
    "EnumerationResult",
    "enumerate_families_fixed_jacobi",
    "gauge_distance",
    "phase_shift_coeffs",
    "same_family",
]


# --------------------------------------------------------------------------
# Gauge-invariant distance (final classification -- see module docstring #2)
# --------------------------------------------------------------------------


def phase_shift_coeffs(
    cosc: NDArray[np.float64], sinc: NDArray[np.float64], delta_theta: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``(cos, sin)`` harmonic coefficients of the SAME loop with its
    time origin shifted by ``delta_theta`` radians of phase angle.

    Derivation: ``f(theta) = dc + sum_k [c_k cos(k theta) + s_k sin(k theta)]``.
    Substituting ``theta -> theta + delta`` and expanding
    ``cos(k theta + k delta)`` / ``sin(k theta + k delta)`` via the angle-sum
    identities and collecting the ``cos(k theta)`` / ``sin(k theta)``
    coefficients gives a per-harmonic 2D ROTATION by angle ``k * delta``:

        new_c_k =  c_k cos(k delta) + s_k sin(k delta)
        new_s_k = -c_k sin(k delta) + s_k cos(k delta)

    ``cosc``/``sinc`` may be a single coordinate's ``(n_harmonics,)`` vector
    or the full ``(3, n_harmonics)`` array (broadcasts over the leading axis
    identically for every coordinate, since a single shared time-origin
    shift applies to all three coordinates of one loop together).
    """
    n = cosc.shape[-1]
    k = np.arange(1, n + 1, dtype=np.float64)
    c = np.cos(k * delta_theta)
    s = np.sin(k * delta_theta)
    new_cos = cosc * c + sinc * s
    new_sin = -cosc * s + sinc * c
    return new_cos, new_sin


def _sample_loop(
    dc: NDArray[np.float64],
    cosc: NDArray[np.float64],
    sinc: NDArray[np.float64],
    n_samples: int,
) -> NDArray[np.float64]:
    """Dense ``(n_samples, 3)`` (x, y, z) sample of one Fourier loop over ``theta in [0, 2pi)``."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_samples, endpoint=False)
    cols = []
    for c in range(3):
        f, _fp, _fpp = _eval_series(theta, float(dc[c]), cosc[c], sinc[c])
        cols.append(f)
    return np.stack(cols, axis=1)


def gauge_distance(
    dc_a: NDArray[np.float64],
    cosc_a: NDArray[np.float64],
    sinc_a: NDArray[np.float64],
    period_a: float,
    dc_b: NDArray[np.float64],
    cosc_b: NDArray[np.float64],
    sinc_b: NDArray[np.float64],
    period_b: float,
    *,
    n_samples: int = 256,
) -> tuple[float, dict[str, float]]:
    """Phase-minimized (gauge-invariant) distance between two Fourier loops.

    Samples both loops on a dense ``theta`` grid, finds the SINGLE best
    circular shift aligning loop B onto loop A via FFT cross-correlation
    (``ifft(fft(a) * conj(fft(b)))``, summed across the three coordinates so
    the shift is chosen jointly, not per-coordinate -- a relative-phase
    difference between coordinates is exactly what this catches and the
    magnitude-fingerprint invariant used for in-loop deflation cannot), then
    reports:

    * ``shape_dist``: RMS distance between loop A and the best-aligned loop B
      (near zero for a genuine phase-shifted duplicate of the SAME orbit,
      per the self-consistency test).
    * ``period_abs_diff`` / ``period_rel_diff``: the period difference (a
      time-shift never changes an orbit's OWN period, so any nonzero period
      difference is real evidence of a physically different orbit, not a
      gauge artifact).

    Returns ``(combined, info)`` where ``combined = sqrt(shape_dist**2 +
    period_abs_diff**2)`` (both terms are O(1) in these nondimensional CR3BP
    units, so no further rescaling is applied) -- a single smooth-ish scalar
    convenient for reporting; :func:`same_family` uses ``info`` directly for
    the sharper two-part accept/reject decision.

    Two-stage alignment: a coarse DISCRETE circular shift via FFT cross-
    correlation (grid resolution ``2*pi/n_samples``) locates the right
    vicinity cheaply, then a bounded 1D refinement (``scipy.optimize.
    minimize_scalar``) applies :func:`phase_shift_coeffs` -- an EXACT,
    continuous, non-interpolated rotation of the harmonic coefficients, not
    a resampling trick -- to polish the shift to near machine precision
    within one grid cell of the coarse estimate. Without this refinement the
    reported distance for a genuine phase-shifted duplicate floors at the
    sampling grid's OWN resolution error (~1e-2 at ``n_samples=256`` for a
    typical few-percent-amplitude orbit) rather than at machine precision --
    confirmed by :func:`test_gauge_distance_self_consistency_phase_shift`.
    """
    sa = _sample_loop(dc_a, cosc_a, sinc_a, n_samples)
    sb = _sample_loop(dc_b, cosc_b, sinc_b, n_samples)
    corr = np.zeros(n_samples)
    for c in range(3):
        fa = np.fft.rfft(sa[:, c])
        fb = np.fft.rfft(sb[:, c])
        corr += np.fft.irfft(fa * np.conj(fb), n=n_samples)
    shift_idx = int(np.argmax(corr))
    grid_step = 2.0 * np.pi / n_samples
    # ``np.roll(sb, +shift_idx)`` (a SAMPLE-INDEX roll) is the correlation-
    # peak convention, but :func:`phase_shift_coeffs` shifts the FUNCTION
    # ARGUMENT (f(theta) -> f(theta + delta)), which is the OPPOSITE sense
    # (rolling samples forward by m corresponds to evaluating the original
    # function at an EARLIER argument at each fixed sample point) -- hence
    # the sign flip. Verified directly: for a known delta, the analytic
    # refinement below only converges from ``-shift_idx * grid_step``, never
    # ``+shift_idx * grid_step`` (see test_deflated_variational_periodic_orbit.py).
    delta0 = -shift_idx * grid_step

    def _obj(delta: float) -> float:
        cos_shift, sin_shift = phase_shift_coeffs(cosc_b, sinc_b, float(delta))
        sb_shift = _sample_loop(dc_b, cos_shift, sin_shift, n_samples)
        return float(np.mean(np.sum((sa - sb_shift) ** 2, axis=1)))

    refined = minimize_scalar(
        _obj,
        bounds=(delta0 - grid_step, delta0 + grid_step),
        method="bounded",
        options={"xatol": 1e-13},
    )
    shape_dist = float(np.sqrt(refined.fun))
    period_abs_diff = abs(float(period_a) - float(period_b))
    denom = max(float(period_a), float(period_b), 1e-9)
    period_rel_diff = period_abs_diff / denom
    combined = math.sqrt(shape_dist**2 + period_abs_diff**2)
    return combined, {
        "shape_dist": shape_dist,
        "period_abs_diff": period_abs_diff,
        "period_rel_diff": period_rel_diff,
        "shift_samples": float(shift_idx),
        "shift_delta_theta": float(refined.x),
    }


def gauge_distance_z(
    z_a: NDArray[np.float64], z_b: NDArray[np.float64], n_harmonics: int, *, n_samples: int = 256
) -> tuple[float, dict[str, float]]:
    """:func:`gauge_distance` on two raw free-variable vectors (anchors all free)."""
    dc_a, cosc_a, sinc_a, period_a = _unpack(z_a, n_harmonics, None, None, None)
    dc_b, cosc_b, sinc_b, period_b = _unpack(z_b, n_harmonics, None, None, None)
    return gauge_distance(
        dc_a, cosc_a, sinc_a, period_a, dc_b, cosc_b, sinc_b, period_b, n_samples=n_samples
    )


def same_family(
    z_a: NDArray[np.float64],
    z_b: NDArray[np.float64],
    n_harmonics: int,
    *,
    shape_tol: float = 1e-3,
    period_rel_tol: float = 1e-3,
    n_samples: int = 256,
) -> tuple[bool, dict[str, float]]:
    """Strict gauge-invariant same-family test used for final de-duplication.

    ``True`` iff BOTH the phase-aligned shape distance is below ``shape_tol``
    AND the relative period difference is below ``period_rel_tol`` -- period
    alone is a cheap, exact discriminator (unaffected by phase gauge), shape
    catches loops with matching period but a genuinely different geometry.
    """
    _combined, info = gauge_distance_z(z_a, z_b, n_harmonics, n_samples=n_samples)
    is_same = info["shape_dist"] < shape_tol and info["period_rel_diff"] < period_rel_tol
    return is_same, info


# --------------------------------------------------------------------------
# In-loop deflation (cheap, smooth, exact invariant -- see module docstring #1)
# --------------------------------------------------------------------------


def _magnitude_fingerprint(z: NDArray[np.float64], n_harmonics: int) -> NDArray[np.float64]:
    """Exactly phase-shift-INVARIANT descriptor of a loop: offsets, period, and
    per-harmonic-per-coordinate complex magnitudes ``sqrt(c_k^2 + s_k^2)``.

    A time-origin shift rotates each ``(c_k, s_k)`` pair (see
    :func:`phase_shift_coeffs`), which leaves its magnitude UNCHANGED --
    so this fingerprint is bit-for-bit identical (up to the loop's own
    self-consistency) for any two phase-shifted representations of the SAME
    orbit, with no search/alignment needed. NOT a full gauge-invariant
    distance on its own (see module docstring); used only to weight the
    in-loop Farrell deflation multiplier via `deflated_newton.
    deflation_factor`, never for the final accept/reject decision (that is
    :func:`same_family`'s job).
    """
    dc, cosc, sinc, period = _unpack(z, n_harmonics, None, None, None)
    mag = np.sqrt(cosc**2 + sinc**2)
    return np.concatenate([dc, np.array([period]), mag.ravel()])


def _deflated_residual_fixed_jacobi(
    z: NDArray[np.float64],
    theta: NDArray[np.float64],
    mu: float,
    n_harmonics: int,
    jacobi_target: float,
    jacobi_weight: float,
    known_root_z: list[NDArray[np.float64]],
    p: float,
    shift: float,
) -> NDArray[np.float64]:
    """The #606 fixed-Jacobi residual, Farrell-deflated (#524) against
    ``known_root_z`` via the gauge-invariant magnitude fingerprint."""
    base = _residual_fixed_jacobi(z, theta, mu, n_harmonics, jacobi_target, jacobi_weight)
    if not known_root_z:
        return base
    fp_z = _magnitude_fingerprint(z, n_harmonics)
    fp_known = [_magnitude_fingerprint(u, n_harmonics) for u in known_root_z]
    m, _grad = deflation_factor(fp_z, fp_known, p=p, shift=shift)
    if not math.isfinite(m):
        # Exactly on a known root's fingerprint (measure-zero in float
        # arithmetic) -- return a large-but-finite residual so
        # scipy.optimize.least_squares (which requires a finite ndarray, not
        # None) is repelled without crashing.
        return np.full_like(base, 1.0e8)
    return m * base


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class EnumeratedFamily:
    """One accepted, Radau-cross-checked, gauge-deduplicated periodic-orbit family."""

    result: VariationalOrbitResult
    radau_closure_residual: float
    radau_ok: bool
    min_gauge_distance_to_prior: float
    """The minimum :func:`gauge_distance` to every PREVIOUSLY accepted family
    (``inf`` for the first family found) -- an audit trail confirming this
    was a genuinely distinct find, not a borderline dedup call."""


@dataclass(frozen=True)
class EnumerationResult:
    """Outcome of one :func:`enumerate_families_fixed_jacobi` run."""

    families: list[EnumeratedFamily]
    n_attempts: int
    n_converged_raw: int
    n_rejected_ghost: int
    n_rejected_duplicate: int
    jacobi_target: float
    system_label: str = field(default="")


def enumerate_families_fixed_jacobi(
    system: cr3bp.CR3BPSystem,
    jacobi_target: float,
    *,
    n_harmonics: int = 8,
    n_collocation: int | None = None,
    n_restarts: int = 60,
    center_guess: tuple[float, float, float] = (0.8, 0.0, 0.0),
    period_guess: float = 3.0,
    period_guess_range: tuple[float, float] = (0.5, 1.8),
    coefficient_noise: float = 0.05,
    jacobi_weight: float = 50.0,
    n_continuation_steps: int = 8,
    tol: float = 1e-8,
    jacobi_tol: float = 1e-6,
    max_nfev: int = 15000,
    p: float = 2.0,
    shift: float = 1.0,
    radau_closure_floor: float = INDEPENDENT_CLOSURE_FLOOR_NONDIM,
    dedup_shape_tol: float = 1e-3,
    dedup_period_rel_tol: float = 1e-3,
    rng: np.random.Generator | None = None,
) -> EnumerationResult:
    """Deflated enumeration of distinct periodic-orbit families at fixed Jacobi (#648).

    From ``n_restarts`` independent random cold starts, deflated-solves
    (see module docstring) :func:`~cyclerfinder.search.
    variational_periodic_orbit.discover_periodic_orbit_fixed_jacobi`'s
    residual (Jacobi-fixing row, all amplitude anchors free) via
    ``scipy.optimize.least_squares`` (Levenberg-Marquardt, the SAME robust
    solver #606/#611/#612 already validated for this residual -- deflation
    is baked into the residual function it optimizes, not routed through a
    hand-rolled undamped Newton loop, which would be materially less
    reliable on this stiff a cold-start problem). Every candidate that
    clears the (undeflated) residual/Jacobi gates is Radau-cross-checked
    (#620 ghost-minima discipline) and gauge-deduplicated
    (:func:`same_family`) against every family already accepted before being
    reported as genuinely new.

    A seed whose deflated solve diverges, stalls, times out, or reconverges
    to an already-known family (rejected by the dedup check) contributes
    nothing to ``families`` -- exactly the same "expected behaviour, not an
    error" semantics as `deflated_newton.enumerate_roots`.
    """
    if n_harmonics < 1:
        raise ValueError(
            f"enumerate_families_fixed_jacobi: n_harmonics must be >= 1, got {n_harmonics}"
        )
    mu = system.mu
    n_coll = n_collocation if n_collocation is not None else 6 * n_harmonics
    theta = np.linspace(0.0, 2.0 * np.pi, n_coll, endpoint=False)
    n_free = _n_free_actual(n_harmonics, None, None, None)
    gen = rng if rng is not None else np.random.default_rng()

    known_root_z: list[NDArray[np.float64]] = []
    families: list[EnumeratedFamily] = []
    n_converged_raw = 0
    n_rejected_ghost = 0
    n_rejected_duplicate = 0

    for _attempt in range(n_restarts):
        z0 = np.zeros(n_free)
        z0[0] = center_guess[0] + gen.normal(scale=0.05)
        z0[1] = center_guess[1] + gen.normal(scale=0.02)
        z0[2] = center_guess[2] + gen.normal(scale=0.02)
        z0[3:-1] = gen.normal(scale=coefficient_noise, size=n_free - 4)
        z0[-1] = np.log(period_guess * float(gen.uniform(*period_guess_range)))

        # Natural-parameter continuation in the target Jacobi constant (see
        # variational_periodic_orbit._fixed_jacobi_continuation_solve's
        # docstring for why a single cold-start solve is unreliable) --
        # ``known_root_z`` (hence the deflation multiplier) is held FIXED
        # across the continuation steps of one attempt; only the target C
        # schedule and the warm-started coefficients change.
        state0_0, _period0 = _reconstruct_state0(z0, n_harmonics, None, None, None)
        c0 = cr3bp.jacobi_constant(state0_0, mu)
        schedule = np.linspace(c0, jacobi_target, max(2, n_continuation_steps))[1:]
        z = z0.copy()
        sol = None
        for c_step in schedule:
            sol = least_squares(
                _deflated_residual_fixed_jacobi,
                z,
                args=(
                    theta,
                    mu,
                    n_harmonics,
                    float(c_step),
                    jacobi_weight,
                    known_root_z,
                    p,
                    shift,
                ),
                method="lm",
                xtol=1e-15,
                ftol=1e-15,
                gtol=1e-15,
                max_nfev=max_nfev,
            )
            z = sol.x
        assert sol is not None

        state0, period = _reconstruct_state0(z, n_harmonics, None, None, None)
        if not (0.05 < period < 50.0):
            continue
        base_res = _harmonic_balance_residual(z, theta, mu, n_harmonics, None, None, None)
        residual_rms = float(np.sqrt(np.mean(base_res**2)))
        jacobi_actual = cr3bp.jacobi_constant(state0, mu)
        jacobi_err = abs(jacobi_actual - jacobi_target)
        if not (residual_rms < tol and jacobi_err < jacobi_tol):
            continue
        n_converged_raw += 1

        radau_ok, radau_info = check_periodic_orbit_closure(
            state0, period, system, floor=radau_closure_floor
        )
        if not radau_ok:
            n_rejected_ghost += 1
            continue

        min_dist = math.inf
        is_dup = False
        for u in known_root_z:
            is_same, info = same_family(
                z,
                u,
                n_harmonics,
                shape_tol=dedup_shape_tol,
                period_rel_tol=dedup_period_rel_tol,
            )
            combined = math.sqrt(info["shape_dist"] ** 2 + info["period_abs_diff"] ** 2)
            min_dist = min(min_dist, combined)
            if is_same:
                is_dup = True
        if is_dup:
            n_rejected_duplicate += 1
            continue

        vresult = VariationalOrbitResult(
            state0=state0,
            period=period,
            jacobi=jacobi_actual,
            converged=True,
            residual_rms=residual_rms,
            closure_residual=radau_info["closure_residual"],
            n_harmonics=n_harmonics,
            n_collocation=n_coll,
            n_restarts_tried=1,
            cost=float(sol.cost),
            raw_coeffs=z.copy(),
        )
        families.append(
            EnumeratedFamily(
                result=vresult,
                radau_closure_residual=radau_info["closure_residual"],
                radau_ok=True,
                min_gauge_distance_to_prior=min_dist,
            )
        )
        known_root_z.append(z.copy())

    return EnumerationResult(
        families=families,
        n_attempts=n_restarts,
        n_converged_raw=n_converged_raw,
        n_rejected_ghost=n_rejected_ghost,
        n_rejected_duplicate=n_rejected_duplicate,
        jacobi_target=jacobi_target,
        system_label=f"{system.primary}-{system.secondary}",
    )
