"""Impulsive primer-vector optimality diagnostic (READ-ONLY).

Classical first-order necessary conditions for the optimality of an
*impulsive* N-burn trajectory in a two-body (Kepler) field, after Lawden's
primer vector theory and the Lion & Handelsman fixed-time diagnostic.

Theory (citations)
------------------
* Lawden, D. F. (1963). *Optimal Trajectories for Space Navigation.*
  Butterworths. — defines the primer vector ``p(t)`` as the velocity
  adjoint (costate) of the impulsive minimum-ΔV problem.
* Lion, P. M., & Handelsman, M. (1968). *Primer Vector on Fixed-Time
  Impulsive Trajectories.* AIAA Journal 6(1), 127-132. DOI 10.2514/3.4452.
  — the *diagnostic* used here: propagate ``p(t)`` for a GIVEN (possibly
  non-optimal) fixed-time impulse schedule; if ``|p(t)| > 1`` on any coast
  the schedule is non-optimal and an added/relocated midcourse impulse
  lowers the total ΔV.
* Prussing, J. E. (2010). *Primer Vector Theory and Applications*, Ch. 2 in
  Conway (ed.), *Spacecraft Trajectory Optimization*, Cambridge UP. —
  standard textbook treatment, including the coplanar-circular Hohmann vs
  bi-elliptic threshold used as the sourced golden gate (the Hohmann
  transfer satisfies the primer necessary conditions iff the radius ratio
  is below ~11.94; above ~15.58 bi-elliptic is strictly better).

The primer vector
-----------------
On a ballistic (coast) arc the primer satisfies the same linear variational
equation as a position perturbation in the two-body field::

    p̈ = G(r(t)) · p ,   G(r) = (μ / r³) (3 r̂ r̂ᵀ - I)

i.e. ``G`` is the gravity-gradient (tidal) matrix of the Kepler field along
the reference arc ``r(t)``. Because ``(p, ṗ)`` obeys the *same* linear ODE as
the state perturbation, the primer over a coast propagates with the arc's
state-transition matrix (STM)::

    [p(t); ṗ(t)] = Φ(t, t0) · [p(t0); ṗ(t0)]

Boundary conditions (Lawden): at each impulse ``i`` the primer equals the
unit ΔV direction, ``p(t_i) = Δv_i / |Δv_i|`` (so ``|p(t_i)| = 1`` by
construction). For an interior coast bounded by two impulses we therefore
have a two-point BVP: ``p`` is pinned (as a 3-vector) at both ends, and the
unknown ``ṗ(t0)`` is recovered from the STM by inverting the upper-right
3x3 block ``Φ_rv``.

Caveat — degenerate (zero-magnitude) bounding impulse: the BC
``p(t_i) = Δv_i / |Δv_i|`` (Guzman 2002, Eq. 33) is **ill-posed as
``|Δv_i| → 0``** — the direction is undefined in the limit. Our Aldrin coast-0
Earth departure is *exactly* this degenerate case (``|Δv_0| = 0`` by
construction), so the endpoint primer direction supplied there is a fallback
unit vector, not a physical ``Δv/|Δv|``. The resulting ``|p|`` bulge on such a
coast is partly an artifact of the imposed (non-physical) endpoint direction
and should be read with that in mind; the survey assumes non-degenerate
impulses and does not cover this boundary condition.

Diagnostic verdict
------------------
Necessary (NOT sufficient) conditions for an optimal impulse schedule:

1. ``|p(t)| ≤ 1`` throughout every coast.
2. ``|p(t_i)| = 1`` at each impulse (satisfied here by construction of the
   BCs).

If ``max|p| > 1`` on a coast the schedule fails condition 1: the published
interpretation (Lion & Handelsman) is that an added or relocated impulse
near the peak reduces total ΔV. This module reports the per-coast
``max|p|`` and its time, and a verdict enum. It does **not** re-optimise.

Caveat (Guzman 2002, NOT yet acquired): linearised primer theory has
singularities and degrades on long multi-revolution arcs (Guzman, Mailhe,
Schiff, Hughes & Folta 2002, IAC-02-A.6.09 / NTRS 20030032208). Long
multi-rev cycler legs are exactly such a case; treat results on them as
DIAGNOSTIC / PROVISIONAL until that survey is in hand and the STM
propagation is validated against its cases.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.constants import MU_SUN_KM3_S2

Vec3 = NDArray[np.float64]


class PrimerVerdict(enum.Enum):
    """Per-coast (and overall) primer optimality verdict.

    The check is a *necessary*-conditions test only (Lawden / Lion &
    Handelsman); it can refute optimality but never prove it.
    """

    OPTIMAL_NECESSARY_CONDITIONS_MET = "optimal_necessary_conditions_met"
    """``|p(t)| ≤ 1`` throughout the coast (within tolerance). The schedule
    satisfies the first-order necessary conditions on this coast; this is NOT
    a proof of optimality (sufficiency is not checked)."""

    IMPROVABLE_ADD_IMPULSE = "improvable_add_impulse"
    """``max|p| > 1`` on the coast: the schedule is non-optimal here and the
    published interpretation is that an added / relocated impulse near the
    peak reduces total ΔV."""

    INDETERMINATE_ILL_CONDITIONED = "indeterminate_ill_conditioned"
    """The two-point primer BVP is singular / ill-conditioned on this coast —
    the regime of a near-integer-revolution ("phasing-orbit") coast, where the
    coast STM tends to the identity and its position-velocity block ``Φ_rv``
    becomes rank-deficient (Şaloğlu, Taheri & Landau 2023, Sec. III.F; Glandorf
    1969). The endpoint primer directions cannot be independently enforced on
    such an arc, so the two-point inversion (Eq. 28) returns garbage. We instead
    propagate the primer by continuity (truncated-SVD minimum-norm ``ṗ(0)``,
    dropping the unobservable null direction of ``Φ_rv``), which keeps ``|p(t)|``
    finite, and report the coast as INDETERMINATE rather than emitting a
    spurious IMPROVABLE verdict. See ``primer_on_coast`` for the citation and
    :data:`_ILL_CONDITIONED_CAVEAT`."""


_NECESSARY_CONDITIONS_CAVEAT: str = (
    "NECESSARY-conditions check only (Lawden / Lion & Handelsman 1968): "
    "|p| <= 1 can refute optimality but does not prove it (sufficiency not "
    "checked). Linearised primer theory degrades on long multi-rev arcs "
    "(Guzman 2002, not yet acquired) -- treat long-arc results as provisional."
)

_ILL_CONDITIONED_CAVEAT: str = (
    "At least one coast is a near-integer-revolution (phasing-orbit) arc whose "
    "two-point primer BVP is singular (Phi_rv rank-deficient; Saloglu, Taheri & "
    "Landau 2023 Sec. III.F, Glandorf 1969). The endpoint primer directions "
    "cannot be enforced there; the primer was propagated by continuity "
    "(truncated-SVD minimum-norm pdot0) and the coast verdict is INDETERMINATE, "
    "NOT a true IMPROVABLE/OPTIMAL result. On such multi-rev coasts an interior "
    "|p|->1 touch with pdot~0 is an iso-DeltaV impulse degeneracy consistent "
    "with optimality, not an add-an-impulse signal."
)

# Reciprocal-condition-number floor below which Phi_rv is treated as singular, so
# the continuity (truncated-SVD) primer-rate solve is used instead of a direct
# inversion. Picked well above float64 round-off (~1e-16) yet far below the
# conditioning of any well-separated Kepler coast (rcond ~ 1e-2 in the empirical
# sweep). NOTE: Phi_rv is singular at EVERY half-integer revolution (transfer
# angle = k*180 deg), which includes the legitimate 180-deg Hohmann transfer --
# so a singular Phi_rv alone is NOT a defect signal. The defect (the false
# IMPROVABLE the direct inversion produces on multi-rev coasts) is detected by
# the BVP RESIDUAL below: a singular arc is only INDETERMINATE when the endpoint
# primer direction also lies OUTSIDE the range of Phi_rv (so the two-point BC is
# genuinely unenforceable). Hohmann's BC happens to lie in-range -> residual ~0
# -> trustworthy; an integer-rev phasing coast with mismatched directions has a
# residual ~O(1) -> INDETERMINATE.
_PHI_RV_RCOND_FLOOR: float = 1.0e-8

# Tolerance on the two-point BVP residual ||Phi_rr p0 + Phi_rv pdot0 - p1|| used,
# ONLY on a rank-deficient (singular) Phi_rv, to decide whether the prescribed
# endpoint primer direction is actually reachable. Well above the ~1e-14 residual
# of an in-range (e.g. Hohmann) BC and far below the ~1e-1..1e0 residual of a
# genuinely unenforceable integer-rev BC.
_PRIMER_BVP_RESIDUAL_TOL: float = 1.0e-6


def gravity_gradient(r: Vec3, mu: float = MU_SUN_KM3_S2) -> NDArray[np.float64]:
    """Gravity-gradient (tidal) matrix ``G(r) = (μ/r³)(3 r̂ r̂ᵀ - I)``.

    This is the Jacobian ``∂(-μ r / r³)/∂r`` of the two-body acceleration;
    the primer obeys ``p̈ = G p`` along a coast. ``G`` is symmetric and
    traceless by construction (``tr(3 r̂ r̂ᵀ - I) = 3 - 3 = 0``).

    Parameters
    ----------
    r:
        Heliocentric position, km, ``(3,)``.
    mu:
        Central-body gravitational parameter, km³/s². Defaults to the Sun.

    Returns
    -------
    (3, 3) float64 symmetric, traceless matrix (units 1/s²).
    """
    r_arr = np.asarray(r, dtype=np.float64)
    r_n = float(np.linalg.norm(r_arr))
    r_hat = r_arr / r_n
    return (mu / r_n**3) * (3.0 * np.outer(r_hat, r_hat) - np.eye(3))


def _coast_stm(
    r0: Vec3,
    v0: Vec3,
    duration_s: float,
    mu: float,
    *,
    n_samples: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Integrate the 6x6 STM of the variational equations along a Kepler arc.

    Integrates the reference state ``(r, v)`` and the state-transition matrix
    ``Φ(t, t0)`` of the linear variational system

        d/dt [δr; δv] = [[0, I], [G(r(t)), 0]] [δr; δv]

    over ``[0, duration_s]``. The position-perturbation block of ``Φ`` is the
    same operator that maps an initial ``(p, ṗ)`` to a later ``(p, ṗ)`` (the
    primer obeys the identical linear ODE), so this STM is reused to solve the
    primer two-point BVP and to evaluate ``|p(t)|`` on a sample grid.

    Returns
    -------
    (times, ref_states, stms):
        ``times`` ``(n_samples,)`` from 0 to ``duration_s``; ``ref_states``
        ``(n_samples, 6)`` the propagated reference ``(r, v)``; ``stms``
        ``(n_samples, 6, 6)`` the STM ``Φ(t_k, 0)``.
    """
    y0 = np.empty(6 + 36, dtype=np.float64)
    y0[:3] = np.asarray(r0, dtype=np.float64)
    y0[3:6] = np.asarray(v0, dtype=np.float64)
    y0[6:] = np.eye(6, dtype=np.float64).ravel()

    def rhs(_t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        r = y[:3]
        v = y[3:6]
        phi = y[6:].reshape(6, 6)
        g = gravity_gradient(r, mu)
        a = np.zeros((6, 6), dtype=np.float64)
        a[:3, 3:] = np.eye(3)
        a[3:, :3] = g
        dy = np.empty_like(y)
        dy[:3] = v
        # Two-body acceleration on the reference arc.
        r_n = float(np.linalg.norm(r))
        dy[3:6] = -mu * r / r_n**3
        dy[6:] = (a @ phi).ravel()
        return dy

    t_eval = np.linspace(0.0, duration_s, n_samples)
    sol = solve_ivp(
        rhs,
        (0.0, duration_s),
        y0,
        t_eval=t_eval,
        method="DOP853",
        rtol=1.0e-10,
        atol=1.0e-10,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"primer STM integration failed: {sol.message}")
    ys = sol.y.T  # (n_samples, 42)
    ref_states = ys[:, :6].copy()
    stms = ys[:, 6:].reshape(-1, 6, 6).copy()
    return sol.t.copy(), ref_states, stms


def _solve_primer_rate(
    phi_rr: NDArray[np.float64],
    phi_rv: NDArray[np.float64],
    p0: Vec3,
    p1: Vec3,
) -> tuple[NDArray[np.float64], bool, float]:
    """Recover ``ṗ(0)`` from the coast STM, robust to near-integer-rev coasts.

    The primer two-point BVP on a coast bounded by two impulses fixes ``p(0)``
    and ``p(T)`` (the unit ΔV directions); the unknown initial rate is recovered
    from the STM end-block decomposition ``p(T) = Φ_rr p(0) + Φ_rv ṗ(0)``, i.e.
    ``ṗ(0) = Φ_rv⁻¹ (p(T) - Φ_rr p(0))`` (Glandorf 1969, Eq. 28 of Şaloğlu,
    Taheri & Landau 2023).

    Over a (near-)integer-revolution **phasing orbit** the coast STM tends to
    the identity and ``Φ_rv`` becomes rank-deficient (Saloglu 2023 Sec. III.F):
    a direct inversion then amplifies round-off by ``~1/sigma_min ≈ 10¹¹`` and the
    diagnostic silently returns a primer of magnitude ``~10⁹`` — a false
    IMPROVABLE on exactly the multi-rev cycler-leg regime we care about. The
    *continuity* fix the paper prescribes is to propagate the primer through the
    phasing orbit rather than re-solve the singular two-point problem; at the
    single-coast level that is the **truncated-SVD minimum-norm** solution,
    which discards the unobservable null direction of ``Φ_rv`` (along which the
    endpoint BC carries no information) and so keeps ``ṗ(0)`` — and hence
    ``|p(t)|`` — finite and continuous. On a full-rank coast the truncation
    drops nothing and the result is bit-for-bit the direct solve.

    ``Φ_rv`` is singular at *every* half-integer revolution (transfer angle a
    multiple of 180°), which includes the legitimate 180° Hohmann transfer — so
    a singular ``Φ_rv`` alone is not a defect. The coast is reported
    ``ill_conditioned`` (⇒ INDETERMINATE) only when ``Φ_rv`` is rank-deficient
    **and** the prescribed endpoint direction lies outside its range, i.e. the
    two-point BC is genuinely unenforceable, measured by the BVP residual
    ``‖Φ_rr p0 + Φ_rv ṗ0 - p1‖`` exceeding :data:`_PRIMER_BVP_RESIDUAL_TOL`.
    Hohmann's BC is in-range (residual ~1e-15) ⇒ trustworthy; an integer-rev
    phasing coast with mismatched impulse directions has a residual ~O(1) ⇒
    INDETERMINATE.

    Returns
    -------
    (pdot0, ill_conditioned, rcond):
        the recovered initial primer rate; a flag that ``Φ_rv`` was singular
        *and* the endpoint BC unenforceable (so the verdict must be reported
        INDETERMINATE); and the reciprocal condition number ``sigma_min/sigma_max``
        of ``Φ_rv``.
    """
    rhs = p1 - phi_rr @ p0
    u, s, vt = np.linalg.svd(phi_rv)
    smax = float(s[0])
    smin = float(s[-1])
    rcond = smin / smax if smax > 0.0 else 0.0
    singular = rcond < _PHI_RV_RCOND_FLOOR
    # Truncated-SVD (minimum-norm) pseudo-inverse: invert only the singular
    # values above the floor; this is the continuity-propagation primer on the
    # singular arc and is identical to the direct solve when full rank.
    s_inv = np.array(
        [1.0 / sv if sv > _PHI_RV_RCOND_FLOOR * smax else 0.0 for sv in s],
        dtype=np.float64,
    )
    pdot0 = (vt.T * s_inv) @ (u.T @ rhs)
    # On a singular Φ_rv the truncated solve cannot represent any BC component in
    # the null direction; a large residual there means the endpoint primer
    # direction is genuinely unreachable and the verdict must be INDETERMINATE.
    bvp_residual = float(np.linalg.norm(phi_rv @ pdot0 - rhs))
    ill_conditioned = singular and bvp_residual > _PRIMER_BVP_RESIDUAL_TOL
    return pdot0, ill_conditioned, rcond


@dataclass(frozen=True)
class CoastPrimerResult:
    """Primer-vector profile and verdict for a single coast arc.

    Attributes
    ----------
    coast_index:
        Zero-based index of the coast (= leg index in the impulse schedule).
    max_primer_magnitude:
        ``max_t |p(t)|`` over the coast (dimensionless).
    time_of_max_s:
        Time of the maximum, seconds from the coast start.
    duration_s:
        Coast duration, seconds.
    verdict:
        :class:`PrimerVerdict` for this coast.
    endpoint_magnitudes:
        ``(|p(t0)|, |p(t1)|)`` — both ≈ 1 by construction on a well-posed
        coast; reported so a caller can confirm the BCs were honoured. On an
        ``INDETERMINATE_ILL_CONDITIONED`` (near-integer-rev) coast the second
        entry may differ from 1: the endpoint BC is unenforceable there and the
        primer was propagated by continuity (see ``phi_rv_rcond``).
    phi_rv_rcond:
        Reciprocal condition number ``sigma_min/sigma_max`` of the coast STM
        block ``Φ_rv``. A value at/near ``0`` flags a near-integer-revolution
        (phasing-orbit) coast where the two-point primer inversion is singular
        (Saloglu, Taheri & Landau 2023 Sec. III.F).
    ill_conditioned:
        ``True`` when ``phi_rv_rcond`` fell below :data:`_PHI_RV_RCOND_FLOOR`
        and the continuity (truncated-SVD) primer-rate fallback was used; the
        verdict is then ``INDETERMINATE_ILL_CONDITIONED``.
    """

    coast_index: int
    max_primer_magnitude: float
    time_of_max_s: float
    duration_s: float
    verdict: PrimerVerdict
    endpoint_magnitudes: tuple[float, float]
    phi_rv_rcond: float = 1.0
    ill_conditioned: bool = False


@dataclass(frozen=True)
class PrimerDiagnostic:
    """Overall primer diagnostic for an impulsive N-burn schedule.

    Attributes
    ----------
    coasts:
        Per-coast results in schedule order.
    overall_verdict:
        Aggregated verdict. ``INDETERMINATE_ILL_CONDITIONED`` if any coast is
        ill-conditioned (a near-integer-rev coast whose verdict cannot be
        trusted dominates the schedule); else ``IMPROVABLE_ADD_IMPULSE`` if any
        coast is improvable; else ``OPTIMAL_NECESSARY_CONDITIONS_MET``.
    max_primer_magnitude:
        The largest per-coast ``max|p|`` across the whole schedule. NB: an
        ill-conditioned coast's ``max|p|`` is a continuity-fallback value, not a
        trustworthy optimality magnitude — read it together with
        ``any_ill_conditioned``.
    caveat:
        Standard necessary-conditions / multi-rev caveat string; the
        ill-conditioned (phasing-orbit) caveat is appended when any coast tripped
        the continuity fallback.
    any_ill_conditioned:
        ``True`` when at least one coast was a near-integer-rev singular arc.
    """

    coasts: tuple[CoastPrimerResult, ...]
    overall_verdict: PrimerVerdict
    max_primer_magnitude: float
    caveat: str
    any_ill_conditioned: bool = False


def primer_on_coast(
    r0: Vec3,
    v0: Vec3,
    p0_hat: Vec3,
    p1_hat: Vec3,
    duration_s: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    coast_index: int = 0,
    n_samples: int = 200,
    tol: float = 1.0e-6,
) -> CoastPrimerResult:
    """Primer profile on one coast pinned to unit ΔV directions at both ends.

    Solves the two-point BVP ``p(0) = p0_hat``, ``p(T) = p1_hat`` for the
    primer along the Kepler coast that starts at ``(r0, v0)`` and lasts
    ``duration_s``. The initial primer rate ``ṗ(0)`` is recovered from the
    arc STM (Lawden; the primer shares the variational equation), then
    ``|p(t)|`` is evaluated on a sample grid.

    Parameters
    ----------
    r0, v0:
        Heliocentric state at the start of the coast (just after the
        departure impulse), km and km/s.
    p0_hat, p1_hat:
        Unit ΔV directions at the bounding impulses (``Δv/|Δv|``). Normalised
        defensively here.
    duration_s:
        Coast duration, seconds (> 0).
    mu:
        Central-body gravitational parameter, km³/s².
    coast_index:
        Index recorded on the result.
    n_samples:
        Number of grid points for the ``|p(t)|`` scan (endpoints included).
    tol:
        Slack on the unit-magnitude test: a coast is ``IMPROVABLE`` only when
        ``max|p| > 1 + tol``.

    Returns
    -------
    CoastPrimerResult

    Notes
    -----
    On a near-integer-revolution (phasing-orbit) coast the position-velocity
    STM block ``Φ_rv`` is rank-deficient and the two-point inversion (Eq. 28 of
    Şaloğlu, Taheri & Landau 2023) is singular; ``ṗ(0)`` is then recovered by
    continuity (truncated-SVD minimum-norm; see :func:`_solve_primer_rate`) and
    the verdict is reported as ``INDETERMINATE_ILL_CONDITIONED`` — never a
    spurious IMPROVABLE — with the coast's ``Φ_rv`` reciprocal condition number
    on the result.
    """
    if duration_s <= 0.0:
        raise ValueError(f"coast duration must be positive, got {duration_s}")
    p0 = np.asarray(p0_hat, dtype=np.float64)
    p1 = np.asarray(p1_hat, dtype=np.float64)
    p0 = p0 / np.linalg.norm(p0)
    p1 = p1 / np.linalg.norm(p1)

    times, _ref, stms = _coast_stm(r0, v0, duration_s, mu, n_samples=n_samples)

    # End-state STM maps [p(0); ṗ(0)] -> [p(T); ṗ(T)].
    # p(T) = Φ_rr p(0) + Φ_rv ṗ(0)  =>  ṗ(0) = Φ_rv⁻¹ (p1 - Φ_rr p0).
    # Robust to near-integer-rev singularity of Φ_rv (continuity fallback).
    phi_end = stms[-1]
    phi_rr = phi_end[:3, :3]
    phi_rv = phi_end[:3, 3:]
    pdot0, ill_conditioned, rcond = _solve_primer_rate(phi_rr, phi_rv, p0, p1)

    state0 = np.concatenate([p0, pdot0])
    mags = np.empty(times.shape[0], dtype=np.float64)
    for k in range(times.shape[0]):
        pk = stms[k][:3, :] @ state0
        mags[k] = float(np.linalg.norm(pk))

    k_max = int(np.argmax(mags))
    max_mag = float(mags[k_max])
    if ill_conditioned:
        # Near-integer-rev coast: the endpoint BC is unenforceable, so a |p| > 1
        # excursion here is NOT a valid add-an-impulse signal (it may be the
        # iso-ΔV phasing-orbit degeneracy of Saloglu 2023 Sec. III.F). Report
        # INDETERMINATE rather than a false OPTIMAL/IMPROVABLE verdict.
        verdict = PrimerVerdict.INDETERMINATE_ILL_CONDITIONED
    elif max_mag > 1.0 + tol:
        verdict = PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    else:
        verdict = PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET
    return CoastPrimerResult(
        coast_index=coast_index,
        max_primer_magnitude=max_mag,
        time_of_max_s=float(times[k_max]),
        duration_s=float(duration_s),
        verdict=verdict,
        endpoint_magnitudes=(float(mags[0]), float(mags[-1])),
        phi_rv_rcond=float(rcond),
        ill_conditioned=bool(ill_conditioned),
    )


def diagnose_impulse_schedule(
    coast_states: list[tuple[Vec3, Vec3, float]],
    impulse_dirs: list[Vec3],
    *,
    mu: float = MU_SUN_KM3_S2,
    n_samples: int = 200,
    tol: float = 1.0e-6,
) -> PrimerDiagnostic:
    """Primer diagnostic for a fixed-time N-impulse schedule.

    Parameters
    ----------
    coast_states:
        One ``(r0, v0, duration_s)`` per coast, in schedule order. ``(r0,
        v0)`` is the heliocentric state at the *start* of the coast (just
        after that coast's departure impulse).
    impulse_dirs:
        Unit ΔV directions ``Δv_i/|Δv_i|`` at each impulse, length
        ``len(coast_states) + 1`` (an impulse bounds each coast end).
    mu, n_samples, tol:
        Passed to :func:`primer_on_coast`.

    Returns
    -------
    PrimerDiagnostic
    """
    n_coasts = len(coast_states)
    if len(impulse_dirs) != n_coasts + 1:
        raise ValueError(
            f"need {n_coasts + 1} impulse directions for {n_coasts} coasts, got {len(impulse_dirs)}"
        )
    coasts: list[CoastPrimerResult] = []
    for i, (r0, v0, dur) in enumerate(coast_states):
        coasts.append(
            primer_on_coast(
                r0,
                v0,
                impulse_dirs[i],
                impulse_dirs[i + 1],
                dur,
                mu=mu,
                coast_index=i,
                n_samples=n_samples,
                tol=tol,
            )
        )
    max_mag = max(c.max_primer_magnitude for c in coasts)
    any_ill = any(c.ill_conditioned for c in coasts)
    if any_ill:
        # A near-integer-rev coast whose BVP is singular cannot be ruled
        # OPTIMAL or IMPROVABLE; its uncertainty dominates the schedule verdict.
        overall = PrimerVerdict.INDETERMINATE_ILL_CONDITIONED
    elif any(c.verdict is PrimerVerdict.IMPROVABLE_ADD_IMPULSE for c in coasts):
        overall = PrimerVerdict.IMPROVABLE_ADD_IMPULSE
    else:
        overall = PrimerVerdict.OPTIMAL_NECESSARY_CONDITIONS_MET
    caveat = _NECESSARY_CONDITIONS_CAVEAT
    if any_ill:
        caveat = f"{caveat} {_ILL_CONDITIONED_CAVEAT}"
    return PrimerDiagnostic(
        coasts=tuple(coasts),
        overall_verdict=overall,
        max_primer_magnitude=max_mag,
        caveat=caveat,
        any_ill_conditioned=any_ill,
    )


def hohmann_primer_diagnostic(
    r1_km: float,
    r2_km: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    n_samples: int = 400,
    tol: float = 1.0e-6,
) -> PrimerDiagnostic:
    """Primer diagnostic for a coplanar-circular Hohmann transfer.

    Two tangential impulses connect circular orbits of radii ``r1_km`` and
    ``r2_km`` via the transfer ellipse (periapsis ``r1``, apoapsis ``r2``).
    Both ΔVs are tangential (along the local velocity), so the primer unit
    directions at the two impulses are the (anti-parallel across the
    half-revolution) velocity directions at periapsis and apoapsis.

    This is the SOURCED golden gate: the Hohmann transfer satisfies the
    primer necessary conditions iff the radius ratio is below the published
    Lawden/Marchal threshold (~11.94); above it ``max|p| > 1`` and a
    three-impulse (bi-elliptic) schedule is implied (Prussing & Conway).

    Parameters
    ----------
    r1_km, r2_km:
        Inner and outer circular radii, km. ``r2_km > r1_km`` (raise).
    mu, n_samples, tol:
        As :func:`primer_on_coast`.

    Returns
    -------
    PrimerDiagnostic with a single coast (the transfer arc).
    """
    if r2_km <= r1_km:
        raise ValueError(f"expected r2 > r1, got r1={r1_km}, r2={r2_km}")
    a_t = 0.5 * (r1_km + r2_km)
    # Transfer-ellipse periapsis state (start), in the ecliptic plane.
    r0 = np.array([r1_km, 0.0, 0.0], dtype=np.float64)
    v_peri = float(np.sqrt(mu * (2.0 / r1_km - 1.0 / a_t)))
    v0 = np.array([0.0, v_peri, 0.0], dtype=np.float64)
    # Half the ellipse period is the transfer time.
    tof = float(np.pi * np.sqrt(a_t**3 / mu))

    # Tangential ΔVs: at periapsis the departure burn raises speed in +y;
    # at apoapsis the arrival burn raises speed in -y (the spacecraft is on
    # the far side travelling in -y). Both impulses are along the local
    # velocity, so the primer end-directions are +y and -y.
    p0_hat = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    p1_hat = np.array([0.0, -1.0, 0.0], dtype=np.float64)

    return diagnose_impulse_schedule(
        [(r0, v0, tof)],
        [p0_hat, p1_hat],
        mu=mu,
        n_samples=n_samples,
        tol=tol,
    )


__all__ = [
    "CoastPrimerResult",
    "PrimerDiagnostic",
    "PrimerVerdict",
    "diagnose_impulse_schedule",
    "gravity_gradient",
    "hohmann_primer_diagnostic",
    "primer_on_coast",
]
