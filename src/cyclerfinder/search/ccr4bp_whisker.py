"""One-period amplification + segment-anchored discrete-QR/CLV whisker-direction
extraction for CCR4BP tori (`#691`, Stage B sub-task 3 of `#686`'s CCR4BP plan).

Scope (per `#691`'s own dispatch, two ordered steps)
------------------------------------------------------
1. **One-period amplification.** Cheap post-hoc state+STM extraction off an
   already-converged `#690` CCR4BP torus (:mod:`cyclerfinder.search.variational_ccr4bp_torus`):
   evaluate the torus's planar state at a phase ``(theta1, theta2)``, embed it into
   `#689`'s full 6-state ``(x, y, z, vx, vy, vz)`` (``z = vz = 0`` -- the CCR4BP torus
   is planar), and propagate state+STM for ONE full ``theta1`` cycle (one Ganymede
   synodic period, ``torus.period``) via :func:`cyclerfinder.core.ccr4bp.ccr4bp_stm_eom`.
   No new solver, no new corrector -- reuses `#689`'s STM equations verbatim
   (:func:`one_period_stm`), then reports eigenvalues/singular values of the
   resulting one-period map (:func:`amplification_diagnostics`).
2. **Perturbation-robustness diagnostic.** `#619` built a one-shot manifold-
   direction extractor for the QBCP torus (`search.variational_qbcp_torus.
   manifold_state_vec_pseudospectral`) that was later shown (`#646`) to be
   perturbation-UNSTABLE at the violently-unstable EM-L2 torus (a sub-invariance-
   residual base-point offset swings the extracted unstable eigenvector up to 88
   degrees, because the flow amplifies that offset ~2e4x over one period and the
   STM ends up linearized along a diverged trajectory). `#646`'s fix -- a
   segment-anchored discrete-QR / covariant-Lyapunov-vector (CLV) extraction
   (Benettin 1980 / Ginelli et al. 2007 / Dieci-Van Vleck) -- splits the one-period
   integration into ``n_segments`` short segments, RE-PROJECTS the base state onto
   the torus's own explicit functional representation
   (:func:`variational_ccr4bp_torus.evaluate_torus_state`) at every segment
   boundary (so the base-point error never compounds across the whole period, only
   across one short segment), composes the segment STMs with QR
   re-orthonormalization, and eigendecomposes the reconstructed one-period STM
   for the manifold-tangent EIGENvector (not the leading singular vector -- the
   correctly-computed one-period map can be strongly non-normal at some phases;
   see `#646`'s own module-level comment in ``variational_qbcp_torus.py`` for the
   full argument). This module adapts that SAME logic (:func:`manifold_direction_segmented_clv`)
   -- only the state-dimension/RHS backend changes: the CCR4BP's 6-state Cartesian
   rotating frame needs NO position/momentum (PV<->PM) conversion at all (simpler
   than the QBCP adapter `#646` built), so :func:`_segment_period_stm` is a
   directly simplified analogue of ``variational_qbcp_torus._segment_period_stm``.

   Both :func:`manifold_direction_one_shot` and :func:`manifold_direction_segmented_clv`
   accept an optional ``base_perturbation`` (a raw 6-vector offset applied to the
   extraction base point -- one-shot: the whole initial state; segmented: ONLY
   segment 0's initial condition, since every later segment boundary re-projects
   exactly onto the torus regardless of what happened during the previous
   segment). :func:`perturbation_swing_deg` runs the exact `#619`/`#646`
   perturbation-robustness test: extract the branch direction unperturbed and
   under a small random offset, report the angular swing between the two (a
   large swing under a sub-residual offset means the direction is NOT
   trustworthy at that phase; `#646` found <0.01 degrees for the segmented method
   vs up to 88 degrees for the one-shot method on the EM-L2 QBCP torus).

Mandatory positive control (this module's own tests)
------------------------------------------------------
Per this project's "verify against a known-good case before trusting a novel
one" discipline, the tests in ``tests/search/test_ccr4bp_whisker.py`` FIRST run
this exact one-period-STM + QR/CLV pipeline on the UNFORCED (``mu_gan = 0``)
Jupiter-Europa 3:4 resonant CR3BP periodic orbit -- where `#689` already proved
the CCR4BP EOM/STM reduce to the plain CR3BP exactly -- and cross-check the
result against a standard, independent CR3BP monodromy-matrix computation
(:func:`cyclerfinder.core.cr3bp.propagate` with ``with_stm=True`` over the
orbit's own natural period, NOT going through any of this module's or `#690`'s
code) before trusting the diagnostic on the physical-mass, genuinely forced
CCR4BP torus.

Discipline / scope
-------------------
Pure post-hoc diagnostic module: does NOT modify `#689`'s ``core/ccr4bp.py`` or
`#690`'s ``search/variational_ccr4bp_torus.py``. Explicitly NOT the
mesh-intersection heteroclinic search or manifold globalization (later `#686`
Stage-B sub-tasks), and NOT the full Kumar-Anderson-de la Llave-Gunter
bundle-solve whisker construction (only indicated if this module's own
perturbation diagnostic finds the cheap QR/CLV direction untrustworthy -- see
the `#691` dispatch note). No catalogue writeback -- a capability
characterization build, not a discovery result.

References
----------
* Kumar, B., Anderson, R. L., de la Llave, R., Gunter, B. (2021). AAS 21-651 /
  arXiv:2109.14815 (the CCR4BP whiskered-torus object class).
* Benettin, G., et al. (1980). Lyapunov exponents computation.
* Ginelli, F., et al. (2007). Covariant Lyapunov vectors.
* Dieci, L., & Van Vleck, E. S. Discrete-QR Lyapunov-exponent algorithms.
* `#646`'s ``search/variational_qbcp_torus.py`` segment-anchored CLV extraction
  (the method this module adapts, EOM-swapped onto the CCR4BP).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.ccr4bp as ccr4bp
from cyclerfinder.genome.qp_torus_manifold import local_stability
from cyclerfinder.search.variational_ccr4bp_torus import (
    CCR4BPTorusVariationalResult,
    evaluate_torus_state,
)

_N_FULL = 6  # embedded (x, y, z, vx, vy, vz) state the CCR4BP STM propagates.


def _embed6(u4: NDArray[np.float64]) -> NDArray[np.float64]:
    """Embed the CCR4BP torus's planar 4-state ``(x, y, vx, vy)`` into the full
    6-state ``(x, y, z, vx, vy, vz)`` with ``z = vz = 0`` -- `#689`'s planar-orbit
    convention (the CCR4BP torus corrector only ever represents planar motion)."""
    return np.array([u4[0], u4[1], 0.0, u4[2], u4[3], 0.0], dtype=np.float64)


# ---------------------------------------------------------------------------
# Step 1: one-shot one-period state + STM extraction and its characterization.
# ---------------------------------------------------------------------------


def one_period_stm(
    torus: CCR4BPTorusVariationalResult,
    theta1: float,
    theta2: float,
    *,
    base_perturbation: NDArray[np.float64] | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """One-shot state + one-period STM at a CCR4BP torus phase.

    Evaluates the torus's planar state at ``(theta1, theta2)``, embeds it into
    the full 6-state, optionally adds ``base_perturbation`` (a raw 6-vector
    offset -- the `#619`/`#646` perturbation-robustness-test mechanism: a small
    offset applied to the extraction point, below the torus's own invariance
    residual), and propagates state + STM FORWARD for one full ``torus.period``
    (one Ganymede-synodic ``theta1`` cycle) via `#689`'s
    :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp`, starting at absolute time
    ``t0 = (theta1 mod 2*pi) / torus.omega1``.

    Returns ``(state6_at_t0, stm_6x6)``. No new solver -- pure post-hoc
    integration off an already-converged torus (`#691` step 1).
    """
    t0 = (float(theta1) % (2.0 * np.pi)) / torus.omega1
    u4 = evaluate_torus_state(torus, float(theta1), float(theta2))
    state0 = _embed6(u4)
    if base_perturbation is not None:
        state0 = state0 + np.asarray(base_perturbation, dtype=np.float64)
    arc = ccr4bp.propagate_ccr4bp(
        torus.system, state0, torus.period, with_stm=True, t0=t0, rtol=rtol, atol=atol
    )
    assert arc.stm is not None
    return state0, arc.stm


@dataclass(frozen=True)
class OnePeriodAmplification:
    """Eigen/singular-value characterization of a torus point's one-period STM.

    ``eigenvalues`` are ALL 6 (possibly complex) eigenvalues of the one-period
    STM, unsorted (as returned by ``numpy.linalg.eigvals``); ``singular_values``
    are its 6 singular values, descending. ``lam_u`` (``> 1 + hyperbolicity_tol``)
    / ``lam_s`` (``< 1 - hyperbolicity_tol``) are the real hyperbolic eigenpair
    picked by :func:`cyclerfinder.genome.qp_torus_manifold.local_stability`
    (``None`` if the point is locally non-hyperbolic under this measure).
    ``separation_ratio`` is ``|largest eigenvalue| / |second-largest eigenvalue|``
    (by magnitude) -- how well the dominant (generically unstable) direction is
    separated from the rest of the spectrum; large means a clean, well-isolated
    unstable direction, near-1 means the spectrum is crowded near the dominant
    magnitude (e.g. two comparably-unstable real eigenvalues, or a complex pair
    with magnitude close to the real one).
    """

    state6: NDArray[np.float64]
    stm: NDArray[np.float64]
    eigenvalues: NDArray[np.complex128]
    singular_values: NDArray[np.float64]
    lam_u: float | None
    lam_s: float | None
    vec_u: NDArray[np.float64] | None
    vec_s: NDArray[np.float64] | None
    separation_ratio: float


def amplification_diagnostics(
    state6: NDArray[np.float64],
    stm: NDArray[np.float64],
    *,
    hyperbolicity_tol: float = 1e-4,
) -> OnePeriodAmplification:
    """Characterize a one-period STM's amplification (`#691` step 1 report).

    Reuses :func:`cyclerfinder.genome.qp_torus_manifold.local_stability`
    (system-agnostic: it only needs ``(state, stm)``) for the real hyperbolic
    eigenpair, and adds the full eigenvalue/singular-value spectrum plus a
    dominant-direction separation ratio.
    """
    eigvals = np.asarray(np.linalg.eigvals(stm), dtype=np.complex128)
    svals = np.asarray(np.linalg.svd(stm, compute_uv=False), dtype=np.float64)
    stab = local_stability(state6, stm, hyperbolicity_tol=hyperbolicity_tol)
    mags = np.sort(np.abs(eigvals))[::-1]
    separation = float(mags[0] / mags[1]) if mags[1] > 0.0 else float("inf")
    return OnePeriodAmplification(
        state6=state6,
        stm=stm,
        eigenvalues=eigvals,
        singular_values=svals,
        lam_u=stab.lam_u,
        lam_s=stab.lam_s,
        vec_u=stab.vec_u,
        vec_s=stab.vec_s,
        separation_ratio=separation,
    )


def manifold_direction_one_shot(
    torus: CCR4BPTorusVariationalResult,
    theta1: float,
    theta2: float,
    branch: str,
    ref_vec: NDArray[np.float64] | None,
    *,
    base_perturbation: NDArray[np.float64] | None = None,
    hyperbolicity_tol: float = 1e-4,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    """One-shot (GMOS/`#619`-convention) manifold direction at a torus phase.

    The "naive" extraction :func:`manifold_direction_segmented_clv` is checked
    against: a SINGLE one-period STM eigendecomposition, no segmentation, no
    re-projection. ``branch`` selects ``vec_u`` (unstable) or ``vec_s`` (stable);
    ``ref_vec`` fixes the sign by dot-product continuity. Returns ``None`` if no
    real hyperbolic eigenpair of the requested branch exists at this phase.
    """
    if branch not in ("unstable", "stable"):
        raise ValueError(f"branch must be 'unstable' or 'stable', got {branch!r}")
    state0, stm = one_period_stm(
        torus, theta1, theta2, base_perturbation=base_perturbation, rtol=rtol, atol=atol
    )
    stab = local_stability(state0, stm, hyperbolicity_tol=hyperbolicity_tol)
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None
    vec = np.asarray(vec, dtype=np.float64)
    vec = vec / np.linalg.norm(vec)
    if ref_vec is not None and float(np.dot(vec, ref_vec)) < 0.0:
        vec = -vec
    return state0, vec


# ---------------------------------------------------------------------------
# Step 2: segment-anchored discrete-QR / CLV extraction (`#646` adapted).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SegmentedCLVDiagnostics:
    """Conditioning diagnostics for a segment-anchored CLV extraction (`#646`
    adapted to the CCR4BP; see ``variational_qbcp_torus.SegmentedCLVDiagnostics``
    for the full field-by-field rationale, identical here).

    ``n_segments`` segments were used; ``lyapunov_exponents`` are the discrete-QR
    Lyapunov exponents (log growth / period, descending, 6 of them);
    ``per_segment_leading_growth`` is the leading orthonormal-frame growth factor
    per segment (should stay modest, NOT match the full one-period amplification
    -- that is the whole point of segmenting); ``lam_u``/``lam_s`` are the
    reconstructed one-period STM's real hyperbolic eigenvalues;
    ``nonnormality_ratio`` is ``sigma_max(P) / |lam_u|``; ``qr_vs_direct_stm_relerr``
    is the relative Frobenius error between the QR-reconstructed and directly-
    composed one-period STM (self-consistency check on the QR accumulation).
    """

    n_segments: int
    lyapunov_exponents: NDArray[np.float64]
    per_segment_leading_growth: NDArray[np.float64]
    lam_u: float | None
    lam_s: float | None
    nonnormality_ratio: float
    qr_vs_direct_stm_relerr: float


def _reproject6(
    torus: CCR4BPTorusVariationalResult, theta1: float, theta2: float
) -> NDArray[np.float64]:
    """Evaluate the torus's planar state at ``(theta1, theta2)`` and embed it into
    the full 6-state -- the `#646` re-projection, adapted (no PM/PV conversion
    needed for the CCR4BP's plain Cartesian rotating-frame state)."""
    return _embed6(evaluate_torus_state(torus, float(theta1), float(theta2)))


def _segment_period_stm(
    torus: CCR4BPTorusVariationalResult,
    theta1: float,
    theta2: float,
    n_segments: int,
    rtol: float,
    atol: float,
    *,
    base_perturbation: NDArray[np.float64] | None = None,
) -> (
    tuple[
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
    ]
    | None
):
    """Segment-anchored, torus-re-projected composition of the one-period STM.

    Direct simplification of ``variational_qbcp_torus._segment_period_stm``: the
    CCR4BP's Cartesian rotating-frame state needs no PV<->PM conversion, so each
    segment is a single :func:`cyclerfinder.core.ccr4bp.propagate_ccr4bp` call
    with ``with_stm=True``. Returns ``(base_state_t0, p_qr, p_direct,
    per_segment_leading_growth, log_growth)``, or ``None`` if any segment
    propagation fails. ``base_perturbation`` (if given) is added ONLY to segment
    0's initial condition -- every later segment boundary re-projects exactly
    onto the torus regardless (the `#619`/`#646` perturbation-robustness-test
    mechanism: does a sub-residual offset at the very start still get amplified
    to the same degree, given the periodic re-projection resets the trajectory
    back onto the torus after each short segment).
    """
    omega1, omega2 = torus.omega1, torus.omega2
    period = torus.period
    t0 = (float(theta1) % (2.0 * np.pi)) / omega1
    dt = period / n_segments

    q_acc = np.eye(_N_FULL)
    r_prod = np.eye(_N_FULL)
    p_direct = np.eye(_N_FULL)
    log_growth = np.zeros(_N_FULL)
    seg_lead = np.empty(n_segments, dtype=np.float64)
    base_state_t0 = np.empty(_N_FULL)

    for k in range(n_segments):
        sk = t0 + k * dt
        tl_k = theta1 + omega1 * (sk - t0)
        tt_k = theta2 + omega2 * (sk - t0)
        state_k = _reproject6(torus, tl_k, tt_k)
        if k == 0:
            if base_perturbation is not None:
                state_k = state_k + np.asarray(base_perturbation, dtype=np.float64)
            base_state_t0 = state_k
        try:
            arc = ccr4bp.propagate_ccr4bp(
                torus.system, state_k, dt, with_stm=True, t0=sk, rtol=rtol, atol=atol
            )
        except RuntimeError:
            return None
        assert arc.stm is not None
        m_k = arc.stm
        p_direct = m_k @ p_direct
        # QR re-orthonormalization of the running product (Dieci-Van Vleck):
        # Y = M_k Q_{k-1} = Q_k R_k; R_prod <- R_k R_prod so P = Q_K R_prod.
        y = m_k @ q_acc
        q_acc, r_k = np.linalg.qr(y)
        sign = np.sign(np.diag(r_k))
        sign[sign == 0.0] = 1.0
        q_acc = q_acc * sign
        r_k = (r_k.T * sign).T
        r_prod = r_k @ r_prod
        seg_lead[k] = abs(r_k[0, 0])
        log_growth += np.log(np.abs(np.diag(r_k)))

    p_qr = q_acc @ r_prod
    return base_state_t0, p_qr, p_direct, seg_lead, log_growth


def manifold_direction_segmented_clv(
    torus: CCR4BPTorusVariationalResult,
    theta1: float,
    theta2: float,
    branch: str,
    ref_vec: NDArray[np.float64] | None,
    *,
    n_segments: int = 16,
    hyperbolicity_tol: float = 1e-4,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    base_perturbation: NDArray[np.float64] | None = None,
    return_diagnostics: bool = False,
) -> (
    tuple[NDArray[np.float64], NDArray[np.float64]]
    | tuple[NDArray[np.float64], NDArray[np.float64], SegmentedCLVDiagnostics]
    | None
):
    """Segment-anchored discrete-QR / CLV manifold direction at a CCR4BP torus
    phase (`#646` extraction, adapted to the CCR4BP -- `#691` step 2).

    Drop-in analogue of :func:`manifold_direction_one_shot` with the same
    ``(state6, unit_eigenvector)`` contract (add ``return_diagnostics=True`` for a
    trailing :class:`SegmentedCLVDiagnostics`). See the module docstring and
    ``variational_qbcp_torus``'s own module-level comment for the full rationale.
    ``branch`` is ``"unstable"`` or ``"stable"``; ``ref_vec`` fixes the sign by
    dot-product continuity. Returns ``None`` if a segment propagation fails or no
    real hyperbolic eigenpair of the requested branch exists.
    """
    if branch not in ("unstable", "stable"):
        raise ValueError(f"branch must be 'unstable' or 'stable', got {branch!r}")
    if n_segments < 1:
        raise ValueError(f"n_segments must be >= 1, got {n_segments}")

    composed = _segment_period_stm(
        torus, theta1, theta2, n_segments, rtol, atol, base_perturbation=base_perturbation
    )
    if composed is None:
        return None
    base_state6, p_qr, p_direct, seg_lead, log_growth = composed

    stab = local_stability(base_state6, p_qr, hyperbolicity_tol=hyperbolicity_tol)
    vec = stab.vec_u if branch == "unstable" else stab.vec_s
    if vec is None:
        return None
    vec = np.asarray(vec, dtype=np.float64)
    vec = vec / np.linalg.norm(vec)
    if ref_vec is not None and float(np.dot(vec, ref_vec)) < 0.0:
        vec = -vec

    if not return_diagnostics:
        return base_state6, vec

    sigma_max = float(np.linalg.svd(p_qr, compute_uv=False)[0])
    nonnormality = sigma_max / abs(stab.lam_u) if stab.lam_u else float("inf")
    relerr = float(np.linalg.norm(p_qr - p_direct) / (np.linalg.norm(p_direct) + 1e-300))
    diag = SegmentedCLVDiagnostics(
        n_segments=n_segments,
        lyapunov_exponents=np.sort(log_growth / torus.period)[::-1],
        per_segment_leading_growth=seg_lead,
        lam_u=stab.lam_u,
        lam_s=stab.lam_s,
        nonnormality_ratio=nonnormality,
        qr_vs_direct_stm_relerr=relerr,
    )
    return base_state6, vec, diag


# ---------------------------------------------------------------------------
# Perturbation-robustness diagnostic (`#619`/`#646`'s go/no-go test).
# ---------------------------------------------------------------------------


def perturbation_swing_deg(
    torus: CCR4BPTorusVariationalResult,
    theta1: float,
    theta2: float,
    branch: str,
    *,
    method: str,
    perturbation_size: float = 1e-4,
    n_segments: int = 16,
    seed: int = 0,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> float | None:
    """Angular swing (degrees) of the extracted branch direction under a small
    random base-point perturbation -- the exact `#619`/`#646` perturbation-
    robustness test.

    Draws a random unit direction (seeded by ``seed``) scaled to
    ``perturbation_size`` and adds it as ``base_perturbation`` to
    :func:`manifold_direction_one_shot` (``method="one_shot"``) or
    :func:`manifold_direction_segmented_clv` (``method="segmented"``); compares
    the perturbed direction to the unperturbed one via ``arccos(|dot|)`` (the
    absolute value avoids a trivial 180-degree sign-flip artifact -- both
    extractors resolve eigenvector sign only via an optional ``ref_vec``, which
    is not supplied here). Returns ``None`` if either extraction fails (no real
    hyperbolic eigenpair at this phase).
    """
    rng = np.random.default_rng(seed)
    delta = rng.normal(size=_N_FULL)
    delta = delta / np.linalg.norm(delta) * perturbation_size

    base_vec: NDArray[np.float64] | None
    pert_vec: NDArray[np.float64] | None
    if method == "one_shot":
        base_one_shot = manifold_direction_one_shot(
            torus, theta1, theta2, branch, None, rtol=rtol, atol=atol
        )
        pert_one_shot = manifold_direction_one_shot(
            torus, theta1, theta2, branch, None, base_perturbation=delta, rtol=rtol, atol=atol
        )
        base_vec = None if base_one_shot is None else base_one_shot[1]
        pert_vec = None if pert_one_shot is None else pert_one_shot[1]
    elif method == "segmented":
        base_seg = manifold_direction_segmented_clv(
            torus, theta1, theta2, branch, None, n_segments=n_segments, rtol=rtol, atol=atol
        )
        pert_seg = manifold_direction_segmented_clv(
            torus,
            theta1,
            theta2,
            branch,
            None,
            n_segments=n_segments,
            base_perturbation=delta,
            rtol=rtol,
            atol=atol,
        )
        base_vec = None if base_seg is None else base_seg[1]
        pert_vec = None if pert_seg is None else pert_seg[1]
    else:
        raise ValueError(f"method must be 'one_shot' or 'segmented', got {method!r}")

    if base_vec is None or pert_vec is None:
        return None
    dot = abs(float(np.dot(base_vec, pert_vec)))
    return float(np.degrees(np.arccos(min(1.0, dot))))
