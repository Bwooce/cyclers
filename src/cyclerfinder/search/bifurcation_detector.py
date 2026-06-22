"""Period-multiplying and saddle-center bifurcation detector for CR3BP families
(#266 Phase 2, extended #347 Phase 1).

When a one-parameter family of periodic orbits is continued through parameter
space, the Floquet multipliers (eigenvalues of the monodromy matrix) sweep
through the complex unit circle. **A bifurcation point is a parameter value
where a multiplier crosses a primitive k-th root of unity**: at that point a
new family with k times the parent period branches off. This module identifies
such crossings in an ordered list of converged family members.

The detector is the *minimum-viable* substrate for tulip-orbit family hunting
(#266 Phase 2): given a continuation along the NRHO family it flags every
adjacent pair where a multiplier crossed a root-of-unity, narrowing the search
for new petal-count families to a bracket the corrector can later refine.

The k=1 saddle-center case (extension #347 Phase 1, RTR2026 + Hamiltonian
bifurcation theory) is the bifurcation at lambda=+1 with codimension 1: a
non-trivial Floquet multiplier coalesces with the trivial pair at +1. In the
CR3BP Hamiltonian flow this manifests as a complex-conjugate pair on the unit
circle coalescing on the real axis at +1 and splitting into a real reciprocal
pair (one branch slightly above +1, one below). The pre-existing
:func:`detect_period_multiplying` excludes k=1 (the trivial-pair degeneracy
swamps the signal); :func:`detect_saddle_center_bracket` is the k=1
specialisation. It identifies the bracket where the two "secondary" non-
trivial eigenvalues transition from complex-on-unit-circle (pre-bifurcation)
to real-near-+1 (post-bifurcation).

Discipline:

  * The detector returns **brackets**, never asserts a single bifurcation
    point. Refining the bracket is the corrector's job (Phase 3).
  * "Crossing" is identified by a sign flip of ``dist - tol`` across adjacent
    family members; the detector itself does not interpolate.
  * Only adjacent-member crossings are reported. A multiplier that crosses a
    root and then crosses back inside a single tolerance band is reported as
    two brackets (the discrete signal); the caller resolves.

Style cross-reference: the topological-flag pattern follows
:func:`cyclerfinder.search.binary_star_search.winding_topology` -- compute a
labelling invariant per member, then derive the family-level signal from the
sequence of labels.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp

# ---------------------------------------------------------------------------
# Public types.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FamilyMember:
    """One converged orbit along a continued family.

    Attributes
    ----------
    label :
        Caller-supplied identifier (e.g. continuation index, parameter value).
        Used to label the bracket in :class:`BifurcationPoint`.
    state0 :
        IC at the family member.
    period :
        Full nondim period of the member.
    mu :
        Mass parameter the member is referenced to. Stored per-member because
        :func:`monodromy` accepts a CR3BPSystem at compute time and the family
        may be continued in mu.
    parameter :
        The continuation parameter value at this member (e.g. perilune
        altitude, Jacobi constant). Optional -- when ``None`` the family is
        labelled by ``label`` only.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    mu: float
    parameter: float | None = None


@dataclass(frozen=True)
class BifurcationPoint:
    """A bracketed period-multiplying bifurcation between two family members.

    Attributes
    ----------
    k :
        The period-multiplying integer (``k = 2`` -> period doubling,
        ``k = 3`` -> tripling, ...).
    members :
        The two adjacent FamilyMembers between which a multiplier crossed a
        primitive k-th root of unity.
    eig_before, eig_after :
        The multiplier most closely matching a primitive k-th root of unity at
        each member.
    dist_before, dist_after :
        The associated distances to the nearest such root. The crossing is
        identified by ``(dist_before - tol)`` and ``(dist_after - tol)``
        having opposite signs.
    tol :
        The tolerance used when judging proximity to a root of unity.
    """

    k: int
    members: tuple[FamilyMember, FamilyMember]
    eig_before: complex
    eig_after: complex
    dist_before: float
    dist_after: float
    tol: float
    # Extra diagnostic context -- kept structural rather than free-form for the
    # caller to log without parsing a string.
    extras: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Monodromy + Floquet multipliers.
# ---------------------------------------------------------------------------


def monodromy(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    stm_mode: cr3bp.StmMode = "variable",
) -> NDArray[np.float64]:
    """Return the 6x6 monodromy matrix for one full period.

    Integrates the standard (state, STM) augmented variational EOM (see
    :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom`) for time ``period`` from
    ``state0``. The returned matrix is the STM ``Phi(T)`` at the end of one
    period; this is the monodromy whose eigenvalues are the Floquet
    multipliers.

    Parameters
    ----------
    system :
        CR3BP system; only ``system.mu`` is read.
    state0 :
        6-vector IC at the start of the period.
    period :
        Full nondim period.
    rtol, atol :
        Integrator tolerances. Defaults match the rest of the cr3bp module.
    stm_mode :
        Forwarded to :func:`cyclerfinder.core.cr3bp.propagate`. ``"variable"``
        (default; backward-compatible) is the legacy variable-step variational
        path. ``"fixed_path"`` is the Pellegrini-Russell 2016 mitigation
        (Pellegrini & Russell, JGCD, DOI 10.2514/1.G001920, §III.A.2 + Conc. 2)
        recommended for monodromy at clustered eigenvalues (saddle-center
        cluster point, λ approaching primitive roots of unity), where the
        eq. 17 step-size IC-dependence contaminates the multipliers under the
        variable-step variational path. See the ``propagate`` docstring for
        the full derivation and the cost ratio.

    Returns
    -------
    NDArray (6, 6) :
        The monodromy matrix.

    Raises
    ------
    RuntimeError
        Propagated from the integrator if it fails.
    ValueError
        If ``state0`` is not a 6-vector, or ``stm_mode`` is not supported.

    Reference
    ---------
    Pellegrini, E. & Russell, R.P. (2016), "On the Computation and Accuracy of
    Trajectory State Transition Matrices", *Journal of Guidance, Control, and
    Dynamics*, DOI 10.2514/1.G001920.
    """
    s0 = np.asarray(state0, dtype=np.float64)
    if s0.shape != (6,):
        raise ValueError(f"monodromy: state0 must be a 6-vector, got shape {s0.shape}")
    arc = cr3bp.propagate(
        system, s0, float(period), with_stm=True, rtol=rtol, atol=atol, stm_mode=stm_mode
    )
    assert arc.stm is not None
    return np.asarray(arc.stm, dtype=np.float64)


def floquet_multipliers(monodromy_matrix: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Return the monodromy's eigenvalues sorted by magnitude descending.

    Uses Schur decomposition to re-orthonormalize the eigenvector basis
    (Gram-Schmidt) so clustered/degenerate multipliers don't collapse the basis,
    fixing numerical conditioning per Doedel 1991.

    The CR3BP monodromy has reciprocal-pair structure (Liouville / time-reversal
    symmetry): if ``lambda`` is an eigenvalue so is ``1/lambda``. One pair is
    always ``(1, 1)`` (energy/time-translation). The detector below uses these
    multipliers, sorted, without separating the trivial pair -- the
    period-multiplying check is on every multiplier individually.
    """
    import scipy.linalg

    t_mat, _ = scipy.linalg.schur(monodromy_matrix, output="complex")
    eigs = np.diag(t_mat)
    order = np.argsort(-np.abs(eigs))
    return np.asarray(eigs[order], dtype=np.complex128)


# ---------------------------------------------------------------------------
# Period-multiplying detection.
# ---------------------------------------------------------------------------


def _primitive_roots_of_unity(k: int) -> list[complex]:
    """Return the primitive k-th roots of unity (j coprime with k)."""
    return [
        complex(math.cos(2 * math.pi * j / k), math.sin(2 * math.pi * j / k))
        for j in range(1, k)
        if math.gcd(j, k) == 1
    ]


def _nearest_kth_root_distance(eig: complex, k: int) -> tuple[float, complex]:
    """Return (min distance from eig to a primitive k-th root, the root itself)."""
    roots = _primitive_roots_of_unity(k)
    if not roots:  # k = 1 is excluded externally; guard anyway
        return float("inf"), complex(1.0)
    dists = [abs(eig - r) for r in roots]
    idx = int(np.argmin(dists))
    return float(dists[idx]), roots[idx]


def detect_period_multiplying(
    eigs: NDArray[np.complex128],
    *,
    k_max: int = 6,
    tol: float = 1e-3,
) -> list[tuple[int, complex, float]]:
    """For each multiplier, report the primitive k-th-root-of-unity it is nearest to.

    For each eigenvalue and each k in ``[2, k_max]``, compute the distance to
    the nearest primitive k-th root of unity. Emit ``(k, lambda, distance)``
    whenever ``distance < tol`` -- a multiplier sitting on (within tolerance)
    a primitive k-th root of unity signals a period-k bifurcation **at this
    family member**.

    The list is sorted by distance ascending: the closest match first. The
    same multiplier may appear under multiple ``k`` if it happens to lie near
    several roots; the caller's responsibility to pick the relevant one (the
    primitive smallest-k root is the standard convention).

    Parameters
    ----------
    eigs :
        Array of Floquet multipliers (complex, length N).
    k_max :
        Largest period-multiplying integer to check. Defaults to ``6`` --
        the genome literature rarely catalogues k > 6 bifurcations.
    tol :
        Distance tolerance below which a multiplier is "on" a root.

    Returns
    -------
    list of (k, lambda, distance) :
        Sorted by distance ascending. Empty when no multiplier is within
        ``tol`` of any primitive k-th root for ``2 <= k <= k_max``.
    """
    if k_max < 2:
        raise ValueError(f"detect_period_multiplying: k_max must be >= 2, got {k_max}")
    if tol <= 0.0 or not math.isfinite(tol):
        raise ValueError(f"detect_period_multiplying: tol must be positive finite, got {tol}")
    out: list[tuple[int, complex, float]] = []
    for eig in eigs:
        for k in range(2, k_max + 1):
            dist, _ = _nearest_kth_root_distance(complex(eig), k)
            if dist < tol:
                out.append((k, complex(eig), dist))
    out.sort(key=lambda x: x[2])
    return out


def scan_family_for_bifurcations(
    seeds: Sequence[FamilyMember],
    *,
    primary: str = "Earth",
    secondary: str = "Moon",
    k_max: int = 6,
    tol: float = 1e-2,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[BifurcationPoint]:
    """Scan an ordered family for period-multiplying bifurcations.

    For each consecutive pair ``(member_i, member_{i+1})`` and each ``k`` in
    ``[2, k_max]``:

      1. Compute the Floquet multipliers at both members.
      2. For each multiplier at each member, find the nearest primitive k-th
         root of unity and its distance ``d``.
      3. Take the minimum distance ``d_min`` over multipliers at each member.
      4. If ``d_min`` crosses ``tol`` between the two members (one above, one
         below), emit a :class:`BifurcationPoint` bracket for this ``k``.

    This is the discrete crossing signal; the bracket can be refined later by
    a corrector that interpolates parameter -> multiplier and zeroes
    ``d(parameter) - tol`` (Phase 3 follow-up).

    Parameters
    ----------
    seeds :
        Ordered family members. **Order matters**: the bracket is between
        adjacent members in the sequence the caller supplied.
    primary, secondary :
        Labels passed through to the CR3BPSystem; informational only --
        ``monodromy`` only reads ``system.mu``, which is set from each
        member's ``mu`` attribute.
    k_max, tol, rtol, atol :
        Forwarded to the detection routine and integrators.

    Returns
    -------
    list of BifurcationPoint :
        Adjacent-member brackets, one per (k, pair) signal. May be empty.
    """
    if len(seeds) < 2:
        return []

    # Per-member eigenvalue + (k -> nearest-root, distance) cache.
    per_member_eigs: list[NDArray[np.complex128]] = []
    per_member_k_summary: list[dict[int, tuple[complex, float]]] = []
    for mem in seeds:
        sysm = cr3bp.CR3BPSystem(
            mu=mem.mu,
            primary=primary,
            secondary=secondary,
            l_km=1.0,
            t_s=1.0,  # placeholders; monodromy() only reads mu
        )
        mat = monodromy(sysm, mem.state0, mem.period, rtol=rtol, atol=atol)
        eigs = floquet_multipliers(mat)
        per_member_eigs.append(eigs)
        k_summary: dict[int, tuple[complex, float]] = {}
        for k in range(2, k_max + 1):
            # Best (closest) match across the 6 multipliers at this k.
            best_dist = float("inf")
            best_root = complex(1.0)
            best_eig = complex(1.0)
            for eig in eigs:
                d, root = _nearest_kth_root_distance(complex(eig), k)
                if d < best_dist:
                    best_dist = d
                    best_root = root
                    best_eig = complex(eig)
            k_summary[k] = (best_eig, best_dist)
            del best_root  # not currently emitted -- kept available for future
        per_member_k_summary.append(k_summary)

    brackets: list[BifurcationPoint] = []
    for i in range(len(seeds) - 1):
        a, b = seeds[i], seeds[i + 1]
        for k in range(2, k_max + 1):
            eig_a, dist_a = per_member_k_summary[i][k]
            eig_b, dist_b = per_member_k_summary[i + 1][k]
            # Signed proximity to the tolerance band; sign flip => crossing.
            sa = dist_a - tol
            sb = dist_b - tol
            # Sign flip across the tolerance band -> a crossing. Exact-on-tol
            # is treated as a crossing of the boundary.
            crossed = (sa == 0.0 or sb == 0.0) or (sa > 0.0) != (sb > 0.0)
            if crossed:
                pa = a.parameter if a.parameter is not None else float("nan")
                pb = b.parameter if b.parameter is not None else float("nan")
                brackets.append(
                    BifurcationPoint(
                        k=k,
                        members=(a, b),
                        eig_before=eig_a,
                        eig_after=eig_b,
                        dist_before=dist_a,
                        dist_after=dist_b,
                        tol=tol,
                        extras={"param_before": pa, "param_after": pb},
                    )
                )
    return brackets


# ---------------------------------------------------------------------------
# Saddle-center / pitchfork (k=1) detection (#347 Phase 1).
# ---------------------------------------------------------------------------


def _classify_secondary_pair(eigs: NDArray[np.complex128]) -> tuple[str, complex, complex]:
    """Classify the two "secondary" Floquet multipliers of a planar CR3BP orbit.

    The CR3BP monodromy has 6 eigenvalues with reciprocal-pair structure:

      * one trivial pair near (1, 1) (energy + time-translation; split only by
        integrator round-off);
      * one "primary" real reciprocal pair (lambda_max, 1/lambda_max) carrying
        the orbit's dominant unstable manifold (e.g. ~ 2.5e5 / 4e-6 for the
        (3,2) Earth-Moon C32 anchor);
      * one "secondary" pair: either a complex conjugate pair on the unit
        circle (pre-bifurcation) OR a real reciprocal pair near +1 (post-
        bifurcation). This is the pair the saddle-center detector tracks.

    The classification picks the secondary pair by excluding (a) the 2
    eigenvalues closest to +1 (trivial), and (b) the eigenvalues with the
    largest |log|lambda|| (primary saddle pair). The remaining 2 are the
    secondary pair.

    Returns
    -------
    (kind, lam_a, lam_b) :
        ``kind`` is one of "complex_unit_circle", "real_near_one", or
        "real_far". ``(lam_a, lam_b)`` are the two secondary eigenvalues; for
        a complex pair, ``lam_a`` has the positive imaginary component.

    Notes
    -----
    The classification is local-state-only — it makes no continuation
    inference. The bifurcation BRACKET is produced by comparing the
    classification across adjacent family members.

    "real_near_one" means both secondary eigenvalues are real AND within 0.1
    of +1 (the post-saddle-center signature; the bifurcated pair starts at
    (+1, +1) and slowly separates as the family advances). "real_far" is for
    orbits where the secondary pair is real but well separated — that can
    happen far past the bifurcation. "complex_unit_circle" requires both
    eigenvalues complex with non-zero imaginary and |lambda - 1| > 1e-6.
    """
    eigs_c = np.asarray(eigs, dtype=np.complex128)
    if eigs_c.shape[0] < 4:
        raise ValueError(f"_classify_secondary_pair: need >= 4 eigenvalues, got {eigs_c.shape[0]}")
    # Step 1: exclude the 2 trivial-pair eigenvalues (closest to +1 by |lam - 1|).
    dists_to_one = np.abs(eigs_c - 1.0)
    triv_idx = np.argsort(dists_to_one)[:2]
    triv_set = set(int(i) for i in triv_idx)
    remaining = [(i, eigs_c[i]) for i in range(eigs_c.shape[0]) if i not in triv_set]
    # Step 2: among remaining, exclude the primary saddle pair by largest |log|lambda||.
    log_mags = [(i, abs(math.log(max(abs(e), 1e-300)))) for i, e in remaining]
    log_mags.sort(key=lambda t: -t[1])
    primary_indices = {log_mags[0][0], log_mags[1][0]}
    secondary = [(i, e) for i, e in remaining if i not in primary_indices]
    if len(secondary) != 2:
        # Fallback: take the 2 remaining whose |lambda| is closest to 1.
        secondary = sorted(remaining, key=lambda t: abs(abs(t[1]) - 1.0))[:2]
    lam_a = complex(secondary[0][1])
    lam_b = complex(secondary[1][1])
    # Conventional ordering: positive-imag first.
    if lam_a.imag < lam_b.imag:
        lam_a, lam_b = lam_b, lam_a
    # Classify.
    imag_a = abs(lam_a.imag)
    imag_b = abs(lam_b.imag)
    is_complex = imag_a > 1e-6 and imag_b > 1e-6
    if is_complex:
        return "complex_unit_circle", lam_a, lam_b
    near_one_a = abs(lam_a - 1.0) < 0.1
    near_one_b = abs(lam_b - 1.0) < 0.1
    if near_one_a and near_one_b:
        return "real_near_one", lam_a, lam_b
    return "real_far", lam_a, lam_b


def _deflated_determinant(eigs: NDArray[np.complex128]) -> float:
    """Compute the deflated determinant of (M - I) as a scalar test function.

    The determinant of the augmented Jacobian (Doedel 1991 Part I Eq. 2.9) is
    proportional to the product of (lambda - 1) for the non-trivial eigenvalues.
    This scalar test function changes sign across any true bifurcation but does
    not change sign across a simple fold, providing a robust dual-detector check
    for the saddle-center (k=1) bracket.
    """
    dists = np.abs(eigs - 1.0)
    # Exclude the 2 eigenvalues closest to +1 (the trivial pair)
    triv_idx = set(int(i) for i in np.argsort(dists)[:2])
    non_triv = [eigs[i] for i in range(len(eigs)) if i not in triv_idx]
    p = np.prod([lam - 1.0 for lam in non_triv])
    return float(p.real)


def detect_saddle_center_bracket(
    seeds: Sequence[FamilyMember],
    *,
    primary: str = "Earth",
    secondary: str = "Moon",
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[BifurcationPoint]:
    """Scan an ordered family for the saddle-center / pitchfork (k=1) bifurcation.

    For each consecutive pair ``(member_i, member_{i+1})``:

      1. Compute the monodromy + Floquet multipliers at both members.
      2. Classify the secondary eigenvalue pair (see
         :func:`_classify_secondary_pair`).
      3. If the classification TRANSITIONS from ``complex_unit_circle`` at
         member_i to ``real_near_one`` at member_{i+1} (or vice versa), emit a
         :class:`BifurcationPoint` bracket with ``k=1``.

    The signal is unambiguous in eigenvalue terms: a complex conjugate pair on
    the unit circle has coalesced on the real axis at +1 and split into a
    real reciprocal pair. This is the saddle-center / pitchfork at lambda=+1
    (codimension 1 Hamiltonian-with-symmetry bifurcation, RTR2026 Sec. p.5).

    Parameters
    ----------
    seeds :
        Ordered family members. Order matters: the bracket is between
        adjacent members in the sequence.
    primary, secondary :
        CR3BPSystem labels passed through; only ``system.mu`` is read from
        each member.
    rtol, atol :
        Integrator tolerances for the monodromy. Defaults match the rest of
        the module (1e-12 / 1e-12).

    Returns
    -------
    list of BifurcationPoint :
        Adjacent-member brackets with ``k=1`` and
        ``extras["classification_before"]`` / ``extras["classification_after"]``
        carrying the secondary-pair labels at each end. May be empty.

    Notes
    -----
    The detector is robust to the discrete signal flipping more than once: a
    family that crosses the bifurcation, then crosses back, will emit two
    brackets in canonical adjacency order. The caller resolves.
    """
    if len(seeds) < 2:
        return []
    per_member_kind: list[str] = []
    per_member_eigs: list[NDArray[np.complex128]] = []
    per_member_pair: list[tuple[complex, complex]] = []
    per_member_det: list[float] = []
    for mem in seeds:
        sysm = cr3bp.CR3BPSystem(
            mu=mem.mu,
            primary=primary,
            secondary=secondary,
            l_km=1.0,
            t_s=1.0,
        )
        mat = monodromy(sysm, mem.state0, mem.period, rtol=rtol, atol=atol)
        eigs = floquet_multipliers(mat)
        kind, lam_a, lam_b = _classify_secondary_pair(eigs)
        per_member_kind.append(kind)
        per_member_eigs.append(eigs)
        per_member_pair.append((lam_a, lam_b))
        per_member_det.append(_deflated_determinant(eigs))

    brackets: list[BifurcationPoint] = []
    transition_kinds = {"complex_unit_circle", "real_near_one"}
    for i in range(len(seeds) - 1):
        k_a = per_member_kind[i]
        k_b = per_member_kind[i + 1]

        # Detector 1: Eigenvalue classification transition
        class_transition = (k_a != k_b) and (k_a in transition_kinds) and (k_b in transition_kinds)

        # Detector 2: Scalar test function sign change (augmented Jacobian det)
        det_a = per_member_det[i]
        det_b = per_member_det[i + 1]
        det_sign_change = det_a * det_b < 0.0

        if not (class_transition or det_sign_change):
            continue

        a, b = seeds[i], seeds[i + 1]
        lam_a, _ = per_member_pair[i]
        lam_b, _ = per_member_pair[i + 1]
        pa = a.parameter if a.parameter is not None else float("nan")
        pb = b.parameter if b.parameter is not None else float("nan")

        ca = (
            0.0
            if k_a == "complex_unit_circle"
            else (1.0 if k_a == "real_near_one" else float("nan"))
        )
        cb = (
            0.0
            if k_b == "complex_unit_circle"
            else (1.0 if k_b == "real_near_one" else float("nan"))
        )

        brackets.append(
            BifurcationPoint(
                k=1,
                members=(a, b),
                eig_before=lam_a,
                eig_after=lam_b,
                dist_before=float(abs(lam_a - 1.0)),
                dist_after=float(abs(lam_b - 1.0)),
                tol=0.1,
                extras={
                    "param_before": pa,
                    "param_after": pb,
                    "classification_before": ca,
                    "classification_after": cb,
                },
            )
        )
    return brackets
