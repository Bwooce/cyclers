"""Period-multiplying bifurcation detector for CR3BP families (#266 Phase 2).

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

    Returns
    -------
    NDArray (6, 6) :
        The monodromy matrix.

    Raises
    ------
    RuntimeError
        Propagated from the integrator if it fails.
    ValueError
        If ``state0`` is not a 6-vector.
    """
    s0 = np.asarray(state0, dtype=np.float64)
    if s0.shape != (6,):
        raise ValueError(f"monodromy: state0 must be a 6-vector, got shape {s0.shape}")
    arc = cr3bp.propagate(system, s0, float(period), with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    return np.asarray(arc.stm, dtype=np.float64)


def floquet_multipliers(monodromy_matrix: NDArray[np.float64]) -> NDArray[np.complex128]:
    """Return the monodromy's eigenvalues sorted by magnitude descending.

    The CR3BP monodromy has reciprocal-pair structure (Liouville / time-reversal
    symmetry): if ``lambda`` is an eigenvalue so is ``1/lambda``. One pair is
    always ``(1, 1)`` (energy/time-translation). The detector below uses these
    multipliers, sorted, without separating the trivial pair -- the
    period-multiplying check is on every multiplier individually.
    """
    eigs = np.linalg.eigvals(monodromy_matrix)
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
