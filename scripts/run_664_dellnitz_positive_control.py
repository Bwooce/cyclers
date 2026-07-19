"""#664 Dellnitz et al. 2005 quasi-Hilda transport positive control (`#661`
shortlist item 2).

Mandatory positive control for the new set-oriented transfer-operator search
paradigm (``search/set_oriented_transfer_operator.py`` -- the generic
GAIO-style box/transition-matrix/spectral-almost-invariant-set machinery --
built on top of ``search/quasi_hilda_positive_control.py``, the Sun-Jupiter
PCR3BP-specific glue). See both modules' own docstrings for the full method
description, the escape/leakage convention, the eigenvector-vs-graph-
partition method-lineage note, and exactly how the box-covering domain and
the R (quasi-Hilda) / Q (Mars-crossing) regions were reconstructed from first
principles (the corpus digest,
``docs/notes/2026-06-30-digest-dellnitz2005-mars-transport-gaio.md``, has no
machine-readable Fig. 1 pixel data, so this is a from-scratch, physically
self-consistent reconstruction at the SAME model/section/energy, not a
literal pixel-for-pixel reproduction).

What this script does
----------------------
1. Build the box grid over the section domain (x in [-0.885, -0.03], xdot in
   [-1.0, 1.0] -- the accessible zero-velocity-curve band closest to the Sun
   at E=-1.52, i.e. the paper's own "interior realm"), masking out boxes
   whose center is off the energy manifold.
2. Build the sparse transition matrix via short (one Poincare-map iterate)
   propagations of Monte Carlo samples per box, using this project's own
   validated ``core.cr3bp`` DOP853 propagator (no new integrator).
3. Extract almost-invariant sets from the leading eigenvectors of the
   transition matrix (three requested, matching the paper's own "three
   almost-invariant sets" finding -- via a DIFFERENT decomposition method,
   see the module docstring).
4. Report which discovered cluster(s) overlap the geometrically-constructed
   quasi-Hilda region R, and each cluster's almost-invariance ratio.
5. Compute the R -> Q (Mars-crossing) transport probability over up to 200
   iterates and compare its order of magnitude against the paper's sourced
   "~6% after 200 iterates" result.

Runtime: a real compute job (tens of thousands of short CR3BP propagations),
NOT a fast unit test -- expect several minutes, run synchronously.
"""

from __future__ import annotations

import time
from collections.abc import Callable

import numpy as np

from cyclerfinder.search.quasi_hilda_positive_control import (
    DELLNITZ_2005_ENERGY,
    a_mars_nondim,
    energy_to_jacobi_constant,
    mars_crossing_indicator,
    quasi_hilda_region_indicator,
    section_map_xxdot,
    sun_jupiter_system,
    zero_velocity_v,
)
from cyclerfinder.search.set_oriented_transfer_operator import (
    BoxGrid,
    almost_invariance_ratio,
    almost_invariant_sets_spectral,
    build_transition_matrix,
    transport_probability,
)

X_MIN, X_MAX = -0.885, -0.03
XDOT_MIN, XDOT_MAX = -1.0, 1.0
NX, NXDOT = 36, 36
N_SAMPLES_PER_BOX = 25
N_EIGVECS = 3
N_TRANSPORT_ITERS = 200
SEED = 664


def build_mask(grid: BoxGrid, mu: float, c_target: float) -> np.ndarray:
    mask = np.zeros(grid.n_boxes, dtype=bool)
    for flat in range(grid.n_boxes):
        center = grid.box_center(flat)
        x, xdot = float(center[0]), float(center[1])
        if zero_velocity_v(x, mu, c_target) - xdot * xdot >= 0.0:
            mask[flat] = True
    return mask


def region_box_indices(
    grid: BoxGrid,
    mask: np.ndarray,
    indicator: Callable[[float, float, float, float], bool | None],
    mu: float,
    c_target: float,
) -> np.ndarray:
    idx = []
    for flat in range(grid.n_boxes):
        if not mask[flat]:
            continue
        center = grid.box_center(flat)
        x, xdot = float(center[0]), float(center[1])
        result = indicator(x, xdot, mu, c_target)
        if result:
            idx.append(flat)
    return np.array(idx, dtype=np.intp)


def main() -> None:
    t_start = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] #664 Dellnitz 2005 quasi-Hilda positive control starting")

    system = sun_jupiter_system()
    mu = system.mu
    c_target = energy_to_jacobi_constant(DELLNITZ_2005_ENERGY, mu)
    print(f"  Sun-Jupiter mu = {mu:.10e}, E = {DELLNITZ_2005_ENERGY}, C_target = {c_target:.10f}")
    print(f"  a_Mars (nondim, Jupiter-SMA units) = {a_mars_nondim():.6f}")

    grid = BoxGrid(
        lower=np.array([X_MIN, XDOT_MIN]),
        upper=np.array([X_MAX, XDOT_MAX]),
        counts=(NX, NXDOT),
    )
    mask = build_mask(grid, mu, c_target)
    n_valid = int(mask.sum())
    print(
        f"  Box grid: {NX}x{NXDOT} = {grid.n_boxes} boxes, {n_valid} masked-in (on energy manifold)"
    )

    r_boxes = region_box_indices(grid, mask, quasi_hilda_region_indicator, mu, c_target)
    q_boxes = region_box_indices(grid, mask, mars_crossing_indicator, mu, c_target)
    overlap = np.intersect1d(r_boxes, q_boxes)
    print(
        f"  R (quasi-Hilda) region: {r_boxes.size} boxes; "
        f"Q (Mars-crossing) region: {q_boxes.size} boxes; "
        f"R & Q overlap: {overlap.size} boxes (must be 0 -- geometrically disjoint by construction)"
    )

    rng = np.random.default_rng(SEED)

    def map_fn(pt: np.ndarray) -> np.ndarray | None:
        return section_map_xxdot(float(pt[0]), float(pt[1]), mu, c_target)

    t0 = time.time()
    print(
        f"[{time.strftime('%H:%M:%S')}] building transition matrix: "
        f"{n_valid} boxes x {N_SAMPLES_PER_BOX} samples = "
        f"{n_valid * N_SAMPLES_PER_BOX} map evaluations..."
    )
    result = build_transition_matrix(grid, mask, map_fn, N_SAMPLES_PER_BOX, rng)
    dt = time.time() - t0
    print(
        f"[{time.strftime('%H:%M:%S')}] transition matrix built in {dt:.1f}s "
        f"({result.p.nnz} nonzeros)"
    )
    mean_escape = float(result.escape_fraction[mask].mean())
    print(f"  Mean per-box escape fraction (masked-in boxes): {mean_escape:.4f}")

    print(
        f"[{time.strftime('%H:%M:%S')}] extracting almost-invariant sets "
        f"(top {N_EIGVECS} eigenvectors)..."
    )
    decomp = almost_invariant_sets_spectral(result, n_eigvecs=N_EIGVECS)
    print("  Leading eigenvalues (real parts):", [f"{v.real:.6f}" for v in decomp.eigenvalues])
    n_clusters_show = min(6, len(decomp.cluster_sizes))
    print(f"  Cluster sizes (largest {n_clusters_show} of {len(decomp.cluster_sizes)}):")
    for rank, size in enumerate(decomp.cluster_sizes[:6]):
        cluster_id = np.where(decomp.labels == rank)[0]
        ratio = almost_invariance_ratio(result, cluster_id) if cluster_id.size else 0.0
        r_overlap = np.intersect1d(cluster_id, r_boxes).size
        q_overlap = np.intersect1d(cluster_id, q_boxes).size
        print(
            f"    cluster {rank}: {size} boxes, almost-invariance ratio={ratio:.4f}, "
            f"overlap with R={r_overlap} ({100 * r_overlap / max(1, r_boxes.size):.1f}% of R), "
            f"overlap with Q={q_overlap} ({100 * q_overlap / max(1, q_boxes.size):.1f}% of Q)"
        )

    print(
        f"[{time.strftime('%H:%M:%S')}] computing R->Q transport probability "
        f"over {N_TRANSPORT_ITERS} iterates..."
    )
    probs = transport_probability(result, r_boxes, q_boxes, N_TRANSPORT_ITERS)
    checkpoints = [1, 5, 10, 25, 50, 100, 150, 200]
    for n in checkpoints:
        if n <= N_TRANSPORT_ITERS:
            print(f"    p_R,Q({n}) = {probs[n]:.6f} ({100 * probs[n]:.4f}%)")

    total_dt = time.time() - t_start
    print(f"\n[{time.strftime('%H:%M:%S')}] total runtime: {total_dt:.1f}s")
    print("=" * 78)
    print("Sourced comparison (Dellnitz et al. 2005, Fig. 4): ~6% after 200 iterates")
    print(f"This reconstruction:  p_R,Q(200) = {100 * probs[200]:.4f}%")
    print("=" * 78)


if __name__ == "__main__":
    main()
