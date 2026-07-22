"""#685 Sun-Earth PCR3BP set-oriented (GAIO) Earth-Mars transport search
(`#661` shortlist item 2 lineage; own-system application of the #664 pipeline).

Applies the #664 validated set-oriented transfer-operator machinery
(``search/set_oriented_transfer_operator.py`` -- generic GAIO primitives,
reused verbatim) to THIS project's own flagship Earth-Mars transport domain
for the first time, via the Sun-Earth PCR3BP glue in
``search/sun_earth_mars_transport.py``. See that module's docstring and #685's
``data/OUTSTANDING.md`` bullet for the full model/section/energy/region
rationale and the honest framing of what a positive vs. negative result means.

This is a genuinely novel SEARCH (not a reproduction), so it runs an explicit
CALIBRATION pass before trusting any spectral result, exactly as #664 did for
its Dellnitz positive control:

  C1. Energy self-consistency: every masked-in box centre's constructed section
      state reproduces the target Jacobi constant under core.cr3bp.jacobi_constant.
  C2. Map conservation: the exterior Poincare return map conserves the Jacobi
      constant on hand-picked trajectories (a real dynamical cross-check).
  C3. Region partition sanity: R (Earth-neighbourhood) and Q (Mars-reaching)
      are disjoint and both non-empty on the actual box grid.
  C4. Escape/leaky-mass accounting: a CORRECTED map wrapper (see below) so
      off-manifold sample ICs are redrawn while genuine dynamical escapes are
      counted as escape; the mean per-box escape fraction is reported and must
      be sane (neither ~0 -- an artificially sealed domain -- nor ~1 -- a
      domain that empties in one iterate, making the operator meaningless).

Only then does it build the transfer operator (two independent seeds),
extract almost-invariant sets from the leading eigenvectors, and measure the
R->Q transport probability.

Escape-accounting note (correction over the raw section map)
------------------------------------------------------------
``sun_earth_mars_transport.section_map_xxdot`` returns ``None`` for BOTH an
off-manifold initial condition AND a genuine dynamical escape (no qualifying
x>0 return). ``build_transition_matrix``'s retry logic is designed to REDRAW a
None sample (treating it as an invalid IC) up to ``max_point_retries`` times
before counting escape -- correct for off-manifold ICs, but it would silently
convert genuine escapes into retained mass, biasing the operator toward
FALSE-POSITIVE (over-trapping) structure. The wrapper below separates the two:
an off-manifold IC returns ``None`` (redrawn, as intended), while a genuine
escape returns a far-field SENTINEL point that lands outside the grid and is
therefore counted as escape WITHOUT a redraw. This keeps escape accounting
honest and biases conservatively (against, not toward, finding structure).

Runtime: a real compute job (tens of thousands of short CR3BP propagations per
seed), NOT a fast unit test -- expect a few minutes; run synchronously.
"""

from __future__ import annotations

import time

import numpy as np

from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.search.set_oriented_transfer_operator import (
    BoxGrid,
    almost_invariance_ratio,
    almost_invariant_sets_spectral,
    build_transition_matrix,
    transport_probability,
)
from cyclerfinder.search.sun_earth_mars_transport import (
    SUN_EARTH_MARS_TRANSFER_C,
    a_mars_nondim,
    earth_neighborhood_region_indicator,
    mars_reaching_indicator,
    osculating_elements_at_section,
    poincare_first_return_exterior,
    section_state6,
    sun_earth_system,
    zero_velocity_v,
)

X_MIN, X_MAX = 1.03, 1.60
XDOT_MIN, XDOT_MAX = -0.85, 0.85
NX, NXDOT = 36, 36
N_SAMPLES_PER_BOX = 25
N_EIGVECS = 3
N_TRANSPORT_ITERS = 200
SEEDS = (685, 686)
# Far outside the grid (x = 10 >> X_MAX): lands in no box -> counted as escape.
ESCAPE_SENTINEL = np.array([10.0, 0.0], dtype=np.float64)


def build_mask(grid: BoxGrid, mu: float, c_target: float) -> np.ndarray:
    mask = np.zeros(grid.n_boxes, dtype=bool)
    for flat in range(grid.n_boxes):
        center = grid.box_center(flat)
        x, xdot = float(center[0]), float(center[1])
        if zero_velocity_v(x, mu, c_target) - xdot * xdot >= 0.0:
            mask[flat] = True
    return mask


def region_box_indices(grid, mask, indicator, mu, c_target) -> np.ndarray:  # type: ignore[no-untyped-def]
    idx = []
    for flat in range(grid.n_boxes):
        if not mask[flat]:
            continue
        center = grid.box_center(flat)
        x, xdot = float(center[0]), float(center[1])
        if indicator(x, xdot, mu, c_target):
            idx.append(flat)
    return np.array(idx, dtype=np.intp)


def make_map_fn(mu: float, c_target: float):  # type: ignore[no-untyped-def]
    """Escape-correct map wrapper (see module docstring): None for off-manifold
    ICs (redrawn), far-field sentinel for genuine escapes (counted as escape)."""

    def map_fn(pt: np.ndarray) -> np.ndarray | None:
        state0 = section_state6(float(pt[0]), float(pt[1]), mu, c_target)
        if state0 is None:
            return None  # off-manifold IC -> redraw
        result = poincare_first_return_exterior(state0, mu)
        if result is None:
            return ESCAPE_SENTINEL  # genuine escape -> counted as escape, no redraw
        return np.array([result[0], result[3]], dtype=np.float64)

    return map_fn


def calibration(grid: BoxGrid, mask: np.ndarray, mu: float, c_target: float) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] --- CALIBRATION ---")
    # C1: energy self-consistency of section states over masked-in box centres.
    max_c_err = 0.0
    for flat in np.where(mask)[0]:
        c = grid.box_center(flat)
        s = section_state6(float(c[0]), float(c[1]), mu, c_target)
        if s is not None:
            max_c_err = max(max_c_err, abs(jacobi_constant(s, mu) - c_target))
    print(f"  C1 section-state Jacobi self-consistency: max |C-C_target| = {max_c_err:.2e}")

    # C2: map conserves Jacobi on hand-picked trajectories.
    worst = 0.0
    for x, xdot in [(1.10, 0.05), (1.30, 0.20), (1.45, -0.10), (1.55, 0.30)]:
        s = section_state6(x, xdot, mu, c_target)
        if s is None:
            continue
        c0 = jacobi_constant(s, mu)
        r = poincare_first_return_exterior(s, mu)
        if r is not None:
            worst = max(worst, abs(jacobi_constant(r, mu) - c0))
    print(f"  C2 exterior map Jacobi conservation: worst |dC| = {worst:.2e}")

    # C3: hand-check an apsis crossing (xdot=0 at x is an osculating apside).
    el = osculating_elements_at_section(1.20, 0.0, mu, c_target)
    assert el is not None
    print(
        f"  C3 osculating hand-check @ (1.20, 0): a={el.a:.4f} e={el.e:.4f} "
        f"r_p={el.r_p:.4f} r_a={el.r_a:.4f} (r_a should ~= 1.20, the crossing apsis)"
    )


def run_seed(
    grid: BoxGrid,
    mask: np.ndarray,
    mu: float,
    c_target: float,
    r_boxes: np.ndarray,
    q_boxes: np.ndarray,
    seed: int,
) -> None:
    n_valid = int(mask.sum())
    rng = np.random.default_rng(seed)
    map_fn = make_map_fn(mu, c_target)
    t0 = time.time()
    print(
        f"[{time.strftime('%H:%M:%S')}] seed {seed}: building transition matrix "
        f"({n_valid} boxes x {N_SAMPLES_PER_BOX} samples = {n_valid * N_SAMPLES_PER_BOX} evals)..."
    )
    result = build_transition_matrix(grid, mask, map_fn, N_SAMPLES_PER_BOX, rng)
    dt = time.time() - t0
    mean_escape = float(result.escape_fraction[mask].mean())
    print(
        f"[{time.strftime('%H:%M:%S')}] seed {seed}: built in {dt:.1f}s, "
        f"nnz={result.p.nnz}, mean per-box escape fraction = {mean_escape:.4f} "
        f"(C4: sane if not ~0 and not ~1)"
    )

    decomp = almost_invariant_sets_spectral(result, n_eigvecs=N_EIGVECS)
    print(
        f"  seed {seed}: leading eigenvalues (real):",
        [f"{v.real:.6f}" for v in decomp.eigenvalues],
    )
    n_show = min(6, len(decomp.cluster_sizes))
    print(f"  seed {seed}: cluster sizes (largest {n_show} of {len(decomp.cluster_sizes)}):")
    for rank, size in enumerate(decomp.cluster_sizes[:6]):
        cluster = np.where(decomp.labels == rank)[0]
        ratio = almost_invariance_ratio(result, cluster) if cluster.size else 0.0
        r_ov = np.intersect1d(cluster, r_boxes).size
        q_ov = np.intersect1d(cluster, q_boxes).size
        print(
            f"    cluster {rank}: {size} boxes, almost-invariance ratio={ratio:.4f}, "
            f"overlap R={r_ov} ({100 * r_ov / max(1, r_boxes.size):.1f}% of R), "
            f"overlap Q={q_ov} ({100 * q_ov / max(1, q_boxes.size):.1f}% of Q)"
        )

    # Does R itself trap mass? (Its own almost-invariance ratio.)
    r_ratio = almost_invariance_ratio(result, r_boxes)
    q_ratio = almost_invariance_ratio(result, q_boxes)
    print(f"  seed {seed}: R self almost-invariance ratio = {r_ratio:.4f}; Q = {q_ratio:.4f}")

    probs = transport_probability(result, r_boxes, q_boxes, N_TRANSPORT_ITERS)
    for n in [1, 5, 10, 25, 50, 100, 150, 200]:
        if n <= N_TRANSPORT_ITERS:
            print(f"    p_R,Q({n}) = {probs[n]:.6f} ({100 * probs[n]:.4f}%)")


def main() -> None:
    t_start = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] #685 Sun-Earth GAIO Earth-Mars transport search starting")
    system = sun_earth_system()
    mu = system.mu
    c_target = SUN_EARTH_MARS_TRANSFER_C
    print(f"  Sun-Earth mu = {mu:.10e}, C_target = {c_target:.10f} (Earth-Mars Hohmann Tisserand)")
    print(f"  a_Mars (nondim, Earth-SMA units) = {a_mars_nondim():.6f}")

    grid = BoxGrid(
        lower=np.array([X_MIN, XDOT_MIN]),
        upper=np.array([X_MAX, XDOT_MAX]),
        counts=(NX, NXDOT),
    )
    mask = build_mask(grid, mu, c_target)
    print(f"  Box grid: {NX}x{NXDOT} = {grid.n_boxes} boxes, {int(mask.sum())} masked-in")

    r_boxes = region_box_indices(grid, mask, earth_neighborhood_region_indicator, mu, c_target)
    q_boxes = region_box_indices(grid, mask, mars_reaching_indicator, mu, c_target)
    overlap = np.intersect1d(r_boxes, q_boxes)
    print(
        f"  R (Earth-neighbourhood): {r_boxes.size} boxes; "
        f"Q (Mars-reaching): {q_boxes.size} boxes; overlap: {overlap.size} (must be 0)"
    )
    assert overlap.size == 0, "R and Q must be geometrically disjoint"

    calibration(grid, mask, mu, c_target)

    for seed in SEEDS:
        run_seed(grid, mask, mu, c_target, r_boxes, q_boxes, seed)

    print(f"\n[{time.strftime('%H:%M:%S')}] total runtime: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
