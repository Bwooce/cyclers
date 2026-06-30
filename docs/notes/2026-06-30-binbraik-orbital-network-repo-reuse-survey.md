# Reuse survey — binbraik/cislunar_orbital_network (Braik-Ross 2026 orbital-networks code)

**Date:** 2026-06-30. Repo: https://github.com/binbraik/cislunar_orbital_network — the MATLAB
implementation behind Braik & Ross 2026 "Orbital Networks in the Three-Body Problem"
(arXiv:2605.31543, already mined #249). Surveyed for reusable code/algorithms.

## License — CLEAR
**MIT** (© 2026 Abdullah Braik & Shane D. Ross). Unrestricted use/modify/redistribute, no copyleft.

## Language
100% MATLAB (R2022b+, optional Parallel Computing Toolbox). No Python/Julia/C → nothing drop-in;
reuse = adopt their DATA as goldens, or reimplement their ALGORITHM design in Python.

## Three actionable reuse items (prioritized)

### 1. (HIGH) Exact corrected ICs — `src/cr3bp_family_ic.m` → unblock #249 recovery gate
The paper's own corrected ICs (x0, ydot0, period) for all 13 families at the Earth-Moon planar
CR3BP energy, **including C11a and C21 — the two cyclers our #249 1-DOF corrector FAILED to
recover**. Critically the repo prints **C21 at CJ = 3.129389531054557 (15 sig figs)**, not the
published rounded 3.1294. THIS IS THE KNOWN ISSUE in [[feedback_published_rounded_values_are_display]]
("C21 was 1.05e-5 outside the family at literal 3.1294") — the unrounded value is exactly what that
memory said we needed. Adopt as **sourced golden ICs** (independently computed via the repo's
`correct_po_to_cj313_v2.m`, verified against Table 2 instability rates/periods).
→ Task #495.

### 2. (HIGH) Paper-validated ΔV/network goldens — `example_outputs/*.csv`
- `dc_refined_summary_mps.csv` — DC-refined ΔV + TOF for all 78 family pairs (independent of our code).
- `DVmatrix_mps.csv` — the 13×13 proxy ΔV matrix (m/s), the paper's primary result.
- `snapshot_summary.csv` — 325-row sweep over (DV_cap, Tmax) with edges_kept / LCC / centrality.
Adopt as sourced golden data to validate the #249 reachable_network scorer at the PAIR level
(bypasses the 6-node-truncation issue that blocked the full-network gate). → Task #495.

### 3. (HIGH) Two-phase corrector design — `src/traj_diffcorr.m` → fix #405/#411 Newton stall
Their inter-family transfer corrector uses **feasibility-first → optimization-second**:
1. SQP (fmincon) from warm start; 2. on constraint-recovery failure, run `lsqnonlin` on ‖r(z)‖²
ALONE (drop the ΔV cost) to reach feasibility; 3. restart SQP from the feasible point for ΔV
optimality. Plus a two-tolerance ODE strategy (loose during FD-gradient, tight on final arc) +
persistent eval cache. **Our #405/#411 cross-system Newton STALLS at |R|=0.59 rad** — almost
certainly the same "initial residual too large for single-phase Newton" failure. Reimplement in
Python: `scipy.optimize.least_squares` (phase-1 feasibility) → `scipy.optimize.minimize(method=
'SLSQP')` (phase-2). Algorithm design is transferable (not MATLAB-specific). → Task #496.

## Lower-value / not reusable
- `atlas_hits_log_from_traj.m` (vectorized voxel-segment logger) — ~4× speedup for #249 atlas if it
  ever bottlenecks; a Python port, not blocking.
- `overlap_pair.m` (theta-mirror LUT for BRS=R(FRS)) — #318/#494 tube-intersection cross-check; we
  use Poincaré sections not voxel bitmaps → complementary, not drop-in.
- All MATLAB infra (cache/parfor/graph/viz) — not reusable.
- #314 heteroclinic: their family-to-family heading-fan + SQP is a DIFFERENT problem from our
  manifold-tube Floquet seeding → no new capability.

## Defensive cross-check (bonus)
Their `cr3bp_potential.m` carries a "BUG FIX: was r^7" comment for the Uxx/Uyy/Uxy denominators.
**Our `core/cr3bp.py` already uses r^5 correctly** (verified) — we do not have this bug.

## Net
MIT-clear. Best wins: (#495) adopt the exact C21/C11a ICs + CSV ΔV goldens to unblock the #249
recovery gate and the rounded-3.1294 issue; (#496) port the two-phase feasibility-first corrector
to break the #411 cross-system stall. Both are HIGH-value, both are reimplement/adopt-data (no code
copy needed). Clone retained at scratchpad/binbraik_repo for extraction.
