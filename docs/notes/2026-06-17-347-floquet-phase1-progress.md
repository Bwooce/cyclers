# #347 Phase 1 — Floquet bifurcation framework progress note

**Status:** in progress, working on `main`.
**Date opened:** 2026-06-17 AET.
**Predecessor:** Phase 0 design doc `docs/notes/2026-06-17-347-floquet-bifurcation-framework-phase0-design.md` (commit `b02073d`).
**Exit criterion (verbatim):** One JSONL row showing a branched orbit with topology distinct from parent (k1, k2) = (3, 2), corrector residual < 1e-10, ≤ 5 days wall-clock.

This note tracks per-sub-task progress for the substrate build. Each P1.x sub-task appends a section below as its commit lands.

---

## P1.1 — Reproduce-gate the (3,2) Earth-Moon symmetric anchor

**Sourced anchor.** Braik-Ross 2026 Table 2 row C32: P=78.613 days, σd=0.1583 day⁻¹, CJ=3.1294. The seed IC (`x0=-0.2752115`, `ydot0_sign=-1`, `half_crossings=6`) lives in `src/cyclerfinder/search/reachable_representatives.py` at `_CYCLER_SEEDS["C32"]`. This was set up at #262 from the Ross & Roberts-Tsoukkas 2025 AAS-25-621 family seed region. The full recovery pipeline `recover_all_cyclers_braik_ross` already produces a `Representative` whose period matches Braik-Ross to 0.0005%.

**Reproduction-gate test.** New test `tests/genome/test_floquet_phase1_anchor.py` consumes the existing recovery + adds three orthogonal cross-checks:

  1. **Period gate (vs sourced)**: recovered period in days must match Braik-Ross's 78.613 d within 1% (the design-doc threshold).
  2. **Topology gate (independent classifier)**: `winding_topology` on the converged orbit yields exactly `(k1, k2) = (3, 2)` with prograde windings.
  3. **Floquet σ gate (independent diagnostic vs sourced)**: per Braik-Ross eq. (20), σd = ln(|λ_max nontrivial|) / T (day⁻¹). Recovered σd must match 0.1583 within 5% (loose tolerance to absorb integrator + the nondim/day conversion to 6 sig figs).

**Anchor-stage numerics (this run).** All three gates pass:

  * Period: 78.612611 d vs sourced 78.613 d. Match to 0.0005%.
  * Jacobi: 3.1294000000 (literal, by construction of the fixed-jacobi corrector).
  * Topology: (k1=3, k2=2), w1=3.0000, w2=2.0000, prograde=True.
  * Floquet σ (TU⁻¹): 0.6882 vs sourced 0.6886 (0.06%); σd: 0.1583 vs sourced 0.1583 (exact to 4 sig fig).
  * Monodromy eigenvalues at the anchor (sorted by |λ|):
    * λ_max = 2.531366e+05 (strongly unstable saddle direction; reciprocal pair λ_min = 3.474447e-06)
    * Two real eigenvalues near +1: 1.0147, 0.9855 — these are the trivial pair (energy/time-translation) split by integrator precision. They are *not* sitting exactly at the saddle-center bifurcation point yet (a saddle-center would require a SECOND eigenvalue arriving at +1, i.e. one of the unit-circle pair).
    * A complex-conjugate pair on the unit circle: 0.6448 ± 0.7643j (|λ|=1.0000). These are the Neimark-Sacker / center-direction eigenvalues whose Floquet exponent will drive σd.

**Observations driving P1.2 design.**

  * The (3,2) Earth-Moon anchor at CJ=3.1294 is NOT itself a saddle-center bifurcation point. RTR2026 Fig. 4 (lower panel, C - C_cr on the y-axis) shows the saddle-center sits at the family's C-maximum and the published representative is below it. The continuation in P1.2 must walk UP in C (toward C_cr).
  * The two non-trivial unit-circle eigenvalues (the 0.6448 ± 0.7643j complex pair) are the candidates for the saddle-center coalescence: at the bifurcation point this complex pair collides on the real axis at +1, then splits into two real eigenvalues (one stable branch, one unstable). The detector in P1.3 must look for this complex pair coalescing on the real axis at +1.
  * The dominant unstable direction (λ_max ≈ 2.5e5) is the σ-driver — this is the (3,2) unstable manifold structure already documented in Braik-Ross; not the bifurcation direction.

**P1.1 gate.** PASS. The anchor is reproduced to all three gates. Proceed to P1.2 (CJ-direction continuation along the (3,2) family).

**P1.1 commit:** (filled at commit time)

---
