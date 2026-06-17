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

**P1.1 commit:** `4b840a3`.

---

## P1.2 — Natural-parameter continuation of (3,2) C32 family in CJ

**Method.** `cyclerfinder.search.cr3bp_continuation.continue_family` walks from the C32 anchor at CJ=3.1294 upward in C by dC=1e-4, using a secant predictor in (C, x0) and the project's `correct_symmetric_fixed_jacobi` Newton corrector. Each member runs the full gauntlet (period bounds, equilibrium, Jacobi conservation, dedup, independent Radau cross-check). The radau_closure_tol is loosened to 5e-2 (vs the default 1e-3) because the (3,2) family is very unstable (|λ_max| ~ 2.5e5 at the anchor; the orbit's chaos amplifies round-trip integrator drift).

**Driver script.** `scripts/floquet_phase1_p1_2_walk.py`. Writes `data/floquet_phase1_c32_family.jsonl` with one header row + one member row per kept family member. Each member row carries: corrected (x0, ydot0, period, jacobi); the 6 Floquet multipliers from `monodromy()` and `floquet_multipliers()`; a diagnostic decomposition splitting the eigenvalues into "trivial pair" (closest to +1) vs "non-trivial" (others); the unit-circle complex pair's |λ| and |arg|; Floquet σ per Braik-Ross eq. 20; topology cross-check.

**Run.** 250 steps, walk time 633.6s. The branch ran to `max_steps` (no fold reached at the corrector layer within the budget).

**Saddle-center bifurcation FOUND between steps 123 and 124 (C ∈ [3.14170, 3.14180]).** The signature is unambiguous in the eigenvalue evolution:

  * Step 120 (C=3.14140): non-trivial complex pair at 0.9943 ± 0.1064j on the unit circle; trivial pair (1, 1)-cluster split by integrator into 0.9999 ± 0.0060j.
  * Step 123 (C=3.14170): non-trivial complex pair has migrated to 0.9989 ± 0.0473j — nearly coalescing on the real axis at +1.
  * Step 124 (C=3.14180): **all 6 eigenvalues are real.** The bifurcating pair has SPLIT on the real axis to (0.9742, 1.0264) — a new real reciprocal pair. Trivial pair sits at (0.9965, 1.0035).
  * Step 130 (C=3.14240): the bifurcated real pair has separated to (0.8778, 1.1392) — a clear new unstable direction; trivial pair is complex again at 0.9999 ± 0.0019j.

This is a saddle-center / pitchfork at λ=+1 — the classical Hamiltonian-with-symmetry bifurcation where a unit-circle eigenvalue pair coalesces at +1 then bifurcates into a real reciprocal pair, spawning a NEW family branch tangent to the marginal eigenvector.

**Member-level numerics around the bifurcation:**
  * C* ∈ (3.14170, 3.14180); a finer-step refinement bracket would tighten this, but for P1.3 the discrete signal is sufficient.
  * Period at C* ≈ 75.94 days (vs 78.61 at the anchor; ~3.4% shorter at the bifurcation).
  * x0 at C* ≈ -0.2843 (vs -0.2752 at the anchor).
  * |λ_max| at the bifurcation ≈ 3.0e4 (the "secondary saddle" direction; the primary saddle direction stays at ~3e4).

**Difference from the RTR2026 Fig. 4 picture.** RTR2026 Fig. 4 lower plot shows the bifurcation as a SADDLE-NODE FOLD in (C, x0) — the (3,2) family curves back at C* with two branches (stable + unstable orange) meeting at it. What we ACTUALLY see in the data is a SUPERCRITICAL-PITCHFORK-TYPE bifurcation where the family stays single-valued in C and continues past C* with a new pair of real eigenvalues, NOT a fold-back. The signature in eigenvalue terms is consistent — both saddle-node and pitchfork have a unit-circle pair coalescing at +1 — but the global topology differs. This is one of two possibilities:

  1. The (3,2) family Braik-Ross uses is a SLICE of a more complex 2-parameter family, and the saddle-center bifurcation in RTR2026 Fig. 4 happens in a different parameter direction (e.g. in x0 at fixed C).
  2. Our continuation parameterization moves SMOOTHLY through the bifurcation point because the marginal eigenvector is transverse to the (C, x0) continuation direction — the family keeps existing for C > C*, just with a new unstable manifold. Past C* there are two NEW families branching off perpendicular to the (C, x0) curve.

Either interpretation is consistent with the eigenvalue evidence. For P1.3-P1.4 the operational fact that matters: at C ∈ (3.14170, 3.14180) the monodromy has a non-trivial multiplier transiting +1, the bifurcation eigenvector exists, and the marginal direction can be used to perturb off the symmetric family onto the new (asymmetric or branched) family.

**P1.2 gate.** PASS. Continuation runs cleanly to 250 members covering C ∈ [3.1294, 3.1544]; the saddle-center signature is clearly visible in the data. Proceed to P1.3 (saddle-center detector).

**P1.2 commits:** `scripts/floquet_phase1_p1_2_walk.py` + `data/floquet_phase1_c32_family.jsonl` artifact (251 rows: 1 header + 250 members).

---
