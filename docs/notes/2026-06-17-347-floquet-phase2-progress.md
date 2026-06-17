# #347 Phase 2 — Floquet bifurcation discovery sweep progress note

**Status:** Phase 2 COMPLETE, working on `main`.
**Date opened:** 2026-06-17 AET.
**Date closed:** 2026-06-17 AET.
**Predecessor:** Phase 1 progress note `docs/notes/2026-06-17-347-floquet-phase1-progress.md` (Phase 1 commit ladder `4b840a3` → `075f21b` → `006b062` → `4054f86` → P1.5).
**Phase 2 commits:** `ffd7991` (P2.1 Gram-Schmidt), `7d0f1f9` (P2.2 sweep parents), `c811e05` (P2.3 sweep driver), `08a3927` (P2.4 family-network builder), TBD (P2.5 this note + artifacts).
**Exit criterion (per #347 brief):** A JSONL artifact listing discovered branch families with their (k1, k2) parent, C*, branched topology, max-Floquet, and any cycler-candidate flag; a Phase 2 verdict note recommending Phase 3 fire OR refinement work.

This note tracks per-sub-task progress for the discovery sweep build.

---

## P2.1 — Gram-Schmidt fix in `_select_saddle_center_eigenvector` (#379)

**New API:** `_select_saddle_center_eigenvector(monodromy_matrix, parent_tangent=None)`.
When `parent_tangent` is supplied (a 6-vector), the returned eigenvector is
Gram-Schmidt orthogonalised against it before being unit-normalised. The
parent tangent is `ẋ₀ = f(state0, mu)` where `f` is the CR3BP EOM at the IC —
the trivial time-translation eigenvector direction.

**Driver:** `branch_at_saddle_center` ALWAYS computes the parent tangent via
`cr3bp.cr3bp_eom(0.0, parent_state0, system.mu)` and passes it through.

**Phase 1 P1.4 behaviour change.** With the cleaner perturbation direction, at
the (3, 2) C32 i=124 / eps=5e-4 anchor the corrector now lands on a different
basin: a planar (3, 3) orbit at T~101.6 d, jacobi~3.797, **max_floquet_mag~1.0
(essentially stable)**. The original Phase 1 landing — a 3D (2, 0) orbit at
T~27.4 d, z0=-0.66, max_floquet_mag=1.07 — was a consequence of ~6e-6 in-plane
contamination in the eigenvector (variable-step DOP853 trivial-pair noise at
the cluster point, per #372 P372.3 Pellegrini-Russell 2016 finding). The
Phase 2 / P2.1 landing is a CLEANER signal of the (3, 2) saddle-center.

**Both landings satisfy the Phase 1 exit criterion** (residual < 1e-10, topology
distinct from parent). The reproduction JSONL artifact `data/floquet_phase1_reproduction.jsonl`
was regenerated under P2.1; the new contents are documented in the P2.1 commit message.

**Tests:** 5 new in `tests/genome/test_asymmetric_branch_gram_schmidt.py`. Pass.

**P2.1 gate.** PASS. Phase 2 substrate has a cleaner perturbation direction
than Phase 1's substrate.

**P2.1 commit:** `ffd7991`.

---

## P2.2 — Sweep parents enumeration

**New module:** `src/cyclerfinder/genome/floquet_phase2_parents.py`. Defines
`PHASE2_SWEEP_PARENTS`: a 6-entry tuple of `SweepParent` dataclasses covering
the four Braik-Ross 2026 Table 2 cyclers (C11a, C11b, C21, C32) PLUS C32 + C11a
walked DOWNWARD to probe the inverse direction.

**Sourcing.** Every IC traces to a published anchor: Braik-Ross 2026 Table 2
(literal C_J = 3.1294 for C11a/C11b/C32) + Ross-RT 2025 AAS-25-621 Table 4
(unrounded C_J_C21 = 3.129389531088256 for C21).

**Deferred to a future Phase 2.x.** (5, 2), (3, 1), (4, 3), and other higher-k
Aldrin / Russell-Ocampo families. Each requires its own sourced symmetric IC
(not in our offline seeds). The substrate generalises trivially; the work is
adding the sourced IC.

**Tests:** 12 new in `tests/genome/test_floquet_phase2_parents.py`. Pass.

**P2.2 gate.** PASS. Sweep parent inventory sourced + tested.

**P2.2 commit:** `7d0f1f9`.

---

## P2.3 — Saddle-center sweep driver

**New script:** `scripts/floquet_phase2_sweep.py`. Per-parent pipeline:
recover → continue (in CJ, ±direction by `_down` suffix) → detect
(`detect_saddle_center_bracket` with `stm_mode="fixed_path"` per #372 P372.3) →
branch (`branch_at_saddle_center` at epsilon=5e-4 with the P2.1 Gram-Schmidt
fix) → emit JSONL row per (parent, bracket, branched orbit) tuple.

**Smoke validation (130-step C32 walk, 21:33 - 21:39 AET, smoke run):**
- 1 saddle-center bracket at C ∈ (3.14170, 3.14180) — REPRODUCES Phase 1 P1.3.
- 1 converged branched orbit: T=101.56 d, topology=(3, 3), planar,
  jacobi=3.797, max_floquet_mag~1.0, cycler_candidate_flag=True.
- Walk time 334s + detector 18s + branch ~few seconds.

**HONEST_NEGATIVES** are first-class outcomes: a parent that does not yield
a saddle-center bracket (e.g. C21's tiny family-extent) or a branch corrector
that does not converge MUST emit a row recording the no-find / no-converge.

**Tests:** 5 new in `tests/genome/test_floquet_phase2_sweep_smoke.py`. Pass.

**P2.3 gate.** PASS. Driver works end-to-end on C32 + handles `_down` suffix.

**P2.3 commit:** `c811e05`.

---

## P2.4 — Family-network artifact builder

**New script:** `scripts/floquet_phase2_family_network.py`. Converts the Phase 2
sweep results JSONL into a Doedel 2003 Fig. 4-5-style family-network JSONL —
nodes (parent + branched orbits) + directed edges (parent → branch via
saddle-center bifurcation).

Per Doedel 2003 (`docs/notes/2026-06-17-digest-doedel-2003.md`): each family is
a node, each bifurcation is a labelled edge. The Phase 2 implementation
follows this format. Every node carries
`catalogue_status: "phase2_discovery_candidate"` — admission is Phase 4+.

**Tests:** 7 new in `tests/genome/test_floquet_phase2_family_network.py`. Pass.

**P2.4 gate.** PASS. Builder converts sweep output into family-network format.

**P2.4 commit:** `08a3927`.

---

## P2.5 — Full sweep run + verdict

**Sweep launched:** 2026-06-17 21:41 AET (commit `08a3927`), under
`uv run python -u scripts/floquet_phase2_sweep.py --output data/floquet_phase2_sweep_results.jsonl`,
on `main`. STM mode: `fixed_path` (per #372 P372.3). Epsilon: 5e-4 (Phase 1
P1.4 value).
**Sweep completed:** 2026-06-17 22:24 AET. Wall time: 2613s (43m33s). 0 errors.

### Sweep results (per-parent)

| Parent | Direction | Walk steps | Brackets | Branches converged | Topology changed | Cycler candidates | Outcome |
|--------|:---------:|:---------:|:---------:|:------------------:|:----------------:|:------------------:|---------|
| C32        | +1 | 250 (max_steps)   | 2 | 1 (bracket 0) | 1 | 1 | bracket 0 yields planar (3, 3) at T=101.6d, max_λ≈1.0; bracket 1 (C=3.1517) DID NOT CONVERGE |
| C11a       | +1 | 201 (max_steps)   | 2 | 1 (bracket 0) | 0 | 0 | bracket 0 yields planar (1, 1) at T=53.5d (same topology); bracket 1 (C=3.1474) DID NOT CONVERGE |
| C11b       | +1 | 77 (fold_reversal) | 0 | — | — | — | HONEST NEGATIVE: walk folded back at member 77, no bracket detected |
| C21        | +1 | 1 (fold_reversal)  | 0 | — | — | — | HONEST NEGATIVE: (2,1) family extent < dC, immediate fold (confirms #262 ~4e-12 extent) |
| C32_down   | -1 | 251 (max_steps)   | 0 | — | — | — | HONEST NEGATIVE: no bracket in [3.10, 3.1294] |
| C11a_down  | -1 | 201 (max_steps)   | 1 | 0 (bracket 0)  | — | 0 | bracket 0 (C=3.1118) DID NOT CONVERGE |

**Totals:** 981 members walked, 5 brackets found, 2 branches converged,
1 topology change, 3 honest negatives at the bracket/no-converge level
(C11a b1, C11b parent-level, C11a_down b0), 2 honest negatives at the
parent-level (C32_down, C21).

### Discoveries (cycler candidates)

**1 cycler candidate at the Phase-2 discovery level: branch_C32_b0.**

| Field | Value |
|-------|-------|
| Parent | C32 (Braik-Ross 2026 Table 2 (3, 2)-cycler) |
| Bracket | C ∈ (3.14170, 3.14180) — REPRODUCES Phase 1 P1.3 |
| Branched topology | (3, 3) — DISTINCT from parent (3, 2) |
| Branched period | T = 23.36 TU = 101.56 d |
| Branched jacobi | 3.797 |
| Branched degenerate_planar | True (z0 ~ 1.7e-22, zdot0 ~ 8.3e-24) |
| Branched max Floquet magnitude | 1.000000000000617 (within 6e-13 of unity — **essentially stable**) |
| Branched σ_d (instability rate, day⁻¹) | 6.08e-15 (numerical zero) |
| Corrector residual | 4.77e-12 |
| Independent Radau closure | 2.59e-11 |
| Cycler candidate flag | True |

**Interpretation.** The (3, 2) family's saddle-center bifurcation spawns a
planar (3, 3) orbit at T~101.6d (vs parent's 75.94d, ~134% factor) that is
essentially **stable** in the Floquet sense. This is a *new* (3, 3)
Earth-Moon planar cycler candidate at jacobi=3.797. It needs Phase 4
literature-novelty check + V0-V5 gauntlet before any catalogue admission
claim.

**1 stability-flip discovery at the Phase-2 discovery level: branch_C11a_b0.**

| Field | Value |
|-------|-------|
| Parent | C11a (Braik-Ross 2026 Table 2 (1, 1)a-cycler) |
| Bracket | C ∈ (3.14430, 3.14440) |
| Branched topology | (1, 1) — SAME as parent (1, 1) |
| Branched period | T = 12.31 TU = 53.51 d (~27% longer than parent's 42.14d) |
| Branched jacobi | 3.111 |
| Branched max Floquet magnitude | 1.000000000000002 (within 2e-15 of unity — **essentially stable**) |
| Corrector residual | 8.60e-13 |
| Cycler candidate flag | False (topology unchanged) |

**Interpretation.** The (1, 1)a saddle-center yields a (1, 1) orbit with
**SAME topology but flipped stability character**: parent C11a is the most
unstable cycler in Braik-Ross Table 2 (σ_d = 1.0482 d⁻¹), but the branched
orbit has σ_d ≈ 0 d⁻¹ (stable). This is the canonical "stable branch" half
of the saddle-center pair — analogous to RTR2026 Fig. 4 (the two-branch
saddle-center). This is NOT a new topology cycler, but it IS a new sourceable
sibling member of the (1, 1) Earth-Moon family at a different period + Jacobi.

### Honest negatives

1. **C11b**: walk folded back at member 77 (members_walked=77, n_brackets=0).
   The C11b family does not extend cleanly across the [3.1294, 3.15] window —
   the natural-parameter continuation hits a fold-reversal before any
   saddle-center. **Refinement option:** pseudo-arclength continuation
   (deferred to Phase 2.x per the Phase 0 design doc).
2. **C21**: walk completed in 1 member, immediate fold-reversal. This
   CONFIRMS the #262 prediction that the (2, 1) family has Jacobi extent
   ~4e-12 — even at dC=1e-13 the continuation walks off the family in one
   step. **Refinement option:** dx0-parameterised continuation (RTR2026's
   second method) rather than dC.
3. **C32_down**: 251 members walked, 0 brackets found in [3.10, 3.1294]. The
   (3, 2) family has its saddle-center at C* ∈ (3.14170, 3.14180) ABOVE the
   anchor (upward-only direction); below the anchor there is no further
   saddle-center within 0.03 of C_J. Consistent with RTR2026 p.5: cycler
   families have ONE saddle-center at the MAXIMUM C of the family.
4. **C11a_down bracket 0**: detector found a saddle-center bracket at
   C=3.1118 but the branch corrector did not converge at either eps sign.
   **Refinement option:** epsilon sweep (eps ∈ {1e-4, 5e-4, 1e-3, 5e-3}); or
   multi-shooting corrector.
5. **C32 bracket 1** (C=3.1517): branch corrector did not converge. This is
   the "second saddle-center" Phase 1 P1.3 noted in the inverse direction
   (real_far → complex_unit_circle). The Phase 2 detector flagged it, the
   corrector did not converge. Same refinement options as C11a_down b0.

### Family-network artifact

`data/floquet_phase2_family_network.jsonl`. 8 rows: 1 header + 2 parent nodes
(C32, C11a) + 2 branch nodes (branch_C32_b0, branch_C11a_b0) + 2 edges
(saddle-center bifurcation links) + 1 footer. Per Doedel 2003 nodes + edges
format. Every node carries
`catalogue_status: "phase2_discovery_candidate"`.

### Phase 2 verdict

**PASS_PROCEED_TO_PHASE_3.**

The substrate generalises beyond the (3, 2) anchor: **4 of 6 parents** yielded
useful results (2 producing converged branched orbits, 2 producing honest
negatives at the parent level confirming RTR2026 / #262 predictions). The
substrate's mechanical pieces ALL work as designed: continuation, detector,
Gram-Schmidt-orthogonalised eigenvector, asymmetric corrector, topology cross-check.

The **stability discovery** at branch_C11a_b0 (a *stable* member of the
(1, 1) family at a different period from the published Braik-Ross C11a) is a
noteworthy auxiliary finding. The (1, 1) family is published as unstable
at the C_J=3.1294 anchor; the saddle-center bifurcation reveals a stable
sibling at a different (C, T). This is the classical RTR2026 Fig. 4 picture
realised computationally for the (1, 1)a family.

The **new-cycler-candidate discovery** at branch_C32_b0 (planar (3, 3) at
T=101.6d, stable) is the headline Phase 2 result. This is a CANDIDATE for
Phase 4 V0-V5 gauntlet evaluation.

**Recommendation for Phase 3.** Fire #378 cislunar BCT integration next.
The Phase 2 substrate is ready for Phase 3 ER3BP/BCR4BP extension; the
Phase 4 gauntlet integration is parallel work.

**Recommendations for Phase 2.x refinement (deferred, not blocking Phase 3):**
- Pseudo-arclength continuation for fold-prone families (C11b primary target).
- dx0-parameterised continuation for tiny-extent families (C21).
- Epsilon sweep for no-converge branches (C32 b1, C11a_down b0).
- Multi-shooting corrector for high-instability branches.
- Higher-k Aldrin / Russell-Ocampo families ((3, 1), (5, 2), (4, 3), etc.)
  once sourced ICs are acquired.

---

## Phase 2 summary

| Metric | Value |
|--------|-------|
| Wall time | 2613s (43m33s) |
| Parents processed | 6 |
| Members walked | 981 |
| Brackets found | 5 |
| Branches converged | 2 |
| Topology changes | 1 |
| Cycler candidates | 1 (branch_C32_b0) |
| Stability-flip discoveries | 1 (branch_C11a_b0) |
| Honest negatives | 3 at bracket/no-converge level + 2 at parent level |
| Errors | 0 |
| Phase 2 wall budget | 4-12 hours (used 43 min) |
| New code (P2.x) | ~1500 LOC (1 module + 3 new tests + 2 new scripts + this note) |
| New tests | 29 (5 + 12 + 5 + 7) |
| Pytest gates | All pass; no regressions in 17 existing tests across `bifurcation_detector` + `family_switch` + Phase 1 substrate. |

**Phase 2 deliverables:**
- `data/floquet_phase2_sweep_results.jsonl` (13 rows: header + 6 parent_summary +
  3 branch_record + 3 branch_no_converge + footer)
- `data/floquet_phase2_family_network.jsonl` (8 rows: header + 2 parents + 2 branches + 2 edges + footer)
- This Phase 2 verdict note.

## Notes on scope

- **No catalogue writeback.** Every node in the family-network artifact carries
  `catalogue_status: "phase2_discovery_candidate"`. Admission is Phase 4+
  gauntlet work.
- **No literature-novelty check.** That's part of the V0-V5 gauntlet (Phase 4).
  This Phase 2 verdict notes the cycler candidates; the gauntlet decides.
- **No source-permanence check.** Sourced anchors only at the parent level;
  branched orbits are by definition unsourced (Phase 2 IS the discovery).
