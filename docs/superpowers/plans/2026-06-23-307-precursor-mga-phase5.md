# #307 — #289 Phase 5: precursor-MGA closure upgrades (implementation plan)

> **For agentic workers:** execute task-by-task with TDD. Steps use `- [ ]`.

**Goal:** Upgrade the epoch-locked precursor-MGA closure driver so it closes
Aldrin/S1L1 precursor insertions to publication grade (≤1 km/s closure, ≤0.1 km/s
flyby continuity), unlocking a literature-fresh precursor-MGA discovery pass.

**Architecture:** Three additive upgrades to `close_epoch_locked` in
`src/cyclerfinder/genome/epoch_aware_genome.py` (the #289 Phase-1 closure framework),
each behind a parameter so existing callers are byte-unchanged at the defaults:
multi-rev Lambert branch selection, automated DSM placement (Vasile-Conway 2006 §3.2),
and cycler-cadence terminal-phase targeting. Then re-run the #302 precursor matcher.

**Tech stack:** existing `core/lambert.py` (`max_revs>0` already supported),
`scipy.optimize.differential_evolution` (DSM outer loop), `verify/literature_check`.

**Source anchors:** Phase-5 recommendation in
`docs/notes/2026-06-16-302-289-phase4-precursor-matcher.md` §"Phase 5 path";
Vasile-Conway 2006 MGA-DSM transcription (§3.2).

---

### Task 1: Multi-rev Lambert branch selection in `close_epoch_locked`

**Files:** Modify `src/cyclerfinder/genome/epoch_aware_genome.py:372` (close_epoch_locked);
Test `tests/genome/test_epoch_aware_multirev.py`.

Design decision (the non-trivial part): each leg's Lambert branch sets v1 (→ V∞ at
encounter k) AND v2 (→ V∞ at k+1), and intermediate encounters average inbound/outbound,
so branch choice is NOT per-leg-separable. Use **bounded best-combination enumeration**:
for `max_revs=R`, each leg has ≤ `1+2R` candidate solutions; enumerate the Cartesian
product (legs are few: 2–4 → ≤ 3^4=81 combos at R=1), close each combination, return the
one with the lowest max-V∞-residual. The single-rev solution is always in the candidate
set, so the result is **never worse** than `max_revs=0` (the clean invariant test).

- [ ] **Step 1: failing test — `max_revs` param exists + never-worse invariant.**
  Build a 3-body EpochLockedTrajectory (analytic-circular ephemeris, per docstring),
  close at `max_revs=0` and `max_revs=1`; assert `closure(max_revs=1).max_vinf_residual
  <= closure(max_revs=0).max_vinf_residual + 1e-9`. Run: FAIL (no `max_revs` kwarg).
- [ ] **Step 2: implement** — add `max_revs: int = 0`; refactor the per-leg loop (lines
  500–541) to collect `lambert(..., max_revs=max_revs)` per non-DSM leg, enumerate
  combinations, evaluate the existing residual machinery per combination, keep the best.
  DSM legs keep their single-rev reference seed (DSM placement is Task 2).
- [ ] **Step 3: test passes.** Also assert `max_revs=0` is byte-identical to today
  (regression-pin the existing `tests/genome/` closure goldens).
- [ ] **Step 4: ruff + mypy + commit.**

### Task 2: Automated DSM placement (Vasile-Conway 2006 §3.2)

**Files:** new `src/cyclerfinder/search/mga_dsm_placement.py`; wire into
close_epoch_locked via an opt-in `optimize_dsms: bool`; Test
`tests/search/test_mga_dsm_placement.py`.

- [ ] Step 1: failing test — a leg with a known continuity deficit closes to ≤ target
  after DSM optimisation (free DSM `fraction_along_leg` ∈ (0,1) + `delta_v_kms` vector
  as continuous decision vars; `differential_evolution` minimising closure residual +
  λ·|ΔV|).
- [ ] Step 2: implement the transcription + outer loop (reuse the existing `DSMSpec`
  two-arc executor already in close_epoch_locked).
- [ ] Step 3: test passes; ΔV is the reported actuator cost (golden discipline: the
  optimiser's ΔV is `computed`, never asserted as a sourced value).
- [ ] Step 4: ruff + mypy + commit.

### Task 3: Cycler-cadence terminal-phase targeting

**Files:** extend `find_cycler_precursors` (`search/precursor_matcher.py:312`) to pass a
target terminal-Earth phase window from the catalogue cycler row; Test
`tests/search/test_precursor_phase_targeting.py`.

- [ ] Step 1: failing test — `epoch_alignment_score > 0` when a target phase window is
  supplied (today it is hardcoded 0 because no window is passed).
- [ ] Step 2: derive the window from the cycler row's period + encounter cadence; score
  the precursor's terminal-Earth arrival against it.
- [ ] Step 3: test passes. Step 4: commit.

### Task 4: Re-run the Aldrin + S1L1 precursor scans through the upgraded substrate

- [x] Run the #302 matcher with Tasks 1–3 enabled; record closure quality + flyby
  continuity + literature-check verdict in a results note.
  → `docs/superpowers/plans/2026-06-23-307-task4-multirev-precursor-verdict.md`.
- [x] Expected (per the Phase-4 note): closure ≤1 km/s, continuity ≤0.1 km/s, and the
  lit-check still flags 0 novel (the structural test) — OR, if a literature-fresh
  precursor survives, route it through the V0–V5 gauntlet + ML flagger + #256 FP guard
  before any catalogue claim (closure discipline: independent cross-check mandatory).
  → OUTCOME: honest negative. `max_revs=2` shifts the validated distribution
  (median continuity −2 km/s) but 0/394 cross the 0.10 km/s gate on either target;
  min continuity stays 4.15 (Aldrin) / 3.95 (S1L1) km/s; 0 literature-fresh. The
  expectation that the upgraded substrate would close the gap is falsified — the
  coplanar/circular-body model wall holds. Both rows stay V0.

## Self-review
- Spec coverage: Tasks 1–3 = the three Phase-5 deficiencies in the #302 note; Task 4 =
  the re-run + honest-negative-or-promote outcome. ✓
- Golden discipline: every optimiser/closure ΔV is `computed`, never asserted sourced;
  any promotion goes through the full gauntlet. ✓
- Backward-compat: all three upgrades are opt-in params; `max_revs=0` + no DSM-opt +
  no phase window = today's behaviour (regression-pinned). ✓
