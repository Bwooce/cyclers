# Multi-arc seed/basin investigation (#248) — partial finding, NOT yet converged

**Status: PARTIAL. Structured-epoch seeding is the right direction (residual 40 →
~1.0 km/s) but full convergence (≤0.1 km/s) was NOT achieved; the last M→M
resonant DSM leg remains the binding blocker.** #245 (FBS default flip) stays
HELD. This note preserves the finding; the production wiring was not completed.

## The problem (from #244)
`close_row_dsm` seeds the whole multi-arc chain at `t0 = 0` (J2000) with the
coplanar short-way transit ToF plus a generic slack for the resonant legs. On the
real `dsm_chain_correct` lane NEITHER the FBS-analytic nor the Lambert+FD gradient
lane converges from that seed on the multi-arc E-E-M-M / E-E-E-M-M rows — the
single, badly-placed charged seed is outside the basin of any feasible solution, so
gradient quality is irrelevant (the #244 conclusion).

## What was tried (#248) and what it bought
Experiments live in untracked `scripts/scratch_*.py` (structured-seed, structured
seed-through-MBH, resonant-ToF sweep, baseline diagnostic) — reproducibility aids
for the follow-up, not production code.

1. **Structured-epoch seed.** Instead of `t0=0`, sweep the transit-leg epoch ±~1
   synodic period (via `self_seeding.joint_epoch_tof_close`, centred ~2027-01-01 —
   the era the #181 ToF-fix closer used; the seed is insensitive to the exact
   centre). **Result: `mcconaghy-2006-em-k2` residual 40 → ~1.0 km/s; several other
   rows to ~6–14 km/s.** Major progress — confirms the blocker is seed placement,
   not the gradient or genome.
2. **Structured seed → MBH multi-start.** Driving the structured seed through the
   existing MBH wrapper did not close the remaining ~1 km/s on the hard rows.
3. **Resonant-leg ToF sweep.** The ~1.0 km/s residual floor is **dominated by the
   LAST (M→M) resonant DSM leg** — the structured-epoch seed fixes the transit
   legs but the resonant leg still gets generic slack. Hypothesis (untested to
   closure): seed the resonant-leg ToF from the row's descriptor resonant period
   rather than generic slack.

## Precise remaining blocker / next step
The transit-leg seed is solved; **the resonant (M→M) leg ToF seeding is the last
mile.** A focused follow-up should: seed each resonant leg's ToF from its
descriptor resonant period (not generic slack), combine with the structured-epoch
transit seed, and re-run — then re-run `scripts/fbs_optimizer_adoption_parity.py`
on whatever row first converges so the #245 convergence-parity comparison finally
becomes meaningful. Until a row converges on EITHER lane, #245 stays HELD.

## Update 2 (2026-06-13, resonant-ToF attempt — IMPROVED but NOT converged)
A second pass added resonant-leg ToF seeding (seed each same-body/resonant leg's ToF
from small-N resonant returns of its body + the descriptor value, with the
structured-epoch transit seed fixed). Scratch harness (untracked, for the resumer):
`scripts/scratch_reso_mcconaghy.py` (coarse) + `scripts/scratch_reso_refine.py`
(fine) + the run-1 `scratch_structured_*`.

- Coarse sweep *claimed* `mcconaghy-2006-em-k2` → res 0.163 km/s at combo
  [E-reso 1096d, transit 189d, M-reso 1718d] — just over the 0.1 gate.
- **That sub-gate hit did NOT reproduce.** Running the fine refine inline (main loop,
  verified — not trusting the scratch comment) the best was **res 2.0643 km/s,
  conv=False**, combo [1144,189,1684], dv_total 3.22. The two scratch scripts also
  appear to use INCONSISTENT residual metrics (coarse "res" vs refine "res" +
  dv_total), so the 0.163 was likely a metric mismatch / non-robust point — a classic
  "it closed!" false signal (see [[feedback_orbit_closure_discipline]]).

**Honest status: PARTIAL.** Structured-epoch + resonant-ToF seeding takes the
multi-arc rows from 40 km/s (non-converging) down to the ~1–2 km/s range — the
seed/basin IS the lever and this is real progress — but **no row cleanly crosses the
0.1 km/s gate**, and the one sub-gate claim was not reproducible.

## Update 3 (2026-06-14, inline probe — INCONCLUSIVE, the scratch numbers are unreliable)
Took it inline to settle the metric. The canonical gate metric is
`dsm_chain_correct(...).max_residual_kms` with `tol_kms=0.1` (and `.converged`).
A quick canonical probe on `mcconaghy-2006-em-k2` (lambert lane, structured seed
with the *base* descriptor resonant ToF = 534 d) returned **58 km/s, non-converged**,
and IDENTICAL across a ±15% resonant-ToF sweep — an identical-across-seeds result is
a bug signature (the resonant seed wasn't affecting the outcome). This also flatly
contradicts the run-1 note's "structured seed → 1.0 km/s," so AT LEAST ONE of those
numbers is wrong/metric-confused. Re-testing with the coarse's actual per-leg
resonant-return values (≈1096 d / 1718 d) **crashed** with `jplephem OutOfRangeError`
— the larger resonant ToFs drove the corrector's epoch search outside DE440's range
(the evaluator needs epoch-range-safe guarding).

**Conclusion: sub-gate multi-arc convergence is NOT demonstrated, and quick inline
probing is unreliable** (inconsistent metrics, an identical-result bug, an
unguarded epoch-range crash). #245 stays HELD. This is a genuine, careful
engineering problem — not solvable by chat-turn thrashing or one-shot agents.

## What a proper #248 solve needs (the real unit of work)
A CLEAN multi-arc closure harness (a real module + test, not scratch):
1. ONE canonical metric throughout (`max_residual_kms` / `converged`), no ad-hoc "res".
2. Epoch-range-SAFE evaluation (catch/penalise `OutOfRangeError` instead of crashing).
3. Correct resonant-return ENUMERATION per leg (discrete small-N returns of the leg's
   body, not a ±% scale of one value) for the seed.
4. Structured-epoch transit seed (done) + MBH multi-start over the resonant-return
   grid, driven through `dsm_chain_correct`.
5. Verify on `mcconaghy-2006-em-k2` first (closest row), then the 6 russell rows.
This belongs inline (deliberate) or in the Track-C discovery daemon (its natural host
for long solver compute) — NOT a one-shot agent (3 hangs) and NOT quick probes.

## Precise next step (NOT a one-shot agent)
#248 has now hung THREE one-shot agents — it is slow iterative solver compute that
makes agents background-and-poll until killed. Do it **inline (main loop, bounded,
ONE consistent residual metric)** or fold it into the **Track-C discovery daemon**.
First fix the metric inconsistency, then a single bounded grid on `mcconaghy-2006-em-k2`
(the closest row), lambert lane only, to settle whether it crosses 0.1. Until a row
converges on EITHER lane with a consistent metric, **#245 stays HELD**.

## Update 4 (2026-06-14, the clean harness WORKS — 0.1044 km/s, at the gate edge)
The `search/multiarc_closure.py` harness (epoch-safe canonical metric + discrete
resonant-return seeds + multi-start) was built and the coordinator ran the bounded
campaign: `scripts/multiarc_closure_run.py --rows mcconaghy-2006-em-k2 --n-starts 25`
→ **best canonical max_residual_kms = 0.1044, converged=False** (25 seeds, ~126 s,
lambert lane). This is the TRUSTWORTHY number (the 58 / 2.06 / 0.163 scratch values
were buggy/metric-confused). It sits a hair over the 0.1 km/s gate — the multi-arc
row is essentially AT convergence. A refine pass (fbs-analytic lane, 40 starts,
finer resonant grid) is running to settle whether it crosses 0.1. Either way the
seed/basin lever + clean harness is the answer; #245 unblocks the moment a row
crosses (or we accept 0.1044 ≈ gate as effective convergence and re-run the parity).

## Discipline
No catalogue writeback — a converged optimization is not a validated cycler (the
V0–V5 gauntlet still governs). The incomplete production scaffold in
`dsm_descriptor_seed.py` was reverted (unused imports / not wired); findings are
preserved here for the follow-up.
