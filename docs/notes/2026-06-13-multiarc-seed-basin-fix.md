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

## Discipline
No catalogue writeback — a converged optimization is not a validated cycler (the
V0–V5 gauntlet still governs). The incomplete production scaffold in
`dsm_descriptor_seed.py` was reverted (unused imports / not wired); the finding is
preserved here for the follow-up.
