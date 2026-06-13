# Multi-arc closure harness — design (#248)

**Goal:** a clean, reliable harness that answers — for a multi-arc E-E-M-M /
E-E-E-M-M catalogue row — *does it converge below the 0.1 km/s gate on either
optimizer lane, robustly?* Replaces the throwaway `scratch_*` scripts whose
inconsistent metrics produced the false "0.163 km/s" signal. Output of a converged
row is a **SILVER candidate → gauntlet**, never a catalogue writeback.

**Why a harness, not more scratch:** three one-shot agents hung on this (slow
solver compute they backgrounded + polled), and inline probing hit a metric
inconsistency, an identical-across-seeds bug, and an unguarded `jplephem
OutOfRangeError` crash. The failures were all *infrastructure*, not the physics —
so the fix is a disciplined harness with the failure modes designed out.

## The four things that must be right (each was wrong in scratch)

### 1. ONE canonical metric
Everything reports `dsm_chain_correct(...).max_residual_kms` and `.converged`
against `tol_kms=0.1`. No ad-hoc "res", no `dv_total` masquerading as residual.
The coarse-vs-refine disagreement (0.163 vs 2.06) was a metric mismatch; banning
ad-hoc metrics kills that class of bug.

### 2. Epoch-range-SAFE evaluation
Wrap every corrector call so a `jplephem OutOfRangeError` (or any ephemeris-range
error — DE440 covers 1549–2650) returns a large PENALTY residual instead of
crashing the run. The optimizer/grid must be free to explore wild epochs without
aborting. (This crashed the inline probe.)

### 3. Correct resonant-leg seeding (discrete returns, not a ±% scale)
The transit leg is seeded via `self_seeding.joint_epoch_tof_close` (structured
epoch — the working part). Each RESONANT leg's ToF is seeded by ENUMERATING the
discrete small-N resonant returns of that leg's body: descriptor resonant period ×
{1,2,3,…N} (plus the descriptor value itself). A ±15% scale around one value (what
the inline probe did) never reaches the right return and gives identical junk.

### 4. Multi-start, not one charged seed
The #244 root cause is single-charged-seed basin selection. Drive the
structured-epoch × resonant-return seed grid through the existing MBH wrapper
(`search/mbh.py`) — basin hopping over the discrete seed set — so we sample basins
systematically. Lambert lane is the convergence reference (#245 needs lambert-lane
convergence parity); fbs-analytic optional as the cross-check.

## Architecture
- **New module** `src/cyclerfinder/search/multiarc_closure.py`:
  - `safe_chain_residual(x, seq, eph, bounds, gradient) -> (max_residual_kms, converged, result|None)` — the epoch-safe canonical-metric wrapper (point 2 + 1).
  - `resonant_return_seeds(row) -> Iterator[seed_vector]` — structured transit epoch (joint_epoch_tof_close) × discrete resonant-return grid per resonant leg (point 3).
  - `close_multiarc_row(row, eph, *, n_starts, gradient="lambert", budget) -> ClosureReport` — MBH multi-start over the seeds, returns best canonical residual + converged + the seed that achieved it (point 4).
- **Reuse:** `dsm_descriptor_seed.seed_dsm_chain_from_descriptor`, `self_seeding.joint_epoch_tof_close`, `dsm_leg.dsm_chain_correct`, `search/mbh.py`, `core/ephemeris`.
- **Tests** `tests/search/test_multiarc_closure.py`:
  - epoch-safe wrapper returns a penalty (not a crash) on an out-of-range epoch;
  - resonant-return enumerator yields the descriptor value + the integer multiples;
  - a fast smoke (single start, bounded) runs end-to-end on one row.
- **Driver** `scripts/multiarc_closure_run.py` — per-row run + a compact canonical table.

## Build sequence (bite-sized, commit each)
1. `safe_chain_residual` wrapper + test (penalty-not-crash). COMMIT.
2. `resonant_return_seeds` enumerator + test. COMMIT.
3. `close_multiarc_row` (MBH multi-start, canonical metric) + smoke test. COMMIT.
4. Run on **`mcconaghy-2006-em-k2`** (closest, E-E-M-M): does best canonical
   `max_residual_kms` cross 0.1, robustly (re-run reproduces)? Record real numbers.
5. If yes → run the 6 russell descriptor rows; re-run `fbs_optimizer_adoption_parity.py`
   on a converged row so the **#245** convergence-parity decision becomes meaningful.
6. If no → a clean NEGATIVE: report the residual floor + dominant leg per row; this
   feeds the genome track (the topology may not close in this genome).

## Execution & the daemon link
The harness is written to run BOTH ways: inline-bounded (a few rows, capped starts)
for the #248 answer now, and as a **Track-C discovery-daemon worker** (resumable,
many starts, checkpointed) later — same `close_multiarc_row` entry point. This is
deliberately the first concrete piece of the Track-C daemon.

## Discipline
- NO catalogue writeback — a converged optimization is a SILVER candidate, not a
  validated cycler (V0–V5 gauntlet governs). Same-model goldens only.
- Not a one-shot agent (3 hangs); build it inline/deliberately, commit per step.
- Honest negative is a valid outcome (report the residual floor, don't force 0.1).

## Success criteria
At least one real multi-arc row's canonical `max_residual_kms < 0.1` on the lambert
lane, **reproducibly** (not a single fluke point), with the epoch-safe wrapper and
discrete-resonant-return multi-start — OR a clean, per-row negative with residual
floors that tells us the genome, not the seeding, is the wall.
