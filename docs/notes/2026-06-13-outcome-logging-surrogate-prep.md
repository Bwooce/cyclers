# Passive solver-outcome logging — ANN surrogate prep (#210)

Date: 2026-06-13

## What this is

A passive, opt-in, near-zero-cost logger that appends one `(inputs -> outcome)`
tuple to a training log every time a solver decides an outcome. It is the
**costless on-ramp** to a future ANN cycler-search surrogate (#210): we begin
accumulating the `(genome -> outcome)` corpus now, as a byproduct of normal runs,
so the data exists when (if) a surrogate is trained.

Module: `src/cyclerfinder/search/outcome_log.py` — `log_outcome(...)`.

## HARD BOUNDARY — training capture, never validation

This log is **only ever written, never read back to validate a cycler.** An ANN /
surrogate output can never be a golden; the V0–V5 validation gauntlet still governs
every catalogue claim. Nothing in the validation path reads this file. This boundary
is restated in the module docstring and is the whole reason the warning lives next to
the code.

## How to enable

Logging is a NO-OP unless the environment variable is set:

```sh
export CYCLERFINDER_OUTCOME_LOG=/path/outside/the/repo/outcomes.jsonl
```

- Unset (the default for all normal runs and CI): zero cost — one env lookup, then
  return. No file is created; behaviour is byte-identical to before.
- Set: each wired solver appends one JSON line at the point its outcome is decided.
- The path is written **outside the repo by convention**; it is env-var-driven and
  never committed. (If a default path under the repo is ever introduced, add it to
  `.gitignore`.)

## Record envelope (every record)

```json
{
  "schema_version": 1,
  "counter": 0,                // process-local monotonic, deterministic-friendly
  "wall_time": 1.7e9,          // convenience side field, NOT the ordering key
  "solver": "<call-site name>",
  "inputs":  { ... },          // the genome / decision inputs (JSON-coerced)
  "outcome": { ... },          // the decided outcome scalars
  "meta":    { ... } | null    // optional caller context
}
```

numpy scalars/arrays are coerced to float/int/bool/list via `to_jsonable`. The call
is wrapped in try/except: a logging failure is swallowed (warned to stdlib logging),
never raised into the caller — capturing data must never break a solve.

## Schema per wired solver (the `inputs`/`outcome` payloads)

| solver string | inputs | outcome |
|---|---|---|
| `cr3bp.correct_periodic` | `state0_guess[6]`, `period_guess`, `mu` | `converged`, `residual`, `period`, `jacobi` |
| `cr3bp.correct_symmetric_fixed_jacobi` | `x0_guess`, `jacobi`, `period_guess`, `ydot0_sign`, `half_crossings`, `mu` | `converged`, `crossing_residual`, `x0`, `ydot0`, `period`, `n_iter` |
| `optimise_cell_idealized` | `cell_id`, `bodies`, `sequence`, `period_k`, `vinf_cap`, `n_starts`, `seed`, `use_de` | `converged`, `constraints_satisfied`, `closure_residual_kms`, `max_vinf_kms`, `taxi_cost_kms` |
| `mbh.hop` | `x_seed`, `rng_seed`, `hop` | `x`, `objective`, `feasible` |

`meta` carries `primary`/`secondary` (CR3BP) or `perturbation` (MBH).

### Wiring (one guarded `log_outcome` call each — behaviour-preserving)

- `src/cyclerfinder/search/cr3bp_periodic.py` — `correct_periodic` and
  `correct_symmetric_fixed_jacobi` (the latter also captures **every per-step
  continuation result**, because `cr3bp_continuation.continue_family` calls that
  corrector once per step — so the continuation campaign is covered without touching
  `cr3bp_continuation.py`).
- `src/cyclerfinder/search/optimize.py` — `optimise_cell_idealized` (the heliocentric
  per-cell search).
- `src/cyclerfinder/search/mbh.py` — `mbh()` (one record per basin-hop trial: the
  seed solve and each perturbed hop).

## Label sources for a future surrogate

Two complementary `(genome -> outcome)` label sources now exist:

1. **The JPL oracle** — `src/cyclerfinder/verify/jpl_periodic_orbits.py` provides
   sourced periodic-orbit ground truth (the positive-label anchor set). *(Not touched
   by this work; named here only as the authoritative label source.)*
2. **Continuation / corrector logs** — the JSONL captured here: dense
   converged/not-converged labels across `(C, x0)` and across heliocentric cells,
   including the **negatives** (non-converged, off-basin), which a surrogate needs as
   much as the positives.

## Data-floor context (Ozaki 2022)

Ozaki et al. 2022 (arXiv:2111.11858) — the reference ANN trajectory surrogate
architecture — needs on the order of **7e6 samples** to train. Passive capture is the
only realistic way to approach that floor at zero marginal cost; a dedicated sampling
campaign for 7e6 points would be prohibitive. This is exactly why the logger is opt-in
and free when off: turn it on during normal large sweeps and the corpus accrues.

## #210 status

- **Costless prep: DELIVERED.** Logging is implemented, wired at the four high-value
  solver outcome points, tested, and on (opt-in via the env var) with zero behaviour
  change to any solver.
- **Still deferred:** the KKT pseudo-body amplifier (the active sample-amplification
  half of #210) — that remains future work; nothing here trains or runs a surrogate.
