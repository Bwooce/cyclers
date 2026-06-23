# Unified Global MGA-DSM Precursor Search Engine — Design

**Date:** 2026-06-23
**Issue lineage:** follows #307 (multi-rev precursor re-run, honest search-method-limited negative). This build attacks the *search-method* limit #307 exposed: a circular-body Tisserand seed + a **local** Nelder-Mead optimiser could not reach a real-ephemeris ballistic basin into the classic Earth-Mars cyclers, even with multi-rev Lambert.

## Goal

Replace the precursor matcher's local-optimiser path with a **global** search engine that hunts ballistic→powered precursors into a target cycler, using **eccentric-body-aware seeds**, and ranks survivors by total ΔV / `dv_band`. Optimised to MAXIMISE the chance of surfacing a real precursor (hunt-the-positive); any survivor is routed through the full V0-V5 gauntlet before any catalogue claim.

## Decisions (locked in brainstorming, 2026-06-23)

1. **Direction:** eccentric-body Tisserand seeds + global search (a targeted upgrade to the existing real-ephemeris matcher — distinct from the ER3BP genome #293, which makes the rotating-frame *primaries* eccentric).
2. **Goal:** hunt the positive (coverage + restart-richness over provability).
3. **Precursor type:** ANY `dv_band`, ranked by cost — ballistic if it exists, powered (DSM) otherwise. Wires in #307 Task-2's `evaluate_dsm_leg`.
4. **Architecture:** full integrated rebuild — a unified global MGA-DSM engine replacing the matcher's optimiser path, with #302/#307 reproducibility protected (the closure primitives are reused unchanged; `find_cycler_precursors` delegates, signature preserved).

## Non-goals (YAGNI)

- NOT an eccentric *primaries* 3-body model (that is ER3BP, #293 — already built, separate question).
- NOT a global guarantee of non-existence (that was the *other* #307 fork, explicitly not chosen).
- NOT a new closure model — `close_epoch_locked` already runs on real DE440; we reuse it verbatim.
- NOT a 3D/broken-plane search (the v∞ sphere is #388 Phase 2 / #414; out of scope here — coplanar real-eph as today).

## Architecture & module boundaries

New module `src/cyclerfinder/search/global_precursor_engine.py` owns the search. Reused, untouched primitives:

- `cyclerfinder.core.lambert.lambert` (multi-rev branch enumeration).
- `cyclerfinder.search.mga_dsm_placement.evaluate_dsm_leg` (#307 Task 2, Vasile-Conway DSM transcription).
- `cyclerfinder.genome.epoch_aware_genome.close_epoch_locked` (real DE440 closure, the residual oracle).
- `cyclerfinder.core.ephemeris.Ephemeris` (DE440 backend).

`cyclerfinder.search.precursor_matcher.find_cycler_precursors` is rewired to **delegate** to the engine. Its public signature is preserved (callers `run_302_*`, `run_307_*`, `tests/search/test_precursor_matcher.py` keep working). The old local-Nelder-Mead `optimise_chain_tofs` is retained as a private fallback (used by the delegated path only if the engine is explicitly disabled), so no caller breaks.

**#302/#307 reproducibility guard:** a regression test pins `close_epoch_locked`'s output (closure residual + flyby continuity + per-encounter V∞) on a fixed, hard-coded candidate. Because the engine reuses that primitive unchanged, this test proves the closure physics did not move; only the *search over* it changed.

## Components

### (a) Eccentric-body Tisserand-Poincaré seeder
`global_precursor_engine.py::eccentric_tp_seeds(...)`. Extends the linkability + resonance-TOF logic in `tisserand_mga_window.py` to use each body's **actual heliocentric state at the encounter epoch** (DE440 `r`, `v`) rather than `PLANETS[body].sma_au` (mean, circular). The Tisserand-Poincaré graph drawn at the body's actual heliocentric radius `r_p` (true-anomaly dependent) shifts/widens the linkable contour set. Sourced from **Campagnola-Russell 2009 Endgame Part B** (multibody T-P graph) and **Strange-Russell-Buffington 2007** ("Mapping the V∞ globe"), both in corpus (`docs/notes/2026-06-05-endgame-tisserand-mining.md`). Output: a set of candidate chains (body sequence, V∞ bins, seed TOFs, seed epochs) — reuses the `MGAChainCandidate` type — that become the differential_evolution init population.

### (b) Unified parameter vector + global optimiser
`global_precursor_engine.py::search_precursors(...)`. Decision vector for an N-leg chain:
`x = [epoch_offset_days, tof_1..tof_N, (eta_i, dvx_i, dvy_i, dvz_i for i=1..N)]`.
Optimiser: `scipy.optimize.differential_evolution` over a bounded box (the load-bearing change vs #307's local Nelder-Mead), with `init` set to a population that includes the (a) eccentric seeds, `polish=True` for a final local refine. Bounds: epoch ±`epoch_half_width_days`; TOFs inside `tof_box_days_per_leg`; eta ∈ [0,1]; DSM components ∈ [−`dsm_max_kms`, +`dsm_max_kms`]. Hunt-the-positive knobs exposed: `popsize`, `maxiter`, `n_restarts`, `seed` list (multi-start).

**Discrete-vs-continuous split:** the body sequence (e.g. E-V-E) is *discrete* and is fixed per DE run; the continuous vector `x` is optimised *within* a fixed sequence. The seeder (a) enumerates the discrete candidate sequences (with their seed `x` vectors); the engine runs one differential_evolution per candidate sequence, then aggregates and ranks (d) across all sequences. Per-sequence DE runs are independent and parallelisable (joblib, per the `#321` substrate).

### (c) Per-leg DSM (Vasile-Conway)
Each leg evaluation: if `eta_i` and `|dv_i|` are effectively zero → ballistic leg (multi-rev Lambert to the next body); else → `evaluate_dsm_leg` (arc1 ballistic propagate with the DSM applied at `eta_i`, arc2 Lambert to arrival). The DSM Δv is thus a first-class decision variable, not a post-hoc repair — the optimiser trades it against closure/continuity directly.

### (d) Cost-ranked output
Objective = `closure_residual + w_cont * flyby_continuity + w_dsm * total_dsm_dv`, with `w_dsm` set so a ballistic solution always scores below an otherwise-equal powered one (ballistic preferred). Survivors ranked by **total ΔV ascending**, each tagged with `dv_band` (strictly_ballistic / essentially_ballistic / low_maintenance / powered_dsm), closure residual, flyby continuity, per-leg DSM tuple, and literature check. Emitted as JSONL via an extended `precursor_match_to_jsonl_record` (adds `per_leg_dsm_kms` + `total_dsm_dv_kms` + `dv_band`).

## Data flow

```
cycler row (seed V∞ + first body)
  → eccentric_tp_seeds()                    # (a) DE440-aware seeds
  → differential_evolution(init=seeds)      # (b) global search
        per trial vector x:
          build N legs (multi-rev Lambert + optional DSM)   # (c)
          close_epoch_locked(real DE440)                    # reused oracle
          cost = closure + w_cont*continuity + w_dsm*ΣdsmΔv  # (d)
  → ranked survivors (by total ΔV, dv_band-tagged)
  → check_literature()
  → JSONL
  → (any survivor) V0-V5 gauntlet + ML flagger + #256 FP guard
```

## Error handling

- Kepler/Lambert non-convergence on a trial vector → that evaluation returns a large **finite** penalty (`_COST_FLOOR`, e.g. 1e6), never `inf`, so differential_evolution's mutation/selection stays well-behaved (mirrors `evaluate_dsm_leg`'s try/except → `pos_res=inf` handling, but capped finite for the optimiser).
- Infeasible geometry (negative/zero TOF, hyperbolic arc where elliptic required, eta out of [0,1] after clamp) → same penalty floor.
- A chain whose seeder produces no feasible legs at any V∞ is dropped from the candidate set with a logged count (no silent truncation — per the negative-results discipline).

## Testing (all goldens sourced; `feedback_golden_tests_sourced_only`)

1. **#302/#307 reproducibility guard** — `close_epoch_locked` output (closure, continuity, per-encounter V∞) byte-stable on a pinned candidate.
2. **Eccentric-T-P seeder golden** — linkability/contour value at a sourced Campagnola-Russell 2009 case (value traced to the paper, not self-computed); plus a reduces-to-circular check (at e→0 the seeder matches the mean-`a` enumerator).
3. **Engine recovers a known ballistic solution** — on a synthetic real-eph E→M direct transfer with a known ballistic answer, `search_precursors` converges to ≈0 total DSM and the correct TOF.
4. **DSM reduces-to-ballistic** — η/Δv→0 leg evaluation matches the no-DSM multi-rev Lambert closure.
5. **dv_band ranking** — unit tests on the ranker (ordering + band thresholds at 1/10/300 m/s boundaries).
6. **Determinism** — fixed `seed` ⇒ reproducible ranked output (differential_evolution is seeded).

## Deliverable run

Re-run Aldrin (`aldrin-classic-em-k1-outbound`) + S1L1 (`s1l1-2syn-em-cpom`) through the engine; rank survivors by `dv_band`. Route any sub-gate ballistic OR low-ΔV powered survivor through the V0-V5 gauntlet + ML flagger + #256 FP guard before any catalogue claim. Honest-negative-or-promote; register the method-versioned outcome in `data/empty_regions.jsonl` (this engine capability-subsumes the #307 multi-rev negative). Compare residual distribution vs the #307 baseline in a results note.

## Decomposition (for writing-plans)

~7 tasks: (1) eccentric-T-P seeder + golden; (2) parameter-vector leg builder (multi-rev + DSM) + reduces-to-ballistic test; (3) cost function + dv_band ranker + tests; (4) `search_precursors` differential_evolution driver + known-ballistic-recovery test; (5) `find_cycler_precursors` delegation + #302/#307 reproducibility guard; (6) extended JSONL record; (7) Aldrin/S1L1 deliverable run script + results note + registry entry.

Estimated ~1000–1400 LOC across the new engine + matcher rewire + tests.
