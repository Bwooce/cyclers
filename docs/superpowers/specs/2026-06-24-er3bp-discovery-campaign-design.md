# ER3BP Discovery Campaign (#293 Phase 4 / #432) — Design

**Date:** 2026-06-24
**Lineage:** First *discovery* use of the #293 ER3BP genome (Phase 1: `core/er3bp.py`,
`genome/er3bp_periodic.py`, `genome/er3bp_continuation.py` — built + golden-tested,
never run for search). Motivated by the #307/#430 model-wall findings pointing at
eccentricity, and by the validation-ceiling finding that novel-cycler discovery is
method-gated, not iteration-gated. The ER3BP is the one capable frontier genome with
zero discovery campaign run (`0` refs in `scripts/`).

## Goal

Continue rotating-frame cycler families from `e=0` into target eccentricity in the
Elliptic Restricted 3-Body Problem, monitoring Floquet stability to (a) flag where
known cyclers **survive vs. die** at a system's real `e`, and (b) flag **bifurcations**
(eigenvalue / unit-circle crossings) as candidate births of *e>0-only* families — the
genuinely open, literature-unanswered question: *do ER3BP cycler families exist that do
not continue back to CR3BP as e→0?* Report-only; no catalogue writeback.

## Decisions (locked in brainstorming, 2026-06-24)

1. **Seeds & systems:** BOTH, sequential. Phase A = catalogue CR3BP-class cyclers
   continued into their native system's real `e` (golden-anchored where data exists).
   Phase B = standard Lyapunov / DRO / resonant families across high-`e` systems
   (Sun-Mercury 0.206, Sun-Mars 0.093, Earth-Moon 0.0549).
2. **Discovery mechanism:** survival + Floquet-flag. Continue each seed to target `e`
   while computing the full-period (`f = 2πn`) monodromy via the existing
   `er3bp_stm_eom`; tag stability and flag unit-circle crossings (bifurcations) and
   survival/death. NO branch-switching in this campaign.
3. **Disposition:** report-only. JSONL traces + verdict note + method-versioned
   `empty_regions` entry. NO catalogue writeback (no ER3BP V0–V5 gauntlet exists —
   matches the #430/#307 discipline). Novel survivors / promising bifurcations spawn
   follow-on tasks (branch-switching; an ER3BP gauntlet), neither built here.
4. **Purity:** ER3BP only. No low-thrust / QP / multi-rev angle (YAGNI).

## Confirmed feasibility constraint (shapes Phase A)

The 12 CR3BP-class catalogue rows (Arenstorf, Genova-Aldrin 3-petal, Wittal, the four
Ross-Tsoukkas RT cyclers, the three Braik-Ross C-cyclers, Russell-Strange Saturnian)
store **no rotating-frame initial conditions** in `catalogue.yaml` (verified:
`ic-ish: []` for all). So Phase-A seeds cannot be read from the catalogue directly;
each must be **reconstructed** from (a) ICs already encoded in repo code constants
(e.g. tulip/NRHO tables, the #255 binary-star figure-read ICs, any Floquet-discovered
branch ICs), or (b) source-paper table/figure ICs. **The guaranteed-available seed
floor is the standard Lyapunov / DRO / resonant families the CR3BP genome generates
from scratch** (Earth-Moon, where Broucke-1969 / Peng-Xu-2015 golden data exists).
Phase A is therefore *best-effort on recoverable catalogue-cycler ICs*; the campaign's
floor does not depend on it.

## Architecture & module boundaries

New module `src/cyclerfinder/search/er3bp_discovery.py` owns the campaign. Reused,
unchanged: `core/er3bp.py` (`ER3BPSystem`, `er3bp_eom`, `er3bp_stm_eom`,
`propagate_er3bp`), `genome/er3bp_continuation.py` (`continue_er3bp_family_in_e`),
`genome/er3bp_periodic.py` (`ER3BPPeriodicOrbit`, `correct_er3bp_periodic`), the CR3BP
genome for generating standard seed families, `search/literature_check.check_literature`,
and the #256 ML false-positive flagger. The Floquet eigenvalue extraction reuses the
#347 monodromy/Floquet conventions.

## Components

### (a) Seed registry — `er3bp_seeds.py` (or a section of the campaign module)
`Er3bpSeed` dataclass: `label`, `system` (primary, secondary, μ), `state0`
(rotating-frame IC at e=0), `period_f` (true-anomaly period, multiple of 2π),
`is_half_period_residual`, `target_e`, `source` (provenance string). Two providers:
- `standard_family_seeds()` — generate Lyapunov / DRO / low-order resonant ICs from the
  CR3BP genome for the target systems (the guaranteed floor; includes the Broucke
  Earth-Moon golden family as a self-check seed).
- `catalogue_cr3bp_seeds()` — best-effort reconstruct ICs for the CR3BP-class catalogue
  rows from recoverable code constants / sourced paper values; SKIP (with a logged
  count, no silent drop) any row whose IC cannot be recovered.

### (b) Floquet monitor — `er3bp_floquet.py`
`er3bp_monodromy(orbit, system) -> NDArray` (full-period STM via `er3bp_stm_eom` over
`f ∈ [0, period_f]`) and `floquet_classify(monodromy) -> (eigenvalues, stability_tag,
on_unit_circle: bool)`. Stability tags mirror the #347 conventions
(stable / unstable / marginal). Sourced golden: a Broucke family's known stability.

### (c) Per-seed continuation driver — in `er3bp_discovery.py`
`continue_and_monitor(seed, n_steps) -> Er3bpContinuationTrace`. Steps `e=0→target_e`
via `continue_er3bp_family_in_e`; at each converged step records `(e, residual,
eigenvalues, stability_tag)`. Classifies the seed outcome:
- **survives** — closes at `target_e` with bounded residual;
- **dies** — continuation fails / folds before `target_e` (records `e_max_reached`);
- **bifurcates** — a unit-circle crossing / stability flip occurs at some `e* < target_e`
  (records `e_star`, the crossing eigenvalue) → candidate e>0-only family birth.

### (d) Survivor adjudication
On survivors and bifurcation-flagged seeds: run `check_literature` (structural
fingerprint) and the #256 ML FP flagger. Tag literature verdict honestly
(not-found / published / inconclusive).

### (e) Campaign runner — `scripts/run_432_er3bp_discovery.py`
Phase A then Phase B. Incremental timestamped progress logging (per
`feedback_incremental_progress_reports`); detached-run friendly (per
`feedback_long_runs_acceptable`). Emits `data/er3bp_discovery_{phaseA,phaseB}.jsonl`
(one record per seed: trace summary + outcome + flags + literature verdict).

## Data flow

```
seed registry (standard families + recoverable catalogue cyclers)
  → for each seed: continue_er3bp_family_in_e (e=0 → target_e)
        per step: er3bp_monodromy → floquet_classify → record
     → classify outcome {survives | dies | bifurcates(e*)}
  → check_literature + #256 ML flag on survivors / bifurcation-flagged
  → JSONL trace + dv-free outcome record
  → verdict note + method-versioned empty_regions entry
  → (any novel survivor / promising bifurcation) → follow-on task (branch-switch + gauntlet)
```

## Error handling

- Continuation non-convergence at a step → record the step as the death point
  (`e_max_reached`), classify **dies**, continue to the next seed (never abort the
  campaign on one seed).
- Monodromy integration failure → record `stability_tag="unknown"` for that step;
  does not crash the seed's trace.
- A catalogue seed whose IC cannot be recovered → SKIP with a logged count (no silent
  truncation — per the negative-results discipline).

## Testing (all goldens sourced)

1. **Floquet golden** — `er3bp_monodromy` + `floquet_classify` on a Broucke-1969
   Earth-Moon family reproduce the published stability character (reuses the existing
   ER3BP golden fixture).
2. **Survives-classification** — a seed known to continue (Broucke Orbit 1 → 59, the
   existing continuation golden) classifies as **survives** with finite per-step
   residuals.
3. **Dies-classification** — a seed driven to an infeasible target `e` classifies as
   **dies** with a recorded `e_max_reached` (no crash).
4. **Seed-registry floor** — `standard_family_seeds()` returns ≥1 usable seed per target
   system; `catalogue_cr3bp_seeds()` returns a (possibly empty) list and logs the
   skip count — the campaign floor never depends on catalogue IC recoverability.
5. **Reduces-to-CR3BP** — at `e=0` a continued orbit matches its CR3BP seed (guards the
   genome reuse; the ER3BP EOM already has this golden).
6. **Determinism** — a fixed seed-set + `n_steps` produces a reproducible JSONL ordering.

## Deliverable run

Detached Phase-A then Phase-B run via `scripts/run_432_er3bp_discovery.py`. Harvest →
verdict note comparing outcomes (survives / dies / bifurcates counts per system) →
method-versioned `empty_regions` entry (first ER3BP discovery capability). Honest
outcome: most likely a structural map (which known families survive to real `e`, where
they bifurcate) with bifurcation flags as the candidate novel-family leads; a clean
"no e>0-only family flagged" is a registry-grade negative. Any genuinely-novel flagged
family → follow-on branch-switching + an ER3BP gauntlet (new tasks; not built here).

## Decomposition (for writing-plans)

~7 tasks: (1) `Er3bpSeed` model + `standard_family_seeds()` + floor test;
(2) `catalogue_cr3bp_seeds()` best-effort IC recovery + skip-logging test;
(3) `er3bp_monodromy` + `floquet_classify` + Broucke stability golden;
(4) `continue_and_monitor` driver + survives/dies/bifurcates classification tests;
(5) survivor adjudication (literature + #256) wiring; (6) campaign runner script +
smoke; (7) deliverable run → verdict → registry.

Estimated ~900–1200 LOC across the new modules + tests.
