# #432 ER3BP discovery campaign — verdict

**Date:** 2026-06-24
**Status:** COMPLETE — capability delivered + golden-validated; structural NEGATIVE
(no novel *e>0-only* cycler at the probed seeds/systems). Conditional on the
seed/system coverage below. No catalogue writeback.

## What was run

First discovery use of the #293 ER3BP genome (Phase 4). The campaign
(`scripts/run_432_er3bp_discovery.py`) continues rotating-frame CR3BP cycler
families from e=0 into a target eccentricity, computing the full-period Floquet
monodromy at each step and classifying each seed **survives / dies / bifurcates**
(a *bifurcation* = an elliptic↔hyperbolic stability transition along the
continuation — the candidate birth of an e>0-only family). Survivors get an
offline KNOWN_CORPUS literature check.

- **Phase A** (Earth-Moon, real e=0.0549, n_steps=60): 3 seeds — the sourced
  Broucke-1969 7P family + 2 Koblick-2023 NRHOs (Np=1 Gateway 9:2, Np=2 butterfly).
- **Phase B** (high-e probes): Earth-Moon Broucke family pushed to synthetic
  e=0.10 and e=0.15. Sun-Mercury (e=0.206) and Sun-Mars (e=0.093) were
  **seed-limited and skipped** — the Broucke IC is Earth-Moon-μ-specific and no
  CR3BP seed IC exists at those mass ratios; **no ICs were fabricated**.

Outputs: `data/er3bp_discovery_{phaseA,phaseB}.jsonl`.

## Result

| phase | seeds | survives | dies | bifurcates | literature-fresh |
|---|---|---|---|---|---|
| A (EM e=0.0549) | 3 | 3 | 0 | 0 | 0 |
| B (EM e=0.10, 0.15) | 2 | 2 | 0 | 0 | 0 |

**Every seeded CR3BP rotating-frame family continues smoothly into the ER3BP
past the real Earth-Moon eccentricity (and on to e=0.15) staying in its
stability regime (all 61 steps "unstable"/hyperbolic) with NO bifurcation
(`e_star=None`).** No new family branches off; all match the published corpus.
So at the probed seeds/systems there is **no novel e>0-only cycler**.

## Reading (honest scope of the negative)

This is a clean, valid discovery-campaign outcome — *and* a delivered, golden-
validated capability (the ER3BP discovery pipeline: Floquet monitor + Floquet
golden vs Broucke stability, continuation driver with survives/dies/bifurcates
classification, seed registry, literature adjudication). But the negative is
**conditional**, and the conditions matter:

1. **Seed-limited to Earth-Moon.** The systems where the ER3BP departs MOST from
   the CR3BP — high-e Sun-Mercury (0.206), Sun-Mars (0.093), Sun-Pluto (0.249) —
   were not probed because we have no CR3BP seed ICs at those μ. The most
   promising regime for e>0-only structure was therefore *not reached*. (This is
   the scoping note's prediction: ER3BP cycler-class existence is "thin/
   speculative", and the EM system is the well-studied, least-likely-novel one.)
2. **Small seed set** (3 distinct CR3BP families). A broader EM seed set
   (Lyapunov / DRO / low-order resonant families) would widen coverage.
3. **Method = bifurcation-along-continuation.** A genuine e>0-only family that
   has *no CR3BP limit at all* would be found by **direct e>0 seeding** (guess
   ICs at the target e, confirm no e→0 limit) — the alternative mechanism
   deliberately out of scope here. And a flagged bifurcation would need
   **branch-switching** to actually trace the new family. Neither was built.

So the precise claim: *the (continuation + Floquet-transition) method finds no
novel e>0-only cycler branching off the Broucke/Koblick Earth-Moon families up
to e=0.15.* It does NOT establish that ER3BP-only cyclers don't exist — the
high-e Sun-planet regime and the direct-seeding method are untested.

## Two correctness traps caught before the run (would have invalidated it)

1. **`period_f` half-vs-full-period contract** (caught by the Task-3 implementer):
   the corrector treats `period_f` as the half-period integration span; the seed
   stores the full 2π. Passing it straight made even the e=0 seed fail to converge
   → spurious *all-dies*. The driver now halves it; the Broucke seed converges to
   5e-15.
2. **Bifurcation detector** (caught via smoke test): the per-orbit "eigenvalue on
   the unit circle" flag is TRUE for every stable Hamiltonian orbit, so it fired
   `bifurcates@e=0` for everything → false-positive discoveries. Re-keyed to a
   stability *transition* across steps; the Broucke seed now correctly classifies
   `survives`.

## Disposition

- No catalogue writeback (no ER3BP V0–V5 gauntlet; matches #430/#307 discipline).
- Method-versioned negative registered in `data/empty_regions.jsonl`.
- Follow-ons (the genuine remaining frontier; logged as tasks): (a) generate
  CR3BP Lyapunov/DRO/resonant seed ICs at high-e Sun-planet μ (Mercury/Mars/Pluto)
  so Phase B can actually probe the high-departure regime; (b) direct-e>0 seeding
  for families with no CR3BP limit; (c) branch-switching to trace any flagged
  bifurcation. The pipeline is ready to consume all three.
