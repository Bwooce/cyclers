# #436 direct-e>0 ER3BP seeding + branch-switching — verdict

**Date:** 2026-06-24
**Status:** COMPLETE — capability delivered; **NO discovery**. The campaign's
`e_only_candidate` classification is a continuation artifact (demonstrated), not
evidence of novel e>0-only families. Registry-grade METHODOLOGICAL negative. No
catalogue writeback.

## What was run

`scripts/run_436_direct_er3bp.py` (parallelised across 12 cores). The last
untested ER3BP discovery mechanism: a **CR3BP-independent** blind symmetric-IC
grid `[x0, 0, 0, 0, ydot0, 0]` placed directly at the target eccentricity,
forward-converged, then **reverse-continued toward e=0** to test for a circular
limit. A family that "dies" before e=0 has no CR3BP ancestor → candidate
e>0-only family. Plus an ER3BP fixed-period branch-switcher (contingent infra).

- Grid 12×12 over x0∈(0.1,0.95), ydot0∈(−4,4), both period_f∈{2π, π}.
- Systems: Earth-Moon (0.0549), Sun-Mercury (0.206), Sun-Mars (0.093),
  Sun-Pluto (0.249). 1152 seeds total. Wall 4188 s.

Output: `data/er3bp_direct_436.jsonl` (200 converged records).

## Raw tally (BEFORE scrutiny)

| system | grid | converged | cr3bp_continuous | e_only_candidate |
|---|---|---|---|---|
| Earth-Moon | 288 | 69 | 0 | 69 |
| Sun-Mercury | 288 | 36 | 1 | 35 |
| Sun-Mars | 288 | 48 | 0 | 48 |
| Sun-Pluto | 288 | 47 | 10 | 36 |
| **total** | 1152 | 200 | 11 | **188** |

188/200 converged seeds flagged `e_only_candidate` (death_e median 0.047, max
0.23). **At face value this would be ~188 novel e>0-only families — which is the
"it closed!" danger signal, not a jackpot.** Two red flags forced scrutiny:
(1) a 94% flag rate; (2) a suspicious μ-ordering — 0/69 continuous at Earth-Moon
(largest μ) vs 10/11 of the continuous ones at Sun-Pluto (smallest μ).

## Scrutiny — the classification is a step-size artifact (DECISIVE)

Re-classifying ONE Earth-Moon candidate (`x1-yd2`, x0=0.177, ydot0=−2.545), the
identical converged orbit reverse-continues to:

| n_steps | classification | death_e |
|---|---|---|
| 30 | **cr3bp_continuous** (reaches e≈0) | — |
| 60 | **e_only_candidate** | 0.054 |
| 90 | **cr3bp_continuous** (reaches e≈0) | — |

The verdict **flips cr3bp_continuous ↔ e_only_candidate purely on the step
count.** Worse, the campaign recorded this seed as `e_only_candidate` at
n_steps=30 while a fresh n_steps=30 run returns `cr3bp_continuous` — so the
classification is also **run-to-run nondeterministic** (BLAS reduction-order
noise at the margin, the same class as the #347 / FBS flakes).

**Conclusion:** the "death" detected by the reverse continuation is a fragility
of the secant-predictor + 2-variable symmetric corrector tracking an arbitrary
blind-grid orbit across a wide e-range — NOT a physical family boundary. The
`e_only_candidate` flag does not reliably distinguish a no-CR3BP-limit family
from a CR3BP-continuous one. **The 188 flags are not evidence of any novel
e>0-only family.** The μ-ordering corroborates: tracking robustness degrades with
μ, so larger-μ Earth-Moon over-flags and tiny-μ Pluto (near-2-body, easy to
track) under-flags.

## Reading (honest negative)

Combined with #432 (continuation-from-CR3BP) and #435 (high-e Sun-planet seeds),
all three ER3BP discovery mechanisms are now exercised and **none has produced a
reliable novel e>0-only cycler:**
- #432/#435: known CR3BP families persist smoothly into e>0, no bifurcation.
- #436: the direct-e>0 + reverse-continuation discriminator is **not reliable**
  — its candidate flag is dominated by continuation step-size/BLAS artifacts.

This is the correct, integrity-preserving outcome: the scrutiny **prevented a
false-positive discovery claim** of ~188 "families" that are numerical artifacts.

The branch-switcher fired (converged) on 5 candidates, but since the candidates
themselves are artifacts these are not meaningful bifurcations; and #432/#435
found zero genuine bifurcations, so the switcher remains contingent infra with no
validated target.

## Disposition

- **No catalogue writeback. No discovery claim.** (The 188 candidates are
  artifacts; matches the closure discipline — "it closed!" is the danger signal.)
- Method-versioned negative registered in `data/empty_regions.jsonl`.
- **The genuine fix (follow-on):** replace the secant reverse-continuation +
  2-var symmetric corrector with a **fold-aware pseudo-arclength family
  continuation** (the robust tracker #434 already uses for 3D families) before
  the no-CR3BP-limit discriminator can be trusted. Only then is a direct-e>0
  campaign capable of a defensible discovery. The seed grid + branch-switcher are
  built and ready to feed a robust tracker.
