# #433 Jupiter Galilean repeated-moon quasi-cycler sweep — verdict

**Date:** 2026-06-24
**Status:** COMPLETE — registry-grade NEGATIVE for a *novel* admissible quasi-cycler.
No catalogue writeback.

## What was run

`scripts/scan_433_jupiter_galilean.py` — the proven #320/#344 epoch-aware
repeated-moon "moontour" sweep (the machinery that produced the Umbriel
catalogue row), configured for the **Jupiter Galilean + Amalthea** system at the
two levers that paid off on Saturn (#344): a finer **48×48** phase grid and a
wider per-leg **n_rev ∈ {0..5}**.

- Moons: Io, Europa, Ganymede, Callisto, Amalthea.
- Length-3 repeated-moon cycles (`X-Y-X`), 20 cycle templates.
- 700 cells evaluated, **558 resonance-screened** (the Galilean Laplace
  resonance 1:2:4 floods the grid with phase-locked periodics; the
  rotation-number irrationality / resonance screen removed them as designed).
- Wall: ~1005 s.

Output: `data/scan_433_jupiter_galilean.jsonl`.

## Result — ranked best closure per cycle

| # | sequence | n_rev | closure (km/s) | dv_band | physical gate | lit anchors | lit-fresh |
|---|---|---|---|---|---|---|---|
| 1 | Callisto-Ganymede-Callisto | (1,1) | **0.0083** | silver_gate | **PASS** | 6 | **No** |
| 2 | Io-Amalthea-Io | (1,1) | 0.0340 | silver_gate | FAIL | 0 | Yes |
| 3 | Ganymede-Io-Ganymede | (1,1) | 0.0373 | silver_gate | FAIL | 5 | No |
| 4 | Io-Europa-Io | (1,1) | 0.0675 | sub_0.1 | PASS | 5 | No |
| 5 | Ganymede-Callisto-Ganymede | (1,1) | 0.0891 | sub_0.1 | PASS | 6 | No |
| … | (12 more cycles, all ≥0.09 km/s) | | | | | | |

(1 SILVER-band cell total; 50 near-miss-band cells; full ranking in the run log.)

## Reading (honest scope of the negative)

The discovery gate is the conjunction **sub-0.05 km/s AND physical-max-bend PASS
AND literature-fresh**. No cycle satisfies all three:

1. The single strongest closure — **Callisto-Ganymede-Callisto at 8.3 m/s**
   (silver band, physical PASS) — is **not literature-fresh** (6 KNOWN_CORPUS
   anchors). It is a re-derivation of the well-studied outer-Galilean resonant
   tour region (Ganymede/Callisto are the canonical Galileo/JUICE/Europa-Clipper
   leveraging pair), not a novel structure.
2. Every **literature-fresh** closure involves **Amalthea** (Io-Amalthea-Io
   3.4 cm... 34 m/s; Amalthea-Europa-Amalthea; etc., all 0 anchors) and every
   one **fails the physical max-bend gate**. Amalthea is a tiny inner moon deep
   in Jupiter's gravity well; its flyby cannot bend a heliocentric-class leg
   enough to close the cycle ballistically. The freshness is real (nobody designs
   Amalthea tours) but the geometry is infeasible — the gate correctly rejects it.

So the precise claim: **the repeated-moon quasi-cycler method finds no novel,
physically-feasible, low-Δv (sub-0.05 km/s) quasi-cycler in the Jupiter Galilean
+ Amalthea system at 48×48 phase resolution and n_rev ≤ 5.** The feasible
low-Δv structures are all known; the fresh structures are all infeasible.

This mirrors the Saturn/Neptune/Uranus moontour outcomes: the method reliably
*re-derives* the literature's resonant tours but the unstudied corners
(tiny/inner moons) are unstudied *because* they are dynamically infeasible.

## Disposition

- **No catalogue writeback** (matches #320/#341/#344 discipline; the one
  feasible strong closure is not novel, so there is no admission to make).
- Method-versioned negative registered in `data/empty_regions.jsonl` (Jupiter
  Galilean was unswept by #320 — this fills that gap).
- **Optional follow-on (NOT this task):** Callisto-Ganymede-Callisto (8.3 m/s,
  physical PASS) is a clean *known* resonant quasi-cycler; if catalogue coverage
  of the Jovian system is later desired it is the natural V0/known-attribution
  admission candidate (the Umbriel #339/#340 template, but with a literature
  anchor rather than a novelty claim). Not pursued here — the discovery goal was
  novelty.
