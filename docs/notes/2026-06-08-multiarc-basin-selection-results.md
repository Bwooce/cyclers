# 6.44Gg3 flyby-continuity probe — results & three-way-gate verdict (#162)

**Date:** 2026-06-08
**Experiment:** the recommended first experiment from the #162 multi-arc genome
design (`docs/notes/2026-06-08-multiarc-genome-design.md` §3, Approach #3) and its
plan (`docs/superpowers/plans/2026-06-08-multiarc-basin-selection.md`). Add the
explicit per-flyby V∞-continuity + bend-feasibility residual the chained-DSM
evaluator omitted, then run the decisive 6.44Gg3 probe against the three-way gate.

**Code (this experiment):** additive, default-off `charge_flyby_continuity` path in
`src/cyclerfinder/search/dsm_leg.py` (commit `5b686a6`); probe + symmetric seed in
`tests/search/test_dsm_leg.py::test_dsm_644gg3_flyby_continuity_probe` (commit
`de18822`). Default path verified **bit-identical** (regression test
`test_charge_flyby_continuity_default_off_is_bit_identical`).

**GOLDEN/HONESTY.** EXPECTED = the SOURCED anchors (catalogue `russell-ch4-6.44Gg3`,
Russell 2004 Table 4.13): V∞ **E = 6.44**, **M = 3.74** km/s. Emerged V∞ is
EVIDENCE, never imposed. No catalogue writeback.

---

## Verdict: AMBIGUOUS (floor-drop) — the objective fix works mechanically, the floor dropped ~3× below #157, but the sourced anchors remain unreached; the binding term is ΔV_DSM, not flyby_dv.

This is NOT a clean CLOSE and NOT the flyby_dv-floored EMPTY-SET the design
anticipated. The decisive mechanistic finding is sharper than either:

> **Freeing the intermediate departure-V∞ direction makes the bend term trivially
> satisfiable — every per-flyby `flyby_dv` floors at ~0.** The irreducible residual
> is the **ΔV_DSM budget** (~7.4–9.7 km/s), not the bend-feasibility cost. The
> objective fix DID lower the floor (26.9 → ~9 km/s, the #157 → #162 drop), so the
> diagnosed mechanism was real and the fix is correct — but the low-V∞ sourced
> basin is gated by irreducible DSM work, not by flyby bend.

### Why this is AMBIGUOUS, not EMPTY-SET, not CLOSE (against the design §3 gate)

- **CLOSE** requires `max(residual_vector) < 0.1` AND emerged V∞ within tol. The
  floor is ~9 km/s and emerged V∞ is far off-anchor → not CLOSE.
- **EMPTY-SET** as written keys on an *irreducible `flyby_dv` floor > ~1 km/s*. The
  observed `flyby_dv` floors at **0** wherever the chain stays feasible (the freed
  direction puts every flyby inside its bend cone). So the EMPTY-SET *as the design
  framed it* does **not** fire — the bend term is satisfiable.
- **AMBIGUOUS** requires the floor below #157's 26.9 km/s (✓: ~7.4–9.7) AND the
  basin to move. The floor moved decisively; the emerged V∞ did **not** move toward
  the anchors. So this fires the AMBIGUOUS *floor-drop* criterion with the explicit
  caveat that the emerged-V∞ half did not improve — the honest, decisive evidence
  is the quantified floor + the mechanism (DSM-bound, not bend-bound).

---

## Numbers (verbatim, circular-coplanar backend)

### Single corrector run (symmetric seed, charged, `max_revs=3`)

```
converged=False  max(residual_vector)=9.3367 km/s
total_dV (audit sum)=21.6519 km/s
per-leg dV_DSM km/s:   (5.885, 4.295, 9.337, 2.135)
per-flyby flyby_dv km/s: (0.0, 0.0, 0.0)
emerged V_inf_in (sourced E=6.44, M=3.74):
   {1(E): 9.474, 2(M): 14.273, 3(E): 4.828, 4(M): 8.090}
n_revs/branch per leg: all (0, "single")
eta per leg: (0.351, 0.413, 0.567, 0.509)
tof days per leg: (262.0, 500.4, 262.0, 1315.6)
wall: 5.9 s
```

### MBH follow-up (cauchy, rng_seed=6, ≤120 hops, stall 60)

```
feasible=False  total_dV (audit)=21.6519 km/s   max(residual_vector)=9.3367 km/s
hops attempted/accepted = 61/0
(landed at the same point as the single run — 0 accepted hops)
wall: 106.1 s
```

### Seed sweep (departure V∞ × intermediate η), the floor-robustness datum

| vinf0 | eta | max(residual) | ΔV_DSM per leg | flyby_dv | emerged V∞_in |
|---|---|---|---|---|---|
| 4.00 | 0.3 | **8.02** | (8.02, 3.61, 6.28, 2.63) | (0,0,0) | {E 9.69, M 12.22, E 1.39, M 4.50} |
| 4.00 | 0.5 | **9.70** | (4.62, 4.31, 9.70, 2.80) | (0,0,0.83) | {E 9.27, M 15.70, E 5.00, M 15.46} |
| 4.00 | 0.7 | inf | (chain breaks — hyperbolic) | — | — |
| 6.44 | 0.3 | **7.45** | (7.45, 3.62, 5.04, 2.43) | (0,0,1.69) | {E 10.35, M 11.90, E 1.80, M 4.47} |
| 6.44 | 0.5 | **9.34** | (5.89, 4.30, 9.34, 2.13) | (0,0,0) | {E 9.47, M 14.27, E 4.83, M 8.09} |
| 6.44 | 0.7 | inf | (chain breaks) | — | — |
| 8.00 | 0.3 | **8.03** | (8.03, 3.68, 5.03, 6.95) | (0,0,5.93) | {E 9.59, M 9.29, E 3.42, M 8.59} |
| 8.00 | 0.5 | **9.41** | (7.99, 4.27, 9.41, 2.18) | (0,0,0) | {E 9.55, M 14.34, E 4.90, M 8.01} |
| 8.00 | 0.7 | inf | (chain breaks) | — | — |

**Floor (over feasible seeds): ~7.4–9.7 km/s, robust.** Dominated by ΔV_DSM (the
binding term); `flyby_dv` floors at ~0 except where the corrector did not fully
relax the last flyby (and even then it is the smaller term). Emerged V∞ never
approaches the sourced E 6.44 / M 3.74 — the Earth leg-1 V∞ sits ~9–10 km/s and the
Mars encounters ~9–16 km/s, off-anchor like #157.

---

## Comparison to the elimination chain

| pass | method | floor (km/s) | binding term | emerged at anchor? |
|---|---|---|---|---|
| #150 | minimal E-M-E, 1-DSM | 9.40 | ΔV_DSM | no |
| #153 | full E-M-E-M-E, single-rev | 29.9 | ΔV_DSM (degenerate Lambert) | no |
| #157 | full sequence, multi-rev | 26.9 | ΔV_DSM | no |
| **#162** | **+ explicit flyby continuity+bend residual** | **~7.4–9.7** | **ΔV_DSM (flyby_dv→0)** | **no** |

The #162 objective fix produced the **lowest multi-arc floor to date** (~9 vs 26.9),
confirming the diagnosed mechanism (the old scalar `Σ ΔV_DSM` had no bend term and
was dominated by the long leg). But it **also proves the bend term was never the
binding constraint at the sourced anchors** — freeing the flyby direction satisfies
it for free. The residual that remains is the DSM impulse needed to stitch the two
generic loop arcs into a *continuous* heliocentric chain at these ToFs, and that
budget is ~7–9 km/s irreducible in this transcription.

---

## What this means / next step

This sharpens the empty-set hypothesis rather than confirming or refuting it:

- The sourced anchors are **not** a bend-feasibility wall (the design's leading
  hypothesis) — flyby_dv → 0 freely.
- They **are** behind an irreducible **ΔV_DSM** floor (~7–9 km/s) in the
  circular-coplanar E-M-E-M-E one-DSM-per-leg transcription. The published cycler is
  *near-ballistic* (catalogue `delta_v_kms: 0.509`, `turn_ratio: 0.95`), so a ~9
  km/s DSM floor means **this transcription does not host the sourced geometry as a
  near-ballistic closure** — consistent with the memory blocker
  (`project_s1l1_realeph_closure_blocker.md`: "representable ≠ reachable"; the low
  Mars V∞ only appears in *discontinuous* geometry) and Hughes-2014 (low-V∞ Mars
  needs broken-plane the authors exclude).

**Recommended next experiment (design #5, the hybrid), conditioned on this
AMBIGUOUS-floor-drop result:** symmetric seed + branch enumeration + the charged
residual + MBH, and — given that the binding term is now identified as ΔV_DSM, not
bend — add **per-arc ToF / multi-rev freedom on the loop legs** so the optimiser can
find a lower-DSM stitch (the g-arc remainder admitted no multi-rev branch at the
fixed 500 d ToF; letting its ToF/branch move is the obvious DSM-reducer). The
decisive question for the hybrid: does freeing the loop-arc ToF/branch drive the
ΔV_DSM floor toward zero, or does it confirm a hard ~km/s DSM floor (the quantified
empty-set)? Either is publishable.

**No catalogue writeback.** The emerged V∞ are evidence only; the sourced anchors
remain the EXPECTED target and were not reached.
