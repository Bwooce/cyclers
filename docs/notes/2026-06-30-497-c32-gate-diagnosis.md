# #497 — #249 C32-dominance gate: diagnosed, not flippable by a cap value

**Date:** 2026-06-30. **Verdict:** the `DV_CAP_MS≈51 m/s` recalibration does NOT flip
`test_validation_gate_c32_dominant`; the gate failure is PRECISELY isolated to OUR proxy-ΔV
fidelity. A NEW passing golden gate validates the centrality math against the published matrix.

## What #495 suggested vs what the data shows
#495 noted "C32 wins all three at DVcap≈51 m/s in Braik's budget sweep" → hypothesis: recalibrate
`DV_CAP_MS` 409.3 → 51 m/s to flip the xfail. **Tested directly — the hypothesis fails:**

- **Our proxy network at 51 m/s = EMPTY** (all centralities 0). Our heading-fan proxy ΔVs all
  exceed 51 m/s because the proxy OVERESTIMATES ΔV (#495: dc_refined < proxy in every pair). So a
  51 m/s cap on our (inflated) proxy drops every edge.
- **Our proxy at 409.3 m/s**: C32 does NOT dominate (betweenness argmax = R21-U; C32 zero
  betweenness) — a ranking error from the proxy ΔV values.

So no cap value flips the gate from OUR proxy: the obstacle is the proxy-ΔV scale + ranking, not the
budget.

## The clean separation (the #497 deliverable)
Fed Braik & Ross's OWN published proxy-ΔV matrix (`data/golden/braik_ross_2026_dvmatrix_mps.csv`,
adopted #495), our `normalized_centralities` **reproduces Table 4 to print precision**:

| metric | our value | Braik Table 4 |
|---|---|---|
| C32 strength | 0.2850 | 0.2850 |
| C32 harmonic | 0.2891 | 0.2891 |
| C32 betweenness | 0.5000 | 0.5000 |

C32 is the argmax of all three. → **the centrality scorer is CORRECT** (validated against published
data). New passing test: `test_centrality_scorer_reproduces_braik_ross_table4` (fast, sourced golden).

## Net
- ADDED: a passing, sourced golden gate — the centrality math reproduces Braik Table 4 exactly.
- UPDATED: `test_validation_gate_c32_dominant` xfail reason from "not yet diagnosed" → the precise
  diagnosis (OUR proxy-ΔV overestimates; empties at 51 m/s, mis-ranks at 409.3). The xfail STANDS
  (honestly — our proxy genuinely doesn't reproduce it); it is NOT tuned to pass.
- The real fix (open, the remaining #249/#497 work): rebuild the gate network from the adopted
  `dc_refined` golden (Braik's DC-refined ΔV, at the right scale) OR calibrate the proxy→refined
  ΔV scale (#495 found a strong ρ=0.96 monotone relation, so a calibration factor is plausible).
  Not attempted here. Parameters were NOT tuned to manufacture a pass.
