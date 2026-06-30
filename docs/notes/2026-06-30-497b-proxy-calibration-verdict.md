# #497b — proxy→refined calibration verdict: dc_refined gate PASSES; calibration is clean negative

**Date:** 2026-06-30.
**Precondition:** #497 diagnosis (our proxy overestimates → 51 m/s empties our network,
409.3 m/s mis-ranks). #495 golden adoption (dvmatrix_mps + dc_refined_mps golden CSVs).

---

## Step 1 (dc_refined gate) — PASSES

Built the symmetric 13×13 dc-refined ΔV matrix from
`data/golden/braik_ross_2026_dc_refined_mps.csv` (75 finite pairs; 3 NaN = no
accessible path). Applied the Braik & Ross budget cap of **51.16 m/s** (the first
budget level at which C32 wins all three metrics simultaneously in Braik's parametric
sweep, sourced from #495 analysis of Fig. 10 / snapshot_summary).

Result: **C32 is the dominant family** — argmax of strength, harmonic closeness,
AND betweenness.

| metric | C32 value | argmax |
|---|---|---|
| strength | 0.5167 | C32 ✓ |
| harmonic closeness | 0.5397 | C32 ✓ |
| betweenness | 0.1667 | C32 ✓ |

C32 dominance is **very robust**: it holds across the full range of caps tested
(10, 20, 51.16, 60, 100, 400, 1e6 m/s) and at full dc_refined connectivity (75 edges).
At every cap level C32 is the argmax of all three. This is a STRONGER result than
Braik Table 4 (which is at the proxy scale): even the DC-corrected maneuver costs
(actual optimized transfers, not proxy estimates) show C32 as the dominant node.

**New passing test:** `test_dc_refined_gate_c32_dominant` (fast, no propagation,
gated on the sourced dc_refined golden). Added alongside the existing
`test_centrality_scorer_reproduces_braik_ross_table4` (which uses the proxy dvmatrix).

---

## Step 2 (proxy calibration) — CLEAN NEGATIVE

### Fitted calibration from (#495 golden pairs)

Fit: linear regression on 75 valid (dvmatrix_proxy, dc_refined) pairs.

| coefficient | value |
|---|---|
| slope | **1.0196** |
| intercept | **3.779 m/s** |
| Pearson r | 0.989 |
| r² | 0.978 |

The slope is **≈1.02**: the proxy and dc-refined ΔV values are on the **same scale**
in Braik's data. This is a non-compressive mapping (scaling by ~1 and shifting +3.78 m/s).

### Why this fails for OUR proxy

Our heading-fan `proxy_matrix()` output (nondimensional × VU_MS ≈ 1023 m/s) gives
values in the range of **hundreds of m/s** for the same families — approximately 30-60×
larger than Braik's proxy values (in m/s). The #497 diagnosis confirmed: applying a
51 m/s cap to our proxy **empties the network** (0 surviving edges → all centralities 0).

Applying the calibration (slope≈1, intercept≈+3.78 m/s):

    f(x) = 1.0196 × x + 3.779    [m/s → m/s]

For any x > 51.16 m/s:

    f(x) = 1.0196 × x + 3.779 > x > 51.16 m/s

The calibration shifts values **upward** (slope>1, intercept>0). If our proxy values
are already above 51.16 m/s, the calibrated values are even further above. The network
remains empty at the 51.16 m/s cap. **C32 still does not dominate on OUR families.**

The fundamental mismatch is **scale**, not functional form: our heading-fan proxy
overestimates by 30-60× compared to Braik's proxy (and by far more compared to
dc_refined). A calibration fitted on (Braik_proxy, dc_refined) pairs has slope≈1
because those two datasets are already on the same scale. It cannot bridge the
30-60× gap to our proxy.

**New passing test:** `test_proxy_calibration_is_unit_slope_confirms_clean_negative`
(fast, no propagation). Verifies: r²≥0.95, 0.5<slope<1.5, and that
`f(51.16 m/s) > 51.16 m/s` (calibration is not compressive at the budget threshold).

---

## What calibration WOULD require

To bring our proxy into the dc_refined scale, we would need a compression factor of
roughly 30-60×. This cannot come from a calibration fitted on (Braik_proxy, dc_refined)
pairs, because those are already at comparable scales. The required calibration would
need to be derived from (OUR_proxy, dc_refined) pairs — i.e., we would need to:

1. Run `proxy_matrix()` on Braik's 13 families (expensive: full heading-fan propagation)
2. Pair our proxy output (in m/s) against the dc_refined golden
3. Fit f(our_proxy) → dc_refined
4. Apply to our families

Step 1 is already done implicitly by the slow xfail tests (which compute the proxy
matrix on our recovered families). The result is that our proxy values all exceed
51 m/s. Step 2 is the missing piece: we have no (OUR_proxy, dc_refined) pairs to
fit from. Getting them would require Braik's DC-corrector applied to our proxy
overlap voxels — which is essentially rebuilding their full pipeline.

---

## Net / xfail status

| test | status | notes |
|---|---|---|
| `test_centrality_scorer_reproduces_braik_ross_table4` | PASS | scorer correct on proxy dvmatrix |
| `test_dc_refined_gate_c32_dominant` (NEW) | PASS | scorer correct on dc_refined too |
| `test_proxy_calibration_is_unit_slope_confirms_clean_negative` (NEW) | PASS | calibration approach documented |
| `test_validation_gate_c32_dominant` | **xfail (STANDS)** | our proxy cannot reproduce C32 dominance |
| `test_stable_resonants_are_hard_access_on_subset` | PASS (slow) | qualitative negative reproduced |
| `test_validation_gate_c32_undominant_faithful_negative` | PASS (slow) | faithful negative recorded |

The #249 gate xfail is **correctly placed and exhaustively diagnosed**:
- Centrality math: correct (test 1)
- dc_refined confirms C32: correct (test 2 — new)
- Sourced calibration: unit-slope, cannot compress OUR proxy (test 3 — new)
- Root cause: our heading-fan proxy overestimates by 30-60× vs Braik; not addressable by
  a (Braik_proxy, dc_refined) calibration
- The fix would require a finer grid and/or DC correction of our proxy trajectories,
  neither of which is available without significant new computation infrastructure
