# #497 — proxy-fidelity rebuild: the 30-60× overestimate FIXED; C32 xfail is a missing-node effect

**Date:** 2026-07-01. **Verdict:** the 30-60× proxy-ΔV overestimate is a real,
localized **bug** in the heading-mismatch patch term — now fixed (median our/Braik
edge ratio 30-60× → **~1.9×**, Spearman ρ=0.84). The C32-dominance gate still does
not flip, but the corrected root cause is **not** proxy fidelity: our 12-node
source-confirmable set **excludes R52-U**, and C32-dominance is provably contingent
on R52-U **in Braik's own published matrix**. The xfail STANDS, re-diagnosed. No
parameters were tuned to force C32.

Supersedes the proxy-scale claims in `2026-06-30-497-c32-gate-diagnosis.md` /
`2026-06-30-497b-proxy-calibration-verdict.md` (which measured the *buggy* proxy).

---

## 1. Root cause of the 30-60× overestimate — the constant-pedestal patch (FIXED)

The pair proxy is `min over shared voxels of (src_turn + tgt_turn + patch)`. The
`src`/`tgt` turn machinery was already correct (hard-access R21-S rows matched Braik
to ratio ~1). The defect was `_patch_cost`:

```python
# OLD (buggy):  a CONSTANT, using a UNIT reference speed and half the coarse voxel
return dv_turn(1.0, 0.5 * grid.dtheta)   # = 89.3 m/s at dtheta=10°, ADDED to every pair
```

Two compounding errors vs Braik & Ross Sec. 4 (whose `DVpatch` column in
`dc_refined_mps.csv` is 1.5–9 m/s and **voxel-varying**):

1. **Unit reference speed** `v=1.0` instead of the arc's true local rotating-frame
   speed at the voxel (~0.2–0.5 nd on the C_J manifold) → ~2–5× too big.
2. **Fixed half-voxel heading** `dtheta/2` at the coarse 10° grid instead of the
   **actual** heading mismatch between the two arcs (≤ one voxel, usually far less)
   → ~10× too big at 10°.

Because the term used `v_ref=1.0` and `dtheta` only, it was **constant across all
voxels and all pairs** — an 89.3 m/s pedestal added to every edge. Empirical
decomposition (12 recovered members, gate grid) confirmed the mechanism exactly:

| pair | our ΔV (m/s) | min(src+tgt) | patch | Braik proxy | ratio |
|---|---|---|---|---|---|
| C21–C32   | 89.3 | **0.0** | 89.3 | 3.26 | 27× |
| C11a–C32  | 89.3 | **0.0** | 89.3 | 3.40 | 26× |
| C11a–C21  | 89.3 | **0.0** | 89.3 | 2.86 | 31× |
| R21-U–C11a| 89.3 | **0.0** | 89.3 | 2.13 | 42× |
| C11b–C32  | 155.8| 66.5    | 89.3 | 3.27 | 48× |
| R21-S–C21 | 326.8| 237.5   | 89.3 | 327.84 | 1.0 |
| R31-U–R52-S| 121.6| 32.4   | 89.3 | 123.70 | 1.0 |

For the low-cost cycler pairs `min(src+tgt)` is **exactly 0** (the natural,
un-maneuvered orbits already overlap) — the physics is right — so the entire 26–48×
overestimate is the pedestal. For the hard pairs the turn cost dominates and the
ratio is already ~1. This is a pure **metric-definition/units bug**, not a units
error in `VU_MS` (verified 1023.16 m/s, correct for Earth-Moon) and not (primarily)
a grid-resolution effect.

### The fix (default, correctness — [[feedback_bugfix_invalidates_past_searches]])

`ReachableSet` now stores the actual heading and speed of the cheapest arc per
voxel; `mirror_reachable_set` time-reverses the heading (π−θ) and preserves speed;
`pair_proxy` computes the physical patch:

```python
v_local  = forward_a.speeds[idx]
mismatch = angular_diff(forward_a.headings[idx], backward_b.headings[idx])
patch    = dv_turn(v_local, mismatch)     # local speed × ACTUAL mismatch
```

**Result (gate grid dx=0.02, dθ=10°, 52 common finite edges):**

| | before fix | after fix |
|---|---|---|
| median our/Braik edge ratio | ~19–48× | **1.94×** (p25 0.83, p75 3.89) |
| Spearman ρ (our vs Braik proxy) | — | **0.844** (p=4e-15) |

The proxy now tracks Braik's published proxy in both scale and rank.

---

## 2. Why the C32 gate STILL does not flip — the missing R52-U node (not proxy)

After the fix, C32 is still **not** the dominant node on our 12-set — **C11a** is
the strength/closeness hub, at every grid resolution tested (dx swept 0.02→0.001,
dθ 10°→1°; C32 best strength rank = 2). Grid resolution is therefore **not** the
residual cause.

The cause is the **node set**. Our source-confirmable set excludes the 5:2 unstable
resonant **R52-U** (never recovered — σ=0.37 collapses onto spurious orbits with the
available single-/free-period correctors). In Braik's published matrix, **C32–R52-U
= 0.62 m/s is the single smallest edge in the whole 13×13** — R52-U's strongest
connection is to C32, and it is a large share of C32's hub strength.

**Proven on Braik's OWN published matrix** (fast, sourced, no propagation —
`test_braik_matrix_c32_dominance_requires_r52u`):

| node set | strength argmax | harmonic argmax | C32 strength rank |
|---|---|---|---|
| full 13 (published) | **C32** | **C32** | 1 |
| drop R52-U → our 12 | **C11a** | **C11a** | **3** |

Removing R52-U from Braik's own data reproduces our scorer's ranking (C11a hub, C32
demoted) exactly. So the negative is a **node-set incompleteness**, not a proxy
defect. Our fixed proxy on the 12-set gives the same hub (C11a) and the same relay
(C32 betweenness) as Braik's matrix restricted to the same 12 nodes — a positive
control that the proxy is now faithful.

---

## 3. Disposition

| test | status | note |
|---|---|---|
| `test_centrality_scorer_reproduces_braik_ross_table4` | PASS | scorer correct (Table 4) |
| `test_dc_refined_gate_c32_dominant` | PASS | C32 dominant on dc_refined |
| `test_braik_matrix_c32_dominance_requires_r52u` (NEW) | PASS | root cause pinned: R52-U |
| `test_proxy_calibration_is_unit_slope_confirms_clean_negative` | PASS | (historical, buggy-proxy era) |
| unit tests (`test_reachable_network.py`) | PASS | physical patch asserted |
| `test_validation_gate_c32_dominant` | **xfail (STANDS)** | re-diagnosed: missing R52-U |
| `test_validation_gate_c32_undominant_faithful_negative` | PASS (slow) | C11a hub, matches Braik-12 |

**xfail stays.** The proxy bug is fixed and the scorer now tracks Braik (ρ=0.84,
1.9×), but the published C32-dominance is contingent on R52-U, which we cannot
source-confirm. Flipping the gate faithfully requires **recovering R52-U** with a
robust Jacobi-constrained multiple-shooting corrector — a separate open task — NOT a
parameter tune ([[feedback_golden_tests_sourced_only]]). No catalogue edits.

## 4. Follow-on

Recover R52-U (5:2 unstable resonant) at C_J=3.1294 to its Table-2 period (56.436 d,
σ=0.37) via a Jacobi-constrained multiple-shooter; add it to `_recover_subset`; then
re-evaluate the gate. With the proxy now faithful (§1) and the mechanism pinned
(§2), R52-U recovery is the single remaining blocker to a faithful C32-dominance
reproduction on OUR members.
