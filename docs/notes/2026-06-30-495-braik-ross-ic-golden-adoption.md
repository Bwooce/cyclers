# Task #495: Braik-Ross IC golden adoption — verdict note

**Date:** 2026-06-30
**Task:** Adopt sourced ICs from Braik & Ross 2026 (MIT-licensed repo) to unblock #249 gate.
**Source:** Abdullah Braik & Shane D. Ross, *Orbital Networks in the Three-Body Problem*,
arXiv:2605.31543 (2026). Repo: https://github.com/BinBraik/cislunar-orbital-network (MIT).

---

## 1. What was extracted

All 13 family ICs from `src/cr3bp_family_ic.m` (function `cr3bp_family_ic`) extracted
verbatim into `data/golden/braik_ross_2026_em_family_ics.yaml`. These are independently
computed (via `correct_po_to_cj313_v2.m`, verified against Table 2 instability rates and
periods) — NOT derived from any code in this repository.

CSV goldens committed:
- `data/golden/braik_ross_2026_dvmatrix_mps.csv` — 13×13 proxy ΔV matrix (m/s)
- `data/golden/braik_ross_2026_dc_refined_mps.csv` — 78-pair DC-refined ΔV+TOF
- `data/golden/braik_ross_2026_tofmatrix_days.csv` — 13×13 proxy TOF matrix (days)

Jacobi self-consistency check passed for all 13 families (CJ computed from stored
(x0, ydot0) matches stored CJ to ≤ 2e-13, essentially machine epsilon).

---

## 2. C11a and C21 recovery at exact Braik ICs

### C11a (CJ = 3.1294, x0 = -0.8116406668238326)

Corrector (`correct_symmetric_fixed_jacobi`, half_crossings=3):
- **Converged: YES**
- x0 match vs Braik IC: **exact** (0.00 difference)
- ydot0 match: **4.2e-17** (sub-machine-eps)
- Period: **42.140 d** (sourced 42.140 d, diff < 0.001 d)
- **Result: C11a recovers AT the Braik IC.** The corrector treats the exact IC as a
  fixed point — it is already a solution to the boundary-value problem.

### C21 (CJ = 3.129389531054557, x0 = 7.237366530581342e-01)

Corrector at **exact Braik CJ** (`correct_symmetric_fixed_jacobi`, half_crossings=4):
- **Converged: YES**
- x0 shift from Braik seed: 1.10e-07 (tiny; family is nearly degenerate in C)
- Period: **84.533 d** (sourced 84.533 d, diff 0.0002 d)
- **Result: C21 recovers at the exact Braik CJ.** Period match is well within 0.5 d gate.

Corrector at the **rounded CJ=3.1294** (the original #236 bug):
- The corrector at 3.1294 converges to a DIFFERENT orbit — period deviates > 5 d from
  the sourced 84.533 d (it slides onto an off-family basin). This confirms that the
  original #236 failure was exactly what the memory note describes: the 1.05e-5 gap
  between the printed 3.1294 and the family's CJ-max placed the seed off-family.

### CJ sourcing comparison

| Source | C21 CJ |
|---|---|
| Braik & Ross 2026 repo (MIT) | **3.129389531054557** |
| Ross-RT 2025 AAS-25-621 Table 4 | 3.129389531088256 |
| Difference | 3.3e-11 (12th decimal place) |

Both sources agree to 11 significant figures. The Braik-repo value is now the PRIMARY
source (MIT-licensed, independently computed). `C_J_C21_AAS` retained as a cross-reference
alias.

---

## 3. Proxy ΔV cross-check (DVmatrix vs dc_refined)

Cross-checking Braik's own output data at the pair level:

| Metric | Result |
|---|---|
| Pairs checked | 75 (of 78 in dc_refined; 3 NaN refined skipped) |
| Pairs where dc_refined > proxy (DVmatrix) | **33 / 75 (44%)** |
| Pearson r (proxy vs refined) | **0.989** |
| Spearman rho (proxy vs refined) | **0.962** |
| Worst violation | R31-U ↔ R52-U: proxy=35.1 m/s, refined=100.8 m/s (+65.7 m/s) |

**Conclusion: the proxy ΔV (DVmatrix) is NOT a guaranteed upper bound on the DC-refined ΔV.**
This is not a defect in our implementation — it is inherent to the voxel-overlap proxy:
the heading-turn + approximate patch cost estimated from the voxel does not bound the
differential-correction patch burn needed to connect the actual trajectory segments.

The proxy IS however a very strong SCREENING TOOL: Pearson r=0.99 and Spearman rho=0.96
show that the proxy ranks pair accessibility reliably. Small proxy → small refined (in the
ordering sense), making the proxy useful for family prioritization and pair pre-screening.

The worst violations cluster in geometrically "hard" pairs (stable resonant orbits:
R52-S, R21-S, R31-S with C32, R52-U, etc.) where the voxel-level heading match
significantly underestimates the actual trajectory correction needed.

**Implication for #249 gate**: the proxy ΔV from our `proxy_matrix` function should be
understood as a screening estimate, not a conservative bound. The gate test
`test_validation_gate_c32_dominant` verifies structural network properties (hub/relay
roles), not absolute ΔV bounds.

---

## 4. Budget-cap finding (C32-dominant reproduction)

From cross-analysis of `snapshot_summary.csv` (not committed; Braik repo only):

- C32 wins **strength, harmonic closeness, AND betweenness** in 500/625 budget-parameter
  combinations (80% of the 25×25 DVcap × Tmax sweep).
- C32 wins all three simultaneously at di≥1, dj≥6, i.e., DVcap=51.16 m/s, Tmax≥11.1 d.
- At di=1, dj=1 (the tightest budget: DVcap=51.16 m/s, Tmax=6.83 d), C11a wins
  strength/closeness and C32 wins only betweenness.

**Root cause of our `xfail` gate (test_validation_gate_c32_dominant):**
Our gate uses `DV_CAP_MS=409.3 m/s` (full-connectivity reference). At this uncapped budget
nearly all edges exist, the network is near-complete, and betweenness carries no signal
(every node reaches every other directly). C32's gateway/relay role only emerges under the
BUDGET-CONSTRAINED formulation where some direct edges drop out and relay routing appears.

To reproduce the C32-dominant paper result, the scorer needs:
1. Budget cap ~51 m/s (the Braik reference scale visible in DVmatrix_mps.csv)
2. Tmax ≥ 11 d (the reference orbital period used in the parametric sweep)

The gate stays `xfail` until `DV_CAP_MS` is recalibrated to the Braik reference scale.

---

## 5. What #249 work this unblocks

**C11a and C21 both now recover** at the exact Braik ICs. Combined with the #262 recovery
of all 4 cyclers (C11a, C11b, C21, C32) and 8 offline nodes (LL1, LL2, DPO, R21-S, R21-U,
R31-S, R31-U, R52-S), the 12-node source-confirmable gate set is available.

The **remaining blocker** for the full 13-node gate (and C32-dominant reproduction) is:
1. R52-U is still not recovered (period-confirmed offline recovery fails at CJ=3.1294;
   the Braik IC for R52-U is now in the YAML golden at x0=-0.2719428329684943).
2. The DV_CAP_MS parameter in `reachable_network.py` needs recalibration to ~51 m/s.

**V-tier recommendation (for adjudication, not enacted here):**
- The 12 recovered families (all except R52-U) are confirmed to sourced Table-2 periods
  and sigmas via independent corrector runs. No catalogue writeback recommended — these
  are accessibility-network members at a non-catalogue common energy, not new catalogue rows.
- If the C32-dominant gate passes after DV_CAP_MS recalibration, the scorer is licensed
  for family-selection prioritization on our catalogue continuation seeds.

---

## 6. Files changed / committed

- `data/golden/braik_ross_2026_em_family_ics.yaml` — 13-family IC golden (NEW)
- `data/golden/braik_ross_2026_dvmatrix_mps.csv` — proxy ΔV matrix (NEW)
- `data/golden/braik_ross_2026_dc_refined_mps.csv` — DC-refined ΔV+TOF (NEW)
- `data/golden/braik_ross_2026_tofmatrix_days.csv` — proxy TOF matrix (NEW)
- `tests/search/test_braik_ross_ic_golden.py` — 9 golden tests (NEW, all pass)
- `src/cyclerfinder/search/reachable_representatives.py`:
  - `C_J_C21` updated from AAS-25-621 value (3.129389531088256) to Braik-repo value
    (3.129389531054557); old value retained as `C_J_C21_AAS` alias
  - `ROSS_C21_JACOBI` aliased to `C_J_C21` (was redundant cross-reference)
  - C11a x0 seed updated from -0.81164067 to -0.8116406668238326 (Braik exact)
  - C21 x0 seed updated from +0.7237335857 to 7.237366530581342e-01 (Braik exact)
  - Several OFFLINE_SEEDS x0 seeds updated to Braik exact values
