# Ross & Roberts-Tsoukkas 2025 adoption — CR3BP-lane results (#212 Part B)

**Date:** 2026-06-12
**Source:** S. D. Ross & M. Roberts-Tsoukkas, "Stable, Low-Energy Prograde
Earth-Moon Cycler Orbits," AAS 25-621 (2025). Mining note:
`docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md`.
**Code:** `src/cyclerfinder/search/cr3bp_periodic.py`
(`ydot0_from_jacobi`, `correct_symmetric_fixed_jacobi`, `barden_stability`).
**Reproduction:** `scripts/cr3bp_ross_reproduce.py`.
**Tests:** `tests/search/test_cr3bp_ross_families.py`.
**Status:** PROPOSED for review. NO catalogue writeback (that is #216).

This pass adopts the paper's construction into the CR3BP lane:
- a **fixed-Jacobi symmetric corrector** (`correct_symmetric_fixed_jacobi`,
  Ross Eqs. 9-12) — single-shooting that holds the Jacobi constant fixed and
  solves the 1-D problem for `x0` so the orbit re-crosses the x-axis
  perpendicularly at the half period. The half-period crossing is selected by
  its **index** among the many x-axis crossings of a multi-rev cycler (the
  crossing nearest `T/2`), held fixed across Newton iterations — this is the
  step that lets the single-shooter land on the cycler branch;
- **Barden half-period monodromy stability** (`barden_stability`, Eqs. 13-15):
  `M = G inv(Phi(T/2)) G Phi(T/2)`, `nu = 1/2(lambda + 1/lambda)`, stable iff
  `|nu| < 1`, evaluated on the **4x4 planar sub-STM** (indices 0,1,3,4) so the
  trivial out-of-plane z eigenpair cannot be mistaken for the nontrivial pair.

Model: **pure planar CR3BP (PCR3BP)**, `mu = 1.2150584270572e-2` (paper p. 3).
The corrector and stability index are the engine the CR3BP continuation campaign
(`docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md`)
was gated on.

## 1. Five-family reproduction (published vs reproduced)

Acceptance run: `scripts/cr3bp_ross_reproduce.py`. EXPECTED values are Ross's
PRINTED `C^stable` / `T^stable` (Table 3, p. 11) and the stability verdict; the
recovered `x0` / `ydot0` are DERIVED (the 1-D solve of the symmetric-orbit
structure, mining note §5), never goldens. `C` is enforced algebraically so it
matches to machine epsilon; the binding acceptance checks are the **period** and
the **stability verdict**.

| family | C^stable (pub = reproduced) | T^stable pub (TU) | T reproduced (TU) | dT (TU) | nu | verdict |
|--------|------------------------------|-------------------|-------------------|---------|------|---------|
| (1,1) | 3.151175879508174 | 10.29206921007976 | 10.2920692262528 | +1.6e-08 | -0.00334 | **STABLE** |
| (2,1) | 3.129389531088256 | 19.44043166795154 | 19.4401604299533 | -2.7e-04 | +0.05007 | **STABLE** |
| (3,1) | 3.161784147013429 | 14.78849241668140 | 14.7882679408449 | -2.2e-04 | +0.01545 | **STABLE** |
| (3,2) | 3.182762663084288 | 17.90058010350006 | 17.9005801151158 | +1.2e-08 | -0.01175 | **STABLE** |
| (3,3) | 3.177224018696528 | 18.14546057589189 | 18.1454605742454 | -1.6e-09 | +0.06001 | **STABLE** |

**All five** published families reproduce, every one **linearly STABLE**
(`|nu| < 1`) and landing near the `nu = 0` midpoint, exactly the Table-3
definition (the published member is the midpoint of the largest stable
subfamily, p. 11). These are the **first `|nu| < 1` verdicts the CR3BP lane has
ever produced** — all prior lane discoveries (the 14 SILVER Saturnian
Lyapunovs) were unstable, error-amplifying orbits. A regression test pins the
direction (`barden_stability` returns `|nu| > 1` for the unstable Arenstorf
figure-eight) so a stable-only suite cannot pass vacuously.

Reproduced derived states (nd, perpendicular x-axis-crossing IC
`(x0, 0, 0, 0, ydot0, 0)` — DERIVED, stored with provenance, NOT goldens):

| family | x0 | ydot0 | T (days) |
|--------|-----|--------|----------|
| (1,1) | -0.7682140805 | -0.2568154857 | 44.753801 |
| (2,1) | +0.7237335857 | +0.4137707374 | 84.533154 |
| (3,1) | -0.3209891696 | -1.8322134075 | 64.304970 |
| (3,2) | -0.3214486013 | -1.8239555099 | 77.838478 |
| (3,3) | -0.3217380626 | -1.8238865086 | 78.903311 |

(1 TU = 27.321661 d / 2pi = 4.348377... d, paper p. 3.)

**Residual provenance.** `dT` for (1,1), (3,2) and (3,3) is at the 1e-8/1e-9
level (wide stable windows). For the razor-thin windows (2,1) (`Delta_p_m`
4.23 km) and (3,1) (narrow inner band) `dT ~ 2e-4`: the exact `nu = 0` `x0`
sits a hair off the published 16-digit `C^stable`, so the corrected member is
marginally displaced along the family and its period lands ~2e-4 TU (~1e-3 d)
from the printed `T^stable`. The verdict (STABLE, `nu ~ 0`) and the enforced
`C` are exact. Independent **Radau** cross-check (different integrator) closes
the full period and conserves Jacobi to `dJ < 1e-12` for all five.

## 2. (3,2): the half-period crossing is the 6th x-axis crossing

The (3,2) Table-3 member (`C^stable = 3.182762663084288`,
`T^stable = 17.90058010350006` TU) is the `nu = 0` midpoint of the **larger of
two** stability windows (p. 13). Its perpendicular-crossing branch is a
high-velocity Earth-side family (`|ydot0| ~ 1.82`) that makes many x-axis
crossings per period; the binding subtlety is **crossing-index selection** — the
half-period perpendicular crossing is the **6th** x-axis crossing
(`half_crossings=6`), not the index a naive nearest-`T/2` guess picks. Seeded
from the same `x0 ~ -0.321` region as (3,1)/(3,3) with `ydot0_sign = -1` and the
correct index, the corrector lands `T = 17.9005801151` (`dT = +1.2e-8`),
`nu = -0.01175` (STABLE), and the independent Radau cross-check closes with
`dJ = 5.3e-13`. (An earlier pass deferred this family as
"integrator-cost-prohibitive"; fixing the crossing index resolves it — the
deferral is retracted.)

## 3. EM CR3BP backfill re-run with corrected mu (#212a)

Part A fixed the `cr3bp_system("Earth","Moon")` mu double-count
(commit 3cec84c). Re-running `scripts/cr3bp_backfill.py` (runlog
`2026-06-12T06:07:57Z`) against the three EM catalogue rows:

- **`arenstorf-em-figure8-1963`** — CONVERGED, and the non-dimensional outputs
  are **UNCHANGED** by the mu fix: `jacobi = 2.85641252`,
  `period_nd = 17.0652165601594`, `state0_nd` identical (closure residual
  `7.749e-11`). This is expected and verified: the Arenstorf row uses the
  explicit test-problem `mu = 0.012277471`, NOT the registry EM mu, so the fix
  cannot touch its nd quantities.
- The only delta is the **dimensional EM characteristic time** reported with the
  row (it is derived from the corrected `G(m1+m2)`):

  | quantity | pre-fix | corrected | rel. delta |
  |----------|---------|-----------|------------|
  | EM mu (registry) | 0.0120047200 | 0.0121505844 | +1.215e-2 |
  | EM t_s (s) | 372931.4 | 375190.3 | +6.06e-3 |
  | lunit_km | 384400.0 | 384400.0 | 0 |

  i.e. the Arenstorf row's reported `tunit_s` shifts `372931 -> 375190`; its
  `jacobi_constant`, `period_nd`, `state_nd`, `lunit_km` do not change. (The
  corrected EM mu matches Ross's printed `1.2150584270572e-2` to 1.0e-8
  relative.)
- **`genova-aldrin-2015-em-3petal-cycler`** and **`wittal-2022-em-cycler-family`**
  remain NO_SOURCED_IC (no fabricated ICs) — unaffected by mu.

The committed `docs/notes/2026-06-10-cr3bp-backfill-results.md` (with its
post-mining Genova-Aldrin permanent-out-of-model addendum) is left intact; the
mu-correction delta is recorded here rather than by overwriting that note. NO
catalogue writeback; the `tunit_s` correction is a review item for #216.

## 4. Proposed (state0, T, C, nu) tuples for #216

Five Tier-2 EM rows (mining note §8), one per family, common fields:
`model_assumption: cr3bp`; `mass_ratio: 1.2150584270572e-2` (p. 3); frame
Earth-Moon rotating barycentric planar; `center`: Earth-Moon barycenter;
`source_ephemeris: n/a` (pure CR3BP); `lunit_km: 384400.0`;
`tunit_s: 375699.8`; `tof_days_bounds: [period, period]`. Source for every
number: Ross & Roberts-Tsoukkas 2025, Table 3 (p. 11) + per-family text.

| proposed id | C^stable | period (TU) | state0 (nd, DERIVED) | nu | validation |
|---|---|---|---|---|---|
| ross-rt-em-cycler-11-2025 | 3.151175879508174 | 10.29206921007976 | [-0.7682140805,0,0,0,-0.2568154857,0] | -0.00334 (stable) | reproduced same-model |
| ross-rt-em-cycler-21-2025 | 3.129389531088256 | 19.44043166795154 | [0.7237335857,0,0,0,0.4137707374,0] | +0.05007 (stable) | reproduced same-model |
| ross-rt-em-cycler-31-2025 | 3.161784147013429 | 14.78849241668140 | [-0.3209891696,0,0,0,-1.8322134075,0] | +0.01545 (stable) | reproduced same-model |
| ross-rt-em-cycler-32-2025 | 3.182762663084288 | 17.90058010350006 | [-0.3214486013,0,0,0,-1.8239555099,0] | -0.01175 (stable) | reproduced same-model |
| ross-rt-em-cycler-33-2025 | 3.177224018696528 | 18.14546057589189 | [-0.3217380626,0,0,0,-1.8238865086,0] | +0.06001 (stable) | reproduced same-model |

Row extras to carry: `stability: linearly stable (|nu|<1, nu~0 midpoint;
nu = 1/2(lambda+1/lambda))`; the family C^max values (Table 3); Case 2/Case 3
energy-regime flag (`C2 = 3.172160450399808`). Promotion guidance: all five
reproduced members closed in-model on `(C^stable, T^stable)` with the published
stability verdict and an independent Radau cross-check — eligible to promote
above family-seed per the ladder.

## 5. Data gaps carried forward (do NOT silently resolve)

Two `C_(k1,k2)` **bound**-column values in Table 3 are internally inconsistent
(mining note §4; re-verify against the 2026 journal). Pinned by tests
(`test_data_gap_c21_bound_inconsistency`, `test_data_gap_c31_bound_inconsistency`)
so the gaps stay visible:

1. `C_(2,1)`: Table 3 `3.1297495000000` vs Eq. 8 min `3.129751730201047`
   (Delta 2.24e-6).
2. `C_(3,1)`: Table 3 `3.1833333078762` vs Table-4-implied `3.1756140`
   (Delta 7.7e-3).

These are the bound (`C_(k1,k2)`) columns only; the `C^stable` / `T^stable`
columns used throughout this pass are internally consistent and are the
catalogue-grade numbers. The bound columns are NOT used as goldens here.

## 6. Honest caveats

- **PCR3BP only.** Everything here is the pure planar circular RTBP. No
  bicircular / solar-gravity / lunar-eccentricity / ephemeris refinement — those
  are the paper's stated future work, and any cross-fidelity persistence claim
  must wait for the 2026 journal (mining note §10).
- **`nu ~ 0`, not `nu = 0`.** The reproduced members sit close to (not exactly
  at) the published `nu = 0` midpoints because the corrector single-shoots the
  fixed published `C`; the exact midpoint is found by continuation, not by this
  single-shot. The STABLE verdict is robust (`|nu|` well below 1 in all five).
- **Derived states are not goldens.** Only Ross's printed `(mu, C, T,
  stability)` are EXPECTED; the `(x0, ydot0)` we recover are derived and carry
  provenance.
- **No catalogue writeback.** All rows/deltas above are proposals for #216.
