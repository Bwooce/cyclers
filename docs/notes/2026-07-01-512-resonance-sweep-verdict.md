# #512 — (n_em, n_se) integer resonance sweep for the #411 cross-system cycle: verdict

**Date:** 2026-07-01.
**Status:** DRAFT — pending final solve results.

---

## 1. Background

#496 (docs/notes/2026-07-01-496-feasibility-first-verdict.md) confirmed a phase-closure
wall for the #411 planar cross-system (Sun-Earth <-> Earth-Moon) cycle at n_em=1, n_se=1:
both legs converge geometrically (EM-L2 -> SE-L2 -> EM-L2) at the seed c_em=3.150,
c_se=3.00086, but the phase residual R2 = -0.72 rad cannot be driven to zero within the
accessible (c_em, c_se) parameter space (EM-L2 family upper bound ~3.152-3.153; SE-L2
manifold convergence threshold ~3.000854, just below the Canalias bifurcation).

#496's next-steps section observed that incrementing n_em or n_se (extra revolutions on
the EM-L2 / SE-L2 orbit before executing the connecting leg) shifts R2 / R1 by a fixed
phase offset per step, and proposed a systematic (n_em, n_se) sweep as the natural next
move. This note reports that sweep.

---

## 2. Method — the sweep decouples exactly, no re-solving needed to rank

`correct_cross_cycle`'s phase residuals (module Task-5 note, `cross_system_cycle.py`
line ~940) are:

```
R1 = wrap[ theta_ret - theta_fwd - omega_rel*(t_fwd + n_se*T_se) ]
R2 = wrap[ theta_fwd - theta_ret - omega_rel*(t_ret + n_em*T_em) ]
```

`n_em`, `n_se` enter ONLY the phase bookkeeping (`_resid`), not the leg construction
(`_solve`): `theta_fwd`, `theta_ret`, `t_fwd`, `t_ret`, `T_em`, `T_se` are unchanged by
`n_em`/`n_se` at a fixed `(c_em, c_se)`. Since `wrap[x - c] = wrap[wrap(x) - c]` for any
`x, c` (they differ by an integer multiple of 2*pi), the sweep decouples exactly:

```
R1(n_se) = wrap[ R1(1) - omega_rel*(n_se-1)*T_se ]   (independent of n_em)
R2(n_em) = wrap[ R2(1) - omega_rel*(n_em-1)*T_em ]   (independent of n_se)
```

This is EXACT, not an approximation — it lets us rank all `(n_em, n_se)` pairs from a
SINGLE cheap seed evaluation (one `_solve` call, ~40s) before spending any solver
compute. `scripts/sweep_411_resonance.py` implements this: Step 1 prints the analytic
table (fast); Step 2 (`--solve`) re-runs `bounded_ls` (the #496 fix) at the top-ranked
candidates to check whether the amplitude corrector can actually close the residual to
zero after re-optimizing `(c_em, c_se)`.

---

## 3. Step 1 — analytic seed baseline and wrap-prediction table

Seed: c_em=3.150, c_se=3.00086 (EM-L2, SE-L2), `return_scan_n=8, return_scan_n_tau=3`
(the #496 scan-resolution fix).

```
R1(n_se=1) = +0.08579 rad     R2(n_em=1) = -0.70760 rad
omega_rel = 2.466216e-06 rad/s
T_em = 1,283,364.4 s (14.854 d)      T_se = 15,432,561.0 s (178.618 d)
R1 step = -0.36092 rad per n_se increment
R2 step = +3.11813 rad per n_em increment (close to pi -> odd/even n_em alternate sign)
```

(These freshly-computed R1/R2 at the seed differ slightly from #496's reported
R1=+0.0859, R2=-0.7246 — 0.0858 vs 0.0859 agrees to 4 digits; R2 differs by ~0.017 rad,
attributable to which of the 4 return-leg branch variants the (re-run) seed solve
selects as best. Both values place the seed in the same ~0.7-0.73 rad |R| regime; the
#512 sweep uses its own freshly-computed, self-consistent seed throughout.)

### 3.1 Practical grid: n_em = 1..8, n_se = 1..4 (32 pairs)

Full ranked table (best 10 shown; see `runlogs/sweep_411_resonance_analytic.log` for
all 32):

| rank | n_em | n_se | R1      | R2      | \|R\| (rad) |
|------|------|------|---------|---------|-------------|
| 1    | 1    | 1    | +0.0858 | -0.7076 | 0.7128      |
| 2    | 1    | 2    | -0.2751 | -0.7076 | 0.7592      |
| 3    | 3    | 1    | +0.0858 | -0.7545 | 0.7594      |
| 4    | 3    | 2    | -0.2751 | -0.7545 | 0.8031      |
| 5    | 5    | 1    | +0.0858 | -0.8014 | 0.8060      |
| 6    | 5    | 2    | -0.2751 | -0.8014 | 0.8474      |
| 7    | 7    | 1    | +0.0858 | -0.8484 | 0.8527      |
| 8    | 7    | 2    | -0.2751 | -0.8484 | 0.8919      |
| 9    | 1    | 3    | -0.6361 | -0.7076 | 0.9515      |
| ...  |      |      |         |         |             |
| 32   | 2    | 4    | -0.9970 | +2.4105 | 2.6086      |

**Key finding: the n_em=1, n_se=1 seed is the analytic BEST of all 32 practical-grid
pairs, by a clear margin.** Every other pair in the swept range is strictly worse.
This is because:
- R1(n_se) walks monotonically away from its near-zero seed value (+0.086) as n_se
  increases over 1-4 (down to -0.997 at n_se=4) — the seed happens to sit almost
  exactly at a local R1 minimum within the practical n_se range.
- R2(n_em) is best (least negative) at n_em=1 among the ODD values (1,3,5,7 -> -0.71,
  -0.75, -0.80, -0.85, monotonically worse) and far worse for any EVEN n_em (+2.2 to
  +2.4 rad, since the per-step shift 3.118 rad is close to pi).

### 3.2 Wide diagnostic scan (informational, not solved)

Extending the analytic (free) scan to locate the true nearest-to-zero crossings for
each dimension independently:

| best n_se for R1≈0 | \|R1\| (rad) | SE-dwell (days) |
|---------------------|--------------|------------------|
| 36                   | 0.0198       | 6430 (17.6 yr)   |
| 1 (seed)             | 0.0858       | 178.6            |
| 19                   | 0.1277       | 3394 (9.3 yr)    |

| best n_em for R2≈0 | \|R2\| (rad) | EM-dwell (days) |
|---------------------|--------------|------------------|
| 104                  | 0.0174       | 1545 (4.2 yr)    |
| 106                  | 0.0295       | 1575 (4.3 yr)    |
| 102                  | 0.0643       | 1515 (4.1 yr)    |

Combining the best of each (n_em=104, n_se=36) would give a predicted joint
`|R| = hypot(0.0174, 0.0198) ≈ 0.026 rad` — closer to zero than the practical-grid
optimum, but at the cost of 104 EM-L2 revolutions (~4.2 years just for that dwell) and
36 SE-L2 revolutions (~17.6 years just for that dwell) before the connecting legs
execute. Total cycle period would be ~20+ years. This is NOT solved by this driver:
it is a mathematical curiosity documented for completeness, not a practical cycler
candidate, and even at 0.026 rad the residual is still ~2.6x the closure tolerance
(`theta_tol_rad=1e-2`) before any amplitude re-solving.

---

## 4. Step 2 — bounded_ls re-solve at the top practical-grid candidates

(1,1) is already fully characterized by #496 (`feasibility_ls`/`bounded_ls`,
best `|R| = 0.517-0.518 rad`, `closed=False`) — not re-run here (existing result
reused).

Top NEW candidates solved: (n_em, n_se) = (1,2), (3,1), (1,3), using the same #496
bounds (`c_em_bounds=(3.112, 3.152)`, `c_se_bounds=(3.00050, 3.00086)`) and seed
(c_em0=3.150, c_se0=3.00086).

<!-- RESULTS TABLE FILLED IN BELOW ONCE THE DRIVER COMPLETES -->

---

## 5. Verdict

<!-- FILLED IN BELOW -->

---

## 6. Code deliverables

- `scripts/sweep_411_resonance.py`: analytic wrap-prediction table (Step 1, exact,
  cheap) + `--solve` driver re-running `bounded_ls` at the top candidates (Step 2).
- `runlogs/sweep_411_resonance_analytic.log`: Step-1-only run.
- `runlogs/sweep_411_resonance.log`: full run including Step 2 solves.
