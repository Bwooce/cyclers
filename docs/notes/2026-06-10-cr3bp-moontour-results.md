# CR3BP Saturnian Midsize-Moon Discovery Results (degeneracy-gated, v2)

**Run timestamp:** 2026-06-10T04:18:57Z
**Script:** `scripts/cr3bp_moontour_run.py`
**Method:** `cr3bp-lyapunov-corrector-v2`  (git `f166e1c`)

**Status:** NO catalogue writeback.  SILVER entries in `data/cr3bp_silver.jsonl`
(review-gated).  EMPTY notes (if any) in `data/empty_regions.jsonl`.

## Why v2

The v1 run (`cr3bp-lyapunov-corrector-v1`, 11 entries, never committed) was
contaminated: the min-norm full-state corrector "converged" onto the L1/L2/L4/L5
libration points themselves (rotating-frame speed ~1e-12, closed for any period),
onto a period-collapse solution (T=4.5e-5), and onto the same equilibrium twice
with different reported periods.  That output was deleted and regenerated from
scratch under the v2 degeneracy gate:

1. **Equilibrium rejection** -- max |v| over the propagated period >= 1e-6 nondim
   AND position amplitude max|r(t)-r(0)| >= 1e-6 nondim.
2. **Period floor** -- period >= 0.1 nondim.
3. **Dedup** -- state0 positions within 1e-9 are one orbit (best residual kept).
4. **Independent cross-check** -- `crosscheck_periodic` (Radau, vs the corrector's
   DOP853) must re-close within 1e-8 with Jacobi drift < 1e-8; dJ recorded.

v2 also fixes the seeding: the v1 draft's linear frequency formula was wrong
(it substituted `1+2c2` where the standard formula takes `c2`, with a sign slip),
and its constant-Jacobi seeding placed initial velocities off the Lyapunov family.
v2 seeds the linearised centre eigensolution at Newton-solved collinear points:
`x0 = x_L + Ax`, `vy0 = -k nu Ax`, `T_guess = 2 pi / nu`.

## Run report (verbatim)

```text
cr3bp_moontour_run  2026-06-10T04:18:57.954431+00:00  method=cr3bp-lyapunov-corrector-v2  git=f166e1c  elapsed=6.2s
==============================================================================
pair: Saturn/Mimas  mu=6.599e-08
  seeds tried        : 10  (L1+L2 x amplitude fracs [0.0005, 0.001, 0.002, 0.003, 0.005])
  converged          : 5
  rejected degenerate: 0  [equilibrium=0 period_floor=0 duplicate=0 crosscheck=0]
  crosscheck passed  : 5 of 5 checked
  SILVER written     : 5
    L2 Ax/gamma=0.001  J=3.0000703462  T=3.03902654  resid=1.12e-12  max|v|=6.700e-06  amp=4.146e-06  dJ=0.00e+00
    L1 Ax/gamma=0.002  J=3.0000704300  T=3.02704118  resid=3.77e-12  max|v|=7.345e-05  amp=5.984e-05  dJ=0.00e+00
    L2 Ax/gamma=0.0005  J=3.0000703462  T=3.03902658  resid=1.70e-11  max|v|=7.124e-06  amp=3.733e-06  dJ=0.00e+00
    L1 Ax/gamma=0.0005  J=3.0000704342  T=3.02700485  resid=6.26e-11  max|v|=1.072e-05  amp=5.438e-06  dJ=0.00e+00
    L1 Ax/gamma=0.001  J=3.0000704338  T=3.02700785  resid=7.89e-11  max|v|=2.347e-05  amp=1.219e-05  dJ=4.44e-16

pair: Saturn/Enceladus  mu=1.901e-07
  seeds tried        : 10  (L1+L2 x amplitude fracs [0.0005, 0.001, 0.002, 0.003, 0.005])
  converged          : 4
  rejected degenerate: 0  [equilibrium=0 period_floor=0 duplicate=0 crosscheck=0]
  crosscheck passed  : 4 of 4 checked
  SILVER written     : 4
    L2 Ax/gamma=0.0005  J=3.0001421651  T=3.04156452  resid=9.25e-13  max|v|=1.040e-05  amp=5.423e-06  dJ=0.00e+00
    L1 Ax/gamma=0.001  J=3.0001424177  T=3.02446170  resid=2.21e-12  max|v|=3.431e-05  amp=1.768e-05  dJ=4.44e-16
    L2 Ax/gamma=0.001  J=3.0001421651  T=3.04156446  resid=7.34e-12  max|v|=9.366e-06  amp=5.943e-06  dJ=4.44e-16
    L1 Ax/gamma=0.0005  J=3.0001424184  T=3.02445848  resid=1.38e-11  max|v|=1.524e-05  amp=7.722e-06  dJ=4.44e-16

pair: Saturn/Tethys  mu=1.086e-06
  seeds tried        : 10  (L1+L2 x amplitude fracs [0.0005, 0.001, 0.002, 0.003, 0.005])
  converged          : 5
  rejected degenerate: 0  [equilibrium=0 period_floor=0 duplicate=0 crosscheck=0]
  crosscheck passed  : 5 of 5 checked
  SILVER written     : 5
    L1 Ax/gamma=0.002  J=3.0004536092  T=3.01774702  resid=6.77e-14  max|v|=2.068e-04  amp=1.186e-04  dJ=4.44e-16
    L2 Ax/gamma=0.0005  J=3.0004521943  T=3.04828531  resid=6.24e-13  max|v|=1.961e-05  amp=1.014e-05  dJ=4.44e-16
    L1 Ax/gamma=0.0005  J=3.0004536427  T=3.01770182  resid=1.53e-12  max|v|=2.643e-05  amp=1.344e-05  dJ=0.00e+00
    L2 Ax/gamma=0.001  J=3.0004521941  T=3.04828561  resid=3.06e-12  max|v|=2.586e-05  amp=1.360e-05  dJ=4.44e-16
    L1 Ax/gamma=0.001  J=3.0004536403  T=3.01770509  resid=3.70e-11  max|v|=6.117e-05  amp=3.139e-05  dJ=0.00e+00

NO catalogue writeback performed.
```

## Per-pair summary

| Pair | Seeds | Converged | Rejected (equilibrium/period/dup/xcheck) | SILVER |
|---|---|---|---|---|
| Saturn/Mimas | 10 | 5 | 0/0/0/0 | 5 |
| Saturn/Enceladus | 10 | 4 | 0/0/0/0 | 4 |
| Saturn/Tethys | 10 | 5 | 0/0/0/0 | 5 |

**Total SILVER candidates:** 14

## Honest caveat

Small collinear-point (L1/L2) planar Lyapunov families are mathematically
guaranteed to exist in any CR3BP (Lyapunov centre theorem) and are NOT novel
discoveries in the literature sense -- such families are tabulated for many
systems.  Their value here is exercising and validating the seed -> STM-corrector
-> degeneracy-gate -> independent-crosscheck pipeline on new (Saturn, midsize-moon)
mass ratios.  They route to review (SILVER) regardless; promotion past SILVER
would require a sourced anchor.  **NO writeback to `data/catalogue.yaml`.**

All survivors are SMALL-amplitude members (Ax ~ 1e-6..1e-5 nondim, i.e. of order
a kilometre at these moons): the full-period single-shooting corrector's
convergence basin on these strongly unstable orbits is tiny, and larger-amplitude
seeds either diverge or are captured by the L4/L5 equilibria (gate-rejected).
Extending the families to useful amplitudes needs a symmetry-exploiting
half-period or multiple-shooting corrector -- a method-capability gap, recorded
here so a future sweep can supersede this one.

## Addendum (run forensics, recorded manually)

**Why every rejected-degenerate count is 0 in this run:** with basin-aware seed
amplitudes, the out-of-basin seeds (fracs 0.002-0.005 at most points) *diverged*
(`converged=False`, counted as non-converged) instead of converging degenerately.
The gate itself was verified directly against the archived v1 contamination:
feeding the v1 L1-equilibrium entry (state0 x=0.99720, |v|~1e-12, "T"=1.5539)
through `degeneracy_gate` returns `("equilibrium", max|v|=1.0e-11, amp=2.7e-12)`,
and the v1 period-collapse entry (T=4.49e-5 at the L4 point) returns
`("period_floor", ...)`.  Both v1 failure classes are rejected.  At the larger
(now unseeded) amplitudes the corrector demonstrably converges onto L4/L5
(e.g. Mimas L1 seed Ax/gamma=0.05 -> x0=0.4999..., |v0|=1.7e-15), which the
equilibrium gate also rejects -- observed during basin measurement.

**Did the v1 Tethys pass ever run?** No completed pass: the v1 (uncommitted,
now-deleted) `data/cr3bp_silver.jsonl` ends with Enceladus entries
(t_added 03:45:09Z); there is no Tethys SILVER entry, no Tethys EMPTY record in
`data/empty_regions.jsonl`, and no results note was written (the v1 script wrote
its note only at the end of `main`).  The v1 run was interrupted at/before the
Tethys pair and never reached its reporting stage.  In v2 all three pairs ran to
completion (per-pair tallies above).

**v1 contamination tally (for the record, 11 entries, all discarded):**
7 libration-point equilibrium entries (L1/L2 at x~0.996-1.004, L4/L5 at
(0.5, +-0.866); rotating-frame speeds 1e-12..1e-15), two of which were the SAME
Enceladus L1 equilibrium under different reported "periods"; 1 period collapse
(T=4.5e-5); and 3 genuine small Lyapunov orbits.  None carried an independent
cross-check, so all 11 were deleted and the survivors regenerated under v2.
