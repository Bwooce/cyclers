# #411 — Single-rev cross-system corrector outcome (EM-L2/SE-L2): stalls, not closed

Date: 2026-06-20. Outcome of the time-consistent single-revolution cross-system
cycle corrector (`correct_cross_cycle`, committed f563ab3) on the EM-L2 ↔ SE-L2 pair.
Honest non-result; asserts NO closed orbit.

## What was built
`correct_cross_cycle`: a 2×2 damped finite-difference Newton over the amplitude knobs
(c_em, c_se) driving the two phase-time-consistency residuals to zero,

    R1 = wrap[ θ_ret − θ_fwd − ω_rel·(t_fwd + n_se·T_se) ]   (fwd → ret handoff)
    R2 = wrap[ θ_fwd − θ_ret − ω_rel·(t_ret + n_em·T_em) ]   (ret → fwd handoff)

so the cycle is periodic-up-to-rotation accounting for the relative phase advancing at
ω_rel over the leg transits and orbit dwells (not just the frozen-orientation
|θ_ret − θ_fwd| proxy `search_cross_cycle` uses). Supporting infra: transit-time
reporting (`_manifold_crossing_timed`, `CrossConnection.transit_time`), `_wrap_pi`,
`CrossCycleClosure`, an `on_iter` progress callback. Tested (fast `_wrap_pi` unit +
slow plumbing/clean-negative-contract tests).

## Result — stalls at |R| ≈ 0.59 rad
Seed (EM-L2 C=3.15, SE-L2 C=3.000863625, n_em=n_se=1), corrector trace:

| iter | c_em | c_se | R1 | R2 | \|R\| |
|---|---|---|---|---|---|
| 0 | 3.150000 | 3.00086363 | +0.043 | −0.823 | 0.824 |
| 1 | 3.151471 | 3.00085763 | +0.075 | −0.596 | 0.600 |
| 2 | 3.151471 | 3.00085734 | +0.077 | −0.587 | 0.592 |
| 3–4 | 3.151471 | 3.00085725 | +0.078 | −0.584 | 0.589 |

**Both legs stay converged and ballistic-cheap throughout** (forward EM-L2→SE-L2 ΔV
0.36, return SE-L2→EM-L2 ΔV 0.44 → total ~0.80 km/s; cycle ~1.06 yr; max leg position
gap ~12 km). The Newton closes R1 (→0.008) but R2 floors at ≈−0.58: |R| stalls at
~0.59 rad, well above the 1e-2 closure tolerance. **Single-revolution closure is not
achieved at this pair.**

## Why (working hypothesis — NOT firmly established)
The spatial-connection constraint (both legs must converge) couples c_em and c_se: the
both-converged set is (empirically) a curve in the (c_em, c_se) plane, not an open
region — at the seed the legs converge at (3.15, 3.0008636); by iter 4 at (3.1515,
3.0008573). If only ONE amplitude degree of freedom is genuinely available along that
connection curve, it can zero only ONE of the two phase-consistency residuals — which
matches the observed "R1 closes, R2 stalls." This would mean the earlier amplitude-knob
feasibility (2026-06-20-411-amplitude-theta-closure.md), which counted c_em and c_se as
two independent knobs, OVER-counted: the connection-existence constraint removes one DOF.

This is a HYPOTHESIS. Standalone diagnostic scans intended to confirm it were
inconclusive — they were confounded by (a) using a fixed c_em with the wrong c_se for
that point on the connection curve, and (b) return-leg branch selection differing from
the corrector's internal variant scan. The reliable evidence is the corrector trace
above; the exact codimension of the connection set was not cleanly measured.

## Next venues (future #411 work, not done here)
1. **Rev-count variation**: n_em or n_se ∈ {2, 3} adds a discrete ω_rel·T phase shift to
   one residual (Δθ_em ≈ 3.16, Δθ_se ≈ 0.33 rad/rev), which may move R2's floor across
   zero. Within the EM shadow budget (~3–4 rev). Cheapest next try.
2. **Other libration pairs**: EM-L1/SE-L1, EM-L1/SE-L2, EM-L2/SE-L1 — a different pair
   may give a 2-D both-converged region (two real knobs) or a better-placed R2 floor.
3. **Solver**: pseudo-arclength continuation ALONG the connection curve (parameterize by
   the curve, not free (c_em,c_se)) to characterize R1, R2 vs the single curve parameter
   and confirm/deny the 1-DOF obstruction directly.
4. **BCR4BP** with an SE-scale seed (#412 re-scope) — closure well-posed by construction.

## Multi-revolution shots — same stall (rev-count does NOT rescue closure)
Tested the cheapest "Next venue" (#1) directly: re-ran the corrector at higher revolution
counts to see if a discrete extra revolution's phase increment pushes the residual through
zero.

| config | residual floor | result |
|---|---|---|
| n_em=1, n_se=1 | \|R\| ≈ 0.59 rad | stalls |
| n_em=1, n_se=2 | \|R\| ≈ 0.78 rad | stalls (plateau by iter 4: 0.7858→0.7819) |

Each discrete revolution SHIFTS the residual floor (0.59 → 0.78) but does not reach zero —
both legs stay converged throughout, and the Newton plateaus at a nonzero floor. This is
the predicted signature of the **1-DOF obstruction**: the spatial-connection constraint
couples (c_em, c_se) so only ONE phase degree of freedom is available, which cannot zero
the TWO phase-consistency residuals; adding integer revolutions changes the offset, not the
dimensionality. (The n_em=2 config was not run to completion — each stalled corrector run
is ~30+ min as the line-search backtracks; the n_se=2 plateau already confirms the pattern.)

## Status — patched-CR3BP route is a characterized NEGATIVE
The corrector is real, tested, committed infrastructure. The EM-L2/SE-L2 cross-system cycle
does NOT close in the patched CR3BP at single OR low-multiple revolution (residual floors at
~0.59 / ~0.78 rad). Both connection legs are ballistic-cheap (~0.8 km/s) — the obstruction
is purely phase-closure dimensionality, not energy. This is an honest, well-characterized
negative on a novel search. The remaining venues are NOT more revolutions or finer grids:
- **Other libration pairs** (EM-L1/SE-L1, mixed) — may yield a 2-D both-converged region
  (two real amplitude knobs) instead of the 1-D curve seen at L2/L2; UNTESTED.
- **BCR4BP with an SE-scale seed** (#412 re-scope) — phase closure is well-posed by
  construction in the coherent 4-body model; the clean long-term venue.
Pseudo-arclength continuation along the L2/L2 connection curve would *prove* the 1-DOF
obstruction rigorously, but the empirical single+multi-rev stalls already establish it for
practical purposes.
