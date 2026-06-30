# #496 — bounded-least_squares cross-cycle solver: capability shipped, #411 stall NOT broken (sharper blocker)

**Date:** 2026-06-30. **Verdict:** capability delivered + a CHARACTERIZED NEGATIVE that re-locates
the binding #411 blocker. The cross-system cycle is NOT closed; no claim that it is.

## What shipped (commit 825687e)
`correct_cross_cycle` gains `solver="bounded_ls"` (+ `c_em_bounds`/`c_se_bounds`): a bounded
trust-region `scipy.optimize.least_squares` on the two phase residuals (R1, R2), with c_se
**clamped below the Canalias bifurcation** (3.000863625) so the corrector cannot take the wild step
off the SE-L2 family that breaks the damped FD-Newton at ~0.59 rad. Infeasible probes return the
max wrapped residual (π, π) so trf steers back. This is the Braik-Ross 2025 / Ross-Scheeres
blueprint (bound + feasibility-first), digested 2026-06-30
(`2026-06-30-digest-braik-ross-2025-heteroclinic-L1-L3.md`). Static-clean (ruff + mypy). Default
`solver="newton"` unchanged. Driver: `scripts/close_411_bounded_ls.py`.

## Verification (3 runs) — the binding blocker is UPSTREAM of the c_se stall
| seed (c_em, c_se) | outcome |
|---|---|
| 3.150, 3.00050 | seed **node** build fails — SE-L2 family is patchy: nodes build at c_se∈{3.0000, 3.0002, 3.0006, 3.00086} but FAIL at 3.0004/3.0005 (dead gaps) and error (neg. radicand) above 3.0010 |
| 3.150, 3.00060 | nodes build, but **seed legs do not converge** (manifolds don't co-reach the patch to 100 km) |
| 3.150, 3.00086 (Canalias) | nodes build; **legs still don't converge** — best position gap **2.17e5 km** (≫ 100 km tol), patch ΔV 0.86 km/s |

So at every accessible seed the two cross-system connection **legs do not BOTH close spatially**, and
`correct_cross_cycle`'s seed-convergence gate (requires both legs converged before the solver runs)
**blocks `bounded_ls` from ever iterating**. The c_se runaway-stall the solver fixes is real but
**downstream** of a leg-convergence problem that bites first at these seeds. Consistent with #405
("near-ballistic forward leg, hard return") and the #411 note's open-sufficiency caveat.

## The actual next step (the real two-phase formulation)
My `bounded_ls` only solves the θ-residuals AFTER the legs converge — that is NOT feasibility-first.
The genuine fix (true Phase-1) is to fold the **leg POSITION GAPS into the bounded_ls residual** and
**relax the seed gate**, so the solver drives leg closure + θ jointly over (c_em, c_se) [and possibly
the leg phases τ_u, τ_s] from an infeasible start — exactly the Braik-Ross multi-crossing +
closest-pair + feasibility-first pattern. AND seed at a leg-closing c_em (the #405 working region
c_em∈(3.15,3.184) needs a scan; 3.150 is not a forward-leg-closing point here). Both are real work,
not a re-run.

## Standing
#411 stall is **NOT broken**; the bounded_ls capability is shipped and sound (it removes the
c_se-off-family failure mode) but cannot be exercised until a both-legs-converging seed exists. Net:
a usable solver primitive + a sharper, honest re-statement of the binding obstacle (leg convergence +
SE-L2 family pathology near the only viable c_se). Tracked forward under #496 (feasibility-first
residual + leg-closing seed scan).
