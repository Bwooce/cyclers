# #249 — recover the off-stable / unstable reachable-set members (plan + scoping)

**Date:** 2026-06-14
**Status:** SCOPED. Not a seeding gap — a fold-turning continuation build.

## The problem (confirmed)
At the Braik-Ross common energy **C_J = 3.1294**, #247 recovered 6 of the 13
network nodes (LL1/LL2/DPO/R21-S exact, C11b exact, C32 ~1.1%). The gap:
- **C11a (42.140 d)** and **C21 (84.533 d)** — not recovered at all.
- **C32 (78.613 d)** — only ~1.1% (off-stable perpendicular-crossing branch).
- **R21-U / R31-S/U / R52-S/U** resonant branches — not recovered at the common C.

With only 6 nodes the C32-dominance gate can't be tested, so the scorer stays
GATED (faithful negative).

## Root cause (diagnosed this task)
C11a/C11b are the two branches of the (1,1) family meeting at a **saddle-center
fold** in C (Ross-RT Fig. 4 shows the (3,2) version). Two tools were tried and
both structurally cannot cross the fold:
1. the 1-DOF perpendicular-x-crossing symmetric corrector lands on the stable
   branch (C11b) for any seed near it; #247's seed/half-crossing grid scan
   already established C11a/C21 are not within 1.5 d anywhere on it.
2. `cr3bp_continuation.continue_family` is **natural-parameter** Jacobi
   continuation (step C, re-solve x0). It diverges exactly at the fold
   (dx0/dC -> infinity), so it cannot walk C11b -> C11a.

## The build (tractable — machinery mostly exists)
**Pseudo-arclength continuation in Jacobi C at fixed mu.** `mu_continuation.py`
already implements fold-turning pseudo-arclength in z = (x0, C, mu): predict along
the tangent (null vector of dr/dz), correct back onto r = xdot(t_half) = 0 with the
arclength constraint, adaptive step, fold-capable. It currently arclengths in mu.
Adapt it to **freeze mu and arclength in (x0, C)** — the identical scheme, one
fewer free variable — to follow the (1,1) family curve from the recoverable C11b,
**turn the saddle-center fold**, and reach the C11a branch; then a fixed-C landing
projection (the `_land_at_mu` analogue at fixed C) drops exactly onto C = 3.1294.

### Tasks
1. `continue_in_jacobi_pseudoarclength(seed, *, mu, half_crossings, ydot0_sign,
   c_target, ...)` in a new `cr3bp_jacobi_arclength.py` (or extend
   `mu_continuation`): reuse `_residual_jac` / `_tangent` / `_correct` with the
   tangent/step taken in (x0, C) at fixed mu; `_land_at_c` projects onto C_target.
   TEST: reproduce C11b from itself (zero-length), and recover a *known* JPL
   family member across a fold as the reproduce-before-trust gate.
2. Recovery driver: seed from C11b (and the C21 stable member if recoverable, else
   from a JPL (2,1)-adjacent seed), arclength to the fold, turn it, land at
   C=3.1294, confirm period vs SOURCED_PERIODS_DAYS (C11a 42.140, C21 84.533, C32
   78.613) within tol. Resonant U-branches: seed from the JPL DB unstable branch
   (`_jpl_seed_near_cj(..., stable=False)`), which is already available.
3. Re-run the C_J=3.1294 gate with the enlarged node set; re-test C32-dominance
   (the `xfail` in `test_reachable_network_gate.py`). Flip xfail->pass ONLY if it
   genuinely reproduces; else keep the faithful negative and record which nodes
   recovered.
4. Update the #247 results note + the scorer's gated/ungated status.

## Honest scope
This is a real numerical build (~1 focused session). It is **discovery
infrastructure** (the prioritizer that says WHERE to seed), not a discovery
itself. Highest leverage given three consecutive blind-search negatives
(repeated-moon x3) argue we are searching without a map. Recommend executing as a
focused plan (subagent-driven or a fresh session) rather than rushed — the
foreground corrector scans are slow (~0.5 s/seed) and timed out repeatedly when
attempted inline.
