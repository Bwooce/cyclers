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

---

## Progress update (2026-06-14, same day)

**Task 1 DONE + the fold-turn capability is PROVEN.**
`src/cyclerfinder/search/cr3bp_jacobi_arclength.py` + tests built and committed
(`e34edb2`); pseudo-arclength continuation in (x0, C) at fixed mu, 5 tests
(reproduce-before-trust: self-consistency, period continuity, unit-null tangent,
fixed-C landing).

Controller-side validation against the real (1,1) saddle-center fold:
**the continuation turns the fold** — seeded from the x0=-0.768 / C=3.1294 /
T=55.995 d member, it walked up to the fold apex **C=3.15117** and back down to
C=3.10060 (75 members, stop=step_underflow). `continue_family` (natural-parameter)
could not do this; the new arclength engine does.

**C11a not yet recovered — two concrete follow-ups identified:**
1. The x0=-0.768, T=55.995 d member is **unstable (nu=+7.5e4)** — it is NOT the
   stable C11b branch despite matching the sourced C11b period. The genuine stable
   C11b (|nu|<1) seed must be identified first (scan x0 at C=3.1294 for the
   |nu|<1 (1,1) member), then continued.
2. The return branch (post-fold, C descending through 3.1294) exists but the
   coarse member scan (tol 5e-4) did not extract its C=3.1294 member — use
   `land_at_jacobi(c_target=3.1294)` from a return-branch member to pull the
   candidate C11a (~42.140 d) exactly, then confirm period + Barden + Radau.

Next session: (1) identify stable-branch seeds for (1,1) and (2,1); (2) land both
fold branches exactly at C=3.1294; (3) confirm C11a/C21/C32 periods vs sourced;
(4) add resonant U-branches via JPL unstable seeds; (5) re-run the C32-dominance
gate with the enlarged node set.

---

## FINAL (2026-06-14): C11a/C21 NOT recovered — blocked on unpublished ICs

Ran the fold-turn extraction: continued the (1,1) seed in C past the apex
(C in [3.119, 3.183], 99 members), then landed every branch straddling C=3.1294.
Members recovered at the common energy:
- T=54.802 d (nu=-1, marginal),  T=55.995 d (nu=+7.5e4, unstable; the C11b period),
- **T=79.499 d (nu=-0.344, STABLE)** — the C32 member (matches #247's ~79.50 d),
  now recovered as a genuine stable member.
- **NO 42.140 d member (C11a)** and **no 84.533 d (C21)** anywhere on the
  fold-connected branch from this seed.

**Verdict:** the fold-turning engine WORKS (built, tested, proven to turn the
saddle-center fold; commit e34edb2) and adds the stable C32 node, but **C11a and
C21 are not reachable from the available seeds** — they live on a different
branch/family whose seed Braik-Ross does not publish (only the period is given,
which does not pin an orbit in a multi-branch family). This is the SAME
missing-published-IC blocker as the RT-Ross binary-star case (#255) and the
paywalled acquisitions (#116): re-derivation from a single period is
under-determined. The scored network grows 6 -> 7 (adds stable C32) but C11a/C21
remain missing, so the C32-dominance gate stays a faithful PARTIAL/negative.

**Disposition:** engine deliverable DONE; full member recovery BLOCKED on
unpublished state vectors (re-file as source-blocked, like #116). Not claimed
recovered.
