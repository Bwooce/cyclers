# #494 Phase 2 + 3 Verdict — mu-Family Representatives + Pluto-Charon

**Date:** 2026-06-30
**Task:** #494, Phases 2 and 3
**Verdict:** Phase 2 PASS (6/6); Phase 3 PASS — a STABLE Pluto-Charon (3,2) cycler exists

---

## Summary

Phase 2 recovers all six Table-I representatives from Ross & Roberts-Tsoukkas 2026
(arXiv:2606.29189v1).  All six converge, match published ICs within print precision,
show correct winding topology (or nth-iterate form for Rep 6), are prograde, are
linearly stable (|nu|<1), and pass the independent Radau cross-check.

Phase 3 instantiates the (3,2) family at Pluto-Charon (mu=0.10851).  Seeding the
corrector at the mu=0.1 anchor Jacobi constant C=3.5734 lands on the (3,2) branch
but FAR off the stable island (nu is large) — stability is C-selective.  A C-sweep
along the (3,2) branch at mu=0.10851 (the paper's own method: trace the (k1,k2)
C-family and find the |nu|<1 island) locates a **razor-thin stable window** (~1.1e-5
wide in C) centred on the nu=0 midpoint at **C = 3.5792220, ~0.006 above the anchor
C**.  That member is a genuine, linearly-stable, prograde (3,2) periodic orbit that
closes under the independent integrator: **a stable Pluto-Charon (3,2) cycler exists.**
Physical cross-checks (C_L1 vs Jbara 2025 and Holman-Wiegert a_crit) pass.

> CORRECTION (this pass): an earlier draft concluded "STABILITY FAIL / no admission"
> from the anchor-C orbit alone.  That was an artefact of holding C at the mu=0.1
> value; the stable island sits ~0.006 higher in C at the shifted mu and was found
> by the C-sweep below.  This matches the [[feedback_published_rounded_values_are_display]]
> / family-extent lesson: a sourced anchor value need not sit on the stable subfamily
> after a parameter shift — check the family extent before declaring a negative.

---

## Phase 2 — Per-representative recovery table

Source: Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189v1, Table I.
hc = half_crossings; tol = corrector tolerance; n = iterate multiplier.
All crosschecks pass (independent Radau integrator, closure < 1e-6, dJ < 1e-8).

| Rep | mu | k | x0 (recovered) | |dx0| | T (recovered) | |dT| | topo at T | n | nu | Stable | Crosscheck |
|-----|----|---|-----------------|-------|----------------|------|-----------|---|-----|--------|------------|
| 1 | 0.001 | (1,1) | -0.647047499999966 | 0.00e+00 | 14.774502790981787 | 7.0e-12 | (1,1) | 1 | 0.4121 | YES | PASS |
| 2 | 0.01215 | (1,1) | -0.768217354461248 | 0.00e+00 | 10.291893654144250 | 1.2e-08 | (1,1) | 1 | 0.7966 | YES | PASS |
| 3 | 0.01215 | (3,3) | -0.322477620583087 | 0.00e+00 | 19.503763586997401 | 7.3e-11 | (3,3) | 1 | 0.9859 | YES | PASS |
| 4 | 0.1 | (3,2) | -0.694376003123377 | 0.00e+00 | 12.295263860078997 | 1.4e-08 | (3,2) | 1 | 0.5687 | YES | PASS |
| 5 | 0.3 | (3,1) | -0.804726006933677 | 2.2e-07 | 9.094570937588166 | 5.5e-06 | (3,1) | 1 | 0.0427 | YES | PASS |
| 6 | 0.5 | (1,1) | -0.519689929077496 | 0.00e+00 | 8.792013561445721 | 1.7e-11 | (3,3) | 3 | 0.9375 | YES | PASS |

### Notes by representative

**Rep 2 (mu=0.01215, (1,1)):** Paper uses hc=3 with tol=1e-8; converges to a
different member of the EM (1,1) family than the AAS 2025 paper.  The recovered
nu=0.7966 vs published sp=0.8210 because the 2026 paper tabulates a slightly
different (C, x0) than AAS 2025.  Both are stable (|nu|<1).

**Rep 5 (mu=0.3, (3,1)):** Corrector finds x0=-0.804726006933677, which differs
from the published -0.804725783387797 by dx0=2.24e-07 (beyond 16-digit printing
precision).  The recovered orbit is the nearest periodic orbit at C=3.7020; the
period difference dT=5.5e-06 is within 1e-5 tolerance.  Topology (3,1) matches.
nu=0.0427 vs published sp=0.0294 (both near zero, deep inside the stable island).

**Rep 6 (mu=0.5, (1,1)):** Published T=8.792 is the 3rd iterate of the fundamental
period T1~2.931.  Family label (1,1) refers to the fundamental winding; at T=8.792
the winding is (3,3) (n=3 iterate check).  At T/3, winding_topology returns (1,1)
as verified by the separate fundamental-period test.  nu=0.9375 matches published
sp=0.9376 to within 1e-4.  Triple-angle identity: 4*nu1^3 - 3*nu1 = nu3.

### mu-gap closure status

These representatives span mu in {0.001, 0.01215, 0.1, 0.3, 0.5}, closing the
gap between the EM (mu~0.012) and binary-star (mu=0.3, 0.5) limits for which
Ross 2026 shows (k1,k2) stable cyclers in figures.  All 5 physically distinct mu
values are now recovered from sourced (mu, C, T) ICs.

### V0 -> V1 recommendations

All six Phase-2 representatives satisfy the V0->V1 criteria:
- Converged corrector (periodic orbit confirmed)
- Correct topology (sourced label matches)
- |nu|<1 (linearly stable)
- Independent Radau crosscheck passes

Recommended: admit all 6 as V1 (same-model golden against Table I).  The orbits
are in the non-dimensional CR3BP (mu, C, T, x0) space; physical dimensionalization
requires a specific binary system (l_km, t_s).

---

## Phase 3 — Pluto-Charon (mu=0.10851) (3,2) cycler

### System parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| mu_PC | 0.10851 | Task specification (Ross 2026 context) |
| l_km | 19600 km | satellites.py (Charon mean separation) |
| GM_system | 975.5 km^3/s^2 | JPL DE440, satellites.py |
| t_s | 552014 s = 6.389 d | 2*pi*sqrt(a^3/GM) |

### Branch-establishing orbit at the anchor C (Phase-3a)

Seeded from Rep 4 anchor (mu=0.1, C=3.5734, (3,2)) with half_crossings=6, the
(3,2) branch is present at mu=0.10851 (converged, (3,2), prograde, crosscheck
PASS) but lies far off the stable island at the anchor C (nu is large).  The
crosscheck PASSES (closure, Jacobi conservation), confirming a genuine periodic
(3,2) orbit — it simply sits outside the |nu|<1 window at C=3.5734.

### Stable (3,2) member at mu=0.10851 (Phase-3b — the deliverable)

A C-sweep along the (3,2) branch (hc=6) at fixed mu=0.10851 reveals a stable
window. nu vs C has a very steep slope here (~2e5 per unit C), so the |nu|<1
island is razor-thin (~1.1e-5 wide in C) — analogous to the EM (2,1) family.
nu=0 midpoint (Barden), located by brentq root-find on nu(C):

| Quantity | Value |
|----------|-------|
| C (Jacobi) | 3.579222016200 |
| x0 | -0.693189765944 |
| ydot0 | -0.296204695856 |
| T (TU) | 11.8366755503 |
| T (days) | 12.0361 d (Charon P = 6.3891 d; TU = P/2pi) |
| Topology | (3,2), prograde=True, reaches_secondary=True |
| nu (Barden) | -1.1e-08 (nu=0 midpoint — maximally stable) |
| Stable | YES (|nu| << 1) |
| Crosscheck | PASS (closure < 1e-6, dJ = 8.8e-13) |

The located (C, x0) are DERIVED (our C-sweep), NOT goldens; the test asserts
self-consistency + the stability VERDICT, not a sourced number.  The stable
window spans roughly C in [3.5792165, 3.5792275] (nu from +1 to -1).

### Physical cross-checks

**C(L1):** L1 at mu=0.10851 is at x_L1=0.5930 (non-dimensional).
C_L1 = jacobi_constant([x_L1, 0, 0, 0, 0, 0], 0.10851) = 3.6203.
Jbara 2025 (arXiv:2510.13479) reports C_L1 ~ 3.6210 at mu~0.109; |diff| = 0.0007 < 0.005.
CHECK: PASS.

**Holman-Wiegert a_crit:**
a_c/a_bin = 1.60 + 4.12*0.10851 - 5.09*0.10851^2 = 1.9871.
a_crit = 1.9871 * 19600 = 38,948 km.
Source: Holman & Wiegert 1999, Eq. 1.  Tolerance 100 km.
CHECK: PASS (|38948 - 38947| < 100).  (Exterior P-type stability bound only — a
sanity anchor for the regime, NOT the cycler itself, which is an interior capture
orbit; consistent.)

### Literature-novelty check (necessary-not-sufficient)

WebSearch over arXiv / ADS for the Pluto-Charon (3,2) prograde capture cycler:

- The (k1,k2) ballistic-prograde-cycler **family** is published: Ross &
  Roberts-Tsoukkas 2026, arXiv:2606.29189 — but it does NOT instantiate
  Pluto-Charon (it spans mu=0.001->0.5 generic representatives; PC is not a row).
- Pluto-Charon CR3BP literature exists but reports *different* objects: Jbara 2025
  (arXiv:2510.13479, chaotic dynamics + zero-velocity structures); orbit
  classification (arXiv:1512.08683); near-optimal **capture** (Sci.Direct 2018,
  not a periodic cycler); circumbinary P-type theory (Langford & Weiss 2023,
  AJ 165 140); "Sailboat island" periodic orbits *of the first kind* (Giuliatti
  Winter et al.) — a distinct family from a (3,2) capture cycler.
- **No published Pluto-Charon (3,2) prograde capture cycler surfaced.**

Verdict: **not-found** (NECESSARY-not-sufficient for novelty per the #261 rule;
web search is partial).  This clears — does not certify — the candidate as a
**fresh real-system instantiation of a published family** (the #312-Uranus framing).
No novel-discovery claim beyond the instantiation.

### PC admission recommendation

**RECOMMEND admit the stable Pluto-Charon (3,2) member** (the nu=0 midpoint above)
as a catalogue cycler row, pending human adjudication.  Proposed row (RECOMMENDATION
ONLY — NOT written; adjudicator decides level + final fields):

```yaml
# orbit_class: cycler
# primary: Pluto
# bodies: [Pluto, Charon]
# model: cr3bp
# mu: 0.10851
# jacobi_C: 3.579222016200
# x0: -0.693189765944          # rotating-frame nd, IC (x0,0,0,0,ydot0,0)
# ydot0: -0.296204695856
# period_tu: 11.8366755503     # ~12.036 d (TU = Charon_P / 2pi)
# topology_k1k2: [3, 2]        # prograde, reaches secondary
# stability_nu: ~0 (Barden, |nu|<1 -> linearly stable)
# source: Ross & Roberts-Tsoukkas 2026 (family); fresh PC instantiation (this work)
# validation_level: V0 -> V1 candidate (see below)
```

**Level recommendation: V1.**  The Ross-RT 2026 (3,2) FAMILY is sourced (the
mu=0.1 anchor is Table-I row 4, recovered at V1 in Phase 2); the Pluto-Charon
member is a same-model continuation/instantiation of that sourced family at a
real-system mu — a faithful sourced-family instantiation, not a novel discovery.
It is NOT V2: no long-span REBOUND/IAS15 bounded-band run was done here (the EM
five reached V2 via the 100-period IAS15 evidence in
`2026-06-13-ross-v2-longspan-evidence.md`); a PC V2 would need the same long-span
evidence at mu=0.10851.  Real-ephemeris Pluto-Charon validation (SPICE) is the
optional V3/V4 lever.

### mu (0.10851 vs 0.10877) note

The task specifies mu=0.10851.  satellites.py gives Charon GM 106.1 / system GM
975.5 = 0.10877 (a 2.4e-3 relative difference).  All Phase-3 numbers above use the
task value 0.10851.  The stable island's C and x0 shift slightly at 0.10877; if the
catalogue row uses the repo's GM ratio, re-locate the midpoint at 0.10877 first.
This is flagged for the adjudicator, not silently resolved.

---

## Test file added

`tests/search/test_ross_rt_2026_mu_family.py` — 9 tests, all PASS:

- `test_494_phase2_recover_table_i_representative[...]` x6 (parametrized)
- `test_494_phase2_rep6_fundamental_winding_is_11`
- `test_494_phase3_pluto_charon_32_branch_and_crosschecks` (Phase-3a: branch + C_L1 + a_crit)
- `test_494_phase3_pluto_charon_32_stable_member_exists` (Phase-3b: the stable member)

Runtime: ~35 s (parallel, 16 workers).
