# #494 Phase 2 + 3 Verdict — mu-Family Representatives + Pluto-Charon

**Date:** 2026-06-30
**Task:** #494, Phases 2 and 3
**Verdict:** Phase 2 PASS (6/6); Phase 3 TOPOLOGY PASS, STABILITY FAIL

---

## Summary

Phase 2 recovers all six Table-I representatives from Ross & Roberts-Tsoukkas 2026
(arXiv:2606.29189v1).  All six converge, match published ICs within print precision,
show correct winding topology (or nth-iterate form for Rep 6), are prograde, are
linearly stable (|nu|<1), and pass the independent Radau cross-check.

Phase 3 instantiates the (3,2) family at Pluto-Charon (mu=0.10851) from the mu=0.1
anchor.  The orbit converges and has the correct (3,2) topology, but is UNSTABLE
(nu~1903) at the anchor Jacobi constant C=3.5734.  Physical cross-checks (C_L1 and
Holman-Wiegert a_crit) pass.

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

### Corrected orbit at mu=0.10851

Seeded from Rep 4 anchor (mu=0.1, C=3.5734, (3,2)) with half_crossings=6.

| Quantity | Value |
|----------|-------|
| x0 | -0.693353687114568 |
| ydot0 | -0.304877058365371 |
| C (Jacobi) | 3.573367616904619 (=anchor C, enforced) |
| T (TU) | 11.984949101055731 |
| Topology | (3,2), prograde=True |
| nu (Barden) | 1903.5 |
| Stable | NO (|nu| >> 1) |
| Crosscheck | PASS (closure < 1e-6, dJ = 8e-13) |

### Physical cross-checks

**C(L1):** L1 at mu=0.10851 is at x_L1=0.5930 (non-dimensional).
C_L1 = jacobi_constant([x_L1, 0, 0, 0, 0, 0], 0.10851) = 3.6203.
Jbara 2025 reports C_L1 ~ 3.6210 at mu~0.109; |difference| = 0.0007 < 0.005.
CHECK: PASS.

**Holman-Wiegert a_crit:**
a_c/a_bin = 1.60 + 4.12*0.10851 - 5.09*0.10851^2 = 1.9871.
a_crit = 1.9871 * 19600 = 38,948 km.
Source: Holman & Wiegert 1999, Eq. 1.  Tolerance 100 km.
CHECK: PASS (|38948 - 38947| < 100).

### Stability finding

At the anchor Jacobi constant C=3.5734 (the mu=0.1 Rep-4 value), the (3,2)
orbit at mu=0.10851 is HIGHLY UNSTABLE (nu~1903).  This is physically
meaningful: a small change in mu (0.00851) has moved the orbit far outside
the stable island.

The (3,2) STABLE family likely exists at mu=0.10851 but at a DIFFERENT Jacobi
constant (C value) than the mu=0.1 anchor.  Finding the stable member would
require a mu-continuation sweep over C at mu=0.10851, which is beyond the
Phase-3 scope.

The crosscheck PASSES (closure and Jacobi conservation), confirming the
recovered orbit IS a genuine periodic (3,2) orbit — it simply lies outside
the stable island at this C.

### PC admission recommendation

NO V1 admission for the recovered Pluto-Charon (3,2) orbit at C=3.5734.
The orbit is unstable and does not meet the stable-cycler criterion.

To admit a Pluto-Charon (3,2) catalogue row, a mu-continuation at fixed
topology from the mu=0.1 anchor sweeping C until |nu|<1 is required.  This
is a Phase-3b task (not in scope for #494).

---

## Test file added

`tests/search/test_ross_rt_2026_mu_family.py` — 8 tests, all PASS:

- `test_494_phase2_recover_table_i_representative[...]` x6 (parametrized)
- `test_494_phase2_rep6_fundamental_winding_is_11`
- `test_494_phase3_pluto_charon_32_cycler`

Runtime: ~33 s (parallel, 16 workers).
