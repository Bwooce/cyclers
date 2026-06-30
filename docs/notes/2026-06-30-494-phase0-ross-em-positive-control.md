# #494 Phase-0 Positive Control — Ross EM (k₁,k₂)-Cycler μ-Family Lane

**Date:** 2026-06-30
**Task:** #494, Phase 0 (positive-control gate before any μ-extension)
**Verdict:** GO — all 5/5 Ross EM families recovered; period, topology, stability, and
independent closure all pass.

---

## Purpose

Before trusting Phase 2 (μ-extension from EM to Pluto-Charon), the recover-from-(μ,C,T)
machinery must faithfully reproduce all five catalogued Ross EM stable (k₁,k₂)-cycler
members. A 0/5 failure here would signal a broken lane (and invalidate Phase 2's all-negative
or all-positive as meaningless). This note records the gate result.

---

## Phase-0 gate tests

The complete Phase-0 gate is two existing test files:

| File | What it covers for all 5 families |
|---|---|
| `tests/search/test_cr3bp_ross_families.py` | corrector convergence, Jacobi = C^stable, period ≈ T^stable (tolerance column below), Barden stability \|ν\|<1 and ν≈0, independent Radau closure |
| `tests/search/test_binary_star_search.py` | `winding_topology` returns published (k₁,k₂), prograde, reaches_secondary — **all 5** (PR #494 Phase-0 added (2,1),(3,2),(3,3); (3,1) and (1,1) were pre-existing) |

Both files pass green:

```
tests/search/test_cr3bp_ross_families.py  19/19 passed
tests/search/test_binary_star_search.py   11/11 passed
```

---

## Per-member results table

Sources: Ross & Roberts-Tsoukkas 2025, AAS 25-621.
C^stable and T^stable are from Table 3 (p. 11); (k₁,k₂) from Definition 1 (p. 5).
μ = 1.2150584270572 × 10⁻² (p. 3).

| Family | C^stable (sourced) | T^stable TU (sourced) | T^recovered TU | \|ΔT\| TU | ν (Barden) | Stable? | Topology | Prograde | Xcheck |
|--------|-------------------|----------------------|----------------|-----------|-----------|---------|----------|----------|--------|
| (1,1)  | 3.151175879508174 | 10.2920692100797595  | 10.2920692262527513 | 1.6e-08 | −0.00334 | YES | (1,1) MATCH | YES | PASS |
| (2,1)  | 3.129389531088256 | 19.4404316679515396  | 19.4401604299533481 | 2.7e-04 | +0.05007 | YES | (2,1) MATCH | YES | PASS |
| (3,1)  | 3.161784147013429 | 14.7884924166814002  | 14.7882679438847049 | 2.2e-04 | +0.01545 | YES | (3,1) MATCH | YES | PASS |
| (3,2)  | 3.182762663084288 | 17.9005801035000616  | 17.9005801012602142 | 2.2e-09 | −0.01174 | YES | (3,2) MATCH | YES | PASS |
| (3,3)  | 3.177224018696528 | 18.1454605758918888  | 18.1454605759149708 | 2.3e-11 | +0.06001 | YES | (3,3) MATCH | YES | PASS |

### Period tolerance discussion

The corrector enforces C^stable algebraically (Jacobi matches to ≤ 1e-12 for all members).
The period residual |ΔT| is controlled by the width of the stable subfamily:

- **(1,1)**: very wide window → |ΔT| = 1.6e-08 TU (below 5e-8; tighter than the 1e-5 test tolerance)
- **(2,1)**: razor-thin window (Δ_pm = 4.23 km) → |ΔT| = 2.7e-04 TU (within 5e-4 test tolerance)
- **(3,1)**: narrow inner band (Δ_pm = 253.70 km, narrow in x₀ space) → |ΔT| = 2.2e-04 TU (within 5e-4)
- **(3,2)**: correct half-period crossing (6th x-axis crossing) → |ΔT| = 2.2e-09 TU (well below 1e-5)
- **(3,3)**: widest window (Δ_pm = 2041.34 km) → |ΔT| = 2.3e-11 TU (machine-level)

The razor-thin (2,1) and narrow (3,1) families' larger |ΔT| values reflect the physical width
of the stable subfamily: the nu=0 midpoint sits at a C slightly different from the published
15-digit C^stable, so the corrected period at the printed C lands a hair off the nu=0 exact
midpoint. This is expected and documented in the catalogue `data_gaps` entries.

### What the Phase-0 gate confirms

For every member recovered from the sourced (μ, C^stable, T^stable):

1. **Period**: T^recovered ≈ T^stable to the tolerances above — the corrector is on the right
   branch.
2. **Topology**: `winding_topology` returns the published (k₁,k₂) label — the U1-/U2+ winding
   count matches.
3. **Stability**: Barden |ν|<1 — every member is linearly stable, matching the paper's
   "stable" verdict.
4. **Closure**: independent Radau cross-check (`crosscheck_periodic`) passes — the orbit is a
   genuine periodic solution, not an artefact of the single-shooter.

---

## Verdict: GO

5/5 Ross EM (k₁,k₂)-cyclers recovered. The recover-from-(μ,C,T) machinery is verified.
Phase 2 (μ-extension: Pluto-Charon instantiation) may proceed. Any all-negative or
all-positive result from Phase 2 is now meaningful (this gate eliminates the broken-lane
alternative explanation).

---

## Test additions in this pass

`tests/search/test_binary_star_search.py`:
- New parametrized test `test_494_phase0_winding_topology_all_5_em_families` covering all 5
  families with sourced (C^stable, T^stable) as inputs. Previously the binary_star test only
  checked (3,1) and (1,1) topology; this adds (2,1), (3,2), (3,3).

---

## Notes

- The Jacobi constant is enforced algebraically by `correct_symmetric_fixed_jacobi` (dJ < 1e-12
  for all members; not tabulated above since it is exact by construction, not a measured tolerance).
- Barden ν values are DERIVED quantities (our computation), not sourced goldens. The sourced
  VERDICT ("stable") is the golden; our ν values self-consistently reproduce it.
- The winding-topology `prograde` flag requires both w₁ > 0 and w₂ > 0 (both primaries wound
  counterclockwise in the rotating frame); all 5 Ross families satisfy this.
- `reaches_secondary` (x_max > x_L1) is true for all 5, confirming each orbit crosses the L1
  neck into the lunar realm as expected for a cycler.
