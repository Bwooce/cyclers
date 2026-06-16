# 325 Topology-gate audit — defensive sweep post-#322

**Date.** 2026-06-16.
**Scope.** #325 Phase 1 (Parts A — defensive audit; B — independent-topology
sweep harness; C — minimum-necessary gate fixes; D — verdict + Phase 2).
**Outcome.** **Sibling-bug count: 0**. The audit confirms #322 was the only
#322-class topology-gate bug in the discovery-corpus JSONLs (24 files, 552
findings, 0 discrepancies after the honesty filter).

## #322 lesson

#322 (commit `c2a77c7`) found that `find_tulip_at_system::petal_count`
misclassified PLANAR Np-petal orbits as 3D tulips at extreme small-mu regimes
(e.g. Mars-Phobos, μ ~ 1.65e-8): the symmetric corrector drove z0 → 0 and
the in-plane projection still had Np petals so the count-only gate fired.
The fix added `cyclerfinder.genome.tulip.is_three_dimensional`, a
complementary check on max|z(t)| against the Koblick-sourced floor.

**Generalised bug pattern.** Any topology gate that
1. Counts a feature (petals, encounters, revolutions) without checking the
   FULL state dimensionality;
2. Was designed for a specific regime (typical μ, typical amplitudes) without
   bounds-testing extreme regimes;
3. Lacks a complementary "is this the topology I think it is" check.

## Part A — per-module audit verdicts

Seven modules audited. Verdict per module:

| # | Module | Verdict | Notes |
| - | ------ | ------- | ----- |
| 1 | `search/cr3bp_general_periodic.py` | **CLEAN** | Convergence gate is the 2x2 return-map residual + `closure_residual` (independent Radau re-propagation). Orbit is z=0 by construction; no #322-class collapse possible. |
| 2 | `search/nrho_continuation.py` | **VULNERABLE (low severity)** | The convergence gate is `residual < tol`. At z0 → 0 the residual is trivially satisfied. No topology gate at the corrector level. **Mitigated**: the downstream caller `find_tulip_at_system` runs `is_three_dimensional` post-#322. |
| 3 | `search/bifurcation_detector.py` | **CLEAN** | Returns brackets, not topology claims. Distance from `eig` to primitive k-th roots of unity. Real -1 IS period-doubling (correct k=2). |
| 4 | `search/tisserand_mga_window.py::linkable_3d` | **CLEAN** | The Galileo Jupiter V∞=11 km/s pin is the real ballistic-Tisserand threshold; Galileo's multi-planet leveraging is out of scope for this predicate. |
| 5 | `five_tier_prioritizer.legs_from_repeated_moon_candidate` | **N/A** | Planar by construction in the planet frame. No 3D-topology claim to falsify. |
| 6 | `genome/multi_shooting.py::multi_shoot_periodic` | **VULNERABLE (low severity, same as #2)** | Same z=0 collapse pattern as #2. **Mitigated**: only production caller is `_multi_shoot_switch` which flows through `switch_family` (Part C `verify_three_dimensional`). |
| 7 | `genome/bcr4bp_genome.py::correct_bcr4bp_periodic` | **CLEAN** | Independent Radau cross-check already in place; `sun_phase_drift` diagnostic for caller visibility. |

**Summary**: 4 clean, 2 vulnerable-but-mitigated, 1 N/A.

## Part B — independent-topology audit harness

New module: `src/cyclerfinder/search/topology_audit.py`. Three sourced
independent checkers:

1. **`check_tulip_topology`** — re-runs the `#322` `is_three_dimensional`
   gate. Sourced to Koblick 2023 Table 4.
2. **`check_periodic_orbit_closure`** — independent Radau full-period
   re-propagation; closure residual must be below 1e-6 nondim. Sourced to
   `correct_bcr4bp_periodic`'s `independent_tol` (1e-6).
3. **`check_floquet_neimark_sacker`** — classifies a Floquet eigenvalue
   pair as Neimark-Sacker iff BOTH `|λ| - 1` within 1e-2 AND `|Im(λ)|`
   strictly above 1e-3. Sourced to Gomez et al. 2001 §3.

The top-level `audit_topology(jsonl_path)` walks a discovery JSONL and
dispatches each row to applicable checker(s). Rows tagged
`degenerate_planar=True` by the genome are skipped — the family tracer is
honestly reporting a planar terminal member, not making a 3D claim.

Output: `TopologyAuditFinding` dataclass per (row, checker) pair.

New script: `scripts/run_325_topology_audit.py`. Sweeps the harness across
24 session JSONLs (#296, #299, #312, #313, #285, #298, #309, #311, #302).

Output: `data/topology_audit_325.jsonl` (header + 552 findings + summary;
24/24 files audited).

Tests: `tests/search/test_topology_audit.py` — 12 tests, all pass.

## Part C — fixes for vulnerable gates

Single minimum-necessary edit: **`genome/family_switch.switch_family`**.

The original gate (line 297-304) is `petal_count == k`. The fix adds an
opt-in `verify_three_dimensional` parameter (default `False` to preserve
behaviour) that runs `is_three_dimensional` on the switched member's IC
and rejects sub-floor z0 / max|z|.

Why opt-in: the canonical Earth-Moon discovery path
(`find_tulip_via_continuation` → `find_tulip_at_system`) already runs the
3D gate downstream, AND the NRHO test suite exercises legitimate small-z0
ICs at Earth-Moon mu. Recommended ON for multi-system / extreme-mu callers
where z0 → 0 collapse is the dominant failure mode.

Modules NOT modified:
* `nrho_continuation.py` — only callers gate downstream; adding the check
  would break the small-z0 NRHO test contract.
* `multi_shooting.py` — same: caller flows through `switch_family`.

Tests added: 3 in `tests/genome/test_family_switch.py`:
1. **Positive control**: Earth-Moon Np=2 still admitted with the flag on.
2. **Backward-compat**: default `verify_three_dimensional=False` produces
   bit-identical results to the pre-#325 call.
3. **Synthetic regression**: `is_three_dimensional` rejects sub-floor z0.

## Part D — re-verification yield + Phase 2

### Yield (re-verification of past JSONL rows)

* **Files audited**: 24/24.
* **Total findings**: 552 (across multiple checkers per row).
* **Discrepancies before honesty filter**: 1.
* **Discrepancies after honesty filter**: 0.

The single pre-filter discrepancy was the terminal member of the #296
Earth-Moon 3D family (`step_index=138`, `state_nd[2] ~ -2.5e-10`). The
family tracer had ALREADY tagged this row as `degenerate_planar=True` —
the genome correctly identified its own family-forward termination at the
planar limit. The audit's honesty filter (introduced after observing this
single case) skips rows that the genome already honestly labeled.

**Net yield**: 0 past JSONL rows downgraded.

### Honest assessment

This audit confirms hygiene rather than catching a sibling bug. #322 was
isolated. The reasons are clear in retrospect:

* The other genome topology gates are either residual-based (no topology
  claim), or already paired with an independent cross-check.
* The family tracer modules honestly label their own degeneracies.
* The single risk-bearing path (`switch_family` → potential planar
  collapse) is already filtered by `find_tulip_at_system` downstream.
  Part C hardens it at the `switch_family` layer too.

### Phase 2 recommendation

Adopt a **convention** rather than a separate audit pass:

> Every new topology gate must come with a complementary check that asks
> the orthogonal question. Example: a petal counter must be paired with a
> 3D-amplitude check; a Floquet root-of-unity classifier must check that
> the eigenvalue is genuinely complex; a Lambert closure gate must check
> V∞ preservation AND physical realisability (the #324 physical-sanity
> gate covers the latter).

The `topology_audit` harness is the standing tool for re-running this
audit on demand.

## Files

* Harness: `src/cyclerfinder/search/topology_audit.py`
* Tests: `tests/search/test_topology_audit.py` (12 tests)
* Script: `scripts/run_325_topology_audit.py`
* Output: `data/topology_audit_325.jsonl` (552 findings, 0 discrepancies)
* Fix: `src/cyclerfinder/genome/family_switch.py` (verify_three_dimensional)
* Fix tests: `tests/genome/test_family_switch.py` (+3 regression tests)
* Doc: `docs/notes/2026-06-16-325-topology-gate-audit.md` (this file)
