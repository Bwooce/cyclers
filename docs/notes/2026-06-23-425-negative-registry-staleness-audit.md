# #425 — negative-results registry staleness audit (method-subsumption + bug-fix)

**2026-06-23, HEAD 4daa977.** Periodic audit required by
`feedback_bugfix_invalidates_past_searches`: on every correctness fix, scan the
recorded negatives for any whose **binding constraint** the fix loosens, and re-run
only those. A negative is invalidated **only** if a fix *loosens a binding
constraint* — for a flyby floor that means a *lower* floor (more bend capability);
for an epoch correction it means a phase shift large enough to change a survival
verdict. A *tighter* constraint (raised floor) can never resurrect a negative.

## Scope: correctness changes since the last audit

The empty-region sweeps were last reverified on **2026-06-14** against the **#205
Lambert dT/dz bug** (`f6a0460`) — all stamped "no flip" (gaps stayed 7.8–21.2 km/s
vs the 6.0 floor). Since then two correctness changes have landed that could in
principle move a recorded negative:

- **F1 — #198 epoch fix** (`439d279`): UTC↔TDB 63 s correction on **real-ephemeris**
  evaluation. In scope only for negatives whose verdict comes from a real-eph
  (DE440) integration.
- **F2 — #426/#428 flyby-floor fix** (`fbe7e56` + `4daa977`): lowered **design**
  floors — Earth/Mars 300→200 km, Mercury 1000→200 km. (Callisto 100→200 km is a
  *raise* = tightening, cannot invalidate anything.) In scope only for negatives
  whose binding constraint is an Earth / Mars / Mercury flyby-bend budget.

The #411 cross-system corrector-direction bug-fix and the #205 Lambert fix are
already reflected in the registry (the #405 entry's resweep text; the empty-region
reverification stamps) and are out of scope here.

## Registry contents audited

- `data/negative_results.yaml` — 10 characterized structural negatives.
- `data/empty_regions.jsonl` — 21 bounded discovery-sweep "empty" reports.

**31 negatives total**, grouped into 4 classes by what could possibly flip them.

## Class-by-class verdict

### Class A — Floquet real-eph V4 negatives (4)
`branch_C32_b0` (#389), `branch_C32_C_3.1774` (#392), `branch_C11a_b0` (#393),
`em_node_floquet_branches` (#392).

These *are* real-ephemeris, so **F1 is nominally in scope**. But the failure margin
is structural and enormous: EM-only control drift ~7×10⁵ km vs full N-body drift
~1.4×10⁹ km — a four-order-of-magnitude exponential blow-up driven by solar tides
(~30 % of Earth gravity at far amplitude) and Floquet-multiplier instability
(~37 at low amplitude) amplifying real lunar eccentricity. A 63 s epoch shift
perturbs the initial phase by <10⁻⁴ of an orbit; it cannot convert a 1.4×10⁹ km
exponential escape into a bounded orbit. These are libration/cycler-node families
with **no flyby**, so F2 is irrelevant. **NOT invalidated.**

### Class B — Earth–Mars single-arc ballistic V1 negatives (4)
`russell-ocampo-3.1.2+1`, `-4.3.1-5`, `-4.5.2-2` (#365), `mcconaghy-2006-em-k2`
(S1L1, #365).

These involve Earth + Mars flybys, so **F2 is nominally in scope**. But the failure
is a **topology mismatch at the closure step**: a single radial-crossing free-return
ellipse fundamentally cannot represent a multi-arc geometry, breaking V∞ continuity
by *tens of km/s*. The flyby floor only sets the max-bend budget at the *linkage*
step, which is downstream of — and irrelevant to — a closure that never forms; and a
bend-*angle* budget cannot repair a velocity-*magnitude* discontinuity of tens of
km/s. The closers carry no real-eph epoch, so F1 is irrelevant. Their resweep gates
(per-arc free-return descriptor acquisition / multi-arc DSM closure infrastructure,
#388/#404) are untouched by either fix. **NOT invalidated.**

> Note: `mcconaghy-2006-em-k2` is exactly the S1L1 row the M7 maintenance work
> (#423) reproduced strictly-ballistically at the new 200 km floor. That is a
> *maintenance-ΔV* result on the constructed cycler, not a closure of the
> single-arc free-return — it does **not** promote this V1 closure negative, whose
> resweep remains the multi-arc closer.

### Class C — patched/coupled CR3BP & BCR4BP structural negatives (2)
`cross_system_se_em_L2_patched_cr3bp` (#405/#411),
`bcr4bp_phase_b_em_libration_seed` (#412).

Autonomous CR3BP / BCR4BP models — **no real-ephemeris epoch**, so F1 is irrelevant.
Libration manifolds with **no flyby**, so F2 is irrelevant. The failures are
phase-closure *dimensionality* (1-DOF obstruction, confirmed at single and
low-multiple revolution) and amplitude-*reach* structure (EM-libration seed stays
EM-bounded). Neither fix touches these. **NOT invalidated.**

### Class D — empty-region discovery sweeps (21)
Jovian/Saturnian VILM (4), CR3BP Jacobi-continuation (9), repeated-moon
multi-rev (5), DA/HOTM new-family band (1), binary-star μ-continuation (2).

- **F1 (epoch):** every entry's `ephem_model` is `circular` or `cr3bp` — **none use
  real DE440**, so the UTC/TDB correction touches no empty-region verdict.
- **F2 (floor):** flyby bodies are outer-planet moons (Io/Europa/Ganymede/Callisto;
  Dione/Rhea/Titan/Enceladus/Tethys/Mimas; Uranian/Neptunian moons) and
  Mars-Phobos/Deimos — **none use Earth, Mars, or Mercury as a flyby body**, so the
  lowered E/M/Mercury floors are not a binding constraint of any entry. The single
  body touched (Callisto 100→200) is a *raise* (less bend), which can only tighten;
  and the Jovian/Saturnian gaps are 7–21 km/s in **V∞ magnitude**, not bend-limited.

The 2026-06-14 #205-Lambert reverification already stamped the four VILM entries
"no flip". **NOT invalidated.**

## Conclusion

**0 of 31 negatives are invalidated** by the #198 epoch fix or the #426/#428
flyby-floor fix. In every case the binding constraint is provably untouched:
real-eph negatives fail by margins ~10⁹ km that 63 s cannot move; ballistic
negatives fail at a closure/topology step upstream of any flyby-bend budget; the
floor change is a *raise* (Callisto) or touches a body (E/M/Mercury) that is not a
flyby target of any recorded negative; and no discovery sweep used real ephemeris.

Because each binding constraint is provably unaffected, **no re-run is warranted** —
re-running deterministic analytic/CR3BP sweeps whose gating quantity cannot change
would only reproduce the same verdict. This is the correct application of the
bug-fix-invalidation rule, not an exception to it.

Each entry is stamped with a `reverification` (empty_regions) / `audit_425`
(negative_results) record pointing here, so a future correctness-fix audit can see
this registry was already checked against F1+F2 at HEAD 4daa977.
