# CR3BP stable-family continuation search -- results (engine-validation campaign)

**Run timestamp:** 2026-06-12T11:01:27Z   elapsed 106.9s
**Script:** `scripts/cr3bp_family_search.py`
**Method:** `cr3bp-jacobi-continuation-v1`  (git `f99ba30`)
**Inertial gate:** `cr3bp-inertial-rebound-ias15-v1` (REBOUND/IAS15, independent frame + code path)
**Spec:** `docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md`
**Status:** SILVER-only. NO catalogue writeback.

## Headline (honest framing -- read this first)

- **Distinct NEW families found: 0.** This is the only count that
  can support a discovery claim, and it is ZERO. Every seed in this campaign is a
  *published* family representative (Ross Table 3 / Arenstorf / a sourced #182
  Lyapunov). Continuing such a seed in the Jacobi constant walks ALONG the already-
  published family -- each stepped member is the SAME family at a neighbouring
  energy, NOT a new orbit. A clean ZERO with the engine validated is a SUCCESS, not
  a failure: it is exactly what the orbit-closure discipline expects when every
  seed is sourced and the JPL 3-body catalogue (#116) is not yet acquired.
- **This corrects an earlier inflated count.** A prior pass labelled every non-
  point-reproduction member 'NOVEL-SILVER' (peaking at 83 'novel' / 30 'stable
  novel'), because the dedup compared distance to a single published *member* while
  Ross prints one representative per family -- so every other point on the same
  family curve counted as a discovery. The classification now keys on FAMILY
  identity / branch continuity, not point-distance:

- Seeds run: **9**
- Members kept (gauntlet-passing, across both directions): **57**
- REPRODUCTION (published representative re-found): **8**
- KNOWN-FAMILY-CONTINUATION (same sourced family, stepped in C; NOT a discovery): **49**
  -- of which linearly STABLE (|nu|<1): 3 (additional members of
  Ross's *named* stable windows -- a same-model cross-check, not new stable families).
- NOVEL-SILVER (member on a branch NOT tracing to a sourced seed = distinct new family): **0**
- Review-queue rows written: 49  (`data/cr3bp_continuation_review_queue.jsonl`)
  -- all KNOWN-FAMILY-CONTINUATION, flagged same-family (NOT discoveries).
- EMPTY-region records written: 9  (`data/empty_regions.jsonl`)

## Model & scope (caveats)

- **Planar CR3BP (PCR3BP) only.** The richest stable families (halos, DROs) are
  3D and need a 3D corrector (deferred).
- **Novelty scope = 'not in our sourced set AND not a continuation point of a
  sourced family'.** A genuinely-novel member would have to sit on a branch that
  does not trace back to a published seed. The JPL Three-Body Periodic Orbit
  Catalog is NOT acquired (#116); even a future genuinely-novel member is 'new vs
  our sourced set', NOT a literature-novelty claim.
- **Reproduction / same-family continuation is the engine-validation SUCCESS.**
  Re-deriving Ross's families (and stepping smoothly along them with nu varying
  continuously and Jacobi conserved) validates the continuation engine against
  sourced truth. That -- not a discovery count -- is what this campaign delivers.

## Per-seed census

| seed | system | mu | members | repro | known-fam-cont | stable-cont | new families | branches (dir: stop @ steps) |
|---|---|---|---|---|---|---|---|---|
| ross-(1,1) | Earth/Moon | 1.215058e-02 | 1 | 1 | 0 | 0 | 0 | +1: fold_radicand@0; -1: topology_jump@1 |
| ross-(2,1) | Earth/Moon | 1.215058e-02 | 1 | 1 | 0 | 0 | 0 | +1: topology_jump@1; -1: topology_jump@1 |
| ross-(3,1) | Earth/Moon | 1.215058e-02 | 17 | 1 | 16 | 0 | 0 | +1: topology_jump@1; -1: fold_reversal@17 |
| ross-(3,2) | Earth/Moon | 1.215058e-02 | 3 | 1 | 2 | 0 | 0 | +1: topology_jump@1; -1: fold_reversal@3 |
| ross-(3,3) | Earth/Moon | 1.215058e-02 | 32 | 1 | 31 | 3 | 0 | +1: fold_reversal@8; -1: jacobi_bound@24 |
| arenstorf-1963 | Earth/Moon | 1.227747e-02 | 0 | 0 | 0 | 0 | 0 | seed off-family: corrected T=5.032046 vs published T=17.065217 (>2%); crossing index / sign do not match |
| saturn-mimas-lyapunov | Saturn/Mimas | 6.598788e-08 | 1 | 1 | 0 | 0 | 0 | +1: fold_radicand@0; -1: jacobi_bound@0 |
| saturn-enceladus-lyapunov | Saturn/Enceladus | 1.901073e-07 | 1 | 1 | 0 | 0 | 0 | +1: fold_radicand@0; -1: topology_jump@1 |
| saturn-tethys-lyapunov | Saturn/Tethys | 1.086441e-06 | 1 | 1 | 0 | 0 | 0 | +1: fold_radicand@0; -1: no_converge@1 |

## Inertial cross-check (false-consensus independence gate)

Each kept member is re-propagated in the inertial REBOUND/IAS15 harness (a
different code path AND frame than the rotating-frame DOP853 corrector) and
graded against LINEAR-INSTABILITY THEORY -- the per-period amplification is the
monodromy eigenvalue magnitude `|lambda| = |nu + sqrt(nu^2-1)|` from the member's
Barden nu. R1 (one-period recurrence <= 0.1A), R2 (Jacobi drift <= 1e-9 over
the bounded span), R3 (observed departure no EARLIER than the linear prediction
`t_dep = T + ln(3A/delta1)/ln|lambda|`). The band is NEVER loosened. A member
that fails this gate is NOT routed (classified CONTINUATION-REJECTED).

**Honest limit of this gate (do not overread a PASS).** The corrector lands
every member to a ~1e-10 perpendicular-crossing residual, so the inertial
one-period recurrence delta1 is ~1e-9..1e-8 -- and even a strongly unstable
member (|nu|~10^2, |lambda|~10^2..10^3) does not amplify that tiny seed past 3A
within the 5-period span, so it stays numerically bounded and R3 records
departure LATER than predicted (a PASS). The inertial gate therefore confirms
each member is a GENUINE periodic orbit through an independent integrator+frame
(it re-closes and conserves Jacobi) -- it is the false-consensus *consistency*
gate. It does NOT independently certify STABILITY: stability is the Barden nu
(reported per member), which IS discriminating here (nu spans -0.0 stable to
+360 wildly unstable across these branches). Read 'inertial PASS' as 'a real
orbit, cross-checked', and the STABLE verdict as 'from nu', never conflated.

## STABLE same-family continuation members (NOT discoveries, NOT catalogued)

**These are NOT new families.** Every stable member below is a continuation
member of one of Ross's *named* (k1,k2) families, at a Jacobi constant a few
steps off the published nu=0 midpoint but INSIDE that family's finite stable
window (Ross Table 3 gives the window widths; e.g. (3,3) spans ~2041 km in
perilune). They are additional members of the published stable WINDOW -- a same-
model cross-check, a SUCCESS -- not new stable families. No literature-novelty
is claimed (JPL 3-body catalog not acquired, #116). NEVER catalogue.

| family | C | T (nd) | nu | x0 | ydot0 | inertial |
|---|---|---|---|---|---|---|
| ross-(3,3) | 3.17782402 | 18.134746 | -0.55589 | -0.32077990 | -1.82898125 | PASS |
| ross-(3,3) | 3.17792402 | 18.139980 | +0.01877 | -0.32077035 | -1.82900645 | PASS |
| ross-(3,3) | 3.17692402 | 18.429612 | -0.37121 | -0.32520754 | -1.80508496 | PASS |

## Distinct new families

**None.** By construction: every seed is a sourced family representative, so
continuation can only re-derive or step along a published family. No branch in
this campaign traces to anything other than a sourced seed. Recorded as method-
versioned EMPTY-for-novelty (re-sweepable when the JPL catalog (#116) or a 3D
corrector ships, or when an UN-sourced seed is introduced).

## Discipline statement

Sourced goldens only (Ross Table 3 / Arenstorf); same-model validation; novelty
keyed on FAMILY identity / branch continuity (a continuation point of a sourced
family is the SAME family, not a discovery), not point-distance; SILVER-only with
NO catalogue writeback; every Jacobi range/step/stop-reason logged (no silent caps).

## Audit appendix (#219 adversarial re-exam of the #218 claims)

Manually appended after the run above (the sections above are script-generated).

- **The "30 stable novel" claim is RETRACTED.** 25 of them were ross-(3,1)
  continuation steps at dC=1e-4 along the ONE published (3,1) family (nu drifting
  smoothly) -- points on a known family curve. Under the corrected family-keyed
  accounting, (3,1) contributes 16 KNOWN-FAMILY-CONTINUATION members (0 stable)
  and 0 discoveries. Honest census: 0 distinct new families.
- **The "arenstorf-seeded stable T=5.032 / nu=-0.55" member: false provenance,
  now gated.** Direct re-check (this audit): with the Arenstorf IC (x0=0.994,
  C~2.8564) and half_crossings=1, the fixed-Jacobi corrector converges to
  x0=+1.31087, ydot0=-0.65647, T=5.032046 -- a GENUINE in-model periodic orbit
  (crossing residual 8e-15, Radau dJ 3e-14, amplitude 0.70 nd; not degenerate, and
  NOT a subharmonic of Arenstorf: T_pub/T = 3.391, no integer multiple matches).
  But it is NOT the Arenstorf figure-eight (canonically UNSTABLE, T=17.0652; our
  golden pins |nu|>1), so continuing it under the 'arenstorf' label would
  fabricate provenance. The SEED_PERIOD_TOL_FRAC=2% off-family guard now abandons
  the seed (0 members emitted) -- exercised by
  `test_offfamily_seed_is_abandoned_not_fabricated`. The T=5.032 orbit itself is
  left unrouted: it surfaced from a failed seed correction, not a provenanced
  search, and carries no source anchor.
- **The (1,1)-branch "stable member": eliminated by the period-continuity gate.**
  The (1,1) walk now stops at topology_jump@1 (the corrector landed on a
  different perpendicular-crossing topology, period discontinuous in C); only the
  seed reproduction survives. The earlier 'stable (1,1) novel' was that
  off-family jump member.
- **The 3 (3,3)-branch stable members: REAL distinct orbits, KNOWN family.** They
  re-converge identically in this re-run (same C/T/nu to print precision), pass
  the full gauntlet + independent Radau + inertial IAS15 gates, and sit INSIDE
  Ross's published (3,3) finite stable window -- additional members of a sourced
  stable window (engine validation), not new stable families.
- **Review-queue routing now agrees with the census**: 49 rows written = 49
  KNOWN-FAMILY-CONTINUATION members (+ 0 novel); the old run wrote rows labelled
  NOVEL-SILVER (and the original #218 pass wrote 0 rows against 83 claimed
  novels). The routing test now exercises a real queue write and the count
  identity (it was vacuous when novel==0).
