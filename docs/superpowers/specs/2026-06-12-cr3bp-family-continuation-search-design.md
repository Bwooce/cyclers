# CR3BP Stable-Family Continuation Search — Design Spec

**Date:** 2026-06-12
**Status:** Draft (design). Execution gated on #212b (Ross adoption) landing the
fixed-Jacobi symmetric corrector + Barden stability index.
**Track:** The project's FIRST genuine novel-discovery campaign on a frontier with
real headroom. Follows CR3BP Tier-2 (#182, the propagator + corrector + 14 SILVER
pipeline) and the Ross & Roberts-Tsoukkas 2025 adoption (#212).

## Goal

Discover new CR3BP periodic-orbit family members — especially STABLE ones — by
natural-parameter continuation along the Jacobi constant, seeded from sourced
published members, across the Earth-Moon and midsize-moon (Saturnian + others)
systems. Validate each discovery in-model (Jacobi conservation + independent
Radau + inertial n-body cross-check) and route every new orbit to the review
queue as SILVER. Record empty regions in the method-versioned registry.

## Why this, why now

Every prior discovery lane is **ceiling-bound or empty**:
- The heliocentric patched-conic novelty sweep (Forge Phase 6) came back
  essentially empty with the current genome; the catalogue is at its
  data-limited validation ceiling and past it is new-input-gated.
- The 14 SILVER Saturnian Lyapunovs (#182) proved the discovery→gauntlet→
  cross-check pipeline works end-to-end, but they were tiny-amplitude and
  **unstable** (error amplification ~2000×/period) — pipeline-validation, not
  literature-novel.

The CR3BP periodic-orbit space is the frontier that is NOT exhausted, and Ross &
Roberts-Tsoukkas 2025 (AAS 25-621, mined 2026-06-11) just supplied the two things
the lane lacked: a **systematic construction method** (the fixed-Jacobi symmetric
single-shooting corrector, their Eq. 12, ~5 iterations to ε=1e-8) and **sourced
stable seeds** (5 (k1,k2) Earth-Moon cycler families with 15-digit C/T, μ printed
to 1e-8, and stable subfamilies — the lane's first chance at a |ν|<1 verdict).
Continuation along Jacobi from those seeds walks each family and harvests members
the single-shot corrector alone cannot reach. This is a real search, in a model
where discoveries get same-model validation (unlike the off-family patched-conic
lanes).

## Dependencies (hard)

- **#212b** must land: the Ross fixed-Jacobi symmetric corrector
  (`correct_periodic` variant or a sibling) and the Barden half-period monodromy
  stability index ν = ½(λ + 1/λ). This campaign continues families with that
  corrector; without it there is no engine. The EM μ double-count fix (#212a,
  done) is also prerequisite — continuation needs the correct μ.
- The existing Tier-2 core (`core/cr3bp.py`: EOM, Jacobi, 42-state STM,
  `cr3bp_system`; `search/cr3bp_periodic.py`: corrector + `crosscheck_periodic`),
  the review queue, and the inertial REBOUND/IAS15 cross-check
  (`scripts/cr3bp_silver_inertial_crosscheck.py`, #182).

## Binding constraints (orbit-closure-discipline + false-consensus defence)

1. **Sourced goldens, never our own computation.** The corrector/continuation is
   validated by reproducing PUBLISHED members: the Arenstorf orbit (Tier-2
   golden) and Ross's 5 families' published (C, T) to their printed precision
   before any new member is trusted. EXPECTED values cite the publication.
2. **Same-model validation only.** A discovered orbit's Jacobi/period/stability
   are compared against in-CR3BP invariants, NEVER against a different-fidelity
   value (the cross-fidelity bug class). State μ and the system for every result.
3. **All binding constraints in the convergence test.** A member is converged
   only with: closure residual < tol AND period in bounds (no collapse) AND not
   an equilibrium (max|v| gate) AND dedup-distinct AND Jacobi conserved ≤1e-10
   over the period AND `crosscheck_periodic` (independent Radau) agrees. The
   "it closed!" danger signal applies — a libration point trivially closes.
4. **Independent inertial cross-check (per the false-consensus doctrine).** Each
   kept member is re-propagated in the INERTIAL n-body harness (different code
   path AND frame), graded against linear-instability theory (the #182 method:
   unstable members can't survive literal N-period boundedness in double
   precision — grade departure-time-consistency, never loosen the band).
   This is the campaign's independence gate; the corrector's own residual is a
   consistency check only.
5. **Stability is a first-class output, honestly reported.** ν computed from the
   half-period monodromy; |ν|<1 = stable, |ν|>1 = unstable. Do NOT advertise a
   member as stable without the ν verdict. The campaign's headline interest is
   the stable subfamilies, but unstable members are recorded too.
6. **No catalogue writeback.** Every new member → `review_queue.jsonl` (SILVER),
   never auto-promoted. Dedup against Ross's published members, Arenstorf, the
   Genova-Aldrin 3-petal (bicircular, out-of-model), and the JPL Three-Body
   Periodic Orbit Catalog if acquired (#116). A member matching a published one
   is a REPRODUCTION (sourced cross-check), not a discovery — label it so.
7. **Empty regions are method-versioned negatives.** A Jacobi range that yields
   no new family member under this method+version is recorded in the registry,
   re-sweepable only when a subsuming capability ships. No silent truncation —
   log every bound (Jacobi range scanned, step, families seeded).

## Architecture — continuation on the Tier-2 foundation

### Phase 0 — corrector/stability acceptance (blocks everything)
Confirm the #212b corrector reproduces ALL 5 Ross families' published (C, T) to
printed precision AND returns the published stable/unstable ν verdict for each.
This is the same-model golden gate; if it fails, stop (the engine is wrong).

### Phase 1 — natural-parameter continuation driver (`search/cr3bp_continuation.py`)
```python
def continue_family(system, seed: PeriodicOrbit, *, dC, n_steps,
                    corrector, max_jacobi, min_jacobi) -> FamilyBranch
```
Pseudo-arclength or natural-parameter (Jacobi) continuation: predict next member
(tangent/secant from the last STM), correct with the fixed-Jacobi symmetric
corrector, run the convergence gauntlet (constraint 3), compute ν (constraint 5),
stop at a fold/turning point or the Jacobi bound. Returns the ordered branch with
per-member (state0, T, C, ν, residuals). Golden: walking from Ross's seed
reproduces the published stable-window edges (their perilune-width values).

### Phase 2 — multi-system seeded campaign (`scripts/cr3bp_family_search.py`)
Seed continuation from: (a) Ross's 5 EM families; (b) Arenstorf; (c) the Tier-2
Saturnian Lyapunov seeds; (d) JPL-catalogue ICs if #116 lands. Walk each family
both directions in Jacobi, run the gauntlet + inertial cross-check per member,
dedup, classify each as REPRODUCTION (matches a published member) / NOVEL-SILVER
(new, validated) / EMPTY (registry). Parallelise across families/systems on the
16-core box (the #182 process-pool pattern).

### Phase 3 — disposition + registry
NOVEL-SILVER members → review queue with full provenance (system, μ, state0, T, C,
ν, all gate residuals, the inertial-crosscheck verdict). Empty Jacobi ranges →
method-versioned registry. A results note with the family census (reproduced vs
novel vs empty), the stable-member count (the headline), and honest caveats.

## What counts as a discovery

A NOVEL-SILVER member is one that (a) passes the full gauntlet + independent
inertial cross-check, (b) is NOT within dedup tolerance of any published member
(Ross, Arenstorf, JPL catalogue, existing catalogue rows), and (c) carries a
complete in-model tuple (μ, state0, T, C, ν). Stable novel members (|ν|<1) are
the campaign's primary target — the project has zero stable discoveries to date.
A member matching a published orbit is a sourced REPRODUCTION (valuable as a
same-model cross-check, but not new) and is labelled as such, never as discovery.

## Out of scope (YAGNI / deferred)

- **Bicircular / full-ephemeris refinement** of CR3BP members (Tier-3; the
  Genova-Aldrin 3-petal lives here and is permanently out of pure-CR3BP scope).
- **3D / halo families** beyond what Ross's planar method + the corrector's
  current planar capability reach — the corrector is planar today; 3D halo
  continuation (e.g. from the Cuevas L1 halo seed, #190) is a follow-on once a
  3D corrector exists.
- **The multi-arc / repeated-moon genome** (the Liang Jovian re-sweep) — a
  separate sub-project; this campaign is CR3BP periodic orbits, not patched-conic
  resonant moon tours.
- **Auto-promotion** of any discovery — always review-gated.

## Testing (TDD)

- `cr3bp_continuation.py`: walking from each Ross seed reproduces the published
  (C, T) of an adjacent member and the published stable-window edge; a seed at a
  known fold turns the branch (no run-past); convergence gauntlet rejects an
  equilibrium / period-collapse injected mid-branch (fault-injection).
- campaign: a seeded family yields ≥1 member OR a clean EMPTY record; a member
  matching a Ross-published one is labelled REPRODUCTION; a novel member routes
  to the review queue, never the catalogue; the inertial cross-check runs on
  every kept member.
- regression: full suite green; the Arenstorf + Ross same-model goldens hold.

## References

- `docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md` — the 5 families,
  construction method, 15-digit (C, T), stability windows, the μ-bug find.
- `docs/superpowers/specs/2026-06-10-cr3bp-tier2-design.md` — the Tier-2
  foundation this builds on.
- `docs/notes/2026-06-10-cr3bp-moontour-results.md` — the 14 SILVER pipeline
  (degeneracy gates, `crosscheck_periodic`); `...-cr3bp-silver-inertial-crosscheck.md`
  — the inertial-cross-check method (linear-instability grading).
- Memory: `orbit-closure-discipline` (incl. the false-consensus defence),
  `validation-ceiling`, `negative-results-registry`, `gmat-install`.
- Acquisition leads (#116): Ross 2026 journal version (3D/non-symmetric families,
  resolves the C-bound table inconsistencies); JPL Three-Body Periodic Orbit
  Catalog (free — sourced halo/DRO seeds).

## Risks / honest caveats

- **Planar-only until a 3D corrector exists** — the richest stable families
  (halos, DROs) are 3D; this campaign harvests planar families first. State this
  in every result; it is not "the CR3BP searched," it is "planar CR3BP families
  from these seeds."
- **Reproduction-heavy first pass is a SUCCESS, not a failure** — if the first
  campaign mostly reproduces Ross/Arenstorf/JPL members, that validates the
  continuation engine against sourced truth (and is the same-model golden the
  discipline wants); genuine novel members are the upside, not the bar.
- **Dedup completeness gates the novelty claim** — a "discovery" that is actually
  a known orbit we failed to dedup is a false positive. Acquire the JPL catalogue
  (#116) before claiming heliocentric-of-the-3-body novelty; until then, scope
  novelty claims to "not in our sourced set" and say so.
