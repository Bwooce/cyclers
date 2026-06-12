# Lambert completeness-fix blast-radius re-run — RESULTS (task #206)

**Date:** 2026-06-12
**Trigger:** the #205 Lambert completeness fixes (commit `f6a0460`,
`core/lambert.py`): (A) the `_dt_dz` √C algebra error that made Newton oscillate
and false-raise `LambertConvergenceError` on long-way (dnu>π) single-rev
transfers, and (B) the log-compressed Illinois residual that rescues multi-rev
high branches the raw-residual stall silently dropped. Both are FALSE-NEGATIVE
generators: any recorded EMPTY / OFF-FAMILY / NO-CLOSE verdict whose pipeline
traversed a long-way or multi-rev Lambert leg may have been wrong.
**Question:** does any recorded negative actually FLIP under the fixed solver?
**Writeback:** NONE. No catalogue change. A flip would go to REVIEW; none found.
**New runlogs (method-versioned, fixed-Lambert HEAD `b27f251`):**
`data/runs/self-seeding-triage-20260612T0150Z.jsonl`,
`data/runs/self-seeding-reachable-20260612T0150Z.jsonl` (APPENDED, the
2026-06-09 originals are untouched).

## HEADLINE

**No recorded negative flips.** The Lambert fix is real and observable — it
recovers conics (long-way solutions that previously false-raised; multi-rev
branches that were dropped) and measurably changes the surveys (more legs close,
one VEM floor dropped 18.54→16.34 km/s) — but in every flagged survey the
recovered conics live in the SAME high-|v∞| basin. The Jones gap stays a
|v∞|-MAGNITUDE gap (floors ~16-18 km/s, ZERO bend-feasible, vs sourced 2.4-7.0);
the self-seeding / dsm verdicts reproduce bit-identically. The fix changes the
conic SET, not the reachable FAMILY.

## 0. The fix is genuine and quantified (OLD `f6a0460^` vs NEW)

A dense OLD-vs-NEW solver sweep (11 232 cases: transfer angle 20-350°, radius
ratio 0.6-2.5, tof 80-1300 d, max_revs 0-2) confirms the conic set changed:

- **318 / 11232 (2.8%) feasibility flips** — long-way transfers (dnu≈200-300°)
  where OLD raised `LambertConvergenceError` and NEW returns a valid single-rev
  solution (defect A; matches the ~4% false-raise estimate).
- **226 branch-set changes** — at long tof OLD dropped a multi-rev branch
  (typically `(1,'low')`) that NEW recovers (defect B).
- **0 velocity deltas** — where BOTH solve, |v1| is bit-identical. The fix only
  ADDS previously-missing solutions; it never perturbs an existing one.

So the re-runs below are a fair test: the conic set genuinely changed in exactly
the long-way / multi-rev regimes these surveys traverse.

## 1. VEM ballistic rediscovery survey

- **Test:** `tests/test_vem_rediscovery.py::test_jones_vem_ballistic_rediscovers_sourced_multiset` (slow xfail).
- **Driver:** `scripts/hunt_vem_ballistic.py` (parallel epoch×branch scan; topologies
  include 1-rev low/high per leg — the defect-B path).
- **Recorded floor (xfail reason, hunt 2026-06-06, 256 epochs):** EMEVVE best
  max-V∞ 17.86; MEEVEM 18.49; ZERO bend-feasible. Gap to sourced 2.42-7.0 km/s.
- **Re-run config (cheapest decisive):** `hunt_vem_ballistic.py 16 astropy magnitude`,
  OLD `f6a0460^` lambert vs NEW, identical grid (176 pts/member).

| member | OLD floor | NEW floor | OLD closed | NEW closed | bend-feasible |
|---|---|---|---|---|---|
| EMEVVE | 18.54 | **16.34** | 51/176 | 70/176 | 0 / 0 |
| MEEVEM | 18.58 | 18.58 | 64/176 | 95/176 | 0 / 0 |

- **New verdict:** the fix recovers conics (closed count up 51→70, 64→95) and
  lowers EMEVVE's floor 18.54→16.34, but the floor stays ~16-18 km/s with ZERO
  bend-feasible, vs sourced 2.42-5.16. **FLIP? NO.** Magnitude gap unmoved.

## 2. Jones n-body shooter survey

- **Test:** `tests/nbody/test_shooter_jones_gate.py::test_jones_vem_nbody_rediscovers_sourced_multiset` (slow xfail).
- **Driver:** `scripts/hunt_vem_nbody_shooter.py`; the Lambert-dependent stage is
  the `near_miss_survey` (single-rev/single-branch default → defect A applies to
  every long-way leg). It gates whether the shoot gets a low-V∞ seed at all.
- **Recorded floor (xfail reason):** near-miss surfaces only the high-V∞
  Lambert-chain basin, max-V∞ ~29-33 km/s, 0 bend-feasible vs Jones 2.5-5.2.
- **Re-run config:** `hunt_vem_nbody_shooter.py 64 0` (near-miss survey only, 0
  shoots — the survey IS the decisive Lambert-dependent stage), OLD vs NEW.

| member | OLD seeds | NEW seeds | OLD best seed | NEW best seed | bend-feasible |
|---|---|---|---|---|---|
| EMEVVE | 30 | 44 | 19.09 | **18.28** | 0 / 0 |
| MEEVEM | 34 | 46 | 18.57 | 18.57 | 0 / 0 |

- **New verdict:** more near-miss seeds recovered (30→44, 34→46), EMEVVE floor
  19.09→18.28, but still ~18 km/s, ZERO bend-feasible, vs sourced 2.5-7.0. No
  low-V∞ seed exists to shoot into the Jones basin. **FLIP? NO.**

## 3. #177 self-seeding triage (OFF-FAMILY / NO-CLOSE + multi-rev G-arc)

- **Drivers:** `scripts/triage_self_seeding.py` (cheap transit-match gate,
  no n-body) and `scripts/validate_self_seeding_reachable.py` (full tail:
  multi-rev best branch → longitude rendezvous → REBOUND/IAS15 confirm). Both use
  the multi-rev `g_arc_branches` / `_refine_lambert` path (defects A+B reachable).
- **Recorded (`docs/notes/2026-06-08-self-seeding-triage-results.md`,
  runlogs `...-20260609T0205Z`/`...-20260609T0210Z`):** 6 REACHABLE / 2 NO-CLOSE
  / 204 NO-DESCRIPTOR; all 6 REACHABLE → OFF-FAMILY-AT-ANCHOR-VINF (0 V3).
- **Re-run:** full triage (212 rows) + reachable validator (6 rows), NEW lambert,
  fresh timestamp `20260612T0150Z` (appended, not overwritten).

| stage | recorded | re-run (fixed Lambert) |
|---|---|---|
| triage distribution | 6 REACHABLE / 2 NO-CLOSE / 204 NO-DESCRIPTOR | **identical** |
| 6 REACHABLE best branches/Δd | 9.353Gg2 short −8 … 6.44Gg3 **long** +30 | **identical** |
| reachable validation | 6× OFF-FAMILY-AT-ANCHOR-VINF | **6× OFF-FAMILY-AT-ANCHOR-VINF, identical emerged v∞** |

  6.44Gg3 long branch still emerges v∞_M 7.83 (anc 3.74); 9.353Gg2 short 21.01
  (anc 10.52); etc. — every value matches the 2026-06-09 run to the printed digit.

- **New verdict:** bit-identical. The 2 NO-CLOSE rows (`5.30gGf3`, `5.75ggF3`)
  stay NO-CLOSE — confirming their non-reach is pure coplanar geometry (aphelion
  sub-Mars), NOT Lambert convergence. **FLIP? NO.**

  > **Cross-reference (not a #205 effect):** the descriptor-row OFF-FAMILY verdict
  > was ALREADY corrected on 2026-06-10
  > (`docs/notes/2026-06-10-dsm-tof-artifact-correction.md`): the Stage-B closer
  > used the coplanar branch ToF (`shape.tof_g_days`) instead of the row's
  > signature transit, inflating emerged Mars v∞ ~1.6-2.1×; a free (epoch, ToF)
  > Lambert at the signature ToF reproduces BOTH anchors to 0.1-0.3 km/s on ~6
  > rows. That ToF-selection fix is the LIVE flip for these rows and is a SEPARATE
  > concern — the #205 algebra fix does not touch it, which is exactly why this
  > re-run reproduces the (still-buggy-ToF) OFF-FAMILY numbers identically. Do not
  > conflate the two: #205 = no flip; the 2026-06-10 ToF correction = the real,
  > already-documented flip awaiting its corrected Spec-2 implementation.

## 4. MBH / dsm multi-rev long-way runs

- **Run:** the #157 multi-rev probe (`docs/notes/2026-06-07-dsm-multirev-probe.md`)
  — `dsm_leg(..., max_revs=3)` selecting multi-rev branches on the 1315.6-d
  (6.44Gg3) and 764.2-d (S1L1 4.991gG2) G-arc M→E loop arcs, the longest /
  most-multi-rev legs in the dsm lane and the prime defect-B candidates.
- **Recorded:** 6.44Gg3 best multi-rev |v1| 15.51 km/s (max_revs=3), floor 26.88
  km/s off-basin, 0/61 hops; S1L1 4.991gG2 |v1| 21.28, 38.88 km/s off-basin.
  (The dsm-multiarc "off-family" closure run, `2026-06-10-dsm-multiarc-closure-results.md`,
  is already SUPERSEDED/WITHDRAWN by the same 2026-06-10 ToF correction above.)
- **Re-run:** OLD vs NEW lambert on the exact #157 loop-arc geometries,
  max_revs 0-3.

| arc | max_revs | OLD branch set / best |v1| | NEW |
|---|---|---|---|
| 6.44Gg3 (1315.6 d) | 3 | 7 branches, 15.508 | **identical** |
| S1L1 4.991gG2 (764.2 d) | 3 | 3 branches, 21.279 | **identical** |

  (max_revs 0/1/2 also bit-identical for both arcs.)

- **New verdict:** the multi-rev branch SET and best |v1| are bit-identical
  OLD↔NEW on both loop arcs — these legs sit comfortably inside the Illinois
  bracket and never hit the stall regime defect B touches. The #157 floors would
  reproduce exactly. **FLIP? NO.**

## 5. Honest yield / interpretation

| survey | recorded floor / verdict | re-run floor / verdict | FLIP |
|---|---|---|---|
| VEM ballistic | ~17.9 / 18.5, 0 bend-feasible | 16.34 / 18.58, 0 bend-feasible | NO |
| Jones n-body near-miss | ~29-33, 0 bend-feasible | 18.28 / 18.57, 0 bend-feasible | NO |
| #177 self-seeding (212+6) | 6 REACH→0 V3, 2 NO-CLOSE | identical | NO |
| #157 dsm multi-rev loop | 15.51 / 21.28 |v1|, off-basin | identical | NO |

The #205 fix is correct and its effect is observable (recovered conics, more
closures, a couple of lowered floors), but it changes the conic SET WITHIN the
high-|v∞| basin — it does not bridge the |v∞|-magnitude gap to the Jones family,
and it does not perturb the dsm/self-seeding branch selections that already
converged correctly. This is the expected and valuable negative: a mirror-free
DE440 control independently floored at 18.16, and the gap is magnitude (15-21 vs
2.4-7.0 km/s), not a missing conic. **No verdict goes to review on #205 grounds.**

## 6. Registry / discipline

- New runlogs APPENDED under fixed-Lambert HEAD `b27f251` (`code_version` field);
  the 2026-06-09 originals (`b0683f4`) are untouched — method-versioned per the
  negative-results-registry rule ("empty" is conditional on the method; the
  fixed-Lambert sweep is a new method version of the same region).
- The xfail markers were NOT edited (triage's separate call). Both VEM xfails
  remain correctly xfail: no member converges within `VEM_VINF_TOL_KMS`.
- No catalogue writeback. Nothing promoted on the search's say-so.
