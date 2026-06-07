# Multi-rev back-arc Lambert in `dsm_leg` + re-probes — task #157

**Date:** 2026-06-07
**Code:** `src/cyclerfinder/search/dsm_leg.py` (additive), `tests/search/test_dsm_leg.py`
(additive). `core/lambert.py`, `search/mbh.py` read-only (imported, never edited).
**Context read (in order):** `docs/notes/2026-06-07-dsm-full-sequence-probe.md`
(the #153 diagnosis: the back-arc Lambert is single-rev, the loop arcs are 1.5-3
revs; single-rev `|v1|`=28.2 vs best multi-rev 15.5 km/s on the 1315.5-d G-arc),
`src/cyclerfinder/search/dsm_leg.py`, `src/cyclerfinder/core/lambert.py` (the
`max_revs` API + `n_revs`/`branch` conventions from M4),
`tests/search/test_dsm_leg.py`, `docs/notes/2026-06-07-dsm-leg-genome.md`,
`docs/notes/multi-arc-classification.md` §7/§12 (the three-rows / two-framings
finding — use each row's OWN anchors, never mix framings).

## What this builds

`dsm_leg(..., max_revs=0, rev_branch=None)` — the back-arc Lambert can now take a
multi-revolution branch:

* `max_revs=0` (default) → **bit-identical** to the prior single-rev path
  (`sols[0]`; regression `test_dsm_leg_max_revs_zero_is_bit_identical_to_default`).
* `max_revs>0` → enumerate the single-rev branch plus every feasible multi-rev
  `low`/`high` branch up to `max_revs` (via `core/lambert.lambert(max_revs=...)`),
  and select the one **minimising the DSM impulse** `||v21 - v12||` (the leg
  objective). The chosen `(n_revs, branch)` are recorded in new audit fields
  `DsmLegResult.n_revs_chosen` / `.branch_chosen`.
* `rev_branch=(n, "low"|"high")` → force an explicit branch (raises
  `LambertError` if that branch is infeasible for the geometry/ToF).

`max_revs` threads through `evaluate_dsm_chain` / `dsm_chain_correct` /
`make_dsm_chain_step`; per-leg `n_revs_per_leg` / `branch_per_leg` are carried on
`DsmChainResult` and surfaced in the MBH-step `info` for audit.

The selection criterion is dV-minimisation, not branch-count maximisation — so the
optimiser is free to keep single-rev wherever that is genuinely cheaper. That is
the honest behaviour: multi-rev is *available*, never *imposed*.

## Regression

`tests/search/test_dsm_leg.py` — the 9 prior tests (6 fast + 3 slow from
#153/#150) stay **green, unchanged**. Full suite (slow + not slow): **14 passed**
(9 fast incl. 3 new mechanics gates; 5 slow incl. 2 new re-probes). ruff check +
ruff format --check on both files: clean. mypy on `dsm_leg.py`: clean (the test
file's 8 mypy notes — import-untyped + two pre-existing `type: ignore` — are
identical on HEAD; not introduced here).

## Multi-rev mechanics gate (CONSTRUCTED) — verbatim

`test_dsm_leg_max_revs_recovers_multirev_improvement_on_g_arc`: on the 1315.5-d
G-arc M→E return (the loop arc that floors the single-rev full sequence), a
near-pure-Lambert leg (η=0.01, modest departure V∞), `dsm_leg`:

```
single-rev (max_revs=0):  |v1| = 27.97 km/s   n_revs=0 branch=single   dV_DSM = 18.35
multi-rev  (max_revs=3):  |v1| = 15.14 km/s   n_revs=3 branch=high     dV_DSM = 10.17
```

This reproduces #153's measured Lambert-family numbers (28.18 → 15.51 km/s) **at
the `dsm_leg` primitive level**: |v1| drops ~13 km/s, the DSM impulse drops ~8
km/s, and a genuine multi-rev branch is selected. Reference = the leg's own
single-rev result + the multi-rev Lambert family (golden-rule mechanics, not a
golden test). Label: **mechanics.**

## Re-probe (a): 6.44Gg3 full E-M-E-M-E sequence, `max_revs=3` — verbatim

(MBH cauchy, rng_seed=6, ≤120 hops, stall 60. Same seeding as the #153 single-rev
baseline; only `max_revs=0→3` changed. EXPECTED = the row's sourced anchors
V∞ E=6.44 / M=3.74; emerged is EVIDENCE.)

```
feasible=False   total_dV = 26.8793 km/s
hops attempted/accepted = 61/0   (stopped on stall)
emerged V_inf_in  {1(E→M):11.16, 2(M→E):15.66, 3(E→M):3.26, 4(M→E):5.85}   sourced E=6.44, M=3.74
per-leg dV_DSM = (7.37, 7.00, 9.12, 3.39) km/s
n_revs per leg = (0, 0, 0, 1)   branch per leg = (single, single, single, low)
eta per leg    = (0.452, 0.476, 0.476, 0.442)
tof days/leg   = (263.9, 502.6, 262.1, 1315.8)
wall = 61.6 s
```

**Outcome: NEGATIVE, but improved and now mechanically representable.** total ΔV
fell from the #153 single-rev **29.94 → 26.88 km/s**; the long G-arc leg 4
selected a genuine multi-rev branch `(1, "low")` and its emerged V∞ dropped (7.55
→ 5.85). The g-arc leg 2 stayed single-rev (the optimiser found single-rev cheaper
there at the landed η/ToF). The floor is lower than the single-rev baseline — the
loop arc is now representable — but the basin still does not close: emerged
encounter V∞ (3-16 km/s) remain off the sourced 6.44/3.74, 0 of 61 hops accepted.
The multi-rev branch removed the *specific* degeneracy #153 identified on the long
leg, but the full multi-arc basin is not reached from this seed under this genome.

## Re-probe (b): S1L1 `russell-ch4-4.991gG2` two-arc, `max_revs=3` — verbatim

**Anchors used are this row's OWN** (catalogue `russell-ch4-4.991gG2`, Russell 2004
Table 4.9 free-return framing): g 1.4612 yr + G 2.8096 yr, transit out/in 150 d,
V∞ E=**4.99** / M=**5.10** km/s, aphelion 1.64 AU, period 4.27 yr. **NOT** the
`s1l1-2syn-em-cpom` 5.65/3.05 framing (a different idealisation of the same
physical cycler; multi-arc-classification §7/§12 — framings are never mixed). Same
E-M-E-M-E unrolling, MBH cauchy seed 6.

```
feasible=False   total_dV = 38.8809 km/s
hops attempted/accepted = 61/0
emerged V_inf_in  {1:11.80, 2:9.94, 3:18.20, 4:37.14}   sourced E=4.99, M=5.10
per-leg dV_DSM = (13.10, 13.58, 6.59, 5.61) km/s
n_revs per leg = (0, 0, 0, 0)   branch per leg = (single, single, single, single)
eta per leg    = (0.389, 0.484, 0.257, 0.340)
tof days/leg   = (152.0, 388.0, 149.0, 866.4)
wall = 55.8 s
```

**Outcome: NEGATIVE, off-basin.** The dV-minimising selection kept all four legs
single-rev at the landed genome (the arcs here are shorter — g 1.46 yr, G 2.81 yr
— so the multi-rev minimum-time threshold is closer, and at the converged η/ToF
single-rev was the cheaper branch). total ΔV floors at 38.88 km/s, emerged V∞
(10-37 km/s) far from the sourced 4.99/5.10, 0/61 accepted. **No close-and-match;
not the first multi-arc closure.** The multi-rev primitive is exercised and
available, but this seed/genome does not reach the 4.991gG2 basin.

## Verdict

The #153 blocker — *the back-arc Lambert cannot represent the multi-rev loop arcs*
— is **fixed and exercised**: `dsm_leg` now selects multi-rev branches, the
mechanics gate reproduces #153's 28→15 km/s improvement at the primitive level,
and re-probe (a) shows the long G-arc leg taking a real `(1, "low")` branch with a
lower floor than the single-rev baseline. But **neither re-probe closes**: both
6.44Gg3 (26.88 km/s, improved) and S1L1 4.991gG2 (38.88 km/s) floor off-basin with
0 accepted hops. Making the topology *representable* was necessary but not
sufficient — the remaining gap is basin/seed selection for the full multi-arc
sequence (the η/ToF/branch combination that ties the two arcs together at the
sourced encounter V∞), not a missing mechanic.

**No catalogue writeback** (per task constraint; both rows stay multi-arc with
their sourced anchors).

## Follow-up

The frontier is now seed/basin selection, not mechanics. Candidate next steps
(not done here):

1. **Branch as a genome coordinate.** The dV-minimising per-leg selection is
   greedy/local; the globally-best closure may need a *fixed* `rev_branch` per leg
   (e.g. force the g arc to `(1,·)` and the G arc to `(2,·)`) explored as a
   discrete MBH coordinate, rather than re-minimised inside each evaluation.
2. **Symmetric / period-locked seeding.** Russell's generic arcs are symmetric
   (t_in = t_out); seeding the loop arcs at their true half-arc geometry rather
   than the (transit, arc−transit) split may land nearer the basin.
3. **Coupled V∞ continuity at the Mars/Earth flybys** as an explicit residual
   (the chain currently charges only the DSM impulses + terminal arrival; the
   flyby V∞-match that *defines* a ballistic cycler is the scorer's job and is not
   in this objective).
