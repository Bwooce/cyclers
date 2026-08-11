# #820 — russell12 lambert genome re-posed per #794's designated-arc semantics: truth is NOW a closure point for 8/12 rows

**Date:** 2026-08-11
**Task:** `#820` (registered by `#813`, dispatched 2026-08-11) — re-pose
`scripts/campaign_russell12.py::build_genome` so the DESIGNATED (uppercase-letter) arc in
Russell's own leg-descriptor notation is the Mars-transit leg regardless of its position in
`free_return_arcs[]` (per `#794`'s primary-source verification: Russell 2004 §4.8 pp. 125-127;
McConaghy/Russell/Longuski 2005, JSR 42(4), DOI 10.2514/1.8123), seed the E-E loops from
`#794`'s written-back `loop-ee-*` segments, re-check the `3.64gGg3`/`8.049gGf2`-class rows
`#794` flagged, and re-run the campaign.

## Headline

**The #135 "sourced geometry is not a residual-zero point of this genome" diagnosis is
OVERTURNED — it was an artifact of the mis-posing, not a property of the rows.** Under the
corrected posing, the seed-at-truth probe converges to a genuine residual-zero closure within
0.2-3.2 days of the sourced geometry on **8 of 12 rows**, with the EMERGED per-body V∞ matching
the independently sourced anchors to ≤ 0.1 km/s (documented tolerance 0.5); one row
(`5.30ggF3`) passes even the strict at-seed STAYED-AT-TRUTH criterion, and the blind epoch-grid
campaign produced its **first-ever CLOSE-AND-MATCH verdicts** (`9.353Gg2`, `3.78Gg3` — every
anchor check OK). The 4 remaining rows fail for a now precisely-identified structural reason
(full-rev resonant arcs in the LOOP set are not Lambert-representable — registered `#825`),
not a diffuse "family-selection wall". This is the V1-grade multi-arc closure evidence beyond
the single-ellipse free-return genome that `#820`'s registration hoped for; adjudication /
writeback deliberately NOT done here (registered `#826`).

## 1. What was mis-posed (three coupled defects, all fixed)

The pre-#820 `build_genome` assumed `arcs[0]` is the Mars free-return arc and `arcs[1:]` are
the E-E loops, with the two sourced taxi transits (`out-em`/`ret-me` segment ToFs) as the E→M
and M→E leg ToFs. Against `#794`'s primary-source semantics this was wrong three ways:

1. **Designated-arc position.** The designated (uppercase) arc is `arcs[0]` for only 3 of the
   12 rows. Census by `raw_descriptor` case (now the code's identification rule):

   | designated position | rows |
   |---|---|
   | `arcs[0]` (old assumption correct) | 9.353Gg2, 3.78Gg3, 9.94Gg3 |
   | `arcs[1]` (G) | mcconaghy-2006-em-k2, 4.991gG2, **8.049gGf2**, **3.64gGg3**, 5.30gGf3, 6.44Gg3 |
   | `arcs[2]` (F) | 3.66gfF3, 5.30ggF3, 5.75ggF3 |

   So 9/12 rows assigned a designated arc to the loop set and dropped a real loop — including
   the two 1:1-resonance rows `#794` flagged for re-check (`3.64gGg3`, `8.049gGf2`: both
   designated at `arcs[1]`, both were mis-posed; the flag was warranted).
2. **Loop seeds.** Loops are now the non-designated arcs in cyclic itinerary order after the
   designated leg, seeded (ToF, n_revs, branch) directly from `#794`'s written-back
   `loop-ee-*` segments (primary-source-derived; `resonant` branch label mapped to the
   corrector's multi-rev `low`). A per-position cross-check (descriptor arc ToF vs segment
   `tof_days`, tol 1 d) raises if the mapping ever drifts; all 12 rows pass.
3. **The beyond-Mars remainder (the deepest defect).** The M→E leg was seeded with the sourced
   t_in — but t_in is the INBOUND-crossing taxi transit (ridden on another cycle's leg), not a
   per-cycle leg ToF. The designated leg's own beyond-Mars time (printed full ToF minus t_out;
   e.g. 1026.2 − 153 = 873.2 d for mcconaghy vs the 153 d the old genome used) had no leg to
   live in, so the seeds could NEVER tile the period and the slack-leg elimination
   reconstructed a garbage loop ToF (1259.9 d where the sourced loop is 533.7 d). Leg 1 is now
   the designated-arc REMAINDER after the Mars encounter at t_out; the full seed vector now
   tiles the sourced period to ~1 d (residual tiling error is purely `period.years`'s
   2-decimal print rounding, absorbed by the slack leg).

Leg 1's `(n_revs, branch)` is not printed anywhere: for f/F designated arcs it is exactly
N−1 revs (the E→M transit is a fraction of one conic period); for g/G arcs the printed
transfer angle bounds it to floor(θ/360) or one less, and the low/high Lambert branch is
undetermined. The small candidate set is enumerated (`_leg1_topology_candidates`) and
discriminated by **residual-at-truth** (`select_leg1_topology` — the same
sourced-geometry-as-discriminator principle `#794` used for Lambert-root selection; losing
candidates are reported alongside, nothing is loosened).

`_classify`'s transit check was re-posed to match: achieved t_out (leg 0) vs the sourced
`transit_times_days[0]`, plus the achieved designated-leg total (legs 0+1) vs the printed
descriptor ToF; t_in has no per-cycle leg to compare against under the corrected posing.

## 2. Re-run (12 rows, circular model, matching #813's scale)

`uv run python scripts/campaign_russell12.py --model circular --epochs 256 --workers 8
--probe-at-truth --phase-epochs 256 --runlog-timestamp 20260811T-820-reposed` (run in two
`--rows` batches, same runlog). Runlog:
`data/runs/russell12-circular-20260811T-820-reposed.jsonl` (12 records, appended
method-versioned; the June-07 and `#813` records untouched). Full per-row probe JSON captured
at 256 phase epochs (values below).

| row | grid outcome | truth res (km/s) | solved res | ToF drift (d) | emerged V∞ E/M | sourced E/M |
|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2 | CLOSE-OFF-ANCHOR | 0.307 | **0.000** | 3.14 | 5.008/5.107 | 4.70/5.00 (McConaghy-flavored; matches the Russell 4.99/5.10 twin — the known #794 model caveat) |
| 4.991gG2 | CLOSE-OFF-ANCHOR | 0.267 | **0.000** | **0.22** | 5.008/5.107 | 4.99/5.10 ✓ |
| 8.049gGf2 | CLOSE-OFF-ANCHOR | 35.22 | 5.079 | 291.2 | — | 8.05/10.02 ✗ (#825 wall) |
| 9.353Gg2 | **CLOSE-AND-MATCH** | 0.182 | **0.000** | 0.71 | 9.386/10.532 | 9.35/10.52 ✓ |
| 3.64gGg3 | CLOSE-OFF-ANCHOR | 39.21 | 0.000* | 204.0 | 13.6/12.6 | 3.64/4.59 ✗ (#825 wall; *converged far off-truth) |
| 3.78Gg3 | **CLOSE-AND-MATCH** | 0.367 | **0.000** | 0.97 | 3.748/4.622 | 3.78/4.63 ✓ |
| 5.30gGf3 | CLOSE-OFF-ANCHOR | 1.972 | 1.913 | 2.10 | — | 5.30/9.17 ✗ (#825 wall, inclined f(3:2) loop) |
| 9.94Gg3 | CLOSE-OFF-ANCHOR | 1.315 | **0.000** | 0.99 | 9.838/10.718 | 9.94/10.76 ✓ |
| 3.66gfF3 | CLOSE-OFF-ANCHOR | 39.87 | 3.957 | 123.9 | — | 3.66/4.66 ✗ (#825 wall) |
| 5.30ggF3 | CLOSE-OFF-ANCHOR | **0.074 (STAYED-AT-TRUTH)** | **0.000** | 0.66 | 5.270/5.419 | 5.30/5.44 ✓ |
| 5.75ggF3 | CLOSE-OFF-ANCHOR | 0.349 | **0.000** | 1.59 | 5.697/9.337 | 5.75/9.36 ✓ |
| 6.44Gg3 | CLOSE-OFF-ANCHOR | 0.213 | **0.000** | 2.18 | 6.497/3.795 | 6.44/3.74 ✓ |

Reading: "truth res" = min-over-phase residual evaluated AT the (rounded) sourced seed;
"solved res"/"drift" = the corrector run seeded exactly there. For the 8 ✓ rows the corrector
converges to an exact closure (residual 0.000, period-constrained, ballistic) within
0.2-3.2 d of the sourced geometry whose emerged V∞ reproduces the printed anchors to
0.01-0.10 km/s — V∞ is never imposed in this genome, so this is genuine evidence, not
imposition. Versus `#813`'s post-M:N-fix baselines on its 4 rows: 5.30ggF3 26.11 → 0.074,
5.75ggF3 17.00 → 0.349, 5.30gGf3 6.89 → 1.97 (residual left = the out-of-ecliptic f(3:2)
loop, see below), 3.66gfF3 26.69 → 39.87 (now dominated by the degenerate f(1:1) loop leg —
a different, sharper failure, see below).

### Honest caveats (nothing loosened)

- The strict STAYED-AT-TRUTH verdict (at-seed residual ≤ 0.1 km/s) fires only on `5.30ggF3`.
  For the other ✓ rows the at-seed residual sits at 0.18-0.37 km/s — dominated by seed print
  rounding: `period.years` alone (4.27 vs the exact 4.2708 = 2×2.1354-yr synodic) injects
  ~1.3 d of slack-leg error. The closure claim rests on the corrector-converged results under
  the documented tolerances (drift ≤ TOL_TRANSIT_DAYS = 5 d; V∞ within TOL_VINF = 0.5 km/s),
  not on any adjusted threshold.
- The blind 256-epoch grid found the anchor basin on only 2 rows (its ~6-9 d t0 spacing steps
  over a basin the phase-scan + corrector resolves); the CLOSE-OFF-ANCHOR grid labels on ✓
  rows reflect grid sampling, not absence of the anchor solution.
- mcconaghy's probe closure emerges at Russell's 4.99/5.10 anchor, not the row's own
  McConaghy-2006 4.7/5.0 — exactly the model caveat `#794` recorded on that row.

## 3. The remaining wall is now precisely identified (#825)

The 4 ✗ rows are EXACTLY the rows whose **loop set** contains a full-rev resonant `f` arc:

- `8.049gGf2`, `3.64gGg3`: in-ecliptic f(1:1) loop (365.25 d, 1 rev, e ≠ 0). A 1:1 loop leg
  departs and re-meets Earth at the SAME point after one Earth year — a 360°-transfer-angle
  Lambert problem, which is DEGENERATE (the (a,e) family through one point with fixed period;
  the printed loop e = 0.2678/0.1220 is not recoverable from endpoints+ToF). The 35-39 km/s
  truth residuals are this degeneracy, not geometry.
- `3.66gfF3`: same f(1:1) degeneracy, plus its loop is genuinely inclined (λ=87.388 ⇒
  i = 7.038°, from `#794`).
- `5.30gGf3`: f(3:2) loop (2 revs, 720° — same degeneracy class) that is also genuinely
  out-of-ecliptic (λ=118.851 ⇒ i = 6.589°); the residual floor ~1.9 km/s it converges to is
  the out-of-plane component a coplanar corrector cannot zero.

Note the F-DESIGNATED rows (`5.30ggF3`, `5.75ggF3`) close fine: the designated arc is split at
the Mars encounter into E→M / M→E legs, which are NOT degenerate. The wall is specifically
resonant arcs in the loop position. Fix path registered as `#825` (resonant-loop-aware leg:
pin the loop conic from the segment's (a,e[,i]) instead of solving a Lambert).

## 4. Code changes

`scripts/campaign_russell12.py`:
- `_designated_arc_index` (uppercase `raw_descriptor` = designated; raises unless exactly one),
  `_g_arc_transfer_angle_deg`, `_leg1_topology_candidates`, `_genome_with_leg1`,
  `select_leg1_topology` (residual-at-truth discrimination, candidates reported), `_t0_center`.
- `build_genome` re-posed as §1; raises (no silent pad) on arc/segment count or ToF mismatch;
  genome now carries `designated_index/designated_raw/designated_tof_days/leg1_candidates`.
- `run_row` and `probe_at_truth` select leg-1 topology before solving; verdict/probe dicts
  carry the selection + losing candidates. `_classify` transit check per §1.

`tests/search/test_russell12_likeforlike_probe.py`: the two lambert-genome pins re-pointed at
the NEW finding on `russell-ch4-4.991gG2` (converges from truth, emerged V∞ matches the
sourced 4.99/5.10 anchors; genome splits the designated `arcs[1]` and tiles the period). The
old pins guarded the mis-posed genome's #135 diagnosis and are superseded by exactly the bug
they were guarding against regressing into.

No catalogue edits in this task. `tests/scripts` preflight AST ratchet not applicable (script
is not `scripts/run_*.py`; suite run anyway, green).

## 5. Follow-ups registered

- **`#825`** — resonant-loop-aware genome leg for the 4 `#825`-wall rows (§3).
- **`#826`** — adjudicate the validation-evidence implications of the 8 truth-region closures
  (per-row: which V-level gate this satisfies, canonical `close_row_dsm` framing per the #388
  gate discipline, mcconaghy model caveat) and do any writeback as its own task — deliberately
  NOT done here per the dispatch.

## 6. Verification

- `uv run ruff check .` / `ruff format --check .` — clean; `uv run mypy src tests` — clean
  (843 files).
- `uv run pytest tests/data tests/search -q` — green (run in per-letter chunks under heavy
  CPU contention: a concurrent self-hosted CI job (12 workers) + `#822`'s concurrent
  screening run; zero FAILED/ERROR across all chunks). `tests/scripts -q` — green.
- `tests/search/test_russell12_likeforlike_probe.py -m 'slow or not slow'` — 4 passed
  (includes the two rewritten slow pins).
