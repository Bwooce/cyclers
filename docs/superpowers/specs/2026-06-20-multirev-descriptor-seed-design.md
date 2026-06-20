# Multi-rev wiring into the DSM descriptor-seed lane — Design Spec

**Date:** 2026-06-20
**Status:** Design approved (brainstorming), pending implementation plan.
**Amends:** `docs/superpowers/specs/2026-06-10-dsm-multiarc-closure-lane-design.md`
(the closure lane whose Component-4 batch run on 2026-06-20 produced a clean
structural negative — see `docs/notes/2026-06-20-dsm-closure-batch-results.md`).
**Tracks:** #307 / #388.

## Goal

Make the DSM descriptor-seed closure lane able to represent the **multi-revolution
resonant same-body legs** (E→E, M→M) that the Russell/McConaghy cyclers are built
from, by (a) seeding each leg at its **published** arc time-of-flight and (b)
enabling the *already-built* multi-rev Lambert branch selection the lane currently
bypasses. No new genome code — the leg-level multi-rev machinery exists and is
tested; only the descriptor→seed bridge is single-rev.

## Why — the 2026-06-20 batch negative, diagnosed

The closure-lane batch converged **0 / 9** descriptor-bearing rows (residuals
20–77 km/s; DSM dV 50–150 km/s), including validated V3 rows. Root cause: the
sequences are `E-E-M-M`; the same-body legs are ~2.4-revolution resonant returns,
but `close_row_dsm` calls `dsm_chain_correct` with the default `max_revs=0`
(single-rev), forcing each resonant leg onto the degenerate near-radial single-rev
Lambert solution. `dsm_leg.py` already documents this exact failure ("the single-rev
branch on a >1-period leg is degenerate — the #153 diagnosis") and already solves it
via `max_revs` + the dV-minimising branch selection. The bridge simply never passes
`max_revs`.

A second, larger defect surfaced during the literature correlation: the current seed
does not even use the **published** arc ToFs. For `russell-ch4-4.991gG2` /
`mcconaghy-2006-em-k2` it produces `[883, 143, 883]` d from a slack heuristic
(`big_G − transit`), when the sourced values are `[533.7, 150, 1026]` d. The lane was
seeding the wrong geometry; the multi-rev cap alone would not have fixed it.

## Literature correlation (the two parameters are sourced, not chosen)

Source paper for these cyclers: Russell, R. P., "Global Search and Optimization of
Free-Fall Cycler Trajectories," Ph.D. dissertation, UT Austin, 2004. Genome source:
Takao 2025 (arXiv:2501.06586), Appendix Eqs. A.1–A.3. Digests:
`docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md`,
`docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md`.

| Parameter | Sourced value | Citation |
|---|---|---|
| **Per-leg rev cap** | `ceil(arc_ToF_yr / body_period_yr) + 1`, globally ≤ Russell's **6-body-period** generic-return ceiling. → `max_revs = 2` for the E-E-M-M rows. | Russell §2.1 (ToF capped at 6 body periods; solutions binned by N); §2.4 (`p.h.s.i`, `i` = signed rev index, `i = −N_MAX..N_MAX`). |
| **Per-leg seed ToF** | The descriptor's **published** arc ToF: g/G `tof_years`×365.25 for the resonant legs; Russell `t_out = t_in = 150 d` (Table 4.9) for the E→M transit. | Catalogue `free_return_arcs[*].tof_years`; Russell 2004 Table 4.9. ToFs sum to 4.27 yr = Russell's 2-synodic repeat (Table 5.2). |
| **ToF box bounds** | Takao A.2/A.3 window **anchored at the published ToF** (bracket the sourced value, not a generic `0.5·P_body`). | Takao Appendix A.2 (`τ ∈ [30 d, P_s + P_H]`) / A.3; anchored at the Russell published arc ToF. |
| **"Auto / let it morph" rev selection** (chosen, not pinned) | Pass a `max_revs` cap; let `dsm_leg`'s existing dV-minimising branch selection pick per leg. | Russell §5.4.2: rejects pinning integer rev structure for gradient solvers ("the solution structure is free to morph"); names S1L1 as the favorable case. |

## Architecture

Three focused changes in `src/cyclerfinder/search/dsm_descriptor_seed.py`, plus a
batch re-run. The genome (`dsm_leg.py`), the corrector, the honesty gates, and the
golden-anchor separation are all unchanged.

### Change 1 — seed each leg at its published arc ToF

In `seed_dsm_chain_from_descriptor`, replace the slack heuristic
(`slack_days = max(30.0, big_G_tof_days − transit_days)`) with the **published**
per-leg ToFs:

- A same-body resonant leg (A == B, e.g. E→E, M→M) takes the published arc ToF for
  that arc: `g_tof_yr × 365.25` for the first resonant arc, `big_g_tof_yr × 365.25`
  for the second. The arc→leg assignment follows the existing g/G order in
  `free_return_arcs`.
- A cross-body transit leg (E→M / M→E) takes the row's **sourced** transit ToF from
  the existing structured field `invariants.transit_times_days` (required on every
  multi-arc row; `[150, 150]` d for `russell-ch4-4.991gG2` = Russell's tabulated
  `t_out = t_in = 150 d`), NOT the computed `arc.tof_g_days` (≈143 d). Transit legs
  map to `transit_times_days` entries in sequence order. No schema change or backfill
  — the field already exists; the seed simply reads it.

Add an explicit `per_leg_tof_days: tuple[float, ...]` to `DsmChainSeed` carrying the
sourced ToFs (audit/evidence), and use it to build `x0`.

### Change 2 — compute and thread the rev cap

In `seed_dsm_chain_from_descriptor`, add `max_revs: int` to `DsmChainSeed`:

```python
def _leg_rev_cap(arc_tof_days: float, body: str) -> int:
    period_days = 2*pi*sqrt((PLANETS[body].sma_au*AU_KM)**3 / MU_SUN_KM3_S2)/DAY_S
    return floor(arc_tof_days / period_days) + 1         # max complete revs in ToF + Russell §2.1 fast/slow headroom

RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP = 6               # Russell §2.1 ToF ceiling
max_revs = min(max(per_leg_rev_caps), RUSSELL_GENERIC_RETURN_BODY_PERIOD_CAP)
```

For a same-body leg the relevant `body` is that body; for a transit leg use the
inner body's period (shorter → larger N → the conservative cap). Take the max across
legs (a single global `max_revs` for the chain; the per-leg dV-min selection still
picks single-rev for the transit leg). In `close_row_dsm`, pass
`max_revs=seed.max_revs` to `dsm_chain_correct`. One added kwarg.

### Change 3 — ToF bounds bracket the published value

In the bounds built by `seed_dsm_chain_from_descriptor`, for each resonant same-body
leg set the lower/upper ToF bound to bracket the **published** arc ToF (e.g.
`[0.7×, 1.3×]` the sourced value, intersected with Takao's A.2 window where finite)
so the corrector cannot collapse the resonant leg back to the degenerate
near-zero-ToF single-rev solution. The transit leg keeps the existing
sequence-keyed bound. The infinite-`P_s` same-body cap (#217) the code already
handles stays.

### Re-run

Re-run `scripts/dsm_closure_batch.py` on real DE440. No script change beyond what it
already reads from the seed (it already records `n_revs_per_leg` via the result). Add
the emerged `n_revs_per_leg` to the runlog/summary line so the multi-rev branch
selection is visible.

## Data flow (unchanged shape; sourced seed)

row (descriptor + published arc ToFs + Russell V∞ anchor)
→ `seed_dsm_chain_from_descriptor` (published ToFs + rev cap + bracketed bounds)
→ `close_row_dsm` (`dsm_chain_correct` on DE440, `max_revs=seed.max_revs`)
→ emerged V∞ + DSM ΔV + emerged `n_revs_per_leg`
→ single-arc-degenerate guard → V1 anchor-match → [held for review].

## Success criterion (the golden cross-check)

The validated regression rows reconverge: a correctly-modeled V3/V1 `ch4` row closes
to `< 0.1 km/s` residual with its emerged V∞ within `V1_TOLERANCE_KMS` (0.5 km/s) of
its **sourced Russell anchor**, and its emerged `n_revs_per_leg` is > 0 on the
resonant legs / 0 on the transit leg. This proves the lane now represents the
geometry. Only after the regression rows pass is any statement made about a V0
promotion.

**S1L1 caveat:** the two rows that *are* S1L1 (`mcconaghy-2006-em-k2` =
`russell-ch4-4.991gG2`) may carry a mis-modeled `E-E-M-M` sequence (memory
`project_s1l1_nomenclature`: S/L are Earth-resonant intervals). They are therefore
NOT the regression test — the other validated `ch4` rows are. The S1L1 promotion
target stays held regardless of its convergence.

## Error handling

- Non-convergence after the sourced fix → recorded negative, stays V0 (a clean
  negative is success). If even the validated regression rows still fail to
  reconverge, **stop** — that is a deeper finding (genome / model mismatch beyond
  rev count), not a tolerance to loosen.
- A descriptor arc that does not reach the body → `seed_dsm_chain_from_descriptor`
  returns `None` (existing OFF-FAMILY-NO-CLOSE contract, unchanged).
- `hyperbolic_impossible` or DSM ΔV above a physical cap → loud flag, reject.

## Honesty gates (from the 2026-06-10 spec — unchanged, non-negotiable)

1. Golden V∞ = the published Russell-table cell; never self-computed.
2. single-arc-degenerate guard rejects rows a single ellipse already closes.
3. No catalogue writeback; any V0→V1 promotion is **held** for session review.
4. No tolerance/budget/cap loosening to force a closure.
5. n-body (V3) confirmation is a separate downstream step, not in this lane.

## Testing (TDD)

- **seed ToF:** `seed_dsm_chain_from_descriptor` returns `per_leg_tof_days` equal to
  the published arc ToFs (`[533.7, 150, 1026]` d for `russell-ch4-4.991gG2`), not the
  old `[883, 143, 883]` slack values.
- **rev cap:** `seed.max_revs == 2` for an `E-E-M-M` two-synodic row; `== 0` for a
  pure single-rev sub-Hohmann sequence; never exceeds the 6-body-period ceiling.
- **closer:** a correctly-modeled regression row converges and its emerged V∞ is
  within `V1_TOLERANCE_KMS` of the **sourced** anchor (golden), with
  `n_revs_per_leg > 0` on the resonant legs.
- **bounds:** the resonant-leg lower ToF bound brackets the published value (the
  corrector cannot reach the degenerate near-zero-ToF region).
- **negative preserved:** a NO-CLOSE / non-reaching descriptor still returns `None`;
  the single-arc-degenerate guard still rejects single-ellipse-closable rows.

## Out of scope (YAGNI)

- The 204 descriptor-less ocampo rows (publication gap; no seed possible).
- Out-of-plane / broken-plane DOF (coplanar seed → real-eph close, as #170).
- Pinning exact `rev_branch_per_leg` (rejected per Russell §5.4.2 — let it morph).
- Re-modeling the S1L1 `E-E-M-M` sequence (a separate, known issue —
  `project_s1l1_nomenclature` / `project_s1l1_realeph_closure_blocker`).
- Any catalogue writeback (the 2026-06-10 spec's Component 5, post-review only).

## References

- `docs/notes/2026-06-20-dsm-closure-batch-results.md` — the negative this fixes.
- `docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md` — `p.h.s.i`,
  6-body-period cap, §5.4.2 morph rationale.
- `docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md` — Appendix A bounds.
- `src/cyclerfinder/search/dsm_leg.py` — `dsm_leg(max_revs=…)`, `dsm_chain_correct`,
  `n_revs_per_leg` (the existing multi-rev machinery).
- `src/cyclerfinder/search/dsm_descriptor_seed.py` — the bridge being amended.
- Memory: `project_s1l1_nomenclature`, `project_s1l1_realeph_closure_blocker`,
  `golden-tests-sourced-only`, `orbit-closure-discipline`,
  `published-rounded-values-are-display`.
