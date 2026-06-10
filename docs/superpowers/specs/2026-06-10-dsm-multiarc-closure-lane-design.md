# DSM Multi-Arc Closure Lane — Design Spec

**Date:** 2026-06-10
**Status:** Design approved (brainstorming), pending implementation plan.
**Sub-project:** 2 of 2 in the "genome upgrade" track. Spec 1 (VILM-leveraging
endgame solver, #179) is complete. This spec is the catalogue-lift sub-project.

## Goal

Wire the existing Takao one-DSM-per-leg genome (`dsm_chain_correct`) into a
catalogue-row **closure-and-validation lane**, so the descriptor-bearing
off-family Russell rows — which close geometrically but whose single-ellipse
genome has no sourced basin — can be re-derived **on-family** as two-arc
trajectories and promoted off V0 where the closure holds.

## Scope reality (read first — this bounds the payoff)

The #177 triage of all 212 unsourced rows
(`docs/notes/2026-06-08-self-seeding-triage-results.md`) found:

| outcome | count | meaning |
|---|---|---|
| REACHABLE | 6 | a coplanar G-arc branch transit lands within 30 d of the tabulated transit |
| OFF-FAMILY-NO-CLOSE | 2 | descriptor exists but the shape does not reach Mars (`5.30gGf3`, `5.75ggF3`) |
| OFF-FAMILY-NO-DESCRIPTOR | 204 | `russell-ocampo` rows: no per-arc g/G descriptor (the `n.m.k` summary format) |

**A multi-arc genome needs a descriptor to seed from. The 204 ocampo rows have
none** — they are *descriptor-gated* (the publication gap, see memory
`validation-ceiling`), NOT genome-gated. No genome closes a row it cannot seed.
So this lane's entire reachable payoff is the **8 descriptor-bearing rows**
(6 REACHABLE + 2 NO-CLOSE). The 204 are explicitly OUT OF SCOPE and their
exclusion is a recorded fact, not a failure of this work. The build's integration
is reusable; the catalogue lift is bounded to ≤ 8 rows.

## Why this, why now

The 6 REACHABLE rows carry a 2-arc g/G free-return descriptor and a tabulated
transit that a coplanar G-arc branch reproduces, but they sit at V0 because a
single-ellipse genome has no residual-zero basin at their sourced V∞ anchors —
they are genuinely two generic-return arcs (the S1L1 lesson, #167). The Takao
η-DSM genome (`dsm_leg.py`, #150) can represent a leg that follows a different
ballistic arc on its front fraction than its back — exactly the two-arc geometry.
It already exists and is tested (used by `free_return_chain`, `continuation_chain`,
`self_seeding`); it is simply not wired into a catalogue-row validation pipeline.
No paper acquisition is required — pure engineering.

## Binding constraints (orbit-closure-discipline, non-negotiable)

Baked into the gates, not optional:
1. **Same-model sourced golden:** the EXPECTED V∞ is each row's **published
   Russell-table summary cell**, never a value our own code computed. The DSM
   closure EMERGES V∞; it is never imposed (golden-rule separation).
2. **Independent cross-check:** lamberthub izzo2015 vs gooding1990 per-leg
   agreement (two Lambert methods — the spec §14 V1 mechanism), PLUS the
   single-arc-degenerate guard (a different genome as second opinion).
3. **single-arc-degenerate guard:** if a plain single-ellipse already closes the
   row to V1, the DSM result is REJECTED — no manufacturing a "multi-arc" row
   that did not need the extra DOF.
4. **n-body confirmation** before any V3 claim.
5. **No catalogue writeback** until n-body passes AND the result is
   session-reviewed (held, as in #170).
6. **No tolerance/budget/cap loosening** to force a closure. A row that will not
   close on-family is a recorded NEGATIVE (stays V0) — a clean negative is a
   success.

## Architecture

A DSM closure lane reusing `dsm_chain_correct`, mirroring the #170 App-C batch.
Four new focused units + a separate held writeback.

### Component 1 — `src/cyclerfinder/search/dsm_descriptor_seed.py`

```python
@dataclass(frozen=True)
class DsmChainSeed:
    sequence: tuple[str, ...]
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    eta_init: tuple[float, ...]      # 0.0 per leg — starts ballistic
    arc_a_au: float                  # coplanar descriptor arc shape (a, e)
    arc_e: float
    transit_branch: str              # short | long | k-rev (from g_arc_branches)
    vinf_anchor_kms: float           # the row's sourced Russell-table V∞ cell

def seed_dsm_chain_from_descriptor(row: dict) -> DsmChainSeed | None
```
Parse the row's 2-arc g/G descriptor; reuse `self_seeding.g_arc_branches` to get
the coplanar arc `(a, e)` and the transit branch matching the tabulated transit;
assemble the `dsm_chain_correct` inputs and the per-leg bounds via
`dsm_leg.sequence_keyed_bounds`. Returns `None` for a row with no usable
descriptor (the 204 ocampo rows — they never enter the lane).

### Component 2 — per-row closer (in `dsm_descriptor_seed.py` or a sibling)

```python
@dataclass(frozen=True)
class DsmClosureResult:
    converged: bool
    max_residual_kms: float          # CONSTRAINED
    dsm_dv_kms: tuple[float, ...]    # CONSTRAINED — per-leg interior impulse
    vinf_per_encounter_kms: tuple[float, ...]  # EMERGED
    vinf_anchor_kms: float
    anchor_match: bool               # |emerged - anchor| <= V1 tol
    hyperbolic_impossible: bool      # loud physics flag (never a filter)

def close_row_dsm(row: dict, ephem: Ephemeris) -> DsmClosureResult
```
Seeds via Component 1, runs `dsm_chain_correct` on the real DE440 ephemeris,
extracts the emerged per-encounter V∞ and the η interior-impulse ΔV. Converged by
residual magnitude (the corrector's own criterion). Never raises — an infeasible
seed returns `converged=False`.

### Component 3 — validation gates (full ladder, capped per-row by data)

In order, each blocking the next:
1. **single-arc-degenerate guard** — call `free_return_chain.single_arc_degenerate`
   (or the single-ellipse corrector); if it already closes the row to V1, REJECT.
2. **V1** — `anchor_match` (emerged V∞ within `V1_TOLERANCE_MPS` of the Russell
   cell) AND lamberthub izzo+gooding per-leg agreement `< V1_TOLERANCE_MPS` AND
   Kepler forward-reprop residual pass. Reuse `verify/crosscheck`.
3. **V3** — n-body horizon-TCM over 3–5 laps `<=` the §14 V3 threshold. Reuse the
   REBOUND harness (the moon-tour/heliocentric variant as appropriate).

### Component 4 — batch driver `scripts/dsm_closure_batch.py`

Iterate the 8 descriptor-bearing rows; per-row seed→close→guard→V1→V3; write a
results note (`docs/notes/YYYY-MM-DD-dsm-closure-batch-results.md`) + a runlog
(`data/runs/dsm-closure-<ts>.jsonl`) with per-row outcome and, for passes, the
**proposed** `_LEVEL_EVIDENCE` text. NO catalogue writeback.
`export PATH="$HOME/.local/bin:$PATH"`; cap workers; honest EMPTY/negative rows.

### Component 5 — writeback (separate, post-review)

After session review of the batch results: add the confirmed `_LEVEL_EVIDENCE`
entries to `src/cyclerfinder/data/validate.py` and bump `validation_level` in
`data/catalogue.yaml` for each passing row. Per-row, sourced-evidence-tagged,
exactly as #168/#170. Not part of the batch run.

## Data flow

row (descriptor + Russell V∞ anchor)
→ `seed_dsm_chain_from_descriptor`
→ `close_row_dsm` (`dsm_chain_correct` on DE440) → emerged V∞ + DSM ΔV
→ single-arc-degenerate guard
→ V1 gate (anchor match + izzo/gooding + Kepler reprop)
→ V3 gate (n-body horizon-TCM)
→ runlog + proposed `_LEVEL_EVIDENCE`
→ [session review]
→ writeback (validation_level + evidence registry).

## Error handling

- Non-convergence / infeasible seed → `converged=False` → recorded negative, no
  promotion.
- The 2 NO-CLOSE rows: the interior impulse MAY bend the arc to reach Mars; if it
  cannot, a legitimate recorded negative (stays off-family).
- single-arc-degenerate → reject (would be a false multi-arc claim).
- `hyperbolic_impossible` or DSM ΔV above a physical cap → loud flag, reject.
- A row with no descriptor → `seed_dsm_chain_from_descriptor` returns `None`; the
  batch records it OUT-OF-SCOPE (the 204; should never reach the closer).

## Testing (TDD)

- **adapter:** `seed_dsm_chain_from_descriptor` builds the expected seed
  (sequence/arcs/bounds/anchor) for `russell-ch4-9.353Gg2`; returns `None` for a
  descriptor-less ocampo row.
- **closer:** a REACHABLE row converges and its emerged V∞ is within
  `V1_TOLERANCE_MPS` of the **Russell-table anchor** (golden — sourced, not
  self-computed).
- **single-arc-degenerate guard:** a row that the single-ellipse genome already
  closes is rejected from the DSM lane.
- **V1 crosscheck:** izzo/gooding per-leg agreement gate reused and asserted.
- **n-body confirm:** marked `@pytest.mark.slow`.
- **negative:** a NO-CLOSE row → `converged=False` or `anchor_match=False` →
  recorded negative, no promotion (no fabricated closure).

## Out of scope (YAGNI / deferred)

- The **204 descriptor-less ocampo rows** — no seed possible; descriptor-gated
  publication gap. Recorded as out-of-scope, not retried.
- Real-eph **inclination / broken-plane** — the seed is the coplanar descriptor
  arc; the corrector closes on DE440 exactly as #170 (coplanar seed → real-eph
  close). Out-of-plane DOF is a later rung.
- DSM as a `residual_mode` inside `ballistic_correct` (Approach B) and grafting
  DSM onto `free_return_chain` (Approach C) — both rejected for redundancy /
  regression risk.

## References

- `docs/notes/2026-06-08-self-seeding-triage-results.md` — the 6/2/204 triage.
- `src/cyclerfinder/search/dsm_leg.py` — `dsm_chain_correct`,
  `make_dsm_chain_step`, `sequence_keyed_bounds` (the Takao genome).
- `src/cyclerfinder/search/self_seeding.py` — `g_arc_branches` (descriptor → arc).
- `src/cyclerfinder/search/appc_corrected.py` — the #170 per-row close precedent.
- `src/cyclerfinder/data/validate.py` — `_LEVEL_EVIDENCE` over-claim registry.
- Memory: `validation-ceiling`, `orbit-closure-discipline`,
  `golden-tests-sourced-only`.
