# #318 Phase 2 — real-ephemeris n-body multi-axis joint search (DESIGN / START)

**Date:** 2026-06-30. Status: Phase-2 STARTED — this is the design for the unified build the
Phase-2 reframe called for (`docs/notes/2026-06-30-318-phase2-blocker-diagnosis-reframe.md`).
The Phase-1 substrate's "joint sweep" can never be jointly 3D because axes A/B/D (heliocentric
Lambert) and axis C (CR3BP) don't share a state model. The fix: do the joint search in the
ONE model where all four axes co-vary naturally — the **real-ephemeris n-body trajectory**.

## The unified joint-cell
A candidate is a single real-eph n-body trajectory; the four axes are its free dimensions:

| Axis | In the real-eph n-body cell |
|---|---|
| **C — 3D / broken-plane** | INTRINSIC: the 6-D state has out-of-plane components; no CR3BP needed |
| **D — epoch-locked** | the anchor/launch epoch sets the real ephemeris geometry (the dominant phasing knob, per the #480 EIGE work) |
| **A — powered maintenance** | per-flyby maneuver ΔV (the `chain_cycles` per-cycle re-target pattern, #223-validated) |
| **B — multi-rev** | the Lambert branch / rev-count of each leg's SEED (then shot in n-body) |

So a cell = `(sequence, epoch, 3D-state-seed, per-leg (n_rev,branch), ΔV-budget)` → build the
seed → shoot in real-eph n-body → return (closure residual, total ΔV, feasibility, V∞ regime).

## Reuse (do NOT rebuild)
- `nbody/shooter.py::shoot` (heliocentric n-body multiple-shooter) + `nbody/jovian.py`
  (`jovian_shoot`, `chain_cycles`, `JovianRestrictedNBody`) — the real-eph propagate+correct.
- `nbody/jovian_stm.py` — the analytic STM Jacobian (block-bidiagonal, ~40× faster than FD).
- `core/lambert.py` (multi-rev seeds, Axis B) + the real-eph ephemerides (DE440 / JUP365 / URA111).
- `search/literature_check.py` for the post-hoc novelty gate (V∞ + sequence fingerprint).

## The crux risk — COMPUTE (carry the #480 lesson)
[[project_dsm_closure_modeljump_blocker]]: real-eph n-body shooting is EXPENSIVE — a single
multi-rev multi-year cycler shoot is >400 s with FD; the analytic STM helps but the corrector
can PLATEAU (the #480 EGGIE/EIGE level-3 wall, `2026-06-30-480-level3-approach-c-verdict.md`).
A naive Sobol sweep with an n-body shoot per cell is compute-infeasible at scale. The design
MUST be compute-aware:
1. **Short cyclers only** (≤ ~1-2 revs, ≤ ~weeks): the #480 EIGE 1-rev case shot in <1 s; the
   multi-year EM/EGGIE cases did not. Restrict the joint search to short-period moon-tour /
   inner-system cyclers where the shoot is cheap.
2. **Analytic-STM corrector** (`jovian_stm`) not FD.
3. **Patched-conic pre-filter** (the existing `close_powered_cycle` / Lambert close) as a cheap
   surrogate to rank cells BEFORE any n-body shoot — only shoot the top-K survivors.
4. **Small Sobol budget** over the joint dims (not Cartesian; Phase-1 showed Cartesian wastes 96%).

## Build order
- **Phase 2a (the bounded first increment — next):** the joint-cell EVALUATOR + a POSITIVE
  CONTROL. Wrap (build-seed → patched-conic surrogate → optional n-body shoot) into one
  `evaluate_joint_cell(...)` returning the verdict tuple. Positive control: reproduce a
  catalogued real-eph member's closure as a joint-cell (e.g. the #312 Uranus V4 quasi-cycler,
  or a short Jovian moon-tour) — the evaluator must re-find it BEFORE any search is trusted
  ([[feedback_verify_gauntlet_with_positive_control]]).
- **Phase 2b:** the compute-aware Sobol search on a target short-cycler system; post-hoc
  lit-novelty on any survivor.

## Honest prior + stop criterion
Strong-prior-EMPTY: Phase-1 EM collapsed to 2 distinct dV (the axes were redundant), and every
2026-06-30 discovery probe (#465, #320 Saturn/Pluto/QP-tori) resolved to V0-known. The realistic
Phase-2 outcome is a **compute-bounded empty-region map**, not a discovery — log what was
searched + the compute cap honestly ([[project_negative_results_registry]]); a clean negative on
this genuinely-novel joint manifold is a legitimate result. STOP if the positive control can't
be reproduced cheaply (the compute wall bites before the search even starts) — that itself is
the finding (the joint manifold is real-eph-n-body-compute-gated, needs a surrogate not a shoot).

## Scope honesty
Multi-week. Phase 2a (evaluator + positive control) is the bounded next deliverable; 2b (search)
follows only if 2a's positive control reproduces cheaply. No catalogue writeback; no novelty
claim until `literature_check` clears a survivor against the published record.
