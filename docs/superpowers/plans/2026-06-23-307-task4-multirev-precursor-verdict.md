# #307 Task 4 — multi-rev precursor re-run verdict (Aldrin / S1L1)

**Date:** 2026-06-23
**Status:** COMPLETE — honest negative (the model wall holds).

## What was run

The #289-Phase-5 substrate (Tasks 1–3) added three capabilities to the
epoch-locked precursor matcher:

1. multi-rev Lambert branch selection in `close_epoch_locked` (`max_revs`),
2. automated Vasile–Conway DSM placement (standalone module),
3. cycler-cadence terminal-phase-window auto-derivation.

Task 4 wired `max_revs` through the full matcher closure stack
(`find_cycler_precursors → optimise_chain_tofs → _epoch_locked_loss →
close_epoch_locked`, commit `3775579`) and re-ran the #302 Aldrin + S1L1
precursor scans with `max_revs=2` as the **only** changed variable vs the
single-rev #302 baseline (`scripts/run_307_precursor_multirev.py`). Identical
launch window (2030–2034), V∞ grid, TOF box (80–500 d), gates, and
400-candidate validation cap.

Outputs: `data/precursor_307_{aldrin,s1l1}_multirev.jsonl` (394 survivors each).

## Result

| target | metric | #302 single-rev | #307 multi-rev (max_revs=2) |
|---|---|---|---|
| Aldrin | min closure residual | 0.4547 km/s | 0.4547 km/s |
| Aldrin | min flyby continuity | 4.1549 km/s | 4.1549 km/s |
| Aldrin | median flyby continuity | 17.48 km/s | **15.40 km/s** |
| S1L1 | min closure residual | 1.7916 km/s | 1.7916 km/s |
| S1L1 | min flyby continuity | 3.9535 km/s | 3.9535 km/s |
| S1L1 | median flyby continuity | 17.59 km/s | **15.40 km/s** |

Ballistic gate = 0.10 km/s flyby continuity. **Candidates passing the gate:
0 / 394 in both runs. Candidates even under 1.0 km/s continuity: 0 / 394.**
Literature-fresh survivors: 0 (every survivor matches a published anchor).

## Reading

- **`max_revs` is genuinely active.** The validated distributions are *not*
  identical (`identical_dist=False`); the median flyby-continuity drops ~2 km/s
  on both targets, i.e. multi-rev branch selection does find different,
  lower-continuity heliocentric arcs for many candidates.

- **But the wall does not fall — it moves between constraints.** The two binding
  residuals never co-minimise on one trajectory. The best *joint* candidate
  (minimising `max(closure, continuity)`) is an E–V–E chain with
  closure ≈ 3.05 km/s **and** flyby-continuity ≈ 4.0–4.2 km/s on both targets —
  both an order of magnitude above gate. A leg can be driven near
  ballistic-continuous (the smoke probe hit `flyby_cont ≈ 0`) only at the cost
  of a large terminal V∞-mismatch (~4 km/s); pushing one residual down lifts the
  other.

- **This is the publication-gap result, re-confirmed with a stronger tool.**
  Ballistic precursor-MGA insertion into the published classic Earth–Mars
  cyclers is not reachable in the coplanar/circular-body epoch-locked model,
  even with multi-rev Lambert branches. The negative is consistent with the
  #302 Phase-4 note (the precursor-MGA space into these cyclers is itself
  published) and with the #388/S1L1 family-selection findings (the real geometry
  needs body eccentricity — the model, not the search, is the limiter).

## Disposition

- `aldrin-classic-em-k1-outbound` and `s1l1-2syn-em-cpom` precursor scans stay
  **V0**. No promotion; no novelty claim (0 literature-fresh).
- The next lever is *not* more branches but a more faithful model: an
  eccentric-body Tisserand enumerator and/or DSM post-repair on the
  near-continuous arcs to absorb the residual V∞-mismatch (Task 2's standalone
  optimiser is the building block, deliberately left un-wired here — wiring it is
  a separate, sourced build, not part of the #307 plan scope).
- #307 Tasks 1–4 complete. The eccentric-Earth Tisserand named in the issue
  title was never decomposed into a plan task and is logged here as the
  motivated follow-on.
