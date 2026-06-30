# #465 — at-scale powered moon-tour campaign: VERDICT (capability-breadth, all V0-known)

**Date:** 2026-06-30. The explicit follow-on to the multi-rev leveraging releg verdict
(`docs/superpowers/plans/2026-06-26-465-multirev-leveraging-verdict.md`): sweep moon-tour
skeletons × ToF-scale × phasing with `MultiRevLeveragingReleg` across the Galilean +
Saturnian systems. **Result: the chain closes in-band for 8 distinct skeletons across 2
systems — and all 8 are V0-KNOWN reproductions of the published resonant-leveraging
endgame/capture-tour family.** No novel cyclers; no catalogue self-admission.

## The sweep (`scripts/releg_moontour_campaign_465.py`)
160 combos = (Jovian Io/Europa/Ganymede/Callisto + Saturnian Enceladus/Tethys/Dione/Rhea/
Titan) contiguous 3- and 4-moon cyclic skeletons × ToF-scale {1.0,1.2,1.5,2.0,2.5} ×
phase-step {0.5,1.0,1.5,2.0}. Each run = one `close_powered_cycle`. Results:
`data/releg_moontour_campaign_465.jsonl`. **71 in-band runs / 89 honest chain-stalls**
(the reachability ceiling V∞≲~2·V_M biting at bad phasing/scale — no fabricated bridges).

## The 8 distinct in-band closures (best ΔV per skeleton; continuity verified)
| system | skeleton | ΔV/cycle (km/s) | continuity (km/s) |
|---|---|---|---|
| Saturn | Dione-Rhea-Titan-Dione | 0.070 | 1.3e-13 |
| Saturn | Tethys-Dione-Rhea-Tethys | 0.091 | 6.6e-14 |
| Jupiter | Europa-Ganymede-Callisto-Europa | 0.163 | 6.8e-4 |
| Saturn | Enceladus-Tethys-Dione-Enceladus | 0.469 | 1.0e-13 |
| Jupiter | Io-Europa-Ganymede-Callisto-Io (4-moon) | 0.570 | 1.2e-13 |
| Jupiter | Io-Europa-Ganymede-Io (the #465 base case) | 0.671 | 5.8e-4 |
| Saturn | Enceladus-Tethys-Dione-Rhea-Enceladus (4-moon) | 0.737 | 2.2e-13 |
| Saturn | Tethys-Dione-Rhea-Titan-Tethys | 0.829 | 1.0e-13 |

All continuity-closed (V∞-continuous by construction), all inside the 3.5 km/s/cycle
powered ceiling, all `powered_dsm` band. "It closed" verified non-degenerate
([[feedback_orbit_closure_discipline]]): continuity residual ≤ 7e-4 km/s on every hit.

## Lit-novelty (the mandatory baseline — all V0-known, lit-grounded)
Per [[feedback_literature_novelty_check_baseline]], grounded against the published record
(not asserted):
- **Jovian** Io-Europa-Ganymede(-Callisto) resonant-leveraging tours are published — the
  GTOC6 winning solution (Izzo et al., resonance-hopping Galilean tour) and Europa-Orbiter
  Io-gravity-assist tours; the #465 decomposition golden already anchors to Campagnola-
  Russell "The Endgame Problem" Part 1.
- **Saturnian** Enceladus-Tethys-Dione-Rhea-Titan leveraging tours are the **Strange-
  Campagnola-Russell** capture tour verbatim ("a leveraging tour... combining gravity assists
  from Titan, Rhea, Dione, and Tethys with VILMs" to capture at Enceladus) — exactly the
  campaign's Saturnian skeletons; the #465 Rhea-Dione decomposition golden anchors to the
  same paper (Part 2).
- **Verdict: all 8 skeletons are V0-KNOWN** reproductions of the published resonant-
  leveraging endgame/capture-tour family. None is SILVER-novel. The honest expectation of
  the #465 verdict ("V0-known, a reproduction of a published tour") is confirmed at breadth.

## Standing / discipline
- **Capability-breadth result, not a discovery.** The campaign confirms the multi-rev
  leveraging chain reproduces the published moon-tour family across 8 skeletons in 2 systems
  (up from the 1 demonstrated in the #465 verdict), and characterizes the reachable region
  (8 closing skeletons) vs the honest reachability stalls (89 runs).
- **No catalogue self-admission.** The 8 are V0-known candidates flagged for the human
  gauntlet (a human-admitted reproduction is V0-known at most, never self-admitted). The
  per-skeleton citation pinning (which exact published tour each reproduces) + the V2
  moontour gauntlet is the human step.
- **No empty-region re-stamp change** beyond #465's existing `multi-rev-leveraging` stamp:
  the 89 stalls are phasing/ToF reachability misses within bridgeable contours, not method-
  level empties (the prefilter found 0 disjoint-contour skips in these contiguous rings).
