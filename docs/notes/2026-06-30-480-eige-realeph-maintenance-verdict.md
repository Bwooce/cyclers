# #480 — EIGE real-ephemeris maintenance positive control: CHARACTERIZED NEGATIVE

**Date:** 2026-06-30. Status: the real-ephemeris EIGE maintenance-ΔV positive control
**cannot be completed with the patched-conic chain lane** — the same 3-D B-plane wall
documented for EGGIE, but TIGHTER (EIGE yields zero feasible real-eph members where
EGGIE yielded nine). The ideal-model EIGE construction (the prior commit) stands; this
is the real-ephemeris side.

## Goal
Reproduce the paper's Fig-5 EIGE maintenance example — ballistic first cycle, then
~30 m/s over 10 repeat cycles (AAS 17-608 pp.10-11) — as the positive control for the
per-cycle maintenance method (`nbody.jovian.chain_cycles`) before quoting an EGGIE number.

## What was built + run (all foreground, instrumented)
1. **Maintenance lane driver** `scripts/eige_maintenance_480.py`: drives `chain_cycles`
   on `sequence=EIGE_SEQUENCE`, 3-leg single-rev branch plan, seeded by
   `eige_ballistic.eige_tof_seed_days()` in the real DE440/JUP365 ephemeris. Works given
   a feasible seed; reusable for the EGGIE Approach-A lane.
2. **Phase-match epoch search**: `optimize_cycle` only tweaks ToFs (±bound), so a blind
   epoch scan never lands the EIGE moon configuration. The ideal member fixes the
   RELATIVE phases (Io−Europa = φ_Io = −0.194 rad, Ganymede−Europa = φ_Gan = −3.876 rad).
   Scanning 2020 for those, the best real-ephemeris alignment is **~0.65 rad (≈37°) off**,
   recurring every ~7.05 d (Ganymede synodic) and never closer over the window — the real
   eccentric/inclined moons do not reach the ideal circular-coplanar EIGE configuration.
3. **Real-eph ballistic corrector** over free (epoch, 3 ToFs): DOES find ballistic
   closures (equal-in/out |V∞| to ~1e-10, total flyby ΔV = 0) in the correct
   low-excess-speed regime (V∞ Europa 6.3-9.0, Io 5.2-6.0, Ganymede 5.9-7.1 km/s) — but
   **every closure has the Io flyby ~−1746 km and Europa ~−1517 km (deeply sub-surface,
   near-180° turns)**.
4. **EGGIE-style feasibility-first discovery** (minimize total ballistic flyby ΔV +
   altitude-window penalty, V∞ an output — the exact objective that found 9 feasible
   EGGIE members, `2026-06-30-480-eggie-realeph-unguided-discovery.md`): from 40
   phase-match seeds across 2020-21, **0 feasible EIGE members**. The optimiser drives ΔV
   to 0 (ballistic) but cannot lift the Io/Europa flybys above surface.

## Why EIGE fails where EGGIE succeeded (the physical reason)
EGGIE is 4-synodic / 5-rev / 4-leg with a repeated Ganymede — enough geometric slack to
arrange feasible flybys; the unguided search found 9 feasible members. **EIGE is
1-synodic / 1-rev / 3-leg** — the single revolution pins the Io and Europa crossings to
points where the required turn is near 180°, i.e. sub-surface, in a patched-conic Lambert
reconstruction (which fixes each leg's V∞ DIRECTION and so exposes no independent B-plane
orientation freedom). The paper's feasible Fig-5 EIGE (Io 2,817 / Europa 470 km) uses its
full 3-D B-plane + maneuver formulation (Eqs 6-7) and Monte-Carlo/NLP optimisation — the
DOF our chain lane lacks. This is the **same wall** as the EGGIE Table-4 member
([[project_dsm_closure_modeljump_blocker]] leading cause: "coplanar reconstruction deletes
the paper's 3-D B-plane flyby orientation DOF").

The ideal coplanar model DID admit a feasible ballistic EIGE (prior commit: Io 2,817 /
Ganymede 13,180 / Europa 1,323 km) because there the circular moons can be placed exactly
on the conic crossings; the ~37° real-ephemeris phase residual is what pushes the real-eph
flybys sub-surface.

## Verdict / honesty
- **No catalogue impact.** No exact-member reproduction claim (consistent with the EGGIE
  verdict). The EIGE V∞ figures are reported construction outputs (the paper prints none).
- **The maintenance METHOD remains validated** via the Liang Member D chained ΔV (#223) —
  the scope doc's stated alternative positive control. The EIGE-specific control is the
  nice-to-have that hits this wall.
- **Not attempted:** the paper's full 3-D B-plane / SNOPT-style NLP (scope Approach C,
  "weeks, last resort"). The 1-rev EIGE infeasibility in the patched-conic lane is the
  characterized boundary; breaking it needs the B-plane-DOF optimiser, not a re-run.

## Bankable
- `scripts/eige_maintenance_480.py` — EIGE/arbitrary-sequence maintenance lane (reusable).
- Ideal EIGE ballistic construction `search/eige_ballistic.py` (committed, golden-gated).
- Lesson reinforced: a too-narrow search FORMULATION (pure equal-|V∞| residual, no
  feasibility) hid the result; the EGGIE-style feasibility-first objective is the right one
  — and even it returns a clean 0 feasible for EIGE ([[feedback_verify_gauntlet_with_positive_control]]).
