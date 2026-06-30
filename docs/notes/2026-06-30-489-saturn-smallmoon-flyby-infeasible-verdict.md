# #489 — small-moon-flyby Saturn #320 candidates: physically INFEASIBLE (R-S premise confirmed)

**Date:** 2026-06-30. Re-evaluated the #320 Saturn two-moon candidates whose novelty was
REOPENED earlier today after acquiring Russell-Strange 2009 (its Saturnian census is
Titan→Enceladus only, so the small-moon-flyby pairs aren't covered). **Verdict: they are
physically INFEASIBLE — the small moons cannot provide the gravity assist — which CONFIRMS
Russell-Strange's premise rather than reopening novelty.** No novel Saturn candidate; no new row.

## The check (the #320 data already recorded it; grounded, not re-derived)
A repeated-moon cycler X-Y-X uses X as the gravity-assist (flyby) body. `scan_320_epoch_aware_saturn.jsonl`
carries `physical_gate_passed` + `max_bend_deg_per_enc` (the #324 max-ballistic-deflection gate):

| candidate | flyby body | μ (km³/s²) | best residual | flyby max-bend | #324 gate |
|---|---|---|---|---|---|
| **Titan-Rhea-Titan** (the SILVER) | Titan | 8978 | 0.0316 | **49.4° / 50.5°** | **PASS** |
| Tethys-Enceladus-Tethys | Tethys | 41 | 0.0258 | **0.44° / 0.44°** | FAIL |
| Dione-Tethys-Dione | Dione | 73 | 0.0387 | **3.1° / 3.2°** | FAIL |

EVERY non-Titan-flyby Saturn candidate has `physical_gate_passed: FALSE` (flyby max-bend
0.07–5°). Only the Titan-flyby cycler clears the gate (Titan turns ~50°). The low residual on
the failed ones is only V∞-MAGNITUDE continuity — the small moon physically cannot REDIRECT the
trajectory, so it is not a functional flyby (a "closure" with a sub-degree turn is geometric,
not dynamical).

## My low-V∞ hypothesis was WRONG (the honest correction of the correction)
I'd flagged "low V∞ (1–4 km/s) → bend ∝ 1/V∞² is large → small moons might be feasible." Wrong:
even at V∞ 2–4 km/s, Dione/Tethys (μ 41–73, ~200× lighter than Titan) give only ~0.4–5°. The
mass deficit dominates the V∞ factor. **Russell-Strange's "only Titan can provide the Saturnian
gravity assists" premise holds**, and the #320 physical-sanity gate already encoded it
correctly — the acquisition confirmed the gate's physics, it did not overturn it.

## The twice-corrected arc (landed right)
1. Morning #320 verdict: Saturn = V0-known (cited R-S without holding it).
2. After acquiring R-S: corrected — R-S only does Titan-Enceladus, so small-moon pairs "reopened."
3. **#489 (this):** the small-moon pairs are physically INFEASIBLE (the gate the #320 pipeline
   already applied) → not novel, not even valid cyclers. **Titan-Rhea-Titan stays V0-known**
   (R-S Titan-flyby cycler targeting Rhea, a documented R-S target moon).

## Standing
The #320 Saturn thread is fully closed: no novel candidate. The lit-check is moot (an
infeasible "cycler" is not a candidate regardless of novelty). The #320 pipeline's #324
physical-sanity gate is validated against Russell-Strange's independent physical argument. No
catalogue change. Net of the whole #320 discovery arc: the only genuinely-novel admitted hit
remains the #312 Uranus quasi-cycler (V4-catalogued); everything else is V0-known, infeasible,
or a published class.
