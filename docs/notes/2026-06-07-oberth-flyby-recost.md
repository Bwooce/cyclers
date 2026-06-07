# Oberth periapsis powered-flyby re-cost — Aldrin maintenance ΔV (2026-06-07)

**Status: DIAGNOSTIC / PROVISIONAL. Single-solver. NOT for the site or the
catalogue.** Task #151. Adds an Oberth-credited periapsis powered-flyby cost
model alongside the existing asymptote-rotation surrogate, and re-costs the
classic Aldrin E-M maintenance ΔV under both.

Context read first (in order): the Takao 2025 MPGA-1DSM mine
(`docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md`, Eq. 11 — periapsis
P-FB ΔV with Oberth credit), the #148 primer-refine finding
(`docs/notes/2026-06-07-primer-refine-recoverable-dv.md` — the 2.9138 km/s is
*entirely* an Earth-flyby turn-deficit charge, heliocentric coasts near-
ballistic, so flyby geometry is the only lever), and the Russell 2004 mine
(`docs/notes/2026-06-07-russell-2004-dissertation-method-mining.md`, Eq. 5.5 —
the independent powered-SOI conditional Δv).

## What was built

* `core/flyby.py::dv_powered_flyby_periapsis(vinf, delta_required, delta_max,
  mu_planet, rp_min)` — Oberth-credited periapsis cost for the equal-|V∞| turn
  deficit. `dv_from_turn_deficit` is **untouched** (the published-comparison
  baseline and the regression anchor).
* `search/maintain.py` — opt-in `flyby_cost_model={"asymptote",
  "oberth_periapsis"}` threaded through `idealized_flyby_turn_deficit`,
  `optimise_maintenance_dv`, and `optimise_aldrin_maintenance_dv`. **Default
  unchanged = `"asymptote"`.** The optimiser objective (which pins the orbital
  anchors via `flyby_dv`/`dv_from_turn_deficit`) is unchanged, so the recovered
  cycler — a, e, t0, leg ToFs — is bit-identical across models; only the
  **reported** `maintenance_dv_kms` changes.

## The two models

Asymptote rotation (`dv_from_turn_deficit`, baseline) charges the deficit as a
rotation of the V∞ asymptote **at infinity** (speed V∞):

    ΔV_asym = 2·V∞·sin(deficit/2),   deficit = δ_req − δ_max.

Oberth periapsis (`dv_powered_flyby_periapsis`, new) supplies the residual turn
by a **tangential impulse at periapsis**, deep in the well. To deliver the full
δ_req ballistically the flyby must run at a lower excess speed `V∞_target` whose
cone opens to δ_req:

    sin(δ_req/2) = 1/(1 + rp·V∞_target²/μ)
      ⇒ V∞_target = √( (μ/rp)·(1/sin(δ_req/2) − 1) ).

Charge the Oberth speed difference (Takao Eq. 11 form) twice — slow into the
widened cone, restore |V∞| on the way out so the closure magnitude is preserved:

    ΔV_Oberth = 2·| √(V∞² + 2μ/rp) − √(V∞_target² + 2μ/rp) |.

## Mathematical sanity properties (label-mechanics, asserted in tests)

`tests/core/test_flyby_oberth.py`:

* (a) zero deficit → exactly 0.0 (within / at the cone edge). ✔
* (c) monotone non-decreasing in δ_req. ✔
* (b) Oberth < asymptote (strict) for nonzero deficit **in the deep-well
  regime** — asserted over a bounded sweep (Earth 200 km flyby, V∞ ≤ 6.5 km/s,
  deficits 1–40°; Mars, V∞ ≤ 2 km/s). ✔
* (d) units/frames consistent with `dv_from_turn_deficit` (km/s, km, km³/s²,
  rad). ✔

### Honest boundary of property (b): NOT universal

Property (b) is **not** a universal guarantee for these two specific models.
The Oberth periapsis maneuver beats asymptote rotation only when the well
dominates (`2μ/rp ≫ V∞²`). For an Earth 200 km flyby the crossover deficit (where
the two costs are equal) is:

| V∞ (km/s) | cone δ_max (deg) | crossover deficit (deg) |
|---|---|---|
| 1.0 | 159.3 | 18.7 |
| 3.0 | 121.1 | 56.9 |
| 5.0 | 90.1 | 87.9 |
| 6.0 | 77.7 | 100.3 |
| **6.86 (Aldrin)** | **68.5** | **109.5** |
| 8.0 | 58.2 | 0 (Oberth never wins) |
| 12.0 | 34.5 | 0 |

At high V∞ the ballistic cone is already so narrow that opening it to δ_req
demands a large magnitude excursion (×2), and the periapsis maneuver *exceeds*
the asymptote rotation. This is asserted explicitly
(`test_oberth_advantage_not_universal_high_vinf`, V∞=12) so the boundary is a
documented property, not a hidden regression. **The Aldrin V∞ ≈ 6.86 km/s sits
comfortably inside the favorable regime** (deficit 24.5° ≪ crossover 109.5°), so
the Oberth credit is real and robust for this schedule.

## Takao Eq. 11 vs Russell Eq. 5.5 — agreement verdict

**Verdict: AGREEMENT, with an important qualification on which case is which.**

Both papers split the powered flyby into a within-cone case and a deficit case,
and both place the cheap (Oberth) maneuver at periapsis and the turn-deficit
maneuver where magnitude/turn cannot be reconciled ballistically:

* **Within-cone (δ_req ≤ δ̄):** Takao Eq. 11 charges the periapsis Oberth speed
  change `|√(V∞in²+2μ/rp) − √(V∞out²+2μ/rp)|`. Russell Eq. 5.5 Case A charges
  `|v∞+ − v∞−|` *at the SOI* (zero-radius, no Oberth credit). These **agree in
  structure** (a magnitude change) but Takao does it at periapsis (Oberth-aware,
  strictly cheaper) and Russell at the SOI (conservative). For our equal-|V∞|
  closure this case costs ~0 (no magnitude change needed).
* **Turn deficit (δ_req > δ̄):** Takao Eq. 13 and Russell Eq. 5.5 Case B are the
  **same law of cosines on the residual turn** at the SOI:
  `√(v∞+² + v∞−² − 2 v∞+ v∞− cos(residual))`. With equal magnitudes both reduce
  **exactly** to `2·V∞·sin(residual/2)` — i.e. to our existing
  `dv_from_turn_deficit`. **Russell Eq. 5.5 Case B ≡ our asymptote baseline.**

The decisive finding: **neither published formula gives an Oberth credit for a
*pure equal-|V∞| turn deficit*.** Both deliberately maneuver at the SOI for the
deficit case (Russell explicitly excludes burns inside the hyperbola as too
risky). The Oberth credit in Takao Eq. 11 is for the *magnitude-change* (within-
cone) case. Our `dv_powered_flyby_periapsis` is therefore a **physically distinct
extension** — it converts the turn deficit into a periapsis magnitude excursion
(widen-cone-then-restore) so the Oberth lever applies — not a verbatim
transcription of either Eq. 11 or Eq. 5.5. It is consistent with their physics
(same hyperbolic mechanics, same μ/rp/V∞ relations) but is OUR construction, and
is labelled DIAGNOSTIC accordingly.

## Re-cost result — Aldrin E-M maintenance ΔV

In-family DE440 solve (V∞-anchor seeding via `real_window_priority_date =
1985-10-28`, the same path the BVP solver uses), `n_starts=5, seed=0`. Recovered
anchors (identical across models): a = 1.5878 AU, e = 0.3932, V∞_E = 6.86 km/s,
Earth turn 93.01° required vs 68.49° achievable (200 km flyby) → deficit 24.52°.

| model | reported maintenance ΔV (km/s) |
|---|---|
| **asymptote (baseline, regression anchor)** | **2.913836** |
| oberth_periapsis (OUR computation, DIAGNOSTIC) | **1.933552** |

Reduction: **0.9803 km/s (33.6 %)**.

* The asymptote value **reproduces ≈2.9138 km/s** to 1e-3 — the established
  regression anchor that cross-checks the `solve_powered_periodic_cycler` BVP
  path. ✔ (`test_aldrin_de440_recost_under_both_models`, slow.)
* The Oberth value **1.9336 km/s** is OUR computation. **Single-solver,
  PROVISIONAL.**

## Cross-check against the independent BVP solver — NOT AVAILABLE

The "independent" powered-cycler BVP solver
(`search/bvp.py::solve_powered_periodic_cycler`) **delegates** to
`optimise_aldrin_maintenance_dv` and charges the deficit through the *same*
`idealized_flyby_turn_deficit` → `dv_from_turn_deficit` path. It does **not**
carry an independent periapsis/Oberth cost model. Therefore:

* The 2.9138 asymptote value IS cross-checked by the BVP path (same code, the
  established evidence chain).
* The 1.9336 Oberth value is **single-solver** — there is no independent solver
  that supports a periapsis model to confirm it. **Honestly: provisional, NOT
  site-ready.**

## Plausibility gate

Both values pass `verify/plausibility.check_publishable(MAINTENANCE_DV_KMS, …)`
(engineering bar 3.0 km/s): asymptote 2.914 (within bar), Oberth 1.934 (within
bar, more margin). ✔ (Asserted in the slow test.)

## Site-update recommendation

The Oberth re-cost lowers the headline Aldrin maintenance ΔV by ~34 % (2.914 →
1.934 km/s) and is physically better-motivated (the maneuver is paid deep in the
well where it is cheapest, which is what a real mission would do). **However it
is single-solver and DIAGNOSTIC.** Recommendation to the main session:

1. Do **not** update the site/catalogue/windows.json on this task (per the task
   rules — left for the main session to decide).
2. Before publishing the lower number, obtain an **independent confirmation** —
   either a second solver carrying a genuine periapsis cost model, or a direct
   integration of the powered flyby — since the BVP path cannot cross-check it.
3. If the main session does update the site, keep the asymptote value as the
   documented conservative upper bound and present the Oberth value as the
   Oberth-optimal lower bound, with this note's caveats.

## Files touched

* `src/cyclerfinder/core/flyby.py` — `dv_powered_flyby_periapsis` (new).
* `src/cyclerfinder/search/maintain.py` — `FlybyCostModel` selector threaded
  (default `"asymptote"`).
* `tests/core/test_flyby_oberth.py` — math sanity properties (new, 12 tests).
* `tests/search/test_maintain.py` — model-selector + slow DE440 re-cost tests.
* `docs/notes/2026-06-07-oberth-flyby-recost.md` — this note.
