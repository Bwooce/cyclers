# Oberth periapsis flyby cost — independent confirmation (2026-06-07)

**VERDICT: CONFIRMED.** Task #154. The Oberth-credited periapsis powered-flyby
cost model `core/flyby.py::dv_powered_flyby_periapsis` (task #151) is confirmed
by direct numerical integration of the two-body hyperbola — a method fully
independent of the analytic formula under test. At the Aldrin Earth-flyby point
the integrated Δv matches the formula to **1 part in 1e9**.

## What was checked and how (option a — direct integration)

New module `src/cyclerfinder/verify/flyby_integrate.py`. The formula computes
the periapsis maneuver from the closed-form hyperbolic bend relation. The check
instead **numerically integrates** the planar two-body equations of motion (the
flyby body's `mu` only) with `scipy.integrate.solve_ivp` (DOP853, rtol 1e-12),
extracts the realised asymptote-to-asymptote turn by propagating to a far field
(default 1e5·rp) up- and down-stream, and **root-solves** (Brent) over the
periapsis speed for the maneuver that delivers the required turn. The reported
Δv is `2·|vp_in − vp_target_integrated|` — a function only of integrated
trajectories and a 1-D root-find. It never calls `dv_powered_flyby_periapsis`
inside the integrated path (the formula is invoked only to populate the
comparison field). Same physics, different method — the SPICE-cross-check
pattern.

The physical model being confirmed (note `2026-06-07-oberth-flyby-recost.md`):
the residual turn is supplied at periapsis by slowing into a widened ballistic
cone, then restoring `|V∞|` on the way out — two tangential periapsis impulses,
`ΔV = 2·|vp_in − vp_target|`. A slower hyperbola bends more, so a turn deficit
(required > cone) is met by `vp_target < vp_in`, found here by integration.

## Harness self-validation (task §3 — validates the check before it judges)

* **Integrated cone vs analytic `max_bend`:** integrating the *unpowered*
  incoming hyperbola and reading its turn reproduces `max_bend` to **2.3e-10
  rad** at the Aldrin point (≤8.3e-9 rad worst-case across the whole sweep).
  The asymptote-extraction error is therefore negligible vs the 1% Δv gate.
* **Asymptote baseline reconstruction vs `dv_from_turn_deficit`:** building the
  asymptote-rotation baseline from explicit V∞ vectors and a subtraction
  reproduces `dv_from_turn_deficit` to **machine precision** (≤2.2e-16). This
  validates that the harness's notion of "rotate the asymptote at infinity" is
  the same physics the baseline formula encodes.

## Aldrin-point result (the headline)

Geometry: V∞ = 6.86 km/s, Earth `rp_min` = 6678.137 km, `mu_E` = 398600.4355,
required turn 93.01° (ballistic cone `max_bend` = 67.99° here → deficit 25.02°;
the note's 68.49°/24.52° uses a 6578 km / 200 km-altitude periapsis — same
geometry class, the agreement check is method-vs-method on identical inputs).

| quantity | value |
|---|---|
| **integrated Δv** | **1.9720337236 km/s** |
| **formula Δv** (`dv_powered_flyby_periapsis`) | **1.9720337215 km/s** |
| relative Δv error | **1.05e-9** |
| turn achieved (root-solve) | 93.0100° (target 93.01°) |
| asymptote-extraction residual | 2.3e-10 rad |
| asymptote baseline (`dv_from_turn_deficit`) | 2.9134 km/s |

The integrated Oberth cost (1.972) sits well below the asymptote baseline
(2.913), confirming the Aldrin point is firmly in the favorable regime. (The
per-flyby ≈1.94–1.97 here is the dominant term of the note's 1.9336 km/s
*schedule total*; the schedule total under the asymptote model is 2.9138 —
both consistent with the per-point pair above.)

## Sweep + crossover (slow gate)

Grid V∞ ∈ {2,4,6.86,9,12} km/s × deficit ∈ {5°,15°,24.52°,40°} at Earth
`rp_min`, 20 points. Every point agrees to ≤1% relative; **max observed
rel_dv_error = 3.4e-5** (one near-escape steep-turn corner, V∞=2/deficit=40°),
all others ≤2e-7. Max asymptote residual across the sweep = 8.3e-9 rad.

**Crossover CONFIRMED to exist** (the #151 note asserts Oberth loses its
advantage at high V∞). Integrated Oberth cost vs asymptote baseline:

| V∞ (km/s) | deficit | dv_integrated | asymptote | Oberth wins? |
|---|---|---|---|---|
| 2 | 24.52° | 0.3074 | 0.8494 | yes |
| 4 | 24.52° | 0.8274 | 1.6988 | yes |
| **6.86 (Aldrin)** | **24.52°** | **1.9435** | **2.9134** | **yes** |
| 9 | 5° | 0.8514 | 0.7851 | **no (crossover)** |
| 12 | 5° | 1.5937 | 1.0469 | **no (crossover)** |
| 12 | 15° | 3.9340 | 3.1326 | **no (crossover)** |
| 12 | 24.52° | 5.5019 | 5.0964 | **no (crossover)** |
| 12 | 40° | 7.2564 | 8.2085 | yes |

The crossover is real and lands near V∞ ≈ 9–12 as the #151 note predicted: for
narrow deficits Oberth already loses by V∞=9; for the 40° deficit even V∞=12 is
still favorable. The Aldrin point (V∞=6.86, deficit ~24.5°) is unambiguously on
the favorable side.

## Test results (verbatim)

Fast layer (`-n0`, default selection):

```
collected 37 items / 21 deselected / 16 selected
tests/verify/test_flyby_integrate.py ................                    [100%]
====================== 16 passed, 21 deselected in 0.76s =======================
```

Full suite incl. slow sweep (`-m "slow or not slow" -n0`):

```
collected 37 items
tests/verify/test_flyby_integrate.py ................................... [ 94%]
..                                                                       [100%]
============================== 37 passed in 6.11s ==============================
```

## Provenance / publication recommendation

With this confirmation the **1.9336 km/s Aldrin maintenance schedule total is
now dual-method** — the analytic periapsis formula AND an independent direct
integration agree on the per-flyby cost to 1e-9 — but it remains **single
OPTIMISER** (one in-family DE440 solve produced the schedule; the integration
confirms the *cost model*, not the orbit-fit / family selection).

Recommendation to the main session (final decision deferred there):

* If publishing, present the maintenance ΔV as a **bound pair with model
  provenance**: **2.914 km/s (conservative, asymptote-rotation)** and **1.934
  km/s (Oberth-optimal, periapsis maneuver)**.
* State the cost model is dual-method-confirmed (formula + integration) but
  single-optimiser; an independent *orbit* solver carrying a periapsis cost
  model would still strengthen the family-selection side.
* Keep the asymptote value as the documented upper bound; the Oberth value is
  the physically better-motivated lower bound (the maneuver is paid deep in the
  well, as a real mission would).

## Files (task #154 — only these three touched)

* `src/cyclerfinder/verify/flyby_integrate.py` — integration + root-solve check.
* `tests/verify/test_flyby_integrate.py` — self-validation, Aldrin, sweep,
  crossover gates.
* `docs/notes/2026-06-07-oberth-independent-confirmation.md` — this note.
