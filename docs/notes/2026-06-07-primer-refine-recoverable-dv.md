# Add-an-impulse refinement — recoverable maintenance ΔV (2026-06-07)

**Status: DIAGNOSTIC / PROVISIONAL. NOT for the site or the catalogue.** Task
#148. Quantifies the recoverable ΔV implied by the #144 primer-vector verdict,
via the Lion & Handelsman (1968) add-an-impulse step. Implemented in
`src/cyclerfinder/verify/primer_refine.py`; tests in
`tests/verify/test_primer_refine.py`. Reads on from the #144 diagnostic note
`docs/notes/2026-06-07-primer-vector-diagnostic.md`.

> **Addendum (2026-06-07, per Guzman 2002 mining):** the "DIAGNOSTIC /
> PROVISIONAL, multi-rev caveats apply" hedge is **upgraded to CONFIRMED per
> Guzman 2002 pp.5-6**. The now-acquired survey corroborates the zero-recovery
> result from its own equations (Eq. 43 δJ ∝ `|Δv₀| = 0`, and the Eq. 40/41
> first-order-gradient-vs-realised-magnitude separation), and the multi-rev
> caveat does not bind (sub-one-rev arc). See
> `docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`. The
> "NOT for the site/catalogue" honesty framing is unchanged.

## What it computes

A single coast is a two-impulse, fixed-endpoint transfer from boundary state A
(`r_A` at `t_A`, imposed pre-departure velocity `v_A_before`) to boundary state
B (`r_B` at `t_B`, imposed post-arrival velocity `v_B_after`). The original cost
is the departure + arrival ΔV of the single Lambert arc A→B:

    ΔV_orig = |v_dep − v_A_before| + |v_B_after − v_arr|

The add-an-impulse step (Lion & Handelsman 1968, AIAA J. 6(1):127-132, DOI
10.2514/3.4452; textbook treatment Prussing & Conway Ch. 2) inserts a midcourse
impulse at an interior time `t_m` and position `r_m`, splitting the coast into
Lambert arcs A→m and m→B:

    ΔV_refined = |v1_Am − v_A_before| + |v1_mB − v2_Am| + |v_B_after − v2_mB|

The **boundary states (endpoints, encounter times, and the imposed boundary
velocities) are held fixed**, so original vs refined is like-for-like. The free
parameters are `(τ, δr)`: the fractional midcourse time `τ = (t_m−t_A)/(t_B−t_A)`
and the 3-D offset `δr` of `r_m` from the original-arc position at that time.
`δr = 0` reproduces the original arc to machine precision (the midcourse burn
vanishes), so the optimum is never worse than the original. The step is **seeded
from the primer peak** — `τ* = t*/T` and `δr ∝ p̂(t*)` from the #144 primer BVP
— and scipy (Nelder-Mead; the |·| burns give a kinked objective) minimises the
4-vector. **One add-an-impulse step only; not iterated to exhaustion.**

## Mechanics gate (constructed, sourced-discipline)

The same fixed-time "long-way" (transfer angle > 180°) Lambert transfer that
#144 flags IMPROVABLE (interior `|p| > 1`) must refine strictly downward after
one step — qualitative teeth, no invented magnitude asserted. Canonical units
(μ = 1), boundary velocities = local circular velocities at each end:

| case (r1→r2, angle, ToF×Hohmann) | ΔV_orig | ΔV_refined | recovered |
|---|---|---|---|
| 1→4, 200°, 1.0× | 0.49290 | 0.44933 | 0.04357 (8.84%) |
| 1→4, 250°, 1.5× | 0.75609 | 0.46522 | 0.29087 (38.47%) |

The optimiser demonstrably recovers ΔV when a genuine heliocentric interior
bulge exists. The tests assert only `refined < original` and `recovered > 0`.

## Application — Aldrin E-M-E maintenance schedule, coast 0 (E→M)

Real ephemeris (`astropy`, real-window priority 1985-10-28), same schedule as
#144 (`optimise_aldrin_maintenance_dv`, total maintenance ΔV 2.9138 km/s).
Coast 0 is the 131.9-day Earth→Mars leg with #144 grid-converged `max|p| =
1.1223` at `t/T = 0.355`.

**Result (gated `@pytest.mark.slow`, ~19 s):**

```
schedule maintenance dv (km/s, our value) = 2.9138
original heliocentric coast-0 dv (km/s)   = 0.003238
refined  heliocentric coast-0 dv (km/s)   = 0.003238
realised recoverable (our optimiser)      = 0.000000 km/s (0.00% of coast-0 cost)
seed primer-peak t/T = 0.355; optimised midcourse t/T = 0.355
```

The seed primer-peak `t/T = 0.355` reproduces #144 exactly (cross-check on the
primer machinery).

### Original vs refined ΔV → recoverable ΔV (the headline numbers)

- **Original heliocentric coast-0 two-impulse ΔV: 0.003238 km/s.**
- **Refined (one add-an-impulse step): 0.003238 km/s.**
- **Realised recoverable ΔV (our optimiser): 0.000000 km/s (0.00%).**

### Honest interpretation (the crux)

The **heliocentric** burns on coast 0 are essentially nil. The Earth departure
burn `|v_depart − (v_planet + v∞_in)|` is **exactly 0.0** (the leg departs at
the encounter velocity by construction), and the Mars-end heliocentric
discontinuity is only **0.003238 km/s**. So the heliocentric coast is already
near-ballistic, and an interior heliocentric impulse has almost nothing to
recover — the realised recoverable ΔV is 0.

This is the load-bearing finding: **the schedule's 2.9138 km/s does not live on
the heliocentric coast impulses at all.** It is the Earth **flyby turn-deficit**
(`idealized_flyby_turn_deficit` in `search/maintain.py`) — the geocentric turn
the Earth flyby cannot deliver ballistically (≈84° required vs ≈72° achievable),
a purely *geometric* cost charged at the flyby, not a heliocentric midcourse
burn. **An interior heliocentric add-an-impulse cannot recover a flyby turn
deficit** — different actuator, different place in the trajectory.

The #144 primer IMPROVABLE flag (`max|p| = 1.122`) is a genuine necessary-
condition violation **of the heliocentric arc's optimality**, but the
heliocentric ΔV it implies is recoverable is ~0 here, because the heliocentric
burns are already ~0. The primer bulge reflects the endpoint primer directions
(set from the heliocentric ΔV directions, one of which is the degenerate
zero-burn Earth departure → fallback direction), not a large recoverable
heliocentric impulse. First-order primer theory predicted *some* heliocentric
improvement exists; the realised magnitude from our optimiser is negligible.

**Bottom line (DIAGNOSTIC, PROVISIONAL):** the recoverable maintenance ΔV on
the Aldrin schedule via heliocentric add-an-impulse is **~0 km/s** — the
maintenance cost is flyby-geometric and is not addressable by this step. To
reduce the 2.9138 km/s one must change the *flyby geometry* (encounter epochs /
leg ToFs / the flyby altitude), which is the maintenance optimiser's job, not a
heliocentric primer refinement. The mechanics gate confirms the refinement
itself is correct and has teeth; it is simply the wrong lever for this schedule.

## Honesty framing (binding)

- The recoverable-ΔV numbers here are **OUR computation** — realised improvement
  from our optimiser, **not primer-predicted**. First-order primer theory only
  asserts that an improvement *exists* on a coast with `|p| > 1`; the magnitude
  comes from the optimiser.
- **CONFIRMED per Guzman 2002 pp.5-6.** The zero-recovery finding is the
  *expected* outcome of a degenerate-endpoint (`|Δv₀| = 0`) coast: the survey
  separates the Eq. 40/41 first-order gradient (improvable) from the realised
  impulse magnitude (which "might not produce an optimal trajectory in the sense
  of Lawden"), and its endpoint-time formula Eq. 43 (δJ ∝ `|Δv₀|`) gives zero
  recoverable when `|Δv₀| = 0` exactly — both predict #148's ~0. The multi-rev
  caveat does **not** bind (sub-one-rev arc). See
  `docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`.
- These numbers must **NOT** be published to the site or written to the
  catalogue. They live in this note only.

## What is NOT done (deliberate scope limits)

- One add-an-impulse step only; not iterated, and no second-bulge re-check
  needed (coast-0 recoverable is 0).
- Coast 1 (M→E, marginal `max|p| = 1.00008` at the endpoint, #144) is not
  refined: its violation is at the noise floor and at an endpoint, not an
  interior bulge — nothing to add an interior impulse for.
- No change to the flyby-geometric maintenance cost (out of scope; that is the
  maintenance optimiser, `search/maintain.py`, which is read-only here).
