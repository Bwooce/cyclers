# #554: Neptune/Amalthea empty-region retrograde-correction stamp

**Date:** 2026-07-10
**Task:** #554 (formalizes the #552 scoping pass's back-of-envelope flyby-bend check)
**Subject:** `data/empty_regions.jsonl` rows `jupiter-galilean-amalthea-repeated-moon-quasicycler-2026-06-24`,
`repeated-moon-neptune-sweep`, `uranus-neptune-regular-moon-endgame-vilm-2026-06-23`,
`neptune-proteus-triton-proteus-multirev-leveraging-2026-06-26`.

## What this checks

#552 scoped a proposed 3D/inclined-releg moontour genome extension and found the payoff claim
("this might reopen the Amalthea and Neptune/Triton empty-region verdicts") does not survive a
flyby-bend physics check. This note reruns that check from the codebase's own sourced constants
and its own hyperbolic-flyby formula, rather than re-asserting #552's numbers, and records the
result as `reverification` entries on the four affected rows.

**Formula** (`src/cyclerfinder/core/flyby.py::max_bend`, unchanged, already the codebase's
physical-max-bend gate used by `search/correct.py`, `search/physical_sanity.py`,
`search/moon_prune.py`, etc.):

```
bend_max = 2 * arcsin( 1 / (1 + r_p * v_inf^2 / mu) )
```

where `mu` is the flyby body's GM (km^3/s^2), `r_p` is the periapsis radius from the body's
center (km) — `radius_eq_km + safe_alt_km` in this codebase's convention (see e.g.
`search/appc_corrected.py:304`) — and `v_inf` is the hyperbolic excess speed (km/s). This formula
has **no inclination, node, or geometry term at all**: it is a pure function of mass, minimum
approach altitude, and approach speed. Whatever a 3D/inclined-releg genome extension changes
(the *geometry* by which a spacecraft reaches a given `v_inf` at a given body), it cannot change
this ceiling for a fixed `(mu, r_p, v_inf)`.

## Amalthea: mass-limited, not geometry-limited

Sourced constants (`src/cyclerfinder/core/satellites.py:205`, JPL SSD `phys_par` ref JUP365):
`mu = 0.16456 km^3/s^2`, `radius_eq_km = 83.5`, `safe_alt_km = 10.0` (the codebase's engineering
floor for tiny irregulars, same convention used for Amalthea/Hyperion/Nix/Hydra) -> `r_p = 93.5 km`.

At the generous end (`v_inf = 0.5 km/s`, near-co-orbital, an optimistic low-V_inf case):

```
bend_max = 2*arcsin(1/(1 + 93.5*0.5^2/0.16456)) = 0.8011 deg
```

confirming #552's 0.80 deg figure to 4 significant figures. The bend collapses further, not
better, at any realistic *higher* V_inf — and the codebase's own #433 sweep
(`data/scan_433_jupiter_galilean.jsonl`) shows the actual V_inf realized at Amalthea encounters
in real closures is far higher than the generous 0.5 km/s case: e.g. idx 105/111/251 report
Amalthea-encounter `v_inf` of 5.15, 22.31, 37.24 km/s (Amalthea sits deep in Jupiter's well —
its own circular speed about Jupiter is ~26 km/s — so any transfer connecting it to an outer
moon inherits a large relative speed, not a small one):

| v_inf (km/s) | source | bend_max (deg) |
|---|---|---|
| 0.5 | generous back-of-envelope | 0.8011 |
| 5.1 | scan_433 idx 111, Io-Amalthea-Io | 0.007753 |
| 22.3 | scan_433 idx 105, Io-Amalthea-Io | 0.000406 |
| 37.2 | scan_433 idx 251, Europa-Amalthea-Europa | 0.000146 |

The real search data are *more* pessimistic than the back-of-envelope case, not less. Since the
formula carries no geometry term, no 3D/inclined-releg capability can rescue this — the ceiling is
set by `mu` and `r_p` alone, both already at their sourced/engineering-floor values.

**Geometry doesn't even differ much here anyway.** Amalthea's actual inclination to Jupiter's
equator is ~0.37 deg (JPL SSD mean orbital elements; not currently a coded field in
`satellites.py`, which has no inclination/node column at all). A 3D/inclined-releg extension would
have essentially nothing to work with for this body regardless of the bend-ceiling argument above.

**Conclusion:** Amalthea's registered empty verdict (`jupiter-galilean-amalthea-repeated-moon-
quasicycler-2026-06-24`) already correctly identifies the low-thrust releg (which can pay a
turn-angle deficit as `Delta V`, cf. `core/flyby.py::dv_from_turn_deficit`) as its live re-sweep
condition; that stays. The 3D/inclined-releg condition never held — it is removed.

## Triton: the code's coplanar-prograde model is generous, not conservative

Sourced constants (`src/cyclerfinder/core/satellites.py:187`): `mu = 1428.49546 km^3/s^2`,
`radius_eq_km = 1352.6`, `safe_alt_km = 100.0` -> `r_p = 1452.6 km`. Neptune system
`mu = 6.836527100580e6 km^3/s^2` (`satellites.py:62`, `PRIMARIES["Neptune"]`).

**The code models Triton as prograde.** `search/discovery_campaign.py::_moon_state` (lines
242-250) computes every moon's planet-frame position/velocity from a *positive* mean motion and
`z=0`:

```python
def _moon_state(theta0, n_rad_day, t_days, sma_km, mu):
    theta = theta0 + n_rad_day * t_days
    v_circ = math.sqrt(mu / sma_km)
    pos = np.array([sma_km*cos(theta), sma_km*sin(theta), 0.0])
    vel = np.array([-v_circ*sin(theta), v_circ*cos(theta), 0.0])
    return pos, vel
```

`n_rad_day` is derived from `SatelliteData.mean_motion_deg_day`, which is *always positive*
(`mean_motion_deg_day_about`, Kepler III, `satellites.py:32-35` — no sign/direction input exists).
There is no retrograde flag anywhere in `core/satellites.py` or `search/discovery_campaign.py`;
Triton's real retrograde motion is recorded only in comments (`satellites.py:113`, `181-182`;
`genome/bcr4bp_systems.py:348,355` carries the identical caveat for the BCR4BP lane). This
confirms #552's finding in code, not just doc claims: every moontour sweep that has ever touched
Triton — including the three Neptune "structurally empty" rows below — used a coplanar model
whose Triton is moving in the *wrong rotational direction*.

**Why the wrong direction is generous, not conservative.** Triton's own circular speed about
Neptune:

```
v_moon = sqrt(mu_Neptune / a_Triton) = sqrt(6.836527100580e6 / 354800) = 4.3896 km/s
```

For a spacecraft arriving at Triton's orbit with local speed `v_sc` (patched-conic, planet frame),
the encounter `v_inf` is the vector difference `|v_sc - v_moon|`. Model Triton's assumed
(WRONG, prograde) velocity and its TRUE (retrograde) velocity as the same vector reversed:
`v_moon_true = -v_moon_prograde`. So:

```
v_inf(model)  = |v_sc - v_moon|
v_inf(true)   = |v_sc - (-v_moon)| = |v_sc + v_moon|
```

For the near-co-orbital case that gives the model's most generous (lowest) `v_inf` — `v_sc` and
`v_moon` close in magnitude (~4.39 km/s) and separated by a small angle `phi` — the chord identity
for two equal-magnitude vectors gives:

```
v_inf(model) = 2*v_moon*sin(phi/2)      (small when phi is small)
v_inf(true)  = 2*v_moon*cos(phi/2)      (near 2*v_moon when phi is small)
```

Solving for `phi` at the model's generous `v_inf` range (1.0-2.0 km/s, the range #552 flagged as
giving 23-59 deg bend — reproduced below) and evaluating the corrected `v_inf(true)`:

| v_inf(model), prograde (km/s) | phi (deg) | v_inf(true), retrograde (km/s) | bend_max(model) (deg) | bend_max(true) (deg) |
|---|---|---|---|---|
| 1.0 | 13.08 | 8.722 | 59.447 | 1.462 |
| 1.5 | 19.68 | 8.650 | 35.413 | 1.487 |
| 2.0 | 26.34 | 8.548 | 22.762 | 1.522 |

(`bend_max` via the same `max_bend` formula with Triton's `mu = 1428.49546`, `r_p = 1452.6 km`.)
This reproduces #552's claimed ranges almost exactly: model `v_inf` 1-2 km/s -> bend 23-59 deg;
corrected `v_inf` converges to ~8.5-8.7 km/s (within the claimed "6-9+ km/s") -> bend collapses to
~1.46-1.52 deg (within the claimed "1-2 deg"). The correction is *self-consistent regardless of
`phi`* over this whole generous-model range: correcting the model's rotation direction roughly
doubles the encounter speed and collapses the achievable bend by more than an order of magnitude.

**Conclusion:** the three Neptune/Triton-adjacent "structurally empty" rows
(`repeated-moon-neptune-sweep`, `uranus-neptune-regular-moon-endgame-vilm-2026-06-23`,
`neptune-proteus-triton-proteus-multirev-leveraging-2026-06-26`) were all produced by a Triton
model that is, if anything, MORE favorable to a Triton-involving closure than reality. Correcting
the rotation direction can only raise the achievable encounter `v_inf` and shrink the achievable
bend — it strengthens these negatives, it cannot reopen them. A genuinely correct 3D/retrograde
Triton model is not a re-sweep condition worth carrying on these rows.

## Registry action

`reverification` entries appended to all four rows in `data/empty_regions.jsonl` recording this
check (see the entries themselves for exact wording; `src/cyclerfinder/data/empty_regions.py`
tolerates the extra key — `_from_payload` only reads the fields on `EmptyRegionReport`, the same
pattern already used by the `repeated-moon-neptune-sweep` row's pre-existing `#425` staleness
`reverification` entry). Net effect:

- Amalthea row: live re-sweep condition narrows to **low-thrust relegs only** (the "3D/inclined
  relegs" clause is retracted as false).
- The 3 Neptune/Triton rows: verdicts marked **strengthened** under a correctly-modeled retrograde
  Triton (no re-sweep condition added; a future correctly-modeled retrograde/inclined Triton lane
  would need to justify itself on some OTHER basis, not on reopening these).

No `data/catalogue.yaml` edit, no new search, no new genome capability. Pure registry-accuracy
work per #554's scope.
