# 2026-06-16 — NEA-augmented cycler search Phase 1 (#308)

## Motivation: where #302 said the fresh ground is

Task #302's structural conclusion (commit `c3d433a`) closed the door on
ballistic Earth-Mars cycler insertion with V/E intermediates — it is fully
published, and exhaustive re-sweeping does not find new families. The
conclusion explicitly identified three fresh-ground directions:

1. **Asteroid leveraging** (NEAs as gravity-assist nodes, not just flyby
   anchors).
2. **Low-thrust** (powered cyclers, opening the bend-deficit gap).
3. **Non-Earth-Mars** (Jovian moon tours, Uranian moon tours, etc.).

This task targets (1) — the sparsest of the three by the
`literature_check.KNOWN_CORPUS` audit. None of the corpus's ~41 anchors lists
a NEA as a cycler node; the Petropoulos-Longuski pump-tour combinatorics
include some asteroid work but as flyby anchors, not as primary structural
elements.

## What Phase 1 built

Two new files:

* `src/cyclerfinder/search/asteroid_leveraging.py` — the
  `NEAEphemeris` value type, the built-in 10-NEA pool (`LARGEST_NEAS`,
  sourced JPL SBDB orbital elements + mission/radar-sourced masses + radii),
  and the `search_nea_augmented_cyclers` driver that enumerates length-3
  heliocentric chains (E -> NEA -> M and analogues) with single-rev prograde
  Lambert per leg + ballistic flyby continuity + the #324 physical-sanity
  max-bend gate at every encounter.
* `tests/search/test_asteroid_leveraging.py` — 27 tests locking the sourced
  data, the formula sanity, the search shape, the gate behaviour, and the
  override hooks.

The scan script `scripts/scan_308_em_nea_m.py` exercises the driver over the
default V_inf grid (3-7 km/s) x 6x6 TOFs (60-500 d) = 1800 cells per NEA x
10 NEAs.

### NEA pool — sourcing

Ten anchors spanning the three Earth-crossing dynamical groups
(Aten/Apollo/Amor). All orbital elements from JPL SBDB (epoch
JD 2461200.5 ≈ 2026-06-30). Mass sources cited per-record in the module:

| NEA        | Designation       | SMA (AU) | Mass (kg)   | R (km) | Source class |
|------------|-------------------|----------|-------------|--------|--------------|
| Eros       | 433 Eros          | 1.4582   | 6.687e15    | 8.42   | spacecraft   |
| Ganymed    | 1036 Ganymed      | 2.6622   | 7.75e16     | 19.00  | density-mod  |
| Apophis    | 99942 Apophis     | 0.9224   | 6.1e10      | 0.170  | radar+orbit  |
| Bennu      | 101955 Bennu      | 1.1264   | 7.329e10    | 0.245  | spacecraft   |
| Itokawa    | 25143 Itokawa     | 1.3243   | 3.51e10     | 0.165  | spacecraft   |
| Ryugu      | 162173 Ryugu      | 1.1896   | 4.50e11     | 0.448  | spacecraft   |
| Didymos    | 65803 Didymos     | 1.6442   | 5.28e11     | 0.390  | spacecraft   |
| Toutatis   | 4179 Toutatis     | 2.5418   | 5.05e13     | 1.34   | spacecraft   |
| Geographos | 1620 Geographos   | 1.2455   | 2.5e13      | 1.30   | density-mod  |
| Castalia   | 4769 Castalia     | 1.0633   | 3.6e12      | 0.70   | density-mod  |

Seven have spacecraft-tracked or radar-fit masses (the highest precision);
three (Ganymed, Geographos, Castalia) are density-modelled from radar
shape + assumed S-type density (lower precision; flagged per-entry).

### Ephemeris fidelity

Phase 1 uses **circular-coplanar** Keplerian states for the NEAs (state at
SMA, prograde, theta = 0 at t = 0 — matching the existing `_CircularBackend`
planet model). Eccentricity and inclination are recorded but NOT consumed by
`NEAEphemeris.state(t_sec)`; they are available for Phase 2+ to switch to an
inclined-eccentric Kepler propagator without changing the record schema.

This matches the fidelity rung of every existing CR3BP / planet
circular-Lambert search in the tree; Phase 2 (real DE440 + epoch-locked
NEAs) is a deferred extension.

## The honest Phase 1 verdict

A NEA flyby's max ballistic bend is

```
delta_max = 2 * arcsin( mu_NEA / (mu_NEA + r_p * V_inf^2) )
```

with `mu_NEA = G * M`. The NEAs in this pool have masses 6e10 - 8e16 kg, so
mu sits at 4e-9 - 5e-3 km^3/s^2 — six to ten orders of magnitude below the
smallest planet in `PLANETS` (Mercury at 2.2e4). Concretely, at the
minimum-safe periapsis radius `r_p = R_NEA + safe_alt_km` (5 km for sub-km
NEAs; 50 km for the larger ones):

| NEA      | v_escape (m/s) | bend @ v_esc | bend @ 1 km/s     | bend @ 5 km/s    |
|----------|----------------|--------------|-------------------|------------------|
| Eros     | 10.30          | 7.71 deg     | 0.00088 deg       | 3.5e-5 deg       |
| Ganymed  | 23.33          | 13.90 deg    | 0.0086 deg        | 3.4e-4 deg       |
| Bennu    | 0.20           | 2.62 deg     | 1.07e-7 deg       | ~ 4e-9 deg       |
| Apophis  | 0.22           | 1.85 deg     | 9.02e-8 deg       | ~ 4e-9 deg       |

At the operational cycler V_inf grid (3-7 km/s) **no NEA in the pool
clears the 5 deg useful-bend floor**. The largest NEA (Ganymed) only clears
the floor at V_inf comparable to its surface escape speed (~23 m/s),
which is six orders of magnitude below the cycler regime.

### Scan results (commit `e6ac690`, 1800 cells)

```
Enumerated cells:             1800
Lambert+closure passers:      45
   per-NEA: Eros 15, Ganymed/Bennu/Itokawa/Didymos/Toutatis/Castalia 5 each,
            Apophis/Ryugu/Geographos 0 (TOF box did not bracket their
            single-rev geometry)
Physical-sanity gate passers: 0
```

The zero gate-passing count IS the Phase 1 answer. A ballistic NEA flyby at
cycler V_inf is geometrically vacuous — the spacecraft sails past with
negligible deflection. The gate is doing its job when it rejects these.

The 45 Lambert-converging pre-gate candidates are kept in the JSONL as a
diagnostic — they map the **geometric** (Lambert + closure) region that Phase
2 should target, even though none of them have any usable ballistic bend.

## Reading the JSONL

`data/scan_308_em_nea_m.jsonl` has three row kinds:

* **Leading `_meta` config row** — sweep parameters, NEA pool, discipline
  note.
* **Per-candidate rows** — none in this run (`post_gate` count was 0); when
  Phase 2 wires powered slingshots, gate-passing rows will appear here first.
* **Pre-gate diagnostic rows** — capped at 200 (sampled), `pass:
  pre_gate_diagnostic_sample`. The 45 rows from this run map the geometric
  feasibility region per NEA.
* **Trailing `_meta` summary row** — per-NEA breakdown of pre-gate and
  post-gate counts, the Phase 1 verdict, the Phase 2 recommendation.

## Phase 2 path (NOT in scope for this task)

The bend-deficit gap that Phase 1's gate exposes is exactly the gap
**low-thrust** opens: a powered slingshot at the NEA encounter modifies the
heliocentric V_inf vector independent of the NEA's gravity, so the
chain's V_inf-shell continuity decouples from the NEA's vacuous bend.
Concretely, Phase 2 should:

1. Re-load the 45 pre-gate diagnostic rows from this run as initial guesses.
2. Wrap each with a Sims-Flanagan low-thrust segment in the NEA-encounter
   leg (the existing `cyclerfinder.search.low_thrust_cycler_search` machinery
   is the natural home; reused from #309).
3. Re-test the gate AFTER the low-thrust segment: if the spacecraft V_inf at
   the NEA encounter is brought low enough by the powered burn that the
   patched-conic bend clears the floor, the NEA flyby contributes geometric
   bending the burn can leverage.

A complementary Phase 2 path: pick the **handful of larger NEAs** (asteroidal
mass > 1e16 kg) and broaden the V_inf grid down to 0.01-0.5 km/s — within
those bounds the largest NEAs in the pool DO pass the gate, and the scan
becomes a useful-NEA-flyby V_inf map (not a cycler, but groundwork for
multi-NEA tours).

## Sourcing discipline

* All orbital elements from JPL SBDB (queried 2026-06-16).
* Mass sources are spacecraft-tracking papers (Konopliv et al. 2002 for
  Eros; Scheeres et al. 2019 for Bennu; Watanabe et al. 2019 for Ryugu;
  Brozovic et al. 2018 for Apophis; etc.) or peer-reviewed radar-shape +
  density papers (Hudson & Ostro 1994/1999 for Castalia/Geographos; Krasinsky
  et al. 2002 for Ganymed). Per-record citations are inline in the module.
* The 5 deg max-bend floor is the same #324 judgment threshold used for the
  Umbriel-Oberon-Umbriel rejection. NOT a sourced physical constant.

## What was NOT done

* **No catalogue writeback** — Phase 1 output is JSONL + this doc.
* **No literature check** — the `literature_check.py` baseline gate (per
  `feedback_literature_novelty_check_baseline`) is downstream and mandatory
  before any "novel" framing. Phase 1 produces 0 gate-passing rows so the
  question doesn't arise; if Phase 2 produces survivors, the lit-check runs
  before any anchor.
* **No real ephemeris** — circular-coplanar Keplerian NEAs; DE440 epoch-locked
  is Phase 2+.
* **No multi-NEA chains** — `max_nea_in_chain=1` hardcap; multi-NEA tours
  are Phase 2+.

## Discipline notes for future runs

* `src/cyclerfinder/search/physical_sanity.py` (#324) — read-only on this
  task. The gate is the right shape; it correctly rejects every NEA at
  cycler V_inf.
* `LARGEST_NEAS` is the discoverable surface — extending the NEA pool is a
  straightforward addition (the API takes any iterable of `NEAEphemeris`),
  but Phase 1 deliberately caps at the 10 largest measured-mass NEAs so the
  zero-survivor verdict is unambiguous.
