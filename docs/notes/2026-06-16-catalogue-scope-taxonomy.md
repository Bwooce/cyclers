# Catalogue scope taxonomy (schema v4.7, 2026-06-15)

**Decision (2026-06-15):** the catalogue's scope expanded from cyclers-only to a
four-class taxonomy. Driven by the #286 frontier-scoping finding that the
literature has a mature class of *epoch-locked* trajectories — Galileo VEEGA,
Cassini VVEJGA, Petropoulos pump tours, Tito 2018 Mars free-return — that the
prior scope could not represent. A catalogue restricted to strict cyclers
describes orbits no spacecraft can actually reach: a cycler is reached *via* an
insertion trajectory, and that insertion is epoch-locked.

## The four classes

| Class | Period? | Epoch-locked? | Returns | Type case |
|---|---|---|---|---|
| `cycler` | strictly periodic | NO | ∞ | Aldrin Earth-Mars; Russell-Ocampo S1L1; Braik-Ross C11a/C11b/C32 |
| `quasi_cycler` | closes-up-to-rotation | YES (10–15 system-periods) | 3–15 | cyclers-of-opportunity inside a planetary-alignment window |
| `precursor_mga` | non-repeating | YES (launch window) | 1 (insertion) | one-shot MGA chain inserting a spacecraft into a steady-state cycler |
| `mga_tour` | non-repeating | YES (launch window) | 1 (terminal) | Galileo VEEGA; Cassini VVEJGA; Tito 2018 |
| `resonant_po` (v4.9, #453) | strictly periodic | NO | ∞ | stable resonant/libration PO, NO transport utility (em-cycler-21-3d-spatial-2026) |

The `cycler` class is the original scope and remains the gold standard. The
fifth class, `resonant_po` (added in schema v4.9, task #453), is NOT
mission-actionable: it is a stable resonant/libration periodic orbit that
shares the `cycler` reachability/periodicity invariants but never encounters
the secondary (no transport). It exists so the catalogue can carry a
known-class-member resonant PO without mislabelling it a transport cycler — see
the C21 3D spatial extension, recharacterized from `cycler` after Region B
(#447) measured its closest lunar pass at 122,628 km vs the 66,183 km lunar SOI.
The
other three classes are admitted because they are *mission-actionable*,
structurally *searchable* (closure equations, not arbitrary flyby chains), and
have **strong existence priors in literature already in the corpus**.

### `quasi_cycler` window is system-relative, not a fixed 10-15 years (corrected 2026-07-11, #557)

The "10-15 yr" figure was an Earth-calibrated example, never actually enforced
anywhere (`catalogue.schema.json`'s `validity_window` is just an ISO-8601
`{start, end}` pair with no duration constraint, and
`data/validate.py`'s semantic gate only checks `epoch_locked=true` + a finite
`n_returns >= 1` — no window-length check exists in code). It was also already
silently wrong: the catalogue's one existing `quasi_cycler` row (task #312,
Umbriel-Oberon-Umbriel Uranian cycler) has an 83-year `validity_window`
(2000-06-21 to 2083-06-21), 5.5x the stated ceiling, admitted through the
full V1-V4-strict gauntlet with every ratchet green.

**Corrected reading:** the return-separation floor and the window length are
both stated in units of **one rotating-frame period of the CR3BP system under
study** — floor = 1 system-period, window = 10-15 system-periods. This is
what #535's Earth calibration (1 yr floor, 10-15 yr window) already was,
just written in absolute years because Earth's system-period happens to be
1 year. For a Sun-planet quasi_cycler this is the planet's heliocentric
period (e.g. Sun-Jupiter: ~11.86 yr floor, ~119-178 yr window); for a
**moon-system** quasi_cycler like #312 itself, the relevant system-period is
the MOON-PAIR's synodic period (~6 days for Umbriel-Oberon), NOT the
primary's heliocentric period (~84 years for Uranus) — using the wrong
period for a moon-pair case would be wrong by a factor of ~5000. `n_returns`
(3-15) and the bounded-geometry ratio (3.0) are already dimensionless and
unchanged.

Two real consequences of this correction (not just a units relabeling):

1. **`validity_window` has two coexisting semantics that any new criterion
   note must name explicitly.** #312's row uses it as a LAUNCH-EPOCH span
   (any June-21 within the 83-year range is a valid start; the 10 returns
   themselves accumulate over ~150 days at 14.94 d/leg). #535's Earth
   criterion used it as an ENCOUNTER-ACCUMULATION span (all 3-15 returns
   must fall inside the window itself). Both are legitimate but distinct;
   a criterion note must say which one it means.
2. **V2 (long-span bounded drift) and V4 (real-ephemeris) validation budgets
   scale in real calendar years and can exceed kernel coverage** at
   giant-planet-heliocentric scale: DE440 spans ~1550-2650 (~1100 yr) —
   Jupiter-scale windows (119-178 yr) fit easily, Saturn-scale (294-442 yr)
   fits, but Uranus-heliocentric-scale (840-1260 yr) is marginal and
   Neptune-heliocentric-scale (1650-2472 yr) EXCEEDS DE440 entirely; moon
   satellite kernels are far shorter still (URA111: 1900-2099 only — #312's
   own V4 failures already cluster in the last 15 years of that coverage,
   `data/validate.py` line ~1152). Any long-window V2/V4 validation attempt
   must evaluate over `min(validity_window, kernel_coverage)` with the
   truncation explicitly recorded, not silently assume full-window coverage
   is achievable.

**Explicitly NOT changed by this correction:** no re-validation of past
negatives is triggered. #535's Earth-system numbers already ARE the
dimensionless criterion evaluated at a 1-year system-period; nothing about
them was wrong. Any future giant-planet-scale quasi_cycler screen (e.g.
Sun-Jupiter quasi-Hilda, task #557) that produces a negative under its own
system-period-scaled criterion must be registered in the empty-region
registry under a criterion-version tag distinct from any Earth-scale
negative, so a future reader doesn't compare across incompatible window
definitions.

## Schema additions (v4.6 → v4.7)

Three additive optional row fields (defaults preserve v4.6 behaviour — every
pre-v4.7 row reads as a strict cycler):

- `orbit_class`: enum {`cycler` | `quasi_cycler` | `precursor_mga` | `mga_tour` | `resonant_po`}, default `cycler` (`resonant_po` added v4.9, #453)
- `epoch_locked`: bool, default `false`
- `n_returns`: integer ≥ 1 or string `"infinite"`, default `"infinite"`

Three optional epoch-locked fields (null/absent for `cycler`):

- `validity_window`: `{start, end}` ISO-8601 (UTC) — when the trajectory closes
- `launch_epoch`: ISO-8601 (UTC) — for `mga_tour` / `precursor_mga`
- `inserts_into`: catalogue ID — required for `precursor_mga`; must resolve to an extant `cycler`

Class invariants (enforced by the Python semantic gate, not the JSON Schema —
the schema cannot express cross-row referential integrity):

- `cycler` / `resonant_po` ⇒ `epoch_locked=false`, `n_returns="infinite"` (pinned by `tests/data/test_schema_v47_orbit_class.py`)
- `quasi_cycler` / `precursor_mga` / `mga_tour` ⇒ `epoch_locked=true`, `n_returns` is a finite integer ≥ 1
- `precursor_mga` ⇒ `inserts_into` resolves to an existing `cycler` row

## `precursor_mga` sub-classification (additive convention, task #386)

The `precursor_mga` class covers any non-repeating chain that inserts a
spacecraft *into* a steady-state cycler. Per Pontani & Conway 2018 ("Optimal
Trajectories for Hyperbolic Rendezvous with Earth-Mars Cycling Spacecraft",
*J. Guid. Control Dyn.*, DOI 10.2514/1.G002984), that insertion can be
mechanistically distinct, so the class carries an **additive
`sub_classification` annotation** (task #386, 2026-06-17 digest wave; see
`docs/notes/2026-06-17-digest-pontani-2018.md` Sec 6.2). This is a free-text
controlled vocabulary, **not a fifth `orbit_class` enum** — Pontani 2018 is a
methodology paper and is *not* admitted as a row.

| `sub_classification` | Insertion mechanism | Type case |
|---|---|---|
| `mga_sequence` | multi-flyby V-infinity-leveraging gravity-assist chain | Rogers 2015 establishment trajectories (the 20 `*-establishment*` rows) |
| `rendezvous_only` | pure hyperbolic rendezvous, no gravity assists | Pontani 2018 impulsive taxi (methodology only; not catalogued) |
| `low_thrust_rendezvous` | continuous-thrust (electric) rendezvous taxi | Pontani 2018 low-thrust taxi (methodology only; not catalogued) |

The field is optional and additive (the schema's `additionalProperties: true`
permits it); it changes no class/tier/validation census. The 20 Rogers 2015
establishment rows were tagged `mga_sequence` under #386. If a second/third
taxi-side paper is ever catalogued, schema v4.x may promote the vocabulary to a
JSON-Schema enum (Pontani 2018 Sec 6.2 Option A).

## Migration (one-time, 2026-06-15)

`scripts/migrate_catalogue_scope_2026-06-15.py` annotates every pre-v4.7 row
with `orbit_class: cycler`, `epoch_locked: false`, `n_returns: "infinite"`,
preserving comments and formatting via ruamel.yaml. **No row counts changed**;
the frozen census ratchet
(`tests/test_catalogue_rediscovery.py::EXPECTED_COVERAGE`) was unchanged by the
migration itself.

The Tito 2018 Mars free-return was admitted as the first `mga_tour` row in a
follow-up commit, bumping `MULTI_ENCOUNTER_SEQUENCE` 223 → 224 (the v1 gauntlet
has no `mga_tour` lane yet, so the classifier files E-M-E under the
multi-encounter bucket).

## V0-V5 gauntlet extension (per-class semantics)

The gauntlet keeps the same V0-V5 ladder shape but the closure definition
becomes class-dependent:

- **V0** sourced + reproduced — unchanged.
- **V1** same-model <1 m/s residual — for epoch-locked classes, evaluated *within* `validity_window`.
- **V2** long-span bounded-drift — for `quasi_cycler`, "long-span" means the full `validity_window`, not infinite; for `mga_tour` / `precursor_mga`, V2 collapses to "from `launch_epoch` to terminal-epoch only" — no drift requirement past the window.
- **V3** independent-model — unchanged, evaluated within window.
- **V4** HFEM real-ephemeris (DE440) — unchanged; epoch-locked classes shine here because they were designed against real ephemeris in the first place.
- **V5** mission-quality — unchanged.

A `precursor_mga` row that does not `inserts_into` an extant `cycler` row fails
V0 by definition.

## Literature-novelty corpus expansion

`src/cyclerfinder/search/literature_check.py::KNOWN_CORPUS` gained 12 anchors
for the new classes (commit `568d8a4`):
Petropoulos-Longuski 2000, Strange-Russell 2007, Heaton-Strange-Longuski 2002,
Vasile-Conway 2006, Hughes-Edelman-Longuski 2014, Genova-Aldrin 2015,
McConaghy 2004, Ceriotti 2010, Vasile-Campagnola 2009,
Diehl-Belbruno-Roberts 1986 (Galileo VEEGA), Tito-MacCallum 2018, plus
Antoniadou-Voyatzis 2018 (the spatial-CR3BP corpus the #287 3D-Aldrin spike
likely rediscovered). A `mga_tour` candidate is not "novel" until it clears
this corpus, same as cyclers clear the cycler corpus.

## Website filter UI (`cyclers.space`, schema v5 site-side)

The site's catalogue filter bar (commits `c2bb4ee` + `c620492` on
`Bwooce/cyclers.space:main`) gained three filters and a class-colored badge:

- **Class** dropdown: All / Cyclers / Quasi-cyclers / Precursors / Tours
- **Epoch window** dropdown: All / Open now / Future / Past (auto-greyed when Class=Cyclers)
- **n_returns** numeric range (auto-greyed when Class=Cyclers)
- Per-row class badge (green / blue / amber / purple)
- Validity-window column rendered conditionally (only when ≥1 epoch-locked row exists)
- `inserts_into` rendered as a clickable cycler link on `precursor_mga` rows

The site fetches `data/catalogue.yaml` from the main repo at predev/prebuild
(`scripts/sync-catalogue.mjs`), so the schema v4.7 fields reach the site
automatically on the next build.

## Tito 2018 disposition (revised)

Was withheld permanently under the prior cyclers-only scope (see
`project_catalogue_scope_cyclers_only` memory, now superseded). Under v4.7 it
qualifies as `mga_tour` with `n_returns: 1`, `launch_epoch: 2018-01-05T07:00:00Z`,
`validity_window: 2018-01-05 → 2019-05-21`. The DE440 reproduction (<1.5% to
Tito's published DE421 Tables III/IV; flyby ballistic continuity 33.4° vs 34.2°
cone) carries directly as V0 evidence. Catalogue ID:
`tito-2018-mars-free-return`.

## Out-of-scope (still)

- **Strict one-shot trajectories with no published existence anchor** — admission
  still requires a sourced reference (literature, mission record, peer-reviewed
  preprint). The four-class taxonomy is about *which mature literature classes
  the catalogue admits*, not about lowering the source bar.
- **Tito-2018-class one-shots without a usable reproduction** — V0 still requires
  the row to reproduce on real ephemeris. The Tito row passes; a row whose
  reproduction fails would be a `validation_artifact:` only (the existing
  pattern), never a graded `mga_tour`.

## What this enables (and what it does not)

This expansion **unlocks** literature gating + catalogue admission for the
quasi-cycler / precursor / tour classes. It **does not** automatically build
the search-side genome that would *generate* discoveries in those classes — that
is task #289 (quasi-cycler + precursor MGA genome build, ~1200 LOC / 2-3 wk,
Track A).

## References

- Memory: `project_catalogue_scope_expanded_2026-06-15` (current);
  `project_catalogue_scope_cyclers_only` (superseded).
- Tito 2018 reproduction: `docs/notes/2026-06-13-tito-maccallum-2018-free-return-reproduction.md`.
- Frontier scoping that surfaced the gap: `docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md`.
- Discovery program spec: `docs/notes/2026-06-13-discovery-program-spec.md`.
