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
| `quasi_cycler` | closes-up-to-rotation | YES (10–15 yr window) | 3–15 | cyclers-of-opportunity inside a planetary-alignment window |
| `precursor_mga` | non-repeating | YES (launch window) | 1 (insertion) | one-shot MGA chain inserting a spacecraft into a steady-state cycler |
| `mga_tour` | non-repeating | YES (launch window) | 1 (terminal) | Galileo VEEGA; Cassini VVEJGA; Tito 2018 |

The `cycler` class is the original scope and remains the gold standard. The
other three classes are admitted because they are *mission-actionable*,
structurally *searchable* (closure equations, not arbitrary flyby chains), and
have **strong existence priors in literature already in the corpus**.

## Schema additions (v4.6 → v4.7)

Three additive optional row fields (defaults preserve v4.6 behaviour — every
pre-v4.7 row reads as a strict cycler):

- `orbit_class`: enum {`cycler` | `quasi_cycler` | `precursor_mga` | `mga_tour`}, default `cycler`
- `epoch_locked`: bool, default `false`
- `n_returns`: integer ≥ 1 or string `"infinite"`, default `"infinite"`

Three optional epoch-locked fields (null/absent for `cycler`):

- `validity_window`: `{start, end}` ISO-8601 (UTC) — when the trajectory closes
- `launch_epoch`: ISO-8601 (UTC) — for `mga_tour` / `precursor_mga`
- `inserts_into`: catalogue ID — required for `precursor_mga`; must resolve to an extant `cycler`

Class invariants (enforced by the Python semantic gate, not the JSON Schema —
the schema cannot express cross-row referential integrity):

- `cycler` ⇒ `epoch_locked=false`, `n_returns="infinite"`
- `quasi_cycler` / `precursor_mga` / `mga_tour` ⇒ `epoch_locked=true`, `n_returns` is a finite integer ≥ 1
- `precursor_mga` ⇒ `inserts_into` resolves to an existing `cycler` row

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
