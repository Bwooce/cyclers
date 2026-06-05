# data/

Top-level project data files. Consumed by:
- M3 golden tests (when those land) — reads `catalogue.yaml`.
- The `cyclers.space` static site repo — pulls `catalogue.yaml` directly.
- M7's catalogue loader (`src/cyclerfinder/data/catalog.py`, future) — reads `catalogue.yaml`.

Files
-----

- `catalogue.yaml` — published-cycler seed catalogue with full attribution per spec.md §16. **Sole source of truth.** Edits go through the same process as code: change values, change source quotes, commit.
- `OUTSTANDING.md` — long-form research questions / source-access gaps / paradigm-mismatch flags log. The YAML's per-entry `notes:` field handles short-form caveats; OUTSTANDING handles the discussion threads.

Conventions
-----------

Three rules — non-negotiable:

1. **Single source of truth.** The YAML at `catalogue.yaml` is the canonical data record. There is no parallel markdown copy of the parameters. The cross-reference table that humans skim is *generated on demand* from the YAML (see "Cross-reference" below). If you want to compare a number against its source, the `source_quotes:` field on each YAML entry has the verbatim citation.

2. **No fabrication.** Every numerical value in the YAML is either backed by a `source_quotes:` entry citing a paper / table / page, OR is `null` with the gap explained in the entry's `notes:` field. If a paper was not accessible, the entry carries the citation but `null` numerics rather than secondary-source-derived guesses.

3. **Provenance is mandatory on edit.** When you change a value, you must also update its `source_quotes:` entry to cite the new source. Edits without a source-quote update should be rejected by review.

Attribution policy
------------------

Per spec.md §16.4:

- **Earliest priority date wins.** A cycler's `priority_date` is the earliest published date in any of its `first_published` or `corroborating_sources`. Attribution goes to those authors regardless of who later re-derives the cycler.
- **Never claim a published cycler as novel.** If `first_published` is set, the cycler is `known-reproduction` even when independently rediscovered by our search code.
- **Retroactive corrections allowed.** If a finder result is tagged `candidate-novel` and later-ingested literature matches its canonical signature, the pipeline auto-downgrades to `known-reproduction` and attaches the new citation.
- **Novelty claims require all four:** validation level V5 + no exact catalogue match + no unresolved probable-match + a documented literature review returning nothing.

Schema (extended) — 2026-05-31 non-heliocentric expansion
---------------------------------------------------------

The base schema in spec.md §16.1 was originally written for heliocentric
(Sun-primary) cyclers; all entries 1-37 in this file are heliocentric and
have no `primary:` field. On 2026-05-31 the catalogue was extended to
carry lunar (Earth-primary) and Jovian (Jupiter-primary) cyclers, with
hooks for future Saturnian (Saturn-primary) and other planet-centric
work. To keep the existing entries valid, a new optional top-level field
was added:

- `primary:` — optional string naming the gravitational primary the
  cycler orbits. Allowed values include `"Sun"` (heliocentric — default
  when the field is absent), `"Earth"` (lunar cyclers), `"Jupiter"`
  (Jovian / Galilean cyclers), `"Saturn"` (Saturnian cyclers; future),
  `"Mars"` (Martian-moon cyclers; future). Consumers MUST treat a
  missing `primary:` field as `"Sun"` to preserve backwards
  compatibility with entries 1-37.

Body codes are scoped to the primary:

- Heliocentric (`primary: "Sun"` or absent): `"V"`, `"E"`, `"M"` for
  Venus, Earth, Mars (as in the existing 15+ entries). Mercury would be
  `"Me"` if ever added.
- Earth-Moon (`primary: "Earth"`): `"E"` is the primary anchor / Earth;
  `"Moon"` is the orbiting natural satellite. Most lunar cyclers in the
  literature are CR3BP periodic orbits; the patched-conic + V_inf
  abstraction the rest of the catalogue uses is a poor fit for them,
  but the catalogue still records the citation + qualitative
  geometry. See entry notes for the model-mismatch flag.
- Jovian (`primary: "Jupiter"`): `"Io"`, `"Europa"`, `"Ganymede"`,
  `"Callisto"`. Full moon names are used (not single letters) to avoid
  collision with heliocentric planet codes.
- Saturnian (`primary: "Saturn"`; reserved): `"Mimas"`, `"Enceladus"`,
  `"Tethys"`, `"Dione"`, `"Rhea"`, `"Titan"`, etc.
- Martian (`primary: "Mars"`; reserved): `"Phobos"`, `"Deimos"`.

Some non-heliocentric entries are deliberately family-seed records (one
entry covering hundreds of catalogued members from a single paper) rather
than per-member records. These are flagged in `notes:` with phrasing
such as "family seed entry" and have `null` numeric fields with the
attribution citation preserved.

### `trajectory_regime:` field — 2026-06-01

A second optional top-level field, `trajectory_regime:`, was added on
2026-06-01 to make the *trajectory class* explicit:

- `ballistic` — Keplerian arcs + impulsive flybys; no deep-space thrust
  required to close the cycle. All but 2 of the current 235 entries are ballistic (the 2 powered Aldrin establishment variants excepted).
- `low-thrust` — requires continuous propulsion (Solar Electric, ion,
  nuclear electric, solar sail, etc.) over the transit legs to close.
  Mathematically: Lambert's problem no longer applies; the trajectory
  is solved with optimal control (Sims-Flanagan transcription, etc.).
  Not in catalogue today; reserved for v2 ingest.
- `manifold` — CR3BP invariant-manifold connections (Lyapunov / halo
  manifold hopping). The conserved quantity is the Jacobi constant,
  not V∞. Not in catalogue today; reserved for v2 ingest.

Consumers MUST treat a missing `trajectory_regime:` field as
`ballistic` (the v1-era default). Backfilled on all 43 existing
entries on 2026-06-01 so the field is now mandatory by convention even
though the default makes it formally optional.

The field is **trajectory geometry**, not propulsion hardware. A
`low-thrust` regime trajectory can be realised by ion / SEP / NEP /
solar-sail propulsion; the choice of propulsion is implementation, not
classification. (See OUTSTANDING.md §H for the rejected
`propulsion:` field discussion.)

When mixed-regime cataloguing begins (v2), use this field plus the
`primary:` field to filter searches and constrain matcher behaviour
(e.g. only match `ballistic` finder hits against `ballistic` catalogue
entries; signatures across regimes are incommensurable).

Schema v2 (2026-06-01) — six additive optional field categories
---------------------------------------------------------------

The 2026-06-01 schema rev adds six additive optional field categories
to every catalogue entry. All are **additive optional**: pre-v2
consumers that omit any of these fields treat the entry as if
`model_assumption = "circular-coplanar"` and every other v2 field is
`null`. None of the new fields participate in the §16.2 canonical
signature (so adding them does not invalidate existing matches against
literature). They are descriptive metadata for the M7 matcher's *pool
filter* stage (e.g. only match `cr3bp` hits against `cr3bp` literature)
and for downstream consumers (mission planners, public site
visualisations).

The convention `null` means **"not in source / not yet derived"**, NOT
"doesn't apply". When a field genuinely doesn't apply (e.g. CR3BP
entries have no Keplerian RAAN/ω/ν), the YAML field is still `null`
and the per-entry `notes:` block carries the model-mismatch flag —
matching the catalogue's existing convention from Arenstorf and
related entries.

### `model_assumption:` (top-level, string)

The model under which the entry's numeric values (V_inf, ToF, a, e,
peri, apo) were computed. Allowed values:

- `circular-coplanar` — default. The catalogue's main convention:
  Earth on a 1.000-AU circle, Mars on a 1.524-AU circle, both in the
  ecliptic plane, both unperturbed. Russell 2004 Tables 3.4/3.9-3.11,
  McConaghy 2002/2006 abstracts, S1L1 CPOM, Jones 2017 VEM-triple
  family seed, the entire 5-/6-synodic Russell catalogue. Every
  cycler whose published values are explicitly a "simplified model"
  reference.
- `analytic-ephemeris` — values come from STOUR or equivalent analytic
  ephemeris solver, not the circular-coplanar baseline. The two
  Rogers/Hughes/Longuski/Aldrin 2012 establishment variants (Aldrin
  4:3(2)- and 3:2(1)-) are explicitly tagged this way; their V_inf
  and ToF reflect a real-ephemeris reference epoch rather than the
  abstract repeating model.
- `cr3bp` — values come from a Circular Restricted Three-Body Problem
  formulation. Arenstorf periodic orbits, Genova/Aldrin 2015 lunar
  cycler, Wittal 2022 Earth-Moon family, Russell/Strange 2009
  Saturnian-multi-moon family. These entries' numeric V_inf and
  Keplerian orbital-element fields are mostly `null` because patched-
  conic abstractions don't apply; the entry exists as a citation +
  model-mismatch flag.

**Why this matters.** M7 should only match finder candidates against
literature with the *same* `model_assumption` — comparing a CR3BP
periodic orbit against a patched-conic cycler's signature is
meaningless, and would create false novelty claims. Use this field as
the matcher's pool prefilter (before the per-signature distance
calculation).

**Default when absent:** `circular-coplanar`. Pre-v2 entries that
predate this field (none currently — backfilled on all 235) should be
read as circular-coplanar.

**Backfill stats (initial v2 rev, 2026-06-01):**

| value                | count | source of classification                                |
| -------------------- | ----: | ------------------------------------------------------- |
| `circular-coplanar`  |   213 | default (Russell 2004 Tables, McConaghy, S1L1, Niehoff…)|
| `analytic-ephemeris` |     2 | Rogers 2012 Aldrin 4:3(2)- and 3:2(1)- establishment    |
| `cr3bp`              |     4 | Arenstorf 1963 + Genova/Aldrin 2015 + Wittal 2022 + Russell/Strange 2009 Saturnian |

### `delta_v_kms:` (top-level, float or null)

Per-cycle maintenance ΔV (km/s) required to close the cycle in the
entry's reference model. `0.0` for strict-ballistic cyclers (per
Russell's footnotes a + b: Aphelion Ratio ≥ 1.0 AND Turn Ratio ≥ 1.0;
also for entries the literature flags as "requires no propulsive
maneuvers", e.g. the McConaghy 2006 S1L1 abstract). Positive value for
near-ballistic entries (Russell's wider net, ARMIN ≥ 0.9 / TRMIN ≥
0.85) — derived per row from Russell's tabulated AR and TR via the M2
powered-flyby surrogate (`src/cyclerfinder/core/flyby.py`):
ΔV ≈ (V∞_E + V∞_M) × max(0, 1 − TR) + (1 − AR) × (V∞_E + V∞_M) × 0.025
(AR-correction only when AR < 1). `null` for entries where this number
isn't extractable from the source.

**Default when absent:** `null` (i.e. "we don't know"). Do not assume
`0.0` for entries that omit the field.

**Backfill stats (post-2026-06-01 Russell ΔV extraction):**

| value          | count | breakdown                                                                                        |
| -------------- | ----: | ------------------------------------------------------------------------------------------------ |
| `0.0`          |   132 | strict ballistic (Russell Tables 3.9–3.11 + 4.9–4.12 strict rows, McConaghy S1L1, Niehoff, …)    |
| positive       |    67 | Russell near-ballistic (Tables 4.x with AR ≥ 0.9 / TR ≥ 0.85); min 0.001, median 0.716, max 3.015 km/s |
| `null`         |    20 | 17 non-Russell entries (Aldrin, McConaghy, Jones, Hollister, CR3BP, etc.) + 3 Russell broad-class seeds |

### `v_infinity_leveraging_dv_kms:` (top-level, float or null)

Establishment-only ΔV (km/s) for the V-infinity leveraging manoeuvre
that gets the spacecraft *onto* the steady-state cycler. Specific to
Rogers 2012's Aldrin 4:3(2)- and 3:2(1)- variants and the broader
SnLm establishment context. `null` on all current entries pending
extraction of Rogers 2012 Tables 3/4 numbers.

**Default when absent:** `null`.

### `fleet_size:` (top-level, integer or null)

The minimum number of spacecraft required on the cycle to maintain the
specified cadence (typically: for a 2-synodic-period cycler one
vehicle gives a 4.27-yr cadence, but a 4-vehicle fleet gives a
~1-yr cadence). `null` on every current entry — the catalogue's
existing sources don't tabulate this systematically; it's typically
derivable from the period but not stated in the published numbers.

**Default when absent:** `null`.

### `flyby_mechanics:` (top-level, list of dicts or null)

Per-encounter geometry parallel to the existing
`vinf_kms_at_encounters:` list. Each element:

```yaml
- body: "E"
  turning_angle_deg: 60     # gravity-assist bend angle delta
  min_altitude_km: 200      # closest-approach altitude above body surface (convention: 200 km per Russell 2004 / Aldrin)
  rp_km: 6578               # closest-approach radius from body center (= body radius + min altitude)
```

`null` (or omitted) on all current 235 entries: while some Russell
entries cite turning-angle multisets like `[93, 93]` deg in their
notes, deriving per-encounter `turning_angle_deg` + `min_altitude_km`
mechanically requires the bend formula `sin(δ/2) = 1 / (1 + r_p ·
V_inf² / μ_planet)` and an `r_p` convention. Manual entry per
encounter is preferred over a derivation that may misalign with
Russell's own tabulated values. Backfill is **deferred**.

**Default when absent:** `null`.

### `orbit_elements.{periapse_km, apoapse_km}` (nested, float or null)

Parallel to the existing `orbit_elements.perihelion_au` /
`.aphelion_au`, but expressed in km from the gravitational primary.
Useful for non-heliocentric entries (Earth-Moon, Jovian, Saturnian)
whose natural unit is km, not AU. `null` for heliocentric entries
that already populate the `_au` fields; `null` for current
non-heliocentric entries because none of the catalogue's existing
sources tabulate them (Arenstorf, Genova, Wittal, Russell/Strange
Jovian/Saturnian all carry `null` orbit_elements per the
model-mismatch flag).

**Default when absent:** `null`.

### `orbit_elements.{raan_deg, arg_periapsis_deg, true_anomaly_deg, epoch_iso8601}` (nested, float / ISO-string or null)

The remaining 3D orbital orientation: Right Ascension of Ascending
Node (Ω), Argument of Periapsis (ω), True Anomaly at epoch (ν), and
the ISO-8601 epoch the anomaly is referenced to. Together with the
existing `a_au`, `e`, `inclination_deg` these complete the six-element
state plus epoch needed for ephemeris-mode propagation.

`null` on every current entry: none of the catalogue's literature
sources publish these (the circular-coplanar model has no need for
Ω/ω, and the published values are phase-independent). M6+ ephemeris-
mode optimisation will populate them naturally for entries that get
re-derived against the real ephemeris.

**Default when absent:** `null`.

Schema v3 (2026-06-01) — `trajectory{}` (OCM-aligned), `family{}`, `data_gaps[]`
--------------------------------------------------------------------------------

The 2026-06-01 v3 rev adds three additive optional top-level fields. The
full rationale — including the comparison to CCSDS OCM/OEM, TLE/OMM, JPL
SBDB/Horizons, MPC, the JPL CR3BP periodic-orbit catalog, and AstDyS
proper elements, and the decision to borrow OCM's trajectory model rather
than migrate to OCM — is in **spec.md §16.6**. This section is the
field-level reference.

As with v2, all three fields are **additive optional** and **none
participate in the §16.2 canonical signature** (`segments` contributes
only the per-leg `(a, e)` multiset that the legacy `legs[]` already
contributed, deduped to the same value; `family{}`, `maneuvers[]`, and
`data_gaps[]` are pure metadata). Existing matches against literature are
unaffected.

### `trajectory:` (top-level, dict) — OCM-aligned, supersedes flat `legs[]`

Reshapes the flat top-level `legs[]` into a self-describing trajectory
block mirroring CCSDS OCM's `TRAJ` (segments) + `MAN` (maneuvers)
decomposition, with an OCM-style metadata header:

```yaml
trajectory:
  center: "Sun"             # OCM CENTER_NAME — gravitational primary (mirrors top-level primary:)
  ref_frame: "ECLIPJ2000"   # OCM REF_FRAME
  time_system: "TDB"        # OCM TIME_SYSTEM; null for idealized circular-coplanar
  epoch_tzero: null         # OCM EPOCH_TZERO; null = epoch-free family (the common case here)
  segments:                 # === OCM TRAJ === ordered conic arcs, one per leg
    - id: "out-em"
      from: "E"
      to: "M"
      traj_type: "keplerian-arc"   # keplerian-arc | cartesian-state
      tof_days: 154
      n_revs: 0
      branch: "single"             # single | low | high  (multi-rev branch selector)
      a_au: 1.30                   # optional per-segment osculating arc elements
      e: 0.257
  maneuvers:                # === OCM MAN === per-encounter flyby / ΔV mechanics
    - at_segment_boundary: ["out-em", "ret-me"]   # flyby between these two segments
      body: "M"
      type: "flyby-ballistic"      # flyby-ballistic | flyby-powered | launch | arrival
      dv_kms: 0.0
      turning_angle_deg: null      # absorbs the v2 flyby_mechanics fields
      periapsis_alt_km: null
```

- **Migration is lazy** (see "Backfill" below): the loader reads
  `trajectory.segments` when present and falls back to the legacy top-level
  `legs[]` otherwise, producing an identical signature either way. An
  entry on either form is valid on disk during the transition.
- `maneuvers[]` absorbs the v2 `flyby_mechanics:` and top-level
  `delta_v_kms:` placeholders for migrated entries (the v2 fields remain
  valid for un-migrated entries).
- **OCM export.** Entries with `model_assumption: analytic-ephemeris`
  (epoch-anchored) can be projected to canonical CCSDS OCM (KVN/XML) by a
  future exporter; idealized circular-coplanar families cannot emit a
  meaningful OCM and are not export targets.

**Default when absent:** consumers fall back to the legacy `legs[]`.

### `family:` (top-level, dict or null) — CR3BP-catalog-style grouping

First-class family linkage, replacing prose-only `notes:` cross-references
(e.g. between the McConaghy "notable" variant and the S1L1 cycler):

```yaml
family:
  id: "s1l1-em"                  # stable family slug
  name: "S1L1 Earth-Mars 2-synodic"
  nomenclature: "Russell-McConaghy SnLm"
  continuation_param: {name: "k_synodic", value: 2}   # along-family parameter (cf. Jacobi constant)
```

Not a signature input (members are told apart by their own signatures); a
grouping/navigation aid for the matcher's human-review stage and the
public site. Missing `family{}` = "ungrouped."

**Default when absent:** `null` (ungrouped).

### `data_gaps:` (top-level, list of dicts) — "we don't know it yet" register

The explicit, machine-readable, **sweepable** register of *known-unknowns*
— fields whose value is expected to exist but has not been filled. This is
the **significant distinction** from a bare `null`: a `null` may be either
"not applicable" or "not yet derived" (the legacy convention above, lines
124-129), whereas a `data_gaps[]` entry asserts specifically "a value
exists; we owe it; here is where to find it."

```yaml
data_gaps:
  - path: "trajectory.segments[ret-me].tof_days"   # dotted/keyed path into this entry
    kind: "unknown"            # unknown | uncertain | derive
    note: "return M->E leg ToF not tabulated in the abstract"
    source_hint: "McConaghy/Longuski/Byrnes 2002, AIAA 2002-4420, Table 2"
    todo_ref: "#54-backfill"
```

- `kind`: `unknown` (no value on hand), `uncertain` (provisional /
  single-source value present), `derive` (computable once a dependency
  lands — e.g. a leg `(a,e)` once the multi-rev solver closes the arc).
- **This does NOT reinterpret existing nulls.** A bare `null` with no
  `data_gaps[]` entry keeps its legacy meaning; genuinely not-applicable
  cases are still flagged in `notes:` (the CR3BP/coplanar pattern).
  `data_gaps[]` is an explicit register layered on top — promoting a vague
  null to a tracked, actionable known-unknown.
- **Sweep:** `cyclerfinder.data.catalog.find_data_gaps(catalog)` and its
  CLI enumerate every gap across the catalogue, so the gap inventory is a
  query, not a manual audit.
- **Lazy fill:** when a value arrives, populate the field, update its
  `source_quotes:` per the provenance rule (rule 3 above), and remove the
  matching `data_gaps[]` entry. The no-fabrication rule (rule 2) is
  unchanged.

**Default when absent:** `[]` (no tracked known-unknowns).

### Backfill (`legs[]` → `trajectory{}`)

Lazy and sweep-driven, **not** a big-bang rewrite (most entries have only
partial leg data; a bulk rewrite would write mostly `null`s and obscure
which nulls are genuine gaps). Plan in spec.md §16.6.5. In short: (1)
loader supports both forms; (2) migrate `s1l1-2syn-em-cpom` first as the
worked exemplar (known E→M segment populated, unknown return + Earth-loop
segments marked in `data_gaps[]`); (3) sweep to inventory remaining
`legs[]`-only entries; (4) migrate the rest opportunistically, source-
gated; (5) "done" when no entry retains top-level `legs[]`.

Schema v4 (2026-06-03) — `cycler_class`, frame-tagged elements, `invariants{}`, `cr3bp{}`, `period.basis`
---------------------------------------------------------------------------------------------------------

The v4 rev (2026-06-03) makes the record *structurally honest* about the kind
of orbit it holds, removing the implicit assumption that every entry is a single
repeating heliocentric Kepler ellipse. **All v4 fields are additive, optional,
and none participate in the §16.2 canonical signature.** Full rationale and
field-level reference in **spec.md §16.7**.

- **`cycler_class`** — the structural kind of orbit:
  `single-ellipse` (default, one repeating Kepler ellipse; e.g. S1L1, Aldrin),
  `multi-arc` (different ellipse per leg — the bulk of the Russell catalogue),
  `non-keplerian` (CR3BP / rotating-frame; Arenstorf, lunar/Jovian entries).
  A `multi-arc` or `non-keplerian` entry MUST NOT carry a non-null top-level
  `orbit_elements.a_au`/`e`.

- **`orbit_elements.reference_frame` + `.center`** — frame/units-tagged elements
  following JPL SBDB conventions. Heliocentric entries keep
  `frame: heliocentric-inertial` (default); planet-centric entries use
  `planetcentric-inertial` with `center` naming the primary.

- **`invariants{}`** — present on `multi-arc` entries only; carries the Russell
  cycle-level descriptors (`aphelion_ratio`, `transit_times_days`, `turn_ratio`)
  promoted from prose `notes` to first-class fields the validator can assert.

- **`orbit_elements.cr3bp{}`** — present on `non-keplerian` entries; mirrors the
  JPL three-body periodic-orbit catalog: `(jacobi_constant, period_nd,
  stability_index)` plus `state_nd` + `mass_ratio` + `lunit_km`/`tunit_s`.

- **`period.basis`** — optional list of `{pair,k}` beat relations for n-body
  (VEM+) cyclers whose period is the beat of several synodic pairs; the legacy
  flat `{pair,k,years}` remains valid as the 2-body special case.

Backfill is lazy and source-gated: tag `cycler_class` first (mechanical sweep),
then populate `invariants{}`/`cr3bp{}` opportunistically as sources are read.
2026-06-04 backfill status: `cycler_class` tagged on all 235 entries
(single-ellipse 28, multi-arc 201, non-keplerian 6); `invariants{}`/`cr3bp{}`
populated on entries where Russell 2004 Tables 4.9–4.13 provide the descriptors.

Schema v4.1 (2026-06-03) — `free_return_arcs[]`
------------------------------------------------

The v4.1 sub-rev adds `free_return_arcs[]` to record Russell's arc-type
decomposition separately from the OCM encounter-segment decomposition:

```yaml
free_return_arcs:   # null when no descriptor is available
  - arc_type: generic          # generic | half-rev | full-rev
    resonance: null            # M:N string for full-rev arcs; null for generic/half-rev
    tof_years: 1.4612          # TOF in years for generic/half-rev; null for full-rev
    raw_descriptor: "g(1.4612,526.02,Ll)"   # verbatim Russell token
```

These are Russell's Earth-to-Earth free-return arcs — a *different* decomposition
of the same orbit from the encounter-leg segments. The two MUST NOT be conflated
(no 1-to-1 correspondence). Full field semantics in **spec.md §16.7.7**. The
`data/catalogue.schema.json` JSON Schema enforces the shape; a `check-jsonschema`
pre-commit hook validates every YAML change against it.

Backfill coverage (2026-06-03): 12 entries with explicit descriptors from Russell
2004 Tables 4.9–4.13; 3 entries with incomplete descriptors gapped; all
`russell-ocampo-*` entries gapped (Russell Ch3 tables carry AR/TR summary, not
leg descriptors).

Schema v4.2 — segment center, TOF bounds, source ephemeris
----------------------------------------------------------

The v4.2 sub-rev adds three additive, optional, non-signature fields:
`trajectory.segments[].center` (free string, absent ⇒ `Sun`, for planet-centric
moon-tour segments — issue #76), `trajectory.segments[].tof_days_bounds`
(a published `[min, max]` day range for sources that state a range), and
top-level `source_ephemeris` (the ephemeris model — `DE405`, `DE430`, `STOUR
ephemeris` — the source's numbers were computed against).

```yaml
source_ephemeris: "DE430"
trajectory:
  segments:
    - id: out
      tof_days: 146           # point idealization
      tof_days_bounds: [161, 172]   # published range; NOT required to contain tof_days
      center: Sun
```

`tof_days_bounds` is deliberately **not** required to contain `tof_days`:
different sourced model framings of the same leg (e.g. Aldrin's circular-coplanar
146 d vs Rogers et al. 2012 STOUR 161–172 d) are both valid. Full field
semantics in **spec.md §16.7.9**. The JSON Schema enforces the shape; the
Python semantic gate enforces `min <= max`, non-empty `source_ephemeris`, and
the non-containment rule. No rows are backfilled in this rev (structure only).

Out-of-paradigm work
--------------------

Cyclers obtained via CR3BP invariant manifolds (Lyapunov / halo
manifolds) or low-thrust / solar-sail trajectories are NOT in scope for
the current patched-conic + gravity-assist matcher and are NOT added to
this YAML. They are flagged in `OUTSTANDING.md` §H so the matcher does
not falsely tag finder hits against them. Re-evaluate ingestion when /
if the project adopts those modelling paradigms (cf. spec §2 stretch
goals).

Cross-reference
---------------

A skimmable cross-reference table of all current catalogue entries can
be regenerated at any time by running:

```sh
uv run --with pyyaml python scripts/render-catalogue.py        # markdown
uv run --with pyyaml python scripts/render-catalogue.py --csv  # spreadsheet ingest
```

The script reads `catalogue.yaml` directly; there is no committed
cross-reference file (it would only drift).

Adding a new published cycler
-----------------------------

1. Read the source paper. Capture verbatim quotes for every numerical value.
2. Add an entry to `catalogue.yaml` following the existing schema (spec.md §16.1 `signature_fields` shape, plus the `primary:` extension above if non-heliocentric).
3. If the entry surfaces a contradiction with an existing entry, or has gaps requiring future work, add a section to `OUTSTANDING.md`.
4. Commit as `data: add <cycler name> (<citation>)`.

Editing an existing entry
-------------------------

Edit the value in `catalogue.yaml` AND its parallel entry in `source_quotes:`. Edits without a source-quote update should be rejected by review.

Recording a discovery
---------------------

When the finder produces a closing trajectory that survives the V0-V3
auto-gauntlet (per spec §14), it writes a new entry with
`source: this-project`. Full lifecycle is in spec.md §16.5. The
practical shape of a project-discovered entry:

```yaml
- id: cyc-vem-novel-0001              # auto-assigned; "cyc-" prefix marks project-discovered
  name: "VEM triple cycler — cyclerfinder discovery"
  source: this-project                # vs "literature" (seed entries) or "both" (re-derived literature)
  primary: "Sun"
  trajectory_regime: ballistic
  bodies: ["V", "E", "M"]
  sequence_canonical: "E-V-M-E-M"
  sense: "outbound"
  period: {pair: "VEM-beat", k: 1, years: 6.41}
  vinf_kms_at_encounters: [...]
  orbit_elements: {a_au: ..., e: ..., ...}
  legs: [...]
  source_quotes: {...}                # cite the discovery_run, not external literature

  # Attribution lifecycle
  first_published: null               # populated when published
  priority_date: "2026-06-15"         # discovery date; revised on publication
  our_status: candidate-novel         # candidate-novel → verified-novel (V5) → published

  # Auto-validation results
  validation:
    level: V3
    gates: {V0: {pass: true}, V1: {max_diff_mps: 0.0}, V2: {max_drift_km: 1.2e4}, V3: {horizon_tcm_mps: 120}}

  # Discovery provenance (separate from reproducibility — see spec §16.5)
  discovery_run:
    finder_version: "0.4.2"
    config_hash: "sha256:..."
    seed: 42
    run_id: "cyc-runs-2026-06-15-1430"
    cell_id: "VEM|E-V-M-E-M|k1|r00000|blllll"
    launch_epoch: "2032-11-02"
    discovery_date: "2026-06-15T14:32:18Z"

  # How to re-derive (may be overwritten by cleaner re-runs)
  reproducibility:
    finder_version: "0.4.2"
    config_hash: "sha256:..."
    seed: 42
    run_id: "cyc-runs-2026-06-15-1430"
    cell_id: "VEM|E-V-M-E-M|k1|r00000|blllll"
```

Transitions:

- `candidate-novel` → `verified-novel` requires V5 (human expert review;
  documented literature search returns nothing). V4 is high-fidelity
  external (GMAT / Tudat / pykep) and is a prerequisite for V5.
- `candidate-novel` → `known-reproduction` happens automatically if a
  later-ingested literature entry matches the canonical signature
  (per spec §16.3 matcher); the literature citation is attached and
  the entry is retagged.
- On publication: populate `first_published` with the human authors and
  publication metadata; revise `priority_date` to the publication date.
  This locks in priority — subsequent literature finds of the same
  cycler bow to our publication.

**Authorship rule (per project commit policy):** `first_published.authors`
on project-discovered entries lists **human contributors** with
substantive engineering / scientific input. No AI co-authorship.
Project self-citation (for people citing cyclerfinder as a tool) lives
in `CITATION.cff` at the repo root.

Audit
-----

The seed catalogue was compiled May–June 2026 from:

- The Russell 2004 dissertation (full PDF, UT Austin handle 2152/1253) —
  primary source for the 24-ballistic-cycler taxonomy and the Aldrin
  cycler's energy parameters. Tables 3.4–3.8 and 4.7–4.13 are the most
  detailed accessible primary source.
- Rogers/Hughes/Longuski/Aldrin 2012 (AIAA 2012-4746, cached PDF from
  engineering.purdue.edu) — primary source for the orbital elements
  table (Aldrin, VISIT-1, VISIT-2, S1L1, plus several other variants).
- Web search snippets quoting AIAA abstracts (McConaghy 2006, Jones 2017)
  — used where direct access was blocked by HTTP 403.
- NTRS metadata records (Hernandez 2017 Jovian, Genova/Aldrin 2015,
  Wittal 2022, Arenstorf 1963 reprint).
- Wikipedia's "Mars cycler" article and the spaceflighthistory blog
  summary of Niehoff 1985 — secondary corroboration for several values.

**Source-access caveats encountered:**

- All AIAA-hosted PDFs and abstract pages (arc.aiaa.org) returned **HTTP 403
  Forbidden** to the web-fetch tool used to compile this catalogue. Quoted
  AIAA abstracts come from secondary sources (search snippets, ResearchGate
  "Request PDF" landing pages) and from cached PDFs hosted at academic
  institutions rather than at AIAA.
- The Russell 2004 dissertation was successfully downloaded from UT
  Austin's open-access repository (<http://hdl.handle.net/2152/1253>).
- Niehoff's original 1985 SAIC presentation and the early Niehoff 1986
  AAS Paper 86-172 were **not accessible** in any digitised form; values
  attributed to "Niehoff 1985" come from Rogers et al. 2012 Table 1 and
  the spaceflighthistory blog summary, with the original Niehoff documents
  cited but not consulted directly.
- The McConaghy 2006 JSR paper (10.2514/1.15215) was inaccessible beyond
  its abstract; orbital elements (a, e, peri, apo) for the "Notable
  Two-Synodic" cycler are therefore `null` in the YAML. The McConaghy
  2005 Purdue PhD dissertation (e-Pubs AAI3166673) is the open-access
  alternative source for the broader SnLm taxonomy and is queued for
  future ingest.
- The Jones/Hernandez/Jesick 2017 VEM paper full text was inaccessible
  beyond the abstract; the entry is a family-seed pending member-level
  ingestion when the paper becomes accessible.

The compiler has NOT personally consulted: the Niehoff 1985 SAIC
presentation, the Niehoff 1986 AAS 86-172 paper, the McConaghy 2006 full
paper (only the abstract), the Russell-Ocampo 2004 JGCD or 2005 JSR
papers (only abstracts; the dissertation is the comprehensive treatment),
or the Jones/Hernandez/Jesick 2017 full paper (only the abstract). Every
numerical value in the YAML is grounded in one of the above sources;
gaps are flagged in `OUTSTANDING.md`.
