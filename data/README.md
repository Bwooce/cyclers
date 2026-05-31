# data/

Top-level project data files. Consumed by:
- M3 golden tests (when those land) — reads `seed_cyclers.yaml`.
- The `cyclers.space` static site repo — pulls `seed_cyclers.yaml` directly.
- M7's catalogue loader (`src/cyclerfinder/data/catalog.py`, future) — reads `seed_cyclers.yaml`.

Files
-----

- `seed_cyclers.yaml` — published-cycler seed catalogue with full attribution per spec.md §16. **Sole source of truth.** Edits go through the same process as code: change values, change source quotes, commit.
- `OUTSTANDING.md` — long-form research questions / source-access gaps / paradigm-mismatch flags log. The YAML's per-entry `notes:` field handles short-form caveats; OUTSTANDING handles the discussion threads.

Conventions
-----------

Three rules — non-negotiable:

1. **Single source of truth.** The YAML at `seed_cyclers.yaml` is the canonical data record. There is no parallel markdown copy of the parameters. The cross-reference table that humans skim is *generated on demand* from the YAML (see "Cross-reference" below). If you want to compare a number against its source, the `source_quotes:` field on each YAML entry has the verbatim citation.

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
  required to close the cycle. All current 43 entries are ballistic.
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
predate this field (none currently — backfilled on all 219) should be
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
maneuvers", e.g. the McConaghy 2006 S1L1 abstract). `null` for
near-ballistic entries (Russell's wider net, ARMIN ≥ 0.9 / TRMIN ≥
0.85) — actual values require row-by-row Russell ΔV extraction which
is deferred. `null` also for entries where this number isn't
extractable from the source.

**Default when absent:** `null` (i.e. "we don't know"). Do not assume
`0.0` for entries that omit the field.

**Backfill stats (initial v2 rev, 2026-06-01):** 117 entries `0.0`
(strict ballistic), 102 `null` (near-ballistic or undetermined).

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

`null` (or omitted) on all current 219 entries: while some Russell
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

The script reads `seed_cyclers.yaml` directly; there is no committed
cross-reference file (it would only drift).

Adding a new published cycler
-----------------------------

1. Read the source paper. Capture verbatim quotes for every numerical value.
2. Add an entry to `seed_cyclers.yaml` following the existing schema (spec.md §16.1 `signature_fields` shape, plus the `primary:` extension above if non-heliocentric).
3. If the entry surfaces a contradiction with an existing entry, or has gaps requiring future work, add a section to `OUTSTANDING.md`.
4. Commit as `data: add <cycler name> (<citation>)`.

Editing an existing entry
-------------------------

Edit the value in `seed_cyclers.yaml` AND its parallel entry in `source_quotes:`. Edits without a source-quote update should be rejected by review.

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
