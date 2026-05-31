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
