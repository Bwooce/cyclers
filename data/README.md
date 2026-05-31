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
