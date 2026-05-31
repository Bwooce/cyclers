# data/

Top-level project data files. Consumed by:
- M3 golden tests (when those land) — reads `seed_cyclers.yaml`.
- The `cyclers.space` static site repo — pulls `seed_cyclers.yaml` directly.
- M7's catalogue loader (`src/cyclerfinder/data/catalog.py`, future) — reads `seed_cyclers.yaml`.

Files
-----

- `seed_cyclers.yaml` — published-cycler seed catalogue with full attribution per spec.md §16. Edits go through the same process as code: change values, change source quotes, commit. See `docs/known-cyclers.md` for the human-readable companion document.

Adding a new published cycler
-----------------------------

1. Read the source paper. Capture verbatim quotes for every numerical value.
2. Add an entry to `seed_cyclers.yaml` following the existing schema (spec.md §16.1 `signature_fields` shape).
3. Add a section to `docs/known-cyclers.md` describing the entry's provenance and any ambiguities.
4. Update the cross-reference table in `docs/known-cyclers.md` §7.
5. Commit as `data: add <cycler name> (<citation>)`.

Editing an existing entry
-------------------------

Edit the value in `seed_cyclers.yaml` AND its parallel entry in `source_quotes:`. Edits without a source quote update are rejected by review.

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
this YAML. They are flagged in `docs/known-cyclers.md` under "Outstanding
questions" so the matcher does not falsely tag finder hits against
them.
