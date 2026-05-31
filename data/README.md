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
