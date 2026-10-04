# Everything in this directory was computed in a model whose Sun turns the wrong way

Task `#891` (2026-10-04): until that date `core/bcr4bp.py` and `search/sun_forced_periodic_884.py`
advanced the Sun counter-clockwise in the Earth-Moon rotating frame. The real Sun moves clockwise
there. Every file here (controls, families at nonzero Sun mass, Melnikov scans, continuations,
verification, summary) belongs to that non-physical model and must not be used or cited as a
result about the Sun-Earth-Moon system. The three-body family walks under `families/` do not
involve the Sun and are unaffected.

The files are kept as the record of what was computed. A rerun with the corrected sense will be
written to a new directory. See `docs/notes/2026-10-04-884-sun-forced-em-cyclers.md` (notice at the
top) and `docs/notes/2026-10-04-884-adversarial-review.md`.
