# Digest — Harch, Carcich et al. 2016, "Accommodating Navigation Uncertainties in the Pluto Encounter Sequence Design"

**SpaceOps 2016 Conference** (Daejeon, Korea, 16-20 May 2016), AIAA. DOI
`10.2514/6.2016-2623`. (The book-chapter version is *Space Operations: Contributions
from the Global Community*, Springer 2017, DOI `10.1007/978-3-319-51941-8_21`.)
SwRI + JHU/APL + KinetX + JPL + NASA Ames; many authors incl. Finley, Olkin, Stern,
Young. Corpus file (private `cyclers_pdf/papers/`):
`harch-bhaskaran-2017-accommodating-navigation-uncertainties-pluto-encounter-sequence-design-space-operations-springer.pdf`.
Processed 2026-06-23 (#429 follow-on, user-supplied).

## What it is
A **navigation/operations methodology** paper on the New Horizons **Pluto flyby**
(2015-07-14): how the encounter command sequence was built to absorb navigation
targeting uncertainty (the "Late Update" pointing/timing-correction process, B-plane
delivery vs knowledge uncertainty, error modeling). It is **not** an orbital-tour or
trajectory-design paper — NH was a single flyby, not an orbiter. (This is the paper
Crossref returns as the only other Finley-coauthored Pluto item; it does NOT carry
Charon gravity-assist tour periapsis numbers — those would have been in the unpublished
Finley "Orbital Tour" paper. See `2026-06-15-279-finley-paper-disposition.md`.)

## Values relevant to us (#429 — observed-flown only, NOT design floors)

From **Fig. 1 (Nominal Pluto Flyby Trajectory)**, the flown New Horizons encounter
geometry (2015-07-14):

- **Pluto closest approach: 11:50:00 UTC, 13,695 km from center.** Minus Pluto mean
  radius 1188.3 km → **≈12,507 km altitude**. (Paper: "minimum distance to the center
  of Pluto was chosen to be 13695 km" — a single designed/flown C/A, set for imaging
  resolution + atmospheric/occultation geometry, NOT a physical floor.)
- **Charon closest approach: 12:04:00 UTC, 29,432 km from center.** Minus Charon mean
  radius 606 km → **≈28,826 km altitude.**
- Geometry context: flyby speed ~14 km/s; diametric Earth occultation by Pluto; solar +
  Earth occultations of both Pluto and Charon; Charon orbit-normal 47° to trajectory.

**Discipline:** both are `observed-flown` (single realised C/A, carries nav/science
margin) — the same class as the Voyager-2 values, **NOT design floors**. They do not
change any `safe_alt_km`. Charon had no reference-table entry at all, so this is its
first sourced datapoint (provenance only). Pluto's sourced *design* floor stays 100 km
(Stern 2020); this NH value is a flown cross-check well above it.

## Catalogue / task implications
- **#429:** Charon + Pluto now carry an `observed-flown` NH datapoint in
  `flyby_altitude_references.yaml`. The Charon/small-Pluto-moon **design-floor** gap
  REMAINS (no published design floor; the Finley tour paper that would have it was never
  published — #279).
- Methodology content (Late Update, B-plane error modeling) is not catalogue-relevant.
  Not an mga_tour/cycler row (single flyby).
