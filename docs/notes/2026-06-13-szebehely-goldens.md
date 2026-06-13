# Szebehely 1967 libration-point / Jacobi goldens wired (#241)

**Date:** 2026-06-13
**Source:** V. Szebehely, *Theory of Orbits: The Restricted Problem of Three
Bodies*, Academic Press, 1967 (cited by book page / appendix / table).
**Mining note:** `docs/notes/2026-06-13-szebehely-1967-theory-of-orbits-mining.md` (#185).
**Test file:** `tests/core/test_szebehely_goldens.py`.

Wires SOURCED, non-JPL CR3BP goldens into the core suite as INDEPENDENT
cross-check anchors for `core/cr3bp.py` (`jacobi_constant`). No libration-point
*solver* exists in the codebase, so the anchors compare the **Jacobi constant
at the printed libration-point location** (the location is the published input;
C is the published expected output). No catalogue writeback.

## Convention conversions applied (from the mining note "Convention trap")

1. **Frame mirror.** Szebehely puts the larger mass 1-mu at +mu and the smaller
   at mu-1; we put the larger at -mu and the smaller at 1-mu. Mirror is `x -> -x`.
   C is mirror-invariant; collinear x-coordinates are negated before use.
2. **Two Jacobi constants.** Our `jacobi_constant()` at v=0 equals C_bar
   (= 2*Omega_bar, no mu(1-mu) term). The Appendix I-III tables print the
   STANDARD C, so we compare via `C_bar = C - mu*(1-mu)`. Table III (p.457)
   prints both C_bar and C and needs no conversion.

## What was wired (real numbers; our value vs Szebehely)

### Triangular points (exact-location anchors)
- **Table III (p.457), mu=1/2:** L4/L5 at x=0, y=+-sqrt(3)/2. Our C_bar =
  2.7500000000 vs printed C_bar = 2.7500 (and C = 3.0000). Exact.
- **p.451 identity C(L4,5) = 3 for all mu:** L4 at our x=1/2-mu, y=sqrt(3)/2
  (r1=r2=1). Our C_bar + mu(1-mu) = **3.0000000000** for every tested mu
  (0.01215068, 0.012, 0.1, 0.2, 0.3, 0.5) -> matches to <1e-12.
  - NOTE: the first draft mistakenly placed L4 at x=0 for all mu and saw C=3.0616
    at mu=0.1 -- that was MY location bug, not a stack defect; L4 sits at the
    equilateral apex x=1/2-mu, not the origin (except mu=1/2). Fixed.

### Collinear point, mu=1/2 (Table III, p.457; x to 4 sig figs)
- **L2 (inner), x=0 exact:** our C_bar = 4.0000000000 vs 4.0000 (C 4.2500). Exact.
- **L1/L3 (outer), his x=-+1.1984:** our C_bar = 3.4568 vs 3.4568 (C 3.7068),
  to the 4 printed figures (location-limited; asserted abs=2e-3).

### Collinear points, Earth-Moon range (App I.D/II.D/III.D, pp.216/220/224)
High-precision: location printed to ~13 sig figs, standard C printed.
Compared as C_bar = C - mu(1-mu) at the mirror-mapped location.

| App | mu | our C_bar vs printed (after conversion) | (ours - printed C) | status |
|-----|------|------------------------------------------|--------------------|--------|
| I.D   | 0.0120 | match to ~3e-11 | +2.9e-11 | TIGHT (abs=1e-9) |
| I.D   | 0.0123 | match to ~1e-11 | -9.2e-12 | TIGHT (abs=1e-9) |
| I.D   | 0.0121 | off by ~3.0e-7   | +3.0e-7  | LOOSE (abs=5e-7) |
| II.D  | 0.0120 | match to ~6e-12 | +6.2e-12 | TIGHT (abs=1e-9) |
| II.D  | 0.0121 | off by ~1.1e-7   | +1.1e-7  | LOOSE (abs=5e-7) |
| III.D | 0.0120 | match to ~4e-11 | -3.8e-11 | TIGHT (abs=1e-9) |
| III.D | 0.0121 | off by ~2.0e-7   | -2.0e-7  | LOOSE (abs=5e-7) |

## Discrepancies found (honest, not fudged)

### A. The mu=0.0121 rows are internally inconsistent (~1e-7)
All three mu=0.0121 collinear rows (App I/II/III) disagree with our code by
1e-7..3e-7, while EVERY mu=0.0120 and mu=0.0123 row agrees to ~1e-11. The
decisive test: recomputing the standard C **directly from the printed x in
Szebehely's own frame** (C = x^2 + 2(1-mu)/r1 + 2mu/r2 + mu(1-mu)) reproduces the
printed C to ~1e-11 for the 0.0120/0.0123 rows but gives the SAME ~1e-7-off value
our code does for the 0.0121 rows. So for mu=0.0121 the printed x and printed C
**do not satisfy Szebehely's own definition of C** -- they cannot both be exact.
- **Read:** most likely a last-digits transcription slip in the mining-note
  capture of the mu=0.0121 rows (the surrounding rows are perfect). Could also be
  original-print rounding at that single mu. Either way our stack is vindicated by
  the two bracketing mu values matching to 1e-11.
- **Disposition:** kept the mu=0.0121 rows as a LOOSE cross-check at abs=5e-7 with
  the EXPECTED value left at Szebehely's printed C (never adjusted toward ours);
  the relaxed tolerance documents the gap. **Re-check against the PDF** (App I.D
  p.216 / II.D p.220 / III.D p.224, mu=0.0121 row) to decide note-typo vs
  print-rounding.

### B. Appendix IV C(x) sample could not be reproduced -- test DROPPED
The single App IV sample in the mining note (mu=0.1, x=-1.0 -> C = 4.97029_60396)
does not reproduce in either frame: the standard C at x=-1.0, mu=0.1 in
Szebehely's own frame is 4.7263636..., and our code matches 4.7263636... after
the mirror. 4.97029 is ~0.24 away and the algebra forces 4.72636, so the sample
is a suspected note transcription error. A golden cannot rest on an unverifiable
number, so no App IV smoke test was wired. **Re-mine App IV (pp.226-229) before
adding a C(x) zero-velocity anchor.**

## Outcome

11 sourced goldens passing. The CR3BP `jacobi_constant` is now cross-checked to
~1e-11 against a 1967 classic, independent of the JPL periodic-orbit oracle, at
five Earth-Moon-range collinear points and the mu=1/2 Copenhagen L1-L5 set, plus
the exact C(L4,5)=3 identity across six mu. Two documented transcription
discrepancies (the mu=0.0121 rows; the App IV sample) are flagged for a PDF
re-check and are NOT treated as stack defects -- the bracketing exact rows show
our equations are correct.
