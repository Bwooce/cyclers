# Digest: Gobetz and Doll 1968, "A Survey of Impulsive Trajectories", UARL Report G-910557-11 (corpus policy; map digest)

F. W. Gobetz and J. R. Doll, United Aircraft Research Laboratories (East Hartford, Connecticut), Report
G-910557-11, "A Survey of Impulsive Trajectories", final report under NASA contract NAS8-21091
(Impulsive Transfer Study). NTRS 19680026430. The title page gives the date as June 1968 and 124 pages. It
is signed by F. W. Gobetz and J. R. Doll and approved by J. L. Cooley (Acting Manager, Systems Analysis).
The foreword thanks A. Schwaniger and R. Burns of NASA Marshall and T. N. Edelbaum (Analytical Mechanics
Associates), who wrote Appendix III.
- Filed as `cyclers_pdf/papers/gobetz-doll-1968-survey-of-impulsive-trajectories-uarl-g-910557-11-ntrs-19680026430.pdf`
  (the NTRS scan as received, 129 PDF pages, md5 dbac4cc329400964d64f6a8f8c71904c). The OCR copy is
  `cyclers_pdf/papers/gobetz-doll-1968-survey-of-impulsive-trajectories-uarl-g-910557-11-ntrs-19680026430-ocr.pdf`.
- The NTRS scan is not text-free. It carries a letter-spaced text layer with many misreads. The new OCR
  layer was made with `ocrmypdf --force-ocr`, so that layer replaced it.
- Journal form: Gobetz and Doll, AIAA Journal 7(5):801-834, May 1969, DOI 10.2514/3.5231 (the DOI was
  confirmed against the Crossref record on 2026-10-11, which lists the same two authors and the title
  "A survey of impulsive trajectories."). The project does not hold the journal form. It is proposed as a
  wanted row in sec. 7.
- Tables transcribed to `data/sources/gobetz-doll-1968-tables.yaml`.
- I read the whole body (printed pp. 1-68), the three appendices, the glossary, the symbol list, the three
  tables, the figure pages and the bibliography text. I did not digitise any curve.

## 0. Verdict

**A survey of planet-centred orbit transfer, intercept and rendezvous with impulses. It contains no
cycler, no repeating or periodic trajectory, no free return and no swingby round trip.**
- The authors list interplanetary applications among the subjects not covered in the survey text
  (printed p.7, with computational techniques, correction maneuvers and terminal-phase rendezvous). The only heliocentric material is a few lines on Mars-Earth three-impulse
  rendezvous (ref 52) and the interplanetary list in Table I.
- The one idea that touches gravity assist is the hyperbola-to-hyperbola "free transfer". It is stated as
  a pure Keplerian fact about a single hyperbolic passage (equal incoming and outgoing speeds at infinity
  are joined with no propulsion, any turning angle by choice of passage radius). It is not applied to a
  planetary swingby in this report.
- Value to the project: the bibliography (Table I interplanetary list, and refs 21, 52, 190, 248, 297, 298)
  points at period-era round-trip and double-conic papers. Two of them are not on the wanted list.

## 1. Map of the survey (printed page numbers; PDF page = printed page + 5)

The report classifies by boundary condition (free orbit, fixed orbit, terminal, subterminal) and by
trajectory type (intercept, transfer, rendezvous). All results are time-open unless stated.

| Section (printed p.) | Trajectory classes | Key results | References it rests on |
|---|---|---|---|
| Summary, conclusions (1) | whole survey | Fixed-time trajectories, optimal multi-impulse trajectories and optimal fixed-time rendezvous are the three under-studied areas | |
| Introduction, method (2-6) | whole survey | 314 bibliography entries, fewer than 5 percent before 1958; 58 percent since 1963; about 10 percent never obtained; existing surveys are entries 38, 53, 179, 180, 148, 195, 198 | surveys 195 and 198 (Marchal and Marec, not in English) |
| Intercept: general (9-10) | one-impulse intercept of a radius, a point | Hohmann-type 180 deg grazing intercept; a quartic gives the point intercept; a two-impulse plane-change leg helps near 180 deg central angle | 22, 147, 270, 263, 249, 71, 11, 58, 82 |
| Intercept: disorbit (10-13) | one-impulse and two-impulse disorbit to a radius, entry speed, entry angle; noncoplanar | Grazing entry is optimal when angle is open; two regimes by entry speed; off-apogee disorbit can save up to 10 percent from near-circular orbits | 40, 186, 75, 270, 286, 266, 114, 221, 74, 274, 275, 8 |
| Transfer: free-orientation (14-16) | orbit-to-orbit, orbit-to-subterminal, subterminal-to-subterminal, with unconstrained orientation | Optimum inclination is zero; Hohmann type is the global optimum for two or more impulses (proved in ref 20); four or more impulses never optimal | 277, 276, 20, 193, 123-125, 114, 283, 219, 220, 112 |
| Transfer: unbounded orbits (16-21) | parabola and hyperbola links, one-, two-, three-impulse escape, transfer through infinity | Tangential pericenter escape; three-impulse escape beats two-impulse for large V-infinity; **hyperbola-to-hyperbola free transfer for equal V-infinity; up to six impulses when passage radius is bounded** | 199, 195, 121, 230, 155, 58, 217, 22, 92, 294, 32, 194, 83, 45, 296, 312 |
| Transfer: circle-to-circle (22-27) | coplanar and noncoplanar, 2- and 3-impulse | Hohmann proved globally optimal (20); bi-parabolic beats Hohmann above radius ratio 11.93876; bi-elliptic always beats it above 15.58172; noncoplanar: three modes (generalized Hohmann, three-impulse, through infinity) with a single map (Fig. 7); through infinity optimal above 60.185 deg | 109, 20, 6, 143, 290, 207, 176, 293, 97, 258, 107, 58, 101, 61, 115, 9, 241, 104, 17, 251, 195, 196, 239, 295, 291, 200, 81, 62 |
| Transfer: circle-to-ellipse (27-30) | coplanar (three configurations) and noncoplanar | Tangential apsidal impulses; Hohmann or through infinity; noncoplanar three-impulse exists | 215, 237, 50, 113, 213, 299, 56, 37, 35 |
| Transfer: coaxial ellipses (30-36) | coplanar and noncoplanar, axes aligned or opposed | Mode tests via Marchal's maps (Figs. 11-13); generalized Hohmann, finite three-impulse, through infinity; Table II and III summarise Fig. 13 | 237, 170, 192, 193, 195, 196, 299, 300, 301, 105, 106 |
| Transfer: special cases (36-45) | congruent orbits; neighboring orbits (linearised); cotangential; ascent | Congruent ellipses: symmetric two-impulse optimal below e = 0.53533, through infinity above; linear theory: impulses never exceed the number of end conditions; Contensou spool; ascent: two-impulse 90 deg or 180 deg transfer angle, 70.529 deg when r1 = r2 | 152, 237, 170, 56, 54, 55, 227, 224, 265, 45, 81, 205, 59, 66, 143, 290, 238, 147, 200, 62, 301, 302, 35, 153, 214, 288, 287, 296, 12, 114, 178, 303, 42, 314 |
| Transfer: terminal-to-terminal (45-47) | two terminal states, time-free and time-fixed | Octic (eighth-order) polynomial; two relative minima; hyperbolic transfer orbits can occur; fixed-time: intermediate impulses help in an inverse-square field | 84, 6, 272, 7, 174, 311, 150, 90, 70, 183 |
| Transfer: ellipse-to-ellipse (47-52) | coplanar time-open and noncoplanar | "Useful angle" under 12.5 deg; switching rules; modes F, FF, FR, FFR, FRR, through infinity; the three-impulse mode needs e1 + e2 above 1.712 | 192, 193, 195, 196, 212, 47, 81, 45, 34, 202, 203 |
| Rendezvous (52-61) | time-open, time-limited (circle-to-circle: Hohmann, bi-elliptic, semitangential; general), time-fixed (neighboring orbits, general orbits) | Time-open rendezvous costs the same as time-open transfer plus waiting (at most one synodic period for Hohmann); fixed-time four-impulse solutions exist for 1 < t < 2.5 reference periods (Fig. 22); **three-impulse Mars-Earth rendezvous beats two-impulse for long-time long-angle transfers (ref 52)** | 25, 57, 282, 250, 31, 231, 207, 142, 29, 308, 238, 302, 297, 298, 52, 202, 183, 313, 129 |
| Discussion of results (61-63) | | Recommends fixed-time trajectories, multi-impulse modes, optimal rendezvous; lists space rescue, near-Earth operations, disorbit, interplanetary probes and landers, and abort as the applications | |
| Appendix I: methods (64-65) | | Calculus of maxima, hodograph, calculus of variations, adaptive steepest descent, impulse-as-independent-variable, primer vector, sufficiency tests | 272, 170, 154, 54, 235, 202, 315, 45, 34, 156, 172, 183, 129, 313 |
| Appendix II: impulsive approximation (66-67) | finite versus impulsive thrust | Finite-thrust penalty from gravity gradient and thrust-direction change; a generalized impulsive approximation (ref 243) is advised for detailed analysis | 243, 163, 184, 292, 306, 87, 96 |
| Appendix III: singular arcs, by Edelbaum (68-69) | | Singular arcs are minimizing for some fixed-time coplanar cases (Robbins); Lawden spiral is not minimizing | 244, 225, 224, 265, 137, 168, 169, 172, 242 |
| Glossary, symbols (70-73) | | 24 defined terms; symbol list | |
| Tables I-III (74-76), Figs. 1-22 (77-98), Bibliography (B-1 to B-26) | | see secs. 2-3 | |

Classes asked about, one line each:
- One-impulse: intercept and disorbit sections, plus the cases noted under transfer where orbits intersect.
- Two-impulse: Hohmann and generalized Hohmann (sec. circle-to-circle onward).
- Three-impulse: bi-elliptic, noncoplanar three-impulse, three-impulse escape, fixed-time rendezvous.
- Rendezvous, plane changes, escape: covered as above.
- Interplanetary: excluded by scope, except ref 52 and the Table I list.
- Round trips, swingbys, free returns, periodic or repeating trajectories: **not treated.** The word "free
  transfer" appears only for hyperbola-to-hyperbola links (p.20). The only periodic-type object is the
  Lawden spiral in Appendix III, which is a singular arc and not a periodic orbit. Ref 315 (McCue and
  Bender) is a three-body two-impulse quasi-linearisation study (p.64), not a periodic-orbit paper.

## 2. Tables

Three tables, all in `data/sources/gobetz-doll-1968-tables.yaml`.
- **Table I (p.74)**: lists of bibliography entry numbers by use (not acquired 25 entries; reviewed but not
  referenced 53 entries; not in translation 9; terminal-phase rendezvous 13; correction maneuvers 11;
  interplanetary applications 10; computational techniques 7). This is a map of the bibliography, so it
  matters for sec. 7.
- **Table II (p.75)**: coordinates of named points in Fig. 13 (pericenter-ratio rho or sigma against
  inclination i2). Kept as printed.
- **Table III (p.76)**: which optimal transfer mode holds in each Fig. 13 zone for M and N.
- Quoted constants from the running text (11.93876, 15.58172, 0.08376, 60.185 deg, 0.53533, 70.529 deg) are
  in the same file under their own key.
- There are no numeric ephemeris, state or orbit-element tables. Nothing here feeds a catalogue row.

## 3. Figures (22 plates on printed pp.77-98; described from the plates and the text, not digitised)

1. Fig. 1: coplanar circle-to-point intercept, normalised delta-V against radius ratio r2/r1 (0 to 10) for
   central angles from 0 to 330 deg, with elliptic, nearly critical, parabolic and hyperbolic regions.
2. Fig. 2: disorbit to a small entry angle, minimum delta-V against radius ratio, with optimum thrust angle.
3. Fig. 3: one-, two- and three-impulse escape sketches from a circular orbit.
4. Fig. 4: three-impulse escape delta-V/Vc against V-infinity/Vc for several radius ratios (from the text,
   p.18).
5. Fig. 5: summary of Hohmann transfers: total delta-V and the two impulses against radius ratio, 0.01 to
   100 (plate read; maximum of the second impulse near radius ratio 10 is visible).
6. Fig. 6: bi-elliptic delta-V against conjunction ratio for several radius ratios, exterior and interior
   conjunctions, Hohmann as the dashed reference.
7. Fig. 7: complete map of circle-to-circle optimal modes in the plane of radius ratio and inclination,
   with contours of delta-V and of conjunction ratio.
8. Fig. 8: three circle-to-ellipse configurations (ellipse inside circle, circle inside ellipse,
   intersecting).
9. Fig. 9: circle-to-ellipse summary diagram from ref 299, with shaded through-infinity regions.
10. Fig. 10: the two candidate two-impulse coplanar coaxial transfers, axes aligned and opposed.
11. Fig. 11: coplanar aligned coaxial orbits, Hohmann against through-infinity boundary at pericenter ratio
    11.938.
12. Fig. 12: non-intersecting coaxial ellipses with axes opposed, five zones (Marchal).
13. Fig. 13: map of optimal mode between aligned coaxial orbits (generalized Hohmann, finite three-impulse,
    through infinity); its points are tabulated in Table II.
14. Fig. 14: coplanar congruent-orbit transfers, eccentricity against angle between pericenters, with the
    boundary between symmetric two-impulse and through infinity.
15. Fig. 15: three configurations of noncoplanar congruent orbits (line of nodes perpendicular to, along, or
    inclined to the axis of symmetry).
16. Fig. 16: the "useful angle" band of thrust directions between local horizontal and local tangent.
17. Fig. 17: test for transfers between intersecting coplanar ellipses (eccentricity e1 against e2).
18. Fig. 18: the same for non-intersecting coplanar ellipses; the three-impulse mode is accelerate,
    accelerate, brake.
19. Fig. 19: phase angles for rendezvous between coplanar circular orbits.
20. Fig. 20: bi-elliptic transfer radius required for rendezvous against phase angle, for constant radius
    ratio.
21. Fig. 21: equal-time rendezvous, delta-V of modified Hohmann against outer bi-elliptic, by radius ratio,
    phase angle and plane change.
22. Fig. 22: reachable final-state map for fixed-time rendezvous between neighbouring circular orbits, with
    zones for two-, three-, four-impulse and terminal-coast solutions (ref 238).
- Plates 4, 6, 7, 14 and 22 have curve content that the text layer does not carry. I confirmed their titles
  from the text and the plate headers only.

## 4. Relevance to this project

- **Cycler, periodic, repeating, free-return or Venus-swingby Mars round trip: none treated.** The survey is
  scoped to a single central body. The project catalogue is therefore not touched by it. No Earth-Mars or
  Earth-Venus row gains a source, a number or a check from this report.
- **Closest physical content**: hyperbola-to-hyperbola free transfer (p.20, ref 194, Marchal 1966). A
  hyperbolic passage with equal V-infinity before and after needs no propulsion and turns by an angle set by
  the passage radius. This is the Keplerian kernel of an unpowered flyby. It is stated with a planetary
  radius floor ("rayon planetaire non nul" in the title of ref 194) but without any swingby application.
  Its use for our swingby legs would need the original French paper, which is not held.
- **Lower-order help for the Earth-Mars rendezvous legs**: ref 52 (Doll and Gobetz 1966) reports
  three-impulse Mars-Earth rendezvous trajectories that beat two-impulse ones for long-time long-angle
  transfers and are less sensitive to launch date. That is a statement about impulsive intermediate burns,
  not ballistic cyclers. It is a possible source for the project's powered-leg comparisons only.
- **Which cited sources the project holds** (grep of CORPUS_INDEX.md and the wanted list, 2026-10-11):
  - Held: ref 36, Breakwell, Gillespie and Ross 1961 (ARS J 31(2):201-207).
  - On the wanted list: ref 248 Ross, "A Systematic Approach to the Study of Nonstop Interplanetary Round
    Trips" (wanted row 4); ref 21 Battin 1959 "The Determination of Round-Trip Planetary Reconnaissance
    Trajectories" (inside the catch-all wanted row 37 as "Battin 1959/1999").
  - Neither: everything else flagged in sec. 7.
- **New identifier for wanted row 4.** The wanted row for Ross says the AAS paper number was not found. This
  survey prints the venue as Advances in the Astronautical Sciences, Vol. 13, American Astronautical
  Society, New York, 1963. That is a verifiable proceedings volume to search. The lead may want to add it.
- **Value for a constraint check**: none. The report has no numbers that the V-gauntlet or the catalogue
  could use. Its critical constants (11.93876, 15.58172, 60.185 deg) are Hohmann/bi-elliptic facts for a
  single central body.
- **Possible use as a negative-results note.** The authors' own conclusion (p.61) says that fixed-time
  multi-impulse interplanetary problems have only been studied numerically case by case. That supports the
  project's view that no general analytic theory of such legs existed in 1968.

## 5. Notes on the text

- The introduction says the bibliography has 314 entries. The printed list runs to entry 315 (McCue and
  Bender, NASA CR-77878). I did not resolve the one-entry difference.
- Refs 195 and 198 (Marchal, and Marchal, Marec and Winn) are the two recent surveys of impulsive transfer
  that the authors single out. Neither was available in English translation to them.

## 6. OCR and three-witness record

- Witnesses: (1) the NTRS text layer, letter-spaced and read after joining the spaced letters; (2) a new
  tesseract layer from `ocrmypdf -l eng --force-ocr --rotate-pages --deskew` on the scan (the scan already
  had the NTRS layer, so plain OCR refused; the copy was then recompressed with Ghostscript, keeping the
  text, from 347 MB to 40 MB); (3) the page render, read by eye (110 dpi for Tables I-III, 250 dpi for
  Table II, 130 dpi for the constants checked below).
- Table I (170 numbers in 7 lists): the render and both text layers agree on every list length. The only
  misreads were in single numbers.

  | Number | NTRS layer | New OCR | Render |
  |---|---|---|---|
  | 120 | "EO" | 120 | 120 |
  | 133 | "1.33" | 133 | 133 |
  | 140 | 140 | "lho" | 140 |
  | 190 | "1-90" | 190 | 190 |

- Table II (21 values): each value is correct in at least one text layer and on the render; none is
  contradicted by a layer except by an evident garble.

  | Value | NTRS layer | New OCR | Render |
  |---|---|---|---|
  | 0.33333 | "0033333" | 0.33333 | 0.33333 |
  | 0.28942 | 0.28942 | "0,289h2" | 0.28942 |
  | 55.6 | "55 96" | 55.6 | 55.6 |
  | 0.447 | 0.447 | "O.447" | 0.447 |
  | 0.707 | "0 707" | 0.707 | 0.707 |
  | 60.185 | "60. a85" | 60.185 | 60.185 |

- Table III: the zone strings and the mode names agree on the render. Both layers garble the primes and
  subscripts on E1, Q, O and R, so the zone strings in the yaml are taken from the render.
- Constants quoted in the text: 11.93876, 15.58172, 0.08376, 60.185, 0.53533, 0.1505, 1.712, 8.7967,
  1.07067 appear in both layers (0.08376 as "0,08376" in the NTRS layer; 11.93876 once as "1143876" there).
  70.529 and the penalty 0.001615 with its squared lambda factor were read on the render because one layer
  misses them.
- Not three-witnessed: the bibliography entries in sec. 7. They were read from the two text layers and the
  entry texts agree between layers except where noted below. Page ranges and volume numbers should be
  checked against the render before a wanted row is searched on them.
  - Entry 285 (Vertregt, BIS Journal, vol. 16, March-April 1958, pp. 326-354): the issue number reads 79 in
    the NTRS layer and 19 in the new OCR. It is unresolved.
  - Entry numbers in the new OCR are garbled for some digits ("el" for 21, "ell" for 211), so the entries
    were matched by author and title, not by number alone.

## 7. Citation mining

The bibliography has 315 printed entries. Only entries that bear on cyclers, round trips, swingbys, free
returns or periodic orbits are listed. Held means in CORPUS_INDEX.md; wanted means on
`2026-10-05-960-wanted-papers.md`. Identifiers are as printed in the report and are not checked unless
stated.

| Ref | Item as printed | Status | Why it matters |
|---|---|---|---|
| 36 | Breakwell, Gillespie and Ross, "Researches in Interplanetary Transfer", ARS Journal 31(2), Feb 1961, pp. 201-207 | Held | Already digested |
| 248 | Ross, S. E., "A Systematic Approach to the Study of Nonstop Interplanetary Round Trips", Advances in the Astronautical Sciences, Vol. 13, American Astronautical Society, New York, 1963 | Wanted (row 4) | Nonstop round trips, directly on topic. This report gives a concrete venue that row 4 lacks |
| 21 | Battin, R. H., "The Determination of Round-Trip Planetary Reconnaissance Trajectories", Journal of the Aerospace Sciences 26, Sept 1959, pp. 545-567; also IAS 27th Annual Meeting preprint 59-31, Jan 1959 | Wanted only as "Battin 1959/1999" inside catch-all row 37; row 37 does not say which Battin 1959 item | Round-trip free-return type geometry from 1959; a candidate for its own row |
| 190 | MacKay, J. S. and Rossa, L. G., "A Variational Method for the Optimization of Interplanetary Round-Trip Trajectories", NASA TN D-1660, 1963 | Neither | Round-trip optimisation; NASA TN, so probably on NTRS |
| 247 | Ross, S. E. et al., "Space Flight Handbook, Vol. 3, Planetary Flight Handbook", NASA SP-35, 1963 | Neither | Handbook of planetary trajectory data from the same period as the round-trip papers |
| 52 | Doll, J. R. and Gobetz, F. W., "Three-Impulse Interplanetary Rendezvous Trajectories", Proceedings of the Southeastern Symposium on Missiles and Aerospace Vehicle Sciences, Huntsville, Vol. II, Dec 1966, pp. 55-1 to 55-14 | Neither | Mars-Earth three-impulse rendezvous, the only heliocentric result in the survey |
| 297, 298 | Willis, E. A., "New Class of Optimal Interplanetary Trajectories with Specified Trip Time", AIAA Paper 65-66, Jan 1965; and "Optimization of Double Conic Interplanetary Trajectories", NASA TN D-3184, Jan 1966 | Neither | Patched-conic trajectory optimisation with fixed trip time; the method behind entry 52 |
| 67 | Escobal, P. R. and Crichton, G. A., "Minimum Time Interplanetary Trajectories", Journal of the Astronautical Sciences 14(5), Sept-Oct 1967, pp. 212-217 | Neither | Interplanetary transfer family, low priority |
| 255 | Schechter, H. B., "Impulsive Interplanetary Transfers for Prescribed Launch Date", ARS Journal 32(11), Nov 1962, pp. 1716-1723 | Neither | Launch-date constrained transfers, low priority |
| 194 | Marchal, C., "Transferts Optimaux Entre Orbites Hyperboliques (Rayon Planetaire Non Nul)", 17th International Astronautical Congress, Madrid, Oct 1966; "to appear as a NASA translation" | Neither | The free hyperbola-to-hyperbola transfer, the Keplerian kernel of a flyby; French, low priority |
| 211, 132, 285, 32 | Moeckel, NASA TN D-80 (1959); Karrenberg and Arthur, ARS Preprint 870-59 (1959); Vertregt, "Interplanetary Orbits", BIS Journal 16, Mar-Apr 1958, pp. 326-354; Bossart, "Techniques for Departure and Return in Interplanetary Flight", Aerospace Engineering 17(10), Oct 1958, pp. 44-52 | Neither | Early interplanetary orbit papers, low priority |
| 315 | McCue, G. A. and Bender, D. F., "Satellite Orbital Transfer Studies", NASA CR-77878, Aug 1966 | Neither | Three-body two-impulse study by quasi-linearisation (p.64); not a periodic-orbit paper |
| 21, 36, 190, 211, 247, 248 | Table I "Interplanetary Applications" list (36, 67, 117, 132, 162, 190, 211, 247, 248, 285) | See rows above | Entries 117 (Il'in, Cosmic Research 2(5), 1964) and 162 (Lawden, "Transfer Between Circular Orbits", Jet Propulsion 26, July 1956) are planet-centred or coplanar-circular and are not proposed |

No entry in the bibliography is on periodic orbits of the three-body problem, free returns, swingby round
trips or repeating Earth-Mars or Earth-Venus trajectories. The 315 entries are impulsive-transfer papers.

### Proposed wanted rows (the lead decides which go on the list)

| Priority | Item | Identifier |
|---|---|---|
| 1 | Gobetz, F. W. and Doll, J. R. (1969), "A survey of impulsive trajectories", AIAA Journal 7(5):801-834 (journal form of this report; not held) | DOI 10.2514/3.5231, Crossref-confirmed 2026-10-11 (title, both authors, volume, issue, pages, May 1969) |
| 2 | MacKay, J. S. and Rossa, L. G. (1963), "A Variational Method for the Optimization of Interplanetary Round-Trip Trajectories", NASA TN D-1660 | Report number as printed; NTRS not checked |
| 3 | Ross, S. E. et al. (1963), "Space Flight Handbook, Vol. 3, Planetary Flight Handbook", NASA SP-35 | Report number as printed; NTRS not checked |
| 4 | Battin, R. H. (1959), "The Determination of Round-Trip Planetary Reconnaissance Trajectories", J. Aerospace Sciences 26:545-567 | Not checked (needs its own row, split from row 37) |
| 5 | Willis, E. A. (1966), "Optimization of Double Conic Interplanetary Trajectories", NASA TN D-3184 | Report number as printed; NTRS not checked |
| 6 | Doll, J. R. and Gobetz, F. W. (1966), "Three-Impulse Interplanetary Rendezvous Trajectories", Southeastern Symposium on Missiles and Aerospace Vehicle Sciences, Huntsville, Vol. II, pp. 55-1 to 55-14 | None found in the report |
| - | Wanted row 4 (Ross 1963): add venue "Advances in the Astronautical Sciences, Vol. 13, AAS, New York, 1963" | As printed in this report's entry 248 |

Catalogue implications: none. This report has no cycler, free-return or swingby content.
