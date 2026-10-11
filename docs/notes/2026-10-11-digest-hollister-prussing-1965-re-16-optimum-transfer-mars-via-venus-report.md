# Digest: Hollister and Prussing 1965, "Optimum Transfer to Mars via Venus", MIT Experimental Astronomy Laboratory report RE-16 (#960 wanted-list row 8, report form)

Walter M. Hollister (Deputy Director) and John E. Prussing (Research Assistant, per the title-page footnotes; the page also carries a middle line "Assistant Professor, Department of Aeronautics and Astronautics" for the
second author as printed in the text layer, which of the three titles belongs to whom is not clear from the OCR), Experimental Astronomy Laboratory, MIT, Cambridge, Massachusetts 02139.
Report RE-16, dated April 1965 on the abstract page; prepared under DSR Project 5007 with NASA research grant NsG 254-62. NTRS 19660006110.
Filed as `cyclers_pdf/papers/hollister-prussing-1965-optimum-transfer-mars-via-venus-mit-eal-re-16-ntrs-19660006110.pdf`
(md5 `8d7f904f1c9425c21a4114563d2140b4`, 42 pp., a 300 dpi scan with an Acrobat Capture 3.0 text layer; original upload not modified).
The wanted-list row's date is "1966, Astronautica Acta 12(2)" for the journal form; the report is the April 1965 source of the same work (the conference form AIAA 65-700, held, is dated September 1965).

Text layer: the brief said "almost no text layer". It has 8,061 words, but spaced-letter OCR errors (for example "t r a n s f e r") and garbled equations. An ocrmypdf force-OCR copy
(`--force-ocr -O 1 -l eng`, md5 2af2e0a34c70c8663d9f1d40010a53f5, 30 MB, 6,552 words) read worse overall, so it was NOT filed; the original, 1.4 MB, is the corpus copy. The Tesseract copy served as a
second text witness for the Fig. 4 numbers only (data file).
Status: READ in full (text layer; equations 2 to 7 read on page renders; Figs. 1, 4, 5, 6, 7 viewed; Fig. 8 and 9 seen at low resolution only). Numbers: `data/sources/hollister-prussing-1965-re16-figure-values.yaml`.
COMPUTED = our arithmetic 2026-10-11 (script and output beside the PDF in cyclers_pdf/papers as `...-check_fig4.py` and `.out`; method and results also in the data file `check:` block). INFERRED = our reading.
Page numbers are the report's printed numbers.
Related digests, not repeated: `2026-10-05-digest-hollister-prussing-1965-optimum-transfer-mars-via-venus.md` (the held AIAA 65-700 form),
`2026-10-11-digest-hollister-1966-predicting-launch-dates-mars-via-venus-jsr.md`, `2026-10-07-digest-hollister-1963-mission-manned-expedition-mars-scd-thesis.md`,
`2026-10-05-digest-hollister-1969-periodic-orbits-interplanetary-flight.md`.

## 0. Verdict

**RE-16 is the same text as the held AIAA 65-700 paper.** Every section, the seven numbered equations plus the common-peripoint equations, the six conclusions, the fifteen references and the nine figures are the
same. The held AIAA copy is a re-set reprint of this report (its abstract still reads "This report describes ..." and it carries the same Experimental Astronomy Laboratory cover; 15 landscape pages against RE-16's
42). The report adds only: a title page, an acknowledgment of the NASA grant, a table of contents, a list of figures, and the nine figures at full size (the Fig. 5 fold-out is on four 8.5 by 11 in sheets). It has no tables,
no appendices, no extra launch years and no additional derivations beyond the AIAA text. So the earlier digest's statement that the paper has no tables is confirmed, and the wanted-list suspicion that the report
might carry more tables and launch years does not hold. It is a one-shot mission study. Nothing in it is periodic.

What the report does add for the project is the numbers printed on its figures, which the held AIAA text layer lacks (the figure annotations were not in that layer): the Fig. 4 optimum and common-peripoint delta-V,
and the via-Venus minimum launch velocities of Fig. 5. They are in the data file. The Fig. 4 optimum reproduces independently (sec. 3).

## 1. Sections (pp.1-28)

- Introduction (p.1): Venus flybys enroute to a Mars stopover mission save energy over direct flights (refs 5, 6, 7); a thrust impulse at Venus may save more (refs 8, 9, 10). The study compares three Earth-Mars
  free-fall transfers: direct, Venus flyby without thrust, and a flyby with a thrust impulse at Venus.
- When Venus is available (pp.2-3, Fig. 1): a Venus encounter can change the heliocentric velocity by up to the circular satellite velocity at periapsis, above 20,000 ft/s (about 6 km/s). Two principles:
  (i) bodies in near-circular orbits inside 1 AU get ahead of Earth and outside 1 AU fall behind, so a fast Mars round trip must spend time inside Earth's orbit; (ii) low-energy transfers have the planetary alignment
  about half way between launch and arrival. So the Sun-Venus-Mars alignment, plotted relative to Earth, predicts whether Venus is available (Fig. 1 viewed: a ring diagram of alignments Feb 1971 to Apr 1981 with arcs
  "Venus available enroute to Mars" and "returning from Mars"). Venus is available on one leg at every Mars opposition and both ways at every third.
- Optimisation of dates (pp.4-5): a flyby links two trajectories (launch planet to flyby planet, flyby planet to destination). A thrusted flyby is needed if the inbound and outbound hyperbolic speeds differ or the pure
  flyby would go below the surface. The authors compute every flyby date for each launch and arrival date directly, not by a gradient optimiser, because there is more than one local minimum.
- General description and details (pp.6-18): the inbound and outbound hyperbolic excess velocity vectors are fixed by the launch, flyby and arrival dates. The impulse optimisation is planar (perfect guidance assumed).
  Lengths are in planet radii and speeds in circular surface speeds. Eq. (1): the peripoint radius in terms of the angular momentum h and V_H (not transcribed; the equation is garbled in the text layer and was not read on
  an image). Eq. (2): Battin's hyperbolic injection velocity, V = (1/2)(sqrt(1 + 4/(V_H^2 r (1 + cos theta))) + 1) V_H + (V_H/2r)(sqrt(...) - 1) r. Eq. (3): the same for the outbound hyperbola. Eq. (4): the inbound hyperbola
  with signs changed (the velocity -V_I achieves -V_HI). Then delta-V = V_O - V_I with theta_O = pi + A - theta_I. Eqs. (5) and (6): x and y components of the impulse in a frame with its origin at the planet and
  x parallel to the inbound asymptote; Eq. (7): delta-V = sqrt(dVx^2 + dVy^2). (Equations 2 to 7 read on page renders; an independent implementation of Eqs. 2 to 4 and 7 is in sec. 3.)
- Computation (pp.19-23): a MAD program for the IBM 7094 maps delta-V over the plane (Fig. 4). The common peripoint (periapsis of both hyperbolas) gives delta-V = difference of the two periapsis speeds, eqs. (8) and (9)
  (read on the AIAA copy's render; the report's printed layout is the same), which the text says is under 3 percent above the minimum for a typical case, so the common-peripoint value serves as the optimum when its
  radius exceeds the planet radius. If the common peripoint is below the surface, the constrained optimum is approximated along the less energetic boundary hyperbola (the one skimming the surface). The argument
  given (pp.22-23): contours of constant delta-V are nested closed curves, so the constrained optimum lies on the boundary of the permissible region, on the hyperbola with the smaller V_H because a lower speed gives a larger turn.
- Results (pp.24-26, Figs. 5 to 9): launch velocity contours of 0.2 EMOS for direct and via-Venus trips, launch years 1970 to 1990. Venus both ways in 1971, outbound only 1973, homebound only 1975, both again 1978.
  Fig. 6 (arrival at Mars JD 2442860, 1976): a pure flyby would go below the surface and a neighbouring direct trip is cheaper. Fig. 7 (JD 2441120, 1971): via-Venus clearly beats direct; pure flyby possible.
  Fig. 8: varying the Venus date with Earth and Mars dates fixed; the pure-flyby date lies near the overall optimum because the delta-V at Venus changes faster with date than that at Earth or Mars; the thrust saving is
  "a few hundred feet per second". Fig. 9 (1984 opposition): Venus both ways allows a Mars stay over a hundred days, or two vehicles (one via Venus, one direct) to meet at Mars and return via Venus.
- Conclusions (pp.27-28), six: (1) the method yields optimum via-Venus transfers for 1970 to 1990; (2) the common-peripoint solution is a few percent above optimum; (3) when a pure flyby would pass below the surface,
  neighbouring direct trips are cheaper; (4) the saving from thrust at Venus is a few hundred ft/s; (5) the figures verify the Fig. 1 availability rule; (6) the results apply only to the trip types studied (other authors,
  refs 12 to 14, found large savings from thrust at Mars flybys).
- References 1 to 15 (p.40): the same as the AIAA copy (Crocco 1956; Lockheed 1962; Ford Aeronutronic EMPIRE 1962; Battin 1964; Hollister 1963; Sohn 1964; Deerwester 1965; Hollister 1964; Gobetz 1963;
  Ross 1963; Leitmann 1962; Lee 1964; Luidens and Kappraff 1965; Titus 1965; Fimple 1962).

## 2. Figures (one line each)

1 alignment ring (viewed). 2 the planetary encounter and 3 the flyby parameters (geometry sketches, not read in detail). 4 loci of constant delta-V (viewed; numbers in the data file: A = 60 deg, V_HI = 0.741, V_HO = 1.07,
contours .19 to .25, common peripoint 0.194, true optimum 0.189). 5 launch-velocity contour fold-out (viewed on four sheets; via-Venus minima .123, .142, .131, .161, .197, .109, .149, .107, .143, .132, .170, .129, .124, .176, .123
read once each). 6 and 7 total delta-V against trip time (viewed; approximate readings in the data file). 8 delta-V against date at Venus (low resolution). 9 1984 opposition contours (low resolution).

## 3. Independent check of Fig. 4 (ours)

A two-body implementation of Eqs. (2) to (4) and (7), with the Fig. 4 inputs (V_HI = 0.741, V_HO = 1.07, A = 60 deg), minimised over the plane outside the planet: the true optimum delta-V = 0.1889 at 1.194 planet radii
(printed 0.189, agrees to the printed digits). The common peripoint (shared periapsis, e = 1 + r V_H^2, asin(1/e_I) + asin(1/e_O) = A) is at r = 1.276 with delta-V = 0.1921, against the printed 0.194 (1 percent higher). The
gap of 0.002 is not explained: the report prints no peripoint radius; the figure's common-peripoint marker sits near the surface at about 1.3 radii, consistent with our 1.276. The optimum reproduces, so the formulation (Eqs. 2 to 7 and the unit system) is confirmed; the common-peripoint
number would be the one to quote with a caution. Energy check of the injection formula (2): speed squared equals V_H^2 + 2/r to machine precision.

## 4. Relevance to this project

- **Venus-swingby Mars round trips versus Earth-Venus cyclers.** The paper is the one-shot Earth-Venus-Mars stopover case: Earth launch, one Venus flyby on the outbound or homebound leg, Mars stopover, no repetition.
  It contains no periodic orbit, no periodic-orbit statement and no cycler statement. The Hollister periodic-orbit lineage (R1(b): Hollister 1967 to 1969, Hollister-Menning 1970 Earth-Venus periodic swingby orbits)
  is separate and later; this report is the stopover-mission predecessor from the same group (Hollister 1963 thesis, ref 5, is its source). It is background to R1(c) (Venus-Mars geometry), not a precursor in the
  periodic sense.
- **Collision with the R1(c) gate.** None new: no planet constants, ephemeris or element set are stated (data file `planet_constants`); the held 65-700 digest's conclusion stands.
- **Useful content.** (a) The availability rule (Sun-Venus-Mars alignment about half way, alignment relative to Earth) restates the geometry that a V-M leg of an Earth-Venus-Mars construction must satisfy; the
  predicted availability years (1971 both ways, 1973 outbound, 1975 homebound, 1978 both) can be a cross-check on any Venus-Mars leg generator. (b) The Fig. 4 optimum reproduces independently, so Fig. 4 is
  a usable sourced control for a planar flyby-impulse optimiser (inputs and optimum 0.189 printed on the figure; the common-peripoint 0.194 is the number we do not reproduce). (c) Fig. 5 via-Venus minimum
  launch velocities (0.107 to 0.197 EMOS, about 3.2 to 5.9 km/s at EMOS = 29.78 km/s, COMPUTED) are sourced single-witness values; the Jan 1978 and Dec 1975 minima (.107, .109) are the lowest.
- **Gate label.** There is no gate anchor for this paper in literature_check.py (grep 2026-10-11). If one is added, the content supports {mga-tour} or {transfer} (a one-shot Venus-Mars gravity-assist
  round trip with an optional impulse), single flyby body Venus, and no cycler, repeated-moon or resonant label.

## 5. What the report does not contain

No tables. No launch-year table (only contours). No periodic-orbit statement. No derivation of eq. (1) beyond the statement of the peripoint-angular-momentum relation. No planet constants. No numeric delta-V listing
for any specific trip other than the Fig. 4 example.

## 6. Cautions

- Author titles: the title page footnotes were read from a spaced-letter text layer ("Deputy Director", "Assistant Professor", "Research Assistant"); the pairing with the two authors was not confirmed on an image.
- The 1965 report was unnumbered "RE-16" in the title; the NTRS record title is the same. The journal form (Astronautica Acta 12(2):169-179 per the held-form digest, 1966) has not been read; any difference from this
  text is unverified (the conference form is identical).
- Fig. 5 via-Venus labels have a single witness.

## 7. Citation mining

All 15 references are the same as the held AIAA 65-700 digest's list (already mined there). New against that digest: none. Items still not held and already noted there: Hollister 1964 AIAA 64-647 (the JSR
engineering note of 1966 is held, the preprint is not), Deerwester AIAA 65-89, Lee AIAA/ION 1964, Luidens and Kappraff NASA TN D-2605, Titus AIAA 65-515 (wanted row 54), Fimple 1962, Gobetz 1963, Ross 1963
(wanted list). No new wanted-list rows proposed.
