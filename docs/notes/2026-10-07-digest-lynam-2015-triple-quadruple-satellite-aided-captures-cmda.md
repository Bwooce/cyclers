# Digest: Lynam 2015, "Broad-search algorithms for finding triple- and quadruple-satellite-aided captures at Jupiter from 2020 to 2080" (CMDA) (#960, #943)

A. E. Lynam (West Virginia University), "Broad-search algorithms for finding triple- and quadruple-satellite-aided
captures at Jupiter from 2020 to 2080", Celestial Mechanics and Dynamical Astronomy 121(4):347-363, April 2015
(online 4 February 2015). doi 10.1007/s10569-015-9602-y (printed on p.1; Crossref title, author, volume, issue and pages
agree). Received 11 November 2014, revised 2 January 2015, accepted 8 January 2015.
- File given: `1d105f69-lynam2015_1.pdf`, 17 pages, text layer (online-first PDF, no printed page numbers).
  md5 d3664f1a90ecca28b661f3e5ec64a354. I cite PDF pages; journal page = 346 + PDF page.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/lynam-2015-broad-search-triple-quadruple-satellite-aided-capture-jupiter-2020-2080-cmda-121-347-doi-10.1007-s10569-015-9602-y.pdf`
- **Wanted list:** not listed.
- **Supplement (not held):** `TripleQuadImportantDates.mat` (online supplementary material; p.14). It holds every first,
  second and third flyby date, perijove date and asymptote date for each sequence type, 2020-2080. Worth fetching for X1.
- **How I read it.**
  - I read the whole text layer. I read pp.7 (Table 2, Fig. 3), 8 (Table 3, Fig. 4, the 7.05-day statement),
    12 (Table 4, Fig. 5), 13 (Tables 5-6, Fig. 6), 14 (Fig. 7, supplement) and 15 (Table 7) on 130-dpi renders.
    Every number quoted from those pages was read on the image.
  - Arithmetic checks: `checks_lynam.py` secs. 1 and 4, output in `checks_lynam.out`.
  - The #943 collision check is in `collision.md` (written first). The verdict is **no collision**.

## 0. Verdict

**A fast, ephemeris-light census of when three or four Galilean moons line up for one capture pass.**
- The method reduces moon motion to two linear phase angles, Delta-lambda(Ca,Ga) and Delta-lambda(Ga,Eu). Io follows from
  the Laplace relation. Each sequence type is a contour in that plane. The "dynamics line" crosses the contours at the
  flyby epochs. The whole 60-year search runs in 50 s.
- **#943 collision: no collision** with gc-1, gc-2, ge-1, ge-2 or ge-3. Callisto-Ganymede-Europa sequences are counted
  (Table 5), so Ganymede-to-Europa legs appear, but only as one inbound pass on a capture hyperbola (rp 6.4-10 R_J,
  V_inf 2.0-6.0 km/s, Table 2). No moon v_inf is printed. My estimate on such conics is 9.2-13.0 km/s at Ganymede
  (INFERRED, `checks_lynam.out` sec. 5); the candidates have 1.4-3.9 km/s. There is no repeat.
- **Related (phasing only):** the Laplacian G-E-I triple "occurs once every 7.05 days" (p.8). That is the Europa-Ganymede
  synodic period (7.051 d, my check). ge-2 (k = 3) spans three of them.
- **What it gives the project:** alignment rates and patterns for X1 (sec. 2), the linear phase-angle machinery (sec. 1),
  and the supplement's date list.
- **Catalogue implication (PROPOSAL only):** none.
- **Lynam 2012 PhD (wanted row 14; the brief says 15, the current list at 52d51c6f has 14):**
  - Post-thesis work (received November 2014, from WVU; the thesis is dated May 2012).
  - It generalises thesis secs. 2.7-2.8 (triple and quadruple capture, the Laplace phase-angle analysis 2.7.1-2.7.2 and
    the Callistan near-resonance 2.7.4) and carries out future-work item 6.1 (preview TOC p.vi).
  - The taxonomy and the Laplace phase-angle conversions are "repeated here for completeness" from Lynam, Kloster &
    Longuski 2011 (p.2-3), which is thesis Ch. 2 (thesis sec. 1.2, refs [10], [11]).
  - It does not cover Ch. 4 (navigation) or Ch. 5 (cyclers). It is not among thesis refs [10]-[16].

## 1. Method (READ, pp.2-9)

- **Taxonomy (Table 1, p.4).** Letters C, G, E, I for flybys; P = perijove without JOI, J = perijove with JOI. Letters
  before P/J are inbound flybys, after are outbound. Four triple categories (C-G-I, C-G-E, C-E-I, Laplacian G-E-I) and
  quadruples. Only 8 of the 16 Laplacian triple permutations and 16 of the 32 quadruple permutations exist; the 1:2:4
  resonance forbids the rest (Lynam et al. 2011). C-E-I is dropped as impractical.
- **Laplacian perijoves:** the resonance fixes perijove at 2.1 R_J for half the types and about 1.1 R_J for the other half (p.3).
- **Phase angles.** Delta-lambda(Ca,Ga) = lambda_Ca - lambda_Ga (eq. 1; the opposite sign to Lynam et al. 2011).
  Ephemeris form via a signed arccos (eqs. 4-6). Laplace relation 180 deg = 2 lambda_Ga - 3 lambda_Eu + lambda_Io (eq. 7)
  and its conversions (eqs. 8-13), for example Delta-lambda(Ga,Io) = 3 Delta-lambda(Ga,Eu) + 180 deg (eq. 10).
  So the C-G-I contours repeat every 120 deg in Delta-lambda(Ga,Eu) (p.8).
- **Input ranges (Table 2, p.7), image-read:**

| category | rp, R_J | V_inf-, km/s (low-thrust / chemical) | flyby altitudes h_p1..h_p4, km |
|---|---|---|---|
| C-G-I | 3.6-6.4 | 2.0-5.5 / 5.5-6.0 | 100, 1000, 300, - |
| C-G-E | 6.4-10 | 2.0-5.5 / 5.5-6.0 | 100, 1000, 300, - |
| G-E-I | 2.1 | 2.0-5.5 / 5.5-6.0 | 100, 1000, 300, - |
| quadruple | 2.1 | 2.0-5.5 / 5.5-6.0 | 100, 1000, 1000, 300 |

- **Contours.** C-G-I and C-G-E types give quasi-rectangular 2-D contours (Fig. 3). Laplacian triples give 1-D ranges of
  Delta-lambda(Ga,Eu) (Table 3, p.8: GIPE 200.17-201.22 deg, GIJE 199.94-200.15, EPIG 208.06-208.89, EJIG 208.08-208.26).
  Quadruples give short 1-D segments (Fig. 4). Contours are discretised as line segments (eqs. 18-19).
- **Dynamics line (eqs. 20-24).** Both phase angles drift linearly at (n_Ca - n_Ga) and (n_Ga - n_Eu), so their ratio
  E1 = (n_Ca - n_Ga)/(n_Ga - n_Eu) gives a straight line. Intersections with the 360-deg-shifted contour segments are
  closed-form (eqs. 25-30). Three samples (25, 50, 75 %) are kept per contour crossing; only crossings are counted.
- **Sun filter (eq. 31).** Sun-asymptote phase angle within 90 +/- 10 deg (chemical) or 90 +/- 25 deg (low thrust).
- **Initial states.** A 3-D patched-conic back-propagation from the first flyby (Rodrigues rotation, eq. 32) gives a state
  for GMAT or Mystic. The author says these initial states "have not yet been tested in a high-fidelity model" (p.14).

## 2. Results (READ; image-checked) - moon-alignment data

- **Table 4 (p.12), C-G-I, 2020-2080, low-thrust / chemical:** PIGC/JIGC 13/6, IPGC/IJGC 99/30, GPIC/GJIC 196/68,
  GIPC/GIJC 147/48, CPIG/CJIG 150/52, CIPG/CIJG 201/66, CGPI/CGJI 104/33, CGIP/CGIJ 17/6. Totals 927/309; 15.5/5.2 per year.
- **Table 5 (p.13), C-G-E:** PEGC/JEGC 4/2, EPGC/EJGC 42/16, GPEC/GJEC 100/36, GEPC/GEJC 52/16, CPEG/CJEG 59/17,
  CEPG/CEJG 109/35, CGPE/CGJE 58/22, CGEP/CGEJ 8/3. Totals 432/147; 7.2/2.5 per year.
- **Table 6 (p.13), Laplacian and quadruple:** CGIPE/CGIJE 0/0, EPIGC/EJIGC 0/0, CEPIG/CEJIG 4/0, GIPEC/GIJEC 2/0,
  EPIG/EJIG 459/174, GIPE/GIJE 470/175. Totals 935/349; 15.6/5.8 per year.
- **Patterns:** CGPI comes "in patterns of 5 or 6 in a two-year period and then ... an approximately two-year gap" (p.12,
  Fig. 5). Quadruples, about 1 per 10 years, "seem to occur in clusters of two separated by about 2.25 years" (p.12).
  Fig. 7 (p.14, figure reading) puts CEPIG near 2026, 2028, 2054, 2056 and GIPEC near 2059, 2061.
- **Independent confirmation (p.14):** the method recovers the 2026 CGPI optimised in MALTO by Patrick & Lynam 2014 and
  the CIJG confirmed in GMAT by Didion & Lynam 2014.
- **Table 7 (p.15), 200-day capture plus perijove raise (PJR), image-read:**

| sequence | 1st PJ | JOI, m/s | JOI + PJR to 9 R_J | JOI + PJR to 14 R_J |
|---|---|---|---|---|
| CGIJ | 5 R_J | 208 | 413 | 563 |
| CGEJ | 9 R_J | 360 | 416 | 569 |
| CGJ | 13 R_J | 568 | - | 661 |
| GIJ | 5 R_J | 349 | 565 | 716 |
| GJ | 13 R_J | 864 | - | 966 |
| IJ | 5 R_J | 555 | 780 | 930 |

- Radiation scale used (p.15, after Garrett et al. 2012): worst below 5 R_J; 1/2 at Io (5.9 R_J), 1/7 at Europa (9.4),
  1/60 at Ganymede (15), negligible at Callisto (26.3).
- Navigation (p.14): an uncorrected error in the first flyby of a CGIP grows "more than 38 times" by the third flyby.

## 3. Checks (`checks_lynam.out`)

- All column sums of Tables 4-6 agree with the printed totals. All per-year rates agree after rounding
  (927/60 = 15.45 printed 15.5; 309/60 = 5.15 printed 5.2; 147/60 = 2.45 printed 2.5).
- All triples: 1236 + 579 + 1278 = 3093; the abstract says "approximately 3,100". Quadruples 4 + 2 = 6 (abstract 6).
- Callisto triples per year: (927 + 432)/60 = 22.65 low thrust and (309 + 147)/60 = 7.60 chemical; the conclusions say
  23 and 8.
- "Once every 7.05 days" = Europa-Ganymede synodic period, 7.051 d.
- Runtime scaling: 167 s x 60/16 x 11 = 6888.8 s (printed 6,889 s); 6889/50 = 137.8 (printed "approximately 138x").
- Table 7 is identical, cell for cell, to Table 1 of Lynam 2015/2016 (doi 10.1007/s10569-015-9649-9), as that paper says.
- **Internal inconsistency (not resolved):** the introduction says "22 in this paper vs. only 2" sequence types (p.2), and
  "11 times the number of flyby types" (p.13) fits 22/2. But the same page says "The entire 60-year search for all 26 types
  took 50 s". The table rows count 8 + 8 + 6 = 22 type pairs. I cannot tell what the 26 counts.

## 4. Citation mining (26 references, pp.16-17)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| work | status |
|---|---|
| Lynam 2014a, Part I, Acta Astronautica 94:246-252 (doi 10.1016/j.actaastro.2013.07.018) | not held; not on the wanted list. **New candidate (low priority);** this paper supersedes its method. |
| Lynam 2014b, Part II, Acta Astronautica 94:253-261 | in batch 35 (this agent; `digest-lynam-2014-acta-part2.md`) |
| Lynam, Kloster & Longuski 2011, CMDA 109:59-84 | manuscript copy in batch 35 with a sibling agent (`5d107f14-MSAC_Laplace_Resonance.pdf`) |
| Lynam & Longuski 2011, JGCD 34:1485-1494 | in batch 35 with a sibling agent (`05b40581-lynam2011_1.pdf`) |
| Didion & Lynam 2014, AIAA 2014-4106 (Earth to Callisto-Io-Ganymede triple capture) | in batch 35 with a sibling agent (`26289618-ImpulsiveTrajectories2014.pdf`) |
| Lynam & Longuski 2012, Acta Astronautica 70:33-43 (navigation) | not held; not on the wanted list. Low priority (navigation). |
| Patrick & Lynam 2014, AIAA 2014-4218 (SEP to triple capture, MALTO) | not held. Low priority. |
| Izzo, Simoes, Maertens, de Croon, Heritier & Yam 2013, "Search for a grand tour of the Jupiter Galilean moons", GECCO 2013 | not held (no izzo hit); not on the wanted list. A multi-moon tour search; possible X1 background. Low priority. |
| Schadegg, Russell & Lantoine 2014, AAS 14-448 (tether capture) | not held. Not relevant. |
| Buffington 2014, AIAA 2014-4105 | HELD (`buffington-2014-trajectory-design-europa-clipper-mission-concept-aiaa-2014-4105-doi-10.2514-6.2014-4105.pdf`) |
| Strange et al. 2012, AIAA 2012-4518; Landau, Strange & Lam 2010 | not held. Low priority. |
| Longman 1968; Cline 1979; Nock & Uphoff 1979; Johannesen & D'Amario 1999; Yam 2008 PhD | not held. Capture history. |
| Davis, Patterson & Howell 2007, AAS 07-275 (solar perturbations, Cassini) | not held (the held Davis papers are other works). |
| Garrett et al. 2012 (JPL 12-9); Riedel et al. 2006 (AutoNav); Sims et al. 2006 (MALTO); Whiffen & Sims 2001 (Mystic); Hughes 2008 (GMAT); Vallado 2007; Laplace 1809; Wilson et al. 1997 | none needed |

- New candidates: Part I (low). Optional: Izzo et al. 2013 for X1 background. Nothing for #943.

*Wanted-list row numbers in this digest are those of the list at commit 52d51c6f.*

*Filed as `cyclers_pdf/papers/lynam-2015-broad-search-triple-quadruple-satellite-aided-captures-jupiter-2020-2080-cmda-121-347-doi-10.1007-s10569-015-9602-y.pdf`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
