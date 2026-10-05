# Digest: Russell & Strange 2007, "Planetary Moon Cycler Trajectories", AAS 07-118 (#960)

R. P. Russell and N. J. Strange, "Planetary Moon Cycler Trajectories", AAS 07-118, AAS/AIAA Space Flight
Mechanics Meeting, Sedona, Arizona, 29 Jan - 1 Feb 2007 (both authors at JPL). No DOI (AAS conference
paper). Filed as `cyclers_pdf/papers/russell-strange-2007-planetary-moon-cycler-trajectories-AAS-07-118.pdf`
(20 pages, text layer, md5 95fef963421bad7f26a0a8d6311fffb6).

It is ref. [26] of Russell & Strange 2009 (JGCD 32(1):143-157, doi 10.2514/1.36610; digest
`docs/notes/2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md`). The 2009 paper is the
journal version without the heliocentric results.

All 20 pages were read from the text layer. Tables 3, 5 (pp.9, 12) and 8 (p.18) were read as page images.

Evidence tags: READ = printed (PDF page cited). COMPUTED = my own check, 2026-10-05. INFERRED = my reading.

Why it was acquired: the `#938` review gate for `#942` (R1) and `#943` (X1). R-S 2009 p.5 says that
Earth-Mars, Earth-Venus, Venus-Mars and Venus-Mercury "heliocentric calibration" searches were run and
that "a few of the resulting previously undocumented heliocentric cyclers are recorded in [26]".

## 0. Verdict in one paragraph

This is the conference parent of R-S 2009. It runs the generalised Russell-Ocampo free-return cycler
search (one working body, massless target) for:
- three heliocentric pairs: Earth-Mars, Venus-Mercury and Venus-Mars.
- four Jovian pairs (five directions).
- Titan-Enceladus.

**There is no Earth-Venus set.** Table 1 lists none, and no Earth-Venus cycler is printed anywhere. So the
R-S 2009 statement that [26] holds Earth-Venus results is not borne out.

The paper does publish:
- one Venus-Mars cycler, VenMar#45. It stays ballistic in a patched-conic ephemeris model.
- three Venus-Mercury cyclers, ideal model only. The ephemeris was abandoned because of Mercury's
  eccentricity.

The Jovian and Saturnian rows are the same rows as in R-S 2009 (the strings match exactly), plus EurGan#93,
which R-S 2009 drops. So 2007 is the first publication of all 30 catalogued `russell-strange-2009-*`
member rows.

Every ideal-model cycler here has a massless target. The authors name "removing the massless assumption
of the target body" as future work, "the Jovian system cyclers in particular" (p.18).

## 1. Method (READ pp.1-6)

- Two steps (p.2):
  1. an exactly periodic ideal model (circular, coplanar).
  2. evaluation and optimisation over multiple cycles in (a) a zero-radius-SOI patched-conic ephemeris
     model and (b) an integrated n-body plus oblateness model (Mystic).
- "Free-return cyclers ... where one of the two orbiting celestial bodies in the ideal model is
  considered massless" (p.2). The flyby body hosts every return. The target is met once per cycle on one
  leg.
- The v_inf globe (Fig. 1, p.3) shows three types of return:
  - full-rev circles (even n pi, resonant).
  - half-rev x's and o's (odd n pi).
  - generic dots (non-n pi, in plane).
- The cycler search enumerates 1..n-leg patterns whose total TOF is an integer number of synodic periods,
  as in Russell-Ocampo 2005 (ref. [7]) (p.3).
- New non-resonant return solver (Fig. 5 algorithm, Eqs. 1-7, pp.5-6):
  - a 1-D root solve on the transfer period T for each (M, N, inbound/outbound).
  - Units: flyby-body orbit radius = 1 LU and primary mu = 1.
  - Example: 220 direct non-resonant returns at v_inf = 0.5 LU/TU for M <= 9 (p.6).
- Retrograde returns are left out (p.6). Eqs. (8)-(9) give the modification.
- The homotopy to the ephemeris linearly interpolates body positions from the ideal model to the
  ephemeris. It replaces the mean-element homotopy of ref. [9] (p.6).
- Rule of thumb (p.11): periods up to 8 flyby-body periods, about 10 legs at most. An exhaustive search
  takes "less than a few days on one modern processor".
- Odd n pi returns were sought in every search, but "ballistic solutions that include odd-n pi free
  returns are much less common for the planet-centered cyclers" (footnote, p.14).

## 2. Models and constants (READ p.7, Tables 1-2)

Table 1, ideal models (primary: flyby body -> target):
- Sun: Earth -> Mars, Venus -> Mercury, Venus -> Mars.
- Jupiter: Ganymede -> Io, Ganymede -> Europa, Ganymede -> Callisto, Europa -> Ganymede.
- Saturn: Titan -> Enceladus.

Earth-Mars "is used to test and calibrate the improved methods noting that solutions to this system are
well documented in [7,9]" (p.7). No Earth-Mars results are printed.

Table 2 (p.8), the ideal-model circular period in seconds, which is needed to reproduce the nomenclature
strings:

| Body | Period (s) |
|---|---|
| Mercury | 7,600,552 |
| Venus | 19,414,153 |
| Mars | 59,354,429 |
| Io | 152,854 |
| Europa | 306,822 |
| Ganymede | 618,153 |
| Callisto | 1,441,931 |
| Titan | 1,377,684 |
| Enceladus | 118,387 |

Gravitational parameters (km^3/s^2) and radii (km):

| Body | mu (km^3/s^2) | Radius (km) |
|---|---|---|
| Sun | 1.3271244e11 | 696,000 |
| Jupiter | 126,686,535 | 71,492 |
| Saturn | 37,931,208 | 60,268 |
| Mercury | 22,321 | 2,440 |
| Venus | 324,860 | 6,052 |
| Mars | 42,828.3 | 3,399 |
| Io | 5,959.92 | 1,827 |
| Europa | 3,202.74 | 1,561 |
| Ganymede | 9,887.83 | 2,634 |
| Callisto | 7,179.29 | 2,408 |
| Titan | 8,978.14 | 2,575 |
| Enceladus | 6.95 | 256.3 |

- Earth's row is absent from Table 2 (INFERRED: Earth-Mars was calibration only).
- Mars at 59,354,429 s is 1.8808 yr (COMPUTED). The project's `RussellModel` uses 1.875 yr, so an exact
  string reproduction needs the Table 2 values.
- Ephemerides: DE414, jup230 and sat242 (p.8 footnote). Flyby minimum altitude: 1000 km at Titan (p.4).

## 3. Results (READ pp.9-14, Tables 3-6)

Table 3 (p.9, image-read) has the heliocentric and Jovian rows. Body A is the flyby body and Body B is the
massless target.

| ID | Synodic period (d) | v_inf A / B (km/s) | Legs | Period (d) | Petal (yr) | Min flyby alt. at A (km) |
|---|---|---|---|---|---|---|
| VenMar#45 | 333.9 | 8.22 / 12.96 | 1 | 667.8 | -65.68 | 19,784 |
| VenMer#22 | 144.6 | 6.62 / 8.61 | 1 | 433.7 | -16.99 | 3,322 |
| VenMer#69 | 144.6 | 8.03 / 10.68 | 2 | 722.8 | 9.13 | 6,111 |
| VenMer#75 | 144.6 | 10.59 / 14.57 | 3 | 1,156.5 | 21.54 | 5,108 |
| EurGan#93 | 7.05 | 2.37 / 4.10 | 3 | 28.2 | -1.33 | 1,293 |
| EurGan#131 | 7.05 | 2.40 / 4.10 | 2 | 21.2 | -1.33 | 1,113 |
| EurGan#159 | 7.05 | 2.45 / 4.11 | 3 | 28.2 | -1.33 | 1,261 |
| GanCal#1 | 12.52 | 3.18 / 3.26 | 3 | 37.6 | 0.41 | 247 |
| GanCal#5 | 12.52 | 3.24 / 3.34 | 2 | 37.6 | 0.41 | 328 |
| GanEur#5 | 7.05 | 1.66 / 2.57 | 1 | 35.3 | -1.33 | 1,819 |
| GanEur#43 | 7.05 | 1.87 / 3.89 | 1 | 14.1 | -1.33 | 8,861 |
| GanEur#316 | 7.05 | 3.20 / 3.81 | 4 | 49.4 | -1.33 | 1,447 |
| GanIo#53 | 2.35 | 3.90 / 9.85 | 2 | 21.2 | -1.33 | 518 |
| GanIo#185 | 2.35 | 3.97 / 9.90 | 6 | 49.4 | -1.33 | 603 |
| GanIo#403 | 2.35 | 4.29 / 4.34 | 2 | 56.4 | -1.33 | 540 |

Table 3 also gives the minimum and maximum distance to the primary (km) and the A->B and B->A transit
times (d). For example, VenMar#45 has 108,067,501 / 341,571,371 km and transits of 554.89 / 112.96 d.

Table 5 (p.12, image-read) gives the formal nomenclature (McConaghy-Russell-Longuski 2005). Capital
letters mark the target-encounter leg.

- VenMar#45: G(2.97216,349.97729,U)
- VenMer#22: G(1.93012,1054.84284,U)
- VenMer#69: g(1.33364,480.11127,Ls) G(1.88322,1037.96013,U)
- VenMer#75: g(1.31055,471.79942,U) G(1.83643,1021.11483,U) f(2:3,74.20508,0.00960)
- EurGan#93: g(1.97030,349.30891,U) G(3.97176,709.83427,U) f(2:1,88.69348,0.57296)
- EurGan#131: G(3.95655,704.35739,U) f(2:1,87.95239,90.00000)
- EurGan#159: G(3.94206,699.14318,U) f(2:1,87.24509,120.26932) f(2:1,87.24509,59.73068)
- GanCal#1: G(1.74871,269.53421,U) g(1.50246,540.88534,L) f(2:1,77.40130,0.03291)
- GanCal#5: g(1.50425,541.53130,L) G(3.74691,628.88825,U)
- GanEur#5: G(4.92758,2493.92898,U)
- GanEur#43: G(1.97103,1069.57159,U)
- GanEur#316: g(1.31322,472.76044,U) h(1.5,540.0,L,-3.98557) g(1.31322,472.76044,U) G(2.77217,1357.97970,U)
- GanIo#53: g(0.97232,710.03419,U) G(1.98423,1434.32320,U)
- GanIo#185: g(0.96288,706.63724,U) f(1:2,84.97911,89.99999) g(0.96288,706.63724,U) f(1:2,84.97911,57.73932) G(1.97285,1430.22609,U) f(1:2,84.97911,90.00890)
- GanIo#403: g(3.72334,1700.40403,Ll) G(4.16078,2217.88234,L)

Tables 4 and 6 (pp.10, 13) give the 20 Titan-Enceladus rows: IDs 37, 145, 183, 207, 217, 227, 231, 235,
314, 370, 492, 510, 539, 552, 572, 586, 594, 602, 624 and 631.
- The synodic period is 1.50 d.
- Example: TitEnc#235 is g(0.88468,678.48383,U) g(1.22599,441.35506,U) f(2:3,55.18988,179.99996)
  F(1:2,57.76202,180.0) f(1:2,57.76202,154.65857), v_inf T/E 3.18/6.04 km/s, period 97.4 d.

Cross-check against R-S 2009 Table 5 (READ, held PDF): the EurGan#131, #159, GanCal#1, #5 and GanEur#5,
#43, #316 strings are identical. R-S 2009 Table 5 omits EurGan#93 and all Ven* rows.

Ephemeris and high-fidelity results (READ pp.14-18):
- **VenMar#45** "easily converges to ballistic" in the patched-conic ephemeris model and is "similar to the
  Earth-Mars Aldrin cycler" because it is a single multi-revolution non-resonant return (Fig. 9a, p.14).
- **Venus-Mercury:** "no ephemeris solutions are presented because the large eccentricity of Mercury
  proved too high for multiple cycles to remain ballistic in the ephemeris model" (p.14). The text names
  "#22 and #45" as Venus-Mercury cyclers. #45 is VenMar, so this is probably a slip for #69 or #75
  (INFERRED).
- GanCal#1 is ballistic over 10 cycles (Fig. 9b). GanEur#316 and EurGan#131 are ballistic (Fig. 10).
  TitEnc#37, #235, #572 and #586 appear in Figs. 11-13.
- Table 7 (p.17) gives the patched-conic ephemeris TitEnc#235. Legs 1-27 are shown, and the start is
  8774.549 d after J2000 (10 Jan 2024).
- Table 8 (p.18) gives single cycles optimised in a high-fidelity n-body plus oblateness model:
  - TitEnc#183: 32 m/s.
  - TitEnc#235: 63 m/s. A second cycle needed 40 m/s.
  - TitEnc#586: 11 m/s.
  - GanEur#316: 121 m/s. Five Ganymede encounters and a targeted Europa flyby at 495 km.
  - EurGan#131 (not tabulated): 58 m/s.

## 4. Gate answers for `#942` (R1) and `#943` (X1)

- **Heliocentric "calibration" sets that were run:** Earth-Mars (calibration, no rows printed),
  Venus-Mercury (3 rows) and Venus-Mars (1 row). **Earth-Venus: not run in this paper and not printed.**
  The R-S 2009 p.5 claim that it is recorded in [26] is not borne out. This is a citation slip in the 2009
  paper, offered with respect.
- **One-working-body Earth-Venus catalogue: no.** There is no Earth-Venus content.
- **One-working-body Venus-Mars catalogue: one documented member, VenMar#45.** The paper says the "complete
  trajectories are archived for future use" (p.9). That is not a publication.
- **Jovian pairs:**
  - Ganymede-Callisto: GanCal#1 and #5. Ganymede is the flyby body and Callisto is massless.
  - Ganymede-Europa in both directions: GanEur#5, #43 and #316 (Ganymede flyby, Europa massless), and
    EurGan#93, #131 and #159 (Europa flyby, Ganymede massless).
  - All ideal-model targets are massless (p.2, p.7).
  - Table 8(d) GanEur#316 is an n-body single-cycle optimisation with both moons physically present. It
    has a 495 km Europa flyby and costs 121 m/s, so it is powered, not ballistic. It is not an ideal-model
    two-working-body cycler.
- **Does anything make X1 a reproduction? No.** p.18: "improved ideal model (circular-coplanar) cyclers
  could be sought by removing the massless assumption of the target body. The Jovian system cyclers in
  particular would benefit from such a change." This is author-named future work, which supports novelty
  policy (i) or (ii).

Per-cell verdicts (combining this paper with the Pisarevsky 2008 digest):

| Cell | Verdict | What is published | What is left |
|---|---|---|---|
| R1(a) Earth-Mars, Mars turns | PARTIAL | Pisarevsky 2008: method, Table 4 (one spatial class III member), Fig. 14 (graphical class I.1 candidates) | enumeration with turn gating, numeric coplanar class I and II members, ephemeris |
| R1(b) Earth-Venus, Venus returns | OPEN | nothing in either paper. AAS 07-118 ran no Earth-Venus set. Precedents: Hollister & Menning 1970 (15 rows), Pisarevsky method | the whole cell |
| R1(c) Venus-Mars | PARTIAL | VenMar#45 (one-working-body, Venus hosts, Mars massless, ballistic in patched-conic ephemeris). The `#938` phrase "never revisited" since Rall is superseded | the two-working-body Venus-Mars cell, and the rest of the one-body catalogue (archived, unpublished) |
| X1 Ganymede-Callisto | OPEN | only the Ganymede-hosted one-body limit (GanCal#1, #5) | both moons bending |
| X1 Ganymede-Europa | OPEN | both one-body limits (GanEur Ganymede-hosted, EurGan Europa-hosted), and a powered n-body single cycle (Table 8(d)) | a ballistic both-moons-bending ideal-model cycler |

"OPEN" means open in these two papers plus the held corpus. Post-2007 citers of R-S 2007 and 2009 were not
searched here. `#943`'s own literature gate should do that.

Side effect on `#938` R5 (Venus-Mercury): the ideal-model Venus-Mercury cyclers are published (VenMer#22,
#69, #75). A patched-conic ephemeris attempt is published as a negative (p.14). R5's ideal-model step is a
reproduction. Only the real-ephemeris step, with a better method, stays open.

## 5. Positive controls with sourced numbers

All of these are 5-decimal Table 5 strings plus the Table 3 metrics. Reproduce them with the Table 2
periods.

1. **GanCal#1**: G(1.74871,269.53421,U) g(1.50246,540.88534,L) f(2:1,77.40130,0.03291). v_inf 3.18 / 3.26,
   period 37.6 d, min Ganymede altitude 247 km. This is the X1 one-body limit at Ganymede-Callisto.
2. **GanCal#5**: g(1.50425,541.53130,L) G(3.74691,628.88825,U). 3.24 / 3.34, 37.6 d, 328 km.
3. **EurGan#131** (Europa-hosted) and **GanEur#43** (Ganymede-hosted). These are the two ends of an X1
   Ganymede-Europa mass continuation.
4. **VenMar#45**: G(2.97216,349.97729,U). 8.22 / 12.96 km/s, 667.8 d = 2 x 333.9 d. This is the R1(c)
   one-body control.
5. **VenMer#22**: G(1.93012,1054.84284,U). 433.7 d = 3 x 144.6 d. This is the R5 control.
6. **TitEnc#235** together with Table 7 (ephemeris legs, start 10 Jan 2024). This is the
   patched-conic-ephemeris control.

## 6. Errata candidates (respectful framing)

- Ref. [6], Russell & Ocampo 2005 JSR "Geometric Analysis ...", is printed as 42(1):694-698. Crossref gives
  42(1):138-152, doi 10.2514/1.5571. The 694-698 pages belong to ref. [8], the nomenclature paper.
- p.14 "Venus-Mercury ideal model cyclers such as #22 and #45": #45 is the Venus-Mars ID.
- R-S 2009 p.5 says [26] records Earth-Venus cyclers. This paper does not.

## 7. Citation mining (policy step 4)

The references [1]-[22] (pp.19-20) were checked against `CORPUS_INDEX.md`, the digest bodies and the
filenames on 2026-10-05.

Held:
- [3] Byrnes-Longuski-Aldrin 1993.
- [8] McConaghy-Russell-Longuski 2005 nomenclature.
- [9] Russell & Ocampo 2006.
- [11] Byrnes-McConaghy-Longuski 2002.
- [13] Chen et al. AIAA 2002-4421.
- [5] Russell & Ocampo 2004 JGCD: covered by the AAS 03-145 preprint and the Russell 2004 dissertation.
- [7] the 2005 Global Search: covered by the dissertation (INFERRED).
- [10] AAS 03-509: covered by the McConaghy 2004 PhD and McConaghy-Landau-Yam 2006.
- [1] and [2]: partly covered by Hollister-Rall 1970 NASA CR and Hollister-Menning 1970.

Not held, in priority order. The DOI was checked through the Crossref API on title, authors, volume and
pages. AAS papers have no DOI.

1. **Turner, A. (2007), "Low Road to Mars: The Venus-Mars Cycler", AAS 07-175.** No DOI (AAS paper).
   Unlocks: the direct literal-collision check for R1(c). It is the study that prompted R-S's Venus-Mars
   search (p.7).
2. Hollister, W. M. (1969), "Periodic Orbits for Interplanetary Flight", JSR 6(4):366-369, doi
   10.2514/3.29664 (CONFIRMED).
3. Russell, R. P. & Ocampo, C. A. (2005), "Geometric Analysis of Free-Return Trajectories Following a
   Gravity-Assisted Flyby", JSR 42(1):138-152, doi 10.2514/1.5571 (CONFIRMED). The v_inf-globe geometry
   that the search rests on.
4. Russell, R. P. & Ocampo, C. A. (2005), "Global Search for Idealized Free-Return Earth-Mars Cyclers",
   JGCD 28(2):194-208, doi 10.2514/1.8696 (CONFIRMED).
5. Russell, R. P. & Ocampo, C. A. (2004), "Systematic Method for Constructing Earth-Mars Cyclers Using
   Free-Return Trajectories", JGCD 27(3):321-335, doi 10.2514/1.1011 (CONFIRMED).
6. Strange, N. J. & Sims, J. A. (2001), "Methods for the Design of V-Infinity Leveraging Maneuvers", AAS
   01-437. No DOI (AAS paper).
7. Landau, D. F. & Longuski, J. M. (2006), "Guidance Strategy for Hyperbolic Rendezvous", AIAA 2006-6299,
   doi 10.2514/6.2006-6299 (CONFIRMED).
8. Nock, K. T. et al. (2003), "An Interplanetary Rapid Transit System Between Earth and Mars", STAIF 2003,
   AIP Conf. Proc. 654:1075-1086, doi 10.1063/1.1541404 (CONFIRMED; R-S print pp.1074-1086).
9. Niehoff, J. (1986), "Pathways to Mars: New Trajectory Opportunities", AAS 86-172. No DOI (AAS paper).
10. Prussing, J. E. (2000), "A Class of Optimal Two-Impulse Rendezvous Using Multiple-Revolution Lambert
    Solutions", J. Astronautical Sciences 48(2-3):131-148, doi 10.1007/BF03546273 (CONFIRMED; Crossref
    misspells the author as "Trussing").
11. Shen, H. & Tsiotras, P. (2003), "Using Battin's Method to Obtain Multiple-Revolution Lambert's
    Solutions", AAS 03-568. No DOI (AAS paper). The related journal paper, JGCD 2003, doi 10.2514/2.5014,
    has a different title.
12. Whiffen, G. J. & Sims, J. A. (2002), AAS 02-208 (SDC optimal control). No DOI (AAS paper).
13. Rall, C. S. (1969), PhD, MIT. No DOI (thesis).
14. Seidelmann, P. K. et al. (2002), IAU/IAG WG report 2000, CMDA 82(1):83-111, doi
    10.1023/A:1013939327465 (resolves at doi.org; Crossref record not returned, so title UNCONFIRMED).
15. Prussing & Conway (1993), "Orbital Mechanics", OUP; Wiesel (1997), "Spaceflight Dynamics". Textbooks.
    Low priority.

## 8. Catalogue and code implications (proposals only; no catalogue writes)

- **Priority date.** The 30 `russell-strange-2009-*` member rows (ganio, ganeur, gancal, eurgan, titenc)
  first appeared here (2007), with identical nomenclature strings.
  - Proposal: set `first_published` to AAS 07-118 (2007, no DOI) and keep R-S 2009 as
    `corroborating_sources`.
  - This needs the ratchets that check first_published DOIs to accept a DOI-less AAS source. Check
    `tests/data` before editing.
- **Label error in the existing rows.** The catalogued notes read, for example, "GanCal#1 ... min Callisto
  flyby altitude 247 km". Table 3's column is "Min flyby alt. at Body A", the FLYBY body (Ganymede for
  GanCal, Europa for EurGan). R-S 2009's header says the same. The notes name the target body. This is a
  wording fix in the notes; the values are correct.
- **New rows available (V0, sourced strings):**
  - EurGan#93 (dropped from R-S 2009).
  - VenMar#45, which has a published patched-conic ephemeris ballistic result, so possibly more than V0
    after review.
  - VenMer#22, #69 and #75.
- `literature_check.py`: added KNOWN_CORPUS anchors for Venus-Mars (`russell-strange-2007-venmar`) and
  Venus-Mercury (`russell-strange-2007-venmer`) in this task. The Jovian and Saturnian body pairs stay
  covered by the R-S 2009 anchors.
