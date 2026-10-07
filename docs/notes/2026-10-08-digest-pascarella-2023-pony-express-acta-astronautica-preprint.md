# Digest: Pascarella et al. 2023, "Low-thrust trajectory optimization for the Solar System Pony Express" (Acta Astronautica preprint): diff against AAS 22-015

A. Pascarella, R. Woollands, E. Pellegrini, M. Sanchez Net, H. Xie & J. Vander Hook, "Low-thrust trajectory
optimization for the solar system pony express", *Acta Astronautica* 203 (Feb 2023) 280-290,
**DOI 10.1016/j.actaastro.2022.11.046** (Crossref: title, six authors in this order, vol. 203, pp. 280-290,
issued 2023-02, record created 2022-11-26).

- Source file: upload `9e80088c-CL22_6828.pdf` (JPL clearance CL#22-6828 in the file name only; no stamp on
  the pages), 28 pp, md5 29ad8c37536aed2aea9bd32441439e0f. pdfTeX 1.40.24, created 2022-11-29 AEDT.
  Elsevier preprint class, double-spaced, line numbers. Footer p.1: "Preprint submitted to Acta Astronautica,
  November 28, 2022". So this is the **accepted author preprint**, not the typeset article. Digital text layer;
  no OCR needed. Figures are raster images (no text in them).
- Proposed filename:
  `cyclers_pdf/papers/pascarella-woollands-pellegrini-sanchez-net-xie-vander-hook-2023-low-thrust-trajectory-optimization-solar-system-pony-express-acta-astronautica-203-280-doi-10.1016-j.actaastro.2022.11.046-preprint.pdf`
- Compared against the held conference form `AAS-22-015-pascarella-pony-express.pdf` (17 pp). Crossref gives
  its published form as 44th Annual AAS Guidance, Navigation and Control Conference 2022, *Advances in the
  Astronautical Sciences* (Springer, 2024), pp. 45-61, DOI 10.1007/978-3-031-51928-4_4. OUTSTANDING H.3 calls
  the venue "AAS/AIAA Space Flight Mechanics Meeting"; that is wrong. The held AAS form has no prior digest
  file; OUTSTANDING H.3 is its digest. Companion paper Sanchez Net et al. 2022 JSR 59(3):861-870 is held
  and mined in `docs/notes/s1l1-target-topology-mining.md` §2.9.
- How I read it: full pdftotext of both versions; word-level and number-level diff of the two text layers;
  every page of both rendered at 110 dpi and all figure pages looked at; Table 1 (p.17) read at 300 dpi. Table 1
  has two witnesses (page image and the pdfTeX text layer); all 21 cells agree (`witness-comparison.tsv`).
  Arithmetic in `check_table1.py` / `check_table1.out`.

## 0. Verdict

Same study, same headline solution, same method. The journal is **not content-identical** to AAS 22-015. It adds
one new table and a new section, and it changes several numbers:

- **New Table 1 (p.17)**: three fully converged solutions A, B, C, with departure and four flyby dates each,
  and injection and targeting propellant in kg. The AAS form gives only one solution and round numbers.
  Solution A is the AAS solution (same Figs., same mass curves).
- **Propellant for Solution A is now 33.23 kg (COI) + 1.65 kg (targeting) = 34.88 kg.** The text and the
  conclusion still say 36 kg and "a further 2 kg". The journal is not self-consistent here.
- Maintenance bound tightened from "< 5 kg" (AAS p.15) to "< 2 kg" (p.17).
- "278 of STAR's patched conic trajectories" becomes "278 out of 2225" (12.5 %), and "several" fully converged
  becomes "3".
- The data goal changes from "1 Petabit per year" with ">8000 Tbits" over 8 flybys (AAS) to "1 Petabit per
  flyby", and the crosslink figure now shows 2 flybys in Pbit.
- New Section 5 "Cross-Link Analysis" with DSOC terminal numbers (22 cm, 4 W, > 2.5 Gbps assumed), SOLT, and a
  new distance plot (Fig. 16).

What it gives the project:
- No catalogue numbers. No catalogue row uses Pascarella; the two Pony Express rows cite Sanchez Net 2022 JSR
  and none of their values occur in either Pascarella version (`catalogue-check.md` section 2).
- It gives a better source for the low-thrust maintenance figure used in code: 1.65 kg from 466.77 kg at
  Isp 4155 s = **144.3 m/s** for Solution A over about 6.3 yr. B: 182.3 m/s, C: 147.1 m/s. The project uses
  163 m/s (2 kg from 500 kg), which uses the wrong starting mass.
- It gives three dated real-ephemeris low-thrust cycler solutions. With dates only (no v_inf, no altitudes),
  they are weak anchors, but Solution A's dates fix the positive-control epoch better than the 2032-05-11 epoch
  that `tests/verify/test_pascarella_low_thrust.py` uses with the Sanchez Net cycler-1 row.

PROPOSALS (do not apply without the lead):
1. `dv_band_acceptance.py` lines 40 and 97: re-source the low-thrust figure to this journal, Table 1 A,
   144 m/s over about 6.3 yr (or the range 144-182 m/s over A-C). The 500 m/s ceiling keeps its margin.
2. `real_closure.py` lines 26 and 91: the "100,000-200,000 km drift" claim is not in either version. The
   plotted deviation (patched-conic vs low-thrust) peaks at about 0.017 AU, about 2.5 million km, and is not a
   cycle-to-cycle drift. Remove or re-source that sentence. "Pascarella 2024" should be 2023 (journal) or
   2022 (AAS).
3. `test_pascarella_low_thrust.py`: the test is not a Pascarella reproduction. Re-scope the name, or rebuild it
   on Solution A dates.
4. CORPUS_INDEX line 180: "mined-by-catalogue (Pony Express rows)" is wrong; no row uses it. OUTSTANDING H.3:
   fix the venue and add the journal DOI and Table 1.
5. File this preprint under the proposed name and index it next to AAS 22-015.

## 1. Section map (AAS -> journal)

| AAS 22-015 | journal preprint | change |
|---|---|---|
| Abstract | Abstract + keywords | "We present a high-fidelity candidate solution ... 36 kg ... further 2 kg ... eight subsequent flybys over a period of six years" becomes "We present three high-fidelity candidate solutions ... In all three cases ... feasible". The numbers left the abstract. "courier" becomes "data mule" throughout. "very feasible" becomes "feasible". Keywords added. |
| Introduction | 1. Introduction | Adds "> 1 Gbps have been demonstrated in an optical link between the Moon and Earth, ... at least one order of magnitude larger than with RF" (p.3). Otherwise same. |
| Indirect optimal control | 2. Indirect Optimal Control Formulation, 2.1-2.3 | Same equations. Adds Fig. 2 (LVLH frame diagram) and refs to Pontani 2021 and Biggs & Maclean for the RTN frame. Smoothing refs move from end of paragraph to [6, 7]. |
| Trajectory design, Assumptions | 3, 3.1 | Same numbers (500 kg, NEXT, Isp 4155 s, 0.235 N, altitude <= 25,000 km, >= 300 km Mars, >= 1000 km Earth). The no-thrust-during-flyby rationale is rewritten as a constraint ("greatly simplifies the attitude and pointing requirements"). |
| Methodology, Cycler injection | 3.2, 3.3 | Same: RK 9(8), Matlab fsolve multiple shooting; mid-June 2035 C3 minimum, about 190 days; target Mars flyby at least 800 days after departure; target point at 3 x Mars SOI. |
| Cycler orbiting and targeting, Steps 1-5 | 3.4 | Same five steps. Adds Fig. 5 (STAR flyby vs Keplerian flyby after Step 2). Step 4 ends "present three candidate solutions" instead of the single solution. |
| Results | 4. Results | New Table 1; new numbers (section 0). Adds Fig. 14 (impulsive vs low-thrust deviation). |
| (crosslink paragraph inside Results) | 5. Cross-Link Analysis | New section: Fig. 16 distance plot, DSOC terminal, PPM, CCSDS HPE, SOLT, Fig. 17 redrawn. |
| Future work | 6. Future Work | Last sentence removed: "we also intend to improve the efficiency of the optimization procedure in order to more quickly converge low-thrust trajectories". Opening sentence on smaller spacecraft classes added. |
| Conclusion | 7. Conclusion | Same text (36 kg, 2 kg, eight flybys, six years), "very feasible" becomes "feasible". |
| References (15) | References (17) | Adds Pontani 2021 (JOTA 191) and Biggs & Maclean (IEEE TAES 50). DOIs dropped. AAS ref 15 "Small Satellites on Cycler Orbits with Optical Communication Enable Regular High Volume Data Transfer from Mars, JSR 2022" becomes ref [17] "Cycler orbits and the Solar System Pony Express, JSR" (title of the held JSR 59(3) paper). Hua Xie added as author. |

## 2. Figure map (AAS -> journal)

| AAS | content | journal | change |
|---|---|---|---|
| Fig. 1 | ConOps | Fig. 1 | caption adds [1] |
| - | LVLH frame | **Fig. 2** | new |
| Fig. 2 | C3 contours, Earth-Mars direct | Fig. 3 | same |
| Fig. 3 | STAR dataset, Mars-Earth < 12 months | Fig. 4 | same |
| - | STAR flyby vs Keplerian flyby | **Fig. 5** | new |
| Fig. 4 | Step 3 multiple shooting | Fig. 6 | redrawn, before/after panels |
| Fig. 5 | Step 4 multiple shooting | Fig. 7 | redrawn |
| Fig. 6 | COI trajectory | Fig. 8 | same data, restyled |
| Fig. 7 | COI thrust and mass (500 -> about 467 kg) | Fig. 9 | same data |
| Fig. 8 | STAR solution | Fig. 10 | same |
| Fig. 9 | impulsive ephemeris | Fig. 11 | same |
| Fig. 10 | low-thrust ephemeris | Fig. 12 | same |
| Fig. 11 | deviation patched-conic vs low-thrust, linear signed axis, peak about +0.0173 AU (dy, about 2040.8, read at 200 dpi); all three components return to about 0 near the flybys | Fig. 13 | same data on log abs axis |
| - | deviation impulsive vs low-thrust, 1e-14 to 1e-3 AU | **Fig. 14** | new |
| Fig. 12 | targeting thrust and mass (about 466.8 -> 465.1 kg) | Fig. 15 | same data |
| - | distance to Earth and Mars, 2040-2046 | **Fig. 16** | new |
| Fig. 13 | crosslink, 8 flybys in two groups of 4, Tbit (per-flyby totals about 850-1500 Tbit) | Fig. 17 | **replaced**: 2 flybys, "Sep 2039" and "Jun 2043", Pbit (about 0.8 and 1.38 Pbit), 15 and 26 days |

## 3. Every number that differs

| quantity | AAS 22-015 | journal | journal page |
|---|---|---|---|
| fully converged solutions | "several", one presented | 3 (A, B, C) | 17 |
| STAR set size | not given ("278 of STAR's") | 2225 | 17 |
| maintenance propellant bound | < 5 kg (p.15) | < 2 kg | 17 |
| Solution A COI propellant | 36 kg for the entire mission | 33.23 kg (Table 1); text still 36 kg "for the entire mission" (p.18) and "36 kg" insertion (p.26) | 17, 18, 26 |
| Solution A targeting propellant | about 2 kg | 1.65 kg (Table 1); "about 2 kg" (p.18) | 17, 18 |
| Solution B | - | dep 2028-12-20; flybys 2031-11-11, 2034-12-06, 2038-09-26, 2039-11-08; 65.36 / 1.94 kg | 17 |
| Solution C | - | dep 2028-12-21; flybys 2032-04-04, 2034-08-26, 2039-12-31, 2041-09-15; 83.88 / 1.50 kg | 17 |
| Solution A dates | "mid-2037 to the beginning of 2046" | dep 2037-08-09; flybys 2039-10-06, 2041-08-13, 2043-07-02, 2044-08-23 | 17 |
| flybys of Solution A | 8 (abstract, p.15, conclusion) | 8 (conclusion p.26); Table 1 lists 4; "Two Earth flybys and three Mars flybys" (p.23, matching Fig. 16) | 17, 23, 26 |
| data goal | about 1 Pbit per year; > 8000 Tbits over 8 flybys | about 1 Pbit per flyby; no total | 24 |
| optical terminal | not given | DSOC heritage, 22 cm aperture, 4 W laser, upgraded to > 2.5 Gbps, PPM, CCSDS HPE | 23 |
| optical demo | not given | > 1 Gbps Moon-Earth | 3 |
| authors | 5 | 6 (Hua Xie, JPL) | 1 |

Unchanged numbers (checked in both): 500 kg; NEXT Isp 4155 s, 0.235 N; altitude <= 25,000 km, >= 300 km
(Mars), >= 1000 km (Earth); DSN 120 deg spacing, 34 m and 70 m dishes; mid-June 2035 C3 minimum, about 190 days;
>= 800 days to the target Mars flyby; 3 x Mars SOI; four COI thrust arcs of weeks to months; targeting arcs of
10-20 h; 1-3 Pbit per flyby; 100 kg spacecraft with space tug; NIAC 80NM0018D0004; 278 converged to Step 3.

## 4. Checks (our arithmetic, `check_table1.py`)

- Propellant to delta-V (c = 4155 x 9.80665 m/s): A COI 2802 m/s, targeting 144.3 m/s; B 5708 / 182.3 m/s;
  C 7483 / 147.1 m/s. Fig. 15 mass (about 466.8 -> 465.1 kg) agrees with 1.65 kg. Fig. 9 end mass (about
  467 kg) agrees with 33.23 kg, not with 36 kg.
- 33.23 + 1.65 = 34.88 kg, not 36 kg. Internal inconsistency.
- Departure to Flyby #1: A 788 d, B 1056 d, C 1200 d. A is below the "at least 800 days" rule (p.11) by 12 d.
  Not explained in the paper.
- Flyby spacing for A: 677, 688, 418 d. Fig. 16 (p.24) gives the order Mars, Earth, Mars, Earth (then Mars
  about 2046.1). So both Mars -> Earth intervals are longer than 12 months: 677 d (22.2 months,
  2039-10-06 -> 2041-08-13) and 418 d (13.7 months, 2043-07-02 -> 2044-08-23). The STAR dataset in Fig. 4
  (p.13) is "constrained to be shorter than 12 months" for the Mars-Earth transit. The paper does not say
  whether Solution A came from that filtered set. Observation only.
- Fig. 17 flyby labels (Sep 2039, Jun 2043) are 2-4 weeks before the Table 1 Mars dates (2039-10-06,
  2043-07-02). The plots cover a 15-day and a 26-day approach, so the labels likely mark the window start.

## 5. Citation mining (17 references)

Held/not-held checked with `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX.

| # | reference | held? |
|---|---|---|
| 1 | Sanchez Net, Pellegrini, Vander Hook, "Cycler orbits and the Solar System Pony Express", IEEE Aerospace 2020 (DOI 10.1109/AERO47225.2020.9172342 per AAS ref 1) | not held; not on wanted list. Low priority (concept paper; the JSR 2022 paper supersedes it). |
| 2 | Hemmati & Caplan 2013, optical comms book chapter | not held; out of scope |
| 3 | Byrnes, Longuski, Aldrin 1993, JSR 30(3) | held |
| 4 | Hollister & Menning 1970, JSR 7(10) | held |
| 5 | Sanchez Net et al., "Solar system data mules: analysis for Mars and Jupiter", IEEE Aerospace 2021 (DOI 10.1109/AERO50100.2021.9438463 per AAS ref 5) | not held; not on wanted list. Possible interest: Jupiter data-mule cyclers. Crossref confirms DOI 10.1109/aero50100.2021.9438463. Medium-low priority. |
| 6 | Taheri & Junkins 2018, JGCD 41(11) | not held; method only |
| 7 | Woollands, Taheri, Junkins 2019, JAS 67(2) | held |
| 8 | Walker, Ireland, Owens 1985, MEE | not held; method only |
| 9 | Pontani 2021, JOTA 191 (new in journal) | not held; method only |
| 10 | Biggs & Maclean, IEEE TAES 50 (new in journal) | not held; out of scope |
| 11 | Battin 1999 | not held; textbook |
| 12 | Lawden 1963 | not held; textbook |
| 13 | Maly et al. 2000, ESPA | not held; out of scope |
| 14 | Fisher 2020, NEXT-C | not held; out of scope |
| 15 | Russell & Ocampo 2006, JGCD 29(2) | held |
| 16 | Landau, Campagnola, Pellegrini, "STAR searches for patched-conic trajectories", JAS (submitted in 2022) | not held (Landau 2018 JGCD held is a different paper). Crossref: JAS 69 (2022), DOI 10.1007/s40295-022-00350-y. Worth acquiring: STAR generated the 2225-cycler set behind both Pony Express papers. |
| 17 | Sanchez Net et al., "Cycler orbits and the Solar System Pony Express", JSR | held (`sanchez-net-2022-cycler-orbits-solar-system-pony-express-JSR.pdf`) |

*Filed as `cyclers_pdf/papers/pascarella-woollands-pellegrini-sanchez-net-xie-vander-hook-2023-low-thrust-trajectory-optimization-solar-system-pony-express-acta-astronautica-203-280-doi-10.1016-j.actaastro.2022.11.046-preprint.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*
