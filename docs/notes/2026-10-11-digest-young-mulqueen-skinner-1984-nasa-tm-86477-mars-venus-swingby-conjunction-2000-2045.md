# Digest: Young, Mulqueen and Skinner 1984, "Mars Exploration, Venus Swingby and Conjunction Class Mission Modes, Time Period 2000 to 2045", NASA TM-86477 (#960, from the owner's Mac)

Archie C. Young, John A. Mulqueen and James E. Skinner (Program Development, NASA George C. Marshall Space Flight Center), NASA Technical Memorandum TM-86477,
31 August 1984 (title page and report-date box), 82 pp., NTRS 19850003653 (N85-11961). Approved by William R. Marshall, Director, Program Development (p.82). The first author is Young, as guessed: the same
Archie C. Young as the held 1966 TM X-53511. The author line on the cover page, the standard title page and the approval page all read the same three names.
Filed as `cyclers_pdf/papers/young-1984-mars-exploration-venus-swingby-conjunction-class-mission-modes-2000-2045-nasa-tm-86477-ntrs-19850003653.pdf`
(md5 `dc299ad48e5d125788ac277cb7929930`, 3,809,391 bytes, 82 pp.: `ocrmypdf --force-ocr -O1 -l eng` of the received scan, source md5 `b45edfdcc893e36666f1b78818ff7729`, 3,732,763 bytes, garbled
text layer of 8,912 words; the OCR copy has 28,578 words, but the pages that hold Tables 1 to 31 are printed sideways and their OCR text is unusable, so every table was read on page images and is in the data file). No
Ghostscript recompression: the file did not grow. A second run with `--rotate-pages` did not fix the sideways pages and was discarded.
Corpus check before filing (2026-10-11): `grep` of CORPUS_INDEX.md and `docs/notes/*.md` for "19850003653", "TM-86477", "2000 to 2045" found nothing. Held and related: Young 1966 TM X-53511
(`2026-10-06-digest-young-1966-multiple-planet-flyby-venus-mars.md`), NASA SP-35 1968 Part 6 (`2026-10-11-digest-nasa-sp-35-1968-part-6-mars-stopover-venus-swingbys-douglas.md`, digested in parallel by corpus-titus-sonnet), Sohn 1964,
Hollister and Prussing (AIAA 65-700 and RE-16).
Page numbers below are the report's printed numbers; PDF page = printed page + 8 for pp.9 to 74.
Status: READ in full (OCR text for the body, page renders for every figure and table; all 31 tables, Figs. 2, 7, 14, 15, 19, 24, 25 on renders, the other figures as thumbnails). Numbers:
`data/sources/young-mulqueen-skinner-1984-tm-86477-tables.yaml`. Checks: `cyclers_pdf/papers/...-check.py` and `...-check.out` beside the PDF (output summarised in sec. 7).

## 0. Verdict

**A set of mission-opportunity tables for ONE-SHOT manned Mars round trips, 2000 to 2045. No trajectory in it repeats and none is called periodic.** Two mission families:
- Opposition class with a 60-day Mars stay and ONE or TWO Venus swingbys: 15 launch opportunities 2001 to 2031 (Table 1; 6 inbound, 7 outbound, 2 double) with Earth-launch dates, hyperbolic excess speeds at every body, Venus closest-approach radius, stage sizing for chemical LOX/LH2 vehicles.
- Conjunction class (long stay, 340 to 550 d) with no swingby: 7 launch opportunities 2031 to 2044 (Table 22), covering one Earth-Mars synodic period.
It extends Young 1966 (1975 to 1980) and SP-35 Part 6 (1975 to 1999) to 2001 to 2031 and adds the double-swingby mode (SP-35 excludes it) and the conjunction mode. The only repetition statements are
about patterns of launch opportunities (a 15-year or seven-opposition cycle of the Mars-eccentricity profile; conjunction requirements "approximately repetitive every seven opportunities"). The numbers in it are consistent
at the 0.1 percent level between independent tables (sec. 7), with a few listed exceptions.

## 1. Introduction, assumptions, opportunities (pp.1-4)

- Purpose: performance and operation requirements and their sensitivity for the Venus-swingby and conjunction-class stopover missions; initial mass in Earth orbit per opportunity; "a precursor mission for the long stopover time, up to 550 days".
  Direct opposition-class round trips (stay up to 60 d) need very high propulsive energy; a Venus swingby (accelerating outbound, decelerating inbound) reduces it; the conjunction class gets low energy by optimising the stay.
- Assumptions (Fig. 1, p.11, and p.1-2): vehicle assembled in and departing from a 270 n.mi circular Space Station orbit; Mars capture and escape orbit 24-hr ellipse with 270 n.mi periapsis; Earth return orbit 24 hr, 270 n.mi perigee; three
  propulsion stages (Earth escape; Mars braking and escape; Earth braking), each mass fraction 0.90, LOX/LH2, Isp 480 s ("relatively low" mass fraction to cover micrometeoroid protection and boil-off); Venus closest approach not below 0.1 planet radius altitude (330 n.mi);
  conjunction class uses type I or II transfers with the stay optimised for minimum initial weight. Weights: Venus-swingby mission module 83,000 lb, probes 42,000 lb, Mars excursion module 95,000 lb; conjunction mission module 117,000 lb, same excursion module, Pioneer Mars base 57,000 lb.
- Opportunities (p.2, Fig. 2 viewed): standard opposition-class missions lie 90 to 180 days before opposition, oppositions every 26 months on average; "the cyclic pattern of mission profile variation repeats every 15 years or every 7 oppositions [ref 1]" because of Mars' eccentricity; Fig. 2 plots the opposition positions around Mars' orbit for two cycles 2000 to 2031 (labels 2001, 2003 at Mars perihelion, ..., 2031).
  Why the swingby: a near-Hohmann outbound leg arrives at Mars after opposition so a near-Hohmann inbound leg is impossible; the heliocentric transfer must exceed 180 degrees; Venus either shifts the outbound geometry or accelerates the inbound.
- Availability test (pp.3-4, Figs. 3 to 6): Venus is on the vehicle's way on one leg; Venus's angular rate is high so it is usually available on one leg; the Sun-Venus-Mars alignment about 200 days before an Earth-Mars opposition (and Earth-Venus conjunction 200 days earlier) selects outbound; Figs. 4 to 6 plot the dates of Earth-Venus alignment, Venus-Mars alignment and
  opposition continuously 2000 to 2031, with a YES / NO / MAYBE table per opposition for outbound and inbound availability. The symmetric alignment about an opposition "occurs every third opposition".

## 2. Mission characteristics (pp.4-5, Table 1, Fig. 7, Table 22, Fig. 19)

- **Table 1 (data file; 15 rows)**: for each opposition-class mission: Earth-Mars opposition day, launch, Venus pass (outbound), Mars arrival, Mars leave, Venus pass (inbound), Earth arrival, each with V_inf in EMOS in brackets, and total mission time. Heading "J.D. 2450000 or 2460000"; the numbers behave as MJD (data file header: Figs. 14 and 24 dates are reproduced exactly by 1858-11-17 + day number).
  Modes by opposition: 2001 in, 2003 out, 2005 out, 2007 in, 2010 double, 2012 out, 2014 in, 2016 in, 2018 out, 2020 in, 2022 out, 2025 out, 2027 in, 2029 double, 2031 out. Mars stay 60 d in every row. Trip times 558 to 737 d (shortest 558 d, September 2007, inbound; longest 737 d, March 2028, double).
  Earth launch precedes opposition by 3 to 6 months (inbound) or 12 to 19 months (outbound). Minimum V_inf at Earth and Mars about 0.10 EMOS for near-Hohmann legs; the 2001 inbound is cheap at Earth departure and Mars arrival, 0.2050 EMOS at Mars departure, 0.1418 at Earth arrival.
- **Fig. 7 (p.17)**: the same 15 missions listed by Earth launch month and trip time: March 2001 606; August 2002 610; June 2004 659; September 2007 558; January 2009 736; November 2010 650; November 2013 634; November 2015 577; April 2017 638; June 2020 594; October 2021 636; September 2023 614; November 2026 570; March 2028 737; January 2030 654. NOTE: Table 1 labels a row by the opposition year, Fig. 7 by the launch.
- **Table 22 (data file; 7 rows)** and **Fig. 19 (p.29)**: conjunction class, opposition years 2031, 2033, 2035, 2037, 2040, 2042, 2044. Outbound 200 to 356 d, stay 340 to 550 d, inbound 200 to 350 d, total 950 to 1004 d (shortest April 2033 launch, stay 550; longest June 2035, stay 530). V_inf about 0.08 to 0.14 EMOS everywhere. The
  stay of 340 d is the lower limit in four of seven rows (2037 to 2044). "Covers the envelope for the conjunction class since it covers one synodic period."
- Text gloss on the longest and shortest trips of Table 22 ("shortest trip time ... 200 days" for the 2033 legs; "longest ... 356 days" for the 2037 outbound leg) agrees with the table.

## 3. Space vehicle, trajectory program, stage sizing (pp.5-8)

- Spacecraft (p.5-6, Fig. 8): crew of six, Mars excursion module takes three to the surface for a 30-day stay (2 to 5 days of orbital survey first), ballute entry, ascent module; mission module 7 psia O2/N2 shirtsleeve; probes and orbiters.
- Trajectory program (p.7): restricted two-body, patched-conic; "while within the sphere of influence of Venus the vehicle is on a free-flight hyperbola ... no change of energy with respect to Venus ... magnitude of the velocity on leaving must equal that on entering; if the required bend angle is too large to achieve by constraining the periapsis altitude to one tenth of a planet radius, a propulsive maneuver is effected". Tables 2 to 18 are the
  launch-window trajectory printouts for the Venus-swingby missions and Tables 23 to 29 for the conjunction missions; the appendix (pp.67-73) defines the printout columns.
- Stage sizing (Figs. 9 to 12, 20 to 22; Tables 19 to 21 and 30, 31; data file `stage_sizing`): the Venus-swingby WO (initial weight in low Earth orbit) ranges 1,470,725 lb (2001 inbound) to 4,064,865 lb (2007 inbound); with the unfavourable opportunities set aside it does not exceed 3.0 million lb. First stage 869,021 to 2,719,261 lb, second 352,989 to 1,093,887 lb. A third-stage range of "24,000 to 102,000 lb" is stated (p.8, Fig. 15) while Table 21 gives 17,513 lb for 2025 (sec. 7).
  Venus-swingby missions are "two favourable out of three" opportunities (conclusion). Double swingbys "available every third opportunity" and "generally higher energy; exception 2007" (p.8); Tables 1 and 19 to 21 show only two double missions (2010, 2029), and 2007 inbound, not 2010, is the most expensive.
  The 2025 outbound vehicle (Fig. 13): total length 276.5 ft, WO 2.08 million lb (matches Table 21: 2,084,630 lb). Conjunction (Figs. 21, 22, 23): first stage 600,000 to 725,000 lb, second 112,000 to 169,000, third 25,000 to 36,500 lb; WO 1.01 to 1.17 million lb; 2031 vehicle length 211.5 ft, 1.17 million lb. Earth departure C3 8.9 to 16.2 km2/s2 (conjunction), 9.2 to 37.0 km2/s2 (Venus swingby, Fig. 15).
- Representative profiles (pp.8-9, Figs. 14, 24): the 2014 inbound swingby (Earth launch 11/22/13, Mars arrival 8/19/14, Mars departure 10/16/14, Venus passage 3/2/15, Earth arrival 8/19/15; Earth about 45 degrees ahead of Mars at Mars arrival; near-tangential departure and arrival, 634 d), the 2025 outbound swingby (launch 9/9/23, Venus 2/18/24, Mars arrival 7/3/24, departure 9/1/24, Earth 5/15/25;
  transfer angle Earth to Venus over 180 deg, Venus to Mars under 180 deg, total just over 360 deg; 614 d), and the 2031 conjunction mission (launch 12/25/30, Mars arrival 10/3/31 after 282 d, stay 500 d, departure 2/14/33, Earth 9/18/33 after 216 d; 998 d).
- Conclusions (pp.8-9, Figs. 15, 25): an outbound or inbound swingby is possible at every opposition; two out of three should be favourable; mission time increases 20 to 50 percent over direct flight; the 60-day stay "may be extended to over 100 days" using both swingbys in favourable opportunities; a swingby vehicle and a direct vehicle could rendezvous at Mars; conjunction requirements "approximately repetitive every seven opportunities" (the 2029 opportunity resembles 2044).

## 4. The printout tables (Tables 2 to 18 and 23 to 29; data file `printout_records`)

Each table is a launch-window sweep of typically 7 records at 2-day launch steps around the optimum date (the optimum sits at a different position in each sweep): 5 or 6 lines per record for the swingby missions (three JD lines for inbound or outbound legs, a Venus passage line, communication distances, DV and weight-ratio line)
and 4 lines for the conjunction missions. Columns are defined in the appendix (DEPART, PASS, ARRIVE day numbers, C3-DD, RA, DECL, I1, V1, PSI1, ECCEN, SMA, THET1, THET2, PERIH, APHEL, V2, PSI2, hyperbolic excess speed, RCP closest-approach radius in planet radii, KAPPA bend angle, UP pericentre speed, AH, EH, ANMAX, INH, RAP/DECP, RAS/DECS, ETA, lighting condition, DVCP power maneuver at the swingby planet, CORD1 to 3 communication distance in AU, DV1 to DV4, GAMA, ETA).
The optimum-date record of each of the 24 tables is transcribed (one witness: the render read by eye; Tesseract as a second read on the digits; the dates and C3 columns were also checked against Table 1 and Table 22, sec. 7). **The other six records per table (about 144 records) are not transcribed**: they are sampling of the launch window, about 6 numbers per cell by 100 cells, and add
nothing the project needs beyond the optimum; this is a scoping decision, stated here, not an omission by oversight. The table gives, for the optimum launch date, the Venus passage radius, bend angle and pericentre speed, in sec. 5.

## 5. Venus passage and V_inf at the optimum dates (Table 1 and printouts; checks in sec. 7)

Closest-approach radius RCP (planet radii; altitude = RCP - 1 times 6051.8 km), bend angle, pericentre speed UP (EMOS), V_inf at Venus (EMOS, in the table; arrival side), powered maneuver DVCP (km/s). Dates are from the day numbers read as MJD (sec. 0 note).

| Mission (opposition, mode) | Launch (MJD) | RCP | Altitude, km | Bend, deg | UP, EMOS | V_inf at Venus, EMOS | DVCP |
|---|---|---|---|---|---|---|---|
| 2001 inbound | 2001-03-28 | 1.106 | 641 | 53.6 | .4141 | .2531 | 0 |
| 2003 outbound | 2002-08-24 | 1.143 | 865 | 75.7 | .3683 | .1781 | 0 |
| 2005 outbound | 2004-06-08 | 2.261 | 7631 | 22.9 | .3968 | .3259 | 0 |
| 2007 inbound | 2007-09-07 | 1.107 | 648 | 38.0 | .4659 | .3308 | 0 |
| 2010 double, outbound leg | 2009-01-23 | 1.100 | 605 | 56.5 | .4124 | .2463 | .21559 |
| 2010 double, inbound leg | (Mars leave 2010-03-01) | 1.271 | 1640 | 56.6 | .3818 | .2370 | 0 |
| 2012 outbound | 2010-11-28 | 1.211 | 1277 | 32.2 | .4734 | .3570 | 0 |
| 2014 inbound | 2013-11-22 | 1.959 | 5804 | 24.7 | .4136 | .3321 | 0 |
| 2016 inbound | 2015-11-10 | 1.665 | 4024 | 21.2 | .4769 | .3958 | 0 |
| 2018 outbound | 2017-04-11 | 1.706 | 4273 | 27.5 | .4240 | .3345 | 0 |
| 2020 inbound | 2020-06-10 | 1.507 | 3068 | 39.2 | .3942 | .2764 | 0 |
| 2022 outbound | 2021-10-31 | 2.676 | 10143 | 40.3 | .2933 | .2012 | 0 |
| 2025 outbound | 2023-09-09 | 1.764 | 4624 | 22.3 | .4542 | .3758 | 0 |
| 2027 inbound | 2026-11-02 | 1.212 | 1283 | 33.7 | .4650 | .3440 | 0 |
| 2029 double, outbound leg | 2028-03-26 | 1.741 | 4484 | 41.0 | .3613 | .2460 | 0 |
| 2029 double, inbound leg | (Mars leave 2029-05-14) | 1.471 | 2850 | 49.9 | .3675 | .2374 | 0 |
| 2031 outbound | 2030-01-29 | 1.226 | 1368 | 32.9 | .4683 | .3516 | 0 |

Range over the 15 missions: Venus V_inf .1781 to .3958 EMOS (5.3 to 11.8 km/s at 29.78 km/s); closest approach 605 to 10,143 km altitude (RCP 1.100 to 2.676); bend 21 to 76 degrees. The 1.1 limit is reached by the 2010 outbound record (RCP 1.100) and nearly reached by 2001 and 2007 (1.106, 1.107); only the 2010 outbound record has a powered maneuver (0.2156 km/s) printed.
Earth-launch V_inf .0947 to .1957 EMOS (launch C3 9.2 to 34.0 km2/s2; Fig. 15's "37.0" is not an Earth-launch C3 of any optimum record, sec. 7).

## 6. Figures (one line each)

1A assumptions block (viewed as a rotated image; text in sec. 1). 2 opposition positions around Mars' orbit (viewed). 3 outbound and inbound conditions for Venus availability, plus the direct stopover profile (thumbnail). 4, 5, 6 availability plots, Earth-Venus alignment, Venus-Mars alignment and opposition dates against Julian date with YES / MAYBE / NO tables (thumbnails; not read). 7 mission list (viewed, transcribed). 8 interplanetary spacecraft (thumbnail). 9 to 12 stage weights and weight in Earth orbit as bar charts per opposition (thumbnails; values equal Tables 19 to 21). 13 vehicle sized for 2025. 14 the 2014 inbound and 2025 outbound profiles (viewed; dates in sec. 3). 15 summary bullets (viewed).
16 conjunction assumptions. 17 oppositions 2030 to 2061 (thumbnail). 18 the direct inverted stopover profile (thumbnail). 19 conjunction mission list (viewed, transcribed). 20 spacecraft. 21, 22 conjunction stage weights and weight in orbit (thumbnails). 23 vehicle sized for 2031. 24 the 2031 conjunction profile (viewed). 25 summary bullets (viewed). 26 transfer geometry (appendix, p.68).

## 7. Independent checks (ours; `-check.out`)

- **Table 1 arithmetic:** launch + total trip = Earth arrival and Mars stay = 60 d in every row, except 2029 (arrival 62592 minus launch 61856 is 736 d; Table 1 and Fig. 7 print 737; the printout Table 17 prints 62592.0 as well). The launch months of Table 1 equal Fig. 7's for all 15 rows (months computed from MJD).
- **Table 22:** the legs add up (outbound + stay + inbound = total) in all 7 rows; the launch months and opposition months equal Fig. 19's.
- **V_inf of Table 1 against the printout C3 columns:** 65 comparisons (Earth launch, Venus, Mars arrival and departure, Earth arrival, for the 13 non-double rows) give an implied EMOS of 29.780 km/s (mean; range 29.55 to 29.80). One outlier: 2018 outbound Mars arrival, Table 1 prints .1865 but the optimum-date record has C3 30.37 (V_inf 5.51 km/s = .1850 EMOS); the Table 1 digit may be 5 for 0 or the record may differ; unresolved.
- **Venus hyperbola:** in 15 of the 17 records AH (EH - 1) equals RCP to 0.002 and the bend angle equals 2 asin(1/EH) to 0.1 degree. The two exceptions are the two double-swingby legs of 2010 (bend 56.5 and 56.6 against 54.8 and 55.6 degrees from EH): the 2010 outbound record is the powered one, the inbound shows unequal arrival and departure speeds (below).
- **Venus speed match:** the text says the flyby conserves speed (V_inf in = V_inf out) and applies a propulsive maneuver only when the 0.1-radius altitude limit binds. In the printed records V_inf at arrival and departure differ by 0.005 to 0.41 km/s (up to 11 percent in C3), for example 2003 outbound 5.304 vs 5.428 km/s, 2022 outbound 5.992 vs 6.206, 2010 inbound 7.060 vs 6.654, with DVCP printed 0 in 16 of 17. The report does not explain the mismatch (SP-35 Part 6's program matched the speeds to 3 m/s). It may reflect the program's tolerance or definitions; not resolved. A reproduction from these tables should not assume an exactly ballistic flyby.
- **Stage sizing (Tables 19 to 21, 30, 31):** in all 22 blocks WSTG = WF + WBO to 2 lb, WBO = 0.1111 WF to 2.5 lb, WO = the stages plus payloads (83,000 + 95,000 + 42,000 lb for swingby; 117,000 + 95,000 + 57,000 lb for conjunction) to 2 lb, and DV1 equals Isp 480 s times g0 times ln(WO/(WO - WF1)) to 0.01 km/s. Digits reconstructed from these relations (smudged on the scan): 2003 WSTG3 70203 (printed with a slashed zero that reads as 8), 2007 WF1/WBO1, 2012 WSTG1, 2020 WO and DV4, 2040 WF2 (data file note). One range statement fails: the stated minimum third stage 24,000 lb (p.8, Fig. 15) is below Table 21's 17,513 lb (2025, WF3 15,762, Earth-braking DV4 0.803 km/s); the maximum 102,000 is Table 21's 101,967 (2031) and 100,164 (2022). Fig. 15's "9.2 to 37.0" Earth-departure C3 has an upper bound that no optimum record reaches (largest 33.98, 2016).
- **Dates:** the dates printed on Figs. 14 and 24 are reproduced exactly by day number + 1858-11-17 for the 2025 outbound mission and the 2031 conjunction mission (nine dates), and, in the 2014 inbound figure, for the launch, Mars arrival and Venus passage; its Mars departure (10/16/14) and Earth arrival (8/19/15) are 2 d and 1 d from the Table 1 day numbers (10/18/14, 8/18/15). The day numbers therefore behave as MJD, not as JD as headed.
- **The 6.4-year configuration repeat (SP-35 Part 6, Young 1966):** shifting each opposition of Table 1 by 2338 d lands within 55 d of a later opposition in 12 cases, and the swingby side (outbound or inbound; a double counts as both) is kept in 11 of the 12 (the exception: 2016 inbound to 2022 outbound). The report itself does not use the 6.4-year number; it states the 15-year, seven-opposition cycle of Mars-eccentricity effects. The 15-year shift does not keep the swingby side (2005 outbound to 2020 inbound, 2007 inbound to 2022 outbound).

## 8. Relevance to this project

- **The 2000-2045 opportunity table against the Earth-Venus and Earth-Mars rows.** Table 1 and Table 22 give dates and V_inf at every body for one-shot E-V-M-E (outbound), E-M-V-E (inbound) and E-V-M-V-E (double) sequences with a 60-day Mars stay, and E-M-E conjunction round trips. Venus passage distances are given (sec. 5; RCP 1.100 to 2.676). No row is a repeating trajectory, so none is a collision for a cycler candidate; they are prior art for the *encounter order* of a Venus-swingby Mars round trip in 2001 to 2031 (R1(c), same
  class as Young 1966's 1977 E-V-M-V-E, 1978 E-V-M-E and 1975 E-M-V-E). V_inf at Earth launch 0.10 to 0.20 EMOS and at Earth arrival 0.09 to 0.28, at Venus 0.18 to 0.40, at Mars arrival 0.12 to 0.22 and departure 0.12 to 0.24 EMOS; the Earth-Mars conjunction legs sit at 0.08 to 0.14 EMOS, the Hohmann-like end of the range.
- **Does any mode repeat or is any called periodic?** No trajectory. The text says the profile pattern of Earth-Mars oppositions "repeats every 15 years or every 7 oppositions" (Mars eccentricity; ref 1), Venus availability is "every opposition, favourable two out of three" and conjunction-class requirements are "approximately repetitive every seven opportunities". These are patterns of opportunities, as SP-35 Part 6 and Young 1966 say for the 6.4-year repeat. No periodic orbit, no cycler, no multi-lap statement.
- **How it extends Young 1966 and SP-35 Part 6.** Young 1966: 1975 to 1980, four missions with three-body sequences, constraints v_inf Earth and Mars at most 0.30 EMOS, Venus at most 0.50, total at most 750 d. SP-35 Part 6 (1968, Douglas): 16 opportunities 1975 to 1999, outbound or inbound only (no double), charts rather than tables (its supplements, on NTRS, carry the values), patched conic with speed match to 3 m/s. TM-86477 (1984) is the next
  slice in time: 15 opportunities 2001 to 2031 with a table of V_inf, Venus passages and sized stages, the double swingby re-introduced (2010 and 2029; Young 1966's 1977 triple is the only earlier one), and conjunction missions. Its Mars stay is fixed at 60 d, where SP-35 charts 0, 30 and 60 d. A cross-check to do when needed: SP-35 Part 6's 1997 and 1999 outbound and 1995 inbound opportunities map by the 2338-day shift onto TM-86477's 2003 outbound, 2005 outbound and 2001 inbound
  (opposition-day gaps 17, 52 and 25 d; sec. 7), so the SP-35 tabulated Venus passage distances and speeds for those three are a direct comparison for Tables 3, 4 and 2 here. Not done here (the SP-35 values are in its supplements, not in the held PDF).
- **Gate label.** There is no gate anchor for this report in literature_check.py. If one is added, the content supports {mga-tour} (one-shot Venus-swingby Mars round trip) for the opposition class and {transfer} for the conjunction class, with Venus (and Mars) as flyby bodies; no cycler, repeated-moon or resonant label.
- **Real-ephemeris use.** The optimum dates and V_inf of Table 1 and the Venus passage radii can seed a real-ephemeris recomputation of the 2001 to 2031 Venus-swingby rows, but the program is patched-conic with the Table 2-1-style fixed ellipses (not stated here); the Venus speed mismatch (sec. 7) means the printed Venus V_inf is not a ballistic value to the 5 m/s level.

## 9. Cautions

- Day numbers are MJD in behaviour though headed JD; Table 22's headings read "J.D. 2460000" and "J.D. 246000".
- Table 1 labels by opposition year, Fig. 7 by launch year; the same mission has two years.
- A few Table 1 digits were smudged on the scan (2003 Venus pass, 2022 launch, 2029 Earth arrival); they are fixed from the printouts and arithmetic, as stated.
- The printouts' unreadable columns on Table 5 (a stripe of the scan) are marked [ill].
- Only the optimum-date record of each printout table is transcribed (sec. 4).

## 10. Citation mining

The reference list has five items (p.10): (1) Sohn 1964 (held, digest); (2) NASA SP-35, Space Flight Handbook Vol. 3, Planetary Flight Handbook, prepared for MSFC, 1963, cited for single-leg trajectory data ("Examination of single leg trajectory data [2] indicates that if the outbound and inbound legs could be optimised separately ..." for the conjunction class) (the 1963 Part 1; SP-35 Part 6 1968 is held; Part 1 is not, see that digest);
(3) Hollister and Prussing, AIAA/ION conference 16 September 1965 (held, AIAA 65-700 and RE-16); (4) Deerwester AIAA 65-89 (wanted-list row 41, not held); (5) Integrated Manned Interplanetary Spacecraft Concept Definition, Vol. 1, The Boeing Company, NASA CR-66558, contract NAS1-6774, January 1968 (the source of the spacecraft element weights; not held, not on the list; a Boeing study of Mars and Venus manned missions with flyby). No new wanted-list rows are
proposed by this digest, as the lead said: the Deerwester row exists, the SP-35 Part 1 question belongs to the SP-35 digest, and CR-66558 is spacecraft design, not trajectory.
