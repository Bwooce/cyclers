# Digest: NASA TM X-53049 (1964), "Proceedings of the Symposium on Manned Planetary Missions, 1963/1964 Status" (#960 batch 32)

NASA George C. Marshall Space Flight Center, Future Projects Office (Research and Development Operations),
"Proceedings of the Symposium on Manned Planetary Missions, 1963/1964 Status", NASA TM X-53049, 12 June 1964.
Symposium held at MSFC, Huntsville, 28-30 January 1964 (chairman J. N. Smith; panel moderator H. H. Koelle).
NTRS 19640017065 (accession N64-26 9xx on the cover; the digits are not legible). No DOI.
- Source file: `19640017065.pdf` (from macpro), 750 PDF pages, 61,278,066 bytes, md5 094e3558a8ed7357f1e36f603de1d8e3.
  Producer "Pdf.Capture 5.2"; it carries the NTRS OCR text layer.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/nasa-msfc-1964-proceedings-symposium-manned-planetary-missions-1963-1964-status-nasa-tm-x-53049-ntrs-19640017065.pdf`
  (Sohn's Part 5 is the item other papers cite. If the filer prefers an author key:
  `sohn-ehricke-et-al-1964-...`. I suggest the corporate key, because 17 parts have 15 authors.)
- **How I read it:**
  - Front matter (TOC, foreword) on the page images, PDF pp. 2-5.
  - Search on `pdftotext -layout` (40,668 lines; deleted after use, re-make with `pdftotext -layout`) for Venus, swingby, swing-by, round trip,
    free return, periodic, symmetric, repeat, cycle, Minovitch, Ross, Gillespie, Ehricke, Hollister, Crocco,
    "32-464", fly-by, perihelion, opposition, conjunction, shuttle. Hit list in `hits.txt`.
  - Every hit page that touches trajectories was read on the page image (100-450 dpi). Each number below is
    from an image, except where a line says "(text layer)". Pages read on images: PDF 2-5, 16-17, 21-22, 26-28,
    38, 82-83, 114, 155-156, 166-168, 218, 230-231, 241, 270, 272, 274, 286, 615, 625. Tables 1-1, 1-2, 1-3
    and Figs. 9 and 16 were re-read at 400-450 dpi crops.
  - Arithmetic checks: `checks_tmx53049.py`, output `checks_tmx53049.out`. All disagreements are listed in sec. 5.
  - Parts 7, 10-14 and 16 (entry bodies, landers, subsystems, programme models, electric propulsion) were skimmed
    on the text layer only. None has trajectory content of cycler interest.
- **PDF page offset is not constant.** The scan drops blank pages, so PDF page = printed page + 5 at the start
  and printed page - 13 at the end. Always quote both.

## 0. Verdict

**A 1964 conference record of the US industry manned-Mars studies. It holds no cycler, no periodic or repeating
trajectory, and no citation of Minovitch.** Its value to the project is history for `#942` R1(c) (Venus-Mars
round trips):
- **Part 5 (R. L. Sohn, STL) is the "NASA TM-53049, Pt. 5" that Gillespie & Ross 1967 cite as their ref. 2.**
  It is the fuller, earlier version of the Venus swingby return mode that Sohn published as the JSR 1(5) note
  later in 1964 (held). It adds Fig. 9: direct and swingby Earth-entry speeds and trip times for all 14 Mars
  opportunities 1971-1999. See sec. 2.
- **Part 2 (K. A. Ehricke, General Dynamics/Astronautics, contract NAS8-5026) is the Ehricke study in this TM.**
  It has "powered fly-by" (PFB) at Venus on the Earth-bound leg of a Mars capture mission (Mission IV, 1975)
  and on the outbound leg (Mission VI, 1977), plus "bi-planet" Mars-Venus missions. Ehricke says bi-planet
  capture missions "will represent the actual future mode of interplanetary transportation". That is a
  one-sentence vision, not a periodic orbit. See sec. 3.
- Part 15 (Ehricke again, BPTSM, NAS8-11084) has a mission-class code that includes tri-planet missions and
  a "084 Heliocentric Orbit Installation" activity class. It gives no orbit for it. See sec. 4.
- Lockheed (Part 3, EMPIRE follow-on), North American (Part 6) and Ames (Part 8) add single-flyby and
  Venus-encounter data points (sec. 4).
- **Ross and Gillespie appear only in the external distribution list** (PDF p.744: "Dr. S. Ross, Code MT-2";
  "R. W. Gillespie, Code MTG", both NASA HQ). They are not authors here.
- **Minovitch, JPL TR 32-464, Hollister, Crocco: zero hits** on the text layer. The OCR is fair on typeset text
  (sec. 6), and the reference lists are few and short (sec. 7), so I judge this a true negative.
- **Catalogue implication: none.** Nothing here is periodic. PROPOSAL only: when the R1(c) history note is next
  touched, cite Part 5 as the earliest full statement of the swingby mode (January 1964 talk; TM June 1964),
  ahead of the JSR note (September-October 1964).

## 1. Table of contents (read on the image, PDF pp.4-5, printed pp. iii-iv)

| Part | Author, affiliation (contract from the part title page) | Title | Printed p. | PDF pp. |
|---|---|---|---|---|
| 1 | W. von Braun, Director, MSFC | Welcome to MSFC | 1 | 6-10 |
| 2 | K. A. Ehricke, General Dynamics/Astronautics (NAS8-5026) | A Study of Manned Interplanetary Missions | 7 | 11-77 |
| 3 | B. P. Martin, Lockheed Missiles and Space Co. (NAS8-5024) | Manned Interplanetary Missions | 75 | 78-143 |
| 4 | C. A. Syvertson, Ames Research Center | Ames Research Center Mars Mission Studies Introduction | 141 | 144-150 |
| 5 | R. L. Sohn, Space Technology Laboratories (NAS2-1409, per p.157 footnote) | Summary of Manned Mars Mission Study | 149 | 151-221 |
| 6 | A. L. Jones, North American Aviation (NAS2-1408) | Manned Mars Landing and Return Mission Study | 221 | 222-244 |
| 7 | R. N. Worth, Northrop Corp. (NAS2-1411) | Maneuverable Descent Systems for Mars Landing | 245 | 245-266 |
| 8 | H. Hornby, Ames Research Center | Ames Research Center Mars Mission Studies Summary | 269 | 267-275 |
| 9 | R. N. Austin, General Dynamics/Fort Worth (NAS8-11004) | A Study of Manned Mars Exploration in the Unfavorable Time Period (1975-1985) | 279 | 276-342 |
| 10 | R. L. Gervais, Douglas Aircraft Co. (NAS8-11005) | Manned Mars Exploration in the Unfavorable (1975-1985) Time Period | 347 | 343-399 |
| 11 | A. L. Jones, North American Aviation (NAS9-1748) | Subsystems Requirements for a Mars Mission Module | 405 | 400-435 |
| 12 | F. P. Dixon, Philco Corp. (Aeronutronic) (NAS9-1608) | Study of a Manned Mars Excursion Module | 441 | 436-517 |
| 13 | D. J. Shapland, Lockheed Missiles and Space Co. (NAS9-1702) | Preliminary Design of a Mars-Mission Earth Reentry Module | 525 | 518-568 |
| 14 | G. W. Morgenthaler, Martin Co. (NAS8-11057) | A Planetary Transportation Systems Model and Its Application to Space Program Planning | 577 | 569-600 |
| 15 | K. A. Ehricke, General Dynamics/Astronautics (NAS8-11084) | A Study of the Development of a Basic Planetary Transportation Systems Model | 611 | 601-674 |
| 16 | B. Pinkel and D. L. Trapp, The RAND Corp. (NAS8-11081) | A Study of Electrical Propulsion in the Space Program | 687 | 675-723 |
| 17 | Panel, moderated by H. H. Koelle, Director, Future Projects Office, MSFC | Panel Discussion | 737 | 724-740 |
| - | - | Approval page; internal and external distribution | - | 741-750 |

- Part titles are as in the TOC. Part 14's own title page drops "Systems" ("A Planetary Transportation Model...").
- Session chairmen (foreword, PDF p.2): H. O. Ruppe (MSFC), C. A. Syvertson (Ames), J. N. Smith (MSFC),
  C. R. Darwin (MSC), V. Gradecak (MSFC).

### One-line summaries

1. Von Braun: welcome; a Saturn I scrub story as a lesson that planetary launch windows are unforgiving.
2. Ehricke: second phase of the GD/A early manned planetary mission study for MSFC. Mission analysis
   (Venus elliptic capture, Mars capture, perihelion braking, powered fly-by, synodic and bi-planet missions),
   nuclear interplanetary vehicles, structures, life support, radiation, entry, operations, cost. (Sec. 3.)
3. Martin (Lockheed): EMPIRE follow-on summary. Three-man single-planet flybys (1974 Venus, 1975 Mars) on
   Saturn V class boosters; radiation, guidance, reconnaissance, power, spacecraft design. (Sec. 4.)
4. Syvertson: introduces the three Ames contracts (STL, NAA, Northrop).
5. Sohn (STL): Mars stopover mission trades 1971-2000; the Venus swingby return mode; aero entry, navigation,
   launch holds, artificial gravity; Appendix A on analytic optimisation. (Sec. 2.)
6. Jones (NAA): Mars landing-and-return mission and module design; a 1970 Venus-encounter case. (Sec. 4.)
7. Worth (Northrop): maneuverable descent and landing systems for Mars and Earth.
8. Hornby (Ames): summary of the Ames studies; radiation, aerobraking, and the STL Venus flyby result. (Sec. 4.)
9. Austin (GD/FW): Mars missions in the "unfavourable" 1975-1985 period; mission maps; the 16-year
   opposition cycle; plane-change trajectories. (Sec. 4.)
10. Gervais (Douglas): the parallel MSFC study of the same 1975-1985 period; vehicles and operations.
11. Jones (NAA): subsystem design criteria for the Mars mission module.
12. Dixon (Philco): Mars Excursion Module preliminary design.
13. Shapland (Lockheed): Earth reentry module for Mars return speeds.
14. Morgenthaler (Martin): a planetary transportation systems model for programme planning.
15. Ehricke (GD/A): Basic Planetary Transportation Systems Model (BPTSM): structure, codes, sub-models. (Sec. 4.)
16. Pinkel and Trapp (RAND): electric propulsion in the space programme (lunar transport, probes, Mars).
17. Panel: consensus that a manned flyby is the sensible first planetary step; unmanned data are needed.

## 2. Part 5, Sohn (STL): the Venus swingby return mode

Read on images PDF pp.155-156, 166-168, 218 (printed 153-154, 164-166, 216).

- **Summary (printed p.154):** "In a typical unfavorable year, 1975, Earth entry velocity is reduced from 66 to
  46,000 fps." "From 1971 to 1999, which covers two cycles of oppositions, the Venus return mode is achievable in
  9 of 14 opportunities." Mars moves about four times slower than Venus, so a short Mars wait gives a Venus
  phasing. An outbound swingby (Venus "gravitational whip effect") covers the rest. "The Venus swingby mode can
  be used in all opportunities during the thirty-year period", with entry below 45-50,000 fps.
- **1975 example (printed p.165, Fig. 8):** direct return crosses inside the Venus orbit; at Earth the heliocentric
  flight path angle relative to Earth is 31 deg; entry 65,600 fps. A dark-side Venus pass at about 3300 km
  altitude cuts heliocentric speed by 15,000 fps; the Earth approach angle is 14 deg; entry 44,000 fps.
  - Fig. 8 shows the 31 deg and 14 deg triangles. Its two speed labels are too blurred to read; I quote none.
- **Coverage (printed p.166):** "about 75 percent of the missions" can use the return swingby; the rest use the
  "reversed" mode (long >180 deg transfer out, Venus whip, short <180 deg back), entry below 50,000 fps.
  "A Venus swingby mode can be achieved in all 14 of the oppositions to 1999." Venus encounter before the
  spacecraft's perihelion stretches the trip to about 500 days; after perihelion, little or no stretch.
- **Fig. 9 (PDF p.218, 450 dpi):** bars of direct and swingby Earth entry speed for 1971, 73, 75, 78, 80, 82, 84,
  86, 88, 90, 93, 95, 97, 99. Gridlines at 40, 50, 60, 70 kfps, with an "Apollo" level below 40.
  - Outbound-swingby dots over **1973, 1980, 1986, 1993, 1999** (five). No dot is visible over 1978.
  - Trip durations under the bars, two rows (the parenthesised row is, by the text, the swingby trip; the
    figure does not say so):

    | Year | 1971 | 73 | 75 | 78 | 80 | 82 | 84 | 86 | 88 | 90 | 93 | 95 | 97 | 99 |
    |---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
    | days | 400 | 425 | 430 | 430 | 440 | 430 | 432 | 382 | 420 | 432 | 425 | 415 | 418 | 407 |
    | (days) | (440) | (490) | (490) | (435) | (520) | (520) | (458) | (470) | (518) | (442) | (500) | (516) | (450) | (490) |

  - Bar heights (read against gridlines, about +-1 kfps): direct entry tops out near 69 kfps (1978) and
    68 kfps (1993); swingby entry is 40-50 kfps in every year, highest about 50 kfps (1997).
    Control: the 1975 bars read about 66 kfps direct and about 44 kfps swingby, which matches the text
    (65,600 and 44,000 fps).
  - Five outbound dots + 9 return swingbys = 14. So Fig. 9 supports "9 of 14" (p.154), not "about 75 percent".
- **Fig. 10:** advantages (cuts entry speed, Venus inspection, "breaks monotony of return trip", general
  availability) and disadvantages ("trip times increased by 20% in some years", more complex profile).
- **Fig. 7 (1975 transit trajectory):** 15.7 Dec 1975 opposition, trip time 430 days, 10-day stopover.
- **Repetition:** the only "cycle" in Part 5 is the Earth-Mars opposition cycle ("two cycles of oppositions",
  1971-1999). Nothing is periodic.
- **Relation to the JSR note** (`sohn-1964-venus-swingby-mode-manned-mars-missions-jsr-1-5-565-...`, held):
  same 65,600 / 3300 km / 15,000 fps / 31 deg / 14 deg / 44,000 fps example and the same 14 opportunities.
  The TM version adds Fig. 9 (per-year bars and trip times), which the JSR digest found unreadable in its own
  Figs. 2, 5, 6. **Real difference between versions:** the JSR text layer (`pdftotext`, line 323) has
  "increases the gross weight from 1.43 to 1.70 X 10^6"; this TM prints 1.77 on both p.153 and p.164
  (images). The JSR digest is right for the JSR. Retro to 60,000 fps (1.61) and to 50,000 fps (2.05) agree.
- Sohn's own reference: R. L. Sohn et al., "Feasibility Studies of Manned Mars Mission", STL Memo 9860.6-1,
  25 March 1963 (footnote 2, printed p.157).

## 3. Part 2, Ehricke (GD/A): Venus capture, powered fly-by and bi-planet missions

Read on images PDF pp.16-17, 27-28 (printed 13-14, 24-25).

- **Windows (printed p.13):** Earth-Venus short and medium transfers (120-150 d) from May to mid-July 1975 and
  November 1976 to end of February 1977; Venus-Earth returns from June 1975 (to November for 200-250 d legs)
  and from February 1977. Best Venus profile: short out, long return with aphelion beyond Earth's orbit.
- **Reference Mission I (Venus 1975, elliptic capture n = 8), Tab. 1-1 (printed p.25):**

  | | Dep. Ea | Transf. | Ar. Ve | Cpt. | Dep. Ve | Transf. | Arr. Ea |
  |---|---|---|---|---|---|---|---|
  | date | 5-20-75 | | 10-7-75 | | 10-27-75 | | 6-24-76 |
  | days | | 140 | | 20 | | 240 | |
  | v_inf* (EMOS) | 0.1327 | | 0.1400 | | 0.2455 | | 0.2824 |
  | principal dv (km/s) | 4.24 | | 1.5 | | 2 x 1.31 + "3 x 24" | | 14 (entry speed) |
  | principal dv (ft/s) | 12,900 | | 4,900 | | 2 x 4300 + 10,200 | | 46,000 |

  - EMOS = Earth mean orbital speed. The printed pairs give 1 EMOS = 29.7-29.9 km/s = 97,400-97,900 ft/s.
  - Capture ellipse: periapsis 1.3 Venus radii, period 19 h, apoapsis circular orbit period 49 h, 23 or 24
    revolutions in the 20-day stay (printed p.13).
  - Text (printed p.14) says the Venus arrival date "was held constant to 11-5-75"; the table prints 10-7-75.
- **Mars reference missions, Tab. 1-2 (printed p.24; 400 dpi image):**

  | Mission | II | III | IV | V | VI |
  |---|---|---|---|---|---|
  | target | Mars | Mars | Ma/Ve | Mars | Ve/Ma |
  | dep. Earth | 9-5-75 | as II | 10-15-75 | 8-31-75 | 1-27-77 |
  | first leg (d) / mode | 160 / CC | / CC/SE | 160 / CC | 150 / PFB | 150 / PFB |
  | capture (d), window (d) | 30, 20 | | 20, 0 | 0, - | 0, - |
  | second leg (d) / mode | - | - | 200 / PFB (Venus) | - | 200 / CC (Mars) |
  | capture (d), window (d) | - | - | 0, - | 0, - | 20, 0 |
  | to Earth (d) | 220 | | 200 | 250 | 220 |
  | mission period (d) | 440 | | 580 | 400 | 590 |
  | Earth dep. dv1* (EMOS) | .153 | | .164 | .166 | .1605 |
  | target arr. dv2* | .1735 | | .1163 | - | - |
  | PFB dv* (first planet) | - | | - | .0359 | .0171 |
  | target dep. dv3* | .200 | | .1725 | - | - |
  | target arr. dv4* | - | | - | - | .106 |
  | PFB dv* (Venus, return) | - | | .00505 | - | - |
  | target dep. dv5* | - | | - | - | .169 |
  | unbraked Earth entry (EMOS) | .74 | | .415 | .592 | .547 |
  | Earth retro to v_E = .512 (dv6*) | .236 | | 0 | .0825 | .036 |
  | mission velocity (EMOS / ft/s / km/s) | .7625 / 74,500 / 22.7 | | .4578 / 44,600 / 13.6 | .2844 / 27,800 / 8.47 | .3886 / 38,000 / 11.6 |

  - Mission IV: Mars capture, then a powered Venus fly-by on the way home. It removes the Earth retro burn
    (entry .415 EMOS = about 12.4 km/s, below the .512 limit) for +160 days over Mission II.
  - Mission V: Mars powered fly-by round trip. Mission VI: Venus PFB on the way to Mars (1977).
  - Text (printed p.35): Venus fly-by plus a small powered manoeuvre cuts vehicle weight by about
    1.2 x 10^6 lb (547 t), at +160 days.
- **Tab. 1-3 (printed p.25)** gives Missions IIB, IIF (perihelion brake 0.08 EMOS) and IVB with calendar dates.
  IVB: dep. Earth 10-15-75, arr. Mars 3-23-76, dep. Mars 4-12-76, Venus PFB 10-29-76 (0.00505 EMOS),
  arr. Earth 6-17-77; 580 d; 0.45785 EMOS = 44,700 ft/s = 13.65 km/s.
- **Tab. 1-4** compares single-planet Mars missions with Mars capture / Venus PFB missions (9-5-75 and 10-5-75
  departures). Several cells are blurred on the scan; I transcribe none.
- **Synodic missions (printed pp.18-19; images PDF 21-22):** stay until the next favourable return window. Long-transfer
  ("conjunction") missions: about one year each way, 7-8 months at Mars, 800-1000 d total, mission velocity
  0.3-0.4 EMOS (9-12 km/s), Earth entry about 0.389 EMOS (11.6 km/s). Short-transfer: 170-250 d legs,
  400-500 d stay. A "synodic base" for longer surface stays is named.
- **Bi-planet missions (printed p.23; image PDF 26):** because Venus and Mars differ a lot in angular rate, "it is
  practically always possible" to meet Venus on the return from Mars, or to go Venus first. Classes: bi-planet
  capture, capture/PFB, PFB/PFB. Ehricke: such missions, "rather than single-planet missions, will represent the
  actual future mode of interplanetary transportation" (needs nuclear pulse or gas-core engines).
- **Gate relevance (`#942` R1(c)):** these are one-shot E-M-V-E and E-V-M-E round trips with powered fly-bys.
  None repeats. No collision.

## 4. Other cycler-adjacent content

- **Part 3, Martin (Lockheed), printed pp.77-80 (PDF 80-83):** the first phase studied a stopover, a single-launch
  "two planet flyby" passing Mars and Venus, and single-planet flybys. Phase two focused on the 1974 Venus and
  1975 Mars flybys. Charts assume a light-side pass at "500 nm" closest approach, "1.15 for Venus and 1.20 for
  Mars" planet radii. Fig. 1 (Venus): departures Oct-Dec 1973, total trip about 330-370 d. Fig. 2 (Mars):
  departures 10 Sep-10 Oct 1975, total trip about 650-680 d. (Curve readings, about +-5 d.)
  Earth entry (printed p.111): 14.02 km/s (46,000 fps) Venus, 14.63 km/s (48,000 fps) Mars.
- **Part 6, Jones (NAA), printed pp.230-231 and Fig. 16 (PDF 241):** a 1970 launch with a Venus encounter on
  the way to Mars versus a favourable direct May 1971 launch. With encounter: 310 d to Mars, 30 d stay, 230 d to
  Earth, 570 d; dv3 = 9500 fps; V_e = 42,000 fps; 3.1 M lb in Earth orbit. No encounter: 195 + 30 + 265 d,
  printed "480 DAYS"; dv3 = 16,350 fps; V_e = 56,300 fps; 1.5 M lb. Arrows mark 1 Aug 1970 (JD 2440800) and
  1 May 1971 (JD 2441070). With aerobraking at both planets the totals are 21,500 fps (Aug 1970 with Venus)
  versus 28,300 fps (May 1971 direct). Lunar encounter gave at most 200 fps.
- **Part 8, Hornby (Ames), printed p.273 and Fig. 4 (PDF 270, 274):** quotes the STL 1975 result as a cut of
  "almost 23,000 fps from the optimum direct value of some 68,000 fps to ... about 45,000 fps". Hornby: the
  method "can be used on the outbound leg" and might cut propulsive dv; with aerobraking and Venus flyby "the
  differences between 'good' and 'bad' years may be slight".
- **Part 9, Austin (GD/FW), printed p.289:** "The true anomaly of Mars at opposition exhibits a cyclical
  variation with a period of approximately 16 years." Austin's conclusions (PDF p.341, text layer) repeat
  "the approximate 16-year cycle of repeated opposition locations". This is the Earth-Mars opposition cycle (7 synodic periods =
  14.9 yr), the same cycle as Sohn's "two cycles" and the JSR note's "15 yr". It is not a cycler.
- **Part 15, Ehricke (BPTSM):** mission classes A-F (Appendix A, PDF p.672): one-planet fly-by, bi-planet
  fly-by, one-planet capture, bi-planet capture/fly-by, bi-planet capture, tri-planet combinations; manoeuvre
  codes including PFB(#), FB(#), perihelion brake and aphelion brake. Activity "084 Heliocentric Orbit
  Installation, 10-100 persons" (Tab. 4-3). Fig. 4-5 case (1) is a Mars capture with return via Venus; case (5)
  a Saturn capture with an unpowered Jupiter fly-by. The text says Earth "no longer plays a preferred role" in a
  mature, heliocentrically oriented transport system. No trajectory is given for a heliocentric installation.
- **Part 17, panel:** M. Faget and others favour a manned flyby as the first planetary step (PDF pp.734-739,
  text layer). No swingby or repeat-orbit discussion.

## 5. Arithmetic checks (`checks_tmx53049.py` -> `checks_tmx53049.out`)

All cells below were re-read on 400-450 dpi crops (Tab. 1-1, 1-2, 1-3: PDF pp.27-28; NAA Fig. 16: PDF p.241).
They stand as printed ("6-17-77", "2-22-76", "195", "265", "480" are clear at 450 dpi). I keep the printed values.
- **Tab. 1-2, Mission VI:** the dv column sums to .4886 EMOS; printed total .3886 (and 38,000 ft/s, 11.6 km/s,
  which match .3886). One cell or the total is off by exactly 0.1000. Columns II, IV, V sum correctly.
- **Tab. 1-2, Mission II period:** legs 160 + 30 + 20 + 220 = 430 d; printed 440. IV, V, VI add up.
- **Tab. 1-3, IIB and IIF:** 9-5-75 to 2-22-76 is 170 calendar days (printed 160) and 2-22-76 to 3-13-76 is 20
  (printed 30). End-to-end (420 and 450 d) agree. The arrival date "2-22-76" is probably a misprint for 2-12-76.
- **Tab. 1-3, IVB:** 10-29-76 to 6-17-77 is 231 days (printed 200); end-to-end 611 d versus printed 580.
  "6-17-77" is probably a misprint for 5-17-77 (= 10-29-76 + 200 d).
- **Tab. 1-1:** 10-27-75 to 6-24-76 is 241 days (printed 240). Text gives Venus arrival 11-5-75, table 10-7-75.
- **Tab. 1-1 unit pairs:** 4.24 km/s = 13,911 ft/s, printed 12,900; 3.24 km/s = 10,630 ft/s, printed 10,200.
  12,900 ft/s equals v_inf* 0.1327 EMOS (12,965 ft/s), so the km/s cell may be the burn and the ft/s cell the
  v_inf, or one is a misprint. The Dep. Ve cell prints "3 x 24" where the text has 3.24 km/s. Other pairs agree.
- **Venus capture ellipse:** periapsis 1.3 R_V, apoapsis/periapsis 8, modern GM: period 20.4 h (printed 19);
  apoapsis circular orbit 48.4 h (printed 49). Close; their constants differ from mine.
- **Lockheed 500 n mi:** Venus (6052 + 926)/6052 = 1.153 (printed 1.15). Mars (3390 + 926)/3390 = 1.273,
  printed 1.20. The Mars figure does not match 500 n mi with any Mars radius near 3400 km.
- **NAA Fig. 16:** 195 + 30 + 265 = 490 d; printed 480. JD 2440800 = 1 Aug 1970 and JD 2441070 = 28 Apr 1971,
  consistent with the arrows (1 Aug 1970, 1 May 1971).
- **Sohn Part 5 internal differences (images p.153, p.154, p.164, p.165):**
  - 1975 direct entry: 66,500 fps (p.153), "66" kfps (p.154), 65,600 fps (pp.164-165).
  - 1975 swingby entry: 46,000 fps (p.154) versus 44,000 fps (p.165).
  - Propulsive braking at Mars: gross weight 1.43 to 7.15 M lb (p.153) versus 3.25 M lb (p.164).
  - Earth hold: 50 days (p.153) versus 45 days (p.164); "30-day" versus "45-day" hold assumption.
  - Coverage: "9 of 14" return-swingby opportunities (p.154) versus "about 75 percent" (p.166, about 10.5).
    Fig. 9's five outbound-swingby dots support 9 of 14.
  - p.164 prints "reduced to below 5,000 fps" where p.154 has 50,000 fps: a dropped zero.
  These are the author's own inconsistencies (summary written from a different). I quote the body text
  (pp.164-166) as primary.

## 6. Size policy and OCR test

- **`ocrmypdf --skip-text --optimize 2 -j 4` (ocrmypdf 17.13.0) -> `tmx53049-opt.pdf`: 54,170,094 bytes
  (54.2 MB), a saving of 11.4% (image optimisation 11.9%). Still over the 50 MB guideline.** Run time about 3 min.
  - 750 pages kept. Output is PDF/A-2b; the input XMP metadata (all empty or generic) was dropped.
  - Renders: PDF pp.27, 218 and 512 at 150 dpi grey from both files are pixel-identical (difference bbox None).
  - Text layer: kept. `pdftotext` word lists are identical on 630 of 750 pages; on 120 pages only the
    extraction order or token joins differ (for example "DynamicsAstronautics" split differently), not the OCR.
  - `-O3` (lossy JBIG2) would save more but changes the images; I did not run it. The saving from `-O2` is real
    but small. Either file the original (61 MB, over guideline) or the optimised copy (54 MB, still over);
    a further cut needs a lossy setting or splitting by part.
- **Re-OCR test (`--force-ocr`, tesseract via ocrmypdf 17.13.0), two 20-page excerpts made with `qpdf --pages`:**
  - Excerpts made with `qpdf 19640017065.pdf --pages 19640017065.pdf 153-172 -- excerpt20.pdf` and
    `... 15-34 -- excerptB.pdf`, then `ocrmypdf --force-ocr -j 4`. Excerpts and text dumps deleted after use.
  - Excerpt A, PDF pp.153-172 (Sohn, clean typeset text): dictionary-word fraction 0.813 (NTRS layer) versus
    0.814 (new). Key strings (65,600; 3300; 44,000; 1999; "75 percent") found equally; one more "44,000" hit.
  - Excerpt B, PDF pp.15-34 (Ehricke, tables): 0.829 versus 0.812. Both layers garble Tab. 1-2 cells
    (NTRS ".425", new "418", image .415). New layer drops ".3886".
  - **Verdict: a full re-OCR would not help.** The NTRS layer is as good as tesseract on typeset text, and both
    fail on tables. Tables need image reading in any case. (`ocrq.py`, `ocrq.out`.)

## 7. Citation mining (cycler-relevant references only)

The TM has very few reference lists. Only Part 8 (one item) and Part 13 (radiation and guidance) have one;
Part 5 has one footnote. Ehricke cites his own Vol. I Condensed Summary only.

Not held (checked with `ls cyclers_pdf/papers | grep -i` and `grep -i CORPUS_INDEX.md`):
1. Sohn, R. L., et al. (1963), "Feasibility Studies of Manned Mars Mission", STL Memo 9860.6-1, 25 March 1963
   (Part 5 footnote 2, printed p.157). The STL company study behind the swingby mode. Not on the wanted list.
   Not held. Low priority; probably never public.
2. Wong, T. J. & Anderson, J. L. (1964), "A Preliminary Study of Mission and Spacecraft Requirements for Manned
   Mars Orbiting and Landing", NASA TM X-54032 (Part 8 ref. 1, printed p.275; read on the image). Not held; not on
   the wanted list. Background only (Ames in-house study); not proposed.
3. Breakwell, J. V., Helgostam, L. F. & Krop, M. A. (1963), "Guidance Phenomena for Mars Missions", AAS Symposium
   on Exploration of Mars, Denver (Part 13 refs., printed p.576, text layer). Not held. Guidance, not trajectory
   design; not proposed. (Do not confuse with wanted-list row 48, Breakwell-Gillespie-Ross 1959/1961.)

Works this TM is cited for (for the filer):
- Gillespie & Ross 1967 (HELD, `gillespie-ross-1967-venus-swingby-mission-mode-manned-exploration-mars-jsr-4-2-170-...`)
  ref. 2 = Part 5 (Sohn) of this TM. The Gillespie-Ross digest lists it as "not held"; filing this TM closes it.
- TM X-53049 itself is not a row on the wanted list (`grep -i 53049` finds nothing there). The sibling MSFC
  report TM X-53055 (Wood et al. 1964) is row 65 and is not in this volume.
- The task brief mentions "Ehricke's study cited as TM-53049". In this TM, Ehricke has Part 2 (NAS8-5026,
  manned interplanetary missions) and Part 15 (NAS8-11084, BPTSM). I found no repo note that cites Ehricke with
  this number; the only repo citation of TM-53049 is the Gillespie-Ross one, and that one is Sohn. If an outside
  source cites "Ehricke, TM X-53049", it means Part 2 (the trajectory part) unless it is about programme models.

Absent: no citation of Minovitch (JPL TR 32-464 or any other), Ross 1963 (wanted row 5), Crocco 1956 (row 11),
Hollister 1963, Battin, or Breakwell-Gillespie-Ross (row 48).

*Filed as `cyclers_pdf/papers/nasa-msfc-1964-proceedings-symposium-manned-planetary-missions-1963-1964-status-nasa-tm-x-53049-ntrs-19640017065.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
