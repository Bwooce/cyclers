# Digest: Sturms 1967, "Trajectory Analysis of an Earth-Venus-Mercury Mission in 1973", JPL Technical Report 32-1062 (#960 batch 30; wanted row 61; #951 R5)

F. M. Sturms, Jr., Jet Propulsion Laboratory Technical Report 32-1062, 1 January 1967 (NASA contract NAS 7-100).
NTRS 19670008498. Cover shows accession N67-17827; a handwritten NASA CR number is on the cover
(it looks like 81636; I do not rely on it). 62 PDF pages: cover, title, abstract (p. x), 52 printed pages, 10 tables, 81 figures.
- The report number is **32-1062**. The NTRS text layer says "32-7062". The cover image and the running head on every page I
  rendered (cover, abstract p. x, pp. 5 and 21) read 32-1062. The title is "in 1973": the cover image says so. (The OCR layer once reads "1978"; that is an OCR error.)
- Original NTRS file: supplied file `ntrs_19670008498.pdf`, md5 fb99cbdcd47bdf3d869a467beba5fffa.
- OCR copy (filed): `ntrs_19670008498-ocr.pdf`, md5 90596b73ee97b823642972cbfdbec1c6.
- Filed as `cyclers_pdf/papers/sturms-1967-trajectory-analysis-earth-venus-mercury-mission-1973-jpl-tr-32-1062-ntrs-19670008498.pdf`.
- NOT HELD at the start of this batch (`ls cyclers_pdf/papers | grep -i sturms` and `grep -i sturms CORPUS_INDEX.md`: no hit). Wanted row 61.
- How I read it: the OCR text for the prose and the reference list; page images for every number below.
  - Image-checked: cover, abstract (p. x), printed p. 5 (Tables 1-2, Fig. 4) at 110 dpi; printed p. 21 (Table 3) at 110 dpi; printed p. 23 (Table 4, rotated page) at 200 dpi.
  - Not image-checked: the 81 figures (labels only, from the text layer), and Tables 5-10 (guidance covariances). The intro prose and the reference list come from the text layer. No number in this digest depends on those, except the guidance numbers, which are in the image-checked abstract.

## 0. Verdict

**This is a design study of ONE one-way Earth-Venus-Mercury mission, the 1973 opportunity (the Mariner 10 geometry). It contains no repetition-cycle analysis, no cycler and no free return. Its value for #951 R5 is (a) a sourced
survey statement of which 1970s E-V-Me opportunities exist, and (b) the pre-flight conic data for the 1973 mission.**
- **VanderVeen's claims about this report are CONFIRMED.** VanderVeen 1969 (held) cites it as ref. 4,
  "Trajectory Analysis of an Earth-Venus-Mercury Mission in 1973, JPL TR 32-1062, 1 Jan 1967". Date, number and subject match. He also says "Ref 3 (Sturms 1966) found no acceptable 1977 or 1978 trajectories below its launch-energy limit." Sturms 1967
  repeats that survey result: "two (1977 and 1978) result in no trajectories to Mercury within the constraints of the survey", where the constraint is C3 below 21 km^2/s^2. VanderVeen's 1978 mission has C3 = 30.2 km^2/s^2, above that limit. So there is no conflict (sec. 2).
- **VanderVeen's "E-V-Me opportunities recur every 1.6 yr (refs 2-6)" is only loosely supported here.** Sturms
  counts "six Venus opportunities in the 1970's" (about one per 1.6 yr), but never states a recurrence period or an 8-yr cycle. Sturms says nothing on a 13-yr or 8-yr cycle.
- **Catalogue (PROPOSAL, not an edit):** `mariner-10-venus-mercury` (data/catalogue.yaml, line 50698; the altitude is on line 50699) is "not a cycler"
  and stays so. Sturms Table 4 gives a Venus closest-approach altitude of 5796.5 km for the 31 Oct 1973 launch (conic); the catalogue row has 5794 km for the flown trajectory (flyby_altitudes). That is a useful independent pre-flight corroboration of the order of 0.05 per cent, but Sturms is a design, not the flown trajectory. The flown launch was 3 Nov 1973 05:45 UTC (catalogue line 50716), Venus 5 Feb 1974 (catalogue notes) and Mercury-I 29 Mar 1974 (lines 50749, 50764), i.e. a 146 d flight. Sturms's Oct 31 case is 146.02 d with Mercury arrival 26 Mar 1974. I take the 3-day launch offset as the cause. I did not verify this against the held Giberson-Cunningham 1975 paper or Dunne-Burgess SP-424.

## 1. Method (text layer)

- Method of the first step: patched conics (SPARC program, Ref. 4 = Joseph & Richard, EPD-406), 1-day grid in launch date and Venus arrival date, analytic mean elements for planet positions. Then integrated trajectories (SEARCH driving SPACE) for five launch dates 5 days apart.
- Three constraints define the "region of possible trajectories": (1) energy match at Venus (arrival speed equals departure speed on the Venus hyperbola), (2) Venus periapsis radius above the surface (Venus radius taken as 6200 km), (3) C3 < 21 km^2/s^2 (about the Atlas-Centaur limit).
- Types and classes: eight combinations of Type I/II Earth-Venus, Type I/II Venus-Mercury, Class I/II. Only Type I + Type I gives positive Venus altitude, in both classes. Class II has shorter regions, lower altitude and longer flight time, so Sturms drops it.
- The 1973 range of Venus arrival dates in the region is 2-7 Feb 1974.
- Midcourse: three corrections about 6 d after injection, 6 d before Venus and 8 d after Venus; about 120 m/s total; final rms miss at Mercury 1400-2900 km (abstract, image-checked).

## 2. Opportunity survey (prose) and the 1977-78 statement

- Sturms Sec. I.B: "Of the six Venus opportunities in the 1970's, two (1970 and 1973) result in favorable trajectories to Mercury; two (1972 and 1975) result in unfavorable trajectories to Mercury because of very low altitudes at Venus closest approach; and two (1977 and 1978) result in no trajectories to Mercury within the constraints of the survey." The source is the 1966 SPS 37-39 note (ref. 3, not held).
- My check of the six opportunities against the Earth-Venus synodic period (583.92 d, `check_cycles.out`): 18 Aug 1970 to 4 Nov 1973 is 1174 d = 2.01 synodic periods. Launch dates from VanderVeen's Table 2 (18 Aug 70, 1 Apr 72, 4 Nov 73) and his 10 Aug 78 are consistent with six launches in 1970-1978 spaced about 584 d, with an uneven spacing of 548-591 d from eccentricity. This confirms the "about every 1.6 yr" count, but it is my arithmetic.
- The 1970 and 1978 launches are 8 yr apart (Earth-Venus repeat) but Mercury is about 78 deg out of phase after 8 yr (33.22 Mercury revolutions). VanderVeen's results (1978 needs C3 = 30.2) are what that phase error would predict qualitatively. This explanation is mine, not in either report.
- Sturms says the 1973 mission has, relative to the 1970 design, "(1) higher launch energies, (2) larger altitudes at closest approach to Venus, (3) higher Venus arrival speeds, (4) smaller deflection angle at Venus, (5) smaller total flight times, (6) smaller communication distance at Mercury encounter."

## 3. Tables read on the page images

**Table 1, launch period against maximum launch energy, Venus arrival date free** (max C3 km^2/s^2: range of launch dates 1973, launch period days): 21.0: 12 Oct-21 Nov, 41 | 20.5: 16 Oct-20 Nov, 36 | 20.0: 19 Oct-18 Nov, 31 | 19.5: 20 Oct-17 Nov, 29 | 19.0: 25 Oct-15 Nov, 22 | 18.5: 31 Oct-8 Nov and 12-13 Nov, 11 | 18.15: 1-2 Nov, 2. The minimum C3 is about 18.15 km^2/s^2 (text).

**Table 2, same, Venus arrival date fixed.** 21.0: arrival 3 Feb, 17 Oct-10 Nov, 25 d; 4 Feb, 12 Oct-14 Nov, 34 d; 5 Feb, 20 Oct-17 Nov, 29 d | 20.5: 3 Feb, 22 Oct-6 Nov, 16 d; 4 Feb, 16 Oct-12 Nov, 28 d; 5 Feb, 20 Oct-15 Nov, 27 d | 20.0: 4 Feb, 20 Oct-9 Nov, 21 d; 5 Feb, 20 Oct-13 Nov, 25 d | 19.5: 5 Feb, 20 Oct-11 Nov, 23 d; 6 Feb, 31 Oct-14 Nov, 15 d | 19.0: 5 Feb, 25 Oct-7 Nov, 14 d; 6 Feb, 31 Oct-12 Nov, 13 d | 18.5: 6 Feb, 31 Oct-8 Nov, 9 d | 18.15: 6 Feb, 1-2 Nov, 2 d.
- Arithmetic check: all 5 + 12 launch-period counts equal the inclusive day count of the printed date ranges (`check_cycles.out`). The abstract's 23 d from 20 Oct to 11 Nov also checks.
- The selected design: Venus arrival 5 Feb 1974, C3 limit 19.5, launch 20 Oct-11 Nov 1973 (23 d); a second period at C3 limit 19.0 is 25 Oct-7 Nov (14 d).

**Table 3, launch window parameters** (parking orbit 90 nmi; launch azimuth sector 90-114 deg; coast 1072-1728 s; injection latitude about -21 to -29 deg). Five launch dates, three azimuths each. One row is checked: 26 Oct 1973, azimuth 90 deg: launch 05 35 34 GMT, coast 1704 s, injection 06 15 13. The gap is 2379 s, and 1704 s coast + 570 s ascent + 105 s burn (the report's assumptions 3 and 5) = 2379 s. It matches exactly. Not important for R5.

**Table 4, conic (SPARC) against integrated (SPACE), five launch dates** (read at 200 dpi on a rotated page; first value is the SPARC conic; SPACE values differ in the 3rd-4th digit):

| Launch 1973 | C3 km^2/s^2 | V_inf Venus km/s | Venus closest-approach altitude km | E-V TOF d | V-M TOF d | Total TOF d | Mercury arrival | V_inf Mercury km/s |
|---|---|---|---|---|---|---|---|---|
| 21 Oct | 19.3417 | 7.8544 | 4531.93 | 106.73 | 55.2348 | 161.96 | 1 Apr 1974 | 10.3924 |
| 26 Oct | 18.9385 | 8.0329 | 5481.01 | 101.74 | 51.9277 | 153.67 | 28 Mar 1974 | 10.5883 |
| 31 Oct | 18.7623 | 8.2211 | 5796.52 | 96.75 | 49.2729 | 146.02 | 26 Mar 1974 | 11.2869 |
| 5 Nov | 18.8733 | 8.4179 | 5657.68 | 91.76 | 46.9316 | 138.69 | 23 Mar 1974 | 12.3469 |
| 10 Nov | 19.3620 | 8.6230 | 5171.54 | 86.78 | 44.8850 | 131.66 | 21 Mar 1974 | 13.6346 |

- Venus encounter is fixed at 5 Feb 1974 0h in every SPARC column. Checks: E-V TOF plus V-M TOF equals total TOF in all five columns (to 0.005 d), and 5 Feb minus launch date gives 107, 102, 97, 92, 87 d, matching E-V TOF.
- Mercury arrival speed rises from 10.4 to 13.6 km/s as launch moves later, since the Venus-Mercury leg shortens (55 to 45 d). The lowest C3 (18.76) is at 31 Oct.
- The label "h_ca, km" is the Venus closest-approach altitude (Fig. 3 of the report is titled "Altitude of closest approach at Venus"). I take it as altitude, not radius, because the 6200 km Venus radius constraint would otherwise put it below the surface.
- The rotated page was read column by column. The Mercury arrival hours in the same table are not reproduced here.

## 4. Repetition structure, cycler, free return

- **None stated.** The report treats a single launch opportunity. Section I.B is the only place with an opportunity count.
- **Not a cycler, not a free return.** One E-V-Me one-way mission; the spacecraft ends at Mercury. Mercury encounter is a flyby, with an Earth-occultation or Sun-side aiming point (guidance Sec. V).
- Comparison, for R5 only: the Mariner 10 mission that followed (not discussed in this report) added a 2:1 Mercury resonance of 176 d (Mercury-I 29 Mar 1974, -II 21 Sep 1974, -III 16 Mar 1975, from the catalogue row). That is the only repetition structure in this lane. It comes from a post-Mercury-I Mercury flyby, not from the Earth-Venus-Mercury geometry.

## 5. Citation mining

| Ref | Item | Status |
|---|---|---|
| 1 | Minovitch, JPL TR 32-464, 31 Oct 1963 | not held; wanted row 1 |
| 2 | Sturms & Cutting, AIAA Paper 65-90, 25-27 Jan 1965 (1970 Mercury mission via Venus) | not held; row 61 (as JPL TM 312-505) |
| 3 | Sturms, SPS 37-39 Vol. IV, pp. 1-5, 30 Jun 1966 (E-V-Me opportunities in the 1970s) | not held; row 61. This is the survey behind the 1970-1978 table in Sec. I.B and the most relevant to the opportunity question. |
| 4 | Joseph & Richard, EPD-406, Jul 1966 (SPARC conic program) | not held; not on the wanted list, no cycler content |
| 5 | Clarke et al., JPL TR 32-77, 16 Jan 1963 (Design Parameters for Ballistic Interplanetary Trajectories, Part 1) | not held; not on the wanted list. Possible source for E-V opportunity contour plots (low priority). Also Clarke et al., TM 33-99 (Earth-Venus 1970, 1964), not held. |
| 6, 7, 8 | Kizner EP 674 (1959); Newhall TM 33-203 (1965); White et al. TM 33-198 (1965) | not held; tools, not relevant |
| 9, 11 | Sturms SPS 37-37 (Feb 1966) and 37-27 (Jun 1964), midcourse and error analysis of multi-planet trajectories | not held; not relevant to cyclers |
| 10 | Gates, JPL TR 32-504 (1963), midcourse execution errors | not held; not relevant |

Catalogue: no row change. One PROPOSAL: note in the `mariner-10-venus-mercury` source list that Sturms TR 32-1062 is the pre-flight conic design for the same Venus date. No new wanted rows except possibly adding AIAA 65-90 as an alias of row 61's TM 312-505 item.

*`check_cycles.out` is filed beside the Manning 1967 PDF as `cyclers_pdf/papers/manning-1967-minimal-energy-ballistic-trajectories-manned-unmanned-missions-mercury-nasa-tn-d-3900-ntrs-19670013962-check_cycles.out`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
