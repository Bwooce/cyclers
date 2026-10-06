# Digest: Manning 1967, "Minimal Energy Ballistic Trajectories for Manned and Unmanned Missions to Mercury", NASA TN D-3900 (#960 batch 30; wanted row 61; #951 R5)

L. A. Manning (NASA Headquarters, Mission Analysis Division, Moffett Field), NASA Technical Note D-3900,
April 1967. NTRS 19670013962. 33 PDF pages (cover, title page, summary, 20 pp. of text, references,
symbols and appendix, Tables I-IX at printed pp. 23-31, Figures 1-13 interleaved in the text).
- Original NTRS file: supplied file `ntrs_19670013962.pdf`, md5 e1b677815dd0dfb61a68a87551b47cd6.
- OCR copy (filed): `ntrs_19670013962-ocr.pdf`, md5 a8fc2868f35218fbf12f1edb77190671.
- Filed as `cyclers_pdf/papers/manning-1967-minimal-energy-ballistic-trajectories-manned-unmanned-missions-mercury-nasa-tn-d-3900-ntrs-19670013962.pdf`.
- NOT HELD at the start of this batch (`ls cyclers_pdf/papers | grep -i manning` and `grep -i manning CORPUS_INDEX.md`: no hit). Wanted row 61.
- How I read it: the OCR text for the prose; the page images for every number below.
  - Image-checked at 130-200 dpi: Tables I, II, III, IV, V (all printed pp. 23-26) and Table VI (p. 27, rotated; re-read upright at 200 dpi).
  - Not image-checked: Tables VII-IX (stopover with long stay, mode comparison), Figures 1-13 (labels and shapes only seen as text layer). No number in this digest depends on them, except the "6-1/2-year" and "6-7 year" phrases, which are prose.

## 0. Verdict

**This is a survey of one-way and round-trip Earth-Mercury missions for launches 1980-1999. It is not a cycler paper.
It gives the repetition structure of the direct Earth-Mercury geometry: a 13-year cycle of 4750 days, plus a
weaker 6-7 year cycle. Both check with synodic arithmetic (sec. 2).**
- **VanderVeen's claim about this report is CONFIRMED, with one limit.** VanderVeen 1969 (held) says Manning
  found "a 13-yr direct-Mercury cycle (41 x 116-d synodic periods)". Manning says: "The pertinent
  Earth-Mercury geometry repeats itself every 4750 days (13 Earth years)" (Results and Discussion, "Direct flybys"; repeated for the return legs: "repeats every 4750 days (13 years)"). 41 synodic periods of 115.88 d = 4751 d.
  - The limit: Manning states the cycle for the DIRECT mode only. He never writes "41 synodic periods". He
    never discusses an exact Earth-Venus-Mercury cycle. VanderVeen's "104 yr" (= 8 x 13) is his own step,
    not in Manning.
  - Manning does NOT show that the Venus-swingby results repeat every 13 years (sec. 3).
- **No cycler, no free return.** Every trajectory is one-way or a round trip with a stay at Mercury (Table VI: stays 63-89 d, total about 1 yr, total delta-V 21-24 km/s). Nothing repeats by itself.
- **Use for #951 R5:** a sourced statement that Earth-Mercury phasing repeats at 4750 d, and the first
  published tables of Venus-swingby-to-Mercury opportunities over a full 13-yr cycle (Tables III, IV).
  Both give opportunity dates only. Neither gives Venus-Mercury repetition (sec. 4).
- **Catalogue:** none. These are not cyclers. PROPOSAL: cite as the source of the 13-year Mercury cycle if the project ever writes a Venus-Mercury lane note.

## 1. What the report does

- Trajectories are patched two-body conics (planet motion from two-body equations, appendix). The data are
  "compiled from references 1-6".
- Missions: unmanned flyby (minimise Earth departure delta-V), unmanned orbiter (departure delta-V plus
  Mercury capture to a 1000 km orbit), manned stopover (outbound plus return). 1980-1999.
- Three modes: direct; unpowered Venus swingby; "modified pericenter" Venus swingby (a powered pass at
  250 km altitude, bounded below at 250 km, with a delta-V at Venus).
- Venus swingby needs an Earth-Venus leg first. Earth-Venus opportunities "occur every Earth-Venus conjunction (about every 17 months)" (Results and Discussion, unpowered swingby), "plotted in reference 6" = NASA SP-35. Manning then searches Venus-Mercury legs that match the Venus arrival speed to 15 m/s. He did not compute every trajectory.

## 2. The repetition arithmetic (my computation, `check_cycles.py` and `check_cycles.out` filed beside the PDF)

Mean sidereal periods: Earth 365.256 d, Venus 224.701 d, Mercury 87.969 d. Synodic periods: E-V 583.92 d, E-Me 115.88 d, V-Me 144.57 d.

| Interval | Earth phase error | Mercury phase error | Comment |
|---|---|---|---|
| 4750 d (Manning) = 13.005 yr | 0.0 rev (13.0046) | 53.996 rev, -1.4 deg | 40.99 E-Me synodic periods. This is Manning's cycle. |
| 13 Julian yr = 4748.25 d | 0 | 53.976 rev, -8.6 deg | Slightly worse than 4750 d. |
| 6.5 yr | 180 deg (6.5 rev) | 26.99 rev, -4.3 deg | Mercury returns to the same longitude, Earth is on the other side. Not a true repeat of a conjunction. |
| 7 yr | 0 | 29.06 rev, +23 deg | Mercury is 23 deg ahead. Weaker repeat. |
| 6 yr | 0 | 24.91 rev, -32 deg | Worse. |

- Manning's phrase "An apparent 6-7 year cycle also exists, but is much less exact than the 13-year period"
  is consistent with the 7-yr row (+23 deg). The "6-1/2-year" figure for the delta-V pattern is consistent
  with the 6.5-yr row: Mercury lands at about the same place, and a conjunction still occurs because
  20.49 E-Me synodic periods fit in 6.5 yr, but the phase flips by half a synodic period. I did not
  verify this physical reading against any figure.
- **Table VI repeats to the day.** Earth departure dates for 1980-1986 and 1993-1999 differ by exactly 4750 d in all seven pairs I checked (4366->9116, 4716->9466, 5109->9859, 5468->10218, 5828->10578, 6196->10946, 6562->11312). The delta-V and trip columns repeat to 0.1 km/s or less for the pairs I compared (e.g. 1980 total 22.0 vs 1993 21.9; 1981 23.5 vs 1994 23.4; 1982 23.4 vs 1995 23.5; 1983 22.8 vs 1996 23.0). So the 4750 d cycle is real in his own data. The tables for direct mode contain exactly one 13-year cycle (1980-1992 in Tables I-II, 1980-1999 in Table VI), so Tables I and II cannot show the repeat themselves; the later years of Table VI do.
- The table dates in this report are Julian Date minus 2 440 000; the footnote says 2440000 = 23 May 1968. My conversion script gives 1968-05-23 (`jd_dates.py`).

## 3. Tables (read on page images)

**Table I, direct flyby, yearly minima** (Earth departure date as JD-2440000, delta-V_Earth km/s, trip days, HBEV at Mercury km/s):
1980 4570 5.0 105 15.2 | 1981 4912 5.0 114 13.8 | 1982 5264 5.0 114 12.8 | 1983 5614 5.2 116 13.0 | 1984 5964 5.4 115 13.2 | 1985 6422 5.2 105 16.7 | 1986 6776 5.1 100 15.9 | 1987 7116 5.0 110 15.1 | 1988 7464 5.0 114 12.7 | 1989 7814 5.1 114 13.0 | 1990 8166 5.2 114 13.1 | 1991 8512 5.5 116 13.2 | 1992 8974 5.1 102 16.7.
- Departure-to-departure spacing in days: 342 352 350 350 458 354 340 348 350 352 346 462. The two long gaps (458, 462) are where the yearly minimum jumps to a different conjunction. The pattern HBEV 16.7 (1985) and 16.7 (1992) is a 7-year repeat in the table, not 6.5.
- Manning's lower bound for direct transfer: 4.8 km/s, 115 d (Hohmann to Mercury aphelion radius in the ecliptic plane).

**Table II, direct orbiter, three local minima per year.** Columns: departure JD-2440000, trip days, delta-V_total, delta-V_Earth, delta-V_Mercury km/s. Lowest total in each year, read on the image:
1980 4366 90 13.4 6.8 6.5 | 1981 4716 85 14.5 6.9 7.6 | 1982 5270 110 15.4 5.1 10.3 | 1983 5468 140 15.4 8.9 6.5 | 1984 5828 130 14.4 8.1 6.3 | 1985 6196 110 13.3 7.0 6.3 | 1986 6562 95 12.9 6.7 6.3 | 1987 6918 85 13.9 6.9 6.9 | 1988 7264 85 15.0 6.8 8.2 (ties with 7470, 110 d, 15.0) | 1989 7815 115 15.7 5.1 10.6 | 1990 8022 135 14.9 8.6 6.4 | 1991 8386 120 13.9 7.6 6.3 | 1992 8752 105 13.0 6.8 6.1.
- The lowest orbiter total is 12.9 km/s (1986). The yearly minimum is lowest in 1986 (12.9) and 1992 (13.0), 6 yr apart in the table, and highest in 1982-83 and 1989 (15.4-15.7), 6-7 yr apart. This matches Manning's "about 6-1/2-year" variation (Fig. 2). I did not read Fig. 2 itself.
- The other two opportunities in a year are usually far worse (up to 18.2 km/s).

**Table III, unpowered Venus swingby (flyby and orbiter)**. Read at 200 dpi, sums checked (all 13 rows pass: trip = Mercury arrival minus Earth departure; delta-V_total = delta-V_Earth + delta-V_Mercury).
Conjunction year / Earth departure / delta-V_E / Venus passage / pericenter altitude km / Mercury arrival / delta-V_Me / total trip d / total delta-V (JD-2440000):
1980 4370 5.3 4490 340 4594 7.6 224 12.9 | 1982 5000 4.2 5168 1980 5305 8.6 305 12.8 | 1983 5484 4.8 5569 270 5628 9.4 144 14.2 | 1985 6210 4.7 6361 1150 6480 10.2 270 14.9 | 1986 6634 6.1 6722 250 6770 7.6 136 13.7 | 1988 7160 8.2 7423 320 7595 10.0 435 18.2 | 1990 7680 6.3 7865 1190 7922 8.3 242 14.6 | 1991 8412 5.4 8556 300 8650 10.0 238 15.4 | 1993 8840 6.7 9026 2960 9075 10.0 235 16.7 | 1994 9650 4.2 9824 1360 9940 7.6 290 11.8 | 1996 (JD 2450000+) 0005 7.4 0192 3420 0235 9.1 230 16.5 | 1998 0810 4.0 0976 3620 1100 7.9 290 11.9 | 1999 1210 7.8 1361 600 1400 7.5 190 15.3.
- Flyby note: the delta-V_Me column is the Mercury hyperbolic excess speed for flyby (9.4 km/s for 1980, 11.8 for 1985 per the text).
- Some years are skipped (no 1981, 1984, 1987, 1989, 1992, 1995, 1997). Those conjunctions have no low-energy swingby.

**Table IV, modified pericenter swingby** (Venus altitude fixed at 250 km; delta-V at Venus; sums checked, all 6 rows pass). Conjunction / departure / delta-V_E / Venus date / delta-V_V / Mercury arrival / delta-V_Me / trip d / total:
1980 4336 4.3 4412 0.9 4460 6.4 124 11.6 | 1985 6200 4.7 6363 0.1 6460 9.6 260 14.4 | 1988 7330 5.6 7519 0.9 7630 4.7 300 11.2 | 1991 8440 7.0 8657 1.5 8780 5.4 340 13.9 | 1996 0180 4.5 0305 3.0 0390 6.8 210 14.3 | 1999 1310 3.9 1443 3.2 1540 7.1 230 14.2.

**Calendar dates** (converted from the JD columns, `jd_dates.out`): Table III 1980: depart 10 May 1980, Venus 7 Sep 1980, Mercury 20 Dec 1980. Table III 1993: 5 Aug 1992, 7 Feb 1993, 28 Mar 1993. Table IV 1980: 6 Apr 1980, 21 Jun 1980, 8 Aug 1980.

**Table V, direct return** (three per year; Mercury departure date, trip, delta-V_Mercury, Earth entry speed). First group (1980): 4538 200 8.7 19.2 | 4634 100 6.2 14.5 | 4754 80 10.8 14.3. Last group: 8936 185 8.4 18.4 | 9036 80 6.3 15.6 | 9162 140 9.9 14.2. These are the return opportunities for a 13-yr cycle; the dates span 4538-9162, i.e. 4624 d.

**Table VI, manned stopover, stay under 1 Mercury year** (direct out, direct back). Selected rows read on the 200 dpi image (launch year, departure JD, delta-V_E, outbound trip, delta-V_arrival, stay d, delta-V_departure, return trip d, Earth entry km/s, total delta-V, mission length d):
1980 4366 6.8 90 6.5 82 8.7 200 19.2 22.0 372 | 1986 6562 6.7 95 6.3 83 8.6 200 19.3 21.5 376 | 1992 8752 6.8 105 6.1 79 8.4 185 18.4 21.4 369 | 1993 9116 6.8 90 6.4 82 8.7 200 19.2 21.9 372 | 1999 1312 6.6 95 6.3 81 8.5 195 19.0 21.4 371. Stay times are 63-89 d: about one Mercury year (88 d), because arrival and departure both occur near Mercury's ascending node.
- Rounding check (all 20 rows): delta-V and mission-length columns sum to within 0.1 km/s and 2 d, except 1986 where 95+83+200 = 378 against printed 376 (trip times are rounded to 5 d). Not an error of reading.
- Manning also lists (text) a longer stay of about two Mercury years (176 d) that lowers the total delta-V to about 19 km/s (Table VII, not image-checked).

## 4. Does anything here repeat as a cycler or free return?

- **Direct Earth-Mercury:** the geometry repeats at 4750 d. That is a repeat of launch opportunity dates and delta-V, not of a trajectory that continues by itself. No flyby of Earth or Mercury is used to reshape the orbit.
- **Venus swingby:** Manning gives opportunities only for conjunctions that work. He does not claim that the swingby results repeat at 13 or 8 years. In Table III the 1980 and 1993 rows are not 13 years apart in launch date (4370 vs 8840 = 4470 d), and total delta-V is 12.9 vs 16.7 km/s. So Table III does not show the 4750 d repeat. The reason is mechanical: the Earth-Venus phase does not close at 13 yr (Venus 21.13 revolutions, 47 deg off; sec. 5).
- **Round trip with stay:** a stopover, not a free return. Earth entry speeds of 13-22 km/s (Tables V-VI).
- **Mercury-Mercury resonance from the spacecraft side** (Mariner-10-style 176 d orbit, two Mercury years) is hinted only by Manning's "stay of about two Mercury years (176 days)". The report does not treat it as a cycler.

## 5. Venus and joint E-V-Me repeats (my arithmetic, same script)

- 8 yr (VanderVeen's E-V cycle): Venus 13.00 rev (+1.5 deg), Mercury 33.22 rev (+78 deg). So Earth and Venus repeat, Mercury does not. VanderVeen's "25-day Mercury error" is 78 deg x 88 d / 360 = 19 d by this arithmetic (he gives 25 d; the difference is the rounding in his "within six days" estimate).
- 13 yr: Venus 21.13 rev (+47 deg), Mercury 53.98 rev (-8 deg). So 13 yr repeats E-Me, not E-V.
- 104 yr (VanderVeen): Venus 169.05 rev (+19.6 deg), Mercury 431.82 rev (-66 deg). Not an exact repeat either; his "104 yr" counts only the commensurability of the two periods in days, not a three-body phase repeat. Best whole-year joint repeat under 120 yr: 72 yr (Venus +13.6 deg, Mercury -17.9 deg). I did not check whether 72 yr is meaningful.
- Venus-Mercury alone: best near-resonances (continued fraction of the period ratio 2.5543): 9 Venus orbits (5.54 yr) = 23 Mercury orbits, drift -4.0 deg, 14.0 V-Me synodic periods; 83 Venus orbits (51.1 yr) = 212 Mercury orbits, +2.7 deg; 92 Venus orbits (56.6 yr) = 235 Mercury orbits, -1.3 deg. The 9:23 relation is a short V-Me repeat, but Earth does not return to the same place after 5.54 yr, so it does not repeat the E-V-Me opportunity. Not stated in either report.

## 6. Citation mining

| Ref | Item | Status |
|---|---|---|
| 1 | Minovitch, JPL TR 32-464, 31 Oct 1963 (NASA CR-53033) | not held; wanted row 1. Manning gives the NASA CR number 53033, which the wanted-list row lacks. |
| 2 | Sturms & Cutting, AIAA Paper 65-90 (1970 Mercury mission via Venus) | not held; wanted row 61 (listed there as JPL TM 312-505, probably the JPL memorandum form of the same work). The title Manning and Sturms 1967 both cite is "Trajectory Analysis of a 1970 Mission to Mercury Via a Close Encounter with Venus". |
| 3 | Niehoff, IITRI Rep. T-12, May 1965 | HELD (`niehoff-1965-analysis-gravity-assisted-trajectories-ecliptic-plane-iitri-report-t-12-n65-34460-ntrs-19650024859.pdf`) |
| 4 | General Dynamics, Nuclear Pulse Space Vehicle Study, Vol. IV, NAS 8-11053, Feb 1966 | not held; not a cycler item |
| 5 | Hollister & Prussing, AIAA 65-700, Optimum Transfer to Mars via Venus | not held (the held Hollister items are `hollister-1969-periodic-orbits-interplanetary-flight-...` and `hollister-menning-1970-periodic-swingby-earth-venus-JSR-7-10.pdf`, different papers). Not on the wanted list as a row; low priority. |
| 6 | NASA SP-35, Space Flight Handbook Vol. 3, Planetary Flight Handbook Part 1 (1963) | not held; source of the Earth-Venus opportunity plots |

Catalogue: none (no cycler content). No new wanted rows needed; consider adding NASA CR-53033 to row 1.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
