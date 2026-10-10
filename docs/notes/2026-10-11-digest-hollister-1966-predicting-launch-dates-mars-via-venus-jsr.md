# Digest: Hollister 1966, "Predicting Launch Dates for Mars Transfer via Venus" (#960 batch 43)

W. M. Hollister (MIT, Dept. of Aeronautics and Astronautics), J. Spacecraft and Rockets 3(6):925-927 (Engineering Note,
June 1966), doi 10.2514/3.28562 (Crossref checked 2026-10-11: title, volume 3, issue 6, pp. 925-927, June 1966).
Presented as AIAA Preprint 64-647 (AIAA/ION Astrodynamics Guidance and Control Conference, Los Angeles, 24-26 Aug 1964;
revision received 28 Sep 1965). The preprint, doi 10.2514/6.1964-647 (Crossref: record exists, 16 Aug 1964), is not held.
- Proposed filename: `cyclers_pdf/papers/hollister-1966-predicting-launch-dates-mars-transfer-via-venus-jsr-3-6-925-doi-10.2514-3.28562.pdf`.
  - File as given: 3 pp., md5 4aedfe92c28051b78e6d7fcbaf0137d8. Publisher typesetting with a text layer. No OCR needed.
  - Page 3 (p.927) carries only the reference list of this note (top left) and an unrelated note by Goldstein, Braun and
    Rose (thermionic diode). Nothing in it is by Hollister apart from refs 1-6.
- **How it was read:** pp.1-2 on page images at 110 dpi (all numbers below). The text layer was the second witness
  and agrees on every number. p.3 (references) was also checked on the page image: refs 1-6 read as printed in sec. 5.
  Figures 1-5 were read as images and are not digitised.
- Companion files in my folder: `hollister-1966-venus-swingby-trips.yaml` (the printed trip numbers),
  `hollister_checks.py/.out` (date and alignment checks), `lambert_check.py/.out` (patched-conic check of the printed v_inf
  values). All values in the scripts are ours.

## 0. Verdict

**A two-page qualitative note plus two worked 1970-72 itineraries. It has no equations and no periodic orbit.**
- The "prediction rule" is geometric. A Venus-assisted Mars trip needs the Earth-Venus, Venus-Mars and Mars-Earth
  alignments in that order, about 140 days apart (minimum energy) or about 200 days apart (longer trips that dip inside
  Venus's orbit). This repeats about every six years, "prior to every third opposition of Mars".
- **I reproduced the rule and all the printed trip numbers** (sec. 2, 3). The Sun-Venus-Mars alignment sits +73 degrees
  from the Earth-Sun line in March 1971 and -81 degrees in February 1972 (printed: about 70 ahead, about 70 behind).
- **R1(c) (`#942`): no collision.** Both trips are one-shot. The same itinerary is in the held 1963 thesis. Nothing is
  repeated except the launch opportunity. Same finding as the thesis digest and the Hollister-Prussing 1965 digest.
- **New against the 1963 thesis digest: the Mars arrival speed differs.** This note prints 0.343 EMOS. The thesis
  (Table 11.1, p.80, re-read on the page image: ".284") prints 0.284. My independent Lambert check of the printed dates gives
  0.341, so 0.343 reproduces and 0.284 does not. The thesis's own mass ratio does not settle it (sec. 3.3).
  **UNRESOLVED**, listed as an open item for the thesis digest. PROPOSAL only.
- Use: a sourced, checked statement of the E-V-M opportunity cycle (about 6.4 years) and a dated 1970-72 positive
  control for a patched-conic Venus-swingby code. Not a catalogue input.

## 1. Content (READ, pp.925-926)

**The rule (p.925, text).**
- Treat the trip as two single transfers, Earth to Venus then Venus to Mars. Two terminal dates fix each one.
- "Minimum energy transfers between planets are characterized by the fact that alignment of the launch and arrival planet
  occurs about midway through the transfer."
- For low-energy Earth-Venus, Venus-Mars and Mars-Earth legs, the three alignment dates "must occur in that order, with a
  spacing of about 140 days between them". This "is seen to occur about every six years or prior to every third
  opposition of Mars". Examples in Fig. 2: oppositions of 1965, 1971 and 1978. 1971 is also near Mars perihelion.
- When the orientation suits Earth-Venus-Mars before an opposition, it also suits Earth-Mars-Venus after it.
- Longer trips: let the path go inside Venus's orbit. The Earth-Venus leg is then not minimum energy. The spacing becomes
  "about 200 days". Outbound via Venus is then possible before the oppositions of 1967, 1973 and 1980. Return via Venus is
  possible after 1969, 1975 and 1982. "Every opposition period holds the possibility of a transfer via Venus."
- **Geometric form (Fig. 3, p.925-926):** in the frame where the Earth-Sun line is fixed, "Venus is continually getting
  ahead of Earth at a rate slightly greater than 1/2 degree per day. Mars is continually falling behind Earth at a rate
  slightly less than 1/2 degree per day." The Mars-Venus-Sun alignment direction shows the spacing of the other two
  alignments.
  - Min-energy Earth-Venus followed by min-energy Venus-Mars: alignment about **70 degrees ahead** of the Earth-Sun line.
  - Min-energy Mars-Venus followed by min-energy Venus-Earth: about **70 degrees behind**.
  - Longer trips: a range, shown by arrows, "somewhat arbitrary". Outside it the trips are too long, or too costly in
    launch energy.
- Fig. 3 labels the dated alignment directions: FEB 71, FEB 72, DEC 72, DEC 73, OCT 74, OCT 75, AUG 76, AUG 77, JUNE 78,
  JUNE 79, APR 80, APR 81 (read on the image; month and two-digit year only).
- Fig. 1a/1b give dials of Mars oppositions 1965-1978 and Venus conjunctions 1970-1977. Fig. 2 gives alignment dates
  1964-1980 on a Julian-date strip (ticks 2439000 to 2444000). Neither is tabulated in the text. Not digitised.

**The two worked trips (p.926, all numbers read on the image).** Units are EMOS, "approximately 18.5 miles per sec".

| Item | Trip 1 (Earth-Venus-Mars-Earth, 370 d, Fig. 4) | Trip 2 (Earth-Mars-Venus-Earth, 400 d, Fig. 5) |
|---|---|---|
| Earth departure | 1 Sep 1970, longitude 338 deg E of Aries, parking orbit 1.1 Re, burn 0.126, v_inf 0.126 | 18 Jun 1971 |
| Mars | arrival 19 May 1971, v_inf 0.343 (drag brake possible), 10-day stay | arrival 26 Sep 1971, 10-day stay, departure from 1.3 Mars radii, Mars longitude 353 deg, burn 0.141, v_inf 0.196 |
| Venus | 31 Dec 1970, v_inf in 0.174, passes 1.70 planet radii south, sunlit side, increment +0.022, v_inf out 0.211, turn 57 deg | 3 Apr 1972, v_inf in 0.244, 1.2 planet radii south, sunlit, increment -0.033, v_inf out 0.186 |
| Earth arrival | 6 Sep 1971 | 22 Jul 1972, v_inf 0.230, direct entry |
| Link to 1971 opposition | arrival 84 d before it | departure 56 d after it |

- A drag brake at Venus is "conceivable (but not recommended)".
- The Earth-Venus-Mars trip lets Mars be reached early enough for an economic return. The Mars-Venus-Earth trip allows a late
  Mars departure without a large Earth arrival speed.
- Ref. 6 (Hollister and Prussing, AIAA 65-700) "also presents computed data for transfers via Venus from 1970 to 1990. The
  data confirm the conclusions argued here."
- The note gives no planet constants, no ephemeris source and no flyby-geometry formula.

## 2. Checks of the date rule (OURS, `hollister_checks.py/.out`)

Mean J2000 longitudes with mean motions, circular-orbit alignments. Accuracy is tens of days, because Mars is eccentric.
1. **Leg lengths of the printed dates.** Trip 1: 121 + 139 + 100 d (after the 10-d stay), total 370 d. Printed 370.
   Trip 2: 100 d to Mars, 180 d from Mars (after the stay) to Venus, 110 d to Earth, total 400 d. Printed 400.
2. **Opposition link.** 19 May 1971 + 84 d = 11 Aug 1971. 6 Oct 1971 - 56 d = 11 Aug 1971. The two statements agree.
   Meeus gives 10 Aug 1971 for the opposition (external value, recalled, not checked here). One day apart.
3. **Synodic periods** (sidereal periods 365.2564, 224.7008, 686.9796 d): E-V 583.92 d, E-M 779.94 d, V-M 333.92 d.
   Three Earth-Mars synodic periods = 2339.8 d = 6.41 yr: "about every six years". Venus gains 0.6165 deg/d on Earth
   (printed "slightly greater than 1/2"). Mars falls back 0.4616 deg/d (printed "slightly less than 1/2").
4. **The 70-degree rule.** Each V-M conjunction steps by +205.9 deg relative to the Earth-Sun line per 333.9 d. Angle
   of the V-M alignment from the Earth-Sun line, mean model:
   - 9 Mar 1971: +73 deg. 6 Feb 1972: -81 deg. 1 Aug 1977: +73 deg. 1 Jul 1978: -81 deg. These are the 1971 and 1978
     windows. Printed: about +70 and about -70.
   - Spacing for the 1971 trip 1: Earth-Venus alignment (Venus inferior conjunction 11 Nov 1970) to V-M conjunction
     (9 Mar 1971) is 118 d. V-M to Mars opposition (14 Aug 1971 in the mean model, 10 Aug true) is about 158 d. The mean
     is about 138 d. Printed "about 140".
   - The midpoints of the printed legs fall on the alignments: trip 1 Earth-Venus midpoint 31 Oct 1970 (conjunction 11 Nov);
     Venus-Mars midpoint 10 Mar 1971 (conjunction 9 Mar).
5. **Window lists.** V-M conjunction before and after each opposition (mean model):
   - 1965, 1971, 1978: before by 156-161 d at +71 to +73 deg, after by 173-178 d at -81 to -83 deg. These match the 140-d
     favourable windows (Venus first, then Mars, then Earth).
   - 1967, 1973, 1980: before by 268-273 d at +123 to +125 deg. These are the printed outbound-via-Venus windows with the
     longer spacing.
   - 1969, 1975, 1982: after by 284-288 d at -131 to -133 deg. These are the printed return-via-Venus windows.
6. **Fig. 3 labels against the mean model.** Chain A (FEB 71, DEC 72, OCT 74, AUG 76, JUNE 78, APR 80) is exactly 22
   months per step. Chain B (FEB 72, DEC 73, OCT 75, AUG 77, JUNE 79, APR 81) is also 22 months per step. The model gives
   the 12 conjunctions 11 months apart, within -17 to +22 days of the labels, with alternating sign. INFERRED: the
   alternation is Mars's eccentricity (consecutive conjunctions are 175 deg apart on the Mars orbit). Labels are month
   precision, so I take the figure as consistent.

## 3. Checks of the worked trips (OURS, `lambert_check.py/.out`)

Standish approximate mean elements (1800-2050) for Venus, Earth and Mars, Sun-only Lambert, EMOS = 29.7847 km/s.

| Quantity | Printed | Ours |
|---|---|---|
| Trip 1 Earth v_inf | 0.126 | 0.127 |
| Trip 1 Venus in / out | 0.174 / 0.211 | 0.175 / 0.211 (Venus on 30 Dec) |
| Trip 1 Venus turn | 57 deg | 59.0 deg (Venus on 30 Dec; 60.8 on 31 Dec) |
| Trip 1 Mars arrival v_inf | 0.343 | 0.341 |
| Trip 1 Venus increment | +0.022 | +0.0216 from the printed v_inf pair, rp 1.70 Rv |
| Trip 2 Mars departure | 0.196 | 0.197 |
| Trip 2 Venus in / out | 0.244 / 0.186 | 0.243 / 0.184 |
| Trip 2 Earth arrival | 0.230 | 0.228 |
| Trip 2 Venus increment | -0.033 | -0.0325 from the printed pair, rp 1.2 Rv |

1. All printed speeds reproduce to 0.003 EMOS or better with an approximate ephemeris. The turn angle is 2 deg off
   (59 against 57), which is inside the sensitivity to the Venus date and the ephemeris.
2. **The Venus date.** The thesis digest found JD 2440950 (the thesis value) is 30 Dec 1970 (0 h on JD 2440950.5), while
   the text prints "31 Dec". This note also prints 31 Dec. The 30 Dec case fits the printed v_inf slightly better
   (Earth 0.127 against 0.131 on 31 Dec). INFERRED, not decisive: both dates are within the model error.
3. **Mars arrival: this note and the 1963 thesis differ.** Thesis Table 11.1 prints Mars v_inf 0.284 for the same date
   (19 May 1971). This note prints 0.343 and my Lambert check gives 0.341. Other cells differ slightly too:

   | Cell | 1963 thesis | 1966 note |
   |---|---|---|
   | Venus passage, trip 1 | 1.74 Rv, in 0.171, impulse 0.024 | 1.70 Rv, in 0.174, impulse 0.022 |
   | Venus in, trip 2 | 0.265, impulse 0.045 | 0.244, impulse 0.033 |
   | Earth arrival, trip 2 | 0.226 | 0.230 |
   | Mars arrival, trip 1 | 0.284 | 0.343 |

   Out v_inf, both Mars departure cells and all dates are the same in the two sources. The 1966 values are the later set
   (revision of Sep 1965). They reproduce. The thesis Venus-in value for trip 2 (0.265) and Mars arrival (0.284) do not.
   I re-read thesis Table 11.1 (PDF p.88, printed p.80) on the page image: it prints ".284" at Mars 19 May 1971 and
   ".171 - .211" at Venus. So the digest transcribed the thesis correctly. Two further internal tests of the thesis
   (ours, `mars_arrival_mr_check.out`, using the thesis eq. 9.1 as coded in the held `checks.py`, Mars surface speed 0.121 EMOS):
   - Eq. 9.1 mass ratio 1.9 (V_pi - 0.04)^(1/4): v_inf 0.284 gives 1.39-1.40, v_inf 0.343 gives 1.45. The printed 1.42 lies
     between them (it needs about 0.31). So the mass ratio does not discriminate.
   - The thesis entry limit for Mars (p.67, V_pi 0.34) corresponds to v_inf 0.294. A value of 0.284 sits just under that
     limit; 0.343 would sit above it. INFERRED: the thesis table may list a limit-respecting design value, not the
     Lambert v_inf of the dates. This is a guess.
   **Verdict: UNRESOLVED.** The Lambert check supports 0.343 for the stated dates and approximate ephemeris. The thesis
   internal checks do not decide. PROPOSAL: do not cite the thesis 0.284 as a Mars arrival speed without this caveat.

## 4. R1(c) history (the project's Earth-Venus cycler rows)

| Work | What it gives for R1(c) | Collision |
|---|---|---|
| Hollister 1963 Sc.D. thesis (held) | first Venus-swingby Mars trips; bi-elliptic idea; E-V-M every ~6 yr; 1970-72 itinerary | none (one-shot) |
| **Hollister 1964/66 (this note)** | the alignment rule and two worked 1970-72 trips; same itinerary as the thesis, later numbers | none (inherited from the thesis `...-evm-repeat-check.md`, which found no collision with ev-A, ev-B, ev-C or vm2-1; I did not re-test it) |
| Hollister and Prussing 1965 (held, AIAA 65-700) | Venus both ways 1971, outbound only 1973, homebound only 1975, both again 1978; thrusted flyby optimum | none |
| Hollister 1969 and Hollister-Menning 1970 (held) | the periodic Earth-Venus orbits (the periodic-orbit turn of the group's work) | not assessed here (digested separately) |
| Rall and Hollister 1971 (held) | periodic swing-by orbits Earth-Mars | not assessed here (digested separately) |

- This note is the bridge between the one-shot window studies (1963-66) and the periodic-orbit papers (1969-71). It
  states the geometry (the 6-year cycle and the alignment spacing) that the periodic orbits later repeat.
- Consistency with Hollister-Prussing 1965: that paper lists Venus available both ways in 1971 and 1978, outbound only in
  1973, homebound only in 1975. This note's windows (outbound 1967, 1973, 1980; return 1969, 1975, 1982; both ways
  1965, 1971, 1978) agree for the years that overlap.
- **Project use (PROPOSALS only):**
  - A 1970-72 patched-conic positive control for any Venus-swingby code: printed dates and v_inf values (yaml). Tolerance
    0.005 EMOS with an approximate ephemeris. The 0.343 Mars arrival is the cell that discriminates a correct code.
  - A test for any "next E-V-M opportunity" function: V-M conjunction at +73 deg on 9 Mar 1971 and 1 Aug 1977, 6.41-yr
    period (three E-M synodic periods), 333.9-d V-M synodic period.
  - Nothing for the catalogue.

## 5. Citation mining

References printed (p.927 top): [1] Ross, ed., Planetary Flight Handbook, NASA SP-35 (1963). [2] Carlson, N. A., Eng.
thesis, MIT (Feb 1965). [3] Hollister, Sc.D. thesis (1963). [4] Sohn, JSR 1:565-567 (1964). [5] Deerwester, AIAA 65-89
(Jan 1965). [6] Hollister and Prussing, AIAA 65-700 (Sep 1965).

Checked 2026-10-11 with `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX.md.
- Held: [3] Hollister 1963 thesis (`hollister-1963-mission-manned-expedition-mars-mit-scd-thesis-dspace-1721.1-15698.pdf`,
  digested 2026-10-07). [4] Sohn 1964 (`sohn-1964-venus-swingby-mode-manned-mars-missions-jsr-1-5-565-doi-10.2514-3.27702.pdf`).
  [6] Hollister and Prussing 1965 (`hollister-prussing-1965-optimum-transfer-mars-via-venus-aiaa-65-700-doi-10.2514-6.1965-700.pdf`,
  digested 2026-10-05).
- Not held: [1] Planetary Flight Handbook NASA SP-35 (wanted-list row exists as "Ross 1963" per the Hollister-Prussing
  digest). [2] Carlson 1965 MIT Engineer thesis (low-thrust mission geometry; no wanted-list row; low priority, it is
  low-thrust). [5] Deerwester AIAA 65-89, doi 10.2514/6.1965-89 CONFIRMED: wanted-list row 45.
- **Wanted-list row 41 is this file's conference form.** Row 41 reads "Hollister (1964), Predicting launch dates for Mars
  transfer via Venus, AIAA 64-647, doi 10.2514/6.1964-647 CONFIRMED, R1(c) background". With this journal form filed,
  PROPOSAL: mark row 41 "content held in another form (JSR 3(6):925-927, filed batch 43); acquire the preprint only for
  attribution". Same treatment as row 8 (Hollister and Prussing).
- New proposal: none. The note adds no new prior art.

*Filed as `cyclers_pdf/papers/hollister-1966-predicting-launch-dates-mars-transfer-via-venus-jsr-3-6-925-doi-10.2514-3.28562.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/hollister-1966-venus-swingby-trips.yaml`.*
