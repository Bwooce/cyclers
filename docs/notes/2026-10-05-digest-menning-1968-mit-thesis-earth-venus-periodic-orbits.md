# Digest: Menning 1968, "Free-Fall Periodic Orbits Connecting Earth and Venus", MIT MS thesis (#960)

M. D. Menning, "Free-Fall Periodic Orbits Connecting Earth and Venus", S.M. thesis, Department of
Aeronautics and Astronautics, MIT. Submitted 30 July 1968; supervisor W. M. Hollister.
- 72 pages. The upload was image-only (0 words, md5 4d18add1f960f158ea993a06370a4457).
- I OCRed it with `ocrmypdf --skip-text` (tesseract) and filed the result as
  `cyclers_pdf/papers/menning-1968-free-fall-periodic-orbits-connecting-earth-venus-mit-ms-thesis-ocr.pdf`
  (md5 926edd41a0e02bb7b4455f511f8d293e).
- The OCR is poor on the typewriter digit 4, which comes out as h, L or u, and it scrambles the
  double-spaced body text. Per the policy's table rule, every number below was read from page images or
  checked against them.
- No DOI (thesis). It is ref. [2] of Pisarevsky 2008, and the source of Hollister & Menning 1970 (JSR
  7(10):1193-1199, doi 10.2514/3.30134).

Page convention: the thesis's own page numbers are given, with the PDF page alongside. The Appendix pages
A-1 to A-17 are PDF pp.54-70.

## 0. Verdict

This thesis is the primary source of the 15 Earth-Venus periodic orbits that Hollister & Menning 1970
published as Table 3. Appendix A tables A-3 to A-17 hold the same numbers as JSR Table 3, row for row.
The thesis adds:
- the flyby dates of Hollister's three circular-coplanar orbits 1H, 2H and 3H (p.A-2).
- the orbit type key (p.A-1).
- the method in full: direct-return types, the symmetric-return Lambert modification, minimax turn-angle
  selection by "cone vectors", and the iteration history.
- the estimate of the total number of orbits (p.41-42).

As a `#942` positive control it settles one JSR print error (orbit 13) and shows that the JSR date print
errors were already in the thesis.

## 1. Method (READ; text pages cited)

- **Definition (p.1-2):** a periodic orbit "recurrently" flies by a planet sequence. The first and last
  flybys are at the same planet "with identical spacecraft velocities and absolute planet orientations", so
  the period is an integer multiple of the time for the planets to repeat their absolute orientation.
- **Model (p.3-4):** patched conics with zero-radius spheres. Ellipses run "from planet center to planet
  center" and time on the hyperbolas is neglected. The equal-|v_inf| condition at each flyby is reached by
  iterating the flyby dates. Flybys must rotate the v_inf vector without passing below the surface; passes
  beyond 1.1 planet radii are acceptable (p.6).
- **Direct-return types (ch. 2, pp.7-10):**
  - full-revolution return: one planet period. Its v_inf tips lie on a circle, a "double infinity". The
    half-revolution return is a special case.
  - symmetric return: coplanar, time of flight between one and two planet periods (Fig. 2.1).
- **Symmetric-return semi-major axis (ch. 3, pp.12-16):** a modified Lambert time-of-flight equation, Eq.
  (3.6), with one extra revolution. It has two roots, one of them the planet's own orbit, and branch logic
  over curves A, B and C (Fig. 3.1) is used to reach the transfer root.
- **Turn-angle selection (ch. 4, pp.18-29):**
  - For a full-revolution return, the free return direction is a "cone vector" on a right-circular cone
    (Fig. 4.2).
  - The rule: choose the cone vector(s) that minimise the largest turn angle. This is a minimax over the
    single or two consecutive full-revolution returns (secs. 4.21-4.22, Fig. 4.3).
  - "Approximately two seconds are required to compute the turn angles at 16 planet encounters" (p.21).
- **Iteration (ch. 5, pp.30-39):**
  - Initial approximations are Hollister's three circular-coplanar orbits 1H, 2H and 3H (p.30-31). "Since
    orbit 1H contains no symmetric returns, the solution for this orbit is rigorous in the patched conic
    sense" (p.31).
  - Steepest descent reduced the summed |delta v| to 0.1 EMOS, then Newton-Raphson to an assumed
    convergence of 0.005 EMOS (p.33). This worked for 1H and 2H (orbits 1 and 2).
  - 3H stalled at 0.11 EMOS, a "ravine" (p.34). Davidon and conjugate-gradient methods (Fletcher-Powell,
    Fletcher-Reeves) were tried (pp.34-36).
  - 3H was solved by first replacing symmetric returns with full-revolution returns and then restoring them
    one at a time (p.37), which gave orbits 3-8.
  - The run time per orbit was a few minutes on the 1968 MIT Computation Center (p.40).
- **Number of orbits (p.41-42):** "Each orbit contains five encounters at Earth and five encounters at
  Venus. At each encounter either one of two direct return orbits is available ... a minimum of 1024
  acceptable orbits may exist". Excluding orbits with six or more symmetric returns, which would pass below
  a surface, "reduces the total number periodic orbits with acceptable flybys to 648".
- **Variations (p.42):** half-revolution, full-revolution, half-revolution sequences at Venus, or
  reversing the full-revolution and symmetric order, give more orbits.
- **Repeat (p.43):** 16-year cycle. Earth and Venus repeat every 8 years, so two spacecraft fly each orbit,
  "30 spacecraft" for the 15 orbits.
- **Further study (p.43):** "While Mars is a particularly interesting planet for further study, its
  relatively low gravitational attraction ... make flybys above the planet surface extremely difficult."
  This is a stated reason Mars-hosted returns were not pursued, which bears on R1(a)/(c).

## 2. Appendix A (READ, page images)

**A-1 key (p.A-1)**, in the order the returns occur:

| Orbits | Earth | Venus |
|---|---|---|
| 1 and 1H | 5FR | 5TFR |
| 2 and 2H | 5FR | 5FRSY |
| 3 and 3H | 5SY | 5TFR |
| 4 | 2FR, SY, 2FR | 5TFR |
| 5 | FR, 2SY, 2FR | 5TFR |
| 6 | 2FR, 2SY, FR | 5TFR |
| 7 | FR, 3SY, FR | 5TFR |
| 8 | FR, 4SY | 5TFR |
| 9 | 5FR | 2TFR, FRSY, 2TFR |
| 10 | 5FR | 2TFR, 2FRSY, TFR |
| 11 | 5FR | FRSY, TFR, 2FRSY, TFR |
| 12 | 5FR | FRSY, TFR, 3FRSY |
| 13 | 2FR, SY, 2FR | FRSY, 4TFR |
| 14 | 2FR, SY, 2FR | FRSY, 2TFR, FRSY, TFR |
| 15 | 2FR, SY, FR, SY | FRSY, 2TFR, FRSY, TFR |

FR = full-revolution, SY = symmetric, FRSY = full-revolution then symmetric, TFR = two consecutive
full-revolution returns.

**A-2: Hollister's orbits 1H, 2H, 3H**, flyby dates (JD - 2440000), blocks of E E V V V, closing E:
- 1H: 441 806 971 1196 1421 / 1592 1957 2125 2350 2575 / 2797 3163 3316 3541 3765 / 3935 4300 4471 4696
  4921 / 5077 5442 5664 5889 6114 / 6285
- 2H: 417 782 914 1139 1470 / 1612 1977 2086 2311 2642 / 2763 3128 3253 3478 3809 / 3927 4293 4427 4642
  4953 / 5107 5472 5591 5816 6149 / 6261
- 3H: 352 852 970 1195 1420 / 1542 2042 2142 2367 2592 / 2697 3197 3297 3522 3747 / 3853 4353 4477 4702
  4927 / 5038 5538 5644 5869 6094 / 6196
- "Repeating after 16 years." Each spans 5844 d (COMPUTED). Only dates are given; there are no V_r, turn or
  Rmin values for these three.

**A-3 to A-17: orbits 1-15.** Columns: planet; date (JD - 2440000); relative velocity (EMOS); turn angle
(deg); minimum distance to planet centre (planet radii). 26 rows each.

Checked against the Hollister & Menning 1970 Table 3, using the clean text-layer copy filed 2026-10-05 as
corrected (see `docs/notes/2026-10-05-hollister-menning-1970-table3-recheck.md`):
- Orbits 1, 2, 3, 5, 10, 11 and 13 were checked visually against page images.
- The other orbits were checked by digit-whitelisted tesseract at 400 dpi: 1417 well-formed cells. Cells
  whose 4 was lost were excluded.
- The thesis equals the JSR table EXCEPT:
  - **Orbit 13, row 1 V_r: thesis 0.124 (p.A-15), JSR 0.129.** The thesis is right: row 2 has 0.124 with the
    same turn (74.6 deg) and Rmin (2.99).
  - **Orbit 3, row 1 Rmin: thesis 1.78 (p.A-5), JSR 2.78.** The thesis typo is fixed in JSR. The closing row
    prints 2.78, and the Tisserand-hyperbola check from 0.162 EMOS and 58.8 deg gives 2.78 (COMPUTED).
  - **Orbit 5, row 12 Rmin (3192, 0.171, 37.7 deg): thesis 5.08 (p.A-7), JSR 5.03.** The check gives
    about 5.05, so this one is unresolved.
- The JSR date print errors are already in the thesis, so the thesis cannot settle them:
  - orbit 2: 5715.
  - orbit 5: 5877.
  - orbit 6: 2585 and 2810.
  - orbit 8: 4870.
- Also identical in both: orbit 1 row 12 is planet E at 3163, a full-revolution return from 2798, and orbit
  6 row 3 is 995.

## 3. Gate relevance

- **`#942` R1(b), Earth-Venus:** the thesis is the original, inclined-elliptic, patched-conic computation of
  the 15 published Earth-Venus cyclers.
  - Earth and Venus both host returns.
  - Only the 5FR/5TFR/FRSY/SY families above, with five E and five V encounters per 16-year cycle, were
    computed.
  - The "1024 / 648" count is an estimate of the size of that same family set. It is not a published
    catalogue.
  - Anything outside these 15 itineraries is not in this thesis. Examples: other return types, Earth-Venus
    with Venus n-pi returns of the Russell type, and other periods.
- **R1(a)/(c), Mars:** the author declines Mars on the grounds of a weak turn (p.43). There is no Mars
  computation.

## 4. Positive controls (sourced)

- The 15 orbit tables, A-3 to A-17. Use the corrected JSR Table 3 with the orbit 13 row 1 fix above.
- The 1H-3H dates (A-2). 1H is "rigorous in the patched conic sense".
- The convergence criterion: summed |delta V_r| = 0.005 EMOS (p.33).

## 5. Citation mining (policy step 4; references pp.71-72, PDF pp.71-72)

Held:
- [11] Breakwell & Perko 1965: AIAA 65-689 is not held, but its matching results are in the held
  Breakwell & Perko 1974 and Perko papers.
- [8] Broucke 1964 JPL SPS: not held. The related Broucke 1968 TR 32-1168 is held.
- Everything else in the list is not held.

Not held, in priority order. DOIs are left as UNCONFIRMED where Crossref was not checked here.
1. Hollister, W. M. (1968), "Periodic Orbits for Interplanetary Flight", MIT Experimental Astronomy
   Laboratory report RE-36. No DOI (report). The JSR 1969 version is doi 10.2514/3.29664 (CONFIRMED in the
   batch 1 list).
2. Minovitch, M. A. (1964), "The Determination and Characteristics of Ballistic Interplanetary Trajectories
   Under the Influence of Multiple Planetary Attractions", JPL TR 32-464. No DOI (NTRS). It is the
   multi-flyby Earth-Mars-Venus round trip that the thesis cites.
3. Hickman, D. E. (1968), "Guidance Requirements for Periodic Orbits", S.M. thesis, MIT. No DOI. Guidance
   for these Earth-Venus orbits.
4. Gillespie, R. W. & Ross, S. (1966), "Venus Swingby Mission Mode ...", AIAA Paper 66-37. The JSR 1967
   version is JSR 4(2):170-175 (cited by VanderVeen 1969).
5. Crocco, G. A. (1956), "One Year Exploration Trip Earth-Mars-Venus-Earth", Proc. VII IAC, Rome. No DOI.
6. Ross, S. (1963), "A Systematic Approach to the Study of Nonstop Interplanetary Round Trips", AAS 9th
   Annual Meeting. No DOI.
7. Battin, R. H. (1959), "The Determination of Round-Trip Planetary Reconnaissance Trajectories", J.
   Aero/Space Sci. 26(9):545-567.
8. Lockheed 1962 study 3-17-62-1; STL 1964 Manned Mars Mission Study 8572-6009-RU-000; Broucke 1964 JPL
   SPS. Reports.
9. Optimisation (low priority): Davidon 1959, ANL-5990; Fletcher & Powell 1963; Fletcher & Reeves 1964.
