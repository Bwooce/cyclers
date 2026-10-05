# Digest: Rall & Hollister 1971, "Periodic Swing-By Orbits Connecting Earth and Mars" (#960)

C. S. Rall (Bellcomm) & W. M. Hollister (MIT), "Periodic Swing-By Orbits Connecting Earth and Mars",
J. Spacecraft and Rockets 8(10):1017-1020 (October 1971), doi 10.2514/3.59763.
- CONFIRMED with `scripts/crossref_check.py`; the DOI is also printed on every page.
- Presented as AIAA Paper 71-92 (9th Aerospace Sciences Meeting, January 1971; Crossref
  10.2514/6.1971-92, titled "Free-fall periodic orbits connecting earth and Mars").
- "Based on an October 1969 MIT Sc.D. thesis under W. M. Hollister", which is held as
  `cyclers_pdf/papers/rall-1969-free-fall-periodic-orbits-connecting-earth-mars-mit-scd-thesis-msl-te-34-ntrs-19700017824.pdf`.
- Filed as `cyclers_pdf/papers/rall-hollister-1971-periodic-swing-by-orbits-earth-mars-jsr-8-10-1017-doi-10.2514-3.59763.pdf`.
  4 pages, text layer, md5 895b9f9248c49da6709941f0bbcdbb6e. Page 1 opens with the tail of an unrelated
  lunar-escape paper.
- I read all pages from the text layer. Table 1 (p.1019) is scrambled there, so I read it on the page
  image.

## 0. Verdict

This is the journal summary of Rall's thesis. It publishes ballistic Earth-Mars periodic orbits (M4-1,
M5-1, M5-2, M6-1) built from two reciprocal Earth-Mars-Earth round trips (Ross 1963) and two series of
direct returns at Earth.
- **`#942` R1(a), Earth-Mars where Mars also gives gravity assists:** relevant prior art. Table 1 shows
  that the Mars swing-bys turn the v_inf by 2.3-4.3 deg in the circular-coplanar case, and up to
  11.1-13.6 deg in the eccentric-inclined case.
  - So Mars is a (weak) working body in these 1969-1971 orbits.
  - But the method avoids direct returns at Mars, because "the small mass of Mars means that much less
    change in velocity occurs" (p.1017).
  - Whether this changes the R1(a) PARTIAL verdict is for the gate owner. The thesis is fuller (see the
    2026-06-07 mining note and its 2026-10-05 addendum).
- **`#942` R1(c), Venus-Mars:** not addressed here. The thesis sec. 4.4 records the failed Mars-Venus
  search.
- **Planet model (for the H-M control):** p.1018 says "Several periodic orbits have been computed for a
  solar system model that included accurate values for the planets' eccentricity, relative inclination,
  period, and location of the ascending node". No values are given; the thesis listing (p.136) gives the
  1960 mean elements.

## 1. Content (READ)

- **Definition (p.1017):** a periodic swing-by orbit is a free-fall trajectory that "visits one or more
  planets and revisits these same planets repeatedly for an indefinite period of time".
- **Method (pp.1017-1018):**
  - patched conic, with the sphere of influence neglected;
  - direct returns: half-revolution, full-revolution and symmetric (in Ross's sense, symmetric about the
    line of apsides);
  - reciprocal round trips, with dates the negatives of each other relative to opposition;
  - "The numerical techniques used are basically those of Menning", extended to half-revolution returns
    and multi-revolution legs.
- **Labels (p.1018):** Mn-m. n is the number of Earth-Mars synodic periods before the circular-coplanar
  pattern repeats (4, 5 or 6); each orbit makes two round trips to Mars per pattern. There are n versions
  (phases) of each orbit.
- **Eccentric-inclined period (p.1019):** lcm of n and 15 synodic periods (32 yr).
  - So M5 orbits repeat in 15 synodic periods (32 yr).
  - M4-1 repeats in 60 synodic periods (128 yr), "seventy-five" independent dates.
- **Table 1 (p.1019, page image):** speeds in EMOS, passing distances in planet radii, turn angles in deg.
  - A = Earth next to the short transfers, B = Mars, C = Earth next to the long transfers.
  - Row a) = circular coplanar; rows b)-d) = average, highest and lowest in the eccentric-inclined case.

| Orbit, region | v_inf a / b / c / d | Pass dist. a / b / c / d | Turn a / b / c / d |
|---|---|---|---|
| M4-1 A | 0.257 / 0.260 / 0.270 / 0.250 | 1.54 / 1.40 / 1.63 / 1.21 | 48.3 / 51.0 / 57.8 / 44.6 |
| M4-1 B (Mars) | 0.314 / 0.324 / 0.405 / 0.245 | 3.77 / 3.04 / 7.63 / 1.12 | 4.3 / 6.1 / 11.2 / 2.9 |
| M4-1 C | 0.181 / 0.195 / 0.229 / 0.178 | 1.30 / 3.60 / 38.7 / 1.16 | 77.4 / 61.1 / 83.0 / 6.3 |
| M5-1 A | 0.249 / 0.276 / 0.371 / 0.188 | 1.78 / 1.65 / 2.52 / 1.07 | 46.0 / 45.7 / 80.7 / 37.2 |
| M5-1 B (Mars) | 0.316 / 0.322 / 0.390 / 0.252 | 7.30 / 2.59 / 6.86 / 1.00 | 2.3 / 7.8 / 11.1 / 3.3 |
| M5-1 C | 0.211 / 0.216 / 0.260 / 0.171 | 1.55 / 1.53 / 1.76 / 1.31 | 60.7 / 60.8 / 80.8 / 51.2 |
| M5-2 A1/A2 | 0.245 / 0.244 / 0.283 / 0.212 | 1.42/2.06, 1.32/2.06, 1.73/2.41, 1.02/1.81 | 54.1/42.7, 57.9/43.6, 67.5/54.8, 41.4/33.4 |
| M5-2 B (Mars) | 0.314 / 0.325 / 0.399 / 0.246 | 4.79 / 5.57 / 14.2 / 1.31 | 3.4 / 4.4 / 13.6 / 1.4 |
| M5-2 C | 0.183 / 0.195 / 0.230 / 0.172 | 1.37 / 2.64 / 23.5 / 1.11 | 74.8 / 62.3 / 86.2 / 9.9 |

  - M4-1 change in encounter dates from the circular to the eccentric case, A / B / C, in days:
    average -0.8 / -0.3 / -0.8; rms 17.7 / 28.2 / 5.9; mean absolute 15.2 / 25.3 / 4.7.
  - M5-1 Mars lowest passing distance is 1.00 radius, which grazes the surface.
- **Transportation system (p.1020):** n vehicles (4, 5 or 6), one per version, give a short transfer
  each way at every opposition. "As few as four spacecraft" (abstract).

## 2. Gate relevance

- R1(a): see the verdict. The record of Mars turn angles above is sourced (Table 1, page image).
- The H-M control (#942): this paper only says that "accurate values" of e, i, period and node were used
  in the eccentric case. It confirms that the Hollister group computed in an elliptic-inclined model; the
  values are in the thesis listing.

## 3. Positive controls

- Table 1 circular-coplanar rows (a) for M4-1, M5-1 and M5-2: v_inf, passing distance and turn angle at
  Earth and Mars. The encounter dates are in the thesis appendices (transcribed 2026-06-07).

## 4. Citation mining (references 1-5)

Held: [1] Hollister 1969 JSR; [2] Hollister & Menning 1970 JSR; [4] Rall 1969 thesis.

Not held:
1. [3] Ross, S. (1963), "A Systematic Approach to the Study of Non-stop Interplanetary Round Trips", AAS
   9th Annual Meeting. Already in the wanted list.
2. [5] Hickman, D. E. (1968), "Guidance Requirements for Periodic Orbits", S.M. thesis, MIT. Already in
   the wanted list (not on DSpace).
