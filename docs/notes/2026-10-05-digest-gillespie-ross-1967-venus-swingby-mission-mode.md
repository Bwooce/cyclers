# Digest: Gillespie & Ross 1967, "The Venus-Swingby Mission Mode and Its Role in the Manned Exploration of Mars" (#960)

R. W. Gillespie & S. Ross (NASA Office of Manned Space Flight), "Venus-Swingby Mission Mode and Its Role
in the Manned Exploration of Mars", J. Spacecraft and Rockets 4(2):170-175 (February 1967),
doi 10.2514/3.28830.
- DOI confirmed with `scripts/crossref_check.py` (journal and volume match), 2026-10-05. The conference
  version is AIAA Preprint 66-37 (3rd Aerospace Sciences Meeting, January 1966), doi 10.2514/6.1966-37.
- Filed as `cyclers_pdf/papers/gillespie-ross-1967-venus-swingby-mission-mode-manned-exploration-mars-jsr-4-2-170-doi-10.2514-3.28830.pdf`.
  6 pages, text layer, md5 8a874ef00a9c7631a00be44f987f7a40. Former wanted-list row 10.
- I read all pages from the text layer. Table 1 is scrambled in the text layer, so I read it on the
  page image (p.171).

## 0. Verdict

These are one-shot Earth-Mars round trips with ONE Venus swingby, outbound or homebound. Nothing is
periodic, and no trajectory is chained into a repeating path. So there is no literal collision for
`#942` R1(c), a Venus-Mars cycler.

The paper is the source of two things that later work reuses:
- the swingby type numbers #1-#7, one per Mars-Venus alignment in a syzygistic period, which VanderVeen
  1969 uses;
- the E-V-M commensurability: a 2338-day "syzygistic period" (6.4 yr), and 5 periods = a 32-yr cycle
  after which the three planets return to the same absolute positions (p.170).

## 1. Content (READ)

- **Commensurability (p.170):** "relative configurations among the three bodies ... repeat fairly closely
  every (approximately) 2338 days". My check with mean synodic periods:
  - V-M 333.9 d x 7 = 2337 d;
  - E-M 779.9 d x 3 = 2340 d;
  - E-V 583.9 d x 4 = 2336 d;
  - 5 x 2338 d = 11,690 d = 32.0 yr.
- **Symmetry principle (p.170, from Ross 1963):**
  - On 24 Aug 1987 (JD 2447032) the Sun and the three planets are on one line, 4 deg from the Mars line
    of apsides.
  - A swingby trajectory with dates D1..Dn relative to that date mirrors a trajectory flown in the
    opposite direction with dates -Dn..-D1, with departure and arrival speeds swapped.
  - So outbound results follow from homebound ones (coplanar assumption).
- **Types (pp.170-171, Fig. 2):**
  - Mars and Venus align 7 times per syzygistic period. Each alignment gives one group (#1-#7) of easy
    Mars-to-Venus segments, and likewise outbound.
  - #2, #4 and #6 have no Venus-to-Earth continuation. #7 is not competitive. #1 is marginal. #5 is of
    occasional interest. #3 is promising.
- **Table 1 (p.171, page image), 1981-1987:**

| Type | Syzygy date | Trip, d | Stopover, d | LV Earth | AR Mars | LV Mars | AR Earth (EMOS) |
|---|---|---|---|---|---|---|---|
| 1 homebound | 5 Apr 1981 | - | - | no engineering window | | | |
| 3 homebound | 31 Jan 1983 | 567 | 28 | 0.15 | 0.11 | 0.21 | 0.17 |
| 5 outbound | 24 Dec 1983 | 446 | 20 | 0.25 | 0.26 | 0.15 | 0.30 |
| 5 homebound | 27 Nov 1984 | 464 | 18 | 0.14 | 0.30 | 0.18 | 0.13 |
| 3 outbound | 23 Oct 1985 | 560 | 29 | 0.14 | 0.20 | 0.11 | 0.11 |
| 1 outbound | 24 Aug 1987 | - | - | no engineering window | | | |

  EMOS = Earth mean orbital speed (about 29.8 km/s). Contour charts (Figs. 3-12) are in units of 0.1 EMOS.
- **Study results (pp.173-174):**
  - No #2, #4 or #6 trips exist, and no #1 trip is practical. The exception is a small 1974 outbound
    area at 0.26 EMOS, which "repeats once during each 32-yr syzygistic cycle".
  - #3 trips last about 550 d and recur every 3.2 yr, alternately outbound and homebound, "regularly and
    for decades into the future". They arrive at Mars near conjunction.
  - #5 trips come in outbound/homebound pairs every 6.4 yr, last about 460 d, and arrive at Mars near
    opposition.
- **Repeat caveat (p.175):** because Mars is eccentric and the planets' longitudes differ between
  periods, "the events in one period are not exactly reproduced during the period following". The 6.4-yr
  period is a guide, not a true period.
- **Direct-mission background (p.174):** opposition-class and conjunction-class trips repeat every
  2.13-yr synodic period, with a 15- or 17-yr cycle for Mars eccentricity effects.

## 2. Gate relevance

- **`#942` R1(c), Venus-Mars:** there is no literal collision. The paper is the first published statement
  (to my knowledge, INFERRED) of the 2338-d / 32-yr E-V-M commensurability that any repeating Venus-Mars
  or E-V-M path would sit on.
  - Russell-Strange VenMar#45 has period 667.8 d = 2 V-M synodic periods (R-S digest, sec. 4). A
    real-ephemeris check of a V-M cycler should test closure over 7 V-M synodic periods (one syzygistic
    period) and 35 (one 32-yr cycle) (INFERRED).
  - The #3 "every 3.2 yr, alternately outbound and homebound" pattern is a repeating opportunity, not a
    repeating trajectory.
- **Cite with:** VanderVeen 1969 (the E-V-M-V-E triples built on these types) and Sohn 1964 (the
  swingby origin).

## 3. Positive controls

- Table 1 dates and v_inf, for a one-shot patched-conic E-V-M or M-V-E Lambert check. The dates are
  syzygy dates, not encounter dates; the encounter dates are only on the contour charts (Figs. 6, 8-12).
  Read them on the page images before use.

## 4. Citation mining (references 1-7, p.175)

Held:
- none of them. (Hollister 1963 is not held; Hollister 1969 and Hollister-Menning 1970 are.)

Not held:
1. Hollister, W. M. (1963), "The Mission for a Manned Expedition to Mars", Sc.D. thesis, MIT. Already
   on the wanted list (VanderVeen digest).
2. Sohn, R. L. (1964), "Summary of manned Mars mission study", NASA TM-53049, Pt. 5. No DOI.
3. Sohn, R. L. (1966), "Manned trips using Venus flyby modes", JSR 3:161-169. DOI not checked. The
   follow-up to Sohn 1964 (filed today).
4. Deerwester, J. M. (1965), "Initial mass savings associated with the Venus swingby mode of Mars round
   trips", AIAA Paper 65-89. DOI not checked.
5. Ross, S. (1963), "A systematic approach to the study of nonstop interplanetary round trips", Adv.
   Astronaut. Sci. 13. Already in the wanted list. The symmetry principle.
6. Ross, S. (ed.) (1963), Planetary Flight Handbook, NASA SP-35. Background.
