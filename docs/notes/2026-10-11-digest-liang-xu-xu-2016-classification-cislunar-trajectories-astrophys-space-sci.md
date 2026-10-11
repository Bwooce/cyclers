# Digest: Liang, Xu and Xu 2016, "The classification of cislunar trajectories and its applications in the Earth-Moon system" (Astrophys Space Sci 361:5)

**Digested:** 2026-10-11, owner upload, read in full from the text layer (two extractions) with the pages that hold
the table, the figures with printed numbers and the application results rendered and read.
**Citation:** Yuying Liang, Ming Xu and Shijie Xu (School of Astronautics, Beihang University), Astrophysics and Space
Science 361:5 (2016), DOI 10.1007/s10509-015-2577-z. Received 19 May 2015, accepted 16 November 2015, published online
2 December 2015.
**Filed as** `cyclers_pdf/papers/liang-xu-xu-2016-classification-cislunar-trajectories-applications-earth-moon-astrophys-space-sci-361-5-doi-10.1007-s10509-015-2577-z.pdf`
(19 pages, text layer, md5 a19f655b0be5d7350f7ccf9e6c73fb45, 22 MB). Table transcription:
`data/sources/liang-xu-xu-2016-tables.yaml`.

## Verdict

This is a transfer-classification paper, not a periodic-orbit paper. It tabulates no p:q resonant periodic orbit with
initial conditions. The one resonant orbit it draws is a Casoliva class member shown only as a picture, with a Jacobi
constant and a time of flight in a figure label. Nothing in it belongs in the task #1047 literal table as a numeric
row. The classification scheme does not name any of the project's families A, B, D or the nine symmetric families.

## What the paper does

Planar Earth-Moon circular restricted problem, mu = 0.0121516, C1 = 3.1884 and C2 = 3.1722 printed.
- **Two classes of transfer.** Direct: C far below C1, flight 2 to 6 days (Apollo-like, 3.13 days at C 2.08).
  Low-energy "SMART-like": C just below C1 with the L1 neck barely open, two-body energy E with respect to the Moon
  negative at perilune, flight 70 to 120 days (citing Parker and Anderson 2013).
- **One special class.** "Resonance cyclers", drawn as a 1:2 example (see below).
- **Method.** A grid of initial conditions on the section x = x_L1 = 0.83691257 with xdot > 0, y0 in [-0.50, 0.50],
  C in [1.99, 3.19], velocity angle in [-pi/2, pi/2] (120 x 360 x 100 grid points by the printed steps), propagated
  backward to the first Earth perigee (section Sigma2) and forward to the first lunar perilune (section Sigma3).
  Limits: 28 days in the Moon's region, 84 days in the Earth's, perigee radius 0.35 or less, perilune radius 0.15 or less.
  Perigee and perilune radii are matched to 1e-5 (3 km).
- **New parameters.** E (two-body energy about the Moon at perilune), beta1 and beta2 (phase angles of perigee and
  perilune).
- **Stated laws.**
  - Minimum time of flight rises in steps with C (Table 1).
  - In the Earth's region, smaller C gives larger perigee radius. The paper says direct transfers have perigee below
    0.10 and low-energy ones above 0.15.
  - In the Moon's region some perilune phase angles are unreachable for small C. Successive perilune passes differ by
    about 11 days (about 2.5 in normalised time).
- **Applications.** A geosynchronous (GEO) deorbiting plan with low-thrust arcs, and a lunar relay satellite (LRSS) on a
  "periodic" orbit.

## Table 1 (transcribed; three witnesses agree)

Minimum time of flight in days, by Jacobi constant C. Same values in `data/sources/liang-xu-xu-2016-tables.yaml`.

| C | 3.18 | 3.08 | 2.98 | 2.88 | 2.78 | 2.68 | 2.58 | 2.48 | 2.38 | 2.28 | 2.18 | 2.08 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TOF (d) | 67.4244 | 41.8138 | 45.9003 | 30.8022 | 4.7950 | 4.3102 | 4.3084 | 4.2464 | 3.4042 | 3.3654 | 3.2545 | 3.1319 |

The column is not monotone between C 3.08 and 2.98. The page prints it so. It is the only table in the paper.

## Orbits with numbers that the paper prints

| Item | Printed numbers | p:q, sense, initial conditions |
|---|---|---|
| 1:2 resonance cycler (Fig. 2 red curve, Fig. 10 four samples, Fig. 13 red dots) | C = 2.25 and TOF = 14 days (Fig. 2 label); "periods about 14 days" (text, citing Casoliva 2010) | label 1:2; sense not printed; no initial conditions |
| LRSS "periodic orbit" (Fig. 18) | period about 77.6 d (29.53 d at the Moon, 48 d at the Earth); perilune radius 1946 km; perigee 6513.4 km, apogee 6710.8 km | not a resonant orbit; no C, no initial conditions, no p:q |
| GEO deorbiting free arc | perigee a 1.8317e5 km, e 0.6145; perilune a 2.5213e5 km, e 0.7801 | a transfer arc, not periodic |

Nothing else with numbers is a closed orbit. The paper never prints a state vector, a period in the rotating frame, or
a stability index for any orbit.

### The resonance label is reversed against Casoliva

The paper defines a p:q resonance as p T_M = q T_S/C, with a_S/C = (p/q)^(2/3). Casoliva 2010 Eq. 9 reads
q T_moon = p T_spacecraft, and the project follows Casoliva. The paper's "1:2" cycler has a period of about 14 days,
which is half the lunar period. So it carries two spacecraft revolutions per lunar revolution. That is the project's
and Casoliva's 2:1, not 1:2. This is INFERRED from the printed period and the printed formula, not from a figure
count. The paper's remark that "1:2n" cyclers can be found from Fig. 10 is then a statement about the project's 2n:1.

## Relevance to this project

**Task #1047 literal table.** The question is whether this paper holds an Earth-Moon p:q resonant periodic orbit that
the table should carry as a numeric row.
- **No numeric row.** The only resonant orbit is the 1:2 (project 2:1) cycler. Its numbers are C = 2.25 and a
  14-day time of flight, taken from figure labels and a text remark. There are no initial conditions, no sense and no
  period to the printed precision. The paper credits the orbit to Casoliva 2010, so it is not a new orbit.
- **Class-level note only.** If the table wants one line, it is: same p:q as the project's 2:1 pro family (by the
  inference above), C 2.25 printed, sense not printed. The project's 2:1 pro symmetric family spans C 2.382 to 2.443
  (literal table 2026-10-11). C = 2.25 is 0.13 below the range, so even a sense-agnostic row fails the match rule.
  The row would not collide with any project family.
- **The LRSS orbit is not a collision candidate.** It has no p:q, no C and no state. Its printed Earth perigee and
  apogee (altitude 113 to 311 km) contradict its own Fig. 18, where the Earth-region loops reach roughly 0.2 to 0.8
  of the Earth-Moon distance. Treat those two numbers as a misprint or a different definition. Do not use them.
- **Classification names.** The scheme has three labels: direct, low-energy SMART-like, resonance cycler. It does not
  name families A, B, D or any symmetric family, and it does not use winding numbers.
- **Other Liang papers.** It is the earliest of the three held Liang, Xu and Xu cislunar papers. The 2017 paper
  (digest 2026-10-06) builds the polygonal-like orbits, and the 2020 paper (digest 2026-10-05) builds the p:q cycler
  infrastructure. Neither digest records that this 2016 paper uses the reverse p:q convention, so a table that merges
  the three Liang sources must convert this one.
- **#875 and the literature gate.** The paper adds one published class-level record of resonant cislunar cyclers,
  with no numbers to compare. No candidate-novel row changes.

## Figures

- **Fig. 1 (page 4).** Earth-Moon rotating frame, x_r from -1.5 to 1.5. Five libration points, an L1 Lyapunov orbit
  with its stable (black) and unstable (blue) manifolds, and the green Hill-region boundary (labelled 3.17). The
  Earth's region, the neck and the Moon's region are labelled.
- **Fig. 2 (page 4).** The three classes in one frame. A dense blue SMART-like curve (C = 3.18, TOF over 70 days)
  inside the green Hill boundary. A black direct curve (C = 2.26, TOF 3.2 days) from near the Earth to the Moon. A red
  1:2 curve (C = 2.25, TOF 14 days): an Earth-centred inner loop plus a large outer loop that reaches x_r near -0.75
  and passes the Moon.
- **Fig. 3 (page 6).** Definition of beta1 (perigee phase) and beta2 (perilune phase) on a feasible transfer.
- **Fig. 4 (page 6).** Earth segments of all low-energy transfers: a blue annulus bounded inside by the red minimum
  parking orbit of radius 0.20.
- **Fig. 5 (page 8).** Six panels. Stable manifolds of an L1 Lyapunov orbit mapped to their perigees at C = 3.08 and
  C = 3.18, against the minimum parking orbit. At 3.08 some mapped points fall inside it. At 3.18 none do at the
  first perigee, and they do at the sixth. This is the paper's explanation of the stepped time of flight.
- **Figs. 6 and 7 (pages 9 and 10).** Time of flight against beta1 and beta2 for chosen perigee and perilune radii,
  coloured by C band (below 2.2, 2.2 to 2.6, 2.6 to 3.1, above 3.1). Four parallel lines of slope 1.34 in Fig. 6.
  Tick-like lines separated by unreachable beta2 intervals in Fig. 7.
- **Fig. 8 (page 11).** Trajectories ended at the first, second and third perilune, explaining the 11-day spacing.
- **Fig. 9 (page 12).** Unstable manifolds at C = 3.18 twisting around the Moon, with three transfer segments.
- **Fig. 10 (page 12).** Four 1:2 cyclers in the rotating frame: multi-loop figure-eight curves around the Earth and
  the Moon. No axis values beyond the frame.
- **Figs. 11 and 12 (pages 13 and 14).** Perilune radius against E for a set of C values, near and far side branches.
- **Figs. 13 and 14 (pages 15 and 16).** Maps of perilune and perigee positions, grey scale by C, with "mountain
  ridge" structure and regions B1 to B5. Red dots in Fig. 13 mark four 1:2 cyclers.
- **Figs. 15 and 16 (page 17).** Candidate perilune and perigee regions for the GEO deorbiting arc, and the resulting
  trajectory.
- **Fig. 17 (page 17).** Candidate perilune sector for the relay orbit.
- **Fig. 18 (page 18).** The relay orbit: a three-lobed curve around the Earth joined to a small loop at the Moon,
  period 77.6 days printed.

## Internal inconsistencies noted (the paper's own)

- Worked insertion example: the text says C = 2.13, the Fig. 12 caption says C = 3.13, both with E = -0.09.
- Relay orbit Earth perigee and apogee against Fig. 18, above.
- The 1:2 label against the printed period, above.

## Citation mining (references checked against CORPUS_INDEX and the wanted list)

Held: Casoliva et al. 2010 (JGCD 33(5)), Szebehely 1967, the Koon, Lo, Marsden and Ross book (held as the 2006
file; the paper cites it as 2007).
Not held and not wanted: Belbruno and Marsden 1997 (the 2017 digest already records it).
Not held, not on the wanted list before today:
- Xu, M., Wei, Y. and Xu, S. (2013), "On the construction of low-energy cislunar and translunar transfers based on the
  libration points", Astrophys Space Sci 348(1):65-88. Same group, Poincare-map transfer survey. Low priority.
- Parker, J. S. and Anderson, R. L. (2013), "Low-Energy Lunar Trajectory Design", Wiley (the book). The corpus holds
  Parker's 2007 thesis. Parker, Anderson and Peterson (2013), JGCD 36(5):1501-1511, a survey of ballistic transfers.
- Li, W., Huang, H. and Peng, F. (2015), Adv. Space Res. 56:273-280, transit-orbit classification by support vector
  machine. Ren, Y. and Shan, J. (2014), CNSNS 19(3):554-569, spatial transit orbits.
- Quarta and Mengali 2005, Betts and Erb 2003, Mingotti et al. 2011, Pontani and Teofilatto 2013, Ivashkin 2001,
  Macau 2000, Uesugi 1993, Sousa-Silva and Terra 2014, Topputo 2013: transfer-design background, not wanted.
Already on the wanted list: Xu and Liang, IAC-13-A3.P.57 (the "Cislunar Bus System"), cited here as 2013b. It is the
Liang-group cycler paper that precedes this one in use and is the best remaining source for a numeric row.
