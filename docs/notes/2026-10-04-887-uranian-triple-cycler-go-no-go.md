# #887 stage 1: Uranian triple cycler, requirements and GO / NO-GO (reading and arithmetic only)

**Date:** 2026-10-04. **Scope:** establish what the published triple-cycler constructions
require of three moons' orbital periods, then test every Uranian triple against that. No search,
no project code, no trajectory. All numbers below come from one standalone script
(Appendix); every claim about a paper is a quotation with a page, or "not stated".

**Sources and page conventions.** All four are filed in the private paper corpus as:
`liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf` (the author's DRAFT manuscript,
header "Draft, Journal of Guidance, Control, and Dynamics", 22 PDF pages, so "p.N" is the
draft-PDF page, not the printed pp.146-155);
`lynam-longuski-2011-laplace-resonant-triple-cyclers-jupiter-acta-astronautica-69-doi-10.1016-j.actaastro.2011.03.011.pdf`
(journal pages; PDF page 1 is p.158);
`hernandez-jones-jesick-2017-one-class-io-europa-ganymede-triple-cyclers-AAS-17-608.pdf` (PDF
pages; PDF page 2 is conference p.974);
`russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf`
(PDF page; journal page = PDF page + 142).

## Verdict in brief

1. Liang's method needs a near-commensurability of two synodic periods that share a moon, about
   7:4 for Callisto-Ganymede-Europa with a leftover of 0.7365 d per 49.4 d (1.5 percent), and it
   absorbs the leftover by alternating two double cyclers (each closing the hub moon and one
   partner exactly, while the third moon drifts) and by multi-revolution phasing arcs. Decisive
   quote (p.2): the earlier methods "are not applicable because the relative phases of the three
   moons cannot be reset to their initial values. The proposed strategy ... mitigates the phase
   mismatch by switching between two double-cyclers". It states no numeric tolerance.
2. Yardsticks: Liang's triple has eps 13.6 degrees and M 0.7365 d per 49.4 d, and a largest
   inertial shift W of 12.4 degrees at 49.8 d; Io-Europa-Ganymede has eps about 0 and an
   inertial shift of 5.2 degrees per 7.05 d. On the relative-phase condition (Liang Eq. 1-2) all
   ten Uranian triples reach Liang's grade, so the arithmetic discriminates weakly; the best
   are Miranda-Ariel-Umbriel (eps 0.19 degrees, M 0.005 d at a 6.4 d repeat), Miranda-Ariel-Oberon
   (2.0 degrees, 83.7 d) and Miranda-Umbriel-Oberon (2.3 degrees, 113.7 d), the same top three
   on the stricter inertial figure W. Miranda-Ariel-Umbriel is far better than Liang's triple,
   but a best of ten: a 2.3 percent per-triple chance coincidence is about 20 percent for some
   triple of the ten. The GO rests on the margin over the yardstick, not on significance.
3. Miranda cannot bend a trajectory (0.4 degrees at 2 km/s) and sits 4.2 degrees off the
   equatorial plane of the other four.
4. **GO** (conditional, a narrow and cheap first test) for Miranda-Ariel-Umbriel. **NO-GO** for
   Miranda-Ariel-Oberon and Miranda-Umbriel-Oberon (dominated by it; the Oberon-reaching orbits
   need at least about 3.2 km/s at the middle moon). The non-Miranda triples
   (Ariel-Umbriel-Titania, Ariel-Titania-Oberon in particular) are NOT excluded by the arithmetic:
   they are decided by one test on the existing Liang reproduction.

## 1. What each published method requires

### 1.1 Russell and Strange 2009 (two moons; for comparison)

- Two-moon repeat: "total combined flight time that is an integer multiple of the synodic
  period of the two bodies of interest" (PDF p.3, journal p.145); "The cycler period is the
  combined time of flight (TOF) of each leg and is constrained to be an integer multiple of the
  system synodic period" (PDF p.8, journal p.150).
- Brief synodic periods "lead to an abundance of cycler solutions with short repeat periods and
  frequent initiation opportunities" (PDF p.1, journal p.143).
- Mass: "Although the target body in the idealized model is assumed to be massless, a
  free-return cycler connecting two massive bodies is possible as long as the flyby altitudes at
  the target body are sufficiently high" (PDF p.4, journal p.146). So a very light moon is
  admissible as a passive target body in the ideal model; mass alone does not disqualify
  Miranda.
- Consequence (my inference, arithmetic): for two moons the relative phase returns after every
  synodic period by definition, so no commensurability of orbital periods is needed. A third
  moon adds a second relative phase, and a repeat now needs both to reset together. That is the
  first genuine requirement on three periods, and it is what the three-moon papers address.

### 1.2 Lynam and Longuski 2011 and Hernandez, Jones and Jesick 2017 (Laplace-resonant)

- Lynam and Longuski, p.158: "Because of the dynamics of the Laplace resonance in the Jupiter
  system, cycler trajectories that periodically return to three bodies are possible"; the
  Laplace resonance is "a 1:2:4 orbital resonance among three moons or planets" and "the only
  Laplace resonance in the Solar System" among Io, Europa and Ganymede (p.158); "The
  synchronicity of the Laplace resonance allows three-moon cyclers to have periods that are
  commensurate with the period of the Laplace resonance" (pp.158-159).
- Their constraint, p.160 (Eq. 14): "n(T_Lap) = T_Eu,Io + T_Io,Ga + T_Ga,Eu where n is a
  positive integer and T_Lap is the time required for Ganymede, Europa, and Io to complete one
  cycle of the Laplace resonance". Also p.160: "the Laplace resonance does not have exactly the
  same period as Ganymede's orbit (7.05 days vs 7.15 days). This lack of periodicity requires
  that the argument of perijove of the triple cycler's trajectory precess by about 5.2 degrees
  every Laplace resonance period ... the 5.2-degree offset must be corrected using the
  gravity-assist effect of all three flybys." And: "in order to repeat indefinitely, the
  triple-cyclers must have the same orbital periods and perijoves for every cycle."
- Hernandez et al., PDF p.2: "The three inner Galilean moons are special in that their orbit
  periods exhibit what is known as a Laplace resonance, that is (almost) perfect integer
  multiples 1:2:4. Therefore, the synodic period is one Ganymede orbit, two Europa orbits, or
  four Io orbits, equivalent to 7.05 days. After one synodic period, all the moons return to
  their initial relative configuration; however, an angular shift in their inertial location of
  5.2 degrees occurs." Abstract (PDF p.1): the resonance "allows for trajectories that
  periodically fly by the three bodies, and under idealized assumptions repeat indefinitely".
  PDF p.2: ballistic repeatability of the high-fidelity solution "in general, will only last for
  a few cycles".
- Checkable conditions these methods impose: three periods in the ratio 4:2:1, so that
  S(inner,mid) : S(mid,outer) is 1:2 exactly, the relative configuration returns every 7.05 d,
  and the residual inertial shift is only about 5 degrees per repeat, which three flybys per
  repeat can correct.
- Why this does not transfer (arithmetic, section 3): the Uranian moons have no such relation.
  The angle n1 - 3 n2 + 2 n3 is 3e-5 deg/day for Io-Europa-Ganymede (circulation period about
  3e4 years, effectively locked; period set B, section 2) and the same angle for
  Miranda-Ariel-Umbriel is -0.0785 to -0.0884 deg/day, a circulation period of 11 to 13 years.
  Miranda-Ariel-Umbriel has the same (2,1,3) structure as the Laplace triple, but it circulates,
  and the whole configuration is rotated by -160 degrees (not 5 degrees) after each basic 6.44 d
  repeat. The Laplace "reset" argument does not apply to any Uranian triple.

### 1.3 Liang et al. (Callisto-Ganymede-Europa) 2024/25: the key method

**Why the earlier methods fail for a near-resonance (p.2):** "However, as for the near-resonance
system, the previous methods in [9, 10] are not applicable because the relative phases of the
three moons cannot be reset to their initial values. The proposed strategy in this Note
mitigates the phase mismatch by switching between two double-cyclers, through which the
dependence on the dynamics of Laplace-resonance system is avoided."

**Repeat condition (p.3, Eq. 1-3).** Callisto-Ganymede and Ganymede-Io have "an approximately
16:3 integer ratio for their synodic periods". The condition used is on two synodic periods that
share the middle moon: S(C,G)/S(G,E) = 12.5232/7.0509 = 1.7761, close to 7/4. "There is a
mismatch in the imperfect resonance per resonance period: M = 4 S(C,G) - 7 S(G,E) = 50.0928 -
49.3563 = 0.7365 days ... It means that the angular alignment of the three Galilean moons will
repeat approximately every 50 days." Eq. 3: the near synodic period T_S satisfies
T_S/T_C approx 3, T_S/T_G approx 7, T_S/T_E approx 14, "This means that every time the three
Galilean moons are synodic, their inertial location will not change greatly." And the premise:
"if the mismatch of the imperfect resonance can be handled through gravity-assist, it is
possible to design triple cyclers using the near-resonances."

**Choice of repeat period.** It is read off the commensurability, not searched: the synodic
alignment recurs about every 50 d, "approximately equal to 4 C-G synodic periods or 11 C-E
synodic periods" (p.10; 11 = 4 + 7). Conditions on the spacecraft orbit (pp.6-7): (1) "the near
synodic period is an integer multiple of the orbit period of the spacecraft"; (2) "the orbit
period of the starting moon is an integer ratio to the spacecraft's orbit period"; (3) the
period must put the apojove beyond the outermost moon. Perijove is between the planet radius and
the innermost moon's orbit radius (p.7, Eq. 15). In the example the spacecraft period equals
Callisto's (p.11).

**Tolerance on near-commensurability.** Not stated as a number. What is stated: the mismatch is
0.7365 d per about 49.4 to 50.1 d (1.5 percent; my division); "In about 150 days, the phase
differences will be about 100 degrees" if nothing absorbs it (p.9), which my arithmetic
reproduces (3 repeats = 2.21 d = 113 degrees of Ganymede-Europa relative phase, Appendix
section 2); so I use the same conversion for every triple: degrees = 360 x mismatch days / that
pair's synodic period. The existing project operator (`analyze_near_resonance`, note
`2026-07-03-alternating-double-cycler-operator.md`) records a scale-free 5 percent default
tolerance and a preference for the lowest-order match; that is the project's choice, not
Liang's.

**How the leftover is absorbed (pp.9, 18-19).** "the triple cycler is divided into two double
cyclers, one double cycler is completed in one cycle period, which is about 50 days and then
another double cycler is switched to and completed in one cycle period. The repeat period of the
triple cycler is now 2 cycle periods, which is about 100 days." (p.9). Switching is done "by
mainly changing the moon encounter epochs of the spacecraft rather than mainly rotating the
arch line of the spacecraft orbit", via the spacecraft period, "maneuvering at the perijove or
apojove and the moon gravity-assist can both be adopted" (p.9). Transfers are multi-revolution
Lambert arcs: "the multi-revolution transfer can also be seen as a process to reshape the phase
of the spacecraft" (p.9); "the gravity-assist is mainly used to modify the semi-major axis of the
spacecraft's orbit rather than mainly shifting the orbit's argument of periapsis ... the
spacecraft ... could 'wait' or 'chase' the targeting body" (pp.18-19). Epochs come from an
optimiser (IMODE, p.13) seeded by Eq. 16; flyby epochs may be shifted by "one or several synodic
period(s)" (p.12); transfer times are varied by 2 to 6 days (p.13). Not stated: whether the
inertial repeat of Eq. 3 is a requirement or an observed property of the chosen 50 d period. The
evidence in his own results points to the latter: "the phases of the Callisto and Ganymede at
each flyby epoch are similar, while there are significant differences among the phases of the
Europa" (p.15). My arithmetic on his Table 1 mean motions agrees: after one C-G half cycle
(4 S_CG = 50.09 d) Callisto and Ganymede are back (+0.6 degrees each) while Europa has moved
+38.3 degrees; after one C-E half cycle (11 S_CE = 49.62 d) Callisto and Europa shift -9.5
degrees and Ganymede -23.1. So each half cycle closes the hub moon (Callisto) and one partner,
and the third moon drifts freely until its own half. The two halves differ in length by
|4 S_CG - 11 S_CE| = 0.472 d (0.95 percent), the hub-specific mismatch; Eq. 2's 0.7365 d is the
Ganymede-hub one.

**Sequence pattern (pp.9-10).** Switched-double-cycler sequence CGCEC, base moon Callisto
visited three times per 100 d repeat; initial conic guess structures 1-1-1 and 1-1-0 (p.8, results p.13).
The sequences share a base moon, so the two doubles need the same cycle duration (about 50 d).

**Numbers reported.** Mean motions (p.4, Table 1): Europa 1.7693, Ganymede 0.8782, Callisto
0.3765 rad/day. Idealised model: ballistic triple cyclers "with a maintenance time of more than
1000 days" (p.13). Excess speeds (Tables 3, 5, 7, pp.15-17): 4.49 to 6.99 km/s for the two
high-perijove cases and 7.64 to 12.02 km/s for the low-perijove case. Minimum flyby altitude
constraint 50 km (p.11). Cycle TOF varies by "not more than 0.35 days" in the idealised model
(p.18) and by "almost 1 day" in the ephemeris model (p.19), "still about 100 days". Ephemeris
model: 10 cycles, 25 Sep 2033 to 22 Jun 2036, maximum Delta-v 1.0383e-7 m/s, all altitudes above
100 km (p.19). Not guaranteed indefinite: "the cyclers are not guaranteed to repeat
indefinitely, but the duration is enough for the potential Jupiter missions" (p.20).

**Liang's tables use very little bending (my arithmetic, his Eq. 7).** The altitudes in Tables
3, 5, 7 are "calculated with the Delta-v provided by the gravity-assist ... and the gravitation
of Jupiter is not taken into consideration" (p.18), so they are only an indicator. Re-evaluated
with Eq. 7, the Europa flybys bend 0.5 to 2.1 degrees (Delta-v 102 to 173 m/s, 0.7 to 1.3 percent
of Europa's circular speed) and the Callisto flybys 0.4 to 1.1 degrees (40 to 123 m/s, 0.5 to 1.5
percent of Callisto's circular speed). Appendix section 7.

### 1.4 Requirements stated as checkable conditions

| # | Condition | Source | Applies to Uranus? |
|---|---|---|---|
| R1 | Exact 1:2:4 mean-motion relation, so S(a,b) = 2 S(b,c) exactly and the relative configuration returns every 7.05 d | Lynam and Longuski p.158-160; Hernandez PDF p.2 | No: no Uranian triple has it |
| R2 | Residual inertial shift per repeat small enough (about 5 degrees) for the flybys to correct by turning the apse line | Lynam and Longuski p.160; Hernandez PDF p.2 | Only if a triple has it (section 3) |
| R3 | A near-integer relation n_ab S(a,b) approx n_bc S(b,c), leftover M about 1.5 percent of the repeat (0.95 percent with Callisto as hub) (Liang), lowest order available | Liang p.3, p.10 | Testable on periods |
| R4 | Near-integer relation of the repeat to each moon's own period (3:7:14 for CGE): inertial positions recur | Liang p.3 Eq. 3 | Testable; Liang's p.15 results show the unvisited moon drifts (38 degrees per half cycle), so this is probably NOT required; not stated |
| R5 | Spacecraft period T_sc: repeat is an integer multiple of T_sc, T_sc in integer ratio to the start moon's period, apojove beyond the outermost moon, perijove inside the innermost | Liang pp.6-7 | Testable on periods and radii |
| R6 | Two hub-sharing double cyclers with equal cycle duration, alternated | Liang pp.9-10 | Needs two real doubles per triple |
| R7 | Moons circular and coplanar in the construction model; real ephemerides tolerated only over about 1000 d | Liang p.4, p.13, p.20; Hernandez PDF p.2 | Inclination matters (Miranda) |

## 2. Moon data used

Orbital periods are NOT recorded in `src/cyclerfinder/core/satellites.py`; it records semi-major
axis, GM, mean radius and a flyby altitude floor, from JPL SSD (URA111 physical, mean elements,
accessed 2026-06-14 per its comments), and a system GM for Uranus (5.7945564e6 km^3/s^2, DE440)
and Jupiter (1.26686534e8). Primary periods are Keplerian, P = 2 pi sqrt(a^3 / GM_system); this
is what the project's own six two-moon rows use (their `synodic_period_days` reproduce to 4
decimals for five of the six rows; the sixth is quoted as 5.99). Eccentricity and inclination are not in `satellites.py`; they come from Pergola 2007
Table 1 (citing Seidelmann 2006, referred to the Uranus equatorial plane; digest
`2026-10-03-digest-pergola-2007-iepc-305-uranus-moons-manifold-electric-propulsion.md`).

| Moon | a (km) | GM (km^3/s^2) | R (km) | project flyby floor (km) | P Keplerian (d) | e | i (deg) | v_circ (km/s) |
|---|---|---|---|---|---|---|---|---|
| Miranda | 129846 | 4.3 | 235.8 | 100 | 1.41351 | 0.0013 | 4.232 | 6.68 |
| Ariel | 190929 | 83.5 | 578.9 | 50 | 2.52037 | 0.0012 | 0.260 | 5.51 |
| Umbriel | 265986 | 85.1 | 584.7 | 50 | 4.14423 | 0.0039 | 0.205 | 4.67 |
| Titania | 436298 | 226.9 | 788.9 | 50 | 8.70624 | 0.0011 | 0.340 | 3.64 |
| Oberon | 583511 | 205.3 | 761.4 | 50 | 13.46572 | 0.0016 | 0.700 | 3.15 |

The `#693` note (JPL SSD, relative to the same plane) gives Miranda 4.4 degrees and differences of
at most 0.1 degrees among the other four; Pergola's Seidelmann values for those four span 0.205
to 0.700 degrees, so the two sources disagree at the 0.5 degree level among the large moons. Both
are quoted; neither changes the Miranda conclusion.

**Second period set (sensitivity only, NOT verified in this session).** Sidereal periods typed
from standard tables: Miranda 1.413479, Ariel 2.520379, Umbriel 4.144176, Titania 8.705867,
Oberon 13.463234 d; Io 1.769138, Europa 3.551181, Ganymede 7.154553 d. Keplerian periods differ
from these by up to 0.0025 d (Oberon). For every entry within T <= 120 d the quality figure W
(section 3) changes by at most 0.4 degrees between the two sets and the best T is identical; for
long-T entries W changes by up to 2.5 degrees (for example Miranda-Titania-Oberon at 714 d: 3.7
against 6.2 degrees). Resonance-level quantities need the better set: the Laplace angle of
Io-Europa-Ganymede is -0.026 deg/day with Keplerian periods and -0.00003 with set B, so
Keplerian periods are not good enough to judge a true resonance, only to rank near-commensurabilities
at the several-degree level.

**Yardstick data.** Liang's triple: his Table 1 mean motions (reproduces his own S(C,G) =
12.5232, S(G,E) = 7.0509, ratio 1.7761 and M = 0.7365 d to 3 decimals: 12.5238, 7.0510, 1.7762,
0.7379 d, the small difference from the 4-decimal mean motions). Io-Europa-Ganymede: Keplerian
periods from the project's a and Jupiter GM; control: Hernandez's 5.2 degree inertial shift is
reproduced as -5.1 to -5.3 degrees after one synodic period (Appendix section 3a).

## 3. Arithmetic

Definitions. S(x,y) = 1/|1/P_x - 1/P_y|. For a triple (a,b,c) ordered by period, the
relative-phase layer (Liang Eq. 1-2) uses integers n_ab, n_bc with repeat T = (n_ab + n_bc) S(a,c)
(the third pair resets exactly; the other two are off by +eps and -eps), leftover
M = n_ab S(a,b) - n_bc S(b,c) days and eps = 360 (T/S(a,b) - n_ab) degrees. The inertial layer
(Liang Eq. 3) is W(T) = max over the three moons of 360 |frac(T/P_moon)| degrees, the largest
inertial shift of any moon after T; relative errors are at most twice W. Liang's p.15 results
show his construction does not close all three moons at once (section 1.3), so eps (with M) is
the primary figure and W is a stricter secondary figure, reported because the Laplace papers
do need an inertial shift small enough for flybys to correct (R2). Repeat-time cap: T <= 120 d for the primary
comparison (Liang's triple repeats every about 100 d, two cycles), and T <= 730 d (2 years) as the
outer limit: Liang's ballistic lifetime is about 1000 d (p.13, p.19), so a repeat longer than
about 300 d gives fewer than 3 repeats, and the published Uranian tours last 424 to 700 d
(digest `2026-10-03-digest-liang-2026-ice-giant-trajectory-design-review.md`). Matches that appear
only at several hundred days are what chance produces (the null below).

### 3.1 Pairwise (ten moon pairs)

| Pair | P ratio | S (d) | | Pair | P ratio | S (d) |
|---|---|---|---|---|---|---|
| Miranda-Ariel | 1.7831 | 3.2186 | | Ariel-Titania | 3.4544 | 3.5473 |
| Ariel-Umbriel | 1.6443 | 6.4322 | | Ariel-Oberon | 5.3428 | 3.1007 |
| Umbriel-Titania | 2.1008 | 7.9089 | | Umbriel-Oberon | 3.2493 | 5.9867 |
| Titania-Oberon | 1.5467 | 24.6321 | | Miranda-Umbriel | 2.9319 | 2.1452 |
| Miranda-Titania | 6.1593 | 1.6875 | | Miranda-Oberon | 9.5264 | 1.5793 |

(The project's rows quote 5.99, 24.6321, 3.1007, 7.9089, 3.5473 and 6.4322 for the six
non-Miranda pairs that are in the catalogue; all agree.)

### 3.2 Ten triples and the two Jovian yardsticks

Best figures per T limit, Keplerian periods (full output: run the Appendix script).
eps columns: best relative-phase eps over repeats with T <= 60 d and T <= 120 d, with
(n_ab,n_bc), T and M/T. W columns: best W within T <= 120 d, <= 365 d, <= 730 d (T in days). The last
column normalises as Liang does, T <= 3.1 x the outer moon's period (his repeat is 2.98 Callisto
periods).

| Triple | eps, T <= 60 d (n_ab,n_bc; T; M; M/T) | eps, T <= 120 d | W <= 120 d | W <= 365 d | W <= 730 d | W, T <= 3.1 P_out |
|---|---|---|---|---|---|---|
| **Liang C-G-E (his Table 1)** | 13.6 deg (4,7); 49.6 d; 0.738 d; 1.49 % | 1.9 deg (9,16); 112.8 d; -0.103 d | 12.4 @ 49.8 | 12.4 @ 49.8 | 3.4 @ 501 | 12.4 @ 49.8 |
| **Io-Europa-Ganymede** | 0.06 deg (2,1); 7.1 d; 0.002 d; 0.03 % | same | 3.2 @ 7.1 | 3.2 @ 7.1 | 2.5 @ 508 | 3.2 @ 7.1 |
| Miranda-Ariel-Umbriel | 0.19 deg (2,1); 6.4 d; 0.005 d; 0.08 % | same | 4.3 @ 58.0 | 4.3 @ 58.0 | 4.3 @ 58.0 | 18.9 @ 12.6 |
| Miranda-Ariel-Oberon | 6.7 (1,1); 3.2 d; 0.118 d; 3.7 % | 2.0 (26,27); 83.7 d; -0.035 d | 7.5 @ 80.6 | 7.5 @ 80.6 | 5.2 @ 512 | 19.5 @ 12.7 |
| Miranda-Umbriel-Oberon | 4.4 (14,5); 30.0 d; 0.099 d; 0.33 % | 2.3 (53,19); 113.7 d; -0.053 d | 10.6 @ 53.8 | 10.6 @ 53.8 | 2.8 @ 485 | 21.1 @ 12.7 |
| Ariel-Umbriel-Titania | 13.2 (5,4); 31.9 d; 0.525 d; 1.64 % | 2.5 (16,13); 102.9 d; 0.099 d | 24.9 @ 95.6 | 8.3 @ 174 | 8.3 @ 174 | 35.2 @ 25.3 |
| Ariel-Titania-Oberon | 2.5 (7,1); 24.8 d; 0.199 d; 0.80 % | same | 26.4 @ 52.9 | 9.4 @ 270 | 8.5 @ 592 | 38.9 @ 25.5 |
| Miranda-Ariel-Titania | 3.6 (11,10); 35.4 d; -0.068 d; 0.19 % | same | 18.3 @ 35.3 | 9.0 @ 209 | 3.5 @ 383 | 29.4 @ 25.4 |
| Miranda-Titania-Oberon | 4.5 (29,2); 49.0 d; -0.327 d; 0.67 % | same | 26.3 @ 26.8 | 11.6 @ 122 | 3.7 @ 714 | 26.3 @ 26.8 |
| Miranda-Umbriel-Titania | 4.7 (11,3); 23.6 d; -0.130 d; 0.55 % | same | 12.6 @ 8.4 | 12.6 @ 8.4 | 7.7 @ 609 | 12.6 @ 8.4 |
| Ariel-Umbriel-Oberon | 12.9 (1,1); 6.2 d; 0.445 d; 7.2 % | 5.7 (13,14); 83.7 d; -0.196 d | 21.1 @ 12.7 | 16.6 @ 270 | 2.0 @ 431 | 21.1 @ 12.7 |
| Umbriel-Titania-Oberon | 10.0 (3,1); 23.9 d; -0.905 d; 3.8 % | same | 24.5 @ 95.2 | 14.6 @ 270 | 14.6 @ 270 | 42.2 @ 25.4 |

(Rows are in order of eps at T <= 120 d, the primary figure, then W.)

How to read it. The yardstick is two-sided: Liang's triple (eps 13.6 deg, M 0.7365 d per 49.4 d,
W 12.4 deg at about 50 d) is the least that has been made to work with the double-cycler
switching, and Io-Europa-Ganymede (eps 0.06, W 3.2 to 5.2 deg at 7 d) is the Laplace case.

- **Relative-phase layer (Liang Eq. 1-2), the stated condition.** All ten Uranian triples have
  eps at or below Liang's 13.6 degrees at some T <= 60 d (Ariel-Umbriel-Titania, 13.2, is
  marginal; Ariel-Umbriel-Oberon passes only with the trivial (1,1) repeat at 6.2 d and
  M/T = 7.2 percent, which is not a commensurability). At T <= 120 d only Umbriel-Titania-Oberon
  (10.0 degrees, M/T 3.8 percent) misses Liang's M/T of 1.5 percent by a clear margin. The
  chance level is high: for a random synodic ratio the best order-11 fit has a median error of
  11.7 degrees and 58 percent of random triples are at or below 13.6 degrees (own Monte Carlo,
  Appendix section 4). So this layer discriminates weakly and cannot by itself justify a GO.
  By eps at T <= 120 d the order is Miranda-Ariel-Umbriel 0.19, Miranda-Ariel-Oberon 2.0,
  Miranda-Umbriel-Oberon 2.3, Ariel-Umbriel-Titania 2.5, Ariel-Titania-Oberon 2.5 (this last at
  a 24.8 d repeat, half of Liang's 50 d), then 3.6 to 5.7 for four more.
- **Inertial layer W (stricter, secondary).** At T <= 120 d only Miranda-Ariel-Umbriel (4.3),
  Miranda-Ariel-Oberon (7.5), Miranda-Umbriel-Oberon (10.6) and Miranda-Umbriel-Titania (12.6)
  are at or below Liang's 12.4: the same top three as on eps. The best non-Miranda triple is
  Ariel-Umbriel-Oberon at 21.1 degrees, at 12.7 d, under one Oberon period, so not a
  multi-revolution cycle; the next three non-Miranda triples (Umbriel-Titania-Oberon,
  Ariel-Umbriel-Titania, Ariel-Titania-Oberon) are 24 to 26 degrees, twice Liang's, and reach 8
  to 9 degrees only at 174 to 270 d (3.5 to 5.4 times Liang's 50 d cycle). Because Liang's own
  construction lets the unvisited moon drift by 38 degrees per half cycle, a W of 24 to 26
  degrees for the non-Miranda triples is not a failure of Liang's condition; it is a failure of
  the stricter Laplace-style one.
- **Normalised to Liang's 3.1 outer periods** no Uranian triple beats Liang's W (the best, at
  12.6 degrees, is Miranda-Umbriel-Titania at a single Titania period). Miranda-Ariel-Umbriel
  reaches 4.3 degrees only at 14 Umbriel periods (58 d); in absolute time that is within
  Liang's 100 d scale. Which normalisation is right is a judgement; I rank on absolute time
  because mission time is the constraint, and say so.
- **Chance baseline (my own construction, not from any paper); a best-of-ten correction.**
  Random period triples log-uniform in each system's range, best W with T <= 60 d: Uranian
  range, P(W <= 4.4) = 2.3 percent, P(W <= 12.3) = 17.5 percent; Jovian range, P(W <= 12.3) =
  9.3 percent (Liang's own triple is a moderately rare 7 to 9 percent case). The relative-phase
  figure for Miranda-Ariel-Umbriel (eps 0.19) is a 0.9 percent case. These are per-triple
  probabilities; with ten triples tried, the chance that some Uranian triple does this well is
  about 1 - 0.977^10 = 21 percent (W) or 1 - 0.991^10 = 9 percent (eps). So
  Miranda-Ariel-Umbriel is unusual but not statistically conclusive, and any other
  Liang-grade result is expected by chance; a match first appearing after about 300 d is chance
  and is not counted.
- **Miranda-Ariel-Umbriel in detail.** 2 S(M,A) = 6.4373 d against S(A,U) = 6.4322 d: leftover
  0.005 d per 6.44 d (0.08 percent, about 19 times tighter than Liang's 1.49 percent; 0.19 degrees
  of phase against his 13.6). That gives the relative configuration back every 6.44 d, but the
  inertial rotation is -160 degrees per basic repeat, so the strict usable repeat is 9 of them,
  57.97 d: Miranda 41.01, Ariel 23.00, Umbriel 13.99 revolutions, inertial shifts +4.0, +0.2,
  -4.3 degrees. This is the Uranian counterpart of Liang's 3:7:14 and of the Laplace relation,
  but it is not locked (circulation period 11 to 13 years), so over the 58 d repeat it holds to
  the quoted figures and drifts slowly: relative leftover 0.19 degrees per 6.44 d is about 11
  degrees per year, against Liang's 13.6 degrees per 49.6 d, about 100 degrees per year.

## 4. Physical feasibility at a glance

All figures use the project's circular-coplanar constants; they are lower bounds on difficulty.

**Bending (Eq. 7 of Liang, minimum altitude = the project's floor).** Maximum bend in degrees and
maximum Delta-v (2 v_inf sin(delta/2)), as a percentage of the moon's circular speed:

| Moon | v_inf 1 km/s | 2 km/s | 3 km/s | 4 km/s |
|---|---|---|---|---|
| Miranda | 1.4 deg, 25 m/s (0.4 %) | 0.4 deg, 13 m/s (0.2 %) | 0.2 deg, 9 m/s | 0.1 deg, 6 m/s |
| Ariel | 13.5 deg, 234 m/s (4.3 %) | 3.7 deg, 129 m/s (2.3 %) | 1.7 deg, 87 m/s | 0.9 deg, 66 m/s |
| Umbriel | 13.6 deg, 236 m/s (5.1 %) | 3.7 deg, 130 m/s (2.8 %) | 1.7 deg, 88 m/s | 1.0 deg, 66 m/s |
| Titania | 24.6 deg, 426 m/s (11.7 %) | 7.3 deg, 253 m/s (7.0 %) | 3.3 deg, 175 m/s | 1.9 deg, 133 m/s |
| Oberon | 23.3 deg, 404 m/s (12.8 %) | 6.8 deg, 238 m/s (7.6 %) | 3.1 deg, 164 m/s | 1.8 deg, 125 m/s |

Compared with what Liang's flybys use (0.5 to 1.5 percent of circular speed, 0.4 to 2.1 degrees),
Ariel, Umbriel, Titania and Oberon offer 2.3 to 7.6 percent at 2 km/s and still 1.6 to 5.2
percent at 3 km/s; Miranda offers 0.4 percent at 1 km/s and less above, below the lower end of
Liang's use at every v_inf from 1 km/s up. The project's two-moon rows (0.9 to 2.2 km/s) sit well inside the capable range of the
four large moons.

**Spacecraft orbit (Liang conditions 1-3).** Checked with the same vis-viva construction (circular
coplanar); the positive control reproduces Liang's Table 3 first row (orbit with perijove 660988 km
and period equal to Callisto's: Callisto 5.673 km/s against 5.6730 printed, Ganymede 6.994
against 6.9919). Condition 2 reads "integer ratio", so spacecraft periods T_sc = (q/p) P_moon
with p, q <= 4 were allowed, with the repeat T equal to an integer number of T_sc (within 0.08),
perijove at or inside the inner moon and apojove at or beyond the outer one, the lowest-energy
orbit for each period. Lowest options by largest V-inf (km/s at the three moons in the order of
the triple name; full lists in the Appendix output):

| Triple (repeat T) | spacecraft period, revolutions | V-inf at the three moons |
|---|---|---|
| Miranda-Ariel-Umbriel (58 d) | 2/3 P_Umbriel = 2.763 d, 21 revs | 1.11, 2.03, 1.14 |
| | P_Ariel = 2.520 d, 23 revs | 2.00, 2.21, 1.03 |
| Miranda-Ariel-Oberon (80.6 d) | 1/2 P_Oberon = 6.733 d, 12 revs | 1.89, 3.32, 1.42 |
| Miranda-Umbriel-Oberon (53.8 d) | 1/2 P_Oberon = 6.733 d, 8 revs | 1.89, 3.24, 1.42 |
| Ariel-Umbriel-Titania (95.6 d) | 2 P_Ariel = 5.041 d, 19 revs | 1.80, 2.11, 0.92 |
| Ariel-Titania-Oberon (52.9 d) | 3 P_Ariel = 7.561 d, 7 revs | 1.28, 1.90, 1.09 |

Each of these keeps all three V-inf at or below the project's 2.2 km/s band except the two
Oberon-reaching Miranda triples, whose best options leave 3.2 to 3.3 km/s at Ariel or Umbriel
(the absolute floor for any orbit spanning Miranda and Oberon is 3.19 and 3.28 km/s at Umbriel
and Ariel, perijove at Miranda, apojove at Oberon). At 3.3 km/s Ariel and Umbriel bend about 1.4
degrees, and Miranda 0.2 degrees. Miranda-Ariel-Umbriel can stay at 1.0 to 2.2 km/s with two
good bending moons.

**Miranda's inclination.** 4.232 degrees (Pergola, Seidelmann 2006) or 4.4 degrees (JPL via the
#693 note). Out-of-plane extent a sin i = 9,580 km, 41 Miranda radii; velocity out of plane at a
node v_circ sin i = 0.49 km/s, which is half to a fifth of the V-inf the triple would use at
Miranda (1.0 to 2.2 km/s). A spacecraft in the equatorial plane of the other four moons reaches
Miranda within 2,000 km of its plane only within 12 degrees of a node, 13 percent of an orbit.
The basic 6.44 d repeat rotates the configuration by -160 degrees, so Miranda's argument of
latitude changes by that much at every repeat and the node window would be missed; the 58 d
repeat shifts it by 4.0 degrees, which keeps the encounter at a fixed point of Miranda's orbit.
That is the first thing a build must check (not computed here: Miranda's nodal regression
under Uranus J2). The other four moons are tilted by 0.2 to 0.7 degrees (z up to 7,100 km at
Oberon); the project's two-moon rows pass the real-ephemeris V4 gate with those tilts.

**Shared-moon V-inf from the existing two-moon rows (my inference, not Liang's text; each row is
one representative of a 30-member family, so indicative only).** A ballistic flyby preserves
|V-inf|, so the two double cyclers joined at the shared hub moon must have similar |V-inf|
there, unless other flybys bridge the gap. Values from `vinf_kms_at_encounters`:

| Moon | rows | values (km/s) | difference |
|---|---|---|---|
| Umbriel | Ariel-Umbriel / Umbriel-Titania | 1.3004 / 1.2296 | 0.07 |
| Oberon | Titania-Oberon / Ariel-Oberon | 1.968 / 1.8286 | 0.14 |
| Ariel | Ariel-Umbriel / Ariel-Titania | 0.979 / 1.2306 | 0.25 |
| Ariel | Ariel-Titania / Ariel-Oberon | 1.2306 / 1.5209 | 0.29 |
| Titania | Ariel-Titania / Titania-Oberon | 1.7186 / 2.1618 | 0.44 |
| Titania | Umbriel-Titania / Ariel-Titania | 1.0058 / 1.7186 | 0.71 |

Ariel-Umbriel-Titania could share Umbriel (0.07) but not Titania; Ariel-Titania-Oberon could
share Oberon (0.14). No catalogue row has Miranda, so Miranda-Ariel-Umbriel has no row-based
check at all; that is a build task.

## 5. Ranking and verdicts

Ranking: by eps within T <= 120 d (Liang's stated condition, with M), W within 120 d as the
secondary figure, then the physical flags of section 4. The top three, Miranda-Ariel-Umbriel,
Miranda-Ariel-Oberon, Miranda-Umbriel-Oberon, are the same in that order on both eps and W, so
they do not depend on the reading of Liang's Eq. 3. All three contain Miranda; that is where the
arithmetic ends and the physics decides.

| Rank | Triple | eps, T <= 120 d (M) | W, T <= 120 d | Flags |
|---|---|---|---|---|
| 1 | Miranda-Ariel-Umbriel | 0.19 deg @ 6.4 d (0.005 d) | 4.3 deg @ 58 d | Miranda passive, i = 4.2 deg; V-inf 1.0 to 2.2 km/s possible |
| 2 | Miranda-Ariel-Oberon | 2.0 @ 83.7 d (0.035 d) | 7.5 @ 81 d | V-inf at least 3.3 km/s at Ariel |
| 3 | Miranda-Umbriel-Oberon | 2.3 @ 113.7 d (0.053 d) | 10.6 @ 54 d | V-inf at least 3.2 km/s at Umbriel |
| 4 | Ariel-Umbriel-Titania | 2.5 @ 102.9 d (0.099 d) | 24.9 @ 96 d (8.3 @ 174 d) | non-Miranda; Umbriel hub V-inf matches |
| 5 | Ariel-Titania-Oberon | 2.5 @ 24.8 d (0.199 d) | 26.4 @ 53 d (9.4 @ 270 d) | non-Miranda; Oberon hub V-inf matches |
| 6 | Miranda-Ariel-Titania | 3.6 @ 35.4 d (0.068 d) | 18.3 @ 35 d | Miranda |
| 7 | Miranda-Titania-Oberon | 4.5 @ 49.0 d (0.327 d) | 26.3 @ 27 d | Miranda |
| 8 | Miranda-Umbriel-Titania | 4.7 @ 23.6 d (0.130 d) | 12.6 @ 8.4 d | Miranda; W repeat under one Titania period |
| 9 | Ariel-Umbriel-Oberon | 5.7 @ 83.7 d (0.196 d) | 21.1 @ 12.7 d | non-Miranda; no matching hub V-inf |
| 10 | Umbriel-Titania-Oberon | 10.0 @ 23.9 d (0.905 d) | 24.5 @ 95 d | non-Miranda; M/T 3.8 percent |

Yardstick: Liang's triple eps 13.6 deg, M 0.7365 d per 49.4 d (0.472 d per 49.6 d with Callisto
as hub), W 12.4 deg at 49.8 d; Io-Europa-Ganymede eps 0 and an inertial shift of 5.2 degrees per
7.05 d (W 3.2 deg at the best repeat). Two-moon tolerance the project already lives with: for the
six rows, 2 x leg time against the synodic period gives integers to within 0.0002 synodic
periods (0.1 degrees) for five rows and 4.991 against 5 (-3.3 degrees, 0.054 d per 29.9 d cycle)
for Umbriel-Oberon. That is by construction (the leg time is free) and the cycle is 1, 3 or 5
synodic periods, so it measures what a free leg time can do, not what a repeat needs; but it
shows the project already accepts about 3 degrees of relative phase per cycle in the ideal
model, with (for the Umbriel-Oberon row) bounded drift of 86,000 to 530,000 km in the real
ephemeris, and windowed duty cycles of 62 to 79 percent across the six rows.

**1. Miranda-Ariel-Umbriel: GO (conditional).** Arithmetic better than Liang's by a factor of 3
in W and 70 in eps (M of 0.005 d against his 0.472 to 0.737 d), though a per-triple 2.3 percent
chance coincidence is about 20 percent for the best of ten, so the GO rests on the margin over
the yardstick, not on significance. Spacecraft orbits satisfying Liang's conditions exist at
V-inf of 1.0 to 2.2 km/s (2/3 P_Umbriel, 21 revolutions: 1.11, 2.03, 1.14 km/s), with two good
bending moons (Ariel, Umbriel) and Miranda used as a passive target body, which Russell and
Strange's ideal model permits (PDF p.4). Conditions: it is the only Uranian triple that clears
the yardstick by a margin that justifies paying the Miranda penalty (no bending, 4.2 degree
tilt); and the first test is the Miranda node geometry, not a search.
**What a build does next:** (a) positive control: `analyze_near_resonance` on Liang's Table 1
(documented to return 7/4 and 0.7365 d) and then on Miranda-Ariel-Umbriel; (b) build the two
hub-sharing doubles with the existing alternating-double-cycler operator, Ariel-Umbriel (an
existing catalogue row, leg 3.216 d, V-inf 0.98/1.30 km/s) and Miranda-Ariel / Miranda-Umbriel
(no rows exist; construct as Russell-Strange free-return doubles with Miranda as the massless
target) in the circular-coplanar model with the spacecraft periods above and T = 58 d (or the
6.44 d basic repeat if the inertial layer is not needed); (c) gate: Miranda encounter distance
from its orbital plane at every repeat, with Miranda's nodal regression included; (d) only then
3D and the windowed V0-V4 gauntlet the two-moon rows passed.
**What would change it to NO-GO:** the Miranda node window cannot be met at every repeat, or the
Miranda doubles do not exist at V-inf of at most 2.2 km/s.

**2. Miranda-Ariel-Oberon: NO-GO.** It passes the yardstick (eps 2.0 degrees, W 7.5 at 81 d) but
it is dominated: the same Miranda burden, a larger and faster orbit (the lowest option found is
half Oberon's period, 12 revolutions, V-inf 1.89, 3.32, 1.42 km/s; the floor for any orbit
spanning Miranda and Oberon is 3.28 km/s at Ariel, where Ariel bends about 1.4 degrees), and a
longer repeat. **What would change it:** a failure of Miranda-Ariel-Umbriel that is specific to
Umbriel (not Miranda), or a Miranda-Oberon double cycler appearing at V-inf below 2.2 km/s.

**3. Miranda-Umbriel-Oberon: NO-GO.** Same reasoning (lowest option 1/2 P_Oberon, 8
revolutions, 1.89, 3.24, 1.42 km/s; floor 3.19 km/s at Umbriel; repeat 113.7 d for the best eps,
53.8 d for the best W). Same change criteria.

**Non-Miranda triples (Ariel-Umbriel-Titania, Ariel-Titania-Oberon; ranks 4 and 5): not excluded
by the arithmetic; decided by one test.** On Liang's stated condition they match
Miranda-Umbriel-Oberon: Ariel-Titania-Oberon has S(T,O) = 24.64 d against 7 S(A,T) = 24.83 d,
leftover 0.199 d (0.80 percent, about half of Liang's 1.49) at a 24.8 d repeat, half of Liang's
50 d, and Ariel-Umbriel-Titania has 16 S(A,U) = 102.92 d against 13 S(U,T) = 102.82 d (0.099 d,
0.1 percent) at 102.9 d, about twice Liang's 50 d half cycle. All their moons are good benders,
both constituent doubles already exist at V4 (hub V-inf agreeing to 0.07 at Umbriel and 0.14 at
Oberon), and the least-energy V-inf are 0.8 to 2.2 km/s. The open points: (a) on the stricter
inertial reading they are 21 to 26 degrees at T <= 120 d, twice Liang's; (b) spacecraft orbits
meeting Liang's conditions at low V-inf exist at their W-optimal repeats (Ariel-Umbriel-Titania,
95.6 d: 1.80, 2.11, 0.92 km/s; Ariel-Titania-Oberon, 52.9 d: 1.28, 1.90, 1.09 km/s) but at the
eps-optimal repeats (102.9 d and 24.8 d) the only orbits found need 3.1 to 3.2 km/s at the middle
moon. **The deciding test:** on the existing Liang reproduction (`search/cge_scaffold.py` members
A/B/C, model in `2026-06-13-liang-abc-reproduction.md`) measure the inertial rotation per
cycle of the spacecraft apse line and the phase drift of the unvisited moon that his cyclers
actually absorb. If the construction tolerates about 25 degrees per half cycle the
non-Miranda triples become the better GO candidates (no Miranda penalty); if it tolerates only
about 12 degrees or less they stay behind Miranda-Ariel-Umbriel. Liang's p.15 results (38 degrees
of drift of the unvisited moon per half cycle) already suggest the former.

## 6. What this study does not establish

- It is arithmetic and reading, not a trajectory. No Lambert arc, no flyby, no ballistic
  closure was computed for any Uranian triple.
- Liang's method has no stated numeric tolerance; the yardstick is his one worked system
  (eps 13.6 deg, W 12.4 deg), taken from his own mean motions and verified against his printed
  S(C,G), S(G,E) and M.
- The 120 d and 730 d caps, the choice of eps with M as the primary figure and W as a secondary
  one, the 3.1 P_out normalisation, the spacecraft-period search and the chance baselines are my
  own choices and constructions; none is in a paper.
- The Keplerian periods are not the observed mean motions; the second set is unverified in this
  session. Matches that rely on the last digits (anything beyond about 120 d) are not secure.
- Eccentricities (0.001 to 0.004) and the 0.2 to 0.7 degree tilts of the four large moons were
  not propagated; Miranda's nodal regression was not computed; the shared-hub V-inf comparison
  uses one representative per pair.
- It does not show a Miranda-Ariel-Umbriel cycler exists; it shows only that its period
  arithmetic is better than Liang's and that nothing in the arithmetic or the V-inf
  estimate forbids a first test.

## Appendix: the arithmetic script

```python
# Standalone arithmetic for #887 stage 1. No repo imports. python3 appendix.py
import itertools, math, random
import numpy as np

# ---- constants (src/cyclerfinder/core/satellites.py: a [km], GM [km^3/s^2], R [km]; PRIMARIES system GM) ----
MU_U, MU_J = 5.7945564e6, 1.26686534e8
A_U = {"Miranda":129846., "Ariel":190929., "Umbriel":265986., "Titania":436298., "Oberon":583511.}
A_J = {"Io":421800., "Europa":671100., "Ganymede":1070400., "Callisto":1882700.}
GM = {"Miranda":4.3,"Ariel":83.5,"Umbriel":85.1,"Titania":226.9,"Oberon":205.3,
      "Io":5959.916,"Europa":3202.739,"Ganymede":9887.834,"Callisto":7179.289}
RAD = {"Miranda":235.8,"Ariel":578.9,"Umbriel":584.7,"Titania":788.9,"Oberon":761.4,
       "Io":1821.49,"Europa":1560.8,"Ganymede":2631.2,"Callisto":2410.3}
SAFE = {"Miranda":100,"Ariel":50,"Umbriel":50,"Titania":50,"Oberon":50}
kep = lambda a, mu: 2*math.pi*math.sqrt(a**3/mu)/86400.0
PU = {k: kep(v, MU_U) for k, v in A_U.items()}          # primary period set (Keplerian from project a)
PJ = {k: kep(v, MU_J) for k, v in A_J.items()}
# sensitivity set B: sidereal periods typed from standard tables, NOT verified in this session
PU_B = {"Miranda":1.413479,"Ariel":2.520379,"Umbriel":4.144176,"Titania":8.705867,"Oberon":13.463234}
# Liang Table 1 mean motions (rad/day), p.4 of the draft PDF
NL = {"Europa":1.7693, "Ganymede":0.8782, "Callisto":0.3765}
PL = {k: 2*math.pi/v for k, v in NL.items()}
PJ_B = {"Io":1.769138,"Europa":3.551181,"Ganymede":7.154553}   # typed from standard tables, NOT verified here

syn = lambda p, q: 1.0/abs(1.0/p - 1.0/q)
wrap = lambda x: (x + 0.5) % 1.0 - 0.5

def Wbest(P, tri, Tc, step=0.004):
    """Return best (T, W deg, shifts) with T<=Tc; W = max over the 3 moons of |inertial shift| after T."""
    Ps = [P[n] for n in tri]
    T = np.arange(0.9*max(Ps), Tc, step)
    sh = np.array([wrap(T/p)*360 for p in Ps]); W = np.max(abs(sh), axis=0)
    idx = [i for i in range(1, len(T)-1) if W[i] <= W[i-1] and W[i] <= W[i+1]]
    i = min(idx, key=lambda i: W[i]) if idx else int(np.argmin(W)); return float(T[i]), float(W[i]), [round(float(sh[k][i]),1) for k in range(3)]

def eps_best(P, tri, Tc):
    """Relative-phase layer: integers (n_ab,n_bc), balanced repeat T=(n_ab+n_bc)*S_ac; eps = phase error of pair ab (= minus that of bc)."""
    a, b, c = tri; Sab, Sbc, Sac = syn(P[a],P[b]), syn(P[b],P[c]), syn(P[a],P[c]); best = None
    for nab in range(1, 400):
        for nbc in range(1, 400):
            T = (nab+nbc)*Sac
            if T > Tc: break
            e = abs(360*(T/Sab - nab)); M = nab*Sab - nbc*Sbc
            if best is None or e < best[0]: best = (e, nab, nbc, T, M)
    return best

def vinf_on_orbit(rp, ra, r, mu):
    a = (rp+ra)/2; v2 = mu*(2/r - 1/a); h = math.sqrt(mu*a*(1-((ra-rp)/(ra+rp))**2))
    vt = h/r; vr = math.sqrt(max(v2-vt*vt, 0)); vc = math.sqrt(mu/r); return math.hypot(vr, vt-vc)

def bend(m, vinf, h):
    s = GM[m]/((RAD[m]+h)*vinf**2 + GM[m]); return 2*math.degrees(math.asin(s)), 2*vinf*s

if __name__ == "__main__":
    print("== 1. periods (Keplerian from project a and system GM) [d]")
    for k, v in PU.items(): print(f"  {k:9s} {v:9.5f}   set B {PU_B[k]:9.6f}   n={360/v:9.4f} deg/d")
    for k, v in PJ.items(): print(f"  {k:9s} {v:9.5f}")
    print("== 2. control: Liang Table 1 (CGE)")
    SCG, SGE, SCE = syn(PL["Callisto"],PL["Ganymede"]), syn(PL["Ganymede"],PL["Europa"]), syn(PL["Callisto"],PL["Europa"])
    print(f"  S_CG={SCG:.4f} S_GE={SGE:.4f} S_CE={SCE:.4f}  ratio={SCG/SGE:.4f}  M=4S_CG-7S_GE={4*SCG-7*SGE:.4f} d (paper: 12.5232, 7.0509, 1.7761, 0.7365)")
    print(f"  M/T = {(4*SCG-7*SGE)/(7*SGE)*100:.2f} %; 3 repeats = {3*(4*SCG-7*SGE):.3f} d = {360*3*(4*SCG-7*SGE)/SGE:.0f} deg of G-E phase (paper p.9: 'about 100 degrees' in about 150 d)")
    print(f"  4*S_CG={4*SCG:.3f}  11*S_CE={11*SCE:.3f}  7*S_GE={7*SGE:.3f}  3P_C={3*PL['Callisto']:.3f} 7P_G={7*PL['Ganymede']:.3f} 14P_E={14*PL['Europa']:.3f}")
    for lab, T_ in (("4 S_CG", 4*SCG), ("11 S_CE", 11*SCE), ("7 S_GE", 7*SGE)):
        print(f"  Liang half-cycle T={lab}={T_:.3f} d: inertial shifts [deg] " + ", ".join(f"{m} {360*wrap(T_/PL[m]):+.1f}" for m in ("Callisto","Ganymede","Europa")))
    print(f"  hub mismatches (Liang hub = Callisto): |4 S_CG - 11 S_CE| = {abs(4*SCG-11*SCE):.3f} d ({abs(4*SCG-11*SCE)/(11*SCE)*100:.2f} %)")
    print("== 3. yardsticks and ten Uranian triples (primary period set)")
    print("  name | P ratios | synodics S_ab S_bc S_ac | eps(T<=60) | eps(T<=120) | W(T<=60) | W(T<=120) | W(T<=365) | W(T<=730) | Wset B (<=120, <=730)")
    rows = [("Liang CGE (Table 1)", PL, ("Callisto","Ganymede","Europa")), ("Io-Europa-Ganymede", PJ, ("Io","Europa","Ganymede"))]
    rows += [("-".join(t), PU, t) for t in itertools.combinations(sorted(PU, key=PU.get), 3)]
    for nm, P, tri in rows:
        a, b, c = tri; S = (syn(P[a],P[b]), syn(P[b],P[c]), syn(P[a],P[c]))
        e60, e120 = eps_best(P, tri, 60), eps_best(P, tri, 120)
        w = [Wbest(P, tri, Tc) for Tc in (60,120,365,730)]
        wb = ""
        if P is PU:
            wb = " | B: " + " ".join(f"{Wbest(PU_B, tri, Tc)[1]:.1f}@{Wbest(PU_B, tri, Tc)[0]:.0f}" for Tc in (120,730))
        print(f"  {nm} | {P[b]/P[a]:.4f} {P[c]/P[b]:.4f} | {S[0]:.4f} {S[1]:.4f} {S[2]:.4f} | e {e60[0]:.2f}deg n=({e60[1]},{e60[2]}) T={e60[3]:.1f} M={e60[4]:.3f}d M/T={abs(e60[4])/e60[3]*100:.2f}% | e {e120[0]:.2f} n=({e120[1]},{e120[2]}) T={e120[3]:.1f} M={e120[4]:.3f}"
              f" | " + " | ".join(f"W {x[1]:.1f}@{x[0]:.1f}d T/Pout={x[0]/max(P[n] for n in tri):.2f}" for x in w) + wb)
    T0 = syn(PJ["Europa"], PJ["Ganymede"])
    print(f"== 3a. control (Hernandez p.2: 7.05 d, 5.2 deg): at T=S_Europa-Ganymede={T0:.4f} d inertial shifts [deg]: " + ", ".join(f"{m} {360*wrap(T0/PJ[m]):.1f}" for m in ("Io","Europa","Ganymede")))
    T1 = 2*syn(PU["Miranda"], PU["Ariel"])
    print(f"== 3a'. Miranda-Ariel-Umbriel basic (2,1,3) repeat T={T1:.4f} d: inertial shifts [deg] " + ", ".join(f"{m} {360*wrap(T1/PU[m]):.1f}" for m in ("Miranda","Ariel","Umbriel")))
    print("== 3b. Liang-style normalisation: best W for T <= 3.1 x P_outer (Liang's repeat is 2.98 Callisto periods)")
    for nm, P, tri in rows:
        po = max(P[n] for n in tri); t = Wbest(P, tri, 3.1*po); print(f"  {nm}: W {t[1]:.1f} deg at T={t[0]:.1f} d (T/P_out={t[0]/po:.2f})")
    print("== 3c. shifts [deg] of each moon (listed in the order of the triple) at the W-optimal T, T<=120 and T<=730")
    for nm, P, tri in rows:
        a1, a2 = Wbest(P, tri, 120), Wbest(P, tri, 730)
        print(f"  {nm}: <=120: T={a1[0]:.1f} shifts {a1[2]} (max pair diff {max(a1[2])-min(a1[2]):.1f}); <=730: T={a2[0]:.1f} shifts {a2[2]}")
    print("== 4. chance baseline (own construction): random sorted log-uniform period triples, best W for T<=60 d")
    random.seed(7)
    for lab, lo, hi, ref in (("Uranian range 1.41-13.46 d", 1.41, 13.46, (4.4, 10.5, 12.3, 18.3, 24.9)), ("Jovian range 3.55-16.69 d", 3.55, 16.69, (12.3,))):
        v = []
        for _ in range(1200):
            ps = sorted(math.exp(random.uniform(math.log(lo), math.log(hi))) for _ in range(3))
            v.append(Wbest({"a":ps[0],"b":ps[1],"c":ps[2]}, ("a","b","c"), 60, step=0.01)[1])
        v3 = []
        for _ in range(1200):
            ps = sorted(math.exp(random.uniform(math.log(lo), math.log(hi))) for _ in range(3))
            v3.append(Wbest({"a":ps[0],"b":ps[1],"c":ps[2]}, ("a","b","c"), 3.1*ps[2], step=0.002*ps[2])[1])
        v3.sort(); print(f"  {lab}, T<=3.1 P_out: median {v3[600]:.1f}, 10th {v3[120]:.1f}; P(W<=12.4)={sum(1 for y in v3 if y<=12.4)/len(v3):.3f}; P(W<=4.3)={sum(1 for y in v3 if y<=4.3)/len(v3):.3f}")
        v.sort(); print(f"  {lab}: median {v[600]:.1f}, 10th pct {v[120]:.1f}; " + "; ".join(f"P(W<={x})={sum(1 for y in v if y<=x)/len(v):.3f}" for x in ref))
    def best_eps_rho(rho, N):
        return min(abs(360*(nbc - nab*rho)/(1+rho)) for nab in range(1, N) for nbc in range(1, N-nab+1))
    ve = sorted(best_eps_rho(math.exp(random.uniform(math.log(0.1), math.log(10))), 11) for _ in range(5000))
    print(f"  relative-phase eps, order n_ab+n_bc <= 11, random synodic ratio: median {ve[2500]:.1f} deg; P(eps<=13.6)={sum(1 for y in ve if y<=13.6)/len(ve):.2f}; P(eps<=0.19)={sum(1 for y in ve if y<=0.19)/len(ve):.3f}")
    print("== 5. Laplace-type angle n_1-3n_2+2n_3 [deg/d]")
    for lab, J, U in (("Keplerian (project a)", PJ, PU), ("set B (unverified tables)", PJ_B, PU_B)):
        nI = {k: 360/v for k, v in J.items()}; nU = {k: 360/v for k, v in U.items()}
        x = nU['Miranda']-3*nU['Ariel']+2*nU['Umbriel']
        print(f"  {lab}: Io,Europa,Ganymede {nI['Io']-3*nI['Europa']+2*nI['Ganymede']:.5f} ; Miranda,Ariel,Umbriel {x:.4f}  circulation period {360/abs(x)/365.25:.1f} yr")
    print("== 6. spacecraft orbit with T_sc = outer-moon period (Liang condition 1-2 style), apojove side fixed by period; rp=a_inner and rp=0.5 a_inner")
    def orb(P, A, mu, tri, T):
        k = round(T/P[tri[2]]); Tsc = T/k; a = (Tsc*86400*math.sqrt(mu)/(2*math.pi))**(2/3.0)
        out = []
        for f in (1.0, 0.5):
            rp = f*A[tri[0]]; ra = 2*a-rp
            out.append((f, rp, ra, [round(vinf_on_orbit(rp, ra, A[m], mu),2) for m in tri]))
        return k, Tsc, a, out
    print("  control (Liang Table 3 first row: Callisto 5.6730, Ganymede 6.9919 km/s; perijove r_Eu-10000 = 660988 km, T_sc=T_Callisto):")
    PLk = {"Europa":PJ["Europa"],"Ganymede":PJ["Ganymede"],"Callisto":PJ["Callisto"]}
    a = A_J["Callisto"]; rp = 660988.0; ra = 2*a-rp
    print("   ", {m: round(vinf_on_orbit(rp, ra, A_J[m], MU_J),3) for m in ("Callisto","Ganymede","Europa")})
    for nm, tri, T in (("M-A-U",("Miranda","Ariel","Umbriel"),57.97),("M-A-O",("Miranda","Ariel","Oberon"),80.6),("M-U-O",("Miranda","Umbriel","Oberon"),53.75),("A-U-T",("Ariel","Umbriel","Titania"),95.6),("A-T-O",("Ariel","Titania","Oberon"),52.9)):
        k, Tsc, a, out = orb(PU, A_U, MU_U, tri, T)
        print(f"  {nm}: T={T} d, k={k} (T/P_out={T/PU[tri[2]]:.2f}), T_sc={Tsc:.3f} d, a={a:.0f} km")
        for f, rp, ra, v in out: print(f"      rp={rp:.0f} ra={ra:.0f}: Vinf {dict(zip(tri, v))}")
    print("  admissible T_sc = (q/p) P_moon, p,q<=4 (Liang condition 2: integer ratio), T = k T_sc within 0.08 of an integer k; apses just reach the outer moon; ranked by largest V-inf:")
    cases = (("M-A-U",("Miranda","Ariel","Umbriel"),(58.0,)),("M-A-O",("Miranda","Ariel","Oberon"),(80.6,83.7)),("M-U-O",("Miranda","Umbriel","Oberon"),(53.8,113.7)),("A-U-T",("Ariel","Umbriel","Titania"),(95.6,102.9)),("A-T-O",("Ariel","Titania","Oberon"),(24.8,52.9)))
    from math import gcd
    for nm, tri, Ts in cases:
        for T in Ts:
            opts = []
            for m in tri:
                for p in range(1, 5):
                    for q in range(1, 5):
                        if gcd(p, q) > 1: continue
                        Tsc = PU[m]*q/p; k = round(T/Tsc)
                        if k < 1 or abs(T/Tsc-k) > 0.08: continue
                        a = (Tsc*86400*math.sqrt(MU_U)/(2*math.pi))**(2/3.0); rp = min(A_U[tri[0]], 2*a-A_U[tri[2]]); ra = 2*a-rp
                        if rp < 25559+500 or ra < A_U[tri[2]]-1: continue
                        v = [vinf_on_orbit(rp, ra, A_U[x], MU_U) for x in tri]
                        opts.append((max(v), f"T_sc={q}/{p} P_{m}={Tsc:.3f} d, {k} revs, rp={rp:.0f} ra={ra:.0f}: Vinf {[round(x,2) for x in v]}"))
            opts.sort()
            print(f"    {nm} T={T}: " + (f"{len(opts)} options; lowest max-Vinf:" if opts else "none"))
            for o in opts[:3]: print("        " + o[1])
    print("  least-energy orbit (rp=a_in, ra=a_out), all ten triples: lower bound on Vinf")
    for tri in itertools.combinations(sorted(A_U, key=A_U.get), 3):
        rp, ra = A_U[tri[0]], A_U[tri[2]]
        print("   ", "-".join(tri), {m: round(vinf_on_orbit(rp, ra, A_U[m], MU_U),2) for m in tri}, f"T_orb={kep((rp+ra)/2, MU_U):.2f} d")
    print("== 7. maximum bend [deg] / max dv [m/s] at project safe altitude")
    for m in ("Miranda","Ariel","Umbriel","Titania","Oberon"):
        vc = math.sqrt(MU_U/A_U[m]); print(f"  {m}: v_circ {vc:.2f} km/s; " + "; ".join(f"v_inf {v}: {bend(m,v,SAFE[m])[0]:.1f} deg / {1000*bend(m,v,SAFE[m])[1]:.0f} m/s ({bend(m,v,SAFE[m])[1]/vc*100:.1f}% of v_circ)" for v in (1,2,3,4)))
    print("  Liang Tables 3/5/7 flybys re-evaluated with his Eq. 7 (his altitudes ignore Jupiter gravity, p.18):")
    for tab, rows_ in (("T3",[("Callisto",5.6698,33839),("Europa",4.6685,6241),("Callisto",5.8721,19765)]),("T5",[("Callisto",5.6698,60877),("Europa",4.4853,10825),("Callisto",5.7914,36258)]),("T7",[("Callisto",7.6409,27438),("Europa",12.0213,3636),("Callisto",7.7838,12516)])):
        for m, v, h in rows_:
            d, dv = bend(m, v, h); vc = math.sqrt(MU_J/A_J[m]); print(f"    {tab} {m} v_inf {v} h {h} km: bend {d:.2f} deg, dv {1000*dv:.0f} m/s = {dv/vc*100:.2f}% of v_circ")
    print("== 8. Miranda geometry: i=4.232 (Pergola Table 1) ; z_max = a sin i; out-of-plane speed at node = v_circ sin i")
    i = math.radians(4.232); vc = math.sqrt(MU_U/A_U['Miranda']); print(f"  z_max {A_U['Miranda']*math.sin(i):.0f} km ({A_U['Miranda']*math.sin(i)/RAD['Miranda']:.0f} Miranda radii); v_circ sin i = {vc*math.sin(i):.3f} km/s; window |z|<2000 km: |u|<{math.degrees(math.asin(2000/(A_U['Miranda']*math.sin(i)))):.1f} deg of a node")
```
