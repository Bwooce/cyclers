# Digest: Campagnola, Buffington, Lam, Petropoulos & Pellegrini 2019, "Tour Design Techniques for the Europa Clipper Mission" (#960)

S. Campagnola, B. B. Buffington, T. Lam, A. E. Petropoulos and E. Pellegrini (JPL), "Tour Design Techniques
for the Europa Clipper Mission", Journal of Guidance, Control, and Dynamics 42(12):2615-2626 (2019), DOI
10.2514/1.G004309. Crossref-confirmed 2026-10-05: title, five authors, 42(12), pp. 2615-2626.

Filed as `cyclers_pdf/papers/campagnola-buffington-lam-petropoulos-pellegrini-2019-tour-design-techniques-europa-clipper-JGCD-42-12-2615-doi-10.2514-1.G004309.pdf`.
- The PDF is the 12-page "Article in Advance" version: text layer, pages numbered 1-12, md5
  da2d99abfa33f9707f6e1f765f5beb75. Page numbers below are Article-in-Advance pages.
- I read all pages from the text layer. Pages 9 and 10 (Figs. 11-13, Table 3) were also read as page
  images.
- The large tour tables (Table 1, pp.3-6) were read only for the structure noted below. Their text layer is
  scrambled into columns.

Evidence tags: READ = printed (page cited). INFERRED = my reading. Figure-derived values are flagged.

Why it was acquired: `#938` sec. 5 item 4 is the literal-collision check for:
- `#943` X1 (two-working-body Ganymede-Callisto).
- `#945` R2 (Io-containing Jovian triples).
- `#953` R7 (Laplace-locked triple cycler).
It is also the source of the "Callisto pi-transfer".

## 0. Verdict

This is a mission-design paper about the Europa Clipper candidate tour 18F17. It does contain one item that
bears directly on `#943`:
- **A published two-working-body Ganymede-Callisto cycler, "GCGC" (sec. III.C, Fig. 11, p.9).** It is
  given as an original contribution, but with no numbers beyond Table 3's summary row.
- So the both-moons-bending Ganymede-Callisto class is no longer unpublished. **X1 Ganymede-Callisto moves
  from OPEN to PARTIAL.** What is left: any numeric or ideal-model periodic member, and its states.
- The paper has no Ganymede-Europa cycler, no Io and no repeating triple. So X1 Ganymede-Europa, R2 and R7
  stay as they were.

## 1. Model and vocabulary (READ pp.2-3)

- Design model: "patched conics ... patched by instantaneous Delta-Vs, to simulate each satellite flyby", with
  moons on circular coplanar orbits and direct orbits around Jupiter only (p.2).
- Final tours are converged in high fidelity: jTOP, then JPL's COSMIC (p.3).
- Flyby variables: v_inf, pump angle alpha and crank angle kappa on the v_inf sphere (Fig. 2).
- Transfer types (p.2):
  - resonant n:m. "n is the number of moon revolutions and m is the number of spacecraft revolutions"; in
    some literature the ratio is m:n.
  - nonresonant n:m+ (long) or n:m- (short), which are in-plane.
  - n-odd-pi transfers, typically inclined.
- Sequences:
  - resonant hopping.
  - petal rotation: alternating pump-up and pump-down nonresonant pairs n1:m1(+-)/n2:m2(+-) (ref. [16],
    Anderson et al. 2018).
  - crank-over-the-top (COT) sequences.
  - pi-transfer sequences: half-COT, a pi transfer, then half-COT.
  - switch-flip: a pi transfer at a different moon (p.3).

## 2. Tour 18F17 (READ pp.3-8)

- 45 Europa, 6 Ganymede and 9 Callisto flybys over about 4 years.
- Flyby naming is xxByy: xx is the spacecraft revolution count, B the body and yy the body count (p.3).
- Phases (p.3):
  - Ganymede pump-down.
  - COT-1 to COT-3: 4:1 resonances at Europa.
  - nonresonant 7:2 transfers to avoid solar conjunction.
  - 4:1/5:1- petal rotation.
  - **Transition to Europa Campaign 2 by a low-dose Callisto petal rotation with Ganymede leveraging**: 06G05
    leverages down v_inf at Callisto, then the Callisto petals, then the reverse.
  - a leading-edge flyby E29: outgoing 5:1/4:1-, 25 km altitude.
  - COT-4 and COT-5 at 5:1, then COT-6 at 4:1.
  - Ganymede-hopping disposal.
- Table 2 compares 18F17 with the PDR baseline 17F12 on Delta-V, dose, TOF and eclipses.
- Callisto-Ganymede-Callisto transfers are used to reach a 1:1 Callisto orbit from low Ganymede v_inf, "both
  at the beginning (34C01-35G05-37C02) and at the end (46C07-48G06-49C08)" (p.8).
- The (rho, v_inf) graphs are Figs. 7, 8 and 10, at Ganymede, Europa and Callisto. Example: Ganymede 2:1 to
  4:1 hopping is possible only for 3.5 < v_inf < 4.5 km/s, and 4:1 to 10:1 has v_inf,max about 6.5 km/s
  (p.8).

## 3. Sec. III.C: 180-degree rotations by cycler, pi-transfer or petals (READ pp.8-10)

**(a) Ganymede-Callisto cyclers (p.9, Fig. 11).** Verbatim:
- "The first technique is Ganymede-Callisto cyclers. Cyclers are trajectories repeatedly transferring a
  spacecraft between two bodies [24]. Their application to fast line-of-apses rotation is an original
  contribution of this paper, and is inspired by the Clipper tour 15F09."
- "These types of cyclers can be modeled with n parameters, which are found solving a system of n nonlinear
  algebraic equations. Cyclers do not exist in families but as individual solutions, which can be found
  numerically."
- "The details of the cycler technique for fast line-of apses rotations will be discussed in future papers;
  here we only present one example used for the Europa Clipper application."
- "The example is shown in Fig. 11, optimized in a high-fidelity model (jTOP), with two cycles each providing
  ~90 deg rotation. ... First guess solutions are designed with the help of Star [25]."
- Ref. [24] is Russell & Strange 2009. Ref. [25] is Landau, AAS 15-585 (STAR).
- Fig. 11 caption: "The 180 deg rotation of the spacecraft orbit and of the flyby location, from G1 to G9,
  repeating two GCGC cycles (G1-G5 and G5-G9)."
- Fig. 11 plots, in km, labelled flybys G1, C2, G3, C4, G5 (C6, G7, ...) and G9. Both moons host flybys, and
  each cycle is G-C-G-C-G.

Table 3 (p.10), "rotation only" row for the cycler: TID 300 krad, TOF 5.5 months, v_inf 3.5/4.5 km/s. Its
qualitative rows:
- Phasing: "Fixed v_inf and phasing".
- Untargeted flybys: "Can disrupt the phasing invalidating the trajectory".
- Solar conjunction: same.
- Eclipses: "Few chances of mid eclipses".
- Footnote: "The cycler approach was not implemented into an end-to-end tour, and so the full cost of the
  transition to EC2 is not available."

**(b) Callisto pi-transfer sequence with 1:1 resonant transfers (p.9).** "The number of flybys in the
sequence ... is function of the v_inf". Table 3: TID 25 krad, TOF 4 months, v_inf 2.7 km/s, "Mainly
out-of-plane orbits with no untargeted flybys". In the 17F13 tour it is used as a switch-flip, and the dose
is accrued mostly in Europa flybys E26-E33 (pp.9-10).

**(c) Callisto petal rotation (p.9-10, Figs. 12-13).**
- Fig. 12 gives TOF and TID for 180 deg as functions of Callisto v_inf (2.5-4.5 km/s) for these families:
  1:1+/2:1-, 1:1+/3:3-, 4:4-/3:2+, 3:3-/3:2+, 2:2-/3:2+ and 1:1+/2:2-.
- The 18F17 design point (a) and the 17F12 point (b) are marked.
- "The period for one n1:m1/n2:m2 petal is about n1 + n2 Callisto revolutions" (p.9).
- Table 3: TID 60 krad, TOF 5 months, v_inf 3.6 km/s.
- 18F17 uses 1:1+/3:3-/1:1+/2:2-/2:2+ (C02-C07 in the text, C02-C08 in the Fig. 13 caption).
- "The trajectory shadows three periodic orbits and their heteroclinic connections in the Callisto-Jupiter
  circular, restricted, three-body problems (CR3BP). Each periodic orbit is associated to single
  n1:n1/n2:n2 pair" (p.10). Fig. 13 shows the three as shapes only, with no initial conditions.

**(d) Leading-point flybys (pp.10-11).** These give the altitude over the leading or trailing point as a
function of the mean pump angle, Eqs. (1)-(2), with families in Fig. 16. They are not related to cyclers.

## 4. Gate answers

- **`#943` X1 Ganymede-Callisto: PARTIAL.**
  - The two-working-body G-C-G-C cycler class is published. It is described as solving n nonlinear equations
    for isolated solutions, and one example was converged in a high-fidelity model.
  - The published numbers are only Table 3's v_inf 3.5/4.5 km/s and the 180 deg totals.
  - Not published: states, per-flyby v_inf or altitude, period, cost, or any ideal-model periodic member. The
    "future papers" have not been located (follow-up below).
  - Policy consequence: a `#943` Ganymede-Callisto closure is "first numeric / first ideal-model periodic
    member of a published class". Attribute Campagnola et al. 2019. It is not "first of class". It must be
    compared against Fig. 11, which can only be digitised, so any comparison is figure-derived.
  - INFERRED: the "v_inf 3.5/4.5" pair is most likely Ganymede/Callisto in that order. The table does not
    say which is which.
- **`#943` X1 Ganymede-Europa: still OPEN.** The paper has no Ganymede-Europa cycler.
- **`#945` R2 (Io triples) and `#953` R7 (Laplace triple cycler): no literal collision.** No Io flyby, no
  Io-Europa-Ganymede repetition and no triple cycler appear.
- **The "Callisto pi-transfer"** is a tour-segment technique: a COT plus a pi transfer at Callisto with 1:1
  resonances, not a cycler. It is not prior art for a cycler row.
- **Callisto petals** are near-periodic "shadowing" of three Callisto-Jupiter CR3BP periodic orbits, with no
  numbers. That is prior art only for a symmetric n:n resonant family claim at Callisto.

## 5. Positive controls

There are no sourced numeric states. Usable sourced numbers are limited to:
- Table 3: cycler v_inf 3.5/4.5 km/s; 5.5 months and 300 krad for 180 deg in two GCGC cycles.
- Fig. 11 flyby order G1, C2, G3, C4, G5, with about 90 deg apse rotation per cycle (READ, caption and text).

They are useful as shape and energy targets for a `#943` Ganymede-Callisto member. They are not a
reproduction control.

## 6. Citation mining (policy step 4)

The references [1]-[25] were checked against `CORPUS_INDEX.md`, the digests and the filenames on
2026-10-05.

Held:
- [2] Wolf & Smith 1995.
- [8] Campagnola et al. 2015 Neptune-Triton.
- [13] Campagnola & Russell 2010 Endgame part 2 (held as AAS 09-227).
- [24] Russell & Strange 2009.

Not held, in priority order. The DOI was checked through the Crossref API on title, authors, volume and
pages.

1. Landau, D. F., "Efficient Maneuver Placement for Automated Trajectory Design", AAS 15-585 (2015). The
   journal version is JGCD 41:1531-1541 (2018), doi 10.2514/1.G003172 (CONFIRMED). This is STAR, which
   made the GCGC first guesses. It is the most direct follow-up for the `#943` collision check.
2. Anderson, R. L., Campagnola, S. & Buffington, B. B., "Analysis of Petal Rotation Trajectory
   Characteristics", JGCD 41(4):827-840 (2018), doi 10.2514/1.G002571 (CONFIRMED). Petal families as CR3BP
   periodic orbits.
3. Lam, T., Buffington, B. B. & Campagnola, S., "A Robust Mission Tour for NASA's Planned Europa Clipper
   Mission", AIAA 2018-0202, doi 10.2514/6.2018-0202 (CONFIRMED). Tours 17F12 and 17F13.
4. Buffington, B. B., Campagnola, S. & Petropoulos, A. E., "Europa Multiple-Flyby Trajectory Design", AIAA
   2012-5069, doi 10.2514/6.2012-5069 (CONFIRMED). COT definitions.
5. Campagnola, S. & Kawakatsu, Y., "Three-Dimensional Resonant Hopping Strategies and the Jupiter
   Magnetospheric Orbiter", JGCD 35(1):340-344 (2012), doi 10.2514/1.53334 (CONFIRMED).
6. Uphoff, C., Roberts, P. H. & Friedman, L. D., "Orbit Design Concepts for Jupiter Orbiter Missions", JSR
   13(6):348-355 (1976), doi 10.2514/3.57096 (CONFIRMED). The origin of the n-odd-pi transfer.
7. Buffington, B. B., "Trajectory Design for the Europa Clipper Mission Concept", AIAA 2014-4105, doi
   10.2514/6.2014-4105 (CONFIRMED; the Crossref title reads "Trajectory Design Concept for the Proposed Europa
   Clipper Mission").
8. Strange, N. J., Campagnola, S. & Russell, R. P., "Leveraging Flybys of Low Mass Moons to Enable an
   Enceladus Orbiter", Adv. Astronaut. Sci. 135:2207-2225 (2010; AAS 09-208). No DOI (AAS paper).
9. Tardivel, S., Klesh, A. & Campagnola, S., "Average Daily Radiation Dose ...", JSR 54(6):1367-1375 (2017),
   doi 10.2514/1.A33790 (CONFIRMED). Radiation, low priority.
10. Lower priority, no DOI (AAS, IAC, ISSFD papers and JPL reports):
    - Buffington, Strange & Campagnola, ISSFD 2012.
    - Boutonnet & Schoenmaekers, AAS 12-207.
    - Lam, Arrieta-Camacho & Buffington, AAS 15-657.
    - Scott et al., AAS 17-437.
    - Buffington et al., IAC 2017.
    - Smith & Buffington, Adv. Astronaut. Sci. 135.
    - Campagnola et al., ISSFD 2015.
    - Evans, JPL IOM 2018.
    - Garrett et al., JPL Publication.
    - Jordan, NOVICE 1982.
    - Campagnola, Takashima & Kawakatsu, JAXA.
    - Hand 2009, Nature (not checked).
11. Not located: the "future papers" on the cycler technique for fast apse rotation (sec. III.C). A forward
    citation search on this paper is the follow-up.

## 7. Novelty-gate anchor

A `KNOWN_CORPUS` anchor has been added in `src/cyclerfinder/search/literature_check.py`:
- key: `campagnola-2019-europa-clipper-gcgc`
- primary: Jupiter
- body_set: {Ganymede, Callisto}
- topology: repeated-moon

A Ganymede-Callisto cycler candidate now meets this paper as well as the Russell & Strange 2009 GanCal
anchor.

## 8. Addendum: the sourced Ganymede-Callisto sequence numbers (READ, Table 1 continued, p.6, page image at 220 dpi)

Transition to Europa Campaign 2 in tour 18F17. This is a high-fidelity tour, so these are flown-design
values, not ideal-model invariants. n and m are the numerically measured revolution counts of the leg after
the flyby (p.3).

| Flyby | Date (ET) | Alt, km | v_inf, km/s | n | m |
|---|---|---|---|---|---|
| 33E26 | 21-May-2027 | 1932 | 3.77 | 4.53 | 1.35 |
| 34C01 | 06-Jun-2027 | 135 | 4.84 | 0.92 | 0.52 |
| 35G05 | 22-Jun-2027 | 1302 | 4.71 | 3.30 | 1.82 |
| 37C02 | 15-Jul-2027 | 25 | 3.62 | 1.56 | 1.56 |
| 38C03 | 10-Aug-2027 | 1614 | 3.61 | 2.30 | 2.30 |
| 41C04 | 18-Sep-2027 | 1633 | 3.62 | 1.56 | 1.56 |
| 42C05 | 14-Oct-2027 | 543 | 3.63 | 1.26 | 1.26 |
| 44C06 | 04-Nov-2027 | 1292 | 3.63 | 2.59 | 2.59 |
| 46C07 | 17-Dec-2027 | 2772 | 3.63 | 1.59 | 1.54 |
| 48G06 | 12-Jan-2028 | 4657 | 3.20 | 3.13 | 1.39 |
| 49C08 | 04-Feb-2028 | 5644 | 4.37 | 1.00 | 1.00 |
| 50C09 | 21-Feb-2028 | 176 | 4.38 | 1.25 | 1.54 |

- The two Callisto-Ganymede-Callisto transfers are 34C01-35G05-37C02 and 46C07-48G06-49C08 (p.8).
- In these, Ganymede v_inf is 4.71 and 3.20 km/s, and Callisto v_inf is 4.84 / 3.62 and 3.63 / 4.37
  km/s.
- These are tour legs, not a repeating cycler. They are the only sourced per-flyby numbers in the paper
  for a G-C sequence.
- Caution: their v_inf range (3.2-4.8 km/s) overlaps the GCGC cycler's 3.5/4.5 km/s.

Follow-up for 15F09, which "inspired" the GCGC cycler: the 2015 Clipper trade-study paper Lam,
Arrieta-Camacho & Buffington, AAS 15-657, is a likely published description of 15F09. This is unverified.
It moves up the missing-papers list, to just after Landau's STAR paper.
- **Resolved 2026-10-05 (batch 9):** AAS 15-657 is filed and digested
  (`2026-10-05-digest-lam-arrieta-camacho-buffington-2015-europa-mission-flyby-trades.md`). Its tour
  15F9-A22 (INFERRED to be 15F09) rotates 180 deg with a Callisto-only switch flip (E23, C1-C7 with a
  C4-C5 pi-transfer, E24). It has no G-C chain. The only G-C-G-C chain in that paper is the non-resonant
  13F7-A21 pump-down. X1 Ganymede-Callisto stays PARTIAL; there is no literal collision.
