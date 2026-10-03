# Digest — Pergola, Geurts, Casaregola, Andrenucci (2007), "Three Body Invariant Manifold Transition with Electric Propulsion" (IEPC-2007-305)

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, 13 pages; Tables 1-3 and Figures 5-10
checked against the rendered pages as well as the text layer). **Purpose:** the project must know exactly
what in-system three-body trajectory work at Uranus is published. This is the only paper in the corpus that
builds planar circular restricted three-body (PCR3BP) invariant-manifold trajectories for Uranus-moon pairs.
This note is facts only; no novelty verdict is written here. Where a number is not printed it is marked
"not stated"; arithmetic of mine is labelled as such.

## Citation
P. Pergola, K. Geurts, C. Casaregola, M. Andrenucci, "Three Body Invariant Manifold Transition with Electric
Propulsion", IEPC-2007-305, 30th International Electric Propulsion Conference, Florence, Italy, 17-20
September 2007 (header as printed: "Presented at the 30th International Electric Propulsion Conference,
Florence, Italy, September 17-20, 2007").

- Affiliations as printed: Pergola, Geurts, Casaregola (Ph.D. students, Dep. of Aerospace Engineering,
  University of Pisa; research engineers, Alta S.p.A., Pisa); Andrenucci (CEO Alta S.p.A.; Professor,
  University of Pisa). Open access from the conference archive.
- Filed in the private paper corpus as
  `pergola-geurts-casaregola-andrenucci-2007-three-body-invariant-manifold-transition-electric-propulsion-uranus-moons-IEPC-2007-305.pdf`,
  md5 `a9ad6d010cafe869a9a8b00b2f71f4f6`.

## What the paper is
Abstract (quoted): "The work presented in this paper was aimed at analyzing the possibility to successfully
exploit the characteristics of three-body systems using different manifolds, combined with Electric Propulsion
within the Uranus planetary system. Electric Propulsion is used to link the low energy paths between stable
and unstable invariant manifolds computed for the different Uranus - moon systems. A tour of the planetary
system involving orbiting phases around the five main moons was designed. This mission requires only a
minimum of onboard resources and can be carried out with an acceptable transfer duration."

Sections as printed: I Introduction (A Planetary System); II Mathematical Models (A System Dynamics; B Design
Approach; 1 Non Linear Programming); III Results; IV Conclusion; References (12). Three tables, ten figures.
One tour is designed and one set of totals is printed; the paper is a single-case study, not a survey.

Stated aim (quoted): "The main purpose of the proposed study is to develop a tour within the Uranus system
that visits each one of the selected moons, including a temporary capture obtained by a ballistic arrival and
departure arc." The tour (quoted): "The trajectory starts with a ballistic approach arc towards the outer
moon, Oberon, and ends with a stable orbit around Uranus with a radius smaller than the inner moon Miranda's
distance. During the complete transfer several closed orbits around Oberon, Titania, Umbriel, Ariel and
Miranda are executed. This would enable scientific studies of the main moons for a considerably longer
duration when compared with fly-by's." The approach to Uranus "has been studied in previous works" (its Refs.
5 and 6, both AIAA-2007 papers by the same authors); not treated here.

## Model (exactly what is computed)
- Model, quoted: "The tour is computed within the planar circular restricted three body model, based on which
  the stable and unstable manifolds associated with the libration points of that system are computed." The
  moons are "considered to move in circular, equatorial orbits. Therefore application of different, coupled
  Planar Restricted Three Body Problems (PCR3BP) is considered a valid approximation." All trajectories are
  referred to the Uranus equatorial plane. Equations (1)-(4): planar CR3BP equations of motion, effective
  potential, Jacobi constant C = -(xdot^2 + ydot^2) + 2 Omega, and "E = -C/2".
- Each moon defines its own Uranus-moon three-body system; the systems are "coupled" only by a coordinate
  transformation: "The manifolds associated with the libration points of each moon are computed in the
  relative synodic frame and subsequently translated and scaled to the Uranus - Oberon system, which is chosen
  as the main system of reference for the tour construction. The transformation between the two systems takes
  into account the initial phase of the moons and the associated non-autonomous phase difference during the
  entire transfer. ... In the main reference frame the manifolds of Oberon are time independent, whereas the
  manifolds associated with the other moons are time dependent (periodic)."
- At no time is more than one moon acting on the spacecraft in a single integration. Each arc is integrated
  in one Uranus-moon PCR3BP (the moon of that stage), with a thrust acceleration added during powered arcs.
- Computed objects per system: the collinear points L1 and L2 ("L1 and L2 are studied in this work") and the
  one-dimensional stable and unstable manifolds of their linearisation (Eqs. 7-8: x0 +/- d v, eigenvector
  perturbation propagated backward for stable, forward for unstable; "two legs for each manifold"). The
  manifold of L2 (stable) is used for capture, the manifold of L1 (unstable) for departure. Energy level used:
  "The lowest energy value that permits a transit from the outer realm, traversing the moon and going to the
  inner realm is the energy value associated with the second libration point L2, therefore the manifolds
  computed in this study are based on this libration point, thus energy level."
- The perturbation size d is "an arbitrary small value that assures the manifolds existence for all values of
  mu considered in the tour"; its numerical value is **not stated.** Quoted effect: "it seems to affect only
  the number of closed orbits around the moons, which increase with the diminution of d, and the minimum
  altitude above them. For some values of the perturbation even impact trajectories are possible".
- Jacobi constants or energies per moon at the L2 level: **not stated** in text or tables. Figure 2 (Hill
  regions) is "produced for an energy value that assures an opening being present even for the smallest mu,
  Miranda"; the value is not printed. Figure 10 plots C and E along the whole tour in the Uranus-Oberon frame
  (C roughly 3 at the start rising to roughly 5.5 at the end, E roughly -1.5 to -2.7, read by eye from the
  figure; no numbers printed).
- Periodic orbits: none computed. L1 and L2 are fixed points only. Future work (quoted): "Further development
  of the present work could deal with the use of the manifolds associated with the periodic orbits around the
  libration points. This provides more freedom in the intersection of the different manifolds and also the
  possibility to control the size and the number of the closed orbits around the moons."
- Heteroclinic connection (quoted): "the stable manifold associated with the second libration point of the
  specific moon, when propagated for a ballistic time greater than tmani, performs various closed orbits around
  that moon after which it passes onto the unstable manifold of L1 of the same moon. This transition is called
  a heteroclinic connection and is used to obtain the starting conditions for the next passage." Figures 8-9
  are captioned "Ballistic Captures and Escapes", with text: "The ballistic arcs that correspond to heteroclinic
  connections between the stable manifold of L2 and the unstable one of L1 of the same moon are shown in Fig.
  8." These are connections of one moon's own L2 stable and L1 unstable manifolds (the "closed orbits around
  the moon" are temporary-capture arcs, not computed periodic orbits).
- Resonant orbits, quasi-periodic orbits, tori, bicircular or four-body (two-moon) models: none computed (see
  the explicit statements below).

### Table 1, exact transcription (page 3; checked against the rendered page)
Caption as printed: "Table 1. Uranus Moons characteristics". Inclinations "are referred to the Uranus
equatorial plane".

| Uranus Moon | Mass [kg] | Semi-major Axis [km] | Eccentricity | Inclination [deg] |
|---|---|---|---|---|
| Oberon | 3.014 x 10^21 | 583519 | 0.0016 | 0.700 |
| Titania | 3.526 x 10^21 | 435910 | 0.0011 | 0.340 |
| Umbriel | 1.200 x 10^21 | 266300 | 0.0039 | 0.205 |
| Ariel | 1.350 x 10^21 | 190020 | 0.0012 | 0.260 |
| Miranda | 6.590 x 10^19 | 129390 | 0.0013 | 4.232 |

Source of Table 1 as cited: Ref. 7 (Seidelmann, Explanatory Supplement to the Astronomical Almanac, 2006).

### Table 2, exact transcription (page 4; checked against the rendered page)
Caption as printed: "Table 2. Uranus-moon identification parameters".

| Primaries | mu | Distance Unit [km] | Time Unit [sec] | Mass Unit [kg] |
|---|---|---|---|---|
| Uranus-Oberon | 3.4792 x 10^-5 | 583519 | 1.8539 x 10^5 | 8.6628 x 10^25 |
| Uranus-Titania | 4.0703 x 10^-5 | 435910 | 1.1970 x 10^5 | 8.6629 x 10^25 |
| Uranus-Umbriel | 1.3853 x 10^-5 | 266300 | 5.7157 x 10^4 | 8.6626 x 10^25 |
| Uranus-Ariel | 1.5584 x 10^-5 | 190020 | 3.4452 x 10^4 | 8.6626 x 10^25 |
| Uranus-Miranda | 7.6075 x 10^-7 | 129390 | 1.9358 x 10^4 | 8.6625 x 10^25 |

Text: "Due to the low mass parameters of the different systems [Table 2], no intersections in the position space
are present (Fig. 3). Thus, the only approaches feasible to perform the tour are either using a multi-burn
strategy or a continuous thrust, which modifies the spacecraft energy during propelled arcs."

### Table 3, exact transcription (page 8; checked against the rendered page)
Caption as printed: "Table 3. Fixed Parameters" (spacecraft "taken from previous studies").

| Parameter | Value |
|---|---|
| Initial Mass [kg] | 500 |
| Specific Impulse [sec] | 3200 |
| Thruster Power [W] | 1000 |
| Thrust Efficiency [-] | 0.5 |

Thrust modulus (Eq. 11): |Thr| = 2 eta P / (Isp g0); the value of thrust in newtons is **not stated.** (My
arithmetic from Table 3 and Eq. 11: about 31.9 mN.)

### Per-moon summary (what the paper prints for each system)
| System | mu (Table 2) | Objects computed | Stage in tour |
|---|---|---|---|
| Uranus-Oberon | 3.4792e-5 | L1 and L2 stable/unstable manifolds; main reference frame; ballistic capture shown in Fig. 5 | first: start fixed "on the intersection of the x-axis of the stable L2 manifolds of Oberon", t0 = 0 |
| Uranus-Titania | 4.0703e-5 | same; ballistic capture/escape in Figs. 8-9 | second; Oberon to Titania propulsion about 11.5 days |
| Uranus-Umbriel | 1.3853e-5 | same | third; propulsion arc to Umbriel about 150 days |
| Uranus-Ariel | 1.5584e-5 | same | fourth; propulsion duration not stated |
| Uranus-Miranda | 7.6075e-7 | same; "the PCR3BP approximation for Miranda does not hold as well due to its higher inclination (Table 1)" | fifth and last; propulsion arc to Miranda about 151 days |

Libration points L3, L4, L5 are mentioned only as part of the generic CR3BP description; no object is computed
at them. Only the planetary distance is scaled in Fig. 3 ("the manifolds are shown as seen each in their
rotating reference frame").

## How transfers between moons are made
- Method, quoted: "The connection of two states on an unstable and stable manifold, respectively, enables the
  passage from one system to the other. The connection of two states requires an energy change provided by the
  Electric Propulsion, where the thrust direction and modulus are two parameters of optimization." Each
  transfer is an electric-propulsion arc (powered, in the Uranus-Oberon rotating frame) from the unstable L1
  manifold of the previous moon to the initial state of the stable L2 manifold of the next moon, followed by a
  ballistic stable-manifold arc and capture-like loops around the moon.
- In-plane thrust: acceleration a cos(alpha) T + a sin(alpha) N, with scaling tau in [0,1], mass-flow equation
  in Eq. 10; "The thrust angle alpha is measured counterclockwise from the velocity direction." Initial guess
  alpha "always near 180 because the tour is computed from the outer to the inner moon and an anti-tangential
  thrust (in a synodic frame) assures a diminution of the spacecraft energy."
- Control vector per passage (Eq. 12): u = {t0, tel, tmani, alpha, tau, theta}: exit time from the previous
  moon's L1 unstable manifold, electric thrust duration, duration on the stable L2 manifold, thrust angle,
  thrust scaling, and the initial angle of the next moon. Thrust law discretised on N mesh points, so u has
  dimension 2N + 4; N = 10 ("only doubled for long propulsion phases").
- Optimisation: nonlinear programming (sequential quadratic programming, quasi-Newton Hessian update) with
  equality constraints on the final state, minimising propellant mass; then a simplex refinement minimising
  phase-space distance to the insertion state, with "a tolerance of 5% ... on the phase-distance of the
  conjunction states". Text: "the chaotic dynamics of the model lead to completely different solutions even
  for very similar initial conditions."
- Chaining: the heteroclinic L2-stable to L1-unstable ballistic arc of a moon provides the t0 of the next
  passage.
- Propulsion durations printed: Oberon to Titania "approximately 11.5 days"; passages to Umbriel and Miranda
  "approximately 150 and 151 days respectively". The Titania-to-Umbriel and Umbriel-to-Ariel and Ariel-to-
  Miranda split is not tabulated; the sentence names only "Umbriel and Miranda". Delta-V, propellant or
  duration per passage: **not tabulated** (no table of per-passage results). Fig. 6 shows green powered arcs
  after the Titania stage, around Umbriel to Ariel, and ending at Miranda (read by eye).
- Time spent in each moon's ballistic capture: "for the tour studied the time spent around each moon ranges
  from several days to almost a month". Number of closed orbits per moon: not stated. Capture altitudes:
  "no control on the altitude of the capture orbit around the generic moon, therefore also impact trajectories
  can be present." Closing to a chosen orbit "can be obtained with two very small Delta-Vs, in the order of
  some [m/s]" (stated, not computed).
- Does the tour repeat or return? One-way. The trajectory "starts with a ballistic approach arc towards the
  outer moon, Oberon, and ends with a stable orbit around Uranus with a radius smaller than the inner moon
  Miranda's distance"; order Oberon, Titania, Umbriel, Ariel, Miranda, each visited once. "The total time
  required is approximately 930 dys to lead the spacecraft beyond the inner moon." Each moon is visited once;
  no return to an earlier moon is described.

## Totals as printed (Conclusion, page 12)
"The propellant mass fraction required for the entire tour is approximately 0.062879 that corresponds to a
Delta-V of 2.0387 [km/s]. This value is very small and can easily be included in the launch mass budget of
the spacecraft. The total time required is approximately 930 dys to lead the spacecraft beyond the inner
moon. This time must be added to the transfer time towards Uranus that remains the most stringent constraint
to the feasibility of this kind of missions."

- My arithmetic: 0.062879 x 500 kg = 31.4 kg propellant; Delta-V = Isp g0 ln(1/(1 - 0.062879)) = 2.039 km/s
  with Isp 3200 s, which matches the printed 2.0387 km/s. Fig. 6 mass panel runs from 500 kg to about 469 kg
  (read by eye).
- No chemical-propulsion delta-V total, no per-moon table, no Umbriel-specific delta-V.

## Everything the paper says about Umbriel (complete)
Umbriel never gets its own paragraph. Every occurrence is a table row, a figure label, or part of a list of the
five moons. The sentences that name it:
- "During the complete transfer several closed orbits around Oberon, Titania, Umbriel, Ariel and Miranda are
  executed." (Introduction, page 2)
- "On the other hand the passage to Umbriel and Miranda require very long propulsion arcs (approximately 150
  and 151 days respectively), as they are physically very different from the previous moon and with relatively
  large radial difference." (page 9)
- Table 1 row Umbriel: mass 1.200 x 10^21 kg, semi-major axis 266300 km, e 0.0039, inclination 0.205 deg
  (largest eccentricity of the five). Table 2 row Uranus-Umbriel: mu 1.3853 x 10^-5, distance unit 266300 km,
  time unit 5.7157 x 10^4 s, mass unit 8.6626 x 10^25 kg.
- Figures: Fig. 2 and Fig. 3 (Hill region and L1/L2 manifolds, labelled UL1, UL2), Fig. 5 and Fig. 7 (tour),
  Fig. 6 (position/velocity/mass versus time; Umbriel stage sits in the middle of the tour, roughly 450 to 600
  days read by eye), Figs. 8-9 ("Ballistic Capture [Uranus-Umbriel rotating frame]" with an L1 and L2 close-up).
- Umbriel is the third moon of five. It is neither first nor last, has no per-moon result, and is never
  paired with another named moon in text except in the arc-length sentence above. There is no Umbriel-Titania
  or Umbriel-Oberon connection computed other than the sequential chain Titania, Umbriel, Ariel in the tour.

## Text-search counts (whole words, case-insensitive, lines joined, ligatures normalised; includes figure labels and the reference list)
| Term | Paper A | Term | Paper A |
|---|---|---|---|
| cycler / cycling | 0 / 0 | Umbriel | 14 |
| periodic | 3 | Titania | 17 |
| repeating | 0 | Oberon | 29 |
| free-return | 0 | Ariel | 13 |
| resonant | 0 | Miranda | 17 |
| resonance | 1 | Triton | 0 |
| heteroclinic | 3 | Proteus | 0 |
| homoclinic | 0 | Nereid | 0 |
| quasi-periodic | 0 | torus / tori (whole word) | 0 / 0 |
| four-body | 0 | bicircular | 0 |

Quoted sentences for the dynamical terms:
- "periodic" (3): "the manifolds associated with the other moons are time dependent (periodic)" (describing the dependence
  of a non-Oberon moon's manifold on time in the Oberon frame, not a periodic orbit); "Further development of
  the present work could deal with the use of the manifolds associated with the periodic orbits around the
  libration points."; and the title of Ref. 11 ("Heteroclinic Connections between Periodic Orbits and Resonance
  Transitions in Celestial Mechanics").
- "resonance" (1): only in the title of Ref. 11.
- "heteroclinic" (3): "This transition is called a heteroclinic connection and is used to obtain the starting
  conditions for the next passage."; "The ballistic arcs that correspond to heteroclinic connections between
  the stable manifold of L2 and the unstable one of L1 of the same moon are shown in Fig. 8."; and the title of
  Ref. 11.
- Closest words to "four-body": "a complete propagation in a multi-body model, such as superposition of
  bi-circular models" (Conclusion, future work: "bi-circular" hyphenated, one occurrence; the word "bicircular"
  as one word does not appear).

## Explicit answers
- **Trajectory returning to the same two moons repeatedly:** no. The tour is one-way, Oberon, Titania, Umbriel,
  Ariel, Miranda, each once.
- **Resonant periodic orbit:** none. No m:n resonance with any moon appears; the only occurrence of the word
  "resonance" is a reference title.
- **Orbit or connection computed with two moons acting at once:** none. Each arc is in a single Uranus-moon
  PCR3BP. Two moons never act simultaneously. The paper suggests "superposition of bi-circular models" only as
  future work.
- **Quasi-periodic orbit or torus:** none.
- **Closed or periodic orbit around a moon, computed:** none as an object; "closed orbits" around each moon are
  temporary-capture loops of the L2 stable manifold, with no periodicity condition and no stability analysis.

## What it does NOT contain
- No periodic orbits (Lyapunov, halo, resonant or other) of any system; only L1 and L2 fixed points and their
  linear (eigenvector) manifolds.
- No Jacobi constant or energy values printed per moon.
- No per-passage delta-V, propellant, thrust-arc or duration table; only the Oberon-Titania 11.5 days, the
  Umbriel and Miranda "approximately 150 and 151 days", and the totals 0.062879 and 2.0387 km/s and 930 days.
- No thrust value in newtons (only power, efficiency and Isp).
- No value of the manifold perturbation d.
- No resonant, quasi-periodic, torus, four-body or two-moon model; no spatial (3D) model (the Conclusion lists
  "a three dimensional approach" as future work).
- No tour that visits a moon more than once; no cycler, free-return or repeating trajectory.
- No treatment of the Uranus approach leg (deferred to the authors' earlier papers, Refs. 5-6).
- No handling of Miranda's 4.232 deg inclination beyond the remark that the PCR3BP "does not hold as well".
- No ephemeris, no perturbations (no J2 of Uranus, no non-coplanar effects).

## Relevant bibliography (as printed in this paper; all 12 references)
1. S. D. Ross, W. S. Koon, M. W. Lo and J. E. Marsden, "Design of a multi-moon orbiter", 13th AAS/AIAA Space
   Flight Mechanics Meeting, Ponce, Puerto Rico, Paper No. AAS 03-143, 2003.
2. S. D. Ross, W. S. Koon, M. W. Lo and J. E. Marsden, "Constructing a Low Energy Transfer Between Jovian
   Moons", Contemporary Mathematics 292, 129-145, 2002.
3. G. Gomez, W. S. Koon, M. W. Lo, J. E. Marsden, J. Masdemont and S. D. Ross, "Invariant Manifolds, the
   Spatial Three-Body Problem and Space Mission Design", Advances in the Astronautical Sciences, volume 109,
   part 1, p. 3-22, AAS 01-301, 2001.
4. F. Topputo, "Low-Thrust Non-Keplerian Orbits: Analysis, Design and Control", Ph.D Thesis, 2007.
5. C. Casaregola, K. Geurts, P. Pergola and M. Andrenucci, "Exploitation of Three-Body Dynamics by Electric
   Propulsion for Outer Planetary Missions", AIAA-2007-5228, 43th Joint Propulsion Conference, Cincinnati,
   Ohio, July 2007.
6. C. Casaregola, K. Geurts, P. Pergola and M. Andrenucci, "Radioisotope Low-Power Electric Propulsion
   Missions to the Outer Planets", AIAA-2007-5234, 43th Joint Propulsion Conference, Cincinnati, Ohio, July
   2007.
7. P. K. Seidelmann, "Explanatory Supplement to the Astronomic Almanac", University Science Books, California,
   2006.
8. E. A. Belbruno and J. K. Miller, "Sun-Perturbed Earth-to-Moon Transfers with Ballistic Capture", Journal of
   Guidance, Control and Dynamics, Vol. 16, No. 4, pp. 770-775, 1993.
9. V. Szebehely, "Theory of Orbits, The Restricted Problem of Three Bodies", Academic Press Inc., New York,
   1967.
10. F. B. Zazzera, F. Topputo, M. Massari, "Assessment of Mission design Including Utilization of Libration
    Points and Weak Stability Boundaries", Ariadna Study id: 03/4103.
11. S. D. Ross, W. S. Koon, M. W. Lo and J. E. Marsden, "Heteroclinic Connections between Periodic Orbits and
    Resonance Transitions in Celestial Mechanics", Chaos 10(2), 427-469, 2000.
12. J. T. Betts, "Survey of Numerical Methods for Trajectory Optimization", Journal of Guidance, Control and
    Dynamics, Vol. 21, No. 2, pp. 193-207, 1998.

Refs. 1, 2, 3 and 11 are the multi-moon-orbiter and Jovian-tour three-body papers on which the approach is
modelled (the Introduction says "An approach similar to the mentioned Jovian tour is implemented, however, in
this study Electric Propulsion is considered for the execution of the required manifold transitions", and "A
similar strategy was investigated by Topputo, developing a low thrust assisted version of the Jovian tour").

## Suggested topology label
**halo** is the nearest of the project's labels (libration-point manifold transfers between systems), with a
caveat: the paper computes no periodic or halo orbits, only the L1 and L2 manifolds. Justification: "Electric
Propulsion is used to link the low energy paths between stable and unstable invariant manifolds computed for
the different Uranus - moon systems". It is not "repeated-moon" or "resonant" (each moon visited once, no
resonance), not "pump-tour" and not "mga-tour" (no flybys; powered manifold transfers with capture loops).

## Suggested corpus index description
Pergola et al. 2007 (IEPC-2007-305): one-way electric-propulsion tour of the five main Uranian moons
(Oberon, Titania, Umbriel, Ariel, Miranda) built from planar CR3BP L1/L2 manifolds of each Uranus-moon system;
2.0387 km/s, about 930 days; no periodic, resonant or two-moon orbits.
