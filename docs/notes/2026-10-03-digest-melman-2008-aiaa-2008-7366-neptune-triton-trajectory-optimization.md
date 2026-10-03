# Digest — Melman, Orlando, Safipour, Mooij, Noomen (2008), "Trajectory Optimization for a Mission to Neptune and Triton" (AIAA 2008-7366)

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, 21 pages; Tables 3-5, Figures 12-17 and the
Triton-transfer pages 16-20 checked against the rendered pages as well as the text layer). **Purpose:** the
project is about to work on Neptune-Triton resonant-orbit objects and must know exactly what in-system
three-body trajectory work at Neptune is published. About two-thirds of this paper is interplanetary cruise
and Neptune capture; the Neptune-system phase is the last third. This note is facts only; no novelty verdict is
written here. Where a number is not printed it is marked "not stated"; arithmetic of mine is labelled as such.

## Citation
J. Melman, G. Orlando, E. Safipour, E. Mooij, R. Noomen, "Trajectory Optimization for a Mission to Neptune and
Triton", AIAA 2008-7366, AIAA/AAS Astrodynamics Specialist Conference and Exhibit, 18-21 August 2008,
Honolulu, Hawaii.

- Affiliation as printed: Delft University of Technology, Faculty of Aerospace Engineering, Kluyverweg 1,
  2629 HS Delft, The Netherlands (Melman, PhD student; Orlando, MSc student; Safipour, space software engineer,
  Logica; Mooij and Noomen, assistant professors). Copyright 2008 by Delft University of Technology.
- Filed in the private paper corpus as
  `melman-orlando-safipour-mooij-noomen-2008-trajectory-optimization-mission-neptune-triton-aiaa-2008-7366-doi-10.2514-6.2008-7366.pdf`,
  md5 `d5e8d88926593825ae3a14dabd0e414c`.

## What the paper is
Abstract (quoted): "This paper treats a mission to Neptune and Triton, focusing on the design and optimization
of its trajectory. All trajectory parts are dealt with: the interplanetary cruise to reach Neptune, capture at
Neptune, and the reconnaissance of the Neptunian system, including both a Neptune and a Triton orbiter. The
spacecraft uses five swing-bys to arrive at Neptune: Venus twice, Earth, Jupiter and Saturn. At Neptune an
aerocapture takes place. At the first apoapsis after aerocapture, the two orbiters separate from each other.
The Neptune orbiter goes on a trajectory that results in up to 160 flybys of Neptune's inner moons within one
month. This trajectory has been determined by a multi-objective optimization implementing Genetic Algorithms.
The Triton orbiter performs a transfer that has been determined using the Circular Restricted Three-Body
Problem. It enters a near-polar orbit around Triton, which ensures good coverage. Finally, two orbiters of
about 600 kg will orbit both Neptune and Triton."

Sections as printed: I Introduction; II Models (A Interplanetary Trajectory; B Orbit Insertion at Neptune; C
Transfer to Triton; D Optimization); III Results (A Interplanetary; B Orbit Insertion at Neptune; C Transfer
to Triton); IV Conclusions and Recommendations; References (22). Five tables, 17 figures. Three MSc-thesis
projects ("performed in the framework of three MSc thesis projects"); the Triton part is Orlando's (Ref. 14,
Master's thesis, Delft, September 2008).

## Interplanetary phase (summary only; not the subject of this digest)
Patched conics with Gooding Lambert; 77 swing-by sequences (Table 1); best VVEJS, launch 2012/05/25, Neptune
arrival 2029/12/18, total Delta-V 5.819 km/s, time of flight 17.6 years; Neptune arrival C3 = 149.823 km^2/s^2,
periapsis velocity 25.033 km/s at 4000 km (Table 2). Neptune delivered mass 4092 kg. (My arithmetic: C3 gives a
hyperbolic excess speed of about 12.24 km/s.) Launcher Atlas V 551 with Star 48V; Isp 320 s.

## The Neptune-system phase

### Capture options and what is said about moons (section II-B and III-B)
- Capture by moon gravity assist is rejected (quoted): "An examination of the gravity-assisted capture strategy
  learns that the moons of Neptune are not massive enough to initiate capture. The maximum achievable change in
  energy can be obtained by using Triton, since it is Neptune's most massive moon. The orbital energy of the
  spacecraft arriving at Neptune is about 75 km^2 s^2 ... The achievable change in energy due to Triton is just
  0.7515 km^2 s^2, which is negligible with respect to the spacecraft energy. Therefore, it is concluded that
  the moons of Neptune will not be used for capturing into an orbit around Neptune." (units as printed.)
- Strategies analysed: chemical capture (Eq. 3, Isp 320 s) and aerocapture (Galileo-probe geometry, lift
  coefficient 1.4, 22 g and 16 kW/cm^2 limits, peak 3300 W/cm^2, 710 s passage, TPS 1119 kg, deceleration under
  3.85 g). Aerocapture is the preferred strategy ("because of its high performance when compared to the chemical
  capture").
- Triton is the destination of a second orbiter; the Neptune orbiter is the one that flies the inner moons.
  Triton's inclination is printed as 157.3 deg ("The Triton orbiter targets a trajectory with the same
  inclination as Triton (157.3 deg)"). The aerocapture entry is at 157.3 deg and "The aerocapture maneuver
  targets an apoapsis radius that is equal to the semi-major axis of Triton's orbit". Triton's semi-major axis in
  km: not stated.
- Neptune orbiter science orbit: inclination 60 to 120 deg, periapsis 26164 to 29764 km, apoapsis 118000 to
  130000 km. "the inner moon furthest away from Neptune, Proteus, has a circular orbit with a radius of 118000
  km". Inner moons (Fig. 8 and Fig. 10 labels): Naiad, Despina, Thalassa, Galatea, Larissa, Proteus, with
  "Inner moons (except Triton) 41,000 - 117,000 km". Definition of a flyby: "a flyby takes place once the
  distance between the spacecraft and a moon is less than 15000 km". Perturbations included in the Neptune
  analysis: "the J2-term of Neptune's gravity, the third-body effect of Triton, and the aerodynamic forces once
  in Neptune's atmosphere" (those above 1e-4 m/s^2).
- Flyby count is a number of encounters within 15000 km in one month of observation, optimised with a
  multi-objective genetic algorithm (GAVaPS; 5000 initial individuals, 50 generations). No individual flyby
  (moon, date, distance, V-infinity) is listed; no moon sequence, resonance with any moon, periodic orbit or
  repeating pattern is described for the Neptune orbiter. The 15000 km flyby is a distance criterion, not a
  gravity-assist design.

### Table 3, exact transcription (page 16; checked against the rendered page)
Caption as printed: "Table 3. Comparison of chemical and aerocapture strategies. For available mass and number
of flybys a range is indicated, which is intrinsic to the nature of multi-objective optimization."

| Option | Total Neptune Orbiter's Mass except TPS | Payload Mass (total minus fuel and TPS) | Maximum number of flybys (observation period of 1 month) | Mission complexity and risks |
|---|---|---|---|---|
| Chemical capture | 2046 kg | 347-383 kg | 44-176 | Standard |
| Aerocapture | 1486.5 kg * | 581-964 kg | 41-161 | High |

Footnote as printed: "* Half of the mass that results from subtracting the TPS mass (1119 kg) from the
spacecraft mass that arrives at the SOI of Neptune (4092 kg)."

Consistency notes from the text: chemical case "the maximum number of flybys increased from 145 in the first
generation to 176 in the 50th generation"; aerocapture "The number of flybys ranges from 41 to 161 during one
month of observation time". The abstract and Conclusion say "up to 160 flybys" and "performs 160 flybys of the
inner moons"; Table 5 says "161 moon flybys". The paper prints 160 and 161 in different places (not
reconciled). Aerocapture delta-V range for the Neptune orbiter: "1.36 to 2.95 km/s" (inclination change and
periapsis raise, with an extra 0.2 km/s contingency stated as included); chemical total "5.26 to 5.59 km/s".

### The Triton transfer (section II-C and III-C)

**Model, quoted:** "The transfer to Triton is designed by making use of low-energy capture at Triton in the
Neptune-Triton Circular Restricted Three-Body Problem. The stability and instability properties near the
collinear Lagrangian points become thus crucial in the process of the trajectory design." and "To solve the
problem of low-energy transfer and temporary capture, a method has been implemented that involves a new
technique called Energy Surface Reconstruction (ESR). The method is based on nested manifolds, with the main
difference to consider several values of the Jacobi Energy and to allow building 'windows' of viable capture
trajectories." The model is spatial: the Conclusion says "To capture into a high-inclination orbit, it proved
to be very useful to look at the Circular Restricted Three-Body Problem in a three-dimensional way." Frame:
dimensionless Neptune-Triton Rotating Reference Frame (NTRRF).

**The eight-step ESR method, as printed (condensed, all wording from the paper):**
1. Backwards integration of near-circular polar orbits (an orbit phi, larger semi-major axis than the target,
   hence unstable) until a section Sigma at the L2 point, parallel to the yz-plane. The set of escape
   trajectories, excluding impacts, is S = W(phi)_{C,L2} intersect Sigma ("subscript C stands for capture
   (forwards in time) and L2 indicates that such trajectories pass in the neighborhood of the L2 point").
2. Pattern recognition on the surface S: "ordered structures and disordered regions"; points in ordered
   structures "remain in close proximity to each other after terminating on the osculating capture trajectory."
3. Selection of dense subsets of interest, called sigma.
4. Enlargement of sigma through ESR: neighbouring points on Sigma are assigned a velocity with the same Jacobi
   energy as the original point (partial-derivative construction for ydot and zdot, then xdot from the Jacobi
   constant); the enlarged set is sigma-tilde; "All the points on this set lead to temporary capture when
   integrated forward in time."
5. Further backwards integration to a more suitable surface Sigma2, "close to the xy-plane again" (Triton's
   orbital plane), widening the target window.
6. Bounds of the window by Overlapping Convex Hulls (OCH).
7. Design of a trajectory targeting a point in the window.
8. Matching velocity at that point by a further ESR; the Delta-V is the velocity difference.

**Relation to periodic orbits, quoted:** "The method presents analogies with the methods envisaged in
literature, where temporary capture trajectories are designed as intersections between two stable manifolds:
the stable manifold of the Lyapunov orbit at L2 and the stable manifold of an unstable periodic orbit around
the planetary satellite." and "It is actually possible that each of these highly perturbed Keplerian orbits is a
perturbed orbit originating from a local unstable periodic orbit of the CR3BP". The paper itself computes no
Lyapunov, halo or other periodic orbit and no invariant manifold of one; it integrates escape trajectories of a
polar Keplerian-like orbit backwards. Future work, quoted: "Further research can also focus on the possibility
to use periodic orbits of the Circular Restricted 3-Body Problem and their manifolds to find the optimal set of
maneuvers leading to the final target orbit."

**Inputs and numbers printed (Section III-C):**
- Target orbit, quoted: "a nearly circular target orbit around Triton with a semi-major axis of 1900 km and an
  inclination of 70 deg, in analogy to the target orbit used in the work by Russell" (Ref. 12). "Such an orbit
  requires a maintenance Delta-V of 31 m/s for a mission time of 6 months". (The 1900 km is a semi-major axis,
  not an altitude.)
- Initial conditions: "the exiting conditions from the aerocapture and lie in Triton's orbital plane".
- Orbit phi: "set equal to a nearly circular polar orbit with a semi-major axis of 8000 km"; "For a circular
  polar orbit with a semi-major axis of 8000 km, the time to escape is in the order of 10 days." Two thousand
  escape trajectories "for different values of the initial true anomaly".
- Sigma is x = x_L2; Sigma2 is "the set of points x = 0.25125 in the dimensionless Neptune-Triton Rotating
  Reference Frame". Window: "an about 2500 km wide and 100 km high window in the configuration space"
  (Fig. 14 axes: y about 4.39 to 4.415 x 10^5 km, z about +/-125 km, read by eye).
- Optimisation: genetic algorithm (35 generations, 300 individuals), local Monte Carlo, BFGS quasi-Newton
  refinement with mixed quadratic and cubic line search.
- Six manoeuvres, quoted: Delta-V1 "raises the periapsis of the starting orbit in order to avoid a possible
  reentry into the atmosphere of Neptune and inserts the spacecraft into a phasing orbit"; Delta-V2 "changes
  the phasing orbit to target the window on the Sigma2 surface in the NTRRF"; Delta-V3 "inserts the spacecraft
  into a temporarily captured trajectory by use of the ESR method"; Delta-V4 "is executed after the passage of
  the L2 point and changes the orbital parameters in order to approach the ones from the target orbit";
  Delta-V5 "lowers the altitude of the phi orbit in order to move on an eccentric orbit whose pericenter
  coincides with the pericenter of the target orbit"; Delta-V6 "circularizes the transfer orbit resulting from
  the previous step". Arcs: 1-2 phasing (green), 2-3 aim at the window (red), 3-4 temporary capture (blue),
  4-5 "orbits Triton on an unstable near-polar orbit (black)", 5-6 final transfer to the target orbit.
- Total, quoted: "The total Delta-V (including all six Delta-Vs mentioned before) is 2850 m/s." Maximum
  excursion from Triton's plane "about 2 x 10^4 km". "the mass delivered at Triton is me = 599.6 kg" from
  m0 = 1486.5 kg and Isp = 320 s.

### Table 4, exact transcription (page 19; checked against the rendered page)
Caption as printed: "Table 4. Delta-Vs for capture trajectory in m/s."

| Delta-V1 | Delta-V2 | Delta-V3 | Delta-V4 | Delta-V5 | Delta-V6 |
|---|---|---|---|---|---|
| 1229 | 34 | 1251 | 20 | 31 | 284 |

My arithmetic: the six values sum to 2849 m/s (printed total 2850 m/s; rounding). Delta-V1 plus Delta-V3 is
2480 m/s, about 87% of the total (the paper says "the main contributions are due to Delta-V1 and Delta-V3").

### Table 5, exact transcription (page 20; checked against the rendered page)
Caption as printed: "Table 5. Overview of the entire mission to Neptune and Triton. The numbers for the mass
indicate the spacecraft mass after the mentioned mission part has taken place."

| Mission Part | Mass [kg] | Delta-V [km/s] | Note |
|---|---|---|---|
| Launch | 4706 | 1.58 | C3 = 18.94 km^2/s^2 |
| Swing-bys | 4092 | 0.439 | VVEJS |
| Aerocapture | 2973 * | - | TPS mass subtracted |
| Neptune orbit insertion | 581 ** | 2.95 | 161 moon flybys |
| Triton orbit insertion | 600 | 2.85 | target orbit |

Footnotes as printed: "* This mass will be divided evenly over the two orbiters. ** This corresponds to the
solution on the Pareto front with the largest number of moon flybys."

## Verification of the project's records
- "Circular restricted three-body treatment with energy-surface reconstruction": confirmed (the Neptune-Triton
  CR3BP, spatial, with the ESR method; "Energy Surface Reconstruction (ESR)" is named by the authors).
- "6 manoeuvres": confirmed (Delta-V1 to Delta-V6, Table 4).
- "2850 m/s": confirmed as the printed total (the six table values sum to 2849 m/s).
- "1900 km 70-degree orbit": confirmed, with a refinement: 1900 km is the semi-major axis of a "nearly
  circular" target orbit, inclination 70 deg, taken "in analogy to the target orbit used in the work by Russell"
  (Ref. 12, a Europa paper), with a 31 m/s per 6 months maintenance cost.
- Not in the records but printed: the capture is a temporary capture via L2 passage (arc 3-4 "passes the L2
  point"), the starting orbit phi is an 8000 km polar orbit integrated backwards, the transfer is powered by
  chemical impulses (Isp 320 s), and only the Neptune-Triton pair is a three-body system in the paper.

## Text-search counts (whole words, case-insensitive, lines joined, ligatures normalised; includes figure labels and the reference list)
| Term | Paper B | Term | Paper B |
|---|---|---|---|
| cycler / cycling | 0 / 0 | Umbriel | 0 |
| periodic | 5 | Titania | 0 |
| repeating | 1 | Oberon | 0 |
| free-return | 0 | Ariel | 0 |
| resonant | 0 | Miranda | 0 |
| resonance | 0 | Triton | 70 |
| heteroclinic | 0 | Proteus | 5 |
| homoclinic | 0 | Nereid | 0 |
| quasi-periodic | 0 | torus / tori (whole word) | 0 / 0 |
| four-body | 0 | bicircular | 0 |

Quoted sentences for the dynamical terms:
- "periodic" (5): "the stable manifold of the Lyapunov orbit at L2 and the stable manifold of an unstable periodic orbit
  around the planetary satellite."; "It is actually possible that each of these highly perturbed Keplerian
  orbits is a perturbed orbit originating from a local unstable periodic orbit of the CR3BP"; "In fact, the
  unstable periodic orbits applied in one of these papers are highly perturbed polar near-circular orbits.";
  "Further research can also focus on the possibility to use periodic orbits of the Circular Restricted 3-Body
  Problem and their manifolds"; and the title of Ref. 12 ("Designing Ephemeris Capture Trajectories at Europa
  Using Unstable Periodic Orbits").
- "repeating" (1): "Reconstructing the neighborhood of a generic point i consists of repeating such a process for
  different values of delta-y and delta-z." (a computational step, not a trajectory property).
- Proteus (5) appears only in moon lists and figure labels (Fig. 8, Fig. 10, and the Proteus orbit-radius
  sentence). "Lyapunov" (1) appears only in the first sentence quoted above.

## Explicit answers
- **Trajectory returning to the same two moons repeatedly:** none. The Triton orbiter flies one transfer to one
  target orbit; the Neptune orbiter's inner-moon flybys are counted, not listed, and no two-moon pattern is
  described.
- **Resonant periodic orbit:** none; no resonance (with Triton or any other moon) is mentioned.
- **Orbit or connection computed with two moons acting at once:** none. The Neptune-orbiter analysis includes
  Triton as a perturbation (the third-body effect) alongside J2 and drag. The Triton transfer is a
  Neptune-Triton CR3BP with no other moon.
- **Quasi-periodic orbit, torus, heteroclinic or homoclinic connection:** none.
- **V-infinity at Triton or at any Neptunian moon:** not stated. The only printed arrival numbers are at
  Neptune (C3 149.823 km^2/s^2; periapsis velocity 25.033 km/s).
- **Flyby sequence of Triton:** none; Triton is the capture target, not flown by. The only Triton gravity numbers
  are the 0.7515 km^2/s^2 energy change (rejected capture) and the temporary-capture arc 3-4.
- **Use of other Neptunian moons:** only as the Neptune orbiter's counted inner-moon flybys (Naiad, Thalassa,
  Despina, Galatea, Larissa, Proteus, 41000 to 117000 km) at a 15000 km criterion.

## What it does NOT contain
- No periodic orbit (Lyapunov, halo, resonant, other) is computed; no manifold of one.
- No resonance with Triton or any moon; no V-infinity values at any moon; no flyby list for the Neptune orbiter.
- No Jacobi constant or mass parameter of the Neptune-Triton system is printed (the value of mu is not stated,
  nor the Jacobi energy of the capture).
- No second moon in any three-body model; no four-body or bicircular model.
- No repeating or cycling trajectory of any kind; one-way mission.
- No ephemeris-model validation of the Triton capture (the transfer is in the CR3BP; the perturbations listed are
  for the Neptune orbiter).
- No numerical value of Triton's semi-major axis, period or mass in text.
- Nereid and the irregular moons are not mentioned.

## Relevant bibliography (as printed in this paper; the references bearing on the Neptune-system phase, with the rest listed for completeness)
- 1. W. Carroll and A. Ostlie, "An Introduction to Modern Astrophysics", Addison-Wesley, Reading, 1996 (moon data).
- 2. M. Noca and R. Bailey, "Mission Trades for Aerocapture at Neptune", AIAA-2004-3843, 40th
  AIAA/ASME/SAE/ASEE Joint Propulsion Conference and Exhibit, Fort Lauderdale, Florida, July 11-14, 2004.
- 11. S. Campagnola and M. Lo, "BepiColombo Gravitational Capture and the Elliptic Restricted Three-Body Problem",
  Proceedings in Applied Mathematics and Mechanics, November 2007.
- 12. R. Russell and T. Lam, "Designing Ephemeris Capture Trajectories at Europa Using Unstable Periodic Orbits",
  Journal of Guidance, Control, and Dynamics, Vol. 30, No. 2, March-April 2007.
- 13. G. Orlando and R. Noomen, "Temporary Capture Method based on Energy Surface Reconstruction, Manifolds
  Dynamics and Patterns Recognition", to be published.
- 14. G. Orlando, "Trajectory Optimization for a Mission to Neptune and Triton: release of a Triton orbiter and
  capture at Triton", Master's thesis, Delft University of Technology, The Netherlands, September 2008.
- 15. V. Szebehely, "The Theory of Orbits", Academic Press, 1967.
- 22. G. Orlando, E. Mooij and R. Noomen, "Optimal Orbital Stability around Planetary Satellites as Optimization
  Problem", to be published.
- Other references (interplanetary, aerocapture and optimiser methods): 3 Gooding 1990 (Lambert); 4 Cornelisse,
  Schoyer and Wakker, Rocket Propulsion and Spaceflight Dynamics, 1979; 5 Sponnick and Jensen, Atlas Launch System
  Mission Planner's Guide, 2004; 6 Larson and Wertz, SMAD, 3rd ed., 1999; 7 Miner and Wessen, Neptune: The Planet,
  Rings and Satellites, 2002; 8 Vasile, Ann. N. Y. Acad. Sci. 1065, 2005; 9 Regan and Anandakrishnan, Dynamics of
  Atmospheric Re-entry, 1993; 10 Lockwood, "Neptune Aerocapture Systems Analysis", AIAA-2004-4951; 16 Goldberg,
  1989; 17 Michalewicz, 1996; 18-21 Broyden 1970, Fletcher 1970, Goldfarb 1970, Shanno 1970.

## Suggested topology label
**none.** Justification: the Neptune-system content is a one-way, powered transfer into a temporary capture at
Triton through the L2 neighbourhood ("The Triton orbiter performs a transfer that has been determined using the
Circular Restricted Three-Body Problem"), with no libration-point orbit, no resonance and no repeated
encounter; the inner-moon flybys are counted, not designed. (The nearest label, "halo", does not fit: no
libration-point orbit or its manifold is computed.)

## Suggested corpus index description
Melman et al. 2008 (AIAA 2008-7366): Neptune and Triton mission trajectory study (VVEJS cruise, aerocapture);
Triton low-energy capture by spatial Neptune-Triton CR3BP "Energy Surface Reconstruction", 6 impulses,
2850 m/s into a 1900 km, 70 deg orbit; no periodic, resonant or multi-moon orbits.
