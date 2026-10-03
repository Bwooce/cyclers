# Digest — Strange, Landau, Longuski (2013), "Design of Initial Inclination Reduction Sequence for Uranian Gravity-Assist Tours" (AAS 13-801)

**Digested:** 2026-10-03 (17 pages; text layer has doubled (overprinted) text in the abstract and in some
figure labels, so page 1 and pages 13 to 15 (Figures 8 to 11) were viewed as rendered images, and Figure 8
was re-rendered at 250 dpi to read its labels; the equations were read from the text layer only). **Purpose:**
the Uranian moon-tour literature check for the six catalogued Uranian two-moon quasi-cyclers. This paper is
Ref. [10] of Landau, Davis & Karimi 2023 (AAS 23-460). This note is facts only; no novelty verdict is written
here. Statements marked "by my arithmetic" or "by eye" are mine; everything else is quoted or paraphrased.

## Citation
Nathan J. Strange, Damon F. Landau and James M. Longuski, "DESIGN OF INITIAL INCLINATION REDUCTION SEQUENCE
FOR URANIAN GRAVITY-ASSIST TOURS".

- Header line as printed (page 1): "AAS 13-801". Printed page numbers run 1469 to 1485.
- Affiliations as printed (page-1 footnotes): Strange, Systems Engineer, Jet Propulsion Laboratory,
  California Institute of Technology, 4800 Oak Grove Dr., Pasadena, CA, and Ph.D. Candidate, Purdue
  University, School of Aeronautics and Astronautics, 701 W. Stadium Ave. West Lafayette, IN; Landau,
  Systems Engineer, Jet Propulsion Laboratory, California Institute of Technology, 4800 Oak Grove Dr.,
  Pasadena, CA; Longuski, Professor, Purdue University, School of Aeronautics and Astronautics, 701 W.
  Stadium Ave., West Lafayette, IN, AAS Member, AIAA Associate Fellow.
- Venue, place, date, volume: **not stated.** The paper prints no conference name, city, month or year of
  publication. The only dated items are mission dates (launch 2025) and reference dates. (The 2023 Landau
  et al. paper cites it as "AAS/AIAA Spaceflight Mechanics Meeting, American Astronautical Society, 2014,
  pp. 1469-1485. AAS 13-801"; this paper does not say so itself.) The page numbers 1469-1485 match.
- Acknowledgment printed: "This research was carried out at Purdue University and the Jet Propulsion
  Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space
  Administration."
- Filed in the private paper corpus as
  `strange-landau-longuski-2013-initial-inclination-reduction-sequence-uranian-gravity-assist-tours-AAS-13-801.pdf`,
  md5 `7a5a17d68b42f87c57f05eddc7f57a7f`.

## What the paper is
Abstract (quoted from the rendered page): "Although a gravity-assist tour of the Uranian moons would be
desirable component of a Uranus mission, such tours are especially challenging due to its distance from the
Sun and the planet's very high obliquity (97.77 deg). The high obliquity means that the initial orbits at
Uranus tend to be very highly inclined (60 deg - 80 deg), except in the rare case of arrival during the
Uranian equinox (the equinoxes occur every 42 years, with the next one in 2049). The long flight time to
Uranus means that there may be precious little time left in a mission for inclination reduction flybys to
reach the Equatorial plane. This paper presents a method for the design of the initial capture orbit
maneuvers to target a satellite v-infinity that allows for an efficient gravity assist inclination
reduction sequence. We also provide an example case for a 2025 mission with a 13-year trajectory to Uranus,
a 1-year inclination reduction sequence, and a 2.5 km/s total mission delta-V."

It is a method paper on the initial inclination-reduction phase only. Quoted scope: "In this paper, we will
focus on the design of initial inclination reduction sequence needed to set up an equatorial multiple
satellite Uranian tour. We will not address other aspects of potential Uranian mission designs such as the
delivery of atmospheric probes and the design of a science tour of the Uranian satellites." Sections as
printed: Introduction; Analysis (Parameterization of Maneuvers in Terms of Twist and Incline Angles;
Location of Inclination Change Maneuvers; Arrival Conditions; Changing the Orbit with Flybys; Pump Angle;
Crank Angle; The Tisserand Invariant; Flyby Bending Angle; Finding Inclination Reduction Sequences);
Example Mission (Tisserand Graph for Uranian Satellites; delta-V and Flight Time); Conclusion;
Acknowledgment; Notation; References. Eleven figures, no numbered tables. Equations (1) to (64) are
about two-thirds of the paper.

## The design method
- Model, quoted: "In the theory of zero sphere of influence patched conics, the gravity assist is treated as
  an instantaneous rotation of the v-infinity vector. A flyby cannot change the v-infinity magnitude."
  "By the patched-conic assumption, the spacecraft's velocity relative the central body is given by the
  vector sum of the gravity-assist body velocity and the spacecraft's v-infinity with respect to the
  gravity-assist body." Ephemeris for the moons: not stated. The gravity-assist body orbit is taken as
  circular in the Tisserand relation (Eq. 49: "we will assume the gravity-assist body to be in a circular
  orbit").
- Geometry: orbit insertion and apoapsis manoeuvres parameterised by twist angle (xi) and incline angle
  (zeta); flybys decomposed into pump angle (alpha) and crank angle (kappa) changes using the bending angle
  delta (Eq. 58, sin(delta/2) = mu_ga / (mu_ga + r_p,fb v_inf^2)); Tisserand invariant relation
  v_inf^2 / v_c^2 = 3 - C_Tiss (Eq. 49).
- Procedure, quoted: "We begin our analysis in the middle at the first flyby of the inclination reduction
  sequence and work our way out. ... For a given target moon and v-infinity, we can optimize the free
  parameters subject to the node crossing constraints to find the lowest delta-V solution. For the example
  mission in the following section this was done using MATLAB's fmincon function."
- Constant V-infinity assumption, quoted: "The inclination reduction sequence cannot change the
  v-infinity with respect to the flyby moon without the addition of v-infinity leveraging maneuvers or
  backflip transfers with other moons. In this analysis, we will neglect those cases and assume the
  v-infinity is constant. This means that as inclination is reduced, the eccentricity will increase (by
  Eq. 49) and the periapsis will decrease. Placing periapsis at the outer edge of the ring system then
  provides a lower bound on orbit period. We will target the minimum integer resonance (i.e. 1:1, 2:1,
  3:1, etc.) that avoids the rings. Until that period is achieved, the flybys will pump down to the lowest
  integer resonance achievable and use any excess bending to reduce crank. After we reach that period,
  flybys will reduce crank until inclination is zero."
- Transfers between moons: not designed. Only flybys of one target moon (Titania or Oberon) are used for
  the inclination reduction in the example; the paper says it does not treat v-infinity leveraging or
  "backflip transfers with other moons".
- Manoeuvres: orbit insertion plus periapsis-raise at Uranus ("Satellite tours such as those flown by
  Galileo and Cassini typically start with an orbit insertion maneuver and a periapsis raise maneuver that
  will together target a v-infinity at a moon. We will use this same approach at Uranus.") and a
  deep-space manoeuvre in cruise. The flybys themselves are treated as ballistic bends; no per-flyby
  manoeuvre is printed.
- Example constraints, quoted: "We will limit the Uranian periapsis altitude to no lower than 5,000 km and
  the ring plane crossing to no lower than 51,140 km to avoid the rings. In this example, our first
  Uranian moon flyby will be 200 km and subsequent flybys will limited to 50 km or higher." The Figure 9
  caption text: "The v-infinity contours have tick marks denoting the spacing between 50 km flybys of each
  moon."

## The tour and flyby lists
There is **no moon-flyby sequence table or list** in this paper: no flyby numbers, moon-by-moon dates,
altitudes, or per-flyby V-infinity for the inclination-reduction sequence. The only dated list is the
interplanetary cruise of the example, in Figure 8 ("Example 2025 EESU trajectory"; the body text calls it
"a 13 year delta-V-EGA with a Saturn flyby, launching in 2025"). Transcribed from the rendered page at
250 dpi (the date of event 3 is overprinted in the print; it reads 07/23/2028):

| Event | Printed label |
|---|---|
| 1: Earth | 06/02/2025, C3 = 51.2 km^2/s^2, Dec. = -7.8 deg |
| 2: DSM | 10/26/2026, delta-V = 616 m/s |
| 3: Earth | 07/23/2028, V-inf = 12.558 km/s, Alt = 500 km |
| 4: Saturn | 11/12/2031, V-inf = 8.842 km/s, Alt = 1250000 km |
| 5: Uranus | 05/30/2038, V-inf = 6.524 km/s, Dec. = 75.6 deg |

Text for the same trajectory: "It arrives with a v-infinity of 6.5 km/s and a declination of 75.6 deg."
Date format taken as mm/dd/yyyy (05/30/2038 and 10/26/2026 parse only that way).

Other numbers the paper prints for the inclination-reduction sequence (text and figures):
- "Our goal will be to find a 1-year reduction sequence, which would then leave us enough time for a 2-year
  planar tour in a 15 year mission."
- "We see that a Titania 3:1 resonance (26 days) is very close to an Oberon 2:1 resonance (27 days). To
  avoid the rings at this resonance we need a v-infinity of 4.5 km/s or less at Titania or a v-infinity of
  less than about 3.9 km/s at Oberon. We also note that a v-infinity of 2 km/s is the lowest v-infinity at
  either body to allow periapses near the rings."
- "We see that to achieve a 1 year inclination reduction sequence we will need a v-infinity of about 2.2
  km/s at Titania or about 2 km/s at Oberon. This corresponds to an initial orbit dv of 1.86 km/s for
  Titania and 1.9 km/s for Oberon."
- Conclusion: "the initial inclination can be reduced from 76 deg to 0 in 1 year for a total mission
  deterministic delta-V of 2.5 km/s (including orbit insertion and deep space maneuvers)."
- Figure 9, a Tisserand (r_a, r_p) graph for Oberon (black) and Titania (blue) at zero inclination, labelled
  in the image with v-infinity contours 0.5 to about 4.5 km/s (Titania) and 0.5 to 4.0 km/s (Oberon) and
  integer resonance lines (1:1 up to 10:1 for Oberon and 1:1 up to 15 for Titania as legible on the
  page). Figure 10 plots "Insertion + Raise delta-V" (about 1.27 to 1.9 km/s) against satellite target
  v-infinity (2 to 6.5 km/s); Figure 11 plots flight time of the inclination-reduction tour (0.5 to 6
  years) against target v-infinity (2 to 5 km/s). Values read from these figures are approximate (by eye):
  at 2 km/s, Figure 10 gives about 1.89 km/s for Oberon and about 1.86 km/s for Titania, agreeing with the
  text; Figure 11 shows about 1.2 years for Oberon at 2 km/s and about 0.9 years for Titania at 2.2 km/s,
  not exactly the "1 year" of the text.

## Answers to the specific questions
**(4) Alternation of two different moons.** None. The paper designs flybys of one moon at a time (Titania or
Oberon, as alternative choices, shown as two curves in Figures 9 to 11), not a sequence in which the two
alternate. It remarks that Titania 3:1 (26 days) is close to Oberon 2:1 (27 days) in period, and that
v-infinity leveraging and "backflip transfers with other moons" are neglected.

**Repeated flybys of the same moon and the resonance stated.** The method "will target the minimum integer
resonance (i.e. 1:1, 2:1, 3:1, etc.) that avoids the rings", with "flybys" (plural) of the one target moon
pumping the orbit down to that resonance and then reducing crank. No specific resonance is stated for the
example sequence, and no number of flybys or intervals between flybys is printed. The only specific
resonances named are the Titania 3:1 (26 days) and Oberon 2:1 (27 days) in the Tisserand-graph passage.

**(5) Word search** (full text layer, lines joined, ligatures normalised, case-insensitive; figure-image
labels such as "1:1 Oberon Resonance" are not in the count): "cycler" 0; "cycling" 0; "periodic" 0;
"repeating" 0; "repeat" 0; "free-return" 0; "resonant" 0; "petal" 0; "resonance" 6 (the passages
quoted above: "minimum integer resonance", "lowest integer resonance achievable", "the integer resonances
with each satellite", "a Titania 3:1 resonance", "an Oberon 2:1 resonance", "at this resonance"). No
trajectory is described as returning to the same two moons repeatedly or indefinitely.

**(6) V-infinity, duration, delta-v.** No single lowest or highest moon-flyby V-infinity is printed for a
designed tour. The lowest value the paper calls feasible is 2 km/s ("the lowest v-infinity at either body
to allow periapses near the rings"); Figure 9 shows contour labels down to 0.5 km/s, and Figure 10 starts
at 2 km/s. The highest values: the introduction says initial high-inclination capture orbits begin "with
high v-infinity relative to the satellites (5-8 km/s)"; Figure 10 extends to 6.5 km/s; the Uranus arrival
V-infinity of the example is 6.5 km/s (6.524). The example targets about 2 to 2.2 km/s at the moon for a
1-year sequence. Duration and delta-v of the example: 13-year cruise, 1-year inclination reduction
(76 deg to 0), 15-year mission with a 2-year planar tour; total deterministic delta-V 2.5 km/s including
orbit insertion and deep-space manoeuvres (the DSM is 616 m/s; the initial orbit dv is 1.86 km/s for
Titania or 1.9 km/s for Oberon).

## What it does NOT contain
- No flyby-by-flyby table: no flyby count, dates, altitudes (beyond the 200 km first and 50 km later
  limits), V-infinity, resonance ratio or delta-v for the inclination-reduction sequence itself.
- No two-moon alternation, no cycler, periodic, repeating or free-return trajectory; the words do not occur.
- No design of the equatorial science tour that follows ("We will not address ... the design of a science
  tour of the Uranian satellites"), no Miranda, Ariel or Umbriel flybys, and no moon-to-moon transfers
  (v-infinity leveraging and backflip transfers are explicitly neglected).
- No moon ephemeris, no n-body or numerical integration; no verification of the example in a higher-fidelity
  model.
- No statement of venue or date for the paper itself.
- The inclination-reduction sequence is not itemised; the example is described only through its
  total delta-V, duration and the v-infinity values above.

## Relevant bibliography (as printed in this paper)
- [1] Vision and Voyages for Planetary Science in the Decade 2013-2022. National Academies Press, 2011.
- [2] D. Landau, T. Lam, and N. Strange, "Broad Search and Optimization of Solar Electric Propulsion
  Trajectories to Uranus and Neptune," AAS/AIAA Astrodynamics Specialists' Conference, Aug. 2009. AAS Paper
  09-428.
- [3] A. Heaton and J. Longuski, "Feasibility of a Galileo-Style Tour of the Uranian Satellites," Journal of
  Spacecraft and Rockets, Vol. 40, No. 4, 2003, pp. 591-596.
- [4] J. McAdams, C. Scott, Y. Guo, J. Dankanich, and R. Russell, "Conceptual Mission Design of a Polar
  Uranus Orbiter and Satellite Tour," Space Flight Mechanics Conference, Feb. 2011. AAS Paper 11-188.
  (Cited for the Decadal Survey study that "was not able to achieve any significant orbit shaping with
  gravity-assists of the Uranian Satellites".)
- [5] J. Deerwester, J. McLaughlin, and J. Wolfe, "Earth-Departure Plane Change and Launch Window
  Considerations for Interplanetary Missions," Journal of Spacecraft, Vol. 3, No. 2, 1966, pp. 169-174.
- [6] R. Luidens and B. Miller, "Efficient Planetary Parking Orbits with Examples for Mars," NASA TN D-3220,
  Jan. 1966.
- [7] D. Landau, J. Longuski, and P. Penzo, "Method for Parking-Orbit Reorientation for Human Missions to
  Mars," Journal of Spacecraft and Rockets, Vol. 42, No. 3, 2005, pp. 517-522.
- [9] C. Uphoff, P. H. Roberts, and L. D. Friedman, "Orbit Design Concepts for Jupiter Orbiter Missions,"
  AIAA Mechanics and Control Conference, Aug. 1974. AIAA Paper 74-781.
- [13] N. Strange and J. Longuski, "Graphical Method for Gravity-Assist Trajectory Design," Journal of
  Spacecraft and Rockets, Vol. 39, No. 1, 2001, pp. 9-16.
- [14] S. Campagnola and R. Russell, "The Endgame Problem Part 2: The Multi-Body Technique and the T-P
  Graph," Journal of Guidance, Control and Dynamics, Vol. 33, No. 2, 2010, pp. 476-486.
- Listed, not about tour design: [8] Bowditch, The American Practical Navigator (2002); [10] Tisserand,
  Traite de Mecanique Celeste, Vol. 4 (1896); [11] Brouke, "On the History of the Slingshot Effect and
  Cometary Orbits," AAS Paper 01-435; [12] Battin, An Introduction to the Mathematics and Methods of
  Astrodynamics (1999).
