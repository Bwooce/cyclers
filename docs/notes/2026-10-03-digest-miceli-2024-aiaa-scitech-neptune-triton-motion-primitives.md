# Digest: Miceli, Bosanac, Stuart & Alibay (2024), "Motion Primitive Approach to Spacecraft Trajectory Design in the Neptune-Triton System"

**Citation as printed on the PDF.** G. E. Miceli (University of Colorado Boulder), N. Bosanac
(University of Colorado Boulder), J. R. Stuart and F. Alibay (Jet Propulsion Laboratory),
"Motion Primitive Approach to Spacecraft Trajectory Design in the Neptune-Triton System".
The text layer of this author version prints NO venue, NO paper number and NO date: there is no
"AIAA SciTech" string, no "2024" and no DOI anywhere in the 20 pages (text search finds none).
The venue, the paper number (AIAA 2024-1280) and the DOI (10.2514/6.2024-1280) come from the
task brief and CrossRef, not from the PDF. The 2026 journal digest records this paper as its
reference [62], "the AIAA-2024 precursor". The only dates printed in the PDF are in the
reference list and body: reference [50] is a private communication of July 2023, references [22],
[23], [49] carry retrieval dates of November 2023, and the initial-state epoch is
October 2, 2045 (p. 10). The PDF file metadata (not printed text) gives a creation date of
2023-11-14.

Filed in the private paper corpus as
`miceli-bosanac-stuart-alibay-2024-motion-primitive-approach-trajectory-design-neptune-triton-aiaa-scitech-2024-1280-author-version.pdf`
(md5 `f0652877af6a3db44d8e2b16035a3ff8`, 20 pages, text layer; a full read of every page
including figures 4 and 5 rendered as images).

Context: companion to `docs/notes/2026-08-01-776-miceli-bosanac-2026-neptune-triton-digest.md`
(the 2026 JAS paper, Miceli and Bosanac). Purpose of this digest: establish exactly which
periodic orbits and families this conference paper computes, so that the project's two
symmetric periodic orbits near the 4:5 saddle homoclinics can be checked against it.

## 1. Model and constants (all as printed)

* Planar-and-spatial formulation written for the Neptune-Triton CR3BP, rotating frame with the
  system barycentre at the origin, x-hat from Neptune to Triton (p. 2). The transfer design
  itself is planar (p. 10, p. 11).
* Characteristic quantities (p. 3): "m* ~= 1.024569 x 10^26 kg is the sum of the masses of
  Neptune and Triton; l* = 354, 760 km is set equal to the average distance between Neptune and
  Triton, and t* ~= 8.081353 x 10^4 s sets the mean motion of the primary system to unity".
* Mass ratio (p. 3): "mu = M2/(M1 + M2) ~= 0.00020895" (rounded display value; no more digits
  printed). Primaries' average eccentricity quoted as 0.000016 (p. 2, from the NASA fact sheet).
* Potential (Eq. 2, p. 3): U* = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2, with r1 = sqrt((x + mu)^2 + ...)
  and r2 = sqrt((x - 1 + mu)^2 + ...) (Eq. 3). Neptune sits at x = -mu, Triton at x = 1 - mu.
* Jacobi constant (p. 3): "C_J = 2U* - xdot^2 - ydot^2 - zdot^2". This is the 2U* - v^2 convention.
* Stability indices (p. 4): s1, s2 are the sum of each non-trivial eigenvalue pair of the
  monodromy matrix; a pair with index in [2, -2] means bounded motion nearby, otherwise the
  orbit admits stable and unstable manifolds.
* Differential correction tolerance for periodic orbits 1e-10 (multiple shooting, p. 3-4);
  collocation constraint tolerance 1e-12 (p. 5).

## 2. Periodic orbits and families the paper computes or uses

Overall: the paper prints NO initial condition and NO period for any periodic orbit, and prints
a Jacobi constant only for the items listed below. Families are shown in figures (Figs. 3-5)
and summarised by Jacobi-constant ranges only graphically (Fig. 4).

Family list as printed (p. 12, Sec. IV.C), "planar, periodic orbit families and hyperbolic
invariant manifolds":

* L3 Lyapunov family
* 1:2 prograde resonant family, periapsis on the -x-axis
* 1:3 prograde resonant family, periapsis on the -x-axis
* 1:4 prograde resonant family, periapsis on the -x-axis
* 1:5 resonant family, periapsis on the -x-axis
* 3:1 retrograde resonant family, periapsis on the +x-axis
* 3:1 retrograde resonant family, periapsis on the -x-axis
* 4:1 retrograde resonant family
* 2:3 prograde resonant family
* 3:4 prograde resonant family, periapsis on the -x-axis
* 3:5 prograde resonant family
* 4:5 resonant family, periapsis on the -x-axis
* manifolds (listed separately, see Section 3).

The Fig. 5 graph labels (p. 15) additionally show a "1:1 Res S" family tag that is not in the
printed text list (figure only). "S stands for symmetric, R for retrograde and M for
manifolds" (Fig. 5 caption). Resonance convention (p. 3): a p:q orbit completes p revolutions
about Neptune in about the time Triton completes q revolutions; "When p > q, the orbit is
labeled an interior resonance ... when p < q ... an exterior resonance". Initial guesses come
from a two-body resonant ellipse with the eccentricity varied (p. 3), then multiple shooting and
pseudo-arclength continuation (p. 4).

Per-orbit data actually printed:

* Target orbit (p. 11, Sec. IV.B): "a prograde, planar, 3:4 resonant orbit at C_J,f = 1.75598".
  Chosen "because it satisfies the requirements with the lowest value of the Jacobi constant";
  the requirements (from Cochrane et al.) are a 300 - 350 km altitude pass of Triton with
  closest passages every 3/2 period of Triton's orbit (p. 11). Intro (p. 2) says it "performs two
  flybys with Triton at an altitude between 300 and 350 km". Period of the target: not stated.
  Initial condition: not stated. Shown in Fig. 2 (p. 12).
* Initial state (p. 10-11): periapsis epoch 2045-10-02 11:52:51 UTC, periapsis altitude 2460.11 km
  above Neptune's surface, v-infinity 11.5252 km/s, declination 8.3778 deg; spatial Jacobi
  constant 0.950382, inclination 20.64 deg; the planar modification has C_J,0 = 0.963141.
* Table 1 (p. 13), periodic orbits whose manifolds are generated, transcribed exactly:

  | Orbit family | Jacobi constant of selected periodic orbits |
  |---|---|
  | 1:2 prograde resonant orbit family | C_J = [1.50, 1.56, 1.61, 1.65, 1.70, 1.75, 1.81] |
  | 1:3 prograde resonant orbit family | C_J = [1.22, 1.30, 1.35, 1.40, 1.45, 1.51] |
  | 1:4 prograde resonant orbit family | C_J = [1.15, 1.17, 1.20, 1.25] |
  | 1:5 prograde resonant orbit family | C_J = [1.09, 1.10, 1.12, 1.13, 1.15] |
  | 3:4 prograde resonant orbit family | C_J = 1.75598 |

  Caption: "Periodic orbits used to generate stable and unstable manifolds". Stability indices,
  periods and initial conditions are not tabulated. The text says only that "The resonant orbits
  used to generate stable and unstable manifolds also tend to have in-plane stability indices
  that are close to 2" (p. 12).
* Fig. 4 (p. 14) plots, per family, the C_J span of each family's primitives (log-style axis with
  ticks 0.96, 1.5, 1.755, 2, 2.5, 3). Read off the plot only, no numbers printed: the 4:5 Res
  row spans roughly C_J 1.6 to 2.7. These are plot readings, not text.
* Library window (p. 12-13): only primitives with periodic-orbit C_J between C_J,0 = 0.963141
  and C_J,f = 1.75598 are used in the graph.

### The 4:5 orbit, role and quotations

A 4:5 resonant family appears in exactly two roles, both as a library family, never as a
manifold-bearing family:

1. In the family list (p. 12): "4:5 resonant family with periapsis on the -x-axis" (quoted).
2. In Fig. 4c (p. 14, row "4:5 Res") and in Fig. 5 (p. 15, tag "4:5 Res" in the second
   high-level group of primitives).

Table 1 has no 4:5 row, and the manifold bullets (p. 12) cover only the 1:2, 1:3, 1:4, 1:5
and 3:4 families. A text search for "4:5" finds the single line in the family list. No 4:5
orbit is singled out by Jacobi constant, period or initial condition.

## 3. Invariant manifolds

Quoted list (p. 12): "Stable and unstable invariant manifolds associated with 7 members of the
1:2 prograde resonant orbit family", "... 6 members of the 1:3 ...", "... 4 members of the 1:4
...", "... 5 members of the 1:5 prograde resonant orbit family" and "Stable invariant manifolds
associated with the target 3:4 prograde resonant orbit". (Fig. 3d caption on p. 13 and the text
on p. 13 describe the 3:4 manifold as "stable"; Fig. 3 caption says "3:4 target resonant orbit
unstable manifolds" for panel d, an inconsistency inside the paper; Fig. 5 tags it "3:4 Res SM".)
Manifold method (p. 4): half-manifolds from a perturbation along the stable or unstable
eigenvector, propagated backward or forward. Segments are stored only "once they sufficiently
depart the vicinity of the resonant orbits" (p. 12), split at curvature extrema, sampled with 10
nodes and clustered. They are used as arcs (primitives) in the graph.

Search of the text for "homoclinic" and "heteroclinic": zero matches for each. No sentence in
the paper constructs or discusses a connection of any orbit to itself or to another orbit as a
stand-alone object. The word "connect" appears only for graph connectivity of primitives (for
example "A graph that captures their potential connectivity is then searched", p. 1) and
"intersecting one of the primaries" is a continuation termination (p. 4). Whether manifold
arcs happen to link in a graph path is a matter of the numerical q metric (Eq. 17), not a
computed manifold intersection.

## 4. Design scenario

* Start (p. 10-11): the planar Neptune orbit insertion state above, propagated 7.5 days either
  way; "the first maneuver [may] occur at any state within 3.7 days before or after the initial
  state" (p. 11).
* End: the 3:4 resonant target orbit at C_J,f = 1.75598.
* Library (p. 12-14): the families above, clustered by WEAC (n_clust = [3, 11], ward linkage,
  p. 12), curvature-extrema sampling, medoid primitives with exemplars.
* Graph (p. 13-15, Fig. 5): four high-level components (initial arc, lower-C_J primitives,
  primitives near C_J,f, target). alpha_pos = 5, alpha_vel = 1 (Eq. 17).
* Search (p. 15): Dijkstra and Yen (top 1,000, path length up to 50 primitives) plus depth-first
  search (sequences of 5 primitives). Yen paths all had at least 10 primitives (p. 15).
* Correction (p. 15-16): collocation with IPOPT, objective J = w_geo (dr)^2 + w_man sum(dv_i)^2
  (Eq. 18), continuation from [0.9, 0.1] to [0.1, 0.9] in typically 10 steps. Maneuvers at
  apses, curvature extrema and segment joins.
* Results (p. 16-17): of the five initial guesses (Fig. 6), c) hits Neptune and a) gives an
  "excessively high total delta-v"; Fig. 7 shows corrected b), d), e). After continuation:
  TOF 36.21 days with 4.05 km/s; TOF 12.60 days with 1.64 km/s; TOF 12.93 days with 1.87 km/s
  (p. 17). The lowest, 1.64 km/s, is "61% of the available delta-v = 2.708 km/s" of the Uranus
  Orbiter and Probe study [51] (p. 17); Uranus orbit insertion there is 1.0867 km/s. Per-trajectory
  numbers in Fig. 7 and Fig. 9 are in figure images and were not transcribed.

## 5. Reproduction targets

The paper tabulates almost nothing numerical that is usable as a reproduction target. Usable
values: mu = 0.00020895, l*, t*, m*, the target C_J,f = 1.75598, C_J,0 = 0.963141, Table 1
Jacobi constants (a Jacobi-constant list per family, to two decimals except the 3:4 orbit), and
the delta-v and flight-time figures above. There is no initial-condition table, no period table
and no eigenvalue table. Compare the 2026 paper, which ships machine-readable supplementary
data (see that digest). This conference paper gives nothing at the 1e-6 level except the 3:4
C_J and the two C_J values for the initial state.

## 6. Differences from the 2026 JAS paper (as far as the two digests allow)

* Target: here a 3:4 exterior resonant orbit, C_J = 1.75598 (a capture-then-science orbit,
  from a high-energy arrival). The 2026 digest describes Scenario 1 (1:7 resonant target,
  C_J = 1.8, period 41.14 days) and Scenario 2 (3:2 resonant orbit to a Triton-region low
  prograde orbit). Neither of those targets is in this paper.
* Scope: one planar scenario here; two scenarios in 2026.
* Families: this paper lists the 12 families of Section 2. The 2026 digest does not enumerate
  its families, so a family-by-family comparison cannot be made from the digests; not stated.
* Data: no supplementary files here; the 2026 paper ships four.
* Method: the same Smith and Bosanac motion-primitive method; this paper's improvements are the
  curvature-extrema sampling and the Yen search (p. 6, p. 9).

## 7. What the paper does NOT contain (text searches on the extracted text)

* "homoclinic": 0 matches. "heteroclinic": 0 matches.
* "cycler": 0 matches. No repeated-flyby periodic trajectory is designed or named as such. The
  only repeated Triton passage is the 3:4 target orbit itself, with passes every 3/2 Triton
  period (p. 11).
* Any periodic orbit near a homoclinic connection: none. No orbit is described as lying near a
  manifold intersection.
* Any orbit with period near 68.75 or 86.90 nondimensional units: a search for "68." finds
  only the URL of reference [46]; "86." finds nothing. No period is printed for any orbit at
  all, so a period comparison is impossible from the text. The only way the 4:5 family could
  hold such orbits is through its family curve (Fig. 4c), not through any value printed.
* Any initial condition x = 1.16933872 or x = -1.38561105: no initial condition of any kind is
  printed.
* A 4:5 orbit's Jacobi constant, period, stability or manifolds: not stated (Section 2).

Checking caveat: Figs. 1-9 are raster/vector images; their internal labels were inspected only
for Figs. 4 and 5. Figures 3, 6, 7, 8, 9 were not decoded beyond their captions.

## 8. Bibliography entries about Neptune and Triton trajectory design (transcribed as printed)

* [1] National Academies of Sciences, Engineering, and Medicine, "NASA 2023 Decadal Survey,"
  Technical report, National Academies Press, Washington, D.C., 2023.
* [2] Swenson, B., "Neptune atmospheric probe mission," Guidance, Navigation and Control
  Conference, American Institute of Aeronautics and Astronautics, 1992.
  https://doi.org/10.2514/6.1992-4371.
* [3] Masters, A., Achilleos, N., Agnor, C., Campagnola, S., Charnoz, S., Christophe, B.,
  Coates, A., Fletcher, L., Jones, G., Lamy, L., Marzari, F., Nettelmann, N., Ruiz, J., Ambrosi,
  R., Andre, N., Bhardwaj, A., Fortney, J., Hansen, C., Helled, R., Moragas-Klostermeyer, G.,
  Orton, G., Ray, L., Reynaud, S., Sergis, N., Srama, R., and Volwerk, M., "Neptune and Triton:
  Essential pieces of the Solar System puzzle," Planetary and Space Science, Vol. 104, 2014,
  pp. 108-121. https://doi.org/10.1016/j.pss.2014.05.008.
  (Campagnola appears only as the fourth author of this entry; there is no separate Campagnola
  trajectory-design reference.)
* [4] Rymer, A. M. et al., "Neptune Odyssey: A Flagship Concept for the Exploration of the
  Neptune-Triton System," Planet. Sci. J., Vol. 2, No. 5, 2021, p. 184.
  https://doi.org/10.3847/PSJ/abf654. (Long author list abbreviated here; the PDF prints it in
  full.)
* [6] Marley, M. e. a., "Planetary Science Decadal Survey JPL Rapid Mission Architecture
  Neptune-Triton-KBO Study Final Report," Technical report, NASA Jet Propulsion Laboratory,
  Pasadena, CA, 2010.
* [7] Melman, J., Orlando, G., Safipour, E., Mooij, E., and Noomen, R., "Trajectory Optimization
  for a Mission to Neptune and Triton," AIAA/AAS Astrodynamics Specialist Conference and
  Exhibit, 2007, p. 7366.
* [20] Cochrane, C. J., Persinger, R. R., Vance, S. D., Midkiff, E. L., Castillo-Rogez, J.,
  Luspay-Kuti, A., Liuzzo, L., Paty, C., Mitchell, K. L., and Prockter, L. M., "Single- and
  Multi-Pass Magnetometric Subsurface Ocean Detection and Characterization in Icy Worlds Using
  Principal Component Analysis (PCA): Application to Triton," Earth and Space Science, Vol. 9,
  No. 2, 2022. https://doi.org/10.1029/2021EA002034. (Source of the 300 - 350 km, 3/2-period
  magnetometry requirement.)
* [51] Simon, A., Nimmo, F., and Anderson, R. C., "Journey to an Ice Giant System: Uranus
  Orbiter and Probe," Planetary Mission Concept for the 2023-2032 Planetary Science Decadal
  Survey, June 7 2021. Retrieved May 1, 2022. (Uranus, used for the delta-v comparison.)
* [50] "Private communication with Dr. Reza Karimi," July 2023, NASA Jet Propulsion Laboratory
  (source of the arrival state).

Method-lineage references used (not Neptune-specific): [9] Smith and Bosanac, CMDA 134(1), 2022,
p. 7; [10] Smith and Bosanac, AAS/AIAA Astrodynamics Specialists Conference, 2022; [25] Vaquero
Escribano, Purdue PhD thesis, 2013 (resonant-orbit guess construction, p. 3); [31] Smith and
Bosanac, JAS 70(34), 2023. Mention in text: "Marley et al. [6] and Masters [3]" for patched-conic
multiple-Triton-flyby studies, and "Melman et al. computed a transfer ... from a Neptune-centered
orbit to a circular near-polar orbit around Triton in the CR3BP [7]" (p. 2). No Campagnola
trajectory paper is in the list.

## Verdict (for the corpus index)

One-line description: Miceli, Bosanac, Stuart and Alibay (AIAA SciTech 2024 author version),
planar Neptune-Triton CR3BP motion-primitive transfer from a 2045 interplanetary arrival
(C_J 0.963141) to a 3:4 resonant target (C_J 1.75598), 12.60 days and 1.64 km/s; 12 resonant and
L3-Lyapunov library families (including a 4:5 family) and manifolds of the 1:2, 1:3, 1:4, 1:5 and
3:4 orbits; no initial conditions, periods, homoclinics or cyclers printed.

Does not contain the project's two near-homoclinic 4:5 orbits (x = 1.16933872, x = -1.38561105;
periods 68.7485 and 86.898): no initial conditions or periods are printed for any orbit, and the
words homoclinic and heteroclinic do not occur.
