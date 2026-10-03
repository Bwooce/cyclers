# Digest: Owen & Baresi 2023 (AAS 23-110), conference version of the knot-theory heteroclinic paper

**Paper:** "Detecting Heteroclinic Connections Between Quasi-Periodic Invariant Tori Using Knot Theory"
**Authors:** Danny Owen (PhD student), Nicola Baresi (Lecturer), Surrey Space Centre, University of Surrey
**Number printed on the paper:** AAS 23-110 (page 1 header).
**Venue and date:** the paper body does NOT print a meeting name, place or month. The only dating in the file is
the repository cover sheet: "Owen, D., & Baresi, N. (2023)" and "@2023 The Authors", with "Downloaded On 2026/10/03".
The identification as the 33rd AAS/AIAA Space Flight Mechanics Meeting, Austin, January 2023, comes from the task
brief, not from the paper's text.
**Filed:** in the private paper corpus as
`owen-baresi-2023-detecting-heteroclinic-connections-quasi-periodic-tori-knot-theory-AAS-23-110-surrey-open-research.pdf`
(md5 `097c8bd686509364213fa49cbc8ef638`; 19 pages: a Surrey Open Research cover sheet, then 18 numbered pages).
Page references below use the paper's printed page numbers; the PDF page is the printed page plus one.
**Length and structure:** 17 pages of text and figures plus references. Sections: Introduction, Background
(linking number, CR3BP, quasi-periodic tori, heteroclinic connections), Methodology (manifold generation, level
curves, linking number, tracking, theta interpolation, differential correction), Results (one short section),
Discussion, Conclusion, References (24).

Companion record of the journal version: [[2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic]].
Context for why this matters: [[2026-07-10-postmortem-548-linking-number-pipeline]] (Addendum 2026-10-03, #883).

All figure readings below are by eye from rendered pages and are marked "(by eye)". The paper prints no table of
numeric results.

## 1. Every number and convention printed for the Earth-Moon examples

### The energy sentence

Results section, p. 14, verbatim, with the sentence before and after:

> "For each of the following applications, the perturbation value for the initialisation of the manifold was set to 1 x 10^-6.
> Connections were successfully found for combinations of quasi-halo to quasi-halo (Fig 12), quasi-halo to Lissajous (Fig 13), and Lissajous to Lissajous (Fig 14). **The energy level of these orbits is given by 3.1460717.**
> Heteroclinic connections that arrive or depart from a Lissajous orbit were produced successfully in this work."

Which orbits: "these orbits" follows directly the list of the three pairings, so it covers all the orbits in Figs 12, 13
and 14: the quasi-halo and Lissajous tori of all three Earth-Moon examples. One value, 3.1460717, is given for all of
them. The paper does not name the quantity as "Jacobi constant" in that sentence ("energy level"), but the Background
section (p. 5) defines the Jacobi integral C (Eq 5) and states that connections exist "only between orbits that are of
the same Jacobi integral, also known as isoenergetic orbits", and the Conclusion calls the tori "isoenergetic". The
digit count (7 decimals) is the only place any energy is printed in the whole paper; the Introduction, Background and
captions print no other value of C.

### Everything else printed

| Item | What the paper prints | Page |
|---|---|---|
| System | "In this work we will be considering the Earth-Moon system" | 4 |
| Mass ratio | defined symbolically only: mu = m/(M+m), Eq 1. No numeric value is printed | 4 |
| Units | distance = Earth-Moon separation, time = reciprocal of the Moon's angular velocity, both 1 | 4 |
| Frame | barycentric rotating frame; Earth at [-mu,0,0], Moon at [1-mu,0,0]; y parallel to the Moon's velocity | 4 |
| Equations of motion | Eqs 2-4 (standard CR3BP, q1,q2,q3); Jacobi integral C = Omega - (qdot1^2+qdot2^2+qdot3^2)/2, Eq 5 (Omega not defined in text) | 4-5 |
| Surface of section | x = 1 - mu, "the x position of the Moon"; the section the manifolds are propagated to | 7 |
| Manifold offset | epsilon = 1 x 10^-6, "for each of the following applications" | 14 |
| Offset form | manifold function U~(theta) = U(theta) + epsilon dU(theta); eigenvectors "appropriately scaled" (normalisation not stated) | 5, 7 |
| Torus method | GMOS (ref 23, Olikara and Scheeres 2012) "were implemented" | 5 |
| Torus angles | theta0 "longitudinal", theta1 "latitudinal"; both plotted over 0 to 2 pi | 5, 6, 13 |
| Grid | "N 6x1 state vectors" on an invariant circle (equally spaced latitude), M invariant circles at equally spaced longitude. N and M are NOT given numerically | 7 |
| Number of manifold trajectories | N x M per manifold; the value is not stated. Figs 9 and 12-14 show only a handful of corrected trajectories | 7 |
| Quasi-halo / Lissajous | the pairings quasi-halo to quasi-halo, quasi-halo to Lissajous, Lissajous to Lissajous are demonstrated | 14-16 |
| L1 / L2, north / south | the tori sit on either side of the Moon in the figures (see section 3); the paper never writes "L1", "L2", "northern" or "southern" about its own examples | 15-16 |
| Frequencies, rotation numbers | not stated | - |
| Periods, amplitudes | not stated (see by-eye extents in section 3) | - |
| Initial conditions | not stated | - |
| Fourier mode counts | not stated | - |
| Integration tolerances | not stated | - |
| Correction solver | MATLAB "fsolve" (named as an example: "such as") | 14 |
| Unknowns of the correction | the four torus angles (theta_s, theta_u); propagation times are NOT unknowns | 14 |
| Numbers of connections | "indicating 4 heteroclinic connections" (Fig 8 caption, p. 11); no counts printed for Lissajous cases | 11 |

## 2. The algorithm as the conference paper states it

1. **Tori.** GMOS torus function U(theta) = X, with the Floquet (monodromy) matrix of each invariant circle as a
   by-product (p. 5). "N 6x1 state vectors that lie on an invariant circle of the torus, each with equal longitude and
   with equally spaced latitude angles. We then produce M additional invariant circle sets at equally spaced longitude
   angles" (p. 7). Each circle's monodromy matrix gives dU_s(theta) and dU_u(theta).
2. **Manifold initialisation.** "The eigenvectors extracted are 6N x 1 vectors, which are appropriately scaled and
   applied to U(theta) to produce states" (p. 7), giving two N x M sets of states (stable, unstable). Torus maps built
   from the state variables let any point of the torus be initialised by interpolating the angles.
3. **Propagation.** "the stable manifold is propagated in backward time and the unstable manifold is propagated in
   forward time. They are propagated until they intersect the surface of section x = 1 - mu" (p. 7). The result sets
   are W_s and W_u. Which crossing is taken (first or second) is chosen in advance (p. 16).
4. **Torus maps on the section.** "We can map the values of y, z, ydot, and zdot at the surface of section to the
   surface of the torus by considering the value of theta at which each manifold trajectory was initialised,
   providing a torus map for each variable" (p. 7). x is fixed by the section; xdot is dropped (recovered later from C).
5. **Scanning variable.** z. "Consider the two sets of z-values at the surface of section, Z_s and Z_u ... the set
   Z_s intersect Z_u. This set contains every value of z which a heteroclinic connection could have as it passes the
   surface of section" (pp. 7-8). Level curves of the z torus map at each scanned z give theta curves.
6. **Linking curves.** "The values of theta0 and theta1 that correspond to these curves can be used to interpolate
   the values of y, ydot, and zdot for the same positions on their equivalent torus maps ... plotting the interpolated
   values of y, ydot, and zdot provides one or more closed curves in a 3-dimensional phase space" (p. 8). One set of
   curves from the stable map, one from the unstable map. "we have accounted for x via the surface of section and z as
   the curves are unique to a z value. That means an intersection in these curves would constitute a heteroclinic
   connection, as 5 state variables would be found to be equal" (p. 8). Axes of the curve plots are y, ydot, zdot
   (Figs 5, 7).
7. **Linking number of sets.** The paper says "while technically an inaccurate definition of the linking number" it
   uses L = sum over all pairs K(a_i, b_j) of the stable-curve set a and unstable-curve set b (Eq 6, p. 8).
8. **Computing K for a pair (p. 9, Fig 6).** Each curve is a polyline of n points. One curve is the "mesh curve":
   Q = (1/n) sum of its points (the centroid); n triangles [p_j, p_{j+1}, Q]; normal R = [p_{j+1} - p_j] x [Q - p_{j+1}].
   For each segment of the other curve, find the point P where the segment's line meets the triangle's plane; reject if
   the distance from p_j to P exceeds the segment length or is negative along the segment; otherwise "we take three dot
   products ... If all three resulting dot products are positive, the line segment must intersect the interior of the
   triangle" (ref 24, Amanatides and Choi ray tracing). Direction is the sign of the dot product of the segment with
   R; "positive intersections are given a value of 1, and negative intersections -1. The sum total of all
   intersections is taken to be the linking number."
9. **Scan.** "By repeating this process over a dense range of z values given by the intersection of the stable and
   unstable manifolds at the surface of section, we can track the linking number as z ... evolve (Fig 7, 8)" (p. 10).
10. **Initial guess in z.** "Any time a change in the linking number occurs, an average of the z values before and
    after the change is calculated and taken to be an initial guess for the z value of the heteroclinic connection"
    (p. 10).
11. **Initial guess state.** "the stable and unstable curves are interpolated for this value and a simple KNN search
    is performed to find the points in the curves which are most close to the other set. The average of these points
    is then taken, providing the y, ydot, and zdot values. The Jacobi integral is then used to calculate xdot"
    (p. 10). The paper does not state the sign convention for xdot.
12. **Recovering torus angles (pp. 11-13, Figs 10-11).** Level curves on the four torus maps (y, z, ydot, zdot) at the
    guess-state values X_i. Intersection classes: ab = (y,z), bc = (z,ydot), ca = (ydot,y), d = zdot with any of the
    others. "we construct every possible triangle that has one vertex from each of these three sets of intersections.
    We then compare each triangle with each intersection point in d, calculating the distance from each triangle vertex
    to the d intersection. These three distances are summed. The triangle with the smallest minimum summed distance is
    taken ... the mean of the positions of the triangle's vertices is calculated, providing a value of theta." Done for
    both tori. The paper calls the method "limited" and says "a more robust method of interpolating these values is in
    development".
13. **Correction (p. 14).** Psi(theta) = Psi(theta0, theta1) = X maps torus angles to the section state (Eq 7). The
    system to solve: Psi_s(theta_s) - Psi_u(theta_u) = X_s - X_u = 0 (Eq 8), "a boundary value problem ... solved using
    freely available nonlinear equation solvers from existing MATLAB libraries, such as 'fsolve'". "Unlike Henry and
    Scheeres, who included the propagation time to the surface of section, we are choosing to end the propagation at
    the surface of section, removing the need for these two additional variables". So 4 unknowns (theta0,theta1 for
    each torus) against the section-state difference. Which components of X enter the residual (4 or 5 or 6) is not
    stated; X_s and X_u are the section states "Psi" outputs.

### Stated here but not recorded in the journal digest

- The exact form of the intersection-triangle recipe in step 12 and its four classes (the journal digest only says
  "interpolate the torus angles").
- The correction unknowns are the four angles only and the propagation to the section is not a variable (the journal
  digest records the initial-guess steps only up to the 6D guess).
- The statement that the correction uses "fsolve" with Eqs 7-8.
- The statement that Fig 3 is a z map for a quasi-halo torus and the "Z_s intersect Z_u" definition of the scan
  range.
- The pseudo-code-level detail of the linking-number triangle test (Q = centroid, R definition, three dot products).
  The journal digest records this at the summary level only (Eq 19-20 there).
- The Discussion's statement of the miss condition: a connection is missed only if "the linking number changes twice
  within a very short range of z or two or more intersections occur near-simultaneously" (p. 14).

### Differences from the journal digest in the algorithm

- Scanning variable naming: conference uses z and curve axes (y, ydot, zdot); the journal digest uses a generic D and
  (A,B,C). Same content for the Earth-Moon example.
- The conference calls the extra work for Lissajous "adjusted to construct the torus maps using the second
  intersection" (section 4 below). The journal digest attributes the second-crossing requirement to the
  quasi-halo to Lissajous case only; the journal text itself (p. 13, 4.1.2) also applies it to Lissajous to Lissajous
  (checked against the journal PDF).
- Regularisation: the conference only mentions it in the Discussion (p. 16, citing Celletti, ref 22) and prints no
  equations. The journal has the regularisation section.

## 3. Results and figures

Connections: the Results section prints no counts and no table. Only the Fig 8 caption gives a count: "indicating 4
heteroclinic connections". Fig 9 is titled "Four initial guess heteroclinic connections propagated from the surface of
section". Figs 12, 13, 14 show corrected connections for quasi-halo to quasi-halo, quasi-halo to Lissajous and
Lissajous to Lissajous respectively. The paper never states which of the three pairings Figs 8, 9, 10, 11 belong
to; Fig 9 has the same geometry as Fig 12, and Fig 3 is stated to be a quasi-halo torus map.

| Figure (printed page) | What it shows | Readings (all by eye) |
|---|---|---|
| Fig 3 (p. 6) | z torus map on a quasi-halo torus, theta0 and theta1 each 0 to 6.28, colour bar x10^-3 | colour bar ticks -5 to 25 (x1e-3); map spans roughly -7e-3 to +27e-3; maximum near (theta0, theta1) = (3.1, 5.3); minimum near (5.8, 4.5) |
| Fig 4 (p. 8) | same map with z = 0 level curves in red | two red curves, near theta0 about 1.6 and about 4.7, nearly vertical with mild waviness |
| Fig 5 (p. 9) | linked stable (green) and unstable (red) curves in y-ydot-zdot | y axis 0 to about 0.08; ydot about -0.8 to 0.2; zdot about -0.5 to 1 |
| Fig 7 (p. 11) | curves at four z values: z = -0.0075837, -0.00088665, 0.0058104, 0.012507 | y axis 0 to 0.1 in all four; at z = -0.00088665 the ydot axis is about -2 to 2 and zdot about -4 to 0 (extreme near-Moon values); other panels ydot about -1 to 1, zdot about -1 to 1 |
| Fig 8 (p. 11) | linking number against z (x-axis in units of 1e-3) | linking number takes values 0, +1, 0, -1, 0 in order of increasing z, so four changes. Scanned z runs from about -8e-3 to about +24e-3 (axis ends). Changes at about z = -7.5e-3 (0 to +1), about -1.2e-3 (+1 to 0), about 0 (0 to -1), about +9e-3 (-1 to 0). The zero-valued gap between the second and third change is about 1e-3 wide (uncertainty of a few tenths of 1e-3) |
| Fig 9 (p. 12) | four initial-guess connections propagated from the section | x axis 0.8 to about 1.15; y -0.05 to 0.1; z -0.05 to 0.05; Moon marked as a black dot near x = 1 |
| Fig 10 (p. 13) | y, z, ydot, zdot level curves on the theta0-theta1 plane, 0 to 6.28 | initial-guess point near (theta0, theta1) = (2.2, 4.0) |
| Fig 11 (p. 13) | zoom of Fig 10 | theta0 2.17 to 2.26, theta1 3.945 to 4.04; guess near (2.23, 3.995); the four curves pass within about 0.01 rad of each other but do not meet at one point |
| Figs 12-14 (pp. 15-16) | corrected connections | axes x 0.85 to 1.1+, y -0.1 to 0.1, z -0.1 to 0.1; Moon dot near x = 1 |

Consistency of Fig 7 with Fig 8 (by eye): the four z values in Fig 7 fall at one in each of the intervals of Fig 8:
-0.0075837 near the first step up, -0.00088665 in the short zero gap, 0.0058104 in the -1 plateau, 0.012507 in the final
zero region.

Torus extents (by eye, Figs 9 and 12-14): the left torus lies at about x = 0.82 to 0.88 and the right torus at about
x = 1.10 to 1.17, one each side of the Moon dot near x = 1. Quasi-halo tori in Fig 12 appear as loops of roughly
0.05 to 0.1 in y and z; the Lissajous tori in Figs 13-14 appear flatter and denser. These figures give extent only to
the nearest 0.05 on axes labelled in steps of 0.1 and cannot be turned into amplitudes. The paper does not label the
left and right tori as L1 and L2.

Colour code (Fig 5 legend): green = stable curve, red = unstable curve. Figs 9 and 12-14 have no legend, but the
trajectories at the left torus are green and those at the right torus are red, which under the Fig 5 code would mean the
stable manifold attaches to the left torus and the unstable manifold to the right. The paper does not say so.

## 4. The Lissajous sentence

Results, p. 14, verbatim:

> "Heteroclinic connections that arrive or depart from a Lissajous orbit were produced successfully in this work.
> However, it has been shown that such connections do not exist for single passes of the Moon,6 so the process was
> adjusted to construct the torus maps using the second intersection with the surface of section instead of the first.
> This demonstrates a need to generalise the process to consider multiple revolutions of the Moon. This could be
> achieved by constructing the sets of stable and unstable curves based on the direction through which they pass the
> surface of section."

What was changed: for any pairing in which a Lissajous torus is the arrival or departure object (quasi-halo to
Lissajous, Lissajous to Lissajous), the torus maps (section 2, step 4) are built from the SECOND crossing of x = 1 - mu
by each manifold trajectory instead of the FIRST. Reference 6 is Gomez and Masdemont 2000 ("Some zero cost transfers
between libration point orbits"). The Discussion (p. 16) adds: "the number of crossings of the surface of section must
be designated before the process is implemented", with the proposed generalisation of storing multiple passes and
splitting the curve sets by crossing direction. The paper says nothing more about the second-crossing construction:
no definition of whether the second crossing is counted in one direction or both, and no treatment of trajectories
that never return.

The conference text does not say that quasi-halo to quasi-halo needs only one crossing. It only says the
adjustment applies to Lissajous connections, which by exclusion leaves the quasi-halo pair at the first crossing.

## 5. How the tori were computed and the energy question

- Method and reference: GMOS, ref 23, Z. P. Olikara and D. J. Scheeres, "Numerical method for computing
  quasi-periodic orbits and their stability in the restricted three-body problem", Advances in the Astronautical
  Sciences 145 (2012) 911-930 (p. 5, p. 18). "numerical methods ... developed by Gomez, Mondelo, Olikara, and Scheeres
  (GMOS) were implemented". No code, software, library or solver settings are named for the torus computation (the only
  named software is MATLAB's fsolve, p. 14). Olikara's thesis (ref 15) is cited elsewhere for his connection method,
  not as the source of the tori.
- Tori exist "in 2-parameter families, with variable energy and fundamental frequencies. Setting the Jacobi integral of
  the tori as constant, they then exist in 1-parameter families defined by their unique fundamental frequencies" (p. 5).
  The paper never says which frequency or amplitude was selected for each example.
- Shared energy: the energy sentence (section 1) gives one number for the quasi-halo and Lissajous tori alike, and the
  text requires isoenergetic pairs ("heteroclinic connections can only exist between orbits that are of the same Jacobi
  integral", p. 5). So within every pair the two tori share the printed energy, and the text assigns the same number to
  the quasi-halo and the Lissajous tori. The paper does not say how the equal energy was imposed (fixing C versus
  matching after the fact) and does not print the residual of the Jacobi constant.

## 6. Conference against journal

Journal values are from the journal PDF (Astrodynamics 8(4):577-595, 2024), pp. 13-14, and from the project's journal
digest where noted.

| Item | Conference AAS 23-110 | Journal 2024 |
|---|---|---|
| Earth-Moon Jacobi constant | "3.1460717" (p. 14), all three pairings | "3.15" ("All example connections in the Earth-Moon system presented in this paper have a Jacobi integral of 3.15", p. 13) |
| Mass ratio mu | symbolic only, no number | 0.012153643 (p. 13) |
| Tori used | quasi-halo to quasi-halo (Fig 12), quasi-halo to Lissajous (Fig 13), Lissajous to Lissajous (Fig 14) | same three pairings (4.1.1, 4.1.2, 4.1.3) |
| Frequencies | none printed | QH-QH latitudinal 0.2739 (L1) and 0.02163 (L2); L-L 0.3226 and 0.3578 (p. 13) |
| Frequencies for QH-L pair | none | not recorded in the journal's Earth-Moon section (journal digest does not list them) |
| Section | x = 1 - mu | x = 1 - mu |
| Offset epsilon | 1e-6 | 1e-6 (journal p. 5 and p. 13 on Jacobi drift) |
| Connections | 4 stated (Fig 8 caption); none stated for Lissajous cases | QH-QH 4; L-L 8 (p. 13); QH-L not counted in digest |
| Crossings | second intersection for Lissajous cases; QH-QH implied first | QH-QH first only; L-L second ("to demonstrate the method's robustness, torus maps are constructed using ... the second passing"); QH-L "often require an additional passing" |
| Systems | Earth-Moon only | Earth-Moon, Sun-Earth, Jupiter-Ganymede |
| Regularisation | one sentence, cites Celletti, no equations | full section (journal digest records KS-type regularisation Eqs 7-18) |
| Scanning variable | z (explicitly) | "D", chosen as z or zdot, with the symmetric-Lissajous reason given |
| Lissajous to Lissajous claim | demonstrated in Results (Fig 14) but Conclusion says "We are confident that it should prove successful" (internally inconsistent) | demonstrated, 8 connections |
| Linking-number scan range, Fig 8 against journal Fig 15 | scanned z about -8e-3 to +24e-3 (by eye); changes about -7.5e-3, -1.2e-3, 0, +9e-3 (by eye) | scanned z about -7e-3 to +20e-3 (axis -5 to 20; by eye); changes at about -7e-3, -1.5e-3, +1e-3, +9.5e-3 (by eye). Same 0, +1, 0, -1, 0 sequence |
| Linking-number scan range for Lissajous pair | not shown | journal Fig 17: z axis -0.06 to 0.06 (x1e-3), i.e. a range of order 1e-4, not 1e-2 |
| Correction | fsolve on 4 angles, Eq 8 | same ("fsolve", Eq 22) |

Correction to the project's journal digest (by eye): it records the QH-QH scan as "z in [-6e-3, +7e-3]". The journal's
own Fig 15 axis runs -5e-3 to 20e-3 with the linking number changing at about -7e-3 and +9.5e-3, and the conference Fig 8
runs about -8e-3 to +24e-3 with changes at about -7.5e-3 and +9e-3, so the digest's range does not match either figure.
The same quantity as a sign-change interval is roughly [-7.5e-3, +9.5e-3].

Journal Fig 15 and conference Fig 8 show the same linking-number history, which suggests they are the same computation
at the same energy; the paper itself does not say so. The journal's rounded "3.15" is consistent with 3.1460717
rounded to two decimals.

## 7. Bibliography (full, as printed in the paper)

1. W. S. Koon, M. W. Lo, J. E. Marsden, and S. D. Ross, "The Genesis trajectory and heteroclinic connections," AAS/AIAA Astrodynamics Specialist Conference, Citeseer, 1999, pp. 99-451.
2. V. Angelopoulos, "The ARTEMIS mission," The ARTEMIS mission, pp. 3-25, Springer, 2010.
3. S. B. Broschart, M.-K. J. Chung, S. J. Hatch, J. H. Ma, T. H. Sweetser, S. S. Weinstein-Weiss, and V. Angelopoulos, "Preliminary trajectory design for the Artemis lunar mission," Advances in the Astronautical Sciences, Vol. 135, No. 2, 2009, pp. 1329-1343.
4. M. Smith, D. Craig, N. Herrmann, E. Mahoney, J. Krezel, N. McIntyre, and K. Goodliff, "The Artemis program: an overview of NASA's activities to return humans to the moon," 2020 IEEE Aerospace Conference, IEEE, 2020, pp. 1-10.
5. G. Gomez, J. Llibre, and J. Masdemont, "Homoclinic and heteroclinic solutions in the restricted three-body problem," Celestial mechanics, Vol. 44, No. 3, 1988, pp. 239-259.
6. G. Gomez and J. Masdemont, "Some zero cost transfers between libration point orbits," Advances in the Astronautical Sciences, Vol. 105, No. 2, 2000, pp. 1199-1216.
7. L. Arona and J. J. Masdemont, "Computation of heteroclinic orbits between normally hyperbolic invariant 3-spheres foliated by 2-dimensional invariant Tori in Hill's problem," Conference Publications, Vol. 2007, American Institute of Mathematical Sciences, 2007, p. 64.
8. R. L. Anderson and M. W. Lo, "Flyby design using heteroclinic and homoclinic connections of unstable resonant orbits," 2011.
9. R. L. Anderson, "Tour design using resonant-orbit invariant manifolds in patched circular restricted three-body problems," Journal of Guidance, Control, and Dynamics, Vol. 44, No. 1, 2021, pp. 106-119.
10. B. Kumar, R. L. Anderson, and R. d. l. Llave, "High-order resonant orbit manifold expansions for mission design in the planar circular restricted 3-body problem," Communications in Nonlinear Science and Numerical Simulation, Vol. 97, 2021, p. 105691.
11. A. F. Haapala and K. C. Howell, "Representations of higher-dimensional Poincare maps with applications to spacecraft trajectory design," Acta Astronautica, Vol. 96, 2014, pp. 23-41.
12. A. F. Haapala and K. C. Howell, "Trajectory design using periapse Poincare maps and invariant manifolds," AAS/AIAA Space Flight Mechanics Meeting, AAS, 2011, pp. 11-131.
13. S. Bonasera and N. Bosanac, "Transitions between quasi-periodic orbits near resonances in the circular restricted three-body problem," AAS/AIAA Astrodynamics Specialist Conference, Lake Tahoe, CA (Virtual), 2020.
14. S. De Smet and D. J. Scheeres, "Identifying heteroclinic connections using artificial neural networks," Acta Astronautica, Vol. 161, 2019, pp. 192-199.
15. Z. P. Olikara, Computation of quasi-periodic tori and heteroclinic connections in astrodynamics using collocation techniques. PhD thesis, University of Colorado at Boulder, 2016.
16. R. C. Calleja, E. J. Doedel, A. R. Humphries, A. Lemus-Rodriguez, and E. Oldeman, "Boundary-value problem formulations for computing invariant manifolds and connecting orbits in the circular restricted three body problem," Celestial Mechanics and Dynamical Astronomy, Vol. 114, No. 1, 2012, pp. 77-106.
17. D. B. Henry and D. J. Scheeres, "A Survey of Heteroclinic Connections in the Earth-Moon System," 73rd International Astronautical Congress (IAC 2022), 2022, pp. C1,9,3,x70258.
18. C. Shonkwiler and D. Vela-Vick, "Higher-dimensional linking integrals," Proceedings of the American Mathematical Society, Vol. 139, No. 4, 2011, pp. 1511-1519.
19. R. L. Ricca and B. Nipoti, "GAUSS' LINKING NUMBER REVISITED," Journal of Knot Theory and Its Ramifications, Vol. 20, No. 10, 2011, pp. 1325-1343.
20. Wikipedia contributors, "Linking number - Wikipedia, The Free Encyclopedia," (no URL or date printed).
21. W. S. Koon, M. W. Lo, J. E. Marsden, and S. D. Ross, "Dynamical systems, the three-body problem and space mission design," Equadiff 99: (In 2 Volumes), pp. 1167-1181, World Scientific, 2000.
22. A. Celletti, "Basics of regularization theory," Chaotic worlds: from order to disorder in gravitational N-body dynamical systems, pp. 203-230, Springer, 2006.
23. Z. P. Olikara and D. J. Scheeres, "Numerical method for computing quasi-periodic orbits and their stability in the restricted three-body problem," Advances in the Astronautical Sciences, Vol. 145, No. 911-930, 2012, pp. 911-930.
24. J. Amanatides and K. Choi, "Ray tracing triangular meshes," Proceedings of the Eighth Western Computer Graphics Symposium, Vol. 43, 1997.

(Accents dropped for plain text; titles as printed, including the oddities in refs 20 and 23.)

## What it does NOT contain

- No numeric mass ratio. No frequencies, rotation numbers, periods, amplitudes (Az, Ax) or rotation of any torus.
- No initial conditions for any torus or connection; no table of connection states; no torus Fourier coefficients.
- No grid sizes N and M, no number of Fourier modes, no number of manifold trajectories, no integration method or
  tolerances, no torus-invariance residual, no corrector tolerance, no fsolve options.
- No connection counts for the Lissajous cases, no numeric scan windows, no z values of the sign changes (only the
  four z values of Fig 7, which are curve snapshots, not connection locations).
- No statement of which libration point each torus belongs to, nor north or south.
- No Sun-Earth or Jupiter-Ganymede cases; no discussion of the quasi-halo to Lissajous discontinuity problem
  (the journal's Section 5); no regularisation equations; no data or code release statement.
- No statement of venue, meeting, location or date in the body.
- No explanation of how xdot is signed after being recovered from the Jacobi integral, and no definition of Omega.

## Inputs this gives a reproduction attempt

- Energy: C = 3.1460717 for every Earth-Moon torus, quasi-halo and Lissajous alike (p. 14). This is 7 digits where the
  journal prints 3.15.
- Manifold offset epsilon = 1e-6 (p. 14); section x = 1 - mu, propagate unstable forward and stable backward (p. 7).
- Torus method: GMOS, Olikara and Scheeres 2012 (p. 5).
- The scan variable is z; curves are in (y, ydot, zdot); xdot from C (pp. 8, 10).
- Quasi-halo to quasi-halo: first crossing, 4 connections, linking number sequence 0, +1, 0, -1, 0, with sign changes
  at about z = -7.5e-3, -1.2e-3, 0, +9e-3 and scanned range about -8e-3 to +24e-3 (Fig 8, p. 11, by eye).
- Fig 7 (p. 11) z snapshots: -0.0075837, -0.00088665, 0.0058104, 0.012507, useful as probes of the y, ydot, zdot
  curves.
- Section torus maps: z map spans about -7e-3 to +27e-3 on a quasi-halo torus, with z = 0 at theta0 near 1.6 and 4.7
  (Figs 3-4, pp. 6, 8, by eye).
- Tori location: one torus at x about 0.82 to 0.88 and one at x about 1.10 to 1.17 (Figs 9, 12, by eye).
- Lissajous connections: use the second crossing for both manifolds (p. 14).
- Initial guess for angles: the triangle-of-intersections recipe (p. 12); the correction is fsolve on four angles with
  Eq 8 (p. 14). An example target for a check: theta guess near (2.23, 3.995) on one torus (Fig 11, by eye).
- Cross-reference: the frequencies 0.2739 and 0.02163 and mass ratio 0.012153643 are in the journal, not here.
