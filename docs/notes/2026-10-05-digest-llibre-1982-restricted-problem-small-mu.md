# Digest: Llibre 1982, on the restricted three-body problem when the mass parameter is small

Date: 2026-10-05 (Sydney). Reading, and an independent numerical check of every closed-form value printed in the paper's Figure 1
(Section 5); no project code was changed.

Source: Jaume Llibre, "On the restricted three-body problem when the mass parameter is small", Celestial Mechanics 28:83-105 (1982),
DOI 10.1007/BF01230662 (Universitat Autonoma de Barcelona; paper presented at the 1981 Oberwolfach conference on mathematical methods in
celestial mechanics). Filed in the private paper corpus as
`llibre-1982-restricted-three-body-problem-mass-parameter-small-celest-mech-28-83-doi-10.1007-BF01230662.pdf`. A 23-page scan with an OCR
text layer that garbles equations; I read the equations of Section 3 from the text layer with the regularised system cross-checked against
the 1985 spatial paper (whose equations I viewed on page images), and viewed Figure 1b on a page image.

## 1. Verdict

The paper is the planar predecessor of Llibre and Martinez Alfaro 1985 (digest
`2026-10-05-digest-llibre-martinez-alfaro-1985-ejection-collision-spatial-rtbp.md`) and is directly about the circular planar restricted
problem at small mu. It has two halves:

- Sections 1 to 3: the global flow of the two-body rotating problem (mu = 0), and from it the ejection-collision, collision-to-other-primary
  and parabolic orbits for small mu > 0 (Theorems A, B, C). A McGehee-type blow-up (collision manifold a torus) is the regularisation.
- Section 4: Theorem D, that for small mu all of the bounded phase space except a set of arbitrarily small Lebesgue measure lies in the
  closure of the symmetric periodic orbits (ten classes by pairs of right-angle crossings). Twist-map and KAM argument; of limited use to
  the project except as a statement of the structure that persists from mu = 0.

For the project: the useful items are (a) the exact mu = 0 skeleton for each Jacobi constant (a table of the sets I_hC and ring boundaries in
the surface of section), with printed numerical radii and energies that I verified to the last printed digit (a rare sourced-test source);
(b) the regularisation written out; (c) existence of two m1-to-m1 orbits and a unique m1-to-m2 orbit for small mu.
It prints no orbit initial condition and no value of mu_0(C).

## 2. Model and conventions (page 83; Section 4 changes them)

Section 1 (and the collision analysis): synodic frame, rotation frequency 1, larger primary m1 (mass 1 - mu) at the origin and the smaller
primary m2 (mass mu) at e_2 = (-1, 0). Hamiltonian (eq. 1.1), as in the 1985 paper:

    H = |p|^2/2 + q2 p1 - q1 p2 - 1/|q| + mu (1/|q| - 1/|q - e_2| - mu/2)   (the constant is read from the scan and from eq. 3.3)

with C = -2H, and the paper's Jacobi constant differs from the usual one by mu (1 - mu) (printed; same offset as the 1985 paper).

Section 4 (symmetric periodic orbits) changes the frame: m1 = 1 - mu at e_1 = (mu, 0), m2 = mu at e_2 = (mu - 1, 0), with
H = |p|^2 + q2 p1 - q1 p2 - (1 - mu)/|q - e_1| - mu/|q - e_2| as printed (the kinetic term has no factor 1/2 in the scan). The two
conventions are not the same frame; do not mix them.

Symmetry (page 84): (q1, q2, p1, p2, t) -> (q1, -q2, -p1, p2, -t). In (r, v, theta, u, tau) coordinates (Section 3) it reads
S : (r, v, theta, u, tau) -> (r, -v, -theta, u, -tau).

## 3. The two-body rotating problem, mu = 0 (Section 2)

Integrals: sidereal energy h = |p|^2/2 - 1/|q| and sidereal angular momentum M = q1 p2 - q2 p1, with H = h - M = -C/2. In polar
coordinates (q1 + i q2 = r exp(i theta)): for M not zero, r = M^2 [1 + (1 + 2 h M^2)^(1/2) cos(theta - theta')]^(-1) (eq. 2.1);
for M = 0 (collision orbits) r-dot^2 = 2/r + 2h (eq. 2.2). The sets I_hC are plotted in the half-plane (r, r-dot) for seven regimes of C
(Fig. 1a to 1g), reducing to five qualitative pictures: C > 3, C = 3, 0 < C < 3, C = 0, C < 0. Circular orbits have radii r_i, the positive
roots of the polynomial (eq. 2.3, citing Birkhoff)

    4 r^3 - C^2 r^2 + 2 C r - 1 = 0,

with h_i = -1/(2 r_i). The collision cylinder has h_c = -C/2 and r_c = -1/h_c = 2/C when h_c < 0; the parabolic cylinder meets r-dot = 0
at r_p = C^2/8; the boundary radii r_a, r_b (named r_i and r_o in the figures) are the two positive roots of r^3 - C r + 2 = 0. For C > 3 the
set I_C has two components, I^1_C (h in [h_1, h_2]) and I^2_C (h in [h_3, infinity)); for C below 3 the cubic has only one positive root and I_C has one
component. Topologies (Tables I to V): circular retrograde orbit (S^1) at h = h_1, tori of retrograde ellipses (S^1 x S^1), the cylinder of elliptic
collision orbits at h = h_c (S^1 x R), tori of direct ellipses, a circular direct orbit; for C > 0 the stable and unstable collision sets coincide (elliptic collision
orbits), for C = 0 two copies of S^1 x R of parabolic collision orbits, for C < 0 two copies of hyperbolic collision orbits (Tables IV and V). Levi-Civita regularisation
of the binary collision identifies the boundary of the collision cylinder, turning the open solid torus I^1_C into a 3-sphere for C > 3 (and for
0 < C < 3 an open solid torus). All collision orbits at mu = 0 are retrograde.

Table VI (page 100): the rings of the surface of section Sigma (r-dot = 0, r-ddot > 0), with inner and outer radius, by regime:

| orbits | C | inner radius | outer radius |
| --- | --- | --- | --- |
| retrograde | (-infinity, 0) | C^2/8 parabolic | r_1 circular |
| retrograde | 0 | 0 parabolic collision | r_1 circular |
| retrograde | (0, infinity) | 0 elliptic collision | r_1 circular |
| direct | (0, 3) | 0 elliptic collision | C^2/8 parabolic |
| direct | 3 | 0 elliptic collision | r_2 = r_3 = 1 circular |
| direct | 3 | r_2 = r_3 = 1 circular | C^2/8 parabolic |
| direct | (3, 32^(1/3)) | 0 elliptic collision | r_2 circular |
| direct | (3, 32^(1/3)) | r_0 elliptic | C^2/8 parabolic |
| direct | (3, 32^(1/3)) | r_0 elliptic | r_3 circular |
| direct | [32^(1/3), infinity) | 0 elliptic collision | r_2 circular |
| direct | [32^(1/3), infinity) | C^2/8 parabolic | r_3 circular |

(The Table VI layout is read from the text layer; verify on a page image before building a test on it. Here 32^(1/3) is the printed threshold
(32)^(1/3) = 2^(5/3), about 3.1748, at which r_0 meets r_p = C^2/8.)

## 4. Collision orbits, regularisation and the small-mu results (Section 3)

### 4.1 Regularisation (McGehee-type blow-up, planar)

The paper says the usual McGehee variables fail because H is not of the form |p|^2/2 + V(q) with V homogeneous, but the ideas of McGehee
work. Polar coordinates (Q1 = r, Q2 = theta, momenta P1, P2); the Hamiltonian becomes

    H = (P1^2 + P2^2 Q1^(-2))/2 - P2 - 1/Q1 + mu [1/Q1 - (Q1^2 + 1 + 2 Q1 cos Q2)^(-1/2) - P1 sin Q2 - P2 Q1^(-1) cos Q2],

(the printed first form carries a mu-dependent bracket; read from the scan). Velocity variables y = r-dot, x = r theta-dot; system (3.2):

    r' = y,  y' = x^2/r + 2x + r - r^(-2) + mu [cos(theta) + r^(-2) - (r + cos(theta)) (r^2 + 1 + 2 r cos(theta))^(-3/2)],
    theta' = x/r,  x' = -x y/r - 2y + mu sin(theta) [(r^2 + 1 + 2 r cos(theta))^(-3/2) - 1],

energy relation (3.3): H = (x^2 + y^2 - r^2)/2 - 1/r - mu [mu/2 + r cos(theta) + (r^2 + 1 + 2 r cos(theta))^(-1/2) - 1/r].
Scaled velocities u = r^(1/2) x, v = r^(1/2) y, time dt/dtau = r^(3/2). The regularised system (3.5):

    r' = r v
    v' = v^2/2 + u^2 - 1 + mu + 2 u r^(3/2) + r^3 + mu r^2 [cos(theta) - (r + cos(theta)) (r^2 + 1 + 2 r cos(theta))^(-3/2)]
    theta' = u
    u' = -u v/2 - 2 r^(3/2) v + mu r^2 sin(theta) [(r^2 + 1 + 2 r cos(theta))^(-3/2) - 1]

with energy relation (3.6): (u^2 + v^2)/2 - 1 + mu = r H + r^3/2 + mu r [mu/2 + r cos(theta) + (r^2 + 1 + 2 r cos(theta))^(-1/2)].
(The scan loses some terms in the first line of v'; I compared with the 1985 paper's page-image equations, where the spatial line reduces to this one at phi = 0, u2 = 0.
The exact grouping of the mu terms is as read, not verified.)

Collision manifold at r = 0 (3.7): (u^2 + v^2)/2 = 1 - mu with theta arbitrary, a two-dimensional torus; vector field on it (3.8):
v' = u^2/2, theta' = u, u' = -u v/2. Equilibria: two circles u = 0, v = +-[2(1 - mu)]^(1/2), theta arbitrary; every other orbit goes from the lower
circle to the upper one (gradient-like in v), the same picture as the Kepler problem (Devaney). Hence (Theorem A) the set of orbits ending
(beginning) at collision with m1 is a cylinder, and with m2 the same by symmetry.

At mu = 0 the collision orbits have M = 0, the energy relation reduces to v^2/2 - 1 = r H, and u = -r^(3/2) (minus sign: collision orbits are
retrograde). The cylinder of collision orbits for H < 0 is plotted in Fig. 6 in cylindrical coordinates (r, theta, v).

### 4.2 Results

Theorem A: for each Jacobi constant, the set of orbits that end or begin at collision with m1 or m2 is topologically a cylinder. (The 1985
paper generalises this to S^2 x R in space.)

Theorem B: for mu small enough,
(i) for each C > 0 there is mu_0 = mu_0(C) such that for every mu in (0, mu_0] there are at least two orbits that leave the collision with
    mass 1 - mu, cross the surface r-dot = 0 one time and go to collision with mass 1 - mu again;
(ii) for each C < 2 there is mu_0(C) such that for every mu in (0, mu_0] there is one and only one orbit that leaves the collision with mass 1 - mu
     and goes to collision with the mass mu without crossing r-dot = 0; by the symmetry, one and only one from mass mu to mass 1 - mu.
(The "surface r = 0" of the 1985 paper's Theorem B is r-dot = 0 here, which resolves the oddity I flagged in that digest: the 1985 statement
inherits a typo or OCR loss, and the surface is the section r-dot = 0.) Note that (ii) here says "one and only one", while the 1985 spatial paper
says "at least one".

Proof idea (page 97 to 98): for H < 0 fixed, let gamma^u and gamma^s be the first cuts of the unstable and stable collision cylinders with the
plane v = 0. At mu = 0 they coincide and form a circle of radius -H^(-1) = 2/C in the (r, theta, u) space. For small positive mu both are real
analytic simple curves, not equal, intersecting at the points theta = 0 and theta = pi nontangentially (from the symmetry S, S(gamma^u) = gamma^s). Each
intersection is an ejection-collision orbit with m1; that gives (i). For (ii): if C < 2 then 2/C > 1, so the circle of radius 2/C encloses the circle
r = 1 on which the m2 position lies; for mu small, gamma^u and gamma^s remain near that circle, so the orbit through the corresponding point reaches m2.
For C <= 0, W^u goes to infinity without crossing r-dot = 0 at mu = 0, and for small mu its projection crosses the circle r = 1.
The paper adds that, with the Poincare map on the section v = 0, "we have put the Bernoulli shift as a subsystem" and given a geometrical
interpretation, in a companion not yet published ([7], Llibre and Pinyol, "Collision orbits in the planar circular restricted three-body problem",
to appear); I have not seen it.

Theorem C (parabolic orbits, McGehee 1973 type result: parabolic orbits form a cylinder): for mu small,
(i) for each epsilon > 0 and all C outside (-epsilon, epsilon) there is mu_0(C) such that for every mu in (0, mu_0] at least two parabolic orbits leave
    infinity, cross r-dot = 0 once, and reach infinity;
(ii) for each small epsilon > 0 and all C in (-8^(1/2) + epsilon, 8^(1/2) - epsilon), there is mu_0(C) such that exactly one parabolic orbit leaves infinity and
     goes to collision with m2 without crossing r-dot = 0 (and, by symmetry, one from m2 to infinity). The bound 8^(1/2) is the value where
     the parabolic circle r_p = C^2/8 reaches radius 1 (the paper notes that for C in (-8^(1/2), 8^(1/2)) the parabolic cylinder cuts r-dot = 0 in a circle of
     radius less than 1). For C large and mu small the parabolic first cuts intersect in only two points, nontangentially (Llibre-Simo 1980,
     Theorem 5.1 of that paper).

## 5. Printed numbers usable as sourced tests, and my verification

All of Figure 1 values below were reproduced by my independent computation (roots of the quartic of eq. 2.3 and of r^3 - C r + 2; the
sidereal energy h = r^2/2 - 1/r at the boundary radii, where the radial velocity vanishes and M = r^2; consistency C = 2(M - h)). The
agreement is to the last printed digit, which also tests the transcription of eq. 2.3 and of the Jacobi-constant convention C = 2(M - h) = -2H.

| Fig. | C | printed values | my check |
| --- | --- | --- | --- |
| 1a | 3.25 | r_1 = 0.2367865, r_2 = 0.5783759, r_3 = 1.825462 (circular radii); r_i = 0.7401394, r_o = 1.3149066; h_i = -1.0771937, h_0 = 0.1039795 | roots of 4r^3 - C^2 r^2 + 2 C r - 1: 0.23678647, 0.57837587, 1.82546266; roots of r^3 - C r + 2: 0.7401394, 1.31490664; h at r_i, r_o: -1.0771937, 0.1039794 |
| 1b | 32^(1/3) = 2^(5/3) | r_1 = 2^(-5/3)(3 - 5^(1/2)), r_2 = 2^(-2/3), r_3 = 2^(-5/3)(3 + 5^(1/2)), r_i = 2^(-2/3)(5^(1/2) - 1), r_0 = 2^(1/3), h_i = 2^(-1/3)(1 - 5^(1/2)), h_0 = 2^(2/3) (as printed) | all radii and h_i reproduce (0.240624, 0.629961, 1.649258, 0.778674, 1.259921, -0.981068); r_0 = r_p = C^2/8 = 2^(1/3); at r_0 my h is 0, not 2^(2/3): the printed h_0 for this case does not reproduce (see below) |
| 1c | 3.1 | r_1 = 0.2445555, r_2 = 0.7022517, r_3 = 1.4556927, r_i = 0.8288294, r_0 = 1.1933108, h_i = -0.8630418, h_0 = -0.1260092 | 0.24455552, 0.70225175, 1.45569273; 0.82882939, 1.19331083; h: -0.86304182, -0.12600932 |
| 1d | 3 | r_1 = 0.25, r_2 = r_3 = r_i = r_0 = 1, h_i = h_0 = -0.5 | double root at 1 and r_1 = 0.25 (roots 0.25 and 1, 1); h(1) = -1/2 |
| 1e | 1 | r_1 = 0.4320408 | root 0.4320408 |
| 1f | 0 | r_1 = 2^(-2/3) | 0.629961 |
| 1g | -1 | r_1 = 1 | root 1 |

My reading of the h_0 discrepancy in Fig. 1b: at the threshold C = 2^(5/3) the outer boundary radius r_0 = 2^(1/3) coincides with the parabolic
radius r_p = C^2/8, where h = 0 (verified: r^2/2 - 1/r = 2^(-1/3) - 2^(-1/3) = 0). The printed "h_0 = 2^(2/3)" is, by my computation,
the sidereal angular momentum M = r_0^2 = 2^(2/3), or a slip for 0; I did not find the page image clear enough to settle which. Do not use
the printed h_0 of Fig. 1b as a test value; use the other entries.

Further sourced statements: thresholds in C at which the topology changes (C = 3, C = 32^(1/3), C = 0, and for Theorem C the values +-8^(1/2)); apocentre of
the mu = 0 collision cylinder r = 2/C (printed as -H^(-1)); the collision-manifold energy radius (u^2 + v^2)/2 = 1 - mu and equilibrium speeds
+-[2(1 - mu)]^(1/2) (the planar analogue of the 1985 values); the Jacobi-constant offset mu (1 - mu); the closed form for the twist map of Section 4
(Section 6 below).

Counts: Theorem B(i) at least two orbits, B(ii) exactly one, Theorem C(i) at least two, C(ii) exactly one. No mu_0(C) printed.

## 6. Theorem D and the twist-map construction (Section 4)

Theorem D: for any fixed C and any epsilon > 0 there is mu_0 > 0 such that, for mu in [0, mu_0], the set of bounded orbits not in the closure
of one of the ten classes C_ij of symmetric periodic orbits has Lebesgue measure less than epsilon. The classes are the pairs of right-angle crossings of
the q1 axis among four types: C_1 lower passage at opposition (x > mu, y = 0, r-dot = 0, r-ddot > 0), C_2 higher passage at opposition (r-ddot < 0),
C_3 lower and C_4 higher passage at conjunction (x < mu); ten combinations C_ij = C_i intersect C_j. The topology on the surface of section is that
induced from (R^2 minus {0, e_2}) x R^2 at fixed C. The surface of section Sigma = {r-dot = 0, r-ddot > 0} has one connected component for C <= 3, two for C > 3.

At mu = 0 on a ring R of Table VI the Poincare map is a twist map (eq. 4.1): r' = r, theta' = theta + f(r), with f(r) = -2 pi a^(3/2) (the
mean motion in the rotating frame, as printed) where the semimajor axis a = a(r) is given implicitly by (eq. 4.2)
r = a {1 - [1 - (1 - C a)^2 (4 a^3)^(-1)]^(1/2)}, with |df/dr| > 0 on the ring (printed). Lower crossings are the segments theta = 0 and theta = pi
on the ring; higher crossings are the arcs gamma_2: theta = pi (1 - a^(3/2)) and gamma_4: theta = -pi a^(3/2), with a = a(r), per Birkhoff. For mu > 0 the map is
r' = r + mu g_1, theta' = theta + f(r) + mu g_2 (analytic, periodic) and area-preserving; the Kolmogorov-Arnold-Moser theorem (Moser's twist form, s = 5) gives
invariant curves leaving out arbitrarily small measure (Arnold, Theorem 4.4), and a density lemma (Lemma 4.6, Theorem 4.7) puts symmetric periodic points
arbitrarily close to every invariant curve. The proof at mu = 0 uses Jacobi's theorem on irrational rotations of the circle (Lemma 4.1 and 4.2).

## 7. Relation to the 1985 spatial paper, to Henon, to Broucke, and to the Rodriguez del Rio thesis

- Llibre and Martinez Alfaro 1985 (spatial): Theorem A, S^2 x R in place of the cylinder; Theorem B(i) and (ii) the same structure with "at least" in (ii); the
  collision manifold S^2 x S^2 replaces the torus; the same two-sphere (here two-circle) equilibria with v = +-[2(1 - mu)]^(1/2); the apocentre r = -1/H = 2/C is the same skeleton.
  The planar paper is the one with the exact mu = 0 phase portrait (Figure 1 and Tables I to VI), which the spatial one only summarises (its Table II and
  Fig. 3, C thresholds 3 and 0). The note on the "r = 0" typo above applies.
- Henon's second-species arcs (digest `2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`): the m1-to-m2 orbit of Theorem B(ii) is the
  radial (zero angular momentum about m1) member of the same objects Henon's arcs generalise: Henon's arcs leave and return to m2 along Keplerian ellipses about m1
  with nonzero angular momentum, whereas Llibre's mu = 0 skeleton is only the M = 0 collision cylinder at h = -C/2. So this paper does not cover generic second-species
  arcs; it covers the degenerate (rectilinear) end of the family. At mu = 0 and C < 2 the apocentre 2/C exceeds 1, which is the radial condition for reaching m2.
- Broucke 1969 (digests `2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-a.md` and part b): I did not re-read them; the structural link is
  the symmetric periodic orbit classes of Theorem D (right-angle crossings, symmetric under S), the same construction that Broucke's symmetric periodic
  orbits use in the elliptic problem.
- Rodriguez del Rio 2021 thesis (digest `2026-10-05-digest-rodriguez-del-rio-2021-thesis-part-b.md`; part A is not yet present as I write): the thesis
  lists this paper in its bibliography ("Not in the index" at that time) and treats the planar n-ejection-collision orbits much more completely, with
  a quantitative numerical method and with the Hill problem. Its part B digest also notes that the 1985 spatial paper's analysis ignored the z axis and
  that its claim of families of spatial symmetric ECOs is contradicted by the thesis (the planar symmetric ECOs are isolated and do not continue to spatial
  families). The thesis, not these two 1980s papers, is the source to use for any numerical count (four planar ECOs at mu = 0.1, C = 4, per its digest). The existence
  of four ECOs for every mu is attributed in the thesis to Chenciner and Llibre 1988 (Ergodic Theory Dyn. Syst. 8:63-72), which sharpens Theorem B(i) here
  from "at least two" in the small-mu regime to "four for any mu".

## 8. Techniques applicable to the project's problems

### 8.1 `#928` (collision-orbit test case for the regularised propagator)

- Regularisation comparison: the collision manifold is the torus (u^2 + v^2)/2 = 1 - mu, equilibria on two circles, and every orbit passes from the
  lower to the upper one; the planar Levi-Civita regularisation (the paper cites Stiefel and Scheifele) identifies the ends of the collision cylinder and gives
  a smooth passage with no blow-up. The two are complementary, as concluded for the 1985 paper: blow-up for the asymptotic classification, Levi-Civita or KS for
  the integration through a pass. The paper makes no numerical comparison.
- Exact controls with printed numbers (a better source than the 1985 paper because the planar mu = 0 flow has closed forms): for C > 0 the zero-angular-momentum
  orbit about the origin has apocentre r = 2/C in the rotating frame; for C = 3.25 this is 0.6153846 (my arithmetic, 2/3.25), which lies between the printed
  r_2 = 0.5783759 and r_i = 0.7401394 in Fig. 1a, as drawn (r_c between r_2 and r_i). A mu = 0 regularised propagator started on the radial ejection ray at the
  collision cylinder should return to r = 2/C at v = 0 (u = -r^(3/2) at the apocentre) and complete a retrograde loop back to collision; this tests the Coriolis term
  and the regularisation together, beyond the project's nonrotating radial-fall closed form.
- Printed Fig. 1 radii and energies give circular-orbit controls for the rotating-frame mu = 0 propagator: a circular orbit of radius r_1 at the printed C is a
  fixed point; e.g. C = 3.1 with r_2 = 0.7022517 (direct) and r_1 = 0.2445555 (retrograde). These need no regularisation, so they are weak controls for the pass
  but valid controls for the frame and the Jacobi-constant convention (offset mu (1 - mu) at mu > 0).
- At mu > 0, no printed number: Theorem B gives existence only.

### 8.2 `#899` (ejection-collision orbits as seeds or skeleton)

- The seed for an m1-to-m1 orbit at small mu is a point on the symmetric half-line theta = 0 or theta = pi of the intersection of the first cuts gamma^u and gamma^s with the plane v = 0,
  starting from the mu = 0 circle of radius 2/C. This is a one-dimensional correction in theta (at theta = 0 and pi by symmetry), cheaper than
  the generic two-angle search of the spatial case, and it is the same bracket-by-symmetry the thesis uses for planar ECOs.
- The m1-to-m2 orbit (Theorem B(ii)) is the zero-angular-momentum seed for a second-species-like arc that starts at m1 collision and ends at m2 collision
  for C < 2; unique for small mu. It is the unique "radial" member and is a natural first seed for continuation in mu or C, then relaxed to nonzero angular momentum
  by the Henon families.
- Do not use this paper for counts beyond "at least two" (use the thesis and Chenciner and Llibre: four for any mu), nor for mu_0(C).

## 9. Cautions and follow-ups

- OCR: the equations in Sections 3 and 4 are read from a degraded text layer; the regularised system was rebuilt by comparison with the 1985 paper's
  page images. Table VI and Figures 1a to 1g (other than 1b) were not viewed on page images; the numbers in Section 5 were verified by computation, which
  checks the numbers but not the figure-to-regime assignment.
- The paper presents Theorem B(i)/(ii) with "r-dot = 0" in the scan; my reading of the 1985 "r = 0" as the same section is an inference.
- Fig. 1b h_0 discrepancy: see Section 5.

Follow-ups (no task numbers registered):
- Llibre and Pinol, "Collision orbits in the planar circular restricted three-body problem" (cited as to appear in this paper; the Bernoulli shift for the
  first-cut Poincare map): not located and no DOI printed; the later Llibre-Pinol 1990 paper on the elliptic problem is listed in the thesis digest.
- Llibre and Simo 1980, "Oscillatory solutions in the planar restricted three-body problem", Math. Ann. 248:153-184 (cited here for the parabolic first cuts and
  the density of periodic orbits, Theorem 5.1): not held.
- Gomez and Llibre 1981, "A note on a conjecture of Poincare", Celest. Mech. 24:335-343, the earlier density theorem that Theorem D improves: not held.
- McGehee 1973, "A stable manifold theorem for degenerate fixed points with applications to celestial mechanics", J. Differential Equations 14:70-88:
  the source of the parabolic-orbit cylinder: not held; no DOI printed here.
- Chenciner and Llibre 1988 (Ergodic Theory Dyn. Syst. 8:63-72) and Lacomba and Llibre 1988 (J. Differential Equations 74:69-85): flagged in the thesis digest; not held.
- Devaney 1981 (Singularities in classical mechanical systems) is the source for the collision manifold used here; not held.
