# Digest: Llibre and Martinez Alfaro 1985, ejection and collision orbits of the spatial restricted three-body problem

Date: 2026-10-05 (Sydney). Reading, reasoning and short algebraic checks; no project code was changed.

Source: Jaume Llibre and J. Martinez Alfaro, "Ejection and collision orbits of the spatial restricted three-body problem",
Celestial Mechanics 35:113-128 (1985), DOI 10.1007/BF01227665, received March 1984, accepted July 1984 (Universitat Autonoma de
Barcelona and Universitat de Valencia). Filed in the private paper corpus as
`llibre-martinez-alfaro-1985-ejection-collision-orbits-spatial-restricted-three-body-celest-mech-35-113-doi-10.1007-BF01227665.pdf`.
A 16-page scan with an OCR text layer that garbles the equations (spacing, sub- and superscripts, mu rendered as D, p, U and
similar). I read the text layer for the structure and checked the equations of the regularised system (eqs. 3.1 to 3.6) on page
images (journal pp. 123 and 124). The remaining equations were read from the text layer and cross-checked algebraically (Section 6).

## 1. Relevance verdict

This paper IS about the restricted problem: the circular spatial restricted three-body problem for every mass parameter mu in [0, 1),
with collisions with either primary. That makes it the first held paper whose framework applies to collision with a primary, unlike
the 2019 and 2021 papers of Alvarez-Ramirez et al. (see their digests, which exclude it).

What it gives is qualitative and topological: for each fixed Jacobi constant the set of ejection orbits (and of collision orbits) is a
copy of S^2 x R, and for small mu some of them connect (ejection from one primary to collision with either primary). It prints no
numerical orbit, period, count or tolerance. The usable content is: the regularisation (a McGehee-type blow-up in spherical
coordinates, with the collision manifold written out), the structural statements, the exact mu = 0 skeleton (radial Kepler orbits with
apocentre r = 2/C), and the Jacobi-constant thresholds C = 2 (and C = 0, C = 3 in the figure) at which the picture changes. It is
a result for the small-mass-parameter limit, with the usual perturbative caveat that "mu small enough" is not quantified.

## 2. The model (page 113 to 114)

Spatial circular restricted problem in the synodic frame (rotation frequency 1). The larger primary m1 (mass 1 - mu) is at the origin
and the smaller primary m2 (mass mu) at e_2 = (-1, 0, 0). The Hamiltonian for the massless particle (eq. 1.1), with p conjugate to q:

    H = |p|^2/2 + q2 p1 - q1 p2 - 1/|q| + mu (1/|q| - 1/|q - e_2| - mu/2)

(the constant term is read from the scan and from the form in eq. 3.2, where it appears as -mu/2 inside the mu bracket; the OCR is
ambiguous, so treat the constant as read, not verified). The Jacobi integral is C = -2H. The paper states that its Jacobi constant
differs from the usual one by the constant mu (1 - mu). That is the same offset as in the Barrabes-Mondelo-Olle 2009b digest; it is
relevant to anyone moving a value of C between this paper and `core.cr3bp.jacobi_constant`. The planar problem is the restriction
q3 = p3 = 0. The problem is non-integrable for mu in (0, 1) (cited Moser and Llibre-Simo 1980).

Definitions (page 114): a solution has a collision (ejection) with m1 at t0 if the distance to m1 tends to zero as t tends to t0 from
below (from above); likewise with m2. An orbit on (t-, t+) with an ejection at t- and a collision at t+ is an ejection-collision
(or bi-collision) orbit.

## 3. Main results as printed

Theorem A. For every mu in [0, 1) and every Jacobi constant C, the set of ejection orbits (resp. collision orbits) with m1 or m2 (only
with m1 if mu = 0) of the spatial restricted problem is diffeomorphic to S^2 x R. (The proof says they need not be embedded submanifolds.)

Theorem B (statements as printed):
(i) For every C > 0 with C not equal to 2 there is mu_0 = mu_0(C) such that for every mu in (0, mu_0] there are at least two circles of
    orbits that start in an ejection from the mass 1 - mu, cross one time the hypersurface r = 0 (with r = |q|), and end in a
    collision with the mass 1 - mu. (As printed the hypersurface is r = 0; I cannot tell from the scan whether this is the paper's
    intended hypersurface, since r = 0 is the collision with m1 itself. See the caution in Section 9.)
(ii) For every C < 2 there is mu_0 = mu_0(C) such that for every mu in (0, mu_0] there is at least one orbit that starts in an ejection
     from the mass 1 - mu and ends in a collision with the mass mu without crossing r = 0. By the symmetry S (Section 5) there is also one
     starting from the mass mu and ending at 1 - mu.

They extend the planar results of Llibre (1982, Celest. Mech. 28:83-105), which is their reference [6].

## 4. The regularisation: a McGehee-type blow-up in spherical coordinates (Section 3)

The paper says the binary collision with m1 can be regularised by McGehee variables (Devaney 1981) or by Levi-Civita (Stiefel and
Scheifele 1971). For the spatial problem the usual McGehee variables "do not work", so they use the idea of McGehee as in the planar
paper, with a spherical-coordinate version. This is a blow-up, not Levi-Civita and not KS.

Spherical coordinates (eq. 2.9): q1 = r cos(phi) cos(theta), q2 = r cos(phi) sin(theta), q3 = r sin(phi), with r > 0, theta in [0, 2 pi),
phi in (-pi/2, pi/2), and the conjugate momenta (p_r, p_theta, p_phi). In these coordinates

    H = [p_r^2 + p_phi^2 r^(-2) + p_theta^2 (r cos(phi))^(-2)]/2 - p_theta - 1/r + mu [1/r - (r^2 + 2 r cos(phi) cos(theta) + 1)^(-1/2) - mu/2].

Velocity variables y = r', x1 = r theta' cos(phi), x2 = r phi' (spherical components of the velocity, the paper writes y = r-dot). The
system (3.1) is then no longer Hamiltonian (the energy relation (3.2) defines an invariant hypersurface), with

    r' = y,   theta' = (r cos(phi))^(-1) x1,   phi' = r^(-1) x2,
    y' = r^(-1)(x1^2 + x2^2) + 2 x1 cos(phi) + r cos^2(phi) - r^(-2) + mu [cos(phi) cos(theta) + r^(-2) - (r + cos(phi) cos(theta)) (r^2 + 2 r cos(phi) cos(theta) + 1)^(-3/2)]
    (x1', x2' as printed on page 122-123, with Coriolis terms and mu terms in sin(theta), sin(phi))

and (3.2): H = (x1^2 + x2^2 + y^2 - r^2 cos^2(phi))/2 - 1/r + mu[1/r - r cos(phi) cos(theta) - (r^2 + 2 r cos(phi) cos(theta) + 1)^(-1/2) - mu/2].

Scaled variables v = r^(1/2) y, u1 = r^(1/2) x1, u2 = r^(1/2) x2 and time d t / d tau = r^(3/2). The regularised system (eq. 3.3,
read from the page image; r = 0 is invariant):

    r' = r v
    theta' = (cos(phi))^(-1) u1
    phi' = u2
    v' = v^2/2 + u1^2 + u2^2 - 1 + 2 u1 r^(3/2) cos(phi) + r^3 cos^2(phi) + mu r^2 [cos(phi) cos(theta) + r^(-2) - (r + cos(phi) cos(theta)) (r^2 + 2 r cos(phi) cos(theta) + 1)^(-3/2)]
    u1' = -u1 v/2 - 2 v r^(3/2) cos(phi) + u1 u2 tan(phi) + 2 u2 r^(3/2) sin(phi) + mu r^2 sin(theta) [(r^2 + 2 r cos(phi) cos(theta) + 1)^(-3/2) - 1]
    u2' = -u2 v/2 - r tan(phi) (r^(-1) u1^2 + r^2 cos^2(phi) + 2 u1 r^(1/2) cos(phi)) + mu r^2 sin(phi) cos(phi) [(r^2 + 2 r cos(phi) cos(theta) + 1)^(-3/2) - 1]

(the first mu term in v' contains r^(-2) times r^2, so it is finite at r = 0 when the factors are combined; the printed grouping
is as in the scan). Energy relation: (u1^2 + u2^2 + v^2)/2 - 1 + mu = r H + r^3 cos^2(phi)/2 + mu r [mu/2 + r cos(phi) cos(theta) + (r^2 + 2 cos(phi) cos(theta) + 1)^(-1/2)] (page 123; the last bracket is as printed).

Collision manifold at r = 0 (eq. 3.5): Lambda = {(u1^2 + u2^2 + v^2)/2 = 1 - mu, (theta, phi) in S^2}, a four-dimensional manifold
diffeomorphic to S^2 x S^2. Vector field on Lambda (eq. 3.6):

    theta' = u1/cos(phi),  phi' = u2,  v' = v^2/2 + u1^2 + u2^2 - 1 + mu = (u1^2 + u2^2)/2,
    u1' = -u1 v/2 + u1 u2 tan(phi),   u2' = -u2 v/2 - u1^2 tan(phi).

The flow on Lambda is gradient-like in v: v increases along every non-equilibrium orbit. Equilibria: u1 = u2 = 0, v = +-[2(1 - mu)]^(1/2),
with (theta, phi) arbitrary. These are two spheres of equilibria S^2_+- (not isolated points, unlike the 2019 and 2021 papers),
the physical picture being a radial approach or departure along any direction in space.
- Proposition 3.1: there is an invariant manifold of the equilibrium (theta*, phi*, v = -[2(1 - mu)]^(1/2), u1 = u2 = 0) that is also the
  invariant manifold of (theta*, phi*, v = +[2(1 - mu)]^(1/2), u1 = u2 = 0): the orbits on Lambda connect the lower sphere to the upper one,
  with theta* and phi* related through the first integral cos(alpha) cos(phi) = A (3.8) and solutions (3.9) to (3.10); the arrival direction
  is read off from limits as time tends to plus and minus infinity.
- Proposition 3.3: S^2_+ and S^2_- are normally hyperbolic spheres of equilibria. The linearisation on S^2_+ has eigenvalues lambda,
  lambda, -lambda/2, -lambda/2, 0, 0 with lambda = [2(1 - mu)]^(1/2); on S^2_- the signs are reversed. The zeros are the sphere directions.
  The paper states dim W^s(S^2_+) = dim W^u(S^2_-) = 4 and dim W^u(S^2_+) = dim W^s(S^2_-) = 3 once v is omitted by using the Jacobi
  integral (my reading of the printed sentence; the dimension bookkeeping in the scan is not clear).
- Proof of Theorem A: normal hyperbolicity (Hirsch, Pugh and Shub 1977; a theorem of Devaney type) gives stable and unstable invariant
  manifolds of S^2_+ and S^2_-, each diffeomorphic to S^2 x R; the stable manifold of the right sphere consists of the collision orbits,
  the unstable one of the ejection orbits. For m2, swap the primaries by symmetry.

## 5. Global flow of the integrable case mu = 0 and the small-mu argument (Section 2 and proof of B)

At mu = 0 the problem is the spatial two-body problem in a rotating frame. Integrals (eq. 2.8): C = 2(M3 - h) with h = |p|^2/2 - 1/|q| (sidereal
energy, in the rotating-frame p of eq. 2.5), F = M1^2 + M2^2 = (q2 p3 - q3 p2)^2 + (q3 p1 - q1 p3)^2, and ||M||^2 = F + M3^2.
Printed bifurcation set in the space (C, h, F) (page 121): F = 0, h = 0 and F = -(h + C/2)^2 - 1/(2h), corresponding to the planar
problem, parabolic orbits and circular orbits. Topology of the invariant sets in Tables I (M3 not zero) and II (M3 = 0); Table II, for
the case M = 0 (zero sidereal angular momentum, so all orbits are ejection-collision for C > 0): for F = 0 the set I_Ch0 is S^2 x R and its
orbits are ejection-collision if C > 0; for C = 0 it is two copies of S^2 x R (ejection-parabolic and parabolic-collision); for C < 0
two copies (ejection-hyperbolic, hyperbolic-collision). The ejection orbits and the collision orbits coincide at mu = 0 for C > 0 because the
zero-angular-momentum Kepler orbits go out and back along the same ray.

Fixed C > 0, mu = 0: the first cut of the collision set with the hypersurface v = 0 is, from (3.11) and (3.12),

    {r = -1/H, v = 0, u1 = -r^(3/2) cos(phi), u2 = 0, (theta, phi) arbitrary}  (eq. 3.13),

diffeomorphic to S^2. With H = -C/2 this is r = 2/C, the apocentre of the radial Kepler orbit; v^2/2 - 1 = r H (eq. 3.12) is the radial
energy relation.

Proof of B, as printed: for C > 0 and mu small and positive the cuts gamma^s and gamma^u remain diffeomorphic to S^2 and close to
(3.13). The symmetry S(r, theta, phi, v, u1, u2, t) = (r, -theta, phi, -v, u1, -u2, -t) maps gamma^u onto gamma^s, so for small mu
the two meet at least at points with theta = 0, v = 0, u2 = 0, with u1^2/2 = 1 - mu + rH + r^3 cos^2(phi)/2 + mu r[...] (eq. 3.14); at
mu = 0 this reduces to u1^2 = r^2 cos^2(phi), and gamma^s intersect gamma^u contains at least two circles (the paper's reason for circles is that the intersections of two close spheres are circles; see the caution in Section 9). Because the small primary sits at
(r = 1, theta = pi, phi = 0), the case C near 2 (apocentre r = 2/C = 1) must be avoided in (i). Part (ii): if C < 0 (hyperbolic) or
C = 0 (parabolic) at mu = 0 the ejection orbit W^u goes to infinity, crossing the position (-1, 0, 0) of m2 for the right direction; for
0 < C < 2 one has -1/H = 2/C > 1, so the apocentre lies beyond the m2 position and the radial orbit in the direction of m2
passes through it. So the mechanism is: a radial Kepler orbit about m1 in the direction of m2, long enough to reach m2, persists for small mu
as an orbit from ejection at m1 to collision at m2.

## 6. Printed numbers usable as sourced tests

The paper prints no numerical values of orbits; the usable items are exact statements, each stated or directly derived from printed text:

1. Collision manifold radius: (u1^2 + u2^2 + v^2)/2 = 1 - mu at r = 0; equilibrium speeds v = +-[2(1 - mu)]^(1/2); eigenvalues lambda, lambda,
   -lambda/2, -lambda/2, 0, 0 with lambda = [2(1 - mu)]^(1/2) on S^2_+ (reversed on S^2_-).
2. Monotonicity: v' = (u1^2 + u2^2)/2 >= 0 on r = 0; independent algebra check: with (3.5), v^2/2 + u1^2 + u2^2 - 1 + mu = (u1^2 + u2^2)/2.
3. mu = 0 apocentre of a radial (zero angular momentum) ejection-collision orbit: r = 2/C for C > 0, from r = -1/H and H = -C/2.
   Consistent with the radial Kepler orbit: energy h = -C/2 gives apocentre -1/h = 2/C. This is an exact analytic control for a
   regularised propagator in the rotating frame.
4. The threshold C = 2: apocentre 2/C equals 1, the distance to m2; Theorem B(i) excludes C = 2, Theorem B(ii) applies for C < 2.
5. Bifurcation relation F = -(h + C/2)^2 - 1/(2h): my independent derivation from the integrals (2.8) and the circular orbit condition
   (circular iff h = -1/(2 ||M||^2)) gives ||M||^2 = -1/(2h), so F = ||M||^2 - M3^2 = -1/(2h) - (h + C/2)^2 with M3 = h + C/2. This agrees with
   the printed formula.
6. Figure 3 (C values at which the topology changes): C > 3, C = 3, 0 < C < 3, C = 0, C < 0. I derived C = 3 at mu = 0 as the circular
   orbit of radius 1 with M3 = 1 and h = -1/2, using C = 2(M3 - h) = 3 (my derivation; the paper prints only the figure labels).
7. The Jacobi-constant offset mu (1 - mu) relative to the usual constant (printed statement on page 113).
8. Counts: Theorem B(i) at least two circles, B(ii) at least one orbit for each admissible (C, mu in (0, mu_0(C)]). The values mu_0(C) are
   not printed and not computed.

No numerical value of mu_0, no initial condition, no period.

## 7. Relation to second-species arcs and to Broucke's collision orbits

Henon 1968 (digest `2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`): second-species arcs run between consecutive
collisions with the small primary, and for small mu at mu = 0 they become Keplerian arcs about the large primary that pass through the
position of the small one. The Llibre and Martinez Alfaro mu = 0 skeleton is narrower: ejection-collision orbits with respect to the LARGE
primary only, with zero sidereal angular momentum (radial arcs). Theorem B(ii) is the nearest overlap: a radial arc that leaves m1 and
reaches m2's position (when 2/C > 1) persists as an m1-to-m2 collision orbit. A generic Henon arc has nonzero angular momentum about
m1, so is not in this skeleton, but the same machinery (a collision with the small primary regularised by a blow-up, and a continuation
from mu = 0) applies to it in principle; the paper does not treat it. Related held digests: Gomez-Olle 1991 (second species circular and elliptic),
Hitzl-Henon 1977 and 1977b, Perko 1976 to 1981, Font-Nunes-Simo 2002 and 2009 (consecutive quasi-collisions), Brjuno 1978.

Broucke 1969 (digests `2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-a.md` and part b): I did not re-read them for this
digest; the periodic collision orbits there belong to the elliptic restricted problem. The structural link is only the same:
collision orbits as a skeleton from which families of periodic orbits are built (a one-parameter family of ejection orbits, closed by
symmetry S). The symmetry S of this paper, (theta, v, u2, t) -> (-theta, -v, -u2, -t), is the usual time-reversal symmetry that makes symmetric
periodic collision orbits come in families.

## 8. Techniques applicable to the project's problems

### 8.1 `#928` (regularised propagator for close passes)

- Alternative or complement: complement. The paper regularises collision with a primary by a McGehee-type blow-up, and states that
  Levi-Civita is the other known way (and that for the spatial problem the usual McGehee variables fail). It does not compare them numerically.
  For a propagator the blow-up gives a flow on and near r = 0 with spheres of equilibria; Levi-Civita (planar) or KS (spatial) gives a smooth
  passage with no equilibria. The two serve different tasks: the blow-up classifies what ejects and collides, KS or Levi-Civita integrates through
  a pass.
- Collision-orbit test case with an exact answer, no numerics needed from the paper: at mu = 0 with zero angular momentum, an orbit ejected
  from m1 reaches apocentre r = 2/C and falls back (in the rotating frame, with the Coriolis term); a regularised propagator started on the radial
  ejection direction should reproduce the apocentre r = 2/C for C > 0 and a symmetric return, to the integrator's tolerance. This is the
  same family as the project's radial-fall closed-form control, but with the Jacobi constant normalisation of the paper (C = -2H, offset mu (1 - mu)).
  At mu > 0 there is no printed value to test against.
- Equation (3.3) is a possible alternative formulation to compare with `core/cr3bp_regularized.py` (Sundman only) if a McGehee-type
  integrator were ever wanted, but the spherical-coordinate singularity at phi = +-pi/2 (the tan(phi) terms) and the r^(-2) mu-term grouping
  make it unattractive as a propagator.

### 8.2 `#899` (ejection-collision orbits as seeds or skeleton)

- The paper's skeleton is exactly what `#899` needs at mu = 0: for each C > 0 the zero-angular-momentum radial orbits form the sphere S^2 of
  directions (apocentre 2/C). For small mu the orbits connecting ejection from one primary to collision with the other persist (Theorem B(ii)), starting from a
  radial orbit in the direction of the small primary with 2/C > 1.
- Seed construction implied (my reading, not in the paper): for a chosen C < 2 and small mu, take the radial Kepler orbit from m1 along the
  direction of m2 and continue the connecting orbit in mu by correcting the ejection direction; the symmetry S gives the converse orbit.
- Limit of the paper: only existence for mu in (0, mu_0(C)], with mu_0 not quantified, and nothing for the generic (nonzero-angular-momentum) second-species arcs
  of Henon 1968. For the latter, the restricted-problem papers of Olle, Rodriguez and Soler (listed in the `#899` entry; not held) are the
  continuation of this line; this paper is cited there as the spatial predecessor (Llibre and Martinez Alfaro are in the `#899` entry's reference list).

### 8.3 `#906` (demanded-turn gate)

A collision orbit has zero angular momentum about the primary it collides with, so the flyby geometry that the demanded-turn gate assumes
(a defined incoming and outgoing asymptote with an impact parameter) degenerates: a collision with m2 has zero periapsis and a defined direction but
the hyperbolic turn is 180 degrees on a radial pass or undefined, in line with the existing rule that a demanded turn of 0 or 180 degrees
is not an encounter and that an undefined turn is "indeterminate". This paper adds one concrete item: at mu = 0 and C < 2 the m1-m2 connecting orbit
is radial, and the whole sphere S^2 of directions is a degenerate set of such cases, so any seed generated from the ejection skeleton must be routed
to "indeterminate" by the gate, not rejected. No number is supplied.

## 9. Cautions and follow-ups

Cautions:
- OCR: equations on pages 114 to 117 were not viewed on images; the topology tables (I, II) and Fig. 3 are read from the text layer, whose figure labels
  are partly garbled. The key equations (3.1 to 3.6) were checked on images.
- Theorem B(i), "cross one time the hypersurface r = 0", is odd as printed because r = 0 is the collision with m1 itself. It may denote another hypersurface
  in the original; do not build on this clause.
- The proof that gamma^s and gamma^u intersect in at least two circles is a sketch (two nearby spheres in S^2-like cuts); the paper does not
  discuss transversality.
- The normal-hyperbolicity statement and the dimension counts of W^s and W^u are given as printed; the scan's bookkeeping is not fully clear.
- Reference numbering in the printed list is alphabetical without numbers in the scan; my attribution of [3], [6], [10] in the text is by position.

Follow-ups (no task numbers registered):
- Llibre, "On the restricted three-body problem when the mass parameter is small", Celest. Mech. 28:83-105 (1982): the planar predecessor; not in the
  corpus per the CORPUS_INDEX grep (no Llibre 1982 row); DOI not printed in this paper.
- Llibre and Simo 1980, "Oscillatory solutions in the planar restricted three-body problem", Math. Ann. 248:153-184; Devaney 1981, "Singularities
  in classical mechanical systems" (Birkhauser, Ergodic Theory and Dynamical Systems I, 211-333), the source of the McGehee blow-up used here; Devaney 1978,
  Invent. Math. 45:221-251 (collision orbits in the anisotropic Kepler problem). None is held; no DOIs are printed in this paper.
- The Olle, Rodriguez and Soler papers (DOIs in the `#899` entry) for the restricted-problem continuation and the McGehee-based spatial paper.
- Stiefel and Scheifele 1971 is held (digest `2026-10-04-digest-stiefel-scheifele-1971-linear-regular-celestial-mechanics.md`); Kustaanheimo-Stiefel 1965 is held and
  digested (`2026-10-05-digest-kustaanheimo-stiefel-1965-ks-regularization.md`). Together they cover the Levi-Civita and KS side of the comparison in 8.1.

## 10. Note added 2026-10-05: the planar predecessor (Llibre 1982)

Llibre 1982 (Celest. Mech. 28:83-105, DOI 10.1007/BF01230662; digest `2026-10-05-digest-llibre-1982-restricted-problem-small-mu.md`) is the planar
version of this paper. Its Theorem B prints the section as r-dot = 0 (cross one time the surface r-dot = 0), which settles the oddity flagged in
Section 9 above: the "r = 0" of Theorem B(i) here is most likely r-dot = 0 (the scan or the printing lost the dot), so those orbits cross the section
r-dot = 0, not the collision. Its Theorem B(ii) says "one and only one" orbit from m1 to m2, where this paper says "at least one". The planar paper carries
the exact mu = 0 phase portrait with printed, verified numbers (circular radii, boundary radii and energies at C = 3.25, 3.1, 3, 2^(5/3), 1, 0, -1) that
this paper only summarises. The Rodriguez del Rio thesis (part B digest) contradicts this paper's claim of families of spatial symmetric orbits.
