# Digest: Barrabes & Gomez 2002, "Spatial p-q resonant orbits of the RTBP" (#920)

Celestial Mechanics and Dynamical Astronomy 84:387-407 (2002), received 14 Aug 2001, accepted 19 Feb 2002.
Filed in the private paper corpus as `barrabes-gomez-2002-spatial-pq-resonant-orbits-rtbp-cmda-84-387-doi-10.1023-A1021137127909.pdf`
(21 PDF pages; PDF page n is journal page 386+n; text layer present, `pdftotext | wc -w` = 7921).
Status of every statement below: READ (from the page images, all 21 pages read) unless marked COMPUTED (my own short
arithmetic, 2026-10-04) or INFERRED. Page references are journal pages (p.387 to p.407).

An earlier, shorter treatment exists in `docs/notes/2026-07-28-742-anderson-barrabes-franz-earth-moon-jovian-digest.md`
section 3; this note does not repeat it. This note is for programmers: every formula needed for the OUT half of a seed, as printed.
The companion paper (Barrabes & Gomez 2003, digest `docs/notes/2026-10-04-digest-barrabes-gomez-2003-pq-resonant-second-species.md`)
supplies the IN map and the matching; the full seed recipe is in `docs/notes/2026-10-04-digest-casoliva-2008-aiaa-families-cycler-trajectories-seeds.md`.

## 1. What the paper does and does not do

READ (abstract, p.387): "The purpose of this paper is to extend the study of the so called p-q resonant orbits of the
planar restricted three-body problem to the spatial case. ... If E, M and P denote the larger primary, the smaller one and the
infinitesimal body, respectively, then p and q are the number of revolutions that P gives around M and M around E,
respectively, between two consecutive close approaches. For fixed values of p and q and suitable initial conditions on a sphere
of radius mu^alpha around the smaller primary, we will derive expressions for the final position and velocity on this sphere
for the orbits under consideration."

Two points to take from that sentence.
1. The abstract sentence is garbled as printed: "p revolutions that P gives around M and M around E" is not what Definition 1.1 and
   Eq. 26 say. The equations say p = revolutions of P around E (the LARGE primary), q = revolutions of M around the origin
   (see section 3 below). Follow the equations, not the abstract. READ + COMPUTED (section 6).
2. The paper gives ONLY the outer (out) map and the necessary conditions for a p-q resonant orbit (Theorem 3.1). It does not
   give the inner solution, the in-map, or the matching; the conclusion (p.406) says "In a forthcoming paper this study will be
   used to get solutions of the RTBP close to periodic second species solutions" (that is the 2003 paper).
   So the 2002 paper alone does NOT produce a seed; it produces the admissible set of (C_J, direction) values and the outer map.

## 2. Model and conventions (section 1, pp.387-388)

- Spatial circular RTBP. Masses 1-mu (E, big) and mu (M, small); mu in [0,1]; primaries complete one inertial revolution in 2 pi;
  separation 1; G = 1. mu << 1 - mu is assumed throughout.
- Sidereal system OXYZ (inertial, origin at the barycentre; E at radius mu, M at radius 1-mu): Eq. 1,
  R'' = -(1-mu)(R - R_E)/R_1^3 - mu (R - R_M)/R_2^3.
- Synodic system Oxyz: E at r_E = (mu, 0, 0)^T, M at r_M = (mu - 1, 0, 0)^T (Moon on the left). R = G(t) r with G(t) a planar
  rotation by angle t (matrix printed in the 2003 paper, Eq. 14 there: [[cos t, -sin t, 0],[sin t, cos t, 0],[0,0,1]]).
- Eq. 2: r'' + A3 r' = grad Omega(r), A3 = [[0,-2,0],[2,0,0],[0,0,0]],
  Omega(x,y,z) = (1/2)(x^2 + y^2) + (1-mu)/r_1 + mu/r_2, r_1^2 = (x-mu)^2 + y^2 + z^2, r_2^2 = (x-mu+1)^2 + y^2 + z^2.
- Eq. 3: |r'|^2 = 2 Omega(x,y,z) - C_J. (Jacobi constant; NOT the half-sign convention.)
- Notation (p.389): subscript i = initial conditions (in both systems); subscript 1 / 2 on a position vector = measured from E / from M,
  r_1 = r - r_E, r_2 = r - r_M.

### Definition 1.1 (p.388), the p-q resonant orbit (verbatim, key parts)
"Let R(t) be the solution of (1) with initial conditions R(t_1) = R_i, R'(t_1) = R'_i and R_TB(t) the solution of R'' = -R/R^3 with the
same initial conditions. Denoting by B(M, mu^alpha) the sphere of center M and radius mu^alpha, 0 < alpha < 1, we will assume that
R_i is on B(M, mu^alpha) and that for a certain t_2 > t_1, R(t_2) belongs to B(M, mu^alpha) too, but that if t in (t_1, t_2) then
|R(t) - R_M(t)| > mu^alpha. Under these conditions, we say that the orbit is p-q resonant if

  t_2 - t_1 = 2 pi q + eps mu^alpha + O(mu^(2 alpha)) = 2 pi p tau + delta mu^alpha + O(mu^(2 alpha))     (4)

where p, q in N are relatively prime and 2 pi tau is the period of R_TB."
Consequences (READ): tau = q/p at leading order; so the two-body (heliocentric about E) orbit has period 2 pi q/p, p inertial Kepler
revolutions in the time the Moon makes q. p = spacecraft (P) revolutions about E, q = M revolutions. This is the labelling used in
all of Barrabes & Gomez and in the Casoliva Table 2 designations (checked, section 6).

## 3. The outer solution (section 2, pp.390-396)

Setting (p.390): at t = t_1 the orbit is on the sphere B(M, mu^alpha), moving away from M; t_2 is the first return to the sphere:
R_2(t_2) = mu^alpha, R_2(t) > mu^alpha on (t_1, t_2). Eq. 5: with a the angle between r_2i and r'_i,
cos a = <r_2i, r'_i>/(v_i mu^alpha) >= 0 (leaving the sphere). Eq. 6: v_i^2 = 3 - C_J + 2 mu^(1-alpha) + O(mu^(2 alpha)).
READ (p.389): "This equation implies that C_J <= 3. This is a logical restriction, since it is known that there are no zero velocity curves
in the planar case for C_J <= 3 - mu(1-mu). In the spatial case, there are zero velocity surfaces for values of C_J > -mu(1-mu), but for
values of the Jacobi constant between -mu(1-mu) and 3 - mu(1-mu) these surfaces do not intersect the z = 0 plane."

Eq. 7-9: q = (R, R') with q' = G(q,mu) + F(q,mu); G is the Kepler (E-only) field, F = (0, -mu (R - R_M)/R_2^3) the Moon term.

THEOREM 2.1 (p.391): "With the preceding notations and hypothesis and assuming, furthermore, that R = |R| cannot be arbitrarily small, we
have that q(t) = q_TB(t) + O(mu^(1-2 alpha)) for all t in [t_1, t_2]." (Gronwall, uses |F| <= mu/mu^(2 alpha).)
Remark (p.392): to make q_TB an approximation one needs alpha < 1/2.
LEMMA 2.2 (p.393): if |R_i - R_M(t_1)| = mu^alpha and cos a >= eps > 0 then there exist k >= eps' > 0 and t_(alpha/2) with R_2(t_(alpha/2)) = mu^(alpha/2)
and R_2(t) >= mu^alpha + k (t - t_1) on [t_1, t_(alpha/2)]. (Proof pp.394-396, a lower bound on R_2'' using Jacobi: K_alpha, p(k) parabola;
the proof requires v_i^2 cos^2 a - 4 K_alpha (mu^(alpha/2) - mu^alpha) >= 0, shown positive from (3 - C_J) cos^2 a - 8(1 + sqrt(3-C_J)) mu^(alpha/2) - ... > 0.)
THEOREM 2.3 (p.393): "If there exists t_(alpha/2) and t'_(alpha/2) such that t_1 < t_(alpha/2) < t'_(alpha/2) < t_2, R_2(t_(alpha/2)) = R_2(t'_(alpha/2)) = mu^(alpha/2),
R_2(t) > mu^(alpha/2) for t in (t_(alpha/2), t'_(alpha/2)), and the initial conditions verify cos a >= eps > 0, then q(t) = q_TB(t) + O(mu^(1-alpha)) for all t in [t_1, t_2]."
Choice of alpha (p.394), verbatim: "from now on we will impose that mu^(2 alpha) = O(mu^(1-alpha)) or, equivalently, alpha > 1/3. Thus, we have alpha in (1/3, 1/2)."
So: alpha in (1/3, 1/2); the orbit leaves with cos a >= eps > 0 (not tangentially).
Figure 1 (p.390): sample outer solution, synodic planar, mu = 1e-3, alpha = 0.4. The only numerical value of alpha printed anywhere in the paper is 0.4 (also the
q_max table, below).

## 4. Restrictions on the initial conditions (section 3, pp.396-403)

### 4.1 Spherical coordinates (Eq. 16, 17, 18; p.397), exactly as printed
  r_i = ( mu - 1 + mu^alpha cos(varphi) cos(theta),  mu^alpha cos(varphi) sin(theta),  mu^alpha sin(varphi) )^T,
  r'_i = v_i ( cos(phi) cos(psi),  cos(phi) sin(psi),  sin(phi) )^T.                                              (16)
Here varphi, theta are latitude/longitude of the position on the sphere about M (theta measured from +x); phi, psi are elevation/azimuth
of the velocity (psi measured from +x). (INFERRED reading of the symbols from the matrices; the paper does not name them.) Expansions:
  phi = phi_0 + Dphi mu^alpha + O(mu^(2 alpha)),   psi = psi_0 + Dpsi mu^alpha + O(mu^(2 alpha)).                  (17)
  cos a = cos a_0 + DC mu^alpha + O(mu^(2 alpha)),                                                               (18)
  cos a_0 = cos(varphi) cos(phi_0) cos(theta - psi_0) + sin(varphi) sin(phi_0),
  DC = Lambda_0 Dpsi + Dphi ( sin(varphi) cos(phi_0) - cos(varphi) sin(phi_0) cos(psi_0 - theta) ),
  Lambda_0 = cos(varphi) cos(phi_0) sin(theta - psi_0).
The subscript 0 denotes the zero-order term in mu^alpha. Planar case: varphi = phi = 0 (and so phi_0 = 0), then cos a_0 = cos(theta - psi_0).

### 4.2 Energy and angular momentum of the Kepler ellipse (p.397-399)
Kepler energy h = |R'_i|^2/2 - 1/|R_i|. Using |R_i|^2 = 1 - 2 mu^alpha cos(varphi) cos(theta) + O(mu^(2 alpha)) and |R'_i|^2 expanded:
  h = h_0 + Dh mu^alpha + O(mu^(1-alpha)),                                                                      (19)
  h_0 = 1 - C_J/2 - sqrt(3 - C_J) cos(phi_0) sin(psi_0),
  Dh = sqrt(3 - C_J) ( Dphi sin(phi_0) sin(psi_0) - Dpsi cos(phi_0) cos(psi_0) + cos(varphi) cos(phi_0) sin(psi_0 - theta) ) - 2 cos(varphi) cos(theta).
Ellipse needs h_0 < 0:
  sqrt(3 - C_J) < cos(phi_0) sin(psi_0) + sqrt( 1 + cos^2(phi_0) sin^2(psi_0) ),                                   (20)
"From this inequality we can deduce that C_J > -2 sqrt 2." Angular momentum c = R_i x R'_i not equal to zero (Kepler orbit not a collision):
  (3 - C_J) sin^2(phi_0) + ( 1 - sqrt(3 - C_J) cos(phi_0) sin(psi_0) )^2 not equal to 0,
  hence  phi_0 not equal to 0,  or  sin(psi_0) not equal to 1/sqrt(3 - C_J),  provided C_J < 3.                    (21)

### 4.3 Return to the sphere (Eq. 22, 23; p.399-400)
|R(t_2) - R_M(t_2)| = mu^alpha |w| + O(mu^(1-alpha)), with
  |w|^2 = 1 + (eps - delta)^2 + (3 - C_J) delta^2 + 2 (eps - delta) delta sqrt(3 - C_J) cos(phi) sin(psi)
          + 2 delta sqrt(3 - C_J) ( cos(varphi) cos(phi) cos(psi - theta) + sin(varphi) sin(psi) ) + 2 (eps - delta) cos(varphi) sin(theta).      (22)
The final term inside the second bracket is printed "sin varphi sin psi"; by analogy with cos a_0 in (18) one expects "sin varphi sin phi", so
INFERRED: a typo in the printed (22) (the zero-order condition (23) below is printed with cos a_0, which has sin(phi_0) there). It does not affect the planar case.
Requiring |w_0| = 1 gives the ellipse in the (delta, eps - delta) plane:
  (eps - delta)^2 + delta^2 (3 - C_J) + 2 (eps - delta) cos(varphi) sin(theta) + 2 delta sqrt(3 - C_J) cos(a_0)
     + 2 (eps - delta) delta sqrt(3 - C_J) cos(phi_0) sin(psi_0) = 0.                                              (23)
It is an ellipse unless 1 = cos^2(phi_0) sin^2(psi_0) (excluded). Figure 3 (p.401) draws these ellipses for p = 1, q = 2.

### 4.4 The two matching conditions between the ellipse and the Kepler period (Eq. 24-29; pp.400-402)
Period: tau^(2/3) = 1/(2|h|) (semi-major axis a = 1/(2|h|) with mu_E = 1), expanded
  1/(2|h|) = (1/(2|h_0|)) ( 1 + (Dh/|h_0|) mu^alpha + O(mu^(1-alpha)) ),                                          (24)
  tau = q/p + ((eps - delta)/(2 pi p)) mu^alpha + O(mu^(2 alpha)).                                                 (25)
Order zero gives the resonance relation (READ, Eq. 26):
  (q/p)^(2/3) = 1/(2|h_0|),   i.e.   C_J - 2 + 2 sqrt(3 - C_J) cos(phi_0) sin(psi_0) = (p/q)^(2/3).                (26)
Interval of C_J for given p/q (Eq. 27, p.400): from (26), (p/q)^(2/3) <= C_J - 2 + 2 sqrt(3-C_J) <= 2, hence "p/q <= 2 sqrt 2", and
  C_J1 = (p/q)^(2/3) - 2 sqrt( 2 - (p/q)^(2/3) )  <=  C_J  <=  (p/q)^(2/3) + 2 sqrt( 2 - (p/q)^(2/3) ) = C_J2.   (27)
The admissible C_J interval is [C_J1, C_J2] contained in (-2 sqrt 2, 3). Figure 3 caption (p.401): for p = 1, q = 2, "C_J can take any value on the interval
(-1.711013183, 2.970934233)". The endpoints themselves must be excluded (p.403, and 2003 p.149: cos^2 phi_0 sin^2 psi_0 = 1 there, not acceptable).
Collision-free Kepler orbit (Eq. 28): phi_0 not equal to 0, or C_J not equal to (p/q)^(2/3). Figure 4 (p.402): the admissible band as a function of p/q;
the dotted curve C_J = (p/q)^(2/3) is excluded when phi_0 = 0.

Size limit of the construction in terms of mu (p.402), verbatim: "Finally, |h_0| must not be too close to zero. For this, it will be enough that
(1/2) mu^(2 alpha)/h_0 = O(mu^(1-alpha)) or, equivalently, (p/q)^(2/3) > mu^(3 alpha - 1)  <=>  p > q mu^((3 alpha - 1) 3/2). This relation restricts
the range of values of p and q as a function of mu. For example, if we take alpha = 0.4, the maxim values for q when p = 1 are"

  | mu | 1e-3 | 1e-4 | 1e-5 | 1e-6 |
  | q_max | 7.943 | 15.84 | 31.6 | 63.09 |

COMPUTED check (p = 1, alpha = 0.4, so 3 alpha - 1 = 0.2): q_max = mu^(-0.3) = 7.943, 15.849, 31.623, 63.096. The printed 7.943 is exact; 15.84, 31.6 and 63.09
are the same numbers truncated (not rounded) to the printed digits. The table agrees with the printed inequality (exponent (3 alpha - 1) 3/2 = 0.3). READ+COMPUTED: no error.

Second condition (Eq. 29, p.402), the order-mu^alpha time condition, with Dh from (19):
  eps - delta = 6 pi q (q/p)^(2/3) ( sqrt(3 - C_J) ( cos(varphi) cos(phi_0) sin(psi_0 - theta) + Dphi sin(phi_0) sin(psi_0) - Dpsi cos(phi_0) cos(psi_0) ) - 2 cos(varphi) cos(theta) ).   (29)
COMPUTED derivation check: from (24), (25) and (26): (q/p)^(2/3) (Dh/|h_0|) = (1/3)(p/q)^(1/3) (eps - delta)/(pi p), and 1/|h_0| = 2 (q/p)^(2/3), giving
eps - delta = 6 pi q (q/p)^(2/3) Dh, which is (29) with Dh substituted. Consistent.
That (eps - delta) must lie between the roots k_1 < 0 < k_2 of the ellipse (23) (p.403):
  k_1,2 = ( cos a_0 cos phi_0 sin psi_0 - cos(varphi) sin theta  -/+  sqrt(Dk) ) / ( 1 - cos^2 phi_0 sin^2 psi_0 ),
  Dk = ( cos a_0 - 2 cos(varphi) cos phi_0 sin theta sin psi_0 ) cos a_0 + cos^2(varphi) sin^2 theta.
Writing lambda = 6 pi q (q/p)^(2/3), A, B, C as functions (A has the factor lambda; B = (cos a_0 cos phi_0 sin psi_0 - cos(varphi) sin theta)/(1 - cos^2 phi_0 sin^2 psi_0);
C = sqrt(Dk)/(1 - cos^2 phi_0 sin^2 psi_0)), the condition k_1 < eps - delta < k_2 becomes |A - B| <= C (Eq. 30).

THEOREM 3.1 (p.403), verbatim: "Let r(t) be a solution of (2) with initial conditions (r_i, r'_i) on the sphere B(M, mu^alpha). If it is a p-q resonant orbit, then the Jacobi
constant C_J in (C_J1, C_J2) and the variables varphi, theta, phi and psi, defined by (16) and (17), must verify C_J - 2 + 2 sqrt(3 - C_J) cos phi_0 sin psi_0 = (p/q)^(2/3),
phi_0 not equal to 0 or C_J not equal to (p/q)^(2/3), |A(varphi, theta, phi_0, psi_0, Dphi, Dpsi) - B(varphi, theta, phi_0, psi_0)| <= C(varphi, theta, phi_0, psi_0)."
This is a NECESSARY condition only ("If it is a p-q resonant orbit, then ..."): the paper does not prove existence; existence of orbits near the
second-species solutions is the subject of the 2003 paper and, for the periodic ones, of the continuation in Casoliva et al.

## 5. The outer map (section 4, pp.404-406), exact as printed
Result of the outer solution: R(t_2) = R_i + R'_i delta mu^alpha + O(mu^(1-alpha)), R'(t_2) = R'_i - R_i delta mu^alpha + O(mu^(1-alpha)) (sidereal, Eq. 31, using |R| = 1 + O(mu^alpha)).
Synodic form (Eq. 32), I the 3x3 identity, A3 as in Eq. 2:
  r_e  = ( I - ((eps - delta)/2) mu^alpha A3 ) r_i + delta mu^alpha r'_i,                                           (32a)
  r'_e = ( I - ((eps + delta)/2) mu^alpha A3 ) r'_i.                                                                (32b)
Spherical coordinates of the exit point (p.405): r_e = (mu - 1 + r_2e cos(varphi_e) cos(theta_e), r_2e cos(varphi_e) sin(theta_e), r_2e sin(varphi_e))^T,
r'_e = v_e (cos(phi_e) cos(psi_e), cos(phi_e) sin(psi_e), sin(phi_e))^T. The exit distance is r_2e = mu^alpha ( 1 + mu^alpha Dr_e + O(mu^(1-alpha)) ) with
  Dr_e = delta sqrt(3 - C_J) [ Lambda_0 (eps - delta) + Dpsi ( Lambda_0 + (eps - delta) cos phi_0 cos psi_0 )
          + Dphi ( sin(varphi) cos phi_0 - cos(varphi) sin phi_0 cos(psi_0 - theta) - (eps - delta) sin phi_0 sin psi_0 ) ] - (eps - delta)^2 cos(varphi) cos(theta).      (33)
Direction cosines of the exit point (Eq. 34a-c, pp.405-406), with s = sqrt(3 - C_J):
  cos(varphi_e) cos(theta_e) = cos(varphi) cos(theta) + delta s cos(phi_0) cos(psi_0)
        - mu^alpha [ delta s ( Dphi sin phi_0 cos psi_0 + Dpsi sin psi_0 cos phi_0 ) - (eps - delta) cos(varphi) sin(theta)
                    + Dr_e ( cos(varphi) cos(theta) + delta s cos phi_0 cos psi_0 ) ] + O(mu^(1-alpha)),            (34a)
  cos(varphi_e) sin(theta_e) = cos(varphi) sin(theta) + delta s cos phi_0 sin psi_0 + (eps - delta)
        + mu^alpha [ delta s ( Dpsi cos phi_0 cos psi_0 - Dphi sin psi_0 sin phi_0 ) - (eps - delta) cos(varphi) cos(theta)
                    - Dr_e ( cos(varphi) sin(theta) + (eps - delta) + delta s cos phi_0 sin psi_0 ) ] + O(mu^(1-alpha)),       (34b)
  sin(varphi_e) = sin(varphi) + delta s sin phi_0 + mu^alpha [ delta s Dphi cos phi_0 - Dr_e ( sin(varphi) + delta s sin phi_0 ) ] + O(mu^(1-alpha)).    (34c)
Velocity (Eq. 35): v_e = v_i + O(mu^(2 alpha)); if cos(varphi) not equal to 0,
  phi_e = phi + O(mu^(2 alpha)),   psi_e = psi_0 + ( Dpsi - (eps + delta) ) mu^alpha + O(mu^(2 alpha)).                (35a, b)
(The symbol "Delta r_e" in (33) is called "Delta r_e"; the 2003 paper prints the same expression with the Dphi bracket written
sin(varphi) cos phi_0 - sin phi_0 ( cos(varphi) cos(psi_0 - theta) + (eps - delta) sin psi_0 ), which is algebraically the same. COMPUTED by expanding.)

## 6. Checks from printed numbers (all COMPUTED, 2026-10-04, Python, double precision)
1. Eq. 27 for p = 1, q = 2: C_J1 = -1.711013183, C_J2 = 2.970934233; printed on p.401: (-1.711013183, 2.970934233). Agree to all printed digits.
   Figure 3 value C_J = 0.6299605255 equals (1/2)^(2/3) = 0.629960525 (the excluded value for phi_0 = 0, Eq. 28).
2. The q_max table (section 4.4): agrees with p > q mu^((3 alpha - 1) 3/2), values truncated rather than rounded.
3. Eq. 29 follows from Eq. 24, 25, 26 (section 4.4).
4. Planar bound p/q <= 2 sqrt 2 = 2.8284271 (the same constant that Casoliva et al. use for the "forbidden resonances").
5. Labelling: for the Casoliva 2008 Table 2 rows (p-q designation, mu = 1e-6) the semi-major axis computed from the printed state is (q/p)^(2/3)
   to within 6e-3 relative, and T is within 3.1e-2 of 2 pi q, for all ten rows (see the Casoliva 2008 digest). So p = spacecraft revolutions,
   q = Moon revolutions, in this paper's Eq. 26 and in those designations.
6. Cross-paper check, Bolotin 2006 (p.238, Remark 3): its Jacobi integral h = H - G is for h in (-3/2, sqrt 2). With C = -2h (the 2002 convention for C_J relates to h
   by C_J = -2h; INFERRED from Casoliva Eq. 3 and 5), h in (-3/2, sqrt 2) maps to C_J in (-2.828, 3.0), the same range as
   (-2 sqrt 2, 3) here (Eq. 27 with its limits). COMPUTED endpoints: -2 sqrt 2 = -2.8284271, 3.

## 7. Failure modes and limits (what the paper says, and what follows)
- alpha in (1/3, 1/2) only; the error is O(mu^(1-alpha)) and no constant is given. "mu_1" below which the construction holds is never quantified.
- cos a >= eps > 0: orbits that leave the sphere tangentially (cos a -> 0) are outside the theory; the Lemma 2.2 constants degenerate.
- (3 - C_J) cos^2 a must beat 8 (1 + sqrt(3 - C_J)) mu^(alpha/2) (Lemma 2.2 proof, p.396): a lower bound on the leaving angle that grows with mu^(alpha/2); at
  mu = 1e-6 and alpha = 0.4, mu^(alpha/2) = 1e-6^0.2 = 0.063 (COMPUTED), so the constants are not small unless cos a is of order one.
- C_J strictly inside (C_J1, C_J2); the endpoints are cos^2 phi_0 sin^2 psi_0 = 1.
- The theorem is for the OUTER solution only; a seed also needs the inner solution and the matching (2003).
- The sign convention hazard: the abstract's "p revolutions around M, M around E" is inverted relative to Eq. 26 (section 1 above).
- Not in the paper: any numerical example of a corrected orbit, any table of seeds, any continuation in mu. Those are in Casoliva et al.

## 8. How this applies to the project
- The in/out-map technique is the generator behind `#899` (P2 in `docs/notes/2026-10-04-897-technique-synthesis-papers-to-problems.md`). Everything needed
  from this paper for the PLANAR seed is: Eq. 16 to 18 (coordinates), Eq. 26 and 27 (C_J interval and psi), Eq. 32 (outer map, only needed to verify), and the list of
  restrictions in section 4. The planar specialisation (varphi = phi_0 = 0) of the matching is in the 2003 paper, Eq. 44 to 46.
- Positive controls that this paper alone gives: Eq. 27 for p = 1, q = 2 reproduces (-1.711013183, 2.970934233); the q_max table; the admissibility bound p/q <= 2 sqrt 2.
  A programmer can implement these in a minute and they test the labelling convention (p and q not swapped) before anything else is built.
- A mismatch of labelling is the most likely silent error: with p and q swapped the interval for (1, 2) becomes the interval for (2, 1), (0.3027, 2.8721).

## 9. Proposed corpus-index row
| barrabes-gomez-2002-spatial-pq-resonant-orbits-rtbp-cmda-84-387-doi-10.1023-A1021137127909.pdf | 2026-10-04-digest-barrabes-gomez-2002-spatial-pq-resonant-orbits.md | Spatial p-q resonant orbits of the CR3BP at small mu: necessary conditions on C_J (interval Eq. 27), psi (Eq. 26), the period/ellipse matching (Eq. 23, 29, 30) and the explicit OUTER map (Eq. 32-35) on a sphere of radius mu^alpha, alpha in (1/3, 1/2); no inner map or matching (that is the 2003 paper); text layer (7921 words); full digest with checks. | DIGESTED (#920) |
