# Digest: Font, Nunes and Simo 2002, consecutive quasi-collisions in the planar circular RTBP

Date 2026-10-04. Purpose: establish what is proved and computed for ONE small secondary about orbits that make
repeated close encounters with it, so that the #890 result (a periodic orbit of the planar concentric circular
restricted four-body problem, Uranus plus Titania plus Oberon, one close flyby of each moon per cycle; see
`docs/notes/2026-10-04-890-titania-oberon-candidate.md` and `docs/notes/2026-10-04-890-literature-check.md`) is
worded correctly. Companion digests: `docs/notes/2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`
and `docs/notes/2026-10-04-digest-font-nunes-simo-2009-second-species-numerical-study.md`.

Citation: J. Font, A. Nunes and C. Simo, "Consecutive quasi-collisions in the planar circular RTBP",
Nonlinearity 15:115-142 (2002), DOI 10.1088/0951-7715/15/1/306. 28 journal pages (the PDF has an IOP cover
sheet, so PDF page n is journal page n+113). Received 6 November 2000, published 12 December 2001.
Filed in the private paper corpus as
font-nunes-simo-2002-consecutive-quasi-collisions-planar-circular-RTBP-nonlinearity-15-115-doi-10.1088-0951-7715-15-1-306.pdf.

Evidence labels: READ means read from the printed page in this session (all pages were read as page images; the
long digit strings of the initial conditions were re-read at 220 dpi). INFERRED means my reading of what the
print implies. Page numbers are journal pages. Formulas are transcribed from images; anything marked [unclear]
could not be read with confidence.

## 0. Summary (READ)

Abstract (p. 115): "In this paper, we consider the planar circular restricted three-body problem and, in
particular, the existence of orbits which undergo consecutive close encounters with the small primary. The
number of revolutions of the small bodies around the larger one between successive encounters can be chosen to
be two arbitrary sequences of natural numbers, with constraints depending on the Jacobi constant. We prove
that such orbits exist as a consequence of the fact that, when the mass parameter mu is small, the first
return map defined on a region of phase space whose projection is a circle around the small primary is a
'horseshoe' map. The proof is constructive, in the sense that it is based on the computation of an
approximate expression for this return map. When mu is small, the approximate return map contains the
essential information about the dynamics from the quantitative as well as qualitative point of view. Using
this information, we have been able to carry out a numerical study of this problem for mu up to 10^-3."

The method is a perturbative construction (resonant strips, OUT-map, IN-map, return map), not a variational
one. The theorem (Theorem 5.1) is stated with an explicit smallness condition mu(N) << N^-10 (eq. (57)), so
unlike Bolotin and MacKay it carries a stated scaling of the admissible mass, though still no numerical value.

## 1. Setting (pp. 115-117, READ)

Model (p. 115 to 116): "Consider the planar circular restricted three-body problem, choosing, as usual,
appropriate mass, length and time units so that the masses of the two primaries are m_1 = 1 - mu, m_2 = mu,
mu in (0, 1/2), the angular velocity of their motion around the fixed centre of mass is 1 and their positions
in the plane are given in synodic coordinates (x, y) by (mu, 0) and (mu - 1, 0), respectively." Equations of
motion (1):

  xddot - 2 ydot = dOmega/dx = x - (1 - mu)/r1^3 (x - mu) - mu/r2^3 (x + 1 - mu),
  yddot + 2 xdot = dOmega/dy = y - (1 - mu)/r1^3 y - mu/r2^3 y,

with Omega = (1 - mu)/r1 + mu/r2 + (x^2 + y^2)/2 + mu(1 - mu)/2, r1^2 = (x - mu)^2 + y^2,
r2^2 = (x - mu + 1)^2 + y^2. Jacobi function (2): C_J = 2 Omega - (xdot^2 + ydot^2). "[The equations] are
singular at the points (x, y) = (mu, 0) and (x, y) = (mu - 1, 0), which correspond to collision of the massless
body Z with the big primary E and with the small primary M, respectively." "We recall that C_J = 2(M_sid -
E_sid), twice the difference between the angular momentum and the energy of Z with respect to the centre of
mass." Z is the massless body, E the big primary, M the small primary. The frame is the synodic (rotating)
frame with M at (mu - 1, 0), to the left of the origin.

Regularisation (p. 126 to 127): the IN-map uses rescaled Hill-type variables a_1 = mu^(-1/3)(x + 1 - mu),
a_2 = mu^(-1/3) y, b_1 = mu^(-1/3) xdot - a_2, b_2 = mu^(-1/3) ydot + a_1 (eq. (35)), a time-dependent rotation
(38), then Levi-Civita (41): "a~_1 + i a~_2 = z~^2, b~_1 + i b~_2 = w~/conj(z~), dt/d(rho) = 2|z~|^2", then a scaling
(44) with epsilon = mu^(alpha - 1/3) that maps "the disk of radius mu^alpha around (-1 + mu, 0) in synodic
coordinates on the unit ball in the complex z plane". Numerically (p. 133): "All the orbits are obtained by
numerical integration of the equations of the RTBP with regularization of collision with the small primary,
using the Runge-Kutta-Fehlberg algorithm of orders 7-8."

Definition 1 (p. 116), quoted: "We shall say that a solution gamma(t), t in [t_1, t_2], of equations (1) is
p-q resonant, p, q in N and coprime, (p, q) = 1, if both gamma(t_1) and gamma(t_2) are on a circle C of radius
mu^alpha centred at M, with alpha in (1/3, 1/2) a real parameter whose value will be specified later, if
gamma(t) is outside the disk D of radius mu^alpha centred at M for all t in (t_1, t_2), and if moreover the time
interval t_2 - t_1 between these successive approaches to M is such that
(1/2pi)(t_1 - t_2) = q + O(mu^alpha) = p tau + O(mu^alpha) holds, where 2 pi tau is the period of the orbit of
Z as a two-body problem with respect to E with the same initial conditions of gamma at t_1." (The printed
left-hand side is "(t_1 - t_2)", presumably a sign slip for t_2 - t_1; the same relation is written
t_2 - t_1 = 2 pi q + O(mu^(alpha)) on p. 123.) So q is the number of full turns of the synodic frame and p the
number of Kepler revolutions about E between two consecutive approaches.

Definition 2 (p. 116): "A solution of equations (1) will be called a consecutive quasi-collision orbit if it
can be obtained by concatenation of resonant orbits and, for each such consecutive couple of such orbits
having, say, related time intervals [t_1, t_2] and [t_3, t_4], the solution remains inside D for
t in (t_2, t_3)."

Encounter size: the encounter circle has radius mu^alpha with alpha in (1/3, 1/2) (eq. (5), p. 118); the
numerics use alpha = 2/5 (p. 133: "We take the parameter alpha = 2/5 since the errors involved in the analytic
computations of the OUT- and IN-maps match for this value of alpha, and global error estimate for the return
map is then of order mu^(1/5)"). Since alpha > 1/3 the circle is smaller than the Hill radius scale mu^(1/3)
(INFERRED: the Hill radius is (mu/3)^(1/3); alpha = 2/5 gives radius 0.0251 at mu = 1e-4 and 0.063 at mu = 1e-3,
the figure captions show the circle of that radius).

Why alpha is pinned at 2/5 (p. 129): "The expansions in (50) are in powers of mu^alpha and mu^(1-alpha). We note
that the lowest-order congruence between these exponents, taking into account the restriction (5), is
obtained for 2(1 - alpha) = 3 alpha and then alpha = 2/5."

Jacobi constant range (p. 118, 122): Delta_J = sqrt(3 - C_J) > 0 (Remark 2.1), "3 - C_J > 0 is bounded away
from zero", i.e. orbits that can reach the small primary's neighbourhood from the large primary's side with C_J
below 3. Admissibility of a pair (p, q) at a given C_J, eq. (21): -2 Delta_J <= 2 - C_J + (p/q)^(2/3) <=
2 Delta_J, and p/q < 2 sqrt 2; "Given C_J, we shall say that a pair (p, q) is admissible if (21) is verified
strictly (see figure 1). Note that this excludes the orbits which are close to orbits of the two-body problem
with perpendicular crossings of the x-axis." Figure 1 (p. 123) plots the admissible region in the plane
(C_J horizontal, (p/q)^(2/3) vertical); the plotted C_J axis runs from -3 to 4 and the dashed admissible region
spans roughly C_J from -3 to 3 (reading from the figure; INFERRED for the exact extent).

Avoiding intermediate close encounters (pp. 116-117): "At this point we show that these encounters can be
avoided by choosing C_J outside a finite union of small intervals, provided mu is small enough. The number of
these intervals depends on p and q." The argument is the Kepler collision condition (3), (4) at the second
intersection point B of the orbits of M and Z. Eq. (3):
2 pi (j_M + delta_M) +- 2 arctan( sqrt(1 - e^2) sin E / (cos E - e) ) = (q/p)(2 pi (j_Z + delta_Z) - 2(E - e sin E)),
with constraint (4) 1 = a(1 - e cos E). "Notice that this argument is essentially contained in figure 4 of the
pioneering work of Henon [2]."

## 2. The OUT-map and IN-map (pp. 117-130, READ)

Proof outline (p. 117): "(a) In section 2 we derive an approximate expression for the resonant strips, i.e.
for the set of initial conditions on the circle C of radius mu^alpha around M which correspond to p-q
resonant orbits. (b) In section 3 we obtain an approximate expression for the OUT-map, i.e. the map sending a
point on a p-q resonant strip to the position and velocity of the orbit of that point at the time of its first
intersection with the circle C. (c) In section 4 we obtain an approximate expression for the IN-map, i.e. the
map sending the initial conditions of an orbit entering the disk D to the position and velocity of that orbit
when it crosses again the circle C leaving the disk D. (d) Finally, in section 5 we select the optimal alpha to
minimize the errors made in the previous approximations and show that the composition of the OUT- and IN-maps
defined on the resonant strips is a 'horseshoe' map."

Coordinates on C (p. 121): the triad (phi, v_s, psi_s), "phi is the angular coordinate on the circle defined
above, v_s is the modulus of the synodic velocity and psi_s is the direction of the velocity of the particle in
synodic coordinates". The synodic speed on C (6): v_s^2 = -C_J + 3 + 2 mu^(1-alpha) + O(mu^(2 alpha)).

Lemma 2.1 (p. 120): approximation of a resonant orbit between encounters by the two-body orbit about E, with
error O(mu^(1-alpha)) provided 3 - C_J > 0, |cos a_i| bounded away from 0, and the orbit crosses the circle of
radius mu^(alpha/2) only twice (so no intermediate close encounter).

Resonant strips (Definition 3, p. 123, eqs. (22), (23)): "For a given value of the Jacobi constant C_J, a p-q
resonant strip is the set of values (phi, psi_s) in the annulus {(phi, psi_s) in S^1 x S^1 : |phi - psi_s| <
pi/2} such that psi_s = psi~_s - mu^alpha (B(phi) + xi C(phi))/K", with
v_0^2 = 2 - (p/q)^(2/3), e = (2 - v_0^2)^(5/2)/(6 pi p),
sin psi~_s = (4 - C_J - v_0^2)/(2 Delta_J), cos psi~_s = sigma sqrt(4 v_0^2 - (2 - C_J - v_0^2)^2)/(2 Delta_J),
B(phi) = (1/(Delta_J - sin psi~_s)) (2 cos phi + (Delta_J - e/cos psi~_s) sin(phi - psi~_s)),
C(phi) = e/((Delta_J - sin psi~_s) cos psi~_s), K = 2 Delta_J^2 cos psi~_s/(2 - C_J + v_0^2),
sigma = sign(cos psi~_s), xi in [-xi_m, xi_m] subset (-1, 1). "The condition xi_m < 1 ensures that the resonant
orbits enter the disk D away from tangency." The strips have width of order mu^alpha in the (phi, psi) torus
(p. 156 of the companion paper repeats this).

OUT-map (p. 124 to 126), result Lemma 3.1 (34): phi_e = psi~_s + pi - arcsin(xi) + O(mu^(1-2 alpha)),
psi_e = psi~_s - mu^alpha (k_1 cos phi + k_2 sin phi + k_3 xi + k_4 sqrt(1 - xi^2)) + O(mu^(1-alpha)), with
v_e given by (33): v_e = Delta_J + (2 mu^(1-alpha) + 3 cos^2(phi_e) mu^(2 alpha))/(2 Delta_J) + O(mu^gamma),
gamma = min{1, 3 alpha}, and
k_1 = 2e(4 - C_J - v_0^2)/(4 v_0^2 - (2 - C_J - v_0^2)^2), k_2 = -2 sigma e/sqrt(4 v_0^2 - (2 - C_J - v_0^2)^2),
k_3 = 4 e Delta_J/(4 v_0^2 - (2 - C_J - v_0^2)^2) + sigma (2 v_0^2 - 2)/(Delta_J sqrt(4 v_0^2 - (2 - C_J - v_0^2)^2)),
k_4 = -2/Delta_J. (The placement of the squares in the k_1, k_3 denominators is as read from the image;
treat as [unclear] if used for computation.)

IN-map (pp. 126-130). Away from collision, Lemma 4.1 (p. 129), valid when |phi_e - psi_e - pi| >> mu^alpha:
  phi_s = 2 psi_e - phi_e + pi + 4 mu^alpha cos(phi_e - psi_e)/Delta_J + mu^(1-alpha) 2 cos(phi_e - psi_e)/(Delta_J^2 sin(phi_e - psi_e)) + O(mu^(2 alpha)),
  psi_s = psi_e + 4 mu^alpha cos(phi_e - psi_e)/Delta_J + mu^(1-alpha) 2 cos(phi_e - psi_e)/(Delta_J^2 sin(phi_e - psi_e)) + O(mu^(2 alpha)),
  v_s = sqrt(3 - C_J - 3 mu) + (2 mu^(1-alpha) + 3 cos^2(phi_s) mu^(2 alpha))/(2 Delta_J) + O(mu^(3 alpha)).
So, away from the collision direction, the exit velocity direction psi_s differs from the entry direction psi_e by
O(mu^alpha) (and the mu^(1-alpha) term grows as sin(phi_e - psi_e) -> 0, i.e. toward the collision direction).
Near collision, Lemma 4.2 (p. 129): "Assume that phi_e - psi_e = pi - Delta_J^(-1) mu^alpha - eta Delta_J^(-2)
mu^(1-alpha). Then ... phi_s = phi_e + Delta phi - ..., psi_s = phi_s - mu^alpha/Delta_J + O(mu^(4 alpha - 1)),
where cos Delta phi = (1 - eta^2)/(1 + eta^2), sin Delta phi = 2 eta/(1 + eta^2)". The parameter eta sweeps the
deflection through a full turn; the transition region is of order mu^(2 alpha - 1) in eta.

Time inside D (p. 128): crossing time T = T_0 + O(mu^(beta + alpha) log mu), T_0 = -2 mu^alpha Delta_J^(-1)
cos(phi - psi), with beta = 2 alpha (resp. alpha) away from (resp. close to) collision.

## 3. Main results (pp. 130-133, READ)

Return map (Lemma 5.1, p. 130), far from collision (|xi - xi_0| >> mu^alpha):
phi_s = psi~_s + arcsin(xi) + O(mu^(1-2 alpha)), psi_s = psi~_s + mu^alpha g(phi, xi) + O(mu^(1-alpha)),
g = -k_1 cos phi - k_2 sin phi - k_3 xi - k_4 sqrt(1 - xi^2). "The behaviour of this return map is essentially the
following: it maps each resonant strip into another strip of width O(mu^alpha), close to the resonant strip,
after contracting it by a factor of order mu^(1-2 alpha) along the horizontal direction (i.e. the direction of the
phi axis) and expanding it by a factor of order mu^(-alpha) along the vertical direction (the direction of the
psi-axis)." Eq. (53): |dF/dphi| = O(mu^(1-2 alpha)), |(dF/dpsi)^(-1)| = O(mu^alpha). "Notice that this hyperbolic
behaviour is entirely due to the OUT-map". Lemma 5.2 (p. 130): the map near collision.

At alpha = 2/5 (p. 131): "For this value of alpha, the error in the return map is of order mu^(1/5), and all the
expansions are in integer powers of mu^(1/5)." The stretching factor is mu^(-2/5), the contraction mu^(1/5).

Theorem 5.1 (p. 131), quoted exactly: "For any integer N, there exists mu(N) such that any admissible sequence
(p_n, q_n) with q_n < N, n in Z, is realized by a consecutive quasi-collision orbit of the RTBP with
0 < mu < mu(N)."

Proof outline (pp. 131-133): annulus A = {|phi - psi| < c mu^beta}, beta = 1 - 2 alpha; countable set of
(p, q, sigma) strips, finitely many with q < N; "the distance between consecutive strips has a lower bound of
order 1/(N^2 - N)"; Conley-Moser conditions on horizontal and vertical strips H_i, V_i give the Smale-Moser
horseshoe and a conjugacy to a subshift on the alphabet A = {1, ..., M} with transition matrix M - Id
(excluding repeating the same strip). The two constraints on N (p. 133): the period of the two-body orbit must
be small enough that the error bounds hold: "N mu must be at most of order mu^(1-alpha) i.e. N must be much
smaller than 1/mu^alpha. For alpha = 2/5 this is assured by N << mu^(-1/10)", and the strips must not overlap
(distance between strips 1/(N^2 - N) exceeds mu^(1-2 alpha)). Result (57): "mu(N) << N^(-10)." "Note however
that the first constraint, mu << N^(-2.5), is much less restrictive."

Periodic orbits: "In particular, our result implies the existence of infinitely many consecutive quasi-collision
periodic orbits, corresponding to periodic sequences (p_n, q_n)_(n in Z). These are what Poincare called periodic
orbits of the second species (POSS) [7]. The existence of symmetric (with respect to the line joining the two
primaries in synodic coordinates) POSS was proved by Perko [6], and by Henrard [3] and Marco and Niederman [4]
by a different method. Most of the POSS that we find, however, are asymmetric. The existence of POSS in a more
general setting, which includes the RTBP, was proved by Bolotin and MacKay [1] through variational methods."
(p. 117 to 118.)

Instability: the printed stability parameters of the computed orbits are 2.338645E+7, 1.132360E+6 and
-4.228471E+9 at mu = 1e-4 (Section 4 below; the paper uses the phrase "the orbit is unstable and the stability
parameter is near ..." and does not define the stability parameter in the text I read). The analytic scaling
statement is the mu^(-alpha) vertical stretching per encounter (eq. (53)) against mu^(1-2 alpha) contraction.

### What the paper itself says about the size of angle changes (the Bolotin and MacKay remark)

The phrase "small angle changes" is Bolotin and MacKay's (2006, p. 434: "limited to orbits with small angle
changes at collisions"); Font, Nunes and Simo do not use it. What they print (READ):
- The horseshoe is built on the annulus |phi - psi| < c mu^beta, beta = 1 - 2 alpha, of encounter geometries
  in which the orbit passes the small primary on a non-collision geometry; for these, Lemma 4.1 (eq. (50)) gives
  a velocity-direction change of order mu^alpha (psi_s - psi_e = 4 mu^alpha cos(phi_e - psi_e)/Delta_J + ...).
- Orbits aimed close to collision are covered by Lemma 4.2 and "are mapped by the IN-map on ejections with
  arbitrary exit angle", but only in "a small strip of width of order mu" (p. 131), and the horseshoe argument
  uses them only through the gap in the annulus, of length of order mu^(1 - 2 alpha) (p. 132), not as symbols.
- The conclusion (p. 141) says the result is for the case "when the quasi-collisions take place, essentially,
  around the same location in sidereal coordinates", i.e. the same circle C is revisited.
INFERRED: this is the content behind "small angle changes": the symbols of the shift are resonant arcs that
leave and re-enter the same small circle after p Kepler revolutions, with the encounter bending the velocity by
O(mu^alpha) at the distances used; the large turns of a real flyby (tens of degrees) arise only for encounter
distances of order mu^(1/3) or less or near the collision direction, which the paper treats only as a thin
transition set. The companion paper (2009) adds the S-arcs (collision-to-collision arcs) and recovers
arbitrary exit angles after a near-collision encounter; see its digest.

## 4. Numerical results (pp. 133-141, READ)

Parameters: alpha = 2/5 throughout. Jacobi constants used: C_J = 2.8 (most figures) and C_J = 1.8 (Figure 7).
Mass parameters: mu = 1e-6, 1e-4, 1e-3 (Figures 5, 6, 7), and the continuation 1e-4 to 1.5e-3 for the (7,4,+1)
strip (Figures 11 to 14). The largest value printed is mu = 1.5e-3 (Figure 14(c)). The Earth-Moon value
(about 0.0122) is not reached. No Jacobi constant above 2.8 or any value near 3 is used.

Numerical procedure (p. 134, items (a) to (h)): fix mu, C_J, sigma; find on the circle C an initial condition
(phi, psi) whose orbit, after about 2 pi q time units, is tangent to the circle (a boundary point of the strip);
continue in phi along the strip boundary; the strip is bounded by the two tangency lines |phi - psi| = pi/2 and
the two tangency curves. "Due to the high stretching rate of the map, it is necessary to take a small step size
on the lines that correspond to lines b and d of figure 3(a)." The strips and images (Figures 4, 5, 6) are
integrations of the regularised equations.

Admissible pairs at C_J = 2.8 (p. 134 to 136): bound q <= 7 "because of the restriction q << mu^(-alpha)
(see (57)) and because we want to consider values of mu as high as 1e-3. Note that for this value of mu, the
bound N = 7 is not even small compared to mu^(-2/5)." Pairs removed because of intermediate close encounters,
found "by an algorithm based on Henon's work" at mu = 1e-3 with Delta t < (1e-3)^(1/5): eq. (58)
{(5,4), (4,5), (6,5), (7,5), (8,5), (11,5), (5,6), (11,6), (13,6), (5,7), (10,7), (11,7), (12,7), (13,7), (15,7)}.

Figures 4 to 7 (pp. 134-137): resonant strips and images for (2,1,-1) at mu = 1e-4, C_J = 2.8 (Fig. 4); for
(5,3,-1) at mu = 1e-6, 1e-4, 1e-3 (Fig. 5); the full set for q <= 7 at the three mu values with strips with q <= 3
labelled (Fig. 6); and for C_J = 1.8, (p,q,sigma) = (3,2,1) at mu = 1e-6, 1e-4, 1e-3, compared with the analytic
strip (22) and image (34), (49) (Fig. 7). Text (p. 137): "numerical curves are well approximated by the
analytical curves, even for mu = 1e-3. ... the separation has been checked to be O(mu^gamma), gamma >= 2/3, and
hence of the order of the error terms of the analytical expressions." No numerical tables of strip positions
are printed; the content is in the figures.

Periodic orbits (pp. 137-139), all at mu = 1e-4, C_J = 2.8. "We use extended precision to compute the periodic
orbits and the stability parameter." The printed data (digits as read; (phi, psi) are the angle on the circle
of radius mu^alpha around M and the synodic velocity direction; "the initial conditions are (phi, psi)"):

| Figure | symbolic sequence (p_n, q_n, sigma_n) | phi | psi | stability parameter | symmetry |
|---|---|---|---|---|---|
| 8 | (4,3,1) for n odd, (2,1,1) for n even | 0.4840458051695259 | 0.3657582929540782 | 2.338645E+7 | asymmetric (not symmetric about the synodic x axis) |
| 9 | (2,1,1) n odd, (2,1,-1) n even | 1.060031115468069 | 0.9672594502381595 | 1.132360E+6 | symmetric about the x axis |
| 10 | period three: (2,1,1), (2,1,-1), (3,2,1) | 0.91622091785178612688 | 0.94492534193201301129 | -4.228471E+9 | asymmetric |

Captions give the "near" value of the stability parameter, so the figures are as printed, not more precise.
No period (orbital period) values are printed in this paper for these orbits. (Compare the 2009 companion,
which does print periods.) Figure 8 caption: arcs 1 to 4, "the dotted line, from 1 to 2, is the part of the
orbit associated with the (p,q,sigma) = (4,3,1) resonance, and the full line, from 3 to 4, ... (2,1,1)".

Symmetric orbits recipe (p. 139): "it is possible to obtain symmetric periodic orbits, for instance, by taking a
sequence of the form (p,q,sigma) if n odd and (p,q,-sigma) if n even."

Loss of strips with growing mu (pp. 139-141):
- "as the parameter mu increases from 1e-4 to 1e-3, several resonant strips have been lost. This may be due to
  two different causes. One of them is illustrated by the disappearance of the (9,7,+-1) strips. Although
  (9,7) is not in the set (58), if we define a close encounter in a slightly less restrictive way, say by
  requiring Delta t to be smaller than 1.6 x (1e-3)^(1/5), we find that the corresponding orbit has two close
  encounters between resonances. The resulting accumulated perturbation destroys the resonant strip."
- "The other mechanism is more subtle. In this case, the resonant orbits still exist, but the locus of the
  corresponding initial conditions is strongly deformed. It ceases to be close to a horizontal strip, the
  analytical theory cannot be applied (mu is too large) and it can no longer be found by our standard numerical
  procedure."
- The (7,4,+1) strip followed from mu = 1e-4 to 1.5e-3 (Figures 11-14), together with four sets S_a^+, S_b^+,
  S_a^-, S_b^- of orbits returning after "approximately 2 pi q time units" (return times differ by about 5 percent
  from 2 pi q): "S_b^- and S_b^+ merge for mu around 5.25e-4 to form a unique connected component S. This set S
  will in turn merge with the (7,4,+1) strip for mu around 6.74e-4, to form the set S'. As mu increases, the phi
  range of the set S' gets smaller and S' merges with the set S_a^+, producing a connected set S'', for mu
  around 1.4e-3, while (7,4,-1) and S_a^- get closer (see figure 14(b)). For mu still larger, these three sets
  merge to form a unique connected set S''' (see figure 14(c))." Figure 14(c) is at mu = 1.5e-3. Figures 12 and
  13 caption values: mu = 5e-4, 5.25e-4, 5.5e-4; Figure 14: 6.74e-4, 1.4e-3, 1.5e-3.

## 5. Conclusions and statements on other settings (pp. 141-142, READ)

Quoted in full from p. 141 to 142: "The existence of orbits having consecutive passages very close to the
smallest of the primaries (quasi-collisions) in the RTBP has been shown. Furthermore, the passages can be
selected with some 'randomness' due to the freedom in selecting the involved resonances between the mean
motions. Concerning this problem existence results are available in the literature, but we have taken a direct
and more constructive approach which allows us to make accurate predictions of the initial conditions. This is
illustrated with some highly unstable periodic orbits, either symmetric or not.

For conciseness, we have considered the case when the quasi-collisions take place, essentially, around the same
location in sidereal coordinates. The same ideas can be used when quasi-collisions taking place near different
locations along the orbit of the smallest primary, occur.

Beyond theoretical considerations, these orbits are relevant when it is interesting to perform successive
fly-byes of the small primary by a spacecraft. This can help put the spacecraft on a desired nominal orbit
around the largest primary at low cost."

Elliptic problem (p. 118): "We remark that although we have considered the equations of motion for the circular
RTBP, essentially the same steps can be carried out in the restricted elliptic case." Not carried out in the
paper (READ: no elliptic computation or theorem appears).

More than one small body: not mentioned anywhere in the paper (READ, by absence: the model is the three-body
problem with ONE small primary; the "different locations along the orbit of the smallest primary" sentence
concerns the same single small primary at different phases of its own orbit, not a second secondary).
Time-dependent problem beyond the elliptic remark: not mentioned. Spacecraft applications: only the closing
sentence above, with no mission analysis.

## 6. What this does and does not cover for #890

(a) Established for ONE small secondary. Analytically: Theorem 5.1, for fixed Jacobi constant C_J < 3 outside a
finite union of small intervals, a stated smallness mu(N) << N^(-10) for symbol sequences with q_n < N; the
encounter distance is of order mu^alpha with alpha in (1/3, 1/2), alpha = 2/5 in the numerics; the unstable
direction is stretched by about mu^(-alpha) per encounter. Numerically: the approximate return map is compared
with integrated strips for mu = 1e-6, 1e-4, 1e-3 and agrees (stated "even for mu = 1e-3"), periodic orbits are
computed at mu = 1e-4, C_J = 2.8, and strips begin to merge or disappear between mu = 5e-4 and 1.5e-3. The
largest mass parameter printed is 1.5e-3.

(b) Two secondaries with different periods: nothing is said. The model has one small primary on a circular
orbit; the paper's only remark outside it is the restricted elliptic case (stated as possible, not done).

(c) Reproducible in `src/cyclerfinder/core/cr3bp.py` as a check of the flyby-chain continuation machinery
(INFERRED, my recipe, not tested): the three periodic orbits of Figures 8, 9, 10. Take mu = 1e-4, C_J = 2.8,
alpha = 0.4 (circle radius mu^0.4, about 0.0251). Initial state in synodic coordinates (INFERRED convention:
phi measured from the +x axis about M, psi the synodic velocity direction from the +x axis, which matches the
leaving condition |phi - psi| < pi/2): x = mu - 1 + mu^alpha cos(phi), y = mu^alpha sin(phi), speed
v_s = sqrt(2 Omega - C_J), (xdot, ydot) = v_s (cos psi, sin psi). Integrate (regularised if the code can; the
orbits stay outside the circle between encounters but pass through it). The orbit should close (Figure 9, with the
x-axis mirror symmetry, is the easiest first test; Figure 8 and 10 are asymmetric). The printed digits of phi
and psi are given to 15 to 20 places, so closure should be tested to the integrator's precision, not to the
digits; the stability parameters (1e6 to 1e9) are the expected sensitivity scale. The Figure 5, Figure 7 and
Figure 11 strips are figures only, with no tabulated points. The orbital periods are not printed in this
paper (they are in the 2009 companion, Tables 3 and 4). Caveat: the definition of the stability parameter and
the exact orientation convention of psi were not found in the printed text, so closure is the check, not the
stability value.

(d) Wording supported: "For the planar circular restricted three-body problem with one small secondary,
Font, Nunes and Simo (2002) proved, for sufficiently small mass parameter, that sequences of resonant arcs
joined by encounters at distance of order mu^alpha define a horseshoe, hence periodic and chaotic orbits of
consecutive quasi-collisions (their Theorem 5.1), and computed examples at mu = 1e-4, with the analytic
approximation checked up to mu = 1e-3." Also supported: the instability per encounter grows as mu^(-alpha) in
their scaling, and the unstable periodic orbits they print have stability parameters of 1e6 to 1e9, the same
order as the 8.4e5 per cycle of #890 (INFERRED; the quantities are not defined identically, so this is a
magnitude comparison only). The encounter circle of #890 is also consistent: Titania's flyby periapsis of 2,766
km against mu^(2/5) times 436,298 km (using mu of about 4e-5 for Titania from GM 226.9 over a Uranus GM of
about 5.8e6, my recollection, not read from the project note, giving a circle of order 7,600 km) puts the
flyby inside the paper's encounter circle (INFERRED, order of magnitude only).
Not supported: any statement that the paper covers two moons, two different periods, or a four-body problem;
that it proves or computes a periodic orbit at mu of order 1e-2 or above; that it addresses flybys with large
turn angles (the proved symbols are the small-deflection resonant arcs; see section 3); or that it gives a
numerical admissible mu for any N. Whether #890's two-moon orbit is a second species orbit of a four-body
problem is not a statement these authors make.

## 7. References cited by the paper for second species orbits (p. 142, READ; full list)

[1] Bolotin S V and MacKay R S 2000 Periodic and chaotic trajectories of the second species for the n-centre
problem, Celest. Mech. (submitted; published as Celest. Mech. Dyn. Astron. 77, 49-75).
[2] Henon M 1966 Sur les Orbites Interplanetaires qui Rencontrent Deux Fois la Terre, Bull. Astronomique Paris
1, 377-402.
[3] Henrard J 1980 On Poincare's second species solutions, Celest. Mech. 25, 83-97.
[4] Marco J P and Niederman L 1995 Sur la construction des solutions de seconde espece dans le probleme plan
restreint des trois corps, Ann. Inst. H Poincare Phys. Theor. 62, 211-49.
[5] Moser J 1973 Stable and Random Motions in Dynamical Systems (Princeton: Princeton University Press).
[6] Perko L 1976 Second species periodic solutions with an O(mu) near-moon passage, Celest. Mech. 14, 395-427.
[7] Poincare H 1899 Les Methodes Nouvelles de la Mecanique Celeste Tome III (Paris: Gauthier-Villars).
[8] Szebehely V 1967 Theory of Orbits (New York: Academic Press).
(Note: the volume of Henrard in this list is 25; the Bolotin and MacKay 2006 digest cites it as 21, 83-97
(1980), and the 2009 companion lists it as Celest. Mech. 21. The two papers by the same authors differ; not
resolved here.)
