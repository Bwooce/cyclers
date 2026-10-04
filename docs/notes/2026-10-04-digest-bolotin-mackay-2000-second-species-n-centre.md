# Digest: Bolotin and MacKay 2000, second species periodic and chaotic trajectories for the n-centre problem

Date 2026-10-04. Purpose: establish exactly what the general n-centre theorem requires, so the project can word
how far published "second species" theory covers the #890 object (a periodic orbit of the planar concentric
circular restricted four-body problem, Uranus plus Titania plus Oberon on circular orbits with different
periods, one close flyby of each moon per cycle, found by continuing a patched-conic flyby chain in the moon
mass; see `docs/notes/2026-10-04-890-titania-oberon-candidate.md` and
`docs/notes/2026-10-04-890-literature-check.md`). This digest covers the part the 2006 follow-up only cites;
for the 2006 paper see `docs/notes/2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`, which is
not repeated here.

Citation: S. V. Bolotin and R. S. MacKay, "Periodic and chaotic trajectories of the second species for the
n-centre problem", Celestial Mechanics and Dynamical Astronomy 77:49-75 (2000), DOI
10.1023/A:1008393706818. 27 pages (journal pages 49 to 75; PDF page n is journal page n+48). Filed in the
private paper corpus as
bolotin-mackay-2000-periodic-chaotic-trajectories-second-species-n-centre-problem-cmda-77-49-doi-10.1023-A-1008393706818.pdf.

Evidence labels: READ means read from the printed page in this session. INFERRED means my reading of what the
printed text implies. All 27 pages were read as page images. Formulas are transcribed from the images; page
numbers are journal pages. The paper has no closing "comments" or "extensions" section: it ends with Section 6
(regularisation) and the references. All statements about extensions are the two Remarks after Theorem 2.1 and
the remarks inside Section 1.

## 0. Summary (READ)

Abstract (p. 49): "For the n-centre problem of one particle moving in the potential of attracting centres of
small mass fixed in an arbitrary smooth potential and magnetic field, we prove the existence of periodic and
chaotic trajectories shadowing sequences of collision orbits. In particular, we obtain large subshifts of
solutions of this type for the circular restricted 3-body problem of celestial mechanics. Poincare had
conjectured existence of the periodic ones and given them the name 'second species solutions'."

Structure: Section 1 states the general theorem (Theorem 1.1, planar or spatial, n fixed centres). Sections 2
and 3 apply it to the planar circular restricted three-body problem (Theorem 2.1 there, with Lemmas 3.1 to
3.3). Sections 4 to 6 prove Theorem 1.1 (Lemma 4.1, Theorem 5.1, regularisation Lemmas 6.1 to 6.4, Proposition
6.1). The word "fixed" in the abstract is load-bearing: the centres are fixed points of the configuration
manifold.

## 1. The general setting (pp. 49-52, READ)

### 1.1 The Lagrangian and the centres (p. 49)

"Let Q be a smooth manifold and P = {p_1, ..., p_n} a finite set in Q. Consider a Lagrangian system (L_eps)
with configuration space Q \ P and Lagrangian

  L_eps(q, qdot) = L(q, qdot) - eps V(q).   (1.1)

We assume that L is C^4 on TQ and quadratic in the velocity:

  L(q, qdot) = T(q, qdot) + <omega(q), qdot> - W(q),   (1.2)

where T = <A(q) qdot, qdot>/2 is a Riemannian metric on Q, and that V is a C^4 function on Q \ P having
Newtonian singularities at the points p_i in P. This means that in a neighbourhood U_i of any point p_i in P,

  V(q) = - f_i(q) / dist(q, p_i),   (1.3)

where f_i is a C^4 function on U_i with f_i(p_i) > 0, and the distance is defined by means of the Riemannian
metric T. We call the system (L_eps) with the Lagrangian (1.1) the n-centre problem. We will study this problem
for small eps > 0 as a perturbation of the Lagrangian system (L) with the Lagrangian (1.2). Examples are given
by the classical integrable 2-centre problem of Euler and by the restricted circular 3-body problem in the
rotating frame, when one of the masses is small. For physical reasons (and also to simplify regularisation of
singularities) we assume that Q is two or three-dimensional."

Energy (p. 50): "Let H_eps = H + eps V, H(q, qdot) = T(q, qdot) + W(q), be the energy integral. We fix E such
that the domain D = {q in Q | W(q) < E} contains the set P and study the system (L_eps) on the energy level
{H_eps = E} in T(Q \ P)."

Answers to the setting questions:
- How many centres: a finite number n, any n. Corollary 1.1 (p. 52) is stated "for any n >= 2".
- Fixed: yes. P is a fixed finite subset of Q and L, V are functions of (q, qdot) and q only. The abstract
  says "centres of small mass fixed". Nothing in Section 1 gives a centre a motion in time.
- May they move in time: the Lagrangian (1.1) to (1.3) has no time argument and H_eps is "the energy integral"
  (p. 50). No remark in the paper allows explicit time dependence of L or V. (READ for the form of the
  Lagrangian; "no remark allows" is a statement about absence, checked across all 27 pages.)
- Several centres with different small masses: V is one function carrying all the singularities, multiplied
  by one eps. The functions f_i may differ (only f_i(p_i) > 0 is required), so centres need not have equal
  strength; they scale together with the single eps. (READ for (1.3); the remark about relative strengths is
  INFERRED from f_i being arbitrary positive C^4 functions.) The theorem's constants depend on the set of
  collision orbits, not on any ratio between the f_i.
- Role of eps: sole small parameter, multiplying the whole singular potential. Remark 2 (p. 52): "One can
  allow the nonsingular part L of the Lagrangian L_eps also to depend on eps. Then all the results remain true
  with L replaced by L|_{eps=0}. Also V may depend on eps." That is dependence on the small parameter, not on
  time.
- Dimension: Q two or three-dimensional (p. 49). Remark after Lemma 4.1 (p. 64): "The proof of Lemma 4.1 is
  the only place in the proof of Theorem 1.1 where the assumption dim Q = 2 or 3 is used. ... Hence Theorem 1.1
  holds also in the multidimensional case".

### 1.2 Definitions (pp. 50-52)

Collision trajectory (p. 50): "Suppose that the system (L) with the Lagrangian (1.2) has nondegenerate
(definition to be recalled shortly) collision trajectories gamma_k: [0, tau_k] -> D, k in K (a finite set),
of energy E such that gamma_k(0), gamma_k(tau_k) in P and gamma_k does not pass through other points of the
set P. More precisely, gamma_k(0) = p_{alpha_k}, gamma_k(tau_k) = p_{beta_k}, and gamma_k((0, tau_k)) in
D \ P. A sequence (k_i in K)_{i in Z} will be called a chain if gamma_{k_i}(tau_i) = gamma_{k_{i+1}}(0) and
gammadot_{k_i}(tau_{k_i}) != +- gammadot_{k_{i+1}}(0) for all i in Z."

Nondegeneracy (p. 51), "expressed in five equivalent ways": (1) the end points are not conjugate along the
collision orbit on the energy level {H = E}; (2) rank dg(gammadot_k(0), tau_k) is maximal, equal to dim Q,
for g(v, tau) = pi(phi_tau(v)) = p_beta on the sphere S_alpha = T_{p_alpha}Q intersect {H = E}; (3) the
collision orbit is a nondegenerate critical point of the action F(u, tau) = integral_0^tau (L + E) dt (1.7);
(4) the same for the Maupertuis-Jacobi functional J(gamma) = integral_0^1 g_E dt with g_E(q, qdot) =
2 sqrt((E - W(q)) T(q, qdot)) + <omega(q), qdot>, "(g_E is called the Jacobi metric)"; (5) the rank of the
system q(lambda, 0) = p_alpha, q(lambda, tau) = p_beta, h(lambda) = E (1.8) in the variables lambda, tau is
maximal, "equals 2 dim Q + 1" (p. 52). All five use the conserved energy H or the fixed level {H = E}.

### 1.3 Theorem 1.1, quoted in full (p. 50)

"THEOREM 1.1. There exists eps_0 > 0 such that for all eps in (0, eps_0] and any chain (gamma_{k_i})_{i in Z}:
- there exists a trajectory gamma: R -> D \ P of energy E for the system (L_eps) shadowing the chain
  (gamma_{k_i})_{i in Z}, and it is unique (up to a time shift) if the W_k are chosen small enough.
- the orbit gamma converges to the chain of collision orbits as eps -> 0: there exists an increasing sequence
  (t_i)_{i in Z} such that
    max_{t_i <= t <= t_{i+1}} dist(gamma(t), gamma_{k_i}([0, tau_{k_i}])) <= C eps,   (1.4)
  where the constant C > 0 depends only on the set K of collision orbits.
- the orbit gamma avoids collision by a distance of order eps: there exists a constant c in (0, C), depending
  only on K such that
    c eps <= min_{t_{i-1} <= t <= t_i} dist(gamma(t), p_{alpha_{k_i}}).   (1.5)"

Hypotheses actually used (READ): the same energy E for all collision orbits and for the perturbed trajectory
("trajectory ... of energy E"); a fixed finite set P; the Lagrangian (1.1) to (1.3) as above; a FINITE set K of
nondegenerate collision orbits; consecutive orbits in a chain meet at the same centre (beta_{k_i} = alpha_{k_{i+1}},
graph (1.6)); direction change gammadot_{k_i}(tau_{k_i}) != +- gammadot_{k_{i+1}}(0).

Consequences stated (p. 50): "Theorem 1.1 implies that there is an invariant subset in {H_eps = E} on which the
system is a suspension of a subshift of finite type. The topological entropy is positive provided the graph with
the set of vertices K and the set of edges G = {(k, l) in K^2 | beta_k = alpha_l, gammadot_k(tau_k) != +-
gammadot_l(0)} (1.6) has a connected branched subgraph. Note that in the case of a periodic sequence
(k_i)_{i in Z}, uniqueness of the trajectory gamma implies that it is also periodic, with the minimal period
consistent with shadowing the periodic sequence."

Precise version: Theorem 5.1 (p. 66), same hypotheses, "There exists eps_0 > 0 such that for any eps in (0,
eps_0] and any chain (k_i in K)_{i in Z} of collision orbits there exists a unique (up to a time shift)
trajectory gamma ... of energy E for the system (L_eps)". It adds: "The constant eps_0 > 0 depends only on the
set {gamma_k}_{k in K} of collision orbits and is independent of the sequence (k_i in K). Thus gamma(t) is
O(eps)-close to a chain of collision orbits. By (4.8) the trajectory gamma|[b_i, a_{i+1}] avoids p_j by a
distance of order eps." Proof strategy (p. 66): "continuation from the case eps = 0, which in the language of
[4] is an 'anti-integrable limit'." The proof applies the implicit function theorem to the gradient of the
formal action F_eps on a product of neighbourhoods and requires "eps_0^{-1} > C max_{(k,l) in G} max_{B_k x A_l}
||s''||", where C is the bound in (5.1) on the inverse of the second derivative of g_k (p. 67).

Remarks (p. 52, READ):
1. Nondegeneracy may be replaced by a minimality condition (the action of any other curve in W_k from p_alpha
   to p_beta exceeds that of gamma_k), or by an isolated minimax critical point; "We will prove the
   nondegenerate version of Theorem 1.1." The uniqueness assertion is lost in the degenerate case.
2. Dependence of L, V on eps allowed (quoted above).
3. If shadowing orbits may come arbitrarily close to collision, or pass through regularised collisions, then
   the direction-change condition "gammadot_{k_i}(tau_{k_i}) != -gammadot_{k_{i+1}}(0)" (the paper writes only
   the minus-sign case) "can be skipped. Theorem 1.1 remains true with the exception of the last statement
   (inequality (1.5))."
After the remarks (p. 52): "Note that the conditions of Theorem 1.1 are formulated in terms of the Lagrangian L
and the set P only, not involving the potential V provided it has Newtonian singularities at P."

Corollary 1.1 (p. 52): "Suppose Q is a closed manifold and E > -min_{q, qdot} L(q, qdot) = min_q (1/2 <A^{-1}(q)
omega(q), omega(q)> + W(q)). Then for any n >= 2 and almost all points p_1, ..., p_n in Q, for small eps > 0,
there exist chaotic trajectories of energy E close to chains of collision orbits." The proof uses that g_E is
then a positive definite Finsler metric and Morse theory. Example given: n = 2, Q = S^2 with the standard
metric, chaotic trajectories exist "provided p_1 and p_2 are not antipodal points of S^2".

## 2. Application to the planar circular restricted three-body problem (pp. 53-62, READ)

### 2.1 Reduction to one fixed centre (pp. 53-54)

"Consider the planar restricted circular 3-body problem (Sun, Jupiter, and Asteroid) and suppose that the mass
of Jupiter is small compared to the Sun. We suppose the masses to be normalised to 1 - eps (Sun), eps
(Jupiter), and 0 (Asteroid), with the centre of mass O stationary and the first two masses in circular orbits
about O, having separation and angular frequency both normalised to 1. To apply Theorem 1.1, consider the
motion of the Asteroid in the frame rotating at angular frequency 1 about the Sun." (p. 53)

L = (1/2)|qdot|^2 + x ydot - y xdot + (1/2)|q|^2 + 1/|q|   (2.1)
V(q) = 1/|q| - 1/|q - j| + x   (2.2)

"Here, q = (x, y) in R^2 and j = (1, 0) is the position of Jupiter (the Sun is at (0, 0)). Hence L has the form
(1.2), V has the form (1.3), and Q = R^2 \ {0}." The Coriolis term x ydot - y xdot is the magnetic term
<omega, qdot> of (1.2). The energy integral is the Jacobi integral H_eps = (1/2)|qdot|^2 - (1/2)|q|^2 -
(1 - eps)/|q| - eps/|q - j| + eps x, with "The energy constant H_eps = E is traditionally replaced by the Jacobi
constant C = -2E" and C = 2(h - script-E) (script-E the energy in the fixed frame, h the angular momentum about
O) (pp. 53-54).

So the three-body problem is a ONE-centre problem: the single set P is {j}, Jupiter, stationary in the frame
rotating with it. The Sun is not a centre of the perturbation; its 1/|q| sits in W and the point 0 is removed
from Q. Sign note: comparing (2.1) with (1.2), the W of (1.2) is -(1/2)|q|^2 - 1/|q| (INFERRED); (2.1) is transcribed
as printed.

### 2.2 Theorem 2.1 (p. 54), quoted in full

"THEOREM 2.1. In the planar circular restricted 3-body problem with masses 1 - eps, eps, 0, for all values of
the Jacobi constant C in (-sqrt(8), +3) there exists a dense subset S_C of rationals in the set A_C of allowed
frequencies for Kepler ellipses crossing the unit circle, such that for all finite subsets T subset S_C there
exists eps_0 > 0 such that for any sequence sigma = (Omega_n = m_n/k_n)_{n in Z} in T and 0 < eps < eps_0 there
is a unique trajectory of Jacobi constant C near to a chain of collision trajectories formed by transforming
ellipses of frequencies Omega_n traversed m_n times to the rotating frame, and it converges to the chain as
eps -> 0."

Arcs used (pp. 55-56): Kepler ellipses about the Sun, with Jacobi constant C = a^{-1} +- 2 sqrt(a(1 - e^2))
(3.1) (plus for prograde, minus for retrograde), mean motion Omega = a^{-3/2} (3.3), crossing the unit circle
when (C - a^{-1})^2 < 8 - 4 a^{-1}. "Thus given C in (-sqrt(8), +3) we obtain a non-empty set A_C of frequencies
of ellipses that cross the unit circle, which consists of one interval for C in [2, 3) or (-sqrt(8), 0] and two
intervals if C in (0, 2)." Whole revolutions: "choose a so that Omega is rational, say m/k in lowest terms,
start Jupiter and the Asteroid at either of the two intersection points, and let Jupiter make k complete
revolutions and the Asteroid make m complete revolutions ..., after which they will collide again (at exactly
the same point)." These are class-1 arcs in the later paper's terminology: whole numbers of revolutions of a
coplanar orbit, each starting and ending at Jupiter. The paper's own comment: "The only catch is that they might
collide earlier, at the other intersection point of the ellipse with the unit circle, which we call an 'early
collision'."

Lemma 3.1 (p. 56): "For all C in (-sqrt(8), +3) there is a dense set S_C of rational Omega in the allowed set
A_C such that there is no early collision." Proof (pp. 56-61) by analyticity of G_C(Omega) = (M - f Omega)/pi
and case analysis.
Lemma 3.2 (p. 62): "For C in (-sqrt(8), +3), the collision orbit at eps = 0 in the rotating frame,
corresponding to a whole number m of revolutions of an ellipse with rational frequency m/k in A_C starting and
ending at collision with Jupiter, is non-degenerate." Proof via formulation 5 of nondegeneracy; the key
derivative is dC/d(SF) = +- SF / (2a sqrt(a - SF^2/4a)) != 0. The paper remarks: "It is fortunate that we want
nondegeneracy for fixed Jacobi constant rather than fixed time (or fixed energy), because any point P on a Kepler
ellipse is equal-time self-conjugate (and equal-energy self-conjugate)".
Lemma 3.3 (p. 62): "Given any sequence (Omega_n)_{n in Z} in S_C there exists an 'orbite a chocs' consisting of
ellipses of frequency Omega_n connecting collisions at Jupiter, which leave each collision in neither the same
nor the opposite direction as they arrive." This is the direction-change lemma; see section 5 below.

Range of Jacobi constant: C in (-sqrt(8), +3) for the planar problem (the 2006 spatial paper gives (-2, +3)).
Arcs: only elliptic, whole-revolution, coplanar arcs of rational frequency Omega = m/k, each beginning and ending
at the single small body.

### 2.3 Remarks on the application (p. 55, READ)

Remark 1: "Given a choice of subset T in Theorem 2.1, (1.5) shows that there exists c > 0 such that the
constructed orbits remain at least c eps from collision, so the result can even be applied to an extended
spherically symmetric Jupiter (its external field is equivalent to that of a point mass), provided its radius is
less than c eps. Similarly, the Sun can be an extended body provided the eccentricity of the ellipses with
frequencies in T is bounded away sufficiently from 1 to prevent penetration of the Asteroid into the Sun.
Spherical symmetry can even be abandoned if the Sun and Jupiter maintain the same orientations relative to each
other (as for the Earth-Moon system). The potential V still has the form (1.3) outside a sphere S_J containing
Jupiter in the concentric sphere of half the radius, and then V can be artificially extended to satisfy (1.3)
inside S_J too. As already mentioned, Theorem 1.1 can be modified to allow V to depend on eps. The resulting
solutions avoid S_J if its radius is less than c eps. The effect of asphericity of the Sun simply modifies the
potential W, so requires no further analysis."

Remark 2: "It would be interesting to see whether our Theorem 1.1 could also be used to prove existence of second
species orbits not confined to the plane of rotation of the Sun and Jupiter, which Poincare also claimed to
exist. ..." (the nonplanar case, taken up in the 2006 paper).

## 3. Statements on more than one small mass, time-periodic, elliptic, full three-body problem (READ)

Everything the paper says on these topics, exhaustively (searched across all pages):

(a) The unrestricted (two small masses) problem, p. 54, immediately after Theorem 2.1: "For periodic sequences
sigma, these solutions are a special case of the 'second species solutions' whose existence Poincare claimed to
prove in chapter XXXII of [18] (published 100 years ago). 'Second species solutions' are ones which in the
limit eps -> 0 tend to chains of collision orbits, which he called 'orbites a chocs'. He claimed to treat the
'unrestricted case' of masses 1, eps m_1, eps m_2 and to construct non-planar solutions too. His sketch proof of
all of these, however, leaves a lot to be desired (see Levy's comments on p. 629 of [19])."
The paper proves nothing for the unrestricted case; it describes Poincare's claim only.

(b) The two-centre problem appears only as an example of a system of the form (1.1): "Examples are given by the
classical integrable 2-centre problem of Euler and by the restricted circular 3-body problem in the rotating
frame, when one of the masses is small." (p. 49) and as a sphere example in Corollary 1.1 (p. 53). In both the
centres are fixed points of Q (Euler's problem has two fixed gravitating centres). There is no statement about
two orbiting small bodies.

(c) Time-periodic problems, the elliptic restricted problem, and quasi-periodic forcing: NOT mentioned anywhere
in the paper. The paper's reference list (pp. 74-75) contains Gomez and Olle (1991), "Second species solutions
in the circular and elliptic restricted three body problem, I and II" (ref. 10), cited only inside the general
list of earlier analyses on p. 54. No statement of an extension to the elliptic case is made; the later paper
(2006, p. 434) attributes the elliptic extension to Bolotin (2005, 2006).

(d) Parabolic and hyperbolic arcs, infinite sets of arcs, unbounded orbits: not mentioned in this paper (those
are the 2006 paper's "probably" statements).

(e) Regarding the regime of validity of earlier work (p. 54): "Several of the above proofs, however, are
incomplete from a mathematical viewpoint. The paper [15] gives a mathematical proof of the existence of
two-segment second species periodic orbits, meaning ones whose orbite a chocs consists of alternation between
two collision trajectories; it was not checked there that these collision trajectories do not suffer additional
collisions between the intended ones, but Marco filled this gap in a private communication to us. The phenomenon
of additional collisions really occurs for many other collision trajectories as the Jacobi constant is varied,
so is a genuine problem." (This is the "early collision" issue that Lemma 3.1 handles.) Also: "Aperiodic
sequences of collision orbits were treated in the collinear case by Saari and Xia [20], but their shadowing
orbits are constructed for a regularised system and collide just as often as the given sequence of collision
orbits, so do not generate solutions of the 3-body problem."

(f) Closing statements about extensions: none. The only forward-looking statements are the Remark 2 on p. 55
(nonplanar case, "It would be interesting to see") and the Remark on p. 64 (arbitrary dimension: "we are not
aware of any physical application of this generalisation").

## 4. Quantitative content (READ)

- Closest approach: (1.5) and (4.8): dist(gamma(t), p_i) >= c eps, with c depending only on the set K of
  collision orbits. The mechanism: the transit time near a centre is -log eps + bounded (Proposition 6.1, p. 70:
  "tau_eps(a, b) = -log eps + mu(a, b, eps)", mu uniformly bounded) and, in regularised variables, "(6.15)
  min |x_{a,b}^eps(t)|^2 = 2 eps(|v^+(a)||v^-(b)| - <v^+(a), v^-(b)>) + o(eps)" (p. 71). Because the regularising
  map squares the coordinate, the regularised distance of order sqrt(eps) is a physical distance of order eps; the
  paper writes "Hence x^eps_{a,b}(t) avoids 0 by a distance of order sqrt(eps) provided that v^+(a) != v^-(b)"
  for the regularised variable.
- Shadowing distance: (1.4), at most C eps in position over each arc, with C depending only on K; transit time
  of each arc within order eps of the collision-orbit time (Theorem 5.1: b_i - a_i -> tau_{k_i} - tau^-(x_{k_i})
  - tau^+(y_{k_i}) as eps -> 0; Lemma 4.1: tau(a, b, eps) -> tau^+(a) + tau^-(b)). The O(eps) transit-time bound
  stated explicitly in the 2006 paper's Theorem 2.1 does not appear in this paper's Theorem 1.1 as printed
  (only the limit statements above). (READ; compared against the 2006 digest.)
- Instability: the 2000 paper establishes the subshift-of-finite-type structure (p. 50) and says the
  topological entropy is positive under the branched-subgraph condition. It does NOT state the Lyapunov
  exponent estimate of order log eps^{-1} or uniform hyperbolicity as a theorem: the 2006 paper says
  "Theorem 2.2 can be deduced from the proof in Bolotin and MacKay (2000) of Theorem 2.1, but it was not proved
  there". The only hint of strong instability here is the transit time -log eps in Proposition 6.1. No Lyapunov
  exponent, Floquet multiplier or hyperbolicity statement appears. (READ.)
- Numerical bound on eps: none. "There exists eps_0 > 0" in Theorem 1.1, Theorem 2.1 and Theorem 5.1; the only
  characterisation of eps_0 is the implicit condition eps_0^{-1} > C max ||s''|| on p. 67, in terms of
  unspecified constants. The paper contains no numerical example, no mass ratio of any real system and no
  numerical value of c, C or eps_0. (READ; checked across all pages.)
- Jacobi constant range (planar CR3BP): C in (-sqrt(8), +3).

## 5. The erratum issue: direction change in the inertial frame (READ for the quotes, INFERRED for the identification)

The 2006 paper's footnote 3 (p. 444 there) says: "At the analogous point in Bolotin and MacKay (2000) we
mistakenly studied the direction change in the inertial frame; this error was corrected in MacKay (2005)." The
2006 introduction (footnote 1) adds: "See also MacKay (2005) for a summary, some minor additions and a
correction".

The passages in the 2000 paper that this refers to (INFERRED identification; the 2000 paper itself does not flag
any error, and I have not read MacKay 2005):

1. p. 54, last sentence of Section 2: "In the next section we construct a set of collision trajectories with the
   given value of C for the case eps = 0 (Lemma 3.1), check their non-degeneracy (Lemma 3.2) and construct a
   non-trivial set of chains of collision trajectories which change direction at each collision (Lemma 3.3). These
   three things are most easily done in the non-rotating frame, to which we shall now revert."
2. p. 62, Lemma 3.3 and its proof, in full: "LEMMA 3.3. Given any sequence (Omega_n)_{n in Z} in S_C there exists
   an 'orbite a chocs' consisting of ellipses of frequency Omega_n connecting collisions at Jupiter, which leave
   each collision in neither the same nor the opposite direction as they arrive. Proof. For each C in
   (-sqrt(8), 3) and Omega in A_C and position of Jupiter on the unit circle, there are two possible ellipses
   leaving it with frequency Omega, and they leave in different and non-opposite directions (see Figure 3). Thus
   at least one of them leaves in a different and non-opposite direction from that in which the Asteroid arrived
   along the previous collision trajectory." Figure 3 (p. 63) caption: "The two ellipses which leave the same
   point of the unit circle with the same values of (a^{-1}, C)."

Why this is a frame issue (INFERRED, simple geometry): Theorem 1.1's direction-change condition
gammadot_{k_i}(tau_{k_i}) != +- gammadot_{k_{i+1}}(0) is a condition on the velocities of the Lagrangian system
(L), which in Section 2 is the system in the rotating frame (Lagrangian (2.1)). Lemma 3.3 is proved for
ellipses in the non-rotating frame (Figure 3 draws the ellipses about the Sun with Jupiter's orbit as the unit
circle). The rotating-frame velocity is the inertial velocity minus a rotation of the position; at a collision
the position is the same for the arriving and leaving arcs, so the correction is the same vector for both, but
parallel or antiparallel in one frame does not imply parallel or antiparallel in the other. The 2006 paper,
Section 4.5, redoes the check in the rotating frame (using the relative velocity to Jupiter, with
v_phi = G_z = cos(iota)) and uses the sign of the mirror ellipse; its statement of the result is the 2006
digest's section 3.3.

Consequence for readers: the planar Theorem 2.1 here is stated with a proof of its direction-change step that
the authors later said was done in the wrong frame. What the correction in MacKay (2005) concludes about the
planar Theorem 2.1 (whether it stands unchanged, stands with an extra condition, or is repaired) is not
established by the 2000 or 2006 papers as I read them; "correction" in 2006 footnote 1 and "error was
corrected" in footnote 3 are all that is printed. Anyone citing the 2000 planar theorem should cite it together
with MacKay (2005) (Equadiff 2003, pp. 59-72), which the project does not hold (see the 2006 digest, section 8,
item 4).

## 6. What this does and does not cover for #890

### 6.1 Hypotheses of Theorem 1.1 against the #890 setting

The #890 model: Uranus fixed (concentric model), Titania and Oberon on concentric circular orbits with
different periods, a massless particle (planar). The time dependence in any rotating frame is that of the
second moon unless the periods are equal.

(a) Can the planar problem of one large mass and TWO small masses on circular orbits with different periods be
put into the paper's setting (n fixed centres in a time-independent Lagrangian with conserved energy) by any
change of frame? From the printed hypotheses, no (INFERRED from READ hypotheses). The paper's own reduction in
Section 2 works because one small body on a circular orbit about a fixed large one is stationary in the frame
that rotates with it; that gives a Lagrangian (2.1) with no time argument, hence a conserved energy (the Jacobi
integral) and a single fixed centre. With two moons of different periods no uniformly rotating frame makes both
stationary: in the frame of one, the other circulates at the difference of the angular rates, so the potential
depends explicitly on time (periodically if the periods are commensurate, quasi-periodically otherwise). The
Lagrangian (1.1) to (1.3) does not admit that, there is no energy integral H_eps, and so "the same energy E",
the Jacobi metric g_E, nondegeneracy formulations 2 to 5 (pp. 51-52), the set D = {W < E} and Theorem 1.1 itself
are not defined. The only escape written in the paper, Remark 2 (p. 52), lets L and V depend on eps, not on
time. If the two moons had the SAME period (a co-orbital pair, or any configuration of the large and small
masses rotating rigidly), then (INFERRED) there would be a frame in which both moons are fixed points, the
Lagrangian would be time-independent with the Coriolis and centrifugal terms in the form (1.2), the set P would
have n = 2 points with f_1, f_2 reflecting the two masses, Q would be the plane minus the large body, and the
n-centre setting of Theorem 1.1 would apply in form, subject to the remaining hypotheses (a finite set of
nondegenerate collision arcs between the two centres at one common Jacobi constant, with direction change, and
closeness to the singular limit eps -> 0 with both moon masses scaled by a single eps). That equal-period case is
not the #890 case (Titania and Oberon have different periods), and the paper does not discuss it.

(b) Wording the paper supports and does not support for #890 (all statements about the paper, not about #890):

Supported:
- "Bolotin and MacKay (2000) prove, for n fixed centres in a time-independent Lagrangian with a conserved energy
  and for the planar circular restricted three-body problem with one small secondary, that for sufficiently small
  eps there exist periodic and chaotic trajectories that shadow chains of nondegenerate collision arcs of one
  common energy (Jacobi constant), with closest approach of order eps."
- "The theorem requires nothing about the size of eps beyond existence of some eps_0; the paper gives no
  numerical eps_0 and no real-system mass ratio."
- "The paper's n-centre theorem is stated for fixed centres; its restricted three-body application has one
  centre."

Not supported:
- That the n-centre theorem covers several small bodies moving on circular orbits with different periods. "n
  centres" means n fixed points of the configuration space Q, not n moving secondaries.
- That the paper proves, or conjectures, or leaves open, a two-secondary (restricted four-body) version, a
  time-periodic or quasi-periodic version, or an elliptic version. It states none of these. The only mention of
  more than one small mass is Poincare's unrestricted claim of masses 1, eps m_1, eps m_2, reported as unproved
  (p. 54).
- That the #890 orbit "exists by Bolotin and MacKay". The #890 orbit was computed by continuation, not obtained
  from the theorem. In addition, the theorem needs a common energy and the #890 model has no conserved energy.
- That eps = the Titania or Oberon mass ratio is within the theorem's range. No eps_0 is given.
- Any statement about instability (Lyapunov exponent of order log eps^{-1}) from this paper alone: that
  statement is the 2006 paper's Theorem 2.2, not proved here.
- Any use of the planar Theorem 2.1 without the MacKay (2005) correction (section 5 of this digest).

Safe one-line wording (INFERRED from the above): "Bolotin and MacKay (2000) proved the existence of periodic and
chaotic orbits that shadow chains of collision arcs for n FIXED centres in a time-independent Lagrangian with a
conserved energy, and for the circular restricted three-body problem with one small secondary; a two-moon
problem with different orbital periods has no conserved energy in any frame and is outside the printed
hypotheses."

## 7. Related papers not held (checked against the 2006 digest's section 8)

Bolotin, Second species periodic orbits of the elliptic 3-body problem, CMDA 93 (2005) 345-373 (the
time-dependent extension; the nearest published non-autonomous analogue); Bolotin, Shadowing chains of
collision orbits, DCDS 14 (2006) 235-260; MacKay, Chaos in three physical systems, Equadiff 2003 (2005) 59-72
(the correction); Marco and Niederman, Ann. IHP Phys. Theor. 62 (1995) 211-249 (ref. 15 here); Gomez and Olle,
CMDA 52 (1991) 107-146, 147-166 (ref. 10 here); Henrard, CMDA 21 (1980) 83-97 (ref. 12); Perko, CMDA 24 (1981)
155-171 (ref. 17); Bruno, CMDA 24 (1981) 255-268 (ref. 7); Henon, Generating families in the restricted three
body problem, Springer (1997) (ref. 11, "for a nice survey" of numerics, p. 54); Poincare, Methodes Nouvelles
III (1899) chapter XXXII (ref. 18) and Levy's note, Oeuvres VII (1952) p. 629 (ref. 19); Saari and Xia, J. Diff.
Eq. 82 (1989) 342-355 (ref. 20, collinear case). New in this list relative to the 2006 digest: none beyond the
Henon survey, which the 2000 paper recommends for the planar numerics.
