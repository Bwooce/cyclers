# Digest: Dimare 2010, chaotic quasi-collision trajectories in the 3-centre problem

Date 2026-10-04. Purpose: establish exactly what this paper treats, so the project can word how far published
"second species" theory covers the #890 object (a periodic orbit of the planar concentric circular restricted
four-body problem, Uranus plus Titania plus Oberon on circular orbits with different periods, one close flyby of
each moon per cycle; see `docs/notes/2026-10-04-890-titania-oberon-candidate.md` and
`docs/notes/2026-10-04-890-literature-check.md`, section 5 of which says only that "Dimare's 3-centre result
points that way"). Companion digests: `2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`,
`2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`,
`2026-10-04-digest-marco-niederman-1995-seconde-espece-plan-restreint.md`.

Citation: L. Dimare, "Chaotic quasi-collision trajectories in the 3-centre problem", arXiv:0911.3557v2 [math-ph]
(4 May 2010), 22 pages; Dipartimento di Matematica, Universita di Roma "La Sapienza". Journal reference (from a
web search result, not from the PDF, which carries no journal data): Celestial Mechanics and Dynamical Astronomy
107:427-449 (2010), DOI 10.1007/s10569-010-9284-4. Filed in the private paper corpus as
dimare-2010-chaotic-quasi-collision-trajectories-3-centre-problem-arxiv-0911.3557v2.pdf. Text layer present
(born-digital, about 12,990 words by `pdftotext | wc -w`).

Evidence labels: READ means read from the printed page in this session. INFERRED means my reading of what the
printed text implies. All 22 pages were read as page images. Page numbers are the arXiv PDF page numbers (the
paper's own numbering is the same). The paper has no conclusion or outlook section: it ends with Section 5.3,
the acknowledgement, Appendix A and the references. All statements about extensions are therefore the passages
quoted in section 5 below.

## 0. Summary (READ)

Abstract (p. 1): "We study a particular kind of chaotic dynamics for the planar 3-centre problem on small
negative energy level sets. We know that chaotic motions exist, if we make the assumption that one of the
centres is far away from the other two (see [7]): this result has been obtained by the use of the
Poincare-Melnikov theory. Here we change the assumption on the third centre: we do not make any hypothesis on
its position, and we obtain a perturbation of the 2-centre problem by assuming its intensity to be very small.
Then, for a dense subset of possible positions of the perturbing centre in R^2, we prove the existence of
uniformly hyperbolic invariant sets of periodic and chaotic almost collision orbits by the use of a general
result of Bolotin and MacKay (see [4], [5]). To apply it, we must preliminarily construct chains of collision
arcs in a proper way. We succeed in doing that by the classical regularisation of the 2-centre problem and the
use of the periodic orbits of the regularised problem passing through the third centre."

In one sentence: two large fixed centres of equal intensity plus ONE small fixed centre C, energy small and
negative, planar; the unperturbed arcs are periodic orbits of the integrable (Euler) 2-centre problem that pass
through C; chains of such arcs, all starting and ending at the single small centre C, are shadowed (Bolotin and
MacKay 2000, 2006) by periodic and chaotic orbits making repeated close passages to C. Nothing in the paper
alternates between two different small centres.

## 1. The setting exactly (READ)

- Dimension and centres (p. 1): "We consider the motion of a particle in the plane, under the gravitational
  action of three point masses at fixed positions (the planar restricted 3-centre problem). We fix a Cartesian
  reference system Oxy on the plane and choose suitable dimensionless coordinates such that the two centres
  with greater masses occupy the positions C1 = (1,0), C2 = (-1,0). Following the common terminology, we refer
  to these as primaries. We suppose for simplicity that the primaries have equal intensities a1 = a2 = a > 0
  (symmetric problem), and we always consider all the three centres having positive intensities."
- Third centre (p. 1): "Let C = (x0, y0) in R^2 minus {C1, C2} be the position of the third centre and
  epsilon > 0 be its intensity. We assume epsilon to be very small and consider the limit epsilon to 0: this
  means that epsilon is a perturbation parameter. In other words, we make the hypothesis that the conditions are
  such that we can consider the problem as a one-parameter perturbation of an integrable one: the 2-centre
  problem."
- So: 3 centres, all FIXED in an inertial plane; 2 large (equal intensity a, at (+1,0) and (-1,0)), 1 small
  (intensity epsilon, at arbitrary fixed position C). The particle is massless and does not act on the centres.
- Autonomous with conserved energy: yes. Lagrangian (p. 2) L_eps = L0 + eps / sqrt((x-x0)^2 + (y-y0)^2), with
  L0 = (xdot^2 + ydot^2)/2 + a/sqrt((x+1)^2 + y^2) + a/sqrt((x-1)^2 + y^2); Hamiltonian H_eps = H0 + eps V,
  H0 = (px^2 + py^2)/2 + W, W the potential of the primaries and V = -1/sqrt((x-x0)^2 + (y-y0)^2). The energy
  H_eps = E is the fixed level.
- Energy regime (p. 2): "For negative values of the energy H_eps = E < 0, with E to 0, we prove the existence of
  chaotic motions passing arbitrarily close to the perturbing centre C." On p. 3: "We too study the case of
  small negative energies, but our point of view is different, in that the position of the third centre does not
  go to infinity. This allows us to study the orbits which undergo close encounters with the perturbing centre
  C." The energy must lie in (-E0, 0) with E0 > 0 depending on C and the rational set I (Theorem 3).
- Planar only. The paper is entirely in R^2.
- Small parameters: there are two linked ones. (i) epsilon, the intensity of C (Theorem 3: epsilon in
  (0, epsilon0)). (ii) the energy parameter beta = |E|/E1 (with E1 the regularised separation energy; section
  3.2, p. 8) or equivalently |E|; the construction needs beta in (0, beta0), "to require beta small enough is
  equivalent to ask for the absolute value of the energy E to be small enough" (p. 21). The order of limits is:
  fix C, fix a finite set of rationals, choose E close to zero, then take epsilon small enough (Theorem 3).
- Nondimensionalisation: the primaries' intensity a is kept as a free parameter in the formulas; the primaries'
  separation is 2.
- Contrast with Bolotin and Negrini [7] (p. 3): "The authors study the restricted 3-centre problem on the plane,
  when the third centre is very far from the other two and consider small negative energies E, in the limit
  E to 0. Then they have a two-parameter perturbation of the 2-centre problem on the zero-energy level. They
  succeed in applying the Poincare-Melnikov theory, thus proving the existence of chaotic motions." Dimare
  instead keeps the third centre at a finite, generic position and makes its intensity small.

## 2. Main theorems (READ, transcribed)

### Theorems 1 and 2 (p. 5), quoted from Bolotin and MacKay, used as black boxes

Definitions (p. 4): a solution gamma: [0,T] to M of fixed energy E for the unperturbed system L0 "will be called
a collision arc if gamma(0) = gamma(T) = C, and gamma(t) is not equal to C for any t in (0,T). In particular, the
latter condition means that there are no early collisions." Fix E < 0 with C in D = {W < E}; collision arcs are
critical points of the Maupertuis-Jacobi functional J_E(gamma) = integral of g_E(gamma, gammadot) dt with
g_E = sqrt(2(E - W)(xdot^2 + ydot^2)) on the space Omega of W^{1,2} curves in D starting and ending at C. An arc
is nondegenerate if it is a nondegenerate critical point of J_E. A finite set of nondegenerate collision arcs
gamma_k: [0,T_k] to D, k in K, with the same energy E; a sequence (k_i), i in Z, is a collision chain if
"direction change" holds: gammadot_{k_i}(T_{k_i}) is not equal to +- gammadot_{k_{i+1}}(0) for any i
(equation (1) defines the graph Gamma on K).

Theorem 1 (Bolotin-MacKay, 2000), p. 5: "Given a finite set K of nondegenerate collision arcs with the same
energy E, there exists epsilon_0 > 0 such that for all epsilon in (0, epsilon_0] and any collision chain
(gamma_{k_i})_{i in Z}, k_i in K, there exists a unique (up to a time shift) trajectory gamma: R to D minus {C}
of energy E of system (L_eps), which shadows the chain (gamma_{k_i})_{i in Z} within order epsilon. More
precisely, there exist constants B, B' > 0, independent of epsilon and the collision chain, and a sequence of
times (t_i), i in Z, such that |t_{i+1} - t_i - T_{k_i}| <= B epsilon, dist(gamma(t), gamma_{k_i}([0, T_{k_i}]))
<= B epsilon, for t_i <= t <= t_{i+1}, and dist(gamma(t), C) >= B' epsilon."

Consequence stated by Dimare: an invariant subset Lambda_eps of the energy shell {H_eps = E} on which (L_eps) is
"a suspension of a subshift of finite type", the subshift being the shift on paths of the graph Gamma; "The
important fact about the invariant set Lambda_eps is that it is uniformly hyperbolic."

Theorem 2 (Bolotin-MacKay, 2006), p. 5: "There exists a cross-section N in {H_eps = E}, such that the
corresponding invariant set M_eps = Lambda_eps intersect N of the Poincare map is uniformly hyperbolic with
Lyapunov exponents of order log eps^{-1}. In particular, the set Lambda_eps is uniformly hyperbolic as a
suspension of a hyperbolic invariant set with bounded transition times."

### Theorem 3 (p. 6), the paper's main result

"Let I a subset of Q+ be a finite set of positive rationals. There exists a dense open subset X_I of M of
possible positions for the third centre C, such that, fixed C in X_I, the following is true. There is a small
value E_0 > 0, depending on C and I, such that, fixed an energy value E in (-E_0, 0), we have that there is
epsilon_0 > 0, such that for any epsilon in (0, epsilon_0) and any sequence (q_k), k in Z, q_k in I, there exists a
trajectory of the planar restricted 3-centre problem (L_eps) on the energy shell {H_eps = E}, which avoids
collision with the third centre C by order epsilon and is within order epsilon a concatenation of pieces of
periodic orbits for the planar restricted 2-centre problem (L0), passing through C and of classes q_k, k in Z
(see Subsection 3.3 for notation). The resulting invariant set formed by these orbits is uniformly hyperbolic."

Follow-up text (p. 6): "In particular, as it follows from Theorem 2, fixed a small enough energy value E < 0 and
epsilon > 0, there is a cross section in the energy shell, such that the associated Poincare map has a chaotic
invariant set with Lyapunov exponents of order log epsilon^{-1} and it contains infinitely many periodic orbits,
corresponding to periodic collision chains. The topological entropy of the Poincare map is O(epsilon)-close to
that of the topological Markov chain associated to the graph Gamma, defined by (1). It's easy to verify that the
finite set of collision arcs that we construct in the proof of Theorem 3 determine positive topological entropy
(see Appendix A)." And: "Theorem 3 is still true if we substitute in the statement the set X_I with a set X,
which is dense in M and is independent of the set of rationals I. In this case, we do not know if the set X is
open or not, the only thing that we can say about it is that it is dense (see Remark 6)."

Answers to the questions in the brief:
- Kinds of orbit: periodic orbits (infinitely many, from periodic collision chains) and chaotic orbits
  (subshift of finite type), uniformly hyperbolic, on the energy shell E in (-E_0, 0).
- Do they pass close to ONE small centre repeatedly or alternate between centres? ONE. Every collision arc starts
  and ends at the single small centre C; the "chain" is a sequence of such arcs of classes q_k joined at C. The
  letters of the alphabet are arcs (class q and one of four transverse velocity choices), not centres. There is
  no second small centre in the model.
- Symbolic dynamics: shift on paths of the graph Gamma whose vertices are the finite set of collision arcs and
  whose edges are the direction-change-admissible pairs (equation (1)). Appendix A (p. 22) shows h_top >= log 2.
- Scaling: shadowing within order epsilon (Theorem 1: B epsilon), closest approach to C at least B' epsilon (and,
  by shadowing a collision arc, of order epsilon; the paper states "avoids collision with the third centre C by
  order epsilon"), Lyapunov exponents of order log(1/epsilon) (Theorem 2). Transition-time error at most B epsilon.
- Numerical bound: none. epsilon_0, E_0, beta_0, B, B' are existence constants; no value is computed, estimated
  or illustrated anywhere. READ: no table, no worked numerical case with an epsilon value.

## 3. Method (READ)

Strategy (pp. 3-4): apply Bolotin-MacKay (Theorem 1 and 2) with the unperturbed problem being the 2-centre problem
and the singular perturbation being the Newtonian potential of C. Dimare's own content is constructing the finite
set K of nondegenerate collision arcs and verifying direction change; the shadowing and hyperbolicity come from
[4], [5]. He does not use Bolotin and Negrini's Melnikov argument (that is the contrast case [7]) and does not
use a first-return-map horseshoe (the Font-Nunes-Simo route [11], [12], discussed on pp. 3-4: "in this manner we
don't get an horseshoe map, but we still have a symbolic dynamics on the set of quasi-collision orbits found").

Step A, regularisation (section 3.1, pp. 6-7). Elliptic coordinates x + i y = cosh(xi + i phi), cylinder
R x S^1, primaries at (xi, phi) = (0,0) and (0,pi). L_eps = (xidot^2 + phidot^2)/2 (cosh^2 xi - cos^2 phi)
- W(xi,phi) - eps V(xi,phi), with W = -2a cosh(xi)/(cosh^2 xi - cos^2 phi) and
V = -1/sqrt(cosh^2 xi - sin^2 phi + (x0^2 + y0^2) - 2(x0 cosh xi cos phi + y0 sinh xi sin phi)) (printed form;
the second term under the root is as I read it from the image). Time change tau with dtau = dt/(cosh^2 xi -
cos^2 phi) (equation (2)); regularised Hamiltonian
calH_eps = (p_xi^2 + p_phi^2)/2 - 2a cosh xi - [E - eps V](cosh^2 xi - cos^2 phi) on the zero level.

Step B, separated problem (section 3.2, pp. 7-8). At eps = 0 the regularised system separates (equation (3)):
(xi')^2/2 - 2a cosh xi + |E| cosh^2 xi = -E1 and (phi')^2/2 - |E| cos^2 phi = E1. Assumptions (4):
|E| < a, E1 > 0, |E| + E1 < 2a. Scaled: A1 = E1/(2a), beta = |E|/E1, giving (5): A1, beta > 0, 2 beta A1 < 1,
A1 < 1/(1 + beta), and (6): (xi')^2/(4a) = cosh xi - beta A1 cosh^2 xi - A1,
(phi')^2/(4a) = beta A1 cos^2 phi + A1. The xi-motion is periodic (two inversion points xi_- and xi_+, with
cosh(xi_+-) = -(a/E)(1 +- sqrt(1 + E E1/a^2))); the phi-motion is a pendulum rotation.

Lemma 1 (p. 9): periods T1 = 2 sqrt(2 a^{-1}) K(kappa1) / (1 - 4 beta A1^2)^{1/4} and
T2 = 2 sqrt(a^{-1}) K(kappa2) / sqrt(A1 (1 + beta)), with kappa1^2 = (A1(1 - beta) + sqrt(1 - 4 beta A1^2)) /
(2 sqrt(1 - 4 beta A1^2)), kappa2^2 = beta/(1 + beta), K the complete elliptic integral of the first kind.
Limits (8) as A1 to 0 and to 1/(1+beta) (T2 infinite at A1 to 0, T1 infinite at A1 to 1/(1+beta)); T2 strictly
decreasing and T1 strictly increasing in A1.

Step C, periodic orbits (section 3.3, pp. 9-10). Proposition 1: for beta in (0,1) and each positive rational q
there is a unique A1_hat(beta, q) in (0, 1/(1+beta)) with q T1 = T2 (equation (9)); the orbit has energy
E = -2 a beta A1_hat. Proposition 2: A1_hat(beta, .) is strictly decreasing in q, with limits 1/(1+beta) at
q to 0+ and 0 at q to infinity. The class q = m/n (coprime) records the ratio of oscillations in xi to
revolutions in phi. Different classes have different energies at fixed beta, which is why a finite set of classes
needs the energy chosen small enough afterwards (p. 21).

Step D, orbits through C (section 4.1, p. 11). Proposition 3: for fixed q there is beta_0 > 0 such that for
beta in (0, beta_0) there is a periodic orbit of class q through C (C lies inside the ellipse xi = xi_+ as
beta to 0 because cosh xi_+ to infinity), and xi_0 is not an inversion point. Two orbits per class, differing by
the sign of the velocity (xi'_0, phi'_0) at C (the reversed-direction pair counts as the same orbit).

Step E, avoiding the primaries (section 4.2, pp. 11-16), "the central result of the paper" (p. 12): Proposition 4
(orbit through a primary: n odd gives both primaries in a period at half-period separation, n even gives one),
Proposition 5 and Remarks 2, 3 (q = 1 case), Lemma 2 (A1_hat smooth in beta, extendable smoothly to beta = 0,
A1_hat(0,q) in (0,1)), and Theorem 4 (p. 14): "Let q in Q+ a fixed positive rational number. There is a dense
open subset X'_q of R x S^1, such that fixed (xi_0, phi_0) in X'_q, there is beta_0 > 0, such that for each
beta in (0, beta_0), the periodic orbits through (xi_0, phi_0), associated with A1_hat(beta, q), do not pass
through the primaries C1, C2. In particular, after scaling time, they are orbits of the not-regularised 2-centre
problem with Lagrangian L0, and they have energy E = -2 a beta A1_hat." Proof: a finite set S of rationals and
smooth functions G+- (xi_0, phi_0, beta) = (P +- Q)/T1 with passage through a primary requiring G+- in S; the
excluded positions are those with G+-(xi_0, phi_0, 0) in S, an open dense set of full Lebesgue measure remains.
Remark 4 (p. 15): for q = 1 at least one of the two orbits through C misses the primaries. Remark 5 (p. 16):
the theorem "would be improved" if one showed that the derivatives of G+- vanish only at isolated points or never
vanish; "At present, we don't have any of these improvements." Corollary 1 (p. 16): finite I, dense open X'_I.
Remark 6 (p. 16): by Baire, X = intersection over all rational q of X'_q is dense, independent of I, "we don't
know if the set X is open"; "Note that the choice of a sufficiently small beta cannot be independent of I,
instead."

Step F, early collisions (section 4.3, pp. 16-18). An orbit of the 2-centre problem starting from C can pass
through C again before a period; "In general, we can only say that early collisions cannot be excluded: we don't
have any knowledge about the conditions that determine them. Anyway, it is not a problem for the proof of
Theorem 3" (p. 18): the collision arc is then the partial arc from C to its first next return to C. Fully
characterised only for q = 1 (p. 17): "if q = 1, we can say that there is an early collision at C if and only if
C lies on the y-axis". Proposition 6 (p. 18): an early return forces the two orbits through C to coincide up to
reversal, with a transverse autointersection at C. The velocity map from elliptic to Cartesian coordinates is
v = U(xi, phi)(xi', phi')^T with U = [[sinh xi cos phi, -cosh xi sin phi], [cosh xi sin phi, sinh xi cos phi]]
(equation (10)), invertible away from the primaries; U(-xi,-phi) = -U(xi,phi).

Step G, nondegeneracy check (section 5.1, pp. 19-20). Sufficient condition from p. 5: collision arc with energy E
corresponds to a solution (v0, T) of f(C, v0, T) = C, H0(C, v0) = E, and is nondegenerate if the Jacobian of this
system is nonzero (the statement "actually this is the characterisation" is Dimare's). He reduces to variables
(beta, A1): F(beta, A1) = [q T1 - T2](beta, A1) = 0 and -2 a beta A1 = E, with Jacobian matrix
J = [[dF/dbeta, dF/dA1], [-2 a A1, -2 a beta]]. Since dF/dA1 > 0 (strictly, from the monotonicity of T1, T2) and
A1_hat in (0, 1/(1+beta)), the determinant is nonzero at beta = 0, A1 = A1_hat(0,q) and by regularity for small
beta. Conclusion (p. 20): "We conclude that the collision arc gamma is nondegenerate for beta small enough."
The verification is therefore for small beta only, via continuity from beta = 0; no explicit beta_0.

Step H, direction change (section 5.2, pp. 20-21). Two orbits through C of the same A1_hat with velocities
(xi'_0, phi'_0) and (-xi'_0, phi'_0) are transverse since xi'_0 is not zero; transversality survives the time
reparametrisation and the invertible matrix U. At the other elliptic image (-xi_0, -phi_0) the possible
velocities are the same. Hence exactly four collision arcs per (beta, q), in two pairs of opposite initial
velocities, "each arc in one pair has transverse initial velocity to each arc in the other" (Proposition 7,
pp. 20-21; I read it for X'_I and finite I).

Step I, proof of Theorem 3 (section 5.3, p. 21). E(beta, q) = -2 a beta A1_hat(beta, q) is strictly decreasing in
beta for beta < beta_0 with E(0,q) = 0, E(beta_0, q) = -E_0 < 0. Quote (p. 21): "For beta fixed, the energy
increase with the class q (see Proposition 2). This means that we cannot state a general result valid for a fixed
energy E and any q in Q+. On the other hand, we don't need such a result, because to apply Theorem 1 we only want
a finite number of collision arcs. What is certainly true is that for any finite set of classes I of cardinality
i, we can choose a small enough energy value E < 0 to form a set of 4i collision arcs of energy E through the
centre C." Then: arcs of different classes cannot share A1_hat, so (xi'_0, phi'_0) differ, but "we can't state
that collision arcs of different classes determine a different set of directions at the point C. What we can
certainly assure is that if we fix an arc of class q in I, then for each class q' in I, not necessarily different
from q, we can always choose a pair of arcs of class q', which start at C with directions transverse to the
velocity with which the arc of class q has arrived at C." This gives infinite collision chains with each gamma_k a
piece of a periodic 2-centre orbit of class q_k.

Appendix A (p. 22), positivity of entropy: P_n = number of periodic collision chains of period n; case (i) two
transverse periodic orbits through C (arcs are entire orbits): P_{2m+1} = 0 and P_{2m} = 2^{2m+1}; case (ii) one
orbit with a transverse autointersection at C (arcs are parts between passages): P_n = 2^{n+1}; in both
p = limsup (1/n) log P_n = log 2. "We conclude that the entropy corresponding to any choice of I in Q+ is
h_top >= log 2."

Reuse assessment (INFERRED): the construction is specific to the 2-centre problem (separability, elliptic
coordinates, Charlier's regularisation, rational resonance of the xi and phi periods). The transferable part is
Theorem 1 and 2 as a black box plus the checklist (finite set of nondegenerate collision arcs of equal energy,
direction change) and the Jacobian nondegeneracy test (11). The collision-arc construction itself does not
transfer to a Kepler-plus-moons or rotating-frame setting.

## 4. Explicit examples, numbers, figures (READ)

There are no tables and no numerical evaluation of the theorem constants. The only numbers are illustrative
parameter values for pictures of 2-centre periodic orbits:
- Figures 1 and 2 (p. 8): potential energy in xi and in phi for a = 1 and E = -0.5. (From the text: for xi the
  potential has local maximum at xi = 0 of value |E| - 2a, minima at cosh(+-xi_m) = a/|E| where it equals
  -a^2/|E|.)
- Figure 3 (p. 17): "Periodic orbits in Cartesian coordinates when a = 1, q = 1, beta = 1/7 and the orbit through
  the primaries": 18 periodic orbits, the two colliding with the primaries drawn in green, primaries as red
  asterisks, plus the orbit xi = xi_+ (different energy).
- Figure 4 (p. 17): "Periodic orbits through (xi_0, phi_0) = (2/3 xi_+, 0) in Cartesian coordinates, when a = 1,
  q = 1, beta = 1/7" (my reading of the fraction in the caption; the image is small).
- Figures 5 and 6 (p. 18): same point, a = 1, q = 2, beta = 1/7, and an enlargement showing autointersections
  (axis range of Figure 5 about +-80, so these orbits extend far from the primaries).
None of the figures includes the small centre's perturbed dynamics; they are all eps = 0 orbits.

## 5. Statements about extensions (READ)

The paper makes no claim about moving centres, the restricted three- or four-body problem treated by this method,
time-periodic problems, or several small masses as results of its own. The relevant passages, quoted:

- p. 2, on non-integrability literature (Bolotin 1985 [3]): "In particular, we have analytic non-integrability for
  the restricted circular many-body problem, in which a particle moves in a rotating plane, under the
  gravitational attraction of n centres fixed on this plane, when n > 2. Nevertheless, this generalisation does
  not add any additional information about the n-centre problem on a fixed plane, in which case it reduces to the
  result given in [2]." This is about "n centres fixed in a rotating plane", i.e. all centres at rest in one
  rotating frame, and is a non-integrability remark for positive energy, not a quasi-collision result.
- p. 3, on Bolotin and MacKay [4] applied to the circular restricted three-body problem: "The result of [4] has
  been applied by the two authors to the planar circular restricted problem of three bodies, in which a massless
  particle moves under the gravitational attraction of the other two bodies, the primaries, and the latters are
  supposed on a circular orbit about their centre of mass. The second primary is supposed to have small mass
  epsilon. The result is the existence of periodic and chaotic orbits which undergo consecutive close
  encounters with the smaller primary and which approach, as epsilon to 0, arcs of Kepler ellipses around the
  first primary, starting and ending at collision with the second."
- p. 3-4, on the alternative route [11], [12] (Font, Nunes and Simo): "A similar existence result has been
  obtained also in [11] by a completely different method ... first return map ... horseshoe like ... orbits with
  consecutive infinite close approaches with the small primary. A complete numerical study of these orbits is
  carried out in [12]."
- p. 3, on the Bolotin-MacKay general result: "valid for Lagrangian systems with Newtonian singularities on a
  d-dimensional Riemannian manifold, with d = 2, 3: they show that if the singular part of the Lagrangian is a
  perturbing term, then for each chain formed by arcs of the unperturbed problem which start and end at a
  singular point (collision arcs), there is a unique orbit with the same energy which shadows it."
  (READ: the phrase "a singular point" is singular in the printed text here, and Dimare's own definitions use one
  point C. Whether the original Bolotin-MacKay theorem allows chains through several singular points is
  documented in the 2000 digest, not here.)
- p. 2 literature on positive energy: [6] and [13] (spatial n-centre, topological entropy positive for n >= 3 and
  E >= 0; no analytic independent integral for E > E_th). Not about close passages.

Nothing else on extensions: no conclusion section, no list of open problems, no mention of Uranus, moons,
cyclers, rotating frames as a model for his result, or mass ratios other than epsilon of one centre.

## 6. What this does and does not cover for #890

(a) Orbits passing close to TWO different small centres in sequence. READ: no. The model has exactly one small
centre. All collision arcs begin and end at that one point C. The alphabet of the symbolic dynamics is a finite
set of arcs of that one centre (4 per class q). Any sequence of the form "close to moon A, then close to moon B"
is outside the printed result. INFERRED: the Theorem 1 mechanism (chains of collision arcs shadowed within order
epsilon) is stated for collision arcs ending at singular points of the perturbation; if its original statement
admits several singular points, a two-moon chain would be a different application of [4], not something Dimare
proves. His paper does not do it.

(b) Fixed centres and autonomy, and two moons with different periods. READ: the three centres are fixed in an
inertial plane and the problem is autonomous with conserved energy E; the energy level is the whole organising
structure (arcs of equal energy, E negative and tending to zero). Two moons with different orbital periods
cannot both be at rest in one rotating frame, and in an inertial frame they move; so a Titania-Oberon system is
not in this setting, and the paper offers nothing that gets round this: it never introduces a rotating frame,
time-dependent centres, a Jacobi-type integral, or a centre moving on a circle. The only rotating-plane
remark (p. 2) is Bolotin's non-integrability remark for n centres fixed on the rotating plane, which again has
all centres at rest in the rotating frame. A further mismatch (INFERRED from the setting): the unperturbed
problem here is the integrable Euler two-fixed-centre problem with two EQUAL large primaries at distance 2, and
the arcs are periodic 2-centre orbits through C with rational resonance q between the xi and phi oscillations;
the #890 unperturbed problem is a single large primary (Uranus), whose arcs are Kepler arcs.

(c) Wording. Supported by the printed text: "Dimare (2010, CMDA 107:427) proves, for the planar problem of a
massless particle in the field of two fixed equal primaries and one fixed small third centre of intensity
epsilon, at small negative energy and for a dense set of positions of the small centre, the existence of
uniformly hyperbolic invariant sets of periodic and chaotic orbits that shadow chains of collision arcs through
the small centre, with closest approach and shadowing error of order epsilon and Lyapunov exponents of order
log(1/epsilon), via the Bolotin-MacKay theorems"; and, as context, "the second-species mechanism has been shown
to work with non-Kepler unperturbed arcs (periodic orbits of the 2-centre problem)". Not supported, and should
not be attributed to this paper: any statement that periodic orbits with close flybys of two different small
bodies are proved or predicted by it; any statement about moving moons, a rotating frame, circular restricted
four-body or concentric circular models, or mass ratios of order 1e-5 to 1e-4 (no numerical bound on
epsilon_0, E_0, beta_0 is given); any statement that the #890 orbit is a member of a family it constructs.
INFERRED consequence for the literature check: the sentence "Dimare's 3-centre result points that way" in
`docs/notes/2026-10-04-890-literature-check.md` section 5 overstates this paper. It covers repeated close
passages to one small centre among fixed centres; it does not point at a two-moon alternating chain. The "first
predicted" question for #890 therefore still rests on whether Bolotin and MacKay 2000 itself covers several small
singular points in a time-dependent setting (see its digest), not on Dimare.

## 7. References cited by Dimare on second species and n-centre problems (p. 22)

Corpus status was checked by searching `docs/notes/CORPUS_INDEX.md` for the author names; "held" means a row
exists, "not held" means no row was found by that search.

1. [1] S. Aubry, R. S. MacKay, C. Baesens, "Equivalence of uniform hyperbolicity for symplectic twist maps and
   phonon gap for Frenkel-Kontorova models", Physica D 56 (1992), 123-134. Not held. (Cited for the equivalence
   of nondegeneracy and hyperbolicity.)
2. [2] S. V. Bolotin, "Nonintegrability of the n-center problem for n > 2", Moscow Univ. Mech. Bull. 39, No. 3,
   24-28 (1984); translation from Vestnik Mosk. Gos. Univ., Ser. I, math. mekh. 3 (1984), 65-68. Not held.
3. [3] S. V. Bolotin, "Influence of singularities of the potential energy on the integrability of dynamical
   systems", J. of Applied Math. and Mech. 48 (1985) No. 3, 255-260. Not held.
4. [4] S. V. Bolotin and R. S. MacKay, "Periodic and Chaotic Trajectories of the Second Species for the n-Centre
   Problem", Celest. Mech. & Dyn. Astr. 77 (2000), 49-75. HELD (digest
   `2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`).
5. [5] S. V. Bolotin and R. S. MacKay, "Non-planar second species periodic and chaotic trajectories for the
   circular restricted three-body problem", Celest. Mech. & Dyn. Astr. 94 (2006), No. 4, 433-449. HELD (digest
   `2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`).
6. [6] S. V. Bolotin and P. Negrini, "Regularization and topological entropy for the spatial n-center problem",
   Ergod. Th. & Dynam. Sys. 21 (2001), 383-399. Not held.
7. [7] S. V. Bolotin and P. Negrini, "Chaotic behavior in the 3-center problem", J. Differential Equations 190
   (2003), 539-558. Not held. (The far-away-third-centre, Melnikov-method predecessor.)
8. [8] C. L. Charlier, "Die Mechanik des Himmels", Bd. I, II, Verlag Von Veit & Comp., Leipzig, 1902. Not held.
9. [9] Y. Duan and J. Yuan, "Periodic orbits of the hydrogen molecular ion", Eur. Phys. J. D 6, 0 (1999),
   319-326. Not held.
10. [10] Y. Duan, J. Yuan and C. Bao, "Periodic orbits of the hydrogen molecular ion and their quantization",
    Phys. Rev. A 52, 5 (1995), 3497-3502. Not held.
11. [11] J. Font, A. Nunes, C. Simo, "Consecutive quasi collisions in the planar circular RTBP", Nonlinearity 15
    (2002), 115-142. HELD (index row `font-nunes-simo-2002-...`, digest `2026-10-04-digest-font-nun...`).
12. [12] J. Font, A. Nunes, C. Simo, "A numerical study of the orbits of second species of the planar circular
    RTBP", Celest. Mech. & Dyn. Astr. 103 (2009), 143-162. HELD (index row `font-nunes-simo-2009-...`).
13. [13] A. Knauf and I. A. Taimanov, "On the integrability of the n-centre problem", Math. Ann. 331 (2005),
    631-649. Not held.

The index does not list Marco and Niederman 1995 among Dimare's references; Dimare does not cite it (READ:
reference list above is complete).
