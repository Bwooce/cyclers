# Digest: Bolotin and MacKay 2006, nonplanar second species periodic and chaotic trajectories

Date 2026-10-04. Purpose: establish exactly how far the published "second species" theory covers the
#890 object (a periodic orbit of the planar concentric circular restricted four-body problem, Uranus plus
Titania plus Oberon on circular orbits with different periods, one close flyby of each moon per cycle,
found by continuing a patched-conic flyby chain in the moon mass; see
`docs/notes/2026-10-04-890-titania-oberon-candidate.md` and
`docs/notes/2026-10-04-890-literature-check.md`).

Citation: S. Bolotin and R. S. MacKay, "Nonplanar second species periodic and chaotic trajectories for the
circular restricted three-body problem", Celestial Mechanics and Dynamical Astronomy 94:433-449 (2006),
DOI 10.1007/s10569-006-9006-0. 17 pages (journal pages 433 to 449; PDF page n is journal page n+432).
Filed in the private paper corpus as
bolotin-mackay-2006-nonplanar-second-species-periodic-chaotic-trajectories-cr3bp-cmda-94-433-doi-10.1007-s10569-006-9006-0.pdf.

Evidence labels: READ means read from the printed page in this session. INFERRED means my reading of what
the printed text implies. The whole paper (all 17 pages) was read as page images. Formulas are transcribed
from the images; anything marked [unclear] could not be read with confidence. Page numbers are journal pages.

## 0. One-paragraph summary (READ)

Abstract (p. 433): "For the circular restricted three-body problem of celestial mechanics with small
secondary mass, we prove the existence of uniformly hyperbolic invariant sets of non-planar periodic and
chaotic almost collision orbits. Poincare conjectured existence of periodic ones and gave them the name
'second species solutions'. We obtain large subshifts of finite type containing solutions of this type."
The paper has two layers. Section 2 states an abstract theorem (Theorem 2.1, proved in Bolotin and MacKay
2000, not here) for a time-independent Lagrangian system with a finite set of fixed singular points, a
conserved energy, and a small parameter epsilon multiplying the singular potential. Section 3 applies it to
the spatial circular restricted three-body problem with ONE small mass, in the frame rotating with that mass,
using nonplanar Kepler ellipses as collision arcs. Appendix A proves the hyperbolicity statement
(Theorem 2.2), the only new abstract result.

## 1. Section 2, general setting (pp. 434-436, READ)

### 1.1 The Lagrangian system (p. 434)

"Let P = {p_1, ..., p_n} be a finite set in a 3D manifold Q. Consider a Lagrangian system (L_eps) with
configuration space Q \ P and Lagrangian

  L_eps(q, qdot) = L_0(q, qdot) - eps V(q).   (2.1)

We assume that L_0 is C^4 everywhere in Q and quadratic in the velocity:

  L_0(q, qdot) = T(q, qdot) + <w(q) . qdot> - W(q),   (2.2)

where the kinetic energy T(q, qdot) is a positive definite quadratic form in qdot, and w(q) is a covector
field on Q. Let V be a C^4 function on Q \ P having Newtonian singularities on P. This means that in a
neighborhood U_alpha of any point p_alpha in P,

  V(q) = - f_alpha(q) / dist(q, p_alpha),   f_alpha(p_alpha) > 0,   (2.3)

where f_alpha is a C^4 function on U_alpha, and the distance is defined by means of the Riemannian metric T.
We study system (L_eps) for small eps > 0. Then it is a singular perturbation of system (L_0)."

Energy (p. 434): "Let H_eps = H_0 + eps V, H_0(q, qdot) = T(q, qdot) + W(q) (2.4) be the energy integral. We fix
E such that the domain D = {q in Q | W(q) < E} contains the set P and study system (L_eps) on the energy
level {H_eps = E}."

Answers to the specific questions on the setting:
- Number of centres: a finite set, any number n. The text calls the problem "the general n-centre problem" in
  the introduction (p. 434), referring to Bolotin and MacKay 2000 (title: "Periodic and chaotic trajectories
  of the second species for the n-centre problem").
- Are the centres fixed points of configuration space: yes. P is a fixed finite subset of the manifold Q, and
  the Lagrangian is a function of (q, qdot) only. Nothing in Section 2 lets a centre move. (READ for the
  setting as printed; the absence of any mechanism for moving centres is INFERRED from there being no time
  argument anywhere in (2.1) to (2.4).)
- Autonomous with conserved energy: yes. L_eps has no explicit time dependence, H_eps is called "the energy
  integral", and the whole analysis is on a fixed level set {H_eps = E}. The theorem statement itself
  requires "the same energy E" for all arcs (see below).
- Same epsilon for all centres: V is a single function carrying all the singularities, multiplied by one
  epsilon. All centres therefore scale together; the positive functions f_alpha may differ, so the
  relative strengths of the centres are not tied to be equal. (READ for the form (2.3); the remark about
  relative strengths is INFERRED from f_alpha being an arbitrary positive C^4 function.)
- In the application of Section 3 the CR3BP has exactly one centre: "The singular set consists of one point
  P." (p. 436), and the Sun's singularity is not in V but in W, with Q = R^3 \ {0} so that L_0 is C^4 on Q.

### 1.2 Definitions (p. 435)

Collision arc: "We say that a solution gamma: [0, tau] -> D of system (L_0) is a collision arc if gamma(0),
gamma(tau) in P and has
- No early collisions: gamma(t) not in P for 0 < t < tau."

Energy, Jacobi functional: "Let E be the energy of gamma. Then gamma is a critical point of the
Maupertuis-Jacobi functional (see e.g. Arnold et al. 1989) J_E on the set Omega of nonparameterized curves
in D with end points in P:

  J_E(gamma) = integral_0^tau g_E(gamma(t), gammadot(t)) dt,
  g_E(q, qdot) = 2 sqrt((E - W(q)) T(q, qdot)) + <w(q) . qdot>,

where g_E is the Jacobi metric. Since W|_D < E, the functional J_E is well defined on Omega. We say that the
collision arc gamma is
- Nondegenerate if it is a nondegenerate critical point of J_E."

Practical nondegeneracy (p. 435): "Represent the general solution of system (L_0) as q(t) = f(q_0, v_0, t),
where q_0, v_0 in R^3 are initial position and velocity. Then collision arcs with energy E connecting p_alpha
to p_beta correspond to solutions of the system of four equations

  f(p_alpha, v, tau) = p_beta,   H_0(p_alpha, v) = E   (2.5)

in four variables v, tau. The nondegeneracy condition is that the Jacobian at the solution is nonzero."

Collision chain: "Suppose that system (L_0) has nondegenerate collision arcs gamma_k: [0, tau_k] -> D, k in K
(a finite set) with the same energy E connecting the points p_alpha_k and p_beta_k. A sequence (gamma_{k_i})
for i in Z is called a collision chain (Poincare called them 'orbites a chocs') if beta_{k_i} = alpha_{k_{i+1}}
and satisfies
- Direction change: gammadot_{k_i}(tau_{k_i}) is not equal to +- gammadot_{k_{i+1}}(0) for all i."

So consecutive arcs must meet at the same centre (the arc ends at the point where the next begins), and the
outgoing velocity of the next arc must be neither parallel nor antiparallel to the incoming velocity of the
previous one. Graph (p. 435): "Collision chains correspond to paths in the graph Gamma with the set of
vertices K and the set of edges Gamma = {(k, k') in K^2 | beta_k = alpha_k', gammadot_k(tau_k) not equal to
+- gammadot_k'(0)}. (2.6)"

### 1.3 Theorems (pp. 435-436), quoted in full

"Theorem 2.1. Given a finite set K of nondegenerate collision arcs with the same energy E, there exists
eps_0 > 0 such that for all eps in (0, eps_0] and any collision chain (gamma_{k_i}), k_i in K, there exists a
unique (up to a time shift) trajectory gamma: R -> D \ P of energy E of system (L_eps), which shadows the
chain (gamma_{k_i}) within order eps. More precisely, there exist c, C > 0, independent of eps and the
collision chain, and a sequence (t_i) such that |t_{i+1} - t_i - tau_{k_i}| <= C eps, dist(gamma(t),
gamma_{k_i}(t - t_i)) <= C eps for t_i <= t <= t_{i+1}, and dist(gamma(t), P) >= c eps."

Immediately after (p. 435): "Hence there is an invariant subset Lambda_eps in {H_eps = E} on which (L_eps) is
a suspension of a subshift of finite type. It can be proved to be uniformly hyperbolic, and strongly so."
Attribution: "The following result is proved in Bolotin and MacKay (2000)." (So Theorem 2.1 is imported, not
proved here.)

"Theorem 2.2. There exists a cross-section N in {H_eps = E} such that the corresponding invariant set
M_eps = Lambda_eps intersect N of the Poincare map is uniformly hyperbolic with Lyapunov exponents of order
log eps^{-1}."

"Corollary 2.1. The set Lambda_eps is uniformly hyperbolic, as a suspension of a hyperbolic invariant set with
bounded transition times."

Following text (p. 436): "Theorem 2.2 can be deduced from the proof in Bolotin and MacKay (2000) of
Theorem 2.1, but it was not proved there. Thus we prove it here in Appendix A.

The topological entropy of the Poincare map on M_eps is positive provided the graph Gamma has a connected
branched sub-graph. In fact the topological entropy is O(eps)-close to that of the topological Markov chain
determined by the graph Gamma. In the case of a periodic sequence (k_i), i in Z, local uniqueness of the
trajectory gamma implies that it is also periodic, closing after one cycle of the sequence."

"Remarks. One can allow the nonsingular part L_0 of the Lagrangian L_eps also to depend on eps. Then all the
results remain true with L_0 replaced by L_0|_{eps=0}.

The result can be extended to some L_eps, which are not quadratic in the velocity; a case like this arises
for the reduction of the motion of two charges in a uniform magnetic field with respect to Euclidean
symmetry."

Note that the dependence of L_0 on eps is dependence on the small parameter, not on time. No remark in the
paper allows explicit time dependence. (READ for what is printed; the statement that none is allowed is a
statement about the absence of such a remark.)

### 1.4 The one quantitative ingredient hidden in the theorem (READ)

Appendix A, p. 448, eq. (5.17) region: the Maupertuis action of the connecting orbit near a centre is
J_E = S_alpha(a, b, eps) = S_+(a) + S_-(b) + eps s_alpha(a, b, eps) - c_alpha eps log eps, with
c_alpha = f_alpha(p_alpha) from (2.3), and s_alpha(a, b, 0) = c_alpha log || v_+(a) x v_-(b) || (5.17).
The proof that connecting orbits have time intervals of order log eps^{-1} is cited to Bolotin and MacKay
(2000) (p. 448 to 449). This is the origin of the log eps^{-1} Lyapunov exponent.

## 2. Section 3, the application (pp. 436-437, READ)

Setting (p. 436): "Consider the spatial circular restricted three-body problem (Sun, Jupiter, and Asteroid,
with the Sun and Jupiter moving in circles around their center of mass and the Asteroid of zero mass free to
move in 3D) and we suppose that the mass of Jupiter is small with respect to the mass of the Sun. We
normalize the masses to 1 - eps (Sun), eps (Jupiter), and 0 (Asteroid), with the center of mass stationary
and the first two masses in circular orbits about it, having separation and angular frequency both normalized
to 1."

Frame: "the motion of the Asteroid in the frame Oxyz rotating anti-clockwise about the z-axis through the Sun
at angular frequency 1 with respect to the Sun. Then the Sun is at O = (0, 0, 0), and Jupiter can be chosen
at P = (1, 0, 0)." The Lagrangian is of the form (2.1) with

  L_0(q, qdot) = (1/2)|qdot|^2 + x ydot - y xdot + W(q),   W(q) = (1/2)|q|^2 + 1/|q|   (3.7)
  V(q) = 1/|q| - 1/|q - P| + x   (3.8)

"Hence L_0 has the form (2.2), V has the form (2.3), and Q = R^3 \ {0}. The singular set consists of one
point P." The energy integral in the rotating frame is the Jacobi integral
H_eps = (1/2)|qdot|^2 - (1/2)|q|^2 - (1 - eps)/|q| - eps/|q - P| + eps x, with value -C/2, and
C = -2 E_fixed + 2 G_z (E_fixed the energy in the fixed frame, G_z the z-component of angular momentum about O).

So what plays the role of the centre P is the single small secondary, which is stationary in the frame
rotating with it. The Sun's 1/|q| is moved into the unperturbed part W (the base problem L_0 is the Kepler
problem of the Sun in the rotating frame, p. 437). Note that the Sun is excluded from Q rather than treated as
a collision centre, and the arcs are Kepler ellipses about the Sun.

Jacobi constant range (p. 437): at eps = 0 the orbits are Kepler ellipses, with angular frequency
Omega = a^{-3/2} and Jacobi constant C = a^{-1} + 2 sqrt(a(1 - e^2)) cos iota, iota the inclination to
Jupiter's orbital plane. "Given C in R we define the set A_C of allowed frequencies of Kepler ellipses to be
(0, 1) if C in [-1, +2]; (0, (2 + C)^{3/2}) if C in (-2, -1); ((3 - C)^{3/2}, 1) if C in (2, 3); and empty if
C not in (-2, +3)."

"Theorem 3.1. For any C in (-2, +3) there exists a dense subset S of the set A_C of allowed frequencies, such
that for any finite set Lambda subset S there exists eps_0 > 0 such that for any sequence (Omega_n), n in Z,
in Lambda and eps in (0, eps_0) there is a trajectory of the spatial circular restricted three-body problem
with Jacobi constant C, which avoids collisions by order eps and in the rotating frame is within order eps of
a concatenation of collision orbits formed from arcs of Kepler ellipses of frequencies Omega_n with
inclinations iota_n satisfying cos iota_n = C/2 - Omega_n^{2/3}. The resulting invariant set is uniformly
hyperbolic.

In particular, the Poincare map for given Jacobi constant has a chaotic invariant set with Lyapunov exponents
of order log eps^{-1}, and it contains infinitely many nonplanar periodic second species orbits."

The inclination relation follows from C = Omega^{2/3} + 2 cos iota (p. 443, Lemma 4.1).

## 3. Section 4, construction of the collision arcs (pp. 437-444, READ)

### 3.1 Which arcs (pp. 437-438)

Four classes of Kepler arc between two intersections with a circle (the unit circle, Jupiter's orbit),
following Poincare 1899 (p. 437): "1. a whole number of revolutions of a coplanar orbit; 2. a segment of
coplanar orbit between distinct intersection points; 3. a whole number of revolutions of a noncoplanar orbit;
4. a segment of a noncoplanar orbit between points at opposite ends of a straight line through the Sun."
(p. 438) "Here we consider only the last class of orbit, because construction of subshifts from the first
class was already done in Bolotin and MacKay (2000), construction from the third class looks problematic to
us, and the orbit of Marco and Niederman (1995) was generated from two arcs in the second class so it is not
as virgin territory as the fourth class (though still merits treating one day). An interesting question that
we will address at the end of the paper is whether one can make subshifts using arcs from a combination of
classes."

So this paper's arcs are class 4: ellipses with inclination iota in (0, pi) (planar orbits excluded), crossing
Jupiter's orbit at the ends of a diameter. Collision arc duration (4.12): tau = 2 pi k + pi, with Jupiter's
period 2 pi. Condition for crossing at +-j: a(1 - e^2) = 1 (4.11), with a = sin^{-2} eta, Omega = sin^3 eta,
C = sin^2 eta + 2 cos iota (Section 4.2, p. 440). The collision arcs are parametrised by integers (k, m, l)
through eq. (4.13), pi(2k+1) sin^3 eta - pi(2m+1) + sigma g(eta) = 0, g(eta) = pi - 2 eta + sin 2 eta, and
Lemma 4.1 (p. 443) gives a dense set Lambda_sigma in (0, 1) of attainable frequencies.

The important point for #890: the collision arcs are heliocentric (Sun-centred) Kepler ellipses that start and
end at the single small body's instantaneous position, which is fixed in the rotating frame. The "collision"
is a collision with Jupiter. Both "ends" of a collision arc are at Jupiter.

### 3.2 Nondegeneracy (Section 4.4, pp. 443-444)

A collision arc is nondegenerate if the derivative of (D, tau, C) with respect to the initial velocity v is
invertible, where D is the distance from the origin at which the nearby trajectory re-pierces the horizontal
plane, tau the transit time and C the Jacobi constant. Replacing v by the position F of the second focus, the
three moves (F around a circle perpendicular to Pi; F radially in Pi; F parallel to the intersection of Pi
with the horizontal plane) change C, tau and D respectively at nonzero rate, giving a triangular derivative
with nonzero diagonal: "Thus the derivative of (D, tau, C) with respect to F is triangular with nonzero
diagonal entries, so invertible." The key nonvanishing quantity is dD/d alpha = -+ 2e not equal to 0 at
alpha = +-pi/2 (p. 444). Section 4.2 separately proves uniqueness and nondegeneracy of the solutions eta of
(4.13); the case l = 0 (eta = pi/2) is nondegenerate but ignored, and l = k (eta = 0) is degenerate and
excluded (p. 441).

### 3.3 Direction change (Section 4.5, p. 444) and footnote 3

"We make collision chains by connecting collision arcs with the same Jacobi constant, but we must be sure
that they satisfy the 'direction change' condition. This requires the velocity in the rotating frame just
after each collision to be neither parallel nor opposite to the velocity just before the collision.^3"

Footnote 3 (p. 444): "At the analogous point in Bolotin and MacKay (2000) we mistakenly studied the
direction change in the inertial frame; this error was corrected in MacKay (2005)."

The relative velocity to Jupiter is v = v_z e_z + (v_phi - 1) e_phi + v_rho e_rho in the rotating frame with
v_phi = G_z = cos iota. The condition fails only if cos iota = cos iota' or cos iota + cos iota' = 2; the
second is impossible for nonplanar orbits. For Omega not equal to Omega' the first is excluded by
conservation of C. For Omega = Omega' the paper switches the sign lambda of the perihelion side (the
mirror-image ellipse) so the condition always holds, because the orbits are nonplanar: "So by choosing to
switch sign of lambda we satisfy the changing direction condition." The authors conclude the proof of
Theorem 3.1.

Relevance (INFERRED): the direction-change condition is formulated in the frame in which the small body is
stationary. For one small body on a circular orbit that frame is unambiguous. With two moons on circles of
different periods there is no frame in which both are stationary, so the condition has no frame-independent
meaning unless one chooses one moon's frame and treats the other as moving.

## 4. Quantitative statements (READ)

- Closest approach: Theorem 2.1, "dist(gamma(t), P) >= c eps", with c independent of eps and of the chain.
  Theorem 3.1: "avoids collisions by order eps". The p. 445 comments say the extended-mass replacement of
  Jupiter works for a sphere "of radius c eps" and that deviation from spherical symmetry must decay "to much
  less than 1/c^2 eps at this radius" (the printed fragment is "1/c^2 eps"; transcription of the exponent
  and placement of eps is [unclear] at that point of the page).
- Shadowing: within order eps in position, within C eps in transit times (Theorem 2.1).
- Lyapunov exponents: order log eps^{-1} (Theorem 2.2, Theorem 3.1, p. 445 comment: "highly unstable, with
  Lyapunov exponents of order log eps^{-1}. This strong instability implies strong controllability, a key
  fact long recognized by the designers of solar system exploration missions using flyby.")
- Topological entropy: O(eps)-close to that of the topological Markov chain of the graph Gamma, positive
  if Gamma has a connected branched sub-graph (p. 436).
- Size of eps_0: not given anywhere. Theorem 2.1 and Theorem 3.1 are pure existence statements, "there exists
  eps_0 > 0", with c and C unspecified. The paper contains no numerical value of eps_0, no mass ratio for any
  real system, and no numerical example. This was checked across all pages.
- Jacobi constant range: C in (-2, +3), with the allowed frequency sets A_C (Section 2 of this digest).

## 5. Comments and open questions (Section 5, p. 444-445), quoted

"For small enough ratio eps of secondary to primary mass in the circular restricted three-body problem and any
value of Jacobi constant in the range for which there exist nonplanar Kepler ellipses crossing the unit circle
twice, we have proved existence of arbitrarily large uniformly hyperbolic subshifts of finite type consisting
of non-planar second species orbits. This result is a 3D analog of the planar result proved in Bolotin and
MacKay (2000)."

Proved extensions:
"As in the planar case, the result remains true if the Sun is replaced by an extended mass distribution
provided it is constant in the rotating frame, because the only effect is to make a small change to the
potential W. Similarly, Jupiter can be replaced by a spherically symmetric mass distribution confined to a
sphere of radius c eps, because it produces the same field as a point mass at its center and the constructed
orbits avoid collision by at least c eps. In fact, one can also replace Jupiter by any nonspherically
symmetric mass distribution provided it is constant in the rotating frame, contained within a radius less
than c eps about its center of mass, and the effect on the gravitational field of deviation from spherical
symmetry has decayed to much less than 1/c^2 eps at this radius. This allows Jupiter a significant oblateness
(J_2 component) for example."

"We have also proved in Appendix A a general result which implies that both in the planar and nonplanar
cases, the resulting second species orbits are highly unstable, with Lyapunov exponents of order log eps^{-1}."

Said to be probable ("Probably", "Presumably"; not proved):
- "Probably we could also make subshifts of finite type using some parabolic and hyperbolic Kepler arcs in
  addition to the elliptical ones used here."
- "Probably we could make subshifts using infinitely many collision arcs (also in the planar case), by
  restricting to sequences for which the direction change is bounded away from 0 and pi, but this would need
  more careful control on the nondegeneracy, and uniform hyperbolicity for the flow (though perhaps not the
  map) would be lost because the durations of the collision arcs would be unbounded. Probably, we could make
  unbounded orbits too, ... for any C in (-2 sqrt 2, +2 sqrt 2) for the planar case."
- "Presumably we could extend the result of Bolotin and MacKay (2000) to make planar subshifts using class 2
  arcs ..., like Marco and Niederman's orbit. Presumably we could combine them with the class 1 arcs used in
  Bolotin and MacKay (2000), to make even bigger planar subshifts."

Open questions: "An interesting question is whether one could construct subshifts based on sequences of both
planar and nonplanar collision arcs. Their existence does not follow from our analysis because although the
planar arcs used in Bolotin and MacKay (2000) are nondegenerate with respect to variations in the plane, they
are degenerate with respect to 3D variations."

Difficulty: "We recall from Bolotin and MacKay (2000), however, a problem with using class 3 arcs, namely that
they are all degenerate. So more delicate analysis would be required to make trajectories to shadow sequences
of them. Existence is not impossible, but might require an analog of the method of Bolotin (2006)."

What the closing comments do NOT mention (READ, by absence): more than one small body, two secondaries with
different periods, a time-dependent (elliptic or otherwise) restricted problem as a proved extension, the full
three-body problem (small massive third body). Section 1 mentions only that the elliptic extension has been
made in other work (Bolotin 2005, 2006), and Poincare's original claim concerned "second and third masses m,
mu small compared to the primary mass M", that is, two small masses in the full problem; this paper restricts
to mu = 0.

## 6. The history paragraph (pp. 433-434, READ)

Quoted from the introduction: "In chapter XXXII of Poincare 1899, he proposed that there are periodic
solutions of the three-body problem of celestial mechanics with second and third masses m, mu small compared
to the primary mass M, which as m, mu -> 0 converge to pairs of segments of Kepler orbit joined at
collisions. He derived several necessary conditions on the sequences of collision arcs which occur as the
limits and sketched an argument that these are sufficient for existence of nearby second species orbits when
m, mu are small enough. He christened them 'second species' orbits.

It is agreed (Levy 1952), however, that Poincare did not provide a proof and that the result is not true in
the full generality that he claimed. Despite many analyses (e.g., Alexeev 1970; Henrard 1980; Bruno 1981;
Perko 1981; Gomez and Olle 1991), it is only recently that any complete proofs have been written, and so far
they are only for the 'restricted' case where mu = 0. Marco and Niederman (1995) proved the existence of a
periodic second species orbit with two collisions per period for the planar circular restricted case. For the
same case, we proved in Bolotin and MacKay (2000) the existence of large sets of second species orbits
[footnote 1: See also MacKay (2005) for a summary, some minor additions and a correction], including
aperiodic analogues, to which we proposed to extend the same name. They form uniformly hyperbolic subshifts.
A similar result was subsequently obtained by Font et al. (2002) by a different method but it is limited to
orbits with small angle changes at collisions. In contrast, the angle changes are large in Bolotin and MacKay
(2000) and in the present paper. In recent work (Bolotin 2005, 2006), an extension has been made to the
slightly elliptic case.

In the present paper, we extend our analysis of the circular restricted three-body problem to prove existence
of uniformly hyperbolic subshifts of nonplanar second species orbits. The method is the same as in Bolotin
and MacKay (2000). Indeed we already prepared the ground there by allowing for the 3D case in our analysis of
the general n-center problem and in remarking that it was likely that existence of one of Poincare's classes
of nonplanar second species orbits could be proved by this method."

Who proved what, from that paragraph:
- Poincare 1899: proposed (conjectured) the existence, in the full problem with two small masses; named them
  second species; sketched an argument; no proof.
- Levy 1952 (the editorial notes in Poincare's Oeuvres): agreed that it was not proved and not true in the
  generality claimed.
- Marco and Niederman 1995: one periodic orbit with two collisions per period, planar circular restricted.
- Bolotin and MacKay 2000: large sets, including aperiodic, uniformly hyperbolic subshifts, planar circular
  restricted; large angle changes; n-centre general setting.
- Font, Nunes and Simo 2002: similar result by a different method, small angle changes only (the full
  title is "Consecutive quasi-collisions in the planar circular RTBP").
- Bolotin 2005 and 2006: extension to the slightly elliptic case (from the text; I have not read those
  papers, so what exactly each proves is not established here).
- MacKay 2005: summary of Bolotin and MacKay 2000, minor additions, a correction (the frame in which the
  direction change is measured, footnote 3).

Full citations, from the reference list (p. 449, READ):
- Alexeyev, V.M.: Sur l'allure finale du mouvement dans le probleme des trois corps. Actes du Congres Int.
  Math. 2, 893-907 (1970).
- Arnold, V.I., Kozlov, V.V., Neishtadt, A.I.: Mathematical Aspects of Classical and Celestial Mechanics,
  Encyclopedia of Mathematical Sciences, vol. 3, Springer-Verlag, Berlin (1989).
- Aubry, S., MacKay, R.S., Baesens, C.: Equivalence of uniform hyperbolicity for symplectic twist maps and
  phonon gap for Frenkel-Kontorova models. Physica D 56, 123-134 (1992).
- Bolotin, S.: Second species periodic orbits of the elliptic 3-body problem. Celest. Mech. Dyn. Astron. 93,
  345-373 (2005).
- Bolotin, S.: Shadowing chains of collision orbits. Discr. Conts. Dyn. Syst. 14, 235-260 (2006).
- Bolotin, S.V., MacKay, R.S.: Periodic and chaotic trajectories of the second species for the n-centre
  problem. Celest. Mech. Dyn. Astron. 77, 49-75 (2000).
- Bruno, A.D.: On periodic flybys to the Moon. Celest. Mech. 24, 255-268 (1981).
- Font, J., Nunes, A., Simo, C.: Consecutive quasi-collisions in the planar circular RTBP. Nonlinearity 15,
  115-142 (2002).
- Gomez, G., Olle, M.: Second species solutions in the circular and elliptic restricted three body problem,
  I and II. Celest. Mech. Dyn. Astron. 52, 107-146, 147-166 (1991).
- Henon, M.R.: Generating Families in the Restricted Three Body Problem, Lect. Notes in Phys. Monographs, 52,
  Springer, Berlin (1997).
- Henrard, J.: On Poincare second species solutions. Celest. Mech. 21, 83-97 (1980).
- Levy, P.: Editorial notes. In: Poincare, H. Oeuvres Tome VII, p. 629, Gauthiers-Villars, Paris (1952).
- MacKay, R.S.: Chaos in three physical systems. In: Dumortier, F., Broer, H., Mawhin, J., Vanderbauwhede,
  A., Verduyn Lunel S. (eds.) Equadiff 2003, pp. 59-72, World Scientific, Singapore (2005).
- Marco, J.-P., Niederman, L.: Sur la construction des solutions de seconde espece dans le probleme plan
  restreint des trois corps, Ann. Inst. H. Poincare Phys. Theor. 62, 211-249 (1995).
- Perko, L.M.: Second species solutions with an O(mu^nu), 0 < nu < 1 near-Moon passage. Celest. Mech. 24,
  155-171 (1981).
- Poincare, H.: Les Methodes Nouvelles de la Mecanique Celeste, Tome III, Gauthiers-Villars, Paris (1899).
- Xia, Z.: Arnold diffusion and oscillatory solutions in the planar three-body problem. J. Diff. Eq. 110,
  289-321 (1994).

## 7. What this does and does not cover for #890

### 7.1 Hypotheses of Theorem 2.1 against the #890 setting

The #890 object lives in a planar restricted four-body model: Uranus (primary), Titania and Oberon on
concentric circular orbits with different periods (the model of the project note, with the particle of zero
mass). In any rotating frame at most one moon is stationary, so the problem is non-autonomous (time periodic
if the two periods are commensurate, quasi-periodic otherwise) and has no energy integral or Jacobi integral.

Hypotheses of Theorem 2.1 as printed (READ), and their status for #890:
1. The singular set P is a finite set of points of a fixed manifold Q. FAILS for #890 in any single frame:
   one moon can sit at a fixed point, the other moves on a circle in that frame.
2. The Lagrangian L_eps = L_0 - eps V is time-independent, so H_eps = E is a conserved energy and all
   collision arcs must share "the same energy E". FAILS: with two moons of different periods there is no
   conserved Jacobi-type energy, so "the same energy E" is not defined and the Maupertuis-Jacobi variational
   characterisation (the definition of nondegenerate arc through J_E) does not apply as written. (A
   restricted conservation in the single-moon limit is not a property of the two-moon model.)
3. V = sum of Newtonian singularities, all multiplied by one eps, with f_alpha positive C^4. HOLDS in
   spirit if both moons are scaled together; a mass ratio between the moons is allowed through the f_alpha
   only if the f_alpha are fixed functions on Q, which again presupposes fixed centres.
4. L_0 is C^4, quadratic in velocity, and may depend on eps (Remarks p. 436). The remark does not allow
   dependence on time. Uranus plus the particle with one fixed frame would be a legitimate L_0, but the
   moving second moon cannot be put into W (a time-independent potential), and it cannot be treated as a
   perturbation because it is a singular centre of the same strength as the first.
5. Collision arcs are nondegenerate critical points of J_E with the same energy; a chain requires
   beta_{k_i} = alpha_{k_{i+1}}, i.e. each arc ends at the centre where the next begins. For a chain of the
   #890 type the arcs would alternate between two centres (Titania then Oberon then Titania). That is
   permissible in the n-centre setting as stated (P may have n > 1 points), but only when the centres are
   fixed, and the arcs are travelled at fixed energy with the transit time determined by the arc, not by the
   moons' phasing. In #890 the Oberon encounter is an event in time, set by the phase of Oberon; that
   phase condition has no counterpart in Theorem 2.1.
6. The direction-change condition, Section 4.5, is formulated in the frame in which the single centre is
   stationary (footnote 3). It has no frame-independent version for two moving centres.

Verdict on (a), from hypotheses 1, 2 and 4: Theorem 2.1 as printed does not apply to one large and TWO small
masses on circular orbits with different periods. The failing hypotheses are the fixed finite singular set P
and the autonomous, energy-conserving form of L_eps. The paper never claims such an extension.

A further scope point (READ): Theorem 2.1 is stated for a 3D manifold Q. The #890 model is planar. The
planar single-small-mass case is Bolotin and MacKay 2000, not this paper; this paper's Theorem 3.1 requires
nonplanar arcs (iota in (0, pi)) and explicitly says planar arcs are degenerate with respect to 3D
variations. A planar statement would rest on the 2000 paper (not read here), so even the one-moon planar
analogue is, for the project, a statement taken from a paper not yet held.

### 7.2 Does the paper cite or mention any result for two or more small secondaries (b)

No (READ). The only places that touch more than one small body are Poincare's original two-small-masses
proposal (m and mu both small, p. 433), which the authors report as unproven in full generality and "only
recently" proved in the restricted case mu = 0, and the "n-centre problem" of the 2000 paper (fixed centres,
not moving moons). Neither is a result for two orbiting secondaries with different periods. The reference
list contains nothing on a restricted four-body problem, a bicircular model, or any problem with two
secondaries, and Section 5 states no extension in that direction, even as a probable one. The reference
list does carry Gomez and Olle 1991 (circular and elliptic restricted problem, which is still one
secondary).

### 7.3 Wording the paper supports and does not support for #890 (c)

Supported (each is a statement about this paper, not about #890):
- "For the circular restricted three-body problem with a small secondary, Bolotin and MacKay proved, in the
  planar case (2000) and the spatial case (2006), that for sufficiently small mass ratio there are
  uniformly hyperbolic families of periodic and chaotic orbits that shadow chains of collision arcs, with
  closest approach of order eps and Lyapunov exponents of order log(1/eps)." (Theorem 3.1 for the spatial
  case; the planar case is attributed to the 2000 paper.)
- "The #890 object is the two-moon analogue of such orbits. The one-moon version is a theorem; a two-moon
  version with different periods is not covered by Theorem 2.1 as stated, whose hypotheses (fixed singular
  points, time-independent Lagrangian, conserved energy) it does not satisfy."
- "The theorem is an existence statement for sufficiently small eps with no stated eps_0, so it does not say
  whether the physical Titania and Oberon masses are within its range, even for one moon."
- The mass-continuation heuristic is of the same character as the theory: both are statements about the
  singular small-mass limit in which the arc between flybys is a Kepler arc about the primary. That is
  consistent with section 2.9 of the #890 note, but the paper does not prove or test the project's
  continuation in eps from 0.01 to 1. (INFERRED.)

Not supported:
- Any statement that the #890 orbit "exists by Bolotin and MacKay", "is a second species orbit proved to
  exist", or "is covered by the published theorem". The project's orbit was computed, not obtained from the
  theorem.
- Any statement that the paper proves the persistence of a two-moon chain at the physical masses, or that it
  gives a range of validity (no numerical eps_0, no real-system mass ratio).
- Any statement that a two-secondary version is "known", "expected", or "left open by the authors"; the
  paper mentions none of these. The authors' probable extensions are parabolic and hyperbolic arcs,
  infinitely many arcs and unbounded orbits, class 2 arcs, and combined planar and nonplanar subshifts.
- Any claim of novelty or non-novelty of the #890 object resting on this paper alone. It bears on the
  THEORY class (spec 16.4 (ii) "system"/class language in the literature-check note, section 5), not on the
  body set.
- Any assumption that the elliptic extension (Bolotin 2005, 2006) covers a two-moon problem. Those papers
  treat the slightly elliptic single-secondary problem, which is time dependent, so they are the nearest
  published non-autonomous analogue; what exactly they require has not been read (INFERRED from the title
  and the sentence on p. 434).

## 8. Related papers the project should obtain

None of these is in the corpus index or held (checked by grep of `docs/notes/CORPUS_INDEX.md` and the
directory listing, 2026-10-04).

1. Bolotin, S.V., MacKay, R.S.: Periodic and chaotic trajectories of the second species for the n-centre
   problem. Celest. Mech. Dyn. Astron. 77, 49-75 (2000). Source of Theorem 2.1, the planar case, the n-centre
   setting; needed to state the planar analogue and the exact hypotheses.
2. Bolotin, S.: Second species periodic orbits of the elliptic 3-body problem. Celest. Mech. Dyn. Astron. 93,
   345-373 (2005). The time-dependent (slightly elliptic) extension; the nearest published non-autonomous
   case.
3. Bolotin, S.: Shadowing chains of collision orbits. Discr. Conts. Dyn. Syst. 14, 235-260 (2006). The
   method that handles degenerate arcs and the elliptic case.
4. MacKay, R.S.: Chaos in three physical systems. In: Dumortier, F., Broer, H., Mawhin, J., Vanderbauwhede,
   A., Verduyn Lunel S. (eds.) Equadiff 2003, pp. 59-72, World Scientific, Singapore (2005). Summary with
   the frame-of-direction-change correction.
5. Font, J., Nunes, A., Simo, C.: Consecutive quasi-collisions in the planar circular RTBP. Nonlinearity 15,
   115-142 (2002). Independent method, small-angle encounters, and likely numerical content on the size of
   the mass ratio.
6. Marco, J.-P., Niederman, L.: Sur la construction des solutions de seconde espece dans le probleme plan
   restreint des trois corps, Ann. Inst. H. Poincare Phys. Theor. 62, 211-249 (1995).
7. Gomez, G., Olle, M.: Second species solutions in the circular and elliptic restricted three body problem,
   I and II. Celest. Mech. Dyn. Astron. 52, 107-146, 147-166 (1991). Numerical continuation of second
   species orbits; closest to the project's computational approach.
8. Poincare, H.: Les Methodes Nouvelles de la Mecanique Celeste, Tome III, Gauthiers-Villars, Paris (1899),
   chapter XXXII (second species) and Levy's editorial note, Oeuvres VII (1952), p. 629.
9. Aubry, S., MacKay, R.S., Baesens, C.: Equivalence of uniform hyperbolicity for symplectic twist maps and
   phonon gap for Frenkel-Kontorova models. Physica D 56, 123-134 (1992). Only if the hyperbolicity proof
   itself is ever needed.

Not in this paper's reference list but worth a search for the two-secondary question (not cited here, so
no citation is claimed): work on second species or collision-chain orbits in the bicircular or restricted
four-body problem and in the full three-body problem with two small masses.


## Coordinator's addendum (2026-10-04, later the same day)

This digest says MacKay (2005) is not held. It is now held and digested:
`docs/notes/2026-10-04-digest-mackay-2005-chaos-in-three-physical-systems.md`. Its footnote a
(p. 3) is the correction referred to here: the direction change at a collision must be measured
in the rotating frame (with velocities relative to the secondary); it "holds in the rotating
frame for all such sequences except those containing two consecutive Omega_n satisfying" a
printed exceptional relation, and MacKay states that the planar theorem stands. The exceptional
relation as printed could not be rederived by the digest's author and is recorded as unresolved
(`#910`).
