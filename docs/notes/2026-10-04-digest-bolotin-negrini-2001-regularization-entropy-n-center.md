# Digest: Bolotin and Negrini 2001, regularisation and topological entropy for the spatial n-centre problem

Date 2026-10-04. Purpose: establish exactly what this paper proves, because the project's #890 literature check
(`docs/notes/2026-10-04-890-literature-check.md`, section 5) said that "the Bolotin-Negrini snippet ... point[s]
that way" for a theory of the Titania-Oberon flyby orbit, and the coordinator's correction at the end of that note
withdrew the claim for the other second-species papers without the present paper having been read. It has now been
read in full. Headline: this is NOT a second-species or collision-chain paper. It proves that the n-centre
problem (n fixed Newtonian centres, any masses, energy not below zero) has positive topological entropy, by
regularising the collisions and then applying a topological theorem (Gromov, Paternain) to the geodesic flow of
the Maupertuis metric. It contains no small parameter, no collision arcs, no chains, no shadowing and no
direction-change condition. Its relation to the second-species theory of Bolotin and MacKay is one of shared
authorship and shared object (n fixed centres), not shared method. Related digests, not repeated here:
`docs/notes/2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`,
`docs/notes/2026-10-04-digest-bolotin-mackay-2006-nonplanar-second-species.md`,
`docs/notes/2026-10-04-digest-dimare-2010-quasi-collision-3-centre.md`,
`docs/notes/2026-10-04-digest-marco-niederman-1995-seconde-espece-plan-restreint.md`; the companion digest for
MacKay (2005) is `docs/notes/2026-10-04-digest-mackay-2005-chaos-in-three-physical-systems.md`.

Citation: S. V. Bolotin and P. Negrini, "Regularization and topological entropy for the spatial n-center
problem", Ergodic Theory and Dynamical Systems 21(2):383-399 (2001), DOI 10.1017/S0143385701001195. Received 16
July 1999, accepted in revised form 24 September 1999. 18 pages in the held file: PDF page 1 is the Cambridge
cover sheet and PDF page n is journal page n + 381 (so journal pages 383 to 399 are PDF pages 2 to 18). Filed in
the private paper corpus as
bolotin-negrini-2001-regularization-topological-entropy-spatial-n-center-problem-etds-21-383-doi-10.1017-S0143385701001195.pdf.
Text layer: present, about 9,800 words by `pdftotext | wc -w` (born digital; the formulas were read from the page
images, not from the text layer).

Evidence labels: READ means read from the printed page in this session. INFERRED means my reading of what the
printed text implies, or a comparison with project work. All 18 pages were read as page images. Page numbers
below are journal pages. Quotation marks mean verbatim from the page; mathematical notation is transcribed into
plain text and "not equal" is written in words.

## 0. Summary (READ)

Abstract (p. 383): "We show that the n-center problem in R^3 has positive topological entropy for n >= 3. The
proof is based on global regularization of singularities and the results of Gromov and Paternain on the
topological entropy of geodesic flows. The n-center problem in S^3 is also studied."

Structure (p. 386, "The paper is organized as follows"): Sections 2 to 4, global regularisation of Newtonian
singularities in a three-dimensional manifold (a global version of the Kustaanheimo-Stiefel regularisation);
Section 5, the topology of the regularised configuration space; Section 6, topological entropy of the regularised
flow and the S^3 theorem; Section 7, the non-compact case R^3 and the proofs of Theorems 1.1 and 1.2.

## 1. The setting and the main theorems (pp. 383-386, READ)

### 1.1 The system (p. 383)

"Let P = {p_1, ..., p_n} be a finite set in R^3. The n-center problem is a Lagrangian system

  q'' = -grad V(q), q in R^3 \ P,   (1.1)

with the Lagrangian

  L(q, q') = (1/2)|q'|^2 - V(q), V(q) = - sum_{i=1}^n mu_i / |q - p_i| + phi(q), mu_i > 0,   (1.2)

where phi in C^infinity(R^3). Hence the potential energy V has n Newtonian singularities. We assume that V < 0
everywhere and there exists R > 0 such that

  V(q) + (1/2)<q, grad V(q)> < 0, for |q| >= R.   (1.3)

For the classical n-center problem phi = 0 and condition (1.3) is obviously satisfied. The classical two-center
problem is a completely integrable system as shown by Euler. For Lagrange's two-center problem, phi =
-|q|^2/2."

Energy: "Let H(q, q') = (1/2)|q'|^2 + V(q) be the energy integral. Fix an energy level Sigma_E = {(q, q') in
(R^3 \ P) x R^3 | H(q, q') = E}, E >= 0." Consequence (1.4), "By Lagrange's equality and condition (1.3) for any
solution of equation (1.1) of energy E >= 0, d^2/dt^2 |q|^2 = 4(E - V - (1/2)<q, grad V(q)>) > 0, for |q| >= R.
Hence all trajectories on Sigma_E exiting the R-ball B escape to infinity, and non-trivial dynamics occurs in
Sigma_E intersect {q in B} only."

Answers to the setting questions:
- Number of centres: any finite n; the theorems need n >= 3 (R^3) or n >= 5 (S^3).
- Fixed or moving: fixed. P is a fixed finite subset of R^3 and V is a function of q alone.
- Masses: mu_i is any positive number. No mass is small; there is no small parameter in the paper. The
  theorems hold for every choice of the mu_i.
- Energy range: E >= 0 only (p. 384 Remark 3 below says this is necessary). With V < 0 everywhere, E >= 0 means
  the particle is unbound relative to the centres: the dynamics studied is scattering-type motion through the
  region around the centres, not bound motion.
- Time dependence: none. L has no time argument, H is conserved, and the whole method works on the fixed level
  Sigma_E.
- Lagrangian form: only kinetic minus potential. No velocity-linear (Coriolis or magnetic) term appears anywhere
  in the paper; I looked for one on every page. (Compare the Bolotin-MacKay Lagrangian (1.2) in the 2000 paper,
  which has a <omega(q), q'> term and is what makes the rotating-frame restricted three-body problem fit.)
  Remark 4 below allows a conformally Euclidean kinetic metric but not a gyroscopic term.
- Dimension: R^3 (spatial) or S^3, and in Section 2 a general connected three-dimensional Riemannian manifold Q.
  The planar R^2 problem is not proved here; it is cited as known (Remark 1 below).

### 1.2 Theorem 1.1 and the definition of entropy (pp. 384)

"THEOREM 1.1. For n >= 3 and E >= 0 the topological entropy h_top of the n-center problem on Sigma_E is
positive."

How entropy is defined for an incomplete flow (p. 384): "Note that no natural definition of the topological
entropy for incomplete flows with singularities seems to exist." The trajectories of energy E are mapped by the
Maupertuis principle to geodesics of g_E(q, q') = 2(E - V(q))|q'|^2 on R^3 \ P, with the time change ds =
sqrt(2(E - V)) dt. The geodesic flow g^t on Sigma = {g_E = 1} is not complete because of collisions, "However, it
is well known that collisions can be regularized, so that g^t is extended to a smooth flow g~^t : Sigma~ ->
Sigma~ without singularities on a manifold Sigma~ containing Sigma. Topologically, Sigma is obtained from
Sigma~ by removing two-dimensional spheres S_i^2 corresponding to collisions with p_i. The trajectories of the
flow g~^t cross S_i^2 transversally." The definition adopted: "We define the topological entropy of the geodesic
flow g^t as the topological entropy of the regularized flow h_top(Sigma, g^t) = h_top(Sigma~, g~^t)." Sigma~ is
non-compact, but by (1.4) the non-wandering set Lambda lies in q in B, so the entropy is that of the restriction
to the compact invariant set Lambda, and "the positiveness of the topological entropy is invariant under strictly
monotone time reparametrization", hence independent of which regularisation is used.

Remarks (p. 384, READ):
1. "It was proved in [2, 8] that the n-center problem in R^2 is non-integrable for n >= 3, and has positive
   topological entropy. If phi = 0 and all the centers in R^3 lie in one plane, then the n-center problem in R^3
   contains the n-center problem in R^2 as a subsystem and, hence, this also has positive topological entropy.
   However, this argument does not hold if the centers do not lie in a plane." ([2] is Bolotin 1984; [8] is
   Klein and Knauf 1992.)
2. "Under the conditions of Theorem 1.1 the system (1.1) can be partly integrable. For example, when phi = 0 and
   the centers lie on a line, there is an integral of angular momentum."
3. "The condition E >= 0 is necessary. For any n it is possible to construct an n-center problem satisfying all
   conditions of Theorem 1.1 that is integrable on Sigma_E for some E < 0." (p. 384 to 385: for example V such
   that for large negative E the set {V <= E} is a union of balls U_i around p_i with V = -1/|q - p_i| there.)
4. "Theorem 1.1 holds also if the Lagrangian is replaced by L = T(q, q') - V(q), where the Riemannian metric T
   (kinetic energy) is conformally Euclidean near p_i and tends to the Euclidean metric at infinity, so that
   Lagrange's inequality holds for large |q|."

### 1.3 Theorems 1.2 and 1.3, Proposition 1.1 (pp. 385-386)

"THEOREM 1.2. For n >= 3 and E >= 0 there exists c > 0 such that N_lambda(p, q) >= e^{c lambda} for almost all
p, q in B \ P." Here N_lambda(p, q) is "the number of such trajectories gamma : [a, b] -> B \ P with the
Maupertuis action integral_gamma sqrt(g_E) dt = integral_a^b |gamma'(t)|^2 dt <= lambda" connecting p and q
(trajectories of energy E; by Lagrange's inequality any such trajectory lies in B). "Thus the number of
trajectories connecting two generic points grows exponentially with the action. Hence the geodesic entropy h_geod
of the metric g_E on B \ P is positive." Geodesic entropy is defined by (1.5): h_geod = limsup_{lambda -> infinity}
(1/lambda) log of the integral over M x M of N_lambda(p, q) with respect to the volume measure. The paper notes
that Mane and Paternain-Paternain proved h_top = h_geod on closed manifolds, that no such equality holds in
general for non-compact manifolds, and that it will show h_top >= h_geod for M = B \ P "due to a special
behavior of the metric at infinity". The constant c is not given a value anywhere.

"THEOREM 1.3. For n >= 5 and E >= 0 the topological entropy of the system on Sigma_E and the geodesic entropy of
the Maupertuis metric g_E on S^3 \ P are both positive." (the n-centre problem on the sphere S^3, with a
Riemannian metric on S^3 conformally Euclidean near P.) "There exist counter examples for any n <= 4."

"PROPOSITION 1.1. For n <= 4 there exist a C^infinity Riemannian metric T on S^3 conformally Euclidean in a
neighborhood of P and a C^infinity potential V < 0 having Newtonian singularities at P such that the n-center
problem on S^3 with L = T(q, q') - V(q) is integrable on the energy level Sigma_E for some energy E > 0, and the
topological entropy on Sigma_E is zero." The proof of Theorem 1.3 is written out only for even n; "The case of
odd n is a little different and will be treated elsewhere" (p. 386).

Remark (p. 386): "Similar results hold for an n-center problem in an arbitrary two-dimensional manifold Q [3, 7,
8]. Then the topological entropy is always positive for n greater than twice the Euler characteristics of Q. We
are unable to treat an n-center problem in an arbitrary three-dimensional manifold Q, because we are using the
methods of Gromov and Paternain, and so have to assume that Q is simply connected. Similar results hold if the
fundamental group of Q is large enough. For example, using the theorem of Taimanov [16] it is easy to show that
the n-center problem on a closed manifold Q is non-integrable if dim H_1(Q, R) >= 5."

## 2. The method (pp. 386-398, READ)

### 2.1 Global regularisation, Theorem 2.1 (pp. 386-388)

Setting: Q a connected three-dimensional Riemannian manifold, kinetic energy T conformally Euclidean near each
centre (2.1) T = a_i(q)|q'|^2/2, potential V bounded above with Newtonian singularities V = -f_i(q)/|q| (2.2),
V < 0, H = T + V, and the Maupertuis-Jacobi metric g_E = 4(E - V) T (2.3). Two cases: Q non-compact or with
boundary; or Q closed and n even (for closed Q and n odd "the global regularization of singularities is more
complicated and will be discussed in another paper").

Theorem 2.1 (p. 387), in outline: for n even there exist a four-dimensional manifold M, a smooth circle action
Phi_t on M, a smooth surjective map phi : M -> Q and a Riemannian metric G on M such that phi is invariant under
Phi_t; the preimage of each centre p_i is a single point s_i; the action is free on M minus S and phi restricted
to M minus S is a fibre bundle over Q minus P with structure group S^1; G is Phi_t-invariant; and the norm of
phi'(x) x' in g_E equals the norm of x' in G for every x' orthogonal (in G) to the Killing field v generating Phi_t.
Corollary 2.1: the momentum integral F(x, x') = <v(x), x'>_G of the geodesic flow of G, restricted to the zero
level F = 0, maps geodesics of G on M minus S into geodesics of g_E on Q minus P.

The regularised energy level (p. 388): "the geodesic flow of the metric g_E on Q minus P is obtained from the
geodesic flow of the metric G on M by means of a reduction with respect to the symmetry group Phi_t [1]. The
geodesics of G passing through one of the points s_i in S are projected to geodesics of g_E that collide with p_i
and are reflected in the opposite direction." Sigma~ is the quotient Gamma / Phi'_t of the zero level of the
momentum integral, a smooth manifold; Sigma~ minus Sigma is a union of n two-spheres (the unit tangent spheres at
the s_i, Hopf-fibred). READ: this sentence is the entire treatment of what happens at a collision in this paper:
the regularised flow reflects the collision orbit back along itself.

Remarks (p. 388): for Q two-dimensional there is "an analog of Theorem 2.1 ... a global version of the
Levi-Civita regularization [3]", with phi : M -> Q a two-sheeted ramified covering. "There are several other ways
to regularize double collisions. The most famous is Moser's regularization [12]. ... However, other known
regularizations do not preserve the Lagrangian structure of the system (for example, Moser's regularization
interchanges coordinates and momenta), and hence are not suitable for our purposes."

### 2.2 The local model: the Kustaanheimo-Stiefel (KS) regularisation (pp. 388-389)

For one centre in R^3 the paper reformulates KS geometrically. The Hopf map h : R^4 -> R^3 is quadratic with |h(x)|
= |x|^2, h(e^{Jt} x) = h(x) with J skew-symmetric and J^2 = -1, and if <Jx, x'> = 0 then |h'(x) x'| = 2|x||x'|.
In quaternions h(x) = x i x-bar and Jx = x i. The metric G(x, x') = u(h(x))|x'|^2, u(q) = 8(E|q| + f(q))a(q)
(3.2), has the isometry group x -> e^{Jt} x with first integral F = u(x)<x', Jx>; the zero level of F projects to
geodesics of g_E. G has no singularity at x = 0, so the geodesic flow of G on {F = 0} regularises that of g_E.
The connection form of the Hopf fibration, lambda_x = <Jx, dx>/|x|^2 (3.3), and the curvature (monopole) form
Omega = (q x dq) wedge dq / (2|q|^3) (3.4) with the integral of Omega over S^2 equal to 2 pi, are what the global
construction patches together.

### 2.3 Topology of the regularised space and the entropy argument (Sections 4 to 7)

- Section 4 (pp. 390-391): the circle bundle over B minus P is built from its Chern class, using closed 2-forms
  omega = sum (-1)^i Omega_i (4.1) with alternating signs at the centres. For n odd in a closed manifold the signs
  cannot sum to zero, which is why n must be even (Stokes: 0 = 2 pi sum of plus or minus 1).
- Section 5 (pp. 392-394): the regularised manifold M_n is computed. For one centre in R^3, M = R^4. Proposition
  5.1: "For n >= 1 the manifold M_{n+2} is a connected sum of M_n and S^2 x S^2: M_{n+2} = M_n # (S^2 x S^2)",
  by surgery. Consequences, (5.1) to (5.3): for Q = S^3 and k >= 2, M_{2k} is a connected sum of k - 1 copies of
  S^2 x S^2; "For n >= 1 centers in R^3, M_n is homotopically equivalent to S^2 v ... v S^2 (n - 1 times)."
  Corollary 5.1: for n >= 2 even centres in S^3, H_2(M, Z) = Z^{n-2}; for n >= 1 centres in R^3, H_2(M, Z) =
  Z^{n-1}. Proposition 1.1 (integrable, zero entropy for n <= 4 on S^3) is proved here by exhibiting S^2 x S^2
  with the product of two rotation actions as the regularised space of a four-centre problem.
- Section 6 (pp. 394-396): topological entropy of the regularised flow is bounded below by the exponential growth
  of the number n_N(x, lambda) of geodesics of length at most lambda leaving a fibre N = phi^{-1}(q) orthogonally and
  ending at x (inequality (6.5), Gromov and Paternain); Lemma 6.1 bounds the sum of Betti numbers of the path
  space Omega(N, x) by n_N(x, Ck); equation (6.8), b_i(Omega(N, x)) = b_i(Omega(M)) + b_{i-1}(Omega(M)). A simply
  connected M is rationally hyperbolic if the sum of Betti numbers of its loop space grows exponentially.
  Corollary 6.1: "If M is rationally hyperbolic, the topological entropy h_top(Sigma, g^t) is positive." For
  n = 2k centres in S^3, "dim(pi_2(M) tensor Q) = dim H_2(M, Q) = 2(k - 1)" and by (6.9) M is rationally
  hyperbolic for k >= 3, that is for even n >= 6 (the proof is written for even n only, p. 386; n = 4 is the
  integrable counterexample of Proposition 1.1). Lemma 6.2: for almost all q, p the connecting geodesics in M project to trajectories that avoid
  P, so N_lambda(q, p) = n_N(x, lambda).
- Section 7 (pp. 397-398), the non-compact case Q = R^3: Lagrange's inequality makes the ball B geodesically
  convex in the Maupertuis metric (Lemma 7.1 builds G so that the preimage W = phi^{-1}(B) is geodesically convex
  in G as well). The Gromov-Paternain machinery needs closed manifolds and Yomdin's theorem; the authors replace
  the flow by one with compact support, extend it to a closed manifold by the identity, and prove (7.2), "h_top(N
  perp, G~) >= limsup (1/lambda) log integral_W n_N(x, lambda) d mu(x)". Theorem 1.2 and then Theorem 1.1 follow
  from (6.6), (7.1) and (7.2).

The method is thus purely topological and global: positive entropy comes from the second homology of the
regularised configuration space, which has n - 1 independent two-spheres for n centres in R^3. Nothing in it
selects, constructs or enumerates individual orbits. The only "arcs" in the paper are the geodesics counted by
n_N(x, lambda), and they are shown to avoid collisions for almost every pair of endpoints (Lemma 6.2).

Reference list (p. 399): Arnold-Kozlov-Neishtadt; Bolotin 1984 (two papers); Friedlander-Halperin; Gromov 1986;
Katok-Hasselblatt; Knauf 1987; Klein-Knauf 1992; Kobayashi-Nomizu; Mane 1997; Milnor-Stasheff; Moser 1970;
Paternain 1992 and 1997; Paternain-Paternain 1994; Taimanov 1988; Kustaanheimo-Stiefel 1965; Yomdin 1987. It
cites neither Poincare nor any second-species, collision-orbit or shadowing literature, and it does not cite
Bolotin and MacKay (2000), which appeared after this paper was accepted.

## 3. Statements on extensions (READ, exhaustive)

The paper makes no statement about: moving centres; the restricted three-body or four-body problem; two small
bodies; time-periodic, quasi-periodic or otherwise time-dependent potentials; negative energy (other than Remark
3, which says that E >= 0 is necessary because integrable examples exist for some E < 0); gyroscopic (Coriolis)
terms; shadowing of collision chains; a small parameter. I read every page for these. The only generalisations it
states are the ones quoted above: a kinetic metric conformally Euclidean near the centres (Remark 4, p. 385);
other two-dimensional manifolds Q (p. 386 Remark; Levi-Civita analogue p. 388); other three-dimensional
manifolds when simply connected or with large enough fundamental group (p. 386 Remark); and the sphere S^3.

## 4. What it does and does not cover for the two-moon object (printed hypotheses only)

The #890 model is a planar periodic orbit of a massless particle moving in the field of Uranus with Titania and
Oberon on concentric circular orbits of different periods, alternating close flybys of the two moons, bound to
Uranus (negative Kepler energy relative to Uranus), found by continuation in the moon mass.

Against the printed hypotheses of Theorem 1.1 (all READ; the conclusion that the theorem does not apply is
INFERRED from them):
1. Centres fixed: fails. The moons move, with different periods, so no frame makes both stationary.
2. Autonomous, conserved energy: fails in any frame (the potential depends on time through two incommensurate or
   commensurate moon phases). The method works on the level set Sigma_E of a time-independent Hamiltonian.
3. No gyroscopic term: the rotating frame of one moon introduces a Coriolis term not in the Lagrangian (1.2); the
   inertial frame has no Coriolis term but the moons are then not fixed.
4. Energy E >= 0 with V < 0 everywhere: the #890 orbit is bound to Uranus. Remark 3 (p. 384) says the theorem is
   false in general for E < 0.
5. Dimension: planar. The theorem is for R^3, and R^2 is cited from other work (Remark 1).
6. Smallness: not required, which is a difference in the other direction: the theorem holds for the physical
   moon masses if its other hypotheses hold. But it gives no information about any individual orbit.

Conclusion (INFERRED from READ): this paper proves nothing about the #890 object and does not "point that way".
What it does show, at the level of mathematical background: for n >= 3 fixed Newtonian centres in space, at
non-negative energy and for any masses, there is chaos in the sense of positive topological entropy, and the
number of collision-free connecting trajectories grows exponentially with the action. Anyone wishing to cite it
for the project should cite it for those two facts and no others. Safe wording (INFERRED): "Bolotin and Negrini
(2001) proved that the problem of n >= 3 fixed Newtonian centres in R^3 at energy not below zero has positive
topological entropy, for all positive masses; it is a statement about fixed centres at non-negative energy and
does not cover moving centres, bound orbits, or the existence of any specific flyby orbit."

## 5. Techniques applicable to the project's problems

This paper supplies little that is directly constructive. The constructive material (collision arcs, chains,
admissibility, shadowing) is in the 2000 and 2006 Bolotin-MacKay papers and in MacKay (2005); the full
enumeration recipe is written out in section 5(a) of the MacKay digest and only summarised here, to avoid two
copies drifting apart. What this paper adds is listed first.

What this paper adds (READ for the facts, INFERRED for the uses):
- Collision regularisation as a numerical device. The paper's regularisations are KS in space and Levi-Civita in
  the plane (p. 388 Remark 1), chosen because they preserve the Lagrangian structure (p. 388 Remark 2). A
  regularised moon-centred integration through a close approach lets a numerical search continue through, or very
  near, an encounter without step collapse, and lets a continuation in the moon's mass go to values where the
  periapsis is of order epsilon times the arc size. The fixed-centre theorem does not transfer to a moving moon,
  but the local change of variables is a smooth, time-rescaled substitution and does not itself need the moon to
  be fixed (INFERRED; the project would have to check that a time-dependent moon position does not spoil the
  regularised form).
- The statement that at a regularised collision the flow "reflects" the orbit back along itself (p. 388) is the
  reason the Bolotin-MacKay direction-change condition excludes the reversal case: a chain edge whose arrival
  velocity is the exact opposite of the departure velocity is a regularised bounce, not a flyby with a
  collision-avoiding periapsis (INFERRED, connecting this sentence to the 2000 paper's Remark 3).
- Counting. Theorem 1.2 says the number of collision-free connecting trajectories grows like e^{c lambda} in the
  Maupertuis action, and Lemma 6.2 says the collision orbits are a measure-zero nuisance. For a search this
  means the number of candidate chains of bounded action grows exponentially with the allowed action (or
  time), so any enumeration must be indexed by an action or time cap and by a pruning rule; c is not given.

(a) Graph-of-collision-arcs enumeration: not in this paper. See section 5(a) of the MacKay digest for the
construction, the three admissibility tests (equal energy or Jacobi constant, nondegeneracy, direction change in
the relative-velocity frame) and how to test each for Kepler arcs between moon encounters.

(b) Direction change versus the demanded-turn gate (`src/cyclerfinder/verify/turn_gate.py`): this paper has no
direction-change condition and no small-mass scaling, so it says nothing about the gate. The only relevant fact
is the one above, that exact reversal of direction is a regularised collision (a periapsis of zero), which the
gate registers as a demanded turn of 180 degrees with a required altitude of minus the body radius, that is,
infeasible. See section 5(b) of the MacKay digest.

(c) Hyperbolicity, Lyapunov exponents, navigation cost, conditioning: this paper gives positive topological
entropy and an existence-only exponential count with an unquantified constant c. It states no Lyapunov exponent,
no hyperbolicity, and no dependence on mass. For navigation cost it therefore gives nothing quantitative. The
quantitative instability statement (exponent of order log of one over epsilon, per unit of the secondary's
period) is MacKay's (2005, p. 2) and the 2006 paper's Theorem 2.2; see the MacKay digest section 5(c). One
qualitative inference from this paper (INFERRED): because the entropy is global and topological, it survives
changes of the metric and of the potential that keep the singularity type (Remarks 3 and 4 and "the positiveness
... is independent of the construction of the regularized flow"), so some chaotic structure is expected to be
robust to model changes; this is a statement about existence, not about the size of any exponent.

(d) Extension to two moons with different periods: this paper's closest relevant statement is negative. Its
hypotheses (fixed centres, autonomous, E >= 0) exclude the object, and its method (a global topological
argument on the regularised configuration space of a time-independent system) has no obvious time-dependent
version, because with time dependence there is no energy level set to carry a geodesic flow. (INFERRED.) A
related point: two fixed centres alone (n = 2) is integrable (Euler; here also M = S^2 x R^2 has loop-space
Betti numbers of linear growth, so it is rationally elliptic), and chaos starts at n = 3 in R^3. The moons are
perturbations of the Uranus Kepler problem, not extra fixed centres of comparable strength, so the n = 3 theorem
is not the relevant regime even if the moons were frozen. What would have to be proved or computed is listed in
section 5(d) of the MacKay digest and is not repeated here.

(e) Frame conventions as a warning: this paper has no frames issue because its centres are fixed in the inertial
frame, which is itself the warning. The theorem is stated for fixed centres in one fixed frame; as soon as the
centres move, the frame in which "fixed" and "direction" are judged must be stated, and the 2000 paper's error
was exactly a direction test made in the wrong frame. The check of which frame the project's turn gate uses is in
section 5(e) of the MacKay digest: it measures the turn between velocities relative to the flyby body.
