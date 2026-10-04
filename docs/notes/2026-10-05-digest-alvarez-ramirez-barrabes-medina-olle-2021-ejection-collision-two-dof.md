# Digest: Alvarez-Ramirez, Barrabes, Medina and Olle 2021, ejection-collision orbits in two-degree-of-freedom problems

Date: 2026-10-05 (Sydney). Reading, reasoning and consistency checks against the 2019 paper; no project code was changed.

Source: M. Alvarez-Ramirez, E. Barrabes, M. Medina and M. Olle, "Ejection-Collision Orbits in Two Degrees of Freedom
Problems in Celestial Mechanics", J. Nonlinear Sci. 31:68 (2021), DOI 10.1007/s00332-021-09721-5, received 14 August 2020,
accepted 13 May 2021 (communicated by M. Leok). Filed in the private paper corpus as
`alvarez-ramirez-barrabes-medina-olle-2021-ejection-collision-orbits-two-degrees-of-freedom-jns-31-68-doi-10.1007-s00332-021-09721-5.pdf`.
A 33-page PDF with a good text layer. I read the whole text; the equations below were checked against each other
(Section 3 lists the two internal consistency checks I ran). I did not view page images, so the arrow-diagram figures in
Theorems 2 to 5 are given from the text layer only (the arrows are lost there and are described in words).

Companion: the same group's 2019 numerical paper, digested in
`2026-10-05-digest-alvarez-ramirez-barrabes-medina-olle-2019-ejection-collision-four-body.md`. This 2021 paper is its
purely analytic generalisation. It prints no numerical result of any kind (no table, no orbit, no tolerance).

## 1. Verdict in brief

- Scope: a general Hamiltonian with two degrees of freedom, a diagonal constant kinetic matrix and a potential homogeneous of
  degree -1 with exactly two partial-collision rays. It covers the collinear three-body problem, the rectangular,
  rhomboidal and symmetric collinear four-body problems, and a symmetric planar 2N-body problem.
- The planar restricted three-body problem (rotating frame) and the Hill problem are NOT covered (Section 3). The paper does not
  mention either, nor any restricted problem, nor the Olle, Rodriguez and Soler papers (it cites none of them; the reference
  list is in Section 9).
- It gives no numerical method (the 2019 paper does). What it adds is a theorem-level classification of which collision
  sequences are guaranteed to occur, in terms of the shape of the one-dimensional manifolds on the collision manifold.
- Use for the project: a template for how ejection-collision orbits are organised (a Poincare map on the partial-collision
  section and backwards iteration of arcs of the stable manifold), and a counting argument whose first cases the 2019 numbers
  reproduce exactly. Not a source of a test for a planar-pass propagator.

## 2. The general setting

System (eq. 1) and Hamiltonian (eq. 2):

    q' = A^(-1) p,   p' = grad U(q),   H(q, p) = (1/2) p^T A^(-1) p - U(q),

with q in R^2 minus a singular set Delta, p in R^2, A = diag(a1, a2) with a1, a2 > 0 constant, and U homogeneous of degree -1
on R^2 minus Delta (singular on Delta, which stands for all partial collisions; q = 0 is total collision). Energy fixed at
H = h < 0. Proposition 1: bounded motion needs h < 0 (from the Lagrange-Jacobi relation I'' = U + 2h with I = (1/2) q^T A q).

Hypotheses on the potential, stated on the shape function V(theta) = r U(q) (Proposition 2):

    V(theta) = c_b/sin(theta_b - theta) + c_a/sin(theta - theta_a) + Vtilde(theta),   theta in (theta_a, theta_b),

with 0 < theta_b - theta_a <= pi; c_a > 0 and c_b >= 0 constants, with c_b = 0 if and only if theta_b - theta_a = pi;
Vtilde > 0 smooth and bounded on [theta_a, theta_b]; V has exactly one critical point theta_c in (theta_a, theta_b), non-degenerate,
a minimum. Under these, the regularised system has two hyperbolic saddle equilibria E+- (below), with (on a fixed energy level)
dim W^u(E-) = 1, dim W^s(E-) = 2, dim W^u(E+) = 2, dim W^s(E+) = 1. The proof is deferred to Martinez (2012). The transversality
of the stable and unstable manifolds of the total-collision and total-ejection homothetic orbits follows from the non-degenerate
critical point, by Simo and Llibre (1981).

Worked examples given (V in theta):
- Rectangular four-body, theta in (0, pi/2): V = 2 + 2/cos(theta) + 2/sin(theta). Fits with c_a = c_b = 2 and Vtilde = 2.
- Rhomboidal four-body, theta in (0, pi/2): V = 1/(sqrt(2) cos(theta)) + alpha^(5/2)/(sqrt(2) sin(theta)) + 4 sqrt(2) alpha^(3/2)/(alpha cos^2(theta) + sin^2(theta)).
- Symmetric collinear four-body, theta in [theta_alpha, pi/2], theta_alpha = arctan(sqrt(alpha)): the same V as eq. 4 of the 2019 paper.
  The 2021 text prints it without the 2019 factor layout but with the same four terms.
- Collinear three-body (McGehee): V(xi) with xi in (-1, 1), singularities at xi = -+1, lambda a mass-dependent constant.
- Symmetric planar 2N-body problem with 2N equal masses (Martinez 2013): V = 1/sin(pi/N - theta) + 1/sin(pi/N + theta) + Vtilde,
  theta in [-pi/N, pi/N]. The printed text writes the singular terms with an index "n" and the section range with N; I read it as
  one parameter.

## 3. Does the restricted problem satisfy the hypotheses?

No, on four independent grounds (my analysis; the paper says nothing on the point):

1. The Hamiltonian must be H = (1/2) p^T A^(-1) p - U(q) with U homogeneous of degree -1. The planar circular restricted problem in
   a rotating frame has H = (1/2)|p|^2 + (p_x y - p_y x) (Coriolis, linear in p) - (1 - mu)/r1 - mu/r2 - correction terms; the gyroscopic
   term is outside the class and the effective potential is not homogeneous. Hill's problem adds the tidal terms with the same defect.
2. Total collision (q = 0 with all bodies coincident) is not a feature of the restricted problem. Collision with the small
   primary is a BINARY collision with the other bodies far away; there is no polar shape variable theta at that point with two
   distinguished partial-collision rays.
3. In the two-body (Kepler) limit relevant to collision with a primary, U = 1/|q| gives V(theta) = 1, constant. There is no
   unique non-degenerate minimum (every angle is critical), and no c_a/sin singular terms. So the paper's equilibrium E+- is
   degenerate there: a Kepler radial ejection or fall exists at every angle (a circle of homothetic orbits), which is the
   familiar reason angular momentum has to be zero for a collision.
4. In the limit mu = 0 the problem is Kepler in a rotating frame; the collision orbits there are the radial Keplerian ones, with
   no hyperbolic saddle structure of this type.

What the restricted problem does offer instead (not from this paper): a one-dimensional family of ejection orbits from the
small primary, parametrised by the ejection angle at a given Jacobi constant, with Levi-Civita regularisation. The Olle,
Rodriguez and Soler papers listed in the `#899` entry of `data/OUTSTANDING.md` (CNSNS 55:298; 90:105294; 94:105550; 111:106410) are the
restricted-problem ones; they are not held, and this paper does not cite them. The 2021 paper should therefore not be quoted as
covering the restricted problem. If the restricted-problem papers define their regularisation as a McGehee-type blow-up of collision
with the primary, as the OUTSTANDING note suggests for the 2022 spatial one ("spatial, McGehee"), their framework differs from
this one: the shape manifold there is not an interval (theta_a, theta_b) with two partial collisions.

## 4. The McGehee-type coordinates, written out

Step 1, polar-like blow-up. With r^2 = q^T A q and q = r w, w = (A^(-1))^(1/2) (cos(theta), sin(theta))^T (so w^T A w = 1),
the transverse unit is u = (A^(-1))^(1/2) (-sin(theta), cos(theta))^T (u^T A u = 1, w^T A u = 0); radial speed r' = w^T p;
q' = r' w + r theta' u.

Step 2, scaled velocities and time. v = r^(1/2) r', u = r^(3/2) theta' (the symbol u is reused for the unit vector and for this
coordinate in the printed text), p = r^(-1/2) A (v w + u u), new time d tau = r^(-3/2) d tau-tilde.

System (eq. 3), r = 0 invariant:

    dr/dtau = r v
    dv/dtau = v^2/2 + u^2 - V(theta)
    dtheta/dtau = u
    du/dtau = -v u + V'(theta)/2

Energy relation (eq. 4): h r = (v^2 + u^2)/2 - V(theta). Zero-velocity curve (eq. 5): V(theta) + h r = 0.
Equilibria E+- = (r, v, theta, u) = (0, +-v_c, theta_c, 0) with v_c^2 = 2 V(theta_c) (note the factor 2, which differs from the
2019 paper's v_c^2 = V(theta_c); the 2019 paper's variables carry a different scaling of p and of the energy, so the two
statements describe different normalisations; I did not reconcile them). Total collision manifold (eq. 6):
C = {r = 0, theta_a < theta < theta_b, u^2 + v^2 = 2 V(theta)}. On C: dv/dtau = h r + u^2/2 = u^2/2 >= 0 (gradient-like in v).

Step 3, removal of the two partial collisions simultaneously (Sundman type). With f(theta) = sin(theta - theta_a) sin(theta_b - theta)
if theta_b - theta_a != pi, else f = sin(theta_b - theta); W(theta) = f(theta) V(theta); F(theta) = f/sqrt(W);
w = F(theta) u; d tau = F(theta) dt. Regularised system (eq. 7):

    dr/dt = r v F(theta)
    dv/dt = F(theta) (2 h r - v^2/2) + sqrt(W(theta))
    dtheta/dt = w
    dw/dt = -F(theta) v w/2 + (W'(theta)/W(theta)) (f(theta) - w^2/2) + f'(theta) (1 + (f(theta)/W(theta)) (2 h r - v^2))

Energy relation (eq. 8): W(theta) w^2 + f(theta)^2 v^2 = 2 f(theta)^2 h r + 2 W(theta) f(theta).
Symmetry (eq. 9, 10): (r, v, theta, w, t) -> (r, -v, theta, -w, -t). Scaling: (lambda r, v, theta, w) with energy lambda h, so one
negative h suffices. Phase space F = [0, infinity) x R x [theta_a, theta_b] x R.

Relation to the 2019 eqs. 5 to 7 (my comparison): the 2019 f is cos(theta)(sin(theta) - sqrt(alpha) cos(theta)), which equals the 2021
f(theta) = sin(theta - theta_a) sin(theta_b - theta) multiplied by the constant 1/cos(theta_alpha) = sqrt(1 + alpha) (for
theta_a = theta_alpha, theta_b = pi/2). A constant rescaling of f rescales the regularised time t by a constant, so the two sets of equations
agree up to that. The 2019 paper's v and dv/ds line "sqrt(W) + (2rh - v^2/2) Delta" matches the 2021 dv/dt line
"F (2hr - v^2/2) + sqrt(W)" term for term, which independently supports my reading of the 2019 eq. 6 (page-image reading).

Internal consistency checks run on the 2021 text (independent algebra, not a numerical test):
- Using eq. 3 and eq. 4, dv/dtau = v^2/2 + u^2 - V = (h r + u^2/2) as the paper says: (v^2+u^2)/2 - V + u^2/2 = v^2/2 + u^2 - V. Holds.
- Eq. 8 follows from eq. 4 with u = w sqrt(W)/f: multiply h r = (v^2 + u^2)/2 - V by 2 f^2 and use F^2 = f^2/W. Holds.
- Proposition 4 proof: at theta = theta_b, w = 0, f(theta_b) = 0, the last bracket of dw/dt reduces to f'(theta_b) = -sin(theta_b - theta_a),
  as printed. Holds for theta_b - theta_a != pi.

## 5. The dynamics on the collision manifold

- Lemma 1: on C every maximum or minimum of theta has w = 0 and either theta = theta_a,b or v^2 = 2 V(theta).
- Proposition 4: if at t0 either theta > theta_c with w > 0, or theta < theta_c with w < 0, the orbit reaches the section Sigma
  (below) at least once. Proof: theta cannot reach a maximum at theta < theta_b (there dw/dt > 0 would make it a minimum), and
  it cannot approach theta_b asymptotically (dw/dt = -sin(theta_b - theta_a) != 0 there).
- Proposition 5: on C every orbit oscillates between maxima and minima on theta = theta_a,b and/or the curve v^2 = 2 V(theta), w = 0.
  Consequence (Remark): forwards in time theta oscillates infinitely while v -> infinity, or oscillates finitely and tends to E+.
  Any heteroclinic orbit E- to E+ on C makes a finite number of oscillations between theta_a and theta_b. v -> +-infinity along the
  right arm (theta > theta_c) or left arm (theta < theta_c) of C.
- Proposition 3: for every h < 0 there is a homothetic orbit with theta = theta_c, w = 0, r(t) -> 0 as t -> +-infinity. It is an
  ECO with no partial collision.

Poincare section (Definition 2): Sigma = Sigma_a union Sigma_b = {r >= 0, w = 0, theta = theta_a} union {r >= 0, w = 0, theta = theta_b};
w = 0 is forced on it by eq. 8. Map P : Sigma -> Sigma is the first return forward in time (eq. 11), with P^(-1) backwards. The
section Sigma is preferred to {theta = theta_c} for ECOs, as in the 2019 paper.

### 5.1 One-dimensional manifolds, four non-degenerate types plus degenerate ones

W^u(E-) lies in C and has two branches W^u_+ (w > 0, goes toward theta_b) and W^u_- (w < 0, toward theta_a). Their successive
intersections with Sigma are u^+-_j (increasing sequences), and those of W^s(E+) are s^+-_j with s^-_j = -u^+_j and s^+_j = -u^-_j (symmetry).
Code I+ : W^u(E-) -> S records the type, a or b, of each intersection. A branch makes a full turn when theta goes from
theta_c out through theta_a and theta_b once and back to theta_c. Lacomba (1983) names three cases: non-degenerate (no heteroclinic
connection on C), symmetric degenerate (both branches of W^u(E-) coincide with both of W^s(E+)), non-symmetric degenerate (one does).

Table 1 (non-degenerate; n >= 1 turns):

| Type | I+(W^u_+(E-)) | I+(W^u_-(E-)) | arms |
| --- | --- | --- | --- |
| I | (b,a)^n, b, a, b, b, b, ... | (a,b)^n, a, b, a, a, a, ... | both branches n full turns, escape by different arms |
| II | (b,a)^(n+1), a, a, a, ... | (a,b)^(n+1), b, b, b, ... | n turns and a half, different arms |
| III | (b,a)^n, b, a, b, b, b, ... | (a,b)^(n+1), b, b, b, ... | n and n+1/2 turns, both escape by the right arm |
| IV | (b,a)^(n+1), a, a, a, ... | (a,b)^n, a, b, a, a, ... | n+1/2 and n turns, both escape by the left arm |

(The printed table typesetting is garbled in the text layer; the rows above are rebuilt from the Section 3.1 text, which gives each
sequence in full. Check Table 1 on a page image before relying on the exact repeat counts.) The ordering of the u and s points
along C intersect Sigma_b and Sigma_a is printed for each type (I: s^-_(2n+3) < s^-_(2n+2) < s^-_(2n+1) < u^+_1 < s^-_(2n) < u^+_2 < ...). It is the
tool the proofs use to show that backward images of arcs hit given arcs; I do not reproduce the full chains here.

### 5.2 Two-dimensional manifolds

W^u(E+) and W^s(E-) meet C along the arcs Gamma^+-_u and Gamma^+-_s (eq. 12): each is an orbit escaping forwards (resp. backwards)
to v -> +-infinity through one of the legs of C. Their partial-collision codes are I+(Gamma^+_u) = (b,b,b,...), I+(Gamma^-_u) =
(a,a,a,...), I-(Gamma^+_s) = (a,a,a,...), I-(Gamma^-_s) = (b,b,b,...); intersections P^+_j (increasing) on Sigma_b and P^-_j on
Sigma_a for W^u(E+); Q^+-_j (decreasing, q^+-_j = -p^-+_j) for W^s(E-). Always (eq. 13): q^-_1 < u^+_1, s^-_1 < p^+_1,
q^+_1 < u^-_1, s^+_1 < p^-_1. The homothetic ECO is in W^u(E+) intersect W^s(E-) and does not meet Sigma.
For fixed theta the motion is on an ellipse in the (v, w) plane whose semiaxes are largest at r = 0 (W w^2 + f^2 v^2 <= 2 W f).

## 6. The analytic results: existence and counting of ejection-collision orbits

An ECO of type sigma = (sigma_1, ..., sigma_m) has these partial collisions forwards in time. Definition 1: an ejection orbit
lies in W^u(E+), a collision orbit in W^s(E-), an ECO in the intersection (r -> 0 at both ends). Proposition 6: the time-reversed
orbit of an ECO of type sigma is an ECO of type the reversed sequence.

- Lemma 2: the first intersection of W^u(E+) with Sigma_b is an arc J^b with endpoints u^+_1 and p^+_1; with Sigma_a an arc J^a with
  endpoints u^-_1 and p^-_1; of W^s(E-) the arcs K^b (q^-_1, s^-_1) and K^a (q^+_1, s^+_1). Proof idea: a semicircle of initial conditions
  phi in [0, pi] on W^u(E+) near E+; the ends phi = 0, pi lie on Gamma^+_u, Gamma^-_u, the midpoint phi = pi/2 on the homothetic orbit; orbits
  near the homothetic one pass close to E- and then follow its unstable branch.
- Theorem 1 (any h < 0, any non-degenerate or degenerate case): for every m >= 1 there are ECOs of type (a)^m and of type (b)^m.
  Proof: J^b meets K^b (using eq. 13); then iterate P^(-1) on arcs of K^b, each time obtaining the arc (q^-_(m+1), s^-_1) which meets J^b.
- Theorem 2 (type I, n >= 1 full turns): (a) there are ECOs with 2n+1 collisions of type (b,a)^n, b and (a,b)^n, a (that is, (b,a,...,b,a,b)
  and (a,b,...,a,b,a)); (b) every sequence reachable on the printed directed graph exists, whose vertices are (b,a)^n b, (b), (a), (a,b)^n a and
  whose arrows (all proved, reversals by Proposition 6) give: (b,a)^n b followed by (b)^m, any m; (a)^m followed by (b,a)^n b;
  and (b,a)^n b linked to (a,b)^n a (the proof builds (a,b)^n a (b,a)^n b ... as written in the proof as (a,b,...,a,b,a,b,a,...,a,b);
  the text layer loses which concatenations are printed, so see the page).
- Theorem 3 (type II, n and a half turns): the same with 2(n+1) collisions, (b,a)^(n+1) and (a,b)^(n+1), and the analogous graph.
- Theorem 4 (types III and IV): for every integer k >= 1 there are ECOs of types (b,a)^(k(n+1)), b linked to (b) and to (b,a)^(k(n+1)-1) b
  (type III), and the same with a and b interchanged (type IV). Proof by induction on k, applying the Theorem 2 arc argument to successive arcs.
- Theorem 5 (non-symmetric degenerate, types D1 and D2, coincident branches of n and a half turns): for every k >= 1, ECOs of types
  (a), (a,b)^(k(n+1)) and (b) for D1; (a), (b,a)^(k(n+1)) and (b) for D2. In the symmetric degenerate case only Theorem 1 ECOs are guaranteed.

Notes on counting: the theorems give existence, not a count; the paper states that in the figures arcs are drawn as if they crossed
only once, and that in general they need not ("in the proofs the reader will find points like E and E-tilde that in the figures seem to be the same
one, but they are not in general"). It also notes that other ECOs with fewer partial collisions may exist but cannot be ensured.

### 6.1 Cross-check against the 2019 numerical counts (my comparison)

At alpha = 1 the 2019 paper reports, by order n (number of binary collisions), 2, 2, 2, 2, 4, 8, 12 ECOs (mirror partners
excluded; with 1 = single collision = b = theta = pi/2 and 2 = double collision = a). Theorem 1 gives the all-1 and all-2 sequences at every
order: that is the 2 per order for n <= 4. At order 5 the 2019 paper finds four, all symmetric: (1,1,1,1,1), (2,2,2,2,2), (1,2,1,2,1),
(2,1,2,1,2). Theorem 2(a) with n = 2 (type I, 2 full turns) guarantees exactly the last two, (b,a,b,a,b) and (a,b,a,b,a), at order 2n + 1 = 5.
So the numerical counts up to order 5 coincide with the guaranteed ones, and the first non-guaranteed ones appear at order 6. This is
consistent with alpha = 1 being a type-I case with n = 2: the 2019 paper places alpha = 1 in (alpha_3, alpha_4) and quotes
Lacomba and Medina for the order-5 ECO (1,2,1,2,1) there; the 2021 paper's Fig. 3 top left is described as "two full turns" and the figure
set is stated to use SC4BP samples, but the text does not say which alpha belongs to which panel, so the identification of alpha = 1 with
type I, n = 2 is my inference. Also the 2019 observation that the orbit set is closed under the swap 1 <-> 2 at alpha = 1 is not predicted by
the theorems here.

## 7. Printed numbers usable as sourced tests

The paper prints no numerical values. What is testable is structural, and each item below is stated by the paper:

1. dv/dtau = h r + u^2/2 >= 0 on r = 0 (so v is non-decreasing on C): a monotonicity test for any implementation of eq. 3 or eq. 7.
2. Eq. 4 and eq. 8 (energy relations) must be conserved along any integration of eqs. 3 and 7 respectively.
3. The scaling symmetry: a solution with energy h gives a solution (lambda r, v, theta, w) with energy lambda h.
4. The time-reversal symmetry (r, v, theta, w, t) -> (r, -v, theta, -w, -t) (eq. 9, 10).
5. Equilibria E+- = (0, +-v_c, theta_c, 0) with v_c^2 = 2 V(theta_c) for eq. 3 in its own normalisation; hyperbolic saddles with the
   stated stable and unstable dimensions (1 and 2 on a fixed energy level).
6. At theta = theta_b, w = 0: dw/dt = -sin(theta_b - theta_a) (when theta_b - theta_a != pi).
7. Counts of guaranteed ECOs (Theorem 1: all (a)^m, (b)^m, m >= 1; Theorem 2: 2n+1 collisions; Theorem 3: 2(n+1) collisions), to be compared with
   a numerical search of the 2019 kind. The first cases line up with the 2019 numbers (Section 6.1); the identification of the
   mass parameter with the type and n needs a printed alpha (not in the text) or a computation of the one-dimensional manifold.

No number can be taken from this paper as an expected value for an orbit, an eigenvalue or an alpha_k.

## 8. Techniques applicable to the project's problems

### 8.1 `#928` (regularised propagator for close passes)

- The paper does not help for a planar pass: the restricted problem is outside its class (Section 3), and its partial-collision
  regularisation is a one-dimensional Sundman time change (the f(theta) factor), valid because the collision shape-manifold is an interval.
  It confirms what the 2019 digest concluded: McGehee blow-up is a complement for the asymptotics of ejection and collision, not a
  replacement for Levi-Civita or KS.
- The part with a general lesson: when the shape space is an interval, one multiplicative time factor f = sin(theta - theta_a) sin(theta_b - theta)
  removes both end singularities at once and gives a smooth vector field on a compact interval; the same product form is how a
  two-centre collision shape manifold would be regularised in the three-body setting.
- No test case for the `#928` controls (radial fall, Kepler ellipse, symplectic check, finite-difference Jacobian). The structural
  checks of Section 7 items 1 to 4 apply only if a McGehee-class propagator is ever built.

### 8.2 `#899` (ejection-collision orbits of the restricted problem as seeds)

- What carries over is the proof architecture, not the equations: take the one-parameter family of ejection orbits (an arc of
  the 2-dimensional unstable manifold, parametrised by an angle), follow it to a section built from the partial-collision set, and use
  backward iterates of the stable arc to find intersections one return at a time. Numerically this is the 2019 bracketing scheme (sign
  change of a continuous function of the angle at the n-th section crossing). The 2021 paper supplies the reason the bracket must exist:
  the end arcs J^a, J^b and K^a, K^b have ordered endpoints (eq. 13) and are mapped onto ordered arcs by P^(-1).
- A concrete guide for the restricted problem: define the section as the set of collisions with either primary (analogue of Sigma_a and Sigma_b
  in the problem with two partial-collision types), code each ejection orbit by the sequence of such collisions, and treat the
  existence of ECOs of type (a)^m and (b)^m for every m (Theorem 1 analogue) as the baseline family, with mixed sequences needing the one-dimensional
  manifold structure. The restricted problem has two collision types by construction (the two primaries), which is the structural
  parallel; the break is that there is no total collision to blow up.
- Do not use this paper for the choice of regularisation or for the mu -> 0 continuation: both are in the Olle, Rodriguez and Soler papers,
  which are not held. See the follow-ups.

## 9. Cited references (with printed DOIs) and follow-ups

Follow-ups (no task numbers registered):
- Acquire and digest the restricted-problem papers named in the `#899` entry of `data/OUTSTANDING.md`: Olle, Rodriguez and Soler, CNSNS 55:298
  (2018), DOI 10.1016/j.cnsns.2017.07.013; CNSNS 90:105294 (2020), DOI 10.1016/j.cnsns.2020.105294; CNSNS 94:105550 (2021),
  DOI 10.1016/j.cnsns.2020.105550; CNSNS 111:106410 (2022), DOI 10.1016/j.cnsns.2022.106410. This paper cites none of them.
- Re-read Table 1 and the Theorem 2 to 5 diagrams on a page image before writing a test on them.
- Test the claimed mapping alpha = 1 to type I, n = 2 by computing the one-dimensional manifold on C for the symmetric collinear four-body problem
  and counting turns; only worth doing if the sequence-counting approach is pursued for `#899`.
- Kaplan 1999 (symbolic dynamics of the collinear three-body problem): the closest source for the Poincare-map coding with the same section;
  Devaney 1980 and Moeckel 1987 for spiralling manifolds (complex eigenvalues), which is the case this paper excludes by assuming real saddles.

References printed with DOIs (none held in the corpus unless noted in CORPUS_INDEX; I did not re-grep each):
McGehee 1974, Invent. Math. 27:191-227, 10.1007/BF01390175; Devaney 1980, Invent. Math. 60:249-267, 10.1007/BF01390017;
Kaplan 1999, Contemp. Math. 246:143-162, 10.1090/conm/246/03781; Lacomba and Simo 1982, Celest. Mech. 28:37-48 (boundary manifolds for energy
surfaces), 10.1007/BF01230658; Simo and Lacomba 1982, Celest. Mech. 28:49-62, 10.1007/BF01230659; Lacomba 1983, Celest. Mech. 31:23-41,
10.1007/BF01272558; Lacomba and Medina 2004, Qual. Theory Dyn. Syst. 5:75-100, 10.1007/BF02968131; Lacomba and Medina 2008, DCDS-S 1:557-587,
10.3934/dcdss.2008.1.557; Lacomba and Perez-Chavela 1993, CMDA 57:411-437, 10.1007/BF00695713; Delgado Fernandez and Perez-Chavela 1991,
10.1007/978-1-4613-9725-0_8; Martinez 2012, DCDS-B 17:943-975, 10.3934/dcdsb.2012.17.943; Martinez 2013, CMDA 117:217-243,
10.1007/s10569-013-9509-4; Simo and Llibre 1981, Arch. Rational Mech. Anal. 77:189-198, 10.1007/BF00250623; Moeckel 1989, Arch. Rational
Mech. Anal. 107:37-69, 10.1007/BF00251426; Tanikawa and Mikkola 2000a, CMDA 76:23-34, 10.1023/A:1008397420958; 2000b, Chaos 10:649-657,
10.1063/1.1287064; Alvarez-Ramirez and Medina 2020, Qual. Theory Dyn. Syst. 19:10-15, 10.1007/s12346-020-00342-z; Alvarez-Ramirez, Medina and Vidal 2015,
Astrophys. Space Sci. 358:1-17, 10.1007/s10509-015-2416-2; Saari 2005, CBMS 104, 10.1090/cbms/104. Without a printed DOI: Broucke 1979 (Astron.
Astrophys. 73:303-313), Simo 1981 (Dekker), Moeckel 1981, 1983, 1984, 1987, 2007, Sekiguchi and Tanikawa 2004.

## 10. Cautions

- The 2021 normalisation v_c^2 = 2 V(theta_c) differs from the 2019 v_c^2 = V(theta_c); I did not reconcile the two scalings of the
  momentum and the energy, so numbers must not be moved between the two papers.
- Statements in Section 6 about which concatenations appear on the Theorem 2 and 3 graphs are partly lost in the text layer (arrows).
- The identification of the 2019 alpha = 1 case with type I, n = 2 is my inference, not a printed fact.
- The assertion that the restricted and Hill problems are outside the hypotheses is my analysis of the stated conditions, not a statement of the paper.
