# Digest: Alvarez-Ramirez, Barrabes, Medina and Olle 2019, ejection-collision orbits in the symmetric collinear four-body problem

Date: 2026-10-05 (Sydney). Reading, reasoning and two small independent computations; no project code was changed.

Source: M. Alvarez-Ramirez, E. Barrabes, M. Medina and M. Olle, "Ejection-Collision orbits in the symmetric collinear
four-body problem", Commun. Nonlinear Sci. Numer. Simulat. 71:82-100 (2019), DOI 10.1016/j.cnsns.2018.10.026, received 3 May
2018, accepted 31 October 2018. Filed in the private paper corpus as
`alvarez-ramirez-barrabes-medina-olle-2019-ejection-collision-orbits-symmetric-collinear-four-body-cnsns-71-82-doi-10.1016-j.cnsns.2018.10.026.pdf`.
A 19-page PDF with a good text layer except equations (eqs. 6, 18, 24 and the Proposition 1 inequalities), which I checked on
page images (journal pp. 85, 94-97); the four "Tables" are images of orbit plots labelled by their collision sequence.

## 1. Relevance verdict

Direct physical relevance to the project: none. The model is a one-dimensional (collinear), symmetric four-body problem with
masses (alpha, 1, 1, alpha); there is no rotating frame, no restricted problem and no spacecraft. Nothing in it is a catalogue
candidate or a tolerance source.

Methodological relevance: moderate, and mostly conceptual. The paper is a clean worked example of three things the project
needs for `#899`, `#928` and the Broucke-1969 collision family (`#933`):

- a regularisation by McGehee blow-up (the collision is replaced by an invariant manifold with a flow on it), which is a
  different object from Levi-Civita or KS (Section 4);
- ejection-collision orbits defined as heteroclinic connections between equilibria on the collision manifold, and counted by
  a symbol code (Section 3);
- a root-finding scheme that locates every such connection by tracking a continuous sign function, not a discontinuous symbol
  (Section 3.3).

It is NOT a replacement for KS or Levi-Civita for a planar close pass (Section 4.1). It provides no test case that our code
can reproduce without first implementing the symmetric collinear four-body problem, which has no project use.

## 2. The model and its regularisation

### 2.1 Equations (page 83-85)

Masses m1 = m2 = 1 at positions +-x; m3 = m4 = alpha at positions +-y/sqrt(alpha); order m3, m1, m2, m4 along the line;
motion symmetric about the centre of mass at the origin. Hamiltonian (eq. 1), with p_x = 2 dx/dt, p_y = 2 dy/dt:

    H = p_x^2/4 + p_y^2/4 - U(x,y)
    U = 1/(2x) + alpha^(5/2)/(2y) + 2 alpha^(3/2)/(y - sqrt(alpha) x) + 2 alpha^(3/2)/(y + sqrt(alpha) x)        (2)

Domain 0 < sqrt(alpha) x < y. Three singularities: single binary collision SBC (x = 0, m1 with m2), double binary collision
DBC (y = sqrt(alpha) x, m1 with m3 and simultaneously m2 with m4), quadruple collision (x = y = 0). Hill region
R_h = {U >= -h}; bounded motion needs h < 0. The paper works at alpha = 1, h = -1.

### 2.2 Blow-up and binary regularisation (page 85, eqs. 4-7)

1. Polar variables: x = r cos(theta)/sqrt(2), y = r sin(theta)/sqrt(2), with
   p_x = sqrt(2) p_r cos(theta) - sqrt(2) p_theta sin(theta)/r (the printed layout puts the 1/r on the p_theta term), and
   theta in (theta_alpha, pi/2), theta_alpha = arctan(sqrt(alpha)) the DBC ray.
2. Potential on the unit shape sphere (eq. 4):
   V(theta) = r U = 1/(sqrt(2) cos(theta)) + alpha^(5/2)/(sqrt(2) sin(theta))
   + 2 sqrt(2) alpha^(3/2)/(sin(theta) - sqrt(alpha) cos(theta)) + 2 sqrt(2) alpha^(3/2)/(sin(theta) + sqrt(alpha) cos(theta)).
3. McGehee coordinates: p_r = r^(-1/2) v, p_theta = r^(1/2) u, time dt = r^(3/2) d tau.
4. Remove both binary collisions at once with the regularised potential
   W(theta) = V(theta) cos(theta) (sin(theta) - sqrt(alpha) cos(theta)),
   positive and real analytic on [theta_alpha, pi/2], and a second time change with
   Delta(theta) = d tau/ds = cos(theta) (sin(theta) - sqrt(alpha) cos(theta))/sqrt(W(theta)), w = Delta(theta) u (eq. 5).

Regularised system, time s (eq. 6, read from the page image):

    dr/ds = r v Delta(theta)
    dv/ds = sqrt(W(theta)) + (2 r h - v^2/2) Delta(theta)
    dtheta/ds = w
    dw/ds = -(v w/2) Delta(theta) + (cos 2 theta + sqrt(alpha) sin 2 theta) ((2 r h - v^2)/sqrt(W) Delta(theta) + 1)
            + (W'(theta)/W(theta)) (cos(theta)(sin(theta) - sqrt(alpha) cos(theta)) - w^2/2)

Energy relation (eq. 7): w^2 = (2 r h - v^2) Delta(theta)^2 + 2 cos(theta)(sin(theta) - sqrt(alpha) cos(theta)).
(The text layer of the printed eq. 6 is partly garbled; the layout above is from the page image. The factor in the dv/ds
line is printed as sqrt(W) + (2rh - v^2/2) Delta and I have not re-derived it; treat the coefficient placement as read, not
verified.)

The field is analytic on F = [0, infinity) x R x [theta_alpha, pi/2] x R, including r = 0. Symmetry (eq. 8):
L1 : (r, v, theta, w, s) -> (r, -v, theta, -w, -s), so a solution reversed in s traces the same path in (r, theta). A solution
is symmetric (Definition 1) if q(s0 + s) = q(s0 - s) for some s0; v = w = 0 at a point means a binary collision or a point on
the zero-velocity curve.

### 2.3 Total collision manifold, equilibria, invariant manifolds (pages 86-88)

- C = {r = 0, w^2 = -v^2 Delta^2 + 2 cos(theta)(sin(theta) - sqrt(alpha) cos(theta))}: a 2-D invariant manifold, a sphere
  minus four points, independent of h; the flow on C is gradient-like in v (dv/ds >= 0).
- Two equilibria E+- = (0, +-v_c, theta_c, 0) with V'(theta_c) = 0 and v_c^2 = V(theta_c) (printed in the Appendix). The
  homothetic solution theta = theta_c divides the plane into the DBC region (theta < theta_c) and the SBC region (theta > theta_c).
- Hyperbolic for every alpha > 0: four distinct real eigenvalues, two positive and two negative (Lemma 1, proved with
  polynomials in w = tan(theta_c)/sqrt(alpha); see Section 5 below for the numerical check I did).
- On the energy surface: dim W^u(E+) = dim W^s(E-) = 2; dim W^u(E-) = dim W^s(E+) = 1, the latter two lying inside C.
- Eigenvectors (Appendix): sigma_1 = (-v_c, -+h, 0, 0) (slow direction, equals the homothetic solution), sigma_3 = (0, 0, 1, lambda_3),
  sigma_4 = (0, 0, 1, lambda_4). Eigenvalues (eq. 18) lambda_1(E+-) = +-lambda, lambda_2 = -lambda_1, lambda_3,4 =
  (lambda/4)(-+1 +- sqrt(1 + 8 V''(theta_c)/V(theta_c))) with lambda = 2 cos(theta_c)(sin(theta_c) - sqrt(alpha) cos(theta_c)) > 0.
  (The sign pattern of the -+ and +- in the printed eq. 18 is hard to read in the text layer; see the caveat in Section 5.)
- Order-one parametrisation of W^u(E+) (eqs. 11, 24):
  Gamma_1^+(xi, phi) = E+ + xi (cos(2 pi phi) sigma_1 + sin(2 pi phi) sigma_3), with unit-normalised sigma_i. The r component of the correction is
  proportional to -cos(2 pi phi) (the first entry of eq. 24), so r > 0 needs phi in (1/4, 3/4); phi = 1/4, 3/4 are
  orbits inside C (fast direction), phi = 1/2 is the homothetic orbit (slow direction). phi in (1/4, 1/2) starts in the
  SBC region, phi in (1/2, 3/4) in the DBC region.
- Higher-order parametrisations Gamma_m^+ (m <= 8) built with the parameterisation method of Haro, Canadell, Figueras, Luque and
  Mondelo (2016) (filed in the corpus; digested in `2026-07-27-733-olikara-thesis-haro-parameterization-book-digest.md`).

## 3. Ejection-collision orbits: definition, computation, counting

### 3.1 Definitions

- ECO (Definition 3): a solution with r(s) -> 0 as s -> +-infinity, that is a heteroclinic connection E+ -> E- in
  W^u(E+) intersect W^s(E-). It starts and ends at quadruple collision.
- Proposition 2: a symmetric solution in W^u(E+) is an ECO; the reversed solution of an ECO is an ECO.
- Section 2.3 of the paper: Sekiguchi and Tanikawa (2004) use the section {theta = theta_c}, on which different ECOs share boundary
  points; the paper also states that their Theorem 1 (every orbit crosses {theta = theta_c} at least once) is false, because their
  proof omits W^u(E+) and W^s(E-). The paper's own Theorem 1 supplies the counter-example.
- Section Sigma_c (eq. 12): Sigma_c = {w = 0, theta = theta_alpha or theta = pi/2}, the set of binary collisions of either type.
- Symbol map P (eq. 16): P(Gamma) = (p_1, p_2, ...), p_j = 1 if the j-th crossing of Sigma_c is an SBC, 2 if a DBC. B is the set
  of sequences over {1, 2}.
- Proposition 3: P(Gamma) is finite iff Gamma is an ECO (proof: infinitely many binary collisions would need an unbounded
  analytic theta(s)); P of the reversed ECO is the reversed sequence; symmetric ECO iff symmetric (palindromic) sequence.
- ECO of order n: P(Gamma) = (p_1, ..., p_n). Proposition 4: a symmetric ECO of order n touches the zero-velocity curve iff n is
  even, between the (n/2)-th and (n/2 + 1)-th binary collision.
- Lacomba and Medina (2004) [their ref. 20], quoted: for any alpha and any n there is an ECO with exactly n SBC (or n DBC); for
  alpha in (alpha_3, alpha_4) an ECO with P = (1,2,1,2,1) exists, in particular at alpha = 1. The sequence {alpha_k} is from
  Simo and Lacomba (1982); alpha = 1 lies in (alpha_3, alpha_4) and alpha = 2 in (alpha_4, alpha_5) (Fig. 3). The values alpha_k
  are not printed in this paper.

### 3.2 Direct escape orbits (Section 3)

Escape type 1: y -> infinity with x bounded (only SBC); type 2: x, y -> infinity (only DBC). Proposition 1 gives sufficient
escape inequalities in (v, theta, w) (a positive two-body energy bound on the outer body, or on a Jacobi coordinate). Theorem 1:
orbits ejecting from quadruple collision and escaping with binary collisions of only one type exist, proved by linearising
about E+ and showing the sign of DF(p0).p1 at phi = 1/4 + epsilon (type 1) or 3/4 - epsilon (type 2). Numerically (alpha = 1,
h = -1, xi = 0.01, order m = 8): all phi in (1/4, 0.49998) escape directly with SBC only, and all phi in (0.500006, 3/4) escape
directly with DBC only; the paper prints "alpha = -1" for this run, which is a slip for alpha = 1 (its own captions say alpha = 1).

### 3.3 Computing ECOs (Section 4.1)

1. Start on W^u(E+) with the order-m parametrisation Gamma_m^+(xi, phi) at fixed xi; each orbit is labelled by phi in (1/4, 3/4).
2. Integrate system (6) forward to the first n+1 crossings of Sigma_c; define P_n(phi) = (p_1, ..., p_n).
3. Proposition 5 (by continuity): if P_{n+1}(phi_1) and P_{n+1}(phi_2) agree in the first n symbols and differ in the (n+1)-th,
   there is an orbit between them with P = (p_1, ..., p_n), that is an ECO of order n.
4. Instead of the discontinuous theta_{n+1}(phi) (jumps between pi/2 and theta_alpha), use the continuous function
   F_{n+1}(phi) = r (theta - theta_c) evaluated at the (n+1)-th crossing: r > 0, so its sign is the sign of theta - theta_c.
   Every zero of F_{n+1} is an ECO of order j <= n. Track the sign of F_{n+1} over a phi grid and refine each bracket by an
   iterative root finder.
5. Observed (not proved): every zero found is transversal, so Proposition 5 would be an "if and only if".

Numerical controls reported:

- ECOs lie in phi in (1/2 - epsilon, 1/2 + epsilon): epsilon about 1e-4 for xi about 1e-2 at alpha = 1. The higher the order n,
  the narrower the interval containing zeros of F_{n+1}.
- Order-one parametrisation needs xi < 1e-6 for an orbit accuracy about 1e-12, and then epsilon about 1e-9 is needed to see the
  structure; hence the high-order parametrisation, which allows a large xi.
- Orbit-error accuracy (Fig. 4, alpha = 1, h = -1): below 1e-8 with m = 5 needs xi up to about 0.015, with m = 8 up to about
  0.05. The run used m = 8 and xi about 1e-2 (xi = 0.05 for Fig. 7). Results were repeated at xi = 0.001, 0.01, 0.05 with the
  same ECOs.
- Fig. 7 (alpha = 1, h = -1, xi = 0.05, phi < 1/2): for n = 4 the function F_5 has four zeros (four ECOs of order <= 4); for
  n = 5, F_6 has six zeros (the same four plus two new of order 5); for n = 6, F_7 shows more (five are plotted in a detail).

## 4. How this relates to KS and Levi-Civita (`#928`)

### 4.1 Alternative or complement?

A complement, and for the planar close-pass problem of `#928` only a weak one.

- What the paper regularises by time scaling is a one-dimensional binary collision. In the collinear problem the two bodies meet
  and reverse (the orbit returns along the same path); the factor cos(theta)(sin(theta) - sqrt(alpha) cos(theta)) in Delta,
  together with the rescaled w, removes the 1/distance singularity. This is the Sundman-type multiplication, which suffices in
  one dimension. It does not turn a planar near-collision into a smooth rotation: the planar case needs Levi-Civita (the u^2
  map), which this paper never uses. So for passes through or near the small primary in the planar restricted problem, KS or
  Levi-Civita remains the tool; this paper does not replace it.
- What the paper adds is McGehee blow-up of the total (quadruple) collision: r -> 0 is replaced by a boundary manifold C with its
  own flow, with hyperbolic equilibria whose stable and unstable manifolds are the ejection and collision orbits. This is
  complementary: Levi-Civita and KS make a collision orbit a regular passage; McGehee makes the asymptotics of an
  ejection or collision orbit an invariant-manifold question.
- In the restricted problem the analogue of "quadruple collision" is the collision with the small primary (a binary collision, so
  Levi-Civita already regularises it). A McGehee blow-up of that binary collision is possible (the equilibria are then the
  radial-fall states) but gives nothing that Levi-Civita does not, for ballistic passes. The McGehee view would matter for ECOs
  at the mu = 0 skeleton only as a classification device.
- Combination worth noting: the paper uses a regularised time in which every binary collision is a regular point at which the
  orbit reflects, then reads the collision sequence off a section through those points. A two-stage structure (Levi-Civita for
  the pass, section at the collision for the symbol) is the natural planar analogue.

### 4.2 Test case for `#928`

None directly. The model is not the CR3BP, so no printed number can check a Levi-Civita CR3BP propagator. The only
transferable regression is a method check, not a number: on the closed-form radial fall, a regularised integrator should keep
the collision as a regular point, and the ejection angle (the analogue of phi) bracket structure is a topological check (see 5.2).

## 5. Printed numbers usable as sourced tests

Honest summary: the paper prints very few numerical values. There is no orbit initial condition, no period, no phi value at an
ECO, no eigenvalue, no theta_c and no v_c. What is printed is a set of counts and sequences, one polynomial identity chain and
a few tolerances.

### 5.1 Counts and symbol sequences (alpha = 1, h = -1; tables are plot images read from the page)

| order n | ECOs found | symmetric | sequences shown in the table (mirror partners omitted) |
| --- | --- | --- | --- |
| 1 | 2 | (all) | (1), (2) |
| 2 | 2 | (all) | (1,1), (2,2) |
| 3 | 2 | (all) | (1,1,1), (2,2,2) |
| 4 | 2 | (all) | (1,1,1,1), (2,2,2,2) |
| 5 | 4 | all 4 | (1,1,1,1,1), (2,2,2,2,2), (1,2,1,2,1), (2,1,2,1,2) |
| 6 | 8 | 2 | symmetric (1^6), (2^6); non-symmetric (1,2,1,2,1,2), (1,1,2,1,2,1), (2,2,1,2,1,2); their mirrors (2,1,2,1,2,1), (1,2,1,2,1,1), (2,1,2,1,2,2) have the same (x, y) projection |
| 7 | 12 | 4 | symmetric (1^7), (2^7), (1,1,2,1,2,1,1), (2,2,1,2,1,2,2); non-symmetric (1,1,1,2,1,2,1), (1,1,2,1,2,1,2), (1,2,1,2,1,2,2), (2,2,2,1,2,1,2); their mirrors (1,2,1,2,1,1,1), (2,1,2,1,2,1,1), (2,2,1,2,1,2,1), (2,1,2,1,2,2,2) |

Further printed statements: for n <= 4 only the all-1 and all-2 sequences exist; the first non-symmetric ECOs appear at n = 6;
if a sequence (p_1, ..., p_n) is an ECO then so is (3 - p_1, ..., 3 - p_n) at alpha = 1 (observed, no proof; the
alpha = 2 manifolds of Fig. 3 differ, so this is not general). The sequences of Tables 3 and 4 that I read at 90 dpi should be
rechecked on a full-resolution page image before they are used as a test, particularly the two non-symmetric order-7 labels
(1,1,1,2,1,2,1) and (2,2,2,1,2,1,2).

Fig. 8: a directed graph of allowed consecutive collision-type transitions for alpha = 1, from Lacomba and Medina. The paper's
only stated consequence: no ECO of the form (2,...,2,1,1) or (1,...,1,2,2) can exist, because the partial sequences (2,1,1)
and (1,2,2) cannot occur. I did not reconstruct the graph's adjacency from the image, and I did not check the table sequences
against it.

### 5.2 Tolerances and sizes (alpha = 1, h = -1)

| quantity | printed value |
| --- | --- |
| orbit accuracy below 1e-8, order m = 5 | xi up to about 0.015 |
| orbit accuracy below 1e-8, order m = 8 | xi up to about 0.05 |
| order-one parametrisation accuracy about 1e-12 | xi below 1e-6, with epsilon about 1e-9 |
| working choice | m = 8, xi about 1e-2 (0.05 in Fig. 7), epsilon about 1e-4 |
| direct escape, SBC only | phi in (1/4, 0.49998), xi = 0.01, m = 8 |
| direct escape, DBC only | phi in (0.500006, 3/4), xi = 0.01, m = 8 |
| F_5 zeros for phi < 1/2 | 4 |
| F_6 zeros for phi < 1/2 | 6 |

### 5.3 Analytic identity chain (Appendix, Lemma 1) and my numerical check

Printed: with w = tan(theta_c)/sqrt(alpha) the condition V'(theta_c) = 0 becomes
w^7 - 2 w^5 - (8 + 17 alpha) w^4 + w^3 + (2 alpha - 8) w^2 - alpha = 0, linear in alpha, so
alpha = w^2 (w^5 - 2 w^3 - 8 w^2 + w - 8) / (17 w^4 - 2 w^2 + 1) (eq. 23); alpha > 0 requires w above the unique positive
root w_bar of w^5 - 2 w^3 - 8 w^2 + w - 8 = 0, with w_bar in [11/5, 12/5].

My computation (independent; not printed in the paper). I coded V(theta) of eq. 4 at alpha = 1, minimised it, and tested eq. 23
and the printed quintic:

| quantity | value |
| --- | --- |
| theta_c | 1.264505 rad |
| tan(theta_c) = w (alpha = 1) | 3.16212 |
| V(theta_c) = v_c^2 | 9.67900 |
| V''(theta_c) (finite difference) | about 90.1 |
| eq. 23 evaluated at that w | 1.0000000001 (alpha = 1 recovered) |
| printed septic residual at that w | 1.8e-7 (finite-difference noise) |
| positive real root of the printed quintic | 2.39681 (inside the printed [11/5, 12/5]) |

This confirms my transcription of eq. 4 and the Appendix polynomials against each other (the formula chain is internally
consistent with V'(theta_c) = 0 at alpha = 1). It does not test eqs. 6 or 18.
Caveat on the eigenvalues: reading eq. 18 as lambda_3 = (lambda/4)(-1 + sqrt(1 + 8V''/V)) for E+ gives, from my numbers,
lambda = 0.393 and lambda_3, lambda_4 = 0.756, -0.952; these are not printed and I cannot confirm the equation sign pattern
from the text layer, so they are not a sourced test. The paper's printed theorem (Lemma 1) is that lambda_1 < lambda_3, which
my reading satisfies.

## 6. Techniques applicable to the project's problems

### 6.1 For `#928` (regularised propagator for close passes)

- Verdict: not a replacement for KS or Levi-Civita. See Section 4.1 above. The paper's regularisation is a McGehee blow-up plus a
  one-dimensional Sundman scaling; none of it applies to planar or spatial near-collisions beyond the collinear fall.
- The useful carry-over for the collinear limit: a radial fall in the planar CR3BP is the collinear collision case. The same
  Delta-scaling (multiply time by a factor that vanishes like the separation) is what `core/cr3bp_regularized.py` already does as
  its Sundman transformation, so the paper confirms that a Sundman-only integrator is sufficient on a radial fall and
  insufficient off it; that matches the existing `#928` plan to test against the radial-fall closed form.
- Sourced test: none that the project can run today. The only test the paper supports is topological (bracket structure on a
  phi-like parameter, Section 6.2), not numerical.

### 6.2 For `#899` (second-species orbits; ejection-collision orbits of the restricted problem as seeds)

- Direct analogue. In the restricted problem a collision orbit with the small primary is determined by a one-parameter
  ejection angle at a given Jacobi constant, in the same way that W^u(E+) is a one-parameter family in phi at fixed energy.
  The paper's recipe carries over: shoot a fine grid of the ejection angle, follow each ejected orbit to the n-th crossing of
  a chosen section, and treat the section itself as the seed of the symbol sequence.
- Use a continuous bracketing function, not the discontinuous symbol. The paper's function F_{n+1} = r (theta - theta_c) is
  continuous, changes sign exactly where the symbol jumps, and is therefore bisectable; a code that bisects on the discrete
  "which crossing came next" label has no gradient. For the restricted problem, an analogue is a signed distance to the
  separatrix between two kinds of next pass (inside or outside a lunar-pass section), multiplied by a positive scale, to make
  every sign change a root.
- Why this finds all of them: the symbol changes only across an ECO, so a tabulated grid with bracketing refines every order-j
  ECO with j below the depth; the paper reports the zeros become transversal in every case it checked, which is the condition
  under which a sign-change finder never misses a root.
- Sensitivity warning for the seed: ECOs concentrate near the slow (homothetic) direction (epsilon about 1e-4 at xi about 1e-2,
  about 1e-9 at the linear level). The analogue is that near-collision seeds cluster close to the collision-ejection direction
  and need a refined grid in the ejection angle there; start at a larger distance from the collision with a higher-order
  local expansion, or the structure is lost in rounding.
- Cross-check against symbol-sequence counting: counts by order n (Section 5.1) are a template for the sort of "number of
  ejection-collision orbits of order n" table `#899` should produce at mu = 0 and at small mu.
- Related corpus items already digested: Gomez-Olle 1991 (`2026-10-04-digest-gomez-olle-1991-second-species-circular-elliptic-I-II.md`),
  Font-Nunes-Simo 2002 and 2009 (`2026-10-04-digest-font-nunes-simo-2002-consecutive-quasi-collisions.md`, `...2009-second-species-numerical-study.md`),
  Perko 1976-1981 second species. This paper's contribution relative to those is the invariant-manifold-based completeness of
  the search, not a new orbit family for the restricted problem.

### 6.3 For Broucke's collision-orbit family (`#933`)

Only as classification vocabulary: ECO of order n and the symmetric-sequence criterion (palindromic P, with the zero-velocity
touching rule of Proposition 4 for even n) give a way to label and count Broucke-type symmetric periodic collision orbits. No
numerical overlap.

## 7. Errata and cautions

- Section 3 numeric run is printed "alpha = -1"; the captions and the rest of the paper use alpha = 1.
- The paper states that Sekiguchi and Tanikawa's Theorem 1 is false because it omits W^u(E+) and W^s(E-). This is the authors'
  claim; I have not checked the 2004 proof.
- Eq. 6 (dv/ds and dw/ds placement) and eq. 18 (sign pattern) are legible on the page image but I did not re-derive them; no
  project code should be written from them without re-deriving from eqs. 1-5.
- Tables 1 to 4 are plot images whose sequence labels I read at 90 dpi.

## 8. Follow-ups (no task numbers registered)

- Full-resolution reread of the Tables 3 and 4 labels and of eq. 6 and eq. 18 before any test is built on them.
- Haro et al. 2016 parameterisation method (corpus, digested under `#733`): the order-m expansion of the local manifold is the
  numerical engine here; reuse the digest if a McGehee-type local expansion is ever wanted at the mu = 0 skeleton of `#899`.
- Lacomba and Medina, "Symbolic dynamics in the symmetric collinear four-body problem", Qual. Theory Dyn. Syst. 5(1):75-100
  (2004), DOI 10.1007/BF02968131 (their ref. 20): the proofs of existence of ECOs and the Fig. 8 transition graph. Not in the
  corpus (grep of CORPUS_INDEX for Lacomba returned no row); acquire if the sequence-counting approach is pursued.
- Simo and Lacomba, "Analysis of some degenerate quadruple collisions", Celest. Mech. 28(1-2):49-62 (1982), DOI 10.1007/BF01230659:
  the sequence {alpha_k} and the collision-manifold flow. Not in the corpus.
- McGehee, "Triple collision in the collinear three-body problem", Invent. Math. 27:191-227 (1974), DOI 10.1007/BF01390175: the
  blow-up used here; the original source of the technique. Not in the corpus index (no row matched McGehee).
- Sekiguchi and Tanikawa, Publ. Astron. Soc. Jpn 56:235-251 (2004): the symbol-sequence classification of the equal-mass case
  (no DOI printed in this paper).
- Restricted three-body ejection-collision papers: this paper does not cite one. The group's restricted-problem work cited by
  nearby digests is Barrabes, Mondelo and Olle (2009, 2009b; in the corpus) and Barrabes and Gomez (2002, 2003; in the corpus), and
  Gomez and Olle (1986, 1991; in the corpus), all about periodic or homoclinic orbits, not ejection-collision orbits. I recall,
  without a source in hand and therefore unverified, that the Barcelona group also has a restricted-problem ejection-collision
  paper; check by a literature search before listing it, and give a DOI only after reading it. Also in the corpus and relevant
  to collision orbits: Bolotin 2006 (shadowing chains of collision orbits, DCDS 14(2):235, DOI 10.3934/dcds.2006.14.235) and
  Font-Nunes-Simo 2002/2009.
- Also cited here and relevant for four-body background only: Alvarez-Ramirez, Medina and Vidal, "The trapezoidal collinear
  four-body problem", Astrophys. Space Sci. 358:1-17 (2015) (no DOI printed); Sweatman, Celest. Mech. Dyn. Astron. 82(2):179-201
  (2002), DOI 10.1023/A:1014599918133; Sweatman 2006, CMDA 94(1):37-65, DOI 10.1007/s10569-005-2289-8; Saari, "Collisions, rings,
  and other Newtonian N-body problems", CBMS 104 (2005), DOI 10.1090/cbms/104.

## 9. Note added 2026-10-05: the 2021 analytic generalisation

The same four authors' 2021 paper (J. Nonlinear Sci. 31:68, DOI 10.1007/s00332-021-09721-5) is the analytic version of this one, digested in
`2026-10-05-digest-alvarez-ramirez-barrabes-medina-olle-2021-ejection-collision-two-dof.md`. It prints no numerics. Its Theorems 1 and 2 reproduce
the counts of Section 5.1 above for orders 1 to 5 (two per order for n <= 4 from the all-1 and all-2 sequences; the four at order 5 are those two plus
(1,2,1,2,1) and (2,1,2,1,2), Theorem 2 with two full turns), and the first ECOs not guaranteed by the theorems are those of order 6.
Its eq. 7 is the 2019 eq. 6 up to a constant time rescaling, with the 2019 dv/ds line matching term for term. Its normalisation has
v_c^2 = 2 V(theta_c), against v_c^2 = V(theta_c) here, so numbers must not be moved between the two. The planar restricted problem is outside the
2021 hypotheses (rotating-frame Coriolis term, no total collision); the restricted-problem papers (Olle, Rodriguez and Soler, listed in the `#899` entry of
`data/OUTSTANDING.md`) are not cited there and are not held.
