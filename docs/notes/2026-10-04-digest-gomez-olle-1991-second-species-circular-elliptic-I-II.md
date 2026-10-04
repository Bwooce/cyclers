# Digest: Gomez and Olle 1991, second-species solutions in the circular and elliptic restricted three-body problem, I and II

Date 2026-10-04. Purpose: the classical analytical (Part I) and numerical (Part II) treatment of periodic second-species
solutions for ONE small secondary, circular AND elliptic problem, read so that the project can use second-species theory as a
GENERATOR of cycler families (enumerate generating chains, decide which survive at the physical mass, continue them) and so
that the demanded-turn gate (#888) and the symmetric closures (#563, #890) can be cross-checked against an independent
asymptotic statement of the same physics. Companion digests, which this one builds on and does not repeat:
`docs/notes/2026-10-04-digest-marco-niederman-1995-seconde-espece-plan-restreint.md`,
`docs/notes/2026-10-04-digest-bolotin-mackay-2000-second-species-n-centre.md`,
`docs/notes/2026-10-04-digest-font-nunes-simo-2009-second-species-numerical-study.md`,
`docs/notes/2026-10-04-digest-font-nunes-simo-2002-consecutive-quasi-collisions.md`,
`docs/notes/2026-10-04-digest-bradley-russell-2014-patched-conics-to-full-gravity-continuation.md`.

Citations:
- Part I: G. Gomez and M. Olle, "Second-species solutions in the circular and elliptic restricted three-body problem. I.
  Existence and asymptotic approximation", Celestial Mechanics and Dynamical Astronomy 52:107-146 (1991), DOI
  10.1007/BF00049446. 40 pages (journal pages 107-146; PDF page n is journal page n+106). Received 21 February 1990,
  accepted 15 August 1991. Filed in the private paper corpus as
  gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-I-existence-asymptotic-cmda-52-107-doi-10.1007-BF00049446.pdf.
- Part II: same authors, "II. Numerical explorations", CMDA 52:147-166 (1991), DOI 10.1007/BF00049447. 20 pages (PDF page n is
  journal page n+146). Filed in the private paper corpus as
  gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-II-numerical-explorations-cmda-52-147-doi-10.1007-BF00049447.pdf.

Text layer: `pdftotext <file> - | wc -w` gives 13261 words (Part I) and 5772 words (Part II); the layer is OCR of an old scan,
formulas are garbled, so every formula here was read from the page images (all 60 pages), and the printed numerical entries
were re-read from 300 dpi crops (Fig. 8 captions, Fig. 14, 16 and 18 captions, the three starting-orbit lines on p. 160, Table III).

Evidence labels: READ = read from the printed page in this session. COMPUTED = my own short arithmetic, stated with its inputs.
INFERRED = my reading of what the print or the computation implies. Journal page numbers throughout. [unclear] = could not be
read with confidence. Nothing in a table below was computed or filled in by me unless labelled COMPUTED.

## 0. Summary (what the two papers are, and are not)

1. What they ARE. Part I proves, for the planar restricted problem with the secondary on a circular or ELLIPTIC orbit, the
   existence and gives first-order asymptotic formulas (small mass parameter mu) for symmetric periodic second-species solutions
   (SPSSS) whose generating (mu = 0) orbit is a RECTILINEAR Kepler ellipse (circular case) or a NEARLY rectilinear ellipse
   (elliptic case), by Perko's O(mu) matching theory (inner hyperbola near the secondary, outer Keplerian orbit, matching in a
   boundary layer). Part II computes characteristic curves of the families these produce, in the circular problem at
   mu = 1e-6 (families A0, A1, B1, B2, C12) with their bifurcation orbits (Tables I-III), and follows three of them in the
   eccentricity e_m of the primaries from 0 up to about 0.9 (Figs. 13-18).
2. What they are NOT. (a) They do not use Henon's S-arc/T-arc vocabulary or the Font-Nunes-Simo "alphabet"; the generating arcs
   are "orbits with consecutive collisions" (OCC), labelled by families A_j, B_j, C_ij of Gomez-Olle 1986 and by Henon 1968's
   transcendental equation (Part II eq. (1), p. 149). (b) Part II fixes mu = 1e-6 and never computes at an Earth-Moon or
   Sun-Jupiter mass; Guillaume's 1e-2 is mentioned only as a citation (p. 155). (c) Neither part describes a corrector, an initial-
   guess construction for general orbits, or any continuation in mu; "a detailed description of the evolution of all the families
   can be found in [10]" (Part II p. 153, the Olle 1989 thesis, not held). (d) No stability or Floquet result appears anywhere
   (READ: none in either part). (e) Nothing about two secondaries. (f) Neither reference list cites Bruno or Hitzl.
3. Most useful for the project: (i) the exact asymptotic formulas for the flyby (Delta = mu K0, periapsis time shift
   mu ln mu/|v1|^3), which reduce to the gate's own periapsis relation (section 8e: COMPUTED/INFERRED, the strongest original
   content of this digest); (ii) the elliptic result: symmetric periodic solutions in the elliptic problem have period 2 k pi and
   start at pericentre or apocentre of the primaries, which I read as the same kind of structural requirement as the #890 perturber-on-axis condition (INFERRED);
   (iii) printed 16-digit initial conditions at mu = 1e-6 that the project's CR3BP and ER3BP code can reproduce, and which I
   verified are internally consistent (section 9).

## 1. Setting, definitions and the object studied (Part I, pp. 107-111)

### 1.1 Introduction statements (Part I p. 107-108, READ, quoted)

Abstract: "An analytical proof of the existence of some kinds of periodic orbits of second species of Poincare, both in the
Circular and Elliptic Restricted three-body problem, for small values of the mass parameter. The proof uses the asymptotic
approximations for the solutions and the matching theory developed by Breakwell and Perko. In the paper their results are
extended to the Elliptic problem and applied to prove the existence of second-species solutions generated by rectilinear
ellipses in the Circular case and nearly-rectilinear ones in the Elliptic case."

p. 107: "When mu -> 0, the established periodic orbits approach two pieces of keplerian ellipses joined at a corner point. This
type of solutions was called periodic second-species solutions, PSSS, by Poincare [15]." And: "For mu > 0, the PSSS are two
perturbed pieces of ellipse, except when the distance to the small primary becomes small: then, the orbit is a perturbed
hyperbolic arc (arc that reduces to a corner point when mu -> 0)."

p. 107-108: the matching theory (inner solution near the small primary, outer solution away, asymptotic matching) "applied to the
Circular RTBP, as well as the error estimates for this boundary layer approximation were developed by Perko [10], assuming that
the angle between the relative incoming and outgoing velocities, for mu = 0, was different from zero and the close passage from
the infinitesimal body to the small primary was O(mu) (we shall call it O(mu)-matching theory)."

p. 108: "A piece of orbit which begins and ends at a collision point with the small primary is said to be an orbit with
consecutive collisions; and a generating ellipse is simple or n-multiple according to it is formed by one or n arcs of orbits
with consecutive collisions." Guillaume [3-6] extends the matching to bifurcation orbits (relative collision velocities
parallel) and "establishes the existence of families of SPSSS generated by bifurcation orbits, as well as the existence of SSS
with the two orthogonal crossings far from the small primary, but with different close passages to it during a period."
Henrard [8] gives an existence proof "using different tools based on the topological equivalence between differential systems in
the neighbourhood of an equilibrium point."

p. 108, scope statement: "In the Circular RTBP, we give an analytical proof of the existence and the asymptotic approximation of
SPSSS generated by rectilinear ellipses (e0 = 1). This result completes the analytical study of those SPSSS which approach
simple generating ellipses in the Circular RTBP and which have an O(mu)-close passage to the small primary during one period.
In the Elliptic RTBP, we prove analytically the existence and the asymptotic approximation of SPSSS 2k pi-periodic, generated by
nearly-rectilinear orbits with consecutive collisions."

### 1.2 Equations and notation (pp. 108-109, READ)

Equation (1) (p. 109): x'' = -x/|x|^3 - mu [ (x - x_m)/|x - x_m|^3 + x_m/|x_m|^3 - x/|x|^3 ], x(t, mu) in R^2 the position of
the particle P2 in an INERTIAL (sidereal) system centred at the big primary P0; x_m(t) the position of the small primary P1,
which satisfies the two-body problem; mu = m1/(m0 + m1). The orbital eccentricities of P1 and P2 are e_m in [0,1) and e0 <= 1.

Synodic frame (from Part II p. 148, p. 152-153 and the numbers in section 9): P0 at the origin of the sidereal frame; the
rotating (circular) or rotating-pulsating (elliptic) frame has the primaries on the x-axis. "In the whole paper we assume that,
in sidereal coordinates, the big primary P0 is located at the origin and the small one, P1, describes a direct ellipse of
eccentricity e_m with mean angular velocity equal to unity, and for t = 0 it is at the pericenter of its orbit" (Part II p. 148).

Generating ellipse parameters (Part I p. 109): a0, e0, and three signs epsilon, epsilon', epsilon'' "that determine respectively
the sign of the pericenter of its orbit, the sense of motion and the point, pericenter or apocenter, where the motion begins"
(these definitions are restated in Part II p. 149: epsilon = +1 if the pericentre of P2 is positive (on the positive axis),
-1 negative; epsilon' = +1 direct, -1 retrograde; epsilon'' = +1 if P2 begins its motion at the pericentre, -1 at the
apocentre). Equivalent variables (tau, eta, epsilon, epsilon', epsilon''): tau = value of the true anomaly f_m of P1 at the
collision, eta = eccentric anomaly of P2 at the collision (p. 109, p. 113, p. 124).

Condition for periodic, symmetric orbits (p. 109, READ, quoted): "Due to the symmetry of synodical equations (see [16]), we know
that in the Circular RTBP, we have a symmetric periodic orbit if it crosses orthogonally twice the synodical x-axis. While in
the Elliptic RTBP, we must add the condition that both perpendicular crosses occur when the small primary is located at its
pericenter or apocenter; that implies, in particular, that symmetric periodic orbits have a period equal to 2 k pi."

Generating family (p. 109): "the candidates to generating ellipses of SPSSS 2k pi-periodic will be those orbits with
consecutive collisions (OCC) with a collision instant t1 near to k pi." Family x0(t; t1) with t1 in J = [k pi - delta, k pi + delta],
0 < delta < pi/2, "with the restriction eta < pi. That is, we consider families of generating ellipses with rectilinear
solutions (e0 = 1). From the diagram of OCC, see Figure 1 ([1], Figure 6), we take those intervals corresponding to families
A0, A1, B1 with t1 in J and eta < pi. For such solutions we have epsilon'' = -1, i.e., the particle begins its motion at the
apocenter. The restriction eta < pi guarantees that, the collision between the primary P1 and the particle, when mu = 0,
happens before the collision between P0 and P2, that is, x0(t; t1) /= 0, for all t in [0, t1]. On the other hand, we impose
0 < delta < pi/2 since t1 = k pi + pi/2 may correspond to a bifurcation generating orbit."

Fig. 1 (p. 110): the (tau/pi, eta/pi) diagram of families of generating ellipses for e_m = 0.5, near rectilinear (t1 in
[k pi - delta, k pi + delta]); curves labelled A0, A1, A2, A3, B1, B2, B3 with sign triplets (epsilon, epsilon', epsilon''); the
dashed curve = degenerate solutions where P1 and P2 coincide at every moment on the same orbit. Fig. 2 (p. 111): regular
generating ellipse in sidereal and synodic coordinates, and a bifurcation generating ellipse (v1 parallel to v1').

Perko's hypothesis (p. 109): there is t1 > 0 with x0(t1) = x_m(t1), x0(t; t1) /= 0 for t in [0, t1], and the relative incoming and
outgoing collision velocities v1 = x0'(t1) - x_m'(t1), v1' are NOT parallel (|v1| > 0). The sufficient condition printed on
p. 112 (the lower bound M0(t1), formula garbled in the scan [unclear layout]) reduces, p. 114: "for e_m = 0, we recuperate the
condition obtained by Perko ([11], p. 206), |v1| > |epsilon' sqrt(a0 (1 - e0^2)) - 1|."

## 2. The matching construction (Part I, Appendix pp. 142-146, and Sections 2-3)

### 2.1 Outer (Keplerian) solution (Appendix, p. 142-143, READ)

Equations x'' = g(x) + mu f(x, t) with g(x) = -x/|x|^3 and perturbation f(x,t) = -[(x - x_m)/|x - x_m|^3 + x_m/|x_m|^3 - x/|x|^3].
Initial conditions x(0, mu) = x0(0; .) + dx0, x'(0, mu) = x0'(0; .) + dx0'. Outer expansion (A.1), valid for t in
[0, t1 - mu^(1/2) k0], any k0 > 0, t1 in J:

(x, x') = (x0, x0') + mu (x1, x1') + mu^2 (x2s, x2s') + O(mu^2 ln^2 mu, mu^(3/2) ln^2 mu),

with x1 the first-order perturbation (variational equation (A.2) with fundamental matrix Phi(t, t0; .) of the two-body problem
along x0, G(t; t1) = [[0, I2],[g_x(x0(t; t1)), 0]]) driven by the forcing f[x0(t'; t1), t'], and x2s the singular part of the
second order. Estimates (A.3): |x0(t) - x_m(t)| > A1 (t1 - t), |x1(t)| < B1 |ln(t1 - t)|, |x2(t)| < B2 |ln(t1 - t)|/(t1 - t).

### 2.2 Inner (hyperbolic) solution (Appendix p. 143-144, READ)

In y = x - x_m (position relative to the small primary), inner region t in [t1 - mu^(1/2) k1, t1 + mu^(1 - eps) k0] [the exponent
is printed with the glyph epsilon; the lemmas use nu; I read it as nu], k1 = k0 + O(mu^(1/2)):

(y, y') = (y0, y0') + (y1, y1') + O(mu^(2 - nu), mu^(1 - nu)), with y0'' = -mu y0/|y0|^3 (A.5),

a HYPERBOLA about the small primary whose parameters (v_inf, Delta, t_p) are fixed by matching. Estimates (A.6):
|y0| <= mu^(1/2) k_bar, |y1| = O(mu^(2 - nu)), |y1'| = O(mu^(1 - nu)). Fig. 8 (p. 144): the hyperbola with v_inf (asymptote
velocity at t = -infinity), the impact parameter Delta (distance from the small primary to the asymptote) and the pericentre
time t_p.

### 2.3 Matching (A.7)-(A.17), READ

Outer and inner expansions agree to O(mu^(2 - nu)) in the matching region t = t1 - mu^(1/2) k0 + O(mu) if the hyperbola's
constants are (i = v1/|v1| unit vector, j perpendicular, a[i,j] = pi/2):

- (A.7) t_p = t1 - dt1 + O(mu^(2 - nu)), with
- (A.8) |v1| dt1 = i . [Phi_rr(t1,0;.) dx0 + Phi_rv(t1,0;.) dx0'] + mu integral_0^{t1} i . b(t;.) dt + (mu/|v1|^2) [ ln(2 |v1|^3 t1/(mu e1)) - 2 ].
- (A.9) Delta = j . [Phi_rr(t1,0;.) dx0 + Phi_rv(t1,0;.) dx0'] + mu integral_0^{t1} j . b(t;.) dt + O(mu^(2 - nu)).
- (A.10) |v_inf| = |v1| + i . [Phi_vr dx0 + Phi_vv dx0'] - mu/(|v1|^2 t1) + mu integral_0^{t1} i . b1 dt + O(mu^(2 - nu)).
- (A.11) a[v1, v_inf] = (1/|v1|) j . (Phi_vr dx0 + Phi_vv dx0') + (mu/|v1|) j . integral_0^{t1} Phi_vv(t1,0;.) f(x0(t;.), t) dt + O(mu^(2 - nu)).
- (A.12) b(t;.) = Phi_rv(t1, t;.) f(x0(t;.), t) - i/(|v1|^2 (t1 - t)) for t /= t1, b(t1;.) = (0,0)^T; (A.13) b1 the analogous
  expression with Phi_vv and an explicit term -1/|v1|^2 ((1 - 3 cos^2 alpha0), -2 sin alpha0 cos alpha0)^T / (6 |x_m(t1)|^3)
  [layout partly garbled], alpha0 = a[v1, x_m(t1)].
- (A.14) e1 = [1 + (v_inf^2 Delta/mu)^2]^(1/2) (eccentricity of the hyperbola).
- (A.15) |y0(t_p)| = |Delta| ((e1 - 1)/(e1 + 1))^(1/2); (A.16) |y0'(t_p)|^2 = v_inf^2 + 2 mu/|y0(t_p)|.
- (A.17) "t_p = t1 + (mu ln mu)/|v1|^3 = t1 + O(mu^(1 - nu)), for all nu > 0" (dominant term).

Remark on orders (p. 142, 146): Perko proved the boundary-layer approximation is uniformly valid to O(mu^(2 - nu)), all nu > 0, for
all mu in (0, mu1(nu)) and t in [0, t1 + O(mu^(1 - nu))], all t1 in J. "We remark that since Delta = O(mu) and e1 > 1, it follows
|y0(t_p)| = O(mu)" (p. 146).

Relations among hyperbola quantities (Lemma 6, eq. (9), p. 122): r_p = |Delta| ((e1 - 1)/(e1 + 1))^(1/2); v_p = |v_inf| ((e1 + 1)/(e1 - 1))^(1/2);
e1 = [1 + (v_inf^2 |Delta|/mu)^2]^(1/2); |Delta| = r_p v_p/|v_inf|. Under hypothesis O(mu) one has r_p >= mu d1 > 0.

## 3. Existence results with their hypotheses (Part I, Section 2 and 3)

"Hypothesis O(mu)" (p. 115): |dx0| <= mu k0, |dx0'| <= mu k0, mu in (0, mu1), t1 in J = [n pi - delta, n pi + delta], 0 < delta < pi/2.

- Lemma 1 (p. 111): if, in the circular problem, x . x' = 0 and x ^ x_m = 0 at an instant t (in the elliptic problem also
  t = k pi), the orbit crosses the synodic axis orthogonally at t. "The proof is trivial taking into account the change of
  variables from sidereal to the rotating-pulsating (synodical) ones."
- Conditions (4)/(5) (p. 112): a second orthogonal crossing exists at t* iff [x(t*) - x_m(t*)] . [x'(t*) - x_m'(t*)] = 0 (5a), and
  x_m(t*) [x'(t*) - x_m'(t*)] = 0, i.e. relative velocity perpendicular to x_m (5b), and t* = k pi if e_m /= 0 (5c). (5a) means closest
  approach to the small primary (Corollary 5, p. 120: |y| attains its minimum at the unique t* of Lemma 4); (5b) means the
  relative velocity at closest approach is perpendicular to the primaries' line.
- Lemma 2 (p. 114): existence, uniqueness and analyticity of the solution on [0, t1 + mu^(1 - nu) k0] provided the hypothesis (6)
  |j . Phi_rr(t1,0;.) dx0 + j . Phi_rv(t1,0;.) dx0' + mu K1(.)| >= mu d0 holds (the particle does not collide with the small
  primary), K1 = integral_0^{t1} j b dt.
- Lemma 3 (p. 117): an intermediate-value lemma (Bolzano) used to locate zeros of F(t, mu) = F0(t, mu) + O(mu^(2 - nu)).
- Lemma 4 (p. 117): under Lemma 2, for every nu > 0 there is a unique t* in [t1 - mu^(1 - nu) k0, t1 + mu^(1 - nu) k0] satisfying (5a),
  analytic in its variables, and t*(dx0, dx0', mu; t1) = t1 - dt1 + O(mu^(2 - nu)) with
  dt1 = (i/|v1|)[Phi_rr dx0 + Phi_rv dx0'] + (mu/|v1|) integral_0^{t1} i . b dt + (mu/|v1|^3)[ln(2|v1|^3 t1/(mu e1)) - 2]
  (this is (A.8) divided by |v1|, consistent with the Appendix).
- Lemma 6 (p. 121): Delta(dx0, dx0', mu; t1) and |v_inf| are analytic functions of their variables.
- Lemma 7 (pp. 122-123): under Lemma 2, there is Delta = Delta*(dx0, dx0', mu; t1) such that a[x_m(t*), y'(t*, Delta, v_inf, mu, t1)] = pi/2
  (mod pi) (condition (5b)); it is analytic and Delta* = mu [K0 + O(mu^(1 - nu))], with
  K0 = C0 / ( v1^2 sqrt( |x_m(t1)|^2 |v1|^2 - C0^2 ) ), C0 = C' - epsilon' sqrt(a0 (1 - e0^2)), C' = sqrt(1 - e_m^2).
  "If e0 = 1, C0 = sqrt(1 - e_m^2)." K0 well defined and positive requires C0 > 0 and |x_m|^2 |v1|^2 - C0^2 = (x_m . v1)^2 > 0.
  The proof (pp. 126-129) reduces (5b) to the equation arctan( C0 / sqrt(|x_m(t1)|^2 |v1|^2 - C0^2) ) = arccos(1/e1) (p. 127) and
  shows it is solved by Delta = mu K0'.
- Lemma 8 (p. 130): the second initial-condition freedom: if a21(t1) /= 0 there is dx0 = dx0*(dv, mu; t1) of order mu with
  dx0' = (0, dv)^T so that (6) holds and Delta(dx0*, dx0') = Delta*; approximated, in the (i-hat, j-hat) frame of the generating
  ellipse, by dx0* = ( [mu (K0 - K1) - a24 dv]/a21 + O(mu^(2 - nu)), 0 ). If a24 /= 0, symmetrically with dx0 = (dr, 0) and
  dx0'* = (0, [mu (K0 - K1) - a21 dr]/a24 + O(mu^(2 - nu))). Here a_pq are the components of Phi(t1, 0;.) R, R the rotation to
  the (i-hat, j-hat) frame (pp. 130, 133); a21, a24 are the (j-component of x(t1)) derivatives with respect to the i-component
  of x(0) and the j-component of x'(0) respectively [my reading of the index layout, INFERRED].
- Lemma 9 (p. 134): for a rectilinear ellipse with t1 = k pi,
  a[x0(0), v1] = pi + arctan( sqrt(1 - e_m^2) / ( |(-1)^k - e_m| sqrt( 2/(1 - (-1)^k e_m) - 1/a0 ) ) ).
- The nondegeneracy a21(k pi) /= 0 or a24(k pi) /= 0 is verified NUMERICALLY only "for some particular values of k, i.e., k = 1, 2, 3, 4"
  (p. 133), over e_m in [0,1); Fig. 5 (p. 135) plots a21 and a24 against e_m for k = 1..4 (no table; both are nonzero on the
  plotted ranges: a21 positive, growing without bound as e_m -> 1 for k = 2 and 4, small and decaying for k = 1 and 3; a24 negative,
  values between about -6.4 and about 0 across the four curves: approximate, read off a plot, not exact). The generating-ellipse condition
  (p. 133, from [1] eq. (17)): eta solves (1 - (-1)^k e_m)^(3/2) (1 - epsilon epsilon'' (-1)^k cos eta)^(3/2) [eta + epsilon epsilon'' (-1)^k sin eta]
  - k pi sin^3 eta = 0, and for families A0, A1: epsilon epsilon'' = -1, k odd; family B1: epsilon epsilon'' = +1, k even. For e_m = 0 this
  becomes (1 - cos eta)^(3/2) [eta + sin eta] - k pi sin^3 eta = 0 (p. 133).

### 3.1 Theorem 10 (circular problem, p. 136-137, READ, quoted in full)

Define K0 = +/- 1/( |v1|^2 sqrt(|v1|^2 - 1) ), K1 = integral_0^{k pi} j b(t) dt, K2 = integral_0^{k pi} i b(t) dt + (1/|v1|^2)[ln(|v1|^3 2 k pi) - 2],
|v1|^2 = 3 - 1/a0, e1 = [1 + |v1|^4 K0^2]^(1/2) = |v1|/sqrt(|v1|^2 - 1), b(t) as in (A.12).

"Theorem 10. Let x0(t) be a rectilinear generating ellipse with epsilon'' = -1, t1 = k pi and eta < pi such that a21(k pi) /= 0 or
a24(k pi) /= 0. Then, given k0 > 0 and epsilon > 0, there exists mu1 > 0 such that for all mu in (0, mu1) there exists a
uniparametric family of second species solutions of the Circular RTBP, x(t, dx0, dx0', mu), symmetric and periodic in synodical
coordinates, determined by the initial conditions x(0, mu) = x0(0) + dx0, x'(0, mu) = dx0', where dx0' = (0, dv)^T with |dv| <= mu k0
and dx0 = dx0*(dv, mu) defined in Lemma 8 and approximated by dx0* = ([mu (K0 - K1) - a24 dv]/a21 + O(mu^(2 - nu)), 0) with respect to
(i-hat, j-hat); or dx0 = (dr, 0)^T with |dr| <= mu k0 and dx0' = dx0'*(dr, mu) approximated by (0, [mu (K0 - K1) - a21 dr]/a24 + O(mu^(2 - nu)));
The synodical period is T = 2 t*(dx0, dx0', mu) = 2 k pi - 2 dt1(dx0, dx0', mu) + O(mu^(2 - nu)), where for dx0 = (dr, 0)^T and dx0' = (0, dv)^T
dt1(dx0, dx0', mu) = (1/|v1|) ( a11 dr + a14 dv - (mu/|v1|^2) ln(mu e1) + mu K2 )."
[The statement prints "epsilon > 0" where the lemmas use nu > 0; read as nu.] NOTATION WARNING: in Theorem 10 the "synodical
period" T = 2 t* is the FULL period; in Part II, "T: half-period of the SPS" (p. 154). The printed T in Part II Fig. 8 and in the
starting-orbit lines on p. 160 is the HALF period.

Remark p. 137 (READ): "Therefore, for all t in [0, t1] x(t, mu) -> x0(t) when mu -> 0, so, in inertial coordinates and during a
period, the periodic solution approaches two pieces of the rectilinear generating ellipse x0(t) joined at a corner occurring at time
t1; and in synodical or rotating coordinates, the periodic solution approaches an orbit with a single corner at the position of the
small primary at time t1 (see Figure 6)." In the circular case conditions (5a) and (5b) alone guarantee symmetry and periodicity
(p. 136); there is a one-parameter freedom (dv or dr), the family parameter. Hypotheses: rectilinear generating ellipse with
epsilon'' = -1 (start at apocentre), eta < pi, t1 = k pi, a21 or a24 nonzero. mu1 exists but is NOT estimated anywhere (READ: no
numerical bound or sufficient-smallness estimate for mu1 in either part).

### 3.2 The generating-orbit identities I checked (COMPUTED)

For e0 = 1, e_m = 0: the relative speed at the collision is |v1|^2 = (radial speed)^2 + (moon speed 1)^2 = (2 - 1/a0) + 1 = 3 - 1/a0
(printed p. 136); Jacobi constant of the rectilinear generating orbit C = 1/a0 = 3 - |v1|^2 (COMPUTED from printed
ICs, section 9). Hypothesis |v1| > 1 (the condition printed p. 114 with e0 = 1) means a0 > 1/2, i.e. C < 2: Theorem 10 does not cover
rectilinear seeds with C >= 2 (INFERRED from the printed hypothesis). Collision of the particle with P0 before the P1 encounter
is excluded by eta < pi, i.e. t1 < pi a0^(3/2) for the radial orbit; for t1 = k pi this needs a0^(3/2) > k (COMPUTED check on the
printed starting orbits: a0^(3/2) = 1.177, 3.162 (k=3), 2.166 (k=2) respectively).

## 4. The asymptotic formulas, collected (Part I)

All in units of the primaries' separation and period/(2 pi) (GM_total = 1, mean motion of P1 equal to 1).

| quantity | formula (READ) | where |
|---|---|---|
| closest-approach impact parameter | Delta* = mu [K0 + O(mu^(1 - nu))] | Lemma 7, p. 123 |
| K0, circular, rectilinear seed | K0 = 1/( |v1|^2 sqrt(|v1|^2 - 1) ), |v1|^2 = 3 - 1/a0 | Theorem 10, p. 136 |
| K0, general | K0 = C0/( v1^2 sqrt(|x_m(t1)|^2 |v1|^2 - C0^2) ), C0 = sqrt(1 - e_m^2) - epsilon' sqrt(a0(1 - e0^2)) | p. 123 |
| hyperbola eccentricity | e1 = [1 + (v_inf^2 Delta/mu)^2]^(1/2) = |v1|/sqrt(|v1|^2 - 1) (circular, rectilinear) | (A.14), p. 136 |
| closest-approach distance | |y0(t_p)| = |Delta| ((e1 - 1)/(e1 + 1))^(1/2) = O(mu) | (A.15), p. 145-146 |
| periapsis time | t_p = t1 + mu ln mu/|v1|^3 + O(mu) ; t* = t1 + mu ln mu/|v1|^3 + (mu/|v1|){ a14 (K1 - K0)/a24 + ln(e1)/|v1|^2 - K2 } + O(mu^(2 - nu)) | (A.17) p. 146; p. 139 |
| half period (second orthogonal crossing) | t*(t1, mu) as above; synodic full period 2 t* = 2 k pi - 2 dt1 + O(mu^(2-nu)) | p. 139, Theorem 10 |
| initial conditions, circular, rectilinear | x(0) = x0(0) + dx0, x'(0) = dx0'; dx0*, dx0'* as in Theorem 10 (O(mu) corrections, free parameter dv or dr) | Theorem 10 |
| initial conditions, elliptic, p. 141 | see section 5 (Theorem 11) | Theorem 11 |
| error | O(mu^(2 - nu)) on the position, O(mu^(1 - nu)) on the velocity, any nu > 0; relative error on the periapsis distance O(mu^(1 - nu)) | Lemma 2, (A.4), (A.6) |

The closest-approach distance is NOT a free parameter once the generating ellipse is fixed: the matching condition (5b) fixes
Delta and hence r_p (this is the content of Lemma 7). The family parameter (dv or dr) moves the initial point along the
characteristic curve at O(mu) and does not change r_p to leading order (INFERRED from the printed formulas: dv and dr do not
appear in K0 or e1).

## 5. Existence statement for the elliptic problem (Part I, Section 3.2, pp. 137-141)

Problem (p. 137): for every e_m in [0,1) the solutions of Section 2 satisfy (5a) and (5b); one still needs (5c), t* = k pi, so the
solution is 2 k pi-periodic and symmetric. Method (pp. 138-140): choose dr = 0 and use the family parameter t1 itself:
"Each SSS is generated by the ellipse defined by t1, eta(t1), epsilon, epsilon', epsilon'' = -1 or by a0(t1), e0(t1), epsilon, epsilon',
epsilon'' = -1" with initial conditions x(0, mu; t1) = x0(0, t1) + dx0(t1), x'(0, mu; t1) = x0'(0, t1) + dx0'(t1),
x0(0; t1) = (-epsilon a0 (1 + e0), 0)^T, x0'(0; t1) = (0, epsilon' ((1 - e0)/(a0 (1 + e0)))^(1/2))^T, dx0 = (0,0)^T,
dx0' = (0, [mu (K0(t1) - K1(t1))]/a24 + O(mu^(2 - nu)))^T. Then (p. 139)

t*(t1, mu) = t1 - (a14/(|v1| a24)) mu (K0 - K1) + (mu/|v1|^3) ln(mu e1) - mu K2/|v1| + O(mu^(2 - nu))
= t1 + (mu ln mu)/|v1|^3 + (mu/|v1|) { a14 (K1 - K0)/a24 + ln(e1)/|v1|^2 - K2 } + O(mu^(2 - nu)),

so t* = t1 + (mu ln mu)/|v1|^3 + O(mu). Continuity argument (pp. 139-140, Fig. 7): on J = [k pi - mu1^(1/2), k pi + mu1^(1/2)] the map
t1 -> t*(t1) is continuous and |v1| >= M > 0, |mu ln mu| <= mu1*^(1/2) with mu1*^(1/2) = min(mu1^(1/2), (1/4) mu1^(1/2) M^3), so
t*(k pi - mu1^(1/2)) < k pi < t*(k pi + mu1^(1/2)), and by the intermediate value theorem some t1-bar in J gives t* = k pi
exactly (the image interval brackets k pi, Fig. 7).

### Theorem 11 (p. 141, READ, quoted)

"Theorem 11. For all e_m in [0, 1), for all k in N such that a24(k pi) /= 0 (or a21(k pi) /= 0), there exists mu1* > 0 and nu > 0,
0 < mu1* <= mu1, such that for each mu in (0, mu1*), there exists t1-bar in [k pi - mu1^(1/2), k pi + mu1^(1/2)] (which corresponds to a
generating ellipse, defined by a0, e0, epsilon, epsilon'), that yields a SPSSS of the Elliptic RTBP of period T = 2 k pi, and with the
initial conditions approximated by x(0, mu) = (r0, 0)^T, x'(0, mu) = (0, v0 + mu (K0 - K1)/a24 + O(mu^(2 - nu)))^T with respect to the
(i-hat, j-hat) reference system, if a24(k pi) /= 0; or by x(0, mu) = (r0 + mu (K0 - K1)/a21 + O(mu^(2 - nu)), 0)^T, x'(0, mu) = (0, v0)^T
if a21(k pi) /= 0, and where r0 = -epsilon a0 (1 + e0), v0 = epsilon' [ (1 - e0)/(a0 (1 + e0)) ]^(1/2)."

Remarks (p. 141): (1) "When mu -> 0 and for all t in [0, t1], x(t, mu; t1) -> x0(t; t1), i.e., the second-species solution approaches its
generating ellipse." (2) The synodic initial conditions x_s^0 = (x_s^01, 0), x_s'^0 = (0, x_s'^02) satisfy x_s'^02 ~ -x_s^01, "since the
generating ellipses of SPSSS in the Elliptic RTBP are nearly rectilinear; and the synodical initial conditions of a rectilinear
solution are x_s^0 = (r0/(1 - e_m), 0), x_s'^0 = (0, -r0/(1 - e_m))."

Which solutions survive eccentricity: the NEARLY rectilinear ones, t1 within mu^(1/2) of k pi, for every e_m in [0,1), for the k with
a24(k pi) /= 0 or a21(k pi) /= 0 (numerically checked for k = 1..4 only). The resonance condition: the orbit's synodic period is
2 k pi, i.e. k full revolutions of the primaries (INTEGER commensuration with the primaries' period), and both orthogonal
crossings occur with the secondary at pericentre or apocentre (t* = k pi), which is why t1 must be tuned to the narrow window
J. No other generating solutions are treated: general (non-rectilinear) generating ellipses in the elliptic problem, and any
generating chain with more than one close passage per period ("n-multiple"), are outside the theorem (READ: not proven here;
stated as existing in the Guillaume/Perko lineage for the circular problem).

## 6. Part II: what was computed (pp. 147-166)

Scope (Part II abstract, quoted): "The characteristic curves of symmetric periodic second-species solutions in the Circular and
Elliptic RTBP are given for small mu > 0. The behaviour in the neighbourhood of the bifurcation orbits, in the Circular case, is
described." Introduction (p. 147): "The aim of this paper is to obtain numerically the characteristic curves of the symmetric
periodic second-species solutions (SPSSS), both in the Circular and Elliptic planar RTBP and analyse the bifurcations which appear
when doing the continuation of the orbits with consecutive collisions (OCC) for mu = 0 to SPSSS for small values of the mass
parameter. ... We shall 'check' that the regular behaviour is the one predicted analytically: the initial conditions of a SPSSS
differ from those of their generating OCC by O(mu), and the distance from the particle to the small primary at the second
orthogonal crossing (the one close to the primary) is also O(mu). We also show the bifurcations: the characteristic curve near a
bifurcation orbit breaks into two different families whose SPSSS have the following behaviour: the variations in the initial
conditions as well as the distance to the small primary are O(mu^nu), 0 < nu < 1. That result can be found in the analytical researches
done by Guillaume [4], [5], [7] and Perko [12], [13]."

Mass ratios actually computed (READ): mu = 10^-6 only. p. 155: "for mu > 0 (mu = 10^-6), Figure 5 shows the corresponding one of
SPSSS in coordinates (T, LC)"; p. 160 (Section 4): "for any given (but small) value of the mass parameter mu > 0 (in what follows
mu = 10^-6)". Guillaume's mu = 10^-2 appears only as "This behaviour was already described by Guillaume for mu = 10^-2 in [6]"
(p. 155). Elliptic primaries: e_m from 0 to about 0.9 (axis range of Figs. 13, 15, 17), printed orbits at e_m = 0.5. The mu used
for the printed elliptic orbits is not restated at Figs. 14, 16, 18; I take it to be the same 10^-6 (INFERRED from p. 160). The
Earth-Moon (0.01215) and Sun-Jupiter (9.5e-4) values are NOT reached.

### 6.1 Orbits with consecutive collisions (OCC) at mu = 0 (Section 2, pp. 148-149)

Setting: P2 on an elliptic orbit (e, a, major axis on the sidereal X-axis); P and Q the points of the first and second collisions
with P1; "We denote by 2 tau the angle between both collisions and we take the origin of angles in such a way that the collisions
happen at f_m = tau and f_m = -tau". The family of symmetric arcs with consecutive collisions is described by (Gomez-Olle 1986 [2]):

r_m cos tau = epsilon a (epsilon'' cos eta - e)
r_m sin tau = epsilon epsilon' a sqrt(1 - e^2) epsilon'' sin eta
2 arctan( k_m tan(tau/2) ) - e_m sin[ 2 arctan( k_m tan(tau/2) ) ] = a^(3/2) (eta - epsilon'' e sin eta)

with eta the particle's eccentric anomaly at collision, r_m the distance between the primaries, k_m = [(1 - e_m)/(1 + e_m)]^(1/2)
(the third equation is Kepler's equation for P1 at true anomaly tau, equated to the particle's time to collision). Eliminating a
and e gives the single transcendental relation between tau and eta, eq. (1) p. 149 (READ; the right-hand-side factor layout
is the least certain part of the print):

( 2 arctan(k_m tan(tau/2)) - e_m sin[2 arctan(k_m tan(tau/2))] ) |sin eta|^3
 = ( (1 - e_m^2)/(1 + e_m cos tau) )^(3/2) (1 - epsilon epsilon'' cos eta cos tau)^(1/2) [ eta (1 - epsilon epsilon'' cos eta cos tau) - sin eta (cos eta - epsilon epsilon'' cos tau) ].

"For any value of e_m there appear an infinity of families of solutions. As an example, we show in Figure 2 the diagram of solutions
for e_m = 0.5 and tau <= 5.5 pi, eta <= 7 pi. For the particular case e_m = 0, the Circular RTBP, Equation (1) becomes the one obtained
by Henon [8]." Fig. 2 (p. 150): the (tau/pi, eta/pi) diagram, curves A0, A1, A2, A3, A4, B1, B2, B3, B4, B5, loops C23, C24, C25, C26,
C27, C34, C45, C46, C47, C56, C57, with the sign triplets (epsilon, epsilon', epsilon'') on each branch (the sign triplets are
too small to transcribe reliably from the scan [unclear]); dashed curve = degenerate solutions where P1 and P2 follow the same orbit.
Fig. 3 (p. 153), the same diagram at e_m = 0 with the six Henon-Broucke types of SPSSS marked by line style.

### 6.2 Bifurcation orbits and the Henon bifurcation diagram (Section 3, pp. 151-153)

Characteristic curves in the (x, C) plane (x = synodic abscissa of the initial condition, C = Jacobi constant). Periodicity
condition F(x0, C0, mu) = 0; points where two characteristics intersect are "bifurcation points" where the linearisation degenerates.
Guillaume's classification of bifurcation orbits (Part II p. 151, READ): "an OCC is of bifurcation if the incoming and outgoing arcs
at a passage through the small primary match tangentially. He shows there are three basic cases of bifurcation orbits for mu = 0:
case 1, the generating orbit is an ellipse of rational mean motion which intersects the unit circle in two different points;
actually this is a 'composite' generating orbit; case 2, the generating orbit is an ellipse of rational mean motion which is
tangent to the unit circle at a point (we shall say this orbit belongs to S cap E_ij, where S is the family of symmetric OCC and
E_ij the family of second-kind orbits with rational mean motion n = j/i); case 3, the generating orbit is a retrograde circular
orbit of radius one. We denote by I_r the family of circular retrograde orbits, we shall say this orbit belongs to S cap I_r. A case 4
could be defined: the generating orbit is a direct circular orbit of radius one, i.e., at every instant, the orbit coincides with that
of the primaries; but it seems that singular perturbation methods used in cases 1, 2 and 3 do not apply here." Existence: types
(case) 1 in Perko [11] (1977); cases 2 and 3 in Perko [13] (1981); first-species to second-species bifurcation orbits asymptotics
in Guillaume [6], [7]. Near a bifurcation orbit of case 2 or 3 the characteristic "behaves as a hyperbola with the tangents to the
basic characteristics (S, E_ij or I_r) as asymptotes", so it breaks into two families with O(mu^nu) variation.

Henon bifurcation diagram (p. 152): Broucke's principle [1] (quoted): "during the continuous evolution of a family, the relative
position of the orthogonal crossings P' and P'' with respect to the primaries never changes, even when the family encounters
collision orbits", giving six types of symmetric periodic orbits (Types 1-6: six sketches of the order of P', P'' and the primaries
P1, P2 along the x-axis; the printed ordering marks are too small to transcribe reliably [unclear]; Broucke 1968 is held in the
corpus). By this principle "Henon gave a first analysis of simple bifurcation orbits. In a neighbourhood of the synodical x-axis,
once at P' far from P1, and again at P'' near it. So, the generating bifurcation orbit has only one orthogonal crossing at a point P
close to P'; its abscissa given by x(P) = (epsilon epsilon'' + cos tau)/(1 + cos eta) and sign(x(P')) = sign(x(P)) for mu small
enough. On the other hand, the second orthogonal crossing of the SPSSS verifies sign(x(P'') - 1 + mu) = sign(epsilon'' sin eta)."
Fig. 3 (p. 153) draws the Henon bifurcation diagram: all families A_j, B_j, C_ij coloured by type.

### 6.3 Families at mu = 1e-6 (Section 3.3, pp. 153-160)

Notation (p. 153): tau, eta of the OCC at mu = 0; x, ydot synodic initial conditions for t = 0 with (x, 0, 0, ydot) defining a
symmetric periodic solution; C Jacobi integral; T the HALF-period; LC the value of the Levi-Civita coordinate u or v at the second
orthogonal crossing, "since |x + 1 - mu| = u^2 or v^2" (the printed expression conflicts with the other x-conventions of the
paper: see section 11). Remarks (p. 154): (1) "As pointed out by Henrard [9], when mu -> 0, C must be < 3. We just take into account
that the possible motion regions are delimited by simple closed ovals around the primaries when C >= C2(mu) (C2(mu) is the value of
C at the equilibrium solution L2) and these ovals collapse to the primaries when mu -> 0. So, for mu > 0 and small, the condition
C < C2(mu) is required for a second-species solution." (2) "when extending the families of OCC in order to obtain the SPSSS, there
appear collision orbits with T = pi/2 + k pi, k in N (time of collision with the small primary for mu > 0). As noted in [10], the
analytical condition which guarantees that along the extended second-species solution P2 doesn't collide with P1 is that
C0 = |epsilon' sqrt(a0(1 - e0^2)) - 1| /= 0. In particular, for those generating ellipses with epsilon' = +1, tau = pi/2 + k pi, we have
a0 = 1/sin^2 eta, e0 = epsilon'' cos eta and then C0 = 0; so this generating ellipse extends to a collision orbit for mu > 0."
[COMPUTED: for these a0, e0, a0(1 - e0^2) = 1, so C0 = 0; confirms the statement.]

Family A0 (Fig. 4, Fig. 5, Fig. 6, Fig. 8, Table I, pp. 154-157): points of the (tau/pi, eta/pi) curve labelled 1-4 (1 and 3
bifurcation orbits, 2 and 4 extend to collision orbits). Fig. 5: the Levi-Civita distance LC at the second orthogonal crossing
against T for mu = 1e-6: LC changes sign and blows up near each bifurcation orbit (Fig. 5 caption, quoted: "The distance, LC, to the small
primary at the second orthogonal crossing grows when we approach a bifurcation orbit. Around each discontinuity families of different types are
obtained. We label with the same number the extended orbits for mu > 0."). Point 1, tau = eta = pi/2: case 3 (A0 cap I_r); for mu > 0 two branches A0^{mu,1} (type 3: the particle passes on the
left of the small primary at its closest passage) and A0^{mu,2} (type 2: passes on its right), Fig. 6; as tau -> pi/2 the distance
from P2 to P1 at the second orthogonal crossing is O(mu^nu) not O(mu). Point 3, tau = 2 pi, eta = pi: case 3 (also called a first-
species/second-species bifurcation) of A0 cap E21; two branches A0^{mu,3} of type 2 and A0^{mu,4} of type 3 (Fig. 7; the text
calls both points 1 and 3 "case 3-bifurcation orbit of first-species-second-species", and the table calls point 3 case 2; I follow the
table).

Family A1 (Fig. 9, Table II): points 1-7; Family B1 (Fig. 10, Table III), B2 (Fig. 11), C12 (Fig. 12): see below. "A detailed description of
the evolution of all the families can be found in [10]."

Family A1 text (p. 156): "We remark that last point tau = 2 pi, pi < eta < 2 pi, extends to a SPSSS such that, for t = 0, P2 ejects from collision with P0"
(the big primary). Fig. 10 text (p. 158): "We remark point 2, tau = eta = pi. In such orbit, for t = 0, P2 ejects from
collision with P1. Thus, for mu > 0 small, the first orthogonal crossing is very close to P1. Also, we notice that around point 3,
tau = eta = 4.4199..., and in the characteristic curve for mu > 0, two families are obtained: one of type 4 which evolves to orbits
around the Lagrangian equilibrium solution L2, and one of type 6 whose orbits evolve to solutions around L1. Finally, point 9,
tau = pi, pi < eta < 2 pi, extends to an orbit such that, at t = 0, P2 ejects from collision with P0." Family B2 (p. 159): "The
bifurcation point tau = eta = 2 pi, epsilon' = -1 belongs to B2 cap {x = 1}; it is a sidereal circular retrograde orbit of radius one, so the particle
ejects from collision with P1 at tau = 0 and collides again at tau = pi and tau = 2 pi. In synodical coordinates, it is a circular orbit of period
pi which describes two revolutions in the interval [0, 2 pi]. For mu > 0, we obtain two different families of SPS (of type 5 for x > 1 - mu, and
of type 4 for x < 1 - mu) which begin very close to P1." Family C12 (p. 159): "two (1 and 2) bifurcation points in (x, C) coordinates,
which correspond to the double point tau = pi, eta = 2 pi. They belong to C12 cap E12, so they are case 2-bifurcation orbits of first-species-
second-species; and there appear (around each bifurcation point, for mu > 0) two different families of types 2 and 3."

### 6.4 Tables I, II, III (pp. 157, 163, READ, transcribed as printed)

Table I (p. 157). Family A0. "By local behaviour we mean how the SP solutions (in a neighbourhood of the bifurcation or relevant point)
of A0^mu, according to the distance between P1 and P2 at the second orthogonal crossing."

| bifurcation and relevant points | C0 | bifurcation case (Guillaume), intersecting families | type of bifurcated branch (Broucke) | local behaviour |
|---|---|---|---|---|
| tau = eta = pi/2 | /= 0 | 3, A0 cap I_r | 3 if tau < pi/2; 2 if tau > pi/2 | O(mu^nu) |
| tau = 3 pi/2, eta < pi | 0 | regular | - | exists a collision orbit with P1 |
| tau = 2 pi, eta = pi | /= 0 | 2, A0 cap E21 | 2 if tau < 2 pi; 3 if tau > pi/2 [as printed; presumably 2 pi] | O(mu^nu) |
| tau = 5 pi/2, pi < eta < 2 pi | /= 0 [as printed; see below] | regular | - | exists a collision orbit with P1 |

Table II (p. 157). Family A1.

| bifurcation and relevant points | C0 | bifurcation case, families | type | local behaviour |
|---|---|---|---|---|
| tau = eta = 3 pi/2 | /= 0 | 3, A1 cap I_r | 3 if tau < 3 pi/2; 2 if tau > 3 pi/2 | O(mu^nu) |
| tau = 2 pi, eta = pi | /= 0 | 2, A1 cap E21 | 3 if tau < 2 pi; 2 if tau > 2 pi | O(mu^nu) |
| tau = 3 pi, eta = 2 pi | /= 0 | 2, A1 cap E32 | 2 if tau < 3 pi; 3 if tau > 3 pi | O(mu^nu) |
| tau = 5 pi/2, pi < eta < 2 pi | 0 | regular | - | exists a collision orbit with P1 |
| tau = 5 pi/2, 0 < eta < pi | /= 0 | regular | - | O(mu) SPSSS |
| tau = 3 pi/2, pi < eta < 2 pi | /= 0 | regular | - | O(mu) SPSSS |
| tau = 2 pi, pi < eta < 2 pi | /= 0 | regular | - | O(mu) SPSSS |

Table III (p. 163). Family B1. (x = 3.16017... and C = 2.94591... printed on the first row only.)

| bifurcation and relevant points | C0 | bifurcation case, families | type | local behaviour |
|---|---|---|---|---|
| tau = 3 pi, eta = pi; x = 3.16017...; C = 2.94591... | /= 0 | 2, B1 cap E31 | 5 if tau < 3 pi; 6 if tau > 3 pi | O(mu^nu) |
| tau = 3 pi, eta = pi | /= 0 | 2, B1 cap E31 | 6 if tau < 3 pi; 5 if tau > 3 pi | O(mu^nu) |
| tau = eta = pi | /= 0 | 2, B1 cap E11 cap {x = 1} | 4 if tau < pi; 5 if tau > pi | O(mu^nu), t = 0, P2 near P1 |
| tau = eta = 4.4199..., x = 1, C = 3 | /= 0 | 4 | 4 if x < 0.993081... (L2); 6 if x > 1.006948... (L1) | O(mu^nu) |
| tau = pi + pi/2, pi < eta < 2 pi | 0 | regular | - | exists a collision orbit with P1 |
| tau = 2 pi, pi < eta < 2 pi | /= 0 | regular | - | O(mu), P2 collides with P1 before second orthogonal crossing |
| tau = pi + pi/2, 0 < eta < pi | /= 0 | regular | - | O(mu) SPSSS |
| tau = 2 pi + pi/2, 0 < eta < pi | 0 | regular | - | exists a collision orbit with P1 |
| tau = 2 pi + pi/2, pi < eta < 2 pi | /= 0 | regular | - | O(mu) SPSSS |
| tau = pi, pi < eta < 2 pi | /= 0 | regular | - | O(mu) SPSSS |

Figure captions that tie the rows to points (READ): Fig. 4 "Family A0 of OCC. Points 1 and 3 are bifurcation orbits; points 2 and 4
extend to collision orbits"; Fig. 9 "Family A1 of OCC. Points 1, 2 and 3 are bifurcation orbits; point 4 extends to a collision orbit.
Points 5, 6 and 7 extend to SPSSS"; Fig. 10 "Family B1 of OCC. Points 1, 2 and 3 are bifurcation orbits; points 4 and 7 extend to collision
orbits and points 5, 6, 8 and 9 extend to SPSSS". The mapping of Table III rows to the Fig. 10 point numbers is not printed; my mapping
(rows 1-2 = point 1, row 3 = point 2, row 4 = point 3, row 5 = point 4, ..., row 10 = point 9) is INFERRED and fits the p. 158 text
(point 2 is tau = eta = pi, point 3 is tau = eta = 4.4199..., point 9 is tau = pi, pi < eta < 2 pi). Inconsistencies in the printed
tables (kept, not corrected): Table I row 4 prints C0 /= 0 for a collision-orbit point while rows with C0 = 0 are the collision
points elsewhere and p. 155 says points 2 and 4 "extend to collision orbits since C0 = 0"; Table III row 6 prints "SPSSS-regular"
with "P2 collides with P1 before second orthogonal crossing" whereas Fig. 10's caption lists its point as extending to a SPSSS;
Table III labels x < 0.993081 as L2 and x > 1.006948 as L1.

### 6.5 Printed initial conditions at mu = 1e-6 (READ; digits re-read at 300 dpi)

Circular problem, Fig. 8 (p. 159), Levi-Civita-coordinate plots with the following captions (T = half period):

| item | T | x | ydot | C |
|---|---|---|---|---|
| Fig. 8(a): SPSSS of type 2 belonging to A0^{mu,3} | 6.273357979998948 | -2.176283412498242 | 1.63874359406601 | 2.969728009802066 |
| Fig. 8(b): SPSSS of type 3 belonging to A0^{mu,4} | 6.294070698585960 | -2.173709228554383 | 1.634970253275049 | 2.971971477141059 |

Starting orbits followed in eccentricity (Section 4, p. 160, quoted): "So, for any given (but small) value of mu > 0 (in what follows
mu = 10^-6), we pick the initial conditions of those SPSSS in the Circular RTBP, e_m = 0, belonging to families A0^mu, A1^mu, B1^mu and with
half period T = k pi. In particular, we choose from:

| family | x | ydot | T (as printed) |
|---|---|---|---|
| A0^mu | -2.22979002878 | 2.22978812061 | pi |
| A1^mu | -4.30865819030 | 4.30865753635 | 3 pi |
| B1^mu | 3.348331642398 | -3.3483307149 | 2 pi |

and we follow the corresponding family when e_m varies from 0 to 1. We denote by A_{0,1}^{mu,e}, A_{1,3}^{mu,e}, B_{1,2}^{mu,e} those new
families obtained in the Elliptic RTBP, where the second subindex j gives the value of the half-period (T = j pi) of the orbits that belong
to the same family."

Elliptic problem, e_m = 0.5 (pp. 162, 164, 165; READ):

| item | x | ydot |
|---|---|---|
| Fig. 14: SPSSS belonging to A_{0,1}^{mu,e} | -4.898301112912180 | 4.898298973960655 |
| Fig. 16: SPSSS belonging to A_{1,3}^{mu,e} | -8.884269706188229 | 8.884269157978011 |
| Fig. 18: SPSSS belonging to B_{1,2}^{mu,e} | 6.467038249660822 | -6.467037966687732 (second line of the Fig. 18 caption read at page-image resolution, not at 300 dpi) |

No eccentricity-continued orbit has a printed C (there is none in a time-dependent problem), a printed period, or a printed
periapsis distance; the Fig. 13, 15 and 17 curves (e against x, and ydot against x) are plots only: A_{0,1}^{mu,e}: e from 0 to about 0.9 as x runs
from about -2.2 to about -25 (ydot to about +25, nearly linear: ydot ~ -x); A_{1,3}^{mu,e}: x from about -4.3 to about -50;
B_{1,2}^{mu,e}: x from about 3.3 to about 40, e rising steeply then flattening near 0.9, ydot ~ -x (approximate, read off the plots).

Remarks after Fig. 16 (p. 164, READ): "(1) In the diagrams of the characteristic curves in synodical coordinates (x, ydot), we can see that
x ~ ydot. The reason is that the generating ellipse of each SPSSS is nearly rectilinear; and for any generating rectilinear arc, (x0, ydot0),
we have x0 = ydot0. (2) When e_m -> 1, the value of |x| -> infinity. For a given rectilinear generating arc defined by tau = k pi, eta < pi,
epsilon'' = -1, we have in sidereal coordinates that P2 begins its motion at point (-2 epsilon a0, 0); so, in synodical coordinates
x0 = -2 epsilon a0/(1 - e_m) = -2 epsilon (1 + e_m)/(1 + e_m cos tau) . (1 + epsilon cos tau cos eta)/sin^2 eta [layout of the second equality
partly garbled]. For families A0, A1 then epsilon epsilon'' = -1, tau = k pi, k odd, so |x0| = 2 (1 + e_m)/(1 - e_m) . 1/(1 + cos eta) -> infinity when
e_m -> 1 (since 0 < eta < pi). For family B1, then epsilon epsilon'' = 1, tau = k pi, k even, so |x0| = 2/(1 + cos eta) -> infinity when e_m -> 1
(since eta -> pi when e_m -> 1, see [2])." Note that remark (1) writes x0 = ydot0 while every printed pair above has x = -ydot (A0, A1) or
x = -ydot (B1: 3.3483... against -3.3483...): a sign slip in the remark, not in the data (see section 11).

What was NOT done, READ: no stability or Floquet analysis anywhere; no continuation in mu; no case beyond mu = 1e-6; no Earth-Moon or Sun-Jupiter
mass; no statement of where an eccentricity-continued family ends (the plots stop near e_m = 0.9 and the paper only says that |x| -> infinity as
e_m -> 1); no bifurcations of the elliptic families are reported.

## 7. Verification I performed on the printed numbers (COMPUTED)

Inputs: mu = 1e-6, big primary at x = -mu, small primary at x = 1 - mu, Jacobi constant C = 2 Omega - v^2 with
Omega = (x^2 + y^2)/2 + (1 - mu)/r0 + mu/r1 + mu (1 - mu)/2 (the standard Szebehely form, including the additive constant).

| check | result |
|---|---|
| Fig. 8(a) C from printed x, ydot | 2.9697280098020755 against printed 2.969728009802066 (difference 9e-15) |
| Fig. 8(b) C from printed x, ydot | 2.9719714771410595 against printed 2.971971477141059 (difference 5e-16) |
| Fig. 8(a) osculating elements about the big primary (inertial speed = ydot + x) | a = 1.58718, e = 0.37116; Tisserand C = 1/a + 2 sqrt(a(1 - e^2)) = 2.969727; mean motion a^(-3/2) = 0.50011 (close to 1/2, a = 2^(2/3) = 1.58740) |
| Fig. 8(b) osculating semimajor axis (inertial speed = ydot + x) | a = 1.58769; so the two Fig. 8 orbits straddle the E21 value a = 2^(2/3) = 1.58740 (n = 1/2), consistent with their being the two branches around the A0 cap E21 bifurcation orbit |
| E31 point of Table III: a = 3^(2/3) = 2.080084, x = a (1 + e) = 3.16017 gives e = 0.51925; C = 1/a + 2 sqrt(a(1 - e^2)) | 2.945905 against printed 2.94591... |
| A0, A1, B1 starting lines: C = 1/a0 with a0 = |x|/2 | 0.896945, 0.464182, 0.597312 (and C recomputed from the printed x, ydot: 0.896955, 0.464188, 0.597320) |
| same lines: t1 = a0^(3/2) (eta + sin eta) with cos eta = 1/a0 - 1 (Kepler time to the collision, e0 = 1, epsilon'' = -1, e_m = 0), eta/pi | A0: 3.141610 (k pi = 3.141593), eta = 0.5329 pi; A1: 9.424812 (9.424778), eta = 0.6800 pi; B1: 6.283235 (6.283185), eta = 0.6319 pi |
| x + ydot on the starting lines | -1.9e-6, -6.5e-7, +9.3e-7 (the O(mu) departure from exact rectilinearity) |

Reading: the frame is the project's own frame (big primary at -mu, small at 1 - mu), and the paper's C includes the constant
mu (1 - mu): the project's `core.cr3bp.jacobi_constant` (C = x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2 - v^2, no additive constant) returns the printed C minus
mu (1 - mu): 2.96972701 and 2.97197048 for Fig. 8(a) and (b) at mu = 1e-6. At the Earth-Moon mu = 0.01215 the offset is 0.012, large
compared with the C windows quoted below. The three starting orbits generate t1 = k pi to about 2e-5 to 5e-5 (not to 1e-6); the departure
is of the order of the O(mu) initial-condition offsets (x + ydot above) times the sensitivity, which I did not test. The Fig. 8 orbits are
NOT rectilinear; they are the non-rectilinear orbits of the A0 family at the 2:1 bifurcation point (a about 1.587, e about 0.371).

## 8. Techniques applicable to the project's problems

### 8a. Recipe: second-species periodic orbits for one small secondary at a given energy (circular problem), from these papers

1. Fix the energy C < 3 (Henrard's necessary condition for mu -> 0, p. 154; at mu > 0 small the condition is C < C2(mu), the Jacobi value at the paper's L2, which is the collinear point BETWEEN the primaries, p. 154; the paper's C includes the constant mu(1 - mu), so C_L4 = 3 exactly in its convention (COMPUTED: Omega at L4 = 3/2); for Earth-Moon C2 = 3.1883 in the project's convention, 3.2003 in the paper's, COMPUTED with the project's mu = 0.0121506)
   and the period class T = k pi (half period), k = 1, 2, 3, 4 (the only k for which the nondegeneracy a21 or a24 /= 0 is checked, p. 133).
2. Enumerate generating arcs (OCC): the solutions (tau, eta) of Part II eq. (1) (for e_m = 0, Henon's equation) with the sign triplet
   (epsilon, epsilon', epsilon''), giving the families A_j, B_j, C_ij; each OCC carries a0, e0 and, by Tisserand, C = 1/a0 + 2 epsilon'
   sqrt(a0(1 - e0^2)) [the Tisserand form is standard, not printed in these papers; the rectilinear case C = 1/a0 is verified in section 7]. Keep those at
   the target C. For rectilinear seeds the condition is eta < pi, t1 = k pi, epsilon'' = -1, and Theorem 10's |v1| > 1 (C < 2).
3. Admissibility tests (all printed): (i) Perko: v1 and v1' not parallel, C0 = |epsilon' sqrt(a0(1 - e0^2)) - 1| /= 0 (equivalently
   |v1| > |epsilon' sqrt(a0(1 - e0^2)) - 1|); C0 = 0 is the collision-orbit case (tau = pi/2 + k pi, epsilon' = +1); (ii) x0 /= 0 on [0, t1]
   (no collision with the big primary, eta < pi); (iii) not a bifurcation orbit: the OCC must not lie on S cap E_ij, S cap I_r (v1 parallel to v1');
   near these the generating arc is replaced by two families with O(mu^nu) departures (Tables I-III).
4. Initial guess at small mu: x(0) = x0(0) + dx0 with x0(0) = (-epsilon a0 (1 + e0), 0), x0'(0) = (0, epsilon' sqrt((1 - e0)/(a0 (1 + e0)))) and the
   O(mu) shifts of Theorem 10: dx0* = ([mu (K0 - K1) - a24 dv]/a21, 0) or dx0'* = (0, [mu (K0 - K1) - a21 dr]/a24); in practice
   the zeroth-order seed (dx0 = 0) is already within O(mu) and the corrector below absorbs the shift. K0 = 1/(|v1|^2 sqrt(|v1|^2 - 1)); K1 requires a
   quadrature of the Moon's perturbing force along the generating ellipse (b(t), eq. (A.12)); a_pq = components of the two-body state-transition matrix
   over [0, t1] rotated to the (i-hat, j-hat) frame.
5. Correct: the SPSSS is symmetric and periodic iff it crosses the x-axis perpendicularly at t = 0 and again at t* = t1 + mu ln mu/|v1|^3 + O(mu)
   (conditions (5a), (5b)); in the circular case that is all (p. 136). Unknown: x0 at fixed C (ydot0 from the Jacobi constant) with the single residual
   xdot(t*) = 0, or (x0, ydot0) at fixed T. The paper's own characteristic curves are in (x, C).
6. Continue: Part II follows regular families by continuation in the Jacobi constant (characteristic curves, (x, C)); near bifurcation orbits the
   curve splits in two and the step must be refined (variations O(mu^nu)). The papers do NOT describe a mu-continuation; the printed facts that
   bound it are: regular members move by O(mu) in the initial conditions, the closest-approach distance scales as mu K0 and the periapsis time as
   mu ln mu/|v1|^3, and a family ceases to be second-species (a) when C0 -> 0 (the extension is a collision orbit with the small primary, tau = pi/2 + k pi with
   epsilon' = +1), (b) when |v1| -> 1 for rectilinear seeds (K0 -> infinity: the passage becomes tangent), (c) at a bifurcation orbit (E_ij tangency or I_r, where it splits into
   two branches with O(mu^nu) variations), (d) when eta -> pi (collision with the big primary first), and (e) when C reaches C2(mu).
7. At the physical mu check: closest approach r_p = mu (e1 - 1)/|v1|^2 (section 8e) against the secondary's radius plus margin, and that the
   hyperbolic leg stays inside the secondary's Hill sphere (the formulas are first order in mu; the theorem's mu1 is not estimated). VALIDITY TESTS for
   using the first-order formulas at all (COMPUTED from the printed orders, my reading, not a printed criterion): mu |ln mu|/|v1|^3 << 1 (the periapsis time
   shift, (A.17), must be small against the encounter's own time scale) and r_p << mu^(1/2) (the inner region of (A.4), (A.6) has radius mu^(1/2) k_bar). For the
   Earth-Moon rectilinear rows of section 8e the first is 0.010 to 0.035 and r_p/(mu^(1/2)) is at most 0.083 (r_p <= 0.0091 against mu^(1/2) = 0.110), so they pass.

### 8b. What the project already has, and what it lacks

Has (names verified by grep in this session):
- CR3BP frame and Jacobi constant: `core.cr3bp.cr3bp_eom`, `jacobi_constant` (note the additive-constant difference above), `propagate` with STM;
  regularised propagation through the encounter: `core.cr3bp_regularized.propagate_regularized`, `extract_perilune_distance`.
- Symmetric fixed-Jacobi perpendicular-crossing corrector: `search.cr3bp_periodic.correct_symmetric_fixed_jacobi`; `correct_periodic`; Barden stability
  `barden_stability`; Jacobi natural-parameter continuation with the gauntlet: `search.cr3bp_continuation.continue_family`.
- Pseudo-arclength continuation in mu of a perpendicular-crossing cycler: `search.mu_continuation.continue_in_mu` (planar, in (x0, C, mu)); multiple
  shooting `search.cr3bp_multiple_shooting.correct_multiple_shooting`; bifurcation tools `search.bifurcation_detector.floquet_multipliers`,
  `detect_period_multiplying`.
- Two-body resonant seeds `search.jovian_resonant_families.two_body_resonant_seed` and the naive-seeding lineage report
  `search.earth_moon_resonant_families.two_body_seed_lineage_check` (seeds at the physical mu directly); the Casoliva Class 1 table
  `earth_moon_resonant_families.TABLE3_ROWS`.
- Two-body STM: `core.kepler_stm.shepperd_stm` (the matrix needed for a21, a24 over [0, t1]).
- Elliptic: `core.er3bp.er3bp_eom` (state derivatives with respect to f; primaries at -mu and 1 - mu; pulsating frame), `er3bp_stm_eom`, `propagate_er3bp`,
  `search.er3bp_periodic`, `search.er3bp_floquet`, `genome.er3bp_continuation`.
- Gate: `verify.turn_gate.demanded_turn_gate`, `required_periapsis_alt_km`, `available_bend_rad`; `search.two_moon_periodic_890` (symmetric multiple-
  shooting continued in the moon-mass scale; `osculating_flyby`, `flyby_hyperbola`).

Lacks (a grep of `src` and `scripts` for second-species and consecutive-collision terms hits only `search.literature_check`, `search.earth_moon_class1_resonant_connections` and
`search.earth_moon_resonant_families`, which name the topic; I did not read them for an enumerator and found none by their function lists):
1. An enumerator of the generating arcs: solving Part II eq. (1) / Henon's equation over (tau, eta, epsilon, epsilon', epsilon'') and listing (a0, e0, C, family
   label A_j / B_j / C_ij), with the admissibility tests of 8a step 3. The project's #563 symmetric-closure enumeration is the TWO-moon patched-conic analogue;
   there is no one-moon OCC list to compare with.
2. The first-order seed of Theorem 10/11: K0, K1, K2, a_pq. Nothing in the repository computes K1 (the Moon-force quadrature along the Kepler arc) or the
   O(mu ln mu) period/time-of-periapsis shift.
3. A diagnostic that extracts (Delta, v_inf, e1, r_p) from a computed encounter and tests them against the formulas of section 4 (`osculating_flyby` extracts
   osculating elements; it does not compare with the asymptotic prediction).
4. Bifurcation-orbit classification by the Guillaume cases 1-3 (E_ij tangency, I_r) and Henon-Broucke type 1-6 labelling of the crossing order.
5. A symmetric ELLIPTIC corrector that enforces the two conditions of Part I (perpendicular crossing at f = 0 and at f = k pi, secondary at pericentre
   and apocentre/pericentre) and continues in e_m from the circular orbit: `er3bp_periodic` and `genome.er3bp_continuation` exist; whether they impose
   exactly these symmetric conditions was not checked.
6. Continuation in mu from 1e-6 that watches r_p/mu, C0 and |v1| (the printed termination mechanisms): `continue_in_mu` follows the branch and has no
   second-species-specific monitors.

### 8c. Application to the catalogued Earth-Moon cycler families

Catalogue facts (read from `data/catalogue.yaml`, `orbit_elements.cr3bp`): the ross-rt and braik-ross planar rows sit at C = 3.129 to 3.183 (periods
6.1 to 19.4 nondimensional); the braik-ross 3D corridors at C = 2.15 to 3.06; the casoliva rows (casoliva-1-2c ... casoliva-7-3c) at C = 0.489 to 2.763
with periods 6.2832, 12.5664 and 18.85 (2 pi, 4 pi, 6 pi: integer multiples of the primaries' period). The mass ratio is 0.012150584...

Consequences, from these papers' own hypotheses:
- The second-species/matched-asymptotics regime (generating arcs are Kepler arcs meeting the secondary's orbit with relative speed v^2 = 3 - C minus
  an angular-momentum term, i.e. C < 3 for mu -> 0, with "C < C2(mu)" at small mu > 0) contains the casoliva rows and the lower-C braik-ross 3D corridors.
  The ross-rt (1,1), (2,1), (3,1), (3,2), (3,3) rows and braik-ross planar rows with catalogue C = 3.13 to 3.18 (paper convention 3.14 to 3.19, adding mu(1 - mu) = 0.0120,
  assuming the catalogue stores the project convention, not checked) are OUTSIDE the C < 3 limit regime (INFERRED): their generating arcs at their own C would have no
  real encounter with the secondary's orbit, so the O(mu)-matching formulas are not defined for them at that C. They are below the inner collinear value C2 (3.1883 in
  the project's convention), so the neck is open and the secondary is reachable, but that is a finite-mu effect. Whether a family reaches them by continuing in C (and in mu)
  from an admissible family is NOT excluded: Table III row 4 has B1 branches near C = 3 evolving to orbits around the collinear points, and Casoliva's strategy (digest
  `docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md`) is continuation in mu, then in C. That connection has not been tested (INFERRED; open).
- The casoliva Class 1 periods 2 k pi map onto the 2 k pi class of Theorem 10 (circular case), and Part I's starting families have synodic periods
  2 T = 2 pi (A0), 6 pi (A1), 4 pi (B1) and, near the 2:1 point, 4 pi (the Fig. 8 orbits: 2 T = 12.547 and 12.588), so the comparison is direct by period:
  casoliva-2-1a/2-1b (T = 2 pi, C = 0.489, 1.196) against A0 (C = 0.897 at mu = 1e-6, k = 1); casoliva-7-3a/b/c (T = 6 pi, C = 1.02 to 1.07) against A1
  (C = 0.464); the T = 4 pi rows (1-2c, 1-2d, 1-2e, 3-2c: C = 1.57, 2.58, 2.76, 0.709) against B1 (C = 0.597), the Fig. 8 branches (C about 2.97, which
  lies above the 1-2e value 2.76) and the A0 family beyond the 2:1 bifurcation. These pairings are by period only (INFERRED); they have not been
  tested. Casoliva's Class 1 uses a second-species-type seed at mu = 1e-6 continued in mu (digest `docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md`),
  which is the same lineage (Olle is a co-author of both).
- What one would compare: (i) per catalogued casoliva row, the OCC family it descends from and the continuation path in mu (does the row sit on the
  O(mu) regular continuation of the A0/A1/B1/B2/C12 family at the same period, or past a bifurcation orbit); (ii) the printed first-order prediction
  of the closest approach, r_p = mu (e1 - 1)/|v1|^2 (8e), against the row's actual perilune distance (the casoliva-7-3 rows are the tight cyclers; the table
  at 8e gives the scale); (iii) the printed IC offsets, O(mu), against the actual displacement from the mu = 1e-6 member.

### 8d. Two secondaries with different periods: what carries over

From these papers' own hypotheses (READ) the theory has ONE secondary and a symmetric-periodicity argument that rests on a single reflection
symmetry of the rotating (circular) or rotating-pulsating (elliptic) frame (Lemma 1, conditions (4)-(5)). Nothing is said about two secondaries
(READ: no statement). INFERRED consequences:
- Carries over: the local analysis at one encounter. The inner hyperbola, the matching (A.7)-(A.17), the impact parameter Delta = mu K0 and the
  periapsis r_p = mu (e1 - 1)/v_inf^2 and the time shift mu ln mu/|v1|^3 are statements about a single passage by a body of mass mu on a given
  Kepler orbit, with O(mu^(2 - nu)) errors, with no use of the symmetry. In a chain with several moons, each passage can be matched separately,
  provided the passages are separated and the "outer" arcs between them are the unperturbed Kepler arcs; the first-order timing shift per passage
  is the new bookkeeping item. The VALIDITY of the first-order theory for the #890 chain is the point to check, and it fails (COMPUTED): per Titania passage
  mu_T = GM_Titania/GM_Uranus = 3.94e-5, so mu ln mu = -4.0e-4 in units of 1/n_Titania (8.706 d / 2 pi = 1.386 d). The #890 note
  (`docs/notes/2026-10-04-890-titania-oberon-candidate.md`) prints an osculating V-infinity of 0.261 km/s at the Titania flyby; Titania's circular speed is 3.646 km/s,
  so |v1| = 0.0716 and |v1|^3 = 3.67e-4. Then mu |ln mu|/|v1|^3 = 1.09, i.e. about 1.5 days: the periapsis time shift is of order one, not a small correction, and
  the flyby radius (2,126 km, the note's Titania value) is the same size as mu^(1/2) a_Titania = 2,736 km, the edge of the inner region (A.6), not O(mu). The
  one-moon first-order matching is therefore OUTSIDE its range of validity for #890 (a slow, near-tangent encounter, C close to 3), and these formulas can neither
  support nor reject #890's explanation of the phase slip. What does carry over for such low-speed encounters is the structure (hyperbola matched to the outer
  arcs; r_p = mu (e1 - 1)/v_inf^2 as the gate's relation, which is exact for a two-body hyperbola whatever the matching accuracy), not the O(mu) error claims.
- Does not carry over: (i) "symmetric periodic" in a synodic frame: with two moons of different periods there is no common rotating frame; periodicity
  becomes a condition on the pair of phases at the ends of the chain (as in #890, where the perturber must lie on the axis at tau = T). (ii) the one-parameter
  family structure of Theorem 10 (dv or dr free) depends on the circular problem's integral C; the concentric circular four-body model has no Jacobi
  constant. (iii) the tabulated nondegeneracies a21, a24 and the bifurcation tables, which are for one secondary.

### 8e. What the asymptotic formulas predict for the flyby distance at the physical mass, and how this is the gate's own relation

Derivation (all steps from the printed formulas; the identification at the end is mine, labelled COMPUTED/INFERRED):
1. From (A.14) and (A.15): e1 = [1 + (v_inf^2 Delta/mu)^2]^(1/2), r_p = |Delta| sqrt((e1 - 1)/(e1 + 1)). Eliminating Delta: Delta = mu sqrt(e1^2 - 1)/v_inf^2,
   so r_p = mu (e1 - 1)/v_inf^2 (COMPUTED).
2. For a hyperbola with eccentricity e1 the bend (deflection) angle delta satisfies sin(delta/2) = 1/e1 (standard), so e1 - 1 = 1/sin(delta/2) - 1 and
   r_p = (mu/v_inf^2) (1/sin(delta/2) - 1). The project gate states r_p = GM/v^2 (1/sin(demanded/2) - 1) (`verify.turn_gate`, docstring
   "required_alt_km"). With GM = mu in the papers' units (GM_total = 1), v = v_inf = |v1|, the two are THE SAME EQUATION (COMPUTED, 8e step 2).
3. In the papers the turn is not an input: it is the turn the OUTER (Keplerian) geometry demands, and the matching condition (5b), in the form
   arctan( C0/sqrt(|x_m|^2 |v1|^2 - C0^2) ) = arccos(1/e1) (p. 127), is exactly the statement that the hyperbola supplies it. Reading C0 = |x_m x v1| =
   |x_m||v1| sin(theta), theta the angle between the moon's position vector and v1 (p. 123-124 identity |x_m|^2|v1|^2 - C0^2 = (x_m . v1)^2): the left side is
   theta, so cos(theta) = 1/e1 = sin(delta/2), delta = pi - 2 theta, i.e. v_in and v_out are the mirror images of each other about the primaries'
   tangent line, and the demanded turn is delta = pi - 2 theta (COMPUTED/INFERRED from p. 127, 123-124; consistent with (5b), periapsis velocity
   perpendicular to x_m).
   For the circular rectilinear seed: sin(theta) = 1/|v1|, so delta = 2 arccos(1/|v1|) and e1 = |v1|/sqrt(|v1|^2 - 1) (agrees with p. 136; COMPUTED).
   Hence the theory's r_p is the periapsis radius at which an unpowered flyby turns by exactly the demanded angle, i.e. the gate's `required_periapsis_alt_km` (plus the
   body radius) IS the leading-order second-species periapsis, and the gate's pass condition (demanded turn <= bend available at the floor) is the statement that the
   theory's r_p is at least the body radius plus the floor altitude. The asymptotics add one thing the gate does not have: it
   says WHICH r_p the chain produces (r_p = mu K0 sqrt((e1 - 1)/(e1 + 1)), uniquely fixed by the generating arc), so a closure at a different r_p is not on
   the second-species branch.
4. Numbers, circular rectilinear seeds, K0 = 1/(|v1|^2 sqrt(|v1|^2 - 1)), |v1|^2 = 3 - 1/a0 = 3 - C (all COMPUTED; leading order in mu; the theorem's mu1 is
   not estimated, so the extrapolation to the physical mu is NOT a theorem):

   | a0 | C = 1/a0 | |v1|^2 | e1 | turn 2 arccos(1/|v1|) (deg) | r_p/mu | r_p, Earth-Moon (km) | r_p, Sun-Jupiter (km) |
   |---|---|---|---|---|---|---|---|
   | 0.6 | 1.667 | 1.333 | 2.000 | 60.0 | 0.750 | 3503 | 556900 |
   | 1.0 | 1.000 | 2.000 | 1.414 | 90.0 | 0.207 | 967 | 153800 |
   | 2.0 | 0.500 | 2.500 | 1.291 | 101.5 | 0.116 | 544 | 86400 |
   | 5.0 | 0.200 | 2.800 | 1.247 | 106.6 | 0.088 | 412 | 65600 |
   | 100 | 0.010 | 2.990 | 1.226 | 109.3 | 0.0755 | 353 | 56100 |

   (Earth-Moon: mu = 0.0121506, lunar distance 384400 km. Sun-Jupiter: mu = 9.5368e-4, a = 778.57e6 km.) Moon radius 1737.4 km is r_p/mu = 0.372 in these
   units: rectilinear seeds with C < 2 give r_p/mu between 0.0755 (C -> 0) and 0.75 (a0 = 0.6), crossing 0.372 near a0 = 0.73, so at the physical
   Earth-Moon mass the rectilinear second-species solutions with C below about 1.37 pass INSIDE the Moon in this first-order estimate (r_p < 1737 km),
   i.e. they are physical only for a0 between 0.5 and about 0.73 (C between about 1.37 and 2), or after non-first-order corrections. For Jupiter
   (radius 71492 km, threshold r_p/mu = 0.0963) the first-order r_p falls below the planet's radius for a0 above about 3.4 (C below about 0.29; 65600 km at a0 = 5 in the
   table). COMPUTED by hand, labelled estimate.
   Validity of this table (COMPUTED, my reading): the first-order formulas need mu |ln mu|/|v1|^3 << 1 and r_p << mu^(1/2); at Earth-Moon mu the five rows give
   0.035, 0.019, 0.014, 0.011, 0.010 and r_p = 0.0091, 0.0025, 0.0014, 0.0011, 0.0009 against mu^(1/2) = 0.110, so they pass; for #890 they do not (section 8d).
   Regime mismatch to say plainly: the rectilinear seeds have C <= 2 (C = 0.46 to 0.90 for the printed starting orbits) while the catalogue's ross-rt and
   braik-ross planar rows are at C about 3.13 to 3.18; only the Fig. 8 orbits (C about 2.97) lie in that band, and they are not rectilinear, for which the general
   K0 (C0 /= 0 form, p. 123) applies and the table above does not.
5. Predictions for the printed starting orbits at mu = 1e-6 (COMPUTED, leading order; each is a falsifiable check of the project's CR3BP integration):
   | orbit | |v1|^2 | e1 | K0 | demanded turn (deg) | predicted closest approach r_p (units of the separation) | time shift mu ln mu/|v1|^3 |
   |---|---|---|---|---|---|---|
   | A0 start (k = 1) | 2.10305 | 1.38079 | 0.45274 | 92.81 | 1.811e-7 | -4.53e-6 |
   | A1 start (k = 3) | 2.53582 | 1.28496 | 0.31821 | 102.20 | 1.124e-7 | -3.42e-6 |
   | B1 start (k = 2) | 2.40269 | 1.30878 | 0.35142 | 99.65 | 1.285e-7 | -3.71e-6 |
   The sign of the closest approach: sign(x(P'') - 1 + mu) = sign(epsilon'' sin eta) (p. 152) with epsilon'' = -1, sin eta > 0: the particle passes the small
   primary on the side x < 1 - mu at the second orthogonal crossing (READ + INFERRED application).

### 8f. The elliptic part and carrying a circular-problem cycler to the elliptic problem

- Statement (Theorem 11, p. 141; p. 109): for e_m in [0, 1) the symmetric second-species solutions that persist are 2 k pi-periodic in the rotating-pulsating
  frame, start (and end half way) at a pericentre or apocentre configuration of the primaries, and are generated by nearly rectilinear ellipses with t1
  within mu^(1/2) of k pi. The mechanism: the circular-case freedom (dv or dr) is spent on the single extra condition t* = k pi, by tuning the generating ellipse
  (t1 -> t1-bar, Fig. 7). Printed continuations: A_{0,1}, A_{1,3}, B_{1,2} families from e_m = 0 to about 0.9, |x| growing without bound as e_m -> 1 (the pulsating
  frame's unit is the instantaneous separation 1 - e_m at t = 0).
- Consequences for a circular-problem cycler (INFERRED, from the printed period requirement): (i) a circular cycler can be carried to the elliptic problem as
  a PERIODIC orbit only if its synodic period is an integer multiple of 2 pi (k primaries' revolutions), and only if both symmetric crossings can be arranged at
  pericentre and apocentre; the catalogue's casoliva rows (periods 2 pi, 4 pi, 6 pi) qualify; the ross-rt rows do not (periods 10.292, 19.440, 14.788, 17.901, 18.145
  nondimensional = 1.638, 3.094, 2.354, 2.849, 2.888 primary periods), so for them the elliptic perturbation produces tori (the project's quasi-periodic genome),
  not periodic orbits. (ii) The pulsating-frame ICs: the printed initial conditions have ydot = -x (to within O(mu)), which is what a point at rest in the INERTIAL frame
  at pericentre has when the derivative is with respect to the true anomaly f (COMPUTED: d/df of a point fixed in the inertial frame is -J x at pericentre where
  r'(f) = 0); a derivative with respect to time would carry an extra factor sqrt(1 + e_m)/(1 - e_m)^(3/2) = 3.46 at e_m = 0.5. So the paper's ydot is a derivative with
  respect to f (INFERRED); `core.er3bp.er3bp_eom` also uses f. (iii) Take the circular family member at mu = 1e-6, set e_m = 0.5 in steps, and re-solve the
  perpendicular-crossing conditions at f = 0 and f = k pi (a 2 x 2 or 3 x 3 Newton problem); the starting orbits are at pericentre (f = 0 is pericentre, confirmed by p. 148).
  The family parameter lost is the Jacobi constant: the elliptic problem has none (p. 109: only the integer-period symmetric solutions survive).
- Structural parallel to #890 (INFERRED): the elliptic condition (5c) t* = k pi (secondary at the extremal point of its orbit at the second orthogonal
  crossing) is the same kind of requirement as the #890 condition that the perturber sit on the axis at tau = T (2.5 forcing periods); in both cases symmetry
  of the full time-dependent problem needs the time-dependent bodies in a reflection-symmetric configuration at the two crossings.

## 9. Positive controls for `core/cr3bp.py` and `core/er3bp.py`

Common conventions (COMPUTED match, section 7): mu = 1e-6; big primary at (-mu, 0), small primary at (1 - mu, 0); the paper's C = 2 Omega - v^2 includes
the additive constant mu (1 - mu), the project's `jacobi_constant` does not (subtract 9.99999e-7). The state is (x, 0, 0, ydot) on the x-axis; T is the HALF period.

CR3BP (circular), the half period T and x, ydot as printed:
1. Fig. 8(a): x = -2.176283412498242, ydot = 1.63874359406601, half period T = 6.273357979998948, C(paper) = 2.969728009802066 (project convention
   2.96972701). Expect: perpendicular x-axis crossing (y = 0, xdot = 0) at t = T, passing the small primary closely at that crossing, type 2 (particle on
   the right of the small primary).
2. Fig. 8(b): x = -2.173709228554383, ydot = 1.634970253275049, T = 6.294070698585960, C = 2.971971477141059 (project 2.97197048), type 3.
3. A0 start: x = -2.22979002878, ydot = 2.22978812061, T = pi; A1 start: x = -4.30865819030, ydot = 4.30865753635, T = 3 pi; B1 start: x = 3.348331642398,
   ydot = -3.3483307149 (ten decimals as printed, do not pad), T = 2 pi. Expect the second perpendicular crossing at about T (to within the O(mu ln mu) shift,
   about -4e-6, plus the offset tests in section 7 and the sign of the closest approach) at distance about 1.8e-7, 1.1e-7, 1.3e-7 from the small primary
   (8e): a test needs the regularised integrator, since the miss distance is 1e-7 in a unit where the Moon distance is 1.
   Positive-control recipe: integrate with `propagate_regularized`, find the second axis crossing, compare xdot there with 0, the crossing time with T, the
   closest approach with the 8e prediction (leading order; the sensitivity a21, a24 of order 1 to 6 means 11-digit ICs reproduce the miss distance to a few
   percent at best, INFERRED).
ER3BP (elliptic, e_m = 0.5, mu = 1e-6 assumed), independent variable f (true anomaly), pericentre at f = 0, derivatives with respect to f:
4. Fig. 14 (A_{0,1}^{mu,e}): x = -4.898301112912180, y' = 4.898298973960655, expect symmetric perpendicular crossing at f = pi and period 2 pi.
5. Fig. 16 (A_{1,3}^{mu,e}): x = -8.884269706188229, y' = 8.884269157978011, crossing at f = 3 pi, period 6 pi.
6. Fig. 18 (B_{1,2}^{mu,e}): x = 6.467038249660822, y' = -6.467037966687732, crossing at f = 2 pi, period 4 pi.
   For each, the first-order estimate of the generating ellipse: a0 = |x| (1 - e_m)/2 = 1.2246, 2.2211, 1.6168 (COMPUTED, rectilinear approximation).
   Success criterion: the second crossing must be perpendicular (y = 0, x' = 0) and the closest approach to the small primary O(mu) = O(1e-6) in the pulsating frame, with
   the primaries at f = k pi at their pericentre or apocentre.

Best candidates: the Fig. 8 orbits (16-digit state, an independent published C to 1e-14, a half period, both branch types), then the three starting lines (the
generating-orbit identity t1 = k pi in section 7 is an independent check), then the elliptic Figs. 14, 16, 18.

## 10. What the papers leave open (quoted)

- p. 133 (Part I): the nondegeneracy hypothesis is established only numerically and only for k = 1..4: "These computations have been done numerically for some
  particular values of k, i.e., k = 1, 2, 3, 4."
- p. 108 (Part I): scope limit, "SPSSS which approach simple generating ellipses in the Circular RTBP and which have an O(mu)-close passage to the small
  primary during one period"; n-multiple generating ellipses and orbits with two orthogonal crossings far from the small primary are cited to Guillaume, not proven here.
- p. 151 (Part II): "A case 4 could be defined: the generating orbit is a direct circular orbit of radius one ... but it seems that singular perturbation methods used in
  cases 1, 2 and 3 do not apply here."
- p. 153 (Part II): "A detailed description of the evolution of all the families can be found in [10]" (the thesis); only family A0 is described in detail here, the others by tables.
- Not addressed (READ: absent): stability, mu values above 1e-6, the size of mu1, a corrector or continuation algorithm, the eccentricity at which the elliptic
  families end, any elliptic bifurcation.

## 11. Reliability notes and internal inconsistencies of the print

1. Frame: Part II's x conventions disagree among themselves: p. 154 "|x + 1 - mu| = u^2 or v^2" places the small primary at -(1 - mu); p. 152 ("sign(x(P'') - 1 + mu)"),
   p. 159 ("x > 1 - mu", "x < 1 - mu") and p. 163 (x = 1, C = 3 near the small primary) place it at +(1 - mu). My Jacobi recomputation (section 7) reproduces the printed C
   to 1e-14 only with the small primary at +(1 - mu) and the big at -mu; Font-Nunes-Simo 2009 place the small primary at the opposite side (mu - 1), so ICs from
   that paper cannot be mixed with these without a reflection x -> -x (y -> -y).
2. Lagrange labels: Table III prints (L2) for x < 0.993081 (the point between the primaries) and (L1) for x > 1.006948 (beyond the small primary); the p. 154 remark
   (ovals around the primaries for C >= C2) is consistent with C2 being the value at the point between the primaries, so the paper uses inner = L2, outer = L1
   consistently. This is a naming convention, not an inconsistency; the project's own labelling (check `search.periapse_map.collinear_l1_l2`, `search.binary_star_search.collinear_lpoints`)
   may differ, so use positions (0.99309 and 1.00695 at mu = 1e-6) not labels.
3. Remark (1) p. 164 states x0 = ydot0; the data have x = -ydot. Sign slip in the remark.
4. T notation: Theorem 10 p. 137 uses T = 2 t* (full period); Part II uses T for the half period.
5. Table I row 4 prints C0 /= 0 for a collision-orbit point (text says C0 = 0); Table III row 6 and Fig. 10's caption disagree on whether the point extends to a SPSSS; Table I row 3 prints
   "3 if tau > pi/2" where 2 pi is presumably meant.
6. The statement of Theorem 10 and Lemma 2 prints "epsilon > 0" where nu > 0 is used in the proofs; the inner-region exponent in the Appendix prints the glyph epsilon for nu.
7. Reference years: Part I [11] prints "Perko 1974 ... SIAM J. Appl. Math. 41, 200-237" (volume 41 is 1981; the year is probably a slip, INFERRED).
8. The reproducibility of my section 8e/9 predictions rests on the first-order theory; the actual error terms are O(mu^(2 - nu)) with unestimated constants and unestimated mu1.

## 12. References on second-species solutions and generating orbits, as printed

Corpus status: "held" = a file is in the corpus directory (checked by listing); "digested" = a digest exists. The corpus index (`docs/notes/CORPUS_INDEX.md`) was grepped; no row exists for the
works marked "not held". The digests in `docs/notes/` cite Perko and Guillaume only by reference.

Part I list (pp. 146):
- [1] G. Gomez and M. Olle, 1986, "A Note on the Elliptic Restricted Three-Body Problem", Celest. Mech. 39, 33-55. Not held.
- [2] G. Gomez and M. Olle, "Second-Species Solutions in the Circular and Elliptic Restricted Three-Body Problem II: Numerical Explorations", Celest. Mech. 52, 147-166 (this issue). Held (this digest).
- [3] P. Guillaume, 1973, "Periodic Symmetric Solutions of the Restricted Problem", Celest. Mech. 8, 199-206. Not held.
- [4] P. Guillaume, 1973, "A Linear Description of the Second-Species Solutions", in B.D. Tapley and V.G. Szebehely (eds.), Recent Advances in Dynamical Astronomy, 161-174. Not held.
- [5] P. Guillaume, 1975, "Linear Analysis of One Type of Second-Species Solutions", Celest. Mech. 11, 213-254. Not held.
- [6] P. Guillaume, 1975, "An Extension of Breakwell-Perko's Matching Theory", Celest. Mech. 11, 449-468. Not held.
- [7] M. Henon, 1968, "Sur les Orbites Interplanetaires qui Rencontrent Deux Fois la Terre", Bull. Astron., Serie 3, 337-402. Not held.
- [8] J. Henrard, 1980, "On Poincare's Second-Species Solutions", Celest. Mech. 21, 83-97. Not held.
- [9] M. Olle, 1989, "Solucions Periodiques de Segona Especie en el Problema Restringit de Tres Cossos", Doctoral dissertation, Universitat Autonoma de Barcelona. Not held.
- [10] L.M. Perko, 1964, "Asymptotic Matching in the Restricted Three-Body Problem", Doctoral dissertation, Stanford University. Not held.
- [11] L.M. Perko, 1974 [as printed], "Periodic Orbits in the Restricted Problem: Existence and Asymptotic Approximation", SIAM J. Appl. Math. 41, 200-237. Not held.
- [12] L.M. Perko, 1976, "Second-Species Periodic Solutions with an O(mu) Near-Moon Passage", Celest. Mech. 14, 395-427. Not held.
- [13] L.M. Perko, 1981, "Periodic Orbits in the Restricted Problem: An Analysis in the Neighbourhood of a 1st Species-2nd Species Bifurcation", SIAM J. Appl. Math. 41, 181-202. Not held.
- [14] L.M. Perko, 1981, "Second-Species Solutions with an O(mu^nu) 0 < nu < 1 Near-Moon Passage", Celest. Mech. 24, 155-171. Not held.
- [15] H. Poincare, 1899, Les Methodes Nouvelles de la Mecanique Celeste, Tome III, Gauthier-Villars, Paris. Not held.
- [16] V. Szebehely, 1967, Theory of Orbits, Academic Press, New York. Held (file present in the corpus directory).

Part II list (p. 166):
- [1] R.A. Broucke, 1968, "Periodic Orbits in the Restricted Three-Body Problem with Earth-Moon Masses", JPL Tech. Report 32-1168. Held and digested (`docs/notes/2026-07-28-744-broucke-leiva-barrabes-earth-moon-lineage-digest.md`).
- [2] G. Gomez and M. Olle, 1986, "A Note on the Elliptic Restricted Three-Body Problem", Celest. Mech. 39, 33-55. Not held.
- [3] G. Gomez and M. Olle, Part I (this issue). Held.
- [4] P. Guillaume, 1971, "Solutions Periodiques Simetriques du Probleme Restreint des Trois Corps pour de Faibles Valeurs du Rapport des Masses", Ph.D. dissertation, Universite de Liege. Not held.
- [5] P. Guillaume, 1973, "A Linear Description of the Second-Species Solutions", in Tapley and Szebehely (eds.), Recent Advances in Dynamical Astronomy, 161-174. Not held.
- [6] P. Guillaume, 1973, "Periodic Symmetric Solutions of the Restricted Problem", Celest. Mech. 8, 199-206. Not held.
- [7] P. Guillaume, 1975, "An Extension of Breakwell-Perko's Matching Theory", Celest. Mech. 11, 449-468. Not held.
- [8] M. Henon, 1968, "Sur les Orbites Interplanetaires qui Rencontrent Deux Fois la Terre", Bull. Astron., Serie 3, 337-402. Not held.
- [9] J. Henrard, 1980, "On Poincare's Second-Species Solutions", Celest. Mech. 21, 83-97. Not held.
- [10] M. Olle, 1989, doctoral dissertation (as Part I [9]). Not held.
- [11] L.M. Perko, 1977, "Second-Species Solutions with an O(mu^nu), 1/3 < nu < 1 Near-Moon Passage", Celest. Mech. 16, 275-290. Not held.
- [12] L.M. Perko, 1981, "Second-Species Solutions with an O(mu^nu), 0 < nu < 1 Near-Moon Passage", Celest. Mech. 24, 155-171. Not held.
- [13] L.M. Perko, 1981, "Periodic Orbits in the Restricted Problem: An Analysis in the Neighbourhood of a First-Species-Second-Species Bifurcation", SIAM J. Appl. Math. 41, 181-202. Not held.

Neither list cites Bruno, Hitzl, or the Hitzl-Henon work, nor Henon's later book "Generating Families in the Restricted Three-Body Problem" (Springer 1997; cited by Casoliva et al., listed in the acquisition backlog, not held).
The papers' lineage for the project: Guillaume and Perko (existence by matching) to Gomez-Olle (circular, elliptic, numerics) to Barrabes-Gomez 2002/2003 (in/out maps; held, digested) to Casoliva et al. 2008/2010 (Earth-Moon, held, digested).

## 13. Proposed `CORPUS_INDEX.md` rows (for the coordinator)

- gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-I-existence-asymptotic-cmda-52-107-doi-10.1007-BF00049446.pdf | 2026-10-04-digest-gomez-olle-1991-second-species-circular-elliptic-I-II.md | Part I: existence and first-order asymptotics (Perko O(mu) matching) of symmetric periodic second-species solutions generated by rectilinear ellipses (circular) and nearly rectilinear ones (elliptic, period 2 k pi); Delta = mu K0, periapsis time shift mu ln mu/|v1|^3; matching formulas reduce to the gate's periapsis relation | digested | text-layer (13261 words; OCR, formulas read from images)
- gomez-olle-1991-second-species-solutions-circular-elliptic-restricted-three-body-problem-II-numerical-explorations-cmda-52-147-doi-10.1007-BF00049447.pdf | 2026-10-04-digest-gomez-olle-1991-second-species-circular-elliptic-I-II.md | Part II: characteristic curves of second-species families A0, A1, B1, B2, C12 at mu = 1e-6 with bifurcation tables, printed 16-digit ICs (circular Fig. 8; elliptic e_m = 0.5 Figs. 14, 16, 18); no corrector or mu-continuation described | digested | text-layer (5772 words; OCR)
