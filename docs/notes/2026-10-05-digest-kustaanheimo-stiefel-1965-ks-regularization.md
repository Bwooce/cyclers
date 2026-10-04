# Digest: Kustaanheimo and Stiefel 1965, "Perturbation theory of Kepler motion based on spinor regularization"

P. Kustaanheimo and E. Stiefel, "Perturbation theory of Kepler motion based on spinor regularization", Journal fuer die reine und angewandte Mathematik 218:204-219 (1965), DOI 10.1515/crll.1965.218.204 (received 24 July 1964; written at the ETH Zurich in summer 1964; IBM sponsored). Filed in the private paper corpus as `kustaanheimo-stiefel-1965-perturbation-theory-kepler-motion-spinor-regularization-crelle-218-204-doi-10.1515-crll.1965.218.204.pdf` (16 PDF pages; PDF page n is journal page 203 + n, PDF 1 is the title/abstract page and the paper starts on journal p.204; a scan with an OCR layer whose equations are mostly lost).
Which pages I read: the text layer of all pages for prose, and page IMAGES (110 dpi) of PDF 2 to 10 and 15 and 16 (journal pp.205 to 213, 218, 219), which hold every equation from (4) to (53) and (70), (71) and the end matter. Equations (54) to (69) (journal pp.214 to 217: Stumpff form of Kepler's equation, first-order perturbation integrals, the non-conservative time law) were read on the text layer only, where they are partly garbled, and are summarised, not transcribed, and marked so. No digit in a transcribed equation was illegible.
Evidence tags: READ (journal page) is what the page image shows; COMPUTED is my own check of 2026-10-05 (scratch, not committed); INFERRED is my reading across sources.
Companions: `docs/notes/2026-10-04-digest-stiefel-scheifele-1971-linear-regular-celestial-mechanics.md`, `docs/notes/2026-10-04-digest-aarseth-zare-1974-regularization-three-body.md`, `docs/notes/2026-10-04-digest-peters-1968-numerical-regularization.md`, `docs/notes/2026-10-04-digest-aarseth-1971-direct-integration-n-body.md`.

## 0. What the paper is

The original KS paper. Abstract (p.204, READ): a regularisation of Kepler motion in R3 built from a simple map of R4 onto R3; in R4 the equations of any undisturbed Kepler motion are linear with constant coefficients and so stay regular at the centre; no knowledge of spinors is needed; the stated drawback is that the regularised methods need more integrations. It prints NO numerical example: the last sentence (p.219) says an account of numerical experiments by M. Roessler "will appear in ZAMP" (not held, not identified). It prints no table and no figure. Its value for #928 is as the reference for the KS map, the bilinear condition, the generalised forces and the regularised equations in the authors' own signs and factors.

## 1. Levi-Civita and the Hurwitz argument (pp.204 to 206, READ)

Levi-Civita: z = w^2 (1), x1 = u1^2 - u2^2, x2 = 2 u1 u2; dx1 = 2(u1 du1 - u2 du2), dx2 = 2(u2 du1 + u1 du2) (2); the 2 x 2 matrix [[u1, -u2], [u2, u1]] has linear elements (property 1) and orthogonal rows of norm u1^2 + u2^2 (property 2), giving dx1^2 + dx2^2 = 4(u1^2 + u2^2)(du1^2 + du2^2) (4). Generalisation to n dimensions needs an n x n matrix with linearly homogeneous elements and orthogonal rows of norm sum u_i^2: by Hurwitz this exists only for n = 1, 2, 4, 8 ("the case n = 3 is out of question; therefore there is no easy generalization of Levi-Civita's transformation in R3"), so they take n = 4. The matrix (READ p.205), "intimately connected with the rule of multiplication of quaternions":
 A = [[u1, -u2, -u3, u4], [u2, u1, -u4, -u3], [u3, u4, u1, u2], [u4, -u3, u2, -u1]].
Multiplying A into the column du gives four forms (5); "a disappointment": only the first three are complete differentials, giving by integration
 x1 = u1^2 - u2^2 - u3^2 + u4^2, x2 = 2(u1 u2 - u3 u4), x3 = 2(u1 u3 + u2 u4) (6) (READ p.205),
the Kustaanheimo spinor map. The fourth form is annihilated by imposing, on the differentials, the side condition (READ p.206)
 u4 du1 - u3 du2 + u2 du3 - u1 du4 = 0 (7),
after which dx = 2 A3 du with the 3 x 4 matrix A3 = the first three rows of A, and the last row 0 (8); and x = A u with the fourth component identically 0 (9). "We did not succeed in defining a Levi-Civita transformation mapping two 4-dimensional spaces onto each other. But we found a mapping of R4 onto R3 having the desired behaviour."

## 2. Geometry (pp.206 to 208, READ)

- The (u1, u2) plane maps to the (x1, x2) plane by Levi-Civita's rules: the map continues Levi-Civita's.
- dx1^2 + dx2^2 + dx3^2 = 4(u1^2 + u2^2 + u3^2 + u4^2)(du1^2 + ... + du4^2) (10); x1^2 + x2^2 + x3^2 = (sum u_i^2)^2 (11), r = u1^2 + u2^2 + u3^2 + u4^2 (12) (r is the square of the distance in R4); hence dx^2 = 4 r du^2 (13), valid only when (7) holds.
- Inverse and kernel (14): all points v with v1 = u1 cos(phi) - u4 sin(phi), v2 = u2 cos(phi) + u3 sin(phi), v4 = u1 sin(phi) + u4 cos(phi), v3 = -u2 sin(phi) + u3 cos(phi) map to the same x; the preimage is a circle of radius sqrt(r) with tangent (-u4, u3, -u2, u1) (15), so (7) says du is orthogonal to that circle; the rotation group (14) is the "kernel". This is the fibre, with the same signs as Stiefel and Scheifele chapter XI (43,7), (43,8) and Waldvogel (7).
- Inverse map (16) (READ p.207): u1^2 + u4^2 = (x1 + r)/2, u2 = (x2 u1 + x3 u4)/(x1 + r), u3 = (x3 u1 - x2 u4)/(x1 + r); alternative for negative x1: u2^2 + u3^2 = (-x1 + r)/2, u1 = (x2 u2 + x3 u3)/(-x1 + r), u4 = (x3 u2 - x2 u3)/(-x1 + r). "The formulae in the first line determine two of the u_i up to a rotation (14); the remaining u_k follow from the formulae below."
- Differential inverse (17) (READ p.207): du = (1/(2r)) A3^T dx, i.e. rows (u1, u2, u3), (-u2, u1, u4), (-u3, -u4, u1), (u4, -u3, u2) acting on (dx1, dx2, dx3); (7) is then satisfied.
- Planes (18) to (21): for two vectors u, v with u4 v1 - u3 v2 + u2 v3 - u1 v4 = 0 (18) (the bilinear relation, in the same sign as (7), and as the S and S (9,34)) the plane they span maps conformally onto a plane in Levi-Civita style; image x = rho^2 [a cos(2 theta) + (b cos(omega) + c sin(omega)) sin(2 theta)] (20), a, b, c the columns of the rotation matrix (21) (a Cayley parametrisation); conics centred at the origin of such a plane go to conics with a focus at the origin, a straight line to a parabola.

## 3. The equations of motion (pp.208 to 212, READ)

Particle of mass m, force P_k: m xddot_k = P_k (22). The corresponding motion u_i(t) must satisfy the non-holonomic condition u4 udot1 - u3 udot2 + u2 udot3 - u1 udot4 = 0 (23). Lagrange: v^2 = 4 r (udot1^2 + ... + udot4^2) (24), r = sum u_i^2 (25), T = (1/2) m v^2 = 2 m r sum udot_i^2 (26); generalised forces from sum P_k dx_k = sum Q_i du_i (27) give, by the transpose of (8),
 Q = 2 A3^T P: (Q1, Q2, Q3, Q4)^T = 2 [[u1, u2, u3], [-u2, u1, u4], [-u3, -u4, u1], [u4, -u3, u2]] (P1, P2, P3)^T (28) (READ p.209).
Lagrange's equations (29) give 4 m [ d/dt(r udot_i) - u_i sum udot_j^2 ] = Q_i (30). Initial values from x0 by (16) and velocity (READ p.209)
 udot1^0 = (1/(2 r0))(u1 xdot1 + u2 xdot2 + u3 xdot3), udot2^0 = (1/(2 r0))(-u2 xdot1 + u1 xdot2 + u4 xdot3), udot3^0 = (1/(2 r0))(-u3 xdot1 - u4 xdot2 + u1 xdot3), udot4^0 = (1/(2 r0))(u4 xdot1 - u3 xdot2 + u2 xdot3) (31),
which satisfies (23) at t = 0; multiplying (30) by (u4, -u3, u2, -u1) and adding shows r(u4 udot1 - u3 udot2 + u2 udot3 - u1 udot4) is constant (using (28)), hence zero for all t: the side condition is a first integral (their Theorem 6 equivalent in the book). "The final steps of the proof are left to the reader."
Regularisation (pp.210 to 212): the work W (32), the kinetic energy equation 2 m r sum udot^2 - (1/2) m v0^2 = W (33) (v0 the initial velocity) turns (30) into 4 m d/dt(r udot_i) - (1/r)(2W + m v0^2) u_i = Q_i (34); regularising time s = integral dt/r, d/dt = (1/r) d/ds (35), and
 4 m d^2u_i/ds^2 - (2W + m v0^2) u_i = r Q_i (36),
regular at the origin. Initial velocities in s (READ p.211): du1/ds = (1/2)(u1 xdot1 + u2 xdot2 + u3 xdot3), du2/ds = (1/2)(-u2 xdot1 + u1 xdot2 + u4 xdot3), du3/ds = (1/2)(-u3 xdot1 - u4 xdot2 + u1 xdot3), du4/ds = (1/2)(u4 xdot1 - u3 xdot2 + u2 xdot3) (37), that is du/ds = (1/2) A3^T xdot. Conservative case: W = V0 - V and Q_i = -dV/du_i (38). Non-conservative: restrict to the time change alone (39), (40): m [4 d^2u_i/ds^2 - v^2 u_i] = r Q_i with v^2 = (4/r) sum (du_j/ds)^2 (41).
Kepler plus perturbation (section 6, READ p.212): central potential -mM/r gives Q_i = -(2 m M/r^2) u_i + Q_i' (42), W = m M (1/r - 1/r0) + W' (43), and
 m [ 4 d^2u_i/ds^2 + (2M/r0 - v0^2) u_i ] = r Q_i' + 2 W' u_i (44), 1/a0 = 2/r0 - v0^2/M (45) (a0 the osculating semi-major axis at s = 0),
 d^2u_i/ds^2 + (M/(4 a0)) u_i = (1/(4m))(r Q_i' + 2 W' u_i) (46),
"completely regular at the origin; even if a trajectory collides with the central mass it is computable". The non-conservative form (47): m[4 u'' + (2M/r - v^2) u] = r Q_i, whose coefficient of u is not constant.

## 4. The initial value problem and first-order perturbations (pp.213 to 217)

Unperturbed motion (READ p.213): d^2u/ds^2 + (M/(4 a0)) u = 0 (48), "the image point moves as if connected with the origin by an elastic string of rigidity M/4a0", path a conic centred at the origin; for a0 > 0 a frequency omega^2 = M/(4 a0) (49), u_i = alpha_i cos(omega s) + beta_i sin(omega s) (50), alpha = u at s = 0, omega beta = du/ds there; sum alpha^2 = r0, sum beta^2 = r0 v0^2/(4 omega^2) (51), sum alpha_i beta_i = (r0 v0/(2 omega)) cos(theta) (52) (theta the angle between the initial position and velocity), r(s) = r0 cos^2(omega s) + (r0 v0^2/(4 omega^2)) sin^2(omega s) + (r0 v0/omega) cos(theta) sin(omega s) cos(omega s), and the Kepler time (53) t = r0 s - (r0/(4 omega))(1 - v0^2/(4 omega^2))(2 omega s - sin(2 omega s)) + (r0 v0/(4 omega^2)) cos(theta)(1 - cos(2 omega s)); the Stumpff form is "Stumpff's principal equation" (54) (text layer), with lambda = 2 omega s the eccentric-anomaly increment, to be solved by Newton's method. The paper notes the R4 has served only "as a tool for establishing that equation".
Perturbations (text layer for (55) to (65)): write the right side of (46) as F_i (55), (56) a forced oscillator; the particular integral is (1/omega) integral of F_i sin(omega(s - sigma)) d sigma (57); the perturbed motion is the Kepler motion plus Delta u_i = (1/omega) integral F_i sin(omega(s - sigma)) d sigma (59), velocity perturbation (60), element perturbations Delta alpha_i = -(1/omega) integral F_i sin(omega sigma) d sigma, Delta beta_i = (1/omega^2) integral F_i cos(omega sigma) d sigma-type forms (61) to (63), "alpha_i, beta_i as elements of the motion in the parametric space"; the physical time is perturbed: Delta t = integral Delta r ds (64) with Delta r = 2 (u + u_K) . Delta u (the (u_j^2 - u_jK^2) expression), so time must be integrated along (65). General perturbations: Fourier-expand F in s (eccentric anomaly), Hansen's technique; planar basic motion (u1, u2 plane) needs no u3, u4 for the perturbations of the plane motion; parabolic and hyperbolic modifications "are not too difficult". These sections (pp.214 to 217) were read on the text layer only.

## 5. The non-conservative case, constant-frequency time (pp.217 to 219; READ p.218 and 219)

With the perturbing force per unit mass written q_i = (1/m)(...) (66) and a new time s1 with s = f(s1), the equations (67) have the coefficient of u_i made a constant k by choosing f' from (68) (k = (1/r)((M/2) f'^2 - sum u_j'^2) READ p.218); then
 u_i'' + k u_i = (1/(2M)) (k r + sum u_j'^2)(r q_i + (u_i'/k) sum q_j u_j') (70) (READ p.218, primes d/ds1),
whose right side vanishes without perturbation, so the unperturbed motion is again a harmonic oscillation. Constraints (READ p.218): k has the sign of a0 and an arbitrary magnitude; the method breaks down for parabolic initial data because a0 is infinite and k would be 0 (k appears in a denominator of (70)); for elliptic motion a good choice is k = M/(4 a0), whereupon u_i' = du_i/ds at the start and the Kepler-time rules (53), (54) hold in s1; the time perturbation is Delta t = integral Delta r f' ds1 = sqrt(2/M) integral Delta r sqrt(k r + sum u_j'^2) ds1 (READ p.219). This is the authors' own answer to "evolve the frequency": it is the ancestor of the constant-frequency E-time (generalised eccentric anomaly) of the book's section 19. It is elliptic-only, as the book says.

## 6. Sign and convention comparison for the KS map (what #928 plans to test against)

All sources share ONE underlying matrix. With L(u) = the 4 x 4 KS matrix [[u1,-u2,-u3,u4],[u2,u1,-u4,-u3],[u3,u4,u1,u2],[u4,-u3,u2,-u1]] (KS 1965 eq. (5); Stiefel and Scheifele (9,27)) and L3 its first three rows (3 x 4):

| Source | Object | Relation to L3 | Where |
|---|---|---|---|
| KS 1965 | x = A u, dx = 2 A3 du; du = A3^T dx/(2r); Q = 2 A3^T P; du/ds = (1/2) A3^T xdot | A = L | (8), (9), (17), (28), (37) |
| Stiefel and Scheifele 1971 | x = L(u) u; x' = 2 L u'; u' = (1/2) L^T xdot (s-derivative, x' = dx/ds); Q = (|u|^2/2) L^T (-dV/dx + P) | L = L; "(L^T P)" is the 4-vector from the transpose | (9,28), (9,43), (9,71), (9,38) |
| Peters 1968 | A* is the 4 x 3 matrix [[u1,u2,u3],[-u2,u1,u4],[-u3,-u4,u1],[u4,-u3,u2]]; P_u = 2 A* P | A* = L3^T | digest of Peters |
| Aarseth 1971 | L^T (3 x 4 transposed to 4 x 3) [[u1,u2,u3],[-u2,u1,u4],[-u3,-u4,u1],[u4,-u3,u2]] (eq. 34); R = L u; P = 2 mu L u'/R | L^T (his) = L3^T | digest of Aarseth 1971 |
| Aarseth and Zare 1974 | A_1 = 2 [[Q1,Q2,Q3],[-Q2,Q1,Q4],[-Q3,-Q4,Q1],[Q4,-Q3,Q2]] (eq. 51); q = (1/2) A^T Q (52); P_K = A p (50); p = A^T P/(4R) (53) | A_1 = 2 L3^T | digest of A and Z |

Findings:
1. NO difference in sign. KS 1965 (8), (9), (17), (28), (31), (37) agree term by term with Stiefel and Scheifele (9,27), (9,69) to (9,71), Peters's A*, Aarseth's L^T and one half of A_1 of Aarseth and Zare; the forward map (6), the inverse (16) and (9,69)/(9,70) are identical (including the x1 < 0 branch and the arbitrary rotation); the bilinear condition (7) (and (18), (23)) is u4 v1 - u3 v2 + u2 v3 - u1 v4 = 0 in the same sign as Stiefel and Scheifele (9,34); Peters's control "-u4 P_u1 + u3 P_u2 - u2 P_u3 + u1 P_u4 = 0" and Aarseth's (35) are the same equation multiplied by -1.
2. Differences are scale factors only, all documented: (a) Aarseth and Zare's A_1 and Peters's momentum map P_u = 2 A* P carry a factor 2 relative to L3^T (these are the matrices that map the physical momentum into the momentum conjugate to u in the Hamiltonian forms; the factor 2 is the same 2 as in dx = 2 L3 du); (b) the book and KS use the time derivative in s (dt = r ds, so du/ds = (1/2) L^T xdot), Aarseth 1971 and Peters use dtau with dt = R dtau; the velocity formula "xdot = 2 L u'/|u|^2" is the same with u' = du/ds; (c) KS 1965 carries the mass m and the energy through a0 (u'' + (M/(4 a0)) u = ...) whereas the book uses h = M/(2 a) so that M/(4 a0) = h/2 (COMPUTED consistent; the book's h is the NEGATIVE energy, the KS paper's a0 is the semi-major axis); Peters's h has the other sign convention (digest); (d) Aarseth 1971 uses R = r_k - r_l and Peters's printed definitions are reversed relative to his own gradient formula (digest of Aarseth 1971, footnote p.126), a relative-vector sign, not a KS-map sign.
3. KS 1965 prints no hyperbolic or Stumpff generalisation of the oscillator equation (it says hyperbolic and parabolic modifications "are not too difficult"); it prints the Stumpff principal equation for elliptic motion; the full uniform treatment is in the 1971 book (c-functions, (11,34)).
4. The perturbation transformed force: KS Q = 2 A3^T P (28) against the book's Q = (|u|^2/2) L^T P: the factor (|u|^2/2) = r/2 is absorbed because KS puts r on the right of (36), (44): KS's equations read u'' + (M/(4 a0)) u = (r/(4m)) Q' + ... = (r/(2m)) A3^T P' + ..., which equals the book's (r/2) L^T P per unit mass COMPUTED consistent.

COMPUTED test (2026-10-05, scratch, five random draws of u with rng seed 1965, normal entries): maximum absolute error over all items listed below 1.8e-15.
- x from the matrix, x = L3 u, against the polynomial form (6): 0 to rounding; |x| = u.u: e.g. u = (-0.20699972, -0.42939569, -2.60940368, 1.61001564) gives x = (-4.35836898, 8.58013103, -0.30237589) and r = 9.628367449 (both lengths equal to 1.8e-15).
- Peters's A* equals L3^T, Aarseth's L^T equals L3^T, Aarseth and Zare's A_1 equals 2 L3^T: exact (the same entries).
- L3 L3^T = r E3 and L^T L = r E4 (orthogonality, (4), (10)).
- For any dx: du = L3^T dx/(2r) satisfies 2 L3 du = dx, and the side condition (7) holds (value 0 to 1e-17) for du/ds = (1/2) L3^T xdot (the vector (37)).
- Q = 2 A* P (Peters) equals KS (28) 2 L3^T P: exact, and Q . (du/ds) = 2 P . (L3 du/ds): work invariance (27).

## 7. Which convention #928 should adopt, and the exact test

Adopt the KS 1965 / Stiefel and Scheifele 1971 convention: the 4 x 4 matrix L(u) = [[u1,-u2,-u3,u4],[u2,u1,-u4,-u3],[u3,u4,u1,u2],[u4,-u3,u2,-u1]] with x = L(u) u (fourth component identically 0), the 3 x 4 block L3 for dx = 2 L3 du, the velocity u' = (1/2) L3^T xdot in the regularising time dt = r ds (so xdot = 2 L3 u'/|u|^2), the inverse (16) / (9,69)/(9,70) with the x1 >= 0 and x1 < 0 branches, and the bilinear control u4 u1' - u3 u2' + u2 u3' - u1 u4' = 0. Reasons: it is the source's own sign, the book that holds the hyperbolic and perturbed equations uses it, and every other source reduces to it by a documented factor. Keep the book's h as the NEGATIVE energy only if the code follows the book's equations; otherwise document the sign (the project's physical energy is -h).
The test (before any dynamics test): for random u (rng draws of four normal numbers, avoid |u| < 1e-3) and a random physical velocity xdot,
1. x = L3 u; check x1 = u1^2 - u2^2 - u3^2 + u4^2, x2 = 2(u1 u2 - u3 u4), x3 = 2(u1 u3 + u2 u4), and |x| = u.u, to 1e-14 relative;
2. build the code's 3 x 4 matrix M3(u) used for the map; require M3 = L3 elementwise; require M3^T equals Peters's A* and Aarseth 1971's L^T entries as printed above, and equals (1/2) of A_1 of Aarseth and Zare eq. 51 transposed; require M3 M3^T = |u|^2 E3;
3. velocity: u' = (1/2) M3^T xdot; require the bilinear control u4 u1' - u3 u2' + u2 u3' - u1 u4' = 0 (to 1e-15 relative) and the round trip xdot_back = 2 M3 u'/|u|^2 equals xdot (1e-14);
4. inverse maps (16) for x1 >= 0 and x1 < 0 reproduce x from the forward map for a (u1, u2, u3, u4) on the same fibre (the recovered u need not equal the original u: require only that the forward map of the recovered u returns x), and that rotating u by an arbitrary angle along the fibre (14) leaves x unchanged and the control at 0 for the transformed velocity;
5. forces: Q = 2 M3^T P (KS 28) equals Peters's P_u = 2 A* P and the work identity Q . u' = 2 P . (M3 u') to 1e-14; and the right-hand side of the book's u-equation, (r/2) M3^T P, equals r Q/4 with Q = 2 M3^T P (per unit mass), which is KS (46) with W' = 0.
All sources then agree up to the documented factors listed in section 6. Tolerance: rounding (1e-14 relative); my run gave 1.8e-15 absolute on entries of size up to 10.

## 8. Printed numbers

None (no example, no table). Structural identities usable as exact tests: (10)-(13), (7), (16), (17), (28), (31), (37), (45), (46), (51) (sum alpha^2 = r0, sum beta^2 = r0 v0^2/(4 omega^2)), (52) and (53), the radius of the fibre circle sqrt(r) (14).

## 9. Reconciliation with project code

No KS map exists in `src/`; `core/cr3bp_regularized.py` is Sundman-time only (its dt = r ds is the KS time, with scale 1); `core/kepler.py` (universal variables) is the natural exact control for (53), (54) (the KS Kepler time is Stumpff's principal equation). `scripts/_validated_taylor_integrator.py` (planar Levi-Civita jet) is the planar case x = L(u) u with u3 = u4 = 0. No project constant corresponds to anything in this paper.

## 10. Techniques applicable to the project's problems

- #928: implement the map of section 7 with test items 1 to 5 first (this is the "three-source map test" the Peters and Aarseth digests recommend; this paper adds the original source and the exact entries). The KS regularised equation to code is (46) or, for the project's forces, the book's total-energy form (9,53), not (36); the work-W form (44) needs the perturbation work as an extra integral; the book's total-energy form avoids it. Controls drawn from this paper: the Kepler time (53) (against `core/kepler.py` over a full revolution), the invariants (51), (52), the side condition (23) as a monitor, the radial fall (alpha, beta collinear: sum alpha_i beta_i = sqrt(...) cos(theta) with theta = 0 gives a collision orbit, r reaches 0 smoothly), and the planar subcase (u3 = u4 = 0, Levi-Civita, same code).
- Constant-frequency time (section 5): useful only for a long elliptic arc around the moon (a captured orbit); not for a hyperbolic flyby (k fixed by a0 > 0).
- #899, #924, #929: nothing new beyond the Stiefel and Scheifele digest, except that the perturbation-of-elements forms (59) to (63) give an analytic first-order check on a weakly perturbed ellipse about the moon.

## 11. Follow-ups (not registered)

1. Add the section 7 map test to the #928 test plan verbatim.
2. Obtain M. Roessler's ZAMP account of the numerical experiments (cited as to appear) and Kustaanheimo 1964 (Ann. Univ. Turkuens. A I 73); check the corpus index first.
3. Read journal pp.214 to 217 on images (equations (54) to (69)) if the Stumpff form or first-order perturbations are implemented.
