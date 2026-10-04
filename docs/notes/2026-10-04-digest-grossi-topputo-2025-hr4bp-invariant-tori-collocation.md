# Digest: Grossi & Topputo (2025), "Computation of invariant tori in the Hill restricted four-body problem via collocation"

AAS 25-427, 2025 (Politecnico di Milano). 20 pages. No DOI is printed on the paper and none is claimed here; the PDF metadata (LaTeX, created 6 March 2025) carries no title, venue
or DOI, and the paper does not name the meeting. Treat it as "AAS 25-427" only.
Filed in the private paper corpus as `grossi-topputo-2025-invariant-tori-hill-restricted-four-body-problem-collocation-aas-25-427.pdf`.

Digested 2026-10-04. Status: the text layer was read in full; pages 4 (equations 1 to 3) and 10 (equation 12) were also read as images and agree with the text layer.
The figure pages (2 to 7 and 14 to 18) were read by caption and axis text only, not as images, so any number taken from a plot is marked "axis" (a printed tick value of a colour bar or axis) and is not a measurement.
READ = seen on the page. COMPUTED = my own arithmetic on 2026-10-04. INFERRED = my reading or comparison with project code.
Related digests (not repeated here): `2026-10-04-digest-sanaga-howell-2025-hill-restricted-four-body-ephemeris-transition.md` (Hill model derivation, periodic orbits),
`2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md` (HR4BP periodic orbits), `2026-07-27-733-olikara-thesis-haro-parameterization-book-digest.md` (the collocation method's source thesis),
`2026-10-04-digest-rosales-jorba-jorba-cusco-2021-bcp-halo-like-tori-l2.md` (BCP tori).

## 0. What the paper is, and is not

READ (abstract, p1): a collocation version of the GMOS algorithm (Gomez-Mondelo, Olikara-Scheeres) for two-dimensional invariant tori of the Hill restricted four-body problem (HR4BP), with
Fourier discretisation of the invariant curve of the stroboscopic map and Gauss-Legendre collocation in time, pseudo-arclength continuation, and stability as a by-product of the Newton Jacobian.
The method is Olikara's (2016 PhD thesis, ref 28) and Baresi's (2017, ref 27); the HR4BP application is the new part. Families computed: planar and vertical quasi-Lyapunov at the dynamical equivalents
of L1 and L2, quasi-halos at L1 (via a bifurcation of the planar family) and L2 (from a resonant pi-periodic halo), quasi-DROs (planar and vertical), and the stable and unstable manifolds of one quasi-halo.

READ: THE PAPER PRINTS NO TABLE. There is no table of initial conditions, frequencies, rotation numbers, multipliers or Jacobi constants anywhere in it. Everything numeric is either a model
constant (section 1), a discretisation setting (section 3), or a figure axis. There is no published numerical control for a torus in this paper. This matters for section 5 (it is a method
reference, not a reproduction target) and is the paper's main limit for the project.

## 1. The HR4BP as printed (pp3 to 5)

### 1.1 Assumptions and frame (p3)

READ. P0 = Sun, P1 = Earth, P2 = Moon, P3 = particle of negligible mass. Assumptions: (I) m0 much larger than m1, m2; (II) particle mass negligible; (III) the three small bodies are much
closer to each other than to P0. The Sun and the Earth-Moon barycentre (EMB) follow Keplerian orbits; the Sun's orbit about the EMB is circular; the Earth-Moon relative motion is a member of the Hill variation-orbit family,
one parameter m, "the synodic period of the Earth-Moon orbit in years" (as printed; this is the ratio of the Sun's mean motion to the Earth-Moon synodic mean motion, see section 1.4). READ: m = 0.0808.
The frame is centred on the EMB and rotates about z at the constant absolute rate Omega = 1 + 1/m, "where 1 is the normalised rotation rate of the EMB around the Sun, and 1/m is the average rotation rate of the Earth-Moon orbit"
(this is in the year-based unit; the equations below use the second unit, section 1.4). x points from the EMB to the average Moon position, z along the angular momentum of the two primaries.
Length unit: the mean Earth-Moon distance, 384,400 km; time unit: one synodic month (29.53 d) equals 2 pi time units. Earth and Moon are NOT stationary in the frame; they oscillate about mean positions.

### 1.2 Variation orbit (eq. 1, p4)

READ (image): r_v.o.(t) = ( 1 + sum_{n>=1} (a_n/a_0 + a_{-n}/a_0) cos 2t , sum_{n>=1} (a_n/a_0 - a_{-n}/a_0) sin 2t , 0 ), with a_{-n}(m), ..., a_0(m), ..., a_n(m) polynomials in m "computed ... following Reference 18",
used "up to order 6 in m, as in Reference 16" (Scheeres 1998). The coefficients themselves are NOT printed here. (The Sanaga & Howell digest records the order m^9 in their paper.)
Earth r1(t) = -mu r_v.o.(t); Moon r2(t) = (1 - mu) r_v.o.(t), with mu = m2/(m1 + m2). READ: mu = 0.01215, M = (m1 + m2)/m0 = 3.0404e-6.
Sun: r0(t)^T = M^(-1/3) a0(m) (cos t, -sin t, 0), as printed.

### 1.3 Equations of motion (eqs. 2 and 3, p4)

READ (image), state x = (r, v), r = (x, y, z):

    xdot = f(x, t) = ( v ; dV/dr (r, t) + 2(1 + m) v x e_z )

    V(r, t) = (1/2)(1 + 2m + (3/2) m^2)(x^2 + y^2) - (1/2) m^2 z^2
              + (3/4) m^2 ((x^2 - y^2) cos 2t - 2 x y sin 2t)
              + (m^2 / a0(m)^3) ( (1 - mu)/r13(t) + mu/r23(t) )

with r_i3 = |r - r_i|. a0(m) = m^(2/3) (1 - 2m/3 + O(m^2)). At m = 0 the system reduces to the CR3BP (the Sun at infinity). The vector field is pi-periodic in t: the Sun's tide has half-month
periodicity because the quadrupole is symmetric. READ: "periodic Hamiltonian system".

### 1.4 Checks (COMPUTED / INFERRED)

- COMPUTED: the Sun tide is the quadrupole of a Sun at angle theta0 - t (regressing). The printed term (3/4) m^2 ((x^2 - y^2) cos 2t - 2xy sin 2t) is (3/2) m^2 (x cos t - y sin t)^2 up to a multiple of x^2 + y^2,
  which is the quadrupole of the direction (cos t, -sin t). This is the same sense as the project's Sun after #891 and the same as `bcr4bp._sun_position` (theta = theta_sun0 - omega_sun t). It agrees with the Sanaga & Howell digest (T1 there).
- INFERRED, a unit reading that reconciles two printed statements: Omega = 1 + 1/m is in a time unit in which the Sun's mean motion about the EMB is 1 (one year = 2 pi); eq. 2's Coriolis term 2(1 + m) is in the synodic-month unit (multiply by m). The two are the same rate.
  The frame rate in the Sun-year unit, 1 + 1/m = 13.376 at m = 0.0808 (COMPUTED), times m, is 1 + m. This is a unit change, not a conflict, but the paper never states it.
- COMPUTED, a probable printing slip in the Sun position: the paper prints r0 = M^(-1/3) a0(m) (cos t, -sin t, 0). With M = 3.0404e-6 and m = 0.0808: M^(-1/3) = 69.028; a0(m) to first order = 0.17683;
  the product is 12.2 length units, but the Sun is at 388.81 Earth-Moon distances (the project's `_ANDREU_A_S`). M^(-1/3) / a0(m) = 390.4, within 0.4 percent of 388.81, which is the size of the O(m^2) term dropped from a0(m).
  So the printed formula is probably M^(-1/3) / a0(m). This is consistent with the potential carrying m^2 / a0^3 in front of the Earth and Moon terms. The Sun's position does not appear in the equations of motion (eqs. 2 and 3), so the slip does not affect the dynamics; it would only matter for drawing the Sun.
  Flag: not confirmed against Scheeres 1998 (ref 16), which is not held in the corpus under that name (grep of the index on 2026-10-04 found no entry for it).
- COMPUTED: 1 + 2m + (3/2) m^2 = 1.1714 at m = 0.0808.
- COMPUTED relation to the project's bicircular constants: `core/bcr4bp.py` takes omega_sun = 0.925195985520347 (Rosales-Jorba 2023 Table 3, time unit 1/n_sidereal-like); 1/omega_sun - 1 = 0.08085, so m = 0.0808 is the same ratio (the Sanaga digest got 0.080852).
  This ties the Hill time unit (synodic month = 2 pi) to the bicircular one (Sun period 2 pi/0.9252 = 6.7911): the ratio of the two units is 1 + m = 1.0808.

### 1.5 Comparison with the project's bicircular and quasi-bicircular models

| Item | HR4BP (this paper) | `core/bcr4bp.py` (BCP) | `core/qbcp.py` (QBCP) |
|---|---|---|---|
| Earth-Moon motion | Hill variation orbit (coherent), Fourier in t, order 6 in m, pulsating, eccentric-like | circular, rigid, Earth at -mu, Moon at 1 - mu | quasi-circular, Fourier alpha_i(theta_S) from the coupled three-body solution |
| Sun | tidal quadrupole only (Sun at infinity); direction (cos t, -sin t); pi-periodic | full point mass at 388.81 EM units, circular, direct plus indirect term; period 2 pi/omega_S | coherent, same scale |
| Forcing period | pi (in synodic-month units) | 2 pi/omega_S = 6.791 (in sidereal-like units); the BCP's fundamental is the full Sun period, no pi-periodicity | Sun synodic period |
| Sun sense in the rotating frame | regresses, angle theta0 - t | regresses, theta0 - omega_S t (#891) | regresses |
| Sun phase | enters only through t (pi-periodic) | theta_sun0 free | as BCP |
| Mass parameters | mu = 0.01215, M = 3.0404e-6, m = 0.0808 | mu = 0.0121505816, mu_S = 328900.5423, a_S = 388.8111, omega_S = 0.925196 | as BCP |
| Free parameter to switch the Sun off | m -> 0 (Sun to infinity) | mu_sun = 0 | none, structural |

INFERRED consequences. (1) The HR4BP forcing period is a half synodic month, the BCP's is a full Sun synodic period: a periodic orbit of period 2 pi k in the Hill unit is the same object as a BCP orbit closing after k Sun periods,
the commensurability condition the #884/#905 code and its check in `OUTSTANDING.md` already state. (2) The Hill model carries the Earth-Moon eccentricity-like pulsation that the BCP lacks. (3) The Hill model drops terms of relative order 1/a_S (about 0.26 percent at 388.8 EM units);
the paper does not quantify the difference and neither does this digest. (4) The project has NO HR4BP code: `grep -ri "hill\|hr4bp" src/cyclerfinder/core` finds only unrelated Mars-Hill-sphere text in `core/sunmars_wsb.py` and `core/wsb.py`. An HR4BP RHS would be new, but it is 30 lines on top of the CR3BP RHS if the variation-orbit coefficients are available
(they are not printed here; Brown 2024 and Sanaga & Howell 2025 give the model, the coefficients come from Henon-Petit 1986 or Wintner, refs 18 and 19).

## 2. Dynamical equivalents and the periodic orbits used as seeds (pp5 to 6)

READ. The CR3BP equilibria become small pi-periodic orbits ("dynamical equivalents", DE) found by continuing in m from 0 to 0.0808 with multiple shooting. At L2 a broken pitchfork bifurcation produces several pi-periodic planar Lyapunov-like orbits
(Fig. 2 labels HR4BP L1, HR4BP L2, pi-POa L2, pi-POb L2). Periodic orbits of the CR3BP survive in the HR4BP only as isolated solutions when their period is commensurate with pi (Fig. 3: DRO a/b, DPO a/b, vertical L1 a/b, Lyapunov L1 a/b, L2 halo north and south).
Solutions labelled (a) and (b) come from the same CR3BP orbit and differ in the point of the orbit used at t = 0. The L1 and L2 DE are centre x centre x saddle, the same as the CR3BP equilibria. READ: this normal behaviour implies two (planar and vertical)
one-parameter Cantorian families of 2D quasi-periodic orbits and two-parameter families of 3D ones around each pi-periodic orbit.
No initial conditions, periods or multipliers are printed for any of them (Figures 2 to 4 only).

## 3. The method as printed (pp7 to 15)

### 3.1 Torus, invariance equation, stroboscopic map (eqs. 4 to 9)

READ. Angles on the unit circle (theta in [0,1]). Torus function v(theta1, theta2): T^2 -> R^6 with the invariance equation
omega1 dv/dtheta1 + omega2 dv/dtheta2 = f(v, t) (eq. 5). The HR4BP is forced with period T = pi, so tau = t/T, omega1 = 1/T1 = 1, and theta1 is time. Fixing theta1 and integrating one forcing period
(Delta tau = 1) shifts theta2 by the rotation number rho = T1 omega2 = omega2/omega1 (irrational for non-resonance). The invariant-curve equation of the stroboscopic map phi_T1 is w(theta2 + rho) = phi_T1(w(theta2)) (eq. 6).

### 3.2 Cylinder function and the unfolding parameter (eqs. 7 to 12)

READ. The cylinder function u(tau; theta2) on [0,1] x T^1 relates to the torus by v(theta1, theta2) = u(tau; theta2 - rho tau), tau = theta1 (eq. 8); closure is u(0; theta2) = u(1; theta2 - rho) (eq. 9).
The dynamics is the augmented field (eqs. 10, 11):

    (1/T) du/dtau = f_a(u, tau, theta2; lambda) = f(u, tau) + lambda J du/dtheta2

with unfolding parameter lambda (a free scalar that makes the system square, Olikara 2016) and J the constant matrix of the gradient of a local integral (the torus action) (eq. 12, image):
with m~ = 2(m + 1). Read from the page image: the upper-left 3x3 block has +m~ at (1,2) and -m~ at (2,1); the upper-right 3x3 block is -I (negative identity); the lower-left 3x3 block is +I; the lower-right block is zero.
So J = [[A, -I],[I, 0]] with A = [[0, m~, 0], [-m~, 0, 0], [0, 0, 0]], m~ = 2(1 + m), i.e. the Coriolis coefficient of eq. 2 (the Coriolis rotation 2(1+m) v x e_z appears as the antisymmetric block).
INFERRED: this is the Jacobian of the symplectic-structure gradient of the Hamiltonian-like quantity in the rotating frame (the momentum shift p = v + (1+m) e_z x r) and matches Olikara's choice for the CR3BP with 2 replaced by 2(1+m).
Not checked against Olikara's thesis in this digest.

### 3.3 Discretisation (eqs. 13 to 23)

READ. Fourier in theta2 truncated at M_F, N_F = 2 M_F + 1 samples, c_k(tau) complex coefficients (eq. 13); the collected state X(tau) = {x^(k)(tau) = u(tau; k/N_F)}, k = 0..N_F-1 (eq. 14), and field
F(X, tau; lambda) = {f_a(x^(k), tau, k/N_F + rho tau; lambda)} (eq. 15); X' = T F (eq. 16). Boundary condition X(0) = R(-rho) X(1) (eq. 17), R built from DFT, phase shift -rho and inverse DFT.
Time: N intervals with N+1 evenly spaced nodes, M Gauss-Legendre collocation points per interval (roots of the degree-M Legendre polynomial rescaled to (0,1)); each trajectory is a degree-M Lagrange polynomial
X(tau) = sum_{j=0..M} X_{i,j} l_j(tau^), tau^ = (tau - tau_i0)/(tau_{i+1,0} - tau_i0) (eqs. 18, 19). Continuity: sum_j X_{i,j} l_j(1) = X_{i+1,0} (eq. 20). Collocation, N M n* constraints (eq. 21):
sum_m X_{i,m} dl_m/dtau^ (tau^_j) = (T/N) F(X(tau_{i,j}), tau_{i,j}; lambda), i = 0..N-1, j = 1..M (the factor T/N comes from d/dtau = N d/dtau^). Quasi-periodic closure X_{0,0} = R(-rho) X_{N,0} (eq. 22).
Unknown vector X in R^(n* N (M+1) + n*) with n* = n N_F (eq. 23); n = 6.

### 3.4 Continuation (eqs. 24 to 30)

READ. One extra scalar phase condition anchors the theta2 phase: (X_{0,0} - X~_{0,0})^T dX~_{0,0}/dtheta2 = 0 (eq. 24), the derivative by DFT (the "minimise phase difference to the previous curve" condition of McCarthy-Howell, Olikara-Scheeres).
Unknowns Y = (X, rho, lambda), dimension N_Y = n* N (M+1) + n* + 2 (eq. 25), the "+2" being rho and lambda. Pseudo-arclength constraint (eq. 26): (1/(N1 N_F)) (X - X~)^T X~' + (rho - rho~) rho~' - Delta s = 0, N1 = N(M+1)+1.
The constraint function H: R^NY -> R^NY collects continuity, collocation, quasi-periodicity, phase and arclength. The tangent is the null vector of the Jacobian with the lambda column and the collocation row removed, normalised so that (1/(N1 N_F)) X'^T X' + rho'^2 = 1 (eq. 27).
Bifurcations are flagged when the null space grows; direction is kept with the flipping condition of Dahlke and Bettinger (ref 37). Predictor Y(0) = (X~ + Ds X'~, rho~ + Ds rho~', lambda~) (eq. 29); corrector Newton, Y(n) = Y(n-1) - DH^-1 H, tolerance ||H|| < 1e-12 (eq. 30, printed with a plus sign: see the note below).
The unknowns count equals the constraint count because lambda is free, so DH is square and plain Newton with sparse LU applies. READ: because lambda is exponentially small, O(lambda) terms are dropped from the collocation Jacobian df_a/dX_{i,j} to keep it sparse "without affecting convergence".

NOTE on eq. 30: the text layer prints Y(n) = Y(n-1) + DH^-1 H. A Newton step is Y - DH^-1 H. INFERRED typesetting slip; not read as an image.

### 3.5 Initial guess (eq. 28)

READ. From the monodromy matrix of the underlying periodic orbit, the eigenvector v_lambda of a unit-modulus complex eigenvalue gives the initial curve
u(0)(0; theta2) = Re[v_lambda] cos(2 pi theta2) - Im[v_lambda] sin(2 pi theta2), and the initial rotation number is "rho(0) = tan^-1(lambda_I / lambda_R)".
INFERRED: the printed rho is in turns (theta in [0,1]) and the monodromy matrix is over one forcing period T = pi, so the correct estimate is atan2(lambda_I, lambda_R)/(2 pi), not the bare arctangent in radians. The printed form is a typesetting shortcut or an error;
the project's own `genome/qp_tori._seed_invariant_circle` and `_canonicalize_ns_eigenpair` already handle this convention.

### 3.6 Discretisation settings and error control (pp14 to 15)

READ. "N = 20 and M = 7" typically per family; initial N_F = 10; N_F is raised by 5 whenever the invariance error on a finer grid (each point rotated by -rho and propagated for T = pi under the true flow) exceeds the tolerance;
the step size is reduced if N_F rises more than 10 times between two solutions or Newton fails. Fig. 8 example sparsity: N = 20, M = 8, M_F = 10 gives nnz = 324744 (axis of Fig. 8: matrix size about 2 x 10^4 per side, the printed axis shows scale 10^4 and 2). These are the only discretisation numbers printed.
COMPUTED: the text says the initial N_F is 10 and N_F grows by 5 (10, 15, 20, ...), but N_F = 2 M_F + 1 is odd by definition (eq. 13). INFERRED: either the implementation does not require N_F odd, or the text means M_F = 10 (N_F = 21).
Fig. 8 says "Fourier order M_F = 10". The unknown count for N = 20, M = 8, M_F = 10 (N_F = 21, n* = 126) is n* N (M+1) + n* + 2 = 126*20*9 + 126 + 2 = 22,808 (COMPUTED). With N_F = 10 (n* = 60) it would be 10,862. The sparsity plot's axis labels run to 2 on a 10^4 scale, which fits 22.8 thousand and not 10.9 thousand (axis text only, plot not viewed), so M_F = 10 with N_F = 21 is the likelier reading.

### 3.7 Stability (p15)

READ. Linear stability follows Olikara 2016: assumes the tori are reducible (Jorba 2001, ref 38); the Floquet matrix eigenvalues come "directly from the Jacobian DH"; unstable eigenvectors are first computed on the initial invariant curve and propagated over the whole torus by the flow.
Perturbing along them and integrating forward and backward builds the hyperbolic manifolds (Fig. 18). The paper prints no eigenvalue or Floquet exponent.

## 4. Results as printed (pp15 to 18; all from figures)

READ (captions and axis text):
- Fig. 9 planar quasi-Lyapunov at L1: invariant curves at tau = theta1 = 0 span about x = 0.75 to 0.95, y = -0.15 to 0.15 (axis); colour bar T2 = 2.5 to 2.9 (axis). L2 planar: x = 1.12 to 1.16, y = -0.06 to 0.06; T2 about 3.12 to 3.14 (axis), approaching T1 = pi = 3.14159, a 1:1 resonance.
- Text (p15): the L1 planar family is bounded by the pi-periodic Lyapunov orbits of Fig. 2; the L2 planar family is constrained by a 1:1 resonance with T2 tending to pi, last members in a nearly-resonant region from pi-POa and pi-POb (Olikara-Gomez-Masdemont 2016, ref 35).
- Fig. 11 vertical family at L1: T2 about 1.42 to 1.52 (axis).
- A bifurcation on the planar L1 family produces northern and southern quasi-halo families (ref 23). Fig. 12 L1 quasi-halo curves: T2 about 2.54 to 2.72 (axis). Fig. 14 L2 quasi-halo: T2 about 2.6 to 3.0 (axis); these start from the resonant pi-periodic northern halo of Fig. 3.
- Figs. 16, 17: quasi-DROs from a planar centre and a vertical centre of a pi-periodic DRO; Fig. 18: stable and unstable manifolds from two invariant curves of a quasi-halo, and the image of the whole torus under the hyperbolic flow (a 2-torus).
- Conclusion (p18): future work on boundaries of tori (dynamical limits) and higher-dimensional tori.
- Claim (p2, p15): the collocation version "demonstrated superior efficiency" to shooting, but NO run time, iteration count or comparison is printed in this paper (the comparison is referred to Baresi 2017, Chapter 3.5).

The values marked axis are plot ticks, good to a few units in the last digit, and are not suitable as test expectations. T2 = pi/rho (p15: T2 = T1/rho = pi/rho) so rho can be inferred from them (COMPUTED example: T2 = 2.7 gives rho = pi/2.7 = 1.16 turns per forcing period, i.e. the rotation number defined mod 1 is 0.16, which shows
the paper's T2 and rho are not mod-1 equivalent; see section 7).

## 5. Printed numbers usable as sourced tests

No printed number is a torus or periodic-orbit result. The usable ones are model constants and discretisation settings. Page and equation given.

| ID | Quantity | Printed value | Where | Use |
|---|---|---|---|---|
| G1 | Hill parameter m for Sun-Earth-Moon | 0.0808 | p3 | check 1/omega_sun - 1 = 0.08085 against `bcr4bp._ANDREU_OMEGA_S` (3 s.f.; the printed m is rounded, so tolerance 1e-3 relative) |
| G2 | mu = m2/(m1+m2) | 0.01215 | p4 | against `_ANDREU_MU_EM` = 0.0121505816 (rounded to 4 s.f.) |
| G3 | M = (m1+m2)/m0 | 3.0404e-6 | p4 | reciprocal 328899 against `_ANDREU_MU_S` = 328900.54 (Sun mass in Earth+Moon units; 1/3.0404e-6 = 328,904: COMPUTED, agrees to 1.1e-5 relative, i.e. to the printed digits) |
| G4 | Length unit | 384,400 km | p3 | |
| G5 | Synodic month | 29.53 d = 2 pi time units | p3 | time unit 4.700 d (COMPUTED) |
| G6 | Forcing period T | pi | pp4, 8 | HR4BP pi-periodicity test; the BCP counterpart is 2 pi/omega_S with the Sun sense fixed |
| G7 | Maximum inclination neglected | about 5 degrees | p4 | |
| G8 | Eq. 3 coefficients | 1 + 2m + 3/2 m^2, -1/2 m^2 z^2, 3/4 m^2 [(x^2-y^2) cos 2t - 2xy sin 2t], m^2/a0^3 | p4 | structural identity test (when an HR4BP RHS is built): at m -> 0 reduces to the CR3BP; the tide equals the quadrupole of a Sun at angle theta0 - t (shared with the Sanaga digest T1) |
| G9 | a0(m) | m^(2/3) (1 - 2m/3 + O(m^2)) | p4 | |
| G10 | Coriolis coefficient | 2(1+m); J entry m~ = 2(m+1) | pp4, 10 | 2.1616 at m = 0.0808 (COMPUTED) |
| G11 | Discretisation | N = 20, M = 7 typical; initial N_F = 10; N_F steps of 5; Newton tolerance 1e-12 | pp14 to 15 | for a corrector test configuration, not a result |
| G12 | Sparsity | N = 20, M = 8, M_F = 10, nnz = 324744 | p14 Fig. 8 | structural: a sparse-Jacobian builder with these settings should have 324,744 nonzeros if built with the same block structure (an implementation-specific number; useful only if the same stencil is used; see 3.6 for an unknown-count check) |
| G13 | L1, L2 planar families bounded by | T2 -> pi (L2) | p15 | qualitative: 1:1 resonance wall at the pi-periodic orbit |
| G14 | Normal behaviour of the L1 and L2 DE | centre x centre x saddle | p6 | a Floquet-spectrum check once an HR4BP periodic-orbit finder exists |

Not printed, so no test possible from this paper: any periodic-orbit initial condition, period, multiplier, rotation number, Jacobi constant or torus amplitude.

## 6. Reconciliation with project code

Existing (grep and reads of 2026-10-04):

- `core/bcr4bp.py`: BCR4BP RHS, STM, `propagate_bcr4bp`, `sun_commensurate_period`; Sun at theta_sun0 - omega_S t (#891). It is the Sun-point-mass model; no Hill model.
- `core/qbcp.py`: QBCP RHS in the (p, m) canonical frame with the eight alpha_i(theta_S), STM, `propagate_qbcp_pv`.
- `genome/qp_tori.py`: the Olikara-Scheeres GMOS corrector, SHOOTING variant: `correct_qp_torus`, `_gmos_residual` (stroboscopic flow integration), `_seed_invariant_circle`, `_canonicalize_ns_eigenpair`. No collocation, no unfolding parameter lambda, no pseudo-arclength in this file.
- `genome/bcr4bp_torus.py` (276 lines): BCR4BP torus corrector built on `qp_tori.evaluate_invariant_circle` (shooting-GMOS lineage); `genome/qbcp_torus.py`: the same for the QBCP (works for Sun-Earth L2, hopeless for Earth-Moon L1/L2 because of the one-period amplification).
- `search/variational_qp_torus.py` (CR3BP, 767 lines), `search/variational_qbcp_torus.py` (QBCP, 1277 lines, #617), `search/variational_ccr4bp_torus.py` (#690), `search/variational_crnbp_torus.py` (#720): a SEEDLESS 2D pseudospectral (Fourier in both angles, algebraic residual on a grid, no integration in the residual) torus corrector with free rotation numbers. This is the harmonic-balance relative of the collocation method: two Fourier angles instead of Fourier x piecewise Gauss-Legendre polynomials.
- `search/pseudo_arclength.py`, `search/qp_torus_fixed_jacobi_continuation.py`, `genome/qp_tori_energy_walk.py`: continuation machinery; `search/er3bp_floquet.py` and `search/variational_periodic_orbit*.py`: Floquet and periodic-orbit code.
- `data/validation/v1_qp.py`, `v2_qp.py`, `v3_qp.py`: the QP-torus validation tiers (v2 is the long-span invariance residual test).

Absent: an HR4BP RHS; Gauss-Legendre piecewise collocation of the cylinder function; the unfolding parameter lambda and the matrix J; a Newton Jacobian whose by-product is the torus Floquet spectrum; a sparse LU solve of the collocation system.
Stability of tori in the project is by other means (the GMOS corrector's invariant-curve Floquet analysis in `qp_tori`; hyperbolic-whisker code in `search/ccr4bp_whisker.py`). The 2026-07-27 Olikara-thesis digest already records that `qp_tori.py` has no collocation counterpart (its lines 140 to 153).

Conflicts: none with the project's Sun sense (section 1.4). Two points need a decision rather than being conflicts: (1) the paper's rho convention is mod 1 in turns but the T2 axis of its colour bars is pi/rho, not mod 1 (section 7);
(2) the project's pseudospectral torus correctors use a fixed Fourier grid in both angles and the BCP/QBCP forcing period, while this paper discretises time by collocation on the HR4BP's pi-period. Both treat the forcing angle as a torus angle; nothing in the paper contradicts them.

What the paper's method would add: (a) a COLLOCATION torus corrector whose Jacobian is sparse and square (unfolding parameter), which converges from poor guesses with a wide basin ("significant numerical stability", p2, via ref 27); like the existing seedless pseudospectral correctors it avoids one-period amplification of unstable flows; in addition it yields a continuous torus, arclength continuation and bifurcation detection (null-space growth) in one framework; (b) the floquet spectrum of a torus as a by-product; (c) adaptive Fourier order along a family. What it does NOT add beyond the existing `variational_*_torus` modules: the ability to cross the shooting wall (already crossed by #612/#617).

## 7. Errata and ambiguities flagged (benefit of the doubt: typesetting)

1. Sun position r0 = M^(-1/3) a0(m)(cos t, -sin t, 0) gives 12.2 EM units; M^(-1/3)/a0(m) gives 390.4 against the true 388.8 (section 1.4). Probably a printed slip; does not enter the EOM.
2. Eq. 30 prints Y(n) = Y(n-1) + DH^-1 H; Newton is a minus.
3. Eq. 28 text: rho(0) = atan(lambda_I/lambda_R) lacks the 2 pi (and should be atan2) if rho is in turns.
4. T2 = T1/rho = pi/rho in the results while rho is an angle fraction defined mod 1: for the L1 vertical family T2 about 1.42 to 1.52 means rho = pi/T2 = 2.07 to 2.21 turns, which is not mod-1 equivalent to the rotation number of the monodromy eigenvalue. Either rho is carried unreduced along the family (continuous in the parameter), or T2 is reported modulo a branch. Not settled by the text.
5. "initial N_F = 10" against N_F = 2 M_F + 1 (section 3.6).
6. The statement that m "is the synodic period of the Earth-Moon orbit in years" is imprecise; m = 0.0808 is dimensionless and equals (1 year / synodic month)^-1 = 29.53/365.24 = 0.0808 (COMPUTED: 29.5306/365.2422 = 0.08085). So m is the synodic month expressed in years, which is correct as printed. No slip.

## 8. Techniques applicable to the project's problems

Scope rule from the coordinator: #890 and #895 are the Titania-Oberon (Uranian) orbit and its real-ephemeris arcs, not Earth-Moon solutions, so only the METHOD is mapped to them. Earth-Moon Sun-forced work is #884, #905 and #902.

### 8.1 #902 (rebuild the bicircular validation tiers on a real orbit)

The paper has no orbit to use as a control (no table). Its contribution is a different one: the HR4BP is the Sun-at-infinity limit of the BCP, so the L1 dynamical-equivalent periodic orbit of the Hill model (Fig. 2) and the BCP L1 replacement orbit reproduced from Jorba et al. 2020 should differ by O(1/a_S) (about 0.3 percent in position, INFERRED). That gives a cross-model sanity band for the rebuilt L1 orbit, not a published control: a test asserting the BCP orbit sits within a stated fraction of the Hill-model orbit would be self-computed on both sides and must not be labelled a golden value. Use it only after Brown 2024 or Henry et al. prints an L1 DE initial condition. The L1 DE of the paper is centre x centre x saddle; the same classification test on the BCP L1 replacement orbit is a published-type structural control (G14) and costs nothing.

### 8.2 #884 / #905 (Sun-forced Earth-Moon search rerun in the corrected models)

1. Resonance structure: the paper states that periodic orbits persist only if their period is commensurate with pi (Hill) and that, off resonance, a CR3BP periodic orbit becomes a quasi-periodic 2-torus with the forcing frequency as one basic frequency (p5). Under #905 the closure-over-p-Sun-periods check is the same statement; the torus view says the orbits found by the search are the isolated resonant members of a family of tori, and the nearby non-resonant members are free (no closure needed) objects. A cycler needs a closed orbit, but a quasi-periodic neighbour at the same Jacobi constant gives an independent test of "survives the Sun": the invariant-curve closure residual of the stroboscopic map.
2. Bifurcation detection in the family walk (the #905 review item "bifurcation detection in the family walk"): the paper's rule, null space of the reduced Jacobian grows (p13), is for tori; the periodic-orbit analogue is a multiplier crossing 1 or the null space of (DPhi - I) growing, which the family walk can monitor with the existing STM. The paper's L2 broken pitchfork (Fig. 2) is a concrete warning that a walk from L2 in m yields several branches.
3. Flip condition (Dahlke and Bettinger, ref 37) for keeping orientation along a family: directly applicable to the family walk, given the tangent sign.
4. The Hill model as a fast intermediate model: pi-periodic forcing halves the closure time and removes the Sun's phase degree of freedom; useful as a cheap pre-screen only if HR4BP code is added, and the screen would need the BCP/QBCP re-check (the Sanaga digest argues the same for the L2 halo transition region).
5. Do not use Fig. 3's orbits as seeds: no numbers are printed.

### 8.3 #890 / #895 (Titania-Oberon orbit and its real-ephemeris arcs; METHOD ONLY)

The Uranian problem is a planar two-moon model with several commensurate frequencies. The method maps in two ways. (1) The #890 orbit as a one-frequency forced periodic orbit and the extra frequency (Oberon's eccentricity) as a second torus angle: that is exactly #922's torus question, and this paper's cylinder-function / stroboscopic-map formulation (eqs. 6 to 11) is the formulation to use; the unfolding-parameter trick makes the Newton system square. (2) The paper's adaptive N_F rule (rotate by -rho, propagate for one period, compare) is a ready-made invariance gauge for the `#895` real-ephemeris arcs: take points of the model torus, propagate them in the real-ephemeris model for one forcing period, and measure the distance to the rotated curve. This replaces a statement of the form "the arcs are near the orbit" by the invariance residual of the paper's own error-control step. The physics (Hill approximation, a Sun at infinity) does not transfer: the Uranian moons are not in a Hill limit and the paper's quadrupole Sun term has no analogue.

### 8.4 #913 to #918 and #921 to #924

- #913 (77 Earth-Mars real-ephemeris cycler solutions), #914, #915, #918: the paper has nothing for heliocentric cyclers; no mapping. A statement of absence rather than a stretch.
- #916 (persistence conjecture of Ross and Roberts-Tsoukkas, stable near-commensurate cyclers as invariant curves of the one-period map): this is the problem the paper's method solves. The invariant curve of the one-period map, its rotation number rho, and its stability from the Jacobian are the objects. A cycler's one-period map in the corrected BCP/QBCP is the stroboscopic map; the stable cyclers' near-commensurate orbit is the periodic point whose neighbours form the curve. The corrector needs: the stroboscopic map by integration (existing STM propagators), the Fourier curve (existing `genome/qp_tori.py`), and a continuation in the perturbation strength (Sun mass, eccentricity), with the rotation number as an observable. Collocation is worth adding when the cyclers are unstable (large one-period amplification), which is the same reason the project built the seedless pseudospectral correctors.
- #922: as in 8.3.
- #917 (Casoliva rows in the elliptic problem): the ER3BP is a periodic perturbation with the same structure (forcing angle = true anomaly); the paper's method applies verbatim with f = the ER3BP field and T = 2 pi.
- #921, #923: no mapping (navigation budget and taxi delta-v). The torus Floquet spectrum could give a stable-direction growth rate as one input to a maintenance budget, but the paper prints no such quantity.
- #924 (FLI, MEGNO screen): a torus family gives the regular-motion reference; the paper's invariance-error control is itself a check against chaos near resonant tongues (the L2 planar family ends at a 1:1 resonance). It is a cross-check for the screen's regular case, not a substitute.

## 9. Recommended follow-ups (not registered; no task numbers assigned)

1. Add an HR4BP right-hand side and STM (eqs. 2 and 3) with the m^2 expansion of the variation-orbit coefficients from a source that prints them (Henon-Petit 1986 or Wintner; check the corpus first), tested by: reduction to `cr3bp` at m = 0; the Sun tide equals the quadrupole of the project's regressing Sun; and a published periodic-orbit control (Brown 2024 or Henry 2023 print initial conditions if any; check the Brown digest).
2. Decide, with an experiment, whether a Gauss-Legendre collocation corrector adds anything over `variational_qbcp_torus` and `genome/qp_tori`: take the #617 torus that the pseudospectral corrector finds and compare Newton iterations, basin and wall time with a collocation implementation. The paper claims faster than shooting and prints no numbers; the pseudospectral correctors already beat shooting in the project, so the expected gain is small (a judgement, not measured).
3. Resolve eq. 28 and the T2 versus rho convention (section 7, items 3 and 4) against Olikara 2016 (held, digested 2026-07-27) before reusing the rho handling.
4. For #916: build the invariant-curve test of the one-period map for the printed stable cycler members, reusing `genome/qp_tori.py`, and record rho and the Floquet spectrum from the Jacobian.
5. For #905: add the null-space and multiplier-crossing bifurcation flag and the Dahlke-Bettinger orientation flip to the family walk.
6. For #922: set the problem up in the cylinder-function form (section 3.2); use the paper's adaptive invariance check (3.6) as the comparison metric for the #895 arcs.
7. Check the Sun-position slip (section 1.4) against Scheeres 1998 if that paper can be obtained, and file the result.
8. Obtain the three conference references this paper leans on that the corpus lacks: Henry et al. 2023 (AIAA, ref 23) and Henry et al. ISTS 2023 (ref 14), which carry the HR4BP torus results with tables of rotation numbers, if they exist; check `CORPUS_INDEX.md` first.
