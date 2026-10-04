# Digest: Jorba and Olle 2004, invariant curves near Hamiltonian-Hopf bifurcations of four-dimensional symplectic maps

Date: 2026-10-05 (Sydney). Reading, reasoning and independent computations; no project code or tests were changed.

Source: Angel Jorba and Merce Olle, "Invariant curves near Hamiltonian-Hopf bifurcations of four-dimensional symplectic maps",
Nonlinearity 17:691-710 (2004), DOI 10.1088/0951-7715/17/2/019, received 29 January 2003, in final form 7 October 2003, published
27 January 2004 (Universitat de Barcelona and Universitat Politecnica de Catalunya). Filed in the private paper corpus as
`jorba-olle-2004-invariant-curves-hamiltonian-hopf-bifurcations-4d-symplectic-maps-nonlinearity-17-691-doi-10.1088-0951-7715-17-2-019.pdf`.
A 21-page PDF with a publisher cover page and a good text layer (the model maps and tables are cleanly printed). I read the whole text.

Evidence tags: READ (p.N) is the printed page. COMPUTED is my own computation (scratch scripts, not committed; recipes in section 4). INFERRED is my reasoning.

Companions: `docs/notes/2026-10-05-digest-olle-pacha-villanueva-2004-l4-vertical-hopf.md` (the restricted-problem application, which cites this paper for the numerical method),
`docs/notes/2026-10-05-digest-olle-pacha-villanueva-2005-normal-form-1-1-resonance.md` (the analysis).

## 0. Answer first

A numerical paper on two exactly specified model maps with closed forms, which makes it the best sourced test source of the three Hopf papers. READ: the abstract states a "numerical description of an extended neighbourhood of a fixed point of a
symplectic map undergoing an irrational transition from linear stability to complex instability ... both the direct and inverse cases". All the quantities that need a formula (Jacobian, stability parameters, critical parameter, critical eigenvalue) are
printed in closed form, and Table 1 (the rotation numbers at the two ends of the Lyapunov families) and Fig. 2 (two invariant curves with 11- and 12-digit data) are printed to the digit. I reproduced all of them (section 4): Table 1 to 14 digits, omega_crit = arccos(3/4) exactly, and
both Fig. 2 curves (the rotation numbers 0.6008342478 and 0.798309349426 from the printed x~2 values) with an independent Fourier-Newton invariant-curve solver.

## 1. The model maps (READ pp.693-695, section 2)

For x = (x1, x2, x3, x4) in R^4 and parameters K1, K2, L1, L2 (eq. in section 2, p.693):

    Ts(x) = ( x1 + K1 sin(x1 + x2) + L1 sin(x1 + x2 + x3 + x4),  x1 + x2,  x3 + K2 sin(x3 + x4) + L2 sin(x1 + x2 + x3 + x4),  x3 + x4 )
    Tt(x) = the same with sin(x1 + x2 + x3 + x4) replaced by tan(x1 + x2 + x3 + x4)

generalisations of the standard map (used by Froeschle 1972 with L1 = L2, and by Pfenniger 1985a and Olle and Pfenniger with L1 + L2 = 0). Symmetry: T(-x) = -T(x). The origin is a fixed point for all parameters; both maps have the same Jacobian

    A = [[1 + K1 + L1, K1 + L1, L1, L1], [1, 1, 0, 0], [L2, L2, 1 + K2 + L2, K2 + L2], [0, 0, 1, 1]].

For L1 + L2 = 0 the maps are symplectic with conjugate pairs (x1, x2) and (x4, x3) (COMPUTED check: A^T W A = W with the corresponding W, residual 0 at K = -1, L = 0.24; also det A = 1 and A against a finite-difference Jacobian of Ts to 1.3e-15).
Characteristic polynomial z^4 + alpha z^3 + beta z^2 + alpha z + 1 = (z^2 + b1 z + 1)(z^2 + b2 z + 1), with b1 = -(lambda + 1/lambda), b2 = -(mu + 1/mu), Delta = (b1 - b2)^2, alpha = b1 + b2, beta = 2 + b1 b2 (the same parameters as the Hadjidemetriou 1975b digest,
section 2.3, where Delta was defined as alpha^2 - 4(beta - 2); the two are equal: (b1 - b2)^2 = (b1 + b2)^2 - 4 b1 b2). READ (p.694): the origin is linearly stable if |b1| <= 2, |b2| <= 2 and Delta >= 0, complex unstable when Delta < 0. Closed forms:
alpha = -(4 + K1 + K2), beta = 6 + 2(K1 + K2) + K1 L2 + K2 L1 + K1 K2, Delta = (K1 - K2 + L1 - L2)^2 + 4 L1 L2.
Restricting to K2 = 0, L1 = -L2 = L, K1 = K: Delta = K (K + 4L); for -8 < K < 0 the stability-to-complex-instability transition (Delta = 0) in L is at L_crit = -K/4. At L_crit the Jacobian A_crit has the double eigenvalue
lambda = 1 + K/4 + i (1/4) sqrt(8|K| - K^2) (a conjugate pair, modulus 1), with eigenvector u + i v = ((K + i sqrt(8|K| - K^2))/4, 1, (K + i sqrt(8|K| - K^2))/4, 1)^T and invariant plane x1 = x3, x2 = x4 (p.694).
The paper fixes K = -1 (similar results for other -8 < K < 0): L_crit = 1/4, Ts being the direct case and Tt the inverse case. The closed form inverse of Ts is printed (p.694); the closed form inverse of Tt is printed with auxiliary functions f, g (eqs 1, 2, p.695) and D = cos(x2 + x4); Tt is a local diffeomorphism near the origin, not one-to-one on R^4.

## 2. The paper's results (READ pp.695-706)

### 2.1 Direct case, Ts with K = -1 (section 3)
- For L < L_crit: two Lyapunov families of invariant curves (2D tori of the flow) born at the origin with rotation numbers tending to the linear angles omega_i of the eigenvalues exp(+-i omega_i); the families are Cantorian (holes below double-precision resolution, exponentially small for Diophantine frequencies); both exist for L = 0.24, 0.245, 0.249 (Fig. 1, p.695) and meet at omega_crit when L -> L_crit. Normal behaviour of the curves near the origin is that of the elliptic point.
- For L > L_crit: the origin is complex unstable; its 2D unstable and stable manifolds are almost coincident (a small "loop") and trajectories are trapped nearby; for L = 0.26 the iterates of the initial curve on the unstable manifold (k = 1, 200, 270, 300; Fig. 3, p.698) and of the stable manifold (Fig. 4) go far away and come back chaotically (splitting of separatrices; the manifolds coincide for an integrable system); slices by x1 = 0 (Fig. 5); for L = 0.28 and 0.30 the manifolds grow and are no longer almost coincident (Fig. 6), and as L decreases to L_crit they shrink to the origin.
- Hopf unfolding (section 3.3, Fig. 7, p.701): the two Lyapunov families that exist for L < L_crit meet at L_crit and detach from the origin as one family for L > L_crit (continued to L = 0.26, 0.27, 0.29, 0.31), remaining elliptic (normal eigenvalues 1, 1, exp(+-i nu); Fig. 8 for L = 0.24, 0.25, 0.27): so tori coexist with the invariant manifolds of the unstable origin.
- Confinement (section 3.4): 2D tori through the initial condition x_i = 0.001 (i = 1..4) before (L = 0.249), at (0.25) and after (0.251) the transition (Fig. 9): the plot size grows after the transition because of the manifolds; at L = 0.25 the (x1, x3) projection of the tori collapses to a plane (the invariant plane x1 = x3). Fig. 10: slices x1 = 0 of tori from x1 = x2 = x3 = x4 = 0.1, 0.07, 0.001 at L = 0.249, 0.25, 0.252, 0.254: topology change near the origin as L crosses L_crit; invariant curves and tori coexist with the manifolds and the manifolds do not separate regions of the 4D space.
### 2.2 Inverse case, Tt with K = -1 (section 4)
- For L < L_crit the origin is surrounded by 2D invariant tori. The two Lyapunov families born at the origin meet at some distance from it and form one global family that begins and ends at the origin; near the origin the curves are elliptic, farther they are hyperbolic, with two parabolic transition curves per family (Fig. 11, p.704; L = 0.24, 0.245, 0.249; Fig. 12: rotation number against the modulus of the normal eigenvalues). As L -> L_crit this connecting loop collapses to the origin, the hyperbolic curves and their nearly connected manifolds merge with the origin and it becomes parabolic at L = L_crit. For L > L_crit no invariant curves remain, the origin is complex unstable and its manifolds are far from connected (Fig. 13 for L = 0.24, 0.249, 0.2501, 0.251): "all the nearby trajectories escape" (the dynamics after the bifurcation is trivial). Compare Meyer 1998 and Pacha 2002.
### 2.3 Methods (section 5, pp.706-709)
- Invariant curves: x(theta) = a0 + sum_{k=1..N} (a_k cos k theta + b_k sin k theta) with rotation number omega assumed known and irrational, F(x)(theta) = f(x(theta)) - x(theta + omega) discretised on 2N + 1 equispaced points; Newton with the analytic differential (chain rule); one coordinate fixed at theta = 0 to remove the phase freedom (an extra row; Gaussian elimination with pivoting leaves it as 0 = 0); N chosen automatically so that E(x, omega) = sup |f(x(theta)) - x(theta + omega)| < 1e-12, tabulated on a mesh finer than the discretisation (p.707).
- Normal behaviour: the linear skew product x' = A(theta) x, theta' = theta + omega, A(theta) = D f(x(theta)), analysed through the operator (L_omega psi)(theta) = A(theta - omega) psi(theta - omega); Proposition 5.1: if lambda is an eigenvalue so is lambda exp(i k omega); Definition 5.1 (omega-unrelated eigenvalues); Proposition 5.2: n omega-unrelated eigenvalues give the reduction to B = diag(lambda_1..lambda_n); the spectrum of the truncation (dimension n (2N + 1)) is computed by shifted QR; all curves in the paper are reducible to 1e-12 (p.708).
- Manifolds: a closed curve sigma(s) = h (cos s u + sin s v)/norm on the linear approximation spanned by u + i v (eigenvector of r2 exp(i omega), 0 < r1 < 1 < r2), with h <= 1e-5 giving an O(h^2) error below 1e-15 (p.709); iterates remeshed by the Fourier transform to equal arc parameter; curves for h between h0 and r2 h0 cover a fundamental interval; slices by a coordinate hyperplane from the Fourier coefficients; stable manifold by the inverse map (closed form here); repeated in quadruple precision with no visible difference.

## 3. Printed numbers usable as tests

| item | value | page | status |
|---|---|---|---|
| maps and Jacobian | Ts, Tt, A (section 1) | 693-695 | exact; Jacobian reproduced (1.3e-15) |
| alpha, beta, Delta (general) | alpha = -(4 + K1 + K2); beta = 6 + 2(K1 + K2) + K1 L2 + K2 L1 + K1 K2; Delta = (K1 - K2 + L1 - L2)^2 + 4 L1 L2 | 694 | verified (K = -1: alpha = -3, beta = 4 + L, Delta = 1 - 4L) |
| restricted: Delta, L_crit | Delta = K(K + 4L); L_crit = -K/4 (K = -1: 1/4) | 694 | verified |
| critical eigenvalue | lambda = 1 + K/4 + (i/4) sqrt(8 abs(K) - K^2) | 694 | verified: eigenvector residual 0, modulus 1 |
| omega_crit | arccos(3/4) = 0.722 734 247 813 42 (K = -1) | 696 | exact: computed 0.7227342478134157 |
| Table 1: L = 0.24 | 0.643 501 108 793 29 and 0.795 398 830 184 14 | 697 | reproduced (0.64350110879328, 0.79539883018414: differences in the last printed digit) |
| Table 1: L = 0.245 | 0.667 526 397 108 77 and 0.774 680 351 224 54 | 697 | reproduced to all printed digits |
| Table 1: L = 0.249 | 0.698 494 191 441 32 and 0.746 325 490 506 20 | 697 | reproduced |
| Table 1: L = 0.25 | 0.722 734 247 813 42 (double) | 697 | exact (the numerical eigenvalues of A_crit give 0.722734245 and 0.722734251: the double-eigenvalue splitting, section 4) |
| Fig. 2 left | x~2 = -0.227 291 802 38, omega = 0.600 834 247 8 (L = 0.24) | 697 | reproduced: omega = 0.600834247816 |
| Fig. 2 right | x~2 = -1.071 257 417 2, omega = 0.798 309 349 426 (L = 0.24) | 697 | reproduced: omega = 0.798309349454 (2.8e-11) |
| settings | K = -1; L = 0.24, 0.245, 0.249, 0.25; 0.26 (manifolds, k = 1, 200, 270, 300); 0.28, 0.30; 0.27, 0.29, 0.31 (Fig. 7); 0.249, 0.25, 0.251 (tori from x_i = 0.001); 0.2501, 0.251 (inverse); error threshold 1e-12; h <= 1e-5 | 696-709 | settings |

## 4. Independent reproduction (COMPUTED)

1. Table 1 and the critical angle. Eigenvalue angles of the closed-form Jacobian A(K = -1, L): L = 0.24: 0.64350110879328, 0.79539883018414; L = 0.245: 0.66752639710877, 0.77468035122454; L = 0.249: 0.69849419144132, 0.74632549050620 (all equal to the printed 14-digit values, the
   printed ...793 29 against computed ...793 28 being the last-digit rounding); arccos(3/4) = 0.7227342478134157 equals the printed critical value. At L = 0.25 the numerical eigenvalues of A_crit are at 0.722734245 and 0.722734251: the double eigenvalue with a Jordan block splits to 3e-9, the sqrt of machine precision, as expected.
2. The critical orbit structure (closed form, so exact): A_crit has eigenvector residual 0 for the printed u + i v; singular values of A_crit - lambda I: 1.732, 1.379, 0.314 and 1.5e-16 (one null direction only, not two: a single Jordan block per eigenvalue); condition number of the eigenvector matrix 2.1e8.
3. Delta and multipliers across the transition, K = -1: Delta = 1 - 4L exactly (linear in L), with multiplier moduli (computed): L = 0.2501: 0.99247 and 1.007588; 0.251: 0.97639 and 1.024181; 0.26: 0.927512 and 1.078153; 0.27: 0.899471 and 1.111765, angle about 0.7228 to 0.7290.
   Just above the transition max |lambda| - 1 = 0.7588 sqrt(L - L_crit) (L = 0.2501), so the modulus excess scales like the square root of the distance while Delta scales linearly. Random perturbations of size eta added to A_crit move the eigenvalues off the unit circle by about 0.34 sqrt(eta) (median of 200 trials: eta = 1e-14: 3.4e-8; 1e-12: 3.5e-7; 1e-10: 3.7e-6).
4. Invariant curves (Fig. 2). Independent solver, not the authors' code: Fourier series with N = 30 (61 mesh points), unknowns the Fourier coefficients of all four components and omega, equations F(x)(theta_j) = Ts(x(theta_j)) - x(theta_j + omega) = 0 (244 equations) plus the two constraints x1(0) = 0 and x2(0) = x~2 (printed value imposed),
   Levenberg-Marquardt least squares, continuation in x~2 from 1e-3 amplitude along the linear eigenvector of lambda = exp(i 0.6435) (60 steps) for the first curve and of lambda = exp(i 0.7954) for the second. Results: the first curve reaches omega = 0.600834247816 with residual 4.9e-16 (printed 0.600 834 247 8), x~ = (0, -0.22729180, 0, -0.41453329);
   the second reaches omega = 0.798309349454 with residual 8.9e-16 (printed 0.798 309 349 426: difference 2.8e-11), x~ = (0, -1.07125742, 0, -0.41664847). The x~4 values are mine (not printed). So the printed 12-digit data of Fig. 2 are reproduced from the printed map and parameters, which also confirms the printed Ts.
   (A first attempt with omega fixed collapsed to the trivial curve x = 0 because the trivial solution satisfies F = 0 for any omega; the amplitude (x~2) must be the fixed quantity and omega the unknown; recorded because the same trap applies to any invariant-curve code.)

## 5. Reconciliation with the project

No project module computes invariant curves of a map or the Hopf unfolding. The project's Poincare-map machinery (`search/cr3bp_*`, `genome/qp_tori_energy_walk.py`, `search/qp_torus_fixed_jacobi_continuation.py`) is for the circular problem and quasi-periodic tori; I did not read them for this digest, so any overlap with the Fourier-Newton curve method above is not assessed. The maps Ts and Tt are self-contained and need no ephemeris: they are suitable as unit-test fixtures for any invariant-curve solver the project builds (`#916`, `#899` tori) and for the `#931` classifier (closed-form alpha, beta, Delta).

## 6. Techniques applicable to the project's problems

### 6.1 `#931`: classifier and critical tag, an exactly solvable sourced case
- The map Ts at K = -1 is an exactly solvable Hamiltonian-Hopf test: alpha = -3, beta = 4 + L, Delta = 1 - 4L (all printed), critical parameter L = 1/4, critical multiplier exp(+-i arccos(3/4)) (printed to 14 digits). Unit tests: the classifier's Delta sign change at L = 1/4 within 1e-12 (Delta is linear in L, so bisection converges exactly); the critical tag at L = 1/4 (Delta = 0 to rounding, |alpha/2| = 1.5 < 2, one null singular value of A - lambda I, eigenvector condition about 2e8); the printed Table 1 as angle tests on both sides (L = 0.24, 0.245, 0.249: both pairs on the unit circle with the printed angles to 1e-13).
- Tolerance near the transition, from section 4: the multiplier modulus excess is about 0.76 sqrt(L - L_crit) and rounding noise of the Jacobian alone is about 0.34 sqrt(eta); a modulus test with tolerance tol is blind for L - L_crit < (tol/0.76)^2 (about 1.7e-6 at tol = 1e-3), while Delta changes sign cleanly (1 - 4L). Use Delta for detection and the modulus only as a secondary check.
- Both the direct and the inverse maps share this Jacobian: the classifier cannot tell the two cases apart from the linear data; the side is decided by the nonlinear terms (Ts and Tt differ only by sin against tan). So a "critical" tag must carry "direct/inverse undetermined" unless a nonlinear test (stability of the bifurcating tori, as in the 2004 L4 paper) is run.
- The RTBP companion (`docs/notes/2026-10-05-digest-olle-pacha-villanueva-2004-l4-vertical-hopf.md`) gives the 3D restricted-problem test of the same detector; the two together cover a closed-form map (this paper) and a flow (that paper).

### 6.2 `#916`: persistence of stable members as invariant curves; tori born at the transition
- What this paper establishes numerically for the stable side (L < L_crit): two Lyapunov families of invariant curves exist for every L tested, Cantorian with holes below double precision, and near the origin the curves are elliptic. For the direct case these two families merge at L_crit and the merged family detaches and persists on the unstable side (elliptic, coexisting with the manifolds); for the inverse case the global family has hyperbolic parts that collapse to the origin at L_crit and no curves survive on the unstable side. For `#916`: a stable cycler member near a Delta = 0 event might lose all its nearby curves if the transition is inverse and keep a detached family if it is direct; determine which from a nonlinear test.
- Practical numbers for a project solver: Fourier truncation chosen automatically to reach an invariant-curve residual below 1e-12 (the 2004 L4 paper uses N up to 50); the printed Fig. 2 curves need N about 30 in my reproduction (residual 5e-16 to 9e-16); a first guess along the linear eigenvector with the amplitude as the continuation parameter works (section 4); the reducibility test by the spectrum of L_omega gives the normal eigenvalues.
- Positive control for an invariant-curve solver: the two Fig. 2 curves (map Ts, K = -1, L = 0.24, x1(0) = 0, x2(0) = -0.22729180238 gives omega = 0.6008342478; x2(0) = -1.0712574172 gives omega = 0.798309349426), tolerance 1e-9 on omega. Negative control: with omega fixed rather than the amplitude the solver returns the trivial curve (section 4); a solver should refuse a result whose Fourier amplitude is zero.
- Born tori: the paper's statement for the direct case is that the two Lyapunov families of elliptic 2D tori become one detached family, with normal eigenvalues exp(+-i nu) kept elliptic (Fig. 8); in the inverse case the stable-side tori are destroyed (all trajectories escape). The 2D-torus claims for the flow are in the L4 paper (stable 3D tori near the orbit change topology across the transition).

## 7. Recommended follow-ups (no task numbers registered)
1. Build a small test-only module with Ts and the Jacobian (closed form) and pin: Delta = 1 - 4L, alpha = -3, beta = 4 + L, Table 1, omega_crit = arccos(3/4), and the critical-tag test. No integrator needed.
2. Add the Fig. 2 invariant-curve reproduction (section 4, item 4) as the positive control for any Fourier-Newton solver, with the amplitude-fixed continuation recipe and the trivial-solution guard.
3. Reproduce a manifold iterate or the Lyapunov-family detachment (Fig. 7) on Ts to test a manifold or family-continuation routine; not done here.
4. Acquire Pfenniger 1985 (Astron. Astrophys. 150:97 and 112), Olle 2000 and Olle and Pfenniger 1999 for the earlier numerical Hopf studies on these maps (the L1 + L2 = 0 case), and Jorba 2001 (Nonlinearity 14:943) for the operator-spectrum method.

## 8. Summary
- The paper: numerical Hamiltonian-Hopf study on two explicit 4D symplectic maps Ts (direct) and Tt (inverse) with K = -1 and L_crit = 1/4, based on Fourier-Newton invariant curves, normal behaviour by operator spectra, and manifolds grown from the linear approximation.
- Closed forms printed and verified: Jacobian, alpha = -(4 + K1 + K2), beta, Delta = K(K + 4L), L_crit = -K/4, critical eigenvalue and eigenvector, omega_crit = arccos(3/4) = 0.72273424781342.
- Reproduced: Table 1 to 14 digits and both Fig. 2 invariant curves (omega = 0.600834247816 and 0.798309349454 against the printed 0.6008342478 and 0.798309349426) with an independent solver; the critical Jacobian has a single Jordan block per eigenvalue; Delta is linear and multiplier moduli scale like the square root of the distance to the transition, noise-amplified to about 0.34 sqrt(eta).
- Project use: an exactly solvable, fully printed test case for the `#931` Delta detector and critical tag, a positive control for invariant-curve solvers (`#916`), and the direct/inverse distinction that the linear data cannot resolve.
