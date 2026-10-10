# Digest: Jorba 2001, "Numerical computation of the normal behaviour of invariant curves of n-dimensional maps" (#960 batch 43)

Angel Jorba (Departament de Matematica Aplicada i Analisi, Universitat de Barcelona), Nonlinearity 14(5):943-976 (2001),
doi 10.1088/0951-7715/14/5/303 (Crossref checked 2026-10-11: title, vol. 14, issue 5, pp. 943-976, published 2 July 2001).
PII S0951-7715(01)18300-3. Received 25 Oct 2000, final form 3 May 2001. Recommended by L. Bunimovich. MSC 37C15.
- Proposed filename: `cyclers_pdf/papers/jorba-2001-numerical-computation-normal-behaviour-invariant-curves-n-dimensional-maps-nonlinearity-14-943-doi-10.1088-0951-7715-14-5-303.pdf`.
  - File as given: 35 pp. (IOPscience cover page, then journal pp. 943-976 in 34 pages), md5
    974de013e766688c7c2dc0db55bd3797. Publisher typesetting with a clean text layer. No OCR needed. PDF page n = journal page 942 + n.
- **How it was read:** the whole text through the text layer. The five printed tables (pp.956, 958, 962, 963, 972) and
  Figs. 1 and 4 were read on page images at 120 dpi. The images decide. The text layer agrees on every table cell.
  Table transcriptions are in `jorba-2001-normal-behaviour-tables.yaml` (216 cells, each also found in the text layer
  by script, `jorba-2001-yaml-verify.out`; that is a presence check, the row order was checked by eye on the images).
- Companion files in my folder: `jorba_checks.py/.out`, `jorba_compare.py/.out` (independent reproduction, sec. 3). All
  values in those files are ours.
- Held companion: `docs/notes/2026-10-05-digest-jorba-olle-2004-invariant-curves-hamiltonian-hopf.md`. Its follow-up item 4
  asked for this paper. Its sec. 2.3 summary of the method agrees with this paper (checked, sec. 4).

## 0. Verdict

**The method paper behind the normal-behaviour (stability) step of every "invariant curve of a Poincare map" computation
in the Barcelona school.** It is a numerical-analysis paper with proofs of the theory it needs and two worked examples.
- It computes the linear normal behaviour of an invariant curve of a diffeomorphism of R^n. The input is the curve as
  a Fourier series, its rotation number omega, and a routine for the map and its differential. The output is the normal
  eigenvalues (elliptic or hyperbolic), their eigenfunctions (the reducing transformation, and the first-order stable,
  unstable and centre directions), and a warning when the system is not reducible.
- It needs no symplectic structure and no closeness to integrable. It works on Poincare sections of ODEs (the section of
  the CR3BP or BCP at one period) because it only needs the map and its Jacobian.
- **For the torus lanes (`#916`, `#899`, `#996`) it gives:**
  - the complete algorithm for the stability of a stable cycler member viewed as an invariant curve of the one-period map
    (sec. 1);
  - the rule for picking the eigenvalue argument (it is defined only modulo multiples of omega), which any
    continuation along a family must handle (sec. 1.4);
  - a hard test of reducibility: eigenvectors whose Fourier coefficients do not decay (sec. 1.5);
  - pinned published numbers we can test against (sec. 2) and a reproduction of the Hill-equation ones to 1e-14 (sec. 3).
- **What it does not give:** a BCP test case we can pin. The BCP parameters (mu, m_S, a_S, omega_S) are not printed (sec. 2.3).
  Its examples are one linear skew-product (Hill with two forcing frequencies) and the BCP curves.
- **Project gap (checked, sec. 4):** no project module computes normal behaviour of an invariant curve. The torus code
  cites this paper in a docstring (`src/cyclerfinder/search/variational_qp_torus.py`, References) but uses it for the
  Fourier-collocation corrector idea only.

## 1. The method (READ)

### 1.1 Setting (sec. 1, pp.943-946)
- Map x_bar = f(x), x in R^n (n >= 2). Invariant curve x(theta), theta in T^1 of length 2 pi, with irrational rotation
  number omega: f(x(theta)) = x(theta + omega) (eq. 1). A non-autonomous map x_bar = f(x, theta), theta_bar = theta + omega
  has f(x(theta), theta) = x(theta + omega) (eq. 2).
- Curves without such a parametrisation exist (Broer, Simo and Tatjer 1998, ref [5]) and are excluded.
- **Computing the curve itself (sec. 1.1, brief).** Real Fourier series x(theta) = a0 + sum_{k=1..N} (a_k cos k theta + b_k sin k theta)
  (eq. 3). Mesh theta_j = 2 pi j/(2N+1), j = 0..2N (eq. 4). F(x)(theta) = f(x(theta)) - x(theta + omega) is sampled on
  the mesh and driven to zero by Newton, with the differential from the chain rule. The shift symmetry (x(theta + phi) is
  also a solution) makes the Jacobian singular. The cure is one extra linear row fixing a component at theta = 0. Gaussian
  elimination with row pivoting then leaves the last row as 0 = 0, which is dropped.
- **Linear normal behaviour (sec. 1.2).** x + h near the curve gives the quasi-periodic skew product x_bar = A(theta) x,
  theta_bar = theta + omega with A(theta) = D f(x(theta)) (eq. 5). No special matrix group is assumed.
- **Reducible (Definition 1.1):** a (maybe complex) change x = C(theta) y gives y_bar = B y with B = C^-1(theta + omega) A(theta) C(theta)
  independent of theta (eq. 6). The dynamics then follow the eigenvalues of B. This is the quasi-periodic Floquet theorem.
  The lack of uniqueness is the same as for periodic Floquet: arguments are defined only modulo multiples of omega (here)
  or 2 pi/T (periodic case).

### 1.2 Characterisation of reducibility by an eigenproblem (sec. 2, pp.946-950)
- Generalised eigenproblem A(theta) psi(theta) = lambda T_omega psi(theta), (T_omega psi)(theta) = psi(theta + omega) (eq. 7).
- Prop. 2.1: for an autonomous map, 1 is an eigenvalue with eigenfunction x'(theta) (the tangent to the curve).
- **Prop. 2.2:** if lambda is an eigenvalue so is lambda exp(i k omega) for every integer k (eigenfunction exp(-i k theta) psi).
  The closure of the spectrum is a union of circles.
- Prop. 2.3: in the reducible case, eigenvalues of B are eigenvalues of (7), and every eigenvalue of (7) is exp(-i k omega)
  times an eigenvalue of B for some k. So the set of eigenvalues is never empty.
- **Definition 2.1:** lambda1, lambda2 are omega-unrelated if lambda1 is not exp(i k omega) lambda2 for any integer k.
  Otherwise omega-related.
- **Prop. 2.4:** if (7) has n omega-unrelated eigenvalues, the system is reducible to B = diag(lambda_1, ..., lambda_n).
  Corollary 2.1: it cannot have more than n unrelated eigenvalues. So reducibility is equivalent to finding n unrelated
  eigenvalues with continuous eigenfunctions. The proof builds C(theta) = (psi_1, ..., psi_n) and shows it is regular
  (three steps, pp.947-948; the unrelatedness is what keeps psi_j(theta_0) independent at a point).
- log|lambda_j| are the Lyapunov exponents. The arguments measure the rotation about the curve, modulo omega. The
  reducing transformation is not unique: exp(-i k theta) C(theta) reduces to exp(i k omega) B.
- sec. 2.1: tangent versus normal spectrum. The tangent direction is always the eigenvalue 1 for an autonomous map. The
  author did not remove it from the eigenproblem (the saving, one dimension, "is not big enough to justify" the extra
  code). Mather 1968 (ref [31]) studies the full spectrum for Anosov maps. De la Llave (refs [29, 30]) shows the spectrum
  depends on the function space.
- sec. 2.2: eigenfunctions psi_nu = exp(-i nu theta) psi with non-integer nu are discontinuous periodic extensions. They
  would give any argument for the eigenvalue. Only continuous transformations are considered, so those are excluded.

### 1.3 Numerical approximation (sec. 3, pp.950-954)
- Discretise T_{-omega} o A(theta) as a matrix. For the Fourier coefficient vector x_N (2N+1 vector coefficients) compute
  the values f(x(theta_j)), then the Fourier coefficients of f(x(theta)) (eqs. 12-13). The differential of the first step is
  D_x f(x(theta_j)) times the identity times cos(k theta_j) or sin (chain rule). The Jacobian is (2N+1) x (2N+1) boxes of
  n x n blocks. T_omega is block diagonal with 2 x 2 rotations by k omega on each {cos k theta, sin k theta} pair. T_{-omega}
  is its transpose. The product is the matrix. Its eigenvalues and eigenvectors come from a shifted QR (Golub and Van Loan;
  Wilkinson). The matrix dimension is n (2N+1).
- A finite matrix always has eigenvalues, so a discretised spectrum can exist when the true operator has none. This is the
  non-reducible danger (sec. 5).
- **Accuracy of eigenvalues and which to keep (sec. 3.2, eq. 14).**
  - If (5) is reducible with eigenvalue mu_0 of B, the operator has all mu_k = mu_0 exp(i k omega). The discretisation
    contains only finitely many of them, with different accuracy.
  - Use the norm ||psi||^(p) = sum_j |psi_j| |j|^p (p = 1 is "enough"). The truncation error is
    TE(psi, N) = sum_{|j| > N} |psi_j| |j|^p. For exp(-i k theta) psi, TE is small only for a small set of |k|. So sort the
    eigenvectors by ||psi||^(p) and keep the smallest ones. These are the accurate eigenvalues. One keeps about 25 % (sec. 3.7).
  - A large norm of the representative means a larger N is needed.
- **Hyperbolic directions (sec. 3.3).** Eigenvalues with modulus not 1. Real eigenvalues: only one eigenvalue of each
  related class is real. The real eigenfunctions span the stable or unstable linear directions of the curve (the first-order
  invariant manifolds). Complex saddle (Hamiltonian): four eigenvalues lambda, 1/lambda, conjugate, 1/conjugate with |lambda|
  not 1. One picks a lambda and checks that the other three are eigenvalues.
- **Elliptic directions: the normal frequencies (sec. 3.4).** Modulus 1. Any member of a related class can serve. Choose the
  one with the minimum-norm eigenvector, which gives the simplest C. Alternatively use a priori information on the rotation
  number (the frequency-halving case in sec. 4.1).
- **Procedure (sec. 3.5).** Take the set of accurate eigenvalues. Identify the hyperbolic ones, choose a representative of
  each class, and delete the whole class (all related eigenvalues). Then do the elliptic ones the same way. Stop when the
  set is empty.
- **Resonances.** There is no non-resonance assumption. The paper says that, unlike the method of Jorba, Ramirez-Ros and
  Villanueva 1997 ([19], valid only close to constant coefficients and "not too close to resonances"), this method "can deal
  with these two situations" (p.955). The Prop. 2.4 hypothesis, n unrelated eigenvalues, is the only requirement. When the
  normal eigenvalues are omega-related (a normal-internal resonance) that hypothesis can fail. The paper does not treat
  that case (it cites Broer, Hanssmann, Jorba, Villanueva and Wagener 2001, ref [3], in preparation). INFERRED.
- **Symplectic case (sec. 3.6).** The method does not use symplecticity. For an autonomous symplectic map of R^6 under generic
  conditions the eigenvalue 1 has a 2 x 2 Jordan block. psi_1 is the tangent. psi_2 is the tangent to the one-parameter
  family of invariant curves. The off-diagonal entry delta is the derivative of the frequency along the family. A(theta) psi_2 =
  delta psi_1(theta + omega) + psi_2(theta + omega). Other eigenvalues come in pairs (lambda real or modulus 1) or fours.
  The author uses these relations only as accuracy checks.
- **Implementation (sec. 3.7).** C, double precision. Read the Fourier coefficients, build the matrix, shifted QR, normalise
  eigenvectors to sup norm 1, sort by ||.||^(p), drop large norms (keep about 25 %), check the last Fourier modes as the
  accuracy estimate and restart with more modes if too big, until the accuracy is reached or N is too large (the program
  then stops with a message).
- The discretisation the eigenproblem solves is the same as the one that computes the curve (sec. 1.1), so the Fourier
  coefficients of the curve and of A(theta) share one mesh.

### 1.4 The argument is ambiguous, and what to do about it (secs. 3.4, 4.1)
- The argument of a normal eigenvalue is defined only modulo multiples of omega. In a family the selected argument "does
  not have to be continuous with respect to the continuation parameter". Fig. 2 (bottom) and Fig. 6 show jumps. They are a
  change of determination, not physics. Adding a suitable multiple of omega (or of 2 pi omega2/omega1 in the Hill example)
  repairs the plot (Fig. 3 magnifies one jump near a^2 = 0.4).
- The argument is also "not necessarily the rotation number" (footnote 3, p.959). The rotation number is that value plus k omega.
- For a single torus this does not matter. For a family it does. Any project code that tracks a normal frequency along a
  family must unwrap by the family's omega.

### 1.5 Detecting non-reducibility (sec. 5, pp.966-972)
- Example (18): x_bar = A(theta) x, theta_bar = theta + omega, A = rotation by theta (Rychlik 1992 [34]). Prop. 5.1: not
  reducible by any continuous C. Proof by Fourier coefficients: the 1,1 entry of the reducing equation forces c_1,1^(k)
  e^{i k omega} lambda = c_1,1^(k-1). If a coefficient is non-zero then all are, and the series diverges or fails to decay
  unless |lambda_1| = 1, which then contradicts continuity. A geometric reason: the rotation after n steps is O(n^2), not O(n),
  so the rotation number of x0 is not defined (eq. 19).
- Prop. 5.2: if lambda is in the non-eigenvalue spectrum, so is every mu with |mu| = |lambda|. Prop. 5.3: the spectrum of
  (18) is the whole unit circle, with no eigenvalues. So the failure is "no eigenvalue", visible only through the
  eigenfunctions.
- **Numerical signature (sec. 5.3, Table 5).** N = 16 gives 66 eigenvalues with modulus 1 to 1e-15 and arguments within
  1e-15 of integer multiples of 2 pi/(2N+1) = 0.190 399 554 763 017 772 02. The argument is then "totally undefined" and is
  chosen by the truncation. The eigenvector Fourier coefficients are all the same size, as if divergent. Use: **a computed
  normal "eigenvalue" whose eigenvector has non-decaying Fourier coefficients, or whose argument sits on the lattice
  k 2 pi/(2N+1), is not a true eigenvalue.**
- sec. 5.4 and Fig. 10: the Hill equation far from constant coefficients (a^2 = 0.408, b = 0.6) is non-reducible by Broer and
  Simo's earlier computation. The computed eigenvalues cluster on circles. The outer radius matches exp(Lyapunov exponent) =
  exp(0.2967086) = 1.3454232 (printed, from Broer-Simo), so "the spectral radius of the operator seems to coincide with the
  maximal Lyapunov exponent" even when not reducible. The paper says a deeper study is in progress. Figure 10 uses N = 32,
  64, 128 and 256.

### 1.6 Error estimates, in one place
- Truncation error of an eigenfunction: TE(psi, N) (eq. 14), estimated in practice by the size of the last Fourier
  coefficients. Table 1 and Table 2 print this as the "Error" column (the norm of the last two Fourier coefficients).
- Automatic N: the program starts at N = 8 and raises N until the estimate reaches the requested accuracy (10^-12 in the
  examples).
- Eigenvalues of an eigenfunction with large norm are poor. Table 1 shows errors from 1.3e-14 (norm 10.7) to 2.6 (norm 91).
  Rows with error of order 1 (the arg 0 rows, and the last two) are not true eigenvalues.
- Self-consistency: the symplectic relations are not imposed, so their residual measures accuracy: Table 3 gives
  |lambda1 lambda2 - 1| about 4e-15; Table 4 gives |(alpha1 + 2 omega) - alpha3| about 1e-14.
- Manifold step: with the linear direction h psi_1, the second-order error is estimated by
  E(h) = max_theta || f(x(theta) + h psi_1(theta)) - x(theta + omega) - h lambda_1 psi_1(theta + omega) ||_2. For the BCP
  curve h = 10^-7 gives E(h) < 10^-13 (p.964).
- Cost: for N = 14 the Hill scan needs 29 evaluations of A(theta2) and a dimension-58 eigenproblem (p.958). For higher
  tori the number of harmonics grows roughly exponentially in the torus dimension r (sec. 6, footnote 5). The
  discretisation is fully parallel. The tangent eigendirections can be removed when the torus has low codimension.

## 2. The worked examples and every printed number

### 2.1 Quasi-periodic Hill equation (sec. 4.1, pp.954-959)
- System (eq. 15): x' = y, y' = -(a^2 + b p(theta1, theta2)) x, theta1' = omega1, theta2' = omega2, with
  p = cos theta1 + cos theta2. It is the normal behaviour of the 2-torus x = y = 0, which exists for all parameters.
- Section theta1 = 0. A(theta2) is the 2 x 2 fundamental matrix at time 2 pi/omega1. The cocycle is z_bar = A(theta2) z,
  theta2_bar = theta2 + 2 pi omega2/omega1 (eq. 16). The map is linear but non-autonomous, so Prop. 2.1 does not apply
  and there is no eigenvalue 1. A(theta2) came from a 16th-order Taylor integrator, local error 1e-16 (ref [25]).
- **Example 1 (p.955).** a^2 = 1, b = 0.1, omega1 = sqrt 2, omega2 = sqrt 3 (the case of Jorba, Ramirez-Ros and Villanueva
  1997, ref [19]). Requested accuracy 1e-12 forces N = 8 (17 values of theta2, dimension-34 eigenproblem). N = 8 gives error
  below 1e-14. Result: reduced matrix has two eigenvalues of modulus 1 with argument nu = 1.857 169 546 872 597. The flow
  frequency from [19] is alpha = 0.996 203 551 143 221. Relation alpha 2 pi/omega1 + nu = 2 pi (mod 2 pi) holds to about 0.8e-15
  (printed). **Ours: the printed 15-digit values give 2.7e-15 for that residual** (rounding of the printed digits).
- **Example 2 (p.956, Tables 1 and 2, Fig. 1).** omega1 = 1, omega2 = gamma = (1 + sqrt 5)/2. The cocycle frequency is
  omega = 4 pi/(1 + sqrt 5) = 3.883 222 077 450 9332 (printed). b = 0.2, a^2 = 0.7. N = 8 gives error estimate 5.6e-12, above
  the target, so the program repeats with N = 10 and reaches 1.3e-14. Table 1 (N = 10) lists the full spectrum. The leading
  pair has modulus 1.421 399 533 740 337 and 0.703 531 960 059 501, argument 1.199 981 614 864 327, which "seems to coincide with"
  pi - omega/2 = 1.199 981 614 864 3267. Here modulus is not 1: this is hyperbolic (INFERRED label: "instability pocket" is my wording, from the
  Broer-Simo title), Lyapunov exponent
  log 1.4214 = 0.3516 (ours). Fig. 1 plots the eigenvalues in the plane for N = 10 (top) and N = 64 (bottom). Two circles of
  inverse radius appear, and a few outlying points are inaccurate eigenvalues.
- **Frequency halving (Table 2, eq. 17).** The argument pi - omega/2 shows the system cannot be made real on T. Work on the
  2-torus: A(2 eta) psi(eta) = lambda psi(eta + omega/2). N = 16 is needed. The first two eigenvalues are now real and negative
  (modulus 0.703 531 960 059 501 and 1.421 399 533 740 336, argument 3.141 592 653 589 793), so the system reduces to a real
  constant matrix.
- **Scan (p.958, Fig. 2).** b = 0.2, a^2 from 0.25 to 1.00, N between 8 and 14 automatically. The top plot (log of the unstable
  eigenvalue modulus) is "identical" to Fig. 4 (left) of Broer and Simo 1998. The bottom plot is the argument on [0, pi]. A
  jump near a^2 = 0.4 is a change of determination (Fig. 3).
- Printed numbers: all in the yaml (Tables 1 and 2, example 1, golden frequency).

### 2.2 Bicircular problem (sec. 4.2, pp.959-966)
- BCP: Earth and Moon on circular orbits about their barycentre, the barycentre on a circular orbit about the Sun, all coplanar.
  Rotating frame with Earth and Moon fixed. Hamiltonian (p.961): H = (px^2 + py^2 + pz^2)/2 + y px - x py - (1-mu)/r_PE - mu/r_PM -
  m_S/r_PS - (m_S/a_S^2)(y sin theta - x cos theta), r_PE^2 = (x-mu)^2 + y^2 + z^2, r_PM^2 = (x-mu+1)^2 + y^2 + z^2,
  r_PS^2 = (x-x_S)^2 + (y-y_S)^2 + z^2, x_S = a_S cos theta, y_S = -a_S sin theta, theta = omega_S t. First introduced by
  Cronin, Richards and Russell 1964 ([8]). **The values of mu, m_S, a_S and omega_S are not printed in this paper.**
- Near each triangular point the BCP has three periodic orbits (period of the Sun). The out-of-plane linear oscillation
  continues to a Cantorian one-parameter family of two-frequency tori. The map is the time-of-the-Sun section P, a symplectic
  map of R^6. Its differential is obtained by integrating the variational equations.
- Fig. 4: rotation number against zdot (value of zdot at the point of the curve with z = 0, zdot > 0). Read on the image: axes
  zdot 0 to 0.9, rotation number 0.505 to 0.55. VF1 and VF2 meet at a turning point near zdot about 0.2 to 0.25. VF3
  continues to zdot about 0.9 where the rotation number falls to about 0.505. The holes of the Cantor structure are too small to
  see. These are approximate readings, not goldens.
- VF1 contains hyperbolic tori, VF2 and VF3 are mostly elliptic. Notation of Castella and Jorba 2000 ([6]).
- **Table 3 (p.962):** unstable curve of VF1 near the VF1-VF2 transition: omega = 0.535 033 339 385 478, zdot = 0.080 508 698 608 030.
  Normal eigenvalues (the two neutral directions, tangent and conjugate, skipped): saddle x centre.
  - lambda1: modulus 1.091 942 641 437 887, arg 0.
  - lambda2: modulus 0.915 799 019 152 856, arg 0.
  - lambda3, lambda4: modulus 0.999 999 999 999 985, arg +/-2.035 517 841 801 725.
  Printed |lambda1 lambda2 - 1| about 4e-15. Ours from the printed digits: -4.33e-15.
- Fig. 5: a pair of eigenvalues near the VF1-VF2 transition. It goes from hyperbolic to elliptic by a collision at 1 ("an
  analytical study of some of these transitions" is ref [3]). On Cantor sets, so the effect is invisible at the plot scale except
  for low-order resonances.
- **Table 4 (p.963):** unstable curve of VF3 (narrow hyperbolic set): omega = 0.505 512 051 568 360, zdot = 0.895 865 023 909 962.
  - lambda1, lambda2: modulus 1.000 103 145 715 284, arg +/-0.436 023 756 155 351.
  - lambda3, lambda4: modulus 0.999 896 864 922 655, arg +/-1.447 047 859 292 059.
  Printed |(alpha1 + 2 omega) - alpha3| about 1e-14: a complex saddle (the four eigenvalues lambda, 1/lambda, conjugates).
  Ours: 1.2e-14, and |lambda1||lambda3| = 0.9999999999999977.
- Fig. 6: arguments of the two elliptic eigenvalues along VF3, with the unsmoothed jumps shown on purpose.
- **Unstable manifold (sec. 4.2.1, Figs. 7-9).** The curve of Table 3. Linear unstable manifold x(theta) + h psi_1(theta), h = 10^-7
  (E(h) < 1e-13). C_1 is its image, C_j = f(C_{j-1}). Fig. 8: (x, y) projections of the invariant curve with C_j for j = 1..100
  (a); j = 120..240 (b); C_330 (c); C_400 (d). The manifold bends and comes back near the curve (as for the periodic-orbit case close
  to a stability change). The stable and unstable manifolds "do not coincide but intersect". By symmetry there are four
  homoclinic connections, not two (Fig. 9: distance between the curve and C_330 against the curve parameter). Iterates are
  trustworthy for j up to 400. Errors show (spikes, mesh non-uniformity) only after the second passage near the curve.
  The paper describes this as a crude manifold method, for illustration, and a building block for heteroclinic chains and
  Arnold diffusion between tori.

## 3. Independent reproduction (OURS, `jorba_checks.py/.out`, `jorba_compare.py/.out`)

Own implementation: monodromy matrices by DOP853 (rtol 1e-13), complex Fourier collocation of the operator, numpy eig.
Not the author's code. The norm column and error column use my own normalisation and do not reproduce. Moduli and
arguments do.
1. **Table 1 (golden-mean Hill, a^2 = 0.7, b = 0.2, N = 10): reproduced.** Ours modulus 1.4213995337403 and 0.70353196005949 against
   printed 1.421 399 533 740 337 and 0.703 531 960 059 501 (max difference 1.4e-14 over six eigenvalues). Arguments, ours against
   printed: 1.199981614864328 (diff 6.7e-16), 2.683240462586605 (-1.8e-15), 0.283277232857952 (-1.1e-16), 2.116685996870701
   (1.3e-15), 1.766536080580230 (-2.2e-15), 0.633427149148422 (2.4e-15). The last two have printed error estimates of 1.8e-8
   and 5.8e-7, so agreement far below their error is expected when the discretisation is the same as the author's.
2. **Table 2 (N = 16, eq. 17): reproduced.** Argument pi, 1.199981614864, 0.741629423861, 2.683240462587, moduli 0.70353 and 1.42140. The
   half-frequency problem needs the reduced cocycle frequency 3.8832 (not 2 pi gamma = 10.1664) because the choice of
   omega representative changes the sign of lambda in the halved problem. This is a trap for any implementation: with 2 pi
   gamma I got argument 0 instead of pi.
3. **Example 1: reproduced.** nu = 1.857169546872602 (+5e-15 from printed). The partner arguments are 3.0139 and 0.4451 (ours).
4. **Table 5 pattern: reproduced.** For system (18) the discretised operator, with any irrational omega, gives moduli 1 +/- 4e-15
   and arguments within 1.7e-15 of multiples of 2 pi/33. Omega for (18) is not printed, so the individual digits (k 2 pi/33) are
   an arithmetic fact: 2 pi/33 = 0.19039955476301776, and the printed table values equal k 2 pi/33 to 1e-15.
5. **Arithmetic relations in the printed numbers:** omega = 4 pi/(1 + sqrt 5) = 3.8832220774509327 (printed ...9332); pi - omega/2 =
   1.1999816148643268 (printed ...3267); 1.421399533740337 x 0.703531960059501 = 1.0000000000000000 (the pair is lambda, 1/lambda).
6. **Not reproduced: Tables 3 and 4.** They need the BCP constants (not printed). Only the internal relations were checked (sec. 2.2).
7. Cost note (ours, arithmetic only): the dense eigenproblem has dimension n (2N+1). For the BCP section (n = 6) with N = 30 to 50
   that is 366 to 606, which is small for a dense eigensolver today. (The paper does not print the N used for Tables 3 and 4.)

## 4. Reconciliation with the project and the Jorba-Olle digest

- **Jorba-Olle 2004 digest (held):** its sec. 2.3 summarises Proposition 5.1 (eigenvalue lambda gives lambda exp(i k omega)), Definition 5.1
  (omega-unrelated), Proposition 5.2 (reduction to a diagonal B), the spectrum of the truncation by shifted QR, and
  reducibility to 1e-12. These are this paper's Prop. 2.2, Def. 2.1 and Prop. 2.4 (the 2004 paper recaps them under its own numbers). The substance agrees. This paper adds
  the proofs, the selection of eigenvalues by eigenvector norm (the 2004 digest does not mention it), the 25 % rule, the
  frequency-halving device, and the non-reducibility test. Its follow-up item 4 ("acquire Jorba 2001 for the operator-spectrum
  method") is now met.
- **Code:** `variational_qp_torus.py` cites this paper (References) and builds the Fourier-collocation corrector (the sec. 1.1
  idea). A grep for "normal", "reducib", "eigen" in that file and in `variational_crnbp_torus.py`, `crnbp_torus_ghost_guard.py`
  finds no normal-behaviour solver. `crnbp_torus_ghost_guard.py` says the torus case has "no eigenvector direction". So the
  stability step of a torus (the normal eigenvalues, the reducibility test, the omega-unrelated bookkeeping) is not in the project.
  I did not read the code beyond those greps.
- Reusing the Jorba 2001 normal-behaviour step on top of the project's torus correctors needs only the Jacobian of the map along the
  curve. The project's variational torus modules already carry variational vectors (`#924`), so the input exists.

## 5. Techniques for the project (PROPOSALS only, no task numbers registered)

1. **`#916` (persistence of stable cycler members as invariant curves).** Use the method as the stability test of each curve: the
   curve is stable if all normal eigenvalues have |lambda| = 1, after removing the tangent and the
   conjugate Jordan pair. Two checks stand in for the missing symplectic imposition: the product lambda1 lambda2 = 1 residual and
   the pair structure, as in Tables 3 and 4 (4e-15 and 1e-14 as the achieved level).
2. **Selecting the elliptic representative.** Use the minimum eigenvector norm ||psi||^(1) with sup-norm scaling, and unwrap the
   argument against omega along a family before comparing to a resonance condition. Never compare a raw argument with
   k omega across two family members.
3. **Reducibility guard.** Refuse to report a normal frequency when the eigenvector Fourier coefficients do not decay, or when
   the argument lands on the lattice k 2 pi/(2N+1). This is the paper's own warning (sec. 5.3).
4. **Test fixtures (published values only).** (a) Example 1: a^2 = 1, b = 0.1, omega1 = sqrt 2, omega2 = sqrt 3, nu = 1.857169546872597
   (tolerance 1e-12). (b) Golden Hill, a^2 = 0.7, b = 0.2: leading modulus pair 1.421399533740337 and 0.703531960059501, argument
   1.199981614864327 (tolerance 1e-11 for N = 10; the printed digits beyond 1e-13 depend on the discretisation). (c) The N = 16
   half-frequency case: first two real negative eigenvalues of modulus 0.703531960059501 and 1.421399533740336. (d) The Table 5
   lattice property as a negative control (non-reducible cocycle). Our `jorba_checks.py` is a reference implementation for (a) to (d).
   The BCP Tables 3 and 4 cannot be fixtures until the BCP constants are found (Castella and Jorba 2000 [6], not held).
5. **Positive control for any invariant-curve plus normal-behaviour pipeline** combining this paper and Jorba-Olle 2004: the
   Hill example here for the normal-behaviour step, and the Ts map of Jorba-Olle for the curve step (held digest sec. 3, 4).
6. **Manifold seeds.** The linear unstable direction h psi_1 with h chosen so E(h) < 1e-13 is the seed for tube or manifold
   computations from a torus (`#899`, `#996`). Expect the pattern in sec. 2.2: usable iterates until the second return near the
   curve.

## 6. Citation mining (references 1-38, p.975-976)

Checked 2026-10-11 with `ls cyclers_pdf/papers | grep -i` for each author and with CORPUS_INDEX.md and the wanted list
(row 31 is this paper). The `Not held` entries below carry the printed reference-list form; DOIs are given only where Crossref returned the record.
- **Held:**
  - [22] Jorba and Villanueva 1997, J. Nonlinear Sci. 7:427-473 (`jorba-villanueva-1997-persistence-lower-dimensional-invariant-tori-...pdf`).
  - Jorba-Olle 2004 (the follow-up paper; `jorba-olle-2004-...pdf`).
- **Not held, relevant to the torus lanes (PROPOSALS for the wanted list, in priority order):**
  1. [6] Castella and Jorba 2000, "On the vertical families of two-dimensional tori near the triangular points of the bicircular
     problem", Celest. Mech. Dyn. Astron. 76:35-54, doi 10.1023/a:1008321605028 (Crossref verified 2026-10-11). Source of the BCP constants and the VF1/VF2/VF3 notation. It would let Tables 3
     and 4 be reproduced. No wanted-list row found.
  2. [21] Jorba and Villanueva 1997, "On the normal behaviour of partially elliptic lower dimensional tori of Hamiltonian systems",
     Nonlinearity 10:783-822, doi 10.1088/0951-7715/10/4/001 (Crossref verified 2026-10-11). The theory of the normal behaviour that this paper makes numerical. No row found.
  3. [3] Broer, Hanssmann, Jorba, Villanueva and Wagener, "Normal-internal resonances in quasiperiodically forced oscillators: a
     conservative approach", Nonlinearity 16:1751-1791 (2003), doi 10.1088/0951-7715/16/5/312 (Crossref verified 2026-10-11; the
     "in preparation" paper of ref. [3]). It treats the omega-related resonant case this paper does not.
  4. [4] Broer and Simo 1998, "Hill's equation with quasi-periodic forcing: resonance tongues, instability pockets and global
     phenomena", Bol. Soc. Bras. Mat. 29:253-293 (1998), doi 10.1007/bf01237651 (Crossref verified 2026-10-11). Source of the Hill-equation reference values (the Lyapunov exponent 0.2967086).
  5. [17] Jorba 2000, "A numerical study on the existence of stable motions near the triangular points of the real Earth-Moon
     system", Astron. Astrophys. 364:327-338 (as printed in this paper's reference list; a Crossref bibliographic query
     did not return it, so no DOI is verified), and [16] Jorba 2000 (HAMSYS-98 proceedings, 197-213, as printed). Practical-stability
     application of the BCP tori.
  6. [37] Simo, Gomez, Jorba and Masdemont 1995 (as printed), "The bicircular model near the triangular libration points of the RTBP", in From
     Newton to Chaos (Plenum), 343-370.
  7. [15] Jorba 1999 (as printed), "A methodology for the numerical computation of normal forms, centre manifolds and first integrals of
     Hamiltonian systems", Exp. Math. 8:155-195, and [23] Jorba and Villanueva 1998, Physica D 114:197-229 (the normal-form side).
- Not held, lower priority (background): [19] Jorba, Ramirez-Ros and Villanueva 1997, SIAM J. Math. Anal. 28:178-188 (effective
  reducibility); [20] Jorba and Simo 1992, J. Diff. Eq. 98:111-124; [8] Cronin, Richards and Russell 1964, Icarus 3:423-428 (origin of
  the BCP); [36] Simo 1998 CISM 387 (effective computations in celestial mechanics); [35] Simo 1999 (figure-eight; high-precision
  Poincare map expansion); [10]-[12] Eliasson; [26]-[28] Krikorian (reducibility on compact groups); [29]-[31] de la Llave, Mather
  (spectrum of the operator); [18] Jorba, de la Llave and Zou (Lindstedt series); [5] Broer, Simo and Tatjer 1998 (curves without a
  parametrisation); [1], [2], [7], [9], [13], [14], [24], [25], [32], [33], [34], [38].
- The wanted list holds a row for this paper (row 31). With it filed, PROPOSAL: remove row 31.
- A grep of CORPUS_INDEX.md for Castella, Broer, Eliasson, Cronin, Simo (C.) found none of the above. The corpus index has no
  Jorba-Ollé companion rows for them either.

*Filed as `cyclers_pdf/papers/jorba-2001-numerical-computation-normal-behaviour-invariant-curves-n-dimensional-maps-nonlinearity-14-943-doi-10.1088-0951-7715-14-5-303.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/jorba-2001-normal-behaviour-tables.yaml`.*
