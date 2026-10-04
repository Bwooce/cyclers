# Digest: Olle, Pacha and Villanueva 2005, normal form around a non-semi-simple 1:-1 resonant periodic orbit

Date: 2026-10-05 (Sydney). Reading, reasoning and small independent computations; no project code was changed.

Source: Merce Olle, Juan R. Pacha and Jordi Villanueva, "Quantitative estimates on the normal form around a non-semi-simple
1:-1 resonant periodic orbit", Nonlinearity 18:1141-1172 (2005), DOI 10.1088/0951-7715/18/3/012, received 9 July 2004, in final
form 9 November 2004 (Universitat Politecnica de Catalunya). Filed in the private paper corpus as
`olle-pacha-villanueva-2005-normal-form-non-semi-simple-1-1-resonant-periodic-orbit-nonlinearity-18-1141-doi-10.1088-0951-7715-18-3-012.pdf`.
A 33-page PDF (a publisher cover page, then journal pp.1141-1172) with a good text layer. I read the whole text layer, and
checked the theorem page (p.1143), the constants page (p.1169) and the figure page (p.1171) on page images, because the text layer
loses the Greek epsilon (the invariant sign) and the subscript Omega_s.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my own computation (scratch scripts, not committed; the recipes are in
section 6). INFERRED is my reasoning.

## 0. The answer to the questions that were asked, first

1. There is NO example system, no table, no critical orbit, no multipliers, no energies and no mu in this paper. It is a pure
   analysis paper. READ: its only numerical statements are the value of one proof constant (lambda-tilde, p.1169, remark A.6) and the
   plot of one auxiliary sequence (Figure A1, p.1171). The restricted three-body problem appears only in the reference list: Bruno 1994
   (The Restricted 3-Body Problem), Jorba and Villanueva 1998 (Physica D 114, numerical normal forms around periodic orbits of the
   restricted three-body problem, cited as the numerical implementation of an elliptic-case normal form, p.1143), Giorgilli et al. 1988,
   and Olle, Pacha and Villanueva 2004 (Cel. Mech. Dyn. Astron. 90:89, "Motion close to the Hopf bifurcation of the vertical family of
   periodic orbits of L4", the application to the restricted problem, cited as [17] and in the Introduction, not held). So the answer to
   "is the example the restricted three-body problem" is: the paper has no example; its application elsewhere is the vertical family of
   periodic orbits of L4 in the spatial restricted problem (by the title of reference [17]), mu not given in this paper.
2. Hence there is no printed critical orbit to use as a sourced test. The sourced numbers are few (section 4). A usable test case for the
   Delta < 0 transition has to come from elsewhere: the project's own Hadjidemetriou 1975b Table I has a collision of two unit-circle pairs
   between rows 10 and 11 (section 5.1), and a synthetic normal-form case can be built from this paper's own normal form (section 6).
3. What the paper does give that matters for the `#931` classifier is the exact algebraic character of the transition (non-diagonalizable
   normal block, multipliers {lambda, lambda, 1/lambda, 1/lambda}) and a precise statement that near it the normal form cannot be
   improved to exponential smallness (the remainder is only R^(r_opt/2) with r_opt growing like a Lambert function): the analytic expression of why
   there is no sharp stability radius around such an orbit.

## 1. Setting (READ pp.1141-1143)

A real analytic three-degree-of-freedom Hamiltonian, a periodic orbit with four non-trivial characteristic multipliers (those not equal to 1)
that collide pairwise on the unit circle at two conjugate points different from +-1, equivalently two normal frequencies that are equal.
The Introduction places this as the transition orbit of a one-parameter family (families of periodic orbits are parametrised by the energy
when no non-trivial multiplier equals 1, citing Siegel and Moser 1971): from "stability (the four non-trivial characteristic multipliers
are different and of modulus one) to complex instability (one non-trivial characteristic multiplier is a complex number of modulus different
from one and the other three are the complex conjugate number and the corresponding inverse numbers) through a pairwise collision of its
non-trivial characteristic multipliers in the unit circle (at two conjugate points different from +-1)", "usually referred to as the
quasi-periodic Hamiltonian Andronov-Hopf bifurcation" (p.1142; reference [18], Pacha's thesis, for a broad study). In the project's
vocabulary (Hadjidemetriou 1975b digest, section 2.3): two unit-circle pairs with b1 = b2 and |b| < 2, that is Delta = 0 with
|b| < 2, the boundary between "Delta > 0, |b1|, |b2| < 2" (stable) and "Delta < 0" (complex quartet).
"Common in several mathematical models of science"; the paper deliberately does not name the contexts.

Two classifications (READ p.1142): semi-simple against non-semi-simple (whether the 4 x 4 box of the monodromy matrix on the normal
directions is diagonalizable); the generic case is non-semi-simple and it is the one treated. Rational against irrational
omega1/omega2 (the resonance); the paper takes it irrational and, for quantitative work, Diophantine:
|<k, omega>| >= gamma |k|_1^(-tau) for all k in Z^2 minus {0}, gamma > 0, tau >= 1 (eq. 3).

Coordinates: (theta, x, I, y) in T^1 x R^2 x R x R^2 with 2-form dtheta^dI + dx^dy; theta describes the orbit, which is the circle
{I = 0, x = y = 0}. After a Floquet transformation (symplectic, linear in (x, y), 2 pi periodic in theta) the Hamiltonian is (eq. 1, p.1142)

    H(theta, x, I, y) = omega1 I + omega2 (y1 x2 - y2 x1) + epsilon (y1^2 + y2^2)/2 + Hhat(theta, x, I, y),

omega1 the orbit's angular frequency, omega2 its (only) normal frequency, so the non-trivial multipliers are {lambda, lambda, 1/lambda,
1/lambda} with lambda = exp(2 pi i omega2/omega1). The sign epsilon = +-1 is an invariant of the collision; a change of time t -> -t
exchanges the two contexts, and the paper takes epsilon = +1. The term epsilon(y1^2 + y2^2)/2 is the nilpotent part that makes the normal
variational equations non-diagonalizable. (The text layer drops the epsilon in eq. 1; I read it on the page image of p.1142-1143 and in
the later formula (20).)

## 2. The normal form and the main theorem (READ pp.1143-1144)

Domain (eq. 2): D(rho0, R0) = {|Im theta| <= rho0, |I| <= R0^2, |(x, y)| <= R0}. The adapted degree of a monomial I^l q^m p^n exp(ik theta)
is 2l + |m|_1 + |n|_1 (I counts twice), the weighted norm (eq. 9) is ||f||_{rho,R} = sum |f_klmn| R^(2l+|m|+|n|) exp(|k| rho).

Theorem 1.1 (p.1143-1144), epsilon = +1, H analytic in D(rho0, R0) with finite weighted norm and omega Diophantine (3). Given any
epsilon' > 0 and sigma > 1 there is R* in (0, 1), depending on rho0, R0, ||H||, |omega1|, |omega2|, gamma, tau, epsilon' and sigma, such that for
0 < R <= R* there is a real analytic canonical diffeomorphism Phi^(R) with
(i) defined on D(sigma^-2 rho0/2, R), mapping into D(rho0/2, sigma R);
(ii) Phi - Id has components 2 pi periodic in theta with sup bounds: (1 - sigma^-2) rho0/2 for the theta component, (sigma^2 - 1) R^2 for I,
     (sigma - 1) R for each of X_j and Y_j;
(iii) H o Phi = Z^(R)(x, I, y) + R^(R)(theta, x, I, y), with the integrable normal form

     Z^(R) = omega1 I + omega2 (y1 x2 - y2 x1) + (y1^2 + y2^2)/2 + Zbar^(R)( I, -(x1^2 + x2^2)/2, (y1 x2 - y2 x1)/2 ),

     Zbar a polynomial of degree at most floor(r_opt(R)/2) in its three arguments starting at degree two;
(iv) r_opt(R) := 2 + exp( W( log( 1 / R^(1/(tau + 1 + epsilon')) ) ) ), W the Lambert function, W(z) exp(W(z)) = z, W: (0, inf) -> (0, inf);
(v) ||Z^(R)||_R <= ||H||_{rho0, R0} and ||R^(R)||_{sigma^-2 rho0/2, R} <= R^(r_opt(R)/2);
(vi) R^(R) goes to zero faster than any power of R.
Remark 1.2: the Lebesgue measure of the omega in R^2 failing (3) for some gamma is zero for any tau > 1 (Lochak and Meunier 1988, appendix 4).

So the normal form is a function of the three combinations I, x1^2 + x2^2 and y1 x2 - y2 x1, plus the fixed linear part with the
nilpotent term: it is integrable. The remainder is smaller than any power of R but NOT exponentially small (comparison below).

Contrast with the elliptic case (READ pp.1144-1145): for omega = (omega1, omega2, omega3) Diophantine with tau >= 2 a normal form
Z^(r)(I, (x1^2 + y1^2)/2, (x2^2 + y2^2)/2) exists to any order r, and with r = r_opt = b/R^(2/(tau+1)) the remainder is at most
a exp(-b/R^(2/(tau+1))) (Jorba and Villanueva 1997, reference [11]); remark 8.2 (p.1165) gets exp(-(c/R)^(1/(tau+1))) by the
Giorgilli-Galgani route, matching Delshams and Gutierrez 1996. The reason for the difference (p.1145): in the non-semi-simple case the
homological equations are lower triangular with a Jordan-type block of size growing with the order, with the small divisor on its diagonal,
so the solution at order s carries the factor 1/Delta^(M+N+1) (with fact(M+N) in front); in the elliptic case only a single power of one small
divisor appears. The authors "do not claim that the estimates are optimal" but are convinced they cannot be strongly improved with the
standard approach (remarks 7.2, A.9).

## 3. How the normal form is constructed and where the estimates come from (READ sections 3-9)

Method: Lie series by the Giorgilli-Galgani algorithm T_G (definition 3.1, eq. 13), generating function G = sum G_s, homological equation
L_{H2} G_s + Z_s = F_s (eq. 16) with F_s from the recursion (15). Complex coordinates (eq. 17): x1 = (q1 - p2)/sqrt2, x2 = i(q1 + p2)/sqrt2,
y1 = (q2 + p1)/sqrt2, y2 = i(q2 - p1)/sqrt2 turn the quadratic part into H2 = omega1 I + i omega2 (q1 p1 + q2 p2) + epsilon q2 p1 (eq. 20); the
transformed Hamiltonian is complex analytic with the symmetry S (eq. 12), preserved by Poisson brackets, real values corresponding to
qbar1 = -p2, qbar2 = p1.

Homological operator (eqs 24-25): L_{H2} acting on a monomial alpha = I^l q^m p^n exp(i k theta) gives (Delta + m1 q2/q1 - n2 p1/p2) alpha ... with the
scalar Delta = Delta_{k,|m|,|n|} = i omega1 k + i omega2 (|m|_1 - |n|_1) (the part of the printed formula that the text layer garbles reads: L_{H2} alpha =
(Delta + epsilon (m1 q2/q1 - n2 p1/p2)) alpha ... I checked the structure, not every sign, on the page image). Consequences:
- Delta = 0 exactly for k = 0 and |m|_1 = |n|_1 (irrationality of omega1/omega2); those monomials are the possible non-removable (resonant) terms, only in even degree s.
- For the others (Delta not equal to 0, space E^+), the equation is a block lower-triangular system with Delta on the diagonal (eq. 27, matrix of size (M+1)(N+1), blocks
  Id times j and D_N = Delta Id - P_N with P_N nilpotent); explicit inverse (eq. 43-44), and Lemma 5.3: |Lambda^-1|_1 <= (1 + 1/|Delta|)^(M+N) fact(M+N)/|Delta|, "quite sharp"
  (remark 5.4: a summand fact(M+N)/|Delta|^(M+N+1) occurs).
- On the resonant block (Delta = 0, k = 0, M = N) L has trivial kernel but non-trivial cokernel; the non-removable terms Zhat_i xi1^i xi3^(M-i) are given in closed form (eq. 38) in the variables
  xi1 = q1 p2, xi2 = q2 p1, xi3 = i(q1 p1 + q2 p2)/2, xi4 = (q1 p1 - q2 p2)/2 with {xi1, H2} = -2 xi4, {xi2, H2} = 0, {xi3, H2} = 0, {xi4, H2} = xi2 (p.1152); this is what gives the normal form in the three combinations of (23).
  The change of variables to xi costs at most a factor 2^M in the norms (lemma A.4).
Bounds (proposition 5.1, eq. 41-42): ||G_s||_rho <= (2^s fact(s)/Omega_s^(s+1)) ||F_s||_rho, ||Z_s|| <= 2^(s/2) ||F_s||_rho, with
Omega_s := min over k in Z^2 minus {0}, |k2| <= s of min{|<k, omega>|, 1} (the minimum over k up to |k2| <= s only, not over |k|_1). Under the Diophantine
condition (lemma A.7) Omega_s >= gamma-tilde s^(-tau) with gamma-tilde = min{(3/2 + |omega2/omega1|)^(-tau) gamma, |omega1|, 1}. Remark 5.2: in the non-semi-simple
case these bounds, even with the Diophantine condition, do NOT give convergence of G = sum G_s (a strong difference from the semi-simple case).
Bounds for the process (proposition 6.1, lemma 6.2): ||G_s||_{3 rho0/4} <= (1/2)fact(s-1) beta_3...beta_s lambda^(s-3) c^(s-2) Delta_s^(s-3) / R0^(s-2) (read on the image of p.1158; Delta_s here is the homogenisation factor of eq. 50, Delta_j = 4 + (4/(e rho0)) sum_{l=3}^{j-1} 1/l, not the small divisor) with
beta_j = 2^j fact(j)/Omega_j^(j+1) (this is the printed beta_j; (50) in the text layer shows the exponent j+1 on Omega_j) and a universal constant lambda = 3 lambda-tilde, from recurrences (55)-(56) with a_3 = 1 and b_{l,0} = 1; lemma A.5: a_k <= lambda-tilde^(k-3), b_{l,m} <= lambda-tilde^m.
The key product (lemma A.8, eq. A3): beta_3...beta_s <= d_eps (s + 1)^((tau + 1 + eps)(s+1)^2/2), the source of the Lambert-function order, with remark A.9 saying numerics suggest the exponent cannot be lowered.
Proposition 7.1: with R-hat^(r) = ((sigma - 1)/sigma' ... ) R0 / (r-2)^((tau+1+eps')(r-2)/2) the remainder is <= c R0^2 (R/R-hat)^(r+1) for R <= R-hat^(r). Proposition 8.1 optimises r: minimising
h(r) = log(c R0^2) + (tau+1+eps)(r-2)(r+1) log(r-2)/2 + (r+1) log R gives, from its dominant part (tau + 1 + eps)(r-2) log(r-2) + log R = 0, the choice r_opt, with (r-2)^((tau+1+eps)(r-2)) <= 1/R, hence the remainder bound
R^(r_opt/2).

## 4. Every number in the paper that is usable as a test (READ page, tag)

1. lambda-tilde = 20.362 07... (remark A.6, p.1169). It is defined (lemma A.5, p.1169) as lambda-tilde := max{2^(3/2) C/(1 - D), A, B}/6 with
   A = (3/2) sum_{j>=1} fact(j+2)fact(j+3)/fact(2j+2), B = 2 sum_{j>=1} (fact(j+2))^2/fact(2j+1), C = max{(2/3)A, (3/4)B}, D = 9 - 6 e^(1/3).
   COMPUTED (mpmath, 30 digits, by direct summation of those series): A = 20.7484407, B = 32.2813582, C = 24.2110186, D = 0.6263254, and the printed
   formula gives lambda-tilde = 30.5431120. The printed numerical value 20.36207 is exactly 2/3 of that: 30.5431120/1.5 = 20.3620747. So the printed number and the printed definition
   disagree by a factor 3/2, to eight digits. A consistent reading is that the printed value corresponds to C = B/2 (that is, "(1/2)B" in place of "(3/4)B" in the printed
   definition of C), or to a factor 3/2 lost elsewhere in the recurrence; I cannot tell which from the paper. INFERRED: a slip in one of the two, harmless for the theorem (lambda-tilde is only used as "some constant greater than 1"). A test of a normal-form
   code must not pin lambda-tilde to either value without this caveat. The factor is the same to 1e-8, so it is not rounding.
2. D = 9 - 6 e^(1/3) = f(1/3) with f(x) = ((x-1)e^x + 1)/x^2 (p.1169): 0.62632545 (COMPUTED, and < 1 as the proof requires).
3. Figure A1 (p.1171): the sequence alpha_s = log(beta_3...beta_s)/(s^2 log s) for the golden-mean frequency vector omega = (1, (sqrt5 - 1)/2) (tau = 1), plotted for s from 5e4 to 1e6; the curve is a sawtooth between 0.9652 and 0.9738 on the axis (axis range 0.964 to 0.974), starting near 0.9655 at s = 5e4 and ending near 0.9738 at s = 1e6;
   the text concludes lim sup alpha_s = alpha = 1 (the limit is approached very slowly). COMPUTED (reproduction): with Omega_j = ||F_k omega2|| = omega2^k for the largest Fibonacci F_k <= j, beta_j = 2^j fact(j)/Omega_j^(j+1), I get
   alpha_s = 0.96552 (s = 5e4), 0.96842 (1e5), 0.96842 (2e5), 0.97196 (4e5), 0.97263 (6e5), 0.97183 (8e5), 0.97376 (1e6), matching the figure's endpoints (0.9655 and 0.9738) and its sawtooth. This confirms the printed definition of beta_j (Omega_j with exponent j+1) and of Omega_s.
   This is a reproduction of a plotted curve, not of printed digits; I did not digitise the curve beyond reading its endpoints and axis.
4. r_opt(R) (Theorem 1.1 iv; proposition 8.1). The relation (r-2)^((tau+1+eps)(r-2)) = 1/R holds exactly for the real r_opt (before taking the integer part): from W(z) e^(W(z)) = z with z = log(R^(-1/(tau+1+eps))). A direct test of any implementation. COMPUTED values (eps = 0.01):
   R = 1e-2: r_opt = 4.500 (tau = 1), 4.084 (tau = 2); R = 1e-3: 5.067, 4.502; R = 1e-4: 5.587, 4.887; R = 1e-6: 6.542, 5.591; R = 1e-8: 7.422, 6.238; R = 1e-12: 9.043, 7.427.
   Remainder bound R^(r_opt/2): for tau = 1, R = 1e-2: 3.2e-5; R = 1e-4: 6.7e-12; R = 1e-8: 2.1e-30; R = 1e-12: 5.6e-55. The normal form degree (floor(r_opt/2) in Zbar) is only 2 to 4 over R from 1e-2 to 1e-12: very little normalisation is justified by this theorem.
5. Orders: Zbar has degree at most floor(r_opt/2), starting at 2; the remainder Taylor expansion starts at degree r + 1 in the adapted degree (p.1163).
6. Not printed anywhere: R*, d-tilde, c, rho0, any value of omega1, omega2, tau, gamma for a real system; any multiplier; any energy; any orbit; any mu; any normal-form coefficient (Zbar is not given explicitly beyond its structure).

## 5. Reconciliation with the project

### 5.1 Where the transition appears in what the project holds

- Hadjidemetriou 1975b Table I (equal masses; `docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md`): rows 1-10 stable with printed b1, b2 real and |b| < 2; rows 11-15 have no b (Delta < 0). The change between rows 10 and 11 is a Delta = 0 event with |b| < 2: a collision of two unit-circle pairs, the transition this paper treats (a 1:-1 collision), provided the pair signatures are opposite (the paper's epsilon invariant; the Krein signature, which I did not compute). COMPUTED from my recomputed monodromy eigenvalues of that digest: row 10 has multipliers
  -0.8604 +- 0.5099 i and -0.5885 +- 0.8083 i (angles about 2.60 and 2.20 rad), row 11 has a complex quartet with moduli about 1.096 and 0.912 and both members at angle about 2.38 rad; so the collision is at an angle of about 2.4 rad (about 137 degrees), omega2/omega1 about 0.38 (not close to a low-order rational), b1 = b2 about 1.45 at the collision (alpha/2 between 1.449 at row 10 and 1.457 at row 11, from the computed alpha values 2.8977 and 2.9140). That is a sourced-in-the-project bracket (Delta changes sign between two printed rows, one of which has printed b values) but not a printed critical orbit.
- No project code computes a normal form. grep of `src/` for "normal form", "normal_form", "lie series", "giorgilli" finds only `genome/bcr4bp_continuation.py` and a docstring mention in `search/cr3bp_3d_family_tracer.py`, neither a normal-form implementation.
- Classifiers: `search/er3bp_floquet.floquet_classify` (unit-circle tolerance 1e-3 on |lambda|; per the `#931` entry in `data/OUTSTANDING.md` it labels a quartet unstable by modulus and cannot return "stable"), `search/cr3bp_3d_family_tracer._classify_floquet` (unit_tol = 1e-3, picks the trivial pair as the two eigenvalues nearest +1, in a 6 x 6 monodromy, and "unstable" only for |lambda| > 1 + unit_tol), `search/cr3bp_periodic.barden_stability` (one pair only), `search/er3bp_periodic.monodromy_eigenstructure` (saddle plus centre; see `#931`).
  The 1e-3 tolerance is the quantity that matters for a 1:-1 transition, and section 6 shows why.

## 6. Synthetic test built from the paper's own normal form (COMPUTED, not sourced)

To get a controlled case I took the quadratic part of the normal form, adding one symmetric unfolding term: H = omega2 (y1 x2 - y2 x1) + (y1^2 + y2^2)/2 + (mu/2)(x1^2 + x2^2), omega1 = 1, omega2 = (sqrt5 - 1)/2 (golden mean, so that omega1/omega2 is irrational and Diophantine with tau = 1), period T = 2 pi/omega1. The linear flow matrix A (state x1, x2, y1, y2: xdot = dH/dy, ydot = -dH/dx), the monodromy M = expm(T A). (This unfolding is my construction; the paper treats only mu = 0, the transition orbit itself.) Results:
- mu = 0: the Hamiltonian-matrix eigenvalues are +- i omega2 each doubled; M has the double multiplier lambda = exp(2 pi i omega2/omega1) = -0.73737 - 0.67549 i (and its conjugate), and rank(M - lambda I) = 3 for the 4 x 4 block: not diagonalizable, one Jordan block of size 2 per eigenvalue (the non-semi-simple case, confirmed). Condition number of the eigenvector matrix about 2e16.
- mu > 0 (with epsilon = +1 and this sign of the unfolding term): complex instability. mu = 1e-4: Hamiltonian eigenvalues +-0.01 +- 0.61803 i, |lambda| = 1.064848 and 0.939101 (each twice), Delta = alpha^2 - 4(beta - 2) = -2.886e-2.
- mu < 0: stable. mu = -1e-4: eigenvalues +- 0.62803 i and +- 0.60803 i, all |lambda| = 1, Delta = +2.878e-2, two distinct unit-circle pairs.
- Scalings (COMPUTED): Delta/mu = -288.22 for all small mu tested (1e-4, 1e-6, 1e-8), so Delta is linear in mu and changes sign at the transition; the distance of the multiplier from the unit circle is max|lambda| - 1 = T sqrt(mu) to four digits (mu = 1e-6: 6.303e-3 against T sqrt(mu) = 6.283e-3; mu = 1e-8: 6.285e-4 against 6.283e-4). So the multiplier moduli move like sqrt(mu) while Delta moves like mu.
- Sensitivity of a computed monodromy at the transition (COMPUTED): at mu = 0 the roundoff-only error of max ||lambda| - 1| is 2.9e-15, but adding random noise of size eta to the matrix entries (the effect of an integration or finite-difference error in the monodromy matrix) makes max ||lambda| - 1| about 1.3 sqrt(eta) (median of 200 trials: eta = 1e-14 gives 1.4e-7; 1e-12 gives 1.5e-6; 1e-10 gives 1.3e-5; 1e-8 gives 1.4e-4), the square-root sensitivity of a double eigenvalue with a Jordan block.
Consequences for classification near the transition:
- The sign of Delta is a far better indicator of the transition than a modulus test: Delta responds linearly in mu and its noise is of the order of eta (not sqrt(eta)). The `#931` entry's plan to report alpha, beta, Delta and decide by modulus is sound away from the transition; within it, report both.
- A modulus test with tolerance tol (the project's classifiers use 1e-3) labels the quartet as "stable/marginal" whenever T sqrt(mu) < tol, here for mu < (1e-3/(2 pi))^2 = 2.5e-8, that is for |Delta| below about 7e-6 (COMPUTED, using Delta = -288.2 mu); with a monodromy error eta = 1e-8 the noise floor 1.4e-4 is below the 1e-3 tolerance, so the tolerance, not the numerical noise, sets the blind window. The width of this window in the family parameter depends on how fast mu varies along the family and on T; it must be measured per family, not assumed.
- Which side is unstable is not decided by the double multiplier. Here epsilon = +1 and the sign of the unfolding give instability for mu > 0. The paper states epsilon is an invariant of the collision (and that t -> -t exchanges the two contexts) but does not discuss the unfolding. The side on which the quartet appears is determined by the Krein signature of the colliding pairs and the family's direction of motion; a classifier should report the sign of Delta along the family and not infer the side. (INFERRED from the printed epsilon and my synthetic case.)
- A sign-change detector across family members (alpha, beta, Delta) is the right detector, with a stored bracket [Delta_i, Delta_{i+1}] and a refinement by bisection in the family parameter, then a double-multiplier test at the bracket centre: rank(M - lambda I) = 3 numerically (singular values 4 of M - lambda I with one tiny and the others not), non-diagonalizability, angle theta with lambda = e^(i theta), theta not equal to 0 or pi.

## 7. Techniques applicable to the project's problems

Entries read from `data/OUTSTANDING.md`: `#931` (b) to (d), `#899`, `#916`.

### 7.1 `#931`: detecting and classifying the 1:-1 transition along a family

1. Detector. Along a continued family (for example the `#899` mass or Jacobi-constant continuation, or the `#931` (a) mass continuation) store alpha, beta, Delta (and b1, b2 when Delta >= 0) at every step. A transition of the kind treated here is a sign change of Delta with |b| < 2 at the collision (equivalently |alpha/2| < 2 at Delta = 0) and with the colliding pairs not at +-1. Discriminate it from (i) a period-doubling or tangent event (b crossing -2 or +2 with Delta > 0, one pair at +-1), and (ii) a collision at +-1 (a fold-like degeneracy). Test, from the printed data: Hadjidemetriou 1975b Table I rows 10 and 11 (Delta > 0 with b1 = 1.7211, b2 = 1.1766 at row 10; Delta < 0 at row 11, none printed) must be flagged as a Delta sign change with |alpha/2| about 1.45 < 2. A negative control from the same table: rows 1 to 4 have b1, b2 close to -2 with Delta > 0 (approach to a period-doubling-type |b| = 2 boundary, no Delta sign change) and must not be flagged. This needs the `#931` (a) integrator to reproduce the rows, as the `#931` entry already plans.
2. Tolerance near the transition. Use Delta's sign as the primary indicator and the modulus test only as a secondary one, because section 6 shows the modulus test is blind to a window of width |Delta| below about 288 (tol/T)^2 (for the synthetic case, 7e-6 at tol = 1e-3) and the eigenvalue noise at the double multiplier is 1.3 sqrt(eta). Concretely: for a monodromy computed to relative error eta, do not label |lambda| - 1 < 2 sqrt(eta) as "on the unit circle" evidence of stability, and do not treat |lambda| - 1 of that size as evidence of instability. For eta = 1e-12 the threshold is about 2e-6, far below the 1e-3 tolerance in `floquet_classify` and `_classify_floquet`, so the 1e-3 tolerance dominates and the blind window is set by it.
3. A sourced test case does not exist in this paper (section 0). Best available: (a) Hadjidemetriou 1975b Table I rows 10-11 (a Delta sign change between printed rows, collision angle about 2.4 rad computed); (b) the synthetic normal-form case of section 6 as a unit test of the classifier. I did not derive a closed form for its eigenvalues (only computed them: for the sign convention used, mu < 0 gives two distinct imaginary pairs and mu > 0 a quartet), so the unit test should compare against the eigenvalues of expm(T A(mu)) computed independently (for example by a different matrix-exponential routine), not against a formula. Controls: the exact double multiplier at mu = 0 with rank 3, Delta/mu = -288.2 at small mu, max|lambda| - 1 = T sqrt(mu) for mu > 0.
4. Non-diagonalizability as a flag. At the transition the monodromy 4 x 4 block has the Jordan structure above; a numerical test is the condition number of the eigenvector matrix (about 2e16 at mu = 0, 1e2 at |mu| = 1e-4 in the synthetic case) or the rank test. A classifier should output a "critical (1:-1)" tag when Delta is within noise of zero and |b| < 2 and not the stable tag.
5. The paper's Diophantine requirement. Away from the transition the resonance is irrational and generic; in a family the angle theta of the colliding pair is whatever the family gives; there are infinitely many other (rational) resonances at which the monodromy has multipliers that are roots of unity. If theta/2 pi is close to a low-order rational at the collision (the 1975b bracket, about 0.38, should be checked against low-order rationals before it is called generic), a further resonance interacts; this paper's theory excludes that.

### 7.2 `#899` and `#916`: what happens to families and tori born at a complex-instability transition

What this paper does and does not say (READ): it treats only the transition orbit. It does not treat the nearby periodic orbits of the family, nor the invariant objects that exist on either side; it cites for those the "quasi-periodic Hamiltonian Andronov-Hopf bifurcation" (Introduction, p.1142; reference [18]), the dynamics of the normal form (reference [16], in press), invariant curves near Hamiltonian-Hopf bifurcations of four-dimensional symplectic mappings (Jorba and Olle 2004, Nonlinearity 17:691, reference [10]) and the vertical L4 application (reference [17]); none of these is held. So the following is mapping and inference, with the limits stated.
- What can be taken from this paper for the stable side (INFERRED from the printed normal form and the cited program, not printed): on the stable side of the transition the normal form has two real normal frequencies and invariant tori are expected around the periodic orbit (as in the elliptic case); near the transition they are governed by the non-semi-simple normal form of the theorem, whose remainder is only R^(r_opt/2). So any stability or persistence claim for a stable member that is close to a Delta = 0 event rests on a weaker estimate than for a member away from it: with the numbers of section 4 item 4, r_opt is 4 to 9 for R from 1e-2 to 1e-12, so the integrable approximation holds to a power of R and is not exponentially small. For `#916`'s test of the printed persistence conjecture on stable members as invariant curves of the one-period map, this means: members whose Delta is small (close to the transition) are the worst test candidates; add Delta (or the angular separation of the two normal frequencies) as a column, and expect the invariant curves to be unreliable or fractured there. Control: the Hadjidemetriou 1975b orbit a (eq. 73, stable, iterated orbits form invariant curves with about 42.7 intersections per circuit, computed) has Delta = -5.5e-5 within noise (two nearly equal rotation angles 0.1466 and 0.1474), that is, it is itself very close to a 1:-1 collision of two unit-circle pairs (computed in the 1975b digest, section 5(d)): a published case where invariant curves were still found to be closed or nearly closed (Fig. 3 and 4 of that paper) for perturbations of 0.02 to 0.1, so proximity to the collision did not destroy them over 40 to 1200 iterations. That is evidence against over-interpreting the weak remainder bound, not a contradiction: the theorem is an upper bound.
- The unstable side: a complex-quartet orbit has a two-dimensional unstable manifold and no invariant tori from that orbit alone; the paper says nothing about the unstable side. For `#899` (continuation through near-collision seeds and families): when a continued orbit passes into Delta < 0, do not read the loss of "stable" as the end of the family or as a fold. The Introduction (p.1142) states that a periodic orbit whose monodromy has no non-trivial eigenvalue equal to 1 lies in a one-parameter family parametrised by the energy; at the Hopf collision the double multiplier lambda is not 1, so the family continues through the transition and a fixed-Jacobi-constant or fixed-period corrector keeps a nonsingular Jacobian. A positive control for the `#899` and `#931` (d) logging: det(M - I) of the normal block is |lambda - 1|^4 > 0 at the collision (COMPUTED for the synthetic case, golden-mean omega2: det(M - I) = 12.07 = (2 - 2 cos(2 pi omega2))^2), unlike at a fold or at a period-doubling, where an eigenvalue +1 or -1 appears. A sourced test: Newton on the periodic orbit converges across Hadjidemetriou 1975b Table I rows 10 to 11 with the same corrector, no loss of rank (both sides were computed by the author).
- Where the families born at the transition would be treated (the quasi-periodic Hopf bifurcation: invariant tori or invariant curves of the Poincare map on one side): not in this paper; acquire reference [17] (the L4 vertical family, a restricted three-body application of exactly this Hopf transition, with periodic orbits and invariant tori computed) and reference [10] before using the theory for `#916`.

### 7.3 Controls and what not to do

- Do not use this paper to justify a stability radius: it gives no R*, and its remainder bound is polynomial in R.
- Do not pin lambda-tilde or any proof constant in a test (section 4 item 1).
- A positive control for any numerical normal-form code (should the project ever build one): the structure, not the numbers: Zbar depends only on I, x1^2 + x2^2 and y1 x2 - y2 x1; odd-degree terms vanish (Z_s = 0 for odd s, proposition 4.1); non-removable terms exist only for k = 0 and |m| = |n|; and the closed-form homological solutions of eq. (37) and (38) are algebraic enough to test against finite-order manual examples.

## 8. Recommended follow-ups (no task numbers registered)

1. Acquire Olle, Pacha and Villanueva 2004 (Cel. Mech. Dyn. Astron. 90:89, the vertical family of L4: a restricted three-body Hopf transition with a computed critical orbit, multipliers and tori) and Jorba and Olle 2004 (Nonlinearity 17:691); either would give the sourced critical orbit and the unfolding that this paper does not.
2. For the `#931` classifier: add the synthetic test of section 6 (compare against the numerically computed eigenvalues, Delta linear in mu, max|lambda| - 1 about T sqrt(mu), rank 3 at mu = 0), a Delta sign-change detector with bisection, and a "critical (1:-1)" tag; keep a modulus test only as a secondary indicator and report its blind window.
3. Run the detector over Hadjidemetriou 1975b Table I (rows 10-11) and any CR3BP 3D families the project holds, to get the first sourced Delta bracket; compute the Krein signature of the colliding pairs there before calling it a Hamiltonian Hopf (a collision of equal-signature pairs is not a bifurcation).
4. Treat the printed lambda-tilde (20.36207) against the printed definition (30.54311) as an erratum candidate to report if the paper's constants are ever relied on.

## 9. Summary

- Paper: pure analysis of the normal form around a non-semi-simple 1:-1 resonant periodic orbit of a 3-degree-of-freedom Hamiltonian (the Hamiltonian-Hopf transition orbit, where a stable orbit's four multipliers collide pairwise on the unit circle and leave it as a quartet). Normal form (Theorem 1.1): Z = omega1 I + omega2 (y1 x2 - y2 x1) + (y1^2 + y2^2)/2 + Zbar(I, -(x1^2 + x2^2)/2, (y1 x2 - y2 x1)/2), remainder at most R^(r_opt/2) with r_opt = 2 + exp(W(log(1/R^(1/(tau+1+eps)))), polynomially not exponentially small because the homological equations have growing Jordan blocks with small divisors (factor fact(M+N)/Delta^(M+N+1)).
- No example, no table, no critical orbit, no multipliers, no mu: the restricted problem enters only through references. Sourced numbers: lambda-tilde = 20.36207 (does not match its own printed definition, which gives 30.54311; factor exactly 3/2) and Figure A1 (alpha_s from 0.9655 to 0.9738 for the golden mean, reproduced).
- For `#931`: Delta = 0 with |b| < 2 is this transition; Delta changes linearly across it while multipliers move like sqrt(mu), monodromy noise is amplified to 1.3 sqrt(eta) at the double multiplier, so decide by Delta's sign and treat modulus tolerances as a blind window; Hadjidemetriou 1975b Table I rows 10-11 bracket one such event (collision angle about 2.4 rad, computed).
- For `#899` and `#916`: the paper is silent about nearby orbits and tori; the family continues through the transition, tori on the stable side are less reliable near it, and the sourced treatment is in two cited papers not held.
