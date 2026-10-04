# Digest: Rhouma & Chicone 2000, "On the continuation of periodic orbits"

Citation: M. B. H. Rhouma and C. Chicone, "On the continuation of periodic orbits", *Methods and
Applications of Analysis* 7(1):85-104 (March 2000), International Press, DOI
`10.4310/maa.2000.v7.n1.a5`. Revised 16 July 1999. 20 PDF pages (journal pp. 85-104). Filed in the
private paper corpus as
`rhouma-chicone-2000-continuation-periodic-orbits-forced-oscillators-maa-7-85-doi-10.4310-maa.2000.v7.n1.a5.pdf`
(md5 `c12838a9f63a8f308a0ce80a981ca057`). The PDF has a text layer, but its equations are badly
garbled; every equation below was read from page images (130 dpi), and the numbers in the
examples were re-checked (section 5).

Wanted for: `#916` and `#922` (persistence of periodic orbits and tori when a periodic forcing such as
the Sun is switched on), and as the exact statement behind the Melnikov step used in `#884` and `#905`.

How to read the markers: PRINTED means taken from the paper; DERIVED means worked out here from the
paper's statement and not printed in it; COMPUTED means a number I computed with a script this
session (scratch only, nothing committed); PROJECT means a statement about this repository.

## 1. What the paper does

Setting (PRINTED, Eq. 1.1 to 1.3, p. 85). The system is

    y' = f(y) + eps g(t, y, eps),     y in R^n,  eps >= 0,

with g periodic in t of period T > 0 for every fixed (y, eps). The unperturbed system y' = f(y) has a
*resonant period manifold* Z: a manifold (dimension m) consisting entirely of T0-periodic orbits,
with relatively prime positive integers M and N such that

    M T0 = N T.                                                                   (1.3)

The question is persistence: for small eps > 0, is there a periodic solution of period N T near each
unperturbed orbit in Z. The method is Poincare's: a zero of the displacement function
v -> y(NT, v, eps) - v is an NT-periodic solution. The unperturbed solution through zeta in Z is
NT-periodic because NT = M T0 (the orbit is traversed M times).

Persistence is defined (p. 87) in a weak sense: zeta in Z persists if there is eps0 > 0 and a
CONTINUOUS function beta: [0, eps0) -> R^n with beta(0) = zeta and y(NT, beta(eps), eps) - beta(eps) = 0.
The paper deliberately allows non-smooth continuations beta(eps) = zeta + eps^(1/p) gamma(eps^(1/p))
(gamma smooth) so that the implicit function theorem can be applied to a rescaled problem
(eps = mu^p).

Notation (PRINTED, Section 2).

- Phi(t, zeta) = y_v(t, zeta, 0): principal fundamental matrix at t = 0 of the first variational
  equation w' = Df(y(t, zeta, 0)) w along the unperturbed solution.
- Infinitesimal displacement: Phi(NT, zeta) - I.
- K(zeta) = kernel of Phi(NT, zeta) - I. The tangent space T_zeta Z is always contained in K(zeta)
  (differentiate y(T0, sigma(s), 0) - sigma(s) = 0 along a curve in Z).
- Z is *normally nondegenerate* if dim K(zeta) = dim T_zeta Z = m for every zeta in Z. Z is
  *l-degenerate* if dim K(zeta) - dim T_zeta Z has the constant value l on Z.
- y_eps(NT, zeta, 0) = Phi(NT, zeta) integral_0^NT Phi^-1(t, zeta) g(t, y(t, zeta, 0), 0) dt      (2.8)
- y_vv(NT, zeta, 0) = Phi(NT, zeta) integral_0^NT Phi^-1(t, zeta) D^2 f(y(t, zeta, 0)) (Phi(t, zeta), Phi(t, zeta)) dt   (2.7)
  (I read Eq. 2.7 as the variation-of-constants formula for the second variation; the garbled print
  shows the integrand D^2 f applied to the pair of columns of Phi. DERIVED reading, consistent with the
  line above Eq. 2.7.)
- Proposition 2.1 (smooth bases). If M(z) is a smooth family of n x n matrices of constant rank r, there
  are smooth families K(z) (columns: a basis of a complement of the kernel, then a basis of the kernel)
  and R(z) (columns: a basis of the range, then a basis of a complement of the range) with
  R^-1(z) M(z) K(z) = diag(I_r, 0). Complements may be orthocomplements. (Used to put the displacement
  equation in normal form; the constant-rank hypothesis is Remark 4.4.)

## 2. The theorems, stated exactly

pi_R(zeta) and pi_C(zeta) denote the projections onto the range of Phi(NT, zeta) - I and onto a complement
of that range, corresponding to R^-1(zeta) (first n-k components, last k components, where k = dim K).

**Theorem 3.1 (PRINTED, p. 92). Normally nondegenerate case.** Suppose the unperturbed system has a
normally nondegenerate resonant period manifold Z, Phi(t, zeta) is the principal fundamental matrix of the
first variational equation at zeta in Z, and pi_C(zeta) is the projection onto the complement of the range
of Phi(NT, zeta) - I. If the *subharmonic Melnikov function*

    zeta -> pi_C(zeta) y_eps(NT, zeta, 0)                                          (3.3)

has a simple zero zeta0 in Z, then zeta0 persists. Conclusion in the Introduction (p. 86): there is a
smooth curve eps -> gamma(eps) on an open set N0 containing 0 such that beta(eps) = zeta0 + eps gamma(eps)
is the initial point of a periodic solution of the forced system for each eps in N0 (a smooth,
two-sided continuation).

The bifurcation function (their F) has the two components 0 = nu + pi_R y_eps + O(mu) and
0 = pi_C y_eps + O(mu) (Eq. 3.2) in unknowns (nu, zeta, mu). The leading displacement at the solution
is nu0 = -pi_R(zeta0) y_eps(NT, zeta0, 0) (Eq. 3.4).

Nondegeneracy of the zero, restated (PRINTED, p. 93). zeta0 is a simple zero of (3.3) if and only if

1. y_eps(NT, zeta0, 0) lies in the range of Phi(NT, zeta0) - I (this is "zeta0 is a zero"), and
2. the composition pi_C(zeta0) y_{v eps}(NT, zeta0, 0) has rank m (the paper words it as the columns of
   y_{v eps} generating a basis of a complement of the range; the coordinate argument on the same page
   proves the rank-m form: the derivative of (3.3) along Z is pi_C y_{v eps} DT, with T a coordinate map).

If condition 1 fails at every point of Z, no point of Z has a SMOOTH continuation (p. 93); non-smooth
continuations (eps^(1/p), p > 1) remain possible and are not excluded.
Two degeneracies the paper names and leaves open: Z normally degenerate (treated next), and a zero of
(3.3) that is not simple ("analysis near a singularity of the Melnikov function is beyond the scope").

**Theorem 4.1 (PRINTED, p. 94). Normally degenerate case.** Let Z be an l-degenerate resonant period
manifold, k = m + l. Choose a smooth basis K_1..K_{n-k}, K_{n-k+1}..K_{n-m}, K_{n-m+1}..K_n of R^n
adapted to the constant-rank hypothesis: the first n-k vectors complement the kernel, the last k span
the kernel, and the last m of those are tangent to Z; the l vectors K_i that complement T Z inside the
kernel carry the new variable kappa in R^l. The *generalized subharmonic Melnikov function*
M: Z x R^l -> R^k is

    M(zeta, kappa) = pi_C(zeta) ( y_eps(NT, zeta, 0) + (1/2) y_vv(NT, zeta, 0) (K(zeta) V(0, kappa))^2 ),     (4.2)

with K(zeta)V(0, kappa) = sum_{i=1..l} kappa_i K_i(zeta)      (4.3). If (zeta0, kappa0) in Z x R^l is a simple
zero of M, then zeta0 persists, and there are TWO continuation curves, both defined for eps >= 0
(they come from mu = +sqrt(eps) and mu = -sqrt(eps) in the rescaled equation, eps = mu^2). For l = 0 this
reduces to (3.3). Only the first term of M depends on the perturbation g.

The second-order term appears because with eps = mu^2 and xi = mu K V(nu, kappa), the displacement
expansion (2.2) divided by mu^2 becomes (4.1): 0 = (nu, 0) + mu R^-1 (y_eps + (1/2) y_vv (K V)^2) + O(mu^2).

**Remark 4.2 (PRINTED).** Negative eps is handled by the substitution eps -> -eps with g(t, y, eps) -> -g(t, y, -eps),
which gives y_eps -> -y_eps and the function M_-(zeta, kappa) = pi_C(-y_eps + (1/2) y_vv (K V)^2). A
simple zero of M continues to eps >= 0, a simple zero of M_- to eps <= 0, and "it is very likely that
the simple zeros of M and M_- are different" (the Section 7 example shows this).

**Remark 4.3 (PRINTED).** Phi(NT, zeta) must be "known". For planar autonomous systems Diliberto's theorem
gives the variational equation by quadratures along the orbit; "we do not know how to obtain a substitute
for Diliberto's theorem in R^n for n > 2". (Numerically this is no obstacle: Phi is integrated.)

**Remark 4.4 (PRINTED).** The constant-dimension assumption on K(zeta) is needed for the smooth bases. Example
(Eq. 4.4): x' = y - x z (x^2 + y^2 - 1), y' = -x - y z (x^2 + y^2 - 1), z' = 0. The cylinder
x^2 + y^2 = 1 is a 2 pi-period manifold; dim K = 2 for z > 0 (normally nondegenerate, m = 2) and dim K = 3 at
z = 0 (2-degenerate, m = 1), so three separate period manifolds (z > 0, z = 0, z < 0) are analysed.

**Section 5 (PRINTED): Cronin's theorem.** Cronin's 1996 continuation result (Meth. Appl. Anal. 3:370)
for y' = f(y) + mu^2 g(t, mu) with eigenvalue 1 of Phi(T, zeta) of algebraic multiplicity one yields a
bifurcation equation A y^2 + B = 0, but Proposition 5.1 shows A = 0 identically, so Cronin's theorem is
logically correct and vacuous; it reduces to B = 0, in agreement with Theorem 3.1. Proposition 5.1: for
zeta in Z and u in the kernel, with u = u_nor + u_tan (u_tan tangent to Z),
pi_C y_vv(NT, zeta)(u_tan, u_tan) = 0 and pi_C y_vv(NT, zeta)(u_tan, u_nor) = 0. Consequence for us:
second-variation terms along the orbit direction never contribute to the bifurcation function.

**Section 6 (PRINTED): other exponents and the linear case.** With eps = mu^p: for 1 < p < 2 the equation
reduces to the p = 1 continuations; for p > 2 the leading function (zeta, kappa) -> pi_C (y_vv (K V kappa)^2) is
homogeneous in kappa, so it cannot have a simple zero. For an unperturbed LINEAR system
y' = A y + eps g, all derivatives y_vv and higher vanish and, when Phi(NT) - I = 0 (the planar oscillator
at resonance), it suffices that zeta -> -y_eps(NT, zeta, 0) has a simple zero (the classical result).
Remark 6.1: the "controllably periodic" perturbations of Farkas (period of g allowed to vary with a
parameter) require eigenvalue 1 of Phi(T, zeta) to be simple; this paper extends that to eigenvalue 1
that is not simple.

## 3. Examples and every printed number

**Example 1 (Eq. 7.1, p. 100-101), one orbit, 1-degenerate.**

    x' = y - x (x^2 + y^2 - 1)^2 - eps cos t
    y' = -x - y (x^2 + y^2 - 1)^2 + eps sin t
    z' = z (x^2 + y^2) + eps sin t.

Period manifold Z = {x^2 + y^2 = 1, z = 0}, a single 2 pi-periodic orbit (T0 = 2 pi, forcing period
T = 2 pi, M = N = 1). Orbit coordinate theta: t -> (cos(t - theta), -sin(t - theta), 0). PRINTED:

- Phi(t) = [[cos t, sin t, 0], [-sin t, cos t, 0], [0, 0, e^t]], independent of theta.
- Phi(2 pi, theta) - I = diag(0, 0, e^(2 pi) - 1); kernel = the plane z = 0 (dim 2 against m = 1, so l = 1).
- K(theta) = [[0, cos theta, -sin theta], [0, sin theta, cos theta], [1, 0, 0]],
  R(theta) = [[0, 0, 1], [0, 1, 0], [e^(2 pi) - 1, 0, 0]].
- y_eps(2 pi, v, 0) = (-2 pi, 0, e^(2 pi) integral_0^(2 pi) e^(-t) sin t dt); after pi_C its contribution is (0, -2 pi).
  pi_C maps (u, v, w) to (v, u).
- (1/2) y_vv term projected: (-8 pi kappa^2 sin theta, -8 pi kappa^2 cos theta).
- M(theta, kappa) = (-8 pi kappa^2 sin theta, -2 pi (1 + 4 kappa^2 cos theta));
  M_-(theta, kappa) = (-8 pi kappa^2 sin theta, -2 pi (1 - 4 kappa^2 cos theta)).
- Simple zeros: (theta, kappa) = (pi, +-1/2) for M (continuations with eps >= 0) and (0, +-1/2) for M_-
  (eps <= 0). "For each simple zero there are two continuation curves. However for this simple case, these
  continuation curves likely match to form one smooth curve" (stated as "likely", not proved).

COMPUTED check (scratch script, DOP853 at rtol 1e-13, finite differences): at theta = 0.7, kappa = 0.3,
the projected (1/2) y_vv term is (-1.45725, -1.73011) against the printed formula (-1.45719, -1.73003)
(agreement 6e-5, the size of the finite-difference error), and y_eps(2 pi, zeta, 0) = (-6.28319, -5e-9, 267.2458)
against (-2 pi, 0, (e^(2 pi) - 1)/2 = 267.2458). So the Example 1 numbers are reproduced and are usable
as a test of any generic Melnikov-function code (see section 6).

**Example 2 (Eq. 7.2 to 7.4, p. 102-103), m uncoupled-limit-cycle oscillators, m-degenerate.** The j-th
oscillator is x_j' = mu_j y_j - x_j (x_j^2 + y_j^2 - lambda_j^2)^2 + eps g_{1,j}, y_j' = -mu_j x_j - y_j (...)^2 + eps g_{2,j}.
Z is an m-torus (product of limit cycles x_j = lambda_j cos(mu_j t - theta_j), y_j = -lambda_j sin(mu_j t - theta_j)),
l = m, pi_C = identity, kernel = all of R^(2m). PRINTED:

    M_+-(Theta, kappa) = -4 N T (lambda_1^3 kappa_1^2 cos theta_1, lambda_1^3 kappa_1^2 sin theta_1, ..., lambda_m^3 kappa_m^2 cos theta_m, lambda_m^3 kappa_m^2 sin theta_m)^T +- (p_{1,1}, p_{1,2}, ..., p_{m,1}, p_{m,2})^T,

    (p_{j,1}, p_{j,2})^T = integral_0^NT [[cos mu_j s, -sin mu_j s], [sin mu_j s, cos mu_j s]] (g_{j,1}(s, X0(s), Y0(s), 0), g_{j,2}(s, X0(s), Y0(s), 0))^T ds.

(The printed (7.2) shows lambda_j^3 in the first group and the m-th line as lambda_m^3; the text print of the
exponent on the second line of each pair is garbled in places, so I take lambda_j^3 kappa_j^2 for every row.)
Up to 2^m continuation curves at each continuation point.

Specific case m = 2, mu_j = lambda_j = 1, g_{1,1} = 0, g_{1,2} = x_2, g_{2,1} = 0, g_{2,2} = -x_1 + a cos t (a constant):
M(Theta, kappa) for eps >= 0 (Eq. 7.4):

    ( -8 pi kappa_1^2 cos theta_1 - pi sin theta_2,
      -8 pi kappa_1^2 sin theta_1 + pi cos theta_2,
      -8 pi kappa_2^2 cos theta_2 + pi sin theta_1,
      -8 pi kappa_2^2 sin theta_2 - pi cos theta_1 + pi a ).

For 0 < |a| < 1, eight simple zeros (two continuation points on the torus). PRINTED:

    (Theta, kappa) = (0, -pi/2, +-(1/4) sqrt 2, +-(1/8) sqrt(1 - a))   and   (pi, pi/2, +-(1/4) sqrt 2, +-(1/8) sqrt(1 + a)).

M_- (eps <= 0) has eight other zeros: (0, pi/2, +-(1/4) sqrt 2, +-(1/8) sqrt(1 - a)) and
(pi, -pi/2, +-(1/4) sqrt 2, +-(1/8) sqrt(1 + a)). At the point (theta_1, theta_2) = (0, -pi/2) continuation
curves exist along the four directions (+-(1/4) sqrt 2, 0, 0, +-(1/8) sqrt(1 - a)) in R^4.

PRINTED-NUMBER DISCREPANCY (COMPUTED, flag): substituting the printed Theta = (0, -pi/2) into the printed
(7.4) gives kappa_1^2 = 1/8 (so kappa_1 = +-(1/4) sqrt 2, matching) and, from the fourth component,
8 pi kappa_2^2 - pi + pi a = 0, i.e. kappa_2 = +-sqrt((1 - a)/8) = +-(1/(2 sqrt 2)) sqrt(1 - a), NOT the printed
+-(1/8) sqrt(1 - a) (which squares to (1 - a)/64). The same factor affects all eight zeros of each function and
the quoted continuation directions. I recomputed (7.4) from (7.2) and (7.3) (p_{1,1} = -pi sin theta_2,
p_{1,2} = pi cos theta_2, p_{2,1} = pi sin theta_1, p_{2,2} = -pi cos theta_1 + pi a) and it matches the
printed (7.4), so the slip is in the quoted zeros, not in the function. Treat the zero locations in kappa_2 as
a typesetting slip in the paper (the theta values, kappa_1 and (7.4) itself are consistent), and test against
the function (7.4), not the quoted kappa_2.

## 4. What the theorem says about the Earth-Moon problem (`#884`, `#902`, `#905`, `#916`, `#922`)

The setting maps as follows. The Earth-Moon synodic frame with the Sun as a time-periodic
fourth body (`core/bcr4bp.py`) is y' = f(y) + eps g(t, y): f is the CR3BP field (`core/cr3bp.py`,
autonomous, state in R^4 planar or R^6 spatial), and g is the Sun's direct plus indirect acceleration
with the Sun phase advancing at the synodic rate. Scaling the Sun mass, mu_S -> eps mu_S, makes g(t, y, eps)
exactly linear in eps, with g(t, y, 0) the physical solar acceleration, and eps = 1 the physical model.
The Sun's synodic period is T = 2 pi / omega_S = 6.79 TU (about 29.5 d). A CR3BP periodic orbit of synodic
period T0 satisfies the resonance condition if T0 = (N/M) T, the persisted orbit having period N T = M T0.
(PROJECT: this is what `sun_commensurate_period` and the `#884` stored orbits do, with their "n" = N and
their "a" = M.)

**DERIVED: the bifurcation function of the CR3BP is the Jacobi-constant work integral.** For an
autonomous Hamiltonian f the eigenvalue 1 of Phi(T0) has algebraic multiplicity at least 2 (the flow
direction f and the family direction), and gradient-of-C is a LEFT eigenvector of Phi with eigenvalue 1
(conservation of C gives grad C(y(t))^T Phi(t) = grad C(zeta)^T). So it spans the complement of the
range of Phi(NT) - I when the kernel is one dimensional, and pi_C y_eps = grad C(zeta) . y_eps. With (2.8)
and grad C(zeta)^T Phi^-1(t) = grad C(y(t))^T this gives

    Mel(zeta) = grad C(zeta) . y_eps(NT, zeta, 0) = integral_0^NT grad C(y(t)) . g(t, y(t), 0) dt = -2 integral_0^NT v(t) . a_Sun(t, r(t)) dt,

the first-order change of the Jacobi constant over NT per unit eps (the work of the solar force along the
unperturbed orbit). This is the same function that the `#884` module
`src/cyclerfinder/search/sun_forced_periodic_884.py` already computes ("Mel(theta0)", docstring lines
42-58; `melnikov_eval`, `melnikov_variational`, `find_zeros`) and that Brown et al. (2025) call their
Eq. 2.7; it was derived there from the same Lyapunov-Schmidt argument. So Theorem 3.1 is the statement
that makes that code's use of the function a theorem (sufficiency), and the paper adds what that code does
not state: the exact nondegeneracy hypotheses below, the degenerate case, and the M_- companion.

Because Z is the one orbit (the unperturbed period T0 differs along a family, so Z at fixed T0 is a single
orbit, m = 1), "the point zeta" is the phase along the orbit, equivalent to the phase theta0 of the Sun at
t = 0. The Melnikov function is therefore a function on a circle. DERIVED: its mean over the circle is
zero (the k = 0 term of the solar potential does no net work around a closed orbit, and every k != 0 harmonic
averages out), so unless it vanishes identically it has at least two sign changes. Combined with Theorem 3.1:
a nondegenerate commensurate orbit whose Melnikov function is not identically zero persists, at each simple
zero, as a forced periodic orbit for small eps, in at least two phase-distinct copies. (The standard
Poincare-Birkhoff picture, with the two copies of opposite stability type, is not in this paper; it is the
usual picture for a twist map and is stated here as background, not as a result of the paper.)

**Hypotheses that must be verified for the persistence claim (DERIVED consequences of Definitions in section 2).**

1. Resonance: T0 = (N/M) T_S exactly. The theorem says nothing about T0 not commensurate with T_S.
   DERIVED necessary condition: if a forced family of periodic orbits of bounded period N' T_S converges as
   eps -> 0, its limit is a CR3BP orbit whose period divides N' T_S, so a non-commensurate orbit cannot persist
   as a periodic orbit of bounded order. This is consistent with Rosales 2021 (generic-period orbits become
   2-tori, commensurate ones stay periodic), but the TORI themselves are not produced by this paper's method
   (they need KAM or normal-hyperbolicity arguments; section 5).
2. Normal nondegeneracy of the M-fold repeated orbit: dim ker(Phi(T0)^M - I) = 1 (planar and spatial alike, the
   orbit being the only point of Z). In the CR3BP this holds exactly when (a) the period is not stationary
   along the family at this orbit, dT0/dC != 0 (otherwise the Jordan block of the eigenvalue 1 splits and the
   kernel grows to 2, a fold of the period function: the paper's n = 2 remark on p. 86 that a periodic orbit
   in an annulus is nondegenerate iff the derivative of the period function does not vanish), and (b) none of the
   other Floquet multipliers lambda of the one-lap monodromy matrix satisfies lambda^M = 1 (a period-M-multiplying
   or tangent bifurcation of the repeated orbit). (a) and (b) fail at exactly the bifurcation points that
   `#905`'s "bifurcation detection in the family walk" must catch; at them Theorem 4.1 (l = 1 or more) applies
   instead and the continuation is in sqrt(eps).
3. A simple zero of Mel as a function of the Sun phase: Mel(theta0) = 0 and dMel/dtheta0 != 0 (the rank-1
   form of condition 2 of Theorem 3.1). Zeros where the slope also vanishes (the orbit is symmetric and the
   forcing symmetric: the `#884` review found zeros at the symmetric phases) are fine as long as the slope is
   nonzero; double zeros need the higher-order Melnikov functions of Brown et al. 2025 (`2026-10-04-digest-brown-peterson-henry-scheeres-2025-jas-high-order-resonance.md`),
   which the paper does not provide.
4. Mel not identically zero. If Mel is identically zero (the harmonic content of the orbit does not meet the
   forcing: in the project's model the quadrupole solar term has Sun-angle harmonics 0 and 2, the octupole 1 and 3,
   so only harmonics that are multiples of M contribute), the theorem is silent; this is Brown et al.'s "Melnikov
   function identically zero for h2" situation.
5. eps = 1 is not small. The theorem is local in eps with no radius (the paper cites Farkas 1978 [15] for
   existence regions and gives none). The physical solar coupling is mu_S / a_S^3 = 328900.5423 / 388.8111^3 =
   5.60e-3 per unit squared distance (COMPUTED from the constants in `core/bcr4bp.py`), small but not
   infinitesimal against the CR3BP's near-Moon sensitivities. Existence at eps = 1 must be established by
   continuing eps from 0 to 1 and reporting folds, not by invoking the theorem.

**Where the paper's setting does NOT match the CR3BP setting (stated plainly).**

- The paper's examples are dissipative (limit cycles, Phi has eigenvalue 1 only through the flow direction,
  or a z-direction that is unrelated to the flow). The CR3BP is conservative: eigenvalue 1 is always at least
  double, the nondegeneracy condition becomes the twist condition dT0/dC != 0 (item 2 above), and the cokernel
  is spanned by grad C, a structure the paper does not use and does not mention.
- The paper has no example of a Hamiltonian system, of n = 4 or n = 6, or of a family with a first integral.
  Remark 4.3's closed-form route (Diliberto) is available only for n = 2; for the CR3BP Phi is integrated.
- Theorems 3.1 and 4.1 need a single forcing period. The Sun alone provides one (T_S). Adding Oberon's or
  Titania's eccentricity (`#922`) gives two incommensurate frequencies: in the Titania frame the Titania-Oberon
  synodic period is 24.637 d and Oberon's anomalistic period 13.463 d, a ratio 1.830 (COMPUTED from the
  registry's printed periods 8.705869 d and 13.463237 d), not a ratio of small integers, so the forcing is
  quasi-periodic and the theorem does not apply. It applies to the existing planar circular two-moon model
  (period 5 synodic periods, `#890`) as a persistence test for forcing terms of that one period.
- Existence only. The theorem gives no stability, no count beyond the zero count, no size of the resonance
  zone, and no information on orbits in a resonance overlap region; nothing about the "near-commensurate"
  invariant-curve statement of `#916` (Ross & Roberts-Tsoukkas) beyond the exactly commensurate case.
- Smoothness: g must be smooth along the orbit. Close approaches to Earth or Moon are not a problem for g (the
  Sun term is smooth), but Phi(NT) can have multipliers up to 2.5e14 (the `#884` review table), so the rank
  decisions of the theorem (kernel dimension) cannot be made from singular values of Phi - I numerically. Use
  grad C as the known left null vector and test the right kernel through the Jordan structure at the family level.
- The Sun in `core/bcr4bp.py` is planar (z = 0), the orbit may be spatial. A planar orbit viewed in the spatial
  model has the vertical monodromy block's multipliers in addition; condition 2(b) must be checked for them too.

## 5. How to compute the bifurcation function numerically from the project's models

For an orbit given as a CR3BP periodic orbit (`core/cr3bp.py`, `cr3bp.propagate(system, state, T0, with_stm=True)`;
state vector 6, planar orbits have z = vz = 0) with T0 = (N/M) T_S:

1. Nondegeneracy screen on the unperturbed orbit. Monodromy Phi1 = Phi(T0). Eigenvalues of Phi1^M: exactly two near 1 (the
   Jordan pair), the others satisfying |lambda^M - 1| bounded away from 0 (report the minimum). Estimate dT0/dC from the
   neighbouring family members (`search/cr3bp_continuation.py`, `search/cr3bp_jacobi_arclength.py`); refuse or flag if it is
   within the continuation error of zero.
2. Bifurcation function by quadrature (cheapest; the form already in `sun_forced_periodic_884.melnikov_scan`):
   Mel(theta0) = -2 integral_0^NT v(t) . a_Sun(t; r(t), theta0) dt, with a_Sun from `core/bcr4bp.py`
   (`_sun_acceleration(x, y, z, t, system)` with `BCR4BPSystem.theta_sun0 = theta0`), integrating one lap and extending
   periodically (unstable orbits must not be integrated for M laps in one shot).
3. Cross-check by variational equation: augment the state with w' = Df(y) w + g(t, y(t), 0), w(0) = 0 (g = the Sun
   acceleration at the physical mu_S), then Mel = grad C(y(NT)) . w(NT). `melnikov_variational` does this.
   Cross-check by finite differences: y_eps = [y(NT; +eps) - y(NT; -eps)] / (2 eps) with
   `BCR4BPSystem(mu_sun = eps mu_S, ...)` through `propagate_bcr4bp(..., t0=0)`; also
   [C(y(NT; eps)) - C(zeta)] / eps -> Mel as eps -> 0. (Negative mu_sun is mathematically fine here: it is Remark 4.2's M_-.)
4. Scan theta0 over its period (2 pi / M, per the `#884` module), find sign changes, and require a nonzero slope at each
   zero (condition 3). Mirror-symmetric orbits have zeros at the symmetric phases by reversibility, so a zero there is
   not evidence of anything (the `#884` review's point about the positive control).
5. Predictor for the forced orbit, straight from Eq. 3.1/3.4: solve (Phi(NT) - I) xi = -y_eps(NT, zeta0, 0) in the
   least-squares sense with the component along the flow direction left free; set v(eps) ~ zeta0 + eps xi.
6. Corrector: Newton on y(NT; v, eps) - v = 0 using `propagate_bcr4bp(..., with_stm=True)`; at eps = 0 the Jacobian
   Phi - I is singular, so continue in eps with a bordered system or pseudo-arclength in (v, eps), starting from
   (zeta0 + eps xi, eps) at small eps. Report any fold in eps before reaching eps = 1.
7. For l >= 1 (orbit at a fold of the period function or a multiplier-1 point of the repeated orbit): the generalized
   function (4.2) needs the second sensitivity y_vv (second variational equations or a finite difference of the STM) and the
   kernel basis K_i; unknowns (zeta, kappa), k equations; the continuation starts at eps = mu^2 with the two signs of mu.
   Not implemented anywhere in the repository.

The same recipe applies unchanged to `core/ccr4bp.py` (one extra moon on a concentric circle, forcing period
2 pi / |omega_gan|): g is `_ganymede_acceleration`.

## 6. Techniques applicable to the project's problems

- `#905` (rerun `#884` in the corrected models). Add the two nondegeneracy screens of section 5 step 1 to the family walk, with
  the paper as the source of why they matter: at dT0/dC = 0 or lambda^M = 1 the continuation is in sqrt(eps) with two
  branches (Theorem 4.1) and a naive eps-walk will either miss one branch or jump. Report, for every row, the Melnikov
  zero count and slopes, not only existence. The rule "closure over N Sun periods with the Sun back at its starting phase"
  in `#905` is exactly M T0 = N T of the paper.
- `#884` (existing record). The Melnikov function in `sun_forced_periodic_884.py` is the paper's Eq. 3.3 specialised
  (section 4 of this note); the paper supports quoting Theorem 3.1 as the sufficiency statement. What it does not
  support: any claim that a non-simple zero or a degenerate orbit persists, and any count of forced orbits beyond
  the zeros found.
- `#902` (rebuild the bicircular validation tiers on a real orbit). The paper supplies two controls that depend on no
  project code, which a generic Melnikov or continuation routine (taking f, g and a period as inputs) can be tested
  against: Example 1 (the printed M and M_- and zeros (pi, +-1/2), (0, +-1/2)) and Example 2 (Eq. 7.4, with the corrected
  kappa_2 of section 3). The Example 1 numbers are reproduced by a finite-difference computation (section 3).
- `#916` (test the printed persistence conjecture of Ross & Roberts-Tsoukkas as invariant curves of the one-period
  map). For orbits that are exactly commensurate the paper gives the exact statement: the circle of fixed points of
  the stroboscopic map at eps = 0 breaks to isolated fixed points at the simple zeros of Mel. A direct test: for
  each stored commensurate orbit, count and locate the fixed points of F^N (the N-forcing-period map) near the
  unperturbed circle in the corrected model and compare with the Mel zeros. A mismatch is a bug in one of the two
  computations. For near-commensurate stable orbits (the conjecture proper) the paper gives nothing; the invariant-curve
  existence needs KAM-type arguments, and the project's `search/pertbp_strob_889.py` rotation-number machinery is the
  relevant tool.
- `#922` (`#890` orbit as an invariant two-torus). Theorem 3.1 applies only while the forcing has one period. Use it
  as the first rung: perturb the existing planar circular model (period 5 synodic periods) by one
  commensurate extra forcing, get the Melnikov zeros, then add the incommensurate frequency (ratio 1.830
  above) and treat the result as a torus problem (Rosales 2021, `search/pertbp_strob_889.py`). The paper's
  two-oscillator example (Example 2) is a toy version of a torus Z, which shows the structure that arises when the
  unperturbed period manifold is itself a torus (m-degenerate, up to 2^m continuation curves), not a single orbit.
- `#884` and the CCR4BP lanes. The same bifurcation function, with the same cokernel (grad C), applies to
  `core/ccr4bp.py` (one moon on a concentric circle), so it can be added once and reused.

## 7. Follow-ups (task numbers not registered here)

1. Implement `bifurcation_function(f, g, period, ...)` generically (any f, g, T0 = (N/M) T) and test it on Example 1
   (M, M_-, zeros) and Example 2 (Eq. 7.4 with the paper's printed zeros corrected as in section 3), as sourced
   expected values. Then specialise to the CR3BP-plus-Sun case and test agreement with the existing
   `melnikov_variational` on the stored `#884` orbits (the two should agree to quadrature error).
2. Add the nondegeneracy screens (dT0/dC, lambda^M) to the family-walk data of `#905` and record them per row.
3. Implement the generalized Melnikov function (4.2) for l = 1 and test on Example 1 before using it on a fold point.
4. For `#916`: the fixed-point count of F^N against the Mel zero count, per commensurate stored orbit.
5. Ask whether the paper's higher-order or square-root continuation (Theorem 4.1) explains any of the `#884` family
   members that the corrector found only at a particular Sun phase.
6. Cite correctly: the continuation theorem is Theorem 3.1 (nondegenerate) and 4.1 (degenerate) of this paper;
   the Melnikov-integral form for Hamiltonian f in section 4 is derived here and previously in the `#884` module, not
   printed in the paper.
7. Corpus: Cronin 1996 (Meth. Appl. Anal. 3:370), Chicone 1994 (J. Diff. Eqns. 112:407) and Farkas 1978 (SIAM J.
   Math. Anal. 9:876, existence regions) are the nearest sources for an existence radius in eps; not held.

## 8. Reading quality

Text layer present, equations garbled; I read pp. 92, 94, 100-103 from page images and cross-checked the
printed example numbers numerically. Digits I could not read with certainty: none in the equations used above.
The printed kappa_2 in Example 2 is inconsistent with the paper's own (7.4) (section 3).
