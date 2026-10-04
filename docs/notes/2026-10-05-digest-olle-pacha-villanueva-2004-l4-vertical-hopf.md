# Digest: Olle, Pacha and Villanueva 2004, motion close to the Hopf bifurcation of the vertical family of periodic orbits of L4

Date: 2026-10-05 (Sydney). Reading, reasoning and independent computations with the project's `core.cr3bp`; no project code or tests were changed.

Source: Merce Olle, Joan R. Pacha and Jordi Villanueva, "Motion close to the Hopf bifurcation of the vertical family of periodic
orbits of L4", Celestial Mechanics and Dynamical Astronomy 90 (2004), DOI 10.1007/s10569-004-1592-0, Universitat Politecnica de
Catalunya. Filed in the private paper corpus as
`olle-pacha-villanueva-2004-motion-close-hopf-bifurcation-vertical-family-l4-cmda-90-87-doi-10.1007-s10569-004-1592-0.pdf`.
Page numbers: the file name and the coordinator give 90:87-107; the scan of the paper prints "Celestial Mechanics and Dynamical
Astronomy 90: 89-109, 2004" on its first page and its running heads number the 21 pages 89 to 109, so I cite printed pages 89-109
below and flag that the file name uses the other numbering. Printed dates: received 31 October 2003, revised 26 March 2004, accepted
26 February 2004 (acceptance printed before revision; as printed). I read the whole text layer, and the key page (printed p.99) on its page image.
The text layer drops minus signs and some Greek letters (mu appears as "l", kappa and omega as "j" and "x"); I checked the signs below
against the page image of p.99 and against ordering constraints stated in the text.

Evidence tags: READ (p.N) is the printed page. COMPUTED is my own computation with the project's `core.cr3bp` (scratch scripts, recipe in
section 4). INFERRED is my reasoning.

Companions: `docs/notes/2026-10-05-digest-olle-pacha-villanueva-2005-normal-form-1-1-resonance.md` (the analysis paper with no example),
`docs/notes/2026-10-05-digest-jorba-olle-2004-invariant-curves-hamiltonian-hopf.md`,
`docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md`.

## 0. Answer first

This is the example the 2005 normal-form paper lacked. It is a numerical paper in the spatial circular restricted three-body problem
(RTBP): the vertical family of periodic orbits of L4 for mass parameters mu slightly above Routh's value (0.03852), where the family
passes from complex instability to stability as its vertical velocity grows, through a collision of two pairs of multipliers on the unit
circle: the Hamiltonian-Hopf transition. The critical orbit is printed for mu = 0.04 by its energy, rotation angle and one coordinate,
not by a full state. I reproduced it independently from `core.cr3bp` (section 4):
- h_crit = -1.4571360 (printed 1.45714146, with the sign lost in print, and 5.4e-6 apart), the critical rotation angle 1.83262879
  (printed 1.8326287), and the printed coordinate y = 0.86386304 of the critical orbit at x = -0.462 reproduce to all eight printed digits;
- the monodromy at the critical orbit is non-diagonalizable (one Jordan block per eigenvalue), the quantity Delta = alpha^2 - 4(beta - 2)
  changes sign linearly across it, and the family continues through it without loss of nonsingularity.
So the project has, for the first time, a sourced and independently reproduced test case for a sign-of-Delta detector and a "critical"
tag in a 3D restricted problem (section 6).

## 1. Setting (READ pp.91-93)

RTBP with masses mu in (0, 1/2] and 1 - mu at (mu - 1, 0, 0) and (mu, 0, 0), synodical frame, momenta p_x = xdot - y, p_y = ydot + x,
p_z = zdot, Hamiltonian (eq. 1)
H = (p_x^2 + p_y^2 + p_z^2)/2 + y p_x - x p_y - (1 - mu)/r1 - mu/r2,
r1^2 = (x - mu)^2 + y^2 + z^2, r2^2 = (x - mu + 1)^2 + y^2 + z^2. So the paper places the small primary at mu - 1 (negative x) and the
large one at mu; this is the project's convention (`core/cr3bp.py`: large at -mu, small at 1 - mu) reflected in x with time reversed
(COMPUTED/INFERRED: the map (x, t) -> (-x, -t) takes one form of the equations to the other, and takes the paper's L4 = (mu - 1/2, +sqrt3/2) to the project's L4;
monodromy multipliers are inverted, so Delta and the multiplier set are unchanged). H = -C/2 where C is the Jacobi constant
(`core.cr3bp.jacobi_constant`): COMPUTED H(L4) = -(3 - mu(1 - mu))/2 = -1.4808 for mu = 0.04, so all the paper's printed energies are negative (their signs were lost in the text layer and, in one place, in print: section 2.3 below).
Triangular points L4,5 = (mu - 1/2, +-sqrt3/2, 0); characteristic exponents (eq. 2) +-i (vertical oscillation, frequency 1) and the planar pair
sqrt(-1/2 +- (1/2) sqrt(1 - 27 mu (1 - mu))) (reconstructed from the standard result; the printed line is garbled in the text layer),
purely imaginary and different for 0 < mu < mu_R = (1/2)(1 - sqrt(23/27)) = 0.03852 (Routh), complex unstable for mu_R < mu <= 1/2.

Vertical family: the Lyapunov family of the vertical oscillation about L4 (frequency 1), parametrised by its vertical amplitude or by zdot at the crossing
of z = 0 (zdot = 0 is L4 itself); computed by arc-length continuation (Belbruno et al. 1994; Simo 1998; Gomez et al. 2001) as fixed points of the 4D Poincare map
P_h on R_h = {z = 0} intersect {H = h}, with P_h(p) = q the SECOND passage (pz obtained from H, and positive).

Stability parameters (Broucke 1969): characteristic polynomial of the monodromy p(z) = (z - 1)^2 (z^4 + a z^3 + b z^2 + a z + 1), with
a = -(lambda + 1/lambda + sigma + 1/sigma), b = 2 + (lambda + 1/lambda)(sigma + 1/sigma) (printed as a and b; these are the project's alpha and beta). Fig. 1
(p.93) plots, over (mu, zdot), the regions S (stable), Cu (complex unstable), Du (double unstable) and Su (semi-unstable) of the vertical family: for mu greater than and
close to mu_R, small-amplitude orbits are complex unstable (as L4 is) and "when decreasing zdot" the family goes from stable to complex unstable: a transition orbit exists
where the four multipliers collide in pairs on the unit circle, lambda = 1/lambda-bar = exp(i 2 pi kappa).

## 2. The paper's results (READ pp.93-107)

### 2.1 Rational collision (section 3, pp.94-95)
mu = 0.069608... (digits truncated in print) gives kappa = 1/4 at the transition: in Broucke's (a, b) diagram the path of the family C passes through
the marked point P4 (rational collisions kappa = 1/n at P1 ... P6). Two families of 4-periodic orbits QP1 and QP2 bifurcate (period about 4 T);
QP1 starts stable and later becomes unstable, all of QP2 is unstable. Printed (Fig. 3, p.95, minus signs restored): C member h = -1.116234, T = 6.335791;
QP1 member h = -1.124786, T = 25.337930; QP2 member h = -1.1411396, T = 25.309643.

### 2.2 Irrational collision (section 4, pp.95-107), mu = 0.04
- kappa = 0.291678572... "irrational", with 8167/28000 as a good rational approximation (p.96). Broucke diagram of the family, Fig. 4 (p.96).
- Summary (p.96): for h > h_crit the vertical orbit is linearly stable; around it there are Cantorian families of 3D KAM tori and the two Lyapunov families of
  elliptic 2D tori (invariant curves of P_h) born at the periodic orbit, one for each normal frequency; at h = h_crit the two families become one and for h < h_crit it
  detaches from the orbit as a single family of elliptic 2D tori at a finite distance from the (now complex unstable) periodic orbit. This is the DIRECT
  Hamiltonian-Hopf bifurcation; "for the vertical family the bifurcation has always direct character" (p.89-90): stable objects bifurcate on the unstable side. The direct
  character is shown numerically, by the stability of the bifurcated tori; a rigorous identification would need the nonlinear normal form (p.91).
- Invariant curves (section 4.1, p.98): parametrised X(theta) in R_h with P_h(X(theta)) = X(theta + omega), omega = 2 pi omega1/omega2 for the flow's two frequencies; Fourier
  series with N up to 50, Newton on the discretised F(X)(theta) = P_h(X(theta)) - X(theta + omega) on 2N + 1 mesh points, with one coordinate fixed at theta = 0 to remove the phase freedom.
- Section 4.2 (p.99): the curves are labelled by X(0) = (x~, y~, p~x, p~y) with x~ = -0.462 fixed (close to L4's x = mu - 1/2 = -0.46). Printed data (image of p.99 checked):
  h = -1.457018 (> h_crit): the periodic orbit has y~ = 0.86385432 and normal frequencies omega1 = 1.8028458, omega2 = 1.8625327; as h decreases to h_crit the two frequencies collide;
  at the critical orbit y~_crit = 0.86386304 and omega = omega_crit = 2 pi kappa = 1.8326287 (the "x" in the text layer is omega). The text prints "h_crit = 1.45714146" with no minus sign:
  with h = -1.457018 > h_crit this must be -1.45714146 (printed sign slip; the computation below gives -1.4571360).
- Fig. 6 (p.97): the two Lyapunov families in the (omega, y~) plane for h = -1.457018, and the detached family for h = -1.458308 < h_crit; Fig. 7 (p.100): invariant curves with omega = 1.90052495 (stable
  side) and omega = 1.79602495 (complex unstable side) at the same x~ (p.100).
- Linear normal behaviour (section 4.3, pp.100-103): reducibility of the linear skew product X' = A(theta) X, theta' = theta + omega, via the spectrum of the operator (L_omega w)(theta) = A(theta - omega) w(theta - omega)
  (Proposition 1: if lambda is an eigenvalue so is lambda exp(i k omega); Proposition 2: n omega-unrelated eigenvalues give a diagonal reduction); a curve is linearly stable if the spectrum lies on the unit circle and
  unstable if its closure is three circles of radius 1, 0 < c < 1 and 1/c; Fig. 8 (p.102) shows the spectrum of the discretised operator on the unit circle for N = 8 and N = 50 at h = -1.457018.
  Result: all Lyapunov families computed are linearly stable and persist, still linearly stable, on the complex unstable side (p.102-103), so the irrational collision gives the direct Hopf pattern.
- Dynamics close to the orbits (section 4.4, pp.103-107): (i) the 3D invariant manifolds of the complex unstable orbits, grown as 2D manifolds of the fixed point of P_h from the eigenvector of lambda_2 = r2 exp(i omega) (r1 < 1 < r2),
  initial curve of size c on the linear approximation, iterated by P_h; for h = -1.4577796221 < h_crit the 100th, 200th and 300th iterates (Fig. 9) go far away and come back repeatedly, and the slice x = x0 (Fig. 10)
  of the stable and unstable manifolds look coincident but are not (splitting of separatrices, p.105); the truncated integrable normal form has coinciding manifolds. (ii) confinement: a point near the manifolds stays near them for many
  iterates; Fig. 11 (p.106): 5000 iterates of a point on a 2D torus near the stable orbit h = -1.44908826 > h_crit, the same at h = h_crit = -1.45714146, and 20,000 iterates near the hyperbolic orbit at h = -1.45777962 with offset
  epsilon = 1e-6 (chaotic along the manifolds) and epsilon = 1e-4 (stays on what looks like a 2D torus). (iii) 3D KAM tori close to the periodic orbit change topology after the transition (secondary 3D tori from the bifurcated 2D tori), those far enough survive; Nekhoroshev-type effective
  stability is expected (Olle and Pfenniger 2000) but not for exponentially long times because the orbit is near a resonance.

## 3. Printed numbers usable as tests (READ page, with minus signs restored by the argument in section 1)

| item | value | page | status |
|---|---|---|---|
| mu_R | 0.03852 (from 0.5 (1 - sqrt(23/27))) | 92 | reproduces (COMPUTED 0.038520) |
| rational case mass parameter | 0.069608... (truncated), kappa = 1/4 | 94 | consistent with section 4 |
| C member | h = -1.116234, T = 6.335791 (mu = 0.069608) | 95 | T reproduced to 5e-6 by interpolation (section 4.3) |
| QP1 member | h = -1.124786, T = 25.337930 | 95 | not reproduced (needs the bifurcated branch) |
| QP2 member | h = -1.1411396, T = 25.309643 | 95 | not reproduced |
| irrational case mass parameter | mu = 0.04 | 96 | exact |
| critical rotation number | kappa = 0.291678572, 8167/28000 | 96 | printed value inconsistent with the printed omega_crit (section 4.2) |
| x~ (section coordinate of the curves) | -0.462 | 99 | used |
| h (stable side) | -1.457018: y~ = 0.86385432, omega1 = 1.8028458, omega2 = 1.8625327 | 99 | y~ to 1e-8; omega_i to 6e-6 (explained by the rounding of h) |
| h_crit | 1.45714146 as printed (sign lost), i.e. -1.45714146 | 99, 106 | computed -1.4571360299, difference 5.4e-6 |
| critical y~ | 0.86386304 | 99 | reproduces to all eight digits (0.86386304) |
| critical rotation angle | omega_crit = 2 pi kappa = 1.8326287 | 99 | computed 1.8326287865 |
| unstable-side h used for figures | -1.458308 (Fig. 6), -1.4577796221 (Fig. 9), -1.45777962 (Fig. 11), also h = -1.44908826 (stable) | 97, 104, 106 | multipliers in section 4.3 |
| curve rotation numbers | 1.90052495 (stable side), 1.79602495 (unstable side) | 100 | not reproduced (needs the curve solver) |
| numerical settings | N <= 50 Fourier modes; operator test N = 8 and 50; iterates 100, 200, 300 (manifolds), 5000 and 20,000 (tori); epsilon 1e-6, 1e-4 | 98, 102, 104, 106 | settings only |
| Fig. 1 (mu, zdot) map of S, Cu, Du, Su | graph only | 93 | not tabulated |

## 4. Independent reproduction with `core.cr3bp` (COMPUTED)

Recipe. Project frame (large primary at -mu), mu = 0.04. Start at L4 = (0.5 - mu, sqrt3/2), momenta from rest (px = -y, py = x), energy h. Fixed points of the 4D Poincare map on z = 0 at fixed
H = h, second passage, found by Newton on u = (x, y, px, py) (analytic 6 x 6 state transition matrix from `core.cr3bp.cr3bp_stm_eom`, DOP853 at rtol 1e-13, atol 1e-14; the Newton Jacobian was a central finite
difference of the map at 1e-7); continuation in h from H(L4) + 1e-4 upward, 30 steps to h = -1.4572. Monodromy M = the 6 x 6 state transition matrix over the period; the four non-trivial multipliers are the ones other than the pair
nearest 1; alpha, beta from their characteristic polynomial; Delta = alpha^2 - 4(beta - 2); the critical h by Brent root-finding of Delta (to 1e-14 in h). Residual of the fixed-point equations at every step: 1e-15 to 3e-14.
The paper's frame is the project's frame reflected in x (and time reversed), so the paper's x~ = -0.462 is the project's x = +0.462, y unchanged.

### 4.1 Critical orbit, mu = 0.04
- h_crit = -1.4571360299, period T = 6.2852869399, Jacobi constant C = 2.9142720597 (computed from -2 h), fixed point (x, y, px, py) = (0.46524528, 0.86271844, -0.84269655, 0.45444607) (project frame),
  state at the section z = 0: (x, y, z, vx, vy, vz) = (0.46524528, 0.86271844, 0, 0.0200219, -0.01079922, 0.2163659); z amplitude 0.21636644.
- alpha = 1.03540396, beta = 2.26801534, Delta = 2.9e-12 (zero to the root-finding tolerance); eigenvalues of the 6 x 6 monodromy: the trivial pair exactly 1 and the four others
  -0.25885147 +- 0.96591714 i and -0.25885051 +- 0.96591739 i (equal to 1e-6, the sqrt of the 1e-12 root tolerance, as expected for a double eigenvalue with a Jordan block); angle arccos(-alpha/4) = 1.8326287865 rad.
  Printed omega_crit = 1.8326287 (7 decimals): difference 8.7e-8 (rounding of the printed value). Rotation number computed omega/(2 pi) = 0.29167193.
- Non-diagonalizable: the singular values of M - lambda I (lambda = exp(i 1.8326287865)) are 34.17, 9.81, 1.587, 1.586, 0.315 and 8.5e-15: exactly one null direction, not two (a semi-simple eigenvalue of
  multiplicity 2 would give two). Condition number of the eigenvector matrix 8.8e7 (against 1e2 to 1e3 away from the critical orbit in the synthetic case of the 2005 digest).
- The orbit's point at the plane x = 0.462 (project frame; the paper's x~ = -0.462): y = 0.86386304 at z = +-0.03277, with the section coordinates printed by the paper reproduced to all digits: printed y~_crit = 0.86386304.
  (So the paper's invariant-curve coordinate X(0) is not on the plane z = 0 for the periodic orbit; INFERRED: the representative point is taken on the plane x = x~.)
- Sign change of Delta (the detector): Delta(h) from the scan, with Delta = 0 at h_crit and slope d(Delta)/dh = +112.6 per unit h (finite difference: Delta(h_crit +- 1e-5) = +-1.126e-3):
  h = -1.4573: Delta = -1.847e-2 (max |lambda| 1.0358); -1.4572: -7.205e-3 (1.0222); -1.4571: +4.057e-3 (1.0000, all four on the unit circle); -1.4570: +1.531e-2. Both sides at the same step 1e-4.
- Smooth continuation through the transition: the corrector converged at every step on both sides with residual 1e-15 to 1e-14, with the same corrector and no change of method; det(M - I) of the non-trivial block is |lambda - 1|^4 > 0 at lambda = exp(i 1.83), so the Jacobian stays nonsingular
  (INFERRED from the algebra; the computed convergence confirms it).
- Sensitivity: the multiplier modulus excess near the transition is max|lambda| - 1 = c sqrt(h_crit - h) with c about 2.78 (from h = -1.4572 and -1.4573: 0.02221 and 0.035796 against sqrt(6.4e-5) and sqrt(1.64e-4)); a modulus test with tolerance 1e-3 is therefore blind only within
  h_crit - h < (1e-3/2.78)^2 = 1.3e-7, that is |Delta| below about 1.5e-5, and the double eigenvalue itself is known only to about sqrt(root tolerance).

### 4.2 Consistency of the printed rotation number
The printed kappa = 0.291678572 and its rational approximant 8167/28000 = 0.291678571 disagree with the printed omega_crit = 1.8326287: 2 pi x 0.291678572 = 1.8326705, whereas 1.8326287 / 2 pi = 0.29167192 and my
computed value is 0.29167193. So the printed kappa (and the approximant built from it) appears to be a transcription slip of about 6.6e-6 (INFERRED: digits "...67857" for "...67192"); omega_crit = 1.8326287 is correct. Irrelevant to the paper's argument (any irrational value works), relevant to anyone using kappa as a test value.

### 4.3 Other printed energies, mu = 0.04 (stable side, hyperbolic side)
| h | T | multipliers (moduli, angles) | Delta | alpha, beta | note |
|---|---|---|---|---|---|
| -1.44908826 (Fig. 11, stable) | 6.286004009 | all four on the unit circle, angles 1.5906655 and 2.0827435 | +0.8836 | 1.019487, 2.038931 | stable, two distinct pairs |
| -1.457018 | 6.285297449 | unit circle, angles 1.8028516, 1.8625268 | + | | printed omega1 = 1.8028458, omega2 = 1.8625327 |
| -1.4577796221 / -1.45777962 | 6.285229641 | moduli 0.9326889, 1.0721688, angle 1.8322988 | -0.07262 | 1.036641, 2.286810 | complex unstable, fixed point (0.465105, 0.862807, -0.843329, 0.454603) |
| -1.458308 (Fig. 6) | 6.285182604 | moduli 0.910244, 1.098606, angle 1.8320274 | - | | complex unstable |
The printed omega1, omega2 at h = -1.457018 differ from the computed angles by +-5.8e-6 with the same sum (3.6653785 against 3.6653784): the splitting 0.0596869 printed against 0.0596752 computed. That is consistent with
the rounding of the printed energy: the splitting scales like sqrt(h - h_crit), and a change of 4e-7 in h (printed to six decimals) changes it by 0.2 percent as observed. The y coordinate of the orbit at x = 0.462 is 0.86385433 (printed 0.86385432).

### 4.4 The rational case, mu = 0.069608 (COMPUTED, coarse)
Same continuation for mu = 0.069608 (H(L4) = -1.4676). Delta changes sign between h = -1.12662 (Delta = -0.299, max |lambda| 1.146, alpha = +0.0013) and h = -1.12262 (Delta = +0.256, all on the unit circle,
alpha = -0.0011): the Delta = 0 point is at about h = -1.1245 with alpha close to 0 (the collision with lambda = +-i, kappa = 1/4 has alpha = 0, beta = 2); the printed mu is truncated, so alpha is zero only to about 1e-3 there.
The period of the family at the printed energy h = -1.116234, interpolated from the h = -1.118619 (T = 6.335454) and -1.114619 (T = 6.336028) steps, is 6.335796 against the printed 6.335791 (5e-6, within the interpolation). So the printed C member is on the stable side, just above the transition, as the text says (the QP orbits at -1.124786 and -1.1411396 lie on the unstable side).

## 5. Reconciliation with the project

- `core/cr3bp.py` (3D equations and 6 x 6 state transition matrix) reproduces every printed number I could test; no change is needed. I did not find a module that continues the vertical L4 family from L4 at fixed energy on the Poincare section z = 0 (files in `src/cyclerfinder/search` and `genome` that mention "vertical" include `jpl_family_census.py`, `earth_moon_resonant_families.py`, `variational_periodic_orbit.py`, `deflated_variational_periodic_orbit.py`, `genome/qp_tori_energy_walk.py` and `genome/known_corpus_3d.py`; I did not read them for this, so a continuation may exist that I missed); the 3D family tracer `search/cr3bp_3d_family_tracer.py` classifies a 6 x 6 monodromy with `_classify_floquet` (unit_tol = 1e-3). I wrote the scratch continuation myself.
- The project's `_classify_floquet` picks the trivial pair as the two eigenvalues nearest +1: at the critical orbit the other four are at angle 1.83 rad, so there is no ambiguity here; it reports "stable" within 1e-3 of the unit circle, which at h = -1.4572 (max |lambda| = 1.0222) correctly gives "unstable" but would call h_crit - h < 1.3e-7 "stable" (section 4.1).
- Conventions to carry: the paper's frame is the project's reflected in x (section 1); energies are H = -C/2 and negative.

## 6. Techniques applicable to the project's problems

### 6.1 `#931`: a sourced test case for the sign-of-Delta detector and the critical tag
This paper supplies the first independently reproduced Hamiltonian-Hopf critical orbit in a restricted problem. Concrete test, using `core.cr3bp` only:
1. Detector. Along the mu = 0.04 vertical L4 family continued in h (the recipe of section 4), Delta must change sign between h = -1.4572 and h = -1.4571 (computed -7.2e-3 and +4.1e-3, step 1e-4), with |alpha/2| = 0.5177 < 2 at the crossing, and the root at h = -1.457136 (printed -1.45714146, 5.4e-6 apart; use a tolerance of 2e-5 in h, three times the observed difference, because the printed value is the author's own computation to unknown accuracy). Tolerance on the critical rotation angle: arccos(-alpha/4) = 1.8326287 to 1e-6.
2. Critical tag. At the root: Delta within 1e-9 of zero, |alpha/2| < 2, one non-trivial singular value of (M - lambda I) near zero and the next at least 0.1 (computed 8.5e-15 and 0.315), eigenvector-matrix condition number above 1e6 (computed 8.8e7). A semi-simple double eigenvalue would give two small singular values; this is the discriminator the 2005 digest proposed, now on a real orbit.
3. Side. With the paper's family the complex-unstable side is the small-amplitude side (h < h_crit): the detector must report "stable above, complex unstable below" and the quartet moduli at h = -1.4577796221 must be 0.93269 and 1.07217 (angle 1.83230), at h = -1.458308: 0.91024 and 1.09861. These are my values, not printed; the printed anchors are the energies and the critical numbers.
4. Blind window of a modulus test (section 4.1): a classifier with tolerance 1e-3 on |lambda| calls the quartet "stable" for h_crit - h < 1.3e-7; a unit test should show Delta (sign) detecting the transition while the modulus test lags by that amount, which is invisible on any practical grid (step 1e-4) but matters when bisecting.
5. Negative controls from the same paper: the rational collision at mu = 0.069608 (alpha crossing 0) is also a Delta sign change with |alpha/2| < 2; the detector flags it as a Hopf-type collision with kappa = 1/4 (angle pi/2), where the theory predicts bifurcating 4-periodic orbits, not tori: the classifier should report the angle and mark the collision as rational when theta/(2 pi) is within a stated tolerance of a low-order rational (1/4 here). For the irrational case the computed kappa = 0.2916719 is 5e-6 from 7/24 = 0.2916667, so whether a decimal is called rational depends on the tolerance chosen; the printed 8167/28000 makes the same point.
Also useful: the family goes into and out of the Hopf-type collision smoothly in a fixed-energy Poincare-map corrector (`#899`/`#931` (d) logging).

### 6.2 `#916`: persistence of stable members as invariant curves
- The paper shows (numerically, mu = 0.04) that for stable vertical orbits there are two Lyapunov families of invariant curves (2D tori of the flow), that these persist and stay linearly stable on the complex unstable side, and that 3D KAM tori close to the orbit change topology across the transition while those far away survive (p.107). For `#916`'s test of a persistence conjecture on stable members: a stable member close to a Delta = 0 event must be expected to have its nearest tori reorganised when the family parameter crosses the transition; use Delta (or the distance in h to the transition) as a column and treat members close to the transition as a separate class.
- Control: the Jorba-Olle map digest gives exactly solvable test maps for the same transition (Table 1 and the printed curves reproduced), so the invariant-curve solver can be validated on the map before being applied to a restricted-problem section.
- Direct against inverse: for the vertical L4 family the transition is always direct (stable objects bifurcate on the unstable side); a family with an inverse transition (the map T_t of Jorba-Olle) has no stable objects on the unstable side and trajectories escape quickly: the project cannot assume the direct case for cyclers; it must be determined per family (here numerically from the stability of the bifurcated tori; analytically from the nonlinear normal form).
- The paper's Nekhoroshev remark (p.107): effective stability is expected for long but not exponential times near a resonant orbit, matching the 2005 digest's finding that the normal-form remainder is only R^(r_opt/2).

## 7. Recommended follow-ups (no task numbers registered)
1. Add the mu = 0.04 critical-orbit test of section 6.1 to the `#931` classifier tests (it needs a small scratch-to-`src/` fixed-point corrector on the section z = 0; or pin the continuation at 3 values of h: -1.4572, -1.4571 and the root).
2. Quote kappa as 0.29167193 (computed) and not 0.291678572 (printed) in any test; record the printed value as an erratum candidate.
3. Reproduce the rational case fully: the QP1 and QP2 branches (T = 25.337930 at h = -1.124786 and 25.309643 at h = -1.1411396, mu = 0.069608) need a period-4 branch-switching corrector; not done here.
4. Reproduce the invariant curves of Fig. 6 and 7 (omega = 1.90052495 and 1.79602495) with the Fourier-Newton method used on the map in the Jorba-Olle digest; the periodic-orbit inputs of section 4 are enough.

## 8. Summary
- The paper: numerical study of the direct Hamiltonian-Hopf transition of the vertical L4 family of the spatial RTBP at mu = 0.04 (irrational collision) and mu = 0.069608 (rational, kappa = 1/4), with 2D tori, their normal behaviour by operator spectra, invariant manifolds of the complex unstable orbits, confinement and secondary 3D tori; no table, a handful of printed numbers.
- Reproduced (project `core.cr3bp`): h_crit = -1.4571360 (printed 1.45714146 with the sign lost), omega_crit 1.8326287865 (printed 1.8326287), the critical coordinate y = 0.86386304 (printed, eight digits), y at h = -1.457018 to 1e-8, the rational-case period 6.335796 against 6.335791; non-diagonalizable monodromy at the critical orbit, Delta linear across it.
- Errata candidates: the sign of h_crit in the printed text (p.99, p.106), and kappa = 0.291678572 (inconsistent with the printed omega_crit by 6.6e-6); the file name's page number (87) differs from the printed pages (89-109).
- Project use: the first sourced, reproduced test case for the `#931` sign-of-Delta detector and critical tag, in 3D CR3BP, with the printed energies at the stable and unstable sides.
