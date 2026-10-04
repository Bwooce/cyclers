# Digest: Peng & Xu (2015), "Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted three-body problem"

Single-paper digest, written 2026-10-04 (Sydney time) from a read of all 25 pages of the complete copy (the earlier held copy was a
truncated 19 kB file; it was replaced the same day). Filed in the private paper corpus as
`peng-xu-2015-stability-two-groups-multirev-elliptic-halo-orbits-ERTBP-cmda-123-279-doi-10.1007-s10569-015-9635-2.pdf`.
The file has a text layer; the tables below were transcribed from the text layer and checked against the page images.
Statements are marked READ (printed in the paper), COMPUTED (arithmetic shown) or INFERRED (reasoning, not printed).
Journal page numbers are used throughout (PDF page n is journal page 278 + n).

## 0. Headline

- Peng, H. and Xu, S., Celestial Mechanics and Dynamical Astronomy 123:279-303 (2015), DOI 10.1007/s10569-015-9635-2. Authors:
  Hao Peng and Shijie Xu, Beihang University, Beijing. Received 27 April 2014, accepted 30 June 2015, online 23 July 2015.
- It is the method-and-stability paper behind every later multi-revolution elliptic halo (ME-Halo) result the project holds:
  Peng, Bai & Xu 2017 (same first author, Sun-Mercury), Neelakantan & Ramanan 2022, Park & Howell 2024, Singh, Park & Howell 2026.
- It prints NO initial condition, NO period table and NO orbit state anywhere. Its five tables (Tables 1 to 5) are monodromy
  eigenvalues only. So it supplies stability controls, not closure controls (section 7). This is the single most important fact
  for the project's checks.
- Model: Earth-Moon-type mass ratio sweep mu in [0.001, 0.020], eccentricity e in [0, 0.21], one resonance (M5N2) around L1,
  two orbit groups (periapsis and apoapsis), north halo only.
- Method: multi-segment (multiple-shooting) optimisation with an interior-point NLP solver, natural-parameter continuation in mu
  then in e at fixed step 0.001 (not arclength). Continuation fails above a critical eccentricity for small mu; the authors cannot
  tell whether this is a bifurcation or a fold ("either there is a bifurcation or the characteristic curve is not monotonic",
  p292). Peng, Bai & Xu 2017 later resolve exactly this with tangential continuation.
- Results the project can use: (i) the two groups have structurally different stability from the first order in e (a pair at +1
  splits into a real pair in one group and a unit-circle pair in the other, with nearly equal magnitude; section 6(c)); (ii) eleven
  eigenvalue-distribution types are observed, with collisions and bifurcations inside intervals as narrow as 1e-5 in e; (iii) the
  Earth-Moon orbits (mu 0.0122, e 0.0554) have largest multipliers of about 1.54e6 over 4 pi; (iv) printed eigenvalues at e = 0 for
  mu = 0.009 and 0.015 are an independent check on a circular-problem corrector (section 7).

## 1. The model exactly as printed (Section 2.1, pp281-283)

READ. Restricted problem, third body of negligible mass. Primaries m1 (larger) and m2, mu = m2/(m1+m2).

- Circular problem (Eq. 1, 2, p282): synodic frame, origin at the barycentre, x from m1 to m2, z along the angular momentum, y
  completing the right-handed set. Length unit = primary separation r12, mass unit = m1+m2, time unit = 1/n. m1 at x1 = -mu, m2 at
  x2 = 1-mu.
  x'' - 2y' = Omega_x, y'' + 2x' = Omega_y, z'' = Omega_z (dots are d/dt);
  Omega = (1/2)(x^2+y^2) + (1-mu)/r1 + mu/r2 + (1/2) mu (1-mu), r1 = sqrt((x+mu)^2+y^2+z^2), r2 = sqrt((x-1+mu)^2+y^2+z^2).
- Elliptic problem (Eq. 3 to 7, pp282-283): primaries on a Kepler ellipse of eccentricity e.
  r12(f) = a12 (1-e^2) / (1 + e cos f)                                                   (Eq. 3)
  d/dt = (df/dt)(d/df) = sqrt( G(m1+m2) / (a^3 (1-e^2)^3) ) (1 + e cos f)^2 d/df          (Eq. 4)
  x'' - 2y' = omega_x, y'' + 2x' = omega_y, z'' = omega_z (primes are d/df)               (Eq. 5)
  omega(x,y,z,f) = (1 + e cos f)^(-1) Omega_tilde(x,y,z)                                  (Eq. 6)
  Omega_tilde = Omega - (1/2) e cos f z^2                                                (Eq. 7)
  The synodic frame is "instantaneously normalized by r12(f), the total primary mass (m1+m2) and the reciprocal of the mean motion
  n-bar", so it pulsates and rotates non-uniformly. The independent variable is the true anomaly f. "the epoch when primaries are at
  their periapsis is set to be f0 = 0" (Fig. 1 caption and p283).
- No Jacobi integral: "Therefore, as a well-known fact, there does not exist the Jacobi integral in the ERTBP" (p283). The system is
  non-autonomous and 2 pi periodic in f.
- Libration points (Eq. 8, p283): omega_x = omega_y = omega_z = 0 gives five fixed points at the same place as in the circular
  problem; they are "geometrical libration points but not dynamical ones", the collinear ones oscillating along the x-axis in the
  pulsating frame.
- Symmetry (Eq. 9, p284, after Sarris 1989): (f; x, y, z, x', y', z') -> (-f; x, -y, z, -x', y', -z'), i.e. symmetry about the x-z
  plane combined with time reversal. The authors say Sarris discussed other symmetries but this one is the one used.
- Parameters used: mu in [0.001, 0.020] (Earth-Moon 0.0122 marked in every surface plot), e in [0, 0.210] (Earth-Moon 0.0554), steps
  delta mu = delta e = 0.001 (p290-291). Only the north halo branch is studied (p289).

COMPUTED check against the project. `core/er3bp.py:er3bp_eom` has x'' = 2y' + s (x - grav_x), y'' = -2x' + s (y - grav_y),
z'' = -s (e cos f z + grav_z), s = 1/(1 + e cos f). From Eq. 5 to 7: omega_x = Omega_x / (1+e cos f); omega_z = (Omega_z - e cos f z)
/ (1+e cos f) with Omega_z = -grav_z. These agree term by term. READ for the project's code, COMPUTED for the match. The constant
(1/2) mu (1-mu) in Omega drops out of the equations of motion.

INFERRED consequence of Eq. 6: the equations depend on f only through e cos f, so replacing f by f + pi is the same as replacing e
by -e. The paper does not state this; Peng, Bai & Xu 2017 do (their Eq. 6, see the 2017 digest), and the project's `er3bp_eom` takes
any real e, so an apoapsis-group orbit can be run from f = 0 with a negative e (not tested here).

## 2. The orbits

### 2.1 Definition (Section 2.2, pp283-284)

READ (p284): "For an orbit to be periodic it is sufficient that it has two perpendicular crossings with the syzygy-axis, and that
the crossings happen at moments when the two primaries are at an apse, (i.e., at maximum or minimum elongation, or apoapsis and
periapsis)." (Moulton 1920, Broucke 1969; planar case.) Spatial version after Campagnola (2010): "it is sufficient that it has two
perpendicular crossing with either the normal plane or the syzygy axis, or both of them, when the primaries are at apse."

Commensurability (Eq. 10, p284): T_E = N 2 pi = M T_C, with M, N positive integers, T_E the ERTBP period and T_C the CRTBP period,
both measured in true anomaly (equal to scaled time in the circular problem). "M indicates the revolution number of the third body
around the libration region and N indicates the revolution of the primaries within one orbit period. So actually it is an M : N
resonant orbit." Campagnola first generated these (Sun-Mercury and Earth-Moon, Campagnola et al. 2008) and called them elliptic
halo orbits; here they are called ME-Halo orbits. The reason for M > 1: the circular halo period at small mu is too narrow to reach
2 pi N with M = 1 (Fig. 2, p285), so the orbit must close only after M revolutions.

The authors keep M and N small: p285 chooses M = 5, N = 2 ("M5N2"). "too large M and N will cause numerical difficulties for stability
study" (p285; quote in section 4.1 below).

Fig. 2 (p285), Earth-Moon mu = 0.0122: halo families around L1 (left plot) and L2 (right plot) as z-amplitude Az against period
T_C, with vertical dashed lines at commensurate periods, labelled M3N1, M5N2, a third label read as M8N3 (low image resolution,
unclear) and M7N3 on the L1 plot, and M7N3, M9N4, M2N1 on the L2 plot. COMPUTED periods T_C = 2 N pi / M: M3N1 2.0944, M8N3 2.3562
(if the label is right), M5N2 2.5133, M7N3 2.6928, M9N4 2.7925, M2N1 3.1416. These fall in the plotted axes (L1 1.8 to 2.8, L2 2.6 to
3.4), which supports the readings. No Az values are printed in text.

### 2.2 The groups (Section 2.4, pp288-289)

Section title: "Four groups of ME-Halo orbits". Because the continuation is in mu and e at a fixed (M, N), the orbits "can be continuously
parameterized by mu and e", so the word "group" rather than "family" is used (the circular-problem halo orbits form a family
parameterised by Az; the elliptic ones exist only at discrete (mu, e) for each (M, N)).

READ (p288). "According to the periodicity criterion, ME-Halo orbits perpendicularly cross the x-z plane twice and these two crosses
can occur whether primaries are at periapsis or apoapsis."

- M odd: "after half a period, two crosses occur at two sides of ME-Halo, so the orbit is different whether it starts at f0 = 0 or
  f0 = pi, which indicates the primaries are at their periapsis or apoapsis. Let the orbit starts from the left side (in the x-y
  projection) of the ME-Halo orbit, define Periapsis Group to start at f0 = 0 and Apoapsis Group to start at f0 = pi. If the orbit
  starts from the right side, it will coincide with one of the two groups depending on that N is even or odd." Fig. 4 (p288):
  Earth-Moon L1 Periapsis and Apoapsis ME-Halo, M5N2, with Earth-Moon eccentricity "adopted as 0.0554"; the perpendicular crossing is
  marked with a red circle, on the outer circle of the y-z projection for periapsis and on the inner circle for apoapsis.
- M even: "two crosses occur at the same side of ME-Halo orbit and the orbit can be different whether the cross is on the left or
  right side of the orbit. Let the orbit starts at f0 = 0, define Left Group to start from the left side and Right Group from the
  right side. If the orbit starts from f0 = pi, it will coincide with one of the two groups as well." Fig. 5 (p289): Earth-Moon L2
  Left and Right ME-Halo, M2N1; "The most obvious difference is the position of the bifurcation of the orbit. For Left Group it is
  at the top and for Right Group the bottom."
- North and south halos are mirror images; only the north branch is used.

INFERRED parity bookkeeping (not printed). The half period is T_E/2 = N pi, so the second perpendicular crossing is at f0 + N pi.
For N even both crossings are at the same apse (M5N2: both at periapsis for the periapsis group, both at apoapsis for the apoapsis
group); for N odd the two crossings are at opposite apses (M3N1, M2N1). Starting from the other crossing of a given orbit relabels
it, which is the paper's "will coincide with one of the two groups depending on N". The classification is stated for coprime (M, N);
INFERRED that M4N2 (not coprime, T_C = pi, the same circular orbit as M2N1) is the 2-fold cover of the M2N1 circular orbit and the
paper's two-group count does not obviously apply to it (section 6(e)).

## 3. Method: generating an ME-Halo orbit (Sections 2.2 to 2.3, pp284-288)

### 3.1 Step 1: a circular halo of the right period (p285)

READ. "The very first step of the continuation is to find a halo orbit with precise period T_C satisfying Eq. (9)". (The printed
cross-reference is to Eq. 9, the symmetry map; the commensurability condition is Eq. 10. INFERRED cross-reference slip.) "an orbit is numerically extracted from the whole
halo orbit family by a dichotomy of Az". The families are made "by continuations along Az with simple differential correction method".
Then T_C = 2 N pi / M at e = 0, the initial condition X0(e = 0) is taken, and e is increased by delta e with X0(e + delta e)
obtained from the multi-segment optimisation. Continuation in mu works the same way. "There are infinite sets of (M, N) close to any
desired period since the rational number is dense."

### 3.2 Circular-problem differential corrector (Eq. 11 to 15, pp285-286)

READ. State X = [x, y, z, x', y', z'] (dots in the circular problem). Symmetric halo: X0 = [x0, 0, z0, 0, y0', 0]^T, half-period
state must be [x^, 0, z^, 0, y^', 0]. Error at half period delta X_T/2 = X^_T/2 - X_T/2; first-order update
0 = X^_T/2 - X_T/2 = (dX_T/2/dX0) delta X0 + (dX_T/2/dT) delta T + O(.)            (Eq. 11)
with Phi the state-transition matrix, Phi' = A_C Phi, Phi(t0, t0) = I6 (Eq. 12), A_C = [[0, I3],[H_C, K_C]], H_C the Hessian of
Omega, K_C = [[0, 2, 0], [-2, 0, 0], [0, 0, 0]] (Eq. 13). The reduced system (Eq. 14, 15, with the propagation stopped at y1 = 0 so
that delta y = 0):
[delta x1'; delta z1'] = [[phi41, phi43, phi45],[phi61, phi63, phi65]] [delta x0; delta z0; delta y0']
                         - (1/y1') [x1''; z1''] [phi21, phi23, phi25] [delta x0; delta z0; delta y0'].
The continuation is carried out "with z-axis amplitude Az", so z0 is fixed and the correction is [delta x0, 0, 0, 0, delta y0', 0]
(Eq. 15). The authors say this method has a small convergence domain, will fail for a long integration, and can converge to a
retrograde distant orbit because it "does not control the boundary of the correction". It is used only for the circular families.

### 3.3 Multi-segment optimisation for the elliptic problem (Eq. 16 to 20, pp286-287; Fig. 3)

READ. This is the paper's key technique.
- Initial state at epoch f0 (Eq. 16): X0 = X(f0) = [x0, 0, z0, 0, y0', 0]^T (primes are d/df). State after half the ERTBP period
  (Eq. 17): X_T/2 = X(f0 + T_E/2) = [x, y, z, x', y', z']^T at f0 + T_E/2. "in the ERTBP the independent variable at the same time
  indicates the phase angle of the primaries, so not only the duration T_E/2 but also the initial epoch f0 must be specified". f0 = 0
  for the periapsis group, f0 = pi for the apoapsis group (p289).
- Terminal condition (Eq. 18): y = x' = z' = 0 at T_E/2, i.e. perpendicular crossing of the x-z plane.
- "Borrowing the idea of the multiple shooting method introduced by Howell and Pernicka", the orbit is broken into n segments;
  X_i and f_i (i = 0 ... n-1) are the starting state and epoch of segment i+1. Continuity (Eq. 19): X^-_{i+1}(X_i, f_i) = X_i for
  i = 0 ... n-2, where X^-_{i+1} is the end of segment i integrated from f_i to f_{i+1}; that gives 6(n-1) nonlinear constraints.
  X_n = [x_n, y_n, z_n, x_n', y_n', z_n'] is the end point at f0 + T_E/2; the remaining three components are free.
- Cost (Eq. 20): min J(X0, f0; ...; X_i, f_i; ...; X_{n-1}, f_{n-1}) = sqrt( y_n^2 + x_n'^2 + z_n'^2 ). "The minimum of this problem
  should be zero and it implies a periodic orbit." Fig. 3 (p287) illustrates n = 5 segments with an example step from (mu0 = 0.0122,
  e0 = 0.038) to (mu1 = 0.0122, e1 = 0.088) in the x-y projection (the figure shows the nodes X0 to X5 and epochs f1 to f4; the parameter values are read from the figure legend at low resolution).
- Solver: Matlab fmincon, interior-point method by default (tolerances can be set separately for position and velocity); SQP "is also
  utilized in some studies, for example when trying to re-close up a continuation result with higher accuracy" and "appears faster
  than interior point method when a good initial guess is available". "Extrapolation techniques for initial guesses show some
  improvements but not that much efficient." (p287-288). The previous ME-Halo orbit is the default initial guess.
- Number of segments (p288): "The convergence of the algorithm increases as n increases, but the time cost increases as well. So
  after a trial and error process, we recommend M/2 <= n <= 2M for continuations along e and n >= 3M for continuations along mu."
  In the stability sweep (p291): about 32 segments for one step in mu, about eight for one step in e (M = 5); about 170 s per mu step
  and 60 s per e step on a 3.2 GHz PC, serial.
- Tolerances (p291): Matlab ode45, absolute and relative 3e-14 (normalised units); nonlinear constraint tolerance 1e-10, "which is
  the upper boundary of all discontinuities at the connection points between two segments"; "After obtain databases with thousands of
  ME-Halo orbits, the eigenvalues ... are calculated". Continuation: first along mu at e = 0, then along e for each mu, in parallel.
- Why segments (p286): the differential corrector "has a convergence domain small for a nonlinear problem" and "a long-time
  integration of nonlinear equations will almost certainly fail the method even a good initial guess is available". The authors recast
  shooting "as an optimization problem and solve it with mature solvers" instead of building a more advanced method.
- Limits the authors state about their own continuation (p291, p292): "all successful continuations can converge with any accuracy no
  higher than the integration accuracy, but any continuation covering e^(mu) cannot converge to the tolerance"; the step is simply
  reduced until it is below the integration accuracy.

Not printed: the number of unknowns assembled explicitly, whether the interior epochs f_i are truly free in the solver (the text lists
(X_i, f_i) as arguments of J, so they appear to be), the Jacobian method (fmincon default, finite differences, INFERRED; no STM
Jacobian is mentioned for the corrector), and any initial states.

## 4. Stability (Section 2.5 pp289-290; Sections 3.1 to 3.4 pp290-301)

### 4.1 Definition and computation

READ (p289-290).
- Monodromy matrix = state transition matrix over one period. In the ERTBP Phi'(f, f0) = A_E(X, f) Phi(f, f0), Phi(f0, f0) = I6, with
  A_E = [[0, I3], [H_E, K_E]], H_E the Hessian of omega, K_E as in Eq. 13 (Eq. 21).
- "The periodic orbit in a periodic system is stable if and only if all eigenvalues of the monodromy matrix Psi(f) = Phi(f + T_E, f)
  have modules smaller than one (Bittanti and Colaneri 2009)." Eigenvalues are invariant along the orbit, "they can be calculated at
  any convenient point". Psi(f0) is the STM from f0 to f0 + T_E with f0 = 0 for the periapsis, left and right groups and f0 = pi for
  the apoapsis group. (INFERRED comment: with every eigenvalue of a symplectic monodromy having a reciprocal, "all modules smaller
  than one" cannot hold, so the operative criterion is "all on the unit circle", as the stability index below makes explicit.)
- Multi-segment product (Eq. 22): Psi(f0) = Phi(f0 + T_E, f0) = Phi(f0 + T_E, f_{n-1}) ... Phi(f2, f1) Phi(f1, f0). "If directly
  integrate for one period T_E from the very first initial condition, the initial deviation will be magnified exponentially, and the
  deviation it caused to the monodromy matrix can affect its eigen-structure. But if Psi is calculated separately according to Eq.
  (22), the numerical process appears more steady." In the sweep each orbit is first refined to four segments, with the same
  tolerances (p291).
- The stated limit, verbatim (p290): "However, as M and N grow larger, the total integration time increases, and this causes the
  largest eigenvalue become larger exponentially and the smallest become smaller. Very quickly they will differ more than 20 orders,
  and this will cause the Psi intrinsically not able to be numerically obtained will enough accuracy. So this is a crucial
  restriction that the ME-Halo orbits with N = 2 are studied, even though we could generate ME-Halo orbits with a larger N. The
  authors suppose a more profound insight into the stability of these orbit requires more advanced mathematic methods dealing with
  the nonlinearity."
- Eigenvalue structure: "In the CRTPB, the monodromy matrix of a halo orbit usually has a pair of eigenvalues equal to one because
  the system is Hamiltonian and autonomous. But in the ERTBP this is not the truth because of the appearance of eccentricity e and
  the non-autonomous features (Broucke 1969)." Eigenvalues come in reciprocal pairs lambda1, 1/lambda1, lambda2, 1/lambda2, lambda3,
  1/lambda3 (Campagnola 2010).
- Stability index (Broucke 1969, Sarris 1989), k_i = lambda_i + 1/lambda_i, i = 1, 2, 3, with the modification (Eq. 23, p290)
  k_i = 2 max(|lambda_i|, |1/lambda_i|) if lambda_i is complex and |lambda_i| not equal to 1, and k_i = lambda_i + 1/lambda_i otherwise.
  Reason: "The only exception is that if there are two pairs of real reciprocal complex eigenvalues which are conjugate but not on the
  unit circle, they will give complex k_i." With the modification "two pairs of conjugate complex eigenvalues not on the unit circle
  will give the same index larger than 2". The text says the orbit is unstable "if any k_i > 2"; the figures draw the lines k = 2 and
  k = -2 (Fig. 7 and 8). INFERRED: the operative criterion is |k_i| > 2 (negative real eigenvalues give k below -2). The authors note
  that the index is discontinuous where eigenvalues leave the unit circle and that Broucke 1969 observed the same "complex instability"
  in the planar problem.

### 4.2 Results, periapsis group, L1, M5N2 (Section 3.1, pp291-296)

READ. Characteristic curves of the orbits in the subspace (x0, z0, y0') (Fig. 6, p291; axes read as x0 about 0.8 to 0.95, z0 about
0.1 to 0.25, y0' about 0.1 to 0.3; the figure prints no numbers beyond the axes). For each mu the curve is parametrised by e from 0
(green) to 0.21 (red). "At mu^ = 0.012 the curve is the shortest, and it separates the two trends. For mu > mu^, y0' of the curves
keep descending, but for smaller mu < mu^ the curves bend in the middle, y0' of the curves start to go upward and stop at a critical
eccentricity e^(mu) < 0.21. The continuation algorithm fails to converge at e^(mu)." "The ME-Halo orbits later displayed in Fig. 10
do not show any obvious mechanisms that might cause this problem, like a close flyby, so it can be inferred that either there is a
bifurcation or the characteristic curve is not monotonic with respect to the continuation parameter" (p292). "the stability
analysis later shows that the eigenvalues of the ME-Halo orbit near e^(mu) have great changes."

Stability index surfaces (Figs. 7, 8, 9, p292; text p293): "There is a clear and big great gap at mu^ = 0.012 in all figures, and the
trends on two sides of the gap are different. When mu < mu^, k1 decreases, k2 increases and k3 increases with e, but when mu > mu^,
k1, k2 and k3 all increase with e." For mu <= mu^ the curves k1(e) drop below -2 at e^(mu) (nearly vertically, "a tiny incensement of
e will cause very large variation of k1"); k3 is always below 2 (unit-circle pair), and for mu = mu^ the k3 curve "stops almost exactly
at k3 = 2", i.e. an eigenvalue almost on the unit circle at the end.

Two branches are chosen as representatives, "lighter branch" (mu = 0.009 < mu^) and "heavier branch" (mu = 0.015 > mu^). Shapes
(Figs. 10 and 11, pp293-294) in the pulsating (red) and non-pulsating (blue) frames: the pulsating-frame loops shrink with e; in the
non-pulsating frame the orbits are alike up to e = 0.12 and differ after that.

- Lighter branch, Table 1 (mu = 0.009, e = 0 to 0.16, Fig. 10). READ (p294): "the lighter branch shows a big change from e = 0.12 to
  0.16, where lambda1 starts to decrease and finally falls below zero and lambda2 becomes complex at the same time." A fine scan with
  delta e at least 2.5e-6 finds that "the eigenvalues collide and bifurcate three times" (Fig. 12, p295): intervals read from the
  figure (small print) e in [0.15818, 0.15819] (a), [0.15819, 0.15820] (b), [0.15834, 0.15835] (c). Between (a) and (b) "the ME-Halo
  orbit has all eigenvalues complex but two pairs of them conjugate with respect to the unit circle", between (b) and (c) "two pairs
  of negative eigenvalues". The orbits before and after do not look different; "We suggest this phenomenon is caused by the
  interaction between mu and e, which cannot likely be explained by numerical studies."
- Heavier branch, Table 2 (mu = 0.015, e = 0 to 0.21, Fig. 11): lambda1 grows steadily, lambda2 grows; no collision; "when mu is
  larger, the effect of the eccentricity is relatively smaller."
- Critical orbit (Fig. 13, p296): mu^ = 0.012 and e^(mu^) = 0.143 are the centre of the figure; at mu = 0.012 the continuation fails
  at about e = 0.143 (the orbits at the continuation end are shown). For small mu "the loops in the middle shrink greatly compared
  with their beginnings", for larger mu they keep their size.

### 4.3 Results, apoapsis group, L1, M5N2 (Section 3.2, pp296-300)

READ. All curves extend to e = 0.21 successfully, unlike the periapsis group. The y0' curves keep descending for mu >= 0.011 and bend
and rise for mu <= 0.010, so a separator mu^ in (0.010, 0.011) "separating these two trends as well. But it was not directly
detected by the discrete grids." k1 is always above 2 (lambda1 and 1/lambda1 real positive) (Fig. 15, p297); the k1 curves for
mu < mu^ "fall down to as small as 500, but turn upward after a certain eccentricity e^"; k2 and k3 (Fig. 16, p298) are separated for
mu < mu^ and mu > mu^, coincide in intervals, and show discontinuities where lambda2 and lambda3 become complex but not on the unit
circle. Three kinds of stability evolution: lighter (mu = 0.004, Table 3, Fig. 17), medium (mu = 0.009, Table 4, Fig. 18), heavier
(mu = 0.015, Table 5, Fig. 19). Lighter and medium branches have complex-off-circle intervals: lighter only about e = 0.10 to 0.14;
medium around e = 0.12 and e = 0.14 to 0.16; then both change to a pair of negative real eigenvalues. "The heavier branches has three
pairs of eigenvalues, and they all show increasing trends ... which indicates that the ME-Halo orbit is very unstable with these
parameters, all states near in the phase space will eventually leaving it along the three-dimensional unstable manifold quickly."

### 4.4 Earth-Moon orbits (Section 3.3, p300)

READ. mu = 0.0122, e = 0.0554. "The L1 Periapsis ME-Halo orbit has two pairs of real eigenvalues, (lambda1, 1/lambda1) and (lambda2,
1/lambda2), and a pair of complex unit eigenvalues (lambda3, 1/lambda3). lambda1 is very large, about 1.5427 x 10^6, but lambda2 =
1.0086 is only slightly bigger than 1. The L1 Apoapsis ME-Halo orbit has only one pair of real eigenvalue (lambda1, 1/lambda1), and
two pairs of complex unit eigenvalues (lambda2, 1/lambda2) and (lambda3, 1/lambda3), where lambda2 = 1.5431 x 10^6 is larger than
that of the Periapsis." "Although the ME-Halo orbit will not be precisely periodic in the ephemeris model, it can be seen that the
appearance of eccentricity has dynamically changed the stability properties. Generally speaking, the orbits in the ERTBP appear more
unstable than that in the CRTBP."
INFERRED slip: the 1.5431e6 value must be the single real pair, i.e. lambda1 (the text labels it lambda2, but a pair of complex unit
eigenvalues cannot be 1.5e6).
COMPUTED plausibility: 1/1.5427e6 = 6.4821e-7; 1/1.0086 = 0.99147. Interpolating Tables 1 and 2 linearly in mu at e = 0.06 gives
lambda1 about 1.3615e6 + (0.0032/0.006)(1.6845e6 - 1.3615e6) = 1.534e6, and the lambda2 - 1 = 0.0098 (mu = 0.015) at e = 0.06
against 0.0086 printed at e = 0.0554: scaling the two tabulated values 0.0119 (mu = 0.009) and 0.0098 (mu = 0.015) by
(0.0554/0.06)^2.5 = 0.8192 (the e^(5/2) behaviour of section 6(c)) gives 0.00975 and 0.00803, and interpolating to mu = 0.0122 gives
0.0088, within 2 percent of 0.0086. Both are consistent with the printed Earth-Moon numbers.

### 4.5 Eigenvalue types, collisions and bifurcations (Section 3.4, p301; Fig. 20)

READ. Fig. 20 draws eleven eigenvalue distributions on the complex plane ("Type 1" to "Type 11", plus the circular halo). "It can be
observed in other plots of Fig. 20 that an ME-Halo orbit can have as many as three pairs of real eigenvalues." Observed sequences
(read from the page; the order of numbers is as printed, the figure itself shows only small sketches):

| Branch | Type sequence as e grows |
|---|---|
| L1 Periapsis, lighter | 1 -> 2 -> 3 -> 4 -> 5 -> 6 |
| L1 Periapsis, heavier | 1 |
| L1 Apoapsis, lighter | 8 -> 9 -> 8 -> 10 -> 11 -> 10 |
| L1 Apoapsis, medium | 8 -> 9 -> 8 -> 9 -> 8 -> 10 -> 11 |
| L1 Apoapsis, heavier | 8 -> 9 -> 7 |

(The two longer apoapsis sequences were read from the page image; digits unclear at the low resolution of that line, so treat the exact
order as approximate. The periapsis rows are clear.) Also READ (p301): "Campagnola had observed that the stability of left and right
ME-Halo orbits with M2N1 in the Earth-Moon system bifurcates at e = 0 (Campagnola et al. 2008; Campagnola 2010)." "Generally the
heavier branches show less complexity. On the other hand, the eccentricity e usually has relatively greater effects on the stability
property of ME-Halo orbits when mu is small. This can be preliminarily explained by the Legendre polynomial expansion of the
equations of motion (Lei et al. 2013). The eccentricity e only arises from the second order terms, it will have less effect, but when
it is comparable with effect of the mass ratio mu, for example if e is 10 times of mu it might show more severe effect. However, a
quantified result can hardly be drawn from the present numerical study."

### 4.6 Conclusions (pp301-302), verbatim

"ME-Halo orbits are strictly periodic orbits in the nonautonomous ERTBP model which revolves M circles around the libration point
region in one period when primaries revolve N circles." "The multi-segment optimization method is proposed to accomplish the
continuation because the traditional differential correction method diverges." "Totally eleven distribution types of eigenvalues
are observed in this study." "According to the numerical exploration, the Periapsis Group has two different branches, and the
apoapsis group has three different branches. The emergence of the eccentricity in the system introduces great complexity. An ME-Halo
orbit can have as much as three pairs of real eigenvalues, one or two pairs of negative eigenvalues, and two pairs of complex
eigenvalues out of the unit circle. Also, eigenvalues of particular orbits are negative. These are all very different from halo
orbits in the CRTBP." "In the Periapsis Group, a continuation barrier arises for small mass ratio, which seems to be caused by the
change of eigenvalues from positive to negative but still needs more analytical studies in the future." "Since the ME-Halo orbit
captures more natural dynamics, its specific stability features can provide potentially practical applications. These properties will
be helpful especially in fast systems like the Earth-Moon system or very eccentric systems like the Sun-Mercury system."

## 5. Tables, transcribed from the text layer (digit by digit)

Conventions. The paper's columns are e, lambda1, 1/lambda1, lambda2, 1/lambda2, lambda3, 1/lambda3, with the eigenvalues of the
monodromy matrix Psi over T_E = 4 pi (M5N2), L1, north halo, units are nondimensional (eigenvalues are pure numbers). e = eccentricity
of the primaries. No row prints an initial condition. Numbers are copied as printed ("e+06" is the paper's notation for a power of
ten; a minus sign before a number is the paper's U+2212). For complex pairs the paper prints the pair and its conjugate (or, for off-circle
pairs, the reciprocal), and the columns are copied in the same way. The authors' own caveat (p294): at e = 0 the unit eigenvalue is
printed as 1.0002 (Table 1) and 1.0000 + 0.0001i (Table 2), "caused by numerical errors during long-time integrations, and the data
suggest that there are at least four significant digits."

### Table 1 (p294): periapsis group, lighter branch, mu = 0.009 (the caption does not print mu; Fig. 10 gives mu = 0.009)

| e | lambda1 | 1/lambda1 | lambda2 | 1/lambda2 | lambda3 | 1/lambda3 |
|---|---|---|---|---|---|---|
| 0 | 1.3112e+06 | 7.6258e-07 | 1.0002 | 0.9998 | 0.9587 + 0.2843i | 0.9587 - 0.2843i |
| 0.02 | 1.3169e+06 | 7.5937e-07 | 1.0008 | 0.9993 | 0.9607 + 0.2777i | 0.9607 - 0.2777i |
| 0.04 | 1.3338e+06 | 7.4995e-07 | 1.0043 | 0.9957 | 0.9662 + 0.2577i | 0.9662 - 0.2577i |
| 0.06 | 1.3615e+06 | 7.3453e-07 | 1.0119 | 0.9882 | 0.9745 + 0.2245i | 0.9745 - 0.2245i |
| 0.08 | 1.3978e+06 | 7.1573e-07 | 1.0250 | 0.9756 | 0.9840 + 0.1782i | 0.9840 - 0.1782i |
| 0.10 | 1.4324e+06 | 6.9841e-07 | 1.0452 | 0.9567 | 0.9927 + 0.1206i | 0.9927 - 0.1206i |
| 0.12 | 1.3515e+06 | 7.3979e-07 | 1.0778 | 0.9278 | 0.9976 + 0.0696i | 0.9976 - 0.0696i |
| 0.14 | 8.0075e+05 | 1.2486e-06 | 1.2073 | 0.8283 | 0.9929 + 0.1190i | 0.9929 - 0.1190i |
| 0.16 | -2.1676e+05 | -4.6133e-06 | 0.8705 + 0.4922i | 0.8705 - 0.4922i | 0.9896 + 0.1440i | 0.9896 - 0.1440i |

### Table 2 (p295): periapsis group, heavier branch, Fig. 11 (mu = 0.015)

| e | lambda1 | 1/lambda1 | lambda2 | 1/lambda2 | lambda3 | 1/lambda3 |
|---|---|---|---|---|---|---|
| 0 | 1.5966e+06 | 6.2662e-07 | 1.0000 + 0.0001i | 1.0000 - 0.0001i | 0.9021 + 0.4316i | 0.9021 - 0.4316i |
| 0.02 | 1.6061e+06 | 6.2247e-07 | 1.0006 | 0.9994 | 0.9052 + 0.4250i | 0.9052 - 0.4250i |
| 0.04 | 1.6351e+06 | 6.1197e-07 | 1.0035 | 0.9965 | 0.9142 + 0.4053i | 0.9142 - 0.4053i |
| 0.06 | 1.6845e+06 | 5.9376e-07 | 1.0098 | 0.9903 | 0.9282 + 0.3721i | 0.9282 - 0.3721i |
| 0.08 | 1.7563e+06 | 5.6919e-07 | 1.0202 | 0.9802 | 0.9457 + 0.3250i | 0.9457 - 0.3250i |
| 0.10 | 1.8536e+06 | 5.3984e-07 | 1.0359 | 0.9654 | 0.9647 + 0.2635i | 0.9647 - 0.2635i |
| 0.12 | 1.9825e+06 | 5.0427e-07 | 1.0594 | 0.9439 | 0.9824 + 0.1870i | 0.9824 - 0.1870i |
| 0.14 | 2.1633e+06 | 4.6247e-07 | 1.1003 | 0.9089 | 0.9948 + 0.1014i | 0.9948 - 0.1014i |
| 0.16 | 2.4712e+06 | 4.0464e-07 | 1.1624 | 0.8603 | 0.9982 + 0.0601i | 0.9982 - 0.0601i |
| 0.18 | 2.7552e+06 | 3.6307e-07 | 1.2374 | 0.8082 | 0.9971 + 0.0764i | 0.9971 - 0.0764i |
| 0.20 | 2.9533e+06 | 3.3824e-07 | 1.3246 | 0.7549 | 0.9954 + 0.0962i | 0.9954 - 0.0962i |
| 0.21 | 3.0388e+06 | 3.2940e-07 | 1.3749 | 0.7273 | 0.9944 + 0.1059i | 0.9944 - 0.1059i |

### Table 3 (p299): apoapsis group, lighter branch, mu = 0.004 (Fig. 17)

| e | lambda1 | 1/lambda1 | lambda2 | 1/lambda2 | lambda3 | 1/lambda3 |
|---|---|---|---|---|---|---|
| 0 | 9.3110e+05 | 1.0740e-06 | 1.0004 | 0.9996 | 0.9021 + 0.4316i | 0.9021 - 0.4316i |
| 0.02 | 9.3063e+05 | 1.0746e-06 | 0.9905 + 0.1378i | 0.9905 - 0.1378i | 0.9052 + 0.4250i | 0.9052 - 0.4250i |
| 0.04 | 9.2647e+05 | 1.0793e-06 | 1.0000 + 0.0060i | 1.0000 - 0.0060i | 0.9142 + 0.4053i | 0.9142 - 0.4053i |
| 0.06 | 9.0444e+05 | 1.1056e-06 | 0.9952 + 0.0977i | 0.9952 - 0.0977i | 0.9282 + 0.3721i | 0.9282 - 0.3721i |
| 0.08 | 8.1800e+05 | 1.2225e-06 | 0.9967 + 0.0808i | 0.9967 - 0.0808i | 0.9457 + 0.3250i | 0.9457 - 0.3250i |
| 0.10 | 6.4486e+05 | 1.5507e-06 | 1.0320 + 0.1017i | 0.9597 - 0.0946i | 0.9647 + 0.2635i | 0.9647 - 0.2635i |
| 0.12 | 4.3834e+05 | 2.2812e-06 | 1.0502 + 0.1743i | 0.9267 - 0.1538i | 0.9824 + 0.1870i | 0.9824 - 0.1870i |
| 0.14 | 2.3927e+05 | 4.1794e-06 | 1.0019 + 0.2994i | 0.9163 - 0.2738i | 0.9948 + 0.1014i | 0.9948 - 0.1014i |
| 0.16 | 8.1398e+04 | 1.2285e-05 | 0.5867 + 0.8098i | 0.5867 - 0.8098i | 0.9982 + 0.0601i | 0.9982 - 0.0601i |
| 0.18 | 2.8907e+03 | 3.4593e-04 | -46.5318 | -0.0215 | 0.9971 + 0.0764i | 0.9971 - 0.0764i |
| 0.20 | 5.0494e+04 | 1.9805e-05 | -2.6794 | -0.3732 | 0.9954 + 0.0962i | 0.9954 - 0.0962i |
| 0.21 | 1.4006e+05 | 7.1400e-06 | -0.1870 + 0.9824i | -0.1870 - 0.9824i | 0.9944 + 0.1059i | 0.9944 - 0.1059i |

### Table 4 (p299): apoapsis group, medium branch, mu = 0.009 (Fig. 18)

| e | lambda1 | 1/lambda1 | lambda2 | 1/lambda2 | lambda3 | 1/lambda3 |
|---|---|---|---|---|---|---|
| 0 | 1.3112e+06 | 7.6272e-07 | 1.0000 + 0.0001i | 1.0000 - 0.0001i | 0.9587 + 0.2843i | 0.9587 - 0.2843i |
| 0.02 | 1.3169e+06 | 7.5949e-07 | 1.0000 + 0.0008i | 1.0000 - 0.0008i | 0.9607 + 0.2777i | 0.9607 - 0.2777i |
| 0.04 | 1.3339e+06 | 7.4968e-07 | 1.0000 + 0.0043i | 1.0000 - 0.0043i | 0.9662 + 0.2577i | 0.9662 - 0.2577i |
| 0.06 | 1.3623e+06 | 7.3434e-07 | 0.9745 + 0.2245i | 0.9745 - 0.2245i | 0.9999 + 0.0119i | 0.9999 - 0.0119i |
| 0.08 | 1.4018e+06 | 7.1325e-07 | 0.9997 + 0.0248i | 0.9997 - 0.0248i | 0.9841 + 0.1778i | 0.9841 - 0.1778i |
| 0.10 | 1.4487e+06 | 6.9018e-07 | 0.9990 + 0.0457i | 0.9990 - 0.0457i | 0.9931 + 0.1169i | 0.9931 - 0.1169i |
| 0.12 | 1.4488e+06 | 6.9014e-07 | 1.0102 + 0.0656i | 0.9858 - 0.0640i | 1.0102 - 0.0656i | 0.9858 + 0.0640i |
| 0.14 | 1.0349e+06 | 9.6637e-07 | 1.0121 + 0.1266i | 0.9728 - 0.1217i | 1.0121 - 0.1266i | 0.9728 + 0.1217i |
| 0.16 | 6.2201e+05 | 1.6077e-06 | 1.0015 + 0.2330i | 0.9472 - 0.2204i | 1.0015 - 0.2330i | 0.9472 + 0.2204i |
| 0.18 | 2.8905e+05 | 3.4596e-06 | 0.8384 + 0.5451i | 0.8384 - 0.5451i | 0.9623 + 0.2718i | 0.9623 - 0.2718i |
| 0.20 | 6.6894e+04 | 1.4949e-05 | -0.4359 + 0.9000i | -0.4359 - 0.9000i | 0.9421 + 0.3354i | 0.9421 - 0.3354i |
| 0.21 | 1.1938e+04 | 8.3769e-05 | -19.7985 | -0.0505 | 0.9288 + 0.3705i | 0.9288 - 0.3705i |

### Table 5 (p300): apoapsis group, heavier branch, mu = 0.015 (Fig. 19)

| e | lambda1 | 1/lambda1 | lambda2 | 1/lambda2 | lambda3 | 1/lambda3 |
|---|---|---|---|---|---|---|
| 0 | 1.5966e+06 | 6.2636e-07 | 1.0002 | 0.9998 | 0.9021 + 0.4316i | 0.9021 - 0.4316i |
| 0.02 | 1.6061e+06 | 6.2266e-07 | 0.9052 + 0.4250i | 0.9052 - 0.4250i | 1.0000 + 0.0006i | 1.0000 - 0.0006i |
| 0.04 | 1.6351e+06 | 6.1147e-07 | 0.9142 + 0.4053i | 0.9142 - 0.4053i | 1.0000 + 0.0035i | 1.0000 - 0.0035i |
| 0.06 | 1.6851e+06 | 5.9339e-07 | 0.9282 + 0.3721i | 0.9282 - 0.3721i | 1.0000 + 0.0097i | 1.0000 - 0.0097i |
| 0.08 | 1.7590e+06 | 5.6862e-07 | 0.9457 + 0.3250i | 0.9457 - 0.3250i | 0.9998 + 0.0200i | 0.9998 - 0.0200i |
| 0.10 | 1.8623e+06 | 5.3682e-07 | 0.9648 + 0.2630i | 0.9648 - 0.2630i | 0.9994 + 0.0353i | 0.9994 - 0.0353i |
| 0.12 | 2.0062e+06 | 4.9841e-07 | 0.9834 + 0.1814i | 0.9834 - 0.1814i | 0.9982 + 0.0600i | 0.9982 - 0.0600i |
| 0.14 | 2.2197e+06 | 4.5037e-07 | 1.0508 + 0.0908i | 0.9446 - 0.0816i | 1.0508 - 0.0908i | 0.9446 + 0.0816i |
| 0.16 | 2.5436e+06 | 3.9332e-07 | 1.0935 + 0.0561i | 0.9121 - 0.0468i | 1.0935 - 0.0561i | 0.9121 + 0.0468i |
| 0.18 | 2.8381e+06 | 3.5188e-07 | 1.1411 + 0.0264i | 0.8759 - 0.0203i | 1.1411 - 0.0264i | 0.8759 + 0.0203i |
| 0.20 | 3.0812e+06 | 3.2454e-07 | 1.2449 | 0.8033 | 1.1417 | 0.8759 |
| 0.21 | 3.1997e+06 | 3.1249e-07 | 1.2950 | 0.7722 | 1.1489 | 0.8704 |

Sanity checks on the transcription (COMPUTED). The product lambda1 x (1/lambda1) over all 57 rows of the five tables lies between
0.9987 and 1.0010 (largest departures: Table 5 e = 0.18, 0.9987; Table 2 e = 0.20, 0.9989; Table 2 e = 0.21, 1.0010), which is the
four-digit precision the authors claim. For the complex unit pairs |lambda3|^2 at e = 0 is 0.9587^2 + 0.2843^2 = 0.99993 (Table 1) and
0.9021^2 + 0.4316^2 = 1.00006 (Table 2). Reciprocal of a complex off-circle number: 1/(1.0102 + 0.0656i) = (1.0102 - 0.0656i)/1.02479 =
0.9858 - 0.0640i, as printed in Table 4 at e = 0.12.

### Printed anomalies in the tables (INFERRED; flagged so no check is built on a possible slip)

1. The lambda3 column of Table 3 (mu = 0.004, apoapsis) is identical, to all four printed digits, to the lambda3 column of Table 2
   (mu = 0.015, periapsis) in all twelve rows (0.9021 + 0.4316i, 0.9052 + 0.4250i, 0.9142 + 0.4053i, 0.9282 + 0.3721i, 0.9457 + 0.3250i,
   0.9647 + 0.2635i, 0.9824 + 0.1870i, 0.9948 + 0.1014i, 0.9982 + 0.0601i, 0.9971 + 0.0764i, 0.9954 + 0.0962i, 0.9944 + 0.1059i). The
   lambda2 column of Table 5 (mu = 0.015, apoapsis) rows e = 0.02 to 0.12 carries the same sequence to within 0.001. Table 1 (mu = 0.009)
   has a different first row (0.9587 + 0.2843i). A unit-circle pair that is the same at e = 0 for mu = 0.004 and 0.015 but different at
   0.009 is not physically plausible; the Table 3 lambda3 column is probably a copy of the Table 2 column. Do not use the Table 3
   lambda3 column or the Table 3 e = 0 row for lambda3 as a control.
2. Column labels lambda2 and lambda3 are not tracked across rows: Table 4 swaps the roles at e = 0.06 (Table 1 e = 0.06 has
   lambda3 = 0.9745 + 0.2245i; Table 4 e = 0.06 has this as lambda2); Table 5 rows e = 0.02 to 0.12 list the near-1 pair as lambda3 and
   the other as lambda2, the reverse of the convention elsewhere. Compare tables as sets of eigenvalues, not column by column.
3. The Earth-Moon text on p300 calls the single real pair of the apoapsis orbit "lambda2 = 1.5431 x 10^6" (section 4.4).
4. Table 1 prints no mu in its caption (the caption says "illustrated in Fig. 10"); the value mu = 0.009 is from the Fig. 10 title and
   the text on p293 ("0.009 < mu^"). Table 2's caption likewise cites Fig. 11 (mu = 0.015 in the Fig. 11 title).
5. The paper writes "CRTPB" for the circular problem on pp284 and 290, "ERTBPB" and similar typographical slips; none affects a number.

## 6. Techniques applicable to the project's problems

### (a) Recipe: a periodic orbit of the elliptic problem from a circular-problem orbit, what the project has, what it lacks

Steps as the paper does them, with the status in the project (names from reading the modules; where I could not tell, I say so).

1. Choose coprime (M, N); target circular period T_C = 2 N pi / M. In the project: CR3BP periodic-orbit and family code exist
   (`search/cr3bp_periodic.py`, the 3D family tracer named in `genome/er3bp_continuation.py`'s docstring); the target-period
   bisection along a family (the paper's dichotomy of Az) was not looked for and is not known to be present.
2. Take the symmetric circular initial state X0 = [x0, 0, z0, 0, y0', 0] of period T_C. At e = 0 the pulsating frame is the rotating
   frame, so the same state integrated for T_E = N 2 pi = M T_C is an exactly periodic ER3BP orbit (M-fold cover). In the project:
   `core/er3bp.py:propagate_er3bp` with e = 0 does this. READ for the method, INFERRED for the code.
3. Phase: periapsis group f0 = 0, apoapsis group f0 = pi (or, equivalently by Eq. 6, f0 = 0 with e -> -e). In the project:
   `genome/er3bp_periodic.py:correct_er3bp_periodic` and `core/er3bp.py:propagate_er3bp` start at f = 0 only in the corrector
   (`f_span=(0.0, period_f)` hard-coded in the corrector), so the apoapsis group is reachable only via a negative e, which no code
   path has been shown to use. Lacking: an explicit f0 argument. The periapsis group is available.
4. Symmetry conditions at the two apses. Start: y = x' = z' = 0 (Eq. 16 gives X0 with y0 = x0' = z0' = 0) at f0. Half period: at
   f0 + N pi, y = x' = z' = 0 (Eq. 18). Unknown X0 = (x0, z0, y0'), three equations, three unknowns, T_E fixed (so, unlike the
   circular halo, no amplitude Az is fixed). In the project: `correct_er3bp_periodic(free_vars, residual_indices)` is generic; its
   docstring names the 3D symmetric half-period residual (IDX_Y, IDX_XDOT, IDX_ZDOT), and `free_vars=(IDX_X, IDX_Z, IDX_YDOT)` is the
   corresponding choice (INFERRED from reading the signature; not run). The continuation drivers in `genome/er3bp_continuation.py`
   are hard-wired to the planar choice, free (x0, ydot0), residual (y, xdot), seed (x0, 0, 0, 0, ydot0, 0), so they do not do 3D halos.
5. Continue in e with a step (the paper: 0.001) using the previous orbit as the initial guess. In the project:
   `continue_er3bp_family_in_e` (secant) and `continue_er3bp_family_in_e_arclength` (pseudo-arclength, walks through folds); both planar.
6. Correct each step with a multi-segment (multiple-shooting) problem: n segments, unknown nodes and epochs, continuity (Eq. 19),
   terminal residual (Eq. 20), tolerance 1e-10 on discontinuities, n between M/2 and 2M. In the project: single shooting only for
   the elliptic problem. Multiple-shooting correctors exist for the circular problem (`genome/multi_shooting.py`,
   `search/cr3bp_multiple_shooting.py`, both stated to be circular-problem) but no elliptic one was found. Lacking: an ER3BP
   multi-segment corrector (the segment propagator `propagate_er3bp(with_stm=True)` that it would reuse does exist).
7. Stability: multiply segment STMs over f0 to f0 + T_E (Eq. 22) in at least four segments, then eigenvalues and indices (Eq. 23).
   In the project: `search/er3bp_floquet.py:er3bp_monodromy` is a single full-period STM from f = 0 (the segmented product the
   authors prefer for accuracy is not implemented); `floquet_classify` returns only "stable", "unstable" or "marginal" from the largest
   |lambda| and a unit-circle flag. Lacking: the stability indices of Eq. 23, the negative-real and complex-off-circle classification,
   the eigenvalue-type bookkeeping of Fig. 20, and an f0 argument. `search/er3bp_periodic.py:monodromy_eigenstructure` assumes a
   Lagrange-orbit structure (one saddle, one unit-circle centre pair, raises if no complex unit-circle pair is found), which several
   of this paper's orbits violate (Table 5 e = 0.20 and 0.21 have three real pairs and no complex eigenvalue, so it would raise;
   an off-circle complex quartet such as Table 4 e = 0.12 can be mistaken for a centre pair because its modulus is within 0.5 of 1;
   INFERRED from reading the function, not run).

Plan for the missing pieces, in order of value: (i) add an f0 argument (or document the negative-e route) and test it with the
identity propagate(e, f0 = pi) = propagate(-e, f0 = 0); (ii) a segmented monodromy per Eq. 22; (iii) the index of Eq. 23; (iv) an
ER3BP multi-segment corrector for N >= 2 and M >= 4, which is where single shooting loses accuracy (see (e)).

### (b) Stability of the project's Earth-Moon cycler families carried to the elliptic problem

Method: for each circular cycler of period T_C = 2 N pi / M (or tuned to it by the family parameter), build both counterparts
(periapsis and apoapsis starts) and continue in e to 0.0554, computing the monodromy over T_E = 2 N pi by the segmented product, and tabulate (lambda1, lambda2,
lambda3) and the types. Reading the results with this paper's yardstick: (i) at e = 0 the multipliers are the M-th powers of the
single-period ones: COMPUTED, lambda1 = 1.3112e6 (mu = 0.009) means a single-period multiplier of exp(ln(1.3112e6)/5) = 16.73;
mu = 0.015 gives 17.40; (ii) the project's cyclers are not halo orbits; INFERRED that the Earth-Moon cycler families (Casoliva et al. 2008, 2010) are
long planar orbits whose multipliers over several revolutions can be as large as the halo ones, in which case the paper's warning
applies: stability numbers should come from a segmented product, and orbits with N >= 3 would need the multi-segment corrector.
Their actual multipliers have not been computed here. (iii) The
index and the type sequence say what kind of persistence to expect: three real pairs (heavier apoapsis branch, Table 5 e = 0.20)
means no centre direction at all; one unit-circle pair (all periapsis rows) means a centre direction, which is what a station-keeping
or quasi-cycler design relies on.

### (c) Two counterparts for every circular orbit whose period is commensurate with the primaries' period

READ: for each (M, N) there are two ME-Halo groups (periapsis and apoapsis for M odd, left and right for M even), with
different shapes (Fig. 4, 5) and different stability (Tables 1 to 5; the paper's "lighter" and "heavier" branches differ per group:
two branches for the periapsis group, three for the apoapsis group). At e = 0 they are the same orbit (INFERRED from Tables 1 and 4,
and 2 and 5: identical e = 0 rows for lambda1 to 4 digits and for the unit pair, apart from the Table 3 problem above).

New observation, COMPUTED from the tables and INFERRED as a rule. The trivial pair at +1 splits at first order in e into a REAL pair
in the periapsis group and a UNIT-CIRCLE pair in the apoapsis group, with nearly the same magnitude delta of departure from 1:

| mu | e | periapsis, lambda2 - 1 (Table 1 or 2) | apoapsis, imaginary part of the near-1 pair (Table 4 or 5) |
|---|---|---|---|
| 0.009 | 0.02 | 0.0008 | 0.0008 |
| 0.009 | 0.04 | 0.0043 | 0.0043 |
| 0.009 | 0.06 | 0.0119 | 0.0119 (printed as lambda3 = 0.9999 + 0.0119i) |
| 0.009 | 0.08 | 0.0250 | 0.0248 |
| 0.009 | 0.10 | 0.0452 | 0.0457 |
| 0.015 | 0.02 | 0.0006 | 0.0006 |
| 0.015 | 0.04 | 0.0035 | 0.0035 |
| 0.015 | 0.06 | 0.0098 | 0.0097 |
| 0.015 | 0.08 | 0.0202 | 0.0200 |
| 0.015 | 0.10 | 0.0359 | 0.0353 |
| 0.015 | 0.12 | 0.0594 | 0.0600 |

The same unit-circle pair (0.9745 +- 0.2245i at mu = 0.009, e = 0.06 and so on) appears in both groups in the other slot. COMPUTED:
delta / e^(5/2) is 13.4, 13.5, 13.8, 14.3 (mu = 0.009, e = 0.04 to 0.10) and 10.9, 11.1, 11.2, 11.4 (mu = 0.015), nearly constant.
So delta is about C e^(5/2) with C of order 11 to 14, hence (lambda - 1)^2 is about +/- C^2 e^5, a function of e of odd power, whose
sign flips under e -> -e. INFERRED: this is exactly the sign flip between the two groups (periapsis start with +e is apoapsis start
with -e, Eq. 6 of the 2017 paper and section 1 above), so one group gets a saddle pair and the other a centre pair from the same
trivial pair. Park and Howell 2024 (Park-Howell digest, section on counterparts) say the same in words: the trivial pair
"bifurcates into a saddle and a centre pair for each counterpart" and counterparts are linked by positive and negative e for odd p.
CONJECTURE, untested: the exponent 5/2 equals M/2 for M = 5, which would give (lambda - 1)^2 proportional to e^M, a sign flip
between counterparts for odd M, and, for M = 2, a splitting linear in e that crosses the unit circle at e = 0, which would be the
"bifurcation at e = 0" that Campagnola found for left and right M2N1. A test: compute both counterparts of a second resonance (M7N3
or M3N1) and fit the exponent.

Consequences for the project: (1) any statement about the elliptic persistence or stability of a resonant circular cycler must
name which counterpart it is about; one counterpart can be saddle-type and the other centre-type in the near-1 pair while sharing
every other multiplier to a percent at small e; (2) a catalogue row or a test for a resonant orbit in the elliptic model should record
f0 (or the sign of e); (3) both counterparts must be run before a "survives" or "fails" verdict is given. The elliptic discovery
result `#435` ("survives with no bifurcation") used planar seeds from f = 0 only, i.e. one counterpart per orbit (INFERRED from the
continuation code above).

### (d) Relation to the folds in eccentricity that Park and Howell use as a difficulty predictor

READ in this paper: the continuation fails at a critical e^(mu) for mu below mu^ in the periapsis group (e^ about 0.143 at mu =
0.012; about 0.16 at mu = 0.009, INFERRED from where Table 1 and Fig. 10 end), never for the apoapsis group (all reach 0.21); the
authors "inferred that either there is a bifurcation or the characteristic curve is not monotonic". That is the fold in e that Park
and Howell (2024) later use: a fold is a turn in e, where the monodromy has an extra unity eigenvalue pair, and "the proliferation of
fold bifurcations" before e = 0.055 predicts transition difficulty. Linking the two (INFERRED):
- this paper's failure of natural-parameter continuation at e^(mu) is the fold signature; Peng, Bai & Xu 2017 confirm that tangential
  continuation is needed to pass such points (2017 digest, section 5); the project's arclength continuator is the right tool;
- the lighter periapsis branch (mu = 0.009) has, just before its end, a cluster of three eigenvalue collisions within about 2e-5 in
  e (Fig. 12, e about 0.1582 to 0.1584, read from the figure) and a sign change of the large real pair. A fold in Park and Howell's
  sense needs an extra eigenvalue pair at +1; the printed tables do not show one (Table 1 has the near-1 pair at 1.2073 and 0.8283
  at e = 0.14, then complex 0.8705 +- 0.4922i at e = 0.16). So the identification of the failure with a fold rests on the authors'
  "either a bifurcation or a non-monotonic characteristic curve" inference, not on an observed unity eigenvalue. (Also INFERRED: a
  real pair cannot pass through zero because the determinant is 1, so the sign change of lambda1 between e = 0.14 and 0.16 must pass
  through the collisions of Fig. 12.)
- the apoapsis group's absence of failure up to 0.21, with a near-1 pair that stays on the unit circle until e about 0.10 to 0.12 (Tables 4 and 5), is
  consistent with Park and Howell's observation that the two counterparts can have different fold behaviour; for the project this
  means that the fold-count predictor of ephemeris-model difficulty should be evaluated per counterpart, not per circular orbit.

### (e) The Neelakantan and Ramanan orbit the project could not reproduce (M4N2 Lyapunov, two primary revolutions)

The row is Table 8 of Neelakantan & Ramanan 2022, "M4N2 Lyapunov" (x0 = 0.804504659626012, ydot0 = 0.31826866733409, period 4 pi,
mu = 0.0122, e = 0.0554); the project's test records closure errors of 4.5 and 2.0 from periapsis and apoapsis, while its sibling
M2N1 (x0 = 0.804125156956177, ydot0 = 0.31182413982453, period 2 pi) closes to 1e-9 (`tests/core/test_er3bp_neelakantan_2022.py`).
This paper does not treat Lyapunov orbits, M4N2 or planar orbits at all, and prints no initial condition, so nothing here can reproduce
the row. What it does bear on:

1. Not coprime. M = 4, N = 2 has T_C = 2 N pi / M = pi, the same circular period as M2N1; the circular orbit of period pi is one
   orbit, so in the e = 0 limit M4N2 is the 2-fold cover of the M2N1 circular Lyapunov orbit. (INFERRED.) The paper's classification
   (section 2.2) is for coprime (M, N); the covering orbit is a solution for every e (the M2N1 orbit repeated twice with the same
   x0), and a different 4-pi orbit exists only after a period-doubling-type bifurcation. The printed x0 differs from the M2N1 value
   by 3.8e-4 (COMPUTED: 0.804504659626012 - 0.804125156956177 = 3.795e-4), so the printed row would be that distinct orbit, not the
   cover. A check that costs nothing: compute the circular (e = 0) Lyapunov orbit of period pi at mu = 0.0122 with the project's
   circular corrector, then see which of the two printed x0 values lies nearer (INFERRED: both should be within O(e) of it).
2. Digits against amplification. The paper's own warning (p290) is that over longer intervals the largest multiplier grows
   exponentially and the largest and smallest differ "more than 20 orders". Crude COMPUTED estimate: the M2N1 row closes to 1.1e-9
   over 2 pi with a state uncertainty e0 of 1e-15 (printed digits) to 1e-13 (integration tolerance used by the test), so the
   amplification over 2 pi is about 1e4 to 1e6; the M4N2 orbit is over 4 pi with the same circular orbit underneath, so the
   amplification is about the square, 1e8 to 1e12, giving an expected closure error of about 1e-5 to 1e-3. The project observed 4.5 and
   2.0. So the printed digits alone do not obviously explain the failure (three to five orders short), but the estimate assumes the
   covering orbit's growth and a linear model; it supports "a different orbit or a slip in the row" more than "too few digits".
   (INFERRED; the actual monodromy of the M2N1 printed state would sharpen it.)
3. Which side. For M even the paper has a Left and a Right group, both started at f0 = 0, distinguished by the side of the orbit at
   the perpendicular crossing. The L1 point of the Earth-Moon system with mu = 0.0122 is at x about 0.837 (COMPUTED, Hill
   approximation, not from the paper), so the printed x0 = 0.8045 is on the Earth side (left) of the libration region. The project
   tested it from periapsis and from apoapsis; the Right-group counterpart (a crossing on the Moon side) was not tried, since the
   table prints one state. INFERRED possibility only: if the orbit of the two-impulse design was a Right-group one and the Left
   side's number was printed, no closure would result. A search for the nearest closed orbit from both sides would show it.
4. The decisive test, in the paper's own terms: re-close the printed state with the multi-segment method (n between M/2 and 2M, so
   2 to 8 segments for M = 4) at e = 0.0554 and compare. If a closed orbit exists within about 1e-3 of the printed state, the printed
   digits are not the issue and the single-shooting closure test is too sensitive; if no closed orbit exists nearby, the row has a
   slip or is a different group or orbit. This needs the ER3BP multi-segment corrector that the project lacks (section (a), step 6),
   and is the main reason to build it. Not run here (no computation was requested).

## 7. Positive controls

Direct controls (an orbit state to close in `core/er3bp.py`): NONE. The paper prints no initial condition, no period beyond the
rational formula T_E = 2 N pi (M5N2: 4 pi = 12.5664), and no orbit state. Figures 6 and 14 give only axis ranges for the family in the
(x0, z0, y0') subspace (x0 about 0.8 to 0.95, z0 about 0.08 to 0.25, y0' about 0.1 to 0.35), which cannot serve as a control.

Eigenvalue controls (the paper's printed numbers as the expected side, which the project can test; none has a test yet):

| Control | Setup in the project | Expected (printed) | Tolerance |
|---|---|---|---|
| C1. Circular limit, mu = 0.009 | Find the L1 north halo of CRTBP period T_C = 4 pi / 5 = 2.51327 at mu = 0.009 (bisect on Az, x-z plane crossing, `search/cr3bp_periodic.py`), single-period monodromy M1, M1^5 | eigenvalues 1.3112e6, 7.6258e-07 (Table 1; 7.6272e-07 in Table 4), the pair at 1 (1.0002, 0.9998 printed; 1 within the authors' error), 0.9587 +- 0.2843i | 4 digits on the large pair, 3 on the unit pair (the authors' stated reliability) |
| C2. Circular limit, mu = 0.015 | Same at mu = 0.015 | 1.5966e6, 6.2662e-07 (Table 2; 6.2636e-07 in Table 5), pair at 1, 0.9021 +- 0.4316i | same |
| C3. Periapsis near-1 pair, mu = 0.015 | With a corrector for the 3D periapsis orbit at e = 0.06, T_E = 4 pi, f0 = 0 | lambda1 1.6845e6, lambda2 1.0098 and 1/lambda2 0.9903, unit pair 0.9282 +- 0.3721i (Table 2) | 3 digits on lambda1, 2 on lambda2 - 1 |
| C4. Earth-Moon, periapsis | mu = 0.0122, e = 0.0554, f0 = 0, T_E = 4 pi; INFERRED candidate state the N&R M5N2 halo [x0, 0, z0, 0, ydot0, 0] = [0.851666641652152, 0, 0.183285539178136, 0, 0.25828972225268, 0] (closes from periapsis to 2.9e-7 in the project's test; the paper does not say it is the same orbit) | lambda1 about 1.5427e6 (so 1/lambda1 = 6.4821e-7), lambda2 = 1.0086 (1/lambda2 = 0.99147), unit pair (lambda3) | 3 digits on lambda1, lambda2 - 1 to a few percent |
| C5. Earth-Moon, apoapsis | Same, f0 = pi (or e = -0.0554 at f0 = 0) | one real pair about 1.5431e6, two unit-circle pairs | 3 digits |
| C6. Splitting law | Compute the near-1 pair at e = 0.02, 0.04, 0.06, 0.08, 0.10 for both counterparts at mu = 0.015 | |lambda - 1| about 0.0006, 0.0035, 0.0098, 0.0202, 0.0359 (periapsis, real) and imaginary 0.0006, 0.0035, 0.0097, 0.0200, 0.0353 (apoapsis) | 5 to 10 percent |

For C4: the circular M5N2 halo state and periapsis f0 = 0 recipe in `core/er3bp.py` is: state0 = [x0, 0, z0, 0, ydot0, 0]; call
`propagate_er3bp(state0, (0, 4 pi), ER3BPSystem(mu=0.0122, e=0.0554, ...), with_stm=True)`, or better, four calls over (0, pi),
(pi, 2 pi), (2 pi, 3 pi), (3 pi, 4 pi), multiplying the four STMs in order (Eq. 22; a single 4-pi call has relative error about 1e-6
from the 1.5e6 multiplier at the project's default tolerance, INFERRED). Expected closure of the state is the project's measured
2.9e-7. Failing C4 would not disprove
anything on its own, since the paper does not say that the N&R state is its own orbit and the Earth-Moon halos have multiple solutions
(N&R p2); it would be reported, not loosened. C1, C2 and C6 do not depend on that identification.

Controls not to use: the Table 3 lambda3 column (anomaly 1 above).

## 8. What the paper leaves open (quoted), and its references

Left open, verbatim:
- "The authors suppose a more profound insight into the stability of these orbit requires more advanced mathematic methods dealing
  with the nonlinearity." (p290)
- "it can be inferred that either there is a bifurcation or the characteristic curve is not monotonic with respect to the continuation
  parameter, which causes the failure of continuation using our algorithm." (p292)
- "We suggest this phenomenon is caused by the interaction between mu and e, which cannot likely be explained by numerical studies."
  (p295)
- "In the Periapsis Group, a continuation barrier arises for small mass ratio, which seems to be caused by the change of eigenvalues
  from positive to negative but still needs more analytical studies in the future." (p302)
- "However, a quantified result can hardly be drawn from the present numerical study." (p301)
- The separator mu^ for the apoapsis group lies in (0.010, 0.011), "But it was not directly detected by the discrete grids." (p296)
- Only N = 2 (M5N2) is studied, L1 and the north branch; no L2, no larger N, no other M ("The authors suppose ..." above).
- Not addressed at all: any planar (Lyapunov) orbit, the stability relation to Floquet theory of the elliptic problem beyond the
  index, the dependence on the sign convention of e, and any initial condition.

References printed (journal pp302-303). Held status uses `docs/notes/CORPUS_INDEX.md` and the file names in the private corpus; the
check was by index text and file name, not by searching the body of every held file. "Held" means the same work; a different work by
the same author is said so.

| Reference as printed | Held? |
|---|---|
| Antoniadou, K.I., Voyatzis, G.: 2/1 Resonant periodic orbits in three dimensional planetary systems. Celest. Mech. Dyn. Astron. 115(2), 161-184 (2013) | held (digested 2026-07-12) |
| Barden, B.T., Howell, K.C., Lo, M.W.: Application of dynamical systems theory to trajectory design for a libration point mission. 268-281 (1996) | not held |
| Belbruno, E.A., Gidea, M., Topputo, F.: Geometry of weak stability boundaries. Qual. Theory Dyn. Syst. (2012) | not held (Belbruno 2004 textbook is a different work, held) |
| Bittanti, S., Colaneri, P.: Periodic Systems, Communications and Control Engineering vol. 36, Springer, London (2009) | not held |
| Broucke, R.A.: Stability of periodic orbits in the elliptic, restricted three-body problem. AIAA J. 7(6), 1003-1009 (1969) | not held (Broucke 1968 JPL TR 32-1168, a different work, is held) |
| Campagnola, S.: New Techniques in Astrodynamics for Moon Systems Exploration. PhD dissertation, University of Southern California (2010) | not held |
| Campagnola, S., Lo, M.W., Newton, P.: Subregions of motion and elliptic halo orbits in the elliptic restricted three-body problem. 18th AAS/AIAA Spaceflight Mechanics Meeting, Galveston (2008) | not held |
| Farquhar, R.W., Kamel, A.A.: Quasi-periodic orbits about the translunar libration point. Celest. Mech. 7(4), 458-473 (1973) | not held |
| Gomez, G., Koon, W.S., Lo, M.W., Marsden, J.E., Masdemont, J.J., Ross, S.D.: Connecting orbits and invariant manifolds in the spatial restricted three-body problem. Nonlinearity 17(5), 1571-1606 (2004) | held |
| Gurfil, P., Kasdin, N.J.: Niching genetic algorithms-based characterization of geocentric orbits in the 3D elliptic restricted three-body problem. Comput. Methods Appl. Mech. Eng. 191(49-50), 5683-5706 (2002) | held |
| Gurfil, P., Meltzer, D.: Semi-analytical method for calculating the elliptic restricted three-body problem monodromy matrix. J. Guid. Control Dyn. 30(1), 266-271 (2007) | not held (Gurfil ed., Modern Astrodynamics, is a different work, held) |
| Heppenheimer, T.A.: Out-of-plane motion about libration points: nonlinearity and eccentricity effects. Celest. Mech. 7(2), 177-194 (1973) | not held |
| Hiday, L.A., Howell, K.C.: Transfers between libration-point orbits in the elliptic restricted problem. Celest. Mech. Dyn. Astron. 58(4), 317-337 (1994) | not held |
| Hou, X.Y., Liu, L.: On motions around the collinear libration points in the elliptic restricted three-body problem. Mon. Not. R. Astron. Soc. 415(4), 3552-3560 (2011) | not held |
| Howell, K.C., Pernicka, H.J.: Numerical determination of Lissajous trajectories in the restricted three-body problem. Celest. Mech. 41(1-4), 107-124 (1987) | held (file dated 1988; digest not rechecked here) |
| Hyeraci, N., Topputo, F.: Method to design ballistic capture in the elliptic restricted three-body problem. J. Guid. Control Dyn. 33(6), 1814-1823 (2010) | not held |
| Hyeraci, N., Topputo, F.: The role of true anomaly in ballistic capture. Celest. Mech. Dyn. Astron. 116(2), 175-193 (2013) | not held |
| Ichtiaroglou, S.: Elliptic Hill's problem: the continuation of periodic orbits. Astron. Astrophys. 92, 139-141 (1980) | not held |
| Ichtiaroglou, S., Michalodimitrakis, M.: Three-body problem: the existence of families of three-dimensional periodic orbits which bifurcate from planar periodic orbits. Astron. Astrophys. 81, 30-32 (1980) | not held |
| Koon, W.S., Lo, M.W., Marsden, J.E., Ross, S.D.: Shoot the moon. Advances in Astronautical Sciences 105, AAS, San Diego, 1017-1030 (2000) | held |
| Koon, W.S., Lo, M.W., Marsden, J.E., Ross, S.D.: Dynamical Systems, the Three-Body Problem and Space Mission Design. Marsden Books (2011) | held in the 2006 printing of the same book |
| Lei, H., Xu, B., Hou, X., Sun, Y.: High-order solutions of invariant manifolds associated with libration point orbits in the elliptic restricted three-body system. Celest. Mech. Dyn. Astron. 117(4), 349-384 (2013) | not held |
| Mahajan, B.: Libration point orbits near small bodies in the elliptic restricted three-body problem. Masters Theses, Paper 7200, Missouri S&T (2013) | not held |
| Mahajan, B., Pernicka, H.J.: Halo orbits near small bodies in the elliptic restricted problem. 1-9 (2012), doi:10.2514/6.2012-4876 | not held |
| Martin, C., Conway, B.A., Ibanez, P., Offin, D.: Optimal low-thrust trajectories to the interior Earth-Moon Lagrange point. In: Space Manifold Dynamics, 161-184, Springer (2010) | not held |
| Meyer, K.R., Hall, G.R., Offin, D.: Introduction to Hamiltonian Dynamical Systems and the N-body Problem, 2nd edn., Springer (2009) | not held |
| Moulton, F.R.: Periodic Orbits. Carnegie Institution of Washington (1920) | not held |
| Parker, J.S., Anderson, R.L.: Low-Energy Lunar Trajectory Design, JPL Deep-Space Communications and Navigation Series, Wiley (2014) | not held (Parker 2007 thesis, a different work, is held) |
| Pernicka, H.J.: The numerical determination of nominal libration point trajectories and development of a station-keeping strategy. Purdue University (1990) | not held |
| Qi, Y., Xu, S.: Lunar capture in the planar restricted three-body problem. Celest. Mech. Dyn. Astron. 120(4), 401-422 (2014) | not held |
| Qi, Y., Xu, S., Qi, R.: Gravitational lunar capture based on bicircular model in restricted four body problem. Celest. Mech. Dyn. Astron. 120(1), 1-17 (2014a) | not held |
| Qi, Y., Xu, S., Qi, R.: Study of the gravitational capture at Mercury in the elliptic restricted three-body problem. Proc. 24th ISSFD (2014b) | not held |
| Richardson, D.L.: Analytic construction of periodic orbits about the collinear points. Celest. Mech. 22(3), 241-253 (1980) | held (digested 2026-07-12) |
| Russell, R.P.: Survey of spacecraft trajectory design in strongly perturbed environments. J. Guid. Control Dyn. 35(3), 705-720 (2012) | held |
| Sarris, E.: Families of symmetric-periodic orbits in the elliptic three-dimensional restricted three-body problem. Astrophys. Space Sci. 162(1), 107-122 (1989) | not held |
| Szebehely, V.G.: Theory of Orbits: The Restricted Problem of Three Bodies. Academic Press, New York (1967) | held |
| Tarrago, P.I.: Study and assessment of low-energy Earth-Moon transfer trajectories. Universite de Liege (2007) | not held |
| Wiggins, S.: Introduction to Applied Nonlinear Dynamical Systems and Chaos, 2nd edn., Springer (2003) | not held |

Highest value to acquire for the project (INFERRED): Campagnola, Lo & Newton 2008 and Campagnola 2010 (the origin of the elliptic
halo orbit classification and of the M2N1 bifurcation at e = 0 stated on p301), Broucke 1969 (the stability index and the complex
instability), Sarris 1989 (the symmetries and families in the spatial elliptic problem), and Gurfil & Meltzer 2007 (a semi-analytic
monodromy matrix for the elliptic problem, which would cross-check the segmented product).

Papers in the project that cite this one, and what they rely on it for: Peng, Bai & Xu 2017 (its reference [19], called the authors'
own earlier paper, "the continuation failures in our early studies" it resolves; the M:N constraint, group names, f0 and -e equivalence
carry over); Neelakantan & Ramanan 2022 (the ME-Halo concept and "constructed MR halo orbits in the Earth-Moon system using the halo
orbit conditions in CRTBP as the initial guess. They divided the MR halo orbit into multiple segments and employed numerical
continuation on eccentricity to generate the design", which is Section 3.3 above; their manifold-transfer numbers are from Peng & Xu
2015b, a different paper); Park & Howell 2024 (cited once for the computation of the two counterparts, section 3 here); Singh, Park &
Howell 2026 (reference list entry only).
