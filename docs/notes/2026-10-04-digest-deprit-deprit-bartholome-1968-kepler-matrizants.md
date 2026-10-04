# Digest: Deprit & Deprit-Bartholome 1968, "The matrizants of the keplerian motions (the two-dimensional case)"

Evidence tags: READ (journal page) is what the printed page says, every equation and table entry below read on a 150 dpi page image (the PDF is a Persee scan; its text layer was not used for any formula). COMPUTED is my own arithmetic or numerical check of 2026-10-05 (scratch only, no test files written). INFERRED is my reading across sources.

Reference: Andre Deprit and Andree Deprit-Bartholome (Boeing Scientific Research Laboratories, Seattle), Bulletin astronomique 3(3):315-339, 1968, DOI 10.3406/bastr.1968.14544 (27 PDF pages: page 1 Persee cover, page 2 multilingual abstract, PDF page n = journal page n + 313 from page 3 on). Filed in the private paper corpus as `deprit-deprit-bartholome-1968-matrizants-keplerian-motions-two-dimensional-bull-astron-3-315-doi-10.3406-bastr.1968.14544.pdf`.

Citation note: Perko 1981b (digest `docs/notes/2026-10-04-digest-perko-1976b-1981b-error-estimates-first-second-species-bifurcation.md`) cites this paper as "Deprit & Deprit-Bartholome 1969, Bull. Astron. 3:315-339". The printed volume is tome 3, fascicule 3, 1968 (Persee), so the year is 1968; the 1969 in Perko is a citation slip or a later issue date (INFERRED which). Same pages, same title, same journal: the same paper.

## 1. What it is

A closed-form construction of the matrizant (state transition matrix) R(t; t0) of the planar two-body variational equations, for elliptic and hyperbolic motion, in any inertial Cartesian frame, by Jacobi's dual last-multiplier theorem. The variational equations are Hamiltonian (so the multiplier is 1), three independent integrals (angular momentum G and the Laplace vector P, Q) supply three adjoint solutions, and one quadrature supplies the fourth. Time is the running variable, so Kepler's equation is never solved: the matrizant is a function of the state (x, y, X, Y) and the elapsed time t only, with no transcendental functions. Excluded by the construction (pp.316, 319, 322, 324): rectilinear motion (G = 0, needed so the three gradients are independent), circular motion in the derivation step (P^2 + Q^2 = 0 divides the quadratures (25)), and parabolic motion (H = 0 divides (32c)). The final tables contain no P^2 + Q^2 denominator; the circular case was tested numerically (section 6) and the tables hold there too.

The paper then derives the same matrizant in the orbital (radial/transverse, rotating) frame, and in Hill's intrinsic (tangential/normal) frame, by constant-structure rotations K and L. No worked numerical example, no table of numbers, no figure: the paper is entirely formulas.

## 2. Variables and frame (READ pp.317-319, 329, 333)

- Inertial Cartesian frame Oxy centred on the attracting body; the phase space is configuration (x, y) times momentum (X, Y), where (X, Y) is the velocity (unit mass). Position (x, y) and velocity (X, Y), NOT the usual (x, y, xdot, ydot) labelling: X = xdot, Y = ydot. mu is the gravitational parameter (strictly positive), r = (x^2 + y^2)^(1/2), V = (X^2 + Y^2)^(1/2).
- Hamiltonian (1): H = (1/2)(X^2 + Y^2) - mu/r. Orbital integrals: energy H (3); angular momentum G = xY - yX (4); Laplace vector P = GY - mu x/r (5), Q = -GX - mu y/r (6); identity (7) P^2 + Q^2 = mu^2 + 2 H G^2. Four identities (8): Px + Qy = G^2 - mu r; Py - Qx = G r rdot; PX + QY = -mu rdot; PY - QX = G (V^2 - mu/r).
- Displacement or variation u = (u, w, U, W) is a solution of the variational equations (10): udot = U, wdot = W, Udot = (mu/r^3)[(3x^2/r^2 - 1) u + 3(xy/r^2) w] with the as printed on p.318 (coefficient +mu/r^3, consistent with the second derivative of -mu/r), Wdot = (mu/r^3)[3(xy/r^2) u + (3y^2/r^2 - 1) w]. Variation vector ordering throughout is (u, w, U, W) = (dx, dy, dX, dY): positions first, then velocities. G is assumed strictly positive after orienting the plane (p.319); the numerical test (section 6) shows the tables also hold for G < 0 (retrograde).
- Time t: the formulas contain the time t explicitly (secular terms 3 mu t and 3 X t). The time origin is arbitrary: the product A(t) B(t0) is unchanged by a shift of t in both factors (COMPUTED: t0 = 0 and t0 = 3 give identical residuals).
- Matrizant R(t; t0) = A(t) B(t0) (42), B(t) = A(t)^-1; the factor A(t) has as columns four independent solutions of (10), B(t) is its inverse, so R(t0; t0) = I.
- Paper equation numbering has a duplicate: "(17)" is used for both the three variational integrals (p.320, top) and the adjoint equations (p.320, section 2); the reader should treat the latter as (17') (my label). Perko 1981b's "eq. (50)" is Perko's own equation, not equation (50) of this paper (the latter is the orbital-frame Hamiltonian system, p.331).

## 3. The matrizant, entry by entry (READ pp.328, Tables I and II)

Notation as printed: a_ij, b_ij with i the row, j the column; here and below e1 = X r^2 - G y + 3 P t and e2 = Y r^2 + G x + 3 Q t. All quantities are evaluated at the time t named (state at t for A, state at t0 and the time t0 for B).

### 3.1 Table I, the factor A(t) (columns are four solutions of the variational equations)

| | column 1 | column 2 | column 3 | column 4 |
|---|---|---|---|---|
| row 1 (u, dx) | a11 = 2x - 3Xt | a12 = -y | a13 = -yY | a14 = -xY + 2yX |
| row 2 (w, dy) | a21 = 2y - 3Yt | a22 = x | a23 = -yX + 2xY | a24 = -xX |
| row 3 (U, dX) | a31 = -X + 3 mu (x/r^3) t | a32 = -Y | a33 = -Y^2 + mu y^2/r^3 | a34 = XY - mu xy/r^3 |
| row 4 (W, dY) | a41 = -Y + 3 mu (y/r^3) t | a42 = X | a43 = XY - mu xy/r^3 | a44 = -X^2 + mu x^2/r^3 |

Columns 2, 3, 4 are the isoenergetic, bounded-in-t solutions (J applied to the gradients of G, P, Q), and column 1 is the secular solution (the one carrying t, from the quadrature). A = J A*, with J as printed on p.327.

### 3.2 Table II, the factor B(t) = A(t)^-1

| | column 1 | column 2 | column 3 | column 4 |
|---|---|---|---|---|
| row 1 | b11 = -(mu/2H) x/r^3 | b12 = -(mu/2H) y/r^3 | b13 = -X/(2H) | b14 = -Y/(2H) |
| row 2 | b21 = (1/G)(-X + 3 mu (x/r^3) t) | b22 = (1/G)(-Y + 3 mu (y/r^3) t) | b23 = -(1/G)(2x - 3Xt) | b24 = -(1/G)(2y - 3Yt) |
| row 3 | b31 = -(mu/2HG^2)(x/r^3) e1 | b32 = -(mu/2HG^2)(y/r^3) e1 + 1/G | b33 = -(X/2HG^2) e1 + x^2/G^2 | b34 = -(Y/2HG^2) e1 + xy/G^2 |
| row 4 | b41 = -(mu/2HG^2)(x/r^3) e2 - 1/G | b42 = -(mu/2HG^2)(y/r^3) e2 | b43 = -(X/2HG^2) e2 + xy/G^2 | b44 = -(Y/2HG^2) e2 + y^2/G^2 |

PRINTED SLIP in b44 (READ p.328, flagged): the page prints the factor in b44 as (Yr^2 - Gy + 3Qt), the same inner combination as b34 with the sign of the G term and the Q term wrong for row 4. The correct inner factor is e2 = (Yr^2 + Gx + 3Qt), as in b41, b42, b43 and the E-components of p.326 (E123 etc. use "Yr^2 + Gx + 3Qt"). COMPUTED: with the printed b44 the product A B differs from the finite-difference Jacobian by 2.4 to 24 (order one); with e2 the residual is 4e-10 relative (section 6). The tables below and the checks use e2.

### 3.3 Intermediate results usable as checks (READ pp.318-327)

- Gradients (12): G_x = Y, G_y = -X, G_X = -y, G_Y = x; P_x = Y^2 - mu y^2/r^3, P_y = -XY + mu xy/r^3, P_X = -yY, P_Y = 2xY - yX; Q_x = -XY + mu xy/r^3, Q_y = X^2 - mu x^2/r^3, Q_X = 2yX - xY, Q_Y = -xX. (Columns 2 to 4 of Table I are J times these.)
- The three-vector grad G ^ grad P ^ grad Q (15): D123 = -mu G^2 y/r^3, D124 = mu G^2 x/r^3, D134 = G^2 Y, D234 = -G^2 X; it vanishes iff G = 0 (rectilinear). Two-vector components A_ij, B_ij, C_ij are printed in (13) and (14).
- Secular adjoint solution (38): S = (X - 3 mu (x/r^3) t, Y - 3 mu (y/r^3) t, 2x - 3Xt, 2y - 3Yt); det A*(t) = -2 H G^2 (p.326; the column vectors are S, grad G, grad P, grad Q in that order).
- The constants: delta = K1 G r^3/x (24a); alpha = K2 - K1 (3 mu t - r^3 X/x) (24b); beta = K3 + K1 (2 mu/G) integral of x dt (24c); gamma = K4 + K1 (2 mu/G) integral of y dt + K1 r^3/x (24d); K1' with 2 H K1' = -mu G K1.
- Quadratures: (31a) integral dt/r = (1/mu)(r rdot - 2 H t); (33) integral r dt = (r rdot/(4 mu H))(G^2 + mu r) - (1/2)(3 mu/(2H) + G^2/mu) t; (35a) integral x dt = (3P/4H) t + (1/4H)(X r^2 - G y); (35b) integral y dt = (3Q/2H) t + (1/4H)(Y r^2 + G x). The authors say (33), (35a), (35b) "appear to be new in the literature of the problem of two bodies". Note (35b) is printed with 3Q/(2H) where (35a) has 3P/(4H): a likely slip (the secular terms should be symmetric, 3P/4H and 3Q/4H), consistent with (34b) which has 3Q/4H; flagged, not used (the tables carry the checked forms; the table entries b-column 2 and rows 3 to 4 use 3Pt and 3Qt without an extra factor).
- The paper states in closed form (p.325): the integral of energy along a variation, its use to identify isoenergetic variations.

## 4. Orbital and intrinsic frames (READ pp.329-337)

### 4.1 Orbital frame (radial and transverse), Tables III and IV

Axes: Ox' from the centre to the particle, Oy' 90 degrees ahead in the sense of motion. theta is the azimuth, x = r cos theta, y = r sin theta; in this frame x' = r, y' = 0, X' = rdot, Y' = r thetadot = G/r (44). Displacement components (45) u = u' cos theta + w' sin theta, w = -u' sin theta + w' cos theta, and the same for (U, W): u' radial position, w' transverse position, U' radial velocity, W' transverse velocity (the canonical, not an angular, transverse velocity: the same rotation applied to velocities of the inertial frame, not the velocity relative to the rotating frame). Variational Hamiltonian (49) and equations (50) in the rotating frame: u'dot = U' + (G/r^2) w', w'dot = W' - (G/r^2) u', U'dot = -(G/r^2) W' + (mu/r^3) u', W'dot = -(G/r^2) U' - (mu/(2r^3)) w' (printed as read; the Hamiltonian (49) has -(mu/2r^3)(2u'^2 - w'^2), which gives these coefficients mu/r^3 and mu/(2r^3) with the signs as printed on the page; the Lagrangian form is (51)).

Matrizant R' = K' R L', with K' = rotation by theta of the positions and of the velocities (matrix with cos, sin, -sin, cos in two blocks), L' = K'^T (printed with cos, -sin, sin, cos): A' = K' A, B' = B L' (52), R'(t; t0) = A'(t) B'(t0) (53). Table III (A') and Table IV (B') entries, READ p.332:

- A': a'11 = 2r - 3 rdot t, a'12 = 0, a'21 = -3 (G/r) t, a'22 = r, a'31 = -rdot + 3 (mu/r^2) t, a'32 = -G/r, a'41 = -G/r, a'42 = rdot, a'13 = G y/r, a'14 = -G x/r, a'23 = rY + G x/r, a'24 = -rX + G y/r, a'33 = -G Y/r, a'34 = G X/r, a'43 = rdot Y - mu y/r^2, a'44 = -rdot X + mu x/r^2.
- B': b'11 = -mu/(2 H r^2), b'12 = 0, b'21 = -(1/G)(rdot - 3 (mu/r^2) t), b'22 = -1/r, b'31 = -(mu/(2HG^2 r^2)) e1 + y/(G r), b'32 = x/(G r), b'41 = -(mu/(2HG^2 r^2)) e2 - x/(G r), b'42 = y/(G r) (printed with a typesetting slip: the page labels it b'_{42} with the subscripts garbled; the entry is row 4, column 2), b'13 = -rdot/(2H), b'14 = -G/(2 H r), b'23 = -(1/G)(2r - 3 rdot t), b'24 = 3 t/r, b'33 = -(rdot/(2HG^2)) e1 + r x/G^2, b'34 = -e1/(2 H G r), b'43 = -(rdot/(2HG^2)) e2 + r y/G^2, b'44 = -e2/(2 H G r).

COMPUTED: Tables III and IV agree with K'A and B L' to 1e-15 at two states and two times (they are free of printed slips as read).

### 4.2 Intrinsic (tangential/normal, Hill) frame, Tables V and VI

Angle phi is the angle of the velocity vector with the inertial x axis: X = V cos phi, Y = V sin phi (54). Displacement components (58): u'' tangential (along track), w'' normal (across track), U'', W'' the corresponding velocity components. Hill-form second-order equation (64) for the normal variation: w''ddot + Theta w'' = -2 (phidot/V) delta H with Theta = Vddot/V + 2 phidot^2 - Pi_xx - Pi_yy (65); for two bodies (67): Theta = (mu/(V^2 r^3)) [V^2 - 3 (G^2/(V^2 r^2)) (V^2 - mu/r)], equivalently (mu/(V^2 r^3)) [V^2 + 3 G (phidot - thetadot)]. Two-body derivatives: Vdot = -mu rdot/(V r^2), phidot = mu G/(V^2 r^3), Vddot = (mu/(V r^3)) [2 rdot^2 - (G^2/(V^2 r^2))(V^2 - mu/r)], Pi_xx + Pi_yy = mu/r^3. The tangential displacement follows by the quadrature (66).

R'' = K'' R L'', K'' with phi in place of theta, L'' = K''^T (printed with cos, -sin, +sin, cos), A'' = K'' A, B'' = B L''. Table V (A''), READ p.336: a''11 = 2 r rdot/V - 3Vt, a''12 = G/V, a''21 = -2G/V, a''22 = r rdot/V, a''31 = -V - 3 Vdot t, a''32 = 0, a''41 = -3 V phidot t, a''42 = V, a''13 = 2 G Y/V, a''14 = -2 G X/V, a''23 = (2GX + V^2 y)/V, a''24 = (2GY - V^2 x)/V, a''33 = -V y phidot, a''34 = V x phidot, a''43 = V Y + Vdot y, a''44 = -V X - Vdot x. COMPUTED: agrees with K'' A to 2e-15.

Table VI (B''), READ p.337: b''11 = Vdot/(2H), b''12 = V phidot/(2H), b''13 = -V/(2H), b''14 = 0; b''21 = -V/G - 3 (Vdot/G) t, b''22 = -3 (V phidot/G) t, b''23 = -2 r rdot/(G V) + 3 (V/G) t, b''24 = 2/V; b''31 = (Vdot/(2HG^2)) e1 + Y/(G V), b''32 = (V phidot/(2HG^2)) e1 + G X/V (AS PRINTED), b''33 = -(V/(2HG^2)) e1 + x r rdot/(V G^2), b''34 = -x/(G V); b''41 = (Vdot/(2HG^2)) e2 - X/(G V), b''42 = (V phidot/(2HG^2)) e2 + G Y/V (AS PRINTED), b''43 = -(V/(2HG^2)) e2 + y r rdot/(V G^2), b''44 = -y/(G V). PRINTED SLIPS in b''32 and b''42 (flagged): the last terms are printed "GX/V" and "GY/V" (G in the numerator); COMPUTED against B L'': the correct terms are X/(G V) and Y/(G V) (the residual after the leading term, divided by X/(GV) and Y/(GV), is 1.0000000000 to 3e-14). With the printed terms Table VI differs from B'' = B L'' by 0.17 and 0.83 in these two entries only; all other entries agree to 1e-15. The authors' intrinsic frame has b''11 with a dotted V (Vdot/(2H)), which I read from the image; the checks passed with that reading.

Intrinsic-frame results (READ pp.336-338): (68) grad H . a''_1 = -2H and grad H . a''_j = 0 (j = 2, 3, 4): columns 2 to 4 of A'' span the isoenergetic variations. The velocity-direction variation (70): u''_o = V, w''_o = 0, U''_o = Vdot, W''_o = -V phidot, equal to alpha_o a''_2 + beta_o a''_3 + gamma_o a''_4 with alpha_o = -4 phidot - 2H/G (as printed) and beta_o, gamma_o as printed on p.337 (not re-derived; the printed beta_o, gamma_o contain Vdot X/(V G phidot), whose reading I did not check). Basic isoenergetic displacements u_1 (normal position offset) and u_2 (normal velocity offset) with initial conditions (71) u''_1(t0) = 0, w''_1(t0) = 1, w''_1dot = 0 and u''_2 = 0, w''_2 = 0, w''_2dot = 1; (72) u''_1dot(t0) = 2 phidot(t0), u''_2dot(t0) = 0; explicit components (73) and (74) printed on p.338 (long formulas; not transcribed here, not checked).

## 5. Mapping to Perko's a21, a24, a41, a44

Perko 1981b, Theorem 4 (digest section D, p.195 of Perko): "a_ij the components of the transition matrix Phi(t1, 0) R in (i, j) coordinates (i along V1, j perpendicular ...)", with d_1 = a44 a21 - a24 a41, and eq. (50) (a21 dr + a24 dv)(a41 dr + a44 dv) = mu/V1.

Findings (COMPUTED structure, INFERRED about Perko's definitions: I did not re-read Perko's page in this task, so the frame conventions below are those the digest of Perko reports):

1. The a_ij of Perko are NOT the a_ij of Deprit's Table I. Deprit's a_ij is the factor A(t) (columns are solutions); Perko's a_ij are entries of the matrizant itself, Phi(t1, 0) R = R(t1; 0) composed with a rotation. Naming clash only; the paper cited as the source of "well-known formulas for the components of the transition matrix" is this one in the sense that A(t1) B(0) is the closed form.
2. Perko's index order is the same as Deprit's: rows and columns 1, 2 are position components (i, j) and 3, 4 are velocity components (i, j), the same as (u, w, U, W) here. So a21 = R_21 is the response of the j-position at t1 to a unit change of the i-position at t = 0; a24 the response of the j-position to a unit change of the j-velocity at t = 0; a41 and a44 the same for the j-velocity at t1.
3. The pair in eq. (50), with columns 1 and 4: a21 dr + a24 dv = the j-position at t1 (the miss distance from the perturbing body's circle, zero for a collision) for an initial perturbation dr of the i-position and dv of the j-velocity at the perpendicular crossing; a41 dr + a44 dv = the j-velocity at t1 (zero for zero deflection). This matches the Perko-digest reading that the first factor is the second-species characteristic (Delta = 0) and the second the first-species one (V_inf parallel to V1).
4. What Deprit supplies for it: the 2 x 2 block of rows (2, 4) and columns (1, 4) of R(t1; 0) = A(t1) B(0) is
   each Perko a_ij is, in Deprit's table notation, the sum over k of (Table I row i, column k, evaluated at the state and time at t1) times (Table II row k, column j, evaluated at the state and time at 0); for example the (row 2, column 1) entry is the sum over k of a_2k(t1) b_k1(0); the extra rotation R in Perko's Phi R is a constant rotation of the input frame (alpha = plus or minus pi/2 at tangency) which, in the case computed in section 7, reduces to relabelling the axes. Deprit's own rotations K, L do the analogous job for the orbital and intrinsic frames and could be used to build Perko's (i along V1) frame, but Perko's frame at t1 is aligned with the relative velocity V1 of the particle and the perturbing body, not with the particle's inertial velocity: the two coincide only when the perturbing body's velocity is zero. In the sidereal frame the perturbing body moves at the unit circular speed, so V1 is the relative velocity and the rotation to the (i, j) frame at t1 is the rotation taking the x axis to V1, not the intrinsic-frame angle phi (INFERRED from the digest, to be confirmed against Perko p.195-198 before coding).
5. Sign conventions that matter: reversing the orientation of j flips the sign of rows 2 and 4 together, so the product (a21 dr + a24 dv)(a41 dr + a44 dv) is unchanged. Reversing the orientation of dr flips the signs of column 1 (a21 and a41), which changes the cross terms and the sign of d_1. The orientation of dr relative to the radius vector is therefore the one convention to check in Perko's text before using the constant.
6. Units: Deprit's mu is the gravitational parameter of the primary. In Perko's normalisation the primary has mass 1 - mu_Perko, which tends to 1 as the mass ratio tends to 0; use mu_Deprit = 1 in the matrizant and keep Perko's mu (the small mass ratio) separate in the right-hand side mu/V1 of (50). Symbol clash to avoid in code.

## 6. Numerical checks run (COMPUTED, scratch only)

Scratch directory (outside the repository): the Table I and II transcription (with b44 corrected), the finite-difference Jacobian, comparison with the project's analytic Shepperd STM, and the orbital/intrinsic table checks. Units: mu = 1, planar states in the z = 0 plane.

A. Table I/II closed form against a central finite-difference Jacobian (step 1e-6) of `cyclerfinder.core.kepler.propagate` (4 x 4 block of x, y, xdot, ydot): 42 cases (7 states: ellipse with e about 0.3, ellipse e about 0.7, hyperbola, circular prograde, circular retrograde, near-circular, retrograde ellipse; times dt = 0.7, 2.5, 9.0 (the last is several revolutions); two time origins t0 = 0 and 3). Worst residual: 4.2e-10 relative and 4.2e-9 absolute (the absolute worst is dt = 9.0 where the Jacobian entries are order 30; the residual is the size of the finite-difference noise, not of a formula error). The same table with the printed b44 as read fails by 1.2 to 24 absolute. A(t0) B(t0) = I to 2.9e-14 worst; det R - 1 within 1e-14.

B. The same closed form against the project's analytic STM `cyclerfinder.core.kepler_stm.shepperd_stm` (4 states, 3 dt each, 4 x 4 sub-block of the 6 x 6): worst |R - Shepperd| = 4.3e-14. So the closed form agrees to machine precision with an independent analytic STM already in the tree; the finite-difference gap is purely finite-difference noise.

C. Circular orbits (both senses) and a near-circular orbit are covered by A and B: the closed form holds at e = 0 even though the derivation divided by P^2 + Q^2. Retrograde (G < 0) holds too.

D. Orbital-frame tables: Table III against K' A and Table IV against B L': 4e-16 to 2e-15 (four cases). Intrinsic: Table V against K'' A: 2e-15; Table VI against B L'': exact except the two printed-slip entries described in section 4.2 (0.17 and 0.83).

## 7. Perko's constant at C = -1, V1 = 2 (COMPUTED example, with the caveats of section 5)

Primary at the origin with gravitational parameter 1, the perturbing body on the unit circle moving in the positive sense, the retrograde unit circular orbit (type 3, C = -1). Perpendicular crossing at t = 0 at (-1, 0) with velocity (0, +1) in the sidereal frame (retrograde: angle decreasing from pi); the body starts at (1, 0). Relative angular rate 2, so the first tangential meeting is at t1 = pi/2, at (0, 1) with particle velocity (1, 0) and body velocity (-1, 0): V1 = 2 (the relative speed, in agreement with the Perko digest's V1 = 2 for this case). The state at t1 from `propagate` is (-1e-16, 1, 1, 1e-16). The matrizant R(t1; 0) from Tables I and II (rows and columns in the order x, y, xdot, ydot; agrees with the finite-difference Jacobian to 5e-10):

    [  2.712389   -1       2       -0.712389 ]
    [ -2            1      -1        2        ]
    [  1           -1       1       -1        ]
    [ -3.712389     1      -2        2.712389 ]

(The irrational entries match 3 pi/2 - 4 = 0.712389, 3 pi/2 - 2 = 2.712389 and 3 pi/2 - 1 = 3.712389 to the six digits shown: recognised numerically, not derived.)

With i = the x direction (the direction of V1, up to the orientation question), j = the y direction, the columns used by Perko (column 1 the position along the x axis, column 4 the velocity along y), and rows 2 and 4 (y-position and y-velocity at t1):

    a21 = -2,  a24 = 2,  a41 = -3.712389,  a44 = 2.712389,
    d_1 = a44 a21 - a24 a41 = 2.0 (equal to V1 for this orbit; whether d_1 = V1 holds in general is not established).

So, with those conventions, eq. (50) reads (-2 dr + 2 dv)(-3.712389 dr + 2.712389 dv) = mu/V1 = mu/2: a hyperbola in (dr, dv) with asymptotes dv = dr (the collision characteristic) and dv = 1.3688 dr (zero deflection; 3.712389/2.712389). The orientation of dr relative to the position (the outward direction at the crossing is -x here, so dr as written is inward) is the convention to settle against Perko (section 5, point 5). These numbers are COMPUTED and are consequences of the matrizant plus my reading of Perko's conventions; they are not printed in either paper. The one input independent of my conventions is the matrizant R(t1; 0), checked to machine precision against `shepperd_stm`.

## 8. Printed numbers usable as sourced tests

Essentially none numeric: the paper has no worked example, no table of values, no figure. What it gives are structural identities, all of which are testable and independent of the project's code:

1. R(t; t0) = A(t) B(t0) with Tables I and II (corrected b44) equals the Jacobian of Kepler propagation (done above, tolerance 1e-9 relative with a finite-difference reference or 1e-13 against `shepperd_stm`).
2. A(t) B(t) = I (Tables I and II as inverses): to 3e-14.
3. det A*(t) = -2 H G^2 (p.326); det R = 1 (Liouville, symplectic): 1e-14.
4. R(t0; t0) = I.
5. Identity (7) P^2 + Q^2 = mu^2 + 2 H G^2 and the identities (8): textbook, testable for any Kepler state.
6. (68): grad H dotted into the columns of A'': equal to -2H for column 1 and 0 for columns 2 to 4 (a check of the energy gradient and of the table).
7. The Hill coefficient Theta (67) for two bodies; the second-order normal-variation equation (64).
These are recomputations of the same mathematics, not an independent numerical source; they sit alongside, not replace, the finite-difference parity test. No circular golden is created: the expected values are the paper's closed forms evaluated at the test state.

## 9. Techniques applicable to the project's problems

**#899 (Perko eq. 50 constant, avoided-crossing fit at C = -1, V1 = 2).** Section 7 gives the Kepler transition matrix R(t1; 0) for the retrograde unit circle and the four entries a21, a24, a41, a44 under stated conventions, and the leading-order hyperbola constant mu/V1 = mu/2 with asymptotes dv = dr and dv = 1.3688 dr. To test it against an avoided-crossing fit: continue the type-3 family through the junction at (x0, C) = (-1, -1) at two or three values of mu below the Earth-Moon value; map the fitted (dx0, dC) coordinates to Perko's (dr, dv) through the Jacobi integral at the crossing (the linear map the Perko digest says is not printed: dC/dydot0 at fixed x0, from C = 2 Omega - v^2, INFERRED: at fixed x0, dC = -2 ydot0 dydot0 at a perpendicular crossing, plus the position derivative of 2 Omega); fit the product constant and compare with mu/2 (slope in mu, exponent 1/2 on the gap). This depends on the dr orientation and the choice of frame in section 5, so first re-read Perko pp.195-198 on the page images and settle those conventions, then code the test. The closed form here gives the matrizant for any other bifurcation orbit (type 2 as well, a different t1 and V1) with no new derivation. `core/kepler_stm.py` (`shepperd_stm`) already supplies the same quantity numerically; the closed form is a second independent route (section 6, B), useful as a cross-check of which entries to use.

**#928 (an analytic Kepler STM as the P = 0 control for the KS propagator's transition matrix).** The tree already holds `core/kepler_stm.py` (Shepperd). This paper adds an independent planar closed form with no universal-anomaly iteration beyond the state propagation: the P = 0 limit of a Levi-Civita/KS transition matrix, converted to Cartesian coordinates and time, must equal A(t) B(t0) from Tables I and II (positions and velocities at the final time against the initial state; the time-regularised STM must be converted using dt = r ds first). Control chain: Table I/II closed form = Shepperd STM = finite-difference Jacobian of the KS integrator at P = 0. The closed form has a useful property for the KS check: it is analytic in the state and the elapsed time, no Kepler solver, so it tests a KS integrator at a given fictitious time only after the physical time is read off. Constraint: the planar closed form excludes rectilinear and parabolic orbits (and, in the derivation, circular ones, though the tables hold there as tested); the radial-fall control for the Levi-Civita code (the `#928` control list) therefore cannot be checked here, since G = 0.

**#906 (demanded-turn gate hardening).** Indirect relevance only: the gate's near-resonant, near-bifurcation cases are the orbits for which Perko's a21 and a24 matter (a21 a24 not zero was stated by Perko, citing this paper). The closed form lets the project evaluate, for any near-tangent first-to-second-species bifurcation seed, the two linear forms a21 dr + a24 dv and a41 dr + a44 dv, and so the distance of a candidate from the bifurcation point in the product form of (50): a candidate whose product is near mu/V1 is on the hyperbola neck and the first-order demanded turn should be reported as "indeterminate" in line with the #906 amendment. It supplies no new tolerance by itself.

Not applicable: #895 (the Titania-Oberon continuation) uses full-model transition matrices; the planar two-body closed form is only a control for a Kepler arc. The paper is planar only (the title says so); the three-dimensional matrizant is announced as future work (the closing paragraph, p.339, hopes for an airborne-computer form "when they are extended to the three dimensional problem"), and the project's 6 x 6 `shepperd_stm` already covers three dimensions.

## 10. Source honesty and limits

- The paper contains no numerical example or table of values; every number in this note except the formulas is COMPUTED by me.
- Three printed slips found by reading against numerical Jacobians (section 3.2 b44; section 4.2 b''32 and b''42; the likely 3Q/2H in (35b), which is not used in the tables): none affects the structure, all are in transcription of a typeset formula. All were found only because the tables were checked numerically; do not code any entry straight from the printed page without the parity test.
- Not checked: the explicit components (73), (74) of the basic isoenergetic displacements (p.338), the printed constants alpha_o, beta_o, gamma_o of (70) (the printed beta_o and gamma_o have a term with phidot in the denominator that I did not re-derive), the lengthy derivation steps (22), (23), (36), (37) (the final tables are verified numerically, which makes those steps indirectly checked).
- The mapping to Perko's a_ij rests on the Perko digest, not on a fresh reading of Perko's frame definitions. Reading Perko pp.195-198 on the images remains a precondition for any code that uses the constant of (50).
- Persee's cover page says the PDF was generated 27/05/2025; the Persee citation is "tome 3, fascicule 3, 1968".

## 11. Follow-ups

1. Re-read Perko 1981b pp.195-198 on the page images for the exact definitions of dr, dv, the rotation R, the orientation of the i and j axes and the sign of V1; then fix the dr orientation (section 5, point 5) and the value of d_1 in section 7. Write the section 7 numbers into the #899 plan as the prediction to test.
2. Derive the linear map (dr, dv) to (dx0, dC) at the type-3 junction (C = -1, x0 = -1) from the Jacobi integral and test the eq. (50) constant mu/V1 = mu/2 against a fitted avoided crossing of the Casoliva-type continuation at two or three mu (also the exponent 1/2 of the gap in mu).
3. Optional: add the closed form A(t) B(t0) (with the b44 correction) as an independent cross-check of `core/kepler_stm.py::shepperd_stm` in the test suite, with the finite-difference parity at the cases of section 6; it checks the planar sub-block only. Decide whether it deserves a place in `src/` or stays a test fixture.
4. For #928: in the KS/Levi-Civita transition-matrix control, add the Table I/II closed form as the P = 0 reference at the physical time reached.
5. Candidate acquisitions if a universal-variable or three-dimensional matrizant derivation is wanted (check CORPUS_INDEX first; none is claimed held): Goodyear 1965 (Astron. J. 70:189) and 1966 (NASA CR-522) for the universal-time matrizant, Danby 1965 (Astron. J. 70:155; AIAA J. 3:769), Sconzo 1963 and 1967, all listed on p.339 of this paper.
6. The printed slips found (b44, the two intrinsic-frame entries, (35b)) are typesetting slips of the 1968 printing, recorded here with the numerical evidence; no other action.

Note 2026-10-05 (Danby 1965 now held): `docs/notes/2026-10-05-digest-danby-1965-matrizant-keplerian-motion.md` digests Danby, AIAA J. 3:769-770, the 3D Kepler matrizant (two closed forms; elliptic and hyperbolic, a finite). It closes follow-up 5 for Danby 1965 only (Goodyear and Danby 1964 remain not held) and supplies the 3D case that this paper left as future work: both of Danby's forms agree with the finite-difference Jacobian of `core/kepler.py` to 4e-10 relative and with `shepperd_stm` to 5e-14, and the formulas as printed need no correction (this paper's tables needed three). Danby's form 2 (radial/transverse/normal components, any in-plane axes) is the natural cross-check of the section 7 entries for Perko's a_ij; it does not settle Perko's frame conventions (follow-up 1 remains open).
