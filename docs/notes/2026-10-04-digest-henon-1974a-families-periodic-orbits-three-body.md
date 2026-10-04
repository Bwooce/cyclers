# Digest: Henon 1974a, "Families of periodic orbits in the three-body problem"

Date: 2026-10-05 (Sydney; file name carries the 2026-10-04 series date). Reading, reasoning and independent re-computation; no project code was changed.

Source: M. Henon, "Families of periodic orbits in the three-body problem", Celestial Mechanics 10:375-388 (1974), DOI
10.1007/BF01586865, received 4 March 1974 (Observatoire de Nice). Filed in the private paper corpus as
`henon-1974a-families-periodic-orbits-three-body-problem-celest-mech-10-375-doi-10.1007-BF01586865.pdf`. A 14-page scan
(PDF page n is journal page 374 + n). Every page was read as an image at 130 dpi; the two tables (journal p.383, printed rotated) were
rotated and re-read at 220 dpi. No digit was illegible.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my own integration or arithmetic (scratch scripts, not
committed; recipe in section 5). INFERRED is my reasoning.

Companions: `docs/notes/2026-10-04-digest-hadjidemetriou-1975-restricted-to-general-continuation.md`,
`docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md`,
`docs/notes/2026-10-04-digest-hadjidemetriou-christides-1975-families-planar-three-body.md`.

## 0. What the paper is

READ (abstract p.375): "We show by a general argument that periodic solutions of the planar problem of three bodies (with
given masses) form one-parameter families. This result is confirmed by numerical investigations: two orbits found earlier by
Standish and Szebehely are shown to belong to continuous one-parameter families of periodic orbits. In general these orbits have a
non-zero angular momentum, and the configuration after one period is rotated with respect to the initial configuration. Similar
general arguments show that in the three-dimensional problem, periodic orbits form also one-parameter families; in the
one-dimensional problem, periodic orbits are isolated."

It is the paper that removed the conjecture (Szebehely and Peters 1967, Standish 1970, Szebehely 1970, Szebehely and Feagin
1973) that general three-body periodic orbits are isolated for fixed masses. The key fact (p.376): "periodic" means the mutual
distances are periodic, so the configuration may be rotated by an angle phi after one period; searches that insisted on
phi = 0 (fixed axes) found only isolated orbits. Both numerical families have masses 3/12, 4/12, 5/12 (the Pythagorean
masses) and are continued from the starting-at-rest orbits of Standish and Szebehely. Those two orbits
have all velocities zero at t = 0; the continued members have none of the symmetries of the starting orbits (p.385, comment c).

## 1. The dimensional argument (READ pp.375-377, section 2)

Fixed masses, planar problem. Phase space is 12-dimensional (two coordinates and two velocities for each of 3 bodies). There
are 6 isolating integrals: the four centre-of-mass coordinates (x0, y0, u0, v0), the energy E and the angular momentum A. A
solution with given integrals lives on a 12 - 6 = 6-dimensional submanifold S. For periodicity the final configuration must equal
the initial one up to a displacement, which is a translation of the centre of mass (components u0 T, v0 T, fixed by the
integrals) and a rotation by an angle phi (mod 2 pi). So there are 8 variables: 6 initial coordinates on S, the period T and
phi; and 6 conditions (the final configuration is in the 6-dimensional S too). Hence 8 - 6 = 2 free parameters. But a periodic
orbit trivially generates a two-parameter set (shift of the time origin; rotation of the initial configuration about the centre
of mass). Excluding these: "for given values of the integrals, there exist only isolated periodic solutions" (p.376).

Now let the 6 integrals vary: 6 parameters, of which 4 (x0, y0, u0, v0) are removed by changing the reference frame and one
(lambda) by the homologous scaling (distances times lambda^2, times lambda^3, velocities lambda^-1; this is the paper's
scaling, so distances lambda^2, time lambda^3). One parameter remains and it appears non-trivial. Main result (p.377, in italics):
"Periodic solutions of the planar three-body problem form one-parameter families (again, for given masses of the three bodies)."
He states plainly that the reasoning is heuristic and not a rigorous proof, strongly supported by the numerics.

Extensions (section 5, p.384): the same counting with (1, 2, 3) dimensions gives phase space (6, 12, 18), integrals
(3, 6, 10), submanifold dimension s = (3, 6, 8), extra variables (1, 2, 2) (T, and phi in 2 and 3 dimensions), free parameters (1, 2, 2) of
which one is the time shift and, for 2 and 3 dimensions, one is the rotation about the angular-momentum axis; the family
parameters are (0, 1, 1) after removing the frame (2, 4, 6), scaling (1) and, in 3 dimensions, the orientation of the angular
momentum (2). Results: periodic orbits of the rectilinear problem are isolated (Schubart 1956 found such); of the planar and
three-dimensional problems form one-parameter families. Free masses (comment b, p.385): (2, 3, 3) parameters in 1, 2, 3
dimensions. The Szebehely 1970 and Szebehely-Feagin 1973 sequences are sections of the three-parameter planar family under
phi = 0 and m3/m2 = const; the Hadjidemetriou and Christides sequence is a section defined by m2/m1 = const and x30 = const (p.385).

Comment (d), p.385: nothing in the argument assumes anything about the masses, so there is no distinction between the general and the
restricted problems. The elliptic restricted problem is the m3 = 0 case of the general problem with the primaries' eccentricity
as the family parameter (Broucke 1969). A diagram (A, phi), one point per orbit and a curve (a "characteristic") per family; in the
circular restricted problem the characteristic is the vertical line A = m1 m2 (eq. 7, p.385); in the elliptic restricted problem
the horizontal line phi = 0 (mod 2 pi) (p.385-386) and the period is constant along the family; A is related to the primaries'
eccentricity by A = m1 m2 (1 - e^2)^(1/2) (eq. 8, p.386). A further integral I = (2A - 2E - 3 m1 m2)/m3 (eq. 9) tends, as m3 -> 0
with the primaries tending to circular motion, to the Jacobi constant C (eq. 10).

Stability (comment e, p.386): computed "by a method to be described in a forthcoming paper" (the 1975b stability companion is
Hadjidemetriou's, not Henon's); "all orbits were found to be rather unstable: the modulus of the largest eigenvalue ... is of
the order of 10 to 20 for orbits of family 1, and 40 to 60 for orbits of family 2." No individual values are printed.

## 2. Method (READ pp.377-378, section 3, and the Appendix pp.386-388)

Normalisation (conventions (a) to (f), Appendix pp.386-388): G = 1; total mass m1 + m2 + m3 = 1 (eq. 11); origin at the centre
of mass (position and velocity zero); scale fixed by the TOTAL ENERGY (the orbits' distances vary in complex ways, so energy is used):
E = -(1/2)(m1 m2 + m2 m3 + m3 m1) (eqs 2, 13), which reduces to the restricted-problem value E = -(1/2) m1 m2 (eq. 12) and makes the
triangular Lagrange solutions equilateral triangles of side 1 rotating at angular velocity 1; time origin at an extremum of the distance
r23 (r23-dot = 0 at t = 0, eq. 14); the x axis parallel to the direction from body 2 to body 3 at t = 0, so that
y2 = y3 and u2 = u3 at t = 0 (eqs 3, 15). The remaining integral A is the family parameter (section 3 p.377; "most convenient to take
A itself as the parameter").

Initial state (eqs 4, 5): the four coordinates of body 1, (x1, y1, u1, v1), are the free initial data. Then
y2 = y3 = -m1 y1/(m2 + m3) and u2 = u3 = -m1 u1/(m2 + m3) (eq. 4); the remaining x2, x3, v2, v3 follow from the centre-of-mass
conditions sum m_i x_i = 0 and sum m_i v_i = 0, the angular momentum sum m_i (x_i v_i - y_i u_i) = A, and the energy
(1/2) sum m_i (u_i^2 + v_i^2) - sum_{i<j} m_i m_j [ (x_j - x_i)^2 + (y_j - y_i)^2 ]^(-1/2) = E (eq. 5): the first three solve for x3, v2, v3 and
the energy gives an implicit equation for x2 solved iteratively, with an approximate x2 supplied because it generally has more than one root.
(The tables give x1, y1, u1, v1, x2, v2; the other six initial coordinates are recovered from eqs 4 and 5, p.379.)

Integration: Bulirsch-Stoer (1966), Waldvogel (1972) regularisation (necessary because many orbits have close approaches or collisions),
16 significant digits; the integrals A and E preserved to better than 1e-14 at the end (p.378). Integration stops at another extremum of r23
(the appropriate one, depending on the orbit). The angle phi between the x axis and the direction from body 2 to body 3 is measured and the
configuration is rotated by -phi; the resulting (x1, y1, u1, v1) are compared with the initial values and the differences driven to zero by
adjusting the four initial coordinates by a standard differential correction (Szebehely and Peters 1967): the integration is repeated with
an increment of 1e-7 in each initial coordinate in turn, giving a 4 x 4 variational matrix; "after 3 or 4 iterations, a limiting
accuracy of the order of 10^-13 is attained" (this depends critically on the stability of the orbit). Continuation: the family parameter A is stepped
and the previous solution gives the start. So the corrector has four unknowns and four conditions (a square system), the energy
and angular momentum being enforced by construction.

## 3. Results as printed (READ pp.378-383)

Start (p.378): Standish's (1970, Fig. 8) periodic orbit, masses m1 = 3/12, m2 = 4/12, m3 = 5/12 (eq. 6), normalised E = -47/288
(COMPUTED check: -(1/2)(12 + 20 + 15)/144 = -47/288 = -0.163194...). "Family 1" is found "without any particular difficulty"; its A = 0
member is Standish's orbit (all initial velocities zero, so A = 0 and phi = 0). The family continues to negative A by reversing the
signs of A, u1, v1, v2, phi in Table I (equivalent to time reversal, p.379-381). Curiosity noted (p.381): "orbits with a positive angular
momentum have rotated, after one period, by a negative angle phi", so there is no general connection between A and phi.
"Family 2" is started from Szebehely's (1970, Fig. 10) orbit, also at rest and with the added peculiarity of a collision at t = T/2; same masses;
its A = 0 member is Szebehely's orbit, "the collision exists only for that particular orbit" (p.382). Tables are given for positive A only.
Figures: 1 to 3 family 1 at A = 0, 0.002, 0.007 (fixed axes; one period; dots are initial positions; end state is the initial state rotated by phi);
4 family 1 at A = 0.007 in axes rotating at omega = phi/T, where the curves are closed (representation not unique as phi is modulo 2 pi);
5 and 6 family 2 at A = 0.025 in fixed and rotating axes (pp.379-382). Fig. 1 labels masses 3/12, 4/12, 5/12 as full, dashed and dotted lines.

### Table I (READ p.383), family 1, masses 3/12, 4/12, 5/12, E = -47/288

Columns: A, x1, y1, u1, v1, x2, v2, T, phi. (The other six initial coordinates follow from eqs 4, 5; y2 = y3 = -(m1/(m2+m3)) y1 = -y1/3,
u2 = u3 = -u1/3.) Eight decimals printed.

| A | x1 | y1 | u1 | v1 | x2 | v2 | T | phi |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.54402539 | 1.79622952 | 0 | 0 | -1.03258028 | 0 | 20.25306432 | 0 |
| 0.001 | 0.54053935 | 1.79816889 | -0.00585553 | -0.01564698 | -1.03063119 | 0.00460907 | 20.25265542 | -0.26835647 |
| 0.002 | 0.53110395 | 1.80293934 | -0.01103456 | -0.02980389 | -1.02554227 | 0.00866815 | 20.25151930 | -0.51726610 |
| 0.003 | 0.51867478 | 1.80837958 | -0.01512090 | -0.04133431 | -1.01914948 | 0.01174645 | 20.24990518 | -0.72933546 |
| 0.004 | 0.50591701 | 1.81318780 | -0.01822663 | -0.05033904 | -1.01284909 | 0.01387159 | 20.24803759 | -0.90428167 |
| 0.005 | 0.49394848 | 1.81711438 | -0.02064588 | -0.05747761 | -1.00711374 | 0.01528495 | 20.24602141 | -1.05095658 |
| 0.006 | 0.48302567 | 1.82026973 | -0.02260728 | -0.06331350 | -1.00199088 | 0.01619576 | 20.24389395 | -1.17746479 |
| 0.007 | 0.47312277 | 1.82280912 | -0.02425644 | -0.06822325 | -0.99741625 | 0.01674488 | 20.24166819 | -1.28938791 |

### Table II (READ p.383), family 2, same masses and energy

| A | x1 | y1 | u1 | v1 | x2 | v2 | T | phi |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.24368035 | 0.64354890 | 0 | 0 | -1.75039533 | 0 | 10.27881780 | 0 |
| 0.005 | 0.24516190 | 0.64443031 | 0.00489610 | -0.01912525 | -1.74990885 | -0.00123263 | 10.27995212 | 0.14763051 |
| 0.010 | 0.24977335 | 0.64716408 | 0.00969552 | -0.03875714 | -1.74841869 | -0.00241317 | 10.28349038 | 0.30110185 |
| 0.015 | 0.25812442 | 0.65209118 | 0.01428654 | -0.05953310 | -1.74581617 | -0.00347977 | 10.28991666 | 0.46826654 |
| 0.020 | 0.27179428 | 0.66015095 | 0.01850316 | -0.08250792 | -1.74183836 | -0.00433877 | 10.30043936 | 0.66357433 |
| 0.025 | 0.29600664 | 0.67463620 | 0.02193970 | -0.11038997 | -1.73571451 | -0.00476788 | 10.31888994 | 0.92905825 |

(The PDF prints digits in groups of two or three; I have joined them. Signs: the minus signs in u1, v1, x2, v2 are as printed; family 1 has
phi < 0 for A > 0, family 2 has phi > 0 for A > 0.)

Period structure: T falls slowly in family 1 (20.253 to 20.242) and rises in family 2 (10.279 to 10.319); family 2's T is almost exactly half of
family 1's (INFERRED observation: 2 x 10.2788 = 20.5576 against 20.2531; not equal, so I do not claim a relation).

## 4. Independent re-computation (COMPUTED)

An inertial Newtonian integrator (DOP853, rtol 1e-13, atol 1e-14, G = 1, masses (3, 4, 5)/12, no regularisation), initial state built from the
printed (x1, y1, u1, v1, x2, v2) with eq. 4 and the two centre-of-mass relations; E, A from the formulas of eq. 5; phi is the angle of the 2-to-3
direction at t = T; "closure" is the largest difference between the initial state and the state at T rotated by -phi (positions and velocities).

All 13 non-collision rows of Tables I and II reproduce:
- E equals -47/288 to within 5.2e-10 (family 1) and 3.6e-10 (family 2); eq. 2/13 is therefore confirmed as the normalisation.
- The angular momentum computed from the printed initial state equals the A of the row to within 3.0e-9 (family 1) and 4.5e-9 (family 2).
- The angle phi at the printed T equals the printed phi to within 3.1e-7 (family 1) and 7.6e-7 (family 2): the sign convention (phi measured as the
  direction 2 to 3, state rotated by -phi) is confirmed, including phi < 0 for family 1 and phi > 0 for family 2.
- The state at T, rotated by -phi, equals the initial state to within 8.7e-7 (family 1) and 6.7e-7 (family 2): the printed T is a full period of the
  rotated configuration, to about the precision the 8-decimal inputs allow on orbits whose largest eigenvalue modulus is 10 to 60 (p.386).
- r23-dot = 0 at both t = 0 (exactly, by construction of eq. 15) and t = T (to 1.1e-6 or better): T is an extremum of r23 as the paper states.
Family 2, row A = 0 (Szebehely's collision orbit): E = -47/288 to 3e-10 and A = 0 reproduce, but the orbit has a collision between bodies 1 and 3 at t = T/2
(COMPUTED: the 1-3 distance is 1.5e-2 at t = 0.4999 T, 6.8e-4 at t = 0.499999 T and falls toward zero at t = T/2), so an unregularised integrator cannot be run through
it and the closure test is not applicable to that row.

Close approaches (COMPUTED, sampled at 2e5 points over T, so lower bounds on the true minima): family 1 minimum separations between bodies 1-2, 1-3, 2-3 are about
1.07, 0.0033 and 0.0018 at A = 0.001 down to 0.92, 0.0018 and 0.00045 at A = 0.007 (so bodies 2 and 3 pass within 5e-4 to 2e-3 of each other); family 2
(A = 0.005 to 0.025): 1-2 about 0.05 to 0.11, 1-3 about 6e-4 to 2.4e-3, 2-3 about 0.8 to 0.92. Orbits of this kind need regularisation or a high-order
integrator with tight tolerances (the paper used Waldvogel's regularisation); my DOP853 at 1e-13 coped with all but the collision row.

## 5. Recipe for the re-computation

State from the table row: y2 = y3 = -m1 y1/(m2 + m3) = -y1/3, u2 = u3 = -u1/3, x3 = -(m1 x1 + m2 x2)/m3, v3 = -(m1 v1 + m2 v2)/m3, with m = (3, 4, 5)/12.
Integrate to the printed T; at T compute phi = atan2 of (r3 - r2); rotate all positions and velocities by -phi and compare with t = 0; check E = -47/288,
A = sum m_i (x_i v_i - y_i u_i) equal to the row's A. Time per row: seconds.

## 6. Connection to the Hadjidemetriou orbits: none that is printed or computable

Answer: no orbit in Henon's two tables is the same as, or connected through anything printed to, any orbit in Hadjidemetriou 1975b Table I or
Hadjidemetriou and Christides 1975 Table I.
1. Masses. Henon's families have m1:m2:m3 = 3:4:5. 1975b is m1 = m2 = m3 = 1/3. Hadjidemetriou and Christides vary m3 with m1 = m2 = (1 - m3)/2
   (equal primaries). No mass ratio in Henon's paper appears in either of theirs (the nearest is the 5:5:3 end of H&C's path, 0.38775:0.38775:0.2245).
   Under Henon's own observation (comment b, p.385) a one-parameter sequence at fixed masses cannot connect to another set of masses without moving in the mass
   parameters; the paper does not attempt it.
2. Conventions. Henon: inertial frame, total mass 1, energy normalised to E = -(1/2) sum m_i m_j (eq. 13), time origin at an r23 extremum, start data
   (x1, y1, u1, v1) with the x axis along 2 to 3. Hadjidemetriou (both): total mass 1, theta'0 = 1 in the frame tied to the P1-P2 line, so E, p vary along
   the family; time origin at the perpendicular crossing of the third body. These are different normalisations (a homologous rescaling relates them, but only
   if the energy scale is known), and the frames differ (inertial against rotating with the binary). For equal masses with m3 = 0 they coincide up to the
   scale: the H&C row 1 orbit has E = -1/8 = -(1/2) m1 m2, exactly eq. 12.
3. A cross-check that is possible and passes (COMPUTED, using only printed quantities of H&C and Henon's eq. 8): the final H&C orbit (their row 19, m3 = 0) is an elliptic-restricted
   periodic orbit with primaries' eccentricity e = 0.292 (computed 0.29221). Rescaling it to Henon's normalisation (semi-major axis 1, E = -1/8) multiplies lengths by
   1/0.689661 and angular momentum by 1/sqrt(0.689661): its angular momentum 0.198553 (computed, from the H&C digest) becomes 0.239088, and Henon's eq. 8 gives
   m1 m2 (1 - e^2)^(1/2) = 0.25 x 0.95634 = 0.239089. The two agree to 2e-7, confirming that the two papers' conventions describe the same physics in the restricted limit. This
   is a consistency check of conventions, not a connection between orbits.
4. Henon himself notes (p.385, comment b, citing the then-submitted Hadjidemetriou and Christides paper) that the Hadjidemetriou-Christides sequence is a section of a
   three-parameter family defined by m2/m1 = const and x30 = const. So the H&C path (fixed x30 = 0.181) and the 1975b equal-mass family (fixed masses,
   x30 varying) are two different sections through the same three-parameter set; Henon's contribution is the dimensional count and that statement. This is the same fact
   the H&C digest reached from the other direction (their path never reaches equal masses).
5. Therefore the 1975b starting orbit and path remain unprinted; Henon's paper does not supply them either.

## 7. Test-ready numbers

Sourced (every value printed):
1. Tables I and II (p.383): 13 non-collision rows plus the two A = 0 rows. For each row: construct the state with eq. 4 and the centre-of-mass relations, then expect (a) E = -47/288 (tolerance 1e-9), (b) A equals the
   row's A (1e-8), (c) at the printed T the angle phi of the 2-to-3 direction equals the printed phi (1e-6), (d) the state at T rotated by -phi equals the initial state (2e-6), (e) r23-dot = 0 at t = 0 and T (1e-5 at T).
   These tolerances pass with the integrator above; the paper's own corrector closes to 1e-13 with the 13-digit values it holds, which are available only from the author.
2. Masses 3/12, 4/12, 5/12 and E = -47/288 (eq. 6, p.378).
3. Unstable-eigenvalue magnitudes: family 1 of order 10 to 20, family 2 of order 40 to 60 (p.386); the largest-eigenvalue modulus, order of magnitude only, so a test can only bracket (for example, 5 < |lambda_max| < 100).
4. Structural identities from eqs 7 to 10: A = m1 m2 for circular primaries, A = m1 m2 (1 - e^2)^(1/2) for elliptic ones, I = (2A - 2E - 3 m1 m2)/m3 -> Jacobi constant C as m3 -> 0; the check of section 6 item 3 (H&C row 19) is an instance.
5. The negative-A reflection rule: reverse the signs of A, u1, v1, v2, phi (p.379). Test: integrating the reflected state gives the same T and a phi of the opposite sign.
Not usable as sourced numbers: the 13-digit values (available on request only), Fig. 1 to 6 coordinates (graph only), per-row stability values (not printed).

## 8. Reconciliation with the project's code

No general three-body code in `src/` (as in the other two digests; searched earlier: no isoenergetic, routh or hadjidemetriou hits; `nbody/` is the restricted REBOUND wrapper and `core/*` are restricted models). Nothing in Henon's tables can be reproduced by a project module. Henon's eq. 8 and eq. 10 are statements about the restricted limit of the project's own models: A = m1 m2 (1 - e^2)^(1/2) is the elliptic restricted problem's primaries' angular momentum at e (in this normalisation) and I -> C is the Jacobi constant `core/cr3bp.py::jacobi_constant`; both could be added as checks to the elliptic-problem code (`core/er3bp.py`) if its primaries' angular momentum is exposed (not checked). The Appendix's conventions (a) to (f) reduce to the standard restricted-problem conventions the project already uses.

## 9. Techniques applicable to the project's problems

Entries read from `data/OUTSTANDING.md` first: `#931` (a) to (c), `#890`, `#895`.

### 9.1 `#931` (a): the equal-mass test integrator and its continuation test

1. Third sourced control set. Henon's Tables I and II give 13 rows (plus Szebehely's collision orbit) of 3:4:5 mass periodic orbits with a tight closure test (the rotated state returns to 7e-7). They need no rotating-frame bookkeeping: inertial integration and the rotation by -phi suffice, so they test the plain Newtonian integrator and the closure logic independently of the Hadjidemetriou frame. The module should take masses as parameters (it already must for the H&C table), so these rows run in the same test file with m = (3, 4, 5)/12.
2. Non-trivial positive control for the closure test: the rotation by phi. A closure test that compares positions at T in fixed axes (phi = 0 assumed) fails on every row except A = 0; passing here requires rotating by -phi. Henon's point (p.376, p.385 comment a) is exactly this, so the control is the expected failure of the unrotated test and the success of the rotated one (COMPUTED: the unrotated closure error, largest difference between the state at T and at 0, is 0.46 to 1.8 for family 1 rows A = 0.001 to 0.007 and 0.26 to 1.3 for family 2 rows A = 0.005 to 0.025, against 9e-7 after rotating by -phi).
3. Continuation test with a non-trivial rotation: continuing in A at fixed masses (Henon's family parameter) is the other method of continuation from the m3 continuation of H&C. Test: from Table I row A = 0.001, step A to 0.002 with a corrector on the four coordinates of body 1 (a 4 x 4 system as in section 2) and require the printed row at 0.002 within 1e-6; likewise 0.002 to 0.003, and Table II. This checks the corrector's handling of phi and the energy constraint. The family at A = 0 has u = v = 0 (a turning point of the reversed direction), so the corrector must also continue through A = 0 to negative A by the reflection rule; a good test is that A = -0.001 is the reflected row 1.
4. Dimension counting as a test of the continuation code's bookkeeping (section 1): a periodic solution has 12 - 6 = 6 free state variables, plus T and phi, minus 6 conditions, minus 2 trivial parameters, giving zero free parameters at fixed integrals; the corrector's Jacobian after removing the time shift and the rotation must be square and nonsingular (4 x 4 in Henon's formulation). Test: the 4 x 4 variational matrix at a table row is nonsingular (determinant not near zero) while the full 8-variable matrix has a two-dimensional null space (time shift, rotation about the centre of mass).
5. Stability magnitude: the order of the largest eigenvalue modulus (10 to 20 for family 1, 40 to 60 for family 2) is a loose cross-check on any monodromy computation run on these orbits, for the `#931` (b) classifier. They are all strongly unstable, so they exercise the real-hyperbolic regime only, complementing the 1975b Table I (stable through quartet to hyperbolic).

### 9.2 `#890` and `#895` (Titania-Oberon; methods only)

1. Henon's counting argument is the useful frame for what `#890` did, but it applies to a free N-body system with the integrals conserved; his paper treats N = 3 only, and I have not derived the count for four bodies. The point that carries over qualitatively: at fixed integrals a periodic orbit is isolated, so a continuation must move a parameter that is not pinned (energy, angular momentum or a mass), and different choices trace different sections through the same higher-dimensional set. `#890` continued in a moon-mass scale; a continuation in energy at fixed masses is another section and may behave differently (fold, termination). That second section needs a model in which the energy is conserved, which `#890`'s prescribed-circular-moons model is not (the moons are driven). It is a design observation, not a recommendation to build.
2. The periodicity definition matters (p.376): "periodic" means mutual distances periodic, the configuration rotated at the end. In a prescribed-circular-moons model the moons have a fixed angular rate and a periodic orbit of the spacecraft is periodic in the rotating frame of the moons; the corresponding statement for a free system is that closure is up to a rotation. The `#895` real-ephemeris arcs have moon positions imposed by the kernel, so a return to Titania is a closure in the inertial frame at a given epoch and no rotation angle phi enters; I note the distinction only so that a future free-moon model does not copy the fixed-phase closure test.
3. Close approaches: Henon needed regularisation (Waldvogel) for orbits with minimum separations of 5e-4 to 3e-3 of the unit scale, which supports the expectation that `#895` flybys at a few 1e-3 of the moon-planet distance are in the same regime; my own run at those separations succeeded with DOP853 at rtol 1e-13 but not through a collision. Control: Szebehely's A = 0 collision orbit (family 2) is a published case where an unregularised integrator fails at T/2, a positive control for any regularised propagator (`#928`).
4. Stability expectation: even for the mild-looking orbits here the largest multiplier is 10 to 60 per period; `#890`'s reported 8.4e5 per cycle is far larger, but the cycle is much longer in revolutions, so the order-of-magnitude comparison is not informative by itself. A control: report the multiplier per flyby as well, so that comparison against published multipliers is on equal footing.

## 10. Recommended follow-ups (no task numbers registered)

1. Add Henon Tables I and II to the `#931` (a) test module (13 rows; E, A, phi, closure, r23-dot at T), with the unrotated-closure expected failure, and the reflection rule for negative A.
2. A continuation-in-A test through A = 0 and onto negative A (section 9.1 item 3), using Table I as the expected output.
3. A matrix-rank test of the periodicity Jacobian (4 x 4 nonsingular, 8-variable null space of dimension 2) on a table row.
4. Obtain Standish 1970 and Szebehely 1970 (Periodic Orbits, Stability and Resonances, pp.375, 382) for the starting orbits and Szebehely and Feagin 1973, and Broucke 1969 for the elliptic families that Henon's comment (d) relies on.
5. Ask whether the 1975b Table I equal-mass family can be re-derived as a Henon-section: equal masses, A or E as the parameter, starting from any one row; since 1975b rows are printed with 8 digits, Newton on (x1, y1, u1, v1) at fixed E and A (Henon's scheme) would test whether they lie on one family, closing the question of section 6 without needing the unprinted H&C path.

## 11. Summary

- Argument: counting (12 - 6 = 6 phase variables per fixed integrals, plus T and phi, minus 6 closure conditions = 2 parameters, both trivial: time shift and rotation) shows periodic orbits are isolated at fixed integrals; letting the integrals vary and removing the frame (4), the scale (1) leaves one parameter. Periodic solutions of the planar and 3-D problems are therefore one-parameter families (heuristic, not rigorous); the rectilinear ones are isolated. The isolation conjecture came from insisting on phi = 0.
- Numerics: two orbits at rest (Standish; Szebehely, with a collision at T/2), masses 3:4:5, E = -47/288, are continued in the angular momentum A with a 4 x 4 Newton corrector (1e-7 increments, 3 to 4 iterations to 1e-13); Tables I (8 rows) and II (6 rows) print x1, y1, u1, v1, x2, v2, T, phi at 8 decimals. Largest eigenvalue modulus about 10 to 20 (family 1) and 40 to 60 (family 2).
- Reproduced: all 13 non-collision rows (E to 5e-10, A to 5e-9, phi to 8e-7, closure to 9e-7); the collision orbit's E and A.
- Connection to Hadjidemetriou: none. Different masses (3:4:5 against 1:1:1 and m1 = m2) and conventions; the only link is a convention check in the restricted limit (H&C row 19 rescaled satisfies Henon's eq. 8 to 2e-7), and the statement that the H&C and 1975b sequences are different sections of the same three-parameter family.
- Project: no general three-body code; the tables give a third sourced control set for `#931` (a), the A-continuation and phi-closure tests, and a published collision orbit as a regularisation control.
