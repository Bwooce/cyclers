# Digest: Broucke 1969, "Periodic Orbits in the Elliptic Restricted Three-Body Problem" (part A: sections I to III and IV.A to IV.D)

Written 2026-10-05 (Sydney) from a read of the scanned report; equations and tables were read from page images and the
tables were also checked by integration (section 8). Statements are marked READ (printed), COMPUTED (arithmetic or an
integration shown here) or INFERRED (reasoning, not printed). Part B of this digest (families 8P, 8A, 11P, 11A, 10P, the
collision-orbit family and the programs, printed pp51 to 124) is in
`docs/notes/2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-b.md`.

## 0. Citation, file, page offset, scope

- R. A. Broucke, "Periodic Orbits in the Elliptic Restricted Three-Body Problem", Jet Propulsion Laboratory Technical Report
  32-1360, California Institute of Technology, Pasadena, 15 July 1969. Prepared under Contract NAS 7-100 for NASA. NASA NTRS
  19700005781. No DOI (a JPL technical report). The preface states the work was done in the Mission Analysis Division, starting
  1 July 1967 (the end year is not legible in the text layer and was not read from the image).
- Filed in the private paper corpus as
  `broucke-1969-periodic-orbits-elliptic-restricted-three-body-problem-jpl-tr-32-1360-ntrs-19700005781.pdf` (135 PDF pages,
  a scan with a poor text layer; the abstract and contents are roman-numbered).
- **Printed page p is PDF page p + 10** (PDF page 11 is printed page 1; checked on pp1, 3, 33, 42, 54). Page numbers below
  are printed pages.
- Problem: the PLANAR elliptic restricted three-body problem (Section II-A: "two-dimensional"); the three-dimensional form
  appears only inside the series programs (Table 1) and "no research has been done on the three-dimensional case" (p19).
- The report: 15 families, 1127 orbits, a theory of stability with seven types, and a collision-orbit family. Part A covers
  everything up to and including Family 7A.

## 1. Headline

1. The planar equations of motion (Eq. 33, p6) are, term for term, the equations `core/er3bp.py` integrates for x and y
   (section 2). The only differences are the third dimension (core has it; the report does not research it) and the fact that
   the report's pulsating form is singular at e = 1 and r = 0, so the e = 1 and collision work uses inertial and Birkhoff
   coordinates instead (core cannot do those).
2. **The report classifies the one-period monodromy of the planar elliptic problem by two invariants (a1, a2) and names seven
   stability regions** (Table 3, Fig. 2, p29; the author calls it "probably the first discussion of the characteristic exponents
   of a nonconservative dynamical system with two degrees of freedom", p2; that priority claim is the author's, not checked).
   a1 and a2 are exactly alpha and beta of Hadjidemetriou 1975 and b = -k (section 5). It therefore supplies the regime taxonomy
   that `#931` asks for.
3. **Symmetric periodic orbits of the elliptic problem have period 2k pi and need two perpendicular crossings of the
   syzygy axis AT AN APSE of the primaries** (p23); there are two ways to continue every circular symmetric orbit,
   starting at periapsis (family letter P) or apoapsis (letter A); the two families need not be the same (p24). This is
   exactly the published method the project's corrector lacks (`#912`), and the apoapsis start is verified here against
   `core.er3bp` on Table 13 (apoapsis, 120 rows).
4. Tables 4 to 13 (printed pp33 to 55) hold 118 + 152 + 13 + 13 + 22 + 12 + 26 + 23 + 131 + 120 = 630 rows (617 distinct rows
   plus the 13 truncated rows of Table 6) of initial and half-period conditions. **Every distinct row was checked against an independent
   integration** (section 8): 7P and 7A
   (251 rows) reproduce in `core.er3bp` from the printed state; 12A (151 elliptic rows) closes to a perpendicular
   crossing at the half period; Family 9P reproduces down to e = 0.988; the rectilinear rows reproduce in a separate
   inertial integration. These are published positive controls for the project's elliptic code, including an apoapsis start
   and nearly collinear primaries (e = 0.999).
5. **Caution, COMPUTED:** the stability type of 7P and 7A sits on the line a2 = -2 a1 - 2 (an eigenvalue pair near +1), where
   the type depends on the fourth decimal of the initial state. With the printed seven-digit states the monodromy
   gives a different region from the one obtained after re-correcting the state, for several rows (section 9). Any
   classifier must be applied to a re-corrected orbit and must carry a tolerance at k = +-2.
6. Printed anomalies found (section 11): Table 6 and Table 7 number orbits 1 and 2 in opposite order; Table 13 stops at
   e = 0.93 while the text says 0.99; Table 11 rows 9 and 10 are out of order; the sign in Eq. 36c,d describes the mirror of the
   frame the tables use.

## 2. The model as printed, compared term by term with `core/er3bp.py`

### 2.1 Underlying two-body problem (Section II-A, pp3 to 4)

READ. Units: a = 1 (semimajor axis of the primaries' relative orbit), n = 1 (mean motion), G included in the masses:
m1 = 1 - mu, m2 = mu <= 1/2. Separation r = 1 - e cos E = p / (1 + e cos v) (Eq. 1), p = 1 - e^2. Barycentric inertial positions
of the primaries (Eq. 2a to d): xi1 = -mu r cos v, eta1 = -mu r sin v; xi2 = (1-mu) r cos v, eta2 = (1-mu) r sin v. Kepler's
equation t + chi = E - e sin E (Eq. 3), with the phase constant chi = 0 (periapsis at t = 0) or chi = pi (apoapsis at t = 0).
Derivative relations (Eq. 4, primes = d/dt): v' = p^(1/2)/r^2, r' = e sin v / p^(1/2), E' = 1/r. Energy of the Kepler motion
(Eq. 5). With v as the independent variable (dots = d/dv): r_dot = e r^2 sin v / p, r_ddot = 2 r_dot^2 / r + r(1 - r/p) (Eq. 6),
and r_ddot = -(2/p) r^3 + (3/p) r^2 - r (Eq. 8, used for the series because it is polynomial; invalid for p = 0 or e = 1).
Transformation of derivatives (Eq. 9): F' = (p^(1/2)/r^2) F_dot, F'' = p (r F_ddot - 2 r_dot F_dot) / r^5.

### 2.2 The equations in each frame

All READ, p4 to p9.

| Frame | Equations | Where |
|---|---|---|
| Inertial barycentric, time | xi'' = -(1-mu)(xi - xi1)/s1^3 - mu (xi - xi2)/s2^3 and the same for eta (Eq. 13); s1^2 = (xi - xi1)^2 + (eta - eta1)^2 (Eq. 14) | p5 |
| Inertial, true anomaly, reduced (xi_bar = xi/r) | L = (1/2)(xi_dot^2 + eta_dot^2) + (1/2)(r/p - 1)(xi^2 + eta^2) + (r/p)(m1/r1 + m2/r2) (Eq. 21); xi_ddot = (r/p - 1) xi - (r/p)(m1 (xi - xi1)/r1^3 + m2 (xi - xi2)/r2^3) (Eq. 23) | p5 |
| Rotating barycentric (rotation angle = v), time | Eq. 27; Lagrangian Eq. 26 | p6 |
| **Rotating-pulsating (x, y) = (x_bar/r, y_bar/r), v independent** | **x_ddot - 2 y_dot = (r/p)(x - m1 (x - x1)/r1^3 - m2 (x - x2)/r2^3), y_ddot + 2 x_dot = (r/p)(y - m1 y/r1^3 - m2 y/r2^3) (Eq. 33)**; x1 = -mu, x2 = 1 - mu, y1 = y2 = 0; Lagrangian Eq. 32; r/p = 1/(1 + e cos v) | p6 |
| Same, time | Eq. 30a,b (r^2 x'' - 2 y' p^(1/2) + 2 r r' x' - x/r = ...) and Eq. 34, t_dot = r^2 / p^(1/2) | p6 to 7 |
| Rectilinear e = 1, inertial, with new variable s (dt = r^2 ds) | Lagrangian Eq. 35; r_ddot = r^2 (Eq. 40), energy Eq. 41, r_ddot = -2 r^3 + 3 r^2 (Eq. 42); coordinates shifted by x0 = (1 - 2 mu)/2 (Eq. 43 to 46) so the primaries sit at +-1/2; equations Eq. 51 | p7 to 8 |
| Alternate form valid for all e including 1 (r = 0 an equilibrium) | dt = r^2 ds (Eq. 52); xi_ddot = (r - p) xi - r(...) (Eq. 54); with rotation by v, v_dot = p^(1/2) (Eq. 56): x_ddot - 2 p^(1/2) y_dot = r(x - m1 (x - x1)/r1^3 - m2 (y - y1)/r2^3) (Eq. 58; as read from the page image the second attraction term of the x equation shows (y - y1)/r2^3 where (x - x2)/r2^3 is meant, apparently a typesetting slip) | p8 to 9 |

### 2.3 Comparison with `core/er3bp.py`

| Item | Broucke 1969 | `core/er3bp.py` | Verdict |
|---|---|---|---|
| Independent variable | true anomaly v; dots are d/dv (Eq. 33) | `f`, true anomaly; primes d/df | same |
| Frame | barycentric, rotating with v, scaled by r (pulsating) (Eq. 10, 28, 64) | pulsating-rotating | same |
| Primaries | m1 = 1 - mu at x1 = -mu (the larger); m2 = mu at x2 = 1 - mu | primary at -mu with mass 1 - mu; secondary at 1 - mu | same |
| mu | m2, with m1 + m2 = 1 (a = n = 1) | `mu` | same |
| Factor | r/p = 1/(1 + e cos v) | `scale = 1/(1 + e cos f)` | same |
| x equation | x_ddot - 2 y_dot = (r/p)(x - m1(x - x1)/r1^3 - m2(x - x2)/r2^3) | `xdd = 2 yd + scale*(x - grav_x)` | identical (x - grav_x is the bracket) |
| y equation | y_ddot + 2 x_dot = (r/p)(y - m1 y/r1^3 - m2 y/r2^3) | `ydd = -2 xd + scale*(y - grav_y)` | identical |
| z (only in the series programs, Table 1 row P6) | z_dot' = -z + (r/p) z A with A = 1 - m1 sigma1 - m2 sigma2 | `zdd = -scale*(e cos f z + grav_z)` | COMPUTED: (r/p) z (1 - m1 sigma1 - m2 sigma2) - z = -(r/p) z (e cos v + m1 sigma1 + m2 sigma2), the same |
| Variational (Eq. 136, 137) | d(dx_dot) = P48 dx + P49 dy + 2 dy_dot with P48 = (r/p)(1 - V3 + 3[m1 (x - x1)^2/r1^5 + m2 (x - x2)^2/r2^5]) | `u_xx = scale*(1 - (1-mu)/r1^3 - mu/r2^3 + 3(...))` | same |
| Time | Eq. 34 t_dot = r^2/p^(1/2); Kepler Eq. 3 | not carried | n/a |
| e = 1 | separate programs (inertial, Birkhoff); the pulsating form fails at e = 1 | `1 + e cos f` vanishes at f = pi | not reachable in core; core reproduces Table 9 down to e = 0.999 |

COMPUTED check that the tabulated states for the pulsating families map to core without a transformation: 7P and 7A rows
integrated in `core.er3bp.propagate_er3bp` from `[x0, 0, 0, 0, ydot0, 0]` at `f0 = 0` (P) or `f0 = pi` (A) reproduce the printed
state at the half revolution (section 8). The mass ratio used in the tables is mu = 0.012155 (6 digits), not the 0.012155099 of
Broucke 1968: with mu = 0.012155 the 7P row 60 half-revolution miss is 3.6e-7 in x1, with 0.012155099 it is 1.9e-5.

## 3. Coordinate changes, equilibria, energy, regularisation (Section II-F to II-K, pp9 to 18)

### 3.1 Coordinate changes (Eqs. 62 to 67, READ pp9 to 10)

- Inertial (xi, eta) to rotating (x_bar, y_bar) (rotation angle v): x_bar = xi cos v + eta sin v, y_bar = -xi sin v + eta cos v;
  x_bar' = y_bar p^(1/2)/r^2 + xi' cos v + eta' sin v, y_bar' = -x_bar p^(1/2)/r^2 - xi' sin v + eta' cos v (Eq. 62; primes d/dt).
  Inverse: Eq. 63.
- Rotating to pulsating: x = x_bar/r, y = y_bar/r, x_dot = (r x_bar' - r' x_bar)/p^(1/2), y_dot = (r y_bar' - r' y_bar)/p^(1/2)
  (Eq. 64); inverse Eq. 65.
- Geocentric xi_G = xi + mu r cos v, eta_G = eta + mu r sin v (Eq. 66); selenocentric xi_S = xi - (1 - mu) r cos v,
  eta_S = eta - (1 - mu) r sin v (Eq. 67). Velocity forms are not given.
- **Derived here (COMPUTED from Eqs. 62 and 64, then verified in `core.er3bp`), for a symmetric start y0 = 0, xi_dot0 = 0 at an apse.**
  Let r_p = 1 - e, r_a = 1 + e, p = 1 - e^2, and eta' = d eta/dt (the printed "ydot" of the inertial tables).
  - Periapsis (v = 0): x = xi0/r_p; y_dot = -xi0/r_p + r_p eta'/p^(1/2).
  - Apoapsis (v = pi): x = -xi0/r_a; y_dot = xi0/r_a - r_a eta'/p^(1/2).
  - Final state from a periapsis start, at v = pi: xi = -r_a x; eta' = -(p^(1/2)/r_a)(y_dot + x).
  The inertial axes are fixed with the x axis along the periapsis direction. At v = pi the primaries therefore sit with m1 at +mu r
  and m2 at -(1-mu) r in the inertial frame.

### 3.2 The five equilibrium points (Section II-G, pp10 to 13)

READ. The five Lagrange points exist in the pulsating frame as fixed points (Eq. 68, the right side of Eq. 33 set to zero);
L4, L5 at x = (1 - 2 mu)/2, y = +-3^(1/2)/2 (Eq. 69); collinear points are the roots of f(x) = -x + (1-mu)(x - x1)/r1^3 + mu (x - x2)/r2^3 = 0,
one in each of the three intervals (Eq. 70 to 72). Proof in inertial coordinates (Eqs. 77 to 84, Fig. 1): L4 is a rigid
Kepler-scaled ellipse (the triangle has sides r(v); angles depend on mu only; the triangle's shape depends on e only);
Eq. 82 shows the net force is along O-L4. The first-order variational equations about the pulsating-frame solution have
periodic coefficients through r/p (Eq. 73 to 75); the formal characteristic equation of Eq. 75 is
lambda^4 + lambda^2 [4 - (r/p)(2 + U_xx + U_yy)] + (r/p)^2 [(1 + U_xx)(1 + U_yy) - U_xy^2] = 0 (Eq. 76), with coefficients periodic in v,
so it gives only a frozen-coefficient picture; Floquet theory is needed.

### 3.3 Energy equation (Section II-H, pp13 to 14)

READ. There is no Jacobi integral. A "median" energy is defined: E = (1/2)(x_dot^2 + y_dot^2) - (r/p) U with
U = -(1/2)(x^2 + y^2) - m1/r1 - m2/r2 (Eq. 86, "U represents the negative of the quantity in brackets"), and its rate is

dE/dt = -(r'/p) U, dE/dv = -(r_dot/p) U = -(e r^2 sin v / p^2) U (Eq. 87).

Eq. 85 gives the general rule dH/dt = partial dH/partial t, dE/dt = -partial dL/partial t. The report integrates Eq. 87 along the
orbit as an accuracy check (the energy-differential equation) and regularises it for the close approaches.
Sign of U (COMPUTED check): with U = (1/2)(x^2 + y^2) + m1/r1 + m2/r2 (the bracket of Lagrangian Eq. 32), E = L2 - L0 gives Eq. 86 and Eq. 87,
and dU/dx = x - m1 (x - x1)/r1^3 - m2 (x - x2)/r2^3 reproduces the right side of Eq. 33, so Eqs. 86, 87, 135, 137 use this U. The sentence under Eq. 86, "U represents the
negative of the quantity in brackets", does not match that (INFERRED: a slip), and the U of Eq. 74 is a different function (m1/r1 + m2/r2 only, without the
(1/2)(x^2 + y^2) term, which appears separately in Eqs. 73 and 75).

### 3.4 Birkhoff regularisation (Section II-I, II-J, pp13 to 17)

READ. Median coordinates (X, Y) = (x - x0, y), x0 = (1 - 2 mu)/2 (Eqs. 88, 89), primaries at X = -1/2 (m1) and +1/2 (m2)
(Eq. 90). Birkhoff variables zeta = xi_B + i eta_B: Z = X + iY = (1/4)(zeta + 1/zeta) (Eq. 92), S = xi_B^2 + eta_B^2 (Eq. 94),
X = (xi_B/4)(1 + 1/S), Y = (eta_B/4)(1 - 1/S) (Eq. 93); Z' = (zeta^2 - 1)/(4 zeta^2) in modulus, Jacobian
J = |dZ/dzeta|^2 = r1 r2 / S (Eq. 97, 99); r1 = ((S+1) + 2 xi_B)/(4 S^(1/2)), r2 = ((S+1) - 2 xi_B)/(4 S^(1/2)) (Eq. 98). New
independent variable s with dv/ds = J = r1 r2 / S (Eq. 100). Lagrangian Eq. 109, L = (1/2)(xi_B_o^2 + eta_B_o^2) + (A_xi xi_B_o + A_eta eta_B_o) + (r/p) J U
(circle = d/ds); equations of motion Eq. 110 with a term in the "median" energy E_m from the change of independent variable
(E_m = E + (r/p) x0^2, Eq. 115); explicit partial derivatives Eqs. 111 to 114; the regularised energy equation
Eq. 116 and the time equations (Eq. 117: v_o = r1 r2 / S, t_o = r1 r2 r^2/(p^(1/2) S)); a seventh-order system. Hamiltonian forms
Eqs. 118 to 123 (H = (1/2)(p_x^2 + p_y^2) + (y p_x - x p_y) - (r/p)[(1/2)(x^2 + y^2) + m1/r1 + m2/r2]). The value of r is obtained
from the two-body formula (Eq. 1), not integrated, in the Birkhoff program (p17).

### 3.5 Ejection orbits (Section II-K, pp17 to 18)

READ. Ejection from m1 (r1 = 0): r2 = 1, xi_B = -1, eta_B = 0, S = D = 1 (Eq. 125); energy Eq. 124 gives
(1/2)(xi_B_o^2 + eta_B_o^2) - m1 r/p = 0 (Eq. 126), so xi_B_o = (2 m1 r/p)^(1/2) cos theta, eta_B_o = (2 m1 r/p)^(1/2) sin theta
(Eq. 127), theta the ejection angle in the Birkhoff plane; ejection from m2: r1 = 1, r2 = 0, xi_B = +1, with m2 in place of m1
(Eqs. 128 to 130). The regularised velocity does not depend on the energy E. An ejection orbit is fixed by only two numbers,
E and theta. Used in Part B for the collision family.

## 4. Series integration and the variational equations (Section III-A, III-B, pp18 to 23)

READ. Recurrent power series (Steffensen) with v as independent variable, after Deprit (Ref. 7). Eighteen parameters
P1 to P18 (Table 1, p19): P1 = x, P2 = y, P3 = z, P4 = x_dot, P5 = y_dot, P6 = z_dot, P7 = r, P8 = r_dot, P9 = t, P10 = x1 (the distance
r1 as the series variable), P11 = r2, P12 = sigma1 = r1^-3, P13 = sigma2 = r2^-3, P14 = A = 1 - m1 sigma1 - m2 sigma2,
P15 = x A + B with B = m1 m2 (sigma2 - sigma1), P16 = y A, P17 = z A, P18 = r^2. For the planar case P3, P6, P17 are zero. Series
P_i = P_i(1) + P_i(2) v + ... + P_i(18) v^17 (Eq. 131, order 17). Table 2 (p20) gives the recurrence relations; three sum
symbols Sigma1, Sigma2, Sigma3 (Eq. 132). The redundant variables give a built-in accuracy check (Eq. 133, 134: |P10 - r1|,
|P11 - r2|, |P12 - r1^-3|, |P13 - r2^-3|, printed average per step). The author notes the method needs a complete
transformation of the equations but that redundancy is a feature.

Variational equations (Section III-B, p20 to 22): Eq. 135 gives the first-order form; Eq. 136 the variational equations
d(dx_dot)/dv = (r/p)(U_xx dx + U_xy dy) + 2 dy_dot, d(dy_dot)/dv = (r/p)(U_yx dx + U_yy dy) - 2 dx_dot with
U_xx = 1 - V3 + 3[m1 (x-x1)^2/r1^5 + m2 (x-x2)^2/r2^5], U_xy = 3y[m1 (x-x1)/r1^5 + m2 (x-x2)/r2^5],
U_yy = 1 - V3 - 3y^2 [m1/r1^5 + m2/r2^5], V3 = m1/r1^3 + m2/r2^3 (Eqs. 137, 138; here U is the bracket (1/2)(x^2 + y^2) + m1/r1 + m2/r2, which is why d2U/dx2 starts with 1).
Two integrators: a predictor-corrector with 20 first-order equations; and the best results from recurrent series with 32 extra
parameters P19 to P50 (Eqs. 139, 141; P35 to P50 are auxiliary products so that no series is ever multiplied with more than one
other; the recurrences for P23 to P34 repeat those of P19 to P22; the summation index is printed "(9)" where "(q)" is meant in several lines of Eq. 141, a
print or scan artefact). Sigma4 is p = 2..n+1, q = n..1 (Eq. 142). Series are evaluated in the order P1-13, P14-18, P19-34,
P35-38, P39-44, P45-47, P48-50.

## 5. Periodicity, families, differential corrections, stability (Section III-C to III-H, pp23 to 29)

### 5.1 Periodicity and symmetry (Section III-C, p23)

READ. The problem is nonautonomous (the independent variable appears even in rotating or pulsating axes). A periodic solution
must therefore have a period that is an integer multiple of the period of cos v, sin v: **2 k pi, k = 1, 2, ...** (also the
period of the rotating axes, so the orbit is periodic in inertial axes too). Only isolated periodic orbits exist at fixed (e, mu),
all "commensurable"; families arise by varying e or mu. Only symmetric orbits (about the syzygy axis, the x axis) were
studied; k = 1 to 5 investigated.

- Circular problem, weak criterion: a perpendicular crossing of the syzygy axis implies symmetry; two crossings give a periodic
  orbit (time of crossing irrelevant).
- **Elliptic problem, strong criterion (Moulton, Ref. 13): the crossing must occur while the primaries are at an apse (periapsis or apoapsis).
  An orbit is periodic if it has two perpendicular crossings of the syzygy axis, and the crossings are at moments when the primaries
  are at an apse.** The time between two apses is a multiple of pi (k pi), so the period is 2 k pi. The criterion is sufficient,
  not necessary (nonsymmetric periodic orbits may exist).

### 5.2 Families and start epochs (Section III-D, p24)

READ. Two parameters (e, mu). The extra parameter e adds one degree of freedom. For each symmetric periodic orbit at e = 0 there are two ways
to continue it: start the integration on the first perpendicular crossing with the primaries at periapsis (v = 0 at t = 0, called a
**periapsis orbit, letter P**) or at apoapsis (v = pi at t = 0, **apoapsis orbit, letter A**). At e = 0 the orbits are identical apart from a
half-revolution phase delay of the primaries; for e > 0 they differ and give different families, which join at e = 0.
In every run the initial state is (x0, y0 = 0, x_dot0 = 0, y_dot0) with free (x0, y_dot0). The families found: one parameter (mu) fixed, e
variable, or the reverse; any parameter can be varied by the programs. 4000 circular orbits (from the earlier circular study, mu = 0.012155) and the orbits of Bartlett
(Ref. 18) and Henon (Ref. 27) (mu = 1/2) were collected on cards; the isolated orbits with period 2 k pi (k = 1 to 5) were extracted by
automatic interpolation (150 resonant orbits). Continuation in e uses extrapolation from the last ten orbits; the corrector varies two of the
initial parameters to meet y = x_dot = 0 at the end of a half orbit, the third is the family parameter (e or mu). For collision orbits the two
corrected parameters are E and e, the family parameter mu.

### 5.3 Differential corrections (Section III-E, pp25 to 26)

READ. General planar case: four variables (x0, y0, x_dot0, y_dot0); x0 + dx0 = x(T, x0 + dx0) (Eq. 143), linearised
[dx(T, x0)/dx0 - I4] dx0 = x0 - x(T, x0) (Eq. 145), with the partials from the variational equations; four to five iterations when the start is good.
In the CIRCULAR problem the matrix is singular (a double +1 eigenvalue); in the elliptic problem the unit eigenvalue does not occur, but the matrix is often near singular for
other reasons (several solutions); a minimum-norm least-squares solution is used (programs by Dr. C. Lawson), "highly satisfactory".
Symmetric case: two variables, x(T, x0 + dx0) = x_F (Eq. 146 to 148); explicitly
(dy/dx0) dx0 + (dy/dy_dot0) dy_dot0 = -y(x0, y_dot0, T), (d x_dot/dx0) dx0 + (d x_dot/dy_dot0) dy_dot0 = -x_dot(x0, y_dot0, T) (Eq. 149), solved by
Cramer's rule (Eq. 150). Here T is the half period k pi. The collision-orbit version corrects (E, e) (Eq. 151).

### 5.4 The fundamental matrix and its eigenvalues (Section III-F, pp26 to 27)

READ. d(dx)/dv = A dx (Eq. 152), A = [[0, I], [a, 2J]] (Eq. 153) with a the symmetric second-derivative matrix and J = [[0, I], [-I, 0]] (Eq. 154, so J here is
the skew block). S = [[0, I], [-I, 2J]] (Eq. 155). A is skew-symplectic with respect to S: S A^T S^-1 = -A (Eq. 156); the principal fundamental matrix R
satisfies dR/dv = A R (157), R^-1 solves dR^-1/dv = -R^-1 A (158), and S R^T S^-1 = R^-1 (159), so R, R^T, R^-1 have the same
eigenvalues and they come in reciprocal pairs **(lambda, 1/lambda, mu, 1/mu)** (Eq. 160), with det = 1. The characteristic exponents are
(alpha, -alpha, beta, -beta) with lambda = e^(alpha T), mu = e^(beta T) (Eq. 161); T = 2 k pi.

### 5.5 The characteristic equation (Section III-G, pp27 to 28)

READ. s^4 + a1 s^3 + a2 s^2 + a1 s + 1 = 0 (Eq. 162): the two outer coefficients are 1 and the two odd ones are equal, and these relations are used to
check the numerical precision ("at least 10^-10"). Stability indices k1 = lambda + 1/lambda, k2 = mu + 1/mu (Eq. 163); factorised
(s - lambda)(s - 1/lambda)(s - mu)(s - 1/mu) = 0 (Eq. 164); **a1 = -(k1 + k2), a2 = 2 + k1 k2** (Eq. 165); k1, k2 are the roots of X^2 + a1 X + (a2 - 2) = 0
(Eq. 166), k = [-a1 +- (a1^2 - 4 a2 + 8)^(1/2)]/2 (Eq. 167); lambda, 1/lambda = [k1 +- (k1^2 - 4)^(1/2)]/2 (Eq. 169) and likewise mu (Eq. 171). From the
four columns (x_i, y_i, x_dot_i, y_dot_i), i = 1..4, of the matrix R: -a1 = x1 + y2 + x_dot3 + y_dot4 (trace) and
a2 = the sum of the six 2 x 2 principal minors (Eq. 172a,b). The CIRCULAR problem has one stability index (two unit eigenvalues, k = lambda + 1/lambda);
the elliptic problem has two, so stability is a PLANE (a1, a2).

### 5.6 The seven types (Section III-H, pp28 to 29, Table 3, Figs. 2 and 3)

READ. Stability if and only if all four eigenvalues lie on the unit circle (all characteristic exponents purely imaginary). Boundaries:
the discriminant of Eq. 166, a2 = a1^2/4 + 2 (parabola, Eq. 173), and the two lines of Eq. 174, a2 = +2 a1 - 2 and a2 = -2 a1 - 2, tangent to the
parabola at (a1, a2) = (+-4, 6); the circular problem (two unit roots, k = 2) lies on the line a2 = -2 a1 - 2.
(COMPUTED from Eq. 165: k = +2 gives a2 = 2 + 2(-a1 - 2) = -2 a1 - 2, and k = -2 gives a2 = 2 a1 - 2.) Seven regions (Table 3, p29; Fig. 2 diagram; Fig. 3 root configurations):

| Region | Name (Table 3) | a1^2 - 4 a2 + 8 | k1^2 - 4 | k2^2 - 4 | Eigenvalues |
|---|---|---|---|---|---|
| 1 | Stability | > 0 | < 0 | < 0 | all four on the unit circle (lambda and mu complex) |
| 2 | Complex instability | < 0 | complex | complex | a complex quartet, off the unit circle |
| 3 | Even-odd instability | > 0 | > 0 | > 0 | lambda, mu real with lambda mu < 0: two positive and two negative real roots |
| 4 | Even-even instability | > 0 | > 0 | > 0 | four real positive roots |
| 5 | Odd-odd instability | > 0 | > 0 | > 0 | four real negative roots |
| 6 | Even-semi-instability | > 0 | > 0 | < 0 | lambda real > 0; mu complex, on the unit circle |
| 7 | Odd-semi-instability | > 0 | < 0 | > 0 | lambda complex on the circle; mu real < 0 |

Fig. 2 (read, label positions): region 2 is labelled at the top of the diagram inside the "a1^2 - 4 a2 + 8 < 0" area (above the parabola); region 1 sits just under the parabola
above the cusp (0, -2); regions 4 and 5 are labelled at the upper left and upper right outside the tangent lines; regions 6 and 7 are labelled at the left and right; region 3 is below the cusp (0, -2).
(INFERRED from Eq. 165 and Table 3: region 4 has k1, k2 > 2, so a1 = -(k1 + k2) < -4; region 5 has k1, k2 < -2, so a1 > 4; region 6 has k1 > 2 and |k2| < 2; region 7 has k1 < -2 and |k2| < 2.)
Only one region is stable and six are unstable. (INFERRED from Table 3's root descriptions) "even" and "odd" refer to the signs of the real eigenvalues: a negative real eigenvalue
(odd) is the period-doubling kind.

### 5.7 Mapping to the project's classifiers (COMPUTED; INFERRED where stated)

| Quantity | Broucke 1969 | Hadjidemetriou 1975b (see its digest) | `search/er3bp_floquet.floquet_classify` | `#931` item (b) |
|---|---|---|---|---|
| Trace term | a1 = -trace R | alpha = -trace A | not exposed | alpha |
| Second invariant | a2 = sum of the six principal 2x2 minors | beta, same sum | not exposed | beta |
| Discriminant | a1^2 - 4 a2 + 8 | Delta = alpha^2 - 4 (beta - 2), identical | not exposed | Delta |
| Index | k = (-a1 +- sqrt(...))/2, k = lambda + 1/lambda | b = (alpha +- sqrt(Delta))/2, **b = -k** | not exposed | b1, b2 |
| Stable | both |k| < 2, Delta > 0 | Delta > 0, |b1| < 2, |b2| < 2 (same) | tag "marginal" (see below) | regime "stable" |
| Complex quartet | Delta < 0 (region 2) | Delta < 0 | tag "unstable" | regime "quartet" |
| Line lambda = +1 | a2 = -2 a1 - 2 | beta = -2 alpha - 2 (same) | n/a | fold / tangent crossing |
| Line lambda = -1 | a2 = 2 a1 - 2 | beta = 2 alpha - 2 (same) | n/a | period doubling |
| Parabola | a2 = a1^2/4 + 2 | beta = alpha^2/4 + 2 (same) | n/a | collision of two pairs |

- The two papers use the same invariants; a classifier built for `#931` can be tested against Table 3 directly: the regime follows from the signs and sizes of
  (k1, k2) as listed in section 5.6, and the region number gives a name for each regime.
- COMPUTED with the project code (the module's docstring reads "stable: all |lambda| <= 1 + tol", the code returns "stable" only when max |lambda| < 1 - tol): `floquet_classify` applied to the monodromy of a re-corrected 7A orbit (row 75, e = 0.5) returns tag
  "marginal" with on_unit_circle True and all |eigenvalue| = 1 (region 1, stable); applied to 7P (row 100, e = 0.21, region 6) returns "unstable" with the
  largest modulus 1.062 (the pair 1.051 and 0.951 is the planar part; the 1.062 pair is out of plane); applied to 7P row 131 (e = 0.5, region 3) returns
  "unstable", largest modulus 2.066. So the tag "stable" is unreachable for a symplectic monodromy (max modulus >= 1 always), and the docstring and code differ. Using it on a region-1 orbit
  gives "marginal". This is a finding for `#931` (b) and is the same class of defect the Hadjidemetriou digest flags.
- The region number and the sign of the real eigenvalues are NOT available from `floquet_classify`; the planar 4 x 4 block is also needed because it takes a 6 x 6 matrix
  and the out-of-plane pair can be unstable on its own (7P row 100: out-of-plane modulus 1.062).

## 6. Section IV: overview (p30)

READ. 1127 orbits in six groups: (1) 13 isolated rectilinear orbits (1P); (2) four short segments with e = 1 or about 1 (3P, 4P, 5P, 9P); (3) two families related to e = 1 (6A, 12A);
(4) two families with the Earth-Moon mass ratio mu = 0.012155 (7P, 7A); (5) five families with equal masses (8P, 8A, 10P, 11P, 11A); (6) one family of periodic collision orbits.
P and A are the periapsis and apoapsis types. Tables 4 to 19 give initial and final conditions; columns with suffix 0 (X0, YDOT0) are initial, suffix 1 (X1,
YDOT1) are final. **Coordinates:** barycentric INERTIAL for families 1P, 3P, 4P, 5P, 9P, 6A, 12A; barycentric ROTATING-PULSATING for the eight other families
(printed p30). Five frames are plotted: barycentric inertial, barycentric rotating, barycentric rotating-pulsating, geocentric inertial (centre m1), selenocentric inertial (centre m2).
"Much effort was made at the beginning to complete each family, but this has turned out to be impossible because, in the investigation of the continuation of
one family, several new families are always discovered."

## 7. Section IV.A to IV.D: families and tables

### 7.1 IV.A Families 12A and 6A (pp30 to 38)

READ. Orbits symmetric about the syzygy axis; strong periodicity criterion; start at apoapsis (maximum elongation), hence "A". The family is a continuation of a period 4 pi orbit of
Stromgren's Class C of symmetric periodic orbits around L1 (equal masses, mu = 0.5; Fig. 4, p30, shows the evolution; orbit 3 is the double-collision orbit; orbit 5 has period 4 pi).
Continuation in e at mu = 0.5 from e = 0.0 up to e about 0.75 with the series programs, then with an inertial Runge-Kutta program (the pulsating equations are singular at e = 1; five to six places required)
to e = 1; the shape in rotating axes is similar for all e (Fig. 5, e = 0.74, p31); in inertial axes (Fig. 6) the orbit shrinks onto a rectilinear oscillation along Oy
for e = 1: the satellite reaches about +-1.69 when the primaries collide and passes between them at maximum elongation. Schubart (1956) had published this rectilinear orbit. **Family 12A:
mu = 0.5, e from 0.0 to 1.0, 152 orbits (Fig. 8, p35; Table 5). Family 6A: e = 1.0, mu from 0.5 to 0.166, 118 orbits (Fig. 7, p32; Table 4); its mu = 0.5 start is the rectilinear
orbit.** For decreasing mu the initial rectilinear orbit bends into a back-and-forth oscillation of period 4 pi along a pseudo-parabolic path, inside an envelope of zero-velocity
points (dashed in Fig. 7). The family 6A is probably extensible to mu = 0 ("a regularization appears to be necessary"). The stability of 6A and 12A "has not been studied" (p32). More
orbits could be computed with fixed e (not 1) and variable mu.

**Table 4, Family 6A, printed pp33 to 34** (118 rows). Header printed: "RECT.ELL.PROBLEM.MU=VARIABLE,E=1.0, (MU=0.5 IS RECTIL.) FAMILY 6"; columns NR, X0, YDOT0, MASS RATIO, ECC.
Barycentric inertial coordinates, state (X0, 0, 0, YDOT0) at t = 0 with the primaries at apoapsis (the e = 1 collision is at t = pi), period 4 pi. All rows have X0 <= 0 (row 1 is exactly 0.0000000).

```
  NR          X0        YDOT0   MASS RATIO         ECC
   1    0.0000000   1.0530847    0.5000000   1.0000000
   2   -0.0050094   1.0530873    0.4990000   1.0000000
   3   -0.0150282   1.0531077    0.4970000   1.0000000
   4   -0.0250467   1.0531486    0.4950000   1.0000000
   5   -0.0350646   1.0532099    0.4930000   1.0000000
   6   -0.0450817   1.0532917    0.4910000   1.0000000
   7   -0.0550978   1.0533940    0.4890000   1.0000000
   8   -0.0651127   1.0535167    0.4870000   1.0000000
   9   -0.0751261   1.0536600    0.4850000   1.0000000
  10   -0.0851378   1.0538237    0.4830000   1.0000000
  11   -0.0951477   1.0540081    0.4810000   1.0000000
  12   -0.1051555   1.0542130    0.4790000   1.0000000
  13   -0.1151609   1.0544385    0.4770000   1.0000000
  14   -0.1251638   1.0546847    0.4750000   1.0000000
  15   -0.1351639   1.0549516    0.4730000   1.0000000
  16   -0.1451611   1.0552392    0.4710000   1.0000000
  17   -0.1551551   1.0555476    0.4690000   1.0000000
  18   -0.1651458   1.0558768    0.4670000   1.0000000
  19   -0.1751328   1.0562268    0.4650000   1.0000000
  20   -0.1851162   1.0565980    0.4630000   1.0000000
  21   -0.1950956   1.0569900    0.4610000   1.0000000
  22   -0.2050708   1.0574032    0.4590000   1.0000000
  23   -0.2150417   1.0578375    0.4570000   1.0000000
  24   -0.2250081   1.0582930    0.4550000   1.0000000
  25   -0.2349699   1.0587699    0.4530000   1.0000000
  26   -0.2449265   1.0592681    0.4510000   1.0000000
  27   -0.2548782   1.0597878    0.4490000   1.0000000
  28   -0.2648248   1.0603291    0.4470000   1.0000000
  29   -0.2747658   1.0608921    0.4450000   1.0000000
  30   -0.2847019   1.0614773    0.4430000   1.0000000
  31   -0.2946312   1.0620836    0.4410000   1.0000000
  32   -0.3045552   1.0627123    0.4390000   1.0000000
  33   -0.3144731   1.0633632    0.4370000   1.0000000
  34   -0.3243849   1.0640364    0.4350000   1.0000000
  35   -0.3342904   1.0647319    0.4330000   1.0000000
  36   -0.3441890   1.0654498    0.4310000   1.0000000
  37   -0.3540818   1.0661909    0.4290000   1.0000000
  38   -0.3639676   1.0669547    0.4270000   1.0000000
  39   -0.3738465   1.0677414    0.4250000   1.0000000
  40   -0.3837184   1.0685514    0.4230000   1.0000000
  41   -0.3935833   1.0693848    0.4210000   1.0000000
  42   -0.4034412   1.0702420    0.4190000   1.0000000
  43   -0.4132914   1.0711226    0.4170000   1.0000000
  44   -0.4231344   1.0720273    0.4150000   1.0000000
  45   -0.4329699   1.0729563    0.4130000   1.0000000
  46   -0.4427978   1.0739097    0.4110000   1.0000000
  47   -0.4526181   1.0748878    0.4090000   1.0000000
  48   -0.4624305   1.0758906    0.4070000   1.0000000
  49   -0.4722352   1.0769190    0.4050000   1.0000000
  50   -0.4820318   1.0779726    0.4030000   1.0000000
  51   -0.4918205   1.0790519    0.4010000   1.0000000
  52   -0.5016011   1.0801573    0.3990000   1.0000000
  53   -0.5113735   1.0812889    0.3970000   1.0000000
  54   -0.5211378   1.0824474    0.3950000   1.0000000
  55   -0.5308937   1.0836324    0.3930000   1.0000000
  56   -0.5406414   1.0848449    0.3910000   1.0000000
  57   -0.5503806   1.0860851    0.3890000   1.0000000
  58   -0.5601115   1.0873532    0.3870000   1.0000000
  59   -0.5698338   1.0886497    0.3850000   1.0000000
  60   -0.5795476   1.0899749    0.3830000   1.0000000
  61   -0.5892529   1.0913295    0.3810000   1.0000000
  62   -0.5989496   1.0927137    0.3790000   1.0000000
  63   -0.6086376   1.0941278    0.3770000   1.0000000
  64   -0.6183170   1.0955725    0.3750000   1.0000000
  65   -0.6279877   1.0970482    0.3730000   1.0000000
  66   -0.6376498   1.0985555    0.3710000   1.0000000
  67   -0.6473030   1.1000946    0.3690000   1.0000000
  68   -0.6569476   1.1016664    0.3670000   1.0000000
  69   -0.6665833   1.1032712    0.3650000   1.0000000
  70   -0.6762103   1.1049098    0.3630000   1.0000000
  71   -0.6858286   1.1065826    0.3610000   1.0000000
  72   -0.6954380   1.1082903    0.3590000   1.0000000
  73   -0.7050386   1.1100337    0.3570000   1.0000000
  74   -0.7146305   1.1118133    0.3550000   1.0000000
  75   -0.7242135   1.1136298    0.3530000   1.0000000
  76   -0.7337878   1.1154840    0.3510000   1.0000000
  77   -0.7385716   1.1164254    0.3500000   1.0000000
  78   -0.7624579   1.1212802    0.3450000   1.0000000
  79   -0.7862893   1.1263905    0.3400000   1.0000000
  80   -0.8100663   1.1317705    0.3350000   1.0000000
  81   -0.8337889   1.1374358    0.3300000   1.0000000
  82   -0.8574577   1.1434038    0.3250000   1.0000000
  83   -0.8810729   1.1496934    0.3200000   1.0000000
  84   -0.9046351   1.1563259    0.3150000   1.0000000
  85   -0.9281448   1.1633246    0.3100000   1.0000000
  86   -0.9516025   1.1707156    0.3050000   1.0000000
  87   -0.9750087   1.1785277    0.3000000   1.0000000
  88   -0.9983642   1.1867933    0.2950000   1.0000000
  89   -1.0216696   1.1955485    0.2900000   1.0000000
  90   -1.0449255   1.2048340    0.2850000   1.0000000
  91   -1.0681326   1.2146951    0.2800000   1.0000000
  92   -1.0912918   1.2251836    0.2750000   1.0000000
  93   -1.1144036   1.2363577    0.2700000   1.0000000
  94   -1.1374689   1.2482838    0.2650000   1.0000000
  95   -1.1604884   1.2610378    0.2600000   1.0000000
  96   -1.1834629   1.2747065    0.2550000   1.0000000
  97   -1.2063931   1.2893905    0.2500000   1.0000000
  98   -1.2292798   1.3052059    0.2450000   1.0000000
  99   -1.2521239   1.3222887    0.2400000   1.0000000
 100   -1.2749260   1.3407982    0.2350000   1.0000000
 101   -1.2976870   1.3609228    0.2300000   1.0000000
 102   -1.3204076   1.3828872    0.2250000   1.0000000
 103   -1.3430885   1.4069615    0.2200000   1.0000000
 104   -1.3657305   1.4334735    0.2150000   1.0000000
 105   -1.3883343   1.4628257    0.2100000   1.0000000
 106   -1.4109007   1.4955178    0.2050000   1.0000000
 107   -1.4334302   1.5321793    0.2000000   1.0000000
 108   -1.4559237   1.5736156    0.1950000   1.0000000
 109   -1.4783821   1.6208768    0.1900000   1.0000000
 110   -1.5008043   1.6753471    0.1850000   1.0000000
 111   -1.5231925   1.7389361    0.1800000   1.0000000
 112   -1.5321383   1.7674878    0.1780000   1.0000000
 113   -1.5410783   1.7981249    0.1760000   1.0000000
 114   -1.5500127   1.8311020    0.1740000   1.0000000
 115   -1.5589428   1.8667278    0.1720000   1.0000000
 116   -1.5678681   1.9053503    0.1700000   1.0000000
 117   -1.5767853   1.9473569    0.1680000   1.0000000
 118   -1.5856986   1.9932978    0.1660000   1.0000000
```


**Table 5, Family 12A, printed pp36 to 38** (152 rows; rows 1 to 60 p36, 61 to 120 p37, 121 to 152 p38). Header: "ELL.PROBLEM.MU=0.5,E=1.0 TO 0.0 (E=1 IS RECT.) FAMILY 12". Barycentric inertial,
apoapsis start, period 4 pi. Row 1 equals Table 4 row 1 (the rectilinear orbit); row 152 is the Stromgren orbit at e = 0.

```
  NR          X0        YDOT0   MASS RATIO         ECC
   1    0.0000000   1.0530847    0.5000000   1.0000000
   2   -0.0055055   1.0536413    0.5000000   0.9990000
   3   -0.0095377   1.0547550    0.5000000   0.9970000
   4   -0.0123140   1.0558693    0.5000000   0.9950000
   5   -0.0145724   1.0569843    0.5000000   0.9930000
   6   -0.0165242   1.0580999    0.5000000   0.9910000
   7   -0.0182709   1.0592159    0.5000000   0.9890000
   8   -0.0198647   1.0603329    0.5000000   0.9870000
   9   -0.0213403   1.0614508    0.5000000   0.9850000
  10   -0.0227208   1.0625691    0.5000000   0.9830000
  11   -0.0240226   1.0636881    0.5000000   0.9810000
  12   -0.0252577   1.0648078    0.5000000   0.9790000
  13   -0.0264356   1.0659282    0.5000000   0.9770000
  14   -0.0275644   1.0670495    0.5000000   0.9750000
  15   -0.0286476   1.0681712    0.5000000   0.9730000
  16   -0.0296924   1.0692937    0.5000000   0.9710000
  17   -0.0307019   1.0704170    0.5000000   0.9690000
  18   -0.0316795   1.0715410    0.5000000   0.9670000
  19   -0.0326282   1.0726657    0.5000000   0.9650000
  20   -0.0335496   1.0737910    0.5000000   0.9630000
  21   -0.0344479   1.0749174    0.5000000   0.9610000
  22   -0.0353229   1.0760444    0.5000000   0.9590000
  23   -0.0361770   1.0771721    0.5000000   0.9570000
  24   -0.0370116   1.0783006    0.5000000   0.9550000
  25   -0.0378280   1.0794298    0.5000000   0.9530000
  26   -0.0386280   1.0805601    0.5000000   0.9510000
  27   -0.0394105   1.0816907    0.5000000   0.9490000
  28   -0.0401786   1.0828223    0.5000000   0.9470000
  29   -0.0409324   1.0839547    0.5000000   0.9450000
  30   -0.0416727   1.0850879    0.5000000   0.9430000
  31   -0.0424003   1.0862219    0.5000000   0.9410000
  32   -0.0431149   1.0873565    0.5000000   0.9390000
  33   -0.0438194   1.0884923    0.5000000   0.9370000
  34   -0.0445121   1.0896288    0.5000000   0.9350000
  35   -0.0451943   1.0907661    0.5000000   0.9330000
  36   -0.0458664   1.0919042    0.5000000   0.9310000
  37   -0.0465289   1.0930431    0.5000000   0.9290000
  38   -0.0471828   1.0941831    0.5000000   0.9270000
  39   -0.0478264   1.0953236    0.5000000   0.9250000
  40   -0.0484623   1.0964651    0.5000000   0.9230000
  41   -0.0490899   1.0976075    0.5000000   0.9210000
  42   -0.0497096   1.0987507    0.5000000   0.9190000
  43   -0.0503217   1.0998948    0.5000000   0.9170000
  44   -0.0509258   1.1010396    0.5000000   0.9150000
  45   -0.0515242   1.1021857    0.5000000   0.9130000
  46   -0.0521150   1.1033324    0.5000000   0.9110000
  47   -0.0526992   1.1044801    0.5000000   0.9090000
  48   -0.0532770   1.1056286    0.5000000   0.9070000
  49   -0.0538486   1.1067781    0.5000000   0.9050000
  50   -0.0544148   1.1079287    0.5000000   0.9030000
  51   -0.0549738   1.1090798    0.5000000   0.9010000
  52   -0.0555278   1.1102320    0.5000000   0.8990000
  53   -0.0560763   1.1113852    0.5000000   0.8970000
  54   -0.0566195   1.1125393    0.5000000   0.8950000
  55   -0.0571574   1.1136944    0.5000000   0.8930000
  56   -0.0576902   1.1148503    0.5000000   0.8910000
  57   -0.0582182   1.1160073    0.5000000   0.8890000
  58   -0.0587413   1.1171652    0.5000000   0.8870000
  59   -0.0592598   1.1183241    0.5000000   0.8850000
  60   -0.0597738   1.1194839    0.5000000   0.8830000
  61   -0.0602834   1.1206447    0.5000000   0.8810000
  62   -0.0607890   1.1218065    0.5000000   0.8790000
  63   -0.0612888   1.1229694    0.5000000   0.8770000
  64   -0.0617854   1.1241332    0.5000000   0.8750000
  65   -0.0622780   1.1252980    0.5000000   0.8730000
  66   -0.0627666   1.1264638    0.5000000   0.8710000
  67   -0.0632513   1.1276306    0.5000000   0.8690000
  68   -0.0637312   1.1287985    0.5000000   0.8670000
  69   -0.0642094   1.1299673    0.5000000   0.8650000
  70   -0.0646830   1.1311372    0.5000000   0.8630000
  71   -0.0651530   1.1323082    0.5000000   0.8610000
  72   -0.0656195   1.1334802    0.5000000   0.8590000
  73   -0.0660825   1.1346532    0.5000000   0.8570000
  74   -0.0665432   1.1358272    0.5000000   0.8550000
  75   -0.0669987   1.1370024    0.5000000   0.8530000
  76   -0.0674518   1.1381786    0.5000000   0.8510000
  77   -0.0679018   1.1393559    0.5000000   0.8490000
  78   -0.0683487   1.1405343    0.5000000   0.8470000
  79   -0.0687926   1.1417138    0.5000000   0.8450000
  80   -0.0692324   1.1428944    0.5000000   0.8430000
  81   -0.0696713   1.1440759    0.5000000   0.8410000
  82   -0.0698891   1.1446672    0.5000000   0.8400000
  83   -0.0709677   1.1476275    0.5000000   0.8350000
  84   -0.0720289   1.1505949    0.5000000   0.8300000
  85   -0.0730738   1.1535693    0.5000000   0.8250000
  86   -0.0741024   1.1565510    0.5000000   0.8200000
  87   -0.0751157   1.1595399    0.5000000   0.8150000
  88   -0.0761143   1.1625362    0.5000000   0.8100000
  89   -0.0770985   1.1655401    0.5000000   0.8050000
  90   -0.0780689   1.1685514    0.5000000   0.8000000
  91   -0.0790258   1.1715706    0.5000000   0.7950000
  92   -0.0799699   1.1745974    0.5000000   0.7900000
  93   -0.0809014   1.1776321    0.5000000   0.7850000
  94   -0.0818206   1.1806748    0.5000000   0.7800000
  95   -0.0827279   1.1837256    0.5000000   0.7750000
  96   -0.0836236   1.1867845    0.5000000   0.7700000
  97   -0.0845082   1.1898517    0.5000000   0.7650000
  98   -0.0853816   1.1929273    0.5000000   0.7600000
  99   -0.0862444   1.1960114    0.5000000   0.7550000
 100   -0.0870968   1.1991041    0.5000000   0.7500000
 101   -0.0879389   1.2022054    0.5000000   0.7450000
 102   -0.0887711   1.2053155    0.5000000   0.7400000
 103   -0.0895934   1.2084345    0.5000000   0.7350000
 104   -0.0904064   1.2115624    0.5000000   0.7300000
 105   -0.0912099   1.2146995    0.5000000   0.7250000
 106   -0.0920042   1.2178458    0.5000000   0.7200000
 107   -0.0935661   1.2241663    0.5000000   0.7100000
 108   -0.0950933   1.2305249    0.5000000   0.7000000
 109   -0.0965871   1.2369225    0.5000000   0.6900000
 110   -0.0980485   1.2433599    0.5000000   0.6800000
 111   -0.0994787   1.2498380    0.5000000   0.6700000
 112   -0.1008790   1.2563577    0.5000000   0.6600000
 113   -0.1022474   1.2629201    0.5000000   0.6500000
 114   -0.1035882   1.2695259    0.5000000   0.6400000
 115   -0.1049006   1.2761761    0.5000000   0.6300000
 116   -0.1061855   1.2828717    0.5000000   0.6200000
 117   -0.1074432   1.2896137    0.5000000   0.6100000
 118   -0.1086732   1.2964032    0.5000000   0.6000000
 119   -0.1098796   1.3032409    0.5000000   0.5900000
 120   -0.1110591   1.3101281    0.5000000   0.5800000
 121   -0.1122135   1.3170658    0.5000000   0.5700000
 122   -0.1133429   1.3240550    0.5000000   0.5600000
 123   -0.1144479   1.3310968    0.5000000   0.5500000
 124   -0.1155300   1.3381924    0.5000000   0.5400000
 125   -0.1165856   1.3453429    0.5000000   0.5300000
 126   -0.1176189   1.3525494    0.5000000   0.5200000
 127   -0.1196155   1.3671351    0.5000000   0.5000000
 128   -0.1215201   1.3819593    0.5000000   0.4800000
 129   -0.1233338   1.3970321    0.5000000   0.4600000
 130   -0.1250580   1.4123639    0.5000000   0.4400000
 131   -0.1266935   1.4279656    0.5000000   0.4200000
 132   -0.1282409   1.4438486    0.5000000   0.4000000
 133   -0.1297006   1.4600249    0.5000000   0.3800000
 134   -0.1310729   1.4765068    0.5000000   0.3600000
 135   -0.1323581   1.4933075    0.5000000   0.3400000
 136   -0.1335557   1.5104406    0.5000000   0.3200000
 137   -0.1346660   1.5279203    0.5000000   0.3000000
 138   -0.1356886   1.5457617    0.5000000   0.2800000
 139   -0.1366232   1.5639806    0.5000000   0.2600000
 140   -0.1374693   1.5825935    0.5000000   0.2400000
 141   -0.1382262   1.6016178    0.5000000   0.2200000
 142   -0.1388937   1.6210719    0.5000000   0.2000000
 143   -0.1394707   1.6409752    0.5000000   0.1800000
 144   -0.1399566   1.6613480    0.5000000   0.1600000
 145   -0.1403505   1.6822119    0.5000000   0.1400000
 146   -0.1406515   1.7035896    0.5000000   0.1200000
 147   -0.1408588   1.7255053    0.5000000   0.1000000
 148   -0.1409711   1.7479844    0.5000000   0.0800000
 149   -0.1409871   1.7710539    0.5000000   0.0600000
 150   -0.1409055   1.7947423    0.5000000   0.0400000
 151   -0.1407288   1.8190812    0.5000000   0.0200000
 152   -0.1404511   1.8441018    0.5000000   0.0000000
```


### 7.2 IV.B Rectilinear problem (pp32 to 47)

READ. By studying the continuation of Stromgren's orbits it is seen that e = 1 "cannot be avoided in the natural prolongation of the families" and that the rectilinear problem plays a special role
(Schubart had proposed e = +1 as the start of a systematic study). Thirteen isolated periodic orbits were found by integrating regularly spaced initial conditions and then correcting; the
integrator is classical Runge-Kutta with variable step 0.005 r1 r2 (about five to six place accuracy); the corrections are the two-dimensional linear ones. Initial conditions have the form (x0, 0, 0, y_dot0) (Eq. 175),
restricted to 0.1 < x0 < 1.0, 0.2 < y_dot0 < 2.2 (Eq. 176) and y_dot0 below the parabolic limit y_dot0 = (2/x0)^(1/2) (Eq. 177) (the shaded region of Fig. 9, p39), swept with step 0.05 in each. The 13 points of Fig. 9 are the orbits found;
eleven of them (3 to 13) form a single sequence, each having one loop around one primary and an increasing number of loops around the other. All 13 have period T = 2 pi. e and mu: these orbits are isolated, but families arise near
e = 1, mu = 0.5. The four short segments: 5P near orbit 1 (22 orbits, e = 1, mu 0.5 to 0.479, Table 8); 9P near orbit 2 (12 orbits, mu = 0.5, e 1.0 to 0.988, Table 9); 4P near orbit 2 (26 orbits, e = 1,
mu 0.5 to 0.458, Table 10); 3P near orbit 3 (23 orbits, e = 1, mu 0.5 to 0.466, Table 11); Fig. 12 pp43 to 45. Figs. 10 and 11 (pp40 to 41) show the 13 orbits in barycentric (inertial) and in geocentric and selenocentric frames; cusps in
the latter come from the accelerated translation of the frame.

**Table 6, initial conditions for 13 periodic orbits, printed p38** (six decimals; compare with Table 7 below: the values are TRUNCATED, not rounded, e.g. 0.3962118 appears as 0.396211, and orbits 1 and 2 are numbered the other way round from Table 7).

```
Orbit        x0         ydot0
    1   0.579152     1.680050
    2   0.574506     0.346101
    3   0.682309     0.981417
    4   0.492110     1.192088
    5   0.396211     1.347039
    6   0.336523     1.473383
    7   0.295148     1.581601
    8   0.264484     1.677079
    9   0.240690     1.763014
   10   0.221596     1.841481
   11   0.205876     1.913916
   12   0.192667     1.981358
   13   0.171617     2.104198
```


**Table 7, initial and final conditions for the 13 periodic orbits, printed p42.** Header "RECTILINEAR ELLIPTIC RESTRICTED THREE-BODY PROBLEM 1 HALF REV. PERI.". Periapsis start (the primaries collide at t = 0), state after half a revolution (v = pi equivalent, t = pi),
X1 and YDOT1. NR 1 here is orbit 2 of Table 6, NR 2 is orbit 1.

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.5745062  0.3461014   2.2888544   0.2107393   0.5000000   1.0000000
   2  0.5791522  1.6800505  -1.6589334  -1.0153199   0.5000000   1.0000000
   3  0.6823097  0.9814177  -1.6020894   0.9510169   0.5000000   1.0000000
   4  0.4921103  1.1920888  -0.5101565  -0.9672414   0.5000000   1.0000000
   5  0.3962118  1.3470394  -1.2921322   1.4104449   0.5000000   1.0000000
   6  0.3365234  1.4733837  -0.6546279  -1.0962909   0.5000000   1.0000000
   7  0.2951483  1.5816018  -1.2049788   1.7100109   0.5000000   1.0000000
   8  0.2644841  1.6770793  -0.7239846  -1.2009274   0.5000000   1.0000000
   9  0.2406901  1.7630140  -1.1620703   1.9393112   0.5000000   1.0000000
  10  0.2215965  1.8414813  -0.7662645  -1.2895975   0.5000000   1.0000000
  11  0.2058762  1.9139168  -1.1360074   2.1282105   0.5000000   1.0000000
  12  0.1926678  1.9813582  -0.7952567  -1.3671887   0.5000000   1.0000000
  13  0.1716180  2.1041985  -0.8166100  -1.4366001   0.5000000   1.0000000
```


**Table 8, family 5P, printed p42.** Header "RECT.ELL.PROBLEM.MU=VARIABLE,E=1.0, FAMILY 5 (PERIAPSIS)".

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.5791522  1.6800505  -1.6589334  -1.0153199   0.5000000   1.0000000
   2  0.5798281  1.6787450  -1.6577717  -1.0153041   0.4990000   1.0000000
   3  0.5805048  1.6774394  -1.6566089  -1.0152883   0.4980000   1.0000000
   4  0.5811814  1.6761353  -1.6554461  -1.0152713   0.4970000   1.0000000
   5  0.5818578  1.6748327  -1.6542834  -1.0152530   0.4960000   1.0000000
   6  0.5825344  1.6735311  -1.6531205  -1.0152337   0.4950000   1.0000000
   7  0.5832109  1.6722311  -1.6519577  -1.0152131   0.4940000   1.0000000
   8  0.5838874  1.6709325  -1.6507949  -1.0151914   0.4930000   1.0000000
   9  0.5845639  1.6696351  -1.6496321  -1.0151686   0.4920000   1.0000000
  10  0.5852403  1.6683392  -1.6484693  -1.0151446   0.4910000   1.0000000
  11  0.5859166  1.6670446  -1.6473065  -1.0151194   0.4900000   1.0000000
  12  0.5865931  1.6657512  -1.6461436  -1.0150931   0.4890000   1.0000000
  13  0.5872694  1.6644593  -1.6449808  -1.0150656   0.4880000   1.0000000
  14  0.5879457  1.6631688  -1.6438180  -1.0150370   0.4870000   1.0000000
  15  0.5886221  1.6618792  -1.6426551  -1.0150074   0.4860000   1.0000000
  16  0.5892983  1.6605914  -1.6414924  -1.0149764   0.4850000   1.0000000
  17  0.5899745  1.6593047  -1.6403296  -1.0149444   0.4840000   1.0000000
  18  0.5906507  1.6580193  -1.6391668  -1.0149113   0.4830000   1.0000000
  19  0.5913268  1.6567353  -1.6380040  -1.0148770   0.4820000   1.0000000
  20  0.5920030  1.6554525  -1.6368412  -1.0148416   0.4810000   1.0000000
  21  0.5926790  1.6541710  -1.6356785  -1.0148051   0.4800000   1.0000000
  22  0.5933551  1.6528909  -1.6345158  -1.0147675   0.4790000   1.0000000
```


**Table 9, family 9P, printed p46.** Header "ELL.RESTR.PROBLEM. C E VARIABLE FROM 1.0 MU=0.5 PERI. 1 HALF REV. FAM 9" (so "C E" printed). Inertial coordinates, periapsis start, mu = 0.5, e from 1 to 0.988. Row 1 repeats Table 7 NR 1.

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.5745062  0.3461014   2.2888544   0.2107393   0.5000000   1.0000000
   2  0.5554261  0.3536558   2.2735915   0.2485901   0.5000000   0.9990000
   3  0.5391346  0.3588893   2.2672086   0.2649357   0.5000000   0.9980000
   4  0.5226942  0.3634967   2.2623948   0.2776792   0.5000000   0.9970000
   5  0.5057545  0.3677720   2.2584234   0.2885861   0.5000000   0.9960000
   6  0.4880889  0.3718801   2.2550035   0.2983681   0.5000000   0.9950000
   7  0.4694702  0.3759493   2.2519789   0.3074152   0.5000000   0.9940000
   8  0.4496181  0.3801073   2.2492501   0.3159891   0.5000000   0.9930000
   9  0.4281412  0.3845098   2.2467439   0.3243015   0.5000000   0.9920000
  10  0.4044309  0.3893839   2.2443955   0.3325661   0.5000000   0.9910000
  11  0.3774115  0.3951273   2.2421329   0.3410634   0.5000000   0.9900000
  12  0.2984375  0.4150494   2.2371638   0.3619545   0.5000000   0.9880000
```


**Table 10, family 4P, printed p46.** Header "RECT.ELL.PROBLEM.MU=VARIABLE,E=1.0, FAMILY 4 (PERIAPSIS)".

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.5745062  0.3461014   2.2888544   0.2107393   0.5000000   1.0000000
   2  0.5674326  0.3542697   2.2889281   0.2144953   0.4990000   1.0000000
   3  0.5603530  0.3624437   2.2889941   0.2181885   0.4980000   1.0000000
   4  0.5532670  0.3706281   2.2890527   0.2218205   0.4970000   1.0000000
   5  0.5461754  0.3788288   2.2891041   0.2253941   0.4960000   1.0000000
   6  0.5390791  0.3870518   2.2891482   0.2289118   0.4950000   1.0000000
   7  0.5319790  0.3953027   2.2891852   0.2323760   0.4940000   1.0000000
   8  0.5248759  0.4035872   2.2892151   0.2357887   0.4930000   1.0000000
   9  0.5177707  0.4119105   2.2892379   0.2391520   0.4920000   1.0000000
  10  0.5035575  0.4286952   2.2892623   0.2457375   0.4900000   1.0000000
  11  0.4893471  0.4456986   2.2892590   0.2521460   0.4880000   1.0000000
  12  0.4751470  0.4629622   2.2892282   0.2583889   0.4860000   1.0000000
  13  0.4609653  0.4805276   2.2891704   0.2644764   0.4840000   1.0000000
  14  0.4468103  0.4984367   2.2890859   0.2704176   0.4820000   1.0000000
  15  0.4326905  0.5167328   2.2889751   0.2762204   0.4800000   1.0000000
  16  0.4186150  0.5354604   2.2888385   0.2818918   0.4780000   1.0000000
  17  0.4045927  0.5546663   2.2886764   0.2874382   0.4760000   1.0000000
  18  0.3906333  0.5743994   2.2884891   0.2928654   0.4740000   1.0000000
  19  0.3767464  0.5947116   2.2882772   0.2981786   0.4720000   1.0000000
  20  0.3629422  0.6156580   2.2880409   0.3033823   0.4700000   1.0000000
  21  0.3492310  0.6372978   2.2877808   0.3084808   0.4680000   1.0000000
  22  0.3356234  0.6596942   2.2874972   0.3134780   0.4660000   1.0000000
  23  0.3221304  0.6829159   2.2871903   0.3183777   0.4640000   1.0000000
  24  0.3087633  0.7070368   2.2868607   0.3231828   0.4620000   1.0000000
  25  0.2955335  0.7321377   2.2865088   0.3278966   0.4600000   1.0000000
  26  0.2824529  0.7583063   2.2861349   0.3325218   0.4580000   1.0000000
```


**Table 11, family 3P, printed p47.** Header "RECT.ELL.PROBLEM.MU=VARIABLE,E=1.0, FAMILY 3 (PERIAPSIS)". Rows 9 and 10 are printed out of order (mass ratio 0.491 before 0.492); transcribed as printed.

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.6823097  0.9814177  -1.6020894   0.9510169   0.5000000   1.0000000
   2  0.6836911  0.9788757  -1.6002956   0.9519198   0.4990000   1.0000000
   3  0.6850722  0.9763395  -1.5985012   0.9528221   0.4980000   1.0000000
   4  0.6864528  0.9738089  -1.5967064   0.9537237   0.4970000   1.0000000
   5  0.6878332  0.9712840  -1.5949112   0.9546246   0.4960000   1.0000000
   6  0.6892132  0.9687647  -1.5931155   0.9555248   0.4950000   1.0000000
   7  0.6905929  0.9662509  -1.5913194   0.9564243   0.4940000   1.0000000
   8  0.6919722  0.9637427  -1.5895228   0.9573231   0.4930000   1.0000000
   9  0.6947297  0.9587429  -1.5859284   0.9591187   0.4910000   1.0000000
  10  0.6933511  0.9612401  -1.5877258   0.9582213   0.4920000   1.0000000
  11  0.6961080  0.9562512  -1.5841305   0.9600155   0.4900000   1.0000000
  12  0.6988633  0.9512840  -1.5805335   0.9618069   0.4880000   1.0000000
  13  0.7016172  0.9463383  -1.5769350   0.9635956   0.4860000   1.0000000
  14  0.7043695  0.9414138  -1.5733348   0.9653814   0.4840000   1.0000000
  15  0.7071203  0.9365104  -1.5697331   0.9671644   0.4820000   1.0000000
  16  0.7098694  0.9316278  -1.5661298   0.9689446   0.4800000   1.0000000
  17  0.7126169  0.9267658  -1.5625252   0.9707218   0.4780000   1.0000000
  18  0.7153627  0.9219241  -1.5589190   0.9724963   0.4760000   1.0000000
  19  0.7181068  0.9171026  -1.5553115   0.9742678   0.4740000   1.0000000
  20  0.7208490  0.9123010  -1.5517027   0.9760364   0.4720000   1.0000000
  21  0.7235895  0.9075192  -1.5480925   0.9778020   0.4700000   1.0000000
  22  0.7263281  0.9027568  -1.5444811   0.9795647   0.4680000   1.0000000
  23  0.7290649  0.8980138  -1.5408685   0.9813245   0.4660000   1.0000000
```


**Convention of the mass column in Tables 8 to 11 (COMPUTED).** Integrating the rectilinear problem in inertial coordinates (satellite under both primaries, primaries on their degenerate Kepler orbits r = 1 - cos E, t = E - sin E, collision at E = 0) reproduces the printed
(X1, YDOT1) of Tables 7, 8, 10, 11 only if the SMALLER primary (mass mu) is on the +x side at distance (1 - mu) r from the barycentre and the larger (1 - mu) is on the -x side at mu r (the Eq. 2 arrangement
at v = 0). Eq. 36c,d (p7) prints xi1 = +m2 r > 0, xi2 = -m1 r < 0 (the larger at +mu r): that is the mirror image. The two differ by a reflection xi -> -xi, which is a symmetry of the rectilinear problem, so the physics is unchanged, but the sign of X0 flips. In
Table 4 (apoapsis start, E from pi) X0 < 0 matches the Eq. 36 arrangement. In both cases the satellite starts on the same side as the SMALLER primary.

### 7.3 IV.C Family 7P (pp39 to 51)

READ. Earth-Moon mass ratio mu = 0.012155, e variable. Begins with an orbit of "family C of retrograde satellite orbits around the smaller primary m2 in the circular problem" described in Ref. 6 (Broucke 1968), the one close to its orbit
87 because that has period 2 pi. Initial conditions in rotating axes: **x0 = 0.15212027, y_dot0 = 3.16076559, mu = 0.012155, e = 0.0** (Eq. 178a to d). About 130 orbits, e from 0 to 0.50; no higher e because there is a collision with the larger primary just above 0.50
(the orbits come closer to m1 than to m2 although they are "satellite orbits" of the circular problem; Fig. 13, p48). **Stability (text, p39 and p51):** good information only for e up to 0.35; all orbits are unstable, belong to region 6, and lie very close to the line a2 = -2 a1 - 2 ("essentially two-body orbits
about the larger primary, which appear in inertial axes as ellipses weakly perturbed by the smaller primary; there is still approximately an integral and a pair of eigenvalues close to +1"). In the circular problem all such orbits lie on a2 = -2 a1 - 2.

**Table 12, family 7P, printed pp49 to 51** (131 rows; rows 1 to 60 p49, 61 to 120 p50, 121 to 131 p51). Header "ELLIPTIC PROBLEM MU=0.012155 PS 1 HALF REV PERIAPSIS FAMILY 7". Barycentric rotating-pulsating coordinates, periapsis start (v = 0), final state at one half
revolution (v = pi). The state is (X0, 0, 0, YDOT0); X1 and YDOT1 are the x and y_dot at v = pi (y and x_dot there vanish).

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.1520965  3.1608994   1.8428499  -1.5495568   0.0121550   0.0001000
   2  0.1520728  3.1610334   1.8427078  -1.5494262   0.0121550   0.0002000
   3  0.1520491  3.1611673   1.8425657  -1.5492956   0.0121550   0.0003000
   4  0.1520253  3.1613014   1.8424237  -1.5491650   0.0121550   0.0004000
   5  0.1520059  3.1613835   1.8422775  -1.5490262   0.0121550   0.0005000
   6  0.1519822  3.1615176   1.8421355  -1.5488957   0.0121550   0.0006000
   7  0.1519584  3.1616519   1.8419936  -1.5487653   0.0121550   0.0007000
   8  0.1519347  3.1617862   1.8418517  -1.5486348   0.0121550   0.0008000
   9  0.1519109  3.1619206   1.8417098  -1.5485044   0.0121550   0.0009000
  10  0.1518872  3.1620550   1.8415679  -1.5483741   0.0121550   0.0010000
  11  0.1518634  3.1621895   1.8414260  -1.5482437   0.0121550   0.0011000
  12  0.1518354  3.1623755   1.8412882  -1.5481215   0.0121550   0.0012000
  13  0.1517879  3.1626453   1.8410047  -1.5478611   0.0121550   0.0014000
  14  0.1517404  3.1629148   1.8407213  -1.5476007   0.0121550   0.0016000
  15  0.1516928  3.1631846   1.8404379  -1.5473404   0.0121550   0.0018000
  16  0.1516453  3.1634547   1.8401546  -1.5470802   0.0121550   0.0020000
  17  0.1515978  3.1637250   1.8398714  -1.5468201   0.0121550   0.0022000
  18  0.1515503  3.1639956   1.8395883  -1.5465602   0.0121550   0.0024000
  19  0.1515028  3.1642664   1.8393054  -1.5463004   0.0121550   0.0026000
  20  0.1514552  3.1645375   1.8390225  -1.5460407   0.0121550   0.0028000
  21  0.1514077  3.1648088   1.8387397  -1.5457812   0.0121550   0.0030000
  22  0.1513602  3.1650804   1.8384570  -1.5455218   0.0121550   0.0032000
  23  0.1513126  3.1653523   1.8381745  -1.5452625   0.0121550   0.0034000
  24  0.1512651  3.1656244   1.8378920  -1.5450033   0.0121550   0.0036000
  25  0.1512175  3.1658967   1.8376096  -1.5447442   0.0121550   0.0038000
  26  0.1511700  3.1661693   1.8373273  -1.5444853   0.0121550   0.0040000
  27  0.1511224  3.1664422   1.8370452  -1.5442264   0.0121550   0.0042000
  28  0.1510748  3.1667153   1.8367631  -1.5439678   0.0121550   0.0044000
  29  0.1510272  3.1669887   1.8364811  -1.5437092   0.0121550   0.0046000
  30  0.1509797  3.1672623   1.8361993  -1.5434507   0.0121550   0.0048000
  31  0.1509321  3.1675362   1.8359175  -1.5431924   0.0121550   0.0050000
  32  0.1508845  3.1678104   1.8356358  -1.5429342   0.0121550   0.0052000
  33  0.1508369  3.1680847   1.8353543  -1.5426761   0.0121550   0.0054000
  34  0.1508131  3.1682220   1.8352135  -1.5425471   0.0121550   0.0055000
  35  0.1507893  3.1683594   1.8350728  -1.5424181   0.0121550   0.0056000
  36  0.1506941  3.1689094   1.8345102  -1.5419026   0.0121550   0.0060000
  37  0.1505751  3.1695985   1.8338074  -1.5412588   0.0121550   0.0065000
  38  0.1504560  3.1702892   1.8331053  -1.5406158   0.0121550   0.0070000
  39  0.1503370  3.1709814   1.8324038  -1.5399736   0.0121550   0.0075000
  40  0.1502179  3.1716753   1.8317030  -1.5393321   0.0121550   0.0080000
  41  0.1500987  3.1723708   1.8310027  -1.5386913   0.0121550   0.0085000
  42  0.1499795  3.1730680   1.8303031  -1.5380514   0.0121550   0.0090000
  43  0.1498604  3.1737667   1.8296041  -1.5374122   0.0121550   0.0095000
  44  0.1497411  3.1744671   1.8289057  -1.5367737   0.0121550   0.0100000
  45  0.1495029  3.1758696   1.8275105  -1.5354985   0.0121550   0.0110000
  46  0.1492642  3.1772816   1.8261181  -1.5342269   0.0121550   0.0120000
  47  0.1490255  3.1787003   1.8247281  -1.5329582   0.0121550   0.0130000
  48  0.1487866  3.1801254   1.8233405  -1.5316926   0.0121550   0.0140000
  49  0.1485477  3.1815572   1.8219554  -1.5304299   0.0121550   0.0150000
  50  0.1483086  3.1829956   1.8205727  -1.5291703   0.0121550   0.0160000
  51  0.1480694  3.1844405   1.8191925  -1.5279136   0.0121550   0.0170000
  52  0.1478301  3.1858921   1.8178146  -1.5266599   0.0121550   0.0180000
  53  0.1475907  3.1873503   1.8164392  -1.5254092   0.0121550   0.0190000
  54  0.1473511  3.1888152   1.8150662  -1.5241614   0.0121550   0.0200000
  55  0.1471115  3.1902868   1.8136956  -1.5229166   0.0121550   0.0210000
  56  0.1468718  3.1917651   1.8123273  -1.5216747   0.0121550   0.0220000
  57  0.1466319  3.1932502   1.8109615  -1.5204358   0.0121550   0.0230000
  58  0.1463920  3.1947419   1.8095980  -1.5191998   0.0121550   0.0240000
  59  0.1461519  3.1962405   1.8082369  -1.5179667   0.0121550   0.0250000
  60  0.1459117  3.1977458   1.8068781  -1.5167366   0.0121550   0.0260000
  61  0.1456714  3.1992580   1.8055217  -1.5155093   0.0121550   0.0270000
  62  0.1454310  3.2007770   1.8041677  -1.5142850   0.0121550   0.0280000
  63  0.1451905  3.2023028   1.8028160  -1.5130636   0.0121550   0.0290000
  64  0.1449498  3.2038355   1.8014666  -1.5118451   0.0121550   0.0300000
  65  0.1437451  3.2116029   1.7947543  -1.5057956   0.0121550   0.0350000
  66  0.1425372  3.2195499   1.7881000  -1.4998187   0.0121550   0.0400000
  67  0.1413268  3.2276708   1.7815017  -1.4939117   0.0121550   0.0450000
  68  0.1401138  3.2359716   1.7749588  -1.4880748   0.0121550   0.0500000
  69  0.1388980  3.2444577   1.7684709  -1.4823074   0.0121550   0.0550000
  70  0.1376795  3.2531297   1.7620370  -1.4766087   0.0121550   0.0600000
  71  0.1364583  3.2619923   1.7556563  -1.4709780   0.0121550   0.0650000
  72  0.1352343  3.2710488   1.7493281  -1.4654149   0.0121550   0.0700000
  73  0.1340076  3.2803021   1.7430516  -1.4599185   0.0121550   0.0750000
  74  0.1327783  3.2897563   1.7368262  -1.4544885   0.0121550   0.0800000
  75  0.1315462  3.2994151   1.7306510  -1.4491240   0.0121550   0.0850000
  76  0.1303114  3.3092823   1.7245254  -1.4438247   0.0121550   0.0900000
  77  0.1290740  3.3193620   1.7184486  -1.4385899   0.0121550   0.0950000
  78  0.1278339  3.3296584   1.7124201  -1.4334192   0.0121550   0.1000000
  79  0.1265911  3.3401755   1.7064393  -1.4283121   0.0121550   0.1050000
  80  0.1253457  3.3509190   1.7005051  -1.4232678   0.0121550   0.1100000
  81  0.1240976  3.3618924   1.6946172  -1.4182861   0.0121550   0.1150000
  82  0.1228469  3.3731008   1.6887750  -1.4133666   0.0121550   0.1200000
  83  0.1215936  3.3845496   1.6829778  -1.4085087   0.0121550   0.1250000
  84  0.1203377  3.3962437   1.6772250  -1.4037120   0.0121550   0.1300000
  85  0.1190792  3.4081887   1.6715159  -1.3989761   0.0121550   0.1350000
  86  0.1178181  3.4203902   1.6658501  -1.3943007   0.0121550   0.1400000
  87  0.1165544  3.4328541   1.6602269  -1.3896852   0.0121550   0.1450000
  88  0.1152882  3.4455858   1.6546459  -1.3851295   0.0121550   0.1500000
  89  0.1140193  3.4585936   1.6491062  -1.3806330   0.0121550   0.1550000
  90  0.1127480  3.4718822   1.6436076  -1.3761955   0.0121550   0.1600000
  91  0.1114742  3.4854586   1.6381494  -1.3718167   0.0121550   0.1650000
  92  0.1101979  3.4993308   1.6327311  -1.3674963   0.0121550   0.1700000
  93  0.1089191  3.5135053   1.6273523  -1.3632339   0.0121550   0.1750000
  94  0.1076379  3.5279901   1.6220123  -1.3590293   0.0121550   0.1800000
  95  0.1063542  3.5427935   1.6167108  -1.3548823   0.0121550   0.1850000
  96  0.1050681  3.5579232   1.6114471  -1.3507926   0.0121550   0.1900000
  97  0.1037797  3.5733881   1.6062209  -1.3467599   0.0121550   0.1950000
  98  0.1024888  3.5891975   1.6010317  -1.3427842   0.0121550   0.2000000
  99  0.1011956  3.6053602   1.5958790  -1.3388650   0.0121550   0.2050000
 100  0.0999002  3.6218855   1.5907624  -1.3350025   0.0121550   0.2100000
 101  0.0986024  3.6387860   1.5856812  -1.3311961   0.0121550   0.2150000
 102  0.0973023  3.6560701   1.5806353  -1.3274460   0.0121550   0.2200000
 103  0.0960001  3.6737485   1.5756241  -1.3237518   0.0121550   0.2250000
 104  0.0946957  3.6918346   1.5706472  -1.3201137   0.0121550   0.2300000
 105  0.0920803  3.7292754   1.5607948  -1.3130047   0.0121550   0.2400000
 106  0.0894569  3.7684904   1.5510745  -1.3061179   0.0121550   0.2500000
 107  0.0868251  3.8096025   1.5414837  -1.2994542   0.0121550   0.2600000
 108  0.0841858  3.8527260   1.5320190  -1.2930122   0.0121550   0.2700000
 109  0.0815395  3.8979911   1.5226774  -1.2867919   0.0121550   0.2800000
 110  0.0788861  3.9455605   1.5134566  -1.2807952   0.0121550   0.2900000
 111  0.0762265  3.9955870   1.5043534  -1.2750217   0.0121550   0.3000000
 112  0.0735613  4.0482459   1.4953653  -1.2694722   0.0121550   0.3100000
 113  0.0708908  4.1037439   1.4864898  -1.2641485   0.0121550   0.3200000
 114  0.0682157  4.1622978   1.4777244  -1.2590522   0.0121550   0.3300000
 115  0.0655371  4.2241443   1.4690667  -1.2541846   0.0121550   0.3400000
 116  0.0628552  4.2895686   1.4605147  -1.2495490   0.0121550   0.3500000
 117  0.0601710  4.3588715   1.4520660  -1.2451478   0.0121550   0.3600000
 118  0.0574860  4.4323779   1.4437183  -1.2409829   0.0121550   0.3700000
 119  0.0548006  4.5104873   1.4354699  -1.2370589   0.0121550   0.3800000
 120  0.0521161  4.5936304   1.4273187  -1.2333795   0.0121550   0.3900000
 121  0.0494341  4.6822788   1.4192626  -1.2299480   0.0121550   0.4000000
 122  0.0467555  4.7770100   1.4113001  -1.2267703   0.0121550   0.4100000
 123  0.0440821  4.8784508   1.4034294  -1.2238513   0.0121550   0.4200000
 124  0.0414156  4.9873240   1.3956485  -1.2211963   0.0121550   0.4300000
 125  0.0387577  5.1044787   1.3879562  -1.2188122   0.0121550   0.4400000
 126  0.0361103  5.2308875   1.3803507  -1.2167060   0.0121550   0.4500000
 127  0.0334761  5.3676731   1.3728305  -1.2148847   0.0121550   0.4600000
 128  0.0308572  5.5161770   1.3653942  -1.2133572   0.0121550   0.4700000
 129  0.0282564  5.6779674   1.3580405  -1.2121326   0.0121550   0.4800000
 130  0.0256766  5.8549100   1.3507681  -1.2112207   0.0121550   0.4900000
 131  0.0231214  6.0492417   1.3435756  -1.2106323   0.0121550   0.5000000
```


### 7.4 IV.D Family 7A (p51 and pp52 to 55)

READ. Begins with the same orbit as 7P (apoapsis start now); about 120 orbits "with eccentricity e from 0.0 to 0.99 (Table 13)"; all computed in rotating-pulsating axes with the Nechville transformation; e = 1 would need inertial axes. Stability (p51): good information only for
e below 0.75; "all of the orbits are stable, and belong to region 1 (Fig. 2). They are on a line close to a2 = -2 a1 - 2. When e increases, the orbits approach the separation point (-4, 6)". (Note the text's "-2a2 - 2" is a typo for a2 = -2 a1 - 2.)

**Table 13, family 7A, printed pp54 to 55** (120 rows; rows 1 to 60 p54, 61 to 120 p55). Header "ELLIPTIC PROBLEM MU=0.012155 E VARIABLE APOAPSIS FAMILY 7". Apoapsis start (v = pi), final state after half a revolution (v = 2 pi). **The table ends at e = 0.93, not 0.99
as the text says;** Fig. 14 (pp52 to 53) shows orbits at e = 0.95 and 0.975 which are not tabulated.

```
  NR         X0       YDOT0          X1       YDOT1  MASS RATIO         ECC
   1  0.1521203  3.1607656    1.8429920   -1.5496874   0.0121550   0.0000000
   2  0.1527133  3.1574396    1.8465530   -1.5529639   0.0121550   0.0025000
   3  0.1530096  3.1557914    1.8483396   -1.5546093   0.0121550   0.0037500
   4  0.1533057  3.1541530    1.8501301   -1.5562597   0.0121550   0.0050000
   5  0.1536016  3.1525243    1.8519247   -1.5579149   0.0121550   0.0062500
   6  0.1538973  3.1509053    1.8537234   -1.5595751   0.0121550   0.0075000
   7  0.1541929  3.1492960    1.8555261   -1.5612401   0.0121550   0.0087500
   8  0.1544883  3.1476964    1.8573329   -1.5629101   0.0121550   0.0100000
   9  0.1547835  3.1461064    1.8591438   -1.5645850   0.0121550   0.0112500
  10  0.1550785  3.1445259    1.8609589   -1.5662649   0.0121550   0.0125000
  11  0.1553734  3.1429551    1.8627780   -1.5679497   0.0121550   0.0137500
  12  0.1556681  3.1413938    1.8646014   -1.5696395   0.0121550   0.0150000
  13  0.1568450  3.1352432    1.8719365   -1.5764490   0.0121550   0.0200000
  14  0.1580192  3.1292428    1.8793394   -1.5833395   0.0121550   0.0250000
  15  0.1591906  3.1233909    1.8868111   -1.5903119   0.0121550   0.0300000
  16  0.1603591  3.1176858    1.8943528   -1.5973674   0.0121550   0.0350000
  17  0.1615248  3.1121257    1.9019655   -1.6045069   0.0121550   0.0400000
  18  0.1626877  3.1067091    1.9096503   -1.6117315   0.0121550   0.0450000
  19  0.1638478  3.1014346    1.9174086   -1.6190425   0.0121550   0.0500000
  20  0.1650050  3.0963007    1.9252413   -1.6264407   0.0121550   0.0550000
  21  0.1661593  3.0913060    1.9331499   -1.6339276   0.0121550   0.0600000
  22  0.1673108  3.0864492    1.9411355   -1.6415041   0.0121550   0.0650000
  23  0.1684595  3.0817291    1.9491993   -1.6491717   0.0121550   0.0700000
  24  0.1696052  3.0771444    1.9573428   -1.6569314   0.0121550   0.0750000
  25  0.1707481  3.0726941    1.9655672   -1.6647846   0.0121550   0.0800000
  26  0.1718881  3.0683771    1.9738739   -1.6727326   0.0121550   0.0850000
  27  0.1730252  3.0641923    1.9822643   -1.6807768   0.0121550   0.0900000
  28  0.1741594  3.0601389    1.9907398   -1.6889186   0.0121550   0.0950000
  29  0.1752907  3.0562158    1.9993019   -1.6971593   0.0121550   0.1000000
  30  0.1764191  3.0524223    2.0079521   -1.7055004   0.0121550   0.1050000
  31  0.1775445  3.0487575    2.0166919   -1.7139435   0.0121550   0.1100000
  32  0.1786670  3.0452207    2.0255228   -1.7224900   0.0121550   0.1150000
  33  0.1797866  3.0418113    2.0344465   -1.7311415   0.0121550   0.1200000
  34  0.1809032  3.0385284    2.0434646   -1.7398996   0.0121550   0.1250000
  35  0.1820169  3.0353717    2.0525787   -1.7487660   0.0121550   0.1300000
  36  0.1831275  3.0323404    2.0617907   -1.7577424   0.0121550   0.1350000
  37  0.1842352  3.0294342    2.0711022   -1.7668305   0.0121550   0.1400000
  38  0.1853399  3.0266526    2.0805150   -1.7760321   0.0121550   0.1450000
  39  0.1864415  3.0239952    2.0900311   -1.7853489   0.0121550   0.1500000
  40  0.1875401  3.0214616    2.0996522   -1.7947830   0.0121550   0.1550000
  41  0.1886357  3.0190516    2.1093803   -1.8043361   0.0121550   0.1600000
  42  0.1908177  3.0146011    2.1291655   -1.8238077   0.0121550   0.1700000
  43  0.1929873  3.0106425    2.1494033   -1.8437800   0.0121550   0.1800000
  44  0.1951445  3.0071750    2.1701108   -1.8642704   0.0121550   0.1900000
  45  0.1972890  3.0041986    2.1913062   -1.8852970   0.0121550   0.2000000
  46  0.1994207  3.0017143    2.2130087   -1.9068790   0.0121550   0.2100000
  47  0.2015395  2.9997237    2.2352385   -1.9290368   0.0121550   0.2200000
  48  0.2036452  2.9982293    2.2580168   -1.9517915   0.0121550   0.2300000
  49  0.2057375  2.9972341    2.2813659   -1.9751659   0.0121550   0.2400000
  50  0.2078163  2.9967422    2.3053095   -1.9991837   0.0121550   0.2500000
  51  0.2098814  2.9967586    2.3298726   -2.0238701   0.0121550   0.2600000
  52  0.2119325  2.9972889    2.3550814   -2.0492517   0.0121550   0.2700000
  53  0.2139692  2.9983396    2.3809640   -2.0753567   0.0121550   0.2800000
  54  0.2159915  2.9999184    2.4075496   -2.1022151   0.0121550   0.2900000
  55  0.2179989  3.0020336    2.4348697   -2.1298583   0.0121550   0.3000000
  56  0.2199911  3.0046948    2.4629574   -2.1583200   0.0121550   0.3100000
  57  0.2219678  3.0079123    2.4918476   -2.1876358   0.0121550   0.3200000
  58  0.2239287  3.0116978    2.5215779   -2.2178435   0.0121550   0.3300000
  59  0.2258732  3.0160639    2.5521878   -2.2489833   0.0121550   0.3400000
  60  0.2278010  3.0210246    2.5837195   -2.2810981   0.0121550   0.3500000
  61  0.2297116  3.0265951    2.6162180   -2.3142334   0.0121550   0.3600000
  62  0.2316046  3.0327919    2.6497310   -2.3484378   0.0121550   0.3700000
  63  0.2334794  3.0396330    2.6843098   -2.3837631   0.0121550   0.3800000
  64  0.2353354  3.0471379    2.7200087   -2.4202647   0.0121550   0.3900000
  65  0.2371720  3.0553277    2.7568863   -2.4580018   0.0121550   0.4000000
  66  0.2389887  3.0642254    2.7950048   -2.4970377   0.0121550   0.4100000
  67  0.2407846  3.0738557    2.8344313   -2.5374404   0.0121550   0.4200000
  68  0.2425591  3.0842454    2.8752375   -2.5792827   0.0121550   0.4300000
  69  0.2443113  3.0954236    2.9175004   -2.6226427   0.0121550   0.4400000
  70  0.2460404  3.1074215    2.9613031   -2.6676046   0.0121550   0.4500000
  71  0.2477456  3.1202732    3.0067347   -2.7142589   0.0121550   0.4600000
  72  0.2494257  3.1340152    3.0538917   -2.7627032   0.0121550   0.4700000
  73  0.2510797  3.1486873    3.1028781   -2.8130430   0.0121550   0.4800000
  74  0.2527066  3.1643325    3.1538062   -2.8653922   0.0121550   0.4900000
  75  0.2543051  3.1809973    3.2067979   -2.9198743   0.0121550   0.5000000
  76  0.2558739  3.1987323    3.2619855   -2.9766229   0.0121550   0.5100000
  77  0.2574116  3.2175921    3.3195124   -3.0357835   0.0121550   0.5200000
  78  0.2589166  3.2376365    3.3795351   -3.0975142   0.0121550   0.5300000
  79  0.2603875  3.2589301    3.4422239   -3.1619875   0.0121550   0.5400000
  80  0.2618223  3.2815433    3.5077649   -3.2293915   0.0121550   0.5500000
  81  0.2632192  3.3055532    3.5763620   -3.2999322   0.0121550   0.5600000
  82  0.2645763  3.3310436    3.6482386   -3.3738353   0.0121550   0.5700000
  83  0.2658912  3.3581063    3.7236403   -3.4513489   0.0121550   0.5800000
  84  0.2671616  3.3868419    3.8028373   -3.5327460   0.0121550   0.5900000
  85  0.2683849  3.4173607    3.8861284   -3.6183279   0.0121550   0.6000000
  86  0.2695584  3.4497839    3.9738437   -3.7084276   0.0121550   0.6100000
  87  0.2706789  3.4842453    4.0663497   -3.8034149   0.0121550   0.6200000
  88  0.2717433  3.5208923    4.1640539   -3.9037004   0.0121550   0.6300000
  89  0.2727479  3.5598884    4.2674108   -4.0097420   0.0121550   0.6400000
  90  0.2736889  3.6014145    4.3769288   -4.1220520   0.0121550   0.6500000
  91  0.2745619  3.6456723    4.4931783   -4.2412046   0.0121550   0.6600000
  92  0.2753624  3.6928864    4.6168014   -4.3678461   0.0121550   0.6700000
  93  0.2760853  3.7433082    4.7485234   -4.5027062   0.0121550   0.6800000
  94  0.2767250  3.7972199    4.8891665   -4.6466119   0.0121550   0.6900000
  95  0.2772754  3.8549392    5.0396665   -4.8005040   0.0121550   0.7000000
  96  0.2777297  3.9168252    5.2010930   -4.9654574   0.0121550   0.7100000
  97  0.2780805  3.9832852    5.3746733   -5.1427053   0.0121550   0.7200000
  98  0.2783196  4.0547833    5.5618231   -5.3336694   0.0121550   0.7300000
  99  0.2784378  4.1318508    5.7641832   -5.5399970   0.0121550   0.7400000
 100  0.2784251  4.2150985    5.9836655   -5.7636074   0.0121550   0.7500000
 101  0.2782701  4.3052322    6.2225111   -6.0067493   0.0121550   0.7600000
 102  0.2779601  4.4030728    6.4833634   -6.2720740   0.0121550   0.7700000
 103  0.2774809  4.5095803    6.7693616   -6.5627298   0.0121550   0.7800000
 104  0.2768164  4.6258848    7.0842611   -6.8824819   0.0121550   0.7900000
 105  0.2759484  4.7533272    7.4325915   -7.2358700   0.0121550   0.8000000
 106  0.2748561  4.8935101    7.8198638   -7.6284169   0.0121550   0.8100000
 107  0.2735159  5.0483671    8.2528496   -8.0669066   0.0121550   0.8200000
 108  0.2719001  5.2202534    8.7399588   -8.5597626   0.0121550   0.8300000
 109  0.2699771  5.4120692    9.2917616   -9.1175708   0.0121550   0.8400000
 110  0.2677096  5.6274301    9.9217234   -9.7538131   0.0121550   0.8500000
 111  0.2650536  5.8709058   10.6472574  -10.4859223   0.0121550   0.8600000
 112  0.2619565  6.1463636   11.4912723  -11.3368286   0.0121550   0.8700000
 113  0.2583550  6.4674731   12.4845035  -12.3372918   0.0121550   0.8800000
 114  0.2541711  6.8384717   13.6691338  -13.5295236   0.0121550   0.8900000
 115  0.2493079  7.2753611   15.1046178  -14.9730126   0.0121550   0.9000000
 116  0.2436419  7.7978572   16.8774426  -16.7542861   0.0121550   0.9100000
 117  0.2404594  8.0997033   17.9287531  -17.8100023   0.0121550   0.9150000
 118  0.2370121  8.4347116   19.1183074  -19.0040927   0.0121550   0.9200000
 119  0.2332711  8.8088615   20.4745557  -20.3650155   0.0121550   0.9250000
 120  0.2292023  9.2297084   22.0342278  -21.9295106   0.0121550   0.9300000
```


## 8. Closure checks in `core.er3bp` (COMPUTED, scratch scripts, not committed)

Method: `propagate_er3bp` (DOP853, rtol = atol = 1e-13), state `[x, y, z, x', y', z']`, mu = 0.012155 for 7P and 7A, mu = 0.5 for 12A and 9P. Pulsating families: start
`[x0, 0, 0, 0, ydot0, 0]` at `f0 = 0` (P) or `f0 = pi` (A). Inertial families: converted with section 3.1. Rectilinear rows: a separate inertial integration (not `core.er3bp`; it cannot take e = 1).

**Convention mapping to `core.er3bp` (CONFIRMED by these closures):**
1. State and frame: the printed (x0, y0 = 0, x_dot0, y_dot0) of the pulsating families IS the `core.er3bp` state (x, y, x', y') with x' = dx/df. m1 (larger) at -mu, m2 at 1 - mu, mu = m2.
2. Start epoch: P families at f = 0, A families at f = pi. The final column (X1, YDOT1) is the state at f0 + pi (the half period, where y = x' = 0).
3. mu for the Earth-Moon tables is 0.012155 (not 0.012155099).
4. Inertial families (12A, 9P, 6A, 5P, 4P, 3P, 1P): inertial barycentric, x axis along the periapsis direction, dots d/dt; conversion in section 3.1; 12A closes only with the apoapsis conversion x = -xi0/r_a, y_dot = xi0/r_a - r_a eta'/p^(1/2) (the variant without the sign flip misses by about 1).

| Orbit | e | Check | Residuals |
|---|---|---|---|
| 7P row 1 (Table 12) | 0.0001 | half revolution from f = 0 to pi | y_end 4.1e-6, x'_end 9.6e-8; x1 miss -5.7e-6, y_dot1 miss +6.6e-6 (printed 1.8428499, -1.5495568) |
| 7P row 60 | 0.026 | same | y_end -2.3e-7, x'_end -5.1e-8; x1 miss 3.6e-7, y_dot1 miss -3.3e-7 |
| 7P row 102 | 0.22 | same | y_end 7.7e-6, x'_end -7.6e-6; x1 miss -1.5e-5, y_dot1 miss 1.7e-5 |
| 7P row 131 | 0.50 | same | y_end -3.6e-5, x'_end 2.7e-4; x1 miss 1.8e-4, y_dot1 miss -1.8e-4 |
| 7A row 1 (Table 13) | 0.0 | half revolution from f = pi to 2 pi | y_end -3.6e-6, x'_end -8.4e-8; x1 miss 4.9e-6, y_dot1 miss -5.7e-6 |
| 7A row 40 | 0.155 | same | y_end 3.2e-6, x'_end 1.2e-6; x1 miss -3.6e-6, y_dot1 miss 4.8e-6 |
| 7A row 75 | 0.50 | same | y_end 1.1e-6, x'_end 8.5e-7; x1 miss -8.3e-7, y_dot1 miss 9.9e-7 |
| Eq. 178 start (8 digits), e = 0 | 0.0 | full 2 pi closure; half revolution | full-period state difference 2.2e-5; half-revolution state 1.8429921, -1.5496876 against Table 13 row 1 (1.8429920, -1.5496874) |
| 12A row 152 (Table 5) | 0.0 | full 4 pi closure; half-period crossing | 1.1e-5; y_end 1.3e-6, x'_end 1.0e-6 |
| 12A row 127 | 0.5 | same | 8.2e-6; y_end 1.6e-6, x'_end 2.2e-6 |
| 12A row 90 | 0.8 | same | 4.5e-6; y_end 8.5e-7, x'_end 7.6e-7 |
| 12A row 60 | 0.883 | same | 4.6e-5; y_end 4.8e-7, x'_end 5.0e-6 |
| 9P rows 2 to 12 (Table 9) | 0.999 to 0.988 | periapsis start, inertial to pulsating, half revolution | x1 miss from 2.5e-6 (e = 0.999) to 6.0e-5 (e = 0.988), y_dot1 miss at most 7.3e-6, y_end at most 2.3e-5 |

**Whole-table checks (every row, from the printed state):**
- Table 12 (7P, 131 rows): relative miss of the printed (x1, y_dot1): median 4.2e-6, maximum 1.5e-4 (rows 131, 130, 126; the rows with the largest multipliers).
- Table 13 (7A, 120 rows): median 8.3e-7, maximum 2.4e-3 (row 112, e = 0.87; the printed orbit there is strongly sensitive); only that one row exceeds 1e-3.
- Table 5 (12A, 151 elliptic rows): residual of the perpendicular crossing at the half period (y and x'), median 8.5e-7, maximum 8.9e-5 (row 2, e = 0.999).
- Table 4 (6A, 118 rows), rectilinear inertial integration from apoapsis to the half period: median 8.4e-7, maximum 4.7e-4 (rows 114 to 117, mu about 0.17).
- Tables 7, 8, 10, 11 (rectilinear, own integration, half revolution from collision): maxima 1.4e-6, 1.0e-6, 3.3e-5 (row 26), 2.6e-7.
- Transcription slips (a dropped minus sign, a decimal point read as a digit) were found by exactly these checks; see section 12 for how the tables were transcribed.

Observation for `#930`: the printed seven-digit states close only to about 1e-5 over a full period; for rows with a large multiplier (7P row 131, e = 0.5) the full-period miss of the PRINTED state is 2.2, while a state corrected in two parameters
(x0 = 0.02312131, y_dot0 = 6.04924543, a change of 3.8e-6 in y_dot0) closes to 1.1e-10. The half-period test (y = x' = 0 at the second apse) is the right test of a printed orbit; the full-period test of a printed state is not.

## 9. Stability recomputed for 7P and 7A (COMPUTED)

Method: for each row, correct (x0, y_dot0) by Newton iteration (fsolve) so that y = x' = 0 at v0 + pi, then integrate the planar variational equations over 2 pi in eight segments (`er3bp_stm_eom`), take the planar 4 x 4 block, a1 = -trace, a2 = the sum of the six principal minors, k = roots of Eq. 166, region by Table 3.

| Family | Row (e) | a1 | a2 | k1, k2 | Region (re-corrected) | Largest planar modulus | Region from the printed seven-digit state |
|---|---|---|---|---|---|---|---|
| 7P | 1 (0.0001) | -2.9416 | 3.8832 | 2.00000, 0.94162 | 6 | 1.0008 | 6 |
| 7P | 30 (0.0048) | -2.9317 | 3.8633 | 2.00003, 0.93164 | 6 | 1.0056 | 1 |
| 7P | 60 (0.026) | -2.8848 | 3.7694 | 2.00018, 0.8846 | 6 | 1.0134 | 1 |
| 7P | 78 (0.100) | -2.6909 | 3.3806 | 2.00085, 0.69 | 6 | 1.0297 | 1 |
| 7P | 100 (0.21) | -2.283 | 2.5616 | 2.00252, 0.28043 | 6 | 1.0515 | 6 |
| 7P | 111 (0.30) | -1.7816 | 1.5525 | 2.0048, -0.22322 | 6 | 1.0717 | 7 |
| 7P | 116 (0.35) | -1.4016 | 0.78599 | 2.0066, -0.60501 | 6 | 1.0846 | 3 |
| 7P | 125 (0.44) | -0.42653 | -1.1873 | 2.01125, -1.58472 | 6 | 1.1118 | 3 |
| 7P | 131 (0.50) | 0.53416 | -3.1403 | -2.54998, 2.01582 | 3 | 2.0659 | 3 |
| 7A | 1 (0.0) | -2.9418 | 3.8837 | 2.0000, 0.94183 | on the line k = 2 (6 to numerical accuracy) | 1.0001 | 1 |
| 7A | 16 (0.035) | -3.0111 | 4.0224 | 1.9998, 1.01132 | 1 | 1 | 6 |
| 7A | 40 (0.155) | -3.197 | 4.3944 | 1.99939, 1.19759 | 1 | 1 | 6 |
| 7A | 75 (0.50) | -3.4772 | 4.9548 | 1.99935, 1.47786 | 1 | 1 | 1 |
| 7A | 100 (0.75) | -3.5958 | 5.1917 | 1.99959, 1.59616 | 1 | 1 | 1 |
| 7A | 120 (0.93) | -3.7489 | 5.4979 | 1.99973, 1.74918 | 1 | 1 | 1 |

Findings:
- 7A: region 1 (stable) for the whole family tested (e = 0.035 to 0.93), as printed (p51: below 0.75, "all stable, region 1"), with a1 and a2 moving towards the separation point (-4, 6) as printed (-3.75, 5.50 at e = 0.93). The orbit sits within 1e-3 of the line k = 2 throughout.
- 7P: region 6 for e = 0.0001 to 0.44 as printed (p39, p51: region 6, close to the line a2 = -2 a1 - 2). The largest planar modulus grows from 1.0008 to 1.11: mild, not the "88" obtained from the unrefined printed state at row 102. **At e = 0.50 (row 131) the orbit
  is region 3** (the second index k passes through -2 somewhere between e = 0.44 and 0.50, a flip, between rows 125 and 131); the printed report gave stability information for e up to 0.35 only. This is a computed extension, not a printed result, and was obtained from
  one re-corrected orbit at each e. The two orbits at either side were not continued to locate the crossing.
- The region from the PRINTED seven-digit state differs from the re-corrected one in nine of the fifteen rows above (counting 7A row 1, which is on the boundary), because k1 - 2 is of order 1e-4 to 1e-3 and the printed state is not quite periodic. This is the practical content of the caution in section 1 item 5.

## 10. Techniques applicable to the project's problems

### 10.1 `#912` (the corrector lacks what the published method needs)

- **Start epoch.** Add a starting true anomaly v0 in {0, pi} to the corrector. Table 13 (120 rows) is a published control for the apoapsis start in `core.er3bp` (section 8: half-revolution miss 1e-6 to 5e-6 at e = 0, 0.155, 0.5).
- **Conditions.** Unknowns (x0, y_dot0) with y0 = x_dot0 = 0 at v0; targets y = x_dot = 0 at v0 + k pi (Eqs. 146 to 149, solved by Cramer's rule Eq. 150); the Jacobian from the 2 x 2 block of the STM. Use event or terminal-condition integration or a fixed span k pi (fixed here, because the epoch is fixed by the apse condition).
- **Near-singular Jacobians:** the report switches to minimum-norm least squares (Lawson) when the matrix is nearly singular, which happens at the near-unit eigenvalue (k near 2, section 9). The corrector should do the same; the circular case at e = 0 is exactly singular.
- **Seeding recipe** (from the circular problem): take the catalogued symmetric circular orbits, interpolate to those whose period is 2 k pi (k = 1 to 5) in the synodic frame; each is a seed with TWO continuations (P and A) that differ for e > 0 and meet at e = 0. This is the same statement as Peng and Xu 2015, applied to the planar problem (they use N/M resonances in 3D). 150 resonant orbits came from 4000 circular ones.
- **Stability classification for an interval of e.** The monodromy is the one-period map (2 k pi), so the planar 4 x 4 block with a1, a2 and the seven regions of Table 3 are available; the segmented product that Peng and Xu use for 3D is a different object (6 x 6) and its planar block can be tested against this section.
- **Both signs.** Because the pulsating equations depend on f only through e cos f, an apoapsis start at (f0 = pi, e) is a periapsis start at (f0 = 0, -e) (already a test in `tests/core/test_er3bp_peng_xu_2015.py`). This report's Table 13 versus Table 12 shows the two families really differ for e > 0 (7P at e = 0.5: region 3; 7A: region 1).

### 10.2 `#925` (elliptic controls that do not reproduce: model or paper?)

- This report gives 251 elliptic rows (7P, 7A) and 151 (12A) that `core.er3bp` reproduces from the printed state (section 8). They show that the pulsating equations, the factor 1/(1 + e cos f), the mu convention and the apse start are right at the level of 1e-5 or better for e = 0 to 0.5 (7P, 7A) and e = 0 to 0.999 (12A, 9P). INFERRED: the two unresolved controls (Mako and Salamon unstable bands; Neelakantan and Ramanan M4N2) are then less likely to come from a defect in the planar equations as such, which leaves the paper, the start epoch, the stability classifier, or the three-dimensional terms. The three-dimensional z equation was shown algebraically identical to Table 1 (section 2.3) but not tested against a printed three-dimensional orbit (the report has none).
- The Mako and Salamon unstable bands are quoted in terms of a largest modulus; the report's Table 3 gives a finer classification: the planar quartet instability (region 2, Delta < 0) has modest growth and is invisible to a real-eigenvalue index. Compute a1, a2, Delta of the planar block at the printed unstable true anomalies before attributing the miss (the same recommendation as the Hadjidemetriou 1975b digest).

### 10.3 `#930` (the corrector can pass a bad orbit silently)

- Broucke's published guard, which the project can copy: the two relations among the characteristic-polynomial coefficients (the s^3 and s coefficients equal; the s^4 and s^0 coefficients 1), "precision of at least 1e-10", plus the series' redundant variables (Eq. 134) averaged each step. These test the INTEGRATION, not the closure. The closure test is the strong criterion itself (y = x_dot = 0 at the second apse; section 5.1); a state satisfying it to 1e-6 can still miss the full period by 1 (section 8, row 131), so a hard failure should be on the full-period miss of the CORRECTED state, not the printed one.
- Control (published numbers, with measured negatives): Table 12 row 60 and Table 13 rows 1, 40, 75 half-revolution states (section 8) accept (7P row 60: y_end -2.3e-7, x'_end -5.1e-8). Negatives on the same 7P row 60 state, measured: started at the wrong epoch f0 = pi (apoapsis) it ends with y = 0.43, x' = -1.06; integrated over a span of 1.5 pi instead of pi, y = -1.44, x' = -1.22; with a deliberately wrong mu (0.0121505) the x1 miss is 8.5e-4 (against 3.6e-7 at mu = 0.012155). The wrong-epoch case is the same idea as the existing test `test_apoapsis_start_does_not_close`.
- The tolerance of a closure test on a printed seven-digit state has to scale with the monodromy size (section 8, row 131: a change of 3.8e-6 in y_dot0 moves the full-period miss from 1.1e-10 to 2.2); do not use one absolute tolerance for all rows.

### 10.4 `#931` (stability classification and general-three-body controls)

- Use Table 3 as the sourced definition of the regime: with a1 = -trace, a2 = sum of principal minors, Delta = a1^2 - 4 a2 + 8, k = (-a1 +- sqrt(Delta))/2, the seven regions follow from (Delta sign; |k1|, |k2| against 2; the signs of the real k). Synthetic spectra for the tests: build a symplectic 4 x 4 matrix from each of the seven region patterns (this is Broucke's Fig. 3).
- Tolerance and refinement at the line k = 2 (section 9): the common case in this family of problems is |k - 2| of order 1e-4 to 1e-3; a classifier needs a tolerance, a statement of which side is "stable" at equality, and a re-corrected state.
- `floquet_classify` can never return "stable" for a symplectic monodromy (section 5.7); a region-1 orbit comes out "marginal". Add the planar 4 x 4 block and the (a1, a2, region) output; do not use `monodromy_eigenstructure` for regimes without a unit-circle pair (its docstring says it raises then; the `#931` entry records that it accepts any complex eigenvalue within 0.5 of the circle).
- First sourced regime controls from an independent paper: 7A (region 1, from e about 0.035 to 0.93, a1 and a2 approaching (-4, 6)) and 7P (region 6 up to e = 0.44) from this report; the printed text only claims these (no a1, a2 table is printed in part A), so the numbers in section 9 are computed here and must not be called printed.
- Relation of the two papers' conventions: a1 = alpha, a2 = beta, b = -k, Delta identical, same boundary lines. The Bray-Goudas sign erratum noted in the Hadjidemetriou digest (Delta < 0 is instability, not stability) is consistent with Table 3 here (region 2 is "complex instability").

### 10.5 `#899` (second-species generators, continuation through close approaches)

- The report's continuation in TWO parameters (e and mu) from a catalogued circular orbit, with the family parameter chosen as whichever of e or mu is convenient and the two corrected parameters (x0, y_dot0), is the same fixed-structure continuation the project's `genome` code does; the new element is the second parameter. Family 6A (e = 1, mu from 0.5 to 0.166) is a published continuation in mass ratio from an equal-mass orbit to a small mass ratio of the elliptic problem, ending at a bending towards a pseudo-parabolic path (p32), and shows that continuation in mu at fixed e is feasible when the base orbit has an apse symmetry.
- The 7P family ends at e about 0.50 because of a collision with the larger primary (p39): a published example of a continuation ending on a collision, the failure mode the `#899` entry expects, and a control (Table 12 rows 120 to 131) for a corrector approaching one. The 7A family continues to e = 0.93 tabulated (0.99 stated).
- Close approaches are handled in the report by Birkhoff regularisation (section 3.4); the regularised equations (Eq. 110) are compact and the full set of partial derivatives is printed (Eqs. 111 to 114), which makes an independent regularised propagator feasible (`#928`).

## 11. Anomalies and slips found while reading

1. Tables 6 and 7 number orbits 1 and 2 in opposite order (Table 7 NR 1 is Table 6 orbit 2). Table 6 truncates to six decimals (Table 7 has seven).
2. Table 13 ends at e = 0.93 (row 120); the text (p51) says e from 0.0 to 0.99; Fig. 14 shows e = 0.95 and 0.975 that are not tabulated.
3. Table 11 rows 9 and 10 are printed out of order (mass ratio 0.491 then 0.492; X0 0.6947297 then 0.6933511).
4. Eq. 36c,d (xi1 = +m2 r, xi2 = -m1 r) is the mirror image of the frame the tables of 3P, 4P, 5P use (smaller primary at +x). Tables 8, 10, 11 have X0 > 0.
5. Eq. 58 (p9): the second attraction term of the x equation reads (y - y1)/r2^3 in the image where (x - x2)/r2^3 is meant (reading of the image, not certain).
6. Eq. 141 prints "(9)" for "(q)" in the sums (scan or print artefact).
7. p51 text: "-2a2 - 2" should read a2 = -2 a1 - 2.
8. The text layer of the scan is poor (decimal points read as letters, dropped minus signs); the tables below were not taken from it.
9. 7P (p39): "orbits of this family belong to the class of satellite orbits ... but they come closer to the larger primary m1 than to m2": the orbits are around m1, so "satellite orbit of m2" in the introduction of the section refers to the circular-problem seed, not to the elliptic orbit.

## 12. How the tables were transcribed (provenance)

Page images were read at 150 dpi for every equation and table page. The 630 table rows were then transcribed by three independent optical-recognition passes at different scales, one cell at a time on the page-image rows, with a majority vote on the seven
decimals; the minus signs, the integer digits and the 33 disagreements were settled from the page image. Every row was then checked by integration against its own printed final state (section 8), and column smoothness was checked (the cubic fit through the four neighbours in e or
mu leaves relative residuals of at most 2e-5 for every column of Table 12, at most 1.2e-4 for Table 13 (largest in rows 111 to 116, e 0.86 to 0.92, where the family is steep), 5.4e-5 for Table 4 and 3e-7 for Table 5 y_dot0 (Table 5 X0 rows 3 to 6 deviate by up to 7e-3 near e = 1, where the curve is not smooth enough for a cubic; those rows close in the integration); no outlier). Errors that this process found and corrected from the images: dropped minus signs, a '9' for '0' in the eccentricity column,
a few single digits (for example Table 12 row 106 y_dot0 3.7684904, row 105 3.7292754; Table 13 row 93 y_dot0 3.7433082, row 86 x1 3.9738437, row 97 x1 5.3746733, row 117 x1 17.9287531), and Table 5 row 150 y_dot0 1.7947423. **A digit error in the sixth or seventh decimal that is
smaller than the integration check could see (about 1e-6 relative for the low multiplier rows) cannot be excluded** for the seven-digit numbers in Tables 4, 5, 12, 13; the integration reproduces them at the level stated in section 8.

## 13. Follow-ups (not registered)

1. A two-parameter apoapsis and periapsis corrector with a test pinned to Tables 12 and 13 (rows 1, 60, 102, 131 and rows 1, 40, 75), written from section 5.2 and 5.3 (`#912`).
2. A planar 4 x 4 stability classifier returning (a1, a2, Delta, k1, k2, region 1 to 7) with a documented tolerance at k = +-2, tested with a synthetic symplectic matrix per region and with the 7A rows (region 1) and 7P rows (region 6); retire the unreachable "stable" tag of `floquet_classify` or document that a symplectic orbit is at best "marginal" (`#931`).
3. Locate the 7P flip between e = 0.44 and e = 0.50 (k2 through -2) by continuing re-corrected orbits and report a1, a2 at the crossing.
4. A rectilinear (e = 1) inertial test integrator written from Eqs. 35 to 42 as a control for any future inertial elliptic code; tests from Tables 7, 8, 10, 11 (13 + 22 + 26 + 23 rows reproduce to 3e-5).
5. Check the `core.er3bp` three-dimensional (z) equation against a published three-dimensional elliptic orbit (this report has none).
6. Read Table 3 against the 8P, 8A, 10P, 11P, 11A families of Part B, which print stability regions as numbers.
7. Obtain Schubart's 1956 rectilinear orbit (Ref. in Section IV-A) and Bartlett (Ref. 18) to cross-check the Stromgren seeds of 12A and 8P.
