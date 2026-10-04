# Digest: Stiefel and Scheifele 1971, "Linear and Regular Celestial Mechanics"

E. L. Stiefel and G. Scheifele, *Linear and Regular Celestial Mechanics: Perturbed Two-body Motion, Numerical Methods, Canonical Theory*, Die Grundlehren der mathematischen Wissenschaften 174, Springer-Verlag, Berlin-Heidelberg-New York, 1971, xi + 301 printed pages, DOI 10.1007/978-3-642-65027-7. Filed in the private paper corpus as `stiefel-scheifele-1971-linear-regular-celestial-mechanics-grundlehren-174-springer-doi-10.1007-978-3-642-65027-7.pdf` (319 PDF pages, a scan with a usable text layer).
Page map (checked 2026-10-04): PDF page = printed page + 11 (printed p.2 is PDF 13, printed p.50 is PDF 61, printed p.140 is PDF 151, printed p.290 is PDF 301). Every page number below is the PRINTED page; a quoted equation "(9,67)" means equation 67 of section 9 (the book's own convention, p.9 footnote), and "Table n" is the book's own table number.
Which pages I read: the text layer of printed pp.1 to 173 (Parts I, chapters I to VII) and of chapter X (pp.229 to 257) in full; chapters VIII, IX and XI by section headings, theorem statements and the passages cited below only (their text layer, not images). Every equation I transcribe and every table digit quoted from sections 9, 15, 17, 23, 27 and 28 was read on a 130 dpi page image (PDF pages 32 to 35, 40 to 46, 72, 73, 81 to 86, 89 to 92, 94, 129 to 136, 150, 151, 164, 165, 170, 184, 246, 250 to 252); where the image and the text layer disagreed the image won (instances listed in section 12). Equations of chapter X other than (38,36), (38,59) to (38,80) were read on the text layer only and are marked "text layer".
Evidence tags: READ (printed page) is what the printed page says; COMPUTED is my own arithmetic or program of 2026-10-04 and 2026-10-05 (scratch programs, not committed; what each one did is listed in section 11); INFERRED is my reading across sources; RECALLED is from memory and not to be relied on.
Companions: `docs/notes/2026-10-04-digest-burdet-1967-regularization-two-body-problem.md`, `docs/notes/2026-10-04-digest-heggie-1974-global-regularisation.md`, `docs/notes/2026-10-04-digest-aarseth-zare-1974-regularization-three-body.md`, `docs/notes/2026-06-19-digest-szebehely-1967.md`, `docs/notes/2026-06-19-digest-bond-allman-2021-modern-astrodynamics.md`.

## 0. What the book is, and why it matters here

A textbook (Introduction pp.1 to 3) that re-expresses perturbed two-body motion as a perturbed harmonic oscillator in four dimensions (the Kustaanheimo-Stiefel, KS, transformation), with the physical time carried as a coordinate. Part I (chapters I to VII, pp.5 to 178) is elementary and numerical; Part II (chapters VIII to X, pp.181 to 267) is canonical theory; Part III (chapter XI) is geometry. The authors say they omit regularisations that do not give linear equations (p.1, citing Szebehely [1] for them) and that sections 26 and 27 (the G-functions) are Scheifele's (p.3).

Why it matters to the project: #928 (a float regularised propagator and transition matrix for close moon passes), #899, #924 and #929 are all perturbed two-body problems about a moon. This is the standard KS reference, the only held source that prints (i) a fully worked 3D regularised initial-value problem, (ii) hyperbolic and elliptic motion in one formalism, (iii) a perturbed-orbit test set with eight-digit reference values, and (iv) the canonical form of KS. It prints NO variational equations and NO transition matrix in KS variables (section 8 of this note); that is the central gap for #928.

Headline findings (details in the sections named):
1. The u-method (9,67), (9,68) holds for ANY sign of the energy; the elliptic-only pieces are the time element (18,39) to (18,42), the regular elements and generalised eccentric anomaly E (section 19), the E-form canonical systems (sections 38 and 41) and the E-form G-method (p.160). A hyperbolic moon pass is inside the u-method (section 9 of this note).
2. Levi-Civita is KS with u3 = u4 = 0 (collection of formulae item 6, p.35), so one KS code handles planar and 3D; COMPUTED check in the rotating CR3BP with a moon-centred origin (section 11).
3. Eight-digit reference values, in a concrete model (J2 plus a circular lunar orbit), are printed for 1 and 50 revolutions (Tables 1 to 2b); I reproduced the one-revolution values to 1e-4 km with an independent Cartesian integration and the 50-revolution oblateness row to 1e-4 km with a KS integration (section 11). That is the strongest published control in the held corpus for a KS propagator.
4. The book says (p.125, Example 6 comment) that for very close approaches to a perturbing moon "the centre of regularization should be shifted from the earth to the moon at an appropriate intermediate time"; the regularisation centre is one body only. That is the situation of #928.
5. The project code docstring `core/cr3bp_regularized.py` cites this book's "Ch. III" for the dt = r1 r2 ds form; chapter III is Kepler motion and the book contains no Lemaitre attribution (section 10).

## 1. Notation and the energy sign convention (read this before using any formula)

K^2 = k^2 (M + m) (p.7); for a massless spacecraft K^2 = k^2 M, the primary's gravitational parameter. Dots are d/dt, primes d/ds, (a, b) is the scalar product (p.9). r = |x|.
Energy: the book's h is the NEGATIVE total energy per unit mass. Kepler energy (-hK) with hK = K^2/r - v^2/2 (1,10) and h = hK - V when a perturbing potential V(t, x) acts (1,14). So h > 0 is elliptic, h = 0 parabolic, h < 0 hyperbolic (p.52, (12,1)); the semi-major axis is a = K^2/(2h) (10,18). The oscillator in the regularised equations is u'' + (h/2) u = ..., frequency w = sqrt(h/2) in s, and h = 2 w^2 (p.37, p.74).
Mapping to the project's physical energy (INFERRED, COMPUTED where stated):
- physical specific energy epsilon = -h. A hyperbolic flyby of a moon has epsilon > 0, so h < 0 and the u-equation is u'' - (|h|/2) u = ...; the oscillator becomes a hyperbolic (exponential) one. A sign slip here silently gives the wrong frequency.
- `scripts/certify_670_wz_oterma_regularized.py` uses H_ENERGY = -C_J/2 (the physical Jacobi energy with C the project Jacobi constant) as the fixed energy of the Levi-Civita Hamiltonian; in the book's convention that corresponds to h = +C_J/2 - c, where c is the additive constant of the perturbing potential (section 7 of this note). COMPUTED in a moon-centred setup (section 11): with V chosen so that V(0) = 0 and the project's `jacobi_constant` C = (x^2+y^2) + 2(1-mu)/r1 + 2mu/r2 - v^2, I found h (book) constant along the arc and equal to C/2 - [(1-mu)^2/2 + (1-mu)].
- the book's K^2 for a moon pass is the moon's gravitational parameter; in the project's nondimensional CR3BP with unit total mass it is mu.

## 2. Chapter I, Preliminaries (pp.6 to 17)

Equation of perturbed motion about the central mass (1,7): x'' + (K^2/r^3) x = P, P any force per unit mass, possibly depending on position, velocity and explicit time (p.8); "perturbing" does not mean small. Energy: d/dt(v^2/2 - K^2/r) = (x', P) (1,9), so the Kepler energy changes only by the work of P; with a potential, h' = -dV/dt (1,15) (p.10): a conservative (time-independent) V keeps the total energy constant.
Section 4 (pp.11 to 13): regularisation means regular DIFFERENTIAL equations, not regular solutions. The example (1 - cos s) x'' - sin s x' + x = 0 is singular at s = 0 but has the regular general solution c1 (1 - cos s) + c2 sin s; no integrator can pass the singular equation (p.12 to 13). "The purpose of regularization is not to obtain regular functions but regular differential equations."
Section 5 (pp.13 to 17), one-dimensional motion: first step is the fictitious time dt = x ds (5,25); for zero energy this turns x = (K^2 t^2 / ... ) branch-point motion into x = K^2 s^2/2, t = K^2 s^3/6 (5,28), "uniformisation". The second step is x = u^2 (5,31): the energy equation becomes u'^2 related and the equation of motion becomes the linear oscillator u'' + (hK/2) u = 0 (5,33); for hK > 0 the collision (ejection) orbit is the cycloid x = a(1 - cos 2ws), t-like variable = a(2ws - sin 2ws), with w^2 = hK/2 and a = K^2/(4 w^2) (5,30), Fig. 3. The radial case is the exact analytic control for #928 (section 12 of this note).
Section 5 also states (p.17) that not every solution of (5,33) represents a Kepler motion: the energy relation (5,32) is a constraint; the general solutions are "used in chapter VII for constructing approximations of perturbed orbits".

## 3. Chapter II, the regularised theory (pp.18 to 35)

### 3.1 Fictitious time and the matrices

Section 6: orthogonal matrix means A^T A = c E with scalar c (p.19). Section 7: dt = r ds (7,7), so the perturbed equation becomes r x'' - r' x' + K^2 x = r^3 (-dV/dx + P), t' = r (7,8), total order 7 (time is the fourth coordinate). It is still singular (a denominator r appears when solving for x''), p.20.

### 3.2 Levi-Civita in the plane (section 8, pp.20 to 23)

READ on image. x1 = u1^2 - u2^2, x2 = 2 u1 u2 (8,13), r = (u,u) = |u|^2 (8,12); angles at the origin are doubled. L(u) = [[u1, -u2], [u2, u1]] (8,15), x = L(u) u (8,17), x' = 2 L(u) u' (8,16), L^T L = (u,u) E (8,18), L^-1 = L^T / (u,u) (8,19), L(u)' = L(u') (8,21), L(u) v = L(v) u (8,22), and (u,u) L(v) v - 2(u,v) L(u) v + (v,v) L(u) u = 0 (8,23).
Result (8,26), p.23:
 u'' + [(K^2/2 - (u',u')) / (u,u)] u = ((u,u)/2) L^T(u) (-dV/dx + P).
By the energy relation (9,48) the bracket equals hK/2 (the printed text layer shows "(u,u')" in the bracket; the image shows (u',u')).

### 3.3 KS in space (section 9, pp.23 to 35)

KS matrix, READ on image p.24:
 L(u) = [[u1, -u2, -u3, u4], [u2, u1, -u4, -u3], [u3, u4, u1, u2], [u4, -u3, u2, -u1]] (9,27),
the Levi-Civita matrix in the upper-left corner. A physical 3-vector is padded with a fourth component 0; x = L(u) u (9,28) gives
 x1 = u1^2 - u2^2 - u3^2 + u4^2, x2 = 2 (u1 u2 - u3 u4), x3 = 2 (u1 u3 + u2 u4), x4 = 0 (9,29); r = (u,u) = |u|^2 (9,31); L^T L = (u,u) E (9,30); L(u)' = L(u') (9,33).
Bilinear relation (9,34): u4 v1 - u3 v2 + u2 v3 - u1 v4 = 0. Theorem 1: if u, v satisfy it then L(u) v = L(v) u. Theorem 2: under it, (u,u) L(v) v - 2(u,v) L(u) v + (v,v) L(u) u = 0 (9,35). The map x <- u is one-to-many (a one-parameter family, a circle of radius sqrt(r), chapter XI Theorem 1), so the equations of motion are postulated and then verified (Theorem 8), not derived. Initial data: u(0) arbitrary on the fibre over x(0), u'(0) from (9,37) (equivalently (9,71)); Theorem 4 shows the initial pair satisfies the bilinear relation; Theorem 6: l(u,u') = u4 u1' - u3 u2' + u2 u3' - u1 u4' is a first integral of (9,40) and, equal to zero initially, stays zero (Theorem 7); then x' = 2 L(u) u' (9,43).
Equations of motion (9,40): u'' + [(K^2 - 2|u'|^2)/(2|u|^2)] u = Q, Q = (|u|^2/2) L^T (-dV/dx + P) (9,38), and Q = (|u|^2/2)( -(1/2) dV/du + L^T P ) (9,46) because L^T(u) dV/dx = (1/2) dV/du (9,45) for a potential V(t, x(u)). Hence if there is only a potential, no L matrix appears in the equations. Kepler energy in KS: hK = (K^2 - 2|u'|^2)/|u|^2 (9,48) (the text layer's "K^2 - 2|u'|" is |u'|^2).
Regular set (9,52), READ on p.30 image (Kepler-energy form, total order ten):
 u'' + (hK/2) u = (|u|^2/2) ( -(1/2) dV/du + L^T P ), hK' = (dV/du, u') - 2 (u', L^T P), t' = (u,u).
Total-energy form (9,53), p.31 image, "normally a much better numerical precision":
 u'' + (h/2) u = -(1/4) d/du (|u|^2 V) + (|u|^2/2) L^T P, h' = -|u|^2 dV/dt - 2 (u', L^T P), t' = (u,u),
with h = hK - V (9,50). For a time-independent V and no P, h is constant and the oscillator frequency is constant (comment 2, p.31). Both sets "remain regular at collision (u = 0)" provided the perturbing forces stay finite (comment 1). A drag-like force P = lambda x' gives L^T P = 2 lambda |u|^2 u', which is regular (comment 3). The authors argue the rise in order (6 to 10) is harmless (comment 4, section 16).
Distance equation (9,54) to (9,59), p.32: r'' + 2 hK r = K^2 + r( -(1/2) dV/du + L^T P, u ), or with the total energy r'' + 2 h r = K^2 - 2 r V + r( -(1/2) dV/du + L^T P, u ) (9,55), equivalently in Cartesian terms r'' + 2 h r = K^2 - 2 r V + r( -dV/dx + P, x ) (9,58); unperturbed r'' + 2 hK r = K^2 (9,59), the generalisation of the one-dimensional (5,29).
Inverse KS, p.33 to 35 (READ on image): for x1 >= 0 choose u1^2 + u4^2 = (r + x1)/2 and set u2 = (x2 u1 + x3 u4)/(r + x1), u3 = (x3 u1 - x2 u4)/(r + x1) (9,69); for x1 < 0 choose u2^2 + u3^2 = (r - x1)/2, u1 = (x2 u2 + x3 u3)/(r - x1), u4 = (x3 u2 - x2 u3)/(r - x1) (9,70). Velocity: u1' = (1/2)(u1 x1' + u2 x2' + u3 x3'), u2' = (1/2)(-u2 x1' + u1 x2' + u4 x3'), u3' = (1/2)(-u3 x1' - u4 x2' + u1 x3'), u4' = (1/2)(u4 x1' - u3 x2' + u2 x3') (9,71) (dots on the right, i.e. physical velocity; this is u' = (1/2) L^T xdot), and h = K^2/r - |xdot|^2/2 - V at t = 0 (9,72). Checks: energy (9,73) h = (K^2 - 2|u'|^2)/|u|^2 - V, bilinear (9,74), and the r-equation (9,75). Item 6 of the collection (p.35): Levi-Civita is u3 = 0, u4 = 0, and for V = 0 (9,67) coincides with (9,52); the transformed force relation L^T dF/dx = (1/2) dF/du may be used for a potential part F.
Physical velocity from the state: xdot = 2 L(u) u' / |u|^2 (9,65); transformed force (L^T P)_i given in (9,66).

## 4. Chapter III, Kepler motion (pp.36 to 51)

Section 10, elliptic motion by Levi-Civita: with the frequency w = sqrt(h/2) (10,2) the plane solution is u1 = alpha cos ws, u2 = beta sin ws, the parametric orbit being an ellipse with the minor axis along u1 (Fig. 4); the eccentric anomaly is E = 2 w s (10,8); a = (alpha^2 + beta^2)/2 (10,9), a e = (beta^2 - alpha^2)/2 (10,10), alpha = sqrt(a(1 - e)), beta = sqrt(a(1 + e)) (10,11), x1 = a(cos E - e), x2 = a sqrt(1 - e^2) sin E (10,12), r = a(1 - e cos E) (10,13), tan(phi/2) = sqrt((1+e)/(1-e)) tan(E/2) (10,16), a = K^2/(2h) (10,18), dE/dt = K/(r sqrt(a)) (10,21), Kepler's equation t = (a^(3/2)/K)(E - e sin E) (10,26), period T = 2 pi a^(3/2)/K (10,27), M = E - e sin E (10,29) with Newton-Raphson correction (M - (E - e sin E))/(1 - e cos E) (10,30).
Section 11, ALL three conic types in one formalism (pp.42 to 50), READ. The single equation u'' + eta u = 0 with Taylor series gives u(s) = u(0) c0(eta s^2) + u'(0) s c1(eta s^2) (11,34) with the Stumpff functions c_n(z) = sum_k (-z)^k/(2k + n)! (11,35); the family s^n c_n(eta s^2) is closed under differentiation and integration (11,37), (11,38); c0(x^2) = cos x, c1(x^2) = sin x / x, c0(-x^2) = cosh x, c1(-x^2) = sinh x / x, c2(x^2) = (1 - cos x)/x^2, c3(x^2) = (x - sin x)/x^3, and so on (p.45); identities c0(z)^2 + z c1(z)^2 = 1, c0(z)^2 - z c1(z)^2 = c0(4z), c1(z)^2 = 2 c2(4z), c0 c1 = c1(4z) (11,45). Cancellation warning: the closed forms lose digits for small z; use the series (p.45).
Orbit with the pericentre on the x1-axis, start at pericentre (s = 0, t = 0): (READ on text layer, spot-checked on image for (11,51) to (11,59))
 hK = (K^2 - 2 u2'(0)^2)/u1(0)^2 (11,47), e = (1/K^2)[2 u2'(0)^2 - h u1(0)^2] (11,48), q = u1(0)^2 = pericentre distance (11,53),
 x1 = q - K^2 s^2 c2(2h s^2), x2 = K sqrt(q(1 + e)) s c1(2h s^2), r = q + K^2 e s^2 c2(2h s^2) (11,55),
 u1 = sqrt(q) c0(h s^2/2), u2 = (K/2) sqrt(1 + e) s c1(h s^2/2) (11,56) (argument (h/2) s^2, which is why the Levi-Civita orbit in u has half the frequency in s of x),
 r + e x1 = q (1 + e), so r = q(1 + e)/(1 + e cos phi), p = q(1 + e) (11,57): "Kepler's second law holds true for a pure Kepler motion of any kind";
 velocity x1' = -K^2 s c1(2h s^2), x2' = K sqrt(q(1 + e)) c0(2h s^2), r' = K^2 e s c1(2h s^2) (11,58), (x, xdot) = r rdot = r' = K^2 e s c1(2h s^2) (11,59), repeated in the collection as (11,70);
 generalised Kepler equation t = q s + K^2 e s^3 c3(2h s^2) (11,60), Newton-Raphson on it with derivative r (11,61). Parabolic (h = 0, e = 1): x1 = q - K^2 s^2/2, x2 = K sqrt(2q) s, t = q s + K^2 s^3/6 (Barker's equation), pp.48 to 49.
Collection of formulae, pp.50 to 51: items (11,64) to (11,73) with the elliptic specialisations in braces (a, E, M, T).

## 5. Chapter IV, the initial value problem (pp.52 to 72)

Section 12 (the ordinary Cartesian initial-value problem, any energy): h from (12,1); eccentricity and s from e c0(2h s^2) = 1 - 2h r/K^2, e s c1(2h s^2) = (x, xdot)/K^2 (12,2), e^2 = (1 - 2h r/K^2)^2 + 2h (x, xdot)^2/K^4 (12,3), pericentre distance q = r - K^2 e s^2 c2(2h s^2) (12,4), pericentre and orbital-plane unit vectors A = (c0/r) x - s c1 xdot, B = [1/sqrt(q(1+e))] [(s c1/r) x + (q - s^2 c2) xdot] for K = 1 (12,6), the printed K-dependent form being in the text layer only; circular orbits leave s undetermined (any point is a pericentre). Elliptic specialisation (12,7).
Section 13 (the KS use), pp.54 to 65:
- the dictionary L(u) u = x, L(u) u' = (r/2) xdot, L(u') u' = (1/2)(x, xdot) xdot - (v^2/4) x (13,15); (u,u) = r, (u,u') = (1/2)(x, xdot), (u',u') = (r/4) v^2 (13,16) (the text layer is garbled for the last; the numbers of section 15 confirm it: 1.111743 x 0.988674 / 4 = 0.274788, printed 0.274788).
- the parametric solution u(s) = u(0) c0((h/2) s^2) + u'(0) s c1((h/2) s^2) (13,20), pericentre data |u(0)| = sqrt(q), |u'(0)| = (K/2) sqrt(1 + e) (13,18; the second follows from (13,16) and (11,50), and I checked it that way);
- "The Problem of the Projectile" (pp.60 to 62), READ on images (PDF 71 to 73): the exact solution at fictitious time s from ANY initial x, xdot (no pericentre detour), valid for every sign of h:
 u(s) = c0((h/2) s^2) u + s c1((h/2) s^2) u' (13,38),
 x(s) = [1 - (K^2/r) s^2 c2(2h s^2)] x + [r s c1(2h s^2) + (x, xdot) s^2 c2(2h s^2)] xdot (13,39),
 r(s) = r c0(2h s^2) + K^2 s^2 c2(2h s^2) + s c1(2h s^2) (x, xdot) (13,40),
 t(s) = r s c1(2h s^2) + K^2 s^3 c3(2h s^2) + s^2 c2(2h s^2) (x, xdot) (13,41) (Stumpff's "principal equation"; derivative of the right side with respect to s is r, so Newton's rule converges), and the variant t = r s + r (v^2/2 - h) s^3 c3(2h s^2) + s^2 c2(2h s^2) (x, xdot) (13,42), where the middle coefficient (v^2/2 - h) is as printed on the image.
 Equation (13,39) is the Lagrange f and g solution in the universal variable; it is the exact analytic control for any pass, elliptic, parabolic or hyperbolic. Note the book's s equals the project's Sundman variable with dt = r ds and scale 1.
- ejection and "nearcentre" (pp.62 to 65): ejection (x = 0) has infinite velocity but the parametric data are finite: u(0) = 0 and u'(0) = (K / sqrt(2)) u0' with u0' the unit vector from (13,48) (13,49); for a start near the centre use the direction vectors x0, xdot0 and u = sqrt(r) u0, u' = sqrt((K^2 - r hK)/2) u0' (13,47) with (13,45), (13,46) for u0, u0'; "if a numerical procedure fails in an exceptional case, then the procedure suffers from loss of accuracy in situations approaching that singular case" (p.62).
Section 14, classical elements (pp.65 to 70): u in terms of node, inclination and argument (14,50) to (14,53); the three classical angles Omega, omega, J are undefined for circular orbits, equatorial orbits and collision (p.69); the book defines regular elements in section 19 instead. Collection (14,56), (14,57), and the angle recovery (14,58), (14,59).
Section 15, numerical examples (pp.70 to 72), K = 1. All values READ on images; the 3 examples are one data set:
 Example 1 (ordinary IVP, (13,30) to (13,33)). Data: x = (0.38623, 1.03156, 0.15060), xdot = (-0.90484, 0.33171, 0.24476), (x, xdot) = 0.029563, r = 1.111742, v^2 = 0.988674, 2/r = 1.798979, 2hK = 0.810305. e c0 = 0.099150, e s c1 = 0.029563, (e s c1)^2 = 0.00087397, e^2 = 0.0105389, e = 0.102659, c0 = 0.965819, s c1 = 0.287973. Iteration s = 0.291300 (iterates 0.291299, 0.291300; first approximation cos sqrt(2hK) s1 = c0, sqrt(2hK) s1 = 15 deg 01.4 min = 0.262206), 2hK s^2 = 0.0687590, c2 = 0.497142, s^2 c2 = 0.0421853, e s^2 c2 = 0.004331, q = 1.107411. A = (0.596104, 0.800638, 0.060349), B = (-0.781710, 0.561568, 0.271245) (c0/r = 0.868744, s c1/r = 0.259029, q - s^2 c2 = 1.065226, sqrt(q(1+e)) = 1.105032). Classical elements: J = 16 deg 08.0 min (three ways), Omega + omega = 53 deg 48.6 min, Omega - omega = 28 deg 43.6 min.
 Example 2 (inverse KS, (9,69), (9,71), u4 = 0): u = (0.865440, 0.595975, 0.087008, 0), u' = (-0.282049, 0.413169, 0.145277, 0.058505); bilinear value -0.12e-6; r = (u,u) = 1.111743, (u,u') = 0.014782, (u',u') = 0.274788.
 Example 3 (parametric IVP, (13,34) to (13,37), elliptic): hK = 0.405151 (the text layer omits it; image), 1/a = 0.810302, 1/sqrt(a) = 0.900168, e cos E = 0.099152, e sin E = 0.026613, e = 0.102661, cos E = 0.965816, sin E = 0.259231, cos(E/2) = 0.991417, sin(E/2) = 0.130738; alpha = (0.893193, 0.447428, 0.041871, -0.016149), beta = (-0.435593, 0.846970, 0.284074, 0.110473); bilinear value alpha4 beta1 - alpha3 beta2 + alpha2 beta3 - alpha1 beta4 = -0.14e-6; A = (0.596111, 0.800633, 0.060347), B = (-0.781705, 0.561573, 0.271245).
COMPUTED 2026-10-04: I recomputed examples 1 to 3 from the printed 5-digit inputs in double precision; every printed number is reproduced to about 2e-6 (u, u', r, (u,u'), (u',u'), e, s, q, A, B, E-quantities), the differences being the 6-digit hand arithmetic of 1971 (two printed routes give e = 0.102659 and 0.102661; mine 0.102660). Use a tolerance near 5e-6 on these. The printed bilinear values (-0.12e-6, -0.14e-6) are rounding of the inputs, not targets: with the input rounded to double precision the bilinear relation holds to 1e-17 in the 4-vector form.

## 6. Chapter V, the fundamental differential equations (pp.73 to 99)

### 6.1 Stability (section 16, pp.73 to 77), READ on images pp.74 and 75
For pure Kepler motion the regularised system is r'' + 2 h r = K^2 and u_j'' + (h/2) u_j = 0 (j = 1 to 4), t' = r (16,1), (16,2): linear, constant coefficients, total order 11 against order 6 for the nonlinear Newtonian equations. Theorem 1 (h > 0): every solution of the regularised system is stable in Lyapunov's sense. Deviation solutions about a reference solution (these are the only printed variational-type solutions in the book; h = 2 w^2): 
 Du_j(s) = Du_j(0) cos ws + [Du_j'(0)/w] sin ws, Du_j'(s) = -Du_j(0) w sin ws + Du_j'(0) cos ws,
 Dr(s) = Dr(0) cos 2ws + [Dr'(0)/(2w)] sin 2ws, Dt(s) = Dt(0) + [Dr(0)/(2w)] sin 2ws + [Dr'(0)/(4 w^2)] (1 - cos 2ws).
These hold at FIXED h; a variation of h is not covered (comment 5 below). Theorem 2: every elliptic solution of the Newtonian equations is unstable (two circular orbits with different radii and therefore periods drift apart; p.75). Comments (pp.75 to 77): the stabilisation comes from carrying the constant energy h into the equations; "introduction of elements into the differential equations may produce stabilization"; redundant equations are often well behaved (order 6 to 11 is no disadvantage); an erroneous value of h destroys the advantage, "the Kepler energy or the major axis must be given with high accuracy" (comment 5, repeated p.65).

### 6.2 Step regulation and truncation error (section 17, pp.77 to 83)
dt = r ds is an "analytical step regulation". Theorem 3 (p.78, READ): for dt = r^mu ds with mu >= 0, a necessary condition for reaching the central mass along a collision orbit is mu < 3/2 (from r = t^(2/3), s = integral of t^(-2 mu/3) dt finite). So dt = r^2 ds (the true anomaly, mu = 2) cannot reach a collision orbit's centre; dt = r ds (mu = 1) can. The authors "never used for our numerical experiments such sophisticated step-regulators". An automatic step-size control was sometimes added on top.
Error comparison, classical Runge-Kutta of order 4 on a circular orbit, K = 1, a = 1, initial data x = (1, 0), xdot = (0, 1) (11), exact x = (cos t, sin t) (12). Newtonian discretisation error after n steps of length h (13), (14), READ on image:
 Dx1 = -(1/2880) h^4 [66 t sin t + (45/2)(cos 2t + 2 cos t - 3)] + O(h^5), Dx2 = (1/2880) h^4 [66 t cos t - (45/2)(sin 2t + 2 sin t)] + O(h^5); leading growth linear in t (15).
Regularised problem (17), (18): u1'' + u1/4 = 0, u2'' + u2/4 = 0, t' = u1^2 + u2^2, u1(0) = 1, u2(0) = 0, u1'(0) = 0, u2'(0) = 1/2, t(0) = 0, exact u = (cos(s/2), sin(s/2)); error estimates (21), (22) (READ on the 260 dpi image):
 Du1 ~ (1/3840)(s/n)^4 s sin(s/2), Du2 ~ -(1/3840)(s/n)^4 s cos(s/2), Dt ~ -(1/768)(s/n)^4 s (the text layer has a spurious minus on Du1);
 explicit coordinate errors Dx1 ~ (1/1920)(t/n)^4 t sin t, Dx2 ~ -(1/1920)(t/n)^4 t cos t (23); the physical-time error Dt acts again on the coordinates (25), (26), giving total errors Dx1 ~ -(1/1280)(t/n)^4 t sin t, Dx2 ~ (1/1280)(t/n)^4 t cos t (27); coefficient ratio (66/2880)/(1/1280) = 88/3, "about 30 times better" than the Newtonian equations.
General elliptic regularised case with a = 1 (21a): Du1 ~ (sqrt(1-e)/3840)(s/n)^4 s sin(s/2), Du2 ~ -(sqrt(1+e)/3840)(s/n)^4 s cos(s/2), Dt ~ -(1/3840)(s/n)^4 [(5 - 2 e cos s) s + 5 e sin s]; conclusion "The accuracy of the numerical integration of the regularized equations is not sensitive to the value of the eccentricity" (the ratio to the circular estimate is at most sqrt(2) after one revolution).
Measured relative error (p.83 table, READ on image), classical RK4, two revolutions, 100 steps per revolution, with automatic step regulation in both columns; the Newtonian error is the relative error of the position after two revolutions:

| e | Newtonian | regularised | ratio |
|---|---|---|---|
| 0 | 0.8113e-5 | 0.01765e-5 | 46 |
| 0.1 | 1.060e-5 | 0.02159e-5 | 49 |
| 0.2 | 1.540e-5 | 0.02664e-5 | 58 |
| 0.3 | 2.239e-5 | 0.03336e-5 | 67 |
| 0.4 | 4.399e-5 | 0.04271e-5 | 103 |
| 0.5 | 9.545e-5 | 0.05652e-5 | 169 |
| 0.6 | 28.65e-5 | 0.07880e-5 | 364 |
| 0.7 | 97.44e-5 | 0.1205e-5 | 809 |
| 0.8 | 540.8e-5 | 0.2316e-5 | 2335 |
| 0.9 | 7697.0e-5 | 0.6332e-5 | 12160 |

The printed text says the e = 0 ratio "already exceeds the value predicted by the theory" (46 against about 30) and that, at e = 0.9, the Newtonian equations "are not recommended". COMPUTED 2026-10-05: plain fixed-step RK4 on the e = 0 case (100 steps per revolution, two revolutions, no automatic step regulation) gives an error of 7.7e-6 for the Newtonian equations and 1.70e-7 for the regularised equations (total error at the computed time), ratio 45.4, against 8.113e-6, 1.765e-7 and 46 printed: within 5 percent per column and 1.3 percent on the ratio; formula (14) predicts 4.5e-6 and (27) 1.53e-7 for the same run. Use only the ratio, not the absolute values, as a regression number for a fixed-step RK4 on the circular case.

### 6.3 The time element (section 18, pp.83 to 87), elliptic only (h > 0)
An "element" is any quantity that is a linear function of the independent variable in pure Kepler motion (constants included); with s the eccentric anomaly is an element (pp.83 to 84). The time element is tau = t + (1/h)(u, u') (18,32), i.e. t = tau - (1/h)(u, u') (18,38), with
 tau' = (1/(2h))(K^2 - 2 r V) - (r/(4h))(u, dV/du - 2 L^T P) - (h'/h^2)(u, u') (18,40) (text layer only),
h' as in (9,68), and the oscillator u'' + (h/2) u = -(1/4) d/du(|u|^2 V) - ... as (18,39). For pure Kepler motion tau' = K^2/(2h) = const, so any step method integrates it without discretisation error and the only remaining time error is the (u, u') term, negligible on a circle (comment 3, p.87). The set is stated "applicable only if h > 0" (p.87).

### 6.4 Regular elements and the generalised eccentric anomaly (section 19, pp.87 to 95), elliptic only
dE/ds = 2 w with w = sqrt(h/2) (19,45) defines the generalised eccentric anomaly E even for perturbed motion; then u = alpha(E) cos(E/2) + beta(E) sin(E/2) (19,50) with variation of constants (19,51), (19,52). Ten scalar elements: the frequency w (energy), the time element tau, and two 4-vectors alpha, beta (bilinear relation alpha4 beta1 - alpha3 beta2 + alpha2 beta3 - alpha1 beta4 = 0 (19,53)); they are "well defined for any pure elliptic Kepler motion even if collision occurs", hence "regular". Element equations (19,61) to (19,63), used in the form (23,49) to (23,52) in section 23. Initial values at E = 0: alpha = u, beta = 2 u*, with u* = du/dE (19,69). Collection, p.90 to 92: u* = -(1/2) alpha sin(E/2) + (1/2) beta cos(E/2) (19,55); x1 = u1^2 - u2^2 - u3^2 + u4^2 etc. (19,56); r = |u|^2 (19,57); physical time t = tau - (1/w)(u, u*) (19,59); checks (19,70), (19,71).
Comment 2 (p.92): if the perturbing force P is switched off at some E0 the elements stay constant and the particle moves on the osculating Kepler orbit. Comment 3 (pp.92 to 94): with a conservative potential V there are two options. (a) Treat V as the force P = -dV/dx, V = 0: w varies and the orbit is the osculating one. (b) Keep V as potential, P = 0: w is constant, the approximating orbit is the harmonic oscillation (19,71) with the modified frequency; "From the computational point of view the option (a) should be rejected" (Example 3 below: five digits lost). "In general it is a fruitful idea to work with a frequency w that is as constant as possible" (p.94).
Comment 6 (p.95): this element set "is restricted to perturbed elliptic motion", h = K^2/r - |xdot|^2/2 - V > 0 required; unrestricted sets are in section 40. The u-method (p.95): the original equations (9,67), (9,68) or (19,46) to (19,48); better than elements for at most one revolution with classical integrators, but "the u-method is however refined" in chapter VII until it "competes successfully with the element method". Near-parabolic initial conditions (oscillator term the same order as the perturbing term): the book recommends the u-method (u'' = 0 for pure parabolic motion; r is quadratic in s so the time quadrature is exact; p.95 to 96).
Appendix (pp.96 to 99): orthogonal elements (alpha-bar, beta-bar with the phase theta, (19,74) to (19,86)), re-orthogonalisation after each step; "did not produce an appreciable gain of precision" (p.99).

## 7. Chapter VI, typical perturbations (pp.100 to 126)

### 7.1 Potentials and the additive constant (sections 20 and 21, pp.100 to 110)
Third-body and body potentials, Legendre expansions (20,8) to (20,13), spheroids (20,18), (20,19) with J_n and R (p.105 to 106). Section 21 (p.107 to 109): the oblateness perturbing potential V = (1 + m/M) k^2 M Sum_{n>=2} J_n (R/r)^n P_n(cos theta)/r (21,23); the distance equation (21,26): r'' + 2 h r = K^2 + K^2 Sum_{n>=2} (n - 1) J_n (R/r)^n P_n(cos theta). The "problem of the additive constant": V is determined only up to a constant c; the regularised equations (9,67) contain V itself (through h and the potential term), so c changes the frequency and the numerics (not the physics). For the oblateness potential the choice is c = 0 (the potential must vanish at infinity so that the r-equation tends to the pure Kepler one; otherwise a term -2 r c grows without bound, p.109). The oblateness potential is infinite at the origin so the regularisation of the total potential "seems to be an unsolved or insoluble problem" (comment 2, p.109): the equations stay singular at the centre because of the J2 term.

### 7.2 Third body attraction (section 22, pp.110 to 117), READ on images pp.113 and 117
Particle m about central mass M at the origin; third body M' on a prescribed orbit a(t), distance rho = |a|, Delta = |x - a|. Perturbing force (22,29) P = -k^2 M' [ (x - a)/Delta^3 + a/rho^3 ] (principal part plus indirect part). Potentials: principal U = -k^2 M'/Delta (30), indirect U' = -(k^2 M'/rho^3)(x, a) (31), total V = -k^2 M' [1/Delta - (x, a)/rho^3] (32). If the third-body positions come from a table (no derivatives) use the force (29); if from formulae use the potential and differentiate it literally.
Interior problem (r < rho): V = -k^2 M'/rho - k^2 M' Sum_{n>=1} P_n(cos theta) r^n / rho^(n+1) (22,35). Rule for adjustment (box, p.113): add c(t) = +k^2 M'/rho to V so that V = -k^2 M' Sum_{n>=2} P_n(cos theta) r^n/rho^(n+1) (22,36): the new potential vanishes at r = 0 and the perturbation in the r-equation is O(r^3) instead of O(r). c depends only on time so the forces are unchanged, but the energy law (9,68) is modified. Exterior problem (r > rho): origin at the centre of mass of M1, M2, with V = -k^2 (M1 + M2 + m)/(M1 + M2) [ M1 (1/Delta1 - 1/r) + M2 (1/Delta2 - 1/r) ] (22,43), vanishing at infinity; no constant should be added (p.116). The adjustment rules "are applicable provided the distances r, rho of the particle and the third body satisfy either r << rho or r >> rho; in other cases no simple rule for adjusting the potential is available, and we recommend the use of the perturbing force (29) without introducing a potential" (comment 2, p.117, READ). Also, "any conservative potential - for instance an oblateness potential - should be used and not be replaced by the corresponding force", because with the potential plus a third-body force P* the energy changes only at the third-body rate h' = -(P*, xdot), whereas treating the oblateness as a force too gives h' = -(P, xdot) - (P*, xdot), a much larger rate (p.117). Hill's potential V = -k^2 M' P2(cos theta) r^2/rho^3 (22,47), the first term of (36).
Relevance (INFERRED): for a moon-centred pass with the planet as third body, r << rho holds on the pass and the interior rule (V(0) = 0) is the book's recommendation; I used it in the COMPUTED CR3BP check (section 11).

### 7.3 The numerical examples of section 23 (pp.118 to 125), with every printed number
Programme (p.118): the element equations (19,61) to (19,63) in the form (23,49) to (23,52) with the auxiliary 4-vector Q = -(1/4) d(rV)/du + (r/2) L^T P (23,48), integrated by fourth-order Runge-Kutta in E; CDC 6500, Fortran; 100 steps of size pi/50 in E is "roughly one revolution". Satellite (p.118 to 119): initial position x = (0, -5888.9727, -3400.0000) km, velocity xdot = (10.691338, 0, 0) km/s, mass m = 0 so K^2 = k^2 M = 3.98601e5 km^3 s^-2; pericentre r = 6800 km; unperturbed e = 0.95; inclination about 30 degrees; u4 = 0 chosen in (19,65). Oblateness: J2 term only, r V = lambda (x3^2/r^4 - (1/3)(1/r^2)), dV/dt = 0 (23,54), lambda = (3/2) K^2 J2 R^2, gradient of (rV) in (23,55), J2 = 1.08265e-3, R = 6.37122e3 km. Moon (interior problem, circular orbit, no third-body potential, force (22,29) used): ephemeris (23,56) READ on image p.119: a1 = rho sin(sigma t), a2 = -(sqrt(3)/2) rho cos(sigma t), a3 = -(1/2) rho cos(sigma t) (an orbit inclined 30 degrees to the x1x2 plane), rho = 384400 km, k^2 M' = 4902.66 km^3 s^-2, sigma = sqrt((K^2 + k^2 M')/rho^3) = 2.6653 1578 0887e-6 s^-1 (Kepler's third law). COMPUTED: sqrt((3.98601e5 + 4902.66)/384400^3) = 2.66531578089e-6, agreeing with the printed digits to 11 places. "Revolution time" is defined arbitrarily as the time at which E = 2 pi was reached in the first run. Time unit "m. solar day" = mean solar day (86400 s, INFERRED; my computed values below confirm it). Energy check column Dh/h (units 1e-8, defined Dh = h(E) - h(0), h(0) = 1.4724 04283 km^2 s^-2) and "Bil" (units 1e-13): u4 u1* - u3 u2* + u2 u3* - u1 u4* divided by |u||u*|.

Table 1 (oblateness only, 100 steps of pi/50 in E; image; Dh/h in units 1e-8, Bil in 1e-13):

| E (rad) | t (days) | x1 (km) | x2 (km) | x3 (km) | Dh/h | Bil |
|---|---|---|---|---|---|---|
| 0.00000000 | 0.00000000 | 0.0000 | -5888.9727 | -3400.0000 | 0 | -0 |
| 0.62831853 | 0.06395592 | 24888.0803 | 16498.7443 | 9544.2529 | -61 | 36 |
| 1.25663706 | 0.32255979 | 40231.8905 | 75105.8062 | 43406.3971 | -20 | 24 |
| 1.88495559 | 0.89613145 | 40172.2876 | 147544.8516 | 85249.0073 | -10 | 24 |
| 2.51327412 | 1.78468651 | 24732.0745 | 206146.6547 | 119089.7882 | -7 | 36 |
| pi = 3.14159265 | 2.86792773 | -191.1140 | 228527.3163 | 132002.7386 | -7 (a mark follows the digit; uncertain) | 74 |
| 3.76991118 | 3.95119440 | -25077.4679 | 206138.1834 | 119055.5554 | -7 | 36 |
| 4.39822972 | 4.83981608 | -40421.2469 | 147531.1447 | 85193.6174 | -10 | 24 |
| 5.02654825 | 5.41347008 | -40361.6506 | 75092.0991 | 43351.0072 | -20 | 24 |
| 5.65486678 | 5.67214058 | -24921.4055 | 16490.2721 | 9510.0204 | -61 | 36 |
| 2 pi = 6.28318531 | 5.73612194 | 3.5019 | -5888.9781 | -3399.9919 | -1 | 150 |
| 6.27904795 (t = 5.7359 3218 d, eq. (58)) | 5.73593218 | -171.7802 | -5888.0345 | -3399.6521 | -74 | 136 |

The last row is the position at the time (58) found by regula falsi; "the correct values corresponding to the time (58)" are x = (-171.7768, -5888.0346, -3399.6521) km; computing time for the table 1.2 s.
Table 1a (more revolutions, 93 steps per revolution with automatic step regulation; "revolution time" = (58) = 5.7359 3218 d; Dh/h in 1e-8, Bil in 1e-13): revolutions 0, 1, 5, 10, 50 at t = 0, 5.73593218, 28.67966090, 57.35932180, 286.79660900 d:

| rev | x1 | x2 | x3 | Dh/h | Bil |
|---|---|---|---|---|---|
| 0 | 0.0000 | -5888.9727 | -3400.0000 | 0 | 0 |
| 1 | -171.7761 | -5888.0347 | -3399.6521 | 2 | 0 |
| 5 | -857.7427 | -5865.5917 | -3391.3374 | 25 | 3 |
| 10 | -1708.4853 | -5796.3121 | -3365.7752 | 52 | 6 |
| 50 | -7745.1020 | -4043.7641 | -2775.1425 | -60 | 17 |

(computing time 1.4 s per revolution.) Table 1b (50 revolutions, automatic step regulation, steps per revolution; "more than about 200 steps do not improve the 8 digits"):

| steps/rev | x1 | x2 | x3 | Dh/h | Bil |
|---|---|---|---|---|---|
| 93 | -7745.1020 | -4043.7641 | -2775.1425 | -60 | 17 |
| 150 | -7745.1060 | -4043.7620 | -2775.1414 | -31 | -0 |
| 236 | -7745.1076 | -4043.7613 | -2775.1410 | 1 | -1 |
| 344 | -7745.1076 | -4043.7613 | -2775.1410 | 0 | -2 |

Example 2 (oblateness plus Moon, 100 steps of pi/50 in E; "revolution time" redefined as 5.7625 537882 d): Table 2 (h in km^2/s^2, Bil in 1e-13):

| E (rad) | t (days) | x1 | x2 | x3 | h | Bil |
|---|---|---|---|---|---|---|
| 0.00000000 | 0.00000000 | 0.0000 | -5888.9727 | -3400.0000 | 1.472404283 | -0 |
| 0.62831853 | 0.06395556 | 24888.0043 | 16498.6297 | 9544.1866 | 1.472401787 | 36 |
| 1.25663706 | 0.32258325 | 40231.6721 | 75111.3640 | 43409.6068 | 1.471938684 | 24 |
| 1.88495559 | 0.89647587 | 40162.1947 | 147602.7196 | 85282.4212 | 1.470773153 | 24 |
| 2.51327412 | 1.78633212 | 24661.1896 | 206365.6506 | 119216.2215 | 1.469486693 | 36 |
| pi | 2.87249227 | -436.8877 | 229017.4901 | 132285.6989 | 1.468724567 | 73 |
| 3.76991118 | 3.96033964 | -25620.9722 | 206893.9912 | 119491.8117 | 1.468564200 | 36 |
| 4.39822972 | 4.85458101 | -41259.4346 | 148374.4153 | 85680.2952 | 1.468622214 | 24 |
| 5.02654825 | 5.43379782 | -41293.6469 | 75750.6175 | 43730.9663 | 1.468668160 | 24 |
| 5.65486678 | 5.69676984 | -25642.9521 | 16743.2636 | 9655.8457 | 1.468769762 | 36 |
| 2 pi | 5.76300515 | -264.1754 | -6111.7102 | -3529.0015 | 1.468837473 | 129 |
| 6.27372868 (t = 5.76255379, eq. 59) | 5.76255379 | -672.6948 | -6099.8128 | -3522.5763 | 1.468835535 | 108 |

The exact values at the time (59) are x = (-672.6854, -6099.8133, -3522.5765) km; "the lunar perturbation is rather strong since the apogee is roughly at 2/3 of the distance of the moon"; 1.5 s. Table 2a (75 steps per revolution, automatic step regulation): rev 1: t = 5.76255379, x = (-672.6885, -6099.8139, -3522.5762), h = 1.468847614, Bil 82; rev 5: t = 28.81276894, x = (-45005.0833, 74170.7695, 42639.5792), h = 1.455353267, Bil 31; rev 10: t = 57.62553788, x = (-49984.4634, 134011.1715, 76967.0416), h = 1.450853872, Bil 33; rev 50: t = 288.12768941, x = (-24218.8545, 227961.9402, 129753.3046), h = 1.451086953, Bil 231 (rev 0 as Table 2 first row). Table 2b (50 revolutions, steps per revolution; "more than 200 steps do not improve"):

| steps/rev | x1 | x2 | x3 | h | Bil |
|---|---|---|---|---|---|
| 75 | -24218.8545 | 227961.9402 | 129753.3046 | 1.451086953 | 231 |
| 129 | -24219.0652 | 227962.1115 | 129753.4445 | 1.451086579 | 0 |
| 189 | -24219.0500 | 227962.1065 | 129753.4423 | 1.451086585 | -6 |
| 300 | -24219.0500 | 227962.1061 | 129753.4421 | 1.451086590 | -11 |
| 498 | -24219.0503 | 227962.1064 | 129753.4424 | 1.451086590 | -19 |

Example 3 (the oblateness orbit of Example 1 with option (a), V = 0 and the J2 term as a force P, 50 revolutions, 117 steps per revolution): t = 286.79660900 d, x = (-7721.8293, -4054.9312, -2780.2035) km; "loss of about five significant digits" against Table 1b; the reason for rejecting option (a).
Example 4 (small eccentricity e = 0.174, oblateness only, same initial point, xdot = (8.3000, 0, 0) km/s, 50 revolutions, t = 4.30899150 d); Table 4 (Dh/h 1e-8, Bil 1e-13), READ on the image (the text layer says 114 steps; the image says 84):

| steps/rev | x1 | x2 | x3 | Dh/h | Bil |
|---|---|---|---|---|---|
| 84 | 817.9032 (the last digit carries a smudge on the image; 2 is most likely) | -5960.7194 | -3175.7308 | 0 | -2 |
| 135 | 817.9034 | -5960.7192 | -3175.7306 | 0 | 6 |

Example 5 (the same orbit with oblateness as potential and Moon as force; element method against the u-method, both fixed step pi/90 in E, results at E = 2 pi), Table 5, x in km, t in mean solar days:

| method | t | x1 | x2 | x3 |
|---|---|---|---|---|
| u-integration | 0.08617 98373 96585 07 | 16.75756 31484 72 | -5889.00276 20052 98 | -3399.90938 73079 83 |
| element integration | 0.08617 98373 85137 70 | 16.75760 21409 89 | -5889.00276 23840 91 | -3399.90938 74068 08 |
| reference orbit | 0.08617 98373 85138 32 | 16.75760 21404 27 | -5889.00276 23840 67 | -3399.90938 74068 48 |

(the element method is better; "in a highly eccentric case the results of the u-method would be much closer to the results of the element method").
Example 6 (the book's only hyperbolic case): ejection from the centre of the Earth, parabolic initial Kepler energy hK = 0, moon force (22,29) only (oblateness not defined at the centre), u-method in the form (9,68), classical RK4. Direction of ejection x0dot = (0.9379 8034 755, 0.3376 9370 427, 0.07845 9095 724), inclined 4.5 degrees to the x1x2 plane; in this example the Moon's circular orbit lies IN the x1x2 plane with a1 = rho cos(sigma t), a2 = rho sin(sigma t), a3 = 0, rho = 384400 km, sigma = 2.6653 1578 0887e-6 s^-1 (READ on image p.124 to 125). Minimum distance from the moon about 48000 km ("this close approach required a drastic shortening of the step-size"); 82 steps to reach t = 5.0 mean solar days with seven significant digits; hK < 0 for all s > 0, so the osculating motion "is at any time of hyperbolic type". Table 6 (km):

| | x1 | x2 | x3 |
|---|---|---|---|
| perturbed | 659516.70 | 255497.87 | 39220.28 |
| unperturbed (rectilinear) | 651278.19 | 234474.58 | 54477.37 |
| difference | 8238.51 | 21023.29 | -15257.09 |

Comment (p.125): "The method thus outlined is not recommended in the case of very close approaches to the moon. If such an event occurs, the centre of regularization should be shifted from the earth to the moon at an appropriate intermediate time." Whether K^2 in Example 6 is the same 3.98601e5 as in the other examples is not restated (INFERRED yes).
Final remark (p.125): "Clearly our statements with respect to the numerical behaviour of the various underlying differential equations of celestial mechanics are not influenced by the numerical integration method which is adopted."
COMPUTED reproduction (2026-10-04 and 2026-10-05, section 11): K^2, J2, R, the Moon ephemeris (23,56) and sigma reproduce h(0) = 1.4724 04283 to 9 digits, the "correct values" at the times (58) and (59) to 1e-4 km, and the 50-revolution oblateness row of Table 1b to 1e-4 km (by a KS integration); the lunar 50-revolution rows of Tables 2a and 2b agree to a few 1e-3 km.

## 8. Chapter VII, refined numerical methods (pp.127 to 178)

Aim (p.127): numerical methods for x' = f(t, x) + eps g(t, x) that integrate the UNPERTURBED equation without discretisation error when eps = 0 (a "modified" method). The Newtonian equations with Runge-Kutta or difference methods, and Encke's method, do not have the property; the element method does (elements constant or linear), and the oscillator u'' + w^2 u = 0 can be given it by modified coefficients.
Section 24 (difference methods, due partly to Bettis): Cowell second-sum formula x(h) = 2 x(0) - x(-h) + h^2 [alpha0 f1 + alpha1 Df + ...]; the classical coefficients from (1/z... ) log expansion (24,4), (24,5); the modified highest two coefficients alpha_{n-1}, alpha_n are chosen so that f = a cos wt + b sin wt is integrated exactly, via the rational sequences S_m(u), R_m(u) of (24,9) with u = 4 sin^2(sigma), sigma = wh/2 (24,12); modified Stormer, Cowell, Adams-Bashforth, Adams-Moulton formulas (24,24) to (24,28) in the collection pp.132 to 134 (coefficients not transcribed here). Section 25: applying them to u'' + w^2 u = eps g; for the u-method the modified coefficients are computed once for the initial w and, for weak non-conservative forces, need not be re-adjusted every step; the time integral t = Int r ds uses Adams-Bashforth and Adams-Moulton with more retained differences than the u integration; the recommended set is (19,46) to (19,48) with constant frequency 1/2, where the modified coefficients are exact (comment 3, p.137). Long-term behaviour (pp.138 to 139): second-difference Cowell with the classical alpha2 = 1/12 integrates e^(i w t) at a drifting angle, sin sigma' = sigma [1 + 4 alpha2 sigma^2]^(-1/2) (25,38); the modified coefficient alpha2 = (1/4)[1/sin^2(sigma) - 1/sigma^2] (READ on image p.139) restores the exact points; fourth-difference classical coefficients spiral inward.
Example 7 (pp.139 to 140): a satellite on an almost circular orbit with a period of about 24 hours, 30 steps per revolution, about 50 revolutions, modified Adams-Bashforth and Adams-Moulton with the first six differences. Printed initial conditions (image p.139): x = (0, -37159.6555, -20000.0) km, xdot labelled (0, 0, 3,0733) km/s. COMPUTED: with the velocity directed along x1 (not x3) the orbit is circular at r = 42200 km as the tables show; with the speed 3.07335669 km/s = sqrt(K^2/r) (the printed 3,0733 is truncated; a speed rounded to 3.0733 moves the satellite about 700 km along the track by the end because of the mean-motion sensitivity) my Cartesian integration reproduces the printed reference rows of Tables 7a and 7b to 1e-4 km in all three coordinates (section 11). Treat "xdot3 = 3,0733" as a typesetting slip. Tables 7a (J2 only) and 7b (J2 plus Moon of (23,56)), READ on image:

| | x1 (km) | x2 (km) | x3 (km) | r (km) | t (days) |
|---|---|---|---|---|---|
| 7a modified | 216.9382 | -37160.1011 | -19997.9955 | 42200.0000 | 49.9260052 |
| 7a classical | 216.9335 | -37160.1188 | -19998.0051 | 42200.0202 | 49.9260171 |
| 7a reference | 216.9382 | -37160.1011 | -19997.9955 | 42200.0000 | 49.9260052 |
| 7b modified | 301.4760 | -37160.3635 | -19996.8650 | 42200.2147 | 49.9249537 |
| 7b classical | 298.3180 | -37160.4118 | -19996.9080 | 42200.2551 | 49.9249774 |
| 7b reference | 301.4918 | -37160.3634 | -19996.8649 | 42200.2146 | 49.9249535 |

Comment (p.141): as w tends to zero the modified and classical coefficients coincide (the parabolic case, where the classical coefficients are exact), so the modification pays more as the semi-major axis shrinks; "Perhaps the most important advantage ... is the numerical stability" over many revolutions.
Section 26 (pp.141 to 150, the G-functions, Scheifele): for x'' + alpha x = eps f(x, x', t) define G_n(t) = t^n c_n(alpha t^2) (26,59), the Stumpff functions; G_n' = G_(n-1), G_n+2 + ... recurrence (26,48); the solution is built from G0, G1 for the homogeneous part and G_(n+2) for polynomial forcing; the Taylor series is re-expressed as a G-series by b0 = a0, b1 = a1, b_(k+2) = a_(k+2) + alpha a_k (26,62), which is "no extra labour" because the b_k are obtained along with the a_k by recurrences (95), (96). Theorem (pp.147 to 148): with R_n(t) and r_n(t) the residuals of the G-method and of the power method, R_n(t) ~ eps c t^(n-1)/(n-1)! while r_n(t) ~ (-alpha a_(n-1) + eps c) t^(n-1)/(n-1)!: the G residual carries the factor eps, so "even if eps = 0 the method of power series expansion produces a truncation error, whereas already the two first terms of the G-expansion integrate the differential equation exactly". Step control by the number of terms (not by shortening the step); never use plain power series for weak perturbations (remarks, pp.149 to 150).
Section 27 (pp.150 to 160): the u-method with oblateness: u_j'' + alpha u_j = ..., alpha = h/2 (27,105); r'' + 4 alpha r = K^2 - 2 r V - (1/2)(dV/du, u) (27,106); t' = r (107); the modified distance rbar = r - K^2/(4 alpha) (110) and Steffensen auxiliary variables Q = 1/r^3, F = 1/r^2 with recursions (105a) to (116a); the time increment Dt = (K^2/(4 alpha)) h + Sum sigma_k/(k+1) (k+1)! G_(k+1)(4 alpha; h) (p.153); for non-conservative forces (third body) the independent variable becomes E of (19,45) with (118), (119) (frequency 1/4 constant), the Moon ephemeris interpolated by a polynomial and re-expressed in E by auxiliary series T_m = t^m (pp.155 to 159). The G-method is "elliptic only" (E independent variable, p.160 comment).

Tables (READ on images), the oblateness satellite of Example 1 at t = 5.73593218 days after 100 steps, with n terms:

Table 8 (G-method): 

| n | x1 | x2 | x3 | r |
|---|---|---|---|---|
| 2 | -171.6919 | -5888.2081 | -3401.2795 | 6801.8018 |
| 3 | -171.8571 | -5888.0339 | -3399.6518 | 6801.1848 |
| 4 | -171.7765 | -5888.0347 | -3399.6520 | 6801.1831 |
| reference | -171.7768 | -5888.0346 | -3399.6521 | (not printed) |

Table 8a (ordinary power series, same data):

| n | x1 | x2 | x3 | r |
|---|---|---|---|---|
| 2 | -22379.8218 | 10486.2036 | 6038.9283 | 12811.1035 |
| 3 | 1499.0814 | -5815.4568 | -3355.6255 | 6878.2359 |
| 4 | -170.6858 | -5888.0230 | -3399.6476 | 6805.3415 |
| reference | -171.7768 | -5888.0346 | -3399.6521 | (not printed) |

Example 9 (50 revolutions of the same satellite, 100 steps per revolution, automatic n; about n = 10 near pericentre and 3 near apocentre), Table 9: t = 286.79660900 d, x = (-7745.1086, -4043.7609, -2775.1406), r = 9167.3396 km. Example 9a (e = 0.174 satellite, 50 revolutions, 50 steps per revolution; the element method of Example 4 needed about 80): Table 9a: t = 4.30899150 d, x = (817.9031, -5960.7194, -3175.7307), r = 6803.2643 km. Note Table 9 differs from the converged Table 1b row by about 1e-3 km; both are the book's own values.
Example 10 (the Example 5 data, one revolution, E = 2 pi, 100 steps; n G-functions), Table 10 (t in days):

| n | t | x1 | x2 | x3 |
|---|---|---|---|---|
| 3 | 0.086179839538 | 16.76028152 | -5889.00278921 | -3399.90944937 |
| 4 | 0.086179837386 | 16.75760225 | -5889.00275412 | -3399.90934248 |
| 5 | 0.086179837385 | 16.75760214 | -5889.00276238 | -3399.90938740 |
| reference | 0.086179837385 | 16.75760214 | -5889.00276238 | -3399.90938741 |

Table 10a (ordinary power series, n powers):

| n | t | x1 | x2 | x3 |
|---|---|---|---|---|
| 3 | 0.086174883677 | 20.95746487 | -5889.04384923 | -3399.91887732 |
| 4 | 0.086179834417 | 16.75822481 | -5888.97883136 | -3399.89557523 |
| 5 | 0.086179838362 | 16.75739457 | -5889.00275991 | -3399.90938682 |
| reference | 0.086179837385 | 16.75760214 | -5889.00276238 | -3399.90938741 |

"The gain of precision by the use of the G-functions is about four digits, but the computing time ... is the same" (p.160). Final remark (p.160): which of the element method or a refined u-method is better depends on the problem and on the machine (cost of sines and cosines).
Section 28 (first-order methods, pp.160 to 178): insert a fixed reference orbit u = alpha(0) cos(E/2) + beta(0) sin(E/2) into the right-hand sides of the element equations; the equations then become independent quadratures; Fourier-expand with respect to E (and for a time-dependent third body with respect to Hansen's auxiliary variable E' of (28,141), double Fourier series); secular terms a00 E and resonances m lambda1 mu' + n = 0 (28,149). Theorem (p.171): to first order the frequency w has no secular term under a third-body perturbation and no quadratic secular term appears in the time element. Convergence of the Fourier expansions of 1/Delta (28,156) to (28,163): rates k (in E) and kappa (in M), kappa >= k; Table 11 (READ on image p.173), rates of convergence of the E-expansion k and the M-expansion kappa against the perturbing-body distance rho:

| e (rho0) | | rho = 1.5 | 2.0 | 10.0 | 100.0 |
|---|---|---|---|---|---|
| 0.0 (inf) | k | 0.667 | 0.500 | 0.100 | 0.010 |
| | kappa | 0.667 | 0.500 | 0.100 | 0.010 |
| 0.5 (3) | k | 0.474 | 0.377 | 0.089 | 0.009 |
| | kappa | 0.713 | 0.666 | 0.637 | 0.637 |
| 0.9 (0.211) | k | 0.310 | 0.254 | 0.066 | 0.007 |
| | kappa | 0.969 | 0.969 | 0.969 | 0.969 |
| 1.0 (0) | k | 0.209 | 0.172 | 0.046 | 0.005 |
| | kappa | 1.000 | 1.000 | 1.000 | 1.000 |

(the critical distance separating the two cases is printed as rho0 = 2 (1/e - e) (163), READ on the image; the table lists rho0 = 3 at e = 0.5 and 0 at e = 1, which agree with the formula, but 0.211 at e = 0.9 where the formula gives 0.422. I did not resolve which carries the slip; the k and kappa entries are unaffected at the tabulated distances). "For larger values of the eccentricity the M-series converges very slowly and it is in fact useless for e > 0.5 whereas the E-series converges well even if collision occurs" (p.173).
