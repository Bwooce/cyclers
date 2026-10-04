# Digest: Rodriguez del Rio 2021 PhD thesis, part A (front matter, chapters 1-3), with Olle, Rodriguez & Soler 2018 compared

Date: 2026-10-05 (Sydney). Reading, two checks run in scratch (a Levi-Civita integrator written from the printed equations, and the
core.cr3bp Jacobi constant and equations of motion); no project code was changed and no test files were written.

Source: O. Rodriguez del Rio, "Ejection-collision orbits in the Restricted Three-Body Problem", PhD dissertation, Universitat
Politecnica de Catalunya, Programa de Doctorat en Matematica Aplicada, Barcelona, June 2021; advisors M. Olle Torner and J. Soler
Villanueva; DOI 10.5821/dissertation-2117-351117; 180 PDF pages. Filed in the private paper corpus as
`rodriguez-del-rio-2021-ejection-collision-orbits-restricted-three-body-problem-phd-thesis-upc-doi-10.5821-dissertation-2117-351117.pdf`.
Page numbering: PDF page = printed page + 14 (printed page 15, the start of chapter 2, is PDF page 29). Pages below are PRINTED pages.
Part B (chapters 4-7: the analytical existence proofs, the Lyapunov-orbit and transit-region material, the 3D case) is in
`docs/notes/2026-10-05-digest-rodriguez-del-rio-2021-thesis-part-b.md`.

Evidence tags: READ (page) means read on the page image; TEXT means read from the PDF text layer only (the layer is reliable for
this thesis except for equation layout, so every equation I rely on was either read on an image or re-derived, see section 6).
Chapters 1-3 (this part) have NO numbered tables (chapters 4-7, part B, do: e.g. its Table 1 of 2-EC angles). Every printed number is in running text or a figure caption; there are no data tables to transcribe.

Companion digests: `docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`,
`docs/notes/2026-10-05-digest-alvarez-ramirez-barrabes-medina-olle-2019-ejection-collision-four-body.md`,
`docs/notes/2026-10-05-digest-alvarez-ramirez-barrabes-medina-olle-2021-ejection-collision-two-dof.md`. No digest of Llibre &
Martinez Alfaro 1985 existed when this was written (the paper is cited by the thesis as [LMA85]).

## 0. What part A is, and which paper each chapter is

| Thesis chapter | Content | Published as |
|---|---|---|
| Summary, Introduction (pp. i-4) | aims, the 7-chapter plan, the integrators used | the Introduction lists four papers (below) |
| Ch. 1 The RTBP (pp. 5-14) | planar circular RTBP, Hill problem, 3D RTBP, conventions | background only |
| Ch. 2 Local regularisation (pp. 15-30) | McGehee blow-up, Levi-Civita, three global regularisations | background; McGehee part is also section 2-3 of ORS 2018 |
| Ch. 3 Numerical n-EC orbits (pp. 31-65) | definition of n-EC, three methods, 1-EC and n-EC families, bifurcations, evolution | ORS 2018 (CNSNS 55:298), extended |
| Ch. 4 (pp. 67-) analytic existence I | four n-EC orbits for small mu | ORS 2020 (CNSNS 90:105294) per the Introduction list; part B |
| Ch. 5 analytic existence II, Hill problem | threshold C-hat(mu, n) for all mu | ORS 2020 / the 2022 paper; part B |
| Ch. 6 transit regions, LPO1, colour diagrams | | ORS 2021 (CNSNS 94:105550); part B |
| Ch. 7 3D case | 3D McGehee | ORS 2022 (CNSNS 111:106410); part B |

The Introduction (TEXT p. 3-4) says the thesis is made up of [ORS18] (CNSNS 55:298), then the 2020, 2021 and 2022 CNSNS papers. The
assignment of chapters 4-7 to 2020/2021/2022 is the other agent's to confirm; for part A I confirmed against the 2018 manuscript
itself that **chapter 3 is the published form of ORS 2018**: sections 2-4 of the paper are the McGehee setting (thesis section 2.1),
the collision manifold (thesis 2.1.1-2.1.2), and the numerical n-EC study (thesis chapter 3). The thesis chapter 3 is a strict
superset: it adds Levi-Civita, the angular-momentum method, the parameterisation method, the range mu in (0, 1) and n up to 100.
See section 8 for the line-by-line comparison.

Introduction facts (TEXT p. 4): all computations in double precision; integrators are the authors' own Runge-Kutta (7)8 (Fehlberg 1968)
and Runge-Kutta (8)9 (Verner 1978) with adaptive step control of Dormand and Prince 2005, plus a Taylor method (Jorba and Zou 2005).

## 1. RTBP conventions (chapter 1, pp. 5-14) compared with `core/cr3bp.py`

### 1.1 What the thesis prints (READ pp. 5-9 for the structure; eq. 1.4 and 1.9 checked on pages)

- Units: m1 + m2 = 1, primary separation 1, primary period 2 pi (so G = 1, mean motion 1). m2 = mu, m1 = 1 - mu, mu in (0, 1).
- Synodic frame: origin at the barycentre; **the primary of mass 1 - mu (P1) is at (mu, 0); the primary of mass mu (P2) is at
  (mu - 1, 0)** (p. 12, eq. 1.21 text; implied by r1 = sqrt((x - mu)^2 + y^2), r2 = sqrt((x - mu + 1)^2 + y^2)).
- Equations (1.3): xddot - 2 ydot = Omega_x, yddot + 2 xdot = Omega_y, with
  Omega = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2 + mu (1 - mu)/2 = ((1 - mu) r1^2 + mu r2^2)/2 + (1 - mu)/r1 + mu/r2 (eq. 1.4).
  The identity of the two forms holds: (1 - mu) r1^2 + mu r2^2 = x^2 + y^2 + mu (1 - mu) (I expanded it).
- Jacobi constant (eq. 1.9): **C = 2 Omega - xdot^2 - ydot^2**, with H = -C/2 (Hamiltonian momenta p_x = xdot - y, p_y = ydot + x).
  Hill region R(C) = {2 Omega >= C} (eq. 1.10).
- Symmetry (1.6): (t, x, y, xdot, ydot) -> (-t, x, -y, -xdot, ydot).
- Lagrange points: x_L2 < mu - 1 < x_L1 < mu < x_L3 (this ordering is stated for the frame above; the naming L1/L2/L3 is relative to
  P2 at mu - 1, see 1.3). Series (eq. 1.7): x_L1(mu) = 1 - (3^(2/3)/3)(1 - mu)^(1/3) + (3^(1/3)/9)(1 - mu)^(2/3) - (26/27)(1 - mu) + O((1 - mu)^(4/3))
  near mu = 1. C_L1(mu) = 3 + 9 ((1 - mu)/3)^(2/3) - 7 (1 - mu)/3 + O((1 - mu)^(4/3)) (eq. 1.11, near mu = 1).
- Routh value (p. 8): the thesis prints mu_R = (1/2)[1 - sqrt(69)/9] (the text layer garbles the radical) and the number **mu_R = 0.03852089650455137181950** (TEXT). I evaluated the closed form: 0.0385209, agrees to the digits I computed by hand; the 23 printed digits were not independently reproduced.
- Lyapunov orbits (p. 10): the monodromy eigenvalues are 1, 1, lambda, 1/lambda; the thesis states the largest real eigenvalue of
  LPO1 is "approximately between 2000 and 4000" for all mu in (0, 1) over C in [C_L2,3, C_L1) (Figure 1.4, a colour map).
- Hill problem (p. 11-12): x = mu + (1 - mu)^(1/3) x_h, y = (1 - mu)^(1/3) y_h; Hill potential Psi = (3/2) x_h^2 + 1/sqrt(x_h^2 + y_h^2)
  (eq. 1.14); K = 2 Psi - xdot_h^2 - ydot_h^2 (1.18); **C = 3 + (1 - mu)^(2/3) K** (1.19); K_L = 3^(4/3) at both equilibria
  at (+-3^(-1/3), 0). So the limit is mu -> 1: P1 is then the SMALL body.
- 3D (pp. 12-14): same frame plus z; Hamiltonian (1.23); symmetries (1.25) about the (x,y) and (x,z) planes; C = 2 Omega - v^2 (1.26).

### 1.2 Term-by-term mapping to `core/cr3bp.py` (frame, units, constants)

`core/cr3bp.py` (lines 3-4, 199-210, 252-259): mu = m2/(m1 + m2); primary of mass 1 - mu at (-mu, 0, 0); secondary of mass mu at
(1 - mu, 0, 0); r1 to the first, r2 to the second; eom ax = x + 2 vy - (1-mu)(x+mu)/r1^3 - mu (x-1+mu)/r2^3;
jacobi_constant = (x^2 + y^2) + 2(1 - mu)/r1 + 2 mu/r2 - v^2, which **omits** the constant mu (1 - mu).

| Item | Thesis (mu_T) | core.cr3bp (mu_c) |
|---|---|---|
| Equation form | xddot - 2 ydot = Omega_x etc. | identical form (same Coriolis sign, same rotation sense) |
| Frame orientation | primary of mass 1 - mu_T at +mu_T, primary of mass mu_T at mu_T - 1 | primary of mass 1 - mu_c at -mu_c, primary of mass mu_c at 1 - mu_c |
| Exact identification | **mu_c = 1 - mu_T, with NO reflection and NO rotation**: then mass 1 - mu_T sits at x = 1 - mu_c = mu_T and mass mu_T sits at x = -mu_c = mu_T - 1. | |
| Which body the thesis ejects from | P1 = the primary of mass 1 - mu_T, which in core terms is the secondary (mass mu_c) at 1 - mu_c | |
| Jacobi constant | C_T = 2 Omega - v^2 including the constant mu_T (1 - mu_T) | omits it: **C_T = C_c + mu (1 - mu)** (symmetric in mu <-> 1 - mu) |
| Hamiltonian | H = -C_T/2 | |
| Units, period 2 pi, primary separation 1 | same | same |

Verified numerically (section 6): for a state integrated in the printed Levi-Civita system at mu_T = 0.1, C_T = 5, the state mapped
into core coordinates has `jacobi_constant(state, mu_c = 0.9)` equal to C_T - mu_T (1 - mu_T) = 4.91 to 1e-14, and the core equations
of motion integrated from that state reproduce the Levi-Civita trajectory to 1e-14 over 0.3 time units.

Second equivalent reading (the 2018 paper's convention, section 8): for mu in (0, 0.5] with P1 the BIG primary, the 2018 paper's
frame (mass 1 - mu at (-mu, 0), mass mu at (1 - mu, 0)) is exactly core's frame with mu_c = mu. The thesis frame is that frame
rotated by pi about the barycentre (x, y) -> (-x, -y), which preserves the dynamics. So for the project's usual "mu_c <= 0.5, eject from
the BIG body" the thesis angles map as theta_thesis = theta_core + pi (McGehee angle), and the Levi-Civita angle theta_0 as pi/2 shifts; for
"eject from the SMALL body" (the Moon) use mu_T = 1 - mu_c = 0.9878494157 for Earth-Moon, no rotation.

### 1.3 Earth-Moon numbers (computed here, not printed in the thesis)

For the Moon as P1: mu_T = 1 - 0.0121505843 = 0.9878494157. C_L1 in core convention 3.188341105672621 (my root of the collinear
equation between the primaries, x_L1 = 0.836915132216393); in the thesis convention C_T,L1 = 3.2003440532737897. The thesis series
(1.11) at this mu_T gives 3.200327722890686: difference 1.6e-5, consistent with an O((1 - mu)^(4/3)) remainder (checked at three
other mu, section 6). Any claim in the thesis of the form "C >= C_L1" must therefore be tested against 3.2003 in thesis convention,
3.1883 in core convention.

### 1.4 Printed defect found on a page image

Eq. (1.5) and eq. (1.23) print the constant term of the Hamiltonian as **+ (1/2) mu (1 - mu)**, while eq. (2.1) (READ p. 15) and the 2018
paper print **- (1/2) mu (1 - mu)**. The minus sign is correct: with p_x = xdot - y, p_y = ydot + x one gets
(p_x^2 + p_y^2)/2 + y p_x - x p_y = (xdot^2 + ydot^2)/2 - (x^2 + y^2)/2, and H = -C/2 with Omega containing +mu(1 - mu)/2 requires
H = ... - mu(1 - mu)/2. So (1.5) and (1.23) have a sign slip; only the constant is affected, no dynamics.

## 2. Local regularisation (chapter 2, pp. 15-30)

### 2.1 McGehee blow-up of the first primary (pp. 15-19; eq. 2.1 READ p. 15; the rest TEXT, cross-checked against the 2018 paper)

Translate P1 to the origin: xbar = x - mu, ybar = y, p_xbar = p_x, p_ybar = p_y. Hamiltonian (2.1) (READ p. 15):

    H = (p_xbar^2 + p_ybar^2)/2 + ybar p_xbar - xbar p_ybar - mu p_ybar - (1 - mu)/r1 - mu/r2 - mu (1 - mu)/2,
    r1 = sqrt(xbar^2 + ybar^2), r2 = sqrt((xbar + 1)^2 + ybar^2).

Polar canonical change: xbar = r cos(vartheta), ybar = r sin(vartheta), p_xbar = p_r cos(vartheta) - (p_vartheta/r) sin(vartheta),
p_ybar = p_r sin(vartheta) + (p_vartheta/r) cos(vartheta). With r2 = sqrt(r^2 + 2 r cos(vartheta) + 1) (eq. 2.2):

    H = (1/2)(p_r^2 + p_vartheta^2/r^2) - p_vartheta - (1 - mu)/r - mu/r2 - mu (p_r sin(vartheta) + (p_vartheta/r) cos(vartheta)) - mu (1 - mu)/2.

Velocity-like variables **v = rdot r^(1/2), u = r^(3/2) vartheta-dot** and the time change **dt/dtau = r^(3/2)** (primes d/dtau) give
the regularised system (2.5) (TEXT p. 16):

    r' = v r
    vartheta' = u
    v' = v^2/2 + u^2 - (1 - mu) + 2 u r^(3/2) + r^3 + mu r^2 cos(vartheta) - mu r^2 (r + cos(vartheta))/r2^3
    u' = -u v/2 - 2 v r^(3/2) - mu r^2 sin(vartheta) (1 - 1/r2^3)

Energy relation (2.6): 0 = -h r + (1/2)(u^2 + v^2) - (1 - mu) - r^3/2 - mu r^2 (1/2 + cos(vartheta)... as printed in the text layer
with a garbled bracket; the 2018 paper prints the same relation as eq. (12) with the mirrored angle, and the invariant-manifold
statement below does not depend on the O(r) terms. Treat the exact grouping of the mu r terms as TEXT, not re-derived.

**Collision manifold** (eq. 2.7-2.8): Lambda = {r = 0, u^2 + v^2 = 2 (1 - mu)}, a torus (vartheta in [0, 2 pi]), independent of h. The flow on it:
vartheta' = u, v' = u^2/2 (using v^2/2 + u^2 - (1-mu) = u^2/2 on Lambda), u' = -u v/2. Two circles of equilibria
S+ = {r = 0, u = 0, v = +v0} and S- = {r = 0, u = 0, v = -v0}, **v0 = sqrt(2 (1 - mu))**. Linearisation M+- (matrix printed, TEXT p. 17) has eigenvalues
lambda = +-v0 (double), 0 (along the circle), -+v0/2, eigenvectors v1 = (0, 0, 1, 0), v2 = (1, 0, 0, 0), v3 = (0, 1, 0, 0),
v4 = (0, -2/v0, 0, 1) in state order (r, vartheta, v, u). So each P+ in S+ has a 2D unstable and a 1D stable manifold; each P- in S- a 2D stable and
a 1D unstable manifold.

**How ejection and collision with P1 are represented (p. 18).** Ejection orbits = W^u(S+): r > 0 for all finite tau and tends to a
point of S+ as tau -> -infinity (ejection takes INFINITE regularised time). Collision orbits = W^s(S-). An ejection-collision (EC)
orbit is a heteroclinic connection W^u(S+) intersect W^s(S-). For mu = 0 the EC condition is M = r^2 (dvartheta/dt + 1) = 0, i.e.
u = -r^(3/2), and the manifold at energy h is r h = v^2/2 - 1 (eq. 2.10): for h < 0 the ejection and collision manifolds coincide, so
every ejection orbit is an EC orbit (the Kepler radial case); for mu not equal to 0 the small primary deforms them apart and EC orbits
become isolated.

### 2.2 Levi-Civita regularisation (pp. 19-24; derivation TEXT, eq. 2.26 and U READ p. 22, polynomial forms re-derived symbolically)

General form: for xddot + B xdot = grad Omega(x), singularity at x = p, change of time **dt = a r ds** (a constant, r = |x - p|), then
**x = L(u) u + p**, L(u) = [[u, -v], [v, u]], so r = u^2 + v^2 = (u, u), x' = 2 L(u) u', and the result (eq. 2.22) is

    u'' + a L^T B L u' = (a^2/4) grad_u [ (u, u) U ],    U = Omega - C/2,

using the Jacobi integral 2 Omega - C = (xdot, xdot) = 4 (u', u')/(a^2 (u, u)). For the RTBP B = [[0, -2], [2, 0]], and the thesis
chooses **a = 4** and **p = (mu, 0)** for P1 (p = (mu - 1, 0) for P2). The 3D analogue is the Kustaanheimo-Stiefel operator L(u)
(4x4, eq. 2.23), cited to Stiefel & Scheifele 1971; the thesis does not use it (chapter 7 uses a 3D McGehee form instead).
Compact complex form (Szebehely): z = f(w), dt/ds = |f'(w)|^2, f(w) = p + w^2 for Levi-Civita.

**The regularised RTBP (eq. 2.26, READ p. 22):**

    u'' - 8 (u^2 + v^2) v' = d/du [4 U (u^2 + v^2)]
         = 4 mu u + 16 mu u^3 + 12 (u^2 + v^2)^2 u + 8 mu u / r2 - 8 mu u (u^2 + v^2)(u^2 + v^2 + 1)/r2^3 - 4 C u
    v'' + 8 (u^2 + v^2) u' = d/dv [4 U (u^2 + v^2)]
         = 4 mu v - 16 mu v^3 + 12 (u^2 + v^2)^2 v + 8 mu v / r2 - 8 mu v (u^2 + v^2)(u^2 + v^2 - 1)/r2^3 - 4 C v

with (READ p. 22)

    U = (1/2) [ (1 - mu)(u^2 + v^2)^2 + mu ((1 + u^2 - v^2)^2 + 4 u^2 v^2) ] + (1 - mu)/(u^2 + v^2) + mu/r2 - C/2,
    r2 = sqrt((1 + u^2 - v^2)^2 + 4 u^2 v^2),    ' = d/ds,  dt = 4 (u^2 + v^2) ds,
    x = u^2 - v^2 + mu,  y = 2 u v.

I verified symbolically (sympy, random parameter draws, residual below 5e-15) that the printed polynomial right-hand sides equal the
gradient of 4 U (u^2 + v^2), so the printed eq. 2.26 is self-consistent, and numerically that it reproduces the core CR3BP
equations (section 6). Here C is the thesis Jacobi constant (includes mu (1 - mu)).

Properties (p. 22-23):
- Jacobi integral (2.28): **u'^2 + v'^2 = 8 (u^2 + v^2) U**; at the primary (u = v = 0): **u'^2 + v'^2 = 8 (1 - mu)** (2.29), independent of C. The
  collision speed in regularised variables is fixed, so only the direction is free.
- Symmetries (2.27): (s, u, v, u', v') -> (-s, u, -v, -u', v') [from (1.6)] and (-s, -u, v, u', -v') [from the double cover].
  The map (u, v) -> (x, y) is two-to-one; every equilibrium is duplicated (collinear points lie on both axes of the (u, v) plane).
- Hill region (2.30): (u^2 + v^2) U >= 0. The system is regular everywhere except at P2 (r2 = 0), so this local regularisation
  suffices only while the Hill region does not reach P2, i.e. **C >= C_L1** (p. 23). For C below C_L1 the thesis (chapter 6) uses two local
  Levi-Civita charts, one per primary, with the global synodic chart between them and "three different times" (t, s1, s2); see Figure 2.5
  (mu = 0.2, C = 3.8).
- Relation between the McGehee and Levi-Civita angles: **vartheta_0 = 2 theta_0** (p. 35).

Hill problem (eq. 2.31-2.36, p. 24): same construction with p = (0, 0); u_h'' - 8(u_h^2 + v_h^2) v_h' = 4 U_h-gradient, U_h = 3 (u_h^2 - v_h^2)^2/2 + 1/(u_h^2 + v_h^2) - K/2 (as
printed in the text layer: "3 (u_h^2 - v_h^2)^2 / 2 + 1/(u_h^2 + v_h^2) ... - K/2", TEXT; I did not derive the Hill
polynomial); collision speed u_h'^2 + v_h'^2 = 8 (2.35); equilibria L1 = (+-3^(-1/6), 0), L2 = (0, +-3^(-1/6)) in (u_h, v_h).

### 2.3 Global regularisations (pp. 25-30; summary only)

Primaries put at (+-1/2, 0) with q = z + 1/2 - mu: q'' structure dt/ds = |f'(w)|^2. Thiele-Burrau f(w) = cos(w)/2 (periodic, r1 =
(cosh v - cos u)/2, r2 = (cosh v + cos u)/2, |f'|^2 = r1 r2/2, eq. 2.40-2.41); Birkhoff f(w) = (2w + 1/(2w))/4 (adds a singularity at w = 0
corresponding to infinity, eq. 2.42-2.44); Lemaitre f(w) = (w^2 + 1/w^2)/4 (eq. 2.45); and a general family f = (h(w)... ) covering Wintner,
Broucke generalisations. Collision speeds (p. 30): |w'|^2 = 2 (1 - mu) |h'(w1)|^2 at P1, 2 mu |h'(w2)|^2 at P2. The thesis chooses NOT to use them
(Birkhoff and Lemaitre polynomials are of higher degree than Levi-Civita) and uses two local Levi-Civita charts for C < C_L1.

## 3. Numerical method for n-EC orbits (chapter 3, pp. 31-41)

### 3.1 Definition (Definition 3.1, p. 31)

An **n-ejection-collision orbit (n-EC)** of a primary is the orbit that the particle describes when it ejects from the primary and reaches
**n times a relative maximum of the distance to that primary** before colliding with it. Equivalent statement (p. 57): it collides at its
n-th relative minimum of distance (counting the collision itself as the n-th minimum, i.e. n - 1 genuine close passages: Figure 3.1 caption,
"for n = 2 (n = 3) there are 1 (2) close passages to collision between ejection and collision"). The definition is local: it fails once
families are continued to low C (section 5.3 below: orbits of the alpha1 and gamma1 families acquire 2 and 3 maxima).

### 3.2 Initial conditions for ejection / collision

**McGehee (p. 32-34).** Ejection takes infinite tau, so start near the equilibrium on the unstable manifold. For P = (0, vartheta_0, v0, 0) in S+,
at energy H = h, the tangent plane of W^u(P) is spanned by (1, 0, 0, 0) and (0, 0, 1, 0); the energy-level normal is n = (-h - 3 mu/2, 0, v0, 0),
so the unit tangent vector is **w1 = (v0, 0, h + 3 mu/2, 0) / sqrt((h + 3 mu/2)^2 + v0^2)** (eq. 3.1) and the linear initial condition is
**(0, vartheta_0, v0, 0) + s w1** with s typically 1e-6 to 1e-8 (eq. 3.2). Collision orbits: the mirror point on S-, integrated backward.
To avoid the long integration time of small s the thesis builds a high-order parameterisation of W^u(P+) (parameterisation method of
Cabre-Fontich-de la Llave): invariance equation F(W(s)) = DW(s) f(s), expansion in powers of **s^(1/2)** (the system does not admit an
integer power series), normal-form choice f(s) = v0 s + ..., first terms printed to order s^(7/2) and s^4 (p. 34), used "usually" at order 5-10.
The 2018 paper used only the linear form with s = 1e-6 (tested 1e-7 to 1e-5).

**Levi-Civita (p. 35, eq. 3.9, READ).** Start exactly at the primary: **(u, v, u', v') = (0, 0, 2 sqrt(2 (1 - mu)) cos(theta_0), 2 sqrt(2 (1 - mu)) sin(theta_0))**,
theta_0 in [0, 2 pi), integrated forward (ejection) or backward (collision); because the (u, v) plane double-covers the configuration plane,
theta_0 in [0, pi) is enough. Ejection and collision occur in FINITE s. This is the key practical advantage over McGehee.
(Implementation note from my check: starting exactly at (0, 0) is a removable singularity of nothing in the printed system, but the
equations contain 1/(u^2 + v^2) only inside U, whose product with (u^2 + v^2) is regular; I started at s0 = 1e-9 along the straight line to avoid a literal 0/0 in a code that
evaluates U separately.)

### 3.3 Sections (p. 35-36)

McGehee: Sigma_M = {v = 0, v' < 0} (local maximum of distance), Sigma_m = {v = 0, v' > 0} (local minimum).
Levi-Civita (eq. 3.10): **Sigma_M = {h = u u' + v v' = 0, h' < 0}, Sigma_m = {h = 0, h' > 0}** (h is half the derivative of the regularised radius u^2 + v^2).
D_k+- (d_k+-) = k-th intersection of the ejection (collision) manifold with Sigma_M (Sigma_m).

### 3.4 The three methods (pp. 35-38) and their comparison (pp. 39-42)

- **Method I, manifold intersection.** Integrate W^u(S+) forward and W^s(S-) backward to a section and intersect the curves. Lemma 1
  (Levi-Civita): |D_i+ intersect D_j-| = 2 x (number of (i + j - 1)-EC orbits); Lemma 2: |d_i+ intersect d_j-| = 2 x (number of (i + j)-EC orbits) (factor
  2 from the double cover). Using symmetry (1.6) only D_(k)+ need be computed: for n = 2k - 1 use D_k+ (and obtain D_k- by symmetry), for n = 2k use d_k+.
  An EC orbit is **symmetric** if the ejection orbit meets the section with (x, y = 0, xdot = 0, ydot), else its mirror image is also an EC orbit.
  In McGehee variables only D_1- is well defined (d_i+- and D_j- for j > 1 hit the heteroclinic connections), so Method I is practical only for n = 1 there.
- **Method II, angular momentum (new in the thesis).** At the n-th crossing of Sigma_m the angular momentum about P1 vanishes iff the orbit
  is heading to collision: M_xy = (x - mu) ydot - y xdot = r^2 vartheta-dot = r^(1/2) u (eq. 3.11); in Levi-Civita **M_LC = u v' - v u'** (eq. 3.12), not
  equal to M_xy but with the same zeros. **n-EC orbits are the zeros of M_n(theta_0) := M at the n-th Sigma_m crossing.** The sign change of M_n
  detects the orbit; bisect.
- **Method III, singularity in time (the 2018 paper's method).** In McGehee variables the EC orbit is a heteroclinic connection that takes infinite time, so
  the time tau to reach the n-th Sigma_m crossing has a vertical asymptote at an EC theta_0; equivalently u changes sign there (u = M_xy / r^(1/2)).
  Singularities in time accumulate with n (the second crossing shows 1-EC and 2-EC orbits, ...): new ones appear only when comparing crossing n - 1 and n;
  practical up to n of about 10 (Figure 3.8 shows the 25th crossing, a dense mess).
- **Comparison (3.3).** Method I is half the integration time (symmetry) but needs an extra intersection step; II and III integrate about twice as far but read
  n-EC orbits off a single scalar. For n = 1 all three agree. **Massive computations use Levi-Civita + Method II** ("much lower computational cost"; finite
  time, start at the collision; the price is that the collision manifold information is lost). McGehee is intuitive (gives polar coordinates, handles
  multiple collision in other problems) but forces dealing with close approaches to the collision manifold.

### 3.5 Root finding, accuracy, search procedure for C-hat

- Scan theta_0 in [0, pi) for fixed (mu, C), bracket sign changes of M_n, refine by bisection (the 2018 paper: Newton to hit the section, bisection on u = 0).
  Remark 3 (p. 51): below C_L1 the particle can be near P2, so the thesis "numerically verified the condition u = v = 0" at the n-th crossing.
- Accuracy statements: double precision throughout; the thesis gives NO tolerance for integrator or root, and NO convergence study (the 2018 paper tested s from 1e-7 to 1e-5).
- Threshold C-hat(mu, n) (p. 53-54): fix a large C_b, find the four expected roots of M_n(theta_0) for each C from C_b down; at the first C at which more than four roots
  appear, refine C: that C is C-hat(mu, n) (the frontier before new families appear). The thesis states that for all mu in (0, 1) and **n from 1 to 100**
  there is a C-hat(mu, n) such that for C >= C-hat there are exactly four n-EC orbits (p. 48), and that C-hat(mu, 1) < C_L1(mu) (so n = 1 is "not considered").
  The analytic proof of the same statement is chapter 5 (part B). Plots, not tables: C-hat(0.1, n), n = 2..20 (Figure 3.26); C-hat(mu, n), mu in (0, 1), n = 2..10 (Figure 3.27).

## 4. Results printed in chapter 3 (pp. 42-65)

### 4.1 1-EC orbits

- Theorem 1 (Chenciner & Llibre 1988, quoted): for each mu in (0, 1) there is C-hat(mu) such that for C >= C-hat there are exactly four 1-EC orbits (p. 42).
  Corollary 3.4.1: two are symmetric about the x axis (called **alpha_1, gamma_1**), two are mirror images of each other (**beta_1, delta_1**).
- Numerical extension: "a complete exploration for mu in (0, 1) and C >= C_L1(mu)" shows four 1-EC orbits for all of it, so the numerical C-hat(mu) = C_L1(mu) (p. 43).
- Limit angles as C -> infinity (p. 43): the Levi-Civita angle **theta_0 -> 0, pi/4, pi/2, 3 pi/4 for gamma_1, delta_1, alpha_1, beta_1**; faster for larger mu (Figure 3.12 at
  mu = 0.1, 0.5, 0.8, 0.999; Figure 3.13 over mu in (0, 1), C in [C_L1, 8]).
- Below C_L1 (mu = 0.5, C in [C_L2, 6], Figures 3.14-3.18): new bifurcated families **eta_1, xi_1**; eta_1 symmetric, born first as C decreases; xi_1 non-symmetric, born
  from eta_1; eight 1-EC orbits exist at mu = 0.5 and C = C_L2. Regions A, B, C of the (theta_0, C, tau) diagram are described at C = 3.835 (Figure 3.15).
  (For mu = 0.5, C_L1 = 4.25 and C_L2 = 3.706796224086153 in thesis convention; those two constants are PRINTED only in the 2018 paper as H values, section 8.)

### 4.2 n-EC orbits

- Four orbits alpha_n, gamma_n (symmetric), beta_n, delta_n (mirror pair) for C >= C-hat(mu, n); same limit angles as n = 1; the limit is reached more slowly as n grows (Figure 3.20, n = 1, 2, 5, C in [5.5, 20]).
- Bifurcation diagrams of M_n(theta_0) over (theta_0, C), C in [C_L1, 8]: mu = 0.1 (Figure 3.21) and mu = 0.8 (Figure 3.22), n = 1..8. At mu = 0.8 C-hat is smaller, and for n = 2, 3, 4 it is below C_L1.
- Two kinds of bifurcation at C-hat (p. 51-53): **(a) n not equal to 3: four to six** n-EC orbits (a tangency of M_n with zero near alpha_n creates a pair, one the mirror of the other); **(b) n = 3: four to eight** (two tangencies, each giving two families).
- **Confluence**: beta_n and delta_n always end by merging into gamma_n for large enough n (Figure 3.28, n = 5, mu = 0.1).

### 4.3 Evolution at low C (section 3.4.3, pp. 56-65), mu varied

- mu = 0.0001 and 0.01 (Figure 3.29): beta_1 and delta_1 converge on gamma_1; alpha_1 and gamma_1 approach but do not merge and approach a periodic orbit (Figure 3.30), so
  the n-EC definition is not stable along a family.
- mu = 0.1 (n = 1, 2, 3, Figure 3.31); mu = 0.2 (Figure 3.32): continuation until collision with P2; delta_1 deformed near C = 2 by L4; gamma_2 deformed by LPO1 (Figure 3.33, C = 3).
- mu = 0.23, 0.235, 0.24 (Figure 3.34), mu = 0.3 (Figures 3.35-3.38): L4, L5 split gamma_1 into gamma_1 and gamma_1'; limit orbits are the heteroclinic chain ejection -> L4 -> L5 -> collision
  (Figures 3.36-3.37). In Figure 3.35 the families alpha_1 and gamma_1' are continued to the sixth collision with P2 (asterisks); the k-th such orbit of alpha_1 looks like the (k-1)-th with one extra exterior revolution (Figure 3.38).
- mu = 0.62, 0.64, 0.66 (Figure 3.39), mu = 0.7 (Figures 3.40-3.42, C = 3): the same splitting for beta_1 and delta_1 through L4/L5-collision connections.

## 5. Every printed number usable as a sourced test

All values are in the THESIS convention (primary of mass 1 - mu at +mu, C including mu (1 - mu)). Convert to core with mu_c = 1 - mu_T and
C_c = C_T - mu_T (1 - mu_T). "Reproduced" means my own check (section 6); "not reproduced" or "not attempted" are stated.

| # | Quantity | Printed value | Where | Status |
|---|---|---|---|---|
| T1 | C-hat(0.1, 2), tangent bifurcation of M_2(theta_0) from alpha_2 | **3.72442505** ("approximately") | p. 52, Fig. 3.24 caption (TEXT, digits plain) | my tangency C = 3.72441077 (theta_0 = 2.029024 rad = 0.645858 pi), difference 1.4e-5 (4e-6 relative); root count 4 at C = 3.76 and 6 at C = 3.69 REPRODUCED (my brackets: 3 roots at 3.7240 and 3.7243 in the window theta_0/pi in [0.640, 0.652], 1 at 3.72442505 and above, grid step 6e-5 rad) |
| T2 | C-hat(0.1, 3), 4 to 8 orbits | **3.80644009** | p. 53, Fig. 3.25 caption (TEXT) | root counts REPRODUCED (4 at C = 3.9, 8 at C = 3.7); the tangency value itself not attempted |
| T3 | Root counts mu = 0.1, n = 2 | 4 at C = 3.76; 6 below C-hat; Fig. 3.24 uses C = 3.76, 3.72442505, 3.69 | p. 52 | REPRODUCED (my theta_0/pi at C = 3.76: 0.142143, 0.397921, 0.64346, 0.880051; at C = 3.69: 0.146612, 0.400764, 0.626981, 0.647728, 0.669379, 0.88601; the printed text gives no angles, so these are my values, not sourced) |
| T4 | Root counts mu = 0.1, n = 3 | 4 at C = 3.9, 8 at C = 3.7 (Fig. 3.25 caption) | p. 53 | REPRODUCED |
| T5 | Four 1-EC orbits exist at C >= C_L1(mu) for every mu in (0, 1); limit theta_0/pi -> 0, 1/4, 1/2, 3/4 | counts and limits | pp. 42-43 | count 4 REPRODUCED at (mu = 0.1, C = 5): theta_0/pi = 0.044676, 0.319449, 0.545514, 0.768041 (offsets from the limits 0.045, 0.069, 0.046, 0.018; my values) |
| T6 | Fig. 3.30, mu = 0.01, theta_0 = 2 rad, four EC orbits at **C = 2.472170770645, 1.970463731686, 1.970412419219, 1.970412407739** | printed to 13 digits | p. 56 caption (READ) | M_1(theta_0 = 2) = 7e-14 at the first and 3.5e-9 at the second (so those two are 1-EC); at the third M_2 = 2.1e-5 (2-EC expected: near zero relative to a typical 1e-2 but not at 1e-9, sensitivity of a near-periodic orbit to the 12th digit of C); at the fourth M_3 = 7e-3, NOT reproduced as a 3-EC. Treat the last two as unconfirmed |
| T7 | C_L1 series | C_L1 = 3 + 9 ((1 - mu)/3)^(2/3) - 7 (1 - mu)/3 + O((1-mu)^(4/3)) | eq. 1.11 | numerical C_L1 minus series = 4.5e-5, 1.4e-5, 8.7e-7 at mu = 0.99, 0.999, 0.9999 (consistent with the stated remainder; sign positive) |
| T8 | K_L (Hill problem) | 3^(4/3); C = 3 + (1 - mu)^(2/3) K | eq. 1.18-1.19 | algebra only, not tested |
| T9 | mu_R (Routh) | 0.03852089650455137181950 | p. 8 | not tested (closed form (1 - sqrt(69)/9)/2) |
| T10 | Collision speeds | u'^2 + v'^2 = 8 (1 - mu) (Levi-Civita), 2 (1 - mu) (McGehee, u^2 + v^2) | eq. 2.29, 2.7 | Levi-Civita value used as the ejection speed in every check |
| T11 | Theorem-level statements usable as qualitative tests | exactly four n-EC orbits for C >= C-hat(mu, n) (n from 1 to 100 numerically); C-hat(mu, n) -> 3 as mu -> 1 (p. 55); C-hat(mu, n) increases with n | pp. 48, 54-55 | the mu -> 1 limit is stated from Figure 3.27, not tabulated |
| T12 | Figure parameters (for reproducing the pictures, not numbers to test) | Fig. 3.1: mu = 0.2, C = 4.25, n = 1, 2, 3; Fig. 3.3-3.5, 3.7: mu = 0.1, C = 5; Fig. 3.10: mu = 0.1, 0.5, 0.8, C = 5; Fig. 3.11: mu = 0.2, C = C_L1, 5, 7; Fig. 3.19: mu = 0.2, C = 10, n = 2, 4, 8; Fig. 3.28: mu = 0.1, n = 5; Fig. 3.14-3.18: mu = 0.5 | various | not tested |

The 2018 paper additionally prints three constants that are not in the thesis (section 8), and they ARE reproducible to 15 digits.
There is no printed initial angle theta_0 for any named orbit except theta_0 = 2 rad in Figure 3.30. A sourced golden for an EC orbit therefore
rests on T1-T6; the angles in T3 and T5 are mine and must not be used as expected values.

## 6. The checks I ran (scratch only)

Setup: a Levi-Civita integrator written from the printed eq. 2.26 (DOP853, rtol 1e-12 to 1e-13, atol 1e-14 to 1e-15), the ejection initial condition of eq. 3.9, section
Sigma_m = {u u' + v v' = 0, increasing}, M_n = u v' - v u' at the n-th crossing, roots by scan plus bisection.

1. **Printed equations self-consistent.** The printed polynomial/rational right-hand sides of eq. 2.26 equal the symbolic gradient of 4 U (u^2 + v^2) with the printed U (random draws of u, v, mu, C; residuals about 1e-15).
2. **Convention mapping with core.cr3bp.** Generic state at mu_T = 0.1, C_T = 5: mapped to (x, y, xdot, ydot) with x = u^2 - v^2 + mu_T, y = 2 u v, xdot = x'/(4 (u^2 + v^2)), mu_c = 0.9;
   `jacobi_constant` (core) = 4.909999999999998 against C_T - mu_T (1 - mu_T) = 4.91 (difference 8.9e-15 to 4.4e-14 over several states). Core `cr3bp_eom`
   (solve_ivp DOP853, rtol 1e-13) against Levi-Civita over regularised intervals 0.01, 0.05, 0.2 (physical time 0.0138, 0.0728, 0.2995): maximum state difference 3.1e-16, 2.4e-15, 1.8e-14.
   (An earlier attempt that straddled the collision point compared core against a trajectory that passes within 1e-27 of the primary and returned a meaningless difference of 4.6e3; that is a property
   of the test interval, not of the equations.)
3. **EC closure.** The four 1-EC angles at (mu = 0.1, C = 5) (T5) each give M_LC = 1e-13 at the first minimum crossing and, for theta_0/pi = 0.319449, |(u, v)| = 4.4e-14 there, i.e. the orbit returns to the primary
   (physical distance about 1e-27) at regularised time s = 0.7226 (about 0.4 time units in the units of the problem).
4. **Counts and bifurcation values** T1-T4, T6 as in the table.
5. **C_L1 series** (T7), and C_L1 for Earth-Moon (section 1.3).

Positive control: the root count 4 -> 6 for n = 2 and 4 -> 8 for n = 3, in the right C ranges and with the right type of change, is the published behaviour, so the integrator,
section and root-finder reproduce what the thesis computed; the one numerical value that is off in the 6th digit (T1, 1.4e-5) is a tangency whose location is poorly conditioned (the thesis prints it as "approximately").
The thesis gives no precision, so this is NOT a demonstrated error on either side.

## 7. Techniques applicable to the project's problems

**#928 (a regularised propagator and transition matrix for close passes).**
- The thesis eq. 2.26 is a complete, self-consistent planar Levi-Civita CR3BP with the right-hand side in closed form (no 1/(u^2 + v^2) singular term once multiplied through, as the printed
  polynomial form shows) and a=4. It is the same construction as the `#670` proof code. It adds the **Jacobi-constant bridge: C_T = C_c + mu (1 - mu)** and the mapping mu_c = 1 - mu_T, both verified (section 6).
- **A collision-orbit control for the regularised propagator** that is independent of the radial-fall closed form: for mu_T = 0.1, C_T = 5 the four 1-EC orbits exist (T5; my angles to 1e-13 in M_LC) and each, started
  at the primary with the eq. 3.9 velocity, must return to |(u, v)| below 1e-12 at its first minimum crossing; and for mu_T = 0.1 the root-count changes of T3 and T4 (4 -> 6 at C below about 3.7244 for n = 2; 4 -> 8 below about 3.806 for n = 3) are a published test of a propagator that survives many close passes
  (the n = 3 orbit passes close to the primary twice). Only the root COUNTS and the C-hat values (T1 to about 1e-5, T2 not attempted) are sourced; my angles are not.
- A tangency test: M_n(theta_0) = 0 and dM_n/dtheta_0 = 0 solved by Newton on (theta_0, C) located C-hat(0.1, 2) in a few evaluations. A transition matrix for the regularised flow would give dM_n/dtheta_0 analytically and make this a sensitive end-to-end gradient check (compare against a finite difference with h = 1e-6, which gave residual 1.7e-10).
- Practical: start at the collision itself (the thesis's own point against McGehee: finite regularised time, no near-singular start); a 3D version would need KS (the thesis cites Stiefel & Scheifele and in chapter 7 deliberately avoids it, see part B).

**#899 (second-species seeds at the Earth-Moon mass; relation to Henon's consecutive-collision arcs).**
- Mapping: Henon 1968 arcs are Kepler conics about M1 (the Sun/Earth) that meet the small body M2 twice (two collisions with M2 at a prescribed interval), pieces of Poincare's second-species orbits in the limit where the
  small body's mass and the distance of approach go to zero together. In the thesis the object is the EC orbit of the primary with mass 1 - mu_T; for a second-species skeleton around the MOON choose mu_T = 0.9878494157 (P1 = Moon, Hill problem as mu_T -> 1, chapter 5). EC orbits of the Moon are orbits that start and end at lunar collision: they are one-collision-to-the-same-collision arcs
  and are the n-EC analogue of Henon's two-collision arcs (Henon's arcs connect two collisions with the small body through a large Kepler ellipse of the big body; an EC orbit connects a collision to itself or its mirror through a Moon-centred loop plus Earth perturbation).
- Regime warning that matters for the use: the four-n-EC-orbit theorem and the families C-hat(mu, n) are for **C >= C_L1** (Hill region closed around P1): for the Moon that is C_T >= 3.2003 (core 3.1883), i.e. orbits that stay captured near the Moon; they do NOT
  contain Earth-Moon transfers. Earth-Moon second-species arcs (which must leave the lunar Hill region) live at C below C_L1, where the thesis shows only numerics (chapter 3 section 3.4.1 and 3.4.3: new bifurcated families, eight 1-EC orbits at mu = 0.5, C = C_L2) and chapter 6 (part B) for the transit structure around LPO1. So for #899 the thesis supplies (i) the regularised start at the Moon, (ii) a count-and-bifurcation control, (iii) the knowledge that the "four orbits" skeleton is lost exactly where transfers appear.
- The scalar M_n(theta_0) (angular momentum at the n-th minimum distance, equivalently the sign of the closest-approach side) is a cheap root-finding function for seeding Moon-collision arcs: zeros are collision seeds; nearby (nonzero) values are near-collision seeds, which is the use #899 needs ("continuation through near-collision seeds"). Continuation in C of the four families (initial angle versus C, Figures 3.12-3.13, 3.20) is the template.
- Henon relation, concretely: Henon's variable pair (half-interval tau, half-turn eta) picks a conic; the EC orbit has the integer n instead of tau, and the continuous parameter is theta_0 (the ejection direction) together with C. Both are one-parameter-per-integer families with two symmetric members (symmetric about the x axis) and mirror pairs.

**#906 (harden the demanded-turn gate).** The thesis is relevant only as a source of gate test cases, not of a criterion:
- Radial collision arcs (every EC orbit) end with an undefined turn, so the gate must return indeterminate for them, not reject; a regularised EC arc started at the collision (eq. 3.9) gives such a case with a sourced existence (T3-T5).
- Near-collision passes: for the n-EC orbits with n >= 2 the particle passes the primary closely n - 1 times (Figure 3.1, p. 31); at (mu_T = 0.1, C = 5) my 1-EC orbits reach the primary to physical distance about 1e-27, which is the limiting case of a zero-radius pass. A flyby chain built from such an orbit has a demanded turn that is 0 or 180 degrees, which the MacKay-2005 amendment (section on #906 in OUTSTANDING) says is not an encounter.
- Convention for any gate input that uses thesis-frame velocities: they are body-relative only if the body is P1 or P2 of the thesis frame; the frame is the project's rotated by pi (or equal to it with mu_c = 1 - mu_T), see 1.2.

## 8. Olle, Rodriguez & Soler 2018 compared with the thesis

Source: M. Olle, O. Rodriguez and J. Soler, "Ejection-collision orbits in the RTBP", Commun. Nonlinear Sci. Numer. Simulat. 55:298-315 (2018), DOI 10.1016/j.cnsns.2017.07.013. Accepted manuscript, 33 pages (received 7 February 2017, revised 10 July 2017,
accepted 18 July 2017; Highlights page plus 32 numbered pages). Filed in the private paper corpus as
`olle-rodriguez-soler-2018-ejection-collision-orbits-rtbp-cnsns-55-298-doi-10.1016-j.cnsns.2017.07.013.pdf`. Page numbers below are the manuscript's. This is the published form of **thesis chapter 3 plus section 2.1**.
The manuscript is the text layer (TEXT); the figures were not viewed. Where the text layer is garbled (equations 9, 11, 12) I report only what is legible.

**Conventions.**
- Mass parameter mu in **(0, 0.5]**, P1 = the BIG primary (mass 1 - mu) at **(-mu, 0)**, P2 (mass mu) at **(1 - mu, 0)** (p. 6). This IS `core/cr3bp.py`'s frame with mu_c = mu, no transformation. The thesis (mu in (0, 1), P1 = mass 1 - mu at +mu) differs by a rotation of pi about the
  barycentre, and by allowing P1 to be the small body. The paper's angle is theta; the thesis's McGehee angle is vartheta = theta + pi.
- Hamiltonian printed with **- mu (1 - mu)/2** (eq. 6), matching thesis eq. 2.1 and contradicting thesis eq. 1.5 (section 1.4). Jacobi C = 2 Omega - v^2 including mu(1 - mu) (eq. 2-3), H = -C/2 (eq. 5). **The paper works in H, not C**: every energy is H = -C/2, e.g. H = -5.25 is C = 10.5.
- The paper's Hill-region restriction is C >= C_L2 (p. 7) and the numerical studies of n-EC orbits use H below H_L1 (p. 24).

**What the paper has that the thesis does not print.**
- Three constants (manuscript pp. 15, 20 and 25, approximate, the text layer has no reliable page marks near them), all in the paper's H: **H_L1(0.5) = -2.125** (C = 4.25), **H_L2(0.5) = -1.853398112043077** (C = 3.706796224086154), **H_L1(0.1) = -1.843476614939948** (C = 3.686953229879896). I recomputed
  them from the collinear equation and core `jacobi_constant` plus mu(1 - mu): -2.125, -1.8533981120430765 and -1.8434766149399473; they agree to 15 digits. These are the only 15-digit sourced Jacobi constants for named (mu, point) pairs in either document, and a clean sourced test of the thesis-to-core convention (C_T = C_c + mu (1 - mu)).
- Energy levels of the figures: Figure 7, 9, 10: mu = 0.5, H = -5.25, -3.25, -2.75 and H_L1; Figure 11 diagram for H up to H_L2; Figure 17: mu = 0.1, H = -5.05, -3.05; Figure 18: mu = 0.1, n = 1..5, H <= H_L1(0.1).
- Regularisation-step parameters: s = 1e-6 typical, tested 1e-7 to 1e-5 "giving rise to the same results" (p. 15); integrator **Runge-Kutta-Fehlberg 8(7)**, double precision (p. 18); a Newton method to land on the Poincare section (p. 14); the grid mu in [0.01, 0.5] for the four-orbit existence study (p. 17) and "mu in [0.01, 0.5]" for n-EC (p. 24).
- Extended wording of introduction and review (pp. 3-4): Llibre 1982 (at least two EC orbits, small mu and large C), Lacomba & Llibre 1988 (no C^1-extensible regular integrals), Chenciner & Llibre 1988 (four EC orbits, mu in (0, 0.5], H small), Llibre & Martinez Alfaro 1985 (spatial), Llibre & Pinyol 1990 and Pinyol 1995 (elliptic), Henon 1965 (mu = 0.5) and 1969 (Hill), Bozis 1970 (16 collision periodic orbits followed in mu), applications (irregular-moon capture, Kuiper-belt binaries, hydrogen in a microwave field, rubble-pile shedding); the thesis only cites a subset ([ARV13], [CL88], [DF89], [LL88], [LMA85], [ML98], [Pin95]).
- Its abstract and conclusion state **n from 1 to 10** ("we study the cases 1 <= n <= 10", p. 2; "for n <= 10", p. 29), while section 4 describes results for 1 <= n <= 5 (the summary line on p. 5). The thesis states n up to 100 and mu in (0, 1).

**What the thesis adds beyond the paper.**
- Levi-Civita regularisation and all its numerics (chapter 2.2, 3.1.2); the Hill-problem version; the global regularisations; Method II (angular momentum) and the Lemma 1-2 counting; the McGehee parameterisation method; mu in (0, 1) (the paper stops at 0.5), n up to 100; the numerical thresholds C-hat(mu, n) (Figures 3.26-3.27) and the two bifurcation types (4 to 6 and 4 to 8) with the C-hat(0.1, 2) and C-hat(0.1, 3) values (T1, T2); the evolution study at mu = 0.0001 to 0.999 with the L4/L5 heteroclinic chains (section 3.4.3).

**Differences in numbers, definitions and wording.**
| Item | 2018 paper | Thesis chapter 3 |
|---|---|---|
| mu range | (0, 0.5] ("mu in [0.01, 0.5]" computed) | (0, 1) (figures at mu = 0.0001 ... 0.999) |
| Energy variable | H (levels such as -5.25) | C = -2 H (levels such as 10.5, 5.5) |
| First primary | the big one, always | mass 1 - mu, so big for mu < 0.5 and small for mu > 0.5 |
| Frame | big primary at (-mu, 0) | big (mass 1 - mu) at (+mu, 0): rotated by pi |
| McGehee angle | collision with P2 at theta = 0 | at vartheta = pi (eq. 2.5 text) |
| Linear eigenvector ordering (collision manifold) | v1 = (0, -2/v0, 0, 1), v2 = (0, 0, 1, 0), v3 = (1, 0, 0, 0), v4 = (0, 1, 0, 0); lambda1 = -v0/2, lambda2 = lambda3 = v0, lambda4 = 0 (listed for M+) | v1 = (0, 0, 1, 0), v2 = (1, 0, 0, 0), v3 = (0, 1, 0, 0), v4 = (0, -2/v0, 0, 1); lambda1 = lambda2 = +-v0, lambda3 = 0, lambda4 = -+v0/2 |
| Tangent vector for the ejection start | w = (1, 0, (H + 3 mu/2)/v0, 0), normalised, s typically 1e-6 | w1 = (v0, 0, h + 3 mu/2, 0)/norm (the same direction, scaled), s typically 1e-6 to 1e-8, or parameterisation order 5-10 |
| Methods | two: manifold intersection with section Sigma = {v = 0}; "singularity in time" T_2n; bisection on u = 0 at the 2n-th crossing of Sigma (both max and min) | three: I (D_k, d_k, Lemmas 1-2, McGehee and Levi-Civita), II (angular momentum at Sigma_m), III (time asymptote); massive runs in Levi-Civita with Method II |
| Section definition | one section {v = 0}; crossing 2n is the n-th minimum | Sigma_M (v' < 0) and Sigma_m (v' > 0) separate |
| Symmetric/non-symmetric classification | by theta_n = 0 or pi at the n-th section crossing | by x-axis crossing at the half-way section, Corollary 3.4.1 |
| Names of the low-C 1-EC families (mu = 0.5, H between H_L1 and H_L2) | **eta_1 and eta_2** (eta_1 symmetric, born at H = H_b1, eta_2 non-symmetric, born at H_b2 on eta_1) | **eta_1 and xi_1** (same roles, new letter xi) |
| Bifurcation description | in H increasing ("as H increases, a new branching point") | in C decreasing |
| Existence statement for n-EC | for all mu in (0, 0.5] and all n a H-hat(mu, n) exists such that H <= H-hat gives four n-EC orbits | for all mu in (0, 1) and n from 1 to 100 a C-hat(mu, n) exists such that C >= C-hat gives four |
| Collapse / bifurcation examples | beta_4 and delta_4 collapse on gamma_4 (Fig. 20) and alpha_4 bifurcates into alpha_41, alpha_42 (Fig. 21), mu = 0.1 | beta_5, delta_5 into gamma_5 (Fig. 3.28) and n = 2 (4 to 6), n = 3 (4 to 8) at mu = 0.1 |
| Reference list | [10] and [13] numbering duplicated (the manuscript lists two [10]) | thesis bibliography |
| Orbit-count language | "exactly four 1-EC orbits for any mu and H small enough" (Chenciner & Llibre) | the same, and "(Theorem 1) C-hat(mu)" |

The comparison turns up no numerical value that the paper prints and the thesis prints differently: they share no printed number except mu = 0.5 and 0.1 as example parameters; the thesis does not print the paper's three H constants.

## 9. Follow-ups (no task numbers registered)

1. A sourced golden for the convention: assert `jacobi_constant(L1 state, mu) + mu (1 - mu)` equals -2 H_L1 for the 2018 paper's three printed H values (mu = 0.5 L1, mu = 0.5 L2, mu = 0.1 L1; section 8); cheap, independent of the thesis.
2. Once #928's Levi-Civita propagator exists: add the n = 2 and n = 3 root-count tests at mu_T = 0.1 (T3, T4) and the C-hat(0.1, 2) tangency to about 1e-4, and the four-orbit counts at (mu_T = 0.1, C = 5).
3. The one off-digit comparison (T1, 3.72441077 against 3.72442505) deserves a tighter reproduction (higher-precision section crossing, Newton on the extended system with an analytic M_n derivative) before it is used as a sourced tolerance; the near-periodic orbit of Figure 3.30 (third and fourth values, T6) could be reproduced by continuing theta_0 = 2 with an arclength corrector rather than by evaluating at the printed 12-digit C.
4. Fix the sign slip in thesis eq. 1.5 and 1.23 in any note that copies it (the constant is minus mu(1 - mu)/2).
5. A Moon-primary run: compute the four 1-EC orbits and C-hat(0.9878494157, n) (Hill-problem limit) for the Earth-Moon mass and compare with the Hill-problem results of chapter 5 (part B); then test whether any four-orbit skeleton persists below C_L1 = 3.2003 (thesis convention) where Earth-Moon transfers live.
6. Add the check in section 6 item 3 (EC closure to |(u, v)| below 1e-12) to a Levi-Civita regularised-propagator test as the collision-orbit control requested for #928.
7. Cross-check part B for the C-hat analytic expression (chapter 5) against the numerical Figure 3.27 pattern the thesis says it matches ("C-hat -> 3 as mu -> 1").
8. A digest of Llibre & Martinez Alfaro 1985 (the spatial EC orbits) and of Chenciner & Llibre 1988 (Theorem 1 here) is still missing from `docs/notes/`; the thesis relies on both.
