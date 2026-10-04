# Digest: Villegas-Pinto, Baresi, Locoche & Hestroffer (2023), "Resonant quasi-periodic near-rectilinear Halo orbits in the Elliptic-Circular Earth-Moon-Sun Problem"

Advances in Space Research 71(1):336-354 (2023), DOI 10.1016/j.asr.2022.08.011, HAL accepted manuscript hal-04011123
(CC BY 4.0), 26 PDF pages (a HAL cover page plus 25 manuscript pages). Received 12 April 2022.
Filed in the private paper corpus as
villegas-pinto-baresi-locoche-hestroffer-2023-resonant-quasi-periodic-NRHO-elliptic-circular-earth-moon-sun-asr-71-336-doi-10.1016-j.asr.2022.08.011-hal-accepted-manuscript.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page, section, equation, figure or table) or
INFERRED (our reading, or a comparison with project code). Page numbers are the manuscript's own printed page numbers (1 to
25); the PDF page is the manuscript page plus 1. Line numbers in the margin of the manuscript are not used. Context: tasks
`#884` (Sun-perturbed Earth-Moon periodic orbits) and `#891` (the Sun's sense in `core/bcr4bp.py`); the project models in
question are `core/bcr4bp.py`, `core/qbcp.py` and `core/er3bp.py`.

## 0. What the paper is

READ (abstract, p1): the paper computes "dynamical substitutes of the Earth-Moon's resonant Near-Rectilinear Halo Orbits
(NRHO) under the Elliptic-Circular Restricted Four-Body Problem formulation of the Earth-Moon-Sun system. This model
considers that the Earth and Moon move in elliptical orbits about each other and that a third body, the Sun, moves in a
circular orbit about the Earth-Moon barycenter." Because of the Moon's eccentricity the periodic NRHOs of the circular
restricted three-body problem (CR3BP) "are hereby replaced by two-dimensional quasi-periodic tori". It treats three
resonances: the 9:2 synodic resonant L2 southern NRHO (the planned Lunar Gateway orbit), the 4:1 synodic and the 4:1 sidereal
resonant NRHOs (p3). It claims "the first attempt at calculating and analyzing quasi-periodic trajectories that take into
account both the perturbation of the Sun and the eccentricity of the Moon's orbit in the context of trajectory and mission
design" (p3).

## 1. The model as printed (section 2, pp3-6)

### 1.1 Setting

READ (pp3-4): the primary (Earth, mass m1) and secondary (Moon, mass m2) move on coplanar Keplerian ellipses about their
barycentre; the third body (Sun) moves on a circular orbit about that barycentre, in the same plane. The paper says (p3):
"just as the BCR4BP, the ECR4BP is not a coherent model. That is, the motion of their massive bodies is not a solution of the
Three-Body Problem since the Sun is considered to not affect the motion of the Earth and the Moon." It names the coherent
quasi-bicircular model (Andreu 1998) and the BiElliptic model as other restricted four-body models and does not use them.
The model reduces to the bicircular problem when `e = 0` and to the ER3BP when `mu_3 = 0` (p6: "when setting the
eccentricity to zero the equations simplify to the Bicircular case and when setting mu_3 to zero they simplify to the ER3BP
model").

### 1.2 Frame, independent variable, units

READ (p4): a co-rotating frame centred on the barycentre of the primary and secondary, "such that the x-axis constantly points
from the primary to the secondary, the z-axis points in the direction of their mutual orbit's angular momentum vector, and the
y-axis completes the orthogonal frame." Because the primaries' separation is not constant, the normalisation is that of the
ER3BP (Szebehely 1967): the frame is PULSATING and the independent variable is the true anomaly `nu` of the primaries'
orbit, not time. Units (eqs. 1, 2 and p5):

- Length unit `l = a (1 - e^2)/(1 + e cos nu)` (eq. 1), the instantaneous distance between the primaries.
- `nu_dot = (G M)^(1/2) (1 + e cos nu)^2 / (a^(3/2) (1 - e^2)^(3/2))` (eq. 2); time unit `[TU] = 1/nu_dot`.
- Mass unit `[MU] = m1 + m2`, `M = m1 + m2`.
- "Due to the normalization used, the primary and secondary are always on the x-axis at -mu and 1 - mu, respectively" (p6).

So the Earth is at x = -mu and the Moon at x = 1 - mu, the same as in `core/bcr4bp.py` (READ, project docstring: "Earth at
(-mu, 0, 0), Moon at (1 - mu, 0, 0)"; INFERRED agreement).

### 1.3 Equations of motion

Dimensional third-body acceleration (eq. 3, p5): `r3_ddot* = -(G m3*/r3*^3) r3* - (G m3*/a3*^3) [x3*, y3*, z3*]^T`, where
the text defines `m3*` ("`mu_3*`" in the print) as the third body's gravitational parameter, `a3*` as the semi-major axis of
the third body's orbit, `x3*, y3*, z3*` as its position in the co-rotating frame, and `r3*` as the spacecraft-Sun distance
(the printed definition of `r3*` is written as a squared sum; transcribed here without the square, as the text has it).
Normalised third-body acceleration, eq. 4 (p5):
`(G M / l^2)( -mu3 r3/r3^3 - (mu3/a3^3) [x3, y3, 0]^T ) = (G M / l^2) r3''`.

Normalised equations in the pulsating frame (eqs. 5a-5c, p5), with `'` denoting `d/d nu`:

    x'' - 2 y' = d psi / d x
    y'' + 2 x' = d psi / d y
    z''        = d psi / d z

    psi = 1/(1 + e cos nu) [ 1/2 (x^2 + y^2 - e z^2 cos nu) + (1 - mu)/r1 + mu/r2 + mu3/r3 - (mu3/a3^3)(x3 x + y3 y + z3 z) ]   (eq. 6)

with `r1 = sqrt((x + mu)^2 + y^2 + z^2)`, `r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2)` and `r3` the spacecraft-Sun distance. Equivalent
form (eqs. 7a-7c and 8, p5): `x'' - 2y' = d psi_tilde/dx`, `y'' + 2x' = d psi_tilde/dy`, `z'' + z = d psi_tilde/dz`, with
`psi_tilde = 1/(1 + e cos nu) [ 1/2 (x^2 + y^2 + z^2) + (1 - mu)/r1 + mu/r2 + mu3/r3 - (mu3/a3^3)(x3 x + y3 y + z3 z) ]`.
Here `mu3 = G m3*/(G M)` is the Sun's mass in units of the Earth-Moon mass and `a3 = a3*/l` is the Sun's distance in the
pulsating unit, so it varies with `nu` (INFERRED from the definitions: `a3*` is constant, `l` is not). The term
`-(mu3/a3^3)(x3 x + y3 y + z3 z)` is the indirect (barycentre acceleration) term (READ as printed; INFERRED identification).

### 1.4 The Sun's position and sense (eqs. 9-14, p6)

READ, eq. 9: `x3 = a3 cos sigma = (a3*/l) cos sigma`, `y3 = -a3 sin sigma = -(a3*/l) sin sigma`, `z3 = 0`, where `sigma` is
"the angular position of the third body in the pulsating frame". Eq. 10: `sigma = (nu - n3 t) + sigma_0`. Eq. 11:
`sigma' = d sigma/d nu = 1 - n3/nu_dot = 1 - n3 (1 - e^2)^(3/2)/(n (1 + e cos nu)^2)`, with `n = sqrt(G M/a^3)`; `n3` is the Sun's mean
motion. Two auxiliary variables `c = cos sigma`, `s = sin sigma` are appended (eqs. 12a, 12b: `c' = -s sigma'`, `s' = c sigma'`; eq. 13
`chi' = sigma' J chi`, `chi = [c, s]^T`, `J = [[0, -1], [1, 0]]`) so that the system is `2 pi`-periodic in the Sun angle: the
state is `X = [r, v, chi]` (eight-dimensional) and `X' = f(X, nu, xi) = f(X, nu + 2 pi, xi)` with parameter vector
`xi = [mu, e, mu3]^T` (eq. 14).

Sense of the Sun's motion, READ: with `y3 = -a3 sin sigma` and `sigma` increasing with `nu` (for `n3 < n`, `sigma' > 0`), the
Sun at (a3 cos sigma, -a3 sin sigma) moves clockwise in this frame. INFERRED: this is the same sense as the corrected
`core/bcr4bp.py` (Sun at `(a_S cos theta, a_S sin theta)` with `theta = theta_S0 - omega_S t`, clockwise) and as Jorba-Cusco,
Farres & Jorba 2018 (digest `2026-10-04-digest-jorba-cusco-farres-jorba-2018-two-periodic-models.md`), and it supports the
`#891` correction. The paper does not remark on the sense in words; it is only visible in the sign of eq. 9 and in Fig. 1,
which draws the Sun in the lower right quadrant (x > 0, y < 0) with `sigma` measured clockwise by an arrow (p4, read from
the figure; INFERRED that the arrow is clockwise from its position).

The initial Sun angle used for the synodic resonant orbits is zero: "when all bodies are aligned and the Sun is in the
positive x-axis" (p14). For the sidereal ones the initial true anomaly is zero (p14). Two phase constraints fix `sigma_0` at 0
or `pi`: `p_c(X) = c_0 - cos(sigma_0)`, `p_s(X) = s_0 - sin(sigma_0)` (eqs. 21, 22, p10).

### 1.5 Constants as printed

READ, parameters in the text (p6, section 2, last paragraph): CR3BP Earth-Moon system `mu = 0.01215`, `e = 0.0`, `mu3 = 0.0`;
`e = 0.0549`, `mu3 = 0.0` for the ER3BP; `e = 0.00`, `mu3 = 328900.55` for the BCR4BP. The final value of the continuation in
eccentricity is `e = 0.0549` and the Sun-mass homotopy parameter `epsilon` runs to 1 "corresponding to the actual Sun mass"
(p16).

READ, Table 1 (p13), "Physical constants and parameters used for the different models. GM refers to the gravitational
parameter, R refers to the body's radius, L to the distance or semi-major axis between two bodies. Sources: Wieczorek et al.
(2006); IAU Division Working Group (2022); Mamajek et al. (2015)":

| Parameter | Value |
| --- | --- |
| GM Earth [m^3 s^-2] | 3.9860044189 x 10^14 |
| GM Moon [m^3 s^-2] | 4.902801076 x 10^12 |
| GM Sun [m^3 s^-2] | 1.32712440018 x 10^20 |
| R Earth [km] | 6378.137 |
| R Moon [km] | 1737.103 |
| R Sun [km] | 6.9551 x 10^5 |
| L Earth-Moon [km] | 384 399 |
| L Sun-Earth [km] | 149.5978707 x 10^6 |

The Sun's normalised distance `a3` and mean motion `n3` are not printed as numbers; Table 1 gives "L Sun-Earth", the
Earth-Sun distance, not the Sun to Earth-Moon-barycentre distance (READ; INFERRED that the paper's `a3*` is the Earth-Sun
distance, as the table label suggests; the paper does not say so in words). The project's BCR4BP uses `a_S` = 388.811... EM
distances (READ, `core/bcr4bp.py` docstring); the two are not compared numerically here.

### 1.6 Differences from the bicircular and quasi-bicircular models

READ (pp3-4): in the bicircular model the primaries are on circles, the system is periodic with the single synodic frequency
`Omega_3 = n - n3`. In the ECR4BP the system "represents a quasi-periodic system with two (incommensurate) frequencies, one
equal to Omega_3, same as for the Bicircular, and another being the rate of the true anomaly of the two primaries, which is
influenced by their eccentricity." One sentence: the ECR4BP differs from the bicircular problem in that the Earth and Moon
move on Keplerian ellipses of eccentricity `e = 0.0549` (so the frame pulsates and time is replaced by the true anomaly)
while the Sun is still on a circle, which makes the system quasi-periodic with two incommensurate frequencies and replaces
periodic orbits by two-dimensional tori. Against the quasi-bicircular problem (Andreu 1998; Jorba-Cusco et al. 2018), the paper
says only that QBCP is a coherent model and "we do not consider" it (p3). The QBCP keeps the system periodic (a single
frequency, the Sun's) and accounts for the Earth-Moon eccentricity and the Sun's pull on the primaries through Fourier
coefficients (INFERRED from the project's `core/qbcp.py` docstring and the 2018 digest, not stated in this paper).
The model's eccentricity `e = 0.0549` is a constant (the mean value), not the time-varying lunar orbit of the real ephemeris;
inclinations of the orbital planes are neglected, as the paper says on p24 ("the Elliptic-Circular model does not take into
account the inclination of each of the bodies' orbital planes").

## 2. What is computed

### 2.1 Resonances and the torus structure (sections 3, 4.2, 5; pp7-8, 12)

READ: "p:q" means p orbital periods in q synodic (Earth-Moon-Sun) or sidereal (Earth-Moon) months (p8: "p is the number of
orbital periods and q is the number of either synodic or sidereal months"). The three cases:

| Resonance | Meaning (p2, p8) | Intermediate periodic model | Torus frequencies (p12) |
| --- | --- | --- | --- |
| 9:2 synodic | 9 revolutions per 2 synodic periods ("The synodic period of the Earth-Moon-Sun is about 29.5 days", p2) | Bicircular, `epsilon` homotopy from the CR3BP | `omega_0 = 2 pi/T`, `T = q T_syn`; `omega_1 = 1` (the true-anomaly rate); `rho = T` |
| 4:1 synodic | 4 revolutions per synodic period | Bicircular | same as above |
| 4:1 sidereal | 4 revolutions per sidereal period | ER3BP, homotopy in `e` from the CR3BP | `omega_0 = 1/q` (`omega_0 = 1` for 4:1, `omega_0 = 0.5` for a 3:2 sidereal); `omega_1 = rho/T = (1/T) int_0^(2 pi) sigma' d theta_0` (eq. 37) |

For synodic resonant tori the angle `theta_1` reflects the Earth-Moon true anomaly and `theta_0` the Sun; for sidereal tori
`theta_0` reflects the true anomaly and `theta_1` the Sun angle (p12). A periodic orbit that is resonant with one of the
two perturbations gives a TWO-dimensional torus; a non-resonant one would give a three-dimensional torus (p8, a statement
the paper uses to motivate the choice).

No rotation number is printed as a number. Eq. 37 and the identity `rho = T` for synodic tori are given in general; the
periods `T` of the multi-revolution orbits (and the NRHO periods in days) are not tabulated. The only time values quoted
are the synodic period of about 29.5 days (p2) and "approximately 177 days" for `3 x 2 T_syn` (p22).

### 2.2 Algorithms (sections 3 and 4, pp7-11)

READ: tori are computed by the GMOS two-point boundary value algorithm (Gomez & Mondelo 2001; Olikara 2016; Baresi et al.
2018; Villegas-Pinto et al. 2021) with the invariance condition `[R_{-rho}] phi_T(X_0) - X_0 = 0` (eq. 16), where `[R_{-rho}]`
is a rotation operator built by the discrete Fourier transform (p7), multiple shooting with `N_0` nodes (eq. 18-19), phase
conditions (eqs. 20-22), the pseudo-arclength condition (eq. 23), full system `F = [G; p; q] = 0` (eq. 24), Newton iteration
to a tolerance of `10^-10` (p10), and the family step `z^{m+1} = z^m + (dz^m/dh) dh` from the null-space tangent (eq. 25).
Partial derivatives with respect to the Sun mass parameter `epsilon` (eq. 28, `B = 1/(1 + e cos nu) r3''`) and with respect to `e`
(eqs. 29-35, including the closed forms for `da3/de` and `d sigma'/de`) are given (p11).

The continuation is a two-step homotopy (pp8-9, 14-16), READ:

1. Start from the CR3BP `p:q` resonant NRHO and repeat it `p` times to get a multi-revolution periodic orbit with period
   `T = p T_res`.
2. Continue it to the intermediate periodic model: Sun mass `epsilon` from 0 to 1 (CR3BP to BCR4BP) for the synodic cases;
   eccentricity `e` from 0 to 0.0549 (CR3BP to ER3BP) for the sidereal case. "By design, synodic orbits remain periodic in
   the Bicircular model, whereas sidereal orbits remain periodic in the ER3BP" (p8). This step uses the same GMOS code with
   `N_1 = 1` and `N_0 = p` shooting nodes, one node at each apolune.
3. Continue the multi-revolution periodic orbit into the ECR4BP: eccentricity `e` from 0 to 0.0549 for the synodic cases;
   `epsilon` from 0 to 1 for the sidereal case, with `N_1` = 50 points on the invariant circle (p15: "We set the value of N_1
   to 50 ... can be set between 30 and 50 (or larger)"; lower-perilune orbits such as the 9:2 need more).

The initial guess for the first family member repeats each multiple-shooting state `N_1` times; the true anomaly (synodic) or the
Sun angle (sidereal) along the initial circle is `2 pi j/N_1` plus node offsets (eqs. 39-43, pp15-16). All nodes are placed
at apolune because perilune is highly sensitive (p13); the first node uses Sun angle zero (synodic) or true anomaly zero
(sidereal). Symmetries of the intermediate models, eqs. 38a, 38b (p14): `(x, y, z, x_dot, y_dot, z_dot, tau) -> (x, -y, -z, -x_dot, y_dot,
-z_dot, -tau)` and `-> (x, y, -z, -x_dot, y_dot, -z_dot, -tau)`, where `tau` is the Sun angle (bicircular) or the true anomaly
(ER3BP); "a sufficient condition for symmetric periodic orbits in these models is that they cross the x-axis or xz-plane
perpendicularly at either tau = 0 or tau = pi". The ECR4BP has no such symmetry used here.

Convergence statement (p12): "while all the resonances studied in this work converge correctly in the full Elliptic-Circular
model ..., many of the other synodic and sidereal resonant orbits fail to do so, even after being computed in the
intermediate models." The explanation offered is a bifurcation along the stability curve of the 3:1 and 5:1 synodic NRHOs in
the CR3BP to bicircular continuation (Boudad et al. 2020), which "seems to later prevent the numerical continuation to the full
Elliptic-Circular model ... or to ephemeris models"; "extensive testing was not performed".

### 2.3 Results printed

There is no table of initial conditions, periods, frequencies, rotation numbers or stability indices in the paper. The only
table is Table 1 (constants, section 1.5). Everything else is in figures. What the figures and text state:

- Figs. 3 and 4 (p14): the CR3BP resonant southern NRHOs and the multi-revolution periodic orbits in the intermediate models;
  Figs. 5 to 7 (pp17-18): the resonant tori at `e = 0.0549` with `epsilon = 1`. Axes in kilometres in the Earth-Moon rotating
  frame centred on the Moon; no numerical values are printed other than axis ticks (z axis about -60000 to 0 km).
- Figs. 8 and 9 (pp18-19): the initial invariant circles along the continuation. Legend values, as printed:
  9:2 synodic (a): `e` = 1.76e-03, 2.98e-02, 5.49e-02; 4:1 synodic (b): `e` = 1.74e-03, 2.96e-02, 5.49e-02; 4:1 sidereal (c):
  `epsilon` = 0.010, 0.520, 1.000 (legend symbol printed as a lunate epsilon in all three panels; for (a) and (b) it is the
  eccentricity and for (c) the Sun mass parameter, by section 5.2, p15; INFERRED). Tick values for the circles in Fig. 8 are in
  normalised units (for example 9:2 synodic: x about 0.016 to 0.0175, z about -0.19 to -0.175) and in Fig. 9 in kilometres
  (for example 9:2 synodic: x about 6000 to 6800, z about -70000 to -69200) (read from plot, approximate).
- Section 5.2 text (p18): the second perturbation produces a "thickening" of the solution around the periodic resonant
  NRHO, "most prevalent in the 9:2 synodic and 4:1 synodic resonant NRHOs"; the torus is seen "thicker" than the converged
  ephemeris solutions of the literature, which the authors say might be different dynamical solutions, "although analyses of
  these quasi-periodic solutions in a full-ephemeris model would be necessary for a more accurate comparison".
- Section 6, stability (pp19-21): finite Lyapunov exponents `phi_i = Re(ln(lambda_i)/T)` (eq. 44) from the Floquet matrix of
  the torus, `B = [R_{-rho}] diag(Phi_0, ..., Phi_{N_1-1})` (eq. 45), `T` the fundamental period (one revolution along
  `theta_0`). Fig. 10 (p21) plots them against the continuation parameter: for each resonance, a central pair at zero and
  a pair of opposite non-zero exponents (read from plot, approximate: about +-0.52 for 9:2 synodic, about +-0.64 for 4:1 synodic,
  about +-0.57 for 4:1 sidereal, in the unit of the figure, which the caption does not state). Text (p20): "the Lyapunov
  exponents of their quasi-periodic counter parts remain approximately the same throughout the continuation procedure ...
  a bifurcation along one of the stable (zero) Lyapunov exponent pairs exists from the start of the perturbation, but does not
  depart significantly from the zero value (around 0.003 for e = 0.0549)" for the 9:2 case. Conclusion (p24): the stability
  "remains very close to the near-stable behavior presented by their periodic counterparts".
- Section 7, torus maps (pp21-23, Figs. 11 to 14): the conic eclipse model of Montenbruck & Gill (2000), maps of Sun
  visibility and of altitude over the Moon surface against the two torus angles. Colour-bar end labels as printed: 9:2 synodic
  altitude 853 to 65187 km (Fig. 12b); 4:1 synodic 4171 to 68697 km (Fig. 13b); 4:1 sidereal 1953 to 66327 km (Fig. 14b)
  (labels of the colour axis; whether they are the exact extrema of the torus is not stated). 9:2 synodic: "only two small eclipse
  regions, close to theta_0 = 3 pi/4 and theta_0 = 7 pi/4"; a path avoiding eclipse for "at least up to three full revolutions of
  theta_0 ... up to 3 x 2 T_syn, approximately 177 days" (p22); 4:1 synodic: no eclipse regions; 4:1 sidereal: "four sets of two
  eclipse regions, spaced approximately by pi/4 along the theta_0 direction and by pi along the theta_1 direction" (p23). Only
  Moon eclipses are detected (p22-23).

## 3. Cyclers, transfers and lunar flybys

READ by absence over all 25 pages: the paper treats NO cycler, no Earth-Moon transfer orbit and no resonant orbit with lunar
flybys of its own. Its subject is quasi-periodic tori that stay in the neighbourhood of the Moon (NRHO perilune altitudes of
the order of 1000 to several thousand kilometres, from the torus-map colour bars). The introduction mentions the work of
others on transfers: McCarthy & Howell 2022 (transfers from Earth vicinity to quasi-periodic Halo orbits of the BCR4BP via
invariant manifolds), Scheuerle & Howell 2021 (low-energy impulsive transfers towards the 9:2 synodic NRHO and a low lunar
orbit in the BCR4BP) and, in section 7, Henry & Scheeres 2020 (transfers between intersecting tori) (READ, pp2, 22). It
suggests phasing manoeuvres between "hyper-points" of a torus (p22). No lunar-flyby sequence, free-return or
Earth-Moon cycler appears.

## 4. What this means for the project

Project side (READ, module docstrings): `core/bcr4bp.py` is the incoherent bicircular model with the Sun at
`(a_S cos theta, a_S sin theta)`, `theta = theta_S0 - omega_S t` (clockwise, `#891`); `core/qbcp.py` is Andreu's QBCP with
the eight Fourier tables; `core/er3bp.py` is the ER3BP in the pulsating-rotating frame with true anomaly as the independent
variable (Szebehely 1967). The grep over `src` and `docs/notes` finds no ECR4BP module; only the `#884` literature check
and one digest mention the Elliptic-Circular model (INFERRED: the project has no ECR4BP).

1. The ECR4BP is, in the project's terms, the ER3BP of `core/er3bp.py` plus the Sun term of eq. 6 (INFERRED, from eqs. 5-8 and
   the er3bp module description, whose right-hand side was not re-derived here). Its limiting cases are checkable by the
   project's own models: with `e = 0` (so `nu` equals time in the unit of the primaries' rate) it must reproduce
   `core/bcr4bp.py` once `a3`, `n3` and the Sun's sense are matched; with `mu3 = 0` it must reproduce `core/er3bp.py`.
   These are identity tests, not published goldens.
2. Published numbers a project model can reproduce: the constants of Table 1 and `mu = 0.01215`, `e = 0.0549`,
   `mu3 = 328900.55`. The paper prints no initial condition, period or rotation number, so there is NO golden orbit in this paper
   to reproduce: a reproduction would be of the method (resonant multi-revolution periodic orbit in the bicircular model,
   continuation in `e` to a torus), checked qualitatively against Figs. 5 to 7, not digit by digit. The 9:2 and 4:1 synodic
   NRHOs in the BCR4BP are the same periodic orbits as in Boudad, Howell & Davis 2020 (digest
   `2026-10-04-digest-boudad-howell-davis-2020-synodic-resonant-nrho-bcr4bp.md`), which is the better source of printed
   bicircular numbers (INFERRED from the paper's own statement that its intermediate bicircular orbits are those of Boudad et
   al., p3, p12).
3. Sense of the Sun: eq. 9 (Sun at `(a3 cos sigma, -a3 sin sigma)`, `sigma` increasing) is a further printed statement, in the
   Colorado and Paris literature, of the clockwise sense that `#891` corrected the project to, and it agrees with Fig. 2 of
   Jorba-Cusco et al. 2018 as recorded in that digest (INFERRED). It is not a positive control with numbers.
4. Overlap with `#884` (Sun-perturbed periodic Earth-Moon orbits, continued in the Sun's strength; READ in the `#884` entry of
   `data/OUTSTANDING.md` and its literature check): the paper uses the same device, a homotopy in the Sun mass from 0 to the
   physical value (`epsilon` from 0 to 1, eq. 28 for the partial derivative), as the continuation in the Sun's strength recorded for `#884` (INFERRED from the Leiva-Briozzo summary in that entry),
   for resonant multi-revolution NRHOs in particular, and discusses folds and bifurcations that stop such continuations (3:1
   and 5:1 synodic resonances, p12 and p20). It does so for the NRHO family, not for the cycler-class periodic orbits of
   `#884`; and it adds the second perturbation (lunar eccentricity), which turns the periodic orbits into tori. Whether a
   periodic member of the `#884` family survives the lunar eccentricity as a torus is a question this paper's method
   (GMOS with `N_1` points on an invariant circle, `e` continuation) answers, but the paper does not do it for those orbits.
5. For the project's own resonant NRHO or cislunar work with `core/er3bp.py` or `core/bcr4bp.py`: the paper's list of
   continuations that fail (resonances other than 9:2, 4:1 synodic and 4:1 sidereal are reported as not converging, p12) is a
   warning that a failure of a continuation is a published occurrence for this family, to be checked against the bifurcation
   explanation before it is read as a negative result about the orbit.

## 5. References cited for earlier work on periodic and quasi-periodic orbits in Sun-perturbed Earth-Moon models

Transcribed as printed in the paper's reference list (pp24-25). Those cited in the text for the bicircular, quasi-bicircular
and quasi-periodic Sun-perturbed Earth-Moon work:

- Andreu, M. A. (1998). The quasi-bicircular problem. Ph.D. thesis.
- Boudad, K. K., Howell, K. C., & Davis, D. C. (2020). Dynamics of synodic resonant near rectilinear halo orbits in the
  bicircular four-body problem. Advances in Space Research, 66(9), 2194-2214. doi:10.1016/j.asr.2020.07.044.
- Castella, E., & Jorba, A. (2000). On the vertical families of two-dimensional tori near the triangular points of the
  bicircular problem. Celestial Mechanics and Dynamical Astronomy, 76(1), 35-54. doi:10.1023/A:1008321605028.
- Castella, E., & Jorba, A. (2003). The Lagrangian points in the real Earth-Moon system. In International Conference on
  Differential Equations (pp. 3-12). Hasselt, Belgium. (the paper's stated first mention of the elliptic-circular model, p5)
- Gomez, G., Simo, C., Llibre, J., & Martinez, R. (2001). Dynamics and mission design near libration points - Vol. II
  fundamentals: the case of the triangular libration points, volume 3. World Scientific.
- Gomez, G., & Mondelo, J. M. (2001). The dynamics around the collinear equilibrium points of the RTBP. Physica D: Nonlinear
  Phenomena, 157(4), 283-321. doi:10.1016/S0167-2789(01)00312-8.
- Jorba-Cusco, M., Farres, A., & Jorba, A. (2018). Two periodic models for the Earth-Moon system. Frontiers in Applied
  Mathematics and Statistics, 4(July), 1-14. doi:10.3389/fams.2018.00032.
- Rosales, J. J., Jorba, A., & Jorba-Cusco, M. (2021). Families of Halo-like invariant tori around L2 in the Earth-Moon
  Bicircular Problem. Celestial Mechanics and Dynamical Astronomy, 133(4), 1-30. doi:10.1007/s10569-021-10012-0.
- McCarthy, B., & Howell, K. C. (2022). Characterization of Families of Low-Energy Transfers to Cislunar Four-Body
  Quasi-Periodic Orbits. In AIAA SciTech Forum (pp. 1-17). San Diego, California. doi:10.2514/6.2022-1889.
- Scheuerle, S. T., & Howell, K. C. (2021). Characteristics and Analysis of Families of Low-Energy Ballistic Lunar Transfers.
  In AAS/AIAA Astrodynamics Specialist Conference (pp. 1-17). Big Sky, Montana (Virtual).
- Davis, D. C., Bhatt, S. A., Howell, K. C., Jang, J. W., Whitley, R. J., Clark, F. D., Guzzetti, D., Zimovan, E. M., & Barton,
  G. H. (2017). Orbit maintenance and navigation of human spacecraft at cislunar near rectilinear halo orbits. Advances in the
  Astronautical Sciences, 160, 2257-2276.
- Zimovan, E. M., Howell, K. C., & Davis, D. C. (2017). Near rectilinear halo orbits and their application in cis-lunar space.
  In 3rd IAA Conference on Dynamics and Control of Space Systems. Moscow, Russia.
- Howell, K. C., & Breakwell, J. V. (1984). Almost rectilinear halo orbits. Celestial Mechanics, 32(1), 29-52.
  doi:10.1007/BF01358402.
- NASA (2019). White Paper: Gateway Destination Orbit Model: A Continuous 15 Year NRHO Reference Trajectory.
- Henry, D. B., & Scheeres, D. J. (2020). Transfers between intersecting quasi-periodic tori. In AIAA SciTech Forum (pp. 1-17).
  Orlando, Florida.
- Scheeres, D. J. (1998). The Restricted Hill Four-Body Problem with Applications to the Earth-Moon-Sun System. Celestial
  Mechanics and Dynamical Astronomy, 70(2), 75-98. doi:10.1023/A:1026498608950.
- Scheeres, D. J., & Bellerose, J. (2005). The Restricted Hill Full 4-Body Problem: Application to spacecraft motion about
  binary asteroids. Dynamical Systems, 20(1), 23-44. doi:10.1080/14689360412331281267 (as printed; the doi digits are
  partly illegible in the scan).
- Assadian, N., & Pourtakdoust, S. H. (2010). On the quasi-equilibria of the BiElliptic four-body problem with non-coplanar
  motion of primaries. Acta Astronautica, 66(1-2), 45-58. doi:10.1016/j.actaastro.2009.05.014.
- Chakraborty, A., & Narayan, A. (2019). BiElliptic Restricted Four Body Problem. Few-Body Systems, 60(1).
  doi:10.1007/s00601-018-1472-x.

Method references (tori and continuation): Olikara, Z. P. (2016), Computation of quasi-periodic tori and heteroclinic
connections in astrodynamics using collocation techniques, Ph.D. thesis, University of Colorado Boulder; Baresi, N., Olikara,
Z. P., & Scheeres, D. J. (2018), Fully numerical methods for continuing families of quasi-periodic invariant tori in
astrodynamics, Journal of the Astronautical Sciences, 65(2), 157-182, doi:10.1007/s40295-017-0124-6; Villegas-Pinto, D.,
Baresi, N., Hestroffer, D., & Canalias, E. (2021), On the numerical computation of quasi-periodic families and applications
to the Martian Moons Exploration mission, in ICATT 2021; Seydel, R. (2009), Practical Bifurcation and Stability Analysis
(3rd ed.), Springer; Szebehely, V. (1967), Theory of orbits, Academic Press; Campagnola, S., Lo, M., & Newton, P. (2008),
Subregions of motion and elliptic halo orbits in the elliptic restricted three-body problem, Advances in the Astronautical
Sciences, 130 PART 2, 1541-1555.
