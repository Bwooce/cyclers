# Digest: Shu & Lin (2025), "Analysis of pitchfork bifurcations and symmetry breaking in the elliptic restricted three-body problem"

Celestial Mechanics and Dynamical Astronomy 137:21 (2025), DOI 10.1007/s10569-025-10252-4, open access (CC BY 4.0), 22 pages
(article text pp1-18, appendix pp19-20, references pp21-22). Haozhe Shu (Advanced Institute for Material Research and
Mathematical Institute, Tohoku University) and Mingpei Lin (Mathematical Institute, Tohoku University). Received 22 March
2025, accepted 3 June 2025, published online 26 June 2025.
Filed in the private paper corpus as
shu-lin-2025-pitchfork-bifurcations-symmetry-breaking-elliptic-restricted-three-body-problem-cmda-137-21-doi-10.1007-s10569-025-10252-4.pdf

Digested 2026-10-04 from all 22 pages. Each statement is marked READ (seen on the page, with page, section, equation or
figure) or INFERRED (our reading, or a comparison with project code). Page numbers are the journal's printed "Page n of 22".

## 0. What the paper is, and a warning about its content

READ (abstract, p1; section 1, p3): a semi-analytical (Lindstedt-Poincare, trigonometric series) treatment of the local
phase space near the collinear libration points L1, L2 (and in the formulation L3) of the spatial elliptic restricted
three-body problem (ERTBP). It introduces a "coupling coefficient" eta and "bifurcation equations" Delta(eta, e, alpha1,
alpha2, alpha3) = 0; a nonzero eta solving Delta = 0 means a pitchfork bifurcation, i.e. symmetry breaking, has occurred.

INFERRED, and important for how this paper is used: it is a local, formal-series paper. It prints no table of numbers. The
only printed numerical values of any kind are the Sun-Earth parameters (p12, Fig. 2 caption) and a handful of amplitudes
and eta values in figure captions (section 5 below). Nothing in it is a periodic-orbit initial condition, a period, a
stability index, or a bifurcation eccentricity. It does not treat resonance between the orbit period and the primaries'
period (no M:N commensurability appears anywhere), so it does not define "multi-revolution" orbits at all. The orbits it
constructs are formal truncated series (up to order 7 for the Fig. 3b examples, p14) about the libration point; no
differential correction, no continuation in e, no Floquet analysis and no integration is reported.

## 1. The model as printed (section 2, pp3-5)

READ (p3): "positioning the origin at the centroid of the two primaries, the orientation of the X-axis is given by the line
that goes from the smaller primary to the larger primary, while the Z-axis has the orientation determined by the angular
motion of the primaries. Y-axis completes the right-handed coordinate system. ... the normalized coordinates for the
smaller and the larger primary are (mu - 1, 0, 0) and (mu, 0, 0), respectively, where mu = m2/(m1 + m2) is the system
parameter". Pulsating-synodic frame, true anomaly f as independent variable, primes are d/df (p4).

Equations of motion, eq. (1), p4:

    X'' - 2Y' = dOmega/dX
    Y'' + 2X' = dOmega/dY
    Z'' + Z  = dOmega/dZ

    Omega(X,Y,Z,f) = 1/(1 + e cos f) * [ (1/2)(X^2 + Y^2 + Z^2) + (1-mu)/r1 + mu/r2 + (1/2) mu (1-mu) ]

    r1^2 = (X - mu)^2 + Y^2 + Z^2,   r2^2 = (X + 1 - mu)^2 + Y^2 + Z^2        (eq. 2)

READ: "When the eccentricity e = 0, this dynamical model reduces to the well-known autonomous CRTBP" (p4). Units: the
standard normalisation (distance between primaries, total mass, period 2 pi); not spelled out beyond "normalized" (p3).

Symmetries, eq. (3), p4 (READ; the sign pattern in the scan is small, and the S2 and S3 sign patterns are the ones expected
for a time reversal combined with a reflection, INFERRED):

    S1: (f, X, Y, Z, X', Y', Z') <-> (f, X, Y, -Z, X', Y', -Z')          reflection about the (X,Y) plane
    S2: (f, X, Y, Z, X', Y', Z') <-> (-f, X, -Y, Z, -X', Y', -Z')        time reversal with a Y reflection
    S3 = S1 composed with S2:   (f, ...) <-> (-f, X, -Y, -Z, -X', Y', Z')

Quote (p4): "Moreover, after time reversal, its reflection about the (X, Z) plane (X(-f), -Y(-f), Z(-f)) also satisfies the
governing differential equations." The paper says there are "five Lagrange points pulsating in the synodic coordinate
system" (p4).

INFERRED comparison with project code (`core/er3bp.py`, pulsating-rotating frame, true anomaly f, primes d/df): the
equations are the same system. Differences are conventions only. (i) Paper: larger primary at x = mu, smaller at x = mu - 1
(the Barcelona-school and Jorba convention, mirror of the Szebehely placement the project uses: larger at -mu, smaller at
1 - mu). (ii) Paper writes the z-equation as Z'' + Z = dOmega/dZ with the full (1/2)Z^2 inside Omega; the project writes
the equivalent z'' = -(e cos f z + grav_z)/(1 + e cos f). These agree: Z'' = -Z + Z/(1+e cos f) + grav term/(1+e cos f)
= -e cos f Z/(1+e cos f) + grav term/(1+e cos f). Only the mirror of the X axis (and the sign of Y accordingly) needs care
if a figure's initial conditions were ever compared.

### Libration-point-centred coordinates (eqs. 4-10, pp4-5)

READ, eq. (4): X = -gamma_i x + mu + a_i, Y = -gamma_i y, Z = gamma_i z for i = 1, 2 (with a_1 = -1 + gamma_1,
a_2 = -1 - gamma_2); X = gamma_i x + mu + gamma_i, Y = gamma_i y, Z = gamma_i z for i = 3. "gamma_i denotes the instantaneous
distance between the libration point L_i and its closest primary", the unique positive root of Euler's quintic, eq. (5),
printed as:

    gamma^5 -/+ (3-mu) gamma^4 + (3-2mu) gamma^3 - mu gamma^2 +/- 2 mu gamma - mu = 0,   i = 1, 2   (upper sign L1, lower L2)
    gamma^5 + (2+mu) gamma^4 + (1+2mu) gamma^3 - (1-mu) gamma^2 - 2(1-mu) gamma - (1-mu) = 0,   i = 3

The motion is then eq. (6), x'' - 2y' = (1/gamma_i^2) dOmega/dx etc., expanded into the Richardson-type recurrence, eq. (7),
with homogeneous polynomials T_n, R_n from eqs. (8) and (9) (T_n = ((2n-1)/n) x T_{n-1} - ((n-1)/n)(x^2+y^2+z^2) T_{n-2};
R_n = ((2n+3)/(n+2)) x R_{n-1} - ((2n+2)/(n+2)) T_n - ((n+1)/(n+2))(x^2+y^2+z^2) R_{n-2}; T_0 = 1, T_1 = x, R_0 = -1,
R_1 = -3x) and coefficients c_n(mu) from eq. (10). Eq. (7) carries the e-dependence as sums of (-e)^i cos^i f. The linear
part, eq. (11), contains the e-terms; the autonomous linear part, eq. (12), gives the frequencies, eq. (13):

    omega_0 = sqrt( (2 - c2 + sqrt(9 c2^2 - 8 c2)) / 2 ),  nu_0 = sqrt(c2),  lambda_0 = sqrt( (c2 - 2 + sqrt(9 c2^2 - 8 c2)) / 2 )
    kappa_1 = -(omega_0^2 + 2 c2 + 1)/(2 omega_0),   kappa_2 = -(lambda_0^2 - 2 c2 - 1)/(2 lambda_0)

with x = alpha1 cos(theta1) + alpha3 cos(theta3), y = kappa1 alpha1 sin(theta1) + sqrt(-1) kappa2 alpha3 sin(theta3),
z = alpha2 cos(theta2), theta1 = omega0 f + phi1, theta2 = nu0 f + phi2, theta3 = sqrt(-1) lambda0 f + phi3 (eq. 13, p6).

## 2. What is computed: the "unified framework" and the bifurcation equations

READ, abstract (p1): "a unified trigonometric series-based framework is proposed to analyze these bifurcated families from
the perspective of coupling-induced bifurcation mechanisms. By introducing a coupling coefficient and various bifurcation
equations into the ERTBP, different symmetry breaking is achieved when the coupling coefficient is nonzero. This unified
semi-analytical framework captures bifurcations of both periodic/quasi-periodic and transit/non-transit orbits.
Furthermore, it reveals that pitchfork bifurcation solutions in the ERTBP fundamentally depend solely on the orbital
eccentricity and three amplitude parameters of the system's degrees of freedom, governing both the elliptic direction and
the hyperbolic one."

What "unified" means here (READ, sections 1 and 3, pp3, 6-9): (a) the hyperbolic part of the solution is written as a
trigonometric function of complex phase with a complex amplitude alpha3, instead of the exponential form of Masdemont (2005)
and Lei et al. (2013); alpha3 real gives non-transit orbits, alpha3 in sqrt(-1) R gives transit orbits (p6: "The solution
with amplitude alpha3 lying on the real axis describes the motion of non-transit orbit, while the imaginary-valued alpha3
corresponds to transit motion"); (b) the in-plane and out-of-plane frequencies are not forced into 1:1 resonance as in
the classical halo construction (Richardson 1980; Lei et al. 2013), the symmetry is broken instead by adding a coupling term
eta*Delta to one equation (Remark 2, p9: "frequencies of both in-plane and out-of-plane motion are preserved in the
coupling-bifurcation computation").

The quantity that characterises the pitchfork (READ): the coupling coefficient eta. Delta = 0 is the "bifurcation equation"
(p8: "For any choice of the quartet (alpha1, alpha2, alpha3, e), if there exist some eta != 0 satisfying the polynomial
bifurcation equation Delta = 0, it indicates the occurrence of a bifurcation. Conversely, no bifurcation occurs if eta = 0.").
Delta is expanded as Delta = sum d_ijkm alpha1^i alpha2^j alpha3^k e^m (p8), frequencies as power series in the same
variables, eq. (17). Reference: "pitchfork bifurcation solutions in the ERTBP fundamentally depend solely on the orbital
eccentricity and three amplitude parameters" (abstract, p1), and "all types of bifurcated solutions in the form of
trigonometric series are shown to depend solely on the eccentricity and three amplitudes corresponding to the system's DOFs,
where the critical conditions are derived explicitly" (p3).

### Coupling constructions (Table 1, p10; section 3.1-3.2, pp7-9)

Table 1 (p10), transcribed:

| Type of symmetry breaking | Coupling direction | Type of bifurcated orbit |
| --- | --- | --- |
| S1 | x -> z | Halo/quasi-halo orbits and their corresponding transit/non-transit orbits |
| S2 | y -> z, z -> y | Axial/quasi-axial orbits and their corresponding transit/non-transit orbits |

READ: S1 breaking (planar Lyapunov to halo, North-South symmetry broken): the z equation of (7) gets + eta*Delta*x, eq. (7.1);
modified linear solution eq. (15) with kappa_3 = (nu0^2 - omega0^2)/(nu0^2 + lambda0^2), d_0000 = nu0^2 - omega0^2.
S2 breaking (vertical Lyapunov to axial): two variants, + eta*Delta*y in the z equation, eq. (7.2), solution eq. (18), with
d_0000 = (nu0^2 - omega0^2)/kappa_1, kappa_3 = (kappa_2/kappa_1)(nu0^2 - omega0^2)/(nu0^2 + lambda0^2); and + eta*Delta*z in
the y equation, eq. (7.3), solution eq. (19), with d_0000 = 1/(2 nu0) - nu0/2, kappa_3 = -1/(2 nu0) - 3 nu0/2. The formal
solutions are eq. (16) (cosine series, S1) and eq. (20) (sine series, S2). Remark 3 (p9): S3 breaking is "a combination of
the first two, meaning that its absence essentially corresponds to two sequential symmetry breaking"; it corresponds to the
W4 and W5 families of the global bifurcation diagram of Doedel et al. (2007), which "is significantly far from the collinear
libration points" so the local Lindstedt-Poincare method does not reach it.

## 3. Bifurcation equations and critical conditions (section 4, pp10-18)

All "READ". Third-order truncation. Printed forms:

- S1 (halo), eq. (21), p10: Delta = d_0000 + d_2000 alpha1^2 + d_0200 alpha2^2 + d_0020 alpha3^2 + d_0002 e^2 = a eta^4 + b eta^2 + c = 0,
  with a = l1 alpha1^2 + l2 alpha3^2, b = l3 alpha1^2 + l4 alpha2^2 + l5 alpha3^2, c = l6 alpha1^2 + l7 alpha2^2 + l8 alpha3^2 + l9 e^2 +
  (nu0^2 - omega0^2). "The coefficients l_i (i = 1..9) depend solely on the system parameter mu" (Fig. 1, p11: plotted against mu in
  (0, 0.5) for L1; values only as curves, not tabulated). A quadratic in eta^2.
  Quote (p10): "Here, the relatively negligible impact of the small orbital eccentricity on the bifurcation equation
  illustrates some significant similarities shared by the non-autonomous ERTBP and its approximated circular model from
  quantitative perspectives."
  Sign facts read from text and Fig. 1: l6 > 0 and l7, l8 < 0 for all mu in (0, 0.5) (p11), so c = 0 is a one-sheet hyperboloid in
  (alpha1, alpha2, alpha3) for the transit case (eq. 22). Cases 1.1 to 1.3 (pp11-12) give critical surfaces c = 0 (eq. 22),
  a = 0 (eq. 24, a = l1 alpha1^2 - l2 alpha3^2 for transit) and the discriminant surface b^2 - 4ac = 0 (eq. 25); feasible eta
  from eq. (23), eta = +/- sqrt((-b +/- sqrt(b^2 - 4ac))/(2a)); four, two or no feasible solutions depending on the region (Fig. 2a, p13).
- Center-manifold reduction (alpha3 = 0), eqs. (26)-(28), pp12-13: the bifurcation curve in the (alpha1, alpha2) plane is the
  hyperbola l6 alpha1^2 + l7 alpha2^2 = omega0^2 - nu0^2 - l9 e^2 (eq. 27); the critical amplitude is
  |alpha1,cri| = sqrt((omega0^2 - nu0^2 - l9 e^2)/l6) (Fig. 3a). For alpha2 = 0 and |alpha1| > |alpha1,cri|, the pair
  eta = +/- sqrt( (-l3 alpha1^2 - sqrt(l3^2 alpha1^4 - 4 l1 alpha1^2 (l6 alpha1^2 + l9 e^2 + (nu0^2 - omega0^2))))/(2 l1 alpha1^2) ),
  eq. (28), is "the northern halo orbits and southern halo orbits, respectively"; alpha2 != 0, eta != 0 gives quasi-halo orbits;
  eta = 0 gives Lyapunov / Lissajous orbits (p13).
- Hyperbolic (alpha1 = alpha2 = 0), eqs. (29)-(30), p14: Delta = l2 alpha3^2 eta^4 + l5 alpha3^2 eta^2 + l8 alpha3^2 + l9 e^2 +
  (nu0^2 - omega0^2) = 0; bifurcated hyperbolic orbits leave the (x,y) plane (Fig. 4b). Transit/non-transit (section 4.1.4, pp14-15):
  transit orbits bifurcate too; for large |eta| they are dominated by the hyperbolic direction and escape the plane quickly,
  for small |eta| they are "transit orbits of halo orbits" (Fig. 5).
- S2, z -> y coupling, eqs. (31)-(34), pp15-16: Delta = h1 alpha2^2 eta^2 + h2 alpha1^2 + h3 alpha2^2 + h4 alpha3^2 + h5 e^2 + 1/(2 nu0) - nu0/2 = 0;
  "In the Sun-Earth system, the coefficients satisfy h2, h3 > 0, h4 < 0" (p15). Critical surface eq. (32); for alpha3 = 0 the
  critical curve is an ellipse in (alpha1, alpha2) (Fig. 7a); inside it eta = +/- sqrt(-(1/(h1 alpha2^2))(h2 alpha1^2 + h3 alpha2^2 + h5 e^2 +
  1/(2 nu0) - nu0/2)), eq. (33), two families of axial orbits bifurcating from vertical Lyapunov orbits, with
  |alpha2| <= sqrt(-(1/h3)(h5 e^2 + 1/(2 nu0) - nu0/2)) (eq. 34, alpha1 = 0).
- S2, y -> z coupling, eqs. (35)-(38), pp17-18, eq. (35) as printed: Delta = (k1 alpha1^2 + k2 alpha3^2) eta^4 +
  (k3 alpha1^2 + k4 alpha2^2 + k5 alpha3^2) eta^2 + (k6 alpha1^2 + k7 alpha2^2 + k8 alpha3^2 + k9 e^2 + (nu0^2 - omega0^2)/kappa_1) = 0.
  Reduced (alpha3 = 0), eq. (36): k1 alpha1^2 eta^4 + (k3 alpha1^2 + k4 alpha2^2) eta^2 + (k6 alpha1^2 + k7 alpha2^2 + k9 e^2 +
  (nu0^2 - omega0^2)/kappa_1) = 0. Critical |alpha1,cri| = sqrt(-(1/k6)(k9 e^2 + (nu0^2 - omega0^2)/kappa_1)), eq. (37); eta from
  eq. (38). Fig. 9: axial orbits with alpha1 = 0.28 and alpha1 = 0.2.
  (The k_i, like l_i and h_i, depend only on mu; no values are printed.)

Role of e (READ and INFERRED): e enters the third-order bifurcation equations only through one extra term (l9 e^2, h5 e^2,
k9 e^2), additive to the circular constant (nu0^2 - omega0^2) etc. READ conclusion (p18): "the emergence of bifurcations is
governed by the existence of feasible solutions (eta != 0) to the bifurcation equation Delta = 0, where solutions are
determined solely by the orbital eccentricity e and three amplitude parameters alpha_i (i = 1, 2, 3)." INFERRED: at this
truncation order e only shifts the critical amplitude (for example |alpha1,cri|^2 = (omega0^2 - nu0^2 - l9 e^2)/l6); it does
not create new bifurcation types and says nothing about the period-commensurability needed for a genuinely periodic
orbit at e > 0.

## 4. Comparison with the circular problem

READ: the formulas reduce to the CRTBP for e = 0 (p4 and p10, quote above). The S1 and S2 pitchfork structure (planar
Lyapunov to halo, vertical Lyapunov to axial) is the classical CRTBP pitchfork, cited to the global bifurcation diagram of
Doedel et al. 2007 (p9, Remark 3). The new content claimed is (i) a uniform treatment extending to transit/non-transit
hyperbolic motion, (ii) the explicit e-dependence at third order. No numerical comparison of ERTBP against CRTBP bifurcation
values (for example a critical amplitude at e = 0 versus e = 0.0167) is printed.

## 5. Every printed number (no tables exist)

| Where | Quantity | Printed value |
| --- | --- | --- |
| p12, Fig. 2 text | Sun-Earth system parameter mu | 3.040423398444176e-6 |
| p12, Fig. 2 text | Sun-Earth orbital eccentricity e | 0.01671022 |
| p13, Fig. 3b caption | halo orbit (red) | alpha1 = 0.15 |
| p13, Fig. 3b caption | quasi-halo orbit (green) | alpha1 = 0.15, alpha2 = 0.04 |
| p15, Fig. 5 caption | transit orbits, branch (a) | alpha1 = 0.15, alpha3 = +/- 0.005 sqrt(-1), eta = 0.7085 |
| p15, Fig. 5 caption | transit orbits, branch (b) | alpha1 = 0.15, alpha3 = +/- 0.005 sqrt(-1), eta = 9.4484 |
| p16, Fig. 7b caption | axial orbit (red) | alpha1 = 0, alpha2 = 1 |
| p16, Fig. 7b caption | quasi-axial orbit (green) | alpha1 = 0.005, alpha2 = 1 |
| p18, Fig. 9b caption | axial orbits | alpha1 = 0.28 (red), alpha1 = 0.2 (blue) |

Fig. 1 (p11) plots l1 to l9 and d_0000 against mu in (0, 0.5) for L1 as curves only. Only signs are usable (l6 > 0, l7 and
l8 < 0, stated in the text); axis magnitudes in the scan are small and are not transcribed. The amplitude values above are
in the paper's own normalised libration-point-centred units. Orbit periods, Floquet multipliers and Jacobi-like quantities
are never printed.

## 6. Method in a few lines

READ: Richardson-style expansion of the right-hand side about the collinear point (eqs. 7-10), a Lindstedt-Poincare
perturbation in the three amplitudes (alpha1, alpha2, alpha3) and e (eqs. 16, 17, 20), with the unknown series coefficients
found order by order by solving the linear systems of the appendix, eqs. (39)-(49) (the resonant singular cases (s,t,u,r) =
(1,0,0,0), (0,1,0,0), (0,0,1,0) are handled by fixing x_1000 = 0 etc., eqs. 43, 48). Frequencies omega, nu, lambda are
amplitude-and-e dependent series. Once the series coefficients are known, the third-order Delta fixes the critical surfaces.
Fig. 3b and others are drawn from the series to order 7 (p14). No numerical validation against differential correction or
against an independent integration is shown; "Data Availability: The datasets generated during and analyzed during the
study are available from the corresponding author upon reasonable request" (p21).

## 7. What this means for the project

INFERRED throughout.

Positive-control candidates: none with printed orbit numbers. The only numeric inputs are Sun-Earth mu = 3.040423398444176e-6
and e = 0.01671022, with no printed output to compare against. Two weak, structure-only checks the project's
`core/er3bp.py` machinery could support, both qualitative: (a) the sign of the planar-to-halo pitchfork exists at e = 0.01671022
for the Sun-Earth L1 Lyapunov family (halo and quasi-halo families emerge beyond a critical in-plane amplitude); (b) the critical
amplitude shifts by a term proportional to e^2 relative to the CR3BP. Neither has a published number to match, so neither is
a positive control in the project's sense (expected value traceable to the source); they would be consistency checks only.
A real control would need the series coefficients l_i, h_i, k_i, which are plotted (Fig. 1) but not tabulated, or the
authors' data by request.

Overlap with open tasks: #435 and #437 (e > 0 continuation and fold-aware continuation of planar and multi-revolution
families). This paper does not contribute to them. It is a local, formal, near-libration-point analysis; it contains no
continuation in e, no commensurability condition, no folds, no isolated branches, and says nothing about the
multi-revolution resonant halo orbits that Peng and Xu 2015 and Peng, Bai and Xu 2017 treat. It does state, in its
introduction (p2), that "Peng and Xu (2015) generated multi-revolution halo orbits via continuation methods and
multi-segments optimization methods" and "Ferrari and Lavagna (2018) succeeded in finding periodic orbits in the ERTBP
through a differential correction algorithm", citations only; no result, table or initial condition from Peng and Xu 2015
is reproduced. It therefore does NOT fill any part of the Peng and Xu 2015 gap.
A modest conceptual relevance: eq. (21) and its siblings say that, near a collinear point, e only perturbs the circular
bifurcation condition at order e^2, so a continuation in e of a halo orbit that bifurcates from a Lyapunov orbit should
carry the pitchfork with it with a small shift; the multi-revolution resonant families of Peng et al. are not of this
local kind and are not covered.

Catalogue classes: nothing here is a cycler, quasi-cycler, or resonant periodic orbit with flybys. It is
libration-point-orbit theory (Lyapunov, halo, axial, Lissajous, transit orbits) with no flybys and no resonance with the
primaries' period. Out of catalogue scope; of background value only.

## 8. References cited by the paper on multi-revolution halo orbits and elliptic-problem continuation (full citations, pp21-22)

Printed reference-list entries relevant to the subject (READ):
- Peng, H., Xu, S.: Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted three-body problem. Celest. Mech. Dyn. Astron. 123(3), 279-303 (2015), https://doi.org/10.1007/s10569-015-9635-2 (cited at p2 for generating multi-revolution halo orbits "via continuation methods and multi-segments optimization methods"; this is the project's corrupt-copy gap paper).
- Peng, H., Bai, X.L., Masdemont, J.J., Gomez, G., Xu, S.J.: Libration transfer design using patched elliptic three-body models and graphics processing units. J. Guid. Control Dyn. 40(12), 3155-3166 (2017), https://doi.org/10.2514/1.G002692
- Ferrari, F., Lavagna, M.: Periodic motion around libration points in the elliptic restricted three-body problem. Nonlinear Dyn. 93(1), 453-462 (2018), https://doi.org/10.1007/s11071-018-4203-4
- Paez, R.I., Guzzo, M.: Transits close to the Lagrangian solutions L1, L2 in the elliptic restricted three-body problem. Nonlinearity 34(9), 6417-6449 (2021), https://doi.org/10.1088/1361-6544/ac13be
- Paez, R.I., Guzzo, M.: On the semi-analytical construction of halo orbits and halo tubes in the elliptic restricted three-body problem. Phys. D Nonlinear Phenom. 439, 133402 (2022), https://doi.org/10.1016/j.physd.2022.133402
- Lei, H.L., Xu, B., Hou, X.Y., Sun, Y.S.: High-order solutions of invariant manifolds associated with libration point orbits in the elliptic restricted three-body problem. Celest. Mech. Dyn. Astron. 117(4), 349-384 (2013), https://doi.org/10.1007/s10569-013-9515-6
- Jorba, A., Nicolas, B., Rodriguez, O.: A dynamical study of Hilda asteroids in the circular and elliptic RTBP. Chaos Interdiscip. J. Nonlinear Sci. (2024), https://doi.org/10.1063/5.0234410 (volume as printed; not transcribed)
- Jorba-Cusco, A., Epenoy, R.: Low-fuel transfers from Mars to quasi-satellite orbits around Phobos exploiting manifolds of tori. Celest. Mech. Dyn. Astron. 133(5), (2021), https://doi.org/10.1007/s10569-021-10017-9
- Celletti, A., Lhotka, C., Pucacco, G.: The dynamics around the collinear equilibrium points of the elliptic three-body problem: a normal form approach. Phys. D Nonlinear Phenom. 468, 134302 (2024), https://doi.org/10.1016/j.physd.2024.134302
- Lin, M., Luo, T., Chiba, H.: Semi-analytical computation of bifurcation of orbits near collinear libration points in the elliptic restricted three-body problem. Phys. D Nonlinear Phenom. 470, 134404 (2024), https://doi.org/10.1016/j.physd.2024.134404
- Lin, M., Chiba, H.: Bifurcation mechanism of Quasi-Halo orbit from Lissajous orbit. J. Guid. Control Dyn. 48(1), 71-83 (2025), https://doi.org/10.2514/1.G008233
- Masdemont, J.J.: High-order expansions of invariant manifolds of libration point orbits with applications to mission design. Dyn. Syst. 20(1), 59-113 (2005), https://doi.org/10.1080/14689360412331304291
- Jorba, A., Masdemont, J.M.: Dynamics in the center manifold of the collinear points of the restricted three body problem. Phys. D Nonlinear Phenom. 132(1-2), 189-213 (1999), https://doi.org/10.1016/S0167-2789(99)00042-1
- Doedel, E.J., Romanov, V.A., Paffenroth, R.C., Keller, H.B., Dichmann, D.J., Galan-Vioque, J., Vanderbauwhede, A.: Elemental periodic orbits associated with the libration points in the circular restricted 3-body problem. Int. J. Bifurc. Chaos 17(08), 2625-2677 (2007), https://doi.org/10.1142/S0218127407018671
- Broucke, R.: Stability of periodic orbits in the elliptic, restricted three-body problem. AIAA J. 7(6), 1003-1009 (1969), https://doi.org/10.2514/3.5267
- Richardson, D.L.: Analytic construction of periodic orbits about the collinear points. Celest. Mech. 20(3), 241-253 (1980), https://doi.org/10.1007/BF01229511 (volume as printed)
- Szebehely, V.: Chapter 10 - Modifications of the restricted problem, Theory of Orbits, Academic Press, 556-652 (1967)
- Antoniadou, K.I., Voyatzis, G.: 2/1 resonant periodic orbits in three dimensional planetary systems. Celest. Mech. Dyn. Astron. 115, 161-184 (2013), https://doi.org/10.1007/s10569-012-9457-4

Of these, the project already holds Antoniadou-Voyatzis-class and Peng-Bai-Xu 2017 material (see CORPUS_INDEX); the items
not yet obviously in the corpus and plausibly worth acquiring are Paez and Guzzo 2021 and 2022 (Floquet-Birkhoff
normalisation of halo orbits and transits in the ERTBP), Ferrari and Lavagna 2018, Lei et al. 2013 and Celletti et al. 2024.
(Corpus status not re-checked here; check CORPUS_INDEX before filing.)

## 9. Slips noted, stated factually

- The paper switches between "alpha3 in sqrt(-1) R" for transit and "alpha3 in R" for non-transit; the first sentence of
  Case 1.1 (p11) says "alpha3 in sqrt(-1) R" for the transit case while the surrounding text describes alpha3 in R as the
  non-transit case and then restates both. The cases are internally consistent once read in order; flagged only because a
  quick read can swap them.
