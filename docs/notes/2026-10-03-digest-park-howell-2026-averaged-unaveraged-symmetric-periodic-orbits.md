# Digest — Park & Howell (2026), "Linking Averaged and Unaveraged Three-Body Dynamics Near Smaller Primaries: Symmetric Periodic Orbits" (arXiv:2606.08485v1)

**Digested:** 2026-10-03. Text-layer PDF, 38 pages, read in full via `pdftotext -layout`; Table 5 (page 17)
checked against the rendered page and the text layer matches. Table 4 and Table 6 use coloured
triangle markers that the text layer drops (only the apse identifiers are transcribed below).
Filed in the private paper corpus as
`park-howell-2026-linking-averaged-unaveraged-three-body-dynamics-smaller-primaries-symmetric-periodic-orbits-arxiv-2606.08485v1.pdf`
(md5 `f7b06dfa55acb502c48062168b514a94`).

## Citation as printed
Beom Park (Apollo 11 Postdoctoral Fellow) and Kathleen C. Howell (Hsu Lo Distinguished Professor), School of
Aeronautics and Astronautics, Purdue University, West Lafayette, IN, USA, 47907 (park1103@purdue.edu,
howell@purdue.edu). "Linking Averaged and Unaveraged Three-Body Dynamics Near Smaller Primaries: Symmetric
Periodic Orbits", arXiv:2606.08485v1 [math.DS], 7 Jun 2026. Preprint; no journal is named in the text.

## What the paper is
A theory-plus-demonstration paper on orbits close to the smaller primary (an "orbiter" regime). It maps the
equilibria of the integrable doubly averaged model (the Lidov-Kozai-type model) to symmetric periodic orbits
(POs) of the unaveraged Hill and circular restricted three-body problems, using a common frequency lattice,
the parity of the resonance integers, and symmetry counting to predict a priori how many branches and which
symmetry types exist. It then traces how those families connect (to planar families, halo orbits, collision
orbits) and draws archetypical bifurcation diagrams. Abstract: "explicitly linking averaged equilibria to
symmetric periodic orbits in the unaveraged three-body systems. A unified frequency framework is introduced to
characterize the mapping of invariant tori across the dynamical models. Leveraging the parity of the
resonance ratio, an initialization scheme is developed to identify admissible apse configurations, enabling
the a priori prediction of solution multiplicity and symmetry types."

## Models, systems, constants
- Three models (Section 2.2): the Doubly Averaged Dynamical Model (DADM; "P1 is the only perturbing force",
  circular primaries, Legendre-truncated perturbation averaged over P1's period and the spacecraft's Keplerian
  period); the Hill restricted three-body problem (HR3BP; second-order Legendre term, "limiting assumption that
  mu -> 0", equation v_H' = -2 z x v_H + grad U_H, U_H = (3 x_H^2 - z_H^2)/2 + 1/r_H); and the circular
  restricted three-body problem (CR3BP, equations of motion in the barycentric rotating frame, pseudo-potential
  U = (x^2 + y^2)/2 + (1-mu)/d + mu/r). All three are spatial (the planar orbits are the i = 0 and 180 degree limits).
- Notation: mu = mu2/(mu1 + mu2); l* constant primary distance; t* = sqrt(l*^3/(mu1 + mu2)); Hill frame
  scaled by mu^(1/3); a_H = a/(l* mu^(1/3)); frames BRF, HRF, PIF.
- Systems and constants: the only numerical mass ratio printed is "mu ~ 0.01215 that represents the
  Earth-Moon system" (used for all CR3BP results and Fig. 3(c),(d); Fig. 3(a),(b) use mu = 0 for the HR3BP
  limit). The paper says "While the model is applicable to general three-body systems, this work specifically
  adopts mu ~ 0.01215". No other planet-moon system's mass ratio or constants are printed. Planet-moon
  systems named in the text: Earth's Moon; "Icy satellites of Jupiter and Saturn" (introduction, motivation);
  Europa (as a cited science-orbit application, ref [15]). Jupiter-Europa and Saturn-Enceladus appear only as
  the title of Moreno et al. (ref [23]); Ganymede only in the JUICE reference title (ref [2]).
- The HR3BP orbit plots use "the radius of the secondary body ... arbitrarily set to 0.05 [nd]".
- DADM range limit: a_H <= 0.45 nd, "corresponding to a short-period frequency nu_s >~ 3.3 [radians/nd]".

## Objects computed and methods
- DADM (Section 3): integrals C_Z = (1-e^2) cos^2 i; C_H = e^2 (2/5 - sin^2 i sin^2 omega); C_a = a (Eqs. 12-14).
  Circular equilibria (e = 0): elliptic for sin^2 i < 2/5, hyperbolic for sin^2 i > 2/5, bifurcation at
  sin^2 i = 2/5. Frozen equilibria (e > 0): k = 0, h = +/- sqrt((5 sin^2 i - 2)/3), i.e. e = sqrt((5 sin^2 i - 2)/3),
  omega = +/- 90 degrees; they exist for sin^2 i > 2/5, i_crit- ~ 39.2 degrees and i_crit+ ~ 140.8 degrees.
  Lidov diagram (Fig. 2).
- Frequency framework: three phase angles phi_s = M, phi_m = t + t0 - Omega, phi_l with frequencies nu_s, nu_m,
  nu_l; dimensional lift T_Dfull from T_Dsep with D_full = D_sep + 2 - r (Eq. 29); resonance integer relation
  m . nu = 0. Focus: T1 (PO) lifted from T0 (equilibria) with r = 1, p nu_m = q nu_s... as printed "p nu_m = q nu_s
  with coprime p, q"; resonance ratio eta := nu_s/nu_m = p/q (Eq. 31), p > q: "the spacecraft completes p cycles in
  phi_s (mean anomaly) while completing q cycles in phi_m (negative of the RAAN within the rotating frames)".
  Period P = p (2 pi/nu_s). nu_s = n/n_1 depends on a only. nu_m circular = 1 + (3/4)(1-mu) cos i (n_1/n) (Eq. 32);
  nu_m frozen = 1 + (1/4)(1-mu)(n_1/n)(3/5)(20 sin^2 i - 5) cos i/|cos i| (Eq. 33), with n_1/n = a_H^3 in Hill scaling
  (the printed expression for the circular formula has the symbol nr in the extracted text; read as n).
- Time-reversing symmetries (Table 2; the HR3BP group has order 8, the CR3BP group order 4): fixed sets
  Sigma_OX: y_H = z_H = xdot_H = 0 (HR3BP and CR3BP); Sigma_XOZ: y_H = xdot_H = zdot_H = 0 (both); Sigma_OY:
  x_H = z_H = ydot_H = 0 (HR3BP only); Sigma_YOZ: x_H = ydot_H = zdot_H = 0 (HR3BP only). Singly symmetric
  (N = 1), doubly symmetric (N = 2); "in the generic spatial domain, triply symmetric orbits do not exist";
  planar orbits are at least singly symmetric (sigma), and may be triply symmetric.
- Apse configurations: 16 configurations (Table 3), identifiers vertaxis with vert in {A, N, D, S} (ascending
  node, north apex, descending node, south apex) and axis in {+x, -x, +y, -y}. Modulo-4 arithmetic of
  (p, q, +/-p - q) gives apsidal sequences per quarter period (Table 4).
- Initialisation (Section 4.3, Fig. 7): fix eta = p/q and either a or i (Eqs. 32-33 give the other), choose
  equilibrium type (e, omega), symmetry pairing and initial apse via Table 6 (u = omega + theta in
  {0, 90, 180, 270} degrees; Omega = chi_H in {0, 90, 180, 270} degrees; t0 = +/- u), transform osculating to
  rotating-frame state, then "differential correctors ... target the symmetric apse conditions at t = P/2 (for
  singly symmetric orbits) or t = P/4 (for doubly symmetric orbits), following standard procedures, e.g., by
  Robin and Markellos [39]"; long-period orbits may be discretised into multiple segments. The osculating
  elements are taken equal to the doubly averaged elements (footnote 6).
- Families and bifurcations (Section 5): continuation by one-parameter families; monodromy-matrix stability
  indices s = (lambda_z + 1/lambda_z)/2 (Eq. 43), s_j for j = 1, 2 (Eq. 50); DADM prediction s_DADM = cos(2 pi/eta_0)
  with eta_0deg = eta - 1 and eta_180deg = eta + 1 (Eqs. 42-47); period-multiplying, pitchfork, "broken"
  pitchfork, fold (saddle-node) and period-doubling bifurcations; Kustaanheimo-Stiefel regularisation for the
  rectilinear collision orbits "adapted from Howell and Breakwell [19]". Fig. 14-17 plots; Fig. 18 archetypical
  bifurcation diagrams.

## Exact tables (transcribed; tables of orbital data do not exist, see below)

Table 1 (Dimensional lift, rows = reduced DADM tori, columns = full-phase-space resonance count r):

| DADM tori (separated) | r = 0 | r = 1 | r = 2 |
|---|---|---|---|
| T0 (Equilibria) | T2 (2D QPO) | T1 (PO) | N/A |
| T1 (PO) | T3 (3D QPO) | T2 (2D QPO) | T1 (PO) |

Table 2 (Time-reversing symmetries and fixed sets):

| Symmetry | Transformation (Sx, Sy, Sz) | Valid models | Fixed set |
|---|---|---|---|
| OX | x-axis rotation (+1, -1, -1) | HR3BP and CR3BP | y_H = z_H = xdot_H = 0 |
| XOZ | xz-plane reflection (+1, -1, +1) | HR3BP and CR3BP | y_H = xdot_H = zdot_H = 0 |
| OY | y-axis rotation (-1, +1, -1) | HR3BP only | x_H = z_H = ydot_H = 0 |
| YOZ | yz-plane reflection (-1, +1, +1) | HR3BP only | x_H = ydot_H = zdot_H = 0 |

Table 3 (Symmetric apse configurations, identifiers; prograde assumed):
- Axial (Z = 0), OX apses: A+x, D+x, A-x, D-x (ascending/descending at 0 or 180 degrees); OY apses: A+y, D+y, A-y, D-y
  (at 90 or 270 degrees).
- Reflectional (Zdot = 0), XOZ apses: N+x, S+x, N-x, S-x; YOZ apses: N+y, S+y, N-y, S-y.
- Apse location chi_H = 0, 90, 180, 270 degrees for +x, +y, -x, -y.

Table 4 (Circular equilibria; evolution of symmetric apsides per quarter period; "Pro / Ret" azimuthal steps; sample
prograde sequence at t = 0, P/4, 2P/4, 3P/4 as printed in the text layer):

| Parity (p:q) | p mod 4 | q mod 4 | +/-p - q mod 4 | Latitudinal step | Azimuthal step (Pro / Ret) | Sample sequence (hz > 0) |
|---|---|---|---|---|---|---|
| odd:odd | 1 | 1 | 0/2 | +90 | 0 / 180 | A+x, N+x, D+x, S+x |
| odd:odd | 1 | 3 | 2/0 | +90 | 180 / 0 | A+x, N-x, D+x, S-x |
| odd:odd | 3 | 1 | 2/0 | -90 | 180 / 0 | A+x, S-x, D+x, N-x |
| odd:odd | 3 | 3 | 0/2 | -90 | 0 / 180 | A+x, S+x, D+x, N+x |
| odd:even | 1 | 0 | 1/3 | +90 | +90 / -90 | A+x, N+y, D-x, S-y |
| odd:even | 1 | 2 | 3/1 | +90 | -90 / +90 | A+x, N-y, D-x, S+y |
| odd:even | 3 | 0 | 3/1 | -90 | -90 / +90 | A+x, S-y, D-x, N+y |
| odd:even | 3 | 2 | 1/3 | -90 | +90 / -90 | A+x, S+y, D-x, N-y |
| even:odd | 0 | 1 | 3/3 | 0 | -90 / -90 | N+x, N-y, N-x, N+y |
| even:odd | 0 | 3 | 1/1 | 0 | +90 / +90 | N+x, N+y, N-x, N-y |
| even:odd | 2 | 1 | 1/1 | 180 | +90 / +90 | N+x, S+y, N-x, S-y |
| even:odd | 2 | 3 | 3/3 | 180 | -90 / -90 | N+x, S-y, N-x, S+y |

(Units of the steps are degrees. The paper prints each sample sequence twice, the second copy as coloured markers.)

Table 5 (Evolution of symmetric PO configurations; numbers in parentheses give the count of distinct symmetric
configurations; DS doubly symmetric, SS singly symmetric, None no symmetric solutions; verified against the page image):

| Parity | Pairing (quarter period) | Circular HR3BP | Circular CR3BP | Frozen HR3BP | Frozen CR3BP |
|---|---|---|---|---|---|
| odd:odd | OX-XOZ | DS: OX/XOZ (2) | DS: OX/XOZ (2) | SS: XOZ (4) | SS: XOZ (4) |
| odd:odd | OY-YOZ | DS: OY/YOZ (2) | None | SS: YOZ (4) | None |
| odd:even | OX-YOZ | DS: OX/YOZ (2) | SS: OX (2) | SS: YOZ (4) | None |
| odd:even | OY-XOZ | DS: OY/XOZ (2) | SS: XOZ (2) | SS: XOZ (4) | SS: XOZ (4) |
| even:odd | OX-OY | DS: OX/OY (2) | SS: OX (2) | None | None |
| even:odd | XOZ-YOZ | DS: XOZ/YOZ (2) | SS: XOZ (2) | DS: XOZ/YOZ (4) | SS: XOZ (4) |

Table 6 (Initialisation lookup, rows u, columns Omega = chi_H): u = 0 (Ascending): A+x, A+y, A-x, A-y at
Omega = 0, 90, 180, 270 degrees; u = 90 (North): N+x, N+y, N-x, N-y; u = 180 (Descending): D+x, D+y, D-x, D-y;
u = 270 (South): S+x, S+y, S-x, S-y. Grey (pruned) columns for the CR3BP: Omega = 90 and 270; grey rows for frozen
equilibria: u = 0 and 180. Initial phase t0 = +/- u.

Individual orbit data actually printed (all examples, no initial-condition vectors):
- Circular example: p = 4, q = 1 (even:odd), i = 50 degrees gives a_H ~ 0.37 nd; two doubly symmetric XOZ/YOZ
  geometries in the HR3BP initialised at (Omega, M) = (0, 90) and (0, 270) degrees (Figs. 8-9); in the CR3BP
  (Earth-Moon) they degenerate to singly symmetric XOZ (Fig. 10).
- Frozen example: p = 5, q = 1, i = 50 degrees gives a_H ~ 0.30 nd; four singly symmetric XOZ configurations
  (omega = -90 and +90 degrees, two apsidal sequences each; Figs. 11-13).
- Other printed values: i_crit- ~ 39.2 and i_crit+ ~ 140.8 degrees; a_H ~ 0.32 nd threshold above which
  circular equilibria "are more appropriately identified as linking to the DPO family within the CR3BP"; period-doubling
  s_1 = -1 at a_H ~ 0.38 nd (collision/L2-halo comparison); collision orbit example z_H = -0.4 nd (Fig. 16(a)), southern
  collision orbit at a_H = 0.2 nd; retrograde asymptotic connection to i -> 90 degrees only for eta = p/q >= 12; p = 6,
  q = 1 CR3BP broken-bifurcation case (Fig. 17); preliminary note "retrograde 10:1 to prograde 5:1" connection by
  period doubling. Limits (Eqs. 48-49): nu_m,90- = 1 + (15/(4 nu_s)) (1-mu) sqrt(3/5); nu_m,90+ = 1 - (15/(4 nu_s)) (1-mu)
  sqrt(3/5).
No table lists initial conditions, periods, Jacobi constants or stability indices of specific orbits; those appear
only as plots (Figs. 14-17).

## Central results (quotes)
- Multiplicity and symmetry prediction: "Leveraging the parity of the resonance ratio, an initialization scheme
  is developed to identify admissible apse configurations, enabling the a priori prediction of solution
  multiplicity and symmetry types." Each circular equilibrium "lifts to four doubly symmetric PO geometries in the
  HR3BP (two per symmetry pairing) that the CR3BP subsequently prunes"; "the odd-odd case retains two doubly
  symmetric geometries (OX/XOZ), whereas the odd-even and even-odd cases are reduced to two singly symmetric
  geometries each". Frozen equilibria: "Axially symmetric configurations cannot be generated from frozen
  equilibria"; they have higher multiplicity (omega = +/-90 degrees doubles the count).
- Planar limits: "nominally, symmetric periodic orbits originating from the circular equilibria with a
  resonance ratio p : q connect to g/LPO families via a period-(p - q)-multiplying bifurcation at i = 0 degrees,
  and connect to f/DRO families via a period-(p + q)-multiplying bifurcation at i = 180 degrees".
- Polar limit: "nominally, symmetric periodic orbits originating from the frozen equilibria with a resonance ratio
  p : q connect to collision (HR3BP) and L2 halo orbits (CR3BP) via a period-p-multiplying bifurcation as i -> 90
  degrees."
- Broken bifurcations: in the HR3BP family g (triply symmetric) undergoes a pitchfork giving family g'; "Since the CR3BP
  does not possess rho_Y-symmetry, the ideal pitchfork bifurcation observed within the HR3BP cannot persist [23].
  Rather, the dynamics demonstrate a 'broken' bifurcation"; "the LPO family within the CR3BP appears to merge the
  members of family g at low a_H values with the 'western' branch of family g'"; the rest forms the DPO family. The
  retrograde family f continues to the DRO family "essentially avoiding the complexity of broken bifurcations".
  Likewise the HR3BP collision-to-halo pitchfork becomes two disconnected L1 and L2 halo branches in the CR3BP, and
  the p = 6, q = 1 circular-to-frozen bifurcation is shown to be broken in the CR3BP (Fig. 17).
- Link to Vertical Self-Resonant (VSR) orbits: "VSR orbits are spatial symmetric POs that originate from the DADM circular
  equilibria" (complementing Robin-Markellos [39], Aydin [24], Aydin-Batkhin [25], Peng et al. [31]).
- Averaged-to-unaveraged agreement: strong at a_H << 1; retrograde families show "remarkable agreement"; prograde
  disparities grow with a_H as orbits deform.

## Answers to the questions
1. Models: doubly averaged (DADM), HR3BP, CR3BP; spatial; Hill and circular restricted, not elliptic. Systems: only
   Earth-Moon (mu ~ 0.01215) is computed in the CR3BP; HR3BP is the mu -> 0 limit. Planet-moon systems mentioned: Earth's
   Moon; icy satellites of Jupiter and Saturn (motivation); Europa (cited science-orbit works). Constants for any system
   other than mu ~ 0.01215 are not stated.
2. Families: circular POs (e = 0) and frozen POs (e > 0, omega = +/-90 degrees) lifted from DADM equilibria; planar limits
   family g/LPO and g' (prograde), f/DRO (retrograde), collision orbits, L1 and L2 halo families. Classification: symmetry
   type (OX, XOZ, OY, YOZ; singly or doubly symmetric; HR3BP and CR3BP differ), multiplicity (Table 5), and resonance
   ratio eta = nu_s/nu_m = p/q of the orbit's short-period (Keplerian) motion to the rotating-frame nodal motion; the
   period is P = p (2 pi/nu_s). The resonance is with the rotation of the primaries' line (the secondary's period in
   the rotating frame), as p revolutions per q cycles of the negative RAAN in the rotating frame; period expressed in
   nondimensional time units.
3. See "Central results".
4. See "Exact tables".
5. Counts and quotes follow below.
6. Methods: DADM analytical lift and frequency mapping; modulo-4 parity analysis for symmetric apse configurations;
   differential correction targeting apse crossings at P/2 or P/4 (Robin-Markellos procedure); one-parameter family
   continuation; monodromy-matrix stability indices; bifurcation detection via stability index (s = +/-1) and
   comparison with DADM-predicted indices; KS regularisation for collision orbits; segment discretisation for long
   periods. Software and step sizes are not stated.

## Text-search counts (ligatures normalised, lines joined)
cycler 0; homoclinic 2; heteroclinic 0; "symmetric periodic" 23; resonan 60; Triton 0; Neptune 0; Europa 7 (introduction
"and at Europa [15]" plus six in bibliography titles/context); Titan 0; Moon 15 (Earth-Moon mu, "Earth's Moon", lunar
citations; "Earth-Moon" about 5 in the text); Uranus 0; torus 2; tori 10.
- Homoclinic (both uses are the averaged Lidov diagram, not three-body connections): "A homoclinic separatrix (black curve)
  partitions the circulating and librating motions of omega" (Fig. 2(c)); "The level set C_H = 0 also contains the homoclinic
  orbit characterized by sin^2 i sin^2 omega = 2/5".
- Heteroclinic: the word does not occur; no heteroclinic connections are discussed. Invariant manifolds are mentioned as
  motivation: "integrability removes the natural transport mechanisms, i.e., stable/unstable manifolds" and "invariant tori,
  whose families and invariant manifolds supply the backbone of the solution space"; none computed.
- Passing close to the secondary / flybys: not stated. The orbits are bounded orbiters near the smaller primary ("orbiters,
  i.e., ballistic trajectories that remain bounded near a smaller primary"); collision (rectilinear) orbits approach the
  secondary's centre by construction, and the periapsis behaviour of high-e frozen orbits is plotted, but no repeated flyby,
  cycler or resonant-flyby framing is made. "cycler" does not occur.
- Neptune, Triton, Uranus's moons: absent. Jupiter's moons: Europa (as cited application) and Ganymede (JUICE reference).
  Saturn: "Icy satellites of Jupiter and Saturn", Enceladus reference [3]. Earth-Moon: the only computed system.

## What it does NOT contain
No tables of orbit initial conditions, periods, Jacobi constants or stability-index values; no planet-moon system other than
Earth-Moon computed; no Neptune-Triton, Uranus or Jovian-moon constants; no homoclinic or heteroclinic connections of
periodic orbits, no invariant-manifold computation, no tori computed (only a dimensional-lift taxonomy); no resonant orbits
in the sense of mean-motion resonance with the secondary's orbit about the primary (here "resonance" is between the
Keplerian revolution about the small primary and the nodal rotation in the rotating frame); no cycler, flyby or transfer
design; no elliptic or ephemeris model; no quasi-periodic orbits ("Future efforts may extend this ... to quasi-periodic orbits
or higher-fidelity ephemeris models").

## Relevant bibliography (in full as printed)
- [19] K. C. Howell and J. V. Breakwell. Almost rectilinear halo orbits. Celestial mechanics, 32(1):29-52, 1984.
- [20] E. M. Zimovan-Spreen, K. C. Howell, and D. C. Davis. Near rectilinear halo orbits and nearby higher-period dynamical
  structures: orbital stability and resonance properties. Celestial Mechanics and Dynamical Astronomy, 132(5):28, 2020.
- [21] R. P. Russell. Global search for planar and three-dimensional periodic orbits near Europa. The Journal of the
  Astronautical Sciences, 54(2):199-226, 2006.
- [22] R. L. Restrepo and R. P. Russell. A database of planar axisymmetric periodic orbits for the solar system. Celestial
  Mechanics and Dynamical Astronomy, 130(7):49, 2018.
- [23] A. Moreno, C. Aydin, O. van Koert, U. Frauenfelder, and D. Koh. Bifurcation graphs for the CR3BP via symplectic
  methods: On the Jupiter-Europa and Saturn-Enceladus systems. The Journal of the Astronautical Sciences, 71(6):51, 2024.
- [24] C. Aydin. Exploration of vertical self-resonant bifurcations from DRO in the Earth-Moon CR3BP. arXiv:2508.17286, 2025.
- [25] C. Aydin and A. Batkhin. Studying network of symmetric periodic orbit families of the Hill problem via symplectic
  invariants. Celestial Mechanics and Dynamical Astronomy, 137(2):1-77, 2025.
- [26] K. C. Howell, D. J. Grebow, and Z. P. Olikara. Design using Gauss' perturbing equations with applications to lunar
  south pole coverage. AAS/AIAA Space Flight Mechanics Meeting, Sedona, Arizona, January 2007.
- [27] R. P. Russell and M. Lara. Long-lifetime lunar repeat ground track orbits. JGCD, 30(4):982-993, 2007.
- [28] M. Lara and J. F. San Juan. Dynamic behavior of an orbiter around Europa. JGCD, 28(2):291-297, 2005.
- [29] M. Lara, R. Russell, and B. Villac. Classification of the distant stability regions at Europa. JGCD, 30(2):409-418, 2007.
- [30] D. C. Koblick and P. Kelly. Novel three-body tulip-shaped orbit families for lunar missions. The Journal of the
  Astronautical Sciences, 72(4):32, 2025.
- [31] L. Peng, P. Shi, Y. Liang, and N. Pushparaj. Analog of lunar Sun-synchronous orbits based on spatial distant retrograde
  orbits. JGCD, 48(6):1439-1448, 2025.
- [32] J. Zhang, X. Jiang, Y. Yuan, and H. Dai. Time-regularized bifurcation framework for constructing tulip orbits in the
  CRTBP. Nonlinear Dynamics, 114(9):620, 2026.
- [38] C. Aydin. The linear symmetries of Hill's lunar problem. Archiv der Mathematik, 120(3):321-330, 2023.
- [39] I. A. Robin and V. V. Markellos. Numerical determination of three-dimensional periodic orbits generated from vertical
  self-resonant satellite orbits. Celestial Mechanics, 21(4):395-434, 1980.
- [40] M. Henon. Generating families in the restricted three-body problem. Springer, 2002.
- [41] E. Stromgren. Connaissance actuelle des orbites dans le probleme des trois corps. Bulletin astronomique, Observatoire
  de Paris, 9(1):87-130, 1933 (English translation by Maxime Murray and J. D. Mireles James, 2018).
- [42] D. Guzzetti, N. Bosanac, A. Haapala, K. C. Howell, and D. C. Folta. Rapid trajectory design in the Earth-Moon ephemeris
  system via an interactive catalog of periodic and quasi-periodic orbits. Acta Astronautica, 126:439-455, 2016.
- [43] C. Joung, D. Koh, and O. van Koert. Bifurcations of highly inclined near halo orbits using Moser regularization.
  arXiv:2512.03849, 2025.
- [16] R. P. Russell and A. T. Brinckerhoff. Circulating eccentric orbits around planetary moons. JGCD, 32(2):424-436, 2009.
- [17] S. McArdle and R. P. Russell. Circulating, eccentric periodic orbits at the Moon. Celestial Mechanics and Dynamical
  Astronomy, 133(4):18, 2021.
- [18] C. J. Franz and R. P. Russell. Database of planar and three-dimensional periodic orbits and families near the Moon.
  The Journal of the Astronautical Sciences, 69(6):1573-1612, 2022.
Remaining references [1]-[15], [33]-[37] are mission overviews, averaging theory and general texts (Gateway status, JUICE,
Enceladus, Szebehely, von Zeipel, Lidov, Kozai, Ito-Ohtsuka, Scheeres et al., Prado, Broucke, Ely, Folta-Quinn, Ely-Lieb,
Lara-Russell, Vallado, Longuski et al., Nie-Gurfil, Arnold, Celletti).

## Suggested topology label
**halo** is the closest fit among the project's labels only for the collision-to-halo and DRO/DPO connections; the paper's
main content (symmetric resonant periodic orbits near the smaller primary and their bifurcation web) is best filed under
**resonant** with a note that the resonance is between the Keplerian revolution and the rotating-frame nodal drift, not a
mean-motion resonance with a moon's orbital period about a planet in the cycler sense. Suggested: resonant.

## Corpus index line
Park and Howell 2026 (arXiv 2606.08485v1): doubly averaged model linked to HR3BP and CR3BP symmetric periodic orbits near
the smaller primary; parity-of-p:q rule for symmetry type and branch count (Table 5), planar g/LPO/DPO/DRO and
collision/halo connections, broken bifurcations; Earth-Moon mu ~ 0.01215 only; no orbit data tables, no connections.
