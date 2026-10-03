# Digest — Guido & Efthymiopoulos (2026), "Arches of chaos, heteroclinic connections of first-order MMRs and the chaotic transport of small bodies in the Sun-Jupiter system" (arXiv:2604.00679v1)

**Digested:** 2026-10-03. Text-layer PDF, 20 pages, read in full via `pdftotext -layout`
(figures not rendered; their captions were read). Filed in the private paper corpus as
`guido-efthymiopoulos-2026-arches-of-chaos-heteroclinic-connections-first-order-mmr-sun-jupiter-arxiv-2604.00679v1.pdf`
(md5 `5c1a96972eb8c67de51cfcb317c10bd5`).

## Citation as printed
Alessia Francesca Guido (1 Dept. of Physics, University of Trento; 2 Dept. of Mathematics,
University of Rome Tor Vergata) and Christos Efthymiopoulos (3 Dept. of Mathematics Tullio
Levi-Civita, University of Padova), "Arches of chaos, heteroclinic connections of first-order MMRs
and the chaotic transport of small bodies in the Sun-Jupiter system", arXiv:2604.00679v1
[astro-ph.EP], 1 Apr 2026. Preprint; no journal is named in the text. Keywords as printed:
"mean motion resonances, chaotic transports, invariant manifolds, CRTBP, FLI maps".

## What the paper is
A numerical celestial-mechanics study in the Sun-Jupiter planar restricted three-body problem. It
computes, on a Poincare section at fixed Jacobi energy, the stable and unstable manifolds of the
unstable periodic orbits of first-order interior (2:1, 3:2) and exterior (2:3) mean-motion
resonances (MMRs) and of the short-period Lyapunov orbit PL3 at L3, overlays them on short-time
Fast Lyapunov Indicator (FLI) maps, and reads off a network of heteroclinic connections between
resonances. It uses these to explain the "arches of chaos" of Todorovic, Wu and Rosengren (Science
Advances 2020, ref [52]) and "resonance hopping" of small bodies. Purpose stated is small-body
chaotic transport (quasi-Hildas, Jupiter-family comets, Centaurs), not spacecraft design.

## Model and constants
- Planar circular restricted three-body problem (CRTBP), Sun-Jupiter, Hamiltonian in barycentric
  synodic cylindrical coordinates (r, phi), Eq. (1). Most results are in this model. The planar
  elliptic RTBP (ERTBP) is used only for FLI maps in Section 5. Spatial motion is not treated
  ("our study is limited to the planar RTBP").
- Constants as printed: M_S = 1 M_sun; M_J = 0.00096 M_sun; a_J = 5.19 AU; G = 4 pi^2 AU^3/(yr^2 M_sun);
  Omega_J = [G (M_S + M_J)/a_J^3]^(1/2). The text also says "Jupiter (mu = 0.001)" in the
  introduction. Mass ratio mu is not given as a separate number for the computations.
- Elliptic model: e_J = 0.0487, Jupiter launched from pericentre at t = 0 (Eq. 6 and following).
- Jacobi relations as printed: T_J = -2 C_J (Eq. 3); C_J = (a_J/(G(M_S+M_J))) H = (a_J/(G(M_S+M_J))) E_J (Eq. 4);
  T_J = a_J/a + 2 sqrt((a/a_J)(1-e^2)) (Eq. 2).
- Energy range studied, as printed: "−1.4 < C_J < −1.525, which corresponds to Tisserand parameters
  2.8 < T_J < 3.05" (the inequality is printed in that order; read as C_J between -1.525 and -1.4).
  Representative levels: E_J = C_J = -1.48 (T_J = 2.96); -1.42 (T_J = 2.84); -1.52 (T_J = 3.04).
  Fig. 5 uses 40 equidistant levels from E_J = -1.52 to -1.40.
- Integration: "standard fourth-order ODE integrator with adjustable timestep in Matlab, with control
  of the precision up to about ten significant figures"; variational equations integrated alongside.
- Resonances treated: first order |q-p| = 1: interior 2:1, 3:2 (4:3 appears in the phase portrait),
  exterior 2:3 (1:2 appears in the phase portrait and as a reached resonance), plus the co-orbital
  1:1 (PL3, PL4/PL5, quasi-satellite, PL1/PL2). Second-order 5:3 and 3:5 appear in Fig. 1 only as
  islands of stability. Convention: "q > p for interior resonances, q < p for exterior resonances".

## Objects computed and methods
- Poincare map: pericentric surface of section P_EJ = {p_r = 0, (phi, p_phi) in T x D_EJ}, with r the
  root of H(phi, r, p_r = 0, p_phi) = E_J with dp_r/dt > 0 (Eq. 7); prograde (p_phi > 0). phi = g - lambda_J
  (argument of pericentre minus Jupiter's mean longitude, rotating frame), p_phi = sqrt(G a (1-e^2)) in
  the text's notation. Said to be "similar to the one used in [29]" (Koon-Lo-Marsden-Ross 2000), with G in
  place of L.
- Periodic orbits P_q:p: "the unstable periodic orbit of the resonance q:p has multiplicity m = q in
  the surface of section P_EJ, hence it corresponds to a fixed point of the q-th iterate of the Poincare
  return map". Found by Newton-Raphson; Jacobian of the return map by central finite differences,
  symplectic to "within eight significant figures"; monodromy matrix A_q:p = product of q one-step
  Jacobians; eigenvalues real with lambda_1 lambda_2 = 1.
- Manifold curves C^S_q:p,EJ and C^U_q:p,EJ: initial conditions along a segment of length
  Delta S = 1e-3 on the unstable (stable) eigendirection, iterated with the forward (backward) return map.
  No parameterisation-method manifolds, no continuation in energy of the periodic orbits themselves
  (energies are sampled), no spatial manifolds.
- Short-time FLI maps: FLI(T) = max_{0<=t<=T} log10 ||xi(t)||, T = 100 yr, 300 x 300 grids of
  (phi in [0, 2 pi], p_phi in [9, 16]); no window function. Forward integration shows stable manifolds,
  backward shows unstable manifolds (method of Guzzo and Lega, ref [24]).
- Resonances whose manifolds are explicitly computed: 2:1, 3:2, 2:3 and PL3. Resonance 1:2 and 4:3
  are seen only in FLI maps.

## Exact numbers printed (no data tables exist in the paper)
The paper contains no numbered tables. No periodic-orbit initial conditions, periods, eigenvalues or
Jacobi constants of individual resonant orbits are printed. Everything quantitative is listed here:
- Energy levels: E_J = -1.48 (T_J = 2.96); -1.42 (2.84); -1.52 (3.04); range -1.52 to -1.40 (Fig. 5).
  Statement "Since E_J > E_J,L3 ..." and "for values E_J < E_J,L3, communication between the interior and
  exterior MMRs is no longer energetically allowed"; the numerical value of E_J,L3 is not printed.
  Fig. 1 text: "C_J = -1.48 ... slightly above the one of the Lagrangian points L1, L2".
- Quasi-Hilda T_J range quoted from Toth (ref [54]): 2.9 to 3.05; evolved JFCs 2.8 < T_J < 3.1; scattered
  disc (Rickman et al. [46], Nice model): 2.4 < T_J < 3.
- Heteroclinic paths named in the figures: PL1D, PL1U, PL2U (classical Koon et al. paths); P2112
  (2:1 to 2:3 and even 1:2, avoiding PL1 and PL2 neighbourhoods; clear at T_J = 2.84, also present at 2.96
  at smaller distance from PL1D); P2132 (2:1 to 3:2, the only prominent connection type at T_J = 3.04).
- Section line phi = pi/3 with ridge points A, B, C, D, E (Fig. 5); pericentric section for the arches is
  l = 0, omega = pi/3 (the original [52] section was l = pi/3, omega = omega_Jupiter).
- Resonance-hopping time scale: "nearly all trajectories ... along the same section line are observed to
  undergo some resonant hopping in time intervals ranging between few years and few decades" (Fig. 6,
  three example trajectories shown as n/n_J = (a_J/a)^(3/2) against time; plateaus identify the resonance).
- Figure 5 mapping: p_phi = sqrt(G M_S a (1-e^2)) (printed with M_S), a = -G M_S/(2 E_K), E_K = p_phi^2/(2 r^2) - G M_S/r.
- Fig. 1 islands visible at C_J = -1.48: interior 2:1, 3:2, 4:3; exterior 2:3, 1:2; second-order 5:3 and 3:5.

## What "arches of chaos" means here and how heteroclinic connections relate
- Definition used: "According to [52], the 'arches of chaos' are manifold structures observed in the
  projection to the plane (a, e) (semimajor axis - eccentricity) of a particular section of the phase
  space", visualised by forward-integrated FLI maps (Todorovic et al. 2020 use l = pi/3, omega = omega_Jupiter).
- Result: the arches are the same invariant structures as the ridges of the (phi, p_phi) FLI maps: each
  arch is the manifold of the corresponding MMR that "locally gives rise to a ridge in the FLI map";
  joining ridge/manifold points across the 40 V-shaped (a, e) curves at successive E_J gives the arches.
  The structures are described as a "unique Lagrangian Coherent Structure" (LCS) of the flow,
  citing Gawlik et al. (ref [22]).
- Heteroclinic role: intersections of stable curves of one resonance with unstable curves of another define
  heteroclinic orbits; chains give "resonance hopping". Direct (not through the co-orbital resonance)
  connections between 2:1 and 2:3 are claimed (property (ii) in the introduction); indirect ones go through PL3
  (2:1 to L3 to 2:3). Abstract: "we observe direct heteroclinic connections between the manifolds of the
  interior with exterior MMRs, which do not involve the manifolds of any of the periodic orbits of the
  co-orbital resonance". The text also claims evidence "of heteroclinic connections between the manifolds of
  the short-period orbits around L3 and of the periodic orbits associated with interior or exterior first
  order MMRs".
- Statement on the section 4.1 (properties (i)-(iii) in Fig. 3 discussion): the union of the manifolds of
  2:1, 3:2, 1:1 (PL3), 2:3 covers "essentially the entire chaotic saddle observed by the numerical FLI maps";
  PL3 manifolds "penetrate deeply in the domains of both the Hilda (3:2) and Hecuba (2:1) interior MMRs".
  Fig. 3 caption prints "Top right" twice (typesetting slip); the last panel is the 2:1 versus 2:3 comparison.
- Energy dependence: direct inner-outer connections become prevalent as E_J rises (T_J decreases below 3).
  At T_J = 3.04 interior and exterior domains are energetically separated and only interior-interior
  connections (e.g. P2132) occur.
- ERTBP (Section 5): the FLI-map structure "very similar" to the circular case; "most heteroclinic paths such
  as the paths PL1U, PL1D marked in Fig. 1 survive with their structure essentially unaltered".

## Answers to the questions
1. Model: planar circular RTBP, Sun-Jupiter (M_J = 0.00096 M_sun, a_J = 5.19 AU), Jacobi energies in
   -1.52 to -1.40 (T_J 2.8 to 3.05); planar ERTBP (e_J = 0.0487) for FLI maps only. First-order MMRs:
   interior 2:1 and 3:2 (4:3 seen), exterior 2:3 (1:2 seen), plus the co-orbital 1:1.
2. Computed: fixed points of the q-th return map for 2:1, 3:2, 2:3 and PL3, their stable and unstable
   one-dimensional curve manifolds on the pericentric section at fixed E_J (Newton-Raphson, finite-difference
   Jacobian, eigendirection segments of length 1e-3); FLI maps; heteroclinic connections between 2:1 and 3:2
   (direct), 2:1 and 2:3 (direct and via PL3), 3:2 / 2:3 via PL1, PL2 (classical), and PL3 with 2:1, 3:2, 2:3,
   3:4. Not by continuation in energy and not by a parameterisation method.
3. Reproduction targets: only the energy levels, constants and path labels above. No periodic-orbit initial
   conditions are tabulated.
4. See previous section.
5. Homoclinic: the word occurs twice, both general. Quote: "Any point of intersection of the curves
   C^S_q:p,EJ, C^U_q:p,EJ defines the initial condition of a homoclinic orbit in the Poincare section P_EJ."
   And: "the manifolds' homoclinic lobes are rather constrained around the resonance" (for the 2:1 resonance
   when E_J < E_J,L3). No homoclinic orbit, and no family of symmetric periodic orbits accumulating on one,
   is computed or tabulated. The phrase "symmetric periodic" does not occur. The orbits P_q:p are described
   only as unstable resonant periodic orbits; their symmetry is not discussed. Tori: "partially hyperbolic
   invariant tori, which generalize the periodic orbits P_q:p of the circular problem" are mentioned as needed
   for the elliptic problem and not computed.
6. Planetary-moon systems: Earth-Moon appears only in citations (counts below): "chaotic transitions through
   MMRs for the cislunar orbits of artificial bodies in the Earth-Moon system ([3],[53],[23],[33],[31])" and
   "heteroclinic connections similar to those here discussed were reported also in the Earth-Moon-particle
   problem ([31]), thus their appearance in the framework of the restricted three body problem appears rather
   generic, for sufficiently large values of the mass parameter mu." Neptune, Triton, Europa, Titan, Uranus:
   absent.

## Text-search counts (ligatures normalised, lines joined)
cycler 0; homoclinic 2; heteroclinic 43; "symmetric periodic" 0; resonan 93; Triton 0; Neptune 0; Europa 0;
Titan 0; Moon 5 (all Earth-Moon: three in the introduction, including "Sun-Jupiter or Earth-Moon RTBP", and two in
bibliography titles [27], [31]); Uranus 0; torus 0; tori 2 (both the sentence about partially hyperbolic invariant tori).

## What it does NOT contain
No spacecraft or cycler trajectories; no flyby or encounter bookkeeping; no tabulated initial conditions,
periods, monodromy eigenvalues or Jacobi constants of the resonant orbits; no continuation of orbit
families; no homoclinic connection or symmetric-orbit family computed; no spatial or planetary-moon
computation; no mass-ratio study (single Sun-Jupiter value); no quantitative transit-time statistics beyond
"few years and few decades"; no journal named; no tables.

## Relevant bibliography (in full as printed, resonant / heteroclinic / manifold works)
- [2] E. Belbruno and B. G. Marsden. Resonance Hopping in Comets. Astronomical Journal, 113:1433, 1997.
- [3] E. Belbruno, F. Topputo, and M. Gidea. Resonance transitions associated to weak capture in the
  restricted three-body problem. Advances in Space Research, 42(8):1330-1351, 2008.
- [9] M. Dellnitz, O. Junge, M. Lo, J. Marsden, K. Padberg-Gehle, R. Preis, S. Ross, and B. Thiere. Transport
  of mars-crossing asteroids from the quasi-hilda region. Physical Review Letters, 94, 2005.
- [16] J. Fitzgerald and S. Ross. Geometry of transit orbits in the periodically-perturbed restricted
  three-body problem. Advances in Space Research, 70, 2022.
- [22] E. Gawlik, J. Marsden, P. Du Toit, and S. Campagnola. Lagrangian coherent structures in the planar
  elliptic restricted three-body problem. Celestial Mechanics and Dynamical Astronomy, 103:227-249, 2009.
- [23] B. R. Gupta and V. Kumar. Characterization of the Phase Space Structure of Circular Restricted
  Three-Body Problem: An Alternative Approach. International Journal of Bifurcation and Chaos,
  26(2):1650029-1351, 2016.
- [24] M. Guzzo and E. Lega. Evolution of the tangent vectors and localization of the stable and unstable
  manifolds of hyperbolic orbits by fast lyapunov indicators. SIAM J. Appl. Math., 74:1058-1086, 2013.
- [25] M. Guzzo and E. Lega. Theory and applications of fast Lyapunov indicators to model problems of
  celestial mechanics. Celestial Mechanics and Dynamical Astronomy, 135(4):37, 2023.
- [27] A. Jorba and Begona N. Transport and invariant manifolds near l3 in the earth-moon bicircular model.
  Communications in Nonlinear Science and Numerical Simulation, 89:105327, 2020.
- [29] W. Koon, M.W. Lo, J.E. Marsden, and S.D. Ross. Heteroclinic connections between periodic orbits and
  resonance transitions in celestial mechanics. Chaos, 10:427-469, 2000.
- [30] W. Koon, M.W. Lo, J.E. Marsden, and S.D. Ross. Resonance and Capture of Jupiter Comets. Celestial
  Mechanics and Dynamical Astronomy, 81:27-38, 2001.
- [31] B. Kumar, A. Rawat, A. Rosengren, and S. Ross. Investigation of interior mean motion resonances and
  heteroclinic connections in the earth-moon system. 2024.
- [33] H. Lei and B. Xu. Resonance transition periodic orbits in the circular restricted three-body problem.
  Astrophysics and Space Science, 363(4):70, 2018.
- [42] R. Paez and M. Guzzo. On the semi-analytical construction of halo orbits and halo tubes in the elliptic
  restricted three-body problem. Physica D: Nonlinear Phenomena, 439:133402, 2022.
- [43] S. Pan and X. Hou. Analysis of resonance transition periodic orbits in the circular restricted
  three-body problem. Applied Sciences, 12, 2022.
- [51] G. Tancredi, M. Lindgren, and H. Rickman. Temporary satellite capture and orbital evolution of Comet
  P/Helin-Roman-Crockett. Astronomy and Astrophysics, 239(1-2):375-380, 1990.
- [52] N. Todorovic, D. Wu, and A. Rosengren. The arches of chaos in the solar system. Science Advances, 6, 2020.
- [53] F. Topputo, E. Belbruno, and M. Gidea. Resonant motion, ballistic escape, and their applications in
  astrodynamics. Advances in Space Research, 42(8):1318-1329, 2008.
- [55] J. C. Urschel and J. R. Galante. Instabilities in the sun-jupiter-asteroid three body problem. Celestial
  Mechanics and Dynamical Astronomy, 115:233-259, 2013.
- [56] X. Wang and R. Malhotra. Mean Motion Resonances at High Eccentricities: The 2:1 and the 3:2 Interior
  Resonances. Astronomical Journal, 154(1):20, 2017.
- [41] R. Paez and C. Efthymiopoulos. Trojan resonant dynamics, stability, and chaotic diffusion, for
  parameters relevant to exoplanetary systems. Celestial Mechanics and Dynamical Astronomy, 121(2):139-170, 2015.
- [57] J. Wisdom. The resonance overlap criterion and the onset of stochastic behavior in the restricted
  three-body problem. Astronomical Journal, 85:1122-1133, 1980.
Other references (asteroid/comet populations, resonance-overlap theory) are small-body astronomy and are
not listed here; they are in the paper's reference list, numbers [1]-[57].

## Suggested topology label
**resonant** (resonant periodic orbits and their manifold-mediated connections): the objects are unstable
resonant periodic orbits of first-order MMRs and heteroclinic connections between them, in a restricted
three-body model, with no encounter itinerary.

## Corpus index line
Guido and Efthymiopoulos 2026 (arXiv 2604.00679v1): planar circular Sun-Jupiter RTBP; explicit stable and
unstable manifolds of the 2:1, 3:2, 2:3 and PL3 unstable periodic orbits on a pericentric Poincare section at
T_J 2.84 to 3.04, overlaid on FLI maps; direct and PL3-mediated heteroclinic connections between interior and
exterior MMRs; explains the arches of chaos; no tables, no moon system, no symmetric-orbit families.
