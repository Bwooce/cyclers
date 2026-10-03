# Digest — Bonasera & Bosanac (2023), "Computing Natural Transitions Between Tori Near Resonances in the Earth-Moon System" (JGCD 46(3))

**Digested:** 2026-10-03 (text-layer PDF; all 12 pages read, pages 3-11 also checked against the
rendered images for equations, figure axes and the printed numbers). **Citation:** Stefano
Bonasera and Natasha Bosanac, University of Colorado Boulder, "Computing Natural Transitions
Between Tori Near Resonances in the Earth-Moon System," Journal of Guidance, Control, and
Dynamics **46(3)**, March 2023, pp. 443-454, **DOI 10.2514/1.G006941**. Received 4 May 2022,
accepted 16 December 2022, published online 31 January 2023. Filed in the private paper corpus
as `bonasera-bosanac-2023-computing-natural-transitions-between-tori-near-resonances-earth-moon-jgcd-46-3-443-doi-10.2514-1.G006941.pdf`
(12 pages, md5 `8edd91688b8c90cc1393c5d96b27f414`). **Corpus index description (one line):**
Heteroclinic natural transfers between spatial 2-tori near Earth-Moon resonances in the CR3BP,
initial guesses from UMAP-projected Poincare crossings, multiple-shooting correction and
two-parameter continuation (3:2 to 1:2 at C_J = 2.73; also 3:1 to 1:3 and 2:3 to 1:5).

The acknowledgements state that the work also appears in the first author's PhD dissertation,
"Incorporating Machine Learning into Trajectory Design Strategies in Multi-Body Systems"
(University of Colorado Boulder, April 2022).

## What it is
A method paper plus worked examples. It builds families of spatial invariant 2-tori near two
different mean-motion resonances at one common Jacobi constant, generates the hyperbolic
stable and unstable manifolds of one torus from each family, records their crossings of a
surface of section (a five-dimensional dataset), projects those crossings to two dimensions
with UMAP (uniform manifold approximation and projection) to pick nearby crossing pairs,
turns each pair into a discontinuous initial guess, corrects it with a multiple-shooting
problem posed as a constrained optimisation, and then continues the result across the two
torus families. The paper uses "heteroclinic connections", "natural transfers" and "natural
transitions" interchangeably (Sec. V.A). There are no tables; all results are figures
(Figs. 1-11) and numbers in the text.

## Sections
I Introduction; II Dynamical Model; III Periodic Orbits Near Mean-Motion Resonances; IV
Quasi-Periodic Trajectories and Invariant 2-Tori; V Using Dimension Reduction to Visualize
Higher-Dimensional Poincare Maps (A Poincare Maps; B Manifold Learning); VI Computing Natural
Transitions Between Spatial Invariant 2-Tori (A Initial Guess Generation; B Trajectory
Correction and Continuation); VII Natural Transitions Between Tori near Distinct Resonances
(A Recovering Point Solutions for Natural Transfers; B Family Continuation; C Recovering
Transfers Between Invariant Tori near Additional Resonances); VIII Conclusions.

## (1) Model, constants, resonances, tori
- **Model:** the circular restricted three-body problem (CR3BP), Earth-Moon. Only constant
  stated: mass parameter mu ~ 0.01215 (Sec. II). Nondimensional rotating frame, origin at the
  barycentre, x from P1 to P2; equations of motion and pseudo-potential U written out
  (Eq. 1); Jacobi constant C_J = 2U - xdot^2 - ydot^2 - zdot^2. The dimensional length, time and
  velocity units are not printed (the paper nevertheless quotes dimensional equivalents: 100
  km, 400 m, 1 mm/s, and periods in days).
- **Resonance definition (Sec. III):** body B completes p orbits about A in the time C
  revolves q times about A; the p:q resonance is interior when p > q and exterior when p < q.
  In the CR3BP B is P3, A is P1 (Earth), C is P2 (Moon). Planar p:q periodic orbit initial
  guesses are built in the two-body problem following Vaquero and Anderson, transformed to the
  rotating frame, corrected by differential correction, and extended by pseudo-arclength
  continuation. The paper notes only a finite number of members of a "p:q family" are exactly
  resonant.
- **Resonance pairs treated:** (a) 3:2 (interior) to 1:2 (exterior) at C_J = 2.73 (the main
  example, Secs. VII.A-B); (b) 3:1 to 1:3 at C_J = 3 (Sec. VII.C, Fig. 10); (c) 2:3 to 1:5 at
  C_J = 2.6 (Sec. VII.C, Fig. 11). For (a) the two generating planar periodic orbits have
  "approximate periods of 55.92 and 50.54 days, respectively" (printed in the order "3:2 and
  1:2 resonances"; how the period relates to the resonant ratio is not stated).
- **Torus computation:** the approach of Jorba, Gomez-Mondelo and Olikara-Scheeres: an invariant
  curve satisfying R_(-rho) phi_T(s(theta1,theta2)) - s(theta1,theta2) = 0 under a stroboscopic
  map of time T = 2 pi / omega1 with rotation angle rho = 2 pi omega2 / omega1 (Eqs. 2-3);
  truncated Fourier series with an odd number N_Q of states equally spaced in theta2; initial
  guess s0 + eps (Re[v_C] cos theta2 + Im[v_C] sin theta2) (Eq. 4) from the complex-unit-
  eigenvalue eigenvector of the periodic orbit's monodromy matrix, T taken as the orbit
  period and rho = Re(-i ln lambda_C); multiple-shooting correction; **one-parameter family at
  fixed average Jacobi constant** (average over the states of the initial invariant curve);
  pseudo-arclength continuation along the family. Stability from the eigenstructure of DS
  (the differential of the invariance condition): "concentric circles" in the Gauss plane;
  non-unit radii give stable and unstable modes.
- **Numbers printed for the main example (Sec. VII.A):** eps = 5e-5 (Eq. 4 step); N_Q = 25
  states on each of M_Q + 1 = 4 invariant curves; 20 tori computed along each family, with
  the 20th "at the boundary of the range of tori that are computed with the selected invariant
  curve discretization, but additional tori may exist further along the family". For Sec.
  VII.C the tori use N_Q = 25 and M_Q = 101.
- **Torus accuracy:** not given as a number for the tori themselves. Stated qualitatively
  (Sec. IV): each torus "is not exactly recovered" because of the discrete Fourier transform,
  numerical integration and a nonzero tolerance; the average Jacobi constant is constrained
  but "states along each invariant curve" have a nonzero spread between minimum and maximum C_J.
  Rotation numbers, fundamental frequencies omega1, omega2 and tolerances of the torus solver
  are **not stated**. The tori are characterised only by their maximum out-of-plane distance at
  apolune (z_apo).
- **Integrator:** Runge-Kutta Prince-Dormand (8,9) in C++ from the GNU Scientific Library,
  absolute tolerance 1e-15, relative tolerance 1e-14 (Sec. VI.A step 1).
- **Torus z_apo range (read from the axes of Figs. 8 and 9, approximate):** departure (3:2)
  tori about 0.02 to about 0.11 (nondimensional); arrival (1:2) tori about 0.05 to about 0.32.
  Not printed as numbers in the text.

## (2) What a natural transition is and exactly how it is computed
**Definition:** an intersection of the stable manifold of one torus with the unstable
manifold of another torus; the resulting arc leaves the departure torus and approaches the
arrival torus naturally (no manoeuvre). In the numerical problem it is "an approximation to a
nearby heteroclinic connection" because of the numerical tolerances.

**Phase 1, initial guess (Sec. VI.A, five steps):**
1. Compute two families of tori near two distinct resonances at the same C_J (the starting
   periodic orbits must have both hyperbolic and centre manifolds).
2. Surface of section **y = 0, no restriction on the sign of any velocity component** (a
   two-sided map). (The preliminary planar illustration in Fig. 1b uses a one-sided map, ydot > 0,
   up to 10 crossings, plotted in (x, xdot).)
3. Choose one torus from each family "to possess similar maximum out-of-plane components"
   (empirically a good start). Generate each manifold with "a small displacement, equivalent
   to 100 km in the configuration space, along the stable and unstable eigenvectors"
   (the eigenvectors of DS), from invariant curves at **101 values of theta1** (Sec. VII.A).
   Whether all N_Q states of each curve are used, and whether both signs of the eigenvector
   are used, is **not stated**. Propagate up to 12 returns to y = 0; the manifolds linger near the
   generating torus for about six revolutions, so **only crossings 7 to 12 are kept**.
   Crossings recorded: **15,114** (3:2 unstable manifold) and **14,925** (1:2 stable manifold)
   (sum 30,039, my arithmetic).
4. Each crossing is the five-dimensional vector [x, z, xdot, ydot, zdot] (y = 0 removes one
   dimension; the Jacobi constant is deliberately not used to remove another, because the
   states on each invariant curve share only the average C_J). UMAP (Python umap-learn)
   projects the combined stable-plus-unstable set to two dimensions with **n_n = 100, m_dist = 0.0,
   n_c = 2**; the parameters were chosen manually, iteratively, to favour global structure and a
   compact embedding. The paper quotes UMAP complexity ~ O(N^1.14) for the graph step and
   ~ O(n_n N) for the optimisation.
5. Inspect the embedding by eye for regions where blue (stable-manifold) and magenta
   (unstable-manifold) points are close; take one crossing from each; propagate the unstable
   one backward and the stable one forward in time to form the arcs; **concatenate three
   revolutions of the associated torus at the beginning and end** of the transfer. The pair is
   a discontinuous initial guess. Four regions (Fig. 4) gave four transfers T1-T4 (Fig. 5); the
   region labelled T2 was re-examined and four more pairs produced four further transfers
   (Fig. 7), two geometry groups.

**Phase 2, correction (Sec. VI.B):**
- Multiple shooting with the free variable vector V = [s_1,...,s_M, t_(1,2),...,t_(M-1,M)]
  in R^(7M-1) (M states, M - 1 arcs, one propagation time per arc), constraint vector
  F(V) in R^(6(M-1)) enforcing full-state continuity (Eqs. 5-6).
- **Objective** f(V) = ||s_1 - s_T1||^2 + ||s_M - s_T2||^2 (Eq. 8): the squared distances from
  the first and last states of the transfer to the nearest points s_T1 and s_T2 on the
  departure and arrival tori. Problem: minimise f subject to F = 0 (Eq. 7).
- **Nearest torus point** recomputed at every optimiser iteration: the invariant curve is
  rotated by tau1 in the longitudinal direction (equivalent to forward propagation) and tau2 in
  the transverse direction, a two-variable single-shooting Newton iteration on
  G(Y) = s_1 - R(tau2) U(theta1 + tau1, theta2)|_1 with Y = (tau1, tau2), initial guess
  Y = (0.01, 0.01), **10 Newton steps** per evaluation; stated to work if the torus is described
  by a sufficiently large number of invariant curves.
- Solver: MATLAB fmincon, interior-point. The paper says an optimisation formulation "was
  observed to produce a less numerically sensitive process" than an equality-only
  multiple-shooting formulation and that mitigating the sensitivity of the latter is future work.
- **Acceptance thresholds:** a transfer is a natural connection if **f(V) <= 1e-12 and
  ||F(V)||_2 <= 1e-12** (Fig. 8 caption and text write "< 1e-12"). The paper states the
  objective threshold "corresponds to a cumulative approximate displacement of 400 m and
  1 mm/s in the Earth-Moon system" from the departure and arrival torus, justified by the
  nonzero spread of C_J along each approximate invariant curve.
- **Run time:** about 15 s per trajectory on an i7-2600K 3.40 GHz processor, "while remaining
  close to the initial guess" (Sec. VII.A).

**Continuation (Sec. VI.B and VII.B):** a grid. Hold the departure torus fixed and step the
arrival torus along its family, seeding each step with the previous solution; stop when the
family ends or no feasible transfer is found; then move to the next departure torus and repeat.
Only **one transfer per torus pair, near the initial guess's geometry**, is sought; the departure
and arrival longitudes are not themselves continued ("continuation is not used to find similar
transfers connecting the tori at various longitudinal and transverse angles"). Two families
were continued: those resembling T1 (Fig. 8) and T4 (Fig. 9) of the 3:2-to-1:2 example.

**Verification as one continuous trajectory:** by construction, the continuity constraint
F(V) = 0 to 1e-12 plus the torus-proximity objective below 1e-12. The corrected transfers'
y = 0 crossings are plotted over the manifold-crossing clouds (Fig. 6, projections onto
(x, xdot) and (x, xdot, z)) to show they sit near the projected manifold intersections.
No independent re-propagation of the full transfer from a single initial state, no
higher-fidelity (ephemeris) check, no residual table, and no time-of-flight or delta-V
statement is printed.

## (3) Autonomous or time-periodic; phase and energy matching
- **Autonomous** (CR3BP, circular primaries). No time-periodic model is used. The paper
  mentions Kumar, Anderson and de la Llave [17] (GPU search for connections between tori in
  periodically perturbed planar CR3BPs) as prior work only.
- **Energy matching:** both tori are constrained to the **same average Jacobi constant** at
  construction (2.73, 3, 2.6 in the three examples). Because only the average over the initial
  invariant curve is constrained, individual states differ slightly in C_J; this is why the
  five-dimensional (not four-dimensional) crossing space is used and why the 1e-12 tolerance is
  interpreted as about 400 m and 1 mm/s.
- **Phase matching:** there is no fixed-phase requirement. The departure and arrival points are
  free on the tori: the objective measures the distance to the nearest point on each torus,
  found by the (tau1, tau2) search above. The paper states departure and arrival locations
  could be varied and might enlarge the set of torus pairs that connect (not done).

## (4) Printed numbers that could be reproduction targets
No table of initial conditions, periods, rotation numbers, transfer times or delta-Vs exists.
Everything printed (transcribed from text and checked against the rendered pages):

| Quantity | Value as printed | Where |
|---|---|---|
| Mass parameter | mu ~ 0.01215 | Sec. II |
| Illustration orbits, C_J | 2.73 (3:2 interior and 1:2 exterior planar periodic orbits) | Sec. V.A, Fig. 1 |
| Periodic-orbit periods at C_J = 2.73 | 55.92 d and 50.54 d ("respectively", 3:2 then 1:2) | Sec. VII.A |
| Torus step eps | 5e-5 | Sec. VII.A |
| N_Q; M_Q + 1 (main example) | 25; 4 | Sec. VII.A |
| N_Q; M_Q (3:1 to 1:3 example) | 25; 101 | Sec. VII.C |
| Tori per family | 20 (the 20th at the edge of the computed range) | Sec. VII.A |
| theta1 samples for manifolds | 101 | Sec. VII.A |
| Manifold offset | 100 km equivalent | Sec. VI.A step 3 |
| Returns propagated; used | 12; 7th to 12th (main example). 18; last 8 (3:1 to 1:3) | Secs. VI.A, VII.A, VII.C |
| Crossings | 15,114 (3:2 unstable), 14,925 (1:2 stable) | Sec. VII.A |
| UMAP parameters | n_n = 100, m_dist = 0.0, n_c = 2 | Sec. VI.A step 4 |
| Concatenated torus revolutions | 3 at each end | Sec. VI.A step 5 |
| Newton steps / start for torus-point search | 10; Y = (0.01, 0.01) | Sec. VI.B |
| Acceptance | f <= 1e-12, ||F||_2 <= 1e-12 (about 400 m, 1 mm/s) | Sec. VI.B |
| Integrator tolerances | abs 1e-15, rel 1e-14, RKPD(8,9), GSL | Sec. VI.A |
| Time per transfer | about 15 s, i7-2600K 3.40 GHz | Sec. VII.A |
| 3:1 to 1:3 example | C_J = 3 (Fig. 10) | Sec. VII.C |
| 2:3 to 1:5 example | C_J = 2.6 (Fig. 11) | Sec. VII.C |
| Fig. 8 / 9 axes | max z_apo of departure torus (about 0.02-0.11) and arrival torus (about 0.05-0.32), black marker = feasible transfer | Figs. 8-9 |

Qualitative results usable as checks: transfers T1-T4 differ in geometry; T4 makes several
revolutions in the Earth vicinity and a final close lunar flyby before approaching the 1:2
torus and "deviates most significantly from the initial and final tori". For the T1-like
family, a transition from a given 3:2 torus exists only to some 1:2 tori and appears tied to the
difference in maximum out-of-plane component: nearly planar departure tori connect only to
nearly planar arrival tori, and the allowed difference grows along the families (Fig. 8).
Sample transfers A1 (20th torus of each family, largest out-of-plane motion), A2, A3 (second
torus of each family, almost planar), A4 (Fig. 8). The T4-like family (Fig. 9, samples
B1-B4) connects a wider range of torus pairs including large out-of-plane differences, which
the paper attributes "likely" to its multiple close flybys. The number of grid points in the
Fig. 8 and 9 scatter plots is not stated; the figure appears to show a roughly 20-by-20 grid
(not counted exactly), plus one marker at the lower-left corner that the text does not discuss.

## (5) Homoclinic transitions
**Not computed.** Every example is heteroclinic between tori of two different resonant
families (3:2 to 1:2, 3:1 to 1:3, 2:3 to 1:5). The word "homoclinic" does not appear. The
paper does say the method "may be used to compute natural transfers between various other
unstable 2-tori", without demonstrating a torus-to-itself connection.

## (6) Limitations and failure cases stated
- Visualising the 5-D crossings through 2-D or 3-D projections is unreliable: crossings close
  in a projection "may not be close in the full five-dimensional phase space"; a fourth
  dimension or extra constraints could help but would complicate analysis or "significantly
  shrink" the design space.
- The first six crossings of each manifold are unusable (the manifold lingers near its torus).
- Equality-constrained multiple shooting was numerically sensitive; the optimisation form is a
  workaround; reducing that sensitivity is "an ongoing effort ... in future work".
- The torus-point search works "if the torus is originally described by a sufficiently large
  number of invariant curves".
- Torus families were computed only to 20 members, the last at the edge of what the chosen
  discretisation supports.
- Continuation seeks one transfer per torus pair near the initial guess; it ends when no
  feasible transfer is found or the family ends; other transfers or other departure/arrival
  longitudes might connect more torus pairs (not explored).
- Transfers between tori are numerical approximations within the stated tolerance, not exact
  heteroclinic connections.
- The 3:1 to 1:3 and 2:3 to 1:5 examples are single point solutions; additional
  geometries "may be generated" from other regions of the embedding (not done).
- The mapping from UMAP-selected pairs to successful connections (hit rate, number of pairs
  tried) is not reported beyond four initial pairs in Fig. 4 and four in Fig. 7.

## What it does NOT contain
- No table; no initial conditions for any torus or transfer; no rotation numbers or
  fundamental frequencies; no Jacobi-constant values for individual transfers; no
  transfer times, time-of-flight, periods of the tori, or delta-V; no length/time/velocity
  unit values.
- No homoclinic (torus-to-itself) connection; no connection between a torus and a periodic
  orbit; no connections between more than two tori (no multi-leg chain).
- No time-periodic, elliptic, bicircular or ephemeris model; no higher-fidelity validation.
- No statement on cyclers, cycler orbits or repeated lunar encounters as a design goal (the
  transfers do include lunar flybys, e.g. T4, but no sequence or V-infinity is analysed).
- No comparison of UMAP against brute-force or other crossing-pair search (no hit-rate or
  speed comparison beyond the qualitative claim of fewer distance computations).
- No resonance-order, island or Poincare-section analysis of the tori beyond the examples.

## Bibliography (all 35 entries, as printed)
1. Carrico, J., Jr., Dichmann, D., Policastri, L., Carrico, J., III, Craychee, T., Ferreira, J., Intelisano, M., Lebois, R., Loucks, M., Schrift, T., and Sherman, R., "Lunar-Resonant Trajectory Design for the Interstellar Boundary Explorer (IBEX) Extended Mission," Advances in the Astronautical Sciences, Vol. 142, Dec. 2011, pp. 771-789.
2. Dichmann, D., Parker, J., Nickel, C., and Lutz, S., "Trajectory Design for the Transiting Exoplanet Survey Satellite," International Symposium on Spacecraft Flight Dynamics, Johns Hopkins Univ., Applied Physics Lab., May 2014.
3. Short, C., Howell, K. C., Haapala, A., and Dichmann, D., "Mode Analysis for Long-Term Behavior in a Resonant Earth-Moon Trajectory," Journal of the Astronautical Sciences, Vol. 64, No. 4, 2017, pp. 156-187. https://doi.org/10.1007/s40295-016-0098-9
4. Anderson, R. L., "Tour Design Using Resonant-Orbit Invariant Manifolds in Patched Circular Restricted Three-Body Problems," Journal of Guidance, Control, and Dynamics, Vol. 44, No. 1, Jan. 2021, pp. 106-119. https://doi.org/10.2514/1.G004999
5. Anderson, R. L., and Lo, M. W., "Role of Invariant Manifolds in Low-Thrust Trajectory Design," Journal of Guidance, Control, and Dynamics, Vol. 32, No. 6, Nov.-Dec. 2009, pp. 1921-1930. https://doi.org/10.2514/1.37516
6. Koon, W. S., Lo, M. W., Marsden, J. E., and Ross, S. D., "Resonance and Capture of Jupiter Comets," Celestial Mechanics and Dynamical Astronomy, Vol. 81, Sept. 2001, pp. 27-38. https://doi.org/10.1023/A:1013398801813
7. Barrabes, E., Mondelo, M. J., and Olle, M., "Numerical Continuation of Families of Heteroclinic Connections Between Periodic Orbits in a Hamiltonian System," Nonlinearit [sic], Vol. 26, No. 10, 2013, pp. 2747-2765. https://doi.org/10.1088/0951-7715/26/10/2747
8. Duarte, G., and Jorba, A., "Using Normal Forms to Study Oterma's Transition in the Planar RTBP," Discrete and Continuous Dynamical Systems-B, Vol. 28, No. 1, 2023, pp. 230-244. https://doi.org/10.3934/dcdsb.2022073
9. Lykawka, P. S., and Mukaia, T., "Resonance Sticking in the Scattered Disk," Icarus, Vol. 192, No. 1, 2007, pp. 238-247. https://doi.org/10.1016/j.icarus.2007.06.007
10. Haapala, A. F., and Howell, K. C., "Trajectory Design Strategies Applied to Temporary Comet Capture Including Poincare Maps and Invariant Manifolds," Celestial Mechanics and Dynamical Astronomy, Vol. 116, No. 3, 2013, pp. 299-323. https://doi.org/10.1007/s10569-013-9490-y
11. Calleja, R. C., Doedel, E. J., Humphries, A. R., Lemus, A., and Oldeman, B. E., "Boundary-Value Problem Formulations for Computing Invariant Manifolds and Connecting Orbits in the Circular Restricted Three Body Problem," Celestial Mechanics and Dynamical Astronomy, Vol. 114, No. 1, 2012, pp. 77-106. https://doi.org/10.1007/s10569-012-9434-y
12. Jorba, A., "Numerical Computation of the Normal Behavior of Invariant Curves of n-Dimensional Maps," Nonlinearity, Vol. 14, No. 5, 2001, pp. 943-976. https://doi.org/10.1088/0951-7715/14/5/303
13. Gomez, G., and Mondelo, J. M., "The Dynamics Around the Collinear Equilibrium Points of the RTBP," Physica D, Vol. 157, No. 4, 2001, pp. 283-321. https://doi.org/10.1016/S0167-2789(01)00312-8
14. Olikara, Z. P., and Scheeres, D. J., "Numerical Method for Computing Quasi-Periodic Orbits and Their Stability in the Restricted Three-Body Problem," IAA Conference on Dynamics and Control of Space Systems, AAS Paper 12-361, Univelt Inc., San Diego, CA, March 2012, pp. 911-930.
15. Olikara, Z., "Computation of Quasi-Periodic Tori and Heteroclinic Connections in Astrodynamics Using Collocation Techniques," Ph.D. Thesis, Ann and H. J. Smead Aerospace Engineering Sciences, Univ. of Colorado, Boulder, CO, 2016.
16. McCarthy, B., "Cislunar Trajectory Design Methodologies Incorporating Quasi-Periodic Structures with Applications," Ph.D. Thesis, School of Aeronautics & Astronautics, Purdue Univ., West Lafayette, IN, 2022.
17. Kumar, B., Anderson, R. L., and de la Llave, R., "Using GPUs and the Parameterization Method for Rapid Search and Refinement of Connections Between Tori in Periodically Perturbed Planar Circular Restricted 3-Body Problems," Proceedings of the 2021 AAS/AIAA Space Flight Mechanics Meeting, Vol. 176, Advances in the Astronautical Sciences Series, Univelt Inc., San Diego, CA, 2021, pp. 2311-2330.
18. Murphy, K. P., Probabilistic Machine Learning: An Introduction, MIT Press, Cambridge, MA, 2022, Chap. 20.
19. McInnes, L., Healy, J., and Melville, J., "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction," ArXiv e-prints, Feb. 2018. ArXiv: 1802.03426.
20. Cao, J., Spielmann, M., Qiu, X., Huang, X., Ibrahim, D. M., Hill, A. J., Zhang, F., Mundlos, S., Christiansen, L., Steemers, F. J., Trapnell, C., and Shendure, J., "The Single-Cell Transcriptional Landscape of Mammalian Organogenesis," Nature, Vol. 566, No. 7745, 2019, pp. 496-502. https://doi.org/10.1038/s41586-019-0969-x
21. Diaz-Papkovicha, A., Anderson-Trocme, L., Ben-Eghan, C., and Gravel, S., "UMAP Reveals Cryptic Population Structure and Phenotype Heterogeneity in Large Genomic Cohorts," PLos Genetics, Vol. 15, No. 11, 2019, Paper e1008432. https://doi.org/10.1371/journal.pgen.1008432
22. Bloch, T., Watt, C., Owens, M., McInnes, L, and Macneil, A. R., "Data-Driven Classification of Coronal Hole and Streamer Belt Solar Wind," Solar Physics, Vol. 295, No. 3, 2020, pp. 1-29. https://doi.org/10.1007/s11207-020-01609-z
23. Szebehely, V., Theory of Orbits: The Restricted Problem of Three Bodies, Academic Press, London, 1967, Chaps. 1, 5.
24. Vaquero, M., "Spacecraft Transfer Trajectory Design Exploiting Resonant Orbits in Multi-Body Environments," Ph.D. Thesis, School of Aeronautics & Astronautics, Purdue Univ., West Lafayette, IN, 2013.
25. Anderson, R. L., and Lo, M. W., "Dynamical Systems Analysis of Planetary Flybys and Approach: Planar Europa Orbiter," Journal of Guidance, Control, and Dynamics, Vol. 33, No. 6, 2010, pp. 1899-1912. https://doi.org/10.2514/1.45060
26. Barrabes, E., and Gomez, G., "Spatial p-q Resonant Orbits of the RTBP," Celestial Mechanics and Dynamical Astronomy, Vol. 84, No. 4, 2002, pp. 387-407. https://doi.org/10.1023/A:1021137127909
27. Koon, W., Lo, M., Marsden, J., and Ross, S., Dynamical Systems, the Three-Body Problem and Space Mission Design, Marsden Books, 1967 [sic], Chap. 4, 2011.
28. Verhulst, F., Nonlinear Differential Equations and Dynamical Systems, Springer-Verlag, Berlin, 1996, Chap. 5. https://doi.org/10.1007/978-3-642-61453-8
29. Haapala, A., "Trajectory Design in the Spatial Circular Restricted Three-Body Problem Exploiting Higher-Dimensional Poincare Maps," Ph.D. Thesis, School of Aeronautics & Astronautics, Purdue Univ., West Lafayette, IN, 2014.
30. Baresi, N., "Spacecraft Formation Flight on Quasi-Periodic Invariant Tori," Ph.D. Thesis, Ann and H. J. Smead Aerospace Engineering Sciences, Univ. of Colorado, Boulder, CO, 2017.
31. Perko, L., Differential Equations and Dynamical Systems, 3rd ed., Springer, New York, 2000, Chap. 3.
32. Gomez, G., Koon, W. S., Lo, M. W., Marsden, J. E., Masdemont, J., and Ross, S. D., "Connecting Orbits and Invariant Manifolds in the Spatial Restricted Three-Body Problem," Nonlinearity, Vol. 17, No. 5, 2004, pp. 1571-1606. https://doi.org/10.1088/0951-7715/17/5/002
33. McInnes, L., Healy, J., and Astels, S., "How UMAP Works," 2016, https://umap-learn.readthedocs.io/en/latest/how_umap_works.html [retrieved Feb. 2022].
34. Galassi, M., Davies, J., Theiler, J., Gough, B., Jungman, G., Alken, P., Booth, M., Rossi, F., and Ulerich, R., GNU Scientific Library Reference Manual, 3rd ed., 2021.
35. "MATLAB," MathWorks, Natick, MA, 2021.

(Diacritics were dropped when transcribing author names; typographical slips marked [sic] are
as printed.)

## Cross-references (facts from this repository's corpus index at digest time)
Of the cited works, the index already lists Olikara 2016 thesis [15], Gomez et al. 2004 [32],
Anderson and Lo 2010 [25], Barrabes and Gomez 2002 [26] and Vaquero 2013 [24]. The Kumar,
Anderson and de la Llave GPU-connections work [17] is represented by a 2025 SIAM J. Appl.
Dyn. Syst. arXiv version and two related 2021/2023 papers. Not found in the index by name:
Calleja et al. 2012 [11], Haapala and Howell 2013 [10], Haapala thesis [29], McCarthy thesis
[16], Baresi thesis [30], Short et al. 2017 [3], Jorba 2001 [12], Olikara-Scheeres 2012 [14],
Gomez-Mondelo 2001 [13] (searched by file-name patterns only; absence from the index is not a
claim about the paper corpus). Bonasera and Bosanac 2023 is itself cited by the Rosengren
et al. 2026 Astrodynamics Primer (see its digest) in a list of works on the Earth-Moon CR3BP
periodic-orbit atlas.
