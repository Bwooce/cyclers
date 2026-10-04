# Digest: Lantoine & Russell 2011, "Near ballistic halo-to-halo transfers between planetary moons"

Citation: G. Lantoine and R. P. Russell, "Near Ballistic Halo-to-Halo Transfers between Planetary Moons",
*Journal of the Astronautical Sciences* 58(3):335-363 (July-September 2011), DOI
`10.1007/BF03321174`; presented at the George H. Born Symposium, Boulder, May 2010. 29 PDF pages (journal
pp. 335-363). Filed in the private paper corpus as
`lantoine-russell-2011-near-ballistic-halo-to-halo-transfers-planetary-moons-jas-58-335-doi-10.1007-BF03321174.pdf`
(md5 `7f798c90194709871d122d6742ce49ad`). Text layer present but the tables and equations are partly
garbled; Tables 1 to 6 and Eqs. 7, 11, 13, 14 were read from page images. (The wanted-list description of this
paper as a flyby-derivatives paper was wrong; it is a moon-to-moon transfer design paper.)

Markers: PRINTED = from the paper (page and table given); COMPUTED = computed this session with project code or a
scratch script (nothing committed); DERIVED = worked out here; PROJECT = a statement about this repository.

## 1. What the paper does

Goal: a systematic method for fuel-optimal, near-ballistic transfers from a halo orbit of one planetary moon
to a halo orbit of another moon of the same planet, by combining resonant-hopping gravity assists with the
invariant manifolds of the halo orbits. Four stages of increasing fidelity (p. 337, abstract):

1. Initial guess from unstable resonant periodic orbits (K:L resonances of the planar CR3BP) and one invariant
   manifold leg falling off or onto the halo orbit.
2. Decomposition into two INDEPENDENT ideal CR3BPs (one per moon), each optimised by multiple shooting.
3. End-to-end refinement in a PATCHED three-body model (phasing of the two moons included).
4. Continuation (homotopy) to a four-body ephemeris model.

Benchmark: Ganymede L1 halo to Europa L2 halo (Jupiter system), delta-v 54.7 m/s and 204.5 days in the
ephemeris model (Table 7). The authors claim (p. 360) the first end-to-end near-ballistic transfer in any continuous
force model between loosely captured states of two different moons.

## 2. The models (equations)

**Ideal CR3BP phase (p. 348, PRINTED description; equations DERIVED in the standard form).** Each phase is a CR3BP in
which the moon of interest is on a prescribed circular coplanar orbit about the planet-moon barycentre. In
nondimensional rotating units (total mass 1, moon-planet distance 1 = 1 DU, angular rate 1; TU = 1 / n) the standard
equations, which the paper does not print, are

    x'' - 2 y' = Omega_x,  y'' + 2 x' = Omega_y,  z'' = Omega_z,
    Omega = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2,  C = 2 Omega - |v|^2   (Jacobi constant, DU^2/TU^2).

Integration is carried out in the NON-rotating frame (p. 348, p. 350): the initial state of each leg is rotated from the
rotating frame to the inertial frame, then integrated with Runge-Kutta 7(8) at tolerance 1e-12. The two phases
are integrated in opposite time senses (outer-moon phase forward, inner-moon phase backward, "forward-backward",
p. 348) so each is treated symmetrically.

**Patched three-body model (p. 338, p. 353).** The 4-body problem is decomposed into two CR3BPs, each under one moon. They are
linked by (a) a common patch point (planar orbit) chosen on the Tisserand-Poincare graph (Eq. 11, 12) and (b) a phase
between the two moons at the patch, Eq. 13. There is no simultaneous two-moon force anywhere in this model; the
moon not under consideration is absent.

**T-P graph patch point (PRINTED, Eq. 11, p. 352).** With a_M1, a_M2 the semimajor axes of the two moons, and (r_a, r_p) the
apoapsis and periapsis radii of the planar patch orbit about the planet:

    C_M1 = 2 a_M1 / (r_a + r_p) + 2 sqrt( 2 r_a r_p / ((r_a + r_p) a_M1) )
    C_M2 = 2 a_M2 / (r_a + r_p) + 2 sqrt( 2 r_a r_p / ((r_a + r_p) a_M2) )

with C_M1 the Jacobi constant of the forward trajectory and C_M2 that of the backward one; the solution (r_a*, r_p*) is the
intersection of the two Tisserand level sets. Each phase then has a one-dimensional final constraint g = r_p/a - r_p/a* (Eq. 12;
forward phase targets r_p*, backward phase targets r_a*). The Tisserand parameter is "almost equivalent" to the Jacobi
constant only when the Poincare section is far from the minor body (footnote 4, p. 351).

**Phasing (PRINTED, Eq. 13, p. 353).** Required phase of Europa with respect to Ganymede:

    theta = arccos( R*_Gan . R*_Eur / (r*_Gan r*_Eur) ) + w_Eur TOF,

with R* the position vectors of the two moons at the patch (no initial phasing), TOF the total time of flight and w_Eur Europa's
angular rate. Reading note: arccos returns an unsigned angle, so the sense of the phase offset is not fixed by the
printed formula; a direction check is needed in any use (DERIVED).

**Ephemeris model by homotopy (PRINTED, Eq. 14, p. 354).** X(lambda) = (1 - lambda) X_CR3BP + lambda X_ephem, lambda in [0, 1],
applied to the states of the planet and moons; stepped from 0 to 1 re-solving each time. The ephemeris is a self-made
"fake" one: Jupiter plus the two moons as point masses, integrated with an n-body propagator, initial conditions from two-body
motion (Table 6), centre-of-mass motion subtracted, cubic-spline interpolated; the epoch must satisfy Eq. 13 once
over the interval. The four-body equations (mutual point-mass gravity, standard) are not printed. Masses are not printed
beyond the mass ratios of Table 1 ("accurate to a few digits of the most recent estimates [49,50]").

**Optimisation (Eqs. 8, 9, 10; PRINTED).** Multiple shooting with nodes at each flyby; decision vector per leg
Z_i = [tau, X0_i, t0_i, tf_i, DeltaV_1..DeltaV_mi] (i = 1; tau the position on the halo) or without tau (i > 1);
minimise sum ||DeltaV_j,i|| subject to boundary constraint Psi_0(Z_1) = 0 (X0_1 on the halo, Eq. 10: X0,1 - X_h(tau) = 0, implemented with a
cubic spline of the halo orbit in tau; the manifold parameter epsilon is deliberately NOT in this constraint so the solver finds the
natural departure with a small first delta-v), continuity constraints g_i(Z_i, Z_i+1) = 0, and the final constraint (Eq. 12). Eight
delta-v impulses per inertial revolution are in the decision vector (314 variables in the Ganymede phase, 963 in the Europa phase,
Table 5). Tools: OPTIFOR (in-house) with SNOPT; STM by complex-step differentiation; constraints satisfied to a normalised 1e-8,
"around 10 m and 0.1 mm/s" (p. 355). Intel Fortran, 2.0 GHz processor.

## 3. The initial-guess construction

- Resonance notation K:L. Table 4 rows label 3:4, 9:7, 4:3, 11:8, 7:5 with a = (K/L)^(2/3) in DU (COMPUTED: (3/4)^(2/3) = 0.8255,
  (9/7)^(2/3) = 1.1824, (4/3)^(2/3) = 1.2114, (11/8)^(2/3) = 1.2365, (7/5)^(2/3) = 1.2515; the printed a column matches all five;
  Table 3 gives 4:5 -> 0.8618 and 6:5 -> 1.1292, also matching). The prose on p. 339 describes K as moon revolutions and L as
  spacecraft revolutions, which is the opposite of what these numbers imply; follow the numbers.
- Resonant path selection rules (p. 339): a monotonic change of semimajor axis consistent with the travel direction; the last
  resonance bounded by the Hohmann orbit (a_B L/K > a_H for the outer portion); limited flight time.
- Halo orbits are produced by continuation in the Jacobi constant from small planar Lyapunov orbits through the vertical bifurcation (p. 340).
- Manifold initial condition (Eqs. 1 to 3, p. 340-341): t_tau = t0 + tau T (tau in [0, 1] on the orbit),
  V_u(X(tau)) = Phi(tau, t0) V_u(X0), V_s likewise, X_u,s(X(tau)) = X(tau) +- epsilon V_u,s / ||V_u,s||.
  The sign picks the "interior" manifold for the outer moon and the "exterior" for the inner moon.
  (tau, epsilon) are two parameters for one manifold trajectory, not independent: many pairs give the same trajectory at a
  different time.
- Contour maps (Figs. 4, 5): the semimajor axis of the planetocentric orbit at the first crossing of the opposite x-axis, on a
  (log epsilon, tau) grid, for C = 3.0069 (Fig. 4) and 3.0059 (Fig. 5), Ganymede L1 halos. Contours are nearly straight lines
  because of Floquet theory (Eq. 4, V_u(X(tau)) = (rho_u)^tau p(tau), p periodic), with slope tan(alpha) = log(rho_u) (Eq. 5).
  PRINTED check (p. 344): for C = 3.0069 the slope read from the contour is alpha = 1.2741, the Floquet prediction alpha = 1.2696;
  the 6:7 contour of Fig. 5 reaches the same departure at (tau, eps) = (0, 1e-6) and about (1, 1e-3). The logarithm is base 10
  (DERIVED: a three-decade rise in epsilon per revolution matches tan(1.27) = 3.2 and the axis range of Fig. 9; the paper writes
  only "log"). alpha versus C is Fig. 8 (alpha about 1.22 to 1.28 for C from 3.0045 to 3.0070; plot only, no table).
- Empirical relation (Eqs. 6, 7, p. 345-346): xi = log(epsilon) sin(alpha) - tau cos(alpha);
  a = (a_max + a_min)/2 + ((a_max - a_min)/2) cos(omega xi + phi). Least-squares fit for C = 3.0069 (PRINTED): a_min = 9.1285e5 km,
  a_max = 9.5769e5 km, omega = 3.6763, phi = 3.9916. For a given a there are at most two xi: "type I" (da/dxi > 0) and "type II"
  (da/dxi < 0); Figs. 10, 11 show both for the 4:5 resonance (a = 0.8618 DU = 9.224e5 km, COMPUTED with LU = 1.070339e6 km).
  Not reproduced here: I could not recover the figure's I and II positions from Eq. 7 and the printed fit to better than my reading
  error of the plot, because the fit applies to C = 3.0069 and the transfer used C = 3.0066.
- The ten-step recipe for one portion is on pp. 346-348.

## 4. Printed numbers usable as sourced tests

All from the Ganymede-Europa benchmark. Units: DU = LU of Table 1 (moon-planet distance), TU = 1/n; the orbital period quoted in Table 1 is
2 pi TU (COMPUTED: 7.154280561 d / 2 pi = 1.1386 d per TU for Jupiter-Ganymede; the GM implied by 2 pi/P and LU, 1.26696e8 and
1.26690e8 km^3/s^2, equals Jupiter's GM to within the moon mass, which fixes the reading).

**Table 1 (p. 350): CR3BP parameters.**

| CR3BP | Mass ratio | Orbital radius LU (km) | Orbital period (days) = 2 pi TU |
|---|---|---|---|
| Jupiter-Ganymede | 7.8037e-5 | 1.070339e6 | 7.154280561 |
| Jupiter-Europa | 2.5280e-5 | 6.709e5 | 3.550439254 |

Registry comparison (COMPUTED from `core/satellites.py`: GM_Jupiter 1.26686534e8, Ganymede GM 9887.834, a 1070400 km; Europa GM 3202.739,
a 671100 km): mass ratio Ganymede 7.8044e-5 (paper 7.8037e-5, 9e-5 relative), Europa 2.5280e-5 (identical), radius Ganymede 1.0704e6 (paper
1.070339e6). So the registry reproduces the paper's Europa numbers and differs by 9e-5 in the Ganymede ratio.

**Table 2 (p. 352): halo orbits, rotating frame, y0 = 0, xdot0 = 0, zdot0 = 0.** (The printed column head "Y0 (DU/TU)" is the velocity ydot0.)

| | x0 (DU) | z0 (DU) | ydot0 (DU/TU) | Period (TU) | C |
|---|---|---|---|---|---|
| Halo 1 (Ganymede L1) | 0.9768297703815 | 0.0067575508137 | -0.0338705588350 | 3.01368039319 | 3.0066 |
| Halo 2 (Europa L2) | 1.0118043920085 | 0.0087547927135 | 0.0357066388227 | 3.06456025428 | 3.0024 |

**Table 3 (p. 352): manifold trajectories.** Manifold I + 4:5 resonance: tau = 0.5, epsilon = 1.7e-6, flight time 30.34682557231 TU, a = 0.8618 DU,
C = 3.0066 (Ganymede). Manifold II + 6:5 resonance: tau = 0.5, epsilon = 1.05e-6, flight time 42.88481483746 TU, a = 1.1292 DU, C = 3.0024 (Europa).

**Table 4 (p. 353): periodic resonant orbits, y0 = z0 = 0, xdot0 = zdot0 = 0.**

| Resonance | x0 (DU) | ydot0 (DU/TU) | Period (TU) | a (DU) | C |
|---|---|---|---|---|---|
| 3:4 | 0.9639250025000 | -0.037537693295765 | 19.1527202833 | 0.8255 | 3.0066 |
| 9:7 | 1.022912512093 | 0.035866768601460 | 56.8415853699 | 1.1824 | 3.0024 |
| 4:3 | 1.025860244947 | 0.038403138969070 | 25.3393083838 | 1.2114 | 3.0024 |
| 11:8 | 1.028261885259 | 0.040854642490345 | 69.268896450 | 1.2365 | 3.0024 |
| 7:5 | 1.029619747357 | 0.042564223924969 | 44.117093502 | 1.2515 | 3.0024 |

(Which mass ratio each row uses is not printed; the C column assigns 3.0066 to the Ganymede system and 3.0024 to the Europa system.)

COMPUTED test of all seven orbits with `core/cr3bp.py` (`jacobi_constant`, `propagate`, DOP853 at 1e-12), mass ratio 7.8037e-5 for the Ganymede rows and 2.5280e-5
for the Europa rows, one printed period, rows as printed:

| Orbit | C computed (printed) | closure of position / velocity after one printed period |
|---|---|---|
| Halo 1 | 3.006603 (3.0066) | 1.4e-5 / 3.9e-5 |
| Halo 2 | 3.002402 (3.0024) | 7.0e-6 / 2.0e-5 |
| 3:4 | 3.006598 (3.0066) | 1.6e-7 / 2.2e-7 |
| 9:7 | 3.002371 (3.0024) | 1.5e-6 / 1.2e-6 |
| 4:3 | 3.002354 (3.0024) | 4.7e-7 / 2.5e-7 |
| 11:8 | 3.002374 (3.0024) | 9.5e-7 / 3.2e-7 |
| 7:5 | 3.002379 (3.0024) | 5.9e-7 / 1.7e-7 |

The printed initial conditions and periods therefore close to 1e-5 or better (nondimensional) at the stated mass ratios. The residual is larger for the halos, which
are unstable and for which the 5-significant-digit mass ratio (relative uncertainty of order 1e-5 of the ratio) is the probable cause (not tested). The Jacobi
constants of the resonant orbits are 3.00235 to 3.00238, which the paper rounds to 3.0024; test them to 5e-5, not to the printed 4 decimals. These are published
orbits with printed numbers and are valid positive controls for `core/cr3bp.py` at Jupiter-Ganymede and Jupiter-Europa mass ratios (suggested tolerances: closure 1e-4 for the halos,
1e-5 for the resonant orbits; C to 5e-5).

**Table 5 (p. 354): the two portions.** Ganymede-dominant: C = 3.0066, targeted apse r_p* = 6.9466e5 km, 2 flybys, 314 variables, 52 constraints. Europa-dominant:
C = 3.0024, r_a* = 1.01763e6 km, 5 flybys, 963 variables, 107 constraints. Resonant paths (p. 356): Ganymede L1 halo (type I manifold) -> 4:5 -> 3:4 -> r_p*;
Europa (backward in time) halo (type II manifold) -> 6:5 -> 9:7 -> 4:3 -> 11:8 -> 7:5 -> r_a*. Seven flybys in all.

COMPUTED check of Eq. 11 with the printed Table 5 apses: C_M1 evaluated with a_M1 = 1.070339e6 km, (r_a, r_p) = (1.01763e6, 6.9466e5) km gives 3.0068 (printed 3.0066); C_M2 with
a_M2 = 6.709e5 km gives 3.00238 (printed 3.0024). So (r_a*, r_p*) reproduces the Europa constant to rounding and the Ganymede constant to 2e-4 (a position of r_a* uncertain by about
3000 km would account for it; not further investigated). A usable independent test of Eq. 11.

**Table 6 (p. 354): inertial initial conditions of the ephemeris model (km, km/s).**
Jupiter position [-13.2225, 10.6217, 0], velocity [-0.2176e-3, -0.2708e-3, 0]. Ganymede position [1.07026e6, 0, 0], velocity [0, 10.8790, 0].
Europa position [5.2303e5, -4.2015e5, 0], velocity [8.6057, 10.7129, 0].
COMPUTED consistency: circular speed sqrt(GM_J / r) = 10.880 km/s at 1.07026e6 km (printed 10.8790); Europa's radius from the printed position is 6.709e5 km and its
speed 13.74 km/s, matching sqrt(GM_J / 6.709e5 km) = 13.742 km/s. Europa starts at longitude -38.8 degrees (atan2 of the printed components).

**Table 7 (p. 355): results.**

| Model | delta-v | TOF | runs | computational time | function calls |
|---|---|---|---|---|---|
| Independent CR3BPs | 40.5 m/s | 204.4 d | 2 | about 50 min | 2200 |
| Patched CR3BP | 42.2 m/s | 204.3 d | (blank in the print) | about 10 min | 370 |
| Ephemeris four-body | 54.7 m/s | 204.5 d | 1000 | about 1 week | about 100,000 |

The 1000 runs are the lambda continuation in steps of 1e-2 ("a small variation of 1e-2 for lambda is necessary to ensure convergence", p. 358); the lowest delta-v
is the independent-phases one, which is not continuous at the patch. Three main impulses, two near the first and last flybys (p. 359, Figs. 21, 23; maximum
impulse axis about 14 m/s on Fig. 21, read from a plot). The authors compare ~50 m/s with a navigation allowance of about 5 m/s per flyby times 7 flybys, ~35 m/s (p. 358).
Epoch: not printed as a date (the phasing of Eq. 13 sets it). Flyby altitudes are not printed. Figures 13 to 23 are plots without data tables.

## 5. Reconciliation with project code

- PROJECT: `core/cr3bp.py` (planar and spatial CR3BP, STM, Jacobi) is the model of the independent phases; the Jovian K:L resonant orbits are built by
  `search/jovian_resonant_families.py` (Jupiter-Europa, targets Anderson & Lo 2011 Table 1) and `search/jovian_resonant_connections.py`; Halo continuation by Jacobi
  constant exists as `search/halo_family_at_jacobi.py` and `search/cr3bp_3d_family_tracer.py`; the Tisserand machinery is `search/tisserand.py`. Not found by grep: any T-P
  patch point solver for two different moons (Eq. 11 and 12), a halo-manifold semimajor-axis contour map, the phasing Eq. 13, or the lambda blend of Eq. 14.
- PROJECT: the project's two-moon models are NOT the paper's patched model. `core/ccr4bp.py` (Europa base, Ganymede as a simultaneous perturber on a concentric circle,
  no moon-moon force), `core/crnbp.py` (N = 5), `search/two_moon_periodic_890.py` (planar circular Titania-Oberon four-body, both moons acting at once) and
  `search/titania_oberon_realeph_895.py` (URA111 kernel, J2 and J4, unsoftened moons, the Sun) all include both moons at all times. The paper's patched model has one moon at a time
  and the phasing as a patch condition; the paper's own four-body ephemeris is a plain Jupiter-plus-two-moons n-body integration, which is also different from the CCR4BP (moons are
  massive and interact; the concentric-circle idealisation is not used).
- PROJECT: `nbody/jovian.py` and `nbody/jovian_ideal.py` give Galilean-moon states from the JUP365 kernel; the paper's "fake" ephemeris is not reproducible from the printed data (no
  masses, no epoch), only its initial conditions (Table 6) are.
- PROJECT: constants. Registry reproduces Table 1 for Europa and to 9e-5 for the Ganymede mass ratio (above). The unit TU is 1/n of each CR3BP (the paper's Table 1 period column is 2 pi TU).
- Related held papers: Kumar, Anderson, de la Llave 2023 (Ganymede-Europa transfers, `2026-07-27-727-kumar-anderson-delallave-2023-ganymede-europa-transfers-digest.md`), Anderson & Lo 2011
  (resonant flybys), Campagnola & Russell 2010 (T-P graph, cited as [17]), Bradley & Russell 2014 (`2026-10-04-digest-bradley-russell-2014-patched-conics-to-full-gravity-continuation.md`),
  whose mass-scaling continuation is a different homotopy from Eq. 14. Lantoine, Russell & Campagnola 2011 (Acta Astronautica 68, the resonant-orbit-boundary predecessor [30]) is not held.

## 6. Techniques applicable to the project's problems

- `#890` and `#895` (the Titania-Oberon orbit; the transfer between two moons' regions is the core problem). The Uranian pair is in the same regime as the paper's
  Jovian pair: mass ratios (COMPUTED from the registry: GM_Uranus 5.7945564e6, Titania 226.9, Oberon 205.3) 3.92e-5 and 3.54e-5, against 7.80e-5 and 2.53e-5 for Ganymede and Europa.
  Directly usable: (a) the T-P patch point, Eq. 11, with a_M1 = 436,298 km (Titania) and a_M2 = 583,511 km (Oberon), as the construction of the connecting planar orbit and of its
  periapsis and apoapsis, to compare with `#888`'s patched-conic enumeration; (b) the resonant-path pruning rules of p. 339 (monotonic semimajor axis, last resonance bounded by the
  Hohmann orbit, here a_H = (436,298 + 583,511)/2 = 509,904 km as a first estimate, DERIVED); (c) the phasing relation Eq. 13, to be checked against the `#890` convention that Titania
  and Oberon are on opposite sides at the middle encounter (2.5 synodic periods per leg); (d) the epsilon-tau contour method and the Floquet slope relation (Eqs. 4 to 7) for halo
  orbits at Titania or Oberon. What does not transfer: the orbit in `#890` is PERIODIC in a two-moon four-body model with hyperbolic flybys of both moons; the paper
  is an open one-way transfer with halo boundaries, nowhere periodic, with a patched (one moon at a time) model, so it is a method source and a cost benchmark (delta-v and
  convergence cost at the same mass ratios), not a result for the periodic orbit. In `#895` (real-ephemeris arcs) the paper's lesson is that the step from the idealised model to the
  ephemeris was "far from trivial" (about 100,000 function calls and a week for 1000 continuation steps at lambda spacing 1e-2); size the `#895` continuation budget accordingly and
  compare lambda blending with the Bradley-Russell mass scaling already used.
- `#893` (two-moon model audit: identity test plus published control). The paper supplies a published control with printed numbers for `core/cr3bp.py` at two mass ratios (section 4, seven
  orbits closing to 1e-5 or better) and a T-P control (Eq. 11 reproduces C_M2 = 3.0024 from the printed apses). It supplies NO control for the project's simultaneous two-moon models
  (`ccr4bp`, `crnbp`, `two_moon_periodic_890`): its four-body ephemeris has no printed masses or epoch, so the 54.7 m/s and 204.5 d cannot be reproduced. State that plainly in the audit.
  The model difference (patched versus simultaneous) is also a useful negative: the patched-versus-ephemeris delta-v differs by 30 percent (42.2 vs 54.7 m/s), a stated size for the error of dropping the other moon.
- `#914` (Russell & Strange moon cyclers with integrated flybys). The multi-stage ramp (independent three-body, patched, ephemeris), delta-v impulses at every node (8 per revolution) and
  optimiser structure are applicable, and the lambda blend of Eq. 14 is an alternative to a mass-scaling homotopy; the paper cites Russell & Strange 2009 [48] as the origin of the blend.
  No cycler content (not periodic), no numbers to reproduce.
- `#929` (integrator controls through one flyby). Printed integrator facts: RK7(8) at 1e-12, nondimensional constraint tolerance 1e-8 (position 10 m, velocity 0.1 mm/s), complex-step STM.
  The paper's flybys are at Ganymede and Europa with unspecified altitudes. A flyby-through reversibility control in the project's Jovian system could use the seven-flyby patched solution only if
  the solution were published; it is not (plots only).
- `#913` to `#918`. `#915` (optimiser as published): the staged ramp and the sensitivity to the final ephemeris stage parallel the eccentricity, inclination, ephemeris ramp; no numbers shared.
  `#916` and `#922` (invariant curves and tori): the paper contributes the manifold-contour idea and the Floquet slope; nothing about invariant tori. `#913`, `#917`, `#918` (heliocentric rows,
  elliptic-problem rows, global search): not applicable.
- Halo and resonant orbit seeds. Table 2 and Table 4 can serve as seeds for `search/halo_family_at_jacobi.py` and the Jovian resonant corrector at the two printed Jacobi constants, with the sourced closure test of section 4.

## 7. Follow-ups (task numbers not registered here)

1. Add the seven printed orbits (Tables 2 and 4) as a published-control test for `core/cr3bp.py` at both mass ratios (closure and Jacobi constant tolerances in section 4).
2. Add a T-P patch-point solver (Eq. 11, 12) with the Table 5 apses as the sourced test, then apply it to Titania-Oberon with the registry constants and compare with the `#888` enumeration.
3. Reproduce Fig. 8 and the epsilon-tau contour of Fig. 4 for a Ganymede L1 halo near C = 3.0069 (alpha = 1.2741 measured, 1.2696 predicted) and check the sinusoidal fit parameters of Eq. 7; then do it for Titania.
4. Resolve the sign of the phasing in Eq. 13 by a direct construction before using it; compare with the `#890` convention.
5. Decide what the project's simultaneous two-moon models owe the patched model: a head-to-head of one transfer in both models, to put a number on the 30 percent patched-versus-ephemeris gap.
6. Find the Lantoine, Russell & Campagnola 2011 resonant-hopping paper (Acta Astronautica 68, doi 10.1016/j.actaastro.2010.09.021) and Campagnola & Russell 2010 (T-P graph, J. Guid. Control Dyn. 33:476) in the corpus index; acquire if missing.
7. If `#895` continuation cost matters, record the paper's cost (1000 runs, one week, 1e5 function calls at lambda step 1e-2) as the only published comparison.

## 8. Reading quality

Tables 1 to 7 and Eqs. 7, 11, 13, 14 read from images; the text layer is garbled for tables and Figures 3 to 11 (no data to extract; figure coordinates are not tabulated). Digits I could not read with
certainty: none in the numbers quoted. The prose definition of K and L (p. 339) contradicts the semimajor axes printed in Tables 3 and 4 (section 3).
