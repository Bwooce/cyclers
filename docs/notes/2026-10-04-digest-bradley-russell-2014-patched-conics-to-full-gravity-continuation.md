# Digest: Bradley & Russell (2014), "A Continuation Method for Converting Trajectories from Patched Conics to Full Gravity Models"

The Journal of the Astronautical Sciences 61:227-254 (2014), DOI 10.1007/s40295-014-0017-x, 28 pages, published online
31 October 2014. Authors: Nicholas Bradley and Ryan P. Russell (University of Texas at Austin).
Filed in the private paper corpus as
bradley-russell-2014-continuation-method-converting-trajectories-patched-conics-to-full-gravity-models-jas-61-227-doi-10.1007-s40295-014-0017-x.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page, section, equation or figure) or INFERRED
(our reading, or a comparison with project work). Page numbers are PDF pages 1 to 28; the journal page is the PDF page
plus 226 (INFERRED from the printed range 227-254; the PDF carries no page numbers on its pages). The paper has no
numbered tables; its numbers are in the text and in figures, and values read off figures are marked "read from plot,
approximate". Context: task `#890` (Titania-Oberon-Titania periodic orbit by continuation in the moons' mass), whose
literature check flagged this paper as possibly already publishing the continuation step.

## 0. What the paper is

READ (abstract, p1): "A method is introduced to transition space trajectories from low fidelity patched conics models to
full-ephemeris n-body dynamics. The algorithm incorporates a continuation method that progressively re-converges solution
trajectories in systems with incremental changes in the dynamics. Continuation is accomplished through the variation of a
control parameter, which is tied to body ephemeris locations and masses, minimum flyby altitudes, and sphere of influence
sizes. The intermediate models provide a continuous and differentiable path between solutions in the simplified and
n-body dynamics, and the boundary values of the control parameter replicate the patched conics and full-ephemeris models
exactly. Each successive step preserves the qualitative properties of the initial guess by successively altering flyby
states and body masses."

It is a trajectory-design paper about open ballistic tours (a heliocentric Earth-Venus-Venus-Earth-Jupiter transfer and a
seven-encounter Jovian moon tour, plus a twelve-encounter chained extension). It treats no periodic orbit, no cycler and
no restricted four-body problem (sections 2 and 3 below).

## 1. The method, step by step

### 1.1 Starting model and final model

- Start: the zero-sphere-of-influence (ZSOI) patched-conic model (READ, p1, p3, p4): each target body "is a point mass that
  instantaneously turns the trajectory", and "the effects of the central body and third bodies are nonexistent for the
  'instantaneous' flyby" (p11). The initial trajectory is built with Lambert arcs between encounters (READ, p9: "the
  trajectory is re-calculated between each encounter ... using Lambert targeting"). The paper says the method "is
  applicable to any simplified starting model; the ZSOI model is chosen because it is a common preliminary model for space
  trajectory designers" (p3), and the abstract adds that "a similar approach ... may be taken with any simplified starting
  guess, such as a restricted three-body model".
- End: "a spacecraft trajectory in an n-body dynamical ephemeris (e.g. JPL's SPICE ephemerides) with dynamics governed by
  bodies of actual mass and size (as opposed to a massless or point-mass assumption)" (READ, p4). Equations of motion about
  the central body, eq. 9 (p12): `r_ddot = -mu_CB r/r^3 - sum_i mu_TB,i ( r_TB,i/r_TB,i^3 + (r - r_TB,i)/|r - r_TB,i|^3 )`.
  Equations about target body i, eq. 10 (p12), include the central body and the other target bodies as perturbers.
  In the examples the bodies are the Sun, Venus, Earth and Jupiter (EVVEJ) and Jupiter, Callisto, Ganymede and Europa
  (moon tour) (READ, p17, p19). Whether further bodies (other planets, other moons, J2) were in the final model is not
  stated; the text of the examples names only the encounter bodies. INFERRED: the final model contains the central body and
  the encounter bodies only.

### 1.2 The continuation parameter

READ (p3, p11): a single scalar `kappa` in [0, 1], incremented in fixed steps `kappa_{i+1} = kappa_i + Delta kappa` (p13).
"A control parameter kappa in [0, 1] is defined that controls the dynamics and ephemeris model used, with kappa = 0
representing a purely ZSOI Keplerian model and kappa = 1 representing an n-body ephemeris dynamical model" (p11).
Quote (p11): "The control parameter kappa is tied to three important aspects of the model: body ephemerides, body radii, and
body masses (which relate to sphere of influence sizes)."

What is scaled with kappa:

| Quantity | Rule | Where |
| --- | --- | --- |
| Body ephemerides | `x_fake(kappa) = (1 - kappa) x_kepler + kappa x_real` (linear blend of a mean-Keplerian ephemeris and the real one) | eq. 1, p5; p11 |
| Target-body GM | `mu_TB = kappa mu_TB,final` | eq. 8, p11 |
| Target-body radius | `rho_TB = kappa rho_TB,final` | eq. 8, p11 |
| Central-body GM, when moving about a target body | `mu_CB = kappa mu_CB,final` | eq. 8, p11 |
| Central-body GM, when moving about the central body | `mu_CB = mu_CB,final` (held fixed) | eq. 8, p11 |
| Sphere-of-influence radius | `r_SOI = r_AB (kappa mu_B / mu_A)^(2/5)` (the GM carries kappa) | eq. 12, p14 |
| Close-approach radius at the terminal bodies | user-defined, "scaled linearly by kappa at each continuation step" | p15 |

Quote (p11): "Around each target body (TB), the mass of all bodies is varied by kappa. Around the central body (CB), the mass of
the central body remains constant, while the masses of the target bodies are tied to kappa." Equation 8 as printed:
`mu_CB = kappa mu_CB,final` "if about target body", `mu_CB,final` "if about central body"; `mu_TB = kappa mu_TB,final`;
`rho_TB = kappa rho_TB,final`.

INFERRED from eq. 8 and eq. 9: the central-body GM is not reduced to compensate the target-body masses, so the total mass
of the system is not held fixed as kappa varies.

Footnote 3 (p11) and the text above eq. 1 (p5): the "fake" ephemeris is a simple linear blend; the footnote notes "a more
robust way ... would be to calculate a new auxiliary ephemeris for every value of kappa (for example, by numerically
propagating body states according to the dynamics at the current value of kappa)".

### 1.3 The auxiliary Keplerian ephemeris

READ (pp5-10): the ZSOI guess is first re-converged on a mean Keplerian ephemeris of the bodies (a two-body ephemeris
about the central body fitted to the real one) so that the initial guess and the target bodies obey the same dynamics.
Mean motion `n_B = Theta/(t_f - t_0)` from the accumulated angle (eq. 2, Algorithm 2, p8); `a_B = (mu_A/n_B^2)^(1/3)`
(eq. 3); the other elements are time means of the osculating ones (eq. 4); the true anomaly at t_0 by eq. 5; then a
batch least-squares correction of the initial Cartesian state against the real ephemeris with the two-body STM (eq. 6,
p9). Fig. 4 (p10): the Keplerian-minus-real position difference for Ganymede and Callisto oscillates, "but there is no
secular growth", with Ganymede up to about 2300 km and Callisto up to about 400 km over 400 days (read from plot,
approximate). The re-convergence of the ZSOI trajectory on this ephemeris uses a Newton-Raphson on the n-2 intermediate
encounter times to match interior V-infinity magnitudes, with Lambert arcs between encounters, tolerance 1 m/s (eq. 7,
p10).

### 1.4 What is held fixed at each step, and the flyby re-seeding

READ, introduction (pp2-3): "Between each step, the flyby radii, body masses, and sphere of influence sizes are
artificially altered such that the turning angles and flyby properties from the previous converged solution are
preserved. This preservation is illustrated schematically in Fig. 1."

Fig. 1 caption (p3): "Schematic of the continuation method preserving flyby geometry and turning angle. To maintain a
constant turn angle for any value of kappa: r_p = kappa r_p,f". The figure shows kappa = 0, m = 0 (a sharp corner, ZSOI),
kappa = 0.1, 0.2, 0.3 with m = 0.1 m_f, 0.2 m_f, 0.3 m_f, and kappa = 1, m = m_f with the periapsis radius r_p,f inside
the sphere of influence.

READ (p13): "At each step in the continuation algorithm, the converged states from the previous kappa value are kept (other
than the flyby periapsis states, which are altered to preserve flyby turning angles), and are assumed to represent a
reasonable first guess to the feasible solution at the new value of the control parameter."

READ (p19, on Fig. 9, the second Venus flyby, Venus-centred frame): "As kappa is increased, the algorithm increases the
flyby periapsis between successive iterations to maintain the same hyperbolic turning angle, and the radius of Venus is
artificially increased in a linear fashion."

Construction of the flyby state from the Lambert solution (Appendix, pp24-26, eqs. 16-25), READ:

- Turn angle and eccentricity from the two V-infinity vectors: `delta = arccos(v_inf^- . v_inf^+ / |v_inf|^2)`,
  `e = 1/sin(delta/2)` (eq. 16, p24).
- Periapsis direction `r_p_hat = (v_inf^- - v_inf^+)/|v_inf^- - v_inf^+|` (eq. 17); orbit normal
  `h_hat = (v_inf^- x v_inf^+)/|...|` (eq. 18); periapsis velocity direction `v_p_hat = h_hat x r_p_hat` (eq. 19).
- `r_p = mu (e - 1)/|v_inf|^2`, `v_p = sqrt(|v_inf|^2 + 2 mu/r_p)` (eq. 20); `r_p = r_p r_p_hat`, `v_p = v_p v_p_hat` (eq. 21).
- Terminal encounters (one V-infinity only): angular-momentum direction from the target body's own orbit plane (eq. 22),
  `e = 1 + r_p |v_inf|^2/mu_B`, `delta = 2 arcsin(1/e)` (eq. 23), rotate v_inf by delta/2 about h_hat, periapsis state by
  eq. 25.

Note that in eq. 20 the periapsis radius is proportional to mu at fixed e and V-infinity, which is why `r_p = kappa r_p,f`
preserves the turn when mu scales with kappa (INFERRED, direct algebra on eq. 20 and the Fig. 1 caption).

Which flyby quantities are the unknowns (eq. 15, p15), READ: for each non-terminal encounter n the six-vector
`x_p,n = [ r_p,n/rho_body,n, alpha_p,n, beta_p,n, v_inf,n, alpha_vp,n, beta_vp,n ]` (periapsis radius over body radius,
right ascension and declination of the periapsis position, V-infinity magnitude, right ascension and declination of the
periapsis velocity). Quote: "This specific formulation of the state vector is chosen because it means the parameters do
not vary significantly between steps of the continuation procedure, regardless of the value of Delta kappa." "The first
quantity in Eq. 15 represents the ratio of the flyby periapsis radius to the body radius at each step of the continuation
method. This quantity remains nearly constant through the continuation process, since both quantities in the ratio change
nearly linearly with the continuation parameter kappa." The periapsis speed is computed from two-body dynamics from the
V-infinity value.

### 1.5 The inner problem: unknowns and constraints

READ (pp14-16): a multiple-shooting problem over the encounter chain.

- Patch points: at the sphere-of-influence crossing before and after each flyby, `t_i^- = t_p,i - Delta_SOI`,
  `t_i^+ = t_p,i + Delta_SOI` (eq. 13, p14), where `Delta_SOI` is the fixed time of flight from periapsis to the SOI,
  recomputed for each kappa because the SOI radius changes (p14: "The patch times are not held constant between continuation
  steps"). An intermediate state is placed halfway in time between consecutive encounters, `t_inter,i = (t_i + t_{i+1})/2`
  (eq. 14), "for n encounters, there exist n - 1 intermediate times, each with six defined state elements" (p15).
- Free parameters: n flyby periapsis times; 6(n-2) full periapsis state elements; 6(n-1) intermediate state elements; 10
  terminal-body state elements. Total "13n - 8" (p16). The encounter times are free and Delta_SOI is fixed per flyby.
- Equality constraints: position and velocity continuity at the patch points, weighted to match "within 1 km and 1 m/s";
  12 per non-terminal encounter, 12 at the terminal bodies, and 2 terminal periapsis constraints (`r_p . v_p = 0`), total
  12n - 10 equalities (p16).
- Inequality constraints: time increases monotonically through each event, and periapsis altitudes of the intermediate
  flybys above a user-prescribed limit (body radii from eq. 8), total 5n - 5 (p16). The limit's value is not printed.
- Solver: SNOPT (SQP) "with no objective function to obtain a feasible solution" (p16), gradients by complex-step
  differentiation (pp16-17). If a problem is infeasible SNOPT switches to "elastic mode", which "allows the continuation
  to proceed through regions of the homotopy where the problem may be infeasible"; "no cases were encountered where the
  problem became infeasible and subsequently returned to a feasible solution at kappa = 1" (p16). The paper notes the
  problem "may easily be cast as an optimization problem rather than a feasibility problem" (p16).
- Dynamics inside the inner loop is the two-body or n-body propagation of eqs. 9-10 with the kappa-scaled masses, with the
  ZSOI-converged state as the kappa = 0 solution. The paper does not say which integrator it used for the n-body
  propagation (not read on any page).

### 1.6 Continuation strategy

READ (pp12-13): a convex homotopy `H(x_p, kappa) = kappa G(x_p, kappa) + (1 - kappa) F(x_p, kappa)` (eq. 11); a constant
increment Delta kappa ("the most basic approach"); the algorithm is "the most basic embedding algorithm". Predictor-corrector
or pseudo-arclength methods are mentioned as improvements ("a relatively simple improvement ... could include implementing
a PC method ... and/or dynamic correction of the step size", p13) and were not used. Algorithm 1 (p4): initialise
`kappa_0 < 1`, set `Delta kappa`, loop `kappa = [kappa_0 : Delta kappa : 1]`: update the flyby states from kappa, converge
with the inner loop. Quote (p13): "Large values of Delta kappa may lead to a solution in a different family of
trajectories, or may prevent the algorithm from converging at all." And: for very simple problems "it may be possible to set
Delta kappa to 1".

## 2. Worked examples

There are no tables in the paper. All numbers printed in the text for the two examples, plus the extension, are below.

### 2.1 Earth-Venus-Venus-Earth-Jupiter transfer (p17, Figs. 6-8)

READ: initial ZSOI trajectory from EXPLORE, "a patched conics tour generating program external to the current work"
(p17). Encounter order: Earth departure, Venus flyby 1, Venus flyby 2, Earth flyby, Jupiter arrival (labels in Figs. 6, 7).
Dates are not printed in the text. "After tuning the continuation parameters, the values kappa_0 = 0.05 and Delta kappa =
0.05 were found to converge." "The specified initial and final radii about Earth and Jupiter are r_p,f = 10,000 km and
r_p,f = 100,000 km, respectively." Total CPU time "332 seconds" (Fortran 90, gfortran, MacBook Pro with a 2.53 GHz Core i5
and 4 GB of RAM, footnote 4). The number of continuation steps is not printed; kappa from 0.05 to 1 in steps of 0.05 is
20 values (INFERRED). Excess speeds, flyby altitudes and final discontinuities are not tabulated; Fig. 8 (p18) plots them
against kappa. Read from the plots, approximate: V-infinity of about 4 to 15 km/s across the five encounters; the Earth flyby
curve is the highest (about 14 to 15 km/s); the Jupiter arrival flight date moves by about 9 days (earlier) over the
continuation; the other encounter dates move by about 1 day or less; `r_p/r_body` between about 1.1 and 2.5 across flybys.
Final discontinuities: not printed; the constraint tolerance is 1 km and 1 m/s (p16).

READ (p19): "Other than the encounter time, the parameters are nearly constant through the entire continuation process ...
The fact that most parameters are constant suggests that the algorithm may converge directly when kappa_0 = 1, although the
variation of the flyby time suggests the necessity of the continuation method."

### 2.2 Seven-encounter tour of three Jovian moons (pp19-21, Figs. 10-13)

READ (p19): "The spacecraft follows an encounter order of: Callisto - Ganymede - Ganymede - Callisto - Ganymede - Europa -
Ganymede (CGGCGEG). This trajectory is complex, involves a large number of intermediate flybys with multiple revolutions about
the central body, and is in a region of fast-changing non-Keplerian dynamics. The result shown is converged using an
initial kappa_0 = 0.02 and Delta kappa = 0.02, with both the initial and final radii at r_p,f = 3000 km. The continuation
process completed in 22 minutes." Encounters are numbered from 0: C0, G1, G2, C3, G4, E5, G6 (Fig. 12, p21). Epoch
`t_0` is JD 2461287.866 (Fig. 12 caption, p21). Number of steps: kappa from 0.02 to 1 in steps of 0.02 is 50 values
(INFERRED; not printed). The tour is not periodic: it begins at Callisto and ends at Ganymede, with Jupiter as central
body (Figs. 10, 11); the spacecraft ends in a different state from the start.

Read from Fig. 12 (p21), approximate: the distance from Jupiter oscillates between about 10 and 29 Jupiter radii; encounter
times (days from `t_0`) at about 0 (C0), 29 (G1), 56 (G2), 63 (C3), 80 (G4), 85 (E5), 86 (G6); total span about 86 days.
Read from Fig. 13 (p21), approximate: V-infinity between about 1.8 and 4.2 km/s; `r_p/r_body` between about 1.1 and 3.9; the
flyby dates deviate from the ZSOI guess by about 0.25 day at most (C0) over the continuation.

READ (p20): "Figure 13 shows the variation of the model parameters through the continuation process. Unlike the EVVEJ
variation in Fig. 8, the model parameters do not remain nearly constant. This behavior signifies the necessity of a
continuation procedure; the scenario does not converge directly for kappa_0 = 1. It is evident that the converged value
for some parameters varies significantly between a near-ZSOI case (kappa << 1) and the full ephemeris model (kappa = 1).
Because of this variation through the continuation, the initial ZSOI guess for the parameters is not a sufficiently good
initial guess for the full model, and using these parameters directly as an initial guess leads to lack of convergence."

### 2.3 Twelve-encounter extension by patching (pp22-23, Figs. 14-16)

READ (p22): "In the Jovian moon tour solution with seven encounters, the flyby state at Europa (E5) may be used as the
initial state in a subsequent tour, and solutions of arbitrary length may be constructed by chaining smaller sequences
together." "An example patched tour is shown in Figs. 15 and 16. This extension continues the CGGCGEG tour by patching an
EGGEGGE tour, effectively adding five more encounters. In the procedure, a small Delta v (<= 1 m/s) is allowed at the patch
time immediately after the patching body (Europa, in this case), while periapsis state and flyby time are held constant."
The labels in Fig. 15 are E5, G6, G7, E8, G9, G10, E11 (12 encounters, numbered 0 to 11). Quote (p23): "Current design
practice is to chain shorter (one or two flyby) solutions together to form longer tours; the algorithm presented here
extends the capability of tour design to include more encounters in a single calculation sequence." Numerical details for
the extension (kappa schedule, run time) are not printed.

## 3. Limits and failure modes

READ (p20, p23):

- "The ability to find a continuous solution may be limited by the problem itself, as a given ZSOI trajectory may not exist
  ballistically in a full ephemeris model without violating constraints. The existence of a solution in a full model is
  dependent on a variety of factors, including the number of encounters and how closely the full ephemeris dynamics are
  modeled as being Keplerian. The ability to find these solutions is sensitive to inner-loop software parameters such as
  convergence tolerance, maximum parameter step size, and the continuation step size Delta kappa. Some cases with a high
  number of flybys proceed with acceptable convergence properties until larger values of kappa are encountered, at which
  point convergence fails and no solution is found. In many of these cases, the solution does not exist with positive flyby
  altitudes, which may indicate that these many-encounter ZSOI tour solutions simply do not admit a continuous ballistic
  solution."
- In homotopy terms the failure "may potentially indicate the nonexistence of a smooth curve connecting the two homotopy
  levels kappa = 0 and kappa = 1 ... Likewise, multiple solutions may exist for the kappa = 1 homotopy level; smooth curves
  may branch during the continuation, and multiple 'zero points' ... are certainly admissible. In addition, the curve may
  simply terminate with no existing solution for some values of kappa, including kappa = 1." (p21)
- Remedies suggested (p21): "including extra margins in the equivalent flyby altitudes and reducing the number of encounters
  in the ZSOI initial guess"; intermediate delta-v manoeuvres on each leg, or casting the problem as a total delta-v
  minimisation.
- Conclusions (pp23-24): "the maximum number of encounters for a Jovian moon tour is typically between five and nine,
  depending on the scenario ... The method is not always guaranteed to find a feasible solution in the n-body ephemeris model.
  For a relatively small number of encounters (two to four), the method reliably finds a feasible solution in the full model.
  For more encounters, especially in more complex dynamical regimes, feasible solutions may simply not exist for the n-body
  ephemeris model. However, the algorithm is capable of finding feasible tours consisting of a relatively high number of
  encounters (five to nine)."
- Step size (p13): "Large values of Delta kappa may lead to a solution in a different family of trajectories, or may
  prevent the algorithm from converging at all."
- Minimum flyby altitude: only that a user-prescribed lower limit is an inequality constraint (p16); its values are not printed
  for either example. The terminal radii are 10,000 km (Earth), 100,000 km (Jupiter) and 3000 km (the Jovian moon tour).
- Low excess speed, long-duration flybys, resonant legs, multi-revolution legs: no dedicated discussion. The moon tour is
  described as "multiple revolutions about the central body" and "fast-changing non-Keplerian dynamics" (p19) and converged;
  the paper does not test resonant (repeating) legs as such, and does not discuss low-V-infinity flybys, where the sphere of
  influence of a moon is a small fraction of its Hill sphere. INFERRED from absence: nothing in the 28 pages addresses them.

## 4. Comparison with the project's `#890` continuation

Project side (READ in `docs/notes/2026-10-04-890-titania-oberon-candidate.md` sections 1.1 to 1.4 and
`data/OUTSTANDING.md` `#890`): planar concentric circular restricted four-body model (Uranus, Titania, Oberon); mass scale
`lam` multiplies both moons' GM and the planet's GM is reduced so the total stays the registry value; continuation from
`lam` = 1e-3 (not zero) to 1; guess at each flyby = the moon-centric hyperbola that turns the patched-conic incoming
V-infinity onto the outgoing one, scaled to `lam`; symmetric multiple shooting with perpendicular crossings at the two
flybys, solving a square 4M+4 system by Newton with the variational equations; the target is a PERIODIC orbit (period five
synodic periods) in the model, not a trajectory in an ephemeris.

| Item | Bradley-Russell 2014 | Project `#890` | Same or different |
| --- | --- | --- | --- |
| Continuation parameter | One scalar `kappa` in [0, 1] scaling target-body GM, body radius, SOI radius, close-approach radius and the blend of the ephemeris (eqs. 1, 8, 12; READ) | One scalar `lam` scaling both moons' GM (and nothing else); moons stay on circles (READ, project note 1.1) | SAME idea for the mass (a scalar on the flyby bodies' GM from a small value to the physical one); DIFFERENT in that the paper also scales radii, SOI and ephemeris, and does not reduce the central body's GM to hold the total fixed (INFERRED from eq. 8) |
| Start of the continuation | `kappa = 0`, the ZSOI solution is the limit; in practice `kappa_0` = 0.05 (EVVEJ) and 0.02 (moon tour) (READ, pp17, 19) | `lam` = 1e-3, because the patched conic is the singular limit and not an exact solution at 0 (READ, project note 1.3) | SAME in effect: both start at a small positive mass fraction and not at zero. The paper does not discuss the singular limit; it solves at `kappa_0` directly from the ZSOI-seeded guess (INFERRED) |
| Seeding of each flyby | Periapsis state from the patched-conic V-infinity pair by the hyperbola construction, eqs. 16-21 (periapsis direction along `v_inf^- - v_inf^+`, velocity along `h x r_p`, `r_p = mu (e-1)/v_inf^2`); between steps the periapsis state is "altered to preserve flyby turning angles", `r_p = kappa r_p,f` (READ, pp3, 13, 24-25) | Hyperbola of the patched-conic turn, periapsis along `v_in_hat - v_out_hat`, velocity along `v_in_hat + v_out_hat`, scaled to `lam` (READ, project note 1.4) | SAME construction: this is already published in the appendix and Fig. 1. The geometry in the project (periapsis along the difference of the V-infinity unit vectors) is the same as eq. 17 (INFERRED by comparing the two descriptions; the project's exact code was not re-read) |
| Target model | N-body ephemeris (SPICE) with real moon motion; dynamics via eqs. 9-10 (READ, pp4, 12) | Planar circular restricted four-body model, moons on circular orbits, no ephemeris (READ, project note 1.1) | DIFFERENT |
| What is imposed at the end | Position and velocity continuity at patch points along an OPEN tour; free initial and final epochs; terminal periapsis radii; no periodicity (READ, p16) | Periodicity of the full orbit and a mirror (time-reversal) symmetry with perpendicular crossings at both flybys; a fixed point of the stroboscopic map (READ, project note 1.2, 1.4) | DIFFERENT: no symmetry, periodicity, monodromy or stability analysis anywhere in this paper (READ by absence over 28 pages) |
| Inner solver | SNOPT feasibility NLP with inequality constraints and elastic mode; complex-step gradients (READ, p16) | Newton on a square system with analytic variational equations and step backtracking; fold detection by the sign of det J (READ, project note 1.4) | DIFFERENT |
| Step control | Fixed `Delta kappa` (0.05 and 0.02); no predictor-corrector, no arclength; no fold handling beyond elastic mode (READ, pp13, 16) | 1e-3 to 1 with step halving on failure to a minimum step of 1e-4, and fold detection (READ, project note 1.4) | DIFFERENT |
| Central-body mass treatment | Held fixed about the central body; scaled by `kappa` in the target-centred equations (READ, eq. 8) | Planet GM reduced so total GM is fixed at every `lam` (READ, project note 1.1) | DIFFERENT (INFERRED that the paper's total mass is not conserved) |

Plain statement of what is already published here (INFERRED from the table above, with the READ items cited):

1. Already published in this paper: continuation in a scalar that multiplies the flyby bodies' masses from near zero to the
   physical value, starting from a patched-conic (ZSOI) solution, re-seeding each flyby with the patched-conic-turn
   hyperbola whose periapsis scales with the mass (`r_p = kappa r_p,f`) so the turning angle is preserved, solving a
   multiple-shooting problem at each step, applied to a Jovian moon tour of seven encounters (and a chained twelve), and the
   starting from a small positive value in place of zero. The project cannot present the mass continuation of the
   patched-conic flyby chain, with hyperbola re-seeding, as its own method.
2. Not in this paper: the restricted four-body (planet plus two moons) circular model as the target; imposing periodicity
   and a mirror symmetry (a periodic orbit rather than an open tour); a flyby chain between two moons that closes on itself;
   the fold and Floquet analysis; the model-versus-ephemeris distinction the project draws (a periodic orbit of a model that
   is not a trajectory of the real system); Uranus.
3. The honest framing for `#890`: the continuation step is Bradley-Russell's method specialised to a different target (a
   periodic orbit of a restricted four-body model); the new content, if any, is the periodic closure in that model, not the
   continuation in the mass.

## 5. References the paper cites that bear on cyclers, moon tours and continuation in mass

Transcribed as printed (paper reference list, pp26-28; numbers are the paper's).

Continuation applied to cyclers and Earth-Moon families, the pair the introduction cites (p2, "computing Earth-Moon
transfer families and cycler trajectories [10, 39]"):

- [10] Casoliva, J., Mondelo, J.M., Villac, B.F., Mease, K.D., Barrabes, E., Olle, M.: Two classes of cycler trajectories in
  the earth-moon system. J. Guid. Control Dyn. 33(5), 1623-1640. doi:10.2514/1.46856
- [39] Yagasaki, K.: Computation of low energy Earth-to-Moon transfers with moderate flight time. Phys. D Nonlinear Phenom.
  197(3-4), 313-331. doi:10.1016/j.physd.2004.07.005

What the paper claims about them: only that continuation methods "have been applied" to that problem (p2). No result, method
detail or comparison is given. On p13 it groups [5, 10, 11, 28, 39] as works that use continuation "to determine families of
trajectories", and [10, 11, 17, 23] as users of predictor-corrector schemes. INFERRED from the titles and that grouping
(neither paper was read): these are continuation in an orbit-family parameter, not a continuation in the flyby body's mass.

Other works the paper cites for cyclers or patched-conic-to-full-model continuation (READ, p2 and the list):

- [30] Russell, R.P., Ocampo, C.A.: Optimization of a broad class of ephemeris Earth-Mars cyclers. J. Guid. Control Dyn.
  29(2), 354-367 (2006). doi:10.2514/1.13652 (cited for "calculating n-body patched conics trajectories from a basic circular
  coplanar model")
- [31] Russell, R.P., Strange, N.J.: Cycler trajectories in planetary moon systems. J. Guid. Control Dyn. 32(1), 143-157
  (2009). doi:10.2514/1.36610 (cited for "transitioning between ideal patched conics to full ephemeris patched conics")
- [22] Lantoine, G., Russell, R.P.: Near ballistic halo-to-halo transfers between planetary moons. J. Astronaut. Sci. 58(3),
  335-363 (2011). (continuation "from low-fidelity three-body to four-body ephemerides"; the paper says its ephemeris
  interpolation "is similar to the methodology used by Lantoine and Russell for intermoon halo-to-halo transfers", p11)

Moon-tour and patched-conic tour design (cited as "ballistic patched conics tour design [9, 13, 18, 19, 32, 33]", p2):

- [9] Campagnola, S., Buffington, B.B., Petropoulos, A.E.: Jovian tour design for orbiter and lander missions to Europa. Acta
  Astronautica 100, 68-81 (2014). doi:10.1016/j.actaastro.2014.02.005
- [13] Englander, J.A., Conway, B.A., Williams, T.: Automated interplanetary mission planning. AIAA/AAS Astrodynamics
  Specialist Conference, August, pp. 4517-4537 (2012)
- [18] Heaton, A.F., Strange, N.J., Longuski, J.M., Bonfiglio, E.P.: Automated design of the Europa orbiter tour. J.
  Spacecr. Rocket. 39(1), 17-22 (2002)
- [19] Izzo, D., Simoes, L.F., Martens, M., de Croon, G.C., Heritier, A., Yam, C.H.: Search for a grand tour of the Jupiter
  Galilean moons. Proceedings of the 15th Annual Conference on Genetic and Evolutionary Computation, pp. 1301-1308. ACM Press,
  New York (2013). doi:10.1145/2463372.2463524
- [32] Spreen, C., Mueterthies, M., Kloster, K.W., Longuski, J.M.: Preliminary analysis of ballistic trajectories to uranus
  using gravity-assists from venus, earth, mars, jupiter, and saturn. Proceedings of the AAS/AIAA Astrodynamics Specialist
  Conference, pp. 3411-3428 (2011) (an interplanetary Uranus-arrival paper, not a Uranian-moon tour)
- [33] Strange, N.J., Longuski, J.M.: Graphical method for gravity-assist trajectory design. J. Spacecr. Rocket. 39(1), 9-16
  (2002). doi:10.2514/2.3800

Continuation in mass or in a model parameter, and other continuation references:

- [29] Picot, G.: Shooting and numerical continuation methods for computing time-minimal and energy-minimal trajectories in the
  Earth-Moon system using low propulsion. Discrete Contin. Dyn. Syst. Ser. B 17(1), 245-269 (2011). doi:10.3934/dcdsb.2012.17.245
- [8] Caillau, J.B., Daoud, B., Gergaud, J.: Minimum fuel control of the planar circular restricted three-body problem.
  Celest. Mech. Dyn. Astron. 114(1-2), 137-150 (2012). doi:10.1007/s10569-012-9443-x
- [5] Biggs, J.D., McInnes, C.: Solar sail formation flying for deep-space remote sensing. J. Spacecr. Rocket. 46(3), 670-678
  (2009). doi:10.2514/1.42404
- [28] Paffenroth, R.C., Doedel, E., Dichmann, D.J.: Continuation of Periodic orbits around lagrange points and AUTO2000. In:
  Advances in the Astronautical Sciences, pp. 41-60 (2001)
- [11] Chow, C.C., Villac, B.F.: Mapping autonomous constellation design spaces using numerical continuation. J. Guid. Control
  Dyn. 35(5), 1426-1434 (2012). doi:10.2514/1.56434
- [16] Griesemer, P.R., Ocampo, C., Cooley, D.S.: Targeting ballistic lunar capture trajectories using periodic orbits. J.
  Guid. Control Dyn. 34(3), 893-902 (2011). doi:10.2514/1.46843 (a three-body solution used as a guess for a multi-body
  one)
- [23] Lara, M., Pelaez, J.: On the numerical continuation of periodic orbits. Astron. Astrophys. 389, 692-701 (2002).
  doi:10.1051/0004-6361
- [1] Allgower, E.L., Georg, K.: Numerical Continuation Methods: An Introduction. Springer, New York (1990); [2] Arkowitz, M.:
  Introduction to Homotopy Theory. Springer, New York (2011); [21] Keller, H.B.: Numerical solution of bifurcation and
  nonlinear equation problems. In: Rabinowitz, P.H. (ed.) Applications of Bifurcation Theory, p. 359 (1977); [34]
  Sundararajan, P., Noah, S.T.: Dynamics of forced nonlinear systems using shooting / arc-length continuation method -
  application to rotor systems. J. Vib. Acoust. 119, 9-20 (1997)
- [6] Buffington, B., Strange, N.: Patched-integrated gravity-assist trajectory design. AIAA/AAS Astrodynamics Specialist
  Conference (2007); [24] Marchand, B.G., Howell, K.C., Wilson, R.S.: Improved corrections process for constrained trajectory
  design in the n-body problem. J. Spacecr. Rocket. 44(4), 884-897 (2007). doi:10.2514/1.27205 (alternatives the paper
  lists for the model-fidelity transition, p2)

Entries [10] and [39] are printed in the paper's list without a year; volume and pages are as printed.
