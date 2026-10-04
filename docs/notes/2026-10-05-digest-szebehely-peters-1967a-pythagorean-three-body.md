# Digest: Szebehely and Peters 1967, "Complete solution of a general problem of three bodies" (the Pythagorean or Burrau problem)

V. Szebehely and C. F. Peters, "Complete solution of a general problem of three bodies", The Astronomical Journal 72(7):876-883 (September 1967), DOI 10.1086/110355 (Yale University Observatory; received 16 May 1967). Filed in the private paper corpus as `szebehely-peters-1967a-complete-solution-general-problem-three-bodies-aj-72-876-doi-10.1086-110355.pdf` (8 pages, ADS-type image scan, journal pages 876 to 883; PDF page n is journal page 875 + n). All pages, both tables and every printed equation were read on page images; Table I was re-read at 220 dpi. In the sibling paper (the 1967b digest below) this one is called "Paper 2".
Evidence tags: READ (journal page) is what the printed page says; COMPUTED is my own integration of 2026-10-05 (scratch programs, not committed; description in section 5); INFERRED is my reading across sources.
Companions: `docs/notes/2026-10-05-digest-szebehely-peters-1967b-periodic-pythagorean.md`, `docs/notes/2026-10-04-digest-peters-1968-numerical-regularization.md`, `docs/notes/2026-10-04-digest-aarseth-zare-1974-regularization-three-body.md`.

## 0. What the paper is

The first complete numerical solution of Burrau's 1913 problem: three bodies of masses 3, 4, 5 released at rest at the apices of a 3-4-5 right triangle, with the opposite sides of lengths 3, 4, 5. Burrau stopped at t = 3.35. Abstract (p.876): the solution "is neither quasi-periodic nor periodic but it assumes the form known in the recent Soviet literature as 'elliptic-hyperbolic'. In the final configuration two of the three participating bodies form a permanent binary while the third body is rejected to infinity. A new method of treating close approaches which allows achieving the solution is also described." The work was done in parallel at Yale (E. M. Standish), the Institute for Space Studies (R. Spinelli, with M. Lecar and V. Szebehely) and ETH Zurich (L. Stanek, under E. Stiefel) (p.876).

## 1. The problem and the units (READ pp.876-877)

- Masses m1 = 3, m2 = 4, m3 = 5; G = 1 (T^2 G M / L^3 = 1). Positions at t = 0 (Fig. 1): P1 = (1, 3) with m1 = 3, P2 = (-2, -1) with m2 = 4, P3 = (1, -1) with m3 = 5; all initial velocities zero, so the motion is planar, the angular momentum is zero, and the centre of mass is at the origin (COMPUTED: the mass-weighted position is exactly zero).
- Energy: "E = V = 769/60 = -12.8166..." (p.876, printed as +769/60 with the sign carried by E = V = -12.8166). COMPUTED: V = -(3*4/5 + 3*5/4 + 4*5/3) = -12.81667 = -769/60 exactly.
- Dimensional examples (p.877): in cgs units the unit of time is 3872 s, about 1.08 h; for three stars of 3, 4, 5 solar masses at 3, 4, 5 pc the unit of time is about 1.43e7 yr and L = 1 pc = 3e18 cm; the first 30 time units are "the first 430 million years" (p.880; COMPUTED: 30 x 1.43e7 = 4.29e8).
- Equations of motion, dimensionless, d^2 r_i/dt^2 = sum_{j != i} m_j (r_j - r_i)/|r_j - r_i|^3 (p.877).

## 2. The motion (READ pp.877-880; Figs. 2 to 8 are plots, no numbers)

Description printed, in the paper's time units:
- The initial motion is toward the centre of mass. At t1 = 1.879 the second and third bodies approach each other, r23 about 1e-2; the third (heaviest) body "performs repeated close approaches to the first and second bodies but these latter bodies do not find themselves near to each other": the heaviest body acts as the agent between the two lighter ones (p.877). There is an instant (Fig. 7) when m2 and m3 are closer to each other than to the first body (r12 > r23) and "at this critical time are already in the process of forming a close binary system that will be their final configuration" (p.877-878).
- At t = 3.35 (the end of Burrau's integration) m1 moves away from the origin while m2 and m3 are on approach trajectories (p.878).
- Close approaches m2-m3 for 0 to 10: 1 < t1 < 2, 3 < t2 < 4, 8 < t3 < 9, r23(t1) about 1e-2, r23(t2) about 6e-2, r23(t3) about 8e-3; "the last close approach occurs at t3 = 8.760 and it is the smallest distance during the interval of time shown in this figure". m1-m3 approaches at t4 about 3, 6 < t5 < 7, t6 about 10 (p.878). For 10 to 20: m1-m3 at t7 between 14 and 15 and at t8 about 17; m2-m3 at 11 < t9 < 12, t10 = 15.8299236, 19 < t11 < 20. "The close approach at t10 is r23(t10) about 4e-4 which is the smallest distance occurring between any bodies at any time during the evolution of this dynamic system between t = 0 and t = infinity." The relative velocity of m2 and m3 at that close approach is "approximately 171 in our units" while the velocity of m1 is v1 = 0.04 (p.878-879). In physical terms (the 3, 4, 5 pc example) this is 1.7e4 solar radii, 80 a.u., and a maximum velocity about 10 km/s; the first body's 0.04 is about 10 km/h in the cgs example (p.879).
- Near periodicity (p.879): a collision of m2 and m3 together with m1 at rest would reverse the solution (periodic); half the period would be t10 = 15.83 = T/2 if periodic; at t12 = 2 t10 = 31.66 = T the three bodies are close to their initial positions with small velocities (Fig. 5); but because the initial conditions are not exactly reproduced "the subsequent motion in these two figures is entirely different" (sensitivity). The existence of a periodic orbit near the Pythagorean conditions "is still being investigated" (answered by the 1967b paper).
- Escape (pp.879-881): at t13 = 47 the smallest-mass body m1 seems to depart in the third quadrant with xdot < 0, ydot < 0, while m2 and m3 form a binary moving away in the opposite direction; but at t14 = 53 m1 reverses (almost zero velocity). Fig. 7 shows the final binary formation and the penetration of the binary by m1: at t15 = 54 m1 approaches the binary with xdot > 0, ydot > 0; at t16 = 59.2 the binary is at periastron and at t17 = 59.4 nearly at apastron, with m1 about on the line joining m2 and m3; then m1 departs fast while m2 and m3 reach periastron again at t18 = 59.7. At t19 = 60 the speed of m1 is about 2.7 and its distance from the centre of mass about 2; treating m1 and the binary (mass m2 + m3) as an artificial two-body system with speed 0.9, the energy is positive: hyperbolic. Fig. 8 (t = 60 to 70): "The binary and the body with mass m1 depart with hyperbolic velocities. The period of the binary is approximately 0.9 time unit and therefore this binary performs 15 revolutions during one complete rotation of the galaxy." Alexeev's (1961) hyperbolic-elliptic conditions are satisfied at t = 69; the integration was continued to t20 = 102 (1.5e9 yr) with no change from Fig. 8. The final configuration: the escaper is body 1 (the lightest, m1 = 3); the binary is bodies 2 and 3.

## 3. Printed tables

Table I (p.881), "Close approaches between t = 0 and t = 30" (time; approximate distance; pair; the printed order is by listing, not by time; entries as read at 220 dpi):
| time | approx. distance | pair |
|---|---|---|
| t1 = 1.879 | 1e-2 | m2, m3 |
| t4 = 3.026 | 0.6 | m1, m3 |
| t2 = 3.801 | 6e-2 | m2, m3 |
| t5 = 6.898 | 0.1 | m1, m3 |
| t3 = 8.760 | 8e-3 | m2, m3 |
| t6 = 9.962 | 0.5 | m1, m3 |
| t9 = 11.611 | 0.2 | m2, m3 |
| t7 = 14.618 | 0.2 | m1, m3 |
| t10 = 15.830 | 4e-4 | m2, m3 |
| t8 = 17.001 | 0.3 | m1, m3 |
| t11 = 19.807 | 0.2 | m2, m3 |
| 21.791 | 0.4 | m1, m3 |
| 22.966 | 2e-2 | m2, m3 |
| 24.537 | 0.1 | m1, m3 |
| 27.780 | 5e-2 | m2, m3 |
| 28.679 | 0.5 | m1, m3 |
| 29.802 | 3e-3 | m2, m3 |
There are 17 rows; the 14 rows with t at most 25 are the "14 close approaches" of the Table II range.

Table II (p.883), "Comparison of different systems", range 0 to 25 where 14 close approaches occur:
| System | steps | time on computer (s) | change in total energy | tolerance |
|---|---|---|---|---|
| R (rectangular) | 17631 | 395 | 7e-7 | 1e-12 |
| V (only time regularised, dtau = V dt, V the potential energy) | 4300 | 280 | 8e-7 | 1e-12 |
| L (Levi-Civita) | 3450 | 155 | 2e-11 | 1e-10 |
| L | 600 | 26 | 7e-8 | 1e-6 |
Text (p.883): "The numbers listed under the heading of tolerance represent the accuracy attained in the first few integration steps"; in R and V the accuracy degenerates at each close encounter and faster than in L; L integrates through collisions; L was used for the majority of the time (when the particles are less than a distance of 2 apart).

## 4. Method (READ pp.880-883)

- Regularisation, section "Method of regularization" (pp.880-882): Sundman's time transformation alone is judged unsatisfactory for computation because terms x_i/r (both tending to zero) appear; the restricted-problem experience (Szebehely 1967, Astron. J. 72:370 and PNAS 57) is that regularising the independent variable alone increases the complexity of the equations. The method used for the general planar problem is Levi-Civita's: new time tau* = integral of dt/r + tau0* and x + i y = (xi + i eta)^2, applied to the relative coordinates x, y of one colliding body with respect to the other. Only isolated binary collisions are regularised: below a limiting separation r_ij of a pair its equations are regularised; if all three separations are below the limit the closest pair is chosen; when the pair separates the equations are transferred back to the original variables, "this procedure requires that regularization be switched on and off quite often", and every change of variables restarts the integration with new initial conditions (p.882).
- Integration: fifth-order Runge-Kutta of Zonneveld (1964, Math. Tracts 8), whose last (fifth-order) term estimates the truncation error; the step is chosen by Ollongren's (1966, private communication) formula h_{i+1} = (E_max/(E_max + E_5) + 0.45) h_i, with the step rejected and repeated if E_5 is not below E_max; about one step in 1000 was rejected (p.882).
- Checks (pp.882-883): total energy was monitored throughout, plus angular momentum and the centre-of-mass integral "at selected times"; "the requirement of keeping the total energy constant proved to be the most sensitive control", but energy is "a necessary but not a sufficient condition for accuracy". In the regularised system the transformed Hamiltonian, identically zero, is the energy control. Reversibility: integrate forward from t0 to t1 and back; "no numerical errors are allowed for"; reversals were executed at t = 32 and t = 62; the initial positions were reconstructed with errors in the tenth and third decimals respectively, with the energy error below 1e-10 in both; the reversal at t = 62 needed about 19000 steps (p.883).
- Machines: 1967 Yale Computer Center (IBM 7094, double precision, 16 digits, per Peters 1968).

## 5. Quick check of the printed numbers (COMPUTED, 2026-10-05; scratch scripts, no test files)

Method: two independent integrators, both DOP853 (scipy) with rtol 1e-12 to 1e-13: (a) Cartesian with the time transformation dt = g ds, g = r12 r13 r23/(r12 r13 + r12 r23 + r13 r23) (regular at binary collisions to the extent needed here), to t = 30 (4897 steps at 1e-12; relative energy error 1.2e-10 at 1e-12 and 2.0e-10 at 1e-13); (b) a Levi-Civita-regularised pair (2,3) with Aarseth's second-order perturbed-KS equations (see the Aarseth 1971 digest) with the binding energy integrated by h' = 2 u'.L^T(F_k - F_l), to t = 16 (935 steps; relative energy error 7.8e-13). Local minima of each pair distance were refined by a bracketed minimisation. The two tolerances (1e-12, 1e-13) agree on all minima to the printed precision (to about 1e-8 in time).
Residuals (computed minus printed; printed distances are one-figure approximations):
| printed t | computed t | computed minus printed | printed distance | computed distance | pair |
|---|---|---|---|---|---|
| 1.879 | 1.879343 | +3e-4 | 1e-2 | 9.700e-3 | 23 |
| 3.026 | 3.024213 | -1.8e-3 | 0.6 | 0.5714 | 13 |
| 3.801 | 3.800505 | -5e-4 | 6e-2 | 6.114e-2 | 23 |
| 6.898 | 6.897696 | -3e-4 | 0.1 | 0.1033 | 13 |
| 8.760 | 8.759755 | -2.5e-4 | 8e-3 | 8.516e-3 | 23 |
| 9.962 | 9.962810 | +8e-4 | 0.5 | 0.4617 | 13 |
| 11.611 | 11.611864 | +8.6e-4 | 0.2 | 0.1656 | 23 |
| 14.618 | 14.617499 | -5e-4 | 0.2 | 0.2263 | 13 |
| 15.830 (text: 15.8299236) | 15.829920 (LC: 15.8299203) | -3.6e-6 against the 8-digit value | 4e-4 | 4.138e-4 | 23 |
| 17.001 | 17.000993 | -7e-6 | 0.3 | 0.2619 | 13 |
| 19.807 | 19.806907 | -9e-5 | 0.2 | 0.2079 | 23 |
| 21.791 | 21.789934 | -1.1e-3 | 0.4 | 0.4188 | 13 |
| 22.966 | 22.965822 | -1.8e-4 | 2e-2 | 1.754e-2 | 23 |
| 24.537 | 24.536809 | -1.9e-4 | 0.1 | 0.1155 | 13 |
| 27.780 | 27.779632 | -3.7e-4 | 5e-2 | 5.002e-2 | 23 |
| 28.679 | 28.677960 | -1.0e-3 | 0.5 | 0.5387 | 13 |
| 29.802 | 29.801522 | -4.8e-4 | 3e-3 | 2.794e-3 | 23 |
All 17 printed events are found, each to within 2e-3 in time and to the printed precision in distance, and no extra close approach (a distance below 0.6) occurs in 0 to 30 that the table omits (the pair 12 distance never falls below 1.79, at t = 2.94). The one 8-digit value, t10 = 15.8299236, differs from my 15.8299203 by 3.3e-6, which I cannot attribute: the energy conservation of my run is better than 1e-12 to t = 16, and two tolerances and two formulations agree to 1e-8; a 1967 fifth-order Runge-Kutta run through several close approaches may simply carry that error (INFERRED).
Other printed statements checked: the relative speed at the t10 approach: computed |v2 - v3| = 208.55 at the minimum distance (this equals sqrt(2 x 9/r23) as expected for a nearly parabolic pass), against the printed "approximately 171": the printed value is not reproduced (it would correspond to r23 about 6.2e-4; INFERRED that it was read at a slightly different instant or from a different quantity); v1 = 0.0404 against the printed 0.04 (reproduced). Position at t = 31.66: m1 (1.032, 2.716), m2 (-2.117, -1.188), m3 (1.074, -0.679) against the initial (1, 3), (-2, -1), (1, -1): "approximately their initial positions" is reproduced. At t = 60: the speed of m1 is 2.730 (printed about 2.7) and its distance from the centre of mass 2.08 (printed about 2). At t = 70: m1 is at distance 21.4 (Fig. 8 axis labelled to 20, consistent) with speed 1.76; the binary has semi-major axis 0.553 and period 0.858 (printed "approximately 0.9"). The escaper is body 1. The energy of the artificial two-body system is positive. Sensitivity (the printed "sensitivity of this dynamical system to initial conditions"): between the two tolerances the positions at t = 62 differ by about 1e-2 and the binary semi-major axis by 3e-3, but at t = 32 only by 1e-8, so a late-time comparison (after t = 50) is not a usable test, while the first 30 time units are.

## 6. Comparison with Peters 1968 and Aarseth and Zare 1974

- Peters 1968, Table 1 (digest `docs/notes/2026-10-04-digest-peters-1968-numerical-regularization.md`): the same four rows of steps, seconds and energy (R 17631, 395, 7e-7; V 4300, 280, 8e-7; L 3450, 155, 2e-11; L 600, 26, 7e-8), "to t = 25", with "fourteen close approaches"; this paper's Table II is the same table with the tolerance column added (1e-12, 1e-12, 1e-10, 1e-6). Peters's text that regularised L is "about 15 times" faster than R compares R's 395 s with L's 26 s (the 600-step row); this paper does not state a factor. COMPUTED against Table I above: 14 close approaches have t at most 25. So Peters 1968 Table 1 is not a new calculation but this table, reprinted.
- Aarseth and Zare 1974, Table I Example I (digest `docs/notes/2026-10-04-digest-aarseth-zare-1974-regularization-three-body.md`): the same initial conditions (masses 3, 4, 5 at (1, 3), (-2, -1), (1, -1) at rest), integrated to t = 16. Their text: "the second and third particles experience four encounters within a distance of 0.01, whereas the closest approach between the first and third particles is about 0.1", and the smallest distance between the first and second is about two length units. COMPUTED: r23 is below 0.01 at three times in 0 to 16 (t = 1.879 with 9.70e-3, t = 8.760 with 8.52e-3, t = 15.830 with 4.1e-4), and a fourth r23 minimum (t = 3.800, 6.1e-2) is not within 0.01; so "four within 0.01" is not reproduced (the printed Table I here lists 6e-2 for that approach too). The other two statements are reproduced: closest m1-m3 distance 0.1033 (at t = 6.898) and smallest r12 1.80 (at t = 2.94). Their Table III uses a different problem (three unit masses). Their t_f = 16 contains the t10 = 15.83 deep approach (4e-4) that the 1967 paper calls the smallest distance for all time.
- Aarseth 1971 cites this paper as "Szebehely and Peters 1967, Astron. J. 72, 876" for "a critical three-body case" that "may also be reproduced reasonably well by the ordinary method" (digest: `docs/notes/2026-10-04-digest-aarseth-1971-direct-integration-n-body.md`, section 8 settles the citation).

## 7. Printed numbers usable as sourced tests

All are sourced (READ); the computed columns above are independent confirmations.
1. Initial state (exact): masses 3, 4, 5; positions (1, 3), (-2, -1), (1, -1); zero velocities; centre of mass at the origin; E = -769/60 = -12.81666...
2. Table I: 17 close-approach times to 3 decimals with the pair and the one-figure distance, usable to tolerance 2e-3 in time and a factor of 1.3 in distance. Two stronger values: t10 = 15.8299236 (text) with r23(t10) about 4e-4; r23(t3 = 8.760) about 8e-3.
3. Qualitative sequence: the pair alternates (m2,m3) and (m1,m3) with m3 in every approach; no (m1,m2) approach below 1.79; the smallest distance for all time is r23 about 4e-4 at t10 with v1 = 0.04; the final state is body 1 escaping hyperbolically with the (2,3) binary of period about 0.9.
4. Table II as the published efficiency comparison for the three formulations (R, time-only V, Levi-Civita L): ratios 5.1 and 29 in steps (R against the two L rows), 4.1 in steps for the 25-body case in Peters's paper.
5. Reversibility: reversal at t = 32 recovers the initial positions to 1e-10 and at t = 62 to 1e-3 with energy error below 1e-10 (about 19000 steps at 62). COMPUTED at 62: forward 10164 steps (DOP853, time-transformed, 1e-12), consistent with that order.
6. Not reproduced: the relative speed 171 at t10 (computed 208.55 at the minimum distance).
7. The first 30 time units are robust to integrator tolerance (agreement 1e-8 at t = 32); beyond about 50 the trajectory depends on tolerance (differences of 1e-2 at t = 62).

## 8. Reconciliation against project code

Searched `src/cyclerfinder` for Burrau, Pythagorean and three-body regularisation on 2026-10-05: no hit. The project has no general (non-restricted) three-body propagator in `src/`; `core/cr3bp_regularized.py` is a Sundman option for the restricted problem. Nothing in the project is sourced from this paper.

## 9. Techniques applicable to the project's problems

All INFERRED; nothing built.

**#931 (a published general three-body control).** This is the standard published general three-body control with a long, well-documented sequence and a published efficiency table: build the problem (masses 3, 4, 5 at rest at the Pythagorean triangle) as the first test of any general three-body integrator: (i) energy and angular momentum and centre-of-mass invariants to the paper's 1e-10; (ii) the 17 close-approach times and pairs of Table I to 2e-3 and the one-figure distances, with the t10 = 15.83 approach at about 4e-4 and the relative-speed and distance scales of the 8-digit value; (iii) the escape of body 1 with the (2,3) binary at period about 0.9 and a hyperbolic relative orbit by t = 60 to 70 (qualitative only: sensitive after t = 50); (iv) the reversal at t = 32 to 1e-10 in positions. The companion periodic orbit with a binary collision (1967b digest) is the sharper test because its answer (period, collision at T/2) is known to 1e-10. Control: my scratch run reproduces Table I; the run description is in section 5.
**#928 (regularisation) and #929 (integrator comparison).** The paper is the origin of the Peters 1968 and Aarseth-Zare test problem, so the project's comparison of plain, Sundman and Levi-Civita or KS integrators on the Pythagorean problem has published expected behaviour: a plain integrator needs thousands more steps (R 17631 against L 3450 at similar or better energy) and degrades with each close approach; a time-only transformation (V) helps in step count but not in energy conservation (8e-7 against 2e-11). My scratch run adds a data point: DOP853 with the product-of-distances time transformation reaches 1.2e-10 in about 4900 steps to t = 30, and the Levi-Civita pair integrator reaches 7.8e-13 in 935 steps to t = 16 (INFERRED comparison across different tolerances). For #929, the reversal test (t = 32 and 62) is the published reversibility control; and it demonstrates that a late-time comparison is chaotic, so compare early. The restricted-problem projects (#928 close lunar passes) are not general three-body problems; the Pythagorean orbit tests the regularisation machinery (pair selection, switching) that the restricted problem does not need.

## 10. Recommended follow-ups (not registered)

1. Add the Pythagorean problem as the first general three-body control of #931, with the Table I times (tolerance 2e-3) and the invariants to 1e-10; include the 1967b periodic orbit.
2. If a Levi-Civita or KS float propagator is built for #928, run this problem with it and compare steps and energy against the DOP853 baseline of section 5 and against Table II.
3. The printed relative speed 171 should be recomputed from a published source if it matters; I did not reproduce it.
4. Szebehely 1967 (PNAS 57, the historical and restricted-problem background, "Paper 1") and Burrau 1913 (A.N. 195:113) are not held.

## 11. Dated note, 2026-10-05 (after digesting the sibling paper)

`docs/notes/2026-10-05-digest-szebehely-peters-1967b-periodic-pythagorean.md` (AJ 72:1187, DOI 10.1086/110398) modifies this paper's initial positions by a differential correction (all velocities zero) to obtain a periodic orbit with a binary collision of bodies 2 and 3 at T/2 = 15.9115, T = 31.8229622453, E = -12.7616527695. Its observation that "a near collision and an almost zero velocity occurring simultaneously" at t = 15.830 suggests a periodic orbit is the observation of section 2 above (r23 about 4e-4 and v1 = 0.04). My integration reproduces that periodic orbit once one printed digit of its Table I is corrected (x_2: -0.0129612126 to -0.0129612186). The same scratch integrator (Levi-Civita pair regularisation) is the one used for the check of section 5 here.
