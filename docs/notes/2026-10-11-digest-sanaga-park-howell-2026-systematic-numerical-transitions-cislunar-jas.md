# Digest: Sanaga, Park and Howell 2026, "Systematic Numerical Transitions of Cislunar Trajectories Across Dynamical Models" (#960 wanted-list row 21)

R. R. Sanaga, B. Park, K. C. Howell (Purdue), The Journal of the Astronautical Sciences 73:10 (2026), 39 pp., DOI 10.1007/s40295-026-00571-5.
Received 26 February 2025, accepted 6 January 2026, published online 10 February 2026. Springer, exclusive licence (not open access).
The brief and the wanted-list row call it "73(1), article 10"; the first page prints "(2026) 73:10", so the issue number is not printed there.
Filed as `cyclers_pdf/papers/sanaga-park-howell-2026-systematic-numerical-transitions-cislunar-trajectories-across-dynamical-models-jas-73-10-doi-10.1007-s40295-026-00571-5.pdf`
(md5 `f604075271d2bf0ddf8507f1cd97333e`, 39 pp., text layer present, Adobe PDF Library; original upload not modified).
Identification checked on page 1: title, three authors, journal, DOI and the three dates agree with the brief. No copy existed before: CORPUS_INDEX.md and
`docs/notes/*.md` were grepped for "Sanaga" and the DOI; the only Sanaga hit is the 2025 Astrodynamics paper below.

Status: READ in full from the text layer (all 39 pages, every equation, both numeric tables, all captions, the appendix and the reference list). Page images were
viewed for pp.10 (Table 2, Fig. 1), 21 (Table 3, Fig. 8), 23 (Fig. 10), 27 (Figs. 13, 14) and 28 (Fig. 15). Figures not viewed are described from their captions and
the running text only, and are marked CAPTION. Numbers are transcribed to `data/sources/sanaga-park-howell-2026-tables.yaml` (three witnesses: text layer read twice,
page render once; all cells agree). COMPUTED = our arithmetic on 2026-10-11; INFERRED = our reading, not printed.
Related digests, cross-referenced and not repeated:
`2026-10-04-digest-sanaga-howell-2025-hill-restricted-four-body-ephemeris-transition.md` (the held companion; Hill model and broken bifurcations in the L2 halo family),
`2026-10-04-digest-park-howell-2024-transition-challenging-regions-er3bp-l2-halo.md`, `2026-10-04-digest-singh-park-howell-2026-aas-26-654-l2-families-intermediary-models.md`,
`2026-10-03-digest-park-howell-2026-averaged-unaveraged-symmetric-periodic-orbits.md`, `2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md`,
`2026-10-04-digest-scheeres-1998-restricted-hill-four-body-problem.md`.

## 0. Verdict

A method paper on carrying Earth-Moon CR3BP trajectories into an ephemeris model. It builds one numerical scheme, the Unified Transition Scheme (UTS), in which the CR3BP,
the ER3BP, the HR4BP, the in-plane Quasi-Hill HR4BP (I-QHR4BP) and a Sun-Earth-Moon ephemeris model (HFEM, DE440) are all written as one equation of motion with 13
time-dependent coefficients in a pulsating-rotating frame. A homotopy parameter epsilon moves between models. Three examples: L2 halo orbits stacked for 15 years (pp.18-21),
an L1 Lyapunov orbit where the number of homotopy steps changes which family the result lies in (pp.21-24), and a lunar free return trajectory (LFRT) swept over epoch in four
models (pp.24-29). No periodic orbit is treated as a cycler or as resonant. There is no table of orbit initial conditions. The only tables are the lunar-frequency
table (Table 2) and a period table (Table 3). See "Relevance to this project" at the end.

## 1. Introduction (pp.1-4)

- Context (p.2): the Moon is the near-term exploration focus (Global Exploration Roadmap, ISECG); the CR3BP is the preliminary model; designs must be validated in an HFEM.
- Two transition routes (p.2): "direct" (discretise the CR3BP orbit into patch points and correct them in the HFEM; used for Gateway baselines, refs 7, 18, 19, 31)
  and "homotopy through intermediary models" (ref 17, Lantoine and Russell 2011; refs 1, 21). A second family uses intermediate Hamiltonian models (ER3BP, (Q)HR4BP;
  refs 3, 6, 28) that keep periodic or quasi-periodic structure, useful when the direct route fails or the result depends on epoch (ref 20).
- Two nuisances (p.3): frames (HFEM inertial J2000 at Earth or Moon; the low-fidelity models rotating, with their own scalings) and independent variables.
- Contributions, three goals (pp.3-4): Goal 1, find the key parameters of a transition (frame, scaling, relation between independent variables). Goal 2, homotopy
  in the pulsating-rotating (pul-rot) frame between the CR3BP and HFEM. Goal 3, use the intermediate models inside the same frame to explain HFEM behaviour as quasi-periodic in epoch.
- Authors state (p.30) that the same UTS was used for invariant tori in ref 26 (Sanaga, Park, Howell, Commun. Nonlinear Sci. Numer. Simul. 152:109323, 2026).

## 2. Pulsating-rotating formulation (pp.4-10)

- Frame (Eq. 1, p.4): origin at the Earth-Moon barycentre; C = [x_hat y_hat z_hat] with x_hat = R_EM / l, z_hat = h / |h|, y_hat = z_hat x x_hat, h = R_EM x V_EM
  the instantaneous angular momentum; l the instantaneous Earth-Moon distance used to scale length so the primaries sit at -mu and 1 - mu on x_hat.
  mu = mu_M / (mu_M + mu_E).
- Independent variable (Table 1, p.5): uniform, dt/dT = t'* = sqrt(mu*/l*^3) with l* constant, or non-uniform, dt/dT = sqrt(mu*/l^3) with the instantaneous l. Table 1 carries
  no numbers. The uniform rate is used with t'* about (4.348 d)^-1 later.
- Models (pp.5-6): CR3BP (circular); HFEM (Earth, Moon, Sun point masses from DE440, coherent in the sense that all bodies obey N-body motion, p.11); ER3BP
  (elliptic Earth-Moon conic, sidereal/anomalistic frequency); HR4BP (coherent Sun, Hill approximation, circular Sun about the barycentre, Earth-Moon as the Lunar Variational Orbit
  of the Hill problem, synodic frequency); I-QHR4BP (Earth-Moon as a 2-D quasi-periodic orbit about the LVO, two frequencies, anomalistic and synodic). All intermediate
  models are coplanar for the Earth-Moon-Sun configuration (p.6).
- Eq. (2) (p.6): rho'' = [b1 b2 b3]^T + [[b4 b5 0], [-b5 b4 b6], [0 -b6 b4]] rho' + [[b7 b9 b8], [-b9 b10 b11], [b8 -b11 b12]] rho + b13 grad(Omega). The common form of all five models.
- Eq. (3): Omega = (1-mu)/rho_Ec + mu/rho_Mc + mu_S/rho_Sc, mu_S = mu~_S / (mu~_M + mu~_E). The Sun term is dropped where the model has no Sun.
- Eq. (4) (p.7): formulas for b1 ... b13 from l, h, the Moon's acceleration A_EM and jerk J_EM relative to the Earth, the barycentre acceleration B'', and the time-scale
  derivatives t', t''. Dimensional inputs come from DE440 for the HFEM. Printed in the paper; no number.
- Eq. (5): CR3BP coefficients b_C = (0, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 1). The b5 = 2 term is the usual Coriolis factor; b7 = b10 = 1 the centrifugal terms; b13 = 1.
- Epoch E is carried as an eighth state (sec. 2.3, App. A) so that every model, autonomous or not, is handled by one multiple-shooting code.
- ER3BP (sec. 2.4, Eq. 6, p.8): e = 0.055, l* is the SEMI-MINOR axis of the lunar ellipse (the mean Earth-Moon distance), a = l*/sqrt(1-e^2); f' = sqrt(mu~_E+mu~_M)(1+e cos f)^2
  (a(1-e^2))^(-1.5); h = l^2 f'. Eq. (6) lists b_E: b1 = b2 = b3 = b6 = b8 = b9 = b11 = 0 and the remaining ones as in Eq. (4) with the planar ER3BP simplifications.
  Epoch variable f with E_dot = f' / t'. Footnote 2 says choosing the semi-minor axis is arbitrary and the semi-major axis is the obvious alternative.
- (Q)HR4BP (sec. 2.5, Eqs. 7, 8, 9, pp.8-10): the Hill problem EOM in the Hill rotating frame (xi, eta, zeta) are re-expressed in t (Eq. 7) via tau' = (1/m) sqrt(mu~_S / l_SB*^3)
  (read as printed from a garbled glyph line; INFERRED reading) with m = 0.080848933808312 and l_SB* = 149598058.09228115 km; Eq. (8) maps velocities xi_dot = xi' tau_dot etc.
  The Earth-Moon state in the HR4BP is read either by integrating Eq. (7) or from a Fourier series in the single angle theta2; in the I-QHR4BP from a 2-D Fourier series in (theta1, theta2)
  (Eq. 9: xi(theta1, theta2) approx sum c_jk exp(i j theta1) exp(i k theta2); the paper writes the sums over 1 <= j <= N_j, a slip for the usual symmetric index range, INFERRED).
  Fig. 1 (viewed, p.10): the 2-D QPO in the Hill rotating frame (xi, eta both within +-0.2 ndim). In panel (a) the fixed-theta1 curves are nearly circles of radius about 0.17, theta1
  sets the osculating eccentricity; in panel (b) the fixed-theta2 curves are small closed loops placed around the ring, theta2 is the Sun-Earth-Moon phase.

## 3. Angular representation for epoch (pp.10-11, Table 2, Eq. 10)

Table 2 (data file): the first three lunar frequencies from DE440 (nu1 = 0.991548233570925, nu2 = 0.925198775804741, nu3 = 1.004021856216999, ndim, reference JD0 = 2.2875165e6, 28 Nov 1550)
and initial phases (theta10 = 3.207895704735595, theta20 = 1.877120527070590, theta30 = 1.619525974729288 rad), with the meanings theta = 0: Moon at perigee, Earth-Moon-Sun collinear,
local maximum of zeta. Eq. (10): theta_i = nu_i k (JD - JD0) + theta_i0. COMPUTED: the three frequencies equal sidereal/anomalistic, sidereal/synodic and sidereal/draconic month ratios
(0.99154811, 0.92519868, 1.00402176) to 1e-7; so the ndim unit here is the sidereal rate. Not printed in the paper; our identification.
The forward map JD to angles is trivial; the inverse (find an epoch with chosen angles) can have no solution within the ephemeris span once two or more angles are targeted (p.11).
Model-to-angle map (p.11): ER3BP to nu1, theta1; HR4BP to nu2, theta2; I-QHR4BP to both; HFEM has all (the paper says five dominant frequencies in Earth-Moon motion, p.27).

## 4. HFEM formulation in different frames (pp.11-14, Eqs. 11-15, Figs. 2, 3)

- N-body EOM about a central body P_q: R_qc'' = -(mu~_c + mu~_q) R_qc / R_qc^3 + sum_j mu~_j (R_cj/R_cj^3 - R_qj/R_qj^3). Propagated in J2000 with SPICE and DE440; Earth, Moon, Sun only.
- Experiment: Eq. (11) an MCI state (R = [20000, 64177, 15938] km, V = 0.02 x + 0.01 y km/s) at 2025 JAN 19, 10 days, propagated in the MCI frame and in the ECI frame (via Eq. 12, R_Ec = R_EM + R_Mc).
  Figs. 2a, b: visually identical. Fig. 2c: the isochronous position difference is not zero. Fig. 2d: norm of the difference between the Earth-Moon relative position from
  Eq. (13) (coherent Earth-Moon motion inside the Sun-Earth-Moon problem) and DE440. CAPTION only; the plotted magnitudes are not read.
- Eq. (13): R_EM'' = -(mu~_M+mu~_E) R_EM / R_EM^3 + mu~_S (R_MS/R_MS^3 - R_ES/R_ES^3). The ECI and MCI EOM agree exactly if and only if the Earth-Moon motion obeys Eq. (13); DE440 does not,
  because it contains more bodies, body shapes and radiation pressure. Hence the HFEM trajectory depends on the central body chosen.
- Pul-rot version (p.13): the same state converts to rho = (0.937090267094130, -0.158289097365981, -0.0427424365372180), rho_dot = (-0.161445581784592, 0.0422875877217298, -0.00349116306688774) ndim.
  Eqs. (14) and (15): rho = l C (R_Ec) + R_BE and rho = l C (R_Mc) + R_BM relate inertial and pul-rot positions (the printed text layer loses a hat or bracket, so the exact grouping is
  INFERRED as rho = (1/l) C^T [R + R_B...]; do not code from this digest, read page 14). Fig. 3 plots the pul-rot difference and the Earth-Moon distance difference.
- Consequence stated (p.14): trajectories built in the pul-rot formulation are not exactly continuous in the inertial formulations used by Copernicus and GMAT; a corrections step fixes it.
  "No single formulation is universally superior."

## 5. Homotopy in the pul-rot frame (pp.15-24)

### 5.1 Multiple shooting (Eqs. 16-22)
- Free variables: n patch-point states X_i (6 each), n epochs E_i, n-1 segment times dt_i, total (8n - 1). Constraints F (Eq. 16): state continuity X_if - X_(i+1) = 0 and epoch continuity
  E_if - E_(i+1) = 0, total 7(n-1). Underdetermined, so the update is minimum-norm: X_(j+1) = X_j - DF^T (DF DF^T)^-1 F (p.16). Tolerance 1e-10 on the infinity norm of F.
- Eq. (17) block Jacobian DF; Eqs. (18)-(22) its entries from the 7 x 7 STM (state plus epoch): dF_X/dX = [Phi(1:6,1:6) -I], dF_X/dE = Phi(1:6,7), dF_X/d dt = X_dot_if, dF_E/dE = [Phi(7,7) -1],
  dF_E/d dt = E_dot_if. A fixed epoch removes its column.
- Minimum-norm choice means "geometrically similar to the initial guess" results, not a unique answer (p.15).

### 5.1b Homotopy (Eqs. 23-26)
- Eq. (23) b_i = (1 - eps) b_C,i + eps b_M,i, Eq. (24) b_i' = eps b'_M,i (the derivatives enter the STM). eps = 0 is the CR3BP, eps = 1 the target model.
- Eq. (25): the solar potential becomes eps mu_S / rho_Sc. With the non-uniform time rate b13 = 1 in every model, so the homotopy is straightforward. With the uniform rate eps also scales b13, so
  Eq. (26) b13 grad Omega = ((1-eps) b_C,13 + eps b_M,13) grad Omega_C + eps b_M,13 grad(mu_S/rho_Sc) keeps the acceleration linear in eps.

### 5.2 Single-step MYTs (pp.18-19, Figs. 5, 6)
Multi-Year Trajectory (MYT): stacked CR3BP revolutions corrected in the target model. L2 halo of period 14.1974 d (about 14.19), apolune start, 3 patch points per rev, 386 revs for 15 y,
epoch 1 Jan 2025 00:00 UTC; epoch scaled so that E1 = 1. Non-uniform rate keeps the xy-plane crossings centred on the CR3BP fixed point but shortens the mean perilune interval to
14.1678 d (uniform: 14.1983 d, CR3BP: 14.1974 d). With uniform rate t'* = (4.3484 d)^-1 the transition is single step (eps 0 to 1). Fig. 6 (CAPTION): crossings in ER3BP, (Q)HR4BP, and the
I-QHR4BP at theta1 about 278 deg, theta2 about 17.67 deg; the I-QHR4BP crossings are closest to the HFEM ones because both use the two dominant frequencies. Both periodic models
give "QPO-like" stacks because the invariant curve is not targeted.

### 5.3 Tuning parameters (pp.20-21, Figs. 7, 8, Table 3)
- The dimensional scale t of the epoch variable (E = T / t) changes the corrected result because it rescales the STM. Fig. 7 (L2 halo, period 14.765 d, 371 revs, epoch 1 Jan 2026): smaller
  E1 stays closer to the CR3BP fixed point; effect larger for the uniform rate. All solutions satisfy the constraints.
- Table 3 (data file; three witnesses agree): mean perilune interval of the 15-y HFEM analog versus t'*. Test 1 (CR3BP 14.1974 d): 14.1405, 14.2006, 14.2571 d at t'* = (4.330, 4.348, 4.366 d)^-1.
  Test 2 (CR3BP 14.50 d): 14.4421, 14.5028, 14.5588 d. COMPUTED: the interval increases by about 0.057 d across each 0.018 d step in the time unit, i.e. nearly proportionally
  (0.018/4.348 = 0.41 percent of 14.2 d = 0.059 d); so the HFEM analog period scales with the chosen time unit, which is what the paper says ("t'* guides the targeting process to locate structures with varying underlying periods").
  The nominal 4.348 d unit gives the closest match, so 4.348 d is the CR3BP characteristic time of the case.

### 5.4 Homotopy benefit and its limit (pp.21-24, Figs. 9, 10)
- L1 Lyapunov orbit, period 17.0855 d, perilune start, 5 patch points per rev, 15 y from 1 Jan 2025, uniform rate: direct (1 step) analog goes strongly non-planar; a 3-step homotopy stays near the planar orbit
  (Fig. 9, CAPTION and text).
- Monodromy eigenvalues (printed) [408.69, 0.002, 1, 1, 1.14, 0.87] (real pair near 1.14 and 0.87) and for a nearby orbit at 17.2594 d [408-ish 393.55, 0.002, 1, 1, 0.99 +- 0.12 i]: the centre pair goes from real to
  complex, a stability change that spawns the axial family (ref 10, Doedel et al. 2007, held). COMPUTED: the 17.0855 d list has product 0.81, not 1, because 0.002 is printed to one digit; consistent with a symplectic spectrum.
- Reversal with the non-uniform rate: direct transition stays planar while a 20-step homotopy ends in the out-of-plane analog (Fig. 10a, viewed: 20-step crossings are scattered in z over about +-0.03 ndim,
  1-step crossings stay within about +-0.003). Fig. 10b (viewed): continuing the direct HFEM result back to eps = 0 with 20 steps does not return to the original Lyapunov orbit; the eps versus initial-y hodograph shows a
  broken bifurcation. The authors read this as a crossing between the out-of-plane centre manifold of the Lyapunov family and the nearby axial family, unconnectable without symmetry (refs 12, 16, 21, 25).
- Lessons (p.23): results depend on number of homotopy steps and time rate; sensitivities tie to broken bifurcations; no single best formulation; computational cost not analysed.

## 6. Sample transfer trajectories: the lunar free return (pp.24-29, Eq. 27, Figs. 11-16)

- LFRT: leaves a 170 km circular Earth orbit by a tangential TLI, passes the Moon, returns to a perigee, no further manoeuvre (ref 15, Jesick and Ocampo 2011). Fig. 11 (CAPTION): the planar CR3BP family at fixed 170 km
  perigee, coloured by perilune altitude; the chosen member has perilune 1005.4 km and tangential TLI 3.1587 km/s. The captions of Figs. 13-16 say perilune 1000.5 km; the two values are not reconciled in the paper (data file).
- Count (p.26): intermediate models: free = (x, y, xdot, ydot post-TLI, ToF, E1) = 6; constraints = 2 perigee conditions + 2 perigee altitudes + 1 perilune altitude = 5; so a one-parameter family in E1. HFEM adds
  z, zdot and z = 0 at start and end (zdot nonzero), keeping the one-parameter family in E1.
- Transition: 5 patch points (X1 and X5 perigees, X3 perilune), 5 homotopy steps; the number of steps and time rate do not change the eps = 1 answer because the solution is isolated at each epoch (p.26).
  Eq. (27): the dimensional altitude constraints rho_Ec^2(E_i) - (rho_Ec,d / l(E_i))^2 = 0 and rho_Mc^2(E_i) - (rho_Mc,d / l(E_i))^2 = 0, since l changes with epoch.
- Results (CAPTION and text; Figs. 13, 14 and 15 viewed): ER3BP, TLI and ToF are approximately sinusoidal in theta1 over 2 pi (Fig. 13); HR4BP repeats every pi in theta2, with smaller swing (Fig. 14).
  The ER3BP swing in the Delta-v versus ToF plane is larger, so the pulsation (theta1) matters more than the Sun (theta2). HFEM (Fig. 15): a solution every 12 h for over 7900 days from 1 Jan 2025;
  the Delta-v versus ToF region is a filled band (Delta-v about 3.150 to 3.165 km/s, axis labels) looking like a sum of the ER3BP and HR4BP curves; the theta1-theta2 colour map is nearly periodic in theta2 with the larger
  variation along theta1. The theta pair with the ER3BP and HR4BP minima is near the HFEM local minimum (magenta star). I-QHR4BP on a 120 x 120 grid reproduces the HFEM map (Fig. 16c, CAPTION and text).
  Conclusion (p.29): a trajectory-correction study can use the ER3BP to see the epoch-dependent perturbation from pulsation.

## 7. Concluding remarks (pp.29-30) and Appendix A (pp.30-37)

Conclusion restates UTS, MYT construction, key parameters, intermediate-model decomposition and refers to ref 26 for invariant tori and broken bifurcations.
Appendix A: Eq. (A1) the 7 x 7 augmented STM (epoch as state); (A2) dPhi/dt = A Phi; (A3) matrix A with Omega second derivatives and the b coefficients; (A4)-(A6) dx''/dT etc.; (A7) dE_dot/dT; (A8) t' and t'' for
the uniform and non-uniform rates. A.1 HFEM: the b_i' (derivatives of b_i) in terms of l, h and their derivatives, with l', l'', l''', h'', B''' and S_EM from DE440; the Sun position in the pul-rot frame.
A.2 ER3BP: the same with ER3BP simplifications and no Sun, homotopy forms (A12), (A13). A.3 (Q)HR4BP: the Hill-approximation solar term and its derivatives, the i_H unit vector, homotopy forms (A14), (A15).
Equations in A.1 to A.3 were read but not transcribed or re-derived; anyone coding the STM should work from the page images (the text layer is partly garbled in the fraction stacks).
Data availability: "No datasets were generated or analysed". Funding: Purdue and NASA JSC (80NSSC18M0122).

## 8. Equation roles (numbered equations only)

| Eq. | Role |
|---|---|
| 1 | direction cosine matrix of the pul-rot frame |
| 2 | common EOM, 13 coefficients, all models |
| 3 | potential including Sun term |
| 4 | b1 to b13 formulas (HFEM and general) |
| 5 | CR3BP coefficients |
| 6 | ER3BP coefficients |
| 7 | Hill problem EOM in t |
| 8 | velocity relation between Hill time and pul-rot time |
| 9 | 2-D Fourier representation of the QPO (I-QHR4BP) |
| 10 | angle versus epoch, uses Table 2 |
| 11 | MCI test state of the frame experiment |
| 12 | MCI to ECI translation |
| 13 | coherent Earth-Moon motion inside Sun-Earth-Moon |
| 14, 15 | inertial to pul-rot position map |
| 16 | multiple-shooting constraints |
| 17 to 22 | Jacobian of the corrector |
| 23, 24 | homotopy of coefficients and derivatives |
| 25 | homotopy of solar potential |
| 26 | b13 grad Omega under uniform time rate |
| 27 | altitude constraints of the LFRT |
| A1 to A15 | STM, A matrix, coefficient time derivatives, homotopy forms in the three model groups |

## 9. Figures (one line each)

1 QPO ring in Hill frame (viewed). 2 and 3 frame-formulation differences for the eq. 11 test (CAPTION). 4 multi-shooting schematic. 5 uniform versus non-uniform rate HFEM crossings (CAPTION).
6 crossings in ER3BP, (Q)HR4BP, I-QHR4BP (CAPTION). 7 epoch scale t (CAPTION). 8 t'* (viewed): crossings are small linear streaks, nearest the yellow CR3BP point at 4.348 d. 9 and 10 L1 Lyapunov, direct versus homotopy (10 viewed).
11 LFRT family in the barycentre rotating frame. 12 LFRT patch-point layout. 13, 14, 15, 16 LFRT characteristics (13, 14, 15 viewed).

## 10. Relevance to this project

**Does the paper treat any periodic orbit as a cycler, repeated-moon or resonant object? No.** Evidence:
- The words "cycler" and "resonance/resonant" occur nowhere in the body (pp.1-37). "Resonant" appears only inside the title of reference 6 (Chujo and Takao, "Synodic resonant halo orbits of solar sails") and
  "9:2 NRHO" inside the title of reference 5. The 3:1 interface region and the 9:2 NRHO of the 2025 companion paper are not mentioned. Grep of the text layer, 2026-10-11.
- The periodic orbits used are an L2 halo (periods 14.1974, 14.19, 14.50, 14.765 d), an L1 Lyapunov (17.0855 d, nearby 17.2594 d) and the 5.4 monodromy spectra. They are stacked 371 to 386 times
  for 15 years as test subjects for corrections convergence (secs. 5.2 to 5.4). They are not families searched for, not counted by revolutions, not classified by synodic or sidereal ratio, and no period ratio is quoted.
- The LFRT (sec. 6) is a one-shot transfer with a fixed 170 km perigee and a fixed perilune altitude; it is a free-return transfer, not repeated.

**Labels the content supports** (topology label vocabulary of the literature gate; I could not confirm the closed set, so these are descriptive):
- `{model-transition}`: supported (title, abstract, secs. 2 to 6). This is the paper's subject.
- Also supportable as secondary content: periodic-orbit (halo and Lyapunov) analogs, ER3BP and HR4BP and I-QHR4BP intermediate models, lunar free-return transfer, broken bifurcation.
- NOT supported: any cycler, quasi-cycler, repeated-moon, resonant (p:q) or low-energy-transfer label; no V-infinity, no resonance order, no encounter sequence.

**Working-body scope.** Sec. 2.2: all models treat bodies as point masses with a restricted spacecraft. CR3BP and ER3BP have two working bodies (Earth, Moon). HR4BP, I-QHR4BP and HFEM add the Sun as a
perturbing third body (p.5, "Earth-Moon-Sun"); Eq. (3) shows the Sun only as an added term in Omega, and in the Hill models the Sun acts through the tide and through the Earth-Moon motion (sec. 2.5). So the
working-body set is {Earth, Moon}, with the Sun a perturbation present only in some models. That fits the anchor's current `primary="Earth"`, `body_set={"Moon"}`. If the gate has a scope field that
counts bodies (the n_bodies_scope of #1052), the paper's content is scope 2 for the CR3BP and ER3BP analogs and scope 3 for the HR4BP, I-QHR4BP and HFEM analogs; I would record 2 as the working-body scope with the Sun noted as a perturber.
Not a repeated-moon paper and not Earth-Moon-only in the n-body count; do not use it to exclude a single-body or two-body-scope candidate under rule U.

**Recommendation on the gate anchor** (src/cyclerfinder/search/literature_check.py lines 1583-1598; not edited here):
1. Keep `topology_label = {"model-transition"}`. The title-only label is confirmed by the content.
2. Replace the comment "title-level only ... no digest in the corpus yet" with a pointer to this digest and the Filed-as copy.
3. The anchor must not act as prior art for any cycler, resonant or repeated-moon signature. Its keywords are about the fidelity transition framework only, which is right.
4. If the label vocabulary allows more than one label, adding a periodic-orbit-analog label (L2 halo and L1 Lyapunov) and a free-return label is accurate, but neither affects novelty of a cycler candidate.
5. The companion Sanaga and Howell 2025 (held) is the one that handles resonant L2 halos (the 3:1 interface and broken bifurcations); if any anchor in the gate needs a resonant label from this author group it is that paper, and even there the digest says "Not a cycler paper".

Other uses: (a) the Eq. (13) incoherence result (DE440 differs from Eq. (13) Earth-Moon motion, so ECI and MCI propagations differ) is a citable caution for any project check that compares a coherent-model orbit with
an ephemeris propagation. (b) Table 3 shows an HFEM analog's period scales with the chosen nondimensional time unit, so any comparison of a catalogued CR3BP period with an ephemeris analog's period
must say which time unit was held. (c) The L1 Lyapunov example is a worked case where the number of homotopy steps selects the family, relevant to #1047 families judged only by a single-step continuation.

## 11. Internal inconsistencies and cautions

- Perilune altitude 1005.4 km (text, Fig. 11 caption) versus 1000.5 km (Figs. 13-16 captions).
- Mean perilune interval for the L2 halo of period 14.1974 d at uniform rate: 14.1983 d (p.18, epoch 2025) versus 14.2006 d (Table 3, epoch 2026); different runs, the paper does not say so.
- p.8 prints the Hill ndim-time rate with a 1/m factor in a garbled line; see the Sanaga and Howell 2025 digest sec. 2.2 for the same "1/m" versus "(1+m)" ambiguity. Do not code from the text layer.
- "Earth-Moon ... planar" for the intermediate models while the LFRT HFEM is non-planar; the comparison in Fig. 15 mixes planar and non-planar families (the paper notes the extra z, zdot variables).
- The paper's "five dominant frequencies" (p.27) versus the three tabulated in Table 2; only two are used.

## 12. Citation mining

References examined (all 31). Held means a file in `cyclers_pdf/papers/` and a CORPUS_INDEX row, checked 2026-10-11 by grep and filename.

Bearing on cyclers, resonant Earth-Moon orbits, model transition, or Earth-Moon periodic-orbit families:

| Ref | Item | Status |
|---|---|---|
| 3 | Boudad, Howell, Davis, AIAA SciTech 2022-1276, "Analogs for Earth-Moon halo orbits and their evolving characteristics in higher-fidelity force models" | NOT held (the held Boudad 2020 ASR 66:2194 is a different paper); not on the wanted list |
| 5 | Brown, McCarthy, Peterson et al., IAC 2024, "Dynamics around the 9:2 NRHO in the Sun-Earth-Moon system: staging, phasing, transitioning" | NOT held; not on the list |
| 6 | Chujo, Takao, JSR 59(6):2129-2147 (2022), DOI 10.2514/1.A35324, "Synodic resonant halo orbits of solar sails in restricted four-body problem" | NOT held; not on the list |
| 9 | Dei Tos, Topputo, Acta Astronaut. 136:236-247 (2017), DOI 10.1016/j.actaastro.2017.02.025, hierarchical structure of astrodynamical models | NOT held |
| 10 | Doedel et al., IJBC 17(8):2625 (2007) | held |
| 11 | Gomez, Masdemont, Mondelo, A&A 390(2):733-749 (2002), DOI 10.1051/0004-6361:20020625, "Solar system models with a selected set of frequencies" | NOT held |
| 12 | Henry, Rosales, Brown et al., AAS 23 (Big Sky, Aug 2023), "Quasi-periodic orbits around Earth-Moon L1 and L2 in the Hill restricted four-body problem" | NOT held |
| 13 | Hill 1878, Am. J. Math. 1(2):129 | classical, not checked here |
| 15 | Jesick, Ocampo, JGCD 34(1):98-106 (2011), DOI 10.2514/1.50550, symmetric lunar free-return generation | NOT held (the held Jesick 2019 is Mars Trojans) |
| 16 | Jorba-Cusco, Farres, Jorba, Front. Appl. Math. Stat. 4:32 (2018) | held |
| 17 | Lantoine, Russell, JAS 58(3):335-363 (2011) | held (two files) |
| 18 | McCarthy, Howell, JAS 70(4):24 (2023), DOI 10.1007/s40295-023-00389-5, heteroclinic connections between QPOs | NOT held |
| 20 | Park, Howell, CMDA 136:6 (2024), DOI 10.1007/s10569-023-10178-9, "Assessment of dynamical models for transitioning from the CR3BP to an ephemeris model" | NOT held (the held Park-Howell 2024 file is the AIAA SciTech paper, a different one) |
| 21 | Park, Howell, Adv. Space Res. (2025), DOI 10.1016/j.asr.2025.01.009, L2 halo analogs with the ER3BP | NOT held |
| 22 | Park, Sanaga, Howell, CMDA 137:5 (2025), DOI 10.1007/s10569-025-10236-4, frequency-based hierarchy of models | NOT held; the central model-hierarchy source |
| 23 | Park, Sanaga, Howell, Acta Astronaut. (2025), DOI 10.1016/j.actaastro.2025.03.033, QHR4BP | NOT held |
| 25 | Sanaga, Howell, AIAA SciTech 2024-1451, DOI 10.2514/6.2024-1451, challenging region in the L2 halo family via HR4BP | NOT held (the held 2025 Astrodynamics paper is its journal extension) |
| 26 | Sanaga, Park, Howell, Commun. Nonlinear Sci. Numer. Simul. 152:109323 (2026), DOI 10.1016/j.cnsns.2025.109323, UTS for invariant tori | NOT held |
| 27 | Scheeres 1998 CMDA 70:75 | held (digest exists) |
| 28 | Scheuerle, Howell, Davis, Adv. Space Res. (2024), DOI 10.1016/j.asr.2024.07.035, energy-informed ballistic lunar transfers | NOT held |
| 31 | Zimovan-Spreen, Scheuerle, McCarthy et al., AAS 23 (Big Sky), baseline orbit generation for NRHOs | NOT held |

Not relevant to the project (not proposed): refs 1 (Almanza-Soto thesis), 2 (Berry-Coppola indirect term), 4 (Brown lunar theory), 7 (Davis 2017 stationkeeping), 8 (Dei Tos thesis), 14 (Hou-Liu triangular points),
19 (Oguri Equuleus), 24 (DE440), 29 (Copernicus), 30 (Wintner). Ref 7 (Davis et al., AAS 17-826) is the "transition region" source cited in the Sanaga-Howell 2025 digest; it was left out here because it concerns stationkeeping.

### Proposed wanted-list rows (not added; the lead decides)

Priority order, by bearing on the project. Identifier status is as printed in this paper's reference list; none was resolved against Crossref by me.

1. Park, B., Sanaga, R. R. & Howell, K. C. (2025), "A frequency-based hierarchy of dynamical models in cislunar space: leveraging periodically and quasi-periodically perturbed models", Celest. Mech. Dyn. Astron. 137:5, DOI 10.1007/s10569-025-10236-4 (printed). Why: the source of the ER3BP, HR4BP and I-QHR4BP frequency assignments (nu1 anomalistic, nu2 synodic) used throughout this paper; the underlying model-hierarchy paper.
2. Park, B. & Howell, K. C. (2024), "Assessment of dynamical models for transitioning from the CR3BP to an ephemeris model with applications", CMDA 136:6, DOI 10.1007/s10569-023-10178-9 (printed). Why: the pul-rot formulation and the direct versus intermediate-model comparison this paper builds on.
3. Gomez, G., Masdemont, J. J. & Mondelo, J. M. (2002), "Solar system models with a selected set of frequencies", A&A 390(2):733-749, DOI 10.1051/0004-6361:20020625 (printed). Why: origin of the frequency-based models and the lunar angles of Table 2.
4. Park, B., Sanaga, R. R. & Howell, K. C. (2025), "Bridging ephemeris transition gaps: quasi-periodic extensions for the Hill restricted four-body problem in cislunar space", Acta Astronaut., DOI 10.1016/j.actaastro.2025.03.033 (printed). Why: the QHR4BP and the 2-D QPO Fourier representation of Eq. (9).
5. Park, B. & Howell, K. C. (2025), "Characterization of Earth-Moon L2 halo analogs in an ephemeris model using the elliptic restricted three-body problem", Adv. Space Res., DOI 10.1016/j.asr.2025.01.009 (printed). Why: ER3BP halo analogs; the held SciTech 2024 paper is the conference form of related work, not this one.
6. Sanaga, R. R., Park, B. & Howell, K. C. (2026), "Unified transition scheme for invariant tori across various models in the cislunar domain", Commun. Nonlinear Sci. Numer. Simul. 152:109323, DOI 10.1016/j.cnsns.2025.109323 (printed). Why: extends the UTS to invariant tori; closest to the project's quasi-periodic-orbit interests.
7. Henry, D. B., Rosales, J. J., Brown, G. M. et al. (2023), "Quasi-periodic orbits around Earth-Moon L1 and L2 in the Hill restricted four-body problem", AAS/AIAA Astrodynamics Specialist Conference, Big Sky, August 2023 (no paper number printed). Why: QPOs in the HR4BP; companion to the held Brown digests.
8. Chujo, T. & Takao, Y. (2022), "Synodic resonant halo orbits of solar sails in restricted four-body problem", JSR 59(6):2129-2147, DOI 10.2514/1.A35324 (printed). Why: the only reference with "resonant" in its title; synodic-resonant halo orbits in a four-body model. Solar-sail, so low priority.
9. Boudad, K. K., Howell, K. C. & Davis, D. C. (2022), AIAA SciTech 2022-1276, DOI 10.2514/6.2022-1276 (printed). Why: halo analogs in higher-fidelity models; low priority (the held Boudad 2020 covers the synodic resonant NRHO).
10. Sanaga & Howell (2024), AIAA SciTech 2024-1451, DOI 10.2514/6.2024-1451 (printed): the conference form of the held 2025 Astrodynamics paper; low priority, the held journal version supersedes it.
11. Dei Tos & Topputo (2017), Acta Astronaut. 136:236, DOI 10.1016/j.actaastro.2017.02.025 (printed): the hierarchy-of-models paper. Low priority.
12. Jesick & Ocampo (2011), JGCD 34(1):98-106, DOI 10.2514/1.50550 (printed): lunar free-return generation. Low priority (not a cycler paper; useful only for LFRT context).
