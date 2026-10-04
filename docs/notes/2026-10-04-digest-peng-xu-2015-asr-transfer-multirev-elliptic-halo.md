# Digest: Peng & Xu (2015), "Transfer to a multi-revolution elliptic halo orbit in Earth-Moon elliptic restricted three-body problem using stable manifold"

Single-paper digest, written 2026-10-04 (Sydney time) from a read of all 13 pages. Filed in the private paper corpus as
`peng-xu-2015-transfer-multirev-elliptic-halo-earth-moon-ER3BP-stable-manifold-asr-55-1015-doi-10.1016-j.asr.2014.11.013.pdf`.
The file has a text layer; equations and tables were read from the page images (130 dpi) where the text layer was ambiguous.
Statements are marked READ (printed), COMPUTED (arithmetic or a run shown) or INFERRED (reasoning, not printed).
Journal page numbers are used throughout (PDF page n is journal page 1014 + n).

## 0. Headline

- Peng, H. and Xu, S., Advances in Space Research 55(4):1015-1027 (2015), DOI 10.1016/j.asr.2014.11.013. Beihang University,
  Beijing. Received 21 May 2014, revised 8 September 2014, accepted 18 November 2014, online 26 November 2014. (This paper appeared
  before its companion, the Celestial Mechanics and Dynamical Astronomy paper, digested in
  `docs/notes/2026-10-04-digest-peng-xu-2015-stability-multi-revolution-elliptic-halo.md`; it cites that work as "Peng and Xu 2014",
  a 24th ISSFD paper, "Numerical stability study of multi-circle elliptic halo orbit in the elliptic restricted three-body problem",
  not held. Neelakantan & Ramanan 2022 call this ASR paper "Peng and Xu 2015b".)
- Subject: direct transfers from a 185 km Earth parking orbit to one Earth-Moon L1 ME-Halo orbit (M5N2, periapsis group) along
  its stable manifold, designed directly in the elliptic problem. First transfer study to a multi-revolution elliptic halo.
- THE ONE ORBIT STATE PRINTED IN THIS PAPER (Table 2, p1020), which the CMDA companion does not print:
  `[0.8516666, 0, 0.1832855, 0, 0.2582897, 0]` in scaled units, mu = 0.0122, e = 0.0554, T_E = 4 pi. It agrees to every printed digit
  with the 15-digit M5N2 row of Neelakantan & Ramanan 2022 Table 8 (section 6). This is a complete, independently published
  M5N2 state, and the companion paper's controls C4 (section 7 of that digest, "INFERRED candidate state") is now confirmed.
- COMPUTED in this session with `core.er3bp` (section 6): the printed 15-digit state closes from periapsis (1.6e-8) and its
  segmented monodromy has largest eigenvalue 1.54273e6, second real eigenvalue 1.00858 and a unit-circle pair 0.94956 +- 0.31359i,
  against the printed 1.5427e6 and 1.0086 (CMDA p300) and 6.481e-7 and 0.992 (this paper p1020). So the project's pulsating-frame
  equations of motion and variational equations reproduce a published eigenvalue of 1.5e6 at the Earth-Moon eccentricity. This is
  evidence on the model-or-paper question of `#925`, with limits (section 8).
- The transfer numbers (total delta-v 3.791 to 4.676 km/s for perigee manifold injection, 3.388 to 3.934 km/s for locally
  optimised injection, transfer times 54 to 86 days) are read off figures and the text; none is tabulated. They are sourced
  numbers for a future transfer test but only to the printed digits and only with the paper's unreproducible optimiser details
  (section 5).
- Nothing in the paper resolves the M4N2 Lyapunov non-closure in `#925` item 2 (section 8).

## 1. The model as printed (Section 2.1, pp1016-1018)

READ. Third body of negligible mass in the field of two primaries m1 (Earth) and m2 (Moon) on Keplerian ellipses of eccentricity e.
Pulsating synodic frame: origin at the barycentre, x from m1 to m2, z along the primaries' angular momentum, y completing a
right-handed set. The frame is "instantaneously normalized by the primary distance r12(f), the total primary mass (m1+m2) and the
reciprocal of the mean motion n". Primaries fixed at x1 = -mu, x2 = 1 - mu. Independent variable true anomaly f; the epoch with the
primaries at periapsis is f0 = 0.

- r12(f) = a12 (1 - e^2) / (1 + e cos f), a12 the primaries' semimajor axis.
- Eq. 1: d/dt = (df/dt) d/df = sqrt( G (m1+m2) / (a^3 (1-e^2)^3) ) (1 + e cos f)^2 d/df.
- Eq. 2 (dots are d/df): x'' - 2 y' = omega_x, y'' + 2 x' = omega_y, z'' = omega_z.
- Eq. 3: omega(x,y,z,f) = (1 + e cos f)^(-1) Omega_tilde(x,y,z).
- Eq. 4: Omega_tilde = Omega - (1/2) e cos f z^2.
- Eq. 5: Omega = (1/2)(x^2 + y^2) + (1-mu)/r1 + mu/r2 + (1/2) mu (1-mu), r1 = sqrt((x+mu)^2 + y^2 + z^2),
  r2 = sqrt((x-1+mu)^2 + y^2 + z^2).
- Eq. 6 and 7: the velocity-squared integral and its decomposition. x'^2 + y'^2 + z'^2 = 2 integral(omega_x dx + omega_y dy + omega_z dz),
  and since d omega = (omega_x dx + omega_y dy + omega_z dz) + omega_f df, the result is
  x'^2 + y'^2 + z'^2 = 2 omega - 2 integral from f0 to f of omega_f ds - C(f0). The text layer drops the integrand subscript and the
  limits; the form above is from the printed Eq. 6 and the printed differential, not read digit by digit from an image. The
  paper's point (p1017): the system is non-autonomous with period 2 pi in f, "There is no Jacobi constant anymore", and the
  integral term is "caused by the pulsating of the system".
- Libration points: same positions as the circular problem in the scaled pulsating frame, "actually pulsating with the primaries
  along x-axis"; the paper calls them "libration point regions".
- Table 2 parameters (p1020, read from the image): Earth-Moon semimajor axis 384,400 km; Earth-Moon period "27.5 days" (the sidereal
  month is 27.32 days; the paper's value is as printed); mass ratio mu = 0.0122; eccentricity e = 0.0554; target L1 ME-Halo M5N2;
  initial state above; manifold transfer direction External; position perturbation "100 m"; height of parking LEO 185 km; maximum
  transfer duration 6 pi (in f); integration tolerance 1e-10; optimisation solver tolerance 1e-7.

COMPUTED check: `core/er3bp.py:er3bp_eom` matches Eqs. 2 to 4 term by term (section 6 of the companion digest records the same
match). Nothing new here; the closure and monodromy runs in section 6 now exercise it numerically.

## 2. The orbit: construction, groups, relation to the CMDA paper

READ (Section 2.2, pp1018-1019, and Section 4.1 Step 1, pp1022-1023).

- Periodicity criterion (Campagnola 2008, spatial form): a symmetric periodic orbit must cross the x-z plane perpendicularly
  twice. The paper words this as "respectively when primaries are at their periapsis and apoapsis" (p1018). INFERRED: that holds
  for N odd; for N even both crossings are at the same apse (T_E/2 = N pi, so M5N2 has both crossings at periapsis, f = 0 and f = 2 pi).
  The CMDA paper words it correctly ("at an apse"). COMPUTED: the printed M5N2 state has its half-period perpendicular crossing at f = 2 pi
  (section 6), which is periapsis, consistent with N even.
- Eq. 8: T_C = T_E / M = 2 N pi / M, M and N positive integers. T_C the period of the CR3BP halo, T_E the ER3BP period. "Planar and
  vertical Lyapunov orbits with large amplitudes have periods large enough to be an integer multiple of 2 pi (Broucke 1969, Sarris 1989).
  But a halo orbit usually has a period less than 2 pi, so it needs to revolve several revolutions".
- Generation (p1019 and Step 1 of the strategy, p1022): take the halo of period T_C = 4 pi / 5 (the M5N2 case) from the whole
  circular family; raise e gradually; adjust the initial condition by the multi-segment optimisation of the companion paper to close
  it. The method is NOT described here ("Details are out of the scope of this paper"), only named: the differential corrector is
  rewritten "into a multi-segment optimization problem and a more advanced optimization solver, interior point method, is adopted".
  "This new method is much easier to program and relatively less sensitive to the initial guess."
- Table 1 (p1018), four groups by starting state (read from the text layer, matches the CMDA classification):

  | Group | M | f0 at start | Primary at | Starting side in x-y |
  |---|---|---|---|---|
  | Periapsis | odd | 0 | periapsis | left |
  | Apoapsis | odd | pi | apoapsis | left |
  | Left | even | 0 | periapsis | left |
  | Right | even | 0 | periapsis | right |

  Only north orbits are considered; south ones are symmetric. "group" is used because the orbits form a family with respect to e,
  distinguished from the circular halo family.
- Fig. 2 (p1019): four examples in the non-pulsating frame, with axes in km: (a) L1 periapsis M5N2 (x about 3.0e5 to 3.8e5 km,
  y within about +-5e4 km, z up to about 6e4 km), (b) L1 apoapsis M5N2, (c) L2 left M2N1 (x about 3.6e5 to 4.6e5 km), (d) L2 right M2N1.
  No numeric values are printed in text. "The periapsis ME-Halo orbit stretches to larger x after starting point because the
  primaries are pulsating away mutually, but the apoapsis ME-Halo orbit stretches to smaller x because the primaries are getting
  closer."
- Which orbits have two real eigenvalue pairs (p1020): "Earth-Moon L1 perigee ME-Halo with M5N2 and L2 Left ME-Halo with M2N1 have
  two real pairs of eigenvalues". The first is the target. (The apoapsis M5N2 orbit has one real pair, per the CMDA paper p300.)
- Relation to the CMDA paper: this paper uses the CMDA paper's method and its Earth-Moon result (the L1 periapsis M5N2 orbit with
  two real pairs, lambda1 about 1.5427e6) and builds the transfer on it. The target orbit's eigenvalues here (lambda1 = 6.481e-7,
  lambda2 = 0.992 for the two eigenvalues below one) are the reciprocals of the CMDA's 1.5427e6 and 1.0086 (checked in section 4).
  No continuation, no stability sweep and no other (M, N) appear in this paper.

## 3. Stable manifold of a periodic orbit of a non-autonomous system (Section 3, pp1019-1022)

READ.

- Monodromy (Eq. 9, 10): Phi'(f, f0) = A_E(X, f) Phi(f, f0), Phi(f0, f0) = I6, A_E = [[0, I3], [H_E, K_E]] evaluated on (X, f), H_E the
  Hessian of omega (printed with entries omega_xx ... omega_zz), K_E = [[0, 2, 0], [-2, 0, 0], [0, 0, 0]]. Psi = Phi over one period T_E
  with six eigenvalues in reciprocal pairs lambda1, 1/lambda1, lambda2, 1/lambda2, lambda3, 1/lambda3.
- Circular halo (p1020): one real pair (lambda1 < 1), a pair of unit eigenvalues (lambda2 = 1), a complex unit pair. Elliptic: "as many as
  three pairs of real eigenvalues". The paper's claim: a halo's stable manifold is 2-dimensional in the CR3BP, but an ME-Halo with a
  second real eigenvalue below one has a 3-dimensional stable manifold. A "redundant stable direction".
- Orbit parameterisation (p1020): tau = f / T_E in [0, 1), tau = 0 the start, tau = 1 the return. For the M5N2 orbit "the length for
  each revolution is Delta tau = 0.2 roughly" (Fig. 3).
- Stable directions: eigenvectors v^s of Psi(tau) for eigenvalues below one. The manifold branch is generated by perturbing the orbit
  state along the position part r^s of the eigenvector by epsilon and the velocity part by epsilon ||u^s|| / ||r^s||, then integrating
  BACKWARD. The target has eigenvalues below one lambda1 = 6.481e-7 (main stable direction v^s_1) and lambda2 = 0.992 (redundant
  stable direction v^s_2). Eigenvalues are invariant along the orbit, the eigenvectors are not, so an orbit injection point tau is
  fixed.
- Skew frame (Fig. 4): the two stable directions form a skew frame with an angle gamma (varying with tau), spanning a 2-D hyperplane
  tangent to the orbit in the 6-D phase space. Four quadrants of the angle beta: (0, gamma), (gamma, 180 deg), (180 deg, gamma + 180 deg),
  (gamma + 180 deg, 0). Quadrant one gives manifolds leading to the Moon backward ("internal"), quadrant three to the Earth backward
  ("external"). Fig. 5 shows the manifolds for the four quadrants at tau = 0 (position perturbation printed as 100 km in the Fig. 5
  caption and 100 m in Table 2: an internal inconsistency, flagged), integrated backward for Delta f of about 1.82 pi. The manifold along
  v^s_2 stays so close to the orbit that it cannot be seen in Fig. 5 because lambda1 is so much smaller than lambda2.
- Three parameters (tau, beta, sigma): tau the orbit injection (OI) point on the target, beta the stable direction, sigma the epoch
  of a point on that bundle counted backward from the OI point (a true anomaly). "a spacecraft can only be located at the point
  corresponding to the epoch f* + 2 k pi on an ME-Halo orbit or its manifolds" (p1022): the transfer is tied to a specific
  Earth-Moon phase, which the circular problem does not require.

## 4. Transfer design (Section 4, pp1022-1026)

READ.

1. Strategy (p1022-1023, Steps 1 to 5): construct the ME-Halo orbit; pick an OI point tau; choose external or internal; perturb by
   epsilon along the stable direction to obtain the OI state; integrate backward to a manifold injection (MI) point X_MI; define the
   manifold-injection maneuver dV_MI collinear with the manifold velocity V_MI; apply it and integrate backward to the parking LEO of
   height 185 km; adjust dV_MI until the arrival at the LEO is tangential; compute the LEO departure maneuver dV_LEO. The last
   adjustment "is solved by an optimization solver in this paper, whose initial guess is provided by a Hohmann transfer in two-body
   problem" (Alessi et al. 2010 solved it by differential correction). Inclination and position on the LEO are unrestricted.
2. Perigee transfer (4.2): MI point at a perigee of the manifold; the state is converted from the pulsating frame to an Earth-centred
   inertial frame for the perigee test. External transfers go to the Earth directly after dV_MI (Fig. 6, OI points tau = 0.06 to
   0.12); internal transfers meet a lunar flyby first (Fig. 7, tau = 0.14 to 0.2) and were found "very chaotic" combined with the
   multi-revolution property, so only external transfers are surveyed. Cost is the sum of norms dV = ||dV_LEO|| + ||dV_MI||.
3. Survey results, perigee transfer (Figs. 8 and 10): total dV from 3.791 km/s to 4.676 km/s; five peaks and troughs in tau (one per
   revolution of the orbit); transfer time from 56 days to 71 days; "generally speaking the larger the dV is, the faster the transfer
   is". A discontinuity near tau = 0.31 to 0.32 is a constraint artifact, not a property of the orbit: the perigee distance is limited
   to 0.7 normalised length units, and when the manifold perigee is further out "the program automatically jumps to the next
   perigee" (Fig. 9).
4. Local optimal MI point (4.3, Figs. 11 to 14): for the OI point tau = 0.22 the MI parameter sigma was scanned from about -14 to -9.5
   (Fig. 12 axis, read from the image); the minimum total cost is at sigma about -12.5 (the text layer prints "12.5" without the
   sign; the image shows the minimum marker at -12.5), "around the apogee of the manifold", with the perigee (dashed line) near
   sigma = -10.5 where dV_MI and the total are highest and dV_LEO lowest. The MI point is limited between the first epoch with x < 0.8
   and the second epoch with x = 0, a window containing at least one perigee and one apogee. Result over all OI points (Fig. 13):
   total dV from 3.388 km/s to 3.934 km/s, "greatly reduced compared with that of the perigee transfer", consistent with
   Alessi et al. 2010 and Li and Zheng 2010; transfer time (Fig. 14) mostly 54 to 86 days, average 70 days. The authors say other
   local optima exist beyond the window because the manifold circles the Earth many times, and do not seek a global optimum.
5. Redundant stable direction (4.4, Eq. 11): v(tau) = -cos(eta) v^s_1(tau) - sin(eta) v^s_2perp(tau), eta in [-5 deg, 5 deg], where
   v^s_2perp is v^s_2 orthogonalised against v^s_1 and the minus sign selects external directions. Eleven directions, one per degree,
   at each OI point, with MI at the manifold apogee (apogee transfers). Fig. 15: the main-direction (blue) curve repeats the form of
   Fig. 13; around tau = 0.1, 0.3, 0.5, 0.7, 0.9 (the part of each revolution nearest the Moon) the main direction costs the most, and
   around tau = 0.2, 0.4, 0.6, 0.8, 1.0 the cost varies little and the main direction can be the minimum. Fig. 16 (tau = 0.44, the
   most affected region) shows MI points separated over a large region for a small angular spread (the text layer reads "10
   separation"; given eta in [-5, 5] deg in 1 deg steps the intended value is not certain, so no number is taken from it).
   The redundant direction adds a dimension to the optimisation, more computing time and more local minima.
6. Conclusion (p1026): the orbit can serve "as a natural formation flying orbit in the Earth-Moon system, utilizing its multi-
   revolution property", and "as a relatively more stable orbit for a libration point space station"; the redundant stable direction
   "frees this restriction to some extent" (the epoch-phase coupling).

No transfer cost is printed in a table. All delta-v numbers above are range endpoints from the text; the paper does not print the
optimum transfer's departure state, epoch or maneuver vectors. Neelakantan & Ramanan 2022 compare directly with these ranges
(their digest, p8 and the Fig. 10 comparison at tau = 126.648 degrees there: 57.42011 days, 4.6121 km/s on the Peng-Xu side).

## 5. Every printed number usable as a test, with source

Values read from images are marked IMG; from the text layer TXT.

| # | Quantity | Value as printed | Where | Use |
|---|---|---|---|---|
| P1 | mu, Earth-Moon | 0.0122 | Table 2, p1020 (IMG) | model constant for any ME-Halo control |
| P2 | e, Earth-Moon | 0.0554 | Table 2 (IMG) | same |
| P3 | Earth-Moon semimajor axis | 384,400 km | Table 2 (IMG) | scale constant |
| P4 | Earth-Moon period | 27.5 days | Table 2 (IMG) | scale only; true sidereal 27.32 d, so not a time test |
| P5 | M5N2 L1 periapsis state at f0 = 0 | [0.8516666, 0, 0.1832855, 0, 0.2582897, 0] | Table 2 (IMG) | closure and monodromy (section 6) |
| P6 | Orbit period | T_E = 4 pi (M = 5, N = 2); T_C = 4 pi / 5 | Eq. 8 and p1020 (TXT) | exact rational identity |
| P7 | Eigenvalues of Psi below one | lambda1 = 6.481e-7, lambda2 = 0.992 | p1020 (IMG) | monodromy test; lambda2 see below |
| P8 | Revolution length in tau | Delta tau about 0.2 | p1020 | arithmetic: 1 / M |
| P9 | Perigee transfer total dV range | 3.791 to 4.676 km/s | p1023 (TXT) | transfer test only with the paper's optimiser |
| P10 | Perigee transfer duration range | 56 to 71 days | p1024 (TXT) | same |
| P11 | Locally optimal total dV range | 3.388 to 3.934 km/s | p1024 (TXT) | same |
| P12 | Locally optimal duration | mostly 54 to 86 days, mean 70 | p1025 (TXT) | same |
| P13 | Maximum perigee distance constraint | 0.7 normalised length | p1024 (TXT) | constraint value |
| P14 | LEO height, maximum transfer duration, tolerances | 185 km; 6 pi in f; 1e-10 integration; 1e-7 solver | Table 2 (IMG) | setup |
| P15 | Stable-manifold backward time in Fig. 5 | Delta f about 1.82 pi | Fig. 5 caption (IMG) | figure setup |

Digits flagged unsure:
- P7 lambda2 = 0.992: the CMDA paper (p300) prints lambda2 = 1.0086 for the same orbit, whose reciprocal is 0.99146, which rounds to
  0.991. COMPUTED from the project (section 6): 0.99149. So 0.992 is either a rounding up or a slip of one in the third decimal;
  a test should use a tolerance of 1.5e-3 on this value, or use the CMDA 1.0086 and its reciprocal.
- P7 lambda1 = 6.481e-7: the reciprocal of the CMDA 1.5427e6 is 6.4821e-7; the project's 15-digit run gives 1/lambda1 = 6.4821e-7 from
  the eigenvalue 1.54273e6 and 6.480e-7 to 6.486e-7 from the smallest eigenvalue directly (the smallest eigenvalue of a matrix with
  entries near 1e6 carries a relative error of about 1e-3 for double precision); the printed value is within 2e-4 of the
  reciprocal. A test should compare 1 / (largest eigenvalue) at 4 digits, not the smallest eigenvalue.
- The Fig. 5 perturbation (caption "100 km", Table 2 "100 m"): contradictory; neither is a test value.
- P14 "Maximum transfer duration 6 pi": the symbol rendered as "Df max 6p" in the text layer; the image reads "6 pi".
- The Fig. 12 sigma values are read from axis ticks; not a test value.

## 6. Reconciliation against the project and the closure check (run 2026-10-04)

### 6.1 The state against the project's existing M5N2 row

`tests/core/test_er3bp_neelakantan_2022.py` holds the 15-digit row "M5N2 halo": x0 = 0.851666641652152, z0 = 0.183285539178136,
ydot0 = 0.25828972225268, period 4 pi, mu = 0.0122, e = 0.0554, source N&R 2022 Table 8. COMPUTED differences from this paper's
Table 2 state (0.8516666, 0.1832855, 0.2582897): 4.2e-8, 3.9e-8 and 2.2e-8. All three are the 15-digit values rounded (not
truncated) to seven decimals. The two papers therefore publish the same orbit, and this paper is a second, independent source for it.
The project's test measured closure 2.9e-7 for the 15-digit row (that test integrates with DOP853, rtol = atol = 1e-13). My run below
gives 1.6e-8 with `propagate_er3bp`; the difference is the integration path and tolerance, both at the level the 1.5e6 multiplier
amplifies (section 6.3).

### 6.2 Runs (script and log in the session scratch directory; no test file written)

`core.er3bp.propagate_er3bp`, DOP853, rtol = atol = 1e-13, mu = 0.0122, e = 0.0554, the monodromy as the product of segment STMs
(the authors' Eq. 22 of the CMDA paper). All numbers COMPUTED.

| Case | Result |
|---|---|
| 15-digit N&R state, periapsis (f0 = 0), one 4 pi step | closure ||X(4 pi) - X0|| = 1.6e-8; half-period residuals at f = 2 pi (y, x', z') = (-2.7e-10, -3.0e-10, 1.4e-9) |
| same, apoapsis start (f0 = pi) | closure 1.18 (does not close), as the apoapsis counterpart is a different orbit |
| 15-digit state, monodromy by 1, 4, 8 and 20 segments | largest eigenvalue 1.54273e6 for all four; second real eigenvalue 1.00858 for all four; unit-circle pair 0.949557 +- 0.31359i for all four; reciprocal pair partner of the second real eigenvalue 0.99149 for all four; smallest eigenvalue 6.4800e-7 to 6.4860e-7 (precision-limited) |
| paper's own 7-digit state (P5), periapsis | closure 3.7e-2; half-period residuals (y, x', z') = (-3.5e-5, -3.8e-5, 1.8e-4); monodromy largest eigenvalue 1.42665e6 and a different unit pair, so the rounded state is not on the periodic orbit to the precision the monodromy needs |
| single-shooting re-correction of the 7-digit state (three unknowns x0, z0, ydot0; three residuals y, x', z' at f = 2 pi; scipy fsolve) | converges (17 residual evaluations) to residual < 1.5e-13; the corrected state differs from the printed one by (4.2e-8, 3.9e-8, 2.2e-8) and from the N&R 15-digit state by (1.0e-13, 2.2e-13, -1.2e-13), i.e. the printed seven digits plus a three-equation solve recover the 15-digit published orbit to 2e-13; closure then 2.7e-7; monodromy 1.54273e6, 1.009764 (the second real eigenvalue differs from the N&R state's 1.00858 at the 3rd decimal: this pair is nearly neutral, so its value is sensitive to a 1e-7 state shift or to the segmentation path), unit pair 0.949557 +- 0.313595i, partner 0.99033 |
| identity: state propagated over f in [pi, 3 pi] with e, versus over f in [0, 2 pi] with -e | agree to 5.9e-15 (the equations depend on f only through e cos f) |

Printed against computed (15-digit state):

| Quantity | Printed | Computed | Agreement |
|---|---|---|---|
| largest eigenvalue | 1.5427e6 (CMDA p300) | 1.54273e6 | 4 digits |
| second real eigenvalue | 1.0086 (CMDA p300) | 1.00858 | 4 digits |
| its reciprocal | 0.992 (ASR p1020) | 0.99149 | 0.992 is off by 5e-4 (see P7) |
| 1 / largest | 6.481e-7 (ASR p1020) | 6.4821e-7 | 2e-4 relative |

The unit-circle pair (0.94956 +- 0.31359i, modulus 1.0000) is not printed in either paper for this orbit; it is a project prediction, not a
control.

### 6.3 Why the printed seven digits do not close

COMPUTED: a state error of about 4e-8 in each of three components, multiplied by an unstable multiplier of 1.5e6 over the 4 pi interval
gives 6e-2 as an upper bound (4e-8 x 1.5e6). The measured closure is 3.7e-2, so the rounding of
the printed state alone explains the miss. This is the same amplification mechanism the CMDA paper warns about (p290: eigenvalues "differ
more than 20 orders" over longer times). A test built on this paper's state must assert the half-period perpendicular-crossing residual
with a tolerance of about 5e-4 (measured 1.8e-4 maximum component), or assert agreement of the re-corrected state with the printed one to
1e-7, not a closure to 1e-6.

### 6.4 Code and constants

- `core/er3bp.py`: `er3bp_eom` agrees with Eqs. 2 to 4 (the paper uses the same forms as the CMDA paper). `propagate_er3bp(with_stm=True)`
  supplies the segment STMs; the first four lines of the table above exercise it.
- `search/er3bp_floquet.py:er3bp_monodromy` is a single full-period STM from f = 0; this run showed that a single 4 pi step and a 20-segment
  product agree to 5 digits on the largest eigenvalue at tolerance 1e-13 for this orbit, so the segmented product matters only at looser
  tolerance or longer periods (the companion digest's estimate of a 1e-6 relative error "at the project's default tolerance" was not run
  here; the project's default is rtol = atol = 1e-12).
- No ER3BP stable-manifold, bridge-segment or perigee-transfer code exists (grep of `core`, `genome` and `search` for "manifold" finds only
  CCR4BP, BCR4BP, QBCP, a CR3BP asymmetric-branch module, a bicircular transfer module (`genome/bct_transfer.py`) and torus modules). The
  transfer half of this paper has nothing to reconcile against.
- Constants: the project's elliptic tests use mu = 0.0122 and e = 0.0554, the papers' values. `search/er3bp_discovery.py` defaults to
  target e = 0.0549 (the mean lunar eccentricity). `core/bcr4bp.py` and `core/qbcp.py` use mu = 0.0121505816 (Andreu-style). The two
  mass ratios differ by 0.4 percent and the two eccentricities by 0.9 percent. With a monodromy of 1.5e6, a printed-value check must use the
  paper's mu and e exactly; the companion digest's Table 2 (mu = 0.015) shows lambda1 changing by about 1.5 percent per 0.01 in e, so a
  0.0005 shift in e changes lambda1 by about 0.08 percent, which is larger than the printed precision (4 digits). INFERRED, not run.
- Fig. 2 and Fig. 3 are plotted in km and in the pulsating frame with no numbers; nothing to reconcile.

## 7. Printed orbit versus the project's elliptic code: what it can and cannot settle

See section 8 for the ledger items. Summary: the state is a periapsis-group orbit, so `core.er3bp` from f = 0 is the right counterpart. The
apoapsis group is not printed here. The paper's own generation method (multi-segment, interior point) is not described here, only named.

## 8. Do the printed orbits resolve `#912` and `#925`?

### `#912` (the elliptic-problem corrector lacks the published method)

Read first: `data/OUTSTANDING.md` `#912` lists the missing pieces as an apoapsis starting anomaly, an elliptic multi-segment corrector,
the segmented monodromy with stability indices and eigenvalue types, and three-dimensional continuation.

What this paper adds, from the runs in section 6:
1. A positive control for a three-dimensional elliptic periodic orbit exists and works. The printed Table 2 state, re-corrected by a
   plain single-shooting solve in `core.er3bp` (three unknowns, three half-period perpendicular-crossing residuals, f0 = 0, free x0, z0,
   ydot0, T_E = 4 pi fixed), converged to a residual of 1e-13 in 17 residual evaluations, moved the state by about 4e-8, and recovered the 15-digit N&R state to 2e-13. So the
   claim in the companion digest that "single-segment correction" was a gap for M5N2 (N = 2, M = 5) does not hold for this orbit at
   least when the seed is within 1e-7 of the solution. The gap `#912` names is real for seeds far from the solution (the continuation
   in e and mu with natural-parameter steps of 0.001 from a circular halo) and for larger N; this paper does not test that, because it
   does not describe the continuation.
2. The segmented monodromy product and the single-step monodromy agree to five digits on this orbit (4, 8 and 20 segments against one), at the
   project's standard tolerance. Not a reason to skip the segmented product in general, but evidence that the missing piece is not
   urgent for N = 2 and multipliers of 1.5e6.
3. The elliptic corrector's free-variable choice `(IDX_X, IDX_Z, IDX_YDOT)` with residuals `(IDX_Y, IDX_XDOT, IDX_ZDOT)` is the right
   formulation for this orbit class: this run solved exactly that system with an off-the-shelf solver.
4. It does not supply anything for the apoapsis group (no state printed), three-dimensional continuation, the stability index or the
   eigenvalue-type classification. The companion paper's tables remain the only eigenvalue sources.

### `#925` (elliptic-problem controls that do not reproduce: model or paper?)

Read first: `#925` item 1 (Mako & Salamon 2025 unstable bands) and item 2 (N&R Table 8 M4N2 Lyapunov does not close).

- Item 2 (M4N2 Lyapunov): NOT resolved. This paper prints no planar orbit, no M4N2 orbit and no Lyapunov orbit. It says only that a
  M2N1 left orbit at L2 exists and has two real pairs. It does not bear on whether the M4N2 row of N&R has a slip.
- What this paper does add to the model-or-paper question: the project's `core.er3bp` reproduces, from a published state (two papers
  agree on it to the printed digits), a published eigenvalue of 1.5427e6 and 1.0086 to four digits at e = 0.0554. This tests, at a
  periodic 3D orbit over 4 pi, the equations of motion (Eqs. 2 to 4) AND the variational equations (`er3bp_stm_eom`) in the project, which the
  four closing Table 8 rows of N&R test only for the state. A model defect in the elliptic term (sign of e cos f, the z term of Eq. 4,
  the pulsating-frame factor) would change the multipliers at the first order in e. So it moves the posterior toward "the M4N2 row is
  the problem" and away from "`core.er3bp` is wrong", for orbits of this class. Limits: the orbit is a periodic orbit near the
  L1 libration region; this does not test the regimes of Mako & Salamon's true-anomaly sweeps.
- Item 1 (Mako & Salamon): untouched; their bands are not an ME-Halo property.

## 9. Techniques applicable to the project's problems

Each ledger entry was read before this mapping. The paper is a transfer-design paper for a halo-type orbit; its methods map unevenly. Where
a mapping is weak this is said.

### `#912` (elliptic corrector and stability tools)

- The orbit-classification and parameter conventions: four groups, f0 = 0 or pi, the phase-epoch coupling (a spacecraft on an ME-Halo
  or its manifolds only exists at epochs f* + 2 k pi). A catalogue row for an elliptic orbit must record f0.
- The sensitivity lesson (section 6.3): a published seed rounded to seven digits misses closure by 4e-2 on an orbit with multiplier
  1.5e6. A corrector's acceptance test must scale with the monodromy norm: use the half-period perpendicular residual (here 1.8e-4)
  or a re-correction distance, not a closure threshold of 1e-6.
- The three-unknown, three-residual single-shooting formulation worked for this orbit (section 8); the multi-segment solver is
  needed for continuation robustness, not for closing a well-seeded orbit.

### `#925` (model or paper for elliptic controls)

- The Table 2 state and the eigenvalue match are a new, independent control (section 8). Recommended as a permanent test (follow-up 1).
  The test's expected side is the printed numbers, not project output: the state P5, eigenvalues 1.5427e6 and 1.0086 (CMDA p300).
- Not a method, but a template: when a published orbit does not close, first re-correct and report the distance moved against the
  multiplier-amplified rounding (section 6.3). For M4N2 this would turn "does not close" into "closest closed orbit is at distance d"
  (the N&R digest section 6(e) item 4 already asks for this).

### `#902` (bicircular validation tiers on a real orbit)

- This paper's setting is a time-periodic model with a phase variable (true anomaly) that must be matched at injection: the manifold
  depends on (tau, beta, sigma) and the epoch f* + 2 k pi. The bicircular model has the same structure with the Sun angle replacing the
  true anomaly. Lesson for the validation tiers: a real, published orbit must be paired with its published phase, and the tier tests
  must pass the phase as an explicit argument; the Table 2 state is a worked example of a complete (state, phase, period) triple
  (periapsis, f0 = 0, 4 pi) that closes to the published eigenvalue.
- The positive-control discipline: close a published orbit, then reproduce a published number that depends on the variational
  equations (a multiplier), as done in section 6. A bicircular orbit that only closes does not yet test the model's partials.
- Weak: the paper has no bicircular or Sun content.

### `#905` (rerun `#884` with bifurcation detection in the family walk)

- The eleven-type eigenvalue bookkeeping and the modified stability index (CMDA digest, section 4.1) are the right frame for "bifurcation
  detection in the family walk": a collision or a swap of eigenvalue type between two continuation steps is the signal.
- The tables' column mislabelling (CMDA digest, printed anomaly 2) is the failure to avoid: track eigenvalues as sets across steps.
- The segmented monodromy product is an accuracy tool for long-period family members (here the single step and the 20-segment product
  agree, so it is optional at 4 pi).
- The paper itself does not address the C32 and C31 members at 5/2 or Leiva and Briozzo's arcs. Nothing here settles that question.

### `#913` to `#915` (real-ephemeris Earth-Mars cycler reproduction, moon-cycler flybys, Russell-Ocampo optimiser)

- The transfer half: backward integration from the target orbit to a free junction state, with the junction optimised against a cost
  and an initial guess from a two-body (Hohmann) transfer, then a tangential-arrival condition at the departure orbit. This is the
  same architecture as a leg-by-leg cycler optimiser: the free junction (the MI point) is where the cost lives, and the paper found the
  best junction at the manifold apogee, not at the intuitive perigee (3.4 to 3.9 km/s against 3.8 to 4.7 km/s). For `#915` the
  analogue is: do not fix the flyby geometry at a "natural" point before optimising it.
- A survey-hygiene lesson: the discontinuity at tau = 0.31 to 0.32 in Fig. 8 is an artifact of a hard constraint (perigee distance
  0.7), diagnosed by re-running that region and plotting the trajectories (Fig. 9). The same check belongs in any `#913` to `#915`
  sweep: when a cost curve jumps, examine whether a constraint switched branch before reporting a physical structure. (Compare the
  project's memory note on isolated sweep flips.)
- No ephemeris, no planetary flyby and no Mars content: the connection is architecture only. `#914` (integrated-flyby moon cycler) could
  borrow the bridge segment and tangential-arrival idea; this paper gives no flyby treatment (its internal transfers have a lunar flyby and the authors
  avoid them as "very chaotic").

### `#916` (printed persistence conjecture, stable near-commensurate cyclers under eccentricity)

- The companion tables show which near-commensurate orbits stay centre-type: the apoapsis group's near-1 pair stays on the unit circle
  for small e (CMDA Tables 4 and 5) while the periapsis group's becomes real. This paper's target (periapsis, lambda2 = 1.0086 real) is
  therefore the saddle counterpart: it has no centre direction beyond the complex unit pair; relevant to which counterpart a
  persistence claim concerns (the companion digest's section 6(c)).
- This paper adds the observation that a second real stable eigenvalue (0.992, nearly neutral) creates a 3-D stable manifold, so near
  such an orbit the dynamics has a nearly neutral stable direction: a very slow dissipative-looking drift, which is a nuisance for an
  invariant-curve persistence test (convergence of an invariant-curve solver may stall along it). INFERRED.

### `#917` (Casoliva rows in the elliptic problem with the fold count)

- The ME-Halo target's two-real-pair structure shows what the elliptic stability of a resonant orbit looks like in practice at the
  Earth-Moon eccentricity: three-dimensional manifolds, multipliers of 1e6, a counterpart-dependent real or unit pair. For a Casoliva
  row carried to e = 0.0554, the first diagnostic is the near-1 pair's behaviour (companion digest section 6(c)).
- The 2-D versus 3-D stable-manifold dimension count is a classification of the outcome for the fold-count study.
- Blocked on `#912` as registered; the closed Table 2 control makes the corrector side testable now for a 3-D periodic orbit.

### `#918` (global dynamic-programming search for zero-radius real-ephemeris cyclers)

- Not applicable. No technique in this paper maps; the survey over (tau, direction) is an exhaustive scan of a three-parameter
  manifold with local optimisation of one junction parameter, not a graph search.

## 10. What the paper leaves open, and gaps found

- Does not print the continuation method, the multi-segment details, an apoapsis orbit state, any M other than 5 for the target, any
  circular-halo state, the unit-circle eigenvalue of the target, any transfer's maneuver vectors or epoch, or any table of transfers.
- "However, it is hard to make a precise and direct comparison since the target ME-Halo orbit only exist in the ERTBP" (p1025): no
  benchmark against a circular-problem transfer.
- Authors' stated limits: internal transfers are too chaotic to survey; no global optimum along the manifold is sought; the redundant
  direction raises computing time and local minima.
- Slips found: "CRTPB" for CRTBP (p1016, p1018, p1020), "ERTPB" (p1020), "Df max" and "Dt" symbol rendering in the text layer; the Fig. 5 caption
  (100 km) against Table 2 (100 m); "respectively when primaries are at their periapsis and apoapsis" (true for N odd only);
  lambda2 = 0.992 against 1.0086's reciprocal 0.9915; "Earth-Moon period 27.5 days".

## 11. Recommended follow-ups (not registered)

1. A permanent test, `tests/core/test_er3bp_peng_xu_2015_asr.py`, with expected values from this paper and the CMDA paper only: (a) the
   Table 2 state re-corrected by single shooting moves by less than 1e-7 and closes the half-period perpendicular crossing; (b) the
   largest monodromy eigenvalue within 1e-3 relative of 1.5427e6 and the second real eigenvalue within 5e-4 of 1.0086; (c) the printed
   7-digit state's half-period residual under 5e-4 (not a closure test); (d) the apoapsis start does not close (control); (e) the
   identity f0 = pi with e against f0 = 0 with -e at 1e-12. Not written here (instruction: no new test files).
2. Compute the apoapsis counterpart of the same orbit by re-correcting from f0 = pi (or e = -0.0554 at f0 = 0) and check the CMDA p300 numbers
   (one real pair about 1.5431e6, two unit-circle pairs); this closes the CMDA digest's C5.
3. Run the CMDA digest's C1, C2 and C6 now that the orbit solve is shown to work (circular limits at mu = 0.009 and 0.015 and the
   e^(5/2) splitting law); a three-unknown shooting solve at e = 0.02 to 0.10 with natural continuation from the circular halo is within
   what was run here.
4. Re-correct the N&R M4N2 Lyapunov row by the same single-shooting solve from both apses and report the distance moved and the nearest
   closed orbit (the N&R digest's item 4); this is the direct `#925` item 2 test and costs minutes.
5. Add `f0` (or a signed e) as an argument to `genome/er3bp_periodic.correct_er3bp_periodic` and to the planar and 3-D continuation drivers
   (the `#912` first missing piece); the identity at 5.9e-15 shows the negative-e route is already exact.
6. If a transfer capability is ever wanted: an ER3BP stable-manifold generator at a fixed (tau, f0) with the two-stable-eigenvalue case
   handled (the existing manifold code is CCR4BP/BCR4BP/QBCP and does not assume this), tested against this paper's printed
   ranges (3.388 to 3.934 km/s) with the paper's perigee-distance constraint 0.7.
7. Check how the project reports the unit-circle pair (0.94956 +- 0.31359i) against `monodromy_eigenstructure`, which the companion digest
   predicts would raise on the Table 5 orbits; the Table 2 orbit has one complex unit pair, so it should pass, and is a cheap positive control
   for that function.
8. Acquire the 24th ISSFD paper "Numerical stability study of multi-circle elliptic halo orbit in the elliptic restricted three-body
   problem" (Peng and Xu 2014), which may print the multi-segment method and orbit states; and Campagnola 2010 (the origin of the
   classification).

## 12. References printed (journal pp1026-1027)

Cited and relevant: Alessi, Gomez & Masdemont 2010 (Adv. Space Res. 45:1276, two-manoeuvre transfers LEO to Lissajous orbits; DOI
10.1016/j.asr.2009.12.010); Broucke 1969 (AIAA J. 7:1003); Campagnola 2010 (USC thesis) and Campagnola, Lo & Newton 2008 (AAS/AIAA
Galveston); Hou & Liu 2011 (MNRAS 415:3552); Howell & Pernicka 1987 (Celest. Mech. 41:107); Hyeraci & Topputo 2010, 2013; Li & Zheng 2010a,b
(Celest. Mech. Dyn. Astron. 108:203; Acta Astronaut. 66:1481); Parker & Anderson 2014 (Wiley, Low-Energy Lunar Trajectory Design);
Peng & Xu 2014 (24th ISSFD); Qi, Xu et al. 2012a,b, 2014; Richardson 1980; Sarris 1989; Szebehely 1967; Topputo, Belbruno & Gidea 2008.
Held status is recorded in the companion digest's reference table; none of the not-held entries was found held in this session beyond what
that table says.
