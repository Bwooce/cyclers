# Digest: Hitzl & Henon 1977, "The stability of second species periodic orbits in the restricted problem (mu = 0)"

Acta Astronautica 4:1019-1039 (1977), DOI 10.1016/0094-5765(77)90004-2, received 29 November 1976, presented at the XXVIIth IAF Congress, Anaheim, 10-16 October 1977.
Filed in the private paper corpus as
`hitzl-henon-1977b-stability-second-species-periodic-orbits-restricted-problem-mu-0-acta-astronaut-4-1019-doi-10.1016-0094-5765-77-90004-2.pdf`
(21 PDF pages; PDF page n is journal page 1018 + n). Every equation, table digit and caption below was read from the 150 and 250 dpi page images (the one table, Table 1, is printed rotated on p.1034 and was read from a rotated 250 dpi crop).
Evidence tags: READ (page) = read on the page image; COMPUTED = my own arithmetic on 2026-10-05 (60-digit and 30-digit mpmath, scratch scripts not kept in the repository); DERIVED = algebra from READ equations; INFERRED = my reading across sources.
Companion digest, which this paper completes ("paper I"): `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md` (Celest. Mech. 15:421). Related: `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`, `docs/notes/2026-10-04-digest-henon-2001-generating-families-II-part-a-type-1.md` and `...-part-b-type-2.md`, `docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md`, `docs/notes/2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-a.md`.

## 0. The paper in eight lines

1. Question: how stable are the second species periodic orbits of the planar restricted problem for small mass ratio mu, and in the limit mu = 0 which generating orbits sit on the stability boundary (READ p.1019-1022).
2. Method: first-order matched asymptotics of Breakwell and Perko (1966). The orbit is a Keplerian arc from one lunar departure to the next, plus a hyperbolic flyby that the Moon-centred hyperbola turns into a 4x4 "recurrence matrix" M with the property that one entry carries a factor 1/mu (READ eqs. 1-19, pp.1024-1028).
3. Stability index (READ eq. 20, p.1029): k = (Tr - 2)/2 with Tr the trace of the 4x4 M; stable means 0 <= Tr <= 4, that is -1 <= k <= 1; critical orbits (stability boundary) have k = +1 or -1.
4. Scaling result (READ eqs. 18, 19, 23, p.1028-1029): Tr = sum T_ii + h1 T34 + (h2/mu) T24, with every T_ij of order 1. Perturbations are amplified by order 1/mu per period, "which explains why second species solutions are in general strongly unstable and become infinitely unstable for mu = 0". Stability at small mu requires the single element T24 = d b1 / d alpha_-1 to be small, of order mu.
5. The necessary condition S = V T24 = 0 (eq. 21), reduced to a closed trigonometric formula in the generating-orbit variables (tau, eta) (eq. 40), is, to within multiplicative factors, identical to the earlier sufficient criticality condition G* = 0 of paper I (eq. 42). Hence the headline (box, p.1038): second species orbits are critical in the limit mu = 0 if and only if C attains an extremal value.
6. At such an orbit the two critical orbits k = +1 and k = -1 "coalesce into one for mu tending to 0", k jumps from large positive to large negative values (Fig. 2).
7. Extra: vertical (out of plane) criticality (pp.1033, 1038): the only second species generating orbits critical both horizontally and vertically are the retrograde circles A_i(0).
8. Table 1 (p.1034) gives seven critical orbits with a, e, x0, x1, C, T (Figs. 6-11 and the limiting orbit m1, m2 = A0(-1)); no orbit at mu > 0 is computed; stability at mu > 0 is deferred to two later papers (p.1039).

## 1. Setting and notation (READ pp.1019-1023)

- Intro (p.1019): "Analysis of the long-term behavior of orbital motion ...". Close-encounter periodic orbits "are highly unstable" in numerical work; "an early theoretical stability analysis (Abraham, 1967, p. 229) was a failure". The Breakwell-Perko second-order matching (1974) was developed "to obtain a stability analysis for these orbits valid through O(mu)" (p.450 of that paper); what was missing was a supply of correct generating solutions. Three possible families have been described previously (Henon and Guyot, 1970); the simplest is Fig. 1.
- Fig. 1 (p.1020): critical orbits of the retrograde family m that enclose both finite masses, drawn at mu = 0, 0.1, 0.2, 0.3, 0.4, 0.5 (Earth M1, Moon M2). At mu = 0 the orbits m1 and m2 coincide (the limiting orbit with one angular point at the Moon, called "important" on p.1021) and m3 is the large circle-like orbit; m1, m2, m3 separate with mu; m2 and m3 are absent at mu = 0.4 and 0.5 (only m1 drawn). Caption, READ p.1020: "These orbits are stable between m1 and m2 and outside m3 for 0 < mu < mu0 = 0.327... and for 1 - mu0 < mu < 1; for mu0 < mu < 1 - mu0, the orbits are stable outside m1. Otherwise, the orbits are unstable. (From Henon and Guyot, 1970, p. 354)." (The caption wording is partly ambiguous in the rotated print; the numeric mu0 = 0.327... is clear.) The text (p.1021) says the limit orbit m1, m2 "provided the motivation for this entire analysis".
- Definitions (READ p.1021): second species (Poincare) is a Keplerian ellipse about the Earth that at some time passes close to the Moon; as mu tends to 0 these approach arcs of Keplerian ellipses joined at "angular points" (corners) at the Moon, called consecutive collision orbits (Henon 1968) or generating orbits (Perko 1974). Critical means stability index k = +-1 (Henon and Guyot 1970; Hitzl 1975a).
- Fig. 2 (READ p.1022): the stability index k against the initial condition x0 along one family at mu > 0, on a log-like vertical axis with ticks 1000, 100, 10, 1, 0, -1, -10, -100, -1000; the stable band |k| <= 1 is shaded; critical orbits are marked where the curve crosses k = +1 and k = -1. Caption: "As mu tends to 0, the maximal values for |k| tend to infinity and the slope dk/dx0 becomes generally steeper. In the limit mu = 0, the entire locus for k is at +-infinity with occasional jumps from one side to the other at isolated critical orbits possessing extremal values for C." (The caption's printed text "the maximal values for |k| -> infinity" is READ as shown.)
- Geometry (READ pp.1022-1023, Figs. 3 and 4): Earth M1 at the origin of an inertial frame, Moon M2 on the unit circle, period 2 pi (n = 1), direct; particle M3 on a coplanar ellipse with perigee S' or apogee S'' at t = 0 when the Moon is at R; Fig. 3 is the generating orbit with collisions at P (t = -tau) and Q (t = +tau), angle tau marked each side of R. Relative to the Moon the passage is a hyperbola (Fig. 4: asymptotes, impact parameter b, angle alpha of the relative velocity, deflection delta). For mu tending to 0 the periapsis r_p tends to 0 while delta stays finite, so the hyperbola degenerates to an intersection with the Moon, an angular point.
- Four times (READ p.1024, Fig. 5 sketch of m1 in inertial axes): t_-1 previous departure from the Moon's orbit; t_0 = 0 perigee passage near the Earth; t_1 = -t_-1 next encounter; t_2 next departure; t_2 tends to t_1 as mu tends to 0, but the four are kept distinct. A complete orbit is t_-1 <= t <= t_2.

## 2. The stability matrix M, equation by equation (READ pp.1024-1029)

### 2.1 Eq. 1 (p.1024): propagation from one lunar departure to the next arrival

First-order asymptotic matching (Breakwell and Perko 1966, p.175 eq. 4) relates the hyperbolic quantities at arrival (t_1) and at the previous departure (t_-1) through the state transition matrix Phi(t_1, t_-1) of the Keplerian arc, plus an integral of the "gross" bias terms (g, g'):
- Left side at t_1 (printed as a two-row block): the position part is {(mu/V_1^2)[2 - ln(2 V_1^3 (t_1 - t_-1)/(mu e_H))] + V_1 (t_1 - t_p1)} i_1 + b_1 j_1 and the velocity part is v_inf1 - V_1 (1 - mu/(V_1^3 (t_1 - t_-1))). The right side is Phi(t_1, t_-1) applied to the same block built at t_-1 (with ln(...) - 2 and t_-1 - t_p-1 signs, as printed), plus the integral of (g, g') dt from t_-1 to t_1. The print of the log arguments is clear; the exact signs inside the braces are as printed and are not needed below.
- Symbols (READ p.1024): V_1 = relative arrival velocity at the massless Moon on the unperturbed conic at t_1; e_H = eccentricity of the hyperbolic passage (e_H > 1); t_p = time of perigee passage of the osculating hyperbola; b = perpendicular distance from the Moon to the asymptote (impact parameter); v_inf = relative velocity at infinity on the osculating Moon-centred hyperbola; i_1 = direction of V_1; j_1 = k x i_1 normal to V_1 (direction of the impact parameter, b_1 = b_1 j_1).
- Eq. 2 (p.1025): [g; g'] = mu [Phi_rv(t_1, t) ; Phi_vv(t_1, t)] f[r_0(t), t] - (mu V_1/V_1^3)[(t_1 - t)^-1; (t_1 - t)^-2] + (mu V_-1/V_-1^3) Phi(t_1, t_-1) [(t_-1 - t)^-1; (t_-1 - t)^-2], with f (eq. 3) the restricted-problem perturbation function f = -[(r - r_M)/|r - r_M|^3 + r_M/|r_M|^3 - r/|r|^3] and the inertial position expanded as r = r_0 + mu rho_1 + mu^2 rho_2 + ... (eq. 4); r_0 is the unperturbed ellipse through the massless Moon at t_1 and t_-1. Bounded as t tends to t_1 and t_-1 (the singular parts are subtracted).

### 2.2 Perturbation vector and its linearisation (READ eqs. 5-10, p.1026)

- Perturbation vector (eq. 5): P = (dt_p, db, dv_inf, d alpha) (time of periapsis, impact parameter, speed at infinity, direction angle alpha of the relative velocity; Fig. 4).
- Reference orbit symmetric: the relative speeds at t_1 and t_-1 are equal, V_1 = V_-1 = V = |V| (eq. 6). Perturbation only in the direction of V is allowed (footnote "see second comment at the end of this section"): dV = V d alpha j (eq. 7).
- Differentiating eq. 1 (eq. 8): the "gross" biases remain effectively constant as the previous departure conditions are varied. The 1/mu terms come from d e_H.
- Hyperbolic relations (eq. 9): r_p v_p = b v_inf; e_H^2 = 1 + v_inf^4 b^2 / mu^2; delta = alpha_2 - alpha_1 = 2 arcsin(1/e_H). Hence (eq. 10) d e_H = (1/(2 e_H)) d e_H^2 = (V^3/e_H)[2 (b/mu)^2 dv_inf + V (b/mu^2) db] (printed as shown; the exponents are as READ).

### 2.3 The 4x4 recurrence, eqs. 11-14 (READ pp.1027)

P_1 = T P_-1 (eq. 11) with the sixteen elements T_ij (eq. 12), written with the 2x2 blocks of the 4x4 Keplerian state transition matrix Phi = [Phi_rr Phi_rv; Phi_vr Phi_vv] (eq. 13) and Q = (b/mu) V / e_H^2 (eq. 14). Complete list as READ (all T_ij are O(1); terms of O(mu) dropped):

| | col 1 | col 2 | col 3 | col 4 |
|---|---|---|---|---|
| row 1 | i1^T Phi_rr i_-1 - V Q j1^T Phi_rr i_-1 | -(1/V) i1^T Phi_rr j_-1 + Q [i1^T Phi_rr i_-1 + j1^T Phi_rr j_-1] - V Q^2 j1^T Phi_rr i_-1 | -(1/V) i1^T Phi_rv i_-1 + Q j1^T Phi_rv i_-1 | i1^T Phi_rv j_-1 - V Q j1^T Phi_rv j_-1 |
| row 2 | -V j1^T Phi_rr i_-1 | j1^T Phi_rr j_-1 - V Q j1^T Phi_rr i_-1 | j1^T Phi_rv i_-1 | -V j1^T Phi_rv j_-1 |
| row 3 | -V i1^T Phi_vr i_-1 | i1^T Phi_vr j_-1 - V Q i1^T Phi_vr i_-1 | i1^T Phi_vv i_-1 | -V i1^T Phi_vv j_-1 |
| row 4 | j1^T Phi_vr i_-1 | -(1/V) j1^T Phi_vr j_-1 + Q j1^T Phi_vr i_-1 | -(1/V) j1^T Phi_vv i_-1 | j1^T Phi_vv j_-1 |

(Transcription caution: T_12 in the printed row 1 has three groups; my table keeps them in the printed order. T_22 prints "j1^T Phi_rr j_-1 - V Q j1^T Phi_rr i_-1" and T_32 prints "i1^T Phi_vr j_-1 - V Q i1^T Phi_vr i_-1"; the subscripts on Phi in the last term of T_22 and T_32 are as READ.) Row 2 gives db_1, row 3 dv_inf1, row 4 d alpha_1, row 1 dt_p1.

### 2.4 One full period: M (READ eqs. 15-19, pp.1028)

- Symmetric orbits (eq. 15): b_1 = b_-1 = b; v_inf1 = v_inf-1 = v_inf = V + O(mu); alpha_2 = -alpha_1.
- Differentiating eq. 9c for the hyperbolic turn (eq. 16): d alpha_2 = d alpha_1 - [4 sqrt(e_H^2 - 1)/(V e_H^2)] dv_inf - [2 V^2/(mu e_H^2)] db_1 = d alpha_1 + h_1 dv_inf1 + (h_2/mu) db_1. So h_1 = -4 sqrt(e_H^2 - 1)/(V e_H^2) and h_2 = -2 V^2/e_H^2 (DERIVED from the READ equation by matching terms; the paper defines h_1, h_2 only through this equation).
- Recurrence from one departure to the next (eq. 17): P_2 = M P_-1. "It is interesting that eqn (16) affects only the fourth row of T" (eq. 18): M_ij = T_ij for i = 1 to 3, and M_4j = T_4j + h_1 T_3j + (h_2/mu) T_2j, j = 1 to 4.
- Physical reading (READ p.1028): "A slight perturbation Delta b_1 in the impact parameter b of the incoming trajectory produces then a much larger perturbation Delta alpha_2 of the outgoing direction; in fact, the ratio Delta alpha_2/Delta b_1 is of order mu^-1. Thus, the close approach has the effect of strongly amplifying the perturbations ... after the encounter, the perturbation is essentially in the direction; conversely, it is essentially the perturbation in the impact parameter before the encounter which matters."
- Trace (eq. 19): Tr = sum M_ii = sum T_ii + h_1 T_34 + (h_2/mu) T_24, with, for stability, 0 <= Tr <= 4. Index (eq. 20): k = (Tr - 2)/2 (Hitzl 1975a p.154), so critical orbits have k = +-1.
- Determinant (eq. 22): |M| = |T| = 1 (volume preservation, Hitzl 1975c p.193).
- INFERRED (the paper does not say it): the 4x4 M contains the two trivial unit eigenvalues of a Hamiltonian flow with an integral (the time-shift and Jacobi-constant directions), so Tr = 2 + (lambda + 1/lambda) for the one non-trivial reciprocal pair, and k = (lambda + 1/lambda)/2. This is the only reading consistent with the paper's own statement that stability means 0 <= Tr <= 4 together with k = (Tr - 2)/2 and |k| <= 1. It is also Henon's z = (a + d)/2 of the 1997 book (companion digest section 7).

## 3. The stability condition S = 0 (READ pp.1029-1032)

### 3.1 Why S = V T_24 (eqs. 21, 23, p.1029)

- As mu tends to 0, Tr is finite only if the one 1/mu coefficient vanishes: T_24 = M_24 = d b_1 / d alpha_-1 must vanish. Define (eq. 21) S = V T_24 = -V^2 j_1^T Phi_rv j_-1 = 0, "a necessary stability condition for second species orbits in the limit mu = 0".
- Interpretation (READ p.1029): velocity perturbations perpendicular to the relative velocity at departure cause no position perturbation perpendicular to the relative velocity at the next encounter; equivalently perturbations in the direction of the departure velocity have no effect on the impact parameter b at the next encounter.
- Instability scaling (eq. 23): over one full period from just before a close approach to just before the next, Delta b_1/Delta b_-2 = (Delta b_1/Delta alpha_-1)(Delta alpha_-1/Delta b_-2) = (Delta b_1/Delta alpha_-1)(Delta alpha_2/Delta b_1) = Delta alpha_2/Delta alpha_-1, "which is of order mu^-1" because the Earth arc is practically independent of mu so Delta b_1/Delta alpha_-1 is of order unity, and Delta alpha_2/Delta b_1 is of order 1/mu.
- Mechanism of the critical orbit (READ p.1029): along a family Delta b_1/Delta alpha_-1 varies continuously and "may happen ... at some particular orbit this quantity goes through zero and changes sign. Consequently, the stability index k quickly jumps from large positive values to large negative values, or vice versa, and therefore there exist two critical orbits with k = +1 and k = -1 (see Fig. 2); these two orbits coalesce into one for mu tending to 0."

### 3.2 Reduction of S to (tau, eta) (READ eqs. 24-40, pp.1030-1032)

- Inertial-axis states of M3 just after the collision at t = -tau (eq. 24, E = -eta): x_-1 = sigma0 a (sigma2 c_eta - e), y_-1 = -sigma0 sigma1 sigma2 a sqrt(1 - e^2) s_eta; switches (eq. 25) sigma0 = +1 if perigee at x > 0, sigma1 = +1 if direct, sigma2 = +1 if M3 at perigee at t = 0 (-1 for the opposite); c_x = cos x, s_x = sin x. Kepler's equation t = a^(3/2)(E - sigma2 e s_E) (eq. 26), 1 = a^(3/2)(1 - sigma2 e c_E) E_dot = sqrt(a) ... as printed "= V(a) E_dot" (eq. 27; paper I eq. 6). Velocity after collision (eq. 28): x_dot_-1 = sigma0 sigma2 sqrt(a) s_eta, y_dot_-1 = sigma0 sigma1 sigma2 sqrt(a(1 - e^2)) c_eta.
- Moon at t = -tau (eq. 29): (c_tau, -s_tau) = (x_-1, y_-1); velocity (-y_-1, x_-1). Relative velocity direction (eq. 30): V j_-1 = -(y_dot_-1 - x_-1) x^ + (x_dot_-1 + y_-1) y^. Arrival at t = +tau (eq. 31): x_1 = x_-1, y_1 = -y_-1, x_dot_1 = -x_dot_-1, y_dot_1 = y_dot_-1; V j_1 = -(y_dot_1 - x_1) x^ + (x_dot_1 + y_1) y^ (eq. 32).
- The needed matrizant Phi_rv (position change at t_1 from velocity change at t_-1) is taken from Danby (1962), coefficients v_ij in Danby's perifocal axes (pericentre on +x', direct) related to ours by x = sigma0 x', y = sigma0 sigma1 y' (eq. 33); (Delta x'_1, Delta y'_1) = [v11 v12; v21 v22](Delta x'_dot_-1, Delta y'_dot_-1) (eq. 34), giving (eq. 35) S = (y_dot_1 - x_1)[v11 (x_-1 - y_dot_-1) + sigma1 v12 (y_-1 + x_dot_-1)] - (x_dot_1 + y_1)[sigma1 v21 (x_-1 - y_dot_-1) + v22 (y_-1 + x_dot_-1)]. Danby's eccentric anomaly E' relates to ours by E = E' for sigma2 = +1 and E = E' + pi for sigma2 = -1 (eq. 36); Danby's symbols in our notation (eq. 37): S0 = -sigma2 sin eta, C0 = sigma2 cos eta, S = sigma2 sin eta, C = sigma2 cos eta, E - E0 = 2 eta. Eq. 38: r = r0 = 1, n = a^(-3/2). Eq. 39 repeats paper I: a = (1 - sigma c_tau c_eta)/s_eta^2, e = (sigma2 c_eta - sigma0 c_tau)/(1 - sigma c_tau c_eta), sigma = sigma0 sigma2.
- Closed result (eq. 40, READ at 300 dpi):

  S = 2 s_2tau - 6 (c_eta - sigma c_tau)^2 tau / s_eta^2 + 2 sigma_eta [ c_eta (2 c_eta^2 + 3 c_tau^2 - c_tau^4) - sigma c_tau (2 + 2 c_eta^2 - c_tau^2 + c_eta^4) ] / ( s_eta^2 sqrt(1 - sigma c_tau c_eta) ),

  with s_2tau = sin 2 tau and sigma_eta = sgn(sin eta). Eq. 41: rho = sigma_eta sqrt(1 - sigma c_tau c_eta) = sqrt(a) s_eta. Eq. 42: -(1/2) sigma rho s_eta^2 S = G* = sigma (rho^4/s_eta^2) G = 0, "to within these multiplicative factors, the two conditions are identical".
- COMPUTED checks of eq. 40 and eq. 42 (30-digit arithmetic):
  - Eq. 42 as an identity: with G* of paper I (eq. 27, companion digest section 2), -(1/2) sigma rho s_eta^2 S(eq. 40) equals G* to all digits at generic (tau, eta, sigma), for example (tau, eta, sigma) = (1, 2.3, -1) and (0.7, 1.9, +1), by algebra: the secular part of G* is rho sigma [3 tau (c_eta - sigma c_tau)^2 - s_2tau s_eta^2] and the algebraic part equals sigma times the bracket of eq. 40. This is exact, not only on F0 = 0.
  - Eq. 40 at the printed critical orbits: S evaluates to (|S|) 1.4e-11 at A0(1), 4.8e-9 at A0(2), 8.9e-9 at A1(-2) and 1.9e-10 at A0(-1) using the 10-decimal tau/pi, eta/pi of Table I of paper I; and, from the 5-decimal tau/pi, eta/pi of this paper's Table 1, to 2.9e-5 (A0(-1)), 1.4e-6 (C23(1)), 2.5e-4 (C24(2)), 5.2e-8 (C35(1)), 2.4e-4 (C36(2)), 9.4e-6 (C37(1)), 1.3e-4 (C37(2)), the size expected from rounding at 5 decimals. My first coding of eq. 40 had c_eta^4 where the print has c_tau^4 in the first bracket term, and S was then nonzero at every row; with the printed c_tau^4 it vanishes. So the printed eq. 40 is consistent with paper I.
- Second comment, p.1032 (READ): "when calculating stability, it is sufficient to consider perturbations which do not change the value of the Jacobi constant C (Darwin, 1911). Since V = sqrt(3 - C), this means that it is sufficient to consider perturbations which do not change the magnitude V of the relative velocity at collision. Thus one has only to consider perturbations in the impact parameter b and in the velocity direction alpha." (That is the justification for keeping only db and d alpha in the entry T_24.) The factors s_eta and rho vanish only at integer (tau/pi, eta/pi) with sigma = (-1)^(i+j), where, as in paper I, G = 0 is meaningless; so no extraneous solutions are introduced (p.1032). The factor V vanishes only in the degenerate case of M3 coinciding with M2 for all time (an entire family of degenerate solutions, p.1033).

### 3.3 Comments (READ p.1033)

- Hierarchy: F0 = 0 restricts to lines (characteristics, families); G* = 0 restricts to points on them.
- The choice of G* is arbitrary. In deriving G* there were secular terms in both tau and eta; inserting F0 = 0 eliminated the secular terms in eta. "Our condition G* = 0 is only one of several - in fact, an infinity of - possible equations, which all intersect F0 = 0 at the critical orbits, but which all behave differently outside of the characteristics." (This is consistent with the identity above, which holds because S was reduced with the same elimination.)
- Direct orbits C_ij(1) occur at maxima of C, retrograde C_ij(2) at minima of C (p.1033).

## 4. The stability scaling result (DERIVED from READ eqs. 16, 18, 19, 20, 21; not printed as a formula by the paper)

Combining eqs. 19-21: k = (Tr - 2)/2 = (h_2/(2 mu)) T_24 + O(1), with T_24 = S/V and h_2 = -2 V^2/e_H^2, so

  k = -(V / (mu e_H^2)) S + O(1),   hence   |lambda_max| of the non-trivial pair = |k| + sqrt(k^2 - 1) is of order 2 V |S| / (mu e_H^2) per close passage.

Consequences:
- Away from S = 0, |k| is of order 1/mu, so the instability multiplier of a single-passage orbit is of order 1/mu (the paper's eq. 23); as mu tends to 0, k tends to infinity.
- Sign: for V > 0, sign(k) = -sign(S) outside the O(1) neighbourhood of S = 0 (INFERRED sign convention: it depends on the paper's orientation of i, j and of Danby's coefficients; must be tested before use). Since S changes sign across a critical generating orbit, k jumps between +large and -large (the paper's Fig. 2 picture).
- Stable window: |k| <= 1 requires |S| of order mu e_H^2 / V or smaller. If S crosses zero linearly along the family (arclength l in the (tau, eta) plane, dS/dl nonzero), the window has width of order mu e_H^2 / (V |dS/dl|) in l, so it shrinks linearly with mu. The O(1) terms of Tr (sum T_ii + h_1 T_34) shift the window centre by an amount of order mu in l and set its exact edges; the leading-order width is therefore an order-of-magnitude estimate only. The paper states only the qualitative content ("very narrow regions of stability", p.1039).
- The paper's own caption for Fig. 2 states the same trend in words: maximal |k| grows and dk/dx0 steepens as mu tends to 0.
- Exact validity: one close approach per period (the eq. 11 recurrence is for one encounter between t_-1 and t_1). For n encounters per period the monodromy is a product of n such factors, so the amplification is of order mu^-n unless the encounters are related by symmetry; the paper does not treat this.
- COMPUTED order-of-magnitude example at the Earth-Moon mass mu = 0.0121529529, with e_H = 1/sin(delta/2) from eq. 9c, delta from the rotating-frame turn recipe of the companion digest (section 5(b)), V = sqrt(3 - C), dS/dl along the characteristic by 1e-8 finite differences on the printed 5-decimal points (so the values near integer points are rough), l in radians of the (tau, eta) plane. The row A0(-1) uses the 10-decimal point of paper I.

| Orbit | V | turn delta (deg) | e_H | dS/dl | mu e_H^2 / (V abs(dS/dl)) |
|---|---|---|---|---|---|
| A0(-1) | 1.8437 | 65.094 | 1.859 | -4.315 | 0.0053 |
| C24(2) | 1.7105 | 21.315 | 5.407 | -6.74 | 0.031 |
| C36(2) | 1.6730 | 15.884 | 7.237 | -5.271 | 0.072 |
| C37(2) | 1.5134 | 13.431 | 8.552 | 3.649 | 0.16 |
| C37(1) | 0.5085 | 4.769 | 24.04 | -0.366 | 38 |
| C35(1) | 0.2289 | 3.216 | 35.64 | -0.0741 | 910 |
| C23(1) | 0.1691 | 4.353 | 26.33 | -0.0404 | 1230 |

  Reading: for the fast retrograde generating orbits (V above 1.5, e_H below 9) the leading-order stable window at the Earth-Moon mass is of order 0.005 to 0.2 rad of family length, and the formula is plausibly in range. For the slow direct orbits (V below 0.6) the leading-order estimate exceeds 1 rad, which means the O(1) terms dominate and the asymptotic description does not constrain stability at this mass; this agrees with the companion digest's finding that the first-order theory is outside its validity for most of Table I at mu = 0.0122 (only 90 of 160 have mu |ln mu| / V^3 below 1). The A0(-1) turn 65.094 degrees agrees with the companion digest's 65.09 degrees, computed independently of this paper (eq. 24-32 geometry versus the companion's Kepler-velocity recipe, same formulae). These are my numbers, not the paper's.

## 5. Critical orbits, three dimensions, A_i(0) (READ pp.1033-1038)

### 5.1 The seven critical orbits of Table 1 and the Figs. 6-11 (READ pp.1033-1037)

Table 1 (p.1034, printed rotated), "Critical second species orbits for mu = 0. Parameters for the 7 orbits shown in Figs. 1 and 6-11". Columns: name, figure number, n*, tau/pi, eta/pi, sigma0, sigma1, sigma2, a, e, x0, x1, C, T. n* is "one half the number of intersections of the orbit with the x axis (rotating frame) during a complete period of the motion" (p.1033). READ (all digits as printed; the signs of sigma_i read from the image; the minus signs are narrow marks, double-checked against paper I's table of the same orbits, agreement in all seven rows):

| Name | Fig. | n* | tau/pi | eta/pi | sigma0 | sigma1 | sigma2 | a | e | x0 | x1 | C | T |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A0(-1) | 1 | 1 | 0.21345 | 0.39333 | - | - | + | 1.41019 | 0.88445 | -0.16294 | 2.65743 | -0.39913 | 1.34113 |
| C23(1) | 6 | 2 | 1.99451 | 2.99244 | - | + | + | 0.76343 | 0.30997 | -0.52679 | 1.00007 | 2.97130 | 12.53191 |
| C24(2) | 7 | 6 | 2.12158 | 3.77226 | - | - | - | 0.69628 | 0.57791 | 1.09867 | -0.29389 | 0.07423 | 13.33030 |
| C35(1) | 8 | 3 | 2.99611 | 4.99402 | + | + | + | 0.71156 | 0.40543 | 0.42307 | -1.00005 | 2.94756 | 18.82511 |
| C36(2) | 9 | 9 | 3.08311 | 5.83932 | + | - | - | 0.66010 | 0.58829 | -1.04843 | 0.27177 | 0.20092 | 19.37177 |
| C37(1) | 10 | 5 | 2.99563 | 6.98822 | + | + | + | 0.56901 | 0.75797 | 0.13772 | -1.00030 | 2.74153 | 18.82208 |
| C37(2) | 11 | 10 | 3.03775 | 6.90030 | + | - | + | 0.58287 | 0.75226 | 0.14440 | -1.02133 | 0.70963 | 19.08677 |

Notes on Table 1:
- Every row agrees, to the printed 5 decimals, with paper I's Table I (10-decimal tau/pi, eta/pi) for the same orbit (C2,3(1) = 1.994 514 662 3, 2.992 442 288 9, signs -++, and so on; the sign string in paper I's appendix is sigma0 sigma1 sigma2), and with paper I's Table II for A0(-1). COMPUTED check by me against the companion appendix on 2026-10-05: no disagreement.
- COMPUTED reproduction from the printed tau/pi, eta/pi and signs by eqs. 39, 21 of paper I (C = 2 sigma s_tau/rho + s_eta^2/rho^2), x0 = sigma0 a (sigma2 - e), x1 = -sigma0 a (sigma2 + e), T = 2 pi tau/pi: every column agrees to within the 5-decimal input rounding. Largest differences: A0(-1) all columns within 2e-5; C24(2), C36(2), C37(2) within 1.1e-4; C35(1) and C37(1) within 1.4e-4; C23(1) a, e, x0 within 6e-4 (this orbit sits within 6e-3 of the integer point (2, 3) where a, e depend steeply on the input digits). A test against the printed row should therefore use tolerance of order 1e-3 unless the 10-decimal paper I values are used as input.
- Resonance reading (INFERRED, from a = (i/j)^(2/3) and T of order 2 pi i): C23(1) has a = 0.763 versus (2/3)^(2/3) = 0.763, so 3 particle revolutions per 2 lunar periods (T = 2 tau, about 4 pi), C35(1) a = 0.712 versus (3/5)^(2/3) = 0.711, C37(1) a = 0.569 versus (3/7)^(2/3) = 0.568; these three are near-resonant orbits with the smallest V in the table.
- The fixed-point mapping to the Barrabes-Gomez (p, q) labelling used in the companion digest is unchanged.

Figures 6-11 (READ, inertial axes, Earth at the centre, particle's orbit through the Moon's point at x = 1.00 marked by a dot, axes +-1.25):
- Fig. 6 (p.1035): C23(1), a direct orbit in fixed axes with three loops (a three-leaf rosette, a loop at each of three points at 0, 120, 240 degrees).
- Fig. 7 (p.1035): C24(2), the second critical orbit, retrograde, a many-times-crossing interlaced curve.
- Fig. 8 (p.1036): C35(1), "the first (direct) critical orbit", five loops. Fig. 9: C36(2), the second (retrograde) critical orbit.
- Fig. 10 (p.1037): C37(1), the first (direct) critical orbit, seven loops with a close Earth passage. Fig. 11: C37(2), second (retrograde); "only half the orbit is drawn for clarity".
- Text (p.1033): these orbits "can also have rather close passages by the Earth (see, for example, Figs. 10 and 11)". The Fig. 10 and 11 close passages are visible near the Earth in the plots (a tight inner loop); r_p(Earth) = a(1 - e) from Table 1: C37(1) 0.138, C37(2) 0.144 (COMPUTED, in Earth-Moon distance units); x0 and x1 are the signed intersections with the x axis.

### 5.2 Three-dimensional critical orbits (READ pp.1033, 1038)

- Reasoning (p.1033): second species orbits are, in general, infinitely unstable vertically too; the transition from +-infinity to -+infinity occurs when an orbit launched from P with a small vertical velocity returns to the xy plane at Q, which needs P, M1, Q aligned. Two cases: (i) P and Q the same point, integer tau/pi and eta/pi (the case tau/pi integer, eta/pi not, with e = 1, does not give vertical critical orbits); (ii) P and Q diametrically opposite, half-integer tau/pi (Henon 1968, pp.390-391, Fig. 10, paragraph c).
- Henon's 1973 conjecture (p.424 of that paper) that "a vertical critical orbit represents the intersection of the family of plane periodic orbits under consideration with a family of three-dimensional periodic orbits" is "verified here". Category (i) corresponds to intersections with a 3D family of orbits all of period 2 pi i, i a positive integer; category (ii) to the 3D family formed by rotation about the y axis of Fig. 3, as P, M1 and Q all lie on that axis.
- Result (p.1038, READ, emphasised in the text): "the only second species orbits which are critical both horizontally and vertically are the orbits A_i(0). These critical orbits A_i(0) are simply retrograde circular orbits travelling in the same plane as the Moon but in the opposite direction." Fig. 12: an inclined retrograde circular orbit belonging to a family emanating from A_i(0) (angle i marked, a ring tilted about a diameter; caption "A three-dimensional retrograde periodic orbit which emanates from the critical second species orbit A_i(0) with i >= 0").
- Cross-check with paper I (INFERRED): A0(0) has a = 1, e = 0, C = -1, T = pi, tau/pi = eta/pi = 0.5 (companion digest Table II).

## 6. Conclusions and future plans (READ pp.1038-1039)

- Box (p.1038): "Second species orbits are critical in the limiting case of mu = 0 iff C attains an extremal value."
- Plans (p.1039): (1) numerically integrate the equations of motion near a few critical generating orbits at small mu > 0 and form the four sensitivities A, B, C, D to get k = AD + BC (Henon and Guyot 1970, p.352); iterate on k to determine critical orbits k = +-1 as a function of mu; "very narrow regions of stability will be established and, especially, the case mu* = 1/82.30 = 0.01215067... corresponding to Earth-Moon will be examined in detail. Of particular interest are the distances of closest approach to both the Moon and the Earth." (2) Use the second-order matched asymptotic equations of Breakwell and Perko (1974) to get analytic approximations for the stability boundaries at small mu and check several boundary points numerically. "Results ... will be presented in two forthcoming papers." (Not held; no later paper is cited in our corpus. A reference search is a follow-up.)
- Reference list (p.1039), for the record, with what each is used for: Abraham (1967) failed analysis; Arenstorf 1963a, 1963b, Perko 1974 (existence); Breakwell and Perko 1966 (first-order matching, eq. 1), 1974 (second order); Danby 1962 (matrizant for Phi_rv); Darwin 1911 (Jacobi-constant-preserving perturbations); Henon 1968 (generating orbits, F0), Henon and Guyot 1970 (critical families, k = AD + BC), Henon 1973 (vertical stability, equal masses); Hitzl 1975a (index k = (Tr - 2)/2), 1975c (volume preservation, p.193); Hitzl and Henon 1976 (paper I, cited "Celes. Mech. (to appear)"); Perko 1972 (computational asymptotic matching). Note: the paper's mu* = 1/82.30 is printed as 0.01215067...; the exact quotient is 0.012150668 (COMPUTED). The Earth-Moon mu used in the companion digest and the project, 0.0121529529, differs from it in the fourth significant digit (relative difference 1.9e-4, COMPUTED); this paper does not say which constants it used, so the cause is not determined here, and a comparison at the stated mu* must use mu* itself.

## 7. How it relates to the companion paper and to Henon 2001

- Paper I (companion): supplies the enumeration (Table I, 179 critical generating orbits), the criticality equation G* = 0 derived from an extremum of C by a Lagrange multiplier. This paper supplies the physical derivation: S = V T_24 from the matrix M, and proves the two coincide through eq. 42. Hence (companion section 2) the claim "an extremum of C along a periodic family is a critical orbit with k = +1 (Henon 1965)" is now independently verified for the generating limit by direct stability analysis: both the "if" (extremum implies critical) and the "only if" (critical implies extremum) hold for mu = 0.
- The companion's INFERRED section 5(b) said that nothing happens to the demanded turn at a C extremum. This paper is consistent with it: S is a smooth function of (tau, eta) and the turn is a smooth function of the same variables; the stability flip is in S, not in the turn. The demanded turn enters the scaling only through e_H = 1/sin(delta/2), which makes the factor 1/e_H^2 = sin^2(delta/2) in the 1/mu coefficient small for the near-resonant (small-turn) orbits (my section 4).
- Henon 2001 (companion digests, part A section 3.6-3.7 and part B): near a first-to-second-species bifurcation the distance of the true orbit from the generating orbit scales as Delta C = O(mu^nu) with 0 < nu < 1 and the fusion at Delta C = O(mu^(1/2)) (all y_i, Delta a_i of order mu^(1/2)); node and antinode arcs fuse. This paper's scaling sits on the other side of that picture: it applies for Delta C large compared with mu^(1/2), where each encounter is a plain flyby with finite turn and e_H = O(1), and then |k| = O(1/mu) away from S = 0. INFERRED: near a bifurcation (the 64 of 160 critical generating orbits with turn below 2 degrees, companion section 5(b)) the same formula would involve e_H of order mu^(-1/2) or larger, since the turn tends to zero there, so that 1/e_H^2 is small and the order-1/mu amplification is correspondingly reduced; this paper does not cover that regime and Henon's 2001 treatment of stability near bifurcation is not in the digests. Henon 2001 part A (digest line "K values ... Hitzl and Henon 1977") lists K values as being in "volume I sect. 8.2.1 and Hitzl and Henon 1977"; this paper contains no K values and no tabulated G1, G2, G3, so that citation must refer to paper I (or to the book), not to this paper (READ: no such quantity appears in pp.1019-1039).
- Henon 1997 (companion digest section 4.4-4.5, 7): critical arc = extremum of C along a family (his eq. 4.80, Table 4.4); critical orbit = z = +1 or -1 where z is half the monodromy trace. This paper's k is the same quantity.

## 8. Reconciliation with the project's stability code (see also `#931`)

### 8.1 Stability indices across sources (conversion table)

| Source | Symbol | Definition | Stable range | Critical |
|---|---|---|---|---|
| This paper (eq. 20) and Henon | k, z | (Tr - 2)/2 of the 4x4 M including the trivial pair; = (lambda + 1/lambda)/2 for the non-trivial pair (INFERRED) | -1 to +1 | +-1 |
| Hadjidemetriou 1975b | b (and nu = -b/2) | factor (lambda^2 + b lambda + 1) | abs(b) < 2 | b = -2 (lambda = 1), +2 (lambda = -1) |
| Broucke 1969 | a1, a2 (his k1, k2) | roots of the reciprocal quartic | one region of seven | abs(k) = 2 |
| Casoliva et al. 2010 / `search/earth_moon_resonant_families.py` | k | lambda + 1/lambda | abs(k) < 2 | abs(k) = 2 |

So Hitzl-Henon k = -b/2 = Casoliva's k / 2 (for the planar problem with one non-trivial pair). Nothing in this paper is in conflict with the previous conversion in the companion digest section 7.

### 8.2 `search/er3bp_floquet.floquet_classify`

READ in the code (2026-10-05): tag "unstable" if max modulus > 1 + tol; "stable" if max modulus < 1 - tol; otherwise "marginal". For a symplectic monodromy the eigenvalues come in reciprocal pairs and the trivial pair (1, 1) is present, so the maximum modulus is never below 1: the "stable" branch is unreachable, a linearly stable orbit lands in "marginal". The docstring ("stable: all |lambda| <= 1 + tol") says the opposite of the code. This paper adds a reason and a datum to the fix:
- A correct classifier should decide by the non-trivial pair's index k = (lambda + 1/lambda)/2 (planar, one pair) or by (alpha, beta, Delta) (Hadjidemetriou 1975b, Broucke 1969): "stable" iff abs(k) < 1 - tol, "unstable" iff abs(k) > 1 + tol, "critical" in between; and it should also report k itself, because the sign of k (positive hyperbolic, lambda > 0, versus inverse hyperbolic, lambda < 0) is what flips at a critical generating orbit, and its size (order 1/mu per close passage) is the diagnostic that distinguishes a second-species orbit from an ordinary unstable one.
- Test inputs from this paper: none at mu > 0 (the paper contains no mu > 0 orbit and no printed monodromy), so the tests must come from Hadjidemetriou 1975b and Broucke 1969 (companion digests). This paper contributes the closed form of the transitional behaviour: along a family, k is proportional to -V S/(mu e_H^2) plus O(1), a synthetic generator for the regime test: k = -c S / mu with S crossing zero produces the sequence unstable (k large) to stable window to unstable (k large of opposite sign).
- The only consumer `search/er3bp_discovery.py` maps stable and marginal alike to "elliptic" (per `#931`), so no stored result is affected. A second INFERRED risk from this paper: a one-sided tolerance on max modulus cannot even resolve a stable window of width order mu e_H^2/(V |dS/dl|) when it is narrow (the 0.005 rad examples of section 4), because nearby samples on a continuation with a larger step jump from k = +1000 to k = -1000 without a sample inside; the classifier should therefore be paired with a sign-change detector of k (or of the 1/mu coefficient) between consecutive family members.

### 8.3 `search/cr3bp_3d_family_tracer._classify_floquet`

READ in the code: with six Floquet multipliers, identifies the trivial pair as the two closest to +1 (rejects if the second closest is more than 0.1 away), then labels "unstable" if any non-trivial multiplier has |lambda| > 1 + unit_tol (and "hyperbolic_pair" if one of those is real), "stable" if all four non-trivial are within unit_tol of the unit circle, else "unstable". Observations:
- Unlike `floquet_classify`, "stable" is reachable. It tests only the modulus, so it cannot separate stable from critical (lambda = +1 or -1 on the unit circle); the final "unstable" fallthrough covers a non-trivial modulus away from 1 by more than unit_tol.
- It treats the trivial pair as "the two closest to +1". For a 3D orbit with a vertical critical orbit (this paper section 5.2: the vertical pair crossing lambda = +1 at the plane-3D family intersections, integer tau/pi, eta/pi) a non-trivial pair also sits at +1 and the identification "two closest to +1" becomes ambiguous (four eigenvalues within 1e-1 of 1). Result: either "degenerate" or a mis-identified pair; no error is raised. This is a case worth a synthetic test at A_i(0) and at the plane-3D intersections once a 3D second-species seed exists. INFERRED; not run.
- It returns a coarse tag and not the numerical index, and the project's #682 census stores the tag. For second-species-type orbits, store log10 abs(k) or lambda_max in addition, because the stable windows are of width order mu and carry no signal in a modulus tolerance.
- For the same unit_tol = 1e-3 (modulus): a real pair lambda = 1 + eps has k - 1 = eps^2/2 (COMPUTED), so the tolerance accepts k up to 1 + 5e-7 as "stable"; adequate, and not sensitive to the narrow-window issue.

### 8.4 `monodromy_eigenstructure` in `search/er3bp_periodic`

Not touched by this paper beyond `#931` (accepts any complex eigenvalue within 0.5 of the unit circle as the centre, would read part of a complex quartet as a centre). The stability of a second-species orbit at mu > 0 has no complex quartet in the paper's analysis (one reciprocal pair; the paper's discussion is for real k), so the quartet issue is unrelated to the present material.

## 9. Techniques applicable to the project's problems

### 9.1 `#899` (continuation through near-collision seeds, Earth-Moon mass)

1. Predict the sign and size of k for a second-species seed before continuing. From a seed's (tau, eta) (or from the generating orbit of its (p, q, C) family, via the recipe in the companion digest section 5(a)), evaluate S (eq. 40) and e_H = 1/sin(delta/2) (eq. 9c) and form k_pred = -V S/(mu e_H^2). COMPUTED example orders: for A0(-1) at mu = 1e-6 the leading magnitude of |k| away from S = 0 is of order V |S| / (mu e_H^2) = 1.8437 |S| / (1e-6 * 3.456) (e_H = 1.859) = 5.3e5 |S|; with |S| of order 1 that is of the order of the 8.4e5 multiplier question in `#931`(c). A Floquet multiplier of order 1e5 to 1e6 at a 1e-6 seed is therefore expected and is not evidence of a faulty monodromy; a multiplier of order 1 at mu = 1e-6 would be suspect for a single-encounter second-species orbit.
2. Where the stable windows are: a continuation in mu starting from a 1e-6 seed with k = O(1/mu) will not be stable except in windows around zeros of S along the family (the critical generating orbits at the C extrema, Table I of paper I). The set of C extrema per (i, j) family is finite (two for C_ij: C_ij(1) maximum, C_ij(2) minimum), so for a Casoliva-type fixed-period continuation the C interval in which a stable window can occur at small mu is the neighbourhood of those two values: within about the leading-order width of section 4 at the target mu.
3. Fixed-period continuation hits lunar impact "in most cases" (Casoliva et al. 2010, quoted in `#899`). The paper supplies a diagnostic: the width of the stable window shrinks as mu e_H^2/(V |dS/dl|) and the critical orbits k = +1 and k = -1 pair up and coalesce at mu = 0, so following one of the pair in mu while holding k = +1 (the Henon-Guyot iteration mentioned on p.1039) is a continuation in mu along a codimension-one object whose existence ranges are quoted (m1, m2: mu below mu0 = 0.327, Fig. 1 caption). A continuation that holds a critical (k = +1) orbit rather than a fixed period is not published here in detail but is the paper's stated method (plan (1), p.1039).
4. Seed bookkeeping: report for each seed (tau/pi, eta/pi), the integers (i, j) of the nearest integer point, whether it lies on a C extremum (S near 0, check by eq. 40 or paper I eq. 27) and the sign of k_pred, so a continuation that crosses a k sign change is flagged as a window event rather than treated as a bifurcation of unknown type.
5. Validity: eq. 1 assumes O(1) relative speed V and b of order mu; at the Earth-Moon mass, section 4 shows the one-passage estimate loses meaning for V below about 0.6 (the near-resonant direct orbits); do not use k_pred as a ranking metric there (consistent with companion section 5(b)).

### 9.2 `#906` (demanded-turn gate)

- The gate's turn comes from the same quantity as this paper's e_H: sin(delta/2) = 1/e_H (eq. 9c), so e_H = V/|V1| in the companion's notation. A small turn (large e_H) means a distant passage with small 1/mu amplification (k of order V |S|/(mu e_H^2)), so a small-turn near-resonant junction is physically a weakly amplifying one, not a non-encounter. This supports the existing `#906` amendment (return "indeterminate" instead of rejecting at small turn).
- A cheap additional field for the gate's report: e_H = 1/sin(delta/2) and the one-passage amplification estimate V/(mu e_H^2) at the actual mu, so a reader can see whether the encounter is amplifying (O(1/mu)) or nearly passive.

### 9.3 `#931` (stability classification and the first general-three-body controls)

- (b) The classifier: use the reciprocal-quartic route (Hadjidemetriou) with a reported index k; add a synthetic generator k(l) = -c S(l)/mu with a sign change to produce (unstable, stable window, unstable of opposite sign) cases and require the classifier to find the two critical orbits (k = +1 and k = -1) on either side of the zero and to return the sign of k.
- (c) For `#890`/`#895` multipliers computed by integration across a close approach: the product of two O(1/mu) factors per period is not expected unless there are two passages; a multiplier of order mu^-n for n passages is the paper's scaling, so a published or computed value should be compared with mu^-n before it is trusted.
- Do not rely on this paper for any mu > 0 control orbit: it has none (check: tables 1 only, mu = 0).

### 9.4 Positive controls this paper offers (all at mu = 0)

1. Table 1 (section 5.1) as expected values for (a, e, x0, x1, C, T) from (tau/pi, eta/pi, sigma0, sigma1, sigma2), tolerance 1e-3 from the 5-decimal inputs (tighter, 1e-5, if the 10-decimal values of paper I are used as inputs and the printed values compared at their own precision), and A0(-1) agrees to the printed digits with the 10-decimal input. All seven rows reproduced (COMPUTED).
2. Eq. 40 / eq. 42: S(tau, eta, sigma) = 0 at every critical orbit of paper I's Table I (to 1e-8 using 10-decimal inputs for rows tested: A0(1), A0(2), A1(-2), A0(-1)), and the identity -(1/2) sigma rho s_eta^2 S = G* at any (tau, eta, sigma). Two independent functions (eq. 27 of paper I and eq. 40 here) that agree exactly at generic points but are transcribed separately make the identity a check on both transcriptions.
3. Demanded turn of A0(-1) = 65.094 degrees, e_H = 1.859, V = 1.8437 (COMPUTED from the Table 1 row via eq. 9c and the rotating-frame recipe; matches the companion digest's value; the source of the expected value is the paper's own geometry, no project code).
4. Fig. 1 caption, mu0 = 0.327... (existence boundary of m2): an expected value for a future continuation of the m family at mu > 0.
5. Plane-3D intersections (section 5.2): the 3D family emanating from A_i(0) has period 2 pi i; a test of the 3D tracer can seed the retrograde circle of radius 1 at A_i(0) and check the period.

## 10. Errata and ambiguities

- Abstract and text are consistent with the body ("necessary stability condition S = 0 ... identical to a previously obtained sufficient condition G* = 0"): READ; no printed slip found in the equations used above (eq. 40 and eq. 42 verified by computation).
- The caption of Fig. 1 prints the stability statement with wording that is hard to parse in the rotated print (the 0 < mu < 1 range for m3 and the sets mu0 < mu < 1 - mu0); only mu0 = 0.327... is certain.
- The mass ratio mu* = 1/82.30 = 0.01215067... (p.1039): the exact quotient is 0.012150668; the printed seven digits 0.01215067 agree to rounding (COMPUTED).
- Eq. 1 prints several nested braces; the t_p subscripts and the signs inside the braces are not needed below and are not asserted.
- Eq. 10 prints exponents on (b/mu) and (b/mu^2); only the combined structure is used.
- Footnote to eq. 7 refers to "second comment at the end of this section", which is the Jacobi-constant remark on p.1032 (READ).
- The text layer of the PDF was not used. If a text layer exists it garbles the equations; all content above is from the images.

## 11. Follow-ups (no task numbers registered here)

1. Implement eq. 40 and eq. 42 (S and G*) in a `second_species_arcs` module with the two identity tests of section 9.4, and Table 1 as sourced expected values (7 rows, tolerance as stated).
2. Add k_pred = -V S/(mu e_H^2) as a seed diagnostic for `#899`, with a test against a published k at a small mass. Candidate comparison set: Casoliva et al. 2010 Table 3 (k for the nine catalogued rows at the Earth-Moon mass) and the 2008 conference seeds at 1e-6; the sign convention of k_pred must be tested before the sign is used (INFERRED in section 4).
3. Fix `floquet_classify` (label logic and docstring) with k, sign(k), and a k = +-1 test; add a family-level sign-change detector so narrow stable windows are not stepped over. Add a synthetic test for the 3D tracer's trivial-pair identification when a vertical critical orbit puts a second eigenvalue at +1 (A_i(0) intersections).
4. Locate the two forthcoming papers (numerical k = +-1 iteration at mu > 0, including the Earth-Moon mu* = 1/82.30 narrow-window results, and the second-order analytic stability boundaries); these would supply the first mu > 0 numbers for stability of second-species orbits. Candidate search: Hitzl, "Stability of second species orbits ... (mu small)" 1977-1980 and Hitzl and Henon 1978 follow-ups; also Henon 1965 (Ann. Astrophys. 28:992) and Henon-Guyot 1970 (the critical-family classification, k = AD + BC), both cited and not held.
5. Test the claim of section 4 that the stable window shrinks like mu: integrate a family near A0(-1) at mu = 1e-3, 1e-4, 1e-5 with a regularised CR3BP integrator and measure the width of the stable window against mu e_H^2/(V |dS/dl|); then the Earth-Moon mu = 0.01215 comparison with the 0.0053 rad estimate. Positive control: Fig. 1 (m1 and m2 for mu = 0.1 to 0.3).
6. Check what happens to the 1/mu amplification near a bifurcation (small turn, e_H large): the formula of section 4 with 1/e_H^2 = sin^2(delta/2) suggests a reduction, which would make near-resonant seeds less unstable than the generic mu^-1; test with Henon 2001 part A's near-bifurcation scalings (Delta C of order mu^nu, nu in (0, 1/2)) and a computed orbit.
7. Extend the one-passage recurrence to n passages per period (product of n matrices M_i) for chains with several lunar encounters (`#899` for more than one flyby; two moons): derive the n-fold scaling and the conditions for cancellation; not treated in this paper.
8. The 3D statement (section 5.2): A_i(0) are the only second species orbits critical in both directions; verify numerically that the 3D family from A_i(0) of Fig. 12 closes with period 2 pi i, as a control for the 3D tracer's treatment of the retrograde circle.
9. Note the small difference between the paper's mu* = 1/82.30 (0.012150668) and the project's 0.0121529529 (relative 1.9e-4) wherever a published Earth-Moon result is compared with project output; the effect on window widths of order mu is of the same relative size, negligible, but a test pinned to a printed value should use the printed mu.
10. Add the `docs/notes/CORPUS_INDEX.md` entry (done with this digest) and a dated note in the companion digest (done with this digest).
