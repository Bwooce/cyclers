# Digest: Barrabes & Gomez 2003, "Three-dimensional p-q resonant orbits close to second species solutions" (#920)

Celestial Mechanics and Dynamical Astronomy 85:145-174 (2003), received 14 Jan 2002, revised 14 May 2002, accepted 17 Jun 2002.
Filed in the private paper corpus as `barrabes-gomez-2003-three-dimensional-pq-resonant-orbits-second-species-solutions-cmda-85-2-doi-10.1023-A1022098510161.pdf`
(30 PDF pages; PDF page n is journal page 144+n; text layer present, `pdftotext | wc -w` = 11102).
Status: READ (all 30 page images) unless marked COMPUTED (my own arithmetic, 2026-10-04) or INFERRED. Page references are journal pages.
Companion: `docs/notes/2026-10-04-digest-barrabes-gomez-2002-spatial-pq-resonant-orbits.md` (necessary conditions and the outer map). An earlier short treatment:
`docs/notes/2026-07-28-744-broucke-leiva-barrabes-earth-moon-lineage-digest.md` section 2 (not repeated). The seed recipe built from both papers:
`docs/notes/2026-10-04-digest-casoliva-2008-aiaa-families-cycler-trajectories-seeds.md`.

## 1. What the paper gives
READ (abstract): "the three-dimensional p-q resonant orbits that are close to periodic second species solutions (SSS) ... The work is based on an analytic study of the in- and out-maps.
These maps are associated to follow, under the flow of the problem, initial conditions on a sphere of radius mu^alpha around the small primary, and consider the images of those initial points on the same
sphere. The out-map is associated to follow the flow forward in time and the in-map backwards. For both mappings we give analytical expressions in powers of the mass parameter. Once these expressions are
obtained, we proceed to the study of the matching equations between both, obtaining initial conditions of orbits that will be 'periodic' with an error of the order mu^(1-alpha), for some alpha in (1/3, 1/2).
Since, as mu -> 0, the inner solution and the outer solution will collide with the small primary, these orbits will be close to SSS."
Deliverable for a programmer: for given (p, q, C_J) a pair (psi_0, theta) (planar), or (phi_0 = +-pi/2, theta, ...) (spatial), plus eps = delta, at which the out-map and in-map agree through
order mu^alpha. This is "periodic" in the sense that the orbit leaves the sphere, returns to it at the matching point, and the in-map composed with the out-map is the identity to that order; it is NOT
a corrected periodic orbit (the corrector is the differential correction in Casoliva et al.).

## 2. Model and conventions (pp.145-146, 148)
Same as the 2002 paper: spatial CRTBP, masses 1-mu (E, big) and mu (M, small), synodic Oxyz with E at (mu,0,0), M at (mu-1,0,0); Eq. 1: r'' + A r' = grad Omega(r), A = [[0,-2,0],[2,0,0],[0,0,0]];
Omega = (x^2+y^2)/2 + (1-mu)/r_1 + mu/r_2; Eq. 3: |r'|^2 = 2 Omega - C_J. Definition (Eq. 4, p.146), verbatim: "these orbits leave at t = t_1 a sphere B in the configuration space, with center at M and
radius mu^alpha and return for the first time to the same sphere at an epoch t_2 such that t_2 - t_1 = 2 pi q + eps mu^alpha + O(mu^(2 alpha)) = 2 pi tau p + delta mu^alpha + O(mu^(2 alpha)), where p, q in N are
relatively prime, eps and delta are suitable constants and 2 pi tau is the period of the elliptic orbit which approaches the motion of P far from M (Font, 2002)." Initial time t_1 = 0 (p.148).
Labelling: p = revolutions of P about E, q = revolutions of M (period 2 pi); tau = q/p; a = (q/p)^(2/3) (see 2002 digest section 2).
The phrase in the abstract of the 2002 paper about "p around M" is not repeated here; the introduction (p.146) says "Approximately, perform p revolutions around E while M does q revolutions around the center of masses".
Sphere B = B(M, mu^alpha), alpha in (1/3, 1/2) (p.146, abstract).

## 3. Summary of known results (section 2, pp.147-149) = the 2002 paper
Spherical coordinates (Eq. 5, p.148):
  r_i = ( mu - 1 + mu^alpha cos(varphi) cos(theta), mu^alpha cos(varphi) sin(theta), mu^alpha sin(varphi) )^T,  r'_i = v_i ( cos(phi) cos(psi), cos(phi) sin(psi), sin(phi) )^T.
The text under Eq. 5 prints (verbatim): "cos a = cos varphi cos phi (theta - psi) + sin varphi sin phi, where a is the angle between r_2i and r_M." [sic]. This is garbled: the angle a is between r_2i (position from M) and
r'_i (velocity), and the correct expression, printed in Eq. 6 of this paper and in Eq. 18 of the 2002 paper, is cos a_0 = cos(varphi) cos(phi_0) cos(theta - psi_0) + sin(varphi) sin(phi_0). INFERRED: typesetting error; follow Eq. 6.
Eq. 6: phi = phi_0 + Dphi mu^alpha + O(mu^(2 alpha)), psi = psi_0 + Dpsi mu^alpha + ..., cos a = cos a_0 + DC mu^alpha + ..., with
cos a_0 = cos(varphi) cos(phi_0) cos(theta - psi_0) + sin(varphi) sin(phi_0), DC = Lambda_0 Dpsi + Dphi ( sin(varphi) cos(phi_0) - cos(varphi) sin(phi_0) cos(psi_0 - theta) ), Lambda_0 = cos(varphi) cos(phi_0) sin(theta - psi_0).
The restriction set (called "Eq. 7, 8, 9", p.148), a summary of the 2002 results:
  (eps - delta)^2 + delta^2 (3 - C_J) + 2 (eps - delta) cos(varphi) sin(theta) + 2 delta sqrt(3 - C_J) cos(a_0) + 2 (eps - delta) delta sqrt(3 - C_J) cos(phi_0) sin(psi_0) = 0,             (7)
  C_J - 2 + 2 sqrt(3 - C_J) cos(phi_0) sin(psi_0) = (p/q)^(2/3),                                                                                          (8)
  6 pi q (q/p)^(2/3) ( sqrt(3 - C_J) ( cos(varphi) cos(phi_0) sin(psi_0 - theta) + Dphi sin(phi_0) sin(psi_0) - Dpsi cos(phi_0) cos(psi_0) ) - 2 cos(varphi) cos(theta) ) = eps - delta.     (9)
Interval (Eq. 10, p.149): C_J1 = (p/q)^(2/3) - 2 sqrt(2 - (p/q)^(2/3)), C_J2 = (p/q)^(2/3) + 2 sqrt(2 - (p/q)^(2/3)); both endpoints excluded ("if C_J was equal C_J1 or C_J2, then cos^2 phi_0 sin^2 psi_0 = 1 which is
not acceptable").
Eq. 11 (p.149), printed verbatim: "psi_0 not equal to 0 or C_J not equal to (p/q)^(2/3)". DISCREPANCY: the 2002 paper's Eq. 21 and Eq. 28 print "phi_0 not equal to 0 or ..." (the vertical angle phi_0, not the azimuth psi_0).
INFERRED: Eq. 11 here is a typo for phi_0, because the angular momentum c = R_i x R'_i vanishes only when phi_0 = 0 AND sin psi_0 = 1/sqrt(3 - C_J) (2002 Eq. 21), and the latter equals C_J = (p/q)^(2/3) by Eq. 8 when phi_0 = 0.
Eq. 12: v_i^2 = 3 - C_J + 2 mu^(1-alpha) + O(mu^(2 alpha)).

## 4. Out-map (section 3, pp.149-152)
Proposition 3.1 (p.150): if there are t_(alpha/2) < t'_(alpha/2) in (t_1, t_2) with R_2(t) >= mu^(alpha/2) between them and cos a >= eps > 0, then R(t) = R_TB(t) + O(mu^(1-alpha)) and R'(t) = R'_TB(t) + O(mu^(1-alpha))
for all t in [t_1, t_2] (Eq. 13, the content of 2002 Theorem 2.3).
Out-map in sidereal coordinates (Eq. 14): R(t_2) = R_i + R'_i delta mu^alpha + O(mu^(1-alpha)), R'(t_2) = R'_i - R_i delta mu^alpha + O(mu^(1-alpha)).
In synodic coordinates (Eq. 15a, b): r_e = ( I - ((eps - delta)/2) mu^alpha A ) r_i + delta mu^alpha r'_i;  r'_e = ( I - ((eps + delta)/2) mu^alpha A ) r'_i.
Exit spherical coordinates (p.151): r_e = (mu - 1, 0, 0)^T + r_2e (cos(varphi_e) cos(theta_e), cos(varphi_e) sin(theta_e), sin(varphi_e))^T; r'_e = v_e (cos(phi_e) cos(psi_e), cos(phi_e) sin(psi_e), sin(phi_e))^T.
r_2e = mu^alpha ( 1 + mu^alpha Dr_e + O(mu^(1-alpha)) ), Dr_e as in 2002 Eq. 33 (printed with the Dphi bracket
"sin(varphi) cos(phi_0) - sin(phi_0)( cos(varphi) cos(psi_0 - theta) + (eps - delta) sin(psi_0) )", algebraically equal to the 2002 form).
Eq. 16a-c (p.152) and Eq. 17a, b are the 2002 Eq. 34a-c and 35a, b; used here at zero order and first order. The printed 16a, b are the same as 2002 34a, b (compared term by term on the page images).
Notation conflict to watch: here "v_e = v_i + O(mu^(2 alpha))" and "phi_e = phi + O(mu^(2 alpha))" (Eq. 17a), "psi_e = psi_0 + (Dpsi - (eps + delta)) mu^alpha + O(mu^(2 alpha))" (Eq. 17b).

## 5. Inner solution (section 4, pp.152-161)
Idea: inside B the Earth acts as a perturbation and the orbit is close to a hyperbola about M; regularise by Kustaanheimo-Stiefel (KS), blow up scale mu^alpha to unit radius, linearise, then correct the nonlinear term of order mu^alpha with one extra step. 
### 5.1 Variables (Eq. 18-24)
P = r - r_M, Q = r' + (1/2) A r_2 (Eq. 18); Jacobi in (P,Q) (Eq. 20); sidereal-aligned tilde variables P~ = G(t) P, Q~ = G(t) Q (Eq. 21, 22); Jacobi (Eq. 23). KS map with scale (Eq. 24):
P~ = mu^alpha L(u) u, dt/ds = |P~|, with L(u) the 4x4 KS matrix [[u1,-u2,-u3,u4],[u2,u1,-u4,-u3],[u3,u4,u1,u2],[u4,-u3,u2,-u1]], so |P~| = mu^alpha |u|^2, |P~| <= mu^alpha iff |u| <= 1, and t = -mu^alpha integral_s^0 |u|^2 d tau
(initial time t = 0 at s = 0 on the sphere |u| = 1). Bilinear relation l(u,u') = u4' u1 - u3' u2 + u3 u2' - u4 u1' = 0 (Eq. 25), Lemma 4.1 (Stiefel-Scheifele).
Equations of motion (Eq. 27 to 29): with c = (3 - C_J)/4,
  u'' = c u + mu^alpha f(u, u') u + O(mu^(2 alpha)),   f(u, u') = -u2 u1' + u1 u2' - u4 u3' + u3 u4'.                                  (29)
Jacobi (Eq. 28, 30): |u'|^2/|u|^2 = c + mu^(1-alpha)/(2|u|^2) + mu^alpha f(u,u') + O(mu^(2 alpha)). Initial data: u_i with P~_i = mu^alpha L(u_i) u_i, u'_i = (1/2) L(u_i)^T Q~_i.
### 5.2 Approximation (Eq. 31-32, section 4.2)
THEOREM 4.3 (p.157): for |u(s)| <= 1, |zeta(s) - zeta_0(s)| <= |s| K mu^alpha e^(c~ |s|), c~ = max(1, c), where zeta_0 solves the LINEAR problem u'' = c u (K constant from the Gronwall Lemma 4.2). Error in synodic coordinates O(mu^alpha) in r_2 and r'.
This is not accurate enough to match the outer solution (error O(mu^(1-alpha))), so one iterates once: LEMMA 4.4 (p.158): with u_0 the solution of u'' = c u, u(0) = u_i, u'(0) = u'_i, f(u_0(s), u_0'(s)) = f(u_i, u'_i) for all s (f is constant along the linear flow;
proof by the explicit cosh/sinh solution). Hence u'' = c_alpha u + O(mu^(2 alpha)) with
  c_alpha = (3 - C_J)/4 + mu^alpha f(u_i, u'_i),
and THEOREM 4.5 (p.159): |zeta(s) - zeta_h(s)| <= |s| K mu^(2 alpha) e^(c~_alpha |s|) (zeta_h solves u'' = c_alpha u with the same initial data; the time inside B^* is assumed bounded).
The hyperbolic approximation (Eq. 32): u_h(s) = cosh(sqrt(c_alpha) s) u_i + (sinh(sqrt(c_alpha) s)/sqrt(c_alpha)) u'_i, u_h'(s) = sqrt(c_alpha) sinh(sqrt(c_alpha) s) u_i + cosh(sqrt(c_alpha) s) u'_i.
Exit time (Eq. 33): with V_i = |u'_i|, m_alpha = 1 + V_i^2/c_alpha, n_alpha = 2 V_i cos(eta_i)/sqrt(c_alpha) (eta_i = angle between u_i and u'_i), the (backward) time to reach |u_h| = 1 is
  s_f = (1/(2 sqrt(c_alpha))) ln( (m_alpha - n_alpha)/(m_alpha + n_alpha) ),   s_f < 0 needs n_alpha > 0 ,
and <u_i, u'_i> = V_i cos eta_i = (1/2) v_i cos a (p.160), so n_alpha > 0 iff cos a > 0. For the expansion one needs 1 - cos a >= eps' > 0 (the orbit must not leave the sphere in the direction of the radius vector, "In this case, when we move back in time we go towards the centre").
Eq. 34: s_f = (1/(2 sqrt(c))) ln( (1 - cos a)/(1 + cos a) ) + O(mu^alpha). Minimum distance (p.161): |u_h(s_m)|^2 = 1 - (1/2)( m_alpha - sqrt(m_alpha^2 - n_alpha^2) ).
So the admissibility window for the leaving angle is 0 < cos a < 1 (strictly), and the inner flyby geometry is fixed by (C_J, a).
### 5.3 In-map, exact formulas (section 5, pp.161-164)
Definition (p.161): (r_f, r'_f) = (r_h(t_f), r'_h(t_f)), t_f = -mu^alpha integral_(s_f)^0 |u(s)|^2 ds. At s_f (Eq. 35), with Delta_alpha = sqrt(m_alpha^2 - n_alpha^2) = sqrt( (1 + V_i^2/c_alpha)^2 - 4 V_i^2 cos^2(eta_i)/c_alpha ):
  u_f = (1/Delta_alpha) ( (1 + V_i^2/c_alpha) u_i - 2 (V_i cos eta_i / c_alpha) u'_i ),   u'_f = (1/Delta_alpha) ( -2 V_i cos(eta_i) u_i + (1 + V_i^2/c_alpha) u'_i ).
Then P~_f = (1/Delta_alpha^2)(T1 P_i + (mu^alpha/2) T2 Q_i), Q~_f = (1/Delta_alpha^2)(2 mu^(-alpha) T3 P_i + T1 Q_i) with
  T1 = (1 + V_i^2/c_alpha)^2 - 4 V_i^4 cos^2(eta_i)/c_alpha^2,
  T2 = 4 (V_i cos eta_i/c_alpha) ( 2 V_i^2 cos^2(eta_i)/c_alpha - 1 - V_i^2/c_alpha ),
  T3 = -2 V_i cos(eta_i) ( 1 - V_i^4/c_alpha^2 ).
Undo the rotation with G(t_f), t_f = -mu^alpha cos a / sqrt(c) + O(mu^(2 alpha)) (p.162); then (Eq. 36)
  P_f = ( I + (cos a/(2 sqrt c)) mu^alpha A ) P~_f + O(mu^(3 alpha)),   Q_f = ( I + (cos a/(2 sqrt c)) mu^alpha A ) Q~_f + O(mu^(2 alpha)).
Expanded (p.163):
  P_f = ( I + (cos a/(2 sqrt c)) mu^alpha A - (cot^2 a/(2c)) mu^(1-alpha) I ) P_i + (cos a/sqrt c) mu^alpha ( -I + ( (f(u_i,u'_i)/c) I - (cos a/(2 sqrt c)) A ) mu^alpha + (cot^2 a/(4c)) mu^(1-alpha) I ) Q_i + O(mu^(3 alpha)),
  Q_f = ( (cos a/(sin^2 a sqrt c)) mu^(1-alpha) ) ( mu^(-alpha) P_i ) + ( I + (cos a/(2 sqrt c)) mu^alpha A - (cot^2 a/(2c)) mu^(1-alpha) I ) Q_i + O(mu^(2 alpha)).
[The printed Q_f first bracket reads "(cos a/(sin^2 a) mu^(1-alpha)/sqrt c)(mu^(-alpha) P_i)"; I keep the printed grouping.] Back to (r, r') (Eq. 37, using |r_2i| = mu^alpha):
  r_2f = ( 1 - (cot^2 a/(2c)) mu^(1-alpha) ) r_2i + (cos a/sqrt c) mu^alpha ( -1 + (f(u_i,u'_i)/c) mu^alpha + (cot^2 a/(4c)) mu^(1-alpha) ) r'_i - (cos^2 a/(2c)) mu^(2 alpha) A r'_i + O(mu^(3 alpha)),
  r'_f = ( (cos a/(sin^2 a sqrt c)) mu^(1-alpha) ) ( mu^(-alpha) r_2i ) + ( I + (cos a/sqrt c) mu^alpha A - (cot^2 a/(2c)) mu^(1-alpha) I ) r'_i + O(mu^(2 alpha)),
 with r_2f = |r_f - r_M| = mu^alpha (1 + O(mu^(2 alpha))) and v_f = 2 sqrt(c) ( 1 + mu^(1-alpha)/(4c) + O(mu^(2 alpha)) ).
Note: v_f = 2 sqrt(c) (1 + ...) = sqrt(3 - C_J)(1 + ...), consistent with v_i^2 = 3 - C_J + 2 mu^(1-alpha) (Eq. 12). COMPUTED: 4c = 3 - C_J; the 1 + mu^(1-alpha)/(4c) term in v_f gives v_f^2 = (3 - C_J) + 2 mu^(1-alpha), same as v_i^2. Agrees.
Direction cosines of the in-map exit point (Eq. 38, p.164) and the velocity angles (Eq. 39). With cos a_0, Lambda_0, DC as in Eq. 6:
  cos(varphi_f) cos(theta_f) = cos(varphi) cos(theta) - 2 cos a_0 cos(phi_0) cos(psi_0) + 2 mu^alpha [ cos a_0 ( Dpsi cos phi_0 sin psi_0 + Dphi sin phi_0 cos psi_0 ) - DC cos phi_0 cos psi_0 + (cos a_0/sqrt c) cos phi_0 ( cos a_0 sin psi_0 - Lambda_0 cos psi_0 ) ] + O(mu^(1-alpha)),
  cos(varphi_f) sin(theta_f) = cos(varphi) sin(theta) - 2 cos a_0 cos(phi_0) sin(psi_0) - 2 mu^alpha [ cos a_0 ( Dpsi cos phi_0 cos psi_0 - Dphi sin phi_0 sin psi_0 ) + DC cos phi_0 sin psi_0 + (cos a_0/sqrt c) cos phi_0 ( cos a_0 cos psi_0 + Lambda_0 sin psi_0 ) ] + O(mu^(1-alpha)),
  sin(varphi_f) = sin(varphi) - 2 cos a_0 sin(phi_0) - 2 mu^alpha [ Dphi cos a_0 cos phi_0 + DC sin phi_0 + (cos a_0/sqrt c) Lambda_0 sin phi_0 ] + O(mu^(1-alpha)),                                      (38)
  v_f = v_i + O(mu^(2 alpha)),   phi_f = phi + O(mu^(2 alpha)),   psi_f = psi + ( Dpsi + 2 cos(a_0)/sqrt(c) ) mu^alpha + O(mu^(2 alpha)).                                     (39)
Zero-order reading of (38) (INFERRED, geometric): position direction after the inner flyby is the entry direction reflected, r_hat_f = r_hat - 2 (cos a_0) v_hat, a hyperbola of vanishing impact parameter passing through the point mass (a collision orbit). The exit velocity direction is unchanged at zero order (39).

## 6. Matching (section 6, pp.164-168): the equations to solve
Conditions: out-map (Eq. 16, 17a) equals in-map (Eq. 38, 39) through order mu^alpha, plus p-q resonance (Eq. 7) and 0 < cos a_0 < 1.
Zero order (p.165), writing C := delta sqrt(3 - C_J):
  cos varphi cos theta + C cos phi_0 cos psi_0 = cos varphi cos theta - 2 cos a_0 cos phi_0 cos psi_0,
  cos varphi sin theta + C cos phi_0 sin psi_0 + eps - delta = cos varphi sin theta - 2 cos a_0 cos phi_0 sin psi_0,
  sin varphi + C sin phi_0 = sin varphi - 2 cos a_0 sin phi_0.
Result (Eq. 40): "From them we infer that eps - delta = 0 and
  eps = delta = - cos(a_0)/sqrt(c)."                                                                                                             (40)
COMPUTED consistency check: with sqrt(3 - C_J) = 2 sqrt(c), C = delta * 2 sqrt(c) = -2 cos a_0, which is exactly the in-map's -2 cos a_0 shift in the three zero-order equations, and the middle equation then gives eps - delta = 0. Consistent.
So the sphere-to-sphere time is t_2 - t_1 = 2 pi q - mu^alpha cos(a_0)/sqrt(c) + O(mu^(2 alpha)) = 2 pi q - 2 mu^alpha cos(a_0)/sqrt(3 - C_J) (COMPUTED substitution, c = (3 - C_J)/4); eps and delta satisfy the resonance ellipse (7) (paper's remark).
Not to be confused with the PERIOD of the closed orbit: this t_2 - t_1 is the time between leaving and returning to the sphere along the OUTER arc; the inner transit inside B adds time of order mu^alpha cos(a)/sqrt(c) of the opposite sign (INFERRED from t_f = -mu^alpha cos a/sqrt(c) at p.162).
Do not write T = 2 pi q + eps mu^alpha for the closed orbit period. Casoliva et al. Eq. 18 (printed T = 2 pi q + O(mu^alpha)) gives no coefficient, consistent with this.
Velocity condition: the order-mu^alpha terms of (17b) and (39) match iff eps + delta = -2 cos(a_0)/sqrt(c), true by (40): "the expressions of the velocities do not bring any new information".
Order-mu^alpha terms of positions (p.165-166): coefficients O_1, O_2, O_3 (out) and I_1, I_2, I_3 (in), with eps = delta from (40). Printed:
  O_1 = 2 cos a_0 ( Dpsi ( sin psi_0 cos phi_0 + Lambda_0 ( cos varphi cos theta - 2 cos a_0 cos phi_0 cos psi_0 ) ) + Dphi ( sin phi_0 cos psi_0 + ( sin varphi cos phi_0 - cos varphi sin phi_0 cos(psi_0 - theta) ) ( cos varphi cos theta - 2 cos a_0 cos phi_0 cos psi_0 ) ) ),
  O_2 = -2 cos a_0 ( Dpsi ( cos psi_0 cos phi_0 - Lambda_0 ( cos varphi sin theta - 2 cos a_0 cos phi_0 sin psi_0 ) ) - Dphi ( sin phi_0 sin psi_0 + ( sin varphi cos phi_0 - cos varphi sin phi_0 cos(psi_0 - theta) ) ( cos varphi sin theta - 2 cos a_0 cos phi_0 sin psi_0 ) ) ),
  O_3 = -2 cos a_0 ( Dphi ( cos phi_0 - ( sin varphi cos phi_0 - cos varphi sin phi_0 cos(psi_0 - theta) ) ( sin varphi - 2 cos a_0 sin phi_0 ) ) - Dpsi Lambda_0 ( sin varphi - 2 cos a_0 sin phi_0 ) ),
  I_1, I_2, I_3 as printed on p.166 (they are linear in Dpsi, Dphi with the extra c^(-1/2) 2 cos a_0 ... terms).
(The printed O_i and I_i are long and I have not re-derived them; the cell-by-cell grouping above follows the page and the parentheses on p.165 are partly ambiguous in print: UNCLEAR in O_1, O_2, O_3. A programmer should obtain O_i - I_i by differentiating (16)-(17) and (38)-(39) symbolically
rather than typing these.) Setting O_i = I_i gives the linear system (Eq. 42, p.166) in the unknowns Dpsi, Dphi:
  (Lambda_0 G_1) Dpsi + (F G_1) Dphi = C_1,  (Lambda_0 G_2) Dpsi + (F G_2) Dphi = C_2,  (Lambda_0 G_3) Dpsi + (F G_3) Dphi = C_3,                           (42)
with F = sin(varphi) cos(phi_0) - cos(varphi) sin(phi_0) cos(psi_0 - theta) and
  G_1 = cos a_0 ( cos varphi cos theta - 2 cos a_0 cos phi_0 cos psi_0 ) + cos phi_0 cos psi_0,
  G_2 = cos a_0 ( cos varphi sin theta - 2 cos a_0 cos phi_0 sin psi_0 ) + cos phi_0 sin psi_0,
  G_3 = cos a_0 ( sin varphi - 2 cos a_0 sin phi_0 ) + sin phi_0,
  C_1 = c^(-1/2) cos a_0 cos phi_0 ( cos a_0 sin psi_0 - Lambda_0 cos psi_0 ),
  C_2 = -c^(-1/2) cos a_0 cos phi_0 ( cos a_0 cos psi_0 + Lambda_0 sin psi_0 ),
  C_3 = -c^(-1/2) cos a_0 Lambda_0 sin phi_0.                                                                                                      (43)
Together with the resonance condition (9) with eps - delta = 0, which becomes (Eq. 41, p.165):
  - cos phi_0 cos psi_0 Dpsi + sin phi_0 sin psi_0 Dphi = cos(varphi) cos(theta)/sqrt(c) + Lambda_0.                                                  (41)
Compatibility (p.167): the matrix M = [Lambda_0 G_i, F G_i] has rank 0 only if F = Lambda_0 = 0, excluded; otherwise rank 1, and (42) is compatible iff all 2x2 determinants Delta_ij = |G_i C_i; G_j C_j| vanish, with Delta_12, Delta_23, Delta_13 printed
(multiplied by Delta = -cos^2 a_0 cos phi_0/sqrt(c)).
- Planar case (varphi = phi_0 = 0): Delta_ij = 0 for all i, j, so (42) has a solution with one degree of freedom. 
- Spatial case: the only acceptable solution (0 < cos a_0 < 1) of Delta_ij = 0 is phi_0 = +-pi/2. Then: if theta - psi_0 not equal to pi/2 + k pi, (42) has a one-parameter solution; if theta - psi_0 = pi/2 + k pi, (42) reduces to 0 = 0.

### 6.1 Planar case (pp.168-170): the formulas that give the seed
With varphi = phi_0 = 0 (Dphi disappears), and cos a_0 = cos(theta - psi_0), (42) reduces to one equation and, with (41):
  sin^2(theta - psi_0) Dpsi = -cos(theta - psi_0)/sqrt(c),   -cos(psi_0) Dpsi = cos(theta)/sqrt(c) + sin(theta - psi_0).                              (44)
The compatibility (eliminating Dpsi) is
  E(theta, psi_0) = (1/sqrt c)( cos(theta) sin^2(theta - psi_0) - cos(psi_0) cos(theta - psi_0) ) + sin^3(theta - psi_0) = 0.                         (45)
Note: the 2008 and 2010 Casoliva papers print (17) as (2/sqrt(3 - C_J)) ( cos theta sin^2(theta - psi) - cos psi cos(theta - psi) ) + sin^3(theta - psi) = 0, which equals (45) since 1/sqrt(c) = 2/sqrt(3 - C_J). COMPUTED identity.
Properties printed (p.168): E(theta + pi, psi_0) = -E(theta, psi_0); at C_J = -1 (c = 1) psi_0 = pi/2 solves (45); (theta, 0) with tan^3 theta = 1/sqrt(c) and (0, psi_0) with sin^3 psi_0 - (2/sqrt c) sin^2 psi_0 + 1/sqrt(c) = 0 are solutions.
The additional conditions are the 2002 Eq. 8 and 11, planar:
  sin(psi_0) = ( 2 - C_J + (p/q)^(2/3) ) / ( 2 sqrt(3 - C_J) ),   C_J not equal to (p/q)^(2/3),   for C_J in (C_J1, C_J2).                              (46)
Case analysis (p.169), verbatim in substance:
- p/q = 1: sin psi_0 = sqrt(3 - C_J)/2 = sqrt(c). Then psi_0 in [0, pi] \ {pi/2}; "there is only one value for C_J for every psi_0".
- p/q < 1: for psi_0 in [0, 2 pi] \ {pi/2, 3 pi/2} one value of C_J, from sqrt(3 - C_J) = 2 sqrt(c) = sin psi_0 + sqrt( sin^2 psi_0 + 1 - (p/q)^(2/3) ).
- p/q > 1: need sin psi_0 >= sqrt( (p/q)^(2/3) - 1 ) and two values: 2 sqrt(c) = sin psi_0 +- sqrt( sin^2 psi_0 + 1 - (p/q)^(2/3) ).
Existence of theta (p.169-170): E(theta, psi_0) = (1/4) r ( 3 cos(theta - sigma_1) + cos(3 theta - sigma_2) ), r = (3 - 2c - (p/q)^(2/3))/(2c), sigma_1, sigma_2 depending on psi_0 and c, so for fixed psi_0 E = 0 has exactly the same number of solutions
as 3 cos x + cos(3x + s) = 0, namely two (theta_1, theta_2 = theta_1 + pi). Concluding: "If p/q = 1 ... the solutions of (45) are theta = psi_0 +- pi/2. But for these values cos a_0 = cos(theta - psi_0) = 0 so they are not admissible. If p/q not equal to 1 ... we get two solutions theta_1 and theta_2 such that
theta_2 - theta_1 = pi. This implies cos(theta_1 - psi_0) cos(theta_2 - psi_0) < 0. Therefore we only get one admissible value for theta (which verifies cos a_0 = cos(theta - psi_0) > 0)."
So for p not equal to q: for each admissible (C_J, psi_0) there is exactly one theta with cos(theta - psi_0) > 0. Hence two seeds per C_J (two psi_0 branches in Eq. 46).
Why p = q = 1 fails (p.170), verbatim: "However, we cannot say the same for p = q = 1. This can be explained as follows: the generating orbits (for mu = 0) associated to p-q resonant orbits are bifurcation orbits of 1st species-2nd species (Henon, 1997). When p = q = 1, a generating orbit of 1st species can only be an ellipse intersecting the
unit circle twice (type 1) or a retrograde circle of radius 1 and period 2 pi (type 3). Henon (1997, p. 102) proves that there cannot exist orbits of type 1 of 1st species-2nd species. If the generating orbit was of type 3 then the third body P and the secondary M would have a previous encounter at a half period, which cannot be possible because we have supposed that
|r(t) - r_M(t)| > mu^alpha for t in (t_1, t_2)."
Figure 3 (p.171): planar p-q resonant orbits close to periodic SSS: 1-2 resonant at C_J = -0.850431 and 1.059752; 2-1 resonant at C_J = 0.303724 and 0.565699; 2-3 and 3-2 resonant orbits (values of C_J not printed in the caption).
COMPUTED: all four printed C_J lie inside the admissible interval (1-2: (-1.711013, 2.970934); 2-1: (0.302724, 2.872078)). The 2-1 value 0.303724 is exactly 1.000e-3 above the lower end C_J1 = 0.302724 (COMPUTED difference 0.001000);
INFERRED: the authors used C_J1 + 0.001 as an extreme-near-endpoint example. Casoliva's 21a (C_J = 0.3044238 at mu = 1e-6, Table 2 of the AIAA 2008 paper) lies 1.7e-3 above the lower end, which the planar theory treats as the near-endpoint regime.

## 7. Spatial case (section 6.2, pp.170-173)
As found, the system (42) has solutions only for phi_0 = sigma pi/2, sigma = +-1. Two sub-cases:
- psi_0 - theta = pi/2 + k pi: (42) reduces to 0 = 0, and (41) gives Dphi = sigma (-1)^k cos(varphi)/sqrt(c) if cos theta not equal to 0 (indefinite if cos theta = 0), Dpsi arbitrary.
- psi_0 - theta not equal to pi/2 + k pi: (42) reduces to F G_i Dphi = 0, so Dphi = 0 (the G_i cannot all vanish), and (41) becomes sigma sin(psi_0) Dphi = cos(varphi) cos(theta)/sqrt(c), i.e. cos(varphi) cos(theta) = 0. With cos varphi not equal to 0 (else cos a_0 = +-1), cos theta = 0, theta = pi/2 + m pi.
Restrictions (11) and (8) then give: C_J = 2 + (p/q)^(2/3) (because cos phi_0 = 0 in Eq. 8), which requires p < q (so that C_J < 3).
Figure 4 (p.172): spatial examples with theta = pi/2, Dphi = 0: 1-2, 1-3, 1-5 and 3-5 resonant (left: phi_0 = pi/2, right: phi_0 = -pi/2). Initial conditions up to order mu^alpha (p.173):
  x = mu - 1, y = mu^alpha cos(varphi), z = mu^alpha sin(varphi);  x' = 0, y' = 0, z' = v sin(phi_0);
integrated from t = 0 to t = 2 pi q + eps mu^alpha. (Reading: theta = pi/2 puts the sphere point at x offset 0; the velocity is purely vertical, so the orbit leaves the Moon along z. The spatial family is the one that is a strictly vertical flyby; INFERRED.)
Conclusions (p.173): "In this paper, we have proved the existence of solutions of the spatial RTBP close to SSS. The study is done analyzing the conditions under which the inner and outer solutions, both on a sphere of radius mu^alpha around the small primary, match up to a certain order in terms of mu^alpha."
(The conclusion says "proved the existence"; the body proves matching of the asymptotic maps to the stated order, not existence of exact periodic orbits. READ; the distinction matters for seeds, which must still be corrected.)

## 8. What fails and why (summary)
- p = q = 1: no solution (bifurcation orbit of first-species/second-species type; Henon 1997 p.102).
- cos a_0 = 0 or cos a_0 = +-1: the inner hyperbola degenerates (tangential leaving, or leaving along the radius).
- C_J at the interval endpoints (cos^2 phi_0 sin^2 psi_0 = 1) and C_J = (p/q)^(2/3) when phi_0 = 0 (collision orbit of the Kepler problem about E): excluded.
- Spatial: phi_0 not equal to +-pi/2 has no solution to the first matching order; the whole spatial family has psi_0 and theta tied by cos theta = 0, and C_J = 2 + (p/q)^(2/3) with p < q.
- The result is valid for alpha in (1/3, 1/2) and mu below an unquantified bound (not stated). Errors O(mu^(1-alpha)). The paper contains no numerical correction of a seed.
- "The generating orbits (for mu = 0) associated to p-q resonant orbits are bifurcation orbits of 1st species-2nd species" (p.170, quoted): the nonexistence statement for p = q = 1 has a general analogue in the Gomez-Olle digest (passages scaling as mu^nu, nu < 1).
- Print slips noted above: Eq. 11 (psi_0 for phi_0), the cos a expression below Eq. 5, and the ambiguous bracketing of O_1 to O_3.

## 9. Checks
1. (45) equals the printed Casoliva (17) (identity 1/sqrt c = 2/sqrt(3 - C_J)). COMPUTED, symbolic.
2. (40) is consistent with the in-map shift (38). COMPUTED, section 6.
3. For the ten Casoliva 2008 Table 2 rows (mu = 1e-6) the two psi branches from (46), computed with the printed C_J and (p, q) in the 2003 labelling:
   12a: 55.3206, 124.6794 deg (corrected-orbit psi at the y = 0 crossing: 124.5372); 54a: 81.7788, 98.2212 (table 98.2195); 23b: 58.8306, 121.1694 (table 121.2811); 54b: 86.8206, 93.1794 (table 93.5930);
   73a: 80.7476, 99.2524 (table 99.8668); 32b: 40.1071, 139.8929 (table 137.8934); 21a: 88.8471, 91.1529 (table 90.0000); 32a: 87.1200, 92.8800 (table 90.0000); 23a and 52a: sin psi > 1 (C_J just below C_J1, by 1.3e-3 and 8.8e-4).
   (COMPUTED. The table state is a corrected orbit's y = 0 crossing at distance r_2 between 1e-4 and 2e-2 from the Moon, not a point on the mu^alpha sphere, so agreement to 0.1 to 2 degrees in psi is what is expected, and theta cannot be compared directly.)
4. Planar theta from (45) for those psi, with cos(theta - psi) > 0: one solution each (e.g. 12a branch psi = 124.6794 -> theta = 63.956; 54a branch psi = 98.2212 -> theta = 154.54), confirming "one admissible theta" (p.170).

## 10. Applying this to the project (see the Casoliva 2008 digest, section "Recipe for the P2 build", for the step list)
- Positive controls from this paper alone: the zero-order consistency (40); the four printed C_J inside the intervals; the identity (45) = Casoliva (17); the one-theta-per-psi count.
- The seed is only good to O(mu^(1-alpha)) and the radius mu^alpha is a free choice with alpha in (1/3, 1/2) (printed example 0.4); a seed generator must carry alpha as an explicit parameter and the corrector must not depend on it.
- For Casoliva's 7-3 cycler the seed (p, q) = (7, 3) has p/q = 2.333 > 1: two C_J values per psi_0, and the condition sin psi_0 >= sqrt((p/q)^(2/3) - 1) = sqrt(0.7592) (COMPUTED: (7/3)^(2/3) = 1.7592) = 0.8713, so psi_0 must lie in [60.6, 119.4] degrees; the 73a state's psi = 99.87 deg is inside.

## 11. Proposed corpus-index row
| barrabes-gomez-2003-three-dimensional-pq-resonant-orbits-second-species-solutions-cmda-85-2-doi-10.1023-A1022098510161.pdf | 2026-10-04-digest-barrabes-gomez-2003-pq-resonant-second-species.md | In-map (KS-regularised hyperbolic inner solution, Eq. 35-39) and matching with the out-map (Eq. 40-46) for p-q resonant orbits near second-species solutions at small mu; planar seed = (C_J, psi_0 from Eq. 46, theta from E = 0 Eq. 45, eps = delta = -cos a_0/sqrt(c)); spatial family phi_0 = +-pi/2, C_J = 2 + (p/q)^(2/3), p < q; p = q = 1 excluded; text layer (11102 words); full digest with checks. | DIGESTED (#920) |
