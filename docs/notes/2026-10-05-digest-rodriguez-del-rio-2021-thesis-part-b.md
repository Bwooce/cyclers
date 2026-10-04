# Digest: Rodriguez del Rio 2021 PhD thesis, part B (chapters 4 to 7, conclusions, bibliography)

Date: 2026-10-05 (Sydney). Reading, three small independent computations (section 9), no project code changed.

Source: O. Rodriguez del Rio, "Ejection-collision orbits in the Restricted Three-Body Problem", PhD dissertation, Universitat
Politecnica de Catalunya, June 2021 (advisors M. Olle and J. Soler), 180 PDF pages. Filed in the private paper corpus as
`rodriguez-del-rio-2021-ejection-collision-orbits-restricted-three-body-problem-phd-thesis-upc-doi-10.5821-dissertation-2117-351117.pdf`.
Page numbers below are the PRINTED thesis pages (PDF page = printed page + 14). Part A (front matter, chapters 1 to 3:
conventions, regularisation, numerical method) is in `docs/notes/2026-10-05-digest-rodriguez-del-rio-2021-thesis-part-a.md`.
Text layer is good; equations were read from the text layer and spot-checked on page images (printed pp. 95, 122, 125).
Two companion papers are compared in sections 11 and 12.

## 1. Which paper is which chapter (thesis p. 4, "this thesis is made up of")

| Thesis | Paper | Corpus status |
|---|---|---|
| Ch. 3 (numerical n-EC orbits) | ORS18: Olle, Rodriguez, Soler, "Ejection-collision orbits in the RTBP", CNSNS 55:298-315 (2018) | not held |
| Ch. 2 (regularisation) | ORS20b: "Regularisation in ejection-collision orbits of the RTBP", Recent Advances in Pure and Applied Mathematics, 35-47 (2020) | not held |
| Ch. 4 (analytical existence I) and part of ch. 3 | ORS20a: "Analytical and numerical results on families of n-ejection-collision orbits in the RTBP", CNSNS 90:105294 (2020), DOI 10.1016/j.cnsns.2020.105294 | held, section 11 |
| Ch. 5 (analytical existence II, Hill problem) | MSORS21: Martinez-Seara, Olle, Rodriguez, Soler, "Generalised analytical results on n-ejection-collision orbits in the RTBP. Analysis of bifurcations", listed as Preprint 2021 | not held (the coordinator's "2022, CNSNS 111:106410" is probably its published form; that identification is not printed in the thesis and is not verified here) |
| Ch. 6 (transit regions) | ORS21b: "Transit regions and ejection/collision orbits in the RTBP", CNSNS 94:105550 (2021), DOI 10.1016/j.cnsns.2020.105550 | held, section 12 |
| Ch. 7 (spatial case) | ORS21a: "McGehee regularization in the 3D restricted three-body problem. Application to ejection-collision orbits", listed as Preprint 2021 | not held |

The thesis bibliography prints no DOIs at all (journal, volume, pages only).

## 2. Conventions to carry over (thesis ch. 1, pp. 12-14, 61; checked)

- Primary P1 has mass 1 - mu and sits at (mu, 0); primary P2 has mass mu and sits at (mu - 1, 0). The ejecting primary is always P1,
  so mu in (0, 0.5] ejects from the big body and mu in [0.5, 1) from the small one. For the Earth-Moon system: Earth ejection is
  mu = 0.01215; Moon ejection is mu = 1 - 0.01215 = 0.98785.
- The project places the big primary at (-mu, 0). The thesis frame is the project frame rotated by pi about the barycentre,
  (x, y, vx, vy) to (-x, -y, -vx, -vy), which preserves time and the sense of rotation; it is not a one-axis mirror.
- Jacobi constant C = 2 Omega - v^2 with Omega = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2 + mu(1 - mu)/2 (eq. 1.22, 2). The thesis C
  therefore equals the project's standard C plus mu(1 - mu). Independent check (section 9): this reproduces the printed
  C_L1 and C_L2 values to 1e-15. At the Earth-Moon mass the thesis C_L1 is 3.2003381 (standard 3.1883357) and C_L2 is
  3.1841582 thesis (3.1721558 standard); the thesis H = -C/2.
- Hill problem: C = 3 + (1 - mu)^(2/3) K (eq. 1.19), K_L = 3^(4/3) = 4.3267487 (p. 12).
- Levi-Civita: x = mu + u^2 - v^2, y = 2uv, dt/ds = 4(u^2 + v^2); ejection initial condition (0, 0, 2 sqrt(2(1 - mu)) cos th0,
  2 sqrt(2(1 - mu)) sin th0), th0 in [0, pi) (double cover). The physical ejection direction is 2 th0 (derived from
  z = w^2, not printed; consistent with the printed statement that th0 tends to 0, pi/4, pi/2, 3pi/4, orbits lying on the axes).
- n-EC orbit: ejects from P1, reaches n relative maxima of distance from P1, then collides with P1 (so it has n - 1 non-collision
  close passages between). Families alpha_n, gamma_n (each symmetric about the x axis) and beta_n, delta_n (mirror images of one another).

## 3. Chapter 4, analytical existence I (pp. 67-85; source ORS20a)

Theorem 2 (p. 68). For C big enough and every n there exists mu_hat(C, n) such that for mu <= mu_hat(C, n) there exist EXACTLY four
n-EC orbits: two symmetric about the x axis, two mirror images of each other; families alpha_n, gamma_n, beta_n, delta_n.

Hypotheses actually proved: C "big enough" with no explicit bound, and mu small enough with no explicit bound; both
constants are existential (Implicit Function Theorem). Remark 1 of the paper says the same: the value of mu cannot be made explicit.

Method (pp. 68-77): expand the Levi-Civita system in mu, u = u0 + mu u1 + O(mu^2); solve the mu = 0 problem in closed form
(Kepler in sidereal Levi-Civita variables, eqs. 4.12-4.20, rotated to the synodic frame, Lemma 5); get u1 by variation of constants
(4.5). Lemma 3 (p. 69): for C large and mu small an ejection orbit is an EC orbit if and only if the angular momentum
M = u v' - v u' vanishes at a distance minimum.

Printed results usable as tests (all in thesis convention):

- Time of the n-th distance minimum at mu = 0: s0* = n pi / (2 sqrt(C)) (4.19), independent of th0; t(s0*) = 2 pi n / C^(3/2) (4.29); the closed form t(s) of (4.28) is not transcribed (text-layer layout ambiguous).
- Ejection-orbit solution at mu = 0: u^(s) = sqrt(2) cos th0 sin(2 sqrt(C) s) / sqrt(C), v^(s) = sqrt(2) sin th0 sin(2 sqrt(C) s)/ sqrt(C) (4.18).
- Result (4.39), p. 76: M_n(th0) = -15 mu n pi sin(4 th0) / C^(7/2) + mu O(C^(-9/2)) + O(mu^2); the roots are
  th0 = pi m / 4 + O(C^(-1), mu C^(7/2)) (4.41), m = 0..3; m = 0, 2 are the x-symmetric orbits, m = 1, 3 the mirror pair.
  Independent check (section 9): the coefficient, the sign and the linear dependence on n are all reproduced.
- Remark 9: the other proof route (curves D_i^+ and D_j^- with i + j = n + 1) needs higher order in C^(-1/2) and is harder; the
  angular-momentum route was chosen. The paper (section 11) uses the D_i^+, D_j^- route instead.

## 4. Chapter 5, analytical existence II and the Hill problem (pp. 87-113; source MSORS21, not held)

Theorem A (p. 87). There exists L_hat such that for every L >= L_hat, every mu in (0, 1), every n and
C = 3 mu + L n^(2/3) (1 - mu)^(2/3), there are EXACTLY four n-EC orbits (same characterisation).

Theorem B (p. 88). For each n there exists K_hat(n) such that for K >= K_hat(n), every mu in (0, 1) and
C = 3 mu + K (1 - mu)^(2/3), there are exactly four n-EC orbits. K_hat is uniform in mu; as mu tends to 1 the threshold C tends to 3,
as C_L1 does. The constants L_hat and K_hat(n) are not given numerically.

Method: scale u = sqrt(2(1 - mu)/(C - 3 mu)) U, v likewise, tau = 2 sqrt(C - 3 mu) s (5.1); expansion parameter
delta = 1 / sqrt(C - 3 mu) (the expansion of C_L1 suggested the scaling); order 6 is the first with a mu-dependent term (the
first perturbation of the two-body solution is O(delta^6)); same M_n = 0 criterion (Lemma 7, which removes the "mu small" proviso
of Lemma 3). Theorem A follows by the further scaling K = L n^(2/3), which makes the 1/n dependence uniform (Lemmas 8-10).

Higher-order angular momentum (eq. 5.26, p. 95), epsilon = 1 / sqrt(K), image-checked:

M_n(th0) = -(15 mu n pi sin 4th0 / 4) eps^6 + (105 mu (1-mu)^(1/3) n pi (sin 2th0 + 5 sin 6th0) / 64) eps^8
 + (15 mu n^2 pi^2 cos 4th0 / 2) eps^9 - (315 mu (1-mu)^(2/3) n pi (2 sin 4th0 + 7 sin 8th0) / 128) eps^10 + O(eps^11),

with the time corrections tau*_7 = 0, tau*_8 = -35 mu (1-mu)^(1/3) n pi cos 2th0 (5 cos^2 2th0 - 3)/4, tau*_9 = 15 mu n^2 pi^2 sin 4th0 / 2,
tau*_10 = 315 mu (1-mu)^(2/3) n pi (13 - 10 cos 4th0 - 35 cos^2 4th0)/256. The thesis reads the first bifurcation of
n-EC families (two new orbits, from alpha_n, for n not equal to 3) off the sin 6th0 term and the n = 3 bifurcation (four new orbits)
off the sin 8th0 term, with the caveat that this is qualitative.

Hill problem (section 5.5, pp. 108-113): the same proofs with mu = 1. Corollary 5.5.1: for K >= K_hat(n) there are exactly four
n-EC orbits (Hill problem has the extra y-axis symmetry, so alpha_n, gamma_n are mirror images over y and beta_n, delta_n are each
y-symmetric). Corollary 5.5.2: K = L n^(2/3), L >= L_hat. Numerics (printed):
- the numerical K_hat(n) follows the curve L n^(2/3) with L = 2^(2/3) (Fig. 5.6, p. 113, a plot, no table);
- successive bifurcation values of K follow L n^(2/3) with L = (2/p)^(2/3), p = 1..10 (Fig. 5.7);
- n = 5: new families born at K = 5.02714993, eight 5-EC orbits at K = 4.86, collapse at K = 4.72835275 (Fig. 5.4, p. 111);
- n = 9: periodic EC orbit (alpha_9, gamma_9) at K = 4.77318771, (beta_9, delta_9) at K = 4.42215362 (Fig. 5.5, p. 112).
The authors flag two open points (Conclusions, p. 161): the radius of convergence of the series, and why the Hill bifurcations
follow (2/p)^(2/3).

## 5. The mu and C ranges, and the Earth-Moon mass (the question asked)

What is PROVED (thesis theorems):

- Theorem 2 (ch. 4): mu <= mu_hat(C, n), C >= C_hat, both existential. The Earth-Moon mu = 0.01215 is small but there is no
  statement that it lies below mu_hat for any given C and n. Not applicable as a proof at the Earth-Moon mass.
- Theorems A and B (ch. 5): ALL mu in (0, 1), so mu = 0.01215 and mu = 0.98785 are inside the stated range, but only for
  C >= 3 mu + L_hat n^(2/3) (1 - mu)^(2/3) or C >= 3 mu + K_hat(n)(1 - mu)^(2/3) with L_hat, K_hat(n) unknown. The theorem is
  therefore silent on whether any physically interesting C (near C_L1 = 3.2003) is covered.

What is NUMERICAL (thesis ch. 3, part A; re-listed here because it decides the answer):

- n = 1: four 1-EC orbits exist numerically for every mu in (0, 1) and every C >= C_L1(mu) (p. 43: "we can extend numerically
  Theorem 1 ... for values of C >= C_hat(mu) = C_L1(mu)"). So for n = 1 the Earth-Moon mass, both ejection from the Earth
  and from the Moon, is covered down to C_L1 = 3.2003381 (thesis) = 3.1883357 (standard) by computation, not proof.
- n >= 2: C_hat(mu, n) is computed only as a plot (Fig. 3.27, mu in (0, 1), n = 2..10; no table) and as two printed values at mu = 0.1:
  C_hat(0.1, 2) = 3.72442505 and C_hat(0.1, 3) = 3.80644009 (pp. 52-53), against C_L1(0.1) = 3.68695322987989, so for n >= 2 the
  four-orbit regime ends ABOVE C_L1 at mu = 0.1. Fig. 3.27 shows C_hat rising with n, and tending to 3 as mu tends to 1.
  No value at mu = 0.01215 is printed.
- Derived estimate (mine, not printed, order of magnitude only): taking the Hill-limit fit K_hat(n) approximately 2^(2/3) n^(2/3)
  at mu = 0.98785 (Moon ejection) gives C_hat(thesis) of 3.0474, 3.0967, 3.1381, 3.1750, 3.2089, 3.2406 for n = 1..6,
  i.e. about 0.012 less in standard C (3.0354 ... 3.2286). Against the thesis C_L1 of 3.2003, this says four n-EC orbits to the Moon
  persist down to C_L1 for n up to 4, and bifurcations begin above C_L1 from n = 5. The Hill fit is an extrapolation from
  mu tending to 1 (my check: the same mapping at K_L gives C = 3.1922 against the true C_L1 of 3.2003, 0.25 percent off), not a tested statement
  at 0.98785. For Earth ejection (mu = 0.01215) no estimate is available from the thesis.

Consequence: the Earth-Moon mass lies inside the proved range of mu, but the proved range of C is unquantified; the usable
statement is the numerical one (n = 1, all C >= C_L1; n >= 2 above a threshold that must be computed).

## 6. Chapter 6, transit regions (pp. 115-141; source ORS21b)

Range: C in [C_L2,3, C_L1) with C_L2,3 = C_L2 for mu in (0, 0.5] and C_L3 for mu in [0.5, 1) (p. 115); bounded motion in both
regions with a neck at L1, so ejection orbits can pass between the P1 and P2 regions. For the Earth-Moon mass this window is
thesis C in [3.1841582, 3.2003381), width 0.016 (standard C 3.1722 to 3.1883); below it L2 opens (the thesis excludes it).

Results:

- Heteroclinic connections P_i - LPO1 (ejection primary to the L1 Lyapunov orbit): found as W^e(P1) intersect W^s(LPO1) on a
  section x = x_L1 + d (d = 0.1 in the figures). For n close passages before the section they come in pairs H_n^1, H_n^2
  bounding a transit interval I_n = (th_n^1, th_n^2); every ejection angle inside I_n gives an orbit that leaves the P1 region for
  the P2 region, outside gives non-transit orbits that return. Found for n = 0, 2, 3 and NEVER for n = 1 (scanned over all mu
  in (0, 1) and all C in [C_L2,3, C_L1), Fig. 6.5, p. 121).
- Existence of I_n depends on (mu, C, n): for given n there is a minimum and maximum mu, and for given mu a window of C, ending
  at tangency of the two curves (a single connection) with none above (Figs. 6.4, 6.5; plots only, no table).
- Chaos: infinitely many homoclinic orbits to LPO1 (two simple ones HO1, HO2 drawn) generate infinitely many P1 - LPO1
  connections near each of the two basic ones, hence an infinity of orbits E_n^i - PO_k - P_l^j (eject, n passes, k turns about LPO1,
  visit P_l, j passes). The first and second classification described around th_0^1 and th_0^2 is the printed example.
- EC orbits through LPO1: the same mechanism gives infinitely many EC orbits to the same primary (two ejection roads times two
  collision roads, k = 0, 1, 2, ... turns) and, with the double spiral of transit orbits, to the OTHER primary (four countable sets;
  Fig. 6.14); coded E_1^0 - PO_k - C_l^j.
- Colour-code diagrams (Figs. 6.18, 6.19; t in [0, 10]): at mu = 0.5 the first transit tongues appear at C = 4 (two thin tongues,
  after six passes then one turn about P2 and back); n = 0 tongue visible at C = 3.85; more at C = C_L2. At C = C_L2 (mu = 0.2 to 0.9)
  there is little transit at mu = 0.2 over t in [0, 10], a main tongue near th0 = pi/2 grows with mu; I_5 is visible at mu = 0.3.
- Time to the k-th minimum (Fig. 6.15, mu = 0.2, C = 10, 5, 4.2, 3.8): ripples grow as C falls and k rises; theta_0 near pi/2
  (velocity toward P2) shows the longest times.

Relevance to Earth-Moon cyclers: this is the only part of the thesis that treats passage between the Earth and Moon regions, and
it does so for ejection (zero-distance start) and for a narrow band of C. Cyclers pass near, not through, the primaries and mostly
sit at C below C_L2 (Casoliva Class 2 spans standard C 2.89 to 3.18 per the project's Casoliva digest), so the printed
transit intervals do not transfer. The transferable object is the mechanism: the stable and unstable manifolds of the L1 Lyapunov orbit
partition orbits leaving a primary region into transit and non-transit sets, with heteroclinic connections as boundaries.

## 7. Chapter 7, the spatial case (pp. 143-159; source ORS21a, not held)

- Spatial McGehee blow-up (spherical coordinates (r, theta, phi, v, u_theta, u_phi)); time rescalings dt = r^(3/2) dtau fails on the
  z axis (as in Llibre and Martinez Alfaro 1985, whose analysis ignored this set), so dt = r^(3/2) cos(phi) dtau-hat is used, which leaves
  the whole set cos(phi) = 0, u_theta = 0 as equilibria and the z axis invariant; three local charts (about the z, x and y axes) cover
  all directions. Kustaanheimo-Stiefel is not used (future work).
- Collision manifold: v^2 + u_theta^2 + u_phi^2 = 2(1 - mu), diffeomorphic to S^2 x S^2; two spheres of equilibria S^2_+- with
  v0 = +-sqrt(2(1 - mu)); linearisation eigenvalues lambda, lambda, -lambda/2, -lambda/2, 0, 0 with lambda = v0 cos(phi) (eq. 7.15);
  dim W^s(S^2_+) = dim W^u(S^2_-) = 4, dim W^u(S^2_+) = dim W^s(S^2_-) = 3 at fixed energy (Proposition, from Llibre and Martinez Alfaro).
- Detection: EC orbits are heteroclinic connections; candidates are the zeros of the angular-momentum VECTOR direction at the first
  section crossing, found as intersections of curves A (jump alpha to alpha + pi) and D (jump delta to -delta). The authors stress
  these are candidates because the planar transversality argument does not carry over.
- Result at mu = 0.1, C = 4 (H = -2): four planar EC orbits (two x-symmetric, two mirror images) as proved in ch. 5, plus two
  non-planar transversal EC orbits (each symmetric about the xz plane, the pair mirror images across the xy plane), plus an apparently
  continuous set of candidates (curves A2 and D1, not shown to coincide). With increasing C the planar angles tend to 0, pi/2, pi, 3pi/2
  and the two spatial ones tend to (theta0, phi0) = (pi, +-pi/2) (Fig. 7.9; mu = 0.1 and 0.5, C = C_L1, 4, 6, 8, 10, 20, 50, 100).
- Contradiction with the literature (p. 158): Llibre and Martinez Alfaro 1985 claim families of spatial symmetric EC orbits
  parametrised by an angle; the thesis finds the planar symmetric orbits isolated and not continuing to spatial families.
  Printed number: C_L1(mu = 0.1) = 3.68695322987989 (p. 158).
- Only n = 1 is treated; n-EC in 3D is future work.

## 8. Conclusions (p. 161)

Four n-EC orbits for C >= C_hat(mu, n), the shape of the analytical threshold matches the numerical one; the method generalises to any
Newtonian-type perturbed two-body problem (regularise the ejecting body, rewrite as two-body plus perturbation); open items: radius of
convergence, the Hill bifurcation pattern, interaction with the other equilibria (L4, L5 already seen to matter, ch. 3) and the
manifold of infinity, periodic EC orbits, n-EC in 3D with KS regularisation.

## 9. Test-ready numbers (thesis convention: C includes mu(1 - mu); P1 mass 1 - mu at +mu)

Status column: "indep" means I reproduced it with a throwaway integrator (not in the repo); "printed" means taken as printed.
Nothing here is yet a repo test.

| # | Quantity | Value | Page | Status |
|---|---|---|---|---|
| 1 | C_L1(mu = 0.5), C_L2 = C_L3 | 4.25 ; 3.7067962240861525 | p. 138 | indep: 4.25 ; 3.706796224086153 |
| 2 | C_L1(mu = 0.1) | 3.68695322987989 | p. 158 | indep: 3.6869532298798946 |
| 3 | K_L (Hill) | 3^(4/3) | p. 12 | exact |
| 4 | Ejection speed in Levi-Civita | |u'|^2 + |v'|^2 = 8(1 - mu) | p. 21 (eq. 2.29) | used in 6 |
| 5 | M_n(th0) leading term | -15 mu n pi sin(4 th0) / C^(7/2) | p. 76 (4.39) | indep: ratio numerical/predicted 0.9879 (C = 200, n = 1), 0.9877 (n = 3), 0.9599 and 0.9569 (C = 50), tending to 1 as C grows; sign and n-linearity correct |
| 6 | Two-body Levi-Civita times | s0* = n pi / (2 sqrt C) | p. 71 | indep: at C = 8, mu = 0.1 the extrema fall at s = 0.2835, 0.5670, 0.8506, 1.1341 against the mu = 0 values 0.2777, 0.5554, 0.8331, 1.1107 (shifted by the mu perturbation, as expected) |
| 7 | 2-EC angles, mu = 0.1, C = 8 (paper Table 1 only, not in the thesis) | alpha2 1.7037895180958; beta2 2.44063260106480; gamma2 0.1327163661903; delta2 0.9650603806965 | ORS20a p. 16 | indep: each gives a collision at the 4th distance extremum, r^2 about 1e-28 and M about 1e-14; a 1e-3 offset in alpha2 gives r^2 1.5e-10 (control) |
| 8 | C_hat(0.1, 2), C_hat(0.1, 3) | 3.72442505 ; 3.80644009 | pp. 52-53 (also ORS20a p. 19) | indep bracket: 4 two-EC orbits at C = 3.76, 3.74; 6 at C = 3.70; exact value not yet located (section 13) |
| 9 | Hill n = 5 | new families at K = 5.02714993, 8 orbits at 4.86, collapse 4.72835275 | p. 111 | printed |
| 10 | Hill n = 9 periodic EC | K = 4.77318771 (alpha, gamma), 4.42215362 (beta, delta) | p. 112 | printed |
| 11 | Hill fit | K_hat(n) near 2^(2/3) n^(2/3); bifurcation curves (2/p)^(2/3) n^(2/3), p = 1..10 | pp. 112-113 | plots only |
| 12 | Transit interval n = 0, mu = 0.5, C = C_L2 | th_0^1 = 1.558674225724; th_0^2 = 1.932752613334 | p. 122 | printed (also ORS21b) |
| 13 | New heteroclinic angles near them | th_0^{1,1} = 1.584919597063; th_0^{1,2} = 1.702359169703; th_0^{2,1} = 1.867534417000; th_0^{2,2} = 1.925408730327 | p. 122 | printed, image-checked; the sentence attaching 1.702359169703 is garbled ("close to th_0^1 = 1.7023...") and 1.7024 is not close to th_0^1 |
| 14 | Homoclinic-mediated angles | 1.553258226788; 1.557995153267 (near th_0^1); 1.932919335638; 1.934043834973 (near th_0^2) | p. 125 | printed, image-checked |
| 15 | n = 1 transit | no P_i - LPO1 connection exists for n = 1, any mu in (0,1), C in [C_L2,3, C_L1) | p. 121 | printed claim, a negative result under a stated scan |
| 16 | Spatial EC count, mu = 0.1, C = 4 | 4 planar + 2 spatial transversal + continuous candidate set | p. 157 | printed |
| 17 | mu = 0.01, n = 1 continuation (part A, Fig. 3.30, p. 56) | C = 2.472170770645, 1.970463731686, 1.970412419219, 1.970412407739 | p. 56 | printed |

Items 12-14 are for mu = 0.5, C = C_L2 and are ejection angles in the Levi-Civita double cover; a repo test would need the Lyapunov-orbit
manifold intersection (the numbers cannot be reproduced without it) and are therefore best used as a late control. Items 5, 7 and 8
need only a Levi-Civita integrator.

## 10. Techniques applicable to the project's problems

#899 (Earth-Moon cyclers, Casoliva families):

- Do the EC orbits give a skeleton for Earth-Moon cyclers? Partly. The mu = 0 skeleton of every second-species flyby orbit is a
  collision (ejection) orbit of the two-body problem; Henon 1968, Bruno 1981 and Casoliva Class 1 are built on exactly that seed. The thesis
  adds three things: (a) for large C, which n-EC orbits survive for mu above 0: four, at theta0 approximately pi m / 4 (so a seed
  list at C large is four orbits per n, not a continuum); (b) the closed expansion (4.39) and (5.26), which gives the first-order
  residual M_n(th0) as a function of mu, C, n and thereby a quantitative seed-survival rule; (c) a full numerical C_hat(mu, n) map
  (not tabulated) saying where the four-orbit regime ends.
- Limits: all the existence results are for C above C_L1, whereas the Casoliva Class 1 cyclers have C_J 0.49 to 2.76 (standard C,
  per the project note), far below C_L2. So the proved EC skeleton does not reach them; it reaches only the seeds near
  C_L1 = 3.188 standard, and the second species orbits are periodic with collision-near passages, not collisions. A direct EC seed
  for them would need continuation in C from C_L1 downwards through the L1 and L2 openings, where the thesis offers only the
  chaotic classification of ch. 6 and no existence theorem.
- Can they seed the Casoliva families? As starting points in mass (the Casoliva method continues in mu): the thesis shows that EC
  orbits at mu small exist only as isolated points (four per n), so a mu-continuation from
  1e-6 starts from those, not from a family. This is consistent with Casoliva's report that continuation at fixed period ends on lunar impact; it
  does not supply a new strategy. The thesis (pp. 55-56, part A) also shows that along C-continuation alpha_n and gamma_n approach a
  periodic orbit without merging, and the definition of n-EC is local; so a seed list should be classified by the orbit's own
  maxima and minima count, not by n alone.
- Strongest concrete use: M_n(th0) at the n-th distance minimum is a one-dimensional residual with simple roots, so a
  shooting-based seed finder for collision orbits can be validated against items 5, 7 and 8 of section 9, which the project can reproduce with
  a Levi-Civita integrator (cf. #928).

#906 (demanded-turn gate):

- An ejection-collision orbit is the zero-angular-momentum limit of a flyby: periapsis 0, demanded turn 180 degrees, so an exact
  collision is not a flyby and must not be scored as one. The thesis' Lemma 3 shows (C large, mu small) that M = 0 at a distance
  minimum is equivalent to collision. The near-collision middle passes of the 2-EC orbits in section 9 item 7 (r^2 about 1e-7, M about
  1e-3 to 1e-4 in Levi-Civita units) are real orbits with a tiny periapsis: a ready test that the gate returns indeterminate there,
  not a rejection, as the #906 amendment says.

#928 (regularised propagator):

- Chapters 2 and 4 are a worked, tested Levi-Civita system for the CR3BP in the rotating frame (eq. 4.1 = 2.26, with the printed
  Jacobi invariant u'^2 + v'^2 = 8(u^2 + v^2) U). I transcribed eq. 4.1 and integrated it: the physical Jacobi constant recovered
  from (u, v, u', v') is conserved to about 5e-13 at C = 8, mu = 0.1, so the transcription (including the 16 u^3 and -16 v^3 terms and the
  r2 = sqrt((1 + u^2 - v^2)^2 + 4u^2 v^2) terms) is correct as printed in the text layer. This gives the "first #928 test": integrate
  the four Table 1 orbits (section 9, item 7) and expect collisions to 1e-28.
- The regularisation is LOCAL (it removes only the ejecting primary), so the Moon would need the second chart; the thesis (and
  ORS21b section 3) integrates in synodic coordinates away from the primaries and switches charts, carrying three times; this matches
  the project's Birkhoff-global plan only partly.
- The Hill-region exclusion at C >= C_L1 is what makes one local chart sufficient; the project's cyclers (below C_L2) do not satisfy it.
- The McGehee 3D chart (ch. 7) is documented as numerically poor (ejection and collision are asymptotic; time rescalings stall at
  the z axis); the authors' own recommendation is Kustaanheimo-Stiefel for 3D, which the project holds (Stiefel and Scheifele 1971;
  Kustaanheimo and Stiefel 1965, both in the corpus).

## 11. Olle, Rodriguez & Soler 2020 compared with the thesis

Source: Olle, Rodriguez and Soler, CNSNS 90 (2020) 105294, DOI 10.1016/j.cnsns.2020.105294, received 27 Nov 2019, accepted 13 Apr 2020,
22 journal pages. Filed in the private paper corpus as
`olle-rodriguez-soler-2020-families-n-ejection-collision-orbits-rtbp-cnsns-90-105294-doi-10.1016-j.cnsns.2020.105294.pdf`.
Text layer good (equation digits in the paper read from the text layer; Table 1 digits read from the text layer, and verified numerically).

Correspondence: sections 2-3 and 4-6 of the paper are thesis chapter 4 (with chapter 2 for the regularisation); paper section 7 is part of
thesis chapter 3 (n-EC numerics, C_hat). The thesis chapter 5 (Theorems A and B, any mu, scaled variables) is NOT in this paper.

Differences:

1. Theorem statement. Paper Theorem 1 says "there exist four n-EC orbits" for mu <= mu_hat(C, n) (no word "exactly"); the thesis Theorem 2
   says "exactly four". The paper's abstract and section 7 do say "exactly four" and the Lemma 6 proof does count "four and only four", so the
   weaker wording in the theorem is a drafting difference, not a different result.
2. Range of mu. Both the paper's Theorem 1 and the thesis Theorem 2 are small-mu results with no explicit mu_hat or C_hat. The paper's section 2
   fixes mu in (0, 0.5] and shows the small-primary case (mu in [0.5, 1)) only numerically in section 7; the thesis states the any-mu result
   (Theorems A and B, ch. 5) separately, which the paper does not contain.
3. Method. The paper's proof uses the intersections of the curves D_i^+ and D_j^- (kth maxima of distance on the section g = u u' + v v' = 0,
   g' < 0; Lemma 1: |D_i^+ cap D_j^-| = 2 (i + j - 1)-EC orbits; Lemma 6) with Lemmas 2-5; the thesis ch. 4 proves it by the angular momentum at the n-th
   MINIMUM (Lemma 3) and says (Remark 9) it is simpler. Hence different printed expansions:
   - paper Lemma 2: s0k = (2k - 1) pi / (4 sqrt C) (time of the k-th maximum); thesis (4.19): s0* = n pi / (2 sqrt C) (time of the n-th minimum).
   - paper Lemma 3: R0k^2 = 2/C (square of the maximum Levi-Civita distance); not in the thesis.
   - paper Lemma 4 (eq. 60): R1k^2(th0) has leading terms of order 1/C and 1/C^2, then terms containing (3 cos 4th0 + 1), 8 cos 2th0 and
     (3 - 5 cos^2 2th0) at successive orders, and its last printed term is 2(2k - 1) pi sin(4 th0) / C^(11/2), remainder O(C^(-13/2)); the
     exact powers and coefficients are not transcribed here (text layer ambiguous). Not in the thesis.
   - paper Lemma 5: th_k = th0 - (2k - 1) pi / (2 C^(3/2)) + O(mu).
   - paper eq. (34): zero-finding function -4 pi (5i + 7j - 6) sin(4 th) / C^(3/2) + O(C^(-5/2)); the thesis (4.39) is -15 mu n pi sin(4 th0)/C^(7/2).
     The two share the factor sin 4th0 (hence the same four roots pi m / 4) but the paper's coefficient depends on the split (i, j) of
     n + 1, the thesis on n only; they are different functions (radius difference versus angular momentum), so there is no inconsistency,
     but the paper's coefficient has not been checked numerically here.
4. Numbers only in the paper (not in the thesis): Table 1 (section 9 item 7), the ejection speed written as 2 sqrt(1.8) for mu = 0.1,
   Figure 4 (mu = 0.1, n = 2, C = 8, curves D_2^+ and D_1^-), Figure 7 (n-EC orbits for mu = 0.1: n = 15 at C = 15, n = 25 at C = 20),
   Figure 6 (mu = 0.4, n = 3, C = 4.5; mu = 0.9; n = 1, 2, 5), Figure 5 (3D plot of D_2^+, D_1^- versus C for mu = 0.1), Figure 1 (Hill region
   mu = 0.2, C = 3.53 in both coordinate systems). The paper states "n = 1..25" for the earlier numerical literature ([18], [19]).
5. Numbers in both: C_hat(0.1, 2) = 3.72442505, C_hat(0.1, 3) = 3.80644009 (paper p. 19, thesis pp. 52-53); the paper plots C_hat(0.1, n) for n = 2..20
   (Fig. 10) and says C_hat(mu, 1) is below C_L1; the thesis adds the full mu map (Fig. 3.27, n = 2..10) and the bifurcation discussion.
6. Not in the paper: Hill problem, Theorems A and B and the scaled proof, the bifurcation expansion (5.26), the spatial case, transit regions.
7. Conventions: identical to the thesis (primary masses 1 - mu at (mu, 0) and mu at (mu - 1, 0); Jacobi constant with the mu(1 - mu)/2 term).

## 12. Olle, Rodriguez & Soler 2021 compared with the thesis

Source: Olle, Rodriguez and Soler, CNSNS 94 (2021) 105550, DOI 10.1016/j.cnsns.2020.105550, received 26 Apr 2020, revised 23 Aug 2020,
accepted 24 Sep 2020; the held copy is the journal pre-proof (41 pages, "Journal Pre-proof" watermark, not the typeset version). Filed in the private
paper corpus as
`olle-rodriguez-soler-2021-transit-regions-ejection-collision-orbits-rtbp-cnsns-94-105550-doi-10.1016-j.cnsns.2020.105550.pdf`.
Text layer good; figures not individually read (figure-heavy).

Correspondence: the paper is thesis chapter 6, with its sections 2-3 (RTBP, Levi-Civita) taking the place of the thesis chapters 1-2 material,
sections 4, 5, 6 equal thesis 6.1, 6.2, 6.3, section 7 the conclusions. Everything the thesis prints in chapter 6 numerically appears in the paper:
I compared the numeric tokens of the two texts (decimals) and the only decimals that differ are cross-reference numbers (thesis section and figure
numbers) and the paper's "3.53" (Figure 3 and Figure 4 captions: mu = 0.2, C = 3.53 Hill region in Levi-Civita coordinates, and mu = 0.2, C = 3.8
integration regions). Identical in both: the angles 1.558674225724, 1.932752613334, 1.584919597063, 1.702359169703, 1.867534417000,
1.925408730327, 1.553258226788, 1.557995153267, 1.932919335638, 1.934043834973; C_L1 = 4.25 and C_L2 = 3.7067962240861525 at mu = 0.5;
the mu = 0.5, 0.7 transition-interval figures; mu = 0.2 time-to-minimum figure with C = 10, 5, 4.2, 3.8; the colour diagrams
(mu = 0.5, C = C_L1, 4.1, 4, 3.85, C_L2; C = C_L2, mu = 0.2, 0.3, 0.7, 0.8, 0.9). The garbled attachment of 1.702359169703 to "theta_0^1" is also in the paper.

Differences:

1. The paper regularises either primary in one setup: ejection from P1 or P2 by initial velocity 2 sqrt(2(1 - b)) with "b = mu or mu - 1
   for P1 or P2" (section 3.2, text-layer reading). For P2 this gives 1 - b = 2 - mu, whereas the collision speed at P2 should involve mu
   (mass mu); this looks like a printed slip (it should read 1 + b for P2 or the mass directly); not checked on a page image. The thesis
   fixes P1 and relies on the symmetry mu goes to 1 - mu.
2. The paper integrates three different charts and three clocks (synodic t, Levi-Civita s about P1, s about P2) and says so explicitly
   (a 6-dimensional state with two additional times); the thesis chapter mentions this less.
3. The paper's introduction cites more of the application literature, and the paper has Figures 3 and 4 (Hill region and the chart regions at
   mu = 0.2) that the thesis chapter lacks. The thesis adds a final colour-diagram discussion of "I_5" (also in the paper).
4. No table in either; the thesis chapter 6 has no number that the paper lacks.

## 13. Bibliography entries the project would need that are not in `CORPUS_INDEX.md`

Checked by grep of the index and the notes directory on 2026-10-05 (presence of a digest or file name). The thesis prints no DOIs; DOIs
below are therefore marked "not printed" and none is supplied from memory.

Directly relevant to #899, #906, #928 and the ejection-collision thread:
- Chenciner and Llibre 1988, "A note on the existence of invariant punctured tori in the planar circular restricted three-body problem", Ergodic Theory Dyn. Syst. 8:63-72 (four EC orbits for any mu; the Levi-Civita plus McGehee proof). Not in the index.
- Lacomba and Llibre 1988, "Transversal ejection-collision orbits for the restricted problem and the Hill's problem with applications", J. Differential Equations 74:69-85. Not in the index.
- Llibre 1982, "On the restricted three-body problem when the mass parameter is small", Celest. Mech. 28:83-105. Not in the index.
- Delgado Fernandez 1988-89, "Transversal ejection-collision orbits in Hill's problem for C >> 1", Celest. Mech. Dyn. Astron. 44:299-307. Not in the index.
- Llibre and Pinol 1990, "On the elliptic restricted three-body problem", Celest. Mech. Dyn. Astron. 48:319-345; Pinol 1995, 61:315-331. Not in the index.
- Maranhao and Llibre 1998; Alvarez-Ramirez and Vidal 2013 (binary collision in the planar restricted (n+1)-body problem, Physica D 254:1-11). Not in the index.
- Henon 1965 ("Exploration numerique du probleme restreint I. Masses egales", Ann. Astrophys. 28:499-511) and Henon 1969 (Hill's case, Astron. Astrophys. 1:223-238): the index's one Henon hit is the 1997 book. Not in the index (the thesis cites the 1965 paper for EC orbits along periodic families at mu = 0.5).
- Bozis 1970, "Sets of collision periodic orbits in the restricted problem" (in Giacaglia, ed., Periodic orbits, stability and resonances, Springer, 176-191). Not in the index.
- Nagler 2004 and 2005, "Crash test for the Copenhagen problem" and "... the restricted three-body problem", Phys. Rev. E 69:066218 and 71:026227 (crash-probability maps). Not in the index.
- Barrabes, Mondelo and Olle 2013, "Numerical continuation of families of heteroclinic connections between periodic orbits in a Hamiltonian system", Nonlinearity 26 (ORS21b's method for P_i - LPO1 connections). Not in the index (the index holds the 2009 and 2009b papers by the same authors).
- Gomez, Llibre and Masdemont 1988-89, "Homoclinic and heteroclinic solutions in the restricted three-body problem", Celest. Mech. Dyn. Astron. 44. Not in the index.
- Paez and Guzzo 2020, "A study of temporary captures and collisions in the circular restricted three-body problem with normalizations of the Levi-Civita Hamiltonian", Int. J. Non-Linear Mech. 120:103417. Not in the index.
- Astakhov, Burbanks, Wiggins and Farrelly 2003, "Chaos-assisted capture of irregular moons", Nature 423:264-267. Not in the index.
- Broucke 1965, "Regularizations of the plane restricted three-body problem", Icarus 4. Not in the index (Broucke 1968 and 1969 are held).
- Conley 1963 ("On some new long periodic solutions of the plane restricted three body problem", CPAM 16): one index hit exists, but it is a mention, not a held file.
- McGehee 1974 (Invent. Math. 27:191-227) and 1978 (ICM Helsinki 827-834): mentioned in index rows, not held.
- Meyer and Schmidt 1982 (Hill's lunar equations); Wintner 1930; Estes and Lancaster 1968 (NASA TM, power-series Thiele-Burrau); Lemaitre 1955 (global regularisation): not held.
- ORS18 (CNSNS 55:298), ORS20b, MSORS21 (the Hill and any-mu theorems of ch. 5) and ORS21a (spatial McGehee): the two papers whose results are used here are the ones NOT held; MSORS21 would be the primary source of Theorems A and B.

## 14. Follow-ups (no task numbers registered)

1. Acquire MSORS21 (Theorems A and B, bifurcation analysis); check whether it states explicit L_hat or K_hat(n); verify the mapping to the CNSNS 111:106410 reference.
2. Build a small planar Levi-Civita integrator in `src/` (this feeds #928) and add the thesis items 1, 2, 5, 7 of section 9 as sourced tests; locate C_hat(0.1, 2) and C_hat(0.1, 3) by bisection on the number of M_n roots to the printed 9 digits (my bracket is 3.70 to 3.74 so far).
3. Compute C_hat(mu, n) at mu = 0.01215 and 0.98785 for n = 1..10 (the thesis gives only a plot, Fig. 3.27); this is the number the Earth-Moon question needs. Test the Hill-limit estimate of section 5 against it.
4. Decide, for #899, whether to use the n-EC seeds at C near C_L1 (standard 3.188) as the high-energy end of a C-continuation toward Casoliva Class 2 (C 2.89 to 3.18); this needs the L2 opening, which the thesis does not treat.
5. Reproduce the transit-interval angles (items 12-14) with the Lyapunov-orbit manifolds at mu = 0.5, C = C_L2, as a late-stage control for any Lyapunov-manifold tooling; then compute the mu = 0.01215 transit interval for n = 0 (and n = 2, 3) over C in [3.1722, 3.1883] standard, the only band the thesis theory covers.
6. Confirm on a page image the ORS21b ejection-speed formula for P2 (b = mu - 1) and report it if it is a slip.
7. Check the ORS20a eq. 34 coefficient 4 pi (5i + 7j - 6) numerically against a Levi-Civita integration.
8. Revisit the 3D claim that planar symmetric EC orbits do not continue to spatial families (Llibre and Martinez Alfaro 1985 digest in the index), since the project's 3D families may care.
9. A note for the part A digest and #906: thesis C includes mu(1 - mu); `core.cr3bp.jacobi_constant` does not; any test lifted from this thesis must convert.
