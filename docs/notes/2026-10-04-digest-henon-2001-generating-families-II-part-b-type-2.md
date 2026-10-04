# Digest: Henon 2001, "Generating Families in the Restricted Three-Body Problem II: Quantitative Study of Bifurcations", Part B (chapters 17 to 23: type 2 bifurcations)

Springer, Berlin, Lecture Notes in Physics Monographs m65, 2001, DOI 10.1007/3-540-44712-1. Filed in the private paper corpus as
`henon-2001-generating-families-restricted-three-body-problem-II-quantitative-study-bifurcations-lnp-m65-springer-doi-10.1007-3-540-44712-1.pdf`
(308 PDF pages, text layer usable for prose, garbled for equations). This note covers chapters 17 to 23, the Index of Notations and the reference list. Front matter and chapters 11 to 16 (definitions, general equations, type 1) are digested in `docs/notes/2026-10-04-digest-henon-2001-generating-families-II-part-a-type-1.md` (written in parallel by another agent). Volume I is digested in `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`; related digests: `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`, `docs/notes/2026-10-04-digest-perko-1976-second-species-O-mu-near-moon.md` (Perko 1976, 1981 and 1977 as Parts A, B, C), `docs/notes/2026-10-04-digest-guillaume-1973-periodic-symmetric-solutions-small-mu.md`, `docs/notes/2026-10-04-digest-guillaume-1975-extension-breakwell-perko-matching.md`, `docs/notes/2026-10-04-digest-guillaume-1975a-linear-analysis-second-species.md`, `docs/notes/2026-10-04-digest-bruno-1981-periodic-flybys-of-the-moon.md`, `docs/notes/2026-10-04-digest-brjuno-1978-periodic-solutions-arcs-mu-0.md`, `docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`.

Evidence tags: READ (printed page; every formula and table I use was also read on the page image); COMPUTED (my own arithmetic or numerical solution on 2026-10-04, scratch scripts not kept in the repository); INFERRED (my reading across sources, not stated by the book). The nine-step general method and the symbols (mu, Delta C, nu, O-notation) are chapter 11, which I read only for notation.

## 0. Page numbers

Printed pages are used below. PDF page minus printed page, COMPUTED from the running heads of every page: 11 for printed 149 to 181 (chapter 17 opens on PDF 160 = printed 149); 10 for printed 182 to 199; 9 for printed 200 to 271; 8 for printed 272 to 283; 7 for printed 284 to the end (Index of Definitions printed 297 = PDF 304). Printed ranges: chapter 17 pp.149-179, 18 pp.181-197, 19 pp.199-224, 20 pp.225-238, 21 pp.239-270, 22 pp.271-282, 23 pp.283-296; Index of Notations pp.299-300; references pp.301-302 (PDF 308-309).

Read in full (prose from the text layer, formulas and tables from images where used): chapter 17 entire; chapter 18 entire; chapter 19 sections 19.1.1 to 19.1.5, 19.2, 19.3 (opening) and 19.3.3, 19.4; chapter 20 sections 20.1, 20.2.2 to 20.2.5, 20.3 (opening); chapter 21 sections 21.1 (opening and 21.1.1), 21.2, 21.3.2.6 to 21.4; chapter 22 sections 22.2 and 22.3; chapter 23 entire. Skimmed: 19.3.1 to 19.3.2 and 21.3.1 (the branch-ordering proofs), 20.2.4 beyond the symmetric solutions, 21.1.2 to 21.1.5 and 21.3.2.1 to 21.3.2.5 (variational equations and the case-by-case junction pictures), 22.1. These skimmed parts feed the junction tables and are summarised only through their stated results.

The book is, as the Preface of volume I says, a heuristic monograph. Chapters 17 to 23 are formal asymptotics: balance of orders of magnitude, implicit function theorem on the scaled equations, then explicit solution or numerical solution of the scaled equations. Existence is argued by a non-vanishing Jacobian of the asymptotic system at each step (a proof sketch, not a theorem), and the junction results for n of 5 or more rest on numerical tracking of branches (22.2.4) and on a topological "positional method".

## 1. The problem treated, and what the symbols mean

(READ chapter 11.4, pp.14-16, for the notation only.) A type 2 bifurcation orbit of the planar circular restricted problem at mu = 0 is a Keplerian ellipse tangent to the Moon's orbit circle, a = (I/J)^(2/3), e = |a - 1|/a (eq. 17.23, p.153), with I, J mutually prime and not both 1, and direction of motion eps' = +1 or -1. It is made of n identical basic arcs (an arc being 2 pi I long in time and J revolutions of the particle). DeltaC is the distance in Jacobi constant from the bifurcation orbit along a family, mu the mass ratio, and **nu is defined by DeltaC = O(mu^nu)** (eq. 11.80), equivalently ln|DeltaC| = nu ln(mu) to leading order. A fixed nu with mu going to 0 selects one asymptotic regime. For type 2 the book finds the sequence (eq. 11.83, p.15) nu = 0, 0 < nu < 1/3, nu = 1/3, 1/3 < nu < 1/2, nu = 1/2, and nothing above (section 17.8, p.179). For type 1 the sequence is nu = 0, (0, 1/2), 1/2, above 1/2 (eq. 11.82; part A digest). Type 3 has more transitions (section 23.4, below).

Basic constants of a type 2 orbit (eqs. 17.24 to 17.27, p.153, from the page image): fixed-axes velocity at a collision V_X = 0, V_Y = eps' sqrt((2a - 1)/a); rotating-axes v_y = V_Y - 1; relative speed v = |V_Y - 1|; Jacobi constant C = 2 eps' sqrt((2a - 1)/a) + 1/a (which equals 3 - v^2). COMPUTED check against volume I: for eps' = +1 and A = a^(3/2) = 2, 1.5, 3, 4, 5, 6 this C is 2.970934, 2.987424, 2.945907, 2.929161, 2.917266, 2.908345, against the volume I Table 4.4 critical-arc values at nearby A of 2.970940/2.970936/2.970935, 2.987426/2.987425, 2.945910/2.945908, 2.929162, 2.917267, 2.908345 (agreement 6e-6 or better; a critical arc is near but not exactly at the tangent ellipse). So the two volumes are consistent, and the type 2 ellipses with eps' = +1 are the A -> integer limit of the Table 4.4 case 5 and 7 arcs (high-C, slow encounters: v = 0.1705 for (I, J) = (2, 1), 0.1121 for (3, 2), 0.0838 for (4, 3), 0.2326 for (3, 1), COMPUTED).

New notations (17.1, p.149): T^f and T^g replace T^i and T^e so results hold for either sign of a - 1: a T^f arc moves out of the region occupied by the bifurcation orbit (relative to the unit circle), T^g into it. Relative side of passage sigma' = sign(a - 1) sign(x_0 - 1) (eq. 17.1). Branches carry a symbolic sign eps' sign(DeltaC) (the sign "+" branches contain only S-arcs; "-" branches contain T-arcs). Sequences are written with integers (S-arc of m basic arcs), f, g (T-arcs), and signs.

## 2. Chapter 17: the fundamental equations and the scaling cascade

### 2.1 The fundamental system (17.2, pp.150-157; READ, equation images for 17.63 to 17.66)

The time of an encounter is defined as the time of conjunction (M3 crossing the x axis), h_i = signed distance from M2 at conjunction, x = 1 + h_i, sign(h_i) = side of passage (17.5, 17.6). Unknowns per basic arc: Delta a_i (semi-major axis change of the intermediate Keplerian arc), radial velocities u'_i (arrival) and u''_i (departure), h_i, plus DeltaC. After expanding the Kepler relations about the bifurcation ellipse and eliminating the eccentricity, anomaly and time variables, the system reduces to four relations per arc (eqs. 17.63 to 17.67; the page 157 image):

- encounter: h_i (u'_{i+1} - u''_i) = -(2/v) mu [1 + o(1)] (the matching relation projected on the line of centres, 17.63, image p.156; the minus sign is the side-of-passage convention, h being signed);
- arc: u''_i - u'_i = (3 pi I (a - 1)/(a^2 v_y)) Delta a_i [1 + ...] (17.64);
- energy: 2 (a - 1) h_{i-1} = a V_Y DeltaC + (v_y/a) Delta a_i + a u'_i^2 + O(DeltaC^2, Delta a_i^2, u'^4, mu) (17.65);
- time: h_i - h_{i-1} = (3 pi I/(2 a v_y)) Delta a_i (u'_i + u''_i) [1 + ...] + O(mu) (17.66);

with h_0 = h_n = O(mu) for a partial bifurcation (17.67). Equations 17.64 to 17.66 were read on the page image (p.157); I did not re-derive them. The case n = 1 (bifurcations 2T1 total and 2P1 partial) has different properties and is treated in chapter 23 (17.2.3).

Hand-algebra consistency (not a program run): (17.63), |h| |u'_{i+1} - u''_i| = 2 mu/v, is the small-angle limit of the project's gate relation tan(delta/2) = mu/(h v^2), with the demanded turn delta of the rotating-frame relative velocity (of speed v, directed along y at conjunction, so a change of radial velocity Delta u is a turn Delta u/v): delta = Delta u/v = 2 mu/(v^2 h). The book's equation is therefore the same physics as the gate and as Perko's and Guillaume's hyperbola; it adds the dependence of h and Delta u on DeltaC.

### 2.2 The case nu = 0 (17.3, pp.158-160) and 0 < nu < 1/3 (17.4, pp.160-166)

At mu = 0 (branches leaving the bifurcation orbit): T-arc has Delta a = 0 and u' = u'' = +/- sqrt(-V_Y DeltaC) (17.76; + for a T^e arc); S-arc of m basic arcs has Delta a = -(a^2 V_Y/v_y) DeltaC (17.84) and u'_{i+j} = (m - 2j + 2) kappa DeltaC, u''_{i+j} = (m - 2j) kappa DeltaC with kappa = 3 pi I (a - 1) V_Y/(2 v_y^2), j = 1..m (17.86), and h_{i+j} = -j (m - j) (9 pi^2 I^2 a (a - 1) V_Y^2/(2 v_y^4)) DeltaC^2 (17.87).

For mu > 0 and 0 < nu < 1/3 the scaled system has a non-vanishing Jacobian (-4 per T-arc, -2m per S-arc of m basic arcs, 17.110 to 17.112, p.163-164), so the mu = 0 branches persist, with error O(DeltaC^(1/2), mu DeltaC^(-3)) (17.113) and the limit nu < 1/3 comes from requiring that error to vanish (17.103). Dominant terms at the encounters (READ, image p.165):

- TT node: h_i = +/- (mu/v) (-V_Y DeltaC)^(-1/2) [1 + O(DeltaC^(1/2))], + for a T^e T^i node, - for T^i T^e (17.125);
- TS or ST node: h_i = +/- (2 mu/v) (-V_Y DeltaC)^(-1/2) (17.126);
- SS node: h_i = -(4 v mu)/(3 pi I (m_a + m_b)(a - 1) V_Y) DeltaC^(-1) [1 + ...] (17.127);
- T-arc: Delta a_i = a (g' + g'') sign(v_y) mu DeltaC^(-1)/(3 pi I V_Y) with g = 2 if the neighbouring arc is an S-arc, 1 if a T-arc, 0 if absent (17.128 to 17.130); g' + g'' = 0 only in 2P1.

Hand-algebra check (not a program run): (17.125) at a T^e T^i node (u''_i = +sqrt(-V_Y DeltaC), u'_{i+1} = -sqrt(-V_Y DeltaC)) gives h (u'_{i+1} - u''_i) = -2 mu/v, and (17.127) with the S-arc radial velocities of 17.123 at an SS node (u'_{i+1} - u''_i = kappa (m_a + m_b) DeltaC, kappa = 3 pi I (a - 1) V_Y/(2 v_y^2), v_y^2 = v^2) gives -2 mu/v again, matching 17.63 including the sign. So the printed coefficients of 17.125 and 17.127 are mutually consistent.

### 2.3 Transition 2.1 at nu = 1/3 (17.5, pp.166-171) and the range 1/3 < nu < 1/2 (17.6, pp.171-176)

At DeltaC = O(mu^(1/3)): S-arc quantities Delta a, u of order mu^(1/3), h of order mu^(2/3); TT, TS, ST node h of order mu^(5/6); SS node h of order mu^(2/3) (17.131 to 17.134). So h inside an S-arc and h at an SS node have fused to the same order: consecutive S-arcs fuse into one R-region (an R-arc if the bifurcation is partial or T-arcs are present, an R-orbit if the whole orbit is S-arcs, so a total transition 2.1 needs a total bifurcation with only S-arcs; Fig. 17.1, p.167). T-arcs are not affected by transition 2.1 (17.145 and 17.158, Delta a of order mu^(2/3), u of order mu^(1/6)). New scaling (17.138 to 17.142, image p.169):

  DeltaC = mu^(1/3) w*, w* = -(v_y^2 L1/(3 pi I V_Y (a - 1))) w, L1 = [4 (a - 1)/(a v)]^(1/3), u' , u'' = -L1 x mu^(1/3), h = mu^(2/3) (2/(v L1)) y,

with the dimensionless unknowns (w, x', x'', y, z) obeying (17.146 and 17.154), for an R-arc of order n-tilde starting at its origin:

  x''_i - x'_i + w = 0; y_i - y_{i-1} + w (x'_i + x''_i) = 0; y_i (x'_{i+1} - x''_i) - 1 = 0, with y_0 = y_n-tilde = 0 (system 19.1, p.199; z_i = -w has been eliminated, 17.153).
  Eliminating x', x'' gives the three-term map y_{i-1} - 2 y_i + y_{i+1} = -2 w/y_i + 2 w^2 (19.3, p.200, image).

The v in L1 (and in L2, L3 below) is printed as a lower-case v; I read it as the relative speed v = |v_y| of 17.26 (the only v defined in the chapter). INFERRED reading, consistent with 17.141g.

For 1/3 < nu < 1/2 (the R-region regime), the scaled system loses w (17.171, 17.181, p.172-174): with x_i the common value of x'_i and x''_i,

  y_i (x_{i+1} - x_i) - 1 = 0 (i in C), y_i - y_{i-1} + x_i = 0 (i in A), y_0 = y_n = 0  (18.1 for an R-arc).

R-region physical quantities (17.190, 17.191): Delta a = -(a^2 V_Y/v_y) DeltaC (as for S-arcs), u', u'' proportional to mu^(1/2) (-eps' DeltaC)^(-1/2) times the numbers x_i of 18.1, and h proportional to mu^(1/2) (-eps' DeltaC)^(1/2) times y_i, with a constant L2 (17.168, p.173) that I did not transcribe. The orders of magnitude are the usable result (section 2.5). The numbers x_i, y_i are the solutions of Tables 18.2 and 18.3.

### 2.4 Transition 2.2 at nu = 1/2 (17.7, pp.176-178) and the absence of nu > 1/2 (17.8, p.179)

At DeltaC = O(mu^(1/2)) everything fuses: T-arcs and R-regions form a single region, with Delta a_i = O(mu^(1/2)), u'_i = u''_i = O(mu^(1/4)), h_i = O(mu^(3/4)) (17.197). Scaling (17.198 to 17.200, image p.177):

  Delta a = mu^(1/2) z*, u' , u'' = mu^(1/4) x*, h = mu^(3/4) y*, DeltaC = mu^(1/2) W*, W* = -(L3^2/V_Y) W, L3 = sign(a - 1)[2v/(3 pi I a)]^(1/4).

The scaled system (17.205, 17.207) is x' = x'' = X_i, X_i^2 - Z_i - W = 0, Y_i - Y_{i-1} - Z_i X_i = 0, Y_i (X_{i+1} - X_i) - 1 = 0, Y_0 = Y_n = 0 (partial, 21.1) or periodic (total, 22.30). It is an area-preserving mapping (21.11) in the (X, Y) plane; Henon's unpublished numerical exploration shows mixed regular and chaotic orbits, hence "cannot be solved explicitly in general".

**There is no regime nu > 1/2 for type 2 (n larger than 1): W never vanishes (21.7, 22.5), and in fact W > 2^(-1/2) for a partial bifurcation (eqs. 21.4 to 21.7, p.240: |X_i| < W^(1/2), |Z_i| < W, |Y_1| < W^(3/2), |Y_1| > 2^(-1) W^(-1/2) gives W > 2^(-1/2)).** So every branch has a minimum of |DeltaC| of order mu^(1/2) along its characteristic, and all branches are joined two by two at nu = 1/2. At nu = 1/2 the third of the three ways to leave the mu = 0 picture ends: branches reconnect.

### 2.5 The scaling table (type 2, n larger than 1; READ Fig. 17.2 p.177 and equations 17.119 to 17.128, 17.157 to 17.160, 17.184 to 17.194, 17.197)

Orders of magnitude in mu and DeltaC. Rows are the encounter type or arc; "S/R antinode" is an encounter inside an S-arc or R-region.

| Quantity | 0 < nu < 1/3 | nu = 1/3 | 1/3 < nu < 1/2 | nu = 1/2 |
|---|---|---|---|---|
| T-arc u', u'' | (-V_Y DeltaC)^(1/2) | mu^(1/6) | (-V_Y DeltaC)^(1/2) | mu^(1/4) |
| T-arc Delta a | mu DeltaC^(-1) | mu^(2/3) | mu DeltaC^(-1) | mu^(1/2) |
| S or R arc Delta a | -(a^2 V_Y/v_y) DeltaC | mu^(1/3) | same | mu^(1/2) |
| S arc u', u'' | proportional to DeltaC (17.123) | mu^(1/3) | R-region: mu^(1/2) DeltaC^(-1/2) | mu^(1/4) |
| h at S or R antinode | DeltaC^2 | mu^(2/3) | mu^(1/2) DeltaC^(1/2) | mu^(3/4) |
| h at SS node (R-node) | mu DeltaC^(-1) | mu^(2/3) | mu DeltaC^(-1/2) | mu^(3/4) |
| h at TT, TS, ST node | mu DeltaC^(-1/2) | mu^(5/6) | mu DeltaC^(-1/2) | mu^(3/4) |

For nu = 1/3 and 1/2 the columns are the order of the quantities at DeltaC = mu^nu.

### 2.6 What the 2T1 and 2P1 chapters add (23.1 to 23.2, pp.283-294)

**Total bifurcation 2T1 (n = 1, one basic arc, symmetric orbits assumed).** There is a single transition at nu = 1/2 (no 1/3 transition; Fig. 23.2, p.289), because no T-arcs or S-arc fusion exist. The asymptotic equations at nu = 1/2 (23.21 to 23.23, image p.287) are y x' + 1 = 0, 2 x' + z = 0, 2 y + 2 w - z = 0 and reduce to **y^2 + w y - 1 = 0**, y = (-w +/- sqrt(w^2 + 4))/2, z = w +/- sqrt(w^2 + 4), x' = (-w -/+ sqrt(w^2 + 4))/2; the Jacobian is -2x' + 2y = +/- 2 sqrt(w^2 + 4), never zero. Physical variables (23.25, image p.287; scaling 23.18, 23.19 on p.286):

  h_0 = a/(4 (a - 1)) [V_Y DeltaC +/- sqrt(V_Y^2 DeltaC^2 + 16 v mu/(3 pi I a))] + O(mu),
  Delta a_1 = (a^2/(2 v_y)) [-V_Y DeltaC +/- sqrt(V_Y^2 DeltaC^2 + 16 v mu/(3 pi I a))] + O(mu),
  u'_1 = -u''_1 = (3 pi I (a - 1)/(4 v^2)) [V_Y DeltaC -/+ sqrt(V_Y^2 DeltaC^2 + 16 v mu/(3 pi I a))] + O(mu).

This is the avoided crossing (hyperbola) between the first-species family (the "E" branches, circle-like) and the second-species family (the "1" branches): the product of the two linear forms is of order mu, and the minimum separation at DeltaC = 0 is sqrt(16 v mu/(3 pi I a)), of order mu^(1/2). Junctions (Table 23.1, p.289): +E joined to -1, and -E to +1 (the first-species branch passes into the second-species branch of the opposite sign). The book states that this agrees with Guillaume's equation (IV-33) of the 1971 thesis (p.288) and with the Broucke-principle result of volume I Table 8.18. For 0 < nu < 1/2 (23.14 to 23.17, read on images pp.285-286): the first-species branch has h_0 = (a V_Y/(2 (a - 1))) DeltaC [1 + O(DeltaC, mu DeltaC^(-2))], u'_1 = -u''_1 = -(2 (a - 1)/(a V_Y v)) mu DeltaC^(-1) [1 + ...] and Delta a_1 = (4 a sign(v_y)/(3 pi I V_Y)) mu DeltaC^(-1) [1 + ...]; the second-species branch has the S-arc forms of 17.122, 17.123, 17.127 with m = 1 (23.16, 23.17). At large |DeltaC| the plus-sign root of 23.25 reduces to the first-species h_0 (COMPUTED: (a/(4 (a - 1))) 2 V_Y DeltaC = a V_Y DeltaC/(2 (a - 1))), a check that 23.25 and 23.14 are mutually consistent. The case nu > 1/2 (23.1.4) reduces to the W = 0 points of Fig. 23.1, which are a subset of the nu = 1/2 solutions, so nothing new.

COMPUTED illustration of the 2T1 minimum impact distance, |h_0| at DeltaC = 0 = (a/(4|a - 1|)) sqrt(16 v mu/(3 pi I a)) at mu = 0.0121529529 (Earth-Moon, the project's value): (I, J, eps') = (2, 1, +1) 0.0225; (3, 2, +1) 0.0256; (4, 3, +1) 0.0271; (3, 1, +1) 0.0134; (1, 2, +1) 0.0461; (2, 1, -1) 0.0802; (3, 2, -1) 0.1111; (1, 2, -1) 0.0987 (units of the Earth-Moon distance; multiply by 384 400 km for kilometres: 8.6e3 km, 9.8e3, 1.04e4, 5.2e3, 1.77e4 km for the first five). These are first-order numbers at a mass ratio where the approximation is marginal (see section 7); they are scale estimates, not predictions.

**Partial bifurcation 2P1 (n = 1).** The only case where the O(mu) intermediate-orbit approximation fails (11.3.3 note, and the use of the reduced three-equation system 23.26). Results: T-arcs (branches -f, -g, the asymmetric pair): the asymptotic solution exists for 0 < nu < 2/3 (23.35) with u'_1, u''_1 = +/- sqrt(-V_Y DeltaC) [1 + ...], Delta a_1 = O(mu DeltaC^(-1/2)) (23.40; "the only case where g' and g'' both vanish"); at nu = 2/3 the O(mu) term becomes the same order as the error and the approach is indeterminate: "it would be necessary to go to a higher-order approximation ... in the range nu >= 2/3" (p.292). S-arcs (branches +1, -1, the symmetric pair, assuming symmetry persists): Delta a_1 = -(a^2 V_Y/v_y) DeltaC [1 + ...], u'_1 = (3 pi I (a - 1) V_Y/(2 v_y^2)) DeltaC [1 + ...] for 0 < nu < 1 (23.45 to 23.51); the approach fails at nu = 1 (p.294). **So the quantitative method does not describe the 4 branches of 2P1 near the bifurcation, and cannot establish their junctions; the book reproduces the symmetry-derived junctions (Table 23.2, p.294): 2P1S: +1 joined to -1; 2P1A: -f joined to -g.** This is the one place in chapters 17 to 23 where the answer to "what happens at distances of order mu or mu^(2/3) from the bifurcation orbit" is "not determined at this order".

## 3. Chapter 18 (pp.181-197): the R-region equations, solved

The R-region system (18.1, 18.16) for an R-arc (y_0 = y_n = 0) or R-orbit (periodic) eliminates x to the three-term recurrence **y_{i-1} - 2 y_i + y_{i+1} + 1/y_i = 0** (18.4, 18.17), an area-preserving map F: (alpha_i, beta_i) = (y_i, y_{i+1}) -> (beta_i, -alpha_i + 2 beta_i - 1/beta_i) (18.61, 18.62), and the central results are:

- **R-arc of order n-tilde**: exactly 2^(n-tilde - 1) solutions, exactly one for each choice of the signs of y_1 .. y_(n-tilde - 1) (Propositions 18.1.1 to 18.1.3; y_i never vanishes). The sign of y_i is the relative side of passage (17.169); an R-arc therefore labels one of the 2^(n-tilde - 1) decompositions of the arc into S-arcs (a "+" is a node between S-arcs, a "-" an antinode inside one). R-arcs are never critical (dy_n/dx_1 never vanishes, 18.1.3).
- **R-orbit of order n**: exactly 2^n - 2 solutions, one for each sign sequence of (y_0, .., y_(n-1)) except all + and all - (which give infinite y); sub-period solutions are included in the count (spurious solutions with sub-period 2, 3 counted separately in 18.2.3.3 and 18.2.3.5). Sum of 1/y_i over the orbit is zero (18.84). R-orbits are always unstable, z larger than 1 (Proposition 18.2.2, proof by the bound u_i larger than 2 in 18.29 to 18.43), and the Jacobian vanishes only at critical orbits of the first kind (Proposition 18.2.1).
- **Equivalence to the baker transformation** (18.3, pp.192-197): the map F is topologically conjugate to the baker transformation (Proposition 18.3.1, from Devaney 1981; a continuous increasing function v = 0.w_1 w_2 ... in binary of the sign sequence, with the "devil's staircase" of Fig. 18.1). Every R-arc and R-orbit is therefore highly unstable (deviation doubles per step); this is the origin of the factor-2 sensitivity and of the instruction to print 8 digits in Table 18.2 "because the R-arcs are strongly unstable: errors are amplified during the computation of the successive y_i".

Closed forms (READ pp.184, 188-191, images): n-tilde = 2: y_1 = 1/sqrt(2) (18.11); n-tilde = 3: y_1 = y_2 = 1 and y_1 = -y_2 = 1/sqrt(3) (18.12, 18.13); n-tilde = 4: closed forms 18.14 and 18.15 exist (not transcribed; see the Table 18.2 values); R-orbits n = 2: y_0 = -y_1 = 1/2; n = 3: y_0 = 1/sqrt(6), y_1 = y_2 = -2/sqrt(6); n = 4 (18.46, 18.47): y_0 = y_1 = 1/sqrt(2), y_2 = y_3 = -1/sqrt(2) and y_0 = (1 + sqrt(3))/2, y_1 = y_3 = 1, y_2 = (1 - sqrt(3))/2 with sign reversals; n = 5 symmetric: 5 y_1^6 - 20 y_1^4 + 17 y_1^2 - 4 = 0 (18.49); n = 6: the closed forms 18.53 to 18.57 (Table 18.3 values). For n = 7 the first sign sequence with no symmetry appears (+++-+--), computed numerically: y_0 = 0.880142094, y_1 = 1.302150709, y_2 = 0.956199, y_3 = -0.435560, y_4 = 0.468576, y_5 = -0.761411, y_6 = -0.678047 (p.191; printed to 9 digits for y_0 and y_1 and 6 digits for the others, a precision difference in the book itself; COMPUTED residual not run).

### Table 18.2 (p.185): R-arcs, y_1 for every sequence beginning R+

(READ on the page image; COMPUTED: I solved the recurrence 18.4 independently by multi-start Newton for every order 2 to 6; all 31 printed y_1 values, with their sign sequences, are reproduced to the printed 8 digits and the solver finds exactly 2^(n-tilde - 1) solutions per order, so the table is complete as printed. The sign-sequence labels were machine-compared for the 28 rows whose text-layer line parsed cleanly and compared by eye for the other three; all 31 y_1 values match. The sequence "R+R-R+R-R" has y_1 = 0.59586158 and "R+R-R-R+R" has 0.54119610.)

| n-tilde | sequence | y_1 |
|---|---|---|
| 2 | R+R | 0.70710678 |
| 3 | R+R+R / R+R-R | 1.00000000 / 0.57735027 |
| 4 | R+R+R+R / R+R+R-R / R+R-R+R / R+R-R-R | 1.17914724 / 0.94010422 / 0.59967641 / 0.53185593 |
| 5 | R+R+R+R+R / R+R+R+R-R / R+R+R-R+R / R+R+R-R-R | 1.30656296 / 1.14267956 / 0.95043145 / 0.91921106 |
| 5 | R+R-R+R+R / R+R-R+R-R / R+R-R-R+R / R+R-R-R-R | 0.60746124 / 0.59586158 / 0.54119610 / 0.50526000 |
| 6 | R+R+R+R+R+R / R+R+R+R+R-R / R+R+R+R-R+R / R+R+R+R-R-R | 1.40472828 / 1.28123317 / 1.14898685 / 1.12992127 |
| 6 | R+R+R-R+R+R / R+R+R-R+R-R / R+R+R-R-R+R / R+R+R-R-R-R | 0.95398640 / 0.94868391 / 0.92357558 / 0.90676764 |
| 6 | R+R-R+R+R+R / R+R-R+R+R-R / R+R-R+R-R+R / R+R-R+R-R-R | 0.61196435 / 0.60587822 / 0.59651693 / 0.59452140 |
| 6 | R+R-R-R+R+R / R+R-R-R+R-R / R+R-R-R-R+R / R+R-R-R-R-R | 0.54444346 / 0.53958978 / 0.51070145 / 0.48686874 |

(Sequences beginning R- are the same values with all signs reversed, property E', 18.3.)

### Table 18.3 (p.189): R-orbits, y_0 .. y_(n-1) (one representative per orbit, shifts of the origin omitted; "sign reversal" rows are the E' images)

(READ on the image; COMPUTED: every printed row satisfies y_(i-1) - 2 y_i + y_(i+1) + 1/y_i = 0 cyclically to a maximum residual of 4e-6, which is the rounding of six printed decimals.)

| n | y_0 | y_1 | y_2 | y_3 | y_4 | y_5 |
|---|---|---|---|---|---|---|
| 2 | 0.500000 | -0.500000 | | | | |
| 3 | 0.408248 | -0.816497 | -0.816497 | | | |
| 3 | -0.408248 | 0.816497 | 0.816497 | | | |
| 4 | 0.707107 | 0.707107 | -0.707107 | -0.707107 | | |
| 4 | 1.366025 | 1.000000 | -0.366025 | 1.000000 | | |
| 4 | -1.366025 | -1.000000 | 0.366025 | -1.000000 | | |
| 5 | 1.712938 | 1.712938 | 1.129146 | -0.340271 | 1.129146 | |
| 5 | -1.712938 | -1.712938 | -1.129146 | 0.340271 | -1.129146 | |
| 5 | 0.799673 | 0.799673 | -0.450837 | 0.516750 | -0.450837 | |
| 5 | -0.799673 | -0.799673 | 0.450837 | -0.516750 | 0.450837 | |
| 5 | 0.652966 | 0.652966 | -0.878507 | -1.271686 | -0.878507 | |
| 5 | -0.652966 | -0.652966 | 0.878507 | 1.271686 | 0.878507 | |
| 6 | 1.224745 | 0.816497 | -0.816497 | -1.224745 | -0.816497 | 0.816497 |
| 6 | 1.345162 | 0.973460 | -0.425507 | 0.525667 | -0.425507 | 0.973460 |
| 6 | -1.345162 | -0.973460 | 0.425507 | -0.525667 | 0.425507 | -0.973460 |
| 6 | 2.193233 | 1.965259 | 1.228446 | -0.322404 | 1.228446 | 1.965259 |
| 6 | -2.193233 | -1.965259 | -1.228446 | 0.322404 | -1.228446 | -1.965259 |
| 6 | 1.618034 | 1.618034 | 1.000000 | -0.618034 | -0.618034 | 1.000000 |
| 6 | -1.618034 | -1.618034 | -1.000000 | 0.618034 | 0.618034 | -1.000000 |
| 6 | -0.720707 | 0.720707 | 0.774597 | -0.462509 | 0.462509 | -0.774597 |
| 6 | 0.720707 | -0.720707 | -0.774597 | 0.462509 | -0.462509 | 0.774597 |

Counts stated by the book and checkable: the number of solutions of (18.18) is at most 2^n - 2 for an R-orbit; for n = 5 there are 30 symmetric solutions (6 listed times 5 shifts) and the equation 18.49 has all of them (so all solutions are found); for n = 6 the 9 listed families give 54 solutions with shifts, plus 2 spurious sub-period-2 and 6 spurious sub-period-3, total 62 = 2^6 - 2. The closed forms 18.53 to 18.57 for n = 6 (the symmetric and the asymmetric 321 orbits, with 5, 2 and 2 solutions in three groups) are not transcribed; use the Table 18.3 values.

## 4. Chapters 19 and 20: partial and total transition 2.1 (nu = 1/3)

These chapters solve the nu = 1/3 system (19.1 for an R-arc, 20.1 for an R-orbit), whose parameter w replaces DeltaC (DeltaC = mu^(1/3) w*, section 2.3). Results (READ):

- The system contains no parameter: the equations are identical for all transitions 2.1, independently of I, J, eps' (19.1 property 1; contrast type 1 where a constant K enters).
- w never vanishes for n-tilde at least 2 (19.1 property 3: w = 0 forces all y_i equal, contradicting y_0 = y_n = 0) and has the sign of -eps' DeltaC.
- Asymptotic branches: for w -> +/- infinity only the decompositions of the R-arc into S-arcs (19.8: y_{i+j} = -j (m - j) w^2 inside an S-arc of m basic arcs, y_i = 2/(m_a + m_b) w^(-1) at an SS node, x'_{i+j} = -(m - 2j + 2) w/2, x''_{i+j} = -(m - 2j) w/2); for w -> 0 only R-arcs (19.34); **no branch with w -> 0 from below exists for n-tilde larger than 1** (19.1.3, proof 19.46 to 19.47), so the plus-sign branches (w negative, S-arcs only) must join each other and the minus-sign branches (w positive) pass through to nu = 1/2.
- Stability: an R-Jacobian vanishes iff the R-arc is critical, dy_n/dx_1 = 0 (Proposition 19.1.1); the Jacobian of the whole bifurcation vanishes iff it contains a critical R-arc (19.1.2). The R-arc variations are amplified by 2 w m per node (19.1.2).

Small n-tilde: the characteristic in the (w, y_1) plane; closed forms (READ images p.209-211, 229-232, 233-235):

| Case | Equation (as printed) | Extremum in w (critical orbit) | Source |
|---|---|---|---|
| 2.1P2 (n-tilde = 2) | y_1^2 + w^2 y_1 - w = 0 (19.57) | w^3 = -4, w = -2^(2/3) = -1.587401 (extremum) | 19.62, p.209 |
| 2.1P3 | y_1 = y_2 = (-w^2 +/- sqrt(w^4 + 2 w))/2... (19.65, text layer garbled in the factor) | first solution w^3 = -2 (extremum), w^3 = -8/3 (extremum and intersection); second solution w^3 = -8/3 (intersection) | 19.67, 19.68, p.211 |
| 2.1T2 | y_0 y_1 = -w/2, y_0 + y_1 = -w^2 (derived on p.229 from 19.3 with the spurious y_0 = y_1 removed) | w^3 = -2, where z = 1 (the book's z formula, a cubic in w, is garbled in the text layer and not used) | p.229 |
| 2.1T3 | symmetric pair of cases (20.30 to 20.32) | (3 w^3 + 4)(2 w^3 + 3)^2 = 0: w^3 = -4/3 extremum; w^3 = -3/2 is an intersection with a 2T1 orbit described three times (z = (-49 -/+ 81)/32 = 1 for the lower sign) | 20.36 to 20.38, p.231-232 |
| 2.1T4, solutions h = 0 and h = 2 | z - 1 = product of three factors (20.42, image p.233) | first factor p^2 + 4 w vanishes for omega_2 = +/-1, w^3 = -1 and for omega_2 = +1, w^3 = -9/4; second for omega_2 = +1, w^3 = -1; third for w^3 = -8/9 (both signs); all extremums | p.233 |

(The exact positions of the minus signs in the omega statement of 2.1T4 are from the image of p.233; the printed first-factor statement reads "omega_2 = +/-1, w^3 = -1 and for omega_2 = +1, w^3 = -9/4".)

**A disagreement with Guillaume (1971) is recorded** (p.233): for the 2T4 subset with origin at a collision, Guillaume's results make +4 join +121 and +-1111- join +-22-; Henon's results make +4 join +-22- and +-1111- join +121. Henon's explanation: Guillaume did not notice that p^2 + 4 w becomes negative for omega_2 = +1 in the interval -9/4 < w^3 < -1. COMPUTED: I did not re-solve this case.

**Junction rules at transition 2.1** (READ 19.3.3, 19.4, 20.2.5, 20.3): (i) every branch with a minus sign (w positive) passes through transition 2.1 unchanged and goes on to nu = 1/2 (Propositions 19.4.1, 19.4.2); this case needs only Broucke's principle (relative sides of passage). (ii) For n at least 2 all branches with a plus sign (w negative, S-arcs only) are joined among themselves at transition 2.1 (Proposition 19.4.4), so for them nu rises from 0 to a maximum of 1/3 and falls back to 0 without reaching higher values. (iii) n = 1 plus branch passes through (Proposition 19.4.3, the 2P1 case, which later fails, section 2.6). The "positional method" (19.3) finds the junctions by the relative order of the asymptotic branches in the (w, y_1) plane (no two characteristics for different n-tilde intersect), and works for all partial cases to n = 6, but only needs the plus-sign cases after the simpler minus-sign rule; for total transition the plus-sign junctions for n above 4 are by numerical computation (20.2.5.2).

### Table 19.1 (p.224): partial bifurcation of type 2, n = 2 to 6, junctions of the plus branches

(READ on the image, each line pair is one joined pair of branches, written as the S-arc decomposition with a plus sign; the header "2P5----S" is the book's group label.) COMPUTED check: for each n from 2 to 6 the listed pairs partition all 2^(n-1) compositions of n exactly once, each pair has the same sum n, as Proposition 19.4.4 demands; this is a transcription and consistency check, not an independent derivation of which compositions pair.

| n | pairs (+ a, + b) |
|---|---|
| 2 | (2, 11) |
| 3 | (3, 111), (21, 12) |
| 4 | (22, 1111), (4, 121), (211, 112), (31, 13) |
| 5 | (131, 5), (11111, 212), (14, 41), (122, 1211), (32, 311), (1112, 2111), (1121, 221), (113, 23) |
| 6 | (141, 6), (1221, 33), (111111, 2112), (11211, 222), (15, 51), (132, 1311), (42, 411), (12111, 1212), (3111, 312), (123, 321), (1113, 213), (11121, 2121), (11112, 21111), (1122, 2211), (1131, 231), (114, 24) |

### Table 20.1 (p.238): total bifurcation of type 2, n = 2 to 6, junctions of the plus branches (n up to 4 transcribed)

(READ on the image; the header "2T4-----1" etc. labels the origin shift; minus signs around a sequence mark which positions are collisions.)

| Group | pair |
|---|---|
| 2T2---0 | (+2, +-11-) |
| 2T3----0 | (+3, +-111-) |
| 2T3----1 | (+-21-, +12) |
| 2T3----2 | (+-12-, +21) |
| 2T4-----0 | (+4, +-22-), (+-1111-, +121) |
| 2T4-----1 | (+-211-, +13) |
| 2T4-----2 | (+-31-, +-13-), (+211, +112) |
| 2T4-----3 | (+31, +-112-) |

(n = 5 and 6 are on the same page; not transcribed.)

## 5. Chapters 21 and 22: partial and total transition 2.2 (nu = 1/2)

Chapter 21 solves (21.1) for a partial bifurcation, chapter 22 solves (22.1) for a total one. Results (READ):

- **Only minus-sign branches (W positive) reach transition 2.2** (21.1 property 3: W > 0 for n at least 2; plus-sign branches closed at 2.1). W has the sign of -eps' DeltaC.
- W > 2^(-1/2) for all partial bifurcations, n at least 2 (21.7).
- Asymptotic branches for W -> +infinity: sequences of T-arcs and R-arcs only (21.12 to 21.19); T-arc X_i = +/- W^(1/2) (T^f +, T^g -, 21.14), Z_i = g W^(-1) /2 (21.12); R-arc X_i = W^(-1/2) x_i, Y_i = W^(1/2) y_i, Z_i = -W (21.15), with x_i, y_i the numbers of Tables 18.2 and 18.3; nodes Y_i = +/- W^(-1/2) (TR or RT, 21.18; TT with opposite signs 21.17). The branches are denoted by sequences with R (an R-arc from 2.1), f, g, and signs (e.g. -R-R-R+R).
- **Every junction at nu = 1/2 is a pairing of two branches that then continue back to small nu, so the whole path of a minus branch is (Table 21.2, p.265): nu from 0 up to 1/3 (S and T arcs), through transition 2.1 into R and T arcs, to nu = 1/2 where it joins another branch, then back down through transition 2.1 (R-arcs split into S-arcs again) and on to nu = 0.** Example, the book's: -f31 -> (at 2.1) -f-R-R-R+R -> (at 2.2, Fig. 21.8) joined to -f-R-R-g+R -> (back through 2.1) -f2g1.
- Extremum and intersection values of W (the minimum of W along a characteristic is the closest approach in DeltaC to the bifurcation orbit):

| Case | Equation | Extremum | Intersection | Source |
|---|---|---|---|---|
| 2.2P2 first solution | X_1 = -X_2 = sqrt((W -/+ sqrt(W^2 - 2))/2)-type, exists for W at least sqrt(2) | W = sqrt(2) | W = 3/2 with the other family | 21.67 to 21.70, p.249-250 |
| 2.2P2 second | Z_1^4 + W Z_1^3 + Z_1^2 + W Z_1 + 1 = 0 (21.71) | double root at Z_1 = -1, W = 3/2 (extremum and intersection) | | p.250 |
| 2.2P3 symmetric | X_1 = -X_3, X_2 = 0, Y_1 = Y_2 (21.78, 21.79) | W = 2 | an intersection with another family at a W given in 21.83 (expression not transcribed) | 21.82, p.252 |
| 2.2T2 first solution (22.15, W at least 2) | Y_0 = -Y_1, ... | W = 2 | W = 3/sqrt(2) | 22.21, p.274 |
| 2.2T3 symmetric, h = 0 | closed form below | minimum W = sqrt(6) at Y_0 = 24^(-1/4) | Y_0 = ((3 + sqrt(33))/144)^(1/4) with a family of asymmetric orbits | 22.27 to 22.29, p.276 |

2.2T3 symmetric closed form (READ image p.276, equation 22.27): Y_0 = sqrt((W +/- sqrt(W^2 - 6))/12), X_1 = -X_3 = sqrt((W -/+ sqrt(W^2 - 6))/2), X_2 = 0, Y_1 = Y_2 = -sqrt((W +/- sqrt(W^2 - 6))/3), Z_1 = Z_3 = (-W -/+ sqrt(W^2 - 6))/2, Z_2 = -W. **COMPUTED: this satisfies all six equations (22.30) (Y_i (X_(i+1) - X_i) = 1 and Y_i - Y_(i-1) = X_i (X_i^2 - W), cyclically) to 3e-15 at W = 2.6, 3.0, 5.0 for both signs**, and the stability zero Y_0^2 = 24^(-1/2) at W = sqrt(6) is consistent (Y_0^2 = sqrt(6)/12 = 0.2041 = 24^(-1/2) = 0.2041). The stability formula z - 1 = (27/2)(1/(24 Y_0^4))(1/(288 Y_0^8) + 1/(8 Y_0^4) - 3) (22.28) was not re-derived.

Numerical computation (22.2.4, p.276-278): for n larger than 2 (total) and larger than 3 (partial) the branches are followed from large W toward smaller W by a relaxation method with the T/R-arc decomposition as initial approximation (Y_i = 0 in nodes), switching to shooting when the decomposition fails; the three roots of the cubic X_i^3 - W X_i - (Y_i - Y_(i-1)) = 0 near 0, +W^(1/2), -W^(1/2) are used (T^f root near +W^(1/2), T^g near -W^(1/2)). All branches were computed and joined up to n = 5 by inspection of print-outs (22.3); the results are in Tables 21.3 and 22.2. In the partial case the positional method plus a minimal-W argument (the figure 21.9 contradiction, p.264) settles n = 6 for all but two groups of 4 branches; n = 7 leaves 40 branches (of 972) undecided. Regularities observed (p.265): if a junction is not in a trident, one symbol changes between the two joined branches and it is not an end symbol; in a trident two symbols change, at symmetric positions (arcs i and n + 1 - i).

### Table 21.3 (p.267) and Table 22.2 (p.279): junctions of the minus branches (selection, READ)

Headers are of the form 2P<n><signs><S or A> in the book's volume I Table 8.12 format; I did not decode the meaning of the + and - signs inside the header (volume I, not read). Entries are pairs of joined branches, each written as a sequence of integers, f, g and the sign "-".

Table 21.3, n = 2 (four groups as printed): 2P2+S: (-11, -gf); 2P2+A: (-1f, -g1); 2P2-S: (-2, -fg); 2P2-A: (-1g, -f1). Table 22.2 (total), n = 2: 2T2+-+0: (-2, -fg); 2T2+-+A: (-f1, -1g); 2T2-+-0: (--11-, -gf); 2T2-+-A: (-1f, -g1). Tables 21.3 and 22.2 run to n = 5 (Table 21.3 pp.267-269, Table 22.2 pp.279-281); the n = 3 and 4 columns were read but not transcribed here because I could not attach each header unambiguously to its pairs from the image at the resolution I used (flagged: re-read the page images before using them).

## 6. Chapter 23 and the end matter

- **Type 3** (23.4, p.295): not completed. Transitions identified at nu = 1/5, 2/9, 1/4, (n - 1)/(4 n - 5), 1/3, (n - 1)/(2 n - 1), 1/2, 2/3 (eight values, eq. 23.53, read on the page image of p.295; an earlier draft of this note had a spurious 2/7 from reading text-layer commas as digits) and "quite possibly others remain". At least 12 new kinds of arcs would be needed. Henon writes "it seems now unlikely that I will be able to fulfill this program" (a third volume) and offers "about 400 pages of manuscript notes, in french, on type 3" to any interested colleague. So type 3 (the retrograde unit circle, C = -1, the Guillaume 1973 bifurcation at C = -1) has no quantitative treatment.
- **Newton-polyhedron approach for type 2** (23.3.1): applicable in principle (4n + 1 variables, 4n - 1 equations, solutions on 2-dimensional manifolds, simplex polyhedra), not carried out.
- **Proving general results** (23.3.2): not attempted for type 2.
- **Index of Notations** (pp.299-300): L1, L2, L3 are the auxiliary constants of 17.5, 17.6, 17.7; n-tilde the number of basic arcs in an R-region; R an R-region; Tf, Tg the new T-arc names; sigma' the relative side of passage; V_X, V_Y fixed-axes velocity at collision, v_x, v_y rotating; u', u'' radial velocities; g', g'', g the switching variable of 12.4.2 and 17.4. B is a subscript for Bruno's variables (Newton approach, chapter 15).
- **References list** (pp.301-302): Bruno 1976 (Moscow preprint 95, = Bruno 1994 chapter 7), Bruno 1994 (de Gruyter, Plane periodic orbits), Bruno 1998/2000 (Power Geometry in Algebraic and Differential Equations), Devaney 1981 (Commun. Math. Phys. 80:465-476, the baker transformation), Graham, Knuth, Patashnik 1989 (Concrete Mathematics, for O notation), Guillaume 1971 (Liege thesis), 1973b (Celest. Mech. 8:199-206), 1975b (Celest. Mech. 11:449-467), Henon 1997 (volume I), Henon & Guyot 1970 (Periodic Orbits, Stability and Resonances, Reidel, 349-374), Hitzl & Henon 1977b (Acta Astronautica 4:1019-1039), Mihalas & Routly 1968 (Galactic Astronomy), Perko 1965 (thesis), Perko 1976b (Celest. Mech. 14:395-427), 1977a (Celest. Mech. 16:275-290), 1981a (SIAM J. Appl. Math. 41:181-202), 1981b (Celest. Mech. 24:155-171). The reference list contains Hitzl & Henon 1977b only (the Acta Astronautica paper), not the Celest. Mech. 1977a paper that this repository holds.

## 7. Consistency with Perko 1977 and 1981, Guillaume, and the Perko "nu above 1/2" question

Different nu, same structure. Perko's nu (1977, 1981) is the exponent of the initial-condition variations delta r, delta v = O(mu^nu) about a type 1 bifurcation ellipse (Perko 1981 p.156: "Type 1: ... intersects the Moon's position on the unit circle in two distinct points. Type 2: ... tangent ..."), so that the near-Moon distance is Delta = O(mu^nu) and the deflection beta = O(mu^(1 - nu)). Henon's nu is the exponent of DeltaC = O(mu^nu) (11.80). These are related (a variation of the initial velocity by delta v moves C linearly), but the Perko papers treat only types AA(1,1) and AA(0,1), which are Henon's type 1 (and, in Perko's Part B, a first-species to second-species type 1 case; the existence theory for first-species/second-species types 2 and 3 is Perko 1981a, SIAM, and Guillaume 1971, 1973). **Perko 1977 and 1981 contain no counterpart of Henon's type 2 results (transition 2.1 at nu = 1/3, 2.2 at nu = 1/2, the fusion of S-arcs into R-regions, the baker-map structure of the R-regions, W larger than 2^(-1/2)).** My chapters mention Perko only in the opening of volume II (11.1: bifurcations "studied quantitatively by Guillaume (1971) ... and by Perko (1977a, 1981a, 1981b), again in the case of symmetric orbits"), in the acknowledgements (Perko read a draft) and the reference list; **they state no disagreement with Perko anywhere in chapters 17 to 23.** The one recorded disagreement is with Guillaume 1971 (the 2T4 junction, above).

Do the exponents and regimes agree?
1. nu = 1/2 is the hyperbola for type 1 in both: Perko 1981 eq. (9), (a21 dr_1 + a24 dv_1)(A21 dr_1 + A24 dv_1) = 2 mu a2+4/V_1, and the Henon 2T1 equation y^2 + w y - 1 = 0 (23.22) are the same kind of object (a hyperbola with a mu^(1/2) neck). Guillaume (1971), cited for 2T1 by Henon (p.288, IV-33), agrees with Henon's equations. The Henon 2T1 minimum separation sqrt(16 v mu/(3 pi I a)) is the type 2 analogue of Perko's 2 mu a2+4/V_1.
2. Henon's range 0 < nu < 1/2 for type 1 (eq. 11.82) matches Perko 1981 case (i), and Perko 1977's lower limit nu = 1/3 (below which the problem is a regular perturbation, Part C of the Perko digest) is not a type 1 transition in Henon (type 1 has none at 1/3); for type 2 Henon finds a real transition at exactly nu = 1/3, but for a different reason (fusion of consecutive S-arcs into R-regions, because h at antinodes DeltaC^2 and h at SS nodes mu DeltaC^(-1) reach the same order at DeltaC = mu^(1/3)). I would treat the shared value 1/3 as a coincidence of exponents unless a source says otherwise (INFERRED).
3. **Perko 1981's "Henon privately told me the 1977 scaling above nu = 1/2 was wrong" (Perko digest Part B, section 3: for 1/2 < nu < 1, delta v_1 = O(mu^(1 - nu)), not mu^nu).** The two nu's must be kept apart: Perko's nu_P is the exponent of the near-Moon miss distance and of the initial-condition variation (Delta = h = O(mu^(nu_P)), deflection beta = O(mu^(1 - nu_P))); Henon's nu_H is the exponent of DeltaC. The book does not restate Perko's law as such in chapters 17 to 23, and for type 2 there is nothing above nu_H = 1/2 (17.8: W is bounded below), but **the correct scaling is contained in the book's matching relation (17.63), |h| |u'_{i+1} - u''_i| = 2 mu/v: with h = O(mu^(nu_P)) the velocity change, hence the deflection, is O(mu^(1 - nu_P)) in every regime, which is Perko 1981's corrected deflection law and is the same in the 1977 paper.** What Perko's correction concerns is the variation delta v_1. On Perko's hyperbola the two linear forms are of orders mu^(nu_P) and mu^(1 - nu_P) (their product is of order mu); DeltaC is linear in the initial-velocity variation and so tracks the larger one, giving nu_H = min(nu_P, 1 - nu_P), at most 1/2 (INFERRED from the two sources, not stated in either). So Perko's range nu_P > 1/2 maps to Henon's nu_H = 1 - nu_P < 1/2 on the other asymptote of the hyperbola, not to nu_H > 1/2, where (for type 1) the book says only that W -> 0 points occur at isolated crossings of the nu = 1/2 characteristics with the axis W = 0, with Y_i and Delta a_i still of order mu^(1/2) (12.6 and 12.116, read on p.37 only; the other agent's chapters 12 to 14 should confirm). The mapping between the two exponents depends on the arc type, so Perko's type 1 relation does not transfer to type 2: at an SS node h ~ mu/DeltaC gives nu_H = 1 - nu_P; at a T node h ~ mu DeltaC^(-1/2) (17.125) gives nu_H = 2 (1 - nu_P); at an S antinode h ~ DeltaC^2 gives nu_H = nu_P/2. Conclusion: the book neither states nor contradicts Perko's mu^(1 - nu) law; it implies it through (17.63) for type 1 and replaces it by arc-type-specific laws for type 2.
4. For type 2 the book gives more than Perko's theory for the T-arc and S-arc exponents: the demanded turn at a T node scales as (DeltaC)^(1/2), at an SS node as DeltaC, and saturates at mu^(1/4) and mu^(1/3) respectively (section 8), which Perko's AA(1,1) analysis, being an O(mu^nu) variation of one near-Moon passage, does not give.

## 8. How the results quantify the small-turn and near-resonance cases that the #906 gate must treat as indeterminate

The earlier amendment (OUTSTANDING `#906`, from the Hitzl & Henon digest) says a near-zero demanded turn must be returned as "indeterminate: near resonance or outside first-order validity", with mu |ln mu|/v^3 and r_p/sqrt(mu) as validity numbers. This book gives the structure behind that, for orbits near a type 2 (tangent) bifurcation ellipse (I, J, eps'), and a candidate scale for the indeterminate region. All of the following is derived from the printed equations; the closed-form turn laws are INFERRED applications (turn = Delta u/v, section 2.1) and the numbers COMPUTED.

First-order demanded turn at an encounter, in the mu = 0 limit DeltaC >> mu^(1/2) (turn from the encounter relation |h| |Delta u| = 2 mu/v, turn = |Delta u|/v):
- T node (TT): turn = 2 sqrt(-V_Y DeltaC)/v (u' = -u'', 17.119, 17.121), independent of mu; TS or ST node: sqrt(-V_Y DeltaC)/v. Valid while mu DeltaC^(-2) << 1 (error term of 17.125).
- SS node: turn = (m_a + m_b) (3 pi I |a - 1| |V_Y| |DeltaC|)/(2 v^3), proportional to DeltaC (from 17.123, 17.127), valid for DeltaC >> mu^(1/3).
- S-arc antinode: radial velocity is continuous at leading order (no encounter at all in the mu = 0 limit, 17.69b), and the turn that does appear is delta = 2 mu/(v^2 |h|) with h = -j (m - j) (9 pi^2 I^2 a (a - 1) V_Y^2/(2 v_y^4)) DeltaC^2 (17.124), i.e. a turn of order mu/DeltaC^2.
- The turns vanish as DeltaC -> 0 (like DeltaC^(1/2) at T nodes, DeltaC at SS nodes): the demanded turn at the bifurcation orbit itself is zero, which is why a gate that sets "demanded turn = 0 means no encounter" is consistent at mu = 0 (section 4 of the volume I digest) and why a small turn alone is not a rejection reason.

Saturation: the floor on the demanded turn near the bifurcation orbit comes from the transitions. At DeltaC of order mu^(1/3) the S-arc and SS node turns are of order mu^(1/3) (u of order mu^(1/3), 17.131); at DeltaC of order mu^(1/2) all encounters have radial-velocity change of order mu^(1/4), turn of order mu^(1/4)/v. So below DeltaC of order mu^(1/2) the demanded turn no longer depends on DeltaC at all (it is set by mu and by which branch is followed), and, for n at least 2 (several basic arcs), there is a minimum of |DeltaC| of mu^(1/2) times a coefficient on every minus-sign characteristic (W larger than 2^(-1/2) partial; W = 2 at the 2.2T2 extremum, W = sqrt(2) at 2.2P2, W = sqrt(6) at 2.2T3, W = 2 at 2.2P3). **For n at least 2: an orbit with |DeltaC| below about k mu^(1/2) W_min, where k = L3^2/|V_Y| (section 2.4), is not a perturbation of one mu = 0 arc sequence: it lies on a junction of two branches, and any "demanded turn" computed from the nearest mu = 0 generating orbit is not meaningful.** Returning "indeterminate" there is the right output; returning a number is wrong.

**For n = 1 (a single flyby per period, the usual case) the statement is different and the book is weaker.** 2T1: there is no minimum of |DeltaC|; on each branch of the hyperbola y^2 + w y - 1 = 0 the variable w passes through 0 (23.1.4), so DeltaC passes through zero smoothly, and the neck is a separation of the first-species and second-species branches in h_0 (and Delta a_1, u_1), of order mu^(1/2), not in DeltaC. 2P1: the first-order description fails outright at nu_eff = 2/3 for the T-arc (asymmetric) branches and at nu_eff = 1 for the S-arc (symmetric) branches (23.2.3, 23.2.5, pp.292, 294: "it would be necessary to go to a higher-order approximation"), i.e. for |DeltaC| of order mu^(2/3) or smaller (T) and of order mu or smaller (S), no coefficients given; the junctions there are known only from symmetry (Table 23.2). These 2P1 limits are the book's only explicit "indeterminate" boundary for a single-flyby chain near a tangent ellipse; at mu = 0.0121529529, mu^(2/3) = 0.0529 and mu = 0.0122.

COMPUTED scales at mu = 0.0121529529 (the project's Earth-Moon value; mu^(1/3) = 0.2299, mu^(1/2) = 0.1102, mu^(1/4) = 0.332, mu^(2/3) = 0.0529), for the tangent ellipses with eps' = +1 (direct, slow, v small) and a few retrograde ones. k12 = L3^2/|V_Y| so that |DeltaC| = k12 mu^(1/2) W; "validity" = v^2/|V_Y|, the DeltaC at which the T-arc radial velocity sqrt(-V_Y DeltaC) equals v (my criterion; the first-order turn law is meaningless above it); DeltaC at the 2.2T2 extremum is 2 k12 mu^(1/2):

| (I, J, eps') | a | v | k12 | validity v^2/\|V_Y\| | DeltaC(2.2T2 extremum, W = 2) | ratio | mu at which that equals the validity limit |
|---|---|---|---|---|---|---|---|
| (2, 1, +1) | 1.5874 | 0.1705 | 0.09120 | 0.02483 | 0.0201 | 0.81 | 1.9e-2 |
| (3, 2, +1) | 1.3104 | 0.1121 | 0.06996 | 0.01131 | 0.0154 | 1.36 | 6.5e-3 |
| (4, 3, +1) | 1.2114 | 0.0838 | 0.05588 | 0.00647 | 0.0123 | 1.90 | 3.4e-3 |
| (3, 1, +1) | 2.0801 | 0.2326 | 0.07215 | 0.04389 | 0.0159 | 0.36 | 9.3e-2 |
| (3, 4, +1) | 0.8255 | 0.1120 | 0.11031 | 0.01412 | 0.0243 | 1.72 | 4.1e-3 |
| (1, 2, +1) | 0.6300 | 0.3577 | 0.54037 | 0.19915 | 0.1191 | 0.60 | 3.4e-2 |
| (2, 1, -1) | 1.5874 | 2.1705 | 0.32541 | 4.0248 | 0.0717 | 0.02 | 38 |
| (3, 2, -1) | 1.3104 | 2.1121 | 0.30362 | 4.0113 | 0.0669 | 0.02 | 44 |
| (1, 2, -1) | 0.6300 | 1.6423 | 1.15795 | 4.1991 | 0.2553 | 0.06 | 3.3 |

Reading: for the slow direct tangent resonances ((2,1), (3,2), (4,3), (3,4) with eps' = +1) the width of the transition-2.2 region at the Earth-Moon mass is comparable to or larger than (ratios 0.81, 1.36, 1.90, 1.72) the range of DeltaC over which the first-order T-arc turn law is valid at all, so the "regular" regime (DeltaC >> mu^(1/2) and DeltaC << v^2/|V_Y|) does not exist at the Earth-Moon mass for these; for (3,1) and (1,2) direct the ratios are 0.36 and 0.60 (a narrow regular window at best); for the retrograde ones it is wide (ratio 0.02 to 0.06). This agrees with the 1977-digest statement that at the Earth-Moon mass the first-order theory is outside its validity for most near-resonant critical generating orbits, and gives it a sharper, orbit-specific form. INFERRED caveat: L3 and k12 depend on my reading of L3 = [2 v/(3 pi I a)]^(1/4) (p.177), and the 2.2T2 extremum (W = 2) refers to n = 2 total bifurcations; other n have other W_min.

Recommended indeterminate-region test for #906, from this book (INFERRED, not coded): for each junction that has a nearby type 2 bifurcation ellipse, report nu_eff = ln|DeltaC|/ln(mu) and the pair (|DeltaC|/(mu^(1/2) k12), |DeltaC| |V_Y|/v^2). For a chain of n at least 2 basic arcs: if the first ratio is below about W_min (2 for the 2.2T2 case, sqrt(2) for 2.2P2; 2^(-1/2) is only the proven bound) or the second ratio is above 1, return "indeterminate" with those numbers; if nu_eff is at least 1/3 the S-arc and SS-node turn laws are outside their validity (R-region regime) and only the T-node law remains. For n = 1: 2T1, report the h_0 neck of section 2.6 rather than a DeltaC threshold; 2P1, return "indeterminate" when nu_eff is at least 2/3 (T-arc branches) or 1 (S-arc branches). Otherwise return the first-order law above as an expected value with its error scale (mu DeltaC^(-2), DeltaC^(1/2)).

## 9. Reconciliation with project code

COMPUTED (grep of `src/cyclerfinder` and `scripts` on 2026-10-04):
- `src/cyclerfinder/verify/turn_gate.py` implements the demanded-turn gate (arccos of the frame-invariant angle, available bend 2 asin(1/(1 + r_p v^2/GM)), `required_periapsis_alt_km` returning +inf for a zero turn). It has no notion of the bifurcation-neighbourhood regimes, of mu, or of DeltaC; the hyperbola relation it uses (gate relation tan(delta/2) = mu/(h v^2)) agrees with Henon's encounter relation (17.11) at small turn (section 2.1), so nothing in this book contradicts the gate's physics. It does not yet return "indeterminate", as the `#906` amendment requires.
- No module named `second_species_*` exists in `src/`; the only files mentioning second-species are `src/cyclerfinder/search/earth_moon_resonant_families.py` and `src/cyclerfinder/search/earth_moon_class1_resonant_connections.py` (docstrings saying the Casoliva Sec. IV elliptical/second-species method was never implemented) and `search/literature_check.py`. `bifurcation` appears in `search/variational_periodic_orbit.py`, `halo_family_at_jacobi.py`, `er3bp_*.py` and `earth_moon_resonant_families.py` (near-bifurcation continuation walls, bifurcation tracking noted as not deducible from one orbit's monodromy), none of which implements a generating-orbit or avoided-crossing analysis. So there is nothing in the code to reconcile against the type 2 formulas; they are all new.
- The project's `#899` plan (OUTSTANDING) names Henon 1968 Tables 1 to 9, Bruno 1981 Tables I and III and Perko 1981's (2,1) resonance bifurcation (C = -0.406767) as controls. This book adds controls that live at the bifurcation neighbourhood (below) and makes the type 2 / type 1 split explicit: the (2,1) resonance orbit of Perko's examples (C = -0.406767, a type 1 orbit, volume I Table 6.2 first row) is NOT a type 2 orbit, so the type 2 formulas here do not apply to it.

## 10. Techniques applicable to the project's problems

### (a) `#899`: continuation of second-species chains from small mu toward the Earth-Moon value; where Casoliva et al. report lunar impact

1. A mu-continuation at fixed |DeltaC| (or approximately fixed period, as Casoliva et al. do) crosses the nu regimes. At fixed DeltaC the effective exponent nu_eff = ln|DeltaC|/ln(mu) rises from about 0.2 at mu = 1e-6 to 0.3 to 0.7 at mu = 0.0122 for |DeltaC| of 0.05 to 0.01, so a family member starting well inside the nu < 1/3 regime passes through transitions 2.1 and 2.2 on the way (if it lies near a type 2 ellipse at all; see the C-value caveat in item 2), which are exactly where branches reconnect, characteristics have minima and S-arcs fuse. INFERRED, and not the same as the lunar-impact mechanism, but it predicts where a mu-continuation can turn back: for a minus-sign branch a fold in mu occurs where |DeltaC| = k mu^(1/2) W_min, i.e. mu = (|DeltaC|/(k W_min))^2; this is what "families do not exist uniformly in mu" (Casoliva et al. 2010, IV.C) looks like in Henon's variables. The table of section 8 gives k and W_min for the 2.2T2 case; the prediction has not been tested.
2. Lunar impact: Casoliva et al. (2010, IV.C) report fixed-period continuation "in most cases" ending on lunar impact; the project's Casoliva digest (`docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md`) gives no Jacobi constant for those Class 1 impacts. The impact energies it does list (line 133: h = -1.4711 for the He1 and Hm1 families and -1.5892 for Hm2, section V.C) belong to the Class 2 L1-homoclinic families, not to the resonant cyclers; with the digest's convention C_J = -2h (its eq. 3, flagged there as read from an OCR of the paper) they are C_J = 2.9422 and 3.1784, in the Hill-region range near C = 3 where the slow direct tangent ellipses of section 8 also lie (C = 2.87 to 2.99), but they are not resonant-ellipse orbits and no relation to a type 2 ellipse is established here. The Class 1 cyclers of 2010 Table 3 have C_J from 0.4887 (2-1a) to 2.7630 (1-2e), well below the direct type 2 values (2.87 to 2.99) and above the retrograde ones (-1.71, -1.46 for (2,1), (3,2); 0.30 for (1,2)); so the type 2 transitions of this part probably do NOT govern the Casoliva Class 1 continuation, which is near type 1 ellipses (volume I Table 6.2, e.g. C = -0.4068 for (2,1,0), and chapters 12 to 14 of this book, other agent's part). INFERRED from C values alone; not tested.
   Earlier text of this item, kept for the argument: Casoliva et al. report fixed-period continuation "in most cases" ending on lunar impact. In Henon's variables the impact is the periapsis of the local hyperbola at an encounter dropping to the Moon's radius: r_p = h-order for small turn (r_p tends to h, which is of order mu^(3/4) at nu = 1/2 and mu DeltaC^(-1/2)/v at a T node). The book says nothing about impact or a finite Moon radius (it is a point-mass theory); but it gives the miss distance as an explicit function of DeltaC: at the minimum of W the miss distance of the closest encounter is the extremal h of 17.208, so a seed with small DeltaC has a predictable miss distance of order mu^(3/4) at nu = 1/2 (coefficient from 17.208 not transcribed). COMPUTED scale: mu^(3/4) = 0.0366 at the Earth-Moon mass, the Moon's radius 1737 km is 0.00452 of the Earth-Moon distance, so the dimensionless coefficient of the closest encounter would have to exceed 0.12 to clear the surface; this is a scale argument, not a result. I cannot say from the book which Casoliva rows are the impact cases.
3. The 2T1 hyperbola (23.25) is a published, closed-form positive control for a continuation driver: the minimum separation of the first-species and second-species branches at a type 2 orbit scales as sqrt(16 v mu/(3 pi I a)). A test of the driver: continue a first-species family past (I, J) = (2, 1) with eps' = +1 at three values of mu (1e-6, 1e-5, 1e-4), measure the closest approach in the h_0 variable (impact distance at conjunction) at DeltaC = 0, and verify the mu^(1/2) slope and the coefficient (a/(4 |a - 1|)) sqrt(16 v/(3 pi I a)) = 2.04 for (2, 1, +1) (COMPUTED: 0.0225/sqrt(0.0121529529)) in the h_0 variable. Not run.
4. The R-region characteristics (y_i, w) are a ready-made hierarchy for the number of sign-labelled branches to expect: 2^(n-1) R-arcs and 2^n - 2 R-orbits (exact counts), a test of any enumerator of second-species chains near a type 2 orbit (all branches are found if and only if the count matches).
5. Seeds for chains with more than one flyby per period: the plus-sign branches (S-arc compositions) are pure combinatorics of the S-arc alphabet; the book gives which compositions of n reconnect (Table 19.1: pairing, for example n = 4: (22, 1111), (4, 121), (211, 112), (31, 13)), so a chain search that follows one composition through a mu-continuation can predict which composition it arrives on after a transition 2.1. INFERRED use; the pairing is for the generating chains at small mu, not a proven statement at mu_EM.

### (b) `#906`

Section 8 above: the explicit first-order turn laws, the saturation scales mu^(1/3), mu^(1/4), the minimum |DeltaC| of order k mu^(1/2) W_min on every characteristic with n at least 2, the separate n = 1 statements (2T1: no minimum in DeltaC; 2P1: indeterminate beyond nu_eff = 2/3 for T-arcs and 1 for S-arcs), the orbit-specific table, and the recommended return values (nu_eff, two ratios). It also reinforces the amendment: zero turn at the bifurcation orbit is a genuine, regular limit of the family (a zero-turn junction between identical T-arcs is excluded for mu above 0 by Proposition 4.3.2 of volume I; 2T1 shows the exclusion as an avoided crossing of width sqrt(mu)), so "zero turn" must be flagged as "bifurcation orbit", not rejected and not accepted as an encounter.

### (c) `#928`

The book offers no regularisation, but chapters 17.2 to 17.7 define the scale on which a regularised propagator must resolve a near-resonant family: the impact distance h of order mu^(3/4) and relative speed v at the Moon, so the time spent within h is of order h/v; at mu = 1e-6 (the Casoliva seed mass) that distance is 3e-5 of the lunar distance (about 12 km) and at 1e-8 it is 1e-6 (0.4 km) (COMPUTED mu^(3/4) values; order only, coefficients not included). A regularised integrator and transition matrix are therefore needed for the small-mu end of a continuation. The exponents here also set the test cases for `#928`'s plain, Sundman and Levi-Civita comparison: a type 2 orbit with DeltaC chosen so that nu_eff is 1/3, 0.4 and 1/2 at mu = 1e-8, 1e-6, 1e-4 has encounter distances following the table of section 2.5, which the regularised integrator must reproduce. Not run.

## 11. Printed numbers usable as sourced tests

All of these are printed values, transcribed from the page images. Items 1 to 3 are verified numerically (COMPUTED) as described.

1. Table 18.2 (p.185, 31 values, 8 digits): the R-arc solutions of y_{i-1} - 2 y_i + y_{i+1} + 1/y_i = 0 with y_0 = y_n-tilde = 0, one for each sign sequence. Independent reproduction by multi-start Newton: all 31 match to 8 digits; there are exactly 2^(n-tilde - 1) solutions. Expected: n-tilde = 2: y_1 = 0.70710678; n-tilde = 3: 1.00000000, 0.57735027; n-tilde = 4: 1.17914724, 0.94010422, 0.59967641, 0.53185593.
2. Table 18.3 (p.189, rows listed in section 3, 6 decimals): R-orbit solutions of the cyclic recurrence; residual of each printed row 4e-6 or less; all 2^n - 2 sign sequences except +++ and --- have exactly one solution (n up to 6, stated by the book; proved by the baker map, 18.3).
3. Closed forms: n-tilde = 2, y_1 = 1/sqrt(2); n-tilde = 3, y_1 = 1 and 1/sqrt(3); the n = 5 R-orbit polynomial 5 y^6 - 20 y^4 + 17 y^2 - 4 = 0, with roots 0.652966, 0.799673, 1.712938 (COMPUTED, matching Table 18.3 rows 5); the n = 6 closed values (Table 18.3), the n = 7 numerical orbit of p.191.
4. Extremum values of w (nu = 1/3): w^3 = -4 (2.1P2, equation y^2 + w^2 y - w = 0, which I verified with sympy gives the vertical tangent dF/dy = 0 at w = -2^(2/3), i.e. w^3 = -4); w^3 = -2 and -8/3 (2.1P3, the first COMPUTED with sympy from y^2 + 2 w^2 y - 2 w = 0 giving w = -2^(1/3), i.e. w^3 = -2; the -8/3 is printed only); w^3 = -2 (2.1T2); w^3 = -4/3 and -3/2 (2.1T3); w^3 = -1, -9/4, -8/9 (2.1T4). Extremum values of W (nu = 1/2): sqrt(2) and 3/2 (2.2P2), 2 (2.2P3), 2 and 3/sqrt(2) (2.2T2), sqrt(6) and ((3 + sqrt(33))/144)^(1/4) in Y_0 (2.2T3).
5. The 2.2T3 symmetric closed form (22.27) satisfying (22.30); verified to 3e-15.
6. 2T1: y^2 + w y - 1 = 0 and equations (23.25); Table 23.1 junctions (+E with -1; -E with +1) and Table 23.2 (+1 with -1; -f with -g).
7. Constants of a type 2 ellipse (17.23 to 17.27) with C = 3 - v^2 and the agreement with volume I Table 4.4 as above (to 6e-6).
8. W larger than 2^(-1/2) (21.7); the type 3 transitions list (23.53: 1/5, 2/9, 1/4, (n-1)/(4n-5), 1/3, (n-1)/(2n-1), 1/2, 2/3); the ranges nu < 2/3 (2P1 T-arcs) and nu < 1 (2P1 S-arcs).
9. Junction tables: Table 19.1 (complete transcription, pairing verified), Table 20.1 (n up to 4 transcribed), Tables 21.3, 22.2 (n = 2 partially), Table 23.1 and 23.2.
10. The disagreement with Guillaume (1971): +4 joins +-22- and +-1111- joins +121 (p.233).
11. Unit of cross-check for the mu = 0 side from volume I: the type 2 ellipse C values above against Table 4.4.

Not usable as tests (no printed numeric table): Figs. 17.2, 19.1 to 19.12, 20.1 to 20.5, 21.1 to 21.9, 22.1, 22.2, 23.1, 23.2 (characteristic sketches; Fig. 22.2 shows Y_0 against W for the 2.2T3 symmetric branches and could be digitised); the book contains no computation at any stated small mu, no comparison with a numerically computed periodic orbit of the full problem, and no table in terms of the Earth-Moon mass.

## 12. Printed values flagged as unreadable or not transcribed

- L2 (17.168, p.173) and the exact factor in 23.14 (p.285), the 21.83 intersection expression (p.252), the n-tilde = 4 closed form of 18.14 (p.184), the 2.1P3 factor of 19.65 (p.211), the stability formula z for 2.1T2 (20.19) and the 18.57 closed forms: text-layer garbled; I did not read them on images and did not transcribe them. The tables above use only the verified values.
- Tables 21.3 and 22.2 beyond n = 2: not transcribed (the headers could not be matched to pairs at the resolution used).
- The exact relation of the "(I, J) of the book" to the project's p:q labelling: volume I section 3 gives (p, q) = (J, I) (q is the Moon's revolutions). Henon's I counts the Moon's revolutions in the basic-arc period 2 pi I; this is consistent with a = (I/J)^(2/3).

## 13. Follow-ups (task numbers not registered)

1. Build the 2T1 avoided-crossing control: compute the first-species and second-species branches near the (2,1) direct tangent ellipse at mu = 1e-6, 1e-5, 1e-4 with the project's corrector and compare the h_0 minimum with the closed form of 23.25 (slope 1/2, coefficient 2.04 times mu^(1/2) for (2,1,+1)).
2. Implement the "indeterminate" return of the gate with the two ratios of section 8 (nu_eff, |DeltaC|/(k mu^(1/2) W_min), |DeltaC| |V_Y|/v^2), after the other digests (Gomez & Olle, Perko) confirm the validity numbers.
3. Add the R-arc and R-orbit solvers (the recurrence 18.4/18.17) and Tables 18.2 and 18.3 as tests of any second-species enumerator (the counts 2^(n-tilde - 1) and 2^n - 2 are exact); this is a small, sourced module with an independent check already done here.
4. Verify the type 2 constants (17.23 to 17.27) in a test against volume I Table 4.4 (done here by hand at six points to 6e-6).
5. Read chapters 12 to 14 of this book (other agent) and settle whether the type 1 nu above 1/2 statement matches Perko 1981's mu^(1 - nu) law; if yes, record it as the answer to Perko 1981's private-communication remark.
6. Decide whether Table 19.1 and the junction tables should be encoded as a test of a branch-pairing predictor (none exists in `src/`), and transcribe the remaining Table 20.1, 21.3 and 22.2 columns after re-reading the page images at higher resolution.
7. Obtain Guillaume 1971 (Liege thesis) and Devaney 1981 if the 2T1 cross-check or the baker-map proof of the R-region counts is wanted; neither is needed for a first build.
8. Digitise Fig. 22.2 (Y_0 against W for the 2.2T3 symmetric branches) if a graphical control is wanted; the closed form is already exact.
9. Quantify, with the project's CR3BP corrector, how far from the mu = 0 generating orbit the Earth-Moon Casoliva rows lie in DeltaC and nu_eff, to see which are inside the indeterminate region of section 8 (needs `#899` seeds first).

## 14. References in this part that the project does not hold (checked against `docs/notes/CORPUS_INDEX.md` on 2026-10-04)

Held, with this book's labels matched to the corpus file names: book Perko 1976b (Celest. Mech. 14:395) = corpus `perko-1976`; book 1977a (Celest. Mech. 16:275) = corpus `perko-1977`; book 1981a (SIAM J. Appl. Math. 41:181) = corpus `perko-1981b-periodic-orbits...` (the labels are swapped); book 1981b (Celest. Mech. 24:155) = corpus `perko-1981-second-species...`; the corpus file `perko-1976b` (Rocky Mountain J. Math. 6:675) is NOT cited by this book. Also held: Guillaume 1973 (= 1973b here), Guillaume 1975 (= 1975b here, Celest. Mech. 11:449) and 1975a (Celest. Mech. 11:213), Hitzl & Henon 1977 (Celest. Mech. 15:421, which this book does not cite), Bruno 1981, Brjuno 1978, Henon 1968, Henon 1997, Marco & Niederman 1995, Casoliva et al. 2008 and 2010, Barrabes & Gomez 2002 and 2003, Broucke 1968.

Not held, and relevant to `#899`:
- Guillaume, P., "Solutions periodiques symetriques du probleme restreint des trois corps pour de faibles valeurs du rapport des masses", Ph.D. thesis, Universite de Liege, 1971 (cited for the 2T1 equation IV-33, the 2.1T2 and T3 characteristics, the 2T4 disagreement and Delta a equal within O(mu^(2/3)) in R-regions): the most direct earlier treatment of what this book generalises.
- Devaney, R. L., "The baker transformation and a mapping associated to the restricted three body problem", Commun. Math. Phys. 80:465-476 (1981): the proof that the R-region map is conjugate to the baker map (not needed to build, needed only to cite the proof).
- Hitzl, D. L. & Henon, M., "The stability of second species periodic orbits in the restricted problem (mu = 0)", Acta Astronautica 4:1019-1039 (1977): stability data for the same generating orbits (not the Celest. Mech. paper that is held); also not held per the volume I digest.
- Bruno, A. D., "Periodic solutions of the second kind in the restricted three-body problem" (Moscow preprint 95, 1976) = chapter 7 of Bruno 1994, The Restricted 3-Body Problem: Plane Periodic Orbits, de Gruyter: the Newton-polyhedron (power geometry) route and the second-kind solutions; and Bruno 2000, Power Geometry in Algebraic and Differential Equations (Elsevier): the method behind chapter 15.
- Henon, M. & Guyot, M., "Stability of periodic orbits in the restricted problem" (1970, Reidel): the stability index conventions used here.
- Perko 1965 (thesis): asymptotic matching; low priority.
- Henon's manuscript notes (about 400 pages in French) on type 3, offered in the book to interested colleagues (23.4); relevant to the C = -1 retrograde-circle bifurcation, not to `#899` directly.

## 15. Note 2026-10-05: Devaney 1981 now held and digested

The baker-map proof cited in chapter 18 (Devaney 1981, Commun. Math. Phys. 80:465) is digested in `docs/notes/2026-10-04-digest-devaney-1981-baker-transformation.md`. It confirms the dictionary used in sections 3 and 4 above: the book's recurrence is Devaney's plane map in the coordinates (x_i, y_i); Corollary B gives exactly 2^n - 2 period-n points, all hyperbolic (the R-orbit count and Proposition 18.2.2); the R-arc count 2^(n-tilde - 1) corresponds to two-sided terminating sequences, whose existence Devaney only sketches (the book proves it independently, 18.1.2). I also verified the R-orbit counts for n = 2 to 7 (2, 6, 14, 30, 62, 126, one per sign sequence) and Henon's n = 7 orbit (residual 2.5e-6). Item 7 of section 14 (Devaney not held) is now out of date; the thesis and the other references there remain not held.
