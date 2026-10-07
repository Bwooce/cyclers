# Digest: Lantoine & Russell 2011, "Complete closed-form solutions of the Stark problem" (#960 batch 37, #999, #864)

G. Lantoine and R. P. Russell (Georgia Institute of Technology), Celestial Mechanics and Dynamical Astronomy
109(4):333-366 (2011), doi 10.1007/s10569-010-9331-1 (Crossref: CMDA 109(4), pp.333-366, 2011-04; title and
authors match the title page). Received 16 Apr 2010, revised 25 Aug 2010, accepted 16 Dec 2010, online
5 Feb 2011 (dates from the text layer).
- Supplied file `7ce4bdca-lantoine2011_2.pdf`, 34 pp. = journal pp.333-366 (PDF page n = p.332 + n). md5
  ce4366b2e54cb1dc10eacc0c11dbfe06. Springer text layer (symbols, roots and fractions garbled).
- **Proposed corpus filename:** the corpus already holds a different Lantoine & Russell 2011 paper
  (`lantoine-russell-2011-near-ballistic-halo-to-halo-transfers-planetary-moons-jas-58-335-...`). Following the
  `perko-1981b` convention:
  `cyclers_pdf/papers/lantoine-russell-2011b-complete-closed-form-solutions-stark-problem-cmda-109-333-doi-10.1007-s10569-010-9331-1.pdf`
- **How I read it:** the whole text layer. On page images (150 dpi; Fig. 11 at 300 and 600 dpi): p.337 (eq. 1),
  p.339 (eqs 8-12), p.349 (Fig. 5, eqs 70a-f, footnote 2), p.350 (Table 1, Fig. 6 caption), p.351 (eqs 71-72),
  p.361 (eqs 127-130), p.363 (Fig. 11 and its caption, eqs 132-134), p.364 (Table 2, all cells). The 3D
  quadrature (eqs 86-125) and the planar root formulas (eqs 13-69) were read in the text layer only; the
  3D cubic (eq. 89a) is checked indirectly by the classification test below.
- **Checks:** `check_stark.py` -> `check_stark.out` (secs A-E); `check_fig11_axis.py` -> `.out` (Fig. 11 axis). DOP853, rtol 1e-13, atol 1e-14, our own
  integration as the reference. No value our code computed is offered as a golden.
- Wanted list: not listed (supplied by the owner for `#999`). Proposal for the propagator and the golden tests:
  `stark-plan.md` (same folder).

## 0. Verdict

**The reference closed-form solution of the Stark problem (Kepler plus a constant inertial acceleration),
planar and spatial, in Jacobi elliptic functions. Exact, but with 7 planar and 3 spatial solution forms, a
Newton inversion for physical time, and no partial derivatives. It publishes no output states, so it gives
classification goldens and physics checks, not state goldens.**
- **Model:** r'' = -mu r/r^3 + eps k-hat (eq. 1 planar, eq. 81 spatial). Separable in parabolic coordinates
  (xi^2 = y + r, eta^2 = -y + r for a field along +y, eq. 3) with the Sundman-type time dt = (xi^2 + eta^2) dtau
  = 2r dtau (eq. 4). Two integrals in the plane: H (eq. 2) and the separation constant c (eqs 8-9). In space a
  third: p_phi, the angular momentum about the field axis (eq. 88).
- **Solutions:** five xi forms (xi1-xi5) and two eta forms (eta1, eta2) give 7 planar orbit types in 6 domains
  of the (H/sqrt(mu eps), c/mu) plane (Fig. 5, Table 1, eqs 70a-f). Only xi1-eta2 is bounded. In space the
  sextic in xi (eq. 89a) is reduced through a real root Y* of a cubic in Y = X^2 to the planar forms; three
  spatial types (xiI, eta), (xiII, eta), (xiIII, eta) (eqs 110-125); (xiI, eta) is the bounded one.
- **Small-eps safe:** the roots are rewritten so the forms tend to the Kepler sin/sinh forms as eps -> 0
  (eqs 24-27, 44, 65). For the time equation (eqs 73-79) the k -> 0 limit needs a Taylor expansion in the
  modulus (eq. 80, order 2 shown; order 4 used in Table 2).
- **Time:** t(tau) is closed form (elliptic integrals of the second kind, eqs 73-79; the spatial case adds
  tau2(tau) through elliptic integrals of the third kind, eqs 132-134). tau(t) needs a Newton iteration
  ("convergence is usually obtained in few iterations", p.354, text layer).
- **Accuracy (Table 2, p.364):** over 20 TU with mu = 1, the analytic position differs from a quad-precision
  RKF7(8) truth by 8.9e-15 to 3.9e-14 (planar) and 5.7e-13 to 2.1e-12 (spatial) relative; analytic
  Hamiltonian error 1.1e-16 to 3.5e-14 (planar), 1.3e-15 to 4.1e-14 (spatial). No timings in this paper.
- **What it gives the project:** (1) the exact integrals H, c, p_phi as conserved-quantity tests for any Stark
  propagator; (2) a classifier (bounded or not) from the initial state alone; (3) ten published initial
  conditions, one per orbit type, as a coverage set; (4) exact periodic displaced circular orbits (eqs
  127-130); (5) the eccentricity-excitation law for inclined orbits (Fig. 11, eq. 72). It does not give a
  ready propagator: Hatten & Russell 2015 found that the full algorithm needs both 2D and 3D formulas
  (the 3D formulas do not reduce to 2D), three kinds of elliptic integral and the t inversion.
- **Two printed slips found (check them before coding from this paper):**
  - **Eq. 9 is not conserved as printed.** With the field along +y, the printed
    c = x'(y x' - x y') - mu y/r - (1/2) eps x^2 drifts by 57 over 20 TU for the footnote-2 state (eps = 0.4).
    Hatten & Russell 2015 eq. 25, c = x'(x y' - y x') + mu y/r - (1/2) eps x^2, is constant to 4e-13 and equals
    the c of eq. 8; the xi side and the eta side of eq. 8 agree to 3e-15 (`check_stark.out` sec. B, three
    states). Algebraically the printed eq. 9 equals -c - eps x^2 (checked: eq. 9 + eps x^2 is constant and equals
    -c). So two one-step repairs give a conserved quantity: flip the sign of the (1/2) eps x^2 term (giving -c),
    or flip the first two terms (giving +c, the sign used in eq. 8 and by Hatten). Either way, use eq. 8 or
    Hatten eq. 25 for c. (Eqs 8, 9 and Hatten 25 read on the page images.)
  - **Eq. 127 prints five entries,** X0 = [rho, 0, z, rho w, 0]. The 6-vector must be [rho, 0, z, 0, rho w, 0]
    (velocity rho w along +y). With that reading one period 2 pi/w returns to the start to 1e-12 (sec. D).
- **One unstated criterion and one axis label:**
  - p.361 (image) says that "these circular orbits" (the class) are the least-energy orbits of the Stark
    problem (Namouni & Guzzo 2007, not held), then gives eq. 128 for z "for a given mu and eps" without saying
    how z is chosen. Eqs 129-130 are an equilibrium for any z (eps r^3 = mu z). On that family H is unbounded
    below as z -> 0 and has one stationary point, a maximum, at z = sqrt(mu/(27 eps)), twice the printed eq. 128
    value (sec. D). Fig. 10's z-axis (0.9623) is consistent with eq. 128 as printed at an inferred eps = 0.01
    (z = 0.962250; Fig. 10 does not print eps). The criterion behind eq. 128's z is not stated; unresolved.
  - Fig. 11 (eps = 0.0103, X0 = [1, 0, 0, 0, 0.866, 0.5]): the integration reproduces e in [0, 0.5001] =
    [0, sin 30 deg] and the periapsis argument sitting at 0 or 180 deg with flips. But our first eccentricity
    maximum is at t = 102.5 TU; the figure shows it at about 51 TU. The figure's horizontal axis is the
    fictitious time tau (dt = 2r dtau), not t: its a(DU) panel shows 16 oscillations between axis 0 and 50,
    spacing 3.14 axis units (`check_fig11_axis.out`, 600 dpi), while our osculating a(t) has 8 maxima in
    t = 0-50 TU, spacing 2 pi (sec. E); with r about 1, one orbit spans pi in tau. Our e(t) is periodic with
    period pi/w_s (our derivation from the integration), w_s =
    3 eps/(2 n a) = 0.01545 rad/TU (eq. 72); its FFT peak is 2 w_s = 0.031 rad/TU, while the paper reports
    "a main angular frequency of 0.015 rad/TU" (p.362, text layer), which is w_s itself.
- **Also:** Fig. 6's caption carries a footnote mark 3 with no footnote text on pp.350-351 (page images).
  Sections 2.1-2.5 derive with the field along +y; Sec. 4, Fig. 6 and footnote 2 use the field along +x
  (sec. A below shows only +x reproduces the labels).
- **Catalogue implication (PROPOSAL only):** none. Code implication: see `stark-plan.md`.

## 1. Content

- **Sec. 1 (pp.333-337).** Motivation (low-thrust arcs discretised into constant-thrust segments; solar
  pressure). History: Lagrange 1788 (quadratures), Jacobi and Liouville (separability in parabolic
  coordinates), Born, Epstein (quantum Stark effect), Isayev & Kunitsyn 1972, Beletsky 2001 (planar Jacobi
  forms, singular as eps -> 0, missing xi4 and xi5), Vinti 1964 (secular only), Kirchgraber 1971 (KS
  variables, no explicit forms), Rufer 1976 (low-thrust optimisation with Kirchgraber's solution),
  Poleshchikov 2004 (KS, partly series), Rauch & Holman 1999, Namouni 2005, Namouni & Guzzo 2007.
- **Sec. 2 (pp.337-354), planar.** Eq. 1 (field along +y, image p.337). H (eq. 2). Parabolic coordinates
  (eq. 3), dt = 2r dtau (eq. 4). Separation (eqs 6-8, image p.339): H xi^2 - (1/2) xi'^2 + mu + (1/2) eps xi^4
  = -H eta^2 + (1/2) eta'^2 - mu + (1/2) eps eta^4 = -c. Quadratures (eqs 11-12): P_xi = eps xi^4 + 2H xi^2 +
  2(c + mu), P_eta = -eps eta^4 + 2H eta^2 - 2(c - mu).
  - xi cases (Bowman 1961 reductions): discriminant Delta_xi = (2H)^2 - 8(c + mu) eps (eq. 14). Case A.1 (both
    roots psi+- > 0): xi1 (bounded, xi0^2 < xi2^2, sn form, eq. 31) or xi2 (xi0^2 > xi1^2, 1/sn form, eq. 35);
    Case A.2 (one positive root): xi3 (1/cn, eq. 40); Case A.3 (no positive root): xi4 (sn/cn, eq. 46);
    Case B (Delta_xi < 0): xi5 (Cayley reduction, eqs 49-57).
  - eta cases: Delta_eta = (2H)^2 - 8(c - mu) eps (eq. 58); eta1 (dn form, eq. 63) or eta2 (cn form, eq. 68);
    A.3' and B' are infeasible.
  - Summary (Sec. 2.4): 7 orbit types (Fig. 6), domains I-VI (Fig. 5); boundary curves (eqs 70a-f, image
    p.349): B1 c/mu = -1 + H^2/(2 eps mu); B2 c/mu = 1 + H^2/(2 eps mu); B3, B4 c/mu = -1 (H > 0, H < 0);
    B5, B6 c/mu = 1 (H < 0, H > 0). Small eps sends |H/sqrt(mu eps)| large (far left or right of Fig. 5).
  - Bounded orbit (xi1-eta2): confined between two parabolas (eqs 71a-b), a precessing ellipse whose
    eccentricity reaches 1 in the plane; Stark frequency w_s = 3 eps/(2 n a) (eq. 72, image p.351; citing Hezel
    et al. 1992, Namouni & Guzzo 2007).
  - Stark equation (Sec. 2.5): integrals of xi^2 dtau and eta^2 dtau for each form (eqs 73-79, via
    Mathematica); not valid at k = 0 or 1 (the Kepler limit): use the modulus Taylor expansion (eq. 80).
- **Sec. 3 (pp.354-363), spatial.** Field along +z (eq. 81); 3D parabolic coordinates x = xi eta cos phi,
  y = xi eta sin phi, z = (xi^2 - eta^2)/2 (eq. 83); two times dt = (xi^2 + eta^2) dtau1 = xi^2 eta^2 dtau2
  (eq. 84). Integrals H, c (eq. 87), p_phi (eq. 88). The quadratures have sextics (eq. 89); Y = X^2 makes them
  cubic; a real root Y* and Z^2 = sign(a)(Y - Y*) reduce them to the planar forms (eqs 90-125). Cubic
  discriminant Delta = 4(b^2 - 3ac)^3 - e^2 with e = 2b^3 - 9abc + 27a^2 d (eqs 94-95).
  - Examples: displaced circular ("sombrero", "static") orbits (eqs 127-131); excited inclined orbits whose
    maximum eccentricity is the sine of the inclination of the force to the angular momentum (Namouni 2005),
    e.g. a 500 km, 97.8 deg terminator orbit has an eccentricity period of about 62 days (p.362, text layer only).
  - tau2(tau1) for the bounded case through incomplete elliptic integrals of the third kind (eqs 132-134).
- **Sec. 4 (pp.363-364), validation.** Field along +x (2D) or +z (3D), 20 TU, mu = 1; RKF7(8) in quad
  precision (tol 1e-21) as truth and in double (tol 1e-16) as a comparison; analytic in double. Table 2 below.
- **Sec. 5.** Future work: equilibria and periodic orbits; low-thrust segments; analytic partials for
  optimisation (not given here).

## 2. Table 2 (p.364, page image; every cell matches the text layer)

| Type | Initial conditions | eps | H err quad RKF | H err dbl RKF | H err analytic | rel. pos. dbl RKF | rel. pos. analytic |
|---|---|---|---|---|---|---|---|
| xi1 eta2 | [1, 0.1, 0.05, 1] | 1e-9 | 3.3e-20 | -3.1e-14 | -4.5e-16 | 4.0e-15 | 3.6e-14 |
| xi2 eta2 | [10, 1, 0, 0.1] | 0.1 | 8.7e-22 | -1.6e-15 | -2.0e-16 | 1.5e-15 | 1.2e-14 |
| xi3 eta2 | [10, 1, 0, 1] | 0.001 | 5.8e-22 | 6.2e-15 | 4.8e-15 | 6.4e-16 | 1.1e-14 |
| xi4 eta2 | [1, 1, 1, 1.4] | 0.001 | 4.5e-22 | 1.0e-15 | 1.7e-16 | 1.3e-16 | 3.9e-14 |
| xi4 eta1 | [0.2, 1, 1, 1.4] | 0.01 | 1.1e-21 | 5.7e-15 | 1.1e-16 | 1.2e-15 | 9.0e-15 |
| xi5 eta2 | [0.2, 1, 0, 1.4] | 0.01 | 1.5e-19 | -4.8e-13 | 9.3e-15 | 4.1e-15 | 8.9e-15 |
| xi5 eta1 | [0.33, 1, 1.01, 1.09] | 0.035 | 6.4e-21 | -5.6e-14 | 3.5e-14 | 3.5e-15 | 1.1e-14 |
| (xiI, eta) | [1, 0, 0, 0, 0.1, 0.1] | 1e-9 | 3.6e-20 | -3.5e-13 | -1.3e-15 | 5.2e-14 | 2.1e-12 |
| (xiII, eta) | [0, 0.8, 1, -0.8, 0, 0] | 0.5 | 1.3e-19 | -8.3e-15 | -2.4e-15 | 3.7e-14 | 5.7e-13 |
| (xiIII, eta) | [0, 0.8, 1, -0.8, 0, 0] | 2 | 2.4e-19 | -1.3e-12 | 4.1e-14 | 9.3e-14 | 7.8e-13 |

- Units: mu = 1 (DU^3/TU^2), eps in DU/TU^2, time of flight 20 TU. No final states are printed.
- **Convention, decided by `check_stark.out` sec. A:** 2D state [x, y, x', y'] with the field along +x; 3D
  state [x, y, z, x', y', z'] with the field along +z. With these, our classifier (eqs 14-15, 58-59, 89a,
  94-101) puts 7/7 planar and 3/3 spatial rows in their printed types. With the field along +y (as in eq. 1)
  only 2/7 planar rows match. Over 20 TU our integration holds H to <= 1.2e-11 and c (from eq. 8 or 87, either
  side) to <= 2e-10.
- **Numerical note for implementers:** for the (xiI, eta) row (eps = 1e-9) eq. 95 evaluated in double gives
  Delta = -8.5e-14 (wrong sign: it says one real root), because eq. 95 is 27 a^2 times the ordinary
  discriminant and is formed as a difference of two O(240) numbers. The ordinary discriminant
  18abcd - 4b^3 d + b^2 c^2 - 4ac^3 - 27a^2 d^2 = +15.37 gives three real roots and the right type. Any
  implementation must not use eq. 95 as printed for small eps.
- **Footnote 2 (p.349, image):** (mu = 1, eps = 0.4, x0 = 1, y0 = 0.1, x0' = 0.05, y0' = 1) and (mu = 1,
  eps = 0.7, x0 = 1, y0 = 0.1505, x0' = 0.1137, y0' = 1) share (c/mu, H/sqrt(mu eps)). With the field along +x
  both give c/mu = -0.0020 and H/sqrt(mu eps) = -1.4132 (our arithmetic) and the types xi1-eta2 and xi2-eta2,
  as the footnote says (sec. C). With the field along +y they differ.

## 3. Golden candidates (published values only; tests built on them are physics checks unless stated)

- **G-L1 classification set:** the 10 Table 2 rows (inputs, eps, printed type), with the convention above.
  Test: the classifier returns the printed type. Published-value test (the type labels are published).
- **G-L2 footnote 2 pair:** equal (c/mu, H/sqrt(mu eps)) and types xi1-eta2 / xi2-eta2.
- **G-L3 integrals:** H (eq. 2/82), c (eq. 8, or Hatten eq. 25; NOT Lantoine eq. 9 as printed), p_phi (eq. 88)
  constant along any Stark arc. A propagator test, not a value golden.
- **G-L4 boundary curves:** eqs 70a-f (Fig. 5) for classification edge tests.
- **G-L5 displaced circular orbit (physics check, truth = our integration):** eqs 127-130 with
  X0 = [rho, 0, z, 0, rho w, 0]: periodic with period 2 pi/w. Fig. 10's z = 0.9623 DU is consistent with
  eq. 128 at an inferred eps = 0.01, mu = 1 (eps is not printed).
- **G-L6 inclined-orbit excitation:** Fig. 11 inputs (mu = 1, eps = 0.0103, X0 = [1, 0, 0, 0, 0.866, 0.5]):
  published: e oscillates between 0 and sin 30 deg (p.362 and the figure), argument of periapsis at 0 or
  180 deg; w_s = 3 eps/(2 n a) (eq. 72). The figure's time axis is tau, not t (see Verdict).
- The Table 2 error columns are not goldens (they describe their implementation).

## 4. Citation mining

| Cited work | Held? | Wanted list |
|---|---|---|
| Beletsky 2001, Essays on the Motion of Celestial Bodies | not held | not listed |
| Biscani & Izzo 2014 (cited in Hatten 2015, not here) | not held | not listed; see Hatten digest |
| Bowman 1961, Introduction to Elliptic Functions | not held | not listed (textbook; low) |
| Cordani 2003, The Kepler Problem | not held | not listed |
| Dankowicz 1994, CMDA 58:353 | not held | not listed |
| Forward 1991, JSR 28:606 ("Statite") | not held | not listed |
| Hezel et al. 1992, Am. J. Phys. 60:324 | not held | not listed |
| Isayev & Kunitsyn 1972, CMDA 6:44 | not held | not listed |
| Kirchgraber 1971, CMDA 4:340 | not held | not listed |
| Lantoine & Russell 2009, ISSFD 21 ("The Stark Model ... low-thrust trajectory optimization") | not held | not listed; **candidate row** (the low-thrust use of this solution) |
| McInnes 1998, JGCD 21:799 (displaced non-Keplerian orbits) | not held | not listed |
| Namouni 2005, AJ 130:280; Namouni & Guzzo 2007, CMDA 99:31 | not held | not listed; Namouni & Guzzo is cited for the least-energy property of the displaced circular orbits and may give the criterion behind eq. 128 (optional) |
| Poleshchikov 2004, Cosmic Res. 42:398 | not held | not listed |
| Rauch & Holman 1999, AJ 117:1087 | not held | not listed |
| Redmond 1964, Phys. Rev. 133:B1352 ("1954" in the text, "1964" in the list; text layer) | not held | not listed |
| Rufer 1976, CMDA 14:91 (low-thrust optimisation with the closed form) | not held | not listed; low |
| Stump 1998, Eur. J. Phys. 19:299 | not held | not listed |
| Vinti 1964, IAU Symp. 25 | not held | not listed |
| Physics refs (Stark 1914, Epstein 1916, Born 1960, Banks & Leopold 1978, Froman 2008, etc.) | not held | not listed; out of scope |

Held status from `ls cyclers_pdf/papers | grep -i <author>` and `grep -i <author> CORPUS_INDEX.md`; every
hit's file name was read to confirm it is not the cited work (the corpus `lantoine-russell-2011` file is the
halo-to-halo JAS paper). Proposed new wanted-list row: Lantoine & Russell 2009 ISSFD (and, from Hatten,
Biscani & Izzo 2014).

*Filed as `cyclers_pdf/papers/lantoine-russell-2011-complete-closed-form-solutions-stark-problem-cmda-109-333-doi-10.1007-s10569-010-9331-1.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*
