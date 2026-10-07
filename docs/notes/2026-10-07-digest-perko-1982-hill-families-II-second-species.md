# Digest: Perko 1982 II, "Families of Symmetric Periodic Solutions of Hill's Problem II: Second Species Periodic Solutions for C << -1" (#960 batch 36)

Lawrence M. Perko (Northern Arizona University), "Families of Symmetric Periodic Solutions of Hill's Problem II:
Second Species Periodic Solutions for C << -1", American Journal of Mathematics 104(2):353-397 (April 1982),
doi 10.2307/2374163 (Crossref: title, 104(2), first page 353, 1982-04 agree). Received 15 October 1980. The
title on the page image (p.353) matches. Acknowledges M. Henon (Nice).
- Supplied file `bfdfdbd3-perko1982_1.pdf`, 46 pp. (1 JSTOR cover page + journal pp.353-397),
  md5 f01b4c7d2cfa0eb172ced13655097ecf. JSTOR scan with a text layer (equations garbled).
- **Proposed corpus filename:**
  `perko-1982b-families-symmetric-periodic-solutions-hill-problem-II-second-species-C-much-less-minus-1-amer-j-math-104-353-doi-10.2307-2374163.pdf`
- **How I read it:** the whole text from the text layer; on 150 dpi page images: p.358-359 (eqs 15-20, b(t)
  of eq. 19, Theorem 1), p.385-389 (eqs 61-62, the 2nd Species Existence Theorem, the a/c/g' identification,
  Figures 2-3, the conjecture), p.391-392 (type B theorem, family g closest approach), p.394, 396 (composite
  theorem, closing remarks). Checks: `check_perko.py` -> `check_perko.out` (secs A, C, D),
  `check_perko_Q0.py` -> `.out`, `check_perko_composite.py` -> `.out`, `check_perko_g_independent.py` -> `.out`, `check_perko_typeB_B1.py` -> `.out`.
- Wanted list: **row 21**. Combined scope note: `theorem-scope.md`.

## 0. Verdict

**The existence proof, by matched asymptotics, of infinitely many second-species symmetric families of Hill's
problem for C << -1: a, c, g', g'', ... (one per root of 3t1 = 4 tan t1, "type A"), Hill's direct family g
("type B", proof only sketched), and a composite family per root (two angular points). No tables.**
- Every theorem is "there exists a C1 < 0 such that for C < C1"; C1 depends on the root t1 and is never given.
  Symmetric only (two perpendicular axis crossings). Asymmetric orbits are not mentioned anywhere.
- Henon 2003's families: {+-1} = a, c and {+-2} = g'+- and {i, e} = g are covered for C << -1;
  {+-3} = He is Perko's g''+- (covered for C << -1 as a type A family); {+1, -1} = g3 is the composite family at
  t1 = 4.41937 (covered). The other Henon 2003 limits (Ha, Hb, Hc, Hd, Hf, Hg: mixed {+-j} and {i, e} arcs,
  or different roots) are "type AB" / "type AA" / n-point composites that Perko says could be done but
  "the work ... would be monumental" (p.396). Not proved. The Table IV asymmetric groups: not covered.
- Two of Perko's conjectures are tested by Henon 2003 (sec. 3). Neither matches as stated.
- **My checks of the printed numbers (sec. 4):**
  - Constants are right: roots 4.41937, 7.68213; 2 theta1 = 17 deg (p.387) and 163 deg (p.394) (17.157 deg and
    162.843 deg); the integral (67) = 1.298417 (printed 1.298416).
  - Family a leading term k1(1 - cos t1)|C|^(1/2) = 1.331162 |C|^(1/2) matches Henon 1969 Table 7 Q1 = 1.33116.
    The first-order term: my integration gives (abs(xi0) - lead)|C| = 0.2818 at C = -400 (still drifting by 1e-3);
    eq. (61) with my evaluation of B1 (eq. 19) gives 0.2670. A 5 % gap, unresolved (my B1 is a reconstruction).
  - Family g: the x'(0) correction -1/(4C) is right (my (x' - |C|^(1/2))|C| = 0.2500 at C = -80). The y(0)
    correction -(1 + B1)/C also agrees: B1 from eq. (19) with j = (0, -1) (right-handed with i along V1 = (-1, 0))
    gives 1 + B1 = 0.0938; my integration gives 0.0975 at C = -160 and is still falling.
  - **Family g closest approach is wrong as printed.** p.392 gives x(t*, C) = -+[j . r1b'(t1)]^2/(2C^4) with
    j . r1b'(t1) = 1.298416. My integration (positive control: Henon 1969 Table 4 at Gamma = -1 and -3 to all
    printed digits; independent re-integration with a second method) gives a closest approach 2.95 times larger at
    C = -10 and 31.5 times larger at C = -320. The data fit x(t*) = (1.5 ln|C| + c)^2/(2C^4), c -> -1.363
    (c = -1.224, -1.307, -1.342, -1.356, -1.361, -1.363 at C = -10, -20, ..., -320). That is an ln(eps) term,
    eps = |C|^(-3/2), in the miss distance, which the printed formula does not have. Perko quotes Henon (private
    communication) for "excellent agreement" down to C = -10; my numbers do not agree with that.
  - Composite family vs g3 (tracked from Henon 2003's g3 at Gamma = -4): (abs(xi0) - lead)|C| tends to about
    -1.04; eqs (61) with d = (1/2) tan t1 and my B1 give +0.092. Unresolved; see sec. 4.
  - So: the existence theorems stand (they are qualitative), but do not use the printed first-order corrections
    as seeds without a check. Use Henon's tables.
- **What it gives R15 / R16:** R15: background only for the symmetric families; does not cover the asymmetric
  groups. R16: does not cover. Details in `theorem-scope.md` sec. 5.
- **Catalogue implication (PROPOSAL only):** none.

## 1. Content

- **Sec. 1 (pp.353-354).** Claims an infinite number of families a, c, g, g', g'', ...; all except g tend to a
  limit orbit with one angular point at the origin; g tends to a limit orbit "with two cusps at the origin". The
  Henon 1969 identification of a, c, g, g' "allows us to also verify the conjectures made by Matukuma in 1952".
- **Sec. 2 (pp.354-389), type A.** Generating solution (7): x = k1(cos t - cos t1), y = k1(-2 sin t + (3/2) t cos t1),
  with 3t1 = 4 tan t1 (eq. 6), k1 = +-(1 - (3/4) cos^2 t1)^(-1/2) (so the arrival speed V1 = 1). Outer expansion
  (eqs 10-20, with ln(tau) singular parts), inner hyperbola near the origin, matching (eqs 33-35), Theorems 1-4,
  then the theorem (p.386, image):
  "For each t1 satisfying 3t1 = 4 tan t1, there exists a C1 < 0 such that for C < C1 there is a unique
  one-parameter family of 2nd species symmetric periodic solutions, r(t, C), of Hill's equations (1) with
  parameter C and initial conditions x'(0, C) = y(0, C) = 0,
  x(0, C) = k1[|C|^(1/2)(1 - cos t1) + Q0/C + O(ln|C|/|C|^(7/4))],
  y'(0, C) = k1[|C|^(1/2)(-2 + (2/t1) sin t1) + Q1/C + O(ln|C|/|C|^(7/4))]",
  Q0, Q1 from (61), (62) with d = -(1/2) cot t1 and B1 = int_0^t1 j . b(t) dt.
  - Eccentricity of the inner hyperbola e0 = (1 + (1/4) cot^2 t1)^(1/2) + O(eps); "agrees with Henon's equation (24)".
  - Identification (p.387-388): first root, k1 > 0 = family a (from L2 (3^(-1/3), 0)); k1 < 0 = c; second root =
    g', whose two symmetric branches g'+- "meet in a common orbit with C = 4.5" (Henon 1969 Table 10: g'1 = g1 at
    4.49999). Roots n = 2, 3, ...: g''+-, g'''+-, ... Conjecture (p.389): g+(n) and g-(n) meet for every n.
    The p.388 text says "second solution of (19)"; this must mean eq. (6) (eq. 19 is b(t)). A slip.
- **Sec. 3 (pp.390-392), type B, family g.** Generating solution (63) x = sin t, y = 2 cos t + 2 (one ellipse
  through the origin), t1 = pi. Theorem (p.391): "There exists a C1 < 0 such that for all C < C1, there is a unique
  one-parameter family of symmetric periodic solutions r(t, C) of Hill's equation (1) with parameter C and initial
  conditions x(0, C) = y'(0, C) = 0, x'(0, C) = +-[|C|^(1/2) - 1/(4C) + O(ln|C|/|C|^(7/4))],
  y(0, C) = +-[4|C|^(1/2) - (1 + B1)/C + O(ln|C|/|C|^(7/4))]", B1 = int_0^pi j . b(t) dt. Stated "without going
  through the details of the proof". Period T = 4 pi + delta t1 + O(eps^2 ln eps). Closest approach (p.392):
  x(t*, C) = -+[j . r1b'(t1)]^2/(2C^4) + O(ln|C|/|C|^(11/2)), j . r1b'(t1) = 1.298416 from eq. (67).
- **Sec. 4 (pp.392-396), composite.** Two type A arcs with the same t1 and k1 = +-; perpendicular y-axis crossing
  at t* = t1 + O(eps |ln eps|). Theorem (p.394): same form as type A with d = (1/2) tan t1; "symmetric with
  respect to both the x and y axes". Conjecture: the continuation reaches orbits of radius ~1/C circling the origin
  three times retrograde with T ~ 6 pi/C^(3/2) for C >> 1 (and 2n + 1 times for the n-th root). Type AA
  (different t1), type AB (type A + type B arcs) and n-point composites: asserted possible, not done.

## 2. Frame and conventions

As paper I: Hill's equations, C = Henon's Gamma, eps = |C|^(-3/2), x -> x/|C|^(1/2), time not scaled.
Type A families start at a perpendicular x-axis crossing far from the origin (t = 0) and reach the origin
side at t ~ t1 (second perpendicular x-axis crossing; half period). Henon's sign conventions differ by mirror
images (Henon tabulates xi0 < 0 for the far crossing of a/c at Gamma < 0).

## 3. Against Henon 2003 (held; Tables III-V, secs 4.6, 4.8)

- Henon's arc labels: {+j} / {-j} = Perko type A with the j-th root of 3t = 4 tan t (Henon Table III: {+1} = a,
  {-1} = c, {+2}/{-2} = g-/g+, {i, e} = g). I take {+3} = He = Perko's g''+ (third root 10.87356).
- **Meeting conjecture (p.389), n = 2:** Henon 2003 sec. 4.6: He starts at {+3}, rises to Gamma max = -8.615196
  and returns to Gamma -> -infinity ending at {+1, e, i}; its mirror He' joins {-3} to {-1, i, e}. So g''+ does
  not meet g''- (on Henon's numerics). The conjecture fails at n = 2, given my label identification.
- **Composite conjecture (p.394):** the composite {+1, -1} family is Henon's g3. g3 has Gamma max = 3.806201 and
  returns to -infinity at {+1, i, -1, e}; it does not reach C >> 1. But it meets "f described three times" at
  Gamma = -1.411618 and 0.015388, and f-thrice is exactly the orbit Perko describes (radius ~1/C, three retrograde
  turns, T ~ 6 pi/C^(3/2)). So the conjecture holds only if one changes branch at that bifurcation.
- Henon 2003 does not cite Perko 1982 or 1983 (reference list checked).

## 4. Checks (scripts and outputs in this folder)

All DOP853, rtol 1e-12, atol 1e-14, shooting on a perpendicularity residual.
- Constants (sec. A): roots of 3t = 4 tan t: 4.419371, 7.682131, 10.873562 (Henon 1969 Table 7: 4.41937,
  7.68213, 10.87356). k1(1 - cos t1) = 1.331162 (t1 = 4.41937) and 0.838237 (t1 = 7.68213) = Henon 1969 Table 7
  Q1 = 1.33116, 0.83824. Eq. (67) = 1.298417 against printed 1.298416 (last digit: quadrature with the
  1/(2(pi - t)) subtraction; agreement to 1e-6).
- **Family a (sec. D).** Positive control: Henon 1969 Table 8 at Gamma = -4: -2.738494 (printed -2.73849); at
  Gamma = -100: -13.314449 (printed -13.31445), near-origin crossing 0.0001215 (printed 0.0001215; sign differs,
  sign convention not resolved). (abs(xi0) - 1.331162 |C|^(1/2)) |C| = 0.3047, 0.2957, 0.2900, 0.2852,
  0.2827, 0.2818 at C = -4, -9, -16, -36, -100, -400.
  - Perko's value: eq. (61), d = -(1/2) cot t1 = -0.150851, with B1 = int_0^t1 j . R1 Phi_rv(t1, t) f[r0(t)] dt
    (the -i/tau of eq. 19 has no j part) computed by quadrature: B1 = -3.861892, Q0 = -0.258489, so the predicted
    limit is -abs(k1) Q0 (sign fixed by x(0) = k1[... + Q0/C] with C < 0) = +0.2670. A 5 % gap to my 0.2818 limit (the data
    look converged to 1e-3; the error term is O(ln|C|/|C|^(3/4))). The orientation of (i, j) and of R1 in my
    reconstruction is the weakest link. Not called a misprint.
- **Family g (sec. C), closest approach.** Positive control: Henon 1969 Table 4 at Gamma = -1: xi = 0.0210367,
  T/2 = 5.112276 (printed 0.02104, 5.11228); Gamma = -3: 0.00452259, 5.709887 (printed 0.004523, 5.70989).

  | C | closest approach (mine) | Perko 0.842944/C^4 | ratio | sqrt(2 x C^4) - 1.5 ln abs C |
  |---|---|---|---|---|
  | -5 | 1.5110e-3 | 1.3487e-3 | 1.12 | -1.040 |
  | -10 | 2.4868e-4 | 8.4294e-5 | 2.95 | -1.224 |
  | -20 | 3.1728e-5 | 5.2684e-6 | 6.02 | -1.307 |
  | -40 | 3.4309e-6 | 3.2928e-7 | 10.4 | -1.342 |
  | -80 | 3.3228e-7 | 2.0580e-8 | 16.1 | -1.356 |
  | -160 | 2.9821e-8 | 1.2862e-9 | 23.2 | -1.361 |
  | -320 | 2.5340e-9 | 8.0389e-11 | 31.5 | -1.363 |

  - Independent check (`check_perko_g_independent.py`, Radau from the top of the orbit): 2.4865e-4 at C = -10,
    3.425e-6 at C = -40 (the 6-decimal y0 limits the agreement to 1e-3 relative).
  - (x'(0) - |C|^(1/2))|C| = 0.2065, 0.2409, 0.2460, 0.2487, 0.2496, 0.2499, 0.2500, 0.2501 (C = -1, -3, -5, -10, ..., -160): the
    printed -1/(4C) term is right.
  - (y(0) - 4|C|^(1/2))|C| = 0.390, 0.300, 0.241, 0.175, 0.134, 0.112, 0.102, 0.098 (C = -1 ... -160). The theorem
    says this tends to 1 + B1. `check_perko_typeB_B1.py`: B1 = int_0^pi j . R Phi_rv(pi, t) f(r0(t)) dt = -+0.906168
    for j = (0, +-1); with j = (0, -1) (right-handed with i along V1 = r0'(pi) = (-1, 0)) 1 + B1 = 0.093832. The data
    (0.0975 at C = -160, still decreasing) agree. This also checks my B1 machinery used for type A below.
- **Composite / g3 (`check_perko_composite.py`).** Shooting xi0 with a perpendicular second y-axis crossing near the
  origin. Positive control from the same shooter: at Gamma = -4 it returns xi0 = -2.490440 (Henon 2003 Table XIII
  g3: -2.490440), second y-axis crossing at y = -0.965692. Then continuation in Gamma, each step seeded from the
  previous root: (abs(xi0) - 1.331162 |C|^(1/2))|C| = -0.6875, -0.7444, -0.8014, -0.8744, -0.9206, -0.9506,
  -0.9850, -1.0025, -1.0232 at C = -4, -6, -9, -16, -25, -36, -64, -100, -256; T/4 -> 4.4177 (t1 = 4.41937).
  The theorem with d = (1/2) tan t1 and my B1 (-3.861892) gives +0.092. Unresolved: the same B1 gives type A to
  5 %, and would have to be -9.71 to fit here; either eq. (61) does not carry over to the composite with only
  d changed, or my reading of the composite set-up is wrong.

## 5. Citation mining

| Cited work | Held? | Wanted list |
|---|---|---|
| [1] Hill 1878, Amer. J. Math. 1:5-26, 129-147, 245-260 | not held | not listed |
| [2] Henon 1969, A&A 1:223 | held (+ tables.txt) | - |
| [3] Lang 1962 (differentiable manifolds) | not needed | - |
| [4] Perko 1964 PhD thesis, Stanford | not held | row 30 |
| [5] Perko 1974, SIAM J. Appl. Math. 27:200 | held, digested | - |
| [6] Perko 1976, Celest. Mech. 14:395 | held, digested | - |
| [7] Perko 1982 I | supplied this batch | row 21 |
| [8] Poincare 1899, Methodes nouvelles v.3 | not checked | - |
| Matukuma 1952 (cited via Henon 1969 p.234) | not held | row 61 (Matukuma 1930-1957) |

No new row proposed.

*Filed as `cyclers_pdf/papers/perko-1982b-families-symmetric-periodic-solutions-hill-problem-II-second-species-amer-j-math-104-353-doi-10.2307-2374163.pdf`. The Perko check scripts and the combined theorem-scope note are filed beside the 1982 I PDF (`perko-1982a-...-<file name>`).*

*Wanted-list row numbers in this digest are the pre-batch-36 numbering; the list was renumbered in batch 36.*
