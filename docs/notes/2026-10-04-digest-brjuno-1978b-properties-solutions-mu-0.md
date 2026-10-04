# Digest: Brjuno 1978 part III, "Researches on the restricted three-body problem III. Properties of the solutions for mu = 0"

A. D. Brjuno (Bruno), Celestial Mechanics 18:51-101 (1978), DOI 10.1007/BF01233090. Institute of Applied Mathematics, Moscow. Received 28 September 1977. English translation by M. Henon (April 1977) of Preprint No. 25 (1973), Inst. Appl. Math. This is the "Bruno 1973" of the Bruno 1981 reference list. Filed in the private paper corpus as
`brjuno-1978b-researches-restricted-three-body-problem-III-properties-solutions-mu-0-celest-mech-18-51-doi-10.1007-BF01233090.pdf`
(51 PDF pages; PDF page n is journal page 50 + n; a text layer exists but loses most equations and all figure labels, so tables and formulas were read on the page images: Tables I and II, Theorems 2.2 and 2.3, eqs. 2.32 to 2.34 and 3.1).
Evidence tags: READ (page) = read on the page image (journal page numbers); COMPUTED = my own arithmetic on 2026-10-05 (double or 40-digit, scratch scripts not kept in the repository); INFERRED = my reading across sources.
Companions: `docs/notes/2026-10-04-digest-brjuno-1978-periodic-solutions-arcs-mu-0.md` (part II, the paper this continues; "II" below), `docs/notes/2026-10-04-digest-henon-1968-consecutive-collision-orbits.md`, `docs/notes/2026-10-04-digest-bruno-1981-periodic-flybys-of-the-moon.md`, `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`, `docs/notes/2026-10-04-digest-henon-2001-generating-families-II-part-a-type-1.md`, `docs/notes/2026-10-04-digest-henon-2001-generating-families-II-part-b-type-2.md`, `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md`.

## 0. What this paper is, and short answers to the questions put

Part III continues part II (families S of symmetric arcs and T_N of asymmetric arcs at mu = 0, their characteristics in the global section). It adds: (a) a gallery of 109 numbered synodic orbits of the families A0 to C35 with a table of their elements (Table I, Figures 1 to 11); (b) a new exact equation of the characteristics (2.13) and, from it, Theorem 2.1 to 2.5: in the variables (N^-1 = a^(3/2), e*) every plane characteristic is a cosine or sine curve (Theorem 2.2), in the variables (x, y) of (2.30) every characteristic of B_k and C_jk is a straight line (2.31); (c) the values of the Jacobi constant on the characteristics (section 3), with Theorem 3.1 (no two non-conjugate "special" points share a Jacobi constant except at N = 1, proved with an elliptic curve) and Theorem 3.2 (where the extremal points of C along a family are); (d) four unsolved arithmetic problems and their conditional answer (Schanuel's hypothesis). It contains no orbit at mu > 0 and no stability computation: "Sections 2 and 3 aim at preparing the ground for the construction of natural families of generating periodic solutions" (p.100); the continuation to small mu is announced for "our next paper" and is not here.

Answers to the questions put:
1. **W question (Bruno 1981 Table IV, sidereal versus synodic speed): not settled here.** W does not occur in this paper. The only related formula is eq. 3.1, C = -2 H0 = 2 eps' sqrt(a (1 - e^2)) + 1/a, and Table II's "psi(x)" which is unrelated to W. The answer found from part II (the printed column is the e = 1 form of the denominator with the true V) stands unchanged.
2. **Type III expansion (1/2 versus 1/12): not settled by a printed formula here, but this paper's own machinery supports 1/12 with a check.** The expansion of part II (p.44) is not repeated. The one similar formula here, eq. (2.33), a - 1 = -(1/4) e^4 - (5/16) e^6 + O(e^8), is for a different curve (the curve f of eq. 2.16', a < 1, eps' = -1) and reproduces to three digits from Table II (section 3.3). It does not contradict 1/12 and does not confirm 1/2. The independent check of 1/12 is my own solution of (3.11), recorded in the part II digest, and the extremal-point equation of section 3.C reproduced the Bruno 1981 maximum-C orbit (section 5), which tests the same machinery.
3. **What it adds for the #899 enumerator:** an analytic, family-indexed description of all arc families by one-dimensional curves, the exact intersection (bifurcation) structure, exact Jacobi-constant bounds, the locations of the extremal points, and the "sum of simpler orbits" structure of arcs near type I points (section 2). No new arc type beyond part II's S and T_N, and nothing about continuation to mu > 0.

## 1. Definitions and equations (READ pp.51-101)

Notation as in part II: families S (Henon's A_j, B_k, C_jk, symmetric arcs) and T_N (asymmetric arcs on a resonant two-body orbit with N = a^(-3/2) = (p + q)/p rational); section Gamma; symmetry plane pi with coordinates (a~, e~) where a~ is a signed and e~ = eps'(1 + eps'' e) (II, eq. 3.5); domains omega_1 to omega_4. Brjuno's eps is Henon's eps times eps'' (II p.36).

### 1.1 The characteristic equation (section 2.A, eqs. 2.1 to 2.13, READ pp.70-71)

The collision system (2.1) is the system (3.4) of II with t0 = theta: cos(tau - theta) = a (cos eta - eps'' e), sin(tau - theta) = a eps' sqrt(1 - e^2) sin eta, tau = N^(-1) (eta - eps'' e sin eta), and (2.2) 1 = a (1 - eps'' e cos eta). Instead of eliminating a, e as Henon did (leading to II eq. 3.11 in (tau, eta)), Brjuno eliminates eta and tau, giving the characteristic as a relation among a, e, theta only:
- cos eta = (1 - 1/a)/(eps'' e) (2.3); eta = k pi + eta', |eta'| <= pi/2, kappa = sgn eta' (2.4); cos eta = (-1)^k cos eta', sin eta = (-1)^k sin eta' = (-1)^k kappa sqrt(1 - (a - 1)^2/(a^2 e^2)) (2.5); (-1)^k = eps'' sgn(a - 1) (2.6).
- X = (a - 1 - a e^2)/e, Y = sqrt(1 - e^2) sqrt(a^2 e^2 - (a - 1)^2)/e (2.7; the scan prints the radicand with a stray 4, read here from the algebra, COMPUTED check X^2 + Y^2 = 1 identically), U = (a - 1)/(a e), V = sqrt(a^2 e^2 - (a - 1)^2)/(a e) (2.9), W = -i sqrt(a) sqrt(a^2 e^2 - (a - 1)^2) (2.11, my reading; the scan shows "-i" and the radical with a stray 4).
- Q1 = {sgn(a - 1)}^|N^-1 - 1| (X + iY)^-1 (U + iV)^(N^-1) exp W (2.12), all factors of modulus 1 (X^2 + Y^2 = 1, U^2 + V^2 = 1, W imaginary), and the characteristic equation (2.13) exp{i theta - i pi k (N^-1 - 1)} = {Q1}^(kappa sgn(a - 1)) = Q (Q depends also on eps' and k).
- Which sign conventions make (2.12) true. The printed text is garbled, so I tested the following reconstruction (COMPUTED): arg Q1 = arg(X + iY) + N^-1 arg(U - i eps' V) + eps' sqrt(a) sqrt(a^2 e^2 - (a - 1)^2) for a > 1; and arg Q1 = arg(X - i eps' Y) + N^-1 arg(U + iV) - sqrt(a) sqrt(a^2 e^2 - (a - 1)^2) + pi |N^-1 - 1| for a < 1, with e* = Re sqrt(Q1) (sign of e* aside). They agree with each other for eps' = -1 and reproduce every test of section 3 below. The printed (2.12) formula for arg Q1, p.85, reads arg Q1 = arg(X - iY) + N^-1 arg(U + iV) - iW + pi|N^-1 - 1|, i.e. the eps' = -1 case with the signs as I found them.

### 1.2 The coordinate e* and Theorems 2.1 to 2.5

- Domain: with e' = eps'(1 - e) the domain delta' = {a > 0, |e'| < 1} is covered four times by delta; the four omega_i are "packed" into omega' (the region a(1 - e) < 1 < a(1 + e); READ pp.71-72, Figure 13). The new coordinates are N^-1 = a^(3/2) and e* = Re sqrt(Q1) (2.15); for a < 1 the map is two-valued and is made single-valued with e* = Re sqrt(Q1) sgn P (2.17), where P = a e^2 + a - 1 + eps' N^-1 sqrt(1 - e^2)(1 - a + a e^2) (2.16, as I read it) vanishes exactly on the curve f: P = 0 (2.16'), on which e' = -1 and a < 1 (the point of the plane where the transformation is singular; Table II).
- **Theorem 2.1 (READ p.77, eq. 2.19).** In the variables N^-1, theta, e* the characteristics of the families S and T_N are e* = +- cos (1/2){theta - pi (N^-1 - 1) k}, eps'' sgn(a - 1) = (-1)^k, with k >= 0 an integer.
- **Theorem 2.2 (READ p.78).** In the variables N^-1, e* the plane characteristics of S (theta = 0 or pi) are: e* = (-1)^j sin{(2[j/2] + 1)(pi/2)(N^-1 - 1)}, 1 <= N^-1, for A_j in omega_3; e* = (-1)^(j-1) sin{[(j + 1)/2] pi (N^-1 - 1)}, 1 <= N^-1, for A_j in omega_4; e* = +- cos{k (pi/2)(N^-1 - 1)}, 1 - 1/k < N^-1, for B_k; e* = +- cos{k (pi/2)(N^-1 - 1)}, (j - 1)/k < N^-1 < (j + 1)/k for C_jk in omega_1 and omega_2; e* = +- sin{k (pi/2)(N^-1 - 1)} on the same interval for C_jk in omega_3 and omega_4 ([x] is the integer part). The bounds for B_k and C_jk are not strict (they are the zeros of the cosine or sine; strict bounds are given by the curves f*). The characteristic of A_0 in omega_4 is the half axis e* = 0, a > 1 (k = 0 only on it).
- **Theorem 2.3 (READ p.78).** If k is not 0 on a characteristic then sgn(d e*/d N^-1) = -kappa sgn(a - 1): the sign of eta' (the half-integer side of eta modulo pi) is the slope of the characteristic in omega*.
- **Theorem 2.4 (READ pp.79-81).** A point (N0^-1, e0*) of omega* on a characteristic of an S family corresponds to four points of pi, one in each omega_i. If N0 is irrational only one of the four lies on a characteristic of S and the family is unique; if N0 = (p + q)/p is rational, exactly two conjugate points lie on S characteristics and both are intersection points of S characteristics; the integer indices obey k1 + k2 = r (p + q) (2.25), and the parity rules (2.27) to (2.29) (q even: l1 = l2 mod 2; p + q even: k1 = k2 mod 2; q and p + q odd: l1 + l2 = k1 + k2 mod 2).
- **Theorem 2.5 (READ pp.81-82).** At an intersection point of the first type (type I of II) the number e is transcendental, proved by Lindemann's theorem (exp W algebraic with W nonzero algebraic is impossible). Consequence: on the curve f there are no intersection points.
- Ordering of local characteristics at a type I or II point: the family with the lowest k has the lowest eta and the smaller slope |d e*/d N^-1|; for k1 = k2 either theta1 not equal theta2 (the eta have opposite signs) or e* = 0 (then the lowest eta belongs to the one with negative n, which by Theorem 2.3 has positive slope) (p.81).
- a is algebraic at every intersection point since a = N^(-2/3) with N rational; points of the second type lie on e* = +- 1 with e = |a - 1|/a (algebraic) (p.81).

### 1.3 a < 1: the straight-line coordinates (section 2.E, eqs. 2.30 to 2.34, READ pp.82-86)

For a < 1 and e* > 0: x = 1/(N - 1) with N = a^(-3/2) the mean motion (verified numerically on Table II) and y = 2 arccos e*/(pi (1 - N^-1)) (2.30). The triangle FGH (Figure 14) goes to the infinite triangle F'G'H': x >= z0 = 1/(2 sqrt 2 - 1), 0 <= y <= x - z0 + phi(x), phi(x) = 0 for x <= z0, phi < 0.2 for moderate x. **The characteristics are straight lines (2.31): y = k for B_k, y = (k - j)|x - j/(k - j)| for C_jk** (two segments meeting on the x axis at x = j/(k - j)); so all the complexity sits in the boundary function phi(x), given in Table II for x < 5000 and by (2.32): phi(x) = 1/2 + z0 + L1 x^(1/4) + L2 x^(-1/4) + O(x^(-3/4)), z0 = 1/(2 sqrt 2 - 1) = 0.5469181606..., L1 = -(16/(9 pi)) (3/8)^(1/4) = -0.4428283507..., L2 = -(39/(45 pi)) (8/3)^(1/4) = -0.3525286384.... Derivation (pp.83-85): the curve f has the parametric form a = (c^2 + 1)/(c^3 - c + 2), 1 - e^2 = c^2/a with c = eps' sqrt(a (1 - e^2)) the area integral, 0 <= c <= 1 (as printed; on f the sign is eps' = -1); near a = 1, e = 0: a - 1 = -(1/4) e^4 - (5/16) e^6 + O(e^8) (2.33); and the curve F'H' has y = 1/(N^-1 - 1)... as printed y = x + 1/2 - (16/(9 pi)) e^(-1) - (14/(45 pi)) e + O(e^3) (2.34, with x as above), e = (8/3)^(1/4) x^(-1/4) (1 - (5/6) (3/8)^(1/2) x^(-1/2)) + O(x^(-5/4)). For x < 50 one may take phi = 0 and the boundary is y = x - z0.
- Example 1 (p.85): for N = 2, x = 1 and y <= 0.7, so only one intersection y = 0; in pi that is four points, two in omega_3 (type II) and two in omega_1 (type II). For N = 3/2, x = 2, y < 1.7: two intersections, y = 0 and y = 1; the first gives two pairs of conjugate type II points, the second two pairs of type I.
- Example 2 (pp.85-86): C12 (y = |x - 1|) and C56 (y = |x - 5|) intersect at x = 3, y = 2, where B_2 (y = 2) also passes; N = 1 + 1/x = 4/3.
- Intersection counts (p.82): for given rational N < 1 there are exactly 2(p + q) + 1 intersection points in omega', hence 4(p + q) + 2 in the four omega_i, two of which lie on the curve P2**.

### 1.4 Spatial characteristics of T_N (section 2.F, READ p.86)

For rational N fixed, the characteristic of T_N in the variables (theta, e*) is (2.19) with theta free: e* = +- cos (1/2){theta - pi (N^-1 - 1) k}. For a > 1 it has two curves with theta from -pi to pi; for a < 1 the band |e*| < f*(a) is excluded and theta runs over [-2 arccos f*(a), 2 arccos f*(a)] on two curves whose end points are equivalent. Characteristics of different k are displaced in theta by k pi q/(p + q); intersections lie on the vertical lines theta = k pi/(2(p + q)).

## 2. The orbit gallery (section 1, Table I, Figures 1 to 11, READ pp.51-65)

Synodic orbits of the families with consecutive collisions at mu = 0 for A0, A1, A2, B1, B2, B3, C12, C23, C24, C34, C35 (Figures 1 to 11); their points in pi are on Figure 12 (every fourth line of Table I). Orbits are on two scales (a <= 1 and a > 1). The text gives these structural facts, each checkable against the table:
- Near a type I intersection point an arc is "close to a sum" of simpler arcs: orbit 9 (A0) is an A0 orbit like orbit 1 taken twice and one A1 orbit like 16; orbit 25 (A1) is orbit 1 taken three times plus a B1 orbit like 42 taken twice; orbit 29 (A2) is one orbit 1 plus orbit 42 twice; orbit 32 (A2, the vertical part of its characteristic) is five orbits like 1; orbit 37 is one 3-like and two 46-like; orbit 47 (B1) is orbit 1 twice plus orbit 42 once; orbit 52 (B2) is orbit 1 twice plus a 63-like; orbit 60 is two orbits 1 and one 55; orbit 62 is orbit 18 twice; orbit 68 is two 3-like plus one 46-like; orbit 74 is two 1-like plus one 69-like; orbit 78 is two 1-like plus three 42-like; orbit 80 is two 1-like plus one 3-like; orbit 96 is two 45-like plus one 79-like; orbit 107 is two 79-like plus one 45-like, completing orbit 96 to a two-body orbit run three times (N = 3, as I read "described three times"). This is the same structure as II (3.12): the arcs of S near a type I point tend to (L'L'')^l L' (L''L')^l.
- Orbit 2 (a~ = e~ = -1; a = 1, e = 0, the circular orbit with the Moon's orbit) is the "triple orbit" between orbits 18 and 20; orbits run on a circle described once to six times (orbit 2 described twice, four, five, six times in B1, B2, A2, B3) appear in the vertical parts of the characteristics.
- Collision orbits (a change of sign of e~) separate consecutive orbits within a family where the text says so (for example 3 and 5, 9 after the collision of A0 with P1, 46 and 47, 69 after 68).
- Orbits 77 and 66 have a~ = 1.40572, N = 3/5 (COMPUTED: 1.40572^1.5 = 1.66666); orbits 104, 109 have a~ = 0.71138 and N^-1 = 0.6000 (N = 5/3).

### 2.1 Table I (READ pp.52-53; 97 rows)

Columns: orbit number, family, a~, e~ (as in the part II digest: e~ = eps'(1 + eps'' e), so |e~| > 1 for an orbit at pericentre and < 1 at apocentre; a~ = eps a), tau/pi. The scan is sharp; every digit below was read on the image. The printed orbit numbers skip some (4, 12, 15, 17, 19, 24, 31, 33, 39, 40, 43, 50, 76: "not represented"); "99*" is printed as such.

```
orbit  family  a~         e~         tau/pi
1      A0       -1.37126   -1.87614   0.21640
2      A0       -1.00000   -1.00000   0.50000
3      A0       -1.06131   -0.04349   0.90000
5      A0       -1.58265    0.62882   1.94000
6      A0       -1.58740    0.62996   2.00000
7      A0       -1.59392    0.62198   2.08000
8      A0       -1.94110    0.02417   2.88000
9      A0       -2.06961   -0.03887   3.16000
10     A0       -2.40018   -0.39725   3.84000
11     A0       -2.51984   -0.39685   4.00000
13     A1       -2.78148    0.03366   4.82000
14     A1       -2.51984    0.39685   4.00000
16     A1       -2.09757   -0.02948   2.86000
18     A1       -1.58740   -0.62996   2.00000
20     A1       -1.01085   -1.77607   1.78000
21     A1       -1.22983    1.54507   2.42000
22     A1       -1.31037    1.23686   3.00000
23     A1       -1.38149    1.57281   3.54000
25     A1       -1.62785   -1.70211   4.42000
26     A1       -1.84202   -1.45712   5.00000
27     A2       -1.84202    1.45712   5.00000
28     A2       -1.75038    1.76829   4.38000
29     A2       -1.59546   -1.94443   3.84000
30     A2       -1.43400   -1.35184   3.24000
32     A2       -1.01887   -1.17035   2.48000
34     A2       -1.11216    0.25002   3.26000
35     A2       -1.21141    0.82548   4.00000
36     A2       -1.21879    0.77252   4.26000
37     A2       -1.33256    0.09164   4.82000
38     A2       -1.39895   -0.06771   5.16000
41     B1        2.08008    0.48075   3.00000
42     B1        1.95711    0.28446   2.48000
44     B1        0.72634   -1.54399   0.85511
45     B1        0.85656    1.54216   1.24000
46     B1        1.29462    0.21700   1.72000
47     B1        1.57091   -0.05684   2.16000
48     B1        1.82052   -0.45988   2.68000
49     B1        2.08008   -0.48075   3.00000
51     B2        1.95339   -1.57480   5.68000
52     B2        1.83933   -1.94297   5.18000
53     B2        1.66594    1.69049   4.56000
54     B2        1.51172    1.64708   3.44000
55     B2        1.32171   -1.92590   2.84000
56     B2        1.20558   -1.27253   2.36000
57     B2        0.81355   -0.57229   1.78677
58     B2        0.88772    0.31730   2.22000
59     B2        1.15655    1.76811   2.74000
60     B2        1.29870   -1.92374   3.16000
61     B2        1.38827   -1.46565   0.56000
62     B2        1.65419   -1.41227   4.14000
63     B2        1.84468   -1.94319   4.82000
64     B2        1.99766    1.78826   5.40000
65     B3        1.50543    0.12991   5.76000
66     B3        1.40572    0.71138   5.00000
67     B3        1.38301    0.57374   4.60000
68     B3        1.28840    0.09683   4.18000
69     B3        1.22289   -0.06513   3.86000
70     B3        1.12893   -0.82437   3.32000
71     B3        0.86207   -1.33211   2.74754
72     B3        0.92722    1.60949   3.26000
73     B3        1.12977    0.12185   3.82000
74     B3        1.20236   -0.08686   4.16000
75     B3        1.25279   -0.57516   4.52000
77     B3        1.40572   -0.71138   5.00000
78     B3        1.59279   -0.05570   5.84000
79     C12      -0.76848    0.09495   0.90027
80     C12      -0.94747   -0.10710   1.14000
81     C12      -0.98748   -0.68029   1.38453
82     C12      -0.62996   -0.41260   1.00000
83     C12      -0.55970   -0.10229   0.96210
84     C12      -0.60100    0.24240   1.06052
85     C12      -0.62996    0.41260   1.00000
86     C23      -0.75484    1.37062   2.10967
87     C23      -0.66966   -1.75109   1.89181
88     C23      -0.76314   -1.31037   2.00000
89     C23      -0.87804   -1.20291   2.21562
90     C23      -0.80487    1.55520   1.79428
91     C24       0.61917    0.32024   2.06194
92     C24       0.58915   -0.19188   1.94389
93     C24       0.72013   -0.30183   2.14089
94     C24       0.65945    0.22368   1.89955
95     C34      -0.80197    0.51279   3.20632
96     C34      -0.76743    0.14142   3.12000
97     C34      -0.73847   -0.36418   2.84303
98     C34      -0.79479   -0.72852   2.92404
99     C34      -0.99980   -0.81364   3.44000
99*    C34      -0.98357   -0.09859   3.14000
100    C34      -0.90219    0.19958   2.82000
101    C34      -0.83212    0.73816   2.82894
102    C35       0.69727    1.59295   3.12490
103    C35       0.65763   -1.76053   2.89925
104    C35       0.71138   -1.40572   3.00000
105    C35       0.74046   -1.37698   3.08124
106    C35       0.80633   -1.77810   3.16000
107    C35       0.75706    1.85110   2.88000
108    C35       0.72905    1.58305   2.84915
109    C35       0.71138    1.40572   3.00000
```

COMPUTED check of Table I (every row, from the printed equations): from a~ and e~ I took eps = sgn a~, eps' = sgn e~, eps'' = + if |e~| > 1 and - otherwise, a = |a~|, e = ||e~| - 1|, then eta from (2.3) and tau = a^(3/2) (eta - eps'' e sin eta) for every branch eta = k pi +- eta'. Row 2 (a = 1, e = 0) is the circular orbit at tau/pi = 0.5, the same as Henon's row. Of the other 96 rows, 82 reproduce the printed tau/pi to 3e-4; 12 more are the tangent rows of type II (tau/pi an integer, |cos eta| = 1 where the root is ill conditioned: orbits 6, 11, 14, 18, 22, 26, 27, 35, 41, 49, 104, 109; differences up to 2.5e-3, which is the rounding of the 5-digit a~ at the tangency), and two rows are misprints:
- **Orbit 53 (B2): a~ is printed 1.66594; the arc with the printed e~ = 1.69049 and tau/pi = 4.56 has a = 1.66394.** COMPUTED: solving the timing equation (II eq. 3.11) at tau/pi = 4.56 with the B2 sign product gives eta/pi = 2.30388, a = 1.663941, e = 0.690489; the printed e~ = 1.69049 matches that e, and the printed a~ fails the closure of eqs. (3.9) by 0.03 in cos tau. One digit (3 for 5) is wrong in a~.
- **Orbit 61 (B2): tau/pi is printed 0.56000; the arc with a~ = 1.38827, e~ = -1.46565 has tau/pi = 3.56002** (the unique sign-consistent root; Henon 1968 Table 6 passes through a = 1.40272, e = 0.43687 at tau/pi = 3.6, and a = 1.36903, e = 0.51919 at 3.5, bracketing it). The leading digit 3 was dropped in the print.
Use the recomputed values in a test and keep the printed ones as strict expected-failure records.

### 2.2 Cross-check of Table I against Henon 1968 and Bruno 1981

Orbit 1 (A0): a~ = -1.37126, e~ = -1.87614, tau/pi = 0.21640 is the Bruno 1981 Table IV row tau/pi = 0.2164 (a = 1.37126, e = 0.87614, eps' = -1, eps'' = +1 gives e~ = -(1 + 0.87614)). Orbit 6 (A0, tau/pi = 2) is Henon Table 2 row tau/pi = 2.0 (a = 1.58740, e = 0.37004, eps' = +1, eps'' = -1: e~ = 1 - 0.37004 = 0.62996 as printed). Orbit 14 is the Table 3 / Table 4 end (a = 4^(2/3) = 2.51984, e = 0.60315, e~ = +- 0.39685 = 1 - 0.60315). The remaining rows are not in any of Henon's printed tables (they are at the unprinted tau values, 0.9, 1.94, 2.88 and so on) and so are new sourced values, closed by the equations as above.

## 3. Tables and numbers beyond Table I

### 3.1 Table II (READ p.77; 28 rows): the curve f, P = 0

Columns a, e, N^-1, e*, x, phi(x), for the curve f of eq. (2.16'), listed for 28 values of the area integral c.

```
a        e        N^-1     e*       x            phi(x)
0.50952  0.99880  0.36370  0.99492  0.57158      0.07623
0.52058  0.99528  0.37561  0.98824  0.60155      0.10189
0.53320  0.98961  0.38934  0.97965  0.63758      0.12004
0.54736  0.98193  0.40496  0.96879  0.68057      0.13433
0.56307  0.97243  0.42252  0.95531  0.73166      0.14607
0.58029  0.96125  0.44205  0.93880  0.79228      0.15588
0.59898  0.94857  0.46358  0.91886  0.86420      0.16408
0.61908  0.93454  0.48710  0.89511  0.94969      0.17083
0.64048  0.91928  0.51257  0.86717  1.05159      0.17619
0.66307  0.90291  0.53994  0.83476  1.17361      0.18018
0.68671  0.88552  0.56906  0.79767  1.32052      0.18275
0.71121  0.86716  0.59978  0.75584  1.49865      0.18381
0.73635  0.84785  0.63187  0.70935  1.71643      0.18320
0.76190  0.82757  0.66503  0.65851  1.98538      0.18073
0.78757  0.80624  0.69893  0.60382  2.32145      0.17613
0.81306  0.78377  0.73314  0.54599  2.74730      0.16907
0.83807  0.75998  0.76722  0.48596  3.29594      0.15908
0.86226  0.73464  0.80068  0.42482  4.01695      0.14555
0.88530  0.70745  0.83298  0.36381  4.98741      0.12763
0.90688  0.67800  0.86362  0.30418  6.33249      0.10413
0.92669  0.64579  0.89208  0.24721  8.26618      0.07333
0.94448  0.61012  0.91789  0.19408  11.17874     0.03259
0.96002  0.57008  0.94063  0.14582  15.84376     -0.02237
0.97312  0.52432  0.95996  0.10333  23.97472     -0.09895
0.98368  0.47082  0.97562  0.06730  40.01225     -0.21153
0.99161  0.40607  0.98744  0.03826  78.60926     -0.39278
0.99690  0.32280  0.99536  0.01669  214.33000    -0.74087
0.99960  0.19802  0.99939  0.00321  1650.00170   -1.83099
```

COMPUTED check (all 28 rows): (i) P(a, e) = a e^2 + a - 1 - a^(3/2) sqrt(1 - e^2)(1 - a + a e^2) is below 3e-5 on every row (that is eq. 2.16' with eps' = -1); (ii) N^-1 = a^(3/2) to the printed 5 digits; (iii) x = 1/(N - 1) with N = a^(-3/2) to 1e-5 relative except the last rows where the print of N^-1 (5 digits) limits it (x = 1650.0017 printed against 1665.83 from the rounded N^-1; the printed x is the better value); (iv) phi(x) = y - (x - z0) with y = 2 arccos(e*)/(pi (1 - N^-1)) reproduces the printed phi to 3e-5 on every row except the last two (rounding again); (v) the parametric form a = (c^2 + 1)/(c^3 - c + 2), 1 - e^2 = c^2/a with c^2 = a (1 - e^2) reproduces a to 1e-6 on the rows tested (0.99960, 0.50952, 0.71121, 0.90688); (vi) (2.32) with z0, L1, L2 gives phi(1650.0017) = -1.83072 against -1.83099 printed (error 3e-4, consistent with the O(x^(-3/4)) remainder), and the constants z0 = 0.5469181607, L1 = -0.4428283507, L2 = -0.3525286384 reproduce to the digits printed; (vii) (2.33) at a = 0.99960, e = 0.19802 gives -0.000403 against -0.000400 (a - 1 from the printed a). The maximum of phi is 0.18381 at x = 1.49865 (a = 0.71121). The table is a sourced control for the e = 1 region (a < 1) of the characteristics.

### 3.2 The e* formula and Theorem 2.2 reproduce Henon's tables (COMPUTED)

With the sign reconstruction of section 1.1 I computed |e*| = |Re sqrt(Q1)| from (a, e, eps') for every row of Henon 1968 Tables 2 to 9 that has a determinate eps' and compared with the right side of Theorem 2.2 for the row's family and N^-1 = a^(3/2): all 202 rows with a > 1 (A0, A1, A2, B1, B2) agree to 3e-3, and all 74 rows with a < 1 (B1, B2, C12, C23, C24) agree. (Eleven further a > 1 rows, at tau/pi = n + 1/2, were excluded because eps' is not fixed by my sign solver there.) This is an independent confirmation that the paper's equation (2.13) and Theorem 2.2 describe the Henon families, and that the e* transformation can be used as a membership test: given (a, e) it identifies the family and the index k. It does not test Theorem 2.2's stated N^-1 intervals (the paper says they are not strict).

### 3.3 Section 3: Jacobi constant (READ pp.86-100)

- C = -2 H0 = 2 eps' sqrt(a (1 - e^2)) + 1/a (3.1), depends only on a, e, eps'. On omega': -2 sqrt 2 < C < 3; C = 3 only at a = 1, e' = 1 (type IV*, a saddle point of the level lines); C = -1 at a = 1, e' = -1 (type III); on the upper boundary of omega' C tends to 2 sqrt 2 and on the lower to -2 sqrt 2 as a tends to infinity. COMPUTED: on the boundary a e = |a - 1| (3.2), a (1 - e^2) = 2 - 1/a, so C = 1/a + 2 eps' sqrt(2 - 1/a), which tends to +- 2 sqrt 2; and at a = 2^(2/3) = 1.58740, eps' = +1 it is 2.97096, equal to Henon Table 2 row tau/pi = 2.0 (2.97093) to 3e-5. Level lines are single-valued e'(a); the level line enters omega' from below for a > 1 if -2 sqrt 2 < C < -1, from below for a < 1 if -1 < C < 2, from above for a < 1 if 2 < C < 3, and leaves only from above for a > 1 if 2 sqrt 2 < C < 3 (p.89).
- Values at special points: type IV*: C = 3; type III: C = -1; type II: C algebraic with -2 sqrt 2 < C < 3, C not -1; type I: C transcendental (Theorem 2.5).
- **Theorem 3.1 (READ pp.89-92).** If two non-conjugate intersection points are such that one of them is of type II, III or IV*, their Jacobi constants differ. The proof for two type II points reduces (with b_i = a_i^(-1/2), so b_i^6 = N_i^2) to N_1^2 = (1 + 2t - t^2)^3 and N_2^2 = (1 - 2t - t^2)^3 (3.3), then by a birational map to the elliptic curve y^2 = x^3 + 4 x (3.10); its rational points are only (0, 0), the point at infinity and (2, +-4) (3.11), found from the tables of Birch and Swinnerton-Dyer (rank 0, torsion subgroup of order 4); so C1 = C2 only for N1 = N2 = 1. COMPUTED: (2, 4) satisfies 16 = 8 + 8 and (0, 0) satisfies 0 = 0.
- **Problem 1** (open; two type I points with rational N and the same C): answered "no" conditional on Schanuel's hypothesis (section 3.D).
- **Theorem 3.2 (READ pp.93-96).** For a characteristic chi_0 of S: (a) to each type II point (a0, e0) one can associate an extremal point of C along the family, a maximum if e0 > 0 and a minimum if e0 < 0; (b) chi_0 A_k has one more maximum near a = 1 (on the segments of type IV); (c) the point a = 1, e' = -1 on chi_0 A_k is a minimum (C = -1); (d) the point a = e' = 1 on chi_0 B_k is a maximum (C = 3). Proof by classifying the segments of the characteristics between boundary points of omega' into six types (Figure 18) and the sign of dC/da at their ends: types I and V have a maximum and a minimum, III and IV a maximum, II and VI possibly none. "Apparently this theorem exhausts all extremal points" (conjecture, unproved); more accurate extremal values are in Henon et al. 1970 (Henon and Guyot).
- Extremal points: (3.15) F = (1/i) ln Q - theta + pi k (N^-1 - 1) = 2 pi s, and the extremal condition (3.16) G = (F_a C_e - F_e C_a)/(...) = 0 with an explicit G* = 2/(3 a e ...) (3.16'; not transcribed: R = 2 eps' sqrt(1 - e^2) a sqrt(a)(a e^2 - a + 1) - 2 a^4 e^4 + a^4 e^2 - 5 a^3 e^2 + a^2 e^2 + a e^2 + a - 1 + a^4 - a^3, read from a garbled text layer, not checked on the image). Every extremal point satisfies (2.13) and (3.16) with the same integers k and kappa. On T_N the extrema of C are exactly the extrema of theta, which are the symmetric points (type II for N not equal 1; types III and IV for N = 1).
- **COMPUTED test of the extremal condition.** The condition is that the gradient of the characteristic function F is parallel to the gradient of C. With my reconstruction of F (k = 0, eps' = -1, a > 1) and C of (3.1), F_a C_e - F_e C_a = -7.8e-6 at (a, e) = (1.41019, 0.88445) against a gradient scale of 3.58, whereas it is -0.50 at (1.2, 0.8) and +0.66 at (1.5, 0.6). (1.41019, 0.88445) is the maximum-C orbit of Bruno 1981 Table IV (C = -0.39913); so the extremal-point condition reproduces that published orbit from the paper's own equation (2e-6 relative, at the rounding of the printed a and e).
- Problems 2 to 4 (open): do extremal points coincide with intersection points; do two extremal points share a C; does an intersection point share a C with an extremal point. All three are answered negatively (so each singular point of the families is isolated in C) conditional on Schanuel's hypothesis, using Lindemann's theorem for the algebraic-independence step: at an extremal point (other than a = 1, e' = +-1) at least one of a, e is transcendental unconditionally, and both are under Schanuel (pp.96-100).

## 4. Consistency with the other digests

- **Henon 1968:** the characteristic structure here is the same set of curves as his Figure 4 / part II Figure 17 (A0, A1, ...; B1, ...; C_jk), and (3.9)-(3.11) are his eqs. 27-30; section 3.2 above shows the e* machinery reproduces his Tables 2 to 9 (276 rows).
- **Bruno 1981:** Table IV's maximum-C orbit is reproduced from the extremal-point equation (section 3.3). The "arrival angle" alpha of part II 4.C is not used here. Table I orbit 1 equals the Bruno row tau/pi = 0.2164.
- **Henon 1997:** the book's S-arc lines Z = beta A - alpha in the (A, Z) plane at fixed C are the analogue of this paper's straight lines y = k, y = |(k - j) x - j| in (x, y); both turn the characteristics into lines with all the complication in one transcendental boundary (here phi(x); in the book the domain boundary D). That the two coordinate systems are related is INFERRED from the digests (the 1997 digest, section 2.1 and relation 4.28), not checked in the book. The book's T-arc existence condition I/J > 2^(-3/2) is the same as this paper's N < 2^(3/2) (II p.32). The book's "critical arcs" (Table 4.4, extrema of C) are the extremal points of Theorem 3.2: the book's critical orbits A0(1) etc. correspond to the type II associated maxima and minima; the book cites Hitzl & Henon 1977 and does not cite this paper (not found in the 1997 digest).
- **Henon 2001 (parts A and B):** the type 2 bifurcation orbit of Part B (a Keplerian ellipse tangent to the Moon circle, a = (I/J)^(2/3), e = |a - 1|/a, eq. 17.23) is exactly the type II point of this paper (|a~| = (p/(p + q))^(2/3), a e = |a - 1| from (3.2)); the Jacobi constant and speed of that orbit given in Part B (eqs. 17.24 to 17.27: C from V_Y = eps' sqrt((2a - 1)/a)) agree with C = 1/a + 2 eps' sqrt(2 - 1/a) of section 3.3 (a(1 - e^2) = 2 - 1/a) (INFERRED from the Part B digest's statement, not re-read in the book). Part A's type 1 and the book's nu-regimes concern mu > 0 and are not in this paper.
- **Part II (Brjuno 1978):** this paper supplies the proofs that part II asserts about the intersection types (Theorems 2.2, 2.4, 2.5) and replaces part II's qualitative Figure 18 by exact curves. The part II type III expansion with the coefficient 1/2 (p.44) is not re-derived here.
- **Hitzl & Henon 1977:** their critical orbits (extremal C) are exactly the extremal points of Theorem 3.2; this paper gives their number and position rule (one near every type II point plus the A_k points near a = 1) but no values.

## 5. Techniques applicable to the project's problems

#890/#895 concern Titania-Oberon; only methods are mapped there.

### 5.1 #899 (second-species continuation to Earth-Moon, Casoliva reproduction, with the Henon enumerator)

What part III adds, concretely:
1. **A one-dimensional enumerator at fixed a.** Instead of solving the (tau, eta) timing equation, take N^-1 = a^(3/2) and find e in (|a - 1|/a, 1) with |e*(a, e, eps')| equal to the Theorem 2.2 value for the family (A_j, B_k, C_jk) and branch; then eta from (2.3), tau from (2.1). My reconstruction of e* (section 1.1) reproduces 276 Henon rows, so the test exists (section 3.2). Build a generator that returns, for a given family label and a (or N^-1), all (e, eps', eps'', eta, tau) and checks them against Henon's tables and Table I.
2. **Exact bifurcation nodes.** The type I to IV* points and their number for rational N (2(p + q) + 1 in omega' for N < 1; Examples 1 and 2) give the junctions of the arc alphabet exactly. For the Earth-Moon resonances of Casoliva (small p, q) the junctions can be listed without a numerical search; the junction Jacobi constants are algebraic for type II (C = 1/a + 2 eps' sqrt(2 - 1/a)), exactly -1 and 3 for III and IV*, and Theorem 3.1 says they are all distinct except at N = 1, so at fixed C below 3 each type II junction is met at most once per family pair and no two coincide.
3. **Extremal points as bifurcation seeds.** Theorem 3.2 locates one extremal point of C next to every type II point, plus the A_k maxima near a = 1; condition (3.16) (gradient of the characteristic function parallel to the gradient of C) reproduced the Bruno 1981 maximum-C orbit (section 3.3). That is a general routine for the critical generating orbits, an independent route to Hitzl & Henon Table I and the 1997 book Table 4.4 (it should reproduce the 1977 digest's recomputed values; not run here).
4. **Composite structure.** Near a type I point an orbit is the sum of simpler arcs (section 2); the enumerator can seed composite chains from the simple ones by the tendency (L'L'')^l L' (L''L')^l of II (3.12), instead of searching.
5. **New arc type?** No. S and T_N are all the arcs at mu = 0 (II section 5). The part III information beyond Henon's enumerator is analytic structure and the Jacobi-constant bounds, not new seeds. In particular, nothing here continues a generating orbit to small mu; for that part III points to the next paper in the series (Bruno 1981's reference list names preprints 103/1978, 51/1980, 148/1980, 149/1980, 25/1981 for it; none held).
6. **Where the Earth-Moon target lies.** With the Jacobi constant of arcs bounded by -2 sqrt 2 < C < 3 (omega'), the Earth-Moon resonant cyclers of Casoliva (Jacobi constants listed in the catalogue rows) can be checked for membership before any continuation: at mu = 0 a Casoliva orbit whose Jacobi constant is outside this interval has no second-species generating arc of this type (the interval is for the mu = 0 arcs; at the Earth-Moon mass the Jacobi constants shift).

### 5.2 #906 (turn gate; do not reject small turns near resonance)

- Type II junctions (tangent ellipses, a e = |a - 1|) have v1 = 0 and so demanded turn exactly 0 (II 4.C), and C = 1/a + 2 eps' sqrt(2 - 1/a): a closed-form flag "resonant tangent junction (p, q) = from a = (p/(p + q))^(2/3)" for the gate to return as "indeterminate". Type IV* (a = 1, e = 1) has C = 3 and V = sqrt(3 - C) = 0: the relative speed vanishes, the arc shrinks to the secondary, the demanded turn is undefined and the first-order periapsis is infinite: another "indeterminate" case. Type III (a = 1, e = 0, C = -1) is the circular orbit coincident with the secondary's orbit, V = 2, turn 0.
- Theorem 3.1 means those special C values are isolated: a gate that sees C within rounding of 3 or -1 or of 1/a + 2 sqrt(2 - 1/a) for a rational resonance can attribute the near-zero turn to the resonance with confidence.

### 5.3 #928 (regularised propagator and transition matrix for close passes)

Little. Controls: the extremal-point equation (3.16) and the e* membership test give analytic relations between a, e and the family that a propagated arc must satisfy at its collision point (a check of the mu -> 0 limit of a regularised pass); and the Table II/Table I rows are sourced (a, e, tau) triples at which an arc's two collisions occur at the secondary. Nothing printed concerns the regularisation itself.

### 5.4 #890/#895 (methods only)

The x, y straight-line coordinates are a method for making a two-parameter matching problem one-dimensional; for two moons the analogous elimination of the anomaly in favour of an angle (eq. 2.13) is a candidate for organising the Titania-Oberon junctions. Not evaluated.

## 6. Recommended follow-ups (task numbers not registered)

1. Implement the e*-based membership test and the one-dimensional arc enumerator (section 5.1) with Table I (96 rows at 3e-4, with orbit 53 and 61 as strict expected failures and tangent rows at 3e-3) and Table II (28 rows, all five checks of section 3.1) as sourced tests; recheck the sign reconstruction of Q1 against the printed (2.12) on the images of pp.71 and 85 (the text layer is garbled there).
2. Add the type II junction closed forms (C = 1/a + 2 eps' sqrt(2 - 1/a), v1 = 0) and the type III/IV* values to the #906 gate flags.
3. Run the extremal-point equation (3.16) over all families and compare with the Hitzl & Henon 1977 Table I recomputation (179 orbits) as an independent check.
4. Transcribe (3.16') (the R polynomial) from the page image and test it; this digest left it unchecked.
5. Obtain the follow-up preprints on generating periodic solutions and their continuation (Bruno 1981 Appendix: 103/1978, 51/1980, 148/1980, 149/1980, 25/1981); none is held, and they hold what part III only prepares.
6. Record the two Table I misprints (orbit 53 a~, orbit 61 tau/pi) in a note beside the data.
