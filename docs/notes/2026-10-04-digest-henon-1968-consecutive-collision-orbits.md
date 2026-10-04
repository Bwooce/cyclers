# Digest: Henon 1968, "Sur les orbites interplanetaires qui rencontrent deux fois la Terre" (orbits with consecutive collisions)

Henon, M., Bulletin astronomique 3(3):377-402 (1968), in French, with abstracts in English, German and Russian. DOI 10.3406/bastr.1968.14547 (Persee). Manuscript received 22 December 1967.
Filed in the private paper corpus as
`henon-1968-orbites-interplanetaires-qui-rencontrent-deux-fois-la-terre-bull-astron-3-377-doi-10.3406-bastr.1968.14547-french.pdf`
(28 PDF pages; PDF page n is journal page 374 + n, so the text starts on PDF p.3, journal p.377, and the nine tables are journal pp.394-402).
Evidence tags: READ (page) means read on the page image (all equations, Figure 4 and the table pages were viewed as images; the PDF text layer has lost every Greek letter and sometimes a minus sign, so it was used only for the digits of the tables and each table row was then checked against the equations); COMPUTED means my own arithmetic on 2026-10-04 (double precision, scratch scripts not kept in the repository); INFERRED means my reading across sources. The paper is translated as read; quotations are my translations.
Companion digests: `docs/notes/2026-10-04-digest-hitzl-henon-1977-critical-generating-orbits.md` (the same timing condition, with the critical orbits), `docs/notes/2026-10-04-digest-henon-1997-generating-families.md` (the book that supersedes this paper's notation), and `docs/notes/2026-10-04-digest-bruno-1981-periodic-flybys-of-the-moon.md` (the application to the Earth-Moon and Sun-Jupiter systems, which uses Tables 5 to 9 of this paper).

## 0. What this paper is, in eight lines

1. Question (READ p.377-378): a probe leaves a body M2 (the Earth, on a circular orbit about M1 the Sun) and meets M2 again later. Which conics about M1 do this with no mid-course correction? Approximations: M2 on a circle, M2 has no size and no mass, other bodies ignored. The orbit of M3 is then one two-body arc.
2. Setting (READ p.379): units with the circle of M2 of radius 1 and angular rate 1 (Sun-Earth: length is 1 AU, time is 1 yr / 2 pi, velocity unit is the Earth's 29.79 km/s). The two meetings (P at time -tau and Q at +tau, so the interval is 2 tau) are "collisions" of the two points.
3. Three cases (READ p.379): (1) P = Q, the orbit period is commensurable with the Earth's (the trivial case, Roy 1963); (2) P and Q distinct and the planes different (then P and Q are diametrically opposite, and the solutions follow from case 3 by rotating the plane of M3 about the line PQ); (3) P and Q distinct and coplanar, "which does not seem to have been studied until now", the subject of the paper.
4. Method: symmetry about the x axis, then the three Kepler-motion equations at the collision point; eliminate a and e to get ONE implicit equation F(tau, eta; signs) = 0 in the two numbers tau (half the time between collisions) and eta (half the change of anomaly), solved numerically (eqs. 10 and 30). Hyperbolic orbits form one family (A0 in its hyperbolic part); elliptic orbits form an infinite set of one-parameter families.
5. Families are named A0, A1, A2, ...; B1, B2, ...; C12, C23, C24, ... on a plane of (tau/pi, eta/pi) (Figure 4, p.384). Characteristic elements are tabulated for A0, A1, A2, B1, B2, C12, C23, C24 up to tau/pi = 4 (Tables 1 to 9, pp.394-402).
6. Practical limits (READ pp.391-392, section III): the orbit must not hit M1 (Sun radius 0.004652, Moon-as-primary case Earth radius 0.01659 in units of the Moon's orbit), the departure speed must be technically feasible (5 km/s is 0.1679 of the Earth's speed; 4.887 in the Earth-Moon unit 1.023 km/s), and the duration 2 tau must be acceptable. At the Earth-Moon scale "tau/pi" counts sidereal months of 27.32 d.
7. Section IV (READ p.393): these arcs are the pieces of Poincare's "second species periodic solutions": let the mass of M2 and the distance of approach go to zero together at constant ratio, and the periodic orbit of the restricted problem (Figure 11a, a strong deviation near M2; Broucke 1962 shows many for the Earth-Moon case) becomes a chain of conic arcs joined at angular points (Figure 11b). "Each of these arcs is one of the orbits with consecutive collisions studied here." The paper declines to pursue this ("the study will possibly be taken up elsewhere").
8. What is absent: no stability analysis, no existence proof (the families are found numerically), no positive-mass orbit, no treatment of the two-body arcs that visit M2 more than twice per period (composite orbits come with Perko 1976 and Bruno 1972).

## 1. The definitions and equations, in the paper's order (READ pp.379-390)

Notation. The paper's epsilon, epsilon', epsilon'' are written here as eps, eps', eps''. Compared with the other digests: eps = sigma0 (side of the perihelion/apse S on the x axis), eps' = sigma1 (direct or retrograde), eps'' = sigma2 (apse at S is perihelion or aphelion) of Hitzl & Henon 1977, and the combination s = eps eps'' is their sigma; in Bruno 1981 the symbols are the same as here. (INFERRED from the identical equations, section 5 below.) The "collision point" Q is at t = +tau, P at t = -tau, symmetric about the x axis; the Moon-position (the paper's M2) is (cos t, sin t); the apse S of M3 on the x axis is passed at t = 0.

### 1.1 Switches (eqs. 3, 4, 25; READ pp.380, 383)

- eps = +1 or -1 if the perihelion S has positive or negative abscissa (eq. 3).
- eps' = +1 or -1 if the orbit of M3 is direct or retrograde (eq. 4).
- eps'' = +1 or -1 if the central crossing S is the perihelion or the aphelion (eq. 25).

### 1.2 Hyperbolic orbits (eqs. 5 to 15; READ pp.380-382)

x = eps a (e - cosh F), y = eps eps' a sqrt(e^2 - 1) sinh F, t = a^(3/2) (e sinh F - F), with a > 0 and e > 1 (eq. 5). At F = eta the point is at Q with t = tau (eq. 6):

- cos tau = eps a (e - cosh eta), sin tau = eps eps' a sqrt(e^2 - 1) sinh eta, tau = a^(3/2) (e sinh eta - eta).
- Three equations in four unknowns (tau, a, e, eta), so a simple infinity of solutions. Combining the first two: 1 = a (e cosh eta - 1) (eq. 7).
- Elimination (eq. 8): a = (1 - eps cos tau cosh eta) / sinh^2 eta, e = (cosh eta - eps cos tau) / (1 - eps cos tau cosh eta).
- Timing equation (eq. 10): (1 - eps cos tau cosh eta)^(1/2) [ sinh eta (cosh eta - eps cos tau) - eta (1 - eps cos tau cosh eta) ] - tau sinh^3 eta = 0.
- "Solutions exist only for eps = -1, eps' = -1" (eq. 11). They form one family, A0 (hyperbolic part), Table 1.
- Periapsis abscissa x0 = eps a (e - 1) (eq. 12); Jacobi constant C = 2 eps' sqrt(a (e^2 - 1)) - 1/a (eq. 13); relative speed at the collision V = sqrt(3 - C) (eq. 14); the Jacobi constant is defined "with a view to application to the restricted problem".
- Limit tau -> 0, eta -> infinity (eq. 15): eta = -ln(tau^2/2), a = tau^2, e = 1 + tau^2/2, x0 = -tau^4/2, C = -1/tau^2 (the orbit is rectilinear: M3 goes from R to M1 and back out with infinite speed).
- The first column of every table is tau/pi, "equal to the time between the two collisions, in years" (because the interval is 2 tau and one year is 2 pi time units).

### 1.3 Parabolic orbit (eqs. 16 to 22; READ p.382)

x = eps (p/2)(1 - s^2), y = eps eps' p s, t = p^(3/2) (s/2 + s^3/6) (eq. 16). At s = sigma (the paper's sigma, not Hitzl-Henon's) the point is at Q (eq. 17). Then p = 1 + eps cos tau, sigma^2 = (1 - eps cos tau)/(1 + eps cos tau) (eq. 19), and the third equation becomes tau = (1/3)(2 + eps cos tau) sqrt(1 - eps cos tau) (eq. 20), with the single solution eps = -1, tau/pi = 0.16393 (eq. 21). The other parameters (eq. 22): eps' = -1, x0 = -p/2 = -0.06485, C = -0.72028, V = 1.92880, sigma = 3.7973. COMPUTED: solving eq. 20 gives tau/pi = 0.163926, p = 0.129702, x0 = -0.064851, sigma = 3.79736, C = -2 sqrt(p) = -0.720283, V = 1.928804. All five printed values reproduce.

### 1.4 Elliptic orbits (eqs. 23 to 30; READ pp.382-383)

M3 crosses the x axis 2n+1 times (n >= 0); the central crossing S (number n+1) is perihelion or aphelion. The two cases are written as one (eq. 26):

- x = eps a (eps'' cos E - e), y = eps eps' a sqrt(1 - e^2) eps'' sin E, t = a^(3/2) (E - eps'' e sin E).
- At E = eta > 0 the point is at Q, time tau (eq. 27): cos tau = eps a (eps'' cos eta - e); sin tau = eps eps' a sqrt(1 - e^2) eps'' sin eta; tau = a^(3/2) (eta - eps'' e sin eta). (Same equations at E = -eta, t = -tau for P.)
- First two: 1 = a (1 - eps'' e cos eta) (eq. 28). Then (eq. 29): a = (1 - eps eps'' cos tau cos eta) / sin^2 eta, e = (eps'' cos eta - eps cos tau) / (1 - eps eps'' cos tau cos eta).
- Timing equation (eq. 30): (1 - eps eps'' cos tau cos eta)^(1/2) [ eta (1 - eps eps'' cos tau cos eta) - sin eta (cos eta - eps eps'' cos tau) ] - tau |sin eta|^3 = 0.
- "The eps enter (30) only through the product eps eps''." Each solution of (30) has a twin with eps and eps'' both reversed; only one is kept because e from (29) must be positive. eps' does not appear in (30) and is fixed by the second of (27), whose two sides must have the same sign.
- The solution set in the (tau, eta) plane is Figure 4: "an infinity of families of solutions".
- Columns of the elliptic tables: tau/pi; eta/pi (the path M3 travels on its ellipse between the collisions, in revolutions: eta/pi = 1/2 is half an ellipse, 1 is a whole ellipse); the three signs eps, eps', eps''; a; e; x0 = eps a (eps'' - e) (abscissa of S); x1 = eps a (-eps'' - e) (abscissa of the other intersection with the x axis, which M3 actually passes through if and only if eta/pi >= 1, eq. 31); C = 2 eps' sqrt(a (1 - e^2)) + 1/a (eq. 33); V = sqrt(3 - C).

### 1.5 The A0 family and the asymptotic forms (eqs. 34 and 35; READ pp.384-388)

- A0 continues from the parabolic orbit into ellipses and goes to tau -> infinity. Its cycle of shapes (READ p.386): the eccentricity falls and is zero at tau/pi = 0.5 (M2 and M3 describe two identical half-circles in opposite directions); then e rises and the orbit is rectilinear at tau/pi = 1; the sense reverses and is direct; at tau/pi = 2 (eta/pi = 1) M3 describes exactly one ellipse of period 2 years tangent to the Earth's orbit; at tau/pi = 3 the orbit is rectilinear through the Sun again ("evidently not achievable"); the sense reverses to retrograde; at tau/pi = 4 an ellipse of period 4 years tangent to the Earth's orbit; "then the cycle repeats with a period of 4 in tau/pi" while a grows continuously.
- Asymptotics for tau -> infinity (eq. 34): eta = pi - (pi/tau)^(1/3) sqrt(2) sin(tau/2), e = 1 - (pi/tau)^(2/3) cos^2(tau/2), a = (tau/pi)^(2/3). So eta(tau) is a damped sinusoid of period 4 pi about eta = pi. For the other families (eq. 35): eta = j pi - (j pi/tau)^(1/3) sqrt(2) sin((tau - k pi)/2), e = i - (j pi/tau)^(2/3) cos^2((tau - k pi)/2) [the leading term of e reads as "i" or "1" in the scan; e tends to 1 as tau grows, so 1 is meant (INFERRED, not confirmed on a second copy)], a = (tau/(j pi))^(2/3), with integers j and k: family A_i has (j = i, k = i + 1, except for i = 0) and (j = i + 1, k = i); family B_i has (j = i, k = i) and (j = i, k = i + 2).
- The product eps eps'' is constant along a family, "a general property" (it is explained by S staying on the same side of M1; Broucke 1962). It is -1 for every A family (so x0 < 0), +1 for every B family, and (-1)^(i+j) for the family C_ij (READ pp.387-388).

### 1.6 The families A, B, C in the (tau/pi, eta/pi) plane (READ pp.384-388, Figure 4)

- A0 is special: it joins the hyperbolic orbits and runs from eta = 0 to infinity in tau. The families A_i (i >= 1) each have tau falling from infinity to a minimum and rising again to infinity. The B_i are a second infinite set interleaved with the A_i with the same behaviour. Curves of consecutive A (or B) families are approximately translates of each other along the diagonal eta = tau. Tables 3 to 6 cover A1, A2, B1, B2 for tau/pi <= 4.
- C families: curves in the shape of a figure eight (Figure 4), with a double point at integer coordinates tau/pi = i, eta/pi = j (the label is C_ij); they exist for all i, j with 1 < j/i < 2 sqrt(2) (eq. 36; the printed condition is "1 < j/i < 2 sqrt 2"; the preceding wording "i < j" holds). C12, C23 and C24 are tabulated (Tables 7, 8, 9); they close on themselves (the last row of each table equals the first, and so does the last drawing of Figure 8).
- Figure 4 labels, reading from the scan: for each branch the three signs eps eps' eps'' are printed beside it (e.g. A0 near eta/pi = 0.5 carries "- - +" and "+ - -" on its two sides, matching Table 2).

### 1.7 Special solutions (eqs. 38 to 50; READ pp.390-391)

- (a) eta = tau with eps eps' = 1 satisfies (30) (eq. 38) and gives a = 1, e = 0, eps' = 1 (eq. 39): M3 has the same orbit as M2, coincident at all times. "Evidently without practical interest." These form the dashed diagonal of Figure 4.
- (b) tau/pi = i, eta/pi = j (integers), eps eps'' = (-1)^(i+j), eps' = +-1 (eq. 40): P and Q coincide at R or at the opposite point R', and M3 is an ellipse tangent to the Earth's orbit there, with a = (i/j)^(2/3), e = |1 - (j/i)^(2/3)| (eq. 41) [the scan reads a = (i/j)^(2/3), e = |1 - (j/i)^(2/3)|]. e <= 1 requires j/i <= 2 sqrt(2) (eq. 42). For given tau and eta, eps' = +1 and -1 are the same ellipse traversed in the two senses, so this is a double point of the plane with two branches: if i < j both branches belong to C_ij; if i = j one belongs to B_j and the other is the diagonal eta = tau; if i > j with i - j even both belong to B_j; if i > j with i - j odd they belong to A_(j-1) and A_j. The slopes of the two branches there (second-order expansion of (30)) are d tau / d eta = +- sqrt(2 (i/j)^(2/3) - 1) (eq. 43); the requirement that the radicand be non-negative reproduces (42).
- (c) tau/pi = i + 1/2, i integer (eq. 44): (30) no longer contains the product eps eps'', so both values are allowed and there is again a double point, one branch always in an A family and the other in a B family or on the diagonal. The two solutions are mirror images across the y axis, on which P and Q lie (Figure 10), and rotating the orbit of M3 about PQ gives a family of non-coplanar solutions (case 2 of the introduction).
- (d) Near tau/pi = eta/pi = i + 1/2 (eq. 45) the families A_i and C_(i,i+1) almost touch and have very sharp bends. Put tau = (i + 1/2) pi + X, eta = (i + 1/2) pi + Y (eq. 46); expansion of (30) to second order gives (eq. 47, partly garbled in the scan) 2Y - (1 + eps eps'') X + ... = 0. For eps eps'' = +1 this becomes (Y - X)[...] + ... = 0 (eq. 48) with two solutions: X = Y (the diagonal) and an expression Y = 4/(3 (i + 1/2) pi)-type (eq. 49) representing B_i approximately, "which has an almost horizontal portion near this point". For eps eps'' = -1: two straight lines of slopes 0 and -1 (eq. 50), approximately the families A_i and C_(i,i+1); the true curve is a hyperbola because of the neglected higher-order terms. (Eqs. 47 to 50 are read from the scan with some uncertainty in the coefficients; the qualitative statements are what is relied on.)

### 1.8 Practical limits and the Earth-Moon numbers (eqs. 51 to 56; READ pp.391-392)

- Sun radius in AU: r1 = 0.004652 (eq. 51). Earth radius in lunar distances: r1 = 0.01659 (eq. 52) [Bruno 1981 uses 0.01665 for the same thing]. Exclude every orbit with |x0| < r1 (eq. 53) and, in the elliptic case, those with |x1| < r1 when M3 passes through x1 (eq. 54). This removes "a good part" of the hyperbolic orbits and only the very flat ellipses.
- Departure speed V <= 5 km/s is V <= 0.1679 in the Earth's speed unit 29.79 km/s (eq. 55); that removes all hyperbolic orbits and most elliptic orbits, and the survivors have eps' = +1 (the same sense as the Earth). In the Earth-Moon case the speed unit is 1.023 km/s and V <= 4.887 (eq. 56), "much less severe": it leaves some hyperbolic orbits and all the elliptic ones.
- Duration: 2 tau is tau/pi years in the Sun-Earth case; the tables stop at tau/pi = 4. Orbits exist for every tau, "whereas for commensurable periods the duration is necessarily a whole number of years". For the Earth-Moon case the unit is the sidereal month, 27.32 d, so much larger tau/pi are possible.
- Real orbits need corrections, the most important being the eccentricity of the Earth's orbit.
- Bibliography (READ p.393): Roy 1963 (Astronautica Acta 9:31-46), Kovalevsky 1963, Broucke 1962 (Louvain dissertation, p.75 and figures, Earth-Moon periodic orbits), Poincare 1892-99 vol. III ch. XXXII.

## 2. The printed numbers

All values are transcribed from the table images (journal pp.394-402). The table columns are as in section 1.4. The first column tau/pi is the elapsed time between collisions in years (Sun-Earth) or sidereal months (Earth-Moon). Where the printed value is wrong against the printed equations it is reported as printed and listed in section 2.3. Rows come in a printed order that follows each family as tau first falls, passes its minimum, and rises again (so the same tau/pi appears twice in a table, with different eta/pi); that is not a typographical duplicate.

### 2.1 Table 1, family A0, hyperbolic orbits (READ p.394; 16 rows)

Columns: tau/pi, eta (hyperbolic anomaly, not divided by pi), a, e, x0, C, V. Here eps = -1, eps' = -1.

```
tau/pi   eta       a        e        x0         C            V
0.01000  7.60129   0.00100  1.00049  -0.00000   -999.87824   31.66825
0.02000  6.18664   0.00412  1.00196  -0.00001   -242.59088   15.67134
0.03000  5.33620   0.00968  1.00437  -0.00004   -103.33590   10.31193
0.04000  4.71177   0.01816  1.00767  -0.00014   -55.08889    7.62161
0.05000  4.20759   0.03030  1.01174  -0.00036   -33.05260    6.00438
0.06000  3.77652   0.04716  1.01644  -0.00078   -21.28156    4.92763
0.07000  3.39334   0.07032  1.02154  -0.00151   -14.33098    4.16305
0.08000  3.04287   0.10217  1.02672  -0.00273   -9.93608     3.59668
0.09000  2.71521   0.14649  1.03157  -0.00462   -7.02035     3.16549
0.10000  2.40341   0.20948  1.03553  -0.00744   -5.01990     2.83194
0.11000  2.10214   0.30205  1.03798  -0.01147   -3.61651     2.57226
0.12000  1.80669   0.44508  1.03823  -0.01701   -2.61923     2.37049
0.13000  1.51162   0.68402  1.03560  -0.02435   -1.90724     2.21523
0.14000  1.20829   1.14007  1.02957  -0.03371   -1.40028     2.09768
0.15000  0.87742   2.28196  1.01982  -0.04522   -1.04266     2.01064
0.16000  0.44359   9.35105  1.00630  -0.05890   -0.79447     1.94794
```

COMPUTED check (2026-10-04) of rows 0.01, 0.05, 0.10, 0.13, 0.14, 0.15, 0.16 from eqs. 8, 12, 13, 14 and the third of 6: a, e, x0, V and tau reproduce to the printed digits; C reproduces to 1e-5 except row 0.01, where C = -999.8827 computed against -999.87824 printed (the row is ill-conditioned: a is 1e-3 and C is the difference of large terms; the printed a has only 3 significant digits). Row 0.16 a = 9.35115 computed against 9.35105 printed (same conditioning). Treat rows 0.01 and 0.16 as controls at a looser tolerance.

### 2.2 Table 2, family A0, elliptic orbits (READ pp.394-395; 42 rows)

First row is the parabolic orbit (a = infinity, x1 = infinity). Sign columns for A0, READ from the images: (eps, eps', eps'') = (-, -, +) from the parabolic orbit to tau/pi = 0.4; at tau/pi = 0.5 (a = 1, e = 0, circular) only eps' = - is printed; (+, -, -) for tau/pi 0.6 to 0.9; at tau/pi = 1.0 (e = 1, rectilinear) eps = + and eps'' = - with eps' blank; (+, +, -) for 1.1 to 2.9; at 3.0 (e = 1) eps = +, eps' blank, eps'' = -; (+, -, -) for 3.1 to 4.0. A blank sign means the orbit is rectilinear or circular and that sign is undefined. eps eps'' = -1 throughout.

(The block of eight columns below is the table; the three sign columns are in the previous paragraph. Values come from the text layer and were verified by recomputing every row from eqs. 27 to 33 at a tolerance of 4e-3 on the closure residual, 2e-4 on x0 and 3e-4 on C and V; the one text-layer defect, a dropped minus sign on x0 at tau/pi = 2.3, was repaired from the image.)


#### Table 2, family A0 (elliptic), p.394-395 (42 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
0.16393 0.00000      inf 1.00000  -0.06485       inf  -0.72028 1.92880
0.17000 0.16734  6.92689 0.98922  -0.07467  13.77911  -0.62644 1.90432
0.18000 0.25967  2.97647 0.96897  -0.09237   5.86057  -0.51697 1.87536
0.20000 0.35543  1.67789 0.92089  -0.13274   3.22303  -0.41394 1.84768
0.30000 0.48269  1.03501 0.62227  -0.39096   1.67907  -0.62661 1.90437
0.40000 0.49875  1.00122 0.31255  -0.68829   1.31416  -0.90219 1.97540
0.50000 0.50000  1.00000 0.00000  -1.00000   1.00000  -1.00000 2.00000
0.60000 0.50059  1.00057 0.31069  -1.31144   0.68971  -0.90214 1.97538
0.70000 0.50358  1.00673 0.59509  -1.60583   0.40763  -0.61940 1.90247
0.80000 0.50967  1.02552 0.81926  -1.86568   0.18535  -0.18627 1.78501
0.90000 0.51924  1.06131 0.95651  -2.07647   0.04616   0.34122 1.63058
1.00000 0.53286  1.11489 1.00000  -2.22978   0.00000   0.89695 1.45019
1.10000 0.55140  1.18349 0.96437  -2.32481   0.04216   1.42056 1.25676
1.20000 0.57584  1.26118 0.87748  -2.36784   0.15452   1.87013 1.06295
1.30000 0.60708  1.34003 0.76872  -2.37014   0.30992   2.22700 0.87920
1.40000 0.64565  1.41212 0.66061  -2.34498   0.47926   2.49238 0.71248
1.50000 0.69157  1.47175 0.56616  -2.30499   0.63851   2.67946 0.56616
1.60000 0.74430  1.51659 0.49058  -2.26059   0.77258   2.80562 0.44088
1.70000 0.80288  1.54743 0.43445  -2.21971   0.87515   2.88709 0.33602
1.80000 0.86603  1.56705 0.39646  -2.18832   0.94578   2.93661 0.25177
1.90000 0.93229  1.57910 0.37518  -2.17155   0.98665   2.96293 0.19255
2.00000 1.00000  1.58740 0.37004  -2.17480   1.00000   2.97093 0.17049
2.10000 1.06719  1.59580 0.38183  -2.20512   0.98648   2.96172 0.19566
2.20000 1.13144  1.60823 0.41290  -2.27228   0.94419   2.93182 0.26112
2.30000 1.18980  1.62884 0.46658  -2.38882   0.86886   2.87159 0.35835
2.40000 1.23917  1.66133 0.54474  -2.56632   0.75633   2.76372 0.48609
2.50000 1.27742  1.70735 0.64366  -2.80630   0.60840   2.58570 0.64366
2.60000 1.30439  1.76490 0.75167  -3.09153   0.43828   2.31900 0.82523
2.70000 1.32166  1.82883 0.85282  -3.38848   0.26917   1.95921 1.02019
2.80000 1.33126  1.89302 0.93299  -3.65920   0.12685   1.51860 1.21713
2.90000 1.33475  1.95242 0.98324  -3.87212   0.03271   1.02161 1.40655
3.00000 1.33297  2.00400 1.00000  -4.00800   0.00000   0.49900 1.58145
3.10000 1.32605  2.04714 0.98427  -4.06208   0.03221  -0.01711 1.73698
3.20000 1.31360  2.08357 0.94097  -4.04415   0.12300  -0.49728 1.87010
3.30000 1.29489  2.11699 0.87835  -3.97645   0.25752  -0.91862 1.97955
3.40000 1.26912  2.15212 0.80698  -3.88884   0.41541  -1.26812 2.06594
3.50000 1.23595  2.19343 0.73763  -3.81137   0.57550  -1.54409 2.13169
3.60000 1.19588  2.24373 0.67884  -3.76685   0.72060  -1.75411 2.18039
3.70000 1.15030  2.30353 0.63541  -3.76721   0.83984  -1.90980 2.21581
3.80000 1.10114  2.37137 0.60878  -3.81501   0.92773  -2.02167 2.24091
3.90000 1.05043  2.44457 0.59843  -3.90746   0.98168  -2.09623 2.25748
4.00000 1.00000  2.51984 0.60315  -4.03968   1.00000  -2.13546 2.26616
```

#### Table 3, family A1, p.395-396 (58 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
4.00000 1.00000  2.51984 0.60315  -4.03968   1.00000   2.92916 0.26616
3.90000 0.95040  2.51043 0.60904  -4.03939   0.98147   2.91169 0.29717
3.80000 0.90209  2.49759 0.62914  -4.06894   0.92625   2.85720 0.37788
3.70000 0.85641  2.47825 0.66279  -4.12080   0.83569   2.76112 0.48876
3.60000 0.81451  2.45014 0.70884  -4.18690   0.71337   2.61635 0.61940
3.50000 0.77734  2.41232 0.76515  -4.25811   0.56652   2.41454 0.76515
3.40000 0.74560  2.36555 0.82790  -4.32399   0.40712   2.14803 0.92302
3.30000 0.71979  2.31241 0.89110  -4.37300   0.25182   1.81263 1.08966
3.20000 0.70018  2.25679 0.94672  -4.39332   0.12025   1.41076 1.26065
3.10000 0.68689  2.20294 0.98570  -4.37439   0.03150   0.95410 1.43035
3.00000 0.68000  2.15432 1.00000  -4.30865   0.00000   0.46418 1.59242
2.90000 0.67958  2.11250 0.98490  -4.19312   0.03189  -0.02982 1.74064
2.80000 0.68580  2.07664 0.94071  -4.03015   0.12313  -0.49610 1.86979
2.70000 0.69895  2.04356 0.87274  -3.82707   0.26005  -0.90638 1.97646
2.60000 0.71935  2.00857 0.78970  -3.59474   0.42240  -1.24106 2.05938
2.50000 0.74730  1.96669 0.70109  -3.34552   0.58786  -1.49153 2.11932
2.40000 0.78299  1.91400 0.61501  -3.09114   0.73687  -1.65932 2.15855
2.30000 0.82639  1.84864 0.53696  -2.84129   0.85599  -1.75307 2.18015
2.20000 0.87730  1.77091 0.46979  -2.60286   0.93896  -1.78484 2.18743
2.10000 0.93534  1.68286 0.41429  -2.38005   0.98566  -1.76714 2.18338
2.00000 1.00000  1.58740 0.37004  -2.17480   1.00000  -1.71101 2.17049
1.90000 1.07079  1.48759 0.33605  -1.98750   0.98769  -1.62525 2.15064
1.80000 1.24723  1.38609 0.31125  -1.81752   0.95467  -1.51623 2.12514
1.70000 1.22898  1.28494 0.29479  -1.66372   0.90615  -1.38811 2.09478
1.60000 1.31584  1.18548 0.28614  -1.52470   0.84627  -1.24301 2.05986
1.50000 1.40803  1.08836 0.28493  -1.39846   0.77825  -1.08118 2.02019
1.45000 1.45675  1.04028 0.28585  -1.33765   0.74292  -0.99350 1.99837
1.43000 1.47766  1.02031 0.28392  -1.31000   0.73063  -0.95698 1.98922
1.42000 1.49124  1.00761 0.27434  -1.28404   0.73118  -0.93812 1.98447
1.42000 1.49527  1.00392 0.26259  -1.26754   0.74030  -0.93750 1.98431
1.43000 1.49877  1.00086 0.22183  -1.22288   0.77884  -0.95187 1.98793
1.45000 1.49977  1.00011 0.15713  -1.15726   0.84296  -0.97538 1.99384
1.50000 1.50000  1.00000 0.00000  -1.00000   1.00000  -1.00000 2.00000
1.60000 1.50039  1.00038 0.31011  -0.69015   1.31061  -0.90213 1.97538
1.70000 1.50202  1.00377 0.59192  -0.40962   1.59792  -0.61878 1.90231
1.80000 1.50517  1.01340 0.81455  -0.18793   1.83887  -0.18115 1.78358
1.90000 1.51018  1.03147 0.95402  -0.04743   2.01552   0.36064 1.62461
2.00000 1.51771  1.05888 1.00000   0.00000   2.11776   0.94440 1.43374
2.10000 1.52882  1.09496 0.95901  -0.04488   2.14503   1.50633 1.22216
2.20000 1.54520  1.13727 0.85289  -0.16731   2.10724   1.99285 1.00357
2.30000 1.56915  1.18157 0.71299  -0.33912   2.02402   2.37069 0.79329
2.40000 1.60320  1.22249 0.57133  -0.52405   1.92094   2.63288 0.60590
2.50000 1.64897  1.25549 0.45111  -0.68913   1.82186   2.79650 0.45111
2.60000 1.70608  1.27882 0.36150  -0.81653   1.74111   2.89072 0.33058
2.70000 1.77236  1.29360 0.30062  -0.90472   1.68249   2.94255 0.23969
2.80000 1.84505  1.30225 0.26260  -0.96028   1.64422   2.97013 0.17284
2.90000 1.92165  1.30717 0.24229  -0.99045   1.62388   2.98351 0.12842
3.00000 2.00000  1.31037 0.23686  -1.00000   1.62074   2.98742 0.11214
3.10000 2.07804  1.31361 0.24610  -0.99033   1.63689   2.98302 0.13031
3.20000 2.15333  1.31872 0.27272  -0.95908   1.67836   2.96796 0.17900
3.30000 2.22258  1.32795 0.32269  -0.89943   1.75646   2.93448 0.25596
3.40000 2.28150  1.34403 0.40387  -0.80122   1.88684   2.86517 0.36720
3.50000 2.32624  1.36903 0.51919  -0.65825   2.07981   2.73045 0.51919
3.60000 2.35617  1.40225 0.65699  -0.48098   2.32352   2.49861 0.70809
3.70000 2.37414  1.44005 0.79335  -0.29759   2.58251   2.15549 0.91897
3.80000 2.38364  1.47807 0.90481  -0.14069   2.81545   1.71191 1.13494
3.90000 2.38708  1.51292 0.97599  -0.03633   2.98952   1.19680 1.34283
4.00000 2.38558  1.54267 1.00000   0.00000   3.08534   0.64823 1.53355
```

#### Table 4, family A2, p.397 (39 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
4.00000 1.62647  1.63120 1.00000   0.00000   3.26240   0.61304 1.54498
3.90000 1.62629  1.60756 0.97804  -0.03530   3.17983   0.09356 1.70483
3.80000 1.63271  1.58805 0.91440  -0.13593   3.04016  -0.39056 1.84135
3.70000 1.64644  1.57057 0.81824  -0.28547   2.85568  -0.80419 1.95043
3.60000 1.66853  1.55197 0.70418  -0.45910   2.64484  -1.12472 2.03094
3.50000 1.70014  1.52882 0.58813  -0.62967   2.42797  -1.34590 2.08468
3.40000 1.74197  1.49858 0.48284  -0.77500   2.22216  -1.47672 2.11583
3.30000 1.79395  1.46042 0.39522  -0.88323   2.03760  -1.53545 2.12966
3.20000 1.85517  1.41514 0.32658  -0.95298   1.87730  -1.54210 2.13122
3.10000 1.92431  1.36450 0.27486  -0.98945   1.73955  -1.51338 2.12447
3.00000 2.00000  1.31037 0.23686  -1.00000   1.62074  -1.46114 2.11214
2.90000 2.08107  1.25436 0.20954  -0.99152   1.51719  -1.39302 2.09595
2.80000 2.16665  1.19764 0.19055  -0.96943   1.42585  -1.31366 2.07694
2.70000 2.25612  1.14103 0.17825  -0.93764   1.34442  -1.22576 2.05567
2.60000 2.34912  1.08501 0.17165  -0.89877   1.27125  -1.13070 2.03241
2.50000 2.44556  1.02983 0.17018  -0.85457   1.20508  -1.02896 2.00723
2.47000 2.47530  1.01338 0.17038  -0.84073   1.18604  -0.99711 1.99928
2.46000 2.48544  1.00784 0.17009  -0.83642   1.17926  -0.98635 1.99658
2.45500 2.49070  1.00498 0.16943  -0.83470   1.17525  -0.98093 1.99523
2.45016 2.49700  1.00156 0.16512  -0.83619   1.16693  -0.97564 1.99390
2.45500 2.49955  1.00020 0.14228  -0.85789   1.14251  -0.98005 1.99501
2.46000 2.49980  1.00008 0.12595  -0.87412   1.12604  -0.98423 1.99605
2.47000 2.49995  1.00001 0.09426  -0.90575   1.09428  -0.99112 1.99778
2.50000 2.50000  1.00000 0.00000  -1.00000   1.00000  -1.00000 2.00000
2.60000 2.50029  1.00028 0.30983  -1.31020   0.69036  -0.90213 1.97538
2.70000 2.50141  1.00262 0.59067  -1.59484   0.41040  -0.61855 1.90225
2.80000 2.50353  1.00909 0.81281  -1.82930   0.18889  -0.17934 1.78307
2.90000 2.50693  1.02118 0.95309  -1.99446   0.04790   0.36752 1.62249
3.00000 2.51213  1.03961 1.00000  -2.07921   0.00000   0.96190 1.42762
3.10000 2.52007  1.06414 0.95673  -2.08223   0.04604   1.54004 1.20829
3.20000 2.53237  1.09341 0.84143  -2.01344   0.17338   2.04465 0.97742
3.30000 2.55177  1.12467 0.68456  -1.89457   0.35477   2.43528 0.75148
3.40000 2.58210  1.15389 0.52287  -1.75723   0.55055   2.69794 0.54960
3.50000 2.62680  1.17710 0.38789  -1.63369   0.72052   2.84954 0.38789
3.60000 2.68608  1.19270 0.29277  -1.54188   0.84352   2.92695 0.27029
3.70000 2.75668  1.20192 0.23275  -1.48166   0.92217   2.96443 0.18861
3.80000 2.83451  1.20697 0.19758  -1.44545   0.96849   2.98245 0.13249
3.90000 2.91639  1.20970 0.17950  -1.42684   0.99255   2.99065 0.09671
4.00000 3.00000  1.21141 0.17452  -1.42283   1.00000   2.99299 0.08375
```

#### Table 5, family B1, p.398-399 (65 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
4.00000 1.29827  2.45202 1.00000   4.90405   0.00000   0.40783 1.61002
3.90000 1.29837  2.40588 0.98722   4.78101  -0.03075   0.91001 1.44568
3.80000 1.29288  2.35360 0.94947   4.58826  -0.11894   1.38793 1.26967
3.70000 1.28122  2.29791 0.89025   4.34363  -0.25220   1.31608 1.08808
3.60000 1.26254  2.24299 0.81649   4.07438  -0.41161   2.17521 0.90818
3.50000 1.23595  2.19343 0.73763   3.81137  -0.57550   2.45591 0.73763
3.40000 1.20106  2.15287 0.66352   3.58135  -0.72439   2.65998 0.58311
3.30000 1.15836  2.12285 0.60190   3.40059  -0.84511   2.79811 0.44932
3.20000 1.10927  2.10256 0.55688   3.27343  -0.93169   2.88437 0.34005
3.10000 1.05578  2.08945 0.52951   3.19584  -0.98306   2.93102 0.26264
3.00000 1.00000  2.08008 0.51925   3.16017  -1.00000   2.94591 0.23258
2.90000 0.94389  2.07075 0.52522   3.15836  -0.98315   2.93201 0.26075
2.80000 0.88918  2.05776 0.54685   3.18304  -0.93248   2.88797 0.33471
2.70000 0.83736  2.03769 0.58382   3.22733  -0.84805   2.80865 0.43744
2.60000 0.78971  2.00781 0.63569   3.28417  -0.73146   2.68569 0.56063
2.50000 0.74730  1.96669 0.70109   3.34552  -0.58786   2.50847 0.70109
2.40000 0.71092  1.91478 0.77658   3.40175  -0.42780   2.26584 0.85683
2.30000 0.68107  1.85478 0.85556   3.44166  -0.26791   1.94927 1.02505
2.20000 0.65796  1.79127 0.92777   3.45316  -0.12938   1.55710 1.20121
2.10000 0.64160  1.72955 0.98022   3.42489  -0.03422   1.09880 1.37884
2.00000 0.63192  1.67416 1.00000   3.34832   0.00000   0.59732 1.55006
1.90000 0.62892  1.62743 0.97843   3.21975  -0.03511   0.08736 1.70665
1.80000 0.63280  1.58868 0.91445   3.04145  -0.13591  -0.39073 1.84139
1.70000 0.64404  1.55409 0.81545   2.82136  -0.28681  -0.79972 1.94929
1.60000 0.66335  1.51744 0.69458   2.57142  -0.46346  -1.11341 2.02815
1.50000 0.69157  1.47175 0.56616   2.30499  -0.63851  -1.32054 2.07859
1.40000 0.72957  1.41125 0.44133   2.03407  -0.78842  -1.42343 2.10319
1.30000 0.77813  1.33299 0.32580   1.76728  -0.89870  -1.43292 2.10545
1.20000 0.83809  1.23711 0.21944   1.50858  -0.96563  -1.36195 2.08853
1.10000 0.91086  1.12564 0.11614   1.25636  -0.99491  -1.21917 2.05406
1.00000 1.00000  1.00000 0.00000   1.00000  -1.00000  -1.00000 2.00000
0.91375 1.10000  0.87591 0.14896   0.74543  -1.00639  -0.70925 1.92594
0.86101 1.20000  0.77250 0.36401   0.49130  -1.05370  -0.34275 1.82832
0.87251 1.30000  0.70088 0.72607   0.19199  -1.20977   0.27545 1.65062
1.00000 1.36836  0.71333 1.00000   0.00000  -1.42666   1.40188 1.26417
1.10000 1.38693  0.76131 0.90143   0.07504  -1.44758   2.06901 0.96488
1.20000 1.39725  0.82653 0.66157   0.27972  -1.37334   2.57337 0.65317
1.30000 1.40325  0.90516 0.35008   0.58828  -1.22204   2.88717 0.33591
1.40000 1.40658  0.99377 0.02166   0.97225  -1.01530   2.99956 0.02094
1.50000 1.40803  1.08836 0.28493   1.39846  -0.77825   2.91882 0.28493
1.60000 1.40799  1.18432 0.54598   1.83094  -0.53771   2.66786 0.57632
1.70000 1.40658  1.27691 0.74961   2.23410  -0.31972   2.27899 0.84912
1.80000 1.40373  1.36177 0.89194   2.57639  -0.14715   1.78962 1.10017
1.90000 1.39915  1.43564 0.97401   2.83397  -0.03732   1.23938 1.32688
2.00000 1.39229  1.49691 1.00000   2.99383   0.00000   0.66804 1.52708
2.10000 1.38222  1.54610 0.97675   3.05625  -0.03595   0.11367 1.69892
2.20000 1.36756  1.58608 0.91425   3.03615  -0.13601  -0.39002 1.84120
2.30000 1.34640  1.62202 0.82641   2.96248  -0.28156  -0.81770 1.95389
2.40000 1.31674  1.66046 0.73056   2.87352  -0.44740  -1.15758 2.03901
2.50000 1.27742  1.70735 0.64366   2.80630  -0.60840  -1.41430 2.10102
2.60000 1.22917  1.76574 0.57682   2.78425  -0.74723  -1.60461 2.14583
2.70000 1.17448  1.83519 0.53321   2.81374  -0.85664  -1.74718 2.17880
2.80000 1.11635  1.91304 0.51104   2.89067  -0.93541  -1.85503 2.20341
2.90000 1.05750  1.99586 0.50722   3.00820  -0.98353  -1.93403 2.22127
3.00000 1.00000  2.08008 0.51925   3.16017  -1.00000  -1.98441 2.23258
3.10000 0.94545  2.16215 0.54549   3.34158  -0.98273  -2.00228 2.23658
3.20000 0.89508  2.23871 0.58479   3.54788  -0.92953  -1.98075 2.23176
3.30000 0.84984  2.30689 0.63598   3.77402  -0.83975  -1.91072 2.21601
3.40000 0.81043  2.36485 0.69717   4.01355  -0.71616  -1.78208 2.18680
3.50000 0.77734  2.41232 0.76515   4.25811  -0.56652  -1.58546 2.14137
3.60000 0.75082  2.45083 0.83502   4.49732  -0.40434  -1.31472 2.07719
3.70000 0.73096  2.48358 0.90023   4.71937  -0.24779  -0.96974 1.99242
3.80000 0.71766  2.51469 0.95346   4.91233  -0.11704  -0.55866 1.88644
3.90000 0.71079  2.54809 0.98811   5.06587  -0.03031  -0.09849 1.76025
4.00000 0.71019  2.58651 1.00000   5.17302   0.00000   0.38662 1.61659
```

#### Table 6, family B2, p.399-400 (47 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
4.00000 2.00000  1.58740 0.37004   1.00000  -2.17480   2.97093 0.17049
3.90000 1.93242  1.58324 0.37684   0.98660  -2.17987   2.96263 0.19333
3.80000 1.86666  1.57716 0.40058   0.94538  -2.20895   2.93542 0.25412
3.70000 1.80474  1.56722 0.44263   0.87352  -2.26092   2.88322 0.34173
3.60000 1.74863  1.55154 0.50490   0.76817  -2.33491   2.79489 0.45289
3.50000 1.70014  1.52882 0.58813   0.62967  -2.42797   2.65410 0.58813
3.40000 1.66046  1.49902 0.68918   0.46593  -2.53211   2.44140 0.74740
3.30000 1.62991  1.46388 0.79842   0.29509  -2.63267   2.14010 0.92731
3.20000 1.60788  1.42664 0.89953   0.14333  -2.70995   1.74452 1.12048
3.10000 1.59331  1.39094 0.97270   0.03797  -2.74392   1.26632 1.31669
3.00000 1.58522  1.35968 1.00000   0.00000  -2.71937   0.73547 1.50484
2.90000 1.58308  1.33426 0.97084   0.03890  -2.62962   0.19568 1.67461
2.80000 1.58702  1.31423 0.88558   0.15038  -2.47809  -0.30410 1.81772
2.70000 1.59800  1.29725 0.75614   0.31635  -2.27816  -0.71985 1.92869
2.60000 1.61784  1.27929 0.60338   0.50740  -2.05117  -1.02225 2.00555
2.50000 1.64897  1.25549 0.45111   0.68913  -1.82186  -1.20350 2.05024
2.40000 1.69346  1.22203 0.31817   0.83322  -1.61084  -1.27770 2.06826
2.30000 1.75189  1.17780 0.21223   0.92783  -1.42777  -1.27204 2.06689
2.20000 1.82323  1.12444 0.13024   0.97799  -1.27089  -1.21340 2.05266
2.10000 1.90606  1.06455 0.06337   0.99708  -1.13201  -1.12003 2.02978
2.00000 2.00000  1.00000 0.00000   1.00000  -1.00000  -1.00000 2.00000
1.90000 2.10810  0.93059 0.07910   1.00420   0.85698  -0.84871 1.96181
1.80000 2.25648  0.84501 0.26484   1.06880  -0.62122  -0.58942 1.89458
1.78733 2.30000  0.82291 0.36612   1.12419  -0.52163  -0.47312 1.86363
1.80000 2.35584  0.79898 0.57500   1.25838  -0.33957  -0.21103 1.79193
1.90000 2.40819  0.79369 0.91382   1.51898  -0.06840   0.53632 1.56961
2.00000 2.42546  0.81166 1.00000   1.62332   0.00000   1.23205 1.32964
2.10000 2.43456  0.84089 0.92687   1.62029  -0.06150   1.87768 1.05940
2.20000 2.43991  0.87913 0.73257   1.52316  -0.23511   2.41394 0.76554
2.30000 2.44310  0.92472 0.45782   1.34807  -0.50137   2.79127 0.45688
2.40000 2.44486  0.97572 0.14435   1.11657  -0.83488   2.97976 0.14225
2.50000 2.44556  1.02983 0.17018   0.85457  -1.20508   2.97104 0.17018
2.60000 2.44537  1.08440 0.45575   0.59019  -1.57862   2.77599 0.47330
2.70000 2.44430  1.13678 0.69115   0.35109  -1.92247   2.42078 0.76106
2.80000 2.44226  1.18452 0.86341   0.16179  -2.20724   1.94238 1.02841
2.90000 2.43897  1.22574 0.96646   0.04111  -2.41038   1.38448 1.27103
3.00000 2.43395  1.25947 1.00000   0.00000  -2.51893   0.79399 1.48527
3.10000 2.42632  1.28584 0.96904   0.03981  -2.53186   0.21774 1.66801
3.20000 2.41457  1.30638 0.88445   0.15096  -2.46180  -0.30123 1.81693
3.30000 2.39621  1.32414 0.76422   0.31221  -2.33608  -0.72912 1.93109
3.40000 2.36777  1.34354 0.63358   0.49230  -2.19478  -1.04926 2.01228
3.50000 2.32624  1.36903 0.51919   0.65825  -2.07981  -1.26955 2.06629
3.60000 2.27175  1.40272 0.43687   0.78991  -2.01553  -1.41783 2.10186
3.70000 2.20793  1.44366 0.38698   0.88499  -2.00234  -1.52314 2.12677
3.80000 2.13918  1.48954 0.36278   0.94916  -2.02992  -1.60329 2.14553
3.90000 2.06904  1.53810 0.35824   0.98709  -2.08911  -1.66563 2.16001
4.00000 2.00000  1.58740 0.37004   1.00000  -2.17480  -1.71101 2.17049
```

#### Table 7, family C12, p.400-401 (26 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
1.00000 2.00000  0.62996 0.58740  -1.00000   0.25992   2.87208 0.35766
1.04586 2.10000  0.61575 0.65614  -1.01977   0.21173   2.80836 0.43777
1.05021 2.20000  0.58186 0.88828  -1.09871   0.06501   2.41935 0.76201
1.00000 2.20824  0.55756 1.00000  -1.11512   0.00000   1.79353 1.09839
0.98624 2.20000  0.55498 0.99118  -1.10506   0.00489   1.60443 1.18134
0.96270 2.10000  0.58086 0.75873  -1.02157   0.14014   0.72867 1.50709
1.00000 2.00000  0.62996 0.58740  -1.00000   0.25992   0.30272 1.64234
1.10000 1.84377  0.72566 0.42865  -1.03672   0.41461  -0.16122 1.77798
1.20000 1.71909  0.81492 0.35753  -1.10628   0.52357  -0.45901 1.85984
1.30000 1.60714  0.90454 0.31954  -1.19358   0.61550  -0.69688 1.92273
1.35000 1.55384  0.95051 0.30930  -1.24450   0.65652  -0.80221 1.94993
1.37000 1.53226  0.96974 0.30835  -1.26877   0.67072  -0.84234 1.96019
1.38000 1.52046  0.98040 0.31125  -1.28554   0.67525  -0.86194 1.96518
1.38406 1.51400  0.98624 0.31723  -1.29911   0.67338  -0.86966 1.96714
1.38000 1.50713  0.99225 0.34859  -1.33814   0.64636  -0.85946 1.96455
1.37000 1.50582  0.99308 0.38165  -1.37208   0.61407  -0.83523 1.95837
1.35000 1.50548  0.99248 0.44021  -1.42938   0.55557  -0.78144 1.94459
1.30000 1.50696  0.98762 0.57329  -1.55381   0.42143  -0.61600 1.90158
1.20000 1.51348  0.96748 0.79387  -1.73553   0.19943  -0.16262 1.77838
1.10000 1.52538  0.93014 0.94283  -1.80710   0.05318   0.43224 1.60242
1.00000 1.54792  0.86958 1.00000  -1.73915   0.00000   1.14998 1.36015
0.90716 1.60000  0.77836 0.92148  -1.49560   0.06112   1.97012 1.01483
0.87673 1.70000  0.69631 0.74201  -1.21298   0.17964   2.55495 0.66712
0.90297 1.80000  0.66074 0.63468  -1.08009   0.24138   2.76978 0.47981
0.94855 1.90000  0.64236 0.58540  -1.01840   0.26632   2.85633 0.37903
1.00000 2.00000  0.62996 0.58740  -1.00000   0.25992   2.87208 0.35766
```

#### Table 8, family C23, p.401-402 (31 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
2.00000 3.00000  0.76314 0.31037  -0.52629   1.00000   2.97125 0.16956
2.07094 3.10000  0.75884 0.33416  -0.50527   1.01241   2.95988 0.20030
2.13186 3.20000  0.75083 0.41019  -0.44285   1.05882   2.91236 0.29604
2.14993 3.30000  0.72761 0.63691  -0.26418   1.19103   2.68958 0.55715
2.10000 3.33569  0.70147 0.86223  -0.09664   1.30630   2.27403 0.85204
2.00000 3.33900  0.67362 1.00000   0.00000   1.34725   1.48451 1.23105
1.91025 3.30000  0.66527 0.85602  -0.09579   1.23475   0.65988 1.52975
1.89079 3.20000  0.68925 0.55727  -0.30515   1.07336   0.07214 1.71110
1.93428 3.10000  0.72407 0.40069  -0.43394   1.01420  -0.17817 1.78274
2.00000 3.00000  0.76314 0.31037  -0.52629   1.00000  -0.35051 1.83044
2.07736 2.90000  0.80523 0.25434  -0.60043   1.01002  -0.49378 1.86917
2.16211 2.80000  0.84996 0.21820  -0.66450   1.03542  -0.62291 1.90340
2.25222 2.70000  0.89728 0.19477  -0.72251   1.07204  -0.74373 1.93487
2.34660 2.60000  0.94722 0.18030  -0.77644   1.11801  -0.85889 1.96441
2.42455 2.52000  0.98916 0.17461  -0.81644   1.16187  -0.94761 1.98686
2.43408 2.51000  0.99452 0.17533  -0.82000   1.16890  -0.95811 1.98950
2.43883 2.50300  0.99829 0.18190  -0.81670   1.17988  -0.96324 1.99079
2.43500 2.50142  0.99912 0.19852  -0.80078   1.19746  -0.95845 1.98958
2.43000 2.50111  0.99925 0.21482  -0.78459   1.21391  -0.95183 1.98792
2.42000 2.50097  0.99925 0.24583  -0.75361   1.24489  -0.93715 1.98423
2.40000 2.50105  0.99899 0.30602  -0.69327   1.30471  -0.90207 1.97537
2.30000 2.50277  0.99496 0.58206  -0.41583   1.57408  -0.61712 1.90187
2.20000 2.50608  0.98491 0.80232  -0.19470   1.77513  -0.16944 1.78029
2.10000 2.51167  0.96645 0.94743  -0.05080   1.88210   0.40562 1.61071
2.00000 2.52128  0.93739 1.00000   0.00000   1.87479   1.06679 1.39040
1.90000 2.53992  0.89503 0.93750  -0.05594   1.73413   1.77572 1.10647
1.80288 2.60000  0.82737 0.67519  -0.26874   1.38600   2.55057 0.67039
1.80791 2.70000  0.78843 0.45653  -0.42849   1.14837   2.84835 0.38942
1.86054 2.80000  0.77397 0.36099  -0.49458   1.05336   2.93291 0.25902
1.92786 2.90000  0.76725 0.31898  -0.52251   1.01198   2.96370 0.19052
2.00000 3.00000  0.76314 0.31037  -0.52629   1.00000   2.97125 0.16956
```

#### Table 9, family C24, p.402 (16 rows)

```
tau/pi  eta/pi    a        e        x0        x1        C         V
2.00000 4.00000  0.62996 0.58740   1.00000  -0.25992   2.87208 0.35766
2.04748 4.10000  0.62312 0.63595   1.01940  -0.22685   2.82320 0.42047
2.07037 4.20000  0.60978 0.79101   1.09212  -0.12744   2.59544 0.63605
2.00000 4.24864  0.58475 1.00000   1.16950   0.00000   1.71013 1.13572
1.94797 4.20000  0.58400 0.88047   1.09820  -0.06980   0.98771 1.41855
1.95681 4.10000  0.60400 0.68913   1.02037  -0.18779   0.52900 1.57194
2.00000 4.00000  0.62996 0.58740   1.00000  -0.25992   0.30272 1.64234
2.05453 3.90000  0.65833 0.54570   1.01758  -0.29908   0.15916 1.68548
2.10870 3.80000  0.68800 0.56053   1.07365  -0.30235   0.07968 1.70889
2.14144 3.70000  0.71701 0.67147   1.19846  -0.23556   0.13971 1.69124
2.10000 3.63950  0.72741 0.88312   1.36979  -0.08502   0.57448 1.55741
2.00000 3.63397  0.70994 1.00000   1.41989   0.00000   1.40856 1.26152
1.90569 3.70000  0.66894 0.84199   1.23218  -0.10570   2.37740 0.78905
1.90964 3.80000  0.64651 0.67584   1.08345  -0.20958   2.73203 0.51766
1.94976 3.90000  0.63635 0.60086   1.01871  -0.25399   2.84677 0.39145
2.00000 4.00000  0.62996 0.58740   1.00000  -0.25992   2.87208 0.35766
```


Sign columns for Tables 3 to 9: the values of eps, eps', eps'' are printed on every row of the table images (journal pp.395-402) and were read to confirm the family rule eps eps'' = -1 (A), +1 (B), (-1)^(i+j) (C_ij), and that eps' is blank on the rectilinear (e = 1) rows and the circular (e = 0) rows. They are NOT transcribed row by row here; the signs enter the numbers only through the closure test, which every row passed (section 2.3) with the printed sign pattern or its equivalent. Anyone building a sign-dependent test must read the image.

### 2.3 Printed values that contradict the paper's own equations (COMPUTED, 2026-10-04)

All rows of Tables 2 to 9 (324 rows) were recomputed from eqs. 27 to 33 (closure residual below 4e-3, x0 within 2e-4, C and V within 3e-4). Everything agrees except these, which are printing or text-layer defects:

| Table, row | Printed | Computed | Verdict |
|---|---|---|---|
| Table 3 (A1), tau/pi = 1.80000 | eta/pi = 1.24723 | 1.14723 (root of eq. 30 is 1.147232) | one wrong digit in the printed eta/pi; a, e, x0, x1, C, V on that row are consistent with 1.14723 and the neighbouring rows (1.07079 at 1.9, 1.22898 at 1.7 as printed) |
| Table 5 (B1), tau/pi = 3.70000 | C = 1.31608 (V = 1.08808) | C = 1.81607 (from a = 2.29791, e = 0.89025, eps' = +1) | the printed V = 1.08808 matches C = 1.81607, so the C digit "3" is a misprint of "8"; the sequence 0.91001, 1.38793, 1.81607, 2.17521 is smooth |
| Table 6 (B2), tau/pi = 1.90000, eta/pi = 2.10810 | x1 = 0.85698 | -0.85698 by eq. 32 | probably a dropped minus sign (not confirmed against a second copy) |

The text layer of the PDF is not reliable for the tables: it drops minus signs (Table 2 tau/pi = 2.3, x0; Table 8 tau/pi = 2.07736, x0) and breaks some numbers across spaces ("2. 4591" for C = 2.94591 at B1 tau/pi 3.0, "2.03201" for 2.93201 at tau/pi 2.9, "2. 36*85" for a at B1 tau/pi 3.4 in the second half, where x0 = 4.01355 = a (1 + e) gives a = 2.36485 by eq. 32 (COMPUTED; that row's a was not read on the image)). The values in the blocks above include these repairs. Anyone building a test from this paper must read the printed table image.

### 2.4 Independent reproductions (COMPUTED, 2026-10-04) suitable as test controls

Each of these is a printed value reproduced by me from the printed equation, so it is a sourced test whose expected side is the paper's number, not our code's:

- Parabolic orbit, eq. 20: tau/pi = 0.16393 (computed 0.163926), x0 = -0.06485 (-0.064851), C = -0.72028 (-0.720283), V = 1.92880 (1.928804), sigma = 3.7973 (3.79736).
- Table 2 / Bruno Table IV row tau/pi = 0.17: root of eq. 30 with eps eps'' = -1 gives eta/pi = 0.167344 (printed 0.16734); at tau/pi = 0.3 gives 0.482687 (printed 0.48269).
- B1 at tau/pi = 1 (Table 5 e = 1 row): root of eq. 30 with eps eps'' = +1 gives eta/pi = 1.368358 (printed 1.36836), a = 0.71333, e = 1 (this is the "interior orbit" Bruno uses at the Moon).
- Eq. 41 tangent ellipses: for tau/pi = 2, eta/pi = 1 (the A0 row at tau/pi = 2.0, eta/pi = 1.0): a = (2/1)^(2/3) = 1.58740, e = |1 - (1/2)^(2/3)| = 0.37004, which equals the printed a = 1.58740, e = 0.37004 (Table 2, row 2.0). For tau/pi = 4, eta/pi = 1: a = 4^(2/3) = 2.51984, e = |1 - 4^(-2/3)| = 0.60315, printed 2.51984 and 0.60315 (Table 2 last row and Table 3 first row). For tau/pi = 2, eta/pi = 3 (C23 first row, i < j): a = (2/3)^(2/3) = 0.76314, e = 0.31037, printed 0.76314, 0.31037. For tau/pi = 1, eta/pi = 2 (C12 first row): a = 0.62996, e = 0.58740, printed.
- Self-consistency of C and V: V = sqrt(3 - C) holds on every row to 3e-4.

Not reproduced independently: the rows' eta/pi digits near turning points (the root of (30) is badly conditioned where the curve is horizontal in the (tau, eta) plane).

## 3. Where this paper sits relative to the others (INFERRED from the digests)

- The timing equation (30) is exactly the equation "F0 = 0" of Hitzl & Henon 1977 (their eq. 11, with sigma = eps eps''), and the collision-geometry equations (27), (29) are their (5), (6), (8): the equation forms are identical term by term (INFERRED from reading both). Hitzl & Henon supply what this paper lacks: the Jacobi constant on the arc in closed form (their eq. 21, C = 2 sigma sin(tau)/rho + sin^2(eta)/rho^2), its extrema (the critical orbits), and the analytic conditions.
- The 1997 book recasts the same arcs as the S-arc families with integer pair (alpha, beta) at fixed Jacobi constant (its eqs. 4.27 and 4.28 are equivalent to relation (30), per the book itself, p.47); this paper's family names A, B, C survive in the book's Table 4.4 (A0, A1, A2, B1, B2, ...).
- Bruno 1981 and Bruno 1972 use this paper's notation without change (same eps, eps', eps''), reuse Figure 4 as his Figure 2, take his Table III from the e = 1 rows of Tables 5, 6, 7, 8, 9 and his Table IV from Table 2 (Henon sent him an unpublished finer table).
- Egorov 1958 had pieces of A0, A1, A2, B1 and B2 (Bruno's Figure 3), which Egorov used to argue that orbits around the Moon with near-Earth passages are impossible; Bruno says three solutions of this paper not found by Egorov make them possible.

## 4. Reconciliation with the project code and ledger

Searched `src/cyclerfinder` for second-species machinery on 2026-10-04: there is none. `grep -i "second.species\|consecutive.collision"` finds only comments in `search/earth_moon_resonant_families.py` and `search/earth_moon_class1_resonant_connections.py` (which say that Casoliva's elliptical/second-species differential-correction method "never" was built here) and a search-query string in `search/literature_check.py`. No code evaluates eq. 30, eq. 29 or the C of eq. 33.

Agreements:
- The `#899` enumerator specification (data/OUTSTANDING.md, from the 1997 book, chapter 4) says to fix C below 3 and enumerate integer-pair arcs. This paper supplies the version at ALL energies in the two-variable (tau, eta) form with the table of the families, and is the fastest way to a first set of seeds. The relation between the two is the book's eq. 4.27/4.28 equivalence. The C values here are the project's Jacobi constant (C = 2 eps' sqrt(a(1-e^2)) + 1/a is the project's Tisserand form too, and V = sqrt(3 - C) is the excess speed at the unit-circle crossing with mu = 0).
- The `#906` amendment to part (b) (do not reject small demanded turns near resonance) is supported from a different direction: this paper's A0 family contains the exact tangent ellipses (eqs. 38 to 41, eta = tau and the i, j integer points), where the two arcs meet the circle tangentially at the Moon and the demanded turn is zero for a reason (a resonance, not a missing encounter). A zero or near-zero turn is a legitimate family member there. That sharpens the Hitzl-Henon recommendation (return "indeterminate", do not reject).
- `#888`/`#890` gate: the gate's turn is computed between V-infinity vectors at an encounter. The arcs here have a turn that follows from eq. (3) of Bruno 1981 (velocities V1, V2 at collision), see the Bruno digest; this paper gives only the speed V = sqrt(3 - C), not the turn. The turn on the two sides of the collision is the angle between (V1, V2) and (-V1, V2), so cos(turn) = (V2^2 - V1^2)/V^2 with V1 = eps'' e sqrt(a) sin(eta) and V2 = eps' sqrt(a(1 - e^2)) - 1 (Bruno's eqs. 3 and 4; the in-and-out velocities differ by the sign of V1: V1' = -V1, V2' = V2). COMPUTED consequence: V1 = 0 exactly when e = 0 or sin(eta) = 0, so the demanded turn of a generating arc is zero there (the tangent ellipses of eq. 41 and the circular orbit at tau/pi = 0.5) and 180 degrees when V2 = 0. These are the cases `#906(b)` must treat as "not an encounter or resonance".

Conflicts or gaps:
- This paper contains no Jacobi-constant extremum condition and no stability result; the project must not quote any stability statement from it.
- The paper's sample conditions (r1 = 0.01659 for the Earth) conflict with Bruno's 0.01665; neither is the project's own Earth radius convention. Use the project's body radii for any gate and cite each paper's value only in a comparison test.
- The three printed errata of section 2.3 mean that any test built from Tables 3 and 5 at those rows will fail by design; use the recomputed values (marked COMPUTED) and keep the printed value as an expected-failure record.

## 5. Techniques applicable to the project's problems

Mapping per the coordinator's 2026-10-04 correction: `#890` and `#895` are the Titania-Oberon periodic orbit and its real-ephemeris arcs, not Earth-Moon solutions; only METHODS are mapped to them. Earth-Moon physics maps to `#899`, `#884`/`#905`.

### 5.1 `#899`: second-species continuation to the Earth-Moon mass (reproduction of Casoliva et al. 2010) with the Henon enumerator

What to build, using this paper:
1. A seed generator that solves eq. 30 (and its hyperbolic twin, eq. 10) for the (tau, eta) curves of each family, then returns a, e, x0, C, V from eqs. 29 and 33. The inputs are the signs (eps eps'' fixed per family: -1 for A, +1 for B, (-1)^(i+j) for C_ij) and a branch label. Use the tabulated rows (section 2) as the golden set; section 2.4 lists the independent checks.
2. Seed selection for the Earth-Moon cyclers of Casoliva et al. (class 1: p-q resonant orbits that pass near the Moon). In this paper's variables the resonant tangent ellipses of eq. 41 are at integer (tau/pi, eta/pi) = (i, j) with a = (i/j)^(2/3); the orbit is the p:q resonance with p = i, q = j (period 2 pi i / j in time units, in whole sidereal months tau/pi). Families A_i, B_i with i >= 1 pass through these points (case b of section 1.7). If a Casoliva p:q is read as p sidereal months over q revolutions (the mapping of Casoliva's p-q convention to (i, j) has to be checked in that paper before use), the tangent ellipse has a = (p/q)^(2/3): 7:3 gives 1.7584, 5:2 gives 1.8420, 3:2 gives 1.3104, 2:1 gives 1.5874 (COMPUTED arithmetic only). e = |1 - (q/p)^(2/3)| is below 1 only for q/p <= 2 sqrt(2). Which A or B branch a given resonance lies on follows the parity of i - j (eq. 40 text); the critical-orbit tables of the 1997 book (Table 4.4) give the extrema along each.
3. Because the printed tables stop at tau/pi = 4 (sidereal months), the 7:3-type seeds (tau/pi = 7 months) are beyond them. The generator must be run by solving eq. 30 beyond tau/pi = 4, using the asymptotic eq. 34/35 as an initial guess (they give eta, e, a to leading order in (pi/tau)^(1/3)).
4. The seed at mu = 1e-6 is then a second-species orbit of the circular restricted problem with the arc of this paper joined at an angular point of near-zero distance (Perko 1976, Bruno 1981 give the first-order correction). Continuing in mass to 0.01215 is the Casoliva three-step strategy. This paper gives ONLY the mu = 0 seed. The turn at the Moon (the angular point's deflection) is given in section 4.
5. Control: the printed rows of section 2 for the mu = 0 arcs, and the nine Casoliva states already in the catalogue (with the 2008 conference table of seeds at mu = 1e-6) for the continued orbits.

### 5.2 `#906`: turn-gate hardening (do not reject small turns near resonance)

- Add the exact-zero-turn locus from section 4: the arcs with V1 = 0 (e = 0, or sin(eta) = 0 i.e. eta/pi an integer). On the family A0 these are tau/pi = 0.5 (e = 0, eta/pi = 0.5 gives V1 = e sqrt(a) sin(eta) = 0) and every integer eta/pi at which the A and B families cross the integer points of eq. 40. The gate should label them "tangent resonance (zero demanded turn by construction)", not "turn feasible", and not reject.
- Add the exact 180-degree locus V2 = 0: eps' sqrt(a (1 - e^2)) = 1.
- Use eq. 43 (slopes of the two branches at a tangent point) to give the local expansion in which the demanded turn grows linearly off the tangent point, so the gate can report "within X of a tangent resonance" with a sourced slope.

### 5.3 `#890` and `#895` (Titania-Oberon periodic orbit and its real-ephemeris arcs): methods only

- The two-variable elimination in this paper (reduce a three-equation, four-unknown matching problem to one implicit equation in two numbers, then draw the solution curves) is the technique that organises candidate itineraries cleanly: characteristic curves are separated in (tau, eta) but "hopelessly entangled" in (e, a) (Figures 4 and 9). For a two-moon chain the analogue is to list the two-body arcs between consecutive encounters as curves in (time-of-flight, anomaly) and read feasible joins from their intersections, as the `#890` closure does implicitly.
- The special points (eqs. 38 to 45, double points at integer and half-integer tau/pi) are the places where a family can switch branches. For continuation (`#895`) these are the bifurcation points to detect; the paper's local expansions (eqs. 46 to 50) give the branch directions. This is a method note, not a statement about the Uranian system.

### 5.4 `#884`/`#905` (Earth-Moon cyclers under the Sun)

- The mu = 0 families with r1 = 0.01659 and V <= 4.887 (eqs. 52, 56) give the screening for which arcs can be second-species skeletons of an Earth-Moon cycler: nearly all elliptic arcs pass the speed limit and all with |x0| and |x1| above the Earth radius pass the surface limit.
- The `#905` rerun should use this paper's tangent-ellipse points (eq. 41) to label the resonance of each recovered orbit (p:q = i:j), which is a sourced and independent label to compare with the resonance Casoliva assigns.

### 5.5 `#913` to `#918`

Little. `#913` (Russell's 77 heliocentric solutions), `#914`, `#915`, `#918` concern heliocentric or moon-tour cyclers; the paper's Sun-Earth case is a ZSOI problem on a circular Earth orbit, so its tables (probe leaves the Earth and returns with no deep-space manoeuvre) are a coarse sanity control for the first leg of an Earth-return trajectory only. `#916` (persistence of stable cyclers under eccentricity and the Sun): this paper's closing remark that the eccentricity of the Earth's orbit is the most important correction is the same question at mu = 0; no numbers. `#917` (Casoliva rows in the elliptic problem): the generating arcs are the same at any eccentricity of the primary orbit only at leading order; nothing from this paper beyond the arcs.

## 6. Recommended follow-ups (task numbers not registered)

1. Implement `search/second_species_arcs.py` (name suggestive): solve eqs. 10 and 30 for any tau, branch and sign set; return a, e, x0, x1, C, V, the collision velocity components (V1, V2) of Bruno's eq. 3 and the demanded turn. Tests: the printed rows of Tables 1 to 9 at the tolerances in section 2.2 and the controls of section 2.4; expected failure (strict) records for the three printed defects of section 2.3.
2. Extend the generator past tau/pi = 4 with the asymptotic forms eqs. 34 and 35 as initial guesses, and tabulate the families to tau/pi = 8 for the 7:3 and 5:2-type seeds; compare with Table 4.4 of the 1997 book (the book gives the critical orbits of the same families).
3. Add the exact zero-turn and 180-degree loci of section 5.2 to the `#906` gate test list.
4. Read the sign columns of Tables 6 to 9 and Figure 4 labels from a second copy and add them to the data (this digest read the numbers but not every sign).
5. Resolve the letter-versus-digit reading of the leading term in eq. 35 and the coefficients of eqs. 47 to 50 from a second copy of the paper (low priority; they are not used in section 5).
6. Record the Earth radius convention difference (0.01659 here, 0.01665 in Bruno 1981) as a note in the body-radius table the gates use.
