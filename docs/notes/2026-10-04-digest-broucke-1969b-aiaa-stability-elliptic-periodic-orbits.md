# Digest: Broucke 1969, "Stability of periodic orbits in the elliptic, restricted three-body problem" (AIAA Journal)

Written 2026-10-05 (Sydney) from all seven pages, read from page images (the text layer has badly formatted numbers). READ = printed;
COMPUTED = my integration in scratch with `core.er3bp` (not committed, no test files); INFERRED = reasoning.

## 0. Citation and scope

- R. Broucke, "Stability of Periodic Orbits in the Elliptic, Restricted Three-Body Problem", AIAA Journal 7(6):1003-1009, June 1969, DOI 10.2514/3.5267.
  Jet Propulsion Laboratory. Presented as AAS Paper 68-086 (Jackson, Wyoming, 3-5 September 1968); revision submitted 16 December 1968. 7 pages (pp.1003 to 1009).
- Filed in the private paper corpus as `broucke-1969b-stability-periodic-orbits-elliptic-restricted-three-body-problem-aiaa-j-7-1003-doi-10.2514-3.5267.pdf`.
- It is the short journal companion of JPL TR 32-1360 (digested in `docs/notes/2026-10-04-digest-broucke-1969-elliptic-periodic-orbits-part-a.md` and `-part-b.md`).
  It contains NO table. Everything numerical is in the text: five initial states and a handful of stability extremes.

## 1. What it adds to or changes from the TR

1. **The stability text for six families, in one place.** TR part A and B state the region of each family in scattered paragraphs; this paper gives them
   together (p1008 to 1009): 7P region 6 (good information to e = 0.35), 7A region 1 (to e = 0.75), 8P region 6 (to e = 0.55), 8A region 4 then 2, region 1 at e = 0.85,
   11P region 4 to e = 0.45, 11A region 6 for e below 0.30 then region 3 to e = 0.80. "We have computed the stability for six families only (about 620 orbits)" (p1008),
   which no table in the TR reflects (the TR prints no per-row a1 or a2).
2. **Three extra numbers on 11P and 11A**: 11P a2 reaches about 13,400 near e = 0.18 and about 1000 at e = 0.45, a1 rises from -5000 to -200 between e = 0 and 0.45;
   11A a1 has a minimum of -5500 at e = 0.18 and a2 a minimum of "-5200" at e = 0.67 (p1009). COMPUTED check in section 4: all agree except the 11A a2 minimum.
3. **8P eccentricity range resolved.** This paper says "About 120 orbits ... with eccentricities ranging from 0.0 to 0.975" (p1009); the TR text says 0.075, a slip that part B had inferred
   from the table. Confirmed here.
4. **Circular-problem comparison stated explicitly:** in the circular problem there are only three classes (stability, even instability, odd instability), because two eigenvalues are +1
   (Eq. (a): lambda, 1/lambda, +1, +1); in the elliptic problem (Eq. (b)) lambda, 1/lambda, mu, 1/mu give one stable and six unstable classes (abstract, p1003 to 1005).
5. Claim of priority (abstract): "probably for the first time, a classification of the multipliers is given for these systems" (author's claim, not checked).
6. Not changed: the equations (Eqs. 1 to 17 equal TR Eqs. 2, 3, 12, 15, 17 in the planar case), the differential correction (Eq. 23 = TR Eq. 149), the strong periodicity criterion,
   P/A start epochs, the (a1, a2) invariants (Eqs. 37 to 45 = TR Eqs. 165 to 174), the seven regions (Fig. 1 = TR Fig. 2; Fig. 2 = TR Fig. 3). The same orbit starts are printed:
   7P, x0 = 0.15212027, ydot0 = 3.16076559, mu = 0.012155, e = 0 (p1008; TR Eq. 178); 8P, x0 = -0.4017933, ydot0 = +3.1437189, mu = 0.5, e = 0 (p1008; TR Eq. 179);
   11P, x0 = -0.07084826, ydot0 = +0.82832745, mu = 0.5, e = 0 (p1009; TR Table 16 row 1 gives -0.0708483, 0.8283275).

## 2. Equations as printed (pp1003 to 1008), compared with the TR

- Inertial barycentric primaries (Eq. 1): xi1 = -mu(cos E - e) = -mu r cos v, eta1 = -mu(1 - e^2)^(1/2) sin E; xi2 = (1-mu)(cos E - e). r = 1 - e cos E = (1 - e^2)/(1 + e cos v) (Eq. 2); a = n = 1; periapsis at t = 0.
- Rotating (Eq. 6), pulsating (Eq. 10), Lagrangians Eqs. 12 and 14, equations of motion **x'' - 2 y' = (r/p)(x - (1-mu)(x - x1)/r1^3 - mu(x - x2)/r2^3, y'' + 2 x' = (r/p)(y - (1-mu) y/r1^3 - mu y/r2^3)** (Eq. 15; primes now d/dv; p = 1 - e^2), x1 = -mu, x2 = 1 - mu (Eq. 11). Identical to TR Eq. 33 and to `core/er3bp.py`.
- Time: t_dot = r^2/p^(1/2) (Eq. 16); r from r'' = 2 r'^2/r + r(1 - r/p) = -(2/p) r^3 + (3/p) r^2 - r (Eq. 17); a seventh-order autonomous system. 15 dependent variables (18 in three dimensions; parameters P1 to P18; 32 more, P19 to P50, for the variational solutions, Eq. 22 lays them out as the matrix R).
- Variational equations Eqs. 19 to 21 (same as TR Eqs. 136 to 138); correction equations Eq. 23 for (Delta x0, Delta ydot0) with y_F = x_dot_F = 0 at k pi; "converges in fewer than about five iterations"; Lawson's least-squares routine for a small determinant.
- Fundamental matrix: A = [[0, I], [a, 2J]] (Eq. 26), J = [[0, 1], [-1, 0]] (Eq. 27), S (Eq. 28), S A^T S^-1 = -A (Eq. 29), R^-1 = S R^T S^-1 (Eq. 32); lambda = e^(alpha T), mu = e^(beta T) (Eq. 33); polynomial s^4 + a1 s^3 + a2 s^2 + a1 s + 1 = 0 (Eq. 34), reciprocity checked to 1e-10.
- Indices k1 = lambda + 1/lambda, k2 = mu + 1/mu (Eq. 35); (s - lambda)(s - 1/lambda)(s - mu)(s - 1/mu) (Eq. 36); a1 = -(k1 + k2), a2 = 2 + k1 k2 (Eq. 37); X^2 + a1 X + (a2 - 2) = 0 (Eq. 38);
  parabola a2 = a1^2/4 + 2 (Eq. 44); lines a2 = 2 a1 - 2 and a2 = -2 a1 - 2, tangent to the parabola at (+-4, 6) (Eq. 45); circular problem on a2 = -2 a1 - 2.

## 3. Disagreements and slips (relative to the TR and internal)

1. **Eq. 39 sign:** printed k = [a1 +- (a1^2 - 4 a2 + 8)^(1/2)]/2; Eq. 38 gives k = [-a1 +- ...]/2, and the TR (Eq. 167) prints -a1. The journal Eq. 39 is a sign slip (COMPUTED: roots of X^2 + a1 X + (a2 - 2) with a1 = -(k1 + k2) return k1, k2 only with -a1). A classifier written from this paper's Eq. 39 would swap k with -k, i.e. confuse region 4 with 5 and 6 with 7.
2. Eq. 36 prints "(beta - mu)" for "(s - mu)". Eq. 18 prints dy_dot/dt for dy_dot/dv. A footnote (p1005) notes the symbol mu in Eq. (b) is an eigenvalue, not the mass ratio.
3. The 8P range 0.075 in the TR text versus 0.975 here (section 1).
4. 11P: "51 orbits ... e from 0 to 0.453" here; the TR Table 16 has 52 rows (the e = 0 orbit plus 51).
5. Text p1009 on 8A: "the family starts in Region Four, passes in Region Two ... up to e = 0.85, except that the last computed orbit e = 0.85 itself had just crossed the parabola ... belongs to Region One": same as the TR (p66).
6. The abstract says "Eleven hundred periodic orbits" (the TR: 1127). 11A: "At this last value of e (0.895) there is a close approach with one of the primaries".
7. Regions are defined by Eqs. 39 to 43 and the seven numbered definitions on pp1007 to 1008 (no disagreement with the TR; listed for the classifier): region 3 has k1^2 > 4, k2^2 > 4 and opposite-sign real eigenvalue pairs, region 4 k1 > 0 and k2 > 0 (both with square above 4), region 5 both negative, region 6 k1 > 2 with k2^2 < 4, region 7 k1^2 < 4 with k2 < -2. These agree with the TR's Table 3.

## 4. Check against `core.er3bp` (COMPUTED, scratch)

mu = 0.5 for 8, 11 families, DOP853 at 1e-13; each row taken from the TR tables transcribed in part B, corrected in (x0, ydot0) so that y = x_dot = 0 at the half period, then the planar 4 x 4 block of the one-period STM, a1 = -trace, a2 = sum of principal minors, regions by k.

| Family, row (e) | a1 | a2 | k1, k2 | Region | Printed claim |
|---|---|---|---|---|---|
| 8P row 30 (0.20) | -92.2 | 93.7 | 1.006, 91.2 | 6 | region 6 |
| 8P row 60 (0.50) | -481.9 | -331.5 | -0.691, 482.6 | 6 | region 6 (to 0.55) |
| 8A row 17 (0.10) | -25.6 | 60.8 | 2.55, 23.0 | 4 | starts in region 4 |
| 8A row 37 (0.30) | -11.5 | 37.2 | complex | 2 | region 2 |
| 8A row 57 (0.50) | -4.61 | 20.3 | complex | 2 | region 2 |
| 8A row 87 (0.80) | 0.409 | 3.875 | complex | 2 | region 2 |
| 8A row 92 (0.85) | 0.759 | 2.036 | -0.708, -0.051 | 1 | region 1 at 0.85 (a1^2/4 + 2 = 2.144 > a2) |
| 11P row 1 (0.00) | -4918.9 | 9835.9 | 2, 4917 | 6 (on k = 2) | a1 about -5000 at e = 0, "even instability of the circular problem" |
| 11P row 22 (0.18) | -3120 | 13325 | 4.28, 3116 | 4 | a2 maximum about 13,400 near 0.18 |
| 11P row 34 (0.30) | -1785 | 10788 | 6.06, 1779 | 4 | region 4 |
| 11P row 49 (0.45) | -313.6 | 1568.6 | 5.08, 308.6 | 4 | a1 about -200 and a2 about 1000 at 0.45 (loosely) |
| 11A row 13 (0.10) | -5503.8 | 4859.6 | 0.883, 5503 | 6 | region 6 below 0.30 |
| 11A row 21 (0.18) | -5612.6 | -434.7 | -0.078, 5613 | 6 | a1 minimum -5500 at 0.18 |
| 11A row 28 (0.25) | -5371 | -6010 | -1.12, 5372 | 6 | region 6 |
| 11A row 43 (0.40) | -3559 | -21333 | -5.99, 3565 | 3 | region 3 from 0.30 |
| 11A e = 0.66 to 0.69 | 2696 to 3132 | -50906 (0.66), -49039 (0.69) | | 3 | a2 minimum "-5200" at 0.67 |
| 11A row 83 (0.80) | 1882 | -12006 | -1889, 6.36 | 3 | region 3 to 0.80 |

- The region statements and the 11P and 11A extremes reproduce (the 11A a1 minimum -5613 against -5500, 2 percent). **The one disagreement is the 11A a2 minimum: computed about -5.09e4 near e = 0.66 to 0.67; printed -5200.** INFERRED: a dropped digit (-52,000 would be within about 2 percent); the paper's own a1 numbers are 2 percent low in the same way. Do not use -5200 as a test value.
- Initial states: 8P (Eq. in p1008) full 2 pi closure at e = 0: 1.0e-4; 11P: 4.0e-5 (printed seven and eight digits, large multipliers). 7P/8P/11P half-revolution results are those of the TR check (part A section 8; part B).
- 8A row 92 (e = 0.85) being region 1 is robust to the 1e-3 level: k1 = -0.708, k2 = -0.051 are far from the boundary values, and a2 lies below the parabola by 0.11.

## 5. Test-ready numbers (all printed, with page; none is a table)

| Quantity | Value | Page | Tolerance advice |
|---|---|---|---|
| 7P start | x0 = 0.15212027, ydot0 = 3.16076559, mu = 0.012155, e = 0 | 1008 | full 2 pi closure 2e-5 (COMPUTED in part A) |
| 8P start | x0 = -0.4017933, ydot0 = 3.1437189, mu = 0.5, e = 0 | 1008 | closure 1e-4 over 2 pi |
| 11P start | x0 = -0.07084826, ydot0 = 0.82832745, mu = 0.5, e = 0 | 1009 | closure 4e-5 over 2 pi |
| Regions | 7P 6, 7A 1 (below 0.75), 8P 6 (to 0.55), 8A 4 then 2 then 1 at 0.85, 11P 4 (to 0.453), 11A 6 (below 0.30) then 3 (to 0.80) | 1008-1009 | region decided from re-corrected orbit; see section 6 |
| 11P a2 maximum | about 13,400 at e about 0.18 (computed 13,325 at 0.18) | 1009 | 2 percent |
| 11P a1 | -5000 (e = 0) to -200 (e = 0.45) (computed -4919, -314) | 1009 | loose, "about" |
| 11A a1 minimum | -5500 at e = 0.18 (computed -5613) | 1009 | 3 percent |
| Convergence limits | 11P last orbit e = 0.453 (10 iterations); e = 0.454 diverged | 1009 | text only |
| Orbit counts | 11P 51 orbits to 0.453; 8A about 100 to 0.85; 8P about 120 to 0.975; 11A about 100 to 0.895; 7P about 130 to 0.50; 7A about 120 to 0.99 | 1008-1009 | TR tables have 52, 92, 118, 93, 131, 120 |

## 6. Techniques applicable to the project's problems

### `#931` (stability classifier, tolerance at k = +-2, seven regions)

- **Definition to implement:** a1 = -trace of the planar 4 x 4 monodromy, a2 = sum of the six principal 2 x 2 minors, Delta = a1^2 - 4 a2 + 8, k = [-a1 +- sqrt(Delta)]/2 (the TR sign; the journal Eq. 39 has a sign slip, section 3), and regions 1 to 7 from Delta and the sizes and signs of k1, k2 (section 3 item 7, Fig. 1). The two Broucke texts and Hadjidemetriou 1975b share (a1, a2) = (alpha, beta), b = -k.
- **Tolerance at k = +-2.** Broucke's own computed families sit on the lines a2 = -2 a1 - 2 (7P, 7A, and every circular-problem orbit; 11P at e = 0 gives k1 = 2.000 exactly). On those lines the region (1 versus 6 for k near +2; 3 versus 7 and 6) is decided in the fourth decimal of k, and part A showed the region flips between the printed seven-digit state and a re-corrected one in 9 of 15 rows. So: compute from a re-corrected orbit, treat |k - 2| below a tolerance (about 1e-3 for the 7A and 7P rows) as "on the boundary" and report it as such; do not label an orbit region 6 or region 1 unless |k1 - 2| exceeds the tolerance; at the exact circular limit (e = 0) the orbit is always on the line.
- **Large-index regimes need relative tolerances.** In regions 3, 4, 6 of the 8 and 11 families |a1| reaches 5600 and a2 5e4 (computed), so a1 and a2 as test numbers need relative, not absolute, tolerances; the region is robust (k1 and k2 are far from +-2 and the sign pattern is stable) but a1 and a2 vary by 2 percent between printed-state and corrected-state values.
- **Test spectra:** a symplectic 4 x 4 matrix per region built from the Fig. 2 root patterns; plus the sourced sequence 8A: region 4 (e = 0.10) to region 2 (e = 0.30 to 0.80) to region 1 (e = 0.85), 11A: region 6 (e = 0.10 to 0.25) to region 3 (e = 0.40 to 0.80), 7A: region 1 (0.035 to 0.93) from the TR. These transitions cross the Delta = 0 parabola (4 to 2, 2 to 1) and the k = -2 line (6 to 3), which gives a sourced example of each transition type that `#931` item (d) wants logged.
- Do not use the printed -5200 for the 11A a2 minimum (section 4).

### `#933` (the Broucke test module)

- This paper adds no table rows, so the module's data come from the TR parts A and B; what this paper supplies is the stability layer for it: a set of (family, e, region) assertions (section 5) to be checked with the classifier above on re-corrected orbits, using the rows named in section 4 (8P rows 30 and 60, 8A rows 17, 37, 57, 87, 92, 11P rows 1, 22, 34, 49, 11A rows 13, 21, 28, 43, 83), each from a printed region statement.
- Record as known text defects, not as test failures: the 11A a2 minimum (-5200 printed, about -5.1e4 computed), the Eq. 39 sign, the 8P range 0.075 in the TR, 11P 51 versus 52 orbits.
- The text-level region claims are computed to hold at the sampled rows; the exact eccentricities of the region changes (8A 4 to 2 near e = 0.3, 11A 6 to 3 near e = 0.30, 7P 6 to 3 between 0.44 and 0.50 from part A) were not located and should not be asserted as tests beyond the sampled ones.

## 7. Follow-ups (not registered)

1. Locate the region change eccentricities (8A 4 to 2, 11A 6 to 3, 7P 6 to 3) by continuation and compare with 0.30 and 0.30 as printed.
2. Confirm the 11A a2 minimum with a finer scan (about -5.1e4 near 0.66 here).
3. The planar classifier of `#931` with the tolerance logic above, tested on the sequence in section 6.
