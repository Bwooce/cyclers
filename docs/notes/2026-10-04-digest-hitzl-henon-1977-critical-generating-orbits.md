# Digest: Hitzl & Henon 1977, "Critical generating orbits for second species periodic solutions of the restricted problem"

Celestial Mechanics 15:421-452 (1977), DOI 10.1007/BF01228610, received 3 November 1975. Filed in the private paper corpus as
`hitzl-henon-1977-critical-generating-orbits-second-species-periodic-solutions-restricted-problem-celest-mech-15-421-doi-10.1007-BF01228610.pdf`
(32 PDF pages; PDF page n is journal page 420 + n; text layer present but equations and the sign columns of Table I are garbled in it, so every formula and
table below was read from the page images).
Evidence tags: READ (page) means read on the page image; COMPUTED means my own arithmetic on 2026-10-04 (60-digit arithmetic, scratch scripts not kept in the repository);
INFERRED means my reading across sources. Page numbers are journal pages. Companion digest of the book that carries the same orbits as a general theory:
`docs/notes/2026-10-04-digest-henon-1997-generating-families.md`.

## 0. What this paper is, in six lines

1. Scope: planar restricted problem in the limit of mass ratio zero. The "generating orbits" (Henon's "consecutive collision orbits") are single Kepler arcs that start
   at the secondary and return to it, in the rotating frame, after one period. They form an infinite set of one-parameter families, and they are "the true limit, for mu > 0, of second species periodic solutions" (abstract).
2. A generating orbit is "critical" when the Jacobi constant C is extremal along its family. Along a family of periodic orbits an extremal C forces stability index k = +1 (or -1), so the orbit sits on the boundary of stability.
3. Mechanism: orbits with consecutive collisions are "infinitely unstable" (k = plus or minus infinity). At a critical generating orbit k jumps from +infinity to -infinity, so for small mu greater than 0 the neighbouring orbits pass through the stable range -1 to +1 over a short interval: "neighboring orbits will then have a finite (but small) region of stability" (abstract).
4. Method: the timing condition F0 = 0 (Henon 1968) in two variables (tau, eta) plus one new analytic condition G* = 0 (their eq. 27). Critical generating orbits are the intersections of the two loci.
5. Result: 179 critical generating orbits, with tau/pi and eta/pi to 10 decimals in Table I (pp.438-442) for 0 to 10 in both variables; ten of them drawn and tabulated in Table II (p.450).
6. What it does not do: no orbit at positive mass is computed. Stability near the critical orbits is promised for a second paper (section 8).

## 1. Definitions, exactly as the paper sets them

### 1.1 Generating orbit and species (READ, pp.421-423, 427)

- Abstract (READ p.421): "The second species periodic solutions of the restricted three body problem are investigated in the limiting case of mu = 0. These orbits, called consecutive collision orbits by Henon and generating orbits by Perko, form an infinite number of continuous one-parameter families and are the true limit, for mu > 0, of second species periodic solutions for mu > 0." (The printed text says "for mu > 0" twice; the second is surely mu tending to 0. INFERRED typographical slip; the body text on p.427 reads: "consecutive collision orbits are the natural limit of second species periodic orbits of the restricted problem as mu tends to 0".)
- The limit is odd, and the paper says so (READ p.427): "in the limit mu = 0, one does not fall back to periodic solutions of the first or second kind of the two-body problem. For this particular type of orbit, as mu tends to zero, the distance of closest approach d to the second body also goes to zero in such a way that, for mu = 0, a finite value of mu/d is obtained. This is somewhat paradoxical for it means that, even when the second body has no mass (mu = 0), the effect of the second body is not negligible."
- Second species is cited to Poincare 1892, vol. III, chapter 32. The paper states the proof situation (READ p.423): "an existence proof of the type previously given for periodic solutions of the second kind (Arenstorf, 1963a; Barrar, 1965) has not as yet been developed for the second species periodic solutions"; an existence proof "based upon approximating expansions in mu of order mu^m (m at most 2) together with error estimates" is Perko 1974.
- The hyperbolic case is dispatched (READ p.424): the hyperbolic consecutive-collision orbits "form (for mu = 0) a single family (Henon, 1968, Table I) along which there is no extremum of the Jacobi constant C. Therefore, we shall consider exclusively the case of elliptic consecutive collision orbits."

### 1.2 Geometry and the three switches (READ pp.423-425, eqs. 1-5, Fig. 2)

Earth (mass 1 - mu) at the origin of a non-rotating frame, the Moon on the unit circle with unit angular rate (period 2 pi), direct sense. At t = 0 the Moon is at R = (1, 0) in the rotating frame sense and the particle is at perigee (S') or apogee (S'') of its Kepler ellipse. The collisions occur at P at t = -tau and at Q at t = +tau, so the interval between collisions is 2 tau.
Three switches (eq. 1): sigma0 = +1 if perigee is at x > 0, -1 if x < 0; sigma1 = +1 if the orbit is direct, -1 if retrograde; sigma2 = +1 if the particle is at perigee at t = 0, -1 if at apogee. Combined switch sigma = sigma0 sigma2 (eq. 8), and sigma_eta = sgn(sin eta).
Fixed-axes position (eqs. 2-4): the Moon is at (cos t, sin t); the particle is at x = sigma0 a (sigma2 cos E - e), y = sigma0 sigma1 sigma2 a sqrt(1 - e^2) sin E, with Kepler's equation t = a^(3/2) (E - sigma2 e sin E). E = 0 is perigee when sigma2 = +1 and apogee when sigma2 = -1.
Collision at Q (t = tau, E = eta), eq. 5 with c and s for cosine and sine: c_tau = sigma0 a (sigma2 c_eta - e), s_tau = sigma0 sigma1 sigma2 a sqrt(1 - e^2) s_eta, tau = a^(3/2) (eta - sigma2 e s_eta). The same system holds at P with t = -tau, E = -eta. "we have a system of three equations in the four unknowns tau, eta, a, e and thus an infinity of solutions" (READ p.425).
Elimination (eqs. 6-10, READ p.425): a = (1 - sigma c_tau c_eta) / s_eta^2, e = (sigma2 c_eta - sigma0 c_tau) / (1 - sigma c_tau c_eta), rho^2 = 1 - sigma c_tau c_eta = a s_eta^2, rho = sqrt(a) s_eta = sigma_eta sqrt(1 - sigma c_tau c_eta).
Timing (periodicity) condition (eq. 11, READ p.425): F0(tau, eta; sigma, sigma_eta) = rho [eta rho^2 - s_eta (c_eta - sigma c_tau)] - tau s_eta^3 = 0. "This equation was derived originally by one of us (Henon, 1968)." Solutions exist for all tau greater than 0.515 00... and eta greater than 0 (p.425).
Because the orbit leaves and returns to the Moon in the rotating frame after 2 tau, the period of the periodic orbit in the rotating frame is T = 2 tau (Table II, column T: A0(0) has tau/pi = 0.5 and T = 3.141 59; READ p.450). INFERRED from that table: one close approach per period, an "angular point" at the Moon: "All orbits are symmetrical with respect to the x axis and have one orthogonal crossing of this axis and one intersection with the Moon at a so-called 'angular point'. (The angle may be zero in particular cases.)" (READ p.445). That last sentence is the only statement in the paper about the turn angle at the collision; see section 4(b).

### 1.3 Characteristics and the naming of families (READ pp.425-426, Fig. 3a-3b)

- "We consider the characteristic of a family of symmetric periodic orbits to be the appropriate locus of initial conditions x0, ydot0 and C plotted in either the x0, C plane or the x0, y0 plane. Consequently, since solutions tau and eta of (11) specify both the orbit (a and e) and the initial conditions (x0, y0, and C), it seems natural to label the solution curves in Figure 3a the characteristics for consecutive collision orbits." (READ p.425.)
- In the (tau, eta) plane the characteristics are "neatly separated and form easily recognizable patterns. On the other hand, in the a, e plane the loci are hopelessly entangled" (READ p.425; Henon 1968 Figure 9 for the latter).
- Families are labelled A0, A1, A2, ...; B1, B2, ...; C12, C23, ... (READ p.425, Fig. 3a). From the figures and Table I (INFERRED from the table content): A_i are the families whose characteristic passes through (tau/pi, eta/pi) = (i + 1/2, i + 1/2) (the trivial critical orbit A_i(0)); B_i pass through the line tau = eta and eta near i; C_ij have critical orbits near (tau/pi, eta/pi) = (i, j) with i smaller than j. Figure 3a plots tau/pi to 8 and eta/pi to 6 roughly.
- Critical orbit naming (READ pp.436-437): "the origin of a numbering system l = 0, +1, +2, ... for the critical orbits A_i(l) which extends to infinity in both directions. For the direction of numbering, it seems natural to take l greater than 0 for the upper branch of a family A_i (i.e., the branch with larger values of eta) and l less than 0 for the lower branch." For B_i: "increasing positive values of l are attained as we move to the right (tau increasing) along the characteristic B_i from the origin B_i(0)." For C_ij only two critical orbits exist per family, labelled (1) and (2) (READ p.442, below Table I): "The first, C_ij(1), is located just below the self-intersection of the characteristic for family C_ij on the branch with positive slope. The second, C_ij(2), is located to the right of the self-intersection on the branch with negative slope." C_ij(1) is "in all cases a direct orbit" and C_ij(2) "always retrograde" (READ p.451).

### 1.4 Jacobi constant on the generating orbits (READ pp.427-429, eqs. 12-22)

- Eq. 12: x^2 + y^2 + 2(1 - mu)/r1 + 2 mu/r2 - v^2 = C. Eq. 14-15: p_x = xdot - y, p_y = ydot + x, C = 2 (p_y x - p_x y) - 2 (T + V) = 2 (h_z - E).
- Two-body limit (eqs. 16-17): E = -1/(2a), h_z = sigma1 sqrt(a (1 - e^2)), so C = 2 sigma1 sqrt(a (1 - e^2)) + 1/a. This is the Tisserand relation and holds "for any elliptical orbit (second species or not)" (READ p.428).
- On a consecutive-collision orbit (eq. 21, READ p.428), the "remarkably simple form": C(tau, eta; sigma, sigma_eta) = 2 sigma s_tau / rho + s_eta^2 / rho^2.
- Periodicity (eq. 22, READ p.429): C is periodic in tau and eta, and invariant under tau goes to tau + pi with sigma going to -sigma, or eta goes to eta + pi with sigma going to -sigma (then sigma_eta flips automatically). So only 0 to 2 pi in each variable and one sign of sigma need be computed.
- Surface plots of C(tau, eta): Figs. 4a-4c, slices Figs. 5a-5d; loci dC/deta = 0 and dC/dtau = 0 (eqs. 23-24, called F1 = 0 and F2 = 0) in Fig. 6; the F2 = 0 locus "is the subset of the F1 = 0 locus composed of the diagonal straight lines" (READ p.434, Fig. 6 caption). Eqs. 23-24 themselves are not transcribed here; they follow by differentiating eq. 21.
- A caution the paper gives (READ p.429): eq. 21 "is meaningful only for a second species solution (i.e., only along a characteristic in Figure 3)", because it uses the collision equations (5).

## 2. What "critical" means and why it matters (READ pp.421-423, 434-435, 451)

- Definition (READ p.421): "The orbits studied here are further restricted to so-called critical orbits which are located on the boundary of stability with the stability index k = +-1 (Henon and Guyot, 1970; Hitzl, 1975a)." The stable range is -1 to +1: "the transitions across the stable region -1 <= k <= 1" (p.422).
- Why generating orbits matter (READ pp.421-422): "In general, orbits with consecutive collisions are 'infinitely unstable' with a stability index k = +- infinity. However, following a family of consecutive collision orbits, we find that the Jacobi constant C varies continuously and occasionally passes through extremal values. Furthermore, previous work on critical orbits for the restricted problem (Henon, 1965) showed that, along a family of periodic orbits, an orbit with an extremal value in C must be a critical orbit with k = +1. Hence, by determining consecutive collision orbits with extrema in C for the limiting case of mu = 0, isolated orbits possessing jumps in the stability index k from +infinity to -infinity are obtained. These critical generating orbits are thus of great interest since, for small mu > 0, k will only reach very large values (rather than infinite values) and the transitions across the stable region -1 <= k <= 1 will now occur over a finite (but small) interval. As a result two critical orbits with k = +1 and k = -1, respectively, will be obtained which separate a small interval of stable orbits from the remaining orbits which are generally strongly unstable."
- Figure 1 (READ p.422): stability index k against C along a family b of periodic orbits at mu = 1/2, with critical orbits labelled b1, b2, ..., b28 (Henon 1965): the picture of k blowing up and returning through -1 to +1.
- Expectation (READ pp.422-423): "a consecutive collision orbit possessing an extremum in C will be the common termination of two families of critical orbits (Henon and Guyot, 1970)."
- "Criticality" is therefore a statement about the periodic-orbit families that end at the generating orbit as mu tends to 0, not a failure of continuation. The paper does not say that continuation from mu = 0 is non-unique or fails at a critical generating orbit. What it shows is that at such an orbit two critical-orbit families terminate, k flips sign through infinity, and a small stable window opens at small positive mu. INFERRED reading; no sentence of the paper addresses non-uniqueness of the continuation.

### The condition that identifies one (READ pp.434-436, eqs. 25-28)

Critical generating orbits are the points of the (tau, eta) plane where F0 = 0 (on a characteristic) and C is extremal along it. By Lagrange multiplier (eq. 25: d/dtau and d/deta of [C + lambda F0] equal zero) and elimination of lambda (eq. 26): G = C_tau F0_eta - C_eta F0_tau = 0, i.e. the gradients of C and F0 are parallel. Worked out (eq. 27, READ p.435):

G*(tau, eta; sigma, sigma_eta) = (rho^4 / s_eta^2) sigma G = rho [T1 c_eta^2 + T2 c_eta + T3] + c_tau c_eta^4 - 2 sigma c_eta^3 + 2 c_tau c_eta^2 - sigma c_tau^2 (3 - c_tau^2) c_eta + c_tau (2 - c_tau^2) = 0,

with the secular terms (eq. 28) T1 = sigma (3 tau + s_2tau), T2 = -6 tau c_tau, T3 = sigma (3 tau c_tau^2 - s_2tau), where s_2tau = sin 2 tau. The square bracket can be written sigma [3 tau (c_eta - sigma c_tau)^2 - s_2tau s_eta^2] (READ p.435). G* is periodic in eta and invariant under eta goes to eta + pi, sigma going to -sigma, sigma_eta going to -sigma_eta (eq. 29). G* = 0 was solved numerically for sigma = -1 over 0 <= tau/pi <= 8, 0 <= eta/pi <= 2 (Figs. 7a, 7b); "the diagonal straight lines belonging to both F1 = 0 and F2 = 0 are indeed solutions"; the critical orbits are the intersections of this locus with the characteristics F0 = 0 (READ p.436). A related equation in terms of a and e was obtained by Bruno 1973 (READ p.435; not held).
Cross-validation (READ p.451, section 8): G* = 0 "has also been obtained independently by both of us from a direct stability analysis of the second species orbits in the limiting case of mu = 0. In particular, the condition G* = 0 comes directly out of first order asymptotic matching (Breakwell and Perko, 1966) as a necessary condition for the trace of a 4 x 4 'stability matrix' M to have a finite value in the limiting case of mu tends to 0." The paper lists three consequences: it confirms the prediction that extrema of C are critical orbits, it is a strong check on the computations, and it leads to further consequences of asymptotic matching.

### Trivial and spurious critical orbits (READ pp.436-437)

- Two infinite trivial sets. First, A_i(0) at tau/pi = eta/pi = i + 1/2 with sigma = -1: "simply retrograde circular orbits travelling in the opposite sense to that of the Moon" (eq. 30). COMPUTED: eqs. 7 and 17 at these points give a = 1, e = 0, C = -1 for every i. Second, B_i(0): the degenerate case at the intersection of B_i with the line tau = eta, where the particle "describes exactly the same orbit as the Moon and is coincident with it for all time". COMPUTED: a = 1, e = 0, sigma1 = +1, C = 3.
- The "something quite subtle" (READ p.437): points tau/pi = i, eta/pi = j (integers, i not equal j, sigma = (-1)^(i+j), eq. 31) satisfy both F0 = 0 and G* = 0 but "are spurious solutions": two branches of F0 = 0 meet there, F0_tau = F0_eta = 0, and eq. 26 is satisfied for any function C. "Nature has in fact played a nasty trick upon us here, because it so happens that there is an actual extremum of C very close to each (i, j) point on the branch with positive slope." Examples: the third extremum of A0 is at tau/pi = 1.997 356, eta/pi = 0.998 207, not at (2, 1). These near-integer extrema are the bulk of Table I. COMPUTED consequence for anyone coding this: Newton from a start more than about 1e-3 away converges to the spurious integer point (A5(1), section 6).
- Existence of the C_ij critical orbits (eqs. 45-47, READ p.450 image). For C_ij(1), tau/pi and eta/pi are close to i and j; from eq. 4: a is approximately (tau/eta)^(2/3), approximately (i/j)^(2/3), less than 1 (eq. 45); from eq. 6: e is approximately |1 - 1/a|, approximately |1 - (j/i)^(2/3)| (eq. 46); and "the condition e < 1 for elliptic generating orbits gives j < 2 sqrt(2) i together with the lower bound 1/2 < a < 1 for the semi-major axis a" (eq. 47). "Consequently, for j/i close to 2 sqrt(2), these orbits have values of a approaching 1/2 from above together with values of e approaching 1 from below." The Figure 16 caption (p.449): C38(1) "Here i/j = 0.375 is quite close to the limiting value i/j > sqrt(2)/4 = 0.353 55... given by (46)." COMPUTED: the Table I C_ij list stops exactly at j less than 2 sqrt(2) i (or at j = 10): C12, C23 to C25, C34 to C38, C45 to C4,10, C56 to C5,10, ...; 30 families, 60 orbits.

## 3. Method of computation (READ pp.442-445, section 6, eqs. 32-42)

Table I took "a battery of methods" because "no single method gives good results in all cases" (READ p.442).
(A) Near-integer orbits (most of the table). Write tau = i pi + X, eta = j pi + Y (eq. 32). Then "all terms of order 0, 1, and 2 vanish" in F0 and "all terms of order 0, 1, 2, 3 vanish" in G*, so the first significant terms are of order 3 and 4 and direct evaluation loses all digits (rho^2 = 1 - 1 in the limit). Change variables to cos X = 1 - u, cos Y = 1 - v (eq. 33); eqs. 34-35 are then of order 3/2 and 2 in u, v and free of the cancellation (READ p.443). Inverse (eq. 36): tau = i pi + (-1)^j sigma1 sigma_eta arctan( sqrt(2u - u^2) / (1 - u) ), eta = j pi + (-1)^j sigma_eta arctan( sqrt(2v - v^2) / (1 - v) ). Starting values from series (eqs. 37-38): with a0 = (i/j)^(2/3), eta approximately j pi - [sqrt(a0) - sigma1 sqrt(2 a0 - 1)]^2 / (3 i pi sqrt(a0) (a0 - 1)^2) and tau approximately i pi - sigma1 sqrt(2 a0 - 1) times the same quantity. For A0(1) they give tau/pi = 1.997 356 89, eta/pi = 0.998 207 72 against the exact 1.997 355 654 8 and 0.998 206 723 0. Accuracy of derived quantities "at least 1e-12". Final residuals of F0 and G* are 1e-14 to 1e-15 (p.438).
(B) B_i(0), i = 1 to 9: solve tau = eta and F0_tau = 0 (eq. 39), which reduces to tan tau = (3/4) tau (eq. 40; "The first 10 solutions of this equation had been previously computed (Henon, 1969)"), with sigma1 = +1, sigma_eta = (-1)^i, sigma = +1.
(C) A_i(-1) and C_{i,i+1}(2), i = 1 to 9: extremely close to each other (Fig. 3b), need very good starts, from an expansion about tau/pi = eta/pi = i + 1/2: tau0 = (i + 1/2) pi (eq. 41), tau approximately tau0 - 4/(3 tau0) -+ (4/(3 tau0))^2 / sqrt(3), eta approximately tau0 +- (4/(3 tau0))^2 / (2 sqrt(3)) (eq. 42; upper sign C_{i,i+1}(2), lower sign A_i(-1)). For i = 9: A9(-1) gives 9.486 146, 9.499 817 and C9,10(2) gives 9.485 413, 9.500 183 (to be compared with Table I: 9.486 121 529 8, 9.499 824 169 7 and 9.485 385 209 0, 9.500 192 484 1).
(D) A0(-1), B_i(-1) (i = 1 to 9), A_i(-2) (i = 2 to 8), B6(1), B7(1), C_{i,i+2}(2) (i = 3 to 8), C6,9(2), C7,10(2): starts from a tabulation of the family; (11) and (27) solved directly.
Check on completeness (READ p.445): "the alternance of maxima and minima in C provides a nice additional check that no critical orbits have been missed. Also we note that for families A_i and C_ij, even values of l correspond to minima in C and odd values to maxima; for families B_i, this is reversed." COMPUTED: this alternation holds for every family in Table I using C from eq. 17 on the printed (tau, eta) (section 6). The paper also says: "We cannot offer a mathematical proof that our list is exhaustive as, in fact, the set of all critical orbits is considered an unsolved theoretical problem (see Bruno, 1973). However, we feel quite confident that our list is complete" (READ p.438).

## 4. Tables

### 4.1 Table II (READ p.450): critical generating orbits drawn in Figs. 8-17, mu = 0

Columns: name, figure, n* (half the number of intersections of the orbit with the x axis in one period: n* = 1 simple-periodic, 2 double-periodic, ...), a, e, x0, x1, C, T. x0 = sigma0 a (sigma2 - e) at E = 0 (eq. 43) and x1 = -sigma0 a (sigma2 + e) at E = pi (eq. 44); "the particle only reaches this location if eta/pi is greater than 1". T = 2 tau. The last column of the second table below is my reproduction from the printed tau, eta of Table I.

| Name | Fig. | n* | a | e | x0 | x1 | C | T |
|---|---|---|---|---|---|---|---|---|
| A0(-1) | 8 | 1 | 1.410 19 | 0.884 45 | -0.162 94 | 2.657 43 | -0.399 13 | 1.341 13 |
| A0(0) | 9 | 1 | 1.0 | 0.0 | -1.0 | 1.0 | -1.0 | 3.141 59 |
| A0(1) | 10 | 2 | 1.587 20 | 0.369 96 | -2.174 40 | 0.999 99 | 2.970 94 | 12.549 76 |
| A1(-1) | 11 | 3 | 1.005 21 | 0.268 16 | -1.274 76 | 0.735 65 | -0.936 94 | 8.919 15 |
| B1(-1) | 12 | 2 | 1.366 17 | 0.370 34 | 1.872 12 | -0.860 23 | -1.439 48 | 8.417 07 |
| C12(1) | 13 | 3 | 0.632 89 | 0.581 82 | -1.001 12 | 0.264 66 | 2.874 12 | 6.203 10 |
| C25(1) | 14 | 4 | 0.544 59 | 0.838 13 | -0.088 15 | 1.001 03 | 2.641 31 | 12.526 52 |
| C25(2) | 15 | 7 | 0.565 60 | 0.823 20 | -0.100 00 | 1.031 20 | 0.914 06 | 12.798 04 |
| C38(1) | 16 | 7 | 0.521 11 | 0.920 48 | -1.000 78 | 0.041 44 | 2.483 18 | 18.826 27 |
| C38(2) | 17 | 11 | 0.526 36 | 0.915 26 | -1.008 12 | 0.044 60 | 1.315 28 | 18.927 11 |

(The text layer has "0.91.4 06" for C25(2)'s C and "2.495" for B2(0)'s tau/pi; both are OCR slips, the image shows 0.914 06 and 2.445 298 131 4.)
COMPUTED reproduction: from the printed Table I (tau, eta) and signs, eqs. 7, 17, 43, 44 return a, e, x0, x1, C and T = 2 tau equal to every entry above to the printed digits (for example A0(1): a = 1.58720, e = 0.36996, x0 = -2.17440, x1 = 0.99999, C = 2.97094, T = 12.54976; C38(2): a = 0.52636, e = 0.91526, C = 1.31528, T = 18.92711). The check needs sigma1 = -1 for A0(-1) and A1(-1) (printed, Table I), which breaks the otherwise alternating pattern of sigma1 along A_i.
Cross-source check (COMPUTED): Henon's 1997 book re-tabulates the critical arcs in new variables with the old names (Table 4.4, pp.71-72) and prints C to six decimals. My C from the printed (tau, eta) of Table I agrees with the book's printed C for all 45 rows that I compared (A0 to A4 and B1 to B4 members with |l| up to 3, and C12 to C45 members 1 and 2), the largest difference being 5e-6, so Table I and the book are consistent with each other and with eq. 17 above.
Also COMPUTED: the relative speed at the collision satisfies v^2 = 3 - C exactly for the incoming and the outgoing rotating-frame vectors (section 5(b)); this is the zero-mass form of the Jacobi integral at the secondary and is the identity behind the C less than 3 requirement.

### 4.2 Table I (READ pp.438-442): the 179 critical orbits with 0 to 10 in tau/pi and eta/pi

Transcribed in the appendix (section 8) with the printed tau/pi and eta/pi. The sign columns sigma0, sigma1, sigma2 are given there as three signs. They follow a pattern that I read off all five page images and that the appendix reproduces: for A_i with i even, l less than 0 has (sigma0, sigma2) = (-, +) and l greater than 0 has (+, -); for i odd it is the reverse; for B_i, sigma0 = sigma2 = (+ for i even, - for i odd); sigma1 for B_i is + for even l and - for odd l; for A_i, sigma1 is + for odd l except that A_i(-1) has sigma1 = - and A_i(0) has sigma1 = -; for C_ij, sigma0 is + for i odd and - for i even, sigma2 is + for j odd and - for j even, sigma1 = + for (1) and - for (2). So sigma = sigma0 sigma2 is -1 on every A orbit, +1 on every B orbit, and for C_ij it is (-1)^(i+j) (COMPUTED fit, equal to the pattern of eq. 31).
Counts (READ and COMPUTED): 179 critical orbits in all. Of these, 10 are the trivial A_i(0) (i = 0 to 9) and 9 are the degenerate B_i(0) (i = 1 to 9), leaving 160 nondegenerate orbits; the table therefore has 169 rows with decimal values (160 plus the nine B_i(0)) and 10 rows A_i(0) with the half-integer values. (Figure 3a's caption counts "the approximate location of 62 critical orbits" in the plotted region only.)

## 5. Techniques applicable to the project

### (a) Enumeration of generating orbits and families (the part of the paper that is directly implementable)

Everything needed is eqs. 5-11 and 27-28 plus the switch table. Recipe, precise enough to code:
1. Unknowns (tau, eta) in R^2; fixed switch sigma in {+1, -1} (eq. 22 lets one fix sigma = -1 for the tau/pi, eta/pi range, with the alternate obtained by shifting eta by pi). sigma_eta = sgn sin eta.
2. For each (tau, eta) compute rho^2 = 1 - sigma c_tau c_eta and, when rho^2 is greater than 0 (and a, e acceptable: a greater than 0, 0 less than e less than 1), a and e from eq. 7.
3. Characteristics are the zero set of F0 (eq. 11) in the (tau, eta) plane; contour or continue it. Solutions exist for tau/pi greater than 0.515 (p.425). For each admissible pair of the remaining switches the pair (a, e) gives x0, xdot0, C (eq. 17 or 21), the Moon-frame relative speed v^2 = 3 - C at the collision, and the period T = 2 tau.
4. Classify: the family labels (A_i, B_i, C_ij) are determined by the characteristic's topology (Fig. 3a): where it crosses tau = eta (B), through (i + 1/2, i + 1/2) (A), or the self-intersecting loops with ends near the integer points (C).
5. Critical orbits: intersect F0 = 0 with G* = 0 (eqs. 27-28) with the u, v substitution (eqs. 33-36) near integer points.
Relation to what `#899` needs. The P2 seeds of the synthesis come from Barrabes & Gomez's matched in/out maps at small mu, for a chosen (p, q) resonance and Jacobi constant. Hitzl and Henon enumerate the zero-mass skeleton, which has no mass dependence at all: one (tau, eta) pair fixes (a, e, C, T). The two descriptions meet at the resonance points. COMPUTED: the Table I orbit near (tau/pi, eta/pi) = (i, j) has a approximately (i/j)^(2/3) and T = 2 tau approximately 2 pi i, so in Barrabes & Gomez's labelling (p revolutions of the particle, q revolutions of the Moon, a = (q/p)^(2/3), period 2 pi q) it is the resonance (p, q) = (j, i). Example: A0(1) is the Casoliva "1-2" resonance (a = 1.5874 = 2^(2/3), T = 12.55 close to 4 pi); C12 is "2-1" (T = 6.20 close to 2 pi); C37 is "7-3". INFERRED mapping, checked by a = (i/j)^(2/3) only. A second INFERRED cross-source link: the bound j less than 2 sqrt(2) i in eqs. 45-47 is the same condition under which Barrabes & Gomez's C_J interval (their eq. 10, C_J1 = (p/q)^(2/3) - 2 sqrt(2 - (p/q)^(2/3)) and C_J2 with the plus sign) exists, because it is e less than 1 at the resonance; I did not re-derive it.
How Henon's enumeration helps `#899`: for every p-q resonance it gives, in closed form, the C interval over which a family of generating orbits exists. COMPUTED: the C_ij characteristic has exactly two critical orbits, with C = C_ij(1) maximum and C_ij(2) minimum, so (INFERRED) the zero-mass generating orbits of that family cover C between them. Examples from the appendix: C12 from -0.870 66 to 2.874 12; C23 from -0.963 28 to 2.971 30; C37 from 0.709 63 to 2.741 53; C25 from 0.914 06 to 2.641 31. Casoliva's Table 3 values in `search/earth_moon_resonant_families.py` fall inside: 2-1 rows C = 0.4887, 1.1964, 1.7352, 1.9522 inside the C12 range; 3-2 rows 0.1259, 0.7089, 1.6506 inside C23; 7-3 rows 1.0216, 1.0688, 1.2892 inside C37. This is a necessary-condition check (a seed (p, q, C) outside the interval has no generating orbit in the C_ij family), not a sufficient one, and it is INFERRED that the C range of a C_ij family is the closed interval between its two critical values.
Hard-won tips: use the u, v variables near integer points; the characteristics are degenerate at integer points; "no single method gives good results in all cases" (p.442).

### (b) The test for critical generating orbits, and the demanded-turn gate

The test: solve F0 = 0 and G* = 0 (section 2). To use it on a candidate you only need (tau, eta, sigma): compute F0 and G* and check both are small relative to the sensitivity there (they are third and fourth order in the distance to an integer point, so a residual near 1e-9 is normal at 10-decimal rounding).
What this says about the gate. Nothing in this paper is about the turn angle, so the relation below is INFERRED and partly COMPUTED.
- Frame. The two velocities at a junction belong to different epochs (t = +tau arriving, t = -tau leaving, which are the same instant of the rotating frame, one period apart). The gate's input must therefore be the body-relative velocity in the rotating frame (the Moon at rest), obtained from the fixed-axes velocities by rotating by -t. This is the same point as `#906`(a) and the correction in MacKay 2005 (direction change measured in the rotating frame). COMPUTED recipe: u_in = R(-tau) (v(E = eta) - v_moon(tau)), u_out = R(+tau) (v(E = -eta) - v_moon(-tau)), with v the Kepler velocity from eq. 3 and v_moon(t) = (-sin t, cos t). The recipe returns |u_in| = |u_out| = sqrt(3 - C) to machine precision for every Table I orbit, and the demanded turn is the angle between u_in and u_out.
- COMPUTED demanded turns (160 non-degenerate orbits): from 65.09 degrees (A0(-1), v^2 = 3.399) down to 0.112 degrees. 64 of the 160 demand under 2 degrees and 23 under 0.5 degrees. Along A_i(1) the turn falls from 1.765 degrees (i = 0) to 0.557 degrees (i = 8) with v^2 falling from 0.0291 to 0.0011. The turn goes to zero as the orbit approaches the integer point (INFERRED from the trend; the integer points themselves are spurious branch crossings, p.437).
- Consequence for `#906`(b). Rule (b) says a demanded turn of 0 or 180 degrees is not an encounter. The Hitzl-Henon critical orbits that sit next to the resonance (bifurcation) points are exactly the ones with a small demanded turn at finite speed. A rule "turn near zero means no encounter" with a loose tolerance, for example 1 degree, would discard 64 of the 160 critical orbits, including the near-resonant seeds of the first-species and second-species bifurcation (Barrabes & Gomez 2003: "the generating orbits ... associated to p-q resonant orbits are bifurcation orbits of 1st species-2nd species", cited in the synthesis). INFERRED risk; the paper never mentions a gate.
- Concrete check to add to `turn_gate.py` (not to be coded under this task). At every junction compute, besides ratio, (i) the first-order periapsis r_p = mu (e_h - 1) / v^2 with e_h = 1 / sin(turn/2) (the gate's own relation) and (ii) the two validity numbers mu |ln mu| / v^3 and r_p / sqrt(mu). When the demanded turn is below the encounter tolerance, or the validity numbers are out of range, return a verdict "indeterminate: near resonance / outside first-order validity" and do not reject; report the nearest p-q ratio of the Kepler period, |a^(3/2) - p/q|, with the p/q of lowest q. Never convert "turn near zero" into "no encounter" without that report. Reason: the gate's turn comes from an integrated trajectory, but the near-resonant junction is a case where the passage is distant (r_p of order 1 at mu = 0.0122 for most of these orbits) and the turn is small for a physical reason.
- Validity at the Earth-Moon mass (COMPUTED, mu = 0.0121529529, thresholds are my choices). Of the 160 nondegenerate orbits only 90 have mu |ln mu| / v^3 below 1, and 85 below 0.1; of those 90, 75 have first-order periapsis below sqrt(mu) = 0.110 and only 1 has it below 0.1 sqrt(mu) = 0.011. So the one-moon asymptotic description cannot be applied at the Earth-Moon mass to most of Table I, and the first-order periapsis relation should not be used to rank them. (Equivalent: the second-species asymptotics need mu far below the Earth-Moon value for the near-resonant orbits, as the Gomez & Olle digest also found for the `#890` chain.)
- Does a critical generating orbit risk a wrong rejection by itself? No for the reason above: the turn at a critical orbit is a smooth function along the family and nothing happens to it at the C extremum (INFERRED; the paper's statement that "the angle may be zero in particular cases" is the only remark). The risk is at the resonance crossings and in the use of generating-orbit theory outside its mass range, not at the C extremum.

### (c) More than one close approach per period; the Earth-Moon mass ratio; two small bodies

- More than one close approach to the secondary per period: each generating orbit has exactly one angular point per period (READ p.445), so chains with several lunar encounters are not treated in this paper. Several close approaches to the primary are: "Note the five close passages by the Earth for each encounter with the Moon" for C25(1) (Fig. 14, p.448); n* = i + j for the retrograde C_ij(2) (Fig. 15); "the first critical orbit C38(1) with eight close passages by the Earth for each lunar encounter" (Fig. 16). The paper's remark (READ p.451): "stable orbits possessing close passages by both the Earth and the Moon together with rather small periods should be possible for mu greater than 0" and "It is conjectured that, for small mu > 0, the region of stability for C_ij(2) will be larger than that for C_ij(1)".
- Earth-Moon mass ratio: no mu greater than 0 orbit is computed. What the paper gives are published existence ranges of the critical periodic families that end at three generating orbits (READ captions, pp.437, 446-447): the retrograde families m1 and m2 that end at A0(-1) exist "for 0 < mu <= 1 and 0 <= mu < 0.327..." and "k jumps from -infinity to +infinity"; the double-periodic h24 and h25 that end at A0(1) exist for "0 < mu <= 0.094... and 0 < mu <= 0.102..."; the retrograde double-periodic f2 and f7 that end at B1(-1) exist for "0 < mu <= 0.156... and 0 < mu <= 0.8...". All these ranges contain mu = 0.0122. INFERRED: the critical families at those three generating orbits persist to the Earth-Moon mass. Family m3 at A0(0) (the retrograde circle surrounding both primaries) has k jumping "from -1 to -infinity" (p.445).
- Two small bodies: nothing. The paper is the restricted problem with one secondary, in the zero-mass limit.

### (d) Printed numbers as positive controls (recipe and expected values)

1. Table I as a test of `F0` and `G*` coded from eqs. 11, 27, 28. Expected: for each printed row, residuals at the level of the 10-digit rounding. COMPUTED on 2026-10-04 with 60-digit arithmetic: all 160 non-degenerate rows polish (Newton on (F0, G*) from the printed values) to shifts below 5.0e-11 in tau/pi and eta/pi and none falls onto the spurious integer points, except for the two rows with printed digit slips described below. Test inputs: only the printed tau/pi, eta/pi and the sign columns; the pass criterion must be stated in terms of the polished shift, not the raw residual, because the residual is a weak test near integer points (third and fourth order in the offset).
2. B_i(0), i = 1 to 9: tau/pi = eta/pi and tan tau = 3 tau / 4. COMPUTED: the printed values satisfy tan tau - 3 tau/4 to between 1e-9 and 6e-8 (the root of eq. 40 reproduces the printed 10 decimals; residual grows with tau because the printed value is rounded to 1e-10 in tau/pi, which multiplies by pi and by the slope). Expected values (printed): 1.406 729 614 4; 2.445 298 131 4; 3.461 162 219 9; 4.469 866 862 7; 5.475 376 062 1; 6.479 179 123 9; 7.481 963 251 2; 8.484 089 938 8; 9.485 767 638 6.
3. Table II as the test of the conversion from (tau, eta, sigma) to a, e, x0, x1, C, T (eqs. 7, 17, 43, 44): ten orbits, all reproduced to the printed digits (section 4.1). These are the expected values with a source independent of any project code.
4. Using the generating orbits to test `core/cr3bp.py` at tiny mu: place the particle at the Moon at the angular point with the incoming rotating-frame velocity from section 5(b), integrate backwards and forwards with the CR3BP equations at mu = 1e-9 or smaller with regularised close passage, and require closure of the arc over T = 2 tau to a tolerance that scales like mu^(1 - alpha). Not run; that is a scoped test, not a result. Expected: the integrated Jacobi constant equals C of Table II (x0, xdot0 with xdot0 = x0-frame from eq. 3) to rounding at any mu away from the encounter.
5. The extremum of C itself: for a sequence A_i(l), the value of C from eq. 17 on the printed (tau, eta) must alternate max/min with parity of l (even l minimum for A and C, maximum for B). COMPUTED to hold for all families in Table I; that is a check on a coded C(tau, eta), with the source of the expected pattern being the paper (READ p.445).

## 6. Errata found in the printed tables (COMPUTED)

Two printed tau/pi entries in Table I do not satisfy F0 = 0 and G* = 0, and for each a single changed digit does. Both were read twice on the page image (the printed digits are the ones below, not an OCR slip).
- A5(-2): printed tau/pi = 6.283 866 904 9, eta/pi = 4.763 717 241 5 (p.439). F0 = -0.152, G* = 0.131 at the printed point. Newton on (F0, G*) from the printed point converges to tau/pi = 6.288 866 904 884, eta/pi = 4.763 717 241 53: eta/pi agrees with the print to every printed digit and tau/pi differs by exactly 0.005 (the digit "3" in the third decimal place for "8"). The neighbours agree: tau/pi for A4(-2) is 5.276 855 614 4 and for A6(-2) is 7.298 243 699 0, so A5(-2) should lie close to their mean 6.2876, not 6.2839.
- A5(1): printed tau/pi = 6.993 912 142 9, eta/pi = 5.999 013 667 1 (p.439). F0 = -4.78e-5 at the printed point. A Newton start from the printed point lands on the spurious integer point (7, 6) (6.999 999 99, 6.000 000 01), which is what the paper's section 5 warns about. With the single digit changed, tau/pi = 6.998 912 142 9, F0 = 1.7e-13 and G* = -4.7e-20; a 2D polish from there gives 6.998 912 142 872, 5.999 013 667 133. The rows around it agree: tau/pi for A4(1) is 5.998 754 665 9 and for A6(1) is 7.999 034 838 7, and the A_i(1) sequence runs 2.997 858, 3.998 260, 4.998 546, 5.998 755, then 6.9989, 7.999 035.
Both are probable typesetting slips (not OCR errors). Anyone using Table I as test expectations should use the corrected values and cite this digest for the correction, or exclude the two rows. All other rows pass (section 5(d)).
Also noted: "The approximate location of 62 critical orbits" in the Fig. 3a caption counts only those plotted. The abstract's "179" matches the table (169 numeric rows plus 10 A_i(0) rows).

## 7. Stability-index conventions (read before comparing with Casoliva's Table 3 or the project's Barden code)

READ: this paper takes the stable range as -1 to +1 and critical orbits at k = +-1 (pp.421-422, 451). The paper itself does not define k; the definition is in the later book, READ there (Henon 1997, section 2.8, p.14; see `docs/notes/2026-10-04-digest-henon-1997-generating-families.md`): the stability index is z = (a + d)/2, half the trace of the linearised return map on a surface of section, stable for |z| less than 1, with critical orbits of the first kind at z = +1 and of the second kind at z = -1. So Henon's index is one half of the trace, equal to (lambda + 1/lambda)/2 for the full-period eigenvalue lambda. Casoliva et al. (2010, their eq. 8) use k = lambda + 1/lambda with stability |k| less than 2 (as recorded in `search/earth_moon_resonant_families.py`), so Casoliva's k equals twice Henon's and Casoliva's boundary |k| = 2 is Henon's critical value |z| = 1. The Barden half-period nu = (lambda_half + 1/lambda_half)/2 used by `search/cr3bp_periodic.py::barden_stability` is built from the half-period monodromy, so it is a different object and not interchangeable with either without an explicit conversion (INFERRED from the docstring; not checked in the code).
Orbits that sit at k = +-1 at the Earth-Moon mass in Casoliva's Table 3 are those with |k| close to 2: 1-2c (1.9996), 1-2e (1.9998), 3-2d (2.0008), 2-1b (2.0374); their closeness to the critical value is not claimed by either paper to be related to a nearby generating orbit. Not tested.

## 8. What the paper leaves open (quoted) and what to obtain

- "We cannot offer a mathematical proof that our list is exhaustive as, in fact, the set of all critical orbits is considered an unsolved theoretical problem (see Bruno, 1973)." (p.438)
- "an existence proof of the type previously given for periodic solutions of the second kind ... has not as yet been developed for the second species periodic solutions." (p.423)
- Second joint paper promised (p.451): "These results will be presented in a forthcoming second joint paper." Also planned: numerical integration of the variational equations for the stability of periodic orbits near a few critical generating orbits at small mu, and "approximate analytical formulas obtained from second-order matched asymptotic expansions (Breakwell and Perko, 1974)". The second joint paper is, READ in the book's reference list (Henon 1997, p.277), Hitzl & Henon, "The stability of second species periodic orbits in the restricted problem (mu = 0)", Acta Astronautica 4:1019-1039 (1977). Not held, not in `docs/notes/CORPUS_INDEX.md`.
- References worth obtaining (from the paper's list, p.451-452), none confirmed in `docs/notes/CORPUS_INDEX.md` by my grep for the author names (Henon, Hitzl, Bruno, Perko, Guillaume, Breakwell: no match for these papers):
  - Henon, M., "Sur les orbites interplanetaires qui rencontrent deux fois la terre", Bull. Astron. 3:377-402 (1968): the origin of F0 and of the characteristics, with the hyperbolic table and Figure 9.
  - Henon, M., "Exploration numerique du probleme restreint. I", Ann. Astrophys. 28:992-1007 (1965): the theorem that an extremal C along a periodic family is a critical orbit with k = +1, and the definition of k.
  - Henon, M. & Guyot, M., "Stability of Periodic Orbits in the Restricted Problem", in Periodic Orbits, Stability and Resonances (Giacaglia ed.), Reidel, 349-374 (1970): the critical families that end at generating orbits.
  - Perko, L. M., SIAM J. Appl. Math. 27:200-237 (1974): existence proof for second species orbits to order mu^m.
  - Breakwell, J. V. & Perko, L. M., Celest. Mech. 9:437-450 (1974): second-order matched asymptotics.
  - Bruno, A. D., Research on the Restricted Three Body Problem III, Inst. Appl. Math. Moscow Preprint 25 (1973), in Russian: the criticality equation in (a, e).
  - Guillaume, P., Celest. Mech. 8:199-206 (1973), 11:213-254 and 11:449-467 (1975): periodic symmetric solutions at larger mu.
  - Henon, M., Astron. Astrophys. 1:223-238 (1969): the first ten solutions of tan x = 3x/4.

## 9. Appendix: Table I (READ pp.438-442) with computed derived columns

Columns tau/pi and eta/pi are printed values (READ, 10 decimals; the two flagged rows are printed as shown and are wrong by one digit, see section 6). The sign column (sigma0 sigma1 sigma2) is READ from the page images. The last three columns (a, e, C) are COMPUTED by me from eqs. 7 and 17 with the corrected values for A5(-2) and A5(1); they are conveniences, NOT sourced values, and must never be used as test expectations (only Table II's printed a, e, C are sourced).

| Orbit | tau/pi | eta/pi | sigma0 sigma1 sigma2 | a (computed) | e (computed) | C (computed) |
|---|---|---|---|---|---|---|
| A0(-1) | 0.213 448 240 7 | 0.393 330 817 4 | --+ | 1.41019 | 0.88445 | -0.39913 |
| A0(0) | 0.5 | 0.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A0(1) | 1.997 355 654 8 | 0.998 206 723 0 | ++- | 1.58720 | 0.36996 | 2.97094 |
| A0(2) | 4.055 382 853 2 | 0.972 767 550 0 | +-- | 2.56117 | 0.61179 | -2.14139 |
| A0(3) | 5.999 581 166 8 | 0.999 823 070 4 | ++- | 3.30189 | 0.69714 | 2.90835 |
| A0(4) | 8.013 300 565 3 | 0.994 978 883 7 | +-- | 4.00778 | 0.75058 | -2.39618 |
| A0(5) | 9.999 821 851 3 | 0.999 938 100 8 | ++- | 4.64157 | 0.78456 | 2.88719 |
| A1(-5) | 7.999 741 236 9 | 0.999 902 196 4 | ++- | 3.99998 | 0.75000 | 2.89575 |
| A1(-4) | 6.023 790 284 4 | 0.989 980 672 8 | +-- | 3.31733 | 0.69890 | -2.30389 |
| A1(-3) | 3.999 173 790 1 | 0.999 588 921 8 | ++- | 2.51977 | 0.60314 | 2.92916 |
| A1(-2) | 2.189 676 963 3 | 0.882 973 214 1 | +-- | 1.76225 | 0.46352 | -1.78510 |
| A1(-1) | 1.419 526 801 3 | 1.493 848 705 5 | +-- | 1.00521 | 0.26816 | -0.93694 |
| A1(0) | 1.5 | 1.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A1(1) | 2.997 857 727 0 | 1.998 317 194 5 | -++ | 1.31031 | 0.23682 | 2.98743 |
| A1(2) | 5.088 519 431 4 | 1.946 970 583 7 | --+ | 1.88126 | 0.47502 | -1.88238 |
| A1(3) | 6.999 483 262 9 | 1.999 728 047 9 | -++ | 2.30520 | 0.56620 | 2.93676 |
| A1(4) | 9.022 471 916 3 | 1.989 371 058 0 | --+ | 2.73376 | 0.63456 | -2.18996 |
| A2(-5) | 8.999 660 926 0 | 1.999 839 287 6 | -++ | 2.72567 | 0.63312 | 2.92275 |
| A2(-4) | 7.039 933 954 1 | 1.979 092 466 9 | --+ | 2.32094 | 0.57037 | -2.07185 |
| A2(-3) | 4.999 091 250 1 | 1.999 445 302 0 | -++ | 1.84198 | 0.45711 | 2.95711 |
| A2(-2) | 3.235 674 264 7 | 1.832 350 236 0 | --+ | 1.43201 | 0.34898 | -1.54455 |
| A2(-1) | 2.450 053 486 8 | 2.497 677 328 7 | --+ | 1.00119 | 0.16338 | -0.97550 |
| A2(0) | 2.5 | 2.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A2(1) | 3.998 260 119 1 | 2.998 541 349 7 | ++- | 1.21139 | 0.17450 | 2.99299 |
| A2(2) | 6.113 619 304 7 | 2.924 742 373 5 | +-- | 1.62395 | 0.39521 | -1.72541 |
| A2(3) | 7.999 456 378 9 | 2.999 677 758 7 | ++- | 1.92298 | 0.47998 | 2.95311 |
| A3(-5) | 9.999 626 140 6 | 2.999 799 095 1 | ++- | 2.23143 | 0.55186 | 2.93961 |
| A3(-4) | 8.053 829 996 9 | 2.968 315 807 6 | +-- | 1.93864 | 0.48658 | -1.91699 |
| A3(-3) | 5.999 118 918 8 | 2.999 402 537 9 | ++- | 1.58738 | 0.37003 | 2.97093 |
| A3(-2) | 4.260 470 423 1 | 2.801 550 424 2 | +-- | 1.30582 | 0.28847 | -1.42250 |
| A3(-1) | 3.463 645 328 7 | 3.498 779 139 7 | +-- | 1.00045 | 0.11775 | -0.98699 |
| A3(0) | 3.5 | 3.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A3(1) | 4.998 546 140 0 | 3.998 734 945 3 | -++ | 1.16038 | 0.13822 | 2.99553 |
| A3(2) | 7.133 558 851 1 | 3.905 774 585 5 | --+ | 1.48612 | 0.34198 | -1.61824 |
| A3(3) | 8.999 457 754 0 | 3.999 652 443 6 | -++ | 1.71706 | 0.41761 | 2.96366 |
| A4(-4) | 9.066 223 004 5 | 3.957 920 505 6 | --+ | 1.73247 | 0.42651 | -1.80381 |
| A4(-3) | 6.999 173 163 4 | 3.999 400 837 5 | -++ | 1.45218 | 0.31138 | 2.97893 |
| A4(-2) | 5.276 855 614 4 | 3.780 005 281 1 | --+ | 1.23804 | 0.24954 | -1.34723 |
| A4(-1) | 4.471 392 338 6 | 4.499 247 223 4 | --+ | 1.00022 | 0.09210 | -0.99193 |
| A4(0) | 4.5 | 4.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A4(1) | 5.998 754 665 9 | 4.998 889 894 5 | ++- | 1.12923 | 0.11445 | 2.99690 |
| A4(2) | 8.149 889 843 4 | 4.889 497 392 1 | +-- | 1.39954 | 0.30359 | -1.53985 |
| A4(3) | 9.999 471 368 9 | 4.999 641 537 0 | ++- | 1.58739 | 0.37004 | 2.97093 |
| A5(-3) | 7.999 231 182 9 | 4.999 416 480 6 | ++- | 1.36797 | 0.26899 | 2.98400 |
| A5(-2) | 6.283 866 904 9 (printed value; see section 6) | 4.763 717 241 5 | +-- | 1.19545 | 0.22186 | -1.29573 |
| A5(-1) | 5.476 408 322 9 | 5.499 489 416 2 | +-- | 1.00012 | 0.07564 | -0.99451 |
| A5(0) | 5.5 | 5.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A5(1) | 6.993 912 142 9 (printed value; see section 6) | 5.999 013 667 1 | -++ | 1.10823 | 0.09766 | 2.99773 |
| A5(2) | 9.163 581 082 1 | 5.875 392 350 5 | --+ | 1.33982 | 0.27439 | -1.47980 |
| A6(-3) | 8.999 286 051 1 | 5.999 439 194 8 | -++ | 1.31036 | 0.23685 | 2.98742 |
| A6(-2) | 7.298 243 699 0 | 5.750 775 493 0 | --+ | 1.16609 | 0.20094 | -1.25809 |
| A6(-1) | 6.479 924 074 2 | 6.499 630 926 8 | --+ | 1.00007 | 0.06418 | -0.99602 |
| A6(0) | 6.5 | 6.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A6(1) | 7.999 034 838 7 | 6.999 113 822 2 | ++- | 1.09310 | 0.08517 | 2.99826 |
| A7(-3) | 9.999 335 859 5 | 6.999 464 273 7 | ++- | 1.26843 | 0.21162 | 2.98985 |
| A7(-2) | 8.305 878 320 6 | 6.740 129 891 3 | +-- | 1.14458 | 0.18445 | -1.22931 |
| A7(-1) | 7.482 526 184 8 | 7.499 720 768 9 | +-- | 1.00005 | 0.05574 | -0.99699 |
| A7(0) | 7.5 | 7.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A7(1) | 8.999 132 955 1 | 7.999 196 135 8 | -++ | 1.08168 | 0.07552 | 2.99863 |
| A8(-2) | 9.312 283 817 4 | 7.731 146 405 8 | --+ | 1.12812 | 0.17103 | -1.20652 |
| A8(-1) | 8.484 530 300 5 | 8.499 781 366 6 | --+ | 1.00003 | 0.04927 | -0.99764 |
| A8(0) | 8.5 | 8.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| A8(1) | 9.999 213 125 1 | 8.999 264 804 7 | ++- | 1.07276 | 0.06783 | 2.99889 |
| A9(-1) | 9.486 121 529 8 | 9.499 824 169 7 | +-- | 1.00002 | 0.04414 | -0.99810 |
| A9(0) | 9.5 | 9.5 | - (undetermined; sigma1 = -) | 1 | 0 | -1 |
| B1(-5) | 9.010 543 163 1 | 0.996 192 135 9 | --- | 4.33267 | 0.76925 | -2.42914 |
| B1(-4) | 6.999 676 438 0 | 0.999 871 279 2 | -+- | 3.65928 | 0.72672 | 2.90137 |
| B1(-3) | 5.034 721 259 0 | 0.984 318 065 4 | --- | 2.94799 | 0.66159 | -2.23579 |
| B1(-2) | 2.998 661 389 1 | 0.999 246 965 4 | -+- | 2.07997 | 0.51922 | 2.94591 |
| B1(-1) | 1.339 617 957 6 | 0.757 577 506 1 | --- | 1.36617 | 0.37034 | -1.43948 |
| B1(0) | 1.406 729 614 4 | 1.406 729 614 4 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B1(1) | 3.098 806 395 2 | 0.946 081 913 3 | --- | 2.16120 | 0.54510 | -2.00228 |
| B1(2) | 4.999 431 580 1 | 0.999 741 838 9 | -+- | 2.92397 | 0.65800 | 2.91727 |
| B1(3) | 7.017 378 075 8 | 0.993 099 303 9 | --- | 3.66996 | 0.72769 | -2.35552 |
| B1(4) | 8.999 787 519 3 | 0.999 923 194 8 | -+- | 4.32673 | 0.76888 | 2.89111 |
| B2(-4) | 7.999 586 921 5 | 1.999 794 476 8 | +++ | 2.51982 | 0.60315 | 2.92916 |
| B2(-3) | 6.057 509 194 7 | 1.967 949 769 7 | +-+ | 2.10401 | 0.52739 | -1.98951 |
| B2(-2) | 3.998 678 292 2 | 1.999 103 737 3 | +++ | 1.58735 | 0.37002 | 2.97094 |
| B2(-1) | 2.360 251 762 9 | 1.715 045 145 0 | +-+ | 1.20569 | 0.27281 | -1.28338 |
| B2(0) | 2.445 298 131 4 | 2.445 298 131 4 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B2(1) | 4.144 354 485 2 | 1.906 138 602 4 | +-+ | 1.65618 | 0.41407 | -1.73904 |
| B2(2) | 5.999 330 782 8 | 1.999 623 542 5 | +++ | 2.08005 | 0.51924 | 2.94591 |
| B2(3) | 8.029 315 229 9 | 1.985 460 833 0 | +-+ | 2.53084 | 0.60551 | -2.13702 |
| B2(4) | 9.999 715 800 3 | 1.999 870 925 1 | +++ | 2.92401 | 0.65800 | 2.91727 |
| B3(-4) | 8.999 553 866 1 | 2.999 749 035 7 | -+- | 2.08007 | 0.51925 | 2.94591 |
| B3(-3) | 7.076 353 232 3 | 2.952 477 145 5 | --- | 1.78249 | 0.44393 | -1.83166 |
| B3(-2) | 4.998 803 424 7 | 2.999 110 934 3 | -+- | 1.40569 | 0.28861 | 2.98173 |
| B3(-1) | 3.370 803 627 3 | 2.692 261 854 3 | --- | 1.14514 | 0.22316 | -1.21299 |
| B3(0) | 3.461 162 219 9 | 3.461 162 219 9 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B3(1) | 5.173 684 910 5 | 2.876 920 476 8 | --- | 1.46489 | 0.34265 | -1.59147 |
| B3(2) | 6.999 319 810 8 | 2.999 571 383 7 | -+- | 1.75919 | 0.43156 | 2.96140 |
| B3(3) | 9.039 782 632 6 | 2.977 718 439 6 | --- | 2.09115 | 0.52308 | -1.98675 |
| B4(-4) | 9.999 545 648 5 | 3.999 722 668 2 | +++ | 1.84201 | 0.45711 | 2.95711 |
| B4(-3) | 8.092 491 004 4 | 3.938 203 084 4 | +-+ | 1.60985 | 0.38608 | -1.71967 |
| B4(-2) | 5.998 929 043 3 | 3.999 158 761 7 | +++ | 1.31035 | 0.23685 | 2.98742 |
| B4(-1) | 4.377 844 835 4 | 3.677 198 550 0 | +-+ | 1.11286 | 0.19193 | -1.17203 |
| B4(0) | 4.469 866 862 7 | 4.469 866 862 7 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B4(1) | 6.194 680 578 0 | 3.854 510 353 7 | +-+ | 1.36225 | 0.29634 | -1.49538 |
| B4(2) | 7.999 339 204 2 | 3.999 551 915 6 | +++ | 1.58739 | 0.37003 | 2.97093 |
| B5(-3) | 9.106 549 736 5 | 4.925 149 069 2 | --- | 1.50129 | 0.34335 | -1.63547 |
| B5(-2) | 6.999 038 142 9 | 4.999 215 407 9 | -+- | 1.25145 | 0.20093 | 2.99081 |
| B5(-1) | 5.383 116 669 4 | 4.666 163 575 9 | --- | 1.09266 | 0.17007 | -1.14495 |
| B5(0) | 5.475 376 062 1 | 5.475 376 062 1 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B5(1) | 7.210 728 717 5 | 4.836 620 756 7 | --- | 1.29769 | 0.26333 | -1.42731 |
| B5(2) | 8.999 369 091 1 | 4.999 549 286 1 | -+- | 1.47972 | 0.32420 | 2.97728 |
| B6(-2) | 7.999 130 142 4 | 5.999 270 754 1 | +++ | 1.21141 | 0.17451 | 2.99299 |
| B6(-1) | 6.387 322 890 1 | 5.657 567 992 5 | +-+ | 1.07878 | 0.15372 | -1.12562 |
| B6(0) | 6.479 179 123 9 | 6.479 179 123 9 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B6(1) | 8.223 558 134 9 | 5.821 891 637 4 | +-+ | 1.25314 | 0.23836 | -1.37636 |
| B6(2) | 9.999 401 747 8 | 5.999 555 498 2 | +++ | 1.40571 | 0.28862 | 2.98173 |
| B7(-2) | 8.999 207 560 9 | 6.999 321 680 9 | -+- | 1.18239 | 0.15426 | 2.99447 |
| B7(-1) | 7.390 816 533 8 | 6.650 592 474 9 | --- | 1.06862 | 0.14092 | -1.11106 |
| B7(0) | 7.481 963 251 2 | 7.481 963 251 2 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B7(1) | 9.234 153 325 9 | 6.809 470 371 1 | --- | 1.22047 | 0.21866 | -1.33667 |
| B8(-2) | 9.999 273 114 2 | 7.999 367 515 8 | +++ | 1.16039 | 0.13822 | 2.99553 |
| B8(-1) | 8.393 799 819 9 | 7.644 763 095 5 | +-+ | 1.06085 | 0.13057 | -1.09967 |
| B8(0) | 8.484 089 938 8 | 8.484 089 938 8 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| B9(-1) | 9.396 399 457 1 | 8.639 782 834 1 | --- | 1.05470 | 0.12199 | -1.09050 |
| B9(0) | 9.485 767 638 6 | 9.485 767 638 6 | (sigma0, sigma2 undetermined; sigma1 = +) | 1 | 0 | 3 |
| C1,2(1) | 0.987 254 706 7 | 1.975 221 027 7 | ++- | 0.63289 | 0.58182 | 2.87412 |
| C1,2(2) | 1.384 687 370 5 | 1.511 649 180 4 | +-- | 0.98836 | 0.32198 | -0.87066 |
| C2,3(1) | 1.994 514 662 3 | 2.992 442 288 9 | -++ | 0.76343 | 0.30997 | 2.97130 |
| C2,3(2) | 2.438 851 706 4 | 2.503 314 987 2 | --+ | 0.99812 | 0.18087 | -0.96328 |
| C2,4(1) | 1.993 619 174 0 | 3.987 512 498 0 | -+- | 0.63070 | 0.58598 | 2.87260 |
| C2,4(2) | 2.121 582 395 1 | 3.772 262 779 5 | --- | 0.69628 | 0.57791 | 0.07423 |
| C2,5(1) | 1.993 657 204 2 | 4.978 633 163 7 | -++ | 0.54459 | 0.83813 | 2.64131 |
| C2,5(2) | 2.036 871 879 1 | 4.882 804 131 7 | --+ | 0.56560 | 0.82320 | 0.91406 |
| C3,4(1) | 2.996 602 057 0 | 3.995 788 939 1 | ++- | 0.82556 | 0.21132 | 2.98747 |
| C3,4(2) | 3.458 087 698 7 | 3.501 566 855 3 | +-- | 0.99938 | 0.12645 | -0.98271 |
| C3,5(1) | 2.996 109 918 3 | 4.994 019 252 3 | +++ | 0.71156 | 0.40543 | 2.94756 |
| C3,5(2) | 3.172 580 121 5 | 4.737 615 897 7 | +-+ | 0.77629 | 0.42436 | -0.30744 |
| C3,6(1) | 2.995 745 049 0 | 5.991 662 489 7 | ++- | 0.63029 | 0.58677 | 2.87231 |
| C3,6(2) | 3.083 113 588 4 | 5.839 320 859 6 | +-- | 0.66010 | 0.58829 | 0.20092 |
| C3,7(1) | 2.995 627 080 4 | 6.988 215 916 0 | +++ | 0.56901 | 0.75797 | 2.74153 |
| C3,7(2) | 3.037 753 377 8 | 6.900 300 030 1 | +-+ | 0.58287 | 0.75226 | 0.70963 |
| C3,8(1) | 2.996 294 347 2 | 7.981 793 913 9 | ++- | 0.52111 | 0.92048 | 2.48318 |
| C3,8(2) | 3.012 343 216 4 | 7.941 476 240 2 | +-- | 0.52636 | 0.91526 | 1.31528 |
| C4,5(1) | 3.997 555 025 7 | 4.997 125 756 5 | -++ | 0.86181 | 0.16036 | 2.99300 |
| C4,5(2) | 4.468 067 666 9 | 4.500 912 723 3 | --+ | 0.99972 | 0.09731 | -0.98995 |
| C4,6(1) | 3.997 257 630 2 | 5.996 220 236 7 | -+- | 0.76321 | 0.31027 | 2.97126 |
| C4,6(2) | 4.205 416 142 3 | 5.716 434 051 6 | --- | 0.82308 | 0.34187 | -0.49019 |
| C4,7(1) | 3.997 003 088 4 | 6.995 121 807 1 | -++ | 0.68873 | 0.45200 | 2.93252 |
| C4,7(2) | 4.114 064 760 8 | 6.814 490 189 5 | --+ | 0.72008 | 0.46560 | -0.11322 |
| C4,8(1) | 3.996 808 503 4 | 7.993 743 560 2 | -+- | 0.63015 | 0.58704 | 2.87221 |
| C4,8(2) | 4.063 831 346 8 | 7.875 592 326 6 | --- | 0.64725 | 0.58944 | 0.24518 |
| C4,9(1) | 3.996 714 308 1 | 8.991 915 888 7 | -++ | 0.58267 | 0.71648 | 2.78126 |
| C4,9(2) | 4.034 970 332 6 | 8.914 791 852 2 | --+ | 0.59217 | 0.71414 | 0.61137 |
| C4,10(1) | 3.996 833 019 9 | 9.989 223 350 5 | -+- | 0.54332 | 0.84103 | 2.63807 |
| C4,10(2) | 4.017 352 127 4 | 9.941 781 209 7 | --- | 0.54836 | 0.83759 | 1.01453 |
| C5,6(1) | 4.998 095 227 3 | 5.997 830 896 1 | ++- | 0.88556 | 0.12923 | 2.99553 |
| C5,6(2) | 5.474 195 111 7 | 5.500 597 419 2 | +-- | 0.99985 | 0.07911 | -0.99343 |
| C5,7(1) | 4.997 897 630 4 | 6.997 281 745 0 | +++ | 0.79910 | 0.25142 | 2.98183 |
| C5,7(2) | 5.229 148 884 5 | 6.701 286 663 5 | +-+ | 0.85391 | 0.28946 | -0.59794 |
| C5,8(1) | 4.997 720 027 8 | 7.996 646 034 1 | ++- | 0.73106 | 0.36789 | 2.95799 |
| C5,8(2) | 5.137 701 773 5 | 7.796 287 568 0 | +-- | 0.76210 | 0.38917 | -0.29617 |
| C5,9(1) | 4.997 566 794 4 | 8.995 897 315 6 | +++ | 0.67588 | 0.47958 | 2.92236 |
| C5,9(2) | 5.084 114 006 5 | 8.858 095 034 5 | +-+ | 0.69414 | 0.48836 | -0.01347 |
| C5,10(1) | 4.997 446 697 3 | 9.994 993 620 3 | ++- | 0.63008 | 0.58717 | 2.87216 |
| C5,10(2) | 5.051 855 673 8 | 9.898 634 339 6 | +-- | 0.64117 | 0.58928 | 0.26578 |
| C6,7(1) | 5.998 441 660 5 | 6.998 262 809 0 | -++ | 0.90235 | 0.10822 | 2.99690 |
| C6,7(2) | 6.478 344 516 4 | 6.500 421 416 2 | --+ | 0.99991 | 0.06666 | -0.99537 |
| C6,8(1) | 5.998 301 193 4 | 7.997 894 505 9 | -+- | 0.82550 | 0.21139 | 2.98746 |
| C6,8(2) | 6.247 420 840 5 | 7.689 619 497 7 | --- | 0.87576 | 0.25285 | -0.66895 |
| C6,9(1) | 5.998 171 790 0 | 8.997 480 044 9 | -++ | 0.76317 | 0.31033 | 2.97125 |
| C6,9(2) | 6.156 784 263 7 | 8.781 801 805 4 | --+ | 0.79328 | 0.33664 | -0.41676 |
| C6,10(1) | 5.998 054 933 2 | 9.997 008 772 1 | -+- | 0.71142 | 0.40565 | 2.94752 |
| C6,10(2) | 6.101 006 501 4 | 9.843 978 982 9 | --- | 0.73003 | 0.41916 | -0.18167 |
| C7,8(1) | 6.998 682 236 0 | 7.998 553 272 6 | ++- | 0.91483 | 0.09310 | 2.99773 |
| C7,8(2) | 7.481 342 125 4 | 7.500 313 197 4 | +-- | 0.99994 | 0.05760 | -0.99657 |
| C7,9(1) | 6.998 577 359 4 | 8.998 289 203 1 | +++ | 0.84575 | 0.18238 | 2.99083 |
| C7,9(2) | 7.262 081 881 2 | 8.680 217 597 4 | +-+ | 0.89204 | 0.22563 | -0.71921 |
| C7,10(1) | 6.998 479 287 0 | 9.997 997 641 2 | ++- | 0.78839 | 0.26841 | 2.97907 |
| C7,10(2) | 7.172 701 451 7 | 9.769 785 485 7 | +-- | 0.81733 | 0.29812 | -0.50242 |
| C8,9(1) | 7.998 858 844 1 | 8.998 761 492 2 | -++ | 0.92449 | 0.08168 | 2.99826 |
| C8,9(2) | 8.483 609 680 0 | 8.500 241 917 7 | --+ | 0.99996 | 0.05071 | -0.99735 |
| C8,10(1) | 7.998 777 593 0 | 9.998 562 929 9 | -+- | 0.86178 | 0.16039 | 2.99299 |
| C8,10(2) | 8.274 198 672 5 | 9.672 402 344 2 | --- | 0.90463 | 0.20451 | -0.75660 |
| C9,10(1) | 8.998 993 912 6 | 9.998 917 837 9 | ++- | 0.93217 | 0.07276 | 2.99863 |
| C9,10(2) | 9.485 385 209 0 | 9.500 192 484 1 | +-- | 0.99997 | 0.04529 | -0.99789 |
