# Digest: Brown, Peterson, Henry & Scheeres 2025, behaviour of symmetric periodic orbits with high-order resonance in the Sun-Earth-Moon system

Date: 2026-10-04. Task: `#884` literature follow-up (see `docs/notes/2026-10-04-884-literature-check.md`, `docs/notes/2026-10-04-884-adversarial-review.md`).

Paper: Gavin M. Brown, Luke T. Peterson, Damennick B. Henry and Daniel J. Scheeres, "Behavior of Symmetric Periodic
Orbits with High-Order Resonance in the Sun-Earth-Moon System", The Journal of the Astronautical Sciences 72:43
(2025), DOI 10.1007/s40295-025-00517-3, 34 pages, open access (CC BY 4.0), accepted 8 July 2025, published online
5 August 2025. A preliminary version was AAS 24-204 (AAS/AIAA Astrodynamics Specialist Conference, Broomfield, 2024).
Filed in the private paper corpus as
brown-peterson-henry-scheeres-2025-behavior-symmetric-periodic-orbits-high-order-resonance-sun-earth-moon-jas-72-43-doi-10.1007-s40295-025-00517-3.pdf.

This is the sequel to the group's 2024 paper (digest: `docs/notes/2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md`;
SIAM J. Appl. Dyn. Syst. 24(1):346-375, 2025). Paper [25] below is that paper.

Marking: READ = stated or printed in the paper (journal page given; the page numbers are the paper's "Page n of 34",
the article number 43 is not repeated). INFERRED = my deduction, not printed. COMPUTED = a short arithmetic check I
ran for this note. Tables were transcribed from the PDF text layer and checked against the page images; they agreed.
Equation numbers are the paper's.

## 0. Summary

- READ (abstract, p1): "The Earth-Moon L2 9:2 Near-Rectilinear Halo Orbit (NRHO) is a periodic orbit in the Circular Restricted 3-Body Problem that is representative of the lunar Gateway's planned orbit. ... In this work, we use Melnikov theory and a continuation algorithm to transition this orbit into the Sun-Earth-Moon (SEM) system using the Hill Restricted 4-Body Problem (HR4BP). The set of periodic orbits corresponding to the 9:2 NRHO in the SEM HR4BP numerically foliate a 2D torus. We hypothesize that this behavior constitutes a limiting case where the destruction of the torus by the 9:2 resonance in the SEM system is not numerically detectable. As this behavior is unexpected, we extend our analysis to other resonant periodic orbits and dynamical models. We find that the SEM HR4BP dynamical equivalents to the L2 7:2 NRHO and 5:2 halo orbit also exhibit similar behavior, and the 9:2 NRHO in the SEM Bicircular Restricted 4-Body Problem exhibits this behavior as well."
- READ: the orbits treated are L2 near-rectilinear halo orbits (9:2, 7:2, 5:2 halo; 3:1 and 4:1 only by reference) and one L4 long-period planar orbit (4:13.5, and 3:11 in the bicircular appendix). No cycler, no Earth-Moon transfer orbit, no orbit defined by repeated lunar flybys, and no low-energy (C near 3) Earth-Moon orbit is treated. The words cycler, flyby, perilune, Jacobi constant and transfer (other than in a reference title) do not occur in the text (COMPUTED: text search of the whole PDF).
- READ: no initial conditions, periods of the continued orbits, Jacobi constants or perilune radii are printed. The only tables are two stability tables for one orbit (Tables 1, 2) and two coefficient tables of the Hill variation orbit (Tables 3, 4). "Data Availability: Not applicable. Code Availability: Not applicable." (p32).
- READ: the useful new content for the project is a rule for the order at which the Melnikov function first stops vanishing, i_n = 2a (Hill problem) and j_n = a - 1 (bicircular), where the resonance is aT* = bTg. Large a means the Sun's first-order effect is weak and the resonance is numerically hard to detect (section 5).
- INFERRED: the paper's bicircular appendix uses exactly the project's frame and Sun sense (Earth at (-mu, 0), Moon at (1 - mu, 0), Sun angle theta_S = -omega alpha + theta_S0 with theta_S0 = pi at the start); see section 1.2. That is an independent published check of the `#891` correction.

## 1. The two models as printed

### 1.1 Hill restricted four-body problem (HR4BP), Sect. 2.1, p3-4 and Eq. 1

READ. Quote (p3): "The equations of motion for the HR4BP are time-periodic with a period of Tg = pi, and were originally presented by Scheeres [23]." P0 is the Sun, P1 the Earth, P2 the Moon, P3 the massless particle. The paper does not repeat the scaling; it refers the reader to Appendix 1 of Brown et al. [25]. From the 2024 digest (READ there): 1 DU is the mean Earth-Moon distance, 1 TU = (1 + m)/n, so "2 pi time units is equal to a synodic month" and the equations are pi-periodic. So in this paper's HR4BP one synodic month is 2 pi and the forcing period Tg = pi is half a synodic month (the quadrupole term has period pi).

Equations (p4, Eq. 1; X is the position and velocity of P3 in the rotating frame B, whose origin is the barycentre of P1 and P2 and which rotates in the direction Omega = k-hat):

- 1a: X' = [ r' ; -2(1 + m) Omega x r' + grad V ]; 1b: [Phi]' = [A][Phi]; 1c: [Psi]' = [A][Psi] + [C], with [Psi] = [dX/dm, dX/dmu].
- 1d: V = (1/2)(1 + 2m + (3/2) m^2)(x^2 + y^2) - (1/2) m^2 z^2 + (3/4) m^2 ((x^2 - y^2) cos(2 tau) - 2 x y sin(2 tau)) + (m^2 / a0^3) [ (1 - mu)/R_(1-mu) + mu/R_mu ].
- True time tau = tau0 + alpha; "tau0 represents the initial phasing of P0, P1, and P2"; alpha is the elapsed time. tau0 = 0 for the entire work (p5).
- Parameters: "The first parameter m represents the effect of P0. Note that the CR3BP equations of motion are obtained as m -> 0 as m^2/a0^3 -> 1 in this case [23]. When m > 0 the positions of P1 and P2 vary with time, unlike in the BCP and CR3BP where P1 and P2 remain fixed. The HR4BP's representation of the SEM system is obtained when m = m_SEM = 0.0808. The second parameter mu is the ratio of the mass of P2 to the combined mass of P1 and P2. Note mu is equivalent to the CR3BP mass-ratio." mu_EM is printed as 0.0122 (p7: "mu_EM ~ 0.0122"). a0 is not given a value here.
- Hill variation orbit: the primaries' motion is the Fourier solution of Hill's variational problem with coefficients d_p (Table 3 here) and c_(n,p) (Table 4 here), computed to maximum order P = 9, N = ceil((P - 1)/2) (Appendix 1, p25), so g_i is accurate only for i <= P (p15).

Sun position and sense (INFERRED, from Eq. 1d and Fig. 1, p4). Equation 1d's solar term (3/4) m^2 ((x^2 - y^2) cos 2 tau - 2 x y sin 2 tau) equals (3/2) m^2 [ (r . s)^2 - r^2/2 ] with s = (cos tau, -sin tau) (COMPUTED by expanding the square), so the Sun-line axis sits at angle -tau in B: the Sun moves retrograde (clockwise) in the rotating frame, as in every other Sun-Earth-Moon model. The quadrupole cannot tell the Sun from the anti-Sun, so the sign of the axis at tau = 0 comes only from Fig. 1, which draws the Sun on the Earth's side of the barycentre (sequence Sun-Earth-Moon at tau = 0, Earth to the left of the Moon). The paper states the same for the bicircular model (next subsection).

### 1.2 Bicircular problem (BCP), Sect. 2.3 and Appendix 2, p7 and p29-31

READ. Frame (Appendix 2, p29): the standard CR3BP frame, "1 DU = 1 LU ~ 3.84 x 10^5 km and the Earth and Moon are fixed in this frame". Eq. 13c,d: Earth at r_(1-mu) = -mu x-hat, Moon at r_mu = (1 - mu) x-hat. Same as the project's convention.

Equations (Eq. 11): X' = [r'; -2 Omega x r' + grad V]; V (Eq. 11d) = (1/2)(x^2 + y^2) + (1 - mu)/R_(1-mu) + mu/R_mu + eps_B mu_S ( 1/R_S - (1/a_S^3) r . r_S ). The last factor includes the indirect term of the Sun (the term in r . r_S). eps_B is the BCP analogue of m: "We obtain motion in the EM CR3BP and SEM BCP when eps_B = 0 and eps_B = eps_B,SEM = 1, respectively" (p7).

Constants (Eq. 12, p30): mu_S = m_S / (m_E + m_M) = 3.2906 x 10^5; a_S = 1 AU / 1 LU = 389.1724; omega = 1 - sqrt((1 + mu_S)/a_S^3) = 0.9252. Quote (p30): "Note 1/omega = 1 + m_SEM = 1.0808 and the perturbation is periodic with Tg = 2 pi/omega = 2 pi (1.0808), where time is scaled similarly to the EM CR3BP such that 2 pi TU is one sidereal month, which is approximately 27.3 days." COMPUTED: with the printed mu_S and a_S, 1 - sqrt(329061/389.1724^3) = 0.92528, consistent with 0.9252.

Sun position (Eq. 13, p30): theta_S = -omega alpha + theta_S0; r_S = a_S (cos theta_S x-hat + sin theta_S y-hat). Quote (p30): "We use theta_S,0 = pi in our analysis of the BCP to be consistent with our selection of tau0 = 0 in our analysis of the HR4BP such that the Earth lies between the Sun and Moon on the x-axis at the initial time." So the Sun is on the -x axis at the start and its angle decreases with time (retrograde): in the project's notation, angle theta_sun0 - omega_sun t with theta_sun0 = pi. This agrees with the corrected `core/bcr4bp.py` (`#891`).

Forcing periods, as the lead asked: HR4BP Tg = pi (half a synodic month, in a time unit with 2 pi = one synodic month); BCP Tg = 2 pi (1.0808) (one synodic month, in a time unit with 2 pi = one sidereal month, 27.3 d). In days (INFERRED from the printed time scalings): the BCP Tg is one synodic month, 29.53 d; the HR4BP Tg is half a synodic month, 14.77 d.

Differences from the project's constants (READ here, project values from `core/qbcp.py` and `core/bcr4bp.py` docstrings): a_S = 389.1724 here (1 AU over a 384,400 km unit) against 388.81114 in Andreu's tables; mu_S = 3.2906 x 10^5 here (five digits) against 328900.54; omega 0.9252 here against 0.925195985520347. These are not the same Sun-distance constants; Brown et al. use the plain AU, Andreu the barycentric distance in the quasi-bicircular solution.

## 2. The Melnikov-type function and the resonance condition

### 2.1 Resonance condition (Sect. 2.2, p5)

READ. "We start with an EM CR3BP periodic orbit that has a period T* and satisfies the appropriate resonance condition aT* = bTg where a and b are relatively prime integers and T* is the period of the orbit in the CR3BP." The perturbed orbit has period T = bTg. "Note we use the traditional convention when describing resonant periodic orbits in the CR3BP (e.g., the CR3BP periodic orbit used to initialize the 9:2 NRHO family in the HR4BP has a period of T* = (2/9)(2 pi) = 4 pi/9, not (2/9)(pi) = 2 pi/9)." So "9:2" counts NRHO revolutions per synodic month pair (9 revolutions in 2 synodic months of 2 pi each) and not in units of Tg. In the HR4BP the 9:2 NRHO has a = 9, b = 4 (9 x 4 pi/9 = 4 pi = 4 Tg); in the BCP a = 9, b = 2 (p20: "b = 2 for the orbits in the BCP and b = 4 for the orbits in the HR4BP").

"High-order resonance" (the title): READ (p15) that large a (not large b) is what matters: "i_n = 2a (where aT* = bTg) in almost all cases." So high-order means large a, the number of orbit revolutions before the orbit's phase and the Sun's phase repeat. The paper never defines "high-order" in one sentence; this is my reading of the results in section 5 below (INFERRED).

### 2.2 Equations of motion expanded in m (Sect. 2.2, p5-6; Eqs. 2-3 and Appendix 1)

READ. f = f_C(X) + eps g(X, tau0 + alpha; Tg, eps) with eps = m (Eq. 2a); eps g = sum_(i>=1) g_i m^i (2b); g_i = g_a,i + g_b,i + g_c,i (3a):

- g_a,i = -sigma_i ( (1 - mu)/R^3_(1-mu,C) R_(1-mu,C) + mu/R^3_(mu,C) R_(mu,C) ), all orders i >= 1 (3b);
- g_b,i = 0 for i = 1, and -sum_(j=0)^(i-2) sigma_(i-j-2) ( (1 - mu)/R^3_(1-mu,C) psi_(1-mu,j+2) + mu/R^3_(mu,C) psi_(mu,j+2) ) for i >= 2 (3c);
- g_c,i = -2 Omega x r' - Omega x Omega x r for i = 1; (3/2) [ x + x cos 2tau - y sin 2tau ; y - y cos 2tau - x sin 2tau ; -(2/3) z ] for i = 2; and 0 for i >= 3 (3d).

So the explicit Sun-driven time dependence (cos 2 tau, sin 2 tau) is in g_c,2 (order m^2) and in the Hill-orbit harmonics that enter g_b,i for i >= 2. "g_a,i appears for all orders (i >= 1), g_b,i for i >= 2 and g_c,i for i <= 2" (p5). Appendix 1 (p25-28) gives the recursions (Sets A to D2, Eqs. 6-10) to arbitrary order K <= P.

### 2.3 The Melnikov function (Sect. 2.2, p6-7; Eq. 4)

READ. "To identify initial continuation points on a CR3BP resonant periodic orbit, we use a Melnikov function, based on the form presented by Cenedese and Haller [29], that represents the leading order term in the expansion of an energy function. The energy function is the energy balance over one period of a periodic orbit, which can be computed by evaluating the work done by the perturbing force on the orbit [29, 30]." Definitions (4a-4c):

eps M(s0, tau0) = integral from 0 to aT* of eps g(*X(s0 + alpha), tau0 + alpha) . *r'(s0 + alpha) d alpha = sum_i m^i M_i(s0, tau0), with M_i = integral from 0 to aT* of g_i(*X(s0 + alpha), tau0 + alpha) . *r'(s0 + alpha) d alpha.

The star marks the unperturbed CR3BP orbit; s0 is the point on that orbit used as initial state when tau = tau0 (p6: "Note that the variable s was used by Brown et al. [25] instead of s0"). Assumption quoted (p6): "we assume that each state X on a perturbed periodic orbit in the HR4BP (with small eps = m) is sufficiently close in phase space to a corresponding state *X on the unperturbed orbit in the CR3BP such that X = *X + O(eps). However, the only requirement for us to reliably use the zeros of M is that the energy function evaluated on the perturbed orbit is accurately represented by its leading order term, and hence by M." Properties from [25] (p6): "M_1(s0, tau0) = 0 for any s0 and tau0. When M_i is identically zero for a particular orbit we compute the next order M_(i+1)."

Continuation then follows from the zeros: "We then initialize a pseudo-arclength continuation scheme starting from points on the CR3BP periodic orbit where M_i = 0 (provided M_i is not identically zero) ... We use tau0 = 0 when evaluating M_i and throughout the continuation. Each periodic orbit in a family is corrected such that the magnitude of the constraint vector F = X_f - X_0 is less than 1e-10, where X_f = phi_T(X_0) is the time-T flow of X_0 and T = bTg" (p7). Bifurcation indicator: sigma_alpha, the smallest singular value of the 6 x 7 corrections Jacobian [DB] = [ [Phi(T,0)] - [I], Psi_m(T,0) ] other than the one that is exactly zero (family direction): "a very small value of sigma_alpha indicates the potential existence of a HR4BP periodic orbit that is both close to the current periodic orbit and does not lie on the family" (p7).

### 2.4 Bicircular Melnikov function (Appendix 2, p30-32; Eqs. 14-17)

READ. In the BCP "eps = eps_B and g = g_1 eps_B + 0 (i.e., g_i = 0 for i >= 2). As a result, the Melnikov function will only have a first order term": g = g_1 = -mu_S ( R_S/R_S^3 + r_S/a_S^3 ) (14a); M_1(s0, theta_S,0) = integral from 0 to aT* of g_1(*X(s0 + alpha), theta_S) . *r'(s0 + alpha) d alpha (14b). Because r_S = a_S is much larger than r near the Moon, 1/R_S^3 is Taylor-expanded in r/a_S (Eqs. 15a-f): g_1 = -(mu_S/a_S^2)( (1/a_S) r + (1/a_S) beta r - beta r-hat_S ), and g_1 = sum_(j>=1) (mu_S/a_S^(2+j)) g_(1,j) (16a), with

- g_(1,1) = -( r - 3 (r^T r-hat_S) r-hat_S ) (16b), the quadrupole (tidal) term;
- g_(1,2) = -(3/2) ( 2 (r^T r-hat_S) r-hat + 2 (r^T r) r-hat_S - 5 (r^T r-hat_S)^2 r-hat_S ) (16c), the octupole term, transcribed as printed (the first vector inside the bracket carries a caret in the printed text, presumably a typesetting slip for r; I did not rederive it).

M_1 = sum_j (mu_S/a_S^(2+j)) M_(1,j) (17a), M_(1,j) the integral of g_(1,j) . r' (17b). The paper notes the contributions of g_(1,j) shrink as 1/a_S^j. Jorba-Cusco et al. [13] used the same expansion for loss of symmetry in the BCP "but that work did not consider its significance to the evaluation of the Melnikov function".

### 2.5 The two rules (statements quoted)

HR4BP (Sect. 3.4, p15): "First, there appears to be a relationship between the number of dynamical equivalents that persist at non-negligible but small values of m and the lowest order of M_(i_n) that is not identically zero. In this way, Delta s0 for M_(i_n) represents the size of the holes in the invariant curve. For most orbits evaluated by Brown et al. [25], i_n = 2 and the zeros of M_(i_n) are separated by Delta s0 = pi/2 ... For the L1 1:0.5 (or 2:1) vertical orbit, i_n = 4 and Delta s0 = pi/4. For the L2 3:1 NRHO, i_n = 6 and Delta s0 = pi/6. For the L4 4:13.5 long-period planar orbit, i_n = 8 and Delta s0 = pi/8. ... While we have not derived an analytical proof guaranteeing the simple zeros of M_(i_n) are separated by Delta s0 = pi/i_n except for i_n = 2 (see Proposition 2 by Brown et al. [25]), such a result would not be surprising. Second, i_n = 2a (where aT* = bTg) in almost all cases. One notable exception to this is the L1 1:0.5 vertical orbit with i_n = 4, but this orbit's double-symmetry is likely the reason for this particular exception [25]. It is noteworthy that i_n has been even for all orbits evaluated by Brown et al. [25] and in this work."

Consistency checks (COMPUTED from the periods printed or implied): L2 3:1 NRHO, T* = 2 pi/3, a = 3, i_n = 6 ok; L4 4:13.5, T* = 27 pi/4, a = 4, i_n = 8 ok; L2 5:2 halo, T* = 4 pi/5 (INFERRED from the 5:2 label), a = 5, i_n = 10 as found with P = 15 (p15: "found i_n = 10 and Delta s = pi/10").

Why harmonics matter (p15): "This spacing Delta s0 may be related to the underlying trigonometric functions present in Eqs. 3d, 7, 9, and 10, which appear at certain orders of the Melnikov function when the relevant Hill variation orbit coefficients c_(n,p) become non-zero (see Table 4 in Appendix 1)." Which harmonics enter at which order (INFERRED from Table 4, section 3 below): the harmonic cos(2 n tau), sin(2 n tau) of the primaries' Hill motion first has a non-zero coefficient at order p = 2|n| (|n| = 1 at p = 2, |n| = 2 at p = 4, |n| = 3 at p = 6, |n| = 4 at p = 8), while the explicit Sun term g_c,2 carries only the second harmonic (n = 1) at i = 2.

BCP (Appendix 2, p32): "Based on evaluations of M_(1,j) for the EM CR3BP L2 5:2 halo orbit, L2 3:1 NRHO, L2 4:1 NRHO, L2 7:2 NRHO, and L4 3:11 long-period planar orbit, we have found two results that we expect to apply for most resonant periodic orbits satisfying aT* = bTg: (1) the lowest degree of the 1st-order BCP Melnikov function that is not identically zero is M_(1,j_n) with j_n = a - 1, and (2) the gaps between the simple zeros of M_(1,j_n) are Delta s0 = (pi/omega)/(j_n + 1)."

## 3. Tables (every table of the paper)

### Table 1 (p9): stability of HR4BP PO 0 (the SEM 9:2 NRHO dynamical equivalent, m = 0.0808)

Caption: "Eigenvalues (lambda) of [Phi(T, 0)] and singular values (sigma) of [DB]_red." Columns are ordered independently (the paper says so, p9).

| lambda | sigma |
|---|---|
| -1.1608 x 10^3 | 1.6696 x 10^3 |
| -8.6150 x 10^-4 | 69.1574 |
| 1 - (9.8820 x 10^-10 - i 4.4191 x 10^-5) | 1.2382 |
| 1 - (9.8820 x 10^-10 + i 4.4191 x 10^-5) | 1.0366 |
| 0.4689 + i 0.8832 | 0.6081 |
| 0.4689 - i 0.8832 | 2.6770 x 10^-11 |

Text (p9): "the monodromy matrix of each periodic orbit in the family has two eigenvalues (lambda_3 and lambda_4) that are close to unity. Furthermore, the maximum value of sigma_alpha for any member in this periodic orbit family is 1.1 x 10^-9." "the smallest singular value of [DB]_red for PO 0 is approximately 2.7 x 10^-11. It is important to note that the nullspace of this reduced corrections Jacobian is technically empty. However, such a small singular value indicates that another periodic orbit in the SEM system may exist close to PO 0." COMPUTED: the real pair multiplies to 1.00004 (symplectic pair) and its stability index is -1160.8; the complex pair has modulus 0.99997 and argument 62.0 degrees, so one oscillatory (stable) pair, one strongly unstable pair, and one pair at unity that is the torus-like direction.

### Table 2 (p12): stability across the whole HR4BP 9:2 family

Caption: "Stability of HR4BP periodic orbits in the SEM 9:2 NRHO family." Columns: eigenvalues of PO 0; for each, the largest magnitude of the difference between PO 0 and every member of the family.

| PO 0 eigenvalue | largest difference over the family |
|---|---|
| -1.1608 x 10^3 | 2.3417 x 10^-7 |
| -8.6150 x 10^-4 | 4.1418 x 10^-13 |
| 1 - (9.8820 x 10^-10 - i 4.4191 x 10^-5) | 1.5310 x 10^-4 |
| 1 - (9.8820 x 10^-10 + i 4.4191 x 10^-5) | 1.5307 x 10^-4 |
| 0.4689 + i 0.8832 | 3.3372 x 10^-10 |
| 0.4689 - i 0.8832 | 3.3372 x 10^-10 |

### Table 3 (p25): Hill variation orbit coefficients d_p (using m)

| p | d_p |
|---|---|
| 0 | 1 |
| 1 | -2/3 |
| 2 | 7/18 |
| 3 | -4/81 |
| 4 | 19565/62208 |
| 5 | -47161/93312 |
| 6 | -2284055/3359232 |
| 7 | -1152145/1259712 |
| 8 | -65557603/1934917632 |
| 9 | 4005971079828870083/13850679916489605120 |

### Table 4 (p26): Hill variation orbit coefficients c_(n,p) (using m)

Zero entries are printed as 0; rows p = 0, 1 are not tabulated (c_(n,p) = 0 for p < 2, p25). Columns n = -4, -3, -2, -1, 1, 2, 3, 4.

| p | n = -4 | -3 | -2 | -1 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|---|---|
| 2 | 0 | 0 | 0 | -19/16 | 3/16 | 0 | 0 | 0 |
| 3 | 0 | 0 | 0 | -5/3 | 1/2 | 0 | 0 | 0 |
| 4 | 0 | 0 | 0 | -43/36 | 7/12 | 25/256 | 0 | 0 |
| 5 | 0 | 0 | 23/640 | -14/27 | 11/36 | 803/1920 | 0 | 0 |
| 6 | 0 | 1/192 | 299/2400 | -7381/82944 | -30749/110592 | 6109/7200 | 833/12288 | 0 |
| 7 | 0 | 7477/215040 | 56339/288000 | 3574153/2488320 | -1010521/829440 | 897599/864000 | 27943/71680 | 0 |
| 8 | 23/6144 | 65239/627200 | 238200053/1105920000 | 55218889/9331200 | -18445871/6220800 | 237203647/368640000 | 12275527/11289600 | 3537/65536 |
| 9 | 795829/28901376 | 2674679587/14224896000 | 146886277/537600000 | 13620153029/1119744000 | -2114557853/373248000 | -11098919887/14515200000 | 27409853579/14224896000 | 18638507/48168960 |

The signs are those of the text layer; the minus signs before the printed fractions in row 6 to 9 were checked against the page image. No cell was computed by me.

## 4. What is treated, in which model, and what happens

READ unless marked. The paper's own sections 3 and 4.

| CR3BP starting orbit | Resonance | Model | What was done and found |
|---|---|---|---|
| L2 9:2 NRHO (northern; the paper says southern is equivalent by symmetry, p8) | T* = 4 pi/9, a = 9, b = 4 (HR4BP); b = 2 (BCP) | HR4BP (Sect. 3, p8-18) | M_i is identically zero for every 1 <= i <= 9 (all orders the Hill coefficients allow), so there are no simple zeros to choose from. The family was started from the point on the xz-plane with z > 0, with y0 = x0' = z0' = 0 enforced ("HR4BP 9:2 NRHO Family 0"), continued in m from 0 to m_SEM = 0.0808 and on to m = 0.5182. The SEM orbit PO 0 is the family member nearest m_SEM, then corrected with m fixed. Continuation at fixed m_SEM with free X0, in the direction of the smallest singular value (2.7e-11), step 2e-4, traces a closed curve Gamma of initial states; PO F is the member whose state at tau = pi matches PO 0 ("a different periodic orbit that contains the same collection of states", its trajectory "time-shifted by Tg"). Gamma splits into b = 4 segments Gamma_k = phi_(k Tg)(Gamma_0). Orbits plotted have nine crossings of the xz-plane (Fig. 4). |F| <= 1e-10 throughout. The monodromy eigenvalues are "virtually the same for all periodic orbits in this family" (Table 2). |
| same | T* = 1.0808 x 4 pi/9, a = 9, b = 2 | BCP (Sect. 4.1, p19-21) | Continued from eps_B = 0 to 1 starting from the same xz-plane point. 93 multiple-shooting nodes ("the continuation of this SEM PO family in the BCP was much more difficult"; |F| <= 5e-10, mostly 5e-11); sigma_alpha rises from eps_B about 0.2 but is below 1e-6 at eps_B = 1. Same closed curve Gamma with b = 2. The authors note the stability is "consistent with the BCP SEM 9:2 NRHO presented in Table 3 of the work of Boudad et al. [43]". |
| L2 7:2 NRHO | a = 7 (INFERRED: T* = 4 pi/7, b = 4, because each orbit shows seven crossings) | HR4BP (Sect. 4.2, p21-23) | Same behaviour: M_i identically zero up to the 9th order; family computed at m = m_SEM (Fig. 13) and at m = 0.1369 (Fig. 14, "this behavior persists at a value of m that is larger than m_SEM"). sigma_alpha steadily increases after m about 0.0843, "slightly larger than m_SEM"; a clear increase starts at sigma_alpha about 1e-10 (for the 9:2 at about 1e-8). |
| L2 5:2 halo | a = 5 (INFERRED: T* = 4 pi/5) | HR4BP (p21) | "this behavior was far weaker for the 5:2"; no figure. With the Hill coefficients recomputed to P = 15, i_n = 10 and Delta s = pi/10 (p15). |
| L4 4:13.5 long-period planar orbit | T* = 27 pi/4, T = 27 pi, a = 4, b = 27 | HR4BP (Sect. 4.3, p23-24, Fig. 15-16) | Planar. M_i identically zero up to 7th order; M_8(s0, 0) has 216 simple zeros for s0 in [0, 27 pi) separated by pi/8 (COMPUTED: 27 pi/(pi/8) = 216). Continued from s0 = 0, tau0 = 0. "At very small values of m (smaller than m about 0.0083) we see a similar pattern ... However, the continuation for this particular orbit reaches a turning point at a relatively small value of m about 0.0253, and the continuation does not reach the SEM system." |
| L2 3:1 and 4:1 NRHO, L1 1:0.5 vertical | a = 3 (T* = 2 pi/3), a = 4, a = 1 | cited only | i_n values (6, and 4) quoted from [25]; in the BCP "the continuation of the 3:1 NRHO and 4:1 NRHO in the BCP (shown in Fig. 24 and Figs. 25-26 in [14])" show "a clear number of dynamical equivalents" (p8). |
| BCP L2 5:2 halo, 3:1 NRHO, 4:1 NRHO, 7:2 NRHO, L4 3:11 planar | resonance as labelled | BCP Melnikov only (Appendix 2) | M_(1,j) evaluated; j_n = a - 1 (section 2.5). |

Orbit and object families (Sect. 3.2, p10-11): "As the system is periodic with Tg = pi, every state at tau = k Tg, k in Z+ on a periodic orbit can be used as the initial state at tau0 = 0 to produce a new periodic orbit containing the same collection of states as the original periodic orbit, but with a state-time history that is shifted in time by k Tg (mod T) compared to the original orbit." Hence the number of distinct orbits in the family is not an integer count but a one-parameter set whose members are phase shifts of each other; "the full family can be characterized by the periodic orbits between PO 0 and PO F".

The torus interpretation (Sect. 3.3, p12-13): "Examining Figs. 4 and 5, this SEM periodic orbit family behaves like a 2D torus with T0 = pi and rho = 9 pi/2. Notice that states on the invariant curve at tau0 experience no net rotation on the invariant curve after tau = 0 (mod 4 pi), so trajectories on such a torus are POs with T = 4 pi" (the rotation number is rho = 2 pi T0/T*, T0 = Tg; COMPUTED: 2 pi x pi/(4 pi/9) = 9 pi/2). Caveat quoted (p13): "all of the statements we will make in this discussion come with the significant caveat that they apply only to the numerical behavior of these orbits. For example, if the initial states traced out a closed curve that was continuous (like the invariant curve of a 2D torus described in the previous paragraph), then the resulting trajectories would be periodic and technically constitute a continuous 1-parameter family of periodic orbits. This behavior would violate our understanding of periodic structures in these types of systems. So, while all of the orbits we computed between PO 0 and PO F satisfy the periodicity constraint |F| <= 1e-10 in our numerical correction algorithm, some of them may not be truly periodic in the mathematical sense, they may not satisfy F = 0. ... we can only state that this family of SEM periodic orbits exists and behaves like a SEM 2D torus numerically."

Persistence in m (Sect. 3.5, p16-17): the 9:2 family was continued to m = 0.5182 ("this value of m would be achieved if the Earth-Moon barycenter was approximately 0.3656 AU from the Sun"); "sigma_alpha begins to steadily increase after m about 0.1688 about 2 m_SEM" and the m versus x0 hodograph starts to flatten at "x0 about 0.9604 and m about 0.2371". The authors decline to give a single critical value m_cr and list three possible endings (orbits drift apart, bifurcation, turning point in m).

Stability (Sect. 3.6, p17): the paper does not classify the SEM orbits as stable or unstable beyond Tables 1, 2 (one large real pair at -1160.8 and -8.6e-4, one oscillatory pair, one pair at unity). It says the modes of lambda_1, lambda_2, lambda_5 and lambda_6 "are the subject of future work".

Ephemeris (Sect. 3.7, p18): not computed. "Zimovan-Spreen et al. [42] identified a number of multi-year trajectories in the ephemeris model that each correspond to the 9:2 NRHO with different solar phasings. These results provide preliminary evidence that, qualitatively, the structure of the 9:2 NRHO is resistant to the effect of solar phasing in the ephemeris model."

Plain statement on the lead's question (item 3, last sentence): no cycler is treated; no Earth-Moon transfer orbit is treated; no orbit is defined or classified by its lunar flybys. The NRHOs are close-lunar-pass orbits by nature (the 9:2 NRHO is the Gateway orbit) but the paper never gives a perilune radius, flyby count or lunar-encounter statistic for any orbit; the only geometric statement is that "the 9:2 NRHO in the SEM system maintains approximately the same configuration relative to the Moon independent of solar phasing" (Fig. 6, p11, orbits plotted relative to the Moon). The orbits here are all in the L2 or L4 regions of the Earth-Moon system, not the low-energy Earth-Moon cycler families of the project.

## 5. Higher-order Melnikov theory when the first-order function is flat

This bears on the project's weak resonances (8:3 and similar).

READ. The paper's statement of the approach is the sentence quoted in section 2.3: "M_1(s0, tau0) = 0 for any s0 and tau0. When M_i is identically zero for a particular orbit we compute the next order M_(i+1)" (p6). In the HR4BP the first-order term vanishes always (it reflects the time scaling; 2024 digest), so the leading Sun effect on any orbit is order m^2 and M_2 is the first function used; when it is identically zero the sequence is continued to M_3, M_4 and so on up to M_9 with the Hill coefficients used here (M_15 with the P = 15 coefficients for the 5:2 halo). This paper extends the 2024 treatment: p3 quotes the earlier work as having "only presented the Melnikov function up to 3rd-order, and it did not comprehensively address the behavior of CR3BP periodic orbits with periods that are not an integer multiple of pi (e.g., the L2 9:2 NRHO has a period of T* = 4 pi/9 in the CR3BP for analysis in the HR4BP)." Results on how the nonzero order scales (quoted in section 2.5): i_n = 2a in the HR4BP and j_n = a - 1 in the BCP, where a is the numerator of the commensurability aT* = bTg, with an explanation by trigonometric harmonics but no proof ("we have not derived an analytical proof").

Interpretation the paper draws: "if a CR3BP resonant periodic orbit has a sufficiently large value of i_n, then we will find this orbit has the potential to numerically behave like a 2D torus in the perturbed system" (p15) and "We also expect that (in most cases) for CR3BP periodic orbits with a sufficiently small value of i_n, m_cr is small enough to be considered negligible, which is the expected behavior for resonant periodic orbits as stated by MacKay et al. [37], '[the critical parameter] k_C(nu)... is typically zero for rational [frequency] nu'" (p16). For the 9:2 NRHO, "If i_n = 2a holds for the 9:2 NRHO, we would expect i_n = 18", out of reach of the P = 9 coefficients (p15). The physical mechanism proposed (p15): "We hypothesize that the perturbation's effect on the energy balance corresponding to the set of states on the initial invariant curve is the physical mechanism that opens these holes. When the perturbation's effect on the energy balance, represented by the Melnikov function, is numerically undetectable, then the holes are not numerically detectable either."

Caveats (INFERRED, not stated by the authors): (1) the higher-order M_i are evaluated on the unperturbed orbit with the unperturbed velocity; once M_1 is flat a genuine second-order Melnikov theory would include the first-order correction to the orbit and to the energy, and the paper does not claim to do that, relying instead on the assumption quoted in section 2.3 that the leading term represents the energy function; (2) i_n = 2a rests on five orbits plus one recomputation (a = 3, 4, 5 and the exceptional a = 1; a = 9 and 7 unresolved), and the BCP rule j_n = a - 1 on five orbits, both stated as "in almost all cases" or "for most".

Translation to the project's notation (INFERRED). For a project period T*/T_sun = p/q with the Sun's synodic period T_sun: in the BCP Tg = T_sun, so a = q and j_n = q - 1: q = 1, 2 are decided by the quadrupole (j = 1), q = 3 first appears at the octupole (j = 2), suppressed by a further factor r/a_S of order 1e-3 relative to the quadrupole (a_S = 389). This agrees with the 2008 Leiva-Briozzo selection (`docs/notes/2026-10-04-digest-leiva-briozzo-2008-rtbp-to-qbcp-periodic-transfer-orbits.md`): their first-order condition (second harmonic of the Sun angle, quadrupole) exists only for q = 1 or 2, and for q = 3 "no condition is obtained at first order". In the HR4BP, Tg = T_sun/2, so a = 1 for q = 1 or 2 (i_n = 2) and a = 3 for q = 3, so i_n = 6 and Delta s0 = pi/6; for the project's 8:3 resonance, T* = (8/3) T_sun = (16/3) pi, so a = 3, b = 16, the predicted zeros are spaced pi/6 apart, giving 96 zeros over the forced period 16 pi (by the same counting that gives 216 for the 4:13.5 orbit; a prediction, not a result in the paper). The 5/2 resonance (a = 1 in the HR4BP and q = 2 in the BCP) is predicted to have the strong quadrupole-driven structure with zeros pi/2 apart, as for "most orbits evaluated by Brown et al. [25]".

## 6. What this means for `#884`

Published here or in the companion 2024 paper [25] (READ, with quotes):

- Sun phases at which equivalents exist are zeros of a Melnikov-type function of the initial phase s0 (equivalently tau0), evaluated on the unperturbed orbit; for the orbits with a = 1 and the half-period symmetry the zeros are at the symmetric phases (2024 digest, Proposition 3: "any points on that orbit that satisfy the half-period symmetry conditions ... are points where A = 0 provided tau0 = k1 pi/2"). The 2025 paper does not restate that, but uses tau0 = 0 throughout. This is the published form of the project's "zeros at the symmetric Sun phases" observation. Not published here: any planar low-energy Earth-Moon orbit and any statement for the Sun's bicircular phases of the project's families.
- Number of equivalents: "we typically find a limited number of HR4BP dynamical equivalents to a resonant periodic orbit in the CR3BP, e.g., two or four" (p12), tied to Delta s0 = pi/i_n (section 2.5). The paper adds that for large a (9:2 NRHO; partly 7:2) the equivalents become numerically a continuum on a closed curve, not an isolated few. For the project's weak 8:3 resonance the rule predicts i_n = 6 (HR4BP) or j_n = 2 (BCP), that is, many zeros and, by the paper's hypothesis, a weak first-order effect; the number of equivalents is not printed for any orbit of the project's kind, and the rule itself is empirical.
- Weak resonances: published qualitatively and with a formula (section 5), including the statement that larger a gives a weaker Sun effect, from first principles of the harmonics. The project's choice of the quadrupole-free resonances (q = 3) as "weak" is consistent with it.
- Folds and turning points: the paper mentions "a turning point in m" as one reason a family stops (p17) and finds one for the L4 4:13.5 orbit at m about 0.0253 (p23). The companion paper has the northern L2 halo family continuing to the southern one at m = 0 via a turning point in m (2024 digest, item 7). So a fold in the Sun parameter that reconnects two three-body orbits at m = 0 is published for the L2 halo family. This paper does not treat folds and reconnections of Earth-Moon cycler-class orbits, and does not give the Sun-mass value of any fold.
- The Sun's sense: the BCP appendix independently fixes theta_S = -omega alpha + theta_S0, theta_S0 = pi, Earth between the Sun and the Moon at the start, the project's convention after `#891` (section 1.2).

Not published here (so not prior art for the project's results):

- Any cycler, Earth-Moon transfer orbit or low-energy resonant orbit with lunar flybys; any planar Earth-Moon family at C of about 3.0 to 3.2; any family at q = 3 or at the 8:3 or 5/2 resonances with a Sun-driven continuation of such orbits.
- Any initial condition, period (beyond T* = 4 pi/9, 27 pi/4, T = 27 pi), Jacobi constant, perilune radius or stability index (beyond Tables 1, 2).
- A positive control for the corrected modules: there is none in this paper. The only printed orbit numbers are the PO 0 eigenvalues (Table 1) and the near-unity differences (Table 2), which would test a 3D Hill-problem code and not `core/bcr4bp.py`; the BCP constants (a_S = 389.1724, mu_S = 3.2906e5) differ from the project's Andreu-based values, so even a BCP 9:2 NRHO run would not reproduce theirs to many digits.
- Second-order Melnikov theory in the strict sense (section 5), a proof of i_n = 2a or j_n = a - 1, or any treatment of the Sun-mass continuation (the project's) as distinct from the m or eps_B continuations here. In the BCP, eps_B multiplies the whole Sun term at fixed omega, which resembles the project's Sun-mass continuation (INFERRED from Eq. 11d); m in the HR4BP changes the primaries' motion as well and does not.

Note for the project's literature record: this paper and the 2024 paper are the best published description of how Sun-phase zeros and equivalents arise, but for L2 NRHOs, L1 and L4 orbits only. For the Earth-Moon cycler-class families the nearest prior art remains Leiva and Briozzo (2005, 2008) in the QBCP.

## 7. References cited for earlier continuation work in these models (as printed, p33-34)

Bicircular and quasi-bicircular models:

- [12] Gomez, G., Jorba, A., Masdemont, J., Simo, C.: Normal form of the bicircular model and related topics. In: Dynamics and Mission Design Near Libration Points, pp. 53-110. World Scientific, Singapore (2001).
- [13] Jorba-Cusco, M., Farres, A., Jorba, A.: Two periodic models for the Earth-Moon system. Front. Appl. Math. Stat. 4, 32 (2018).
- [14] Boudad, K.K., Howell, K.C., Davis, D.C.: Dynamics of synodic resonant near rectilinear halo orbits in the bicircular four-body problem. Adv. Space Res. 66(9), 2194-2214 (2020).
- [15] Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Families of halo-like invariant tori around L2 in the Earth-Moon bicircular problem. Celest. Mech. Dyn. Astron. 133(4), 16 (2021).
- [16] Simo, C., Gomez, G., Jorba, A., Masdemont, J.: The bicircular model near the triangular libration points of the RTBP. In: From Newton to Chaos: Modern Techniques for Understanding and Coping with Chaos in N-Body Dynamical Systems, pp. 343-370. Springer, Boston (1995).
- [17] Castella, E., Jorba, A.: On the vertical families of two-dimensional tori near the triangular points of the Bicircular problem. Celest. Mech. Dyn. Astron. 76, 35-54 (2000).
- [19] Andreu, M.A.: The quasi-bicircular problem. PhD thesis, Universitat de Barcelona (1998).
- [20] Le Bihan, B., Masdemont, J., Gomez, G., Lizy-Destrez, S.: Invariant manifolds of a non-autonomous quasi-bicircular problem computed via the parameterization method. Nonlinearity 30(8), 3040 (2017).
- [21] Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Invariant manifolds near L1 and L2 in the quasi-bicircular problem. Celest. Mech. Dyn. Astron. 135(2), 15 (2023).
- [22] Andreu, M.A.: Dynamics in the center manifold around L2 in the quasi-bicircular problem. Celest. Mech. Dyn. Astron. 84, 105-133 (2002).
- [43] Boudad, K.K., Howell, K.C., Davis, D.C.: Departure and escape dynamics from the near rectilinear halo orbits in the Earth-Moon-Sun system. J. Astronaut. Sci. 69, 1076-1114 (2022).

Hill and elliptic-circular four-body models:

- [23] Scheeres, D.J.: The restricted Hill four-body problem with applications to the Earth-Moon-Sun system. Celest. Mech. Dyn. Astron. 70(2), 75-98 (1998).
- [24] Sanaga, R.R., Howell, K.C.: Synodic resonant near rectilinear halo orbits in the Hill restricted four-body problem. In: 33rd AAS/AIAA Space Flight Mechanics Meeting (2023).
- [25] Brown, G.M., Peterson, L.T., Henry, D.B., Scheeres, D.J.: Structure of periodic orbit families in the Hill restricted 4-body problem. SIAM J. Appl. Dyn. Syst. 24(1), 346-375 (2025).
- [26] Olikara, Z.P., Gomez, G., Masdemont, J.J.: A note on dynamics about the coherent Sun-Earth-Moon collinear libration points. In: Gomez, G., Masdemont, J.J. (eds.) Astrodynamics Network AstroNet-II, pp. 183-192. Springer, Cham (2016).
- [27] Henry, D.B., Rosales, J., Brown, G.M., Scheeres, D.J.: Quasi-periodic orbits near Earth-Moon L1 and L2 in the Hill restricted four-body problem. AAS/AIAA Astrodynamics Specialist Conference (2023).
- [28] Peterson, L.T., Jorba, A., Brown, G.M., Scheeres, D.J.: Dynamics around the Earth-Moon triangular points in the Hill restricted 4-body problem. Celest. Mech. Dyn. Astron. 136(4), 32 (2024).
- [18] Villegas-Pinto, D., Baresi, N., Locoche, S., Hestroffer, D.: Resonant quasi-periodic near-rectilinear halo orbits in the elliptic-circular restricted four-body problem. Adv. Space Res. 71(1), 336-354 (2022).

Other methodological references: [29] Cenedese, M., Haller, G.: How do conservative backbone curves perturb into forced responses? A Melnikov function analysis. Proc. R. Soc. Lond. A 476(2234), 20190494 (2020). [11] Zimovan-Spreen, E.M., Howell, K.C., Davis, D.C.: Dynamical structures nearby NRHOs with applications to transfer design in cislunar space. J. Astronaut. Sci. 69, 718-744 (2022). [42] Zimovan-Spreen, E.M., Scheuerle, S.T., McCarthy, B.P., Davis, D.C., Howell, K.C.: Baseline orbit generation for near rectilinear halo orbits. AAS/AIAA Astrodynamics Specialist Conference (2023). [37] MacKay, R.S., Meiss, J.D., Percival, I.C.: Stochasticity and transport in Hamiltonian systems. Phys. Rev. Lett. 52(9), 697-700 (1984). [33] Mireles James, J.D.: Celestial Mechanics Notes Set 1: Introduction to the N-Body Problem (2007). [10] Lee, D.E.: White Paper: Gateway Destination Orbit Model: A Continuous 15 Year NRHO Reference Trajectory. NASA document 20190030294 (August 2019).

Slip notes (factual): Eq. 16c as printed has a caret over the first vector in the bracket (see section 2.4); Table 3 and Table 4 are labelled "using m" and are the Hill coefficients for the original m, not the modified parameter M of [25], as the appendix text says (p25), so the 2024 paper's Tables 1 and 2 should not be used with this paper's Eq. 3.
