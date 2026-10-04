# Digest: Scheeres (1998), "The restricted Hill four-body problem with applications to the Earth-Moon-Sun system"

Celestial Mechanics and Dynamical Astronomy 70:75-98 (1998), DOI 10.1023/A:1026498608950. Received 14 October 1997, accepted 13 May 1998. D. J. Scheeres, Iowa State University. 24 PDF pages (journal pp. 75-98, one page per PDF page).
Filed in the private paper corpus as `scheeres-1998-restricted-hill-four-body-problem-earth-moon-sun-cmda-70-75-doi-10.1023-A-1026498608950.pdf`.

Digested 2026-10-04. Status: text layer read in full; the journal pages with the key equations were read as images (p82 eqs. 43-46, p84 eq. 61, p88 Fig. 1 and Fig. 2, p97 appendix eqs. 73-78) and agree with the text layer except where noted.
The remaining figures (Figs. 3 to 12) were not read as images; their captions and the surrounding text are used. READ = on the page (journal page given). COMPUTED = my arithmetic or numerical check on 2026-10-04 (scratch scripts, not committed). INFERRED = my reading.
Related digests: `2026-10-04-digest-grossi-topputo-2025-hr4bp-invariant-tori-collocation.md`, `2026-10-04-digest-sanaga-howell-2025-hill-restricted-four-body-ephemeris-transition.md`, `2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md`.

## 0. What the paper is

READ (abstract, p75). A derivation of the restricted Hill four-body problem (HR4BP) from the general four-body problem by the Hill and restricted approximations, two parameters: the mass ratio nu of the restricted three-body problem and the period parameter m of Hill's Variation orbit.
The restricted three-body problem is recovered as m -> 0 (in the rescaled frame) and Hill's problem as nu -> 0. Application: the planar region near the L4 and L5 positions for 0.0115 <= nu <= 0.0130, 0 <= m <= 0.195: a pi-periodic family (the "L4 and L5 fixed points", Earth-Moon-Sun: half-month period, unstable) and a 2 pi family (one month, stable at the Earth-Moon-Sun values), the manifolds of the unstable pi orbit.
This is the origin paper of the HR4BP that Grossi & Topputo, Sanaga & Howell and Brown et al. all build on. It prints no orbit table: periodic orbits appear only as plotted initial-condition curves (Figs. 2-4, 6-7), so there are no initial conditions to use as tests (section 6).

## 1. Derivation as printed (pp76-81)

### 1.1 General four-body problem (section 2, eqs. 1-13)
READ. Bodies M0 (largest, m0 = 1), M1, M2, M3 with m_j = mu alpha_j, sum alpha_j = 1 (eqs. 3-6). Rotating frame with constant angular velocity V; time scaled so G = 1:
R_i'' + V x 2R_i' + V x (V x R_i) = (1/m_i) grad_i U, U = sum_{i<j} m_i m_j/|R_ij| (eqs. 1, 2). Centre of mass of M1..M3: R_c = (1/mu) sum m_j R_j (eq. 8); R0 = -mu R_c (eq. 9). Shift R_j = R_c + r_j (eq. 10) gives eq. 12 (the R_c equation) and eq. 13 (the r_i equations).
(Eqs. 12 and 13 are partly garbled in the text layer; they are not used below.)

### 1.2 Hill approximation (section 3, eqs. 14-25)
READ. Assumptions: mu << 1, |r_ij| << 1, |R0 - Rc| = O(1) (eqs. 14-16); scaling r_j = mu^(1/3) rbar_j (eq. 17), expansion of the third-body term in powers of mu^(1/3) (eq. 19), dropping mu^(2/3) and higher. Results:
- eq. 20: R_c'' + V x [2 R_c' + V x R_c] = -R_c/|R_c|^3 (two-body problem in a rotating frame);
- eq. 21: relative equations r_i'' + V x [2 r_i' + V x r_i] = -r_i/|R_c|^3 + 3 R_c (R_c . r_i)/|R_c|^5 + sum_k' alpha_k r_ik/|r_ik|^3.
- Circular centre-of-mass orbit: R_c = a + mu^(1/3) r_c, a constant in the rotating frame, V x V x a = -a/|a|^3 (eq. 22), so a . V = 0; take |a| = 1 (this sets the length unit, "to ensure |R_c| = O(1)") and then |V| = 1 (the Sun's mean motion about the system is the unit of angular rate). Eq. 24 is the linear equation for r_c. The eccentricity of this "near circular" orbit is O(mu^(1/3)).
- eq. 25, "the Hill four-body problem": r_i'' + V x [2 r_i' + V x r_i] = -r_i + 3 a (a . r_i) + sum_k' alpha_k r_ik/|r_ik|^3. READ: it has linear momentum and energy integrals but no angular momentum integral.

### 1.3 Restricted approximation (section 4, eqs. 26-34)
READ. alpha_3 -> 0, alpha_1 = 1 - nu, alpha_2 = nu (nu <= 1/2), r = r2 - r1: r'' + V x [2 r' + V x r] = -r + 3 a (a . r) - r/|r|^3 (eq. 31), r1 = -nu r, r2 = (1 - nu) r (eqs. 32, 33). Eq. 34: the equation for M3 under the tide of the Sun and the attraction of M1 and M2 at r3 + nu r and r3 - (1 - nu) r.

### 1.4 Hill's problem for the primaries and the Variation orbit (section 5, eqs. 35-42)
READ. Frame (i, j, k): a = i, V = k. Eq. 39 (centre-of-mass linearised equations) and 40 (their solution). Hill's equations for the primaries' relative vector (xi, eta, zeta) (eq. 41):
xi'' - 2 eta' = 3 xi - xi/(xi^2+eta^2+zeta^2)^(3/2), eta'' + 2 xi' = -eta/(...)^(3/2), zeta'' = -zeta - zeta/(...)^(3/2). (The text layer prints a stray "2" in the denominators; it is the exponent 3/2 split across lines, and the coefficient of the Kepler term is 1, as in eq. 31. Checked in 3.3 below by a numerical test of the Variation orbit.)
The Variation orbit is periodic with period 2 pi m in t: xi + i eta = sum_{n=-inf}^{inf} a_n(m) exp(i (2n+1) t/m) (eq. 42, written as cos and sin; zeta = 0). Stability: "stable if retrograde (m < 0) and stable for direct orbits in the interval 0 < m < 0.19510486..." (Henon 1969, p82).

## 2. The restricted Hill four-body problem as used (section 6, eqs. 43-64, pp82-84)

READ (images p82, p84). Transform to a frame rotating at absolute rate (1 + 1/m) about V and rescale: rho3(t) = a0(m) r(t/m) (eq. 43), t = m tau (eq. 44), a0(m) = m^(2/3)(1 - (2/3) m + (7/18) m^2 - (4/81) m^3 + ...) (eq. 45). Then eq. 46:
r'' + (1+m) V x [2 r' + (1+m) V x r] = -m^2 r + 3 m^2 a (a . r) - (m^2/a0^3)(1-nu)(r + nu rho/a0)/|r + nu rho/a0|^3 - (m^2/a0^3) nu (r - (1-nu) rho/a0)/|r - (1-nu) rho/a0|^3.
New axes (i_m, j_m, k): i = cos tau i_m - sin tau j_m, j = sin tau i_m + cos tau j_m (eqs. 47, 48); a = cos tau i_m - sin tau j_m (eq. 51); the primaries' separation rho/a0 = (1 + xibar) i_m + etabar j_m (eq. 52) with
xibar(tau; m) = sum_{n>=1} (a_n/a0 + a_{-n}/a0) cos 2n tau, etabar(tau; m) = sum_{n>=1} (a_n/a0 - a_{-n}/a0) sin 2n tau (eqs. 53, 54).
Scalar equations (eq. 55): x'' - 2(1+m) y' = V_x, y'' + 2(1+m) x' = V_y, z'' = V_z, with (eq. 56)

    V = (1/2)(1 + 2m + (3/2) m^2)(x^2 + y^2) - (1/2) m^2 z^2 + (3/4) m^2 ((x^2 - y^2) cos 2 tau - 2 x y sin 2 tau) + (m^2/a0^3)((1 - nu)/R_{1-nu} + nu/R_nu),
    R_{1-nu} = sqrt([x + nu(1 + xibar)]^2 + [y + nu etabar]^2 + z^2), R_nu = sqrt([x - (1 - nu)(1 + xibar)]^2 + [y - (1 - nu) etabar]^2 + z^2) (eqs. 57, 58).

READ: periodic in tau with period pi; m -> 0 gives the restricted three-body problem; "the distance of the primaries in the Variation orbit is on the order of m^2" (Wintner), "by scaling our distance by a0^3(m) ~ m^2" [sic: the text says a03; the scale is a0 ~ m^(2/3) in length with a0^3 ~ m^2 appearing in the force; the sentence is loose] "the distance of the Sun ... is on the order of 1/m^2" (see section 3.3 for a precise statement).
Symmetry (eqs. 59, 60): (x, y, tau) -> (x, -y, -tau) with (x', y', x'', y'') -> (-x', y', x'', -y''). Consequences: the history of (x0, -y0, -x0', y0') at -tau0 mirrors that of (x0, y0, x0', y0') at tau0; a periodic orbit of period n pi maps to another of the same period and stability; stable and unstable manifolds exchange. The L5 orbits follow from the L4 ones by eq. 67.
Expansion in m (eq. 61, image p84): V = (1 + 2m + (3/2) m^2) V0 - (1/2) m^2 z^2 + m^2 cos 2tau (3/4)(x^2 - y^2) + m^2 cos 2tau nu(1-nu) [(x+nu)/r_{1-nu}^3 - (x-(1-nu))/r_nu^3] - m^2 sin 2tau [(3/2) x y + (11/8) nu(1-nu)(1/r_{1-nu}^3 - 1/r_nu^3) y] + O(m^3), with V0 the restricted three-body potential (eq. 64) and r_{1-nu}, r_nu the fixed-primary distances (eqs. 62, 63).
COMPUTED check of eq. 61 against eq. 56 (nu = 0.0122, three points, the printed a_n/a0 to m^6): the difference is 1.2e-7 to 4.9e-7 at m = 0.02, 1.7e-7 to 5.3e-6 at m = 0.04, 2.0e-5 to 6.7e-5 at m = 0.08, growing as roughly m^3.5, consistent with the stated O(m^3). So eq. 61, including the 11/8 and 3/2, is consistent with eq. 56.

## 3. Parameters, and the Sun's position

### 3.1 Printed values for the Earth-Moon-Sun system
READ. nu = 0.0122 ("approximately", p86; Grossi & Topputo print mu = 0.01215); m = 0.0808 ("of the Earth-Moon-Sun system is approximately 0.0808", p86); parameter box 0.0115 <= nu <= 0.0130 (about 6 percent around the Earth-Moon value), 0 <= m <= 0.195 (the Variation orbit becomes unstable at 0.19510486...). The Sun's actual orbit eccentricity 0.0167 against mu^(1/3) = 0.0145, "so this is a valid assumption" (p85). Inclination of the Earth-Moon plane to the ecliptic about 5 degrees; mu^(1/3) rad is "approximately 1 degree" (p85).
COMPUTED: with mu = (m1+m2)/m0 = 3.0404e-6 (Grossi & Topputo's M), mu^(1/3) = 0.014487 rad = 0.83 degrees, matching the printed 0.0145 and "approximately 1 degree". Scheeres does not print mu's numerical value; his mu is Grossi & Topputo's M and his nu is their mu.

### 3.2 Which Sun term the model has
READ. The Sun enters through (i) the tide 3 m^2 a(a . r) - m^2 r of eq. 46, i.e. the potential terms (1/2)(2m + (3/2) m^2) rho^2 ... and (3/4) m^2 ((x^2 - y^2) cos 2tau - 2xy sin 2tau) of eq. 56, and (ii) the coefficient m^2/a0^3 multiplying the Earth and Moon attraction (which carries the time scaling). The Sun's position vector never appears in the equations of motion.

### 3.3 The Sun's position: the Grossi-Topputo r0 = M^(-1/3) a0 slip is settled
Scheeres prints no formula for the Sun's position. It follows from his definitions (COMPUTED / INFERRED derivation):
- R0 = -mu R_c (eq. 9) and R_c = a + O(mu^(1/3)) with |a| = 1 (eq. 22 and following): the Sun is at distance (1 + mu)|R_c| = 1 (to the order kept) from the three-body barycentre, in the direction -a.
- The Hill coordinates are physical offsets from R_c divided by mu^(1/3) (eq. 17). Hence a physical separation of the primaries of D (in units of the Sun-system distance) corresponds to a Hill-coordinate separation D/mu^(1/3). The primaries' separation in Hill coordinates is |rho| ~ a0(m) (eq. 52: rho/a0 = (1 + xibar, etabar) has unit length at leading order). So D = mu^(1/3) a0(m), and the Sun's distance in units of the primaries' mean separation is 1/D = mu^(-1/3)/a0(m).
So the Sun is at mu^(-1/3)/a0(m) (division, not multiplication), and Grossi & Topputo's printed r0 = M^(-1/3) a0(m) (cos t, -sin t, 0) is a slip. COMPUTED with M = 3.0404e-6, m = 0.0808: M^(-1/3) = 69.0277; a0 by the printed series to m^3 = 0.177301 (first order 0.176832); M^(-1/3)/a0 = 389.32 (390.36 with the first-order a0); M^(-1/3) a0 = 12.24. The project's `_ANDREU_A_S` = 388.8111 is 0.13 percent below 389.32 (0.40 percent below the first-order value). That residual is within what the Hill approximation, the rounding of m and M, and the non-circularity of the Sun's orbit can account for; it is not resolved by the paper. Kepler check (INFERRED): n'^2 a_S^3 = m0, n^2 a^3 = m1 + m2 gives a/a_S = mu^(1/3) (n'/n)^(2/3) = mu^(1/3) m^(2/3)(1 + ...) = mu^(1/3) a0(m), the same relation.
Also INFERRED: the scale statement of p83 ("the distance of the Sun ... is on the order of 1/m^2") is true in the sense of the Hill-force prefactor m^2/a0^3 ~ 1 but is loose as a distance; the exact statement is mu^(-1/3)/a0(m), which is of order mu^(-1/3) m^(-2/3).

### 3.4 Sun sense and phase
READ. a = cos tau i_m - sin tau j_m (eq. 51) and i = cos tau i_m - sin tau j_m (eq. 47): in the frame (i_m, j_m) that rotates with the Earth-Moon mean positions, the Sun-system line i has angle -tau. So the Sun direction regresses (clockwise) in the Earth-Moon rotating frame, rate 1 in tau. The tide is invariant under a -> -a, so only the quadrupole sense matters; the printed term (3/4) m^2 ((x^2 - y^2) cos 2tau - 2 x y sin 2tau) equals the quadrupole (3/2) m^2 (x cos tau - y sin tau)^2 up to a multiple of x^2 + y^2 (COMPUTED), the Sun on the line at angle -tau. This is the same sense as `core/bcr4bp._sun_position` (theta = theta_sun0 - omega_sun t, #891) and as Grossi & Topputo's (cos t, -sin t) direction.
INFERRED: a points from the Sun towards the system barycentre (R_c is measured from the centre of mass of all four bodies, which is nearly the Sun), so the Sun lies in direction -a, at angle pi - tau, not tau-reversed; Grossi & Topputo's (cos t, -sin t) is the same line, displaced by pi, so on the opposite side. The tide cannot tell the two apart (period pi), so the HR4BP Sun has no physical phase beyond tau0 mod pi. At tau = 0 the Variation orbit has xibar about -0.0072, i.e. the Moon is on the +x_m axis (rho = r2 - r1, eq. 30), and the Sun line is along x: full or new Moon at tau = 0 (full Moon if the Sun is on -x, INFERRED).
Frame orientation: the Earth is at x = -nu(1 + xibar) < 0 and the Moon at x = (1 - nu)(1 + xibar) > 0 (eqs. 57, 58): the same orientation as the project's convention (Earth at -mu, Moon at 1 - mu), unlike the Jorba-Rosales papers (rotated by pi).

## 4. The coefficient expansions: do they give the a_n(m) of Sanaga & Howell

READ (appendix, p97-98, image p97): "the Fourier coefficients of the Variation orbit, to the order m^6 [7]" (Wintner 1947, pp379-410):
- a0 = m^(2/3)(1 - (2/3) m + (7/18) m^2 - (4/81) m^3 + ...) (eq. 73) (only to m^3 inside the bracket);
- a1/a0 = (3/16) m^2 + (1/2) m^3 + (7/12) m^4 + (11/36) m^5 - (30749/110592) m^6 - ... (eq. 74);
- a_{-1}/a0 = -(19/16) m^2 - (5/3) m^3 - (43/36) m^4 - (14/27) m^5 - (7381/82944) m^6 + ... (eq. 75);
- a2/a0 = (25/256) m^4 + (803/1920) m^5 + (6109/7200) m^6 + ... (eq. 76);
- a_{-2}/a0 = 0 m^4 + (23/640) m^5 + (299/2400) m^6 + ... (eq. 77);
- a3/a0 = (833/12288) m^6 + ... (eq. 78); a_{-3}/a0 = (1/192) m^6 + ... (eq. 79); |a_{+-n}/a0| < O(m^6) for |n| > 3 (eqs. 80, 81).
The text layer agrees with the image digit for digit.

Answer to the question: Scheeres 1998 gives the a_n(m) to order m^6 only (for a0 only to m^3, then the series is cut by "..."), from Wintner. It does NOT give the m^9 expansion that Sanaga & Howell (their N3) attribute to Olikara & Scheeres 2017. The two are consistent, not in conflict: Grossi & Topputo's "up to order 6 in m, as in Reference 16" is exactly this appendix.
Independent confirmation (COMPUTED, symbolic): the Brown et al. 2024 appendix tables (Tables 1 and 2 of the Brown digest, order m^9, in the variable M = m/(1 - m/3)), re-expanded in m with exact rational arithmetic, reproduce every coefficient above, for all six ratios to m^6 and a0 to m^3, exactly. The next a0 term (not printed by Scheeres) is +(19565/62208) m^4 = 0.3145 m^4 (then -(47161/93312) m^5), from the Brown table. So the project can take the m^9 coefficients from the Brown digest's tables, now cross-validated against this paper to m^6.
Physical consequence test (COMPUTED): the printed a_n/a0 satisfy Hill's equations, eq. 41, as a Fourier series in the complex form z'' + 2i z' = 3 Re z - z/|z|^3 (derivative in t = m s, z = sum a_n exp(i(2n+1)s), n from -3 to 3). The maximum residual relative to the largest term is about 1.5e-7 at m = 0.02 and 3.8e-5 at m = 0.0808, scaling as m^4; the m^4 excess is the missing m^4 term of a0 (fitting a0 gives c4 = 0.27 to 0.32, against 0.3145 from the Brown table). This is a source-independent check on the printed coefficients and on the sign and coefficient 1 of the Kepler term in eq. 41.

## 5. Stability and periodic-orbit numbers (section 7, pp85-96)

READ unless marked.
- The L4/L5 positions are equilibria only at m = 0 (restricted three-body problem), stable for nu < 0.0385... and not equal to 0.0243... and 0.0135... (the Routh limit and the 2:1 and 3:1 resonances; Marchal 1990). COMPUTED from the natural frequencies (eqs. 68-71, ω² = (1/2)[1 -+ sqrt(1 - 27 nu + 27 nu^2)]): Routh limit (1 - sqrt(23/27))/2 = 0.038521; omega2/omega1 = 2 at nu = 0.024294; omega2/omega1 = 3 at nu = 0.013516. All three agree with the printed digits. The asymptotics omega1 ~ (3/2) sqrt(3 nu(1 - nu)) and omega2 ~ 1 - (27/8) nu(1 - nu) (eqs. 69, 71): at nu = 0.0122 the exact values are 0.29887 and 0.95429 against 0.28521 and 0.95933. The printed eqs. 68 and 70 have the outer square root lost in the text layer; the exact form is as given here (it reproduces the printed 0.0385, 0.0243, 0.0135).
- The pi-periodic family continues from the equilibria (frequencies scaled by (1 + m)); the 2:1 resonance (multiplier -1) where omega2 (1 + m) = 1 gives m_cr ~ (27/8) nu (1 - nu) (eq. 72). COMPUTED: 0.04067 at nu = 0.0122. Fig. 1 (image, p88): the stable/unstable boundary is a nearly horizontal line at m about 0.04 (graph read, good to about 0.002) over 0.0115 <= nu <= 0.013, the Earth-Moon-Sun point (m = 0.0808, nu = 0.0122) well inside the UNSTABLE region. The paper says the plotted eq. 72 is "indistinguishable" from the numerical boundary at that resolution; the graph-read 0.04 agrees with 0.0407. (The exact omega2(1+m) = 1 root, 0.0479, is not the criterion; the paper's frequency scaling is first order in m.) Instability is single: a pair of multipliers stays on the unit circle; all members have out-of-plane stability.
- Earth-Moon-Sun values (nu = 0.0122, m = 0.0808): the pi orbit has period of half a month and is unstable; the 2 pi orbit has period of one month and is stable (p91). No numerical multiplier is printed for it.
- 2 pi family: stable for most of the box; members with nu > 0.01191 lose stability through the -1 multiplier as their amplitude regrows; some (nu < 0.01250) regain it before m = 0.195; a stable 4 pi family is born where the 2 pi family loses stability and ends where it regains it (pp90-91). For smaller nu the 2 pi family may terminate on the pi family. The 2 pi orbit's stability change at nu = 0.0122 occurs at m about 0.145 (p95). At small m the family resembles a horseshoe orbit (conjectured link, not demonstrated). At minimum size it is "approximately a 2:1 ellipse" about the primary distance long (p90).
- Manifolds of the unstable pi orbit (nu = 0.0122, m_cr < m < 0.195): the stable and unstable vectors are anti-parallel at m_cr, 90 degrees at m about 0.05, minimum about 30 degrees at m about 0.095, about 50 degrees at m = 0.195 (Fig. 8, graph not read). The multipliers start at -1 at m_cr "and monotonically change to -1.7, -0.59 at m = 0.195" (printed as "characteristic exponents lambda"; COMPUTED: 1.7 x 0.59 = 1.003, so these are a reciprocal pair of multipliers, the printed numbers being multipliers, not exponents). The stable-mode frequency, measured as the angle on the unit circle, rises from 0.31 pi at m_cr to 0.35 pi at m = 0.195 and passes pi/3 at m = 0.10538 (a possible doubly unstable 6 pi family). The manifolds appear trapped for m_cr < m < 0.11 (including the Earth-Moon-Sun value), break down from m about 0.11 and escape rapidly by m = 0.12. These are numerical observations over finite times; the paper says there is no KAM trapping because there is no integral.
- The earlier models (Schechter 1968, Kamel & Breakwell 1970, Kolenkiewicz & Carpenter 1968, which keep the Sun distance as a parameter and have a one-month period): their unstable small orbit corresponds to two small ellipses each traversed in a half month, i.e. the pi orbit, a bifurcation product of the half-month orbit as the length parameter is introduced; their two stable one-month orbits 180 degrees apart correspond to the 2 pi orbit.

## 6. Printed numbers usable as sourced tests

| ID | Quantity | Printed value | Where | Test note |
|---|---|---|---|---|
| S1 | a0(m)/m^(2/3) | 1 - 2m/3 + 7m^2/18 - 4m^3/81 | eq. 45, p82; eq. 73 | exact rationals; next term 19565/62208 not printed here |
| S2 | a1/a0, a_{-1}/a0, a2/a0, a_{-2}/a0, a3/a0, a_{-3}/a0 | eqs. 74-79 in section 4 above | p97-98 | exact rationals to m^6; ALSO cross-checked against Brown 2024 (section 4) |
| S3 | Variation orbit stable | 0 < m < 0.19510486... | p82 (Henon 1969) | monodromy test of the Variation orbit in Hill's problem; no independent value in the project |
| S4 | m for Earth-Moon-Sun | 0.0808 | p86 | against 1/omega_S - 1 = 0.08085 of `bcr4bp._ANDREU_OMEGA_S`; printed to 3 s.f. |
| S5 | nu | 0.0122 (range 0.0115 to 0.0130) | p86 | printed to 3 s.f. (the 0.01215 of Grossi & Topputo is the same quantity) |
| S6 | Restricted-problem stability limits | 0.0385..., 0.0243..., 0.0135... | p87 | reproduced here from eqs. 68-71: 0.038521, 0.024294, 0.013516 (a self-contained test) |
| S7 | m_cr ~ (27/8) nu (1 - nu) | eq. 72 | p87 | 0.04067 at nu = 0.0122; Fig. 1 graph read about 0.04 |
| S8 | pi family at nu = 0.0122, m = 0.0808 | unstable, singly; out-of-plane stable; half-month | pp88, 91 | qualitative check (multiplier classification) |
| S9 | 2 pi family at nu = 0.0122, m = 0.0808 | stable, one month | p91 | qualitative |
| S10 | 2 pi family: loses stability for nu > 0.01191, some regain for nu < 0.01250; 4 pi family born | pp90-91 | qualitative thresholds in nu |
| S11 | Multipliers of the pi orbit | -1 at m_cr; -1.7 and -0.59 at m = 0.195 (nu = 0.0122) | p92 | two digits only; product 1.003 |
| S12 | Stable-mode angle | 0.31 pi at m_cr, 0.35 pi at m = 0.195, pi/3 at m = 0.10538 | p92 | pi/3 is a five-digit root crossing, useful as a continuation target |
| S13 | Manifold angle | 90 degrees at m about 0.05, minimum about 30 degrees at about 0.095, about 50 degrees at 0.195 | p91-92 | graph-level |
| S14 | Symmetry | (x, y, tau) -> (x, -y, -tau); (x', y') -> (-x', y') | eqs. 59, 60, 67 | exact; structural test of any HR4BP RHS |
| S15 | m -> 0 | eq. 55 reduces to the restricted three-body problem | p83 | structural |
| S16 | Neglected | inclination about 5 degrees; eccentricity 0.0167 against mu^(1/3) = 0.0145 | p85 | context |
| S17 | Variation-orbit consistency (COMPUTED) | residual of eq. 41 scales as m^4 with the printed a0 | section 4 | optional, computed by us against a published equation |

Not printed: any periodic-orbit initial condition, period in tau as a number other than pi and 2 pi, Jacobi-type constant (there is none for m > 0), monodromy matrix, or manifold coordinate. The figures (Figs. 2-4, 6-7, 9-12) plot them. The only multiplier digits are S11.

## 7. Reconciliation with project code

- `core/bcr4bp.py` (BCP): point-mass Sun at a_S = 388.8111 EM units, direct and indirect terms, regressing (theta_sun0 - omega_S t, #891), time unit 1/omega_S-based. The HR4BP is the large-a_S limit with a quadrupole tide. Same Sun sense (section 3.4). Same Earth-left, Moon-right orientation. Different time unit (Hill: one synodic month = 2 pi; ratio 1 + m = 1.0808, see the Sanaga digest) and the Hill model has the pulsating Variation-orbit primaries.
- `core/qbcp.py` (QBCP): coherent, primaries from Andreu's coupled solution; no Hill content. The QBCP and HR4BP are both coherent models; the HR4BP drops the 1/a_S terms the QBCP keeps.
- There is no HR4BP code in the project. Eqs. 55-58 with eqs. 45, 53, 54 and the a_n of section 4 are a complete RHS; the Jacobian needs d/dx of eq. 56 only, with the xibar, etabar series as known functions of tau.
- The project's Sun conventions are not contradicted by this paper (sense, orientation); the Sun position formula in section 3.3 can be tested as a constant (mu^(-1/3)/a0) against `_ANDREU_A_S` at the 0.13 percent level.
- Conflicts found: none with project code. With the literature: Grossi & Topputo's Sun position (settled, section 3.3). Sanaga & Howell's `a_n` order m^9 attribution is to a later paper, not this one (section 4).

## 8. Techniques applicable to the project's problems

#926 (Hill four-body follow-ups, registered; its (c) says do NOT build an HR4BP module now and, if wanted, get the coefficients from Olikara & Scheeres 2017, not held):
1. The coefficients are now sourced without that paper: this paper (to m^6, cross-validated) and Brown 2024's tables (m^9, in the Brown digest). (c)'s "get Olikara & Scheeres" can be dropped; the remaining unheld item is only a check of the Brown table's higher orders against an independent source (the m^4 a0 coefficient is confirmed by the Hill-equation fit, section 4).
2. If an HR4BP module is built, the self-contained tests available from this paper alone: S1, S2, S6, S7, S14, S15, the eq. 61 consistency test, the Hill-equation residual test (S17), and the qualitative stability of the pi and 2 pi families (S8, S9, S11, S12). A reproduction of the Routh/resonance values (S6) needs no integration at all.
3. #926(a) (sourced-constants module): add the Sun position mu^(-1/3)/a0(m) as a formula test with the 0.13 percent residual recorded, and the m and omega_S relation.
4. #926(b) (multiplier root-of-unity monitor): this paper gives a ready positive control for the monitor on a Hill-model periodic orbit: the pi orbit's stable pair crosses pi/3 at m = 0.10538 (a 6th root of unity) and the multiplier crosses -1 at m_cr about 0.04 (eq. 72). Both are continuation targets with printed values, in the same family continued in m at nu = 0.0122.

#902 (rebuild the bicircular tiers on a real orbit): the HR4BP has no use as a control for the BCP orbits (different model). It does give a model-independent check of the Sun's strength: the Hill tide strength m^2 and the Sun distance mu^(-1/3)/a0 are consistent with a_S = 388.8 to 0.13 percent, so a BCP orbit whose Sun phase and strength are right should differ from its Hill-model counterpart only at O(1/a_S) (INFERRED). That is a cross-model sanity band, not a golden value.

#905 and #884 (Sun-forced Earth-Moon search rerun): (1) The pi-periodicity and the commensurability of orbit periods with the forcing period (periodic orbits exist only for periods that are multiples of pi in tau) restate the closure-over-p-Sun-periods check; in the Hill model a period of k months in the BCP is 2 pi k in tau and a multiple of the BCP Sun period 2 pi/omega_S, and the half-month tide period pi means that a BCP orbit closing over an odd number of half-Sun-periods is a distinct case. (2) The 2:1 resonance (multiplier -1) at m_cr = (27/8) nu (1 - nu), and the family-walk behaviour in m, are what the family walk's bifurcation detection should find in a Hill-type continuation; the printed 4 pi birth and the 6 pi possibility (pi/3 at 0.10538) are cases for the "one orbit per symmetry class" and "bifurcation detection" items of #905. (3) The symmetry (59) and its use for the mirrored L4/L5 orbits is the analogue of the project's symmetric-orbit classes; it is valid for tau0 = 0 or pi/2 here (the tide is even in tau, the primaries' Variation orbit is symmetric about tau = 0). (4) The paper's own caution (p85-96) that the Variation orbit is only a qualitative Earth-Moon motion and the true system is quasi-periodic when the Earth-Moon eccentricity is added matches the #884 reviewers' point that periodic closures in a periodic model are not physical orbits; the torus formulation of the Grossi digest is the next step.
Neither task is advanced by the L4/L5 results themselves: the project's cyclers are not near the triangular points, and the paper treats only planar L4/L5 families. Do not read them as an Earth-Moon cycler result.

## 9. Recommended follow-ups (not registered)

1. Record the Sun-position answer and the a_n(m) answer in the #926 entry (done in the digests; the ledger entry is the coordinator's).
2. If #926 builds an HR4BP module: take a0 and the a_n from the Brown tables (m^9), pin S1 and S2 as exact rationals against this paper, and add the eq. 61 and Hill-equation residual tests.
3. Reproduce S6, S7 and the pi-orbit multiplier crossings (m_cr, pi/3 at m = 0.10538) as the first quantitative controls of an HR4BP periodic-orbit finder; first digits graph-read from Fig. 1 only to 0.04, so use eq. 72 as the target, not the plot.
4. Read Figs. 2-4 and 6-7 as images and digitise the L4/L5 pi and 2 pi initial conditions if a stricter reproduction is wanted (tolerance of a plot, a few percent; mark as digitised). Not done here.
5. Read Fig. 8 and the multiplier curves if S11 to S13 are to be used beyond two digits.
6. The "mu^(-1/3)/a0 versus 388.81" 0.13 percent: if the Sun strength matters for a Hill-versus-BCP comparison, compute the O(mu^(2/3)) Hill correction before attributing the gap.
