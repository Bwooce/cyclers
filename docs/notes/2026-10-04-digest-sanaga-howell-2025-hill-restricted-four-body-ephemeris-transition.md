# Digest: Sanaga & Howell 2025, "Leveraging the Hill restricted four-body problem to investigate the ephemeris transition characteristics in the Earth-Moon L2 halo orbit region"

Astrodynamics 9(5):785-805 (2025), DOI 10.1007/s42064-024-0250-4. Received 18 June 2024, accepted 31 October 2024. Open access (CC BY 4.0). Purdue University.
Filed in the private paper corpus as `sanaga-howell-2025-hill-restricted-four-body-problem-ephemeris-transition-astrodynamics-9-785-doi-10.1007-s42064-024-0250-4.pdf`
(21 PDF pages = journal pp. 785-805; text layer present, 12,508 words by `pdftotext -layout | wc -w`).
Status: READ in full from the text layer, including every equation, both tables and the caption and axis text of all 21 figures. The figure images themselves were NOT viewed, so no
number below was read off a plot except where the axis text states it; those are marked "axis". Page references are journal pages. COMPUTED = my own arithmetic on 2026-10-04. INFERRED = my reasoning, not printed.
Related digests (not repeated): `2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md` (the printed Hill-model numbers, m=0.0808, mu=0.0122),
`2026-10-04-digest-boudad-howell-davis-2020-synodic-resonant-nrho-bcr4bp.md`, `2026-10-04-digest-park-howell-2024-transition-challenging-regions-er3bp-l2-halo.md`,
`2026-10-04-digest-singh-park-howell-2026-aas-26-654-l2-families-intermediary-models.md`, `2026-10-04-digest-jorba-cusco-farres-jorba-2018-two-periodic-models.md`.

## 1. What the paper is, and is not
READ (abstract, p.785): the Earth-Moon L2 halo family is split into three regions by how its members transition to a higher-fidelity ephemeris model (HFEM): a regular halo region, an "interface region"
(called the "transition region" by Davis et al. 2017, AAS 17-826) around the 3:1 synodic resonant orbit, and the NRHO region (9:2 is the Gateway orbit). Members of the interface region are hard to carry to the
HFEM and the result depends on epoch. The paper does NOT carry anything to the HFEM itself; it builds periodic-orbit families of resonant L2 halos in two intermediate models (the Hill restricted
four-body problem, HR4BP, and a reduced version with the Sun's net acceleration removed, RHR4BP), continues them in the mass-ratio-like parameter m and in a Sun-strength homotopy parameter gamma, and finds that:
(a) broken bifurcations are ubiquitous in the HR4BP families (they link the halo family to nearby higher-period families), (b) removing the Sun's net acceleration (gamma 1 to 0) perturbs these
broken bifurcations and DESTROYS some family branches in the interface region, but not in the regular-halo or NRHO regions, (c) this is offered as a dynamical reason for the HFEM transition difficulty.
Its own caveat (p.803): "Preliminary research by the authors suggests that the RHR4BP provides a better initial guess (epoch selection and initial conditions) for transition to the HFEM"; this is not shown here.
There is NO table of initial conditions, periods, Jacobi constants or multipliers anywhere in the paper (see section 5). It is a method-and-structure paper with figures, not a data paper.
Scope relative to the project: Sun-Earth-Moon only, L2 halos only, no cycler content.

## 2. Models (pp.787-789)
### 2.1 CR3BP (p.787) READ
Standard: characteristic mass m* = M1+M2, length l* = primary separation, time t* = sqrt(l*^3/(G m*)); x-axis M1 to M2, z along angular momentum; nu = M2/(M1+M2);
x'' = 2y' + U*_x, y'' = -2x' + U*_y, z'' = U*_z, U* = (x^2+y^2)/2 + (1-nu)/sqrt((x+nu)^2+y^2+z^2) + nu/sqrt((x-1+nu)^2+y^2+z^2). Earth at (-nu,0,0), Moon at (1-nu,0,0). Matches `core/cr3bp.py` conventions.

### 2.2 HR4BP (pp.787-788, Eqs. 1-4) READ
Assumptions (p.787): (i) M0 (Sun) dominates, M0 >> M1, M2, M3; (ii) M1, M2, M3 are near each other; (iii) their centre of mass is a finite distance from M0; (iv) M3 is infinitesimal. The relative motion of M1, M2 and M0
is the Hill three-body problem (H3BP); the HR4BP is COHERENT because M1-M2 motion is taken from a particular solution of the H3BP, the "variation orbit" (Hill's lunar-theory orbit), written as a Fourier series in a
parameter m ("2 pi m is the period" of the family member in the Wintner/Scheeres parameterisation). Scheeres 1998 (CMDA 70:75) is the source of the model. Frame: the "M1-M2 rotating frame", origin at the M1-M2
barycentre B1, for Sun-Earth-Moon "consistent with the Earth-Moon rotating frame". The text prints the frame's angular velocity as "(1 + 1/m)", while Eq. 1 carries the Coriolis factor (1+m); the two are not the
same expression and the text does not reconcile them. INFERRED: the equation is authoritative and the "1/m" is a typesetting slip, given t* = (1+m)/n (below). Do not code the "1/m" form.

Eq. 1 (p.787): x'' = 2(1+m) y' + V_x, y'' = -2(1+m) x' + V_y, z'' = V_z, with the time-periodic potential (Eq. 2):
  V = (1/2)(1 + 2m + (3/2) m^2)(x^2+y^2) - (1/2) m^2 z^2 + (3/4) m^2 [ (x^2 - y^2) cos 2 tau - 2 x y sin 2 tau ] + (m^2/a0^3) [ (1-nu)/R_{1-nu} + nu/R_nu ],
  R_{1-nu} = sqrt( [x + nu(1+xi_bar)]^2 + [y + nu eta_bar]^2 + z^2 ),  R_nu = sqrt( [x - (1-nu)(1+xi_bar)]^2 + [y - (1-nu) eta_bar]^2 + z^2 ),
  xi_bar(tau; m) = sum_{n>=1} ( a_n(m)/a_0(m) + a_{-n}(m)/a_0(m) ) cos 2 n tau,   eta_bar(tau; m) = sum_{n>=1} ( a_n(m)/a_0(m) - a_{-n}(m)/a_0(m) ) sin 2 n tau   (Eqs. 3, 4).
The coefficients a_n(m) are evaluated to order m^9 (Olikara & Scheeres 2017, ref [19]), "a good assumption when m is small" and valid for Earth-Moon with m about 0.0808. a_0 is not defined in this paper
(it is the leading Fourier coefficient of the variation orbit, defined in Scheeres 1998; see the Brown 2024 digest). The equations are periodic in tau with period 2 pi (p.792, "multiples of 2 pi"), so tau is the synodic phase.
Two free parameters, nu and m. For Sun-Earth-Moon (p.788): nu = M2/(M1+M2) about 0.01215 and m = (n'/n)/(1 - n'/n) about 0.0808, where n' is the mean motion of the Sun about B1 (the Earth-Moon barycentre) and n the Moon's about the Earth.
Characteristic quantities (Eq. 5, p.788): m* = M0 (the Sun's mass), l* = ls * a0 * ((M1+M2+M3)/M0)^(1/3) with ls the B2-B1 distance (B2 the barycentre of all four), t* = m/n' = (1+m)/n.
So one nondimensional time unit is (1+m)/n = 1/(n - n') = 1/(synodic rate), and 2 pi in tau is ONE SYNODIC MONTH (consistent with the paper's 2 x 27.32 sidereal days becoming 2 x 29.5 synodic days at m = 0.0808, p.795).

### 2.3 Where the Sun's sense is in Eq. 2 (COMPUTED, my derivation, 2026-10-04)
The Hill tidal potential of the Sun is (3/2) n'^2 (r . s_hat)^2 minus (1/2) n'^2 z^2, with s_hat the unit Sun direction in the rotating frame and n' t* = m. With s_hat = (cos tau, -sin tau) (a Sun at angle theta0 - tau, i.e.
REGRESSING in the Earth-Moon rotating frame), (r . s_hat)^2 = (1/2)(x^2+y^2) + (1/2)[(x^2-y^2) cos 2 tau - 2 x y sin 2 tau], which reproduces the printed (3/4) m^2 [(x^2-y^2) cos 2 tau - 2 x y sin 2 tau] and the (3/4) m^2 (x^2+y^2) part of the
(3/2) m^2 coefficient. A Sun advancing at +tau would give +2 x y sin 2 tau. So Eq. 2 encodes a regressing Sun, in agreement with the project's #891 fix (`bcr4bp.py`: theta_S = theta_S0 - omega_S t), and with the z coefficient
-(1/2) m^2. This is a check the paper's printed sign passes; it can be made a test (section 5, item T1).

### 2.4 RHR4BP and the gamma homotopy (pp.788-789, Eqs. 6-8) READ
The Sun's effect on the spacecraft splits into the net perturbation (direct plus indirect, Scheeres Eq. 13 terms 1 and 2 of Eq. 6 of this paper) and the perturbation through the primaries (the pulsating Earth-Moon motion).
The RHR4BP removes the net Sun terms on the spacecraft and keeps the Sun's effect on the Earth-Moon motion (the variation orbit). It is "an incoherent three-body model", "not expected to be an accurate representation
of behavior; rather, the pulsation terms are simply introduced first." Common formulation (Eqs. 7, 8): same equations as Eq. 1 with V replaced by W,
  W = (1/2)(1 + 2m + m^2)(x^2+y^2) + gamma (m^2/4)(x^2+y^2) - (1/2) m^2 z^2 + (3/4) m^2 [(x^2 - y^2) cos 2 tau - 2 x y sin 2 tau] + (m^2/a0^3) [(1-nu)/R_{1-nu} + nu/R_nu].
This reading of the garbled text layer (p.789, Eq. 8) passes one check (COMPUTED): at gamma = 1 the (x^2+y^2) coefficient is (1/2)(1 + 2m + m^2) + m^2/4 = (1/2)(1 + 2m + (3/2) m^2), identical to Eq. 2, so gamma = 1 recovers the HR4BP.
It is NOT fully certain: the z^2 term and the cos 2 tau / sin 2 tau tide are printed unscaled by gamma, so in the RHR4BP as printed the tide stays on the spacecraft while only the m^2/4 (x^2+y^2) piece is removed; the paper's prose says all net Sun terms are removed. A page-image read of p.789 settles it before coding Eq. 8.
gamma = 1 is the Sun-Earth-Moon HR4BP, gamma = 0 the RHR4BP, 0 < gamma < 1 a homotopy "that evolves the model" (p.789).
m = 0 recovers the CR3BP (p.792). Figure 10 (p.795) is the map of the three continuations: CR3BP (m = 0) to RHR4BP (gamma = 0, m nonzero) and to HR4BP (gamma = 1, m nonzero) by continuation in m, and RHR4BP to HR4BP by continuation in gamma.

### 2.5 HFEM (p.789) READ
Spacecraft relative to a central body Pq (the Moon, because the southern L2 halos are studied) in J2000 inertial coordinates, with Earth, Moon and Sun as the perturbing bodies, ephemeris from SPICE with DE430.
Equation: r_qi'' = -G (m_i + m_q) r_qi / r_qi^3 + G sum_{j not equal to i, q} m_j ( r_ij / r_ij^3 - r_qj / r_qj^3 ) (the printed form; "G tilde" is the dimensional constant).

### 2.6 Comparison with the project's bicircular and quasi-bicircular models
| Item | HR4BP (this paper, Scheeres 1998) | BCR4BP (`core/bcr4bp.py`) | QBCP (`core/qbcp.py`, Andreu 1998) |
|---|---|---|---|
| Earth-Moon motion | coherent Hill variation orbit, Fourier in tau to order m^9 (pulsation: separation and angle) | circular, rigid | quasi-circular, eight Fourier alpha_i(theta_S) |
| Sun | tidal (quadrupole) field only, Sun at infinity, in the Hill limit; Sun direction -tau | full point mass at a_S = 388.81 EM units, direct plus indirect | coherent, same a_S scale |
| Frame | M1-M2 rotating, origin at B1 | Earth-Moon synodic, Earth at (-mu,0,0) | Earth-Moon synodic, pulsating |
| Time unit | 1/(n - n'), 2 pi = one synodic month | 1/n, Sun period 2 pi/omega_S = 6.79 TU = 29.53 d | 1/n |
| Forcing period | 2 pi in tau (the variation orbit and the tide both repeat each synodic month, tide at pi) | 2 pi/omega_S | 2 pi/omega_S |
| Parameters | nu about 0.01215, m about 0.0808 | mu 0.0121505816, mu_S 328900.54, a_S 388.811, omega_S 0.925196 | as BCR4BP plus alpha tables |
| Periodic orbits exist for | T = 2 pi k in tau | T = 2 pi k/omega_S | same |
The Hill model is the large-Sun-distance limit of the bicircular and coherent models, so it is expected to differ from them at order (a_S)^-1 relative terms; the paper does not quantify this and neither do I.
Key difference for the project: in the HR4BP the Moon moves on a pulsating, eccentric (Hill variation) orbit; in the BCR4BP it does not. The paper's RHR4BP/ER3BP-sidereal comparison (pp.795, 803) concerns exactly that
pulsation. The paper also says (p.803) its P2 RHR4BP counterpart "resembles the BCR4BP counterpart [Boudad et al. 2020]" and the closely bound HFEM geometry (Boudad et al. 2022).

## 3. The ephemeris-transition method (pp.789-791)
READ. "Direct" transition (stacking), as used for the figures:
1. Take a CR3BP halo orbit of resonance p:q. Discretise into patch points along the orbit (the number is "adequate to represent the periodic orbit"; not printed).
2. Stack it repeatedly until the desired interval is reached: the 9:2 orbit has period about 6.55 days, nine stacked orbits give about 59 days (p.789).
3. Choose an epoch. Shift the patch points to the CR3BP Moon-centred rotating-pulsating frame, then to J2000.
4. Differential corrections for continuity in the Sun-Earth-Moon HFEM (a multiple-shooting continuity corrector; the algorithm details and tolerances are not printed).
5. Shift the corrected states back to the Moon-centred rotating-pulsating frame for display.
The rotating-pulsating frame: Moon-centred, variable angular velocity so that the Earth stays at a fixed nondimensional distance of unity although the real Earth-Moon distance varies (p.789).
The initial state for stacking is "arbitrarily" the CR3BP apolune state (p.790).
Findings, with epochs (Figs. 1-4, p.790-791): 9:2 (NRHO region), 2:1 (regular halo) and 3:1 (interface) orbits transitioned over 3 months and over 4 years.
Reference epochs printed in the figure titles: 2024 JUN 02 07:16 (Figs. 1, 2), 2024 JUN 14 13:35 (Fig. 3), 2024 JUN 01 00:00 (Fig. 4). Results: (i) the direct transition for the 3:1 orbit FAILED TO CONVERGE over 4 years (Fig. 2 omits it), "consistently
observed in different orbits from the interface region"; (ii) the 3:1 geometry in the HFEM is strongly epoch-dependent (compare Figs. 1c and 3c) whereas 2:1 and 4:1 are not; (iii) the 9:2 and 2:1 HFEM orbits "seemingly translate their
geometries from the CR3BP over longer time horizons" whereas the 3:1 does not. Park & Howell (AAS 23-118) used a homotopy (not stacking) transition and see the same three-region split, so the failure is not only the algorithm (p.791).
Authors' stance (p.791): case-by-case numerical adjustment "is time-consuming, delivers point designs, and is not based in any dynamical foundations"; the intermediate model is to inform the numerical process.

## 4. Periodic orbits in the (R)HR4BP (pp.791-795)
READ. In a non-autonomous model periodic orbits are isolated with period commensurate with the model's: here multiples of 2 pi (p.792). A CR3BP orbit of minimal period T_c = 2 pi (q/p) (with q/p rational) gives an HR4BP orbit of
period 2 pi q, with p revolutions of the spacecraft in q revolutions of the primaries (sidereal-to-synodic resonance ratio p:q, p/q on the axes of Fig. 5b). (The text prints "2 q pi / p" for T_c and "2 q pi" for the HR4BP orbit; this is the p-revolutions-in-q-periods reading and is internally consistent.)
Phase parameter phi (Eqs. 9, 10): the variation orbit with argument 2n(tau + phi); "phi adjusts the initial location of the Earth and the Moon in their respective orbits" and does NOT change the Sun's net acceleration terms. A nonzero tau0 shifts ALL time-dependent terms, a nonzero phi only the variation orbit.
Perpendicular-crossing symmetry: orbits are symmetric across the xz plane; the states at tau0 and half a period are (x,0,z,0,y',0) (p.792). Four initial-condition classes (Table 1, p.792):
| Condition | tau0 | phi |
|---|---|---|
| C1 | 0 | 0 |
| C2 | pi/2 | 0 |
| C3 | 0 | pi/2 |
| C4 | pi/2 | pi/2 |
C3 and C4 break coherence (Sun, Earth and Moon then vary independently), so only C1 and C2 are used. Two families per CR3BP orbit, P1 and P2 (Table 2, p.794): odd p, P1 from C1 and P2 from C2, either apolune or perilune start (same orbit, different phase);
even p, P1 apolune at C1 and C2 (identical orbit, differing in phase), P2 perilune at C1 and C2. The paper's convention: apolune start for odd p.
Continuation in m (Section 4.1): at fixed gamma, from m = 0 (CR3BP) to m about 0.0808; the orbit period changes with m (see 5:2 below). Continuation in gamma (Section 4.2) at fixed m = 0.0808 from HR4BP to RHR4BP uses the same perpendicular-crossing conditions.
Natural-parameter continuation with step-size reduction to find the "linked" branch (p.799).

## 5. Printed numbers (all of them; none is an orbit state)
The paper prints NO orbit initial-condition table, NO period table, NO Jacobi-constant table and NO multiplier table. Everything below is either a model constant or a number in running text, captions or axes.
| ID | Item | Value | Where | Status |
|---|---|---|---|---|
| N1 | nu = M2/(M1+M2), Earth-Moon | about 0.01215 | p.788 | text |
| N2 | m = (n'/n)/(1 - n'/n) | about 0.0808 | p.788, Figs. 6, 7, 13, 14, 17 | text |
| N3 | a_n(m) order of evaluation | m^9 | p.788 | text |
| N4 | gamma range and meaning | 0 to 1; 0 = RHR4BP, 1 = HR4BP | p.789 | text |
| N5 | forcing period | 2 pi in nondimensional time, orbits periodic at 2 pi k | pp.788, 792 | text |
| N6 | 9:2 synodic resonant halo period | about 6.55 d | p.789 | text |
| N7 | 9 stacked 9:2 orbits | about 59 d | p.789 | text |
| N8 | 5:2 stacked ~5 revs at m = 0 | 2 x 27.32 d (two sidereal months) | p.795 | text |
| N9 | 5:2 stacked at m about 0.0808 | 2 x 29.5 d (two synodic months) | p.795 | text |
| N10 | 5:2 period at general m | 2 x 27.32 x (1 + m) d | p.795 | text |
| N11 | nondimensional 5:2 period at m = 0 / at m | 4 pi/5 / 4 pi (1+m)/5 | p.795 | text |
| N12 | 5:2 broken bifurcation | m about -0.00671; orbit of period about 4 pi (1 - 0.006715)/5 | p.795 | text |
| N13 | multiplier of that orbit | about cos(4 pi/5) + i sin(4 pi/5) (period-5 multiplying bifurcation) | p.795-796 | text |
| N14 | 7:3 broken bifurcation | m about 0.0787; orbit period about 6 pi (1 + 0.0787)/7 | p.796 | text |
| N15 | multiplier at 7:3 | cos(2 n pi/7) + i sin(2 n pi/7), n = 1..7 (period-7 multiplying) | p.796 | text |
| N16 | general rule | period-p multiplying bifurcation at multiplier exp(2 i n pi/p), n = 1..p; underlying halo period in [2 pi q/p, 2 pi q (1+m)/p] as m runs 0 to m | p.796 | text |
| N17 | P5HO5 family members used | period 4 pi (nondim); hodograph Jacobi constant axes about 2.99 to 3.015, period axes 11.5-14.5 | Fig. 12, p.796 | axis only |
| N18 | resonant orbits studied | 9:2, 2:1, 3:1, 4:1 (Figs.1-4); 5:2, 7:3 (Figs.11-13); 12:5, 13:5, 14:5, 20:7, 13:4, 11:3 (Figs.14-17) | text | text |
| N19 | gamma values shown in homotopy figures | 13:4: 0.9499, 0.6517, 0.5412, 0.3383, 0.0646, 0.0071 (Fig. 19); 13:5: 1, 0.8, 0.65, 0.5, 0.2709, 0.2342 (Fig. 20); 12:5: 0.8596, 0.5088, 0.2406, 0.1003 (Fig. 21) | p.801-803 | captions |
| N20 | HFEM epochs | 2024 JUN 02 07:16; 2024 JUN 14 13:35; 2024 JUN 01 00:00 | Figs. 1-4 | captions |
| N21 | HFEM ephemeris | DE430, Earth+Moon+Sun, Moon-centred | p.789 | text |
| N22 | classification of regions | 9:2 NRHO region, 3:1 interface region, 2:1 and 12:5 regular halo; 11:3 NRHO-region behaviour (families intact in gamma); 13:5, 14:5, 20:7, 13:4 interface (inconsistent HR4BP vs RHR4BP geometry) | pp.790, 799-802 | text |
| N23 | perilune-radius axis of the halo hodograph | 0 to about 5 x 10^4 km, p/q from 1.5 to 5 | Fig. 5b | axis |
Internal checks (COMPUTED): 9 x 6.55 = 58.95 d, consistent with 59 d (N6, N7); 2 x 29.53 / 9 = 6.56 d, consistent with N6 as a synodic 9:2 orbit; the project's omega_S = 0.925195985520347 gives n'/n = 1 - omega_S = 0.074804 and
m = 0.074804/0.925196 = 0.080852, matching N2 to its printed 3 s.f.; the Moon's synodic/sidereal ratio 29.53059/27.321661 = 1.08085 gives the same m = 0.08085, confirming N10 and N8-N9; 4 pi (1 - 0.006715)/5 = 2.4964 nondimensional, 6 pi (1.0787)/7 = 2.9047.

### Suitable as sourced tests (T-list)
- T1: the Sun sense in the Hill tide: the printed coefficients of Eq. 2, (3/4) m^2 [(x^2-y^2) cos 2 tau - 2 x y sin 2 tau], equal the quadrupole of a Sun at angle theta0 - tau (regressing) and differ in sign of the cross term from a Sun at theta0 + tau. Independent derivation from the project's own `bcr4bp.py` sun_position (a_S to infinity limit). A non-rotating-frame identity test of the kind the project asks for.
- T2: m from the project's omega_S equals the printed 0.0808 to 3 s.f.; nu from the project's mu equals 0.01215 to 4 s.f. (a constants reconciliation, not a golden).
- T3: 5:2 stack period relation, 2 x 27.32 x (1 + m) days with m = 0.0808 equals 2 synodic months (59.1 d against 2 x 29.53 = 59.06 d). A unit/time-scaling test.
- T4: the Floquet rule N16: at the CR3BP halo whose monodromy has an eigenvalue exp(2 i n pi/p), a period-p multiplying bifurcation sits; a CR3BP test on the known 9:2 NRHO or on the P5HO5 bifurcation orbit (period about 4 pi (1-0.006715)/5 reported here is for the continuation seed, so it is an approximate target, not a golden to 1e-6).
- T5 (qualitative, structural): in the Hill model, orbits continued from the CR3BP exist only at periods T = 2 pi k in tau, equivalently in the BCR4BP at T = 2 pi k/omega_S; assert the project rejects non-commensurate period requests.
NOT test-worthy: any ratio read from figures; the figure axes ranges; the broken-bifurcation m values to more than the printed digits (-0.00671 and -0.006715 are both printed for the same event, p.795).

### Errata and slips in the paper (so nobody trusts them)
- Frame rate "(1 + 1/m)" in prose against (1+m) in Eq. 1 (p.787).
- Section 6 text (p.799) refers to "the 13:5 case (Fig. 16(e))" but Fig. 16(e) is 13:4 per the caption; 13:5 is panel (b); the text then correctly says Fig. 15(b) and 16(b).
- Fig. 16 has only 5 of the 6 listed panels visible in the text layer (panel (d) 20:7 and (f) 11:3 absent from the extracted axes); Fig. 17 has five panels for six ratios (13:5 is explained in text as not reaching m = 0.0808).
- Axis typos in the extraction (1.4 for 1.14 on Fig. 13 and 15; 1.2 for 1.12 on Fig. 6): extraction or print slips, not used.
- "Paramter m" on the figure axes. Fig. 12 caption says "(c) Orbits with period = 4 pi n.d.".

## 6. Reconciliation against project code
Files read: `src/cyclerfinder/core/bcr4bp.py`, `core/qbcp.py`; the CCR4BP modules (`core/ccr4bp*.py`) are moon-pair circular models and have no relation to this paper. No Hill-model code and no Earth-Moon ephemeris-transition code exists under `src/cyclerfinder`: grep for "hill" finds only Hill-sphere radii (wsb, V4 lanes), never the Hill four-body potential. The only circular-to-ephemeris continuation code is `search/continuation.py` (Russell 2004 section 5.4: the PLANET MODEL is the homotopy parameter, heliocentric Earth-Mars), plus `nbody/shooter.py` and `search/multiarc_cycler.py`.
| Item | Paper | Project | Verdict |
|---|---|---|---|
| Sun sense in Earth-Moon rotating frame | Sun at angle theta0 - tau (derived from Eq. 2, section 2.3) | `bcr4bp.py`: theta_S = theta_S0 - omega_S t; QBCP `evaluate_alphas` uses theta = theta_sun0 + omega t | BCR4BP agrees. The QBCP `evaluate_alphas` advances its argument as theta_sun0 + omega t; that is the argument of Andreu's alpha_i Fourier series, whose relation to the Sun's direction in the frame I did not re-derive. Check it once against the #891/#892 non-rotating-frame identity before calling it consistent |
| Mass ratio | nu about 0.01215 | mu = 0.012150581600000 (BCR4BP), 0.012150581623433623 (QBCP) | agree to the printed digits; the two project constants differ at the 1e-11 level, already a known project item |
| Sun parameter | m = (n'/n)/(1-n'/n) about 0.0808 | omega_S = 0.925195985520347 so m = 0.080852 | agree (COMPUTED) |
| Time unit | t* = (1+m)/n = 1/(n - n') (synodic) | 1/n (sidereal) | DIFFERENT. One Hill unit = (1+m) = 1.0808 BCR4BP units (COMPUTED: 29.53 d / 2 pi = 4.70 d against 4.348 d for the EM unit). Lengths are also scaled differently (the HR4BP length unit carries a0 and the Hill factor), so coordinates must not be compared without a conversion the paper does not print. |
| Orbit period labels | 2 pi k in tau (k synodic months), p:q = p revs in q synodic months | `sun_commensurate_period(omega_sun, n) = 2 pi n/omega_sun` | agree; the n Sun revolutions of the code are the q synodic months of the paper. Boudad et al.'s p:q is the same count |
| Sun distance | Sun at infinity (Hill limit), tide only | a_S = 388.811 EM units, direct plus indirect | the HR4BP is a limit, not the same model; do not expect identical orbits |
| Lunar motion | Hill variation orbit, pulsating | circular (BCR4BP) | HR4BP carries the pulsation; BCR4BP does not. QBCP carries a Fourier version of it (alphas) |
| Symmetry classes | tau0 in {0, pi/2} and phi in {0, pi/2}, only C1, C2 coherent | not stated in the code | the BCR4BP symmetry (x, y, t) to (x, -y, -t) holds only for Sun phase theta0 = 0 or pi (the Sun term is odd in time otherwise). The Hill tide has period pi in tau so theta0 = pi/2 is also a symmetric class there. INFERRED, not tested. |
No numerical conflict with any constant in the project was found. The qualitative reconciliations are in section 7.

## 7. Techniques applicable to the project's problems
Honest scope: this paper supplies diagnostics and a continuation strategy, not data. Nothing here can be a V-tier anchor, and the task mapping of #890/#895 below is corrected (see the note under #890/#895).
### #388 (circular-to-ephemeris continuation wall)
1. The paper's central lesson transfers: the failure of a continuation can be STRUCTURAL (a broken bifurcation joining the parent to a nearby higher-period family, destroyed under perturbation) rather than algorithmic. #388 has never tested this. Concrete step: for the converged circular-coplanar parent, compute the multipliers of the linearised return map over k synodic repeats and flag any multiplier near exp(2 i n pi/p) for the repeat count p used. The paper's rule (N16): the bifurcation parameter value is predicted by that eigenvalue condition. If the #388 stalls occur near predicted values of the ramp parameter (eccentricity, inclination), that is a diagnosis, not a tuning problem. (a) 40 percent that a root-of-unity multiplier coincides with a stall; the test is cheap, 2 to 4 agent-hours, and a clean negative is informative. For a heliocentric cycler with flybys the "monodromy" is a composition of Lambert/flyby maps, not a flow, so the multipliers need the analytic STM Jacobian already built for the shooter (see the project memory on the #388 saga).
2. Split the homotopy into physical causes, as the RHR4BP does. The Russell ramp (`search/continuation.py`) moves eccentricity, inclination and ephemeris together in a fixed order. The paper's two-axis picture (pulsation, i.e. eccentricity of the primaries, against net third-body gravity) maps to: eccentricity ramp (planet orbit shape and the pulsation of the synodic geometry) separate from Jupiter/other-planet perturbation and from the planet phase (epoch). Run each alone and in the plane (e, gamma-like Jupiter strength) as a 2-D continuation map, and record where branches disappear. This attributes the wall to one effect, which a single ramp cannot.
3. Epoch dependence: the paper shows the HFEM geometry depends strongly on the epoch only in the interface region. For #388 the corresponding experiment is the same cycler seed over the 21 launch windows already in Russell & Ocampo's procedure, reporting converge/no-converge per window. Cheap positive evidence of the same structure.
4. Direct stacking versus homotopy: the paper's own "direct" transition (stack, shift, correct) fails where a homotopy succeeds (Park & Howell). #388's lane is already a homotopy; the useful corollary is not to spend more effort on direct re-solves.
### #884 / #905 (Sun-forced Earth-Moon search rerun)
1. Closure rule: periodic orbits in a Sun-forced model have T = 2 pi k in tau, equivalently T = 2 pi k/omega_S in the project's BCR4BP. This is the same check #905 already demands ("closure is over p Sun periods"), now sourced to the paper (N5, T5).
2. Symmetry classes: the paper's C1/C2 and P1/P2 table (Tables 1-2) is a ready scheme for #905's "one orbit per symmetry class": the apse start (apolune/perilune) and the Sun phase class give distinct orbits for even p and the same orbit for odd p. In the BCR4BP, search only Sun phases 0 and pi (INFERRED, section 6).
3. Bifurcation detection in the family walk (the #905 review's section 14 item): the paper gives a precise recipe: monitor eigenvalues for exp(2 i n pi/p) along the walk; to find the branch, reduce the step until a continuous curve is seen; to find a broken branch, use a larger step to land on it, then continue. Adopt as the continuation policy, because a single-step-size natural-parameter walk silently jumps branches.
4. Which family to follow at a broken pitchfork: "incoming branch and three outgoing" (p.799). Record the branch identity (not just the family label) in every row, since the geometries and stability change between branches.
5. Cross-model check: run the HR4BP (Brown 2024 and this paper define it) as a third Sun model; agreement of the resonant orbit set across BCR4BP, QBCP and HR4BP is a model-independence check, and the interface-region behaviour (HR4BP versus RHR4BP inconsistency) predicts WHICH Earth-Moon resonances should be fragile. Not needed before the rerun.
### #890 / #895 (note on the task description)
In `data/OUTSTANDING.md`, #890 is the periodic Titania-Oberon-Titania orbit in the planar circular four-body model and #895 the real-ephemeris arcs for it; they are Uranian, not Earth-Moon. The Hill paper does not apply physically (no Sun-like third body). The method transfers in three ways:
1. The #890 failure ("the patched-conic closure misses the Oberon flyby by 967,000 km", energy change from the flyby offset accumulating over 5.5 revolutions) is an orbit isolated in a non-autonomous (commensurate) model that does not survive a perturbation. The homotopy-with-diagnosis scheme applies: introduce a perturbation-strength parameter (Titania and Oberon masses, the moons' eccentricities, Uranus J2) from zero to physical, continue the orbit, and watch multipliers for roots of unity. A branch loss at a root-of-unity multiplier is then a structural explanation to put beside the existing ones.
2. The moons' orbital eccentricity/pulsation is the exact analogue of the Earth-Moon pulsation: the paper's separation of "pulsation" and "net perturbation" suggests ramping eccentricity of the two moon orbits separately from the third moon's gravity.
3. The corrector: the direct stacking-and-correct step (patch points along the orbit, shifted to the inertial frame, corrected for continuity) is the same pattern as the #895 corrector; no new capability.
### #902 (bicircular tiers)
No help for anchoring: the paper prints no orbit state, period or multiplier, so it cannot supply the "real orbit" the rebuilt tiers need. What it adds: (i) the sourced statement that the P2 RHR4BP orbit resembles the BCR4BP counterpart (p.803) which points to Boudad et al. 2020's NRHO set as the better anchor; (ii) the Floquet-root-of-unity rule as a check that a tier orbit is not sitting on a bifurcation (a tier orbit whose multiplier is within the corrector's tolerance of exp(2 i n pi/p) is not isolated and its V3 status is fragile).
### #913 to #918
- #913, #914, #915, #918: heliocentric Russell/Strange solutions and moon cyclers; nothing in this paper applies beyond the generic remark in #388 item 4.
- #916 (printed persistence conjecture, Ross and Scheeres): the paper is direct evidence on a related point: a periodic family that exists in an unperturbed model (CR3BP) is destroyed in a model with a particular perturbation (RHR4BP, gamma to 0) in the interface region but persists elsewhere; persistence depends on the perturbation, the resonance and on nearby higher-period families. Cite as evidence in the #916 test design; not a test of the conjecture.
- #917 (Casoliva rows in the elliptic problem with the fold) and #922 (the #890 orbit as an invariant two-torus): the fold and the broken-bifurcation diagnostics here (hodograph x0 against m, with gaps invisible except at small step) are the same objects; use the step-size policy of section 7 (#884 item 3).

## 8. Recommended follow-ups (not registered; for the owner to number)
1. Constants test: add T1 (Sun-sense identity for the Hill tide), T2 (m and nu reconciliation) and T3 (5:2 period relation) as sourced tests; no code change needed beyond a test module. Cost under 1 agent-hour.
2. Multiplier root-of-unity monitor: a function that, along a family walk or a continuation, reports multipliers within a tolerance of exp(2 i n pi/p), and wire it into the #905 BCR4BP family walk and the #388 ramp. Run first on the CR3BP L2 halo family: the P5HO5 seed near period 4 pi (1 - 0.006715)/5 is the paper's example, and the Boudad NRHO family a second check. This is a positive control (a known bifurcation is found) before any claim.
3. #388 attribution experiment: the 2-D continuation map (eccentricity against other-planet perturbation) and the per-launch-window converge/no-converge table (section 7, #388 items 2-3). Use the existing Russell continuation driver; report where branches vanish.
4. A Hill four-body module (Eqs. 1-4, 7-8) is NOT recommended as a build now: the Brown et al. 2024 digest already carries a fuller account and the project's Sun models are the bicircular and coherent ones. If a third Sun model is wanted for #884/#905, build it then, from a page-image read of Eq. 8 (the text layer is garbled, section 2.4) and the coefficients a_n(m) from Olikara & Scheeres 2017 (not held; check the corpus index first).
5. Read p.789 as an image and the following figures (Figs. 8, 9, 18, 19 hodographs) if a quantitative HR4BP control is wanted; the text layer carries no hodograph numbers, so no x0 values are available.
6. Obtain Davis et al. (AAS 17-826), Park & Howell (AAS 23-118, AAS 22-741), Sanaga & Howell (AAS 23-227) and Olikara & Scheeres (AAS 17) if held-paper work on the Earth-Moon halo interface region proceeds; the printed ICs for the 3:1 orbit, if any, would be in those, not here.
7. Cross-check the QBCP Sun-angle direction once against the #891/#892 non-rotating-frame identity (section 6, first row), so the "theta0 + omega t" in `qbcp.evaluate_alphas` is documented as Andreu's convention rather than left looking like the #891 bug.

## Note added 2026-10-04: what Scheeres 1998 settles (digest `2026-10-04-digest-scheeres-1998-restricted-hill-four-body-problem.md`)

1. N3 (a_n(m) to order m^9, attributed here to Olikara & Scheeres 2017): Scheeres 1998 itself gives the a_n/a0 only to order m^6 (and a0 to m^3) from Wintner; the m^9 values are not in it. The Brown et al. 2024 tables (order m^9, in the Brown digest) re-expanded in m reproduce Scheeres' appendix exactly through m^6 (symbolic check), so the m^9 tables are consistent with the 1998 origin. The Olikara & Scheeres 2017 attribution is not contradicted but is not needed to obtain the coefficients.
2. The model equations here (Eqs. 2 and 3 as read in section 2) agree with Scheeres eqs. 55 and 56 (potential, Coriolis 2(1 + m), period-pi tide); the Sun direction regresses in the Earth-Moon rotating frame, confirming T1 from the source model (a points along angle -tau in the Moon-fixed frame, Scheeres eqs. 47 and 51).
3. T3 and N-values: Scheeres prints m = 0.0808 and nu = 0.0122 (approximately), consistent with N2.
4. Item 4 of the recommendations (a Hill module would need Eq. 8 read from the page image and the a_n coefficients): the coefficients are now available (Scheeres 1998 appendix and the Brown tables). Scheeres' Eq. 61 (expansion in m) was verified here against his Eq. 56 numerically.
