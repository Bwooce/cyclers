# Digest — Brown, Peterson, Henry & Scheeres (2024), "Structure of Periodic Orbit Families in the Hill Restricted 4-Body Problem"

**Digested:** 2026-10-03 (text-layer PDF, no OCR needed; read in full, 27 pages).

**Citation (preprint):** Gavin M. Brown, Luke T. Peterson, Damennick B. Henry, Daniel J. Scheeres
(Smead Aerospace Engineering Sciences, University of Colorado Boulder), "Structure of Periodic
Orbit Families in the Hill Restricted 4-Body Problem", **arXiv:2402.19181v1 [math.DS], 29 Feb 2024**.
**The file we hold is the arXiv v1 preprint, NOT the published version.** CrossRef gives the
published version as SIAM Journal on Applied Dynamical Systems **24(1)**, pp. 346-375 (2025),
**DOI 10.1137/24M1637301**. Page and section numbers below refer to the 27-page preprint; the
published version may differ in content and numbering (not checked). Funding acknowledged: U.S. Air
Force Office of Scientific Research grant FA9550-21-1-0332.

**Filed in the private paper corpus as**
`brown-peterson-henry-scheeres-2024-periodic-orbit-families-hill-restricted-4-body-problem-siads-24-346-arxiv-2402.19181.pdf`
(md5 `4b55047e198ddbfdb517e5602b477517`, 27 pages).

**System:** Sun-Earth-Moon (SEM), Earth-Moon (EM) libration-point orbits, in the Hill restricted
4-body problem (HR4BP). **Acquired for:** what is published about how CR3BP periodic-orbit families
behave when a fourth body's periodic forcing is added (relevant to our BCR4BP / QBCP / CCR4BP
periodic-orbit and torus work).

## What it is
A computational paper on periodic-orbit family structure in a coherent, time-periodic Sun-Earth-Moon
model. Abstract: "The Hill Restricted 4-Body Problem (HR4BP) is a coherent time-periodic model that
can be used to represent motion in the Sun-Earth-Moon (SEM) system. Periodic orbits were computed in
this model to better understand the periodic orbit family structures that exist in these types of
systems. First, periodic orbits in the Circular Restricted 3-Body Problem (CR3BP) representation of
the Earth-Moon (EM) system were identified. A Melnikov-type function was used to identify a set of
candidate points on the EM CR3BP periodic orbits to start a continuation algorithm. A pseudo-arclength
continuation scheme was then used to obtain the corresponding periodic orbit families in the HR4BP
when including the effect of the Sun. Bifurcation points were identified in the computed families to
obtain additional orbit families." The results are presented graphically (hodographs, bifurcation
diagrams); there are no tables of orbit initial conditions or periods (see below).

## The model and its constants (Section 2.1, Appendix A)
**Relation to the CR3BP and the bicircular problem (Section 1).** The paper lists the Elliptic
Restricted 3-Body Problem (ER3BP), the Bicircular Restricted 4-Body Problem (BCP), the Quasi-Bicircular
Model (QBCP) and the HR4BP as increasingly realistic models of cislunar space. It states the BCP "is
incoherent as it does not account for the effect of the Sun on the Earth or the Moon" and that "there
is no accurate dynamical equivalent to L2 which is one major drawback of the BCP [9]". The QBCP
"accounts for the effect of the Sun on the Earth and Moon by modeling their motion as a solution to the
3-body problem [14]". The HR4BP, "developed by Scheeres [16] in 1998, is another coherent time periodic
model describing the motion of a small body (P3) in the presence of three large bodies (P0, P1, and P2).
The model is a higher fidelity model than the CR3BP and BCP, more accurately represents the true
dynamics in the Sun-Earth-Moon (SEM) system, is easier to implement than the QBCP, and has previously
been used to study the Sun's effect on the EM system [16-18]."

**Reduction to the CR3BP (Section 2.1):** "It is important to note that when m -> 0, the HR4BP equations
of motion take on the form of the CR3BP equations of motion. As m increases from zero, the effect of
the more massive body on the system becomes more pronounced, and m_SEM = 0.0808 for the SEM system."

**Bodies and frame.** P0 largest body (Sun); P1, P2 the primaries (Earth, Moon); P3 negligible mass.
B is a rotating frame with constant angular velocity (direction k-hat) with origin at the centre of mass
of P1, P2, P3; P0, P1, P2 lie in the xy-plane. The distance between the Sun and the Earth-Moon barycentre
is taken as d_a, with R_c ~ a = d_a i-hat.

**Parameters (Eq. 1):** mu = M2/(M1+M2); m = (n0/n)/(1 - n0/n); nu = (M1+M2)/M0. n0 is the mean motion of
P0 treating the two primaries as one body at R_c; n is the mean motion of the primaries; "m represents
the relationship between the period of Rc about the total system center of mass and the period of P1 and
P2 about Rc". The paper notes mu and nu are defined differently from Scheeres [16] to match standard
CR3BP notation.

**Constants as printed:** "In the SEM system, m is the ratio of the difference between the synodic and
sidereal months relative to a sidereal month (i.e., m ~ 0.0808), and mu ~ 0.0122." The Conclusion repeats
"m = 0.0808 and mu = 0.0122". Numerical value of nu: not stated. Numerical value of d_a: not stated.

**Scaling (Eq. 2):** 1 MU = M0; 1 DU = a0 d_a nu^(1/3); 1 TU = m/n0 = (1+m)/n. "Time is scaled such that
the relative physical configuration of P0, P1, and P2 repeats every 2 pi time units (i.e., 2 pi time
units is equal to a synodic month). However, the equations of motion are pi-periodic [16]." 1 DU is "the
average distance between P1 and P2". The forcing period is therefore Tg = pi.

**Potential (Eq. 4):** V = (1/2)(1 + 2m + (3/2)m^2)(x^2 + y^2) - (1/2)m^2 z^2 + (3/4)m^2 (x^2 - y^2 cos 2tau
... ) as printed, plus (m^2/a0^3)((1-mu)/R_{1-mu} + mu/R_mu); the full printed form is Eq. 4 on p. 4 (the
text layer scrambles the x^2 - y^2 cos 2 tau - 2xy sin 2 tau term; consult the PDF). Equations of
motion: X' = [r'; -2(1+m) Omega x r' + grad V] (Eq. 5a) with variational equations for the state
transition matrix [Phi] and the parameter sensitivity [Psi] = [dX/dm, dX/dmu] (Eq. 5b, 5c). Full gradient,
Hessian and sensitivity expressions are in Appendix A (Eqs. 17-20).

**Primaries' motion.** The primaries' positions are obtained from the Hill variation orbit (HVO), a
Fourier-series solution of Hill's equations, Eq. 16, with coefficients d_p (Table 1) and c_{n,p} (Table 2)
computed to maximum order P = 9, N = 4 (see the tables below).

**Symmetries (Eq. 6):** S1: (x, y, z, x', y', z', k pi + tau) -> (x, -y, z, -x', y', -z', k pi - tau);
S2: the same with k pi + pi/2. "if the initial condition X(tau = tau0) = [x0, y0, 0, x0', y0', 0]^T
corresponds to a periodic orbit with a period that is an integer multiple of pi, then X(tau = -tau0) =
[x0, -y0, 0, -x0', y0', 0]^T also corresponds to a periodic orbit with the same period."

## Central results (with quotes and section numbers)
The paper is a methods-plus-catalogue paper; it states no theorem about which CR3BP families "survive" or
"break up", and the word "isola" does not appear. What it states:

1. **Resonance with the forcing is a necessary condition (Section 2.2).** "Periodic orbits in this system
   must have a minimal period (T) that is in some resonance with the forcing period (Tg) where Tg = pi for
   the HR4BP (i.e., T = b Tg where b is a positive integer) [35]." A CR3BP orbit of period T* is a
   candidate if "b Tg = a T* where a and b are relatively prime integers", and then the perturbed orbit has
   period T = b Tg.

2. **One CR3BP structure has many HR4BP "dynamical equivalents" (Section 2.2).** "When a time-periodic
   perturbation is added to an autonomous system, a single structure in the autonomous system may have
   multiple 'dynamical equivalents' in the perturbed system... For example, at least four dynamical
   equivalents to the 9:2 NRHO in the BCP have been identified previously [11]. In that system the
   continuation of these orbits followed many possible paths, and we expect a similarly complicated
   behavior in the HR4BP."

3. **Families exist only after fixing period and epoch (Section 3).** "As this system is non-autonomous,
   we do not expect to identify families of periodic orbits as we would in an autonomous system [36].
   However, we do expect to identify periodic orbit 'families' in the HR4BP if we fix the value of T and
   initial time (tau0), and perform a continuation while allowing the initial state (X0 = X(tau = tau0))
   and at least one of the parameters m and/or mu to vary." All results use tau0 = 0 unless stated.

4. **Which CR3BP points continue: zeros of a Melnikov-type function (Sections 2.2, 3.1).** "Provided
   M(s, tau0) is not identically zero for all s in [0, T*), we expect to be able to continue periodic
   orbits from the unperturbed system into the perturbed system at the points s on the unperturbed orbit
   provided the initial time when beginning the integration is tau0." The function used is that of
   Cenedese and Haller [34]: M(s, tau0) = integral over [0, aT*] of g(*X(s+tau), tau0+tau) . r'(s+tau)
   dtau (Eq. 8). The first-order term g1 contributes zero work (it reflects the time scaling), so the
   expansion is carried to h2 (Eq. 11, epsilon = m^2), and to h3, h4 if the function vanishes identically.

5. **Three propositions (Section 3.2, proofs in Appendix B).** Prop. 1: shifting the point along the orbit
   equals shifting tau0 the other way. Prop. 2: M at any point follows from M at two points
   pi/4 apart along the orbit, so "(11) must be integrated only twice"; if M(s, tau0) = 0 then M(s + k pi/2, tau0) = 0,
   and if both vanish "the Melnikov function is identically zero for any other value of s and tau0" on
   that orbit, forcing a higher-order h_j. Prop. 3: for orbits satisfying the half-period symmetry
   conditions (Eq. 13), M(s, tau0) = 2AB with A = sin(2 tau0) if a = 1 and A = 0 otherwise. Quote: "if a CR3BP
   periodic orbit has a period T* = k pi, k in Z+ (i.e., if a = 1), then any points on that orbit that
   satisfy the half-period symmetry conditions presented in (13) are points where A = 0 provided
   tau0 = k1 pi/2. For example, let us consider the CR3BP L1 and L2 planar Lyapunov and halo periodic
   orbits with T* = pi. All points on these orbits that lie on the xz-plane are points where we expect to
   be able to continue a corresponding HR4BP periodic orbit family with tau0 = 0. If T* is not an integer
   multiple of pi (i.e., if a > 1), but there is at least one point on the orbit that satisfies (13), then
   the Melnikov function is identically zero when using h2."

6. **Validation on L4 (Section 4, Fig. 2).** The CR3BP L4 planar orbit with T* = 2 pi (a = 1, b = 2) has
   four Melnikov zeros, "the starting points for four HR4BP orbit families"; these four were already
   found by Scheeres [16] and Peterson et al. [21]. "To validate our techniques, 100 points were selected
   on the initial CR3BP L4 planar 2 pi periodic orbit, and the continuation procedure was attempted
   starting at each of these points. Four of these points produced families in the HR4BP that could be
   continued up to values of m that were not negligible. These four points matched the four points where
   the Melnikov function is zero."

7. **Orbit versus object families (Section 4).** An "orbit" is states with their times tau; an "object"
   is the set of states irrespective of tau. "two points on the CR3BP L2 Northern Halo family member with
   T = pi can be continued (using Proposition 3) to produce two different HR4BP orbit families. At a
   particular value of m, the orbits in these two families contain the same states, with their
   corresponding values of tau shifted by pi/2. Therefore, only one object family in the HR4BP was
   identified corresponding to that particular CR3BP orbit." The continuation in Fig. 3 "starts with m
   increasing (from yellow to magenta), before decreasing (from magenta to dark blue), until the CR3BP L2
   Southern Halo orbit with T = pi is obtained", i.e. the northern halo family continues to the southern
   halo at m = 0 through a turning point in m.

8. **L1 vertical orbit (Section 4, Fig. 4).** The L1 vertical orbit with T* = pi has four points satisfying
   the symmetry conditions, separated by pi/4, so the Melnikov function is identically zero for h2; with
   h4 "four zeros corresponding to the four points on the xz-plane. Each of these points produced a
   different HR4BP orbit family which belonged to one of three different object families... the object
   family on the right side of Figure 4 consists of two different orbit families."

9. **Bifurcation structure (Sections 3.3, 4.1).** Families were continued "until either the maximum m
   value at which the orbit of the primaries is stable was reached (m = 0.19510486) [16], or the
   continuation algorithm failed to identify a new orbit member." Bifurcations are detected as local
   minima of the two smallest non-zero singular values of the modified corrections Jacobian; "many points
   corresponding to these local minima are where multiple families intersect" (Section 4.1). "Note that
   all three families associated with L2 identified in [18] (the A, B, and C families) were identified in
   this work", with "minute differences" near bifurcation points attributed to a different HVO
   representation. "Note new families were only computed at tangent bifurcation points in this study."
   Conclusion: "We expect to find many other families connected to the families identified in this paper
   by analyzing period-multiplying bifurcations. While additional study is needed to completely map out
   the periodic orbit structure in the HR4BP, this work presents techniques that can be used to generate
   resonant periodic orbits in periodically forced systems."

10. **SEM members (Conclusion).** "The set of orbits in these families with m = 0.0808 and mu = 0.0122 are
    periodic orbits in the HR4BP representation of the Sun-Earth-Moon system."

## Families and resonances computed (Section 4.1, Figs. 5-11)
All starting orbits are EM CR3BP orbits whose period is an integer multiple of pi (the forcing period).
Stability: the paper does not give Floquet multipliers or stability indices for any orbit; its only
stability statement is the limit m = 0.19510486 for the primaries' orbit. Periods of the HR4BP members are
those of the starting CR3BP orbit (T = pi, 2 pi or 4 pi), since families are continued at fixed T. No
numerical period, Jacobi constant or amplitude is printed for any individual orbit.

| Region | Starting CR3BP orbit | T | Colour in figure |
|---|---|---|---|
| L2 (Fig. 5, 6) | L2 point itself (the "A" family) | pi | maroon |
| L2 | planar Lyapunov member | 2 pi | dark red, dark orange |
| L2 | vertical family member | 2 pi | orange |
| L2 | northern (and southern) halo member | pi | gold, light green |
| L2 | northern (and southern) butterfly member | pi | green, aquamarine |
| L2 | 9:2 near-rectilinear halo orbit (NRHO) | 4 pi | teal, light blue |
| L1 (Fig. 7, 8) | L1 point | pi | maroon |
| L1 | planar Lyapunov member | pi | red, orange |
| L1 | planar Lyapunov member | 2 pi | gold, green |
| L1 | vertical family member | pi | aquamarine, light blue, blue |
| L3, L4, L5 (Fig. 9, 10) | the L3, L4, L5 points | pi | maroon |
| L4, L5 | planar Lyapunov member | pi | red, orange |
| L4, L5 | planar Lyapunov member | 2 pi | gold, light green, green |
| L3, L4, L5 | vertical family member | pi | aquamarine, light blue, blue |

Families found from bifurcations of these initial families are drawn "magenta through cyan" in each
figure; the paper does not enumerate or name them individually, beyond the L2 A, B, C families of Olikara
et al. [18]. Fig. 11 shows the HR4BP orbits near the EM libration points at the SEM value of m. The
Melnikov-zero starting points are shown for the L4 T* = 2 pi orbit (Fig. 2, four families), the L2 halo
T = pi (Fig. 3) and the L1 vertical T* = pi (Fig. 4, three object families). The northern butterfly is
noted to have a southern counterpart by symmetry (Section 3.3).

Resonance language: the paper treats only resonance of the orbit period with the forcing period (T = b pi;
the 9:2 NRHO is a synodic resonance with period 4 pi, as listed). No other m:n resonances are computed;
no resonant (DRO or lunar-resonant) families are computed.

## Methods (Sections 3, 3.3)
- Expand the HR4BP equations about m = 0 (epsilon g = g1 m + g2 m^2 + g3 m^3 + ...; Eqs. 9, 10) and
  evaluate the Melnikov-type function on CR3BP orbits to find starting points; propositions reduce this to
  two integrations per orbit, or none via Prop. 3 symmetry.
- Pseudo-arclength continuation starting at m = 0, with the state X0 and m varying at fixed T and tau0.
  Continuation runs from m = 0 through m_SEM = 0.0808 to the upper limit 0.19510486 or until failure.
- Bifurcation detection: SVD of a modified corrections Jacobian [DB] = [[Phi(T,0)]^nB - I, [Psi_m(nB T, 0)]]
  (Eq. 14), 6 x 7 so one singular value is zero; local minima of the smallest non-zero singular values
  sigma_alpha, sigma_beta mark candidate bifurcations. nB = 1 for tangent bifurcations, 2 for
  period-doubling, nB for other period-multiplying; "new families were only computed at tangent
  bifurcation points". Branch switching by V0B = [X0A; mA] + Delta s0 * delta V_sigma (Eq. 15), using both
  + and - directions when the null vector is aligned with the symmetries (Eq. 6).
- Symmetry (Eq. 6, 13) used to reduce the starting-point search and to generate mirror families (for
  example southern from northern butterfly).
- Integration scheme, tolerances, step sizes, correction convergence thresholds, software: not stated.
- The HVO coefficients (Tables 1, 2) follow Wintner [38] and are "modified slightly from the form
  presented by Olikara and Scheeres [22]".

## Tables transcribed (Appendix A)
These are the only numeric tables in the paper. Transcribed from the text layer of the PDF and re-checked
against the column layout; coefficients of the Hill variation orbit, not orbit initial conditions.
Eq. 16: M = m/(1 - m/3); a0 = g0 sum_{p=0..P} d_p M^p with g0 = M^(2/3); b_n = a_n/a0 = sum_{p=0..P}
c_{n,p} M^p; rho-bar = [sum_{n=1..N} (b_n + b_-n) cos 2n tau; sum (b_n - b_-n) sin 2n tau; 0]. "We compute
the coefficients up to a maximum order P = 9, so b_n must be determined for integers |n| <= N = 4
excluding n = 0. Note c_{n,p} = 0 for p < 2."

**Table 1: coefficients d_p for the HVO**

| p | d_p |
|---|---|
| 0 | 1 |
| 1 | -8/9 |
| 2 | 133/162 |
| 3 | -1264/2187 |
| 4 | 3319421/5038848 |
| 5 | -13366211/11337408 |
| 6 | 2028830887/2448880128 |
| 7 | -4682845907/5509980288 |
| 8 | 19228022393021/12694994583552 |
| 9 | -5982128249099224247/3119921868853739520 |

**Table 2: coefficients c_{n,p} for the HVO** (columns n = -4, -3, -2, -1, 1, 2, 3, 4; blank = 0)

| p | -4 | -3 | -2 | -1 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|---|---|---|
| 2 | 0 | 0 | 0 | -19/16 | 3/16 | 0 | 0 | 0 |
| 3 | 0 | 0 | 0 | -7/8 | 3/8 | 0 | 0 | 0 |
| 4 | 0 | 0 | 0 | 11/144 | 7/48 | 25/256 | 0 | 0 |
| 5 | 0 | 0 | 23/640 | 5/36 | -1/6 | 553/1920 | 0 | 0 |
| 6 | 0 | 1/192 | 207/3200 | -661/82944 | -34589/110592 | 3743/14400 | 833/12288 | 0 |
| 7 | 0 | 5237/215040 | 1829/288000 | 374797/276480 | -22907/46080 | -28811/864000 | 27337/107520 | 0 |
| 8 | 23/6144 | 263713/7526400 | 124719/40960000 | 98804551/37324800 | -23804639/24883200 | -332659139/1105920000 | 5056291/15052800 | 3537/65536 |
| 9 | 507317/28901376 | 38042489/4741632000 | 48459451/604800000 | 300079583/373248000 | -102469631/124416000 | -4857480211/7257600000 | 472019353/4741632000 | 11705987/48168960 |

For the Melnikov evaluation only (Section 3.1) a different normalisation M = m is used, giving "d0 = 1,
d1 = -2/3, d2 = 7/18, d3 = -4/81, c_{-1,2} = -19/16, c_{1,2} = 3/16, c_{-1,3} = -5/3, and c_{1,3} = 1/2",
the same values as in Scheeres [16]. Note the paper's own caution: "the definition of M given in (16b)
will be used for all calculations except when evaluating the Melnikov function."

Note on the printed text: the Section 3.1 sentence has the "M" definition garbled in the text layer; the
reading above (M = m for the Melnikov expansion, M = m/(1 - m/3) otherwise) follows Eq. 16b and the
accompanying coefficient values. Verify against the PDF before any reuse.

## Reproduction-target data
**None in the form of orbit initial conditions or periods.** The paper tabulates no initial states,
Jacobi-type constants, periods, or stability indices for any HR4BP orbit; families are shown only as
hodographs and singular-value plots, and the 9:2 NRHO, halo, butterfly and Lyapunov members are not
identified by amplitude. What a re-implementation could check: (a) the model constants m = 0.0808,
mu = 0.0122, forcing period pi, m_max = 0.19510486; (b) the HVO coefficient Tables 1 and 2 (these are
printed to full rational precision, and are mutually consistent with the Section 3.1 values by
construction of M = m versus M = m/(1 - m/3), which was not independently verified here); (c) the
qualitative counts: four HR4BP families from the L4 T* = 2 pi planar orbit, three object families from
the L1 T* = pi vertical orbit, the L2 northern halo T = pi continuing to the southern halo, and the
A, B, C L2 families of [18]. Any numerical match to the figures requires the full-size PDF; the figure
axes are not transcribed here.

## What it does NOT contain
- Cyclers, repeated lunar or planetary encounters, flybys, resonant hops or transfers: **none.** The
  words "cycler", "flyby" and "transfer" do not appear in the body. The only trajectory-type content is
  periodic orbits near libration points and the mention of "connections between Sun-Earth and EM
  libration point orbits [22]" in the introduction.
- Distant retrograde orbits and lunar-resonant (m:n) periodic-orbit families: not computed.
- Quasi-periodic orbits and tori: not computed here (cited only as existing work, refs [12, 13, 14, 15,
  19, 20]).
- Floquet multipliers, stability indices or station-keeping information for individual orbits: not given.
- Period-doubling or other period-multiplying branches: not computed ("new families were only computed at
  tangent bifurcation points").
- Ephemeris-model continuation or mission-design cost: not given.
- A statement of which families "break up" or an isola classification: not stated.
- Integration tolerances, software, run times: not stated.
- Comparison with the BCP, QBCP or CCR4BP periodic-orbit families: not made quantitatively.

## Relevance to our project (grounded only in the paper)
- It is a published statement that, in a pi-periodically forced EM model, CR3BP periodic orbits continue
  into the forced model only when their period is in resonance with the forcing (T = b Tg), only at
  specific phase points on the orbit (Melnikov zeros), and that one CR3BP orbit can yield several
  inequivalent forced-model orbits ("dynamical equivalents"). For our BCR4BP and QBCP work the analogous
  forcing period would be the synodic period of the Sun; the paper does not address those models beyond
  the statement that similar propositions are expected ("such as the ER3BP, the BCP, and the QBCP",
  Section 3.2).
- It is not a cycler paper and supplies no reproduce-and-validate target for any catalogue row.

## Relevant bibliography (as printed in the paper; four-body periodic-orbit, quasi-periodic and Melnikov literature it builds on)
- [6] Huang, S.-S.: Very Restricted Four-Body Problem. TN D-501, NASA (Sep 1960).
- [7] Gomez, G., Jorba, A., Masdemont, J., Simo, C.: Normal form of the bicircular model and related
  topics. In: Dynamics and Mission Design Near Libration Points, pp. 53-110. World Scientific (2001).
- [8] Rosales, J.J.: On the effect of the Sun's gravity around the Earth-Moon L1 and L2 libration points.
  PhD thesis, Universitat de Barcelona (2020).
- [9] Jorba-Cusco, M., Farres, A., Jorba, A.: Two Periodic Models for the Earth-Moon System. Frontiers in
  Applied Mathematics and Statistics 4 (2018).
- [10] Simo, C., Gomez, G., Jorba, A., Masdemont, J.: The Bicircular Model Near the Triangular Libration
  Points of the RTBP. In: From Newton to Chaos, pp. 343-370. Springer (1995).
- [11] Boudad, K.K., Howell, K.C., Davis, D.C.: Dynamics of synodic resonant near rectilinear halo orbits
  in the bicircular four-body problem. Advances in Space Research 66(9), 2194-2214 (2020).
- [12] Castella, E., Jorba, A.: On the vertical families of two-dimensional tori near the triangular points
  of the Bicircular problem. Celestial Mechanics and Dynamical Astronomy 76, 35-54 (2000).
- [13] Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Families of Halo-like invariant tori around L2 in the
  Earth-Moon Bicircular Problem. Celestial Mechanics and Dynamical Astronomy 133 (2021).
- [14] Andreu, M.A.: The Quasi-bicircular Problem. PhD thesis, Universitat de Barcelona (1998).
- [15] Rosales, J.J., Jorba, A., Jorba-Cusco, M.: Invariant manifolds near L1 and L2 in the quasi-bicircular
  problem. Celestial Mechanics and Dynamical Astronomy 135 (2023).
- [16] Scheeres, D.J.: The Restricted Hill Four-Body Problem with Applications to the Earth-Moon-Sun
  System. Celestial Mechanics and Dynamical Astronomy 70(2), 75-98 (1998).
- [17] Peterson, L.T., Rosales, J.J., Scheeres, D.J.: The vicinity of Earth-Moon L1 and L2 in the Hill
  restricted 4-body problem. Physica D: Nonlinear Phenomena 455 (2023).
- [18] Olikara, Z.P., Gomez, G., Masdemont, J.J.: A Note on Dynamics About the Coherent Sun-Earth-Moon
  Collinear Libration Points. In: Gomez, G., Masdemont, J.J. (eds.) Astrodynamics Network AstroNet-II,
  pp. 183-192. Springer (2016).
- [19] Henry, D., Rosales, J., Brown, G., Peterson, L., Scheeres, D.: Quasi-Periodic Orbits near Earth-Moon
  L1 in the Hill Restricted Four-Body Problem. In: 34th ISTS (2023).
- [20] Henry, D.B., Rosales, J., Brown, G.M., Scheeres, D.J.: Quasi-Periodic Orbits near Earth-Moon L1 and
  L2 in the Hill Restricted Four-Body Problem. In: AAS/AIAA Astrodynamics Specialist Conference (2023).
- [21] Peterson, L.T., Jorba, A., Brown, G.M., Scheeres, D.J.: Dynamics Around the Earth-Moon Triangular
  Points in the Hill Restricted 4-Body Problem. Communications in Nonlinear Science and Numerical
  Simulation (2024). In Preparation.
- [22] Olikara, Z.P., Scheeres, D.J.: Mapping Connections Between Planar Sun-Earth-Moon Libration Orbits.
  In: 27th AAS/AIAA Space Flight Mechanics Meeting (2017).
- [23] Sanaga, R.R., Howell, K.C.: Synodic Resonant Halo Orbits in the Hill Restricted Four-Body Problem.
  In: 33rd AAS/AIAA Space Flight Mechanics Meeting (2023).
- [2] Whitley, R.J., Davis, D.C., Burke, L.M., McCarthy, B.P., Power, R.J., McGuire, M.L., Howell, K.C.:
  Earth-Moon Near Rectilinear Halo and Butterfly Orbits for Lunar Surface Exploration. In: AAS/AIAA
  Astrodynamics Conference (2018).
- [3] Szebehely, V.: Theory of Orbits: The Restricted Problem of Three Bodies. Academic Press (1967).
- [4] Park, B., Howell, K.C.: Leveraging Intermediate Dynamical Models for Transitioning from the Circular
  Restricted Three-Body Problem to an Ephemeris Model. In: AAS/AIAA Astrodynamics Specialist Conference
  (2022).
- [5] Peng, H., Bai, X.: Natural deep space satellite constellation in the Earth-Moon elliptic system. Acta
  Astronautica 153, 240-258 (2018).
- [24] Melnikov, V.K.: On the stability of the center for time periodic perturbations. Transactions of the
  Moscow Mathematical Society 12, 1-57 (1963).
- [34] Cenedese, M., Haller, G.: How do conservative backbone curves perturb into forced responses? A
  Melnikov function analysis. Proceedings of the Royal Society A 476(2234) (2020).
- [35] Rhouma, M.B.H., Chicone, C.: On the Continuation of Periodic Orbits. Methods and Applications of
  Analysis 7(1), 85-104 (2000).
- [36] Scheeres, D.J.: Orbital Motion in Strongly Perturbed Environments. Springer (2012).
- [38] Wintner, A.: The Analytical Foundations of Celestial Mechanics. Dover (1952).
- Melnikov-theory background only: [25] Greenspan and Holmes (1981); [26] Wiggins (2003); [27]
  Guckenheimer and Holmes (1983); [28] Perko (2001); [29] Haller, Chaos Near Resonance (1999); [30] Guo
  et al. (2022); [31] Veerman and Holmes (1985); [32] Yagasaki (1996); [33] Polcar and Semerak (2019);
  [37] Holmes, Physics Reports 193(3) (1990). [1] is the NASA Strategic Plan 2022.

## Status
Digested from the arXiv v1 preprint; the published SIADS version was not examined. No catalogue row is
affected. Candidate corpus-index description: Brown et al. 2024 (arXiv v1 of SIADS 24(1) 346-375, 2025):
periodic-orbit families of EM libration-point orbits continued from the CR3BP into the Sun-Earth-Moon
Hill restricted 4-body problem (m = 0.0808, mu = 0.0122) using a Melnikov-type function and
pseudo-arclength continuation with SVD bifurcation detection; no cyclers, no initial-condition tables.
