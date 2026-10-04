# Digest: Singh, Park & Howell (2026), "Evolution of L2 orbit families and bifurcations within intermediary-fidelity models"

AAS 26-654 (2026), Purdue University (J. E. Singh, M.S. student; B. Park, Apollo 11 Postdoctoral Fellow; K. C. Howell), 22
pages, conference paper (no journal data printed on the page).
Filed in the private paper corpus as
singh-park-howell-2026-evolution-l2-orbit-families-bifurcations-intermediary-fidelity-models-qbcp-er3bp-AAS-26-654.pdf

Digested 2026-10-04. Each statement is marked READ (seen on the page, with page, section, equation or table) or INFERRED (our
reading or computation). Page numbers are the paper's printed page numbers (1 to 22). All tables were transcribed from the
PDF text layer and checked against the page image (tables 2 and 3 on p12, table 4 on p14, table 1 on p5). Context: tasks #891,
#892 and #884; the paper bears on the quasi-bicircular model (`core/qbcp.py`) and on what is published about periodic and
quasi-periodic orbits in intermediary-fidelity Earth-Moon models.

## 0. What the paper is

READ (abstract, p1): the paper constructs "periodic and quasi-periodic orbits in the Quasi-Bicircular Problem (QBCP) and the
Elliptic Restricted Three-Body Problem (ER3BP) focusing on the L2 region. Orbits within the same CR3BP family differ in
continuation behavior depending on the intermediary-fidelity model. Orbit stability characteristics evolve with the
increased fidelity. Bifurcation analysis into orbit family continuums characterizes a global drift in pitchfork bifurcations
that generate CR3BP structures under more complex dynamics."

READ: the paper prints no initial condition, no Jacobi constant, no stability index, no perilune radius and no table of orbit
states. Its printed numerical results are the three bifurcation tables (tables 2 to 4, section 5 below), a handful of periods
quoted in the text, the model parameters `e = 0.0554` and `omega_S ~ 0.925196`, and the discretisation `N = 51`, `M = 5`. The rest
is in figures (continuation continuums, Lyapunov-exponent plots, geometries) without tabulated data.

## 1. Models, frame, constants, conventions

### 1.1 Frame (READ, p2)

"This investigation adopts the Earth-Moon rotating-pulsating frame (RPF) for all model definitions and trajectory
generation." Origin at the Earth-Moon barycentre, x-hat from the Earth to the Moon, z-hat along the Earth-Moon angular
momentum, y-hat completing the triad; length normalised by "the instantaneous Earth-Moon distance, l_EM, within the
respective model". Positions (eqs. 1 to 4): `r_E = -mu x-hat`, `r_M = (1 - mu) x-hat`, `r_S = x_S x-hat + y_S y-hat + z_S z-hat`.
`mu = mu_M/(mu_M + mu_E)`. This is the same frame as the project's (Earth at -mu, Moon at 1 - mu). The value of mu, the value of
l_EM (nominal) and the dimensional gravitational parameters are not printed. INFERRED (section 5.2 below): the tables'
periods are consistent with a nominal Earth-Moon distance of 384,400 km and the Earth plus Moon gravitational parameter
(about 403,503 km^3/s^2).

### 1.2 CR3BP (READ, p3, eqs. 5, 6)

`x'' = 2y' + dU*/dx, y'' = -2x' + dU*/dy, z'' = dU*/dz` with `U* = (x^2 + y^2)/2 + (1 - mu)/|r - r_E| + mu/|r - r_M|`; time
`t = T sqrt((mu_M + mu_E)/l_EM^3)`. One integral of motion, the Jacobi constant.

### 1.3 HFEM (READ, p3, eqs. 7, 8)

Point masses Earth, Moon, Sun, ephemerides from DE440 ("NASA JPL data, DE440.bsp"), Moon-centred J2000 frame; the RPF state is
obtained by `R = R_B + l_EM C r`. Not used for tables; only in the section on ephemeris transitioning (section 7 below).

### 1.4 ER3BP (READ, p3 to p4, eqs. 9, 10)

True anomaly `t` of the primary orbit as independent variable, so the period is `P_sys = 2 pi`:

    x'' = 2y' + dU_E*/dx,   y'' = -2x' + dU_E*/dy,   z'' + z = dU_E*/dz
    U_E* = (1/(1 + e cos t)) ( (1/2)(x^2 + y^2 + z^2 + mu(1 - mu)) + (1 - mu)/|r - r_E| + mu/|r - r_M| )

(as printed; the parenthesis closes after the three terms). Eccentricity is the continuation parameter, `0 <= e <= 0.0554` ("e =
0.0554 for the Earth-Moon ER3BP", p4; table 1, p5). Reference for the model: Szebehely and Grebenikov 1967 (ref 16).

### 1.5 QBCP (READ, p4, eqs. 11 to 15)

Taken from Andreu's thesis (ref 17, "dynamically coherent motion for the primaries"): "the Sun, Earth, and Moon motions follow
solutions to their mutual planar three-body problem represented as Fourier series". Equations of motion in the RPF, as
printed in the acceleration form (ref 18, Gao, Masdemont, Gomez and Yuan 2022; "with the frame rotated to align the RPF as
described in this investigation"):

    x'' = a1 a4 + (a1'/a1) x' + 2 a3 y' + ( -a1' a2/a1 + a2' + a2^2 + a3^2 ) x + ( -a1' a3/a1 + a3' ) y + a1 a6 dU_Q*/dx   (11)
    y'' = a1 a5 - 2 a3 x' + (a1'/a1) y' + ( a1' a3/a1 - a3' ) x + ( -a1' a2/a1 + a2' + a2^2 + a3^2 ) y + a1 a6 dU_Q*/dy       (12)
    z'' = (a1'/a1) z' + ( -a1' a2/a1 + a2' + a2^2 ) z + a1 a6 dU_Q*/dz                                                       (13)
    U_Q* = (1 - mu)/|r - r_E| + mu/|r - r_M| + mu_S/|r - r_S|,    r_S = -a7 x-hat - a8 y-hat                                (14)

(`a_i` = alpha_i, the primes are derivatives with respect to the model time.) Notes:
- READ: in (11) to (13) the factor `alpha_1 alpha_6` multiplies the gradient of the full `U_Q*`, which contains the Earth, Moon
  and Sun terms. INFERRED: this is the same placement of alpha_6 on all three terms that the project's `#892` correction
  restored (alpha_6 on the whole Newtonian potential, not on the Sun term alone), so this paper is a second published
  source for that placement besides Jorba-Cusco, Farres and Jorba (2018).
- READ: the Sun position in the RPF is `r_S = -alpha_7 x-hat - alpha_8 y-hat`. The sign convention of the Fourier functions
  `alpha_7`, `alpha_8` is not given on the page; the paper says nothing explicit about the sense (clockwise or counter-clockwise)
  of the Sun in the RPF beyond calling omega_S "the nondimensional mean angular velocity of the Sun in the RPF".
- READ: `mu_S` is "the nondimensional gravitational parameter of the Sun based on the characteristic gravitational parameter,
  mu* = mu_M + mu_E". Its value is not printed.
- READ: independent variable `t = T sqrt((mu_M + mu_E)/l_EM^3)` "observing that l_EM is time-varying in this system".
- READ: "the dynamics in the QBCP repeat after one synodic month or a period of P_sys = 2 pi/omega_S where omega_S ~ 0.925196 is
  the nondimensional mean angular velocity of the Sun in the RPF." (p4). Matches the bicircular value 0.925195985518 of the
  Jorba papers to the six digits printed.
- READ: the Fourier coefficients alpha_i "require a semi-analytic procedure to construct the relative motions of the Earth,
  Moon, and Sun in the RPF (see Appendix in Singh and Howell 7)", where reference 7 is J. Singh and K. C. Howell, "Isolating
  Perturbing Cislunar Effects on Orbits Through a Homotopy-based Unified Transitioning Scheme", AAS/AIAA Astrodynamics
  Specialist Conference, Whistler, 2026, paper 26-694. The number of Fourier terms and the numerical alpha tables are NOT in
  this paper and it does not say whether they come from Andreu's tables, Jorba-Cusco et al.'s Table 4 or a fresh computation.
  Say plainly: not stated.
- READ: no continuation parameter as simple as eccentricity exists; an artificial homotopy "linearly blends the CR3BP and QBCP
  accelerations" (eq. 15): `x'' = eps x''_QBCP + (1 - eps) x''_CR3BP` and likewise for y and z, `0 <= eps <= 1`; "Intermediary
  homotopy values serve purely as numerical scaffolding... and alone hold no physical meaning" (p5). Footnote: this differs
  from Andreu's homotopy, which "linearly connects the Hamiltonians of the CR3BP to the QBCP". INFERRED: a continuation
  value reached at `eps < 1` is not a statement about a Sun with a smaller mass.

### 1.6 Apsidal configurations and phase (READ, p5 to p7)

Symmetries in the non-autonomous models: S1 (z reflection, no time reversal, valid at any epoch), S2 and S3 (reverse-time
reflections, valid only at an apse time `tau = j P_sys/2`, `j` a natural number). "These apsidal configurations correspond to
the Moon at perigee/apogee (ER3BP) or the Sun and Moon in conjunction/opposition (QBCP)." Even `j` and odd `j` give the two
configurations; which of conjunction and opposition belongs to even `j` is not stated. The paper gives no explicit phase
convention for the Sun's angle at t = 0 beyond that.

## 2. Orbit families, resonances, selection and continuation

### 2.1 Families (READ, section "Models" p5 to p6, figs. 1, 2)

All are CR3BP orbits about L2: Lyapunov and vertical (doubly symmetric: OX and XOZ), northern halo and northeast axial (singly
symmetric, each with a mirror family by S1). "Each CR3BP orbit family is monotonic in period for the selected ranges." The
NRHO is not treated separately: the paper calls the near-Moon end of the halo family the "NRHO region" (p14). No butterfly,
dragonfly or other bifurcating families, and nothing about L1, are in the paper.
Bifurcation structure in the CR3BP (READ, p8): "two pitchfork bifurcations occur along the L2 Lyapunov family where one
produces the halos and the other results in the axials. The verticals and axials are also linked by a pitchfork bifurcation at
the opposite end of the axial family."

### 2.2 Resonances (READ, p6, p8)

Periodic orbits in the intermediary models exist only when the CR3BP orbit's period is commensurate with the system period:
"the period of the CR3BP orbits must be commensurate with the system period defined by the resonance ratio, p : q, where
p and q are coprime positive integers; the total period of the orbit then results in q P_sys". The resonance is with the
system period (the Moon's sidereal period `2 pi` in true anomaly for the ER3BP, ie about 27.3 days; the Sun-Moon synodic month
`2 pi/omega_S` for the QBCP, ie about 29.5 days), not with the Moon or the Sun as a body one flies by. INFERRED from the tables
(section 5.2) that `P = q P_sys / p` is the CR3BP period of the p:q orbit, so that p counts orbit revolutions per q system
periods. Non-resonant orbits continue as 2D tori.

### 2.3 Selection of members for continuation (READ, p9 to p10, p6)

"Members of the CR3BP L2 libration point orbit families are continued as either POs (Eq. (21)) or QPOs (Eq. (25)) into both
intermediary-fidelity models." The sampled members are indexed by the CR3BP period; how many members, how spaced and their
range are only visible in the figures (Lyapunov roughly 15 to 21 days from fig. 6, halo roughly 6 to 14.8 days from fig. 8; the vertical
and axial ranges are in figs. 7 and 10 and were not read off). The sampling of periodic orbits is by the resonant ratios p:q (a list of ratios such as 24:13,
11:6, 89:61 appears in the tables). INFERRED: they appear to be chosen by scanning ratios across the period range; the
selection rule is not stated.

### 2.4 Continuation methods (READ, p6 to p9)

- Periodic orbits: fixed-time multiple shooting (Newton-Raphson) with `n = 2p + 1` segments for the half-period scheme,
  perpendicular-crossing symmetry constraints (eqs. 19 to 22), continuation variable `c` appended as a free variable, pseudo
  arclength constraint (eq. 23). Initial crossing at an apse time `tau = j P_sys/2`; two counterparts "A" and "B" per resonant
  CR3BP orbit (counterpart B starts at `t1 = P_sys/2` when p is odd, or at the opposite crossing when p is even; after Park
  and Howell, ref 9). Verticals use a quarter-period scheme, and only even-q resonant verticals keep double symmetry.
- Quasi-periodic orbits: Gomez-Mondelo / Olikara-Scheeres (GMOS) invariant-curve method with `N = 51` nodes and `M = 5` maps,
  stroboscopic time = the CR3BP period `P`, rotation angle `rho = 2 pi P/P_sys` (eq. 26 and following). QPOs exist "as
  continuous families in the ER3BP and QBCP, unlike POs, save for resonance gaps".
- Stability: monodromy eigenvalues rebuilt from the multiple-shooting segment STMs with a cyclic block eigenproblem (appendix,
  eqs. 28 to 33, after Jorba and Rosales), with symmetric reduction; for QPOs the eigenvalues of the stroboscopic return
  (eqs. 34 to 37).
- Out-of-plane stability is shown as a Lyapunov exponent `Gamma = ln|lambda| / P` (eq. 27), `P` the orbital period in the model
  (`q P_sys` for p:q POs).
- Continuation variable `c` is `e` for the ER3BP and `eps` for the QBCP (table 1).

### 2.5 What happens to each family (READ, section "Orbit family continuums", pp9 to 17)

The paper names four reasons a continuation stops: numerical challenges of GMOS (NC), resonance gaps around low-order
p:q (D1), a resonance-driven fold (D2) and bifurcation drift (D3).
- Lyapunov (p10 to p12): "All POs in the sampled range successfully reach the target values for the ER3BP and QBCP." Gaps
  around resonances for QPOs. The 1:1 resonant Lyapunov (about 27.3 days) is "strenuous to GMOS". The out-of-plane stability
  sequence is the same in all three models: centre, then saddle at the halo bifurcation, saddle until the axial bifurcation,
  then centre.
- Vertical (p12 to p14): nearly all resonant POs continue; the 19:11 resonant vertical (both OX and XOZ) reaches only
  `eps ~ 0.25` in the QBCP before turning (D2); one small region around the vertical of period 16.1415 days challenges both POs
  and QPOs. "most continuation processes locate a suitable analog in the ER3BP and QBCP, including the vicinity of axial
  bifurcation within the CR3BP (19.2033 days)."
- Halo (p14 to p15, fig. 8): shows all three mechanisms. In the ER3BP (after Park and Howell, fig. 25 of ref 9) three regions: a
  regular region near the Lyapunovs, an NRHO region near the Moon, and an interface region between them where continuation
  fails (D2); the 2:1 halo (`P = 13.6605` days) sits in a resonance gap (D1). In the QBCP "lacks distinct regions throughout the
  evaluated period range... no PO demonstrates a limit in the continuation process anywhere in the family", with gaps around two
  resonant halos (D1). The QBCP therefore does not reproduce the HFEM interface region.
- Axial (p15 to p17, figs. 10 to 12): a narrow period range with no low-q resonant POs, hence no large gaps. A "slant region"
  near the vertical bifurcation where the achieved continuation value trails to zero; the 20:13 axial in the QBCP is shown
  to be a broken bifurcation: the continuation "diverts to a vertical at a specific homotopy parameter whereupon the
  structure eventually reaches the QBCP" (p16). In the regular region the QBCP axials are further out of plane than their
  CR3BP counterparts. The lowest-period axials near the Lyapunov end cannot be reached by the fixed-period continuation;
  continuation in the stroboscopic period at fixed `c` locates them (p17).

## 3. Printed numbers

### 3.1 Table 1 (p5), model parameters

| Model | Perturbing effect | Independent variable | System period | Continuation variable | Range |
| --- | --- | --- | --- | --- | --- |
| ER3BP | Lunar eccentricity | True anomaly | 2 pi | Eccentricity e | 0 <= e <= 0.0554 |
| QBCP | Solar gravity | Non-dimensional time | 2 pi / omega_S | Homotopy parameter eps | 0 <= eps <= 1 |

### 3.2 Table 2 (p12), L2 Lyapunov, halo bifurcation (out-of-plane eigenvalue pair, centre to saddle)

| Model | Orbit type | Period (days) | p:q resonance | Stability change |
| --- | --- | --- | --- | --- |
| CR3BP | PO | 14.8319 | - | Center to Saddle |
| ER3BP | PO | 14.7792 - 14.8825 | 24:13 - 11:6 | Center to Saddle |
| ER3BP | QPO | ~14.8252 | - | Center to Saddle |
| QBCP | PO | 14.7453 - 14.8672 | 2:1 - 121:61 | Center to Saddle |
| QBCP | QPO | ~14.8136 | - | Center to Saddle |

### 3.3 Table 3 (p12), L2 Lyapunov, axial bifurcation (saddle to centre)

| Model | Orbit type | Period (days) | p:q resonance | Stability change |
| --- | --- | --- | --- | --- |
| CR3BP | PO | 18.7183 | - | Saddle to Center |
| ER3BP | PO | 18.7007 - 18.7094 | 89:61 - 35:24 | Saddle to Center |
| ER3BP | QPO | ~18.7047 | - | Saddle to Center |
| QBCP | PO | 18.5370 - 18.5682 | 35:22 - 27:17 | Saddle to Center |
| QBCP | QPO | ~18.5418 | - | Saddle to Center |

### 3.4 Table 4 (p14), L2 vertical, axial bifurcation (centre to saddle)

| Model | Orbit type | Period (days) | p:q resonance | Stability change |
| --- | --- | --- | --- | --- |
| CR3BP | PO | 19.2033 | - | Center to Saddle |
| ER3BP | PO | 19.1730 - 19.3266 | 37:26 - 24:17 | Center to Saddle |
| ER3BP | QPO | ~19.2000 | - | Center to Saddle |
| QBCP | PO | 18.9185 - 19.0262 | 53:34 - 31:20 | Center to Saddle |
| QBCP | QPO | ~19.0029 | - | Center to Saddle |

READ (p11, p12): "As POs exist at discrete periods, the exact bifurcating location is not supplied within the intermediate
models. Rather, the POs supply bounds on the bifurcating period where one end renders one stability type and the opposite end
yields the other... Within these bounds, QPO stability properties are assessed to further refine the locations." So the PO
columns are brackets and the QPO rows are the paper's estimates (marked with a tilde). INFERRED (our observation): the ER3BP
bracket in table 4 (19.1730 - 19.3266) contains the CR3BP value 19.2033; only the QPO estimate (about 19.2000) sits below
it, by 0.0033 days. The ER3BP axial bifurcation (18.7007 - 18.7094) is a bracket about 0.009 day wide at a CR3BP value of
18.7183, outside it. The QBCP bracket in table 2 contains the CR3BP value (14.8319 within 14.7453 - 14.8672); the QBCP brackets in tables 3
and 4 lie wholly below the CR3BP values (18.7183 and 19.2033).

### 3.5 Periods quoted in the text

| Where | Quantity | Value |
| --- | --- | --- |
| p10 | 1:1 resonant Lyapunov (about) | 27.3 days |
| p12 and p13 | CR3BP halo bifurcation; axial bifurcation | 14.8319; 18.7183 days (table 2, 3) |
| p13 | resonant vertical with the problematic QPO and PO region | 16.1415 days |
| p13 | CR3BP vertical at the axial bifurcation | 19.2033 days |
| p13 | 19:11 vertical in QBCP | stops at eps about 0.25 |
| p14 | 2:1 halo (as stated in the text) | 13.6605 days |
| p14 | halos "4:1 and 3:1" respectively | 9.8302 and 7.3727 days |
| p15 | halo near the Lyapunov bifurcation, fig. 9 | 14.8073 days |
| fig. 4 caption | L2 halo continued as a QPO, N = 51 | 12.378 days |
| p18, fig. 13b | axial example | 18.7385 days |
| p18 | lowest-period axial reference | 18.7183 days |
| p18 | fraction of axial family lacking a CR3BP counterpart "if the overall low range in period ... is explored" | roughly 40% |

## 4. Statements on pitchfork bifurcations drifting between models (READ, quoted)

- Abstract: "Bifurcation analysis into orbit family continuums characterizes a global drift in pitchfork bifurcations that
  generate CR3BP structures under more complex dynamics."
- Introduction p2: "Tracking the orbit period across these families, bifurcation analysis uncovers a systematic shift of the
  pitchfork bifurcations that generate the halo and axial families toward lower periods under perturbation, a mechanism
  termed bifurcation drift."
- p12: "In both halo and axial cases, the CR3BP bifurcation occurs at a higher period compared to its intermediary-fidelity
  counterparts. This drift in the bifurcation period is more pronounced in the QBCP as opposed to the ER3BP for both cases.
  The drift is also more notable at the axial bifurcation compared to the halo bifurcation."
- p13 on verticals: "the lower bifurcating period observed in the ER3BP and QBCP (Table 4) contracts the period range
  admitting axial analogs in the intermediate models."
- Conclusion p18: "solar gravity and lunar eccentricity independently decrease the period associated with the bifurcations that
  produce the halo and the axial families. In all cases, solar gravity reduces the bifurcation period more significantly than
  lunar eccentricity."
- INFERRED (our arithmetic on the tables, QPO estimate against CR3BP): halo bifurcation shifts by -0.0067 day (ER3BP) and
  -0.0183 day (QBCP); axial bifurcation by -0.0136 day (ER3BP) and -0.1765 day (QBCP); vertical-axial bifurcation by -0.0033 day
  (ER3BP) and -0.2004 day (QBCP). The relative size of the QBCP axial and vertical shifts, about 0.9 to 1.0 percent, against
  about 0.12 percent for the halo bifurcation, is what the text calls "more notable". On the PO brackets alone the ER3BP
  vertical bracket straddles the CR3BP value; the "lower period" statement for the ER3BP vertical rests on the QPO estimate.
  Reported as an observation.

## 5. Other observations on the printed material

### 5.1 Respectful note on two apparent slips (INFERRED, factual)

- p14: "Some prominent resonance gaps appear that surround the 4:1 and 3:1 (P = 9.8302 days and P = 7.3727 days respectively)
  resonant halo orbits." With the QBCP system period `P_sys = 2 pi/omega_S` of about 29.49 days (section 5.2), a p:1 orbit has
  `P = P_sys/p`: 9.8302 days is `P_sys/3` and 7.3727 days is `P_sys/4`. So the two periods are the 3:1 and the 4:1 orbits
  respectively, the reverse of the order in the sentence. This reading assumes the printed p:q convention of the tables.
- p14: the 2:1 halo "P = 13.6605 days" in the ER3BP is `27.3210/2`; the ER3BP periods implied by the table 2 to 4 brackets
  correspond to `P_sys` = 27.2846 days (13.6423 for 2:1). The two differ by 0.13 percent. This may reflect a different
  constant for the sidereal month between the two places (27.3217 days is the sidereal month itself); it is a rounding or
  constants matter, not a problem for the main results.

### 5.2 Cross-check of the tables' internal consistency (INFERRED, our computation)

Every PO bracket end in tables 2 to 4 satisfies `P = q P_sys / p` with `P_sys(ER3BP) = 27.2846 d` and `P_sys(QBCP) = 29.4906 d`,
to within about 1e-5 relative, for all 12 resonant entries (ER3BP: 24:13 gives P_sys = 27.28468, 11:6 gives 27.28458, 89:61
gives 27.28463, 35:24 gives 27.28454, 37:26 gives 27.28465, 24:17 gives 27.28461; QBCP: 2:1 gives 29.4906, 121:61 gives 29.49068, 35:22
gives 29.49068, 27:17 gives 29.49067, 53:34 gives 29.49060, 31:20 gives 29.49061). This confirms the transcription of the table
digits and the `P = q P_sys/p` convention. The same values follow from `l_EM = 384,400 km` and `mu_E + mu_M = 403,503.2 km^3/s^2`
(time unit 4.34248 days, `2 pi` time units 27.2846 days, `2 pi/omega_S` time units 29.4906 days), so INFERRED the authors use a
nominal rather than ephemeris Earth-Moon distance and a nominal gravitational parameter.

### 5.3 Fig. 13 and the HFEM comparison (p17 to p18)

READ: axial family generated in the HFEM with 20-year multi-year trajectories by the Unified Transition Scheme (Sanaga, Park and
Howell, ref 25); "the HFEM results reflect the QBCP results closely perhaps hinting at the strength of solar gravity in this
region"; "For a given axial period, the CR3BP orbit is the least geometrically similar to the HFEM"; "roughly 40% of the
family lacks a suitable CR3BP counterpart with the same period". The ephemeris epoch and dates are not printed.

### 5.4 Cycler, resonant lunar-flyby or Earth-Moon transfer content: none

READ, stated plainly: the paper does not treat any cycler, any resonant orbit with lunar flybys, any Earth-Moon transfer orbit,
any heteroclinic or homoclinic connection, or any invariant-manifold transfer. It is about the existence, continuation and
stability of orbits near L2 (Lyapunov, vertical, halo including its NRHO end, axial) in the QBCP and ER3BP. The only
resonances in it are between the orbit period and the model's system period (about 27.3 days for the ER3BP, 29.5 days for the QBCP),
which concern the structure of the period-`q P_sys` periodic orbits, not trajectories that return to the Moon. The mention
"orbit chain trajectories" (reference 28, Gomez and Howell, AAS 26-874, 2026, "Anatomical Analysis of Orbit Chain Trajectories
in the Earth-Moon System") is a citation for the multiple-shooting STM decomposition; its content is not described here.

## 6. References cited

Full citations, as printed (reference numbers from the paper):
- [4] E. M. Zimovan-Spreen, K. Howell, and D. C. Davis, "Near rectilinear halo orbits and nearby higher-period dynamical
  structures: orbital stability and resonance properties," Celestial Mechanics and Dynamical Astronomy, Vol. 132, No. 28, 2020,
  10.1007/s10569-020-09968-2.
- [5] D. C. Davis, S. M. Phillips, S. Vutukuri, B. P. McCarthy, and K. C. Howell, "Stationkeeping and Transfer Trajectory Design
  for Spacecraft in Cislunar Space," AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, No. 17-826, 2017.
- [6] K. Boudad, K. C. Howell, and D. Davis, "Analogs for Earth-Moon halo orbits and their evolving characteristics in
  higher-fidelity force models," AIAA SCITECH 2022 Forum.
- [7] J. Singh and K. C. Howell, "Isolating Perturbing Cislunar Effects on Orbits Through a Homotopy-based Unified Transitioning
  Scheme," AAS/AIAA Astrodynamics Specialist Conference, Whistler, BC, No. 26-694, 2026. (Holds the QBCP Fourier construction.)
- [8] G. Gomez, J. Masdemont, and J. Mondelo, "Solar system models with a selected set of frequencies," Astronomy and
  Astrophysics, Vol. 390, No. 2, 2002, pp. 733-749, 10.1051/0004-6361:20020625.
- [9] B. Park and K. C. Howell, "Characterization of Earth-Moon L2 halo analogs in an ephemeris model using the elliptic
  restricted three-body problem," Advances in Space Research, Vol. 75, No. 6, 2025, pp. 5078-5109.
- [10] R. Sanaga and K. C. Howell, "Leveraging the Hill restricted four-body problem to investigate the ephemeris transition
  characteristics in the Earth-Moon L2 halo orbit region," Astrodynamics, Vol. 9, No. 5, 2025, pp. 785-805.
- [11] B. Park, R. Sanaga, and K. C. Howell, "Numerical Assessment of a Frequency-Based Hierarchy for the Cislunar Domain,"
  Journal of Guidance, Control, and Dynamics, Vol. 48, No. 11, 2025, pp. 2462-2479.
- [12] X. Leng and H. Lei, "High-order expansions of multi-revolution elliptic Halo orbits in the elliptic restricted
  three-body problem," CMDA, Vol. 138, No. 4, 2026, 10.1007/s10569-026-10276-4.
- [13] H. Peng and S. Xu, "Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted
  three-body problem," CMDA, Vol. 123, 2015, 10.1007/s10569-015-9635-2.
- [14] R. Broucke, "Stability of periodic orbits in the elliptic, restricted three-body problem," AIAA J., Vol. 7, No. 6, 1969,
  10.2514/3.5267.
- [16] V. Szebehely and E. Grebenikov, "Theory of Orbits - The Restricted Problem of Three Bodies," Soviet Astronomy, Vol. 13,
  1967, p. 364.
- [17] M. Andreu, The Quasi-Bicircular Problem. PhD Dissertation, Universitat de Barcelona, 1998.
- [18] C. Gao, J. J. Masdemont, G. Gomez, and J. Yuan, "The web of resonant periodic orbits in the Earth-Moon Quasi-Bicircular
  Problem including solar radiation pressure," Communications in Nonlinear Science and Numerical Simulation, Vol. 111, 2022,
  p. 106480, https://doi.org/10.1016/j.cnsns.2022.106480.
- [19] S. Campagnola, M. Lo, and P. Newton, "Subregions of motion and elliptic halo orbits in the elliptic restricted
  three-body problem," AAS/AIAA Space Flight Mechanics, Galveston, TX, No. 08-200, 2008.
- [20] G. Brown, L. Peterson, D. B. Henry, and D. Scheeres, "Structure of Periodic Orbit Families in the Hill Restricted 4-Body
  Problem," SIAM J. Applied Dynamical Systems, Vol. 24, No. 1, 2025, pp. 346-375.
- [21] G. Gomez and J. M. Mondelo, "The dynamics around the collinear equilibrium points of the RTBP," Physica D, Vol. 157,
  No. 4, 2001, pp. 283-321.
- [22] Z. P. Olikara and D. J. Scheeres, "Numerical method for computing quasi-periodic orbits and their stability in the
  restricted three-body problem," Advances in the Astronautical Sciences, Vol. 145, 2012, pp. 911-930.
- [23] Z. P. Olikara, G. Gomez, and J. J. Masdemont, "A Note on Dynamics About the Coherent Sun-Earth-Moon Collinear Libration
  Points," Astrodynamics Network AstroNet-II, Astrophysics and Space Science Proceedings, Vol. 44, 2016, pp. 183-192,
  10.1007/978-3-319-23986-6-13.
- [24] J. J. Rosales, A. Jorba, and M. Jorba-Cusco, "Families of Halo-like invariant tori around L2 in the Earth-Moon Bicircular
  Problem," CMDA, Vol. 133, No. 16, 2021, 10.1007/s10569-021-10012-0.
- [25] R. R. Sanaga, B. Park, and K. C. Howell, "Systematic Numerical Transitions of Cislunar Trajectories Across Dynamical
  Models," J. Astronautical Sciences, Vol. 73, No. 1, 2026, p. 10.
- [26] B. Park, A Frequency-Based Model Hierarchy within Cislunar Space. PhD dissertation, Purdue University, 2025.
- [27] A. Jorba, "Numerical computation of the normal behaviour of invariant curves of n-dimensional maps," Nonlinearity, Vol. 15,
  No. 5, 2001, p. 943, 10.1088/0951-7715/14/5/303 (as printed; the volume, issue and DOI differ).
- [28] R. J. Gomez and K. C. Howell, "Anatomical Analysis of Orbit Chain Trajectories in the Earth-Moon System," AAS/AIAA
  Astrodynamics Specialist Conference, Whistler, BC, No. 26-874, 2026.
- [29] B. McCarthy and K. C. Howell, "Leveraging Quasi-Periodic Orbits for Trajectory Design in Cislunar Space," Astrodynamics,
  Vol. 5, No. 2, 2021, pp. 139-165, 10.1007/s42064-020-0094-5.
- [1], [2], [3], [15]: Oguri et al. (EQUULEUS) J. Astronaut. Sci. 67(3), 2020, pp. 950-976; Dei Tos MS thesis (Politecnico di
  Milano, 2014); Dei Tos and Topputo, Adv. Space Res. 59(8), 2017, pp. 2117-2132; Park, Folkner, Williams and Boggs, DE440 and
  DE441, Astron. J. 161(3), 2021, p. 105.

Names the coordinator listed that are NOT cited in this paper: Leiva and Briozzo, Boudad-Howell-Davis 2019 (the NRHO in the
bicircular problem; only the 2022 AIAA analogs paper is cited), and Jorba and co-authors' 2018 and 2020 papers (only Rosales,
Jorba and Jorba-Cusco 2021 is cited, as reference 24). "Brown et al." is reference 20 (Hill four-body, not bicircular).

## 7. What this means for the project

### 7.1 Positive-control candidates (INFERRED)

1. The CR3BP numbers in tables 2 to 4 are controls on `core/cr3bp.py` and its stability and bifurcation code, independent of the
   Sun: Lyapunov family out-of-plane eigenvalue pair crosses +1 at P = 14.8319 d (halo bifurcation) and 18.7183 d (axial
   bifurcation); vertical family crosses at P = 19.2033 d. Tolerance: four decimals are printed, but the period units are
   days, so the nominal conversion constants (384,400 km, mu_E + mu_M = 403,503.2 km^3/s^2 as inferred in section 5.2) must be
   matched and a 1e-4 relative tolerance is a fair allowance. These are published, independent of our code.
2. The system periods 27.2846 d (2 pi) and 29.4906 d (2 pi/omega_S): a unit check on `core/qbcp.py` and `core/bcr4bp.py`
   (the Sun's synodic period in days with a 384,400 km Earth-Moon distance).
3. QBCP PO brackets: a Lyapunov PO of the QBCP at period `q P_sys` with p:q = 35:22 (P = 18.5370 d) must have its
   out-of-plane eigenvalue pair of the opposite stability type to the 27:17 orbit (P = 18.5682 d); the pair 2:1 and 121:61
   brackets the halo bifurcation (14.7453 d and 14.8672 d). A test needs a QBCP whose alpha_i match the authors'; this
   paper does not print them, so the test would have to be made with the Gimeno-Jorba/Jorba-Cusco-Farres coefficients
   already in `core/qbcp.py` and the result is a real test of how close those coefficients are to the authors'. The paper gives no tolerance; a mismatch would not by itself mean our model is wrong, since the authors' alpha_i come from a
   separate semi-analytic construction (ref 7) that is not printed here.
4. The QBCP equations (11) to (14) as printed: the placement of alpha_6 on the whole `U_Q*` (consistent with the #892
   correction), and the Sun at `-alpha_7 x-hat - alpha_8 y-hat`. A symbolic check of `core/qbcp.py`'s acceleration against
   eqs. (11) to (13) is possible without any numerics beyond the alpha_i, since the form is printed. INFERRED that the
   conventions (frame orientation, sign of alpha_7 and alpha_8) agree; check before relying on it.
5. The continuation statement that the periodic Lyapunov out-of-plane stability sequence centre, saddle, centre holds in the QBCP
   (qualitative control; also a check on the PO machinery).
6. The bifurcation drift direction (QBCP bifurcation periods below the CR3BP ones, by roughly 0.018, 0.18 and 0.20 d for the
   halo, axial and vertical bifurcations on the QPO estimates) is a quantitative check on a QBCP model's Sun strength. It needs the
   QPO machinery or a dense scan of resonant POs, so it is a large job.

The three most useful printed numbers: (a) the QBCP axial bifurcation bracket 18.5370 - 18.5682 d with QPO estimate
about 18.5418 d (against the CR3BP 18.7183 d), the strongest effect printed (about 1 percent); (b) the QBCP vertical-axial
bifurcation bracket 18.9185 - 19.0262 d, QPO estimate about 19.0029 d (CR3BP 19.2033 d); (c) the QBCP halo bifurcation
bracket 14.7453 - 14.8672 d (QPO about 14.8136 d; CR3BP 14.8319 d), together with `omega_S ~ 0.925196` and the inferred
system periods.

### 7.2 Overlap with #884

None directly, said plainly (section 5.4): no cycler, no resonant lunar-flyby orbit and no Earth-Moon transfer is treated. The
indirect relevance is as follows (INFERRED): (i) the p:q resonance with the system period is the structural condition for
a periodic orbit under periodic solar forcing, which is the condition for any exactly periodic orbit of a Sun-forced Earth-Moon model
(period commensurate with `2 pi/omega_S`), so the paper's statement that POs "exist at discontinuous period values
within the intermediary-fidelity models" is a published precedent for the project's resonance-gated search; (ii) the `eps`
homotopy and the pair of counterparts A and B at the two apsidal configurations (Sun and Moon in conjunction or in
opposition) correspond to the two phases a Sun-forced search must consider; (iii) the paper's references 18 (Gao et al.,
resonant periodic orbit web in the QBCP), 23 (Olikara, Gomez and Masdemont) and 24 (Rosales, Jorba and Jorba-Cusco, halo-like
tori around L2 in the BCP) are the published literature nearest to a Sun-forced periodic orbit search and should be checked by
the literature review for #884.
