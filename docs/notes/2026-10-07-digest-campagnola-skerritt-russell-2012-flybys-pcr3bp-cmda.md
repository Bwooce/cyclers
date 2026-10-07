# Digest: Campagnola, Skerritt & Russell 2012, "Flybys in the planar, circular, restricted, three-body problem" (#960 batch 37)

S. Campagnola (JPL), P. Skerritt (Caltech), R. P. Russell (UT Austin), "Flybys in the planar, circular, restricted,
three-body problem", Celest. Mech. Dyn. Astron. 113(3):343-368 (2012), doi 10.1007/s10569-012-9427-x
(Crossref, checked 2026-10-07: CMDA 113(3), pp.343-368, published 2012-06-22, authors Campagnola, Skerritt,
Russell). Received 30 Nov 2011, revised 24 Mar 2012, accepted 7 May 2012. Journal form of AAS 11-245.
- Supplied file `8246bb21-campagnola2012.pdf`, 26 pp. = journal pp.343-368 (PDF page n = p.342 + n). Springer
  typeset PDF with a good text layer. md5 ba29997f991772ace3ca31f2bfe9f9b5.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/campagnola-skerritt-russell-2012-flybys-planar-circular-restricted-three-body-problem-cmda-113-343-doi-10.1007-s10569-012-9427-x.pdf`
- **How I read it:** whole text from the text layer. On 150 dpi page images: pp.345-346 (Eqs 1-4, Fig. 1),
  p.347 (Fig. 2 and the alpha_Bmax sentence), p.348 (Eqs 7-10), p.350 (Eqs 11-13), pp.353-356 (Table 1,
  Figs 5-10, Table 2, Eq. 18), p.351 (Eqs 14-16, R2min, footnote 2), p.352 (Fig. 4, Eq. 17), p.357
  (Fig. 11, sec. 4.2), p.358 (Fig. 12, Eq. 19), p.359 (Figs 14-15, the 1.63 / 1.17 km/s baseline),
  pp.360-361 (Figs 16-17, sec. 6.1, Eqs 20-21, footnote 3), p.362 (Eqs 22-27), p.364 (Eqs 30-38), p.367
  (references). At 400 dpi, cropped and pixel-calibrated on the gridlines: Fig. 8a, Fig. 10 (Example 9),
  Fig. 12. Checks: `check_campagnola2012.py` -> `check_campagnola2012.out`.
- Wanted list: **row 22** (this paper). Also read for comparison: Campagnola & Russell 2009 AAS 09-227
  (held, p.6 on the image).

## 0. Verdict

**A methods paper. It defines the "Flyby map", a numerically integrated Poincare map of the planar CR3BP
from a section just after one apse to a section just before the next apse, both far (> 5 R_Hill) from
the moon. It uses the map to (1) put three-body flyby reach on the Tisserand-Poincare (T-P) graph,
(2) find two flyby families, Type I (direct) and Type II (retrograde), and (3) define CWIC ("Conics, When
I Can"): integrate only near the moon, use Kepler conics elsewhere. One worked example: a Callisto-to-
Europa leg with two Ganymede flybys (3:2 then 1:1), all outside the linked-conics domain, Europa orbit
insertion 1.02 km/s instead of 1.17 km/s.**

- **The error number #968 asked for is not in this paper.** No table compares the flyby model with the
  patched conic. There is no sweep in mass ratio, altitude or v_inf. The paper does not report error in
  v_inf, turn angle, periapsis or Tisserand parameter. The only comparison is the post-flyby semi-major
  axis a_B, in two figures, for Europa only:
  - Fig. 8 (p.355): one energy (T = 2.978935, v_inf about 2 km/s), from the 3:2 resonance, H >= 100 km.
  - Fig. 12 (p.358): the T-P graph envelope from 3:2 over a range of T.
  - The qualitative statement (p.355, image): "the linked-conics model underestimates a_B in case of direct
    flybys, and overestimates a_B in case of retrograde flybys." And (p.357): "Type I flybys (direct) are
    more efficient than Type II flybys (retrograde)."
  - Qualifier (from my Fig. 8a readings): the under/over sentence holds on the a-raising branch only. On the
    a-lowering branch the signs flip (Type I gives a lower minimum a_B than LC, Type II a higher one). The
    statement that holds on both branches is: reach of Type I > reach of LC > reach of Type II.
- **Per-moon summary (plain statement for #968):**

  | Moon | What the paper gives | Flyby map vs linked conics |
  |---|---|---|
  | Europa (mu 2.526636e-5) | Figs 5-9, 12; Table 2 Examples 1-2 | Fig. 8 at v_inf about 2 km/s, H = 100 km, from 3:2 (my readings, approximate): maximum a_B, Type I about +2.5 %, Type II about -1.8 % of the LC value; minimum a_B, Type I about -2.4 %, Type II about +2.6 %. As a fraction of the LC change in a: Type I about +26 % (raise) and +11 % (lower); Type II about -18 % and -12 %. Reading uncertainty about +-0.04e5 km (about +-4 % of delta-a). |
  | Ganymede (mu 7.803691e-5) | Fig. 10, Table 2 Examples 3-10, sec. 6 | No LC overlay. At T = 2.9898 (v_inf 1.1 km/s) the example lies wholly outside the LC domain (p.360: "in the linked-conics domain there are no paths forward"). So the LC model has no answer to compare. |
  | Callisto | only the incoming v_inf 1.9 km/s (p.361) | none |
  | Titan, Enceladus | absent | none |

- **The nearest published cost-error number is in the held Endgame Part B (AAS 09-227 p.6, image).**
  "the cost of the VILM endgames can be off up to +-5% when compared to the more accurate CR3BP solutions"
  (Europa). Footnote: "consistent with the +-10% difference observed during the design of the Cassini tour
  when comparing dv costs in the linked conics model with more accurate models (personal communication from
  Nathan Strange)". The worked numbers: direct (long-transfer) endgame 147 m/s in the CR3BP against 154 m/s
  VILM; retrograde (short-transfer) 165 m/s against 155 m/s. Our arithmetic (`.out` sec. K): VILM is +4.8 %
  of the CR3BP cost for the direct case and -6.1 % (-6.5 % relative to VILM) for the retrograde case. So
  the paper's own example exceeds its "up to +-5%". **For #977's VILT cost axis: use about 6.5 % (the
  computed worst case of the published example, Europa) as the sourced band; the quoted +-5 % understates
  it; the +-10 % Cassini figure is unpublished and second-hand.**
- **What it gives the project:**
  - the Flyby map and CWIC definitions with the full coordinate chart (Appendix B), enough to implement;
  - a ready positive control for any future CWIC or flyby-map code: Table 2 initial conditions plus the
    Fig. 8 endpoints and the Fig. 10 shape (section 3);
  - EOI and v_inf numbers for a Ganymede endgame that #977 can use as goldens (section 4).
- **Project code equivalent: none.** `grep -ril "cwic\|flyby_map\|flyby map" src` finds nothing.
  Nearest code: `genome/keplerian_map.py` (Ross-Scheeres Keplerian map; this paper and Lantoine et al.
  2011 say it is accurate only at very low energy), `genome/composed_moon_map.py` (patched Keplerian moon
  map), `search/tisserand.py` (linked-conic Tisserand graph), `search/periapse_map.py` (integrated periapse
  map). The CWIC model is also what Anderson, Campagnola & Lantoine 2016 (held, digested) use for their
  resonant orbits.
- **Error found in project code (PROPOSAL to fix):** the `genome/keplerian_map.py` docstring cites RS07 as
  DOI `10.1137/06065195X` and GR09 as pp.436-443. Crossref (2026-10-07): `10.1137/06065195X` does not
  resolve; Ross & Scheeres SIADS 6(3):576-596 is **10.1137/060663374** (as this paper prints); Grover & Ross
  JGCD 32(2) is **pp.437-444** (10.2514/1.38320).
- **PROPOSALS (code; none for the catalogue):**
  1. #968 control route: a CWIC flyby map in `genome/` (integrate the pcr3bp only inside R2 > 5 R_Hill
     sections, Kepler outside). Positive control: reproduce Fig. 8a's six endpoints and Table 2 Example 1/2
     closest approaches at Europa; then Fig. 10 at Ganymede. Our output would be a check, never a golden.
  2. #977: carry a model-error band on VILT costs: about 6.5 % (Europa; computed from the AAS 09-227 p.6
     example, whose text says "up to +-5%") as the sourced figure; the +-10 % Cassini figure only as a
     labelled second-hand value. Direction matters: VILM was optimistic for the retrograde endgame and
     pessimistic for the direct one.
  3. Fix the two citation slips in `keplerian_map.py` above.
- **Catalogue implication:** none. No cycler, no periodic orbit with published closure data.

## 1. Content

- **Sec. 1 (pp.343-345).** Linked conics fail at low v_inf. The Keplerian map (Ross & Scheeres 2007) is
  accurate only at very low energy (Lantoine et al. 2011). The T-P graph (Campagnola & Russell 2010b) shows
  reachable sets but cannot predict one flyby. Motivation: in a real-ephemeris model direct gravity assists
  were sometimes better than retrograde ones (Campagnola & Russell 2010b). Says the method was used for
  NASA Europa orbiter and ESA JGO endgames.
- **Sec. 2.1 (pp.345-347), linked conics (image).** Units M, a_m, sqrt(a_m^3/GM).
  - Eq. 1: sin(delta/2) = Gm / (Gm + v_inf^2 (r_m + H)), H > H_min.
  - Eq. 3: alpha_Bmin = max(0, alpha_A - delta), alpha_Bmax = min(pi, alpha_A + delta); alpha is the pump
    angle between the moon's velocity and v_inf.
  - Eq. 4 (= Eq. 24, 27): a_B = 1/(1 - 2 v_inf cos alpha_B - v_inf^2);
    e_B = sqrt(1 + (v_inf^2 + 2 v_inf cos alpha - 1)(1 + v_inf cos alpha)^2). Algebra checked (`.out` sec. H).
  - Fig. 1: direct (on the incoming leg) and retrograde (outgoing leg) flybys give the same alpha_B, so the
    Tisserand graph cannot tell them apart.
- **Sec. 2.2 (pp.347-348, image).** pcr3bp, barycentric units. Eq. 7 Jacobi constant includes the
  constant + (1 - mu) mu. Eq. 8: J ~ T = (1 - mu)/a + 2 sqrt(a(1 - e^2)) = 2/(r_a + r_p) +
  2 sqrt(2 r_a r_p/(r_a + r_p)). Eq. 9: T = 3 - v_inf^2. Eq. 10: R_Hill = (mu/(3(1 - mu)))^(1/3).
- **Sec. 3 (pp.349-352), the Flyby map.**
  - Coordinates (a, T, lambda, f): lambda is the rotating-frame longitude of the osculating pericentre
    (a > 1) or apocentre (a < 1); f-bar = 0 (a > 1) or pi (a < 1) (Eq. 11). Osculating elements are about the
    primary's centre with GM = 1 - mu (Appendix B "Alternative formulations").
  - Section A: a != 1, f = -pi + eps (a > 1) or eps (a < 1), R2 > R2min; section B: f = pi - eps or
    2 pi - eps (Eqs 13-16). R2min = 5 R_Hill "to ensure the validity of the Tisserand condition" (p.351).
  - F: A -> B, first crossing. T_A ~ T_B, so T is a parameter and the map is 2-D: (a_A, lambda_A) ->
    (a_B, lambda_B). For mu -> 0 it is the identity.
  - Symmetry (Eq. 17, image): F^-1 = Theta o F o Theta with Theta(a, T, lambda) = (a, T, -lambda).
- **Sec. 4 (pp.352-358).** From the 3:2 resonance (a_A = (3/2)^(2/3)) at Europa:
  - Fig. 5: T = J_L1 = 3.0036678286; only Type I; a_B(lambda_A) has one max (about 8.92e5 km) and one min
    (about 8.65e5 km) (approximate).
  - Fig. 6: T-bar = 2.988, first collision orbit at lambda-bar about 0 deg.
  - Fig. 7: T = 2.978935 (v_inf = 2 km/s "where defined"); Type II appears between collision orbits at
    lambda about +-3 deg.
  - Fig. 8: same with H > 100 km (Eq. 18) and the LC curves overlaid (numbers in section 2).
  - Fig. 10 (Ganymede, T_A = 2.9898, v_inf 1.1 km/s): new multi-loop solutions near lambda_A about -1.4 deg
    and 1.6 deg; shadow unstable bifurcations of the DRO family; R2min = 0.5 here (p.356).
  - Fig. 12: T-P graph envelopes. Type I reaches perijove up to about 7.0e5 km (above Europa's orbit, outside
    the LC domain); Type II stays just below 6.71e5 km.
- **Sec. 5 (pp.358-359), CWIC.** A(i) -F-> B(i) -> A(i+1) = (a_B, T_B, lambda_B - 2 pi sqrt(a_B^3)) (Eq. 19,
  image). Where R2(t) > R2min for a whole revolution, F is replaced by the identity. Fig. 14: the general CWIC
  can add impulses and need not be the pcr3bp.
- **Sec. 6 (pp.359-361), example.** Callisto (v_inf 1.9 km/s) -> Ganymede flyby (backwards from 3:2) ->
  3:2 -> 1:1 (identity map) -> Ganymede Type I flyby, lambda_A^(2) = 1.4784 deg (Example 8) -> a_B^(2) = 0.8076
  -> Europa approach, v_inf 1.38 km/s, EOI 1.02 km/s into a 50 km circular orbit. Literature baseline
  v_inf 1.63 km/s, EOI 1.17 km/s. Footnote 3: the target (a, T) is found on the T-P graph in the Jupiter-
  Europa pcr3bp, not by linked conics. Used for the lander option of the 2011 Europa Habitability Mission
  Study.

## 2. Checks (our arithmetic; `check_campagnola2012.out`)

| Item | Printed (page, image) | Ours | Result |
|---|---|---|---|
| 3:2 a_A | 1.31037 (Table 2) | (3/2)^(2/3) = 1.310371 | agrees |
| v_inf at Europa from T = 2.978935 | 2 km/s (p.354) | 1.9943 km/s (GM_J 1.26686534e8, a 670974 km) | 0.3 % low. T for exactly 2 km/s is 2.978815; gap 1.2e-4 in T. Unit choice (GM vs G(M+m)) does not explain it. Unresolved; small. |
| v_inf at Ganymede from T = 2.9898 | 1.1 km/s (p.356) | 1.0988 km/s | agrees |
| J at L1, mu_Eu | 3.0036678286 (p.353) | 3.0036667439 with Eq. 7's (1-mu)mu term; 3.0036414782 without | 1.1e-6 gap with the paper's convention. Unresolved (perhaps a different mu digit). Note: our `cr3bp.jacobi_constant` omits the (1-mu)mu term, so the paper's J values are larger by (1-mu)mu = 2.53e-5 at Europa. |
| LC max a_B, Fig. 8a | read 9.75e5 km | 974,541 km (Eqs 1, 3, 4; H = 100 km, R_CA 1660.8 km) | agrees: my reading is calibrated |
| LC min a_B, Fig. 8a | read 7.26e5 km | 722,205 km | agrees within reading error |
| LC envelope, Fig. 12 | black circles at r_a 8.3-8.9e5, r_p 5.85-6.5e5 km; upper branch on r_p = 6.71e5 | min side r_a 8.32-8.86e5, r_p 5.77-6.40e5 for v_inf 1.6-2.8 km/s; max side r_p 6.71e5 | agrees |
| Fig. 10 Example 9 peak vs a_A^(0) (Eq. 21 needs a_B > a_A^(0)) | figure | peak read about 1.563e6 km (calibrated on the 1.5e6/1.6e6 gridlines; a_A line reads 1.40e6 vs computed 1,402,539 km); a_A^(0) = 1,559,696 km | consistent, but the margin (about 3,000 km) is inside my reading error. It also supports the dropped-digit reading of "155980 km". |
| 0.8076 a_Ga | 864396 km (p.361) | 864,404 km | agrees (rounding of 0.8076) |
| 1.4572 a_Ga | "155980 km" (p.361) | 1,559,696 km | **print slip: one digit dropped** (about 1559700 km). Fig. 17 puts A(0) at apojove about 2.03e6, perijove about 1.09e6 km; ours 2,026,742 and 1,092,649 km. |
| Callisto v_inf from A(0) (a 1.4572, T_Ga 2.9897) | 1.9 km/s | 1.942 km/s (Tisserand w.r.t. Callisto) | agrees to the printed 2 digits |
| Europa v_inf from B(2) (a 0.8076, T_Ga 2.9897) | 1.38 km/s (Europa approach) | 1.476 km/s by linked-conic Tisserand; B(2) perijove 670,569 km, almost tangent to Europa's orbit. Rounding test (`.out` sec. J): T 2.9897-2.9898 and a 0.80755-0.80765 give 1.464-1.479 km/s | **does not agree, about 0.09-0.10 km/s; not a rounding effect.** Not a contradiction: footnote 3 says the Europa approach was found in the Jupiter-Europa pcr3bp, not by linked conics. It is consistent with a three-body approach doing better than the conic estimate, but the paper does not say this. Unresolved. |
| EOI, v_inf 1.63 km/s, 50 km circular | 1.17 km/s (p.359) | 1.1655 km/s (GM_Eu 3202.739, R 1560.8) | agrees |
| EOI, v_inf 1.38 km/s | 1.02 km/s (p.361) | 1.0150 km/s | agrees |

**Print slips (image-confirmed):**
- p.347: "maximum jump occurs for alpha_B = alpha_Bmax (maximum a_B) and alpha_B = alpha_Bmin (minimum
  a_B)". By Eq. 4 a_B falls as alpha_B rises, so alpha_Bmax gives the MINIMUM a_B (`.out` sec. D). The two
  labels are swapped.
- p.353: "(a_B, T, lambda_B) = F((3/2)^(3/2), T, lambda_A)"; p.361: "(3 : 2)^(3/2)". Both should be 2/3
  (Table 2 and the line above on p.353 use 2/3).
- p.361: "a few zeros of Eq. 21 are found" means Eq. 20. "The free design parameters are lambda_B^(0) and
  lambda_A^(0)" should be lambda_A^(2) (the next sentence chooses lambda_A^(2)).
- p.361: T = 2.9897 "the same used to compute Fig. 10"; Fig. 10 and Table 2 use 2.9898.
- p.348 Eq. 7 prints (Xdot + Ydot) for (Xdot^2 + Ydot^2).
- p.364 Eq. 35 sets h = (T - (1 - mu)/a)/2, which is sqrt(a(1 - e^2)), not h = sqrt((1 - mu) a(1 - e^2));
  Eq. 36 then uses h^2/(a(1 - mu)). Inconsistent at O(mu); harmless for small mu but matters if coded
  literally.
- p.367: Petit & Henon "Icarus 555, 536-555 (1986)"; Crossref: Icarus 66:536-555,
  10.1016/0019-1035(86)90089-8.

## 3. Positive-control data for a future CWIC / flyby-map implementation (PROPOSAL)

Published inputs (Table 2, p.356, image), all with a_A = 1.31037:

| Example | T_A | lambda_A | mu | Use |
|---|---|---|---|---|
| 1 | 2.9789 | -3.39 deg | mu_Eu | Type I, Fig. 8-9, ends at the H = 100 km edge, a_B about 9.99e5 km (Fig. 8a, approx.) |
| 2 | 2.9789 | 2.76 deg | mu_Eu | Type II, a_B about 9.57e5 km (approx.) |
| 3, 4 | 2.9898 | -1.3624, -1.3667 deg | mu_Ga | multi-loop, Fig. 11 |
| 5, 6, 7 | 2.9898 | 1.5893, 1.6411, 1.6568 deg | mu_Ga | multi-loop |
| 8 | 2.9898 | 1.4784 deg | mu_Ga | sec. 6 flyby, a_B = 0.8076 (p.361) |
| 9 | 2.9898 | -2.95 deg | mu_Ga | a_B > a_A^(0) = 1.4572; at the Fig. 10 maximum, about 1.563e6 km (approx.) |
| 10 (a_B, T_B, lambda_B) | 2.9898 | 2.95 deg | mu_Ga | the inverse of 9 by Eq. 17 |

Table 1 constants (p.353): a_Eu = 670974.45148 km, mu_Eu = 2.526636E-5; a_Ga = 1070337.37782 km,
mu_Ga = 7.803691E-5. R2min = 5 R_Hill (p.351), except Fig. 10 R2min = 0.5 (p.356).

**Golden candidates (published values only):** T(J_L1, Europa) = 3.0036678286 (p.353, paper's Eq. 7
convention, which adds (1 - mu) mu to the project's `cr3bp.jacobi_constant`; we reproduce it only to 1.1e-6,
so use a tolerance of 2e-6 or resolve the gap first); Table 2 rows; a_B^(2) = 0.8076 for Example 8 (p.361); EOI 1.02 km/s at v_inf 1.38 km/s and
1.17 km/s at 1.63 km/s, 50 km circular Europa orbit (pp.359, 361); the Endgame Part B costs 147/154 and
165/155 m/s (AAS 09-227 p.6; its quoted +-5 % band is text, the example gives up to 6.5 %). The Fig. 8 and Fig. 12 numbers are figure readings and
must not be used as goldens.

## 4. Comparison with held papers

- **Campagnola & Russell 2009, AAS 09-224 and 09-227 (held; journal forms JGCD 33(2):463-475 and 476-486,
  10.2514/1.44258 and 10.2514/1.44290, not held as separate files).** Part B introduced the T-P graph and
  the direct (long-transfer, L2 side) vs retrograde (short-transfer, L1 side) endgames: CR3BP 147 vs 165 m/s
  for VILM estimates 154 vs 155 m/s. This paper explains that result: Type I (direct) flybys reach farther
  than linked conics predict, Type II less far.
- **Anderson, Campagnola & Lantoine 2016 (held, digested 2026-10-06).** Uses this paper's CWIC flyby map to
  find single-moon unstable resonant orbits (Jupiter-Europa 3:4 ... 7:8; Jupiter-Ganymede C = 2.99). This
  paper is the definition; the 2016 paper is the closure application. For #968, the 2016 orbits are the
  published periodic controls; this paper supplies the one-flyby control.
- **Campagnola, Strange & Russell 2010 and Lantukh, Russell & Campagnola 2015 (held, digested)** are the
  conic VILT side; they carry no three-body error number either.

## 5. Citation mining (selected; all references checked)

| Cited work | Held? | Wanted list |
|---|---|---|
| Anderson 2005, PhD thesis, CU Boulder | not held | not listed |
| Anderson & Lo 2010 JGCD 33(6); 2011 JAS 58(2) | held | - |
| Anderson & Lo 2012, AAS 12-136 (heteroclinic connections of unstable resonant orbits) | not held | not listed; relevant to X1 resonant legs, low priority |
| Belbruno, Topputo & Gidea 2008, Adv. Space Res. 42:1330 | not held (Belbruno 2004 book and Topputo-Belbruno 2015 held) | not listed |
| Campagnola & Lo 2008, PAMM 7 | not held | not listed |
| Campagnola & Russell 2010a, b, JGCD 33(2) | held as AAS 09-224, 09-227 | - |
| Carusi et al. 1982, 1990, 1995; Greenberg et al. 1988 | not held | not listed; Carusi, Kresak & Valsecchi 1995 (Tisserand conservation at close encounters, EM&P 68:71) would give an independent T-vs-J error source; propose as a low-priority row |
| Gawlik et al. 2009 (two papers) | not held | not listed |
| Grover & Ross 2009 JGCD 32(2) | held | - |
| Heaton et al. 2002 JSR 39(1) | not held (Heaton & Longuski 2003 held) | not listed |
| Howell, Davis & Haapala 2012 MPE | held | - |
| Howell, Marchand & Lo 2001 JAS 49(4) | not held | not listed |
| Johannesen & D'Amario 1999 | not held (only Lam, Johannesen et al. 2008 held) | not listed; source of the 3:4-5:6 Europa endgame |
| Kloster, Petropoulos & Longuski 2010 Acta 68:931 | not held (Lynam, Kloster & Longuski 2011 held) | not listed |
| Koon et al. 2000 Chaos 10 | held | - |
| Koon et al. 2001, CMDA 81:63-73 ("Low energy transfer to the Moon") | **not held** (the held `koon-lo-marsden-ross-2001-...-CMDA-81` is a different paper in the same volume) | not listed |
| Koon et al. 2008 book | not held | not listed |
| Labunsky, Papkov & Sukhanov 1998 | not held | row 52 |
| Lantoine & Russell 2012 "to appear" = JAS 58:335 (2011) | held | - |
| **Lantoine, Russell & Campagnola 2011, Acta 68(7-8):1361-1378, 10.1016/j.actaastro.2010.09.021** | **not held** | **not listed. Propose a Tier B row:** it is the quantitative source of the Keplerian-map error at Ganymede, and our `keplerian_map.py` has no validity bound from it. |
| Malyshkin & Tremaine 1999; Petrosky & Broucke 1988 | not held | not listed |
| Perozzi 1986 (Titan encounter) | not held | not listed |
| Petit & Henon 1986, Icarus 66:536 (printed "555") | not held | not listed |
| Ross & Lo 2003 (multi-moon orbiter) | held (`ross-koon-lo-marsden-2003-...AAS-03-143`) | - |
| Ross & Scheeres 2007 SIADS 6(3) | held | - |
| Ross, Jerg & Junge 2009 CNSNS 14 | not held | not listed |
| Russell 2006 JAS 54(2) (Europa periodic orbits); Russell & Lam 2007 JGCD 30(2); Russell & Lara 2009 Acta 65 | not held | not listed |
| Strange & Longuski 2002 JSR 39(1) | held | - |
| Strange, Russell & Buffington 2007 (V-infinity globe) | held | - |
| Sweetser 1993 | not held | row 50 |
| Uphoff, Roberts & Friedman 1976 | not held | row 41 |
| Vaquero & Howell 2011 AAS 11-428 | not held (Vaquero 2013 thesis held) | not listed (row 18 is the different Vaquero & Howell 2014 Acta paper) |
| Vasile & Campagnola 2009 JBIS | held | - |
| Villac & Scheeres 2003, 2004 | not held | not listed |
| Woolley & Scheeres 2009 | not held | row 52 (Woolley 2010 thesis) |

Only new row proposed: Lantoine, Russell & Campagnola 2011 (Acta). Optional: Carusi, Kresak & Valsecchi 1995.

*Filed as `cyclers_pdf/papers/campagnola-skerritt-russell-2012-flybys-planar-circular-restricted-three-body-problem-cmda-113-343-doi-10.1007-s10569-012-9427-x.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the pre-batch-37 numbering; the list was renumbered in batch 37.*
