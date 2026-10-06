# Digest: Breakwell & Perko 1965, "Matched Asymptotic Expansions, Patched Conics and the Computation of Interplanetary Trajectories" (#960 batch 29)

John V. Breakwell (Stanford University) and Lawrence M. Perko (Lockheed Missiles & Space Company, Palo Alto).
AIAA Paper No. 65-689, AIAA/ION Astrodynamics Specialist Conference, Monterey, 16-17 September 1965.
doi 10.2514/6.1965-689 (printed on the page margin of the scan). The abstract page lists the second affiliation as "Lockheed Palo Alto Research Laboratory".
- Filed as `cyclers_pdf/papers/breakwell-perko-1965-matched-asymptotic-expansions-patched-conics-interplanetary-trajectories-AIAA-65-689-doi-10.2514-6.1965-689.pdf`. Original upload, 12 sheets, md5 76d86791770ba12fa8381adf40d6a226.
  It is an AIAA reprint scan (KTH library stamp, 2015) of a typescript. Sheets 2-11 are two printed pages each; printed pages run 1-19, then references, then Figs. 1-5.
- Was wanted-list rank 20; removed in batch 29. (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.)

**How I read it.** The text layer is poor, so I rendered every sheet (130 dpi for the first pass; 220 dpi, cropped per printed page, for
sheets 5-12) and read the images. All formulas below are transcribed from the images (items (1)-(4'), the hyperbola relations, the
flyby map, and the figure axes). Numbers in Figs. 1, 3, 4 are graph readings and are approximate. Where an equation was
checkable I checked it numerically (section 3).

## 0. Verdict

**H5 background: this is the first-order matching result that Perko 1974 / Breakwell-Perko 1974 / Guillaume 1975 build on.**
It states the correction a perturbation makes to a patched-conic (Keplerian arc + planet-centred hyperbola) trajectory, to O(lambda^2)
error, where lambda = m_planet/m_sun.
- **What the patched-conic correction is:** the hyperbola at the planet is not the one a massless-planet patched conic gives.
  The sun and other planets shift it by "gross biases" (constant, computed once by quadrature or by calibration against one accurate run) in
  (1) the speed at infinity, (2) the asymptote direction, (3) the asymptote offset Delta, and (4) the closest-approach time, which also
  carries a "local" logarithmic term in ln(epsilon) of the hyperbola eccentricity.
- **Numerical content is small.** One Earth-Moon plot (Fig. 1), a comparison curve against Lagerstrom-Kevorkian (Fig. 2), an Earth-Mars
  plot of the biases (Figs. 3, 4), and a flow chart for an Earth-Venus-Mars flyby (Fig. 5). No data table. Nothing to key into the catalogue.
- **Use:** `X2`/`X3`/`R4` asymptotics in the wanted-list note. It gives the first-order formulas in the same notation the held 1974 digest
  quotes (Delta, B, B', epsilon). Useful as a test spec: a flyby at small mass ratio, fit the osculating hyperbola at closest approach,
  compare Delta, direction, v_inf and t_p with eqs. (3). The 1974 digest already proposed that test; this paper licenses its first-order half.
- **Catalogue implication (PROPOSAL only):** none. It is a method paper.
- **Venue caution.** The held 1974 digest (`2026-10-04-digest-breakwell-perko-1974-second-order-matching.md`, line 10) cites "Breakwell & Perko 1965
  (Proc. 16th IAF Congress; Progress in Astronautics 17, 1966)". This AIAA preprint (Monterey, September 1965) is the paper the wanted
  list means (its DOI is confirmed). I did not compare the two venues' texts; they may be the same work in different forms. I have not seen the IAF form.

## 1. Method (printed pp. 4-9)

Setup (p. 4). Spaceship of negligible mass in the field of the sun (m0) and planets m_i, i = 1..N. lambda = m_1/m_0 << 1.
Heliocentric equation: r'' = -mu0 r / r^3 + f(r, t), with f = -sum_i mu_i [ (r - r_i)/|r - r_i|^3 + r_i / r_i^3 ], mu_i = G m_i.
Perturbation series r = r^(0) + rho^(1) + rho^(2) + ..., with r^(0) the unperturbed conic. The only restriction is that the approach velocity
relative to the planet is O(1) compared with the planet's orbital speed. Planet orbits may be eccentric and inclined.

First order (p. 5), with the transition matrix Phi(t,t') (blocks Phi_rr, Phi_rv, Phi_vr, Phi_vv, refs Battin and Danby):
[rho^(1); rho^(1)'](t) = Phi(t,t0) [rho^(1); rho^(1)'](t0) + integral from t0 to t of Phi(t,t') [0; f(r^(0)(t'), t')] dt'.
Second order: rho^(2) is the same integral with source (df/dr) rho^(1) + ... Only first order and the singular second-order terms are carried.

Near the arrival planet (p. 6). x-axis along the relative velocity v_1 at arrival time t1, tau = t1 - t:
x^(0) = -v1 tau + (mu0 v1 / (6 r1^3)) (1 - 3 alpha0^2) tau^3 + O(lambda tau^2, tau^4),
y^(0) = -(mu0 v1 / (6 r1^3)) alpha0 beta0 tau^3 + ..., z^(0) = 0 + ...
(alpha0, beta0, 0) are the sun's direction cosines at t1 in the planet frame.
Singular behaviour (p. 7): Phi_rv = (t - t') I_3 + O((t-t')^3), Phi_vv = I_3 + O((t-t')^2). The transported source has a 1/tau'^2 part, so
rho^(1)(t) = i (mu1 / v1^2) ln(tau0 / tau) + rho_b^(1)(t), tau0 = t1 - t0, with rho_b^(1) and its rate bounded.
rho_b^(1)(t) = Phi_rr rho^(1)(t0) + Phi_rv rho^(1)'(t0) - i (mu1/v1^2)(1 - tau/tau0) + integral of b(t,t') dt'.
Second order singular part (p. 8):
rho^(2)(t) = i (mu1^2 / (v1^5 tau)) [ln(tau0/tau) - 3/2] + (mu1 / (2 v1^3 tau)) { 3 [rho_b^(1)(t1) . i] - rho_b^(1)(t1) } + lesser terms.
Footnote p. 8, citing Perko's Stanford thesis: remainder after n terms |R_n| = O[(lambda^2)^n ... |ln tau / tau|^(n-1)] (the printed typescript is hard to read
here; I could not resolve the exact form), so the three-term composite is O(lambda^2) for tau = O(lambda^(1/2)).

Outer (sun-centred) expansion near the planet, eq. (1) (p. 8): x(t), y(t), z(t) as series in tau with terms
x = -v1 tau + [ (mu1/v1^2) ln(tau0/tau) + x_b^(1)(t1) ] + { (mu0 v1/(6 r1^3))(1 - 3 alpha0^2) tau^3 - x_b'^(1)(t1) tau
+ (mu1^2/(v1^5 tau)) [ln(tau0/tau) - 3/2] + (mu1/(v1^3 tau)) x_b^(1)(t1) } + ...;
y = y_b^(1)(t1) - { (mu0 v1/(2 r1^3)) alpha0 beta0 tau^3 + y_b'^(1)(t1) tau + (mu1/(2 v1^3 tau)) y_b^(1)(t1) } + ...; z likewise without the sun term.
Inverted in x (eq. 1', p. 9). Valid for tau = O(lambda^(1/2)), ascending powers of lambda^(1/2) through O(lambda^(3/2)).

Inner (planet-centred) expansion (pp. 9-12). The unperturbed hyperbola has |a| = mu1 / v_inf^2, eccentricity epsilon and pericentre time t_p:
tau* = t_p - t = (mu1 / v_inf^3)(epsilon sinh E - E), r* = (mu1 / v_inf^2)(epsilon cosh E - 1), tanh(E/2) = sqrt((epsilon-1)/(epsilon+1)) tan(theta/2),
deflection delta = 2 arcsin(1/epsilon), and with Delta = (mu1 / v_inf^2) sqrt(epsilon^2 - 1):
x* = -(mu1/v_inf^2)[ (epsilon/2) e^E - 1 + (1/epsilon - epsilon/2) e^(-E) ], y* = Delta (1 - e^(-E)/epsilon).
The sun's tidal term on the planet-centred orbit gives rho* = (mu0 v_inf / (6 r1^3)) { i* - 3 [r1 . i*] r1 / r1^2 } tau*^3 + O(lambda^2).
Large-E forms, eq. (2) (p. 11):
x* = -v_inf tau* + (mu1/v_inf^2)[1 - ln(2 v_inf^3 tau* / (mu1 epsilon))]
   + { -(mu1^2/(v_inf^5 tau*)) [ln(2 v_inf^3 tau*/(mu1 epsilon)) + 1/2] + (mu0 v_inf/(6 r1^3))(1 - 3 alpha0*^2) tau*^3 } + O(lambda^2),
y* = Delta (1 - mu1/(2 v_inf^3 tau*)) - (mu0 v_inf/(2 r1^3)) alpha0* beta0* tau*^3 + O(lambda^2),
z* = -(mu0 v_inf/(2 r1^3)) alpha0* gamma0* tau*^3 + O(lambda^2).

## 2. The matching formulas (printed pp. 12-14)

Equating inner (2') and outer (1') terms gives, to O(lambda^2), eq. (3) (p. 12):
- v_inf = v1 + x_b'^(1)(t1)  (speed at infinity),
- i* = i + j y_b'^(1)(t1)/v1 + k z_b'^(1)(t1)/v1  (asymptote direction),
- j* Delta = j y_b^(1)(t1) + k z_b^(1)(t1)  (offset and orientation of the plane; Delta is the miss distance of the asymptote),
- v_inf (t1 - t_p) = x_b^(1)(t1) + (mu1/v_inf^2)[ ln(2 v_inf^3 tau0 / (mu1 epsilon)) - 1 ]  (time of closest approach),
- epsilon = [1 + (v_inf^2 Delta / mu1)^2]^(1/2).
**The time correction is a gross bias x_b^(1)(t1) plus a local bias (mu1/v_inf^3)[ln(.) - 1] that depends on epsilon.**
Velocity match is only to O(lambda^(3/2)) in the raw expansions; the paper argues no O(lambda^2) term linear in tau can arise, so (3) holds to O(lambda^2).
Compact form, eq. (3') (p. 13): [ (mu1/v1^3)[2 - ln(2 v1^3 (t1-t0)/(mu1 epsilon))] + t1 - t_p ] v1 + Delta ; v_inf - v1 + v1 mu1/(v1^3 (t1-t0)) =
Phi(t1,t0) [r(t0) - r^(0)(t0); r'(t0) - r^(0)'(t0)] + integral of [b(t1,t); b'(t1,t)] dt + O(lambda^2),
with bounded b(t1,t) = Phi_rv(t1,t) f[r^(0)(t),t] - mu1 v1 / (v1^3 (t1-t)) and b' = Phi_vv f - mu1 v1 / (v1^3 (t1-t)^2).
The same sign structure holds for the departure planet (p. 15, with ln(...) - 2 and t0 - t1).

**Two-planet (Earth to Venus) form, eq. (4)-(4') (pp. 14-18).** Consider a Keplerian arc leaving a massless planet 1 at t1 and arriving at a massless planet 2 at
t2, with t0 an intermediate time. The t0 dependence cancels (using Phi(tA,tB) Phi(tB,tC) = Phi(tA,tC), p. 16). The "gross biases" (p. 17), i = 1, 2:
B_i = integral from t1 to t2 of b_i(t_i, t) dt + (mu_i v_i / v_i^3)[ ln(2 v_i^3 (t2 - t1) / mu_i) - 2 ],
B'_1 = integral of b'_1 dt + mu1 v1 / (v1^3 (t2 - t1)), B'_2 = integral of b'_2 dt - mu2 v2 / (v2^3 (t2 - t1)),
and for the other planets i = 3..N: B_i = integral of Phi_rv(t2,t) f_i dt, B'_i = integral of Phi_vv(t2,t) f_i dt.
Abbreviated form, eq. (4') (p. 18): the arrival hyperbola constants (t2 - t_p2 + (mu2/v2^3) ln eps2) v2 + Delta2 and v_inf2 - v2 equal
[B_2; B'_2] + sum over i = 3..N of [B_i; B'_i] + Phi(t2,t1) [B_1; B'_1] + [ (t1 - t_p1 - (mu1/v1^3) ln eps1) v1 + Delta1 ; v_inf1 - v1 ] + O(lambda^2).
The gross biases "can either be computed as definite integrals of bounded functions ... or determined by calibration against an accurate
numerical program. In fact, one accurately computed trajectory will suffice" (p. 19). The sum of the biases through the encounter is then reused.

**Flyby map at the intermediate planet (p. 19).** With k = v_inf2^2 Delta2 / mu2 (so 1 + k^2 = epsilon2^2), arrival (-) to departure (+):
v_inf2(+) = [ (k^2 - 1) v_inf2(-) - (2 v_inf2^3 / mu2) Delta2(-) ] / (1 + k^2),
Delta2(+) = [ (k^2 - 1) Delta2(-) + (2 v_inf2 Delta2^2 / mu2) v_inf2(-) ] / (1 + k^2)
(read on the image; in the second line the right-hand vector carries a small accent mark that I could not resolve and read as v_inf2(-)).

## 3. Numbers and checks

- **COMPUTED check of eq. (2) against the exact hyperbola** (mu1 = 0.01, v_inf = 1, epsilon = 1.8, E = 6, 8, 10). The exact x* from the
  E-parametric form minus the printed eq. (2) through the 1/tau* term: 1.1e-6, 3.8e-8, 1.1e-9. Dropping the 1/tau* term the error is 1.8e-4, 3.2e-5, 5.3e-6.
  So the 1/tau* coefficient -(mu1^2/(v^5 tau*)) [ln(.) + 1/2] and the log term are correct as transcribed. The y* form agrees likewise (E = 10: 0.0149662521 exact vs
  0.0149662519 from eq. 2). The error falls as 1/tau*^2 as expected.
- **COMPUTED check of the flyby map.** With v_inf2(-) along x and Delta2(-) along y, the x- and y- components give cos(delta) = (k^2-1)/(k^2+1) and
  sin(delta) = 2k/(1+k^2). With epsilon = 1.8: 0.3827160494 and 0.9238660214, equal to cos and sin of delta = 2 arcsin(1/epsilon) to 1e-15.
  So v_inf2(+) is v_inf2(-) rotated by the hyperbolic deflection angle, as it must be.
- **Not checked:** eqs. (1), (1'), (3), (4), (4') were not re-derived. The singular second-order coefficient in rho^(2) was not re-derived.
- **Figures (graph readings, approximate).**
  - Fig. 1 (Earth-Moon, coplanar, perigee distance 4050 miles): Delta peaks near +1.9 (units 10^4 mi) at about 22 days and dips to about -0.7 near 10 days; the
    time bias t_p - t1 (units 0.2 day) dips to about -1.2 near 6 days. The curves cross near 14 days. Text: "for transfer times greater than 14 days, the moon has the effect
    of causing the particle to arrive late".
  - Fig. 2: comparison with Lagerstrom and Kevorkian for h0 = 0.01, 0.077, 0.18 (0.18 is the 4050-mile perigee). "The comparison is seen to be very good for h0 sufficiently small."
    Their results assumed h0 <= O(lambda^(1/2)).
  - Fig. 3 (Mars' effect on Earth-Mars trajectories, axis: transfer angle 70-250 deg, t2 - t1 50-550 days): Delta_2 peaks near +2.5 (10^3 mi) at about 275 days; t_p2 - t2 bottoms near -2 (0.1 day) at about
    180 days.
  - Fig. 4 (Earth's effect): the scale is 10^-5 mi for Delta_2 and "t_p2 - t2 (days)", with the curves reaching about -8 (10^-5 mi) near a 220 deg transfer. The text says "The effect of the
    Earth is considerably larger than that of Mars". The printed axis units do not match Fig. 3's scale the way the text implies, and I did not resolve this.
  - Time-of-arrival bias is negative: "Time-of-arrival is typically several hours early because of the attraction toward the destination" (summary, p. 3).
- No numerical example of the full method is printed (no table of a worked flyby).

## 4. Relation to held work

- Perko 1974 (SIAM J. Appl. Math 27:200, HELD, `perko-1974-periodic-orbits-restricted-three-body-...`) and Breakwell-Perko 1974 (Celest. Mech. 9:437, HELD, `breakwell-perko-1974-second-order-matching-...`;
  digest `2026-10-04-digest-breakwell-perko-1974-second-order-matching.md`): the 1974 paper's first-order (2.1.9), (2.1.18) are eqs. above (rho^(1) with the ln(tau0/tau) singular part) and its (b, delta, v_inf)
  are this paper's (Delta, i*, v_inf). The 1974 digest quotes "certain gross biases obtainable by quadrature" for this paper's result, which agrees with the abstract and p. 19 here.
- Perko 1967 (HELD, `perko-1967-method-error-estimation-...`): the error estimates behind the O(lambda^2) claim; this paper cites them only to the thesis (footnote p. 8).
- Guillaume 1975 (HELD, `guillaume-1975-...`) extends the matching to the non-O(1) approach range.

## 5. Citation mining (references, p. 20, read on the image)

1. Lagerstrom, P. A. & Kevorkian, J. (1963), "Earth-to-Moon Trajectories in the Restricted Three Body Problems", Journal de Mécanique 2(2), June 1963. **Not held**, no hit in `ls | grep -i lagerstrom`. On the wanted list, row 70 (as J. Mécanique 2:189). Already a candidate.
2. Perko, L. M. (1964), "Interplanetary Trajectories in the Restricted Three Body Problem", AIAA Journal 2(12), December 1964. **Not held** (held Perko files begin 1967). Wanted row 32 is Perko's 1964 *thesis* ("Asymptotic matching ...", Stanford), not this AIAA Journal paper. **New candidate** (a short journal form of the first-order matching).
3. Battin, G. (1963), "Astronautical Guidance", McGraw-Hill. **Not held** (no hit in `ls | grep -i battin` or in CORPUS_INDEX). Wanted row 47 covers only "Battin 1959/1999". Textbook; low priority; **new candidate only if the 1963 edition is needed** (transition matrices, pp. 38-49 cited for the conic).
4. Danby, J. M. A. (1964), "Matrix Methods in the Calculation and Analysis of Orbits", AIAA Journal 2(1), January 1964. **Not held.** The held Danby file is the 1965 "matrizant" paper (`danby-1965-matrizant-keplerian-motion-...`), a different paper. Wanted row 42 lists "Danby 1964, AIAA J 2:16 and 2:13", so it is already a candidate (row 42).
