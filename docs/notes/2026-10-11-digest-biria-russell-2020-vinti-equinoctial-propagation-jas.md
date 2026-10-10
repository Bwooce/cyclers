# Digest: Biria & Russell 2019/2020, "Analytical Solution to the Vinti Problem in Oblate Spheroidal Equinoctial Orbital Elements" (#960 batch 43)

A. D. Biria and R. P. Russell (University of Texas at Austin), *The Journal of the Astronautical Sciences* **67**(1),
pp. **1-27**. Published online 2019-05-21; issue dated March 2020.
- **DOI 10.1007/s40295-019-00179-y.** Crossref (checked 2026-10-11) confirms the title, the authors Biria and Russell,
  volume 67, issue 1, pages 1-27, online 2019-05-21, print 2020-03.
- Conference form: AAS 17-659, AAS/AIAA Astrodynamics Specialist Conference, Stevenson WA, Aug 2017, Adv. Astronaut.
  Sci. 162:1369-1386 (from the supplement readme). **Not held.**
- File: owner upload `biria2019.pdf`, 27 pp., md5 **af84bca8439220f8dc9cd4b3694c05ba**. Springer publisher PDF
  (LaTeX, Acrobat Distiller 9.0, dated 2019-05-17/21). Digital text layer, good; no OCR needed. PDF page n = journal
  page n.
- Proposed filename:
  `cyclers_pdf/papers/biria-russell-2020-analytical-solution-vinti-problem-oblate-spheroidal-equinoctial-elements-jas-67-1-doi-10.1007-s40295-019-00179-y.pdf`
  (it matches the held supplement zip stem).
- Companion code (held): `cyclers_pdf/supplements/biria-russell-2020-...-code-supplement-v1p1.zip`, md5
  ee709d9bafb96e13dab715ae5ba2488a, v1.10 (Nov 2024), MATLAB, GPL v3 or later. See the batch-41 code notes, sec. 3.
- Companion YAML: `biria-russell-2020-test-cases.yaml` (both published examples, pin levels, conventions block).
- Checks: `checks.py` / `checks.out` (ours; never goldens).
- **How I read it:**
  - I read the full text layer (`pdftotext -layout`).
  - I read these page images at 150 dpi: p.2 (GTOC 9 speed-up), p.4 (Eq. 3), p.18 (Eqs. 92-98), p.19 (Eqs. 99-103),
    p.21 (Fig. 3, start of Examples), p.22 (Example 1 text, Fig. 4), p.23 (Figs. 5, 6), p.24 (Example 2 inputs,
    Fig. 7), p.25 (validation remark, Conclusions).
  - Other equations below are marked "text layer". The text layer drops primes (for example B1' in Eq. 103 reads
    "B1") and garbles fractions.

## 0. Verdict

**This is the propagation half of the Biria-Russell spheroidal equinoctial work. It gives a nonsingular analytic
solution of the J2-only (symmetric) Vinti problem for bounded orbits. It is the paper behind the held MATLAB
supplement. It prints no propagated state, no accuracy number and no runtime for its own propagator. So nothing here
can pin a state.**
- **What it is.** It rewrites Vinti's three kinematic equations (Vinti 1966 with J3 = 0) in the OS equinoctial
  elements {p, q1, q2, p1, p2, L} of the Acta 2018 paper. It derives the equinoctial "integrals of the motion" (the
  secular parts at t_i) from the spheroidal Delaunay variables. It removes the polar singularity in the node rate. It
  then solves an equinoctial generalised Kepler equation by iteration.
- **Bounded only.** Abstract, p.1: "a fully nonsingular analytical solution for bounded orbits". Conclusions, p.25: the
  same. The code prints a message and executes `return` before its outputs are assigned if spheroidal e >= 1 - 1e-4 or
  if the effective e' = gamma1 e >= 1 - 1e-4 (`VintiPropEOE_adb.m` lines 576-583 and 1113-1117). In MATLAB the caller
  then stops with an "output argument not assigned" error. For k = c^2/p^2 >= 0.1716 it only warns.
- **Test cases.** Two worked examples (sec. 3 below; YAML). Neither prints a state.
  - Example 1 (Earth LEO, J2 = 5.08e-2): one printed scalar relation, with incomplete inputs. **Weak golden.** We
    reproduce it (637.6 km against "almost 640 km"). e_K = 0 follows from the printed constant spherical p = 7,300 km
    = rp (1 + e), and I = 0 from "equatorial". Only Re is assumed (6378.137 km); the result does not depend on mu.
  - Example 2 (Saturn orbit insertion): the initial state is complete, but the outputs are figure-level only.
    **Figure-level check.** Our arithmetic shows that the printed inputs fit the printed "about 10 revolutions in a
    little over a year".
  - The supplement's built-in driver state matches **neither** example.
- **What it gives the project.** It is the published specification of the held supplement code, so the code is now
  checkable against a source. It also gives a self-check that needs no golden value: the constants of the motion must
  be preserved to about 15 digits (p.25, sec. 2). The consumers' needs are unchanged from batch 41 (sec. 5): there is
  still no Earth or Jupiter J2 in `src/`.
- **PROPOSALS only:** sec. 5.

## 1. The method (solution steps)

- **Potential** (Eq. 3, p.4, page image): V = -mu rho / (rho^2 + c^2 eta^2), the symmetric Vinti potential.
  - c^2 = Re^2 J2 (Eq. 4, p.5, text layer).
  - It is exact for an oblate spheroid with J4 = -J2^2, J6 = +J2^3, and so on.
  - It includes "roughly 72% of J4" for the Earth (p.5) and "roughly 28% of J4" for Saturn (p.24).
  - **J3 is not included.** Shifting the OS origin for J3 "hamstrings the definition of an equinoctial element set"
    (p.3, text layer).
- **Order:** the secular terms are O(J2^3) in Vinti's form. They become O(J2^4) when u comes from Getchell's
  factorisation, as implemented. The periodic terms are O(J2^2) (pp.6-7, text layer).
- **Elements** (Eq. 21-23, p.8, text layer; the same as the Acta 2018 digest sec. 1): p, q1 = e cos(omega' + K Omega'),
  q2 = e sin(omega' + K Omega'), p1 = tan^K(I/2) cos Omega', p2 = tan^K(I/2) sin Omega', L = f + omega' + K Omega'.
  K = +1 direct, -1 retrograde. This is the *modified* set (p, L), not the standard set (a, lambda). It is "nonsingular
  for the full range of eccentricity and inclination" (p.8).
- **Step A, kinematic equations in equinoctial form** (Eqs. 24-51, pp.9-14, text layer).
  - Define the eccentric longitude F = E + omega' + K Omega' (Eq. 24) and the mean longitude
    lambda = M + omega' + K Omega' (Eq. 25).
  - Add angle combinations to both sides of Vinti's Eqs. 18-20. This gives new constants beta~_j (Eqs. 29-31).
  - Periodic terms: rescale the coefficients A~_jk = A_jk / e^k and B~_jk = B_jk / Q^k (Eqs. 32-33). These must be
    computed by "omission of e^k and Q^k" in the original formulas, never by division.
  - Rewrite e^k sin kf and Q^k sin k psi as polynomials in (q1, q2) and (p1, p2) times sin/cos kL. This uses
    Chebyshev's recursion (Eqs. 37-48).
  - Result: Eqs. 49-51, free of singular angles.
- **Step B, equinoctial integrals of the motion** (Eqs. 52-91, pp.14-18, text layer).
  - Split each angle into secular and periodic parts.
  - The secular parts are the spheroidal Delaunay variables: l = M_s, g = omega'_s, h = Omega'_s (Eqs. 61-63),
    and psi_s = l + g (Eq. 67). The paper notes that omega'_s = g "has never been given explicitly" (p.15).
  - The unknowns omega'_p, Omega'_p and lambda_s at t_i come from Eqs. 89-91. q1s, q2s, p1s, p2s then come by
    rotations (Eqs. 59-60).
- **Polar singularity removed** (Eqs. 83-85, p.16-17, text layer). The node rate h-dot in Eq. 82 has a 1/sqrt(1-S)
  factor. It is replaced by Eq. 84, the rate of Vinti's 1969 slowly varying Omega'. That form has no singularity at
  the pole.
- **Step C, secular propagation** (Eqs. 92-96, p.18, page image).
  - p is constant.
  - (q1s, q2s) rotates by (g-dot + K h-dot) Delta t.
  - (p1s, p2s) rotates by h-dot Delta t.
  - lambda_s(t) = lambda_s(t_i) + (l-dot + g-dot + K h-dot) Delta t.
- **Step D, generalised Kepler equation** (Eqs. 97-107, pp.18-20; Eqs. 97-103 page image).
  - Start with lambda_j = lambda_s and the secular q, p values.
  - Solve "lambda_j = F_j - gamma1 q1_j sin F_j + gamma1 q2_j cos F_j" (Eq. 97), with
    "gamma1 = a/a0 = a/(a + b1) <= 1" (Eq. 98). This is Getchell's notation.
  - The solver is a variable-order Laguerre method (Conway 1986, as used by Der).
  - Then get L_j from sin L_j and cos L_j (Eqs. 99-101) and atan2 (Eq. 102).
  - Update the periodic parts omega'_p, Omega'_p and lambda_p (Eqs. 103-107). Then update q and p with the inverse
    of Eqs. 59-60.
  - Repeat with lambda_j = lambda_s + lambda_pj until converged.
  - The code caps this outer loop at `NitrEOEmax` = 8 with a tolerance of 1e-12 (driver). It prints a message if it
    does not converge (lines 1350 and 1517).
- **ECI output.** The element-to-ECI map is the Acta 2018 algorithm. It is not repeated here.

## 2. Accuracy and runtime claims

**The paper gives no numeric accuracy or runtime for this propagator.** These are the only related statements:
- **400x speed-up (p.2, page image).** This is for "a modified Vinti propagator [14]" (Der & Bonavito 1998) used in
  GTOC 9. It was compared with "an eighth order Runge-Kutta numerical integrator taking 50 steps per revolution" for a
  10-day LEO transfer. The paper says its own solution "would lead to qualitatively similar behavior". **It is not a
  measurement of this propagator.**
- **Accuracy against integration is out of scope.** p.21: "The accuracy of Vinti dynamics relative to low-fidelity
  spherical harmonic representations was assessed in prior work [1, 5] and is not discussed here" (page image).
  [1] is the Biria 2017 thesis and [5] is the CMDA 2018 paper.
- **About 15 digits (p.25, page image).** "The constants are preserved to double precision or about 15 digits only
  when the exact expression for Omega'-dot is used". This is about how well the ECI-to-elements transforms preserve
  the constants of the motion. Traditional approximations do not preserve them. It is not propagation accuracy.
  - Near a pole, the Acta pole-patch approximation "still preserves the constants of the motion".
  - Corner case: for an equatorial orbit, the total angular momentum is also conserved.
  - **This is a self-check that needs no golden value.** At several output times, transform back to elements and
    check that p, e and the other constants keep about 15 digits.
- **Example 1, p.22:** the spherical elements used for comparison were "obtained by numerically integrating Hamilton's
  equations of motion under the Vinti potential". The overlay (Fig. 4) is visual only. No error is printed.
- **Example 2, p.25:** "after one year, the effect of neglecting 72% of J4 will manifest itself as a sizable phase
  error". No number is given.

## 3. Published test cases (two witnesses: page image and text layer agree; full detail in the YAML)

| Case | Page | Inputs printed | Inputs missing | Printed output | Pin level |
|---|---|---|---|---|---|
| Example 1, Earth LEO | pp.21-23, Figs. 3-5 | rp = 7,300 km; "nominally circular equatorial"; J2 = 5.08 x 10^-2; "approximately 18 orbits" | mu, Re, e_K, I, angles, t_i, tf | spherical p constant at 7,300 km; spheroidal p "almost 640 km smaller" (p.22). Fig. 3a axis about 6662 km (figure reading) | weak |
| Example 2, Saturn insertion | pp.23-25, Figs. 6-7 | mu = 3.7931 x 10^7 km^3/s^2; Re = 60,330 km; J2 ≈ 1.6291 x 10^-2; J4 ≈ 9.36 x 10^-4 (as printed); rp = 61,330 km; e_K = 0.99; I_K = 10 deg; Omega_K = 30 deg; omega_K = 11 deg; f_K = 6 deg (osculating Keplerian) | t_i and exact tf | "a little over a year (about 10 revolutions)"; "roughly 28% of J4"; element plots (Fig. 7) | figure-level |

**Our checks** (`checks.py` / `checks.out`; ours, not goldens):
- **Example 1.** e_K = 0 follows from the printed facts (the spherical p is a constant 7,300 km, and p = rp(1 + e) with
  rp = 7,300 km); I = 0 from "equatorial". Only Re is assumed. With Re = 6378.137 km and any mu, the equatorial Vinti turning points give a spheroidal
  p of **6662.36 km**, a gap of **637.6 km**, and a spheroidal e of 0.0691.
  - With mu doubled, the digits are the same. Re = 6378.000 km moves the gap by 0.03 km. Re = 6400 km gives 642.7 km.
  - This matches the printed "almost 640 km", the Fig. 3a axis, and the about 0.07 swing of the spheroidal q1 and q2 in
    Fig. 3a.
  - So the published relation is reproducible. It is a usable positive control for a clean-room spheroidal-element
    layer. It is still not a golden for the propagator.
- **Example 2, initial elements.** Eq. 22 with K = +1 on the printed inputs gives p = 122,046.7 km, q1 = 0.7472,
  q2 = 0.6495, p1 = 0.0758, p2 = 0.0437 and L = 47.00 deg. These match the t = 0 values of p, q1, q2, p1 and p2 in Fig. 7 (figure reading). L at t = 0 sits on the steep
  rise through periapsis and is not legible, so I claim no match for L.
- **Example 2, why a 0.99 orbit makes 10 revolutions in a year.** The Keplerian a_K = 6.133 x 10^6 km has a period of
  179 d. That cannot fit 10 revolutions in a year.
  - At 61,500 km from Saturn's centre, the Vinti potential makes the energy -7.946 km^2/s^2, against the Keplerian
    -3.092. The real orbit is much more tightly bound.
  - The 2-body-equivalent values are a ≈ 2.39 x 10^6 km, period ≈ 43.5 d, and 10 revolutions ≈ 1.19 yr. The apoapsis
    is about 4.7 x 10^6 km, which fits the Fig. 6 extent of about 4.5 x 10^6 km.
  - So the printed inputs are self-consistent and not misprinted. Lesson for consumers: osculating Keplerian elements
    at a low periapsis around an oblate body can misstate the orbit size badly.
- **Supplement driver against the paper.** The driver uses Earth (code values, not published): Re = 6378.137 km,
  mu = 3.986004415e5, J2 = 1.0826360229840e-3, with rp = 6700 km, e = 0.10, I = 40 deg, Omega = 30 deg, omega = 11 deg,
  f = 6 deg, 5 revolutions.
  - It matches **neither** example. It shares only the three angles with Example 2.
  - To run Example 2, set the Saturn constants, rp = 61330, ecc = 0.99, inc = 10 deg and `sims.nOrb = 10`.
  - The v1.10 fix (Nov 2024, "quadrant ambiguity in the true and eccentric longitude, found in rare cases") came after
    the paper's figures. Any rerun uses later code than the paper did.
  - Not run here. MATLAB is not installed, and Octave was not tried.

## 4. How it differs from the two 2018 papers

| | CMDA 2018 (130:23) | Acta 2018 (153:274) | **JAS 2019/2020 (this paper)** |
|---|---|---|---|
| Potential | Vinti 1966, J2 + J3 (origin shift z_delta) | symmetric, J2 only | symmetric, J2 only |
| Elements | classical spheroidal [a, e, S, l0, g0, beta3] | OS equinoctial {p, q1, q2, p1, p2, L} | the same OS equinoctial set |
| Content | Vinti propagator (classical elements) plus an analytic relative-motion STM | point transforms between ECI and elements both ways, pole patch, no propagation | **time propagation**: equinoctial kinematic equations, equinoctial integrals of the motion, generalised Kepler equation |
| Singularities | ECI-to-element map singular at e = 0 and near-equatorial | removed (pole patch for 1 - abs(eta) < 1e-7) | removed throughout, including the polar node rate (Eq. 84) |
| Orbits | bounded | the paper says bounded or unbounded; the code says bounded | bounded only |
| Numeric results | Table 2 inputs; log plots of relative error | Fig. 2 round-trip digits plot | two examples, plots only |
| Code | v3.20 (2018) | v1.10 (2018) | v1.10 (2024) |

- New here: the equinoctial Kepler equation (Eq. 97); the explicit statement omega'_s = g (p.15); and the
  Omega'_s-based h-dot that removes the 1/sqrt(1-S) factor (Eq. 84-85).
- The paper says an analytic STM in these elements is future work: the solution "prescribes an analytical state
  transition matrix (STM) that is nonsingular" (p.25). It does not derive one.

## 5. PROPOSALS for the consumers (no code or catalogue change)

- **The batch-41 finding still holds at 89edd72b.** `grep -rn -i j2 src | grep -v -i j2000` (full output, no
  truncation) and a grep for an Earth or Jupiter value (`1\.0826|14696|earth.*j2|jupiter.*j2`) find no Earth or
  Jupiter J2 in `src/`. `nbody/jovian.py:37` still says "Deliberately EXCLUDED: Jupiter J2".
  - J2 terms that do exist: Uranus J2 + J4 in `search/titania_oberon_realeph_895.py`; Uranus J2 in
    `src/cyclerfinder/data/validation/v4_uranus.py` (+ `_strict`); Saturn J2 in
    `src/cyclerfinder/data/validation/v4_saturn.py` (+ `_strict`).
  - The body-agnostic helper `_j2_acceleration_kms2(r_km, *, mu, j2, r_eq_km)` (`v4_uranus.py:293`, Vallado Eq. 8-37,
    re-exported by `v4_saturn.py`) takes mu, J2 and Re. It is the natural numerical comparator for a Vinti oracle and
    the starting point for an Earth J2 term in #997/#1000.
- **#997/#1000 (Earth-Moon low-perigee legs).** This paper changes nothing in the CMDA digest's sec. 5. Vinti is a
  one-body intermediary and cannot propagate a Moon-shaped leg. The direct route is still a numerical Earth zonal term
  in the existing integrator, after the frame choice is ruled.
  - What this paper adds: a J2-only nonsingular oracle for pure Earth + J2 test arcs, including near-circular,
    near-equatorial and polar arcs. The CMDA J2 + J3 code is singular for those.
  - The oracle should be checked against a numerical integration of the Vinti potential (Eq. 3) itself, so that the
    implied J4 matches. That is how Example 1 builds its spherical comparison (p.22).
  - Example 1's 637.6 km gap is a positive control for any spheroidal-element layer we build.
- **Galilean lanes (Jupiter J2).** Keep the `jovian.py` exclusion. A Vinti leg needs only mu, Re and J2 from a sourced
  set. I give no Jupiter value: none is sourced in the repo.
  - Example 2 shows the regime: a high-e, low-periapsis capture around a strongly oblate planet. There the osculating
    Keplerian a is wrong by a factor of 2.6. That is a warning for any patched-conic Jupiter capture stage with a low
    perijove.
  - Vinti targeting would need the Biria 2020 CMDA "Lagrange coefficients" paper (not held; see the CMDA digest).
- **Clean-room route, if one is ever wanted.** The specification is this paper's Steps A-D, the Acta algorithms and
  Vinti 1966b / Getchell 1970. The tests are the constants-preservation self-check (sec. 2) and the Example 1 gap.
  The GPL code stays an outside oracle, as in batch 41.
- **Corpus bookkeeping (proposal):**
  - File the PDF under the proposed name and add a CORPUS_INDEX paper row next to the supplement row (line 206).
  - Remove wanted row 63, which this file fills.
  - Update the "not held" lines in the Acta digest (sec. 0, sec. 4 and the follow-on note) and in the code notes
    (sec. 3).
  - Optional: file the YAML as `data/sources/biria-russell-2020-vinti-examples.yaml`.

## 6. Print slips and notes (as printed; not corrected)

- **Saturn J4 is printed as "≈ 9.36 x 10^-4", positive** (p.24). In the usual convention Saturn's J4 is negative. The
  "roughly 28%" checks out as J2^2/abs(J4) = 2.654e-4 / 9.36e-4 = 28.4% (my arithmetic). This is a sign-convention
  matter. Jacobson et al. 2006 is not held, so I did not check the source value.
- p.24 says Vinti "captures approximately 72% of J4 for the Earth" and p.25 says "neglecting 72% of J4" for Saturn.
  Both agree with 28% captured at Saturn. This is not a slip.
- Example 1 says "approximately 18 orbits" (p.21), and Fig. 3's time axis ends near 9.2 x 10^4 s. The Keplerian
  period at 7,300 km with Earth mu is about 6,200 s (about 15 orbits). This is not a slip: as in Example 2, the J2
  potential binds the orbit more tightly. The spheroidal a = 6694 km from checks.py gives a period of about 5,450 s,
  so about 17 orbits in 9.2 x 10^4 s (my arithmetic, assuming Earth mu, which the paper does not print). That is
  close to the printed "approximately 18".
- The text layer drops the prime on B1' in Eqs. 103-104 (page image p.19). Read the equations from the page images.

## 7. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` (biria, vinti, getchell, bonavito, broucke-cefola, gim-alfriend,
walker-ireland, conway-1986, jacobson 2006, buffington-strange, wiesel, zurita, wu-tong, der-bonavito, garfinkel,
lara-gurfil, deprit-ferrer, danielson, izsak, biscani), with CORPUS_INDEX and with the `#960` wanted list.

| Ref | Item | Held? |
|---|---|---|
| 1 | Biria 2017, PhD thesis, UT Austin | not held. Candidate (the accuracy study is deferred to it) |
| 2, 5 | Biria & Russell AAS 16-537 / CMDA 130:23 (2018) | CMDA held (batch 41) |
| 3, 4 | Biria & Russell IWSCFF 17-75 / Acta 153:274 (2018) | Acta accepted manuscript held (batch 41) |
| 6 | Biscani & Izzo 2016, MNRAS 455:3480, doi 10.1093/mnras/stv2512 (Euler two-fixed-centres problem) | not held |
| 7 | Broucke & Cefola 1972, Celest. Mech. 5:303, doi 10.1007/BF01228432 | not held |
| 9 | Buffington & Strange 2008, AAS 07-276, patched-integrated gravity assist (Saturn oblateness) | not held (other Buffington papers are held) |
| 10 | Conway 1986, Celest. Mech. 39:199, doi 10.1007/BF01230852 (Laguerre-Kepler solver) | not held |
| 13, 14 | Der 2011 AMOS (superior Lambert); Der & Bonavito 1998, AIAA PA&A vol. 177 | not held |
| 17 | Getchell 1970, JSR 7(4):405, doi 10.2514/3.29954 | not held. Candidate (as in batch 41) |
| 18 | Gim & Alfriend 2005, CMDA 92:295, doi 10.1007/s10569-004-1799-0 | not held |
| 20 | Jacobson et al. 2006, AJ 132:2520, doi 10.1086/508812 (Saturn gravity field) | not held. Source of the Example 2 J2 and J4 |
| 23-29 | Vinti 1959, 1961, 1962, 1963, 1966a, 1966b (J. Res. NBS), 1969 (AJ 74:25) | not held. Free from NIST; 1966b is the base theory |
| 30 | Walden & Watson 1967, NASA TN D-4088 (corrections to the B_j, B_jk coefficients) | not held. Candidate if a clean-room port is built |
| 31 | Walker, Ireland & Owens 1985, Celest. Mech. 36:409, doi 10.1007/BF01227493 (modified equinoctial elements) | not held |
| 32, 33, 35 | Wiesel 2015 JGCD 38:1757; Wright 2016 AFIT thesis; Zurita 2017 AMOS | not held |
| 34 | Wu & Tong 1981, Chin. Astron. Astrophys. 5:192 (third-order Vinti, Omega'_s formula) | not held |
| 8, 11, 12, 15, 16, 19, 21, 22 | Brouwer & Clemence 1961; Danielson et al. 1995; Deprit & Ferrer 1987; Ferrándiz & Floría 1993; Garfinkel & Aksnes 1970; Izsak 1962; Lara & Gurfil 2012; Nayfeh 2004 | not held. Background |

- **Wanted list:** this file fills row 63 (Biria & Russell 2020 JAS). No other reference above is on the list. The
  batch-41 priority suggestions still stand: Vinti 1966b, Bonavito TN D-3562, Getchell 1970, Biria 2017 thesis. Add
  Walden & Watson TN D-4088 only if a clean-room port is started.

*Filed as `cyclers_pdf/papers/biria-russell-2020-analytical-solution-vinti-problem-oblate-spheroidal-equinoctial-elements-jas-67-1-doi-10.1007-s40295-019-00179-y.pdf`. Check scripts, outputs, witness tables and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/biria-russell-2020-test-cases.yaml`.*
