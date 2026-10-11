# Digest: Danby 1964, "Matrix Methods in the Calculation and Analysis of Orbits" and "The Matrizant of Keplerian Motion", AIAA J. 2(1):13-19 (#960)

J. M. A. Danby (Yale University), two papers in AIAA Journal 2(1), January 1964. Both received 28 May 1963, revised 16 October 1963,
supported by Minneapolis-Honeywell Regulator Company. Wanted-list row "Danby (1964) matrix methods ... and the matrizant;
Goodyear and Guseman": both Danby items received; Goodyear and Guseman (1964, IBM/RTCC 12-008) stays wanted.

- **Paper A:** "Matrix Methods in the Calculation and Analysis of Orbits", AIAA J. 2(1):13-16, DOI 10.2514/3.2206.
- **Paper B:** "The Matrizant of Keplerian Motion", AIAA J. 2(1):16-19, DOI 10.2514/3.2207.
- **Filed as**
  - `cyclers_pdf/papers/danby-1964-matrix-methods-calculation-analysis-orbits-aiaa-j-2-1-13-doi-10.2514-3.2206.pdf`
  - `cyclers_pdf/papers/danby-1964-matrizant-keplerian-motion-aiaa-j-2-1-16-doi-10.2514-3.2207.pdf`
  - Table transcription: `data/sources/danby-1964-tables.yaml`.
  - The two files are the received files, unchanged, 4 pdf pages each, native text layers. md5 of the first is
    d4ef9a79e2ad878343dc17949abda3cc (printed pp.13-16) and of the second 1d6649b37dcf152d27502b891c3ca1b2 (printed pp.16-19).
  - **The second file is the companion paper, not a continuation.** Its running heads, title block ("The Matrizant of
    Keplerian Motion", J. M. A. Danby), abstract and DOI watermark (3.2207) are those of a separate paper. Printed p.16 is
    in both files: it ends paper A (section 4 and the references) and starts paper B. Paper A is complete in the first file
    and paper B is complete in the second, so the two files together cover pp.13-19 with p.16 twice.
- **How it was read:** all 8 pages on 150 dpi renders, then every table on 260 to 600 dpi crops. The text layer drops every
  fraction, exponent, tilde and sign in the formulas, so the formulas were read from the render. The check that matters is
  numerical (sec. 3): every table was evaluated at several orbits against the 6 x 6 state-transition matrix integrated from the
  variational equations. There are no figures in either paper.

## 0. Verdict

**Closed-form Keplerian state-transition matrix, tabulated for elliptic, hyperbolic and parabolic reference orbits, plus a
special-perturbation method that uses it. Nothing periodic or repeating appears, and there is no cycler.** The matrizant is the
state-transition matrix (STM) of two-body motion. The tables are exact where they are used, with seven printed slips or inconsistencies listed in sec. 3.
- **Paper B** is the reference. It gives the matrizant Omega(t,t0) = [[L, M],[P, Q]] of Keplerian motion in the orbital
  reference system (x to pericentre, y to true anomaly 90 degrees), with the time origin at a pericentre passage so that
  Omega(t,t0) = Omega(t,T) Omega^-1(t0,T). Four tables: eccentric anomaly (elliptic, hyperbolic by substitution), nearly circular
  (Bower's method, general axes), true anomaly (elliptic and hyperbolic), parabolic (w = tan(f/2)).
- **Paper A** uses it. It integrates special perturbations by quadrature of the matrizant against the disturbing force on a
  Keplerian reference arc, switches to a fresh arc when the quadratic residual (dr)^2 passes the tolerance, and uses a
  "mean orbit" through the true positions at the two ends of an arc. A circumlunar orbit was the test case (eight arcs gave
  an accuracy of 1e-3, about a mile). Afterwards, error analysis, guidance and changes of constants need only closed functions.
- **Cross-reference.** The corpus holds Danby 1965, AIAA J. 3(4):769, "Matrizant of Keplerian Motion" (digest
  2026-10-05-digest-danby-1965-matrizant-keplerian-motion.md). That note gives the same matrizant in terms of Cartesian
  position, velocity and elapsed time, with no anomaly and no Kepler equation (M and N forms). The 1964 tables here are in
  terms of the anomaly. Paper A's reference 2 cites "The matrizant of Keplerian motion, AIAA J. 2, 16-19 (1963)", which
  is paper B (the printed year 1963 is the received year).

## 1. What the papers do

- **Paper A, sec. 1-2 (pp.13-15).** The method of Danby 1962 (Astron. J. 67:287, general perturbations in rectangular
  coordinates) is applied to special perturbations. A dominant -mu r/r^3 force plus a disturbing f(r,t). A reference Keplerian arc
  satisfies the equations with f = 0. Residuals dr, dr' at t0 and t are linked by Eq. 1 (homogeneous part) plus a particular
  integral (Eq. 3). The accuracy of the first-order step is governed by (dr)^2: a new reference orbit is needed when it
  exceeds the tolerance (for example 5e-7 for six decimals). Two ways to continue: restart from an earlier instant with a new
  osculating orbit, or fit a mean orbit through the true positions at t0 and t1 (which is simplest to program and "nearly
  always substantially more effective").
- **Paper A, sec. 3 (pp.15-16).** Practical form: scaled components (tildes), eccentric anomaly E as the independent variable
  (dt = (r/(na)) dE, so Kepler's equation is solved only for output times), repeated Simpson quadrature with a fourth-difference
  step test, hyperbolic orbits by the substitution a -> -a, n -> n1, (1-e) -> (e-1), F for E. The removal of (1-e) from the
  denominators is stated without the algebra. Tables 1 and 2 list the scaled components. The accuracy needed from the
  quadrature is governed by (1-e)^(-3/2). The paper recommends Taylor series about tabulated values for interpolation.
- **Paper A, sec. 4 (p.16).** Analysis of a calculated orbit: for times ta and tb in different arcs, Omega(ta,tb) is the chain
  of the matrizants of the arcs between them (Eq. 6), each built from its own reference-orbit parameters, provided the axes
  keep their orientation in space. Neglected forces and small changes of constants (the mass of the Moon is the example) can be
  put into f in Eq. 3 along the stored reference orbits to give their effect. An unknown force can in principle be recovered
  by treating Eq. 3 as an integral equation.
- **Paper B, sec. 1 (pp.16-17).** The name "matrizant" is preferred as having historical precedence (Peano and Baker; Frazer,
  Duncan and Collar 1938). Bower 1932 (Lick Obs. Bull. 445) gave a method of differential corrections that can compute the matrizant
  for elliptic and hyperbolic orbits. Myachin 1959 (Bull. ITA 7:257) gave explicit components in the orbital reference system.
  The change of axes is M1 = R^T M R (Eq. 4). The derivative relations of Eq. 5 and the symplectic inverse and composition of
  Eqs. 6 and 7 let all other times and axes be built from the pericentre-based tables.
- **Paper B, sec. 2-5 (pp.17-19).** Table 1 (eccentric anomaly) and the hyperbolic substitution; Table 2 (nearly circular
  orbits, Bower's notation: e cosE0 and e sinE0 from r0 and r0', a difference form of Kepler's equation, general axes through
  Phi0); Table 3 (true anomaly f, valid for elliptic and hyperbolic, with pericentre distance q and nu = (mu/q^3)^(1/2));
  Table 4 (parabolic, w = tan(f/2)); nearly parabolic orbits through lambda = (1-e)/(1+e) and a series sigma. The paper warns that
  components become large for e near 1 because (1-e) is in the denominator, which is a property of the orbit and not a
  fault of the formulas.
- **Paper B, sec. 6 (p.19).** Use and limits: first-order theory, so the squares and products of residuals must be negligible
  against the required accuracy; it is dangerous to predict errors near pericentre from errors elsewhere in the orbit; tests on
  nearly circular orbits say nothing about orbits of moderate eccentricity.

## 2. Tables (transcribed in `data/sources/danby-1964-tables.yaml`)

| Table | Content | Components |
|---|---|---|
| Paper B Table 1 | elliptic, eccentric anomaly E | 20 (l, m, p, q at 11, 12, 21, 22, 33) |
| Paper B Table 2 | nearly circular, Bower, 3 x 3 blocks from Phi0 | 22 steps and four block equations |
| Paper B Table 3 | true anomaly, elliptic and hyperbolic | 20 |
| Paper B Table 4 | parabolic, w = tan(f/2) | 20 |
| Paper A Table 1 | elliptic, scaled components for Eq. 4 | 20 |
| Paper A Table 2 | hyperbolic, scaled components for Eq. 4 | 20 |

The printed equations (A: 1-6, B: 1-7) and their roles are in the YAML.

## 3. Independent checks, all ours (scratch scripts, not committed)

The reference is the two-body 6 x 6 STM integrated with the variational equations (DOP853, rtol and atol 1e-13, mu = 1).

- **Paper B Table 1:** all 20 formulas agree to 9 digits at four orbits. Transcription and printed table are both right.
- **Paper B Table 3:** all 20 agree at six orbits: three elliptic (e = 0.3, 0.6, 0.2) and three hyperbolic (e = 1.5, 2.5, 1.2). The same
  formulas serve both.
- **Paper B Table 4 (parabolic): two printed slips.** q12 and q21 do not reproduce the matrix as printed (w^2 leading power, read at
  600 dpi). With w^3 in both (q12 = 4w^3(1 + w^2/5 - w^4/5)/(1+w^2)^3 and q21 = 2w^3(2 + w^2)/(1+w^2)^3) all 20 agree at three orbits.
  The 18 other formulas agree as printed, including q22's 18/5 coefficient.
- **Paper B Table 2 (nearly circular): three printed values do not reproduce the matrix.** For e = 0 only the f' line fails.
  For e > 0, with all three corrections the table reproduces the full 6 x 6 matrix of a rotated elliptic orbit to 1e-12 at five orbits.
  - **f'** is printed "1 - F'/r0". It must be "-F'/r0", the time derivative of f = 1 - F/r0. As printed it adds 1 to the diagonal of P.
  - **J'** is printed "a(k/r) - aG". It must be "a(k/r) - aG'" (G' is defined three lines above). As printed, g' and Q are wrong.
  - **H in L** is printed H = a^(1/2) e sinE, used in L = (1/r)[3J + 2FH + G r0] and, with H' = (k/r) e cosE, in L'. The correct L uses
    the constant H0 = a^(1/2) e sinE0, and L' then has no 2FH' term. With the printed H, L is off by 0.077 at (a, e, E0, E) =
    (1, 0.1, 0.5, 2.0) and by 2.6 at (2, 0.5, 3.0, 5.0). The quantity r' = (k/r)H keeps H at E. This correction set is the one found
    and not proved unique. The paper's section 3 does not say H is taken at E0.
  - The relations of section 3, e cosE0 = r0 r0'^2/mu - 1, e sinE0 = r0 . r0'/(mu a)^(1/2) and n(t - t0) = dE + e sinE0(1 - cos dE) -
    e cosE0 sin dE, are exact (checked to 1e-15).
- **Paper A Tables 1 and 2:** each scaled component equals the paper-B component times a closed factor. The factors are l~ and p~: (r/a)(1-e)^k
  with k = 2, 1, 2, 1 for ij = 11, 12, 21, 22 and 5/2 for 33. For m~ and q~ they are n (r/a)(1-e)^k with k = -1/2, 1/2, -1/2, 1/2 for
  ij = 11, 12, 21, 22 and 1/2 for 33. Hyperbolic: the same with (e-1), n1 and a negative a. The ratios are constant in E (F) at
  (a, e) = (1, 0.3), (2, 0.3), (-1, 1.6) and (-2, 1.6), so all 40 scaled formulas are correct as printed.
- **Paper A Eq. 4:** the first term (response to dr0') reproduces the in-plane blocks of Omega(t,t0)[:, 3:6] to 1e-12 at three
  orbits, with the printed prefactor (a^2/(n r r0))(1-e)^(-3/2). **The z entries (the 33 components) come out low by (1-e)^(3/2):** the printed
  33 components carry powers 5/2 and 1/2, so their products have (1-e)^3 instead of the in-plane (1-e)^(3/2). Eq. 4 as printed
  is exact in the orbit plane and not in z. The paper does not say whether the 33 entries or the prefactor is intended.
- **Paper A Eq. 5:** printed with the 2 x 2 block matrix as [[P, M],[L, Q]], L and P interchanged against the definition used everywhere
  else. The result dr0' = M^-1 dr1 is right because M sits in the right place.
- **Paper B Eq. 5:** L = -dM/dt0, P = -d^2 M/(dt dt0), Q = dM/dt. They follow from the variational equation (our derivation, in the YAML).
  **Eq. 6:** the symplectic inverse [[Q^T, -M^T],[-P^T, L^T]] agrees with the numerical inverse to 1.5e-13.
- **Not checked:** the nearly parabolic series (lambda, sigma and the l21 example in section 5 of paper B), the numerical
  accuracy claims (about a mile, eight arcs, 1e-3), and the Myachin and Bower attributions.

## 4. Relevance to this project

- **Lineage of the project's STM code.** The project has an analytic Kepler STM (`core/kepler_stm.py`, universal variables, Shepperd
  lane) and a self-contained state plus STM co-integrator for the Jupiter-central model (`nbody/jovian_stm.py`, used by the
  EGGIE corrector and, per its docstring, written because REBOUND's variational particles do not differentiate a
  Python force callback). Danby 1964 is the earlier closed-form, anomaly-based route to the same Keplerian STM: the
  elliptic, hyperbolic and parabolic cases in separate tables, with the singular behaviour at e = 1 stated openly. The universal-variable form
  in `kepler_stm.py` avoids those separate tables. **Usable as an independent regression oracle:** Tables 1, 3 and 4 of paper B,
  with the exponent fix for Table 4 q12 and q21, give analytic STM components to compare with `kepler_stm.py` at elliptic,
  hyperbolic and near-parabolic orbits, with the numerical STM above as a second witness. For a unit test the expected values
  must come from the paper's printed formulas with the corrected exponents, not from our code.
- **Singular cases.** Danby's tables show the (1-e)^(-1) and (1-e)^(-2) denominators that make the anomaly form
  ill-conditioned for e near 1. The warning applies to any anomaly-based STM we write and does not apply to the universal-variable form.
- **Arc splitting and mean orbits.** The paper A strategy (reference arcs, restart when (dr)^2 exceeds the tolerance, mean orbit through
  the two end positions, chain of matrizants of Eq. 6) is the ancestor of multi-arc shooting. It has no periodicity or closure condition.
  The relation to the project's multiple-shooting correctors is that of method ancestry, not of shared formulas.
- **Periodic or repeating content:** none. No resonance, no repeated body, no closed orbit and no cycler. The only repeated
  structure is the chain of arcs in Eq. 6.
- **Catalogue and gate:** no catalogue implication and no literature-gate label (a mathematical method paper). Nothing collides with #942, #943 or #973.
- **Risk of misuse.** Do not code Paper B Table 2 or Paper A Eq. 4 from the page without the corrections of sec. 3. Table 4's q12 and q21 are wrong as printed.

## 5. Citation mining (checked against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Held:** Danby 1965, AIAA J. 3(4):769 (the matrizant follow-up, digest 2026-10-05).
- **Already on the wanted list:** Goodyear and Guseman (1964), IBM/RTCC 12-008 (same row).
- **Not in the index and not on the wanted list:**
  - **Danby, J. M. A. (1962), "Integration of the equations of planetary motion in rectangular coordinates", Astron. J. 67:287-299.**
    The method that paper A applies to special perturbations (paper A ref. 1, paper B ref. 5). Medium priority; ADS likely has it.
  - **Bower, E. C. (1932), "Some formulas and tables relating to orbit computation and numerical integration", Lick Obs. Bull. 445**
    (paper B ref. 2). The source of Table 2's method. Low to medium priority.
  - **Myachin, V. F. (1959), "On the estimation of error in the numerical integration of equations in celestial mechanics",
    Bull. Inst. Theor. Astron. (Akad. Nauk SSSR) 7:257-280** (paper B ref. 3). Low priority (Russian).
  - Danby, Fundamentals of Celestial Mechanics (1962), Herget, The Computation of Orbits (1948) and Frazer, Duncan and Collar,
    Elementary Matrices (1938): textbooks, not added.
- **Wanted-list action:** in the Danby row, both Danby items are marked received. Goodyear and Guseman stays.

*Filed as `cyclers_pdf/papers/danby-1964-matrix-methods-calculation-analysis-orbits-aiaa-j-2-1-13-doi-10.2514-3.2206.pdf` and `cyclers_pdf/papers/danby-1964-matrizant-keplerian-motion-aiaa-j-2-1-16-doi-10.2514-3.2207.pdf`. Table transcription: `data/sources/danby-1964-tables.yaml`.*
