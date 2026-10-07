# Digest: Biria & Russell 2018, "A satellite relative motion model including J2 and J3 via Vinti's intermediary" (#960)

A. D. Biria and R. P. Russell (University of Texas at Austin), "A satellite relative motion model including J2 and J3 via
Vinti's intermediary", *Celestial Mechanics and Dynamical Astronomy* **130**(3), article **23**, 40 pp. (2018).
- **DOI 10.1007/s10569-017-9806-4.** Crossref confirms the title, the authors Biria and Russell, volume 130, issue 3,
  article number 23, and publication date 2018-02-23.
- Received 5 Dec 2016, revised 22 Nov 2017, accepted 5 Dec 2017 (p.1).
- Conference form: AAS 16-537, AAS/AIAA Space Flight Mechanics Meeting, Napa CA, Feb 2016, Adv. Astronaut. Sci. 158,
  pp. 3475-3494. Source: the reference list and the code readme. **Not held.**
- File: upload `biria2018.pdf`, 40 pp., md5 **894ca2c08404a888153432fde634240c**. This is the Springer publisher PDF. Its
  metadata reads "Acrobat Distiller 10.1.8, modified using iText 5.3.5", dated 2018-02-23. The text layer is digital and
  good, so no OCR is needed.
- Proposed filename:
  `biria-russell-2018-satellite-relative-motion-model-j2-j3-vinti-intermediary-cmda-130-23-doi-10.1007-s10569-017-9806-4.pdf`
- Companion code: `VintiRelativeMotion_adb.zip`, v3.20 (July 2018), MATLAB, GPL v3 or later. See `code-notes.md`.
- How I read it:
  - I read the full text layer (`pdftotext -layout`).
  - I read these page images at 200 dpi: p.3 (the scope of Vinti and Brouwer theory), p.7 (Eq. 4, c and z_delta, the
    validity condition), p.10 lower half (Getchell speed figures), p.11 (Eq. 13, the O(J2^4) claim), p.21 (the truth model
    and integrator), p.22 lower part (the 1 km statement, the benchmark remark) and p.24 (Table 2, Fig. 6).
  - Table 2 has two witnesses: the page image and the publisher text layer. Its digits agree.
  - I read the figures (Figs. 4-7) only as orders of magnitude on log axes. Any figure value below is marked "read from
    figure".

## 0. Verdict

**It is an analytic first-order state transition matrix (STM) for relative motion. The STM is built on a cleaned-up
Vinti J2+J3 propagator. It gives the project the method, a validity map and GPL MATLAB code. It gives no pinnable
published output values.**
- **What it is.**
  - It revisits Vinti's 1966 spheroidal method. That method captures J2 exactly and J3 exactly through a polar origin
    shift. It also captures the part of J4 equal to -J2^2.
  - It fixes several numerical problems in Bonavito's 1966 procedure:
    - Getchell's quartic factorisation.
    - A cancellation-safe formula for alpha2^2 - alpha3^2.
    - A polar-orbit-safe initialisation of the slowly varying node angle.
    - The b2 = 0 singularity in the partial derivatives.
  - It then derives the analytic partials, so the STM lives in oblate spheroidal (OS) element space.
- **What it is NOT.**
  - Its accuracy plots are errors of a *linearised relative-motion* model: the deputy position from the chief plus the
    STM-mapped offset. They are not the absolute accuracy of the Vinti propagator against an integrator. The paper never
    reports a stand-alone propagator-against-integrator error for one orbit.
  - It is a one-body theory: one oblate primary and nothing else. No third body. Bounded orbits only (Fig. 1 is the
    "bounded Vinti orbit propagator"; the code readme says the same).
- **What it gives the project.**
  1. A description of what a Vinti intermediary captures, and to what order (sec. 1).
  2. The validity limits: J3^2 < 4 J2^3, the focal-circle forbidden zone, and the remaining ECI-to-element singularities
     at e = 0 and near equatorial orbits (sec. 2).
  3. Third-party GPL MATLAB code that could serve as an *outside oracle* for a numerical J2/J3 force term (sec. 5,
     PROPOSAL).
- **No golden values.** Table 2 publishes inputs only: chief elements and relative-orbit parameters. Every result is a
  log-scale plot. Nothing here can be pinned as an expected number.
- **PROPOSALS only (no code or catalogue change):** see sec. 5.

## 1. The Vinti intermediary: what it captures

- Potential (Eq. 4, p.7, page image): V = -mu (rho + eta z_delta) / (rho^2 + c^2 eta^2), in oblate spheroidal coordinates
  (rho, eta, phi). Eq. 3 (p.6) gives the map to ECI: X = sqrt(rho^2 + c^2) sqrt(1 - eta^2) cos phi, Y the same with
  sin phi, Z = rho eta.
- Fitted constants (p.7, page image): c^2 = Re^2 J2 (1 - J3^2 / (4 J2^3)) and z_delta = -Re J3 / (2 J2). "For the Earth,
  c ≈ 210 km and z_delta ≈ 7 km."
- Validity: "The specific condition for model validity is J3^2 < 4 J2^3" (p.7). Otherwise c^2 < 0 and no analytic Vinti
  trajectory exists.
- What it captures:
  - **J2 and J0 exactly.** The potential is fitted to the zeroth and second zonal harmonics exactly (p.2).
  - **J3 exactly**, by the origin shift z_delta (p.2 to p.3: the 1966 potential "fits the J3 harmonic exactly").
  - **J4 in part.** An exact Vinti reference "notionally includes the contributions of J2 + J4 + 2 J6 + ..." (p.2), with
    the implied J4 = -J2^2. On p.3: "approximately 72% of J4 ... the portion of J4 included can be larger or smaller
    depending on the central body, i.e., -J2^2/J4 ≈ 72% for the Earth" (page image).
  - For other bodies, the fraction is -J2^2/J4. Worked example from the repo's own sourced Uranus constants
    (`titania_oberon_realeph_895.py`: J2 = 3510.685e-6, J4 = -34.2e-6 at 25,559 km, Jacobson 2014): 1.2325e-5 / 3.42e-5
    ≈ 36%. This is my arithmetic, not the paper's. I give no Jupiter value here: the repo has no sourced Jupiter J2 or J4,
    and I did not check one. The consumer should compute it from a sourced pair.
- Order of the solution:
  - Vinti's 1966 solution has secular terms through O(J2^3) and periodic terms through O(J2^2) (p.3, page image).
  - This paper computes u = 1/S1 from Getchell's Eq. 13 rather than from Vinti's cubic approximation. That raises the
    secular terms to **O(J2^4)** (p.11, page image).
  - The mean frequencies are limited by the secular coefficients B1, B2 and B3, correct to O(J2^3), O(J2^4) and O(J2^4)
    (p.12, text layer).
- Brouwer for comparison (p.3): secular terms through O(J2^2), periodic terms through O(J2), then J3, J4 and J5 secular
  and long-period terms. Brouwer fails near the critical inclination of 63.4 deg. Vinti 1969 has no singularity there
  (p.3).

## 2. Formulation (pp.5-20)

- Elements:
  - The constant Vinti orbital elements are [a, e, S, l0, g0, beta3] (Eq. 5, p.8).
  - a = (rho1 + rho2)/2 and e = (rho2 - rho1)/(rho2 + rho1), where rho1 and rho2 are roots of the quartic F(rho) (Eqs. 8-9).
  - S = -eta0 eta1 = Q^2 - P^2, from the roots of the quartic G(eta) (Eqs. 10-12). S ≈ sin^2 I.
  - The STM state is the time-varying set [a, e, S, v, psi, Omega'] (p.9), where psi is the true argument of latitude and
    Omega' is Vinti's 1969 slowly varying node angle.
- Quartic factoring (sec. 2.2): the paper uses Getchell's 1970 iteration, with the Getchell-Monuki second-iteration start,
  instead of companion-matrix eigenvalues.
  - Convergence: "typically ... in no more than five iterations each" (p.10).
  - For large J2 (> 1e-1), more iterations are needed, and "the eigenvalue approach may be preferable" (p.10).
- Cancellation (sec. 2.2.3): compute alpha2^2 - alpha3^2 by Eq. 16 for near-equatorial orbits. The paper says this
  problem "is not mentioned in the literature" (p.11).
- Polar orbits (sec. 2.3): a new atan2 formula for Omega', Eq. 19 (p.13). A Walden 1968 differential correction is an
  option on an exact pole.
- The b2 = 0 singularity (sec. 2.4): the authors show it is "artificial" in the partials (p.13).
- STM (sec. 3, Eq. 1 and Fig. 3): Phi_S(t, ti) in OS element space. It maps to ECI or LVLH by linear (Jacobian) or
  nonlinear "bookend" transforms.
  - Partials are checked against complex-step derivatives and central differences (p.5).
  - Singularities left over (abstract, p.1): the linear map from ECI to spheroidal elements is still singular at e = 0
    and for nearly equatorial orbits. The companion Acta paper (the other digest in this batch) removes these for the
    J3 = 0 potential only.

## 3. Accuracy claims (relative-motion STM, not the propagator)

- Truth model (p.21, page image):
  - ECI equations with zonal J2 to J5.
  - Both spacecraft are propagated separately with a 4th-order variable-step Runge-Kutta integrator in double precision,
    over 15 orbits, with an accuracy tolerance of 2.3 x 10^-14.
  - The initial states are formed in quad precision.
  - The comparison metric is the deputy ECI position error.
- Quotable text claims:
  - "The accuracy of the Vinti model can be one to two orders of magnitude higher than that of the GA model" (p.22).
    The GA model is the Gim-Alfriend J2-only first-order STM.
  - "After 15 orbits, the aK = 12,000 km, IK = 63.4 deg, eK = 0.4, delta ri = 200 km case leads to errors on the order of
    1 km" (p.22, page image).
  - Against a two-body plus J2-only truth (Fig. 6, p.24), "the Vinti-based STM is performing slightly worse; its errors
    are still O(J2^2)". The cause is that the secular rates match only to first order.
  - Against a numerically integrated *Vinti* truth (Fig. 7, p.25), "the secular error growth appears to be substantially
    mitigated". Read from the figure: errors stay at about 1e-5 to 1e-4 km over 50 h. This is a figure reading, not a
    published number.
- Error sources the paper names (p.23):
  1. A force-model mismatch: Vinti has J2, J3 and about 72% of J4; the truth has J2 to J5.
  2. The truncation order of the theory.
  3. The linearisation.
- Scale for #997: the paper's chief orbits all have a = 12,000 km (about 5,600 km altitude at e = 0). The paper gives no
  case at a perigee like 177 km.

## 4. Test cases and runtime

- **Table 2 (p.24), two witnesses (page image and text layer):**
  - Chief: aK = 12,000 km; IK = 0, 30, 63.4 or 90 deg; OmegaK = 10 deg; omegaK = 20 deg.
  - Fig. 4 row: eK in [0, 0.4]; vKi = 0; phi_v = 0; delta lambda_v = 0; delta ri = 0.01 km; alpha = 30; beta = 30.
  - Fig. 5 row: eK = 0.4; the same angles; delta ri in [0.01, 100] km; alpha = 30; beta = 30.
  - The relative-orbit parameterisation is from Biria & Russell 2015, JGCD 38(8):1452 (not held).
  - Slip, as printed: Table 2 gives delta ri up to 100 km for Fig. 5, but the text on p.22 and the Fig. 5 legend speak of
    200 km (delta r_avg = 2E+02 km). This is probably a mean size against an initial offset. I left it as printed.
  - Earth constants: none are printed in the evaluation section. The code driver uses Re = 6378.137 km,
    mu = 3.986004415e5 km^3/s^2, J2 = 1.0826360229840e-3 and J3 = -2.5324353457544e-6. These are code defaults, not
    published values.
- **No published output values.** All results are log plots (Figs. 4-7). Nothing can be pinned.
- **Runtime claims (p.10, page image):**
  - In a MATLAB test with J2 = 1.08 x 10^-3, Getchell's factorisation is "roughly 18 times faster for factoring F(rho)
    and 14 times faster for factoring G(eta)" than the eigenvalue method. That example took 7 iterations for F and 4 for
    G.
  - When forced to 7 iterations, G is still about 10x faster (p.11).
  - The paper says the speed of Vinti against Brouwer "has not been benchmarked in the present investigation" (p.22).
    The last comparison it cites is Bonavito et al. 1969 (NASA TN D-5203).

## 5. PROPOSAL: use for #997/#1000 (Earth-Moon) and #968/#1004 (Jovian)

What the repository has today (checked by grep):
- Saturn J2 in the V4 gauntlet: `data/validation/v4_saturn.py`, with the body-agnostic `_j2_acceleration_kms2` (it takes
  mu, j2 and r_eq_km), plus `v4_saturn_strict.py`.
- Uranus J2 in the V4 gauntlet and its strict variant: `v4_uranus.py` and `v4_uranus_strict.py`.
- Uranus J2 + J4 in `search/titania_oberon_realeph_895.py`. It uses `_zonal(n, ...)`, a degree-n zonal term with a
  gradient, so it can feed an STM.
- **No Earth J2 anywhere in `src/`.** The Earth-Moon CR3BP and both regularised cores are point-mass.
- **No Jupiter J2.** `nbody/jovian.py` excludes it on purpose: "Deliberately EXCLUDED: Jupiter J2, solar tide, smaller
  moons" (lines 37-41). The lane exists to measure the gap from patched conics to continuous moon gravity, so adding J2
  would measure a different gap.
- `nbody/forces.py` and `core/` have no zonal term.

Proposals:
1. **Vinti cannot propagate the target trajectories.** It is a one-body intermediary. A Schwaniger-type leg is shaped by
   the Moon as much as by the Earth, and the Galilean lanes are shaped by the moons. For #997/#1000, the direct path to
   an Earth J2/J3 sensitivity check is a numerical zonal acceleration in the existing integrator. Reuse `_zonal` (it has
   a gradient for the STM) or `_j2_acceleration_kms2`.
   - Scale (my arithmetic): at 177 km altitude, the J2 to central-acceleration ratio is about
     1.5 J2 (Re/r)^2 ≈ 1.5e-3, with Earth J2 = 1.0826e-3 from the code driver.
2. **Name the frame choice first.** This is the real question for the Earth-Moon lane.
   - Earth's true pole is tilted from the Moon's orbit plane, by about 23.4 deg ± 5.1 deg. A real-pole J2 in the
     rotating frame breaks the Jacobi integral and the xz-plane mirror symmetry. `correct_symmetric_fixed_jacobi` relies
     on both.
   - A J2 with its pole along the rotating z axis keeps both, but it is a fictitious model.
   - The owner or lead should rule which variant #997/#1000 tests: a symmetric toy, or a real pole that needs an
     asymmetric corrector.
3. **Vinti's useful role is an outside oracle.**
   - Propagate a pure Earth + J2 + J3 arc (no Moon) with the new numerical zonal term and with the authors' MATLAB code,
     run in Octave outside the repo. Compare the two.
   - Fig. 7 of this paper shows the correct oracle setup: integrate the *Vinti potential itself* (Eq. 4) numerically, so
     the J4 = -J2^2 part matches.
   - Label any such numbers "third-party computed (Biria-Russell code v3.20)". They are not published values.
4. **For #968/#1004 (Jupiter):** keep the `jovian.py` exclusion. A Jupiter J2 belongs only on a separate, higher-fidelity
   rung. Vinti could later replace the Kepler legs of a patched-conic or Lambert stage with oblate legs. That needs a Vinti
   Lambert or targeting method, and the 2020 CMDA paper "The Lagrange coefficients of Vinti theory" (Biria, CMDA 132(5),
   article 26, doi 10.1007/s10569-020-09966-4, Crossref-checked, not held) is the lead for it.
   - Getchell's factorisation is fine at Jupiter's J2 (about 1e-2). The paper flags trouble only above 1e-1 (p.10).
5. **Licence.** The code is GPL v3 or later (see `code-notes.md`). The repo has no LICENSE file that I found. A direct
   port into the repo would bring the GPL with it.
   - Recommendation: if a Vinti propagator is ever wanted, build it clean-room from Bonavito TN D-3562 (a NASA report),
     Getchell 1970 and the published equations here.
   - Use the MATLAB code only as an outside oracle.

## 6. Print slips (as printed; not corrected)

- Table 2 gives delta ri in [0.01, 100] km for Fig. 5. The text (p.22) and the Fig. 5 legends say 200 km.
- The PDF's metadata title drops the subscripts: "including and via Vinti's intermediary". This is cosmetic.

## 7. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` (vinti, getchell, bonavito, alfriend, brouwer, kozai, izsak, walden,
wiesel, garfinkel, biria) and with CORPUS_INDEX and the `#960` wanted list. None of them is held or listed.

| Ref | Item | Held? |
|---|---|---|
| Vinti 1959, 1961, 1962, 1963, 1966a,b | J. Res. NBS 63B:105; 65B:169; 66B:5; 67B:191; 70B:1; 70B:17. Free from NIST. The DOIs for 1959, 1961 and both 1966 papers are printed in the Acta paper's reference list (10.6028/jres.063B.012, .065B.017, .070B.002, .070B.003) | not held. **Candidate (free)**: 1966b (J3 inclusion) is the theory used here |
| Vinti 1969 | AJ 74:25-34, doi 10.1086/110770 (as printed in the Acta reference list) | not held |
| Bonavito 1966 | NASA TN D-3562: the computational procedure this paper modifies | not held. **Candidate (free, NTRS)** |
| Bonavito, Watson & Walden 1969 | NASA TN D-5203: Vinti against Brouwer, accuracy and speed | not held. Candidate (free, NTRS) |
| Walden & Watson 1967; Walden 1968 | NASA TN D-4088; AIAA J. 6(7):1305 | not held |
| Getchell 1970 | J. Spacecr. Rockets 7(4):405-408, doi 10.2514/3.29954 (printed in both papers) | not held. **Candidate**: the factorisation used |
| Izsak 1960 | SAO Special Report 52 | not held |
| Lang 1969 | MIT MS thesis: unbounded Vinti orbits | not held |
| Der & Bonavito 1998 | AIAA Progress in Astronautics and Aeronautics vol. 177 (Vinti's collected work) | not held (book) |
| Gim & Alfriend 2003 | JGCD 26(6):956-971, doi 10.2514/2.6924 (as printed) | not held |
| Brouwer 1959; Kozai 1962 | AJ 64:378; AJ 67:446 | not held |
| Wiesel 2015 | numerical action-angle Vinti solution (cited p.3) | not held |
| Biria & Russell 2015 | JGCD 38(8):1452, doi 10.2514/1.G000622: the Table 2 parameterisation | not held |
| Biria 2017 | PhD thesis, UT Austin, "Revisiting Vinti Theory ..." (cited by the Acta paper as ref. 19) | not held. Candidate (likely free from the UT repository; not checked) |

- Wanted-list suggestions (proposal), in priority order, only if a Vinti oracle or propagator is pursued:
  1. Biria & Russell 2020, JAS 67:1-27 (row drafted in `wanted-row.md`).
  2. Vinti 1966b, NBS.
  3. Bonavito TN D-3562.
  4. Getchell 1970.
  5. Biria 2017 thesis.

*Filed as `cyclers_pdf/papers/biria-russell-2018-satellite-relative-motion-model-j2-j3-vinti-intermediary-cmda-130-23-doi-10.1007-s10569-017-9806-4.pdf`. The authors' MATLAB codes (GPL v3 or later) for this paper, the Acta 2018 paper and the JAS 2020 paper are in the private corpus supplements/ directory; code notes are filed beside this PDF.*
