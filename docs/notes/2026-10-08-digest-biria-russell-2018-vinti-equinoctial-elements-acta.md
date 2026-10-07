# Digest: Biria & Russell 2018, "Equinoctial elements for Vinti theory: Generalizations to an oblate spheroidal geometry" (#960)

A. D. Biria and R. P. Russell (University of Texas at Austin), "Equinoctial elements for Vinti theory: Generalizations to
an oblate spheroidal geometry", *Acta Astronautica* **153**, pp. **274-288** (December 2018).
- **DOI 10.1016/j.actaastro.2017.11.013.** Crossref confirms the title, the authors Ashley D. Biria and Ryan P. Russell,
  volume 153, pages 274-288, and issue date 2018-12. Elsevier PII S0094-5765(17)31236-5, article reference AA 6544.
- Received 5 Sep 2017, accepted 13 Nov 2017 (cover page).
- Conference form: IWSCFF 17-75, 9th International Workshop on Satellite Constellations and Formation Flying, Boulder CO,
  June 2017. The manuscript is typeset with the IWSCFF header. **Not held.**
- File: upload `biria2018 (1).pdf`, 29 pp. (1 Elsevier cover page and 28 manuscript pages), md5
  **9e595fa54feeb6525d00b4db187f096c**.
  - This is the **accepted manuscript**, not the typeset version of record. The cover page says it "will undergo
    copyediting, typesetting, and review", and its PDF subject field reads "Accepted manuscript".
  - Page numbers below are manuscript page numbers. The PDF page is one more.
  - The text layer is digital and good, with an "ACCEPTED MANUSCRIPT" watermark that intrudes on some lines. No OCR is
    needed.
- Proposed filename:
  `biria-russell-2018-equinoctial-elements-vinti-theory-oblate-spheroidal-geometry-acta-astronaut-153-274-accepted-manuscript-doi-10.1016-j.actaastro.2017.11.013.pdf`
- Companion code: `VintiEquinoctialElements.zip`, which holds `Vinti_Equinoctial_Online.zip` v1.10 (Aug 2018), MATLAB,
  GPL v3 or later. See `code-notes.md`.
- How I read it:
  - I read the full text layer.
  - I read these page images: manuscript p.3 (the element definition, Eq. 2), p.16 (the Fig. 2 test parameters, at
    200 dpi) and p.24 (step 1 of the element-to-ECI algorithm, c^2 = Re^2 J2).
  - Eq. 3 (the classical-to-equinoctial map) is from the text layer only.

## 0. Verdict

**It is a coordinate-transformation paper. It gives no propagator, no accuracy against integration, no runtime and no
pinnable output. It matters to the project only as the element layer of the later nonsingular Vinti propagator: JAS 2019
or 2020, not held, wanted row drafted.**
- **What it is.**
  - It defines "oblate spheroidal (OS) equinoctial elements" {p, q1, q2, p1, p2, L} for Vinti's problem.
  - It derives the exact point transformations both ways between inertial r, v and these elements.
  - It removes the e = 0 and I = 0 singularities of the ECI-to-element Jacobian, which the CMDA 2018 relative-motion
    paper left open.
  - It removes the long-standing on-pole singularity with a "pole patch" approximation.
- **J2 only.** "The selected flavor of Vinti theory utilizes the symmetric potential wherein J3 = 0 so that the origins
  of the oblate spheroidal (OS) and ECI reference frames coincide" (p.2). Both algorithms start with "Compute
  c^2 = Re^2 J2" (pp.22 and 24; page image for p.24). The J3 version is left to future work (p.2).
  - So the **nonsingular equinoctial route captures J2 and the implied J4 = -J2^2 only**.
  - The **J2 + J3 route** (CMDA 2018) uses classical spheroidal elements. Its map from ECI to elements stays singular at
    e = 0 and near equatorial orbits.
  - The code matches the papers: the JAS supplement propagator has no J3 at all (0 occurrences of "J3" in
    `VintiPropEOE_adb.m`).
- **What it gives the project.** Background only, unless a Vinti propagator is ever built. For a clean-room build, the
  two step-by-step algorithms (pp.22-26) are the specification of the element layer.
- **PROPOSALS only:** see sec. 4. They are the same as the CMDA digest, sec. 5. Vinti is not a propagator for the
  multi-body legs of #997/#1000/#968.

## 1. The Vinti intermediary in this paper

- Potential: the original, symmetric Vinti 1961 potential (z_delta = 0), with c^2 = Re^2 J2 (p.24, page image).
- What it captures: J2 exactly. The higher even zonals follow the oblate-spheroid series J4 = -J2^2, and so on. That
  implied J4 is a body-dependent fraction of the real J4 (about 72% for the Earth; CMDA 2018 p.3). J3 is not included.
- Order of the analytic solution: the paper notes that "the approximate analytical solution ... has been developed to
  the third order in J2 in the literature" (abstract, p.1).
- Geometry. The spheroidal elements are defined by the classical-to-equinoctial formulas, with spheroidal e, I, omega'
  and Omega' (Eq. 3, p.3, text layer):
  - q1 = e cos(omega' + K Omega') and q2 = e sin(omega' + K Omega').
  - p1 = tan^K(I/2) cos Omega' and p2 = tan^K(I/2) sin Omega'.
  - L = f + omega' + K Omega'.
  - K = +1 for direct and -1 for retrograde elements.
- Spheroidal conic equation: rho = p / (1 + q1 cos L + q2 sin L) (Eq. 11, p.4).
- Spheroidal latitude equation for J3 = 0: eta = sin I sin psi (Eq. 13). For J3 not zero: eta = P + Q sin psi (Eq. 14).
- The OS equinoctial frame **rotates**, unlike the fixed frame of the spherical case. Its rate is the rate of the
  spheroidal ascending-node vector (p.3).
- Validity:
  - The transformations are not valid for nearly rectilinear orbits, where the spheroidal semilatus rectum is very small,
    because of the focal-circle forbidden zone (p.2).
  - Eq. (83) for the spheroidal RAAN rate "is invalid" for nearly parabolic, parabolic or hyperbolic orbits. Getchell's
    equations are the route for those (text layer, p.14).
  - The conclusions claim validity "for nearly or exactly zero-energy orbits, while generally maintaining validity for
    other bounded or unbounded orbit regimes" (p.27).
  - Note the mismatch with the code: the shipped code says "Vinti's solution is limited to bounded orbits"
    (`rv2oseoe_adb.m` and `oseoe2rv_adb.m` headers).

## 2. Accuracy claims

- The claim is about **transformation precision**, not propagation accuracy.
  - The test is a forward and back round trip, ECI to OS equinoctial to ECI. "Perfect coordinate transformations would
    give an exactly zero error" (p.16).
  - Fig. 2 (p.15) plots the preserved digits of X, Y, vx and vy against co-latitude, for 1e-15 to 1 deg.
- Text claims (p.16):
  - Near the pole, the exact equations lose X and Y precision to zero for co-latitudes "as large as 10^-8 degrees".
  - With the second-order pole patch, the precision diverges from the exact case at co-latitudes between 10^-1 and 10^-2
    deg (1 - |eta| ≈ 10^-7). The remaining X, Y loss is "an artifact of the equinoctial elements themselves". The
    spherical equinoctial transform has it too, reaching 1 correct digit near 10^-15 deg.
  - The velocities keep "around 16 digits near the poles".
- The pole switch used in the algorithms is 1 - |eta| < 10^-7 (pp.21-26). The conclusions recommend the approximate
  RAAN-rate expression within about 0.01 deg of the pole for Earth (p.27). The approximation is "accurate to O(J2^n) for
  an arbitrary order n, but analytical solutions in the literature do not exceed the third order" (p.27).
- No comparison with numerical integration appears anywhere in this paper.

## 3. Test case and runtime

- **Fig. 2 test inputs (p.16, page image; the text layer agrees):**
  - rp = 7000 km, eK = 0.1, IK = 90 deg, OmegaK = 30 deg, omegaK = 11 deg, with fK varied from 60 to 79 deg.
  - Earth constants: mu_e = 3.986 x 10^5 km^3/s^2, Re = 6378.137 km, J2 = 1.0826 x 10^-3.
  - The shipped code uses mu = 3.986004415e5 and J2 = 1.0826360229840e-3, and its driver uses a different orbit
    (rp = 7000 km, e = 0.01, I = 40 deg). So the code driver does **not** reproduce Fig. 2 as shipped.
  - The output is a digits-preserved plot only. There are **no published output values to pin**. The round-trip test
    itself is a natural self-check for any clean-room port: it needs no golden value.
- **Runtime:** no claim in this paper.

## 4. PROPOSAL (short; full version in the CMDA digest, sec. 5)

- For #997/#1000 (Earth-Moon, e.g. Schwaniger's 177 km perigee) and #968/#1004 (Jovian): this paper adds no propagation
  capability. The J2-only nonsingular propagator is in the JAS paper (not held), and the J2 + J3 one is in CMDA 2018.
- If a Vinti oracle is wanted, use the JAS v1.10 code (J2 only, nonsingular) for near-equatorial or near-circular test
  arcs. Use the CMDA v3.20 code (J2 + J3) where J3 matters and the orbit is away from e = 0 and I = 0.
- For a clean-room implementation, this paper's two algorithm summaries (pp.22-26) plus the forward and back round trip
  are the specification and the self-test for the element layer.

## 5. Print slips and notes (as printed)

- The abstract says the transforms are "exact except near the poles". The conclusions say they are valid for "bounded or
  unbounded" orbits. The code says bounded only. This is a scope difference, not a slip, but a consumer should not assume
  hyperbolic support from the shipped code.
- This is the accepted manuscript. The typeset version of record (Acta 153:274-288) may differ in wording or equation
  numbering. I made no identity check, because the version of record is not held.

## 6. Citation mining

Checked as for the CMDA digest. None of the items below is held or on the `#960` wanted list.

| Ref | Item | Held? |
|---|---|---|
| 1, 2, 4, 17 | Vinti 1959, 1966a, 1966b, 1961 (J. Res. NBS). DOIs as printed: 10.6028/jres.063B.012, .070B.002, .070B.003, .065B.017 | not held; free from NIST |
| 3, 18 | Biria & Russell AAS 16-537 / CMDA 2018 | CMDA form: this batch (digest-biria-russell-2018-relative-motion-vinti.md) |
| 5 | Broucke & Cefola 1972, Celest. Mech. 5:303, doi 10.1007/BF01228432 (equinoctial elements) | not held |
| 6 | Hintz 2008, JGCD 31:785, doi 10.2514/1.32237 (survey of element sets) | not held (only Hintz's 2023 textbook is held) |
| 7 | Danielson et al. 1995, NPS-MA-95-002 (DSST) | not held |
| 8 | Gim & Alfriend 2005, CMDA 92:295, doi 10.1007/s10569-004-1799-0 | not held |
| 9 | Getchell 1970, JSR 7(4):405, doi 10.2514/3.29954 | not held. Candidate |
| 10 | Garfinkel, Hori & Aksnes 1970, AJ 75:651, doi 10.1086/111000 | not held |
| 11 | Vinti 1969, AJ 74:25, doi 10.1086/110770 | not held |
| 12, 13, 15, 16 | Izsak 1960 SAO SR 52; Lang 1969 MIT thesis; Bonavito TN D-3562; Der & Bonavito 1998 | not held |
| 14 | Alfriend et al. 1977, Celest. Mech. 16:441, doi 10.1007/BF01229287 | not held |
| 19 | Biria 2017, PhD thesis, UT Austin | not held. Candidate |
| 20 | Code dataset v1.0 (2017), russell.ae.utexas.edu/index_files/vinti.html | code v1.10 held in this batch's zips |

- Follow-on paper: Biria & Russell, "Analytical Solution to the Vinti Problem in Oblate Spheroidal Equinoctial Orbital
  Elements", JAS 67:1-27 (doi 10.1007/s40295-019-00179-y). This is the "second paper" announced on p.2. Its code is held
  (zip 3), but the paper is **not held**. See `wanted-row.md`.

*Filed as `cyclers_pdf/papers/biria-russell-2018-equinoctial-elements-vinti-theory-oblate-spheroidal-geometry-acta-astronaut-153-274-accepted-manuscript-doi-10.1016-j.actaastro.2017.11.013.pdf`. The authors' MATLAB code (GPL v3 or later) is in the private corpus supplements/ directory.*
