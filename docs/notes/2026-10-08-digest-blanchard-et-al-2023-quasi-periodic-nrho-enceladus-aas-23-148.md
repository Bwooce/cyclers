# Digest: Blanchard, Lo, Restrepo, Landau, Anderson & Elschot 2023, "Quasi-Periodic Near Rectilinear Halo Orbits for Enceladus Jet Targeting"

J. T. Blanchard (Stanford), M. W. Lo, R. L. Restrepo, D. Landau, B. D. Anderson (JPL) & S. Elschot (Stanford), **AAS
23-148** (printed top of p.1). Venue: AAS/AIAA Space Flight Mechanics Meeting, Austin TX, Jan 2023 (the venue is stated
in the paper's own ref [18], "Austin, TX, 2023", and a web search dates the paper 15 Jan 2023). I found **no DOI**
(Crossref returned no match) and I did not find a proceedings DOI.

- Source file: upload `3f7db4af-CL23_4039.pdf` (JPL clearance CL#23-4039), 18 pp, md5
  d35c38b50003d12e0da7f408a9f8f7d0. Digital, good text layer, no OCR needed.
- Proposed filename: `blanchard-lo-restrepo-landau-anderson-elschot-2023-quasi-periodic-nrho-enceladus-jet-targeting-AAS-23-148.pdf`
- Not held (grep of papers/ and CORPUS_INDEX for Blanchard, Enceladus QPO: no hit).
- How I read it: full text layer; the constants block (p.3), Figure 4 (p.5), the results text (pp.14-16) and Figs 8-10
  viewed on page images at 100 dpi. The paper has **no table**. Numbers quoted below are the printed constants and
  the stated results; I cross-checked them by arithmetic (section 3).

## 0. Verdict

A short method-and-demonstration paper. It collects the single-shooting **invariant-circle corrector** for 2D quasi-periodic
tori (Olikara-Scheeres / McCarthy-Howell form, FFT rotation operator) in a form meant to be easy to code, builds
"several dozen" QPOs about northern L2 (and L1) NRHOs of the **Saturn-Enceladus CR3BP**, shows one isoenergetic
four-member set (Fig. 8), and uses an **invariant funnel** for phase targeting (20 m/s per revolution).

What it gives the project:

1. **A clean written-out corrector for the torus lanes (#996, #916):** unknowns, error vector, Jacobian blocks, phase
   constraints and initial guess are all printed (eqs 21-54, section 1). Good as an implementation checklist against
   the held Baresi-Olikara-Scheeres 2018 and the Olikara 2016 thesis.
2. **A Saturn-Enceladus CR3BP constants block** (mu and the three normalising constants, section 2) that checks out
   (section 3).
3. **One concrete base orbit:** an L2 NRHO with rp = 340 km, period 12.8 h, C = 3.0000356, stability index 3.71 (p.14-15).
   That is a seed for a test of any Enceladus torus code, but **no initial state, no rotation number and no QPO
   parameters are printed**, so it is not a golden by itself. The paper points to the JPL periodic-orbit database
   (ssd.jpl.nasa.gov/tools/periodic_orbits.html) for the NRHOs.
4. **For the Enceladus cycler context:** nothing on Titan-Enceladus cyclers or any cycler. The paper is about
   orbits around Enceladus, which are bound orbits, not transfers between moons. The only cross-link is the held
   Davis-Phillips-McCarthy 2018 (Saturnian ocean-worlds orbiters), which is ref [6].

PROPOSALS only: (a) use the paper as the written reference for an invariant-circle corrector test in #996 (seed from
the monodromy eigenvector, alpha = 1e-5; rho from eq. 54); (b) do **not** list it as a source of QPO or rotation-number
values; (c) add the sister paper (ref [18]) and Olikara-Scheeres 2010 to the wanted list (rows 71-72 below).

## 1. Method as printed (pp.6-14, eqs 15-54)

- CR3BP (eqs 1-12, p.2-3): `mu = m2/(m1+m2)`; Omega = (x^2+y^2)/2 + (1-mu)/|r1| + mu/|r2| (**no** `mu(1-mu)/2`
  constant); `C = 2 Omega - v.v`. Primary at -mu, secondary at 1-mu. (Compare: the held Johnston-Lo-Mortari slides
  add the `+ mu(1-mu)` constant to C; with mu = 1.9e-7 the offset is 1.9e-7, which only touches the seventh digit of
  3.0000356.)
- Stability index (eq. 14): `nu = (|lambda_max| + 1/|lambda_max|)/2`; nu = 1 linearly stable. Figure 4 (p.5) plots nu and
  C versus periapsis distance rp for the L1 and L2 halo families: **nu stays below 2 (it peaks near 2 at rp about
  100 km), equals 1 for rp about 230 to 290 km, then rises rapidly (about 7 at rp = 380 km)**. The text calls
  rp < about 290 km "NRHO" (p.5). The dashed vertical line in Fig. 4 at about 250 km is unlabelled; it sits at
  Enceladus's radius (252 km), which I note but the paper does not say.
- Invariant funnels (eqs 15-20, pp.6-7): planar case, a circle dF0 of states at distance d from x*, velocity parallel to
  x*, Jacobi constant C*; spatial case, a 3D closed surface `(x-x*)^T A (x-x*) = 1` with `atan2(ydot,xdot)` fixed and
  C fixed, A = diag(1/a^2, 1/b^2, 1/c^2, 0, 0, 1/d^2) (d much smaller than a, b, c); integrate backward in time.
- Torus (pp.8-14): invariant circle `U` with N points (N odd), `theta_i = 2 pi i / N`; state `u` stored about the halo
  state x0; stroboscopic map `F(U, T)` over one longitudinal period; rotation undone in Fourier space with
  `R(rho) = D^-1 Q(rho) D`, `Q = diag(exp(-i k rho))`; `G = R F`.
  Unknowns `chi = [ {U}; T; rho ]` (size Nn + 2); error `gamma = [ {G} - {U}; C(U) - C*; gamma_psi; gamma_theta ]`
  (size Nn + 3, eq. 47); Jacobian blocks eqs 38-51 (including `(R(rho) kron I6) Phi_tilde`, `dG/dT = f(x0 + g_i)`,
  `dG/drho = D^-1 (dQ/drho) D F`, and the Jacobi-constant gradient eq. 46); phase constraints
  `gamma_psi = <U, dU~/dpsi>`, `gamma_theta = <U, dU~/dtheta>` after Schilder et al. 2005; update `d chi = -J^+ gamma`
  (Moore-Penrose); stop at `||gamma||_inf < eps`, "e.g. 1e-12".
- Initial guess (eqs 52-54): `u0(theta) = alpha (Re v6 cos theta - Im v6 sin theta)` from the eigenvector of the
  complex eigenvalue lambda6 (Im > 0), `alpha` about 1e-5; `T0` = NRHO period; `rho0 = atan(Im lambda6 / Re lambda6)`.
  (A plain `atan` of the ratio loses the quadrant; `atan2` would be safer. My remark, not the paper's.)

## 2. Numbers printed

Read on page images (pp.3, 5, 14-16) and text layer; two witnesses agree on every item below. Machine-readable copy:
`blanchard-2023-saturn-enceladus-constants-and-results.yaml`.

| item | value as printed | page |
|---|---|---|
| mu (Saturn/Enceladus) | 1.9011097359e-7 | 3, eq. 2 |
| R (length unit) | 2.38042e5 km | 3, eq. 3 |
| T (time unit) | 1.88389e4 s | 3, eq. 4 |
| V (velocity unit) | 12.6357 km/s | 3, eq. 5 |
| Enceladus SOI radius | about 488 km (under 2 x equatorial radius) | 1 |
| Europa SOI; Moon SOI | about 9,700 km; about 66,000 km | 1 |
| Saturn-Earth light time | 1.1 to 1.5 h | 1 |
| NRHO / non-NRHO boundary | rp about 290 km | 5 |
| base L2 NRHO for Fig. 8 | rp = 340 km (altitude 88 km), T = 12.8 h, C = 3.0000356, nu = 3.71 | 14-15 |
| Fig. 9 periapsis altitudes (leftmost QPO) | colour scale 60 to 110 km | 15 |
| funnel size (spatial) | max distance 500 m, max velocity difference 1.3 m/s | 16 |
| phase-targeting cost | about 20 m/s per revolution | 16, 17 |

Not printed: any rotation number, any QPO period, any QPO Jacobi constant (the four Fig. 8 QPOs are isoenergetic at
the NRHO's C, p.14), any state vector, any torus-family table, any QPO width in km (Fig. 8 axes: x 1 to 1.004 and
y plus/minus 4e-3 nondimensional, read from the figure), any count other than "several dozen".

Findings stated in words: none of the computed QPOs intersects Enceladus's surface (p.15); wider QPOs spread more
over the south pole but their periapsis distances move away from Enceladus (minimum shrinks, maximum grows, Fig. 8);
L1 NRHOs have mirror families over the left of the south pole (p.15); the funnel found a position intersection for
the Fig. 9 QPO with a ~20 m/s jump (p.16); stable manifolds of a nu = 3.71 NRHO do not spread enough to intersect in one
period (p.15); with nu = 1 there is no stable subspace at all.

Wording issue (not an error): p.15 calls nu = 3.71 a "low stability index" and p.16-17 say NRHOs are "too stable" for
manifolds. nu = 3.71 means |lambda_max| = 7.28 (ours), so the orbit is unstable; the paper's point is that the
manifolds are too slow to be useful in one period.

## 3. Checks (our arithmetic; `checks.py`, `checks.out`)

- **V = R/T:** 2.38042e5 / 1.88389e4 = 12.63566 km/s; printed 12.6357. Agrees.
- **T equals the Enceladus orbital period over 2 pi:** 2 pi x 1.88389e4 s = 1.37000 d (Enceladus sidereal period
  1.370218 d from memory; matches to 4 digits). Not an independent source; plausibility only.
- **mu:** m2/(m1+m2) from rounded JPL GM values (GM Saturn 37931208, Enceladus 7.2112 km3/s2, from memory) gives
  1.90113e-7 against the printed 1.9011097e-7 (1e-5 relative; my rounded GM values are the limit). Agrees.
- **SOIs, R mu^0.4:** Enceladus 487.8 km (printed 488, agrees); Moon 65,859 km (printed 66,000, agrees); Europa
  9,725 km with a = 671,100 km and mu = 2.528e-5 (printed 9,700, agrees).
- **Altitude:** 340 - 252.1 (Enceladus mean radius, from memory) = 87.9 km; printed 88 km. Agrees.
- **nu = 3.71** gives |lambda_max| = 7.283; NRHO T = 12.8 h = 2.446 in nondimensional time, 0.389 of Enceladus's orbital
  period (about 7:18, not a claim of the paper). Fig. 4: at rp = 340 km the L2 curve of nu is on the rising branch
  (about 4 at 340 km by eye), consistent with 3.71.
- Not checked: no QPO can be rebuilt (no state), the corrector was not run.

## 4. Citation mining (27 references)

Held = `ls papers | grep -i` AND CORPUS_INDEX. No not-held item below has a wanted-list row today. Proposed rows from 71.

| ref | work | status |
|---|---|---|
| 1 | "Cassini at Enceladus", Science 311 (Mar 2006) 1422-1425, author list scrambled in the reference | not held; no row; low (pages suggest the Hansen et al. water-vapour paper; not checked) |
| 2 | Spencer et al. 2018, "Plume origins and plumbing", in Enceladus and the Icy Moons of Saturn, pp.163-174 | not held; no row; low |
| 3 | MacKenzie et al. 2020, "Enceladus Orbilander", LPI Contrib. 2547 | not held; no row; low |
| 4 | Szebehely 1967 | **held** (szebehely-1967-theory-of-orbits-...-book.pdf) |
| 5 | Pollard 1966, Mathematical Introduction to Celestial Mechanics | not held; no row; low |
| 6 | Davis, Phillips & McCarthy 2018, Acta Astronaut. 143:16, doi 10.1016/j.actaastro.2017.11.004 | **held** |
| 7 | Zimovan-Spreen, Howell & Davis 2020, CMDA 132, doi 10.1007/s10569-020-09968-2 | **not held, no row.** NRHO stability and resonances; proposed row 73, medium |
| 8 | Zimovan-Spreen, Howell & Davis 2022, JAS 69(3):718, doi 10.1007/s40295-022-00320-4 | not held; proposed row 74, low |
| 9, 10 | Parrish et al. 2020, AIAA 2020-1466 and 2020-1700 | not held; no row; low (NRHO operations) |
| 11, 12 | Davis et al. 2019, NRHO disposal and escape (AAS/AIAA 2019) | not held; no row; low |
| 13 | Boudad, Howell & Davis 2019, IAA-AAS-SciTech-039 | conference version not held; the **2020 ASR journal version is held** (boudad-howell-davis-2020-...-asr-66-2194...). |
| 14 | Phillips & Davis 2020, "Cloud computing methods for NRHO trajectory design", Adv. Astronaut. Sci. 171:3287 | not held; no row; low |
| 15, 16 | Blanchard et al. 2021, "Invariant funnels for resonant landing orbits"; "New tools for tour design: swiss cheese plot, invariant funnel, and resonant encounter map" | not held; no row; proposed rows 75-76, medium (the invariant-funnel lane) |
| 17 | Blanchard & Elschot 2022, MPC with invariant funnels as terminal sets | not held; no row; low |
| 18 | **Blanchard, Lo, Anderson, Landau, Restrepo & Elschot 2023, "Using NRHO invariant funnels to target Enceladus south pole", AAS/AIAA SFM, Austin** | **not held, no row. Sister paper (the funnel half of this study). Proposed row 71, medium-high** |
| 19 | Blanchard et al. 2020, "Low energy capture into high inclination orbits for ocean worlds missions" | not held; no row; low |
| 20 | **Olikara & Scheeres 2010, Adv. Astronaut. Sci. 145:911** (the original torus corrector) | **not held, no row.** The Olikara 2016 thesis and Baresi 2018 are held, which cover the same method. Proposed row 72, medium |
| 21 | McCarthy & Howell 2021, "Leveraging quasi-periodic orbits for trajectory design in cislunar space", Astrodynamics 5(2):139, doi 10.1007/s42064-020-0094-5 | **not held, no row.** Proposed row 77, medium (torus lane) |
| 22 | Baresi, Olikara & Scheeres 2018, JAS 65:157 | **held** |
| 23 | Olikara 2016 PhD thesis | **held** (olikara-2016-...-phd-thesis-colorado.pdf) |
| 24 | McCarthy 2022 PhD thesis, Purdue | not held; no row; low (thesis; check Purdue repository) |
| 25 | Schilder, Osinga & Vogt 2005, SIADS 4:459, doi 10.1137/040611240 | not held; no row; low (phase constraints) |
| 26 | Henry & Scheeres 2022, Proc. ACC, doi 10.23919/ACC53348.2022.9867282 | not held; no row; low |
| 27 | Henry & Scheeres 2021, "Expansion maps", JGCD 44(3):457, doi 10.2514/1.G005492 | not held; no row; low |

Proposed wanted rows:
- 71 | Blanchard, Lo, Anderson, Landau, Restrepo & Elschot (2023), "Using NRHO invariant funnels to target Enceladus south pole", AAS/AIAA SFM Austin | sister paper to AAS 23-148 | medium-high (#996 context)
- 72 | Olikara, Z. P. & Scheeres, D. J. (2010), "Numerical method for computing quasi-periodic orbits and their stability in the restricted three-body problem", Adv. Astronaut. Sci. 145:911-930 | original invariant-circle corrector; thesis and Baresi 2018 held | medium
- 73 | Zimovan-Spreen, Howell & Davis (2020), CMDA 132, doi 10.1007/s10569-020-09968-2 | NRHO stability and resonance, torus seeds | medium
- 74 | Zimovan-Spreen, Howell & Davis (2022), JAS 69(3):718-744, doi 10.1007/s40295-022-00320-4 | low
- 75 | Blanchard et al. (2021), "Invariant funnels for resonant landing orbits", AAS/AIAA SFM | invariant funnel origin | medium
- 76 | Blanchard et al. (2021), "New tools for tour design: swiss cheese plot, invariant funnel, and resonant encounter map", AAS/AIAA Astrodynamics Specialist Conf. | medium
- 77 | McCarthy & Howell (2021), Astrodynamics 5(2):139-165, doi 10.1007/s42064-020-0094-5 | QPO trajectory design, torus lane | medium

*Filed as `cyclers_pdf/papers/blanchard-lo-restrepo-landau-anderson-elschot-2023-quasi-periodic-nrho-enceladus-jet-targeting-aas-23-148-jpl-cl-23-4039.pdf`. Table transcription: `data/sources/blanchard-2023-saturn-enceladus-constants-and-results.yaml`.*
