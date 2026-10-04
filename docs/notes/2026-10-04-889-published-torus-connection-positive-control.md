# #889 — published positive control for the torus-connection lane: Kumar, Anderson & de la Llave (SIAM 2025) reproduced

**Date:** 2026-10-04. **Module:** `src/cyclerfinder/search/pertbp_strob_889.py`. **Tests:**
`tests/search/test_pertbp_strob_889.py` (13 tests, about 5 s). **Driver:**
`scripts/screen_889_published_torus_connection.py` (stages `po`, `tori`, `bundles`, `forward`,
`refine`, `verify`; variants `table1`, `printed`, `table1_hi`). **Outputs:**
`data/found/889_published_torus_connection/`.

**Sources.** The target paper is SIAM J. Applied Dynamical Systems 24(1):219-258 (2025), filed in
the private paper corpus as
`kumar-anderson-delallave-2025-gpu-connections-tori-perturbed-crtbp-siam-ads-24-219-arxiv-2109.14814.pdf`
("SIAM"). The model and torus method come from the companion Celest. Mech. Dyn. Astron. 134:3
(2022), filed as `kumar-anderson-delallave-2021-whiskered-tori-manifolds-cmda-arxiv-2105.11100.pdf`
("CMDA"). The periodic orbits and mass ratio come from arXiv:2109.14800, filed as
`kumar-anderson-delallave-2021-highorder-resonant-manifold-expansions-cnsns-arxiv-2109.14800.pdf`
("CNSNS").

## Bottom line

**Verdict: reproduced within a stated tolerance.** Our refined connection in the Jupiter-Europa
planar elliptic problem is
(x, y, px, py) = (−0.96063509, 0.88783033, −0.51376751, −0.64713768). The printed point is
(−0.96064, 0.88783, −0.51377, −0.64714).

- The differences are (4.9, 0.33, 2.5, 2.3) × 1e-6.
- y, px and py round to the printed five decimals.
- x lies on the rounding boundary. Its difference is 4.9e-6 to 5.1e-6 across our runs, against a
  printed half-unit of 5e-6.
- The tolerance is stated and justified in section 4.

What makes this a reproduction and not a coincidence is listed in section 4. The near-degenerate
geometry means that "the printed point lies on both manifolds" is NOT, on its own, evidence.

## 1. The model as printed, and our mapping

- **Equations of motion.** SIAM Eq. 3.5 (= CMDA Eq. 5), implemented verbatim:
  H = (px² + py²)/2 + n(t)(px y − py x) − (1−µ)/r1 − µ/r2.
  - r1 = |(x + µρ, y)| and r2 = |(x − (1−µ)ρ, y)|.
  - ρ = 1 − e cos E, where E − e sin E = t.
  - n = √(1−e²)/ρ².
- **Coordinates.** Barycentric and rotating with the primaries' line. The frame is NOT pulsating.
  The variables are canonical momenta, px = ẋ − n y and py = ẏ + n x. Time t is the independent
  variable, and periapse is at t = 0.
- **Stroboscopic map.** The time-2π map started at t = 0 mod 2π, so θp = 0 (CMDA Sec. 3.1).
- **e = 0.0094** (SIAM Sec. 3.3).
- **µ = 2.5266448850435028e-5**, as printed in CNSNS Sec. 3.
  - With this µ, our corrector reproduces CNSNS Table 1. For the 3:4 orbit: ẏ to 5e-14, period to
    1e-11, multipliers to 2e-9, C = 3.0024. The 5:6 period agrees to 4e-9 and its multipliers to
    1.1e-7.
  - The other published Europa value, 2.5265115494603433e-5 (AAS 21-651), misses the periods by
    1.2e-5. It is rejected.
- **ω = 4π²/T1** (CMDA Sec. 4.9), where T1 is the parent orbit's period.
  - The Table 1 periods give ωu = 1.5580391955 and ωs = 1.0300114376. These round to the printed
    1.558039 and 1.030011.
  - SIAM's printed 5:6 period, 38.3281, matches Table 1 (38.328135). It does not match
    4π²/1.030011 = 38.32815.
  - The tori are therefore continuations of the Table 1 orbits at C = 3.0024. This is an
    inference, supported by four printed numbers.
- **Not pinned by the text.** Printed ω could be read literally, so we ran it as a second variant
  (`printed`). It is decisively worse (section 4).
- **Independent model check.** At e = 0.2, the Hamiltonian flow agrees with the separately derived
  non-pulsating Lagrangian EOM of `core/er3bp_paper_frame.py` to 1e-9 (test).

## 2. Gates (tests; identities and printed values only)

- **Integrator.** A numba DOP853 port. It matches scipy DOP853 to 1e-11.
- **Reduction at e = 0.**
  - Jacobi conserved to 1e-10 over 10 periods.
  - The vector field equals PCRTBP Eq. 3.2 term by term.
- **Map identities.**
  - Map inverse to 1e-11.
  - 2π-periodicity in t to 1e-12.
  - STM against central differences to 1e-7 relative.
  - Symplectic defect below 1e-10.
  - Time-reversal symmetry M F⁻¹ M = F with M = diag(1, −1, −1, 1), to 1e-11 (CMDA Sec. 5.3).
- **Table 1.** Both orbits reproduced, and the ω/period consistency shown.
- **Unperturbed circle.** Invariant and reversible.

## 3. Tori, bundles, conditioning

**Circles.** Continued in e at fixed ω, from 0 to 0.0094 in about 19 adaptive steps.

| Quantity | 3:4 | 5:6 |
|---|---|---|
| N | 511 | 1023 |
| Residual (on / off grid) | 9e-13 / 3e-11 | 1e-11 / 1e-10 |
| Reversibility defect | 3e-11 | 9e-9 |
| K(0) | (−1.401225, 0, 0, −0.775944) | (−1.240434, 0, 0, −0.852795) |

- K(0) is the perpendicular negative-x-axis crossing. Only y(K(0)) = 0 is imposed; px(K(0)) = 0
  comes out of the solve (6e-13).
- Re-solving at 2N+1 nodes moves the circles by 4e-9 (3:4) and 8e-9 (5:6). This is the circles'
  real accuracy.

**Bundles.** Computed by projective power iteration with constant-multiplier rescaling (CMDA
Sec. 4.8 and 4.11).

| Quantity | 3:4 | 5:6 |
|---|---|---|
| λu | 3.036678 | 2.988406 |
| λs | 0.329307 | 0.334627 |
| λu·λs | 1.000 | 1.000 |

- Off-grid invariance is 2e-9 to 3e-8.
- The reversor maps v_u to v_s to 4e-9 (3:4) and 3e-7 (5:6).
- Second-order terms are included. Their defect at offset 1e-4 is 7e-10 for 3:4 but 2e-8 for 5:6,
  so the charts used are 1e-5 and 1e-6.

**Comparisons with the paper's figures (read by eye).**
- Fig. 9 gives λu ≈ 3.03 at ωu = 1.558. Ours is 3.0367.
- Fig. 15 gives a normalised-column det Df of about 2.5e-6 at ωu = 1.558039. Ours is 2.53e-6.
  This quantity does not depend on the θ or s conventions.

## 4. The comparison

**A. Convention-free location.**
- X* is the printed point.
- The nearest point of our 3:4 unstable manifold is 4.97e-6 away. The nearest point of our 5:6
  stable manifold is 5.02e-6 away. Every component difference is below 5e-6.
- Both searches were seeded from X* alone.
- This is necessary but not discriminating. The manifolds nearly coincide along a curve (a
  "valley"), as they would in the circular problem. Any point in that valley passes this test.

**B. Refinement.**
- The 4×4 damped Gauss-Newton solve converges to a residual of 1.0e-9 at chart radius 1e-5 and
  3.1e-9 at 1e-6.
- The two junctions agree to 2.8e-9.
- The singular values of the normalised columns are (1.99, 0.19, 0.14, 4.8e-5).
- Re-seeding ±3e-3 along the near-null vector returns to the same point within 1e-8.
- At 2N+1 nodes (`table1_hi`) the point moves by at most 1.2e-7. Our x values range from
  −0.96063509 to −0.96063494. The last is the 1e-6 high-resolution run, which stalled at residual
  1.2e-8 and is reported, not dropped.

**C. What discriminates.**
- **Sign change at X\*.** Along the valley (θu fixed, the other three parameters minimised) the
  signed mismatch passes through zero at X*. It is 3e-9 there and grows linearly to 7.4e-7 at
  θu ± 0.003. So X* is the isolated connection, not merely a valley point.
- **The literal ω is wrong.** The `printed` variant (literal ω) puts the zero at
  (−0.96265, 0.88620, −0.51258, −0.64774), about 2e-3 away. ω offsets of 2e-7 and 4e-7 move the
  connection by 2e-3, so the five-decimal agreement is not generic.

**D. Tolerance.**
- Along the valley, the junction moves about 6e3 times the mismatch error (4.5e-3 per 7.4e-7).
- Our floor of 1e-9 to 1e-8 therefore allows 6e-6 to 6e-5 in principle. The measured movement
  across resolution and chart radius is 1.2e-7.
- The paper's stated refinement tolerance (Eq. 7.1 error below 1e-7, Sec. 7) allows its own
  printed point to sit up to about 6e-4 along the valley. The observed agreement, 5e-6, is about
  100 times tighter than the paper guarantees.

**E. Torus parameters.**
- Our 5:6 stable parameters, (θ, s) = (2.39705, −77.7330), match the printed FIRST pair,
  (2.39703, −77.73428).
- Our 3:4 unstable parameters, (1.83709, 203.4469), match the printed SECOND pair,
  (1.83093, 202.62277), after the θ origin is shifted by one mesh step 2π/1024. The shift implies
  rescaling s by |v_u(2π/1024)| = 0.99596, which is not a free parameter. That gives
  (1.83095, 202.6252).
- The printed tuple is labelled (θu, su, θs, ss), with u the 3:4 unstable manifold of Eq. 8.1. The
  values fit with the pairs the other way round. Figs. 12-13 show θu ≈ 2.39 and su ≈ −77 at this
  ω, so the paper is internally consistent and this is not a one-off typo. We cannot tell from the
  paper whether it is a convention or a slip.
- Our pre-registered orientation is v_x(0) > 0 for both bundles. Under this reading, both signs of
  s agree with the paper's.
- **Support for the one-step 3:4 origin:**
  - Fig. 11's x-component of the unit bundle at θ = 0 is about 0.1158. Ours is 0.1236 at the axis
    crossing and 0.1158 one step later.
  - The θ match to 2e-5.
- **Tensions in that reading:**
  - It puts K(0) at y ≈ +0.015. That conflicts with the text's "on the negative x-axis".
  - Our second-order x-coefficient there is 0.3014. Fig. 11 right reads about 0.3005.
- The 5:6 pair needs no shift.
- The s values are determined only to about 1e-3 along the near-null direction (they differ by
  5.6e-4 between our own two chart radii). The remaining s differences, 1.3e-3 and 2.5e-3, are of
  that order.

## 5. The #882 review's verification criteria applied (`verify_table1.json`)

1. **Contraction at the multiplier's rate for at least 8 periods.**
   - **Backward leg from the junction to 3:4: PASS.** 9 consecutive periods are within 2 percent
     of 1/λu (0.3288 to 0.3298 against 0.32931). The distance falls to 1e-7.
   - **Forward leg from the junction to 5:6: FAIL.** Only 4 periods are within 2 percent
     (0.3370 to 0.3397 against 0.33463). The distance then stalls at 2.8e-4.
   - Cause of the forward failure: the junction mismatch of about 3e-9 grows by λu(5:6) ≈ 3 per
     period. 3e-9 × 3^11 ≈ 5e-4, matching the observed stall.
   - Moving the junction only trades forward periods for backward ones. Passing both legs would need
     a junction mismatch of about 1e-13 or less, which the circles' measured accuracy (4e-9 to
     8e-9) rules out.
   - Supplementary only, not counted: the stable leg as computed by the parameterisation gives 11
     clean periods. That is the parameterisation's own leg, the weakness the review identified.
2. **Junction agreement between offsets 1e-5 and 1e-6: PASS** (2.8e-9, criterion 1e-5).
3. **Random 1e-6 perturbations must fail.** We perturbed the junction, not the departure state: a
   single 35-map integration from departure amplifies round-off by about 3^35, and even the true
   trajectory misses (0.17).
   - **Backward: PASS.** True 9.9e-8; perturbed, at least 0.156.
   - **Forward: INCONCLUSIVE.** True 2.8e-4; perturbed, at least 7.1e-4.

## 6. Not done, deviations, notes for the coordinator

- **Step 6 (regenerating the #882 junctions and applying the criteria there) was NOT done.**
  Criterion 1 failed on the forward leg, so the precondition "steps 1 to 5 succeed" is not met.
  - Inference, not run: whether criterion 1 can be met depends on the per-period ratio λu/λs, here
    about 9, and on the junction floor. #882's R3b passed with a ratio of about 2.7. Criterion 1
    as written may be unsatisfiable in double precision for the #882 Jupiter branches too. The
    criterion should be restated before it is adopted as a gate.
- **Corrector bordering.** The fixed-ω circle corrector needs a bordering column J⁻¹DK. The
  fixed-ω invariance operator has a one-dimensional cokernel (the averaged action equation), which
  DK does not span. Bordering with DK alone gave only linear convergence.
- **Paper workflow not replicated.** The GPU mesh search, layers and Whitney continuation were not
  replicated. The connection was located from the printed point by the convention-free search
  above.
- **Reused from #882:** the odd-node Fourier helpers and the PCRTBP symmetric-orbit corrector.
  **Rewritten:** every integrating function, because the #882 interfaces take a `CCR4BPSystem`
  (four-body, velocity coordinates, constant rotation).
- **No new default-suite test for the connection itself.** Building the tori takes minutes. The
  stored outputs are not used as goldens.
