# Digest: Langebartel 1965, "Two-Center Problem Orbits as Intermediate Orbits for the Restricted Three-Body Problem" (NASA TN D-2939) (#960, #948)

R. G. Langebartel (NASA Goddard Space Flight Center, Greenbelt, Md.), "Two-Center Problem Orbits as Intermediate
Orbits for the Restricted Three-Body Problem", **NASA Technical Note D-2939**, NASA, Washington, D.C., September 1965.
NTRS 19650022582 (NTRS record: report number NASA-TN-D-2939; publication date 1965-09-01). Printed "NASA-Langley,
1965 G-636" (p.19). Price $1.00. No DOI (Crossref: no record for the TN; the nearest hit is the author's later journal
paper, see sec. 5).
- TN number read on the cover image (p.i: "NASA TN D-2939", also on the spine strip) and confirmed by the NTRS record.
- Source file: upload `df6b5648-19650022582.pdf`, 25 pp. (cover, title page, summary, contents, blanks, text pp.1-17,
  references p.18, symbols p.19), md5 25f951620e87349dd2f94d9de5f5dd84. AFRL "VSIL Digitizing Team" scan, 1-bit
  CCITT images at 300 dpi, with an old OCR text layer. That layer is usable for prose but breaks words ("t e r m s",
  "Langebmtel", "ORB6TS") and garbles every equation.
- OCR copy (the file to be filed): `langebartel-1965-ocr.pdf` in this folder (`ocrmypdf --force-ocr -l eng`,
  Ghostscript PDF/A, 25 pp.), md5 6f3e5d2ab75406cc1f80bae794337ca5. The prose text is clean (single-letter fragments on
  text p.1: 40 in the old layer, 5 in the new). Equations are still unusable as text; read them on the image.
  - Render check, pp. 1, 9 and 18 at 100 dpi: same pixel size (p.1 is one pixel taller in the copy; compared on the
    common area). Mean grey difference 2.7, 5.0 and 5.1 of 255. A 250-dpi look at the eq. (27) block in the copy is sharp.
  - Size: 14.3 MB against 0.70 MB for the original. `--force-ocr` re-rasterises the 1-bit pages to 400-dpi RGB images.
    A second run with `--output-type pdf` gave 15.9 MB, so the size comes from the rasterising, not from PDF/A. If the
    size matters more than a clean prose layer, the original file is also fit to file (its equations need the image anyway).
- Proposed filename: `cyclers_pdf/papers/langebartel-1965-two-center-problem-orbits-intermediate-orbits-restricted-three-body-problem-nasa-tn-d-2939-ntrs-19650022582.pdf`.
- How I read it: the cover, title and summary pages and every text page (pp.iii, 1-19) on 130-dpi page images, with
  zooms at 250-300 dpi. I checked the equations
  that matter with sympy and scipy (`scripts/check_langebartel.py`, output `scripts/check_langebartel.out`).

## 0. Verdict

**An analytic-method note. It contains NO numbers: no table, no figure, no orbit, and no numerical comparison with the
restricted problem.** The task asked me to transcribe the numerical comparisons. There are none to transcribe. The
concluding remarks (p.17) say this directly: "As far as application to the Apollo project is concerned, work remains to
be done to ascertain the probable range of values for r_i, rho_i for these orbits, so that it can be determined whether
the type of approximation used ... will be reasonable." The perturbation equations are "given in outline" (Summary,
p.iii), and the author does "not carry out the details" (p.15). So there is no `-tables.yaml` and no witness table for
this paper.

What it is:
- It takes Charlier's idea (1902) of using Euler two-fixed-centre orbits as intermediate orbits for the circular
  restricted problem (CR3BP). It writes the CR3BP Hamiltonian in elliptic coordinates (xi, eta) with foci at the two
  primaries.
- It splits H into a separable part H0 and a perturbation H1. H0 keeps both gravity terms. H1 holds the
  rotating-frame (Coriolis and centrifugal) terms that are linear in the momenta. H0 also contains two free functions
  u(xi) and v(eta), and their choice moves part of H1 into H0.
- It solves H0 by Hamilton-Jacobi in Jacobi elliptic functions. It gives a small-modulus expansion in elementary
  functions for slender "lemniscate" orbits that run from near one primary to near the other (the Apollo-type orbits).
  It also writes the variation-of-constants equations, but does not solve them.

What it gives the project:
- **For `#948` (both-primary regularised corrector): the coordinate chart, not the perturbation theory.**
  - The elliptic coordinates are the Thiele-Burrau map, x + iy = c + (1/2) cos(xi - i eta). This is COMPUTED by us; the
    paper never uses the word "regularisation".
  - Near each primary the map is a Levi-Civita squaring: z - z_E ~ (w - pi)^2/4 and z - z_M ~ -w^2/4.
  - Eq. (2) is the full CR3BP Hamiltonian in this chart, including the rotating-frame linear terms. I checked it
    against the barycentric Lagrangian at 20 random points: the maximum difference is 3e-15. A control with the linear
    terms sign-flipped gives 0.60, so the check can fail.
  - Multiplying by D = cosh^2 eta - cos^2 xi = 4 r1 r2 removes both 1/r singularities (COMPUTED). The paper's own time
    law, dt = (1/4)(cosh^2 lambda2 - cos^2 lambda1) d(phi1) (p.11), is exactly dt = r1 r2 d(tau). That is the Sundman
    factor that `core/cr3bp_regularized.py` already uses, but without the coordinate change.
  - So this TN supplies, in one planar formula, the missing half of `#948`: one chart that is regular at both primaries
    at once, with no switching between charts. Proposal in sec. 4.
- The printed partition has a factor-2 inconsistency (sec. 2): the printed H0~ + H1~ does not add back to H. H0~ and
  the intermediate orbit are right; the outlined perturbation equations rest on a wrong H1~.
- The perturbation scheme itself is not worth reusing. Samter (1922) found the perturbations "relatively large" (p.1).
  The terms put into H1 are the Coriolis and centrifugal terms, which are not small at Earth-Moon distances. The author
  stops before any estimate of their size.

## 1. Formulation (pp.1-4, read on the page images)

- **Frame and units (p.1).** "Let the two finite masses ... be indicated by mu and 1 - mu." Circular orbits, barycentric
  rotating coordinates, "the particle with mass mu are (1 - mu, 0) and of the particle with mass 1 - mu are (-mu, 0)".
  Unit distance and unit angular rate are implied by the form of L. This is the same layout as the project frame (Earth
  at -mu, Moon at 1 - mu).
  - L = (1/2)(xdot^2 + ydot^2) + (x ydot - y xdot) + U(x, y)
  - U = (1/2)(x^2 + y^2) + (1 - mu)/|(x + mu)^2 + y^2|^(1/2) + mu/|(x - 1 + mu)^2 + y^2|^(1/2)
  - Cited to Wintner (Reference 1, p.350).
- **Bipolar coordinates (p.2):** r1 and r2 are the distances from (-mu, 0) and (1 - mu, 0). The Hamiltonian in (r1, r2)
  is eq. (1). I did not check eq. (1); it is not needed downstream.
- **Elliptic coordinates (p.2):** "cos xi = r1 - r2, cosh eta = r1 + r2". So
  x = -mu + 1/2 + (1/2) cos xi cosh eta and y = (1/2) sin xi sinh eta.
  - Hence r1 = (cosh eta + cos xi)/2 and r2 = (cosh eta - cos xi)/2 (COMPUTED, sympy, exact).
  - Earth (mass 1 - mu) is at xi = pi, eta = 0. The Moon is at xi = 0, eta = 0.
- **Eq. (2), the CR3BP Hamiltonian in elliptic coordinates.** p1 and p2 are conjugate to xi and eta.
  H = 2 p1^2/(cosh^2 eta - cos^2 xi) + 2 p2^2/(cosh^2 eta - cos^2 xi)
      - sinh eta [cosh eta + (1 - 2 mu) cos xi] p1/(cosh^2 eta - cos^2 xi)
      - sin xi [(1 - 2 mu) cosh eta + cos xi] p2/(cosh^2 eta - cos^2 xi)
      - 2(1 - mu)/(cosh eta + cos xi) - 2 mu/(cosh eta - cos xi)
  - COMPUTED check: a point transformation of H = (1/2)|p|^2 + (y px - x py) - (1 - mu)/r1 - mu/r2 agrees with eq. (2)
    to 3e-15 at 20 random points.
  - There is no centrifugal term in H. With canonical momenta it cancels; only the terms linear in p carry the rotation.
- **Eq. (3), the "two-center" problem (p.3):** L^(2), U^(2) and H^(2) are eq. (2) without the terms linear in p.
  - Author's caveat (p.3): eqs. 3a-3c "are not strictly the formulae for the two fixed centers problem if we adhere to
    our definition of x, y since these represent coordinates in a rotating system whereas in the Newtonian two-center
    problem the two masses are fixed in an inertial frame".
  - Also: "the p_i represent different quantities in H and H^(2)".
- **Hamilton-Jacobi equation (eq. 4, p.3).** "The variables xi, eta do not separate in this partial differential
  equation although they do in Hamilton's Equation formed from H^(2)."

## 2. Partition, separation and the intermediate orbit (pp.3-12)

- **Partition (eq. 5).** H = H0 + H1.
  - H0 = 2(cosh^2 eta - cos^2 xi)^-1 (p1^2 + p2^2) - (cosh^2 eta - cos^2 xi)^-1 [u(xi) p1 + v(eta) p2]
    - 2(1 - mu)/(cosh eta + cos xi) - 2 mu/(cosh eta - cos xi).
  - H1 = -(cosh^2 eta - cos^2 xi)^-1 {(sinh eta [cosh eta + (1 - 2 mu) cos xi] - u(xi)) p1
    + (sin xi [(1 - 2 mu) cosh eta + cos xi] - v(eta)) p2}.
  - "The functions u(xi), v(eta) are at our disposal. If they are taken as zero then H0 is formally the same as the
    two-center Hamiltonian ... the more general H0 ... is still separable" (p.4).
- **Contact transformation (eq. 6):** xi = lambda1, eta = lambda2, p1 = Lambda1 + (1/2) u(lambda1),
  p2 = Lambda2 + (1/2) v(lambda2). This gives eq. (7):
  - H0~ = 2(cosh^2 lambda2 - cos^2 lambda1)^-1 [Lambda1^2 + Lambda2^2 - (1/4)u^2(lambda1) + (1 - 2 mu) cos lambda1
    - (1/4)v^2(lambda2) - cosh lambda2]
  - H1~ = -(cosh^2 lambda2 - cos^2 lambda1)^-1 ({sinh lambda2 [cosh lambda2 + (1 - 2 mu) cos lambda1] - u(lambda1)}
    {Lambda1 + (1/2)u(lambda1)} + {sin lambda1 [(1 - 2 mu) cosh lambda2 + cos lambda1] - v(lambda2)} {Lambda2 + (1/2)v(lambda2)}),
    that is, eq. (5)'s H1 with p_i replaced as in eq. (6). See the inconsistency note below.
  - **Factor-2 inconsistency between eqs. (5), (6) and (7) (COMPUTED; read on 300-dpi zooms of eqs. 5 and 7 and on
    the p.13 restatement; no hidden "2" in the print).**
    - Substituting eq. (6) into eq. (2) and subtracting the printed H0~ + H1~ leaves
      D^-1 [u Lambda1 + v Lambda2 + (u^2 + v^2)/2], with D = cosh^2 lambda2 - cos^2 lambda1 (`check_langebartel.out`,
      item 9). So the printed H0~ + H1~ does not sum back to H.
    - The residual is zero if both H0 and H1 in eq. (5) carry 2u and 2v (that is, "-2D^-1[u p1 + v p2]" in H0 and
      "- 2u", "- 2v" in H1, and the same in H1~). With eq. (5) as printed, the shift in eq. (6) would have to be u/4.
      No single relabelling of u makes eqs. (5), (6) and both parts of (7) consistent as printed.
    - Consequence: the separation and the intermediate orbit (eqs. 8-27) use only H0~, and H0~ is the correct
      separable part for the shift in eq. (6). They are not affected. The perturbation equations (28-41) are built on
      the printed H1~, which is not the true remainder; it needs "- 2u" and "- 2v". Since the perturbation equations are
      only outlined, this is a defect in the outline, not in any number.
- **Separation (eqs. 8-9).** The complete integral is
  W0~ = -h t + Integral^lambda1 ±[(1/4)u^2 - (1 - 2 mu) cos lambda - (1/2) h cos^2 lambda - alpha]^(1/2) d lambda
        + Integral^lambda2 ±[(1/4)v^2 + cosh lambda + (1/2) h cosh^2 lambda + alpha]^(1/2) d lambda.
  "h and alpha are canonical constants." Eq. (10): t - beta1 = dW0~/dh, -beta2 = dW0~/d alpha, and Lambda_i = dW0~/d lambda_i.
- **Which quantities are conserved.**
  - For the intermediate orbit (H0~ alone): h (the value of H0~, an energy-like constant) and alpha (the separation
    constant, the analogue of the Euler problem's second integral). beta1 is the time of passage through the reference
    point, and beta2 is a phase.
  - For the full CR3BP: only H itself is conserved; it is the Jacobi integral in Hamiltonian form, with C = -2H in the
    project convention. h and alpha drift under H1 according to eq. (28).
  - The paper does not state the Jacobi relation. The C = -2H identity is ours.
- **Quadratic choice of u and v (eq. 11).** (1/4)u^2 = u0 + u1 cos lambda + u2 cos^2 lambda and
  (1/4)v^2 = v0 + v1 cosh lambda + v2 cosh^2 lambda, which "keeps the orbits within the general types that appear in
  the two-center problem".
  - Eq. (12) folds these into mu_0 = u0 - alpha, mu_1 = u1 - (1 - 2 mu), mu_2 = u2 - (1/2)h, nu_0 = v0 + alpha,
    nu_1 = v1 + 1 and nu_2 = v2 + (1/2)h.
  - The orbit equations are eqs. (13a), (13b), (14a) and (14b). The motion depends on the roots rho_1, rho_2 of
    mu_2 x^2 + mu_1 x + mu_0 and r_1, r_2 of nu_2 x^2 + nu_1 x + nu_0 (eq. 15).
- **The Apollo-type case (p.6).** The text reads: "A proposed orbit for the Apollo capsule is frequently described ...
  as being first an arc of a two-body ellipse extending from near the earth to the neutral point between the earth and
  the moon and then continuing on another two-body ellipse to near the moon, thereby describing a very slender figure ...
  these are the cases Ib alpha and Ib beta of Charlier (Reference 2, Vol. I, p. 126) which represent quasi-lemniscates
  winding about the two mass points."
  - Assumptions: mu_2 > 0 and nu_2 < 0, "the case, for example, if u2, v2 are zero or are quite small, and h < 0".
  - Further: r1 > 1 > r2, and rho_1, rho_2 either complex or real with rho_1 >= rho_2 > 1.
  - Then 0 <= lambda1 <= 2 pi and 0 <= cosh lambda2 <= r1. The orbit fills the region inside the confocal ellipse
    cosh eta = r1.
- **Reference point (p.7).** lambda1(0) = 0 and lambda2(0) = -r1^ (with cosh r1^ = r1). This is the axis point
  x = -mu + 1/2 + r1/2, y = 0, "on the outside of the point with mass mu". beta1 is the time of passage there.
- **Solution in elliptic functions (eqs. 16-19).** Parameters phi_1 and phi_2 = phi_1 - beta2 (eq. 16) give eq. (17a) for
  cos lambda1 (sn with modulus k1 = [2(rho1 - rho2)/((rho1 - 1)(rho2 + 1))]^(1/2) and argument a1 phi_1,
  a1 = (1/2)[mu_2(rho1 - 1)(rho2 + 1)]^(1/2)). Eq. (17b) gives cosh lambda2, with k2 = [(r1 - 1)(r2 + 1)/(2(r1 - r2))]^(1/2)
  and a2 = (1/2)[-2 nu_2 (r1 - r2)]^(1/2). Eq. (18) gives the momenta, and eq. (19) gives half-angle forms without sn^2.
  - COMPUTED checks against direct quadrature of eq. (16), with scipy's m = k^2:
    - eq. (17a): 3e-15;
    - eq. (19a): 2e-15;
    - eq. (18a): 3e-15;
    - eq. (17b), first form: 5e-13.
  - With m = k instead the error is 0.42. So the paper's "parameter k" is the modulus.
- **Degenerate and small-modulus cases (p.8).** "letting r1 -> 1 gives simply cosh lambda2 = 1 which represents the
  straight line solution running from one of the mass points to the other". "k1 is small if rho1 is close to rho2."
  Small k1 and k2 give "lemniscate-type orbits", expanded in elementary functions.
  - Eq. (20): Fourier series of sn, cn, dn. Eq. (21): the nome q. Eq. (22): K.
  - Eq. (23): sn 2Ku = sin(pi u) + ((1/16) sin(pi u) + (1/16) sin(3 pi u)) k^2 + O(k^4), with similar series for cn and dn.
    COMPUTED: the error is 7.7e-5 at k = 0.2, so the printed signs are right.
  - Eqs. (24a-f) give cos lambda1, cosh lambda2, Lambda1 and Lambda2 to O(k^2).
  - Eq. (25) is a "convergence improvement" form with O(k^8) error. The resulting formulae "are not written out here".
- **Time law (p.11, eq. 26).** dt = (1/4)(cosh^2 lambda2 - cos^2 lambda1) d phi_1. Eq. (27) gives 4(t - beta1) to O(k^2).
  - The author warns (p.13): eq. (27) "is not useful if k2 is small because it is r1 - 1 that is small". The remedy is to
    re-express eq. (24b) in k2 and phi_2 first.

## 3. Perturbation equations (pp.13-17)

- The variation equations use H1~ (eq. 7), restated on p.13 with the same coefficients (so the factor-2 defect of sec. 2
  carries into eqs. 28-41). Sine forms of sin lambda1 and sinh lambda2 are given, and
  u and v are written in terms of phi.
  - "the radicands in u(lambda1) and v(lambda2) are perfect squares if we take u1^2 - 4u0u2 = v1^2 - 4v0v2 = 0".
- Eq. (28): dh/dt = dH1^/d beta1, d alpha/dt = dH1^/d beta2, d beta1/dt = -dH1^/dh and d beta2/dt = -dH1^/d alpha.
- Eqs. (29)-(32) change the independent variable from t to phi_1.
- Eqs. (33)-(41) are the auxiliary partial derivatives of mu_i, nu_i, rho_i, r_i, a_i, k_i and of the phi-t map.
- "We shall not carry out the details of determining explicitly the coefficients in Equation 31 for our particular
  case" (p.15).
- **Concluding remarks (p.17):**
  - "The formulae of the preceding section are lengthy and involved."
  - "For a first order theory it may well be better to deal with the variation equations in t of Equation 28, rather
    than those in phi_1 of Equation 31."
  - Six arbitrary constants remain in u and v, and "an investigation could well be made concerning the feasibility and
    usefulness of some other choice for these functions."

## 4. PROPOSAL for `#948` (not done; for the lead)

`#948` needs a corrector that is regular at both primaries. The `#970` note (sec. 2-3) found that core/ has only the
Sundman time factor (`core/cr3bp_regularized.py`, force law not regularised) and one-centre KS (`core/cr3bp_ks.py`,
Moon only). Its suggested route is a second, Earth-centred KS chart with switching. Langebartel's TN offers a simpler
planar alternative.

1. **Planar global chart (Thiele-Burrau).** Integrate in (xi, eta, p1, p2) with the Hamiltonian eq. (2). Use the new
   time tau, with dt = (cosh^2 eta - cos^2 xi) d tau / 4 (= r1 r2 d tau; the 1/4 is the paper's eq. 26 scaling).
   - Work on the fixed-energy surface K = (D/4)(H - h) with D = cosh^2 eta - cos^2 xi. K is a polynomial in p1 and p2.
     Its coefficients are trigonometric and hyperbolic functions, and it has no 1/r1 or 1/r2 term.
   - The explicit form of D(H - h) is in `scripts/check_langebartel.out`, item 3.
   - Hamilton's equations for K are smooth through both collisions. A Newton corrector on symmetric crossings can work
     in these variables, with the STM from the variational equations of K.
2. **Symmetric-crossing conditions are simple in this chart.** The x-axis (y = 0) is sin xi sinh eta = 0. It consists of
   the segment between the primaries (eta = 0), the ray beyond the Moon (xi = 0) and the ray beyond the Earth (xi = pi).
   Perpendicular crossings become conditions on one momentum. This suits the `#997` F1-F4 continuations, whose floor
   stops (F1 Earth-surface impact at C = 1.09, F3 at C = 1.928, F4 3/7-r at C = 1.355) are where the plain integrator
   stops.
3. **Controls.**
   - (a) Langebartel's own exact limit: r1 -> 1 gives cosh lambda2 = 1, "the straight line solution running from one of
     the mass points to the other" (p.8). This is an intermediate-orbit statement (two fixed centres), not a CR3BP
     orbit. Use it to check the chart, not the corrector.
   - (b) The `#970` Schwaniger table (C, T, perigee 6555 km, periselene 2202.5 km, b_h 525.3, b_v -182.4) must be
     reproduced to its stated tolerances.
   - (c) Consistency between charts: the same orbit propagated in (xi, eta), in physical time and in Moon-KS must agree,
     as in the `#970` three-integrator closure.
4. **Spatial case.** Langebartel is planar only. For 3-D, the held Waldvogel 1967 (spatial Birkhoff) is the matching
   source. Szebehely 1967 (held) sec. 3.6 gives the Thiele-Burrau transformation and the regularised equations in the
   standard form. Use it as the primary formula source and Langebartel eq. (2) as a cross-check of the rotating-frame
   terms.
5. **Do not reuse** the H0/H1 perturbation scheme (eqs. 5-41) for `#948`. It was never carried to numbers, its printed
   H1~ is inconsistent with H by a factor 2 in u and v, Samter found the perturbations large, and a numerical corrector
   needs only the chart and the time law (eq. 2, which we verified).

## 5. Citation mining (4 references)

Checked with `ls cyclers_pdf/papers | grep -i` and `docs/notes/CORPUS_INDEX.md`, and against the wanted list.

| Ref | Item | Status |
|---|---|---|
| 1 | Wintner, A. (1947), The Analytical Foundations of Celestial Mechanics, Princeton UP (p.145 Euler problem; p.350 CR3BP Lagrangian) | not held; textbook; not on the wanted list |
| 2 | Charlier, C. V. L. (1902), Die Mechanik des Himmels, Leipzig: Veit (Vol. I p.126 cases Ib alpha/beta; Vol. II two-centre intermediate orbits) | not held; textbook; not on the wanted list |
| 3 | Samter, H. (1922), "Das Zweizentren-Problem in der Stoerungstheorie", Astron. Nachr. 217(5193-5194):129-152 | not held; **doi 10.1002/asna.19222170902** (Crossref: Astronomische Nachrichten 217(9-10):129-152, 1922); not on the wanted list. Low priority: near-equal-mass case, and the TN reports its outcome |
| 4 | Magnus, W. & Oberhettinger, F. (1948), Formeln und Saetze fuer die speziellen Funktionen der mathematischen Physik, Springer | not held; handbook (elliptic-function identities) |

Forward link found on Crossref: **Langebartel, R. G. (1981), "Formulation of the restricted three-body problem with
two-center problem orbits as intermediate orbits", Astrophysics and Space Science 75(2):437-454, doi
10.1007/bf00648654.** It is not held and not on the wanted list. It is the journal sequel of this TN, and it may carry
the size estimates that the TN leaves open. PROPOSAL: add it to the wanted list at low-to-medium priority, for `#948`
only if the chart in sec. 4 needs a second source.

No reference is a missed cycler source. The TN names "orbits of the type suggested for the Apollo mission" (p.1) but
gives none.

*Filed as `cyclers_pdf/papers/langebartel-1965-two-center-problem-orbits-intermediate-orbits-restricted-three-body-problem-nasa-tn-d-2939-ntrs-19650022582.pdf`.*
