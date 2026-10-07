# Digest: Perko 1983, "Periodic Solutions of the Restricted Problem that are Analytic Continuations of Periodic Solutions of Hill's Problem for Small mu > 0" (#960 batch 36)

L. M. Perko (Northern Arizona University), "Periodic Solutions of the Restricted Problem that are Analytic
Continuations of Periodic Solutions of Hill's Problem for Small mu > 0", Celestial Mechanics 30(2):115-132
(1983), doi 10.1007/BF01234301 (Crossref: Celestial Mechanics 30(2), pp.115-132, 1983-06; the Crossref title
has a broken mu glyph). Received July 1982, accepted October 1982. The title on the page image (p.115) matches.
- Supplied file `7aae3f9c-perko1983.pdf`, 18 pp. = journal pp.115-132 (PDF page n = p.114 + n). The first page
  is the article's title page (title block, abstract, sec. 1), not a separate cover. md5
  ac14495ab1d6380ba659089392c597e6. Springer scan with a text layer (symbols garbled).
- **Proposed corpus filename:**
  `cyclers_pdf/papers/perko-1983-periodic-solutions-restricted-problem-analytic-continuations-hill-problem-small-mu-celest-mech-30-115-doi-10.1007-BF01234301.pdf`
- **How I read it:** the whole text from the text layer; on 150 dpi page images: pp.115, 117 (scaling,
  C = 3 + eps^2 Gamma), 124 (eqs 17-19, Theorem 3), 129-131 (Lemma 2 end, Existence Theorem, Broucke
  comparison, Fig. 1); the Broucke conversion formula at 300 dpi. Checks: `check_perko.py` -> `check_perko.out`
  (secs E, F), `check_perko_C1_folds.py` -> `.out`.
- Wanted list: **row 21**. Combined scope note: `theorem-scope.md`.

## 0. Verdict

**A local continuation theorem: a symmetric periodic Hill orbit with no ejection and a non-zero
non-degeneracy number C1 continues, for mu small enough, to a one-parameter family of symmetric periodic
orbits of the planar circular restricted problem near the small primary ("third species", Henon's term).
No numbers for mu; no tables.**
- Existence Theorem (p.129-130, image): "In any compact subset of the (x0(0), y0'(0))-plane, there are at most
  a finite number of points on any one of the Hill's families a, c, f, g, g', g'', ... for which
  C1(x0(0), y0'(0)), defined by (17), is equal to zero. Given (x0(0), y0'(0)) defining a symmetric periodic
  solution r0(t) of Hill's problem (3) for which C1 != 0 and r0(t) > 0, and given k0 > 0, there exists an
  eps0 > 0 such that for all eps in (0, eps0) and abs(eta0') < k0 there is a unique one-parameter family of
  symmetric periodic solutions of the restricted problem (2), r(t, eps), with parameter eta0' and initial
  conditions x'(0, eps) = y(0, eps) = 0, x(0, eps) = x0(0) + eps(K0 - C4 eta0')/C1 + O(eps^2) and
  y'(0, eps) = y0'(0) + eps eta0'".
  - eps = mu^(1/3); C = 3 + eps^2 Gamma (p.117). eps0 depends on the orbit. No bound on mu anywhere.
  - Symmetric only. "Bifurcation points" (C1^2 + C4^2 = 0, p.129) and ejection orbits (r0(t) = 0) are excluded.
    The paper says regularisation would remove the ejection restriction "however, this will not be carried
    out" (p.130).
  - The family-level clause (finitely many C1 = 0 points) is proved only for a, c, f, g, g', g'', ... and only via
    their ends (Lemma 1: f and g for small period, Hill-Wintner series; Lemma 2: a, c, g', g'', ... for
    Gamma0 << -1, the asymptotics of paper II), plus analyticity.
  - The pointwise part (Theorems 1-3) assumes only "a symmetric periodic solution of Hill's problem ... r0(t) > 0"
    and C1 != 0. So it applies to any symmetric Hill orbit, including Henon 2003's H-families, at each member
    with C1 != 0. Perko does not check this for them (he could not: they were published in 2003).
- Closing claim (p.131): "each symmetric periodic solution of Hill's problem (3) that is not common to two or more
  Hill's families (and does not have an ejection at the origin) can be analytically continued to a local
  one-parameter family"; near a bifurcation the restricted-problem characteristic "will leave one Hill's
  characteristic to parallel another" in an O(eps^(1/2)) neighbourhood (expected, not proved; Broucke 1968,
  mu = 0.012155, Fig. 1).
- **What it gives the project:** the formal basis for "a Hill family seeds a CR3BP family at small mu" and the
  exact scaling map C = 3 + mu^(2/3) Gamma, x_CR3BP = mu^(1/3) x_Hill (origin at the small mass, large mass at
  (-1, 0)). It does not certify any particular mu.
- **Error found:** the Broucke conversion on p.130 prints Gamma = (-2E - 3)/mu^(1/3) + mu^(1/3)(4 + x0^3) +
  O(mu^(2/3)) (300 dpi image). By the paper's own p.117 (C = 3 + eps^2 Gamma, eps = mu^(1/3)) the leading
  denominator must be mu^(2/3) if C = -2E. I did not resolve the O(mu^(1/3)) term (it depends on Broucke's
  E convention and on which x0 is meant). Use Gamma = (C - 3)/mu^(2/3) from p.117.
- **Catalogue implication (PROPOSAL only):** none.

## 1. Content

- **Sec. 1 (pp.115-116).** Hill 1878, Wintner 1926 (f, g for Gamma0 >> 1), Perko I-II (a, c, f, g, g', g'', ...
  for Gamma0 << -1), Moulton 1906, Perron 1937, Siegel 1950-51 (small direct and retrograde orbits in the
  restricted problem). Names: first species (continuations of Kepler motions), second species (angular
  limits as mu -> 0), third species (continuations of Hill orbits), "at the suggestion of M. Henon".
- **Sec. 2 (pp.116-131).** Restricted problem (1) centred at the small mass mu, large mass at (-1, 0); U written
  out and expanded (b_n from (1 + xi)^(-1/2)). Scaling x -> x/eps gives (2): Hill's equations plus eps grad V,
  with lim V = -2 - x^3 + (3/2) x y^2 (p.117, image). Jacobi integral Gamma = -(x'^2 + y'^2) + 3x^2 + 2/r +
  2 eps V(x, y, eps), C = 3 + eps^2 Gamma.
  - First-order perturbation r = r0 + eps r1 + R2 about a symmetric Hill orbit r0 (eq. 7); variational matrix
    g_r; forcing f(r, 0) = (3(y^2 - 2x^2)/2, 3xy) (from V's cubic term).
  - Theorem 1: error bound O(eps^2) on [0, T0] (needs r0 > 0). Theorem 2: the next x-axis crossing time t*.
    Theorem 3 (p.124): if C1 != 0, a unique xi*(eta0', eps) with x'(t*) = 0, and C1 xi* + C4 eta0' = K0 + O(eps).
  - Definitions (17)-(18) (p.124, image): C_i = x0''(t1) a2i - y0'(t1) a3i, i = 1..4, and
    K0 = y0'(t1) int_0^t1 i . Phi_vv(t1, t) f0(t) dt - x0''(t1) int_0^t1 j . Phi_rv(t1, t) f0(t) dt, with
    [aij] = Phi(t1, 0), t1 = T0/2.
  - Lemma 1 (pp.126-127, image): C1 = (21 pi/4) m^(2/3)[1 + O(m)] > 0 on f and g for small period (Hill's
    Fourier series, abs(m) = T0/(2 pi)). The line before it gives (x0''/y0') a21 - a31 = -(21 pi/4) m[1 + O(m)],
    and y0'(0) = (a0/m)(1 + ...) with a0 = m^(2/3)(...) > 0, which multiplies out to -(21 pi/4) a0: the printed
    "> 0" looks like a sign slip. Only C1 != 0 is used, so the lemma stands.
  - Lemma 2 (pp.127-129): C1 = O(Gamma0^2) != 0 on a, c, g', g'', ... for Gamma0 << -1.
  - Existence Theorem (above); Fig. 1 overlays Broucke's mu = 0.012155 families I, G, C, H1, H2 (from a, c, f,
    g, g') on Henon 1969's Fig. 1 characteristics.

## 2. Checks

- **C1 along the families (`check_perko.out` sec. E, `check_perko_C1_folds.out`).** I computed C1 and C4 from
  eq. (17) with the Hill state-transition matrix (DOP853, rtol 1e-12), t1 = the time of Henon's N-th x-axis
  crossing, state order (x, y, x', y'):

  | Orbit (source) | Gamma | N | C1 | C4 |
  |---|---|---|---|---|
  | f (Henon 1969 Table 3) | -10 | 1 | 60.2 | 30.6 |
  | a (Henon 1969 Table 8) | -4 | 1 | 1.64e4 | 8.85e3 |
  | g3 (Henon 2003 XIII) | 1 | 3 | -10.25 | 2.11 |
  | g3, maximum of Gamma | 3.806201 | 3 | 32.4 | 35.9 |
  | g3 | 2 | 3 | 3.83e4 | 6.28e3 |
  | Hg (Henon 2003 X) | 1 | 4 | 10.28 | -6.36 |
  | Hg | 2 | 4 | 145.8 | -32.6 |
  | Hg, maximum of Gamma | 3.836201 | 3 | -4.53e4 | -5.03e4 |
  | Hg = f described four times | 0.823630 | 4 | 2e-6 (0 to print precision) | 2e-5 (0) |

  - C1 is non-zero at every interior member I tried, including the folds (maximum Gamma), where the family is
    parametrised by eta0' at fixed eps and the fold is no obstacle.
  - At Hg's junction with f described four times (Gamma = 0.823630) both C1 and C4 vanish to the precision of
    the printed initial condition. That is Perko's "bifurcation point" and is excluded. It is the lower end of
    Henon's stable interval of Hg [0.823630, 1.914252]. The interior of that interval contains no ejection orbit
    (Henon 2003: Hg ejections at 3.78760, 2.433055 and 3.71730; g3's double vertical ejection at 2.238611), so
    Perko's pointwise theorem applies to the interior stable Hg members (I checked Gamma = 1 only).
  - C1 has opposite signs at g3 Gamma = 1 (first branch) and Gamma = 2 (second branch). Between them lie the
    maximum of Gamma and the double vertical ejection at Gamma = 2.238611, where Phi blows up, so C1 can change
    sign through infinity. A zero of C1 there is not established.
- **Broucke conversion (sec. F):** see the Verdict. mu^(1/3) = 0.22992, mu^(2/3) = 0.05287 at mu = 0.012155: the
  two forms differ by a factor of 4.3.

## 3. What it gives R15, R16, #974, #945

- **R15, asymmetric Henon 2003 families:** does not cover (symmetric only; asymmetric not mentioned).
- **R15, Hg and the symmetric H-families at Europa / Ganymede / Titania:** background only. The theorem says
  each non-ejection member with C1 != 0 continues for mu < mu0(orbit), with the CR3BP characteristic O(mu^(1/3))
  from the Hill one. It gives no mu0, so it does not certify mu = 2.5e-5 (Europa), 7.8e-5 (Ganymede) or
  Titania; the CR3BP correction R15 plans is still needed. It does support seeding with C = 3 + mu^(2/3) Gamma,
  x = mu^(1/3) xi (p.117). That map is leading order only: the Jacobi integral of (2) carries 2 eps V with
  V(x, y, 0) = -2 - x^3 + (3/2) x y^2 (p.117, image), so for a CR3BP state
  C = 3 + mu^(2/3)[Gamma_Hill(state) - mu^(1/3)(4 + 2x^3 - 3xy^2) + O(mu^(2/3))]. The constant part 4 mu^(1/3) is
  0.12 at Europa, 0.17 at Ganymede and 0.25 at Titan in Gamma units, i.e. 11-23 % of the width (1.09) of Hg's
  stable interval. R15 should treat the scaled ends of that interval as uncertain by about 0.2 in Gamma until the
  CR3BP correction is done. (The same term may explain Broucke's "+ mu^(1/3)(4 + x0^3)" on p.130, but the x^3
  coefficient differs, so that is only a possibility.) Hg's stable interval interior is covered pointwise (C1 = 10.3 at Gamma = 1); its lower
  end (f four times) is an excluded bifurcation.
- **R16, Newton-Arenstorf planet-moon cyclers:** does not cover. Third-species orbits stay at O(mu^(1/3)) from the
  moon; the R16 orbits pass both primaries. The existence result for those is Perko 1974 (held) and Arenstorf
  1963 (held).
- **#945 R2 (DRO hubs):** the continuation of f (DRO family) to small mu is covered by Lemma 1 / Theorems 1-3
  except at bifurcation points.

## 4. Citation mining

| Cited work | Held? | Wanted list |
|---|---|---|
| [1] Arenstorf 1963, Amer. J. Math. 85:27 | held | - |
| [2] Birkhoff 1915 | held | - |
| [3] Broucke 1968, JPL TR 32-1168 | held | - |
| [4] Guillaume 1969, A&A 3:57 | not held | not listed; low priority (first/second-species bifurcation background) |
| [5] Henon 1969 | held | - |
| [6] Hill 1878 | not held | not listed |
| [8] Moulton 1906, Trans. AMS 7:537 | not held | not listed |
| [9] Perko 1974 SIAM | held | - |
| [10] Perko 1976 Celest. Mech. 14:395 | held | - |
| [11] Perko 1981 SIAM J. Appl. Math. 41:181 | held (`perko-1981b-...`) | - |
| [12] Perko 1981 Celest. Mech. 24:155 | held | - |
| [13], [14] Perko 1982 I, II | supplied this batch | row 21 |
| [15] Perron 1937, Math. Ann. 113 | not held | not listed |
| [17] Siegel 1950-51, Math. Nachr. 4; [18] Siegel & Moser 1971 | not held | not listed |
| [19] Wintner 1926 ("Winter" in the print), Math. Z. 24; [20] Wintner 1947 | not held | not listed |

No new row proposed. Guillaume 1969 is the only paper that bears on the bifurcation behaviour; optional.

*Filed as `cyclers_pdf/papers/perko-1983-periodic-solutions-restricted-problem-analytic-continuations-hill-problem-small-mu-celest-mech-30-115-doi-10.1007-BF01234301.pdf`. The Perko check scripts and the combined theorem-scope note are filed beside the 1982 I PDF (`perko-1982a-...-<file name>`).*

*Wanted-list row numbers in this digest are the pre-batch-36 numbering; the list was renumbered in batch 36.*
