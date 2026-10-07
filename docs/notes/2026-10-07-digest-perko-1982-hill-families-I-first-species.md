# Digest: Perko 1982 I, "Families of Symmetric Periodic Solutions of Hill's Problem I: First Species Periodic Solutions for C << -1" (#960 batch 36)

Lawrence M. Perko (Northern Arizona University), "Families of Symmetric Periodic Solutions of Hill's Problem I:
First Species Periodic Solutions for C << -1", American Journal of Mathematics 104(2):321-352 (April 1982),
doi 10.2307/2374162 (Crossref: title, volume 104, issue 2, first page 321, 1982-04 all agree). Received
15 October 1980. NSF grant MCS 7703591A01. The title on the page image (p.321) matches.
- Supplied file `2784c018-perko1982.pdf`, 33 pp. (1 JSTOR cover page + journal pp.321-352),
  md5 3a1876f1436feaae63ccb9582c1c0ca6. JSTOR scan with a text layer (equations garbled).
- **Proposed corpus filename:**
  `perko-1982a-families-symmetric-periodic-solutions-hill-problem-I-first-species-C-much-less-minus-1-amer-j-math-104-321-doi-10.2307-2374162.pdf`
- **How I read it:** the whole text from the text layer; on 150 dpi page images: the abstract and sec. 1
  (pp.321-323), Theorem 3 with b(t) (p.343, PDF p.24), the Existence Theorem and its proof end (pp.344-346, PDF
  pp.25-27), the reference list (p.352).
  Checks: `check_perko.py` -> `check_perko.out` (sec. A, B) and `check_perko_I_bt.py` -> `check_perko_I_bt.out`,
  kept in this folder.
- Wanted list: **row 21** (all three Perko 1982-83 papers). This paper is one third of it.
- Combined theorem-scope note for all three papers: `theorem-scope.md` (this folder).

## 0. Verdict

**The existence proof of one first-species family of Hill's problem for C << -1, which is Stromgren's family f
(the retrograde satellite / DRO family), with a first-order asymptotic form. It has no tables.**
- Theorem scope (p.345-346, image): "There is a C1 < 0 such that for all C < C1 there is a unique one-parameter
  family of symmetric periodic solutions ... symmetric with respect to both the x and y axes." C1 is not
  given. Symmetric only. The proof is an existence-plus-error-bound argument (Lemmas 1, 2), plus a Poincare
  continuation proof in the Appendix.
- f is the only first-species family for C << -1 "to first order" (p.348); the text says "it appears that the
  only 1st species family ... is the family f" (p.349). That is an argument, not a theorem.
- **The printed first-order coefficients are wrong.** The theorem prints
  `x(0, C) = |C|^(1/2) + (1 + 4 B1)/|C| + O(|C|^(-5/2))`, `y'(0, C) = -2|C|^(1/2) - 2(1 + 3 B1)/|C| + ...`,
  `B1 = int_0^(pi/2) sin t/(1 + 3 sin^2 t)^(3/2) dt`, and `T = 2 pi + (28/|C|^(3/2)) int_0^(pi/2) dt/(1 + 3 sin^2 t)^(1/2) + O(|C|^(-3))`.
  - B1 = 1/4 exactly (u = cos t), so the theorem says x(0, C) = |C|^(1/2) + 2/|C| and T = 2 pi + 30.19/|C|^(3/2).
  - My integration of family f (which reproduces Henon 1969 Table 3 to all printed digits at Gamma = -6 and
    -10) gives `x(0, C) - |C|^(1/2) = -0.0324 |C|^(-5/2)` (no 1/|C| term at all) and
    `T/2 - pi = -2.1562 |C|^(-3/2)` at C = -300 (sec. 3). The period correction is negative, not +30.19.
    It matches -2K, K = int_0^(pi/2) dt/(1+3 sin^2 t)^(1/2) = 1.078258, to 4 digits.
  - **The cause is a sign error in b(t).** Theorem 3 (p.343, image) defines b(t) = 2 i . Phi_rv(t1, t) f(r0(t)) +
    j . Phi_vv(t1, t) f(r0(t)) and writes it as "= -2 sin t/r0^3(t)", r0^3 = (1 + 3 sin^2 t)^(3/2). Evaluating the
    printed combination with Perko's own Phi (eq. 9; it reproduces a11 = 4, a14 = 2, a41 = -6, a44 = -3 of p.344)
    and f(r) = -r/r^3 gives +2 sin t/r0^3 at every t tried (`check_perko_I_bt.out`). With the corrected sign,
    int b = +1/2, so xi* = -2 int b + 1 = 0 and eta0 = 3 int b - 2 = -1/2: x(0, C) = |C|^(1/2) + O(|C|^(-5/2)) and
    y'(0, C) = -2|C|^(1/2) - 1/(2|C|) + ..., which is what the integration shows. The period formula uses the same
    integral; I did not re-derive it, but the integration gives -4K, not +28K.
  - Henon 1969 eq. (11) (xi0 = -(-Gamma)^(1/2), Table 6) is therefore better than Perko's "improvement" (p.346).
  - The existence statement is not affected. Do not use the printed coefficients as seeds or goldens.
- **What it gives the project:** a cited theorem that the DRO-type family f exists for all C below some C1 and
  tends to the ellipse |C|^(1/2)(cos t, -2 sin t) uniformly (p.346). For R15 and R16: background only (sec. 4).
- **Catalogue implication (PROPOSAL only):** none.

## 1. Content

- **Sec. 1 (pp.321-323).** History: Hill 1878 (families f, g), Wintner 1926 (convergence for C >> 1), Hopf 1929,
  Siegel-Moser (families a, c near the equilibria), Conley 1963. Perko defines first species (limit orbit bounded
  away from the origin) and second species (limit orbit with an "angular point" at the origin) for C << -1. He
  says Henon 1969 [5] showed four families (a, c, g, g') end in second-species limits and that "the result of
  continuing all but these four families is at this time unknown" (p.323).
- **Sec. 2 (pp.323-333).** Hill's equations (1); Jacobi constant C = -(x'^2 + y'^2) + 3x^2 + 2/r (same as Henon's
  Gamma). Symmetry lemma and Birkhoff's theorem (two perpendicular x-axis crossings give an x-symmetric periodic
  orbit); Hill's theorem (a y-axis crossing at t = 0 and an x-axis crossing at T/4 give an orbit symmetric about
  both axes, p.326; it does not hold in the restricted problem). Scaling x -> x/|C|^(1/2), eps = |C|^(-3/2) (eqs 2-3);
  at eps = 0 the general solution is the epicycle x = k1 cos t + k2 sin t + 2k3, y = -2k1 sin t + 2k2 cos t - 3k3 t + k4.
  Lemma 1 (majorising function for error bounds) and Lemma 2 (implicit-function lemma for the second crossing).
- **Sec. 3 (pp.333-349).** Generating solution (7) x = cos t, y = -2 sin t; Theorems 1-3 (error bound, y-axis
  crossing time, perpendicularity), the Existence Theorem (above). Fundamental matrix (9) printed in closed form;
  a11 = 4, a14 = 2, a41 = -6, a44 = -3 at t = pi/2 (p.344, image).
  - The other first-species generating solutions x = cos t, y = -2 sin t + k4 (k4 != 0, +-2) do not generate
    symmetric periodic orbits (pp.348-349). k4 = +-2 pass through the origin: they are the type B generating
    orbits of paper II (family g).
- **Appendix (pp.349-351).** The same existence by Poincare's continuation method; it does not work for second
  species (singular generating orbits), which is why sec. 2's method is needed.

## 2. Frame and conventions

Origin at the small body (Hill's problem, mass parameter scaled out). x towards the far primary's opposite side
as in Henon (x'' = 2y' + 3x - x/r^3). C = Henon's Gamma. Perko's family f has x(0) > 0 with y'(0) < 0
(clockwise, retrograde); Henon tabulates the mirror point xi0 < 0 with eta' > 0. Time is not scaled.

## 3. Checks (`check_perko.py`, `check_perko.out`)

DOP853, rtol 1e-12, atol 1e-14; shoot x(0) so that the first y-axis crossing is perpendicular.

| Gamma = C | xi0 (mine) | T/2 (mine) | Henon 1969 Table 3 | Perko x0 = sqrt(abs C) + 2/abs C | Perko T/2 |
|---|---|---|---|---|---|
| -6 | -2.4492652 | 3.001947 | -2.449265, 3.00195 | 2.782823 | 4.168719 |
| -10 | -3.1621945 | 3.074950 | -3.162194, 3.07495 | 3.362278 | 3.618958 |
| -100 | -9.9999997 | 3.139438 | - | 10.020000 | 3.156688 |
| -300 | -17.3205081 | 3.141178 | - | 17.327175 | 3.144498 |

- Positive control: my shooter matches Henon's printed f values (same criterion: perpendicular crossings). So
  the disagreement is with Perko's coefficients, not with the integration.
- (abs(xi0) - sqrt(abs C)) * C^2 = -0.0081, -0.0083, ..., -0.0032 (C = -100), -0.0019 (C = -300): the
  correction is about -0.0324 |C|^(-5/2). (T/2 - pi) |C|^(3/2) = -2.052, -2.107, ..., -2.1562 (C = -300).
- Read on the image (pp.343-346): b(t) and its closed form in Theorem 3 (p.343), B1 = int b in the proof (p.344),
  the theorem formulas, B1's definition with sin t in the theorem (p.346), "28" in the period formula. I
  re-read p.345-346; the print is clear.

## 4. What it gives R15, R16, #974, #945

- R15 (Henon 2003 asymmetric families; Hg at the moons): background only. f is not one of R15's objects; the
  paper is symmetric-only and does not mention asymmetric orbits.
- R16 (Newton-Arenstorf planet-moon cyclers): does not cover (Hill-limit orbits around the moon only).
- #945 R2 (DRO/f hubs): the existence of f for C << -1 is a theorem here; the numbers to use are Henon 1969's.

## 5. Citation mining

Reference list read on the p.352 image (14 items). Held checks: `ls cyclers_pdf/papers | grep -i` and
CORPUS_INDEX grep for each; wanted list grep.

| Cited work | Held? | Wanted list |
|---|---|---|
| [1] Birkhoff 1915, Rend. Circ. Mat. Palermo 39:1-334 | held (`birkhoff-1915-...`) | - |
| [2] Conley 1963, CPAM 16:449-467 | not held | not listed |
| [3] Gravalos 1941, Harvard PhD (algebraic integrals of Hill's equations) | not held | not listed |
| [4] Graves 1956 (real variables textbook) | not held | not needed |
| [5] Henon 1969, A&A 1:223-238 | held (+ tables.txt) | - |
| [6] Hill 1878, Amer. J. Math. 1:5-26, 129-147, 245-260 | not held | not listed |
| [7] Hopf 1929, Sitzber. Preuss. Akad. Wiss. phys.-math. 1:401-413 | not held | not listed |
| [8] Lang 1962 (differentiable manifolds) | not held | not needed |
| [9] Perko 1982 II ("supply info later") | supplied in this batch | row 21 |
| [10] Poincare, Methodes nouvelles, v.1-3 | not held (no poincare-18*/methodes/nouvelles file; not in CORPUS_INDEX) | not listed |
| [11] Siegel & Moser 1971 | not held | not listed |
| [12] Szebehely 1967 | held | - |
| [13] Wintner 1926, Math. Z. 24:259-265 | not held | not listed |
| [14] Wintner 1928, Math. Z. 28:430 (Hill's orbit of maximum lunation) | not held | not listed |

No new wanted-list row proposed. Conley 1963 is the only candidate (small periodic orbits around the smaller
primary for low energy); low priority.

*Filed as `cyclers_pdf/papers/perko-1982a-families-symmetric-periodic-solutions-hill-problem-I-first-species-amer-j-math-104-321-doi-10.2307-2374162.pdf`. Check scripts, outputs and other files named above are filed beside it (for all three Perko papers, beside the 1982 I PDF) as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the pre-batch-36 numbering; the list was renumbered in batch 36.*
