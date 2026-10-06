# Digest: Hénon 1969, "Numerical Exploration of the Restricted Problem. V. Hill's Case: Periodic Orbits and Their Stability" (#960 batch 21)

M. Hénon (Observatoire de Nice), Astronomy & Astrophysics 1:223-238 (1969), ADS 1969A&A.....1..223H.
Received 15 November 1968.
- Filed as `cyclers_pdf/papers/henon-1969-numerical-exploration-restricted-problem-V-hill-case-periodic-orbits-stability-aa-1-223-ads-1969AA-1-223H.pdf`.
  - This is an `ocrmypdf --force-ocr -l eng` copy of the image-only ADS scan. The source scan is md5
    6fe6bb47ee1b6bd5c45ed8c67eefcc0b, 16 pp.
- **Companion `...-tables.txt`: Tables 1-12 transcribed in full from 400 dpi page images.**
  - Three witnesses: tesseract with a digit whitelist, macOS Vision, and a reading of the image, which
    decides.
  - 1062 printed cells; 282 had at least one witness differing; 7 were real numeric misreads, all
    settled on the image.
  - Internal checks:
    - monotonic trends;
    - Table 2's first row equals 3^(4/3) and 3^(-1/3);
    - the Gamma -> -inf T/2 limits (pi, 2 pi, the roots of tan t = (3/4)t);
    - Table 7 recomputed from scratch (50 cells; one last-digit rounding edge);
    - Tables 6, 8, 9, 11 and 12 recomputed from the asymptotic formulas.
  - Two rows of Table 11 have the root rho_0 solved about 1e-4 inaccurately in the original. The rows are
    internally consistent and are kept as printed.
- I spot-checked Table 2 (Gamma = 4.32675, 0, -100) and Table 3 (Gamma = 4, -3.5) against my own
  reading of the page image. They agree.
- Wanted-list rank 19 before batch 21; removed.

## 0. Verdict

**H5, the Hill-problem (mu -> 0) baseline for the held Hénon 2003 paper (H-series families a, c, f, g,
g').**
- **Family f** consists of the retrograde satellites of the second body. These are the Hill-problem
  distant retrograde orbits (DROs). **All its orbits are stable.** Table 3 gives Gamma, xi, T/2 and the
  stability index a, from Gamma = +inf down to -inf. For example, Gamma = 4: xi = -0.20421,
  T/2 = 0.26837, a = 0.8518. Gamma = -3.5: xi = -1.87053, T/2 = 2.84686, a = 0.2774.
- **Families a and c** start as libration orbits around the Lagrange points. They are all unstable.
  Table 2 starts at Gamma = 3^(4/3) = 4.32675.
- **Families g and g'** are the direct satellites. They are stable in some intervals.
  - Critical orbits (Table 10): g'1 = g1 at Gamma = 4.49999, xi0 = 0.28350, types 1 and 3; g'2 at
    4.27143 / 0.58769, type 6; g'_2 3 at -4.69219 / 0.004302, type 5; g'_2 4 at -4.70479 / 0.004283,
    type 4.
  - Orbits with Gamma > 4.499986 are stable (Fig. 4 caption).
- Asymptotic forms for Gamma -> +-inf are derived and checked (Tables 6-9, 11-12).
- **Use:**
  - `#945` R2 / `#953` R7: the single-moon building blocks in the Hill limit. The DRO stability (f) and
    the direct-satellite stability limit (g) are sourced control values for any Hill-scaled Jovian or
    Saturnian satellite orbit.
  - Pairs with Hénon 2003 (held): there, g3 is the Anderson 2018 petal family.

## 1. Content (READ from the OCR, with tables from the companion)

- Units: Hill's scaling of the restricted problem as mu -> 0 (Hill 1886; Szebehely 1967 p.609), with
  Jacobi constant Gamma.
- Only Matukuma (1930-1957) had searched systematically before, with desk computers. Table 1 compares
  the earlier hand-computed orbits with Hénon's exact values: the differences reach several units in the
  second decimal.
- Families: a and c are mirror images; f is retrograde and symmetric about the origin; g is direct; g'
  branches from g at g1 and is the only family with a double-periodic branch.

## 2. Citation mining

- Hill 1886, Matukuma 1930-1957, Kevorkian 1962, Stumpff 1965, Szebehely 1967 (held), and Hénon
  1965a/b (1965b now held; 1965a on the wanted list) and 1966a/b.
- Matukuma and Hénon 1966a/b (Bull. Astron. 1:57, 1:49) are not held. They are added to one
  low-priority wanted row.
