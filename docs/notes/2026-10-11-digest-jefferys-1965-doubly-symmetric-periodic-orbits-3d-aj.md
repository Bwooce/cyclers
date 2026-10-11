# Digest: Jefferys 1965, "Doubly Symmetric Periodic Orbits in the Three-Dimensional Restricted Problem" (#960 batch 43)

W. H. Jefferys (Yale University Observatory and Van Vleck Observatory), The Astronomical Journal 70(6):393-394,
No. 1331, August 1965. Received 15 March 1965. DOI 10.1086/109752 (Crossref checked 2026-10-11: title, AJ 70, page
393, author Jefferys, issued 1965-08). ADS 1965AJ.....70..393J.
- Proposed filename: `cyclers_pdf/papers/jefferys-1965-doubly-symmetric-periodic-orbits-three-dimensional-restricted-problem-aj-70-393-doi-10.1086-109752-ads-1965AJ-70-393J.pdf`.
  - File the `ocrmypdf --force-ocr -l eng` copy, md5 6a34954cf4e6c6cdf51f82835221f900 (re-made 2026-10-11 with ocrmypdf --force-ocr -O 3; render-checked against the ADS original) (0.64 MB); the source is 0.13 MB,
    because --force-ocr rasterizes the 200 dpi colour JPX pages; `-O 3` at filing would shrink it).
  - Source: ADS scan, 2 pp., md5 eb4d9790cd6a25296b19de4b5121d147, fetched 2026-10-11 from articles.adsabs.harvard.edu.
  - Render check, dark-pixel fraction at 100 dpi, original against OCR copy: p.1 0.042/0.045, p.2 0.015/0.016
    (`hitzl/render_check.out`). Legible, no blank pages.
- **How it was read:** both pages on page images (100 dpi; equations checked by eye). The text layer was not used.

## 0. Verdict

**A two-page existence proof. It has no tables and no numbers.** The brief expected "tables of ICs"; there are none.
So there is no yaml and no integration for this paper.
- It proves, by analytic continuation from mu = 0 in Delaunay-Poincaré variables, that inclined circular orbits of
  the two-body problem continue to doubly symmetric periodic orbits of the 3D circular restricted problem for small mu.
- "Doubly symmetric" means symmetric about both the x-y and the x-z planes: at t = 0, xdot = y = z = 0, and at
  t = T/4, xdot = y = zdot = 0 (p.393).
- Only the first kind (e = 0) is treated. "The question of the existence of periodic orbits of the third kind (e ≠ 0)
  is therefore left open" (p.394).
- **What it gives the project:** an attribution for the existence of 3D doubly symmetric orbits near inclined circular
  orbits (`#994`-type 3D lifts, the 1960s 3D lineage). No positive control.
- **Proposal:** keep wanted-list row 58's 1965 half as received and low priority. The numeric lead is its only
  reference, Goudas 1963 (below), which computed such orbits at large mu.

## 1. Content (page images)

- Frame (p.393): Jupiter (mass mu) and the Sun (1 - mu) fixed on the x axis of axes rotating with unit rate about z.
  Eq. (1): xdd - 2 ydot = dOmega/dx, ydd + 2 xdot = dOmega/dy, zdd = dOmega/dz. Omega depends on y and z only through
  y², z².
- Symmetry (p.393, read on a 300 dpi crop): if x(t), y(t), z(t) solves (1), so does x(-t), -y(-t), ±z(-t). "Therefore, if at t = 0,
  ẋ = y = z = 0, and if at t = ¼T, ẋ = y = ż = 0, then the motion is periodic with period T, and is symmetric with
  respect to both the x-y and x-z planes."
- Eqs. (2)-(4): canonical equations with Hamiltonian F = ½(X1 + X3)^{-2} + X1 - ½(X2² + Y2²) + mu F1, in
  X1 = L - G + H, Y1 = lambda - lambda', X2 = [2(L - G)]^{1/2} cos(lambda' - varpi),
  Y2 = [2(L - G)]^{1/2} sin(lambda' - varpi), X3 = G - H, Y3 = lambda - Omega.
- Periodicity conditions: at t = 0 Y1 = i pi, Y2 = 0, Y3 = j pi; at t = ¼T, Y1 = (i + k) pi, Y2 = 0,
  Y3 = (j + m + ½) pi.
- Eq. (5), mu = 0: Y1 = [(X1 + X3)^{-3} - 1]t + i pi, Y2 = X2(0) sin t, Y3 = (X1 + X3)^{-3} t + j pi.
- Eqs. (6)-(7): X2(0) = 0, ¼(n - 1)T0 = k pi, ¼ n T0 = (m + ½) pi, with n = (X1 + X3)^{-3} = L^{-3}. So
  T0 = 2 pi (2m - 2k + 1) and n/(n - 1) = (2m + 1)/(2k), with n ≠ 1.
- p.394: the Jacobian Delta = ∂(Y1, Y2, Y3)/∂[t, X1(0), X2(0)] at t = ¼T0 equals -3(X1 + X3)^{-4} sin ¼T0 = ±n ≠ 0,
  so the orbit continues for small mu (implicit function theorem). "The existence proof depends on both symmetries."
  If one symmetry is dropped, Delta must be evaluated at ½T0 = l pi, and vanishes.
- p.394: "Goudas (1963) computed doubly symmetric orbits of various types, some of which appear to coincide with the
  present orbits. Since he used large values of mu, it is not unreasonable that the orbits continue to exist over a
  wide range of mu."
- Acknowledgement: V. G. Szebehely suggested the Delaunay-Poincaré elements and the work of Goudas.

## 2. Checks (ours)

- Consistency of eq. (7): from ¼(n - 1)T0 = k pi and ¼ n T0 = (m + ½) pi, dividing gives n/(n - 1) = (2m + 1)/(2k),
  and subtracting gives ¼T0 = (m - k + ½) pi, so T0 = 2 pi (2m - 2k + 1). Both match the print.
- No numeric content to integrate.

## 3. Citation mining

| Reference | Held? |
|---|---|
| Goudas, C. L. 1963, Icarus 2:1 | NOT held (`ls cyclers_pdf/papers | grep -i goudas`: nothing; not in CORPUS_INDEX as a file; not on the wanted list). Crossref: "Three-dimensional periodic orbits and their stability", Icarus 2:1-18, doi 10.1016/0019-1035(63)90003-4. |

- **New lead, proposed as one wanted-list row:** Goudas 1963, Icarus 2:1-18, doi 10.1016/0019-1035(63)90003-4
  (Elsevier, paywalled). It is the numeric source behind this paper: computed 3D doubly symmetric orbits and their
  stability at large mu. It pairs with Bray & Goudas 1967 (row 52) and the Hénon 1973 vertical-stability work (held).
  Priority: low to medium, for 3D symmetric-orbit controls.
- Wanted-list row 58: this is its "Jefferys (1965)" half; the 1971 Atlas half stays not found.
- A Crossref neighbour, not cited here: an AJ 68:534 (1963) abstract, "Existence of Symmetric Periodic Orbits in the
  Restricted Problem of Three Bodies in Three-Dimensional Coordinate System", doi 10.1086/109116. Author not checked.
  Not proposed.

*Filed as `cyclers_pdf/papers/jefferys-1965-doubly-symmetric-periodic-orbits-three-dimensional-restricted-problem-aj-70-393-doi-10.1086-109752-ads-1965AJ-70-393J.pdf`.*
