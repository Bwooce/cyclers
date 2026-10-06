# Digest: Hénon 1970, "Numerical Exploration of the Restricted Problem. VI. Hill's Case: Non-Periodic Orbits" (#960 batch 30)

M. Hénon (Observatoire de Nice), Astronomy & Astrophysics 9:24-36 (1970), ADS 1970A&A.....9...24H.
Received 21 May 1970. 13 pages (pdf page n is journal page 23+n).
- Filed as `cyclers_pdf/papers/henon-1970-numerical-exploration-restricted-problem-VI-hill-case-non-periodic-orbits-aa-9-24-ads-1970AA-9-24H.pdf`.
- OCR copy (filed): `1970AA924H-ocr.pdf`, md5 fe5f678c819d3448439f9c2ae38a5ecd (force-OCR copy), 13 pp.
  Original ADS scan: supplied file `1970AA924H.pdf`, md5 45337173216e5f9e14e667e8639cf7c5
  (13 pp., one 4800x6600 JBIG2 stencil per page, no text layer except the ADS stamp).
- The first OCR copy (`ocrmypdf --redo-ocr`) lost the page images of this ADS JBIG2 scan. The filed copy was re-made with `ocrmypdf --force-ocr` (md5 fe5f678c819d3448439f9c2ae38a5ecd); every page was checked to render.
- How I read it: text from the OCR text layer (two-column layout, interleaved, garbled symbols). Every
  number below that I call "image-checked" I read on a 200 dpi render of the ORIGINAL scan: Table 1
  (p.27), Table 2 (p.32), Table 3 (p.35), the page-28, 30, 31 figures and captions. Scripts and outputs
  are filed beside the PDF.
- Wanted-list row 32 (current numbering of the 2026-10-05 list); removed in batch 30.

## 0. Verdict

**H6: the non-periodic companion to the held H5 (Hénon 1969). It maps where Hill-problem orbits are
bounded (quasi-periodic) and where they escape, with surfaces of section.** For the project it is the
sourced Hill-limit picture behind distant retrograde orbits (DROs), their stability region and the
direct-satellite limit.
- **Retrograde satellites have no upper size bound in Hill's problem.** For every Gamma down to -100
  (checked by Hénon), there is a region of bounded quasi-periodic motion round the retrograde periodic
  orbit f (p.29, p.33, Fig. 13). Exceptions: two resonant Gamma values, 0.017... and -1.41... (the
  stability index of f is -1/2, so the invariant curves vanish there).
- **Direct satellites do have a bound.** Quasi-periodic direct orbits end near Gamma = 4.4. Below
  4.326749 ... the three accessible regions merge and an ergodic orbit escapes. No direct satellite
  exists beyond about 0.29 Hill units (630,000 km for the Earth, p.33).
- **Escape criterion for stars:** Gamma_e(xi -> 0) tends to the Lagrange value 4.326749 ...; but circular
  orbits survive far below it (p.35-36).
- **Catalogue implication (PROPOSAL only):** none for a catalogue row. Use as a cited control for
  `#945` R2 / `#953` R7 (single-moon DRO and direct-satellite limits) and for the Hill-limit
  stability-region claims. I checked no catalogue row against it.
- **Corrects a common belief.** Chebotarev (1968) had concluded that retrograde orbits become unstable near
  the Hill radius 0.693361. Hénon (p.34) shows that conclusion came from the family e = 0 only. Stable
  retrograde orbits continue along family f.

## 1. Setting and method (pp.24-25; image-checked on Fig. 1 and the text)

- Hill's equations as in V: xi'' = 2 eta' + 3 xi - xi/rho^3, eta'' = -2 xi' - eta/rho^3. Jacobi constant
  Gamma = 3 xi^2 + 2/rho - xi'^2 - eta'^2 (the OCR garbles eq. 3; the check script uses this form and
  conserves Gamma to 1e-9 over 200 time units).
- Surface of section: eta = 0, positive crossings (eta' > 0), plotted in the (xi, xi') plane.
- Accessible region on the section: xi'^2 <= 3 xi^2 + 2/|xi| - Gamma. For Gamma > 3^(4/3) = 4.326749 ...
  there are three separate regions (the two outer ones are escape regions); below it there is one.
  I computed 3^(4/3) = 4.3267487 and 3^(-1/3) = 0.6933613 (Hill radius, the Lagrange distance in Hill
  units). Both agree with the printed digits.
- Symmetries: (xi, eta, xi', eta') -> (xi, -eta, -xi', eta') and the point reflection through the origin.
  Hénon uses negative crossings through the second symmetry to double the plotted points.
- Initial conditions for every plotted orbit (Table 3): xi as listed, eta = 0, xi' = 0,
  eta' = (3 xi^2 + 2/|xi| - Gamma)^(1/2); n = number of positive crossings plotted each way.

## 2. Classes of orbit and the boundaries (pp.25-30, Figs. 2-14)

Hénon names three types: quasi-periodic (points on closed invariant curves), ergodic (scattered points
filling a region), escape (distance grows without limit).
- **Gamma >= 5 and down to 4.5 (Figs. 2, 3).** Invariant curves surround f (retrograde) and g (direct).
  Curve size rises with eccentricity. Nothing new happens above 5.
- **Gamma = 4.499986 ... (critical orbit g1, type 3).** g loses stability; two stable g' orbits appear.
  Fig. 4 (Gamma = 4.4) shows the curve region of g splitting in two. The unstable g is surrounded by
  hyperbolic arcs. A set of 7 islands shows. This matches the Gamma > 4.499986 stability limit in the
  held H5 digest.
- **Gamma = 4.35 (Fig. 5).** First ergodic orbit: all scattered points belong to one orbit; two isolated
  curve regions remain round the two g' points.
- **Gamma < 4.326749 ...** The central region opens. Orbits in the ergodic region escape in general.
  "Down to 4.326749 the curve region fills practically the whole accessible space" (p.27).
- **Gamma < 4.27143 ... (critical orbit g'2).** g' become unstable and their curve regions vanish.
- **Retrograde side, family f and the triple-periodic family g3 (new name, p.27).**
  - g3 is a point-symmetric triple-periodic family, computed earlier by Matukuma (1957, his families X and
    Y). By analogy with the mu = 1/2 case (I, Fig. 3; III, Fig. 19) Hénon calls it g3. Hénon did not
    follow its ends or its stability.
  - Table 1 (image-checked) gives xi and T/2 for g3 at 25 values of Gamma, from 3.8 down to -9. Examples:
    Gamma = 3: xi = -0.76939, T/2 = 2.33617. Gamma = 0.5: xi = -0.66969, T/2 = 2.98423. Gamma = -2:
    xi = -1.58227, T/2 = 7.59033. Gamma = -9: xi = -3.90443, T/2 = 8.64646.
  - **My check:** I integrated Hill's equations from (xi, 0, 0, +eta') at T/2 for Gamma = 3, 2, -2, -9.
    The state at T/2 equals the point reflection of the start. The largest residual is 2.4e-3 in eta at
    Gamma = -9 (where the speed is 7.4); it is about 2e-4 at the other three values. This fits the
    5-digit rounding of xi.
    Gamma = 3: xi(T/2) = 0.76946 against 0.76939; eta'(T/2) = -1.17279 against -1.17275. The
    OCR text misread four Table 1 digits (3.19828 for 3.19328; -1.33953 for -1.32953; 7.90873 for 7.90573;
    a row -6 that is -5 on the page). The values above are from the image.
  - g3 limits the size of the curve region round f over a long interval, Gamma from about 1 to -2.5
    (p.28-29). Fig. 6 (Gamma = 4, 3, 2, 1, 0.5) shows the curve region shrink and turn triangular.
  - **Gamma = 0.017 ...:** f meets g3; stability index a = -1/2; rotation angle 2 pi/3 (the 1:3
    resonance); the curves round f disappear (Fig. 7, drawn for Gamma = 0). The curves reappear with the
    triangular shape inverted (Fig. 8, Gamma = -0.5 and -1).
  - **Gamma = -1.41 ...:** a second intersection of f and g3 (a = -1/2 again). For smaller Gamma the curves
    reappear (Fig. 9).
  - **The a = 0 (1:4) resonance** of f gives no singular case. It meets two families of quadruple-periodic
    orbits, which removes the singularity (also seen for mu = 1/2, paper IV).
  - **Gamma < -1.41:** no more intersections with g3 or any triple family, and f is always stable (V).
    Hénon verified a curve region at Gamma = -100.
- **Fig. 11 (Gamma = -16, xi = -4.35).** Hénon describes an elliptic retrograde motion, axis ratio 2:1,
  period of order 2 pi, whose centre librates slowly along a very elongated vertical ellipse.
  - **My check:** I integrated from xi = -4.35, eta' = +8.5573 (Gamma = -16). The orbit stays bounded
    for 200 time units, about 32 periods. xi stays within +-4.350 and eta within -13.70 and +13.57, which
    matches the axes of the figure (xi +-4, eta +-14). The reverse sign of eta' escapes at once. This is
    a short integration, not a proof of boundedness.
- **Boundaries (Fig. 12 and Fig. 13).** The approximate boundaries on the xi axis are plotted against
  Gamma, with accuracy about 0.01 in xi. I do not quote numbers from the figures. Fig. 13 shows the curve
  region round f continuing for Gamma -> -infinity.
- **Fig. 14.** For a start on the xi axis it gives the set of initial velocities that stay bounded. From
  xi = 0 to about 0.65 the curves are nested round the origin. Near xi = 0.65 is the first 1:3 resonance
  along f; between 0.65 and 0.68 there are two little curves; above 0.68 one small curve, which vanishes at
  xi about 1.2 (second 1:3 resonance) and then reappears and grows without limit. For large xi the
  eta' window is narrow and xi' is less critical. The boundaries are interpolated and approximate.

## 3. Applications (pp.31-36)

- **Natural satellites (Table 2, image-checked).** Unit of length mu^(1/3) a' (a' = planet orbit radius),
  unit of time 1/n'. For each satellite Hénon computes the time-averaged Gamma (eq. 9) and the section
  points xi1 = eps a (1 - e), xi2 = eps a (1 + e) (eq. 13), using Allen (1963) data. eps = +-1 is the
  direction of rotation.
  - Examples: Moon (Earth 1) Gamma = 6.5403, xi1 = 0.16769, xi2 = 0.18717. Jupiter 8 (Pasiphae)
    Gamma = 2.4760, xi1 = -0.18416, xi2 = -0.42972. Jupiter 11: 2.5183, -0.23354, -0.35547. Neptune 2
    (Nereid): 30.6425, 0.00789, 0.05789.
  - **My check:** I recomputed Gamma from eq. 9 using a = |xi1 + xi2|/2 and e = (xi2 - xi1)/(xi2 + xi1).
    Moon 6.5405 (6.5403), Jupiter 6 7.4881 (7.4881), Jupiter 8 2.4761 (2.4760), Saturn 9 6.5653 (6.5653).
    Neptune 2 differs by 0.0016 and Jupiter 1 by 0.15, because xi1, xi2 are printed to three significant
    figures and the reconstruction of a and e loses accuracy. These are rounding effects, not misprints.
  - A gap in Gamma: no satellite has 8 < Gamma < 26. The 22 of 31 satellites with Gamma > 26 are Kepler
    orbits with a small perturbation. The outer Jovian satellites (Jupiter 8, 9, 11, 12, Gamma 2.37 to 2.73,
    all with negative xi, that is retrograde) lie below the Lagrange value 4.326749 yet stay bound, because the retrograde region
    reaches far lower Gamma. All 31 lie inside the curve region (Fig. 15), none near its edge.
  - Caveat from Hénon: he ignores inclination (about 30 deg for most) and osculating eccentricity varies
    widely for the four outer Jovian satellites.
- **Artificial satellites.** Hill unit for the Sun-Earth system: 2.168e6 km and 58.13 days. I checked
  mu^(1/3) a' = (3.0e-6)^(1/3) x 1.496e8 km = 2.16e6 km. Direct Earth satellites are quasi-periodic to
  Gamma = 4.4, mean distance 0.29 unit, about 630,000 km (0.29 x 2.168e6 = 6.29e5). Retrograde ones can
  exist at any distance in this model. Real limits come from the Moon, the Earth's finite mass, the
  eccentricity of the Earth's orbit and the planets.
- **Mass-ratio reach.** Hénon and Guyot (1970; HELD) show retrograde orbits are stable for all
  0 < mu < 0.0477 ... (eq. 14), so the Hill result is not only a mu -> 0 limit.
- **Comparison with Chebotarev et al. and Hunter (p.33-34, Fig. 16).** He puts their computed orbits in
  the (Gamma, xi) plane: bounded ones should fall in the curve (or ergodic) region and escaping ones in
  the escape region. Agreement is "quite satisfactory" with a few marginal exceptions. Their
  initial conditions are circular osculating orbits (e = 0), which follow f and g only for high Gamma. His
  inclination cut is under 33 deg. The paper does not give counts, so I have no number for the agreement.
- **Star clusters.** The cluster tidal radius equals the Lagrange distance 0.693361; the stars in
  Fig. 14 stay within xi about 0.68. Encounters launch stars on escape orbits; Gamma_e(0) is near the
  Lagrange value (supported by Wielen 1969). Circular-type orbits exist for much lower Gamma.

## 4. What this gives the project

- A published qualitative map of bounded and escape sets of Hill's problem to test any Hill-limit
  search or census against: the 4.4 and 4.326749 direct limits, the unbounded retrograde side, and the
  two f resonance gaps at Gamma = 0.017 ... and -1.41 ... Both gap Gamma values are printed to three
  decimals only. I did not compare them with the stability table of V.
- Sourced numbers that I could check by integration: Table 1 (g3) at four Gamma values (Section 2).
- No cycler or tour content. Do not cite it for periodic orbits near both primaries.

## 5. Citation mining

Held status checked with `ls cyclers_pdf/papers | grep` and CORPUS_INDEX.
- Hénon 1965a (I) and 1965b (II): HELD (`henon-1965a-...`, `henon-1965b-...`).
- Hénon 1966a, 1966b (III and IV, Bull. Astr. Paris 1): not held. Wanted-list row 65 (also lists
  Matukuma 1930-1957).
- Hénon 1969 (V): HELD (`henon-1969-numerical-exploration-restricted-problem-V-...`).
- Hénon and Guyot 1970 (Sao Paulo symposium): HELD (`henon-guyot-1970-...`).
- Matukuma 1957 (Sendai Astr. Rap. 51): not held; row 65.
- Chebotarev and Bozhkova 1960/1962/1963, Chebotarev and Volkov 1961, Chebotarev 1964/1966/1968:
  not held, not on the wanted list. **New candidate (medium):** Chebotarev 1968, Bull. Inst. Theor. Astr.
  11:341, because it is the "retrograde orbits become unstable at the Hill radius" claim that Hénon
  refutes. Check it before anyone cites that claim in the DRO work.
- Hunter 1967 (MNRAS; printed as vol. 186, p.245, which looks like a misprint; I did not check): not
  held, not listed. Low priority.
- Wielen 1969 (Habilitationsschrift, Heidelberg), Hagihara 1952, Herget 1968, Allen 1963: not held, not
  listed. Low priority; Allen is a reference book.
- Not all the series is held: IV and III (row 65) would complete the mu = 1/2 comparison.
- Hénon's quoted mu -> 0 values (3^(4/3), 3^(-1/3), 4.499986, 4.27143) agree with the H5 digest.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
