# Digest: Newton 1959, "Periodic Orbits of a Planetoid Passing Close to Two Gravitating Masses" (#960 batch 30)

Robert R. Newton (Applied Physics Laboratory, Johns Hopkins University), Smithsonian Contributions to
Astrophysics 3(7):69-78 (1959), doi 10.5479/si.00810231.3-7.69, ADS 1959SCoA....3...69N.
Work supported by the Navy Bureau of Ordnance, contract NOrd 7386. 10 pages (pdf page n is journal page 68+n).
- Filed as `cyclers_pdf/papers/newton-1959-periodic-orbits-planetoid-passing-close-two-gravitating-masses-smithsonian-contrib-astrophys-3-69-doi-10.5479-si.00810231.3-7.69-ads-1959SCoA-3-69N.pdf`.
- OCR copy (filed): `1959SCoA369N-ocr.pdf`, md5 ffd9fa98b6ae3e43d901fbd8a7157f48 (force-OCR copy), 10 pp.
  Original ADS scan: supplied file `1959SCoA369N.pdf`, md5 f4ab8143a3dc7e7398fb9c74ae889d94
  (10 pp., 4352x6120 JBIG2 stencil per page).
- The first OCR copy (`ocrmypdf --redo-ocr`) lost the page images of this ADS JBIG2 scan. The filed copy was re-made with `ocrmypdf --force-ocr` (md5 ffd9fa98b6ae3e43d901fbd8a7157f48); every page was checked to render.
- How I read it: text from the OCR text layer; Table 1 (p.76), eq. 18 and the numerical-method text
  (p.75) read on 170-220 dpi renders of the ORIGINAL scan. Equations 6 to 17 I read from the text layer
  for structure only and do not quote them. The paper is entirely about mu = 1/82.45.
- Wanted-list row 28 in the current numbering of the 2026-10-05 list (the task said 29; row 29 is Ollé
  1989). Remove row 28 on filing.

## 0. Verdict

**The earliest published numerical family of periodic orbits of the Earth-Moon restricted problem that
pass close to both primaries on every period. All 11 orbits in Table 1 reproduce.**
- All orbits are "type 1/2": the planetoid P moves on a near-Kepler ellipse about the heavy mass E with
  sidereal period pi (half the revolution period of E and M), so P goes round E twice per revolution of
  the primaries. The synodic (rotating-frame) period is about 2 pi. P passes the light mass M once per
  period, and reaches its largest distance from E at the same time as it is nearest M (opposition of
  the lobes).
- The 11 initial states are complete (mu, x0, y0 = 0, xdot0 = 0, ydot0 from Table 1), so I integrated all
  of them. **Every orbit closes as a symmetric periodic orbit** (details in section 3).
- **Are they Earth-Moon cycler ancestors? In the weak sense, yes. In the strict sense, not shown.**
  - They are ballistic, closed in the rotating frame, and they return to the Moon every 24 to 34 days.
    The closest lunar approach is 4,660 to 81,500 km from the Moon's centre (my conversion, using
    384,400 km and a time unit of 27.3217/(2 pi) = 4.348 d).
  - Newton does not say "cycler". His aim is a general study of periodic orbits near both masses, and
    he does not compute flyby speeds, transfer cost, stability, or whether the Earth passage is
    low enough for a spacecraft. The smallest Earth distance in Table 1 is 0.035 (13,400 km from the
    centre, from my conversion), so these are not low-Earth-orbit returns. In Table 1 the Earth
    distances run from 0.035 to 0.52 (perigee) and 1.1 to 1.4 (apogee), see the discrepancy in section 3.
  - **Catalogue implication (PROPOSAL only):** none for a row. Cite as the earliest published family for
    `#948` R4; use the 11 orbits as an Earth-Moon check set for a periodic-orbit corrector (section 3).
    I did not compare any catalogue row with this table.

## 1. Zero-order construction (pp.69-71)

- Units: total mass 1, mass of M = mu, mass of E = 1 - mu, E-M distance 1, rotation period 2 pi, G = 1.
  The rotating frame xy is centred on the barycentre, M at (1 - mu, 0), E at (-mu, 0).
- For small mu, E is nearly fixed and P follows a Kepler ellipse about E. To pass close to M it needs
  apocentre near M's distance, so the major axis is about 1.26 (2 a = 1.2599 for a = (1/2)^(2/3), which
  I computed) and Kepler's third law gives T_P = 2 pi a^(3/2) = pi. In general T_P = 2 pi (alpha/beta),
  with alpha, beta coprime. Then 2 pi alpha is the full inertial period and beta is the number of
  ellipses described. In the rotating frame the orbit has beta lobes and a beta-fold symmetry axis when
  mu = 0. For mu > 0 the lobe near M is the most perturbed and the apses precess. The paper treats only
  type 1/2 (the OCR garbles the fraction; the formula T_P = 2 pi (alpha/beta) and the text "T_P = pi" on p.75 settle it as 1/2).
- Two classes at mu = 0, both for each x0 up to 1.26: **direct about E** and **retrograde about E**,
  by the sense of the orbit in the inertial frame XY. P is always retrograde about M in the rotating
  frame at closest approach (otherwise P could not come closer to E than M does).
- Standard initial state: t = 0 at closest approach to M, which is a relative maximum of the distance
  from E; y(0) = 0, xdot(0) = 0, ydot(0) = ydot0(x0). The orbits are symmetric about the x axis
  (x(-t) = x(t), y(-t) = -y(t)), so only half a period needs study; the half-period event is the third
  y = 0 crossing after t = 0.
- Initial states are chosen with x0 the single free parameter; x0 = 1 - mu + rho_M0.

## 2. Classification (pp.72-73, Fig. 3)

For mu = 1/82.45 the paper draws regions in the (rho_M0, ydot0) plane (rho_M0 is the closest-approach
distance). The curves use a crude "instantaneous encounter" estimate of M's kick:
- I: collision with E. II: escape from M. III, IV: escape from the E-M system under energy conservation
  (IV has no use in the range shown). V: zero initial velocity in the inertial frame.
- Retrograde about E, between I and III: one class, retrograde near M too, no subdivision.
- Direct about E, between I and II: three subclasses by the sense near M:
  - **direct** (direct when nearest M);
  - **zero initial velocity** (the orbit has a cusp at the closest approach to M, Fig. 5b);
  - **mixed** (retrograde when nearest M; the cusp opens into a figure-eight loop, Fig. 5a).
- A refined perturbation formula (via a canonical transformation to Kepler elements) gives the period:
  **T = 2 pi + 2 mu P2 / (eps ydot0 rho_M0)** (eq. 18, read on the image; P2 = sqrt(a (1 - e^2)) is the
  Kepler angular momentum after escape from M, with sign). It gives T < 2 pi for orbits direct about E and
  T > 2 pi for retrograde ones, and a precession of the line of apses opposite to the sense of P near E.
  Table 1 agrees with that sign rule in all 11 rows (retrograde 6.302 to 7.8925, all above 2 pi = 6.2832;
  direct, mixed, zero-velocity 5.4292 to 6.192, all below). The paper says the estimate is accurate
  only for rho_M0 not too small (footnote 6, p.75).

## 3. Numerical orbits (pp.75-77, Table 1; every number image-checked)

- Method: numerical integration in xy on a "large-scale digital computer" (type not named). Step =
  the smallest of 0.01, 0.025 rho_E^(3/2), 0.25 rho_M^(3/2). The integral of motion Jacobi constant
  drifted by almost 1 percent in the worst case (x0 = 1.2, rho_E down to 0.035) and typically one part in
  10^6. ydot0 was found by Newton-Raphson on xdot at the third y = 0 crossing, to 1e-4. The orbit of
  zero initial velocity was found by interpolation over several x0.
- Table 1 (mu = 1/82.45); ydot0 = minus the printed column; x0 = 1 - mu + rho_M0, so x0 =
  1.00, 1.05, 1.10, 1.15, 1.20 (the "five values of x0" of the text; x0 = 1.0 is rho_M0 = 0.012129 = mu).

| rho_M0 | -ydot0 | perigee | apogee | period | class |
|---|---|---|---|---|---|
| 0.012129 | 2.3314 | 0.5207 | 1.3911 | 7.8925 | retrograde |
| 0.062129 | 1.7739 | 0.2624 | 1.1022 | 6.5227 | retrograde |
| 0.112129 | 1.6604 | 0.1756 | 1.1531 | 6.3841 | retrograde |
| 0.162129 | 1.5842 | 0.1547 | 1.1968 | 6.3335 | retrograde |
| 0.212129 | 1.498 | 0.052 | 1.242 | 6.302 | retrograde |
| 0.012129 | 1.5364 | 0.1815 | 1.1966 | 5.4292 | mixed |
| 0.03211 | 1.0200 | 0.1945 | 1.1408 | 5.5718 | zero initial velocity |
| 0.062129 | 0.8475 | 0.1732 | 1.1458 | 5.7574 | direct |
| 0.112129 | 0.8303 | 0.1289 | 1.1780 | 5.9750 | direct |
| 0.162129 | 0.9124 | 0.0818 | 1.2206 | 6.1081 | direct |
| 0.212129 | 1.049 | 0.035 | 1.263 | 6.192 | direct |

  Perigee and apogee are called distances from E in the text. The period is the rotating-frame period.
  The paper says the values are believed correct to one unit in the last figure.
- **My reproduction (`check_newton1959.py`, solve_ivp, rtol 1e-12, atol 1e-13; output in
  `check_newton1959.out`).** CR3BP with mu = 1/82.45, start (x0, 0, 0, ydot0):
  - Closest approach to M over the half orbit equals rho_M0 in every row (for example 0.0121 and 0.0321).
  - The third y = 0 crossing has xdot = 0 within 5e-4 in all 11 rows, so the mirror symmetry closes the
    orbit. Twice the time of that crossing equals the printed period: differences are at most 1.0e-3
    (rows 1 to 11: 7.8926, 6.5224, 6.3842, 6.3336, 6.3020, 5.4302, 5.5728, 5.7574, 5.9748, 6.1082,
    6.1918). The two largest differences (1.0e-3) are the mixed and zero-velocity rows.
  - Full-period closure (`check_newton1959_extrema.py`): after the printed period the state error is
    3e-5 (retrograde x0 = 1.15), 1.3e-4 (direct x0 = 1.10) and 7e-3 (retrograde x0 = 1.00, where rho_M
    is 0.012 and ydot0 has five digits). So the tabulated data are good enough to restart a corrector.
  - Jacobi constant (2 Omega - v^2 form) at the start: retrograde -0.483, 0.206, 0.446, 0.663, 0.940;
    mixed 2.592; zero velocity 2.670; direct 2.635, 2.514, 2.340, 2.084. For mu = 1/82.45 the collinear
    points have 3.188 (L1), 3.172 (L2), 3.012 (L3) (`check_newton1959_lagrange.py`), so every orbit lies
    below them and the zero-velocity curve is open at all three necks. Newton does not print Jacobi
    constants.
  - Inertial sense at the start (ydot0 + x0): retrograde rows -1.331, -0.724, ... -0.298; direct rows
    +0.2025 (x0 = 1.05) and +0.151 (x0 = 1.20); the mixed row is -0.536. This matches the class names.
    For the zero-velocity orbit x0 = 1.01998 and ydot0 = -1.0200: the inertial velocity is zero, as it
    should be (the check is exact to the digits printed).
- **Unresolved printing discrepancy: perigee and apogee.** My distances from E do not match the printed
  columns (`check_newton1959_perigee.out`). Distances from the barycentre do match for the direct, mixed
  and zero-velocity rows, to about 0.001 (for example the direct x0 = 1.05 row: 0.1736 and 1.1459
  against printed 0.1732 and 1.1458), and for the retrograde x0 = 1.2 row (1.2424 against 1.242). So the
  printed columns look like distances from the origin, not from E, even though the text says rho_E. The
  retrograde rows with x0 = 1.0 to 1.15 do not fit either way: perigee 0.5207 against 0.4541 (from E) and
  0.4652 (barycentre); 0.1547 against 0.1106 for x0 = 1.15. I do not know if these are misprints or a
  different definition. The periods, rho_M0 and the closure do not depend on it. Do not cite perigee or
  apogee from this table without that caveat.
- Appearance (p.77): in xy there are only two shapes, direct and retrograde, with lobes in the left
  half plane (E side opposite M) larger than the right. In XY, retrograde orbits are precessing ellipses
  that expand and contract on alternate revolutions (Fig. 4, drawn for rho_M0 = 0.112129).

## 4. Comparison with Arenstorf 1963

Compared with `2026-10-06-digest-arenstorf-1963-amer-j-math-...` and the held AIAA J note.
- **The Newton family is an instance of Arenstorf's second-kind continuation.** Arenstorf's generating
  ellipse has a = (m/k)^(2/3) and synodic period 2 pi m. With a^(3/2) = 1/2, Newton's direct orbits have
  m = 1, k = 2 and the retrograde ones k = -2; my reading, since Newton does not use this notation.
  Both are symmetric about the line of the primaries, as in Arenstorf's theorem.
- **Where they differ.** Arenstorf's theorem is for small mu and excludes the finitely many
  eccentricities whose generating ellipse collides with the second mass. For x0 = 1.00 (both rows with
  rho_M0 = mu) the apocentre of the unperturbed ellipse is exactly at M's position at t = 0. That is
  the collision generating orbit. Newton's rows 1 and 6 are therefore at an Arenstorf exceptional value
  and carry the near-collision, second-species character. The rows with larger x0 are away from it and
  inside the theorem's scope for small enough mu. I did not compute how small mu must be.
- Arenstorf's digest says the "Apollo-type" orbits that pass near both masses are only numerical
  motivation in his paper. Newton 1959 is an earlier numerical source for the same idea, for mu = 1/82.45,
  but for an orbit class with a period near 2 pi, not for the figure-eight type of the Arenstorf orbit.
  The Newton orbits have one lunar flyby per period.
- Newton gives no existence proof and no stability result. He cites Poincare (1892, ch. 3) for the
  continuation argument and calls the orbits "periodic orbits of the first kind" in Poincare's terms
  (footnote 3). Newton cites Stromgren 1933 and Darwin 1911 for earlier numerical orbits near both masses
  (of a different kind) and Egorov 1958 for others.

## 5. Use in the open routes

- `#948` R4: earliest published family for the claim "ballistic periodic orbits passing close to both
  Earth and Moon exist at the Earth-Moon mass ratio". Use the 11 rows as a regression set for an Earth-Moon
  periodic-orbit continuation code. My closure test (section 3) shows the printed ydot0 values restart cleanly
  for the rows I closed over a full period. Treat rows 1 and 6 (near-collision) with care.
- Not a cycler paper: it has no transfer cost, no repeated Earth flybys at a specific altitude, no
  stability, no synodic repeat time other than the orbit period.
- Needs a literature caveat: Newton cites Egorov (1958) for other periodic orbits near both masses.

## 6. Citation mining

Held status checked with `ls cyclers_pdf/papers | grep` and CORPUS_INDEX.
- Stromgren 1933 (Bull. Astron. 9:87): HELD (`stromgren-1933-connaissance-actuelle-orbites-probleme-trois-corps-...`);
  digest `2026-10-06-digest-stromgren-1933-...` exists.
- Darwin 1911 (Scientific Papers vol. 4, mass ratio 10 to 1): not held; not on the wanted list. Low priority
  (historical, large mu).
- Egorov 1958 (in Blokhintsev et al., The Russian Literature of Satellites, Pt. 1, p.115): not held; not
  on the wanted list. **New candidate (medium):** it is the only other periodic family near both masses
  that Newton cites, and it is Russian satellite-era work close to the date of Newton's paper.
- Poincare 1892 (Les methodes nouvelles, vol. 1): not held; textbook. Arenstorf's digest already notes
  the same.
- Whittaker 1944 (Analytical Dynamics): not held; textbook.
- Arenstorf 1963: HELD, two files (`arenstorf-1963-existence-periodic-solutions-passing-near-both-masses-...`
  and `arenstorf-1963-periodic-solutions-restricted-three-body-analytic-continuations-...`).
- Nothing else to add to the wanted list.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
