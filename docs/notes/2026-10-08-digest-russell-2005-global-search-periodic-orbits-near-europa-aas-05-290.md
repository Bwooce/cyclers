# Digest: Russell 2005, "Global Search for Planar and Three-Dimensional Periodic Orbits Near Europa" (AAS 05-290) (#960)

R. P. Russell (JPL), "Global Search for Planar and Three-Dimensional Periodic Orbits Near Europa", paper **AAS 05-290**,
2005 AAS/AIAA Astrodynamics Specialist Conference, Lake Tahoe CA (p.1 footnote). The file name carries the JPL clearance
number 05-2188. The conference paper has no DOI.
- Journal form: Russell, R. P., J. Astronaut. Sci. 54(2):199-226 (June 2006), **doi 10.1007/BF03256483** (Crossref:
  title, single author Russell, volume, issue and pages confirmed). **Not held**: `ls cyclers_pdf/papers | grep -i russell`
  and CORPUS_INDEX show no copy. It is row 80 of the `#730` acquisition backlog
  (`docs/notes/2026-07-27-730-acquisition-backlog-master-list.md`). It is not a row of the `#960` wanted list. Because the
  journal form is not held, no identity check was possible. This digest covers the conference form only.
- File: upload `a9dbda5e-05-2188.pdf`, 26 pp., md5 **da06ead940c1e4890b6c7193f393b330**. PDF metadata: Word document
  "TahoePeriodicOrbits3newmarg.doc", Acrobat Distiller 5.0.5, created 2005-08-03. The text layer is digital and good, so
  no OCR is needed.
- Proposed filename: `russell-2005-global-search-planar-three-dimensional-periodic-orbits-near-europa-AAS-05-290.pdf`.
- How I read it:
  - I read the full text layer.
  - I read the page images at 150 dpi for pp.2-3 (Eqs. 1-2, Table 1), p.5 (Eqs. 11-13), p.12 (Table 2, counts, family
    names) and p.25 (Table 3).
  - I read Table 3 row by row on 400-dpi strips.
  - From the held Restrepo-Russell 2018 PDF, I read pp.5, 7, 8 and 12 on page images.
- Deliverables (scratch folder `dg39-russell05/`):
  - `russell-2005-europa-periodic-orbits-tables.yaml`: Table 3, 76 rows, with a conventions block.
  - `witness-comparison.tsv`, built by `compare_witnesses.py`.
  - `check_russell2005.py` and `check_russell2005.out`.

## 0. Verdict

**It is the Europa-only, three-dimensional parent of the Restrepo-Russell 2018 grid-search database. Its tables are
fully reproducible in our CR3BP.**
- **What it is.** A brute-force grid search in the Jupiter-Europa CR3BP for two kinds of symmetric periodic orbit: axi-
  symmetric orbits and doubly symmetric orbits. Both kinds start on the x-axis and perpendicular to it. Each candidate
  found on the grid is corrected with a first-order differential corrector.
  - Size: "over 10 billion grid points", 616,942 converged orbits, about 950 CPU hours (p.12).
  - Planar and 3D results. The paper names five simply periodic (N = 1) planar families and maps them to
    Broucke/Hénon/Robin-Markellos names.
  - Only 76 representative 3D orbits are printed (Table 3). The full 616,942-orbit file was "archived in an electronic
    text file" (p.20). The file is not public, and I found no link to it in the paper.
- **What it gives the project.**
  1. The lineage of the Restrepo-Russell 2018 database. Restrepo-Russell p.7 says their global search "is similar to
     the one described by Russell (2006)". Section 5 below lists what carried over and what changed.
  2. A sourced family-name crosswalk for Europa (sec. 3):
     - Circle-Egg = Broucke H1 = Robin-Markellos g1 = 2018 LPO1.
     - Egg-Diamond = H2 = g2 = 2018 LPO2/DPO.
     - DRO = Broucke C = Hénon/Robin-Markellos f.
     - L1, L2 = 2018 LL1, LL2.

     This is the "Europa family identity" item that consumers ask for. The 2018 side of the crosswalk is from the
     2018 paper p.12 (page image).
  3. Seventy-six 3D Europa periodic orbits with ICs, period, Jacobi constant, symmetry type and stability. All 76
     reproduce in our integrator under the paper's own constants, with the symmetric crossing at T/4 or T/2 as printed
     (sec. 4). They can serve as positive controls for any 3D symmetric-orbit corrector at a small mu (Europa
     2.528e-5). The 2018 database is planar only, so it cannot give 3D controls.
- **Gate relevance.**
  - `#976` (Newton-Arenstorf continuation in mu to Jupiter-Ganymede and Saturn-Titan, gated by the 2018 database): this
    paper adds the method lineage and Europa's planar family names. It holds no Ganymede or Titan data.
  - `#956` R9 is an Earth-Moon exterior-realm cell. This paper is Europa-only and has nothing at Earth-Moon. The
    footnote on p.20 says the solutions "can generally be scaled to be almost-valid for any RTBP system with small mass
    ratios using a ratio of the corresponding L1 distances". That applies only to small mu, and Earth-Moon mu is 0.012.
  - "X4": in the `#938` note, X4 is the DA transfer-map enumerator at Pluto-Charon and Titan. I found no Europa X4.
    Nothing here bears on that cell. The caller should say which X4 is meant.
  - No cycler and no flyby sequence. This is a bound-orbit census, not a cycler source.
- **PROPOSALS only (no catalogue or code change):**
  - (a) A test fixture with 3-5 Table 3 rows (for example 1609237 D N=1, 1480596 D N=2, 1348961 D N=3 retrograde,
    1472520 S N=7, 1606263 S N=4 "squid orbit") as positive controls for a 3D symmetric-orbit corrector.
  - (b) If the project adds the `cr3bp_family` row field proposed in the 2018 digest, record the crosswalk in sec. 3.
  - (c) Ask the author for the archived 616,942-orbit file only if a 3D Europa census is ever needed.

## 1. Model (pp.2-3, page images)

- Jupiter-Europa CR3BP, **centred at Europa**. "The x-axis is fixed opposing the direction to Jupiter". So Jupiter is
  at x = -1 DU, z is along Europa's orbital angular momentum, and y completes the right-handed set.
- Eq. (1): x'' = 2v + (x+1-mu) + kappa(x+1) + mu/rE^3, y'' = -2u + y + kappa y, z'' = kappa z, with
  kappa = -(1-mu)/rJ^3 - mu/rE^3 and mu = GmE/(GmJ+GmE).
- Eq. (2): J = (x+1-mu)^2 + y^2 + 2(1-mu)/rJ + 2mu/rE - u^2 - v^2 - w^2.
- Table 1 (from Bagenal et al. 2004):
  - Jupiter-Europa distance 6.709 x 10^5 km.
  - Europa radius 1560.70 km.
  - GmE 3202.72 km^3/s^2.
  - GmJ 1.2668654 x 10^8 km^3/s^2.
  - Footnote a: mu = 2.528002607976249 x 10^-5.
  - DU = 670900 km and TU = 48822.04433066813 s (p.3).
- **Point-mass CR3BP only.** There is no Jupiter oblateness, no eccentricity and no ephemeris.
- Our checks (`check_russell2005.out` PART 1):
  - mu and TU recomputed from Table 1 reproduce the printed values exactly (relative difference 1e-16; TU to the
    printed digits).
  - The collinear points are 13,744.5 km and -13,559.3 km. Table 2 footnote c prints xL2 = 13,744 km and
    xL1 = -13,559 km.
- Table 3 is printed in **dimensional units**: km, km/s, days, km^2/s^2. Our arithmetic: J_nondim = J_printed / 188.835509.
- Stability (pp.4-6):
  - Monodromy eigenvalues {lambda1, 1/lambda1, lambda2, 1/lambda2, 1, 1}.
  - k_i = lambda_i + 1/lambda_i, from a1 = 2 - trace(Phi) and a2 = {a1^2 + 2 - trace(Phi Phi)}/2 (Eqs. 12-13).
  - Stable means k_i is real and |k_i| <= 2.
  - New instability index, Eq. (11): rho = max(|lambda1|, |1/lambda1|, |lambda2|, |1/lambda2|), so rho = 1 when stable.
  - Eq. (14) builds the full monodromy matrix from the half-period STM (axi-symmetric) or the quarter-period STM (doubly
    symmetric), using the reflection matrices L and K from Robin-Markellos.

## 2. Method (pp.6-11)

- **Symmetries (pp.5-6).**
  - The equations are invariant under {t -> -t, y -> -y}. This is the xz-plane mirror.
  - They are also invariant under {t -> -t, y -> -y, z -> -z}. This is a 180-deg rotation about the x-axis.
  - Axi-symmetric orbit: two perpendicular crossings of the x-axis, separated by T/2.
  - Doubly symmetric orbit: a perpendicular x-axis crossing at t = 0 and a perpendicular xz-plane crossing at T/4. The
    paper notes that it also has the xy-plane symmetry and is axi-symmetric as well.
  - Planar case (Hénon's "N periodic symmetry"): a perpendicular x-axis crossing is also a perpendicular xz-plane
    crossing, so the half period is T/2.
- **Correctors.**
  - Axi-symmetric, Eq. (20): fix x0, vary (v0, w0), and target z = u = 0 at the Nth xz-plane crossing.
  - Doubly symmetric, Eq. (22): fix x0, vary (v0, w0), and target u = w = 0 at the Nth crossing.
  - Planar, Eq. (23): one unknown v0, target u = 0.
  - Each corrector stops the integration on the Nth y = 0 crossing, which removes the time unknown. All three are
    first-order Newton correctors, after Robin & Markellos 1980.
- **Grid (p.8, Fig. 6 algorithm).**
  - Search space {x0, v0, w0, N}, with the start (x0, 0, 0, 0, v0, w0).
  - Only v0 > 0 and w0 > 0 are searched. The w0 < 0 orbits follow from the xy-plane symmetry. The v0 sign reversal
    gives the reciprocal crossing, so it would duplicate orbits (p.9).
  - x0 > 0 starts direct. x0 < 0 starts retrograde.
  - Halo orbits (xz-symmetric but with no x-axis crossing) and asymmetric orbits are excluded (p.9, and the incomplete
    sentence on p.2).
- **Near-solution test (pp.9-10, Figs. 4-5).**
  - On each constant-x0 slice of the (v0, w0) grid, the corrector is called only where one grid step changes the sign
    of both target quantities: (z_f, u_f) for an axi-symmetric orbit, or (u_f, w_f) for a doubly symmetric orbit.
  - Six step directions are tested.
  - The seed is the midpoint between the two linear-interpolated zero crossings on the step with the smallest
    normalised distance.
  - Planar case: a sign change in u_f along v0 at w0 = 0.
- **Uniqueness.** "Check for repeated solutions and false classifications (N and symmetry)" (Fig. 6). The paper gives
  no explicit tolerance for duplicates.
- **Integration and convergence (p.11).**
  - Integrator: variable-step RK7(8) with a stopping condition.
  - Grid shooting: Jacobi constant held to 8 significant digits.
  - Corrector: about 13 digits in J, and the norm of the periodicity constraints converged to **1e-10**.
- **Grid regions, Table 2 (p.12, Nmax = 16).** x0 counts are numbers of equally spaced values.

  | Region | x0 range | x0 values | v0 range (km/s) | v0 values | w0 range (km/s) | w0 values |
  |---|---|---|---|---|---|---|
  | Planar Retrograde | -150,000 km to surface | 2,000 | 0.0-7.0 | 80,000 | 0 | 1 |
  | Planar Direct | surface to L2 | 1,000 | 0.0-2.0 | 40,000 | 0 | 1 |
  | 3D Retrograde I | -150,000 to -50,000 km | 2,000 | 2.0-7.0 | 2,500 | 0.0001-1.5 | 750 |
  | 3D Retrograde II | -50,000 km to surface | 3,000 | 0.0001-2.5 | 1,250 | 0.0001-3.0 | 1,500 |
  | 3D Direct | surface to L2 | 1,000 | 0.0001-2.0 | 1,000 | 0.0001-2.0 | 1,000 |

  - Our arithmetic: the products are 1.6e8 + 4.0e7 + 3.75e9 + 5.625e9 + 1.0e9 = 1.0575e10 (x0, v0, w0) nodes. This is
    consistent with "over 10 billion grid points".
  - With Nmax = 16 crossings to T/2 or T/4, a full orbit can have up to 32 (axi-symmetric) or 64 (doubly symmetric)
    crossings. This matches the "32 times" and "64 times" on p.2, so there is no conflict with Nmax = 16.
- **Counts (p.12).**
  - "over 10 billion grid points are evaluated and 616,942 solutions are found using approximately 950 hours of total
    computer time on Linux machines with 3066 MHz processors".
  - "Approximately 5% of the solutions, or 30,040, are found to be stable in a linear sense and have close approaches
    above Europa's radius. Of those, 19,383 are planar."
  - Subsurface solutions are kept down to 100 km below the mean radius.
  - The abstract says "over 600,000".

## 3. Family nomenclature (pp.12-14, p.12 page image)

- The paper introduces descriptive names for the **five simply periodic (N = 1) planar families**:

  | Russell 2005 | Equivalent named on the page | Restrepo-Russell 2018 (p.12, page image) |
  |---|---|---|
  | Circle-Egg (direct; grazing circular orbit near J = 568.75 km^2/s^2, becomes an egg with its base toward Jupiter; stable almost to impact) | "Brouke's [sic] H1 family and Robin and Markellos' g1 family" | g1 curve = LPO1, plus 2B-LPO at high J (LPO1* where unstable) |
  | Egg-Diamond (direct; grazing egg with its base away from Jupiter, J = 566.18, x0 = 11,500 km; becomes circle, then diamond, then two loops, then impact; stable over most of the egg part up to max J near 566.215) | "Broucke's H2 family or Robin and Markellos' g2 family" | g2 curve = DPO (unstable part, b_h >= 2), LPO2 (stable), LPO2* |
  | L2 family (direct; around the far-side collinear point; unstable throughout) | none | LL2 |
  | L1 family (around the interior collinear point) | none | LL1 |
  | DRO (retrograde; grazing circle, then a vertically aligned near-ellipse; stable "to well beyond 150,000 km") | "Broucke's family C and Henon's and Robin and Markellos' family f"; DRO term from Lam & Whiffen (ref. 24) | DRO |

- Branches named for discussion (p.13): circular-branch and lower egg-branch (Circle-Egg), and upper egg-branch and
  diamond-branch (Egg-Diamond). There is a saddle structure on the circular-branch near J = 567.318 km^2/s^2 (printed
  "567.3 18").
- Hénon comparison (p.13): Hénon's Hill-problem g and g' families intersect. Here, Fig. 7b shows a clear gap between
  Circle-Egg and the Egg-Diamond family. The page prints "Diamond-Circle" at this point, which I take to be the
  Egg-Diamond family. "Henon's g' family is egg shaped on both ends, and the g family transitions from a circle to the
  diamond shapes, opposite from what is seen in the RTBP". This is left to future work. Restrepo-Russell 2018 p.12
  makes the same point: the two prograde curves "connect in the Hill model".
- **There is no family numbering** for the 616,942 orbits and no family ID column. The 3D families are described only
  by how they bifurcate from the planar branches. Most 3D direct solutions lie on the low-J side of the planar saddle
  near 567.25 km^2/s^2. Examples:
  - The "sitting swan" volume plot (Fig. 11).
  - A central, highly inclined stable region near 7,000 km and 70 deg pseudo-inclination.
  - The "squid orbit" near x0 = 11,000 km and i = 87 deg (Fig. 16u). It is stable, but its stability island is tiny.
    This is Table 3 row 1606263: x0 = 1.10342143E4 km, inc 86.7 deg, S, N = 4, rho = 1.
- The Table 3 ID is a running solution number (1216180-1609237). It is not a family label.
- Symmetry labels in Table 3:
  - D = doubly symmetric.
  - S = axi-symmetric. Footnote a says "(A)", but the column prints S. Integration confirms that S means axi-symmetric
    (sec. 4).
- Science result (p.20): highly inclined stable direct 3D orbits are "relatively abundant" compared with the retrograde
  side at similar distances. The paper says this runs against "the general attitude of mission planners that retrograde
  orbits are always more stable".

## 4. Table 3 (p.25): transcription and checks

- 76 rows: 44 stable (rho = 1, increasing x0, Figs. 15-16) and 32 unstable (rho 2.93 to 6.38E6, increasing, Figs.
  17-18). There are 47 D rows and 29 S rows. Every row has w0 > 0, so no row is planar.
- Witnesses (`witness-comparison.tsv`):
  - A: my 400-dpi image reading of every row.
  - B: the publisher text layer.
  - C: raw tesseract psm 6 at 400 dpi.
  - A = B on all 988 cells.
  - C = A raw on 862 cells. The others are glyph noise: "N)" for S, "EO" for E0, a space inside a number, "," for ".".
  - Four cells still differed after normalisation, and I re-read them on the image: 1357937 k1 1.69E0, 1506466 hmin
    4.24E3, 1507698 k1 -4.37E-1, 1609237 k2 -7.48E-1.
- Checks (ours, `check_russell2005.out`; DOP853 rtol 1e-12, atol 1e-13, Eq. (1) model, Table 1 constants):
  - **J**: Eq. (2) from (x0, v0, w0) reproduces the printed J to half a unit of the last digit in **76/76** rows. This is
    an independent check, because J is printed separately from the state.
  - **Pseudo-inclination** atan(w0/v0) matches the printed value to 0.05 deg in **76/76** rows.
  - **rho from k1, k2** (Eq. 12) matches the printed rho within 1 percent in **76/76** rows. Rows 1-44 all have
    |k_i| <= 2 and rho = 1.
  - **Symmetry, N and T** (all 76 rows): the Nth y = 0 crossing falls at T/4 for every D row (47/47) and at T/2 for
    every S row (29/29), within 1e-6 of T.
    - The symmetry residual is below 2.1e-7 nondim: u, w for D rows; z, u for S rows. The other symmetry's residual is
      1e-3 to 0.2, so the test tells the two types apart.
    - This confirms N, the period convention (T = 4 tf for D, T = 2 tf for S) and the meaning of S.
  - **Full-period return** from the printed digits, stable rows: |X(T) - X0| median 3.5e-8 and max 2.0e-5 nondim (max
    at 1449566, N = 15, 29.8 d). Unstable rows grow as expected, up to 1.3e-2 for rho ~ 3.5e5. Jacobi drift is at most
    2e-12.
  - **hmin**: the minimum altitude over the symmetric arc rounds to the printed value in 73/76 rows. For the three
    others the paper prints a higher value: 1261682 (ours 34,139.9 km, printed 3.42E4), 1542144 (525.4, printed 5.26E2)
    and 1606263 (648.2, printed 6.49E2). A minimum sampled at integrator steps would give this. The values are kept as
    printed. (A first pass over the full period gave a false miss on 1313405, rho 6.4e6, because the error grew; the
    symmetric-arc pass removes it.)

## 5. Relation to the Restrepo-Russell 2018 database

Sources: the 2018 paper, pp.5, 7, 8 and 12 (page images; p.4 and the Table 1 row on p.6 from the text layer), and `data/sources/restrepo-russell-2018/README.md`.

| Item | Russell 2005 | Restrepo-Russell 2018 |
|---|---|---|
| Systems | Jupiter-Europa only | 24 systems (Europa is system 502) |
| Dimension | planar and 3D | planar only ("limited to the planar case", p.5) |
| Symmetry | axi-symmetric and doubly symmetric, starting on the x-axis | planar axisymmetric (two perpendicular x-axis crossings) |
| Frame | Europa-centred, Jupiter at x = -1 | same: primary at x = -1, `x(1)` measured from the secondary (README sec. 3 item 6) |
| Units in output | dimensional (km, km/s, days, km^2/s^2) | nondimensional, with `JC-3` = J - 3 |
| Europa constants | DU 670,900 km; GmJ 1.2668654e8; GmE 3202.72; mu 2.528002608e-5 | DU 671,100 km (2018 Table 1, p.6, text layer); GmJ 126,687,000; GmE 3202.74; mu 2.52800922e-5 (README sec. 4, file header) |
| Grid | {x0, v0, w0, N}, v0 > 0, Nmax 16, fixed Table 2 boxes | global {x0, ydot0, N} with x0 within 5 x_L1, ydot0 from 0 to +/- ydot0_max (p.7), Nmax = 10 (p.8), plus local refined searches around LL1, LL2 and DRO (sec. 3.2) |
| Seed test | double sign change on a 2D (v0, w0) slice, six directions, interpolated midpoint | final xdot near 0 (p.7) |
| Corrector | first-order Newton (Robin-Markellos), Eqs. 20, 22, 23; tolerance 1e-10 | full second-order trust region (Conn-Gould-Toint), first- and second-order STTs, forward/backward symmetric propagation (p.7) |
| Stability | k1, k2 = lambda + 1/lambda; rho = max abs(lambda) | b_h, b_v = lambda + 1/lambda (Eq. 6, p.5); same definition, with no factor 1/2. The 2018 file adds `stabA` |
| Names | Circle-Egg, Egg-Diamond, L1, L2, DRO | LPO1/2B-LPO (g1), LPO2/DPO/LPO2* (g2), LL1, LL2, DRO, plus composed families (Hg, Hb, Hm, QDRO, resonant) |
| Counts | 616,942 (19,383 stable planar) | 45,845 rows in Europa_D; 137,961 Europa rows across 9 files (README) |
| Availability | archived text file, not public | public, Apache 2.0 |

- **Carried over:**
  - The idea of a brute-force grid search seeded by sign changes, then a local corrector.
  - Terminating on the Nth axis crossing, which removes the time unknown.
  - The Europa-centred frame with Jupiter at x = -1.
  - The lambda + 1/lambda stability index.
  - The Broucke H1/H2 = g1/g2 identification of the two prograde curves.
- **Changed:**
  - Many systems instead of Europa only.
  - Planar only instead of 3D.
  - Nondimensional output.
  - A second-order trust-region corrector instead of first-order Newton.
  - Local refined searches around the generating families, which find sensitive "connecting" families that a global
    grid misses.
  - A family classification key built from crossing and rotation counters.
  - Updated Europa constants: DU differs by 200 km (0.03 percent) and mu by 2.6e-6 relative.
- **Our arithmetic**: because DU and mu differ, a 2005 orbit and its 2018 counterpart do not share identical
  nondimensional ICs. Any row-level cross-match must convert each orbit with its own paper's constants. There is no
  row-level overlap anyway: 2005 prints only 3D rows, and 2018 has only planar rows.

## 6. Print slips (as printed; not corrected)

- p.2: the sentence "In addition asymmetric solutions" breaks off.
- p.5: "1,5,713,14" is a run-together citation for refs. 1, 5, 7, 13, 14. The page also has "implored" for "employed".
- p.12: Table 2 footnote b says "Europa radius defined in Table 2". It means Table 1.
- p.12: the "# v0" header of Table 2 carries a superscript d, but there is no footnote d.
- p.12: "the grid search is performed in the (x0, y0) space". The grid is in (x0, v0).
- p.12: "Brouke's H1" should be Broucke's.
- p.13: "Diamond-Circle families" is probably the Egg-Diamond family.
- p.21: "Figures A1-A4 and accompanying data in Table A1". These are Figures 15-18 and Table 3.
- p.22: the caption of Figure 16 repeats "Part I".
- p.25: Table 3 footnote a says "(A)", but the column prints S.
- Ref. 8 cites "Jefferys" and the text cites "Jeffreys".

## 7. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX / the `#960` wanted list.

| Ref | Item | Held? |
|---|---|---|
| 1 | Broucke 1969, AIAA J. 7(6):1003-1009 | held (`broucke-1969b-...-aiaa-j-7-1003...`) |
| 2 | Broucke 1968, JPL TR 32-1168 | held |
| 3 | Hénon 2003, CMDA 85:223 | held (tables in `data/sources/henon-2003-hill-families-tables.yaml`) |
| 4 | Hénon 1969, A&A 1:223 | held |
| 5 | Hénon 1974, "Vertical stability ... restricted problem", A&A 30:317 | not held (Hénon 1973 A&A 28:415, part I, is held) |
| 6 | Poincaré, New Methods (1993 ed.) | not held |
| 7 | Szebehely 1967 | held |
| 8, 9 | Jefferys 1965 AJ 70:393; 1966 AJ 71:566 | not held (only Jefferys 1974 held; 1971 Atlas is wanted row 59) |
| 10 | Bray & Goudas 1967, AJ 72:202 | not held (wanted row 53 lists Bray & Goudas 1967 Adv. Astron. Astrophys. 5:71, a different item) |
| 11 | Goudas 1963, Icarus 2:1 | not held |
| 12 | Zagouras & Markellos 1977, A&A 59:79 | not held |
| 13 | Papadakis & Zagouras 1993, Ap&SS 199:241 | not held |
| 14 | Robin & Markellos 1980, Celest. Mech. 21:395 | not held. This is the corrector and the g1/g2 naming source, and it is also cited by Restrepo-Russell 2018. Candidate for the wanted list |
| 15, 16 | Lara & Scheeres 2002 JAS 50:389; Lara & Peláez 2002 A&A 389:692 | not held |
| 17 | Lara, Russell & Villac, "On Parking Solutions Around Europa", AAS 05-384 (same conference) | not held; companion paper (Europa halo and 3D stability). Candidate |
| 18 | Kazantzis & Goudas 1975, Ap&SS 32:95 (3D grid search) | not held. Method ancestor; candidate |
| 19 | Bagenal et al. 2004, "Jupiter" (constants) | not held (book) |
| 20, 21 | Battin 1987; Nayfeh & Balachandran 1995 | not held (textbooks) |
| 22 | Ocampo 1996 PhD thesis | not held |
| 23 | Howell 1983 PhD thesis (halo) | not held |
| 24 | Lam & Whiffen 2005, AAS 05-110 (Europa DROs, "Red Sea" plot) | not held; **wanted row 28** (NTRS record only) |

- Wanted-list suggestions (proposal): Robin & Markellos 1980 and Lara, Russell & Villac AAS 05-384. The journal form of
  AAS 05-384 is probably Lara, Russell & Villac 2007, JGCD 30(2), "Classification of the distant stability regions at
  Europa". I have not checked this on Crossref. Kazantzis & Goudas 1975 is lower priority.
- The journal form of this paper (doi 10.1007/BF03256483) is on the `#730` backlog (row 80). Get it only for attribution
  or for an identity check: Table 3 may differ in the journal.

*Filed as `cyclers_pdf/papers/russell-2005-global-search-planar-three-dimensional-periodic-orbits-near-europa-aas-05-290-jpl-clearance-05-2188.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/russell-2005-europa-periodic-orbits-tables.yaml`.*
