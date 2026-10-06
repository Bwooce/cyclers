# Digest: Hoelker & Winston 1968, "A Comparison of a Class of Earth-Moon Orbits with a Class of Rotating Kepler Orbits" (#960 batch 30)

R. F. Hoelker and B. P. Winston (NASA Electronics Research Center, Cambridge MA), NASA TN D-4903,
November 1968 (text dated August 1968), NTRS 19690001139. No DOI. 65 PDF pages (report pp. 1-56 plus
covers and appendix).
- Filed as `cyclers_pdf/papers/hoelker-winston-1968-comparison-class-earth-moon-orbits-rotating-kepler-orbits-nasa-tn-d-4903-ntrs-19690001139.pdf`.
  - OCR copy (filed): `ntrs_19690001139-ocr.pdf`, md5
    797afda1fe93025ff73d183d197bc658.
  - Original NTRS scan (supplied file `ntrs_19690001139.pdf`): md5 4aaf76c7061c3df7eeade7a4a0381683, 65 pp.
- How I read it:
  - Text: the OCR layer, all pages.
  - Numbers: every figure caption on report pp. 19-45 that carries a velocity or a period was read
    on the 200 dpi page image (PDF pp. 21-22, 26-34, 36, 38-42, 44-47).
  - The OCR is wrong in many captions. Examples: Fig. 23 reads -2.4687 in the OCR and -2.4657 on the
    image; Fig. 121 P/2 reads "3477" in the OCR and 34.77 on the image.
  - Captions of reference and collision orbits on report pp. 21-23, 33, 35, 41 (Figs. 35-40, 78-81,
    83-87, 108-112) are from the OCR only. No number below depends on them, except the collision-orbit
    labels quoted as context.
  - The transcription is in `hw1968_captions.csv` (65 rows, with the printed page of each caption).
  - The check script is `hw1968_check.py`, with output in `hw1968_check_output.txt`. A scan of Figs. 44
    and 82 is in `hw1968_scan_44_82.py`, with output in `hw1968_scan_44_82_output.txt`.
- Wanted-list row 27 ("Same lineage; `#948` R4"). Removed in batch 30.

## 0. Verdict

**This is a tutorial atlas of one planar series at mu = 1/80. It prints 39 three-body initial
conditions with a period or collision time. 36 are periodic: 34 plain, the Fig. 94 libration orbit and the Fig. 48
Earth-collision orbit. I reproduced every one. None of
them is an Earth-Moon cycler in the sense of a close pass at both bodies in every period.**

- **Set-up.** Planar circular restricted problem with mu = 1/80 (p.13).
  - Earth is at x = -0.0125 and the Moon at x = +0.9875.
  - Every orbit starts at x = 1.4875, y = 0, xdot = 0. That is 1.5 units from Earth, half an Earth-Moon
    distance beyond the Moon. Only ydot varies, from -2.47 to +0.3.
  - Each three-body orbit is paired with the Kepler orbit (mu = 0) that starts at x = 1.5 with the
    same ydot.
  - Periodic orbits are drawn to the first perpendicular x-axis crossing. P/2 is printed.
  - The appendix integrates in Arenstorf's regularised variables (ref. 3, AJ 68(8) 1963) by power series.
- **Labels.** Each three-body orbit carries n*, the mean motion of the Kepler orbit it resembles, so
  n = k/l with synodic period 2 pi l. This is exactly Arenstorf's (1963) commensurability a = (m/k)^(2/3)
  (held digest `2026-10-06-digest-arenstorf-1963-amer-j-math-...`). So the periodic orbits are
  numerical members of Arenstorf's Kepler-continuation families at a finite mu, each identified only
  by one point.
- **Reproduction: every periodic caption reproduces at mu = 1/80.** I ran 37 rows (the 35
  periodic ones with a P/2, plus the near-circular Figs. 29 and 117) with a 1-D
  differential correction on ydot: xdot = 0 at the y = 0 crossing nearest the printed P/2, using
  DOP853 with rtol = atol = 1e-12.
  - Corrected minus printed ydot is between 2e-8 and 7.1e-5 in 36 of the 37 rows. The corrected P/2 matches the printed P/2 to its last digit,
    or within 0.01. The libration orbit, Fig. 94 (ydot = -0.847844, P/2 = 4.1952), matches P/2 to
    four decimals.
  - The exception is Fig. 44 (n* = -3/2, ydot = -1.568, P/2 = 6.30). Its perpendicular crossing is a
    pass 0.0057 units (2,176 km) from Earth's centre, at a speed of 18.7 units. At the printed ydot,
    the crossing is perpendicular to within 0.07 deg. The exact root is at ydot = -1.57254 with
    P/2 = 6.308. This is consistent with a 3-decimal caption on a near-collision arc.
  - Fig. 82 (ydot = -0.848798, P/2 = 29.70) is extremely sensitive. xdot at the half-period crossing
    changes by about 5e4 per unit of ydot. The root is at -0.8487999, which is 2e-6 from the printed
    value, with P/2 = 29.685. The paper itself warns that "the correction is concerned with the fifth
    significant digit".
  - **Fig. 89 is printed as ydot = 0.848358, without the minus sign.** -0.848358 reproduces P/2 = 9.973
    exactly (correction -1e-7). +0.848358 would break the ordering of the series: the neighbours, Figs. 88 and 90, are
    -0.848370 and -0.84835. So the missing minus is a printing slip.
  - Fig. 48 is "periodic collision with Earth", ydot = -1.500, n* = +-3/2. It falls radially onto Earth
    at t = 2.0, 6.3 and 10.6. My integration gives closest approaches of 134-170 km to the point
    mass, at those times. The paper prints "2.0; 6.3; 10.6". The Kepler twin, Fig. 51, prints 2.04
    and 6.12.
- **Positive control (Kepler captions, 26 rows).**
  - Vis-viva on each printed ydot (inertial v = ydot + 1.5) gives the captioned n to within 2.5e-4.
  - The integrated perpendicular crossing at mu = 0 falls at pi*l, to the printed digits.
  - The circular captions also check: -2.3165 and -0.68350 (Figs. 32, 119) are within 3.4e-6 of
    v_circ = sqrt(1/1.5).
  - This confirms both my reading of the captions and the sign convention.
- **Cycler verdict.** I took the lunar Hill radius (mu/3)^(1/3) = 0.161 (about 61,900 km at
  L = 384,400 km) as the test for a lunar pass. The numbers below are the minimum distances over one
  period. Distances in km use L = 384,400 km, which the paper does not state.
  - **Orbits with a lunar pass inside the Hill radius, and their closest Earth approach:**
    - Figs. 55, 57: closest Earth approach 0.187 / 0.221 (72,000 / 85,000 km); lunar pass 1,176 / 35,400 km.
    - Figs. 72-73, 75-77: Earth 0.33-0.56 (127,000-215,000 km); lunar pass 12,300-47,600 km.
    - Figs. 82, 89, 94, 96-99: Earth 0.55-0.75 (211,000-287,000 km); lunar pass 21 km to 1,895 km from
      the point mass. Fig. 94's pass is at 0.00162 (624 km), at the half-period crossing. Several pass
      inside the Moon's physical radius of 1,738 km (Figs. 82, 89, 94, 96, 97 and 98, and Fig. 55 at
      1,176 km). The model has point masses.
  - **Earth-grazing orbits:** Fig. 44 (2,176 km, inside the physical Earth) and Fig. 48 (collision).
    Their closest lunar approach is 0.326 and 0.41 units, so they have no lunar encounter.
  - **No orbit combines a close Earth pass with a lunar pass.** The Earth test is mechanical. The
    catalogue's Vaquero rows use a LEO-GEO insertion band of 6,558-42,164 km for the perigee
    (catalogue line 56634). Every orbit with a lunar pass has its closest Earth approach at 71,700 km
    or more, so all fail that band. Several orbits enclose both bodies, or
    alternate between an Earth-satellite phase and a lunar-satellite phase every period. The paper
    discusses this for Figs. 71-77 and calls Fig. 84's neighbour "a periodic orbit that encloses both
    Earth and Moon". These are Earth-Moon periodic orbits that "visit" both bodies. They are not
    perigee-and-flyby cyclers.
  - Also, mu = 1/80 is not the Earth-Moon value (0.01215). Every orbit would need continuation in mu
    before it could be compared with catalogue rows.
- **Use for the project:**
  1. **A sourced control set for a symmetric-orbit corrector at mu = 1/80.** It has 34 periodic
     ICs with P/2, reproduced here. The orbits span retrograde (n* < 0), Earth-collision, lunar-capture,
     L2 libration (Fig. 94, "translunar libration point, designated L1 in reference 1") and direct
     (n* = +5/8 ... +1/4) regimes. It is cheap to adopt.
  2. **Seeds for `#948` R4 (Earth-grazing second-species orbits).**
     - Fig. 48 is a published periodic Earth-collision orbit at finite mu. Fig. 44 is a perpendicular
       crossing 2,176 km from Earth.
     - Both are on the n* = 3/2 Kepler collision line (the inertial-rest start, Kepler Fig. 51). This
       is the same case as Arenstorf's exceptional e values, but for the large primary.
     - **Figs. 89-99 cluster within 1e-3 in ydot of the lunar-collision orbit of Fig. 81
       (ydot = -0.8489, collision at t = 4.18).** Every one passes the Moon at t = 4.2, at 195 km to
       1,895 km. This is a published neighbourhood of a periodic orbit passing near both a lunar
       collision and, separately, the Fig. 48 Earth collision. But no single printed orbit is near both.
  3. **The paper says (p.36) that Fig. 90, a reference orbit that is not periodic as printed, would
     have the topology of a lemniscate (figure-eight) once completed to periodic shape and stripped of
     loops that enclose no mass.** That periodic completion is not printed, and I did not compute it.
     Fig. 89 nearby is periodic but is a different orbit. So the report has no printed finite-mu
     figure-eight. It is only a pointer for the catalogue's `arenstorf-em-figure8-1963` row (line
     9077), at mu = 1/80, with perigees of order 0.6 L, not Earth flybys.
- **Catalogue (PROPOSAL only).**
  - No row cites Hoelker or Winston. There is no grep hit for "hoelker", "winston" or "D-4903" in
    `data/catalogue.yaml`.
  - Arenstorf row (lines 9077-9230): `cr3bp.jacobi_constant`, `period_nd` and `state_nd` are null
    with `data_gaps` (lines 9154-9158). This paper does NOT fill them: its mu is 1/80 and none of its
    orbits is the Apollo-type figure-eight with an Earth perigee.
  - Grep inventory: "earth-moon", "casoliva" and "arenstorf" hit only the Earth-Moon rows. "em-"
    gives 243 hits, mostly Earth-Mars. The Earth-Moon rows are: 9077 Arenstorf; 9235 Genova-Aldrin;
    9469 Wittal; 48219-48756 Ross/Roberts-Tsoukkas; 48445 spatial (2,1); 48901-49164 Braik-Ross;
    55643-56483 Casoliva; 56588-57104 Vaquero.
  - Casoliva rows (lines 55643-56483), Vaquero rows (56588-57104), Ross/Braik-Ross rows (lines
    48219-49164), and Genova-Aldrin and Wittal (9235, 9469): no collision. Those
    are at mu = 0.01215 and have close lunar passes and low perigees. No Hoelker-Winston orbit is a
    member or a duplicate as printed.
  - Proposal: no new row. Optionally, cite this report in the Arenstorf row's `notes` as the first
    published numerical atlas of Arenstorf-type commensurability families at finite mu (n* labels),
    with the caveat that mu = 1/80.

## 1. Content (READ)

- **Kepler orbits in rotating coordinates (pp. 2-12).**
  - For n = k/l, the synodic period is 2 pi l. The synodic mean motion is (k - l)/l. There are |k|
    identical patterns per period, spaced 2 pi (k - l)/|k| apart. Negative-n families have no cusp.
  - n = 1 orbits are libration orbits about the "empty mass point" L at x = 1. They do not enclose
    the mass and do not intersect (Fig. 10; ref. 2, Hoelker 1967).
  - The apsidal velocity diagram (Figs. 15-16, 41) maps (x_R, ydot_R) to n.
- **Series synopsis (Figs. 17-22).** As ydot runs from -2.47 to +0.3, the series goes:
  near-parabolic retrograde, then Earth collision, then lunar collision at ydot = -0.849, then lunar
  satellite orbits, then direct orbits that spiral out.
- **Periodic three-body orbits read on the image** (ydot; P/2; n* if printed). Jacobi constants
  computed from the IC are in the check output.

| Fig | ydot | P/2 | n* | note |
|---|---|---|---|---|
| 23 | -2.4657 | 12.55 | -1/4 | C = -2.5004 |
| 24 | -2.4215 | 9.41 | -1/3 | |
| 25 | -2.3363 | 6.25 | -1/2 | |
| 29 | -2.3138 | 2.02 | | "circle-coordinated" |
| 31 | -2.1005 | 3.21 | -1 | |
| 44 | -1.568 | 6.30 | -3/2 | Earth pass 2,176 km |
| 48 | -1.500 | | +-3/2 | periodic Earth collision |
| 55 | -1.12884 | 6.80 | | lunar satellite |
| 57 | -1.0997 | 14.95 | | general orbit |
| 58 | -1.0604 | 12.21 | | lunar satellite |
| 64 | -0.930243 | 2.973 | +1 | Kepler-libration twin |
| 67 / 68 | -0.89875 / -0.88748 | 20.1 / 8.39 | | 7 and 3 lobes |
| 71-73 | -0.879860 / -0.875970 / -0.874621 | 25.98 / 21.71 / 18.25 | | "losing a loop" |
| 75-77 | -0.859495 / -0.85716 / -0.855166 | 13.49 / 10.26 / 13.33 | | Earth/Moon phase alternation |
| 82 | -0.848798 | 29.70 | | period ~ 9 months, one direct lunar revolution |
| 89 | (-)0.848358 | 9.973 | | sign missing in print |
| 94 | -0.847844 | 4.1952 | | L2-type libration orbit |
| 96-99 | -0.847327 / -0.847252 / -0.847091 / -0.846990 | 20.48 / 17.31 / 14.11 / 20.36 | | |
| 103, 106, 107 | -0.830758 / -0.827640 / -0.815682 | 16.12 / 12.88 / 9.460 | | "pigtail" group |
| 113, 114 | -0.73839 / -0.72640 | 27.0 / 17.4 | +5/8, +3/5 | |
| 117 | -0.71854 | 8.434 | | near-circular |
| 121-123 | -0.640597 / -0.632473 / -0.624343 | 34.77 / 28.48 / 22.26 | +5/11, +4/9, +3/7 | |
| 127-129 | -0.611820 / -0.57621 / -0.52845 | 16.01 / 9.62 / 12.70 | +2/5, +1/3, +1/4 | |

- "Best approach" (not periodic): Figs. 30 (-2.2455, off by 0.3 deg), 42 (-1.810), 50 (-1.21),
  118 (-0.6854, off by 0.8 deg). Lunar collisions: Figs. 39 (-2.2190), 43 (-1.750), 56 (-1.1260),
  81 (-0.8489). Earth collision: Fig. 36 (-2.22012). These are OCR-only for Figs. 36, 39 and 81.
- p.36: Figs. 90 and 91 bracket another lunar collision (second approach).

## 2. Citation mining

- Szebehely 1967, *Theory of Orbits*: HELD (`szebehely-1967-theory-of-orbits-restricted-problem-three-bodies-book.pdf`).
- Hoelker, R. F. (1967), "Numerical Studies of Transitions Between the Restricted Problem of Three
  Bodies and the Problem of Two Fixed Centers and the Kepler Problem", NASA TM X-1465: not held, and
  not on the wanted list. It is a new candidate, low priority. It is the source of the Kepler
  libration family and the "L1" naming used here.
- Arenstorf, R. F. (1963), "New Regularization of the Restricted Problem of Three Bodies", AJ 68(8)
  (Oct 1963): not held. It is not either of the two held Arenstorf 1963 items (the AIAA J 1:238 note
  and Amer. J. Math. 85:27), and it is not wanted-list row 24 (the IAC paper). It is a new candidate,
  low priority: it is a method source for both-primary regularisation, which `#948` R4 says it needs.
- The paper does not cite Huang 1962/1963 (wanted row 25), Newton 1959 (row 28) or Schwaniger 1963.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
