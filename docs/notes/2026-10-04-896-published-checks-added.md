# #896: published checks added from the 2026-10-04 digests and older ones

Date: 2026-10-04. Task: turn every reproducible number printed in the recently digested papers
into a permanent test of the project's models ("use the knowledge from the new papers, always.
add the checks."). New test files only; no model, data or existing test was changed. A check that
fails is held as a strict expected failure (`xfail(strict=True)`) with its numbers in the reason.

Thirteen new files under `tests/core/` (items j and k, Gomez & Olle 1991 and Peng & Xu 2015, were
added later; addenda with the evidence were appended to the thirteen digests concerned). Each exits 0 (strict xfails count as passing); wall times are
with the project's default `-n 6`. `uv run mypy src tests` is clean (893 files), ruff check and
format are clean.

| File | Wall time |
|---|---|
| `test_bcr4bp_oshima_2022.py` | 2 s |
| `test_bcr4bp_rosales_2021_l2.py` | 7 s |
| `test_cr3bp_font_nunes_simo_second_species.py` | 11 s |
| `test_leiva_briozzo_2008_tables.py` | 11 s |
| `test_leiva_briozzo_2005_more.py` | 3 s |
| `test_er3bp_neelakantan_2022_more.py` | 1 s |
| `test_er3bp_mako_salamon_2025.py` | 15 s |
| `test_cr3bp_published_l1_l2_numbers.py` | 6 s |
| `test_ccr4bp_kumar_2021.py` | 25 s (one control alone is 12.9 s, above the 5 s per-test aim) |
| `test_er3bp_peng_2017.py` | 7 s |
| `test_cr3bp_anderson_kumar_2024.py` | under 1 s |
| `test_cr3bp_gomez_olle_1991.py` | 8 s (`-n 0`) |
| `test_cr3bp_peng_xu_2015_stability.py` | 4 s (`-n 0`) |

## Controls

| Paper | Table / page | Model | Asserted | Measured | Result |
|---|---|---|---|---|---|
| Oshima 2022 | Tables 2-5, p1332 (15 rows; Table 4 row 1 was already tested) | bcr4bp, paper's Table 1 constants | each printed y = 0 crossing closes after one Sun period from its printed Sun angle | 9.6e-10 to 2.7e-8 | pass |
| Oshima 2022 | Tables 2-5 | bcr4bp | consecutive crossings are one trajectory; y = 0 crossings fall at the times the Sun angles imply | 1.5e-9 to 6.7e-9; times within 1.6e-9 | pass |
| Oshima 2022 | controls | bcr4bp | Sun reversed, Sun a quarter turn off, Sun off, Andreu mu: all miss | 0.14-0.17, 0.40-0.48, 0.12-0.31, 1.4e-6 to 6.2e-6 | pass |
| Oshima 2022 | Fig. 8, p1331 (read off the figure) | bcr4bp | monodromy moduli: z0 families real pair near 1.06 / 0.94, vz0 families all unit | 1.06045 / 0.94299; others within 2e-13 of 1 | pass |
| Rosales, Jorba & Jorba-Cusco 2021 | Table 2, p7 | bcr4bp (theta_sun0 = pi) | multipliers of the L2 replacement orbit, by pseudo-arclength from L2 through negative Sun mass (turning point near eps = -0.065) | 1e-13 relative on all three, asserted at 2e-12 | pass |
| Rosales et al. 2021 | Fig. 2, Fig. 3, p6 | bcr4bp | the route goes into negative eps; the orbit loops L2 twice | as printed | pass |
| Rosales et al. 2021 | control | bcr4bp, Sun reversed | route starts into positive eps; far-side route gives another orbit | lambda_1 1.0577e6 | pass |
| Font, Nunes & Simo 2002 | Figs. 8-10 captions, pp138-139 | cr3bp | (phi, psi) is the fixed point of the n-encounter return map at mu = 1e-4, C_J = 2.8 | 1.5e-14 to 1.7e-13 | pass |
| Font, Nunes & Simo 2009 | Tables 3, 4, p155 (5 orbits) | cr3bp | fixed point and printed period | fixed point 3.3e-14 to 5.7e-11; period 1.4e-14 to 2.8e-11 | pass |
| Font, Nunes & Simo 2002/2009 | as above | cr3bp | "stability parameter" = trace of the return map (definition INFERRED, not in the papers) | 7e-10 to 3.8e-7 relative; 2002 Fig. 8 off by a factor 10 (exponent) | pass, Fig. 8 strict xfail at printed exponent and pass at E+6 |
| Font, Nunes & Simo | controls | cr3bp, cr3bp_regularized | Sundman-regularized passages agree; four wrong conventions miss by > 1e-3; literal one-period closure for the four least unstable orbits | 1e-11 (1e-9 within 1e-4 of the Moon) | pass |
| Leiva & Briozzo 2008 | Table 1 note, p234 | cr3bp | the section abscissa is L1 at mu = 0.0121505482 | 1.5e-10 (project mu: 1.6e-7) | pass |
| Leiva & Briozzo 2008 | Table 1, p234 (34 orbits) | cr3bp, mu = 0.0121505482 | each returns after its printed period | 3.1e-8 to 4.6e-4 (project mu: 2.9e-5 to 3.4e-2) | pass |
| Leiva & Briozzo 2008 | Tables 3, 4, p240 (24 arcs; 187A_t1 omitted) | qbcp, paper's mu | per the paper's p239 definition, each arc returns in position and velocity after tau = 5/2 or 7/2 Sun periods from t_i | position 7.8e-7 to 1.7e-4 | pass |
| Leiva & Briozzo 2008 | control | qbcp | the same arcs started half a Sun period later miss by 500-12,800 times more | 2.7e-3 to 1.0 | pass |
| Leiva & Briozzo 2008 | Table 5, p241, d_E (35 rows) | qbcp | minimum distance to the Earth's centre x 384,400 km | 26 within 1.4 km, all within budget | pass |
| Leiva & Briozzo 2008 | Table 5, d_M (35 rows) | qbcp | minimum distance to the Moon's centre | 28 within budget | 7 strict xfail |
| Leiva & Briozzo 2008 | control | qbcp | surface distance and pulsating-scale readings fail | ~1737 km; 5.9-860 km | pass |
| Leiva & Briozzo 2005 | Sect. 4.1, p365 | cr3bp | seed orbit returns after one Sun period | 4.7e-5 | pass |
| Leiva & Briozzo 2005 | Sect. 4.1 | cr3bp | localisation orbit h follows from x0, ydot | 1.2e-6 (budget 1.6e-6) | pass |
| Leiva & Briozzo 2005 | Sect. 4.1 | cr3bp | localisation orbit returns after its printed period 6.372441 | 3.5e-2 | strict xfail |
| Leiva & Briozzo 2005 | Sect. 4.3, p368-369 | qbcp | Earth distances 134,588 and 134,381 km | 5.2 and 3.3 km low | 2 strict xfail |
| Leiva & Briozzo 2005 | Sect. 4.3 | qbcp | difference of lunar distances (31 km printed) | 39.2 km | strict xfail |
| Leiva & Briozzo 2005 | Sect. 4.3 | qbcp | closest lunar approach is the start point | yes | pass |
| Leiva & Briozzo 2005 | Sect. 4.3, p369 | qbcp | stability parameters to 3 decimals | orbit 1: 5.5036 / 2.0632 (printed 5.496 / 2.069); orbit 2: 5.5401 (5.531), 2.0697 (2.070) | 3 strict xfail, 1 pass |
| Neelakantan & Ramanan 2022 | p7, p8 | cr3bp | the two circular halo states are periodic at the stated mu = 0.0122 | half-period residuals 0.031, 0.0037 | strict xfail |
| Neelakantan & Ramanan 2022 | p7, p8 | cr3bp | same states periodic at mu = 0.012277471 (found by solving; not printed) | 3.9e-12, 3.3e-12; closure 1.8e-10, 2.7e-11 | pass |
| Neelakantan & Ramanan 2022 | control | cr3bp | not periodic at the project's Earth-Moon mu | 0.054, 0.0061 | pass |
| Neelakantan & Ramanan 2022 | p7, p8 | cr3bp | ratio of the printed periods in days | 1.0522239 printed, 1.0515584 measured | strict xfail |
| Neelakantan & Ramanan 2022 | Table 8, p14, M4N2 Lyapunov | er3bp | closes from periapsis or apoapsis | 4.51, 2.0 | strict xfail |
| Mako & Salamon 2025 | p3-5 Defs. 1, 3; p8 | er3bp | instrument checks: speed scale consistent, low circular orbit returns after its Kepler period, circular start stable, 1.2 v_e unstable | 1.7e-7 | pass |
| Mako & Salamon 2025 | Fig. 7, p13 | er3bp | stable for f0 in [0, 177] and [216, 360) deg | all stable (327-331 deg checked at the lower end of the printed speed's rounding) | pass |
| Mako & Salamon 2025 | Fig. 7, p13 | er3bp | unstable for f0 in [178, 215] deg; transitions at 178 and 215 | stable at all 38 points | strict xfail |
| Mako & Salamon 2025 | Fig. 8, p15 | er3bp | stable for [0, 118] and [260, 360) deg | all stable | pass |
| Mako & Salamon 2025 | Fig. 8, p15 | er3bp | unstable for [119, 259] deg; transitions; launch point one speed step below the boundary | stable at all 141 points; model boundary 0.998 v_e against printed 0.982 v_e | strict xfail |
| Mako & Salamon 2025 | Fig. 7 | er3bp | launch point one speed step below the boundary at f0 = 0 | stable, then unstable | pass |
| Jorba, Jorba-Cusco & Rosales 2020 | Sect. 3.2, p13 | cr3bp linearisation at L1 | planar and vertical frequencies 2.33438585628816, 2.2688310655411 | 1.6e-13, 2e-14 at mu = 1/(1 + 81.300585) | pass |
| Jorba et al. 2020 | Table 1, p3 | cr3bp | at the Table 1 mu both misses are one mass-ratio error | within 2e-11 of the predicted 2.3e-9 | pass |
| Jorba et al. 2020 | control | cr3bp | project default mu misses by about 1.9e-8 | as stated | pass |
| Singh, Park & Howell 2026 | Table 2, p12 | cr3bp | L2 Lyapunov halo bifurcation 14.8319 d | 14.831874 d | pass |
| Singh, Park & Howell 2026 | Table 3, p12 | cr3bp | L2 Lyapunov axial bifurcation 18.7183 d | 18.718299 d | pass |
| Singh, Park & Howell 2026 | Table 4, p14 | cr3bp | L2 vertical axial bifurcation 19.2033 d | 19.203197 d (-1.03e-4 d) | strict xfail |
| Singh, Park & Howell 2026 | Tables 2-4 ER3BP brackets | cr3bp units | the time unit is the project's, to -7.8e-7/+3.2e-7; sidereal-month unit misses by 0.02 d | as stated | pass |
| Kumar, Anderson, de la Llave & Gunter 2021 (AAS 21-651) | Table 1, p8, p12 | ccr4bp constants | GM rows give the printed mu3 and mu_bar2 | 1e-15, only with the two moon GM rows interchanged | pass |
| Kumar et al. 2021 | p8, Fig. 1 | cr3bp limit | 3:4 orbit at the printed rotation number has C = 3.0041 and closest approach 22052 km | 3.0041057; 22051.69 km | pass |
| Kumar et al. 2021 | p8 | ccr4bp | torus at mu3 = 7.804102777055038e-5: invariant, closest approach 18721 km | 3e-10; 18721.37 km | pass |
| Kumar et al. 2021 | Fig. 2, p9 (read off the plot) | ccr4bp | stroboscopic multiplier about 7.18 and 6.915 | 7.1800, 6.9206 | pass |
| Kumar et al. 2021 | control | ccr4bp | Ganymede's synodic rate reversed misses the print | 21774.8 km | pass |
| Peng, Bai & Xu 2017 | Table 1, p4 (6 rows) | er3bp, e = 0 | corrected halo with half period N pi / M is the printed one | within 6.7e-7 | pass |
| Peng et al. 2017 | Table 2, p14 (8 rows) | er3bp, e = 0.2056 | corrected orbit within one printed digit | 3 rows within 1.3e-6 | 5 strict xfail |
| Peng et al. 2017 | Table 3, p14 | er3bp | largest monodromy eigenvalue of the 3 reproduced rows | 7.346e4, 8.777e4, 5.074e4 against -5.7936e5, 74,342, 43,176 | 3 strict xfail |
| Peng et al. 2017 | controls | er3bp | wrong resonance period; f0 = pi orbit started at 0; symplectic monodromy | as stated | pass |
| Kumar & Anderson 2024 (AAS 24-288) | p7 | cr3bp, mu = 3.54326e-5 | Uranus-Oberon C(L1) = 3.00454, C(L2) = 3.00450 within one last decimal below the computed value | 3.0045498, 3.0045025 (the print is the computed value cut, not rounded, for L1) | pass |
| Kumar & Anderson 2024 | control | cr3bp | convention with mu(1 - mu) added misses L1 | 4.5 last-decimal units | pass |
| Gomez & Olle 1991 II | Fig. 8(a), 8(b), p159 | cr3bp, mu = 1e-6 | printed C = jacobi_constant + mu(1 - mu); perpendicular crossing at the printed half period; closure after 2T; Broucke type | C to 1e-14; crossing time 4e-11, xdot 1e-10; closure 1.3e-9, 4.3e-10 | pass |
| Gomez & Olle 1991 II | p163 starting orbits A0, A1, B1 | cr3bp (Sundman) | perpendicular crossing at k pi and passage at O(mu) | xdot -1.1 to -1.26; passage 1e-5 to 3.5e-5 | strict xfail |
| Gomez & Olle 1991 I, II | Lemma 7, Thm 10 | cr3bp | computed k pi members near the printed ones pass at the first-order distance mu K0 | 9e-6 to 3e-5 relative | pass |
| Gomez & Olle 1991 II | Figs. 14, 16, 18, pp162-165 | er3bp, e = 0.5, velocity d/df | periapsis and axis crossing at f = k pi, first-order distance, closure (Fig. 14) | 1e-13 in f; 1e-5 relative; 9.3e-7 | pass |
| Gomez & Olle 1991 | controls | er3bp | d/dt reading and big-primary-at-origin reading fail | no crossing; 0.5-1.2 | pass |
| Gomez & Olle 1991 | Figs. 14, 18 | `genome.er3bp_periodic.correct_er3bp_periodic` | symmetric mode at tol 1e-5 accepts the printed orbit | moves 9e-12, 2e-11 | pass |
| Peng & Xu 2015 | Tables 1-5 e = 0 rows, pp294-300 | cr3bp, L1 halo of period 4 pi/5, monodromy to the fifth power | lambda_1 at mu 0.009, 0.015, 0.004; unit pairs; trivial pair | 1.311204e6, 1.596597e6, 9.31098e5; 0.958732 +- 0.284311i, 0.902069 +- 0.431592i | pass |
| Peng & Xu 2015 | Tables 1, 2, 4, 5 | cr3bp | printed 1/lambda_1 within one last digit | 2.8 to 28.8 units off the exact reciprocal | 4 strict xfail |
| Peng & Xu 2015 | Table 3 lambda_3 | cr3bp | printed unit pair at mu 0.004 | 0.989628 +- 0.143652i against printed 0.9021 +- 0.4316i (copied from Table 2) | strict xfail |
| Peng & Xu 2015 | controls | cr3bp | single-period monodromy, other mu, period 6 pi/7, the other family member of period 4 pi/5 all miss | as stated | pass |

## Findings (every strict expected failure)

1. **Leiva & Briozzo 2008 Table 5, d_M, seven rows.** The printed lunar minimum distance exceeds
   the one computed along the integrated orbit by 1.4 to 7.5 km (budgets 0.8 to 3.5 km), and by
   242 km for 032B_1_t4 (printed 725 km, computed 483.4 km, inside the Moon). In every row where
   the two differ by more than 0.6 km the printed value is the larger. The d_E column agrees in
   all 35 rows, so the length scale and centre-distance reading are right. Diagnosis: not
   determined; a sampled minimum in the paper (missing the true minimum near a fast lunar
   passage) would explain a one-signed excess, INFERRED.
2. **Leiva & Briozzo 2005.** The localisation orbit (x0 = 1.107569, ydot0 = -1.644251) misses its
   start by 3.5e-2 after the printed period 6.372441; the symmetric orbit at that h is at
   x0 = 1.110654 with period 6.370242. The Earth distances are 5.2 and 3.3 km below the printed
   ones; the lunar-distance difference is 39.2 km against 31 km printed; three of four stability
   parameters differ by 0.1 to 0.3 percent, which closure error (4e-6) cannot produce and the
   paper's mass ratio does not remove. Diagnosis: UNDETERMINED. The 2005 orbits do close in the
   module (existing test), so the module and the paper share the orbits but not all their derived
   numbers. Given the `#892` history this should be treated as a possible remaining difference
   between `core.qbcp` and the paper's quasi-bicircular model, not explained away on the paper's
   side. The digest's explanation of the 7 km lunar-distance difference (the pulsating length
   scale) does not hold: that scale gives 32.1 km for the difference but would need a lunar
   radius of 1637 km for orbit 1.
3. **Neelakantan & Ramanan 2022.** The two circular-problem halo states are periodic at
   mu = 0.012277471 (eleven digits the same from both states, independently), not at the paper's
   stated 0.0122 (the Table 8 elliptic rows do close at 0.0122, existing test). The printed day
   periods do not share one time unit (implied 4.3795 and 4.3768 d). "M4N2 Lyapunov" is still not
   reproduced: a differential correction from the printed state reaches only orbits far away, or
   a nearby one of period 2 pi, not 4 pi.
4. **Mako & Salamon 2025.** The paper's weak-stability classifier, implemented on `core.er3bp`
   from Definitions 1 and 3, reproduces the stated stable bands and agrees that both launch points
   lie on the first boundary, but finds no unstable band at [178, 215] deg (Fig. 7) or [119, 259]
   deg (Fig. 8). The model's least stable start in Fig. 7 is near f0 = 330 deg (just before
   perihelion). For Fig. 8 the printed speed (0.982 v_e) is far from the model's boundary
   (0.998 v_e), and the orbits drawn in the figure reach several million km, which the model gives
   near 0.998 v_e only: the printed speed is probably not the one the figure was computed with.
   Four readings of the return section and of the periapsis sign were tried; none changes this.
   The speed scale of the paper is the two-body speed times (1 - e), unexplained in the paper.
5. **Singh, Park & Howell 2026, Table 4.** The L2 vertical-family axial bifurcation is at
   19.203197 d in `core.cr3bp` against the printed 19.2033 d (about one unit in the last printed
   digit). Integrator and tolerance do not move it; the paper's own ER3BP brackets pin its time
   unit to the project's within 1e-6, excluding a unit change; mu would have to change by 4e-6.
   Cause not determined; a coarse bifurcation estimate in the paper's continuation is one
   possibility, INFERRED. The other two periods agree to the printed digits.
6. **Font, Nunes & Simo 2002, Fig. 8 caption.** The stability parameter printed as 2.338645E+7
   is the trace of the return map at exponent 6 (all seven digits agree). An exponent slip in the
   caption; Figs. 9 and 10 agree with their exponents.
7. **Peng, Bai & Xu 2017.** Five of eight Table 2 rows do not converge to the printed orbit with
   a damped single-shooting corrector from the printed state (three stall at residuals 3e-4 to
   5e-3; two converge to other orbits). With mu = 1.6601209e-7 and e = 0.205630 one of them (L1-1)
   converges to within 2e-7 of the print, so rounding of the "approximate" printed constants is a
   likely cause for that row. A better corrector might reproduce the others; this is not proof
   the rows are wrong. Table 3: the measured largest eigenvalues do not match their rows; read
   with the row numbering of the Fig. 17 captions the two L1 values match to 1.2 and 0.9 percent,
   the L2 value still does not. Nothing relabelled is asserted.

8. **Gomez & Olle 1991 II, p163 starting orbits.** From the printed (x, ydot) the axis crossing
   comes 1.35e-6 to 1.46e-6 after k pi with xdot -1.1 to -1.26, and the orbit passes 1e-5 to
   3.5e-5 from the small primary instead of about 1e-7. The computed members whose periapsis lies
   on the axis at t = k pi differ by (dx, dydot) = (+5.89e-6, -5.89e-6), (+9.08e-6, -9.08e-6),
   (-1.537e-5, +1.538e-5), with x + ydot (Theorem 10's family parameter) equal to the printed sum.
   INFERRED: the printed parameter is right and the printed x is off. Those members pass at the
   paper's first-order distance, so the theory reproduces in the project model.
9. **Peng & Xu 2015.** The printed 1/lambda_1 misses the exact reciprocal of the reproduced
   lambda_1 by 2.8 to 28.8 units of the last digit, and the paper's own products
   lambda_1 x (1/lambda_1) are 0.99990, 1.00008, 1.00046 and 1.000046: a numerical error in the
   paper's small eigenvalue. Table 3's lambda_3 column is a copy of Table 2's (confirmed in the
   PDF and now measured: 0.989628 +- 0.143652i at mu 0.004).

## Code findings (not fixed; no source was edited)

- `genome.er3bp_periodic.correct_er3bp_periodic`: at its default tol 1e-10 it raises
  ConvergenceError on all three printed Gomez & Olle elliptic orbits (line search stalls at 5.3e-8,
  9.8e-6, 4.4e-6), because its residual is taken at a fixed f where d(xdot)/dt is about 1e8 near a
  1e-7 periapsis. Its Radau full-period check returns 4.2e-5 and 2.2e-2 against its own 1e-5 bound
  and only logs a warning, so `genome.er3bp_continuation`, which uses it, would accept such orbits
  silently. The plain integrators resolve the 1e-7 passages about as well as a Sundman integration
  (periapsis angle 1e-7 to 5e-6 at 1e-13); the fixed-time residual is what is ill-conditioned.
- `search/er3bp_periodic` holds only a canonical-momentum converter and a monodromy helper.

## Other findings (passing checks that changed what we know)

- **Jorba et al. 2020's RTBP L1 frequencies** match the linearisation to 1e-13 at Earth/Moon
  mass ratio 81.300585, not the 81.300587 behind their Table 1 mu; the 1e-9 miss is exactly that
  mass-ratio difference (INFERRED: the paper does not state which value it used for the RTBP).
- **Leiva & Briozzo 2008 use mu = 0.0121505482**, identified from the printed section abscissa;
  at the project's mu every Table 1 orbit closes 12 to 1100 times worse. The existing Table 2 test
  in `test_qbcp.py` uses the module default mu; two of its rows close worse at the paper's mu,
  nine better.
- **Kumar et al. 2021 Table 1 has the Europa and Ganymede GM values in each other's rows**; read
  swapped, it reproduces the paper's own mass ratios to 16 digits. The registry-built default
  Jupiter-Europa-Ganymede system's Ganymede synodic rate is 1.7e-4 off the printed periods,
  because it uses Kepler's law on the registry semi-major axes (the same class of defect as the
  Uranian mean motions in `#894`). `core/ccr4bp.py` centres Ganymede's circle on the
  Jupiter-Europa barycentre where the paper centres it on Jupiter; the effect on the closest
  approach is 0.015 km. This is the first published control of the model itself (`#893`).
- **Leiva & Briozzo 2008 Table 4 row 187A_t1** is printed with xdot and y both -0.0371102305
  (page image and text layer agree). Solving for y alone with the printed xdot gives -0.0175001
  and a return to 1.4e-5, so the printed y is the misprint. The row is omitted from the tests.
- The 2009 Font-Nunes-Simo "stability parameter" is not defined in either paper; it agrees with
  the trace of the linearised return map to 1e-9..4e-7 (INFERRED definition).

## Digest transcription errors and omissions found against the PDFs

- `2026-10-04-digest-font-nunes-simo-2009-second-species-numerical-study.md`: orbit 12g's phi is
  given as 0.98857090907035220406767336927899; the table (p155) prints
  0.988570907035220406767336927899 (an inserted "09"). The digested value closes only to 2.0e-9,
  the printed one to 1e-12. A test pins the difference.
- The Font-Nunes-Simo digests' (phi, psi) recipe was marked INFERRED; the papers fix it: frame
  rotated by pi from the project's (big primary at +mu), psi is the angle of the synodic velocity
  from the x axis, phi is counter-clockwise from +x on the circle of radius mu^(2/5) about the small
  primary, points are exits from the disk, and C_J includes the constant mu(1 - mu).
- `2026-10-04-digest-rosales-jorba-jorba-cusco-2021-bcp-halo-like-tori-l2.md` section 6 says the
  L2 control discriminates the Sun's phase; it cannot (a phase change is a time shift and the
  continuation starts at time-independent L2; measured: the same multipliers to 1e-13 at
  theta_sun0 = 0 and pi). Only the sense is discriminated. Its "1e-8 relative is a fair target"
  is too cautious: 1e-13 is reached from the cyclic block spectrum (the paper's own method).
- `2026-10-04-digest-jorba-jorba-cusco-rosales-2020-bicircular-l1-vicinity.md` section 2.2
  leaves the 2.3e-9 frequency miss unexplained; it is the Earth/Moon ratio 81.300585 (above).
- `2026-10-04-digest-singh-park-howell-2026-aas-26-654-l2-families-intermediary-models.md`
  section 7.1 proposes a 1e-4 relative tolerance for the bifurcation periods, which would hide
  the Table 4 discrepancy; half a unit of the last printed decimal is used.
- `2026-10-04-digest-neelakantan-ramanan-2022-two-impulse-me-halo-ertbp.md` says the printed
  0.0122 must be used literally for any reproduction; the circular halos need 0.012277471.
- `2026-10-04-digest-leiva-briozzo-2008-rtbp-to-qbcp-periodic-transfer-orbits.md` recommends
  epoch reading (A) first; reading (B), the state at QBCP clock time t_i, is the one that closes
  (already so in `test_qbcp.py`). Its Control B expects about 1e-7 for 013_t3; measured 5e-6.
- `2026-10-04-digest-mako-salamon-2025-weak-stability-transition-region-er3bp.md`: does not flag
  that A7 and the A6 numerator (p19) have the sign opposite to what A4 (p18) requires (a paper
  slip; the printed sign changes the boundary speed by less than 5e-4); presents only the fixed
  return half-line where the proof of Proposition 1 (p6) reads as a rotating one; does not
  mention that p3 says "mean anomaly" while the sweeps are in true anomaly; does not note that
  the Figure 8 speed is inconsistent with the plotted orbit sizes; attributes R_E = 6378 km to
  Remark 1 (it is in the main text, p8). Its 1.22426 km/s variant for Figure 7 escapes at every
  f0, so the p15 v_p remark is the likelier slip, not the printed 1.20446.
- `2026-10-04-digest-oshima-2022-...`: no error; Tables 1 to 5 checked at 250 dpi.
- `2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md` treats e = 0.2056 and mu = 1.660e-7 as exact;
  the paper calls both approximate.
- Not checked against PDFs, flagged during the item (i) inventory: the Restrepo-Russell 2018 digest
  gives the Mars-Sun mu as 3.2271676e-06 (should be about 3.227e-7); Blazevski-Ocampo eq. 15 m2 and
  m3 look swapped; Canales-Howell-Fantino Fig. 20 says 3.0028 where the text says 3.0024;
  Guido-Efthymiopoulos give both M_J = 0.00096 and mu = 0.001.
- Code and test docstrings (not edited): `search/variational_ccr4bp_torus.py` says Kumar et al.
  print no rotation number or energy to reproduce, but they print both and both reproduce;
  `#761`'s 22035.8 km comes from rounding C (22051.69 km with the period from the printed rotation
  number); `tests/core/test_crnbp.py` cites the tri-circular problem paper as AAS 23-257, the
  paper's own number is 23-201.

## Not done

- Oshima 2022 prints no multipliers; Fig. 8 moduli are checked to figure-reading precision.
- Leiva & Briozzo 2005's speeds at closest approach (2023, 2024, 2443, 2446 m/s) are not tested:
  the velocity unit is printed only as "~1024 m/s" and the frame of the speed is not stated;
  neither unit candidate gives all four to the metre per second.
- The item (i) inventory read about 115 of 171 digests in the relevant sections; about 20 mission
  or tour digests were screened by keyword only. Three older papers were implemented (Kumar et al.
  2021 for ccr4bp, Peng et al. 2017 for er3bp, Kumar & Anderson 2024 for cr3bp at the Uranus-Oberon
  mass ratio), the last a small check.
- Process: the repository's pytest `addopts` include `-n 6`, so the five agents of the first wave
  each ran pytest with six workers at times, above the brief's limit of six processes in all.
  mypy on a single test file needs `MYPYPATH=src` (the full `mypy src tests` does not).

## Printed numbers in older digests that still have no test

Four-body moon model (`ccr4bp`):
- Kumar et al. 2021, `2026-07-23-digest-kumar-2021-europa-ganymede-ccr4bp-resonant-orbits.md`,
  p12, Figs. 14-17: Jupiter-Ganymede-frame continuation, omega 3.111756 to 3.116809 at mu_bar2
  1.0015e-5 to 2.506370e-6; 3:2 Ganymede torus to mu_bar2 2.5265e-5. Needs the frame transform
  and a second continuation.
- Kumar et al. 2023 (AAS 23-397), `2026-07-23-digest-kumar-2023-secondary-resonance-overlap-ccr4bp.md`:
  omega range [2.032685, 2.0405], secondary resonances 11/34, 12/37. Needs the 4:3 torus-family
  machinery.
- Kumar & Anderson 2024, `2026-07-27-728-anderson-kumar-2024-oberon-mmr-survey-digest.md`:
  the L1/L2 Jacobi constants are now tested; the Titania-perturbed statements ("no tori below
  C = 3.007714") need torus machinery.
- Aryan & Fitzgerald 2024, `2026-07-26-710-digest-aryan-fitzgerald-2024-jovian-pccfbp.md`,
  Tables 1-2: rotation numbers at C = 3.0034 and 3.0044; no torus state printed.
- Blazevski & Ocampo 2012, `2026-07-27-732-blazevski-negri-baresi-foundational-papers-digest.md`,
  Table 1: ICs in a model with Jupiter fixed (Gm1 = 1), so they cannot close in `ccr4bp.py`.

Many-moon model (`crnbp`):
- Gilliam thesis 2025 Ch. V, `2026-07-26-digest-gilliam-bettinger-2024-crnbp-jovian.md`,
  Tables 5-6: Lagrange-box sizes and displacements for Uranus-Oberon, Saturn-Titan,
  Saturn-Enceladus, Sun-Earth, Sun-Ceres. Perturber phases not printed; window-dependent.
- Baresi, Owen & Scheeres (tri-circular problem), `2026-07-27-722-baresi-owen-scheeres-tri-circular-problem-digest.md`,
  Table 1: constants only; tori and Floquet results figure-only.

Elliptic model (`er3bp`):
- Peng et al. 2017: the five Table 2 rows and Table 3 eigenvalues held above.
- Martinez-Cacho et al. 2025, `2026-06-25-digest-planar-retrograde-ERTBP-2025.md`, Table 6:
  constants only.

Ballistic capture (`wsb`): Topputo & Belbruno 2015, `2026-07-22-digest-topputo-belbruno-2015.md`,
Tables 3-5: Sun-Mars mu = 3.2262081094e-7 and capture dV against periapsis radius; needs their
capture definition.

Bicircular model (`bcr4bp`): Onozaki et al. 2017, `2026-06-30-digest-onozaki2017-tube-4body.md`,
Table 1 constants.

Circular three-body model (`cr3bp`, already well covered):
- Singh et al. 2021 (NRHO and L1 halo digests, 2026-06-17), Tables 1-2: 15-digit ICs with
  stability indices.
- Howell 1984, `2026-06-25-digest-howell-1984-halo-orbits.md`, Tables I-III: halo family at
  mu = 0.04 (ICs, T/2, C, stability index).
- Leiva-Briozzo "Full Atlas" and Broucke 1968, `2026-07-28-744-broucke-leiva-barrabes-earth-moon-lineage-digest.md`:
  16-digit ICs for families 357 and 037; libration points at mu = 0.012155099, family F orbits;
  Barrabes-Gomez resonant C values (sign convention unresolved).
- Frauenfelder-Koh-Moreno 2023 (C = 3.003571774, T0 = 2.1215, eigenvalues), Howell-Davis-Haapala
  2012 (C = 3.068621, 3.000785, 3.17212), both in `2026-07-28-744-symplectic-invariant-periapse-maps-digest.md`.
- Moreno et al. 2024 (`2026-07-27-728-moreno-...-bifurcation-graphs-digest.md`, Appendix B, DPO
  branch at C = 3.00109352); Canales-Howell(-Fantino) 2021/2023 (Jd = 3.00754, Lyapunov
  C = 3.0061, Table 1); Miceli 2024 (C_J = 1.75598, 0.950382, 0.963141); Gomez et al. 2004
  (`2026-06-30-digest-gomez2004-spatial-rtbp.md`, Sun-Earth map C values); Koon-Lo-Marsden-Ross
  2002, Ross-Scheeres 2007, Grover-Ross 2009 (C values); Rawat 2024/2025 (resonance widths, not
  state reproductions).

Constants only: Restrepo-Russell 2018 (Table 1, 24 system mu values; Mars value flagged above),
Pergola 2007 (Uranian mu and e), Rosengren et al. 2026 (lunar resonance a and T), Gurfil 2007,
Bond-Allman 2021.

No model in the repository: Llibre-Martinez-Simo 1985 (Hill problem, N(inf) = 5.1604325,
M(inf) = 2.1330587), Aydin-Batkhin 2025 (Hill family g ICs), Pozzi et al. 2026 (J2-perturbed L1/L2
and halo C ranges).

## Commits

b97188d3, 080b85ea (Oshima); 4d70e952 (Rosales 2021); a18cf2bc, 553daa8b (Font-Nunes-Simo);
0a82441d, d9bfce56 (Neelakantan-Ramanan); ebbb8d47, f0e757e1 (Leiva-Briozzo 2008); 833f6ce7,
0e47d2c7 (Leiva-Briozzo 2005); 299c757e, c97b72fc, f6dabee6 (Jorba 2020, Singh 2026); b752b68b,
fd2eeabf (Mako-Salamon); d6686881, 7fa199e7, 0c503032 (Kumar 2021, Peng 2017); 672660af (Kumar &
Anderson 2024); e08394d8, c02ec829 (Gomez & Olle 1991); 681c9b1a, ac7bf462, 640248a1 (Peng & Xu
2015); 8a8e9f4d, d959d203 and the follow-up (this note, digest addenda, ledger).

Every strict expected failure was run once with `--runxfail`: each fails on its asserted
comparison (an `AssertionError` with the numbers quoted above), none by an exception.
