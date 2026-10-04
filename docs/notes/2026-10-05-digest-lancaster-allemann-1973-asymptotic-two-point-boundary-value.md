# Digest: Lancaster & Allemann 1973, "Numerical analysis of the asymptotic two-point boundary value solution for N-body trajectories"

J. E. Lancaster and R. A. Allemann, AIAA Journal 11(3):259-260 (March 1973), DOI 10.2514/3.50464. A two-page synoptic of AIAA Paper 72-49 (10th Aerospace Sciences Meeting, San Diego, 17-19 January 1972; submitted 1 March 1972, synoptic received 21 September 1972; "Full paper available from AIAA"). Filed in the private paper corpus as `lancaster-allemann-1973-numerical-analysis-asymptotic-two-point-boundary-value-solution-n-body-aiaa-j-11-259-doi-10.2514-3.50464.pdf` (2 pages; the text layer carries the prose and Table 1 in jumbled column order, so the table was read from the page image of p.260 and checked row by row).
Evidence tags: READ (page) is what the printed page says; COMPUTED is my own arithmetic of 2026-10-05 (1 nautical mile = 1.852 km; mu = 0.0121529529; 1 Earth-Moon distance = 384,400 km); INFERRED is my reading across sources.
Related: `docs/notes/2026-10-04-digest-breakwell-perko-1974-second-order-matching.md` (the second-order matching this paper tests), `docs/notes/2026-10-04-digest-perko-1967-error-estimation-singular-perturbation.md` (why the error can only be calibrated).

## 1. What was compared (READ pp.259-260)

- The method: second-order asymptotic (matched expansion) solutions of the two-point boundary value problem for a negligible-mass particle under one primary and N-2 secondaries, "formulated to solve the boundary value problem without iterations", with all unknown parameters as analytic functions of the boundary conditions; "Computation times are comparable to the standard conic approximations with considerable increase in accuracy" (abstract). Outer expansion r = r0 + mu r1 + mu^2 r2 + O(mu^3) (eq. 2), inner expansion about the k-th body with R_k0 the two-body hyperbola and R_k2 a definite integral (eq. 6), matched in the overlap domain; the definite integrals are done by Gaussian quadrature. It builds on Lagerstrom & Kevorkian (1963), Breakwell & Perko (1966), Lancaster (1970 AIAA 70-1060, Moon-to-Earth), Carlson (1969), and the full report Lancaster, MDC G2748 (February 1972, 2 vols; also NASA CR 1973), which "discusses solutions of the following type: Earth-to-moon, interplanetary, midcourse (both lunar and interplanetary), and one- and two-impulse moon-to-Earth".
- The zero-order solution is a Lambert problem between the initial position r(t0) and the Moon's position at the prescribed pericynthion time; the asymptotic solution then gives first- and second-order corrections to the Lambert velocity at t0 (Fig. 1).
- Test: each asymptotic solution (first order and second order) was compared with numerical integration (a CDC 6500 program) of the same boundary-value problem; the printed errors are the differences at pericynthion between the prescribed boundary conditions and what the integrated trajectory started from the asymptotic initial velocity achieves (the sign convention of the errors is not stated). Which model the integration used for the Earth-Moon cases (restricted three-body or with the Sun) is not stated; the title says N-body and the cases range from N = 4 (moon-to-Earth) to N = 7 (Earth-to-Mars), so the Earth-to-Moon runs may include more than the restricted problem (INFERRED; not stated).
- Printed in full for one family only: Earth-to-Moon trajectories, Table 1. For the other types the text gives a few extreme values only (below).

## 2. Table 1, transcribed (READ p.260, image checked)

Prescribed at pericynthion: radius p_M = 1200 nautical miles (COMPUTED: 2,222 km, 0.00578 Earth-Moon distances, 0.476 mu in Earth-Moon units), inclination i_M = -35 degrees ("the negative sign indicating approach from under the moon"). Columns: case, initial radius r(t0) in nautical miles, true anomaly f(t0) in degrees, flight time t_pM in hours, order of solution, and the errors in the time of pericynthion (minutes), pericynthion radius (nautical miles) and inclination (degrees).

| Case | r(t0) (n mi) | f(t0) (deg) | t_pM (hr) | order | error in t_pM (min) | error in p_M (n mi) | error in i_M (deg) |
|---|---|---|---|---|---|---|---|
| 102 | 3,544 | 10 | 80 | 1 | 36 | -92 | -0.05 |
| 102 | | | | 2 | -8 | -299 | -0.33 |
| 111 | 12,850 | 118 | 79 | 1 | 36 | -82 | -0.18 |
| 111 | | | | 2 | 8 | -70 | 0.19 |
| 112 | 40,162 | 149 | 75 | 1 | 36 | -71 | -0.16 |
| 112 | | | | 2 | 13 | -9 | 0.14 |
| 113 | 63,707 | 157 | 70 | 1 | 37 | -62 | -0.14 |
| 113 | | | | 2 | 15 | -2 | 0.14 |
| 114 | 98,489 | 164 | 60 | 1 | 38 | -50 | -0.10 |
| 114 | | | | 2 | 16 | -12 | 0.12 |
| 115 | 125,011 | 167 | 50 | 1 | 38 | -44 | -0.07 |
| 115 | | | | 2 | 17 | -32 | 0.08 |
| 116 | 146,533 | 170 | 40 | 1 | 39 | -40 | -0.04 |
| 116 | | | | 2 | 19 | -54 | 0.06 |
| 117 | 164,463 | 172 | 30 | 1 | 40 | -42 | -0.02 |
| 117 | | | | 2 | 20 | -83 | 0.04 |

Printed prose figures not in the table: for the interplanetary midcourse solution "second-order errors in pericenter time, radius, and inclination were as small as 10^-1 sec, 10^-1 naut mile, and 10^-3 degrees, respectively, starting from midcourse points along a reference 244-day Earth-to-Mars transfer"; the full interplanetary solution and the two-impulse moon-to-Earth solution "resulted in errors somewhat larger than anticipated. In several cases the second-order errors were larger than the corresponding first-order errors (as in Cases 102, 116, and 117 in Table 1)"; computation times "from 1.7 sec for a moon-to-Earth trajectory (N = 4) up to 6.0 sec for a complete Earth-to-Mars trajectory (N = 7)" on a CDC 6500.
The mass ratio is not printed (Earth-Moon, INFERRED), nor the Moon-relative speed at the encounter, nor the number of cases beyond the eight shown (case numbers 102 and 111-117 suggest more were run; "Cases 102 and 111-117 show the effect of moving the initial position away from the Earth").

## 3. Derived values (COMPUTED from the table)

| Case | r(t0) (km; Earth-Moon distances) | first-order radius error (km; percent of 1200 n mi) | second-order radius error (km; percent) | second/first | first / (mu^2 L) | second / (mu^2 L) |
|---|---|---|---|---|---|---|
| 102 | 6,563; 0.017 | 170; 7.7 | 554; 24.9 | 3.25 | 3.00 | 9.75 |
| 111 | 23,798; 0.062 | 152; 6.8 | 130; 5.8 | 0.85 | 2.67 | 2.28 |
| 112 | 74,380; 0.193 | 131; 5.9 | 16.7; 0.8 | 0.13 | 2.32 | 0.29 |
| 113 | 117,985; 0.307 | 115; 5.2 | 3.7; 0.2 | 0.03 | 2.02 | 0.07 |
| 114 | 182,402; 0.475 | 93; 4.2 | 22.2; 1.0 | 0.24 | 1.63 | 0.39 |
| 115 | 231,520; 0.602 | 81; 3.7 | 59.3; 2.7 | 0.73 | 1.44 | 1.04 |
| 116 | 271,379; 0.706 | 74; 3.3 | 100; 4.5 | 1.35 | 1.30 | 1.76 |
| 117 | 304,585; 0.792 | 78; 3.5 | 154; 6.9 | 1.98 | 1.37 | 2.71 |

with mu^2 L = 56.8 km, mu^2 |ln mu| L = 250 km, mu^3 L = 0.69 km, mu^3 |ln mu| L = 3.0 km, mu^3 |ln mu|^2 L = 13.4 km (the scales of the Breakwell & Perko digest, section 9).

## 4. Findings of the paper (READ p.260)

- Dependence on the start position (cases 102 and 111-117): first order: the time-of-flight error is nearly constant (36 to 40 min); the radius and inclination errors improve as the start moves outward "until a point around 150,000 naut miles from Earth is reached. Case 117 then shows a slight increase in the radius error." Second order: the time error increases more rapidly with distance (8 to 20 min); "The pericynthion radius error, however, is initially reduced, dropping from 299 nautical miles for Case 102 down to 2 naut miles for Case 113. It then increases, becoming larger than first-order."
- Mass and boundary-condition dependence: "smaller values of mu should give better results. This was verified by comparing lunar and interplanetary results"; "for a fixed value of mu the accuracy of the asymptotic solution (particularly the second-order) is also dependent on the boundary conditions of the specific problem. Not only do these boundary conditions affect the magnitude of the errors but cases arise where the second-order error is larger than first-order."
- Their explanation: the cases with large second-order errors were those where "the first-order correction to the Lambert velocity at t0 was quite large"; the function r1(t), initially zero, grows rapidly; "In the second-order solution the effect of r1(t) enters quadratically and when r1(t) is large the integration of the quadratic term over the entire trajectory results in a large second-order correction. This breaks down the assumption of uniformity in Eq. (2) and leads to excessive second-order errors. ... asymptotic expansions which are initially convergent may diverge after n terms. The results obtained thus far indicate that the derived form of the asymptotic N-body solution may, for certain choices of the boundary conditions, begin to diverge after two terms."
- Conclusion on use: "the results thus far indicate that determination of midcourse velocity corrections is the best application. In addition, interplanetary applications result in significantly more accurate solutions than lunar applications."

## 5. Does it give the calibration `#899` needs?

Partly, with strong limits.
- What it gives: the only printed measurement in the set of how far a second-order matched solution at the Earth-Moon mass lies from an integrated trajectory, for a flyby at a periapsis of order mu (1200 n mi = 0.48 mu a_moon, the O(mu) regime of the digests of Perko 1976 and of Breakwell & Perko). For first order the radius error is 74 to 170 km (3.3 to 7.7 percent of the 2,222 km radius), equal to 1.3 to 3.0 times mu^2 a_moon, with the time error 36 to 40 min and the inclination error 0.02 to 0.18 degrees. This is of the same scale as the formal O(mu^2) error (mu^2 L = 57 km up to mu^2 |ln mu| L = 250 km), so the first-order scale of the Breakwell & Perko digest is confirmed to within a factor of about 3 (INFERRED comparison; coefficient 1.3 to 3.0 with respect to mu^2 L, or 0.3 to 0.7 with respect to mu^2 |ln mu| L).
- Second order does not show the expected gain: the formal scale is mu^3 |ln mu|^2 L = 13 km, and only cases 112 to 114 are near it (17, 3.7 and 22 km, i.e. 0.07 to 0.39 mu^2 L); cases 102, 116 and 117 are worse than first order (554, 100 and 154 km) and 111 and 115 only slightly better. So at mu_EM and this class of trajectory the second-order result has an error of 4 to 554 km, about five to eight hundred times the mu^3 L scale at the extremes. The reason given by the authors (large first-order Lambert correction, quadratic growth) is a condition on the boundary-value problem, not on the flyby: it applies to the global trajectory from the start point, so for a matched seed in a cycler (whose first-order correction may be small or large) the relevant diagnostic is the size of the first-order correction, which the paper recommends examining.
- Limits of the evidence: one flyby class (an Earth-to-Moon transfer arriving at 1200 n mi, -35 degrees, in 30 to 80 h), one pericynthion radius, one inclination, eight cases, one error per quantity (no scatter), the encounter speed and the number of cases not printed, the model and the sign convention not stated, boundary-value (not initial-value) errors, and no repeated flybys or cyclers. It cannot replace the calibration experiment proposed in the digests of Breakwell & Perko 1974 and Perko 1967 (matching seeds against Casoliva Table 3 flybys at the stated encounter speeds); it can serve as a prior: the first-order seed error at mu_EM is of the order of 1 to 3 mu^2 a_moon (75 to 170 km) for a flyby at 0.5 mu a_moon, and the second-order error is not reliably smaller.
- Recommended test thresholds for `#899`'s seed (INFERRED, to be set before running): first-order periapsis error within 3 mu^2 a_moon; second-order error less than first-order only if the first-order Lambert (initial-velocity) correction is small (the paper's diagnostic); report both.

## 6. The full paper (AIAA 72-49, DOI 10.2514/6.1972-49) and the report

The synoptic says the full paper is available from AIAA and omits the data of the other trajectory types; the full conference paper very likely carries the tables for the interplanetary, midcourse and moon-to-Earth cases and more Earth-to-Moon cases (INFERRED, not seen). The fuller source is Lancaster, "Application of Matched Asymptotic Expansions to Lunar and Interplanetary Trajectories", McDonnell Douglas Report MDC G2748, February 1972 (2 vols; also a NASA Contractor Report, 1973), which the authors say contains "all of the data generated". Obtaining the NASA CR version (likely on NASA's technical reports server; not checked) would give the encounter speeds and the full case list.

## 7. Techniques applicable to the project's problems

- `#899`: the diagnostic and the prior above; a pre-registered pass test for a matched seed: compare the seed's periapsis, time and inclination errors with the first-order scale of 1 to 3 mu^2 a_moon and report whether second order helps; use the size of the first-order correction to the starting velocity as a predictor of second-order failure (the authors' finding).
- `#906`: none beyond noting the 36-min timing error at first order: a turn computed from a matched seed has a timing and radius error that the gate must not treat as physical.
- `#928`: none directly; the paper's computational speed claim (1.7 to 6.0 s on a CDC 6500) is of historical interest only.

## 8. Open points and follow-ups
- Obtain AIAA Paper 72-49 and, better, MDC G2748 / NASA CR (1973) for the full tables and encounter speeds.
- The sign convention, the model (restricted or N-body with the Sun) and the mass ratio are unstated.
